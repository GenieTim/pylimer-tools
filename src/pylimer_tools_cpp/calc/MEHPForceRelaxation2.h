#ifndef MEHP_FORCE_RELAX2_H
#define MEHP_FORCE_RELAX2_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlopt.hpp>
#include <string>
#include <tuple>
#include <vector>
#include <cassert>

namespace pylimer_tools {
namespace calc {
  namespace mehp {
    enum ExitReason
    {
      UNSET,
      TOLERANCE,
      MAX_STEPS,
      GRADIENT_NORM,
      HESSIAN_CONDITION,
      OTHER
    };

    typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

    // improved structures using Eigen
    typedef struct _Network
    {
      double L[3];                /* box sizes */
      double vol;                 /* box volume */
      long int nrOfNodes;         /* number of nodes */
      long int nrOfSprings;       /* number of springs */
      double averageSpringLength; /* average spring length */
      long int nrOfLoops;         /* loops */
      // coordinates & coonectivity
      Eigen::VectorXd coordinates;
      Eigen::ArrayXi springCoordinateIndexA;
      Eigen::ArrayXi springCoordinateIndexB;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
      // interesting properties
      ArrayXb springIsActive;
    } Network;

    typedef struct _AdditionalFunctionParameters
    {
      Network* net;
      double kappa;
      bool is2D;
    } AdditionalFunctionParameters;

    // heavily inspired by Prof. Dr. Andrei Gusev's Code
    class MEHPForceRelaxation2
    {

    public:
      MEHPForceRelaxation2(const pylimer_tools::entities::Universe u,
                           int crosslinkerType = 2)
        : universe(u)
      {
        // interpret network already to be able to give early results
        Network net;
        ConvertNetwork(&net, crosslinkerType);
        this->finalConfig = net;
        this->crosslinkerType = crosslinkerType;
      };

      void configDoOutputSteps(const std::string outputFile,
                               const int outputFreq)
      {
        this->stepOutputFile = outputFile;
        this->stepOutputFrequency = outputFreq;
      }

      void configDoOutputFinalCoordinates(const std::string outputFile)
      {
        this->outputEndNodes = true;
        this->endNodesFile = outputFile;
      }

      double getVolume() { return this->finalConfig.vol; }

      int getNrOfNodes() { return this->finalConfig.nrOfNodes; }

      int getNrOfSprings() { return this->finalConfig.nrOfSprings; }

      int getNrOfActiveNodes() { return this->nrOfActiveNodes; }

      int getNrOfActiveSprings() { return this->nrOfActiveSprings; }

      double getAverageSpringLength() { return sqrt(this->R2Mean); }

      double getFinalPressure() { return this->finalPressure; }

      double getGammaEq() { return this->gammaEq; }

      double getSigmaX() { return this->sigmaX; }

      double getSigmaY() { return this->sigmaY; }

      double getSigmaZ() { return this->sigmaZ; }

      int getNrOfIterations() { return this->nrOfStepsDone; }

      double getNb2() { return this->Nb2; }

      ExitReason getExitReason() { return this->exitReason; }

      void runForceRelaxation(long int maxNrOfSteps = 10000, // default: 10000
                              double tol = 1e-9,             // default: 1e-9
                              double Nb2spec = -1.0,
                              bool is2D = false,
                              double kappa = 1.0,
                              int deltaViolations = 5,
                              double gradientNormTol = 1e-5,
                              double loopTol = 1e-2,
                              const char* algorithm = "LD_MMA")
      {
        long int i, step, Nact, Mact, a, b;
        std::vector<double> r;
        double stress[3][3], G0, r2, G, s20, s2, s20len, s2len;
        Network net = this->finalConfig;
        FILE* fp;
        const int M = this->universe.getMolecules(crosslinkerType).size();
        const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
        const double bM = this->universe.computeMeanBondLength();
        this->Nb2 = Nb2spec > 0 ? Nb2spec : N; // N * bM * bM;
        const int f =
          this->universe.determineFunctionalityPerType()[crosslinkerType];
        this->is2D = is2D;

        /* array allocation */
        r.reserve(3 * net.nrOfNodes);
        std::vector<double> u0 =
          pylimer_tools::utils::initializeWithValue(3 * net.nrOfNodes, 0.0);
        Eigen::VectorXd u = Eigen::VectorXd::Zero(3 * net.nrOfNodes);

        /* initial */
        std::tie(s20, s20len) =
          computeStressAndSquareDistances(&net, u, stress, kappa, loopTol);
        double Fdef = Residual(&net, u, r, kappa);
        double r20 = 0.;
        for (size_t i = 0; i < 3 * net.nrOfNodes; i++) {
          r20 += r[i] * r[i];
        }
        G0 = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;

        /* force relaxation */

        nlopt::opt opt(algorithm, 3 * net.nrOfNodes);

        AdditionalFunctionParameters params;
        params.net = &net;
        params.kappa = kappa;
        params.is2D = this->is2D;
        nlopt::func objectiveF = [](unsigned n,
                                    const double* x,
                                    double* grad,
                                    void* f_data) -> double {
          AdditionalFunctionParameters* fParams =
            static_cast<AdditionalFunctionParameters*>(f_data);
          return MEHPForceRelaxation2::evaluateForceSetGradient(
            fParams->net, fParams->kappa, fParams->is2D, n, x, grad, f_data);
        };
        opt.set_min_objective(objectiveF, &params);
        // set exit conditions
        opt.set_xtol_rel(tol);
        opt.set_ftol_rel(tol);
        opt.set_maxeval(maxNrOfSteps);
        // start/set/run minimization
        double minf;
        opt.optimize(u0, minf);
        // query solution & exit reason
        u = Eigen::Map<Eigen::VectorXd>(u0.data(), u0.size());

        std::cout << "Ran " << opt.get_numevals()
                  << " iterations to a tolerance of " << opt.get_xtol_rel()
                  << ", " << opt.get_ftol_rel() << ". U has norm " << u.norm()
                  << std::endl;

        this->exitReason = opt.get_numevals() >= maxNrOfSteps - 1
                             ? ExitReason::MAX_STEPS
                             : ExitReason::TOLERANCE;
        // TODO: evaluate solution properties properly
        
        Fdef = Residual(&net, u, r, kappa);

        /* acquire equilibrium properties */
        std::tie(s2, s2len) =
          computeStressAndSquareDistances(&net, u, stress, kappa, loopTol);
        G = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;
        this->finalPressure = G;

        /* active springs */
        for (Nact = 0, i = 0; i < net.nrOfSprings; i++) {
          if (net.springIsActive[i] == 1) {
            ++Nact;
          }
        }

        /* active nodes */
        Eigen::VectorXi nrOfActiveNodesConnected =
          Eigen::VectorXi::Zero(net.nrOfNodes);
        for (i = 0; i < net.nrOfSprings; i++) {
          if (net.springIsActive[i] == true) /* active spring */
          {
            a = net.springIndexA[i];
            b = net.springIndexB[i];
            ++(nrOfActiveNodesConnected[a]);
            ++(nrOfActiveNodesConnected[b]);
          }
        }
        for (Mact = 0, i = 0; i < net.nrOfNodes; i++) {
          if (nrOfActiveNodesConnected[i] >= 2) {
            ++Mact;
          }
        }

        /* save results */
        // this->initialG = G0;
        // this->initialS2 = s20;
        // this->finalS2 = s2;
        // this->finalS2len = s2len;
        // this->initialS2len = s20len;
        // this->finalG = G;
        this->gammaEq = s2len; // G;
        this->sigmaX = stress[0][0];
        this->sigmaY = stress[1][1];
        this->sigmaZ = stress[2][2];
        this->finalConfig = net;
        this->nrOfActiveNodes = Mact;
        this->nrOfActiveSprings = Nact;
      };

      /**
       * @brief This function does the step
       *
       * This function is public for testing purposes — TODO: find a better way
       *
       * @param net
       * @param kappa
       * @param is2D
       * @param n
       * @param x
       * @param grad
       * @param f_data
       * @return double
       */
      static double evaluateForceSetGradient(const Network* net,
                                             double kappa,
                                             bool is2D,
                                             unsigned n,
                                             const double* x,
                                             double* grad,
                                             void* f_data)
      {
        Eigen::Map<const Eigen::VectorXd> u =
          Eigen::Map<const Eigen::VectorXd>(x, n);

        double boxHalfs[3];
        boxHalfs[0] = 0.5 * net->L[0];
        boxHalfs[1] = 0.5 * net->L[1];
        boxHalfs[2] = 0.5 * net->L[2];

        assert(u.size() == net->coordinates.size());
        Eigen::VectorXd actualCoordinates = net->coordinates + u;
        // It *could* be more efficient to index u instead of the coordinates
        Eigen::VectorXd coordinatesSpringEndA =
          actualCoordinates(net->springCoordinateIndexA);
        Eigen::VectorXd coordinatesSpringEndB =
          actualCoordinates(net->springCoordinateIndexB);
        Eigen::VectorXd springDistances =
          (coordinatesSpringEndA - coordinatesSpringEndB);
        if (is2D) {
          // springDistances(Eigen::seq(2, Eigen::last, Eigen::fix<3>)) =
          //   Eigen::VectorXd::Zero(net->nrOfSprings / 3);
          for (size_t i = 2; i < net->nrOfSprings; i += 3) {
            springDistances[i] = 0.0;
          }
        }
        assert(springDistances.size() == net->nrOfSprings * 3);
        // Possibly improvable PBC
        for (size_t j = 0; j < 3 * net->nrOfSprings; ++j) {
          while (springDistances[j] > boxHalfs[j % 3]) {
            springDistances[j] -= boxHalfs[j % 3];
          }
          while (springDistances[j] < -boxHalfs[j % 3]) {
            springDistances[j] += boxHalfs[j % 3];
          }
        }
        double s2 = springDistances.norm();
        double constantMultiplier = 0.5 * kappa / (s2);
        if (grad != NULL) {
          for (size_t j = 0; j < net->nrOfSprings; ++j) {
            grad[j] = 0.0;
          }
          for (size_t j = 0; j < net->nrOfSprings; ++j) {
            int a = net->springIndexA[j];
            int b = net->springIndexB[j];
            int nrOfDim = is2D ? 2 : 3;
            for (size_t dir = 0; dir < nrOfDim; ++dir) {
              grad[3 * a + dir] +=
                springDistances[3 * j + dir] * constantMultiplier;
              grad[3 * b + dir] -=
                springDistances[3 * j + dir] * constantMultiplier;
            }
          }
        }
        std::cout << "Evaluated force to " << std::setprecision(15)
                  << 0.5 * kappa * s2 << std::endl;
        return 0.5 * kappa * s2;
      }

    protected:
      /**
       * @brief Convert the universe to a network
       *
       * @param net the target network
       * @param crosslinkerType the atom type of the crosslinker
       * @return true
       * @return false
       */
      bool ConvertNetwork(Network* net, int crosslinkerType = 2)
      {
        pylimer_tools::entities::Universe crosslinkerUniverse =
          this->universe.getNetworkOfCrosslinker(crosslinkerType);
        // crosslinkerUniverse.simplify();
        pylimer_tools::entities::Box box = crosslinkerUniverse.getBox();
        net->L[0] = box.getLx();
        net->L[1] = box.getLy();
        net->L[2] = box.getLz();
        net->nrOfNodes = crosslinkerUniverse.getNrOfAtoms();
        net->nrOfSprings = crosslinkerUniverse.getNrOfBonds();
        net->coordinates = Eigen::VectorXd::Zero(3 * net->nrOfNodes);
        net->springIndexA = Eigen::ArrayXi::Zero(net->nrOfSprings);
        net->springIndexB = Eigen::ArrayXi::Zero(net->nrOfSprings);
        net->springCoordinateIndexA =
          Eigen::ArrayXi::Zero(3 * net->nrOfSprings);
        net->springCoordinateIndexB =
          Eigen::ArrayXi::Zero(3 * net->nrOfSprings);
        net->springIsActive = ArrayXb::Constant(net->nrOfSprings, false);

        // convert beads
        std::vector<pylimer_tools::entities::Atom> allAtoms =
          crosslinkerUniverse.getAtoms();
        std::map<int, int> atomIdToNode;
        for (size_t i = 0; i < allAtoms.size(); ++i) {
          pylimer_tools::entities::Atom atom = allAtoms[i];
          atomIdToNode[atom.getId()] = i;
          net->coordinates[3 * i + 0] = atom.getX();
          net->coordinates[3 * i + 1] = atom.getY();
          net->coordinates[3 * i + 2] = atom.getZ();
        }

        // convert springs
        std::map<std::string, std::vector<long int>> allBonds =
          crosslinkerUniverse.getBonds();
        net->averageSpringLength = 0;
        for (size_t i = 0; i < net->nrOfSprings; ++i) {
          int atomIdFrom = allBonds["bond_from"][i];
          int atomIdTo = allBonds["bond_to"][i];
          net->averageSpringLength +=
            universe.getAtom(atomIdFrom)
              .distanceTo(universe.getAtom(atomIdTo), &box);
          net->springIndexA[i] = atomIdToNode.at(atomIdFrom);
          net->springIndexB[i] = atomIdToNode.at(atomIdTo);
          for (size_t j = 0; j < 3; j++) {
            net->springCoordinateIndexA[3 * i + j] =
              atomIdToNode.at(atomIdFrom) * 3 + j;
            net->springCoordinateIndexB[3 * i + j] =
              atomIdToNode.at(atomIdTo) * 3 + j;
          }
        }

        if (crosslinkerUniverse.getNrOfBonds() != net->nrOfSprings) {
          return false;
        };

        net->averageSpringLength /= net->nrOfSprings;

        /* box volume */
        net->vol = net->L[0] * net->L[1] * net->L[2];
        return true;
      };

      /**
       * @brief Adjust a vector of distances to lie within half the box
       *
       * @param s the distances
       * @param box the box lengths
       */
      void ImposePBC(double s[3], double box[3]) const
      {
        for (int i = 0; i < 3; i++) {
          double half = 0.5 * box[i];
          while (s[i] > half) {
            s[i] -= box[i];
          }
          while (s[i] < -half) {
            s[i] += box[i];
          }
        }

        return;
      }

      /**
       * @brief the residual is minus the assembled force vector
       *
       * @param net the network to compute the residual for
       * @param u the displacement vector
       * @param r the residual
       * @return the definitive force
       */
      double Residual(Network* net,
                      const Eigen::VectorXd& u,
                      std::vector<double>& r,
                      double kappa)
      {
        long int glob[6];
        double s[3], ua[3], ub[3], s2, fele[6];

        /* initial */
        double Fdef = 0.;
        for (size_t i = 0; i < 3 * net->nrOfNodes; i++) {
          r[i] = 0.;
        }

        /* assembly */
        // first, the distances
        assert(u.size() == net->coordinates.size());
        Eigen::VectorXd actualCoordinates = net->coordinates + u;
        // It *could* be more efficient to index u instead of the coordinates
        Eigen::VectorXd coordinatesSpringEndA =
          actualCoordinates(net->springCoordinateIndexA);
        Eigen::VectorXd coordinatesSpringEndB =
          actualCoordinates(net->springCoordinateIndexB);
        Eigen::VectorXd springDistances =
          (coordinatesSpringEndA - coordinatesSpringEndB);

        if (this->is2D) {
          // springDistances(Eigen::seq(2, Eigen::last, Eigen::fix<3>)) =
          //   Eigen::VectorXd::Zero(net->nrOfSprings / 3);
          for (size_t i = 2; i < net->nrOfSprings; i += 3) {
            springDistances[i] = 0.0;
          }
        }

        // then, the residuals
        for (size_t i = 0; i < net->nrOfSprings; i++) {
          double s[3] = { springDistances[3 * i + 0],
                          springDistances[3 * i + 1],
                          springDistances[3 * i + 2] };
          ImposePBC(s, net->L);
          int a = net->springIndexA[i];
          int b = net->springIndexB[i];

          /* energy update */
          s2 = 0.;
          for (size_t j = 0; j < 3; j++) {
            s2 += s[j] * s[j];
          }
          Fdef += 0.5 * kappa * s2;

          /* element force vector */
          for (size_t j = 0; j < 3; j++) {
            fele[j] = kappa * s[j];
            fele[j + 3] = -kappa * s[j];
          }

          /* global numbering scheme */
          for (size_t j = 0; j < 3; j++) {
            glob[j] = 3 * a + j;
            glob[j + 3] = 3 * b + j;
          }

          /* residual update */
          for (size_t j = 0; j < 6; j++) {
            r[glob[j]] -= fele[j];
          }
        };

        return Fdef;
      }

      std::pair<double, double> computeStressAndSquareDistances(
        Network* net,
        const Eigen::VectorXd& u,
        double stress[3][3],
        double kappa,
        double loopTol)
      {
        double s2 = 0, s2len = 0;

        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            stress[j][k] = 0.;
          }
        }

        /* assembly */
        // first, the distances
        assert(u.size() == net->coordinates.size());
        Eigen::VectorXd actualCoordinates = net->coordinates + u;
        // It *could* be more efficient to index u instead of the coordinates
        Eigen::VectorXd coordinatesSpringEndA =
          actualCoordinates(net->springCoordinateIndexA);
        Eigen::VectorXd coordinatesSpringEndB =
          actualCoordinates(net->springCoordinateIndexB);
        Eigen::VectorXd springDistances =
          (coordinatesSpringEndA - coordinatesSpringEndB);

        if (this->is2D) {
          // springDistances(Eigen::seq(2, Eigen::last, Eigen::fix<3>)) =
          //   Eigen::VectorXd::Zero(net->nrOfSprings / 3);
          for (size_t i = 2; i < net->nrOfSprings; i += 3) {
            springDistances[i] = 0.0;
          }
        }

        // then, the stresses
        for (size_t i = 0; i < net->nrOfSprings; ++i) {
          double s[3] = { springDistances[3 * i + 0],
                          springDistances[3 * i + 1],
                          springDistances[3 * i + 2] };
          ImposePBC(s, net->L);
          /* spring contribution to the overall stress tensor */
          for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
              stress[j][k] += kappa * s[j] * s[k];
            }
          }

          double s2local = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);

          /* update */
          s2 += s2local;
          s2len += s2local / std::sqrt(s2local);

          /* loop count */
          if (s2 < loopTol) {
            net->nrOfLoops++;
            net->springIsActive[i] = false;
          } else {
            net->springIsActive[i] = true;
          }
        };

        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            stress[j][k] /= net->vol;
          }
        }

        return std::make_pair(s2 / (double)net->nrOfSprings,
                              s2len / (double)net->nrOfSprings);
      }

    private:
      pylimer_tools::entities::Universe universe;
      bool is2D = false;
      int stepOutputFrequency = 0;
      std::string stepOutputFile;
      bool outputEndNodes = false;
      std::string endNodesFile;
      Network finalConfig;
      int nrOfActiveSprings = 0;
      int nrOfActiveNodes = 0;
      int crosslinkerType;
      int nrOfStepsDone = 0;
      double Nb2 = 0.0;
      double gammaEq = 0.0;
      double finalPressure = 0.0;
      double sigmaX = 0.0;
      double sigmaY = 0.0;
      double sigmaZ = 0.0;
      double R2Mean = 0.0;
      ExitReason exitReason = ExitReason::UNSET;
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
