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

namespace pylimer_tools {
namespace calc {
  namespace mehp {
    enum ExitReason
    {
      UNSET,
      F_TOLERANCE,
      X_TOLERANCE,
      MAX_STEPS,
      HESSIAN_CONDITION,
      OTHER
    };

    typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

    // improved structures using Eigen
    typedef struct Network
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
    };

    typedef struct AdditionalFunctionParameters
    {
      Network* net;
      double kappa;
      bool is2D;
    };

    // heavily inspired by Prof. Dr. Andrei Gusev's Code
    class MEHPForceRelaxation
    {

    public:
      MEHPForceRelaxation(const pylimer_tools::entities::Universe u,
                          int crosslinkerType = 2)
        : universe(u)
      {
        // interpret network already to be able to give early results
        Network net;
        ConvertNetwork(&net, crosslinkerType);
        this->initialConfig = net;
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

      double getVolume() { return this->initialConfig.vol; }

      int getNrOfNodes() { return this->initialConfig.nrOfNodes; }

      int getNrOfSprings() { return this->initialConfig.nrOfSprings; }

      int getNrOfActiveNodes() { return this->nrOfActiveNodes; }

      int getNrOfActiveSprings() { return this->nrOfActiveSprings; }

      double getAverageSpringLength() { return sqrt(this->R2Mean); }

      double getInitialPressure() { return this->initialPressure; }

      double getFinalPressure() { return this->finalPressure; }

      double getInitialResidualNorm() { return this->initialResidualNorm; }

      double getFinalResidualNorm() { return this->finalResidualNorm; }

      double getInitialForce() { return this->initialForce; }

      double getFinalForce() { return this->finalForce; }

      double getGammaEq() { return this->gammaEq; }

      double getSigmaX() { return this->sigmaX; }

      double getSigmaY() { return this->sigmaY; }

      double getSigmaZ() { return this->sigmaZ; }

      int getNrOfIterations() { return this->nrOfStepsDone; }

      double getNb2() { return this->Nb2; }

      ExitReason getExitReason() { return this->exitReason; }

      void runForceRelaxation(bool is2D = false,
                              double Nb2spec = -1.0,
                              const char* algorithm = "LD_MMA",
                              long int maxNrOfSteps = 50000, // default: 10000
                              double xtol = 1e-12,
                              double ftol = 1e-9,
                              double loopTol = 1e-2,
                              double kappa = 1.0)
      {
        long int i, step, Nact, Mact, a, b;
        double stress[3][3], s2, s20, s20len, s2len;

        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            stress[j][k] = 0.;
          }
        }

        Network net = this->initialConfig;
        const int M = this->universe.getMolecules(crosslinkerType).size();
        const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
        const double bM = this->universe.computeMeanBondLength();
        this->Nb2 = Nb2spec > 0 ? Nb2spec : N; // N * bM * bM;
        const int f =
          this->universe.determineFunctionalityPerType()[crosslinkerType];
        this->is2D = is2D;

        /* array allocation */
        std::vector<double> u0 =
          pylimer_tools::utils::initializeWithValue(3 * net.nrOfNodes, 0.0);
        Eigen::VectorXd u = Eigen::VectorXd::Zero(3 * net.nrOfNodes);

        /* initial evaluations */
        double* x0 = new double[3 * net.nrOfNodes];
        double* r0 = new double[3 * net.nrOfNodes];
        for (size_t i = 0; i < net.nrOfNodes * 3; ++i) {
          x0[i] = 0.0;
          r0[i] = 0.0;
        }

        std::tie(s20, s20len) =
          computeStressAndSquareDistances(&net, u, stress, kappa, loopTol);
        double G0 = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;
        double f0 = MEHPForceRelaxation::evaluateForceSetGradient(
          &net, kappa, is2D, 3 * net.nrOfNodes, x0, r0, NULL);
        this->initialForce = f0;
        delete[] x0;

        double r20 = 0.;
        for (size_t i = 0; i < 3 * net.nrOfNodes; i++) {
          r20 += r0[i] * r0[i];
        }
        this->initialResidualNorm = r20;
        delete[](r0);

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
          return MEHPForceRelaxation::evaluateForceSetGradient(
            fParams->net, fParams->kappa, fParams->is2D, n, x, grad, f_data);
        };
        opt.set_min_objective(objectiveF, &params);
        // set constraints to support more algorithms
        // std::vector<double> upperBounds;
        // upperBounds.reserve(3 * net.nrOfNodes);
        // std::vector<double> lowerBounds;
        // lowerBounds.reserve(3 * net.nrOfNodes);
        // for (size_t i = 0; i < net.nrOfNodes; ++i) {
        //   for (size_t dir = 0; dir < 3; ++dir) {
        //     upperBounds.push_back(net.L[dir] * 0.5);
        //     lowerBounds.push_back(-net.L[dir] * 0.5);
        //   }
        // }
        // opt.set_upper_bounds(upperBounds);
        // opt.set_lower_bounds(lowerBounds);
        // set exit conditions
        opt.set_xtol_rel(xtol);
        opt.set_ftol_rel(ftol);
        opt.set_maxeval(maxNrOfSteps);
        // start/set/run minimization
        double minf;
        nlopt::result res = opt.optimize(u0, minf);

        // query solution & exit reason
        assert(u0.size() == 3 * net.nrOfNodes);
        u = Eigen::Map<Eigen::VectorXd>(u0.data(), u0.size());
        this->finalDisplacement = u;
        double* residuals = new double[3 * net.nrOfNodes];
        double fFinal = MEHPForceRelaxation::evaluateForceSetGradient(
          &net, kappa, is2D, 3 * net.nrOfNodes, u, residuals, NULL);
        this->finalForce = fFinal;
        double r2 = 0.0;
        for (int i = 0; i < net.nrOfNodes * 3; i++) {
          r2 += residuals[i] * residuals[i];
        }
        this->finalResidualNorm = r2;
        delete[] residuals;

        std::cout << "Ran " << opt.get_numevals()
                  << " iterations to a tolerance of " << opt.get_xtol_rel()
                  << ", " << opt.get_ftol_rel() << ", exit result " << res
                  << ". U has norm " << u.squaredNorm() << ", minF is "
                  << fFinal << " from " << f0 << std::endl;

        this->exitReason = ExitReason::OTHER;
        if (res == nlopt::result::FTOL_REACHED) {
          this->exitReason = ExitReason::F_TOLERANCE;
        } else if (res == nlopt::result::XTOL_REACHED) {
          this->exitReason = ExitReason::X_TOLERANCE;
        } else if (res == nlopt::result::MAXEVAL_REACHED) {
          this->exitReason = ExitReason::MAX_STEPS;
        }
        // TODO: evaluate solution properties properly
        this->nrOfStepsDone = opt.get_numevals();
        /* acquire equilibrium properties */
        std::tie(s2, s2len) =
          computeStressAndSquareDistances(&net, u, stress, kappa, loopTol);

        double G = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;
        std::cout << "Pressure is " << G << " from " << G0 << ", s2 is " << s2
                  << " from " << s20 << ", s2len is " << s2len << " from "
                  << s20len << ", residual norm is " << r2 << " from " << r20
                  << std::endl;
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
        this->initialPressure = G0;
        // this->initialS2 = s20;
        // this->finalS2 = s2;
        // this->finalS2len = s2len;
        // this->initialS2len = s20len;
        // this->finalG = G;
        this->gammaEq = s2 / this->Nb2; // G;
        this->sigmaX = stress[0][0];
        this->sigmaY = stress[1][1];
        this->sigmaZ = stress[2][2];
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
                                             const double kappa,
                                             const bool is2D,
                                             const unsigned n,
                                             const double* x,
                                             double* grad,
                                             void* f_data)
      {
        Eigen::Map<const Eigen::VectorXd> u =
          Eigen::Map<const Eigen::VectorXd>(x, n);
        return evaluateForceSetGradient(net, kappa, is2D, n, u, grad, f_data);
      }

      /**
       * @brief This function does the step
       *
       * This function is public for testing purposes — TODO: find a better way
       *
       * @param net
       * @param kappa
       * @param is2D
       * @param n
       * @param u
       * @param grad
       * @param f_data
       * @return double
       */
      static double evaluateForceSetGradient(const Network* net,
                                             const double kappa,
                                             const bool is2D,
                                             const unsigned n,
                                             const Eigen::VectorXd& u,
                                             double* grad,
                                             void* f_data)
      {
        assert(n == net->nrOfNodes * 3);
        assert(u.size() == net->coordinates.size());
        Eigen::VectorXd springDistances = MEHPForceRelaxation::evaluateSpringDistances(net, u, is2D);

        double s2 = springDistances.squaredNorm();
        double constantMultiplier = kappa; // * 0.5 / s2;
        if (grad != NULL) {
          for (size_t j = 0; j < n; ++j) {
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
        double s2m = 0.0;
        for (size_t i = 0; i < 3 * net->nrOfSprings; ++i) {
          s2m += springDistances[i] * springDistances[i];
        }
        // std::cout << "Evaluated force to " << std::setprecision(15)
        //           << 0.5 * kappa * s2 << ", whereas manual gives "
        //           << s2m * 0.5 * kappa << std::endl;
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

      double evaluatePressure(Network* net,
                              const Eigen::VectorXd& u,
                              const double kappa)
      {
        auto stressTensor = this->evaluateStressTensor(net, u, kappa, -1);
        return this->evaluatePressure(stressTensor);
      }

      double evaluatePressure(
        const std::array<std::array<double, 3>, 3> stressTensor)
      {
        return (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) /
               3;
      }

      /**
       * @brief Compute the stress tensor (kappa * distance^2 in each direction combination)
       * 
       * @param net 
       * @param u 
       * @param kappa 
       * @param loopTol 
       * @return std::array<std::array<double, 3>, 3> 
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        Eigen::VectorXd springDistances,
        double kappa,
        double volume)
      {
        std::array<std::array<double, 3>, 3> stress;

        for (size_t i = 0; i < springDistances.size() / 3; ++i) {
          double s[3] = { springDistances[3 * i + 0],
                          springDistances[3 * i + 1],
                          springDistances[3 * i + 2] };
          /* spring contribution to the overall stress tensor */
          for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
              stress[j][k] += kappa * s[j] * s[k];
            }
          }
        }

        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            stress[j][k] /= volume;
          }
        }

        return stress;
      }

      /**
       * @brief Compute the stress tensor (kappa * distance^2 in each direction combination)
       * 
       * @param net 
       * @param u 
       * @param kappa 
       * @param loopTol 
       * @return std::array<std::array<double, 3>, 3> 
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        Network* net,
        const Eigen::VectorXd& u,
        const double kappa,
        const double loopTol)
      {
        std::array<std::array<double, 3>, 3> stress;

        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            stress[j][k] = 0.;
          }
        }

        Eigen::VectorXd springDistances = this->evaluateSpringDistances(net, u, this->is2D);

        return this->evaluateStressTensor(springDistances, kappa, net->vol);
      }

      /**
       * @brief Count how many of the springs are active (length > tolerance)
       * 
       * @param springDistances 
       * @param tolerance 
       * @return int 
       */
      int countNrOfActiveSprings(const Eigen::VectorXd springDistances,
                                 const double tolerance = 0.1) const
      {
        return (this->evaluateSpringActiveness(springDistances, tolerance) ==
                true)
          .count();
      }

      /**
       * @brief Iterate all spring distances, mark active ones (length >
       * tolerance)
       *
       * @param springDistances
       * @param tolerance
       * @return ArrayXb
       */
      ArrayXb evaluateSpringActiveness(const Eigen::VectorXd springDistances,
                                       const double tolerance = 0.1) const
      {
        ArrayXb result = ArrayXb::Constant(springDistances.size() / 3, false);
        for (size_t i = 0; i < springDistances.size() / 3; ++i) {
          result[i] =
            sqrt(springDistances[3 * i + 0] * springDistances[3 * i + 0] +
                 springDistances[3 * i + 1] * springDistances[3 * i + 1] +
                 springDistances[3 * i + 2] * springDistances[3 * i + 2]) >
            tolerance;
        }
        return result;
      }

      /**
       * @brief Compute the spring lenghts
       *
       * @param net the network to do the computation for
       * @param u the displacements on top of the network
       * @return Eigen::VectorXd
       */
      static Eigen::VectorXd evaluateSpringDistances(const Network* net,
                                              const Eigen::VectorXd& u, const bool is2D)
      {
        double boxHalfs[3];
        boxHalfs[0] = 0.5 * net->L[0];
        boxHalfs[1] = 0.5 * net->L[1];
        boxHalfs[2] = 0.5 * net->L[2];
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

        if (is2D) {
          // springDistances(Eigen::seq(2, Eigen::last, Eigen::fix<3>)) =
          //   Eigen::VectorXd::Zero(net->nrOfSprings / 3);
          for (size_t i = 2; i < 3 * net->nrOfSprings; i += 3) {
            springDistances[i] = 0.0;
          }
        }
        assert(springDistances.size() == net->nrOfSprings * 3);

        // Possibly improvable PBC
        for (size_t j = 0; j < 3 * net->nrOfSprings; ++j) {
          while (springDistances[j] > boxHalfs[j % 3]) {
            springDistances[j] -= net->L[j % 3];
          }
          while (springDistances[j] < -boxHalfs[j % 3]) {
            springDistances[j] += net->L[j % 3];
          }
        }

        return springDistances;
      }

      std::pair<double, double> computeStressAndSquareDistances(
        Network* net,
        const Eigen::VectorXd& u,
        double (&stress)[3][3],
        const double kappa,
        const double loopTol)
      {
        double s2 = 0., s2len = 0.;

        Eigen::VectorXd springDistances = this->evaluateSpringDistances(net, u, this->is2D);

        // then, the stresses
        for (size_t i = 0; i < net->nrOfSprings; ++i) {
          double s[3] = { springDistances[3 * i + 0],
                          springDistances[3 * i + 1],
                          springDistances[3 * i + 2] };
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
          if (s2local < loopTol) {
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

        return std::make_pair(s2 / static_cast<double>(net->nrOfSprings),
                              s2len / static_cast<double>(net->nrOfSprings));
      }

    private:
      pylimer_tools::entities::Universe universe;
      bool is2D = false;
      int stepOutputFrequency = 0;
      std::string stepOutputFile;
      bool outputEndNodes = false;
      std::string endNodesFile;
      Network initialConfig;
      Eigen::VectorXd finalDisplacement;
      int nrOfActiveSprings = 0;
      int nrOfActiveNodes = 0;
      int crosslinkerType;
      int nrOfStepsDone = 0;
      double Nb2 = 0.0;
      double gammaEq = 0.0;
      double initialPressure = 0.0;
      double finalPressure = 0.0;
      double initialResidualNorm = 0.0;
      double finalResidualNorm = 0.0;
      double initialForce = 0.0;
      double finalForce = 0.0;
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
