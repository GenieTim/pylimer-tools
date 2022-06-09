#ifndef MEHP_FORCE_RELAX2_H
#define MEHP_FORCE_RELAX2_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "cppoptlib/function.h"
#include "cppoptlib/solver/bfgs.h"
#include "cppoptlib/solver/gradient_descent.h"
#include "cppoptlib/solver/newton_descent.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <map>
#include <string>
#include <tuple>

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
      double* nodalDisplacements; /* nodal displacements */
      double* nodalForces;        /* nodal forces */
      double averageSpringLength; /* average spring length */
      long int nrOfLoops;         /* loops */
      // coordinates & coonectivity
      Eigen::VectorXd coordinates;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
      // interesting properties
      ArrayXb springIsActive;
    } Network;

    using FunctionXd = cppoptlib::function::Function<double>;

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

      void runForceRelaxation(long int maxNrOfSteps = 1000, // default: 10000
                              double tol = 1e-6,            // default: 1e-9
                              double Nb2spec = -1.0,
                              bool is2D = false,
                              double kappa = 1.0,
                              int deltaViolations = 5,
                              double gradientNormTol = 1e-5)
      {
        long int i, step, Nact, Mact, a, b;
        double *r, Fdef;
        double stress[3][3], r20, G0, r2, G, s20, s2, s20len, s2len;
        Network net = this->finalConfig;
        FILE* fp;
        const int M = this->universe.getMolecules(crosslinkerType).size();
        const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
        const double bM = this->universe.computeMeanBondLength();
        this->Nb2 = Nb2spec > 0 ? Nb2spec : N; // N * bM * bM;
        const int f =
          this->universe.determineFunctionalityPerType()[crosslinkerType];
        this->is2D = is2D;

        Eigen::VectorXd u = Eigen::VectorXd::Zero(3 * net.nrOfNodes);
        /* array allocation */
        r = (double*)calloc(3 * net.nrOfNodes, sizeof(double));

        /* initial */
        std::tie(s20, s20len) =
          computeStressAndSquareDistances(&net, u, stress, tol);
        Fdef = Residual(&net, u, r);
        for (r20 = 0., i = 0; i < 3 * net.nrOfNodes; i++) {
          r20 += r[i] * r[i];
        }
        G0 = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;

        /* force relaxation */
        double Gamma_eq = 0.0;
        double Rx2_sum = 0.0;
        double Ry2_sum = 0.0;
        double Rz2_sum = 0.0;

        MEHPForce totalForce = MEHPForce(this->is2D, &net, kappa);
        auto state = totalForce.Eval(u);

        // set the exit conditions
        cppoptlib::solver::State<MEHPForce::scalar_t> myStoppingSolverState;
        myStoppingSolverState.num_iterations = maxNrOfSteps;
        myStoppingSolverState.x_delta = static_cast<MEHPForce::scalar_t>(tol);
        myStoppingSolverState.x_delta_violations = deltaViolations;
        myStoppingSolverState.f_delta = static_cast<MEHPForce::scalar_t>(tol);
        myStoppingSolverState.f_delta_violations = deltaViolations;
        myStoppingSolverState.gradient_norm =
          static_cast<MEHPForce::scalar_t>(gradientNormTol);
        myStoppingSolverState.condition_hessian =
          static_cast<MEHPForce::scalar_t>(0);
        myStoppingSolverState.status = cppoptlib::solver::Status::NotStarted;

        cppoptlib::solver::Bfgs<MEHPForce> solver =
          cppoptlib::solver::Bfgs<MEHPForce>(myStoppingSolverState);
        // cppoptlib::solver::NewtonDescent<MEHPForce> solver =
        //   cppoptlib::solver::NewtonDescent<MEHPForce>(myStoppingSolverState);
        // do run solver
        auto [solution, solverState] = solver.Minimize(totalForce, u);
        this->nrOfStepsDone = solverState.num_iterations;

        if (solverState.status == cppoptlib::solver::Status::IterationLimit) {
          this->exitReason = ExitReason::MAX_STEPS;
        } else if (solverState.status ==
                     cppoptlib::solver::Status::XDeltaViolation ||
                   solverState.status ==
                     cppoptlib::solver::Status::FDeltaViolation) {
          this->exitReason = ExitReason::TOLERANCE;
        } else if (solverState.status ==
                   cppoptlib::solver::Status::GradientNormViolation) {
          this->exitReason = ExitReason::GRADIENT_NORM;
        } else if (solverState.status ==
                   cppoptlib::solver::Status::HessianConditionViolation) {
          this->exitReason = ExitReason::HESSIAN_CONDITION;
        } else {
          this->exitReason = ExitReason::OTHER;
        }

        u = solution.x;

        std::cout << "Ran " << solverState.num_iterations
                  << " iterations to a tolerance of " << solverState.x_delta
                  << ", " << solverState.f_delta << std::endl;

        // TODO: evaluate solution properties properly

        Gamma_eq = 0.0;
        Rx2_sum = 0.0;
        Ry2_sum = 0.0;
        Rz2_sum = 0.0;
        Fdef = Residual(&net, u, r);
        for (r2 = 0., i = 0; i < 3 * net.nrOfNodes; i++) {
          r2 += r[i] * r[i];
          if (i % 3 == 0) {
            Rx2_sum += r[i] * r[i];
          } else if ((i + 1) % 3 == 0) {
            Ry2_sum += r[i] * r[i];
          } else {
            assert((i + 2) % 3 == 0);
            Rz2_sum += r[i] * r[i];
          }
        }
        Gamma_eq = r2 / ((double)net.nrOfSprings * Nb2);

        // TODO: this is incorrect.
        double R2_mean = Gamma_eq / (double)net.nrOfSprings;
        double Gamma_x = 3 * Rx2_sum / ((double)net.nrOfSprings * this->Nb2);
        double Gamma_y = 3 * Ry2_sum / ((double)net.nrOfSprings * this->Nb2);
        double Gamma_z = 3 * Rz2_sum / ((double)net.nrOfSprings * this->Nb2);

        /* acquire equilibrium properties */
        std::tie(s2, s2len) =
          computeStressAndSquareDistances(&net, u, stress, tol);
        G = (stress[0][0] + stress[1][1] + stress[2][2]) / 3.;
        this->finalPressure = G;

        /* output */
        // if (this->outputEndNodes) {
        //   fp = fopen(this->endNodesFile.c_str(), "w");
        //   fprintf(fp, "%.10f %.10f %.10f\n", net.L[0], net.L[1], net.L[2]);
        //   fprintf(fp, "%ld #nodes\n", net.nrOfNodes);
        //   fprintf(fp, "%ld #springs\n", net.nrOfSprings);
        //   for (i = 0; i < net.nrOfNodes; i++) {
        //     fprintf(fp,
        //             "%.16f %.16f %.16f\n",
        //             net.nodes[i].x + u[3 * i],
        //             net.nodes[i].y + u[3 * i + 1],
        //             net.nodes[i].z + u[3 * i + 2]);
        //   }
        //   fclose(fp);
        // }

        /* active springs */
        for (Nact = 0, i = 0; i < net.nrOfSprings; i++) {
          if (net.springIsActive[i] == 1) {
            ++Nact;
          }
        }

        /* active nodes */
        Eigen::VectorXi nrOfActiveNodesConnected = Eigen::VectorXi::Zero(net.nrOfNodes);
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

        /** array deallocation */
        free(r);
      };

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
        net->springIsActive = ArrayXb::Constant(net->nrOfSprings, false);

        int usualChainLen =
          this->universe.getMolecules(crosslinkerType)[0].getNrOfAtoms();

        // convert beads
        std::vector<pylimer_tools::entities::Atom> allAtoms =
          crosslinkerUniverse.getAtoms();
        std::map<int, int> atomIdToNode;
        for (int i = 0; i < allAtoms.size(); ++i) {
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
        for (int i = 0; i < net->nrOfSprings; ++i) {
          int atomIdFrom = allBonds["bond_from"][i];
          int atomIdTo = allBonds["bond_to"][i];
          net->averageSpringLength += usualChainLen;
          net->springIndexA[i] = atomIdToNode.at(atomIdFrom);
          net->springIndexB[i] = atomIdToNode.at(atomIdTo);
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
      double Residual(Network* net, const Eigen::VectorXd& u, double* r)
      {
        long int glob[6];
        double kappa, s[3], ua[3], ub[3], s2, fele[6];

        /* initial */
        double Fdef = 0.;
        for (size_t i = 0; i < 3 * net->nrOfNodes; i++) {
          r[i] = 0.;
        }

        /* assembly */

        // first, the distances
        Eigen::VectorXd actualCoordinates = net->coordinates + u;
        // It *could* be more efficient to index u instead of the coordinates
        Eigen::VectorXd coordinatesSpringEndA =
          actualCoordinates(net->springIndexA);
        Eigen::VectorXd coordinatesSpringEndB =
          actualCoordinates(net->springIndexB);
        Eigen::VectorXd springDistances =
          (coordinatesSpringEndA - coordinatesSpringEndB);
        if (this->is2D) {
          springDistances(Eigen::seq(0, Eigen::last, Eigen::fix<3>)) =
            Eigen::VectorXd::Zero(net->nrOfSprings / 3);
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
        double tol)
      {
        double kappa, s2 = 0, s2len = 0;

        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            stress[j][k] = 0.;
          }
        }

        /* assembly */

        // first, the distances
        Eigen::VectorXd actualCoordinates = net->coordinates + u;
        // It *could* be more efficient to index u instead of the coordinates
        Eigen::VectorXd coordinatesSpringEndA =
          actualCoordinates(net->springIndexA);
        Eigen::VectorXd coordinatesSpringEndB =
          actualCoordinates(net->springIndexB);
        Eigen::VectorXd springDistances =
          (coordinatesSpringEndA - coordinatesSpringEndB);
        if (this->is2D) {
          springDistances(Eigen::seq(0, Eigen::last, Eigen::fix<3>)) =
            Eigen::VectorXd::Zero(net->nrOfSprings / 3);
        }

        // then, the stresses
        for (size_t i = 0; i < net->nrOfSprings; i++) {
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
          if (s2 < tol) {
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

      class MEHPForce : public FunctionXd
      {
      private:
        bool is2D;
        Network* net;
        double kappa;
        double boxHalfs[3];
        // Eigen::VectorXd actualCoordinates;
        // Eigen::VectorXd coordinatesSpringEndA;
        // Eigen::VectorXd coordinatesSpringEndB;
        // Eigen::VectorXd springDistances;

      public:
        MEHPForce(bool is2D, Network* net, double kappa)
          : is2D(is2D)
          , net(net)
          , kappa(kappa)
        {
          boxHalfs[0] = 0.5 * net->L[0];
          boxHalfs[1] = 0.5 * net->L[1];
          boxHalfs[2] = 0.5 * net->L[2];

          // actualCoordinates = net->coordinates;
          // coordinatesSpringEndA = actualCoordinates(net->springIndexA);
          // coordinatesSpringEndB = actualCoordinates(net->springIndexB);
          // springDistances = (coordinatesSpringEndA - coordinatesSpringEndB);
        }

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        double operator()(const Eigen::VectorXd& u) const
        {
          Eigen::VectorXd actualCoordinates = net->coordinates + u;
          // It *could* be more efficient to index u instead of the coordinates
          Eigen::VectorXd coordinatesSpringEndA =
            actualCoordinates(net->springIndexA);
          Eigen::VectorXd coordinatesSpringEndB =
            actualCoordinates(net->springIndexB);
          Eigen::VectorXd springDistances =
            (coordinatesSpringEndA - coordinatesSpringEndB);
          if (this->is2D) {
            springDistances(Eigen::seq(0, Eigen::last, Eigen::fix<3>)) =
              Eigen::VectorXd::Zero(net->nrOfSprings / 3);
          }
          // Possibly improvable PBC
          for (size_t j = 0; j < 3 * net->nrOfSprings; j++) {
            while (springDistances[j] > boxHalfs[j % 3]) {
              springDistances[j] -= boxHalfs[j % 3];
            }
            while (springDistances[j] < boxHalfs[j % 3]) {
              springDistances[j] += boxHalfs[j % 3];
            }
          }
          double s2 = springDistances.squaredNorm();
          return 0.5 * kappa * s2;
        }
      };
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
