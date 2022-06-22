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
      OTHER
    };

    typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

    // improved structures using Eigen
    struct Network
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

    struct AdditionalFunctionParameters
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
                          int crosslinkerType = 2,
                          bool is2D = false)
        : universe(u)
      {
        // interpret network already to be able to give early results
        Network net;
        ConvertNetwork(&net, crosslinkerType);
        this->initialConfig = net;
        this->is2D = is2D;
        this->crosslinkerType = crosslinkerType;
        this->currentDisplacements =
          Eigen::VectorXd::Zero(net.coordinates.size());
        this->currentSpringDistances =
          this->evaluateSpringDistances(&net, this->currentDisplacements, is2D);
        this->defaultR0Squared =
          universe.computeMeanSquareEndToEndDistance(crosslinkerType);
        this->defaultNrOfChains =
          universe.getMolecules(this->crosslinkerType).size();
      };

      /**
       * @brief Actually do run the simulation
       *
       * @param algorithm
       * @param maxNrOfSteps
       * @param xtol
       * @param ftol
       * @param loopTol
       * @param kappa
       */
      void runForceRelaxation(const char* algorithm = "LD_MMA",
                              long int maxNrOfSteps = 50000, // default: 10000
                              double xtol = 1e-12,
                              double ftol = 1e-9,
                              double loopTol = 1e-2,
                              double kappa = 1.0)
      {
        this->simulationHasRun = true;
        double stress[3][3];

        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            stress[j][k] = 0.;
          }
        }

        Network net = this->initialConfig;
        const int M = this->universe.getMolecules(crosslinkerType).size();
        const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
        const double bM = this->universe.computeMeanBondLength();
        const int f =
          this->universe.determineFunctionalityPerType()[crosslinkerType];
        bool is2D = this->is2D;

        /* array allocation */
        std::vector<double> u0 =
          pylimer_tools::utils::initializeWithValue(3 * net.nrOfNodes, 0.0);
        Eigen::VectorXd u = Eigen::VectorXd::Zero(3 * net.nrOfNodes);

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
        std::vector<double> upperBounds;
        upperBounds.reserve(3 * net.nrOfNodes);
        std::vector<double> lowerBounds;
        lowerBounds.reserve(3 * net.nrOfNodes);
        for (size_t i = 0; i < net.nrOfNodes; ++i) {
          for (size_t dir = 0; dir < 3; ++dir) {
            upperBounds.push_back(net.L[dir] * 0.5);
            lowerBounds.push_back(-net.L[dir] * 0.5);
          }
        }
        opt.set_upper_bounds(upperBounds);
        opt.set_lower_bounds(lowerBounds);
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
        this->currentDisplacements = u;
        this->currentSpringDistances =
          this->evaluateSpringDistances(&net, this->currentDisplacements, is2D);

        this->exitReason = ExitReason::OTHER;
        if (res == nlopt::result::FTOL_REACHED) {
          this->exitReason = ExitReason::F_TOLERANCE;
        } else if (res == nlopt::result::XTOL_REACHED) {
          this->exitReason = ExitReason::X_TOLERANCE;
        } else if (res == nlopt::result::MAXEVAL_REACHED) {
          this->exitReason = ExitReason::MAX_STEPS;
        }
        this->nrOfStepsDone += opt.get_numevals();
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
        Eigen::VectorXd springDistances =
          MEHPForceRelaxation::evaluateSpringDistances(net, u, is2D);

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

      pylimer_tools::entities::Universe getCrosslinkerVerse(
        int newCrosslinkerType = 2)
      {
        // convert nodes & springs back to a universe
        pylimer_tools::entities::Universe xlinkUniverse =
          pylimer_tools::entities::Universe(this->universe.getBox());
        std::vector<long int> ids;
        std::vector<int> types = pylimer_tools::utils::initializeWithValue(
          this->initialConfig.nrOfNodes, crosslinkerType);
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        std::vector<int> zeros = pylimer_tools::utils::initializeWithValue(
          this->initialConfig.nrOfNodes, 0);
        ids.reserve(this->initialConfig.nrOfNodes);
        x.reserve(this->initialConfig.nrOfNodes);
        y.reserve(this->initialConfig.nrOfNodes);
        z.reserve(this->initialConfig.nrOfNodes);
        for (int i = 0; i < this->initialConfig.nrOfNodes; ++i) {
          x.push_back(this->initialConfig.coordinates[3 * i + 0] +
                      this->currentDisplacements[3 * i + 0]);
          y.push_back(this->initialConfig.coordinates[3 * i + 1] +
                      this->currentDisplacements[3 * i + 1]);
          z.push_back(this->initialConfig.coordinates[3 * i + 2] +
                      this->currentDisplacements[3 * i + 2]);
          ids.push_back(i + 1);
        }
        xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
        std::vector<long int> bondFrom;
        std::vector<long int> bondTo;
        bondFrom.reserve(this->initialConfig.nrOfSprings);
        bondTo.reserve(this->initialConfig.nrOfSprings);
        for (int i = 0; i < this->initialConfig.nrOfSprings; ++i) {
          bondFrom.push_back(this->initialConfig.springIndexA[i] + 1);
          bondTo.push_back(this->initialConfig.springIndexB[i] + 1);
        }
        xlinkUniverse.addBonds(
          bondFrom.size(),
          bondFrom,
          bondTo,
          pylimer_tools::utils::initializeWithValue(bondFrom.size(), 1),
          false,
          false); // disable simplify to keep the self-loops etc.
        return xlinkUniverse;
      }

      int getDefaultNrOfChains() { return this->defaultNrOfChains; }

      double getDefaultR0Square() { return this->defaultR0Squared; }

      double getVolume() { return this->initialConfig.vol; }

      int getNrOfNodes() { return this->initialConfig.nrOfNodes; }

      int getNrOfSprings() { return this->initialConfig.nrOfSprings; }

      /**
       * @brief Get the Nr Of Active Nodes
       *
       * @param tolerance  the tolerance: springs under a certain length are
       * considered inactive
       * @return int
       */
      int getNrOfActiveNodes(double tolerance = 0.1)
      {
        int Mact = 0;

        /* active nodes */
        Eigen::VectorXi nrOfActiveNodesConnected =
          Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
        ArrayXb springIsActive =
          this->findActiveSprings(this->currentSpringDistances, tolerance);
        for (size_t i = 0; i < this->initialConfig.nrOfSprings; i++) {
          if (springIsActive[i] == true) /* active spring */
          {
            int a = this->initialConfig.springIndexA[i];
            int b = this->initialConfig.springIndexB[i];
            ++(nrOfActiveNodesConnected[a]);
            ++(nrOfActiveNodesConnected[b]);
          }
        }
        for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
          if (nrOfActiveNodesConnected[i] >= 2) {
            ++Mact;
          }
        }

        return Mact;
      }

      /**
       * @brief Get the Nr Of Active Springs object
       *
       * @param tol the tolerance: springs under a certain length are considered
       * inactive
       * @return int
       */
      int getNrOfActiveSprings(double tol = 0.1)
      {
        return this->countNrOfActiveSprings(this->currentSpringDistances, tol);
      }

      /**
       * @brief Get the Average Spring Length at the current step
       *
       * @return double
       */
      double getAverageSpringLength()
      {
        double r2 = 0.0;
        for (int i = 0; i < this->initialConfig.nrOfSprings; i++) {
          double r2local = 0.0;
          for (int j = 0; j < 3; ++j) {
            r2local += this->currentSpringDistances[i * 3 + j] *
                       this->currentSpringDistances[i * 3 + j];
          }
          r2 += sqrt(r2local);
        }
        return r2 / this->initialConfig.nrOfSprings;
      }

      std::array<std::array<double, 3>, 3> getStressTensor(double kappa = 1.0)
      {
        return this->evaluateStressTensor(
          this->currentSpringDistances, kappa, this->initialConfig.vol);
      }

      /**
       * @brief Get the Pressure
       *
       * @param kappa the spring constant to use for the force
       * @return double
       */
      double getPressure(double kappa = 1.0)
      {
        return this->evaluatePressure(this->currentSpringDistances, kappa);
      }

      /**
       * @brief Get the Residual Norm at the current step
       *
       * @param kappa the spring constant to use for the force
       * @return double
       */
      double getResidualNorm(double kappa = 1.0)
      {
        double* r = new double[3 * this->initialConfig.nrOfNodes];
        for (size_t i = 0; i < this->initialConfig.nrOfNodes * 3; ++i) {
          r[i] = 0.0;
        }
        this->evaluateForceSetGradient(&this->initialConfig,
                                       kappa,
                                       this->is2D,
                                       3 * this->initialConfig.nrOfNodes,
                                       this->currentDisplacements,
                                       r,
                                       NULL);
        double r2 = 0.;
        for (size_t i = 0; i < 3 * this->initialConfig.nrOfNodes; i++) {
          r2 += r[i] * r[i];
        }

        delete[](r);
        return r2;
      }

      /**
       * @brief Get the Force at the current step
       *
       * @param kappa
       * @return double
       */
      double getForce(double kappa = 1.0)
      {
        return this->evaluateForceSetGradient(&this->initialConfig,
                                              kappa,
                                              this->is2D,
                                              3 * this->initialConfig.nrOfNodes,
                                              this->currentDisplacements,
                                              NULL,
                                              NULL);
      }

      /**
       * @brief Get the Gamma Factor at the current step
       *
       * @param r02 the melt <R_0^2>, for phantom = Nb^2
       * @param nrOfChains the nr of chains to average over (can be different
       * from the nr of springs thanks to omitted free chains or primary loops)
       * @return double
       */
      double getGammaFactor(double r02 = -1.0, int nrOfChains = -1)
      {
        if (r02 < 0) {
          r02 = this->defaultR0Squared;
        }
        if (nrOfChains < 1) {
          nrOfChains = this->defaultNrOfChains;
        }

        return this->evaluateGammaFactor(
          this->currentSpringDistances, r02, nrOfChains);
      }

      int getNrOfIterations() { return this->nrOfStepsDone; }

      ExitReason getExitReason() { return this->exitReason; }

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

        net->averageSpringLength /= net->nrOfSprings;

        // box volume
        net->vol = net->L[0] * net->L[1] * net->L[2];

        return crosslinkerUniverse.getNrOfBonds() == net->nrOfSprings;
      };

      /**
       * @brief Compute the gamma factor from certain spring distances
       *
       * @param springDistances
       * @param r02 the melt <R_0^2>, for phantom = Nb^2
       * @param nrOfChains the nr of chains to average over (can be different
       * from the nr of springs thanks to omitted free chains or primary loops)
       * @return double
       */
      double evaluateGammaFactor(Eigen::VectorXd& springDistances,
                                 double r02,
                                 int nrOfChains) const
      {
        return springDistances.squaredNorm() /
               (static_cast<double>(nrOfChains) * r02);
      }

      /**
       * @brief Evaluate the pressure of the network at specific spring
       * distances
       *
       * @param springDistances the spring distances
       * @param kappa the spring constant for the force factor
       * @return double
       */
      double evaluatePressure(Eigen::VectorXd& springDistances,
                              const double kappa) const
      {
        auto stressTensor = this->evaluateStressTensor(
          springDistances, kappa, this->initialConfig.vol);
        return this->evaluatePressure(stressTensor);
      }

      /**
       * @brief Evaluate the pressure of the network at specific displacements
       *
       * @param net the network to evaluate the pressure for
       * @param u the displacements
       * @param kappa the spring constant for the force factor
       * @return double
       */
      double evaluatePressure(Network* net,
                              const Eigen::VectorXd& u,
                              const double kappa) const
      {
        auto stressTensor = this->evaluateStressTensor(net, u, kappa, -1);
        return this->evaluatePressure(stressTensor);
      }

      /**
       * @brief Evaluate the pressure from the stress tensor
       *
       * @param stressTensor
       * @return double
       */
      double evaluatePressure(
        const std::array<std::array<double, 3>, 3> &stressTensor) const
      {
        return (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) /
               3.0;
      }

      /**
       * @brief Compute the stress tensor (kappa * distance^2 in each direction
       * combination)
       *
       * @param net
       * @param u
       * @param kappa
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        Eigen::VectorXd& springDistances,
        double kappa,
        double volume) const
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
       * @brief Compute the stress tensor (kappa * distance^2 in each direction
       * combination)
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
        const double loopTol) const
      {
        Eigen::VectorXd springDistances =
          this->evaluateSpringDistances(net, u, this->is2D);

        return this->evaluateStressTensor(springDistances, kappa, net->vol);
      }

      /**
       * @brief Count how many of the springs are active (length > tolerance)
       *
       * @param springDistances
       * @param tolerance
       * @return int
       */
      int countNrOfActiveSprings(const Eigen::VectorXd& springDistances,
                                 const double tolerance = 0.1) const
      {
        return (this->findActiveSprings(springDistances, tolerance) == true)
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
      ArrayXb findActiveSprings(const Eigen::VectorXd& springDistances,
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
                                                     const Eigen::VectorXd& u,
                                                     const bool is2D)
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

    private:
      pylimer_tools::entities::Universe universe;
      bool is2D = false;
      bool simulationHasRun = false;
      int stepOutputFrequency = 0;
      int defaultNrOfChains = 0;
      double defaultR0Squared = 0.0;
      std::string stepOutputFile;
      bool outputEndNodes = false;
      std::string endNodesFile;
      Network initialConfig;
      Eigen::VectorXd currentDisplacements;
      Eigen::VectorXd currentSpringDistances;
      int crosslinkerType;
      int nrOfStepsDone = 0;
      ExitReason exitReason = ExitReason::UNSET;
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
