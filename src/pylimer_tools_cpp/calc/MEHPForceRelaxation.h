#ifndef MEHP_FORCE_RELAX2_H
#define MEHP_FORCE_RELAX2_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
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
      FAILURE,
      OTHER
    };

    // heavily inspired by Prof. Dr. Andrei Gusev's Code
    class MEHPForceRelaxation
    {

    public:
      MEHPForceRelaxation(const pylimer_tools::entities::Universe u,
                          int crosslinkerType = 2,
                          bool is2D = false,
                          MEHPForceEvaluator* forceEvaluator = nullptr,
                          double kappa = 1.0)
        : universe(u)
      {
        if (forceEvaluator == nullptr) {
          this->springForceEvaluator = SimpleSpringMEHPForceEvaluator(kappa);
          forceEvaluator = &this->springForceEvaluator;
        }
        this->crosslinkerType = crosslinkerType;
        // interpret network already to be able to give early results
        Network net;
        ConvertNetwork(&net, crosslinkerType);
        this->initialConfig = net;
        this->is2D = is2D;
        this->currentDisplacements =
          Eigen::VectorXd::Zero(net.coordinates.size());
        this->currentSpringDistances =
          this->evaluateSpringDistances(&net, this->currentDisplacements, is2D);
        this->defaultR0Squared =
          universe.computeMeanSquareEndToEndDistance(crosslinkerType);
        this->defaultNrOfChains =
          universe.getMolecules(this->crosslinkerType).size();
        this->setForceEvaluator(forceEvaluator);
      };

      /**
       * @brief Actually do run the simulation
       *
       * @param algorithm
       * @param maxNrOfSteps
       * @param xtol
       * @param ftol
       */
      void runForceRelaxation(const char* algorithm = "LD_MMA",
                              long int maxNrOfSteps = 50000, // default: 10000
                              double xtol = 1e-12,
                              double ftol = 1e-9);

      /**
       * @brief Get the universe consisting of cross-linkers only
       *
       * @param newCrosslinkerType the type to give the cross-linkers
       * @return pylimer_tools::entities::Universe
       */
      pylimer_tools::entities::Universe getCrosslinkerVerse(
        int newCrosslinkerType = 2) const;

      int getDefaultNrOfChains() const { return this->defaultNrOfChains; }

      double getDefaultR0Square() const { return this->defaultR0Squared; }

      double getVolume() const { return this->initialConfig.vol; }

      int getNrOfNodes() const { return this->initialConfig.nrOfNodes; }

      int getNrOfSprings() const { return this->initialConfig.nrOfSprings; }

      // MEHPForceEvaluator getForceEvaluator() const
      // {
      //   return *this->forceEvaluator;
      // }

      void setForceEvaluator(MEHPForceEvaluator* forceEvaluator)
      {
        this->forceEvaluator = forceEvaluator;
        this->forceEvaluator->setNetwork(this->initialConfig);
        this->forceEvaluator->setIs2D(this->is2D);
        this->forceEvaluator->prepareForEvaluations();
      }

      /**
       * @brief Get the Nr Of Active Nodes
       *
       * @param tolerance  the tolerance: springs under a certain length are
       * considered inactive
       * @return int
       */
      int getNrOfActiveNodes(double tolerance = 0.1,
                             int minimumNrOfActiveConnections = 2,
                             int maximumNrOfActiveConnections = -1) const
      {
        return this
          ->getIdsOfActiveNodes(tolerance,
                                minimumNrOfActiveConnections,
                                maximumNrOfActiveConnections)
          .size();
      }

      /**
       * @brief Get the Effective Functionality Of each node
       *
       * Returns the number of active springs connected to each atom, atomId
       * used as index
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return std::unordered_map<long int, int>
       */
      std::unordered_map<long int, int> getEffectiveFunctionalityOfAtoms(
        double tolerance = 0.1) const;

      /**
       * @brief Get the Ids Of active Nodes
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @param minimumNrOfActiveConnections the number of active springs
       * required for this node to qualify as active
       * @return std::vector<long int> the atom ids
       */
      std::vector<long int> getIdsOfActiveNodes(
        double tolerance = 0.1,
        int minimumNrOfActiveConnections = 2,
        int maximumNrOfActiveConnections = -1) const;

      /**
       * @brief Get the Nr Of Active Springs connected to each node
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return Eigen::VectorXi
       */
      Eigen::VectorXi getNrOfActiveSpringsConnected(
        double tolerance = 0.1) const;

      /**
       * @brief Get the Nr Of Active Springs object
       *
       * @param tol the tolerance: springs under a certain length are considered
       * inactive
       * @return int
       */
      int getNrOfActiveSprings(double tol = 0.1) const
      {
        return this->countNrOfActiveSprings(this->currentSpringDistances, tol);
      }

      /**
       * @brief Get the Average Spring Length at the current step
       *
       * @return double
       */
      double getAverageSpringLength() const;

      std::array<std::array<double, 3>, 3> getStressTensor() const;

      /**
       * @brief Get the Pressure
       *
       * @return double
       */
      double getPressure() const
      {
        return this->evaluatePressure(this->currentSpringDistances);
      }

      /**
       * @brief Get the Residual Norm at the current step
       *
       * @return double
       */
      double getResidualNorm() const;

      /**
       * @brief Get the residuals (gradient) at the current step
       *
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd getResiduals() const;

      /**
       * @brief Get the Force at the current step
       *
       * @return double
       */
      double getForce() const;

      /**
       * @brief Get the Gamma Factor at the current step
       *
       * @param r02 the melt <R_0^2>, for phantom = Nb^2
       * @param nrOfChains the nr of chains to average over (can be different
       * from the nr of springs thanks to omitted free chains or primary loops)
       * @return double
       */
      double getGammaFactor(double r02 = -1.0, int nrOfChains = -1) const;

      int getNrOfIterations() const { return this->nrOfStepsDone; }

      ExitReason getExitReason() const { return this->exitReason; }

      /**
       * @brief Compute the spring lenghts
       *
       * @param net the network to do the computation for
       * @param u the displacements on top of the network
       * @return Eigen::VectorXd
       */
      static Eigen::VectorXd evaluateSpringDistances(const Network* net,
                                                     const Eigen::VectorXd& u,
                                                     const bool is2D);

    protected:
      /**
       * @brief Convert the universe to a network
       *
       * @param net the target network
       * @param crosslinkerType the atom type of the crosslinker
       * @return true
       * @return false
       */
      bool ConvertNetwork(Network* net, const int crosslinkerType = 2)
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
        net->oldAtomIds = Eigen::ArrayXi::Zero(net->nrOfNodes);
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
          net->oldAtomIds[i] = atom.getId();
          net->coordinates[3 * i + 0] = atom.getX();
          net->coordinates[3 * i + 1] = atom.getY();
          net->coordinates[3 * i + 2] = atom.getZ();
        }

        // convert springs
        std::map<std::string, std::vector<long int>> allBonds =
          crosslinkerUniverse.getBonds();
        for (size_t i = 0; i < net->nrOfSprings; ++i) {
          int atomIdFrom = allBonds["bond_from"][i];
          int atomIdTo = allBonds["bond_to"][i];
          net->springIndexA[i] = atomIdToNode.at(atomIdFrom);
          net->springIndexB[i] = atomIdToNode.at(atomIdTo);
          for (size_t j = 0; j < 3; j++) {
            net->springCoordinateIndexA[3 * i + j] =
              atomIdToNode.at(atomIdFrom) * 3 + j;
            net->springCoordinateIndexB[3 * i + j] =
              atomIdToNode.at(atomIdTo) * 3 + j;
          }
        }

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
      double evaluateGammaFactor(const Eigen::VectorXd& springDistances,
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
       * @return double
       */
      double evaluatePressure(const Eigen::VectorXd& springDistances) const
      {
        auto stressTensor =
          this->evaluateStressTensor(springDistances, this->initialConfig.vol);
        return this->evaluatePressure(stressTensor);
      }

      /**
       * @brief Evaluate the pressure of the network at specific displacements
       *
       * @param net the network to evaluate the pressure for
       * @param u the displacements
       * @return double
       */
      double evaluatePressure(Network* net, const Eigen::VectorXd& u) const
      {
        auto stressTensor = this->evaluateStressTensor(net, u, -1);
        return this->evaluatePressure(stressTensor);
      }

      /**
       * @brief Evaluate the pressure from the stress tensor
       *
       * @param stressTensor
       * @return double
       */
      double evaluatePressure(
        const std::array<std::array<double, 3>, 3>& stressTensor) const
      {
        return (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) /
               3.0;
      }

      /**
       * @brief Compute the stress tensor
       *
       * @param net
       * @param u
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        const Eigen::VectorXd& springDistances,
        const double volume) const;

      /**
       * @brief Compute the stress tensor
       *
       * @param net
       * @param u
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        Network* net,
        const Eigen::VectorXd& u,
        const double loopTol) const;

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

    private:
      pylimer_tools::entities::Universe universe;
      MEHPForceEvaluator* forceEvaluator;

      SimpleSpringMEHPForceEvaluator
        springForceEvaluator; // helper for memory time
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
