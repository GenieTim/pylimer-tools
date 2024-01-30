#ifndef MEHP_FORCE_RELAX2_H
#define MEHP_FORCE_RELAX2_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include "OutputSupportingSimulation.h"
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
    class MEHPForceRelaxation
      : public pylimer_tools::calc::OutputSupportingSimulation
    {

    public:
      MEHPForceRelaxation(const pylimer_tools::entities::Universe& u,
                          int crosslinkerType = 2,
                          bool is2D = false,
                          MEHPForceEvaluator* forceEvaluator = nullptr,
                          double kappa = 1.0,
                          bool remove2functionalCrosslinkers = false,
                          bool removeDanglingChains = false)
        : universe(u)
      {
        if (forceEvaluator == nullptr) {
          this->springForceEvaluator = SimpleSpringMEHPForceEvaluator(kappa);
          forceEvaluator = &this->springForceEvaluator;
        }
        this->crosslinkerType = crosslinkerType;
        this->defaultR0Squared =
          universe.computeMeanSquareEndToEndDistance(crosslinkerType);
        this->defaultNrOfChains =
          universe.getMolecules(this->crosslinkerType).size();
        // interpret network already to be able to give early results
        Network net;
        ConvertNetwork(&net,
                       crosslinkerType,
                       remove2functionalCrosslinkers,
                       removeDanglingChains);
        this->initialConfig = net;
        this->is2D = is2D;
        this->currentDisplacements =
          Eigen::VectorXd::Zero(net.coordinates.size());
        this->currentForces = Eigen::VectorXd::Zero(net.coordinates.size());
        this->currentVelocities = Eigen::VectorXd::Zero(net.coordinates.size());
        this->currentVelocitiesPlus =
          Eigen::VectorXd::Zero(net.coordinates.size());
        this->currentSpringDistances =
          this->evaluateSpringDistances(&net, this->currentDisplacements, is2D);
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

      void runPhantomSteps(const long int nrOfSteps = 50000,
                           const double dt = 0.01,
                           const double kappa = 1.,
                           const double T = 1.,
                           const double gamma = 0.1);

      /**
       * @brief Get the universe consisting of cross-linkers only
       *
       * @param newCrosslinkerType the type to give the cross-linkers
       * @return pylimer_tools::entities::Universe
       */
      pylimer_tools::entities::Universe getCrosslinkerVerse() const;

      int getDefaultNrOfChains() const { return this->defaultNrOfChains; }

      double getDefaultR0Square() const { return this->defaultR0Squared; }

      double getVolume() override { return this->initialConfig.vol; }

      int getNrOfNodes() const { return this->initialConfig.nrOfNodes; }

      int getNrOfSprings() const { return this->initialConfig.nrOfSprings; }

      size_t getNumBonds() override { return this->getNrOfSprings(); }

      size_t getNumExtraBonds() override { return 0; }

      long int getNumBondsToForm() override { return 0; }

      size_t getNumAtoms() override { return this->getNrOfNodes(); }

      size_t getNumExtraAtoms() override { return 0; }

      Network getNetwork() const { return this->initialConfig; }

      void configAssumeBoxLargeEnough(bool assumption = true)
      {
        this->initialConfig.assumeBoxLargeEnough = assumption;
      }

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
       * @brief Get the Soluble Weight Fraction
       *
       * @param tolerance
       * @return double
       */
      double getSolubleWeightFraction(double tolerance = 0.1)
      {
        return this->computeSolubleWeightFraction(
          &this->initialConfig, this->currentSpringDistances, tolerance);
      }

      /**
       * @brief Get the Dangling Weight Fraction
       *
       * @param tolerance
       * @return double
       */
      double getDanglingWeightFraction(double tolerance = 0.1)
      {
        return this->computeDanglingWeightFraction(
          &this->initialConfig, this->currentSpringDistances, tolerance);
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

      Eigen::VectorXd getCurrentSpringDistances() const
      {
        return this->currentSpringDistances;
      }

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

      double getAverageContourLength() const
      {
        return this->initialConfig.meanSpringContourLength;
      }

      Eigen::VectorXd getSpringContourLength() const
      {
        return this->initialConfig.springsContourLength;
      }

      /**
       * @brief Get the Average Spring Length at the current step
       *
       * @return double
       */
      double getAverageSpringLength() const;

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

      Eigen::VectorXd getCurrentDisplacements() const
      {
        return this->currentDisplacements;
      }

      /**
       * @brief Get the Spring Lengths
       *
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd getSpringLengths() const
      {
        Eigen::VectorXd springDistances = this->getSpringDistances();
        Eigen::VectorXd springLengths =
          Eigen::VectorXd::Zero(this->initialConfig.nrOfSprings);
        for (int i = 0; i < this->initialConfig.nrOfSprings; ++i) {
          springLengths[i] = springDistances.segment(3 * i, 3).norm();
        }
        return springLengths;
      }

      /**
       * @brief Get the Spring Distances
       *
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd getSpringDistances() const
      {
        return this->evaluateSpringDistances(
          &this->initialConfig, this->currentDisplacements, this->is2D);
      }

      /**
       * @brief Compute the spring lengths
       *
       * @param net the network to do the computation for
       * @param u the displacements on top of the network
       * @return Eigen::VectorXd
       */
      static Eigen::VectorXd evaluateSpringDistances(const Network* net,
                                                     const Eigen::VectorXd& u,
                                                     const bool is2D);

      Eigen::Matrix3d getStressTensor() override
      {
        std::array<std::array<double, 3>, 3> stressTensor =
          this->evaluateStressTensor(this->currentSpringDistances,
                                     this->getVolume());
        Eigen::Matrix3d result = Eigen::Matrix3d::Zero();
        for (size_t i = 0; i < 3; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            result(i, j) = stressTensor[i][j];
          }
        }
        return result;
      }

      double getTimestep() override { return this->dt; }

      double getCurrentTime(double currentStep) override
      {
        return this->dt * currentStep;
      }

      void writeRestartFile(std::string& filename) override
      {
        throw std::runtime_error("No restarts allowed here, yet");
      }

      int getNumShifts() override { return 0; }
      int getNumRelocations() override { return 0; }

      Eigen::VectorXd getBondLengths() override
      {
        Eigen::VectorXd lens =
          Eigen::VectorXd::Zero(this->currentSpringDistances.size() / 3);

        for (size_t i = 0; i < this->currentSpringDistances.size() / 3; ++i) {
          double b = lens.segment(3 * i, 3).norm();
          lens[i] = b;
        }

        return lens;
      }

      Eigen::VectorXd getCoordinates() override
      {
        return this->initialConfig.coordinates + this->currentDisplacements;
      }

      double getTemperature() override
      {
        return -1; // TODO: implement?
      }

      size_t getNumParticles() override
      {
        return this->initialConfig.nrOfNodes;
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
      bool ConvertNetwork(Network* net,
                          const int crosslinkerType = 2,
                          bool remove2functionalCrosslinkers = false,
                          bool removeDanglingChains = false)
      {
        std::vector<pylimer_tools::entities::Atom> xlinkers =
          this->universe.getAtomsOfType(crosslinkerType);

        if (remove2functionalCrosslinkers) {
          for (pylimer_tools::entities::Atom xlinker : xlinkers) {
            // change type of cross-linkers with a degree <= 2 to "normal",
            // non-cross-link beads
            size_t vertexId = this->universe.getIdxByAtomId(xlinker.getId());
            if (this->universe.computeFunctionalityForVertex(vertexId) <= 2) {
              this->universe.setPropertyValue(
                vertexId, "type", crosslinkerType - 1);
            }
          }
          xlinkers = this->universe.getAtomsOfType(crosslinkerType);
        }

        std::vector<pylimer_tools::entities::Molecule> crosslinkerChains =
          this->universe.getChainsWithCrosslinker(crosslinkerType);

        // need to include all but dangling and free chains in order to
        // model entanglement
        size_t nrOfSprings = 0;
        size_t omittedChainsAtoms = 0;
        for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
          std::vector<pylimer_tools::entities::Atom> endAtoms =
            crosslinkerChains[i].getAtomsOfType(crosslinkerType);
          RUNTIME_EXP_IFN(crosslinkerChains[i].getType() !=
                            pylimer_tools::entities::MoleculeType::UNDEFINED,
                          "Cross-linker chain's chain type could not be "
                          "detected. Cannot work like that.");
          if (crosslinkerChains[i].getType() ==
                pylimer_tools::entities::MoleculeType::NETWORK_STRAND ||
              crosslinkerChains[i].getType() ==
                pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
            nrOfSprings += 1;
          } else if (!removeDanglingChains &&
                     crosslinkerChains[i].getType() ==
                       pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) {
            std::vector<pylimer_tools::entities::Atom> endAtoms =
              crosslinkerChains[i].getAtomsOfDegree(1);
            RUNTIME_EXP_IFN(endAtoms.size() == 2,
                            "Expected a dangling chain to have two ends, got " +
                              std::to_string(endAtoms.size()) + ".");
            assert(XOR(endAtoms[0].getType() == crosslinkerType,
                       endAtoms[1].getType() == crosslinkerType));

            pylimer_tools::entities::Atom newXlink =
              (endAtoms[0].getType() == crosslinkerType) ? endAtoms[1]
                                                         : endAtoms[0];
            xlinkers.push_back(newXlink);
            nrOfSprings += 1;
          } else {
            // assert(endAtoms.size() == 0); // can also be
            omittedChainsAtoms +=
              (crosslinkerChains[i].getNrOfAtoms() - endAtoms.size());
          }
        }

        size_t nrOfXlinks = xlinkers.size();

        // crosslinkerUniverse.simplify();
        pylimer_tools::entities::Box box = this->universe.getBox();
        net->L[0] = box.getLx();
        net->L[1] = box.getLy();
        net->L[2] = box.getLz();
        net->nrOfNodes = nrOfXlinks;
        net->nrOfSprings = nrOfSprings;
        net->coordinates = Eigen::VectorXd::Zero(3 * net->nrOfNodes);
        net->oldAtomIds = Eigen::ArrayXi::Zero(net->nrOfNodes);
        net->springIndexA = Eigen::ArrayXi::Zero(net->nrOfSprings);
        net->springIndexB = Eigen::ArrayXi::Zero(net->nrOfSprings);
        net->springCoordinateIndexA =
          Eigen::ArrayXi::Zero(3 * net->nrOfSprings);
        net->springBoxOffset = Eigen::VectorXd::Zero(3 * net->nrOfSprings);
        net->springCoordinateIndexB =
          Eigen::ArrayXi::Zero(3 * net->nrOfSprings);
        net->springIsActive = ArrayXb::Constant(net->nrOfSprings, false);
        net->springsContourLength = Eigen::VectorXd::Zero(net->nrOfSprings);
        Eigen::VectorXd targetDistances =
          Eigen::VectorXd::Zero(3 * net->nrOfSprings);

        net->springIndicesOfLinks.reserve(net->nrOfNodes);
        for (size_t i = 0; i < net->nrOfNodes; ++i) {
          net->springIndicesOfLinks.push_back(std::vector<size_t>());
        }

        // convert beads
        std::map<int, int> atomIdToNode;
        for (size_t i = 0; i < xlinkers.size(); ++i) {
          pylimer_tools::entities::Atom atom = xlinkers[i];
          atomIdToNode[atom.getId()] = i;
          net->oldAtomIds[i] = atom.getId();
          net->coordinates[3 * i + 0] = atom.getX();
          net->coordinates[3 * i + 1] = atom.getY();
          net->coordinates[3 * i + 2] = atom.getZ();
        }

        // convert springs
        size_t spring_idx = 0;
        for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
          std::vector<pylimer_tools::entities::Atom> xlinkersOfChain =
            crosslinkerChains[i].getAtomsOfType(crosslinkerType);
          long int nodeIdxFrom;
          long int nodeIdxTo;
          bool addChain = false;
          if (crosslinkerChains[i].getType() ==
              pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
            assert(xlinkersOfChain.size() == 2);
            addChain = true;
            nodeIdxFrom = atomIdToNode.at(xlinkersOfChain[0].getId());
            nodeIdxTo = atomIdToNode.at(xlinkersOfChain[1].getId());

            // spring contour length = nr of bonds between two cross-linkers
            net->springsContourLength[spring_idx] =
              crosslinkerChains[i].getNrOfAtoms() - 1;
          } else if (crosslinkerChains[i].getType() ==
                     pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
            assert(xlinkersOfChain.size() == 1 ||
                   (xlinkersOfChain.size() == 2 &&
                    xlinkersOfChain[0].getId() == xlinkersOfChain[1].getId()));

            nodeIdxFrom = atomIdToNode.at(xlinkersOfChain[0].getId());
            nodeIdxTo = nodeIdxFrom;
            addChain = true;

            net->springsContourLength[spring_idx] =
              crosslinkerChains[i].getNrOfAtoms();
            if (xlinkersOfChain.size() == 2) {
              net->springsContourLength[spring_idx] =
                crosslinkerChains[i].getNrOfAtoms() - 1;
            }
          } else if (crosslinkerChains[i].getType() ==
                       pylimer_tools::entities::MoleculeType::DANGLING_CHAIN &&
                     !removeDanglingChains) {
            // to keep dangling chains, we convert the trailing atom to a
            // cross-link
            assert(xlinkersOfChain.size() == 1);
            std::vector<pylimer_tools::entities::Atom> endsOfChain =
              crosslinkerChains[i].getAtomsOfDegree(1);
            assert(endsOfChain.size() == 2);

            nodeIdxFrom = atomIdToNode.at(endsOfChain[0].getId());
            nodeIdxTo = atomIdToNode.at(endsOfChain[1].getId());
            net->springsContourLength[spring_idx] =
              crosslinkerChains[i].getNrOfAtoms() - 1;
            addChain = true;
          }

          if (addChain) {
            if (nodeIdxFrom > nodeIdxTo) {
              std::swap(nodeIdxFrom, nodeIdxTo);
            }
            std::vector<pylimer_tools::entities::Atom> allChainAtoms =
              crosslinkerChains[i].getAtomsLinedUp();
            targetDistances.segment(3 * spring_idx, 3) =
              crosslinkerChains[i].getOverallBondSum(crosslinkerType);

            if (atomIdToNode.at(allChainAtoms[0].getId()) == nodeIdxTo) {
              targetDistances.segment(3 * spring_idx, 3) *= -1;
            } else {
              assert(atomIdToNode.at(allChainAtoms[0].getId()) == nodeIdxFrom);
            }

            pylimer_tools::utils::addIfNotContained(
              net->springIndicesOfLinks[nodeIdxFrom], spring_idx);
            if (nodeIdxFrom != nodeIdxTo) {
              pylimer_tools::utils::addIfNotContained(
                net->springIndicesOfLinks[nodeIdxTo], spring_idx);
            }

            net->springIndexA[spring_idx] = nodeIdxFrom;
            net->springIndexB[spring_idx] = nodeIdxTo;
            for (size_t j = 0; j < 3; j++) {
              net->springCoordinateIndexA[3 * spring_idx + j] =
                nodeIdxFrom * 3 + j;
              net->springCoordinateIndexB[3 * spring_idx + j] =
                nodeIdxTo * 3 + j;
            }

            spring_idx += 1;
          }
        }

        // box volume
        net->vol = net->L[0] * net->L[1] * net->L[2];
        if (net->nrOfSprings > 0) {
          net->meanSpringContourLength = net->springsContourLength.mean();
        }

        net->springBoxOffset = ((net->coordinates(net->springCoordinateIndexA) -
                                net->coordinates(net->springCoordinateIndexB)) +
                               targetDistances);
        // net->springBoxOffset = this->universe.getBox().getOffset(
        //   net->coordinates(net->springCoordinateIndexB) -
        //   net->coordinates(net->springCoordinateIndexA));

        // check whether spring contour lengths are what we want them to be
        size_t numCrosslinkers = xlinkers.size();
        // this->universe.countPropertyValue<int>("type", crosslinkerType);
        assert((net->springsContourLength.array() -
                Eigen::ArrayXd::Ones(net->nrOfSprings))
                   .sum() +
                 omittedChainsAtoms + numCrosslinkers ==
               this->universe.getNrOfAtoms());

        return true; // crosslinkerUniverse.getNrOfBonds() == net->nrOfSprings;
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
       * @brief Compute the weight fraction of non-active springs
       *
       * @param net
       * @param springDistances
       * @param tolerance
       * @return double
       */
      double computeDanglingWeightFraction(
        Network* net,
        const Eigen::VectorXd& springDistances,
        const double tolerance = 0.1) const
      {
        if (net->nrOfSprings * 3 != springDistances.size()) {
          throw std::invalid_argument(
            "Spring distances and network don't match");
        }
        if (net->nrOfSprings < 1) {
          return 1.;
        }
        // find all active springs
        ArrayXb activeSprings =
          this->findActiveSprings(springDistances, tolerance);
        if (activeSprings.count() == 0) {
          return 1.;
        }
        // as of now, the springsContourLength is equal to the number of bonds
        // from cross-link to cross-link. therefore, the number of atoms of each
        // of these springs is one less
        Eigen::ArrayXd allActiveAtomsPerChains =
          activeSprings.cast<double>() *
          (net->springsContourLength.array() -
           Eigen::ArrayXd::Ones(net->nrOfSprings));
        // finally, normalise by the number of atoms.
        // NOTE: currently, the weight of the atoms is ignored
        return 1. - ((allActiveAtomsPerChains).matrix().sum() +
                     this->getNrOfActiveNodes()) /
                      this->universe.getNrOfAtoms();
      }

      /**
       * @brief Compute the weight fraction of springs connected to active
       * springs (any depth)
       *
       * @param net
       * @param springDistances
       * @param tolerance
       * @return double
       */
      double computeSolubleWeightFraction(
        Network* net,
        const Eigen::VectorXd& springDistances,
        const double tolerance = 0.1) const
      {
        if (net->nrOfSprings * 3 != springDistances.size()) {
          throw std::invalid_argument(
            "Spring distances and network don't match");
        }
        if (net->nrOfSprings < 1) {
          return 1.;
        }
        // find all active springs
        ArrayXb activeSprings =
          this->findActiveSprings(springDistances, tolerance);
        if (activeSprings.count() == 0) {
          return 1.;
        }
        // then, iteratively walk along the springs to mark those as "active"
        // that are connected to active springs
        bool hadChanged = true;
        while (hadChanged) {
          ArrayXb oldActiveSprings = activeSprings;
          for (size_t i = 0; i < net->nrOfNodes; ++i) {
            bool anyActive = false;
            for (size_t spring_idx : net->springIndicesOfLinks[i]) {
              if (activeSprings[spring_idx]) {
                anyActive = true;
                break;
              }
            }

            if (anyActive) {
              for (size_t spring_idx : net->springIndicesOfLinks[i]) {
                activeSprings[spring_idx] = true;
              }
            }
          }
          hadChanged = (oldActiveSprings.count() != activeSprings.count());
        }

        // as of now, the springsContourLength is equal to the number of bonds
        // from cross-link to cross-link. therefore, the number of atoms of each
        // of these springs is one less
        Eigen::ArrayXd allActiveAtomsPerChains =
          activeSprings.cast<double>() *
          (net->springsContourLength.array() -
           Eigen::ArrayXd::Ones(net->nrOfSprings));
        // finally, normalise by the number of atoms.
        // NOTE: currently, the weight of the atoms is ignored
        return 1. - ((allActiveAtomsPerChains).matrix().sum() +
                     this->getNrOfActiveNodes()) /
                      this->universe.getNrOfAtoms();
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
      Eigen::VectorXd currentVelocities;
      Eigen::VectorXd currentVelocitiesPlus;
      Eigen::VectorXd currentForces;
      int crosslinkerType;
      int nrOfStepsDone = 0;
      double dt = 1;
      ExitReason exitReason = ExitReason::UNSET;
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
