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
#ifdef CEREALIZABLE
#include "../utils/CerealUtils.h"
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

namespace pylimer_tools {
namespace sim {
  namespace mehp {
    class MEHPForceRelaxation
      : public pylimer_tools::sim::OutputSupportingSimulation
    {
    private:
#ifdef CEREALIZABLE
      MEHPForceRelaxation() {}; // not exposed to users, only used by Cereal

      friend class cereal::access;
#endif

      // state
      pylimer_tools::entities::Universe universe;
      MEHPForceEvaluator* forceEvaluator;
      bool simulationHasRun = false;
      bool simulationSuggestsRerun = false;
      bool outputEndNodes = false;
      ExitReason exitReason = ExitReason::UNSET;
      int nrOfStepsDone = 0;
      Network forceRelaxationNetwork;
      Eigen::VectorXd currentSpringDistances;
      Eigen::VectorXd currentVelocities;
      Eigen::VectorXd currentVelocitiesPlus;
      Eigen::VectorXd currentForces;
      // config
      SimpleSpringMEHPForceEvaluator
        springForceEvaluator; // helper for memory time
      bool is2D = false;
      int defaultNrOfChains = 0;
      double defaultR0Squared = 0.0;
      double suggestRerunEps = 1e-3;
      int crossLinkerType = 2;
      double dt = 1;

    public:
      MEHPForceRelaxation(const pylimer_tools::entities::Universe& u,
                          int crossLinkerType = 2,
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
        this->crossLinkerType = crossLinkerType;
        // interpret network already to be able to give early results
        Network net;
        ConvertNetwork(&net,
                       crossLinkerType,
                       remove2functionalCrosslinkers,
                       removeDanglingChains);
        // this->defaultR0Squared =
        //   universe.computeMeanSquareEndToEndDistance(crossLinkerType);
        this->defaultNrOfChains = net.springsContourLength.size();
        if (this->defaultNrOfChains > 0) {
          this->defaultR0Squared = net.springsContourLength.mean() *
                                   universe.computeMeanSquaredBondLength();
        }
        this->forceRelaxationNetwork = net;
        this->is2D = is2D;
        this->initializeDefaults();
        this->setForceEvaluator(forceEvaluator);
      };

      MEHPForceRelaxation(const Network& net,
                          bool is2D = false,
                          MEHPForceEvaluator* forceEvaluator = nullptr,
                          double kappa = 1.0)
      {
        this->forceRelaxationNetwork = net;

        if (forceEvaluator == nullptr) {
          this->springForceEvaluator = SimpleSpringMEHPForceEvaluator(kappa);
          forceEvaluator = &this->springForceEvaluator;
        }

        this->initializeDefaults();
        this->setForceEvaluator(forceEvaluator);
      }

#ifdef CEREALIZABLE
      static MEHPForceRelaxation constructFromString(std::string s)
      {
        MEHPForceRelaxation res = MEHPForceRelaxation();
        pylimer_tools::utils::deserializeFromString(res, s);
        return res;
      }
#endif

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

      double getVolume() override { return this->forceRelaxationNetwork.vol; }

      int getNrOfNodes() const
      {
        return this->forceRelaxationNetwork.nrOfNodes;
      }

      int getNrOfSprings() const
      {
        return this->forceRelaxationNetwork.nrOfSprings;
      }

      size_t getNumBonds() override { return this->getNrOfSprings(); }

      size_t getNumExtraBonds() override { return 0; }

      long int getNumBondsToForm() override { return 0; }

      size_t getNumAtoms() override { return this->getNrOfNodes(); }

      size_t getNumExtraAtoms() override { return 0; }

      Network getNetwork() const { return this->forceRelaxationNetwork; }

      void configAssumeBoxLargeEnough(bool assumption = true)
      {
        this->forceRelaxationNetwork.assumeBoxLargeEnough = assumption;
      }

      // MEHPForceEvaluator getForceEvaluator() const
      // {
      //   return *this->forceEvaluator;
      // }

      void setForceEvaluator(MEHPForceEvaluator* forceEvaluator)
      {
        this->forceEvaluator = forceEvaluator;
        this->forceEvaluator->setNetwork(this->forceRelaxationNetwork);
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
      int getNrOfActiveNodes(double tolerance = 0.05,
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
      double getSolubleWeightFraction(double tolerance = 0.05)
      {
        return this->computeSolubleWeightFraction(&this->forceRelaxationNetwork,
                                                  this->currentSpringDistances,
                                                  tolerance);
      }

      /**
       * @brief Count the number of atoms that are in any way connected to an
       * active spring
       *
       * @param tolerance
       * @return double
       */
      double countActiveClusteredAtoms(double tolerance = 0.05)
      {
        return this->countActiveClusteredAtoms(&this->forceRelaxationNetwork,
                                               this->currentSpringDistances,
                                               tolerance);
      }

      /**
       * @brief Get the Dangling Weight Fraction
       *
       * @param tolerance
       * @return double
       */
      double getDanglingWeightFraction(double tolerance = 0.05)
      {
        return this->computeDanglingWeightFraction(
          &this->forceRelaxationNetwork,
          this->currentSpringDistances,
          tolerance);
      }

      /**
       * @brief Get the cross-linker Chains that are active
       *
       * @param tolerance
       * @return std::vector<pylimer_tools::entities::Molecule>
       */
      std::vector<pylimer_tools::entities::Molecule> getActiveChains(
        double tolerance = 0.05) const
      {
        std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
          this->universe.getChainsWithCrosslinker(crossLinkerType);
        std::vector<pylimer_tools::entities::Molecule> resultingChains;
        Eigen::ArrayXb springIsActive =
          this->findActiveSprings(this->currentSpringDistances, tolerance);
        for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
          if (this->forceRelaxationNetwork.moleculeIdxToSpring[i] >= 0 &&
              springIsActive[this->forceRelaxationNetwork
                               .moleculeIdxToSpring[i]]) {
            resultingChains.push_back(crossLinkerChains[i]);
          }
        }
        return resultingChains;
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
        double tolerance = 0.05) const;

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
        double tolerance = 0.05,
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
        double tolerance = 0.05) const;

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
        return this->forceRelaxationNetwork.meanSpringContourLength;
      }

      Eigen::VectorXd getSpringContourLength() const
      {
        return this->forceRelaxationNetwork.springsContourLength;
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
       * @param b02 for the denominator, part of the melt <R_0^2> = b02 *
       * nrOfBondsInSpring
       * @param nrOfChains the nr of chains to average over (can be different
       * from the nr of springs thanks to omitted free chains or primary loops)
       * @return double
       */
      double getGammaFactor(double b02 = -1.0, int nrOfChains = -1) const;

      /**
       * @brief Get all the gamma factors for each spring
       *
       * @param b02 for the denominator, part of the melt <R_0^2> = b02 *
       * nrOfBondsInSpring
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd getGammaFactors(double b2 = 1.) const;

      int getNrOfIterations() const { return this->nrOfStepsDone; }

      ExitReason getExitReason() const { return this->exitReason; }

      /**
       * @brief Get the Spring Lengths
       *
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd getSpringLengths() const
      {
        Eigen::VectorXd springDistances = this->getSpringDistances();
        Eigen::VectorXd springLengths =
          Eigen::VectorXd::Zero(this->forceRelaxationNetwork.nrOfSprings);
        for (int i = 0; i < this->forceRelaxationNetwork.nrOfSprings; ++i) {
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
        return this->evaluateSpringDistances(&this->forceRelaxationNetwork,
                                             this->is2D);
      }

      /**
       * @brief Compute the spring lengths
       *
       * @param net the network to do the computation for
       * @return Eigen::VectorXd
       */
      static Eigen::VectorXd evaluateSpringDistances(const Network* net,
                                                     const bool is2D);
      static Eigen::VectorXd evaluateSpringDistances(
        const Network* net,
        const Eigen::VectorXd& displacement,
        const bool is2D);

      /**
       * @brief Return whether the simulation resulted in offsets close to the
       * limits
       *
       * @return true
       * @return false
       */
      bool suggestsRerun() const
      {
        return this->simulationSuggestsRerun || !this->simulationHasRun;
      }

      void configRerunEps(double eps = 1e-3) { this->suggestRerunEps = eps; }

      /**
       * @brief Get the Stress Tensor
       *
       * @return Eigen::Matrix3d
       */
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

#ifdef CEREALIZABLE
      void writeRestartFile(std::string& filename) override
      {
        throw std::runtime_error("No restarts allowed here, yet");
      }
#endif

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
        return this->forceRelaxationNetwork.coordinates;
      }

      double getTemperature() override
      {
        return -1; // TODO: implement?
      }

      size_t getNumParticles() override
      {
        return this->forceRelaxationNetwork.nrOfNodes;
      }

    protected:
      /**
       * @brief Convert the universe to a network
       *
       * @param net the target network
       * @param crossLinkerType the atom type of the crossLinker
       * @return true
       * @return false
       */
      bool ConvertNetwork(Network* net,
                          const int crossLinkerType = 2,
                          bool remove2functionalCrosslinkers = false,
                          bool removeDanglingChains = false)
      {
        std::vector<pylimer_tools::entities::Atom> springEndAtoms =
          this->universe.getAtomsOfType(crossLinkerType);

        if (remove2functionalCrosslinkers) {
          for (pylimer_tools::entities::Atom xlinker : springEndAtoms) {
            // change type of cross-linkers with a degree <= 2 to "normal",
            // non-cross-link beads
            size_t vertexId = this->universe.getIdxByAtomId(xlinker.getId());
            if (this->universe.computeFunctionalityForVertex(vertexId) <= 2) {
              this->universe.setPropertyValue(
                vertexId, "type", crossLinkerType - 1);
            }
          }
          springEndAtoms = this->universe.getAtomsOfType(crossLinkerType);
        }

        std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
          this->universe.getChainsWithCrosslinker(crossLinkerType);
        net->moleculeIdxToSpring =
          Eigen::VectorXi::Constant(crossLinkerChains.size(), -1);

        // need to include all but dangling and free chains in order to
        // model entanglement
        size_t nrOfSprings = 0;
        size_t omittedChainsAtoms = 0;
        size_t omittedChainsBonds = 0;
        std::vector<bool> vertexAdded =
          pylimer_tools::utils::initializeWithValue(
            this->universe.getNrOfAtoms(), false);
        for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
          std::vector<pylimer_tools::entities::Atom> endAtoms =
            crossLinkerChains[i].getChainEnds(crossLinkerType, true);
          for (pylimer_tools::entities::Atom endAtom : endAtoms) {
            long int endAtomVertexId =
              this->universe.getIdxByAtomId(endAtom.getId());
            if (endAtom.getType() != crossLinkerType &&
                !vertexAdded[endAtomVertexId]) {
              springEndAtoms.push_back(endAtom);
              vertexAdded[endAtomVertexId] = true;
            }
          }
          RUNTIME_EXP_IFN(crossLinkerChains[i].getType() !=
                            pylimer_tools::entities::MoleculeType::UNDEFINED,
                          "Cross-linker chain's chain type could not be "
                          "detected. Cannot work like that.");
          if (crossLinkerChains[i].getType() ==
                pylimer_tools::entities::MoleculeType::NETWORK_STRAND ||
              crossLinkerChains[i].getType() ==
                pylimer_tools::entities::MoleculeType::PRIMARY_LOOP ||
              (!removeDanglingChains &&
               crossLinkerChains[i].getType() ==
                 pylimer_tools::entities::MoleculeType::DANGLING_CHAIN)) {
            net->moleculeIdxToSpring[i] = nrOfSprings;
            nrOfSprings += 1;
          } else {
            // assert(endAtoms.size() == 0); // can also be
            omittedChainsAtoms +=
              (crossLinkerChains[i].getNrOfAtoms() -
               endAtoms.size());
            omittedChainsBonds += crossLinkerChains[i].getNrOfBonds();
          }
        }

        size_t nrOfSpringEnds = springEndAtoms.size();

        // crossLinkerUniverse.simplify();
        pylimer_tools::entities::Box box = this->universe.getBox();
        net->L[0] = box.getLx();
        net->L[1] = box.getLy();
        net->L[2] = box.getLz();
        net->nrOfNodes = nrOfSpringEnds;
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
        net->springIsActive = Eigen::ArrayXb::Constant(net->nrOfSprings, false);
        net->springsContourLength = Eigen::VectorXd::Zero(net->nrOfSprings);
        Eigen::VectorXd targetDistances =
          Eigen::VectorXd::Zero(3 * net->nrOfSprings);

        net->springIndicesOfLinks.reserve(net->nrOfNodes);
        for (size_t i = 0; i < net->nrOfNodes; ++i) {
          net->springIndicesOfLinks.push_back(std::vector<size_t>());
        }

        // convert beads
        std::map<int, int> atomIdToNode;
        for (size_t i = 0; i < springEndAtoms.size(); ++i) {
          pylimer_tools::entities::Atom atom = springEndAtoms[i];
          atomIdToNode[atom.getId()] = i;
          net->oldAtomIds[i] = atom.getId();
          net->coordinates[3 * i + 0] = atom.getX();
          net->coordinates[3 * i + 1] = atom.getY();
          net->coordinates[3 * i + 2] = atom.getZ();
        }

        // convert springs
        size_t spring_idx = 0;
        for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
          std::vector<pylimer_tools::entities::Atom> xlinkersOfChain =
            crossLinkerChains[i].getAtomsOfType(crossLinkerType);
          std::vector<pylimer_tools::entities::Atom> endsOfChain =
            crossLinkerChains[i].getChainEnds(crossLinkerType, true);
          long int nodeIdxFrom = atomIdToNode.at(endsOfChain[0].getId());
          long int nodeIdxTo = atomIdToNode.at(endsOfChain[1].getId());
          bool addChain = false;
          if (crossLinkerChains[i].getType() ==
              pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
            addChain = true;
            // spring contour length = nr of bonds between two cross-linkers
            net->springsContourLength[spring_idx] =
              crossLinkerChains[i].getNrOfBonds();
          } else if (crossLinkerChains[i].getType() ==
                     pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
            addChain = true;

            net->springsContourLength[spring_idx] =
              crossLinkerChains[i].getNrOfBonds();
          } else if (crossLinkerChains[i].getType() ==
                     pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) {
            if (!removeDanglingChains) {
              // to keep dangling chains, we convert the trailing atom to a
              // cross-link
              net->springsContourLength[spring_idx] =
                crossLinkerChains[i].getNrOfBonds();
              addChain = true;
            }
          }

          if (addChain) {
            RUNTIME_EXP_IFN(net->moleculeIdxToSpring[i] == spring_idx,
                            "Expected spring mapping to be consistent.");
            std::vector<pylimer_tools::entities::Atom> allChainAtoms =
              crossLinkerChains[i].getAtomsLinedUp();
            targetDistances.segment(3 * spring_idx, 3) =
              crossLinkerChains[i].getOverallBondSum(crossLinkerType);

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

        RUNTIME_EXP_IFN(spring_idx == net->nrOfSprings,
                        "Expected nr of springs not fulfilled.");

        // box volume
        net->vol = net->L[0] * net->L[1] * net->L[2];
        if (net->nrOfSprings > 0) {
          net->meanSpringContourLength = net->springsContourLength.mean();
        }

        net->springBoxOffset =
          ((net->coordinates(net->springCoordinateIndexA) -
            net->coordinates(net->springCoordinateIndexB)) +
           targetDistances);
        // net->springBoxOffset = this->universe.getBox().getOffset(
        //   net->coordinates(net->springCoordinateIndexB) -
        //   net->coordinates(net->springCoordinateIndexA));

        // check whether spring contour lengths are what we want them to be
        size_t numCrosslinkers = springEndAtoms.size();
        // this->universe.countPropertyValue<int>("type", crossLinkerType);
        size_t contourSum = net->springsContourLength.sum();
        assert(contourSum + omittedChainsBonds ==
               this->universe.getNrOfBonds());
        assert(net->springsContourLength.size() == net->nrOfSprings);
        assert(contourSum - net->nrOfSprings + omittedChainsAtoms +
                 net->nrOfNodes ==
               this->universe.getNrOfAtoms());
        assert(numCrosslinkers == net->nrOfNodes);

        return true; // crossLinkerUniverse.getNrOfBonds() == net->nrOfSprings;
      };

      /**
       * @brief Compute the gamma factor from certain spring distances
       *
       * @param springDistances
       * @param b02 for the denominator, part of the melt <R_0^2> = b02 *
       * nrOfBondsInSpring
       * @param nrOfChains the nr of chains to average over (can be different
       * from the nr of springs thanks to omitted free chains or primary loops)
       * @return double
       */
      double evaluateGammaFactor(const Eigen::VectorXd& springDistances,
                                 double b02,
                                 int nrOfChains) const
      {
        return this->evaluateGammaFactors(springDistances, b02).sum() /
               static_cast<double>(nrOfChains);
      }

      Eigen::VectorXd evaluateGammaFactors(
        const Eigen::VectorXd& springDistances,
        const double b02) const
      {
        INVALIDARG_EXP_IFN(
          springDistances.size() ==
            this->forceRelaxationNetwork.springsContourLength.size() * 3,
          "Invalid sizes.");
        Eigen::VectorXd gammaFactors(springDistances.size() / 3);
        for (size_t i = 0; i < springDistances.size() / 3; ++i) {
          gammaFactors(i) =
            springDistances.segment(3 * i, 3).squaredNorm() /
            (this->forceRelaxationNetwork.springsContourLength(i) * b02);
        }
        return gammaFactors;
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
        auto stressTensor = this->evaluateStressTensor(
          springDistances, this->forceRelaxationNetwork.vol);
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
                                 const double tolerance = 0.05) const
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
        const double tolerance = 0.05) const
      {
        if (net->nrOfSprings * 3 != springDistances.size()) {
          throw std::invalid_argument(
            "Spring distances and network don't match");
        }
        if (net->nrOfSprings < 1) {
          return 1.;
        }
        // find all active springs
        Eigen::ArrayXb activeSprings =
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
       * @brief Count the number of atoms that can be considered part of an
       * active cluster, i.e., are somehow connected to an active spring
       *
       * @param net
       * @param springDistances
       * @param tolerance
       * @return double
       */
      double countActiveClusteredAtoms(Network* net,
                                       const Eigen::VectorXd& springDistances,
                                       const double tolerance = 0.05) const
      {
        INVALIDARG_EXP_IFN(net->nrOfSprings * 3 == springDistances.size(),
                           "Spring distances and network don't match");
        if (net->nrOfSprings < 1) {
          return 0.;
        }
        // find all active springs
        Eigen::ArrayXb activeSprings =
          this->findActiveSprings(springDistances, tolerance);
        if (activeSprings.count() == 0) {
          return 0.;
        }
        // then, iteratively walk along the springs to mark those as "active"
        // that are connected to active springs
        bool hadChanged = true;
        Eigen::ArrayXb nodeIsActive = Eigen::ArrayXb::Zero(net->nrOfNodes);
        while (hadChanged) {
          Eigen::ArrayXb oldActiveSprings = activeSprings;
          for (size_t i = 0; i < net->nrOfNodes; ++i) {
            bool anyActive = false;
            for (size_t spring_idx : net->springIndicesOfLinks[i]) {
              if (activeSprings[spring_idx]) {
                anyActive = true;
                break;
              }
            }

            nodeIsActive(i) = anyActive;
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
        double activeNodes = nodeIsActive.count();

        return ((allActiveAtomsPerChains).matrix().sum() + activeNodes);
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
        const double tolerance = 0.05) const
      {
        double nActiveClusteredAtoms =
          this->countActiveClusteredAtoms(net, springDistances, tolerance);
        // finally, normalise by the number of atoms.
        // NOTE: currently, the weight of the atoms is ignored
        return 1. - nActiveClusteredAtoms / this->universe.getNrOfAtoms();
      }

      /**
       * @brief Iterate all spring distances, mark active ones (length >
       * tolerance)
       *
       * @param springDistances
       * @param tolerance
       * @return Eigen::ArrayXb
       */
      Eigen::ArrayXb findActiveSprings(const Eigen::VectorXd& springDistances,
                                       const double tolerance = 0.05) const
      {
        Eigen::ArrayXb result =
          Eigen::ArrayXb::Constant(springDistances.size() / 3, false);
        for (size_t i = 0; i < springDistances.size() / 3; ++i) {
          result[i] =
            sqrt(springDistances[3 * i + 0] * springDistances[3 * i + 0] +
                 springDistances[3 * i + 1] * springDistances[3 * i + 1] +
                 springDistances[3 * i + 2] * springDistances[3 * i + 2]) >
            tolerance;
        }
        return result;
      }

      void initializeDefaults()
      {
        this->currentForces = Eigen::VectorXd::Zero(
          this->forceRelaxationNetwork.coordinates.size());
        this->currentVelocities = Eigen::VectorXd::Zero(
          this->forceRelaxationNetwork.coordinates.size());
        this->currentVelocitiesPlus = Eigen::VectorXd::Zero(
          this->forceRelaxationNetwork.coordinates.size());
        this->currentSpringDistances = this->evaluateSpringDistances(
          &this->forceRelaxationNetwork, this->is2D);
      }
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
