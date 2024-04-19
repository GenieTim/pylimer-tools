#ifndef MEHP_FORCE_BALANCE_H
#define MEHP_FORCE_BALANCE_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include "EntanglementDetector.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include "OutputSupportingSimulation.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <nlopt.hpp>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {

    class MEHPForceBalance
      : public pylimer_tools::calc::OutputSupportingSimulation
    {

    public:
      MEHPForceBalance(const pylimer_tools::entities::Universe& u,
                       int crosslinkerType = 2,
                       bool is2D = false,
                       double kappa = 1.0,
                       bool remove2functionalCrosslinkers = false,
                       bool removeDanglingChains = false)
        : universe(u)
      {
        this->crosslinkerType = crosslinkerType;
        this->box = u.getBox();
        // interpret network already to be able to give early results
        ForceBalanceNetwork net;
        ConvertNetwork(net,
                       crosslinkerType,
                       remove2functionalCrosslinkers,
                       removeDanglingChains);
        this->initialConfig = net;
        this->is2D = is2D;
        this->currentDisplacements =
          Eigen::VectorXd::Zero(net.coordinates.size());
        this->completeInitialization();
      };

      /**
       * @brief Instantiate this simulator with randomly chosen slip-links.
       *
       * @param universe
       * @param nrOfSliplinksToSample
       * @param cutoff
       * @param minimumNrOfSliplinks
       * @param sameStrandCutoff
       * @param seed
       * @param crosslinkerType
       * @param is2D
       * @param kappa
       * @return MEHPForceBalance
       */
      static MEHPForceBalance constructWithSlipLinks(
        const pylimer_tools::entities::Universe& universe,
        pylimer_tools::calc::entanglement_detection::AtomPairEntanglements
          entanglements,
        int crosslinkerType = 2,
        bool is2D = false,
        double kappa = 1.0);

      /**
       * @brief Instantiate this simulator with randomly chosen slip-links.
       *
       * @param universe
       * @param nrOfSliplinksToSample
       * @param cutoff
       * @param minimumNrOfSliplinks
       * @param sameStrandCutoff
       * @param seed
       * @param crosslinkerType
       * @param is2D
       * @param kappa
       * @return MEHPForceBalance
       */
      static MEHPForceBalance constructWithRandomSlipLinks(
        const pylimer_tools::entities::Universe& universe,
        const size_t nrOfSliplinksToSample,
        const double cutoff,
        const size_t minimumNrOfSliplinks,
        const double sameStrandCutoff,
        const std::string seed = "",
        int crosslinkerType = 2,
        bool is2D = false,
        double kappa = 1.0)
      {
        // sample the "entanglements"
        pylimer_tools::calc::entanglement_detection::AtomPairEntanglements
          entanglements = pylimer_tools::calc::entanglement_detection::
            randomlyFindEntanglements(universe,
                                      nrOfSliplinksToSample,
                                      cutoff,
                                      minimumNrOfSliplinks,
                                      sameStrandCutoff,
                                      seed,
                                      crosslinkerType);

        RUNTIME_EXP_IFN(
          entanglements.pairsOfAtoms.size() >= minimumNrOfSliplinks,
          "Minimum number of slip-links could not be sampled: got " +
            std::to_string(entanglements.pairsOfAtoms.size()) + " instead of " +
            std::to_string(minimumNrOfSliplinks) + ".");

        return MEHPForceBalance::constructWithSlipLinks(
          universe, entanglements, crosslinkerType, is2D, kappa);
      }

      /**
       * @brief Finish initializing some member properties
       *
       */
      void completeInitialization()
      {
        this->currentSpringDistances = this->evaluateSpringVectors(
          this->initialConfig, this->currentDisplacements, is2D);
        this->currentPartialSpringDistances =
          this->evaluatePartialSpringVectors(
            this->initialConfig, this->currentDisplacements, is2D);
        this->defaultR0Squared =
          universe.computeMeanSquareEndToEndDistance(crosslinkerType);
        this->defaultNrOfChains =
          universe.getMolecules(this->crosslinkerType).size();
        this->validateNetwork();
      }

      /**
       * @brief Actually do run the simulation
       *
       * TODO: implement interruptability
       *
       * @param algorithm
       * @param maxNrOfSteps
       * @param xtol
       * @param ftol
       */
      void runForceRelaxation(
        BalanceRunMode mode = BalanceRunMode::ITERATIVE,
        double damping = 1.0,
        long int maxNrOfSteps = 50000, // default: 10000
        double xtol = 1e-9,
        const double initialResidualToUse = -1.0,
        const StructureSimplificationMode simplificationMode =
          StructureSimplificationMode::NO_SIMPLIFICATION,
        const double inactiveRemovalCutoff = -1.0,
        const int outputFrequency = 50,
        bool doInnerIterations = false,
        LinkSwappingMode allowSlipLinksToPassEachOther =
          LinkSwappingMode::NO_SWAPPING,
        const int swappingFrequency = 10,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1)
      {
        this->runForceRelaxation(
          mode,
          damping,
          maxNrOfSteps,
          xtol,
          initialResidualToUse,
          simplificationMode,
          inactiveRemovalCutoff,
          doInnerIterations,
          allowSlipLinksToPassEachOther,
          swappingFrequency,
          oneOverSpringPartitionUpperLimit,
          nrOfCrosslinkSwapsAllowedPerSliplink,
          []() { return false; },
          []() {});
      }

      /**
       * @brief Actually do run the simulation
       *
       * TODO: implement interruptability
       *
       * @param algorithm
       * @param maxNrOfSteps
       * @param xtol
       * @param ftol
       */
      void runForceRelaxation(
        BalanceRunMode mode,
        double damping,
        long int maxNrOfSteps, // default: 10000
        double xtol,
        const double initialResidualToUse,
        const StructureSimplificationMode simplificationMode,
        const double inactiveRemovalCutoff,
        bool doInnerIterations,
        LinkSwappingMode allowSlipLinksToPassEachOther,
        const int swappingFrequency,
        const double oneOverSpringPartitionUpperLimit,
        const int nrOfCrosslinkSwapsAllowedPerSliplink,
        const std::function<bool()>& shouldInterrupt,
        const std::function<void()>& cleanupInterrupt);

      /**
       * @brief Compute the spring update residual
       *
       * @param link_idx
       * @param displacements
       * @param springPartitions
       * @param oneOverSpringPartitionUpperLimit
       * @return double
       */
      double computePartitionUpdateZeroResidual(
        const ForceBalanceNetwork& net,
        const std::vector<size_t>& involvedPartitions,
        const size_t link_idx,
        const Eigen::VectorXd& displacements,
        Eigen::VectorXd& springPartitions,
        double oneOverSpringPartitionUpperLimit = 1.0) const
      {
        // TODO: revise, hard!
        assert(involvedPartitions.size() == 4);
        double firstMeanVal = 0.5 * (springPartitions[involvedPartitions[0]] +
                                     springPartitions[involvedPartitions[1]]);
        double secondMeanVal = 0.5 * (springPartitions[involvedPartitions[2]] +
                                      springPartitions[involvedPartitions[3]]);
        // Eigen::ArrayXi involvedCoordinateIndices = Eigen::ArrayXi(12);
        // for (size_t i = 0; i < 4; ++i) {
        //   involvedCoordinateIndices[3 * i] = 3 * involvedPartitions[i];
        //   involvedCoordinateIndices[3 * i + 1] = 3 * involvedPartitions[i] +
        //   1; involvedCoordinateIndices[3 * i + 2] = 3 * involvedPartitions[i]
        //   + 2;
        // }
        // Eigen::VectorXd displacementsBefore =
        //   displacements(involvedCoordinateIndices);
        Eigen::Vector4d partitionsBefore = springPartitions(involvedPartitions);
        springPartitions[involvedPartitions[0]] = firstMeanVal;
        springPartitions[involvedPartitions[1]] = firstMeanVal;
        springPartitions[involvedPartitions[2]] = secondMeanVal;
        springPartitions[involvedPartitions[3]] = secondMeanVal;
        // this->displaceToMeanPosition(
        //   this->initialConfig, displacements, springPartitions, link_idx);
        double retVal =
          this->updateSpringPartition(net,
                                      displacements,
                                      springPartitions,
                                      link_idx,
                                      oneOverSpringPartitionUpperLimit);
        // it seems to be faster to re-use memory rather than copying the whole
        // vectors
        springPartitions(involvedPartitions) = partitionsBefore;
        // displacements(involvedCoordinateIndices) = displacementsBefore;
        return retVal;
      }

      /**
       * @brief Remove cross-linkers, springs and associated slip-links with the
       * scheme suggested by Andrei
       *
       * @param net
       * @param displacements
       * @param springPartitions
       * @param tolerance
       * @return size_t
       */
      size_t doRemovalAndreisWay(ForceBalanceNetwork& net,
                                 Eigen::VectorXd& displacements,
                                 Eigen::VectorXd& springPartitions,
                                 double tolerance) const;

      /**
       * @brief Remove cross-links which do not have any springs with a certain
       * minimum length
       *
       * @param net
       * @param displacements
       * @param springPartitions
       * @param tolerance
       */
      size_t removeInactiveCrosslinks(ForceBalanceNetwork& net,
                                      Eigen::VectorXd& displacements,
                                      Eigen::VectorXd& springPartitions,
                                      double tolerance) const;

      /**
       * @brief Remove double listed springs from cross-links
       *
       * @param net
       */
      void cleanupPrimaryLoopsInStructure(ForceBalanceNetwork& net);

      /**
       * @brief Remove a spring (and all its parts, incl. slip-links) from the
       * structures
       *
       * @param net
       * @param springPartitions
       */
      void removeSpring(ForceBalanceNetwork& net,
                        Eigen::VectorXd& displacements,
                        Eigen::VectorXd& springPartitions,
                        const size_t springIdx) const;

      /**
       * @brief Remove a certain link from the structures
       *
       * @param net
       * @param displacements
       * @param linkIdx
       */
      void removeLink(ForceBalanceNetwork& net,
                      Eigen::VectorXd& displacements,
                      const size_t linkIdx) const;

      /**
       * @brief Merge two springs around a given cross-link
       *
       * @param net
       * @param springPartitions
       */
      void mergeSprings(ForceBalanceNetwork& net,
                        const Eigen::VectorXd& displacements,
                        Eigen::VectorXd& springPartitions,
                        const size_t removedSpringIdx,
                        const size_t keptSpringIdx,
                        const size_t linkToReduce) const;

      /**
       * @brief Merge two springs around a given cross-link
       *
       * This does not require the resulting network to be valid.
       *
       * @param net
       * @param springPartitions
       */
      void mergePartialSprings(ForceBalanceNetwork& net,
                               const Eigen::VectorXd& displacements,
                               Eigen::VectorXd& springPartitions,
                               const size_t removedSpringIdx,
                               const size_t keptSpringIdx,
                               const size_t linkToReduce,
                               bool skipEigenResize = false) const;

      /**
       * @brief Add a slip-link to a given partial spring
       *
       * This does not require the resulting network to be valid.
       *
       * @param net
       * @param springPartitions
       * @param partialSpringIdx
       * @param slipLinkIdx
       */
      size_t addSlipLinkToPartialSpring(ForceBalanceNetwork& net,
                                        const Eigen::VectorXd& displacements,
                                        Eigen::VectorXd& springPartitions,
                                        const size_t partialSpringIdx,
                                        const size_t slipLinkIdx,
                                        const double alpha) const;

      void relaxationLight(ForceBalanceNetwork& net,
                           Eigen::VectorXd& springPartitions,
                           Eigen::VectorXd& displacements,
                           const size_t linkIdx,
                           const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        Eigen::VectorXd oneOverSpringPartitions = Eigen::VectorXd::Zero(0);
        this->relaxationLight(net,
                              springPartitions,
                              oneOverSpringPartitions,
                              displacements,
                              linkIdx,
                              oneOverSpringPartitionUpperLimit);
      };

      void relaxationLight(ForceBalanceNetwork& net,
                           Eigen::VectorXd& springPartitions,
                           Eigen::VectorXd& oneOverSpringPartitions,
                           Eigen::VectorXd& displacements,
                           const size_t linkIdx,
                           const double oneOverSpringPartitionUpperLimit = 1.0);

      /**
       * @brief Replace the two springs traversinga a two-functional cross-links
       * with a single spring
       *
       * @param net
       * @param displacements
       * @param springPartitions
       */
      size_t removeTwofunctionalCrosslinks(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& displacements,
        Eigen::VectorXd& springPartitions) const;

      /**
       * @brief Add slip-links to this system
       *
       * @param sliplinkDensity
       * @param cutoff
       * @return size_t the nr of actually added slip-links
       */
      size_t randomlyAddSliplinks(const size_t nrOfSliplinksToSample,
                                  const double cutoff = 2.0,
                                  const size_t minimumNrOfSliplinks = 0,
                                  const double sameStrandCutoff = 2.0,
                                  const bool excludeCrosslinks = false,
                                  const int seed = -1);

      /**
       * @brief Add slip-links to this system based on entangled loops
       *
       * @return size_t the nr of actually added slip-links
       */
      size_t addSliplinksBasedOnCycles(const int maxLoopLength = -1);

      /**
       * @brief Deform the system to match the specified box
       *
       * @param box
       */
      void deformTo(pylimer_tools::entities::Box& newBox)
      {
        this->box.adjustCoordinatesTo(this->initialConfig.coordinates, newBox);
        this->box.adjustCoordinatesTo(this->currentDisplacements, newBox);
        this->box = newBox;
        this->universe.setBox(newBox, true);
        this->initialConfig.L[0] = this->box.getLx();
        this->initialConfig.L[1] = this->box.getLy();
        this->initialConfig.L[2] = this->box.getLz();
        this->initialConfig.boxHalfs[0] = 0.5 * this->initialConfig.L[0];
        this->initialConfig.boxHalfs[1] = 0.5 * this->initialConfig.L[1];
        this->initialConfig.boxHalfs[2] = 0.5 * this->initialConfig.L[2];
      }

      /**
       * @brief Investigate one parametrisation optimisation
       *
       * @param link_idx
       * @param displacements
       * @param springPartitions
       * @param innerMaxNrOfSteps
       * @param innerAlphaTol
       * @param distanceBackTolerance
       * @param residualNormSTolerance
       * @param innerMinNrOfSteps
       * @return std::tuple<Eigen::VectorXd, Eigen::VectorXd, size_t, double,
       * double, double>
       */
      std::tuple<Eigen::VectorXd,
                 Eigen::VectorXd,
                 size_t,
                 double,
                 double,
                 double,
                 double>
      inspectParametrisationOptimsationForLink(
        size_t link_idx,
        Eigen::VectorXd& displacements,
        Eigen::VectorXd& springPartitions,
        long int innerMaxNrOfSteps = 500,
        double innerAlphaTol = 1e-9,
        long int innerMinNrOfSteps = 1,
        const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        size_t innerIterationsDone = 0;
        double displacementDone = 0.0;
        double rOverr0 = 0.0;
        double r2 = 0.0;
        double r02 = this->computePartitionUpdateZeroResidual(
          this->initialConfig,
          this->getSpringpartitionIndicesOfSliplink(this->initialConfig,
                                                    link_idx),
          link_idx,
          displacements,
          springPartitions,
          oneOverSpringPartitionUpperLimit);
        do {
          r2 = this->updateSpringPartition(this->initialConfig,
                                           displacements,
                                           springPartitions,
                                           link_idx,
                                           oneOverSpringPartitionUpperLimit);
          rOverr0 = r2 / r02;
          displacementDone =
            this->displaceToMeanPosition(this->initialConfig,
                                         displacements,
                                         springPartitions,
                                         link_idx,
                                         oneOverSpringPartitionUpperLimit);
          innerIterationsDone += 1;
        } while ((innerIterationsDone < innerMaxNrOfSteps &&
                  rOverr0 > innerAlphaTol && std::isfinite(rOverr0)) ||
                 innerIterationsDone < innerMinNrOfSteps);
        return std::make_tuple(displacements,
                               springPartitions,
                               innerIterationsDone,
                               displacementDone,
                               rOverr0,
                               r02,
                               r2);
      }

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

      int getNrOfLinks() const { return this->initialConfig.nrOfLinks; }

      size_t getNumBonds() override { return this->getNrOfSprings(); }

      size_t getNumExtraBonds() override { return 0; }

      long int getNumBondsToForm() override { return 0; }

      size_t getNumAtoms() override { return this->getNrOfNodes(); }

      size_t getNumExtraAtoms() override
      {
        return this->getNrOfLinks() - this->getNrOfNodes();
      }

      int getNrOfSprings() const { return this->initialConfig.nrOfSprings; }

      int getNrOfPartialSprings() const
      {
        return this->initialConfig.nrOfPartialSprings;
      }

      int getNumIntraChainSlipLinks() const;

      Eigen::VectorXd getCurrentDisplacements() const
      {
        return this->currentDisplacements;
      }

      void setCurrentDisplacements(const Eigen::VectorXd displacements)
      {
        this->currentDisplacements = displacements;
      }

      void setSpringContourLengths(const Eigen::VectorXd springsContourLengths)
      {
        INVALIDARG_EXP_IFN(springsContourLengths.size() ==
                             this->initialConfig.springsContourLength.size(),
                           "Contour length must have the correct dimensions.");
        this->initialConfig.springsContourLength = springsContourLengths;
      }

      std::vector<Eigen::ArrayXi> getIndependentCoordinateSets(
        const ForceBalanceNetwork& net) const;

      std::pair<std::vector<Eigen::ArrayXi>, std::vector<Eigen::ArrayXi>>
      getHeuristicallyIndependentCoordinateSets(
        const ForceBalanceNetwork& net) const;

      std::vector<Eigen::ArrayXi> getRandomCoordinateSets(
        const ForceBalanceNetwork& net) const;

      void configAssumeBoxLargeEnough(bool assumption)
      {
        this->assumeBoxLargeEnough = assumption;

        this->currentSpringDistances = this->evaluateSpringVectors(
          this->initialConfig, this->currentDisplacements, is2D);
        this->currentPartialSpringDistances =
          this->evaluatePartialSpringVectors(
            this->initialConfig, this->currentDisplacements, is2D);
      }

      /**
       * @brief Get the Nr Of Active Nodes
       *
       * @param tolerance  the tolerance: springs under a certain length are
       * considered inactive
       * @return int
       */
      int getNrOfActiveNodes(double tolerance = 0.01,
                             int minimumNrOfActiveConnections = 2,
                             int maximumNrOfActiveConnections = -1,
                             bool usePartial = false) const
      {
        return this
          ->getIdsOfActiveNodes(tolerance,
                                minimumNrOfActiveConnections,
                                maximumNrOfActiveConnections,
                                usePartial)
          .size();
      }

      /**
       * @brief Get the Soluble Weight Fraction
       *
       * @param tolerance
       * @return double
       */
      double getSolubleWeightFraction(double tolerance = 0.01)
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
      double getDanglingWeightFraction(double tolerance = 0.01)
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
        double tolerance = 0.01) const;

      /**
       * @brief Compute the weight fraction of non-active springs
       *
       * @param net
       * @param springDistances
       * @param tolerance
       * @return double
       */
      double computeDanglingWeightFraction(
        ForceBalanceNetwork* net,
        const Eigen::VectorXd& springDistances,
        const double tolerance = 0.01) const
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
       * @brief Compute the weight fraction of springs connected to active
       * springs (any depth)
       *
       * @param net
       * @param springDistances
       * @param tolerance
       * @return double
       */
      double computeSolubleWeightFraction(
        ForceBalanceNetwork* net,
        const Eigen::VectorXd& springDistances,
        const double tolerance = 0.01) const
      {
        if (net->nrOfSprings * 3 != springDistances.size()) {
          throw std::invalid_argument(
            "Spring distances and network don't match");
        }
        if (net->nrOfSprings < 1) {
          return 1.;
        }
        std::vector<long int> activeCrosslinkIds =
          this->getIdsOfActiveNodes(tolerance, 2, -1);
        double activeClusterWeightFraction =
          this->universe.computeWeightFractionOfClustersAssociatedWith(
            activeCrosslinkIds);

        return 1. - activeClusterWeightFraction;
      }

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
        double tolerance = 0.01,
        int minimumNrOfActiveConnections = 2,
        int maximumNrOfActiveConnections = -1,
        bool usePartial = false) const;

      Eigen::VectorXd getCurrentSpringDistances() const
      {
        return this->currentSpringDistances;
      }

      std::vector<double> getCurrentSpringLengths() const
      {
        Eigen::VectorXd vecs = this->getCurrentSpringDistances();

        return pylimer_tools::utils::segmentwise_norm(vecs, 3);
      }

      std::vector<double> getOverallSpringLengths() const
      {
        Eigen::VectorXd partialSprings =
          this->getCurrentPartialSpringDistances();
        std::vector<double> results =
          std::vector<double>(this->initialConfig.nrOfSprings, 0.);
        for (size_t i = 0; i < this->initialConfig.nrOfPartialSprings; ++i) {
          results[this->initialConfig.partialToFullSpringIndex[i]] +=
            partialSprings.segment(3 * i, 3).norm();
        }

        return results;
      }

      Eigen::VectorXd getCurrentPartialSpringDistances() const
      {
        return this->evaluatePartialSpringVectors(
          this->initialConfig, this->currentDisplacements, this->is2D);
      }

      std::vector<double> getCurrentPartialSpringLengths() const
      {
        Eigen::VectorXd vecs = this->evaluatePartialSpringVectors(
          this->initialConfig, this->currentDisplacements, this->is2D);

        return pylimer_tools::utils::segmentwise_norm(vecs, 3);
      }

      /**
       * @brief Get the Nr Of Active Springs connected to each node
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return Eigen::VectorXi
       */
      Eigen::VectorXi getNrOfActiveSpringsConnected(
        double tolerance = 0.01) const;

      /**
       * @brief Get the Nr Of Active Springs connected to each node
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return Eigen::VectorXi
       */
      Eigen::VectorXi getNrOfActivePartialSpringsConnected(
        double tolerance = 0.01) const;

      /**
       * @brief Get the Nr Of Active Springs object
       *
       * @param tol the tolerance: springs under a certain length are considered
       * inactive
       * @return int
       */
      int getNrOfActiveSprings(double tolerance = 0.01) const
      {
        return this->countNrOfActiveSprings(this->currentSpringDistances,
                                            tolerance);
      }

      /**
       * @brief Get the Nr Of Active Springs object
       *
       * @param tol the tolerance: springs under a certain length are considered
       * inactive
       * @return int
       */
      int getNrOfActivePartialSprings(double tolerance = 0.01) const
      {
        return this->countNrOfActiveSprings(this->currentPartialSpringDistances,
                                            tolerance);
      }

      /**
       * @brief Get the Average Spring Length at the current step
       *
       * @return double
       */
      double getAverageSpringLength() const;

      Eigen::Matrix3d getStressTensor(
        const double oneOverSpringPartitionUpperLimit) const;

      Eigen::Matrix3d getStressTensorLinkBased(
        const double oneOverSpringPartitionUpperLimit = -1.0,
        const bool xlinksOnly = false) const;

      /**
       * @brief Get the Pressure
       *
       * @return double
       */
      double getPressure() const
      {
        return this->evaluatePressure(this->initialConfig,
                                      this->currentDisplacements,
                                      this->currentSpringPartitionsVec);
      }

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

      void addSlipLinks(const std::vector<size_t>& strandIdx1,
                        const std::vector<size_t>& strandIdx2,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& z)
      {
        std::vector<double> alphas;
        alphas.reserve(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
          alphas.push_back(0.5);
        }
        return this->addSlipLinks(
          strandIdx1, strandIdx2, x, y, z, alphas, alphas);
      }

      void addSlipLinks(const std::vector<size_t>& strandIdx1,
                        const std::vector<size_t>& strandIdx2,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& z,
                        const std::vector<double>& alpha1,
                        const std::vector<double>& alpha2,
                        bool clampAlpha = false)
      {
        std::vector<std::vector<size_t>> loops;
        std::vector<std::vector<size_t>> loopsOfSliplinks;
        return this->addSlipLinks(strandIdx1,
                                  strandIdx2,
                                  x,
                                  y,
                                  z,
                                  alpha1,
                                  alpha2,
                                  loops,
                                  loopsOfSliplinks,
                                  clampAlpha);
      }

      void addSlipLinks(const std::vector<size_t>& strandIdx1,
                        const std::vector<size_t>& strandIdx2,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& z,
                        const std::vector<double>& alpha1,
                        const std::vector<double>& alpha2,
                        std::vector<std::vector<size_t>> loops,
                        std::vector<std::vector<size_t>> loopsOfSliplinks,
                        bool clampAlpha = false);

      /**
       * @brief Compute the spring vectors
       *
       * @param net the network to do the computation for
       * @param u the displacements on top of the network
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd evaluateSpringVectors(const ForceBalanceNetwork& net,
                                            const Eigen::VectorXd& u,
                                            const bool is2D) const;

      /**
       * @brief Compute the spring lenghts
       *
       * @param net the network to do the computation for
       * @param u the displacements on top of the network
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd evaluateSpringLengths(const ForceBalanceNetwork& net,
                                            const Eigen::VectorXd& u,
                                            const bool is2D) const;

      /**
       * @brief Compute one spring length
       *
       * @param net
       * @param springIdx
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluatePartialSpringDistance(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const size_t springIdx,
        bool is2d = false) const
      {
        assert(net.isUpToDate);

        Eigen::Vector3d dist =
          ((net.coordinates.segment(3 * net.springPartIndexB(springIdx), 3) +
            u.segment(3 * net.springPartIndexB(springIdx), 3)) -
           (net.coordinates.segment(3 * net.springPartIndexA(springIdx), 3) +
            u.segment(3 * net.springPartIndexA(springIdx), 3))) +
          net.springPartBoxOffset.segment(3 * springIdx, 3);

        if (this->assumeBoxLargeEnough) {
          this->universe.getBox().handlePBC<Eigen::Vector3d>(dist);
        }

        if (is2d) {
          dist[2] = 0.0;
        }

        return dist;
      }

      Eigen::VectorXd evaluatePartialSpringVectors(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const bool is2D) const;

      size_t getOtherSpringIndex(const ForceBalanceNetwork& net,
                                 const size_t springIdx,
                                 const size_t linkIdx) const
      {
        assert(net.springPartIndexA[springIdx] == linkIdx ||
               net.springPartIndexB[springIdx] == linkIdx);
        return net.springPartIndexA[springIdx] == linkIdx
                 ? net.springPartIndexB[springIdx]
                 : net.springPartIndexA[springIdx];
      }

      /**
       * @brief Query the box offset for a specific spring in a specific
       * direction
       *
       * @param net
       * @param partialSpringIdx
       * @param linkIdx
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d getPartialSpringBoxOffsetTo(
        const ForceBalanceNetwork& net,
        const size_t partialSpringIdx,
        const size_t linkIdx) const
      {
        Eigen::Vector3d dist =
          net.springPartBoxOffset.segment(3 * partialSpringIdx, 3);
        if (net.springPartIndexA(partialSpringIdx) == linkIdx) {
          return -1. * dist;
        }
        return dist;
      }

      Eigen::Vector3d getPartialSpringBoxOffsetFrom(
        const ForceBalanceNetwork& net,
        const size_t partialSpringIdx,
        const size_t linkIdx) const
      {
        return -1. * this->getPartialSpringBoxOffsetTo(
                       net, partialSpringIdx, linkIdx);
      }

      /**
       * @brief Compute one spring length, in a specific direction
       *
       * @param net
       * @param springIdx
       * @param linkIdx the vector "target"
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluatePartialSpringDistanceTo(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const size_t springIdx,
        const size_t linkIdx,
        bool is2d = false) const
      {
        assert(this->isPartOfSpring(net, linkIdx, springIdx));

        Eigen::Vector3d dist =
          this->evaluatePartialSpringDistance(net, u, springIdx, is2d);

        return dist * (net.springPartIndexA(springIdx) == linkIdx ? -1. : 1.);
      }

      /**
       * @brief Compute one spring length, in a specific direction
       *
       * @param net
       * @param springIdx
       * @param linkIdx the vector "source"
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluatePartialSpringDistanceFrom(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const size_t springIdx,
        const size_t linkIdx,
        bool is2d = false) const
      {
        return -1. * this->evaluatePartialSpringDistanceTo(
                       net, u, springIdx, linkIdx, is2d);
      }

      bool validateNetwork() const
      {
        return this->validateNetwork(this->initialConfig,
                                     this->currentDisplacements,
                                     this->currentSpringPartitionsVec);
      }

      bool validateNetwork(const ForceBalanceNetwork& net) const
      {
        return this->validateNetwork(
          net, this->currentDisplacements, this->currentSpringPartitionsVec);
      }

      bool validateNetwork(const ForceBalanceNetwork& net,
                           const Eigen::VectorXd& u,
                           const Eigen::VectorXd& springPartitions) const;

      ForceBalanceNetwork getNetwork() { return this->initialConfig; }

      Eigen::VectorXd getSpringPartitions()
      {
        return this->currentSpringPartitionsVec;
      }

      void setSpringPartitions(const Eigen::VectorXd newSpringPartitionsVec)
      {
        this->currentSpringPartitionsVec = newSpringPartitionsVec;
      }

      /**
       * @brief List all the partial spring indices that are connected to a
       * specified (slip/cross)link
       *
       * @param net
       * @param linkIdx
       * @return std::unordered_set<size_t>
       */
      std::unordered_set<size_t> getPartialSpringIndicesOfLink(
        const ForceBalanceNetwork& net,
        const size_t linkIdx) const
      {
        INVALIDARG_EXP_IFN(linkIdx < net.nrOfLinks,
                           "The requested link does not exist");
        std::unordered_set<size_t> partialSpringIndices;

        std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];

        for (size_t spring_index = 0; spring_index < springIndices.size();
             ++spring_index) {
          std::vector<size_t> springsPartners =
            net.linkIndicesOfSprings[springIndices[spring_index]];
          for (size_t partner_idx = 0; partner_idx < springsPartners.size() - 1;
               ++partner_idx) {
            if (springsPartners[partner_idx] == linkIdx ||
                springsPartners[partner_idx + 1] == linkIdx) {
              size_t globalSpringIndex =
                net.localToGlobalSpringIndex[(springIndices[spring_index])]
                                            [partner_idx];
              partialSpringIndices.insert(globalSpringIndex);
            }
          }
        }
        return partialSpringIndices;
      }

      /**
       * @brief Evaluate the force on one link
       *
       * @param index the link index
       * @param oneOverSpringPartitionUpperLimit
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d getForceOn(
        const size_t index,
        const double oneOverSpringPartitionUpperLimit = 1.0) const
      {
        Eigen::VectorXi debugNrSpringsVisited =
          Eigen::VectorXi::Zero(this->initialConfig.nrOfPartialSprings);
        return this->evaluateForceOnLink(index,
                                         this->initialConfig,
                                         this->currentDisplacements,
                                         this->currentSpringPartitionsVec,
                                         debugNrSpringsVisited,
                                         1.0,
                                         oneOverSpringPartitionUpperLimit);
      }

      /**
       * @brief Evaluate the force on one link
       *
       * @param index the link index
       * @param oneOverSpringPartitionUpperLimit
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d getForceOn(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd&
          springPartitions, /* gives the parametrisation of N */
        const size_t index,
        const double oneOverSpringPartitionUpperLimit = 1.0) const
      {
        Eigen::VectorXi debugNrSpringsVisited = Eigen::VectorXi::Zero(0);
        return this->evaluateForceOnLink(index,
                                         net,
                                         u,
                                         springPartitions,
                                         debugNrSpringsVisited,
                                         1.0,
                                         oneOverSpringPartitionUpperLimit);
      }

      /**
       * @brief Evaluate the current stress on a particulas cross- or slip-link
       *
       * @param linkIdx
       * @param oneOverSpringPartitionUpperLimit
       * @return Eigen::Matrix3d
       */
      Eigen::Matrix3d getStressOn(
        const size_t linkIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0) const
      {
        Eigen::VectorXi debugNrSpringsVisited =
          Eigen::VectorXi::Zero(this->initialConfig.nrOfPartialSprings);
        return this->evaluateStressOnLink(linkIdx,
                                          this->initialConfig,
                                          this->currentDisplacements,
                                          this->currentSpringPartitionsVec,
                                          debugNrSpringsVisited,
                                          1.0,
                                          oneOverSpringPartitionUpperLimit);
      }

      /**
       * @brief Assemble all indices of partial springs for a particular
       * slip-link
       *
       * @param linkIdx
       * @return std::vector<size_t>
       */
      std::vector<size_t> getSpringpartitionIndicesOfSliplink(
        const ForceBalanceNetwork& net,
        const size_t linkIdx) const
      {
        std::vector<size_t> indices =
          pylimer_tools::utils::initializeWithValue<size_t>(4, 0);
        assert(indices.size() == 4);
        this->setSpringpartitionIndicesOfSliplink(indices, net, linkIdx);
        return indices;
      };

      /**
       * @brief Assemble all indices of partial springs for a particular
       * slip-link
       *
       * @param linkIdx
       * @return void
       */
      void setSpringpartitionIndicesOfSliplink(std::vector<size_t>& res_vec,
                                               const ForceBalanceNetwork& net,
                                               const size_t linkIdx) const;

      /**
       * @brief Updates the partition/parametrisation of a spring around one
       * link
       *
       */
      Eigen::VectorXd inspectSpringPartitionUpdate(const size_t linkIdx) const
      {
        Eigen::VectorXd springPartitions = this->currentSpringPartitionsVec;
        this->updateSpringPartition(this->initialConfig,
                                    this->currentDisplacements,
                                    springPartitions,
                                    linkIdx);
        return springPartitions;
      };

      /**
       * @brief Updates the partition/parametrisation of a spring around one
       * link
       *
       */
      double updateSpringPartition(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions, /* gives the parametrisation of N */
        const size_t linkIdx,
        double oneOverSpringPartitionUpperLimit = 1.0,
        bool allowSlipLinksToPassEachOther = false) const
      {
        Eigen::VectorXd oneOverSpringPartitions = Eigen::VectorXd::Zero(0);
        return this->updateSpringPartition(net,
                                           u,
                                           springPartitions,
                                           oneOverSpringPartitions,
                                           linkIdx,
                                           oneOverSpringPartitionUpperLimit,
                                           allowSlipLinksToPassEachOther);
      };

      /**
       * @brief Updates the partition/parametrisation of a spring around one
       * link
       *
       */
      double updateSpringPartition(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions, /* gives the parametrisation of N */
        Eigen::VectorXd&
          oneOverSpringPartitions, /* gives the parametrisation of N */
        const size_t linkIdx,
        double oneOverSpringPartitionUpperLimit = 1.0,
        bool allowSlipLinksToPassEachOther = false) const;

      /**
       * @brief Loop all slip-links and move them if appropriate to other
       * springs
       *
       * @param net
       * @param u
       * @param springPartitions
       * @param oneOverSpringPartitionUpperLimit
       */
      void moveSlipLinksToTheirBestBranch(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const double oneOverSpringPartitionUpperLimit,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true,
        const bool moveAttempt = false);

      /**
       * @brief Move a slip-link if appropriate to other springs
       *
       * @param net
       * @param u
       * @param springPartitions
       * @param oneOverSpringPartitionUpperLimit
       */
      void moveSlipLinkToItsBestBranch(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        size_t slipLinkIdx,
        const double oneOverSpringPartitionUpperLimit,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true,
        const bool moveAttempt = false);

      /**
       * @brief Loop all springs, swap slip-links on them if they are close
       * enough
       *
       * @param net
       */
      void swapSlipLinksInclXlinks(ForceBalanceNetwork& net,
                                   const Eigen::VectorXd& u,
                                   Eigen::VectorXd& springPartitions,
                                   double swappableCutoff,
                                   const bool respectLoops = true);

      /**
       * @brief Loop all springs, swap slip-links on them if they are close
       * enough
       *
       * @param net
       */
      void swapSlipLinks(ForceBalanceNetwork& net,
                         const Eigen::VectorXd& u,
                         Eigen::VectorXd& springPartitions,
                         double swappableCutoff);

      void swapSlipLinks(ForceBalanceNetwork& net,
                         const size_t partialSpringIdx);

      bool swapSlipLinkReversibly(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true,
        const bool moveAttempt = false);

      long int rotateSlipLinkAroundCrosslink(
        ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        double oneOverSpringPartitionUpperLimit = 1.0,
        const bool respectLoops = true);

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @param linkIdx the idx of the link to displace
       * @return double, the distance (squared norm) displaced
       */
      Eigen::VectorXd inspectDisplacementToMeanPositionUpdate(
        const size_t linkIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0) const
      {
        Eigen::VectorXd displacements = this->currentDisplacements;
        this->displaceToMeanPosition(this->initialConfig,
                                     displacements,
                                     this->currentSpringPartitionsVec,
                                     linkIdx,
                                     oneOverSpringPartitionUpperLimit);
        return displacements;
      };

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @param linkIdx the idx of the link to displace
       * @return double, the distance (squared norm) displaced
       */
      Eigen::VectorXd inspectLinkDisplacementToMeanPositionUpdate(
        const size_t linkIdx,
        double damping = 1.0) const
      {
        Eigen::VectorXd displacements = this->currentDisplacements;
        Eigen::VectorXd oneOverSpringPartitions =
          this->assembleOneOverSpringPartition(
            this->initialConfig, this->currentSpringPartitionsVec);
        Eigen::Array3i mask;
        mask << 3 * linkIdx, 3 * linkIdx + 1, 3 * linkIdx + 2;
        this->displaceLinksToMeanPosition(this->initialConfig,
                                          displacements,
                                          oneOverSpringPartitions,
                                          mask,
                                          damping);
        return displacements;
      };

      /**
       * @brief Adjust the two spring's box offsets to work best with the
       * specified slip-link
       *
       * @param net the network to adjust
       * @param slipLinkIdx the slip-link around which to adjust the two springs
       * @param spring1 one of the two partial spring idx
       * @param spring2 the partial spring idx of the other spring
       */
      void reAlignSlipLinkToImages(ForceBalanceNetwork& net,

                                   const Eigen::VectorXd& u,
                                   const size_t slipLinkIdx,
                                   const size_t spring1,
                                   const size_t spring2) const;

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param net the force balance network
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @param linkIdx the idx of the link to displace
       * @return double, the distance (squared norm) displaced
       */
      double displaceToMeanPosition(
        const ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        const size_t linkIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Translate the spring partition vector to its 3*size
       *
       * @param net
       * @param springPartitions0
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd assembleOneOverSpringPartition(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& springPartitions0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      double displaceLinksToMeanPosition(const ForceBalanceNetwork& net,
                                         Eigen::VectorXd& u,
                                         Eigen::VectorXd& springPartitions0,
                                         double damping = 0.5) const;

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param net the force balance network
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @return double, the distance (squared norm) displaced
       */
      double displaceLinksToMeanPosition(
        const ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        const Eigen::VectorXd& oneOverSpringPartitions,
        const Eigen::ArrayXi& resultingCoordinateIndexMask,
        const double damping) const;

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param net the force balance network
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @return double, the distance (squared norm) displaced
       */
      double displaceLinksToMeanPosition(
        const ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        const Eigen::VectorXd& oneOverSpringPartitions,
        const Eigen::ArrayXi& involvedSpringPartCoordinateIndexMask,
        const Eigen::ArrayXi& resultingCoordinateIndexMask,
        const double damping = 0.5) const;

      double getDisplacementResidualNorm(
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      double getDisplacementResidualNormFor(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        const double oneOverSpringPartitionUpperLimit) const;

      double getDisplacementResidualNormFor(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& oneOverSpringPartitions) const;

      /**
       * @brief Get the Link Indices of all neighbours of a specified link
       *
       * @param net
       * @param linkIdx
       * @return std::vector<size_t>
       */
      std::vector<size_t> getNeighbourLinkIndices(
        const ForceBalanceNetwork& net,
        const size_t linkIdx)
      {
        std::vector<size_t> results;
        results.reserve(4);
        for (size_t springIdx : net.springIndicesOfLinks[linkIdx]) {
          for (size_t partialSpringIdx :
               net.localToGlobalSpringIndex[springIdx]) {
            if (net.springPartIndexA[partialSpringIdx] == linkIdx) {
              results.push_back(net.springPartIndexB[partialSpringIdx]);
            } //
            else if (net.springPartIndexB[partialSpringIdx] == linkIdx) {
              results.push_back(net.springPartIndexA[partialSpringIdx]);
            }
          }
        }
        return results;
      }

      void writeRestartFile(std::string& file) override
      {
        throw std::runtime_error("Restart not supported yet");
      }

      double getTimestep() override
      {
        return 1.; // this->dt;
      }

      double getCurrentTime(double currentStep) override
      {
        return this->nrOfStepsDone;
      }

      Eigen::Matrix3d getStressTensor() override
      {
        return this->getStressTensor(-1.0);
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
      bool ConvertNetwork(ForceBalanceNetwork& net,
                          const int crosslinkerType,
                          bool remove2functionalCrosslinkers = false,
                          bool removeDanglingChains = false);

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
       * @brief Evaluate the pressure of the network at specific displacements
       *
       * @param net the network to evaluate the pressure for
       * @param u the displacements
       * @return double
       */
      double evaluatePressure(const ForceBalanceNetwork& net,
                              const Eigen::VectorXd& u,
                              const Eigen::VectorXd& springPartitions) const
      {
        auto stressTensor =
          this->evaluateStressTensor(net, u, springPartitions);
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
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensorLinkBased(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const bool xlinksOnly = false) const;

      /**
       * @brief Compute the stress tensor
       *
       * @param net
       * @param u
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */

      Eigen::Matrix3d evaluateStressTensorForLinks(
        const std::vector<size_t> linkIndices,
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Compute the stress tensor
       *
       * @param net
       * @param u
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Compute the force acting on a slip- or cross-link
       *
       * TODO: use "global" partial distances
       *
       * @param linkIdx
       * @param net
       * @param u
       * @param springPartitions
       * @param kappa0
       * @param minCutoff
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluateForceOnLink(
        const size_t linkIdx,
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        Eigen::VectorXi& debugNrSpringsVisited,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Compute the stress acting on a slip- or cross-link
       *
       * TODO: use "global" partial distances
       *
       * @param linkIdx
       * @param net
       * @param u
       * @param springPartitions
       * @param kappa0
       * @param minCutoff
       * @return Eigen::Vector3d
       */
      Eigen::Matrix3d evaluateStressOnLink(
        const size_t linkIdx,
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        Eigen::VectorXi& debugNrSpringsVisited,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Count how many of the springs are active (length > tolerance)
       *
       * @param springDistances
       * @param tolerance
       * @return int
       */
      int countNrOfActiveSprings(const Eigen::VectorXd& springDistances,
                                 const double tolerance = 0.01) const
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
       * @return Eigen::ArrayXb
       */
      Eigen::ArrayXb findActiveSprings(const Eigen::VectorXd& springDistances,
                                       const double tolerance = 0.01) const
      {
        Eigen::ArrayXb result =
          Eigen::ArrayXb::Constant(springDistances.size() / 3, false);
        for (size_t i = 0; i < springDistances.size() / 3; ++i) {
          result[i] =
            springDistances.segment(3 * i, 3).squaredNorm() > tolerance;
        }
        return result;
      }

      /**
       * @brief Sets the spring into the network
       *
       * @param net
       * @param springIdx
       * @param linkFrom
       * @param linkTo
       */
      void registerPartialSpring(ForceBalanceNetwork& net,
                                 const size_t springIdx,
                                 const size_t linkFrom,
                                 const size_t linkTo)
      {
        net.springPartIndexA[springIdx] = linkFrom;
        net.springPartIndexB[springIdx] = linkTo;
        for (size_t i = 0; i < 3; ++i) {
          net.springPartCoordinateIndexA[3 * springIdx + i] = linkFrom * 3 + i;
          net.springPartCoordinateIndexB[3 * springIdx + i] = linkTo * 3 + i;
        }
      }

      /**
       * @brief Set the Link Properties From two Atom objects
       *
       * @param net
       * @param linkIdx
       * @param atom1
       * @param atom2
       * @param atomType
       */
      void setLinkPropertiesFromAtoms(
        ForceBalanceNetwork& net,
        const size_t linkIdx,
        const pylimer_tools::entities::Atom& atom1,
        const pylimer_tools::entities::Atom& atom2,
        int atomType)
      {
        assert(linkIdx < net.nrOfLinks);
        assert(atom1.getType() != this->crosslinkerType &&
               atom2.getType() != this->crosslinkerType);
        if (atom1.getId() > atom2.getId()) {
          // make sure a second call to this function would result in same
          // result
          setLinkPropertiesFromAtoms(net, linkIdx, atom2, atom1, atomType);
          return;
        }

        Eigen::Vector3d coords = atom1.getCoordinates();
        Eigen::Vector3d dist = atom2.getCoordinates() - coords;
        this->universe.getBox().handlePBC(dist);
        coords += 0.5 * dist;

        net.coordinates.segment(3 * linkIdx, 3) = coords;
        net.linkIsSliplink[linkIdx] = atomType != this->crosslinkerType;
      }

      /**
       * @brief Set the Link Properties From an Atom object
       *
       * @param net
       * @param linkIdx
       * @param atom
       * @param atomType
       */
      void setLinkPropertiesFromAtom(ForceBalanceNetwork& net,
                                     const size_t linkIdx,
                                     const pylimer_tools::entities::Atom& atom,
                                     int atomType = -1)
      {
        if (atomType == -1) {
          atomType = atom.getType();
        }

        net.coordinates.segment(3 * linkIdx, 3) = atom.getCoordinates();
        net.linkIsSliplink[linkIdx] = atomType != this->crosslinkerType;
        if (!net.linkIsSliplink[linkIdx]) {
          net.oldAtomIds[linkIdx] = atom.getId();
        }
      }

      /**
       * @brief Use an existing chain to set the relevant edge properties
       *
       * @param chain
       * @param atom1Idx
       * @param atom2Idx
       * @param edgeId
       * @param chainIdx
       */
      void setPartialSpringPropertiesBasedOnChain(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& partitions,
        const pylimer_tools::entities::Molecule& chain,
        const size_t atom1Idx,
        const size_t atom2Idx,
        const size_t springIdx,
        const size_t partialSpringIdx)
      {
        assert(atom1Idx < atom2Idx);
        if (chain.getType() ==
            pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
          assert(APPROX_WITHIN(atom1Idx, 0, chain.getLength() - 1, 1e-2));
          assert(APPROX_WITHIN(atom2Idx, 1, chain.getLength(), 1e-2));
        } else {
          assert(APPROX_WITHIN(atom1Idx, 0, chain.getLength() - 2, 1e-2));
          assert(APPROX_WITHIN(atom2Idx, 1, chain.getLength() - 1, 1e-2));
        }

        // set the partition, remember the partial spring mapping
        size_t from = net.springPartIndexA[partialSpringIdx];
        size_t to = net.springPartIndexB[partialSpringIdx];
        std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
          chain.getAtomsLinedUp(crosslinkerType, false, true);
        partitions[partialSpringIdx] =
          static_cast<double>(atom2Idx - atom1Idx) /
          static_cast<double>(chain.getNrOfBonds());
        net.partialToFullSpringIndex[partialSpringIdx] = springIdx;
        net.partialSpringIsPartial[partialSpringIdx] = !(
          atom1Idx == 0 &&
          atom2Idx == chain.getLength() -
                        (chain.getType() ==
                             pylimer_tools::entities::MoleculeType::PRIMARY_LOOP
                           ? 0
                           : 1));

        // remember the ids -> springs
        net.oldAtomIdToSpringIndex[linedUpAtoms[atom1Idx].getId()] = springIdx;
        net.oldAtomIdToSpringIndex[linedUpAtoms[atom2Idx].getId()] = springIdx;

        // use the actual position of the vertices!
        Eigen::Vector3d expectedDistance =
          chain.getOverallBondSumFromTo(linedUpAtoms[atom1Idx].getId(),
                                        linedUpAtoms[atom2Idx].getId(),
                                        crosslinkerType);
        Eigen::Vector3d additionalDistance1 =
          linedUpAtoms[atom1Idx].getCoordinates() -
          net.coordinates.segment(3 * from, 3);
        this->universe.getBox().handlePBC(additionalDistance1);
        Eigen::Vector3d additionalDistance2 =
          net.coordinates.segment(3 * to, 3) -
          linedUpAtoms[atom2Idx].getCoordinates();
        this->universe.getBox().handlePBC(additionalDistance2);
        expectedDistance += additionalDistance1 + additionalDistance2;
        if (this->is2D) {
          expectedDistance[2] = 0.0;
        }

        net.springPartBoxOffset.segment(3 * partialSpringIdx, 3) =
          Eigen::Vector3d::Zero();
        Eigen::Vector3d actualDistance = this->evaluatePartialSpringDistance(
          net,
          Eigen::VectorXd::Zero(net.coordinates.size()),
          partialSpringIdx,
          this->is2D);
        net.springPartBoxOffset.segment(3 * partialSpringIdx, 3) =
          expectedDistance - actualDistance;
        assert(this->universe.getBox().isValidOffset(expectedDistance -
                                                     actualDistance));
#ifndef NDEBUG
        Eigen::Vector3d newActualDistance = this->evaluatePartialSpringDistance(
          net,
          Eigen::VectorXd::Zero(net.coordinates.size()),
          partialSpringIdx,
          this->is2D);
        assert(newActualDistance.isApprox(expectedDistance));
#endif
      }

      /**
       * @brief Find the other partial spring connected
       *
       * @param net
       * @param partialSpringIdx
       * @param aroundLinkIdx
       * @param notPartialSpringIdx useful for (double) secondary loops, to
       * distinguish them
       * @return size_t
       */
      size_t getOtherRailPartialSpringIdx(
        const ForceBalanceNetwork& net,
        const size_t partialSpringIdx,
        const size_t aroundLinkIdx,
        const long int notPartialSpringIdx = -1) const
      {
        assert(net.linkIsSliplink[aroundLinkIdx]);
        assert(this->isPartOfSpring(net, aroundLinkIdx, partialSpringIdx));
        size_t fullSpringIdx = net.partialToFullSpringIndex[partialSpringIdx];

        for (size_t i = 0;
             i < net.localToGlobalSpringIndex[fullSpringIdx].size();
             ++i) {
          if (net.localToGlobalSpringIndex[fullSpringIdx][i] ==
              partialSpringIdx) {
            if (i == 0) {
              return net.localToGlobalSpringIndex[fullSpringIdx][i + 1];
            }
            if (i >= (net.localToGlobalSpringIndex[fullSpringIdx].size() - 1)) {
              return net.localToGlobalSpringIndex[fullSpringIdx][i - 1];
            }
            size_t candidate1 =
              net.localToGlobalSpringIndex[fullSpringIdx][i - 1];
            size_t candidate2 =
              net.localToGlobalSpringIndex[fullSpringIdx][i + 1];
            size_t result = candidate1;
            if (candidate1 == notPartialSpringIdx &&
                this->isPartOfSpring(net, aroundLinkIdx, candidate2)) {
              result = candidate2;
            } else if (candidate2 == notPartialSpringIdx &&
                       this->isPartOfSpring(net, aroundLinkIdx, candidate1)) {
              result = candidate1;
            } else if (this->isPartOfSpring(net, aroundLinkIdx, candidate1)) {
              result = candidate1;
            } else {
              result = candidate2;
            }
            assert(this->isPartOfSpring(net, aroundLinkIdx, result));
            return result;
          }
        }

        // alternative:
        /**
         *
        for (size_t inSpringIdx = 1;
           inSpringIdx < net.linkIndicesOfSprings[springIdx].size() - 1;
           ++inSpringIdx) {
        if (net.localToGlobalSpringIndex[springIdx][inSpringIdx] ==
            partialSpringIdx) {
          if (net.linkIndicesOfSprings[springIdx][inSpringIdx] == linkIdx1 &&
              net.linkIndicesOfSprings[springIdx][inSpringIdx + 1] ==
                linkIdx2) {
            RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 == -1,
                            "Expect to find sequence of links only once.");
            otherPartialOfLinkIdx1 =
              net.localToGlobalSpringIndex[springIdx][inSpringIdx - 1];
            otherPartialOfLinkIdx2 =
              net.localToGlobalSpringIndex[springIdx][inSpringIdx + 1];
            firstPositionInSpring = inSpringIdx;
          }
          // else
          if (net.linkIndicesOfSprings[springIdx][inSpringIdx] == linkIdx2 &&
              net.linkIndicesOfSprings[springIdx][inSpringIdx + 1] ==
                linkIdx1) {
            RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 == -1,
                            "Expect to find sequence of links only once.");
            otherPartialOfLinkIdx2 =
              net.localToGlobalSpringIndex[springIdx][inSpringIdx - 1];
            otherPartialOfLinkIdx1 =
              net.localToGlobalSpringIndex[springIdx][inSpringIdx + 1];
            firstPositionInSpring = inSpringIdx;
          }
        }
      }
         */

        throw std::runtime_error("Did not find requested partial spring in "
                                 "full-spring, cannot determine rail.");
      }

      size_t getOtherEnd(const ForceBalanceNetwork& net,
                         size_t partialSpringIdx,
                         size_t linkIdx) const
      {
        assert(this->isPartOfSpring(net, linkIdx, partialSpringIdx));
        return net.springPartIndexA[partialSpringIdx] == linkIdx
                 ? net.springPartIndexB[partialSpringIdx]
                 : net.springPartIndexA[partialSpringIdx];
      }

      bool isPartOfSpring(const ForceBalanceNetwork& net,
                          size_t linkIdx,
                          size_t partialSpringIdx) const
      {
        return (net.springPartIndexA[partialSpringIdx] == linkIdx) ||
               (net.springPartIndexB[partialSpringIdx] == linkIdx);
      }

      bool isPrimaryLoop(const ForceBalanceNetwork& net,
                         size_t partialSpringIdx) const
      {
        return (net.springPartIndexA[partialSpringIdx] ==
                net.springPartIndexB[partialSpringIdx]);
      }

      double getDenominatorOfPartialSpring(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      long int moveSlipLinkFromRailToRail(
        ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const size_t sourcePartialSpringIdx,
        size_t targetPartialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        INVALIDARG_EXP_IFN(sourcePartialSpringIdx != targetPartialSpringIdx,
                           "Source and target must be different when moving "
                           "from one cross-link branch to another.");

        size_t involvedSlipLink =
          net.linkIsSliplink[net.springPartIndexA[sourcePartialSpringIdx]]
            ? net.springPartIndexA[sourcePartialSpringIdx]
            : net.springPartIndexB[sourcePartialSpringIdx];
        assert(net.linkIsSliplink[involvedSlipLink]);
        size_t involvedCrossLink =
          net.linkIsSliplink[net.springPartIndexA[sourcePartialSpringIdx]]
            ? net.springPartIndexB[sourcePartialSpringIdx]
            : net.springPartIndexA[sourcePartialSpringIdx];
        assert(!net.linkIsSliplink[involvedCrossLink]);
        assert(
          (net.springPartIndexA[targetPartialSpringIdx] == involvedCrossLink) ||
          (net.springPartIndexB[targetPartialSpringIdx] == involvedCrossLink));

        size_t otherInvolvedPartialSpring = this->getOtherRailPartialSpringIdx(
          net, sourcePartialSpringIdx, involvedSlipLink);

        // check whether there is space on the target spring, at all.
        size_t targetFullSpringIdx =
          net.partialToFullSpringIndex[targetPartialSpringIdx];
        const double minAlpha =
          (oneOverSpringPartitionUpperLimit > 0.)
            ? 1. / (net.springsContourLength[targetFullSpringIdx] -
                    1. / oneOverSpringPartitionUpperLimit)
            : 1e-9;
        if (minAlpha *
              (net.localToGlobalSpringIndex[targetFullSpringIdx].size() + 1.) >
            1.) {
          //
          return -1;
        }

        // validation: check distances
        Eigen::Vector3d distanceBefore =
          this->evaluatePartialSpringDistanceTo(
            net, u, otherInvolvedPartialSpring, involvedSlipLink, this->is2D) +
          this->evaluatePartialSpringDistanceTo(
            net, u, sourcePartialSpringIdx, involvedCrossLink, this->is2D) +
          this->evaluatePartialSpringDistanceFrom(
            net, u, targetPartialSpringIdx, involvedCrossLink, this->is2D);

        // remove the slip-link from one branch of the x-link
        // but skip resizing the Eigen structures, since the additional rows are
        // still needed
        this->mergePartialSprings(net,
                                  u,
                                  springPartitions,
                                  sourcePartialSpringIdx,
                                  otherInvolvedPartialSpring,
                                  involvedSlipLink,
                                  true);
        // this->validateNetwork(net, u, springPartitions);
        // ... and add it to another
        // assert(currentPartialSpringTargetIdx >= 0);
        // std::cout << "Handling moving link " << involvedSlipLink
        //           << " around cross-link " << involvedCrosslink
        //           << " from partial " << partialSpringIdx << " to "
        //           << targetPartialSpringIdx << std::endl;

        if (targetPartialSpringIdx > sourcePartialSpringIdx) {
          targetPartialSpringIdx -= 1;
        }
        if (otherInvolvedPartialSpring > sourcePartialSpringIdx) {
          otherInvolvedPartialSpring -= 1;
        }

        assert(
          this->isPartOfSpring(net, involvedCrossLink, targetPartialSpringIdx));

        size_t newPartialSpringIdx =
          this->addSlipLinkToPartialSpring(net,
                                           u,
                                           springPartitions,
                                           targetPartialSpringIdx,
                                           involvedSlipLink,
                                           oneOverSpringPartitionUpperLimit);

        // finally, return the idx of the new partial spring
        long int resultingPartialSpringIdx = 0;
        long int remainingPartialSpringIdx = 0;
        if ((net.springPartIndexA[targetPartialSpringIdx] ==
               involvedCrossLink &&
             net.springPartIndexB[targetPartialSpringIdx] ==
               involvedSlipLink) ||
            (net.springPartIndexB[targetPartialSpringIdx] ==
               involvedCrossLink &&
             net.springPartIndexA[targetPartialSpringIdx] ==
               involvedSlipLink)) {
          resultingPartialSpringIdx = targetPartialSpringIdx;
          remainingPartialSpringIdx = newPartialSpringIdx;
        } else {
          RUNTIME_EXP_IFN(
            (net.springPartIndexA[newPartialSpringIdx] == involvedCrossLink &&
             net.springPartIndexB[newPartialSpringIdx] == involvedSlipLink) ||
              (net.springPartIndexB[newPartialSpringIdx] == involvedCrossLink &&
               net.springPartIndexA[newPartialSpringIdx] == involvedSlipLink),
            "Expected to find cross- and slip-link at either partial spring, "
            "but did not.");
          resultingPartialSpringIdx = newPartialSpringIdx;
          remainingPartialSpringIdx = targetPartialSpringIdx;
        }

        // validation: check distances
        if (!this->assumeBoxLargeEnough &&
            !this->isPrimaryLoop(net, resultingPartialSpringIdx) &&
            !this->isPrimaryLoop(net, remainingPartialSpringIdx)) {
          Eigen::Vector3d distanceAfter =
            this->evaluatePartialSpringDistanceTo(net,
                                                  u,
                                                  otherInvolvedPartialSpring,
                                                  involvedCrossLink,
                                                  this->is2D) +
            this->evaluatePartialSpringDistanceTo(
              net, u, resultingPartialSpringIdx, involvedSlipLink, this->is2D) +
            this->evaluatePartialSpringDistanceFrom(
              net, u, remainingPartialSpringIdx, involvedSlipLink, this->is2D);
          assert(distanceAfter.isApprox(distanceBefore));
        }
        return resultingPartialSpringIdx;
      };

      bool swapSlipLinksReversibly(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0);

      bool swapSlipLinkWithXlinkReversibly(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const bool respectLoops = true);

    private:
      pylimer_tools::entities::Universe universe;
      pylimer_tools::entities::Box box;
      bool is2D = false;
      bool assumeBoxLargeEnough = true;
      double kappa = 1.0;
      bool simulationHasRun = false;
      int stepOutputFrequency = 0;
      int defaultNrOfChains = 0;
      double defaultR0Squared = 0.0;
      std::string stepOutputFile;
      bool outputEndNodes = false;
      std::string endNodesFile;
      ForceBalanceNetwork initialConfig;
      Eigen::VectorXd currentDisplacements;
      Eigen::VectorXd currentSpringDistances;
      Eigen::VectorXd currentPartialSpringDistances;
      Eigen::VectorXd
        currentSpringPartitionsVec; /* gives the parametrisation of N */
      int crosslinkerType = 2;
      int sliplinkType = 3;
      int nrOfStepsDone = 0;
      ExitReason exitReason = ExitReason::UNSET;
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
