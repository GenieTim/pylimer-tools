#pragma once

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../sim/MEHPUtilityStructures.h"
#include "../sim/OutputSupportingSimulation.h"
#include "../topo/EntanglementDetector.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <limits>
#include <map>
#include <nlopt.hpp>
#include <random>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>
#ifdef CEREALIZABLE
#include "../utils/CerealUtils.h"
#include <cereal/access.hpp>
#endif

namespace pylimer_tools::sim::mehp {

class MEHPForceBalance final
  : public pylimer_tools::sim::OutputSupportingSimulation
{
private:
#ifdef CEREALIZABLE
  MEHPForceBalance() {}; // not exposed to users, only used by Cereal

  friend class cereal::access;
#endif

  // member properties
  pylimer_tools::entities::Universe universe;
  pylimer_tools::entities::Box box;
  // state
  ExitReason exitReason = ExitReason::UNSET;
  bool simulationHasRun = false;
  ForceBalanceNetwork initialConfig;
  Eigen::VectorXd currentDisplacements;
  Eigen::VectorXd currentSpringVectors;
  Eigen::VectorXd currentPartialSpringVectors;
  Eigen::VectorXd
    currentSpringPartitionsVec; /* gives the parametrisation of N */
  // configuration
  bool is2D = false;
  bool assumeBoxLargeEnough = true;
  double kappa = 1.0;
  int crossLinkerType = 2;
  int sliplinkType = 3;
  int nrOfStepsDone = 0;
  int simplificationFrequency = 10;
  int entanglementType = -999;
  double defaultBondLength = 0.0;
  double springBreakingLength = -1.;

public:
  MEHPForceBalance(const pylimer_tools::entities::Universe& u,
                   const int crossLinkerType = 2,
                   const bool is2D = false,
                   const bool remove2functionalCrosslinkers = false,
                   const bool removeDanglingChains = false)
    : universe(u)
  {
    this->crossLinkerType = crossLinkerType;
    this->box = u.getBox();
    // interpret network already to be able to give early results
    ForceBalanceNetwork net;
    RUNTIME_EXP_IFN(ConvertNetwork(net,
                                   crossLinkerType,
                                   remove2functionalCrosslinkers,
                                   removeDanglingChains),
                    "Failed to convert network.");
    this->initialConfig = net;
    this->is2D = is2D;
    this->currentDisplacements = Eigen::VectorXd::Zero(net.coordinates.size());
    this->completeInitialization();
  };

  MEHPForceBalance(const ForceBalanceNetwork& net,
                   const Eigen::VectorXd& springPartitions,
                   const bool is2D = false)
  {
    this->is2D = is2D;
    this->initialConfig = net;
    this->currentSpringPartitionsVec = springPartitions;
    this->currentDisplacements = Eigen::VectorXd::Zero(net.coordinates.size());
    this->box = pylimer_tools::entities::Box(net.L[0], net.L[1], net.L[2]);
    this->completeInitialization();
  }

  /**
   * @brief Instantiate this simulator with randomly chosen slip-links.
   *
   * @param universe
   * @param entanglements
   * @param crossLinkerType
   * @param is2D
   * @return MEHPForceBalance
   */
  static MEHPForceBalance constructWithSlipLinks(
    const pylimer_tools::entities::Universe& universe,
    pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
      entanglements,
    int crossLinkerType = 2,
    bool is2D = false)
  {
    pylimer_tools::entities::Universe emptyUniverse =
      pylimer_tools::entities::Universe(universe.getBox());
    MEHPForceBalance fb =
      MEHPForceBalance(emptyUniverse, crossLinkerType, is2D, false, false);
    fb.configAssumeBoxLargeEnough(false);
    fb.universe = universe;

    std::vector<std::pair<size_t, size_t>> pairsOfAtoms =
      entanglements.pairsOfAtoms;
    std::vector<long int> pairOfAtom = entanglements.pairOfAtom;

    std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
      universe.getChainsWithCrosslinker(crossLinkerType);

    // add ends of chains
    std::vector<long int> vertexIdToLinkIdx =
      pylimer_tools::utils::initializeWithValue<long int>(
        fb.universe.getNrOfAtoms(), -1);
    size_t currentVertexId = 0;
    size_t numUseableChains = 0;
    for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
      pylimer_tools::entities::Molecule chain = crossLinkerChains[i];
      if (chain.getNrOfBonds() < 1) {
        continue;
      }
      std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
        chain.getAtomsLinedUp(crossLinkerType, false, true);
      if (vertexIdToLinkIdx[fb.universe.getIdxByAtomId(
            linedUpAtoms[0].getId())] == -1) {
        vertexIdToLinkIdx[fb.universe.getIdxByAtomId(linedUpAtoms[0].getId())] =
          currentVertexId;
        currentVertexId += 1;
      }
      if (vertexIdToLinkIdx[fb.universe.getIdxByAtomId(
            pylimer_tools::utils::last(linedUpAtoms).getId())] == -1) {
        vertexIdToLinkIdx[fb.universe.getIdxByAtomId(
          pylimer_tools::utils::last(linedUpAtoms).getId())] = currentVertexId;
        currentVertexId += 1;
      }
      numUseableChains += 1;
    }

    // resize
    // links
    fb.initialConfig.nrOfNodes = currentVertexId;
    fb.initialConfig.nrOfLinks = currentVertexId + pairsOfAtoms.size();
    fb.initialConfig.oldAtomIds.resize(fb.initialConfig.nrOfNodes);
    fb.initialConfig.oldAtomTypes.resize(fb.initialConfig.nrOfNodes);
    fb.currentDisplacements.resize(3 * fb.initialConfig.nrOfLinks);
    fb.currentDisplacements.setZero();
    fb.initialConfig.coordinates.conservativeResize(3 *
                                                    fb.initialConfig.nrOfLinks);
    fb.initialConfig.linkIsSliplink.conservativeResize(
      fb.initialConfig.nrOfLinks);
    fb.initialConfig.springIndicesOfLinks =
      pylimer_tools::utils::initializeWithValue(fb.initialConfig.nrOfLinks,
                                                std::vector<size_t>());
    fb.initialConfig.nrOfCrosslinkSwapsEndured.conservativeResize(
      fb.initialConfig.nrOfLinks - fb.initialConfig.nrOfNodes);
    fb.initialConfig.nrOfCrosslinkSwapsEndured.setZero();

    // springs
    fb.initialConfig.nrOfSprings = numUseableChains;
    fb.initialConfig.springIndexA.conservativeResize(numUseableChains);
    fb.initialConfig.springIndexB.conservativeResize(numUseableChains);
    fb.initialConfig.springCoordinateIndexA.conservativeResize(
      3 * numUseableChains);
    fb.initialConfig.springCoordinateIndexB.conservativeResize(
      3 * numUseableChains);
    fb.initialConfig.springIsActive.conservativeResize(numUseableChains);
    fb.initialConfig.springsContourLength.conservativeResize(numUseableChains);
    fb.initialConfig.springsType.conservativeResize(numUseableChains);
    fb.initialConfig.linkIndicesOfSprings =
      pylimer_tools::utils::initializeWithValue(numUseableChains,
                                                std::vector<size_t>());
    fb.initialConfig.localToGlobalSpringIndex =
      pylimer_tools::utils::initializeWithValue(numUseableChains,
                                                std::vector<size_t>());

    // partial springs
    // we don't know the actual number (yet), but we can over-estimate
    // pretty well, such that we only need to reduce afterwards
    size_t numPartialSpringsEstimate =
      numUseableChains + 2 * pairsOfAtoms.size();
    fb.initialConfig.nrOfPartialSprings = numPartialSpringsEstimate;
    fb.currentSpringPartitionsVec.resize(numPartialSpringsEstimate);
    fb.initialConfig.springPartBoxOffset.conservativeResize(
      3 * numPartialSpringsEstimate);
    fb.initialConfig.springPartCoordinateIndexA.conservativeResize(
      3 * numPartialSpringsEstimate);
    fb.initialConfig.springPartCoordinateIndexB.conservativeResize(
      3 * numPartialSpringsEstimate);
    fb.initialConfig.springPartIndexA.conservativeResize(
      numPartialSpringsEstimate);
    fb.initialConfig.springPartIndexB.conservativeResize(
      numPartialSpringsEstimate);
    fb.initialConfig.partialToFullSpringIndex.conservativeResize(
      numPartialSpringsEstimate);
    fb.initialConfig.partialSpringIsPartial.conservativeResize(
      numPartialSpringsEstimate);

    size_t springIdx = 0;
    size_t partialSpringIdx = 0;
    for (size_t chainIdx = 0; chainIdx < crossLinkerChains.size(); ++chainIdx) {
      pylimer_tools::entities::Molecule chain = crossLinkerChains[chainIdx];
      if (chain.getLength() < 2) {
        continue;
      }

      Eigen::Vector3d overallDistance =
        chain.getOverallBondSum(fb.crossLinkerType);
      fb.initialConfig.springToMoleculeIds.push_back(chainIdx);
      std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
        chain.getAtomsLinedUp(crossLinkerType, false, true);
      size_t previousIdx = 0;

      size_t previousLinkIdx =
        vertexIdToLinkIdx[fb.universe.getIdxByAtomId(linedUpAtoms[0].getId())];
      fb.setLinkPropertiesFromAtom(
        fb.initialConfig, previousLinkIdx, linedUpAtoms[0], fb.crossLinkerType);
      fb.initialConfig.linkIndicesOfSprings[springIdx].push_back(
        previousLinkIdx);
      pylimer_tools::utils::addIfNotContained(
        fb.initialConfig.springIndicesOfLinks[previousLinkIdx], springIdx);
      pylimer_tools::entities::Atom lastAtom =
        pylimer_tools::utils::last(linedUpAtoms);
      size_t lastLinkIdx =
        vertexIdToLinkIdx[fb.universe.getIdxByAtomId(lastAtom.getId())];
      fb.setLinkPropertiesFromAtom(
        fb.initialConfig, lastLinkIdx, lastAtom, fb.crossLinkerType);
      fb.initialConfig.springIndexA[springIdx] = previousLinkIdx;
      assert(linedUpAtoms.size() == chain.getLength() ||
             linedUpAtoms.size() == chain.getLength() + 1);
      if (pairsOfAtoms.size() > 0) {
        for (size_t i = 1; i < linedUpAtoms.size() - 1; i++) {
          pylimer_tools::entities::Atom a = linedUpAtoms[i];
          long int pairIdx = pairOfAtom[universe.getIdxByAtomId(a.getId())];
          if (pairIdx != -1) {
            size_t thisLinkIdx = currentVertexId + pairIdx;
            fb.initialConfig.linkIndicesOfSprings[springIdx].push_back(
              thisLinkIdx);
            pylimer_tools::utils::addIfNotContained(
              fb.initialConfig.springIndicesOfLinks[thisLinkIdx], springIdx);
            // set the mean x,y,z of the two involved atoms
            pylimer_tools::entities::Atom a1 =
              universe.getAtom(pairsOfAtoms[pairIdx].first);
            pylimer_tools::entities::Atom a2 =
              universe.getAtom(pairsOfAtoms[pairIdx].second);
            RUNTIME_EXP_IFN(
              pairOfAtom[universe.getIdxByAtomId(a1.getId())] == pairIdx,
              "Atom 1 does not follow required atom pair pattern");
            RUNTIME_EXP_IFN(
              pairOfAtom[universe.getIdxByAtomId(a2.getId())] == pairIdx,
              "Atom 2 does not follow required atom pair pattern");
            fb.setLinkPropertiesFromAtoms(
              fb.initialConfig, thisLinkIdx, a1, a2, fb.sliplinkType);

            fb.registerPartialSpring(
              fb.initialConfig, partialSpringIdx, previousLinkIdx, thisLinkIdx);
            fb.setPartialSpringPropertiesBasedOnChain(
              fb.initialConfig,
              fb.currentSpringPartitionsVec,
              chain,
              previousIdx,
              i,
              springIdx,
              partialSpringIdx);
            fb.initialConfig.localToGlobalSpringIndex[springIdx].push_back(
              partialSpringIdx);
            //
            previousIdx = i;
            previousLinkIdx = thisLinkIdx;
            partialSpringIdx += 1;
          }
        }
      }

      // close the chain
      fb.registerPartialSpring(
        fb.initialConfig, partialSpringIdx, previousLinkIdx, lastLinkIdx);
      fb.initialConfig.linkIndicesOfSprings[springIdx].push_back(lastLinkIdx);
      pylimer_tools::utils::addIfNotContained(
        fb.initialConfig.springIndicesOfLinks[lastLinkIdx], springIdx);
      fb.setPartialSpringPropertiesBasedOnChain(fb.initialConfig,
                                                fb.currentSpringPartitionsVec,
                                                chain,
                                                previousIdx,
                                                linedUpAtoms.size() - 1,
                                                springIdx,
                                                partialSpringIdx);
      fb.initialConfig.localToGlobalSpringIndex[springIdx].push_back(
        partialSpringIdx);

      fb.initialConfig.springIndexB[springIdx] = lastLinkIdx;
      for (size_t dir = 0; dir < 3; ++dir) {
        fb.initialConfig.springCoordinateIndexA[3 * springIdx + dir] =
          3 * fb.initialConfig.springIndexA[springIdx] + dir;
        fb.initialConfig.springCoordinateIndexB[3 * springIdx + dir] =
          3 * fb.initialConfig.springIndexB[springIdx] + dir;
      }
      fb.initialConfig.springsContourLength[springIdx] = chain.getNrOfBonds();
      auto bondTypes = chain.getBonds()["bond_type"];
      fb.initialConfig.springsType[springIdx] = MEAN(bondTypes);

#ifndef NDEBUG
      Eigen::Vector3d overallDistanceNow = Eigen::Vector3d::Zero();
      for (size_t i = 0;
           i < fb.initialConfig.localToGlobalSpringIndex[springIdx].size();
           ++i) {
        size_t partialSpringIdx =
          fb.initialConfig.localToGlobalSpringIndex[springIdx][i];
        overallDistanceNow +=
          fb.evaluatePartialSpringDistance(fb.initialConfig,
                                           fb.currentDisplacements,
                                           partialSpringIdx,
                                           fb.is2D,
                                           false);
      }
      assert(pylimer_tools::utils::vector_approx_equal(overallDistanceNow,
                                                       overallDistance));
#endif

      partialSpringIdx += 1;
      springIdx += 1;
    }

    if (partialSpringIdx < numPartialSpringsEstimate) {
      // reduce sizes again
      fb.initialConfig.nrOfPartialSprings = partialSpringIdx;
      fb.currentSpringPartitionsVec.conservativeResize(partialSpringIdx);
      fb.initialConfig.springPartBoxOffset.conservativeResize(3 *
                                                              partialSpringIdx);
      fb.initialConfig.springPartCoordinateIndexA.conservativeResize(
        3 * partialSpringIdx);
      fb.initialConfig.springPartCoordinateIndexB.conservativeResize(
        3 * partialSpringIdx);
      fb.initialConfig.springPartIndexA.conservativeResize(partialSpringIdx);
      fb.initialConfig.springPartIndexB.conservativeResize(partialSpringIdx);
      fb.initialConfig.partialToFullSpringIndex.conservativeResize(
        partialSpringIdx);
      fb.initialConfig.partialSpringIsPartial.conservativeResize(
        partialSpringIdx);
    }

    fb.completeInitialization();

    return fb;
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
   * @param crossLinkerType
   * @param is2D
   * @return MEHPForceBalance
   */
  static MEHPForceBalance constructWithRandomSlipLinks(
    const pylimer_tools::entities::Universe& universe,
    const size_t nrOfSliplinksToSample,
    const double upperCutoff,
    const double lowerCutoff = 0.,
    const size_t minimumNrOfSliplinks = 0,
    const double sameStrandCutoff = 3,
    const std::string seed = "",
    const int crossLinkerType = 2,
    const bool is2D = false,
    const bool filterEntanglements = true)
  {
    // sample the "entanglements"
    const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
      entanglements =
        pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
          universe,
          nrOfSliplinksToSample,
          upperCutoff,
          lowerCutoff,
          minimumNrOfSliplinks,
          sameStrandCutoff,
          seed,
          crossLinkerType,
          true,
          filterEntanglements);

    RUNTIME_EXP_IFN(entanglements.pairsOfAtoms.size() >= minimumNrOfSliplinks ||
                      filterEntanglements,
                    "Minimum number of slip-links could not be sampled: got " +
                      std::to_string(entanglements.pairsOfAtoms.size()) +
                      " instead of " + std::to_string(minimumNrOfSliplinks) +
                      ".");

    return MEHPForceBalance::constructWithSlipLinks(
      universe, entanglements, crossLinkerType, is2D);
  }

  /**
   * @brief Finish initializing some member properties
   *
   */
  void completeInitialization()
  {
    this->currentSpringVectors = this->evaluateSpringVectors(
      this->initialConfig, this->currentDisplacements);
    this->currentPartialSpringVectors = this->evaluatePartialSpringVectors(
      this->initialConfig, this->currentDisplacements);
    this->defaultBondLength = universe.computeMeanBondLength();
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
  void runForceRelaxation(const long int maxNrOfSteps = 50000, // default: 10000
                          const double xtol = 1e-9,
                          const double initialResidualToUse = -1.0,
                          const StructureSimplificationMode simplificationMode =
                            StructureSimplificationMode::NO_SIMPLIFICATION,
                          const double inactiveRemovalCutoff = 1e-3,
                          const bool doInnerIterations = false,
                          const LinkSwappingMode allowSlipLinksToPassEachOther =
                            LinkSwappingMode::NO_SWAPPING,
                          const int swappingFrequency = 10,
                          const double oneOverSpringPartitionUpperLimit = 1.0,
                          const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
                          const bool disableSlipping = false)
  {
    this->runForceRelaxation(
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
      disableSlipping,
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
  void runForceRelaxation(long int maxNrOfSteps, // default: 10000
                          double xtol,
                          const double initialResidualToUse,
                          const StructureSimplificationMode simplificationMode,
                          const double inactiveRemovalCutoff,
                          bool doInnerIterations,
                          LinkSwappingMode allowSlipLinksToPassEachOther,
                          const int swappingFrequency,
                          const double oneOverSpringPartitionUpperLimit,
                          const int nrOfCrosslinkSwapsAllowedPerSliplink,
                          const bool disableSlipping,
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
    const double oneOverSpringPartitionUpperLimit = 1.0) const
  {
    // TODO: revise, hard!
    assert(involvedPartitions.size() == 4);
    const double firstMeanVal = 0.5 * (springPartitions[involvedPartitions[0]] +
                                       springPartitions[involvedPartitions[1]]);
    const double secondMeanVal =
      0.5 * (springPartitions[involvedPartitions[2]] +
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
    const Eigen::Vector4d partitionsBefore =
      springPartitions(involvedPartitions);
    springPartitions[involvedPartitions[0]] = firstMeanVal;
    springPartitions[involvedPartitions[1]] = firstMeanVal;
    springPartitions[involvedPartitions[2]] = secondMeanVal;
    springPartitions[involvedPartitions[3]] = secondMeanVal;
    // this->displaceToMeanPosition(
    //   this->initialConfig, displacements, springPartitions, link_idx);
    const double retVal =
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
   * @brief Remove crosslinkers, springs and associated slip-links with the
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
   * @brief Remove crosslinkers which do not have any springs with a certain
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
   * @brief Remove springs that exert a stress higher than
   * `this->springBreakingLength`
   *
   * @param net
   * @param displacements
   * @param springPartitions
   * @return size_t the number of springs broken
   */
  size_t breakTooLongSprings(ForceBalanceNetwork& net,
                             Eigen::VectorXd& displacements,
                             Eigen::VectorXd& springPartitions) const;

  size_t breakTooLongSprings()
  {
    return this->breakTooLongSprings(this->initialConfig,
                                     this->currentDisplacements,
                                     this->currentSpringPartitionsVec);
  }

  /**
   * @brief Remove double listed springs from crosslinkers
   *
   * @param net
   */
  void removeDuplicateListedSpringsFromLinks(ForceBalanceNetwork& net) const;

  void removeDuplicateListedSpringsFromLink(
    ForceBalanceNetwork& net,
    size_t linkIdx,
    bool allowOnEntanglement = false) const;

  /**
   *
   * @param net the network to rid of primary loops
   * @param displacements the corresponding displacements
   * @param springPartitions the corresponding spring partitions
   * @return the number of primary loops removed
   */
  size_t removePrimaryLoops(ForceBalanceNetwork& net,
                            Eigen::VectorXd& displacements,
                            Eigen::VectorXd& springPartitions) const;

  size_t removePrimaryLoops()
  {
    return this->removePrimaryLoops(this->initialConfig,
                                    this->currentDisplacements,
                                    this->currentSpringPartitionsVec);
  };

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
   * @brief break a spring, given its partial spring index
   *
   * @param net
   * @param displacements
   * @param springPartitions
   * @param partialSpringIdx
   */
  void breakPartialSpring(ForceBalanceNetwork& net,
                          Eigen::VectorXd& displacements,
                          Eigen::VectorXd& springPartitions,
                          const size_t partialSpringIdx) const;

  /**
   * @brief Remove a spring, but also all springs that are connected to it
   * and are connected via entanglement links.
   *
   * @param net
   * @param displacements
   * @param springPartitions
   * @param springIdx
   */
  void removeSpringFollowingEntanglementLinks(ForceBalanceNetwork& net,
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
  static void removeLink(ForceBalanceNetwork& net,
                         Eigen::VectorXd& displacements,
                         const size_t linkIdx);

  /**
   * @brief Merge two springs around a given crosslink
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
   * @brief Merge two springs around a given crosslink
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

  void relaxationLight(
    ForceBalanceNetwork& net,
    Eigen::VectorXd& springPartitions,
    Eigen::VectorXd& displacements,
    const size_t linkIdx,
    const double oneOverSpringPartitionUpperLimit = 1.0) const
  {
    Eigen::VectorXd oneOverSpringPartitions = Eigen::VectorXd::Zero(0);
    this->relaxationLight(net,
                          springPartitions,
                          oneOverSpringPartitions,
                          displacements,
                          linkIdx,
                          oneOverSpringPartitionUpperLimit);
  };

  void relaxationLight(
    ForceBalanceNetwork& net,
    Eigen::VectorXd& springPartitions,
    Eigen::VectorXd& oneOverSpringPartitions,
    Eigen::VectorXd& displacements,
    const size_t linkIdx,
    const double oneOverSpringPartitionUpperLimit = 1.0) const;

  /**
   * @brief Replace the two springs traversinga a two-functional
   * crosslinkers with a single spring
   *
   * @param net
   * @param displacements
   * @param springPartitions
   */
  size_t removeTwofunctionalCrosslinks(ForceBalanceNetwork& net,
                                       Eigen::VectorXd& displacements,
                                       Eigen::VectorXd& springPartitions) const;

  /**
   * @brief Add slip-links to this system
   *
   * @param nrOfSliplinksToSample the number of slip-links to sample
   * @param cutoff the distance (norm) between two atoms to be considered as a
   * slip-link
   * @param minimumNrOfSliplinks the minimum number of slip-links to sample
   * @param sameStrandCutoff the number of beads required between two atoms on
   * the same strand
   * @param excludeCrosslinks whether to exclude crosslink atoms when sampling.
   * @param seed
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
    this->box.adjustCoordinatesTo(this->initialConfig.springPartBoxOffset,
                                  newBox);
    this->box = newBox;
    this->universe.setBox(newBox, true);
    for (size_t i = 0; i < 3; ++i) {
      this->initialConfig.L[i] = this->box.getL(i);
      this->initialConfig.boxHalfs[i] = 0.5 * this->initialConfig.L[i];
    }
    this->initialConfig.vol = this->box.getVolume();
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
    const size_t link_idx,
    Eigen::VectorXd& displacements,
    Eigen::VectorXd& springPartitions,
    const long int innerMaxNrOfSteps = 500,
    const double innerAlphaTol = 1e-9,
    const long int innerMinNrOfSteps = 1,
    const double oneOverSpringPartitionUpperLimit = 1.0) const
  {
    size_t innerIterationsDone = 0;
    double displacementDone = 0.0;
    double rOverr0 = 0.0;
    double r2 = 0.0;
    double r02 = this->computePartitionUpdateZeroResidual(
      this->initialConfig,
      this->getSpringpartitionIndicesOfSliplink(this->initialConfig, link_idx),
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
   * @brief Get the universe consisting of crosslinkers only
   *
   * @param newCrosslinkerType the type to give the crosslinkers
   * @return pylimer_tools::entities::Universe
   */
  pylimer_tools::entities::Universe getCrosslinkerVerse() const;

  double getDefaultMeanBondLength() const { return this->defaultBondLength; }

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

  void configAssumeBoxLargeEnough(const bool assumption)
  {
    this->assumeBoxLargeEnough = assumption;

    this->currentSpringVectors = this->evaluateSpringVectors(
      this->initialConfig, this->currentDisplacements);
    this->currentPartialSpringVectors = this->evaluatePartialSpringVectors(
      this->initialConfig, this->currentDisplacements);
  }

  void configMeanBondLength(const double meanBondLength)
  {
    this->defaultBondLength = meanBondLength;
  }

  void configSpringConstant(const double kappa = 1.0) { this->kappa = kappa; }

  void configEntanglementType(const int newEntanglementType = -1)
  {
    this->entanglementType = newEntanglementType;
  }

  void configSimplificationFrequency(const int newRemovalFrequency = 10)
  {
    this->simplificationFrequency = newRemovalFrequency;
  }

  void configSpringBreakingDistance(const double newSpringBreakingForce = -1.)
  {
    this->springBreakingLength = newSpringBreakingForce;
  }

  /**
   * @brief Get the number of active nodes (incl. entanglement nodes, excl.
   * entanglement links)
   *
   * @param tolerance  the tolerance: springs under a certain length are
   * considered inactive
   * @return int
   */
  int getNrOfActiveNodes(const double tolerance = 1e-3) const
  {
    return this->findActiveNodes(tolerance).count();
  }

  /**
   * @brief Get the Soluble Weight Fraction
   *
   * @param tolerance
   * @return double
   */
  double getSolubleWeightFraction(const double tolerance = 1e-3)
  {
    return this->computeSolubleWeightFraction(&this->initialConfig,
                                              this->currentDisplacements,
                                              this->currentSpringPartitionsVec,
                                              tolerance);
  }

  /**
   * @brief Get the Dangling Weight Fraction
   *
   * @param tolerance
   * @return double
   */
  double getDanglingWeightFraction(const double tolerance = 1e-3)
  {
    return this->computeDanglingWeightFraction(&this->initialConfig,
                                               this->currentDisplacements,
                                               this->currentSpringPartitionsVec,
                                               tolerance);
  }

  /**
   * @brief Get the Weight Fraction of Active Springs (atoms)
   *
   * @param tolerance
   * @return double
   */
  double getActiveWeightFraction(const double tolerance = 1e-3)
  {
    return this->computeActiveWeightFraction(&this->initialConfig,
                                             this->currentDisplacements,
                                             this->currentSpringPartitionsVec,
                                             tolerance);
  }

  /**
   * @brief Count the number of atoms that are in any way connected to an
   * active spring
   *
   * @param tolerance
   * @return double
   */
  double countActiveClusteredAtoms(const double tolerance = 1e-3)
  {
    return this->countActiveClusteredAtoms(&this->initialConfig,
                                           this->currentDisplacements,
                                           this->currentSpringPartitionsVec,
                                           tolerance);
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
    double tolerance = 1e-3) const;

  /**
   * @brief Compute the weight fraction of non-active springs
   *
   * We go the full route via active and soluble in order to compensate for
   * removed springs and atoms
   *
   * @param net
   * @param tolerance
   * @return double
   */
  double computeDanglingWeightFraction(ForceBalanceNetwork* net,
                                       const Eigen::VectorXd& u,
                                       const Eigen::VectorXd& springPartitions,
                                       const double tolerance = 1e-3) const
  {
    const double activeWeightFraction =
      this->computeActiveWeightFraction(net, u, springPartitions, tolerance);
    RUNTIME_EXP_IFN(
      APPROX_WITHIN(activeWeightFraction, 0., 1., 1e-6),
      "Expect active weight fraction to be between 0 and 1, got " +
        std::to_string(activeWeightFraction) + ".");
    const double solubleWeightFraction =
      this->computeSolubleWeightFraction(net, u, springPartitions, tolerance);
    RUNTIME_EXP_IFN(
      APPROX_WITHIN(solubleWeightFraction, 0., 1., 1e-6),
      "Expect soluble weight fraction to be between 0 and 1, got " +
        std::to_string(solubleWeightFraction) + ".");
    RUNTIME_EXP_IFN(
      APPROX_WITHIN(activeWeightFraction + solubleWeightFraction, 0., 1., 1e-6),
      "Expect active and soluble weight fraction to add up to maximum 1, "
      "got " +
        std::to_string(activeWeightFraction + solubleWeightFraction) + ".");

    // finally, normalize by the number of atoms.
    // TODO: currently, the weight of the atoms is ignored
    return 1. - activeWeightFraction - solubleWeightFraction;
  }

  /**
   * @brief Compute the weight fraction of active springs
   *
   * @param net
   * @param tolerance
   * @return double
   */
  double computeActiveWeightFraction(ForceBalanceNetwork* net,
                                     const Eigen::VectorXd& u,
                                     const Eigen::VectorXd& springPartitions,
                                     const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(net->nrOfPartialSprings == springPartitions.size(),
                       "Spring partitions and network don't match");
    INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                       "Link displacements and network don't match");
    if (net->nrOfSprings < 1) {
      return 0.;
    }
    // find all active springs
    const Eigen::ArrayXb activeSprings =
      this->findActiveSprings(net, u, springPartitions, tolerance);
    const double nActiveSprings = activeSprings.count();
    if (nActiveSprings == 0) {
      return 0.;
    }

    // as of now, the springsContourLength is equal to the number of bonds
    // from crosslink to crosslink. therefore, the number of atoms of each
    // of these springs is one less
    Eigen::ArrayXd allActiveAtomsPerChains =
      activeSprings.cast<double>() * (net->springsContourLength.array() -
                                      Eigen::ArrayXd::Ones(net->nrOfSprings));

    // TODO: currently, the weight of the atoms is ignored
    // normalize by the number of atoms
    const double nActiveNodes = this->getNrOfActiveNodes(tolerance);
    return (allActiveAtomsPerChains.matrix().sum() + nActiveNodes) /
           (static_cast<double>(this->universe.getNrOfAtoms()));
  }

  /**
   * @brief Find whether springs and nodes are in any way connected to an
   * active spring
   *
   * @param net the network that includes the connectivity
   * @param springDistances the distances used to assert whether the springs
   * is active or not
   * @param tolerance the tolerance for considering springs as active
   * @return std::pair<Eigen::ArrayXb, Eigen::ArrayXb> indices of springs
   * (first) and links (second) connected in any way to active springs
   */
  std::pair<Eigen::ArrayXb, Eigen::ArrayXb> findClusteredToActive(
    const ForceBalanceNetwork* net,
    const Eigen::VectorXd& u,
    const Eigen::VectorXd& springPartitions,
    const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(springPartitions.size() == net->nrOfPartialSprings,
                       "Invalid sizes.");
    INVALIDARG_EXP_IFN(u.size() == net->nrOfLinks * 3, "Invalid sizes.");

    // find all active springs
    Eigen::ArrayXb springIsActive =
      this->findActiveSprings(net, u, springPartitions, tolerance);

    // then, iteratively walk along the springs to mark those as "active"
    // that are connected to active springs
    bool hadChanged = true;
    Eigen::ArrayXb linkIsActive =
      Eigen::ArrayXb::Constant(net->nrOfNodes, false);
    while (hadChanged) {
      hadChanged = false;
      for (size_t i = 0; i < net->nrOfNodes; ++i) {
        if (linkIsActive(i)) {
          continue;
        }
        for (const size_t springIdx : net->springIndicesOfLinks[i]) {
          if (springIsActive[springIdx]) {
            hadChanged = true;
            linkIsActive(i) = true;
            for (const size_t innerSpringIdx : net->springIndicesOfLinks[i]) {
              springIsActive[innerSpringIdx] = true;
            }
            break;
          }
        }
      }
    }

    return std::make_pair(springIsActive, linkIsActive);
  }

  /**
   * @brief Count the number of atoms that can be considered part of an
   * active cluster, i.e., are somehow connected to an active spring
   *
   * @param net
   * @param tolerance
   * @return double
   */
  double countActiveClusteredAtoms(ForceBalanceNetwork* net,
                                   const Eigen::VectorXd& u,
                                   const Eigen::VectorXd& springPartitions,
                                   const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(net->nrOfPartialSprings == springPartitions.size(),
                       "Spring partitions and network don't match");
    INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                       "Link displacements and network don't match");
    if (net->nrOfSprings < 1) {
      return 0.;
    }

    const std::vector<pylimer_tools::entities::Universe> clusters =
      this->universe.getClusters();
    std::vector<long int> atomIdxToClusterIdx(this->universe.getNrOfAtoms());
    for (size_t i = 0; i < clusters.size(); ++i) {
      for (const pylimer_tools::entities::Atom& atom : clusters[i].getAtoms()) {
        atomIdxToClusterIdx[this->universe.getIdxByAtomId(atom.getId())] = i;
      }
    }

    std::vector<bool> clusterIsActive(clusters.size(), false);

    // find active atoms
    const std::vector<long int> activeNodeIndices =
      this->getIndicesOfActiveNodes(net, u, springPartitions, tolerance);

    for (const long int& nodeIdx : activeNodeIndices) {
      const long int universeAtomIdx =
        this->universe.getIdxByAtomId(net->oldAtomIds[nodeIdx]);
      clusterIsActive[atomIdxToClusterIdx[universeAtomIdx]] = true;
    }

    double nClusteredAtoms = 0.;
    for (size_t i = 0; i < clusters.size(); ++i) {
      if (clusterIsActive[i]) {
        nClusteredAtoms += clusters[i].getNrOfAtoms();
      }
    }

    return nClusteredAtoms;
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
  double computeSolubleWeightFraction(ForceBalanceNetwork* net,
                                      const Eigen::VectorXd& u,
                                      const Eigen::VectorXd& springPartitions,
                                      const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(net->nrOfPartialSprings == springPartitions.size(),
                       "Spring partitions and network don't match");
    INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                       "Link displacements and network don't match");
    if (net->nrOfSprings < 1) {
      return 1.;
    }
    const double nActiveClusteredAtoms =
      this->countActiveClusteredAtoms(net, u, springPartitions, tolerance);
    // finally, normalize by the number of atoms.
    // NOTE: currently, the weight of the atoms is ignored
    return 1. - (nActiveClusteredAtoms /
                 (static_cast<double>(this->universe.getNrOfAtoms())));
  }

  /**
   * @brief Get the indices of active Nodes
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return std::vector<long int> the atom ids
   */
  std::vector<long int> getIndicesOfActiveNodes(
    const ForceBalanceNetwork* net,
    const Eigen::VectorXd& u,
    const Eigen::VectorXd& springPartitions,
    double tolerance = 1e-3) const;

  /**
   * @brief Get the Ids of active Nodes
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return std::vector<long int> the atom ids
   */
  std::vector<long int> getIdsOfActiveNodes(double tolerance = 1e-3) const;

  Eigen::VectorXd getCurrentSpringDistances() const
  {
    return this->currentSpringVectors;
  }

  std::vector<double> getCurrentSpringLengths() const
  {
    const Eigen::VectorXd vecs = this->getCurrentSpringDistances();

    return pylimer_tools::utils::segmentwise_norm(vecs, 3);
  }

  std::vector<double> getOverallSpringLengths() const
  {
    const std::vector<double> partialSpringDistances =
      this->getCurrentPartialSpringLengths();
    assert(partialSpringDistances.size() ==
           this->initialConfig.nrOfPartialSprings);
    std::vector<double> results =
      std::vector<double>(this->initialConfig.nrOfSprings, 0.);
    for (size_t i = 0; i < this->initialConfig.nrOfPartialSprings; ++i) {
      results[this->initialConfig.partialToFullSpringIndex[i]] +=
        partialSpringDistances[i];
    }

    return results;
  }

  Eigen::VectorXd getCurrentPartialSpringDistances() const
  {
    Eigen::VectorXd partialSpringVectors = this->evaluatePartialSpringVectors(
      this->initialConfig, this->currentDisplacements);

    return partialSpringVectors;
  }

  std::vector<double> getCurrentPartialSpringLengths() const
  {
    const Eigen::VectorXd vecs = this->evaluatePartialSpringVectors(
      this->initialConfig, this->currentDisplacements);

    return pylimer_tools::utils::segmentwise_norm(vecs, 3);
  }

  /**
   * @brief Get the Nr Of Active Springs connected to each node
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return Eigen::VectorXi
   */
  Eigen::VectorXi getNrOfActiveSpringsConnected(double tolerance = 1e-3) const;

  /**
   * @brief Get the Nr Of Active Springs connected to each node
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return Eigen::VectorXi
   */
  Eigen::VectorXi getNrOfActivePartialSpringsConnected(
    double tolerance = 1e-3) const;

  /**
   * @brief Get the Nr Of Active Springs object
   *
   * @param tol the tolerance: springs under a certain length are considered
   * inactive
   * @return int
   */
  int getNrOfActiveSprings(const double tolerance = 1e-3) const
  {
    return this->countNrOfActiveSprings(tolerance);
  }

  int getNrOfActiveSpringsInDir(const int dir,
                                const double tolerance = 1e-3) const
  {
    return this->countNrOfActiveSpringsInDir(dir, tolerance);
  }

  /**
   * @brief Get the Nr Of Active Springs object
   *
   * @param tol the tolerance: springs under a certain length are considered
   * inactive
   * @return int
   */
  int getNrOfActivePartialSprings(const double tolerance = 1e-3) const
  {
    return this->countNrOfActivePartialSprings(tolerance);
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
   * @brief Get the gamma factor at the current step
   *
   * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
   * computed as phantom = N<b^2>.
   * @param nrOfChains the nr of chains to average over (can be different
   * from the nr of springs thanks to omitted free chains or primary loops)
   * @return double
   */
  double getGammaFactor(double b02 = 0.96,
                        int nrOfChains = -1,
                        double oneOverSpringPartitionUpperLimit = 1.) const;
  double getGamma() override { return this->getGammaFactor(1., -1., 1.); }

  /**
   * @brief Get the per-(partial)-spring gamma factors
   *
   * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
   * computed as phantom = N<b^2>.
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getGammaFactors(
    double b02,
    double oneOverSpringPartitionUpperLimit = 1.) const;

  /**
   * @brief Get the per-(partial)-spring gamma factors
   *
   * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
   * computed as phantom = N<b^2>.
   * @param dir the direction (0=x, 1=y, 2=z)
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getGammaFactorsInDir(
    double b02,
    int dir,
    double oneOverSpringPartitionUpperLimit = 1.) const;

  /**
   * @brief Get the number of force balance iterations done so far
   *
   * @return int
   */
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
    return this->addSlipLinks(strandIdx1, strandIdx2, x, y, z, alphas, alphas);
  }

  void addSlipLinks(const std::vector<size_t>& strandIdx1,
                    const std::vector<size_t>& strandIdx2,
                    const std::vector<double>& x,
                    const std::vector<double>& y,
                    const std::vector<double>& z,
                    const std::vector<double>& alpha1,
                    const std::vector<double>& alpha2,
                    const bool clampAlpha = false)
  {
    const std::vector<std::vector<size_t>> loops;
    const std::vector<std::vector<size_t>> loopsOfSliplinks;
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
                                        const bool is2D,
                                        const bool assumeBoxLarge) const;

  Eigen::VectorXd evaluateSpringVectors(const ForceBalanceNetwork& net,
                                        const Eigen::VectorXd& u) const
  {
    return this->evaluateSpringVectors(
      net, u, this->is2D, this->assumeBoxLargeEnough);
  };

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
  Eigen::Vector3d evaluatePartialSpringDistance(const ForceBalanceNetwork& net,
                                                const Eigen::VectorXd& u,
                                                const size_t springIdx) const
  {
    return this->evaluatePartialSpringDistance(
      net, u, springIdx, this->is2D, this->assumeBoxLargeEnough);
  }

  /**
   * @brief Compute one spring length
   *
   * @param net
   * @param springIdx
   * @param is2D
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d evaluatePartialSpringDistance(const ForceBalanceNetwork& net,
                                                const Eigen::VectorXd& u,
                                                const size_t springIdx,
                                                const bool is2d,
                                                const bool boxLargeEnough) const
  {
    assert(net.isUpToDate);

    Eigen::Vector3d dist =
      ((net.coordinates.segment(3 * net.springPartIndexB(springIdx), 3) +
        u.segment(3 * net.springPartIndexB(springIdx), 3)) -
       (net.coordinates.segment(3 * net.springPartIndexA(springIdx), 3) +
        u.segment(3 * net.springPartIndexA(springIdx), 3))) +
      net.springPartBoxOffset.segment(3 * springIdx, 3);

    if (boxLargeEnough) {
      this->box.handlePBC<Eigen::Vector3d>(dist);
    }

    if (is2d) {
      dist[2] = 0.0;
    }

    return dist;
  }

  Eigen::VectorXd evaluatePartialSpringVectors(const ForceBalanceNetwork& net,
                                               const Eigen::VectorXd& u,
                                               const bool is2D,
                                               const bool assumeLarge) const;

  Eigen::VectorXd evaluatePartialSpringVectors(const ForceBalanceNetwork& net,
                                               const Eigen::VectorXd& u) const
  {
    return this->evaluatePartialSpringVectors(
      net, u, this->is2D, this->assumeBoxLargeEnough);
  };

  /**
   * @brief Sum the partitions up to a given link in a spring
   *
   * @param net
   * @param springPartition
   * @param springIdx
   * @param targetLink
   * @return double
   */
  static double sumToTotalFraction(const ForceBalanceNetwork& net,
                                   Eigen::VectorXd springPartition,
                                   const size_t springIdx,
                                   const size_t targetLink)
  {
    double alpha = 0.;
    for (size_t i = 0; i < net.localToGlobalSpringIndex[springIdx].size();
         ++i) {
      const size_t currentPartialSpringIdx =
        net.localToGlobalSpringIndex[springIdx][i];
      if (net.springPartIndexA[currentPartialSpringIdx] == targetLink) {
        return alpha;
      }
      alpha += springPartition[currentPartialSpringIdx];
      if (net.springPartIndexB[currentPartialSpringIdx] == targetLink) {
        return alpha;
      }
    }
    throw std::runtime_error("Did not find target link in spring.");
  }

  static size_t getOtherSpringIndex(const ForceBalanceNetwork& net,
                                    const size_t springIdx,
                                    const size_t linkIdx)
  {
    assert(net.springPartIndexA[springIdx] == linkIdx ||
           net.springPartIndexB[springIdx] == linkIdx);
    return net.springPartIndexA[springIdx] == linkIdx
             ? net.springPartIndexB[springIdx]
             : net.springPartIndexA[springIdx];
  }

  /**
   * @brief Query the box offset for a specific spring
   *
   * @param net
   * @param partialSpringIdx
   * @param linkIdx
   * @return Eigen::Vector3d
   */
  static Eigen::Vector3d getPartialSpringBoxOffset(
    const ForceBalanceNetwork& net,
    const size_t partialSpringIdx)
  {
    return net.springPartBoxOffset.segment(3 * partialSpringIdx, 3);
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
  Eigen::Vector3d getPartialSpringBoxOffsetTo(const ForceBalanceNetwork& net,
                                              const size_t partialSpringIdx,
                                              const size_t linkIdx) const
  {
    return (net.springPartIndexA(partialSpringIdx) == linkIdx)
             ? (-1. * this->getPartialSpringBoxOffset(net, partialSpringIdx))
             : (this->getPartialSpringBoxOffset(net, partialSpringIdx));
  }

  Eigen::Vector3d getPartialSpringBoxOffsetFrom(const ForceBalanceNetwork& net,
                                                const size_t partialSpringIdx,
                                                const size_t linkIdx) const
  {
    return -1. *
           this->getPartialSpringBoxOffsetTo(net, partialSpringIdx, linkIdx);
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
    const size_t linkIdx) const
  {
    return this->evaluatePartialSpringDistanceTo(
      net, u, springIdx, linkIdx, this->is2D, this->assumeBoxLargeEnough);
  }

  Eigen::Vector3d evaluatePartialSpringDistanceTo(
    const ForceBalanceNetwork& net,
    const Eigen::VectorXd& u,
    const size_t springIdx,
    const size_t linkIdx,
    const bool is2d,
    const bool boxLargeEnough) const
  {
    assert(this->isPartOfSpring(net, linkIdx, springIdx));

    const Eigen::Vector3d dist = this->evaluatePartialSpringDistance(
      net, u, springIdx, is2d, boxLargeEnough);

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
    const size_t linkIdx) const
  {
    return this->evaluatePartialSpringDistanceFrom(
      net, u, springIdx, linkIdx, this->is2D, this->assumeBoxLargeEnough);
  }

  Eigen::Vector3d evaluatePartialSpringDistanceFrom(
    const ForceBalanceNetwork& net,
    const Eigen::VectorXd& u,
    const size_t springIdx,
    const size_t linkIdx,
    const bool is2d,
    const bool boxLargeEnough) const
  {
    return -1. * this->evaluatePartialSpringDistanceTo(
                   net, u, springIdx, linkIdx, is2d, boxLargeEnough);
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
   * @brief Get the Weighted Partial Spring Length for one partial spring
   *
   * @return double
   */
  double getWeightedPartialSpringLength(
    const ForceBalanceNetwork& net,
    const Eigen::VectorXd& u,
    const Eigen::VectorXd& springPartitions,
    size_t partialSpringIdx,
    double oneOverSpringPartitionUpperLimit = 1.) const;

  /**
   * @brief Get the Weighted Partial Spring Lengths
   *
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getWeightedPartialSpringLengths(
    const double oneOverSpringPartitionUpperLimit = 1.) const
  {
    Eigen::VectorXd weightedLengths =
      Eigen::VectorXd(this->initialConfig.nrOfPartialSprings);
    for (size_t i = 0; i < this->initialConfig.nrOfPartialSprings; ++i) {
      weightedLengths(i) =
        this->getWeightedPartialSpringLength(this->initialConfig,
                                             this->currentDisplacements,
                                             this->currentSpringPartitionsVec,
                                             i,
                                             oneOverSpringPartitionUpperLimit);
    }

    return weightedLengths;
  }

  /**
   * @brief List all the partial spring indices that are connected to a
   * specified (slip/cross)link
   *
   * @param net
   * @param linkIdx
   * @return std::vector<size_t>
   */
  std::vector<size_t> getPartialSpringIndicesOfLink(
    const ForceBalanceNetwork& net,
    const size_t linkIdx) const
  {
    INVALIDARG_EXP_IFN(linkIdx < net.nrOfLinks,
                       "The requested link does not exist");
    std::vector<size_t> partialSpringIndices;

    const std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];

    for (const size_t springIdx : springIndices) {
      for (size_t partialSpringIdx : net.localToGlobalSpringIndex[springIdx]) {
        if (this->isPartOfSpring(net, linkIdx, partialSpringIdx)) {
          partialSpringIndices.push_back(partialSpringIdx);
        }
      }
    }
    return partialSpringIndices;
  }

  Eigen::VectorXd getForceMagnitudeVector(
    const double oneOverSpringPartitionUpperLimit = 1.0) const
  {
    Eigen::VectorXd forceMagnitude =
      Eigen::VectorXd::Zero(this->initialConfig.nrOfLinks);
    for (size_t i = 0; i < this->initialConfig.nrOfLinks; ++i) {
      forceMagnitude[i] =
        this->getForceOn(i, oneOverSpringPartitionUpperLimit).norm();
    }
    return forceMagnitude;
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
  static void setSpringpartitionIndicesOfSliplink(
    std::vector<size_t>& res_vec,
    const ForceBalanceNetwork& net,
    const size_t linkIdx);

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
    const double oneOverSpringPartitionUpperLimit = 1.0,
    const bool allowSlipLinksToPassEachOther = false) const
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
    const bool moveAttempt = false) const;

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
    const bool moveAttempt = false) const;

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
                               const bool respectLoops = true) const;

  /**
   * @brief Loop all springs, swap slip-links on them if they are close
   * enough
   *
   * @param net
   */
  void swapSlipLinks(ForceBalanceNetwork& net,
                     const Eigen::VectorXd& u,
                     Eigen::VectorXd& springPartitions,
                     double swappableCutoff) const;

  void swapSlipLinks(ForceBalanceNetwork& net,
                     const size_t partialSpringIdx) const;

  bool swapSlipLinkReversibly(
    ForceBalanceNetwork& net,
    Eigen::VectorXd& u,
    Eigen::VectorXd& springPartitions,
    const size_t partialSpringIdx,
    const double oneOverSpringPartitionUpperLimit = 1.0,
    const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
    const bool respectLoops = true,
    const bool moveAttempt = false) const;

  long int rotateSlipLinkAroundCrosslink(
    ForceBalanceNetwork& net,
    const Eigen::VectorXd& u,
    Eigen::VectorXd& springPartitions,
    const size_t partialSpringIdx,
    double oneOverSpringPartitionUpperLimit = 1.0,
    const bool respectLoops = true) const;

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
   * @brief Displace all links to the mean of all connected neighbours
   *
   * @param net the force balance network
   * @param u the current displacements, wherein the resulting coordinates
   * shall be stored
   * @return double, the distance (squared norm) displaced
   */
  double displaceToMeanPosition(
    const ForceBalanceNetwork& net,
    Eigen::VectorXd& u,
    const Eigen::ArrayXd& oneOverSpringPartitions) const;

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
  static Eigen::VectorXd assembleOneOverSpringPartition(
    const ForceBalanceNetwork& net,
    const Eigen::VectorXd& springPartitions0,
    const double oneOverSpringPartitionUpperLimit = 1.0);

  double getDisplacementResidualNorm(
    const double oneOverSpringPartitionUpperLimit = 1.0) const;

  double getResidual() override
  {
    // this is for the output
    return this->getDisplacementResidualNorm();
  }

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
  static std::vector<size_t> getNeighbourLinkIndices(
    const ForceBalanceNetwork& net,
    const size_t linkIdx)
  {
    std::vector<size_t> results;
    results.reserve(4);
    for (const size_t springIdx : net.springIndicesOfLinks[linkIdx]) {
      for (const size_t partialSpringIdx :
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

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& ar, std::uint32_t const version)
  {
    ar(cereal::virtual_base_class<OutputSupportingSimulation>(this));

    // properties
    ar(universe,
       box,
       exitReason,
       simulationHasRun,
       initialConfig,
       currentDisplacements,
       currentSpringVectors,
       currentPartialSpringVectors,
       currentSpringPartitionsVec);
    // configuration
    ar(is2D,
       assumeBoxLargeEnough,
       kappa,
       crossLinkerType,
       nrOfStepsDone,
       simplificationFrequency,
       entanglementType,
       defaultBondLength,
       springBreakingLength);
  }

  static MEHPForceBalance readRestartFile(std::string& file)
  {
    MEHPForceBalance res = MEHPForceBalance();
    pylimer_tools::utils::deserializeFromFile(res, file);
    return res;
  }

  void writeRestartFile(std::string& file) override
  {
    pylimer_tools::utils::serializeToFile<MEHPForceBalance>(*this, file);
  }

  static MEHPForceBalance constructFromString(std::string s)
  {
    MEHPForceBalance res = MEHPForceBalance();
    pylimer_tools::utils::deserializeFromString(res, s);
    return res;
  }
#endif

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
      Eigen::VectorXd::Zero(this->currentSpringVectors.size() / 3);

    for (size_t i = 0; i < this->currentSpringVectors.size() / 3; ++i) {
      const double b = lens.segment(3 * i, 3).norm();
      lens[i] = b;
    }

    return lens;
  }

  Eigen::VectorXd getCoordinates() override
  {
    return this->initialConfig.coordinates + this->currentDisplacements;
  }

  Eigen::VectorXd getInitialCoordinates()
  {
    return this->initialConfig.coordinates;
  }

  double getTemperature() override
  {
    std::cerr << "Warning: Temperature is not a reasonable metric for this "
                 "type of computation."
              << std::endl;
    return -1; // TODO: implement?
  }

  size_t getNumParticles() override { return this->initialConfig.nrOfNodes; }

  void debugAtomVicinity(const size_t atomId) const
  {
    long int atomIdx = -1;
    for (size_t i = 0; i < this->initialConfig.oldAtomIds.size(); ++i) {
      if (this->initialConfig.oldAtomIds[i] == atomId) {
        atomIdx = i;
        break;
      }
    }
    RUNTIME_EXP_IFN(atomIdx >= 0, "Atom not found.");
    std::cout << "Atom " << atomIdx << " (" << atomId << ")"
              << " connectivity:" << std::endl;
    for (const long int parentSpringIdx :
         this->initialConfig.springIndicesOfLinks[atomIdx]) {
      std::vector<size_t> allSpringIndices = this->getAllFullSpringIndicesAlong(
        this->initialConfig, parentSpringIdx);
      std::string prefix = "";
      for (const size_t springIdx : allSpringIndices) {
        prefix += "\t";
        std::cout << prefix << "Spring " << springIdx << " (";
        std::cout << this->initialConfig.springPartIndexA[springIdx] << " ⟷ "
                  << this->initialConfig.springPartIndexB[springIdx];
        std::cout << ") with distance "
                  << this->currentSpringVectors.segment(springIdx * 3, 3).norm()
                  << std::endl;
        std::cout << prefix << "\t";

        for (const long int linkIdx :
             this->initialConfig.linkIndicesOfSprings[springIdx]) {
          std::cout << linkIdx << " ";
          if (linkIdx < this->initialConfig.nrOfNodes) {
            std::cout << "(" << this->initialConfig.oldAtomIds[linkIdx] << ") ";
          }
        }
        std::cout << std::endl;
      }
    }
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
  bool ConvertNetwork(ForceBalanceNetwork& net,
                      const int crossLinkerType,
                      bool remove2functionalCrosslinkers = false,
                      bool removeDanglingChains = false);

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
    const auto stressTensor =
      this->evaluateStressTensor(net, u, springPartitions);
    return this->evaluatePressure(stressTensor);
  }

  /**
   * @brief Evaluate the pressure from the stress tensor
   *
   * @param stressTensor
   * @return double
   */
  static double evaluatePressure(
    const std::array<std::array<double, 3>, 3>& stressTensor)
  {
    return (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.0;
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
    const double oneOverSpringPartitionUpperLimit = 1.0) const;

  /**
   * @brief Compute the force acting on a slip- or crosslink
   *
   * TODO: use "global" partial distances
   *
   * @param linkIdx
   * @param net
   * @param u
   * @param springPartitions
   * @param minCutoff
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d evaluateForceOnLink(
    const size_t linkIdx,
    const ForceBalanceNetwork& net,
    const Eigen::VectorXd& u,
    const Eigen::VectorXd& springPartitions,
    Eigen::VectorXi& debugNrSpringsVisited,
    const double oneOverSpringPartitionUpperLimit = 1.0) const;

  /**
   * @brief Compute the stress acting on a slip- or crosslink
   *
   * TODO: use "global" partial distances
   *
   * @param linkIdx
   * @param net
   * @param u
   * @param springPartitions
   * @param minCutoff
   * @return Eigen::Vector3d
   */
  Eigen::Matrix3d evaluateStressOnLink(
    const size_t linkIdx,
    const ForceBalanceNetwork& net,
    const Eigen::VectorXd& u,
    const Eigen::VectorXd& springPartitions,
    Eigen::VectorXi& debugNrSpringsVisited,
    const double oneOverSpringPartitionUpperLimit = 1.0) const;

  /**
   * @brief Count how many of the springs are active (length > tolerance)
   *
   * @param tolerance
   * @return int
   */
  int countNrOfActiveSprings(const ForceBalanceNetwork* net,
                             const Eigen::VectorXd& u,
                             const Eigen::VectorXd& springPartitions,
                             const double tolerance = 1e-3) const
  {
    return (this->findActiveSprings(net, u, springPartitions, tolerance))
      .count();
  }

  int countNrOfActiveSprings(const double tolerance = 1e-3) const
  {
    return (this->findActiveSprings(tolerance) == true).count();
  }

  int countNrOfActiveSpringsInDir(const int dir,
                                  const double tolerance = 1e-3) const
  {
    return (this->findActiveSpringsInDir(dir, tolerance) == true).count();
  }

  int countNrOfActivePartialSprings(const double tolerance = 1e-3) const
  {
    return (this->findActivePartialSprings(tolerance) == true).count();
  }

  /**
   * @brief Determine for each spring whether the spring contains at least
   * one partial spring that is considered active (tolerance criterion)
   *
   * @param net
   * @param u
   * @param springPartitions
   * @param tolerance
   * @return Eigen::ArrayXb
   */
  Eigen::ArrayXb findActiveSprings(const ForceBalanceNetwork* net,
                                   const Eigen::VectorXd& u,
                                   const Eigen::VectorXd& springPartitions,
                                   const double tolerance = 1e-3) const
  {
    Eigen::ArrayXb activePartialSprings =
      this->findActivePartialSprings(net, u, springPartitions, tolerance);
    const double nActivePartialSprings = activePartialSprings.count();
    Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net->nrOfSprings, false);

    for (size_t i = 0; i < net->nrOfPartialSprings; ++i) {
      result[net->partialToFullSpringIndex[i]] =
        result[net->partialToFullSpringIndex[i]] || activePartialSprings[i];
    }

    const double nActiveSprings = result.count();
    return result;
  }

  Eigen::ArrayXb findActiveSprings(const double tolerance = 1e-3) const
  {
    return this->findActiveSprings(&this->initialConfig,
                                   this->currentDisplacements,
                                   this->currentSpringPartitionsVec,
                                   tolerance);
  }

  /**
   * @brief Determine for each spring whether the spring contains at least
   * one partial spring that is considered active (tolerance criterion)
   *
   * @param net
   * @param u
   * @param springPartitions
   * @param dir
   * @param tolerance
   * @return Eigen::ArrayXb
   */
  Eigen::ArrayXb findActiveSpringsInDir(const ForceBalanceNetwork* net,
                                        const Eigen::VectorXd& u,
                                        const Eigen::VectorXd& springPartitions,
                                        const int dir,
                                        const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(dir >= 0 && dir < 3, "Invalid direction");
    Eigen::VectorXd partialSpringVectors = this->evaluatePartialSpringVectors(
      *net, u, this->is2D, this->assumeBoxLargeEnough);
    Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net->nrOfSprings, false);

    for (size_t i = 0; i < net->nrOfPartialSprings; ++i) {
      result[net->partialToFullSpringIndex[i]] =
        result[net->partialToFullSpringIndex[i]] ||
        !this->distanceIsWithinTolerance(
          Eigen::Vector3d(partialSpringVectors[3 * i + dir], 0, 0),
          tolerance,
          net->springsContourLength[net->partialToFullSpringIndex[i]],
          springPartitions[i]);
    }

    return result;
  }

  Eigen::ArrayXb findActiveSpringsInDir(const int dir,
                                        const double tolerance = 1e-3) const
  {
    return this->findActiveSpringsInDir(&this->initialConfig,
                                        this->currentDisplacements,
                                        this->currentSpringPartitionsVec,
                                        dir,
                                        tolerance);
  }

  Eigen::ArrayXb findActiveNodes(const double tolerance = 1e-3) const
  {
    Eigen::ArrayXb activeSprings = this->findActiveSprings(tolerance);
    assert(activeSprings.size() == this->initialConfig.nrOfSprings);

    Eigen::ArrayXb activeNodes =
      Eigen::ArrayXb::Constant(this->initialConfig.nrOfNodes, false);
    for (size_t i = 0; i < activeSprings.size(); ++i) {
      if (activeSprings[i]) {
        activeNodes[this->initialConfig.springIndexA[i]] = true;
        activeNodes[this->initialConfig.springIndexB[i]] = true;
      }
    }
    return activeNodes;
  }

  Eigen::ArrayXb findActivePartialSprings(
    const ForceBalanceNetwork* net,
    const Eigen::VectorXd& u,
    const Eigen::VectorXd& springPartitions,
    const double tolerance = 1e-3) const
  {
    Eigen::VectorXd partialSpringVectors = this->evaluatePartialSpringVectors(
      *net, u, this->is2D, this->assumeBoxLargeEnough);
    Eigen::ArrayXb result =
      Eigen::ArrayXb::Constant(net->nrOfPartialSprings, false);

    for (size_t i = 0; i < net->nrOfPartialSprings; ++i) {
      result[i] = !this->distanceIsWithinTolerance(
        partialSpringVectors.segment(3 * i, 3),
        tolerance,
        net->springsContourLength[net->partialToFullSpringIndex[i]],
        springPartitions[i]);
    }

    return result;
  }

  Eigen::ArrayXb findActivePartialSprings(const double tolerance = 1e-3) const
  {
    return this->findActivePartialSprings(&this->initialConfig,
                                          this->currentDisplacements,
                                          this->currentSpringPartitionsVec,
                                          tolerance);
  }

  /**
   * @brief Sets the spring into the network
   *
   * @param net
   * @param springIdx
   * @param linkFrom
   * @param linkTo
   */
  static void registerPartialSpring(ForceBalanceNetwork& net,
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
  void setLinkPropertiesFromAtoms(ForceBalanceNetwork& net,
                                  const size_t linkIdx,
                                  const pylimer_tools::entities::Atom& atom1,
                                  const pylimer_tools::entities::Atom& atom2,
                                  const int atomType)
  {
    assert(linkIdx < net.nrOfLinks);
    assert(atom1.getType() != this->crossLinkerType &&
           atom2.getType() != this->crossLinkerType);
    if (atom1.getId() > atom2.getId()) {
      // make sure a second call to this function would result in same
      // result
      setLinkPropertiesFromAtoms(net, linkIdx, atom2, atom1, atomType);
      return;
    }

    Eigen::Vector3d coords = atom1.getCoordinates();
    Eigen::Vector3d dist = atom2.getCoordinates() - coords;
    this->box.handlePBC(dist);
    coords += 0.5 * dist;

    net.coordinates.segment(3 * linkIdx, 3) = coords;
    net.linkIsSliplink[linkIdx] = atomType != this->crossLinkerType;
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
                                 int atomType = -1) const
  {
    if (atomType == -1) {
      atomType = atom.getType();
    }

    Eigen::Vector3d coords = atom.getCoordinates();
    this->box.handlePBC(coords);
    net.coordinates.segment(3 * linkIdx, 3) = coords;
    net.linkIsSliplink[linkIdx] = atomType != this->crossLinkerType;
    if (!net.linkIsSliplink[linkIdx]) {
      net.oldAtomIds[linkIdx] = atom.getId();
      net.oldAtomTypes[linkIdx] = atom.getType();
    }
  }

  /**
   * @brief Decide whether a distance is within a given tolerance. This is
   * *the* distance criterion for determining whether a spring is active.
   *
   * @param dist
   * @param tolerance
   * @param contourLength
   * @param contourLengthFraction
   * @return true
   * @return false
   */
  static bool distanceIsWithinTolerance(const Eigen::Vector3d& dist,
                                        const double tolerance = 1e-3,
                                        const double contourLength = 1.,
                                        const double contourLengthFraction = 1.)
  {
    return dist.norm() <=
           (tolerance * std::max(contourLengthFraction * contourLength, 1.));
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
    const size_t partialSpringIdx) const
  {
    assert(atom1Idx < atom2Idx);
    if (chain.getType() ==
        pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
      assert(APPROX_WITHIN(atom1Idx, 0, chain.getNrOfAtoms() - 1, 1e-2));
      assert(APPROX_WITHIN(atom2Idx, 1, chain.getNrOfAtoms(), 1e-2));
    } else {
      assert(APPROX_WITHIN(atom1Idx, 0, chain.getNrOfAtoms() - 2, 1e-2));
      assert(APPROX_WITHIN(atom2Idx, 1, chain.getNrOfAtoms() - 1, 1e-2));
    }

    // set the partition, remember the partial spring mapping
    const size_t from = net.springPartIndexA[partialSpringIdx];
    const size_t to = net.springPartIndexB[partialSpringIdx];
    const std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
      chain.getAtomsLinedUp(crossLinkerType, false, true);
    partitions[partialSpringIdx] = static_cast<double>(atom2Idx - atom1Idx) /
                                   static_cast<double>(chain.getNrOfBonds());
    net.partialToFullSpringIndex[partialSpringIdx] = springIdx;
    net.partialSpringIsPartial[partialSpringIdx] =
      !(atom1Idx == 0 &&
        atom2Idx == chain.getNrOfAtoms() -
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
                                    crossLinkerType);
    Eigen::Vector3d additionalDistance1 =
      linedUpAtoms[atom1Idx].getCoordinates() -
      net.coordinates.segment(3 * from, 3);
    this->box.handlePBC(additionalDistance1);
    Eigen::Vector3d additionalDistance2 =
      net.coordinates.segment(3 * to, 3) -
      linedUpAtoms[atom2Idx].getCoordinates();
    this->box.handlePBC(additionalDistance2);
    expectedDistance += additionalDistance1 + additionalDistance2;
    if (this->is2D) {
      expectedDistance[2] = 0.0;
    }

    net.springPartBoxOffset.segment(3 * partialSpringIdx, 3) =
      Eigen::Vector3d::Zero();
    const Eigen::Vector3d actualDistance = this->evaluatePartialSpringDistance(
      net, Eigen::VectorXd::Zero(net.coordinates.size()), partialSpringIdx);
    net.springPartBoxOffset.segment(3 * partialSpringIdx, 3) =
      expectedDistance - actualDistance;
    assert(this->box.isValidOffset(expectedDistance - actualDistance));
#ifndef NDEBUG
    const Eigen::Vector3d newActualDistance =
      this->evaluatePartialSpringDistance(
        net, Eigen::VectorXd::Zero(net.coordinates.size()), partialSpringIdx);
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
    const size_t fullSpringIdx = net.partialToFullSpringIndex[partialSpringIdx];

    for (size_t i = 0; i < net.localToGlobalSpringIndex[fullSpringIdx].size();
         ++i) {
      if (net.localToGlobalSpringIndex[fullSpringIdx][i] == partialSpringIdx) {
        if (i == 0) {
          return net.localToGlobalSpringIndex[fullSpringIdx][i + 1];
        }
        if (i >= (net.localToGlobalSpringIndex[fullSpringIdx].size() - 1)) {
          return net.localToGlobalSpringIndex[fullSpringIdx][i - 1];
        }
        const size_t candidate1 =
          net.localToGlobalSpringIndex[fullSpringIdx][i - 1];
        const size_t candidate2 =
          net.localToGlobalSpringIndex[fullSpringIdx][i + 1];

        if (this->isPartOfSpring(net, aroundLinkIdx, candidate1) &&
            this->isPartOfSpring(net, aroundLinkIdx, candidate2) &&
            // we cannot handle the case where `partialSpringIdx` is a
            // primary loop without `notPartialSpringIdx`
            net.springPartIndexA[partialSpringIdx] !=
              net.springPartIndexB[partialSpringIdx]) {
          // this is a b-a-b-a situation
          // check the ordered linkIndicesOfSprings
          if (net.linkIndicesOfSprings[fullSpringIdx][i] == aroundLinkIdx) {
            return candidate1;
          } else {
            assert(net.linkIndicesOfSprings[fullSpringIdx][i + 1] ==
                   aroundLinkIdx);
            return candidate2;
          }
        }
        assert(this->isPartOfSpring(net, aroundLinkIdx, candidate1) ||
               this->isPartOfSpring(net, aroundLinkIdx, candidate2));
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
                     const size_t partialSpringIdx,
                     const size_t linkIdx) const
  {
    assert(this->isPartOfSpring(net, linkIdx, partialSpringIdx));
    return net.springPartIndexA[partialSpringIdx] == linkIdx
             ? net.springPartIndexB[partialSpringIdx]
             : net.springPartIndexA[partialSpringIdx];
  }

  static bool isPartOfSpring(const ForceBalanceNetwork& net,
                             const size_t linkIdx,
                             const size_t partialSpringIdx)
  {
    return (net.springPartIndexA[partialSpringIdx] == linkIdx) ||
           (net.springPartIndexB[partialSpringIdx] == linkIdx);
  }

  static bool isLoopingSpring(const ForceBalanceNetwork& net,
                              const size_t partialSpringIdx)
  {
    return (net.springPartIndexA[partialSpringIdx] ==
            net.springPartIndexB[partialSpringIdx]);
  }

  /**
   * @brief Check whether a spring is connected to an entanglement bead.
   * Possible to ignore a specific entanglement bead, e.g. to check
   * directionality.
   *
   * @param net
   * @param springIdx
   * @param ignoring
   * @return true
   * @return false
   */
  bool springInvolvesEntanglementBead(const ForceBalanceNetwork& net,
                                      const size_t springIdx,
                                      const long int ignoring = -1) const
  {
    INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                       "Spring index out of range.");
    if (net.springIndexA[springIdx] == net.springIndexB[springIdx]) {
      // "primary loop" e.g. if entanglement is with the same strand
      // cannot use ignore in this case
      return net.oldAtomTypes[net.springIndexA[springIdx]] ==
             this->entanglementType;
    }
    return (net.oldAtomTypes[net.springIndexA[springIdx]] ==
              this->entanglementType &&
            net.springIndexA[springIdx] != ignoring) ||
           (net.oldAtomTypes[net.springIndexB[springIdx]] ==
              this->entanglementType &&
            net.springIndexB[springIdx] != ignoring);
  }

  size_t getInvolvedEntanglementBeadIndex(const ForceBalanceNetwork& net,
                                          const size_t springIdx,
                                          const long int ignoring = -1) const
  {
    INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                       "Spring index out of range.");
    if (net.springIndexB[springIdx] == net.springIndexA[springIdx]) {
      // "primary loop" e.g. if entanglement is with the same strand
      // cannot use ignore in this case
      if (net.oldAtomTypes[net.springIndexA[springIdx]] ==
          this->entanglementType) {
        return net.springIndexA[springIdx];
      }
    }
    if (net.oldAtomTypes[net.springIndexA[springIdx]] ==
          this->entanglementType &&
        net.springIndexA[springIdx] != ignoring) {
      return net.springIndexA[springIdx];
    }
    if (net.oldAtomTypes[net.springIndexB[springIdx]] ==
          this->entanglementType &&
        net.springIndexB[springIdx] != ignoring) {
      return net.springIndexB[springIdx];
    }
    throw std::runtime_error("Did not find involved entanglement bead.");
  }

  /**
   * @brief Get the partial spring indices of all involved springs
   * from one crosslink to another crosslink, jumping all entanglement
   * beads.
   *
   * @return std::vector<size_t>
   */
  std::vector<size_t> getAllPartialSpringIndicesAlong(
    const ForceBalanceNetwork& net,
    const size_t springIdx) const
  {
    INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                       "Spring index out of range.");

    const std::vector<size_t> fullSpringIndices =
      this->getAllFullSpringIndicesAlong(net, springIdx);
    std::vector<size_t> result;
    for (const size_t springIdx : fullSpringIndices) {
      for (size_t partialSpringIdx : net.localToGlobalSpringIndex[springIdx]) {
        result.push_back(partialSpringIdx);
      }
    }

    return result;
  }

  /**
   * @brief Get all "full" spring indices (jumping entanglements) of a
   * particular crosslink.
   *
   * @param net
   * @param nodeIdx
   * @return std::vector<size_t>
   */
  std::vector<size_t> getInvolvedFullSpringIndices(
    const ForceBalanceNetwork& net,
    const size_t nodeIdx) const
  {
    INVALIDARG_EXP_IFN(!net.linkIsSliplink[nodeIdx], "Link must be crosslink.");
    INVALIDARG_EXP_IFN(net.oldAtomTypes[nodeIdx] != this->entanglementType,
                       "Link must be crosslink.");

    std::vector<size_t> result;
    for (const size_t springIdx : net.springIndicesOfLinks[nodeIdx]) {
      std::vector<size_t> subSprings =
        this->getAllFullSpringIndicesAlong(net, springIdx);
      for (size_t subSpringIdx : subSprings) {
        result.push_back(subSpringIdx);
      }
    }

    return result;
  }

  /**
   * @brief Get the full spring indices of all involved springs
   * from one crosslink to another crosslink, jumping all entanglement
   * beads.
   *
   * @return std::vector<size_t>
   */
  std::vector<size_t> getAllFullSpringIndicesAlong(
    const ForceBalanceNetwork& net,
    const size_t springIdx) const
  {
    INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                       "Spring index out of range.");
    INVALIDARG_EXP_IFN(net.springsType[springIdx] != this->entanglementType,
                       "Cannot follow entanglement springs");
    INVALIDARG_EXP_IFN(
      net.oldAtomTypes[net.springIndexA[springIdx]] != this->entanglementType ||
        net.oldAtomTypes[net.springIndexB[springIdx]] != this->entanglementType,
      "Cannot follow springs consisting of entanglement beads");

    std::vector<size_t> result = { springIdx };
    size_t currentSpringIdx = springIdx;
    long int previousEntanglementLinkIdx = -1;
    while (this->springInvolvesEntanglementBead(
      net, currentSpringIdx, previousEntanglementLinkIdx)) {
      const size_t entanglementLinkIdx = this->getInvolvedEntanglementBeadIndex(
        net, currentSpringIdx, previousEntanglementLinkIdx);
      RUNTIME_EXP_IFN(net.oldAtomTypes[entanglementLinkIdx] ==
                        this->entanglementType,
                      "Did not find involved entanglement bead.");
      // cannot assert, since this method might be called in a cleanup loop
      assert(net.springIndicesOfLinks[entanglementLinkIdx].size() == 3 ||
             net.springIndicesOfLinks[entanglementLinkIdx].size() == 2);
      long int nextSpringIdx = -1;
      for (const size_t involvedSpringIdx :
           net.springIndicesOfLinks[entanglementLinkIdx]) {
        if (net.springsType[involvedSpringIdx] == this->entanglementType) {
          continue;
        }
        if (involvedSpringIdx == currentSpringIdx) {
          continue;
        }
        nextSpringIdx = involvedSpringIdx;
      }
      RUNTIME_EXP_IFN(nextSpringIdx != currentSpringIdx,
                      "Did not find continuation spring.");
      RUNTIME_EXP_IFN(nextSpringIdx >= 0, "Spring index not found.");
      assert(net.springsType[nextSpringIdx] != this->entanglementType);
      currentSpringIdx = nextSpringIdx;
      result.push_back(currentSpringIdx);
      previousEntanglementLinkIdx = entanglementLinkIdx;
    }
    assert((net.oldAtomTypes[net.springIndexA[result[0]]] !=
              this->entanglementType ||
            net.oldAtomTypes[net.springIndexB[result[0]]] !=
              this->entanglementType) &&
           (net.oldAtomTypes[net.springIndexA[result.back()]] !=
              this->entanglementType ||
            net.oldAtomTypes[net.springIndexB[result.back()]] !=
              this->entanglementType));
    return result;
  }

  /**
   * @brief Get the indices of all involved entanglement beads
   * from one crosslink to another crosslink, jumping all entanglement
   * beads.
   *
   * @return std::vector<size_t>
   */
  std::vector<size_t> getEntanglementLinkIndicesAlong(
    const ForceBalanceNetwork& net,
    const size_t springIdx) const
  {
    INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                       "Spring index out of range.");
    INVALIDARG_EXP_IFN(net.springsType[springIdx] != this->entanglementType,
                       "Cannot follow entanglement springs");
    INVALIDARG_EXP_IFN(
      net.oldAtomTypes[net.springIndexA[springIdx]] != this->entanglementType ||
        net.oldAtomTypes[net.springIndexB[springIdx]] != this->entanglementType,
      "Cannot follow springs consisting of entanglement beads");

    std::vector<size_t> result = {};
    size_t currentSpringIdx = springIdx;
    long int previousEntanglementLinkIdx = -1;
    while (this->springInvolvesEntanglementBead(
      net, currentSpringIdx, previousEntanglementLinkIdx)) {
      size_t entanglementLinkIdx = this->getInvolvedEntanglementBeadIndex(
        net, currentSpringIdx, previousEntanglementLinkIdx);
      RUNTIME_EXP_IFN(net.oldAtomTypes[entanglementLinkIdx] ==
                        this->entanglementType,
                      "Did not find involved entanglement bead.");
      // cannot assert, since this method might be called in a cleanup loop
      assert(net.springIndicesOfLinks[entanglementLinkIdx].size() == 3 ||
             net.springIndicesOfLinks[entanglementLinkIdx].size() == 2);
      long int nextSpringIdx = -1;
      for (const size_t involvedSpringIdx :
           net.springIndicesOfLinks[entanglementLinkIdx]) {
        if (net.springsType[involvedSpringIdx] == this->entanglementType) {
          continue;
        }
        if (involvedSpringIdx == currentSpringIdx) {
          continue;
        }
        nextSpringIdx = involvedSpringIdx;
      }
      RUNTIME_EXP_IFN(nextSpringIdx != currentSpringIdx,
                      "Did not find continuation spring.");
      RUNTIME_EXP_IFN(nextSpringIdx >= 0, "Spring index not found.");
      currentSpringIdx = nextSpringIdx;
      result.push_back(entanglementLinkIdx);
      previousEntanglementLinkIdx = entanglementLinkIdx;
    }
    return result;
  }

  static double getDenominatorOfPartialSpring(
    const ForceBalanceNetwork& net,
    const Eigen::VectorXd& springPartitions,
    const size_t partialSpringIdx,
    const double oneOverSpringPartitionUpperLimit = 1.0);

  long int moveSlipLinkFromRailToRail(
    ForceBalanceNetwork& net,
    const Eigen::VectorXd& u,
    Eigen::VectorXd& springPartitions,
    const size_t sourcePartialSpringIdx,
    size_t targetPartialSpringIdx,
    const double oneOverSpringPartitionUpperLimit = 1.0) const
  {
#ifndef NDEBUG
    this->validateNetwork(net, u, springPartitions);
#endif
    INVALIDARG_EXP_IFN(sourcePartialSpringIdx != targetPartialSpringIdx,
                       "Source and target must be different when moving "
                       "from one crosslink branch to another.");

    const size_t involvedSlipLink =
      net.linkIsSliplink[net.springPartIndexA[sourcePartialSpringIdx]]
        ? net.springPartIndexA[sourcePartialSpringIdx]
        : net.springPartIndexB[sourcePartialSpringIdx];
    assert(net.linkIsSliplink[involvedSlipLink]);
    const size_t involvedCrossLink =
      net.linkIsSliplink[net.springPartIndexA[sourcePartialSpringIdx]]
        ? net.springPartIndexB[sourcePartialSpringIdx]
        : net.springPartIndexA[sourcePartialSpringIdx];
    assert(!net.linkIsSliplink[involvedCrossLink]);
    assert(
      (net.springPartIndexA[targetPartialSpringIdx] == involvedCrossLink) ||
      (net.springPartIndexB[targetPartialSpringIdx] == involvedCrossLink));
    const size_t targetFullSpringIdx =
      net.partialToFullSpringIndex[targetPartialSpringIdx];
    const size_t sourceFullSpringIdx =
      net.partialToFullSpringIndex[sourcePartialSpringIdx];

    assert(net.localToGlobalSpringIndex[sourceFullSpringIdx][0] ==
             sourcePartialSpringIdx ||
           pylimer_tools::utils::last(
             net.localToGlobalSpringIndex[sourceFullSpringIdx]) ==
             sourcePartialSpringIdx);

    size_t otherInvolvedPartialSpring =
      net.springPartIndexA[sourcePartialSpringIdx] == involvedCrossLink
        ? net.localToGlobalSpringIndex[sourceFullSpringIdx][1]
        : net.localToGlobalSpringIndex
            [sourceFullSpringIdx]
            [net.localToGlobalSpringIndex[sourceFullSpringIdx].size() - 2];

    // check whether there is space on the target spring, at all.
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
    const Eigen::Vector3d distanceBefore =
      this->evaluatePartialSpringDistance(
        net, u, otherInvolvedPartialSpring, this->is2D, false) +
      this->evaluatePartialSpringDistance(
        net, u, sourcePartialSpringIdx, this->is2D, false) +
      this->evaluatePartialSpringDistance(
        net, u, targetPartialSpringIdx, this->is2D, false);

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
    //           << " around crosslink " << involvedCrosslink
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

    const size_t newPartialSpringIdx =
      this->addSlipLinkToPartialSpring(net,
                                       u,
                                       springPartitions,
                                       targetPartialSpringIdx,
                                       involvedSlipLink,
                                       oneOverSpringPartitionUpperLimit);

    // finally, return the idx of the new partial spring
    long int resultingPartialSpringIdx = 0;
    long int remainingPartialSpringIdx = 0;
    if ((net.springPartIndexA[targetPartialSpringIdx] == involvedCrossLink &&
         net.springPartIndexB[targetPartialSpringIdx] == involvedSlipLink) ||
        (net.springPartIndexB[targetPartialSpringIdx] == involvedCrossLink &&
         net.springPartIndexA[targetPartialSpringIdx] == involvedSlipLink)) {
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
    const Eigen::Vector3d distanceAfter =
      this->evaluatePartialSpringDistance(
        net, u, otherInvolvedPartialSpring, this->is2D, false) +
      this->evaluatePartialSpringDistance(
        net, u, resultingPartialSpringIdx, this->is2D, false) +
      this->evaluatePartialSpringDistance(
        net, u, remainingPartialSpringIdx, this->is2D, false);
    assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
      distanceAfter, distanceBefore));

    return resultingPartialSpringIdx;
  };

  bool swapSlipLinksReversibly(
    ForceBalanceNetwork& net,
    Eigen::VectorXd& u,
    Eigen::VectorXd& springPartitions,
    const size_t partialSpringIdx,
    const double oneOverSpringPartitionUpperLimit = 1.0) const;

  bool swapSlipLinkWithXlinkReversibly(
    ForceBalanceNetwork& net,
    Eigen::VectorXd& u,
    Eigen::VectorXd& springPartitions,
    const size_t partialSpringIdx,
    const double oneOverSpringPartitionUpperLimit = 1.0,
    const bool respectLoops = true) const;
};
}

#ifdef CEREALIZABLE
CEREAL_REGISTER_TYPE(pylimer_tools::sim::mehp::MEHPForceBalance);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
  pylimer_tools::sim::OutputSupportingSimulation,
  pylimer_tools::sim::mehp::MEHPForceBalance);
CEREAL_CLASS_VERSION(pylimer_tools::sim::mehp::MEHPForceBalance, 1);
#endif
