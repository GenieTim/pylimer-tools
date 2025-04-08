#pragma once

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include "../sim/MEHPUtilityStructures.h"
#include "../sim/OutputSupportingSimulation.h"
#include "../topo/EntanglementDetector.h"
#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <iostream>
#include <nlopt.hpp>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>
#ifdef CEREALIZABLE
#include "../utils/CerealUtils.h"
#include <cereal/access.hpp>
#endif

namespace pylimer_tools::sim::mehp {
class MEHPForceBalance2 final
  : public pylimer_tools::sim::OutputSupportingSimulation
{
private:
#ifdef CEREALIZABLE
  MEHPForceBalance2()
    : initialConfig() {}; // not exposed to users, only used by Cereal

  friend class cereal::access;
#endif

  // member properties
  pylimer_tools::entities::Universe universe;
  pylimer_tools::entities::Box box;
  // state
  ExitReason exitReason = ExitReason::UNSET;
  bool simulationHasRun = false;
  ForceBalance2Network initialConfig;
  Eigen::VectorXd currentDisplacements;
  // configuration
  bool is2D = false;
  bool assumeBoxLargeEnough = true;
  double kappa = 1.0;
  int crossLinkerType = 2;
  int sliplinkType = 3;
  int nrOfStepsDone = 0;
  int simplificationFrequency = 10;
  double defaultBondLength = 0.0;
  double springBreakingLength = -1.;

public:
  explicit MEHPForceBalance2(const pylimer_tools::entities::Universe& u,
                             int crossLinkerType = 2,
                             bool is2D = false,
                             bool remove2functionalCrosslinkers = false,
                             bool removeDanglingChains = false)
    : universe(u)
  {
    this->crossLinkerType = crossLinkerType;
    this->box = u.getBox();
    // interpret network already to be able to give early results
    ForceBalance2Network net;
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

  MEHPForceBalance2(const ForceBalance2Network& net,

                    bool is2D = false)
  {
    this->is2D = is2D;
    this->initialConfig = net;
    this->currentDisplacements = Eigen::VectorXd::Zero(net.coordinates.size());
    this->box = pylimer_tools::entities::Box(net.L[0], net.L[1], net.L[2]);
    this->completeInitialization();
  }

#ifdef CEREALIZABLE
  static MEHPForceBalance2 constructFromString(std::string s)
  {
    MEHPForceBalance2 res = MEHPForceBalance2();
    pylimer_tools::utils::deserializeFromString(res, s);
    return res;
  }
#endif

  /**
   * @brief Instantiate this simulator with chosen entanglements.
   *
   * @param universe
   * @param entanglements
   * @param crossLinkerType
   * @param is2D
   * @return MEHPForceBalance2
   */
  static MEHPForceBalance2 constructWithEntanglements(
    const pylimer_tools::entities::Universe& universe,
    const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements&
      entanglements,
    int crossLinkerType = 2,
    bool is2D = false)
  {
    pylimer_tools::entities::Universe emptyUniverse =
      pylimer_tools::entities::Universe(universe.getBox());
    MEHPForceBalance2 fb =
      MEHPForceBalance2(emptyUniverse, crossLinkerType, is2D, false, false);
    fb.configAssumeBoxLargeEnough(false);
    fb.universe = universe;

    std::vector<std::pair<size_t, size_t>> pairsOfAtoms =
      entanglements.pairsOfAtoms;
    std::vector<long int> pairOfAtom = entanglements.pairOfAtom;

    // add ends of chains
    size_t currentVertexId = 0;
    size_t numUsableChains = 0;

    // TODO: implement

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
    fb.initialConfig.linkIsEntanglement.conservativeResize(
      fb.initialConfig.nrOfLinks);
    fb.initialConfig.strandIndicesOfLinks =
      pylimer_tools::utils::initializeWithValue(fb.initialConfig.nrOfLinks,
                                                std::vector<size_t>());

    // springs
    fb.initialConfig.nrOfStrands = numUsableChains;
    fb.initialConfig.springsContourLength.conservativeResize(numUsableChains);
    fb.initialConfig.springsType.conservativeResize(numUsableChains);
    fb.initialConfig.linkIndicesOfStrands =
      pylimer_tools::utils::initializeWithValue(numUsableChains,
                                                std::vector<size_t>());
    fb.initialConfig.springIndicesOfStrand =
      pylimer_tools::utils::initializeWithValue(numUsableChains,
                                                std::vector<size_t>());

    // partial springs
    // we don't know the actual number (yet), but we can over-estimate
    // pretty well, such that we only need to reduce afterwards
    size_t numPartialSpringsEstimate =
      numUsableChains + 2 * pairsOfAtoms.size();
    fb.initialConfig.nrOfSprings = numPartialSpringsEstimate;
    fb.initialConfig.springBoxOffset.conservativeResize(
      3 * numPartialSpringsEstimate);
    fb.initialConfig.springCoordinateIndexA.conservativeResize(
      3 * numPartialSpringsEstimate);
    fb.initialConfig.springCoordinateIndexB.conservativeResize(
      3 * numPartialSpringsEstimate);
    fb.initialConfig.springIndexA.conservativeResize(numPartialSpringsEstimate);
    fb.initialConfig.springIndexB.conservativeResize(numPartialSpringsEstimate);
    fb.initialConfig.strandIdxOfSpring.conservativeResize(
      numPartialSpringsEstimate);

    size_t springIdx = 0;
    size_t partialSpringIdx = 0;

    // TODO: implement

    fb.completeInitialization();

    return fb;
  };

  /**
   * @brief Instantiate this simulator with randomly chosen slip-links.
   *
   * @param universe the universe containing the basic atoms and connectivity
   * @param nrOfEntanglementsToSample the number of entanglements to sample
   * @param upperCutoff maximum distance from one sampled bead to its partner
   * @param lowerCutoff minimum distance from one sampled bead to its partner
   * @param minimumNrOfEntanglements the minimum number of entanglements that
   * should be sampled
   * @param sameStrandCutoff distance from one sampled bead to its pair within
   * the same strand
   * @param seed the seed for the random number generator
   * @param crossLinkerType
   * @param is2D
   * @return MEHPForceBalance2
   */
  static MEHPForceBalance2 constructWithRandomEntanglements(
    const pylimer_tools::entities::Universe& universe,
    const size_t nrOfEntanglementsToSample,
    const double upperCutoff,
    const double lowerCutoff = 0.,
    const size_t minimumNrOfEntanglements = 0,
    const double sameStrandCutoff = 3,
    const std::string seed = "",
    int crossLinkerType = 2,
    bool is2D = false,
    bool filterEntanglements = true)
  {
    // sample the "entanglements"
    pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
      entanglements =
        pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
          universe,
          nrOfEntanglementsToSample,
          upperCutoff,
          lowerCutoff,
          minimumNrOfEntanglements,
          sameStrandCutoff,
          seed,
          crossLinkerType,
          true,
          filterEntanglements);

    RUNTIME_EXP_IFN(
      entanglements.pairsOfAtoms.size() >= minimumNrOfEntanglements ||
        filterEntanglements,
      "Minimum number of slip-links could not be sampled: got " +
        std::to_string(entanglements.pairsOfAtoms.size()) + " instead of " +
        std::to_string(minimumNrOfEntanglements) + ".");

    return MEHPForceBalance2::constructWithEntanglements(
      universe, entanglements, crossLinkerType, is2D);
  }

  /**
   * @brief Finish initializing some member properties
   *
   */
  void completeInitialization()
  {
    this->defaultBondLength = universe.computeMeanBondLength();
    RUNTIME_EXP_IFN(this->validateNetwork(), "Invalid internal state");
  }

  /**
   * @brief Actually do run the simulation
   *
   * @param simplificationMode
   * @param inactiveRemovalCutoff
   * @param solver
   */
  void runForceRelaxation(const StructureSimplificationMode simplificationMode =
                            StructureSimplificationMode::NO_SIMPLIFICATION,
                          const double inactiveRemovalCutoff = 1e-3,
                          const SLESolver solver = SLESolver::DEFAULT)
  {
    this->runForceRelaxation(
      simplificationMode,
      inactiveRemovalCutoff,
      solver,
      []() { return false; },
      []() {});
  }

  /**
   * @brief Actually do run the simulation
   *
   * @param simplificationMode
   * @param inactiveRemovalCutoff
   * @param solver
   * @param shouldInterrupt
   * @param cleanupInterrupt
   */
  void runForceRelaxation(const StructureSimplificationMode simplificationMode,
                          const double inactiveRemovalCutoff,
                          const SLESolver solver,
                          const std::function<bool()>& shouldInterrupt,
                          const std::function<void()>& cleanupInterrupt);

  /**
   * @brief Remove cross-linkers, springs and associated slip-links with the
   * scheme suggested by Andrei
   *
   * @param net
   * @param displacements
   * @param tolerance
   * @return size_t
   */
  size_t doRemovalAndreisWay(ForceBalance2Network& net,
                             Eigen::VectorXd& displacements,
                             double tolerance) const;

  /**
   * @brief Remove crosslinkers which do not have any springs with a certain
   * minimum length
   *
   * @param net
   * @param displacements
   * @param tolerance
   */
  size_t removeInactiveCrosslinks(ForceBalance2Network& net,
                                  Eigen::VectorXd& displacements,

                                  double tolerance) const;

  /**
   * @brief Remove springs that exert a stress higher than
   * `this->springBreakingLength`
   *
   * @param net
   * @param displacements
   * @return size_t the number of springs broken
   */
  size_t breakTooLongSprings(ForceBalance2Network& net,
                             Eigen::VectorXd& displacements) const;

  /**
   * @brief Remove double listed springs from crosslinkers
   *
   * @param net
   */
  void removeDuplicateListedSpringsFromLinks(ForceBalance2Network& net) const;

  void removeDuplicateListedSpringsFromLink(
    ForceBalance2Network& net,
    size_t linkIdx,
    bool allowOnEntanglement = false) const;

  size_t removePrimaryLoops(ForceBalance2Network& net,
                            Eigen::VectorXd& displacements) const;

  /**
   * @brief Remove a spring (and all its parts, incl. slip-links) from the
   * structures
   *
   * @param net
   * @param displacements
   */
  void removeSpring(ForceBalance2Network& net,
                    Eigen::VectorXd& displacements,

                    const size_t springIdx) const;

  /**
   * @brief break a spring, given its partial spring index
   *
   * @param net
   * @param displacements

   * @param partialSpringIdx
   */
  void breakPartialSpring(ForceBalance2Network& net,
                          Eigen::VectorXd& displacements,
                          const size_t partialSpringIdx) const;

  /**
   * @brief Remove a spring, but also all springs that are connected to it
   * and are connected via entanglement links.
   *
   * @param net
   * @param displacements
   * @param springIdx
   */
  void removeSpringFollowingEntanglementLinks(ForceBalance2Network& net,
                                              Eigen::VectorXd& displacements,
                                              const size_t springIdx) const;

  /**
   * @brief Remove a certain link from the structures
   *
   * @param net
   * @param displacements
   * @param linkIdx
   */
  void removeLink(ForceBalance2Network& net,
                  Eigen::VectorXd& displacements,
                  const size_t linkIdx) const;

  /**
   * @brief Merge two springs around a given cross-link
   *
   * @param net
   * @param displacements
   * @param removedSpringIdx
   * @param keptSpringIdx
   * @param linkToReduce the index of the link to remove (combine the springs
   * around)
   */
  void mergeSprings(ForceBalance2Network& net,
                    const Eigen::VectorXd& displacements,
                    const size_t removedSpringIdx,
                    const size_t keptSpringIdx,
                    const size_t linkToReduce) const;

  /**
   * @brief Merge two springs around a given cross-link
   *
   * This does not require the resulting network to be valid.
   *
   * @param net
   * @param displacements
   * @param removedSpringIdx
   * @param keptSpringIdx
   * @param linkToReduce
   * @param skipEigenResize
   */
  void mergePartialSprings(ForceBalance2Network& net,
                           const Eigen::VectorXd& displacements,
                           const size_t removedSpringIdx,
                           const size_t keptSpringIdx,
                           const size_t linkToReduce,
                           bool skipEigenResize = false) const;

  /**
   * @brief Replace the two springs traversinga a two-functional
   * crosslinkers with a single spring
   *
   * @param net
   * @param displacements
   */
  size_t removeTwofunctionalCrosslinks(ForceBalance2Network& net,
                                       Eigen::VectorXd& displacements) const;

  /**
   * @brief Deform the system to match the specified box
   *
   * @param newBox the box to deform to
   */
  void deformTo(pylimer_tools::entities::Box& newBox)
  {
    this->box.adjustCoordinatesTo(this->initialConfig.coordinates, newBox);
    this->box.adjustCoordinatesTo(this->currentDisplacements, newBox);
    this->box.adjustCoordinatesTo(this->initialConfig.springBoxOffset, newBox);
    this->box = newBox;
    this->universe.setBox(newBox, true);
    for (size_t i = 0; i < 3; ++i) {
      this->initialConfig.L[i] = this->box.getL(i);
      this->initialConfig.boxHalfs[i] = 0.5 * this->initialConfig.L[i];
    }
  }

  /**
   * @brief Get the universe consisting of cross-linkers only
   *
   * @return pylimer_tools::entities::Universe
   */
  pylimer_tools::entities::Universe getCrosslinkerVerse() const;

  double getDefaultMeanBondLength() const { return this->defaultBondLength; }

  double getVolume() override
  {
    return this->initialConfig.L[0] * this->initialConfig.L[1] *
           this->initialConfig.L[2];
  }

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

  int getNrOfSprings() const { return this->initialConfig.nrOfStrands; }

  int getNrOfPartialSprings() const { return this->initialConfig.nrOfSprings; }

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

  void configAssumeBoxLargeEnough(bool assumption)
  {
    this->assumeBoxLargeEnough = assumption;
  }

  void configMeanBondLength(double meanBondLength)
  {
    this->defaultBondLength = meanBondLength;
  }

  void configSpringConstant(double kappa = 1.0) { this->kappa = kappa; }

  void configSimplificationFrequency(int newRemovalFrequency = 10)
  {
    this->simplificationFrequency = newRemovalFrequency;
  }

  void configSpringBreakingDistance(double newSpringBreakingForce = -1.)
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
  int getNrOfActiveNodes(double tolerance = 1e-3) const
  {
    return this->findActiveNodes(tolerance).count();
  }

  /**
   * @brief Get the Soluble Weight Fraction
   *
   * @param tolerance
   * @return double
   */
  double getSolubleWeightFraction(double tolerance = 1e-3)
  {
    return this->computeSolubleWeightFraction(
      &this->initialConfig, this->currentDisplacements, tolerance);
  }

  /**
   * @brief Get the Dangling Weight Fraction
   *
   * @param tolerance
   * @return double
   */
  double getDanglingWeightFraction(double tolerance = 1e-3)
  {
    return this->computeDanglingWeightFraction(
      &this->initialConfig, this->currentDisplacements, tolerance);
  }

  /**
   * @brief Get the Weight Fraction of Active Springs (atoms)
   *
   * @param tolerance
   * @return double
   */
  double getActiveWeightFraction(double tolerance = 1e-3)
  {
    return this->computeActiveWeightFraction(
      &this->initialConfig, this->currentDisplacements, tolerance);
  }

  /**
   * @brief Count the number of atoms that are in any way connected to an
   * active spring
   *
   * @param tolerance
   * @return double
   */
  double countActiveClusteredAtoms(double tolerance = 1e-3)
  {
    return this->countActiveClusteredAtoms(
      &this->initialConfig, this->currentDisplacements, tolerance);
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
  double computeDanglingWeightFraction(ForceBalance2Network* net,
                                       const Eigen::VectorXd& u,

                                       const double tolerance = 1e-3) const
  {
    double activeWeightFraction =
      this->computeActiveWeightFraction(net, u, tolerance);
    RUNTIME_EXP_IFN(
      APPROX_WITHIN(activeWeightFraction, 0., 1., 1e-6),
      "Expect active weight fraction to be between 0 and 1, got " +
        std::to_string(activeWeightFraction) + ".");
    double solubleWeightFraction =
      this->computeSolubleWeightFraction(net, u, tolerance);
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
  double computeActiveWeightFraction(ForceBalance2Network* net,
                                     const Eigen::VectorXd& u,

                                     const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                       "Link displacements and network don't match");
    if (net->nrOfStrands < 1) {
      return 0.;
    }
    // find all active springs
    Eigen::ArrayXb activeSprings = this->findActiveStrands(net, u, tolerance);
    if (activeSprings.count() == 0) {
      return 0.;
    }
    // as of now, the springsContourLength is equal to the number of bonds
    // from cross-link to cross-link. therefore, the number of atoms of each
    // of these springs is one less
    Eigen::ArrayXd allActiveAtomsPerChains =
      activeSprings.cast<double>() * (net->springsContourLength.array() -
                                      Eigen::ArrayXd::Ones(net->nrOfStrands));

    // TODO: currently, the weight of the atoms is ignored
    // normalize by the number of atoms
    return (allActiveAtomsPerChains.matrix().sum() +
            this->getNrOfActiveNodes(tolerance)) /
           (static_cast<double>(this->universe.getNrOfAtoms()));
  }

  /**
   * @brief Find whether springs and nodes are in any way connected to an
   * active spring
   *
   * @param net the network that includes the connectivity
   * @param u the current displacements of the links
   * @param tolerance the tolerance for considering springs as active
   * @return std::pair<Eigen::ArrayXb, Eigen::ArrayXb> indices of springs
   * (first) and links (second) connected in any way to active springs
   */
  std::pair<Eigen::ArrayXb, Eigen::ArrayXb> findClusteredToActive(
    const ForceBalance2Network* net,
    const Eigen::VectorXd& u,
    const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(u.size() == net->nrOfLinks * 3, "Invalid sizes.");

    // find all active springs
    Eigen::ArrayXb springIsActive = this->findActiveStrands(net, u, tolerance);

    // then, iteratively walk along the springs to mark those as "active"
    // that are connected to active springs
    bool hadChanged = true;
    Eigen::ArrayXb nodeIsActive =
      Eigen::ArrayXb::Constant(net->nrOfNodes, false);
    while (hadChanged) {
      hadChanged = false;
      for (size_t i = 0; i < net->nrOfNodes; ++i) {
        if (nodeIsActive(i)) {
          continue;
        }
        for (size_t springIdx : net->strandIndicesOfLinks[i]) {
          if (springIsActive[springIdx]) {
            hadChanged = true;
            nodeIsActive(i) = true;
            for (size_t innerSpringIdx : net->strandIndicesOfLinks[i]) {
              springIsActive[innerSpringIdx] = true;
            }
            break;
          }
        }
      }
    }

    return std::make_pair(springIsActive, nodeIsActive);
  }

  /**
   * @brief Count the number of atoms that can be considered part of an
   * active cluster, i.e., are somehow connected to an active spring
   *
   * @param net
   * @param tolerance
   * @return double
   */
  double countActiveClusteredAtoms(ForceBalance2Network* net,
                                   const Eigen::VectorXd& u,

                                   const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                       "Link displacements and network don't match");
    if (net->nrOfStrands < 1) {
      return 0.;
    }

    std::vector<pylimer_tools::entities::Universe> clusters =
      this->universe.getClusters();
    std::vector<long int> atomIdxToClusterIdx(this->universe.getNrOfAtoms());
    for (size_t i = 0; i < clusters.size(); ++i) {
      for (const pylimer_tools::entities::Atom& atom : clusters[i].getAtoms()) {
        atomIdxToClusterIdx[this->universe.getIdxByAtomId(atom.getId())] = i;
      }
    }

    std::vector<bool> clusterIsActive(clusters.size(), false);

    // find active atoms
    std::vector<long int> activeNodeIndices =
      this->getIndicesOfActiveNodes(net, u, tolerance);

    for (const long int& nodeIdx : activeNodeIndices) {
      long int universeAtomIdx =
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
   * @param u the current displacements of the links
   * @param tolerance
   * @return double
   */
  double computeSolubleWeightFraction(ForceBalance2Network* net,
                                      const Eigen::VectorXd& u,
                                      const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                       "Link displacements and network don't match");
    if (net->nrOfStrands < 1) {
      return 1.;
    }
    double nActiveClusteredAtoms =
      this->countActiveClusteredAtoms(net, u, tolerance);
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
  std::vector<long int> getIndicesOfActiveNodes(const ForceBalance2Network* net,
                                                const Eigen::VectorXd& u,

                                                double tolerance = 1e-3) const;

  /**
   * @brief Get the Ids of active Nodes
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return std::vector<long int> the atom ids
   */
  std::vector<long int> getIdsOfActiveNodes(double tolerance = 1e-3) const;

  /**
   *
   * @return a vector with the sum of the norms of all the springs per strand
   */
  std::vector<double> getOverallSpringLengths() const
  {
    std::vector<double> partialSpringDistances =
      this->getCurrentPartialSpringLengths();
    assert(partialSpringDistances.size() == this->initialConfig.nrOfSprings);
    std::vector<double> results =
      std::vector<double>(this->initialConfig.nrOfStrands, 0.);
    for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
      results[this->initialConfig.strandIdxOfSpring[i]] +=
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
    Eigen::VectorXd vecs = this->evaluatePartialSpringVectors(
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
   * @param tolerance springs under a certain length are considered inactive
   * @return Eigen::VectorXi
   */
  Eigen::VectorXi getNrOfActivePartialSpringsConnected(
    double tolerance = 1e-3) const;

  /**
   * @brief Get the Nr Of Active Springs object
   *
   * @param tolerance springs under a certain length are considered inactive
   * @return int
   */
  int getNrOfActiveStrands(double tolerance = 1e-3) const
  {
    return this->countNrOfActiveStrands(tolerance);
  }

  int getNrOfActiveStrandsInDir(int dir, double tolerance = 1e-3) const
  {
    return this->countNrOfActiveStrandsInDir(dir, tolerance);
  }

  /**
   * @brief Get the Nr Of Active Springs object
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return int
   */
  int getNrOfActivePartialSprings(double tolerance = 1e-3) const
  {
    return this->countNrOfActivePartialSprings(tolerance);
  }

  /**
   * @brief Get the Average Spring Length at the current step
   *
   * @return double
   */
  double getAverageSpringLength() const;

  Eigen::Matrix3d getStressTensor() override;

  Eigen::Matrix3d getStressTensorLinkBased(const bool xlinksOnly = false) const;

  /**
   * @brief Get the Pressure
   *
   * @return double
   */
  double getPressure() const
  {
    return this->evaluatePressure(this->initialConfig,
                                  this->currentDisplacements);
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
  double getGammaFactor(double b02 = 0.96, int nrOfChains = -1) const;

  double getGamma() override { return this->getGammaFactor(1., -1.); }

  /**
   * @brief Get the per-(partial)-spring gamma factors
   *
   * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
   * computed as phantom = N<b^2>.
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getGammaFactors(double b02) const;

  /**
   * @brief Get the per-(partial)-spring gamma factors
   *
   * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
   * computed as phantom = N<b^2>.
   * @param dir the direction (0=x, 1=y, 2=z)
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getGammaFactorsInDir(double b02, int dir) const;

  /**
   * @brief Get the number of force balance iterations done so far
   *
   * @return int
   */
  int getNrOfIterations() const { return this->nrOfStepsDone; }

  ExitReason getExitReason() const { return this->exitReason; }
  /**
   * @brief Compute one spring length
   *
   * @param net
   * @param u the current displacements of the links
   * @param springIdx
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d evaluatePartialSpringDistance(const ForceBalance2Network& net,
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
   * @param u the current displacements of the links
   * @param springIdx
   * @param is2d
   * @param boxLargeEnough
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d evaluatePartialSpringDistance(const ForceBalance2Network& net,
                                                const Eigen::VectorXd& u,
                                                const size_t springIdx,
                                                bool is2d,
                                                bool boxLargeEnough) const
  {
    Eigen::Vector3d dist =
      ((net.coordinates.segment(3 * net.springIndexB(springIdx), 3) +
        u.segment(3 * net.springIndexB(springIdx), 3)) -
       (net.coordinates.segment(3 * net.springIndexA(springIdx), 3) +
        u.segment(3 * net.springIndexA(springIdx), 3))) +
      net.springBoxOffset.segment(3 * springIdx, 3);

    if (boxLargeEnough) {
      this->box.handlePBC<Eigen::Vector3d>(dist);
    }

    if (is2d) {
      dist[2] = 0.0;
    }

    return dist;
  }

  /**
   *
   * @param net the network
   * @param u the current displacements
   * @param is2D
   * @param assumeLarge
   * @return the vectors of the springs
   */
  Eigen::VectorXd evaluatePartialSpringVectors(const ForceBalance2Network& net,
                                               const Eigen::VectorXd& u,
                                               const bool is2D,
                                               const bool assumeLarge) const;

  Eigen::VectorXd evaluatePartialSpringVectors(const ForceBalance2Network& net,
                                               const Eigen::VectorXd& u) const
  {
    return this->evaluatePartialSpringVectors(
      net, u, this->is2D, this->assumeBoxLargeEnough);
  };

  /**
   *
   * @param net the network
   * @param u the current displacements
   * @return the strand lengths (norm of the strand end-to-end vector)
   */
  Eigen::VectorXd evaluateSpringLengths(const ForceBalance2Network& net,
                                        Eigen::VectorXd u) const;

  /**
   *
   * @param net the network
   * @param u the current displacements
   * @return the strand end-to-end vectors
   */
  Eigen::VectorXd evaluateSpringVectors(const ForceBalance2Network& net,
                                        Eigen::VectorXd u) const;

  /**
   * @brief Sum the partitions up to a given link in a spring
   *
   * @param net
   * @param springPartition
   * @param springIdx
   * @param targetLink
   * @return double
   */
  double sumToTotalFraction(const ForceBalance2Network& net,
                            Eigen::VectorXd springPartition,
                            size_t springIdx,
                            size_t targetLink) const
  {
    double alpha = 0.;
    for (size_t i = 0; i < net.springIndicesOfStrand[springIdx].size(); ++i) {
      size_t currentPartialSpringIdx = net.springIndicesOfStrand[springIdx][i];
      if (net.springIndexA[currentPartialSpringIdx] == targetLink) {
        return alpha;
      }
      alpha += springPartition[currentPartialSpringIdx];
      if (net.springIndexB[currentPartialSpringIdx] == targetLink) {
        return alpha;
      }
    }
    throw std::runtime_error("Did not find target link in spring.");
  }

  size_t getOtherSpringIndex(const ForceBalance2Network& net,
                             const size_t springIdx,
                             const size_t linkIdx) const
  {
    assert(net.springIndexA[springIdx] == linkIdx ||
           net.springIndexB[springIdx] == linkIdx);
    return net.springIndexA[springIdx] == linkIdx ? net.springIndexB[springIdx]
                                                  : net.springIndexA[springIdx];
  }

  /**
   * @brief Query the box offset for a specific spring
   *
   * @param net
   * @param partialSpringIdx
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d getPartialSpringBoxOffset(const ForceBalance2Network& net,
                                            const size_t partialSpringIdx) const
  {
    return net.springBoxOffset.segment(3 * partialSpringIdx, 3);
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
  Eigen::Vector3d getPartialSpringBoxOffsetTo(const ForceBalance2Network& net,
                                              const size_t partialSpringIdx,
                                              const size_t linkIdx) const
  {
    return (net.springIndexA(partialSpringIdx) == linkIdx)
             ? (-1. * this->getPartialSpringBoxOffset(net, partialSpringIdx))
             : (this->getPartialSpringBoxOffset(net, partialSpringIdx));
  }

  Eigen::Vector3d getPartialSpringBoxOffsetFrom(const ForceBalance2Network& net,
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
   * @param u the current displacements
   * @param springIdx
   * @param linkIdx the vector "target"
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d evaluatePartialSpringDistanceTo(
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u,
    const size_t springIdx,
    const size_t linkIdx) const
  {
    return this->evaluatePartialSpringDistanceTo(
      net, u, springIdx, linkIdx, this->is2D, this->assumeBoxLargeEnough);
  }

  Eigen::Vector3d evaluatePartialSpringDistanceTo(
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u,
    const size_t springIdx,
    const size_t linkIdx,
    bool is2d,
    bool boxLargeEnough) const
  {
    assert(this->isPartOfSpring(net, linkIdx, springIdx));

    Eigen::Vector3d dist = this->evaluatePartialSpringDistance(
      net, u, springIdx, is2d, boxLargeEnough);

    return dist * (net.springIndexA(springIdx) == linkIdx ? -1. : 1.);
  }

  /**
   * @brief Compute one spring length, in a specific direction
   *
   * @param net
   * @param u the current displacements
   * @param springIdx
   * @param linkIdx the vector "source"
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d evaluatePartialSpringDistanceFrom(
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u,
    const size_t springIdx,
    const size_t linkIdx) const
  {
    return this->evaluatePartialSpringDistanceFrom(
      net, u, springIdx, linkIdx, this->is2D, this->assumeBoxLargeEnough);
  }

  Eigen::Vector3d evaluatePartialSpringDistanceFrom(
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u,
    const size_t springIdx,
    const size_t linkIdx,
    bool is2d,
    bool boxLargeEnough) const
  {
    return -1. * this->evaluatePartialSpringDistanceTo(
                   net, u, springIdx, linkIdx, is2d, boxLargeEnough);
  }

  ForceBalance2Network getNetwork() { return this->initialConfig; }

  /**
   * @brief Get the Weighted Partial Spring Length for one partial spring
   *
   * @return double
   */
  double getWeightedPartialSpringLength(const ForceBalance2Network& net,
                                        const Eigen::VectorXd& u,
                                        size_t partialSpringIdx) const;

  /**
   * @brief Get the Weighted Partial Spring Lengths
   *
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getWeightedPartialSpringLengths()
  {
    Eigen::VectorXd weightedLengths =
      Eigen::VectorXd(this->initialConfig.nrOfSprings);
    for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
      weightedLengths(i) = this->getWeightedPartialSpringLength(
        this->initialConfig, this->currentDisplacements, i);
    }

    return weightedLengths;
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
    const ForceBalance2Network& net,
    const size_t linkIdx) const
  {
    INVALIDARG_EXP_IFN(linkIdx < net.nrOfLinks,
                       "The requested link does not exist");
    std::unordered_set<size_t> partialSpringIndices;

    std::vector<size_t> springIndices = net.strandIndicesOfLinks[linkIdx];

    for (size_t spring_index = 0; spring_index < springIndices.size();
         ++spring_index) {
      std::vector<size_t> springsPartners =
        net.linkIndicesOfStrands[springIndices[spring_index]];
      for (size_t partner_idx = 0; partner_idx < springsPartners.size() - 1;
           ++partner_idx) {
        if (springsPartners[partner_idx] == linkIdx ||
            springsPartners[partner_idx + 1] == linkIdx) {
          size_t globalSpringIndex =
            net.springIndicesOfStrand[(springIndices[spring_index])]
                                     [partner_idx];
          partialSpringIndices.insert(globalSpringIndex);
        }
      }
    }
    return partialSpringIndices;
  }

  Eigen::VectorXd getForceMagnitudeVector() const
  {
    Eigen::VectorXd forceMagnitude =
      Eigen::VectorXd::Zero(this->initialConfig.nrOfLinks);
    for (size_t i = 0; i < this->initialConfig.nrOfLinks; ++i) {
      forceMagnitude[i] = this->getForceOn(i).norm();
    }
    return forceMagnitude;
  }

  /**
   * @brief Evaluate the force on one link
   *
   * @param index the link index
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d getForceOn(const size_t index) const
  {
    Eigen::VectorXi debugNrSpringsVisited =
      Eigen::VectorXi::Zero(this->initialConfig.nrOfSprings);
    return this->evaluateForceOnLink(index,
                                     this->initialConfig,
                                     this->currentDisplacements,
                                     debugNrSpringsVisited);
  }

  /**
   * @brief Evaluate the force on one link
   *
   * @param index the link index
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d getForceOn(const ForceBalance2Network& net,
                             const Eigen::VectorXd& u,
                             /* gives the parametrisation of N */
                             const size_t index) const
  {
    Eigen::VectorXi debugNrSpringsVisited = Eigen::VectorXi::Zero(0);
    return this->evaluateForceOnLink(index, net, u, debugNrSpringsVisited);
  }

  /**
   * @brief Evaluate the current stress on a particulas cross- or slip-link
   *
   * @param linkIdx
   * @return Eigen::Matrix3d
   */
  Eigen::Matrix3d getStressOn(const size_t linkIdx) const
  {
    Eigen::VectorXi debugNrSpringsVisited =
      Eigen::VectorXi::Zero(this->initialConfig.nrOfSprings);
    return this->evaluateStressOnLink(linkIdx,
                                      this->initialConfig,
                                      this->currentDisplacements,
                                      debugNrSpringsVisited);
  }

  /**
   * @brief Displace one link to the mean of all connected neighbours
   *
   * @param linkIdx the idx of the link to displace
   * @return the displacements afterwards
   */
  Eigen::VectorXd inspectDisplacementToMeanPositionUpdate(
    const size_t linkIdx) const
  {
    Eigen::VectorXd displacements = this->currentDisplacements;
    this->displaceToMeanPosition(this->initialConfig, displacements, linkIdx);
    return displacements;
  };

  /**
   * @brief Adjust the two spring's box offsets to work best with the
   * specified slip-link
   *
   * @param net the network to adjust
   * @param u the current displacements
   * @param slipLinkIdx the slip-link around which to adjust the two springs
   * @param spring1 one of the two partial spring idx
   * @param spring2 the partial spring idx of the other spring
   */
  void reAlignSlipLinkToImages(ForceBalance2Network& net,
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
   * @param oneOverSpringPartitions the one over contour length per direction
   * and spring vector
   * @return double, the distance (squared norm) displaced
   */
  double displaceToMeanPosition(
    const ForceBalance2Network& net,
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
  double displaceToMeanPosition(const ForceBalance2Network& net,
                                Eigen::VectorXd& u,
                                const size_t linkIdx) const;

  /**
   * @brief Translate the spring partition vector to its 3*size
   *
   * @param net
0
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd assembleOneOverSpringPartition(
    const ForceBalance2Network& net) const;

  double getDisplacementResidualNorm() const;

  double getResidual() override
  {
    // this is for the output
    return this->getDisplacementResidualNorm();
  }

  double getDisplacementResidualNormFor(const ForceBalance2Network& net,
                                        const Eigen::VectorXd& u) const;

  double getDisplacementResidualNormFor(
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u,
    const Eigen::VectorXd& oneOverSpringPartitions) const;

  /**
   * @brief Get the Link Indices of all neighbours of a specified link
   *
   * @param net
   * @param linkIdx
   * @return std::vector<size_t>
   */
  std::vector<size_t> getNeighbourLinkIndices(const ForceBalance2Network& net,
                                              const size_t linkIdx) const
  {
    std::vector<size_t> results;
    results.reserve(4);
    for (size_t springIdx : net.strandIndicesOfLinks[linkIdx]) {
      for (size_t partialSpringIdx : net.springIndicesOfStrand[springIdx]) {
        if (net.springIndexA[partialSpringIdx] == linkIdx) {
          results.push_back(net.springIndexB[partialSpringIdx]);
        } //
        else if (net.springIndexB[partialSpringIdx] == linkIdx) {
          results.push_back(net.springIndexA[partialSpringIdx]);
        }
      }
    }
    return results;
  }

#ifdef CEREALIZABLE
  void writeRestartFile(std::string& file) override
  {
    throw std::runtime_error("Restart not supported yet");
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

  int getNumShifts() override { return 0; }
  int getNumRelocations() override { return 0; }

  Eigen::VectorXd getBondLengths() override
  {
    return this->evaluatePartialSpringVectors(this->initialConfig,
                                              this->currentDisplacements);
  }

  Eigen::VectorXd getCoordinates() override
  {
    return this->initialConfig.coordinates + this->currentDisplacements;
  }

  double getTemperature() override
  {
    std::cerr << "Warning: Temperature is not a reasonable metric for this "
                 "type of computation."
              << std::endl;
    return 0;
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
    for (long int parentSpringIdx :
         this->initialConfig.strandIndicesOfLinks[atomIdx]) {
      std::vector<size_t> allSpringIndices =
        this->initialConfig.springIndicesOfStrand[parentSpringIdx];
      std::string prefix = "";
      for (size_t springIdx : allSpringIndices) {
        prefix += "\t";
        std::cout << prefix << "Spring " << springIdx << " (";
        std::cout << this->initialConfig.springIndexA[springIdx] << " ⟷ "
                  << this->initialConfig.springIndexB[springIdx];
        std::cout << prefix << "\t";

        for (long int linkIdx :
             this->initialConfig.linkIndicesOfStrands[springIdx]) {
          std::cout << linkIdx << " ";
          if (linkIdx < this->initialConfig.nrOfNodes) {
            std::cout << "(" << this->initialConfig.oldAtomIds[linkIdx] << ") ";
          }
        }
        std::cout << std::endl;
      }
    }
  }

  bool validateNetwork() const
  {
    return this->validateNetwork(this->initialConfig,
                                 this->currentDisplacements);
  }

  bool validateNetwork(const ForceBalance2Network& net) const
  {
    return this->validateNetwork(net, this->currentDisplacements);
  }

  bool validateNetwork(const ForceBalance2Network& net,
                       const Eigen::VectorXd& u) const;

protected:
  /**
   * @brief Convert the universe to a network
   *
   * @param net the target network
   * @param crossLinkerType the atom type of the crossLinker
   * @return true
   * @return false
   */
  bool ConvertNetwork(ForceBalance2Network& net,
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
  double evaluatePressure(const ForceBalance2Network& net,
                          const Eigen::VectorXd& u) const
  {
    auto stressTensor = this->evaluateStressTensor(net, u);
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
    return (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.0;
  }

  /**
   * @brief Compute the stress tensor
   *
   * @param net
   * @param u
   * @param crosslinksOnly whether to only consider cross-links
   * @return std::array<std::array<double, 3>, 3>
   */
  std::array<std::array<double, 3>, 3> evaluateStressTensorLinkBased(
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u,
    const bool crosslinksOnly = false) const;

  /**
   * @brief Compute the stress tensor
   *
   * @param linkIndices the indices of the links to respect
   * @param net
   * @param u
   * @return std::array<std::array<double, 3>, 3>
   */
  Eigen::Matrix3d evaluateStressTensorForLinks(
    const std::vector<size_t> linkIndices,
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u) const;

  /**
   * @brief Compute the stress tensor
   *
   * @param net
   * @param u
   * @return std::array<std::array<double, 3>, 3>
   */
  std::array<std::array<double, 3>, 3> evaluateStressTensor(
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u) const;

  /**
   * @brief Compute the force acting on a slip- or cross-link
   *
   * @param linkIdx
   * @param net
   * @param u
   * @param debugNrSpringsVisited a vector to keep track of visited springs
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d evaluateForceOnLink(
    const size_t linkIdx,
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u,
    Eigen::VectorXi& debugNrSpringsVisited) const;

  /**
   * @brief Compute the stress acting on a slip- or cross-link
   *
   * @param linkIdx
   * @param net
   * @param u the current link displacements
   * @param debugNrSpringsVisited a vector to keep track of visited springs
   * @return Eigen::Vector3d
   */
  Eigen::Matrix3d evaluateStressOnLink(
    const size_t linkIdx,
    const ForceBalance2Network& net,
    const Eigen::VectorXd& u,
    Eigen::VectorXi& debugNrSpringsVisited) const;

  /**
   * @brief Count how many of the springs are active (length > tolerance)
   *
   * @param net
   * @param u the displacements
   * @param tolerance
   * @return int
   */
  int countNrOfActiveStrands(const ForceBalance2Network* net,
                             const Eigen::VectorXd& u,
                             const double tolerance = 1e-3) const
  {
    return (this->findActiveStrands(net, u, tolerance)).count();
  }

  int countNrOfActiveStrands(const double tolerance = 1e-3) const
  {
    return (this->findActiveStrands(tolerance) == true).count();
  }

  int countNrOfActiveStrandsInDir(const int dir,
                                  const double tolerance = 1e-3) const
  {
    return (this->findActiveStrandsInDir(dir, tolerance) == true).count();
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

   * @param tolerance
   * @return Eigen::ArrayXb
   */
  Eigen::ArrayXb findActiveStrands(const ForceBalance2Network* net,
                                   const Eigen::VectorXd& u,
                                   const double tolerance = 1e-3) const
  {
    Eigen::VectorXd partialSpringVectors = this->evaluatePartialSpringVectors(
      *net, u, this->is2D, this->assumeBoxLargeEnough);
    Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net->nrOfStrands, false);

    for (size_t i = 0; i < net->nrOfSprings; ++i) {
      result[net->strandIdxOfSpring[i]] =
        result[net->strandIdxOfSpring[i]] ||
        !this->distanceIsWithinTolerance(
          partialSpringVectors.segment(3 * i, 3),
          tolerance,
          net->springsContourLength[net->strandIdxOfSpring[i]]);
    }

    return result;
  }

  Eigen::ArrayXb findActiveStrands(const double tolerance = 1e-3) const
  {
    return this->findActiveStrands(&this->initialConfig,
                                   this->currentDisplacements,

                                   tolerance);
  }

  /**
   * @brief Determine for each spring whether the spring contains at least
   * one partial spring that is considered active (tolerance criterion)
   *
   * @param net
   * @param u

   * @param dir
   * @param tolerance
   * @return Eigen::ArrayXb
   */
  Eigen::ArrayXb findActiveStrandsInDir(const ForceBalance2Network* net,
                                        const Eigen::VectorXd& u,

                                        const int dir,
                                        const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(dir >= 0 && dir < 3, "Invalid direction");
    Eigen::VectorXd partialSpringVectors = this->evaluatePartialSpringVectors(
      *net, u, this->is2D, this->assumeBoxLargeEnough);
    Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net->nrOfStrands, false);

    for (size_t i = 0; i < net->nrOfSprings; ++i) {
      result[net->strandIdxOfSpring[i]] =
        result[net->strandIdxOfSpring[i]] ||
        !this->distanceIsWithinTolerance(
          Eigen::Vector3d(partialSpringVectors[3 * i + dir], 0, 0),
          tolerance,
          net->springsContourLength[net->strandIdxOfSpring[i]]);
    }

    return result;
  }

  Eigen::ArrayXb findActiveStrandsInDir(const int dir,
                                        const double tolerance = 1e-3) const
  {
    return this->findActiveStrandsInDir(&this->initialConfig,
                                        this->currentDisplacements,

                                        dir,
                                        tolerance);
  }

  Eigen::ArrayXb findActiveNodes(const double tolerance = 1e-3) const
  {
    Eigen::ArrayXb activeSprings = this->findActiveStrands(tolerance);
    assert(activeSprings.size() == this->initialConfig.nrOfStrands);

    Eigen::ArrayXb activeNodes =
      Eigen::ArrayXb::Constant(this->initialConfig.nrOfNodes, false);
    for (size_t i = 0; i < activeSprings.size(); ++i) {
      if (activeSprings[i]) {
        activeNodes[this->initialConfig.linkIndicesOfStrands[i][0]] = true;
        activeNodes[pylimer_tools::utils::last(
          this->initialConfig.linkIndicesOfStrands[i])] = true;
      }
    }
    return activeNodes;
  }

  Eigen::ArrayXb findActivePartialSprings(const ForceBalance2Network* net,
                                          const Eigen::VectorXd& u,

                                          const double tolerance = 1e-3) const
  {
    Eigen::VectorXd partialSpringVectors = this->evaluatePartialSpringVectors(
      *net, u, this->is2D, this->assumeBoxLargeEnough);
    Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net->nrOfSprings, false);

    for (size_t i = 0; i < net->nrOfSprings; ++i) {
      result[i] = !this->distanceIsWithinTolerance(
        partialSpringVectors.segment(3 * i, 3),
        tolerance,
        net->springsContourLength[net->strandIdxOfSpring[i]]);
    }

    return result;
  }

  Eigen::ArrayXb findActivePartialSprings(const double tolerance = 1e-3) const
  {
    return this->findActivePartialSprings(&this->initialConfig,
                                          this->currentDisplacements,

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
  void registerPartialSpring(ForceBalance2Network& net,
                             const size_t springIdx,
                             const size_t linkFrom,
                             const size_t linkTo)
  {
    net.springIndexA[springIdx] = linkFrom;
    net.springIndexB[springIdx] = linkTo;
    for (size_t i = 0; i < 3; ++i) {
      net.springCoordinateIndexA[3 * springIdx + i] = linkFrom * 3 + i;
      net.springCoordinateIndexB[3 * springIdx + i] = linkTo * 3 + i;
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
  void setLinkPropertiesFromAtoms(ForceBalance2Network& net,
                                  const size_t linkIdx,
                                  const pylimer_tools::entities::Atom& atom1,
                                  const pylimer_tools::entities::Atom& atom2,
                                  int atomType)
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
    net.linkIsEntanglement[linkIdx] = atomType != this->crossLinkerType;
  }

  /**
   * @brief Set the Link Properties From an Atom object
   *
   * @param net
   * @param linkIdx
   * @param atom
   * @param atomType
   */
  void setLinkPropertiesFromAtom(ForceBalance2Network& net,
                                 const size_t linkIdx,
                                 const pylimer_tools::entities::Atom& atom,
                                 int atomType = -1)
  {
    if (atomType == -1) {
      atomType = atom.getType();
    }

    Eigen::Vector3d coords = atom.getCoordinates();
    this->box.handlePBC(coords);
    net.coordinates.segment(3 * linkIdx, 3) = coords;
    net.linkIsEntanglement[linkIdx] = atomType != this->crossLinkerType;
    if (!net.linkIsEntanglement[linkIdx]) {
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
  bool distanceIsWithinTolerance(const Eigen::Vector3d& dist,
                                 double tolerance = 1e-3,
                                 double contourLength = 1.,
                                 double contourLengthFraction = 1.) const
  {
    return dist.norm() <=
           (tolerance * std::max(contourLengthFraction * contourLength, 1.));
  }

  bool isPartOfSpring(const ForceBalance2Network& net,
                      size_t linkIdx,
                      size_t springIdx) const
  {
    return (net.springIndexA[springIdx] == linkIdx) ||
           (net.springIndexB[springIdx] == linkIdx);
  }

  bool isLoopingSpring(const ForceBalance2Network& net, size_t springIdx) const
  {
    return (net.springIndexA[springIdx] == net.springIndexB[springIdx]);
  }

  double getDenominatorOfPartialSpring(const ForceBalance2Network& net,
                                       const size_t partialSpringIdx) const;
};
} // namespace