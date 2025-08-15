#pragma once

#include "../entities/Atom.h"
#include "../entities/Box.h"
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
  double kappa = 1.0;
  int crossLinkerType = 2;
  int entanglementType = 3;
  int nrOfStepsDone = 0;
  double defaultBondLength = 0.0;
  double springBreakingLength = -1.;
  bool assumeNetworkIsComplete = false;

public:
  /**
   * @brief Instantiate this simulator with chosen entanglements.
   *
   * @param universe the universe containing atoms and connectivity
   * @param entanglements the detected entanglements to include in the network
   * @param crossLinkerType the atom type of crosslinkers
   * @param is2D whether to perform 2D simulation (z-coordinates ignored)
   * @param entanglementsAsSprings whether to model the entanglements as merged
   * beads or beads with 1 spring in between
   * @return MEHPForceBalance2
   */
  MEHPForceBalance2(
    const pylimer_tools::entities::Universe& universe,
    const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements&
      entanglements,
    int crossLinkerType = 2,
    bool is2D = false,
    bool entanglementsAsSprings = false);

  explicit MEHPForceBalance2(const pylimer_tools::entities::Universe& u,
                             const int crossLinkerType = 2,
                             const bool is2D = false)
    : MEHPForceBalance2(
        u,
        pylimer_tools::topo::entanglement_detection::AtomPairEntanglements(),
        crossLinkerType,
        is2D,
        false) {};

  /**
   * @brief Constructor from a network and 2D flag
   *
   * @param net the force balance network to use
   * @param is2D whether to perform 2D simulation (z-coordinates ignored)
   */
  MEHPForceBalance2(const ForceBalance2Network& net, const bool is2D = false)
  {
    this->is2D = is2D;
    this->initialConfig = net;
    this->box = pylimer_tools::entities::Box(net.L[0], net.L[1], net.L[2]);
    this->completeInitialization();
  }

  /**
   * @brief Constructor from universe and network
   *
   * @param universe the universe containing atoms and connectivity
   * @param net the force balance network to use
   * @param is2D whether to perform 2D simulation (z-coordinates ignored)
   */
  MEHPForceBalance2(const pylimer_tools::entities::Universe& universe,
                    const ForceBalance2Network& net,
                    const bool is2D = false)
  {
    this->universe = universe;
    this->is2D = is2D;
    this->initialConfig = net;
    this->box = pylimer_tools::entities::Box(net.L[0], net.L[1], net.L[2]);
    this->completeInitialization();
  }

  /**
   * @brief Constructor from universe, old network and spring partitions
   *
   * @param u the universe containing atoms and connectivity
   * @param net1 the force balance network to convert from
   * @param springPartitions the spring partition parameters
   * @param is2D whether to perform 2D simulation (z-coordinates ignored)
   */
  MEHPForceBalance2(const pylimer_tools::entities::Universe& u,
                    const ForceBalanceNetwork& net1,
                    Eigen::VectorXd springPartitions,
                    const bool is2D = false)
  {
    this->universe = u;
    this->is2D = is2D;

    pylimer_tools::sim::mehp::ForceBalance2Network net2 = {};
    for (size_t dir = 0; dir < 3; ++dir) {
      net2.L[dir] = net1.L[dir];
      net2.boxHalfs[dir] = net1.boxHalfs[dir];
    }
    net2.nrOfLinks = net1.nrOfLinks;
    net2.nrOfNodes = net1.nrOfNodes;
    net2.nrOfSprings = net1.nrOfPartialSprings;
    net2.nrOfStrands = net1.nrOfSprings;
    net2.springIndexA = net1.springPartIndexA;
    net2.springIndexB = net1.springPartIndexB;
    net2.springCoordinateIndexA = net1.springPartCoordinateIndexA;
    net2.springCoordinateIndexB = net1.springPartCoordinateIndexB;
    net2.coordinates = net1.coordinates;
    net2.linkIndicesOfStrand = net1.linkIndicesOfSprings;
    net2.springBoxOffset = net1.springPartBoxOffset;
    net2.springIndicesOfStrand = net1.localToGlobalSpringIndex;
    net2.linkIsEntanglement = net1.linkIsSliplink;
    net2.strandIndexOfSpring = net1.partialToFullSpringIndex;
    net2.strandIndicesOfLink = net1.springIndicesOfLinks;

    // what needs a bit more translation
    net2.springIsEntanglement =
      Eigen::ArrayXb::Zero(static_cast<Eigen::Index>(net2.nrOfSprings));
    net2.springContourLength =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(net2.nrOfSprings));
    for (size_t i = 0; i < net2.nrOfSprings; ++i) {
      net2.springContourLength[static_cast<Eigen::Index>(i)] =
        net1.springsContourLength
          [net1.partialToFullSpringIndex[static_cast<Eigen::Index>(i)]] *
        springPartitions[static_cast<Eigen::Index>(i)];
    }
    net2.oldAtomTypes =
      Eigen::VectorXi::Zero(static_cast<Eigen::Index>(net2.nrOfLinks));
    net2.oldAtomIds =
      Eigen::VectorXi::Constant(static_cast<Eigen::Index>(net2.nrOfLinks), -1);
    for (size_t i = 0; i < net2.nrOfNodes; ++i) {
      net2.oldAtomIds[static_cast<Eigen::Index>(i)] =
        net1.oldAtomIds[static_cast<Eigen::Index>(i)];
      net2.oldAtomTypes[static_cast<Eigen::Index>(i)] =
        net1.oldAtomTypes[static_cast<Eigen::Index>(i)];
    }

    this->initialConfig = net2;
    this->box = pylimer_tools::entities::Box(net1.L[0], net1.L[1], net1.L[2]);
    this->completeInitialization();
  }

  /**
   * @brief Constructor from old network and spring partitions
   *
   * @param net1 the force balance network to convert from
   * @param springPartitions the spring partition parameters
   * @param is2D whether to perform 2D simulation (z-coordinates ignored)
   */
  MEHPForceBalance2(const ForceBalanceNetwork& net1,
                    const Eigen::VectorXd& springPartitions,
                    const bool is2D = false)
    : MEHPForceBalance2(
        pylimer_tools::entities::Universe(net1.L[0], net1.L[1], net1.L[2]),
        net1,
        springPartitions,
        is2D) {};

  /**
   * @brief Construct from old force relaxation network
   *
   * @param u
   * @param net
   * @param is2D
   */
  MEHPForceBalance2(const pylimer_tools::entities::Universe& u,
                    const Network& net,
                    const bool is2D = false);

  /**
   * @brief Construct from old force relaxation network
   *
   * @param net
   * @param is2D
   */
  MEHPForceBalance2(const Network& net, const bool is2D = false)
    : MEHPForceBalance2(
        pylimer_tools::entities::Universe(net.L[0], net.L[1], net.L[2]),
        net,
        is2D) {};

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
   * @param filterEntanglements
   * @param entanglementsAsSprings whether to model the entanglements as merged
   * beads or beads with 1 spring in between
   * @return MEHPForceBalance2
   */
  MEHPForceBalance2(const pylimer_tools::entities::Universe& universe,
                    const size_t nrOfEntanglementsToSample,
                    const double upperCutoff,
                    const double lowerCutoff = 0.,
                    const size_t minimumNrOfEntanglements = 0,
                    const double sameStrandCutoff = 3,
                    const std::string seed = "",
                    const int crossLinkerType = 2,
                    const bool is2D = false,
                    const bool filterEntanglements = true,
                    const bool entanglementsAsSprings = false)
    : MEHPForceBalance2(
        universe,
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
          filterEntanglements),
        crossLinkerType,
        is2D,
        entanglementsAsSprings) {};

#ifdef CEREALIZABLE
  static MEHPForceBalance2 constructFromString(std::string s)
  {
    MEHPForceBalance2 res = MEHPForceBalance2();
    pylimer_tools::utils::deserializeFromString(res, s);
    return res;
  }
#endif

  /**
   * @brief Actually do run the simulation
   *
   * @param simplificationMode the mode for simplifying network structure
   * @param inactiveRemovalCutoff tolerance for removing inactive elements
   * @param solver the solver algorithm to use
   * @param residualReduction the target residual reduction
   * @param maxIterations maximum number of iterations allowed
   */
  void runForceRelaxation(const StructureSimplificationMode simplificationMode =
                            StructureSimplificationMode::NO_SIMPLIFICATION,
                          const double inactiveRemovalCutoff = 1e-6,
                          const SLESolver solver = SLESolver::DEFAULT,
                          const double residualReduction = 1e-12,
                          const int maxIterations = 50000)
  {
    this->runForceRelaxation(
      simplificationMode,
      inactiveRemovalCutoff,
      solver,
      residualReduction,
      maxIterations,
      []() { return false; },
      []() {});
  }

  /**
   * @brief Actually do run the simulation
   *
   * @param simplificationMode the mode for simplifying network structure
   * @param inactiveRemovalCutoff tolerance for removing inactive elements
   * @param solverChoice the solver algorithm to use
   * @param residualReduction the target residual reduction
   * @param maxIterations maximum number of iterations allowed
   * @param shouldInterrupt callback function to check if simulation should be
   * interrupted
   * @param cleanupInterrupt callback function to clean up when interrupted
   */
  void runForceRelaxation(StructureSimplificationMode simplificationMode,
                          double inactiveRemovalCutoff,
                          SLESolver solverChoice,
                          double residualReduction,
                          int maxIterations,
                          const std::function<bool()>& shouldInterrupt,
                          const std::function<void()>& cleanupInterrupt);

  /**
   * @brief Remove springs that exert a stress higher than
   * `this->springBreakingLength`
   *
   * @param net the network to modify
   * @param displacements the current displacements (modified in place)
   * @return size_t the number of springs broken
   */
  size_t breakTooLongStrands(ForceBalance2Network& net,
                             Eigen::VectorXd& displacements) const;

  /**
   * @brief Decide for each spring if it should be removed,
   * remove them, then remove orphaned links
   *
   * @param net the network
   * @param displacements the current displacements
   * @param inactiveRemovalCutoff the cut-off for the activity tolerance
   * criterion
   * @return the number of links removed
   */
  size_t removeInactiveLinks(ForceBalance2Network& net,
                             Eigen::VectorXd& displacements,
                             const double inactiveRemovalCutoff) const;

  /**
   * @brief Decide for each spring if it should be removed
   *
   * @param net the network
   * @param displacements the current displacements
   * @param inactiveRemovalCutoff the cut-off for the activity tolerance
   * criterion
   * @return a vector of true/false for each spring, true if it should be
   * removed
   */
  Eigen::ArrayXb markInactiveSpringsToDelete(
    const ForceBalance2Network& net,
    const Eigen::VectorXd& displacements,
    const double inactiveRemovalCutoff) const;

  /**
   * @brief Delete the springs indicated by `toDelete` from the network
   *
   * @param net the network to be modified
   * @param displacements the current displacements (modified in place)
   * @param springsToDeleteIndices a vector of indices of springs to be removed
   */
  void removeSprings(ForceBalance2Network& net,
                     Eigen::VectorXd& displacements,
                     std::vector<size_t>& springsToDeleteIndices) const;

  /**
   * @brief Replace the two springs traversing a two-functional
   * crosslinkers with a single spring
   *
   * @param net the network to be modified
   * @param displacements the current displacements, will be modified
   * @return the number of bifunctional links replaced by 0-functional links
   */
  size_t unlinkBifunctionalLinks(ForceBalance2Network& net,
                                 Eigen::VectorXd& displacements) const;

  /**
   * @brief Remove links without any springs
   *
   * @param net the network to be modified
   * @param displacements the current displacements, will be modified
   * @return the number of 0-functional links removed
   */
  size_t removeZerofunctionalLinks(ForceBalance2Network& net,
                                   Eigen::VectorXd& displacements) const;

  /**
   * @brief Deform the system to match the specified box
   *
   * @param newBox the box to deform to
   */
  void deformTo(const pylimer_tools::entities::Box& newBox);

  /**
   * @brief Get the universe consisting of crosslinkers only
   *
   * @return pylimer_tools::entities::Universe
   */
  [[nodiscard]] pylimer_tools::entities::Universe getCrosslinkerVerse() const;

  /**
   * @brief Get the default mean bond length
   *
   * @return double the default bond length
   */
  double getDefaultMeanBondLength() const { return this->defaultBondLength; }

  double getVolume() override
  {
    return this->initialConfig.L[0] * this->initialConfig.L[1] *
           this->initialConfig.L[2];
  }

  /**
   * @brief Get the number of nodes in the network
   *
   * @return int the number of nodes
   */
  int getNrOfNodes() const { return this->initialConfig.nrOfNodes; }

  /**
   * @brief Get the number of links in the network
   *
   * @return int the number of links
   */
  int getNrOfLinks() const { return this->initialConfig.nrOfLinks; }

  size_t getNumBonds() override
  {
    return static_cast<size_t>(this->getNrOfStrands());
  }

  size_t getNumExtraBonds() override
  {
    return static_cast<size_t>(
      this->initialConfig.springIsEntanglement.count());
  }

  long int getNumBondsToForm() override { return 0; }

  size_t getNumAtoms() override
  {
    return static_cast<size_t>(this->getNrOfNodes());
  }

  size_t getNumExtraAtoms() override
  {
    return static_cast<size_t>(this->getNrOfLinks() - this->getNrOfNodes());
  }

  /**
   * @brief Get the number of strands in the network
   *
   * @return int the number of strands
   */
  int getNrOfStrands() const { return this->initialConfig.nrOfStrands; }

  /**
   * @brief Get the number of springs in the network
   *
   * @return int the number of springs
   */
  int getNrOfSprings() const { return this->initialConfig.nrOfSprings; }

  /**
   * @brief Get the number of intra-chain slip links
   *
   * @return int the number of intra-chain slip links
   */
  int getNumIntraChainSlipLinks() const;

  /**
   * @brief Get the current displacements
   *
   * @return Eigen::VectorXd the current displacements vector
   */
  Eigen::VectorXd getCurrentDisplacements() const
  {
    return this->currentDisplacements;
  }

  /**
   * @brief Set the current displacements
   *
   * @param displacements the new displacements vector
   */
  void setCurrentDisplacements(const Eigen::VectorXd& displacements)
  {
    this->currentDisplacements = displacements;
  }

  /**
   * @brief Set the spring contour lengths
   *
   * @param springsContourLengths the contour lengths for each spring
   */
  void setSpringContourLengths(const Eigen::VectorXd& springsContourLengths)
  {
    INVALIDARG_EXP_IFN(springsContourLengths.size() ==
                         this->initialConfig.springContourLength.size(),
                       "Contour length must have the correct dimensions.");
    this->initialConfig.springContourLength = springsContourLengths;
  }

  /**
   * @brief Configure the mean bond length
   *
   * @param meanBondLength the mean bond length to set
   */
  void configMeanBondLength(const double meanBondLength)
  {
    this->defaultBondLength = meanBondLength;
  }

  /**
   * @brief Configure the spring constant
   *
   * @param kappa the spring constant value
   */
  void configSpringConstant(const double kappa = 1.0) { this->kappa = kappa; }

  /**
   * @brief Configure the spring breaking distance
   *
   * @param newSpringBreakingForce the force threshold for spring breaking
   */
  void configSpringBreakingDistance(const double newSpringBreakingForce = -1.)
  {
    this->springBreakingLength = newSpringBreakingForce;
  }

  /**
   * @brief Configure whether the network is complete
   *
   * If true, the network is assumed to be complete,
   * i.e., no dangling or free chains have been omitted.
   *
   * @param assumeNetworkIsComplete whether to assume the network is complete
   */
  void configAssumeNetworkIsComplete(const bool assumeNetworkIsComplete = false)
  {
    this->assumeNetworkIsComplete = assumeNetworkIsComplete;
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
    if (this->initialConfig.nrOfNodes == 0) {
      return 0;
    }
    const Eigen::ArrayXb activeNodes = this->findActiveNodes(tolerance);
    return activeNodes.count();
  }

  /**
   * @brief Get the Soluble Weight Fraction
   *
   * @param tolerance the tolerance for considering springs as inactive
   * @return double the soluble weight fraction
   */
  double getSolubleWeightFraction(const double tolerance = 1e-3)
  {
    return this->computeSolubleWeightFraction(
      this->initialConfig, this->currentDisplacements, tolerance);
  }

  /**
   * @brief Get the Dangling Weight Fraction
   *
   * @param tolerance the tolerance for considering springs as inactive
   * @return double the dangling weight fraction
   */
  double getDanglingWeightFraction(const double tolerance = 1e-3)
  {
    return this->computeDanglingWeightFraction(
      this->initialConfig, this->currentDisplacements, tolerance);
  }

  /**
   * @brief Get the Weight Fraction of Active Springs (atoms)
   *
   * @param tolerance the tolerance for considering springs as inactive
   * @return double the active weight fraction
   */
  double getActiveWeightFraction(const double tolerance = 1e-3)
  {
    return this->computeActiveWeightFraction(
      this->initialConfig, this->currentDisplacements, tolerance);
  }

  /**
   * @brief Count the number of atoms that are in any way connected to an
   * active spring
   *
   * @param tolerance the tolerance for considering springs as inactive
   * @return double the number of active clustered atoms
   */
  double countActiveClusteredAtoms(const double tolerance = 1e-3)
  {
    return this->countActiveClusteredAtoms(
      this->initialConfig, this->currentDisplacements, tolerance);
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
    const double tolerance = 1e-3) const;

  /**
   * @brief Compute the weight fraction of non-active springs
   *
   * We go the full route via active and soluble in order to compensate for
   * removed springs and atoms
   *
   * @param net
   * @param u the displacements
   * @param tolerance
   * @return double
   */
  double computeDanglingWeightFraction(ForceBalance2Network& net,
                                       const Eigen::VectorXd& u,
                                       const double tolerance = 1e-3) const;

  /**
   * @brief Infer the number of atoms from the network
   *
   * This is a heuristic that counts the number of springs, links and nodes
   * in the network to estimate the number of atoms.
   *
   * @param net the network to infer from
   * @return double the inferred number of atoms
   */
  static double inferNrOfAtomsFromNetwork(const ForceBalance2Network& net)
  {
    return net.springContourLength.sum() -
           static_cast<double>(net.nrOfSprings) +
           static_cast<double>(
             2 * net.nrOfLinks - net.nrOfNodes // count entanglement links twice
           );
  }

  double inferNrOfAtomsFromNetwork() const
  {
    return this->inferNrOfAtomsFromNetwork(this->initialConfig);
  }

  /**
   * @brief Compute the weight fraction of active springs
   *
   * @param net
   * @param u the displacements
   * @param tolerance
   * @return double
   */
  double computeActiveWeightFraction(ForceBalance2Network& net,
                                     const Eigen::VectorXd& u,
                                     const double tolerance = 1e-3) const;

  /**
   * @brief Find whether strands and nodes are in any way connected to an
   * active spring
   *
   * @param net the network that includes the connectivity
   * @param u the current displacements of the links
   * @param tolerance the tolerance for considering springs as active
   * @return std::pair<Eigen::ArrayXb, Eigen::ArrayXb> indices of springs
   * (first) and links (second) connected in any way to active springs
   */
  std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
  findStrandsAndNodesClusteredToActive(const ForceBalance2Network& net,
                                       const Eigen::VectorXd& u,
                                       const double tolerance = 1e-3) const;

  /**
   * @brief Find whether springs and links that are clustered to active springs
   *
   * @param net the network that includes the connectivity
   * @param u the current displacements of the links
   * @param tolerance the tolerance for considering springs as active
   * @return std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
   */
  std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
  findSpringsAndLinksClusteredToActive(const ForceBalance2Network& net,
                                       const Eigen::VectorXd& u,
                                       const double tolerance = 1e-3) const;

  /**
   * @brief Count the number of atoms that can be considered part of an
   * active cluster, i.e., are somehow connected to an active spring
   *
   * @param net
   * @param u the current displacements of the links
   * @param tolerance
   * @return double
   */
  double countActiveClusteredAtoms(ForceBalance2Network& net,
                                   const Eigen::VectorXd& u,
                                   const double tolerance = 1e-3) const;

  /**
   * @brief Compute the weight fraction of springs connected to active
   * springs (any depth)
   *
   * @param net
   * @param u the current displacements of the links
   * @param tolerance
   * @return double
   */
  double computeSolubleWeightFraction(ForceBalance2Network& net,
                                      const Eigen::VectorXd& u,
                                      const double tolerance = 1e-3) const;

  /**
   * @brief Get the indices of active Nodes
   *
   * @param net
   * @param u the current displacements of the links
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return std::vector<long int> the atom ids
   */
  std::vector<int> getIndicesOfActiveNodes(const ForceBalance2Network& net,
                                           const Eigen::VectorXd& u,
                                           const double tolerance = 1e-3) const;

  /**
   * @brief Get the Ids of active Nodes
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return std::vector<long int> the atom ids
   */
  std::vector<long int> getAtomIdsOfActiveNodes(
    const double tolerance = 1e-3) const;

  /**
   *
   * @return a vector with the sum of the norms of all the springs per strand
   */
  std::vector<double> getOverallSpringLengths() const;

  /**
   * @brief Get the current spring distances
   *
   * @return Eigen::VectorXd the current spring distances
   */
  Eigen::VectorXd getCurrentSpringDistances() const;

  /**
   * @brief Get the current spring lengths
   *
   * @return std::vector<double> the current spring lengths
   */
  std::vector<double> getCurrentSpringLengths() const;

  /**
   * @brief Get the Nr Of Active Springs connected to each node
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return Eigen::VectorXi
   */
  Eigen::VectorXi getNrOfActiveStrandsConnected(
    const double tolerance = 1e-3) const;

  /**
   * @brief Get the Nr Of Active Springs connected to each node
   *
   * @param tolerance springs under a certain length are considered inactive
   * @return Eigen::VectorXi
   */
  Eigen::VectorXi getNrOfActiveSpringsConnected(
    const double tolerance = 1e-3) const;

  /**
   * @brief Get the Nr Of Active Springs object
   *
   * @param tolerance springs under a certain length are considered inactive
   * @return int the number of active strands
   */
  int getNrOfActiveStrands(const double tolerance = 1e-3) const
  {
    return this->countNrOfActiveStrands(tolerance);
  }

  /**
   * @brief Get the number of active strands in a specific direction
   *
   * @param dir the direction (0=x, 1=y, 2=z)
   * @param tolerance springs under a certain length are considered inactive
   * @return int the number of active strands in the specified direction
   */
  int getNrOfActiveStrandsInDir(const int dir,
                                const double tolerance = 1e-3) const
  {
    return this->countNrOfActiveStrandsInDir(dir, tolerance);
  }

  /**
   * @brief Get the Nr Of Active Springs object
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return int the number of active springs
   */
  int getNrOfActiveSprings(const double tolerance = 1e-3) const
  {
    return this->countNrOfActiveSprings(tolerance);
  }

  /**
   * @brief Get the Average Spring Length at the current step
   *
   * @return double the average strand length
   */
  double getAverageStrandLength() const;

  /**
   * @brief Get the stress tensor
   *
   * @return Eigen::Matrix3d the stress tensor
   */
  Eigen::Matrix3d getStressTensor() override;

  /**
   * @brief Get the stress tensor based on links
   *
   * @param crosslinksOnly whether to only consider crosslinks
   * @return Eigen::Matrix3d the stress tensor
   */
  Eigen::Matrix3d getStressTensorLinkBased(
    const bool crosslinksOnly = false) const;

  /**
   * @brief Get the Pressure
   *
   * @return double the pressure
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

  /**
   * @brief Get the gamma factor override
   *
   * @return double the gamma factor
   */
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
   * @return int the number of iterations
   */
  int getNrOfIterations() const { return this->nrOfStepsDone; }

  /**
   * @brief Get the exit/stop reason for the simulation
   *
   * @return ExitReason the reason for simulation stop
   */
  ExitReason getExitReason() const { return this->exitReason; }

  /**
   * @brief Compute one spring length
   *
   * @param net
   * @param u the current displacements of the links
   * @param springIdx
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d evaluateSpringVector(const ForceBalance2Network& net,
                                       const Eigen::VectorXd& u,
                                       const size_t springIdx) const
  {
    return this->evaluateSpringVector(net, u, springIdx, this->is2D);
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
  static Eigen::Vector3d evaluateSpringVector(const ForceBalance2Network& net,
                                              const Eigen::VectorXd& u,
                                              const size_t springIdx,
                                              const bool is2d);

  /**
   *
   * @param net the network
   * @param u the current displacements
   * @param is2D
   * @return the vectors of the springs
   */
  static Eigen::VectorXd evaluateSpringVectors(const ForceBalance2Network& net,
                                               const Eigen::VectorXd& u,
                                               const bool is2D);

  Eigen::VectorXd evaluateSpringVectors(const ForceBalance2Network& net,
                                        const Eigen::VectorXd& u) const
  {
    return this->evaluateSpringVectors(net, u, this->is2D);
  };

  /**
   *
   * @param net the network
   * @param u the current displacements
   * @return the strand lengths (norm of the strand end-to-end vector)
   */
  Eigen::VectorXd evaluateStrandLengths(const ForceBalance2Network& net,
                                        Eigen::VectorXd u) const;

  /**
   *
   * @param net the network
   * @param u the current displacements
   * @return the strand end-to-end vectors
   */
  Eigen::Vector3d evaluateStrandVector(const ForceBalance2Network& net,
                                       const Eigen::VectorXd u,
                                       const size_t strandIdx) const;
  /**
   *
   * @param net the network
   * @param u the current displacements
   * @return the strand end-to-end vectors
   */
  Eigen::VectorXd evaluateStrandVectors(const ForceBalance2Network& net,
                                        const Eigen::VectorXd u) const;

  /**
   * @brief Sum the partitions up to a given link in a spring
   *
   * @param net
   * @param springPartition
   * @param springIdx
   * @param targetLink
   * @return double
   */
  static double sumToTotalFraction(const ForceBalance2Network& net,
                                   Eigen::VectorXd springPartition,
                                   size_t springIdx,
                                   size_t targetLink);

  /**
   * @brief Query the link index of the other end of a spring
   *
   * @param net the network with the connectivity
   * @param springIdx the spring to get the ends of
   * @param linkIdx the link which is not searched
   * @return the link index of the other end of the spring
   */
  static size_t getOtherLinkOfSpring(const ForceBalance2Network& net,
                                     const size_t springIdx,
                                     const size_t linkIdx)
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
  static Eigen::Vector3d getSpringBoxOffset(const ForceBalance2Network& net,
                                            const size_t partialSpringIdx)
  {
    return net.springBoxOffset.segment(3 * partialSpringIdx, 3);
  }

  /**
   * @brief Query the box offset for a specific spring in a specific
   * direction
   *
   * @param net
   * @param springIdx
   * @param linkIdx
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d getSpringBoxOffsetTo(const ForceBalance2Network& net,
                                       const size_t springIdx,
                                       const size_t linkIdx) const
  {
    return (net.springIndexA(springIdx) == linkIdx)
             ? (-1. * this->getSpringBoxOffset(net, springIdx))
             : (this->getSpringBoxOffset(net, springIdx));
  }

  Eigen::Vector3d getSpringBoxOffsetFrom(const ForceBalance2Network& net,
                                         const size_t springIdx,
                                         const size_t linkIdx) const
  {
    return -1. * this->getSpringBoxOffsetTo(net, springIdx, linkIdx);
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
  Eigen::Vector3d evaluateSpringVectorTo(const ForceBalance2Network& net,
                                         const Eigen::VectorXd& u,
                                         const size_t springIdx,
                                         const size_t linkIdx) const
  {
    return this->evaluateSpringVectorTo(net, u, springIdx, linkIdx, this->is2D);
  }

  Eigen::Vector3d evaluateSpringVectorTo(const ForceBalance2Network& net,
                                         const Eigen::VectorXd& u,
                                         const size_t springIdx,
                                         const size_t linkIdx,
                                         const bool is2d) const
  {
    assert(this->isPartOfSpring(net, linkIdx, springIdx));

    const Eigen::Vector3d dist =
      this->evaluateSpringVector(net, u, springIdx, is2d);

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
  Eigen::Vector3d evaluateSpringVectorFrom(const ForceBalance2Network& net,
                                           const Eigen::VectorXd& u,
                                           const size_t springIdx,
                                           const size_t linkIdx) const
  {
    return this->evaluateSpringVectorFrom(
      net, u, springIdx, linkIdx, this->is2D);
  }

  Eigen::Vector3d evaluateSpringVectorFrom(const ForceBalance2Network& net,
                                           const Eigen::VectorXd& u,
                                           const size_t springIdx,
                                           const size_t linkIdx,
                                           const bool is2d) const
  {
    return -1. * this->evaluateSpringVectorTo(net, u, springIdx, linkIdx, is2d);
  }

  ForceBalance2Network getNetwork() { return this->initialConfig; }

  /**
   * @brief Get the Weighted Partial Spring Length for one partial spring
   *
   * @return double
   */
  double getWeightedSpringLength(const ForceBalance2Network& net,
                                 const Eigen::VectorXd& u,
                                 size_t partialSpringIdx) const;

  /**
   * @brief Get the Weighted Partial Spring Lengths
   *
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getWeightedSpringLengths() const
  {
    Eigen::VectorXd weightedLengths =
      Eigen::VectorXd(this->initialConfig.nrOfSprings);
    for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
      weightedLengths(i) = this->getWeightedSpringLength(
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
   * @return std::vector<size_t>
   */
  std::vector<size_t> getPartialSpringIndicesOfLink(
    const ForceBalance2Network& net,
    const size_t linkIdx) const;

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
   * @param net
   * @param u
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
   * @brief Evaluate the stress on a particular cross- or slip-link
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
   * @return the displacements afterward
   */
  Eigen::VectorXd inspectDisplacementToMeanPositionUpdate(
    const size_t linkIdx) const
  {
    Eigen::VectorXd displacements = this->currentDisplacements;
    this->displaceToMeanPosition(this->initialConfig, displacements, linkIdx);
    return displacements;
  };

  /**
   * @brief Displace all links to the mean of all connected neighbours
   *
   * @param net the force balance network
   * @param u the current displacements, wherein the resulting coordinates
   * shall be stored
   * @return double, the distance (squared norm) displaced
   */
  double displaceToMeanPosition(const ForceBalance2Network& net,
                                Eigen::VectorXd& u) const;

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
  static Eigen::VectorXd assembleOneOverSpringPartition(
    const ForceBalance2Network& net);

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
   *
   * @param net the network for which the entanglement springs are to be listed
   * @param strandIdx the strand for which the entanglement springs are to be
   * listed
   * @return the indices of the entanglement springs associated with a strand
   */
  std::vector<size_t> getEntanglementSpringsAtStrand(
    const ForceBalance2Network& net,
    size_t strandIdx) const;

  /**
   * @brief Get the Link Indices of all neighbours of a specified link
   *
   * @param net
   * @param linkIdx
   * @return std::vector<size_t>
   */
  static std::vector<size_t> getNeighbourLinkIndices(
    const ForceBalance2Network& net,
    const size_t linkIdx);

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

  Eigen::VectorXd getBondLengths() override;

  Eigen::VectorXd getInitialCoordinates();

  Eigen::VectorXd getCoordinates() override;

  double getTemperature() override;

  size_t getNumParticles() override;

  void debugAtomVicinity(const size_t atomId) const;

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
   * @brief Finish initializing some member properties
   *
   */
  void completeInitialization();

  /**
   * @brief Merge two springs around a given link
   *
   * @param net
   * @param removedSpringIdx
   * @param keptSpringIdx
   * @param linkToReduce the index of the link to combine the springs around
   */
  void mergeSpringsWithoutRemoval(ForceBalance2Network& net,
                                  const Eigen::VectorXd& u,
                                  const size_t removedSpringIdx,
                                  const size_t keptSpringIdx,
                                  const size_t linkToReduce) const;

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
    const auto stressTensor = this->evaluateStressTensor(net, u);
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
   * @param crosslinksOnly whether to only consider crosslinks
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
   * @brief Compute the force acting on a slip- or crosslink
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
   * @brief Compute the stress acting on a slip- or crosslink
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
  [[nodiscard]] int countNrOfActiveStrands(const ForceBalance2Network& net,
                                           const Eigen::VectorXd& u,
                                           const double tolerance = 1e-3) const
  {
    return (this->findActiveStrands(net, u, tolerance)).count();
  }

  [[nodiscard]] int countNrOfActiveStrands(const double tolerance = 1e-3) const
  {
    return (this->findActiveStrands(tolerance) == true).count();
  }

  [[nodiscard]] int countNrOfActiveStrandsInDir(
    const int dir,
    const double tolerance = 1e-3) const
  {
    return (this->findActiveStrandsInDir(dir, tolerance) == true).count();
  }

  [[nodiscard]] int countNrOfActiveSprings(const double tolerance = 1e-3) const
  {
    return (this->findActiveSprings(tolerance) == true).count();
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
  Eigen::ArrayXb findActiveStrands(const ForceBalance2Network& net,
                                   const Eigen::VectorXd& u,
                                   const double tolerance = 1e-3) const
  {
    Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net.nrOfStrands, false);
    Eigen::ArrayXb activeSprings = this->findActiveSprings(net, u, tolerance);
    const size_t nActiveSprings = activeSprings.count();

    for (size_t strandIdx = 0; strandIdx < net.nrOfStrands; ++strandIdx) {
      // without any springs, the strand is inactive
      if (net.springIndicesOfStrand[strandIdx].empty()) {
        continue;
      }
      // with one spring, the strand is active if the spring is active
      if (net.springIndicesOfStrand[strandIdx].size() == 1) {
        result[strandIdx] =
          activeSprings[net.springIndicesOfStrand[strandIdx][0]];
        continue;
      }
      // with more springs, the strand is certainly active,
      // if the first and last spring are active
      // (springs in between may be active due to the other involved strands)
      assert(
        !this->isLoopingSpring(net, net.springIndicesOfStrand[strandIdx][0]));
      bool isActive0 = activeSprings[net.springIndicesOfStrand[strandIdx][0]];
      assert(!this->isLoopingSpring(
        net, net.springIndicesOfStrand[strandIdx].back()));
      bool isActiveN =
        activeSprings[net.springIndicesOfStrand[strandIdx].back()];

      const bool isAnyActive =
        std::ranges::any_of(net.springIndicesOfStrand[strandIdx],
                            [&activeSprings](const long int springIdx) {
                              return activeSprings[springIdx];
                            });
      if (!isAnyActive) {
        continue;
      }

      // however, there is one case where the
      // strand would not yet be marked as active, even though it is:
      // when one end is inactive,
      // because it is bifunctional and in a sandwich with the same entanglement
      // link
      if (!isActive0) {
        const size_t crossLinkIdx0 = net.linkIndicesOfStrand[strandIdx][0];
        const size_t entanglementLinkIdx0 =
          net.linkIndicesOfStrand[strandIdx][1];
        assert(!net.linkIsEntanglement[crossLinkIdx0]);
        assert(net.linkIsEntanglement[entanglementLinkIdx0]);
        for (const size_t strandOfXlink :
             net.strandIndicesOfLink[crossLinkIdx0]) {
          if (strandOfXlink == strandIdx) {
            continue;
          }
          if (net.linkIndicesOfStrand[strandOfXlink][0] == crossLinkIdx0 &&
              net.linkIndicesOfStrand[strandOfXlink][1] ==
                entanglementLinkIdx0) {
            isActive0 = true;
            // break;
          }
          if (net.linkIndicesOfStrand[strandOfXlink].back() == crossLinkIdx0 &&
              net.linkIndicesOfStrand
                  [strandOfXlink]
                  [net.linkIndicesOfStrand[strandOfXlink].size() - 2] ==
                entanglementLinkIdx0) {
            isActive0 = true;
            // break;
          }
        }
      }

      if (!isActiveN) {
        const size_t crossLinkIdxN = net.linkIndicesOfStrand[strandIdx].back();
        const size_t entanglementLinkIdxN =
          net
            .linkIndicesOfStrand[strandIdx]
                                [net.linkIndicesOfStrand[strandIdx].size() - 2];
        assert(!net.linkIsEntanglement[crossLinkIdxN]);
        assert(net.linkIsEntanglement[entanglementLinkIdxN]);
        for (const size_t strandOfXlink :
             net.strandIndicesOfLink[crossLinkIdxN]) {
          if (strandOfXlink == strandIdx) {
            continue;
          }
          if (net.linkIndicesOfStrand[strandOfXlink][0] == crossLinkIdxN &&
              net.linkIndicesOfStrand[strandOfXlink][1] ==
                entanglementLinkIdxN) {
            isActiveN = true;
            // break;
          }
          if (net.linkIndicesOfStrand[strandOfXlink].back() == crossLinkIdxN &&
              net.linkIndicesOfStrand
                  [strandOfXlink]
                  [net.linkIndicesOfStrand[strandOfXlink].size() - 2] ==
                entanglementLinkIdxN) {
            isActiveN = true;
            // break;
          }
        }
      }
      // we can, however, anticipate the activeness of dangling links
      // and mark the strand as inactive in that case
      if (net.linkIndicesOfStrand[strandIdx][0] !=
          net.linkIndicesOfStrand[strandIdx].back()) {
        assert(!net.linkIsEntanglement[net.linkIndicesOfStrand[strandIdx][0]]);
        if (net.strandIndicesOfLink[net.linkIndicesOfStrand[strandIdx][0]]
              .size() == 1) {
          isActive0 = false;
        }
        assert(
          !net.linkIsEntanglement[net.linkIndicesOfStrand[strandIdx].back()]);
        if (net.strandIndicesOfLink[net.linkIndicesOfStrand[strandIdx].back()]
              .size() == 1) {
          isActiveN = false;
        }
      }

      result[strandIdx] = isActive0 && isActiveN;
    }

#ifndef NDEBUG
    const size_t nActiveStrands = result.count();
    assert(nActiveStrands <= nActiveSprings);
#endif

    return result;
  }

  Eigen::ArrayXb findActiveStrands(const double tolerance = 1e-3) const
  {
    return this->findActiveStrands(
      this->initialConfig, this->currentDisplacements, tolerance);
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
  Eigen::ArrayXb findActiveStrandsInDir(const ForceBalance2Network& net,
                                        const Eigen::VectorXd& u,
                                        const int dir,
                                        const double tolerance = 1e-3) const
  {
    INVALIDARG_EXP_IFN(dir >= 0 && dir < 3, "Invalid direction");
    Eigen::VectorXd partialSpringVectors =
      this->evaluateSpringVectors(net, u, this->is2D);
    Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net.nrOfStrands, false);

    for (size_t i = 0; i < net.nrOfSprings; ++i) {
      result[net.strandIndexOfSpring[i]] =
        result[net.strandIndexOfSpring[i]] ||
        !this->distanceIsWithinTolerance(
          Eigen::Vector3d(partialSpringVectors[3 * i + dir], 0, 0),
          tolerance,
          net.springContourLength[net.strandIndexOfSpring[i]]);
    }

    return result;
  }

  Eigen::ArrayXb findActiveStrandsInDir(const int dir,
                                        const double tolerance = 1e-3) const
  {
    return this->findActiveStrandsInDir(
      this->initialConfig, this->currentDisplacements, dir, tolerance);
  }

  /**
   *
   * @param tolerance the tolerance to decide whether a spring is considered
   * active
   * @return for each node whether the node is connected to at least one active
   * spring
   */
  Eigen::ArrayXb findActiveLinks(const double tolerance = 1e-3) const
  {
    Eigen::ArrayXb activeStrands = this->findActiveStrands(tolerance);
    assert(activeStrands.size() == this->initialConfig.nrOfStrands);

    Eigen::ArrayXb activeLinks =
      Eigen::ArrayXb::Constant(this->initialConfig.nrOfLinks, false);
    for (size_t i = 0; i < activeStrands.size(); ++i) {
      if (activeStrands[i]) {
        activeLinks[this->initialConfig.linkIndicesOfStrand[i][0]] = true;
        activeLinks[pylimer_tools::utils::last(
          this->initialConfig.linkIndicesOfStrand[i])] = true;
      }
    }

    return activeLinks;
  }

  /**
   *
   * @param tolerance the tolerance to decide whether a spring is considered
   * active
   * @return for each node whether the node is connected to at least one active
   * spring
   */
  Eigen::ArrayXb findActiveNodes(const double tolerance = 1e-3) const
  {
    Eigen::ArrayXb activeLinks = this->findActiveLinks(tolerance);
    Eigen::ArrayXb activeNodes =
      Eigen::ArrayXb::Constant(this->initialConfig.nrOfNodes, false);
    size_t nodeIdx = 0;
    for (size_t i = 0; i < this->initialConfig.nrOfLinks; ++i) {
      if (!this->initialConfig.linkIsEntanglement[i]) {
        activeNodes[nodeIdx] = activeLinks[i];
        ++nodeIdx;
      }
    }

    return activeNodes;
  }

  /**
   *
   * @param net the network
   * @param u the current displacements
   * @param tolerance the tolerance to decide whether a spring is considered
   * active
   * @return for each spring whether it may be considered active (tolerance
   * criterion)
   */
  Eigen::ArrayXb findActiveSprings(const ForceBalance2Network& net,
                                   const Eigen::VectorXd& u,
                                   const double tolerance = 1e-3) const
  {
    Eigen::VectorXd springVectors =
      this->evaluateSpringVectors(net, u, this->is2D);
    Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net.nrOfSprings, false);

    for (size_t i = 0; i < net.nrOfSprings; ++i) {
      result[i] = !this->distanceIsWithinTolerance(
        springVectors.segment(3 * i, 3), tolerance, net.springContourLength[i]);
    }

    return result;
  }

  Eigen::ArrayXb findActiveSprings(const double tolerance = 1e-3) const
  {
    return this->findActiveSprings(
      this->initialConfig, this->currentDisplacements, tolerance);
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

  static bool isPartOfSpring(const ForceBalance2Network& net,
                             const size_t linkIdx,
                             const size_t springIdx)
  {
    return (net.springIndexA[springIdx] == linkIdx) ||
           (net.springIndexB[springIdx] == linkIdx);
  }

  static bool isLoopingSpring(const ForceBalance2Network& net,
                              const size_t springIdx)
  {
    return (net.springIndexA[springIdx] == net.springIndexB[springIdx]);
  }

  static bool isLoopingStrand(const ForceBalance2Network& net,
                              const size_t strandIdx)
  {
    if (net.linkIndicesOfStrand[strandIdx].size() == 0) {
      return false;
    }
    assert(net.linkIndicesOfStrand[strandIdx].size() >= 2);
    return (net.linkIndicesOfStrand[strandIdx][0] ==
            net.linkIndicesOfStrand[strandIdx].back());
  }

  static double getDenominatorOfPartialSpring(const ForceBalance2Network& net,
                                              const size_t partialSpringIdx);
};
} // namespace
