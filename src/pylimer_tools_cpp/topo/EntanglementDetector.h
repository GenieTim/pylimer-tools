#ifndef ENTANGLEMENT_DETECTOR_H
#define ENTANGLEMENT_DETECTOR_H

#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../sim/MEHPUtilityStructures.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <algorithm>
#include <string>
#include <vector>

namespace pylimer_tools::topo::entanglement_detection {

/**
 * @brief Filtering mode for entanglement sampling.
 */
enum class EntanglementFilteringMode
{
  NONE = 0,          ///< No filtering: all entanglements retained
  OMIT_INACTIVE = 1, ///< Melt density mode: sample from all strands, omit from
                     ///< output if touching dangling/sol
  PRUNE_UPFRONT = 2  ///< Active-only mode: topologically peel leaves upfront,
                     ///< sample strictly between active strands
};

/**
 * @brief Convert string to EntanglementFilteringMode
 */
inline EntanglementFilteringMode
parseEntanglementFilteringMode(const std::string& str)
{
  std::string s = str;
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  if (s == "prune_upfront" || s == "prune" || s == "active_only") {
    return EntanglementFilteringMode::PRUNE_UPFRONT;
  }
  if (s == "none" || s == "unfiltered") {
    return EntanglementFilteringMode::NONE;
  }
  if (s == "omit_inactive" || s == "omit" || s == "filtered") {
    return EntanglementFilteringMode::OMIT_INACTIVE;
  }
  INVALIDARG_EXP_IFN(
    false,
    "Unknown entanglement filtering mode: '" + str +
      "'. Expected 'none', 'omit_inactive', or 'prune_upfront'.");
  return EntanglementFilteringMode::OMIT_INACTIVE;
}

inline std::string
to_string(const EntanglementFilteringMode mode)
{
  switch (mode) {
    case EntanglementFilteringMode::NONE:
      return "none";
    case EntanglementFilteringMode::PRUNE_UPFRONT:
      return "prune_upfront";
    case EntanglementFilteringMode::OMIT_INACTIVE:
    default:
      return "omit_inactive";
  }
}

struct TopologicalLeafPeelingResult
{
  std::vector<bool> strandIsIgnored;
  std::vector<bool> crosslinkerIsIgnored;
};

/**
 * @brief Perform iterative 2-core topological leaf peeling on a
 * CrosslinkerUniverse.
 *
 * @param crosslinkerUniverse The crosslinker & strand representation
 * @return TopologicalLeafPeelingResult containing boolean flags for both
 * strands and crosslinkers
 */
TopologicalLeafPeelingResult
peelTopologicalLeavesDetailed(
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse);

/**
 * @brief Perform iterative 2-core topological leaf peeling on a
 * CrosslinkerUniverse.
 *
 * @param crosslinkerUniverse The crosslinker & strand representation
 * @return std::vector<bool> Boolean vector indicating which strands are peeled
 * (ignored)
 */
std::vector<bool>
peelTopologicalLeaves(
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse);

struct AtomPairEntanglements
{
  /**
   * pairs of atom ids
   */
  std::vector<std::pair<size_t, size_t>> pairsOfAtoms = {};
  /**
   * For each atom, the index of the pair it is associated with, or -1 if
   * none
   */
  std::vector<long int> pairOfAtom = {};

  /**
   * @brief Constructor for AtomPairEntanglements
   *
   * @param pairs Vector of atom pairs (default: empty)
   * @param pairIndices Vector of pair indices for each atom (default: empty)
   */
  AtomPairEntanglements(
    const std::vector<std::pair<size_t, size_t>>& pairs = {},
    const std::vector<long int>& pairIndices = {})
    : pairsOfAtoms(pairs)
    , pairOfAtom(pairIndices)
  {
  }
};

struct StrandPairEntanglements
{
  /**
   * First strand index for each entanglement
   */
  std::vector<size_t> strandIdx1 = {};
  /**
   * Second strand index for each entanglement
   */
  std::vector<size_t> strandIdx2 = {};
  /**
   * Fractional contour position in (0, 1) along strand 1 for each entanglement
   */
  std::vector<double> alpha1 = {};
  /**
   * Fractional contour position in (0, 1) along strand 2 for each entanglement
   */
  std::vector<double> alpha2 = {};
  /**
   * 0-based bead index in strand 1 for each entanglement
   */
  std::vector<size_t> beadIdx1 = {};
  /**
   * 0-based bead index in strand 2 for each entanglement
   */
  std::vector<size_t> beadIdx2 = {};
  /**
   * Sampled x coordinate of each entanglement (midpoint of the two sampled
   * beads in a consistent periodic frame)
   */
  Eigen::VectorXd x;
  /**
   * Sampled y coordinate of each entanglement
   */
  Eigen::VectorXd y;
  /**
   * Sampled z coordinate of each entanglement
   */
  Eigen::VectorXd z;
  /**
   * Sampled 3D coordinates of the bead on strand 1 (stored as 3*N flat vector)
   */
  Eigen::VectorXd position1;
  /**
   * Sampled 3D coordinates of the bead on strand 2 (stored as 3*N flat vector,
   * in the minimum-image frame of position1)
   */
  Eigen::VectorXd position2;

  size_t size() const { return this->strandIdx1.size(); }

  /**
   * @brief Convert strand-level entanglements into AtomPairEntanglements
   * using the deterministic atom/vertex ordering of MCUniverseGenerator.
   *
   * @param crosslinkerUniverse the crosslinker universe
   * @return AtomPairEntanglements
   */
  AtomPairEntanglements toAtomPairEntanglements(
    const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse)
    const;
};

/**
 * @brief Analytically and probabilistically sample spatial entanglements
 * directly on a CrosslinkerUniverse using Brownian bridge / random walk
 * conditional Gaussian distributions and topological filtering, bypassing
 * explicit 3D bead coordinates and Universe / igraph_t construction.
 *
 * @param box the simulation box
 * @param crosslinkerUniverse the crosslinker & strand representation
 * @param nrOfEntanglementsToSample the maximum number of entanglements to
 * sample
 * @param upperCutoff the upper distance cutoff within which two beads can
 * entangle
 * @param lowerCutoff the lower distance cutoff
 * @param minimumNrOfEntanglements the minimum number of entanglements to sample
 * @param sameStrandCutoff the minimum number of bonds between two entangled
 * beads on the same strand
 * @param seed the random seed
 * @param filteringMode filtering mode (NONE, OMIT_INACTIVE, or PRUNE_UPFRONT)
 * @return StrandPairEntanglements
 */
StrandPairEntanglements
randomlyFindEntanglementsAnalytical(
  const pylimer_tools::entities::Box& box,
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse,
  const size_t nrOfEntanglementsToSample,
  const double upperCutoff,
  const double lowerCutoff = 0.,
  const size_t minimumNrOfEntanglements = 0,
  const double sameStrandCutoff = 3,
  const std::string& seed = "",
  EntanglementFilteringMode filteringMode =
    EntanglementFilteringMode::OMIT_INACTIVE);

inline StrandPairEntanglements
randomlyFindEntanglementsAnalytical(
  const pylimer_tools::entities::Box& box,
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse,
  const size_t nrOfEntanglementsToSample,
  const double upperCutoff,
  const double lowerCutoff,
  const size_t minimumNrOfEntanglements,
  const double sameStrandCutoff,
  const std::string& seed,
  const bool filterDanglingAndSoluble)
{
  return randomlyFindEntanglementsAnalytical(
    box,
    crosslinkerUniverse,
    nrOfEntanglementsToSample,
    upperCutoff,
    lowerCutoff,
    minimumNrOfEntanglements,
    sameStrandCutoff,
    seed,
    filterDanglingAndSoluble ? EntanglementFilteringMode::OMIT_INACTIVE
                             : EntanglementFilteringMode::NONE);
}

/**
 * @brief Sample spatial entanglements on a CrosslinkerUniverse using
 * per-strand bounding boxes and on-demand Brownian bridge bead resampling.
 * Only strand seeds and bounding boxes are stored, achieving the O(N_strands)
 * memory scaling of the analytical method while sampling from the exact
 * discrete multi-bead chain conformations.
 *
 * @param box the simulation box
 * @param crosslinkerUniverse the crosslinker & strand representation
 * @param nrOfEntanglementsToSample the maximum number of entanglements to
 * sample
 * @param upperCutoff the upper distance cutoff within which two beads can
 * entangle
 * @param lowerCutoff the lower distance cutoff
 * @param minimumNrOfEntanglements the minimum number of entanglements to sample
 * @param sameStrandCutoff the minimum number of bonds between two entangled
 * beads on the same strand
 * @param seed the random seed
 * @param filteringMode filtering mode (NONE, OMIT_INACTIVE, or PRUNE_UPFRONT)
 * @param nMcSteps number of Monte Carlo equilibration steps for chain
 * resampling (default: 0)
 * @return StrandPairEntanglements
 */
StrandPairEntanglements
randomlyFindEntanglementsOnDemand(
  const pylimer_tools::entities::Box& box,
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse,
  const size_t nrOfEntanglementsToSample,
  const double upperCutoff,
  const double lowerCutoff = 0.,
  const size_t minimumNrOfEntanglements = 0,
  const double sameStrandCutoff = 3,
  const std::string& seed = "",
  EntanglementFilteringMode filteringMode =
    EntanglementFilteringMode::OMIT_INACTIVE,
  int nMcSteps = 0);

inline StrandPairEntanglements
randomlyFindEntanglementsOnDemand(
  const pylimer_tools::entities::Box& box,
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse,
  const size_t nrOfEntanglementsToSample,
  const double upperCutoff,
  const double lowerCutoff,
  const size_t minimumNrOfEntanglements,
  const double sameStrandCutoff,
  const std::string& seed,
  const bool filterDanglingAndSoluble,
  const int nMcSteps = 0)
{
  return randomlyFindEntanglementsOnDemand(
    box,
    crosslinkerUniverse,
    nrOfEntanglementsToSample,
    upperCutoff,
    lowerCutoff,
    minimumNrOfEntanglements,
    sameStrandCutoff,
    seed,
    filterDanglingAndSoluble ? EntanglementFilteringMode::OMIT_INACTIVE
                             : EntanglementFilteringMode::NONE,
    nMcSteps);
}

/**
 * @brief Randomly find pairs of atoms that are close together and could be
 * entanglements
 *
 * @param universe the universe with all the atoms etc.
 * @param nrOfEntanglementsToSample the nr of entanglements to find
 * @param upperCutoff the cut-off distance within which two atoms are
 * considered close enough
 * @param lowerCutoff
 * @param minimumNrOfEntanglements the minimum number of entanglements to
 * find
 * @param sameStrandCutoff the number of beads required between two atoms of
 * the same strand
 * @param seed the random seed
 * @param crossLinkerType the type of the crosslink atoms
 * @param ignoreCrosslinks
 * @param filteringMode filtering mode (NONE, OMIT_INACTIVE, or PRUNE_UPFRONT)
 * @return AtomPairEntanglements
 */
AtomPairEntanglements
randomlyFindEntanglements(const pylimer_tools::entities::Universe& universe,
                          const size_t nrOfEntanglementsToSample,
                          const double upperCutoff,
                          const double lowerCutoff = 0.,
                          const size_t minimumNrOfEntanglements = 0,
                          const double sameStrandCutoff = 3,
                          const std::string& seed = "",
                          int crossLinkerType = 2,
                          bool ignoreCrosslinks = true,
                          EntanglementFilteringMode filteringMode =
                            EntanglementFilteringMode::OMIT_INACTIVE);

inline AtomPairEntanglements
randomlyFindEntanglements(const pylimer_tools::entities::Universe& universe,
                          const size_t nrOfEntanglementsToSample,
                          const double upperCutoff,
                          const double lowerCutoff,
                          const size_t minimumNrOfEntanglements,
                          const double sameStrandCutoff,
                          const std::string& seed,
                          const int crossLinkerType,
                          const bool ignoreCrosslinks,
                          const bool filterDanglingAndSoluble)
{
  return randomlyFindEntanglements(universe,
                                   nrOfEntanglementsToSample,
                                   upperCutoff,
                                   lowerCutoff,
                                   minimumNrOfEntanglements,
                                   sameStrandCutoff,
                                   seed,
                                   crossLinkerType,
                                   ignoreCrosslinks,
                                   filterDanglingAndSoluble
                                     ? EntanglementFilteringMode::OMIT_INACTIVE
                                     : EntanglementFilteringMode::NONE);
}

/**
 * @brief Legacy alias for randomlyFindEntanglements using modern Eigen cell
 * lists and topological leaf peeling.
 */
inline AtomPairEntanglements
randomlyFindEntanglementsV2(const pylimer_tools::entities::Universe& universe,
                            const size_t nrOfEntanglementsToSample,
                            const double upperCutoff,
                            const double lowerCutoff = 0.,
                            const size_t minimumNrOfEntanglements = 0,
                            const double sameStrandCutoff = 3,
                            const std::string& seed = "",
                            int crossLinkerType = 2,
                            bool ignoreCrosslinks = true,
                            EntanglementFilteringMode filteringMode =
                              EntanglementFilteringMode::OMIT_INACTIVE)
{
  return randomlyFindEntanglements(universe,
                                   nrOfEntanglementsToSample,
                                   upperCutoff,
                                   lowerCutoff,
                                   minimumNrOfEntanglements,
                                   sameStrandCutoff,
                                   seed,
                                   crossLinkerType,
                                   ignoreCrosslinks,
                                   filteringMode);
}

inline AtomPairEntanglements
randomlyFindEntanglementsV2(const pylimer_tools::entities::Universe& universe,
                            const size_t nrOfEntanglementsToSample,
                            const double upperCutoff,
                            const double lowerCutoff,
                            const size_t minimumNrOfEntanglements,
                            const double sameStrandCutoff,
                            const std::string& seed,
                            const int crossLinkerType,
                            const bool ignoreCrosslinks,
                            const bool filterDanglingAndSoluble)
{
  return randomlyFindEntanglements(universe,
                                   nrOfEntanglementsToSample,
                                   upperCutoff,
                                   lowerCutoff,
                                   minimumNrOfEntanglements,
                                   sameStrandCutoff,
                                   seed,
                                   crossLinkerType,
                                   ignoreCrosslinks,
                                   filterDanglingAndSoluble
                                     ? EntanglementFilteringMode::OMIT_INACTIVE
                                     : EntanglementFilteringMode::NONE);
}

}

#endif
