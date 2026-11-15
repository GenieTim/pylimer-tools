#pragma once

#include "../entities/Box.h"
#include "../sim/MEHPUtilityStructures.h"
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

namespace pylimer_tools::utils {

/**
 * @brief Topological category of a strand in a CrosslinkerUniverse.
 */
enum class StrandCategory
{
  ALL,
  NETWORK_STRAND,                   ///< Connected to two distinct crosslinkers
  TWO_CROSSLINKERS_INCLUDING_LOOPS, ///< Connected to two crosslinkers (incl.
                                    ///< primary loops)
  DANGLING_CHAIN,                   ///< Connected to exactly one crosslinker
  FREE_CHAIN,                       ///< Unconnected (solvent/free chain)
  PRIMARY_LOOP ///< Both ends connected to the same crosslinker
};

/**
 * @brief Summary statistics for realized bond lengths.
 */
struct BondLengthStats
{
  size_t nrOfBonds = 0;
  double meanBondLength = 0.0;
  double meanSquaredBondLength = 0.0;
  double minBondLength = 0.0;
  double maxBondLength = 0.0;
};

/**
 * @brief Summary statistics for realized end-to-end distances.
 */
struct EndToEndStats
{
  size_t nrOfStrands = 0;
  size_t nrOfExcludedZeroLength = 0;
  double meanEndToEndDistance = 0.0;
  double meanSquaredEndToEndDistance = 0.0;
};

/**
 * @brief Classify a strand in a CrosslinkerUniverse.
 */
StrandCategory
classifyStrand(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
               size_t strandIdx);

/**
 * @brief Check whether a strand matches a requested StrandCategory filter.
 */
bool
strandMatchesCategory(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
                      size_t strandIdx,
                      StrandCategory category);

/**
 * @brief Stream over all strands in a CrosslinkerUniverse, realizing bead
 * coordinates one strand at a time in O(1) strand memory, and compute bond
 * length statistics matching Universe::computeBondLengths().
 */
BondLengthStats
computeBondLengthStats(const pylimer_tools::entities::Box& box,
                       const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
                       int nMcSteps = 0);

/**
 * @brief Stream over strands matching @p category in a CrosslinkerUniverse and
 * compute end-to-end distance statistics matching
 * Universe::computeMeanSquareEndToEndDistance() /
 * Molecule::computeEndToEndDistanceWithDerivedImageFlags().
 */
EndToEndStats
computeEndToEndStats(const pylimer_tools::entities::Box& box,
                     const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
                     StrandCategory category = StrandCategory::NETWORK_STRAND,
                     int nMcSteps = 0);

/**
 * @brief Count atoms per atom type in O(N_xlinks + N_strands) without
 * constructing a Universe.
 */
std::map<int, int>
countAtomTypes(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu);

/**
 * @brief Count atoms per bond degree in O(N_xlinks + N_strands) without
 * constructing a Universe.
 */
std::map<int, size_t>
countAtomsByDegree(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu);

/**
 * @brief Compute effective functionality (mean degree) per atom type in
 * O(N_xlinks + N_strands), matching
 * Universe::determineEffectiveFunctionalityPerType().
 */
std::map<int, double>
determineEffectiveFunctionalityPerType(
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu);

/**
 * @brief Return the number of bonds for each strand matching @p category.
 */
std::vector<int>
computeStrandBondCounts(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
                        StrandCategory category = StrandCategory::ALL);

} // namespace pylimer_tools::utils
