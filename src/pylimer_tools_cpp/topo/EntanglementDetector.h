#ifndef ENTANGLEMENT_DETECTOR_H
#define ENTANGLEMENT_DETECTOR_H

#include "../entities/Universe.h"
#include <string>
#include <vector>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::topo::entanglement_detection {

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
 * @param filterDanglingAndSoluble
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
                          bool filterDanglingAndSoluble = false);

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
 * @param filterDanglingAndSoluble
 * @return AtomPairEntanglements
 */
AtomPairEntanglements
randomlyFindEntanglementsV2(const pylimer_tools::entities::Universe& universe,
                            const size_t nrOfEntanglementsToSample,
                            const double upperCutoff,
                            const double lowerCutoff = 0.,
                            const size_t minimumNrOfEntanglements = 0,
                            const double sameStrandCutoff = 3,
                            const std::string& seed = "",
                            int crossLinkerType = 2,
                            bool ignoreCrosslinks = true,
                            bool filterDanglingAndSoluble = true);

}

#endif
