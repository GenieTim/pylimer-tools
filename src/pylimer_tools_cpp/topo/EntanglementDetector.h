#ifndef ENTANGLEMENT_DETECTOR_H
#define ENTANGLEMENT_DETECTOR_H

#include "../entities/Universe.h"
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools {
namespace topo {

  namespace entanglement_detection {

    struct AtomPairEntanglements
    {
      /**
       * pairs of atom ids
       */
      std::vector<std::pair<size_t, size_t>> pairsOfAtoms;
      /**
       * For each atom, the index of the pair it is associated with, or -1 if
       * none
       */
      std::vector<long int> pairOfAtom;
    };

    /**
     * @brief Randomly find pairs of atoms that are close together and could be
     * entanglements
     *
     * @param universe the universe with all the atoms etc.
     * @param nrOfSliplinksToSample the nr of entanglements to find
     * @param cutoff the cut-off distance within which two atoms are considered
     * close enough
     * @param minimumNrOfSliplinks the minimum number of entanglements to find
     * @param sameStrandCutoff the number of beads required between two atoms of
     * the same strand
     * @param seed the random seed
     * @param crossLinkerType the type of the cross-link atoms
     * @return AtomPairEntanglements
     */
    AtomPairEntanglements randomlyFindEntanglements(
      const pylimer_tools::entities::Universe& universe,
      const size_t nrOfSliplinksToSample,
      const double upperCutoff,
      const double lowerCutoff = 0.,
      const size_t minimumNrOfSliplinks = 0,
      const double sameStrandCutoff = 3,
      const std::string& seed = "",
      int crossLinkerType = 2,
      bool ignoreCrosslinks = true,
      bool filterDanglingAndSoluble = false);

  }

}
}

#endif
