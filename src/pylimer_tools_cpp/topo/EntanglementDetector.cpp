#include "EntanglementDetector.h"
#include "../entities/Molecule.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include <iostream>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace topo {

  namespace entanglement_detection {

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
      const double cutoff,
      const size_t minimumNrOfSliplinks,
      const double sameStrandCutoff,
      const std::string& seed,
      int crossLinkerType,
      bool ignoreCrosslinks)
    {
      INVALIDARG_EXP_IFN(minimumNrOfSliplinks < universe.getNrOfAtoms() / 2,
                         "Minimum number of slip-links must be less than the "
                         "possible number of slip-links to place.");
      INVALIDARG_EXP_IFN(nrOfSliplinksToSample < universe.getNrOfAtoms() / 2,
                         "Number of slip-links to place must be less than "
                         "the possible number of slip-links to place.");
      INVALIDARG_EXP_IFN(nrOfSliplinksToSample >= minimumNrOfSliplinks,
                         "Maximum nr. should be larger than minimum, got " +
                           std::to_string(nrOfSliplinksToSample) + " and " +
                           std::to_string(minimumNrOfSliplinks) + ".");
      INVALIDARG_EXP_IFN(cutoff > 0.0,
                         "Expected a cutoff > 0.0, got " +
                           std::to_string(cutoff) + ".");

      // std::cout << "Randomly finding " << nrOfSliplinksToSample
      //           << " entanglements within cutoff " << cutoff
      //           << " and same strand cutoff " << sameStrandCutoff << "."
      //           << std::endl;

      // initialise some stuff
      std::vector<std::pair<size_t, size_t>> pairsOfAtoms;
      pairsOfAtoms.reserve(nrOfSliplinksToSample);
      std::vector<long int> pairOfAtom =
        pylimer_tools::utils::initializeWithValue<long int>(
          universe.getNrOfAtoms(), -1);

      // assemble some minor performance benefits
      std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
        universe.getChainsWithCrosslinker(crossLinkerType);

      std::unordered_map<size_t, size_t> atomToStrand;
      atomToStrand.reserve(universe.getNrOfAtoms());
      // use long int here to make unproblematic subtractions
      std::unordered_map<size_t, long int> atomIdxInStrand;
      atomIdxInStrand.reserve(universe.getNrOfAtoms());
      // start by setting distribution to sample from
      std::vector<pylimer_tools::entities::Atom> atomsForNeighbourList;
      atomsForNeighbourList.reserve(universe.getNrOfAtoms());
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        std::vector<pylimer_tools::entities::Atom> atoms =
          crossLinkerChains[i].getAtomsLinedUp(crossLinkerType, true, true);
        for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
          pylimer_tools::entities::Atom atom = atoms[atomIdx];
          if (atom.getType() != crossLinkerType || !ignoreCrosslinks) {
            atomsForNeighbourList.push_back(atom);
            atomToStrand.emplace(atom.getId(), i);
            atomIdxInStrand.emplace(atom.getId(), atomIdx);
          }
        }
      }

      // some randomness for placement
      std::mt19937 rng;
      if (seed == "") {
        std::random_device rd{};
        rng = std::mt19937(rd());
      } else {
        std::seed_seq seed2(seed.begin(), seed.end());
        rng = std::mt19937(seed2);
      }
      // std::cout << "Initial sampling rng seed: " << rng << std::endl;
      pylimer_tools::entities::NeighbourList neighbourList =
        pylimer_tools::entities::NeighbourList(
          atomsForNeighbourList, universe.getBox(), cutoff);
      size_t numLinksFoundInIteration = 1;

      // start iteration to sample "entanglements"
      while (pairsOfAtoms.size() < minimumNrOfSliplinks &&
             numLinksFoundInIteration > 0) {
        numLinksFoundInIteration = 0;
        std::shuffle(
          atomsForNeighbourList.begin(), atomsForNeighbourList.end(), rng);
        for (pylimer_tools::entities::Atom a1 : atomsForNeighbourList) {
          size_t atomVertexIdx1 = universe.getIdxByAtomId(a1.getId());
          // make sure this atom does not yet have a pair
          if (pairOfAtom[atomVertexIdx1] != -1) {
            continue;
          }
          // then, find neighbouring atoms (but not from the same strand?!)
          std::vector<pylimer_tools::entities::Atom> neighbours =
            neighbourList.getAtomsCloseTo(a1);
          neighbourList.removeAtom(a1, "After querying neighbours.");
          // filter the neighbours to include only those from other strands
          // NOTE: this skews the whole thing a bit
          std::erase_if(
            neighbours, [&](const pylimer_tools::entities::Atom& a) -> bool {
              return ((
                (atomToStrand.at(a.getId()) == atomToStrand.at(a1.getId())) &&
                (std::abs(static_cast<double>(atomIdxInStrand.at(a.getId()) -
                                              atomIdxInStrand.at(a1.getId()))) <
                 sameStrandCutoff)));
            });
          if (neighbours.size() == 0) {
            // std::cerr << "Not enough close neighbours found." << std::endl;
            continue;
          }
          // then, randomly select one of them
          pylimer_tools::entities::Atom a2 = neighbours[0];
          if (neighbours.size() > 1) {
            size_t randomA2Idx =
              std::uniform_int_distribution<size_t>{ 0, neighbours.size() - 1 }(
                rng);
            a2 = neighbours[randomA2Idx];
          }

          size_t atomVertexIdx2 = universe.getIdxByAtomId(a2.getId());
          RUNTIME_EXP_IFN(
            pairOfAtom[atomVertexIdx2] == -1,
            "Expected not to be able to sample the same atom twice.");
          RUNTIME_EXP_IFN(
            (atomToStrand[a2.getId()] != atomToStrand[a1.getId()]) ||
              ((std::abs(static_cast<double>(atomIdxInStrand[a2.getId()] -
                                             atomIdxInStrand[a1.getId()])) >=
                sameStrandCutoff)),
            "Expected neighbours to have been deleted if too far apart.");
          pairOfAtom[atomVertexIdx2] = pairsOfAtoms.size();
          pairOfAtom[atomVertexIdx1] = pairsOfAtoms.size();
          pairsOfAtoms.push_back(std::make_pair(a1.getId(), a2.getId()));
          numLinksFoundInIteration += 1;
          neighbourList.removeAtom(a2,
                                   "After marking atom as second pair part.");
          // std::cout << "Merging atoms " << a1.getId() << " ("
          //           << atomIdxInStrand[a1.getId()] << "th in "
          //           << atomToStrand[a1.getId()] << ") and " << a2.getId()
          //           << " (" << atomIdxInStrand[a2.getId()] << "th in "
          //           << atomToStrand[a2.getId()] << ") "
          //           << " with distance " << a1.distanceTo(a2, universe.getBox())
          //           << std::endl;
          // if (atomToStrand[a2.getId()] == atomToStrand[a1.getId()]) {
          //   std::cout << "Same strand detected: distance is "
          //             << std::abs(
          //                  static_cast<double>(atomIdxInStrand[a2.getId()] -
          //                                      atomIdxInStrand[a1.getId()]))
          //             << std::endl;
          // }
          if (pairsOfAtoms.size() >= nrOfSliplinksToSample) {
            break;
          }
        }
        if (pairsOfAtoms.size() >= nrOfSliplinksToSample) {
          break;
        }
      }

      AtomPairEntanglements result;
      result.pairsOfAtoms = pairsOfAtoms;
      result.pairOfAtom = pairOfAtom;
      return result;
    }

  }
}
}
