#include "EntanglementDetector.h"
#include "../entities/Molecule.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include <iostream>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace calc {

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
      int crossLinkerType)
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
      std::unordered_map<size_t, size_t> atomIdxInStrand;
      atomIdxInStrand.reserve(universe.getNrOfAtoms());
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        pylimer_tools::entities::Molecule chain = crossLinkerChains[i];
        RUNTIME_EXP_IFN(chain.getType() !=
                          pylimer_tools::entities::MoleculeType::UNDEFINED,
                        "Couldn't determine molecule type.");
        std::vector<pylimer_tools::entities::Atom> atoms =
          crossLinkerChains[i].getAtomsLinedUp(crossLinkerType, true, true);
        for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
          pylimer_tools::entities::Atom atom = atoms[atomIdx];
          if (atom.getType() != crossLinkerType) {
            atomToStrand.emplace(atom.getId(), i);
            atomIdxInStrand.emplace(atom.getId(), atomIdx);
          }
        }
      }

      // start by setting distribution to sample from
      // filter, we don't want cross-links etc. as targets
      std::vector<pylimer_tools::entities::Atom> atomsForNeighbourList =
        universe.getAtomsOfDegree(2);
      atomsForNeighbourList.erase(
        std::remove_if(
          atomsForNeighbourList.begin(),
          atomsForNeighbourList.end(),
          [crossLinkerType](const pylimer_tools::entities::Atom& a) {
            return a.getType() == crossLinkerType;
          }),
        atomsForNeighbourList.end());
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
          neighbourList.removeAtom(
            a1, "After querying neighbours. Impossible case.");
          // filter the neighbours to include only those from other strands
          // NOTE: this skews the whole thing a bit
          neighbours.erase(
            std::remove_if(
              neighbours.begin(),
              neighbours.end(),
              [&](const pylimer_tools::entities::Atom& a) -> bool {
                return (
                  (atomToStrand[a.getId()] ==
                     atomToStrand[a1.getId()] // do not use "at", because not
                                              // all atoms in the neighbours
                                              // have been assigned a strand
                   && (std::abs(static_cast<double>(
                         atomIdxInStrand[a.getId()] -
                         atomIdxInStrand[a1.getId()])) < sameStrandCutoff)));
              }),
            neighbours.end());
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
          assert(pairOfAtom[atomVertexIdx2] == -1);
          pairOfAtom[atomVertexIdx2] = pairsOfAtoms.size();
          pairOfAtom[atomVertexIdx1] = pairsOfAtoms.size();
          pairsOfAtoms.push_back(std::make_pair(a1.getId(), a2.getId()));
          numLinksFoundInIteration += 1;
          neighbourList.removeAtom(a2,
                                   "After marking atom as second pair part.");
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
