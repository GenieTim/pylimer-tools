#include "EntanglementDetector.h"

#include "../entities/EigenNeighbourList.h"
#include "../entities/Molecule.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include "../utils/LazyDistanceMatrix.h"
#include "../utils/StringUtils.h"
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace pylimer_tools::topo::entanglement_detection {
/**
 * @brief Randomly find pairs of atoms that are close together and could be
 * entanglements
 *
 * @param universe the universe with all the atoms etc.
 * @param nrOfEntanglementsToSample the nr of entanglements to find
 * @param cutoff the cut-off distance within which two atoms are considered
 * close enough
 * @param minimumNrOfEntanglements the minimum number of entanglements to
 * find
 * @param sameStrandCutoff the number of beads required between two atoms of
 * the same strand
 * @param seed the random seed
 * @param crossLinkerType the type of the crosslink atoms
 * @param ignoreCrosslinks whether to ignore crosslink atoms when sampling.
 * Careful: if you don't ignore them, the same-strand policy might not work
 * correctly.
 * @param filterDanglingAndSoluble
 * @return AtomPairEntanglements
 */
AtomPairEntanglements
randomlyFindEntanglements(const pylimer_tools::entities::Universe& universe,
                          const size_t nrOfEntanglementsToSample,
                          const double upperCutoff,
                          const double lowerCutoff,
                          const size_t minimumNrOfEntanglements,
                          const double sameStrandCutoff,
                          const std::string& seed,
                          int crossLinkerType,
                          bool ignoreCrosslinks,
                          bool filterDanglingAndSoluble)
{
  INVALIDARG_EXP_IFN(minimumNrOfEntanglements < universe.getNrOfAtoms() / 2,
                     "Minimum number of slip-links must be less than the "
                     "possible number of slip-links to place.");
  INVALIDARG_EXP_IFN(nrOfEntanglementsToSample < universe.getNrOfAtoms() / 2,
                     "Number of slip-links to place must be less than "
                     "the possible number of slip-links to place.");
  INVALIDARG_EXP_IFN(nrOfEntanglementsToSample >= minimumNrOfEntanglements,
                     "Maximum nr. should be larger than minimum, got " +
                       std::to_string(nrOfEntanglementsToSample) + " and " +
                       std::to_string(minimumNrOfEntanglements) + ".");
  INVALIDARG_EXP_IFN(upperCutoff > 0.0,
                     "Expected a cutoff > 0.0, got " +
                       std::to_string(upperCutoff) + ".");
  INVALIDARG_EXP_IFN(
    lowerCutoff < upperCutoff,
    "Expected lower cut-off to be smaller than upper cut-off, got " +
      std::to_string(lowerCutoff) + " ≥ " + std::to_string(upperCutoff) + ".");

  // std::cout << "Randomly finding " << nrOfEntanglementsToSample
  //           << " entanglements within cutoff " << cutoff
  //           << " and same strand cutoff " << sameStrandCutoff << "."
  //           << std::endl;

  // initialise some stuff
  std::vector<std::pair<size_t, size_t>> pairsOfAtoms;
  pairsOfAtoms.reserve(nrOfEntanglementsToSample);
  std::vector<long int> pairOfAtom =
    std::vector<long int>(universe.getNrOfAtoms(), -1);

  if (nrOfEntanglementsToSample == 0) {
    return AtomPairEntanglements(pairsOfAtoms, pairOfAtom);
  }

  // assemble some minor performance benefits
  std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
    universe.getChainsWithCrosslinker(crossLinkerType);
  std::vector<bool> ignoreMolecule =
    std::vector<bool>(crossLinkerChains.size(), filterDanglingAndSoluble);
  if (filterDanglingAndSoluble) {
    std::transform(
      crossLinkerChains.begin(),
      crossLinkerChains.end(),
      ignoreMolecule.begin(),
      [](const pylimer_tools::entities::Molecule& m) {
        return m.getType() ==
                 pylimer_tools::entities::MoleculeType::DANGLING_CHAIN ||
               m.getType() == pylimer_tools::entities::MoleculeType::FREE_CHAIN;
      });
  }

  std::vector<long int> atomToStrand =
    std::vector<long int>(universe.getNrOfAtoms(), -1);
  // use long int here to make unproblematic subtractions
  std::vector<long int> atomIdxInStrand =
    std::vector<long int>(universe.getNrOfAtoms(), -1);
  // start by setting distribution to sample from
  std::vector<pylimer_tools::entities::Atom> atomsForNeighbourList;
  atomsForNeighbourList.reserve(universe.getNrOfAtoms());
  for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
    std::vector<pylimer_tools::entities::Atom> atoms =
      crossLinkerChains[i].getAtomsLinedUp(crossLinkerType, true, true);
    for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
      pylimer_tools::entities::Atom atom = atoms[atomIdx];
      if (((atom.getType() != crossLinkerType) && (atomIdx != 0) &&
           (atomIdx != (atoms.size() - 1))) ||
          !ignoreCrosslinks) {
        igraph_integer_t vertexIdx = universe.getIdxByAtomId(atom.getId());
        atomsForNeighbourList.push_back(atom);
        atomToStrand[vertexIdx] = i;
        atomIdxInStrand[vertexIdx] = atomIdx;
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
      atomsForNeighbourList, universe.getBox(), upperCutoff);
  size_t numLinksFoundInIteration = 1;

  // start iteration to sample "entanglements"
  size_t numEntanglementsSampled = 0;
  do {
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
        neighbourList.getAtomsCloseTo(a1, upperCutoff, lowerCutoff);
      neighbourList.removeAtom(a1, "After querying neighbours.");
      // filter the neighbours to include only those from other strands
      // NOTE: this skews the whole thing a bit
      std::erase_if(
        neighbours, [&](const pylimer_tools::entities::Atom& a) -> bool {
          return (((atomToStrand[universe.getIdxByAtomId(a.getId())] ==
                    atomToStrand[atomVertexIdx1]) &&
                   (std::abs(static_cast<double>(
                      atomIdxInStrand[universe.getIdxByAtomId(a.getId())] -
                      atomIdxInStrand[atomVertexIdx1])) < sameStrandCutoff)));
        });
      if (neighbours.size() == 0) {
        // std::cerr << "Not enough close neighbours found." << std::endl;
        continue;
      }
      // then, randomly select one of them
      pylimer_tools::entities::Atom a2 = neighbours[0];
      if (neighbours.size() > 1) {
        size_t randomA2Idx =
          std::uniform_int_distribution<size_t>{ 0,
                                                 neighbours.size() - 1 }(rng);
        a2 = neighbours[randomA2Idx];
      }

      size_t atomVertexIdx2 = universe.getIdxByAtomId(a2.getId());
      RUNTIME_EXP_IFN(
        pairOfAtom[atomVertexIdx2] == -1,
        "Expected not to be able to sample the same atom "
        "twice. Sampled atom with vertex " +
          std::to_string(atomVertexIdx2) + " and id " +
          std::to_string(a2.getId()) + " again, already present in " +
          std::to_string(pairOfAtom[atomVertexIdx2]) + " being " +
          std::to_string(pairsOfAtoms[pairOfAtom[atomVertexIdx2]]) +
          ". This time sampled with atom " + std::to_string(a1.getId()) + ".");
      RUNTIME_EXP_IFN(
        (atomToStrand[atomVertexIdx2] != atomToStrand[atomVertexIdx1]) ||
          ((std::abs(static_cast<double>(atomIdxInStrand[atomVertexIdx2] -
                                         atomIdxInStrand[atomVertexIdx1])) >=
            sameStrandCutoff)),
        "Expected neighbours to have been deleted if too far apart.");

      if (!filterDanglingAndSoluble ||
          (!ignoreMolecule[atomToStrand[atomVertexIdx2]] &&
           !ignoreMolecule[atomToStrand[atomVertexIdx1]])) {
        pairOfAtom[atomVertexIdx2] = pairsOfAtoms.size();
        pairOfAtom[atomVertexIdx1] = pairsOfAtoms.size();
        pairsOfAtoms.emplace_back(a1.getId(), a2.getId());
      } else {
        pairOfAtom[atomVertexIdx2] = -2;
        pairOfAtom[atomVertexIdx1] = -2;
      }
      numEntanglementsSampled += 1;
      numLinksFoundInIteration += 1;
      neighbourList.removeAtom(a2, "After marking atom as second pair part.");
      // std::cout << "Merging atoms " << a1.getId() << " ("
      //           << atomIdxInStrand[a1.getId()] << "th in "
      //           << atomToStrand[atomVertexIdx1] << ") and " << a2.getId()
      //           << " (" << atomIdxInStrand[a2.getId()] << "th in "
      //           << atomToStrand[atomVertexIdx2] << ") "
      //           << " with distance " << a1.distanceTo(a2,
      //           universe.getBox())
      //           << std::endl;
      // if (atomToStrand[atomVertexIdx2] == atomToStrand[atomVertexIdx1]) {
      //   std::cout << "Same strand detected: distance is "
      //             << std::abs(
      //                  static_cast<double>(atomIdxInStrand[a2.getId()] -
      //                                      atomIdxInStrand[a1.getId()]))
      //             << std::endl;
      // }
      if (numEntanglementsSampled >= nrOfEntanglementsToSample) {
        break;
      }
    }
    if (numEntanglementsSampled >= nrOfEntanglementsToSample) {
      break;
    }
  } while (numEntanglementsSampled < minimumNrOfEntanglements &&
           numLinksFoundInIteration > 0);

  // reset internal thing
  for (long& i : pairOfAtom) {
    if (i == -2) {
      i = -1;
    }
  }

  AtomPairEntanglements result;
  result.pairsOfAtoms = pairsOfAtoms;
  result.pairOfAtom = pairOfAtom;
  return result;
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
 * @param ignoreCrosslinks whether to ignore crosslink atoms when sampling.
 * Careful: if you don't ignore them, the same-strand policy might not work
 * correctly.
 * @param filterDanglingAndSoluble
 * @return AtomPairEntanglements
 */
AtomPairEntanglements
randomlyFindEntanglementsV2(const pylimer_tools::entities::Universe& universe,
                            const size_t nrOfEntanglementsToSample,
                            const double upperCutoff,
                            const double lowerCutoff,
                            const size_t minimumNrOfEntanglements,
                            const double sameStrandCutoff,
                            const std::string& seed,
                            const int crossLinkerType,
                            const bool ignoreCrosslinks,
                            bool filterDanglingAndSoluble)
{
  INVALIDARG_EXP_IFN(minimumNrOfEntanglements < universe.getNrOfAtoms() / 2,
                     "Minimum number of slip-links must be less than the "
                     "possible number of slip-links to place.");
  INVALIDARG_EXP_IFN(nrOfEntanglementsToSample < universe.getNrOfAtoms() / 2,
                     "Number of slip-links to place must be less than "
                     "the possible number of slip-links to place.");
  INVALIDARG_EXP_IFN(nrOfEntanglementsToSample >= minimumNrOfEntanglements,
                     "Maximum nr. should be larger than minimum, got " +
                       std::to_string(nrOfEntanglementsToSample) + " and " +
                       std::to_string(minimumNrOfEntanglements) + ".");
  INVALIDARG_EXP_IFN(upperCutoff > 0.0,
                     "Expected a cutoff > 0.0, got " +
                       std::to_string(upperCutoff) + ".");
  INVALIDARG_EXP_IFN(
    lowerCutoff < upperCutoff,
    "Expected lower cut-off to be smaller than upper cut-off, got " +
      std::to_string(lowerCutoff) + " ≥ " + std::to_string(upperCutoff) + ".");

  // std::cout << "Randomly finding " << nrOfEntanglementsToSample
  //           << " entanglements within cutoff " << cutoff
  //           << " and same strand cutoff " << sameStrandCutoff << "."
  //           << std::endl;

  // initialise some stuff
  std::vector<std::pair<size_t, size_t>> pairsOfAtoms;
  pairsOfAtoms.reserve(nrOfEntanglementsToSample);
  std::vector<long int> pairOfAtom =
    std::vector<long int>(universe.getNrOfAtoms(), -1);

  if (nrOfEntanglementsToSample == 0) {
    return AtomPairEntanglements(pairsOfAtoms, pairOfAtom);
  }

  igraph_lazy_distance_matrix_state_t distance_matrix_computer;
  igraph_t graph = universe.getCopyOfGraph();
  igraph_lazy_distance_matrix_state_init(
    &distance_matrix_computer, &graph, IGRAPH_ALL, sameStrandCutoff);

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
  pylimer_tools::entities::Box box = universe.getBox();
  Eigen::VectorXd coordinates = universe.getUnwrappedVertexCoordinates(box);
  const pylimer_tools::entities::EigenNeighbourList neighbourList =
    pylimer_tools::entities::EigenNeighbourList(coordinates, box, upperCutoff);
  std::vector<bool> vertexIsEligible(universe.getNrOfAtoms(), true);
  const std::vector<int> vertexTypes = universe.getPropertyValues<int>("type");
  const std::vector<int> vertexDegrees = universe.getVertexDegrees();

  std::vector<long int> vertexIdsToSample;
  vertexIdsToSample.reserve(universe.getNrOfAtoms());

  const std::vector<pylimer_tools::entities::MoleculeType>
    moleculeTypePerVertex =
      filterDanglingAndSoluble
        ? universe.identifyObviouslyDanglingAtoms()
        : std::vector({ pylimer_tools::entities::MoleculeType::UNDEFINED });

  // filter out crosslink atoms and soluble atoms

  for (size_t i = 0; i < universe.getNrOfAtoms(); ++i) {
    if (ignoreCrosslinks && vertexTypes[i] == crossLinkerType) {
      vertexIsEligible[i] = false;
    }
    if (ignoreCrosslinks && vertexDegrees[i] > 2) {
      vertexIsEligible[i] = false;
    }
    if (vertexIsEligible[i]) {
      vertexIdsToSample.push_back(i);
    }
  }

  size_t numLinksFoundInIteration = 1;
  Eigen::ArrayXi neighbourIndices = Eigen::ArrayXi::Zero(12);

  // start iteration to sample "entanglements"
  size_t numEntanglementsSampled = 0;
  do {
    numLinksFoundInIteration = 0;
    std::ranges::shuffle(vertexIdsToSample, rng);
    for (const long int atomVertexIdx1 : vertexIdsToSample) {
      // make sure this atom does not yet have a pair
      if (pairOfAtom[atomVertexIdx1] != -1) {
        continue;
      }
      vertexIsEligible[atomVertexIdx1] = false;
      // then, find neighbouring atoms
      const size_t nNeighbours = neighbourList.getIndicesCloseToCoordinates(
        neighbourIndices,
        coordinates.segment(3 * atomVertexIdx1, 3),
        upperCutoff,
        true);
      long int atomVertexIdx2 = -1;

      std::shuffle(
        neighbourIndices.data(), neighbourIndices.data() + nNeighbours, rng);

      for (Eigen::Index idxInNeighbours = 0; idxInNeighbours < nNeighbours;
           ++idxInNeighbours) {
        // filter these neighbours
        if (!vertexIsEligible[neighbourIndices[idxInNeighbours]]) {
          continue;
        }

        Eigen::Vector3d distance =
          coordinates.segment(3 * neighbourIndices[idxInNeighbours], 3) -
          coordinates.segment(3 * atomVertexIdx1, 3);
        box.handlePBC(distance);

        if (distance.norm() < lowerCutoff || distance.norm() > upperCutoff) {
          continue;
        }

        igraph_integer_t pathLength =
          sameStrandCutoff <= 0 ? IGRAPH_INFINITY
                                : igraph_lazy_distance_matrix_path_length(
                                    &distance_matrix_computer,
                                    atomVertexIdx1,
                                    neighbourIndices[idxInNeighbours],
                                    false);
        if (sameStrandCutoff <= 0 || pathLength < sameStrandCutoff) {
          // found the second part of the pair
          atomVertexIdx2 = neighbourIndices[idxInNeighbours];
          break;
        }
      }
      igraph_lazy_distance_matrix_forget_source(&distance_matrix_computer,
                                                atomVertexIdx1);

      if (atomVertexIdx2 < 0) {
        // did not find eligible neighbour, continue to next starting atom
        continue;
      }

      RUNTIME_EXP_IFN(
        pairOfAtom[atomVertexIdx2] == -1,
        "Expected not to be able to sample the same atom "
        "twice. Sampled atom with vertex " +
          std::to_string(atomVertexIdx2) + " again, already present in " +
          std::to_string(pairOfAtom[atomVertexIdx2]) + " being " +
          std::to_string(pairsOfAtoms[pairOfAtom[atomVertexIdx2]]) +
          ". This time sampled with vertex " + std::to_string(atomVertexIdx1) +
          ".");

      if (!filterDanglingAndSoluble ||
          !((moleculeTypePerVertex[atomVertexIdx1] ==
               pylimer_tools::entities::MoleculeType::FREE_CHAIN ||
             moleculeTypePerVertex[atomVertexIdx1] ==
               pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) ||
            (moleculeTypePerVertex[atomVertexIdx2] ==
               pylimer_tools::entities::MoleculeType::FREE_CHAIN ||
             moleculeTypePerVertex[atomVertexIdx2] ==
               pylimer_tools::entities::MoleculeType::DANGLING_CHAIN))) {
        pairOfAtom[atomVertexIdx2] = pairsOfAtoms.size();
        pairOfAtom[atomVertexIdx1] = pairsOfAtoms.size();
        pairsOfAtoms.emplace_back(universe.getAtomIdByIdx(atomVertexIdx1),
                                  universe.getAtomIdByIdx(atomVertexIdx2));
      } else {
        // dangling atoms, entangled but we don't store it
        pairOfAtom[atomVertexIdx2] = -2;
        pairOfAtom[atomVertexIdx1] = -2;
      }
      vertexIsEligible[atomVertexIdx2] = false;

      numEntanglementsSampled += 1;
      numLinksFoundInIteration += 1;

      if (numEntanglementsSampled >= nrOfEntanglementsToSample) {
        break;
      }
    }
    if (numEntanglementsSampled >= nrOfEntanglementsToSample) {
      break;
    }
  } while (numEntanglementsSampled < minimumNrOfEntanglements &&
           numLinksFoundInIteration > 0);

  // reset internal thing
  for (long& i : pairOfAtom) {
    if (i == -2) {
      i = -1;
    }
  }
  igraph_lazy_distance_matrix_state_destroy(&distance_matrix_computer);
  igraph_destroy(&graph);

  AtomPairEntanglements result;
  result.pairsOfAtoms = pairsOfAtoms;
  result.pairOfAtom = pairOfAtom;
  return result;
}
}
