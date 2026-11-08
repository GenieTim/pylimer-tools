#include "EntanglementDetector.h"

#include "../entities/EigenNeighbourList.h"
#include "../entities/Molecule.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include "../utils/LazyDistanceMatrix.h"
#include "../utils/StringUtils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
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
    for (const pylimer_tools::entities::Atom& a1 : atomsForNeighbourList) {
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

        const double dist2 = distance.squaredNorm();
        if (dist2 < lowerCutoff * lowerCutoff ||
            dist2 > upperCutoff * upperCutoff) {
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

AtomPairEntanglements
StrandPairEntanglements::toAtomPairEntanglements(
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse)
  const
{
  const size_t nCrosslinks = crosslinkerUniverse.xlinkX.size();
  const size_t nStrands = crosslinkerUniverse.beadsInStrand.size();
  std::vector<size_t> firstBeadVertexIdx(nStrands, nCrosslinks);
  size_t currentVertexIdx = nCrosslinks;
  for (size_t s = 0; s < nStrands; ++s) {
    firstBeadVertexIdx[s] = currentVertexIdx;
    currentVertexIdx += static_cast<size_t>(
      std::max<int>(0, crosslinkerUniverse.beadsInStrand[s]));
  }
  const size_t totalAtoms = currentVertexIdx;

  AtomPairEntanglements result;
  result.pairsOfAtoms.reserve(this->size());
  result.pairOfAtom.assign(totalAtoms, -1);

  for (size_t e = 0; e < this->size(); ++e) {
    const size_t s1 = this->strandIdx1[e];
    const size_t s2 = this->strandIdx2[e];
    INVALIDARG_EXP_IFN(
      s1 < nStrands && s2 < nStrands,
      "Strand index out of bounds in toAtomPairEntanglements.");
    const size_t k1 = this->beadIdx1[e];
    const size_t k2 = this->beadIdx2[e];
    INVALIDARG_EXP_IFN(
      static_cast<int>(k1) < crosslinkerUniverse.beadsInStrand[s1] &&
        static_cast<int>(k2) < crosslinkerUniverse.beadsInStrand[s2],
      "Bead index out of bounds in toAtomPairEntanglements.");

    const size_t v1 = firstBeadVertexIdx[s1] + k1;
    const size_t v2 = firstBeadVertexIdx[s2] + k2;
    const size_t pairIdx = result.pairsOfAtoms.size();
    result.pairOfAtom[v1] = static_cast<long int>(pairIdx);
    result.pairOfAtom[v2] = static_cast<long int>(pairIdx);
    // Atom IDs in MCUniverseGenerator::getUniverse() are 1-based (vertexIdx +
    // 1)
    result.pairsOfAtoms.emplace_back(v1 + 1, v2 + 1);
  }

  return result;
}

namespace {

  struct PinnedPoint
  {
    double alpha;
    Eigen::Vector3d pos;
  };

  struct StrandModel
  {
    Eigen::Index from = -1;
    Eigen::Index to = -1;
    Eigen::Index nBeads = 0;
    Eigen::Index nBonds = 0;
    Eigen::Index nInternalBeads = 0;
    double sigma2Bond = 1.0 / 3.0;
    double maxSigma = 0.0;
    bool isIgnored = false;
    bool reversedAnchor = false;
    Eigen::Vector3d anchorStart = Eigen::Vector3d::Zero();
    Eigen::Vector3d anchorEnd = Eigen::Vector3d::Zero();
    Eigen::Vector3d midpoint = Eigen::Vector3d::Zero();
    double halfChordLen = 0.0;
    std::vector<PinnedPoint> pinned;

    void evaluateConditionalDistribution(const double alpha,
                                         Eigen::Vector3d& meanOut,
                                         double& sigmaOut) const
    {
      assert(!this->pinned.empty());
      auto it = std::lower_bound(
        this->pinned.begin(),
        this->pinned.end(),
        alpha,
        [](const PinnedPoint& p, const double a) { return p.alpha < a; });

      if (it == this->pinned.begin()) {
        const PinnedPoint& right = *it;
        const double dAlpha = std::max(0.0, right.alpha - alpha);
        meanOut = right.pos;
        sigmaOut = std::sqrt(dAlpha * static_cast<double>(this->nBonds) *
                             this->sigma2Bond);
        return;
      }

      if (it == this->pinned.end()) {
        const PinnedPoint& left = this->pinned.back();
        const double dAlpha = std::max(0.0, alpha - left.alpha);
        meanOut = left.pos;
        sigmaOut = std::sqrt(dAlpha * static_cast<double>(this->nBonds) *
                             this->sigma2Bond);
        return;
      }

      const PinnedPoint& right = *it;
      const PinnedPoint& left = *(it - 1);
      const double span = right.alpha - left.alpha;
      if (span <= 1e-14) {
        meanOut = left.pos;
        sigmaOut = 0.0;
        return;
      }
      const double t = std::clamp((alpha - left.alpha) / span, 0.0, 1.0);
      meanOut = left.pos + t * (right.pos - left.pos);
      const double var = t * (1.0 - t) * span *
                         static_cast<double>(this->nBonds) * this->sigma2Bond;
      sigmaOut = std::sqrt(std::max(0.0, var));
    }

    void insertPinnedPoint(const double alpha, const Eigen::Vector3d& pos)
    {
      auto it = std::lower_bound(
        this->pinned.begin(),
        this->pinned.end(),
        alpha,
        [](const PinnedPoint& p, const double a) { return p.alpha < a; });
      this->pinned.insert(it, PinnedPoint{ alpha, pos });
    }

    void removePinnedPoint(const double alpha)
    {
      for (auto it = this->pinned.begin(); it != this->pinned.end(); ++it) {
        if (std::abs(it->alpha - alpha) <= 1e-14) {
          this->pinned.erase(it);
          return;
        }
      }
    }

    size_t internalIdxToBeadIdx(const size_t m) const
    {
      if (this->from >= 0 && this->to >= 0) {
        return m;
      }
      if (this->from >= 0 || this->to >= 0) {
        if (!this->reversedAnchor) {
          return m;
        }
        return static_cast<size_t>(this->nBeads - 1) - m;
      }
      return m + 1;
    }

    double internalIdxToAlpha(const size_t m) const
    {
      return static_cast<double>(m + 1) / static_cast<double>(this->nBonds);
    }
  };

  double isotropic3DGaussianRadialCDF(const double r,
                                      const double mu,
                                      const double sigma)
  {
    if (r <= 0.0) {
      return 0.0;
    }
    if (sigma <= 1e-12) {
      return (mu <= r) ? 1.0 : 0.0;
    }

    static const double sqrt2 = std::sqrt(2.0);
    static const double twoOverSqrtPi = 2.0 / std::sqrt(M_PI);
    static const double invSqrt2Pi = 1.0 / std::sqrt(2.0 * M_PI);

    const double invSqrt2Sigma = 1.0 / (sqrt2 * sigma);
    if (mu <= 1e-7 * sigma) {
      const double u = r * invSqrt2Sigma;
      const double eps = mu * invSqrt2Sigma;
      const double expMinusU2 = std::exp(-u * u);
      const double f0 = std::erf(u) - twoOverSqrtPi * u * expMinusU2;
      const double corr =
        (2.0 / 3.0) * twoOverSqrtPi * (u * u * u) * expMinusU2 * (eps * eps);
      return std::clamp(f0 - corr, 0.0, 1.0);
    }

    const double a = (r - mu) * invSqrt2Sigma;
    if (a <= -6.0) {
      return 0.0;
    }
    if (a >= 6.0) {
      return 1.0;
    }

    const double b = (r + mu) * invSqrt2Sigma;
    const double erfcB = (b < 6.0) ? std::erfc(b) : 0.0;
    double erfDiff = 0.0;
    if (a < 0.0) {
      erfDiff = 0.5 * (std::erfc(-a) - erfcB);
    } else {
      erfDiff = 1.0 - 0.5 * (std::erfc(a) + erfcB);
    }

    const double b2MinusA2 = (2.0 * r * mu) / (sigma * sigma);
    const double expDiff = (b2MinusA2 < 36.0)
                             ? -std::exp(-a * a) * std::expm1(-b2MinusA2)
                             : std::exp(-a * a);
    const double expTerm = ((sigma * invSqrt2Pi) / mu) * expDiff;
    return std::clamp(erfDiff - expTerm, 0.0, 1.0);
  }

  double shellProbability3D(const double rMin,
                            const double rMax,
                            const double mu,
                            const double sigma)
  {
    if (sigma <= 1e-12) {
      return (mu >= rMin && mu <= rMax) ? 1.0 : 0.0;
    }
    // Quick 6-sigma range rejection (Gaussian tail beyond 6 sigma has P < 1e-9)
    if (mu - 6.0 * sigma > rMax || mu + 6.0 * sigma < rMin) {
      return 0.0;
    }
    const double fMax = isotropic3DGaussianRadialCDF(rMax, mu, sigma);
    if (fMax <= 0.0) {
      return 0.0;
    }
    const double fMin =
      (rMin > 0.0) ? isotropic3DGaussianRadialCDF(rMin, mu, sigma) : 0.0;
    return std::clamp(fMax - fMin, 0.0, 1.0);
  }

  Eigen::Vector3d sampleDisplacementInShell3D(const double rMin,
                                              const double rMax,
                                              const Eigen::Vector3d& deltaMu,
                                              const double sigma,
                                              std::mt19937& rng)
  {
    const double mu = deltaMu.norm();
    if (sigma <= 1e-12) {
      return deltaMu;
    }

    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    const double fMin =
      (rMin > 0.0) ? isotropic3DGaussianRadialCDF(rMin, mu, sigma) : 0.0;
    const double fMax = isotropic3DGaussianRadialCDF(rMax, mu, sigma);

    double rSample = 0.5 * (rMin + rMax);
    if (fMax - fMin <= 1e-14) {
      const double u = uniform01(rng);
      const double rMin3 = rMin * rMin * rMin;
      const double rMax3 = rMax * rMax * rMax;
      rSample = std::cbrt(rMin3 + u * (rMax3 - rMin3));
    } else {
      const double targetCDF = fMin + uniform01(rng) * (fMax - fMin);
      double lo = rMin;
      double hi = rMax;
      for (int iter = 0; iter < 32; ++iter) {
        const double mid = 0.5 * (lo + hi);
        if (isotropic3DGaussianRadialCDF(mid, mu, sigma) < targetCDF) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      rSample = 0.5 * (lo + hi);
    }

    const double u1 = std::clamp(uniform01(rng), 1e-15, 1.0);
    const double u2 = uniform01(rng);
    const double phi = 2.0 * M_PI * u2;

    const double kappa = (mu * rSample) / (sigma * sigma);
    double cosTheta = 1.0 - 2.0 * u1;
    if (kappa > 1e-6 && mu > 1e-12) {
      cosTheta =
        1.0 + std::log(u1 + (1.0 - u1) * std::exp(-2.0 * kappa)) / kappa;
      cosTheta = std::clamp(cosTheta, -1.0, 1.0);
    }
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));

    Eigen::Vector3d e3 =
      (mu > 1e-12) ? (deltaMu / mu) : Eigen::Vector3d(0.0, 0.0, 1.0);
    Eigen::Vector3d helper = (std::abs(e3.x()) < 0.9)
                               ? Eigen::Vector3d(1.0, 0.0, 0.0)
                               : Eigen::Vector3d(0.0, 1.0, 0.0);
    Eigen::Vector3d e1 = e3.cross(helper).normalized();
    Eigen::Vector3d e2 = e3.cross(e1).normalized();

    return rSample * (cosTheta * e3 +
                      sinTheta * (std::cos(phi) * e1 + std::sin(phi) * e2));
  }

} // anonymous namespace

StrandPairEntanglements
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
  const size_t nCrosslinks = crosslinkerUniverse.xlinkX.size();
  const size_t nStrands = crosslinkerUniverse.beadsInStrand.size();
  INVALIDARG_EXP_IFN(crosslinkerUniverse.xlinkY.size() == nCrosslinks &&
                       crosslinkerUniverse.xlinkZ.size() == nCrosslinks,
                     "Crosslinker coordinate vectors xlinkX, xlinkY, xlinkZ "
                     "must have matching sizes.");
  INVALIDARG_EXP_IFN(crosslinkerUniverse.strandFrom.size() == nStrands &&
                       crosslinkerUniverse.strandTo.size() == nStrands,
                     "Strand connectivity vectors strandFrom, strandTo, "
                     "beadsInStrand must have matching sizes.");
  INVALIDARG_EXP_IFN(
    crosslinkerUniverse.meanSquaredBeadDistanceInStrand.empty() ||
      crosslinkerUniverse.meanSquaredBeadDistanceInStrand.size() == nStrands,
    "meanSquaredBeadDistanceInStrand must be empty or match the number of "
    "strands.");

  size_t totalAtoms = nCrosslinks;
  for (size_t s = 0; s < nStrands; ++s) {
    totalAtoms += static_cast<size_t>(
      std::max<int>(0, crosslinkerUniverse.beadsInStrand[s]));
  }

  INVALIDARG_EXP_IFN(minimumNrOfEntanglements < totalAtoms / 2 ||
                       (minimumNrOfEntanglements == 0 && totalAtoms == 0),
                     "Minimum number of slip-links must be less than the "
                     "possible number of slip-links to place.");
  INVALIDARG_EXP_IFN(nrOfEntanglementsToSample < totalAtoms / 2 ||
                       (nrOfEntanglementsToSample == 0 && totalAtoms == 0),
                     "Number of slip-links to place must be less than "
                     "the possible number of slip-links to place.");
  INVALIDARG_EXP_IFN(nrOfEntanglementsToSample >= minimumNrOfEntanglements,
                     "Maximum nr. should be larger than minimum, got " +
                       std::to_string(nrOfEntanglementsToSample) + " and " +
                       std::to_string(minimumNrOfEntanglements) + ".");
  INVALIDARG_EXP_IFN(upperCutoff > 0.0,
                     "Expected a cutoff > 0.0, got " +
                       std::to_string(upperCutoff) + ".");
  INVALIDARG_EXP_IFN(lowerCutoff >= 0.0,
                     "Expected lower cut-off >= 0.0, got " +
                       std::to_string(lowerCutoff) + ".");
  INVALIDARG_EXP_IFN(
    lowerCutoff < upperCutoff,
    "Expected lower cut-off to be smaller than upper cut-off, got " +
      std::to_string(lowerCutoff) + " >= " + std::to_string(upperCutoff) + ".");

  StrandPairEntanglements result;
  if (nrOfEntanglementsToSample == 0 || nStrands == 0) {
    return result;
  }

  std::mt19937 rng;
  if (seed.empty()) {
    std::random_device rd{};
    rng = std::mt19937(rd());
  } else {
    std::seed_seq seedSeq(seed.begin(), seed.end());
    rng = std::mt19937(seedSeq);
  }

  // Topological leaf peeling on CrosslinkerUniverse if filterDanglingAndSoluble
  std::vector<bool> strandIsIgnored(nStrands, false);
  if (filterDanglingAndSoluble) {
    std::vector<std::vector<size_t>> incidentStrands(nCrosslinks);
    std::vector<size_t> activeNonLoopDegree(nCrosslinks, 0);
    std::vector<size_t> primaryLoopCount(nCrosslinks, 0);

    for (size_t s = 0; s < nStrands; ++s) {
      const long int u = crosslinkerUniverse.strandFrom[s];
      const long int v = crosslinkerUniverse.strandTo[s];
      if (u < 0 || v < 0) {
        strandIsIgnored[s] = true;
      } else {
        INVALIDARG_EXP_IFN(
          static_cast<size_t>(u) < nCrosslinks &&
            static_cast<size_t>(v) < nCrosslinks,
          "Crosslinker index out of bounds in strandFrom / strandTo.");
        if (u == v) {
          primaryLoopCount[u] += 1;
        } else {
          incidentStrands[u].push_back(s);
          incidentStrands[v].push_back(s);
          activeNonLoopDegree[u] += 1;
          activeNonLoopDegree[v] += 1;
        }
      }
    }

    std::vector<size_t> peelQueue;
    peelQueue.reserve(nCrosslinks);
    for (size_t xl = 0; xl < nCrosslinks; ++xl) {
      if (activeNonLoopDegree[xl] <= 1 && primaryLoopCount[xl] == 0) {
        peelQueue.push_back(xl);
      }
    }

    size_t qHead = 0;
    while (qHead < peelQueue.size()) {
      const size_t xl = peelQueue[qHead++];
      for (const size_t s : incidentStrands[xl]) {
        if (strandIsIgnored[s]) {
          continue;
        }
        strandIsIgnored[s] = true;
        const size_t u = static_cast<size_t>(crosslinkerUniverse.strandFrom[s]);
        const size_t v = static_cast<size_t>(crosslinkerUniverse.strandTo[s]);
        const size_t other = (u == xl) ? v : u;
        if (activeNonLoopDegree[xl] > 0) {
          activeNonLoopDegree[xl] -= 1;
        }
        if (activeNonLoopDegree[other] > 0) {
          activeNonLoopDegree[other] -= 1;
          if (activeNonLoopDegree[other] <= 1 && primaryLoopCount[other] == 0) {
            peelQueue.push_back(other);
          }
        }
      }
    }
  }

  std::uniform_real_distribution<double> boxXDist(box.getLowX(),
                                                  box.getHighX());
  std::uniform_real_distribution<double> boxYDist(box.getLowY(),
                                                  box.getHighY());
  std::uniform_real_distribution<double> boxZDist(box.getLowZ(),
                                                  box.getHighZ());

  std::vector<StrandModel> strands(nStrands);
  std::vector<size_t> strandInternalOffset(nStrands + 1, 0);
  double maxStrandReach = 0.0;

  for (size_t s = 0; s < nStrands; ++s) {
    StrandModel& sm = strands[s];
    sm.from = crosslinkerUniverse.strandFrom[s];
    sm.to = crosslinkerUniverse.strandTo[s];
    sm.nBeads = std::max<int>(0, crosslinkerUniverse.beadsInStrand[s]);
    sm.isIgnored = strandIsIgnored[s];
    const double b2 =
      (crosslinkerUniverse.meanSquaredBeadDistanceInStrand.size() == nStrands)
        ? crosslinkerUniverse.meanSquaredBeadDistanceInStrand[s]
        : 1.0;
    sm.sigma2Bond = std::max(1e-12, b2 / 3.0);

    if (sm.from >= 0 && sm.to >= 0) {
      sm.nBonds = sm.nBeads + 1;
      sm.nInternalBeads = sm.nBeads;
      Eigen::Vector3d rFrom(crosslinkerUniverse.xlinkX[sm.from],
                            crosslinkerUniverse.xlinkY[sm.from],
                            crosslinkerUniverse.xlinkZ[sm.from]);
      Eigen::Vector3d rTo(crosslinkerUniverse.xlinkX[sm.to],
                          crosslinkerUniverse.xlinkY[sm.to],
                          crosslinkerUniverse.xlinkZ[sm.to]);
      Eigen::Vector3d delta = rTo - rFrom;
      box.handlePBC(delta);
      sm.anchorStart = rFrom;
      sm.anchorEnd = rFrom + delta;
      sm.midpoint = rFrom + 0.5 * delta;
      box.handlePBC(sm.midpoint);
      sm.halfChordLen = 0.5 * delta.norm();
      sm.maxSigma =
        0.5 * std::sqrt(static_cast<double>(sm.nBonds) * sm.sigma2Bond);
      sm.pinned.push_back(PinnedPoint{ 0.0, sm.anchorStart });
      sm.pinned.push_back(PinnedPoint{ 1.0, sm.anchorEnd });
    } else if (sm.from >= 0 || sm.to >= 0) {
      sm.nBonds = sm.nBeads;
      sm.nInternalBeads = std::max<Eigen::Index>(0, sm.nBeads - 1);
      Eigen::Index anchorIdx = sm.from;
      if (sm.from < 0) {
        anchorIdx = sm.to;
        sm.reversedAnchor = true;
      }
      Eigen::Vector3d rAnchor(crosslinkerUniverse.xlinkX[anchorIdx],
                              crosslinkerUniverse.xlinkY[anchorIdx],
                              crosslinkerUniverse.xlinkZ[anchorIdx]);
      sm.anchorStart = rAnchor;
      sm.anchorEnd = rAnchor;
      sm.midpoint = rAnchor;
      box.handlePBC(sm.midpoint);
      sm.halfChordLen = 0.0;
      sm.maxSigma = std::sqrt(static_cast<double>(sm.nBonds) * sm.sigma2Bond);
      sm.pinned.push_back(PinnedPoint{ 0.0, sm.anchorStart });
    } else {
      sm.nBonds = std::max<Eigen::Index>(0, sm.nBeads - 1);
      sm.nInternalBeads = std::max<Eigen::Index>(0, sm.nBeads - 2);
      Eigen::Vector3d rRand(boxXDist(rng), boxYDist(rng), boxZDist(rng));
      sm.anchorStart = rRand;
      sm.anchorEnd = rRand;
      sm.midpoint = rRand;
      sm.halfChordLen = 0.0;
      sm.maxSigma = std::sqrt(static_cast<double>(sm.nBonds) * sm.sigma2Bond);
      sm.pinned.push_back(PinnedPoint{ 0.0, sm.anchorStart });
    }

    strandInternalOffset[s + 1] =
      strandInternalOffset[s] + static_cast<size_t>(sm.nInternalBeads);
    if (sm.nInternalBeads > 0) {
      const double reach = sm.halfChordLen + 4.5 * sm.maxSigma;
      if (reach > maxStrandReach) {
        maxStrandReach = reach;
      }
    }
  }

  const size_t totalInternalBeads = strandInternalOffset[nStrands];
  if (totalInternalBeads < 2) {
    return result;
  }

  // Build spatial grid over strand midpoints
  Eigen::VectorXd strandMidpointCoords(3 * nStrands);
  for (size_t s = 0; s < nStrands; ++s) {
    strandMidpointCoords.segment<3>(3 * s) = strands[s].midpoint;
  }
  const double searchCutoff = maxStrandReach + upperCutoff;
  const pylimer_tools::entities::EigenNeighbourList strandNeighbourList(
    strandMidpointCoords, box, searchCutoff);

  struct BeadRef
  {
    uint32_t strandIdx;
    uint32_t internalIdx;
  };
  std::vector<BeadRef> candidateBeads;
  candidateBeads.reserve(totalInternalBeads);
  for (size_t s = 0; s < nStrands; ++s) {
    const size_t nInt = static_cast<size_t>(strands[s].nInternalBeads);
    for (size_t m = 0; m < nInt; ++m) {
      candidateBeads.push_back(
        BeadRef{ static_cast<uint32_t>(s), static_cast<uint32_t>(m) });
    }
  }

  std::vector<uint8_t> beadIsPaired(totalInternalBeads, 0);
  std::vector<size_t> candidateOrder(totalInternalBeads);
  std::iota(candidateOrder.begin(), candidateOrder.end(), 0);

  // Cached mean positions and standard deviations for all internal beads
  std::vector<Eigen::Vector3d> beadMean(totalInternalBeads);
  std::vector<double> beadSigma(totalInternalBeads);
  const std::function<void(const size_t)> updateStrandBeadCache =
    [&](const size_t s) {
      const StrandModel& sm = strands[s];
      const size_t baseOffset = strandInternalOffset[s];
      const size_t nInt = static_cast<size_t>(sm.nInternalBeads);
      for (size_t m = 0; m < nInt; ++m) {
        const double alpha = sm.internalIdxToAlpha(m);
        sm.evaluateConditionalDistribution(
          alpha, beadMean[baseOffset + m], beadSigma[baseOffset + m]);
      }
    };
  for (size_t s = 0; s < nStrands; ++s) {
    updateStrandBeadCache(s);
  }

  struct CandidatePartner
  {
    size_t globalInternalIdx;
    size_t strandIdx;
    size_t internalIdx;
    double prob;
    double sigma;
    Eigen::Vector3d meanPos;
    Eigen::Vector3d deltaMu;
  };
  std::vector<CandidatePartner> partners;
  partners.reserve(64);
  Eigen::ArrayXi nearbyStrands = Eigen::ArrayXi::Zero(32);

  std::normal_distribution<double> stdNormal(0.0, 1.0);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);

  std::vector<Eigen::Vector3d> sampledPos1;
  std::vector<Eigen::Vector3d> sampledPos2;
  std::vector<Eigen::Vector3d> sampledMidpoints;
  result.strandIdx1.reserve(nrOfEntanglementsToSample);
  result.strandIdx2.reserve(nrOfEntanglementsToSample);
  result.alpha1.reserve(nrOfEntanglementsToSample);
  result.alpha2.reserve(nrOfEntanglementsToSample);
  result.beadIdx1.reserve(nrOfEntanglementsToSample);
  result.beadIdx2.reserve(nrOfEntanglementsToSample);
  sampledPos1.reserve(nrOfEntanglementsToSample);
  sampledPos2.reserve(nrOfEntanglementsToSample);
  sampledMidpoints.reserve(nrOfEntanglementsToSample);

  size_t numEntanglementsSampled = 0;
  size_t numLinksFoundInIteration = 0;

  do {
    numLinksFoundInIteration = 0;
    std::ranges::shuffle(candidateOrder, rng);

    for (const size_t globalIdx1 : candidateOrder) {
      if (beadIsPaired[globalIdx1]) {
        continue;
      }

      const size_t s1 = candidateBeads[globalIdx1].strandIdx;
      const size_t m1 = candidateBeads[globalIdx1].internalIdx;
      StrandModel& sm1 = strands[s1];
      const double alpha1 = sm1.internalIdxToAlpha(m1);

      // 1. Sample 3D position r1 of bead b1 from its conditional Gaussian
      // bridge
      const Eigen::Vector3d& mu1 = beadMean[globalIdx1];
      const double sigma1 = beadSigma[globalIdx1];
      const Eigen::Vector3d r1 =
        mu1 + sigma1 *
                Eigen::Vector3d(stdNormal(rng), stdNormal(rng), stdNormal(rng));

      // Defer pinning (alpha1, r1) on s1 until actually needed:
      // it is only needed if evaluating a same-strand partner (s1 == s2) or
      // upon successfully forming an entanglement.
      bool s1Pinned = false;

      // 2. Query nearby strands around r1
      Eigen::Vector3d r1Wrapped = r1;
      box.handlePBC(r1Wrapped);
      const size_t nNearbyStrands =
        strandNeighbourList.getIndicesCloseToCoordinates(
          nearbyStrands, r1Wrapped, searchCutoff, true);

      partners.clear();
      double survivalProb = 1.0;

      for (size_t nsIdx = 0; nsIdx < nNearbyStrands; ++nsIdx) {
        const size_t s2 = static_cast<size_t>(nearbyStrands[nsIdx]);
        const StrandModel& sm2 = strands[s2];
        const size_t nInt2 = static_cast<size_t>(sm2.nInternalBeads);
        if (nInt2 == 0) {
          continue;
        }

        // Fast strand-level bounding sphere filter
        Eigen::Vector3d midDiff = sm2.midpoint - r1Wrapped;
        box.handlePBC(midDiff);
        const double maxStrandDist =
          sm2.halfChordLen + 4.5 * sm2.maxSigma + upperCutoff;
        if (midDiff.squaredNorm() > maxStrandDist * maxStrandDist) {
          continue;
        }

        const size_t baseOffset2 = strandInternalOffset[s2];
        for (size_t m2 = 0; m2 < nInt2; ++m2) {
          const size_t globalIdx2 = baseOffset2 + m2;
          if (globalIdx2 == globalIdx1 || beadIsPaired[globalIdx2]) {
            continue;
          }

          if (s1 == s2 && sameStrandCutoff > 0.0) {
            const double bondDist =
              std::abs(static_cast<double>(m1) - static_cast<double>(m2));
            // consistent with the discrete samplers: reject only pairs
            // strictly closer than the cutoff
            if (bondDist < sameStrandCutoff) {
              continue;
            }
          }

          Eigen::Vector3d mu2;
          double sigma2 = 0.0;
          if (s1 == s2) {
            if (!s1Pinned) {
              sm1.insertPinnedPoint(alpha1, r1);
              s1Pinned = true;
            }
            const double alpha2 = sm1.internalIdxToAlpha(m2);
            sm1.evaluateConditionalDistribution(alpha2, mu2, sigma2);
          } else {
            mu2 = beadMean[globalIdx2];
            sigma2 = beadSigma[globalIdx2];
          }

          Eigen::Vector3d deltaMu = mu2 - r1;
          box.handlePBC(deltaMu);

          const double maxDist = upperCutoff + 6.0 * sigma2;
          if (deltaMu.squaredNorm() > maxDist * maxDist) {
            continue;
          }
          if (lowerCutoff > 0.0) {
            const double minDist = lowerCutoff - 6.0 * sigma2;
            if (minDist > 0.0 && deltaMu.squaredNorm() < minDist * minDist) {
              continue;
            }
          }
          const double distMu = deltaMu.norm();

          const double p =
            shellProbability3D(lowerCutoff, upperCutoff, distMu, sigma2);
          if (p > 1e-9) {
            partners.push_back(
              CandidatePartner{ globalIdx2, s2, m2, p, sigma2, mu2, deltaMu });
            survivalProb *= (1.0 - p);
          }
        }
      }

      if (partners.empty()) {
        if (s1Pinned) {
          sm1.removePinnedPoint(alpha1);
        }
        continue;
      }

      // 3. Survival probability acceptance test: P(at least one partner in
      // shell)
      const double pAny = 1.0 - survivalProb;
      if (uniform01(rng) >= pAny) {
        if (s1Pinned) {
          sm1.removePinnedPoint(alpha1);
        }
        continue;
      }

      // 4. Select partner b2 with probability proportional to p_c / (1 - 0.5 *
      // p_c)
      double totalWeight = 0.0;
      for (const auto& cand : partners) {
        totalWeight += cand.prob / (1.0 - 0.5 * cand.prob);
      }
      double pick = uniform01(rng) * totalWeight;
      size_t chosenIdx = partners.size() - 1;
      for (size_t idx = 0; idx < partners.size(); ++idx) {
        pick -= partners[idx].prob / (1.0 - 0.5 * partners[idx].prob);
        if (pick <= 0.0) {
          chosenIdx = idx;
          break;
        }
      }

      const CandidatePartner& chosen = partners[chosenIdx];
      const size_t s2 = chosen.strandIdx;
      const size_t m2 = chosen.internalIdx;
      StrandModel& sm2 = strands[s2];
      const double alpha2 = sm2.internalIdxToAlpha(m2);

      // 5. Sample exact displacement d = r2 - r1 within [lowerCutoff,
      // upperCutoff]
      const Eigen::Vector3d d = sampleDisplacementInShell3D(
        lowerCutoff, upperCutoff, chosen.deltaMu, chosen.sigma, rng);

      // Ensure (alpha1, r1) is pinned on s1
      if (!s1Pinned) {
        sm1.insertPinnedPoint(alpha1, r1);
        s1Pinned = true;
      }

      // Pin (alpha2, r2) in strand s2's unwrapped coordinate frame
      const Eigen::Vector3d r2InS2Frame = (chosen.meanPos - chosen.deltaMu) + d;
      sm2.insertPinnedPoint(alpha2, r2InS2Frame);

      beadIsPaired[globalIdx1] = 1;
      beadIsPaired[chosen.globalInternalIdx] = 1;

      updateStrandBeadCache(s1);
      if (s2 != s1) {
        updateStrandBeadCache(s2);
      }

      if (!filterDanglingAndSoluble || (!sm1.isIgnored && !sm2.isIgnored)) {
        result.strandIdx1.push_back(s1);
        result.strandIdx2.push_back(s2);
        result.alpha1.push_back(alpha1);
        result.alpha2.push_back(alpha2);
        result.beadIdx1.push_back(sm1.internalIdxToBeadIdx(m1));
        result.beadIdx2.push_back(sm2.internalIdxToBeadIdx(m2));

        const Eigen::Vector3d pos1 = r1;
        const Eigen::Vector3d pos2 = r1 + d;
        Eigen::Vector3d mid = r1 + 0.5 * d;
        box.handlePBC(mid);
        sampledPos1.push_back(pos1);
        sampledPos2.push_back(pos2);
        sampledMidpoints.push_back(mid);
      }

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

  const size_t nRecorded = result.strandIdx1.size();
  result.x.resize(nRecorded);
  result.y.resize(nRecorded);
  result.z.resize(nRecorded);
  result.position1.resize(3 * nRecorded);
  result.position2.resize(3 * nRecorded);
  for (size_t e = 0; e < nRecorded; ++e) {
    result.x[e] = sampledMidpoints[e].x();
    result.y[e] = sampledMidpoints[e].y();
    result.z[e] = sampledMidpoints[e].z();
    result.position1.segment<3>(3 * e) = sampledPos1[e];
    result.position2.segment<3>(3 * e) = sampledPos2[e];
  }

  return result;
}
}
