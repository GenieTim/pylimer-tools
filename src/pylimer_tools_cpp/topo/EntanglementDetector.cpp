#include "EntanglementDetector.h"

#include "../entities/EigenNeighbourList.h"
#include "../entities/Molecule.h"
#include "../entities/Universe.h"
#include "../utils/RandomWalker.h"
#include "../utils/StringUtils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <list>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace pylimer_tools::topo::entanglement_detection {

namespace {

  void validateEntanglementParameters(const size_t totalAtoms,
                                      const size_t nrOfEntanglementsToSample,
                                      const size_t minimumNrOfEntanglements,
                                      const double upperCutoff,
                                      const double lowerCutoff)
  {
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
        std::to_string(lowerCutoff) + " >= " + std::to_string(upperCutoff) +
        ".");
  }

  std::mt19937 initializeRng(const std::string& seed)
  {
    if (seed.empty()) {
      std::random_device rd{};
      return std::mt19937(rd());
    }
    std::seed_seq seedSeq(seed.begin(), seed.end());
    return std::mt19937(seedSeq);
  }

  struct StrandTopology
  {
    std::vector<int32_t> strandId;
    std::vector<uint32_t> posInStrand;
    std::vector<uint32_t> distToEnd;
    std::vector<uint32_t> strandLength;
    std::vector<uint8_t> strandIsRing;
  };

  /**
   * @brief Label every degree-2 bead with its strand ID, position along the
   * strand, and distance to the nearest non-degree-2 vertex in O(N_atoms).
   */
  StrandTopology buildStrandTopology(const igraph_t& graph,
                                     const std::vector<int>& vertexDegrees)
  {
    const size_t nAtoms = vertexDegrees.size();
    StrandTopology topo;
    topo.strandId.assign(nAtoms, -1);
    topo.posInStrand.assign(nAtoms, 0);
    topo.distToEnd.assign(nAtoms, 0);

    igraph_vector_int_t neighbors;
    igraph_vector_int_init(&neighbors, 8);
    igraph_vector_int_t deg2Neighbors;
    igraph_vector_int_init(&deg2Neighbors, 2);

    // 1. Walk linear strands starting from non-degree-2 vertices
    for (size_t v = 0; v < nAtoms; ++v) {
      if (vertexDegrees[v] == 2) {
        continue;
      }
      REQUIRE_IGRAPH_SUCCESS(igraph_neighbors(&graph,
                                              &neighbors,
                                              static_cast<igraph_integer_t>(v),
                                              IGRAPH_ALL,
                                              IGRAPH_LOOPS_TWICE,
                                              true));
      const igraph_integer_t nDeg = igraph_vector_int_size(&neighbors);
      for (igraph_integer_t k = 0; k < nDeg; ++k) {
        const size_t startNb =
          static_cast<size_t>(igraph_vector_int_get(&neighbors, k));
        if (vertexDegrees[startNb] != 2 || topo.strandId[startNb] != -1) {
          continue;
        }
        const int32_t sid = static_cast<int32_t>(topo.strandLength.size());
        size_t prev = v;
        size_t curr = startNb;
        uint32_t pos = 0;
        while (vertexDegrees[curr] == 2 && topo.strandId[curr] == -1) {
          topo.strandId[curr] = sid;
          topo.posInStrand[curr] = pos++;
          REQUIRE_IGRAPH_SUCCESS(
            igraph_neighbors(&graph,
                             &deg2Neighbors,
                             static_cast<igraph_integer_t>(curr),
                             IGRAPH_ALL,
                             IGRAPH_LOOPS_TWICE,
                             true));
          const size_t n0 =
            static_cast<size_t>(igraph_vector_int_get(&deg2Neighbors, 0));
          const size_t n1 =
            static_cast<size_t>(igraph_vector_int_get(&deg2Neighbors, 1));
          const size_t next = (n0 != prev) ? n0 : n1;
          prev = curr;
          curr = next;
        }
        topo.strandLength.push_back(pos);
        topo.strandIsRing.push_back(0);
      }
    }

    // 2. Walk isolated cyclic rings (components where all vertices have degree
    // 2)
    for (size_t v = 0; v < nAtoms; ++v) {
      if (vertexDegrees[v] != 2 || topo.strandId[v] != -1) {
        continue;
      }
      const int32_t sid = static_cast<int32_t>(topo.strandLength.size());
      size_t prev = v;
      size_t curr = v;
      uint32_t pos = 0;
      while (topo.strandId[curr] == -1) {
        topo.strandId[curr] = sid;
        topo.posInStrand[curr] = pos++;
        REQUIRE_IGRAPH_SUCCESS(
          igraph_neighbors(&graph,
                           &deg2Neighbors,
                           static_cast<igraph_integer_t>(curr),
                           IGRAPH_ALL,
                           IGRAPH_LOOPS_TWICE,
                           true));
        const size_t n0 =
          static_cast<size_t>(igraph_vector_int_get(&deg2Neighbors, 0));
        const size_t n1 =
          static_cast<size_t>(igraph_vector_int_get(&deg2Neighbors, 1));
        const size_t next =
          (pos == 1) ? ((n0 != curr) ? n0 : n1) : ((n0 != prev) ? n0 : n1);
        prev = curr;
        curr = next;
      }
      topo.strandLength.push_back(pos);
      topo.strandIsRing.push_back(1);
    }

    igraph_vector_int_destroy(&deg2Neighbors);
    igraph_vector_int_destroy(&neighbors);

    // 3. Populate distToEnd for every degree-2 vertex
    constexpr uint32_t kRingDist = std::numeric_limits<uint32_t>::max() / 2;
    for (size_t v = 0; v < nAtoms; ++v) {
      const int32_t sid = topo.strandId[v];
      if (sid < 0) {
        continue;
      }
      if (topo.strandIsRing[static_cast<size_t>(sid)] != 0) {
        topo.distToEnd[v] = kRingDist;
      } else {
        const uint32_t pos = topo.posInStrand[v];
        const uint32_t len = topo.strandLength[static_cast<size_t>(sid)];
        topo.distToEnd[v] = std::min(pos + 1, len - pos);
      }
    }

    return topo;
  }

  /**
   * @brief Stamp all vertices within <= maxDepth bonds of source using an
   * epoch counter to avoid clearing bfsVisited between queries.
   */
  void runBoundedBfs(const igraph_t& graph,
                     const size_t source,
                     const size_t maxDepth,
                     std::vector<uint32_t>& bfsVisited,
                     uint32_t& bfsEpoch,
                     std::vector<std::pair<size_t, size_t>>& bfsQueue,
                     igraph_vector_int_t& neighbors)
  {
    ++bfsEpoch;
    if (bfsEpoch == 0) {
      std::fill(bfsVisited.begin(), bfsVisited.end(), 0);
      bfsEpoch = 1;
    }
    bfsQueue.clear();
    bfsVisited[source] = bfsEpoch;
    bfsQueue.emplace_back(source, 0);

    size_t head = 0;
    while (head < bfsQueue.size()) {
      const size_t curr = bfsQueue[head].first;
      const size_t depth = bfsQueue[head].second;
      ++head;
      if (depth >= maxDepth) {
        continue;
      }
      REQUIRE_IGRAPH_SUCCESS(
        igraph_neighbors(&graph,
                         &neighbors,
                         static_cast<igraph_integer_t>(curr),
                         IGRAPH_ALL,
                         IGRAPH_LOOPS_TWICE,
                         true));
      const igraph_integer_t nDeg = igraph_vector_int_size(&neighbors);
      for (igraph_integer_t k = 0; k < nDeg; ++k) {
        const size_t nb =
          static_cast<size_t>(igraph_vector_int_get(&neighbors, k));
        if (bfsVisited[nb] != bfsEpoch) {
          bfsVisited[nb] = bfsEpoch;
          bfsQueue.emplace_back(nb, depth + 1);
        }
      }
    }
  }

} // anonymous namespace

AtomPairEntanglements
randomlyFindEntanglements(const pylimer_tools::entities::Universe& universe,
                          const size_t nrOfEntanglementsToSample,
                          const double upperCutoff,
                          const double lowerCutoff,
                          const size_t minimumNrOfEntanglements,
                          const double sameStrandCutoff,
                          const std::string& seed,
                          const int crossLinkerType,
                          const bool ignoreCrosslinks,
                          const EntanglementFilteringMode filteringMode)
{
  validateEntanglementParameters(universe.getNrOfAtoms(),
                                 nrOfEntanglementsToSample,
                                 minimumNrOfEntanglements,
                                 upperCutoff,
                                 lowerCutoff);

  std::vector<std::pair<size_t, size_t>> pairsOfAtoms;
  pairsOfAtoms.reserve(nrOfEntanglementsToSample);
  std::vector<long int> pairOfAtom(universe.getNrOfAtoms(), -1);

  if (nrOfEntanglementsToSample == 0 || universe.getNrOfAtoms() < 2) {
    return AtomPairEntanglements(pairsOfAtoms, pairOfAtom);
  }

  std::mt19937 rng = initializeRng(seed);

  const pylimer_tools::entities::Box box = universe.getBox();
  const Eigen::VectorXd coordinates = universe.getVertexCoordinates();
  const pylimer_tools::entities::EigenNeighbourList neighbourList(
    coordinates, box, upperCutoff);
  std::vector<uint8_t> vertexIsEligible(universe.getNrOfAtoms(), 1);
  const std::vector<int> vertexTypes = universe.getPropertyValues<int>("type");
  const std::vector<int> vertexDegrees = universe.getVertexDegrees();

  const size_t cutoffCeil = sameStrandCutoff > 0.0
                              ? static_cast<size_t>(std::ceil(sameStrandCutoff))
                              : 0;
  StrandTopology strandTopo;
  std::vector<uint32_t> bfsVisited;
  uint32_t bfsEpoch = 0;
  std::vector<std::pair<size_t, size_t>> bfsQueue;
  igraph_vector_int_t bfsNeighbors;
  if (cutoffCeil > 0) {
    strandTopo = buildStrandTopology(universe.getGraph(), vertexDegrees);
    bfsVisited.assign(universe.getNrOfAtoms(), 0);
    igraph_vector_int_init(&bfsNeighbors, 8);
  }

  std::vector<uint8_t> activeAtoms;
  if (filteringMode != EntanglementFilteringMode::NONE) {
    const std::vector<bool> activeAtomsBool =
      universe.identifyActiveAtoms(crossLinkerType);
    activeAtoms.resize(activeAtomsBool.size());
    for (size_t i = 0; i < activeAtomsBool.size(); ++i) {
      activeAtoms[i] = activeAtomsBool[i] ? 1 : 0;
    }
  }

  std::vector<long int> vertexIdsToSample;
  vertexIdsToSample.reserve(universe.getNrOfAtoms());

  for (size_t i = 0; i < universe.getNrOfAtoms(); ++i) {
    if (vertexDegrees[i] != 2) {
      vertexIsEligible[i] = 0;
    }
    if (ignoreCrosslinks && vertexTypes[i] == crossLinkerType) {
      vertexIsEligible[i] = 0;
    }
    if (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
        !activeAtoms.empty() && activeAtoms[i] == 0) {
      vertexIsEligible[i] = 0;
    }
    if (vertexIsEligible[i] != 0) {
      vertexIdsToSample.push_back(static_cast<long int>(i));
    }
  }

  size_t numLinksFoundInIteration = 1;
  Eigen::ArrayXi neighbourIndices = Eigen::ArrayXi::Zero(12);
  size_t numEntanglementsSampled = 0;

  do {
    numLinksFoundInIteration = 0;
    for (size_t sampleIdx = 0; sampleIdx < vertexIdsToSample.size();
         ++sampleIdx) {
      const size_t randSampleIdx = std::uniform_int_distribution<size_t>(
        sampleIdx, vertexIdsToSample.size() - 1)(rng);
      std::swap(vertexIdsToSample[sampleIdx], vertexIdsToSample[randSampleIdx]);
      const long int atomVertexIdx1 = vertexIdsToSample[sampleIdx];
      const size_t atomIdx1 = static_cast<size_t>(atomVertexIdx1);
      if (pairOfAtom[atomIdx1] != -1) {
        continue;
      }
      vertexIsEligible[atomIdx1] = 0;

      const Eigen::Index nNeighbours =
        static_cast<Eigen::Index>(neighbourList.getIndicesCloseToCoordinates(
          neighbourIndices,
          coordinates.segment<3>(3 * atomVertexIdx1),
          upperCutoff));
      long int atomVertexIdx2 = -1;
      bool hasBfsForSource = false;

      for (Eigen::Index idxInNeighbours = 0; idxInNeighbours < nNeighbours;
           ++idxInNeighbours) {
        const Eigen::Index randNeighbourIdx =
          std::uniform_int_distribution<Eigen::Index>(idxInNeighbours,
                                                      nNeighbours - 1)(rng);
        std::swap(neighbourIndices[idxInNeighbours],
                  neighbourIndices[randNeighbourIdx]);
        const long int candidateVertex = neighbourIndices[idxInNeighbours];
        const size_t candIdx = static_cast<size_t>(candidateVertex);
        if (vertexIsEligible[candIdx] == 0 || pairOfAtom[candIdx] != -1) {
          continue;
        }

        Eigen::Vector3d distance = coordinates.segment<3>(3 * candidateVertex) -
                                   coordinates.segment<3>(3 * atomVertexIdx1);
        box.handlePBC(distance);

        const double dist2 = distance.squaredNorm();
        if (dist2 < lowerCutoff * lowerCutoff ||
            dist2 > upperCutoff * upperCutoff) {
          continue;
        }

        bool pathValid = true;
        if (cutoffCeil > 0) {
          const size_t u = static_cast<size_t>(atomVertexIdx1);
          const size_t v = static_cast<size_t>(candidateVertex);
          const int32_t sidU = strandTopo.strandId[u];
          const int32_t sidV = strandTopo.strandId[v];
          if (sidU == sidV && sidU >= 0) {
            const uint32_t posU = strandTopo.posInStrand[u];
            const uint32_t posV = strandTopo.posInStrand[v];
            const uint32_t diff = (posU > posV) ? (posU - posV) : (posV - posU);
            if (strandTopo.strandIsRing[static_cast<size_t>(sidU)] != 0) {
              const uint32_t len =
                strandTopo.strandLength[static_cast<size_t>(sidU)];
              const uint32_t ringDist = std::min(diff, len - diff);
              pathValid = (ringDist >= cutoffCeil);
            } else if (diff < cutoffCeil) {
              pathValid = false;
            } else if (static_cast<uint64_t>(strandTopo.distToEnd[u]) +
                         static_cast<uint64_t>(strandTopo.distToEnd[v]) <
                       cutoffCeil) {
              if (!hasBfsForSource) {
                runBoundedBfs(universe.getGraph(),
                              u,
                              cutoffCeil - 1,
                              bfsVisited,
                              bfsEpoch,
                              bfsQueue,
                              bfsNeighbors);
                hasBfsForSource = true;
              }
              pathValid = (bfsVisited[v] != bfsEpoch);
            }
          } else if (static_cast<uint64_t>(strandTopo.distToEnd[u]) +
                       static_cast<uint64_t>(strandTopo.distToEnd[v]) <
                     cutoffCeil) {
            if (!hasBfsForSource) {
              runBoundedBfs(universe.getGraph(),
                            u,
                            cutoffCeil - 1,
                            bfsVisited,
                            bfsEpoch,
                            bfsQueue,
                            bfsNeighbors);
              hasBfsForSource = true;
            }
            pathValid = (bfsVisited[v] != bfsEpoch);
          }
        }
        if (pathValid) {
          atomVertexIdx2 = candidateVertex;
          break;
        }
      }

      if (atomVertexIdx2 < 0) {
        continue;
      }
      const size_t atomIdx2 = static_cast<size_t>(atomVertexIdx2);

      RUNTIME_EXP_IFN(pairOfAtom[atomIdx2] == -1,
                      "Expected not to be able to sample the same atom twice.");

      if (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT) {
        pairOfAtom[atomIdx2] = static_cast<long int>(pairsOfAtoms.size());
        pairOfAtom[atomIdx1] = static_cast<long int>(pairsOfAtoms.size());
        pairsOfAtoms.emplace_back(universe.getAtomIdByIdx(atomVertexIdx1),
                                  universe.getAtomIdByIdx(atomVertexIdx2));
      } else if (filteringMode == EntanglementFilteringMode::OMIT_INACTIVE) {
        if (activeAtoms.empty() ||
            (activeAtoms[atomIdx1] != 0 && activeAtoms[atomIdx2] != 0)) {
          pairOfAtom[atomIdx2] = static_cast<long int>(pairsOfAtoms.size());
          pairOfAtom[atomIdx1] = static_cast<long int>(pairsOfAtoms.size());
          pairsOfAtoms.emplace_back(universe.getAtomIdByIdx(atomVertexIdx1),
                                    universe.getAtomIdByIdx(atomVertexIdx2));
        } else {
          pairOfAtom[atomIdx2] = -2;
          pairOfAtom[atomIdx1] = -2;
        }
      } else { // NONE
        pairOfAtom[atomIdx2] = static_cast<long int>(pairsOfAtoms.size());
        pairOfAtom[atomIdx1] = static_cast<long int>(pairsOfAtoms.size());
        pairsOfAtoms.emplace_back(universe.getAtomIdByIdx(atomVertexIdx1),
                                  universe.getAtomIdByIdx(atomVertexIdx2));
      }

      vertexIsEligible[atomIdx2] = 0;
      numEntanglementsSampled += 1;
      numLinksFoundInIteration += 1;

      const size_t targetCount =
        (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT)
          ? pairsOfAtoms.size()
          : numEntanglementsSampled;
      if (targetCount >= nrOfEntanglementsToSample) {
        break;
      }
    }

    const size_t targetCount =
      (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT)
        ? pairsOfAtoms.size()
        : numEntanglementsSampled;
    if (targetCount >= nrOfEntanglementsToSample) {
      break;
    }
  } while (((filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT)
              ? pairsOfAtoms.size() < minimumNrOfEntanglements
              : numEntanglementsSampled < minimumNrOfEntanglements) &&
           numLinksFoundInIteration > 0);

  for (long int& i : pairOfAtom) {
    if (i == -2) {
      i = -1;
    }
  }

  if (cutoffCeil > 0) {
    igraph_vector_int_destroy(&bfsNeighbors);
  }

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

TopologicalLeafPeelingResult
peelTopologicalLeavesDetailed(
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse)
{
  const size_t nCrosslinks = crosslinkerUniverse.xlinkX.size();
  const size_t nStrands = crosslinkerUniverse.beadsInStrand.size();
  TopologicalLeafPeelingResult result;
  result.strandIsIgnored.assign(nStrands, false);
  result.crosslinkerIsIgnored.assign(nCrosslinks, false);
  std::vector<std::vector<size_t>> incidentStrands(nCrosslinks);
  std::vector<size_t> activeNonLoopDegree(nCrosslinks, 0);
  std::vector<size_t> primaryLoopCount(nCrosslinks, 0);

  for (size_t s = 0; s < nStrands; ++s) {
    const long int u = crosslinkerUniverse.strandFrom[s];
    const long int v = crosslinkerUniverse.strandTo[s];
    if (u < 0 || v < 0) {
      result.strandIsIgnored[s] = true;
    } else {
      const size_t uIdx = static_cast<size_t>(u);
      const size_t vIdx = static_cast<size_t>(v);
      INVALIDARG_EXP_IFN(
        uIdx < nCrosslinks && vIdx < nCrosslinks,
        "Crosslinker index out of bounds in strandFrom / strandTo.");
      if (uIdx == vIdx) {
        primaryLoopCount[uIdx] += 1;
      } else {
        incidentStrands[uIdx].push_back(s);
        incidentStrands[vIdx].push_back(s);
        activeNonLoopDegree[uIdx] += 1;
        activeNonLoopDegree[vIdx] += 1;
      }
    }
  }

  std::vector<size_t> peelQueue;
  peelQueue.reserve(nCrosslinks);
  for (size_t xl = 0; xl < nCrosslinks; ++xl) {
    if (activeNonLoopDegree[xl] <= 1 && primaryLoopCount[xl] == 0) {
      peelQueue.push_back(xl);
      result.crosslinkerIsIgnored[xl] = true;
    }
  }

  size_t qHead = 0;
  while (qHead < peelQueue.size()) {
    const size_t xl = peelQueue[qHead++];
    for (const size_t s : incidentStrands[xl]) {
      if (result.strandIsIgnored[s]) {
        continue;
      }
      result.strandIsIgnored[s] = true;
      const size_t u = static_cast<size_t>(crosslinkerUniverse.strandFrom[s]);
      const size_t v = static_cast<size_t>(crosslinkerUniverse.strandTo[s]);
      const size_t other = (u == xl) ? v : u;
      if (activeNonLoopDegree[xl] > 0) {
        activeNonLoopDegree[xl] -= 1;
      }
      if (activeNonLoopDegree[other] > 0) {
        activeNonLoopDegree[other] -= 1;
        if (activeNonLoopDegree[other] <= 1 && primaryLoopCount[other] == 0 &&
            !result.crosslinkerIsIgnored[other]) {
          peelQueue.push_back(other);
          result.crosslinkerIsIgnored[other] = true;
        }
      }
    }
  }
  return result;
}

std::vector<bool>
peelTopologicalLeaves(
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse)
{
  return peelTopologicalLeavesDetailed(crosslinkerUniverse).strandIsIgnored;
}

namespace {

  struct StrandBoundingBox
  {
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    Eigen::Vector3d halfExtents = Eigen::Vector3d::Zero();
    double radius = 0.0;
    uint64_t seed = 0;
    uint32_t nInternalBeads = 0;
    uint32_t firstBeadOffset = 0;
    uint32_t nBonds = 1;
    bool isPrimaryLoop = false;
    bool isIgnored = false;
  };

  inline bool isAABBWithinCutoff(const Eigen::Vector3d& point,
                                 const Eigen::Vector3d& aabbCenter,
                                 const Eigen::Vector3d& aabbHalfExtents,
                                 const pylimer_tools::entities::Box& box,
                                 const double cutoff)
  {
    Eigen::Vector3d delta = aabbCenter - point;
    box.handlePBC(delta);
    double distSq = 0.0;
    for (int d = 0; d < 3; ++d) {
      const double excess = std::abs(delta[d]) - aabbHalfExtents[d];
      if (excess > 0.0) {
        distSq += excess * excess;
      }
    }
    return distSq <= cutoff * cutoff;
  }

  /**
   * @brief Bead-budgeted LRU cache for regenerated strand coordinates.
   */
  class StrandCache
  {
  public:
    explicit StrandCache(const size_t maxBeads = 250000)
      : maxBeads_(std::max<size_t>(1, maxBeads))
    {
    }

    template<typename Generator>
    const Eigen::VectorXd& get(const size_t strandIdx, Generator&& gen)
    {
      const auto it = indexMap_.find(strandIdx);
      if (it != indexMap_.end()) {
        if (it->second != entries_.begin()) {
          entries_.splice(entries_.begin(), entries_, it->second);
        }
        return entries_.front().coords;
      }
      Eigen::VectorXd coords = gen();
      const size_t nBeads = static_cast<size_t>(coords.size() / 3);
      while (cachedBeads_ + nBeads > maxBeads_ && !entries_.empty()) {
        cachedBeads_ -= static_cast<size_t>(entries_.back().coords.size() / 3);
        indexMap_.erase(entries_.back().strandIdx);
        entries_.pop_back();
      }
      cachedBeads_ += nBeads;
      entries_.push_front(Entry{ strandIdx, std::move(coords) });
      indexMap_[strandIdx] = entries_.begin();
      return entries_.front().coords;
    }

    void putIfFits(const size_t strandIdx, Eigen::VectorXd&& coords)
    {
      if (indexMap_.find(strandIdx) != indexMap_.end()) {
        return;
      }
      const size_t nBeads = static_cast<size_t>(coords.size() / 3);
      if (cachedBeads_ + nBeads > maxBeads_) {
        return;
      }
      cachedBeads_ += nBeads;
      entries_.push_front(Entry{ strandIdx, std::move(coords) });
      indexMap_[strandIdx] = entries_.begin();
    }

  private:
    struct Entry
    {
      size_t strandIdx;
      Eigen::VectorXd coords;
    };
    size_t maxBeads_;
    size_t cachedBeads_ = 0;
    std::list<Entry> entries_;
    std::unordered_map<size_t, std::list<Entry>::iterator> indexMap_;
  };

  /**
   * @brief Periodic 3D spatial grid storing strand indices in all cells
   *        overlapped by each strand's padded AABB.
   *
   * Eliminates global max-radius sensitivity and small-box 1x1x1 bucket
   * collapse: a point query only inspects the single cell containing the point.
   */
  class StrandSpatialGrid
  {
  public:
    StrandSpatialGrid(const pylimer_tools::entities::Box& box,
                      const double targetCellSize,
                      const size_t nStrands)
      : boxLow_(box.getLowX(), box.getLowY(), box.getLowZ())
      , boxL_(box.getL())
    {
      const double safeCellSize = std::max(1e-3, targetCellSize);
      for (int d = 0; d < 3; ++d) {
        nCells_[d] = std::max<size_t>(
          1, static_cast<size_t>(std::floor(boxL_[d] / safeCellSize)));
      }
      const size_t totalCells = nCells_[0] * nCells_[1] * nCells_[2];
      const size_t maxCells = std::max<size_t>(4096, 8 * nStrands);
      if (totalCells > maxCells) {
        const double scale = std::cbrt(static_cast<double>(totalCells) /
                                       static_cast<double>(maxCells));
        for (int d = 0; d < 3; ++d) {
          nCells_[d] =
            std::max<size_t>(1,
                             static_cast<size_t>(std::floor(
                               static_cast<double>(nCells_[d]) / scale)));
        }
      }
      for (int d = 0; d < 3; ++d) {
        cellWidth_[d] = boxL_[d] / static_cast<double>(nCells_[d]);
      }
      cells_.resize(nCells_[0] * nCells_[1] * nCells_[2]);
    }

    void insertAABB(const size_t strandIdx,
                    const Eigen::Vector3d& center,
                    const Eigen::Vector3d& halfExtents,
                    const double pad)
    {
      const Eigen::Vector3d effHalf = halfExtents.array() + pad + 1e-6;
      long int lo[3];
      long int hi[3];
      for (int d = 0; d < 3; ++d) {
        const double minRel = center[d] - effHalf[d] - boxLow_[d];
        const double maxRel = center[d] + effHalf[d] - boxLow_[d];
        lo[d] = static_cast<long int>(std::floor(minRel / cellWidth_[d]));
        hi[d] = static_cast<long int>(std::floor(maxRel / cellWidth_[d]));
        if (static_cast<size_t>(hi[d] - lo[d] + 1) >= nCells_[d]) {
          lo[d] = 0;
          hi[d] = static_cast<long int>(nCells_[d]) - 1;
        }
      }
      for (long int ix = lo[0]; ix <= hi[0]; ++ix) {
        const size_t wx = wrapIndex(ix, nCells_[0]);
        for (long int iy = lo[1]; iy <= hi[1]; ++iy) {
          const size_t wy = wrapIndex(iy, nCells_[1]);
          for (long int iz = lo[2]; iz <= hi[2]; ++iz) {
            const size_t wz = wrapIndex(iz, nCells_[2]);
            cells_[flatIndex(wx, wy, wz)].push_back(strandIdx);
          }
        }
      }
    }

    size_t cellIndexForPoint(const Eigen::Vector3d& pt) const
    {
      size_t idx[3];
      for (int d = 0; d < 3; ++d) {
        const long int raw = static_cast<long int>(
          std::floor((pt[d] - boxLow_[d]) / cellWidth_[d]));
        idx[d] = wrapIndex(raw, nCells_[d]);
      }
      return flatIndex(idx[0], idx[1], idx[2]);
    }

    const std::vector<size_t>& queryPoint(const Eigen::Vector3d& pt) const
    {
      return cells_[cellIndexForPoint(pt)];
    }

    size_t numCells() const { return cells_.size(); }

  private:
    static size_t wrapIndex(const long int idx, const size_t n)
    {
      const long int m = static_cast<long int>(n);
      const long int r = idx % m;
      return static_cast<size_t>(r < 0 ? r + m : r);
    }

    size_t flatIndex(const size_t x, const size_t y, const size_t z) const
    {
      return (x * nCells_[1] + y) * nCells_[2] + z;
    }

    Eigen::Vector3d boxLow_;
    Eigen::Vector3d boxL_;
    size_t nCells_[3] = { 1, 1, 1 };
    double cellWidth_[3] = { 1.0, 1.0, 1.0 };
    std::vector<std::vector<size_t>> cells_;
  };

} // anonymous namespace

Eigen::VectorXd
sampleStrandBeadsOnDemand(
  const pylimer_tools::entities::Box& box,
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
  const size_t strandIdx,
  const uint64_t strandSeed,
  const int nMcSteps)
{
  const int nBeads = cu.beadsInStrand[strandIdx];
  if (nBeads <= 0) {
    return Eigen::VectorXd(0);
  }
  const double beadDist = cu.beadDistanceInStrand.size() > strandIdx
                            ? cu.beadDistanceInStrand[strandIdx]
                            : 1.0;
  const double meanSqDist =
    cu.meanSquaredBeadDistanceInStrand.size() > strandIdx
      ? cu.meanSquaredBeadDistanceInStrand[strandIdx]
      : 1.0;

  std::mt19937 strandRng(static_cast<std::mt19937::result_type>(strandSeed));
  const long int from = cu.strandFrom[strandIdx];
  const long int to = cu.strandTo[strandIdx];

  if (from >= 0 && to >= 0) {
    const size_t fromIdx = static_cast<size_t>(from);
    const size_t toIdx = static_cast<size_t>(to);
    const Eigen::Vector3d rFrom(
      cu.xlinkX[fromIdx], cu.xlinkY[fromIdx], cu.xlinkZ[fromIdx]);
    const Eigen::Vector3d rTo(
      cu.xlinkX[toIdx], cu.xlinkY[toIdx], cu.xlinkZ[toIdx]);
    Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChainFromTo(
      box, rFrom, rTo, nBeads, beadDist, meanSqDist, strandRng, true);
    if (nMcSteps > 0) {
      pylimer_tools::sim::equilibrateChainWithMC(positions,
                                                 meanSqDist,
                                                 strandRng,
                                                 true,
                                                 true,
                                                 static_cast<size_t>(nMcSteps));
    }
    return positions.segment(3, 3 * nBeads);
  } else if (from >= 0 || to >= 0) {
    const size_t anchor = static_cast<size_t>(from >= 0 ? from : to);
    const Eigen::Vector3d rAnchor(
      cu.xlinkX[anchor], cu.xlinkY[anchor], cu.xlinkZ[anchor]);
    Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
      nBeads, beadDist, meanSqDist, strandRng);
    if (nMcSteps > 0) {
      pylimer_tools::sim::equilibrateChainWithMC(positions,
                                                 meanSqDist,
                                                 strandRng,
                                                 true,
                                                 false,
                                                 static_cast<size_t>(nMcSteps));
    }
    positions.reshaped(3, nBeads).colwise() += rAnchor;
    return positions;
  } else {
    std::uniform_real_distribution<double> distX(box.getLowX(), box.getHighX());
    std::uniform_real_distribution<double> distY(box.getLowY(), box.getHighY());
    std::uniform_real_distribution<double> distZ(box.getLowZ(), box.getHighZ());
    const Eigen::Vector3d rFrom(
      distX(strandRng), distY(strandRng), distZ(strandRng));
    Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
      nBeads, beadDist, meanSqDist, strandRng);
    if (nMcSteps > 0) {
      pylimer_tools::sim::equilibrateChainWithMC(positions,
                                                 meanSqDist,
                                                 strandRng,
                                                 true,
                                                 false,
                                                 static_cast<size_t>(nMcSteps));
    }
    positions.reshaped(3, nBeads).colwise() += rFrom;
    return positions;
  }
}

StrandPairEntanglements
randomlyFindEntanglementsOnDemand(
  const pylimer_tools::entities::Box& box,
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& crosslinkerUniverse,
  const size_t nrOfEntanglementsToSample,
  const double upperCutoff,
  const double lowerCutoff,
  const size_t minimumNrOfEntanglements,
  const double sameStrandCutoff,
  const std::string& seed,
  const EntanglementFilteringMode filteringMode,
  const int nMcSteps)
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

  validateEntanglementParameters(totalAtoms,
                                 nrOfEntanglementsToSample,
                                 minimumNrOfEntanglements,
                                 upperCutoff,
                                 lowerCutoff);

  StrandPairEntanglements result;
  if (nrOfEntanglementsToSample == 0 || nStrands == 0) {
    return result;
  }

  std::vector<bool> strandIsIgnored(nStrands, false);
  if (filteringMode != EntanglementFilteringMode::NONE) {
    strandIsIgnored = peelTopologicalLeaves(crosslinkerUniverse);
  }

  std::mt19937 masterRng = initializeRng(seed);

  std::vector<uint64_t> strandSeeds(nStrands);
  if (seed.empty() && crosslinkerUniverse.strandSeeds.size() == nStrands) {
    strandSeeds = crosslinkerUniverse.strandSeeds;
  } else {
    std::uniform_int_distribution<uint64_t> seedDist(
      0, std::numeric_limits<uint64_t>::max());
    for (size_t s = 0; s < nStrands; ++s) {
      strandSeeds[s] = seedDist(masterRng);
    }
  }

  std::vector<StrandBoundingBox> bboxes(nStrands);
  std::vector<size_t> strandInternalOffset(nStrands + 1, 0);
  StrandCache cache(250000);

  for (size_t s = 0; s < nStrands; ++s) {
    const int nBeads = crosslinkerUniverse.beadsInStrand[s];
    const long int from = crosslinkerUniverse.strandFrom[s];
    const long int to = crosslinkerUniverse.strandTo[s];

    int nInternalBeads = 0;
    int firstBeadOffset = 0;
    int nBonds = 1;
    if (from >= 0 && to >= 0) {
      nInternalBeads = std::max(0, nBeads);
      firstBeadOffset = 0;
      nBonds = std::max(1, nBeads + 1);
    } else if (from >= 0 || to >= 0) {
      // Last bead is the degree-1 dangling chain end
      nInternalBeads = std::max(0, nBeads - 1);
      firstBeadOffset = 0;
      nBonds = std::max(1, nBeads);
    } else {
      // First and last beads are the degree-1 free chain ends
      nInternalBeads = std::max(0, nBeads - 2);
      firstBeadOffset = 1;
      nBonds = std::max(1, nBeads - 1);
    }

    bboxes[s].seed = strandSeeds[s];
    bboxes[s].nInternalBeads = static_cast<uint32_t>(nInternalBeads);
    bboxes[s].firstBeadOffset = static_cast<uint32_t>(firstBeadOffset);
    bboxes[s].nBonds = static_cast<uint32_t>(nBonds);
    bboxes[s].isPrimaryLoop = (from >= 0 && from == to);
    bboxes[s].isIgnored = strandIsIgnored[s];

    const bool ignoreInOffset =
      (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
       strandIsIgnored[s]) ||
      nInternalBeads <= 0;
    strandInternalOffset[s + 1] =
      strandInternalOffset[s] +
      (ignoreInOffset ? 0 : static_cast<size_t>(nInternalBeads));

    if (ignoreInOffset) {
      continue;
    }

    Eigen::VectorXd coords = sampleStrandBeadsOnDemand(
      box, crosslinkerUniverse, s, strandSeeds[s], nMcSteps);
    if (coords.size() < 3 * nBeads) {
      continue;
    }

    const auto beadBlock =
      coords.segment(3 * firstBeadOffset, 3 * nInternalBeads)
        .reshaped(3, nInternalBeads);
    const Eigen::Vector3d minCoord = beadBlock.rowwise().minCoeff();
    const Eigen::Vector3d maxCoord = beadBlock.rowwise().maxCoeff();

    Eigen::Vector3d unwrappedCenter = 0.5 * (minCoord + maxCoord);
    bboxes[s].halfExtents = 0.5 * (maxCoord - minCoord);
    bboxes[s].center = unwrappedCenter;
    box.handlePBC(bboxes[s].center);
    bboxes[s].radius = bboxes[s].halfExtents.norm();
    cache.putIfFits(s, std::move(coords));
  }

  const size_t totalInternalBeads = strandInternalOffset[nStrands];
  if (totalInternalBeads < 2) {
    return result;
  }

  // Build spatial grid over strand AABBs using median radius for cell sizing
  std::vector<double> activeRadii;
  activeRadii.reserve(nStrands);
  for (size_t s = 0; s < nStrands; ++s) {
    if (bboxes[s].nInternalBeads == 0 ||
        (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
         bboxes[s].isIgnored)) {
      continue;
    }
    activeRadii.push_back(bboxes[s].radius);
  }
  double medianRadius = 0.0;
  if (!activeRadii.empty()) {
    const size_t midIdx = activeRadii.size() / 2;
    std::nth_element(activeRadii.begin(),
                     activeRadii.begin() + static_cast<std::ptrdiff_t>(midIdx),
                     activeRadii.end());
    medianRadius = activeRadii[midIdx];
  }
  StrandSpatialGrid strandGrid(box, 2.0 * medianRadius + upperCutoff, nStrands);
  std::vector<size_t> strandCell(nStrands, 0);
  for (size_t s = 0; s < nStrands; ++s) {
    if (bboxes[s].nInternalBeads == 0 ||
        (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
         bboxes[s].isIgnored)) {
      continue;
    }
    strandGrid.insertAABB(
      s, bboxes[s].center, bboxes[s].halfExtents, upperCutoff);
    strandCell[s] = strandGrid.cellIndexForPoint(bboxes[s].center);
  }

  std::vector<uint8_t> beadIsPaired(totalInternalBeads, 0);
  std::vector<size_t> candidateOrder;
  candidateOrder.reserve(totalInternalBeads);

  struct BatchCandidate
  {
    size_t cellIdx;
    size_t strandIdx;
    size_t beadIdx;
    size_t globalIdx;
  };
  std::vector<BatchCandidate> batchCandidates;

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

  struct PartnerCandidate
  {
    size_t strandIdx;
    size_t beadIdx;
    size_t globalIdx;
    Eigen::Vector3d delta;
  };
  std::vector<PartnerCandidate> candidatePartners;
  candidatePartners.reserve(32);

  size_t numEntanglementsSampled = 0;
  size_t numLinksFoundInIteration = 0;

  do {
    numLinksFoundInIteration = 0;
    // Shuffle all unpaired beads globally so every spatial region is sampled
    // uniformly, then process in batches bounded by remainingNeeded and sorted
    // by spatial cell so nearby strands stay hot in the LRU cache.
    candidateOrder.clear();
    for (size_t g = 0; g < totalInternalBeads; ++g) {
      if (!beadIsPaired[g]) {
        candidateOrder.push_back(g);
      }
    }
    std::ranges::shuffle(candidateOrder, masterRng);

    for (size_t candidatePos = 0; candidatePos < candidateOrder.size();) {
      const size_t currentTargetCount =
        (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT)
          ? result.size()
          : numEntanglementsSampled;
      if (currentTargetCount >= nrOfEntanglementsToSample) {
        break;
      }
      const size_t remainingNeeded =
        nrOfEntanglementsToSample - currentTargetCount;
      const size_t batchEnd =
        std::min(candidateOrder.size(), candidatePos + remainingNeeded);

      batchCandidates.clear();
      batchCandidates.reserve(batchEnd - candidatePos);
      for (size_t idx = candidatePos; idx < batchEnd; ++idx) {
        const size_t g1 = candidateOrder[idx];
        if (beadIsPaired[g1]) {
          continue;
        }
        const std::vector<size_t>::const_iterator it1 =
          std::upper_bound(
            strandInternalOffset.cbegin(), strandInternalOffset.cend(), g1) -
          1;
        const size_t s1 = static_cast<size_t>(
          std::distance(strandInternalOffset.cbegin(), it1));
        const size_t m1 = g1 - *it1;
        batchCandidates.push_back(BatchCandidate{ strandCell[s1], s1, m1, g1 });
      }
      candidatePos = batchEnd;

      std::sort(batchCandidates.begin(),
                batchCandidates.end(),
                [](const BatchCandidate& a, const BatchCandidate& b) {
                  if (a.cellIdx != b.cellIdx) {
                    return a.cellIdx < b.cellIdx;
                  }
                  return a.strandIdx < b.strandIdx;
                });

      for (const BatchCandidate& cand : batchCandidates) {
        const size_t g1 = cand.globalIdx;
        if (beadIsPaired[g1]) {
          continue;
        }
        const size_t s1 = cand.strandIdx;
        const size_t m1 = cand.beadIdx;
        const size_t firstOffset1 =
          static_cast<size_t>(bboxes[s1].firstBeadOffset);
        const size_t k1 = firstOffset1 + m1;
        const size_t baseOffset1 = strandInternalOffset[s1];

        const Eigen::VectorXd& coords1 = cache.get(s1, [&]() {
          return sampleStrandBeadsOnDemand(
            box, crosslinkerUniverse, s1, strandSeeds[s1], nMcSteps);
        });
        const Eigen::Vector3d r1 =
          coords1.segment<3>(static_cast<Eigen::Index>(3 * k1));
        Eigen::Vector3d r1Wrapped = r1;
        box.handlePBC(r1Wrapped);

        candidatePartners.clear();

        // 1. Same-strand candidate beads
        const size_t nInternalBeads1 =
          static_cast<size_t>(bboxes[s1].nInternalBeads);
        for (size_t m2 = 0; m2 < nInternalBeads1; ++m2) {
          if (m1 == m2) {
            continue;
          }
          if (sameStrandCutoff > 0.0) {
            double bondDist =
              std::abs(static_cast<double>(m1) - static_cast<double>(m2));
            if (bboxes[s1].isPrimaryLoop) {
              bondDist = std::min(
                bondDist, static_cast<double>(bboxes[s1].nBonds) - bondDist);
            }
            if (bondDist < sameStrandCutoff) {
              continue;
            }
          }
          const size_t g2 = baseOffset1 + m2;
          if (beadIsPaired[g2]) {
            continue;
          }
          const Eigen::Vector3d r2 = coords1.segment<3>(
            static_cast<Eigen::Index>(3 * (firstOffset1 + m2)));
          Eigen::Vector3d delta = r2 - r1;
          box.handlePBC(delta);
          const double d2 = delta.squaredNorm();
          if (d2 <= upperCutoff * upperCutoff &&
              d2 >= lowerCutoff * lowerCutoff) {
            candidatePartners.push_back({ s1, m2, g2, delta });
          }
        }

        // 2. Different-strand candidate beads
        const std::vector<size_t>& nearbyStrands =
          strandGrid.queryPoint(r1Wrapped);

        for (const size_t s2 : nearbyStrands) {
          if (s2 == s1 ||
              (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
               bboxes[s2].isIgnored) ||
              bboxes[s2].nInternalBeads == 0) {
            continue;
          }
          if (!isAABBWithinCutoff(r1Wrapped,
                                  bboxes[s2].center,
                                  bboxes[s2].halfExtents,
                                  box,
                                  upperCutoff)) {
            continue;
          }

          const Eigen::VectorXd& coords2 = cache.get(s2, [&]() {
            return sampleStrandBeadsOnDemand(
              box, crosslinkerUniverse, s2, strandSeeds[s2], nMcSteps);
          });
          const size_t nInternalBeads2 =
            static_cast<size_t>(bboxes[s2].nInternalBeads);
          const size_t firstOffset2 =
            static_cast<size_t>(bboxes[s2].firstBeadOffset);
          const size_t baseOffset2 = strandInternalOffset[s2];

          for (size_t m2 = 0; m2 < nInternalBeads2; ++m2) {
            const size_t g2 = baseOffset2 + m2;
            if (beadIsPaired[g2]) {
              continue;
            }
            const Eigen::Vector3d r2 = coords2.segment<3>(
              static_cast<Eigen::Index>(3 * (firstOffset2 + m2)));
            Eigen::Vector3d delta = r2 - r1;
            box.handlePBC(delta);
            const double d2 = delta.squaredNorm();
            if (d2 <= upperCutoff * upperCutoff &&
                d2 >= lowerCutoff * lowerCutoff) {
              candidatePartners.push_back({ s2, m2, g2, delta });
            }
          }
        }

        if (candidatePartners.empty()) {
          beadIsPaired[g1] = 1;
          continue;
        }

        const size_t chosenIdx = std::uniform_int_distribution<size_t>{
          0, candidatePartners.size() - 1
        }(masterRng);
        const PartnerCandidate& chosen = candidatePartners[chosenIdx];

        beadIsPaired[g1] = 1;
        beadIsPaired[chosen.globalIdx] = 1;

        if (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT ||
            filteringMode == EntanglementFilteringMode::NONE ||
            (!strandIsIgnored[s1] && !strandIsIgnored[chosen.strandIdx])) {
          const size_t s2 = chosen.strandIdx;
          const size_t m2 = chosen.beadIdx;
          const size_t k2 =
            static_cast<size_t>(bboxes[s2].firstBeadOffset) + m2;
          result.strandIdx1.push_back(s1);
          result.strandIdx2.push_back(s2);
          result.alpha1.push_back(static_cast<double>(m1 + 1) /
                                  static_cast<double>(bboxes[s1].nBonds));
          result.alpha2.push_back(static_cast<double>(m2 + 1) /
                                  static_cast<double>(bboxes[s2].nBonds));
          result.beadIdx1.push_back(k1);
          result.beadIdx2.push_back(k2);

          const Eigen::Vector3d pos1 = r1;
          const Eigen::Vector3d pos2 = r1 + chosen.delta;
          Eigen::Vector3d mid = r1 + 0.5 * chosen.delta;
          box.handlePBC(mid);
          sampledPos1.push_back(pos1);
          sampledPos2.push_back(pos2);
          sampledMidpoints.push_back(mid);
        }

        numEntanglementsSampled += 1;
        numLinksFoundInIteration += 1;

        const size_t targetCount =
          (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT)
            ? result.size()
            : numEntanglementsSampled;
        if (targetCount >= nrOfEntanglementsToSample) {
          break;
        }
      }
    }

    const size_t targetCount =
      (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT)
        ? result.size()
        : numEntanglementsSampled;
    if (targetCount >= nrOfEntanglementsToSample) {
      break;
    }
  } while (((filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT)
              ? result.size() < minimumNrOfEntanglements
              : numEntanglementsSampled < minimumNrOfEntanglements) &&
           numLinksFoundInIteration > 0);

  const size_t nRecorded = result.strandIdx1.size();
  const Eigen::Index nRecordedIdx = static_cast<Eigen::Index>(nRecorded);
  result.x.resize(nRecordedIdx);
  result.y.resize(nRecordedIdx);
  result.z.resize(nRecordedIdx);
  result.position1.resize(3 * nRecordedIdx);
  result.position2.resize(3 * nRecordedIdx);
  for (size_t e = 0; e < nRecorded; ++e) {
    const Eigen::Index eIdx = static_cast<Eigen::Index>(e);
    result.x[eIdx] = sampledMidpoints[e].x();
    result.y[eIdx] = sampledMidpoints[e].y();
    result.z[eIdx] = sampledMidpoints[e].z();
    result.position1.segment<3>(3 * eIdx) = sampledPos1[e];
    result.position2.segment<3>(3 * eIdx) = sampledPos2[e];
  }

  return result;
}
}
