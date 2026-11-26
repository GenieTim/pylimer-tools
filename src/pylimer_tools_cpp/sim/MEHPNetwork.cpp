#include "MEHPNetwork.h"

#include "../entities/Atom.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <stack>

namespace pylimer_tools::sim::mehp {

namespace {

  /**
   * @brief Synchronize box geometry fields (`box`, `L`, `boxHalfs`, `vol`).
   */
  void populateBoxMetadata(MEHPNetwork& net,
                           const pylimer_tools::entities::Box& box)
  {
    net.box = box;
    net.L[0] = net.box.getLx();
    net.L[1] = net.box.getLy();
    net.L[2] = net.box.getLz();
    for (size_t dir = 0; dir < 3; ++dir) {
      net.boxHalfs[dir] = 0.5 * net.L[dir];
    }
    net.vol = net.box.getVolume();
  }

  /**
   * @brief Set `springCoordinateIndexA` and `springCoordinateIndexB` for one
   * spring.
   */
  void setSpringCoordinateIndicesAt(MEHPNetwork& net,
                                    const Eigen::Index sIdx,
                                    const int linkA,
                                    const int linkB)
  {
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      net.springCoordinateIndexA[3 * sIdx + dir] =
        3 * linkA + static_cast<int>(dir);
      net.springCoordinateIndexB[3 * sIdx + dir] =
        3 * linkB + static_cast<int>(dir);
    }
  }

  /**
   * @brief Renumber Eigen vector/array entries while preserving disabled (`<
   * 0`) entries.
   */
  template<typename VecType>
  void safeRenumberWithMapping(VecType& v,
                               const std::vector<long int>& newMapping)
  {
    for (Eigen::Index idx = 0; idx < v.size(); ++idx) {
      if (v[idx] < 0) {
        continue;
      }
      const size_t oldIdx = static_cast<size_t>(v[idx]);
      assert(oldIdx < newMapping.size());
      assert(newMapping[oldIdx] >= 0);
      v[idx] = static_cast<typename VecType::Scalar>(newMapping[oldIdx]);
    }
  }

  /**
   * @brief Recompute `strandContourLength[strandIdx]` and all constituent
   * `springPartitions` from the current `springContourLength` values.
   */
  void recomputeStrandPartitionsFromContourLengths(MEHPNetwork& net,
                                                   const size_t strandIdx)
  {
    const Eigen::Index sIdx = static_cast<Eigen::Index>(strandIdx);
    double totalN = 0.0;
    for (const size_t spIdx : net.springIndicesOfStrand[strandIdx]) {
      totalN += net.springContourLength[static_cast<Eigen::Index>(spIdx)];
    }
    net.strandContourLength[sIdx] = totalN;
    for (const size_t spIdx : net.springIndicesOfStrand[strandIdx]) {
      const Eigen::Index spE = static_cast<Eigen::Index>(spIdx);
      net.springPartitions[spE] =
        (totalN > 0.0) ? (net.springContourLength[spE] / totalN) : 1.0;
    }
  }

  /**
   * @brief Remove all polymer primary loop strands (along with attached
   * entanglement coupling springs), preserving existing crosslink node indices.
   */
  void removePrimaryLoopStrandsInternal(MEHPNetwork& net)
  {
    std::vector<size_t> springsToRemove;
    for (size_t strandIdx = 0; strandIdx < net.nrOfStrands; ++strandIdx) {
      if (net.isEntanglementStrand(strandIdx) ||
          !net.isLoopingStrand(strandIdx)) {
        continue;
      }
      for (const size_t spIdx : net.springIndicesOfStrand[strandIdx]) {
        springsToRemove.push_back(spIdx);
      }
      pylimer_tools::utils::append(
        springsToRemove, net.getEntanglementSpringsAtStrand(strandIdx));
    }
    if (!springsToRemove.empty()) {
      net.removeSprings(springsToRemove, false, nullptr);
    }
  }

  struct StrandEndpointSummary
  {
    size_t nCrosslinks = 0;
    size_t nDanglingEnds = 0;
    size_t nFreeEnds = 0;
    std::vector<size_t> includedStrandIndices;
    std::vector<size_t> danglingEndPartners;
    std::vector<int> danglingEndTypes;
    std::vector<int> freeEndTypes;
  };

  int getStrandBeadType(const CrosslinkerUniverse& cu,
                        const size_t rawStrandIdx)
  {
    if (rawStrandIdx < cu.strandBeadType.size()) {
      return cu.strandBeadType[rawStrandIdx];
    }
    return 1;
  }

  StrandEndpointSummary collectIncludedStrands(const CrosslinkerUniverse& cu,
                                               const bool includeDanglingChains,
                                               const bool includeFreeChains)
  {
    StrandEndpointSummary summary;
    summary.nCrosslinks = cu.xlinkTypes.size();
    const size_t nRawStrands = cu.strandFrom.size();

    INVALIDARG_EXP_IFN(cu.xlinkX.size() == summary.nCrosslinks &&
                         cu.xlinkY.size() == summary.nCrosslinks &&
                         cu.xlinkZ.size() == summary.nCrosslinks,
                       "Crosslinker coordinate arrays in CrosslinkerUniverse "
                       "have mismatched sizes.");
    INVALIDARG_EXP_IFN(
      cu.strandTo.size() == nRawStrands &&
        cu.beadsInStrand.size() == nRawStrands,
      "Strand arrays in CrosslinkerUniverse have mismatched sizes.");

    summary.includedStrandIndices.reserve(nRawStrands);
    for (size_t s = 0; s < nRawStrands; ++s) {
      const long int fromXl = cu.strandFrom[s];
      const long int toXl = cu.strandTo[s];

      if (fromXl >= 0 && toXl >= 0) {
        INVALIDARG_EXP_IFN(static_cast<size_t>(fromXl) < summary.nCrosslinks &&
                             static_cast<size_t>(toXl) < summary.nCrosslinks,
                           "CrosslinkerUniverse strand references out-of-range "
                           "crosslink index.");
        summary.includedStrandIndices.push_back(s);
      } else if (fromXl >= 0 || toXl >= 0) {
        const long int partner = (fromXl >= 0) ? fromXl : toXl;
        INVALIDARG_EXP_IFN(static_cast<size_t>(partner) < summary.nCrosslinks,
                           "CrosslinkerUniverse dangling strand references "
                           "out-of-range crosslink.");
        if (includeDanglingChains) {
          summary.includedStrandIndices.push_back(s);
          summary.danglingEndPartners.push_back(static_cast<size_t>(partner));
          summary.danglingEndTypes.push_back(getStrandBeadType(cu, s));
          summary.nDanglingEnds += 1;
        }
      } else {
        if (includeFreeChains) {
          summary.includedStrandIndices.push_back(s);
          const int bType = getStrandBeadType(cu, s);
          summary.freeEndTypes.push_back(bType);
          summary.freeEndTypes.push_back(bType);
          summary.nFreeEnds += 2;
        }
      }
    }
    return summary;
  }

  /**
   * @brief Return candidate link displacement or zero if displacement vector is
   * omitted.
   */
  Eigen::Vector3d displacementAt(const Eigen::VectorXd* u, const size_t linkIdx)
  {
    if (u == nullptr || u->size() == 0) {
      return Eigen::Vector3d::Zero();
    }
    return u->segment<3>(3 * static_cast<Eigen::Index>(linkIdx));
  }

  /**
   * @brief Split an existing sub-spring cleanly into `existingSpringIdx`
   * `[oldA -> insertedLinkIdx]` and `newSpringIdx` `[insertedLinkIdx -> oldB]`.
   *
   * Updates connectivity arrays, `linkIndicesOfStrand`,
   * `springIndicesOfStrand`, `strandIndicesOfLink`, `springPartitions`,
   * `springContourLength`, and box offsets.
   */
  void insertLinkOntoSubSpring(MEHPNetwork& net,
                               const Eigen::VectorXd& u,
                               const size_t strandIdx,
                               const size_t subSpringPos,
                               const size_t existingSpringIdx,
                               const size_t newSpringIdx,
                               const size_t insertedLinkIdx,
                               const double firstSubSpringPartition,
                               const bool is2D)
  {
    const Eigen::Index dividedE = static_cast<Eigen::Index>(existingSpringIdx);
    const Eigen::Index addedE = static_cast<Eigen::Index>(newSpringIdx);
    const Eigen::Index strandE = static_cast<Eigen::Index>(strandIdx);

    const int partnerAInt = net.springIndexA[dividedE];
    const int partnerBInt = net.springIndexB[dividedE];
    const double oldPartition = net.springPartitions[dividedE];
    const double secondPartition = oldPartition - firstSubSpringPartition;
    RUNTIME_EXP_IFN(firstSubSpringPartition > 0.0 && secondPartition > 0.0,
                    "Split spring partitions must be strictly positive.");

    const Eigen::Vector3d distanceBefore =
      net.evaluateSpringVector(u, existingSpringIdx, is2D);

    net.springIndexB[dividedE] = static_cast<int>(insertedLinkIdx);
    setSpringCoordinateIndicesAt(
      net, dividedE, partnerAInt, static_cast<int>(insertedLinkIdx));

    net.springIndexA[addedE] = static_cast<int>(insertedLinkIdx);
    net.springIndexB[addedE] = partnerBInt;
    setSpringCoordinateIndicesAt(
      net, addedE, static_cast<int>(insertedLinkIdx), partnerBInt);

    net.strandIndexOfSpring[addedE] = static_cast<int>(strandIdx);
    net.springIsEntanglement[addedE] = false;
    net.springBoxOffset.segment<3>(3 * addedE) = Eigen::Vector3d::Zero();

    net.springIndicesOfStrand[strandIdx].insert(
      net.springIndicesOfStrand[strandIdx].begin() +
        static_cast<std::ptrdiff_t>(subSpringPos) + 1,
      newSpringIdx);
    net.linkIndicesOfStrand[strandIdx].insert(
      net.linkIndicesOfStrand[strandIdx].begin() +
        static_cast<std::ptrdiff_t>(subSpringPos) + 1,
      insertedLinkIdx);
    pylimer_tools::utils::addIfNotContained(
      net.strandIndicesOfLink[insertedLinkIdx], strandIdx);

    net.realignSlipLinkToImages(
      u, insertedLinkIdx, existingSpringIdx, newSpringIdx, is2D);

    net.springPartitions[dividedE] = firstSubSpringPartition;
    net.springPartitions[addedE] = secondPartition;
    const double strandN = net.strandContourLength[strandE];
    net.springContourLength[dividedE] = firstSubSpringPartition * strandN;
    net.springContourLength[addedE] = secondPartition * strandN;

    const Eigen::Vector3d distanceAfter =
      net.evaluateSpringVector(u, existingSpringIdx, is2D) +
      net.evaluateSpringVector(u, newSpringIdx, is2D);
    if (net.assumeBoxLargeEnough) {
      Eigen::Vector3d diff = distanceAfter - distanceBefore;
      net.box.handlePBC(diff);
      if (is2D) {
        diff[2] = 0.0;
      }
      RUNTIME_EXP_IFN(
        diff.norm() < 1e-7,
        "Splitting a spring must preserve total end-to-end vector.");
    } else {
      RUNTIME_EXP_IFN(
        pylimer_tools::utils::vector_approx_equal(
          distanceBefore, distanceAfter, 1e-7),
        "Splitting a spring must preserve total end-to-end vector.");
    }
  }

  /**
   * @brief Merge `removedSpringIdx` into `keptSpringIdx` around `linkToReduce`
   * without shrinking array capacities (`removedSpringIdx` is detached with
   * `-1` indices ready for `removeSprings`).
   */
  void mergeTwoSpringsWithoutRemoval(MEHPNetwork& net,
                                     const Eigen::VectorXd& u,
                                     const size_t removedSpringIdx,
                                     const size_t keptSpringIdx,
                                     const size_t linkToReduce)
  {
    const Eigen::Index keptIdx = static_cast<Eigen::Index>(keptSpringIdx);
    const Eigen::Index removedIdx = static_cast<Eigen::Index>(removedSpringIdx);
    const int linkToReduceInt = static_cast<int>(linkToReduce);

    INVALIDARG_EXP_IFN(removedSpringIdx < net.nrOfSprings,
                       "Removed spring index out of bounds.");
    INVALIDARG_EXP_IFN(keptSpringIdx < net.nrOfSprings,
                       "Kept spring index out of bounds.");
    INVALIDARG_EXP_IFN(keptSpringIdx != removedSpringIdx,
                       "Cannot merge a spring with itself.");
    INVALIDARG_EXP_IFN(!net.springIsEntanglement[keptIdx],
                       "Cannot merge entanglement coupling spring.");
    INVALIDARG_EXP_IFN(!net.springIsEntanglement[removedIdx],
                       "Cannot merge entanglement coupling spring.");
    INVALIDARG_EXP_IFN(
      XOR(net.springIndexA[keptIdx] == linkToReduceInt,
          net.springIndexB[keptIdx] == linkToReduceInt),
      "Link to reduce must appear on exactly one end of kept spring.");
    INVALIDARG_EXP_IFN(
      XOR(net.springIndexA[removedIdx] == linkToReduceInt,
          net.springIndexB[removedIdx] == linkToReduceInt),
      "Link to reduce must appear on exactly one end of removed spring.");

    const bool isStrandMerger =
      net.strandIndexOfSpring[keptIdx] != net.strandIndexOfSpring[removedIdx];
    const Eigen::Vector3d totalDistance =
      net.evaluateSpringVectorTo(u, keptSpringIdx, linkToReduce) +
      net.evaluateSpringVectorFrom(u, removedSpringIdx, linkToReduce);

    const bool keptIsA = net.springIndexA[keptIdx] == linkToReduceInt;
    const size_t otherEndKept = static_cast<size_t>(
      keptIsA ? net.springIndexB[keptIdx] : net.springIndexA[keptIdx]);
    const bool removedIsA = net.springIndexA[removedIdx] == linkToReduceInt;
    const size_t otherEndRemoved = static_cast<size_t>(
      removedIsA ? net.springIndexB[removedIdx] : net.springIndexA[removedIdx]);

    const size_t removedStrandIdx =
      static_cast<size_t>(net.strandIndexOfSpring[removedIdx]);
    const size_t keptStrandIdx =
      static_cast<size_t>(net.strandIndexOfSpring[keptIdx]);

    net.springContourLength[keptIdx] += net.springContourLength[removedIdx];
    net.springIndexA[keptIdx] =
      static_cast<int>(keptIsA ? otherEndRemoved : otherEndKept);
    net.springIndexB[keptIdx] =
      static_cast<int>(keptIsA ? otherEndKept : otherEndRemoved);
    net.springBoxOffset.segment<3>(3 * keptIdx) +=
      net.springBoxOffset.segment<3>(3 * removedIdx) *
      (keptIsA == removedIsA ? -1.0 : 1.0);
    setSpringCoordinateIndicesAt(
      net, keptIdx, net.springIndexA[keptIdx], net.springIndexB[keptIdx]);

    if (isStrandMerger) {
      std::vector<size_t> removedStrandsLinks =
        net.linkIndicesOfStrand[removedStrandIdx];
      std::vector<size_t> removedStrandsSprings =
        net.springIndicesOfStrand[removedStrandIdx];

      for (const size_t springIdx : removedStrandsSprings) {
        net.strandIndexOfSpring[static_cast<Eigen::Index>(springIdx)] =
          static_cast<int>(keptStrandIdx);
      }

      std::vector<size_t> uniqueRemovedLinks = removedStrandsLinks;
      pylimer_tools::utils::sort_remove_duplicates(uniqueRemovedLinks);
      for (const size_t lIdx : uniqueRemovedLinks) {
        for (size_t& sIdx : net.strandIndicesOfLink[lIdx]) {
          if (sIdx == removedStrandIdx) {
            sIdx = keptStrandIdx;
            break;
          }
        }
        pylimer_tools::utils::sort_remove_duplicates(
          net.strandIndicesOfLink[lIdx]);
      }

      bool invertSpringDirections = false;
      if (!keptIsA && removedIsA) {
        net.linkIndicesOfStrand[keptStrandIdx].pop_back();
        pylimer_tools::utils::append(
          net.linkIndicesOfStrand[keptStrandIdx], removedStrandsLinks, 1);
        pylimer_tools::utils::append(
          net.springIndicesOfStrand[keptStrandIdx], removedStrandsSprings, 1);
      } else if (keptIsA && !removedIsA) {
        removedStrandsLinks.pop_back();
        net.linkIndicesOfStrand[keptStrandIdx][0] = removedStrandsLinks.back();
        pylimer_tools::utils::prepend(
          net.linkIndicesOfStrand[keptStrandIdx], removedStrandsLinks, 0, 1);
        removedStrandsSprings.pop_back();
        pylimer_tools::utils::prepend(net.springIndicesOfStrand[keptStrandIdx],
                                      removedStrandsSprings);
      } else if (keptIsA && removedIsA) {
        net.linkIndicesOfStrand[keptStrandIdx][0] = removedStrandsLinks[1];
        pylimer_tools::utils::prepend_inverse(
          net.linkIndicesOfStrand[keptStrandIdx], removedStrandsLinks, 0, 2);
        pylimer_tools::utils::prepend_inverse(
          net.springIndicesOfStrand[keptStrandIdx],
          removedStrandsSprings,
          0,
          1);
        invertSpringDirections = true;
      } else {
        removedStrandsLinks.pop_back();
        net.linkIndicesOfStrand[keptStrandIdx].pop_back();
        removedStrandsSprings.pop_back();
        pylimer_tools::utils::append_inverse(
          net.linkIndicesOfStrand[keptStrandIdx], removedStrandsLinks);
        pylimer_tools::utils::append_inverse(
          net.springIndicesOfStrand[keptStrandIdx], removedStrandsSprings);
        invertSpringDirections = true;
      }

      if (invertSpringDirections) {
        for (const size_t springIdx : removedStrandsSprings) {
          if (springIdx == removedSpringIdx) {
            continue;
          }
          const Eigen::Index sIdx = static_cast<Eigen::Index>(springIdx);
          std::swap(net.springIndexA[sIdx], net.springIndexB[sIdx]);
          for (Eigen::Index dir = 0; dir < 3; ++dir) {
            std::swap(net.springCoordinateIndexA[3 * sIdx + dir],
                      net.springCoordinateIndexB[3 * sIdx + dir]);
          }
          net.springBoxOffset.segment<3>(3 * sIdx) *= -1.0;
        }
      }

      net.linkIndicesOfStrand[removedStrandIdx].clear();
      net.springIndicesOfStrand[removedStrandIdx].clear();
      net.strandContourLength[static_cast<Eigen::Index>(removedStrandIdx)] =
        0.0;
    }

    // Detach removed spring
    net.springIndexA[removedIdx] = -1;
    net.springIndexB[removedIdx] = -1;
    net.strandIndexOfSpring[removedIdx] = -1;

    if (!isStrandMerger) {
      pylimer_tools::utils::removeIfContained(
        net.springIndicesOfStrand[keptStrandIdx], removedSpringIdx);
      pylimer_tools::utils::removeIfContained(
        net.linkIndicesOfStrand[keptStrandIdx], linkToReduce);
    }
    net.strandIndicesOfLink[linkToReduce].clear();

    recomputeStrandPartitionsFromContourLengths(net, keptStrandIdx);

    const Eigen::Vector3d distanceAfter =
      (keptIsA ? -1.0 : 1.0) * net.evaluateSpringVector(u, keptSpringIdx);
    const Eigen::Vector3d linkPos =
      net.coordinates.segment<3>(3 * static_cast<Eigen::Index>(linkToReduce)) +
      u.segment<3>(3 * static_cast<Eigen::Index>(linkToReduce));
    const double mergeTol = 1e-7 + 1e-12 * linkPos.norm();
    auto formatDiag = [&]() {
      return "totalDistance=(" + std::to_string(totalDistance[0]) + "," +
             std::to_string(totalDistance[1]) + "," +
             std::to_string(totalDistance[2]) + "), distanceAfter=(" +
             std::to_string(distanceAfter[0]) + "," +
             std::to_string(distanceAfter[1]) + "," +
             std::to_string(distanceAfter[2]) +
             "), kept=" + std::to_string(keptSpringIdx) +
             ", removed=" + std::to_string(removedSpringIdx) +
             ", link=" + std::to_string(linkToReduce) +
             ", strandMerger=" + std::to_string(isStrandMerger) +
             ", keptIsA=" + std::to_string(keptIsA) +
             ", removedIsA=" + std::to_string(removedIsA);
    };
    if (net.assumeBoxLargeEnough) {
      Eigen::Vector3d diff = distanceAfter - totalDistance;
      net.box.handlePBC(diff);
      RUNTIME_EXP_IFN(
        diff.norm() < mergeTol,
        "Merged spring distance mismatch after merging springs: " +
          formatDiag());
    } else {
      RUNTIME_EXP_IFN(
        pylimer_tools::utils::vector_approx_equal(
          totalDistance, distanceAfter, mergeTol),
        "Merged spring distance mismatch after merging springs: " +
          formatDiag());
    }
  }

  /**
   * @brief Collect unique adjacent link neighbors and attached springs for
   * local reversible slip-link swap energy checks.
   */
  std::vector<size_t> collectLocalNeighborsForLinks(
    const MEHPNetwork& net,
    const std::vector<size_t>& centerLinks)
  {
    std::vector<size_t> neighbors = centerLinks;
    for (const size_t linkIdx : centerLinks) {
      const std::vector<size_t> linkNeighbors =
        net.getNeighbourLinkIndices(linkIdx);
      neighbors.insert(
        neighbors.end(), linkNeighbors.begin(), linkNeighbors.end());
    }
    pylimer_tools::utils::sort_remove_duplicates(neighbors);
    return neighbors;
  }

} // namespace

// =============================================================================
// Network Builders
// =============================================================================

MEHPNetwork
MEHPNetwork::fromUniverse(
  const pylimer_tools::entities::Universe& universe,
  const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements&
    entanglements,
  const int crossLinkerType,
  const EntanglementMode entanglementMode,
  const bool keep2FunctionalCrosslinkers,
  const bool removePrimaryLoops,
  const bool assumeBoxLargeEnough)
{
  MEHPNetwork net;
  populateBoxMetadata(net, universe.getBox());
  net.assumeBoxLargeEnough = assumeBoxLargeEnough;
  net.removePrimaryLoops = removePrimaryLoops;
  net.assumeComplete = false;

  const size_t nAtoms = universe.getNrOfAtoms();
  if (nAtoms == 0) {
    net.extractUniverseMetadata(universe);
    net.validate();
    return net;
  }

  std::vector<std::pair<size_t, size_t>> pairsOfAtoms;
  std::vector<long int> pairOfAtom =
    pylimer_tools::utils::initializeWithValue<long int>(nAtoms, -1);

  if (entanglementMode != EntanglementMode::NONE) {
    pairsOfAtoms = entanglements.pairsOfAtoms;
    if (entanglements.pairOfAtom.size() == nAtoms) {
      pairOfAtom = entanglements.pairOfAtom;
    } else if (entanglements.pairOfAtom.empty() && !pairsOfAtoms.empty()) {
      for (size_t pairIdx = 0; pairIdx < pairsOfAtoms.size(); ++pairIdx) {
        const size_t fstV = static_cast<size_t>(universe.getIdxByAtomId(
          static_cast<long int>(pairsOfAtoms[pairIdx].first)));
        const size_t sndV = static_cast<size_t>(universe.getIdxByAtomId(
          static_cast<long int>(pairsOfAtoms[pairIdx].second)));
        pairOfAtom[fstV] = static_cast<long int>(pairIdx);
        pairOfAtom[sndV] = static_cast<long int>(pairIdx);
      }
    } else {
      INVALIDARG_EXP_IFN(
        entanglements.pairOfAtom.size() == nAtoms ||
          (entanglements.pairOfAtom.empty() && pairsOfAtoms.empty()),
        "Entanglement pairOfAtom length must match universe atom count.");
    }
  }

  const igraph_t& graph = universe.getGraph();
  const std::vector<int> vertexDegrees = universe.getVertexDegrees();
  const std::vector<int> atomTypes = universe.getPropertyValues<int>("type");
  const std::vector<long int> atomIds =
    universe.getPropertyValues<long int>("id");
  assert(vertexDegrees.size() == nAtoms);
  assert(atomTypes.size() == nAtoms);
  assert(atomIds.size() == nAtoms);

  // Pre-resolve atom pair IDs -> vertex indices once before DFS
  std::vector<std::pair<size_t, size_t>> pairVertices;
  pairVertices.reserve(pairsOfAtoms.size());
  for (size_t pairIdx = 0; pairIdx < pairsOfAtoms.size(); ++pairIdx) {
    INVALIDARG_EXP_IFN(
      pairsOfAtoms[pairIdx].first != pairsOfAtoms[pairIdx].second,
      "Entanglement atom pairs cannot pair an atom with itself.");
    const size_t fstV = static_cast<size_t>(universe.getIdxByAtomId(
      static_cast<long int>(pairsOfAtoms[pairIdx].first)));
    const size_t sndV = static_cast<size_t>(universe.getIdxByAtomId(
      static_cast<long int>(pairsOfAtoms[pairIdx].second)));
    pairVertices.emplace_back(fstV, sndV);
  }

  for (size_t i = 0; i < pairOfAtom.size(); ++i) {
    if (pairOfAtom[i] >= 0) {
      INVALIDARG_EXP_IFN(
        vertexDegrees[i] == 2,
        "Only degree-2 strand beads are allowed as entanglement targets.");
      INVALIDARG_EXP_IFN(
        !(keep2FunctionalCrosslinkers && atomTypes[i] == crossLinkerType),
        "A 2-functional crosslinker cannot simultaneously serve as an "
        "entanglement bead when keep2FunctionalCrosslinkers is enabled.");
    }
  }

  // Build starting vertex traversal order:
  // 1) degree > 2 junctions
  // 2) optional degree == 2 crosslinkers when keep2FunctionalCrosslinkers
  // 3) degree == 1 chain ends
  // 4) remaining unentangled vertices (handles isolated cycles cleanly)
  std::vector<igraph_integer_t> startingVertices;
  startingVertices.reserve(nAtoms);
  for (size_t vIdx = 0; vIdx < nAtoms; ++vIdx) {
    if (vertexDegrees[vIdx] > 2) {
      startingVertices.push_back(static_cast<igraph_integer_t>(vIdx));
    }
  }
  if (keep2FunctionalCrosslinkers) {
    for (size_t vIdx = 0; vIdx < nAtoms; ++vIdx) {
      if (vertexDegrees[vIdx] == 2 && atomTypes[vIdx] == crossLinkerType &&
          pairOfAtom[vIdx] < 0) {
        startingVertices.push_back(static_cast<igraph_integer_t>(vIdx));
      }
    }
  }
  for (size_t vIdx = 0; vIdx < nAtoms; ++vIdx) {
    if (vertexDegrees[vIdx] == 1 && pairOfAtom[vIdx] < 0) {
      startingVertices.push_back(static_cast<igraph_integer_t>(vIdx));
    }
  }
  for (size_t vIdx = 0; vIdx < nAtoms; ++vIdx) {
    if (vertexDegrees[vIdx] == 2 && pairOfAtom[vIdx] < 0 &&
        !(keep2FunctionalCrosslinkers && atomTypes[vIdx] == crossLinkerType)) {
      startingVertices.push_back(static_cast<igraph_integer_t>(vIdx));
    }
  }

  const Eigen::VectorXd vertexCoordinates =
    universe.getUnwrappedVertexCoordinates(net.box);
  Eigen::VectorXd effectiveCoordinates = vertexCoordinates;

  std::vector<long int> oldVertexIdToNewLinkId =
    pylimer_tools::utils::initializeWithValue<long int>(nAtoms, -1);
  std::vector<bool> vertexVisited =
    pylimer_tools::utils::initializeWithValue(nAtoms, false);
  std::vector<bool> edgeFollowed =
    pylimer_tools::utils::initializeWithValue(universe.getNrOfBonds(), false);

  igraph_adjlist_t adjacencyList;
  igraph_adjlist_init(
    &graph, &adjacencyList, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
  igraph_inclist_t incidenceList;
  igraph_inclist_init(&graph, &incidenceList, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

  std::vector<size_t> springFromVertexIdx;
  std::vector<size_t> springToVertexIdx;
  std::vector<size_t> springContourLen;
  std::vector<Eigen::Vector3d> springBoxOffsets;
  std::vector<Eigen::Vector3d> newLinkPositions;
  std::vector<size_t> newLinkIdxToOldVertexIdx;
  ArrayXArrayXi springIndicesOfStrand;
  ArrayXArrayXi linkIndicesOfStrand;

  long int nextLinkIdx = 0;
  long int currentStrandIdx = 0;
  long int currentSpringContourLength = 0;

  for (const igraph_integer_t vIdx : startingVertices) {
    if (vertexVisited[static_cast<size_t>(vIdx)]) {
      continue;
    }
    vertexVisited[static_cast<size_t>(vIdx)] = true;
    std::stack<igraph_integer_t> vertexStack;
    vertexStack.push(vIdx);
    std::stack<size_t> newLinkIdStack;
    std::stack<size_t> addedVertexIdStack;
    Eigen::Vector3d currentPosition =
      vertexCoordinates.segment<3>(3 * static_cast<Eigen::Index>(vIdx));

    while (!vertexStack.empty()) {
      assert(addedVertexIdStack.size() == newLinkIdStack.size());
      igraph_integer_t currentVertex = vertexStack.top();
      const size_t currentVertexIdx = static_cast<size_t>(currentVertex);

      Eigen::Vector3d distanceToCurrent =
        effectiveCoordinates.segment<3>(
          3 * static_cast<Eigen::Index>(currentVertex)) -
        currentPosition;
      net.box.handlePBC(distanceToCurrent);
      currentPosition += distanceToCurrent;

      const igraph_vector_int_t* neighbors =
        igraph_adjlist_get(&adjacencyList, currentVertex);
      const igraph_vector_int_t* incidentEdges =
        igraph_inclist_get(&incidenceList, currentVertex);

      if (oldVertexIdToNewLinkId[currentVertexIdx] == -1) {
        bool addNewLink = false;
        bool addSpring = false;
        bool addNewStrand = false;

        if (newLinkIdStack.empty()) {
          addNewLink = true;
        } else if (pairOfAtom[currentVertexIdx] >= 0) {
          const size_t pIdx = static_cast<size_t>(pairOfAtom[currentVertexIdx]);
          const size_t fst = pairVertices[pIdx].first;
          const size_t snd = pairVertices[pIdx].second;
          addNewLink = (entanglementMode == EntanglementMode::SPRING) ||
                       (oldVertexIdToNewLinkId[fst] < 0 &&
                        oldVertexIdToNewLinkId[snd] < 0);
          if (!addNewLink) {
            oldVertexIdToNewLinkId[currentVertexIdx] = std::max(
              oldVertexIdToNewLinkId[fst], oldVertexIdToNewLinkId[snd]);
          }
          // In MERGED_LINKS mode, unify both entangled beads at their midpoint.
          // In SPRING mode, each bead keeps its unwrapped position along its
          // own strand and the compliant entanglement spring bridges the two.
          if (entanglementMode == EntanglementMode::MERGED_LINKS) {
            const size_t partnerVertex = (fst == currentVertexIdx) ? snd : fst;
            Eigen::Vector3d dist =
              effectiveCoordinates.segment<3>(
                3 * static_cast<Eigen::Index>(partnerVertex)) -
              currentPosition;
            net.box.handlePBC(dist);
            currentPosition += dist * (addNewLink ? 0.5 : 1.0);
          }
          effectiveCoordinates.segment<3>(
            3 * static_cast<Eigen::Index>(currentVertex)) = currentPosition;
          addSpring = true;
        } else if (vertexDegrees[currentVertexIdx] > 2 ||
                   vertexDegrees[currentVertexIdx] == 1 ||
                   (keep2FunctionalCrosslinkers &&
                    vertexDegrees[currentVertexIdx] == 2 &&
                    atomTypes[currentVertexIdx] == crossLinkerType)) {
          addNewLink = true;
          addSpring = true;
          addNewStrand = true;
        }

        if (addSpring) {
          const size_t strandIdx = static_cast<size_t>(currentStrandIdx);
          if (springIndicesOfStrand.size() <= strandIdx) {
            springIndicesOfStrand.push_back(
              std::vector<size_t>{ springFromVertexIdx.size() });
            linkIndicesOfStrand.push_back(
              std::vector<size_t>{ newLinkIdStack.top() });
          } else {
            springIndicesOfStrand[strandIdx].push_back(
              springFromVertexIdx.size());
          }
          springFromVertexIdx.push_back(addedVertexIdStack.top());

          if (oldVertexIdToNewLinkId[currentVertexIdx] >= 0) {
            effectiveCoordinates.segment<3>(
              3 * static_cast<Eigen::Index>(currentVertex)) = currentPosition;
          }
          springToVertexIdx.push_back(currentVertexIdx);
          linkIndicesOfStrand[strandIdx].push_back(static_cast<size_t>(
            addNewLink ? nextLinkIdx
                       : oldVertexIdToNewLinkId[currentVertexIdx]));

          const Eigen::Vector3d shouldBeDistance =
            currentPosition -
            effectiveCoordinates.segment<3>(
              3 * static_cast<Eigen::Index>(springFromVertexIdx.back()));
          const Eigen::Vector3d targetLinkPos =
            addNewLink ? currentPosition
                       : newLinkPositions[static_cast<size_t>(
                           oldVertexIdToNewLinkId[currentVertexIdx])];
          const Eigen::Vector3d isDistance =
            targetLinkPos -
            newLinkPositions[static_cast<size_t>(
              oldVertexIdToNewLinkId[springFromVertexIdx.back()])];

          const Eigen::Vector3d offset = shouldBeDistance - isDistance;
          springBoxOffsets.emplace_back(offset);
          springContourLen.push_back(
            static_cast<size_t>(currentSpringContourLength));
          currentSpringContourLength = 0;

          if (!addNewLink) {
            newLinkIdStack.push(static_cast<size_t>(
              oldVertexIdToNewLinkId[springToVertexIdx.back()]));
            addedVertexIdStack.push(currentVertexIdx);
          }
        }
        if (addNewStrand) {
          currentStrandIdx += 1;
        }
        if (addNewLink) {
          oldVertexIdToNewLinkId[currentVertexIdx] = nextLinkIdx;
          newLinkIdxToOldVertexIdx.push_back(currentVertexIdx);
          newLinkPositions.push_back(currentPosition);
          effectiveCoordinates.segment<3>(
            3 * static_cast<Eigen::Index>(currentVertex)) = currentPosition;
          newLinkIdStack.push(static_cast<size_t>(nextLinkIdx));
          addedVertexIdStack.push(currentVertexIdx);
          nextLinkIdx += 1;
        }
      }

      bool goDeeper = false;
      for (igraph_integer_t i = 0; i < igraph_vector_int_size(neighbors); ++i) {
        const igraph_integer_t neighbor = igraph_vector_int_get(neighbors, i);
        const igraph_integer_t edge = igraph_vector_int_get(incidentEdges, i);
        const size_t neighborIdx = static_cast<size_t>(neighbor);
        const size_t edgeIdx = static_cast<size_t>(edge);

        if (vertexVisited[neighborIdx] && edgeFollowed[edgeIdx]) {
          continue;
        }
        if (!vertexVisited[neighborIdx]) {
          goDeeper = true;
          currentSpringContourLength += 1;
          vertexVisited[neighborIdx] = true;
          edgeFollowed[edgeIdx] = true;
          vertexStack.push(neighbor);
          break;
        }

        // Loop closure onto an already-visited junction vertex
        const size_t strandIdx = static_cast<size_t>(currentStrandIdx);
        if (springIndicesOfStrand.size() <= strandIdx) {
          springIndicesOfStrand.push_back(
            std::vector<size_t>{ springFromVertexIdx.size() });
          linkIndicesOfStrand.push_back(
            std::vector<size_t>{ newLinkIdStack.top() });
        } else {
          springIndicesOfStrand[strandIdx].push_back(
            springFromVertexIdx.size());
        }
        linkIndicesOfStrand[strandIdx].push_back(
          static_cast<size_t>(oldVertexIdToNewLinkId[neighborIdx]));
        springFromVertexIdx.push_back(addedVertexIdStack.top());
        springToVertexIdx.push_back(neighborIdx);
        springContourLen.push_back(
          static_cast<size_t>(currentSpringContourLength + 1));
        edgeFollowed[edgeIdx] = true;

        const Eigen::Vector3d targetPos = effectiveCoordinates.segment<3>(
          3 * static_cast<Eigen::Index>(neighbor));
        Eigen::Vector3d lastStepToTarget = currentPosition - targetPos;
        net.box.handlePBC(lastStepToTarget);
        const Eigen::Vector3d actualTargetPos =
          currentPosition - lastStepToTarget;

        const Eigen::Vector3d shouldBeDistance =
          actualTargetPos -
          effectiveCoordinates.segment<3>(
            3 * static_cast<Eigen::Index>(springFromVertexIdx.back()));
        const Eigen::Vector3d isDistance =
          newLinkPositions[static_cast<size_t>(
            oldVertexIdToNewLinkId[neighborIdx])] -
          newLinkPositions[static_cast<size_t>(
            oldVertexIdToNewLinkId[springFromVertexIdx.back()])];

        springBoxOffsets.emplace_back(shouldBeDistance - isDistance);
        currentStrandIdx += 1;
        currentSpringContourLength = 0;
      }

      if (!goDeeper) {
        if (newLinkIdStack.empty()) {
          break;
        }
        if (oldVertexIdToNewLinkId[static_cast<size_t>(currentVertex)] ==
            static_cast<long int>(newLinkIdStack.top())) {
          newLinkIdStack.pop();
          addedVertexIdStack.pop();
          if (newLinkIdStack.empty()) {
            break;
          }
        }
        while (oldVertexIdToNewLinkId[static_cast<size_t>(currentVertex)] !=
               static_cast<long int>(newLinkIdStack.top())) {
          vertexStack.pop();
          currentVertex = vertexStack.top();
        }
        currentPosition = newLinkPositions[newLinkIdStack.top()];
      }
    }
  }

  igraph_inclist_destroy(&incidenceList);
  igraph_adjlist_destroy(&adjacencyList);

  const size_t nPolymerSprings = springFromVertexIdx.size();
  if (entanglementMode == EntanglementMode::SPRING) {
    for (const std::pair<size_t, size_t>& pv : pairVertices) {
      const long int newLink1 = oldVertexIdToNewLinkId[pv.first];
      const long int newLink2 = oldVertexIdToNewLinkId[pv.second];
      RUNTIME_EXP_IFN(
        newLink1 >= 0 && newLink2 >= 0,
        "Entanglement bead was not visited during DFS network construction.");
      springFromVertexIdx.push_back(pv.first);
      springToVertexIdx.push_back(pv.second);
      springContourLen.push_back(1);
      springBoxOffsets.emplace_back(
        net.box.getOffset(newLinkPositions[static_cast<size_t>(newLink2)] -
                          newLinkPositions[static_cast<size_t>(newLink1)]));
      linkIndicesOfStrand.push_back(
        { static_cast<size_t>(newLink1), static_cast<size_t>(newLink2) });
      springIndicesOfStrand.push_back({ springFromVertexIdx.size() - 1 });
      currentStrandIdx += 1;
    }
  }

  net.nrOfSprings = springFromVertexIdx.size();
  net.nrOfStrands = springIndicesOfStrand.size();
  net.nrOfLinks = static_cast<size_t>(nextLinkIdx);
  net.springIndicesOfStrand = std::move(springIndicesOfStrand);
  net.linkIndicesOfStrand = std::move(linkIndicesOfStrand);

  const Eigen::Index nSpringsE = static_cast<Eigen::Index>(net.nrOfSprings);
  const Eigen::Index nStrandsE = static_cast<Eigen::Index>(net.nrOfStrands);
  const Eigen::Index nLinksE = static_cast<Eigen::Index>(net.nrOfLinks);

  net.springIndexA = Eigen::ArrayXi(nSpringsE);
  net.springIndexB = Eigen::ArrayXi(nSpringsE);
  net.springCoordinateIndexA = Eigen::ArrayXi(3 * nSpringsE);
  net.springCoordinateIndexB = Eigen::ArrayXi(3 * nSpringsE);
  net.springBoxOffset = Eigen::VectorXd::Zero(3 * nSpringsE);
  net.springContourLength = Eigen::VectorXd(nSpringsE);
  net.springPartitions = Eigen::VectorXd::Zero(nSpringsE);
  net.strandIndexOfSpring = Eigen::ArrayXi::Constant(nSpringsE, -1);
  net.springIsEntanglement = Eigen::ArrayXb::Constant(nSpringsE, false);
  if (entanglementMode == EntanglementMode::SPRING &&
      net.nrOfSprings > nPolymerSprings) {
    net.springIsEntanglement.tail(
      static_cast<Eigen::Index>(net.nrOfSprings - nPolymerSprings)) = true;
  }

  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    const Eigen::Index sIdx = static_cast<Eigen::Index>(i);
    net.springIndexA[sIdx] =
      static_cast<int>(oldVertexIdToNewLinkId[springFromVertexIdx[i]]);
    net.springIndexB[sIdx] =
      static_cast<int>(oldVertexIdToNewLinkId[springToVertexIdx[i]]);
    net.springContourLength[sIdx] = static_cast<double>(springContourLen[i]);
    setSpringCoordinateIndicesAt(
      net, sIdx, net.springIndexA[sIdx], net.springIndexB[sIdx]);
    net.springBoxOffset.segment<3>(3 * sIdx) = springBoxOffsets[i];
  }

  for (size_t s = 0; s < net.nrOfStrands; ++s) {
    for (const size_t spIdx : net.springIndicesOfStrand[s]) {
      net.strandIndexOfSpring[static_cast<Eigen::Index>(spIdx)] =
        static_cast<int>(s);
    }
  }

  net.strandIndicesOfLink = pylimer_tools::utils::initializeWithValue(
    net.nrOfLinks, std::vector<size_t>());
  for (size_t s = 0; s < net.nrOfStrands; ++s) {
    for (const size_t lIdx : net.linkIndicesOfStrand[s]) {
      pylimer_tools::utils::addIfNotContained(net.strandIndicesOfLink[lIdx], s);
    }
  }

  net.coordinates = Eigen::VectorXd::Zero(3 * nLinksE);
  net.linkIsEntanglement = Eigen::ArrayXb::Constant(nLinksE, false);
  net.oldAtomIds = Eigen::ArrayXi::Constant(nLinksE, -1);
  net.oldAtomTypes = Eigen::ArrayXi::Zero(nLinksE);
  net.nrOfCrosslinkSwapsEndured = Eigen::ArrayXi::Zero(nLinksE);

  for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
    const Eigen::Index lIdx = static_cast<Eigen::Index>(linkIdx);
    const size_t oldVertexIdx = newLinkIdxToOldVertexIdx[linkIdx];
    const bool isEnt = (pairOfAtom[oldVertexIdx] >= 0);
    net.linkIsEntanglement[lIdx] = isEnt;
    net.oldAtomIds[lIdx] = isEnt ? -1 : static_cast<int>(atomIds[oldVertexIdx]);
    net.oldAtomTypes[lIdx] = atomTypes[oldVertexIdx];
    net.coordinates.segment<3>(3 * lIdx) = newLinkPositions[linkIdx];
  }
  net.nrOfNodes =
    static_cast<size_t>((net.linkIsEntanglement == false).count());

  net.strandContourLength = Eigen::VectorXd::Zero(nStrandsE);
  for (size_t s = 0; s < net.nrOfStrands; ++s) {
    recomputeStrandPartitionsFromContourLengths(net, s);
  }
  net.recomputeMeanContourLengths();
  net.extractUniverseMetadata(universe);

  if (removePrimaryLoops) {
    removePrimaryLoopStrandsInternal(net);
  }

  net.validate();
  return net;
}

MEHPNetwork
MEHPNetwork::fromUniverseWithRandomEntanglements(
  const pylimer_tools::entities::Universe& universe,
  const size_t nrOfEntanglementsToSample,
  const double upperSamplingCutoff,
  const double lowerSamplingCutoff,
  const size_t minNrOfEntanglements,
  const double sameStrandCutoff,
  const std::string& seed,
  const int crossLinkerType,
  const EntanglementMode entanglementMode,
  const bool withCrosslinkers,
  const bool filterDanglingAndSoluble,
  const bool keep2FunctionalCrosslinkers,
  const bool removePrimaryLoops,
  const bool assumeBoxLargeEnough)
{
  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements;
  if (entanglementMode != EntanglementMode::NONE &&
      nrOfEntanglementsToSample > 0) {
    entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglementsV2(
        universe,
        nrOfEntanglementsToSample,
        upperSamplingCutoff,
        lowerSamplingCutoff,
        minNrOfEntanglements,
        sameStrandCutoff,
        seed,
        crossLinkerType,
        withCrosslinkers,
        filterDanglingAndSoluble);
  }
  return MEHPNetwork::fromUniverse(universe,
                                   entanglements,
                                   crossLinkerType,
                                   entanglementMode,
                                   keep2FunctionalCrosslinkers,
                                   removePrimaryLoops,
                                   assumeBoxLargeEnough);
}

MEHPNetwork
MEHPNetwork::fromCrosslinkerUniverse(
  const pylimer_tools::entities::Box& box,
  const CrosslinkerUniverse& cu,
  const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
    entanglements,
  const int /*crossLinkerType*/,
  const EntanglementMode entanglementMode,
  const bool removeDanglingChains,
  const bool includeFreeChains,
  const bool removePrimaryLoops,
  const bool assumeBoxLargeEnough)
{
  MEHPNetwork net;
  populateBoxMetadata(net, box);
  net.assumeBoxLargeEnough = assumeBoxLargeEnough;
  net.removePrimaryLoops = removePrimaryLoops;
  net.assumeComplete = false;

  const size_t nRawStrands = cu.strandFrom.size();
  const bool includeDanglingChains = !removeDanglingChains;
  const StrandEndpointSummary summary =
    collectIncludedStrands(cu, includeDanglingChains, includeFreeChains);

  const size_t nCrosslinks = summary.nCrosslinks;
  const size_t nDanglingEnds = summary.nDanglingEnds;
  const size_t nFreeEnds = summary.nFreeEnds;
  const size_t nPolymerStrands = summary.includedStrandIndices.size();

  const bool ignoreEntanglements = (entanglementMode == EntanglementMode::NONE);
  const size_t nEntanglements = ignoreEntanglements ? 0 : entanglements.size();
  const bool entanglementsAsSprings =
    (entanglementMode == EntanglementMode::SPRING);

  if (!ignoreEntanglements) {
    INVALIDARG_EXP_IFN(
      entanglements.strandIdx2.size() == nEntanglements &&
        entanglements.alpha1.size() == nEntanglements &&
        entanglements.alpha2.size() == nEntanglements &&
        static_cast<size_t>(entanglements.x.size()) == nEntanglements &&
        static_cast<size_t>(entanglements.y.size()) == nEntanglements &&
        static_cast<size_t>(entanglements.z.size()) == nEntanglements,
      "Inconsistent StrandPairEntanglements vector sizes.");
  }

  std::vector<long int> rawToIncludedStrand(nRawStrands, -1);
  for (size_t idx = 0; idx < nPolymerStrands; ++idx) {
    rawToIncludedStrand[summary.includedStrandIndices[idx]] =
      static_cast<long int>(idx);
  }

  const size_t nNodes = nCrosslinks + nDanglingEnds + nFreeEnds;
  const size_t nEntLinks =
    entanglementsAsSprings ? (2 * nEntanglements) : nEntanglements;
  net.nrOfNodes = nNodes;
  net.nrOfLinks = nNodes + nEntLinks;

  const Eigen::Index nLinksE = static_cast<Eigen::Index>(net.nrOfLinks);
  net.coordinates = Eigen::VectorXd::Zero(3 * nLinksE);
  net.linkIsEntanglement = Eigen::ArrayXb::Constant(nLinksE, false);
  net.oldAtomIds = Eigen::ArrayXi::Constant(nLinksE, -1);
  net.oldAtomTypes = Eigen::ArrayXi::Ones(nLinksE);
  net.nrOfCrosslinkSwapsEndured = Eigen::ArrayXi::Zero(nLinksE);

  for (size_t i = 0; i < nNodes; ++i) {
    net.oldAtomIds[static_cast<Eigen::Index>(i)] = static_cast<int>(i);
  }

  for (Eigen::Index xl = 0; xl < static_cast<Eigen::Index>(nCrosslinks); ++xl) {
    const size_t xlIdx = static_cast<size_t>(xl);
    Eigen::Vector3d pos(cu.xlinkX[xlIdx], cu.xlinkY[xlIdx], cu.xlinkZ[xlIdx]);
    box.handlePBC(pos);
    net.coordinates.segment<3>(3 * xl) = pos;
    net.oldAtomTypes[xl] = cu.xlinkTypes[xlIdx];
  }

  for (size_t dIdx = 0; dIdx < nDanglingEnds; ++dIdx) {
    const Eigen::Index nodeIdx = static_cast<Eigen::Index>(nCrosslinks + dIdx);
    const Eigen::Index partnerXl =
      static_cast<Eigen::Index>(summary.danglingEndPartners[dIdx]);
    net.coordinates.segment<3>(3 * nodeIdx) =
      net.coordinates.segment<3>(3 * partnerXl);
    net.oldAtomTypes[nodeIdx] = summary.danglingEndTypes[dIdx];
  }

  for (size_t fIdx = 0; fIdx < nFreeEnds; ++fIdx) {
    const Eigen::Index nodeIdx =
      static_cast<Eigen::Index>(nCrosslinks + nDanglingEnds + fIdx);
    net.coordinates.segment<3>(3 * nodeIdx).setZero();
    net.oldAtomTypes[nodeIdx] = summary.freeEndTypes[fIdx];
  }

  struct StrandEntanglementPlacement
  {
    double alpha;
    size_t beadIdx;
    size_t entIdx;
    size_t linkIdx;
  };
  std::vector<std::vector<StrandEntanglementPlacement>> entanglementsOnStrand(
    nPolymerStrands);

  const bool hasSeparatePositions =
    (static_cast<size_t>(entanglements.position1.size()) ==
       3 * nEntanglements &&
     static_cast<size_t>(entanglements.position2.size()) == 3 * nEntanglements);

  for (size_t e = 0; e < nEntanglements; ++e) {
    const size_t rawS1 = entanglements.strandIdx1[e];
    const size_t rawS2 = entanglements.strandIdx2[e];
    INVALIDARG_EXP_IFN(rawS1 < nRawStrands && rawS2 < nRawStrands,
                       "Entanglement strand index out of range.");
    const long int s1 = rawToIncludedStrand[rawS1];
    const long int s2 = rawToIncludedStrand[rawS2];
    INVALIDARG_EXP_IFN(s1 >= 0 && s2 >= 0,
                       "Entanglement references an excluded strand.");

    const double a1 = std::clamp(entanglements.alpha1[e], 0.0, 1.0);
    const double a2 = std::clamp(entanglements.alpha2[e], 0.0, 1.0);
    const size_t k1 =
      (e < entanglements.beadIdx1.size()) ? entanglements.beadIdx1[e] : 0;
    const size_t k2 =
      (e < entanglements.beadIdx2.size()) ? entanglements.beadIdx2[e] : 0;
    const Eigen::Index entE = static_cast<Eigen::Index>(e);

    if (!entanglementsAsSprings) {
      const size_t linkIdx = nNodes + e;
      const Eigen::Index linkE = static_cast<Eigen::Index>(linkIdx);
      Eigen::Vector3d mid(
        entanglements.x[entE], entanglements.y[entE], entanglements.z[entE]);
      box.handlePBC(mid);
      net.coordinates.segment<3>(3 * linkE) = mid;
      net.linkIsEntanglement[linkE] = true;
      net.oldAtomIds[linkE] = -1;
      net.oldAtomTypes[linkE] = getStrandBeadType(cu, rawS1);
      entanglementsOnStrand[static_cast<size_t>(s1)].push_back(
        StrandEntanglementPlacement{ a1, k1, e, linkIdx });
      entanglementsOnStrand[static_cast<size_t>(s2)].push_back(
        StrandEntanglementPlacement{ a2, k2, e, linkIdx });
    } else {
      const size_t link1 = nNodes + 2 * e;
      const size_t link2 = nNodes + 2 * e + 1;
      const Eigen::Index link1E = static_cast<Eigen::Index>(link1);
      const Eigen::Index link2E = static_cast<Eigen::Index>(link2);
      Eigen::Vector3d pos1 =
        hasSeparatePositions
          ? Eigen::Vector3d(entanglements.position1.segment<3>(3 * entE))
          : Eigen::Vector3d(entanglements.x[entE],
                            entanglements.y[entE],
                            entanglements.z[entE]);
      Eigen::Vector3d pos2 =
        hasSeparatePositions
          ? Eigen::Vector3d(entanglements.position2.segment<3>(3 * entE))
          : Eigen::Vector3d(entanglements.x[entE],
                            entanglements.y[entE],
                            entanglements.z[entE]);
      box.handlePBC(pos1);
      box.handlePBC(pos2);
      net.coordinates.segment<3>(3 * link1E) = pos1;
      net.coordinates.segment<3>(3 * link2E) = pos2;
      net.linkIsEntanglement[link1E] = true;
      net.linkIsEntanglement[link2E] = true;
      net.oldAtomIds[link1E] = -1;
      net.oldAtomIds[link2E] = -1;
      net.oldAtomTypes[link1E] = getStrandBeadType(cu, rawS1);
      net.oldAtomTypes[link2E] = getStrandBeadType(cu, rawS2);
      entanglementsOnStrand[static_cast<size_t>(s1)].push_back(
        StrandEntanglementPlacement{ a1, k1, e, link1 });
      entanglementsOnStrand[static_cast<size_t>(s2)].push_back(
        StrandEntanglementPlacement{ a2, k2, e, link2 });
    }
  }

  net.nrOfStrands =
    nPolymerStrands + (entanglementsAsSprings ? nEntanglements : 0);
  net.nrOfSprings =
    nPolymerStrands +
    (entanglementsAsSprings ? (3 * nEntanglements) : (2 * nEntanglements));

  const Eigen::Index nStrandsE = static_cast<Eigen::Index>(net.nrOfStrands);
  const Eigen::Index nSpringsE = static_cast<Eigen::Index>(net.nrOfSprings);

  net.strandContourLength = Eigen::VectorXd::Zero(nStrandsE);
  net.springIndicesOfStrand.resize(net.nrOfStrands);
  net.linkIndicesOfStrand.resize(net.nrOfStrands);
  net.strandIndicesOfLink.assign(net.nrOfLinks, std::vector<size_t>());

  net.springIndexA = Eigen::ArrayXi::Zero(nSpringsE);
  net.springIndexB = Eigen::ArrayXi::Zero(nSpringsE);
  net.springCoordinateIndexA = Eigen::ArrayXi::Zero(3 * nSpringsE);
  net.springCoordinateIndexB = Eigen::ArrayXi::Zero(3 * nSpringsE);
  net.springContourLength = Eigen::VectorXd::Zero(nSpringsE);
  net.springPartitions = Eigen::VectorXd::Zero(nSpringsE);
  net.springBoxOffset = Eigen::VectorXd::Zero(3 * nSpringsE);
  net.springIsEntanglement = Eigen::ArrayXb::Constant(nSpringsE, false);
  net.strandIndexOfSpring = Eigen::ArrayXi::Constant(nSpringsE, -1);

  size_t handledDanglingIdx = 0;
  size_t handledFreeIdx = 0;
  size_t nextSpringIdx = 0;
  double sumBondLength = 0.0;
  double sumR02 = 0.0;
  size_t countBondStats = 0;

  for (size_t s = 0; s < nPolymerStrands; ++s) {
    const size_t rawS = summary.includedStrandIndices[s];
    const long int rawFrom = cu.strandFrom[rawS];
    const long int rawTo = cu.strandTo[rawS];
    const int nBeads = cu.beadsInStrand[rawS];

    size_t nodeFrom = 0;
    size_t nodeTo = 0;
    double totalBonds = 1.0;

    if (rawFrom >= 0 && rawTo >= 0) {
      nodeFrom = static_cast<size_t>(rawFrom);
      nodeTo = static_cast<size_t>(rawTo);
      totalBonds = static_cast<double>(nBeads + 1);
    } else if (rawFrom >= 0 && rawTo < 0) {
      nodeFrom = static_cast<size_t>(rawFrom);
      nodeTo = nCrosslinks + handledDanglingIdx++;
      totalBonds = static_cast<double>(std::max(1, nBeads));
    } else if (rawFrom < 0 && rawTo >= 0) {
      nodeFrom = static_cast<size_t>(rawTo);
      nodeTo = nCrosslinks + handledDanglingIdx++;
      totalBonds = static_cast<double>(std::max(1, nBeads));
    } else {
      nodeFrom = nCrosslinks + nDanglingEnds + handledFreeIdx;
      nodeTo = nodeFrom + 1;
      handledFreeIdx += 2;
      totalBonds = static_cast<double>(std::max(1, nBeads - 1));
    }

    if (rawS < cu.beadDistanceInStrand.size() &&
        cu.beadDistanceInStrand[rawS] > 0.0) {
      sumBondLength += cu.beadDistanceInStrand[rawS];
      countBondStats += 1;
    }
    if (rawS < cu.meanSquaredBeadDistanceInStrand.size() &&
        cu.meanSquaredBeadDistanceInStrand[rawS] > 0.0) {
      sumR02 += cu.meanSquaredBeadDistanceInStrand[rawS];
    }

    std::vector<StrandEntanglementPlacement>& placements =
      entanglementsOnStrand[s];
    std::sort(placements.begin(),
              placements.end(),
              [](const StrandEntanglementPlacement& a,
                 const StrandEntanglementPlacement& b) {
                if (std::abs(a.alpha - b.alpha) > 1e-12) {
                  return a.alpha < b.alpha;
                }
                if (a.beadIdx != b.beadIdx) {
                  return a.beadIdx < b.beadIdx;
                }
                return a.entIdx < b.entIdx;
              });

    const size_t kEnt = placements.size();
    const size_t nSubSprings = kEnt + 1;

    net.linkIndicesOfStrand[s].reserve(nSubSprings + 1);
    net.springIndicesOfStrand[s].reserve(nSubSprings);
    net.linkIndicesOfStrand[s].push_back(nodeFrom);
    for (size_t m = 0; m < kEnt; ++m) {
      net.linkIndicesOfStrand[s].push_back(placements[m].linkIdx);
    }
    net.linkIndicesOfStrand[s].push_back(nodeTo);

    for (const size_t lIdx : net.linkIndicesOfStrand[s]) {
      pylimer_tools::utils::addIfNotContained(net.strandIndicesOfLink[lIdx], s);
    }

    const Eigen::Vector3d pStart =
      net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * nodeFrom));
    Eigen::Vector3d strandChord =
      net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * nodeTo)) -
      pStart;
    box.handlePBC(strandChord);
    const Eigen::Vector3d pEnd = pStart + strandChord;

    std::vector<Eigen::Vector3d> unwrappedPositions(nSubSprings + 1);
    unwrappedPositions[0] = pStart;
    unwrappedPositions[nSubSprings] = pEnd;
    for (size_t m = 0; m < kEnt; ++m) {
      const Eigen::Vector3d chordPos =
        pStart + placements[m].alpha * strandChord;
      Eigen::Vector3d diff =
        net.coordinates.segment<3>(
          static_cast<Eigen::Index>(3 * placements[m].linkIdx)) -
        chordPos;
      box.handlePBC(diff);
      unwrappedPositions[m + 1] = chordPos + diff;
    }

    std::vector<double> subContourLengths(nSubSprings, 0.0);
    double prevBondPos = 0.0;
    for (size_t m = 0; m < kEnt; ++m) {
      double bondPos = placements[m].alpha * totalBonds;
      if (std::abs(bondPos - std::round(bondPos)) < 1e-7) {
        bondPos = std::round(bondPos);
      }
      const double minAllowed = prevBondPos + 1e-6;
      const double maxAllowed =
        totalBonds - static_cast<double>(kEnt - m) * 1e-6;
      bondPos = std::clamp(bondPos, minAllowed, maxAllowed);
      subContourLengths[m] = bondPos - prevBondPos;
      prevBondPos = bondPos;
    }
    subContourLengths[kEnt] = std::max(1e-6, totalBonds - prevBondPos);

    for (size_t m = 0; m < nSubSprings; ++m) {
      const size_t spIdx = nextSpringIdx++;
      const Eigen::Index spE = static_cast<Eigen::Index>(spIdx);
      const size_t linkA = net.linkIndicesOfStrand[s][m];
      const size_t linkB = net.linkIndicesOfStrand[s][m + 1];

      net.springIndicesOfStrand[s].push_back(spIdx);
      net.strandIndexOfSpring[spE] = static_cast<int>(s);
      net.springIndexA[spE] = static_cast<int>(linkA);
      net.springIndexB[spE] = static_cast<int>(linkB);
      setSpringCoordinateIndicesAt(
        net, spE, static_cast<int>(linkA), static_cast<int>(linkB));
      net.springContourLength[spE] = subContourLengths[m];

      const Eigen::Vector3d shouldBeDist =
        unwrappedPositions[m + 1] - unwrappedPositions[m];
      const Eigen::Vector3d isDist =
        net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * linkB)) -
        net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * linkA));
      net.springBoxOffset.segment<3>(3 * spE) = shouldBeDist - isDist;
    }

    recomputeStrandPartitionsFromContourLengths(net, s);
  }

  if (entanglementsAsSprings) {
    for (size_t e = 0; e < nEntanglements; ++e) {
      const size_t entStrandIdx = nPolymerStrands + e;
      const size_t spIdx = nextSpringIdx++;
      const Eigen::Index spE = static_cast<Eigen::Index>(spIdx);
      const size_t link1 = nNodes + 2 * e;
      const size_t link2 = nNodes + 2 * e + 1;

      net.linkIndicesOfStrand[entStrandIdx] = { link1, link2 };
      net.springIndicesOfStrand[entStrandIdx] = { spIdx };
      pylimer_tools::utils::addIfNotContained(net.strandIndicesOfLink[link1],
                                              entStrandIdx);
      pylimer_tools::utils::addIfNotContained(net.strandIndicesOfLink[link2],
                                              entStrandIdx);

      net.strandIndexOfSpring[spE] = static_cast<int>(entStrandIdx);
      net.springIsEntanglement[spE] = true;
      net.springContourLength[spE] = 1.0;
      net.springPartitions[spE] = 1.0;
      net.strandContourLength[static_cast<Eigen::Index>(entStrandIdx)] = 1.0;
      net.springIndexA[spE] = static_cast<int>(link1);
      net.springIndexB[spE] = static_cast<int>(link2);
      setSpringCoordinateIndicesAt(
        net, spE, static_cast<int>(link1), static_cast<int>(link2));
      net.springBoxOffset.segment<3>(3 * spE) = box.getOffset(
        net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * link2)) -
        net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * link1)));
    }
  }

  assert(nextSpringIdx == net.nrOfSprings);
  if (countBondStats > 0) {
    net.defaultBondLength = sumBondLength / static_cast<double>(countBondStats);
    net.defaultR02 = (sumR02 > 0.0)
                       ? (sumR02 / static_cast<double>(countBondStats))
                       : (net.defaultBondLength * net.defaultBondLength);
  }

  net.extractSelfMetadata();
  if (removePrimaryLoops) {
    removePrimaryLoopStrandsInternal(net);
  }
  net.recomputeMeanContourLengths();
  net.validate();
  return net;
}

// =============================================================================
// Metadata Extraction & Averages
// =============================================================================

void
MEHPNetwork::extractUniverseMetadata(const pylimer_tools::entities::Universe& u)
{
  this->universeNrOfAtoms = u.getNrOfAtoms();
  this->defaultBondLength =
    u.getNrOfBonds() > 0 ? u.computeMeanBondLength() : 0.0;
  this->defaultR02 =
    u.getNrOfBonds() > 0 ? u.computeMeanSquaredBondLength() : 1.0;
  this->clusterSizes.clear();
  this->atomIdToClusterIdx.clear();

  if (this->universeNrOfAtoms == 0) {
    return;
  }

  igraph_vector_int_t membership;
  igraph_vector_int_t csize;
  igraph_integer_t noClusters = 0;
  igraph_vector_int_init(&membership, 0);
  igraph_vector_int_init(&csize, 0);

  if (igraph_connected_components(
        &u.getGraph(), &membership, &csize, &noClusters, IGRAPH_WEAK) ==
      IGRAPH_SUCCESS) {
    this->clusterSizes.resize(static_cast<size_t>(noClusters));
    for (igraph_integer_t i = 0; i < noClusters; ++i) {
      this->clusterSizes[static_cast<size_t>(i)] =
        static_cast<size_t>(igraph_vector_int_get(&csize, i));
    }
    for (size_t linkIdx = 0; linkIdx < this->nrOfLinks; ++linkIdx) {
      const long int atomId =
        this->oldAtomIds[static_cast<Eigen::Index>(linkIdx)];
      if (atomId >= 0) {
        if (const std::optional<igraph_integer_t> v =
              u.tryGetIdxByAtomId(atomId);
            v.has_value()) {
          this->atomIdToClusterIdx[atomId] =
            static_cast<size_t>(igraph_vector_int_get(&membership, v.value()));
        }
      }
    }
  }

  igraph_vector_int_destroy(&membership);
  igraph_vector_int_destroy(&csize);
}

void
MEHPNetwork::extractSelfMetadata()
{
  const double inferredAtoms = std::round(this->inferNrOfAtomsFromNetwork());
  this->universeNrOfAtoms = static_cast<size_t>(std::max(0.0, inferredAtoms));
  this->clusterSizes.clear();
  this->atomIdToClusterIdx.clear();

  if (this->universeNrOfAtoms == 0 || this->nrOfLinks == 0) {
    return;
  }

  std::vector<bool> linkVisited(this->nrOfLinks, false);
  std::vector<bool> strandVisited(this->nrOfStrands, false);
  std::vector<size_t> queue;
  queue.reserve(this->nrOfLinks);

  for (size_t startLink = 0; startLink < this->nrOfLinks; ++startLink) {
    if (linkVisited[startLink] ||
        this->linkIsEntanglement[static_cast<Eigen::Index>(startLink)]) {
      continue;
    }

    const size_t clusterIdx = this->clusterSizes.size();
    double clusterAtomCount = 0.0;
    queue.clear();
    linkVisited[startLink] = true;
    queue.push_back(startLink);

    size_t qHead = 0;
    while (qHead < queue.size()) {
      const size_t currLink = queue[qHead++];
      const Eigen::Index currLinkE = static_cast<Eigen::Index>(currLink);
      if (!this->linkIsEntanglement[currLinkE]) {
        clusterAtomCount += 1.0;
        if (currLink < static_cast<size_t>(this->oldAtomIds.size()) &&
            this->oldAtomIds[currLinkE] >= 0) {
          this->atomIdToClusterIdx[this->oldAtomIds[currLinkE]] = clusterIdx;
        }
      }

      for (const size_t sIdx : this->strandIndicesOfLink[currLink]) {
        if (sIdx >= this->nrOfStrands || strandVisited[sIdx]) {
          continue;
        }
        const std::vector<size_t>& springs = this->springIndicesOfStrand[sIdx];
        if (springs.empty() || this->isEntanglementStrand(sIdx)) {
          continue;
        }
        strandVisited[sIdx] = true;
        double strandBonds = 0.0;
        for (const size_t spIdx : springs) {
          strandBonds +=
            this->springContourLength[static_cast<Eigen::Index>(spIdx)];
        }
        clusterAtomCount += (strandBonds - 1.0);

        const std::vector<size_t>& strandLinks =
          this->linkIndicesOfStrand[sIdx];
        if (!strandLinks.empty()) {
          const size_t endA = strandLinks.front();
          const size_t endB = strandLinks.back();
          if (endA < this->nrOfLinks && !linkVisited[endA] &&
              !this->linkIsEntanglement[static_cast<Eigen::Index>(endA)]) {
            linkVisited[endA] = true;
            queue.push_back(endA);
          }
          if (endB < this->nrOfLinks && !linkVisited[endB] &&
              !this->linkIsEntanglement[static_cast<Eigen::Index>(endB)]) {
            linkVisited[endB] = true;
            queue.push_back(endB);
          }
        }
      }
    }

    this->clusterSizes.push_back(
      static_cast<size_t>(std::max(0.0, std::round(clusterAtomCount))));
  }
}

double
MEHPNetwork::inferNrOfAtomsFromNetwork() const
{
  double nAtoms = 0.0;
  for (Eigen::Index i = 0; i < this->springContourLength.size(); ++i) {
    if (!this->springIsEntanglement[i]) {
      nAtoms += this->springContourLength[i] - 1.0;
    }
  }
  nAtoms += static_cast<double>(this->nrOfLinks);
  if (this->springIsEntanglement.size() == 0 ||
      !this->springIsEntanglement.any()) {
    nAtoms += static_cast<double>((this->linkIsEntanglement == true).count());
  }
  return nAtoms;
}

void
MEHPNetwork::recomputeMeanContourLengths()
{
  this->meanSpringContourLength =
    (this->nrOfSprings > 0 && this->springContourLength.size() > 0)
      ? this->springContourLength.mean()
      : 0.0;
  this->meanStrandContourLength =
    (this->nrOfStrands > 0 && this->strandContourLength.size() > 0)
      ? this->strandContourLength.mean()
      : 0.0;
}

// =============================================================================
// Slip-Link Insertion & Kinematics
// =============================================================================

double
MEHPNetwork::computeClampedInverseSpringContourLength(
  const bool isPartialSpring,
  const double springContourLen,
  const double strandContourLen,
  const double oneOverSpringPartitionUpperLimit)
{
  if (!isPartialSpring) {
    return (springContourLen > 0.0) ? (1.0 / springContourLen)
                                    : oneOverSpringPartitionUpperLimit;
  }
  if (oneOverSpringPartitionUpperLimit > 0.0) {
    const double upperLimit = oneOverSpringPartitionUpperLimit;
    const double remainingLength = strandContourLen - (1.0 / upperLimit);
    const double lowerLimit = (remainingLength > 0.0)
                                ? std::min(1.0 / remainingLength, upperLimit)
                                : upperLimit;
    const double rawInv =
      (springContourLen > 0.0) ? (1.0 / springContourLen) : upperLimit;
    return std::clamp(rawInv, lowerLimit, upperLimit);
  }
  if (springContourLen > 0.0) {
    const double rawInv = 1.0 / springContourLen;
    if (std::isfinite(rawInv)) {
      return rawInv;
    }
  }
  return std::max(strandContourLen, 1.0);
}

size_t
MEHPNetwork::addSlipLinks(
  const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
    entanglements,
  Eigen::VectorXd* u,
  const bool clampAlpha)
{
  return this->addSlipLinks(entanglements.strandIdx1,
                            entanglements.strandIdx2,
                            entanglements.x,
                            entanglements.y,
                            entanglements.z,
                            entanglements.alpha1,
                            entanglements.alpha2,
                            u,
                            clampAlpha);
}

size_t
MEHPNetwork::addSlipLinks(const std::vector<size_t>& strandIdx1,
                          const std::vector<size_t>& strandIdx2,
                          const Eigen::VectorXd& x,
                          const Eigen::VectorXd& y,
                          const Eigen::VectorXd& z,
                          Eigen::VectorXd* u,
                          const bool clampAlpha)
{
  const std::vector<double> halfAlphas(strandIdx1.size(), 0.5);
  return this->addSlipLinks(
    strandIdx1, strandIdx2, x, y, z, halfAlphas, halfAlphas, u, clampAlpha);
}

size_t
MEHPNetwork::addSlipLinks(const std::vector<size_t>& strandIdx1,
                          const std::vector<size_t>& strandIdx2,
                          const Eigen::VectorXd& x,
                          const Eigen::VectorXd& y,
                          const Eigen::VectorXd& z,
                          const std::vector<double>& alpha1,
                          const std::vector<double>& alpha2,
                          Eigen::VectorXd* u,
                          const bool clampAlpha)
{
  const size_t nNewSlipLinks = strandIdx1.size();
  INVALIDARG_EXP_IFN(strandIdx2.size() == nNewSlipLinks &&
                       static_cast<size_t>(x.size()) == nNewSlipLinks &&
                       static_cast<size_t>(y.size()) == nNewSlipLinks &&
                       static_cast<size_t>(z.size()) == nNewSlipLinks &&
                       alpha1.size() == nNewSlipLinks &&
                       alpha2.size() == nNewSlipLinks,
                     "Input vector sizes to addSlipLinks must match.");
  if (u != nullptr) {
    INVALIDARG_EXP_IFN(
      u->size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
      "Displacement vector u must match existing nrOfLinks before "
      "addSlipLinks.");
  }

  for (size_t i = 0; i < nNewSlipLinks; ++i) {
    INVALIDARG_EXP_IFN(strandIdx1[i] < this->nrOfStrands &&
                         strandIdx2[i] < this->nrOfStrands,
                       "Strand index out of range in addSlipLinks.");
    INVALIDARG_EXP_IFN(!this->isEntanglementStrand(strandIdx1[i]) &&
                         !this->isEntanglementStrand(strandIdx2[i]),
                       "Cannot place slip-link on an entanglement coupling "
                       "strand.");
    if (!clampAlpha) {
      INVALIDARG_EXP_IFN(APPROX_WITHIN(alpha1[i], 0.0, 1.0, 1e-9) &&
                           APPROX_WITHIN(alpha2[i], 0.0, 1.0, 1e-9),
                         "Slip-link alpha fraction must lie within [0, 1].");
    }
  }

  if (nNewSlipLinks == 0) {
    return 0;
  }

  const size_t oldNrOfLinks = this->nrOfLinks;
  const size_t oldNrOfSprings = this->nrOfSprings;
  const Eigen::Index newTotalLinksE =
    static_cast<Eigen::Index>(oldNrOfLinks + nNewSlipLinks);
  const Eigen::Index newTotalSpringsE =
    static_cast<Eigen::Index>(oldNrOfSprings + 2 * nNewSlipLinks);

  this->coordinates.conservativeResize(3 * newTotalLinksE);
  this->linkIsEntanglement.conservativeResize(newTotalLinksE);
  this->oldAtomIds.conservativeResize(newTotalLinksE);
  this->oldAtomTypes.conservativeResize(newTotalLinksE);
  this->nrOfCrosslinkSwapsEndured.conservativeResize(newTotalLinksE);
  this->strandIndicesOfLink.resize(oldNrOfLinks + nNewSlipLinks);

  this->linkIsEntanglement.tail(static_cast<Eigen::Index>(nNewSlipLinks)) =
    true;
  this->oldAtomIds.tail(static_cast<Eigen::Index>(nNewSlipLinks)) = -1;
  this->oldAtomTypes.tail(static_cast<Eigen::Index>(nNewSlipLinks)) = 1;
  this->nrOfCrosslinkSwapsEndured.tail(
    static_cast<Eigen::Index>(nNewSlipLinks)) = 0;

  this->springIndexA.conservativeResize(newTotalSpringsE);
  this->springIndexB.conservativeResize(newTotalSpringsE);
  this->springCoordinateIndexA.conservativeResize(3 * newTotalSpringsE);
  this->springCoordinateIndexB.conservativeResize(3 * newTotalSpringsE);
  this->springBoxOffset.conservativeResize(3 * newTotalSpringsE);
  this->springBoxOffset.tail(static_cast<Eigen::Index>(6 * nNewSlipLinks))
    .setZero();
  this->springContourLength.conservativeResize(newTotalSpringsE);
  this->springPartitions.conservativeResize(newTotalSpringsE);
  this->strandIndexOfSpring.conservativeResize(newTotalSpringsE);
  this->springIsEntanglement.conservativeResize(newTotalSpringsE);
  this->springIsEntanglement.tail(
    static_cast<Eigen::Index>(2 * nNewSlipLinks)) = false;

  Eigen::VectorXd localU;
  Eigen::VectorXd* activeU = u;
  if (activeU != nullptr) {
    activeU->conservativeResize(3 * newTotalLinksE);
    activeU->tail(static_cast<Eigen::Index>(3 * nNewSlipLinks)).setZero();
  } else {
    localU = Eigen::VectorXd::Zero(3 * newTotalLinksE);
    activeU = &localU;
  }

  this->nrOfLinks += nNewSlipLinks;
  this->nrOfSprings += 2 * nNewSlipLinks;

  size_t nextSpringIdx = oldNrOfSprings;
  for (size_t i = 0; i < nNewSlipLinks; ++i) {
    const size_t newLinkIdx = oldNrOfLinks + i;
    const Eigen::Index iE = static_cast<Eigen::Index>(i);
    Eigen::Vector3d linkPos(x[iE], y[iE], z[iE]);
    this->box.handlePBC(linkPos);
    this->coordinates.segment<3>(3 * static_cast<Eigen::Index>(newLinkIdx)) =
      linkPos;

    const std::array<size_t, 2> targetStrands = { strandIdx1[i],
                                                  strandIdx2[i] };
    const std::array<double, 2> targetAlphas = { alpha1[i], alpha2[i] };

    for (size_t leg = 0; leg < 2; ++leg) {
      const size_t strandIdx = targetStrands[leg];
      const std::vector<size_t>& strandSprings =
        this->springIndicesOfStrand[strandIdx];
      RUNTIME_EXP_IFN(!strandSprings.empty(),
                      "Target strand has no springs in addSlipLinks.");

      double alpha = targetAlphas[leg];
      if (clampAlpha) {
        const double totalN =
          this->strandContourLength[static_cast<Eigen::Index>(strandIdx)];
        const double epsAlpha =
          (totalN > 0.0) ? std::min(1e-6, 0.25 / std::max(totalN, 1.0)) : 1e-6;
        alpha = std::clamp(alpha, epsAlpha, 1.0 - epsAlpha);
      } else {
        alpha = std::clamp(alpha, 1e-9, 1.0 - 1e-9);
      }

      double summedAlpha = 0.0;
      size_t subSpringPos = strandSprings.size() - 1;
      for (size_t k = 0; k < strandSprings.size(); ++k) {
        const double segAlpha =
          this->springPartitions[static_cast<Eigen::Index>(strandSprings[k])];
        if (summedAlpha + segAlpha >= alpha || k + 1 == strandSprings.size()) {
          subSpringPos = k;
          break;
        }
        summedAlpha += segAlpha;
      }

      const size_t dividedSpringIdx = strandSprings[subSpringPos];
      const double segAlpha =
        this->springPartitions[static_cast<Eigen::Index>(dividedSpringIdx)];
      double localAlpha = alpha - summedAlpha;
      const double minLocal = std::min(1e-7, 0.25 * segAlpha);
      if (clampAlpha || APPROX_WITHIN(localAlpha, 0.0, segAlpha, 1e-9)) {
        localAlpha = std::clamp(localAlpha, minLocal, segAlpha - minLocal);
      }

      const size_t newSpringIdx = nextSpringIdx++;
      insertLinkOntoSubSpring(*this,
                              *activeU,
                              strandIdx,
                              subSpringPos,
                              dividedSpringIdx,
                              newSpringIdx,
                              newLinkIdx,
                              localAlpha,
                              false);
    }
  }

  assert(nextSpringIdx == this->nrOfSprings);
  this->recomputeMeanContourLengths();
  return nNewSlipLinks;
}

size_t
MEHPNetwork::randomlyAddSlipLinks(const size_t nrOfEntanglementsToSample,
                                  const double cutoff,
                                  const size_t minNrOfEntanglements,
                                  const double sameStrandCutoff,
                                  const int crossLinkerType,
                                  const std::string& seed,
                                  Eigen::VectorXd* u,
                                  const bool clampAlpha)
{
  if (nrOfEntanglementsToSample == 0 || this->nrOfStrands == 0) {
    return 0;
  }
  if (u != nullptr) {
    INVALIDARG_EXP_IFN(
      u->size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
      "Displacement vector u size mismatch in randomlyAddSlipLinks.");
  }

  // Collect polymer strands (non-entanglement strands with endpoints)
  std::vector<size_t> includedNetworkStrands;
  includedNetworkStrands.reserve(this->nrOfStrands);
  std::vector<long int> linkToCuNode(this->nrOfLinks, -1);
  std::vector<size_t> cuNodeToLink;

  for (size_t s = 0; s < this->nrOfStrands; ++s) {
    if (this->isEntanglementStrand(s) ||
        this->linkIndicesOfStrand[s].size() < 2) {
      continue;
    }
    const size_t endA = this->linkIndicesOfStrand[s].front();
    const size_t endB = this->linkIndicesOfStrand[s].back();
    if (this->linkIsEntanglement[static_cast<Eigen::Index>(endA)] ||
        this->linkIsEntanglement[static_cast<Eigen::Index>(endB)]) {
      continue;
    }
    includedNetworkStrands.push_back(s);
    if (linkToCuNode[endA] < 0) {
      linkToCuNode[endA] = static_cast<long int>(cuNodeToLink.size());
      cuNodeToLink.push_back(endA);
    }
    if (linkToCuNode[endB] < 0) {
      linkToCuNode[endB] = static_cast<long int>(cuNodeToLink.size());
      cuNodeToLink.push_back(endB);
    }
  }

  if (includedNetworkStrands.empty()) {
    return 0;
  }

  const size_t nCuNodes = cuNodeToLink.size();
  const size_t nCuStrands = includedNetworkStrands.size();

  CrosslinkerUniverse cu;
  cu.xlinkX.resize(nCuNodes);
  cu.xlinkY.resize(nCuNodes);
  cu.xlinkZ.resize(nCuNodes);
  cu.xlinkTypes.resize(nCuNodes, crossLinkerType);
  std::vector<size_t> nodeDegree(nCuNodes, 0);

  for (size_t cuNode = 0; cuNode < nCuNodes; ++cuNode) {
    const size_t lIdx = cuNodeToLink[cuNode];
    const Eigen::Index lIdxE = static_cast<Eigen::Index>(lIdx);
    const Eigen::Vector3d pos =
      this->coordinates.segment<3>(3 * lIdxE) + displacementAt(u, lIdx);
    cu.xlinkX[cuNode] = pos[0];
    cu.xlinkY[cuNode] = pos[1];
    cu.xlinkZ[cuNode] = pos[2];
    if (lIdxE < this->oldAtomTypes.size()) {
      cu.xlinkTypes[cuNode] = this->oldAtomTypes[lIdxE];
    }
  }

  for (const size_t sIdx : includedNetworkStrands) {
    const size_t cuA = static_cast<size_t>(
      linkToCuNode[this->linkIndicesOfStrand[sIdx].front()]);
    const size_t cuB =
      static_cast<size_t>(linkToCuNode[this->linkIndicesOfStrand[sIdx].back()]);
    nodeDegree[cuA] += 1;
    nodeDegree[cuB] += 1;
  }

  cu.strandFrom.resize(nCuStrands);
  cu.strandTo.resize(nCuStrands);
  cu.beadsInStrand.resize(nCuStrands);
  const double b =
    (this->defaultBondLength > 0.0) ? this->defaultBondLength : 1.0;
  const double b2 = (this->defaultR02 > 0.0) ? this->defaultR02 : (b * b);
  cu.beadDistanceInStrand.assign(nCuStrands, b);
  cu.meanSquaredBeadDistanceInStrand.assign(nCuStrands, b2);

  std::vector<bool> strandIsReversed(nCuStrands, false);
  for (size_t k = 0; k < nCuStrands; ++k) {
    const size_t sIdx = includedNetworkStrands[k];
    const size_t cuA = static_cast<size_t>(
      linkToCuNode[this->linkIndicesOfStrand[sIdx].front()]);
    const size_t cuB =
      static_cast<size_t>(linkToCuNode[this->linkIndicesOfStrand[sIdx].back()]);
    const bool aIsFree =
      (nodeDegree[cuA] == 1 && cu.xlinkTypes[cuA] != crossLinkerType);
    const bool bIsFree =
      (nodeDegree[cuB] == 1 && cu.xlinkTypes[cuB] != crossLinkerType);
    const int nBonds = static_cast<int>(
      std::round(this->strandContourLength[static_cast<Eigen::Index>(sIdx)]));
    const int nFreeEnds = static_cast<int>(aIsFree) + static_cast<int>(bIsFree);
    cu.beadsInStrand[k] = std::max(0, nBonds - 1 + nFreeEnds);
    if (aIsFree && !bIsFree) {
      cu.strandFrom[k] = static_cast<long int>(cuB);
      cu.strandTo[k] = -1;
      strandIsReversed[k] = true;
    } else {
      cu.strandFrom[k] = aIsFree ? -1 : static_cast<long int>(cuA);
      cu.strandTo[k] = bIsFree ? -1 : static_cast<long int>(cuB);
    }
  }

  pylimer_tools::topo::entanglement_detection::StrandPairEntanglements sampled =
    pylimer_tools::topo::entanglement_detection::
      randomlyFindEntanglementsOnDemand(this->box,
                                        cu,
                                        nrOfEntanglementsToSample,
                                        cutoff,
                                        0.0,
                                        minNrOfEntanglements,
                                        sameStrandCutoff,
                                        seed,
                                        false);
  for (size_t e = 0; e < sampled.size(); ++e) {
    const size_t localS1 = sampled.strandIdx1[e];
    const size_t localS2 = sampled.strandIdx2[e];
    if (strandIsReversed[localS1]) {
      sampled.alpha1[e] = 1.0 - sampled.alpha1[e];
    }
    if (strandIsReversed[localS2]) {
      sampled.alpha2[e] = 1.0 - sampled.alpha2[e];
    }
    sampled.strandIdx1[e] = includedNetworkStrands[localS1];
    sampled.strandIdx2[e] = includedNetworkStrands[localS2];
  }

  return this->addSlipLinks(sampled, u, clampAlpha);
}

void
MEHPNetwork::realignSlipLinkToImages(const Eigen::VectorXd& u,
                                     const size_t slipLinkIdx,
                                     const size_t springIdx1,
                                     const size_t springIdx2,
                                     const bool is2D)
{
  const Eigen::Index s1E = static_cast<Eigen::Index>(springIdx1);
  const Eigen::Index s2E = static_cast<Eigen::Index>(springIdx2);
  const Eigen::Index slipE = static_cast<Eigen::Index>(slipLinkIdx);
  INVALIDARG_EXP_IFN(springIdx1 < this->nrOfSprings &&
                       springIdx2 < this->nrOfSprings &&
                       slipLinkIdx < this->nrOfLinks,
                     "Index out of range in realignSlipLinkToImages.");
  INVALIDARG_EXP_IFN(this->springIndexB[s1E] == static_cast<int>(slipLinkIdx) &&
                       this->springIndexA[s2E] == static_cast<int>(slipLinkIdx),
                     "Sub-springs must meet head-to-tail at the slip-link.");

  const Eigen::Vector3d totalOffset =
    this->springBoxOffset.segment<3>(3 * s1E) +
    this->springBoxOffset.segment<3>(3 * s2E);
  const Eigen::Vector3d totalDistanceBefore =
    this->evaluateSpringVector(u, springIdx1, is2D) +
    this->evaluateSpringVector(u, springIdx2, is2D);

  const Eigen::Index srcIdx =
    static_cast<Eigen::Index>(this->springIndexA[s1E]);
  const Eigen::Index dstIdx =
    static_cast<Eigen::Index>(this->springIndexB[s2E]);
  const Eigen::Vector3d sourceCoords =
    this->coordinates.segment<3>(3 * srcIdx) + u.segment<3>(3 * srcIdx);
  const Eigen::Vector3d viaCoords =
    this->coordinates.segment<3>(3 * slipE) + u.segment<3>(3 * slipE);
  const Eigen::Vector3d targetCoords =
    this->coordinates.segment<3>(3 * dstIdx) + u.segment<3>(3 * dstIdx);

  const Eigen::Vector3d aVec = viaCoords - sourceCoords;
  const Eigen::Vector3d bVec = (targetCoords - viaCoords) + totalOffset;

  Eigen::Vector3d bestOffset = Eigen::Vector3d::Zero();
  const size_t nActiveDirs = is2D ? 2 : 3;
  for (size_t dir = 0; dir < nActiveDirs; ++dir) {
    const Eigen::Index d = static_cast<Eigen::Index>(dir);
    const double Ld = this->L[dir];
    const double mStar = (bVec[d] - aVec[d]) / (2.0 * Ld);
    const long int mRound = std::lround(mStar);
    double bestScore = std::numeric_limits<double>::infinity();
    long int bestM = mRound;
    for (long int m = mRound - 1; m <= mRound + 1; ++m) {
      const double off = static_cast<double>(m) * Ld;
      const double v1 = aVec[d] + off;
      const double v2 = bVec[d] - off;
      const double score = v1 * v1 + v2 * v2;
      if (score < bestScore) {
        bestScore = score;
        bestM = m;
      }
    }
    bestOffset[d] = static_cast<double>(bestM) * Ld;
  }

  this->springBoxOffset.segment<3>(3 * s1E) = bestOffset;
  this->springBoxOffset.segment<3>(3 * s2E) = totalOffset - bestOffset;

  const Eigen::Vector3d totalDistanceAfter =
    this->evaluateSpringVector(u, springIdx1, is2D) +
    this->evaluateSpringVector(u, springIdx2, is2D);
  assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
    totalDistanceBefore, totalDistanceAfter, 1e-7));
}

double
MEHPNetwork::updateSpringPartitionAtLink(
  const size_t linkIdx,
  const Eigen::VectorXd& u,
  const double oneOverSpringPartitionUpperLimit,
  const bool is2D)
{
  INVALIDARG_EXP_IFN(linkIdx < this->nrOfLinks,
                     "Link index out of range in updateSpringPartitionAtLink.");
  INVALIDARG_EXP_IFN(
    this->linkIsEntanglement[static_cast<Eigen::Index>(linkIdx)],
    "Only entanglement/slip-links may adjust spring partitions.");

  double residualNorm = 0.0;
  size_t contributions = 0;

  for (const size_t strandIdx : this->strandIndicesOfLink[linkIdx]) {
    const std::vector<size_t>& links = this->linkIndicesOfStrand[strandIdx];
    if (links.size() < 3) {
      continue;
    }
    for (size_t partnerPos = 1; partnerPos + 1 < links.size(); ++partnerPos) {
      if (links[partnerPos] != linkIdx) {
        continue;
      }
      const size_t prevSpringIdx =
        this->springIndicesOfStrand[strandIdx][partnerPos - 1];
      const size_t nextSpringIdx =
        this->springIndicesOfStrand[strandIdx][partnerPos];
      const Eigen::Index prevE = static_cast<Eigen::Index>(prevSpringIdx);
      const Eigen::Index nextE = static_cast<Eigen::Index>(nextSpringIdx);

      const Eigen::Vector3d vecBack =
        this->evaluateSpringVectorFrom(u, prevSpringIdx, linkIdx, is2D);
      const double distanceBack = vecBack.squaredNorm();
      const Eigen::Vector3d vecForward =
        this->evaluateSpringVectorFrom(u, nextSpringIdx, linkIdx, is2D);
      const double distanceForward = vecForward.squaredNorm();

      double idealValue = 0.5;
      if (distanceBack <= 0.0 && distanceForward <= 0.0) {
        idealValue = 0.5;
      } else if (distanceBack <= 0.0) {
        idealValue = 0.0;
      } else if (distanceForward <= 0.0) {
        idealValue = 1.0;
      } else {
        idealValue = 1.0 / (1.0 + std::sqrt(distanceForward / distanceBack));
      }

      const double currentS = this->springPartitions[prevE];
      const double nextS = this->springPartitions[nextE];
      const double totalL = currentS + nextS;
      const double strandN =
        this->strandContourLength[static_cast<Eigen::Index>(strandIdx)];

      if (oneOverSpringPartitionUpperLimit > 0.0) {
        const double limit =
          (totalL * strandN > 0.0)
            ? std::clamp(
                1.0 / (oneOverSpringPartitionUpperLimit * totalL * strandN),
                0.0,
                0.5)
            : 0.5;
        idealValue = std::clamp(idealValue, limit, 1.0 - limit);
      } else {
        idealValue = std::clamp(idealValue, 1e-12, 1.0 - 1e-12);
      }

      const double idealM1 = 1.0 - idealValue;
      const double newS = idealValue * totalL;
      const double compS = idealM1 * totalL;

      double localRes = 0.0;
      if (idealValue > 0.0 && idealValue < 1.0 && strandN * totalL > 0.0) {
        const double iv2 = idealValue * idealValue;
        const double im12 = idealM1 * idealM1;
        localRes = std::fma(distanceBack, im12, -distanceForward * iv2) /
                   (im12 * iv2 * strandN * totalL);
      }
      residualNorm += std::abs(localRes);
      contributions += 2;

      this->springPartitions[prevE] = newS;
      this->springPartitions[nextE] = compS;
      this->springContourLength[prevE] = newS * strandN;
      this->springContourLength[nextE] = compS * strandN;
    }
  }

  return (contributions > 0)
           ? (residualNorm / static_cast<double>(contributions))
           : 0.0;
}

double
MEHPNetwork::displaceToMeanPosition(
  Eigen::VectorXd& u,
  const bool is2D,
  const double oneOverSpringPartitionUpperLimit) const
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
    "Displacement vector u size must match network nrOfLinks.");
  if (this->nrOfLinks == 0 || this->nrOfSprings == 0) {
    return 0.0;
  }

  const Eigen::Index nSpringsE = static_cast<Eigen::Index>(this->nrOfSprings);
  const Eigen::Index nLinksE = static_cast<Eigen::Index>(this->nrOfLinks);
  Eigen::ArrayXd springWeights = Eigen::ArrayXd::Zero(3 * nSpringsE);

  for (size_t i = 0; i < this->nrOfSprings; ++i) {
    const Eigen::Index sIdx = static_cast<Eigen::Index>(i);
    const int strandIdxInt = this->strandIndexOfSpring[sIdx];
    const bool isPartial =
      (strandIdxInt >= 0 && !this->springIsEntanglement[sIdx] &&
       this->springIndicesOfStrand[static_cast<size_t>(strandIdxInt)].size() >
         1);
    const double strandN =
      (strandIdxInt >= 0)
        ? this->strandContourLength[static_cast<Eigen::Index>(strandIdxInt)]
        : this->springContourLength[sIdx];
    const double w = computeClampedInverseSpringContourLength(
      isPartial,
      this->springContourLength[sIdx],
      strandN,
      oneOverSpringPartitionUpperLimit);
    springWeights.segment<3>(3 * sIdx) = w;
  }

  const Eigen::ArrayXd springDistances =
    this->evaluateSpringVectors(u, is2D).array();
  Eigen::ArrayXd objectiveDisp = Eigen::ArrayXd::Zero(3 * nLinksE);
  objectiveDisp(this->springCoordinateIndexA) +=
    springWeights * springDistances;
  objectiveDisp(this->springCoordinateIndexB) -=
    springWeights * springDistances;

  const Eigen::ArrayXd nonLoopMask =
    (this->springCoordinateIndexA != this->springCoordinateIndexB)
      .cast<double>();
  Eigen::ArrayXd weightSum = Eigen::ArrayXd::Zero(3 * nLinksE);
  weightSum(this->springCoordinateIndexA) += springWeights * nonLoopMask;
  weightSum(this->springCoordinateIndexB) += springWeights * nonLoopMask;
  weightSum =
    weightSum.unaryExpr([](const double v) { return v > 0.0 ? v : 1.0; });

  const Eigen::ArrayXd remainingDisp = objectiveDisp / weightSum;

  Eigen::ArrayXd nSpringsPerLink = Eigen::ArrayXd::Zero(3 * nLinksE);
  nSpringsPerLink(this->springCoordinateIndexA) += nonLoopMask;
  nSpringsPerLink(this->springCoordinateIndexB) += nonLoopMask;
  nSpringsPerLink =
    nSpringsPerLink.unaryExpr([](const double v) { return v > 0.0 ? v : 1.0; });

  Eigen::ArrayXd backForthDisp = Eigen::ArrayXd::Zero(3 * nLinksE);
  backForthDisp(this->springCoordinateIndexA) +=
    nonLoopMask * (remainingDisp(this->springCoordinateIndexB) /
                   (nSpringsPerLink(this->springCoordinateIndexA) * 2.0));
  backForthDisp(this->springCoordinateIndexB) +=
    nonLoopMask * (remainingDisp(this->springCoordinateIndexA) /
                   (nSpringsPerLink(this->springCoordinateIndexB) * 2.0));

  const Eigen::VectorXd finalDisp = (remainingDisp + backForthDisp).matrix();
  RUNTIME_EXP_IFN(
    pylimer_tools::utils::all_components_finite(finalDisp),
    "Non-finite displacement computed in displaceToMeanPosition.");
  u += finalDisp;
  return pylimer_tools::utils::segmentwise_norm_max(finalDisp, 3);
}

double
MEHPNetwork::displaceToMeanPosition(
  Eigen::VectorXd& u,
  const size_t linkIdx,
  const bool is2D,
  const double oneOverSpringPartitionUpperLimit) const
{
  const std::vector<size_t> springIndices =
    this->getSpringIndicesOfLink(linkIdx);
  return this->displaceToMeanPosition(
    u, linkIdx, is2D, oneOverSpringPartitionUpperLimit, springIndices);
}

double
MEHPNetwork::displaceToMeanPosition(
  Eigen::VectorXd& u,
  const size_t linkIdx,
  const bool is2D,
  const double oneOverSpringPartitionUpperLimit,
  const std::vector<size_t>& springIndicesOfLink) const
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
    "Displacement vector u size must match network nrOfLinks.");
  INVALIDARG_EXP_IFN(linkIdx < this->nrOfLinks,
                     "Link index out of range in displaceToMeanPosition.");

  Eigen::Vector3d weightedDisp = Eigen::Vector3d::Zero();
  double weightTotal = 0.0;

  for (const size_t spIdx : springIndicesOfLink) {
    INVALIDARG_EXP_IFN(spIdx < this->nrOfSprings,
                       "Spring index out of range in displaceToMeanPosition.");
    if (this->isLoopingSpring(spIdx)) {
      continue;
    }
    const Eigen::Index sE = static_cast<Eigen::Index>(spIdx);
    const int strandIdxInt = this->strandIndexOfSpring[sE];
    const bool isPartial =
      (strandIdxInt >= 0 && !this->springIsEntanglement[sE] &&
       this->springIndicesOfStrand[static_cast<size_t>(strandIdxInt)].size() >
         1);
    const double strandN =
      (strandIdxInt >= 0)
        ? this->strandContourLength[static_cast<Eigen::Index>(strandIdxInt)]
        : this->springContourLength[sE];
    const double w = computeClampedInverseSpringContourLength(
      isPartial,
      this->springContourLength[sE],
      strandN,
      oneOverSpringPartitionUpperLimit);

    weightedDisp += w * this->evaluateSpringVectorFrom(u, spIdx, linkIdx, is2D);
    weightTotal += w;
  }

  if (weightTotal <= 0.0) {
    return 0.0;
  }
  const Eigen::Vector3d delta = weightedDisp / weightTotal;
  u.segment<3>(3 * static_cast<Eigen::Index>(linkIdx)) += delta;
  return delta.norm();
}

bool
MEHPNetwork::mergeSpringsAtLink(const size_t linkIdx,
                                Eigen::VectorXd& u,
                                const bool removeIfZeroFunctional)
{
  INVALIDARG_EXP_IFN(linkIdx < this->nrOfLinks,
                     "Link index out of range in mergeSpringsAtLink.");
  const std::vector<size_t> springIndices =
    this->getSpringIndicesOfLink(linkIdx);
  if (springIndices.size() != 2) {
    return false;
  }
  if (this->isLoopingSpring(springIndices[0]) ||
      this->isLoopingSpring(springIndices[1])) {
    return false;
  }
  const Eigen::Index s0E = static_cast<Eigen::Index>(springIndices[0]);
  const Eigen::Index s1E = static_cast<Eigen::Index>(springIndices[1]);
  if (this->springIsEntanglement[s0E] || this->springIsEntanglement[s1E]) {
    return false;
  }
  if (!this->linkIsEntanglement[static_cast<Eigen::Index>(linkIdx)] &&
      this->strandIndexOfSpring[s0E] == this->strandIndexOfSpring[s1E]) {
    return false;
  }

  mergeTwoSpringsWithoutRemoval(
    *this, u, springIndices[0], springIndices[1], linkIdx);
  this->removeSprings({ springIndices[0] }, removeIfZeroFunctional, &u);
  return true;
}

void
MEHPNetwork::removeLink(const size_t linkIdx, Eigen::VectorXd* u)
{
  INVALIDARG_EXP_IFN(linkIdx < this->nrOfLinks,
                     "Link index out of range in removeLink.");
  INVALIDARG_EXP_IFN(
    this->strandIndicesOfLink[linkIdx].empty(),
    "Connected springs/strands must be removed or unlinked before removeLink.");

  std::vector<size_t> singleLink = { linkIdx };
  std::vector<size_t> linkCoords = { 3 * linkIdx + 2,
                                     3 * linkIdx + 1,
                                     3 * linkIdx };

  pylimer_tools::utils::removeRows(
    this->linkIsEntanglement, singleLink, true, true);
  pylimer_tools::utils::removeRows(
    this->strandIndicesOfLink, singleLink, true, true);
  pylimer_tools::utils::removeRows(this->oldAtomIds, singleLink, true, true);
  pylimer_tools::utils::removeRows(this->oldAtomTypes, singleLink, true, true);
  pylimer_tools::utils::removeRows(
    this->nrOfCrosslinkSwapsEndured, singleLink, true, true);

  pylimer_tools::utils::removeRows(this->coordinates, linkCoords, true, true);
  if (u != nullptr && u->size() > 0) {
    pylimer_tools::utils::removeRows(*u, linkCoords, true, true);
  }

  const std::vector<long int> newMapping =
    pylimer_tools::utils::getMappingForRenumbering(singleLink, this->nrOfLinks);
  safeRenumberWithMapping(this->springIndexA, newMapping);
  safeRenumberWithMapping(this->springIndexB, newMapping);
  pylimer_tools::utils::renumberWithMapping(this->linkIndicesOfStrand,
                                            newMapping);

  std::vector<long int> newCoordMapping;
  newCoordMapping.reserve(newMapping.size() * 3);
  for (size_t i = 0; i < newMapping.size(); ++i) {
    for (long int dir = 0; dir < 3; ++dir) {
      newCoordMapping.push_back(newMapping[i] >= 0 ? (newMapping[i] * 3 + dir)
                                                   : -1);
    }
  }
  safeRenumberWithMapping(this->springCoordinateIndexA, newCoordMapping);
  safeRenumberWithMapping(this->springCoordinateIndexB, newCoordMapping);

  this->nrOfLinks -= 1;
  this->nrOfNodes =
    static_cast<size_t>((this->linkIsEntanglement == false).count());
}

void
MEHPNetwork::removeLink(const size_t linkIdx, Eigen::VectorXd& u)
{
  this->removeLink(linkIdx, &u);
}

void
MEHPNetwork::removeSpring(const size_t springIdx,
                          const bool removeZeroFunctionalLinks,
                          Eigen::VectorXd* u)
{
  INVALIDARG_EXP_IFN(springIdx < this->nrOfSprings,
                     "Spring index out of range in removeSpring.");
  std::vector<size_t> toRemove;
  const int strandIdxInt =
    this->strandIndexOfSpring[static_cast<Eigen::Index>(springIdx)];
  if (strandIdxInt >= 0) {
    const size_t strandIdx = static_cast<size_t>(strandIdxInt);
    toRemove = this->springIndicesOfStrand[strandIdx];
    pylimer_tools::utils::append(
      toRemove, this->getEntanglementSpringsAtStrand(strandIdx));
  } else {
    toRemove.push_back(springIdx);
  }
  this->removeSprings(toRemove, removeZeroFunctionalLinks, u);
}

void
MEHPNetwork::swapSlipLinksOnSpring(const size_t springIdx,
                                   const Eigen::VectorXd& u,
                                   const bool is2D)
{
  INVALIDARG_EXP_IFN(springIdx < this->nrOfSprings,
                     "Spring index out of range in swapSlipLinksOnSpring.");
  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  INVALIDARG_EXP_IFN(!this->springIsEntanglement[sE],
                     "Cannot swap slip-links across an entanglement spring.");
  const size_t link1 = static_cast<size_t>(this->springIndexA[sE]);
  const size_t link2 = static_cast<size_t>(this->springIndexB[sE]);
  INVALIDARG_EXP_IFN(link1 != link2, "Cannot swap slip-link with itself.");
  INVALIDARG_EXP_IFN(
    this->linkIsEntanglement[static_cast<Eigen::Index>(link1)] &&
      this->linkIsEntanglement[static_cast<Eigen::Index>(link2)],
    "Both ends of spring must be slip-links to swap directly.");

  const size_t strandIdx = static_cast<size_t>(this->strandIndexOfSpring[sE]);
  const size_t pos = pylimer_tools::utils::index_of(
    this->springIndicesOfStrand[strandIdx], springIdx);
  RUNTIME_EXP_IFN(
    pos > 0 && pos + 1 < this->springIndicesOfStrand[strandIdx].size(),
    "Swappable sub-spring must lie strictly between two adjacent slip-links.");

  const size_t prevSpringIdx = this->springIndicesOfStrand[strandIdx][pos - 1];
  const size_t nextSpringIdx = this->springIndicesOfStrand[strandIdx][pos + 1];
  const Eigen::Index prevE = static_cast<Eigen::Index>(prevSpringIdx);
  const Eigen::Index nextE = static_cast<Eigen::Index>(nextSpringIdx);

  const Eigen::Vector3d distanceBefore =
    this->evaluateSpringVector(u, prevSpringIdx, is2D) +
    this->evaluateSpringVector(u, springIdx, is2D) +
    this->evaluateSpringVector(u, nextSpringIdx, is2D);

  this->springIndexB[prevE] = static_cast<int>(link2);
  this->springIndexA[nextE] = static_cast<int>(link1);
  setSpringCoordinateIndicesAt(
    *this, prevE, this->springIndexA[prevE], static_cast<int>(link2));
  setSpringCoordinateIndicesAt(
    *this, nextE, static_cast<int>(link1), this->springIndexB[nextE]);

  this->springBoxOffset.segment<3>(3 * prevE) +=
    this->springBoxOffset.segment<3>(3 * sE);
  this->springBoxOffset.segment<3>(3 * nextE) +=
    this->springBoxOffset.segment<3>(3 * sE);

  std::swap(this->springIndexA[sE], this->springIndexB[sE]);
  setSpringCoordinateIndicesAt(
    *this, sE, this->springIndexA[sE], this->springIndexB[sE]);
  this->springBoxOffset.segment<3>(3 * sE) *= -1.0;

  assert(this->linkIndicesOfStrand[strandIdx][pos] == link1);
  assert(this->linkIndicesOfStrand[strandIdx][pos + 1] == link2);
  this->linkIndicesOfStrand[strandIdx][pos] = link2;
  this->linkIndicesOfStrand[strandIdx][pos + 1] = link1;

  const Eigen::Vector3d distanceAfter =
    this->evaluateSpringVector(u, prevSpringIdx, is2D) +
    this->evaluateSpringVector(u, springIdx, is2D) +
    this->evaluateSpringVector(u, nextSpringIdx, is2D);
  assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
    distanceBefore, distanceAfter, 1e-7));
}

long int
MEHPNetwork::moveSlipLinkFromRailToRail(
  const size_t sourceSpringIdx,
  size_t targetSpringIdx,
  Eigen::VectorXd& u,
  const double oneOverSpringPartitionUpperLimit,
  const bool is2D)
{
  INVALIDARG_EXP_IFN(
    sourceSpringIdx != targetSpringIdx && sourceSpringIdx < this->nrOfSprings &&
      targetSpringIdx < this->nrOfSprings,
    "Invalid source/target spring indices in moveSlipLinkFromRailToRail.");

  const Eigen::Index srcE = static_cast<Eigen::Index>(sourceSpringIdx);
  const Eigen::Index tgtE = static_cast<Eigen::Index>(targetSpringIdx);
  const size_t srcA = static_cast<size_t>(this->springIndexA[srcE]);
  const size_t srcB = static_cast<size_t>(this->springIndexB[srcE]);
  const bool aIsSlip =
    this->linkIsEntanglement[static_cast<Eigen::Index>(srcA)];
  const bool bIsSlip =
    this->linkIsEntanglement[static_cast<Eigen::Index>(srcB)];
  INVALIDARG_EXP_IFN(
    XOR(aIsSlip, bIsSlip),
    "Source spring must connect exactly one slip-link and one crosslink.");

  const size_t involvedSlipLink = aIsSlip ? srcA : srcB;
  const size_t involvedCrossLink = aIsSlip ? srcB : srcA;
  INVALIDARG_EXP_IFN(this->isPartOfSpring(involvedCrossLink, targetSpringIdx),
                     "Target spring must touch the same crosslink.");

  const size_t srcStrandIdx =
    static_cast<size_t>(this->strandIndexOfSpring[srcE]);
  const size_t tgtStrandIdx =
    static_cast<size_t>(this->strandIndexOfSpring[tgtE]);
  if (this->isEntanglementStrand(srcStrandIdx) ||
      this->isEntanglementStrand(tgtStrandIdx)) {
    return -1;
  }

  const double tgtN =
    this->strandContourLength[static_cast<Eigen::Index>(tgtStrandIdx)];
  const double minAlpha =
    (oneOverSpringPartitionUpperLimit > 0.0 &&
     tgtN > 1.0 / oneOverSpringPartitionUpperLimit)
      ? 1.0 / (tgtN - 1.0 / oneOverSpringPartitionUpperLimit)
      : 1e-9;
  const size_t tgtSubSpringCount =
    this->springIndicesOfStrand[tgtStrandIdx].size();
  if (minAlpha * static_cast<double>(tgtSubSpringCount + 1) >= 1.0) {
    return -1;
  }

  const std::vector<size_t>& srcSprings =
    this->springIndicesOfStrand[srcStrandIdx];
  const size_t otherInvolvedSpring = (srcSprings.front() == sourceSpringIdx)
                                       ? srcSprings[1]
                                       : srcSprings[srcSprings.size() - 2];

  // Merge sourceSpringIdx into otherInvolvedSpring and remove sourceSpringIdx
  // while keeping link indices intact (`removeZeroFunctionalLinks = false`).
  mergeTwoSpringsWithoutRemoval(
    *this, u, sourceSpringIdx, otherInvolvedSpring, involvedSlipLink);
  this->removeSprings({ sourceSpringIdx }, false, &u);

  if (targetSpringIdx > sourceSpringIdx) {
    targetSpringIdx -= 1;
  }
  const size_t adjustedTgtStrandIdx = static_cast<size_t>(
    this->strandIndexOfSpring[static_cast<Eigen::Index>(targetSpringIdx)]);

  const bool crosslinkIsStartOfTarget =
    (this->springIndexA[static_cast<Eigen::Index>(targetSpringIdx)] ==
     static_cast<int>(involvedCrossLink));
  const size_t subPos = pylimer_tools::utils::index_of(
    this->springIndicesOfStrand[adjustedTgtStrandIdx], targetSpringIdx);

  // Append one new spring slot at the end of the spring arrays
  const size_t newSpringIdx = this->nrOfSprings;
  const Eigen::Index newTotalSpringsE =
    static_cast<Eigen::Index>(newSpringIdx + 1);
  this->nrOfSprings += 1;
  this->springIndexA.conservativeResize(newTotalSpringsE);
  this->springIndexB.conservativeResize(newTotalSpringsE);
  this->springCoordinateIndexA.conservativeResize(3 * newTotalSpringsE);
  this->springCoordinateIndexB.conservativeResize(3 * newTotalSpringsE);
  this->springBoxOffset.conservativeResize(3 * newTotalSpringsE);
  this->springBoxOffset.tail(3).setZero();
  this->springContourLength.conservativeResize(newTotalSpringsE);
  this->springPartitions.conservativeResize(newTotalSpringsE);
  this->strandIndexOfSpring.conservativeResize(newTotalSpringsE);
  this->springIsEntanglement.conservativeResize(newTotalSpringsE);

  const Eigen::Index tgtNewE = static_cast<Eigen::Index>(targetSpringIdx);
  const double oldSegAlpha = this->springPartitions[tgtNewE];
  const double newTgtN =
    this->strandContourLength[static_cast<Eigen::Index>(adjustedTgtStrandIdx)];

  // Ensure space by borrowing `minAlpha` from target strand segments if needed
  const double insertedAlpha = std::min(minAlpha, 0.5 * oldSegAlpha);
  double firstSegAlpha = 0.0;
  if (oldSegAlpha >= 2.0 * minAlpha) {
    firstSegAlpha =
      crosslinkIsStartOfTarget ? minAlpha : (oldSegAlpha - minAlpha);
  } else {
    // Temporarily give the split `insertedAlpha` and rescale remaining segments
    // proportionally so every sub-spring stays >= minAlpha and sums to 1.0.
    firstSegAlpha =
      crosslinkIsStartOfTarget ? insertedAlpha : (oldSegAlpha - insertedAlpha);
  }

  insertLinkOntoSubSpring(*this,
                          u,
                          adjustedTgtStrandIdx,
                          subPos,
                          targetSpringIdx,
                          newSpringIdx,
                          involvedSlipLink,
                          firstSegAlpha,
                          is2D);

  if (oldSegAlpha < 2.0 * minAlpha) {
    // Redistribute surplus above `minAlpha` from larger sibling sub-springs
    const size_t xAdjSpringIdx =
      crosslinkIsStartOfTarget ? targetSpringIdx : newSpringIdx;
    const Eigen::Index xAdjE = static_cast<Eigen::Index>(xAdjSpringIdx);
    const double deficit = minAlpha - this->springPartitions[xAdjE];
    if (deficit > 0.0) {
      double taken = 0.0;
      for (const size_t spIdx :
           this->springIndicesOfStrand[adjustedTgtStrandIdx]) {
        if (spIdx == xAdjSpringIdx) {
          continue;
        }
        const Eigen::Index sE = static_cast<Eigen::Index>(spIdx);
        const double avail = this->springPartitions[sE] - minAlpha;
        if (avail > 0.0) {
          const double take = std::min(avail, deficit - taken);
          this->springPartitions[sE] -= take;
          this->springContourLength[sE] = this->springPartitions[sE] * newTgtN;
          taken += take;
          if (taken >= deficit - 1e-12) {
            break;
          }
        }
      }
      this->springPartitions[xAdjE] += taken;
      this->springContourLength[xAdjE] =
        this->springPartitions[xAdjE] * newTgtN;
    }
  }

  this->recomputeMeanContourLengths();

  if (this->springIndexA[tgtNewE] == static_cast<int>(involvedCrossLink) &&
      this->springIndexB[tgtNewE] == static_cast<int>(involvedSlipLink)) {
    return static_cast<long int>(targetSpringIdx);
  }
  return static_cast<long int>(newSpringIdx);
}

long int
MEHPNetwork::rotateSlipLinkAroundCrosslink(
  const size_t springIdx,
  Eigen::VectorXd& u,
  const double oneOverSpringPartitionUpperLimit,
  const bool respectLoops,
  const bool is2D)
{
  INVALIDARG_EXP_IFN(
    springIdx < this->nrOfSprings,
    "Spring index out of range in rotateSlipLinkAroundCrosslink.");
  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  const size_t linkA = static_cast<size_t>(this->springIndexA[sE]);
  const size_t linkB = static_cast<size_t>(this->springIndexB[sE]);
  if (linkA == linkB) {
    return -1;
  }

  const bool aIsSlip =
    this->linkIsEntanglement[static_cast<Eigen::Index>(linkA)];
  const bool bIsSlip =
    this->linkIsEntanglement[static_cast<Eigen::Index>(linkB)];
  INVALIDARG_EXP_IFN(XOR(aIsSlip, bIsSlip),
                     "Exactly one end of spring must be a slip-link and one a "
                     "crosslink.");

  const size_t involvedCrossLink = aIsSlip ? linkB : linkA;
  const size_t srcStrandIdx =
    static_cast<size_t>(this->strandIndexOfSpring[sE]);

  std::vector<size_t> candidateTargetSprings;
  for (const size_t strandIdx : this->strandIndicesOfLink[involvedCrossLink]) {
    if (this->isEntanglementStrand(strandIdx) ||
        this->springIndicesOfStrand[strandIdx].empty()) {
      continue;
    }
    if (respectLoops && this->isLoopingStrand(strandIdx) &&
        strandIdx != srcStrandIdx) {
      continue;
    }
    const double tgtN =
      this->strandContourLength[static_cast<Eigen::Index>(strandIdx)];
    const double minAlpha =
      (oneOverSpringPartitionUpperLimit > 0.0 &&
       tgtN > 1.0 / oneOverSpringPartitionUpperLimit)
        ? 1.0 / (tgtN - 1.0 / oneOverSpringPartitionUpperLimit)
        : 1e-9;
    if (minAlpha * static_cast<double>(
                     this->springIndicesOfStrand[strandIdx].size() + 1) >=
        1.0) {
      continue;
    }
    const size_t endSpring =
      (this->linkIndicesOfStrand[strandIdx].front() == involvedCrossLink)
        ? this->springIndicesOfStrand[strandIdx].front()
        : this->springIndicesOfStrand[strandIdx].back();
    if (endSpring != springIdx &&
        this->isPartOfSpring(involvedCrossLink, endSpring)) {
      candidateTargetSprings.push_back(endSpring);
    }
  }

  if (candidateTargetSprings.empty()) {
    return -1;
  }

  return this->moveSlipLinkFromRailToRail(springIdx,
                                          candidateTargetSprings.front(),
                                          u,
                                          oneOverSpringPartitionUpperLimit,
                                          is2D);
}

bool
MEHPNetwork::moveSlipLinkAlongStrand(
  const size_t springIdx,
  Eigen::VectorXd& u,
  const LinkSwappingMode mode,
  const double oneOverSpringPartitionUpperLimit,
  const int nrOfCrosslinkSwapsAllowedPerSliplink,
  const bool is2D)
{
  if (mode == LinkSwappingMode::NO_SWAPPING || springIdx >= this->nrOfSprings) {
    return false;
  }

  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  if (this->springIsEntanglement[sE] || this->strandIndexOfSpring[sE] < 0) {
    return false;
  }

  const size_t partnerA = static_cast<size_t>(this->springIndexA[sE]);
  const size_t partnerB = static_cast<size_t>(this->springIndexB[sE]);
  if (partnerA == partnerB) {
    return false;
  }

  const bool aIsSlip =
    this->linkIsEntanglement[static_cast<Eigen::Index>(partnerA)];
  const bool bIsSlip =
    this->linkIsEntanglement[static_cast<Eigen::Index>(partnerB)];
  if (!aIsSlip && !bIsSlip) {
    return false;
  }

  const size_t strandIdx = static_cast<size_t>(this->strandIndexOfSpring[sE]);
  const double strandN =
    this->strandContourLength[static_cast<Eigen::Index>(strandIdx)];
  const double swappableCutoff =
    (oneOverSpringPartitionUpperLimit > 0.0 &&
     strandN > 1.0 / oneOverSpringPartitionUpperLimit)
      ? 1.0 / (strandN - 1.0 / oneOverSpringPartitionUpperLimit)
      : 1e-12;
  if (this->springPartitions[sE] > swappableCutoff) {
    return false;
  }

  const bool involvesCrosslink = (!aIsSlip || !bIsSlip);
  const bool respectLoops = (mode == LinkSwappingMode::ALL_CYCLE ||
                             mode == LinkSwappingMode::ALL_MC_CYCLE ||
                             mode == LinkSwappingMode::ALL_MC_TRY_CYCLE);
  const bool mcMode = (mode == LinkSwappingMode::ALL_MC ||
                       mode == LinkSwappingMode::ALL_MC_CYCLE ||
                       mode == LinkSwappingMode::ALL_MC_TRY ||
                       mode == LinkSwappingMode::ALL_MC_TRY_CYCLE);

  if (!involvesCrosslink) {
    // Two interior slip-links along the same strand
    const size_t pos = pylimer_tools::utils::index_of(
      this->springIndicesOfStrand[strandIdx], springIdx);
    if (pos == 0 || pos + 1 >= this->springIndicesOfStrand[strandIdx].size()) {
      return false;
    }

    if (!mcMode) {
      const size_t prevSpringIdx =
        this->springIndicesOfStrand[strandIdx][pos - 1];
      const size_t nextSpringIdx =
        this->springIndicesOfStrand[strandIdx][pos + 1];
      const Eigen::Index prevE = static_cast<Eigen::Index>(prevSpringIdx);
      const Eigen::Index nextE = static_cast<Eigen::Index>(nextSpringIdx);

      const double denomThis = computeClampedInverseSpringContourLength(
        true,
        this->springContourLength[sE],
        strandN,
        oneOverSpringPartitionUpperLimit);
      const double denomPrev = computeClampedInverseSpringContourLength(
        true,
        this->springContourLength[prevE],
        strandN,
        oneOverSpringPartitionUpperLimit);
      const double denomNext = computeClampedInverseSpringContourLength(
        true,
        this->springContourLength[nextE],
        strandN,
        oneOverSpringPartitionUpperLimit);

      const Eigen::Vector3d vecThis =
        this->evaluateSpringVector(u, springIdx, is2D);
      const Eigen::Vector3d vecPrev =
        this->evaluateSpringVector(u, prevSpringIdx, is2D);
      const Eigen::Vector3d vecNext =
        this->evaluateSpringVector(u, nextSpringIdx, is2D);

      const double forceBefore =
        (-vecPrev * denomPrev + vecThis * denomThis).squaredNorm() +
        (-vecThis * denomThis + vecNext * denomNext).squaredNorm();
      const double forceAfter =
        (vecThis * denomThis + (vecNext + vecThis) * denomNext).squaredNorm() +
        (-vecThis * denomThis - (vecThis + vecPrev) * denomPrev).squaredNorm();
      if (forceAfter > forceBefore) {
        return false;
      }
      this->swapSlipLinksOnSpring(springIdx, u, is2D);
      return true;
    }

    // MC reversible check: test swap, locally relax the pair, keep if stress
    // drops
    const std::vector<size_t> localLinks =
      collectLocalNeighborsForLinks(*this, { partnerA, partnerB });
    const double residualBefore =
      this
        ->evaluateStressTensorForLinks(
          localLinks, u, 1.0, is2D, oneOverSpringPartitionUpperLimit)
        .diagonal()
        .squaredNorm();
    const Eigen::VectorXd uSaved = u;

    this->swapSlipLinksOnSpring(springIdx, u, is2D);
    for (int step = 0; step < 3; ++step) {
      this->displaceToMeanPosition(
        u, partnerA, is2D, oneOverSpringPartitionUpperLimit);
      this->displaceToMeanPosition(
        u, partnerB, is2D, oneOverSpringPartitionUpperLimit);
    }

    const double residualAfter =
      this
        ->evaluateStressTensorForLinks(
          localLinks, u, 1.0, is2D, oneOverSpringPartitionUpperLimit)
        .diagonal()
        .squaredNorm();
    if (residualAfter <= residualBefore) {
      return true;
    }

    // Roll back
    this->swapSlipLinksOnSpring(springIdx, u, is2D);
    u = uSaved;
    return false;
  }

  // Slip-link adjacent to crosslink
  if (mode == LinkSwappingMode::SLIPLINKS_ONLY) {
    return false;
  }

  const size_t slipLinkIdx = aIsSlip ? partnerA : partnerB;
  const size_t xlinkIdx = aIsSlip ? partnerB : partnerA;
  const Eigen::Index slipE = static_cast<Eigen::Index>(slipLinkIdx);
  if (nrOfCrosslinkSwapsAllowedPerSliplink >= 0 &&
      this->nrOfCrosslinkSwapsEndured[slipE] >=
        nrOfCrosslinkSwapsAllowedPerSliplink) {
    return false;
  }

  const std::vector<size_t>& sSprings = this->springIndicesOfStrand[strandIdx];
  if (sSprings.size() < 2) {
    return false;
  }
  const size_t otherRailSpring = (sSprings.front() == springIdx)
                                   ? sSprings[1]
                                   : sSprings[sSprings.size() - 2];

  const bool srcAtFront = (sSprings.front() == springIdx);

  struct CandidateBranch
  {
    size_t strandIdx = 0;
    bool atFront = true;
    size_t currentEndSpring = 0;
  };

  // Collect candidate target branches attached to xlinkIdx
  std::vector<CandidateBranch> candidates;
  for (const size_t candStrandIdx : this->strandIndicesOfLink[xlinkIdx]) {
    if (this->isEntanglementStrand(candStrandIdx) ||
        this->springIndicesOfStrand[candStrandIdx].empty()) {
      continue;
    }
    if (respectLoops && this->isLoopingStrand(candStrandIdx) &&
        candStrandIdx != strandIdx) {
      continue;
    }
    const double tgtN =
      this->strandContourLength[static_cast<Eigen::Index>(candStrandIdx)];
    const double tgtMinAlpha =
      (oneOverSpringPartitionUpperLimit > 0.0 &&
       tgtN > 1.0 / oneOverSpringPartitionUpperLimit)
        ? 1.0 / (tgtN - 1.0 / oneOverSpringPartitionUpperLimit)
        : 1e-9;
    if (tgtMinAlpha *
          static_cast<double>(
            this->springIndicesOfStrand[candStrandIdx].size() + 1) >=
        1.0) {
      continue;
    }
    const bool candAtFront =
      (this->linkIndicesOfStrand[candStrandIdx].front() == xlinkIdx);
    const size_t endSpring =
      candAtFront ? this->springIndicesOfStrand[candStrandIdx].front()
                  : this->springIndicesOfStrand[candStrandIdx].back();
    if (endSpring != springIdx && this->isPartOfSpring(xlinkIdx, endSpring)) {
      candidates.push_back(
        CandidateBranch{ candStrandIdx, candAtFront, endSpring });
    }
  }
  if (candidates.empty()) {
    return false;
  }

  if (!mcMode) {
    Eigen::Vector3d thisVec = this->evaluateSpringVector(u, springIdx, is2D);
    Eigen::Vector3d otherVec =
      this->evaluateSpringVector(u, otherRailSpring, is2D);
    if (slipLinkIdx == partnerA) {
      thisVec *= -1.0;
      otherVec *= -1.0;
    }

    for (const CandidateBranch& cand : candidates) {
      const Eigen::Vector3d candVec = this->evaluateSpringVectorFrom(
        u, cand.currentEndSpring, xlinkIdx, is2D);
      const double forceBefore =
        (-otherVec + thisVec).squaredNorm() + (candVec - thisVec).squaredNorm();
      const double forceAfter =
        candVec.squaredNorm() + (-2.0 * thisVec - otherVec).squaredNorm();
      if (forceAfter < forceBefore) {
        const long int newSp =
          this->moveSlipLinkFromRailToRail(springIdx,
                                           cand.currentEndSpring,
                                           u,
                                           oneOverSpringPartitionUpperLimit,
                                           is2D);
        if (newSp >= 0) {
          this->nrOfCrosslinkSwapsEndured[slipE] += 1;
          return true;
        }
      }
    }
    return false;
  }

  const std::vector<size_t> localLinks =
    collectLocalNeighborsForLinks(*this, { slipLinkIdx, xlinkIdx });
  const double residualBefore =
    this
      ->evaluateStressTensorForLinks(
        localLinks, u, 1.0, is2D, oneOverSpringPartitionUpperLimit)
      .diagonal()
      .squaredNorm();
  const Eigen::VectorXd uSaved = u;
  const auto saveStrandPartitions =
    [this](const size_t sIdx) -> std::vector<double> {
    std::vector<double> saved;
    saved.reserve(this->springIndicesOfStrand[sIdx].size());
    for (const size_t sp : this->springIndicesOfStrand[sIdx]) {
      saved.push_back(this->springPartitions[static_cast<Eigen::Index>(sp)]);
    }
    return saved;
  };
  const auto restoreStrandPartitions =
    [this](const size_t sIdx, const std::vector<double>& saved) {
      const double totalN =
        this->strandContourLength[static_cast<Eigen::Index>(sIdx)];
      assert(saved.size() == this->springIndicesOfStrand[sIdx].size());
      for (size_t k = 0; k < saved.size(); ++k) {
        const Eigen::Index spE =
          static_cast<Eigen::Index>(this->springIndicesOfStrand[sIdx][k]);
        this->springPartitions[spE] = saved[k];
        this->springContourLength[spE] = saved[k] * totalN;
      }
    };
  const std::vector<double> savedSrcPartitions =
    saveStrandPartitions(strandIdx);

  size_t currentSrcSpringIdx = springIdx;
  for (const CandidateBranch& cand : candidates) {
    // Re-resolve endpoint spring index from strandIdx (stable across removals)
    const size_t candSpringIdx =
      cand.atFront ? this->springIndicesOfStrand[cand.strandIdx].front()
                   : this->springIndicesOfStrand[cand.strandIdx].back();
    const std::vector<double> savedCandPartitions =
      saveStrandPartitions(cand.strandIdx);

    const long int newSp =
      this->moveSlipLinkFromRailToRail(currentSrcSpringIdx,
                                       candSpringIdx,
                                       u,
                                       oneOverSpringPartitionUpperLimit,
                                       is2D);
    if (newSp < 0) {
      continue;
    }
    for (int step = 0; step < 3; ++step) {
      this->updateSpringPartitionAtLink(
        slipLinkIdx, u, oneOverSpringPartitionUpperLimit, is2D);
      this->displaceToMeanPosition(
        u, slipLinkIdx, is2D, oneOverSpringPartitionUpperLimit);
      this->displaceToMeanPosition(
        u, xlinkIdx, is2D, oneOverSpringPartitionUpperLimit);
    }
    const double residualAfter =
      this
        ->evaluateStressTensorForLinks(
          localLinks, u, 1.0, is2D, oneOverSpringPartitionUpperLimit)
        .diagonal()
        .squaredNorm();
    if (residualAfter <= residualBefore) {
      this->nrOfCrosslinkSwapsEndured[slipE] += 1;
      return true;
    }

    // Roll back move and restore exact strand partitions and displacements
    const size_t returnTarget =
      srcAtFront ? this->springIndicesOfStrand[strandIdx].front()
                 : this->springIndicesOfStrand[strandIdx].back();
    const long int restored =
      this->moveSlipLinkFromRailToRail(static_cast<size_t>(newSp),
                                       returnTarget,
                                       u,
                                       oneOverSpringPartitionUpperLimit,
                                       is2D);
    RUNTIME_EXP_IFN(restored >= 0,
                    "Failed to roll back rejected slip-link crosslink hop.");
    currentSrcSpringIdx = static_cast<size_t>(restored);
    restoreStrandPartitions(strandIdx, savedSrcPartitions);
    restoreStrandPartitions(cand.strandIdx, savedCandPartitions);
    this->recomputeMeanContourLengths();
    u = uSaved;

    if (mode != LinkSwappingMode::ALL_MC_TRY &&
        mode != LinkSwappingMode::ALL_MC_TRY_CYCLE) {
      return false;
    }
  }
  return false;
}

size_t
MEHPNetwork::performSlipLinkSwaps(
  Eigen::VectorXd& u,
  const LinkSwappingMode mode,
  const double oneOverSpringPartitionUpperLimit,
  const int nrOfCrosslinkSwapsAllowedPerSliplink,
  const bool is2D)
{
  if (mode == LinkSwappingMode::NO_SWAPPING) {
    return 0;
  }
  size_t totalSwaps = 0;
  for (long int spIdx = static_cast<long int>(this->nrOfSprings) - 1;
       spIdx >= 0;
       --spIdx) {
    if (static_cast<size_t>(spIdx) >= this->nrOfSprings) {
      continue;
    }
    if (this->moveSlipLinkAlongStrand(static_cast<size_t>(spIdx),
                                      u,
                                      mode,
                                      oneOverSpringPartitionUpperLimit,
                                      nrOfCrosslinkSwapsAllowedPerSliplink,
                                      is2D)) {
      totalSwaps += 1;
    }
  }
  return totalSwaps;
}

// =============================================================================
// Simplification, Pruning, Deformation & Matrix Update Support
// =============================================================================

std::vector<size_t>
MEHPNetwork::findTooLongStrandSprings(const Eigen::VectorXd& u,
                                      const double breakingLength,
                                      const bool is2D,
                                      size_t* outNumBrokenStrands) const
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
    "Displacement vector u size must match network nrOfLinks.");
  if (outNumBrokenStrands != nullptr) {
    *outNumBrokenStrands = 0;
  }
  if (breakingLength <= 0.0 || this->nrOfStrands == 0) {
    return {};
  }

  std::vector<size_t> springIndicesToDelete;
  size_t nBroken = 0;
  for (size_t strandIdx = 0; strandIdx < this->nrOfStrands; ++strandIdx) {
    for (const size_t springIdx : this->springIndicesOfStrand[strandIdx]) {
      const double weightedLen =
        this->evaluateWeightedSpringLength(u, springIdx, is2D);
      if (weightedLen > breakingLength) {
        nBroken += 1;
        for (const size_t sp : this->springIndicesOfStrand[strandIdx]) {
          springIndicesToDelete.push_back(sp);
        }
        pylimer_tools::utils::append(
          springIndicesToDelete,
          this->getEntanglementSpringsAtStrand(strandIdx));
        break;
      }
    }
  }

  if (outNumBrokenStrands != nullptr) {
    *outNumBrokenStrands = nBroken;
  }
  std::sort(springIndicesToDelete.begin(),
            springIndicesToDelete.end(),
            std::greater<size_t>());
  springIndicesToDelete.erase(
    std::unique(springIndicesToDelete.begin(), springIndicesToDelete.end()),
    springIndicesToDelete.end());
  return springIndicesToDelete;
}

size_t
MEHPNetwork::breakTooLongStrands(const Eigen::VectorXd& u,
                                 const double breakingLength,
                                 const bool is2D)
{
  size_t numBroken = 0;
  const std::vector<size_t> springsToRemove =
    this->findTooLongStrandSprings(u, breakingLength, is2D, &numBroken);
  if (!springsToRemove.empty()) {
    this->removeSprings(springsToRemove, false, nullptr);
  }
  return numBroken;
}

Eigen::ArrayXb
MEHPNetwork::findInactiveSprings(const Eigen::VectorXd& u,
                                 const double tolerance,
                                 const bool is2D,
                                 const bool checkAll) const
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
    "Displacement vector u size must match network nrOfLinks.");
  Eigen::ArrayXb result = Eigen::ArrayXb::Constant(
    static_cast<Eigen::Index>(this->nrOfSprings), false);
  if (this->nrOfSprings == 0) {
    return result;
  }

  const Eigen::ArrayXb activeStrands =
    this->findActiveStrands(u, tolerance, is2D);
  for (size_t strandIdx = 0; strandIdx < this->nrOfStrands; ++strandIdx) {
    if (!activeStrands[static_cast<Eigen::Index>(strandIdx)]) {
      for (const size_t spIdx : this->springIndicesOfStrand[strandIdx]) {
        result[static_cast<Eigen::Index>(spIdx)] = true;
      }
      const std::vector<size_t> entSprings =
        this->getEntanglementSpringsAtStrand(strandIdx);
      for (const size_t spIdx : entSprings) {
        result[static_cast<Eigen::Index>(spIdx)] = true;
      }
    }
  }

  if (checkAll) {
    for (size_t i = 0; i < this->nrOfSprings; ++i) {
      const Eigen::Index sIdx = static_cast<Eigen::Index>(i);
      if (this->strandIndexOfSpring[sIdx] < 0) {
        result[sIdx] = true;
      } else if (this->springIsEntanglement[sIdx]) {
        const Eigen::Vector3d vec = this->evaluateSpringVector(u, i, is2D);
        if (distanceIsWithinTolerance(
              vec, tolerance, this->springContourLength[sIdx])) {
          result[sIdx] = true;
        }
      }
    }
  }

  return result;
}

std::vector<BrokenSpringInfo>
MEHPNetwork::removeSpringsAndCollectMatrixUpdates(
  const std::vector<size_t>& springIndicesToRemove,
  const double kappa,
  Eigen::VectorXd* constants,
  const bool is2D)
{
  std::vector<BrokenSpringInfo> updates;
  if (springIndicesToRemove.empty()) {
    return updates;
  }

  std::vector<size_t> sortedDesc = springIndicesToRemove;
  std::sort(sortedDesc.begin(), sortedDesc.end(), std::greater<size_t>());
  sortedDesc.erase(std::unique(sortedDesc.begin(), sortedDesc.end()),
                   sortedDesc.end());

  if (constants != nullptr) {
    INVALIDARG_EXP_IFN(constants->size() ==
                         static_cast<Eigen::Index>(3 * this->nrOfLinks),
                       "Constants vector size must match 3 * nrOfLinks.");
  }

  updates.reserve(sortedDesc.size());
  for (const size_t spIdx : sortedDesc) {
    INVALIDARG_EXP_IFN(spIdx < this->nrOfSprings,
                       "Spring index out of range in "
                       "removeSpringsAndCollectMatrixUpdates.");
    const Eigen::Index sE = static_cast<Eigen::Index>(spIdx);
    const int idxAInt = this->springIndexA[sE];
    const int idxBInt = this->springIndexB[sE];
    if (idxAInt < 0 || idxBInt < 0 || idxAInt == idxBInt) {
      continue;
    }
    const size_t idxA = static_cast<size_t>(idxAInt);
    const size_t idxB = static_cast<size_t>(idxBInt);
    const double mult = kappa / this->springContourLength[sE];
    updates.push_back(BrokenSpringInfo{ idxA, idxB, mult });

    if (constants != nullptr) {
      const Eigen::Vector3d offset = this->springBoxOffset.segment<3>(3 * sE);
      constants->segment<3>(3 * static_cast<Eigen::Index>(idxA)) +=
        offset * mult;
      constants->segment<3>(3 * static_cast<Eigen::Index>(idxB)) -=
        offset * mult;
      if (is2D) {
        (*constants)[3 * static_cast<Eigen::Index>(idxA) + 2] = 0.0;
        (*constants)[3 * static_cast<Eigen::Index>(idxB) + 2] = 0.0;
      }
    }
  }

  this->removeSprings(sortedDesc, false, nullptr);
  return updates;
}

void
MEHPNetwork::removeSprings(const std::vector<size_t>& springIndicesToRemove,
                           const bool removeZeroFunctionalLinksFlag,
                           Eigen::VectorXd* u)
{
  if (springIndicesToRemove.empty()) {
    return;
  }

  std::vector<size_t> sortedSprings = springIndicesToRemove;
  std::sort(sortedSprings.begin(), sortedSprings.end(), std::greater<size_t>());
  sortedSprings.erase(std::unique(sortedSprings.begin(), sortedSprings.end()),
                      sortedSprings.end());

  INVALIDARG_EXP_IFN(sortedSprings.front() < this->nrOfSprings,
                     "Spring index out of bounds in removeSprings.");

  std::vector<size_t> strandsToDelete;
  for (const size_t spIdx : sortedSprings) {
    const int strandIdxInt =
      this->strandIndexOfSpring[static_cast<Eigen::Index>(spIdx)];
    if (strandIdxInt >= 0) {
      strandsToDelete.push_back(static_cast<size_t>(strandIdxInt));
    }
  }
  for (size_t strandIdx = 0; strandIdx < this->nrOfStrands; ++strandIdx) {
    if (this->linkIndicesOfStrand[strandIdx].empty()) {
      strandsToDelete.push_back(strandIdx);
    }
  }
  std::sort(
    strandsToDelete.begin(), strandsToDelete.end(), std::greater<size_t>());
  strandsToDelete.erase(
    std::unique(strandsToDelete.begin(), strandsToDelete.end()),
    strandsToDelete.end());

  for (const size_t strandIdx : strandsToDelete) {
    for (const size_t spIdx : this->springIndicesOfStrand[strandIdx]) {
      INVALIDARG_EXP_IFN(
        std::binary_search(sortedSprings.begin(),
                           sortedSprings.end(),
                           spIdx,
                           std::greater<size_t>()),
        "Spring " + std::to_string(spIdx) + " of strand " +
          std::to_string(strandIdx) +
          " must be removed alongside all other springs of the strand.");
    }
  }

  for (const size_t strandIdx : strandsToDelete) {
    for (const size_t linkIdx : this->linkIndicesOfStrand[strandIdx]) {
      pylimer_tools::utils::removeIfContained(
        this->strandIndicesOfLink[linkIdx], strandIdx);
    }
  }

  pylimer_tools::utils::removeRows(
    this->springContourLength, sortedSprings, true, true);
  pylimer_tools::utils::removeRows(
    this->springPartitions, sortedSprings, true, true);
  pylimer_tools::utils::removeRows(
    this->strandIndexOfSpring, sortedSprings, true, true);
  pylimer_tools::utils::removeRows(
    this->springIsEntanglement, sortedSprings, true, true);
  pylimer_tools::utils::removeRows(
    this->springIndexA, sortedSprings, true, true);
  pylimer_tools::utils::removeRows(
    this->springIndexB, sortedSprings, true, true);

  std::vector<size_t> toDeleteWithDirs;
  toDeleteWithDirs.reserve(sortedSprings.size() * 3);
  for (const size_t i : sortedSprings) {
    toDeleteWithDirs.push_back(i * 3 + 2);
    toDeleteWithDirs.push_back(i * 3 + 1);
    toDeleteWithDirs.push_back(i * 3);
  }
  pylimer_tools::utils::removeRows(
    this->springBoxOffset, toDeleteWithDirs, true, true);
  pylimer_tools::utils::removeRows(
    this->springCoordinateIndexA, toDeleteWithDirs, true, true);
  pylimer_tools::utils::removeRows(
    this->springCoordinateIndexB, toDeleteWithDirs, true, true);

  pylimer_tools::utils::removeRows(
    this->linkIndicesOfStrand, strandsToDelete, true, true);
  pylimer_tools::utils::removeRows(
    this->springIndicesOfStrand, strandsToDelete, true, true);
  pylimer_tools::utils::removeRows(
    this->strandContourLength, strandsToDelete, true, true);

  const std::vector<long int> newStrandMapping =
    pylimer_tools::utils::getMappingForRenumbering(strandsToDelete,
                                                   this->nrOfStrands);
  pylimer_tools::utils::renumberWithMapping(this->strandIndicesOfLink,
                                            newStrandMapping);
  safeRenumberWithMapping(this->strandIndexOfSpring, newStrandMapping);
  this->nrOfStrands -= strandsToDelete.size();

  const std::vector<long int> newSpringMapping =
    pylimer_tools::utils::getMappingForRenumbering(sortedSprings,
                                                   this->nrOfSprings);
  pylimer_tools::utils::renumberWithMapping(this->springIndicesOfStrand,
                                            newSpringMapping);
  this->nrOfSprings -= sortedSprings.size();

  if (removeZeroFunctionalLinksFlag) {
    this->removeZeroFunctionalLinks(u);
  }
  this->recomputeMeanContourLengths();
}

size_t
MEHPNetwork::removeInactiveLinks(Eigen::VectorXd& u,
                                 const double tolerance,
                                 const bool is2D)
{
  const size_t nLinksBefore = this->nrOfLinks;
  const Eigen::ArrayXb springIsToDelete =
    this->findInactiveSprings(u, tolerance, is2D, false);
  const size_t toDeleteCount = static_cast<size_t>(springIsToDelete.count());
  std::vector<size_t> springIndicesToDelete;
  springIndicesToDelete.reserve(toDeleteCount);
  for (size_t i = static_cast<size_t>(springIsToDelete.size()); i > 0; --i) {
    if (springIsToDelete[static_cast<Eigen::Index>(i - 1)]) {
      springIndicesToDelete.push_back(i - 1);
    }
  }
  this->removeSprings(springIndicesToDelete, false, &u);
  this->unlinkBifunctionalLinks(u, false);
  this->removeZeroFunctionalLinks(u);
  return nLinksBefore - this->nrOfLinks;
}

size_t
MEHPNetwork::unlinkBifunctionalLinks(Eigen::VectorXd& u,
                                     const bool entanglementsOnly)
{
  size_t numRemoved = 0;
  std::vector<size_t> removedSprings;
  for (long int linkIdx = static_cast<long int>(this->nrOfLinks) - 1;
       linkIdx >= 0;
       --linkIdx) {
    const size_t lIdx = static_cast<size_t>(linkIdx);
    const Eigen::Index lIdxE = static_cast<Eigen::Index>(linkIdx);
    if (entanglementsOnly && !this->linkIsEntanglement[lIdxE]) {
      continue;
    }
    const std::vector<size_t> springIndices =
      this->getSpringIndicesOfLink(lIdx);
    if (springIndices.size() == 2 && !this->isLoopingSpring(springIndices[0]) &&
        !this->isLoopingSpring(springIndices[1])) {
      const Eigen::Index s0E = static_cast<Eigen::Index>(springIndices[0]);
      const Eigen::Index s1E = static_cast<Eigen::Index>(springIndices[1]);
      if (this->springIsEntanglement[s0E] || this->springIsEntanglement[s1E]) {
        continue;
      }
      if (!this->linkIsEntanglement[lIdxE] &&
          this->strandIndexOfSpring[s0E] == this->strandIndexOfSpring[s1E]) {
        continue;
      }
      mergeTwoSpringsWithoutRemoval(
        *this, u, springIndices[0], springIndices[1], lIdx);
      removedSprings.push_back(springIndices[0]);
      numRemoved += 1;
    }
  }

  if (!removedSprings.empty()) {
    std::sort(
      removedSprings.begin(), removedSprings.end(), std::greater<size_t>());
    this->removeSprings(removedSprings, true, &u);
  }
  return numRemoved;
}

size_t
MEHPNetwork::removeZeroFunctionalLinks(Eigen::VectorXd& u)
{
  return this->removeZeroFunctionalLinks(&u);
}

size_t
MEHPNetwork::removeZeroFunctionalLinks(Eigen::VectorXd* u)
{
  std::vector<size_t> linksToRemove;
  std::vector<size_t> linksWithDirsToRemove;
  for (long int i = static_cast<long int>(this->nrOfLinks) - 1; i >= 0; --i) {
    const size_t linkIdx = static_cast<size_t>(i);
    if (this->strandIndicesOfLink[linkIdx].empty()) {
      linksToRemove.push_back(linkIdx);
      for (long int dir = 2; dir >= 0; --dir) {
        linksWithDirsToRemove.push_back(linkIdx * 3 + static_cast<size_t>(dir));
      }
    }
  }

  if (linksToRemove.empty()) {
    return 0;
  }

  pylimer_tools::utils::removeRows(
    this->linkIsEntanglement, linksToRemove, true, true);
  pylimer_tools::utils::removeRows(
    this->strandIndicesOfLink, linksToRemove, true, true);
  pylimer_tools::utils::removeRows(this->oldAtomIds, linksToRemove, true, true);
  pylimer_tools::utils::removeRows(
    this->oldAtomTypes, linksToRemove, true, true);
  pylimer_tools::utils::removeRows(
    this->nrOfCrosslinkSwapsEndured, linksToRemove, true, true);

  pylimer_tools::utils::removeRows(
    this->coordinates, linksWithDirsToRemove, true, true);
  if (u != nullptr && u->size() > 0) {
    INVALIDARG_EXP_IFN(
      u->size() == this->coordinates.size() +
                     static_cast<Eigen::Index>(linksWithDirsToRemove.size()),
      "Displacement vector u size mismatch in removeZeroFunctionalLinks.");
    pylimer_tools::utils::removeRows(*u, linksWithDirsToRemove, true, true);
  }

  const std::vector<long int> newIndices =
    pylimer_tools::utils::getMappingForRenumbering(linksToRemove,
                                                   this->nrOfLinks);
  safeRenumberWithMapping(this->springIndexA, newIndices);
  safeRenumberWithMapping(this->springIndexB, newIndices);
  pylimer_tools::utils::renumberWithMapping(this->linkIndicesOfStrand,
                                            newIndices);

  std::vector<long int> newIndicesWithDir;
  newIndicesWithDir.reserve(newIndices.size() * 3);
  for (size_t i = 0; i < newIndices.size(); ++i) {
    for (long int dir = 0; dir < 3; ++dir) {
      newIndicesWithDir.push_back(newIndices[i] >= 0 ? (newIndices[i] * 3 + dir)
                                                     : -1);
    }
  }
  safeRenumberWithMapping(this->springCoordinateIndexA, newIndicesWithDir);
  safeRenumberWithMapping(this->springCoordinateIndexB, newIndicesWithDir);

  this->nrOfLinks -= linksToRemove.size();
  this->nrOfNodes =
    static_cast<size_t>((this->linkIsEntanglement == false).count());
  return linksToRemove.size();
}

size_t
MEHPNetwork::simplify(Eigen::VectorXd& u,
                      const StructureSimplificationMode mode,
                      const double tolerance,
                      const bool is2D,
                      const double breakingLength)
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
    "Displacement vector u size must match network nrOfLinks.");

  size_t totalChanges = 0;
  if (breakingLength > 0.0) {
    totalChanges += this->breakTooLongStrands(u, breakingLength, is2D);
  }

  if (mode == StructureSimplificationMode::NO_SIMPLIFICATION) {
    return totalChanges;
  }
  if (mode == StructureSimplificationMode::X1F_X2F_THEN_INACTIVE) {
    RUNTIME_EXP(
      "StructureSimplificationMode::X1F_X2F_THEN_INACTIVE is not implemented; "
      "use INACTIVE_THEN_X2F instead.");
  }

  size_t removedInLoop = 0;
  do {
    removedInLoop = 0;
    if (mode == StructureSimplificationMode::INACTIVE_ONLY ||
        mode == StructureSimplificationMode::INACTIVE_THEN_X2F) {
      const Eigen::ArrayXb inactiveSprings =
        this->findInactiveSprings(u, tolerance, is2D, false);
      const size_t nInactive = static_cast<size_t>(inactiveSprings.count());
      if (nInactive > 0) {
        std::vector<size_t> toDelete;
        toDelete.reserve(nInactive);
        for (size_t i = static_cast<size_t>(inactiveSprings.size()); i > 0;
             --i) {
          if (inactiveSprings[static_cast<Eigen::Index>(i - 1)]) {
            toDelete.push_back(i - 1);
          }
        }
        this->removeSprings(toDelete, true, &u);
        removedInLoop += nInactive;
      }
      removedInLoop += this->removeZeroFunctionalLinks(u);
    }

    if (mode == StructureSimplificationMode::X2F_ONLY ||
        mode == StructureSimplificationMode::INACTIVE_THEN_X2F) {
      removedInLoop += this->unlinkBifunctionalLinks(u, false);
      removedInLoop += this->removeZeroFunctionalLinks(u);
    }

    totalChanges += removedInLoop;
  } while (removedInLoop > 0);

  return totalChanges;
}

void
MEHPNetwork::deformTo(const pylimer_tools::entities::Box& targetBox,
                      Eigen::VectorXd& u)
{
  this->box.adjustCoordinatesTo(this->coordinates, targetBox);
  if (u.size() > 0) {
    INVALIDARG_EXP_IFN(
      u.size() == this->coordinates.size(),
      "Displacement vector u size must match coordinates size in deformTo.");
    this->box.adjustCoordinatesTo(u, targetBox);
  }
  this->box.adjustCoordinatesTo(this->springBoxOffset, targetBox);
  populateBoxMetadata(*this, targetBox);
}

void
MEHPNetwork::deformTo(const pylimer_tools::entities::Box& targetBox)
{
  Eigen::VectorXd emptyU;
  this->deformTo(targetBox, emptyU);
}

// =============================================================================
// Vector, Length & Stress Evaluation
// =============================================================================

Eigen::VectorXd
MEHPNetwork::evaluateSpringVectors(const Eigen::VectorXd& u,
                                   const bool is2D) const
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
    "Displacement vector u size must match network nrOfLinks.");
  if (this->nrOfSprings == 0) {
    return Eigen::VectorXd::Zero(0);
  }

  const Eigen::VectorXd displacedCoords = this->coordinates + u;
  Eigen::VectorXd vectors = (displacedCoords(this->springCoordinateIndexB) -
                             displacedCoords(this->springCoordinateIndexA)) +
                            this->springBoxOffset;

  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(vectors);
  }
  if (is2D) {
    vectors.reshaped(3, static_cast<Eigen::Index>(this->nrOfSprings))
      .row(2)
      .setZero();
  }
  return vectors;
}

Eigen::Vector3d
MEHPNetwork::evaluateSpringVector(const Eigen::VectorXd& u,
                                  const size_t springIdx,
                                  const bool is2D) const
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
    "Displacement vector u size must match network nrOfLinks.");
  INVALIDARG_EXP_IFN(springIdx < this->nrOfSprings,
                     "Spring index out of range in evaluateSpringVector.");

  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  const Eigen::Index idxA = static_cast<Eigen::Index>(this->springIndexA[sE]);
  const Eigen::Index idxB = static_cast<Eigen::Index>(this->springIndexB[sE]);

  Eigen::Vector3d vec =
    (this->coordinates.segment<3>(3 * idxB) + u.segment<3>(3 * idxB)) -
    (this->coordinates.segment<3>(3 * idxA) + u.segment<3>(3 * idxA)) +
    this->springBoxOffset.segment<3>(3 * sE);

  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(vec);
  }
  if (is2D) {
    vec[2] = 0.0;
  }
  return vec;
}

Eigen::VectorXd
MEHPNetwork::evaluateSpringLengths(const Eigen::VectorXd& u,
                                   const bool is2D) const
{
  const Eigen::VectorXd springVecs = this->evaluateSpringVectors(u, is2D);
  if (this->nrOfSprings == 0) {
    return Eigen::VectorXd::Zero(0);
  }
  return springVecs.reshaped(3, static_cast<Eigen::Index>(this->nrOfSprings))
    .colwise()
    .norm()
    .transpose();
}

Eigen::Vector3d
MEHPNetwork::evaluateSpringVectorFrom(const Eigen::VectorXd& u,
                                      const size_t springIdx,
                                      const size_t linkIdx,
                                      const bool is2D) const
{
  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  const Eigen::Vector3d vec = this->evaluateSpringVector(u, springIdx, is2D);
  if (this->springIndexA[sE] == static_cast<int>(linkIdx)) {
    return vec;
  }
  if (this->springIndexB[sE] == static_cast<int>(linkIdx)) {
    return -vec;
  }
  INVALIDARG_EXP_IFN(false,
                     "Link index is not an endpoint of the specified spring.");
  return vec;
}

Eigen::Vector3d
MEHPNetwork::evaluateSpringVectorTo(const Eigen::VectorXd& u,
                                    const size_t springIdx,
                                    const size_t linkIdx,
                                    const bool is2D) const
{
  return -this->evaluateSpringVectorFrom(u, springIdx, linkIdx, is2D);
}

Eigen::VectorXd
MEHPNetwork::evaluateStrandVectors(const Eigen::VectorXd& u,
                                   const bool is2D) const
{
  const Eigen::VectorXd springVecs = this->evaluateSpringVectors(u, is2D);
  Eigen::VectorXd strandVecs =
    Eigen::VectorXd::Zero(3 * static_cast<Eigen::Index>(this->nrOfStrands));
  for (size_t i = 0; i < this->nrOfSprings; ++i) {
    const Eigen::Index sE = static_cast<Eigen::Index>(i);
    const int strandIdxInt = this->strandIndexOfSpring[sE];
    if (strandIdxInt >= 0) {
      strandVecs.segment<3>(3 * static_cast<Eigen::Index>(strandIdxInt)) +=
        springVecs.segment<3>(3 * sE);
    }
  }
  return strandVecs;
}

Eigen::Vector3d
MEHPNetwork::evaluateStrandVector(const Eigen::VectorXd& u,
                                  const size_t strandIdx,
                                  const bool is2D) const
{
  INVALIDARG_EXP_IFN(strandIdx < this->nrOfStrands,
                     "Strand index out of range in evaluateStrandVector.");
  Eigen::Vector3d total = Eigen::Vector3d::Zero();
  for (const size_t spIdx : this->springIndicesOfStrand[strandIdx]) {
    total += this->evaluateSpringVector(u, spIdx, is2D);
  }
  return total;
}

Eigen::VectorXd
MEHPNetwork::evaluateStrandLengths(const Eigen::VectorXd& u,
                                   const bool is2D) const
{
  const Eigen::VectorXd springLens = this->evaluateSpringLengths(u, is2D);
  Eigen::VectorXd strandLens =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(this->nrOfStrands));
  for (size_t i = 0; i < this->nrOfSprings; ++i) {
    const Eigen::Index sE = static_cast<Eigen::Index>(i);
    const int strandIdxInt = this->strandIndexOfSpring[sE];
    if (strandIdxInt >= 0) {
      strandLens[static_cast<Eigen::Index>(strandIdxInt)] += springLens[sE];
    }
  }
  return strandLens;
}

double
MEHPNetwork::evaluateWeightedSpringLength(const Eigen::VectorXd& u,
                                          const size_t springIdx,
                                          const bool is2D) const
{
  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  return this->evaluateSpringVector(u, springIdx, is2D).norm() /
         this->springContourLength[sE];
}

Eigen::Matrix3d
MEHPNetwork::evaluateStressTensor(
  const Eigen::VectorXd& u,
  const double kappa,
  const bool is2D,
  const double oneOverSpringPartitionUpperLimit) const
{
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();
  const Eigen::VectorXd springVectors = this->evaluateSpringVectors(u, is2D);
  const double oneOverVolume = 1.0 / (this->L[0] * this->L[1] * this->L[2]);

  for (size_t i = 0; i < this->nrOfSprings; ++i) {
    const Eigen::Index sE = static_cast<Eigen::Index>(i);
    const Eigen::Vector3d dist = springVectors.segment<3>(3 * sE);
    const int strandIdxInt = this->strandIndexOfSpring[sE];
    const bool isPartial =
      (strandIdxInt >= 0 && !this->springIsEntanglement[sE] &&
       this->springIndicesOfStrand[static_cast<size_t>(strandIdxInt)].size() >
         1);
    const double strandN =
      (strandIdxInt >= 0)
        ? this->strandContourLength[static_cast<Eigen::Index>(strandIdxInt)]
        : this->springContourLength[sE];
    const double invContour = (oneOverSpringPartitionUpperLimit > 0.0)
                                ? computeClampedInverseSpringContourLength(
                                    isPartial,
                                    this->springContourLength[sE],
                                    strandN,
                                    oneOverSpringPartitionUpperLimit)
                                : (1.0 / this->springContourLength[sE]);

    stress.noalias() += (kappa * invContour) * (dist * dist.transpose());
  }

  stress *= oneOverVolume;
  return stress;
}

Eigen::Matrix3d
MEHPNetwork::evaluateStressTensorForLinks(
  const std::vector<size_t>& linkIndices,
  const Eigen::VectorXd& u,
  const double kappa,
  const bool is2D,
  const double oneOverSpringPartitionUpperLimit) const
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
    "Displacement vector u size must match network nrOfLinks.");
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();
  const double halfOverVolume = 0.5 / (this->L[0] * this->L[1] * this->L[2]);

  for (const size_t linkIdx : linkIndices) {
    if (linkIdx >= this->nrOfLinks) {
      continue;
    }
    const std::vector<size_t> springIndices =
      this->getSpringIndicesOfLink(linkIdx);
    for (const size_t spIdx : springIndices) {
      if (this->isLoopingSpring(spIdx)) {
        continue;
      }
      const Eigen::Index sE = static_cast<Eigen::Index>(spIdx);
      const Eigen::Vector3d dist =
        this->evaluateSpringVectorFrom(u, spIdx, linkIdx, is2D);
      const int strandIdxInt = this->strandIndexOfSpring[sE];
      const bool isPartial =
        (strandIdxInt >= 0 && !this->springIsEntanglement[sE] &&
         this->springIndicesOfStrand[static_cast<size_t>(strandIdxInt)].size() >
           1);
      const double strandN =
        (strandIdxInt >= 0)
          ? this->strandContourLength[static_cast<Eigen::Index>(strandIdxInt)]
          : this->springContourLength[sE];
      const double invContour = (oneOverSpringPartitionUpperLimit > 0.0)
                                  ? computeClampedInverseSpringContourLength(
                                      isPartial,
                                      this->springContourLength[sE],
                                      strandN,
                                      oneOverSpringPartitionUpperLimit)
                                  : (1.0 / this->springContourLength[sE]);
      stress.noalias() += (kappa * invContour) * (dist * dist.transpose());
    }
  }

  return halfOverVolume * stress;
}

// =============================================================================
// Topology Queries
// =============================================================================

std::vector<size_t>
MEHPNetwork::getSpringIndicesOfLink(const size_t linkIdx) const
{
  INVALIDARG_EXP_IFN(linkIdx < this->nrOfLinks,
                     "Link index out of range in getSpringIndicesOfLink.");
  std::vector<size_t> result;
  for (const size_t strandIdx : this->strandIndicesOfLink[linkIdx]) {
    for (const size_t spIdx : this->springIndicesOfStrand[strandIdx]) {
      if (this->isPartOfSpring(linkIdx, spIdx)) {
        result.push_back(spIdx);
      }
    }
  }
  return result;
}

std::vector<size_t>
MEHPNetwork::getNeighbourLinkIndices(const size_t linkIdx) const
{
  const std::vector<size_t> springIndices =
    this->getSpringIndicesOfLink(linkIdx);
  std::vector<size_t> neighbors;
  neighbors.reserve(springIndices.size());
  for (const size_t spIdx : springIndices) {
    neighbors.push_back(this->getOtherLinkOfSpring(spIdx, linkIdx));
  }
  return neighbors;
}

std::vector<size_t>
MEHPNetwork::getEntanglementSpringsAtStrand(const size_t strandIdx) const
{
  INVALIDARG_EXP_IFN(
    strandIdx < this->nrOfStrands,
    "Strand index out of range in getEntanglementSpringsAtStrand.");
  std::vector<size_t> entanglementSprings;
  for (const size_t lIdx : this->linkIndicesOfStrand[strandIdx]) {
    if (!this->linkIsEntanglement[static_cast<Eigen::Index>(lIdx)]) {
      continue;
    }
    for (const size_t otherStrandIdx : this->strandIndicesOfLink[lIdx]) {
      if (otherStrandIdx == strandIdx) {
        continue;
      }
      if (this->isEntanglementStrand(otherStrandIdx)) {
        entanglementSprings.push_back(
          this->springIndicesOfStrand[otherStrandIdx][0]);
      }
    }
  }
  return entanglementSprings;
}

size_t
MEHPNetwork::getOtherLinkOfSpring(const size_t springIdx,
                                  const size_t linkIdx) const
{
  INVALIDARG_EXP_IFN(springIdx < this->nrOfSprings,
                     "Spring index out of range in getOtherLinkOfSpring.");
  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  if (this->springIndexA[sE] == static_cast<int>(linkIdx)) {
    return static_cast<size_t>(this->springIndexB[sE]);
  }
  if (this->springIndexB[sE] == static_cast<int>(linkIdx)) {
    return static_cast<size_t>(this->springIndexA[sE]);
  }
  INVALIDARG_EXP_IFN(false,
                     "Link index is not connected to spring in "
                     "getOtherLinkOfSpring.");
  return linkIdx;
}

bool
MEHPNetwork::isPartOfSpring(const size_t linkIdx, const size_t springIdx) const
{
  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  return (this->springIndexA[sE] == static_cast<int>(linkIdx)) ||
         (this->springIndexB[sE] == static_cast<int>(linkIdx));
}

bool
MEHPNetwork::isLoopingSpring(const size_t springIdx) const
{
  const Eigen::Index sE = static_cast<Eigen::Index>(springIdx);
  return this->springIndexA[sE] == this->springIndexB[sE];
}

bool
MEHPNetwork::isLoopingStrand(const size_t strandIdx) const
{
  if (this->linkIndicesOfStrand[strandIdx].empty()) {
    return false;
  }
  return this->linkIndicesOfStrand[strandIdx].front() ==
         this->linkIndicesOfStrand[strandIdx].back();
}

bool
MEHPNetwork::isEntanglementStrand(const size_t strandIdx) const
{
  if (this->springIndicesOfStrand[strandIdx].size() != 1) {
    return false;
  }
  return this->springIsEntanglement[static_cast<Eigen::Index>(
    this->springIndicesOfStrand[strandIdx][0])];
}

bool
MEHPNetwork::distanceIsWithinTolerance(const Eigen::Vector3d& dist,
                                       const double tolerance,
                                       const double contourLength,
                                       const double contourLengthFraction)
{
  return dist.norm() <=
         (tolerance * std::max(contourLengthFraction * contourLength, 1.0));
}

// =============================================================================
// Network Observables & Topological Analysis
// =============================================================================

Eigen::ArrayXb
MEHPNetwork::findActiveSprings(const Eigen::VectorXd& u,
                               const double tolerance,
                               const bool is2D) const
{
  const Eigen::VectorXd springVectors = this->evaluateSpringVectors(u, is2D);
  Eigen::ArrayXb result = Eigen::ArrayXb::Constant(
    static_cast<Eigen::Index>(this->nrOfSprings), false);
  for (size_t i = 0; i < this->nrOfSprings; ++i) {
    const Eigen::Index sE = static_cast<Eigen::Index>(i);
    result[sE] = !distanceIsWithinTolerance(springVectors.segment<3>(3 * sE),
                                            tolerance,
                                            this->springContourLength[sE]);
  }
  return result;
}

Eigen::ArrayXb
MEHPNetwork::findActiveStrands(const Eigen::VectorXd& u,
                               const double tolerance,
                               const bool is2D) const
{
  Eigen::ArrayXb result = Eigen::ArrayXb::Constant(
    static_cast<Eigen::Index>(this->nrOfStrands), false);
  const Eigen::ArrayXb activeSprings =
    this->findActiveSprings(u, tolerance, is2D);

  auto areSameEntanglement = [this](const size_t linkA,
                                    const size_t linkB) -> bool {
    if (linkA == linkB) {
      return true;
    }
    for (const size_t sIdx : this->strandIndicesOfLink[linkA]) {
      if (this->isEntanglementStrand(sIdx)) {
        const std::vector<size_t>& links = this->linkIndicesOfStrand[sIdx];
        if (links.size() == 2 && (links[0] == linkB || links[1] == linkB)) {
          return true;
        }
      }
    }
    return false;
  };

  for (size_t strandIdx = 0; strandIdx < this->nrOfStrands; ++strandIdx) {
    const std::vector<size_t>& strandSprings =
      this->springIndicesOfStrand[strandIdx];
    const std::vector<size_t>& strandLinks =
      this->linkIndicesOfStrand[strandIdx];
    if (strandSprings.empty()) {
      continue;
    }
    if (strandSprings.size() == 1) {
      result[static_cast<Eigen::Index>(strandIdx)] =
        activeSprings[static_cast<Eigen::Index>(strandSprings[0])];
      continue;
    }

    const bool anyActive =
      std::any_of(strandSprings.begin(),
                  strandSprings.end(),
                  [&activeSprings](const size_t spIdx) {
                    return activeSprings[static_cast<Eigen::Index>(spIdx)];
                  });
    if (!anyActive) {
      continue;
    }

    bool isActive0 =
      activeSprings[static_cast<Eigen::Index>(strandSprings.front())];
    bool isActiveN =
      activeSprings[static_cast<Eigen::Index>(strandSprings.back())];

    if (!isActive0 && strandLinks.size() >= 2) {
      const size_t crossLink0 = strandLinks[0];
      const size_t entLink0 = strandLinks[1];
      if (!this->linkIsEntanglement[static_cast<Eigen::Index>(crossLink0)] &&
          this->linkIsEntanglement[static_cast<Eigen::Index>(entLink0)]) {
        for (const size_t sOther : this->strandIndicesOfLink[crossLink0]) {
          if (sOther == strandIdx ||
              this->linkIndicesOfStrand[sOther].size() < 2) {
            continue;
          }
          if ((this->linkIndicesOfStrand[sOther][0] == crossLink0 &&
               areSameEntanglement(this->linkIndicesOfStrand[sOther][1],
                                   entLink0)) ||
              (this->linkIndicesOfStrand[sOther].back() == crossLink0 &&
               areSameEntanglement(
                 this->linkIndicesOfStrand
                   [sOther][this->linkIndicesOfStrand[sOther].size() - 2],
                 entLink0))) {
            isActive0 = true;
          }
        }
      }
    }

    if (!isActiveN && strandLinks.size() >= 2) {
      const size_t crossLinkN = strandLinks.back();
      const size_t entLinkN = strandLinks[strandLinks.size() - 2];
      if (!this->linkIsEntanglement[static_cast<Eigen::Index>(crossLinkN)] &&
          this->linkIsEntanglement[static_cast<Eigen::Index>(entLinkN)]) {
        for (const size_t sOther : this->strandIndicesOfLink[crossLinkN]) {
          if (sOther == strandIdx ||
              this->linkIndicesOfStrand[sOther].size() < 2) {
            continue;
          }
          if ((this->linkIndicesOfStrand[sOther][0] == crossLinkN &&
               areSameEntanglement(this->linkIndicesOfStrand[sOther][1],
                                   entLinkN)) ||
              (this->linkIndicesOfStrand[sOther].back() == crossLinkN &&
               areSameEntanglement(
                 this->linkIndicesOfStrand
                   [sOther][this->linkIndicesOfStrand[sOther].size() - 2],
                 entLinkN))) {
            isActiveN = true;
          }
        }
      }
    }

    if (!strandLinks.empty() && strandLinks.front() != strandLinks.back()) {
      if (this->strandIndicesOfLink[strandLinks.front()].size() == 1) {
        isActive0 = false;
      }
      if (this->strandIndicesOfLink[strandLinks.back()].size() == 1) {
        isActiveN = false;
      }
    }

    result[static_cast<Eigen::Index>(strandIdx)] = isActive0 && isActiveN;
  }

  return result;
}

Eigen::ArrayXb
MEHPNetwork::findActiveStrandsInDir(const Eigen::VectorXd& u,
                                    const int dir,
                                    const double tolerance,
                                    const bool is2D) const
{
  INVALIDARG_EXP_IFN(dir >= 0 && dir < 3, "Direction must be 0, 1, or 2.");
  const Eigen::VectorXd springVectors = this->evaluateSpringVectors(u, is2D);
  Eigen::ArrayXb result = Eigen::ArrayXb::Constant(
    static_cast<Eigen::Index>(this->nrOfStrands), false);
  for (size_t i = 0; i < this->nrOfSprings; ++i) {
    const Eigen::Index sE = static_cast<Eigen::Index>(i);
    const int strandIdxInt = this->strandIndexOfSpring[sE];
    if (strandIdxInt < 0) {
      continue;
    }
    const Eigen::Index strandE = static_cast<Eigen::Index>(strandIdxInt);
    const Eigen::Vector3d dirVec(
      springVectors[3 * sE + static_cast<Eigen::Index>(dir)], 0.0, 0.0);
    result[strandE] =
      result[strandE] || !distanceIsWithinTolerance(
                           dirVec, tolerance, this->springContourLength[sE]);
  }
  return result;
}

Eigen::ArrayXb
MEHPNetwork::findActiveLinks(const Eigen::VectorXd& u,
                             const double tolerance,
                             const bool is2D) const
{
  const Eigen::ArrayXb activeStrands =
    this->findActiveStrands(u, tolerance, is2D);
  Eigen::ArrayXb activeLinks =
    Eigen::ArrayXb::Constant(static_cast<Eigen::Index>(this->nrOfLinks), false);
  for (size_t i = 0; i < this->nrOfStrands; ++i) {
    if (activeStrands[static_cast<Eigen::Index>(i)] &&
        !this->linkIndicesOfStrand[i].empty()) {
      activeLinks[static_cast<Eigen::Index>(
        this->linkIndicesOfStrand[i].front())] = true;
      activeLinks[static_cast<Eigen::Index>(
        this->linkIndicesOfStrand[i].back())] = true;
    }
  }
  return activeLinks;
}

Eigen::ArrayXb
MEHPNetwork::findActiveNodes(const Eigen::VectorXd& u,
                             const double tolerance,
                             const bool is2D) const
{
  const Eigen::ArrayXb activeLinks = this->findActiveLinks(u, tolerance, is2D);
  Eigen::ArrayXb activeNodes =
    Eigen::ArrayXb::Constant(static_cast<Eigen::Index>(this->nrOfNodes), false);
  size_t nodeIdx = 0;
  for (size_t i = 0; i < this->nrOfLinks; ++i) {
    const Eigen::Index lIdx = static_cast<Eigen::Index>(i);
    if (!this->linkIsEntanglement[lIdx]) {
      activeNodes[static_cast<Eigen::Index>(nodeIdx++)] = activeLinks[lIdx];
    }
  }
  return activeNodes;
}

std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
MEHPNetwork::findSpringsAndLinksClusteredToActive(const Eigen::VectorXd& u,
                                                  const double tolerance,
                                                  const bool is2D) const
{
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(this->nrOfLinks * 3),
    "Displacement vector u size must match network nrOfLinks.");

  Eigen::ArrayXb strandIsActive = this->findActiveStrands(u, tolerance, is2D);
  Eigen::ArrayXb springIsActive = Eigen::ArrayXb::Constant(
    static_cast<Eigen::Index>(this->nrOfSprings), false);
  Eigen::ArrayXb linkIsActive =
    Eigen::ArrayXb::Constant(static_cast<Eigen::Index>(this->nrOfLinks), false);
  std::vector<size_t> linkQueue;
  linkQueue.reserve(this->nrOfLinks);

  for (size_t strandIdx = 0; strandIdx < this->nrOfStrands; ++strandIdx) {
    if (strandIsActive[static_cast<Eigen::Index>(strandIdx)]) {
      for (const size_t spIdx : this->springIndicesOfStrand[strandIdx]) {
        springIsActive[static_cast<Eigen::Index>(spIdx)] = true;
      }
    }
  }

  for (size_t linkIdx = 0; linkIdx < this->nrOfLinks; ++linkIdx) {
    for (const size_t strandIdx : this->strandIndicesOfLink[linkIdx]) {
      if (strandIsActive[static_cast<Eigen::Index>(strandIdx)]) {
        linkIsActive[static_cast<Eigen::Index>(linkIdx)] = true;
        linkQueue.push_back(linkIdx);
        break;
      }
    }
  }

  size_t head = 0;
  while (head < linkQueue.size()) {
    const size_t linkIdx = linkQueue[head++];
    for (const size_t strandIdx : this->strandIndicesOfLink[linkIdx]) {
      if (!strandIsActive[static_cast<Eigen::Index>(strandIdx)]) {
        strandIsActive[static_cast<Eigen::Index>(strandIdx)] = true;
        for (const size_t spIdx : this->springIndicesOfStrand[strandIdx]) {
          springIsActive[static_cast<Eigen::Index>(spIdx)] = true;
        }
      }
      for (const size_t otherLink : this->linkIndicesOfStrand[strandIdx]) {
        if (!linkIsActive[static_cast<Eigen::Index>(otherLink)]) {
          linkIsActive[static_cast<Eigen::Index>(otherLink)] = true;
          linkQueue.push_back(otherLink);
        }
      }
    }
  }

  return std::make_pair(springIsActive, linkIsActive);
}

std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
MEHPNetwork::findStrandsAndNodesClusteredToActive(const Eigen::VectorXd& u,
                                                  const double tolerance,
                                                  const bool is2D) const
{
  const std::pair<Eigen::ArrayXb, Eigen::ArrayXb> clustered =
    this->findSpringsAndLinksClusteredToActive(u, tolerance, is2D);

  Eigen::ArrayXb strandIsActive = Eigen::ArrayXb::Constant(
    static_cast<Eigen::Index>(this->nrOfStrands), false);
  Eigen::ArrayXb nodeIsActive =
    Eigen::ArrayXb::Constant(static_cast<Eigen::Index>(this->nrOfNodes), false);

  for (size_t i = 0; i < this->nrOfSprings; ++i) {
    const Eigen::Index sE = static_cast<Eigen::Index>(i);
    if (clustered.first[sE] && this->strandIndexOfSpring[sE] >= 0) {
      strandIsActive[this->strandIndexOfSpring[sE]] = true;
    }
  }

  size_t xlinkIdx = 0;
  for (size_t i = 0; i < this->nrOfLinks; ++i) {
    const Eigen::Index lIdx = static_cast<Eigen::Index>(i);
    if (this->linkIsEntanglement[lIdx]) {
      continue;
    }
    if (clustered.second[lIdx]) {
      nodeIsActive[static_cast<Eigen::Index>(xlinkIdx)] = true;
    }
    xlinkIdx += 1;
  }

  return std::make_pair(strandIsActive, nodeIsActive);
}

std::vector<int>
MEHPNetwork::getIndicesOfActiveNodes(const Eigen::VectorXd& u,
                                     const double tolerance,
                                     const bool is2D) const
{
  std::vector<int> results;
  results.reserve(this->nrOfNodes);
  const Eigen::ArrayXb strandIsActive =
    this->findActiveStrands(u, tolerance, is2D);

  for (size_t i = 0; i < this->nrOfLinks; ++i) {
    const Eigen::Index lIdx = static_cast<Eigen::Index>(i);
    if (this->linkIsEntanglement[lIdx]) {
      continue;
    }
    for (const size_t strandIdx : this->strandIndicesOfLink[i]) {
      if (strandIsActive[static_cast<Eigen::Index>(strandIdx)]) {
        results.push_back(static_cast<int>(i));
        break;
      }
    }
  }
  return results;
}

double
MEHPNetwork::computeActiveWeightFraction(const Eigen::VectorXd& u,
                                         const double tolerance,
                                         const bool is2D) const
{
  if (this->nrOfStrands < 1) {
    return 0.0;
  }
  const Eigen::ArrayXb activeStrands =
    this->findActiveStrands(u, tolerance, is2D);
  if (activeStrands.count() == 0) {
    return 0.0;
  }

  Eigen::ArrayXb activeLinks =
    Eigen::ArrayXb::Constant(static_cast<Eigen::Index>(this->nrOfLinks), false);
  Eigen::ArrayXb activeSprings = Eigen::ArrayXb::Constant(
    static_cast<Eigen::Index>(this->nrOfSprings), false);

  for (size_t strandIdx = 0; strandIdx < this->nrOfStrands; ++strandIdx) {
    if (!activeStrands[static_cast<Eigen::Index>(strandIdx)]) {
      continue;
    }
    for (const size_t lIdx : this->linkIndicesOfStrand[strandIdx]) {
      activeLinks[static_cast<Eigen::Index>(lIdx)] = true;
    }
    for (const size_t spIdx : this->springIndicesOfStrand[strandIdx]) {
      activeSprings[static_cast<Eigen::Index>(spIdx)] = true;
    }
  }

  const double nActiveAtomsFromSprings =
    ((this->springContourLength.array() - 1.0) * activeSprings.cast<double>())
      .sum();
  const double nActiveCrossLinks = static_cast<double>(
    (activeLinks * (this->linkIsEntanglement == false)).count());
  const double nActiveEntanglementLinks =
    static_cast<double>((activeLinks * this->linkIsEntanglement).count());
  const bool hasEntSprings =
    this->springIsEntanglement.size() > 0 && this->springIsEntanglement.any();

  const double nAtoms = (this->assumeComplete || this->universeNrOfAtoms == 0)
                          ? this->inferNrOfAtomsFromNetwork()
                          : static_cast<double>(this->universeNrOfAtoms);
  if (nAtoms <= 0.0) {
    return 0.0;
  }

  return (nActiveAtomsFromSprings +
          (hasEntSprings ? 1.0 : 2.0) * nActiveEntanglementLinks +
          nActiveCrossLinks) /
         nAtoms;
}

double
MEHPNetwork::countActiveClusteredAtoms(const Eigen::VectorXd& u,
                                       const double tolerance,
                                       const bool is2D) const
{
  INVALIDARG_EXP_IFN(
    static_cast<Eigen::Index>(this->nrOfLinks * 3) == u.size(),
    "Displacement vector u size must match network nrOfLinks.");
  if (this->nrOfStrands < 1) {
    return 0.0;
  }

  if (this->assumeComplete || this->clusterSizes.empty()) {
    const std::pair<Eigen::ArrayXb, Eigen::ArrayXb> clustered =
      this->findSpringsAndLinksClusteredToActive(u, tolerance, is2D);
    const Eigen::ArrayXb& clusteredSprings = clustered.first;
    const Eigen::ArrayXb& clusteredLinks = clustered.second;

    const double nAtomsInActiveStrands =
      (clusteredSprings.cast<double>() *
       (this->springContourLength.array() - 1.0))
        .sum();
    const double nActiveCrossLinks =
      (clusteredLinks.cast<double>() *
       (this->linkIsEntanglement == false).cast<double>())
        .sum();
    const double nActiveEntanglementLinks =
      (clusteredLinks.cast<double>() * this->linkIsEntanglement.cast<double>())
        .sum();
    const bool hasEntSprings =
      this->springIsEntanglement.size() > 0 && this->springIsEntanglement.any();
    return nAtomsInActiveStrands + nActiveCrossLinks +
           (hasEntSprings ? 1.0 : 2.0) * nActiveEntanglementLinks;
  }

  RUNTIME_EXP_IFN(this->universeNrOfAtoms >= this->nrOfNodes,
                  "Expected universeNrOfAtoms >= nrOfNodes, got " +
                    std::to_string(this->universeNrOfAtoms) + " < " +
                    std::to_string(this->nrOfNodes) + ".");

  std::vector<bool> clusterIsActive(this->clusterSizes.size(), false);
  const std::vector<int> activeNodeIndices =
    this->getIndicesOfActiveNodes(u, tolerance, is2D);

  for (const int nodeIdx : activeNodeIndices) {
    const long int atomId =
      this->oldAtomIds[static_cast<Eigen::Index>(nodeIdx)];
    const std::unordered_map<long int, size_t>::const_iterator it =
      this->atomIdToClusterIdx.find(atomId);
    if (it != this->atomIdToClusterIdx.end()) {
      clusterIsActive[it->second] = true;
    }
  }

  double nClusteredAtoms = 0.0;
  for (size_t i = 0; i < this->clusterSizes.size(); ++i) {
    if (clusterIsActive[i]) {
      nClusteredAtoms += static_cast<double>(this->clusterSizes[i]);
    }
  }
  return nClusteredAtoms;
}

double
MEHPNetwork::computeSolubleWeightFraction(const Eigen::VectorXd& u,
                                          const double tolerance,
                                          const bool is2D) const
{
  const double nAtoms = (this->assumeComplete || this->universeNrOfAtoms == 0)
                          ? this->inferNrOfAtomsFromNetwork()
                          : static_cast<double>(this->universeNrOfAtoms);
  RUNTIME_EXP_IFN(nAtoms > 0.0,
                  "Expected number of atoms to be positive, got " +
                    std::to_string(nAtoms) + ".");
  const double nActiveClusteredAtoms =
    this->countActiveClusteredAtoms(u, tolerance, is2D);
  return 1.0 - (nActiveClusteredAtoms / nAtoms);
}

double
MEHPNetwork::computeDanglingWeightFraction(const Eigen::VectorXd& u,
                                           const double tolerance,
                                           const bool is2D) const
{
  const double activeFrac =
    this->computeActiveWeightFraction(u, tolerance, is2D);
  const double solubleFrac =
    this->computeSolubleWeightFraction(u, tolerance, is2D);
  RUNTIME_EXP_IFN(APPROX_WITHIN(activeFrac, 0.0, 1.0, 1e-6),
                  "Active weight fraction out of [0, 1]: " +
                    std::to_string(activeFrac) + ".");
  RUNTIME_EXP_IFN(APPROX_WITHIN(solubleFrac, 0.0, 1.0, 1e-6),
                  "Soluble weight fraction out of [0, 1]: " +
                    std::to_string(solubleFrac) + ".");
  return std::max(0.0, 1.0 - activeFrac - solubleFrac);
}

std::unordered_map<long int, int>
MEHPNetwork::getEffectiveFunctionalityOfAtoms(const Eigen::VectorXd& u,
                                              const double tolerance,
                                              const bool is2D) const
{
  std::unordered_map<long int, int> results;
  results.reserve(this->nrOfNodes);

  const Eigen::ArrayXb activeStrands =
    this->findActiveStrands(u, tolerance, is2D);
  std::vector<int> activeStrandsPerLink(this->nrOfLinks, 0);
  for (size_t s = 0; s < this->nrOfStrands; ++s) {
    if (!activeStrands[static_cast<Eigen::Index>(s)] ||
        this->linkIndicesOfStrand[s].empty() || this->isEntanglementStrand(s)) {
      continue;
    }
    activeStrandsPerLink[this->linkIndicesOfStrand[s].front()] += 1;
    activeStrandsPerLink[this->linkIndicesOfStrand[s].back()] += 1;
  }

  for (size_t i = 0; i < this->nrOfLinks; ++i) {
    const Eigen::Index lIdx = static_cast<Eigen::Index>(i);
    if (this->linkIsEntanglement[lIdx]) {
      continue;
    }
    const long int key = (this->oldAtomIds[lIdx] >= 0)
                           ? static_cast<long int>(this->oldAtomIds[lIdx])
                           : static_cast<long int>(i);
    results.emplace(key, activeStrandsPerLink[i]);
  }
  return results;
}

Eigen::VectorXd
MEHPNetwork::getGammaFactors(const Eigen::VectorXd& u,
                             double b02,
                             const bool is2D) const
{
  if (b02 < 0.0) {
    b02 = (this->defaultBondLength > 0.0)
            ? (this->defaultBondLength * this->defaultBondLength)
            : this->defaultR02;
  }
  INVALIDARG_EXP_IFN(b02 > 0.0, "b02 must be strictly positive.");

  if (this->nrOfSprings == 0) {
    return Eigen::VectorXd::Zero(0);
  }
  const Eigen::VectorXd springVectors = this->evaluateSpringVectors(u, is2D);
  const Eigen::Index nSprings = static_cast<Eigen::Index>(this->nrOfSprings);
  const double commonDenominator = 1.0 / b02;

  Eigen::VectorXd gammaFactors =
    (springVectors.reshaped(3, nSprings)
       .colwise()
       .squaredNorm()
       .transpose()
       .array() *
     commonDenominator / this->springContourLength.array())
      .matrix();

  if (!gammaFactors.allFinite()) {
    for (Eigen::Index i = 0; i < nSprings; ++i) {
      RUNTIME_EXP_IFN(std::isfinite(gammaFactors[i]),
                      "Non-finite gamma factor for spring " +
                        std::to_string(i) + ".");
    }
  }
  return gammaFactors;
}

Eigen::VectorXd
MEHPNetwork::getGammaFactorsInDir(const Eigen::VectorXd& u,
                                  double b02,
                                  const int dir,
                                  const bool is2D) const
{
  INVALIDARG_EXP_IFN(dir >= 0 && dir <= 2, "Direction must be 0, 1, or 2.");
  if (b02 < 0.0) {
    b02 = (this->defaultBondLength > 0.0)
            ? (this->defaultBondLength * this->defaultBondLength)
            : this->defaultR02;
  }
  INVALIDARG_EXP_IFN(b02 > 0.0, "b02 must be strictly positive.");

  if (this->nrOfSprings == 0) {
    return Eigen::VectorXd::Zero(0);
  }

  const Eigen::VectorXd springVectors = this->evaluateSpringVectors(u, is2D);
  const Eigen::Index nSprings = static_cast<Eigen::Index>(this->nrOfSprings);
  const double commonDenominator = 1.0 / b02;

  return (springVectors.reshaped(3, nSprings)
            .row(static_cast<Eigen::Index>(dir))
            .transpose()
            .array()
            .square() *
          commonDenominator / this->springContourLength.array())
    .matrix();
}

double
MEHPNetwork::getGammaFactor(const Eigen::VectorXd& u,
                            double b02,
                            const int nrOfChains,
                            const bool is2D) const
{
  if (this->nrOfSprings == 0) {
    return 0.0;
  }
  const Eigen::VectorXd gammaFactors = this->getGammaFactors(u, b02, is2D);
  const double denom = (nrOfChains < 1) ? static_cast<double>(this->nrOfStrands)
                                        : static_cast<double>(nrOfChains);
  return (denom > 0.0) ? (gammaFactors.sum() / denom) : 0.0;
}

pylimer_tools::entities::Universe
MEHPNetwork::getCrosslinkerVerse(const Eigen::VectorXd& u) const
{
  if (u.size() > 0) {
    INVALIDARG_EXP_IFN(
      u.size() == static_cast<Eigen::Index>(3 * this->nrOfLinks),
      "Displacement vector u size must match network nrOfLinks.");
  }

  pylimer_tools::entities::Universe xlinkUniverse(this->box);
  std::vector<long int> ids;
  std::vector<int> types;
  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<double> zs;
  ids.reserve(this->nrOfNodes);
  types.reserve(this->nrOfNodes);
  xs.reserve(this->nrOfNodes);
  ys.reserve(this->nrOfNodes);
  zs.reserve(this->nrOfNodes);

  std::vector<long int> linkToAtomId(this->nrOfLinks, -1);
  for (size_t i = 0; i < this->nrOfLinks; ++i) {
    const Eigen::Index lIdx = static_cast<Eigen::Index>(i);
    if (this->linkIsEntanglement[lIdx]) {
      continue;
    }
    const double dx = (u.size() > 0) ? u[3 * lIdx + 0] : 0.0;
    const double dy = (u.size() > 0) ? u[3 * lIdx + 1] : 0.0;
    const double dz = (u.size() > 0) ? u[3 * lIdx + 2] : 0.0;
    xs.push_back(this->coordinates[3 * lIdx + 0] + dx);
    ys.push_back(this->coordinates[3 * lIdx + 1] + dy);
    zs.push_back(this->coordinates[3 * lIdx + 2] + dz);

    const long int atomId = (this->oldAtomIds[lIdx] >= 0)
                              ? static_cast<long int>(this->oldAtomIds[lIdx])
                              : static_cast<long int>(i);
    linkToAtomId[i] = atomId;
    ids.push_back(atomId);
    types.push_back(this->oldAtomTypes[lIdx]);
  }

  const std::vector<int> zeros =
    pylimer_tools::utils::initializeWithValue(ids.size(), 0);
  xlinkUniverse.addAtoms(ids, types, xs, ys, zs, zeros, zeros, zeros);

  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  bondFrom.reserve(this->nrOfStrands);
  bondTo.reserve(this->nrOfStrands);
  for (size_t s = 0; s < this->nrOfStrands; ++s) {
    if (this->linkIndicesOfStrand[s].size() < 2 ||
        this->isEntanglementStrand(s)) {
      continue;
    }
    const size_t linkFrom = this->linkIndicesOfStrand[s].front();
    const size_t linkTo = this->linkIndicesOfStrand[s].back();
    if (this->linkIsEntanglement[static_cast<Eigen::Index>(linkFrom)] ||
        this->linkIsEntanglement[static_cast<Eigen::Index>(linkTo)]) {
      continue;
    }
    bondFrom.push_back(linkToAtomId[linkFrom]);
    bondTo.push_back(linkToAtomId[linkTo]);
  }

  xlinkUniverse.addBonds(
    bondFrom.size(),
    bondFrom,
    bondTo,
    pylimer_tools::utils::initializeWithValue(bondFrom.size(), 1),
    false,
    false);
  return xlinkUniverse;
}

size_t
MEHPNetwork::getNumIntraChainSlipLinks() const
{
  size_t count = 0;
  for (size_t s = 0; s < this->nrOfStrands; ++s) {
    if (this->isEntanglementStrand(s) ||
        this->linkIndicesOfStrand[s].size() < 4) {
      continue;
    }
    const std::vector<size_t>& links = this->linkIndicesOfStrand[s];
    std::vector<size_t> interior(links.begin() + 1, links.end() - 1);
    std::sort(interior.begin(), interior.end());
    for (size_t k = 1; k < interior.size(); ++k) {
      if (interior[k] == interior[k - 1]) {
        count += 1;
        while (k + 1 < interior.size() && interior[k + 1] == interior[k]) {
          ++k;
        }
      }
    }
  }
  return count;
}

// =============================================================================
// Validation
// =============================================================================

void
MEHPNetwork::validate(const Eigen::VectorXd& u) const
{
  for (size_t dir = 0; dir < 3; ++dir) {
    RUNTIME_EXP_IFN(std::isfinite(this->L[dir]) && this->L[dir] > 0.0,
                    "Box dimension L[" + std::to_string(dir) +
                      "] must be finite and positive.");
    RUNTIME_EXP_IFN(
      APPROX_EQUAL(this->boxHalfs[dir], 0.5 * this->L[dir], 1e-12),
      "Expected boxHalfs to equal 0.5 * L.");
    RUNTIME_EXP_IFN(
      APPROX_EQUAL(this->L[dir], this->box.getL(static_cast<int>(dir)), 1e-12),
      "MEHPNetwork.L must stay synchronized with MEHPNetwork.box.");
  }
  RUNTIME_EXP_IFN(std::isfinite(this->vol) && this->vol > 0.0,
                  "Box volume must be finite and positive.");

  const Eigen::Index nLinksE = static_cast<Eigen::Index>(this->nrOfLinks);
  const Eigen::Index nStrandsE = static_cast<Eigen::Index>(this->nrOfStrands);
  const Eigen::Index nSpringsE = static_cast<Eigen::Index>(this->nrOfSprings);

  RUNTIME_EXP_IFN(this->coordinates.size() == 3 * nLinksE,
                  "Invalid size of coordinates.");
  if (u.size() > 0) {
    RUNTIME_EXP_IFN(u.size() == 3 * nLinksE,
                    "Invalid size of displacements vector u.");
  }
  RUNTIME_EXP_IFN(this->linkIsEntanglement.size() == nLinksE,
                  "Invalid size of linkIsEntanglement.");
  RUNTIME_EXP_IFN(
    static_cast<size_t>((this->linkIsEntanglement == false).count()) ==
      this->nrOfNodes,
    "nrOfNodes must match count of links with !linkIsEntanglement.");
  RUNTIME_EXP_IFN(this->oldAtomIds.size() == nLinksE,
                  "Invalid size of oldAtomIds.");
  RUNTIME_EXP_IFN(this->oldAtomTypes.size() == nLinksE,
                  "Invalid size of oldAtomTypes.");
  RUNTIME_EXP_IFN(this->strandIndicesOfLink.size() == this->nrOfLinks,
                  "Invalid size of strandIndicesOfLink.");
  RUNTIME_EXP_IFN(this->nrOfCrosslinkSwapsEndured.size() == nLinksE,
                  "nrOfCrosslinkSwapsEndured must be sized to nrOfLinks.");

  RUNTIME_EXP_IFN(this->strandContourLength.size() == nStrandsE,
                  "Invalid size of strandContourLength.");
  RUNTIME_EXP_IFN(this->linkIndicesOfStrand.size() == this->nrOfStrands,
                  "Invalid size of linkIndicesOfStrand.");
  RUNTIME_EXP_IFN(this->springIndicesOfStrand.size() == this->nrOfStrands,
                  "Invalid size of springIndicesOfStrand.");

  RUNTIME_EXP_IFN(this->springIndexA.size() == nSpringsE,
                  "Invalid size of springIndexA.");
  RUNTIME_EXP_IFN(this->springIndexB.size() == nSpringsE,
                  "Invalid size of springIndexB.");
  RUNTIME_EXP_IFN(this->springCoordinateIndexA.size() == 3 * nSpringsE,
                  "Invalid size of springCoordinateIndexA.");
  RUNTIME_EXP_IFN(this->springCoordinateIndexB.size() == 3 * nSpringsE,
                  "Invalid size of springCoordinateIndexB.");
  RUNTIME_EXP_IFN(this->springBoxOffset.size() == 3 * nSpringsE,
                  "Invalid size of springBoxOffset.");
  RUNTIME_EXP_IFN(this->springContourLength.size() == nSpringsE,
                  "Invalid size of springContourLength.");
  RUNTIME_EXP_IFN(this->springPartitions.size() == nSpringsE,
                  "Invalid size of springPartitions.");
  RUNTIME_EXP_IFN(this->strandIndexOfSpring.size() == nSpringsE,
                  "Invalid size of strandIndexOfSpring.");
  RUNTIME_EXP_IFN(this->springIsEntanglement.size() == nSpringsE,
                  "Invalid size of springIsEntanglement.");

  for (size_t linkIdx = 0; linkIdx < this->nrOfLinks; ++linkIdx) {
    const std::vector<size_t>& linkStrands = this->strandIndicesOfLink[linkIdx];
    for (size_t j = 0; j < linkStrands.size(); ++j) {
      const size_t strandIdx = linkStrands[j];
      RUNTIME_EXP_IFN(strandIdx < this->nrOfStrands,
                      "Link references out-of-range strand index.");
      RUNTIME_EXP_IFN(
        std::find(linkStrands.begin() + static_cast<std::ptrdiff_t>(j + 1),
                  linkStrands.end(),
                  strandIdx) == linkStrands.end(),
        "Duplicate strand index in strandIndicesOfLink.");
      RUNTIME_EXP_IFN(std::find(this->linkIndicesOfStrand[strandIdx].begin(),
                                this->linkIndicesOfStrand[strandIdx].end(),
                                linkIdx) !=
                        this->linkIndicesOfStrand[strandIdx].end(),
                      "Strand <-> link mapping inconsistency.");
    }
  }

  for (size_t strandIdx = 0; strandIdx < this->nrOfStrands; ++strandIdx) {
    const std::vector<size_t>& links = this->linkIndicesOfStrand[strandIdx];
    const std::vector<size_t>& springs = this->springIndicesOfStrand[strandIdx];
    RUNTIME_EXP_IFN(
      !links.empty() && links.size() == springs.size() + 1,
      "Each strand must have links.size() == springs.size() + 1.");

    const bool isEntStrand = this->isEntanglementStrand(strandIdx);
    if (isEntStrand) {
      const Eigen::Index sp0 = static_cast<Eigen::Index>(springs[0]);
      RUNTIME_EXP_IFN(links.size() == 2,
                      "Entanglement strand must connect exactly two links.");
      RUNTIME_EXP_IFN(
        this->linkIsEntanglement[static_cast<Eigen::Index>(links[0])] &&
          this->linkIsEntanglement[static_cast<Eigen::Index>(links[1])],
        "Entanglement strand endpoints must be entanglement links.");
      RUNTIME_EXP_IFN(
        APPROX_EQUAL(this->springContourLength[sp0], 1.0, 1e-9) &&
          APPROX_EQUAL(
            this->strandContourLength[static_cast<Eigen::Index>(strandIdx)],
            1.0,
            1e-9),
        "Entanglement strand and spring must have contour length 1.0.");
    }

    double sumSubContour = 0.0;
    double sumPartitions = 0.0;
    const double totalN =
      this->strandContourLength[static_cast<Eigen::Index>(strandIdx)];

    for (size_t m = 0; m < springs.size(); ++m) {
      const size_t spIdx = springs[m];
      RUNTIME_EXP_IFN(spIdx < this->nrOfSprings,
                      "Strand references out-of-range spring index.");
      const Eigen::Index spE = static_cast<Eigen::Index>(spIdx);
      RUNTIME_EXP_IFN(this->strandIndexOfSpring[spE] ==
                        static_cast<int>(strandIdx),
                      "Spring <-> strand index mismatch.");
      RUNTIME_EXP_IFN(this->springIndexA[spE] == static_cast<int>(links[m]) &&
                        this->springIndexB[spE] ==
                          static_cast<int>(links[m + 1]),
                      "Spring endpoint ordering along strand is inconsistent.");
      sumSubContour += this->springContourLength[spE];
      sumPartitions += this->springPartitions[spE];
      RUNTIME_EXP_IFN(APPROX_EQUAL(this->springContourLength[spE],
                                   this->springPartitions[spE] * totalN,
                                   1e-6),
                      "springContourLength must equal springPartitions * "
                      "strandContourLength.");
    }
    RUNTIME_EXP_IFN(
      APPROX_EQUAL(sumSubContour, totalN, 1e-6),
      "Sub-spring contour lengths must sum to strandContourLength.");
    RUNTIME_EXP_IFN(APPROX_EQUAL(sumPartitions, 1.0, 1e-6),
                    "Sub-spring partitions along strand must sum to 1.0.");

    if (!isEntStrand) {
      for (size_t j = 0; j < links.size(); ++j) {
        const Eigen::Index lIdx = static_cast<Eigen::Index>(links[j]);
        const bool expectedEnt = (j != 0) && (j + 1 != links.size());
        RUNTIME_EXP_IFN(this->linkIsEntanglement[lIdx] == expectedEnt,
                        "Polymer strand endpoints must be non-entanglement "
                        "nodes and interior "
                        "links must be entanglements.");
        RUNTIME_EXP_IFN(std::find(this->strandIndicesOfLink[links[j]].begin(),
                                  this->strandIndicesOfLink[links[j]].end(),
                                  strandIdx) !=
                          this->strandIndicesOfLink[links[j]].end(),
                        "Strand missing from strandIndicesOfLink.");
      }
    }
  }

  for (size_t i = 0; i < this->nrOfSprings; ++i) {
    const Eigen::Index sE = static_cast<Eigen::Index>(i);
    RUNTIME_EXP_IFN(std::isfinite(this->springContourLength[sE]) &&
                      this->springContourLength[sE] > 0.0,
                    "Spring contour length must be finite and positive.");
    RUNTIME_EXP_IFN(std::isfinite(this->springPartitions[sE]) &&
                      this->springPartitions[sE] > 0.0,
                    "Spring partition must be finite and positive.");
    const int a = this->springIndexA[sE];
    const int b = this->springIndexB[sE];
    if (a < 0 || b < 0) {
      RUNTIME_EXP_IFN(a == -1 && b == -1 && this->strandIndexOfSpring[sE] == -1,
                      "Disabled spring must have -1 indices and strand.");
    } else {
      RUNTIME_EXP_IFN(a < static_cast<int>(this->nrOfLinks) &&
                        b < static_cast<int>(this->nrOfLinks),
                      "Spring link endpoint index out of range.");
      RUNTIME_EXP_IFN(this->strandIndexOfSpring[sE] >= 0 &&
                        this->strandIndexOfSpring[sE] <
                          static_cast<int>(this->nrOfStrands),
                      "Active spring must map to a valid strand.");
      for (Eigen::Index dir = 0; dir < 3; ++dir) {
        RUNTIME_EXP_IFN(
          this->springCoordinateIndexA[3 * sE + dir] == 3 * a + dir &&
            this->springCoordinateIndexB[3 * sE + dir] == 3 * b + dir,
          "Spring coordinate indices must match 3 * springIndex + dir.");
      }
    }
  }

  for (Eigen::Index idx = 0; idx < this->coordinates.size(); ++idx) {
    RUNTIME_EXP_IFN(std::isfinite(this->coordinates[idx]),
                    "Non-finite link coordinate detected.");
    if (u.size() > 0) {
      RUNTIME_EXP_IFN(std::isfinite(u[idx]),
                      "Non-finite displacement detected.");
    }
  }
  for (Eigen::Index idx = 0; idx < this->springBoxOffset.size(); ++idx) {
    RUNTIME_EXP_IFN(std::isfinite(this->springBoxOffset[idx]),
                    "Non-finite springBoxOffset detected.");
  }
}

} // namespace pylimer_tools::sim::mehp
