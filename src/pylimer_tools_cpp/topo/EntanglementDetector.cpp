#include "EntanglementDetector.h"

#include "../entities/EigenNeighbourList.h"
#include "../entities/Molecule.h"
#include "../entities/Universe.h"
#include "../utils/LazyDistanceMatrix.h"
#include "../utils/RandomWalker.h"
#include "../utils/StringUtils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
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

  igraph_lazy_distance_matrix_state_t distance_matrix_computer;
  igraph_t graph = universe.getCopyOfGraph();
  igraph_lazy_distance_matrix_state_init(
    &distance_matrix_computer,
    &graph,
    IGRAPH_ALL,
    sameStrandCutoff > 0 ? static_cast<igraph_integer_t>(sameStrandCutoff) : 0);

  std::mt19937 rng = initializeRng(seed);

  const pylimer_tools::entities::Box box = universe.getBox();
  const Eigen::VectorXd coordinates =
    universe.getUnwrappedVertexCoordinates(box);
  const pylimer_tools::entities::EigenNeighbourList neighbourList(
    coordinates, box, upperCutoff);
  std::vector<bool> vertexIsEligible(universe.getNrOfAtoms(), true);
  const std::vector<int> vertexTypes = universe.getPropertyValues<int>("type");
  const std::vector<int> vertexDegrees = universe.getVertexDegrees();

  const std::vector<bool> activeAtoms =
    (filteringMode != EntanglementFilteringMode::NONE)
      ? universe.identifyActiveAtoms(crossLinkerType)
      : std::vector<bool>();

  std::vector<long int> vertexIdsToSample;
  vertexIdsToSample.reserve(universe.getNrOfAtoms());

  for (size_t i = 0; i < universe.getNrOfAtoms(); ++i) {
    if (vertexDegrees[i] != 2) {
      vertexIsEligible[i] = false;
    }
    if (ignoreCrosslinks && vertexTypes[i] == crossLinkerType) {
      vertexIsEligible[i] = false;
    }
    if (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
        !activeAtoms.empty() && !activeAtoms[i]) {
      vertexIsEligible[i] = false;
    }
    if (vertexIsEligible[i]) {
      vertexIdsToSample.push_back(static_cast<long int>(i));
    }
  }

  size_t numLinksFoundInIteration = 1;
  Eigen::ArrayXi neighbourIndices = Eigen::ArrayXi::Zero(12);
  size_t numEntanglementsSampled = 0;

  do {
    numLinksFoundInIteration = 0;
    std::ranges::shuffle(vertexIdsToSample, rng);
    for (const long int atomVertexIdx1 : vertexIdsToSample) {
      if (pairOfAtom[atomVertexIdx1] != -1) {
        continue;
      }
      vertexIsEligible[atomVertexIdx1] = false;

      const size_t nNeighbours = neighbourList.getIndicesCloseToCoordinates(
        neighbourIndices,
        coordinates.segment<3>(3 * atomVertexIdx1),
        upperCutoff);
      long int atomVertexIdx2 = -1;

      std::shuffle(
        neighbourIndices.data(), neighbourIndices.data() + nNeighbours, rng);

      for (Eigen::Index idxInNeighbours = 0;
           idxInNeighbours < static_cast<Eigen::Index>(nNeighbours);
           ++idxInNeighbours) {
        const long int candidateVertex = neighbourIndices[idxInNeighbours];
        if (!vertexIsEligible[candidateVertex] ||
            pairOfAtom[candidateVertex] != -1) {
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

        igraph_integer_t pathLength =
          sameStrandCutoff <= 0
            ? IGRAPH_INFINITY
            : igraph_lazy_distance_matrix_path_length(&distance_matrix_computer,
                                                      atomVertexIdx1,
                                                      candidateVertex,
                                                      false);
        if (sameStrandCutoff <= 0 || pathLength >= sameStrandCutoff ||
            pathLength == IGRAPH_INFINITY) {
          atomVertexIdx2 = candidateVertex;
          break;
        }
      }
      igraph_lazy_distance_matrix_forget_source(&distance_matrix_computer,
                                                atomVertexIdx1);

      if (atomVertexIdx2 < 0) {
        continue;
      }

      RUNTIME_EXP_IFN(pairOfAtom[atomVertexIdx2] == -1,
                      "Expected not to be able to sample the same atom twice.");

      if (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT) {
        pairOfAtom[atomVertexIdx2] = static_cast<long int>(pairsOfAtoms.size());
        pairOfAtom[atomVertexIdx1] = static_cast<long int>(pairsOfAtoms.size());
        pairsOfAtoms.emplace_back(universe.getAtomIdByIdx(atomVertexIdx1),
                                  universe.getAtomIdByIdx(atomVertexIdx2));
      } else if (filteringMode == EntanglementFilteringMode::OMIT_INACTIVE) {
        if (activeAtoms.empty() ||
            (activeAtoms[atomVertexIdx1] && activeAtoms[atomVertexIdx2])) {
          pairOfAtom[atomVertexIdx2] =
            static_cast<long int>(pairsOfAtoms.size());
          pairOfAtom[atomVertexIdx1] =
            static_cast<long int>(pairsOfAtoms.size());
          pairsOfAtoms.emplace_back(universe.getAtomIdByIdx(atomVertexIdx1),
                                    universe.getAtomIdByIdx(atomVertexIdx2));
        } else {
          pairOfAtom[atomVertexIdx2] = -2;
          pairOfAtom[atomVertexIdx1] = -2;
        }
      } else { // NONE
        pairOfAtom[atomVertexIdx2] = static_cast<long int>(pairsOfAtoms.size());
        pairOfAtom[atomVertexIdx1] = static_cast<long int>(pairsOfAtoms.size());
        pairsOfAtoms.emplace_back(universe.getAtomIdByIdx(atomVertexIdx1),
                                  universe.getAtomIdByIdx(atomVertexIdx2));
      }

      vertexIsEligible[atomVertexIdx2] = false;
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
    uint32_t nBeads = 0;
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

  Eigen::VectorXd sampleStrandBeadsOnDemand(
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
      Eigen::Vector3d rFrom(cu.xlinkX[from], cu.xlinkY[from], cu.xlinkZ[from]);
      Eigen::Vector3d rTo(cu.xlinkX[to], cu.xlinkY[to], cu.xlinkZ[to]);
      Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChainFromTo(
        box, rFrom, rTo, nBeads, beadDist, meanSqDist, strandRng, true);
      if (nMcSteps > 0) {
        pylimer_tools::sim::equilibrateChainWithMC(
          positions, meanSqDist, strandRng, true, true, nMcSteps);
      }
      return positions.segment(3, 3 * nBeads);
    } else if (from >= 0 || to >= 0) {
      const size_t anchor = static_cast<size_t>(from >= 0 ? from : to);
      Eigen::Vector3d rAnchor(
        cu.xlinkX[anchor], cu.xlinkY[anchor], cu.xlinkZ[anchor]);
      Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
        nBeads, beadDist, meanSqDist, strandRng);
      if (nMcSteps > 0) {
        pylimer_tools::sim::equilibrateChainWithMC(
          positions, meanSqDist, strandRng, true, false, nMcSteps);
      }
      positions += rAnchor.replicate(nBeads, 1);
      return positions;
    } else {
      std::uniform_real_distribution<double> distX(box.getLowX(),
                                                   box.getHighX());
      std::uniform_real_distribution<double> distY(box.getLowY(),
                                                   box.getHighY());
      std::uniform_real_distribution<double> distZ(box.getLowZ(),
                                                   box.getHighZ());
      Eigen::Vector3d rFrom(
        distX(strandRng), distY(strandRng), distZ(strandRng));
      Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
        nBeads, beadDist, meanSqDist, strandRng);
      if (nMcSteps > 0) {
        pylimer_tools::sim::equilibrateChainWithMC(
          positions, meanSqDist, strandRng, true, false, nMcSteps);
      }
      positions += rFrom.replicate(nBeads, 1);
      return positions;
    }
  }

  class StrandCache
  {
  public:
    explicit StrandCache(size_t capacity = 4096)
      : capacity_(std::max<size_t>(1, capacity))
    {
      indexMap_.reserve(capacity_);
    }

    template<typename Generator>
    const Eigen::VectorXd& get(size_t strandIdx, Generator&& gen)
    {
      const auto it = indexMap_.find(strandIdx);
      if (it != indexMap_.end()) {
        if (it->second != entries_.begin()) {
          entries_.splice(entries_.begin(), entries_, it->second);
        }
        return entries_.front().coords;
      }
      Eigen::VectorXd coords = gen();
      if (entries_.size() >= capacity_) {
        indexMap_.erase(entries_.back().strandIdx);
        entries_.pop_back();
      }
      entries_.push_front(Entry{ strandIdx, std::move(coords) });
      indexMap_[strandIdx] = entries_.begin();
      return entries_.front().coords;
    }

  private:
    struct Entry
    {
      size_t strandIdx;
      Eigen::VectorXd coords;
    };
    size_t capacity_;
    std::list<Entry> entries_;
    std::unordered_map<size_t, std::list<Entry>::iterator> indexMap_;
  };

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
  const EntanglementFilteringMode filteringMode)
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

  std::mt19937 rng = initializeRng(seed);

  // Topological leaf peeling on CrosslinkerUniverse if filtering is enabled
  std::vector<bool> strandIsIgnored(nStrands, false);
  if (filteringMode != EntanglementFilteringMode::NONE) {
    strandIsIgnored = peelTopologicalLeaves(crosslinkerUniverse);
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
    if (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
        strandIsIgnored[s]) {
      continue;
    }
    const size_t nInt = static_cast<size_t>(strands[s].nInternalBeads);
    for (size_t m = 0; m < nInt; ++m) {
      candidateBeads.push_back(
        BeadRef{ static_cast<uint32_t>(s), static_cast<uint32_t>(m) });
    }
  }

  const size_t nActiveCandidates = candidateBeads.size();
  if (nActiveCandidates < 2) {
    return result;
  }

  std::vector<uint8_t> beadIsPaired(totalInternalBeads, 0);
  std::vector<size_t> candidateOrder(nActiveCandidates);
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

    for (const size_t candIdx1 : candidateOrder) {
      const size_t s1 = candidateBeads[candIdx1].strandIdx;
      const size_t m1 = candidateBeads[candIdx1].internalIdx;
      const size_t globalIdx1 = strandInternalOffset[s1] + m1;
      if (beadIsPaired[globalIdx1]) {
        continue;
      }

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
        if (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
            sm2.isIgnored) {
          continue;
        }
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

      if (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT ||
          filteringMode == EntanglementFilteringMode::NONE ||
          (!sm1.isIgnored && !sm2.isIgnored)) {
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

      const size_t targetCount =
        (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT)
          ? result.size()
          : numEntanglementsSampled;
      if (targetCount >= nrOfEntanglementsToSample) {
        break;
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
  std::uniform_int_distribution<uint64_t> seedDist(
    0, std::numeric_limits<uint64_t>::max());
  for (size_t s = 0; s < nStrands; ++s) {
    strandSeeds[s] = seedDist(masterRng);
  }

  std::vector<StrandBoundingBox> bboxes(nStrands);
  double maxRadius = 0.0;
  std::vector<size_t> strandInternalOffset(nStrands + 1, 0);

  for (size_t s = 0; s < nStrands; ++s) {
    const int nBeads = crosslinkerUniverse.beadsInStrand[s];
    bboxes[s].seed = strandSeeds[s];
    bboxes[s].nBeads = std::max(0, nBeads);
    bboxes[s].isIgnored = strandIsIgnored[s];

    const bool ignoreInOffset =
      (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
       strandIsIgnored[s]) ||
      nBeads <= 0;
    strandInternalOffset[s + 1] =
      strandInternalOffset[s] +
      (ignoreInOffset ? 0 : static_cast<size_t>(nBeads));

    if (nBeads <= 0 ||
        (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
         strandIsIgnored[s])) {
      continue;
    }

    Eigen::VectorXd coords = sampleStrandBeadsOnDemand(
      box, crosslinkerUniverse, s, strandSeeds[s], nMcSteps);
    if (coords.size() < 3 * nBeads) {
      continue;
    }

    Eigen::Vector3d minCoord = coords.segment<3>(0);
    Eigen::Vector3d maxCoord = coords.segment<3>(0);
    for (int i = 1; i < nBeads; ++i) {
      Eigen::Vector3d pt = coords.segment<3>(3 * i);
      minCoord = minCoord.cwiseMin(pt);
      maxCoord = maxCoord.cwiseMax(pt);
    }

    Eigen::Vector3d unwrappedCenter = 0.5 * (minCoord + maxCoord);
    bboxes[s].halfExtents = 0.5 * (maxCoord - minCoord);
    bboxes[s].center = unwrappedCenter;
    box.handlePBC(bboxes[s].center);
    bboxes[s].radius = bboxes[s].halfExtents.norm();
    if (bboxes[s].radius > maxRadius) {
      maxRadius = bboxes[s].radius;
    }
  }

  const size_t totalInternalBeads = strandInternalOffset[nStrands];
  if (totalInternalBeads < 2) {
    return result;
  }

  Eigen::VectorXd strandCenters(3 * nStrands);
  for (size_t s = 0; s < nStrands; ++s) {
    strandCenters.segment<3>(3 * s) = bboxes[s].center;
  }
  const double searchCutoff = maxRadius + upperCutoff;
  const pylimer_tools::entities::EigenNeighbourList strandNeighbourList(
    strandCenters, box, searchCutoff);

  std::vector<uint8_t> beadIsPaired(totalInternalBeads, 0);
  std::vector<size_t> candidateOrder(totalInternalBeads);
  std::iota(candidateOrder.begin(), candidateOrder.end(), 0);

  StrandCache cache(4096);

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

  Eigen::ArrayXi nearbyStrands = Eigen::ArrayXi::Zero(32);

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
    std::ranges::shuffle(candidateOrder, masterRng);

    for (const size_t g1 : candidateOrder) {
      if (beadIsPaired[g1]) {
        continue;
      }

      auto it1 = std::upper_bound(strandInternalOffset.begin(),
                                  strandInternalOffset.end(),
                                  g1) -
                 1;
      const size_t s1 = std::distance(strandInternalOffset.begin(), it1);
      const size_t m1 = g1 - *it1;

      const Eigen::VectorXd& coords1 = cache.get(s1, [&]() {
        return sampleStrandBeadsOnDemand(
          box, crosslinkerUniverse, s1, strandSeeds[s1], nMcSteps);
      });
      const Eigen::Vector3d r1 = coords1.segment<3>(3 * m1);
      Eigen::Vector3d r1Wrapped = r1;
      box.handlePBC(r1Wrapped);

      candidatePartners.clear();

      // 1. Same-strand candidate beads
      const size_t nBeads1 =
        static_cast<size_t>(crosslinkerUniverse.beadsInStrand[s1]);
      const size_t baseOffset1 = *it1;
      for (size_t m2 = 0; m2 < nBeads1; ++m2) {
        if (m1 == m2) {
          continue;
        }
        if (sameStrandCutoff > 0.0) {
          const double bondDist =
            std::abs(static_cast<double>(m1) - static_cast<double>(m2));
          if (bondDist < sameStrandCutoff) {
            continue;
          }
        }
        const size_t g2 = baseOffset1 + m2;
        if (beadIsPaired[g2]) {
          continue;
        }
        Eigen::Vector3d r2 = coords1.segment<3>(3 * m2);
        Eigen::Vector3d delta = r2 - r1;
        box.handlePBC(delta);
        const double d2 = delta.squaredNorm();
        if (d2 <= upperCutoff * upperCutoff &&
            d2 >= lowerCutoff * lowerCutoff) {
          candidatePartners.push_back({ s1, m2, g2, delta });
        }
      }

      // 2. Different-strand candidate beads
      const size_t nNearby = strandNeighbourList.getIndicesCloseToCoordinates(
        nearbyStrands, r1Wrapped, searchCutoff, true);

      for (size_t nsIdx = 0; nsIdx < nNearby; ++nsIdx) {
        const size_t s2 = static_cast<size_t>(nearbyStrands[nsIdx]);
        if (s2 == s1 ||
            (filteringMode == EntanglementFilteringMode::PRUNE_UPFRONT &&
             bboxes[s2].isIgnored) ||
            bboxes[s2].nBeads == 0) {
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
        const size_t nBeads2 = bboxes[s2].nBeads;
        const size_t baseOffset2 = strandInternalOffset[s2];

        for (size_t m2 = 0; m2 < nBeads2; ++m2) {
          const size_t g2 = baseOffset2 + m2;
          if (beadIsPaired[g2]) {
            continue;
          }
          Eigen::Vector3d r2 = coords2.segment<3>(3 * m2);
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
        result.strandIdx1.push_back(s1);
        result.strandIdx2.push_back(chosen.strandIdx);
        result.alpha1.push_back(
          static_cast<double>(m1 + 1) /
          static_cast<double>(crosslinkerUniverse.beadsInStrand[s1] + 1));
        result.alpha2.push_back(
          static_cast<double>(chosen.beadIdx + 1) /
          static_cast<double>(
            crosslinkerUniverse.beadsInStrand[chosen.strandIdx] + 1));
        result.beadIdx1.push_back(m1);
        result.beadIdx2.push_back(chosen.beadIdx);

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
