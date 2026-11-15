#include "CrosslinkerUniverseStats.h"
#include "../topo/EntanglementDetector.h"
#include "utilityMacros.h"
#include <algorithm>
#include <limits>
#include <random>

namespace pylimer_tools::utils {

namespace {

  std::vector<uint64_t> resolveStrandSeeds(
    const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu)
  {
    const size_t nStrands = cu.strandFrom.size();
    if (cu.strandSeeds.size() == nStrands) {
      return cu.strandSeeds;
    }
    std::mt19937 rng(0);
    std::uniform_int_distribution<uint64_t> seedDist(
      0, std::numeric_limits<uint64_t>::max());
    std::vector<uint64_t> seeds(nStrands);
    for (size_t s = 0; s < nStrands; ++s) {
      seeds[s] = seedDist(rng);
    }
    return seeds;
  }

  double accumulateBond(const pylimer_tools::entities::Box& box,
                        const Eigen::Vector3d& rFrom,
                        const Eigen::Vector3d& rTo,
                        BondLengthStats& stats,
                        double& sumB,
                        double& sumB2)
  {
    Eigen::Vector3d delta = rTo - rFrom;
    box.handlePBC(delta);
    const double len = delta.norm();
    if (stats.nrOfBonds == 0) {
      stats.minBondLength = len;
      stats.maxBondLength = len;
    } else {
      stats.minBondLength = std::min(stats.minBondLength, len);
      stats.maxBondLength = std::max(stats.maxBondLength, len);
    }
    stats.nrOfBonds += 1;
    sumB += len;
    sumB2 += len * len;
    return len;
  }

} // namespace

StrandCategory
classifyStrand(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
               const size_t strandIdx)
{
  INVALIDARG_EXP_IFN(strandIdx < cu.strandFrom.size(),
                     "Strand index out of range.");
  const long int from = cu.strandFrom[strandIdx];
  const long int to = cu.strandTo[strandIdx];
  if (from >= 0 && to >= 0) {
    return (from == to) ? StrandCategory::PRIMARY_LOOP
                        : StrandCategory::NETWORK_STRAND;
  }
  if (from >= 0 || to >= 0) {
    return StrandCategory::DANGLING_CHAIN;
  }
  return StrandCategory::FREE_CHAIN;
}

bool
strandMatchesCategory(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
                      const size_t strandIdx,
                      const StrandCategory category)
{
  if (category == StrandCategory::ALL) {
    return true;
  }
  const StrandCategory actual = classifyStrand(cu, strandIdx);
  if (category == StrandCategory::TWO_CROSSLINKERS_INCLUDING_LOOPS) {
    return actual == StrandCategory::NETWORK_STRAND ||
           actual == StrandCategory::PRIMARY_LOOP;
  }
  return actual == category;
}

BondLengthStats
computeBondLengthStats(const pylimer_tools::entities::Box& box,
                       const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
                       const int nMcSteps)
{
  BondLengthStats stats;
  const size_t nStrands = cu.strandFrom.size();
  if (nStrands == 0) {
    return stats;
  }

  const std::vector<uint64_t> seeds = resolveStrandSeeds(cu);
  double sumB = 0.0;
  double sumB2 = 0.0;

  for (size_t s = 0; s < nStrands; ++s) {
    const long int from = cu.strandFrom[s];
    const long int to = cu.strandTo[s];
    const int nBeads = cu.beadsInStrand[s];

    if (nBeads == 0) {
      if (from >= 0 && to >= 0) {
        const Eigen::Vector3d rFrom(
          cu.xlinkX[from], cu.xlinkY[from], cu.xlinkZ[from]);
        const Eigen::Vector3d rTo(cu.xlinkX[to], cu.xlinkY[to], cu.xlinkZ[to]);
        accumulateBond(box, rFrom, rTo, stats, sumB, sumB2);
      }
      continue;
    }

    const Eigen::VectorXd coords =
      pylimer_tools::topo::entanglement_detection::sampleStrandBeadsOnDemand(
        box, cu, s, seeds[s], nMcSteps);

    if (from >= 0 && to >= 0) {
      const Eigen::Vector3d rFrom(
        cu.xlinkX[from], cu.xlinkY[from], cu.xlinkZ[from]);
      const Eigen::Vector3d rTo(cu.xlinkX[to], cu.xlinkY[to], cu.xlinkZ[to]);
      accumulateBond(box, rFrom, coords.segment<3>(0), stats, sumB, sumB2);
      for (int i = 0; i < nBeads - 1; ++i) {
        accumulateBond(box,
                       coords.segment<3>(3 * i),
                       coords.segment<3>(3 * (i + 1)),
                       stats,
                       sumB,
                       sumB2);
      }
      accumulateBond(
        box, coords.segment<3>(3 * (nBeads - 1)), rTo, stats, sumB, sumB2);
    } else if (from >= 0 || to >= 0) {
      const size_t anchor = static_cast<size_t>(from >= 0 ? from : to);
      const Eigen::Vector3d rAnchor(
        cu.xlinkX[anchor], cu.xlinkY[anchor], cu.xlinkZ[anchor]);
      accumulateBond(box, rAnchor, coords.segment<3>(0), stats, sumB, sumB2);
      for (int i = 0; i < nBeads - 1; ++i) {
        accumulateBond(box,
                       coords.segment<3>(3 * i),
                       coords.segment<3>(3 * (i + 1)),
                       stats,
                       sumB,
                       sumB2);
      }
    } else {
      for (int i = 0; i < nBeads - 1; ++i) {
        accumulateBond(box,
                       coords.segment<3>(3 * i),
                       coords.segment<3>(3 * (i + 1)),
                       stats,
                       sumB,
                       sumB2);
      }
    }
  }

  if (stats.nrOfBonds > 0) {
    const double invN = 1.0 / static_cast<double>(stats.nrOfBonds);
    stats.meanBondLength = sumB * invN;
    stats.meanSquaredBondLength = sumB2 * invN;
  }
  return stats;
}

EndToEndStats
computeEndToEndStats(const pylimer_tools::entities::Box& box,
                     const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
                     const StrandCategory category,
                     const int nMcSteps)
{
  EndToEndStats stats;
  const size_t nStrands = cu.strandFrom.size();
  if (nStrands == 0) {
    return stats;
  }

  const std::vector<uint64_t> seeds = resolveStrandSeeds(cu);
  double sumR = 0.0;
  double sumR2 = 0.0;

  for (size_t s = 0; s < nStrands; ++s) {
    if (!strandMatchesCategory(cu, s, category)) {
      continue;
    }

    const long int from = cu.strandFrom[s];
    const long int to = cu.strandTo[s];
    const int nBeads = cu.beadsInStrand[s];
    double dist = 0.0;

    if (from >= 0 && to >= 0) {
      if (from != to) {
        const Eigen::Vector3d rFrom(
          cu.xlinkX[from], cu.xlinkY[from], cu.xlinkZ[from]);
        const Eigen::Vector3d rTo(cu.xlinkX[to], cu.xlinkY[to], cu.xlinkZ[to]);
        Eigen::Vector3d delta = rTo - rFrom;
        box.handlePBC(delta);
        dist = delta.norm();
      }
    } else if (from >= 0 || to >= 0) {
      if (nBeads >= 1) {
        const size_t anchor = static_cast<size_t>(from >= 0 ? from : to);
        const Eigen::Vector3d rAnchor(
          cu.xlinkX[anchor], cu.xlinkY[anchor], cu.xlinkZ[anchor]);
        const Eigen::VectorXd coords =
          pylimer_tools::topo::entanglement_detection::
            sampleStrandBeadsOnDemand(box, cu, s, seeds[s], nMcSteps);
        Eigen::Vector3d eeVec = coords.segment<3>(0) - rAnchor;
        box.handlePBC(eeVec);
        for (int i = 0; i < nBeads - 1; ++i) {
          Eigen::Vector3d step =
            coords.segment<3>(3 * (i + 1)) - coords.segment<3>(3 * i);
          box.handlePBC(step);
          eeVec += step;
        }
        dist = eeVec.norm();
      }
    } else {
      if (nBeads >= 2) {
        const Eigen::VectorXd coords =
          pylimer_tools::topo::entanglement_detection::
            sampleStrandBeadsOnDemand(box, cu, s, seeds[s], nMcSteps);
        Eigen::Vector3d eeVec = Eigen::Vector3d::Zero();
        for (int i = 0; i < nBeads - 1; ++i) {
          Eigen::Vector3d step =
            coords.segment<3>(3 * (i + 1)) - coords.segment<3>(3 * i);
          box.handlePBC(step);
          eeVec += step;
        }
        dist = eeVec.norm();
      }
    }

    if (dist > 0.0) {
      stats.nrOfStrands += 1;
      sumR += dist;
      sumR2 += dist * dist;
    } else {
      stats.nrOfExcludedZeroLength += 1;
    }
  }

  if (stats.nrOfStrands > 0) {
    const double invN = 1.0 / static_cast<double>(stats.nrOfStrands);
    stats.meanEndToEndDistance = sumR * invN;
    stats.meanSquaredEndToEndDistance = sumR2 * invN;
  }
  return stats;
}

std::map<int, int>
countAtomTypes(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu)
{
  std::map<int, int> counts;
  for (const int xlinkType : cu.xlinkTypes) {
    counts[xlinkType] += 1;
  }
  const size_t nStrands = cu.beadsInStrand.size();
  for (size_t s = 0; s < nStrands; ++s) {
    const int nBeads = cu.beadsInStrand[s];
    if (nBeads > 0) {
      counts[cu.strandBeadType[s]] += nBeads;
    }
  }
  return counts;
}

std::map<int, size_t>
countAtomsByDegree(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu)
{
  std::map<int, size_t> degreeCounts;
  const size_t nCrosslinks = cu.xlinkTypes.size();
  std::vector<int> xlinkDegrees(nCrosslinks, 0);

  const size_t nStrands = cu.strandFrom.size();
  for (size_t s = 0; s < nStrands; ++s) {
    const long int from = cu.strandFrom[s];
    const long int to = cu.strandTo[s];
    const int nBeads = cu.beadsInStrand[s];

    if (from >= 0 && to >= 0) {
      if (!(nBeads == 0 && from == to)) {
        xlinkDegrees[static_cast<size_t>(from)] += 1;
        xlinkDegrees[static_cast<size_t>(to)] += 1;
      }
      if (nBeads > 0) {
        degreeCounts[2] += static_cast<size_t>(nBeads);
      }
    } else if (from >= 0 || to >= 0) {
      const size_t anchor = static_cast<size_t>(from >= 0 ? from : to);
      if (nBeads >= 1) {
        xlinkDegrees[anchor] += 1;
        degreeCounts[1] += 1;
        if (nBeads > 1) {
          degreeCounts[2] += static_cast<size_t>(nBeads - 1);
        }
      }
    } else {
      if (nBeads == 1) {
        degreeCounts[0] += 1;
      } else if (nBeads >= 2) {
        degreeCounts[1] += 2;
        if (nBeads > 2) {
          degreeCounts[2] += static_cast<size_t>(nBeads - 2);
        }
      }
    }
  }

  for (size_t i = 0; i < nCrosslinks; ++i) {
    degreeCounts[xlinkDegrees[i]] += 1;
  }
  return degreeCounts;
}

std::map<int, double>
determineEffectiveFunctionalityPerType(
  const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu)
{
  std::map<int, double> degreeSumByType;
  std::map<int, int> countByType = countAtomTypes(cu);
  for (const auto& [type, count] : countByType) {
    degreeSumByType[type] = 0.0;
  }

  const size_t nCrosslinks = cu.xlinkTypes.size();
  std::vector<int> xlinkDegrees(nCrosslinks, 0);

  const size_t nStrands = cu.strandFrom.size();
  for (size_t s = 0; s < nStrands; ++s) {
    const long int from = cu.strandFrom[s];
    const long int to = cu.strandTo[s];
    const int nBeads = cu.beadsInStrand[s];
    const int beadType = cu.strandBeadType[s];

    if (from >= 0 && to >= 0) {
      if (!(nBeads == 0 && from == to)) {
        xlinkDegrees[static_cast<size_t>(from)] += 1;
        xlinkDegrees[static_cast<size_t>(to)] += 1;
      }
      if (nBeads > 0) {
        degreeSumByType[beadType] += static_cast<double>(2 * nBeads);
      }
    } else if (from >= 0 || to >= 0) {
      const size_t anchor = static_cast<size_t>(from >= 0 ? from : to);
      if (nBeads >= 1) {
        xlinkDegrees[anchor] += 1;
        degreeSumByType[beadType] += static_cast<double>(1 + 2 * (nBeads - 1));
      }
    } else {
      if (nBeads >= 2) {
        degreeSumByType[beadType] += static_cast<double>(2 + 2 * (nBeads - 2));
      }
    }
  }

  for (size_t i = 0; i < nCrosslinks; ++i) {
    degreeSumByType[cu.xlinkTypes[i]] += static_cast<double>(xlinkDegrees[i]);
  }

  for (const auto& [type, count] : countByType) {
    if (count > 0) {
      degreeSumByType[type] /= static_cast<double>(count);
    }
  }
  return degreeSumByType;
}

std::vector<int>
computeStrandBondCounts(const pylimer_tools::sim::mehp::CrosslinkerUniverse& cu,
                        const StrandCategory category)
{
  std::vector<int> bondCounts;
  const size_t nStrands = cu.strandFrom.size();
  bondCounts.reserve(nStrands);

  for (size_t s = 0; s < nStrands; ++s) {
    if (!strandMatchesCategory(cu, s, category)) {
      continue;
    }
    const long int from = cu.strandFrom[s];
    const long int to = cu.strandTo[s];
    const int nBeads = cu.beadsInStrand[s];
    if (from >= 0 && to >= 0) {
      bondCounts.push_back(nBeads + 1);
    } else if (from >= 0 || to >= 0) {
      bondCounts.push_back(nBeads);
    } else {
      bondCounts.push_back(std::max(0, nBeads - 1));
    }
  }
  return bondCounts;
}

} // namespace pylimer_tools::utils
