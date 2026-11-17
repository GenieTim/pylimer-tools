#include "MEHPNetworkBuilders.h"
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace pylimer_tools::sim::mehp {

namespace {

  struct StrandEndpointSummary
  {
    std::vector<size_t> includedStrandIndices;
    std::vector<size_t> danglingEndPartners;
    std::vector<int> danglingEndTypes;
    std::vector<int> freeEndTypes;
    size_t nCrosslinks = 0;
    size_t nDanglingEnds = 0;
    size_t nFreeEnds = 0;
  };

  int getStrandBeadType(const CrosslinkerUniverse& cu, const size_t strandIdx)
  {
    return (strandIdx < cu.strandBeadType.size()) ? cu.strandBeadType[strandIdx]
                                                  : 1;
  }

  StrandEndpointSummary collectIncludedStrands(const CrosslinkerUniverse& cu,
                                               const bool includeDanglingChains,
                                               const bool includeFreeChains)
  {
    StrandEndpointSummary summary;
    summary.nCrosslinks = cu.xlinkTypes.size();
    const size_t nStrands = cu.strandFrom.size();
    summary.includedStrandIndices.reserve(nStrands);

    for (size_t s = 0; s < nStrands; ++s) {
      const long int from = cu.strandFrom[s];
      const long int to = cu.strandTo[s];
      const int beadType = getStrandBeadType(cu, s);

      if (from >= 0 && to >= 0) {
        summary.includedStrandIndices.push_back(s);
      } else if (from >= 0 || to >= 0) {
        if (includeDanglingChains) {
          summary.includedStrandIndices.push_back(s);
          summary.nDanglingEnds += 1;
          summary.danglingEndPartners.push_back(
            static_cast<size_t>(from >= 0 ? from : to));
          summary.danglingEndTypes.push_back(beadType);
        }
      } else {
        if (includeFreeChains && cu.beadsInStrand[s] >= 2) {
          summary.includedStrandIndices.push_back(s);
          summary.nFreeEnds += 2;
          summary.freeEndTypes.push_back(beadType);
          summary.freeEndTypes.push_back(beadType);
        }
      }
    }
    return summary;
  }

  void populateBaseNodeCoordinates(const CrosslinkerUniverse& cu,
                                   const StrandEndpointSummary& summary,
                                   Eigen::VectorXd& coordinates)
  {
    const Eigen::Index nCrosslinks =
      static_cast<Eigen::Index>(summary.nCrosslinks);
    for (Eigen::Index xl = 0; xl < nCrosslinks; ++xl) {
      const size_t xlIdx = static_cast<size_t>(xl);
      const Eigen::Index base = 3 * xl;
      coordinates(base + 0) = cu.xlinkX[xlIdx];
      coordinates(base + 1) = cu.xlinkY[xlIdx];
      coordinates(base + 2) = cu.xlinkZ[xlIdx];
    }
    const Eigen::Index nDanglingEnds =
      static_cast<Eigen::Index>(summary.nDanglingEnds);
    for (Eigen::Index d = 0; d < nDanglingEnds; ++d) {
      const Eigen::Index base = 3 * (nCrosslinks + d);
      const size_t partner =
        summary.danglingEndPartners[static_cast<size_t>(d)];
      coordinates(base + 0) = cu.xlinkX[partner];
      coordinates(base + 1) = cu.xlinkY[partner];
      coordinates(base + 2) = cu.xlinkZ[partner];
    }
    if (summary.nFreeEnds > 0) {
      coordinates
        .segment(3 * (nCrosslinks + nDanglingEnds),
                 static_cast<Eigen::Index>(3 * summary.nFreeEnds))
        .setZero();
    }
  }

  void populateSpringCoordinateIndices(const size_t nSprings,
                                       const Eigen::ArrayXi& springIndexA,
                                       const Eigen::ArrayXi& springIndexB,
                                       Eigen::ArrayXi& springCoordIndexA,
                                       Eigen::ArrayXi& springCoordIndexB)
  {
    springCoordIndexA =
      Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(3 * nSprings));
    springCoordIndexB =
      Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(3 * nSprings));
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(nSprings); ++i) {
      for (Eigen::Index dir = 0; dir < 3; ++dir) {
        springCoordIndexA(3 * i + dir) =
          springIndexA(i) * 3 + static_cast<int>(dir);
        springCoordIndexB(3 * i + dir) =
          springIndexB(i) * 3 + static_cast<int>(dir);
      }
    }
  }

} // namespace

Network
buildForceRelaxationNetwork(const pylimer_tools::entities::Box& box,
                            const CrosslinkerUniverse& cu,
                            const bool includeDanglingChains,
                            const bool includeFreeChains)
{
  Network net;
  net.L[0] = box.getLx();
  net.L[1] = box.getLy();
  net.L[2] = box.getLz();
  net.vol = box.getVolume();

  const StrandEndpointSummary summary =
    collectIncludedStrands(cu, includeDanglingChains, includeFreeChains);
  const size_t nCrosslinks = summary.nCrosslinks;
  const size_t nDanglingEnds = summary.nDanglingEnds;
  const size_t nFreeEnds = summary.nFreeEnds;

  net.nrOfNodes = nCrosslinks + nDanglingEnds + nFreeEnds;
  net.oldAtomIds = net.nrOfNodes > 0
                     ? Eigen::ArrayXi(Eigen::ArrayXi::LinSpaced(
                         static_cast<Eigen::Index>(net.nrOfNodes),
                         0,
                         static_cast<int>(net.nrOfNodes) - 1))
                     : Eigen::ArrayXi(Eigen::ArrayXi::Zero(0));
  net.coordinates = Eigen::VectorXd(net.nrOfNodes * 3);
  populateBaseNodeCoordinates(cu, summary, net.coordinates);

  net.springIndicesOfLinks.assign(net.nrOfNodes, std::vector<size_t>());

  const size_t nSprings = summary.includedStrandIndices.size();
  net.nrOfSprings = nSprings;
  net.springIndexA = Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(nSprings));
  net.springIndexB = Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(nSprings));
  net.springsContourLength =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(nSprings));

  size_t handledDanglingIdx = 0;
  size_t handledFreeIdx = 0;

  for (size_t springIdx = 0; springIdx < nSprings; ++springIdx) {
    const Eigen::Index springEigenIdx = static_cast<Eigen::Index>(springIdx);
    const size_t s = summary.includedStrandIndices[springIdx];
    const long int rawFrom = cu.strandFrom[s];
    const long int rawTo = cu.strandTo[s];
    const int nBeads = cu.beadsInStrand[s];

    size_t nodeFrom = 0;
    size_t nodeTo = 0;

    if (rawFrom >= 0 && rawTo >= 0) {
      nodeFrom = static_cast<size_t>(rawFrom);
      nodeTo = static_cast<size_t>(rawTo);
      net.springsContourLength(springEigenIdx) =
        static_cast<double>(nBeads + 1);
    } else if (rawFrom >= 0 || rawTo >= 0) {
      nodeFrom = static_cast<size_t>(rawFrom >= 0 ? rawFrom : rawTo);
      nodeTo = nCrosslinks + handledDanglingIdx;
      handledDanglingIdx += 1;
      net.springsContourLength(springEigenIdx) = static_cast<double>(nBeads);
    } else {
      nodeFrom = nCrosslinks + nDanglingEnds + handledFreeIdx;
      nodeTo = nodeFrom + 1;
      handledFreeIdx += 2;
      net.springsContourLength(springEigenIdx) =
        static_cast<double>(nBeads - 1);
    }

    net.springIndexA(springEigenIdx) = static_cast<int>(nodeFrom);
    net.springIndexB(springEigenIdx) = static_cast<int>(nodeTo);
    net.springIndicesOfLinks[nodeFrom].push_back(springIdx);
    if (nodeTo != nodeFrom) {
      net.springIndicesOfLinks[nodeTo].push_back(springIdx);
    }
  }

  assert(handledDanglingIdx == nDanglingEnds);
  assert(handledFreeIdx == nFreeEnds);

  populateSpringCoordinateIndices(nSprings,
                                  net.springIndexA,
                                  net.springIndexB,
                                  net.springCoordinateIndexA,
                                  net.springCoordinateIndexB);

  const Eigen::VectorXd actualDistance =
    net.coordinates(net.springCoordinateIndexB) -
    net.coordinates(net.springCoordinateIndexA);
  Eigen::VectorXd expectedDistance = actualDistance;
  box.handlePBC(expectedDistance);
  net.springBoxOffset = expectedDistance - actualDistance;

  net.meanSpringContourLength =
    net.springsContourLength.size() > 0 ? net.springsContourLength.mean() : 1.0;
  net.assumeBoxLargeEnough = true;
  net.assumeComplete = includeDanglingChains && includeFreeChains;

  return net;
}

ForceBalanceNetwork
buildForceBalanceNetwork(const pylimer_tools::entities::Box& box,
                         const CrosslinkerUniverse& cu,
                         const bool includeDanglingChains,
                         const bool includeFreeChains)
{
  ForceBalanceNetwork net;
  net.L[0] = box.getLx();
  net.L[1] = box.getLy();
  net.L[2] = box.getLz();
  for (size_t dir = 0; dir < 3; ++dir) {
    net.boxHalfs[dir] = net.L[dir] * 0.5;
  }
  net.vol = box.getVolume();

  const StrandEndpointSummary summary =
    collectIncludedStrands(cu, includeDanglingChains, includeFreeChains);
  const size_t nCrosslinks = summary.nCrosslinks;
  const size_t nDanglingEnds = summary.nDanglingEnds;
  const size_t nFreeEnds = summary.nFreeEnds;

  net.nrOfNodes = nCrosslinks + nDanglingEnds + nFreeEnds;
  net.nrOfLinks = net.nrOfNodes;
  net.oldAtomIds = net.nrOfNodes > 0
                     ? Eigen::ArrayXi(Eigen::ArrayXi::LinSpaced(
                         static_cast<Eigen::Index>(net.nrOfNodes),
                         0,
                         static_cast<int>(net.nrOfNodes) - 1))
                     : Eigen::ArrayXi(Eigen::ArrayXi::Zero(0));
  net.coordinates = Eigen::VectorXd(net.nrOfNodes * 3);
  populateBaseNodeCoordinates(cu, summary, net.coordinates);

  net.linkIsSliplink =
    Eigen::ArrayXb::Zero(static_cast<Eigen::Index>(net.nrOfLinks));
  net.nrOfCrosslinkSwapsEndured = Eigen::ArrayXi::Zero(0);
  net.oldAtomTypes =
    Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(net.nrOfNodes));
  for (Eigen::Index xl = 0; xl < static_cast<Eigen::Index>(nCrosslinks); ++xl) {
    net.oldAtomTypes(xl) = cu.xlinkTypes[static_cast<size_t>(xl)];
  }
  net.springIndicesOfLinks.assign(net.nrOfNodes, std::vector<size_t>());

  const Eigen::Index nSprings =
    static_cast<Eigen::Index>(summary.includedStrandIndices.size());
  net.springIndexA = Eigen::ArrayXi::Zero(nSprings);
  net.springIndexB = Eigen::ArrayXi::Zero(nSprings);
  net.springsContourLength = Eigen::VectorXd::Zero(nSprings);

  size_t handledDanglingIdx = 0;
  size_t handledFreeIdx = 0;

  for (Eigen::Index springIdx = 0; springIdx < nSprings; ++springIdx) {
    const size_t springSizeIdx = static_cast<size_t>(springIdx);
    const size_t s = summary.includedStrandIndices[springSizeIdx];
    const long int rawFrom = cu.strandFrom[s];
    const long int rawTo = cu.strandTo[s];
    const int nBeads = cu.beadsInStrand[s];

    size_t nodeFrom = 0;
    size_t nodeTo = 0;

    if (rawFrom >= 0 && rawTo >= 0) {
      nodeFrom = static_cast<size_t>(rawFrom);
      nodeTo = static_cast<size_t>(rawTo);
      net.springsContourLength(springIdx) = static_cast<double>(nBeads + 1);
    } else if (rawFrom >= 0 || rawTo >= 0) {
      nodeFrom = static_cast<size_t>(rawFrom >= 0 ? rawFrom : rawTo);
      nodeTo = nCrosslinks + handledDanglingIdx;
      handledDanglingIdx += 1;
      net.springsContourLength(springIdx) = static_cast<double>(nBeads);
    } else {
      nodeFrom = nCrosslinks + nDanglingEnds + handledFreeIdx;
      nodeTo = nodeFrom + 1;
      handledFreeIdx += 2;
      net.springsContourLength(springIdx) = static_cast<double>(nBeads - 1);
    }

    net.springIndexA(springIdx) = static_cast<int>(nodeFrom);
    net.springIndexB(springIdx) = static_cast<int>(nodeTo);
    net.springIndicesOfLinks[nodeFrom].push_back(springSizeIdx);
    if (nodeFrom != nodeTo) {
      net.springIndicesOfLinks[nodeTo].push_back(springSizeIdx);
    }
  }
  assert(handledDanglingIdx == nDanglingEnds);
  assert(handledFreeIdx == nFreeEnds);

  net.nrOfSprings = static_cast<size_t>(nSprings);
  net.springsType = Eigen::ArrayXi::Ones(nSprings);
  net.springIsActive = Eigen::ArrayXb::Ones(nSprings);
  net.partialSpringIsPartial = Eigen::ArrayXb::Zero(nSprings);
  net.nrOfPartialSprings = static_cast<size_t>(nSprings);
  net.nrOfSpringsWithPartition = 0;

  populateSpringCoordinateIndices(static_cast<size_t>(nSprings),
                                  net.springIndexA,
                                  net.springIndexB,
                                  net.springCoordinateIndexA,
                                  net.springCoordinateIndexB);
  net.partialToFullSpringIndex =
    nSprings > 0 ? Eigen::ArrayXi(Eigen::ArrayXi::LinSpaced(
                     nSprings, 0, static_cast<int>(nSprings) - 1))
                 : Eigen::ArrayXi(Eigen::ArrayXi::Zero(0));

  net.linkIndicesOfSprings.reserve(static_cast<size_t>(nSprings));
  net.localToGlobalSpringIndex.reserve(static_cast<size_t>(nSprings));
  net.springToMoleculeIds.reserve(static_cast<size_t>(nSprings));
  for (Eigen::Index i = 0; i < nSprings; ++i) {
    net.linkIndicesOfSprings.push_back(
      { static_cast<size_t>(net.springIndexA(i)),
        static_cast<size_t>(net.springIndexB(i)) });
    net.localToGlobalSpringIndex.push_back({ static_cast<size_t>(i) });
    net.springToMoleculeIds.push_back(static_cast<size_t>(i));
  }

  net.springPartCoordinateIndexA = net.springCoordinateIndexA;
  net.springPartCoordinateIndexB = net.springCoordinateIndexB;
  net.springPartIndexA = net.springIndexA;
  net.springPartIndexB = net.springIndexB;

  const Eigen::VectorXd actualDistance =
    net.coordinates(net.springPartCoordinateIndexB) -
    net.coordinates(net.springPartCoordinateIndexA);
  Eigen::VectorXd expectedDistance = actualDistance;
  box.handlePBC(expectedDistance);
  net.springPartBoxOffset = expectedDistance - actualDistance;

  net.meanSpringContourLength =
    net.springsContourLength.size() > 0 ? net.springsContourLength.mean() : 0.0;

  return net;
}

ForceBalance2Network
buildForceBalance2Network(
  const pylimer_tools::entities::Box& box,
  const CrosslinkerUniverse& cu,
  const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
    entanglements,
  const int /*crossLinkerType*/,
  const bool entanglementsAsSprings,
  const bool includeDanglingChains,
  const bool includeFreeChains)
{
  ForceBalance2Network net;
  net.L[0] = box.getLx();
  net.L[1] = box.getLy();
  net.L[2] = box.getLz();
  for (size_t dir = 0; dir < 3; ++dir) {
    net.boxHalfs[dir] = net.L[dir] * 0.5;
  }

  const size_t nRawStrands = cu.strandFrom.size();
  const size_t nEntanglements = entanglements.size();
  INVALIDARG_EXP_IFN(
    entanglements.strandIdx2.size() == nEntanglements &&
      entanglements.alpha1.size() == nEntanglements &&
      entanglements.alpha2.size() == nEntanglements &&
      static_cast<size_t>(entanglements.x.size()) == nEntanglements &&
      static_cast<size_t>(entanglements.y.size()) == nEntanglements &&
      static_cast<size_t>(entanglements.z.size()) == nEntanglements,
    "Inconsistent StrandPairEntanglements vector sizes.");

  const StrandEndpointSummary summary =
    collectIncludedStrands(cu, includeDanglingChains, includeFreeChains);
  const size_t nCrosslinks = summary.nCrosslinks;
  const size_t nDanglingEnds = summary.nDanglingEnds;
  const size_t nFreeEnds = summary.nFreeEnds;
  const size_t nPolymerStrands = summary.includedStrandIndices.size();

  // Map raw CrosslinkerUniverse strand index -> included strand index (or -1)
  std::vector<long int> rawToIncludedStrand(nRawStrands, -1);
  for (size_t idx = 0; idx < nPolymerStrands; ++idx) {
    rawToIncludedStrand[summary.includedStrandIndices[idx]] =
      static_cast<long int>(idx);
  }

  const size_t nrOfNodes = nCrosslinks + nDanglingEnds + nFreeEnds;
  const size_t nrOfEntanglementLinks =
    entanglementsAsSprings ? (2 * nEntanglements) : nEntanglements;
  const size_t nrOfLinks = nrOfNodes + nrOfEntanglementLinks;

  net.nrOfNodes = nrOfNodes;
  net.nrOfLinks = nrOfLinks;
  net.coordinates =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * nrOfLinks));
  net.linkIsEntanglement =
    Eigen::ArrayXb::Constant(static_cast<Eigen::Index>(nrOfLinks), false);
  net.oldAtomIds = nrOfLinks > 0 ? Eigen::ArrayXi(Eigen::ArrayXi::LinSpaced(
                                     static_cast<Eigen::Index>(nrOfLinks),
                                     0,
                                     static_cast<int>(nrOfLinks) - 1))
                                 : Eigen::ArrayXi(Eigen::ArrayXi::Zero(0));
  net.oldAtomTypes = Eigen::ArrayXi::Ones(static_cast<Eigen::Index>(nrOfLinks));

  for (Eigen::Index xl = 0; xl < static_cast<Eigen::Index>(nCrosslinks); ++xl) {
    const size_t xlIdx = static_cast<size_t>(xl);
    Eigen::Vector3d pos(cu.xlinkX[xlIdx], cu.xlinkY[xlIdx], cu.xlinkZ[xlIdx]);
    box.handlePBC(pos);
    net.coordinates.segment<3>(3 * xl) = pos;
    net.oldAtomTypes(xl) = cu.xlinkTypes[xlIdx];
  }

  for (Eigen::Index dIdx = 0; dIdx < static_cast<Eigen::Index>(nDanglingEnds);
       ++dIdx) {
    const size_t dSizeIdx = static_cast<size_t>(dIdx);
    const Eigen::Index nodeIdx = static_cast<Eigen::Index>(nCrosslinks) + dIdx;
    const Eigen::Index partnerXl =
      static_cast<Eigen::Index>(summary.danglingEndPartners[dSizeIdx]);
    net.coordinates.segment<3>(3 * nodeIdx) =
      net.coordinates.segment<3>(3 * partnerXl);
    net.oldAtomTypes(nodeIdx) = summary.danglingEndTypes[dSizeIdx];
  }

  for (Eigen::Index fIdx = 0; fIdx < static_cast<Eigen::Index>(nFreeEnds);
       ++fIdx) {
    const Eigen::Index nodeIdx =
      static_cast<Eigen::Index>(nCrosslinks + nDanglingEnds) + fIdx;
    net.coordinates.segment<3>(3 * nodeIdx).setZero();
    net.oldAtomTypes(nodeIdx) = summary.freeEndTypes[static_cast<size_t>(fIdx)];
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

    const Eigen::Index entEigenIdx = static_cast<Eigen::Index>(e);
    if (!entanglementsAsSprings) {
      const size_t linkIdx = nrOfNodes + e;
      const Eigen::Index linkEigenIdx = static_cast<Eigen::Index>(linkIdx);
      Eigen::Vector3d mid(entanglements.x[entEigenIdx],
                          entanglements.y[entEigenIdx],
                          entanglements.z[entEigenIdx]);
      box.handlePBC(mid);
      net.coordinates.segment<3>(3 * linkEigenIdx) = mid;
      net.linkIsEntanglement(linkEigenIdx) = true;
      net.oldAtomTypes(linkEigenIdx) = getStrandBeadType(cu, rawS1);

      entanglementsOnStrand[static_cast<size_t>(s1)].push_back(
        StrandEntanglementPlacement{ a1, k1, e, linkIdx });
      entanglementsOnStrand[static_cast<size_t>(s2)].push_back(
        StrandEntanglementPlacement{ a2, k2, e, linkIdx });
    } else {
      const size_t linkIdx1 = nrOfNodes + 2 * e;
      const size_t linkIdx2 = nrOfNodes + 2 * e + 1;
      const Eigen::Index linkEigenIdx1 = static_cast<Eigen::Index>(linkIdx1);
      const Eigen::Index linkEigenIdx2 = static_cast<Eigen::Index>(linkIdx2);
      Eigen::Vector3d pos1 =
        hasSeparatePositions
          ? Eigen::Vector3d(entanglements.position1.segment<3>(3 * entEigenIdx))
          : Eigen::Vector3d(entanglements.x[entEigenIdx],
                            entanglements.y[entEigenIdx],
                            entanglements.z[entEigenIdx]);
      Eigen::Vector3d pos2 =
        hasSeparatePositions
          ? Eigen::Vector3d(entanglements.position2.segment<3>(3 * entEigenIdx))
          : Eigen::Vector3d(entanglements.x[entEigenIdx],
                            entanglements.y[entEigenIdx],
                            entanglements.z[entEigenIdx]);
      box.handlePBC(pos1);
      box.handlePBC(pos2);
      net.coordinates.segment<3>(3 * linkEigenIdx1) = pos1;
      net.coordinates.segment<3>(3 * linkEigenIdx2) = pos2;
      net.linkIsEntanglement(linkEigenIdx1) = true;
      net.linkIsEntanglement(linkEigenIdx2) = true;
      net.oldAtomTypes(linkEigenIdx1) = getStrandBeadType(cu, rawS1);
      net.oldAtomTypes(linkEigenIdx2) = getStrandBeadType(cu, rawS2);

      entanglementsOnStrand[static_cast<size_t>(s1)].push_back(
        StrandEntanglementPlacement{ a1, k1, e, linkIdx1 });
      entanglementsOnStrand[static_cast<size_t>(s2)].push_back(
        StrandEntanglementPlacement{ a2, k2, e, linkIdx2 });
    }
  }

  const size_t nrOfTotalStrands =
    nPolymerStrands + (entanglementsAsSprings ? nEntanglements : 0);
  const size_t nrOfTotalSprings =
    nPolymerStrands +
    (entanglementsAsSprings ? (3 * nEntanglements) : (2 * nEntanglements));
  const Eigen::Index nrOfTotalSpringsEigen =
    static_cast<Eigen::Index>(nrOfTotalSprings);

  net.nrOfStrands = nrOfTotalStrands;
  net.nrOfSprings = nrOfTotalSprings;
  net.springIndicesOfStrand.resize(nrOfTotalStrands);
  net.linkIndicesOfStrand.resize(nrOfTotalStrands);
  net.strandIndicesOfLink.assign(nrOfLinks, std::vector<size_t>());

  net.springIndexA = Eigen::ArrayXi::Zero(nrOfTotalSpringsEigen);
  net.springIndexB = Eigen::ArrayXi::Zero(nrOfTotalSpringsEigen);
  net.springCoordinateIndexA = Eigen::ArrayXi::Zero(3 * nrOfTotalSpringsEigen);
  net.springCoordinateIndexB = Eigen::ArrayXi::Zero(3 * nrOfTotalSpringsEigen);
  net.springContourLength = Eigen::VectorXd::Zero(nrOfTotalSpringsEigen);
  net.springBoxOffset = Eigen::VectorXd::Zero(3 * nrOfTotalSpringsEigen);
  net.springIsEntanglement =
    Eigen::ArrayXb::Constant(nrOfTotalSpringsEigen, false);
  net.strandIndexOfSpring = Eigen::ArrayXi::Constant(nrOfTotalSpringsEigen, -1);

  size_t handledDanglingIdx = 0;
  size_t handledFreeIdx = 0;
  size_t nextSpringIdx = 0;

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

    for (const size_t linkIdx : net.linkIndicesOfStrand[s]) {
      pylimer_tools::utils::addIfNotContained(net.strandIndicesOfLink[linkIdx],
                                              s);
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
      const Eigen::Index spEigenIdx = static_cast<Eigen::Index>(spIdx);
      const size_t linkA = net.linkIndicesOfStrand[s][m];
      const size_t linkB = net.linkIndicesOfStrand[s][m + 1];

      net.springIndicesOfStrand[s].push_back(spIdx);
      net.strandIndexOfSpring(spEigenIdx) = static_cast<int>(s);
      net.springIndexA(spEigenIdx) = static_cast<int>(linkA);
      net.springIndexB(spEigenIdx) = static_cast<int>(linkB);
      for (Eigen::Index dir = 0; dir < 3; ++dir) {
        net.springCoordinateIndexA(3 * spEigenIdx + dir) =
          static_cast<int>(3 * linkA) + static_cast<int>(dir);
        net.springCoordinateIndexB(3 * spEigenIdx + dir) =
          static_cast<int>(3 * linkB) + static_cast<int>(dir);
      }
      net.springContourLength(spEigenIdx) = subContourLengths[m];

      const Eigen::Vector3d shouldBeDist =
        unwrappedPositions[m + 1] - unwrappedPositions[m];
      const Eigen::Vector3d isDist =
        net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * linkB)) -
        net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * linkA));
      net.springBoxOffset.segment<3>(3 * spEigenIdx) = shouldBeDist - isDist;
      assert(box.isValidOffset(net.springBoxOffset.segment<3>(3 * spEigenIdx)));
    }
  }

  if (entanglementsAsSprings) {
    for (size_t e = 0; e < nEntanglements; ++e) {
      const size_t entStrandIdx = nPolymerStrands + e;
      const size_t spIdx = nextSpringIdx++;
      const Eigen::Index spEigenIdx = static_cast<Eigen::Index>(spIdx);
      const size_t link1 = nrOfNodes + 2 * e;
      const size_t link2 = nrOfNodes + 2 * e + 1;

      net.linkIndicesOfStrand[entStrandIdx] = { link1, link2 };
      net.springIndicesOfStrand[entStrandIdx] = { spIdx };
      pylimer_tools::utils::addIfNotContained(net.strandIndicesOfLink[link1],
                                              entStrandIdx);
      pylimer_tools::utils::addIfNotContained(net.strandIndicesOfLink[link2],
                                              entStrandIdx);

      net.strandIndexOfSpring(spEigenIdx) = static_cast<int>(entStrandIdx);
      net.springIsEntanglement(spEigenIdx) = true;
      net.springContourLength(spEigenIdx) = 1.0;
      net.springIndexA(spEigenIdx) = static_cast<int>(link1);
      net.springIndexB(spEigenIdx) = static_cast<int>(link2);
      for (Eigen::Index dir = 0; dir < 3; ++dir) {
        net.springCoordinateIndexA(3 * spEigenIdx + dir) =
          static_cast<int>(3 * link1) + static_cast<int>(dir);
        net.springCoordinateIndexB(3 * spEigenIdx + dir) =
          static_cast<int>(3 * link2) + static_cast<int>(dir);
      }
      net.springBoxOffset.segment<3>(3 * spEigenIdx) = box.getOffset(
        net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * link2)) -
        net.coordinates.segment<3>(static_cast<Eigen::Index>(3 * link1)));
    }
  }

  assert(nextSpringIdx == nrOfTotalSprings);
  return net;
}

} // namespace pylimer_tools::sim::mehp
