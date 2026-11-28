#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPLinearSystemSolver.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPNetwork.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPNonlinearOptimizationSolver.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPSimulator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPSlipLinkRelaxationSolver.h"
#include "../../src/pylimer_tools_cpp/topo/EntanglementDetector.h"
#include "../../src/pylimer_tools_cpp/utils/ExtraEigenTypes.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <random>
#include <vector>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::sim::mehp;

double
roundForOutput(const double val, const double precision = 5)
{
  return std::round(val * std::pow(10, precision)) / std::pow(10, precision);
}

void
outputNetwork(const pcm::MEHPNetwork& net,
              const Eigen::VectorXd& displacements,
              const Eigen::VectorXd& springPartitions)
{
  for (size_t i = 0; i < net.nrOfStrands; ++i) {
    std::cout << "Strand " << i << ", N: "
              << net.strandContourLength[static_cast<Eigen::Index>(i)]
              << std::endl;
    for (size_t j = 0; j < net.linkIndicesOfStrand[i].size(); ++j) {
      const size_t linkIdx = net.linkIndicesOfStrand[i][j];
      std::cout << linkIdx;
      const Eigen::Index linkEigenIdx = static_cast<Eigen::Index>(linkIdx);
      if (!net.linkIsEntanglement[linkEigenIdx]) {
        std::cout << " (" << net.oldAtomIds[linkEigenIdx] << ")";
      }
      std::cout << ": ";
      for (Eigen::Index dir = 0; dir < 3; ++dir) {
        std::cout
          << (net.coordinates[3 * static_cast<Eigen::Index>(linkIdx) + dir] +
              displacements[3 * static_cast<Eigen::Index>(linkIdx) + dir])
          << ", ";
      }
      std::cout << std::endl;
      if (j < net.linkIndicesOfStrand[i].size() - 1) {
        const size_t partialSpringIdx = net.springIndicesOfStrand[i][j];
        std::cout
          << partialSpringIdx << ": "
          << springPartitions[static_cast<Eigen::Index>(partialSpringIdx)]
          << " ("
          << roundForOutput(net.springBoxOffset[3 * static_cast<Eigen::Index>(
                                                      partialSpringIdx) +
                                                0])
          << ", "
          << roundForOutput(net.springBoxOffset[3 * static_cast<Eigen::Index>(
                                                      partialSpringIdx) +
                                                1])
          << ", "
          << roundForOutput(net.springBoxOffset[3 * static_cast<Eigen::Index>(
                                                      partialSpringIdx) +
                                                2])
          << ")" << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

TEST_CASE("Particular slip-link examples", "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"Particular slip-link examples\"" << std::endl;
  const double L = 42.819955007276754;
  const double lmda = 1.2;
  pe::Universe universe = pe::Universe(L, L, L);
  /**
   * Connectivity:
   *
   * 35-(11)-90
   *
   * 10-(12)-1654
   */
  // slip-link 16 in the test-system
  universe.addAtoms({ 35, 90, 10, 1654, 11, 12, 13, 14 },
                    { 2, 2, 2, 2, 1, 1, 1, 1 },
                    { 19.706880857235795,
                      19.288603889563976,
                      10,
                      22.156152142687819,
                      10,
                      10,
                      10,
                      10 },
                    { 1.47224612942217,
                      4.4207926048800461,
                      10,
                      2.821003235624608,
                      10,
                      10,
                      10,
                      10 },
                    { 22.98584649043724,
                      24.328555987207494,
                      10,
                      24.562956857368366,
                      10,
                      10,
                      10,
                      10 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 });
  universe.addBonds({ 35, 11, 10, 12 }, { 11, 90, 12, 1654 });

  universe.setBox(
    pe::Box(lmda * L, L * (1. / sqrt(lmda)), L * (1. / sqrt(lmda))), false);

  pcm::MEHPNetwork net = pcm::MEHPNetwork::fromUniverse(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, true, false, true);
  Eigen::VectorXd displacements = Eigen::VectorXd::Zero(4 * 3);
  net.addSlipLinks(
    { 1, 0 },
    { 1, 1 },
    (Eigen::VectorXd(2) << 22.650980696700437, 15.67228278755659).finished(),
    (Eigen::VectorXd(2) << 38.666887107085628, 3.172504627794056).finished(),
    (Eigen::VectorXd(2) << 20.999752090751294, 27.126973321608833).finished(),
    { 1. - 0.448276, 0.31034482758620685 },
    { 1. - 0.448276, 1. - (0.6896551724137931 * 0.448276) },
    &displacements);

  CHECK(net.springPartitions[0] == Catch::Approx(0.31034482758620685));
  CHECK(net.springPartitions[1] == Catch::Approx(1. - 0.448276));
  CHECK(net.springPartitions[3] == Catch::Approx(1e-7).margin(1e-12));
  CHECK(net.springPartitions[2] == Catch::Approx(0.13912));
  CHECK(net.springPartitions[4] == Catch::Approx(1. - 0.31034482758620685));
  for (int iter = 0; iter < 250; ++iter) {
    const double partRes =
      net.updateSpringPartitionAtLink(5, displacements, 1e10, false);
    const double stepNorm =
      net.displaceToMeanPosition(displacements, 5, false, 1e10);
    if (partRes < 1e-10 && stepNorm < 1e-10) {
      break;
    }
  }
  CHECK(displacements[5 * 3] == Catch::Approx(5.89283));
  CHECK(displacements[5 * 3 + 1] == Catch::Approx(-1.0033));
  CHECK(displacements[5 * 3 + 2] == Catch::Approx(-3.49778));
  CHECK(net.springPartitions[0] == Catch::Approx(0.388968));
  CHECK(net.springPartitions[5] == Catch::Approx(0.112087));
  CHECK(net.springPartitions[3] == Catch::Approx(1e-7).margin(1e-12));
  CHECK(net.springPartitions[1] == Catch::Approx(1. - 0.448276));
  CHECK(net.springPartitions[4] == Catch::Approx(1. - 0.388968));
  CHECK(net.springPartitions[2] == Catch::Approx((1. - 0.25004) * 0.448276));
}

TEST_CASE("MC swap accept and reject work", "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MC swap accept and reject work\"" << std::endl;
  const double L = 42.819955007276754;
  pe::Universe universe = pe::Universe(L, L, L);
  /**
   * The universe looks something like this:
   * Connectivity:
   *
   * 1-(7)-2
   *
   * 3-(8)-4
   *
   * 5-(9)-6
   */
  universe.addAtoms({ 1, 2, 3, 4, 5, 6, 7, 8, 9 },
                    { 2, 2, 2, 2, 2, 2, 1, 1, 1 },
                    { 20., 24., 20., 24., 20., 24., 22., 22., 22. },
                    { 20., 20., 16., 16., 12., 12., 20., 16., 12. },
                    {
                      0.,
                      0.,
                      0.,
                      0.,
                      0.,
                      0.,
                      0.,
                      0.,
                      0.,
                    },
                    { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                    { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                    { 0, 0, 0, 0, 0, 0, 0, 0, 0 });
  universe.addBonds({ 1, 7, 3, 8, 5, 9 }, { 7, 2, 8, 4, 9, 6 });
  pcm::MEHPNetwork baseNet = pcm::MEHPNetwork::fromUniverse(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, true, false, true);
  REQUIRE(baseNet.nrOfStrands == 3);
  baseNet.strandContourLength = Eigen::VectorXd::Constant(3, 30.0);
  baseNet.springContourLength = baseNet.strandContourLength;
  baseNet.recomputeMeanContourLengths();

  SECTION("MC condition accepts as requested")
  {
    pcm::MEHPNetwork net = baseNet;
    Eigen::VectorXd u =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    net.addSlipLinks({ 0, 0 },
                     { 1, 1 },
                     (Eigen::VectorXd(2) << 21.03, 23.01).finished(),
                     (Eigen::VectorXd(2) << 18.02, 18.03).finished(),
                     Eigen::VectorXd::Zero(2),
                     { 0.505, 0.48 },
                     { 0.49, 0.51 },
                     &u);
    const pcm::MEHPNetwork netBefore = net;
    const Eigen::VectorXd uBefore = u;
    CHECK(net.moveSlipLinkAlongStrand(
      5, u, pcm::LinkSwappingMode::ALL_MC_TRY, 1., -1, false));
    const bool netIsIdentical =
      net.springIndexA.isApprox(netBefore.springIndexA) &&
      net.springIndexB.isApprox(netBefore.springIndexB);
    CHECK_FALSE(netIsIdentical);
    CHECK_FALSE(u.isApprox(uBefore));
  }

  SECTION("MC condition rejects as requested")
  {
    pcm::MEHPNetwork net = baseNet;
    Eigen::VectorXd u =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    net.addSlipLinks({ 0, 0 },
                     { 1, 1 },
                     (Eigen::VectorXd(2) << 21., 23.).finished(),
                     (Eigen::VectorXd(2) << 18., 18.).finished(),
                     Eigen::VectorXd::Zero(2),
                     { 0.4, 0.6 },
                     { 0.4, 0.6 },
                     &u);
    const Eigen::VectorXd partitionsBefore = net.springPartitions;
    const Eigen::VectorXd uBefore = u;
    CHECK_FALSE(net.moveSlipLinkAlongStrand(
      5, u, pcm::LinkSwappingMode::ALL_MC_TRY, 1., -1, false));
    CHECK(net.springPartitions.isApprox(partitionsBefore));
    CHECK(u.isApprox(uBefore));
  }
}

TEST_CASE("MC swap accept and reject work with crosslinkers",
          "[analysis][MEHPSimulator]")
{
  std::cout
    << "Running test \"MC swap accept and reject work with crosslinkers\""
    << std::endl;
  const double L = 42.819955007276754;
  pe::Universe universe = pe::Universe(L, L, L);
  /**
   * Connectivity:
   *
   * 1-(7)-2
   *
   * 3-(8)-4-(9)-5
   *      (10)
   *       6
   */
  universe.addAtoms({ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 },
                    { 2, 2, 2, 2, 2, 2, 1, 1, 1, 1 },
                    { 10., 20., 5., 10., 15., 20., 0., 0., 0., 0. },
                    { 10., 20., 5., 5., 15., 5., 0., 0., 0., 0. },
                    { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
                    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
  universe.addBonds({ 1, 7, 3, 8, 4, 9, 4, 10 }, { 7, 2, 8, 4, 9, 5, 10, 6 });

  pcm::MEHPNetwork baseNet = pcm::MEHPNetwork::fromUniverse(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, true, false, true);
  baseNet.strandContourLength = Eigen::VectorXd::Constant(4, 30.);
  baseNet.springContourLength = baseNet.strandContourLength;
  baseNet.recomputeMeanContourLengths();
  CHECK(baseNet.nrOfSprings == baseNet.nrOfStrands);
  CHECK(baseNet.nrOfSprings == 4);

  auto findLinkByOldAtomId = [](const pcm::MEHPNetwork& net, const int atomId) {
    for (size_t i = 0; i < net.nrOfLinks; ++i) {
      if (!net.linkIsEntanglement[static_cast<Eigen::Index>(i)] &&
          net.oldAtomIds[static_cast<Eigen::Index>(i)] == atomId) {
        return i;
      }
    }
    throw std::runtime_error("Link atom ID not found in test");
  };
  auto findStrandBetween =
    [](const pcm::MEHPNetwork& net, const size_t linkA, const size_t linkB) {
      for (size_t s = 0; s < net.nrOfStrands; ++s) {
        const size_t endA = net.linkIndicesOfStrand[s].front();
        const size_t endB = net.linkIndicesOfStrand[s].back();
        if ((endA == linkA && endB == linkB) ||
            (endA == linkB && endB == linkA)) {
          return s;
        }
      }
      throw std::runtime_error("Strand between links not found in test");
    };

  const size_t link1 = findLinkByOldAtomId(baseNet, 1);
  const size_t link2 = findLinkByOldAtomId(baseNet, 2);
  const size_t link3 = findLinkByOldAtomId(baseNet, 3);
  const size_t link4 = findLinkByOldAtomId(baseNet, 4);
  const size_t freeStrand = findStrandBetween(baseNet, link1, link2);
  const size_t swapStrand = findStrandBetween(baseNet, link3, link4);

  SECTION("MC condition accepts as requested")
  {
    pcm::MEHPNetwork net = baseNet;
    Eigen::VectorXd u =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    const double alphaOnSwap =
      (net.linkIndicesOfStrand[swapStrand].front() == link4) ? (1. - 0.966667)
                                                             : 0.966667;
    net.addSlipLinks({ freeStrand },
                     { swapStrand },
                     (Eigen::VectorXd(1) << 10.).finished(),
                     (Eigen::VectorXd(1) << 4.9).finished(),
                     Eigen::VectorXd::Zero(1),
                     { 0.5 },
                     { alphaOnSwap },
                     &u);
    CHECK(net.nrOfSprings != net.nrOfStrands);
    CHECK(net.nrOfStrands == 4);
    CHECK(net.nrOfSprings == 6);
    CHECK(net.nrOfCrosslinkSwapsEndured[6] == 0);
    CHECK(net.linkIndicesOfStrand[swapStrand].size() == 3);
    const size_t swapSpringIdx =
      (net.linkIndicesOfStrand[swapStrand].front() == link4)
        ? net.springIndicesOfStrand[swapStrand].front()
        : net.springIndicesOfStrand[swapStrand].back();
    CHECK_FALSE(net.moveSlipLinkAlongStrand(
      swapSpringIdx, u, pcm::LinkSwappingMode::ALL_MC_TRY, 1., 0, false));
    CHECK(net.moveSlipLinkAlongStrand(
      swapSpringIdx, u, pcm::LinkSwappingMode::ALL_MC_TRY, 1., 5, false));
    // check that the connectivity around the crosslink changed
    CHECK(net.linkIndicesOfStrand[swapStrand].size() == 2);
    CHECK(net.nrOfCrosslinkSwapsEndured[6] == 1);
  }

  SECTION("MC condition rejects as requested")
  {
    pcm::MEHPNetwork net = baseNet;
    Eigen::VectorXd u =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    const double alphaOnSwap =
      (net.linkIndicesOfStrand[swapStrand].front() == link4) ? (1. - 0.966667)
                                                             : 0.966667;
    net.addSlipLinks({ freeStrand },
                     { swapStrand },
                     (Eigen::VectorXd(1) << 10.).finished(),
                     (Eigen::VectorXd(1) << 4.9).finished(),
                     Eigen::VectorXd::Zero(1),
                     { 0.5 },
                     { alphaOnSwap },
                     &u);
    CHECK(net.nrOfSprings != net.nrOfStrands);
    CHECK(net.nrOfStrands == 4);
    CHECK(net.nrOfSprings == 6);
    const size_t freeSpringIdx = net.springIndicesOfStrand[freeStrand].back();
    CHECK_FALSE(net.moveSlipLinkAlongStrand(
      freeSpringIdx, u, pcm::LinkSwappingMode::ALL_MC_CYCLE, 1., -1, false));
    CHECK(net.nrOfCrosslinkSwapsEndured[6] == 0);
  }
}

TEST_CASE("MEHP Force Balance handles slip-links on primary loops",
          "[analysis][MEHPSimulator]")
{
  std::cout
    << "Running test \"MEHP Force Balance handles slip-links on primary loops\""
    << std::endl;
  const pe::Universe universe =
    pe::Universe(42.819955007276754, 42.819955007276754, 42.819955007276754);
  /**
   * Connectivity:
   *          9
   *         / |
   * 4-14-12-3-13
   *        |
   *        11-15-5
   *
   * 6-16-10-7
   */
  pe::Universe mutUniverse = universe;
  mutUniverse.addAtoms({ 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 },
                       { 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                       { 1.6205800871994722,
                         6.515162231365841,
                         12.1,
                         0.9,
                         5.9,
                         5,
                         3,
                         2,
                         6,
                         7,
                         11,
                         14,
                         15,
                         16 },
                       { 7.1412289503058295,
                         6.260972246279709,
                         12.1,
                         0.9,
                         5.9,
                         5,
                         3,
                         2,
                         6,
                         7,
                         11,
                         14,
                         15,
                         16 },
                       { 0.5796829850477182,
                         0.8320529182617298,
                         12.1,
                         0.9,
                         5.9,
                         5,
                         3,
                         2,
                         6,
                         7,
                         11,
                         14,
                         15,
                         16 },
                       { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                       { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                       { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
  mutUniverse.addBonds({ 3, 3, 3, 3, 4, 5, 6, 7, 9, 10, 11, 12 },
                       { 9, 11, 12, 13, 14, 15, 16, 10, 13, 16, 15, 14 });

  SECTION("Unentangled primary loop")
  {
    const pcm::MEHPNetwork net =
      pcm::MEHPNetwork::fromUniverse(mutUniverse,
                                     {},
                                     2,
                                     pcm::EntanglementMode::MERGED_LINKS,
                                     true,
                                     false,
                                     true);
    Eigen::VectorXd displacements =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    CHECK_NOTHROW(net.displaceToMeanPosition(displacements, 0, false, 1.));
    CHECK(displacements[0] == Catch::Approx(7.687).epsilon(1e-5));
    CHECK(displacements[1] == Catch::Approx(2.03926).epsilon(1e-5));
    CHECK(displacements[2] == Catch::Approx(5.88634).epsilon(1e-5));
  }

  SECTION("Entangled primary loop")
  {
    pcm::MEHPNetwork net =
      pcm::MEHPNetwork::fromUniverse(mutUniverse,
                                     {},
                                     2,
                                     pcm::EntanglementMode::MERGED_LINKS,
                                     true,
                                     false,
                                     true);
    net.strandContourLength = Eigen::VectorXd::Constant(
      static_cast<Eigen::Index>(net.nrOfStrands), 15.);
    net.springContourLength = net.strandContourLength;
    net.recomputeMeanContourLengths();
    Eigen::VectorXd displacements =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    size_t loopStrandIdx = 0;
    bool foundLoop = false;
    size_t detachedStrandIdx = 0;
    bool foundDetached = false;
    size_t strandTo5 = 0;
    bool foundStrandTo5 = false;
    for (size_t s = 0; s < net.nrOfStrands; ++s) {
      if (net.isLoopingStrand(s)) {
        loopStrandIdx = s;
        foundLoop = true;
      } else if (net.linkIndicesOfStrand[s].front() != 0 &&
                 net.linkIndicesOfStrand[s].back() != 0) {
        detachedStrandIdx = s;
        foundDetached = true;
      } else if (net.oldAtomIds[static_cast<Eigen::Index>(
                   net.linkIndicesOfStrand[s].front())] == 5 ||
                 net.oldAtomIds[static_cast<Eigen::Index>(
                   net.linkIndicesOfStrand[s].back())] == 5) {
        strandTo5 = s;
        foundStrandTo5 = true;
      }
    }
    REQUIRE(foundLoop);
    REQUIRE(foundDetached);
    REQUIRE(foundStrandTo5);
    const double alphaTo5 = (net.linkIndicesOfStrand[strandTo5].front() == 0)
                              ? 0.9310344827586207
                              : (1.0 - 0.9310344827586207);
    net.addSlipLinks(
      { loopStrandIdx, strandTo5 },
      { detachedStrandIdx, detachedStrandIdx },
      (Eigen::VectorXd(2) << 1.3263401618628183, 42.04664022316877).finished(),
      (Eigen::VectorXd(2) << 6.670300217824844, 7.18272624553976).finished(),
      (Eigen::VectorXd(2) << 41.85951429390015, 0.8578704100544575).finished(),
      { 0.5333333333333333, alphaTo5 },
      { 0.5, 0.5 },
      &displacements);
    CHECK_NOTHROW(net.displaceToMeanPosition(displacements, 0, false, 1.));
    CHECK(displacements[0] == Catch::Approx(0.187321).epsilon(1e-5));
    CHECK(displacements[1] == Catch::Approx(-0.447774).epsilon(1e-5));
    CHECK(displacements[2] == Catch::Approx(-0.925295).epsilon(1e-5));
  }
}

TEST_CASE("MEHPForceBalance phantom with and without removal is the same",
          "[analysis][MEHPSimulator]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  const std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  // a structure with lots of dangling things that can and will be entangled,
  // yet the entanglements removed
  const std::string inputFile =
    suspectedPath +
    "/structure/mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0."
    "0114_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  const pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPSimulator forceBalancerPhantom =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::MERGED_LINKS,
                       nullptr,
                       false,
                       true,
                       false,
                       false);
  pcm::MEHPSimulator forceBalancerPhantomRem =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::MERGED_LINKS,
                       nullptr,
                       false,
                       true,
                       false,
                       false);

  CHECK(forceBalancerPhantom.getNrOfSprings() ==
        forceBalancerPhantomRem.getNrOfSprings());
  CHECK_THAT(
    forceBalancerPhantom.getResidual(),
    Catch::Matchers::WithinRel(forceBalancerPhantomRem.getResidual(), 0.001));
  CHECK_THAT(forceBalancerPhantom.getSolubleWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getSolubleWeightFraction(), 0.001));
  CHECK_THAT(
    forceBalancerPhantom.getGamma(),
    Catch::Matchers::WithinRel(forceBalancerPhantomRem.getGamma(), 0.001));
  CHECK_THAT(forceBalancerPhantom.getDanglingWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getDanglingWeightFraction(), 0.001));
  CHECK((forceBalancerPhantom.getInitialCoordinates() +
         forceBalancerPhantom.getCurrentDisplacements())
          .isApprox(forceBalancerPhantom.getNetwork().coordinates +
                    forceBalancerPhantom.getCurrentDisplacements()));

  // run with and without simplification
  pcm::MEHPSlipLinkSolverConfig phantomConfig;
  phantomConfig.maxNrOfSteps = 10000;
  phantomConfig.xtol = 1e-14;
  phantomConfig.simplificationMode =
    pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
  forceBalancerPhantom.run(phantomConfig);

  pcm::MEHPSlipLinkSolverConfig phantomRemConfig;
  phantomRemConfig.maxNrOfSteps = 10000;
  phantomRemConfig.xtol = 1e-14;
  phantomRemConfig.simplificationMode =
    pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
  phantomRemConfig.inactiveRemovalCutoff = 1e-6;
  forceBalancerPhantomRem.run(phantomRemConfig);

  // compare results
  CHECK(forceBalancerPhantom.getNrOfSprings() >
        forceBalancerPhantomRem.getNrOfSprings());
  CHECK_THAT(
    forceBalancerPhantom.getResidual(),
    Catch::Matchers::WithinAbs(forceBalancerPhantomRem.getResidual(), 1e-4));
  CHECK_THAT(forceBalancerPhantom.getSolubleWeightFraction(),
             Catch::Matchers::WithinAbs(
               forceBalancerPhantomRem.getSolubleWeightFraction(), 5e-3));
  CHECK_THAT(forceBalancerPhantom.getActiveWeightFraction(),
             Catch::Matchers::WithinAbs(
               forceBalancerPhantomRem.getActiveWeightFraction(), 5e-3));
  CHECK_THAT(forceBalancerPhantom.getDanglingWeightFraction(),
             Catch::Matchers::WithinAbs(
               forceBalancerPhantomRem.getDanglingWeightFraction(), 5e-3));
  CHECK_THAT(forceBalancerPhantom.getGammaFactors(1.).sum(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getGammaFactors(1.).sum(), 1e-4));
}

TEST_CASE("MEHP Force Balance handles slip-link convergence correctly",
          "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance handles slip-link "
               "convergence correctly\""
            << std::endl;
  pe::Universe universe =
    pe::Universe(42.819955007276754, 42.819955007276754, 42.819955007276754);
  /**
   * Connectivity:
   *
   * 35-(11)-90
   *
   * 10-(12)-1654
   */
  // slip-link 3 in the test-system
  universe.addAtoms({ 35, 90, 10, 1654, 11, 12, 13, 14 },
                    { 2, 2, 2, 2, 1, 1, 1, 1 },
                    { 12.075848854154861,
                      10.644563425246883,
                      10,
                      14.302483570484272,
                      10,
                      10,
                      10,
                      10 },
                    { 3.574724359149917,
                      5.460837527830988,
                      10,
                      3.718195811318871,
                      10,
                      10,
                      10,
                      10 },
                    { 3.1018436428667284,
                      7.956714096296886,
                      10,
                      4.4007749635446824,
                      10,
                      10,
                      10,
                      10 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 });
  universe.addBonds({ 35, 11, 10, 12 }, { 11, 90, 12, 1654 });

  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::MERGED_LINKS,
                       nullptr,
                       false,
                       true,
                       false,
                       true);
  forceBalancer.addSlipLinks(
    { 1, 0 },
    { 1, 1 },
    (Eigen::VectorXd(2) << 12.650493316819828, 13.197029579176265).finished(),
    (Eigen::VectorXd(2) << 2.8706102036538566, 3.4016980009809297).finished(),
    (Eigen::VectorXd(2) << 8.475863644409664, 5.284899588057222).finished(),
    { 0.13793103448275862, 1. - 0.3103448275862069 },
    { 0.13793103448275862, 0.7931034482758621 },
    false);

  pcm::MEHPNetwork net = forceBalancer.getNetwork();
  Eigen::VectorXd displacements =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
  CHECK(net.springPartitions[2] ==
        Catch::Approx(0.7931034482758621 - 0.13793103448275862));
  // one output step
  net.updateSpringPartitionAtLink(5, displacements, 1.e10, false);
  CHECK(net.springPartitions[0] == Catch::Approx(0.367729));
  CHECK(net.springPartitions[4] == Catch::Approx(1. - 0.367729));
  CHECK(net.springPartitions[1] == Catch::Approx(0.13793103448275862));
  net.displaceToMeanPosition(displacements, 5, false, 1.e10);
  CHECK(displacements[5 * 3] == Catch::Approx(-0.391233));
  CHECK(displacements[5 * 3 + 1] == Catch::Approx(0.413665));
  CHECK(displacements[5 * 3 + 2] == Catch::Approx(0.028894));

  net = forceBalancer.getNetwork();
  displacements.setZero();
  // do many update steps
  for (int iter = 0; iter < 250; ++iter) {
    const double partRes =
      net.updateSpringPartitionAtLink(5, displacements, 1.e10, false);
    const double stepNorm =
      net.displaceToMeanPosition(displacements, 5, false, 1.e10);
    if (partRes < 1.e-10 && stepNorm < 1.e-10) {
      break;
    }
  }
  CHECK(displacements[5 * 3] == Catch::Approx(-0.592091));
  CHECK(displacements[5 * 3 + 1] == Catch::Approx(0.441203));
  CHECK(displacements[5 * 3 + 2] == Catch::Approx(0.44577));
  CHECK(net.springPartitions[1] == Catch::Approx(0.13793103448275862));

  SECTION("Stress tensor computations are equivalent")
  {
    const Eigen::Matrix3d stressTensor1 = forceBalancer.getStressTensor();
    const Eigen::Matrix3d stressTensor2 =
      forceBalancer.getNetwork().evaluateStressTensor(
        forceBalancer.getCurrentDisplacements(), 1.0, false, -1.0);
    for (Eigen::Index i = 0; i < 3; ++i) {
      for (Eigen::Index j = 0; j < 3; ++j) {
        CHECK(stressTensor1(i, j) == Catch::Approx(stressTensor2(i, j)));
      }
    }
  }
}

TEST_CASE("MEHPSimulator slip-link balance runs",
          "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHPSimulator slip-link balance runs\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Balance 3D case")
  {
    std::string largeInputFile =
      suspectedPath +
      "structure/xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      std::cout << "Reading file " << largeInputFile << std::endl;
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      std::cout << "Read file " << largeInputFile << std::endl;

      double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      pcm::MEHPSimulator forceBalancer2 =
        pcm::MEHPSimulator(universe2,
                           {},
                           2,
                           pcm::EntanglementMode::MERGED_LINKS,
                           nullptr,
                           false,
                           false,
                           false,
                           true);

      std::vector<size_t> allLinks(forceBalancer2.getNrOfLinks());
      std::iota(allLinks.begin(), allLinks.end(), size_t{ 0 });
      Eigen::Matrix3d stressTensor1 = forceBalancer2.getStressTensor();
      Eigen::Matrix3d stressTensor2 =
        forceBalancer2.getNetwork().evaluateStressTensorForLinks(
          allLinks, forceBalancer2.getCurrentDisplacements(), 1.0, false, 1.0);
      for (Eigen::Index i = 0; i < 3; ++i) {
        for (Eigen::Index j = 0; j < 3; ++j) {
          CHECK(stressTensor1(i, j) == Catch::Approx(stressTensor2(i, j)));
        }
      }

      SECTION("Actual balance results in correct phantom results")
      {
        std::cout << "Doing phantom force balance" << std::endl;
        pcm::MEHPSimulator forceRelaxer =
          pcm::MEHPSimulator(universe2,
                             {},
                             2,
                             pcm::EntanglementMode::NONE,
                             nullptr,
                             false,
                             false,
                             false,
                             true);
        CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
        CHECK(forceBalancer2.getNrOfIterations() == 0);
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(97.383096 * 97.383096 * 97.383096));
        // initial system values
        CHECK_THAT(
          forceBalancer2.getPressure(),
          Catch::Matchers::WithinAbs(forceRelaxer.getPressure(), 1e-5));
        CHECK(forceBalancer2.getPressure() == Catch::Approx(0.0061105865));
        CHECK_NOTHROW(forceBalancer2.run(pcm::MEHPSlipLinkSolverConfig{}));
        CHECK_NOTHROW(forceBalancer2.validateNetwork());
        CHECK(forceBalancer2.getNrOfSprings() == 9859);
        CHECK(forceBalancer2.getNrOfIterations() > 1);
        CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);

        // conversion factors
        double kb = 1.381e-23; // Boltzmann, J/K
        double T = 300.;       // Temperature, K
        double sigmaToNm = 0.616;
        double sigmaToM = sigmaToNm * 1.e-9;
        double slope = 0.00393 / (sigmaToNm * sigmaToNm); // sigma^2/(g/mol)
        double beadMass = 161.;                           // g/mol
        double Nb = 80.; // nr of beads per strand
        double conversionFactor =
          (forceBalancer2.getNetwork().meanSpringContourLength / Nb) * 3. * kb *
          T / (slope * beadMass); // J/sigma^2
        CHECK(conversionFactor / (sigmaToM * sigmaToM * 79.) ==
              Catch::Approx(0.0002450018));
        double nu =
          nrOfChains / (forceBalancer2.getVolume() * sigmaToM * sigmaToM *
                        sigmaToM); // chain number density, m^-3
        CHECK(nu == Catch::Approx(4.63241e25));

        // final values
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(0.153806 / 79.)); // LJ Units [?]
        CHECK(
          forceBalancer2.getPressure() * conversionFactor /
            (sigmaToM * sigmaToM * sigmaToM) ==
          Catch::Approx(61172.8878)); // shear modulus from the pressure, MPa
        double b02 = slope * beadMass;
        CHECK(forceBalancer2.getGammaFactor(1.) ==
              Catch::Approx(forceBalancer2.getGammaFactors(1.).mean()));
        CHECK(
          forceBalancer2.getGammaFactors(b02).mean() ==
          Catch::Approx((forceBalancer2.getGammaFactorsInDir(0, b02).mean() +
                         forceBalancer2.getGammaFactorsInDir(1, b02).mean() +
                         forceBalancer2.getGammaFactorsInDir(2, b02).mean())));
        CHECK(forceBalancer2.getGammaFactor(b02) ==
              Catch::Approx(forceBalancer2.getGammaFactors(b02).mean()));
        CHECK_THAT(
          forceBalancer2.getGammaFactor(b02) * kb * T * nu,
          Catch::Matchers::WithinRel(61308.3, 0.03)); // ANT shear modulus, Pa
        CHECK_THAT(forceBalancer2.getGammaFactor(b02),
                   Catch::Matchers::WithinRel(
                     forceBalancer2.getGammaFactor(b02, nrOfChains), 0.03));
        CHECK_THAT(forceBalancer2.getGammaFactor(b02, nrOfChains),
                   Catch::Matchers::WithinRel(0.319446, 0.03));
        CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
        CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
        CHECK(forceBalancer2.getNrOfActiveSprings() > 1);
        const pcm::MEHPNetwork& fbNet = forceBalancer2.getNetwork();
        const Eigen::VectorXd& fbU = forceBalancer2.getCurrentDisplacements();
        CHECK(forceBalancer2.getNrOfActiveStrands() <=
              static_cast<size_t>(
                fbNet.findActiveStrandsInDir(fbU, 0, 1e-3, false).count() +
                fbNet.findActiveStrandsInDir(fbU, 1, 1e-3, false).count() +
                fbNet.findActiveStrandsInDir(fbU, 2, 1e-3, false).count()));
        CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
        CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
              forceBalancer2.getNrOfNodes());
      }
      SECTION("Actual balance results in correct slip-link results")
      {
        std::cout << "Doing non-phantom force balance" << std::endl;
        pcm::MEHPSimulator forceRelaxer =
          pcm::MEHPSimulator(universe2,
                             {},
                             2,
                             pcm::EntanglementMode::NONE,
                             nullptr,
                             false,
                             false,
                             false,
                             true);
        CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
        CHECK(forceBalancer2.getNrOfIterations() == 0);
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(97.383096 * 97.383096 * 97.383096));
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(forceRelaxer.getPressure()));
        CHECK_THAT(forceBalancer2.getPressure(),
                   Catch::Matchers::WithinRel(0.0061105865, 1e-3));
        size_t nrOfSprings = forceRelaxer.getNetwork().nrOfSprings;
        forceBalancer2.addSlipLinks(
          { 10, 100, 50, 12, 76, 80, nrOfSprings - 1 },
          { 99, 101, 13, 7, 5, 19, nrOfSprings - 7 },
          Eigen::VectorXd::Zero(7),
          Eigen::VectorXd::Zero(7),
          Eigen::VectorXd::Zero(7));
        CHECK_NOTHROW(forceBalancer2.run(pcm::MEHPSlipLinkSolverConfig{}));
        CHECK(forceBalancer2.getPressure() == Catch::Approx(0.001955022));
      }
    } else {
      std::cout << "Skipping large file PDMS MEHP run" << std::endl;
      CHECK(true);
    }
  }

  SECTION("MEHP Force Balance 2D case")
  {
    CHECK(std::filesystem::exists(suspectedPath));
    universeSeq.initializeFromDataSequence(
      { { suspectedPath + "/structure/"
                          "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
                          "333_2d_t_7500001.structure.out" } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         false,
                         false,
                         true);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::UNSET);
    CHECK(forceBalancer.getNrOfIterations() == 0);
    CHECK(forceBalancer.getVolume() == Catch::Approx(universe.getVolume()));
    pcm::MEHPSlipLinkSolverConfig step5Config;
    step5Config.maxNrOfSteps = 5;
    forceBalancer.run(step5Config);
    CHECK(forceBalancer.getNrOfNodes() != universe.getNrOfAtoms());
    CHECK(forceBalancer.getNrOfIterations() <= 5);
    CHECK(forceBalancer.getNrOfIterations() >= 1);
    CHECK(universe.getAtomsOfType(2).size() == 7200);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::MAX_STEPS);
    CHECK_NOTHROW(forceBalancer.validateNetwork());

    pcm::MEHPSimulator forceBalancer2 =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         false,
                         false,
                         true);
    pcm::MEHPSlipLinkSolverConfig fullConfig;
    fullConfig.maxNrOfSteps = 10000;
    fullConfig.xtol = 1e-10;
    fullConfig.initialResidualToUse = 100.0;
    forceBalancer2.run(fullConfig);
    CHECK(forceBalancer2.getNrOfIterations() > 5);
    CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer2.getGammaFactor(1., forceBalancer2.getNrOfSprings()) ==
          Catch::Approx(1. / 3.).epsilon(0.001));
    Eigen::Matrix3d stressTensor = forceBalancer2.getStressTensor();
    CHECK(forceBalancer2.getPressure() ==
          Catch::Approx(
            (stressTensor(0, 0) + stressTensor(1, 1) + stressTensor(2, 2)) / 3.)
            .epsilon(0.02));
    CHECK_NOTHROW(forceBalancer2.validateNetwork());
    CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
    CHECK(forceBalancer2.getNrOfActiveSprings() > 1);
    CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
    CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
          forceBalancer2.getNrOfNodes());

    pe::Universe universe3 = forceBalancer2.getCrosslinkerVerse();
    CHECK(universe3.getNrOfAtoms() == forceBalancer2.getNrOfNodes());
    CHECK(universe3.getNrOfBonds() == forceBalancer2.getNrOfSprings());
    CHECK(universe3.getAtomsOfType(2).size() == universe3.getNrOfAtoms());
  }
}

TEST_CASE("MEHP Force Balance can adjust the network",
          "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance can adjust the network\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/mini_melt_3_a_10.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  CHECK(std::filesystem::exists(suspectedPath));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  CHECK(universeSeq.getLength() == 1);
  pe::Universe universe = universeSeq.atIndex(0);

  SECTION("Duplicate strand references detected by validate")
  {
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        5,
        5.,
        0.,
        3,
        0.,
        "test_seed",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        false,
        false,
        false);
    CHECK(forceBalancer.getNrOfStrands() == 3);
    pcm::MEHPNetwork net = forceBalancer.getNetwork();
    CHECK_NOTHROW(net.validate());
    const size_t nStrandsBefore = net.strandIndicesOfLink[0].size();
    net.strandIndicesOfLink[0].push_back(net.strandIndicesOfLink[0][0]);
    CHECK(net.strandIndicesOfLink[0].size() == nStrandsBefore + 1);
    CHECK_THROWS_AS(net.validate(), std::runtime_error);
  }

  SECTION("Primary loops")
  {
    // add a primary loop
    universe.addAtoms(
      { 69 }, { 1 }, { 0. }, { 0. }, { 0. }, { 0 }, { 0 }, { 0 });
    universe.addBonds({ 1, 69 }, { 69, 1 });

    pcm::MEHPSimulator withoutRemoval =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        0,
        5.,
        0.,
        0,
        0.,
        "test_seed",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        false,
        false,
        false);
    CHECK(withoutRemoval.getNrOfSprings() == 4);
    CHECK(withoutRemoval.getNumIntraChainSlipLinks() == 0);

    pcm::MEHPSimulator withRemoval =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        0,
        5.,
        0.,
        0,
        0.,
        "test_seed",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        false,
        true,
        false);
    CHECK(withRemoval.getNrOfSprings() == 3);
  }
}

TEST_CASE("MEHP Force Balance can break all springs",
          "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance can break all springs\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/mini_melt_3_a_10.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  CHECK(std::filesystem::exists(suspectedPath));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  CHECK(universeSeq.getLength() == 1);
  pe::Universe universe = universeSeq.atIndex(0);

  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe,
      5,
      5.,
      0.,
      3,
      0.,
      "test_seed",
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      true,
      false,
      nullptr,
      false,
      false,
      false,
      false);
  CHECK(forceBalancer.getNrOfStrands() == 3);
  CHECK(forceBalancer.breakTooLongStrands(1e5) == 0);
  CHECK(forceBalancer.breakTooLongStrands(1e-5) == 3);

  pe::Universe xlinkVerse = forceBalancer.getCrosslinkerVerse();
  CHECK(xlinkVerse.getNrOfBonds() == 0);
  std::unordered_map<long int, int> fs =
    forceBalancer.getEffectiveFunctionalityOfAtoms();
  CHECK(fs.size() == forceBalancer.getNrOfNodes());
  for (const std::pair<const long int, int>& pair : fs) {
    CHECK(pair.second == 0);
  }
}

TEST_CASE(
  "MEHP Force Balance can randomly add slip-links ignoring crosslinkers",
  "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance can randomly add slip-links "
               "ignoring crosslinkers\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  CHECK(std::filesystem::exists(suspectedPath));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  CHECK(universeSeq.getLength() == 1);
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file. " << std::endl;
  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       false,
                       true,
                       false,
                       false);
  size_t nrOfAddedLinks =
    forceBalancer.randomlyAddSlipLinks(250, 2.0, 12, 2.0, "seed1");
  CHECK(nrOfAddedLinks >= 12);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
  std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
  nrOfAddedLinks =
    forceBalancer.randomlyAddSlipLinks(250, 2.0, 25, 2.0, "seed2");
  CHECK(nrOfAddedLinks >= 25);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
  std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

  pcm::MEHPNetwork net = forceBalancer.getNetwork();
  Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
  CHECK(net.nrOfSprings > 0);
  size_t numRemoved = net.unlinkBifunctionalLinks(displacements);
  CHECK_NOTHROW(net.validate(displacements));
  CHECK(numRemoved > 0);
  // remove all springs...
  numRemoved = net.removeInactiveLinks(displacements, 1e5, false);
  CHECK(net.nrOfSprings == 0);
  CHECK(numRemoved > 0);
}

TEST_CASE("MEHP Force Balance can run with empty universe",
          "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance can run with empty universe\""
            << std::endl;
  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator(pe::Universe(),
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       false,
                       true,
                       false,
                       false);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
  CHECK_NOTHROW(forceBalancer.run(pcm::MEHPSlipLinkSolverConfig{}));
  CHECK_NOTHROW(forceBalancer.validateNetwork());
}

TEST_CASE("MEHP Force Balance can be interrupted", "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance can be interrupted\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  universeSeq.initializeFromDataSequence({ { inputFile } });
  CHECK(universeSeq.getLength() == 1);
  pe::Universe universe = universeSeq.atIndex(0);
  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe,
      1000,
      3.,
      0.,
      100,
      0.,
      "test_seed",
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      true,
      false,
      nullptr,
      false,
      true,
      false,
      true);
  pcm::MEHPSlipLinkSolverConfig config;
  config.maxNrOfSteps = 1000;
  config.xtol = 1e-12;
  config.simplificationMode =
    pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
  config.inactiveRemovalCutoff = 1e-3;
  config.doInnerIterations = true;
  config.linkSwappingMode = pcm::LinkSwappingMode::SLIPLINKS_ONLY;
  config.swappingFrequency = 1;
  config.oneOverSpringPartitionUpperLimit = 1.;
  config.nrOfCrosslinkSwapsAllowedPerSlipLink = 1;
  config.slipMode = pcm::SlipMode::SLIDING;
  pcm::MEHPSolveCallbacks cb;
  cb.shouldInterrupt = []() -> bool { return true; };
  cb.cleanupInterrupt = []() -> void {};
  CHECK_NOTHROW(forceBalancer.run(config, cb));
  CHECK(forceBalancer.getExitReason() == pcm::ExitReason::INTERRUPT);
}

TEST_CASE("MEHP Force Balance can run with swapping slip-links",
          "[analysis][MEHPSimulator][long]")
{
  std::cout
    << "Running test \"MEHP Force Balance can run with swapping slip-links\""
    << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         false,
                         false,
                         false);
    size_t nrOfAddedLinks =
      forceBalancer.randomlyAddSlipLinks(250, 2.0, 25, 2.0, "seed_swap");
    CHECK(nrOfAddedLinks >= 25);
    pe::Box oldBox = universe.getBox();
    pe::Box newBox =
      pe::Box(4 * oldBox.getLx(), 0.5 * oldBox.getLy(), 0.5 * oldBox.getLz());
    forceBalancer.deformTo(newBox);

    SECTION("Assuming too small box")
    {
      forceBalancer.configAssumeBoxLargeEnough(false);
      pcm::MEHPSlipLinkSolverConfig cfg1;
      cfg1.maxNrOfSteps = 100;
      cfg1.xtol = 1e-7;
      cfg1.simplificationMode =
        pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
      cfg1.inactiveRemovalCutoff = 0.01;
      cfg1.linkSwappingMode = pcm::LinkSwappingMode::ALL_MC;
      CHECK_NOTHROW(forceBalancer.run(cfg1));
      cfg1.xtol = 1e-8;
      cfg1.linkSwappingMode = pcm::LinkSwappingMode::SLIPLINKS_ONLY;
      CHECK_NOTHROW(forceBalancer.run(cfg1));
      cfg1.xtol = 1e-9;
      cfg1.linkSwappingMode = pcm::LinkSwappingMode::ALL;
      CHECK_NOTHROW(forceBalancer.run(cfg1));
    }

    SECTION("Assuming large enough box")
    {
      forceBalancer.configAssumeBoxLargeEnough(true);
      pcm::MEHPSlipLinkSolverConfig cfg2;
      cfg2.maxNrOfSteps = 100;
      cfg2.xtol = 1e-7;
      cfg2.simplificationMode =
        pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
      cfg2.inactiveRemovalCutoff = 0.01;
      cfg2.linkSwappingMode = pcm::LinkSwappingMode::ALL_MC;
      CHECK_NOTHROW(forceBalancer.run(cfg2));
      cfg2.xtol = 1e-8;
      cfg2.linkSwappingMode = pcm::LinkSwappingMode::SLIPLINKS_ONLY;
      CHECK_NOTHROW(forceBalancer.run(cfg2));
      cfg2.xtol = 1e-9;
      cfg2.linkSwappingMode = pcm::LinkSwappingMode::ALL;
      CHECK_NOTHROW(forceBalancer.run(cfg2));
    }
  }
}

TEST_CASE(
  "MEHP Force Balance can randomly add and remove slip-links with large box",
  "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHP Force Balance can randomly add and remove "
               "slip-links with large box\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         true);
    size_t nrOfAddedLinks =
      forceBalancer.randomlyAddSlipLinks(350, 2.0, 40, 2.0, "large_box_1");
    CHECK(nrOfAddedLinks >= 40);

    pcm::MEHPNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
    size_t numRemoved = net.unlinkBifunctionalLinks(displacements);
    CHECK_NOTHROW(net.validate(displacements));
    CHECK(numRemoved > 0);

    forceBalancer.setCurrentDisplacements(
      Eigen::VectorXd::Zero(3 * forceBalancer.getNrOfLinks()));
    pcm::MEHPSlipLinkSolverConfig slConfig;
    slConfig.maxNrOfSteps = 250;
    forceBalancer.run(slConfig);
    net = forceBalancer.getNetwork();
    displacements = forceBalancer.getCurrentDisplacements();
    size_t numInactiveRemoved =
      net.removeInactiveLinks(displacements, 0.1, false);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numInactiveRemoved > 0);
    CHECK_NOTHROW(net.validate(displacements));

    ////////////////////////////////////////////////////////////////
    forceBalancer = pcm::MEHPSimulator(universe,
                                       {},
                                       2,
                                       pcm::EntanglementMode::NONE,
                                       nullptr,
                                       false,
                                       false,
                                       true,
                                       true);
    nrOfAddedLinks =
      forceBalancer.randomlyAddSlipLinks(1000, 2.0, 100, 2.0, "large_box_2");
    CHECK(nrOfAddedLinks >= 100);

    slConfig.maxNrOfSteps = 50;
    forceBalancer.run(slConfig);
    net = forceBalancer.getNetwork();
    displacements = forceBalancer.getCurrentDisplacements();
    numInactiveRemoved = net.removeInactiveLinks(displacements, 0.1, false);
    CHECK(numInactiveRemoved > 0);
    numInactiveRemoved = net.unlinkBifunctionalLinks(displacements);
    // removeInactiveLinks already unlinked the resulting bifunctional links
    CHECK(numInactiveRemoved == 0);
    CHECK_NOTHROW(net.validate(displacements));
  }
}

TEST_CASE(
  "MEHP Force Balance can randomly add and remove slip-links with small box",
  "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance can randomly add and remove "
               "slip-links with small box\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         false);
    size_t nrOfAddedLinks =
      forceBalancer.randomlyAddSlipLinks(350, 2.0, 42, 2.0, "small_box_1");
    CHECK(nrOfAddedLinks >= 40);

    pcm::MEHPNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
    size_t numRemoved = net.unlinkBifunctionalLinks(displacements);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numRemoved > 0);

    pcm::MEHPSlipLinkSolverConfig slConfig;
    slConfig.maxNrOfSteps = 50;
    forceBalancer.run(slConfig);
    net = forceBalancer.getNetwork();
    displacements = forceBalancer.getCurrentDisplacements();
    size_t numInactiveRemoved =
      net.removeInactiveLinks(displacements, 0.1, false);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numInactiveRemoved > 0);
    CHECK_NOTHROW(net.validate(displacements));

    ////////////////////////////////////////////////////////////////
    forceBalancer = pcm::MEHPSimulator(universe,
                                       {},
                                       2,
                                       pcm::EntanglementMode::NONE,
                                       nullptr,
                                       false,
                                       false,
                                       true,
                                       false);
    nrOfAddedLinks =
      forceBalancer.randomlyAddSlipLinks(350, 2.0, 42, 2.0, "small_box_2");
    CHECK(nrOfAddedLinks >= 40);

    forceBalancer.run(slConfig);
    net = forceBalancer.getNetwork();
    displacements = forceBalancer.getCurrentDisplacements();
    numInactiveRemoved = net.removeInactiveLinks(displacements, 0.1, false);
    CHECK(numInactiveRemoved > 0);
    numInactiveRemoved = net.unlinkBifunctionalLinks(displacements);
    // removeInactiveLinks already unlinked the resulting bifunctional links
    CHECK(numInactiveRemoved == 0);
    CHECK_NOTHROW(net.validate(displacements));
  }
}

TEST_CASE("MEHP Force Balance handles slip-links", "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance handles slip-links\""
            << std::endl;
  pe::Universe universe = pe::Universe(1.0, 1.0, 1.0);
  universe.setBox(pe::Box(-15.0, 15.0, -15.0, 15.0, -15.0, 15.0));
  universe.addAtoms({ { 1, 2, 3, 4, 5, 6, 7, 8, 9 } },     // id
                    { { 2, 2, 2, 2, 1, 1, 1, 1, 1 } },     // type
                    { { -5, 5, 5, -5, 0, 5, -5, -5, 0 } }, // x
                    { { 5, 5, -5, -5, 5, 0, 0, 0, 0 } },   // y
                    { { 2, 2, 2, 2, 2, 2, 2, 2, 2 } },     // z
                    { { 1, 1, 1, 1, 1, 1, 1, 1, 1 } },     // nx
                    { { 1, 1, 1, 1, 1, 1, 1, 1, 1 } },     // ny
                    { { 1, 1, 1, 1, 1, 1, 1, 1, 1 } }      // nz
  );
  universe.addBonds(10,
                    { { 1, 2, 2, 2, 3, 3, 7, 8, 1, 9 } },
                    { { 5, 5, 9, 6, 6, 7, 4, 4, 8, 4 } },
                    { { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } },
                    false,
                    false);

  SECTION("Slip-links are placed where requested")
  {

    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         true);

    pcm::MEHPNetwork baseNet = forceBalancer.getNetwork();
    baseNet.strandContourLength.setConstant(15.0);
    baseNet.springContourLength.setConstant(15.0);
    baseNet.springPartitions.setConstant(1.0);
    baseNet.recomputeMeanContourLengths();
    forceBalancer = pcm::MEHPSimulator(baseNet);

    Eigen::VectorXd zeroX = Eigen::VectorXd::Zero(1);
    forceBalancer.addSlipLinks(
      { 0 }, { 2 }, zeroX, zeroX, zeroX, { 0.62 }, { 0.43 });

    Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
    CHECK(springPartitions[0] == Catch::Approx(0.62));
    CHECK(springPartitions[2] == Catch::Approx(0.43));
    CHECK(springPartitions[6] == Catch::Approx(1. - 0.43));
  }

  SECTION("Relevant slip-links behave correctly")
  {
    universe.addAtoms(
      { 10 }, { 1 }, { 0. }, { 0. }, { 0. }, { 1 }, { 1 }, { 1 });
    universe.addBonds(2, { { 1, 10 } }, { { 10, 3 } });
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         true);
    Eigen::VectorXd sX(1), sY(1), sZ(1);
    sX << 4.2;
    sY << 3.9;
    sZ << 1.2;
    CHECK_NOTHROW(forceBalancer.addSlipLinks({ 4 }, { 5 }, sX, sY, sZ));
    pcm::MEHPNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd displacements =
      Eigen::VectorXd::Zero(forceBalancer.getNrOfLinks() * 3);

    for (int i = 0; i < 5; ++i) {
      net.displaceToMeanPosition(displacements, 4, false, 1.);
      net.updateSpringPartitionAtLink(4, displacements, 1., false);
    }
    CHECK(net.springPartitions[net.springPartitions.size() - 1] ==
          Catch::Approx(0.5));
    CHECK(displacements[3 * 4] == Catch::Approx(-4.2));
    CHECK(displacements[3 * 4 + 1] == Catch::Approx(-3.9));
    CHECK(displacements[3 * 4 + 2] == Catch::Approx(0.8));
  }

  SECTION("Correct rationalisation and displacement")
  {
    // now, construct the force balancer
    pcm::MEHPNetwork initNet = pcm::MEHPNetwork::fromUniverse(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, true, false, false);
    const double N = 23.;
    initNet.strandContourLength.setConstant(N);
    initNet.springContourLength.setConstant(N);
    initNet.recomputeMeanContourLengths();
    pcm::MEHPSimulator forceBalancer2 = pcm::MEHPSimulator(initNet);
    CHECK(forceBalancer2.getNrOfNodes() == forceBalancer2.getNrOfLinks());
    CHECK(forceBalancer2.getNrOfNodes() == 4);
    CHECK(forceBalancer2.getNrOfSprings() == 5);
    pcm::MEHPNetwork net0 = forceBalancer2.getNetwork();
    Eigen::VectorXd springPartitions0 = net0.springPartitions;
    for (size_t i = 0; i < net0.nrOfSprings; i++) {
      std::cout << net0.springIndexA[static_cast<Eigen::Index>(i)] << ", "
                << net0.springIndexB[static_cast<Eigen::Index>(i)] << ": ";
      std::cout << springPartitions0[static_cast<Eigen::Index>(i)] << std::endl;
    }
    std::cout << std::endl;

    SECTION("Irrelevant slip-links are rationalised")
    {
      forceBalancer2.configAssumeBoxLargeEnough(true);
      // check throws when adding slip-links
      CHECK_THROWS(forceBalancer2.addSlipLinks({ 1, 2, 3 },
                                               { 1, 2, 3 },
                                               Eigen::Vector3d(1., 2., 3.),
                                               Eigen::Vector3d(1., 2., 3.),
                                               Eigen::Vector2d(1., 2.)));
      CHECK_THROWS(forceBalancer2.addSlipLinks({ 1, 2, 3 },
                                               { 1, 2 },
                                               Eigen::Vector3d(1., 2., 3.),
                                               Eigen::Vector3d(1., 2., 3.),
                                               Eigen::Vector3d(1., 2., 3.0)));
      // Resolve strands between crosslinkers (1-2, 2-3, 3-4)
      auto findLinkByOldId = [&forceBalancer2](const int atomId) {
        const pcm::MEHPNetwork& n = forceBalancer2.getNetwork();
        for (size_t i = 0; i < n.nrOfLinks; ++i) {
          if (!n.linkIsEntanglement[static_cast<Eigen::Index>(i)] &&
              n.oldAtomIds[static_cast<Eigen::Index>(i)] == atomId) {
            return i;
          }
        }
        throw std::runtime_error("Missing crosslink atom ID");
      };
      auto findStrandBetween = [&forceBalancer2](const size_t a,
                                                 const size_t b) {
        const pcm::MEHPNetwork& n = forceBalancer2.getNetwork();
        for (size_t s = 0; s < n.nrOfStrands; ++s) {
          const size_t f = n.linkIndicesOfStrand[s].front();
          const size_t t = n.linkIndicesOfStrand[s].back();
          if ((f == a && t == b) || (f == b && t == a)) {
            return s;
          }
        }
        throw std::runtime_error("Missing strand");
      };
      const size_t l1 = findLinkByOldId(1);
      const size_t l2 = findLinkByOldId(2);
      const size_t l3 = findLinkByOldId(3);
      const size_t l4 = findLinkByOldId(4);
      const size_t strand12 = findStrandBetween(l1, l2);
      const size_t strand23 = findStrandBetween(l2, l3);
      const size_t strand34 = findStrandBetween(l3, l4);

      // and actually add slip-links
      Eigen::VectorXd slX(2);
      slX << 4.2, 1.3;
      Eigen::VectorXd slY(2);
      slY << 3.9, 1.2;
      Eigen::VectorXd slZ(2);
      slZ << 1.3, 1.1;
      CHECK_NOTHROW(forceBalancer2.addSlipLinks(
        { strand12, strand34 }, { strand23, strand34 }, slX, slY, slZ));
      CHECK(forceBalancer2.getNrOfNodes() + 2 == forceBalancer2.getNrOfLinks());
      // and run with them.
      // Expect the slip-link of between two strands to converge to the central
      // atom Expect the slip-link of two strands to stay at 0.5, 0.5.
      Eigen::VectorXd displacements =
        Eigen::VectorXd::Zero(forceBalancer2.getNrOfLinks() * 3);
      pcm::MEHPNetwork net = forceBalancer2.getNetwork();
      // outputNetwork(net, displacements, net.springPartitions);

      net.displaceToMeanPosition(displacements, 4, false, 1.);
      CHECK(net.coordinates[4 * 3] + displacements[4 * 3] ==
            Catch::Approx(2.5));
      CHECK(net.coordinates[4 * 3 + 1] + displacements[4 * 3 + 1] ==
            Catch::Approx(2.5));
      CHECK(net.coordinates[4 * 3 + 2] + displacements[4 * 3 + 2] ==
            Catch::Approx(2));

      // reset
      displacements.setZero();

      CHECK(net.strandIndicesOfLink.size() == net.nrOfLinks);
      for (int j = 0; j < 125; ++j) {
        net.displaceToMeanPosition(displacements, 5, false, -1.);
        net.updateSpringPartitionAtLink(5, displacements, -1., false);
      }
      for (int i = 0; i < 125; ++i) {
        // do some random 250 steps with these two slip-links
        net.displaceToMeanPosition(displacements, 4, false, 1. * (i + 1));
        net.updateSpringPartitionAtLink(4, displacements, 1. * (i + 1), false);
      }
      // assert expectations are met.
      const Eigen::VectorXd& springPartitions = net.springPartitions;
      for (const size_t targetStrand : { strand12, strand23 }) {
        REQUIRE(net.springIndicesOfStrand[targetStrand].size() == 2);
        const size_t sp0 = net.springIndicesOfStrand[targetStrand][0];
        const size_t sp1 = net.springIndicesOfStrand[targetStrand][1];
        const size_t spAtL2 = net.isPartOfSpring(l2, sp0) ? sp0 : sp1;
        const size_t spOther = (spAtL2 == sp0) ? sp1 : sp0;
        CHECK_THAT(springPartitions[static_cast<Eigen::Index>(spAtL2)],
                   Catch::Matchers::WithinAbs(0., 1e-3));
        CHECK_THAT(springPartitions[static_cast<Eigen::Index>(spOther)],
                   Catch::Matchers::WithinAbs(1., 1e-3));
      }
      // CHECK(springPartitions[8] == Catch::Approx(1.0).margin(1e-6)); // 5-3
      // CHECK(springPartitions[7] + 1e-5 == Catch::Approx(0.0 +
      // 1e-5).margin(1e-6)); // 5-5 CHECK(springPartitions[3] ==
      // Catch::Approx(1.0).margin(1e-6)); // 5-0
    }
  }
}

TEST_CASE("MEHP Force Balance runs with non-network",
          "[analysis][MEHPSimulator][NonGaussianSpringForceEvaluator][long]")
{
  std::cout << "Running test \"MEHP Force Balance runs with non-network\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Balance 3D case")
  {
    std::string largeInputFile =
      suspectedPath +
      "/structure/xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      pcm::MEHPSimulator forceBalancer2 = pcm::MEHPSimulator(
        universe2, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      CHECK_NOTHROW(forceBalancer2.run(pcm::MEHPSlipLinkSolverConfig{}));
      CHECK(forceBalancer2.getNrOfIterations() > 1);
    }
  }
}

TEST_CASE("MEHP Force Balance Free chains collapse",
          "[analysis][MEHPSimulator]")
{
  size_t nrOfBeads = 30;
  size_t nrOfBeadsPerChain = 3;
  pe::Box originalBox =
    pe::Box(nrOfBeads * 10.0, nrOfBeads * 10.0, nrOfBeads * 10.0);
  pe::Universe universe = pe::Universe(originalBox);
  std::vector<double> xPositions;
  xPositions.reserve(nrOfBeads);
  std::vector<double> yPositions;
  yPositions.reserve(nrOfBeads);
  std::vector<double> zPositions;
  zPositions.reserve(nrOfBeads);
  std::vector<long int> atomIds;
  atomIds.reserve(nrOfBeads);
  std::vector<int> atomTypes;
  atomTypes.reserve(nrOfBeads);
  std::vector<int> zeroInts;
  zeroInts.reserve(nrOfBeads);
  std::vector<long int> bondFrom;
  bondFrom.reserve(nrOfBeads - 1);
  std::vector<long int> bondTo;
  bondTo.reserve(nrOfBeads - 1);
  double offset = 10.0;
  for (int i = 0; i < static_cast<int>(nrOfBeads); ++i) {
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i + 1);
    atomTypes.push_back(static_cast<size_t>(i) % nrOfBeadsPerChain == 0 ? 2
                                                                        : 1);
    zeroInts.push_back(0);
    if (i > 0) {
      bondFrom.push_back(i);
      bondTo.push_back(i + 1);
    }
  }
  universe.addAtoms(atomIds,
                    atomTypes,
                    xPositions,
                    yPositions,
                    zPositions,
                    zeroInts,
                    zeroInts,
                    zeroInts);
  universe.addBonds(bondFrom, bondTo);
  CHECK(universe.getNrOfAtoms() == nrOfBeads);
  CHECK(universe.getNrOfBonds() == nrOfBeads - 1);

  pe::Box deformedBox =
    pe::Box(nrOfBeads * 9.0, nrOfBeads * 9.0, nrOfBeads * 9.0);
  universe.setBox(deformedBox, true);

  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       false,
                       true,
                       false,
                       false);
  forceBalancer.deformTo(originalBox);
  CHECK_THAT(forceBalancer.getVolume(),
             Catch::Matchers::WithinRel(originalBox.getVolume()));
  CHECK(forceBalancer.getNrOfSprings() == forceBalancer.getNrOfStrands());
  CHECK(forceBalancer.getNumBonds() == forceBalancer.getNrOfSprings());
  CHECK(forceBalancer.getNumBondsToForm() == 0);
  CHECK(forceBalancer.getNumAtoms() == forceBalancer.getNrOfNodes());
  CHECK(forceBalancer.getNumExtraBonds() == 0);
  CHECK(forceBalancer.getNrOfSprings() == nrOfBeads / nrOfBeadsPerChain);
  CHECK(forceBalancer.getNumShifts() == 0);
  CHECK(forceBalancer.getTimestep() == 1);
  CHECK(forceBalancer.getCurrentTime(12345) == 0.0);
  CHECK(forceBalancer.getNumRelocations() == 0);
  CHECK(forceBalancer.getCoordinates().isApprox(
    forceBalancer.getNetwork().coordinates));
  CHECK(forceBalancer.getCurrentDisplacements().isZero());

  SECTION("Large enough box")
  {
    forceBalancer.configAssumeBoxLargeEnough(true);
    pcm::MEHPSlipLinkSolverConfig cfg;
    cfg.maxNrOfSteps = 50000;
    cfg.xtol = 1e-12;
    CHECK_NOTHROW(forceBalancer.run(cfg));
    CHECK_FALSE(forceBalancer.getCurrentDisplacements().isZero());
    CHECK(forceBalancer.getCoordinates().isApprox(
      forceBalancer.getNetwork().coordinates +
      forceBalancer.getCurrentDisplacements()));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
    CHECK(forceBalancer.getActiveWeightFraction() == 0.);
    CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
               Catch::Matchers::WithinRel(1.));
    CHECK_THAT(forceBalancer.getDanglingWeightFraction(),
               Catch::Matchers::WithinAbs(0., 1e-9));
    CHECK(forceBalancer.getAverageSpringLength() >= 0.0);
    CHECK(forceBalancer.getAverageSpringLength() <= 3e-6);
    CHECK(forceBalancer.getResiduals().isZero(1e-6));
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }

  SECTION("Not large enough box")
  {
    forceBalancer.configAssumeBoxLargeEnough(false);
    pcm::MEHPSlipLinkSolverConfig cfg;
    cfg.maxNrOfSteps = 50000;
    cfg.xtol = 1e-12;
    CHECK_NOTHROW(forceBalancer.run(cfg));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
    CHECK(forceBalancer.getActiveWeightFraction() == 0.);
    CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
               Catch::Matchers::WithinRel(1.));
    CHECK_THAT(forceBalancer.getDanglingWeightFraction(),
               Catch::Matchers::WithinAbs(0., 1e-9));
    CHECK(forceBalancer.getAverageSpringLength() >= 0.0);
    CHECK(forceBalancer.getAverageSpringLength() <= 3e-6);
    CHECK(forceBalancer.getResiduals().isZero(1e-6));
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }
}

TEST_CASE("MEHP Force Balance Entanglement Beads Are Removed",
          "[analysis][MEHPSimulator]")
{
  size_t nrOfBeads = 30;
  size_t nrOfBeadsPerChain = 3;
  pe::Universe universe =
    pe::Universe(nrOfBeads * 10.0, nrOfBeads * 10.0, nrOfBeads * 10.0);
  std::vector<double> xPositions;
  std::vector<double> yPositions;
  std::vector<double> zPositions;
  std::vector<long int> atomIds;
  std::vector<int> atomTypes;
  std::vector<int> zeroInts;
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  std::vector<int> bondTypes;
  double offset = 10.0;
  for (int i = 0; i < static_cast<int>(nrOfBeads); ++i) {
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) + offset);
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset);
    atomIds.push_back(i + 1);
    atomTypes.push_back(static_cast<size_t>(i) % nrOfBeadsPerChain == 0 ? 2
                                                                        : 1);
    zeroInts.push_back(0);
    if (i > 0) {
      bondFrom.push_back(i);
      bondTo.push_back(i + 1);
      bondTypes.push_back(1);
    }
  }
  CHECK(atomTypes[13] == 1);
  CHECK(atomTypes[17] == 1);
  atomTypes[13] = 3;
  atomTypes[17] = 3;
  bondFrom.push_back(13 + 1);
  bondTo.push_back(17 + 1);
  bondTypes.push_back(3);
  universe.addAtoms(atomIds,
                    atomTypes,
                    xPositions,
                    yPositions,
                    zPositions,
                    zeroInts,
                    zeroInts,
                    zeroInts);
  universe.addBonds(bondFrom, bondTo, bondTypes);
  CHECK(universe.getNrOfAtoms() == nrOfBeads);
  CHECK(universe.getNrOfBonds() == bondTo.size());

  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       false,
                       true,
                       false,
                       false);
  CHECK(forceBalancer.getNrOfSprings() == forceBalancer.getNrOfStrands());
  CHECK(forceBalancer.getNrOfSprings() == nrOfBeads / nrOfBeadsPerChain + 3);
  double initialResidual = forceBalancer.getResiduals().squaredNorm();

  pcm::MEHPSlipLinkSolverConfig cfg;
  cfg.maxNrOfSteps = 20000;
  cfg.xtol = 1e-12;
  cfg.initialResidualToUse = initialResidual;
  cfg.simplificationMode = pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;

  SECTION("Large enough box")
  {
    forceBalancer.configAssumeBoxLargeEnough(true);
    CHECK_NOTHROW(forceBalancer.run(cfg));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
    CHECK(forceBalancer.getNrOfSprings() == 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }

  SECTION("Not large enough box")
  {
    forceBalancer.configAssumeBoxLargeEnough(false);
    CHECK_NOTHROW(forceBalancer.run(cfg));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
    CHECK(forceBalancer.getNrOfSprings() == 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }
}

TEST_CASE("MEHP Force Balance fully active chains are fully active",
          "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHP Force Balance fully active chains are "
               "fully active\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  SECTION("MEHP Force Balance 3D case")
  {
    // perfect diamond network = fully connected =>
    // maximum is at perfect crystal structure -> must be all active.
    std::string inputFile =
      suspectedPath +
      "3d-diamond-lattice_10x10x10_a_3_d_0.85_v_0.V-fixed.structure.out";
    if (std::filesystem::exists(inputFile)) {
      std::cout << "Reading file " << inputFile << std::endl;
      universeSeq.initializeFromDataSequence({ { inputFile } });
      pe::Universe universe = universeSeq.atIndex(0);
      std::cout << "Read file " << inputFile << std::endl;

      std::vector<pe::Molecule> molecules =
        universe.getChainsWithCrosslinker(2);
      for (pe::Molecule mol : molecules) {
        CHECK(mol.getType() == pe::NETWORK_STRAND);
      }

      SECTION("Without entanglements")
      {
        pcm::MEHPSimulator forceBalancer =
          pcm::MEHPSimulator(universe,
                             {},
                             2,
                             pcm::EntanglementMode::NONE,
                             nullptr,
                             false,
                             true,
                             false,
                             false);
        size_t initialNSprings = forceBalancer.getNrOfSprings();

        SECTION("For large enough box")
        {
          forceBalancer.configAssumeBoxLargeEnough(true);
          CHECK(forceBalancer.getNrOfActiveSprings() ==
                forceBalancer.getNrOfSprings());
          double initialResidual = forceBalancer.getResiduals().squaredNorm();
          CHECK(std::isfinite(initialResidual));
          pcm::MEHPSlipLinkSolverConfig cfg;
          cfg.maxNrOfSteps = 10000;
          cfg.xtol = 1e-10;
          cfg.initialResidualToUse = initialResidual;
          cfg.simplificationMode =
            pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
          CHECK_NOTHROW(forceBalancer.run(cfg));
          CHECK(forceBalancer.getNrOfIterations() > 0);
          CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
          CHECK(forceBalancer.getNrOfActiveSprings() ==
                forceBalancer.getNrOfSprings());
          CHECK(forceBalancer.getNrOfSprings() == initialNSprings);
          CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                     Catch::Matchers::WithinRel(1.0));
          CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                     Catch::Matchers::WithinAbs(0.0, 1e-9));
          CHECK(initialResidual > forceBalancer.getResiduals().squaredNorm());
        }

        SECTION("For not large enough box")
        {
          forceBalancer.configAssumeBoxLargeEnough(false);
          CHECK(forceBalancer.getNrOfActiveSprings() ==
                forceBalancer.getNrOfSprings());
          double initialResidual = forceBalancer.getResiduals().squaredNorm();
          CHECK(std::isfinite(initialResidual));
          pcm::MEHPSlipLinkSolverConfig cfg;
          cfg.maxNrOfSteps = 10000;
          cfg.xtol = 1e-10;
          cfg.initialResidualToUse = initialResidual;
          cfg.simplificationMode =
            pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
          CHECK_NOTHROW(forceBalancer.run(cfg));
          CHECK(forceBalancer.getNrOfIterations() > 0);
          CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
          CHECK(forceBalancer.getNrOfActiveSprings() ==
                forceBalancer.getNrOfSprings());
          CHECK(forceBalancer.getNrOfSprings() == initialNSprings);
          CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                     Catch::Matchers::WithinRel(1.0));
          CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                     Catch::Matchers::WithinAbs(0.0, 1e-9));
          CHECK(initialResidual > forceBalancer.getResiduals().squaredNorm());
        }
      }

      SECTION("With entanglements")
      {
        pcm::MEHPSimulator forceBalancer =
          pcm::MEHPSimulator::constructWithRandomEntanglements(
            universe,
            400,
            2.0,
            0.0,
            100,
            0.0,
            "a533d",
            2,
            pcm::EntanglementMode::MERGED_LINKS,
            true,
            false,
            nullptr,
            false,
            true,
            false,
            false);
        CHECK(forceBalancer.getNrOfActiveSprings() ==
              forceBalancer.getNrOfSprings());
        double initialResidual = forceBalancer.getResiduals().squaredNorm();
        CHECK(std::isfinite(initialResidual));
        pcm::MEHPSlipLinkSolverConfig cfg;
        cfg.maxNrOfSteps = 10000;
        cfg.xtol = 1e-10;
        cfg.initialResidualToUse = initialResidual;
        cfg.simplificationMode =
          pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
        CHECK_NOTHROW(forceBalancer.run(cfg));
        CHECK(forceBalancer.getNrOfIterations() > 0);
        CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
        CHECK(forceBalancer.getNrOfActiveSprings() ==
              forceBalancer.getNrOfSprings());
        CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                   Catch::Matchers::WithinRel(1.0));
        CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                   Catch::Matchers::WithinAbs(0.0, 1e-9));
        CHECK(initialResidual > forceBalancer.getResiduals().squaredNorm());
      }
    }
  }
}

TEST_CASE("MEHP Force Balance Gives Identical Results for Different PBC "
          "imperfect Diamond Network",
          "[analysis][MEHPSimulator]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath +
    "3d-diamond-lattice_5x5x5_a_3_d_0.85_imperfect.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPSimulator forceBalanceConventional =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        250,
        2.0,
        0.0,
        100,
        2.0,
        "a533d",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        true,
        false,
        true);

    pcm::MEHPSimulator forceBalanceNew =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        250,
        2.0,
        0.0,
        100,
        2.0,
        "a533d",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        true,
        false,
        false);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));

    double initialResidual =
      forceBalanceConventional.getResiduals().squaredNorm();
    CHECK(initialResidual ==
          Catch::Approx(forceBalanceNew.getResiduals().squaredNorm()));

    pcm::MEHPSlipLinkSolverConfig cfg1;
    cfg1.maxNrOfSteps = 250;
    cfg1.xtol = 1e-4;
    forceBalanceConventional.run(cfg1);
    forceBalanceNew.run(cfg1);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getResidualNorm() ==
          Catch::Approx(forceBalanceNew.getResidualNorm()));
    CHECK_THAT(
      forceBalanceConventional.getIndicesOfActiveNodes(),
      Catch::Matchers::Equals(forceBalanceNew.getIndicesOfActiveNodes()));
    CHECK(forceBalanceConventional.getNrOfSprings() ==
          forceBalanceNew.getNrOfSprings());
    CHECK(forceBalanceConventional.getCurrentDisplacements().isApprox(
      forceBalanceNew.getCurrentDisplacements()));

    pcm::MEHPSlipLinkSolverConfig cfg2;
    cfg2.maxNrOfSteps = 250;
    cfg2.xtol = 1e-8;
    cfg2.initialResidualToUse = initialResidual;
    cfg2.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
    cfg2.inactiveRemovalCutoff = 0.01;
    cfg2.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
    cfg2.swappingFrequency = 10;
    cfg2.oneOverSpringPartitionUpperLimit = 1.;
    cfg2.nrOfCrosslinkSwapsAllowedPerSlipLink = 1;
    forceBalanceConventional.run(cfg2);
    forceBalanceNew.run(cfg2);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getResidualNorm() ==
          Catch::Approx(forceBalanceNew.getResidualNorm()));
    CHECK(forceBalanceConventional.getNrOfSprings() ==
          forceBalanceNew.getNrOfSprings());
    CHECK_THAT(
      forceBalanceConventional.getIndicesOfActiveNodes(),
      Catch::Matchers::Equals(forceBalanceNew.getIndicesOfActiveNodes()));
    CHECK(forceBalanceConventional.getCurrentDisplacements().isApprox(
      forceBalanceNew.getCurrentDisplacements()));

    pcm::MEHPSlipLinkSolverConfig cfg3;
    cfg3.maxNrOfSteps = 250;
    cfg3.xtol = 1e-8;
    cfg3.initialResidualToUse = initialResidual;
    cfg3.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg3.inactiveRemovalCutoff = 1.0;
    cfg3.linkSwappingMode = pcm::LinkSwappingMode::ALL_MC;
    cfg3.swappingFrequency = 10;
    cfg3.oneOverSpringPartitionUpperLimit = 1.;
    cfg3.nrOfCrosslinkSwapsAllowedPerSlipLink = 0;
    forceBalanceConventional.run(cfg3);
    forceBalanceNew.run(cfg3);

    CHECK(forceBalanceConventional.getNrOfSprings() ==
          forceBalanceNew.getNrOfSprings());
    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getResidualNorm() ==
          Catch::Approx(forceBalanceNew.getResidualNorm()));
    CHECK(forceBalanceConventional.getCurrentDisplacements().isApprox(
      forceBalanceNew.getCurrentDisplacements()));
    CHECK_THAT(
      forceBalanceConventional.getIndicesOfActiveNodes(),
      Catch::Matchers::Equals(forceBalanceNew.getIndicesOfActiveNodes()));
  }
}

TEST_CASE("MEHP Force Balance Gives Identical Results for Different PBC p = 1",
          "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHP Force Balance Gives Identical Results for "
               "Different PBC p = 1\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPSimulator forceBalanceConventional =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        25,
        2.0,
        0.0,
        20,
        2.0,
        "a533d",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        true,
        false,
        true);

    pcm::MEHPSimulator forceBalanceNew =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        25,
        2.0,
        0.0,
        20,
        2.0,
        "a533d",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        true,
        false,
        false);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));

    double initialResidual =
      forceBalanceConventional.getResiduals().squaredNorm();
    CHECK(initialResidual ==
          Catch::Approx(forceBalanceNew.getResiduals().squaredNorm()));

    pcm::MEHPSlipLinkSolverConfig cfg1;
    cfg1.maxNrOfSteps = 250;
    cfg1.xtol = 1e-4;
    forceBalanceConventional.run(cfg1);
    forceBalanceNew.run(cfg1);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getResidualNorm() ==
          Catch::Approx(forceBalanceNew.getResidualNorm()));
    CHECK_THAT(
      forceBalanceConventional.getIndicesOfActiveNodes(),
      Catch::Matchers::Equals(forceBalanceNew.getIndicesOfActiveNodes()));

    pcm::MEHPSlipLinkSolverConfig cfg2;
    cfg2.maxNrOfSteps = 250;
    cfg2.xtol = 1e-8;
    cfg2.initialResidualToUse = initialResidual;
    cfg2.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
    cfg2.inactiveRemovalCutoff = 0.01;
    cfg2.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
    cfg2.swappingFrequency = 10;
    cfg2.oneOverSpringPartitionUpperLimit = 1.;
    cfg2.nrOfCrosslinkSwapsAllowedPerSlipLink = 1;
    forceBalanceConventional.run(cfg2);
    forceBalanceNew.run(cfg2);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getResidualNorm() ==
          Catch::Approx(forceBalanceNew.getResidualNorm()));
    CHECK_THAT(
      forceBalanceConventional.getIndicesOfActiveNodes(),
      Catch::Matchers::Equals(forceBalanceNew.getIndicesOfActiveNodes()));

    pcm::MEHPSlipLinkSolverConfig cfg3 = cfg2;
    cfg3.linkSwappingMode = pcm::LinkSwappingMode::ALL_MC;
    cfg3.nrOfCrosslinkSwapsAllowedPerSlipLink = 3;
    forceBalanceConventional.run(cfg3);
    forceBalanceNew.run(cfg3);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getResidualNorm() ==
          Catch::Approx(forceBalanceNew.getResidualNorm()));
    CHECK_THAT(
      forceBalanceConventional.getIndicesOfActiveNodes(),
      Catch::Matchers::Equals(forceBalanceNew.getIndicesOfActiveNodes()));
  } else {
    std::cout << "File " << inputFile << " does not exists. Skipping test."
              << std::endl;
  }
}

TEST_CASE("MEHP Force Balance does not collapse", "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance does not collapse\""
            << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  universe.addAtoms(
    { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } },
    { { 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1 } },
    { { 0.,
        2.5,
        5,
        7.5,
        0.1,
        2.5,
        5,
        7.5,
        -0.1,
        5.,
        0.,
        5. } }, // x with slight (0.1) deviation, so we don't start perfect
    { { 0.1, 0., -0.1, 0., 5., 5., 5., 5., 2.5, 2.5, 7.5, 7.5 } }, // y
    { { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. } },        // z
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } });
  universe.addBonds(
    { { 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7 } },
    { { 2, 9, 4, 11, 2, 4, 10, 12, 9, 11, 6, 8, 6, 8, 10, 12 } });

  pcm::MEHPSlipLinkSolverConfig cfg;
  cfg.maxNrOfSteps = 10000;
  cfg.xtol = 1e-15;
  cfg.simplificationMode = pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;

  SECTION("Assuming large enough box")
  {
    pcm::MEHPSimulator forceBalanceConventional =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         true,
                         false,
                         true);
    CHECK_NOTHROW(forceBalanceConventional.run(cfg));
    CHECK(forceBalanceConventional.getNrOfIterations() > 0);
    CHECK(forceBalanceConventional.getExitReason() ==
          pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalanceConventional.getNrOfActiveSprings() ==
          forceBalanceConventional.getNrOfSprings());
    CHECK(forceBalanceConventional.getNrOfActiveSprings() == 8);
    CHECK(forceBalanceConventional.getNrOfActiveNodes() == 4);
    CHECK_THAT(forceBalanceConventional.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.0, 1e-9));
  }
  SECTION("Assuming too small box")
  {
    pcm::MEHPSimulator forceBalanceNew =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         true,
                         false,
                         false);
    CHECK_NOTHROW(forceBalanceNew.run(cfg));
    CHECK(forceBalanceNew.getNrOfIterations() > 0);
    CHECK(forceBalanceNew.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalanceNew.getNrOfActiveSprings() ==
          forceBalanceNew.getNrOfSprings());
    CHECK(forceBalanceNew.getNrOfActiveSprings() == 8);

    const pcm::MEHPNetwork& fbNet = forceBalanceNew.getNetwork();
    const Eigen::VectorXd& fbU = forceBalanceNew.getCurrentDisplacements();
    CHECK(forceBalanceNew.getNrOfActiveStrands() <=
          static_cast<size_t>(
            fbNet.findActiveStrandsInDir(fbU, 0, 1e-3, true).count() +
            fbNet.findActiveStrandsInDir(fbU, 1, 1e-3, true).count() +
            fbNet.findActiveStrandsInDir(fbU, 2, 1e-3, true).count()));
    CHECK(forceBalanceNew.getNrOfActiveNodes() == 4);
    CHECK(forceBalanceNew.getAverageSpringLength() == Catch::Approx(5.0));
    CHECK_THAT(forceBalanceNew.getGammaFactor(1., 100),
               Catch::Matchers::WithinAbs(1.0, 1e-2));
    CHECK_THAT(forceBalanceNew.getPressure(),
               Catch::Matchers::WithinAbs(1. / 30., 1e-3));
    CHECK_THAT(forceBalanceNew.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.0, 1e-9));
  }
};

TEST_CASE(
  "MEHP Force Balance correctly collapses melts even with random slip-links",
  "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance correctly collapses melts "
               "even with random slip-links\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile = suspectedPath + "melt_83_a_100.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.setDataFileAtomStyle(
      { pu::AtomStyle::HYBRID, pu::AtomStyle::BOND, pu::AtomStyle::EDPD });
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        1000,
        2.0,
        0.0,
        100,
        5,
        "my_seed_fb12",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        true,
        false,
        false);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    double initialResidual = forceBalancer.getResidualNorm();
    CHECK(std::isfinite(initialResidual));
    CHECK(forceBalancer.getNumExtraAtoms() > 100);
    CHECK(forceBalancer.getNrOfSprings() > forceBalancer.getNrOfStrands());
    pcm::MEHPSlipLinkSolverConfig cfg;
    cfg.maxNrOfSteps = 50000;
    cfg.xtol = 1e-9;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
    CHECK_NOTHROW(forceBalancer.run(cfg));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
    const pcm::MEHPNetwork& fbNet = forceBalancer.getNetwork();
    const Eigen::VectorXd& fbU = forceBalancer.getCurrentDisplacements();
    CHECK(fbNet.findActiveStrandsInDir(fbU, 0, 1e-3, false).count() == 0);
    CHECK(fbNet.findActiveStrandsInDir(fbU, 1, 1e-3, false).count() == 0);
    CHECK(fbNet.findActiveStrandsInDir(fbU, 2, 1e-3, false).count() == 0);
    CHECK(initialResidual > forceBalancer.getResidualNorm());

    double pressBefore = forceBalancer.getPressure();
    forceBalancer.configAssumeBoxLargeEnough(true);
    CHECK(pressBefore >= forceBalancer.getPressure());
  }
}

TEST_CASE("MEHP Force Balance correctly re-aligns Slip-Links to Images",
          "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance correctly re-aligns "
               "Slip-Links to Images\""
            << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(universe);
  forceBalancer.configAssumeBoxLargeEnough(false);
  CHECK(forceBalancer.getPressure() == Catch::Approx(0.0));

  pcm::MEHPNetwork net;
  net.box = pe::Box(10., 10., 10.);
  net.L = { 10., 10., 10. };
  net.boxHalfs = { 5., 5., 5. };
  net.vol = 1000.;
  net.nrOfLinks = 3;
  net.nrOfNodes = 2;
  net.nrOfStrands = 1;
  net.nrOfSprings = 2;
  net.coordinates = Eigen::VectorXd::Zero(3 * 3);
  net.springIndexA = Eigen::ArrayXi::Zero(2);
  net.springIndexB = Eigen::ArrayXi::Zero(2);
  net.springCoordinateIndexA = Eigen::ArrayXi::Zero(2);
  net.springCoordinateIndexB = Eigen::ArrayXi::Zero(2);
  net.springBoxOffset = Eigen::VectorXd::Zero(2 * 3);
  net.strandIndexOfSpring = Eigen::ArrayXi::Zero(2);
  net.linkIsEntanglement = Eigen::ArrayXb::Constant(3, false);

  SECTION("System 1")
  {
    net.springIndexA << 0, 2;
    net.springIndexB << 2, 1;
    net.linkIsEntanglement[2] = true;

    Eigen::VectorXd u = Eigen::VectorXd::Zero(net.coordinates.size());

    CHECK_NOTHROW(net.realignSlipLinkToImages(u, 2, 0, 1));
    CHECK((net.springBoxOffset.array() == (0.0)).all());

    net.coordinates << 5, 2.5, 0., // crosslink 1,
      5., 7.5, 0.,                 // crosslink 2,
      5., 5., 0.;                  // slip-link

    net.springBoxOffset << 0., 0., 0., // crosslink 1 -> slip-link
      20., 0., 0.;                     // slip-link -> crosslink 2

    CHECK_NOTHROW(net.realignSlipLinkToImages(u, 2, 0, 1));
    CHECK(net.coordinates.segment(3 * 2, 3)[1] == 5.);

    CHECK(net.springBoxOffset[0] == Catch::Approx(10.));
    CHECK(net.springBoxOffset[3] == Catch::Approx(10.));
    for (Eigen::Index i = 0; i < 2 * 3; ++i) {
      if (i % 3 != 0) {
        CHECK(net.springBoxOffset[i] == Catch::Approx(0.));
      }
    }
  }

  SECTION("System 2")
  {
    net.springIndexA << 1, 2;
    net.springIndexB << 2, 0;
    net.linkIsEntanglement[2] = true;

    Eigen::VectorXd u = Eigen::VectorXd::Zero(net.coordinates.size());

    CHECK_NOTHROW(net.realignSlipLinkToImages(u, 2, 0, 1));
    CHECK((net.springBoxOffset.array() == (0.0)).all());

    net.coordinates << 5, 2.5, 0., // crosslink 1,
      5., 7.5, 0.,                 // crosslink 2,
      5., 7.5, 0.;                 // slip-link

    net.springBoxOffset << 20., 0., 0., // crosslink 2 -> slip-link
      0., 0., 0.;                       // slip-link -> crosslink 1

    CHECK_NOTHROW(net.realignSlipLinkToImages(u, 2, 0, 1));

    Eigen::VectorXd expectation = Eigen::VectorXd::Zero(2 * 3);
    expectation(0) = 10.;
    expectation(3) = 10.;

    for (Eigen::Index i = 0; i < 2 * 3; ++i) {
      CHECK(expectation[i] == Catch::Approx(net.springBoxOffset[i]));
    }
  }

  SECTION("System 3")
  {
    net.springIndexA << 2, 0;
    net.springIndexB << 1, 2;
    net.linkIsEntanglement[2] = true;

    Eigen::VectorXd u = Eigen::VectorXd::Zero(net.coordinates.size());

    CHECK_NOTHROW(net.realignSlipLinkToImages(u, 2, 1, 0));
    CHECK((net.springBoxOffset.array() == (0.0)).all());

    net.coordinates << 5, 2.5, 0., // crosslink 1,
      5., 7.5, 0.,                 // crosslink 2,
      5., 2.5, 0.;                 // slip-link

    net.springBoxOffset << 10., 0., 0., // slip-link -> crosslink 2
      -10., 0., 0.;                     // slip-link -> crosslink 1

    CHECK_NOTHROW(net.realignSlipLinkToImages(u, 2, 1, 0));

    Eigen::VectorXd expectation = Eigen::VectorXd::Zero(2 * 3);
    expectation(0) = 0.;
    expectation(3) = 0.;

    for (Eigen::Index i = 0; i < 2 * 3; ++i) {
      CHECK(expectation[i] == Catch::Approx(net.springBoxOffset[i]));
    }
  }
}

TEST_CASE("Particular MEHP Force Balance Example", "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"Particular MEHP Force Balance Example\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure.out-"
    "equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         false,
                         false,
                         true);
    CHECK_NOTHROW(forceBalancer.run(pcm::MEHPSlipLinkSolverConfig{}));

    double pressBefore = forceBalancer.getPressure();
    forceBalancer.configAssumeBoxLargeEnough(false);
    CHECK(pressBefore <= forceBalancer.getPressure());
  } else {
    std::cerr << "File " << inputFile << " does not exist." << std::endl;
  }
}

TEST_CASE("Random sampling example", "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"Random sampling example\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPSimulator forceBalancerLargeOldSampling =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         true);
    forceBalancerLargeOldSampling.randomlyAddSlipLinks(
      1000, 6.0, 900, 3.0, "53467829");

    pcm::MEHPSimulator forceBalancerSmallOldSampling1 =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         false);

    pcm::MEHPSimulator forceBalancerSmallOldSampling2 =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         false);

    CHECK_THAT(forceBalancerSmallOldSampling1.getResidualNorm(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling2.getResidualNorm()));

    forceBalancerSmallOldSampling1.configAssumeBoxLargeEnough(false);
    forceBalancerSmallOldSampling1.randomlyAddSlipLinks(
      1000, 6.0, 900, 3.0, "53467829");
    forceBalancerSmallOldSampling2.randomlyAddSlipLinks(
      1000, 6.0, 900, 3.0, "53467829");
    forceBalancerSmallOldSampling2.configAssumeBoxLargeEnough(false);

    CHECK(forceBalancerLargeOldSampling.getPressure() <=
          forceBalancerSmallOldSampling1.getPressure());
    CHECK_THAT(forceBalancerSmallOldSampling1.getResidualNorm(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling2.getResidualNorm()));
    CHECK_THAT(
      forceBalancerSmallOldSampling1.getPressure(),
      Catch::Matchers::WithinRel(forceBalancerSmallOldSampling2.getPressure()));

    CHECK(forceBalancerLargeOldSampling.getNetwork().springIndexA.isApprox(
      forceBalancerSmallOldSampling1.getNetwork().springIndexA));
    CHECK(forceBalancerSmallOldSampling1.getNetwork().springIndexB.isApprox(
      forceBalancerSmallOldSampling2.getNetwork().springIndexB));

    pcm::MEHPSimulator forceBalancerSmallNewSampling =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        1000,
        6.0,
        0.0,
        900,
        3.0,
        "53467829",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        true,
        false,
        false);

    CHECK_THAT(forceBalancerSmallNewSampling.getNrOfSprings(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling1.getNrOfSprings(), 0.1));
    CHECK_THAT(forceBalancerSmallNewSampling.getPressure(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling1.getPressure(), 0.5));
    CHECK_THAT(forceBalancerSmallNewSampling.getResidualNorm(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling1.getResidualNorm(), 0.5));

    REQUIRE(forceBalancerSmallOldSampling1.getNetwork().nrOfStrands ==
            forceBalancerSmallNewSampling.getNetwork().nrOfStrands);
  } else {
    std::cerr << "File " << inputFile << " does not exist." << std::endl;
  }
}

TEST_CASE("Force Balance random sampling example small",
          "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"Force Balance random sampling example small\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath + "square_lattice_2x2_a_5.2d.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       true,
                       false,
                       false,
                       true);
  forceBalancer.randomlyAddSlipLinks(12, 6.0, 11, 3.0, "86573452");

  pcm::MEHPSimulator forceBalancer2 =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       true,
                       false,
                       false,
                       false);

  pcm::MEHPSimulator forceBalancer3 =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       true,
                       false,
                       false,
                       false);

  CHECK_THAT(forceBalancer2.getResidualNorm(),
             Catch::Matchers::WithinRel(forceBalancer3.getResidualNorm()));

  forceBalancer2.configAssumeBoxLargeEnough(false);
  forceBalancer2.randomlyAddSlipLinks(12, 6.0, 11, 3.0, "86573452");
  forceBalancer3.randomlyAddSlipLinks(12, 6.0, 11, 3.0, "86573452");
  forceBalancer3.configAssumeBoxLargeEnough(false);

  CHECK(forceBalancer.getPressure() <= forceBalancer2.getPressure());
  CHECK_THAT(forceBalancer2.getResidualNorm(),
             Catch::Matchers::WithinRel(forceBalancer3.getResidualNorm()));
  CHECK_THAT(forceBalancer2.getPressure(),
             Catch::Matchers::WithinRel(forceBalancer3.getPressure()));

  CHECK(forceBalancer.getNetwork().springIndexA.isApprox(
    forceBalancer2.getNetwork().springIndexA));
  CHECK(forceBalancer2.getNetwork().springIndexB.isApprox(
    forceBalancer3.getNetwork().springIndexB));
}

TEST_CASE("Yet another sampling example", "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"Yet another sampling example\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    // randomly sample slip-links
    pcm::MEHPSimulator forceBalancerNewSampling =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        1000,
        6.0,
        0.0,
        900,
        3.0,
        "53467829",
        2,
        pcm::EntanglementMode::MERGED_LINKS);
    forceBalancerNewSampling.configAssumeBoxLargeEnough(false);

    // extract these randomly sampled slip-links
    std::vector<size_t> strandIdx1;
    std::vector<size_t> strandIdx2;
    std::vector<double> alpha1;
    std::vector<double> alpha2;

    pcm::MEHPNetwork net = forceBalancerNewSampling.getNetwork();
    Eigen::VectorXd x(static_cast<Eigen::Index>(net.nrOfLinks - net.nrOfNodes));
    Eigen::VectorXd y(static_cast<Eigen::Index>(net.nrOfLinks - net.nrOfNodes));
    Eigen::VectorXd z(static_cast<Eigen::Index>(net.nrOfLinks - net.nrOfNodes));
    auto sumAlphaAlongStrandUpToLink = [&net](const size_t strandIdx,
                                              const size_t targetLinkIdx) {
      double alpha = 0.0;
      const std::vector<size_t>& strandSprings =
        net.springIndicesOfStrand[strandIdx];
      for (const size_t sIdx : strandSprings) {
        const Eigen::Index sIdxE = static_cast<Eigen::Index>(sIdx);
        alpha += net.springPartitions[sIdxE];
        if (static_cast<size_t>(net.springIndexB[sIdxE]) == targetLinkIdx) {
          break;
        }
      }
      return alpha;
    };
    Eigen::Index entIdx = 0;
    for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
      if (!net.linkIsEntanglement[static_cast<Eigen::Index>(linkIdx)]) {
        continue;
      }
      strandIdx1.push_back(net.strandIndicesOfLink[linkIdx][0]);
      strandIdx2.push_back(
        pylimer_tools::utils::last(net.strandIndicesOfLink[linkIdx]));
      x[entIdx] = net.coordinates[3 * static_cast<Eigen::Index>(linkIdx) + 0];
      y[entIdx] = net.coordinates[3 * static_cast<Eigen::Index>(linkIdx) + 1];
      z[entIdx] = net.coordinates[3 * static_cast<Eigen::Index>(linkIdx) + 2];
      alpha1.push_back(sumAlphaAlongStrandUpToLink(
        net.strandIndicesOfLink[linkIdx][0], linkIdx));
      alpha2.push_back(sumAlphaAlongStrandUpToLink(
        pylimer_tools::utils::last(net.strandIndicesOfLink[linkIdx]), linkIdx));
      ++entIdx;
      // TODO: handle these primary loops, then increase required accuracy
      // if (net.strandIndicesOfLink[linkIdx].size() == 1 ||
      //     net.strandIndicesOfLink[linkIdx][0] ==
      //       net.strandIndicesOfLink[linkIdx][1]) {
      //   alpha2[alpha2.size() - 1] +=
      // }
    }

    // generate the same slip-links twice,
    // once for each assumption
    pcm::MEHPSimulator forceBalancerOldSamplingLarge =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::MERGED_LINKS,
                         nullptr,
                         false,
                         false,
                         false,
                         false);
    forceBalancerOldSamplingLarge.configAssumeBoxLargeEnough(true);

    // check that the general conversion is equal
    REQUIRE(forceBalancerOldSamplingLarge.getNetwork().springIndexA.isApprox(
      net.springIndexA));

    forceBalancerOldSamplingLarge.addSlipLinks(
      strandIdx1, strandIdx2, x, y, z, alpha1, alpha2, false);
    CHECK(forceBalancerOldSamplingLarge.getNrOfLinks() ==
          forceBalancerNewSampling.getNrOfLinks());
    CHECK(forceBalancerOldSamplingLarge.getNrOfSprings() ==
          forceBalancerNewSampling.getNrOfSprings());

    // std::cout << "\n\n\nnetwork 4:" << std::endl;
    // outputNetwork(forceBalancer4.getNetwork(),
    //               forceBalancer4.getCurrentDisplacements(),
    //               forceBalancer4.getSpringPartitions());
    // std::cout << "\n\n\nnetwork 1:" << std::endl;
    // outputNetwork(forceBalancer.getNetwork(),
    //               forceBalancer.getCurrentDisplacements(),
    //               forceBalancer.getSpringPartitions());

    pcm::MEHPSimulator forceBalancerOldSamplingSmall =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::MERGED_LINKS,
                         nullptr,
                         false,
                         false,
                         false,
                         false);
    forceBalancerOldSamplingSmall.configAssumeBoxLargeEnough(false);
    forceBalancerOldSamplingSmall.addSlipLinks(
      strandIdx1, strandIdx2, x, y, z, alpha1, alpha2, false);

    auto getDisplacementResidualNorm = [](const pcm::MEHPSimulator& sim,
                                          const double limit) {
      Eigen::VectorXd uTmp = sim.getCurrentDisplacements();
      return sim.getNetwork().displaceToMeanPosition(
        uTmp, sim.getIs2D(), limit);
    };

    // after adding slip-links
    CHECK(forceBalancerOldSamplingLarge.getNetwork().springIndexA.isApprox(
      forceBalancerOldSamplingSmall.getNetwork().springIndexA));
    CHECK_THAT(
      getDisplacementResidualNorm(forceBalancerNewSampling, 1.0),
      Catch::Matchers::WithinRel(
        getDisplacementResidualNorm(forceBalancerOldSamplingLarge, 1.0), 1e-3));
    CHECK_THAT(
      forceBalancerOldSamplingLarge.getPressure(),
      Catch::Matchers::WithinRel(forceBalancerNewSampling.getPressure(), 1e-3));

    CHECK(forceBalancerOldSamplingLarge.getPressure() <=
          forceBalancerOldSamplingSmall.getPressure());
    // it is actually thinkable that the following fails for certain scenarios.
    // however, in general, it should not
    CHECK(getDisplacementResidualNorm(forceBalancerOldSamplingLarge, 1.) <=
          getDisplacementResidualNorm(forceBalancerOldSamplingSmall, 1.));
    CHECK(getDisplacementResidualNorm(forceBalancerOldSamplingLarge, -1.) <=
          getDisplacementResidualNorm(forceBalancerOldSamplingSmall, -1.));

    // check to make sure the slip-links are actually placed identically
    CHECK(forceBalancerOldSamplingLarge.getNetwork().springIndexA.isApprox(
      forceBalancerOldSamplingSmall.getNetwork().springIndexA));
    // cannot compare this: The indices are different depending on sampling
    // or addition method
    // CHECK(forceBalancer2.getNetwork().springPartIndexB.isApprox(
    //   forceBalancer4.getNetwork().springPartIndexB));
  } else {
    std::cerr << "File " << inputFile << " does not exist." << std::endl;
  }
}

TEST_CASE("Conversion of structure is equal for both methods",
          "[analysis][MEHPSimulator][long]")
{
  std::cout
    << "Running test \"Conversion of structure is equal for both methods\""
    << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPSimulator forceBalancerNewConversion =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        0,
        6.0,
        0.0,
        0,
        3.0,
        "",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        true,
        false,
        false);

    pcm::MEHPSimulator forceBalancerOldConversion =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         false);

    CHECK(forceBalancerNewConversion.getNetwork().springIndexA.isApprox(
      forceBalancerOldConversion.getNetwork().springIndexA));
    CHECK(forceBalancerNewConversion.getNetwork().springIndexB.isApprox(
      forceBalancerOldConversion.getNetwork().springIndexB));
    CHECK(forceBalancerNewConversion.getNetwork().coordinates.isApprox(
      forceBalancerOldConversion.getNetwork().coordinates));
    CHECK(forceBalancerNewConversion.getNetwork().springBoxOffset.isApprox(
      forceBalancerOldConversion.getNetwork().springBoxOffset));

    CHECK_THAT(
      forceBalancerOldConversion.getResidualNorm(),
      Catch::Matchers::WithinRel(forceBalancerNewConversion.getResidualNorm()));
    CHECK_THAT(
      forceBalancerOldConversion.getPressure(),
      Catch::Matchers::WithinRel(forceBalancerNewConversion.getPressure()));
  }
}

TEST_CASE("Force Balance can be run without slipping",
          "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"Force Balance can be run without slipping\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe,
      1000,
      6.0,
      0.0,
      900,
      3.0,
      "53467829",
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      true,
      false,
      nullptr,
      false,
      true,
      false,
      false);

  Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
  pcm::MEHPSlipLinkSolverConfig cfg;
  cfg.maxNrOfSteps = 500;
  cfg.xtol = 1e-9;
  cfg.simplificationMode = pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
  cfg.inactiveRemovalCutoff = 0.01;
  cfg.slipMode = pcm::SlipMode::FIXED;
  forceBalancer.run(cfg);
  CHECK(springPartitions.isApprox(forceBalancer.getSpringPartitions()));
}

TEST_CASE("Adding slip-links does not influence other springs",
          "[analysis][MEHPSimulator][long]")
{
  std::cout
    << "Running test \"Adding slip-links does not influence other springs\""
    << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPSimulator forceBalancerOldSamplingSmall =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       false,
                       true,
                       false,
                       false);
  Eigen::VectorXd strandVectorsWithoutSlipLinks =
    forceBalancerOldSamplingSmall.getNetwork().evaluateStrandVectors(
      forceBalancerOldSamplingSmall.getCurrentDisplacements(), false);
  forceBalancerOldSamplingSmall.randomlyAddSlipLinks(
    1000, 6.0, 900, 3.0, "53467829");

  Eigen::VectorXd strandVectorsWithSlipLinks =
    forceBalancerOldSamplingSmall.getNetwork().evaluateStrandVectors(
      forceBalancerOldSamplingSmall.getCurrentDisplacements(), false);

  CHECK(strandVectorsWithoutSlipLinks.isApprox(strandVectorsWithSlipLinks));

  pcm::MEHPSimulator forceBalancer2 =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe,
      1000,
      6.0,
      0.0,
      900,
      3.0,
      "53467829",
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      true,
      false,
      nullptr,
      false,
      true,
      false,
      false);

  pcm::MEHPSimulator forceBalancer2Without =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe,
      0,
      6.0,
      0.0,
      0,
      3.0,
      "53467829",
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      true,
      false,
      nullptr,
      false,
      true,
      false,
      false);

  CHECK(forceBalancer2.getNrOfSprings() >
        forceBalancer2Without.getNrOfSprings());

  CHECK(
    forceBalancer2.getNetwork()
      .evaluateStrandVectors(forceBalancer2.getCurrentDisplacements(), false)
      .isApprox(forceBalancer2Without.getNetwork().evaluateStrandVectors(
        forceBalancer2Without.getCurrentDisplacements(), false)));
}

TEST_CASE(
  "MEHP Force Balance Entanglement Links are just like Entanglement Springs",
  "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHP Force Balance Entanglement Links are just "
               "like Entanglement Springs\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath + "mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0."
                    "0114_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
        universe, 4368, 4., 0., 190, 0, "af1346lhkdsaöf123", 2, true);

  pcm::MEHPSimulator forceBalancerEntanglementSprings =
    pcm::MEHPSimulator(universe,
                       entanglements,
                       2,
                       pcm::EntanglementMode::SPRING,
                       nullptr,
                       false,
                       true,
                       false,
                       true);

  pcm::MEHPSimulator forceBalancerEntanglementLinks =
    pcm::MEHPSimulator(universe,
                       entanglements,
                       2,
                       pcm::EntanglementMode::MERGED_LINKS,
                       nullptr,
                       false,
                       true,
                       false,
                       true);

  REQUIRE(forceBalancerEntanglementSprings.getNrOfSprings() >
          forceBalancerEntanglementLinks.getNrOfSprings());
  REQUIRE(forceBalancerEntanglementLinks.getNrOfStrands() == 464);
  REQUIRE(forceBalancerEntanglementSprings.getNrOfActiveSprings() >
          forceBalancerEntanglementLinks.getNrOfActiveSprings());

  auto getActiveAtomIds = [](const pcm::MEHPSimulator& sim) {
    const pcm::MEHPNetwork& net = sim.getNetwork();
    std::vector<int> activeAtomIds;
    for (const int idx : sim.getIndicesOfActiveNodes()) {
      activeAtomIds.push_back(net.oldAtomIds[static_cast<Eigen::Index>(idx)]);
    }
    std::sort(activeAtomIds.begin(), activeAtomIds.end());
    return activeAtomIds;
  };

  REQUIRE(getActiveAtomIds(forceBalancerEntanglementSprings) ==
          getActiveAtomIds(forceBalancerEntanglementLinks));

  SECTION("With removal of inactive springs")
  {
    pcm::MEHPSlipLinkSolverConfig cfg;
    cfg.maxNrOfSteps = 5000;
    cfg.xtol = 1e-15;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
    cfg.inactiveRemovalCutoff = 1e-3;
    cfg.slipMode = pcm::SlipMode::FIXED;
    forceBalancerEntanglementSprings.run(cfg);
    forceBalancerEntanglementLinks.run(cfg);

    CHECK(forceBalancerEntanglementSprings.getNrOfSprings() >
          forceBalancerEntanglementLinks.getNrOfSprings());
    CHECK(getActiveAtomIds(forceBalancerEntanglementSprings) ==
          getActiveAtomIds(forceBalancerEntanglementLinks));

    CHECK_THAT(
      forceBalancerEntanglementSprings.getSolubleWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getSolubleWeightFraction(), 0.005));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getActiveWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getActiveWeightFraction(), 0.005));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getDanglingWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getDanglingWeightFraction(), 0.005));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getStressTensor().trace(),
      Catch::Matchers::WithinRel(
        forceBalancerEntanglementLinks.getStressTensor().trace(), 0.1));
    CHECK(forceBalancerEntanglementSprings.getStressTensor().trace() <
          forceBalancerEntanglementLinks.getStressTensor().trace());
    CHECK(forceBalancerEntanglementSprings.getGammaFactors(1.).sum() <
          forceBalancerEntanglementLinks.getGammaFactors(1.).sum());
  }

  SECTION("Without simplification of the structure")
  {
    pcm::MEHPSlipLinkSolverConfig cfg;
    cfg.maxNrOfSteps = 5000;
    cfg.xtol = 1e-15;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.slipMode = pcm::SlipMode::FIXED;
    forceBalancerEntanglementSprings.run(cfg);
    forceBalancerEntanglementLinks.run(cfg);

    CHECK(forceBalancerEntanglementSprings.getNrOfSprings() >
          forceBalancerEntanglementLinks.getNrOfSprings());
    CHECK(getActiveAtomIds(forceBalancerEntanglementSprings) ==
          getActiveAtomIds(forceBalancerEntanglementLinks));

    CHECK_THAT(
      forceBalancerEntanglementSprings.getSolubleWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getSolubleWeightFraction(), 0.0025));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getActiveWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getActiveWeightFraction(), 0.005));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getDanglingWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getDanglingWeightFraction(), 0.005));
    CHECK(forceBalancerEntanglementSprings.getStressTensor().trace() <
          forceBalancerEntanglementLinks.getStressTensor().trace());
    CHECK(forceBalancerEntanglementSprings.getGammaFactors(1.).sum() <
          forceBalancerEntanglementLinks.getGammaFactors(1.).sum());
  }
}

TEST_CASE(
  "MEHP Force Balance adding Entanglements vs. Phantom behaves as expected",
  "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHP Force Balance adding Entanglements vs. "
               "Phantom behaves as expected\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath =
    std::string(PYLIMER_TEST_FIXTURES_DIR) + "/structure/";

  std::string inputFile =
    suspectedPath + "mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0."
                    "0114_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
        universe, 4368, 4., 0., 190, 0, "af1346lhkdsaöf123", 2, true);

  pcm::MEHPSimulator forceBalancerEntanglements =
    pcm::MEHPSimulator(universe,
                       entanglements,
                       2,
                       pcm::EntanglementMode::MERGED_LINKS,
                       nullptr,
                       false,
                       true,
                       false,
                       true);

  pcm::MEHPSimulator forceBalancerPhantom =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       nullptr,
                       false,
                       true,
                       false,
                       true);

  double initialResidualE =
    forceBalancerEntanglements.getResiduals().squaredNorm();
  double initialResidualP = forceBalancerPhantom.getResiduals().squaredNorm();

  CHECK(forceBalancerEntanglements.getNrOfSprings() >
        forceBalancerPhantom.getNrOfSprings());
  CHECK(forceBalancerEntanglements.getNrOfActiveSprings() >
        forceBalancerPhantom.getNrOfActiveSprings());
  CHECK(forceBalancerEntanglements.getSolubleWeightFraction() <
        forceBalancerPhantom.getSolubleWeightFraction());

  forceBalancerEntanglements.run(pcm::MEHPSlipLinkSolverConfig{});
  forceBalancerPhantom.run(pcm::MEHPSlipLinkSolverConfig{});

  CHECK(forceBalancerEntanglements.getNrOfSprings() >
        forceBalancerPhantom.getNrOfSprings());
  CHECK(forceBalancerEntanglements.getNrOfActiveSprings() >
        forceBalancerPhantom.getNrOfActiveSprings());
  CHECK(forceBalancerEntanglements.getSolubleWeightFraction() <
        forceBalancerPhantom.getSolubleWeightFraction());

  pcm::MEHPSlipLinkSolverConfig cfgE;
  cfgE.maxNrOfSteps = 5000;
  cfgE.xtol = 1e-10;
  cfgE.initialResidualToUse = initialResidualE;
  cfgE.simplificationMode = pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
  cfgE.inactiveRemovalCutoff = 1e-3;
  cfgE.slipMode = pcm::SlipMode::FIXED;

  pcm::MEHPSlipLinkSolverConfig cfgP = cfgE;
  cfgP.initialResidualToUse = initialResidualP;

  forceBalancerEntanglements.run(cfgE);
  forceBalancerPhantom.run(cfgP);

  CHECK(forceBalancerEntanglements.getNrOfSprings() >
        forceBalancerPhantom.getNrOfSprings());
  CHECK(forceBalancerEntanglements.getNrOfActiveSprings() >
        forceBalancerPhantom.getNrOfActiveSprings());
  CHECK(forceBalancerEntanglements.getSolubleWeightFraction() <
        forceBalancerPhantom.getSolubleWeightFraction());
};

TEST_CASE("Temporary force balance test case",
          "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"Temporary force balance test case\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/tmp/";

  std::string inputFile =
    suspectedPath + "mc_own-si_pdms_crosslinked_melt_2590_a_221_2410_monoa_271_"
                    "r_1.44_wsol_0.278_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPSimulator forceBalancerEntanglements =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe,
      26832,
      3.25,
      0.,
      26000,
      0,
      "",
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      true,
      false,
      nullptr,
      false,
      true,
      false,
      false);

  std::cout << "Sampled entanglements " << std::endl;
  pcm::MEHPSlipLinkSolverConfig cfg;
  cfg.maxNrOfSteps = 50000;
  cfg.xtol = 1e-8;
  cfg.simplificationMode = pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
  cfg.inactiveRemovalCutoff = 1e-3;
  cfg.swappingFrequency = 1000;
  cfg.slipMode = pcm::SlipMode::FIXED;
  CHECK_NOTHROW(forceBalancerEntanglements.run(cfg));
};

#ifdef CEREALIZABLE
TEST_CASE("MEHP Force Balance can be serialized",
          "[MEHPSimulator][CerealUtils]")
{
  std::cout << "Running test \"MEHP Force Balance can be serialized\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  CHECK(std::filesystem::exists(suspectedPath));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  CHECK(universeSeq.getLength() == 1);
  pe::Universe universe = universeSeq.atIndex(0);

  pcm::MEHPSimulator forceBalancerEntanglements =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe,
      1000,
      3.25,
      0.,
      100,
      0,
      "",
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      true,
      false,
      nullptr,
      false,
      true,
      false,
      false);

  SECTION("To binary")
  {
    std::stringstream os;
    cereal::BinaryOutputArchive oarchive(os);
    oarchive(forceBalancerEntanglements);

    std::string binString = os.str();

    pcm::MEHPSimulator forceBalancerLoaded =
      pcm::MEHPSimulator::constructFromString(binString);

    CHECK(forceBalancerEntanglements.getNrOfSprings() ==
          forceBalancerLoaded.getNrOfSprings());
    CHECK(forceBalancerEntanglements.getNrOfActiveSprings() ==
          forceBalancerLoaded.getNrOfActiveSprings());
    CHECK(forceBalancerEntanglements.getResidualNorm() ==
          forceBalancerLoaded.getResidualNorm());
    CHECK(forceBalancerEntanglements.getSolubleWeightFraction() ==
          forceBalancerLoaded.getSolubleWeightFraction());
    CHECK(forceBalancerEntanglements.getActiveWeightFraction() ==
          forceBalancerLoaded.getActiveWeightFraction());
    CHECK(forceBalancerEntanglements.getVolume() ==
          forceBalancerLoaded.getVolume());
    CHECK(forceBalancerEntanglements.getCrosslinkerVerse() ==
          forceBalancerLoaded.getCrosslinkerVerse());
  }

  SECTION("To file")
  {
    std::string outputFile =
      suspectedPath + "/tmp_force_balance_serialization.bin";
    forceBalancerEntanglements.writeRestartFile(outputFile);
    CHECK(std::filesystem::exists(outputFile));

    pcm::MEHPSimulator forceBalancerLoaded =
      pcm::MEHPSimulator::readRestartFile(outputFile);

    std::filesystem::remove(outputFile);

    CHECK(forceBalancerEntanglements.getNrOfSprings() ==
          forceBalancerLoaded.getNrOfSprings());
    CHECK(forceBalancerEntanglements.getNrOfActiveSprings() ==
          forceBalancerLoaded.getNrOfActiveSprings());
    CHECK(forceBalancerEntanglements.getResidualNorm() ==
          forceBalancerLoaded.getResidualNorm());
    CHECK(forceBalancerEntanglements.getSolubleWeightFraction() ==
          forceBalancerLoaded.getSolubleWeightFraction());
    CHECK(forceBalancerEntanglements.getActiveWeightFraction() ==
          forceBalancerLoaded.getActiveWeightFraction());
    CHECK(forceBalancerEntanglements.getVolume() ==
          forceBalancerLoaded.getVolume());
    CHECK(forceBalancerEntanglements.getCrosslinkerVerse() ==
          forceBalancerLoaded.getCrosslinkerVerse());
  }
}
#endif
