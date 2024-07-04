#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
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
roundForOutput(double val, double precision = 5)
{
  return std::round(val * std::pow(10, precision)) / std::pow(10, precision);
}

void
outputNetwork(pcm::ForceBalanceNetwork net,
              Eigen::VectorXd displacements,
              Eigen::VectorXd springPartitions)
{
  for (int i = 0; i < net.nrOfSprings; ++i) {
    std::cout << "Spring " << i << ", N: " << net.springsContourLength[i]
              << std::endl;
    for (int j = 0; j < net.linkIndicesOfSprings[i].size(); ++j) {
      size_t linkIdx = net.linkIndicesOfSprings[i][j];
      std::cout << linkIdx;
      if (!net.linkIsSliplink[linkIdx]) {
        std::cout << " (" << net.oldAtomIds[linkIdx] << ")";
      }
      std::cout << ": ";
      for (int dir = 0; dir < 3; ++dir) {
        std::cout << (net.coordinates[3 * linkIdx + dir] +
                      displacements[3 * linkIdx + dir])
                  << ", ";
      }
      std::cout << std::endl;
      if (j < net.linkIndicesOfSprings[i].size() - 1) {
        size_t partialSpringIdx = net.localToGlobalSpringIndex[i][j];
        std::cout
          << partialSpringIdx << ": " << springPartitions[partialSpringIdx]
          << " ("
          << roundForOutput(net.springPartBoxOffset[3 * partialSpringIdx + 0])
          << ", "
          << roundForOutput(net.springPartBoxOffset[3 * partialSpringIdx + 1])
          << ", "
          << roundForOutput(net.springPartBoxOffset[3 * partialSpringIdx + 2])
          << ")" << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

TEST_CASE("Particular slip-link examples", "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"Particular slip-link examples\"" << std::endl;
  double L = 42.819955007276754;
  double lmda = 1.2;
  pe::Universe universe = pe::Universe(L, L, L);
  /**
   * Connectivity:
   *
   * 35-(11)-90
   *
   * 10-(12)-1654
   */
  // slip-link 16 in the test-system
  universe.addAtoms({ 35, 90, 1654, 10, 11, 12, 13, 14 },
                    { 2, 2, 2, 2, 1, 1, 1, 1 },
                    { 19.706880857235795,
                      19.288603889563976,
                      22.156152142687819,
                      10,
                      10,
                      10,
                      10,
                      10 },
                    { 1.47224612942217,
                      4.4207926048800461,
                      2.821003235624608,
                      10,
                      10,
                      10,
                      10,
                      10 },
                    { 22.98584649043724,
                      24.328555987207494,
                      24.562956857368366,
                      10,
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

  pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
  forceBalancer.addSlipLinks(
    { 1, 0 },
    { 1, 1 },
    { 22.650980696700437, 15.67228278755659 },
    { 38.666887107085628, 3.172504627794056 },
    { 20.999752090751294, 27.126973321608833 },
    { 1. - 0.448276, 0.31034482758620685 },
    { 1. - 0.448276, 1. - (0.6896551724137931 * 0.448276) });

  // do update step
  Eigen::VectorXd displacements = Eigen::VectorXd::Zero(6 * 3);
  Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
  CHECK(springPartitions[0] == Catch::Approx(0.31034482758620685));
  CHECK(springPartitions[1] == Catch::Approx(1. - 0.448276));
  CHECK(springPartitions[2] == Catch::Approx(0.0));
  CHECK(springPartitions[3] == Catch::Approx(0.13912));
  CHECK(springPartitions[4] == Catch::Approx(1. - 0.31034482758620685));
  // outputNetwork(forceBalancer.getNetwork(), displacements, springPartitions);
  /*auto results = */
  forceBalancer.inspectParametrisationOptimsationForLink(
    5,
    displacements,
    springPartitions,
    250,
    1e-10,
    1,
    1e10); // cannot use 1.0 for oneOver... without setting higher contour
           // length fraction
  // outputNetwork(forceBalancer.getNetwork(), displacements, springPartitions);
  CHECK(displacements[5 * 3] == Catch::Approx(5.89283));
  CHECK(displacements[5 * 3 + 1] == Catch::Approx(-1.0033));
  CHECK(displacements[5 * 3 + 2] == Catch::Approx(-3.49778));
  CHECK(springPartitions[0] == Catch::Approx(0.388968));
  CHECK(springPartitions[5] == Catch::Approx(0.112087));
  CHECK(springPartitions[2] == Catch::Approx(0.0));
  CHECK(springPartitions[1] == Catch::Approx(1. - 0.448276));
  CHECK(springPartitions[4] == Catch::Approx(1. - 0.388968));
  CHECK(springPartitions[3] == Catch::Approx((1. - 0.25004) * 0.448276));
}

TEST_CASE("MC swap accept and reject work", "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MC swap accept and reject work\"" << std::endl;
  double L = 42.819955007276754;
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
  pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
  CHECK_THROWS(
    forceBalancer.setSpringContourLengths(Eigen::VectorXd::Constant(2, 30.)));
  forceBalancer.setSpringContourLengths(Eigen::VectorXd::Constant(3, 30.));

  SECTION("MC condition accepts as requested")
  {
    forceBalancer.addSlipLinks({ 0, 0 },
                               { 1, 1 },
                               { 21.03, 23.01 },
                               { 18.02, 18.03 },
                               { 0.0, 0.0 },
                               { 0.601, 0.4 },
                               { 0.401, 0.6 });
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd u = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    // outputNetwork(net, u, partitions);
    CHECK(forceBalancer.swapSlipLinkReversibly(
      net, u, partitions, 5, 1., -1, false, true));
    bool netIsIdentical = net.springPartIndexA.isApprox(
                            forceBalancer.getNetwork().springPartIndexA) &&
                          net.springPartIndexB.isApprox(
                            forceBalancer.getNetwork().springPartIndexB);
    CHECK_FALSE(netIsIdentical);
    // outputNetwork(net, u, partitions);
    CHECK_FALSE(partitions.isApprox(forceBalancer.getSpringPartitions()));
    CHECK_FALSE(u.isApprox(forceBalancer.getCurrentDisplacements()));
  }

  SECTION("MC condition rejects as requested")
  {
    forceBalancer.addSlipLinks({ 0, 0 },
                               { 1, 1 },
                               { 21., 23. },
                               { 18., 18. },
                               { 0.0, 0.0 },
                               { 0.4, 0.6 },
                               { 0.4, 0.6 });
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd u = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    // outputNetwork(net, u, partitions);
    CHECK_FALSE(
      forceBalancer.swapSlipLinkReversibly(net, u, partitions, 5, 1.));
    CHECK(partitions.isApprox(forceBalancer.getSpringPartitions()));
    CHECK(u.isApprox(forceBalancer.getCurrentDisplacements()));
  }
}

TEST_CASE("MC swap accept and reject work with crosslinkers",
          "[analysis][MEHPForceBalance]")
{
  std::cout
    << "Running test \"MC swap accept and reject work with crosslinkers\""
    << std::endl;
  double L = 42.819955007276754;
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

  pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
  forceBalancer.setSpringContourLengths(Eigen::VectorXd::Constant(4, 30.));
  pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
  CHECK(net.nrOfPartialSprings == net.nrOfSprings);
  CHECK(net.nrOfSprings == 4);

  SECTION("MC condition accepts as requested")
  {
    forceBalancer.addSlipLinks(
      { 0 }, { 1 }, { 10. }, { 4.9 }, { 0. }, { 0.5 }, { 0.966667 });
    net = forceBalancer.getNetwork();
    CHECK(net.nrOfPartialSprings != net.nrOfSprings);
    CHECK(net.nrOfSprings == 4);
    CHECK(net.nrOfPartialSprings == 6);
    CHECK(net.nrOfCrosslinkSwapsEndured[0] == 0);
    Eigen::VectorXd u = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    CHECK(net.linkIndicesOfSprings[1].size() == 3);
    // outputNetwork(net, u, partitions);
    CHECK_FALSE(forceBalancer.swapSlipLinkReversibly(
      net, u, partitions, 5, 1., 0, true, false));
    CHECK(forceBalancer.swapSlipLinkReversibly(
      net, u, partitions, 5, 1., 5, false, true));
    // check that the connectivity around the cross-link changed
    // outputNetwork(net, u, partitions);
    CHECK(net.linkIndicesOfSprings[1].size() == 2);
    CHECK(net.nrOfCrosslinkSwapsEndured[0] == 1);
  }

  SECTION("MC condition rejects as requested")
  {
    forceBalancer.addSlipLinks(
      { 0 }, { 2 }, { 10. }, { 4.9 }, { 0. }, { 0.5 }, { 0.966667 });
    net = forceBalancer.getNetwork();
    CHECK(net.nrOfPartialSprings != net.nrOfSprings);
    CHECK(net.nrOfSprings == 4);
    CHECK(net.nrOfPartialSprings == 6);
    Eigen::VectorXd u = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    // outputNetwork(net, u, partitions);
    CHECK_FALSE(
      forceBalancer.swapSlipLinkReversibly(net, u, partitions, 4, 1.));
    CHECK(net.nrOfCrosslinkSwapsEndured[0] == 0);
    // ideally, we could check that the values stayed the same, but they
    // probably did not as of the current implementation. at least we can check
    // that the configuration stayed the same
  }
}

TEST_CASE("MEHP Force Balance handles slip-links on primary loops",
          "[analysis][MEHPForceBalance]")
{
  std::cout
    << "Running test \"MEHP Force Balance handles slip-links on primary loops\""
    << std::endl;
  pe::Universe universe =
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
  universe.addAtoms({ 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 },
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
  universe.addBonds({ 3, 3, 3, 3, 4, 5, 6, 7, 9, 10, 11, 12 },
                    { 9, 11, 12, 13, 14, 15, 16, 10, 13, 16, 15, 14 });

  SECTION("Unentangled primary loop")
  {
    pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    // check unentangled primary loops
    Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    CHECK_NOTHROW(
      forceBalancer.displaceToMeanPosition(net, displacements, partitions, 0));
    CHECK(displacements[0] == Catch::Approx(7.687).epsilon(1e-5));
    CHECK(displacements[1] == Catch::Approx(2.03926).epsilon(1e-5));
    CHECK(displacements[2] == Catch::Approx(5.88634).epsilon(1e-5));
  }

  SECTION("Entangled primary loop")
  {
    pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
    forceBalancer.setSpringContourLengths(
      Eigen::VectorXd::Constant(forceBalancer.getNetwork().nrOfSprings, 15.));
    // entangle primary loop and check again
    // outputNetwork(net, Eigen::VectorXd::Zero(net.nrOfLinks * 3));
    forceBalancer.addSlipLinks({ 0, 2 },
                               { 1, 1 },
                               { 1.3263401618628183, 42.04664022316877 },
                               { 6.670300217824844, 7.18272624553976 },
                               { 41.85951429390015, 0.8578704100544575 },
                               { 0.5333333333333333, 0.9310344827586207 },
                               { 0.5, 0.5 });
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    // outputNetwork(net,
    //               Eigen::VectorXd::Zero(net.nrOfLinks * 3),
    //               forceBalancer.getSpringPartitions());
    Eigen::VectorXd displacements = Eigen::VectorXd::Zero(net.nrOfLinks * 3);
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    CHECK_NOTHROW(forceBalancer.displaceToMeanPosition(
      net, displacements, partitions, 0, 1.));
    CHECK(displacements[0] == Catch::Approx(0.187321).epsilon(1e-5));
    CHECK(displacements[1] == Catch::Approx(-0.447774).epsilon(1e-5));
    CHECK(displacements[2] == Catch::Approx(-0.925295).epsilon(1e-5));
  }
}

TEST_CASE("MEHP Force Balance handles slip-link convergence correctly",
          "[analysis][MEHPForceBalance]")
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
  universe.addAtoms({ 35, 90, 1654, 10, 11, 12, 13, 14 },
                    { 2, 2, 2, 2, 1, 1, 1, 1 },
                    { 12.075848854154861,
                      10.644563425246883,
                      14.302483570484272,
                      10,
                      10,
                      10,
                      10,
                      10 },
                    { 3.574724359149917,
                      5.460837527830988,
                      3.718195811318871,
                      10,
                      10,
                      10,
                      10,
                      10 },
                    { 3.1018436428667284,
                      7.956714096296886,
                      4.4007749635446824,
                      10,
                      10,
                      10,
                      10,
                      10 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 },
                    { 1, 1, 1, 1, 1, 1, 1, 1 });
  universe.addBonds({ 35, 11, 10, 12 }, { 11, 90, 12, 1654 });

  pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
  forceBalancer.configAssumeBoxLargeEnough(true);
  forceBalancer.addSlipLinks({ 1, 0 },
                             { 1, 1 },
                             { 12.650493316819828, 13.197029579176265 },
                             { 2.8706102036538566, 3.4016980009809297 },
                             { 8.475863644409664, 5.284899588057222 },
                             { 0.13793103448275862, 1. - 0.3103448275862069 },
                             { 0.13793103448275862, 0.7931034482758621 },
                             false);

  Eigen::VectorXd displacements = Eigen::VectorXd::Zero(6 * 3);
  Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
  CHECK(springPartitions[3] ==
        Catch::Approx(0.7931034482758621 - 0.13793103448275862));
  // outputNetwork(forceBalancer.getNetwork(), displacements, springPartitions);
  // one output step
  forceBalancer.updateSpringPartition(
    forceBalancer.getNetwork(), displacements, springPartitions, 5, 1.e10);
  CHECK(springPartitions[0] == Catch::Approx(0.367729));
  CHECK(springPartitions[4] == Catch::Approx(1. - 0.367729));
  CHECK(springPartitions[1] == Catch::Approx(0.13793103448275862));
  // TODO: the Mathematica example `exemplaric-sliplink-3` proposes a different
  // solution here – check! CHECK(springPartitions[3] == Catch::Approx(0.655172
  // - 0.306578));
  forceBalancer.displaceToMeanPosition(
    forceBalancer.getNetwork(), displacements, springPartitions, 5, 1.e10);
  CHECK(displacements[5 * 3] == Catch::Approx(-0.391233));
  CHECK(displacements[5 * 3 + 1] == Catch::Approx(0.413665));
  CHECK(displacements[5 * 3 + 2] == Catch::Approx(0.028894));

  displacements.setZero();
  springPartitions = forceBalancer.getSpringPartitions();
  // do many update steps
  /*auto results = */
  forceBalancer.inspectParametrisationOptimsationForLink(
    5,
    displacements,
    springPartitions,
    250,
    1.e-10,
    100,
    1.e10); // cannot use 1.0 for oneOver... without setting higher contour
            // length fraction
  CHECK(displacements[5 * 3] == Catch::Approx(-0.592091));
  CHECK(displacements[5 * 3 + 1] == Catch::Approx(0.441203));
  CHECK(displacements[5 * 3 + 2] == Catch::Approx(0.44577));
  CHECK(springPartitions[1] == Catch::Approx(0.13793103448275862));

  // outputNetwork(forceBalancer.getNetwork(), displacements, springPartitions);

  SECTION("Stress tensor computations are equivalent")
  {
    Eigen::Matrix3d stressTensor1 = forceBalancer.getStressTensor();
    Eigen::Matrix3d stressTensor2 = forceBalancer.getStressTensorLinkBased();
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        CHECK(stressTensor1(i, j) == Catch::Approx(stressTensor2(i, j)));
      }
    }
  }
}

TEST_CASE("MEHP Force Balance runs", "[analysis][MEHPForceBalance][long]")
{
  std::cout << "Running test \"MEHP Force Balance runs\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("MEHP Force Balance 3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      std::cout << "Reading file " << largeInputFile << std::endl;
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      std::cout << "Read file " << largeInputFile << std::endl;

      // BENCHMARK_ADVANCED("MEHP LD_MMA " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceBalance forceBalancer3 =
      //     pcm::MEHPForceBalance(universe2, 2);
      //   meter.measure([&forceBalancer3] {
      //     forceBalancer3.runForceRelaxation("LD_MMA");
      //     return forceBalancer3.getNrOfIterations();
      //   });
      // };
      // BENCHMARK_ADVANCED("MEHP LD_LBFGS " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceBalance forceBalancer3 =
      //     pcm::MEHPForceBalance(universe2, 2);
      //   meter.measure([&forceBalancer3] {
      //     forceBalancer3.runForceRelaxation("LD_LBFGS");
      //     return forceBalancer3.getNrOfIterations();
      //   });
      // };

      double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      pcm::MEHPForceBalance forceBalancer2 =
        pcm::MEHPForceBalance(universe2, 2);

      // SECTION("Stress tensor computations are equivalent")
      // {
      Eigen::Matrix3d stressTensor1 = forceBalancer2.getStressTensor();
      Eigen::Matrix3d stressTensor2 = forceBalancer2.getStressTensorLinkBased();
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          CHECK(stressTensor1(i, j) == Catch::Approx(stressTensor2(i, j)));
        }
      }

      SECTION("Actual balance results in correct phantom results")
      {
        std::cout << "Doing phantom force balance" << std::endl;
        pcm::MEHPForceRelaxation forceRelaxer =
          pcm::MEHPForceRelaxation(universe2, 2);
        // the strands are different -> cannot compare the distances anymore
        // CHECK((forceBalancer2.getCurrentSpringDistances() -
        // forceRelaxer.getCurrentSpringDistances()).isMuchSmallerThan(1e-12));
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
        CHECK_NOTHROW(forceBalancer2.runForceRelaxation());
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
        double nrOfChainCorrection =
          (forceBalancer2.getDefaultNrOfChains() / nrOfChains);
        double expectedNb2 = slope * Nb * beadMass;
        double nb2Correction = 1.;
        // (forceBalancer2.getDefaultR0Square() / (expectedNb2));
        double gammaCorrectionFactor = nrOfChainCorrection * nb2Correction;
        CHECK(
          forceBalancer2.getGammaFactor() //* nrOfChainCorrection * 1.
                                          // forceBalancer2.getDefaultR0Square()
          ==
          Catch::Approx(42.6132)); // as from conversion-less Mathematica script
        CHECK(forceBalancer2.getGammaFactor() * kb *
                T * // gammaCorrectionFactor *
                nu ==
              Catch::Approx(61308.3)); // ANT shear modulus, Pa
        CHECK(forceBalancer2.getGammaFactor() * gammaCorrectionFactor ==
              Catch::Approx(0.319446)); // "correct" gamma factor
        CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
        // TODO: find better, more accurate tests here
        CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
        CHECK(forceBalancer2.getNrOfActiveSprings() > 1);
        CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
        CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
              forceBalancer2.getNrOfNodes());
      }
      // also
      SECTION("Actual balance results in correct slip-link results")
      {
        std::cout << "Doing non-phantom force balance" << std::endl;
        pcm::MEHPForceRelaxation forceRelaxer =
          pcm::MEHPForceRelaxation(universe2, 2);
        // the strands are different -> cannot compare the distances anymore
        // CHECK((forceBalancer2.getCurrentSpringDistances() -
        // forceRelaxer.getCurrentSpringDistances()).isMuchSmallerThan(1e-12));
        CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
        CHECK(forceBalancer2.getNrOfIterations() == 0);
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(97.383096 * 97.383096 * 97.383096));
        // initial system values
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(forceRelaxer.getPressure()));
        CHECK_THAT(forceBalancer2.getPressure(),
                   Catch::Matchers::WithinRel(0.0061105865, 1e-3));
        // add entanglements
        // TODO: these are random values, as are the results... :P
        size_t nrOfSprings = forceRelaxer.getNetwork().nrOfSprings;
        forceBalancer2.addSlipLinks(
          { 10, 100, 50, 12, 76, 80, nrOfSprings - 1 },
          { 99, 101, 13, 7, 5, 19, nrOfSprings - 7 },
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
        CHECK_NOTHROW(forceBalancer2.runForceRelaxation());
        // TODO: replace this value thereafter
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
      { { suspectedPath +
          "structure/equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
          "333_2d_t_7500001.structure.out" } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, true);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::UNSET);
    CHECK(forceBalancer.getNrOfIterations() == 0);
    CHECK(forceBalancer.getVolume() == Catch::Approx(universe.getVolume()));
    forceBalancer.runForceRelaxation(5);
    CHECK(forceBalancer.getNrOfNodes() != universe.getNrOfAtoms());
    CHECK(forceBalancer.getNrOfIterations() <= 5);
    CHECK(forceBalancer.getNrOfIterations() >= 1);
    CHECK(universe.getAtomsOfType(2).size() == 7200);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::MAX_STEPS);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    // run again, this time fully
    pcm::MEHPForceBalance forceBalancer2 =
      pcm::MEHPForceBalance(universe, 2, true);
    forceBalancer2.runForceRelaxation(10000, 1e-10, 100);
    CHECK(forceBalancer2.getNrOfIterations() > 5);
    CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer2.getGammaFactor(1., forceBalancer2.getNrOfSprings()) ==
          Catch::Approx(1. / 3.).epsilon(0.001));
    auto stressTensor = forceBalancer2.getStressTensor();
    CHECK(forceBalancer2.getPressure() ==
          Catch::Approx(
            (stressTensor(0, 0) + stressTensor(1, 1) + stressTensor(2, 2)) / 3.)
            .epsilon(0.02));
    CHECK_NOTHROW(forceBalancer2.validateNetwork());
    // TODO: find better, more accurate tests here
    CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
    CHECK(forceBalancer2.getNrOfActiveSprings() > 1);
    CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
    CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
          forceBalancer2.getNrOfNodes());

    pe::Universe universe3 = forceBalancer2.getCrosslinkerVerse();
    CHECK(universe3.getNrOfAtoms() == forceBalancer2.getNrOfNodes());
    CHECK(universe3.getNrOfBonds() == forceBalancer2.getNrOfSprings());
    CHECK(universe3.getAtomsOfType(2).size() == universe3.getNrOfAtoms());

    // try out different algorithms
    std::vector<std::string> algorithms = { "LD_MMA",
                                            // "LD_TNEWTON_PRECOND_RESTART",
                                            // "GD_STOGO",
                                            "LD_SLSQP",
                                            "GN_DIRECT" };

    // for (std::string algorithm : algorithms) {
    //   pcm::MEHPForceBalance forceBalancerN =
    //     pcm::MEHPForceBalance(universe, 2);
    //   std::cout << "Testing algorithm " << algorithm << std::endl;
    //   auto start = std::chrono::high_resolution_clock::now();
    //   forceBalancerN.runForceRelaxation(true, 15, algorithm.c_str(), 10000);
    //   auto stop = std::chrono::high_resolution_clock::now();
    //   CHECK(forceBalancerN.getGammaEq() ==
    //         Catch::Approx(forceBalancer2.getGammaEq()));
    //   CHECK(forceBalancerN.getFinalPressure() ==
    //         Catch::Approx(forceBalancer2.getFinalPressure()));
    //   auto duration = duration_cast<std::chrono::microseconds>(stop - start);
    //   std::cout << "Took: " << duration.count() << std::endl;
    // }
  }
}

TEST_CASE(
  "MEHP Force Balance can randomly add slip-links ignoring crosslinkers",
  "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHP Force Balance can randomly add slip-links "
               "ignoring crosslinkers\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  std::string inputFile =
    suspectedPath + "structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, false, false, false);
    // TODO: using a seed does not seem to work properly
    size_t nrOfAddedLinks =
      forceBalancer.randomlyAddSliplinks(250, 2.0, 100, 2.0, true, 12);
    CHECK(nrOfAddedLinks >= 50);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
    nrOfAddedLinks =
      forceBalancer.randomlyAddSliplinks(250, 2.0, 100, 2.0, false, 25);
    CHECK(nrOfAddedLinks >= 50);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    size_t numRemoved = forceBalancer.removeTwofunctionalCrosslinks(
      net, displacements, partitions);
    CHECK_NOTHROW(
      forceBalancer.validateNetwork(net, displacements, partitions));
    CHECK(numRemoved > 0);
    // numRemoved = forceBalancer.removeInactiveCrosslinks(net, displacements,
    // partitions, 1e-20); CHECK(numRemoved == 204); // TODO: analyze these
    size_t nrOfSpringsBefore = net.nrOfSprings;
    CHECK_NOTHROW(
      forceBalancer.validateNetwork(net, displacements, partitions));
    // remove all springs...
    numRemoved = forceBalancer.removeInactiveCrosslinks(
      net, displacements, partitions, 1e5);
    CHECK(numRemoved == nrOfSpringsBefore);
  }
}

TEST_CASE("MEHP Force Balance can run with swapping slip-links",
          "[analysis][MEHPForceBalance][long]")
{
  std::cout
    << "Running test \"MEHP Force Balance can run with swapping slip-links\""
    << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  std::string inputFile =
    suspectedPath + "structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, false, true, false);
    size_t nrOfAddedLinks =
      forceBalancer.randomlyAddSliplinks(250, 2.0, 100, 2.0, true, 1209);
    CHECK(nrOfAddedLinks >= 25);
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
    nrOfAddedLinks =
      forceBalancer.randomlyAddSliplinks(250, 2.0, 100, 2.0, false, 1230);
    CHECK(nrOfAddedLinks >= 25);
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
    pe::Box oldBox = universe.getBox();
    pe::Box newBox =
      pe::Box(4 * oldBox.getLx(), 0.5 * oldBox.getLy(), 0.5 * oldBox.getLz());
    forceBalancer.deformTo(newBox);

    SECTION("Assuming too small box")
    {
      forceBalancer.configAssumeBoxLargeEnough(false);
      CHECK_NOTHROW(forceBalancer.runForceRelaxation(
        100,
        1e-7,
        -1.0,
        pcm::StructureSimplificationMode::ALL_TIM,
        0.01,
        false,
        pcm::LinkSwappingMode::ALL_MC));
      CHECK_NOTHROW(forceBalancer.runForceRelaxation(
        100,
        1e-8,
        -1.0,
        pcm::StructureSimplificationMode::ALL_TIM,
        0.01,
        false,
        pcm::LinkSwappingMode::SLIPLINKS_ONLY));
      CHECK_NOTHROW(forceBalancer.runForceRelaxation(
        100,
        1e-9,
        -1.0,
        pcm::StructureSimplificationMode::ALL_TIM,
        0.01,
        false,
        pcm::LinkSwappingMode::ALL));
    }

    SECTION("Assuming large enough box")
    {
      forceBalancer.configAssumeBoxLargeEnough(true);
      CHECK_NOTHROW(forceBalancer.runForceRelaxation(
        100,
        1e-7,
        -1.0,
        pcm::StructureSimplificationMode::ALL_TIM,
        0.01,
        false,
        pcm::LinkSwappingMode::ALL_MC));
      CHECK_NOTHROW(forceBalancer.runForceRelaxation(
        100,
        1e-8,
        -1.0,
        pcm::StructureSimplificationMode::ALL_TIM,
        0.01,
        false,
        pcm::LinkSwappingMode::SLIPLINKS_ONLY));
      CHECK_NOTHROW(forceBalancer.runForceRelaxation(
        100,
        1e-9,
        -1.0,
        pcm::StructureSimplificationMode::ALL_TIM,
        0.01,
        false,
        pcm::LinkSwappingMode::ALL));
    }
  }
}

TEST_CASE(
  "MEHP Force Balance can randomly add and remove slip-links with large box",
  "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHP Force Balance can randomly add and remove "
               "slip-links with large box\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  std::string inputFile =
    suspectedPath + "structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, false, false, false);
    forceBalancer.configAssumeBoxLargeEnough(true);
    size_t nrOfAddedLinks = forceBalancer.randomlyAddSliplinks(1000, 2.0, 100);
    CHECK(nrOfAddedLinks >= 100);
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    size_t numRemoved = forceBalancer.removeTwofunctionalCrosslinks(
      net, displacements, partitions);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numRemoved > 0);

    // run a while to get inactive links
    forceBalancer.runForceRelaxation(100);
    net = forceBalancer.getNetwork();
    displacements = forceBalancer.getCurrentDisplacements();
    partitions = forceBalancer.getSpringPartitions();
    // due to the randomness, it _could_ be one day that actually all strands
    // are active. unlikely, but I can imagine it to be possible.
    size_t numInactiveRemoved = forceBalancer.removeInactiveCrosslinks(
      net, displacements, partitions, 0.1);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numInactiveRemoved > 0);
    CHECK_NOTHROW(
      forceBalancer.validateNetwork(net, displacements, partitions));

    ////////////////////////////////////////////////////////////////
    forceBalancer = pcm::MEHPForceBalance(universe, 2, true, 1.0, true);
    nrOfAddedLinks = forceBalancer.randomlyAddSliplinks(1000, 2.0, 100);
    CHECK(nrOfAddedLinks >= 100);
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
    // check that all f = 2 have already been removed
    // they have not, since more f = 2 are produced by
    // numInactiveRemoved = forceBalancer.removeTwofunctionalCrosslinks(
    //       net, displacements, partitions);
    // CHECK(numInactiveRemoved == 0);

    // run a while to get inactive links
    forceBalancer.runForceRelaxation(100);
    net = forceBalancer.getNetwork();
    displacements = forceBalancer.getCurrentDisplacements();
    partitions = forceBalancer.getSpringPartitions();
    // due to the randomness, it _could_ be one day that actually all strands
    // are active. unlikely, but I can imagine it to be possible.
    numInactiveRemoved = forceBalancer.removeInactiveCrosslinks(
      net, displacements, partitions, 0.1);
    CHECK(numInactiveRemoved > 0);
    numInactiveRemoved = forceBalancer.removeTwofunctionalCrosslinks(
      net, displacements, partitions);
    CHECK(numInactiveRemoved > 0);
    CHECK_NOTHROW(
      forceBalancer.validateNetwork(net, displacements, partitions));
  }
}

TEST_CASE(
  "MEHP Force Balance can randomly add and remove slip-links with small box",
  "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHP Force Balance can randomly add and remove "
               "slip-links with small box\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  std::string inputFile =
    suspectedPath + "structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, false, false, false);
    forceBalancer.configAssumeBoxLargeEnough(false);
    size_t nrOfAddedLinks = forceBalancer.randomlyAddSliplinks(1000, 2.0, 100);
    CHECK(nrOfAddedLinks >= 100);
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    size_t numRemoved = forceBalancer.removeTwofunctionalCrosslinks(
      net, displacements, partitions);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numRemoved > 0);

    // run a while to get inactive links
    forceBalancer.runForceRelaxation(100);
    net = forceBalancer.getNetwork();
    displacements = forceBalancer.getCurrentDisplacements();
    partitions = forceBalancer.getSpringPartitions();
    // due to the randomness, it _could_ be one day that actually all strands
    // are active. unlikely, but I can imagine it to be possible.
    size_t numInactiveRemoved = forceBalancer.removeInactiveCrosslinks(
      net, displacements, partitions, 0.1);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numInactiveRemoved > 0);
    CHECK_NOTHROW(
      forceBalancer.validateNetwork(net, displacements, partitions));

    ////////////////////////////////////////////////////////////////
    forceBalancer = pcm::MEHPForceBalance(universe, 2, true, 1.0, true);
    nrOfAddedLinks = forceBalancer.randomlyAddSliplinks(1000, 2.0, 100);
    CHECK(nrOfAddedLinks >= 100);
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
    // check that all f = 2 have already been removed
    // they have not, since more f = 2 are produced by
    // numInactiveRemoved = forceBalancer.removeTwofunctionalCrosslinks(
    //       net, displacements, partitions);
    // CHECK(numInactiveRemoved == 0);

    // run a while to get inactive links
    forceBalancer.runForceRelaxation(100);
    net = forceBalancer.getNetwork();
    displacements = forceBalancer.getCurrentDisplacements();
    partitions = forceBalancer.getSpringPartitions();
    // due to the randomness, it _could_ be one day that actually all strands
    // are active. unlikely, but I can imagine it to be possible.
    numInactiveRemoved = forceBalancer.removeInactiveCrosslinks(
      net, displacements, partitions, 0.1);
    CHECK(numInactiveRemoved > 0);
    numInactiveRemoved = forceBalancer.removeTwofunctionalCrosslinks(
      net, displacements, partitions);
    CHECK(numInactiveRemoved > 0);
    CHECK_NOTHROW(
      forceBalancer.validateNetwork(net, displacements, partitions));
  }
}

TEST_CASE("MEHP Force Balance handles slip-links",
          "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHP Force Balance handles slip-links\""
            << std::endl;
  // construct an example network
  pe::Universe universe = pe::Universe(1.0, 1.0, 1.0);
  /**
   * The system looks like this (in terms of bonds, not 3D placement):
   *
   * 1-5-2
   * |  /|
   * 8 9 6
   * |/  |
   * 4-7-3
   */
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

    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, false);

    forceBalancer.setSpringContourLengths(
      Eigen::VectorXd::Constant(forceBalancer.getNetwork().nrOfSprings, 15.));

    forceBalancer.addSlipLinks(
      { 0 }, { 2 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.62 }, { 0.43 });

    Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
    CHECK(springPartitions[0] == Catch::Approx(0.62));
    CHECK(springPartitions[2] == Catch::Approx(0.43));
    CHECK(springPartitions[6] == Catch::Approx(1. - 0.43));
  }

  SECTION("Relevant slip-links behave correctly")
  {
    /**
     * This adjusted system looks like this (in terms of bonds, not 3D
     * placement):
     *
     * 1-5-2
     * |\ /|
     * 8 * 6
     * |/ \|
     * 4-7-3
     *
     * where * = 9 & 10, where the latter is connected to 1 & 3,
     * the former to 4 and 2
     */
    universe.addAtoms(
      { 10 }, { 1 }, { 0. }, { 0. }, { 0. }, { 1 }, { 1 }, { 1 });
    universe.addBonds(2, { { 1, 10 } }, { { 10, 3 } });
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, false);
    // add a slip-link between the strands 2-4 & 1-3
    CHECK_NOTHROW(
      forceBalancer.addSlipLinks({ 4 }, { 5 }, { 4.2 }, { 3.9 }, { 1.2 }));
    // ...at the wrong coordinates, to see it converge to the center
    Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
    Eigen::VectorXd displacements =
      Eigen::VectorXd::Zero(forceBalancer.getNrOfLinks() * 3);
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();

    // outputNetwork(net, displacements, forceBalancer.getSpringPartitions());
    for (int i = 0; i < 5; ++i) {
      forceBalancer.displaceToMeanPosition(
        net, displacements, springPartitions, 4);
      forceBalancer.updateSpringPartition(
        net, displacements, springPartitions, 4);
    }
    CHECK(springPartitions[springPartitions.size() - 1] == Catch::Approx(0.5));
    CHECK(displacements[3 * 4] == Catch::Approx(-4.2));
    CHECK(displacements[3 * 4 + 1] == Catch::Approx(-3.9));
    CHECK(displacements[3 * 4 + 2] == Catch::Approx(0.8));
  }

  SECTION("Correct rationalisation and displacement")
  {
    // now, construct the force balancer
    pcm::MEHPForceBalance forceBalancer2 =
      pcm::MEHPForceBalance(universe, 2, false);
    CHECK(forceBalancer2.getNrOfNodes() == forceBalancer2.getNrOfLinks());
    CHECK(forceBalancer2.getNrOfNodes() == 4);
    CHECK(forceBalancer2.getNrOfSprings() == 5);
    double N = 23.;
    forceBalancer2.setSpringContourLengths(
      Eigen::VectorXd::Constant(forceBalancer2.getNetwork().nrOfSprings, N));
    Eigen::VectorXd springPartitions0 = forceBalancer2.getSpringPartitions();
    pcm::ForceBalanceNetwork net0 = forceBalancer2.getNetwork();
    for (int i = 0; i < net0.nrOfPartialSprings; i++) {
      std::cout << net0.springPartIndexA[i] << ", " << net0.springPartIndexB[i]
                << ": ";
      std::cout << springPartitions0[i] << std::endl;
    }
    std::cout << std::endl;

    SECTION("Irrelevant slip-links are rationalised")
    {
      forceBalancer2.configAssumeBoxLargeEnough(true);
      // check throws when adding slip-links
      CHECK_THROWS(forceBalancer2.addSlipLinks({ { 1, 2, 3 } },
                                               { { 1, 2, 3 } },
                                               { { 1., 2., 3. } },
                                               { { 1., 2., 3. } },
                                               { { 1., 2. } }));
      CHECK_THROWS(forceBalancer2.addSlipLinks({ { 1, 2, 3 } },
                                               { { 1, 2 } },
                                               { { 1., 2., 3. } },
                                               { { 1., 2., 3. } },
                                               { { 1., 2., 3.0 } }));
      // and actually add slip-links
      CHECK_NOTHROW(forceBalancer2.addSlipLinks({ { 0, 3 } },
                                                { { 1, 3 } },
                                                { { 4.2, 1.3 } },
                                                { { 3.9, 1.2 } },
                                                { { 1.3, 1.1 } }));
      CHECK(forceBalancer2.getNrOfNodes() + 2 == forceBalancer2.getNrOfLinks());
      // and run with them.
      // Expect the slip-link of between two strands to converge to the central
      // atom Expect the slip-link of two strands to stay at 0.5, 0.5.
      Eigen::VectorXd springPartitions = forceBalancer2.getSpringPartitions();
      Eigen::VectorXd displacements =
        Eigen::VectorXd::Zero(forceBalancer2.getNrOfLinks() * 3);
      pcm::ForceBalanceNetwork net = forceBalancer2.getNetwork();
      // outputNetwork(net, displacements, springPartitions);

      forceBalancer2.displaceToMeanPosition(
        net, displacements, springPartitions, 4);
      CHECK(net.coordinates[4 * 3] + displacements[4 * 3] ==
            Catch::Approx(2.5));
      CHECK(net.coordinates[4 * 3 + 1] + displacements[4 * 3 + 1] ==
            Catch::Approx(2.5));
      CHECK(net.coordinates[4 * 3 + 2] + displacements[4 * 3 + 2] ==
            Catch::Approx(2));

      // reset
      displacements.setZero();

      CHECK(net.springIndicesOfLinks.size() == net.nrOfLinks);
      for (int j = 0; j < 125; ++j) {
        forceBalancer2.displaceToMeanPosition(
          net, displacements, springPartitions, 5, -1.);
        forceBalancer2.updateSpringPartition(
          net, displacements, springPartitions, 5, -1.);
      }
      for (int i = 0; i < 125; ++i) {
        // do some random 250 steps with these two slip-links
        // NOTE: difficulty: finding out which node and spring it is actually
        // after the removal of strand atoms
        forceBalancer2.displaceToMeanPosition(
          net, displacements, springPartitions, 4, 1. / (i + 1));
        forceBalancer2.updateSpringPartition(
          net, displacements, springPartitions, 4, 1. / (i + 1));
      }
      // assert expectations are met.
      // NOTE: difficulty: finding out which spring idx it actually is
      // outputNetwork(net, displacements, springPartitions);
      CHECK_THAT(springPartitions[5],
                 Catch::Matchers::WithinAbs(0., 1e-3)); // 4-1(2)
      CHECK_THAT(springPartitions[6],
                 Catch::Matchers::WithinAbs(0., 1e-3)); // 4-2(3)
      CHECK_THAT(springPartitions[0],
                 Catch::Matchers::WithinAbs(1., 1e-3)); // 0(1)-4
      CHECK_THAT(springPartitions[1],
                 Catch::Matchers::WithinAbs(1., 1e-3)); // 1(2)-4
      // CHECK(springPartitions[8] == Catch::Approx(1.0).margin(1e-6)); // 5-3
      // CHECK(springPartitions[7] + 1e-5 == Catch::Approx(0.0 +
      // 1e-5).margin(1e-6)); // 5-5 CHECK(springPartitions[3] ==
      // Catch::Approx(1.0).margin(1e-6)); // 5-0
    }
  }
}

TEST_CASE("MEHP Force Balance runs with non-network",
          "[analysis][MEHPForceBalance][NonGaussianSpringForceEvaluator][long]")
{
  std::cout << "Running test \"MEHP Force Balance runs with non-network\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("MEHP Force Balance 3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      pcm::MEHPForceBalance forceBalancer2 =
        pcm::MEHPForceBalance(universe2, 2, false);
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      CHECK_NOTHROW(forceBalancer2.runForceRelaxation());
      CHECK(forceBalancer2.getNrOfIterations() > 1);
    }
  }
}

TEST_CASE("MEHP Force Balance Free chains collapse",
          "[analysis][MEHPForceBalance][NonGaussianSpringForceEvaluator]["
          "SimpleSpringMEHPForceEvaluator]")
{
  // TODO: implement the NonGaussianSpringForceEvaluator
  size_t nrOfBeads = 30;
  size_t nrOfBeadsPerChain = 3;
  pe::Universe universe =
    pe::Universe(nrOfBeads * 10.0, nrOfBeads * 10.0, nrOfBeads * 10.0);
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
  for (int i = 0; i < nrOfBeads; ++i) {
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i);
    atomTypes.push_back(i % nrOfBeadsPerChain == 0 ? 2 : 1);
    zeroInts.push_back(0);
    if (i > 0) {
      bondFrom.push_back(i - 1);
      bondTo.push_back(i);
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

  // now, check for every force evaluator, that the maximum entropy is when all
  // these beads overlap first, the gaussian spring one
  pcm::MEHPForceBalance forceBalancer =
    pcm::MEHPForceBalance(universe, 2, false);
  CHECK(forceBalancer.getNrOfSprings() ==
        forceBalancer.getNrOfPartialSprings());
  CHECK(forceBalancer.getNrOfSprings() == nrOfBeads / nrOfBeadsPerChain);

  SECTION("Large enough box")
  {
    forceBalancer.configAssumeBoxLargeEnough(true);
    CHECK_NOTHROW(forceBalancer.runForceRelaxation(50000, 1e-18));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
    CHECK(forceBalancer.getAverageSpringLength() >= 0.0);
    CHECK(forceBalancer.getAverageSpringLength() <= 3e-6);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }

  SECTION("Not large enough box")
  {
    forceBalancer.configAssumeBoxLargeEnough(false);
    CHECK_NOTHROW(forceBalancer.runForceRelaxation(50000, 1e-18));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
    CHECK(forceBalancer.getAverageSpringLength() >= 0.0);
    CHECK(forceBalancer.getAverageSpringLength() <= 3e-6);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }
}

TEST_CASE("MEHP Force Balance Fully active chains are fully active",
          "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHP Force Balance Fully active chains are "
               "fully active\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

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

      pcm::MEHPForceBalance forceBalancer =
        pcm::MEHPForceBalance(universe, 2, false);

      SECTION("For large enough box")
      {
        forceBalancer.configAssumeBoxLargeEnough(true);
        CHECK(forceBalancer.getNrOfActiveSprings(0.001) ==
              forceBalancer.getNrOfSprings());
        double initialResidual = forceBalancer.getDisplacementResidualNorm();
        CHECK(std::isfinite(initialResidual));
        CHECK_NOTHROW(forceBalancer.runForceRelaxation(10000, 1e-12));
        CHECK(forceBalancer.getNrOfIterations() > 0);
        CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
        CHECK(forceBalancer.getNrOfActiveSprings(0.001) ==
              forceBalancer.getNrOfSprings());
        CHECK(initialResidual > forceBalancer.getDisplacementResidualNorm());
      }

      SECTION("For not large enough box")
      {
        forceBalancer.configAssumeBoxLargeEnough(false);
        CHECK(forceBalancer.getNrOfActiveSprings(0.001) ==
              forceBalancer.getNrOfSprings());
        double initialResidual = forceBalancer.getDisplacementResidualNorm();
        CHECK(std::isfinite(initialResidual));
        CHECK_NOTHROW(forceBalancer.runForceRelaxation(10000, 1e-12));
        CHECK(forceBalancer.getNrOfIterations() > 0);
        CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
        CHECK(forceBalancer.getNrOfActiveSprings(0.001) ==
              forceBalancer.getNrOfSprings());
        CHECK(initialResidual > forceBalancer.getDisplacementResidualNorm());
      }
    }
  }
}

TEST_CASE("MEHP Force Balance Gives Identical Results for Different PBC "
          "imperfect Diamond Network",
          "[analysis][MEHPForceBalance]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  // perfect diamond network = fully connected =>
  // maximum is at perfect crystal structure -> must be all active.
  std::string inputFile =
    suspectedPath +
    "3d-diamond-lattice_5x5x5_a_3_d_0.85_imperfect.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPForceBalance forceBalanceConventional =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 250, 2.0, 100, 2.0, "a533d", 2, false);
    forceBalanceConventional.configAssumeBoxLargeEnough(true);

    pcm::MEHPForceBalance forceBalanceNew =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 250, 2.0, 100, 2.0, "a533d", 2, false);
    forceBalanceNew.configAssumeBoxLargeEnough(false);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));

    double initialResidual =
      forceBalanceConventional.getDisplacementResidualNorm();
    CHECK(initialResidual ==
          Catch::Approx(forceBalanceNew.getDisplacementResidualNorm()));

    forceBalanceConventional.runForceRelaxation(250, 1e-4);
    forceBalanceNew.runForceRelaxation(250, 1e-4);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getDisplacementResidualNorm() ==
          Catch::Approx(forceBalanceNew.getDisplacementResidualNorm()));
    CHECK_THAT(forceBalanceConventional.getIdsOfActiveNodes(),
               Catch::Matchers::Equals(forceBalanceNew.getIdsOfActiveNodes()));
    CHECK(forceBalanceConventional.getNrOfPartialSprings() ==
          forceBalanceNew.getNrOfPartialSprings());
    CHECK(forceBalanceConventional.getCurrentDisplacements().isApprox(
      forceBalanceNew.getCurrentDisplacements()));

    forceBalanceConventional.runForceRelaxation(
      250,
      1e-8,
      initialResidual,
      pcm::StructureSimplificationMode::ALL_TIM,
      0.01,
      false,
      pcm::LinkSwappingMode::NO_SWAPPING,
      10,
      1.,
      1);
    forceBalanceNew.runForceRelaxation(
      250,
      1e-8,
      initialResidual,
      pcm::StructureSimplificationMode::ALL_TIM,
      0.01,
      false,
      pcm::LinkSwappingMode::NO_SWAPPING,
      10,
      1.,
      1);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getDisplacementResidualNorm() ==
          Catch::Approx(forceBalanceNew.getDisplacementResidualNorm()));
    CHECK(forceBalanceConventional.getNrOfPartialSprings() ==
          forceBalanceNew.getNrOfPartialSprings());
    CHECK_THAT(forceBalanceConventional.getIdsOfActiveNodes(),
               Catch::Matchers::Equals(forceBalanceNew.getIdsOfActiveNodes()));
    CHECK(forceBalanceConventional.getNrOfPartialSprings() ==
          forceBalanceNew.getNrOfPartialSprings());
    CHECK(forceBalanceConventional.getCurrentDisplacements().isApprox(
      forceBalanceNew.getCurrentDisplacements()));

    forceBalanceConventional.runForceRelaxation(
      250,
      1e-8,
      initialResidual,
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION,
      -1.,
      false,
      pcm::LinkSwappingMode::ALL_MC,
      10,
      1.,
      0);
    forceBalanceNew.runForceRelaxation(
      250,
      1e-8,
      initialResidual,
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION,
      -1.,
      false,
      pcm::LinkSwappingMode::ALL_MC,
      10,
      1.,
      0);

    CHECK(forceBalanceConventional.getNrOfPartialSprings() ==
          forceBalanceNew.getNrOfPartialSprings());
    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getDisplacementResidualNorm() ==
          Catch::Approx(forceBalanceNew.getDisplacementResidualNorm()));
    CHECK(forceBalanceConventional.getCurrentDisplacements().isApprox(
      forceBalanceNew.getCurrentDisplacements()));
    CHECK_THAT(forceBalanceConventional.getIdsOfActiveNodes(),
               Catch::Matchers::Equals(forceBalanceNew.getIdsOfActiveNodes()));
  }
}

TEST_CASE("MEHP Force Balance Gives Identical Results for Different PBC p = 1",
          "[analysis][MEHPForceBalance][long]")
{
  std::cout << "Running test \"MEHP Force Balance Gives Identical Results for "
               "Different PBC p = 1\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  // perfect diamond network = fully connected =>
  // maximum is at perfect crystal structure -> must be all active.
  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPForceBalance forceBalanceConventional =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 25, 2.0, 20, 2.0, "a533d", 2, false);
    forceBalanceConventional.configAssumeBoxLargeEnough(true);

    pcm::MEHPForceBalance forceBalanceNew =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 25, 2.0, 20, 2.0, "a533d", 2, false);
    forceBalanceNew.configAssumeBoxLargeEnough(false);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));

    double initialResidual =
      forceBalanceConventional.getDisplacementResidualNorm();
    CHECK(initialResidual ==
          Catch::Approx(forceBalanceNew.getDisplacementResidualNorm()));

    forceBalanceConventional.runForceRelaxation(250, 1e-4);
    forceBalanceNew.runForceRelaxation(250, 1e-4);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getDisplacementResidualNorm() ==
          Catch::Approx(forceBalanceNew.getDisplacementResidualNorm()));
    CHECK_THAT(forceBalanceConventional.getIdsOfActiveNodes(),
               Catch::Matchers::Equals(forceBalanceNew.getIdsOfActiveNodes()));

    forceBalanceConventional.runForceRelaxation(
      250,
      1e-8,
      initialResidual,
      pcm::StructureSimplificationMode::ALL_TIM,
      0.01,
      false,
      pcm::LinkSwappingMode::NO_SWAPPING,
      10,
      1.,
      1);
    forceBalanceNew.runForceRelaxation(
      250,
      1e-8,
      initialResidual,
      pcm::StructureSimplificationMode::ALL_TIM,
      0.01,
      false,
      pcm::LinkSwappingMode::NO_SWAPPING,
      10,
      1.,
      1);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getDisplacementResidualNorm() ==
          Catch::Approx(forceBalanceNew.getDisplacementResidualNorm()));
    CHECK_THAT(forceBalanceConventional.getIdsOfActiveNodes(),
               Catch::Matchers::Equals(forceBalanceNew.getIdsOfActiveNodes()));

    forceBalanceConventional.runForceRelaxation(
      250,
      1e-8,
      initialResidual,
      pcm::StructureSimplificationMode::ALL_TIM,
      0.01,
      false,
      pcm::LinkSwappingMode::ALL_MC,
      10,
      1.,
      3);
    forceBalanceNew.runForceRelaxation(
      250,
      1e-8,
      initialResidual,
      pcm::StructureSimplificationMode::ALL_TIM,
      0.01,
      false,
      pcm::LinkSwappingMode::ALL_MC,
      10,
      1.,
      3);

    CHECK(forceBalanceConventional.getStressTensor().isApprox(
      forceBalanceNew.getStressTensor()));
    CHECK(forceBalanceConventional.getDisplacementResidualNorm() ==
          Catch::Approx(forceBalanceNew.getDisplacementResidualNorm()));
    CHECK_THAT(forceBalanceConventional.getIdsOfActiveNodes(),
               Catch::Matchers::Equals(forceBalanceNew.getIdsOfActiveNodes()));
  } else {
    std::cout << "File " << inputFile << " does not exists. Skipping test."
              << std::endl;
  }
}

TEST_CASE("MEHP Force Balance does not collapse",
          "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHP Force Balance does not collapse\""
            << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  /**
   * @brief A grid of two rows, each one bead between the two crosslinkers
   *
   */
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

  SECTION("Assuming large enough box")
  {
    pcm::MEHPForceBalance forceBalanceConventional =
      pcm::MEHPForceBalance(universe, 2, true);
    forceBalanceConventional.configAssumeBoxLargeEnough(true);
    CHECK_NOTHROW(forceBalanceConventional.runForceRelaxation(
      10000, 1e-15, -1, pcm::StructureSimplificationMode::ALL_TIM));
    CHECK(forceBalanceConventional.getNrOfIterations() > 0);
    CHECK(forceBalanceConventional.getExitReason() ==
          pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalanceConventional.getNrOfActiveSprings() ==
          forceBalanceConventional.getNrOfSprings());
    // compare to what we expect
    CHECK(forceBalanceConventional.getNrOfActiveSprings() == 8);
    CHECK(forceBalanceConventional.getNrOfActiveNodes() == 4);
    // CHECK(forceBalanceConventional.getAverageSpringLength() ==
    //       Catch::Approx(5.));
    // CHECK_THAT(forceBalanceConventional.getGammaFactor(),
    //            Catch::Matchers::WithinAbs(1.0, 1e-3));
    // outputNetwork(forceBalanceConventional.getNetwork(),
    //               forceBalanceConventional.getCurrentDisplacements(),
    //               forceBalanceConventional.getSpringPartitions());
  }
  SECTION("Assuming too small box")
  {
    pcm::MEHPForceBalance forceBalanceNew =
      pcm::MEHPForceBalance(universe, 2, true);
    forceBalanceNew.configAssumeBoxLargeEnough(false);
    CHECK_NOTHROW(forceBalanceNew.runForceRelaxation(
      10000, 1e-15, -1, pcm::StructureSimplificationMode::ALL_TIM));
    CHECK(forceBalanceNew.getNrOfIterations() > 0);
    CHECK(forceBalanceNew.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalanceNew.getNrOfActiveSprings() ==
          forceBalanceNew.getNrOfSprings());
    // compare to what we expect
    CHECK(forceBalanceNew.getNrOfActiveSprings() == 8);
    CHECK(forceBalanceNew.getNrOfActiveNodes() == 4);
    CHECK(forceBalanceNew.getAverageSpringLength() == Catch::Approx(5.0));
    // forceBalanceNew.setSpringContourLengths(
    //   Eigen::VectorXd::Constant(forceBalanceNew.getNrOfSprings(), 5.));
    CHECK(forceBalanceNew.getDefaultNrOfChains() == 8);
    // TODO: check this again
    CHECK_THAT(forceBalanceNew.getGammaFactor(1., 100),
               Catch::Matchers::WithinAbs(1.0, 1e-2));
    CHECK_THAT(forceBalanceNew.getPressure(),
               Catch::Matchers::WithinAbs(1. / 30., 1e-3));
    // outputNetwork(forceBalanceNew.getNetwork(),
    //               forceBalanceNew.getCurrentDisplacements(),
    //               forceBalanceNew.getSpringPartitions());
  }
};

TEST_CASE(
  "MEHP Force Balance correctly collapses melts even with random slip-links",
  "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHP Force Balance correctly collapses melts "
               "even with random slip-links\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  std::string inputFile = suspectedPath + "melt_83_a_100.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 1000, 2.0, 100, 5, "my_seed_fb12");
    forceBalancer.configAssumeBoxLargeEnough(false);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    double initialResidual = forceBalancer.getDisplacementResidualNorm();
    CHECK(std::isfinite(initialResidual));
    CHECK(forceBalancer.getNumExtraAtoms() > 100);
    CHECK(forceBalancer.getNetwork().nrOfPartialSprings >
          forceBalancer.getNetwork().nrOfSprings);
    CHECK_NOTHROW(forceBalancer.runForceRelaxation(
      50000, 1e-9, -1., pcm::StructureSimplificationMode::ALL_TIM));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
    CHECK(initialResidual > forceBalancer.getDisplacementResidualNorm());

    double pressBefore = forceBalancer.getPressure();
    forceBalancer.configAssumeBoxLargeEnough(true);
    CHECK(pressBefore >= forceBalancer.getPressure());
  }
}

TEST_CASE("MEHP Force Balance correctly re-aligns Slip-Links to Images",
          "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHP Force Balance correctly re-aligns "
               "Slip-Links to Images\""
            << std::endl;

  // it does not really matter with what structure we start...
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);

  pcm::MEHPForceBalance forceBalancer =
    pcm::MEHPForceBalance(universe, 2, false);
  forceBalancer.configAssumeBoxLargeEnough(false);
  CHECK(forceBalancer.getPressure() == Catch::Approx(0.0));

  // ...if we invoke the relevant method with our custom network
  pcm::ForceBalanceNetwork net;
  net.L[0] = 10.;
  net.L[1] = 10.;
  net.L[2] = 10.;
  net.coordinates = Eigen::VectorXd::Zero(3 * 3);
  net.springPartIndexA = Eigen::ArrayXi::Zero(2);
  net.springPartIndexB = Eigen::ArrayXi::Zero(2);
  net.springPartBoxOffset = Eigen::VectorXd::Zero(2 * 3);
  net.partialToFullSpringIndex = Eigen::ArrayXi::Zero(2);
  net.linkIsSliplink = Eigen::ArrayXb::Constant(3, false);

  SECTION("System 1")
  {
    net.springPartIndexA << 0, 2;
    net.springPartIndexB << 2, 1;
    net.linkIsSliplink[2] = true;

    Eigen::VectorXd u = Eigen::VectorXd::Zero(net.coordinates.size());

    CHECK_NOTHROW(forceBalancer.reAlignSlipLinkToImages(net, u, 2, 0, 1));
    CHECK((net.springPartBoxOffset.array() == (0.0)).all());

    net.coordinates << 5, 2.5, 0., // cross-link 1,
      5., 7.5, 0.,                 // cross-link 2,
      5., 5., 0.;                  // slip-link

    net.springPartBoxOffset << 0., 0., 0., // cross-link 1 -> slip-link
      20., 0., 0.;                         // slip-link -> cross-link 2

    CHECK_NOTHROW(forceBalancer.reAlignSlipLinkToImages(net, u, 2, 0, 1));
    CHECK(net.coordinates.segment(3 * 2, 3)[1] == 5.);

    CHECK(net.springPartBoxOffset[0] == Catch::Approx(10.));
    CHECK(net.springPartBoxOffset[3] == Catch::Approx(10.));
    for (size_t i = 0; i < 2 * 3; ++i) {
      if (i % 3 != 0) {
        CHECK(net.springPartBoxOffset[i] == Catch::Approx(0.));
      }
    }
  }

  SECTION("System 2")
  {
    net.springPartIndexA << 1, 2;
    net.springPartIndexB << 2, 0;
    net.linkIsSliplink[2] = true;

    Eigen::VectorXd u = Eigen::VectorXd::Zero(net.coordinates.size());

    CHECK_NOTHROW(forceBalancer.reAlignSlipLinkToImages(net, u, 2, 0, 1));
    CHECK((net.springPartBoxOffset.array() == (0.0)).all());

    net.coordinates << 5, 2.5, 0., // cross-link 1,
      5., 7.5, 0.,                 // cross-link 2,
      5., 7.5, 0.;                 // slip-link

    net.springPartBoxOffset << 20., 0., 0., // cross-link 2 -> slip-link
      0., 0., 0.;                           // slip-link -> cross-link 1

    CHECK_NOTHROW(forceBalancer.reAlignSlipLinkToImages(net, u, 2, 0, 1));

    Eigen::VectorXd expectation = Eigen::VectorXd::Zero(2 * 3);
    expectation(0) = 10.;
    expectation(3) = 10.;

    for (size_t i = 0; i < 2 * 3; ++i) {
      CHECK(expectation[i] == Catch::Approx(net.springPartBoxOffset[i]));
    }
  }

  SECTION("System 3")
  {
    net.springPartIndexA << 2, 0;
    net.springPartIndexB << 1, 2;
    net.linkIsSliplink[2] = true;

    Eigen::VectorXd u = Eigen::VectorXd::Zero(net.coordinates.size());

    CHECK_NOTHROW(forceBalancer.reAlignSlipLinkToImages(net, u, 2, 1, 0));
    CHECK((net.springPartBoxOffset.array() == (0.0)).all());

    net.coordinates << 5, 2.5, 0., // cross-link 1,
      5., 7.5, 0.,                 // cross-link 2,
      5., 2.5, 0.;                 // slip-link

    net.springPartBoxOffset << 10., 0., 0., // slip-link -> cross-link 2
      -10., 0., 0.;                         // slip-link -> cross-link 1

    CHECK_NOTHROW(forceBalancer.reAlignSlipLinkToImages(net, u, 2, 1, 0));

    Eigen::VectorXd expectation = Eigen::VectorXd::Zero(2 * 3);
    expectation(0) = 0.;
    expectation(3) = 0.;

    for (size_t i = 0; i < 2 * 3; ++i) {
      CHECK(expectation[i] == Catch::Approx(net.springPartBoxOffset[i]));
    }
  }
}

TEST_CASE("Particular MEHP Force Balance Example",
          "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"Particular MEHP Force Balance Example\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure.out-"
    "equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, false, true, false);
    forceBalancer.configAssumeBoxLargeEnough(true);
    CHECK_NOTHROW(forceBalancer.runForceRelaxation());

    double pressBefore = forceBalancer.getPressure();
    forceBalancer.configAssumeBoxLargeEnough(false);
    CHECK(pressBefore <= forceBalancer.getPressure());
  } else {
    std::cerr << "File " << inputFile << " does not exist." << std::endl;
  }
}

TEST_CASE("Random sampling example", "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"Random sampling example\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    // generate the same slip-links twice,
    // once for each assumption
    pcm::MEHPForceBalance forceBalancerLargeOldSampling =
      pcm::MEHPForceBalance(universe, 2, false, false, false);
    forceBalancerLargeOldSampling.configAssumeBoxLargeEnough(true);
    forceBalancerLargeOldSampling.randomlyAddSliplinks(
      1000, 6.0, 900, 3.0, false, 53467829);

    pcm::MEHPForceBalance forceBalancerSmallOldSampling1 =
      pcm::MEHPForceBalance(universe, 2, false, false, false);

    pcm::MEHPForceBalance forceBalancerSmallOldSampling2 =
      pcm::MEHPForceBalance(universe, 2, false, false, false);

    // initially
    CHECK_THAT(forceBalancerSmallOldSampling1.getDisplacementResidualNorm(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling2.getDisplacementResidualNorm()));

    // check that the order does not matter
    forceBalancerSmallOldSampling1.configAssumeBoxLargeEnough(false);
    forceBalancerSmallOldSampling1.randomlyAddSliplinks(
      1000, 6.0, 900, 3.0, false, 53467829);
    forceBalancerSmallOldSampling2.randomlyAddSliplinks(
      1000, 6.0, 900, 3.0, false, 53467829);
    forceBalancerSmallOldSampling2.configAssumeBoxLargeEnough(false);

    // after adding slip-links
    CHECK(forceBalancerLargeOldSampling.getPressure() <=
          forceBalancerSmallOldSampling1.getPressure());
    CHECK_THAT(forceBalancerSmallOldSampling1.getDisplacementResidualNorm(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling2.getDisplacementResidualNorm()));
    CHECK_THAT(
      forceBalancerSmallOldSampling1.getPressure(),
      Catch::Matchers::WithinRel(forceBalancerSmallOldSampling2.getPressure()));

    // it is actually thinkable that the following fails for certain scenarios.
    // however, in general, it should not
    CHECK(forceBalancerLargeOldSampling.getDisplacementResidualNorm(1.) <=
          forceBalancerSmallOldSampling1.getDisplacementResidualNorm(1.));
    CHECK(forceBalancerLargeOldSampling.getDisplacementResidualNorm(-1.) <=
          forceBalancerSmallOldSampling1.getDisplacementResidualNorm(-1.));

    // check to make sure the slip-links are actually placed identically
    CHECK(forceBalancerLargeOldSampling.getNetwork().springPartIndexA.isApprox(
      forceBalancerSmallOldSampling1.getNetwork().springPartIndexA));
    CHECK(forceBalancerSmallOldSampling1.getNetwork().springPartIndexB.isApprox(
      forceBalancerSmallOldSampling2.getNetwork().springPartIndexB));

    pcm::MEHPForceBalance forceBalancerSmallNewSampling =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 1000, 6.0, 900, 3.0, "53467829");
    forceBalancerSmallNewSampling.configAssumeBoxLargeEnough(false);

    CHECK_THAT(forceBalancerSmallNewSampling.getNrOfPartialSprings(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling1.getNrOfPartialSprings(), 0.1));
    CHECK_THAT(forceBalancerSmallNewSampling.getPressure(),
               Catch::Matchers::WithinRel(
                 forceBalancerSmallOldSampling1.getPressure(), 0.5));
    CHECK_THAT(
      forceBalancerSmallNewSampling.getDisplacementResidualNorm(),
      Catch::Matchers::WithinRel(
        forceBalancerSmallOldSampling1.getDisplacementResidualNorm(), 0.5));

    // check that the spring vectors are equal
    Eigen::VectorXd springVectorsSmallOldSampling =
      forceBalancerSmallOldSampling1.evaluateSpringVectors(
        forceBalancerSmallOldSampling1.getNetwork(),
        forceBalancerSmallOldSampling1.getCurrentDisplacements());
    Eigen::VectorXd springVectorsSmallNewSampling =
      forceBalancerSmallNewSampling.evaluateSpringVectors(
        forceBalancerSmallNewSampling.getNetwork(),
        forceBalancerSmallNewSampling.getCurrentDisplacements());
    REQUIRE(forceBalancerSmallOldSampling1.getNetwork().nrOfSprings ==
            forceBalancerSmallNewSampling.getNetwork().nrOfSprings);

    CHECK(
      springVectorsSmallOldSampling.isApprox(springVectorsSmallNewSampling));

    forceBalancerSmallNewSampling.configAssumeBoxLargeEnough(true);

    CHECK_THAT(forceBalancerSmallNewSampling.getPressure(),
               Catch::Matchers::WithinRel(
                 forceBalancerLargeOldSampling.getPressure(), 0.5));
    CHECK_THAT(
      forceBalancerSmallNewSampling.getDisplacementResidualNorm(),
      Catch::Matchers::WithinRel(
        forceBalancerLargeOldSampling.getDisplacementResidualNorm(), 0.5));

    CHECK_THAT(
      forceBalancerLargeOldSampling.getForceMagnitudeVector().mean(),
      Catch::Matchers::WithinRel(
        Eigen::median(forceBalancerLargeOldSampling.getForceMagnitudeVector()),
        0.25));
    CHECK_THAT(
      forceBalancerSmallOldSampling1.getForceMagnitudeVector().mean(),
      Catch::Matchers::WithinRel(
        Eigen::median(forceBalancerSmallOldSampling1.getForceMagnitudeVector()),
        0.25));
    CHECK_THAT(
      forceBalancerSmallOldSampling2.getForceMagnitudeVector().mean(),
      Catch::Matchers::WithinRel(
        Eigen::median(forceBalancerSmallOldSampling2.getForceMagnitudeVector()),
        0.25));
    CHECK_THAT(
      forceBalancerSmallNewSampling.getForceMagnitudeVector().mean(),
      Catch::Matchers::WithinRel(
        Eigen::median(forceBalancerSmallNewSampling.getForceMagnitudeVector()),
        0.25));

    // DEBUG:
    // Eigen::VectorXd magnitudeVector =
    //   forceBalancerSmallOldSampling1.getForceMagnitudeVector();
    // size_t maxIndex;
    // double maxValue = magnitudeVector.maxCoeff(&maxIndex);
  } else {
    std::cerr << "File " << inputFile << " does not exist." << std::endl;
  }
}

TEST_CASE("Random sampling example small", "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"Random sampling example small\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  std::string inputFile =
    suspectedPath + "square_lattice_2x2_a_5.2d.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    // generate the same slip-links twice,
    // once for each assumption
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, true, true, false);
    forceBalancer.configAssumeBoxLargeEnough(true);
    forceBalancer.randomlyAddSliplinks(12, 6.0, 11, 3.0, false, 86573452);

    pcm::MEHPForceBalance forceBalancer2 =
      pcm::MEHPForceBalance(universe, 2, true, true, false);

    pcm::MEHPForceBalance forceBalancer3 =
      pcm::MEHPForceBalance(universe, 2, true, true, false);

    // initially
    CHECK_THAT(
      forceBalancer2.getDisplacementResidualNorm(),
      Catch::Matchers::WithinRel(forceBalancer3.getDisplacementResidualNorm()));

    // check that the order does not matter
    forceBalancer2.configAssumeBoxLargeEnough(false);
    forceBalancer2.randomlyAddSliplinks(12, 6.0, 11, 3.0, false, 86573452);
    forceBalancer3.randomlyAddSliplinks(12, 6.0, 11, 3.0, false, 86573452);
    forceBalancer3.configAssumeBoxLargeEnough(false);

    // after adding slip-links
    CHECK(forceBalancer.getPressure() <= forceBalancer2.getPressure());
    CHECK_THAT(
      forceBalancer2.getDisplacementResidualNorm(),
      Catch::Matchers::WithinRel(forceBalancer3.getDisplacementResidualNorm()));
    CHECK_THAT(forceBalancer2.getPressure(),
               Catch::Matchers::WithinRel(forceBalancer3.getPressure()));
    // it is actually thinkable that the following fails for certain scenarios.
    // however, in general, it should not
    CHECK(forceBalancer.getDisplacementResidualNorm(1.) <=
          forceBalancer2.getDisplacementResidualNorm(1.));
    CHECK(forceBalancer.getDisplacementResidualNorm(-1.) <=
          forceBalancer2.getDisplacementResidualNorm(-1.));

    // check to make sure the slip-links are actually placed identically
    CHECK(forceBalancer.getNetwork().springPartIndexA.isApprox(
      forceBalancer2.getNetwork().springPartIndexA));
    CHECK(forceBalancer2.getNetwork().springPartIndexB.isApprox(
      forceBalancer3.getNetwork().springPartIndexB));

    pcm::MEHPForceBalance forceBalancer4 =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 12, 6.0, 11, 3.0, "86573452", 2, true);
    forceBalancer4.configAssumeBoxLargeEnough(false);

    CHECK_THAT(forceBalancer4.getPressure(),
               Catch::Matchers::WithinRel(forceBalancer2.getPressure(), 0.1));
    CHECK_THAT(forceBalancer4.getDisplacementResidualNorm(),
               Catch::Matchers::WithinRel(
                 forceBalancer2.getDisplacementResidualNorm(), 0.1));

    forceBalancer4.configAssumeBoxLargeEnough(true);

    CHECK_THAT(forceBalancer4.getPressure(),
               Catch::Matchers::WithinRel(forceBalancer.getPressure(), 0.1));
    CHECK_THAT(forceBalancer4.getDisplacementResidualNorm(),
               Catch::Matchers::WithinRel(
                 forceBalancer.getDisplacementResidualNorm(), 0.1));

    // std::cout << "FB 1:" << std::endl;
    // outputNetwork(forceBalancer.getNetwork(),
    //               forceBalancer.getCurrentDisplacements(),
    //               forceBalancer.getSpringPartitions());
    // std::cout << "FB 2:" << std::endl;
    // outputNetwork(forceBalancer2.getNetwork(),
    //               forceBalancer2.getCurrentDisplacements(),
    //               forceBalancer2.getSpringPartitions());
    // std::cout << "FB 3:" << std::endl;
    // outputNetwork(forceBalancer3.getNetwork(),
    //               forceBalancer3.getCurrentDisplacements(),
    //               forceBalancer3.getSpringPartitions());
    // std::cout << "FB 4:" << std::endl;
    // outputNetwork(forceBalancer4.getNetwork(),
    //               forceBalancer4.getCurrentDisplacements(),
    //               forceBalancer4.getSpringPartitions());
  } else {
    std::cerr << "File " << inputFile << " does not exist." << std::endl;
  }
}

TEST_CASE("Yet another sampling example", "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"Yet another sampling example\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

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
    pcm::MEHPForceBalance forceBalancerNewSampling =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 1000, 6.0, 900, 3.0, "53467829");
    forceBalancerNewSampling.configAssumeBoxLargeEnough(false);

    // extract these randomly sampled slip-links
    std::vector<size_t> strandIdx1;
    std::vector<size_t> strandIdx2;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> alpha1;
    std::vector<double> alpha2;

    pcm::ForceBalanceNetwork net = forceBalancerNewSampling.getNetwork();
    for (size_t i = 0; i < net.nrOfLinks - net.nrOfNodes; ++i) {
      size_t linkIdx = i + net.nrOfNodes;
      strandIdx1.push_back(net.springIndicesOfLinks[linkIdx][0]);
      strandIdx2.push_back(
        pylimer_tools::utils::last(net.springIndicesOfLinks[linkIdx]));
      x.push_back(net.coordinates[3 * (linkIdx) + 0]);
      y.push_back(net.coordinates[3 * (linkIdx) + 1]);
      z.push_back(net.coordinates[3 * (linkIdx) + 2]);
      alpha1.push_back(forceBalancerNewSampling.sumToTotalFraction(
        net,
        forceBalancerNewSampling.getSpringPartitions(),
        net.springIndicesOfLinks[linkIdx][0],
        linkIdx));
      alpha2.push_back(forceBalancerNewSampling.sumToTotalFraction(
        net,
        forceBalancerNewSampling.getSpringPartitions(),
        pylimer_tools::utils::last(net.springIndicesOfLinks[linkIdx]),
        linkIdx));
      // TODO: handle these primary loops, then increase required accuracy
      // if (net.springIndicesOfLinks[linkIdx].size() == 1 ||
      //     net.springIndicesOfLinks[linkIdx][0] ==
      //       net.springIndicesOfLinks[linkIdx][1]) {
      //   alpha2[alpha2.size() - 1] +=
      // }
    }

    // generate the same slip-links twice,
    // once for each assumption
    pcm::MEHPForceBalance forceBalancerOldSamplingLarge =
      pcm::MEHPForceBalance(universe, 2, false, false, false);
    forceBalancerOldSamplingLarge.configAssumeBoxLargeEnough(true);

    // check that the general conversion is equal
    REQUIRE(forceBalancerOldSamplingLarge.getNetwork().springIndexA.isApprox(
      net.springIndexA));

    forceBalancerOldSamplingLarge.addSlipLinks(
      strandIdx1, strandIdx2, x, y, z, alpha1, alpha2, false);
    CHECK(forceBalancerOldSamplingLarge.getNrOfLinks() ==
          forceBalancerNewSampling.getNrOfLinks());
    CHECK(forceBalancerOldSamplingLarge.getNrOfPartialSprings() ==
          forceBalancerNewSampling.getNrOfPartialSprings());

    // std::cout << "\n\n\nnetwork 4:" << std::endl;
    // outputNetwork(forceBalancer4.getNetwork(),
    //               forceBalancer4.getCurrentDisplacements(),
    //               forceBalancer4.getSpringPartitions());
    // std::cout << "\n\n\nnetwork 1:" << std::endl;
    // outputNetwork(forceBalancer.getNetwork(),
    //               forceBalancer.getCurrentDisplacements(),
    //               forceBalancer.getSpringPartitions());

    pcm::MEHPForceBalance forceBalancerOldSamplingSmall =
      pcm::MEHPForceBalance(universe, 2, false, false, false);
    forceBalancerOldSamplingSmall.configAssumeBoxLargeEnough(false);
    forceBalancerOldSamplingSmall.addSlipLinks(
      strandIdx1, strandIdx2, x, y, z, alpha1, alpha2, false);

    // after adding slip-links
    CHECK(forceBalancerOldSamplingLarge.getNetwork().springPartIndexA.isApprox(
      forceBalancerOldSamplingSmall.getNetwork().springPartIndexA));
    CHECK_THAT(
      forceBalancerNewSampling.getDisplacementResidualNorm(),
      Catch::Matchers::WithinRel(
        forceBalancerOldSamplingLarge.getDisplacementResidualNorm(), 1e-3));
    CHECK_THAT(
      forceBalancerOldSamplingLarge.getPressure(),
      Catch::Matchers::WithinRel(forceBalancerNewSampling.getPressure(), 1e-3));

    CHECK(forceBalancerOldSamplingLarge.getPressure() <=
          forceBalancerOldSamplingSmall.getPressure());
    // it is actually thinkable that the following fails for certain scenarios.
    // however, in general, it should not
    CHECK(forceBalancerOldSamplingLarge.getDisplacementResidualNorm(1.) <=
          forceBalancerOldSamplingSmall.getDisplacementResidualNorm(1.));
    CHECK(forceBalancerOldSamplingLarge.getDisplacementResidualNorm(-1.) <=
          forceBalancerOldSamplingSmall.getDisplacementResidualNorm(-1.));

    // check to make sure the slip-links are actually placed identically
    CHECK(forceBalancerOldSamplingLarge.getNetwork().springPartIndexA.isApprox(
      forceBalancerOldSamplingSmall.getNetwork().springPartIndexA));
    // cannot compare this: the indices are different depending on sampling
    // or addition method
    // CHECK(forceBalancer2.getNetwork().springPartIndexB.isApprox(
    //   forceBalancer4.getNetwork().springPartIndexB));
  } else {
    std::cerr << "File " << inputFile << " does not exist." << std::endl;
  }
}

TEST_CASE("Conversion of structure is equal for both methods",
          "[analysis][MEHPForceBalance]")
{
  std::cout
    << "Running test \"Conversion of structure is equal for both methods\""
    << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    // randomly sample NO slip-links
    pcm::MEHPForceBalance forceBalancerNewConversion =
      pcm::MEHPForceBalance::constructWithRandomSlipLinks(
        universe, 0, 6.0, 0, 3.0, "");
    forceBalancerNewConversion.configAssumeBoxLargeEnough(false);

    pcm::MEHPForceBalance forceBalancerOldConversion =
      pcm::MEHPForceBalance(universe, 2, false, false, false);
    forceBalancerOldConversion.configAssumeBoxLargeEnough(false);

    // check to make sure the links are actually placed identically
    CHECK(forceBalancerNewConversion.getNetwork().springPartIndexA.isApprox(
      forceBalancerOldConversion.getNetwork().springPartIndexA));
    CHECK(forceBalancerNewConversion.getNetwork().springPartIndexB.isApprox(
      forceBalancerOldConversion.getNetwork().springPartIndexB));
    CHECK(forceBalancerNewConversion.getNetwork().coordinates.isApprox(
      forceBalancerOldConversion.getNetwork().coordinates));
    CHECK(forceBalancerNewConversion.getNetwork().springPartBoxOffset.isApprox(
      forceBalancerOldConversion.getNetwork().springPartBoxOffset));

    // and of course the corresponding computations
    CHECK_THAT(forceBalancerOldConversion.getDisplacementResidualNorm(),
               Catch::Matchers::WithinRel(
                 forceBalancerNewConversion.getDisplacementResidualNorm()));
    CHECK_THAT(
      forceBalancerOldConversion.getPressure(),
      Catch::Matchers::WithinRel(forceBalancerNewConversion.getPressure()));

    // std::cout << "FB 1:" << std::endl;
    // outputNetwork(forceBalancerNewConversion.getNetwork(),
    //               forceBalancerNewConversion.getCurrentDisplacements(),
    //               forceBalancerNewConversion.getSpringPartitions());
    // std::cout << "\n\n\nFB 2:" << std::endl;
    // outputNetwork(forceBalancerOldConversion.getNetwork(),
    //               forceBalancerOldConversion.getCurrentDisplacements(),
    //               forceBalancerOldConversion.getSpringPartitions());
  }
}

TEST_CASE("Force Balance can be run without slipping")
{
  std::cout << "Running test \"Force Balance can be run without slipping\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPForceBalance forceBalancer =
    pcm::MEHPForceBalance::constructWithRandomSlipLinks(
      universe, 1000, 6.0, 900, 3.0, "53467829");
  forceBalancer.configAssumeBoxLargeEnough(false);

  Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
  forceBalancer.runForceRelaxation(
    500,
    1e-9,
    -1,
    pcm::StructureSimplificationMode::NO_SIMPLIFICATION,
    0.01,
    false,
    pcm::LinkSwappingMode::NO_SWAPPING,
    10,
    1.0,
    -1,
    true);
  CHECK(springPartitions.isApprox(forceBalancer.getSpringPartitions()));
}

TEST_CASE("Adding slip-links does not influence other springs",
          "[analysis][MEHPForceBalance]")
{
  std::cout
    << "Running test \"Adding slip-links does not influence other springs\""
    << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  std::string inputFile =
    suspectedPath +
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  // generate the same slip-links twice,
  // once for each assumption
  pcm::MEHPForceBalance forceBalancerOldSamplingSmall =
    pcm::MEHPForceBalance(universe, 2, false, false, false);
  forceBalancerOldSamplingSmall.configAssumeBoxLargeEnough(false);
  Eigen::VectorXd springVectorsWithoutSlipLinks =
    forceBalancerOldSamplingSmall.evaluateSpringVectors(
      forceBalancerOldSamplingSmall.getNetwork(),
      forceBalancerOldSamplingSmall.getCurrentDisplacements());
  forceBalancerOldSamplingSmall.randomlyAddSliplinks(
    1000, 6.0, 900, 3.0, false, 53467829);

  Eigen::VectorXd springVectorsWithSlipLinks =
    forceBalancerOldSamplingSmall.evaluateSpringVectors(
      forceBalancerOldSamplingSmall.getNetwork(),
      forceBalancerOldSamplingSmall.getCurrentDisplacements());

  CHECK(springVectorsWithoutSlipLinks.isApprox(springVectorsWithSlipLinks));

  // and the same for sampling method 2
  pcm::MEHPForceBalance forceBalancer2 =
    pcm::MEHPForceBalance::constructWithRandomSlipLinks(
      universe, 1000, 6.0, 900, 3.0, "53467829");
  forceBalancer2.configAssumeBoxLargeEnough(false);

  pcm::MEHPForceBalance forceBalancer2Without =
    pcm::MEHPForceBalance::constructWithRandomSlipLinks(
      universe, 0, 6.0, 0, 3.0, "53467829");
  forceBalancer2Without.configAssumeBoxLargeEnough(false);

  CHECK(forceBalancer2.getNrOfPartialSprings() >
        forceBalancer2Without.getNrOfPartialSprings());

  CHECK(forceBalancer2
          .evaluateSpringVectors(forceBalancer2.getNetwork(),
                                 forceBalancer2.getCurrentDisplacements())
          .isApprox(forceBalancer2Without.evaluateSpringVectors(
            forceBalancer2Without.getNetwork(),
            forceBalancer2Without.getCurrentDisplacements())));
}
