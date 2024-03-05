#include "../../src/pylimer_tools_cpp/calc/MEHPForceBalance2.h"
#include "../../src/pylimer_tools_cpp/calc/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/calc/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <random>
#include <vector>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::calc::mehp;

void
outputNetwork(pcm::ForceBalanceNetwork net, Eigen::VectorXd springPartitions)
{
  for (int i = 0; i < net.nrOfSprings; ++i) {
    std::cout << "Spring " << i << ", N: " << net.springsContourLength[i]
              << std::endl;
    for (int j = 0; j < net.linkIndicesOfSprings[i].size(); ++j) {
      std::cout << net.linkIndicesOfSprings[i][j] << ": ";
      for (int dir = 0; dir < 3; ++dir) {
        std::cout << (net.coordinates[3 * net.linkIndicesOfSprings[i][j] + dir])
                  << ", ";
      }
      std::cout << std::endl;
      if (j < net.linkIndicesOfSprings[i].size() - 1) {
        std::cout << net.localToGlobalSpringIndex.at(i)[j] << ": "
                  << springPartitions[net.localToGlobalSpringIndex.at(i)[j]]
                  << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

TEST_CASE("MEHP Force Balance 2 Particular slip-link examples",
          "[analysis][MEHPForceBalance2]")
{
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
  // slip-link 3 in the test-system
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
    pe::Box(lmda * L, L * (1. / sqrt(lmda)), L * (1. / sqrt(lmda))));

  pcm::MEHPForceBalance2 forceBalancer =
    pcm::MEHPForceBalance2::constructWithSlipLinks(
      universe,
      { 1, 0 },
      { 1, 1 },
      { 22.650980696700437, 15.67228278755659 },
      { 38.666887107085628, 3.172504627794056 },
      { 20.999752090751294, 27.126973321608833 },
      { 0.448276, 0.31034482758620685 },
      { 0.448276, 0.6896551724137931 * (0.448276) });

  // assert placement on the strands
  Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
  CHECK(springPartitions.size() == 6);
  CHECK(springPartitions[0] == Catch::Approx(0.0));
  CHECK(springPartitions[1] == Catch::Approx(1. - 0.448276));
  CHECK(springPartitions[2] == Catch::Approx(0.31034482758620685));
  CHECK(springPartitions[3] == Catch::Approx(1. - 0.31034482758620685));
  // CHECK(springPartitions[4] == Catch::Approx(0.3091558621));
  // CHECK(springPartitions[5] == Catch::Approx(0.1391201379));
  outputNetwork(forceBalancer.getNetwork(), springPartitions);
}

TEST_CASE("MEHP Force Balance 2 MC swap accept and reject work",
          "[analysis][MEHPForceBalance2]")
{
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

  SECTION("MC condition accepts as requested")
  {
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithSlipLinks(universe,
                                                     { 0, 0 },
                                                     { 1, 1 },
                                                     { 21., 23. },
                                                     { 18., 18. },
                                                     { 0.0, 0.0 },
                                                     { 0.6, 0.4 },
                                                     { 0.6, 0.4 });
    CHECK_THROWS(
      forceBalancer.setSpringContourLengths(Eigen::VectorXd::Constant(2, 30.)));
    forceBalancer.setSpringContourLengths(Eigen::VectorXd::Constant(3, 30.));

    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    long int edgeBetweenSlipLinks = -1;
    for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
      if (net.linkIsSliplink[net.springPartIndexA[i]] &&
          net.linkIsSliplink[net.springPartIndexB[i]]) {
        edgeBetweenSlipLinks = i;
      }
    }
    CHECK(edgeBetweenSlipLinks >= 0);
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    outputNetwork(net, partitions);
    CHECK(forceBalancer.swapSlipLinkReversibly(
      net, partitions, edgeBetweenSlipLinks, 1.));
    CHECK_FALSE(net.isUpToDate);
  }

  SECTION("MC condition rejects as requested")
  {
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithSlipLinks(universe,
                                                     { 0, 0 },
                                                     { 1, 1 },
                                                     { 21., 23. },
                                                     { 18., 18. },
                                                     { 0.0, 0.0 },
                                                     { 0.4, 0.6 },
                                                     { 0.4, 0.6 });
    CHECK_THROWS(
      forceBalancer.setSpringContourLengths(Eigen::VectorXd::Constant(2, 30.)));
    forceBalancer.setSpringContourLengths(Eigen::VectorXd::Constant(3, 30.));

    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    long int edgeBetweenSlipLinks = -1;
    for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
      if (net.linkIsSliplink[net.springPartIndexA[i]] &&
          net.linkIsSliplink[net.springPartIndexB[i]]) {
        edgeBetweenSlipLinks = i;
      }
    }
    CHECK(edgeBetweenSlipLinks >= 0);
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    outputNetwork(net, partitions);
    CHECK_FALSE(forceBalancer.swapSlipLinkReversibly(
      net, partitions, edgeBetweenSlipLinks, 1.));
    CHECK(partitions.isApprox(forceBalancer.getSpringPartitions()));
  }
}

TEST_CASE(
  "MEHP Force Balance 2 MC swap accept and reject work with cross-links",
  "[analysis][MEHPForceBalance2]")
{
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

  SECTION("No slip-links, state check")
  {
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2);
    forceBalancer.setSpringContourLengths(Eigen::VectorXd::Constant(4, 30.));
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    CHECK((net.nrOfPartialSprings == net.nrOfSprings));
    CHECK((net.nrOfSprings == 4));
  }

  SECTION("MC condition accepts as requested")
  {
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithSlipLinks(universe,
                                                     { 0 },
                                                     { 1 },
                                                     { 10. },
                                                     { 8.1 },
                                                     { 0. },
                                                     { 0.5 },
                                                     { 0.966667 });
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    CHECK((net.nrOfPartialSprings != net.nrOfSprings));
    CHECK((net.nrOfSprings == 4));
    CHECK((net.nrOfPartialSprings == 6));
    CHECK((net.nrOfCrosslinkSwapsEndured[0] == 0));
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    CHECK(net.linkIndicesOfSprings[1].size() == 3);
    outputNetwork(net, partitions);

    CHECK_FALSE(
      forceBalancer.swapSlipLinkReversibly(net, partitions, 5, 1., 0));
    CHECK(forceBalancer.swapSlipLinkReversibly(net, partitions, 5, 1., 5));
    net = forceBalancer.getNetwork();
    CHECK_FALSE(net.isUpToDate);
    forceBalancer.synchronise();
    net = forceBalancer.getNetwork();
    CHECK(net.isUpToDate);
    // check that the connectivity around the cross-link changed
    outputNetwork(net, partitions);
    CHECK(net.linkIndicesOfSprings[1].size() == 2);
    CHECK((net.nrOfCrosslinkSwapsEndured[0] == 1));
  }

  SECTION("MC condition rejects as requested")
  {
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithSlipLinks(universe,
                                                     { 0 },
                                                     { 2 },
                                                     { 10. },
                                                     { 4.9 },
                                                     { 0. },
                                                     { 0.5 },
                                                     { 0.966667 });
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    CHECK(net.nrOfPartialSprings != net.nrOfSprings);
    CHECK(net.nrOfSprings == 4);
    CHECK(net.nrOfPartialSprings == 6);
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    outputNetwork(net, partitions);
    CHECK_FALSE(forceBalancer.swapSlipLinkReversibly(net, partitions, 4, 1.));
    CHECK(net.nrOfCrosslinkSwapsEndured[0] == 0);
    // ideally, we could check that the values stayed the same, but they
    // probably did not as of the current implementation. at least we can check
    // that the configuration stayed the same
  }
}

TEST_CASE("MEHP Force Balance 2 handles slip-links on primary loops",
          "[analysis][MEHPForceBalance2]")
{
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
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2);
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd originalCoords = net.coordinates;
    // check unentangled primary loops
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    CHECK_NOTHROW(forceBalancer.displaceToMeanPosition(net, partitions, 0));
    Eigen::VectorXd displacements = net.coordinates - originalCoords;
    CHECK(displacements[0] == Catch::Approx(7.687).epsilon(1e-5));
    CHECK(displacements[1] == Catch::Approx(2.03926).epsilon(1e-5));
    CHECK(displacements[2] == Catch::Approx(5.88634).epsilon(1e-5));
  }

  SECTION("Entangled primary loop")
  {
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithSlipLinks(
        universe,
        { 0, 2 },
        { 1, 1 },
        { 1.3263401618628183, 42.04664022316877 },
        { 6.670300217824844, 7.18272624553976 },
        { 41.85951429390015, 0.8578704100544575 },
        { 0.5333333333333333, 0.9310344827586207 },
        { 0.5, 0.5 });
    forceBalancer.setSpringContourLengths(
      Eigen::VectorXd::Constant(forceBalancer.getNetwork().nrOfSprings, 15.));
    // entangle primary loop and check again
    // outputNetwork(net, Eigen::VectorXd::Zero(net.nrOfLinks * 3));
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    outputNetwork(net, forceBalancer.getSpringPartitions());
    Eigen::VectorXd displacements = Eigen::VectorXd::Zero(net.nrOfLinks * 3);
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    Eigen::VectorXd originalCoords = forceBalancer.getNetwork().coordinates;
    CHECK_NOTHROW(forceBalancer.displaceToMeanPosition(net, partitions, 0, 1.));
    displacements = net.coordinates - originalCoords;
    CHECK(displacements[0] == Catch::Approx(0.187321).epsilon(1e-5));
    CHECK(displacements[1] == Catch::Approx(-0.447774).epsilon(1e-5));
    CHECK(displacements[2] == Catch::Approx(-0.925295).epsilon(1e-5));
  }
}

TEST_CASE("MEHP Force Balance 2 handles slip-link convergence correctly",
          "[analysis][MEHPForceBalance2]")
{
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

  pcm::MEHPForceBalance2 forceBalancer =
    pcm::MEHPForceBalance2::constructWithSlipLinks(
      universe,
      { 1, 0 },
      { 1, 1 },
      { 12.650493316819828, 13.197029579176265 },
      { 2.8706102036538566, 3.4016980009809297 },
      { 8.475863644409664, 5.284899588057222 },
      { 1. - 0.13793103448275862, 1. - 0.3103448275862069 },
      { 1. - 0.13793103448275862, 1. - 0.7931034482758621 });

  // do update step
  Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
  outputNetwork(forceBalancer.getNetwork(), springPartitions);
  CHECK(true);

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

TEST_CASE("MEHP Force Balance 2 runs", "[analysis][MEHPForceBalance2][long]")
{
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
      //   pcm::MEHPForceBalance2 forceBalancer3 =
      //     pcm::MEHPForceBalance2(universe2, 2);
      //   meter.measure([&forceBalancer3] {
      //     forceBalancer3.runForceRelaxation("LD_MMA");
      //     return forceBalancer3.getNrOfIterations();
      //   });
      // };
      // BENCHMARK_ADVANCED("MEHP LD_LBFGS " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceBalance2 forceBalancer3 =
      //     pcm::MEHPForceBalance2(universe2, 2);
      //   meter.measure([&forceBalancer3] {
      //     forceBalancer3.runForceRelaxation("LD_LBFGS");
      //     return forceBalancer3.getNrOfIterations();
      //   });
      // };

      double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      pcm::MEHPForceBalance2 forceBalancer2 =
        pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe2, 2);

      // SECTION("Stress tensor computations are equivalent")
      // {
      Eigen::Matrix3d stressTensor1 = forceBalancer2.getStressTensor();
      Eigen::Matrix3d stressTensor2 = forceBalancer2.getStressTensorLinkBased();
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          CHECK(stressTensor1(i, j) == Catch::Approx(stressTensor2(i, j)));
        }
      }
      // }

      // SECTION("Displacement computations are equivalent")
      // {
      //   pcm::ForceBalanceNetwork net = forceBalancer2.getNetwork();
      //   Eigen::VectorXd springPartitions0 =
      //     Eigen::VectorXd::Ones(net.nrOfPartialSprings);
      //   Eigen::VectorXd oneOverSpringPartitions =
      //     forceBalancer2.assembleOneOverSpringPartition(net,
      //     springPartitions0);
      //   CHECK((oneOverSpringPartitions.array() < net.L[0]).all());
      //   Eigen::VectorXd displacements0 =
      //     Eigen::VectorXd::Zero(3 * net.nrOfLinks);
      //   std::vector<Eigen::ArrayXi> vertexSets;
      //   std::vector<Eigen::ArrayXi> springSets;
      //   std::tie(vertexSets, springSets) =
      //     forceBalancer2.getHeuristicallyIndependentCoordinateSets(net);

      //   // SECTION(
      //   //   "HeuristicallyIndependent coordiante sets are unique and
      //   complete")
      //   // {
      //   pcm::ArrayXb vertexSetTest =
      //     pcm::ArrayXb::Constant(3 * net.nrOfLinks, false);
      //   for (int i = 0; i < vertexSets.size(); ++i) {
      //     for (int j = 0; j < vertexSets[i].size(); ++j) {
      //       CHECK_FALSE(vertexSetTest[vertexSets[i][j]]);
      //       vertexSetTest[vertexSets[i][j]] = true;
      //     }
      //   }
      //   for (int i = 0; i < vertexSetTest.size(); ++i) {
      //     CHECK(vertexSetTest[i] == true);
      //   }
      //   // }

      //   for (Eigen::ArrayXi vertexSet : vertexSets) {
      //     forceBalancer2.displaceLinksToMeanPosition(
      //       net, oneOverSpringPartitions, vertexSet, 1.0);
      //     for (size_t i = 0; i < vertexSet.size(); ++i) {
      //       if (i % 3 == 0) {
      //         // NOTE: it is CHECKd, that we process the iterative updates
      //         // also in the order of the vertex sets, as otherwise, the
      //         results
      //         // will not be identical
      //         forceBalancer2.displaceToMeanPosition(
      //           net, springPartitions0, vertexSet[i] / 3);
      //       }
      //     }
      //   }

      //   CHECK(forceBalancer2.validateNetwork(net));
      // }

      SECTION("Actual balance results in correct phantom results")
      {
        std::cout << "Doing phantom force balance" << std::endl;
        pcm::MEHPForceRelaxation forceRelaxer = pcm::MEHPForceRelaxation(
          universe2, 2, false, nullptr, 1.0, false, false);
        CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
        CHECK(forceBalancer2.getNrOfIterations() == 0);
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(97.383096 * 97.383096 * 97.383096));
        CHECK(forceBalancer2.getNrOfSprings() == 10000);
        // remove the inactive ones
        CHECK_NOTHROW(forceBalancer2.removeFreeChains());
        CHECK_NOTHROW(forceBalancer2.synchronise());
        CHECK(forceBalancer2.getNrOfSprings() ==
              9848); // TODO: check that this is reasonable?
        CHECK(forceRelaxer.getNrOfSprings() == 9859);
        // initial system values
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(forceRelaxer.getPressure()));
        CHECK(forceBalancer2.getNrOfSprings() == forceRelaxer.getNrOfSprings());
        CHECK(
          forceRelaxer.getNetwork().meanSpringContourLength ==
          Catch::Approx(forceBalancer2.getNetwork().meanSpringContourLength));
        CHECK(forceBalancer2.getPressure() == Catch::Approx(0.0061105865));
        CHECK_NOTHROW(forceBalancer2.runForceRelaxation());
        CHECK_NOTHROW(forceBalancer2.validateNetwork());
        std::cout << "Removing subfunctional vertices" << std::endl;
        CHECK_NOTHROW(forceBalancer2.removeSubfunctionalVertices());
        CHECK_NOTHROW(forceBalancer2.synchronise());
        CHECK(forceBalancer2.getNrOfSprings() == 9848);
        CHECK_NOTHROW(forceBalancer2.removeInactiveCrosslinks());
        std::cout << "Removing inactive cross-links" << std::endl;
        CHECK_NOTHROW(forceBalancer2.synchronise());
        CHECK(forceBalancer2.getNrOfSprings() == 6693);
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
        double nb2Correction =
          (forceBalancer2.getDefaultR0Square() / (expectedNb2));
        double gammaCorrectionFactor = nrOfChainCorrection * nb2Correction;
        CHECK(
          forceBalancer2.getGammaFactor() * nrOfChainCorrection *
            forceBalancer2.getDefaultR0Square() ==
          Catch::Approx(42.6132)); // as from conversion-less Mathematica script
        CHECK(forceBalancer2.getGammaFactor() * gammaCorrectionFactor * kb * T *
                nu ==
              Catch::Approx(61308.3)); // ANT shear modulus, Pa
        CHECK(forceBalancer2.getGammaFactor() * gammaCorrectionFactor ==
              Catch::Approx(0.319446)); // "correct" gamma factor
        CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
        CHECK_NOTHROW(forceBalancer2.validateNetwork());
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
        CHECK(forceBalancer2.getPressure() == Catch::Approx(0.0061105865));
        // add entanglements
        // TODO: these are random values, as are the results... :P
        size_t nrOfSprings = forceRelaxer.getNetwork().nrOfSprings;
        pcm::MEHPForceBalance2 forceBalancer3 =
          pcm::MEHPForceBalance2::constructWithSlipLinks(
            universe2,
            { 10, 100, 50, 12, 76, 80, nrOfSprings - 1 },
            { 99, 101, 13, 7, 5, 19, nrOfSprings - 7 },
            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
            { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 },
            { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 });
        CHECK_NOTHROW(forceBalancer3.runForceRelaxation());
        // TODO: replace this value thereafter
        CHECK(forceBalancer3.getPressure() == Catch::Approx(0.001955022));
      }
    } else {
      std::cout << "Skipping large file PDMS MEHP run" << std::endl;
      CHECK(true);
    }
  }

  // SECTION("MEHP Force Balance 2 2D case")
  // {
  //   CHECK(std::filesystem::exists(suspectedPath));
  //   universeSeq.initializeFromDataSequence(
  //     { { suspectedPath +
  //         "structure/equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
  //         "333_2d_t_7500001.structure.out" } });
  //   CHECK(universeSeq.getLength() == 1);
  //   pe::Universe universe = universeSeq.atIndex(0);
  //   pcm::MEHPForceBalance2 forceBalancer =
  //     pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2, true);
  //   CHECK(forceBalancer.getExitReason() == pcm::ExitReason::UNSET);
  //   CHECK(forceBalancer.getNrOfIterations() == 0);
  //   CHECK(forceBalancer.getVolume() == Catch::Approx(universe.getVolume()));
  //   forceBalancer.runForceRelaxation(pcm::BalanceRunMode::ITERATIVE, 1.0, 5);
  //   CHECK(forceBalancer.getNrOfNodes() != universe.getNrOfAtoms());
  //   CHECK(forceBalancer.getNrOfIterations() <= 5);
  //   CHECK(forceBalancer.getNrOfIterations() >= 1);
  //   CHECK(universe.getAtomsOfType(2).size() == 7200);
  //   CHECK(forceBalancer.getExitReason() == pcm::ExitReason::MAX_STEPS);
  //   CHECK_NOTHROW(forceBalancer.validateNetwork());
  //   // run again, this time fully
  //   pcm::MEHPForceBalance2 forceBalancer2 =
  //     pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2, true);
  //   forceBalancer2.runForceRelaxation(
  //     pcm::BalanceRunMode::ITERATIVE, 1.0, 10000, 1e-10, 100);
  //   CHECK(forceBalancer2.getNrOfIterations() > 5);
  //   CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  //   CHECK(forceBalancer2.getGammaFactor(25, forceBalancer2.getNrOfSprings())
  //   ==
  //         Catch::Approx(1. / 3.).epsilon(0.001));
  //   auto stressTensor = forceBalancer2.getStressTensor();
  //   CHECK(forceBalancer2.getPressure() ==
  //         Catch::Approx(
  //           (stressTensor(0, 0) + stressTensor(1, 1) + stressTensor(2, 2))
  //           / 3.) .epsilon(0.02));
  //   CHECK_NOTHROW(forceBalancer2.validateNetwork());
  //   // TODO: find better, more accurate tests here
  //   CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
  //   CHECK(forceBalancer2.getNrOfActiveSprings() > 1);
  //   CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
  //   CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
  //         forceBalancer2.getNrOfNodes());

  //   pe::Universe universe3 = forceBalancer2.getCrosslinkerVerse();
  //   CHECK(universe3.getNrOfAtoms() == forceBalancer2.getNrOfNodes());
  //   CHECK(universe3.getNrOfBonds() == forceBalancer2.getNrOfSprings());
  //   CHECK(universe3.getAtomsOfType(2).size() == universe3.getNrOfAtoms());

  //   // try out different algorithms
  //   std::vector<std::string> algorithms = { "LD_MMA",
  //                                           // "LD_TNEWTON_PRECOND_RESTART",
  //                                           // "GD_STOGO",
  //                                           "LD_SLSQP",
  //                                           "GN_DIRECT" };

  //   // for (std::string algorithm : algorithms) {
  //   //   pcm::MEHPForceBalance2 forceBalancerN =
  //   //     pcm::MEHPForceBalance2(universe, 2);
  //   //   std::cout << "Testing algorithm " << algorithm << std::endl;
  //   //   auto start = std::chrono::high_resolution_clock::now();
  //   //   forceBalancerN.runForceRelaxation(true, 15, algorithm.c_str(),
  //   10000);
  //   //   auto stop = std::chrono::high_resolution_clock::now();
  //   //   CHECK(forceBalancerN.getGammaEq() ==
  //   //         Catch::Approx(forceBalancer2.getGammaEq()));
  //   //   CHECK(forceBalancerN.getFinalPressure() ==
  //   //         Catch::Approx(forceBalancer2.getFinalPressure()));
  //   //   auto duration = duration_cast<std::chrono::microseconds>(stop -
  //   start);
  //   //   std::cout << "Took: " << duration.count() << std::endl;
  //   // }
  // }
}

TEST_CASE(
  "MEHP Force Balance 2 can randomly add slip-links ignoring cross-links",
  "[analysis][MEHPForceBalance2]")
{
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
    // TODO: using a seed does not seem to work properly?!?
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithRandomSlipLinks(
        universe, 250, 2.0, 100, 2.0, 12, 2, false, 1.0);
    size_t nrOfAddedLinks = forceBalancer.getNumExtraAtoms();
    CHECK(nrOfAddedLinks >= 50);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

    size_t numRemoved = forceBalancer.removeTwofunctionalLinks();
    CHECK_NOTHROW(forceBalancer.synchronise());
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numRemoved > 0);
    // numRemoved = forceBalancer.removeInactiveCrosslinks(net, displacements,
    // partitions, 1e-20); CHECK(numRemoved == 204); // TODO: analyze these
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    // remove all springs...
    numRemoved = forceBalancer.removeInactiveCrosslinks(1e5);
    // resp., remove all cross-links.
    CHECK(numRemoved <= forceBalancer.getNumAtoms());
    CHECK(forceBalancer.getNumBonds() == 0);
  }
}

TEST_CASE("MEHP Force Balance can run with swapping slip-links",
          "[analysis][MEHPForceBalance2][long]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  std::string inputFile =
    // suspectedPath + "structure/network_p_1_100_a_38_50_xlinks.structure.out";
    suspectedPath + "structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithRandomSlipLinks(
        universe, 250, 2.0, 100, 2.0, -1, 2, false, 1.0);
    size_t nrOfAddedLinks = forceBalancer.getNumExtraAtoms();
    CHECK(nrOfAddedLinks >= 100);
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
    pe::Box oldBox = universe.getBox();
    pe::Box newBox =
      pe::Box(4 * oldBox.getLx(), 0.5 * oldBox.getLy(), 0.5 * oldBox.getLz());
    forceBalancer.deformTo(newBox);
    CHECK_NOTHROW(forceBalancer.runForceRelaxation(
      pcm::BalanceRunMode::ITERATIVE,
      1.0,
      1000,
      1e-9,
      -1.0,
      pcm::StructureSimplificationMode::ALL_TIM,
      -1.0,
      50,
      false,
      pcm::LinkSwappingMode::ALL_MC));
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK_NOTHROW(forceBalancer.runForceRelaxation(
      pcm::BalanceRunMode::ITERATIVE,
      1.0,
      1000,
      1e-9,
      -1.0,
      pcm::StructureSimplificationMode::ALL_TIM,
      -1.0,
      50,
      false,
      pcm::LinkSwappingMode::SLIPLINKS_ONLY));
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK_NOTHROW(forceBalancer.runForceRelaxation(
      pcm::BalanceRunMode::ITERATIVE,
      1.0,
      1000,
      1e-9,
      -1.0,
      pcm::StructureSimplificationMode::ALL_TIM,
      -1.0,
      50,
      false,
      pcm::LinkSwappingMode::ALL));
  }
}

TEST_CASE("MEHP Force Balance 2 can randomly add and remove slip-links",
          "[analysis][MEHPForceBalance2]")
{
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
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithRandomSlipLinks(
        universe, 1000, 2.0, 100, 2.0, 23, 2, false, 1.0);
    size_t nrOfAddedLinks = forceBalancer.getNumExtraAtoms();
    CHECK(nrOfAddedLinks >= 100);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    size_t numRemoved = forceBalancer.removeTwofunctionalLinks();
    CHECK(numRemoved > 0);
    CHECK_NOTHROW(forceBalancer.synchronise());
    CHECK_NOTHROW(forceBalancer.validateNetwork());

    // run a while to get inactive links
    forceBalancer.runForceRelaxation(pcm::BalanceRunMode::ITERATIVE, 1.0, 100);
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    partitions = forceBalancer.getSpringPartitions();
    // due to the randomness, it _could_ be one day that actually all strands
    // are active. unlikely, but I can imagine it to be possible.
    size_t numInactiveRemoved = forceBalancer.removeInactiveCrosslinks(0.1);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(numInactiveRemoved > 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork(net, partitions));

    ////////////////////////////////////////////////////////////////
    forceBalancer = pcm::MEHPForceBalance2::constructWithRandomSlipLinks(
      universe, 1000, 2.0, 100, 2.0, 23, 2, true, 1.0);
    nrOfAddedLinks = forceBalancer.getNumExtraAtoms();
    CHECK(nrOfAddedLinks >= 100);
    // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
    // check that all f = 2 have already been removed
    // they have not, since more f = 2 are produced by
    // numInactiveRemoved = forceBalancer.removeTwofunctionalLinks(
    //       net, displacements, partitions);
    // CHECK(numInactiveRemoved == 0);

    // run a while to get inactive links
    forceBalancer.runForceRelaxation(pcm::BalanceRunMode::ITERATIVE, 1.0, 100);
    net = forceBalancer.getNetwork();
    // due to the randomness, it _could_ be one day that actually all strands
    // are active. unlikely, but I can imagine it to be possible.
    numInactiveRemoved = forceBalancer.removeInactiveCrosslinks(0.1);
    CHECK(numInactiveRemoved > 0);
    numInactiveRemoved = forceBalancer.removeTwofunctionalLinks();
    CHECK(numInactiveRemoved > 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }
}

TEST_CASE("MEHP Force Balance handles slip-links",
          "[analysis][MEHPForceBalance2]")
{
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
  universe.setBox(pe::Box(-50.0, 50.0, -50.0, 50.0, -50.0, 50.0));
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

    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithSlipLinks(universe,
                                                     { 0 },
                                                     { 2 },
                                                     { 0.0 },
                                                     { 0.0 },
                                                     { 0.0 },
                                                     { 0.62 },
                                                     { 0.43 },
                                                     2,
                                                     false);

    forceBalancer.setSpringContourLengths(
      Eigen::VectorXd::Constant(forceBalancer.getNetwork().nrOfSprings, 15.));

    Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
    CHECK(springPartitions[3] == Catch::Approx(0.62));
    CHECK(springPartitions[4] == Catch::Approx(1. - 0.62));
    CHECK(springPartitions[5] == Catch::Approx(0.43));
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
    // add a slip-link between the strands 2-4 & 1-3
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2::constructWithSlipLinks(
        universe, { 4 }, { 5 }, { 4.2 }, { 3.9 }, { 1.2 }, { 0.5 }, { 0.5 });
    CHECK(forceBalancer.getNumExtraAtoms() == 1);
    CHECK(forceBalancer.getNumAtoms() == 4);
    size_t displacedId = 4;

    // ...at the wrong coordinates, to see it converge to the center
    Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    Eigen::VectorXd originalCoords = net.coordinates;

    outputNetwork(net, forceBalancer.getSpringPartitions());
    for (int i = 0; i < 5; ++i) {
      forceBalancer.displaceToMeanPosition(net, springPartitions, displacedId);
      forceBalancer.updateSpringPartition(net, springPartitions, displacedId);
      for (int i = 0; i < net.nrOfPartialSprings; i++) {
        std::cout << net.springPartIndexA[i] << ", " << net.springPartIndexB[i]
                  << ": ";
        std::cout << springPartitions[i] << std::endl;
      }
      std::cout << std::endl;
    }
    Eigen::VectorXd displacements = net.coordinates - originalCoords;
    CHECK(springPartitions[springPartitions.size() - 1] == Catch::Approx(0.5));
    CHECK(displacements[3 * displacedId] == Catch::Approx(-4.2));
    CHECK(displacements[3 * displacedId + 1] == Catch::Approx(-3.9));
    CHECK(displacements[3 * displacedId + 2] == Catch::Approx(0.8));
  }

  SECTION("aye")
  {
    // now, construct the force balancer
    pcm::MEHPForceBalance2 forceBalancer2 =
      pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2, false);
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

    SECTION("displaceLinksToMeanPosition works")
    {
      // test the "displaceLinksToMeanPosition"
      // see also: "Eigen behaves as expected"
      Eigen::VectorXd coords0 = net0.coordinates;
      forceBalancer2.displaceLinksToMeanPosition(net0, springPartitions0, 0.);
      Eigen::VectorXd displacements0 = net0.coordinates - coords0;
      for (int i = 0; i < 3 * net0.nrOfLinks; ++i) {
        CHECK(displacements0[i] + 1e-5 == Catch::Approx(0.0 + 1e-5));
      }
      // repeat the displacement to reach the point where
      // all beads are at 0.0
      for (size_t it = 0; it < 20; ++it) {
        forceBalancer2.displaceLinksToMeanPosition(
          net0, springPartitions0, 0.5);
      }
      displacements0 = net0.coordinates - coords0;
      for (int i = 0; i < net0.nrOfLinks; ++i) {
        CHECK(displacements0[3 * i] + 1e-5 ==
              Catch::Approx(-coords0[3 * i] + 1e-5));
        CHECK(displacements0[3 * i + 1] + 1e-5 ==
              Catch::Approx(-coords0[3 * i + 1] + 1e-5));
        CHECK(displacements0[3 * i + 2] + 1e-5 == Catch::Approx(0.0 + 1e-5));
      }
    }

    SECTION("Irrelevant slip-links are rationalised")
    { // try to add slip-links
      CHECK_THROWS(
        pcm::MEHPForceBalance2::constructWithSlipLinks(universe,
                                                       { { 1, 2, 3 } },
                                                       { { 1, 2, 3 } },
                                                       { { 1., 2., 3. } },
                                                       { { 1., 2., 3. } },
                                                       { { 1., 2. } },
                                                       { { 0.5, 0.5 } },
                                                       { { 0.5, 0.5 } }));
      CHECK_THROWS(
        pcm::MEHPForceBalance2::constructWithSlipLinks(universe,
                                                       { { 1, 2, 3 } },
                                                       { { 1, 2 } },
                                                       { { 1., 2., 3. } },
                                                       { { 1., 2., 3. } },
                                                       { { 1., 2., 3.0 } },
                                                       { { 0.5, 0.5 } },
                                                       { { 0.5, 0.5 } }));
      // and actually add slip-links
      forceBalancer2 =
        pcm::MEHPForceBalance2::constructWithSlipLinks(universe,
                                                       { { 0, 3 } },
                                                       { { 1, 3 } },
                                                       { { 4.2, 1.3 } },
                                                       { { 3.9, 1.2 } },
                                                       { { 1.3, 1.1 } },
                                                       { { 0.5, 0.5 } },
                                                       { { 0.5, 0.5 } });
      CHECK(forceBalancer2.getNrOfNodes() + 2 == forceBalancer2.getNrOfLinks());
      // and run with them.
      // Expect the slip-link of between two strands to converge to the central
      // atom Expect the slip-link of two strands to stay at 0.5, 0.5.
      Eigen::VectorXd springPartitions = forceBalancer2.getSpringPartitions();
      pcm::ForceBalanceNetwork net = forceBalancer2.getNetwork();
      for (int i = 0; i < net.nrOfPartialSprings; i++) {
        std::cout << net.springPartIndexA[i] << ", " << net.springPartIndexB[i]
                  << ": ";
        std::cout << springPartitions[i] << std::endl;
      }
      std::cout << std::endl;

      SECTION("Displacement example 1")
      {
        forceBalancer2.displaceToMeanPosition(net, springPartitions, 4);
        CHECK(net.coordinates[4 * 3] == Catch::Approx(2.5));
        CHECK(net.coordinates[4 * 3 + 1] == Catch::Approx(2.5));
        CHECK(net.coordinates[4 * 3 + 2] == Catch::Approx(2));
      }

      // reset
      SECTION("Displacement Example 2")
      {
        CHECK(net.springIndicesOfLinks.size() == net.nrOfLinks);
        for (int j = 0; j < 5; ++j) {
          forceBalancer2.displaceToMeanPosition(net, springPartitions, 5, -1.);
          forceBalancer2.updateSpringPartition(net, springPartitions, 5, -1.);
          for (int i = 0; i < net.nrOfPartialSprings; i++) {
            std::cout << net.springPartIndexA[i] << ", "
                      << net.springPartIndexB[i] << ": ";
            std::cout << springPartitions[i] << std::endl;
          }
          std::cout << std::endl;
        }
        for (int i = 0; i < 125; ++i) {
          // do some random 125 steps with these two slip-links
          // NOTE: difficulty: finding out which node and spring it is actually
          // after the removal of strand atoms
          forceBalancer2.displaceToMeanPosition(net, springPartitions, 4, -1.);
          // for (int dir = 0; dir < 3; ++dir) {
          //   std::cout << (net.coordinates[3 * 4 + dir] + displacements[3 * 4
          //   + dir])
          //             << ", ";
          // }
          forceBalancer2.updateSpringPartition(net, springPartitions, 4, -1.);
          // std::cout << std::endl;
          // std::cout << springPartitions[0][0] << ", " <<
          // springPartitions[1][0]
          //           << std::endl;
        }
        // assert expectations are met.
        // NOTE: difficulty: finding out which spring idx it actually is
        outputNetwork(net, forceBalancer2.getSpringPartitions());
        for (int i = 0; i < net.nrOfPartialSprings; i++) {
          std::cout << net.springPartIndexA[i] << ", "
                    << net.springPartIndexB[i] << ": ";
          std::cout << springPartitions[i] << std::endl;
        }
        std::cout << std::endl;
        CHECK(springPartitions[1] == Catch::Approx(1.0).epsilon(1e-6));
        CHECK(springPartitions[6] == Catch::Approx(0.0).epsilon(1e-6));
        CHECK(springPartitions[0] == Catch::Approx(1.0).epsilon(1e-6));
        CHECK(springPartitions[5] == Catch::Approx(1.0).epsilon(1e-6));
        // CHECK(springPartitions[8] == Catch::Approx(1.0).margin(1e-6)); // 5-3
        // CHECK(springPartitions[7] + 1e-5 == Catch::Approx(0.0 +
        // 1e-5).margin(1e-6)); // 5-5 CHECK(springPartitions[3] ==
        // Catch::Approx(1.0).margin(1e-6)); // 5-0
      }
    }
  }
}

TEST_CASE(
  "MEHP Force Balance runs with non-network",
  "[analysis][MEHPForceBalance2][NonGaussianSpringForceEvaluator][long]")
{
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
      pcm::MEHPForceBalance2 forceBalancer2 =
        pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe2, 2, false);
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      CHECK_NOTHROW(forceBalancer2.runForceRelaxation());
      CHECK(forceBalancer2.getNrOfIterations() > 1);
    }
  }
}
TEST_CASE("MEHP Force Balance 2 Free chains collapse",
          "[analysis][MEHPForceBalance2][NonGaussianSpringForceEvaluator]["
          "SimpleSpringMEHPForceEvaluator]")
{
  // generate a set of connected chains that must collapse
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
  pcm::MEHPForceBalance2 forceBalancer =
    pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2, false);
  CHECK(forceBalancer.getNrOfSprings() ==
        forceBalancer.getNrOfPartialSprings());
  CHECK(forceBalancer.getNrOfSprings() == nrOfBeads / nrOfBeadsPerChain);
  CHECK_NOTHROW(forceBalancer.runForceRelaxation(
    pcm::BalanceRunMode::ITERATIVE, 1.0, 50000, 1e-18));
  CHECK(forceBalancer.getNrOfIterations() > 0);
  CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  CHECK(forceBalancer.getNrOfActiveSprings() == 0);
  // CHECK(forceBalancer.getAverageSpringLength() ==
  //       Catch::Approx(0.0));
  CHECK(forceBalancer.getAverageSpringLength() >= 0.0);
  CHECK(forceBalancer.getAverageSpringLength() <= 3e-6);
  CHECK_NOTHROW(forceBalancer.validateNetwork());

  std::cout << "Getting cross-linker verse" << std::endl;
  pe::Universe resultingUniverse = forceBalancer.getCrosslinkerVerse();
  // auto distances = resultingUniverse.computeBondLengths();
  // for (auto i : distances) {
  //   std::cout << i << std::endl;
  // }
  // auto residuals = forceBalancer.getResiduals();
  // for (auto i : residuals) {
  //   std::cout << i << " ";
  // }
  // std::cout << std::endl;
  // std::cout << forceBalancer.getNrOfIterations() << ", "
  //           << forceBalancer.getForce() << ", "
  //           << forceBalancer.getResidualNorm() << std::endl;

  Eigen::Matrix3d stressTensorSimpleSpring = forceBalancer.getStressTensor();
  Eigen::Matrix3d stressTensorLangevin = forceBalancer.getStressTensor();
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      CHECK(stressTensorLangevin(i, j) + 1e-5 ==
            Catch::Approx(stressTensorSimpleSpring(i, j) + 1e-5));
    }
  }

  SECTION("Total Removal Works")
  {
    size_t verticesToRemove = forceBalancer.getNrOfNodes();
    CHECK(forceBalancer.removeSubfunctionalVertices() == verticesToRemove);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }

  SECTION("Partial Removal Works")
  {
    size_t verticesToRemove = forceBalancer.getNrOfNodes() - 2;
    CHECK(forceBalancer.removeTwofunctionalLinks() == verticesToRemove);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }
}

TEST_CASE("MEHP Force Balance 2 does not collapse",
          "[analysis][MEHPForceBalance2]")
{
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  /**
   * @brief A grid of two rows, each one bead between the two cross-links
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

  SECTION("Running conventional MEHP")
  {
    pcm::MEHPForceBalance2 forceRelaxerConventional =
      pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2, true);
    forceRelaxerConventional.configAssumeBoxLargeEnough(true);
    REQUIRE_NOTHROW(forceRelaxerConventional.runForceRelaxation());
    REQUIRE(forceRelaxerConventional.getNrOfIterations() > 0);
    CHECK(forceRelaxerConventional.getExitReason() ==
          pcm::ExitReason::X_TOLERANCE);
    CHECK(forceRelaxerConventional.getNrOfActiveSprings() ==
          forceRelaxerConventional.getNrOfSprings());
    // compare to what we expect
    CHECK(forceRelaxerConventional.getNrOfActiveSprings() == 8);
    CHECK(forceRelaxerConventional.getNrOfActiveNodes() == 4);
    // CHECK(forceRelaxerConventional.getAverageSpringLength() ==
    // Catch::Approx(5.)); CHECK_THAT(forceRelaxerConventional.getGammaFactor(),
    //            Catch::Matchers::WithinAbs(1.0, 1e-3));
  }
  SECTION("Running new MEHP")
  {
    pcm::MEHPForceBalance2 forceRelaxerNew =
      pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2, true);
    forceRelaxerNew.configAssumeBoxLargeEnough(false);
    REQUIRE_NOTHROW(forceRelaxerNew.runForceRelaxation());
    REQUIRE(forceRelaxerNew.getNrOfIterations() > 0);
    CHECK(forceRelaxerNew.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceRelaxerNew.getNrOfActiveSprings() ==
          forceRelaxerNew.getNrOfSprings());

    // compare to what we expect
    CHECK(forceRelaxerNew.getNrOfActiveSprings() == 8);
    CHECK(forceRelaxerNew.getNrOfActiveNodes() == 4);
    CHECK(forceRelaxerNew.getAverageSpringLength() == Catch::Approx(5.0));
    CHECK_THAT(forceRelaxerNew.getGammaFactor(),
               Catch::Matchers::WithinAbs(1.0, 1e-2));
  }
};

TEST_CASE("MEHP Force Balance 2 Fully active chains are fully active",
          "[analysis][MEHPForceBalance2]")
{
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

      // validate that this is indeed an infinite structure
      std::vector<pe::Molecule> molecules =
        universe.getChainsWithCrosslinker(2);
      for (pe::Molecule mol : molecules) {
        CHECK(mol.getType() == pe::NETWORK_STRAND);
      }

      pcm::MEHPForceBalance2 forceBalancer =
        pcm::MEHPForceBalance2::constructWithoutSlipLinks(universe, 2, false);
      CHECK(forceBalancer.getNrOfActiveSprings() ==
            forceBalancer.getNrOfSprings());
      CHECK(forceBalancer.getNrOfSprings() == 16000);
      CHECK(forceBalancer.getNrOfPartialSprings() == 16000);
      CHECK_NOTHROW(forceBalancer.runForceRelaxation(
        pcm::BalanceRunMode::ITERATIVE, 1.0, 50000, 1e-12));
      CHECK(forceBalancer.getNrOfIterations() > 0);
      CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK(forceBalancer.getNrOfActiveSprings() ==
            forceBalancer.getNrOfSprings());
      CHECK(forceBalancer.removeSubfunctionalVertices() == 0);
      CHECK(forceBalancer.removeInactiveCrosslinks() == 0);
    } else {
      std::cout << "Diamond structure not found" << std::endl;
    }
  }
}
