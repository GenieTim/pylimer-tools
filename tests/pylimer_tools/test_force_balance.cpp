#include "../../src/pylimer_tools_cpp/calc/MEHPForceBalance.h"
#include "../../src/pylimer_tools_cpp/calc/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/calc/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
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
outputNetwork(pcm::ForceBalanceNetwork net,
              Eigen::VectorXd displacements,
              Eigen::VectorXd springPartitions)
{
  for (int i = 0; i < net.nrOfSprings; ++i) {
    std::cout << "Spring " << i << ", N: " << net.springsContourLength[i]
              << std::endl;
    for (int j = 0; j < net.linkIndicesOfSprings[i].size(); ++j) {
      std::cout << net.linkIndicesOfSprings[i][j] << ": ";
      for (int dir = 0; dir < 3; ++dir) {
        std::cout
          << (net.coordinates[3 * net.linkIndicesOfSprings[i][j] + dir] +
              displacements[3 * net.linkIndicesOfSprings[i][j] + dir])
          << ", ";
      }
      std::cout << std::endl;
      if (j < net.linkIndicesOfSprings[i].size() - 1) {
        std::cout << springPartitions[net.localToGlobalSpringIndex.at(i)[j]]
                  << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

TEST_CASE("Eigen behaves as required", "[analysis][MEHPForceBalance][Eigen]")
{
  SECTION("Summation works with same indices")
  {
    Eigen::VectorXd testVec = Eigen::VectorXd::Zero(10);
    Eigen::ArrayXi testIdx = Eigen::ArrayXi::Zero(5);
    testIdx << 0, 0, 5, 5, 1;
    testVec(testIdx) += Eigen::VectorXd::Ones(5);
    REQUIRE(testVec[5] == Catch::Approx(2.));
    REQUIRE(testVec[0] == Catch::Approx(2.));
    REQUIRE(testVec[1] == Catch::Approx(1.));
    REQUIRE(testVec[2] == 0.0);
  }

  SECTION("Casting bool to double results in 1.0/0.0")
  {
    auto gen = std::bind(std::uniform_int_distribution<>(0, 1),
                         std::default_random_engine());
    Eigen::Array<bool, 1, 100> boolArray;
    for (int i = 0; i < 100; i++) {
      bool b = gen();
      boolArray[i] = b;
    }
    Eigen::ArrayXd castedBoolArray = boolArray.cast<double>();
    for (int i = 0; i < 100; i++) {
      if (boolArray[i]) {
        CHECK(castedBoolArray[i] == 1.0);
      } else {
        CHECK(castedBoolArray[i] + 1e-5 == 1e-5);
      }
    }
  }
}

TEST_CASE("Force Balance Benchmarks", "[analysis][MEHPForceBalance]")
{
  return;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  std::string largeInputFile =
    suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
  if (std::filesystem::exists(largeInputFile)) {
    universeSeq.initializeFromDataSequence({ { largeInputFile } });
    pe::Universe universe2 = universeSeq.atIndex(0);
    //     BENCHMARK_ADVANCED("MEHP Balance Random 1.0 " + largeInputFile)
    //     (Catch::Benchmark::Chronometer meter)
    //     {
    //       pcm::MEHPForceBalance forceBalancer3 =
    //         pcm::MEHPForceBalance(universe2, 2);
    //       meter.measure([&forceBalancer3] {
    //         forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_RANDOM,
    //                                           1.0);
    //         return forceBalancer3.getNrOfIterations();
    //       });
    //     };
    //     BENCHMARK_ADVANCED("MEHP Balance Random 0.75 " + largeInputFile)
    //     (Catch::Benchmark::Chronometer meter)
    //     {
    //       pcm::MEHPForceBalance forceBalancer3 =
    //         pcm::MEHPForceBalance(universe2, 2);
    //       meter.measure([&forceBalancer3] {
    //         forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_RANDOM,
    //                                           0.75);
    //         return forceBalancer3.getNrOfIterations();
    //       });
    //     };
    //     BENCHMARK_ADVANCED("MEHP Balance Strand 1.0 " + largeInputFile)
    //     (Catch::Benchmark::Chronometer meter)
    //     {
    //       pcm::MEHPForceBalance forceBalancer3 =
    //         pcm::MEHPForceBalance(universe2, 2);
    //       meter.measure([&forceBalancer3] {
    //         forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_STRANDS,
    //                                           1.0);
    //         return forceBalancer3.getNrOfIterations();
    //       });
    //     };
    //     BENCHMARK_ADVANCED("MEHP Balance Strand 0.75 " + largeInputFile)
    //     (Catch::Benchmark::Chronometer meter)
    //     {
    //       pcm::MEHPForceBalance forceBalancer3 =
    //         pcm::MEHPForceBalance(universe2, 2);
    //       meter.measure([&forceBalancer3] {
    //         forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_STRANDS,
    //                                           0.75);
    //         return forceBalancer3.getNrOfIterations();
    //       });
    //     };
    //     BENCHMARK_ADVANCED("MEHP Balance All 0.5 " + largeInputFile)
    //     (Catch::Benchmark::Chronometer meter)
    //     {
    //       pcm::MEHPForceBalance forceBalancer3 =
    //         pcm::MEHPForceBalance(universe2, 2);
    //       meter.measure([&forceBalancer3] {
    //         forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_ALL,
    //                                           0.5);
    //         return forceBalancer3.getNrOfIterations();
    //       });
    //     };
    //     BENCHMARK_ADVANCED("MEHP Balance All 0.75 " + largeInputFile)
    //     (Catch::Benchmark::Chronometer meter)
    //     {
    //       pcm::MEHPForceBalance forceBalancer3 =
    //         pcm::MEHPForceBalance(universe2, 2);
    //       meter.measure([&forceBalancer3] {
    //         forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_ALL,
    //                                           0.75);
    //         return forceBalancer3.getNrOfIterations();
    //       });
    //     };
    // BENCHMARK_ADVANCED("MEHP Balance Heuristic 1.0 " + largeInputFile)
    // (Catch::Benchmark::Chronometer meter)
    // {
    //   pcm::MEHPForceBalance forceBalancer3 =
    //     pcm::MEHPForceBalance(universe2, 2);
    //   meter.measure([&forceBalancer3] {
    //     forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_HEURISTIC,
    //                                       1.0);
    //     return forceBalancer3.getNrOfIterations();
    //   });
    // };
    // BENCHMARK_ADVANCED("MEHP Balance Iterative 1.0 " + largeInputFile)
    // (Catch::Benchmark::Chronometer meter)
    // {
    //   pcm::MEHPForceBalance forceBalancer3 =
    //     pcm::MEHPForceBalance(universe2, 2);
    //   meter.measure([&forceBalancer3] {
    //     forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::ITERATIVE, 1.0);
    //     return forceBalancer3.getNrOfIterations();
    //   });
    // };

    BENCHMARK_ADVANCED("Detection of heuristic independent vertices " +
                       largeInputFile)
    (Catch::Benchmark::Chronometer meter)
    {
      pcm::MEHPForceBalance forceBalancer3 =
        pcm::MEHPForceBalance(universe2, 2);
      meter.measure([&forceBalancer3] {
        pcm::ForceBalanceNetwork net = forceBalancer3.getNetwork();
        return forceBalancer3.getHeuristicallyIndependentCoordinateSets(&net);
      });
    };
  }
}

TEST_CASE("MEHP Force Balance handles slip-links on primary loops",
          "[analysis][MEHPForceBalance]")
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
    pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();
    // check unentangled primary loops
    Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    REQUIRE_NOTHROW(
      forceBalancer.displaceToMeanPosition(&net, displacements, partitions, 0));
    CHECK(displacements[0] == Catch::Approx(7.687).epsilon(1e-5));
    CHECK(displacements[1] == Catch::Approx(2.03926).epsilon(1e-5));
    CHECK(displacements[2] == Catch::Approx(5.88634).epsilon(1e-5));
  }

  SECTION("Entangled primary loop")
  {
    pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
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
    outputNetwork(net,
                  Eigen::VectorXd::Zero(net.nrOfLinks * 3),
                  forceBalancer.getSpringPartitions());
    Eigen::VectorXd displacements = Eigen::VectorXd::Zero(net.nrOfLinks * 3);
    Eigen::VectorXd partitions = forceBalancer.getSpringPartitions();
    REQUIRE_NOTHROW(
      forceBalancer.displaceToMeanPosition(&net, displacements, partitions, 0));
    CHECK(displacements[0] == Catch::Approx(0.187321).epsilon(1e-5));
    CHECK(displacements[1] == Catch::Approx(-0.447774).epsilon(1e-5));
    CHECK(displacements[2] == Catch::Approx(-0.925295).epsilon(1e-5));
  }
}

TEST_CASE("MEHP Force Balance handles slip-link convergence correctly",
          "[analysis][MEHPForceBalance]")
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

  pcm::MEHPForceBalance forceBalancer = pcm::MEHPForceBalance(universe, 2);
  forceBalancer.addSlipLinks(
    { 1, 0 },
    { 1, 1 },
    { 12.650493316819828, 13.197029579176265 },
    { 2.8706102036538566, 3.4016980009809297 },
    { 8.475863644409664, 5.284899588057222 },
    { 1. - 0.13793103448275862, 1. - 0.3103448275862069 },
    { 1. - 0.13793103448275862, 1. - 0.7931034482758621 });

  // do update step
  Eigen::VectorXd displacements = Eigen::VectorXd::Zero(6 * 3);
  Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
  outputNetwork(forceBalancer.getNetwork(), displacements, springPartitions);
  /*auto results = */
  forceBalancer.inspectParametrisationOptimsationForLink(
    5,
    displacements,
    springPartitions,
    250,
    1e-10,
    0.0,
    0.0,
    100,
    1e10); // cannot use 1.0 for oneOver... without setting higher contour
           // length fraction
  CHECK(displacements[5 * 3] == Catch::Approx(-0.592091));
  CHECK(displacements[5 * 3 + 1] == Catch::Approx(0.441203));
  CHECK(displacements[5 * 3 + 2] == Catch::Approx(0.44577));
  outputNetwork(forceBalancer.getNetwork(), displacements, springPartitions);

  SECTION("Stress tensor computations are equivalent")
  {
    std::array<std::array<double, 3>, 3> stressTensor1 =
      forceBalancer.getStressTensor();
    std::array<std::array<double, 3>, 3> stressTensor2 =
      forceBalancer.getStressTensorLinkBased();
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        CHECK(stressTensor1[i][j] == Catch::Approx(stressTensor2[i][j]));
      }
    }
  }
}

TEST_CASE("MEHP Force Balance runs", "[analysis][MEHPForceBalance][long]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
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
      std::array<std::array<double, 3>, 3> stressTensor1 =
        forceBalancer2.getStressTensor();
      std::array<std::array<double, 3>, 3> stressTensor2 =
        forceBalancer2.getStressTensorLinkBased();
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          CHECK(stressTensor1[i][j] == Catch::Approx(stressTensor2[i][j]));
        }
      }
      // }

      SECTION("Displacement computations are equivalent")
      {
        pcm::ForceBalanceNetwork net = forceBalancer2.getNetwork();
        Eigen::VectorXd springPartitions0 =
          Eigen::VectorXd::Ones(net.nrOfPartialSprings);
        Eigen::VectorXd oneOverSpringPartitions =
          forceBalancer2.assembleOneOverSpringPartition(&net,
                                                        springPartitions0);
        CHECK((oneOverSpringPartitions.array() < net.L[0]).all());
        Eigen::VectorXd displacements0 =
          Eigen::VectorXd::Zero(3 * net.nrOfLinks);
        std::vector<Eigen::ArrayXi> vertexSets;
        std::vector<Eigen::ArrayXi> springSets;
        std::tie(vertexSets, springSets) =
          forceBalancer2.getHeuristicallyIndependentCoordinateSets(&net);

        // SECTION(
        //   "HeuristicallyIndependent coordiante sets are unique and complete")
        // {
        pcm::ArrayXb vertexSetTest =
          pcm::ArrayXb::Constant(3 * net.nrOfLinks, false);
        for (int i = 0; i < vertexSets.size(); ++i) {
          for (int j = 0; j < vertexSets[i].size(); ++j) {
            CHECK(vertexSetTest[vertexSets[i][j]] == false);
            vertexSetTest[vertexSets[i][j]] = true;
          }
        }
        for (int i = 0; i < vertexSetTest.size(); ++i) {
          CHECK(vertexSetTest[i] == true);
        }
        // }

        Eigen::VectorXd displacements1 =
          Eigen::VectorXd::Zero(3 * net.nrOfLinks);
        for (Eigen::ArrayXi vertexSet : vertexSets) {
          forceBalancer2.displaceLinksToMeanPosition(
            &net, displacements0, oneOverSpringPartitions, vertexSet, 1.0);
          for (size_t i = 0; i < vertexSet.size(); ++i) {
            if (i % 3 == 0) {
              // NOTE: it is required, that we process the iterative updates
              // also in the order of the vertex sets, as otherwise, the results
              // will not be identical
              forceBalancer2.displaceToMeanPosition(
                &net, displacements1, springPartitions0, vertexSet[i] / 3);
            }
            // check that the two displacement vectors are equal, here already
            if (std::isnan(displacements0[vertexSet[i]])) {
              // std::cout << "NaN at " << i << " in " << vertexSet[i]
              //           << std::endl;
              CHECK(std::isnan(displacements1[vertexSet[i]]));
            } else {
              CHECK(displacements0[vertexSet[i]] + 1e-5 ==
                    Catch::Approx(displacements1[vertexSet[i]] + 1e-5));
            }
          }
        }
        // check that the two displacement vectors coincide
        for (size_t i = 0; i < displacements0.size(); ++i) {
          if (std::isnan(displacements0[i])) {
            // std::cout << "NaN at " << i << " in " << vertexSet[i]
            //           << std::endl;
            CHECK(std::isnan(displacements1[i]));
          } else {
            CHECK(displacements0[i] + 1e-5 ==
                  Catch::Approx(displacements1[i] + 1e-5));
          }
        }

        CHECK(forceBalancer2.validateNetwork(&net));
      }

      SECTION("Actual balance results in correct phantom results")
      {
        std::cout << "Doing phantom force balance" << std::endl;
        pcm::MEHPForceRelaxation forceRelaxer =
          pcm::MEHPForceRelaxation(universe2, 2);
        // the strands are different -> cannot compare the distances anymore
        // CHECK((forceBalancer2.getCurrentSpringDistances() -
        // forceRelaxer.getCurrentSpringDistances()).isMuchSmallerThan(1e-12));
        REQUIRE(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
        REQUIRE(forceBalancer2.getNrOfIterations() == 0);
        REQUIRE(forceBalancer2.getVolume() ==
                Catch::Approx(universe2.getVolume()));
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(97.383096 * 97.383096 * 97.383096));
        // initial system values
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(forceRelaxer.getPressure()));
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(0.39911682390778536 / 79.));
        REQUIRE_NOTHROW(forceBalancer2.runForceRelaxation());
        CHECK_NOTHROW(forceBalancer2.validateNetwork());
        CHECK(forceBalancer2.getNrOfSprings() == 8142);
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
              Catch::Approx(0.000245543));
        double nu =
          nrOfChains / (forceBalancer2.getVolume() * sigmaToM * sigmaToM *
                        sigmaToM); // chain number density, m^-3
        CHECK(nu == Catch::Approx(4.63241e25));

        // final values
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(0.153806 / 79.)); // LJ Units [?]
        CHECK(forceBalancer2.getPressure() * conversionFactor /
                (sigmaToM * sigmaToM * sigmaToM) ==
              Catch::Approx(61308.3)); // shear modulus from the pressure, MPa
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
        REQUIRE(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
        REQUIRE(forceBalancer2.getNrOfIterations() == 0);
        REQUIRE(forceBalancer2.getVolume() ==
                Catch::Approx(universe2.getVolume()));
        CHECK(forceBalancer2.getVolume() ==
              Catch::Approx(97.383096 * 97.383096 * 97.383096));
        // initial system values
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(forceRelaxer.getPressure()));
        CHECK(forceBalancer2.getPressure() ==
              Catch::Approx(0.39911682390778536 / 79.));
        // add entanglements
        // TODO: these are random values, as are the results... :P
        size_t nrOfSprings = forceRelaxer.getNetwork().nrOfSprings;
        forceBalancer2.addSlipLinks(
          { 10, 100, 50, 12, 76, 80, nrOfSprings - 1 },
          { 99, 101, 13, 7, 5, 19, nrOfSprings - 7 },
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
        REQUIRE_NOTHROW(forceBalancer2.runForceRelaxation());
        // TODO: replace this value thereafter
        CHECK(forceBalancer2.getPressure() == Catch::Approx(0.0019534759));
      }
    } else {
      std::cout << "Skipping large file PDMS MEHP run" << std::endl;
      REQUIRE(true);
    }
  }

  SECTION("MEHP Force Balance 2D case")
  {
    REQUIRE(std::filesystem::exists(suspectedPath));
    universeSeq.initializeFromDataSequence(
      { { suspectedPath + "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
                          "333_2d_t_7500001.structure.out" } });
    REQUIRE(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, true);
    REQUIRE(forceBalancer.getExitReason() == pcm::ExitReason::UNSET);
    REQUIRE(forceBalancer.getNrOfIterations() == 0);
    REQUIRE(forceBalancer.getVolume() == Catch::Approx(universe.getVolume()));
    forceBalancer.runForceRelaxation(pcm::BalanceRunMode::ITERATIVE, 1.0, 5);
    REQUIRE(forceBalancer.getNrOfNodes() != universe.getNrOfAtoms());
    REQUIRE(forceBalancer.getNrOfIterations() <= 5);
    REQUIRE(forceBalancer.getNrOfIterations() >= 1);
    REQUIRE(universe.getAtomsOfType(2).size() == 7200);
    REQUIRE(forceBalancer.getExitReason() == pcm::ExitReason::MAX_STEPS);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    // run again, this time fully
    pcm::MEHPForceBalance forceBalancer2 =
      pcm::MEHPForceBalance(universe, 2, true);
    forceBalancer2.runForceRelaxation(
      pcm::BalanceRunMode::ITERATIVE, 1.0, 10000, 1e-10, 100, 1e-9);
    REQUIRE(forceBalancer2.getNrOfIterations() > 5);
    CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer2.getGammaFactor(25, forceBalancer2.getNrOfSprings()) ==
          Catch::Approx(1. / 3.).epsilon(0.001));
    auto stressTensor = forceBalancer2.getStressTensor();
    CHECK(forceBalancer2.getPressure() ==
          Catch::Approx(
            (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.)
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

TEST_CASE("MEHP Force Balance can randomly add slip-links",
          "[analysis][MEHPForceBalance]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  std::string largeInputFile =
    suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
  if (std::filesystem::exists(largeInputFile)) {
    REQUIRE(std::filesystem::exists(suspectedPath));
    universeSeq.initializeFromDataSequence(
      { { largeInputFile } });
    REQUIRE(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, true);
    REQUIRE_NOTHROW(forceBalancer.randomlyAddSliplinks(100));
  }
}

TEST_CASE("MEHP Force Balance handles slip-links",
          "[analysis][MEHPForceBalance]")
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
  universe.setBox(pe::Box(-15.0, 15.0, -15.0, 15.0, -15.0, 15.0));
  universe.addAtoms(9,
                    { { 1, 2, 3, 4, 5, 6, 7, 8, 9 } },     // id
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
      1, { 10 }, { 1 }, { 0. }, { 0. }, { 0. }, { 1 }, { 1 }, { 1 });
    universe.addBonds(2, { { 1, 10 } }, { { 10, 3 } });
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, false);
    // add a slip-link between the strands 2-4 & 1-3
    REQUIRE_NOTHROW(
      forceBalancer.addSlipLinks({ 4 }, { 5 }, { 4.2 }, { 3.9 }, { 1.2 }));
    // ...at the wrong coordinates, to see it converge to the center
    Eigen::VectorXd springPartitions = forceBalancer.getSpringPartitions();
    Eigen::VectorXd displacements =
      Eigen::VectorXd::Zero(forceBalancer.getNrOfLinks() * 3);
    pcm::ForceBalanceNetwork net = forceBalancer.getNetwork();

    outputNetwork(net, displacements, forceBalancer.getSpringPartitions());
    for (int i = 0; i < 5; ++i) {
      forceBalancer.displaceToMeanPosition(
        &net, displacements, springPartitions, 4);
      forceBalancer.updateSpringPartition(
        &net, displacements, springPartitions, 4);
      for (int i = 0; i < net.nrOfPartialSprings; i++) {
        std::cout << net.springPartIndexA[i] << ", " << net.springPartIndexB[i]
                  << ": ";
        std::cout << springPartitions[i] << std::endl;
      }
      std::cout << std::endl;
    }
    CHECK(springPartitions[springPartitions.size() - 1] == Catch::Approx(0.5));
    CHECK(displacements[3 * 4] == Catch::Approx(-4.2));
    CHECK(displacements[3 * 4 + 1] == Catch::Approx(-3.9));
    CHECK(displacements[3 * 4 + 2] == Catch::Approx(0.8));
  }

  SECTION("aye")
  {
    // now, construct the force balancer
    pcm::MEHPForceBalance forceBalancer2 =
      pcm::MEHPForceBalance(universe, 2, false);
    REQUIRE(forceBalancer2.getNrOfNodes() == forceBalancer2.getNrOfLinks());
    REQUIRE(forceBalancer2.getNrOfNodes() == 4);
    REQUIRE(forceBalancer2.getNrOfSprings() == 5);
    Eigen::VectorXd springPartitions0 = forceBalancer2.getSpringPartitions();
    pcm::ForceBalanceNetwork net0 = forceBalancer2.getNetwork();
    for (int i = 0; i < net0.nrOfPartialSprings; i++) {
      std::cout << net0.springPartIndexA[i] << ", " << net0.springPartIndexB[i]
                << ": ";
      std::cout << springPartitions0[i] << std::endl;
    }
    std::cout << std::endl;

    SECTION("displaceLinksToMeanPosition = displaceToMeanPosition")
    {
      Eigen::VectorXd oneOverSpringPartitions0 =
        forceBalancer2.assembleOneOverSpringPartition(&net0, springPartitions0);
      Eigen::VectorXd displacements0 =
        Eigen::VectorXd::Zero(forceBalancer2.getNrOfLinks() * 3);
      Eigen::ArrayXi twoIndependentIndices = Eigen::ArrayXi::Zero(2);
      twoIndependentIndices << 0, 2;

      Eigen::ArrayXi coordinateIndicesOfTwoIndependentIndices =
        Eigen::ArrayXi::Zero(3 * twoIndependentIndices.size());

      for (int j = 0; j < twoIndependentIndices.size(); ++j) {
        for (int i = 0; i < 3; ++i) {
          coordinateIndicesOfTwoIndependentIndices[3 * j + i] =
            3 * twoIndependentIndices[j] + i;
        }
      }

      double maxDiff0 = forceBalancer2.displaceLinksToMeanPosition(
        &net0,
        displacements0,
        oneOverSpringPartitions0,
        coordinateIndicesOfTwoIndependentIndices,
        1.0);

      Eigen::VectorXd displacements1 =
        Eigen::VectorXd::Zero(forceBalancer2.getNrOfLinks() * 3);
      double maxDiff1 = 0.;
      for (int i = 0; i < twoIndependentIndices.size(); ++i) {
        maxDiff1 = std::max(
          forceBalancer2.displaceToMeanPosition(
            &net0, displacements1, springPartitions0, twoIndependentIndices[i]),
          maxDiff1);
      }

      for (int i = 0; i < coordinateIndicesOfTwoIndependentIndices.size();
           ++i) {
        CHECK(displacements1[coordinateIndicesOfTwoIndependentIndices[i]] ==
              Catch::Approx(
                displacements0[coordinateIndicesOfTwoIndependentIndices[i]]));
      }

      for (int i = 0; i < displacements1.size(); ++i) {
        CHECK(displacements1[i] == Catch::Approx(displacements0[i]));
      }
      CHECK(maxDiff0 == Catch::Approx(maxDiff1));
    }

    SECTION("displaceLinksToMeanPosition works")
    {
      // test the "displaceLinksToMeanPosition"
      // see also: "Eigen behaves as expected"
      Eigen::VectorXd displacements0 =
        Eigen::VectorXd::Zero(forceBalancer2.getNrOfLinks() * 3);
      forceBalancer2.displaceLinksToMeanPosition(
        &net0, displacements0, springPartitions0, 0.);
      for (int i = 0; i < 3 * net0.nrOfLinks; ++i) {
        CHECK(displacements0[i] + 1e-5 == Catch::Approx(0.0 + 1e-5));
      }
      // repeat the displacement to reach the point where
      // all beads are at 0.0
      for (size_t it = 0; it < 20; ++it) {
        forceBalancer2.displaceLinksToMeanPosition(
          &net0, displacements0, springPartitions0, 0.5);
      }
      for (int i = 0; i < net0.nrOfLinks; ++i) {
        CHECK(displacements0[3 * i] + 1e-5 ==
              Catch::Approx(-net0.coordinates[3 * i] + 1e-5));
        CHECK(displacements0[3 * i + 1] + 1e-5 ==
              Catch::Approx(-net0.coordinates[3 * i + 1] + 1e-5));
        CHECK(displacements0[3 * i + 2] + 1e-5 == Catch::Approx(0.0 + 1e-5));
      }
    }

    SECTION("Irrelevant slip-links are rationalised")
    {
      // and add slip-links
      REQUIRE_THROWS(forceBalancer2.addSlipLinks({ { 1, 2, 3 } },
                                                 { { 1, 2, 3 } },
                                                 { { 1., 2., 3. } },
                                                 { { 1., 2., 3. } },
                                                 { { 1., 2. } }));
      REQUIRE_THROWS(forceBalancer2.addSlipLinks({ { 1, 2, 3 } },
                                                 { { 1, 2 } },
                                                 { { 1., 2., 3. } },
                                                 { { 1., 2., 3. } },
                                                 { { 1., 2., 3.0 } }));
      // and actually add slip-links
      REQUIRE_NOTHROW(forceBalancer2.addSlipLinks({ { 0, 3 } },
                                                  { { 1, 3 } },
                                                  { { 4.2, 1.3 } },
                                                  { { 3.9, 1.2 } },
                                                  { { 1.3, 1.1 } }));
      REQUIRE(forceBalancer2.getNrOfNodes() + 2 ==
              forceBalancer2.getNrOfLinks());
      // and run with them.
      // Expect the slip-link of between two strands to converge to the central
      // atom Expect the slip-link of two strands to stay at 0.5, 0.5.
      Eigen::VectorXd springPartitions = forceBalancer2.getSpringPartitions();
      Eigen::VectorXd displacements =
        Eigen::VectorXd::Zero(forceBalancer2.getNrOfLinks() * 3);
      pcm::ForceBalanceNetwork net = forceBalancer2.getNetwork();
      for (int i = 0; i < net.nrOfPartialSprings; i++) {
        std::cout << net.springPartIndexA[i] << ", " << net.springPartIndexB[i]
                  << ": ";
        std::cout << springPartitions[i] << std::endl;
      }
      std::cout << std::endl;

      forceBalancer2.displaceToMeanPosition(
        &net, displacements, springPartitions, 4);
      CHECK(net.coordinates[4 * 3] + displacements[4 * 3] ==
            Catch::Approx(2.5));
      CHECK(net.coordinates[4 * 3 + 1] + displacements[4 * 3 + 1] ==
            Catch::Approx(2.5));
      CHECK(net.coordinates[4 * 3 + 2] + displacements[4 * 3 + 2] ==
            Catch::Approx(2));

      // reset
      displacements.setZero();

      REQUIRE(net.springIndicesOfLinks.size() == net.nrOfLinks);
      for (int j = 0; j < 5; ++j) {
        forceBalancer2.displaceToMeanPosition(
          &net, displacements, springPartitions, 5, -1.);
        forceBalancer2.updateSpringPartition(
          &net, displacements, springPartitions, 5, -1.);
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
        forceBalancer2.displaceToMeanPosition(
          &net, displacements, springPartitions, 4, -1.);
        // for (int dir = 0; dir < 3; ++dir) {
        //   std::cout << (net.coordinates[3 * 4 + dir] + displacements[3 * 4 +
        //   dir])
        //             << ", ";
        // }
        forceBalancer2.updateSpringPartition(
          &net, displacements, springPartitions, 4, -1.);
        // std::cout << std::endl;
        // std::cout << springPartitions[0][0] << ", " << springPartitions[1][0]
        //           << std::endl;
      }
      // assert expectations are met.
      // NOTE: difficulty: finding out which spring idx it actually is
      outputNetwork(net, displacements, forceBalancer2.getSpringPartitions());
      for (int i = 0; i < net.nrOfPartialSprings; i++) {
        std::cout << net.springPartIndexA[i] << ", " << net.springPartIndexB[i]
                  << ": ";
        std::cout << springPartitions[i] << std::endl;
      }
      std::cout << std::endl;
      CHECK(springPartitions[5] + 1e-2 ==
            Catch::Approx(0.0 + 1e-2).epsilon(1e-6));                 // 4-1
      CHECK(springPartitions[6] == Catch::Approx(1.0).epsilon(1e-6)); // 4-2
      CHECK(springPartitions[0] == Catch::Approx(1.0).epsilon(1e-6)); // 4-0
      CHECK(springPartitions[1] + 1e-2 ==
            Catch::Approx(0.0 + 1e-2).epsilon(1e-6)); // 1-4
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
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
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
      REQUIRE(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE_NOTHROW(forceBalancer2.runForceRelaxation());
      CHECK(forceBalancer2.getNrOfIterations() > 1);
    }
  }
}
TEST_CASE("Free chains collapse",
          "[analysis][MEHPForceBalance][NonGaussianSpringForceEvaluator]["
          "SimpleSpringMEHPForceEvaluator]")
{
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
  REQUIRE(universe.getNrOfAtoms() == nrOfBeads);
  REQUIRE(universe.getNrOfBonds() == nrOfBeads - 1);

  // now, check for every force evaluator, that the maximum entropy is when all
  // these beads overlap first, the gaussian spring one
  pcm::MEHPForceBalance forceBalancer =
    pcm::MEHPForceBalance(universe, 2, false);
  REQUIRE_NOTHROW(forceBalancer.runForceRelaxation(
    pcm::BalanceRunMode::ITERATIVE, 1.0, 50000, 1e-18));
  REQUIRE(forceBalancer.getNrOfIterations() > 0);
  CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  CHECK(forceBalancer.getNrOfActiveSprings() == 0);
  // CHECK(forceBalancer.getAverageSpringLength() ==
  //       Catch::Approx(0.0));
  CHECK(forceBalancer.getAverageSpringLength() >= 0.0);
  CHECK(forceBalancer.getAverageSpringLength() <= 3e-6);
  REQUIRE_NOTHROW(forceBalancer.validateNetwork());

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

  std::array<std::array<double, 3>, 3> stressTensorSimpleSpring =
    forceBalancer.getStressTensor();
  std::array<std::array<double, 3>, 3> stressTensorLangevin =
    forceBalancer.getStressTensor();
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      CHECK(stressTensorLangevin[i][j] + 1e-5 ==
            Catch::Approx(stressTensorSimpleSpring[i][j] + 1e-5));
    }
  }
}
