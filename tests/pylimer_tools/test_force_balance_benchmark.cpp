#include "../../src/pylimer_tools_cpp/calc/MEHPanalysis.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance2.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/utils/MCUniverseGenerator.h"
#include "../../src/pylimer_tools_cpp/utils/StringUtils.h"
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
namespace pcm = pylimer_tools::sim::mehp;

#ifndef PYLIMER_TEST_FIXTURES_DIR
#define PYLIMER_TEST_FIXTURES_DIR "../pylimer_tools/tests/fixtures"
#endif

TEST_CASE("Force Balance Benchmarks",
          "[MEHPForceBalance2][benchmark][long][FBBenchmark]")
{
  std::cout << "Running test \"Force Balance Benchmarks\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::string largeInputFile =
    suspectedPath +
    "/structure/xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
  REQUIRE(std::filesystem::exists(largeInputFile));
  universeSeq.initializeFromDataSequence({ { largeInputFile } });
  pe::Universe universe = universeSeq.atIndex(0);

  pcm::MEHPForceBalance referenceForceBalancer =
    pcm::MEHPForceBalance(universe, 2);
  referenceForceBalancer.configAssumeBoxLargeEnough(false);

  auto start_ref = std::chrono::high_resolution_clock::now();

  referenceForceBalancer.runForceRelaxation();
  auto end_ref = std::chrono::high_resolution_clock::now();
  auto duration_ref =
    std::chrono::duration_cast<std::chrono::microseconds>(end_ref - start_ref);
  std::cout << "Reference Time (FB1) to beat: "
            << std::duration_to_string(duration_ref) << " " << std::endl;

  // BENCHMARK_ADVANCED("MEHP LD_MMA " +
  //                    largeInputFile)(Catch::Benchmark::Chronometer meter)
  // {
  //   pcm::MEHPForceRelaxation forceRelaxer =
  //     pcm::MEHPForceRelaxation(universe, 2);
  //   forceRelaxer.configAssumeBoxLargeEnough(false);

  //   meter.measure([&forceRelaxer, &referenceForceBalancer] {
  //     forceRelaxer.runForceRelaxation("LD_MMA");

  //     CHECK_THAT(
  //       forceRelaxer.getGamma(),
  //       Catch::Matchers::WithinRel(referenceForceBalancer.getGamma(),
  //       1e-2));
  //     return forceRelaxer.getNrOfIterations();
  //   });
  // };

  for (pcm::SLESolver solverChoice :
       // { // pcm::SIMPLICIAL_LLT,
       //   pcm::SIMPLICIAL_LDLT,
       //   //                                    pcm::SPARSE_LU,
       //   //                                    pcm::SPARSE_QR,
       //   pcm::CONJUGATE_GRADIENT,
       //   pcm::CONJUGATE_GRADIENT_IDENTITY,
       //   pcm::LEAST_SQUARES_CONJUGATE_GRADIENT,
       //   pcm::BICGSTAB }
       pylimer_tools::sim::mehp::allSLESolvers) {

    // need to exclude some very slow cases
    if (solverChoice == pcm::SLESolver::SPARSE_QR) {
      continue;
    }

    pcm::MEHPForceBalance2
      forceBalancer = // pcm::MEHPForceBalance2(universe, 2);
      pcm::MEHPForceBalance2(universe,
                             referenceForceBalancer.getNetwork(),
                             referenceForceBalancer.getSpringPartitions(),
                             false);
    auto start = std::chrono::high_resolution_clock::now();

    try {
      forceBalancer.runForceRelaxation(
        pcm::StructureSimplificationMode::NO_SIMPLIFICATION,
        1e-3,
        solverChoice);
    } catch (const std::exception& e) {
      std::cerr << "Exception for solver " << pcm::SLESolverNames[solverChoice]
                << " (" << solverChoice << "): " << e.what() << std::endl;
      continue;
    }

    CHECK_THAT(
      forceBalancer.getGamma(),
      Catch::Matchers::WithinRel(referenceForceBalancer.getGamma(), 1e-4));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Solver (FB2): " << pcm::SLESolverNames[solverChoice] << " ("
              << solverChoice
              << "), Time: " << std::duration_to_string(duration) << std::endl;
  }
}

TEST_CASE("Force Balance Benchmarks randomly functionalized",
          "[MEHPForceBalance][MEHPForceBalance2][benchmark][long][FBBenchmark]")
{
  std::cout
    << "Running test \"Force Balance Benchmarks randomly functionalized\""
    << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(457564875e2);
  generator.setBeadDistance(0.75);
  generator.configNrOfMCSteps(0);

  std::vector<int> chainLengths = pu::initializeWithValue(50, 100);
  generator.addRandomlyFunctionalizedStrands(
    50, chainLengths, 7.2, 1, 2, 1, true);

  CHECK_THAT(generator.getCurrentNrOfAvailableCrosslinkSites(),
             Catch::Matchers::WithinRel(50 * 100 * 7.2, 0.01));
  CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
             Catch::Matchers::WithinAbs(0.0, 1e-10));

  std::vector<int> monofunctionalChainLengths =
    pu::initializeWithValue(200, 10);
  generator.addMonofunctionalStrands(200, monofunctionalChainLengths, 1);

  REQUIRE_NOTHROW(generator.linkStrandsToConversion(
    199. / generator.getCurrentNrOfAvailableCrosslinkSites()));

  pe::Universe universe = generator.getUniverse();

  pcm::MEHPForceBalance referenceForceBalancer = pcm::MEHPForceBalance(
    universe,
    2); // using a different crosslinkerType here makes things faster
  referenceForceBalancer.configAssumeBoxLargeEnough(false);

  auto start_ref = std::chrono::high_resolution_clock::now();

  referenceForceBalancer.runForceRelaxation();
  auto end_ref = std::chrono::high_resolution_clock::now();
  auto duration_ref =
    std::chrono::duration_cast<std::chrono::microseconds>(end_ref - start_ref);
  std::cout << "Reference Time (FB1) to beat: "
            << std::duration_to_string(duration_ref) << " " << std::endl;
  CHECK(referenceForceBalancer.getNrOfActiveNodes(1e-1) == 0);

  // BENCHMARK_ADVANCED("MEHP LD_MMA " +
  //                    largeInputFile)(Catch::Benchmark::Chronometer meter)
  // {
  //   pcm::MEHPForceRelaxation forceRelaxer =
  //     pcm::MEHPForceRelaxation(universe, 2);

  //   meter.measure([&forceRelaxer, &referenceForceBalancer] {
  //     forceRelaxer.runForceRelaxation("LD_MMA");
  //     CHECK(forceRelaxer.getGamma() == referenceForceBalancer.getGamma());
  //     return forceRelaxer.getNrOfIterations();
  //   });
  // };

  for (pcm::SLESolver solverChoice :
       // { // pcm::SIMPLICIAL_LLT,
       //   pcm::SIMPLICIAL_LDLT,
       //   //                                    pcm::SPARSE_LU,
       //   //                                    pcm::SPARSE_QR,
       //   pcm::CONJUGATE_GRADIENT,
       //   pcm::CONJUGATE_GRADIENT_IDENTITY,
       //   pcm::LEAST_SQUARES_CONJUGATE_GRADIENT,
       //   pcm::BICGSTAB }
       pylimer_tools::sim::mehp::allSLESolvers) {

    // need to exclude some very slow cases
    if (solverChoice == pcm::SLESolver::SPARSE_QR) {
      continue;
    }

    pcm::MEHPForceBalance2
      forceBalancer = // pcm::MEHPForceBalance2(universe, 2);
      pcm::MEHPForceBalance2(universe,
                             referenceForceBalancer.getNetwork(),
                             referenceForceBalancer.getSpringPartitions(),
                             false);
    auto start = std::chrono::high_resolution_clock::now();

    try {
      forceBalancer.runForceRelaxation(
        pcm::StructureSimplificationMode::NO_SIMPLIFICATION,
        1e-3,
        solverChoice);
    } catch (const std::exception& e) {
      std::cerr << "Exception for solver " << pcm::SLESolverNames[solverChoice]
                << " (" << solverChoice << "): " << e.what() << std::endl;
      continue;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    CHECK_THAT(
      forceBalancer.getGamma(),
      Catch::Matchers::WithinAbs(referenceForceBalancer.getGamma(), 1e-5));
    CHECK(forceBalancer.getNrOfActiveNodes() == 0);

    std::cout << "Solver (FB2): " << pcm::SLESolverNames[solverChoice] << " ("
              << solverChoice
              << "), Time: " << std::duration_to_string(duration) << std::endl;
  }
}

TEST_CASE("Temporary entanglement sampling benchmark",
          "[EntanglementDetector][benchmark]")
{
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::vector<std::string> files = {
    "3d-diamond-lattice_10x10x10_a_3_d_0.85_imperfect.structure.out",
    "3d-diamond-lattice_10x10x10_a_3_d_0.85_v_0.V-fixed.structure.out",
    "3d-diamond-lattice_3x3x3_a_23_d_3_v_0.structure.out",
    "3d-diamond-lattice_5x5x5_a_3_d_0.85_imperfect.structure.out",
    "3d-diamond-lattice_5x5x5_a_3_d_0.85_v_0.V-fixed.structure.out",
    "crosslinked_M10000_N39_p_0.9.out",
    "crosslinked_p_0.98_melt_100_a_3_50_xlinks_v_14.converted.structure.out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_0.98_melt_100_a_38_50_xlinks_v_22.structure.out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out",
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure."
    "out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_1_1_melt_100_a_3_50_xlinks_v_1.V-fixed.structure.out-"
    "finish_"
    "crosslinking.structure.out",
    "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0.333_2d_t_7500001."
    "structure.out",
    "mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0.0114_f_4_v_1."
    "structure.out",
    "melt_213_a_47_106_xlinks_v_1.structure.out",
    "melt_83_a_100.structure.out",
    "network_100_a_46.structure.out",
    "network_p_1_100_a_38_50_xlinks.structure.out",
    "square_lattice_2x2_a_5.2d.structure.out",
    "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out",
    "xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out"
  };

  for (const std::string& file : files) {
    pe::UniverseSequence universeSeq = pe::UniverseSequence();
    REQUIRE(universeSeq.getLength() == 0);
    std::string inputFile = suspectedPath + "/structure/" + file;
    if (!std::filesystem::exists(inputFile)) {
      std::cerr << "File not found: " << inputFile << std::endl;
      continue;
    }

    std::cout << "Processing file: " << file << std::endl;

    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);

    for (const bool filtered : { true, false }) {
      for (const double sameStrandCutoff : { -1., 2., 5. }) {
        auto start_ref = std::chrono::high_resolution_clock::now();

        pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
          entanglements = pylimer_tools::topo::entanglement_detection::
            randomlyFindEntanglements(universe,
                                      0.1 * universe.getNrOfAtoms(),
                                      2.0,
                                      0.,
                                      0,
                                      sameStrandCutoff,
                                      "",
                                      2,
                                      true,
                                      filtered);

        auto end_ref = std::chrono::high_resolution_clock::now();
        auto duration_ref =
          std::chrono::duration_cast<std::chrono::microseconds>(end_ref -
                                                                start_ref);
        std::cout << "Entanglements v1, " << (filtered ? "" : "un")
                  << "filtered, " << sameStrandCutoff << " same-strand cutoff: "
                  << std::duration_to_string(duration_ref) << " " << std::endl;
        CHECK(entanglements.pairsOfAtoms.size() >=
              0.05 * universe.getNrOfAtoms());
        CHECK(entanglements.pairsOfAtoms.size() <=
              0.11 * universe.getNrOfAtoms());

        auto start_v2 = std::chrono::high_resolution_clock::now();

        pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
          entanglements2 = pylimer_tools::topo::entanglement_detection::
            randomlyFindEntanglementsV2(universe,
                                        0.1 * universe.getNrOfAtoms(),
                                        2.0,
                                        0.,
                                        0,
                                        sameStrandCutoff,
                                        "",
                                        2,
                                        true,
                                        filtered);

        auto end_v2 = std::chrono::high_resolution_clock::now();
        auto duration_v2 =
          std::chrono::duration_cast<std::chrono::microseconds>(end_v2 -
                                                                start_v2);
        std::cout << "Entanglements v2, " << (filtered ? "" : "un")
                  << "filtered, " << sameStrandCutoff << " same-strand cutoff: "
                  << std::duration_to_string(duration_v2) << " " << std::endl;
        CHECK(entanglements2.pairsOfAtoms.size() >=
              0.05 * universe.getNrOfAtoms());
        CHECK(entanglements2.pairsOfAtoms.size() <=
              0.11 * universe.getNrOfAtoms());
      }
    }
  }
}