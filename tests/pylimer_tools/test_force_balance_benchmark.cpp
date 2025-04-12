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

TEST_CASE("Force Balance Benchmarks", "[MEHPForceBalance2][benchmark][long]")
{
  std::cout << "Running test \"Force Balance Benchmarks\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::string largeInputFile =
    suspectedPath +
    "/structure/xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
  if (std::filesystem::exists(largeInputFile)) {
    universeSeq.initializeFromDataSequence({ { largeInputFile } });
    pe::Universe universe = universeSeq.atIndex(0);

    pcm::MEHPForceBalance referenceForceBalancer =
      pcm::MEHPForceBalance(universe, 2);
    referenceForceBalancer.configAssumeBoxLargeEnough(false);

    auto start_ref = std::chrono::high_resolution_clock::now();

    referenceForceBalancer.runForceRelaxation();
    auto end_ref = std::chrono::high_resolution_clock::now();
    auto duration_ref = std::chrono::duration_cast<std::chrono::microseconds>(
      end_ref - start_ref);
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
         { // pcm::SIMPLICIAL_LLT,
           //                                    pcm::SIMPLICIAL_DLT,
           //                                    pcm::SPARSE_LU,
           //                                    pcm::SPARSE_QR,
           pcm::CONJUGATE_GRADIENT,
           pcm::CONJUGATE_GRADIENT_IDENTITY,
           pcm::LEAST_SQUARES_CONJUGATE_GRADIENT,
           pcm::BICGSTAB }) {
      pcm::MEHPForceBalance2 forceBalancer =
        pcm::MEHPForceBalance2(universe, 2);

      auto start = std::chrono::high_resolution_clock::now();

      try {
        forceBalancer.runForceRelaxation(
          pcm::StructureSimplificationMode::NO_SIMPLIFICATION,
          1e-3,
          solverChoice);
      } catch (const std::exception& e) {
        std::cerr << "Exception for solver " << solverChoice << ": " << e.what()
                  << std::endl;
        continue;
      }

      CHECK_THAT(
        forceBalancer.getGamma(),
        Catch::Matchers::WithinRel(referenceForceBalancer.getGamma(), 1e-4));
      auto end = std::chrono::high_resolution_clock::now();
      auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);

      std::cout << "Solver (FB2): " << solverChoice
                << ", Time: " << std::duration_to_string(duration) << std::endl;
    }
  }
}

TEST_CASE("Force Balance Benchmarks randomly functionalized",
          "[MEHPForceBalance][MEHPForceBalance2][benchmark][long]")
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
       { // pcm::SIMPLICIAL_LLT,
         //                                    pcm::SIMPLICIAL_DLT,
         //                                    pcm::SPARSE_LU,
         //                                    pcm::SPARSE_QR,
         pcm::CONJUGATE_GRADIENT,
         pcm::CONJUGATE_GRADIENT_IDENTITY,
         pcm::LEAST_SQUARES_CONJUGATE_GRADIENT,
         pcm::BICGSTAB }) {
    pcm::MEHPForceBalance2 forceBalancer = pcm::MEHPForceBalance2(universe, 2);

    auto start = std::chrono::high_resolution_clock::now();

    try {
      forceBalancer.runForceRelaxation(
        pcm::StructureSimplificationMode::NO_SIMPLIFICATION,
        1e-3,
        solverChoice);
    } catch (const std::exception& e) {
      std::cerr << "Exception for solver " << solverChoice << ": " << e.what()
                << std::endl;
      continue;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    CHECK_THAT(
      forceBalancer.getGamma(),
      Catch::Matchers::WithinAbs(referenceForceBalancer.getGamma(), 1e-5));
    CHECK(forceBalancer.getNrOfActiveNodes() == 0);

    std::cout << "Solver (FB2): " << solverChoice
              << ", Time: " << std::duration_to_string(duration) << std::endl;
  }
}