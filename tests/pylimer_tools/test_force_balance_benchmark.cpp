#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance2.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceRelaxation.h"
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

TEST_CASE("Force Balance Benchmarks", "[analysis][MEHPForceBalance]")
{
  std::cout << "Running test \"Force Balance Benchmarks\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  std::string largeInputFile =
    suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
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
    std::cout << "Reference Time (FB1) to beat: " << duration_ref.count()
              << " microseconds" << std::endl;

    BENCHMARK_ADVANCED("MEHP LD_MMA " +
                       largeInputFile)(Catch::Benchmark::Chronometer meter)
    {
      pcm::MEHPForceRelaxation forceRelaxer =
        pcm::MEHPForceRelaxation(universe, 2);
      forceRelaxer.configAssumeBoxLargeEnough(false);

      meter.measure([&forceRelaxer, &referenceForceBalancer] {
        forceRelaxer.runForceRelaxation("LD_MMA");
        CHECK_THAT(
          forceRelaxer.getGamma(),
          Catch::Matchers::WithinRel(referenceForceBalancer.getGamma()));
        return forceRelaxer.getNrOfIterations();
      });
    };

    for (pcm::SLESolver solverChoice :
         { // pcm::SIMPLICIAL_LLT,
           //                                    pcm::SIMPLICIAL_DLT,
           //                                    pcm::SPARSE_LU,
           //                                    pcm::SPARSE_QR,
           pcm::CONJUGATE_GRADIENT,
           pcm::LEAST_SQUARES_CONJUGATE_GRADIENT,
           pcm::BICGSTAB }) {
      pcm::MEHPForceBalance2 forceBalancer =
        pcm::MEHPForceBalance2(universe, 2);
      forceBalancer.configAssumeBoxLargeEnough(false);

      auto start = std::chrono::high_resolution_clock::now();

      try {
        forceBalancer.runForceRelaxation(
          pcm::StructureSimplificationMode::NO_SIMPLIFICATION,
          1e-3,
          1.0,
          solverChoice);
      } catch (const std::exception& e) {
        std::cerr << "Exception for solver " << solverChoice << ": " << e.what()
                  << std::endl;
        continue;
      }

      CHECK_THAT(forceBalancer.getGamma(),
                 Catch::Matchers::WithinRel(referenceForceBalancer.getGamma()));
      auto end = std::chrono::high_resolution_clock::now();
      auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      std::cout << "Solver (FB2): " << solverChoice
                << ", Time: " << duration.count() << " microseconds"
                << std::endl;
    }
  }
}
