#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceRelaxation.h"
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
        pe::Universe universe2 = universeSeq.atIndex(0);
        BENCHMARK_ADVANCED("MEHP Balance Random 1.0 " + largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_RANDOM,
                                                  1.0);
                return forceBalancer3.getNrOfIterations();
            });
        };
        BENCHMARK_ADVANCED("MEHP Balance Random 0.75 " + largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_RANDOM,
                                                  0.75);
                return forceBalancer3.getNrOfIterations();
            });
        };
        BENCHMARK_ADVANCED("MEHP Balance Strand 1.0 " + largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_STRANDS,
                                                  1.0);
                return forceBalancer3.getNrOfIterations();
            });
        };
        BENCHMARK_ADVANCED("MEHP Balance Strand 0.75 " + largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_STRANDS,
                                                  0.75);
                return forceBalancer3.getNrOfIterations();
            });
        };
        BENCHMARK_ADVANCED("MEHP Balance All 0.5 " + largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_ALL, 0.5);
                return forceBalancer3.getNrOfIterations();
            });
        };
        BENCHMARK_ADVANCED("MEHP Balance All 0.75 " + largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_ALL, 0.75);
                return forceBalancer3.getNrOfIterations();
            });
        };
        BENCHMARK_ADVANCED("MEHP Balance Heuristic 1.0 " + largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::EIGEN_HEURISTIC,
                                                  1.0);
                return forceBalancer3.getNrOfIterations();
            });
        };
        BENCHMARK_ADVANCED("MEHP Balance Iterative 1.0 " + largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                forceBalancer3.runForceRelaxation(pcm::BalanceRunMode::ITERATIVE, 1.0);
                return forceBalancer3.getNrOfIterations();
            });
        };

        BENCHMARK_ADVANCED("Detection of heuristic independent vertices " +
                           largeInputFile)
        (Catch::Benchmark::Chronometer meter) {
            pcm::MEHPForceBalance forceBalancer3 =
                pcm::MEHPForceBalance(universe2, 2);
            meter.measure([&forceBalancer3] {
                pcm::ForceBalanceNetwork net = forceBalancer3.getNetwork();
                return forceBalancer3.getHeuristicallyIndependentCoordinateSets(net);
            });
        };
    }
}
