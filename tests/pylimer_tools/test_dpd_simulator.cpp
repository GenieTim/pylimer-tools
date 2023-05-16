#include "../../src/pylimer_tools_cpp/calc/DPDSimulator.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::calc::mehp;
namespace pcd = pylimer_tools::calc::dpd;

TEST_CASE("DPD Simulator Works", "[analysis][DPDSimulator]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  std::string inputFile = suspectedPath + "melt_83_a_100.structure.out";
  if (std::filesystem::exists(inputFile)) {
    REQUIRE(std::filesystem::exists(suspectedPath));
    universeSeq.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);

    pcd::DPDSimulator simulator = pcd::DPDSimulator(universe, 2, false, "seed");

    // configuration
    REQUIRE_NOTHROW(simulator.configA(25.));
    REQUIRE_NOTHROW(simulator.configSigma(3.));
    REQUIRE_NOTHROW(simulator.configSlipspringLowCutoff(0.5));
    REQUIRE_NOTHROW(simulator.configSlipspringHighCutoff(2.0));
    REQUIRE_THROWS(simulator.configSlipspringLowCutoff(3.0));
    REQUIRE_THROWS(simulator.configSlipspringLowCutoff(2.0));
    REQUIRE_THROWS(simulator.configSlipspringHighCutoff(0.5));

    std::vector<pcd::ComputedValues> outputQuantities = {
      pcd::ComputedValues::STEP,      pcd::ComputedValues::TEMPERATURE,
      pcd::ComputedValues::PRESSURE,  pcd::ComputedValues::STRESS_XX,
      pcd::ComputedValues::STRESS_YY, pcd::ComputedValues::STRESS_ZZ,
      pcd::ComputedValues::MSD
    };

    REQUIRE_NOTHROW(simulator.configValuesToOutput(outputQuantities));

    std::vector<pcd::ComputedValues> averageQuantities = {
      pcd::ComputedValues::TEMPERATURE, pcd::ComputedValues::PRESSURE,
      pcd::ComputedValues::STRESS_XX,   pcd::ComputedValues::STRESS_YY,
      pcd::ComputedValues::STRESS_ZZ,   pcd::ComputedValues::MSD
    };

    REQUIRE_NOTHROW(simulator.configValuesToAverage(averageQuantities));
    std::string averageFile =
      suspectedPath + "melt_83_a_100.structure.avg-out.txt";
    REQUIRE_NOTHROW(simulator.configAveragesFile(averageFile));

    // actual simulation
    REQUIRE_NOTHROW(simulator.runSimulation(75, 0.06, false));
    REQUIRE_NOTHROW(simulator.validateState());
    REQUIRE_NOTHROW(simulator.validateNeighbourlist(1.0));

    CHECK(simulator.getTemperature() == Catch::Approx(1.0).epsilon(0.5));

    simulator.createSlipSprings(100, 2);

    pe::Universe resultUniverse = simulator.getUniverse();
    CHECK(resultUniverse.getNrOfBonds() == 100 + universe.getNrOfBonds());
    CHECK(resultUniverse.getNrOfAtoms() == universe.getNrOfAtoms());

    REQUIRE_NOTHROW(simulator.runSimulation(75, 0.06, true));

    REQUIRE_NOTHROW(simulator.validateState());
    REQUIRE_NOTHROW(simulator.validateNeighbourlist(1.0));

    REQUIRE(std::filesystem::exists(averageFile));
    std::remove(averageFile.c_str());
  }
}
