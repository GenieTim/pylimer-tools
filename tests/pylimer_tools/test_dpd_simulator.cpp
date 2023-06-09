#include "../../src/pylimer_tools_cpp/calc/DPDSimulator.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <filesystem>
#include <string>

namespace pe = pylimer_tools::entities;
namespace pcd = pylimer_tools::calc::dpd;

TEST_CASE("DPD Simulator Works", "[analysis][DPDSimulator]")
{
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile = suspectedPath + "melt_83_a_100.structure.out";
  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator = pcd::DPDSimulator(universe, 2, false, "seed");

    // configuration
    REQUIRE_NOTHROW(simulator.validateState());
    CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
    CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));
    REQUIRE_NOTHROW(simulator.configA(25.));
    REQUIRE_NOTHROW(simulator.configSigma(3.));
    REQUIRE_NOTHROW(simulator.configSlipspringLowCutoff(0.5));
    REQUIRE_NOTHROW(simulator.configSlipspringHighCutoff(2.0));
    REQUIRE_THROWS(simulator.configSlipspringLowCutoff(3.0));
    REQUIRE_THROWS(simulator.configSlipspringLowCutoff(2.0));
    REQUIRE_THROWS(simulator.configSlipspringHighCutoff(0.5));

    std::vector<pcd::ComputedDoubleValues> outputQuantities = {
      pcd::ComputedDoubleValues::TEMPERATURE,
      pcd::ComputedDoubleValues::PRESSURE,  pcd::ComputedDoubleValues::STRESS_XX,
      pcd::ComputedDoubleValues::STRESS_YY, pcd::ComputedDoubleValues::STRESS_ZZ,
      pcd::ComputedDoubleValues::MSD
    };

    pcd::OutputConfiguration config;
    config.filename = "";
    config.outputEvery = 5;
    config.doubleValues = outputQuantities;
    config.intValues = { pcd::ComputedIntValues::STEP };

    REQUIRE_NOTHROW(simulator.configStepOutput({config}));

    std::vector<pcd::ComputedDoubleValues> averageQuantities = {
      pcd::ComputedDoubleValues::TEMPERATURE, pcd::ComputedDoubleValues::PRESSURE,
      pcd::ComputedDoubleValues::STRESS_XX,   pcd::ComputedDoubleValues::STRESS_YY,
      pcd::ComputedDoubleValues::STRESS_ZZ,   pcd::ComputedDoubleValues::MSD
    };

    std::string averageFile =
      suspectedPath + "melt_83_a_100.structure.avg-out.txt";
    pcd::OutputConfiguration avgconfig;
    avgconfig.outputEvery = 20;
    avgconfig.filename = averageFile;
    avgconfig.doubleValues = averageQuantities;

    REQUIRE_NOTHROW(simulator.configAverageOutput({avgconfig}));

    std::vector<size_t> atomIdsForMSD = { 1, 4, 6 };
    REQUIRE_NOTHROW(simulator.startMeasuringMSDForAtoms(atomIdsForMSD));

    // actual simulation
    REQUIRE_NOTHROW(simulator.runSimulation(75, 0.06, false));
    REQUIRE_NOTHROW(simulator.validateState());
    CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
    CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));

    CHECK(simulator.getTemperature() == Catch::Approx(1.0).epsilon(0.5));
    CHECK_NOTHROW(simulator.validateState());

    simulator.createSlipSprings(100, 2);
    CHECK_NOTHROW(simulator.validateState());

    pe::Universe resultUniverse = simulator.getUniverse();
    CHECK(resultUniverse.getNrOfBonds() == 100 + universe.getNrOfBonds());
    CHECK(resultUniverse.getNrOfAtoms() == universe.getNrOfAtoms());
    CHECK_NOTHROW(simulator.validateState());

    CHECK_NOTHROW(simulator.runSimulation(75, 0.06, true));

    CHECK_NOTHROW(simulator.validateState());
    CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));

    REQUIRE(std::filesystem::exists(averageFile));
    std::remove(averageFile.c_str());
  }
};
