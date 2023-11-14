#include "../../src/pylimer_tools_cpp/calc/DPDSimulator.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>

namespace pe = pylimer_tools::entities;
namespace pc = pylimer_tools::calc;
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

    std::vector<pc::ComputedDoubleValues> outputQuantities = {
      pc::ComputedDoubleValues::TEMPERATURE,
      pc::ComputedDoubleValues::PRESSURE,
      pc::ComputedDoubleValues::STRESS_XX,
      pc::ComputedDoubleValues::STRESS_YY,
      pc::ComputedDoubleValues::STRESS_ZZ,
      pc::ComputedDoubleValues::MSD
    };

    pc::OutputConfiguration config;
    config.filename = "";
    config.outputEvery = 5;
    config.doubleValues = outputQuantities;
    config.intValues = { pc::ComputedIntValues::STEP };

    std::vector<pc::OutputConfiguration> configs = { config };
    REQUIRE_NOTHROW(simulator.configStepOutput(configs));

    std::vector<pc::ComputedDoubleValues> averageQuantities = {
      pc::ComputedDoubleValues::TEMPERATURE,
      pc::ComputedDoubleValues::PRESSURE,
      pc::ComputedDoubleValues::STRESS_XX,
      pc::ComputedDoubleValues::STRESS_YY,
      pc::ComputedDoubleValues::STRESS_ZZ,
      pc::ComputedDoubleValues::MSD
    };

    std::string averageFile =
      suspectedPath + "melt_83_a_100.structure.avg-out.txt";
    pc::OutputConfiguration avgconfig;
    avgconfig.outputEvery = 20;
    avgconfig.filename = averageFile;
    avgconfig.doubleValues = averageQuantities;

    std::vector<pc::OutputConfiguration> avgconfigs = { avgconfig };
    REQUIRE_NOTHROW(simulator.configAverageOutput(avgconfigs));

    std::vector<pc::ComputedDoubleValues> autocorrelationQuantities = {
      pc::ComputedDoubleValues::STRESS_XX,
      pc::ComputedDoubleValues::STRESS_YY,
      pc::ComputedDoubleValues::STRESS_ZZ,
      pc::ComputedDoubleValues::STRESS_XY,
      pc::ComputedDoubleValues::STRESS_YZ,
      pc::ComputedDoubleValues::STRESS_XZ,
      pc::ComputedDoubleValues::STRESS_NXY,
      pc::ComputedDoubleValues::STRESS_NYZ,
      pc::ComputedDoubleValues::STRESS_NXZ,
    };
    std::string autocorrFile =
      suspectedPath + "melt_83_a_100.structure.autocorr-out.txt";
    pc::OutputConfiguration autocorrconfig;
    autocorrconfig.outputEvery = 25;
    autocorrconfig.filename = autocorrFile;
    autocorrconfig.doubleValues = autocorrelationQuantities;

    std::vector<pc::OutputConfiguration> autocorrconfigs = { autocorrconfig };
    REQUIRE_NOTHROW(simulator.configAutoCorrelatorOutput(autocorrconfigs));

    std::vector<size_t> atomIdsForMSD = { 1, 4, 6 };
    REQUIRE_NOTHROW(simulator.startMeasuringMSDForAtoms(atomIdsForMSD));

    // restart files
    std::string restartFile = suspectedPath + "dpd_restart_file.bin";
    simulator.configRestartOutput(restartFile, 20);

    // actual simulation
    REQUIRE_NOTHROW(simulator.runSimulation(75, false));
    REQUIRE_NOTHROW(simulator.validateState());
    std::cout << "DPD ran, state validated." << std::endl;
    CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
    CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));

    CHECK(simulator.getTemperature() == Catch::Approx(1.0).epsilon(0.5));
    CHECK_NOTHROW(simulator.validateState());

    simulator.createSlipSprings(100, 2);
    CHECK_NOTHROW(simulator.validateState());
    std::cout << "DPD slip-springs created." << std::endl;

    pe::Universe resultUniverse = simulator.getUniverse();
    CHECK(resultUniverse.getNrOfBonds() == 100 + universe.getNrOfBonds());
    CHECK(resultUniverse.getNrOfAtoms() == universe.getNrOfAtoms());
    CHECK_NOTHROW(simulator.validateState());

    CHECK_NOTHROW(simulator.runSimulation(76, true));

    CHECK_NOTHROW(simulator.validateState());
    std::cout << "DPD ran with slip-springs, state validated." << std::endl;
    CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));

    CHECK(std::filesystem::exists(averageFile));
    std::remove(averageFile.c_str());
    CHECK(std::filesystem::exists(autocorrFile));
    std::remove(autocorrFile.c_str());
    REQUIRE(std::filesystem::exists(restartFile));

    pcd::DPDSimulator sim2 = pcd::DPDSimulator::readRestartFile(restartFile);
    REQUIRE_NOTHROW(sim2.runSimulation(5, false));
    CHECK(sim2.getTemperature() ==
          Catch::Approx(simulator.getTemperature()).epsilon(0.1));
    CHECK_NOTHROW(sim2.validateState());
    std::remove(restartFile.c_str());
  }
};
