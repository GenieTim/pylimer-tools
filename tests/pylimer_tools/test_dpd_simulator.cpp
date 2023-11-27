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

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, false, "12th_seed");

    // configuration
    REQUIRE_NOTHROW(simulator.validateState());
    // CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
    // CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));
    REQUIRE_NOTHROW(simulator.configA(25.));
    REQUIRE_NOTHROW(simulator.configSigma(3.));
    REQUIRE_NOTHROW(simulator.configSpringConstant(1.));
    REQUIRE_NOTHROW(simulator.configSlipspringLowCutoff(0.5));
    REQUIRE_NOTHROW(simulator.configSlipspringHighCutoff(2.0));

    // check that wrong configuration throws
    REQUIRE_THROWS(simulator.configSlipspringLowCutoff(3.0));
    REQUIRE_THROWS(simulator.configSlipspringLowCutoff(2.0));
    REQUIRE_THROWS(simulator.configSlipspringHighCutoff(0.5));

    // initial state – vgl. lammps
    CHECK(simulator.getStressTensor().trace() / 3. == Catch::Approx(23.321285));
    CHECK(simulator.getTemperature() + 1e-2 == Catch::Approx(0. + 1e-2));
    Eigen::Matrix3d initialStressTensor = simulator.getStressTensor();
    CHECK(initialStressTensor(0, 0) == Catch::Approx(23.456778).epsilon(0.75));
    CHECK(initialStressTensor(1, 1) == Catch::Approx(23.374175).epsilon(0.75));
    CHECK(initialStressTensor(2, 2) == Catch::Approx(23.401583).epsilon(0.75));
    CHECK(initialStressTensor(0, 1) ==
          Catch::Approx(-0.0085286852).epsilon(0.002));
    CHECK(initialStressTensor(0, 2) ==
          Catch::Approx(-0.10797027).epsilon(0.02));
    CHECK(initialStressTensor(1, 2) ==
          Catch::Approx(-0.035841973).epsilon(0.02));

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
    // std::string restartFile = suspectedPath + "dpd_restart_file.xml";
    std::string restartFile = suspectedPath + "dpd_restart_file.bin";
    simulator.configRestartOutput(restartFile, 30);

    // actual simulation – without slip-springs yet
    REQUIRE_NOTHROW(simulator.runSimulation(75, false));
    REQUIRE_NOTHROW(simulator.validateState());
    CHECK(simulator.getStressTensor().trace() / 3. ==
          Catch::Approx(20.8).epsilon(1.));
    std::cout << "DPD ran, state validated." << std::endl;
    // CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
    // CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));

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

TEST_CASE("DPD Simulator Can Cross-link", "[analysis][DPDSimulator]")
{
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "melt_213_a_47_106_xlinks_v_1.structure.out";
  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, false, "seed2");

    // configuration
    REQUIRE_NOTHROW(simulator.validateState());
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

    std::unordered_map<int, int> numBondsPerType;
    numBondsPerType[1] = 2;
    numBondsPerType[2] = 4;
    REQUIRE_NOTHROW(
      simulator.configBondFormation(212, numBondsPerType, 2.2, 5));
    size_t numBondsBefore = simulator.getNumBonds();

    // actual simulation
    REQUIRE_NOTHROW(simulator.runSimulation(200, false));
    REQUIRE_NOTHROW(simulator.validateState());
    std::cout << "DPD ran, state validated." << std::endl;

    // check that we actually have formed bonds
    CHECK(simulator.getNumBonds() >= (100 + numBondsBefore));
    pe::Universe universeAfter = simulator.getUniverse();
    std::map<int, int> finalFunctionalityPerType =
      universeAfter.determineFunctionalityPerType();
    CHECK(finalFunctionalityPerType.at(1) == 2);
    CHECK(finalFunctionalityPerType.at(2) <= 4);
  }
};

TEST_CASE("DPD Simulator Computes Correct Forces",
          "[analysis][DPDSimulator]")
{
  // note that the random force might lead to deviations compared to LAMMPS
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile = suspectedPath + "melt_83_a_100.structure.out";
  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, false, "12th_seed");

    // configuration
    REQUIRE_NOTHROW(simulator.validateState());
    // CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
    // CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));
    REQUIRE_NOTHROW(simulator.configA(0.));
    REQUIRE_NOTHROW(simulator.configSigma(0.));
    REQUIRE_NOTHROW(simulator.configSpringConstant(0.));

    REQUIRE_NOTHROW(simulator.refreshCurrentState());
    CHECK(simulator.getStressTensor().trace() / 3. == Catch::Approx(0.));
    // gradually add forces
    REQUIRE_NOTHROW(simulator.configA(25.));
    REQUIRE_NOTHROW(simulator.refreshCurrentState());
    CHECK(simulator.getStressTensor().trace() / 3. == Catch::Approx(25.259583));

    REQUIRE_NOTHROW(simulator.configA(0.));
    REQUIRE_NOTHROW(simulator.configSpringConstant(1.));
    REQUIRE_NOTHROW(simulator.refreshCurrentState());
    CHECK(simulator.getStressTensor().trace() / 3. == Catch::Approx(-1.78695));

    REQUIRE_NOTHROW(simulator.configSigma(3.));
    REQUIRE_NOTHROW(simulator.refreshCurrentState());
    CHECK(simulator.getStressTensor().trace() / 3. ==
          Catch::Approx(-1.8487375));
    // initial state – vgl. lammps
    REQUIRE_NOTHROW(simulator.configA(25.));
    CHECK(simulator.getStressTensor().trace() / 3. == Catch::Approx(23.321285));
    CHECK(simulator.getTemperature() + 1e-2 == Catch::Approx(0. + 1e-2));
    Eigen::Matrix3d initialStressTensor = simulator.getStressTensor();
    CHECK(initialStressTensor(0, 0) == Catch::Approx(23.456778).epsilon(0.75));
    CHECK(initialStressTensor(1, 1) == Catch::Approx(23.374175).epsilon(0.75));
    CHECK(initialStressTensor(2, 2) == Catch::Approx(23.401583).epsilon(0.75));
    CHECK(initialStressTensor(0, 1) ==
          Catch::Approx(-0.0085286852).epsilon(0.002));
    CHECK(initialStressTensor(0, 2) ==
          Catch::Approx(-0.10797027).epsilon(0.02));
    CHECK(initialStressTensor(1, 2) ==
          Catch::Approx(-0.035841973).epsilon(0.02));
  }
}
