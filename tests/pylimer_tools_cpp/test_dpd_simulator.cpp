#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/DPDSimulator.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pe = pylimer_tools::entities;
namespace ps = pylimer_tools::sim;
namespace pu = pylimer_tools::utils;
namespace pcd = pylimer_tools::sim::dpd;

/**
 * @brief Auxiliary helper function: sets different output types
 *
 * @param simulator
 * @param averageFile
 * @param autocorrFile
 */
void
setupAllOutputs(pcd::DPDSimulator& simulator,
                const std::string& averageFile,
                const std::string& autocorrFile)
{
  std::vector<ps::ComputedDoubleValues> outputQuantities = {
    ps::ComputedDoubleValues::TEMPERATURE, ps::ComputedDoubleValues::PRESSURE,
    ps::ComputedDoubleValues::STRESS_XX,   ps::ComputedDoubleValues::STRESS_YY,
    ps::ComputedDoubleValues::STRESS_ZZ,   ps::ComputedDoubleValues::MSD
  };

  ps::OutputConfiguration config;
  config.filename = "";
  config.outputEvery = 1;
  config.doubleValues = outputQuantities;
  config.intValues = { ps::ComputedIntValues::STEP };

  std::vector<ps::OutputConfiguration> configs = { config };
  REQUIRE_NOTHROW(simulator.configStepOutput(configs));

  std::vector<ps::ComputedDoubleValues> averageQuantities = {
    ps::ComputedDoubleValues::TEMPERATURE, ps::ComputedDoubleValues::PRESSURE,
    ps::ComputedDoubleValues::STRESS_XX,   ps::ComputedDoubleValues::STRESS_YY,
    ps::ComputedDoubleValues::STRESS_ZZ,   ps::ComputedDoubleValues::MSD
  };

  ps::OutputConfiguration avgconfig;
  avgconfig.outputEvery = 20;
  avgconfig.filename = averageFile;
  avgconfig.doubleValues = averageQuantities;

  std::vector<ps::OutputConfiguration> avgconfigs = { avgconfig };
  REQUIRE_NOTHROW(simulator.configAverageOutput(avgconfigs));

  std::vector<ps::ComputedDoubleValues> autocorrelationQuantities = {
    ps::ComputedDoubleValues::STRESS_XX,  ps::ComputedDoubleValues::STRESS_YY,
    ps::ComputedDoubleValues::STRESS_ZZ,  ps::ComputedDoubleValues::STRESS_XY,
    ps::ComputedDoubleValues::STRESS_YZ,  ps::ComputedDoubleValues::STRESS_XZ,
    ps::ComputedDoubleValues::STRESS_NXY, ps::ComputedDoubleValues::STRESS_NYZ,
    ps::ComputedDoubleValues::STRESS_NXZ,
  };
  ps::OutputConfiguration autocorrconfig;
  autocorrconfig.outputEvery = 25;
  autocorrconfig.filename = autocorrFile;
  autocorrconfig.doubleValues = autocorrelationQuantities;

  std::vector<ps::OutputConfiguration> autocorrconfigs = { autocorrconfig };
  REQUIRE_NOTHROW(simulator.configAutoCorrelatorOutput(autocorrconfigs));

  std::vector<size_t> atomIdsForMSD = { 1, 4, 6 };
  REQUIRE_NOTHROW(simulator.startMeasuringMSDForAtoms(atomIdsForMSD));
}

TEST_CASE("DPD Simulator Works", "[analysis][DPDSimulator][long]")
{
  std::cout << "Running test \"DPD Simulator Works\"" << std::endl;
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "/structure/melt_83_a_100.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  pe::UniverseSequence universeSequence = pe::UniverseSequence();
  std::vector<pu::AtomStyle> atomStyles = { pu::AtomStyle::HYBRID,
                                            pu::AtomStyle::BOND,
                                            pu::AtomStyle::EDPD };
  universeSequence.setDataFileAtomStyle(atomStyles);
  REQUIRE(universeSequence.getLength() == 0);
  universeSequence.initializeFromDataSequence({ { inputFile } });
  REQUIRE(universeSequence.getLength() == 1);
  pe::Universe universe = universeSequence.atIndex(0);

  pcd::DPDSimulator simulator =
    pcd::DPDSimulator(universe, 2, 9, false, "12th_seed");

  // configuration
  REQUIRE_NOTHROW(simulator.validateState());
  // CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
  // CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));
  REQUIRE_NOTHROW(simulator.configA(25.));
  REQUIRE_NOTHROW(simulator.configSigma(3.));
  REQUIRE_NOTHROW(simulator.configSpringConstant(2.));
  REQUIRE_NOTHROW(simulator.configSlipspringLowCutoff(0.5));
  REQUIRE_NOTHROW(simulator.configSlipspringHighCutoff(2.0));

  // check that wrong configuration throws
  REQUIRE_THROWS(simulator.configSlipspringLowCutoff(3.0));
  REQUIRE_THROWS(simulator.configSlipspringLowCutoff(2.0));
  REQUIRE_THROWS(simulator.configSlipspringHighCutoff(0.5));

  // initial state – vgl. lammps. Again, caution, randomness!
  CHECK_THAT(simulator.getStressTensor().trace() / 3.,
             Catch::Matchers::WithinAbs(23.321285, 0.35));
  CHECK_THAT(simulator.getTemperature() + 1e-2,
             Catch::Matchers::WithinAbs(0. + 1e-2, 1e-15));
  Eigen::Matrix3d initialStressTensor = simulator.getStressTensor();
  CHECK_THAT(initialStressTensor(0, 0),
             Catch::Matchers::WithinAbs(23.456778, 0.75));
  CHECK_THAT(initialStressTensor(1, 1),
             Catch::Matchers::WithinAbs(23.374175, 0.75));
  CHECK_THAT(initialStressTensor(2, 2),
             Catch::Matchers::WithinAbs(23.401583, 0.75));
  CHECK_THAT(initialStressTensor(0, 1),
             Catch::Matchers::WithinAbs(-0.0085286852, 0.2));
  CHECK_THAT(initialStressTensor(0, 2),
             Catch::Matchers::WithinAbs(-0.10797027, 0.2));
  CHECK_THAT(initialStressTensor(1, 2),
             Catch::Matchers::WithinAbs(-0.035841973, 0.2));

  std::string averageFile =
    suspectedPath + "melt_83_a_100.structure.avg-out.txt";
  std::string autocorrFile =
    suspectedPath + "melt_83_a_100.structure.autocorr-out.txt";
  setupAllOutputs(simulator, averageFile, autocorrFile);

  // restart files
  // std::string restartFile = suspectedPath + "dpd_restart_file.xml";
  std::string restartFile = suspectedPath + "dpd_restart_file.bin";
  simulator.configRestartOutput(restartFile, 30);

  // actual simulation – without slip-springs yet
  REQUIRE_NOTHROW(simulator.runSimulation(75, false));
  REQUIRE_NOTHROW(simulator.validateState());
  CHECK_THAT(simulator.getStressTensor().trace() / 3.,
             Catch::Matchers::WithinAbs(20.8, 1.));
  std::cout << "DPD ran, state validated." << std::endl;
  // CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
  // CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));

  CHECK_THAT(simulator.getTemperature(), Catch::Matchers::WithinAbs(1.0, 0.5));
  CHECK_NOTHROW(simulator.validateState());

  simulator.createSlipSprings(100, 2);
  // turn down the DPD & MC steps as we don't do as many total steps
  REQUIRE_NOTHROW(simulator.configNumStepsDPD(25));
  REQUIRE_NOTHROW(simulator.configNumStepsMC(25));
  CHECK_NOTHROW(simulator.validateState());
  std::cout << "DPD slip-springs created." << std::endl;

  pe::Universe resultUniverse = simulator.getUniverse();
  CHECK(resultUniverse.getNrOfBonds() == 100 + universe.getNrOfBonds());
  CHECK(resultUniverse.getNrOfAtoms() == universe.getNrOfAtoms());
  CHECK_NOTHROW(simulator.validateState());

  CHECK_NOTHROW(simulator.configShiftOneAtATime(false));
  CHECK_NOTHROW(simulator.configShiftPossibilityEmpty(true));
  CHECK_NOTHROW(simulator.runSimulation(26, true));
  CHECK_NOTHROW(simulator.configShiftPossibilityEmpty(false));
  CHECK_NOTHROW(simulator.runSimulation(26, true));

  CHECK_NOTHROW(simulator.validateState());
  std::cout << "DPD ran with slip-springs (both), state validated."
            << std::endl;

  CHECK_NOTHROW(simulator.configShiftOneAtATime(true));
  CHECK_NOTHROW(simulator.configShiftPossibilityEmpty(true));
  CHECK_NOTHROW(simulator.runSimulation(26, true));
  CHECK_NOTHROW(simulator.configShiftPossibilityEmpty(false));
  CHECK_NOTHROW(simulator.runSimulation(26, true));

  CHECK_NOTHROW(simulator.validateState());
  std::cout << "DPD ran with slip-springs (single), state validated."
            << std::endl;
  CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));

  CHECK(std::filesystem::exists(averageFile));
  std::remove(averageFile.c_str());
  CHECK(std::filesystem::exists(autocorrFile));
  std::remove(autocorrFile.c_str());

#ifdef CEREALIZABLE
  REQUIRE(std::filesystem::exists(restartFile));

  pcd::DPDSimulator sim2 = pcd::DPDSimulator::readRestartFile(restartFile);
  CHECK_NOTHROW(sim2.validateState());
  std::cout << "DPD read from restart file, state validated." << std::endl;
  REQUIRE_NOTHROW(sim2.runSimulation(5, false));
  CHECK(sim2.getTemperature() ==
        Catch::Approx(simulator.getTemperature()).margin(0.1));
  CHECK_NOTHROW(sim2.validateState());
  std::cout << "DPD ran from restart file, state validated." << std::endl;
  std::remove(restartFile.c_str());
#endif
};

#ifdef CEREALIZABLE
TEST_CASE("DPD Simulator Restart Files Work", "[analysis][DPDSimulator][long]")
{
  std::cout << "Running test \"DPD Simulator Restart Files Work\"" << std::endl;

  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "/structure/melt_83_a_100.structure.out";
  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    std::vector<pu::AtomStyle> atomStyles = { pu::AtomStyle::HYBRID,
                                              pu::AtomStyle::BOND,
                                              pu::AtomStyle::EDPD };
    universeSequence.setDataFileAtomStyle(atomStyles);
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "1st_seed");

    std::string averageFile =
      suspectedPath + "melt_83_a_100.structure.avg-out.txt";
    std::string autocorrFile =
      suspectedPath + "melt_83_a_100.structure.autocorr-out.txt";
    setupAllOutputs(simulator, averageFile, autocorrFile);

    simulator.createSlipSprings(100, 2);
    // turn down the DPD & MC steps as we don't do as many total steps
    REQUIRE_NOTHROW(simulator.configNumStepsDPD(25));
    REQUIRE_NOTHROW(simulator.configNumStepsMC(25));
    CHECK_NOTHROW(simulator.validateState());
    std::cout << "DPD slip-springs created." << std::endl;

    CHECK_NOTHROW(simulator.runSimulation(5, true));
    CHECK_NOTHROW(simulator.validateState());

    std::string restartFile = suspectedPath + "dpd_restart_file.bin";

    simulator.writeRestartFile(restartFile);

    CHECK_NOTHROW(simulator.validateState());

    pcd::DPDSimulator simulator2 =
      pcd::DPDSimulator::readRestartFile(restartFile);

    CHECK_NOTHROW(simulator2.validateState());
    CHECK_NOTHROW(simulator2.runSimulation(5, true));
    CHECK_NOTHROW(simulator2.validateState());
    std::cout << "DPD restart file is read appropriately." << std::endl;

    CHECK(std::filesystem::exists(averageFile));
    std::remove(averageFile.c_str());
    CHECK(std::filesystem::exists(autocorrFile));
    std::remove(autocorrFile.c_str());
    REQUIRE(std::filesystem::exists(restartFile));
  }
}
#endif

TEST_CASE("DPD Simulator Can Crosslink", "[analysis][DPDSimulator][long]")
{
  std::cout << "Running test \"DPD Simulator Can Crosslink\"" << std::endl;
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "/structure/melt_213_a_47_106_xlinks_v_1.structure.out";
  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    std::vector<pu::AtomStyle> atomStyles = { pu::AtomStyle::ANGLE };
    universeSequence.setDataFileAtomStyle(atomStyles);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "seed2");

    // configuration
    REQUIRE_NOTHROW(simulator.validateState());
    REQUIRE_NOTHROW(simulator.configA(25.));
    REQUIRE_NOTHROW(simulator.configSigma(3.));
    // turn down the DPD & MC steps as we don't do as many total steps
    REQUIRE_NOTHROW(simulator.configNumStepsDPD(50));
    REQUIRE_NOTHROW(simulator.configNumStepsMC(50));
    REQUIRE_NOTHROW(simulator.configSlipspringLowCutoff(0.5));
    REQUIRE_NOTHROW(simulator.configSlipspringHighCutoff(2.0));
    REQUIRE_THROWS(simulator.configSlipspringLowCutoff(3.0));
    REQUIRE_THROWS(simulator.configSlipspringLowCutoff(2.0));
    REQUIRE_THROWS(simulator.configSlipspringHighCutoff(0.5));

    simulator.createSlipSprings(100, 2);
    // turn down the DPD & MC steps as we don't do as many total steps
    REQUIRE_NOTHROW(simulator.configNumStepsDPD(25));
    REQUIRE_NOTHROW(simulator.configNumStepsMC(25));

    std::vector<ps::ComputedDoubleValues> outputQuantities = {
      ps::ComputedDoubleValues::TEMPERATURE,
      ps::ComputedDoubleValues::PRESSURE,
      ps::ComputedDoubleValues::VOLUME,
      ps::ComputedDoubleValues::STRESS_XX,
      ps::ComputedDoubleValues::STRESS_YY,
      ps::ComputedDoubleValues::STRESS_ZZ,
      ps::ComputedDoubleValues::MSD
    };

    ps::OutputConfiguration config;
    config.filename = "";
    config.outputEvery = 5;
    config.doubleValues = outputQuantities;
    config.intValues = { ps::ComputedIntValues::STEP };

    std::vector<ps::OutputConfiguration> configs = { config };
    REQUIRE_NOTHROW(simulator.configStepOutput(configs));

    std::unordered_map<int, int> numBondsPerType;
    numBondsPerType[1] = 2;
    numBondsPerType[2] = 4;
    REQUIRE_NOTHROW(
      simulator.configBondFormation(250, numBondsPerType, 2.2, 5));
    size_t numBondsBefore = simulator.getNumBonds();

    // actual simulation
    REQUIRE_NOTHROW(simulator.runSimulation(250, false));
    REQUIRE_NOTHROW(simulator.validateState());
    std::cout << "DPD ran, state validated." << std::endl;

    // check that we actually have formed bonds
    CHECK(simulator.getNumBonds() >= (25 + numBondsBefore));
    pe::Universe universeAfter = simulator.getUniverse(false);
    std::map<int, int> finalFunctionalityPerType =
      universeAfter.determineFunctionalityPerType();
    CHECK(finalFunctionalityPerType.at(1) == 2.);
    CHECK(finalFunctionalityPerType.at(2) <= 4.);
  }
};

TEST_CASE("DPD Simulator Computes Correct Forces",
          "[analysis][DPDSimulator][long]")
{
  std::cout << "Running test \"DPD Simulator Computes Correct Forces\""
            << std::endl;
  // note that the random force might lead to deviations compared to LAMMPS
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "/structure/melt_83_a_100.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  pe::UniverseSequence universeSequence = pe::UniverseSequence();
  REQUIRE(universeSequence.getLength() == 0);
  std::vector<pu::AtomStyle> atomStyles = { pu::AtomStyle::HYBRID,
                                            pu::AtomStyle::BOND,
                                            pu::AtomStyle::EDPD };
  universeSequence.setDataFileAtomStyle(atomStyles);
  universeSequence.initializeFromDataSequence({ { inputFile } });
  REQUIRE(universeSequence.getLength() == 1);
  pe::Universe universe = universeSequence.atIndex(0);

  pcd::DPDSimulator simulator =
    pcd::DPDSimulator(universe, 2, 9, false, "14th_seed");

  // configuration
  REQUIRE_NOTHROW(simulator.validateState());
  CHECK(simulator.getVolume() == Catch::Approx(2766.6667));
  // CHECK_NOTHROW(simulator.validateNeighbourlist(2.0));
  // CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));

  // DISABLE all forces by configuration
  REQUIRE_NOTHROW(simulator.configA(0.));
  REQUIRE_NOTHROW(simulator.configSigma(0.));
  REQUIRE_NOTHROW(simulator.configSpringConstant(0.));

  REQUIRE_NOTHROW(simulator.refreshCurrentState());
  CHECK(simulator.getStressTensor().trace() / 3. == Catch::Approx(0.));
  CHECK(simulator.getStressTensor()(1, 2) == Catch::Approx(0.));
  CHECK(simulator.getStressTensor()(1, 1) == Catch::Approx(0.));

  // gradually add forces
  REQUIRE_NOTHROW(simulator.configA(25.));
  REQUIRE_NOTHROW(simulator.refreshCurrentState());
  CHECK(simulator.getStressTensor().trace() / 3. == Catch::Approx(25.259583));
  CHECK(simulator.getStressTensor()(0, 0) == Catch::Approx(25.394731));
  CHECK(simulator.getStressTensor()(1, 1) == Catch::Approx(25.249986));
  CHECK(simulator.getStressTensor()(1, 2) == Catch::Approx(-0.033026848));

  REQUIRE_NOTHROW(simulator.configA(0.));
  REQUIRE_NOTHROW(simulator.configSpringConstant(2.));
  REQUIRE_NOTHROW(simulator.refreshCurrentState());
  CHECK(simulator.getStressTensor().trace() / 3. == Catch::Approx(-1.78695));
  CHECK(simulator.getStressTensor()(0, 0) == Catch::Approx(-1.8259387));
  CHECK(simulator.getStressTensor()(1, 1) == Catch::Approx(-1.7699664));
  CHECK(simulator.getStressTensor()(1, 2) == Catch::Approx(0.021150486));

  // CAUTION, randomness
  REQUIRE_NOTHROW(simulator.configSigma(3.));
  REQUIRE_NOTHROW(simulator.refreshCurrentState());
  CHECK_THAT(simulator.getStressTensor().trace() / 3.,
             Catch::Matchers::WithinAbs(-1.8487375, 0.5));

  // complete initial state
  REQUIRE_NOTHROW(simulator.configA(25.));
  REQUIRE_NOTHROW(simulator.configSigma(3.));
  REQUIRE_NOTHROW(simulator.configSpringConstant(2.));
  REQUIRE_NOTHROW(simulator.refreshCurrentState());
  CHECK_THAT(simulator.getStressTensor().trace() / 3.,
             Catch::Matchers::WithinAbs(23.321285, 0.5));
  CHECK(simulator.getTemperature() + 1e-2 == Catch::Approx(0. + 1e-2));
  Eigen::Matrix3d initialStressTensor = simulator.getStressTensor();
  CHECK_THAT(initialStressTensor(0, 0),
             Catch::Matchers::WithinAbs(23.456778, 0.75));
  CHECK_THAT(initialStressTensor(1, 1),
             Catch::Matchers::WithinAbs(23.374175, 0.75));
  CHECK_THAT(initialStressTensor(2, 2),
             Catch::Matchers::WithinAbs(23.401583, 0.75));
  CHECK_THAT(initialStressTensor(0, 1),
             Catch::Matchers::WithinAbs(-0.0085286852, 0.15));
  CHECK_THAT(initialStressTensor(0, 2),
             Catch::Matchers::WithinAbs(-0.10797027, 0.15));
  CHECK_THAT(initialStressTensor(1, 2),
             Catch::Matchers::WithinAbs(-0.035841973, 0.15));
  REQUIRE_NOTHROW(simulator.validateState());

  // short run because we can

  std::string averageFile =
    suspectedPath + "melt_83_a_100.structure.avg-out.txt";
  std::string autocorrFile =
    suspectedPath + "melt_83_a_100.structure.autocorr-out.txt";
  setupAllOutputs(simulator, averageFile, autocorrFile);
  REQUIRE_NOTHROW(simulator.runSimulation(5, false));
  REQUIRE_NOTHROW(simulator.runSimulation(5, true));

  CHECK(std::filesystem::exists(averageFile));
  std::remove(averageFile.c_str());
  CHECK(std::filesystem::exists(autocorrFile));
  std::remove(autocorrFile.c_str());
}

TEST_CASE("DPD Simulator Converts Correctly", "[analysis][DPDSimulator]")
{
  std::cout << "Running test \"DPD Simulator Converts Correctly\"" << std::endl;
  // note that the random force might lead to deviations compared to LAMMPS
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "/structure/melt_83_a_100.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  pe::UniverseSequence universeSequence = pe::UniverseSequence();
  REQUIRE(universeSequence.getLength() == 0);
  universeSequence.initializeFromDataSequence({ { inputFile } });
  REQUIRE(universeSequence.getLength() == 1);
  pe::Universe universe = universeSequence.atIndex(0);

  pcd::DPDSimulator simulator =
    pcd::DPDSimulator(universe, 2, 9, false, "12th_seed");

  pe::Universe resultUniverse = simulator.getUniverse();

  CHECK(resultUniverse.getNrOfAtoms() == universe.getNrOfAtoms());
  CHECK(resultUniverse.getNrOfBonds() == universe.getNrOfBonds());
  std::map<std::string, std::vector<long int>> previousEdges =
    universe.getEdges();
  std::map<std::string, std::vector<long int>> newEdges =
    resultUniverse.getEdges();

  for (size_t i = 0; i < resultUniverse.getNrOfBonds(); ++i) {
    CHECK(previousEdges["edge_from"][i] == newEdges["edge_from"][i]);
    CHECK(previousEdges["edge_to"][i] == newEdges["edge_to"][i]);
    CHECK(previousEdges["edge_type"][i] == newEdges["edge_type"][i]);
  }

  CHECK_NOTHROW(simulator.validateNeighbourlist(1.0));
  CHECK_NOTHROW(simulator.attemptBondFormation());
  CHECK_NOTHROW(simulator.relocateSlipSprings());
}

TEST_CASE("DPD can deform box", "[analysis][DPDSimulator][long]")
{
  std::cout << "Running test \"DPD can deform box\"" << std::endl;

  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "/structure/"
                    "crosslinked_p_0.98_melt_100_a_38_50_xlinks_v_22.structure."
                    "out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "25th_seed");

    std::vector<ps::ComputedDoubleValues> outputQuantities = {
      ps::ComputedDoubleValues::TEMPERATURE,
      ps::ComputedDoubleValues::PRESSURE,
      ps::ComputedDoubleValues::VOLUME,
      ps::ComputedDoubleValues::STRESS_XX,
      ps::ComputedDoubleValues::STRESS_YY,
      ps::ComputedDoubleValues::STRESS_ZZ,
      ps::ComputedDoubleValues::MSD
    };

    ps::OutputConfiguration config;
    config.filename = "";
    config.outputEvery = 5;
    config.doubleValues = outputQuantities;
    config.intValues = { ps::ComputedIntValues::STEP };

    std::vector<ps::OutputConfiguration> configs = { config };
    REQUIRE_NOTHROW(simulator.configStepOutput(configs));

    REQUIRE_NOTHROW(simulator.validateState());
    REQUIRE_NOTHROW(simulator.createSlipSprings(100, 2));
    REQUIRE_NOTHROW(simulator.validateState());
    REQUIRE_NOTHROW(simulator.runSimulation(10));
    REQUIRE_NOTHROW(simulator.validateState());

    pe::Box secondBox = pe::Box(1.1 * universe.getBox().getLowX(),
                                1.1 * universe.getBox().getHighX(),
                                universe.getBox().getLowY(),
                                universe.getBox().getHighY(),
                                (1. / 1.1) * universe.getBox().getLowZ(),
                                (1. / 1.1) * universe.getBox().getHighZ());

    CHECK(secondBox.getVolume() == Catch::Approx(simulator.getVolume()));

    SECTION("Deform slowly")
    {
      simulator.configBoxDeformation(secondBox);

      REQUIRE_NOTHROW(simulator.runSimulation(200));
      REQUIRE_NOTHROW(simulator.validateState());

      REQUIRE_NOTHROW(simulator.runSimulation(50));
      REQUIRE_NOTHROW(simulator.validateState());
    }

    SECTION("Deform immediately")
    {
      simulator.deformBoxImmediately(secondBox);

      REQUIRE_NOTHROW(simulator.runSimulation(200));
      REQUIRE_NOTHROW(simulator.validateState());
    }
  }
}

TEST_CASE("New PBC computation is correct", "[analysis][DPDSimulator][1proc]")
{
  std::cout << "Running test \"New PBC computation is correct\"" << std::endl;
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "/structure/"
                    "crosslinked_p_0.98_melt_100_a_3_50_xlinks_v_14.converted."
                    "structure.out-equilibration_do_crosslink.structure.out";
  if (std::filesystem::exists(inputFile)) {
#ifdef OPENMP_FOUND
    // we cannot have more than 1 thread, otherwise the random number generator
    // will not play nicely.
    omp_set_num_threads(1);
#endif

    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "15th_seed");

#ifdef OPENMP_FOUND
    simulator.setNumOmpThreads(1);
#endif
    REQUIRE(simulator.getNumOmpThreads() == 1);
    simulator.createSlipSprings(100, 2);

    // invoke copy-constructor
    pcd::DPDSimulator simulator2 = simulator;
#ifdef OPENMP_FOUND
    simulator2.setNumOmpThreads(1);
#endif
    REQUIRE(simulator.getNumOmpThreads() == 1);

    // switch to "common" PBC
    simulator2.configAssumeBoxLargeEnough();

    simulator.reseedRandomness("15th_seed");
    simulator.refreshCurrentState();
    simulator2.reseedRandomness("15th_seed");
    simulator2.refreshCurrentState();

    CHECK(simulator.getUniformRandBetween0And1() ==
          simulator2.getUniformRandBetween0And1());
    CHECK_THAT(simulator.getTemperature(),
               Catch::Matchers::WithinRel(simulator2.getTemperature()));
    CHECK(simulator.getBondLengths().isApprox(simulator2.getBondLengths()));
    CHECK(simulator.getStressTensor().isApprox(simulator2.getStressTensor()));
    // std::cout << simulator.getStressTensor() << std::endl;
  }
}

TEST_CASE("For large systems the PBC method does not matter",
          "[analysis][DPDSimulator][1proc][long]")
{
  /**
   * @brief Currently, this test fails.
   *
   * This might be ok,
   * demonstrating that the PBC method *does* matter.
   */
  std::cout
    << "Running test \"For large systems the PBC method does not matter\""
    << std::endl;
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile = suspectedPath +
                          "/structure/"
                          "3d-diamond-lattice_3x3x3_a_23_d_3_v_0.structure.out";

#ifdef OPENMP_FOUND
  // we cannot have more than 1 thread, otherwise the random number generator
  // will not play nicely.
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif

  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "19th_seed");
#ifdef OPENMP_FOUND
    simulator.setNumOmpThreads(1);
#endif

    REQUIRE(simulator.getNumOmpThreads() == 1);
    simulator.createSlipSprings(static_cast<int>(0.1 * universe.getNrOfAtoms()),
                                2);

    std::vector<ps::ComputedDoubleValues> outputQuantities = {
      ps::ComputedDoubleValues::TEMPERATURE,
      ps::ComputedDoubleValues::PRESSURE,
      ps::ComputedDoubleValues::MAX_B,
      ps::ComputedDoubleValues::MEAN_B
    };

    ps::OutputConfiguration config;
    config.filename = "";
    config.outputEvery = 1;
    config.doubleValues = outputQuantities;
    config.intValues = { ps::ComputedIntValues::STEP };

    std::vector<ps::OutputConfiguration> configs = { config };
    REQUIRE_NOTHROW(simulator.configStepOutput(configs));

    // invoke copy-constructor
    pcd::DPDSimulator simulator2 = pcd::DPDSimulator(simulator);
#ifdef OPENMP_FOUND
    simulator.setNumOmpThreads(1);
    simulator2.setNumOmpThreads(1);
#endif

    REQUIRE(simulator.getNumOmpThreads() == 1);
    REQUIRE(simulator2.getNumOmpThreads() == 1);
    // switch to "common" PBC
    simulator2.configAssumeBoxLargeEnough();

    REQUIRE(!simulator.assumesBoxLargeEnough());
    REQUIRE(simulator2.assumesBoxLargeEnough());

    simulator.reseedRandomness("20th_seed");
    simulator.refreshCurrentState();
    simulator2.reseedRandomness("20th_seed");
    simulator2.refreshCurrentState();

    REQUIRE(simulator.getUniformRandBetween0And1() ==
            simulator2.getUniformRandBetween0And1());
    CHECK_THAT(simulator.getTemperature(),
               Catch::Matchers::WithinRel(simulator2.getTemperature()));
    REQUIRE(simulator.getBondLengths().isApprox(simulator2.getBondLengths()));
    REQUIRE(simulator.getStressTensor().isApprox(simulator2.getStressTensor()));
    REQUIRE(simulator.getUniformRandBetween0And1() ==
            simulator2.getUniformRandBetween0And1());

    // not sure what value is sensible here...
    simulator.runSimulation(50, true);
    simulator2.runSimulation(50, true);

    REQUIRE(simulator.getNumOmpThreads() == 1);
    REQUIRE(simulator2.getNumOmpThreads() == 1);

    CHECK(simulator.getUniformRandBetween0And1() ==
          simulator2.getUniformRandBetween0And1());
    CHECK_THAT(simulator.getTemperature(),
               Catch::Matchers::WithinRel(simulator2.getTemperature()));
    CHECK(simulator.getBondLengths().isApprox(simulator2.getBondLengths()));
    CHECK(simulator.getStressTensor().isApprox(simulator2.getStressTensor()));
  }
}

TEST_CASE("DPD Simulator's restart files are accurate",
          "[analysis][DPDSimulator][1proc][long]")
{
  // note that the random force might lead to deviations compared to LAMMPS
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile =
    suspectedPath + "/structure/melt_83_a_100.structure.out";
  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "14th_seed");

    simulator.createSlipSprings(100, 2);

    std::string restartFile = "restartFile-for-accuracy-test.bin";
#ifdef CEREALIZABLE
    simulator.writeRestartFile(restartFile);

    REQUIRE(std::filesystem::exists(restartFile));

    pcd::DPDSimulator sim2 = pcd::DPDSimulator::readRestartFile(restartFile);

    std::remove(restartFile.c_str());
#else
    pcd::DPDSimulator sim2 = simulator;
#endif

#ifdef OPENMP_FOUND
    // we cannot have more than 1 thread, otherwise the random number generator
    // will not play nicely.
    omp_set_num_threads(1);
    simulator.setNumOmpThreads(1);
    sim2.setNumOmpThreads(1);
#endif

    CHECK(simulator.getCoordinates().isApprox(sim2.getCoordinates()));
    CHECK(simulator.getBondLengths().isApprox(sim2.getBondLengths()));
    CHECK(simulator.getStressTensor().isApprox(sim2.getStressTensor()));
    CHECK(sim2.getTemperature() == simulator.getTemperature());
    CHECK_NOTHROW(sim2.validateState());

    for (size_t i = 0; i < 1000; ++i) {
      CHECK(sim2.getUniformRandBetween0And1() ==
            simulator.getUniformRandBetween0And1());
      CHECK(sim2.getUniformRandMean0Std1() ==
            simulator.getUniformRandMean0Std1());
    }

    simulator.refreshCurrentState();
    sim2.refreshCurrentState();

    CHECK(simulator.getCoordinates().isApprox(sim2.getCoordinates()));
    CHECK(simulator.getBondLengths().isApprox(sim2.getBondLengths()));
    CHECK(simulator.getStressTensor().isApprox(sim2.getStressTensor()));
    CHECK(sim2.getTemperature() == simulator.getTemperature());

    for (size_t i = 0; i < 10; ++i) {
      CHECK(sim2.getUniformRandBetween0And1() ==
            simulator.getUniformRandBetween0And1());
      CHECK(sim2.getUniformRandMean0Std1() ==
            simulator.getUniformRandMean0Std1());
    }

    simulator.runSimulation(500, false);
    sim2.runSimulation(500, false);

    CHECK(simulator.getCoordinates().isApprox(sim2.getCoordinates()));
    CHECK(simulator.getBondLengths().isApprox(sim2.getBondLengths()));
    CHECK(sim2.getTemperature() == simulator.getTemperature());
  }
}
#ifdef CEREALIZABLE
TEST_CASE("DPD Simulator can be serialized",
          "[DPDSimulator][serialization][1proc][long]")
{
  // note that the random force might lead to deviations compared to LAMMPS
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  for (std::string inputFile :
       { // suspectedPath +
         // "/structure/melt_213_a_47_106_xlinks_v_1.structure.out",
         suspectedPath + "/structure/network_100_a_46.structure.out" }) {
    INFO("Using structure file " << inputFile);
    REQUIRE(std::filesystem::exists(inputFile));
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);
    universe.resampleVelocities(0.0, 2.0, "abcSEED", false);
    CHECK(universe.vertexPropertyExists("vx"));
    CHECK(universe.vertexPropertyExists("vy"));
    CHECK(universe.vertexPropertyExists("vz"));

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "14th_seed");
    CHECK_FALSE(simulator.getVelocities().isZero());

    double T = simulator.getTemperature();
    CHECK(T > 0.0);

    CHECK_THROWS(simulator.configSlipspringHighCutoff(-1));
    CHECK_NOTHROW(simulator.configSlipspringHighCutoff(2.));
    CHECK_THROWS(simulator.configSlipspringLowCutoff(3.));
    CHECK_NOTHROW(simulator.configSlipspringLowCutoff(0.5));
    simulator.createSlipSprings(100, 2);

    simulator.configNumStepsMC(50);
    simulator.configNumStepsDPD(50);
    simulator.configShiftPossibilityEmpty(true);
    simulator.configAssumeBoxLargeEnough();
    CHECK(simulator.getCurrentTime(0.) == 0.);
    CHECK(simulator.getCurrentTime(1000.) == 0.);

    std::string restartFile = "restartFile-for-accuracy-test.bin";
    simulator.writeRestartFile(restartFile);

    REQUIRE(std::filesystem::exists(restartFile));

    pcd::DPDSimulator sim2 = pcd::DPDSimulator::readRestartFile(restartFile);
    std::remove(restartFile.c_str());

    // on 1 proc, expect equal randomness state
    CHECK(simulator.getUniformRandBetween0And1() ==
          sim2.getUniformRandBetween0And1());

    // and equal everything else
    CHECK(simulator.getCoordinates().isApprox(sim2.getCoordinates()));
    CHECK(simulator.getBondLengths().isApprox(sim2.getBondLengths()));
    CHECK(simulator.getStressTensor().isApprox(sim2.getStressTensor()));
    CHECK(sim2.getTemperature() == simulator.getTemperature());
    CHECK(simulator.computeBondLength(2) == sim2.computeBondLength(2));

    // when reseeding the randomness, the state should be different
    CHECK_NOTHROW(sim2.reseedRandomness(""));
    CHECK_NOTHROW(simulator.reseedRandomness(""));
    CHECK(simulator.getUniformRandBetween0And1() !=
          sim2.getUniformRandBetween0And1());

    // finally, some other methods should work as well on the clone,
    // without issues
    CHECK_NOTHROW(sim2.validateNeighbourlist(1.));
    std::unordered_map<int, int> numBondsPerType = {};
    numBondsPerType[0] = 2;
    numBondsPerType[1] = 2;
    numBondsPerType[2] = 4;
    numBondsPerType[3] = 2;
    sim2.configBondFormation(3, numBondsPerType);
    CHECK_NOTHROW(sim2.attemptBondFormation());
    CHECK_NOTHROW(sim2.relocateSlipSprings());
    pe::Box newBox = pe::Box(10, 12, 10);
    CHECK_NOTHROW(sim2.deformBoxImmediately(newBox));
  }
}
#endif

TEST_CASE("OutputSupportingSimulation handles multiple autocorrelation configs "
          "and useEvery > 1",
          "[analysis][DPDSimulator][Correlator]")
{
  pe::Universe u(5.0, 5.0, 5.0);
  u.addAtoms({ 1, 2, 3, 4 },
             { 1, 1, 1, 1 },
             { 0.5, 1.5, 2.5, 3.5 },
             { 1.0, 1.0, 1.0, 1.0 },
             { 1.0, 1.0, 1.0, 1.0 },
             { 0, 0, 0, 0 },
             { 0, 0, 0, 0 },
             { 0, 0, 0, 0 });
  u.addBonds({ 1, 2, 3 }, { 2, 3, 4 });

  pcd::DPDSimulator simulator(u, 2, 9, false, "autocorr_test");

  const std::filesystem::path tmpDir = std::filesystem::temp_directory_path();
  const std::string file1 = (tmpDir / "pylimer_test_autocorr_1.txt").string();
  const std::string file2 = (tmpDir / "pylimer_test_autocorr_2.txt").string();
  std::remove(file1.c_str());
  std::remove(file2.c_str());

  struct TempFileCleanup
  {
    std::string path1;
    std::string path2;
    ~TempFileCleanup()
    {
      std::remove(path1.c_str());
      std::remove(path2.c_str());
    }
  };
  const TempFileCleanup cleanup{ file1, file2 };

  ps::OutputConfiguration ac1;
  ac1.filename = file1;
  ac1.useEvery = 2;
  ac1.outputEvery = 4;
  ac1.doubleValues = { ps::ComputedDoubleValues::STRESS_XX };

  ps::OutputConfiguration ac2;
  ac2.filename = file2;
  ac2.useEvery = 2;
  ac2.outputEvery = 4;
  ac2.doubleValues = { ps::ComputedDoubleValues::PRESSURE };

  REQUIRE_NOTHROW(simulator.configAutoCorrelatorOutput({ ac1, ac2 }, 4, 8, 2));
  REQUIRE_NOTHROW(simulator.runSimulation(4, false));

  REQUIRE(std::filesystem::exists(file1));
  REQUIRE(std::filesystem::exists(file2));

  auto parseDataRows =
    [](const std::string& path) -> std::vector<std::pair<double, double>> {
    std::ifstream in(path);
    std::string line;
    std::vector<std::pair<double, double>> rows;
    while (std::getline(in, line)) {
      if (!line.empty() && line[0] != '#' && line.rfind("Step\t", 0) != 0) {
        const size_t tabPos = line.find('\t');
        REQUIRE(tabPos != std::string::npos);
        const double lag = std::stod(line.substr(0, tabPos));
        const double val = std::stod(line.substr(tabPos + 1));
        rows.emplace_back(lag, val);
      }
    }
    return rows;
  };

  // 4 steps with useEvery = 2 gives 2 samples (steps 2 and 4) -> 2 lag points
  // (0, 1)
  const std::vector<std::pair<double, double>> rows1 = parseDataRows(file1);
  const std::vector<std::pair<double, double>> rows2 = parseDataRows(file2);
  REQUIRE(rows1.size() == 2);
  REQUIRE(rows2.size() == 2);

  for (size_t i = 0; i < 2; ++i) {
    CHECK(std::isfinite(rows1[i].first));
    CHECK(std::isfinite(rows1[i].second));
    CHECK(std::isfinite(rows2[i].first));
    CHECK(std::isfinite(rows2[i].second));
  }

  CHECK(rows1[0].first >= 0.0);
  CHECK(rows1[1].first > rows1[0].first);
  CHECK(rows2[0].first >= 0.0);
  CHECK(rows2[1].first > rows2[0].first);

  const bool distinctValues =
    (std::abs(rows1[0].second - rows2[0].second) > 1e-12) ||
    (std::abs(rows1[1].second - rows2[1].second) > 1e-12);
  CHECK(distinctValues);
}

TEST_CASE(
  "DPD Simulator features, outputs, slipsprings, and serialization work",
  "[analysis][DPDSimulator][serialization]")
{
  // Build a small 2-chain + 1-crosslinker system in a 5x5x5 box
  pe::Universe u(5.0, 5.0, 5.0);
  u.addAtoms(
    { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 },
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2 },
    { 0.5, 1.2, 1.9, 2.6, 3.3, 4.0, 0.5, 1.2, 1.9, 2.6, 3.3, 4.0, 4.2 },
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 1.1 },
    { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
  u.addBonds({ 1, 2, 3, 4, 5, 7, 8, 9, 10, 11 },
             { 2, 3, 4, 5, 6, 8, 9, 10, 11, 12 });
  u.resampleVelocities(0.0, 1.0, "fast_dpd_seed", false);

  pcd::DPDSimulator sim(u, 2, 9, false, "fast_dpd_sim");
  CHECK(sim.getSlipSpringBondType() == 9);
  CHECK(sim.getNumStepsMC() >= 0);
  CHECK(sim.getNumStepsDPD() >= 0);

  sim.configNumStepsMC(2);
  CHECK(sim.getNumStepsMC() == 2);
  sim.configNumStepsDPD(2);
  CHECK(sim.getNumStepsDPD() == 2);
  sim.configShiftPossibilityEmpty(true);
  CHECK(sim.getShiftPossibilityEmpty());
  sim.configShiftOneAtATime(true);
  CHECK(sim.getShiftOneAtATime());

  CHECK_THROWS(sim.configSlipspringHighCutoff(-1.0));
  CHECK_NOTHROW(sim.configSlipspringHighCutoff(2.0));
  CHECK_THROWS(sim.configSlipspringLowCutoff(3.0));
  CHECK_NOTHROW(sim.configSlipspringLowCutoff(0.3));
  CHECK_NOTHROW(sim.configAssumeBoxLargeEnough());

  CHECK(sim.getVolume() == Catch::Approx(125.0));
  CHECK(sim.getTemperature() > 0.0);
  CHECK_FALSE(sim.getVelocities().isZero());
  CHECK(sim.getBondLengths().size() == 10);
  CHECK(sim.computeBondLength(0) > 0.0);
  CHECK(sim.getStressTensor().rows() == 3);
  CHECK_NOTHROW(sim.validateState());
  CHECK_NOTHROW(sim.validateDebugState());
  CHECK_NOTHROW(sim.validateNeighbourlist(1.0));

  // Slip-springs & MSD
  CHECK_NOTHROW(sim.createSlipSprings(4, 9));
  CHECK_NOTHROW(sim.startMeasuringMSDForAtoms({ 1, 2, 3, 4 }));

  // Configure step, average, and restart outputs
  const std::filesystem::path tmpDir = std::filesystem::temp_directory_path();
  const std::string stepFile = (tmpDir / "pylimer_fast_dpd_step.txt").string();
  const std::string avgFile = (tmpDir / "pylimer_fast_dpd_avg.txt").string();
  const std::string restartFile =
    (tmpDir / "pylimer_fast_dpd_restart.bin").string();
  std::remove(stepFile.c_str());
  std::remove(avgFile.c_str());
  std::remove(restartFile.c_str());

  struct Cleanup
  {
    std::vector<std::string> paths;
    ~Cleanup()
    {
      for (const std::string& p : paths) {
        std::remove(p.c_str());
      }
    }
  } cleanup{ { stepFile, avgFile, restartFile } };

  ps::OutputConfiguration stepCfg;
  stepCfg.filename = stepFile;
  stepCfg.outputEvery = 2;
  stepCfg.intValues = {
    ps::ComputedIntValues::STEP,
    ps::ComputedIntValues::NUM_SHIFT,
    ps::ComputedIntValues::NUM_RELOC,
    ps::ComputedIntValues::NUM_ATOMS,
    ps::ComputedIntValues::NUM_EXTRA_ATOMS,
    ps::ComputedIntValues::NUM_BONDS,
    ps::ComputedIntValues::NUM_EXTRA_BONDS,
    ps::ComputedIntValues::NUM_BONDS_TO_FORM,
  };
  stepCfg.doubleValues = {
    ps::ComputedDoubleValues::TIME,        ps::ComputedDoubleValues::VOLUME,
    ps::ComputedDoubleValues::TEMPERATURE, ps::ComputedDoubleValues::PRESSURE,
    ps::ComputedDoubleValues::STRESS_XX,   ps::ComputedDoubleValues::GAMMA,
    ps::ComputedDoubleValues::RESIDUAL,    ps::ComputedDoubleValues::MEAN_B,
    ps::ComputedDoubleValues::MAX_B,       ps::ComputedDoubleValues::MSD,
  };
  std::vector<ps::OutputConfiguration> stepCfgs = { stepCfg };
  REQUIRE_NOTHROW(sim.configStepOutput(stepCfgs));

  ps::OutputConfiguration avgCfg;
  avgCfg.filename = avgFile;
  avgCfg.useEvery = 2;
  avgCfg.outputEvery = 4;
  avgCfg.intValues = { ps::ComputedIntValues::STEP };
  avgCfg.doubleValues = { ps::ComputedDoubleValues::TEMPERATURE,
                          ps::ComputedDoubleValues::MSD };
  REQUIRE_NOTHROW(sim.configAverageOutput({ avgCfg }));

#ifdef CEREALIZABLE
  sim.configRestartOutput(restartFile, 4);
#endif

  // Run 4 steps with MC enabled (exercises DPD integration, slipspring
  // relocation/shifting, step/average/restart output)
  REQUIRE_NOTHROW(sim.runSimulation(4, true));
  CHECK(std::filesystem::exists(stepFile));
  CHECK(std::filesystem::exists(avgFile));

  // Box deformation (immediate & gradual)
  CHECK_NOTHROW(sim.deformBoxImmediately(pe::Box(5.5, 5.5, 5.5)));
  CHECK(sim.getVolume() == Catch::Approx(5.5 * 5.5 * 5.5));
  sim.configBoxDeformation(pe::Box(6.0, 6.0, 6.0));
  CHECK_NOTHROW(sim.runSimulation(2, false));

  // Bond formation
  std::unordered_map<int, int> maxBondsPerType = { { 1, 2 }, { 2, 4 } };
  CHECK_NOTHROW(sim.configBondFormation(1, maxBondsPerType, 1.5, 1, 2, 1));
  CHECK(sim.getNrOfBondsToForm() == 1);
  CHECK_NOTHROW(sim.attemptBondFormation());

  // Reseed & serialization
  CHECK_NOTHROW(sim.reseedRandomness("another_seed"));
#ifdef CEREALIZABLE
  REQUIRE(std::filesystem::exists(restartFile));
  pcd::DPDSimulator sim2 = pcd::DPDSimulator::readRestartFile(restartFile);
  CHECK(sim2.getCoordinates().isApprox(sim2.getCoordinates()));
  CHECK(sim2.getVolume() > 0.0);
#endif
}

TEST_CASE("DPD Simulator slip-spring MC physics and edge cases",
          "[analysis][DPDSimulator]")
{
  SECTION("Defect-free network relocateSlipSprings does not underflow or crash")
  {
    // Closed ring of 6 beads: all atoms have functionality 2 -> chainEndIndices
    // is empty when allowRelocationInNetwork == false.
    pe::Universe u(10.0, 10.0, 10.0);
    u.addAtoms({ 1, 2, 3, 4, 5, 6 },
               { 1, 1, 1, 1, 1, 1 },
               { 1.0, 2.0, 3.0, 3.0, 2.0, 1.0 },
               { 1.0, 1.0, 1.0, 2.0, 2.0, 2.0 },
               { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },
               { 0, 0, 0, 0, 0, 0 },
               { 0, 0, 0, 0, 0, 0 },
               { 0, 0, 0, 0, 0, 0 });
    u.addBonds({ 1, 2, 3, 4, 5, 6 }, { 2, 3, 4, 5, 6, 1 });

    pcd::DPDSimulator sim(u, 2, 9, false, "defect_free_seed");
    sim.configSlipspringLowCutoff(0.5);
    sim.configSlipspringHighCutoff(2.5);
    REQUIRE(sim.createSlipSprings(2) == 2);
    REQUIRE_NOTHROW(sim.validateState());
    CHECK(sim.relocateSlipSprings(1.0) == 0);
    REQUIRE_NOTHROW(sim.validateState());
  }

  SECTION("Zero-velocity step 0 uses target kbT and kbT <= 0 is handled safely")
  {
    pe::Universe u(10.0, 10.0, 10.0);
    u.addAtoms({ 1, 2, 3, 4 },
               { 1, 1, 1, 1 },
               { 1.0, 2.0, 1.0, 2.0 },
               { 1.0, 1.0, 2.0, 3.0 },
               { 1.0, 1.0, 1.0, 1.0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 });
    u.addBonds({ 1, 3 }, { 2, 4 });
    // Add slip-spring (type 9) between atom 1 and atom 3 (distance 1.0)
    u.addBonds(1, { 1 }, { 3 }, { 9 });

    pcd::DPDSimulator sim(u, 2, 9, false, "zero_vel_seed");
    CHECK(sim.getTemperature() == Catch::Approx(0.0));
    CHECK(sim.getTargetKbT() == Catch::Approx(1.0));
    sim.configNumStepsDPD(10);
    sim.configNumStepsMC(5);
    sim.configShiftOneAtATime(false);
    sim.configShiftPossibilityEmpty(false);

    // Negative kbT must throw
    CHECK_THROWS_AS(sim.shiftSlipSprings(-1.0), std::invalid_argument);
    CHECK_THROWS_AS(sim.relocateSlipSprings(-1.0), std::invalid_argument);

    // kbT == 0.0 must not divide by zero or accept uphill moves
    for (int i = 0; i < 50; ++i) {
      CHECK(sim.shiftSlipSprings(0.0) == 0);
    }
    CHECK_NOTHROW(sim.relocateSlipSprings(0.0));

    // Running simulation from zero initial velocities with MC enabled at step 0
    REQUIRE_NOTHROW(sim.runSimulation(2, true));
    REQUIRE_NOTHROW(sim.validateState());
  }

  SECTION("Slip-spring shift forwards kbT and uses 0.5 * k * r^2 energy")
  {
    // Two 2-bead chains:
    //   Chain A: atom 1 at (1, 1, 1) -- atom 2 at (2, 1, 1)
    //   Chain B: atom 3 at (1, 2, 1) -- atom 4 at (2, 3, 1)
    // Slip-spring initially between (1, 3): r0^2 = 1.0
    // Shifting both ends simultaneously always proposes (2, 4): r1^2 = 4.0
    // With k = 2.0: dU = 0.5 * k * (4.0 - 1.0) = 3.0.
    // At kbT = 3.0, acceptance probability = exp(-3.0 / 3.0) = exp(-1) ~=
    // 0.3679 (Without the 0.5 factor, it would be exp(-2) ~= 0.1353.)
    pe::Universe u(10.0, 10.0, 10.0);
    u.addAtoms({ 1, 2, 3, 4 },
               { 1, 1, 1, 1 },
               { 1.0, 2.0, 1.0, 2.0 },
               { 1.0, 1.0, 2.0, 3.0 },
               { 1.0, 1.0, 1.0, 1.0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 });
    u.addBonds({ 1, 3 }, { 2, 4 });
    u.addBonds(1, { 1 }, { 3 }, { 9 });

    pcd::DPDSimulator sim(u, 2, 9, false, "detailed_balance_seed");
    sim.configSpringConstant(2.0);
    sim.configShiftOneAtATime(false);
    sim.configShiftPossibilityEmpty(false);

    const int numTrials = 2000;
    int numAccepted = 0;
    for (int trial = 0; trial < numTrials; ++trial) {
      const int accepted = sim.shiftSlipSprings(3.0);
      numAccepted += accepted;
      if (accepted == 1) {
        // Return move (2, 4) -> (1, 3) is strictly downhill (dU = -3.0)
        // and must always be accepted.
        REQUIRE(sim.shiftSlipSprings(3.0) == 1);
      }
    }
    const double acceptanceRate =
      static_cast<double>(numAccepted) / static_cast<double>(numTrials);
    CHECK_THAT(acceptanceRate,
               Catch::Matchers::WithinAbs(std::exp(-1.0), 0.05));

    // With kbT = 0.0, uphill move must never be accepted (verifies kbT is
    // forwarded rather than defaulting to 1.0)
    int acceptedAtZeroT = 0;
    for (int trial = 0; trial < 200; ++trial) {
      acceptedAtZeroT += sim.shiftSlipSprings(0.0);
    }
    CHECK(acceptedAtZeroT == 0);
  }

  SECTION("shiftOneAtATime handles spring collapse without stale partner state")
  {
    // Single 4-bead chain: 1 - 2 - 3 - 4.
    // Slip-spring between atom 1 (end) and atom 2: shifting partner A (atom 1)
    // deterministically moves it onto atom 2 (collapse: A == B), which must
    // trigger full relocation without operating on stale partner state.
    pe::Universe u(10.0, 10.0, 10.0);
    u.addAtoms({ 1, 2, 3, 4 },
               { 1, 1, 1, 1 },
               { 1.0, 2.0, 3.0, 4.0 },
               { 1.0, 1.0, 1.0, 1.0 },
               { 1.0, 1.0, 1.0, 1.0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 });
    u.addBonds({ 1, 2, 3 }, { 2, 3, 4 });
    u.addBonds(1, { 1 }, { 2 }, { 9 });

    pcd::DPDSimulator sim(u, 2, 9, false, "collapse_seed");
    sim.configSlipspringLowCutoff(0.5);
    sim.configSlipspringHighCutoff(2.5);
    sim.configShiftOneAtATime(true);
    sim.configShiftPossibilityEmpty(false);

    REQUIRE_NOTHROW(sim.shiftSlipSprings(1.0));
    REQUIRE_NOTHROW(sim.validateState());
  }

  SECTION(
    "resetBondDuplicationPenalty covers last atom and bondType=0 sentinel")
  {
    pe::Universe u(10.0, 10.0, 10.0);
    u.addAtoms({ 1, 2, 3 },
               { 1, 1, 1 },
               { 1.0, 2.0, 3.0 },
               { 1.0, 1.0, 1.0 },
               { 1.0, 1.0, 1.0 },
               { 0, 0, 0 },
               { 0, 0, 0 },
               { 0, 0, 0 });
    // Normal bonds (1-2, 2-3) and duplicate slip-spring on last atom (2-3)
    u.addBonds(2, { 1, 2 }, { 2, 3 }, { 1, 1 });
    u.addBonds(1, { 2 }, { 3 }, { 7 });

    pcd::DPDSimulator sim(u, 2, 7, false, "penalty_and_bondtype_seed");
    REQUIRE_NOTHROW(sim.resetBondDuplicationPenalty());
    REQUIRE_NOTHROW(sim.validateState());

    // bondType = 0 sentinel should inherit slipspringBondType (7)
    sim.configSlipspringLowCutoff(0.5);
    sim.configSlipspringHighCutoff(2.5);
    REQUIRE(sim.createSlipSprings(1, 0) == 1);
    pe::Universe uAfterDefault = sim.getUniverse(true);
    std::map<std::string, std::vector<long int>> edgesDefault =
      uAfterDefault.getEdges();
    REQUIRE(edgesDefault["edge_type"].size() == 4);
    CHECK(edgesDefault["edge_type"].back() == 7);

    // Explicit non-zero bondType should be preserved
    REQUIRE(sim.createSlipSprings(1, 5) == 1);
    pe::Universe uAfterExplicit = sim.getUniverse(true);
    std::map<std::string, std::vector<long int>> edgesExplicit =
      uAfterExplicit.getEdges();
    REQUIRE(edgesExplicit["edge_type"].size() == 5);
    CHECK(edgesExplicit["edge_type"].back() == 5);
  }

  SECTION("PBC-crossing slip-spring shifts, partner swaps, and relocations "
          "preserve exact unwrapped bond vectors")
  {
    // Two chains crossing x = 0 / 10 periodic boundary:
    //   Chain A: 1 (9.2) -- 2 (9.8) -- 3 (0.4) -- 4 (1.0) at y = 1.0
    //   Chain B: 5 (9.2) -- 6 (9.8) -- 7 (0.4) -- 8 (1.0) at y = 2.0
    // Slip-spring initially between atom 2 (9.8, 1, 1) and atom 7 (0.4, 2, 1):
    // true minimum-image vector A->B is (+0.6, +1.0, 0.0), length^2 = 1.36.
    pe::Universe u(10.0, 10.0, 10.0);
    u.addAtoms({ 1, 2, 3, 4, 5, 6, 7, 8 },
               { 1, 1, 1, 1, 1, 1, 1, 1 },
               { 9.2, 9.8, 0.4, 1.0, 9.2, 9.8, 0.4, 1.0 },
               { 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0 },
               { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },
               { 0, 0, 0, 0, 0, 0, 0, 0 },
               { 0, 0, 0, 0, 0, 0, 0, 0 },
               { 0, 0, 0, 0, 0, 0, 0, 0 });
    u.addBonds({ 1, 2, 3, 5, 6, 7 }, { 2, 3, 4, 6, 7, 8 });
    u.addBonds(1, { 2 }, { 7 }, { 9 });

    // Compare tracked-offset PBC (default false) vs minimum-image PBC (true)
    // across both 2-end and 1-end shifts and relocations with identical RNG
    // seeds
    for (const bool oneAtATime : { false, true }) {
      pcd::DPDSimulator simTracked(u, 2, 9, false, "pbc_mc_seed", false);
      pcd::DPDSimulator simMic(u, 2, 9, false, "pbc_mc_seed", true);
      CHECK_FALSE(simTracked.assumesBoxLargeEnough());
      CHECK(simMic.assumesBoxLargeEnough());

      simTracked.configSlipspringLowCutoff(0.3);
      simTracked.configSlipspringHighCutoff(2.5);
      simTracked.configShiftOneAtATime(oneAtATime);
      simMic.configSlipspringLowCutoff(0.3);
      simMic.configSlipspringHighCutoff(2.5);
      simMic.configShiftOneAtATime(oneAtATime);

      REQUIRE(simTracked.getBondLengths().isApprox(simMic.getBondLengths()));

      for (int step = 0; step < 40; ++step) {
        const int sTracked = simTracked.shiftSlipSprings(1.5);
        const int sMic = simMic.shiftSlipSprings(1.5);
        CHECK(sTracked == sMic);
        REQUIRE_NOTHROW(simTracked.validateState());
        REQUIRE(
          simTracked.getBondLengths().isApprox(simMic.getBondLengths(), 1e-10));
        CHECK(simTracked.getBondLengths().maxCoeff() < 5.0);

        const int rTracked = simTracked.relocateSlipSprings(1.5);
        const int rMic = simMic.relocateSlipSprings(1.5);
        CHECK(rTracked == rMic);
        REQUIRE_NOTHROW(simTracked.validateState());
        REQUIRE(
          simTracked.getBondLengths().isApprox(simMic.getBondLengths(), 1e-10));
      }
    }

    // Intra-chain PBC-crossing slip-spring where 2-end shift can swap A <-> B
    // or step newPartnerA onto oldPartnerB: chain 1 (9.2) - 2 (9.8) - 3 (0.4) -
    // 4 (1.0)
    pe::Universe uSwap(10.0, 10.0, 10.0);
    uSwap.addAtoms({ 1, 2, 3, 4 },
                   { 1, 1, 1, 1 },
                   { 9.2, 9.8, 0.4, 1.0 },
                   { 1.0, 1.0, 1.0, 1.0 },
                   { 1.0, 1.0, 1.0, 1.0 },
                   { 0, 0, 0, 0 },
                   { 0, 0, 0, 0 },
                   { 0, 0, 0, 0 });
    uSwap.addBonds({ 1, 2, 3 }, { 2, 3, 4 });
    uSwap.addBonds(1, { 2 }, { 3 }, { 9 });

    pcd::DPDSimulator simSwap(uSwap, 2, 9, false, "swap_pbc_seed", false);
    simSwap.configSlipspringLowCutoff(0.3);
    simSwap.configSlipspringHighCutoff(2.0);
    simSwap.configShiftOneAtATime(false);
    for (int step = 0; step < 30; ++step) {
      REQUIRE_NOTHROW(simSwap.shiftSlipSprings(2.0));
      REQUIRE_NOTHROW(simSwap.validateState());
      CHECK(simSwap.computeBondLength(3) < 2.5);
    }
  }

  SECTION("Collapse with no eligible relocation pairs reverts safely and "
          "input validation catches invalid parameters")
  {
    // Two beads far apart (distance 3.0) with highCutoff = 2.0:
    // collapse cannot find any candidate pair within (lowCutoff, highCutoff],
    // so it must revert cleanly rather than leaving a self-loop.
    pe::Universe uSparse(10.0, 10.0, 10.0);
    uSparse.addAtoms({ 1, 2 },
                     { 1, 1 },
                     { 1.0, 4.0 },
                     { 1.0, 1.0 },
                     { 1.0, 1.0 },
                     { 0, 0 },
                     { 0, 0 },
                     { 0, 0 });
    uSparse.addBonds({ 1 }, { 2 });
    uSparse.addBonds(1, { 1 }, { 2 }, { 9 });

    pcd::DPDSimulator simSparse(uSparse, 2, 9, false, "sparse_collapse_seed");
    simSparse.configSlipspringLowCutoff(0.5);
    simSparse.configSlipspringHighCutoff(2.0);
    simSparse.configShiftOneAtATime(true);
    simSparse.configShiftPossibilityEmpty(false);

    CHECK(simSparse.shiftSlipSprings(1.0) == 0);
    REQUIRE_NOTHROW(simSparse.validateState());
    CHECK(simSparse.computeBondLength(1) == Catch::Approx(3.0));

    // Validation & edge cases
    CHECK_THROWS_AS(simSparse.createSlipSprings(-1), std::invalid_argument);
    CHECK_THROWS_AS(simSparse.createSlipSprings(1, -1), std::invalid_argument);
    CHECK(simSparse.createSlipSprings(0) == 0);
    CHECK_THROWS_AS(simSparse.configSlipspringLowCutoff(-0.1),
                    std::invalid_argument);
    CHECK_THROWS_AS(simSparse.configSigma(-1.0), std::invalid_argument);
    CHECK_THROWS_AS(simSparse.configGamma(-1.0), std::invalid_argument);

    simSparse.configSigma(3.0);
    simSparse.configGamma(2.25);
    CHECK(simSparse.getTargetKbT() == Catch::Approx(2.0));

    simSparse.configAssumeBoxLargeEnough(true);
    CHECK(simSparse.assumesBoxLargeEnough());
    simSparse.configAssumeBoxLargeEnough(false);
    CHECK_FALSE(simSparse.assumesBoxLargeEnough());

    CHECK(pcd::DPDSimulator::computeTemperature(Eigen::VectorXd(0)) ==
          Catch::Approx(0.0));
  }
}
