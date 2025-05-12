#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/DPDSimulator.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
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
TEST_CASE("DPD Simulator Restart Files Work", "[analysis][DPDSimulator]")
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

TEST_CASE("DPD Simulator Can Cross-link", "[analysis][DPDSimulator][long]")
{
  std::cout << "Running test \"DPD Simulator Can Cross-link\"" << std::endl;
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

TEST_CASE("DPD Simulator Computes Correct Forces", "[analysis][DPDSimulator]")
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
  REQUIRE_NOTHROW(simulator.runSimulation(10, false));
  REQUIRE_NOTHROW(simulator.runSimulation(10, true));

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
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "15th_seed");

    simulator.createSlipSprings(100, 2);

#ifdef OPENMP_FOUND
    // we cannot have more than 1 thread, otherwise the random number generator
    // will not play nicely.
    omp_set_num_threads(1);
#endif

    // invoke copy-constructor
    pcd::DPDSimulator simulator2 = simulator;

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
  std::cout
    << "Running test \"For large systems the PBC method does not matter\""
    << std::endl;
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  std::string inputFile = suspectedPath +
                          "/structure/"
                          "3d-diamond-lattice_3x3x3_a_23_d_3_v_0.structure.out";
  if (std::filesystem::exists(inputFile)) {
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcd::DPDSimulator simulator =
      pcd::DPDSimulator(universe, 2, 9, false, "19th_seed");

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

#ifdef OPENMP_FOUND
    // we cannot have more than 1 thread, otherwise the random number generator
    // will not play nicely.
    omp_set_num_threads(1);
#endif

    // invoke copy-constructor
    pcd::DPDSimulator simulator2 = pcd::DPDSimulator(simulator);

    // switch to "common" PBC
    simulator2.configAssumeBoxLargeEnough();

    REQUIRE(!simulator.assumesBoxLargeEnough());
    REQUIRE(simulator2.assumesBoxLargeEnough());

    simulator.reseedRandomness("20th_seed");
    simulator.refreshCurrentState();
    simulator2.reseedRandomness("20th_seed");
    simulator2.refreshCurrentState();

    CHECK(simulator.getUniformRandBetween0And1() ==
          simulator2.getUniformRandBetween0And1());
    CHECK_THAT(simulator.getTemperature(),
               Catch::Matchers::WithinRel(simulator2.getTemperature()));
    CHECK(simulator.getBondLengths().isApprox(simulator2.getBondLengths()));
    CHECK(simulator.getStressTensor().isApprox(simulator2.getStressTensor()));

    // not sure what value is sensible here...
    simulator.runSimulation(50, true);
    simulator2.runSimulation(50, true);

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
