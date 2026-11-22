#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPLinearSystemSolver.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPNetwork.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPNonlinearOptimizationSolver.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPSimulator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPSlipLinkRelaxationSolver.h"
#include <array>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::sim::mehp;

void
testGradient(const pcm::MEHPForceEvaluator& forceEvaluator)
{
  pe::Universe universe = pe::Universe(4.0, 4.0, 4.0);
  // how this looks like:
  // 1-2
  // | |
  // 4-3
  universe.addAtoms({ { 1, 2, 3, 4 } },
                    { { 2, 2, 2, 2 } },
                    { { 0.0, 1.0, 2.0, 3.0 } }, // x
                    { { 0.0, 0.0, 0.0, 0.0 } }, // y
                    { { 0.0, 0.0, 0.0, 0.0 } }, // z
                    { { 0, 0, 0, 0 } },
                    { { 0, 0, 0, 0 } },
                    { { 0, 0, 0, 0 } });
  universe.addBonds({ { 1, 2, 3, 4 } }, { { 2, 3, 4, 1 } });
  // assemble network
  pcm::MEHPNetwork net;
  Eigen::VectorXd coordinates = Eigen::VectorXd(12);
  coordinates << 0.0, 0.0, 0.0, // 1
    1.0, 0.0, 0.0,              // 2
    2.0, 0.0, 0.0,              // 3
    3.0, 0.0, 0.0;              // 4
  net.coordinates = coordinates;
  net.springCoordinateIndexA = Eigen::ArrayXi::Zero(12);
  net.springCoordinateIndexB = Eigen::ArrayXi::Zero(12);
  net.springIndexA = Eigen::ArrayXi(4);
  net.springIndexA << 0, 1, 2, 3;
  net.springIndexB = Eigen::ArrayXi(4);
  net.springIndexB << 1, 2, 3, 0;
  for (int i = 0; i < 4; ++i) {
    for (int dir = 0; dir < 3; ++dir) {
      net.springCoordinateIndexA[i * 3 + dir] = net.springIndexA[i] * 3 + dir;
      net.springCoordinateIndexB[i * 3 + dir] = net.springIndexB[i] * 3 + dir;
    }
  }
  net.nrOfLinks = 4;
  net.nrOfNodes = 4;
  net.nrOfStrands = 4;
  net.nrOfSprings = 4;
  net.vol = universe.getVolume();
  net.box = universe.getBox();
  net.L[0] = 4.0;
  net.L[1] = 4.0;
  net.L[2] = 4.0;
  net.boxHalfs[0] = 2.0;
  net.boxHalfs[1] = 2.0;
  net.boxHalfs[2] = 2.0;
  net.springContourLength = Eigen::VectorXd::Constant(4, 1.0);
  net.meanSpringContourLength = 1.0;
  net.meanStrandContourLength = 1.0;
  net.springBoxOffset =
    universe.getBox().getOffset(net.coordinates(net.springCoordinateIndexB) -
                                net.coordinates(net.springCoordinateIndexA));

  REQUIRE(net.coordinates.size() == 12);
  // test that Eigen does as expected
  Eigen::VectorXd coordinatesSpringEndA =
    net.coordinates(net.springCoordinateIndexA);
  Eigen::VectorXd coordinatesSpringEndB =
    net.coordinates(net.springCoordinateIndexB);
  std::vector<double> expectedCoordsA = {
    0.0, 0.0, 0.0, // 1
    1.0, 0.0, 0.0, // 2
    2.0, 0.0, 0.0, // 3
    3.0, 0.0, 0.0, // 4
  };
  std::vector<double> expectedCoordsB = {
    1.0, 0.0, 0.0, // 2
    2.0, 0.0, 0.0, // 3
    3.0, 0.0, 0.0, // 4
    0.0, 0.0, 0.0, // 1
  };
  Eigen::VectorXd springDistances =
    (coordinatesSpringEndB - coordinatesSpringEndA);
  for (size_t i = 0; i < 12; ++i) {
    REQUIRE(coordinatesSpringEndA[static_cast<Eigen::Index>(i)] ==
            expectedCoordsA[i]);
    REQUIRE(coordinatesSpringEndB[static_cast<Eigen::Index>(i)] ==
            expectedCoordsB[i]);
    REQUIRE(springDistances[static_cast<Eigen::Index>(i)] ==
            expectedCoordsB[i] - expectedCoordsA[i]);
  }
  // setup gradient and coordinates
  const double h = 1.e-5;
  std::array<double, 12> grad{};
  Eigen::VectorXd u = Eigen::VectorXd::Zero(12);
  const Eigen::VectorXd springVectors0 = net.evaluateSpringVectors(u, false);
  // actual computation to test gradient
  for (size_t i = 0; i < 12; ++i) {
    // evaluate gradient
    (void)forceEvaluator.evaluateNetworkEnergyAndGradient(
      net, springVectors0, false, grad.data());
    if (i % 3 != 0) {
      // in y and z direction, we expect no spring distance -> 0 gradient
      CHECK(grad[i] == Catch::Approx(0.0));
    } else {
      // test finite difference vs. gradient
      u[static_cast<Eigen::Index>(i)] = -h;
      const double fm = forceEvaluator.evaluateNetworkEnergyAndGradient(
        net, net.evaluateSpringVectors(u, false), false, nullptr);
      u[static_cast<Eigen::Index>(i)] = h;
      const double fp = forceEvaluator.evaluateNetworkEnergyAndGradient(
        net, net.evaluateSpringVectors(u, false), false, nullptr);
      u[static_cast<Eigen::Index>(i)] = 0.0; // reset
      // require gradient to be similar to finite difference
      if (std::abs(grad[i]) == 0.0) {
        CHECK(std::abs(grad[i]) ==
              Catch::Approx(std::abs(fp - fm)).margin(0.000002));
      } else {
        CHECK(grad[i] == Catch::Approx((1.0 / (2.0 * h)) * (fp - fm)));
      }
    }
  }
}

TEST_CASE("MEHP Force Relaxation2 computes correct gradients",
          "[analysis][MEHPSimulator][HookeanForceEvaluator]")
{
  std::cout
    << "Running test \"MEHP Force Relaxation2 computes correct gradients\""
    << std::endl;
  SECTION("Test HookeanForceEvaluator force gradient")
  {
    const pcm::HookeanForceEvaluator forceEvaluatorInstance =
      pcm::HookeanForceEvaluator(1.0);
    testGradient(forceEvaluatorInstance);
  }

  SECTION("Test NonGaussianSpringForceEvaluator force gradient")
  {
    const pcm::NonGaussianSpringForceEvaluator forceEvaluatorInstance =
      pcm::NonGaussianSpringForceEvaluator(1.0, 1.0);
    testGradient(forceEvaluatorInstance);
  }
}

TEST_CASE("MEHP Force Relaxation does not collapse",
          "[analysis][MEHPSimulator][HookeanForceEvaluator]")
{
  std::cout << "Running test \"MEHP Force Relaxation does not collapse\""
            << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  /**
   * @brief A grid of two rows, each one bead between the two crosslinkers
   *
   */
  universe.addAtoms(
    { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } },
    { { 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1 } },
    { { 0.,
        2.5,
        5,
        7.5,
        0.1,
        2.5,
        5,
        7.5,
        -0.1,
        5.,
        0.,
        5. } }, // x with slight (0.1) deviation, so we don't start perfect
    { { 0.1, 0., -0.1, 0., 5., 5., 5., 5., 2.5, 2.5, 7.5, 7.5 } }, // y
    { { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. } },        // z
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } });
  universe.addBonds(
    { { 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7 } },
    { { 2, 9, 4, 11, 2, 4, 10, 12, 9, 11, 6, 8, 6, 8, 10, 12 } });

  SECTION("Running conventional MEHP")
  {
    pcm::MEHPSimulator forceRelaxerConventional =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         true,
                         false,
                         true);
    const pcm::MEHPNonlinearSolverConfig nlConfig;
    REQUIRE_NOTHROW(forceRelaxerConventional.run(nlConfig));
    REQUIRE(forceRelaxerConventional.getNrOfIterations() > 0);
    CHECK(forceRelaxerConventional.getExitReason() ==
          pcm::ExitReason::F_TOLERANCE);
    CHECK(forceRelaxerConventional.getNrOfActiveSprings() ==
          forceRelaxerConventional.getNrOfSprings());
    // compare to what we expect
    CHECK(forceRelaxerConventional.getNrOfActiveSprings() == 8);
    CHECK(forceRelaxerConventional.getNrOfActiveNodes() == 4);
  }
  SECTION("Running new MEHP")
  {
    pcm::MEHPSimulator forceRelaxerNew =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         true,
                         false,
                         false);
    const pcm::MEHPNonlinearSolverConfig nlConfig;
    REQUIRE_NOTHROW(forceRelaxerNew.run(nlConfig));
    REQUIRE(forceRelaxerNew.getNrOfIterations() > 0);
    CHECK(forceRelaxerNew.getExitReason() == pcm::ExitReason::F_TOLERANCE);
    CHECK(forceRelaxerNew.getNrOfActiveSprings() ==
          forceRelaxerNew.getNrOfSprings());

    // compare to what we expect
    CHECK(forceRelaxerNew.getNrOfActiveSprings() == 8);
    CHECK(forceRelaxerNew.getNrOfActiveNodes() == 4);
    CHECK(forceRelaxerNew.getAverageSpringLength() == Catch::Approx(5.0));
    CHECK_THAT(forceRelaxerNew.getGammaFactor(5.),
               Catch::Matchers::WithinAbs(2.5, 1e-2));
    CHECK_THAT(forceRelaxerNew.getGammaFactors(5.).mean(),
               Catch::Matchers::WithinAbs(2.5, 1e-2));
  }
};

TEST_CASE(
  "MEHP Force Relaxation2 runs",
  "[analysis][MEHPSimulator][HookeanForceEvaluator][long]")
{
  std::cout << "Running test \"MEHP Force Relaxation2 runs\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  const std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Relaxation2 3D case")
  {
    const std::string largeInputFile =
      suspectedPath +
      "/structure/xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      const pe::Universe universe2 = universeSeq.atIndex(0);

      const double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      pcm::MEHPNetwork net = pcm::MEHPNetwork::fromCrosslinkerUniverse(
        universe2.getBox(),
        universe2.toCrosslinkerUniverse(2),
        {},
        2,
        pcm::EntanglementMode::NONE,
        true,
        false,
        false,
        false);
      net.universeNrOfAtoms = universe2.getNrOfAtoms();
      pcm::MEHPSimulator forceRelaxer2 = pcm::MEHPSimulator(
        net, std::make_shared<pcm::HookeanForceEvaluator>(1.0), false);
      forceRelaxer2.configMeanBondLength(universe2.computeMeanBondLength());
      forceRelaxer2.configMeanSquaredBondLength(
        universe2.computeMeanSquaredBondLength());
      REQUIRE(forceRelaxer2.getExitReason() == pcm::ExitReason::UNSET);
      CHECK(forceRelaxer2.getNrOfSprings() == 8142);
      REQUIRE(forceRelaxer2.getNrOfIterations() == 0);
      REQUIRE(forceRelaxer2.getNetwork().vol ==
              Catch::Approx(universe2.getVolume()));
      CHECK(forceRelaxer2.getNetwork().vol ==
            Catch::Approx(97.383096 * 97.383096 * 97.383096));
      // initial system values
      CHECK(forceRelaxer2.getPressure() == Catch::Approx(0.0050521117));
      CHECK(forceRelaxer2.getEnergy() * 79 == Catch::Approx(552894.1903005134));
      CHECK(forceRelaxer2.getAverageContourLength() == Catch::Approx(79.0));
      const Eigen::VectorXd contourLengths =
        forceRelaxer2.getSpringContourLength();
      for (int i = 0; i < contourLengths.size(); ++i) {
        CHECK(contourLengths[i] ==
              Catch::Approx(forceRelaxer2.getAverageContourLength()));
      }
      CHECK(forceRelaxer2.getResidualNorm() * 79 ==
            Catch::Approx(1457.465048151));

      // run force relaxation
      const pcm::MEHPNonlinearSolverConfig nlConfig;
      REQUIRE_NOTHROW(forceRelaxer2.run(nlConfig));
      CHECK(forceRelaxer2.getNrOfSprings() == 8142);
      CHECK(forceRelaxer2.getNrOfIterations() > 1);
      CHECK(forceRelaxer2.getSolubleWeightFraction() > 0.);
      CHECK(forceRelaxer2.getSolubleWeightFraction() < 1.);

      // conversion factors
      const double kb = 1.381e-23; // Boltzmann, J/K
      const double T = 300.;       // Temperature, K
      const double sigmaToNm = 0.616;
      const double sigmaToM = sigmaToNm * 1.e-9;
      const double slope = 0.00393 / (sigmaToNm * sigmaToNm); // sigma^2/(g/mol)
      const double beadMass = 161.;                           // g/mol
      const double Nb = 80.; // nr of bonds per strand
      const double conversionFactor =
        3. * kb * T / (slope * beadMass * (Nb / 79.)); // J/sigma^2
      CHECK((conversionFactor / 79.) / (sigmaToM * sigmaToM) ==
            Catch::Approx(0.000245543));
      const double nu =
        nrOfChains / (forceRelaxer2.getNetwork().vol * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      CHECK(forceRelaxer2.getPressure() ==
            Catch::Approx(0.0019468984)); // LJ Units [?]
      CHECK(forceRelaxer2.getPressure() * conversionFactor /
              (sigmaToM * sigmaToM * sigmaToM) ==
            Catch::Approx(
              61308.0398526509)); // shear modulus from the pressure, MPa
      const double expectedNb2 = slope * Nb * beadMass;
      const double bondR02 = forceRelaxer2.getDefaultR0Square();
      const double strandR02 =
        forceRelaxer2.getAverageContourLength() * bondR02;
      const double nb2Correction = strandR02 / expectedNb2;
      const double gammaCorrectionFactor = nb2Correction;
      const double gammaFactor =
        forceRelaxer2.getGammaFactor(bondR02, static_cast<long int>(nrOfChains));
      CHECK(gammaFactor * strandR02 ==
            Catch::Approx(
              42.6130241132)); // as from conversion-less Mathematica script
      CHECK_THAT(gammaFactor * gammaCorrectionFactor * kb * T * nu,
                 Catch::Matchers::WithinRel(61338.8913392133,
                                            1e-3)); // ANT shear modulus, Pa
      CHECK(gammaFactor * gammaCorrectionFactor ==
            Catch::Approx(0.3194444644)); // "correct" gamma factor
      CHECK(forceRelaxer2.getExitReason() == pcm::ExitReason::F_TOLERANCE);
      // TODO: find better, more accurate tests here
      CHECK(forceRelaxer2.getNrOfActiveNodes() > 1);
      CHECK(forceRelaxer2.getNrOfActiveSprings() > 1);
      CHECK(forceRelaxer2.getAverageSpringLength() > 1.0);
      CHECK(forceRelaxer2.getEffectiveFunctionalityOfAtoms().size() ==
            forceRelaxer2.getNrOfNodes());

    } else {
      std::cout << "Skipping large file PDMS MEHP run" << std::endl;
      REQUIRE(true);
    }
  }

  SECTION("2D case")
  {
    REQUIRE(std::filesystem::exists(suspectedPath));
    universeSeq.initializeFromDataSequence(
      { { suspectedPath +
          "/structure/equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
          "333_2d_t_7500001.structure.out" } });
    REQUIRE(universeSeq.getLength() == 1);
    const pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPSimulator forceRelaxer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         true,
                         false,
                         false);
    REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::UNSET);
    REQUIRE(forceRelaxer.getNrOfIterations() == 0);
    REQUIRE(forceRelaxer.getNetwork().vol ==
            Catch::Approx(universe.getVolume()));
    pcm::MEHPNonlinearSolverConfig shortConfig;
    shortConfig.algorithm = pcm::NLoptAlgorithm::LD_LBFGS;
    shortConfig.maxNrOfSteps = 5;
    forceRelaxer.run(shortConfig);
    REQUIRE(forceRelaxer.getNrOfNodes() != universe.getNrOfAtoms());
    REQUIRE(forceRelaxer.getNrOfIterations() <= 5);
    REQUIRE(forceRelaxer.getNrOfIterations() >= 1);
    REQUIRE(universe.getAtomsOfType(2).size() == 7200);
    REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::MAX_STEPS);
    // run again, this time fully
    pcm::MEHPSimulator forceRelaxer2 =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         true,
                         false,
                         false);
    pcm::MEHPNonlinearSolverConfig fullConfig;
    fullConfig.algorithm = pcm::NLoptAlgorithm::LD_LBFGS;
    fullConfig.maxNrOfSteps = 10000;
    fullConfig.xtol = 1e-5;
    fullConfig.ftol = 1e-18;
    forceRelaxer2.run(fullConfig);
    REQUIRE(forceRelaxer2.getNrOfIterations() > 5);
    CHECK(forceRelaxer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    const long int nSprings =
      static_cast<long int>(forceRelaxer2.getNrOfSprings());
    const double gammaFactor2 = forceRelaxer2.getGammaFactor(25., nSprings);
    CHECK(gammaFactor2 == Catch::Approx(2. / 150.).epsilon(0.001));
    const Eigen::Matrix3d stressTensor = forceRelaxer2.getStressTensor();
    CHECK(forceRelaxer2.getPressure() ==
          Catch::Approx(
            (stressTensor(0, 0) + stressTensor(1, 1) + stressTensor(2, 2)) / 3.)
            .epsilon(0.02));
    CHECK(forceRelaxer2.getResidualNorm() > 0.0);
    CHECK(forceRelaxer2.getEnergy() > 0.0);
    // TODO: find better, more accurate tests here
    CHECK(forceRelaxer2.getNrOfActiveNodes() > 1);
    CHECK(forceRelaxer2.getNrOfActiveSprings() > 1);
    CHECK(forceRelaxer2.getAverageSpringLength() > 1.0);
    CHECK(forceRelaxer2.getEffectiveFunctionalityOfAtoms().size() ==
          forceRelaxer2.getNrOfNodes());
    CHECK(forceRelaxer2.getSolubleWeightFraction() < 1.);

    const pe::Universe universe3 = forceRelaxer2.getCrosslinkerVerse();
    CHECK(universe3.getNrOfAtoms() == forceRelaxer2.getNrOfNodes());
    CHECK(universe3.getNrOfBonds() == forceRelaxer2.getNrOfSprings());
    CHECK(universe3.getAtomsOfType(2).size() == universe3.getNrOfAtoms());
  }

  SECTION("Simple melt case")
  {
    const std::string dataFile = suspectedPath + "/lammps_data_file.out";
    universeSeq.initializeFromDataSequence({ { dataFile } });
    const pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPSimulator forceRelaxer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         true,
                         true,
                         false,
                         false);
    // not a proper network -> 0 springs
    CHECK(forceRelaxer.getNrOfSprings() == 0);
    CHECK(forceRelaxer.getSolubleWeightFraction() == 1.);
  }
}

TEST_CASE(
  "MEHP Force Relaxation2 runs with non-gaussian force evaluator",
  "[analysis][MEHPSimulator][NonGaussianSpringForceEvaluator][long]")
{
  std::cout << "Running test \"MEHP Force Relaxation2 runs with non-gaussian "
               "force evaluator\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  const std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("3D case")
  {
    const std::string largeInputFile =
      suspectedPath +
      "/structure/xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      const pe::Universe universe2 = universeSeq.atIndex(0);

      const double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      const std::shared_ptr<pcm::NonGaussianSpringForceEvaluator>
        nonGaussianForceEvaluator =
          std::make_shared<pcm::NonGaussianSpringForceEvaluator>(1.0, 0.98);
      pcm::MEHPNetwork net = pcm::MEHPNetwork::fromCrosslinkerUniverse(
        universe2.getBox(),
        universe2.toCrosslinkerUniverse(2),
        {},
        2,
        pcm::EntanglementMode::NONE,
        true,
        false,
        false,
        false);
      net.universeNrOfAtoms = universe2.getNrOfAtoms();
      pcm::MEHPSimulator forceRelaxer2 =
        pcm::MEHPSimulator(net, nonGaussianForceEvaluator, false);
      forceRelaxer2.configMeanBondLength(universe2.computeMeanBondLength());
      forceRelaxer2.configMeanSquaredBondLength(
        universe2.computeMeanSquaredBondLength());
      REQUIRE(forceRelaxer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE(forceRelaxer2.getNrOfIterations() == 0);
      REQUIRE(forceRelaxer2.getNetwork().vol ==
              Catch::Approx(universe2.getVolume()));
      CHECK(forceRelaxer2.getNetwork().vol ==
            Catch::Approx(97.383096 * 97.383096 * 97.383096));
      // initial system values
      CHECK(forceRelaxer2.getEnergy() ==
            Catch::Approx(22103.6441026747).epsilon(1e-5));
      CHECK(forceRelaxer2.getResidualNorm() ==
            Catch::Approx(58.7568113327).epsilon(1e-5));
      CHECK(forceRelaxer2.getPressure() *
              forceRelaxer2.getNetwork().meanSpringContourLength ==
            Catch::Approx(0.3991168239));
      CHECK(forceRelaxer2.getAverageContourLength() == Catch::Approx(79.0));
      const Eigen::VectorXd contourLengths =
        forceRelaxer2.getSpringContourLength();
      for (int i = 0; i < contourLengths.size(); ++i) {
        CHECK(contourLengths[i] ==
              Catch::Approx(forceRelaxer2.getAverageContourLength()));
      }
      const pcm::MEHPNonlinearSolverConfig nlConfig;
      REQUIRE_NOTHROW(forceRelaxer2.run(nlConfig));
      CHECK(forceRelaxer2.getNrOfSprings() == 8142);
      CHECK(forceRelaxer2.getNrOfIterations() > 1);
      CHECK(forceRelaxer2.getSolubleWeightFraction() > 0.);
      CHECK(forceRelaxer2.getSolubleWeightFraction() < 1.);

      // conversion factors
      double kb = 1.381e-23; // Boltzmann, J/K
      double T = 300.;       // Temperature, K
      double sigmaToNm = 0.616;
      double sigmaToM = sigmaToNm * 1.e-9;
      double slope = 0.00393 / (sigmaToNm * sigmaToNm); // sigma^2/(g/mol)
      double beadMass = 161.;                           // g/mol
      double Nb = 80.; // nr of beads per strand
      double conversionFactor =
        (forceRelaxer2.getNetwork().meanSpringContourLength) * 3. * kb * T /
        (slope * beadMass); // J/sigma^2
      CHECK(conversionFactor / (sigmaToM * sigmaToM * 79. * 79.) ==
            Catch::Approx(0.0002486513));
      const double nu =
        nrOfChains / (forceRelaxer2.getNetwork().vol * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      const Eigen::Matrix3d stressTensor = forceRelaxer2.getStressTensor();
      CHECK(
        forceRelaxer2.getPressure() ==
        Catch::Approx(
          (stressTensor(0, 0) + stressTensor(1, 1) + stressTensor(2, 2)) / 3.)
          .epsilon(0.02));
      CHECK(forceRelaxer2.getPressure() * 79. ==
            Catch::Approx(0.1538073302)); // LJ Units [?]
      CHECK(forceRelaxer2.getPressure() * conversionFactor /
              (sigmaToM * sigmaToM * sigmaToM * 79.) ==
            Catch::Approx(
              62085.0438157862)); // shear modulus from the pressure, MPa
      const double expectedNb2 = slope * Nb * beadMass;
      const double bondR02 = forceRelaxer2.getDefaultR0Square();
      const double strandR02 =
        forceRelaxer2.getAverageContourLength() * bondR02;
      const double nb2Correction = strandR02 / expectedNb2;
      const double gammaCorrectionFactor = nb2Correction;
      const double gammaFactor =
        forceRelaxer2.getGammaFactor(bondR02, static_cast<long int>(nrOfChains));
      CHECK(gammaFactor * strandR02 == Catch::Approx(42.6136781099));
      CHECK(gammaFactor * gammaCorrectionFactor * kb * T * nu ==
            Catch::Approx(61308.980768089)); // ANT shear modulus, Pa
      CHECK(gammaFactor * gammaCorrectionFactor ==
            Catch::Approx(0.319449367)); // "correct" gamma factor
      CHECK(forceRelaxer2.getExitReason() == pcm::ExitReason::F_TOLERANCE);
      // TODO: find better, more accurate tests here
      CHECK(forceRelaxer2.getNrOfActiveNodes() > 1);
      CHECK(forceRelaxer2.getNrOfActiveSprings() > 1);
      CHECK(forceRelaxer2.getAverageSpringLength() > 1.0);
      CHECK(forceRelaxer2.getEffectiveFunctionalityOfAtoms().size() ==
            forceRelaxer2.getNrOfNodes());
    } else {
      std::cout << "Skipping large file PDMS MEHP run" << std::endl;
      REQUIRE(true);
    }
  }
}

TEST_CASE("MEHP Force Relaxation returns the crosslinker universe",
          "[analysis][MEHPSimulator][Universe]")
{
  std::cout << "Running test \"MEHP Force Relaxation returns the crosslinker "
               "universe\""
            << std::endl;
  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  SECTION("Melt universe")
  {
    const std::string inputFile =
      suspectedPath + "/structure/melt_83_a_100.structure.out";
    REQUIRE(std::filesystem::exists(inputFile));
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    const std::vector<pu::AtomStyle> atomStyles = { pu::AtomStyle::HYBRID,
                                                    pu::AtomStyle::BOND,
                                                    pu::AtomStyle::EDPD };
    universeSequence.setDataFileAtomStyle(atomStyles);
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    const pe::Universe universe = universeSequence.atIndex(0);

    pcm::MEHPSimulator forceRelaxer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         false);
    CHECK(forceRelaxer.getNrOfNodes() == 83 * 2);
    CHECK(forceRelaxer.getNrOfSprings() == 83);
    forceRelaxer.run();
    CHECK(forceRelaxer.getNrOfActiveSprings() == 0);
    CHECK(forceRelaxer.getResidualNorm() < 1e-9);
    CHECK(forceRelaxer.getEnergy() < 1e-9);

    const Eigen::Matrix3d stressTensor = forceRelaxer.getStressTensor();
    CHECK(stressTensor.isZero(1e-9));

    const size_t nRemoved = forceRelaxer.removeInactiveLinks();
    CHECK(nRemoved == 83 * 2);
    CHECK(forceRelaxer.getNrOfSprings() == 0);
    const pe::Universe crossLinkUniverse = forceRelaxer.getCrosslinkerVerse();
    CHECK(crossLinkUniverse.getNrOfAtoms() == 0);
    CHECK(crossLinkUniverse.getNrOfBonds() == 0);
  }

  SECTION("Non-empty universe")
  {
    const std::string inputFile =
      suspectedPath +
      "/structure/3d-diamond-lattice_3x3x3_a_23_d_3_v_0.structure.out";
    REQUIRE(std::filesystem::exists(inputFile));
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    const pe::Universe universe = universeSequence.atIndex(0);

    pcm::MEHPSimulator forceRelaxer =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         false);
    CHECK(forceRelaxer.getNrOfNodes() == 3 * 3 * 3 * 40);
    const pe::Universe crossLinkUniverse = forceRelaxer.getCrosslinkerVerse();
    CHECK(crossLinkUniverse.getNrOfAtoms() == forceRelaxer.getNrOfNodes());

    std::unordered_map<long int, int> effectiveFunctionality =
      forceRelaxer.getEffectiveFunctionalityOfAtoms();
    REQUIRE(effectiveFunctionality.size() == forceRelaxer.getNrOfNodes());
    for (size_t i = 0; i < universe.getNrOfAtoms(); ++i) {
      const pe::Atom atom =
        universe.getAtomByVertexIdx(static_cast<igraph_integer_t>(i));
      if (atom.getType() == 2) {
        CHECK(effectiveFunctionality[atom.getId()] == 4);
      }
    }

    CHECK(forceRelaxer.getResidual() > 0.);
    CHECK(forceRelaxer.getEnergy() > 0.);
    CHECK(forceRelaxer.getResidual() == forceRelaxer.getResidualNorm());
    const Eigen::Matrix3d stressTensor = forceRelaxer.getStressTensor();
    CHECK_FALSE(stressTensor.isZero(1e-10));
  }
}

TEST_CASE(
  "MEHP Force Relaxation2 runs with Langevin force evaluator and non-network",
  "[analysis][MEHPSimulator][NonGaussianSpringForceEvaluator]")
{
  std::cout << "Running test \"MEHP Force Relaxation2 runs with Langevin force "
               "evaluator and non-network\""
            << std::endl;
  return;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  const std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("3D case")
  {
    const std::string largeInputFile =
      suspectedPath +
      "/structure/xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      const pe::Universe universe2 = universeSeq.atIndex(0);
      const std::shared_ptr<pcm::NonGaussianSpringForceEvaluator>
        nonGaussianForceEvaluator =
          std::make_shared<pcm::NonGaussianSpringForceEvaluator>(1.0, 0.98);
      pcm::MEHPSimulator forceRelaxer2 =
        pcm::MEHPSimulator(universe2,
                           {},
                           2,
                           pcm::EntanglementMode::NONE,
                           nonGaussianForceEvaluator,
                           false,
                           true,
                           false,
                           false);
      REQUIRE(forceRelaxer2.getExitReason() == pcm::ExitReason::UNSET);
      const pcm::MEHPNonlinearSolverConfig nlConfig;
      REQUIRE_NOTHROW(forceRelaxer2.run(nlConfig));
      CHECK(forceRelaxer2.getNrOfIterations() > 1);
    }
  }
}

TEST_CASE("Inverse Langevin test",
          "[analysis][MEHPSimulator][NonGaussianSpringForceEvaluator]")
{
  std::cout << "Running test \"Inverse Langevin test\"" << std::endl;
  // simple test to make sure the inverse langevin approximation is fine
  CHECK_THAT(pcm::langevin_inv(0.01),
             Catch::Matchers::WithinRel(0.0300018, 0.02));
  CHECK_THAT(pcm::langevin_inv(0.1),
             Catch::Matchers::WithinRel(0.301817, 0.02));
  CHECK_THAT(pcm::langevin_inv(0.25),
             Catch::Matchers::WithinRel(0.779897, 0.02));
  CHECK_THAT(pcm::langevin_inv(0.5), Catch::Matchers::WithinRel(1.79676, 0.02));
  CHECK_THAT(pcm::langevin_inv(0.9999),
             Catch::Matchers::WithinRel(10000, 0.02));
  CHECK_THAT(pcm::langevin_inv(1. / 3.),
             Catch::Matchers::WithinRel(1.07456, 0.02));
  CHECK_THAT(pcm::langevin_inv(1.03078 / 3.),
             Catch::Matchers::WithinRel(1.11306, 0.02));
}

TEST_CASE("Manual NonGaussianSpringForceEvaluator gradient test",
          "[analysis][MEHPSimulator][NonGaussianSpringForceEvaluator]")
{
  std::cout
    << "Running test \"Manual NonGaussianSpringForceEvaluator gradient test\""
    << std::endl;
  // setup network
  pcm::MEHPNetwork net;
  net.nrOfLinks = 2;
  net.nrOfNodes = 2;
  net.nrOfStrands = 1;
  net.nrOfSprings = 1;
  net.springIndexA = Eigen::ArrayXi(1);
  net.springIndexA << 0;
  net.springIndexB = Eigen::ArrayXi(1);
  net.springIndexB << 1;
  net.springContourLength = Eigen::VectorXd::Constant(1, 3.0);
  net.meanSpringContourLength = 3.0;
  net.meanStrandContourLength = 3.0;
  // setup evaluator
  const pcm::NonGaussianSpringForceEvaluator forceEvaluatorInstance =
    pcm::NonGaussianSpringForceEvaluator(1.0, 1.0);
  // setup spring vectors and gradient buffer
  Eigen::VectorXd springDistances =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfSprings));
  std::array<double, 6> r{};
  // first, zero spring vector
  CHECK(forceEvaluatorInstance.evaluateNetworkEnergyAndGradient(
          net, springDistances, false, r.data()) == 0.0);
  for (size_t i = 0; i < r.size(); ++i) {
    CHECK(r[i] == 0.0);
  }
  // sign convention: springVec = r_B - r_A + offset, F = +dE/d(springVec),
  // grad[3*a + d] -= F[d], grad[3*b + d] += F[d]
  springDistances[0] = 1.0;
  CHECK_THAT(forceEvaluatorInstance.evaluateNetworkEnergyAndGradient(
               net, springDistances, false, r.data()),
             Catch::Matchers::WithinRel(0.517942, 0.02));
  const Eigen::Vector3d f1 = forceEvaluatorInstance.evaluateSpringForce(
    Eigen::Vector3d(1.0, 0.0, 0.0), 3.0);
  CHECK_THAT(f1[0], Catch::Matchers::WithinRel(1.07456, 0.02));
  CHECK(r[0] == Catch::Approx(-f1[0]));
  CHECK_THAT(-r[0], Catch::Matchers::WithinRel(1.07456, 0.02));
  CHECK(r[1] == 0.0);
  CHECK(r[2] == 0.0);
  CHECK(r[3] == Catch::Approx(-r[0]));
  CHECK(r[4] == Catch::Approx(-r[1]));
  CHECK(r[5] == Catch::Approx(-r[2]));

  // and again with non-zero x and y
  springDistances[1] = -1.0;
  const double rDist = std::sqrt(2.0);
  CHECK_THAT(forceEvaluatorInstance.evaluateNetworkEnergyAndGradient(
               net, springDistances, false, r.data()),
             Catch::Matchers::WithinRel(1.07797, 0.02));
  CHECK(r[0] == Catch::Approx(-1.6542 * 1.0 / rDist).epsilon(0.02));
  CHECK(r[1] == Catch::Approx(1.6542 * 1.0 / rDist).epsilon(0.02));
  CHECK(r[2] == 0.0);
  CHECK(r[3] == Catch::Approx(-r[0]));
  CHECK(r[4] == Catch::Approx(-r[1]));
  CHECK(r[5] == Catch::Approx(-r[2]));
  // and a 3D stretched spring
  springDistances[2] = 0.25;
  CHECK_THAT(forceEvaluatorInstance.evaluateNetworkEnergyAndGradient(
               net, springDistances, false, r.data()),
             Catch::Matchers::WithinRel(1.11463, 0.02));
  const double oneOverRDist = 1.0 / std::sqrt(2.0 + 0.25 * 0.25);
  CHECK(r[0] == Catch::Approx(-1.68968 * oneOverRDist).epsilon(0.02));
  CHECK(r[1] == Catch::Approx(1.68968 * oneOverRDist).epsilon(0.02));
  CHECK(r[2] == Catch::Approx(-1.68968 * 0.25 * oneOverRDist).epsilon(0.02));
  CHECK(r[3] == Catch::Approx(-r[0]));
  CHECK(r[4] == Catch::Approx(-r[1]));
  CHECK(r[5] == Catch::Approx(-r[2]));
}

TEST_CASE("Force Relaxation free chains collapse",
          "[analysis][MEHPSimulator][NonGaussianSpringForceEvaluator]["
          "HookeanForceEvaluator]")
{
  std::cout << "Running test \"Force Relaxation free chains collapse\""
            << std::endl;
  const size_t nrOfBeads = 30;
  const size_t nrOfBeadsPerChain = 3;
  pe::Universe universe =
    pe::Universe(nrOfBeads * 10.0, nrOfBeads * 10.0, nrOfBeads * 10.0);
  std::vector<double> xPositions;
  xPositions.reserve(nrOfBeads);
  std::vector<double> yPositions;
  yPositions.reserve(nrOfBeads);
  std::vector<double> zPositions;
  zPositions.reserve(nrOfBeads);
  std::vector<long int> atomIds;
  atomIds.reserve(nrOfBeads);
  std::vector<int> atomTypes;
  atomTypes.reserve(nrOfBeads);
  std::vector<int> zeroInts;
  zeroInts.reserve(nrOfBeads);
  std::vector<long int> bondFrom;
  bondFrom.reserve(nrOfBeads - 1);
  std::vector<long int> bondTo;
  bondTo.reserve(nrOfBeads - 1);
  const double offset = 10.0;
  for (int i = 0; i < static_cast<int>(nrOfBeads); ++i) {
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i);
    atomTypes.push_back(static_cast<size_t>(i) % nrOfBeadsPerChain == 0 ? 2
                                                                        : 1);
    zeroInts.push_back(0);
    if (i > 0) {
      bondFrom.push_back(i - 1);
      bondTo.push_back(i);
    }
  }
  universe.addAtoms(atomIds,
                    atomTypes,
                    xPositions,
                    yPositions,
                    zPositions,
                    zeroInts,
                    zeroInts,
                    zeroInts);
  universe.addBonds(bondFrom, bondTo);
  REQUIRE(universe.getNrOfAtoms() == nrOfBeads);
  REQUIRE(universe.getNrOfBonds() == nrOfBeads - 1);
  // add a primary and secondary loop to check that it collapses
  universe.addBonds({ 0, 0 }, { static_cast<long>(nrOfBeadsPerChain), 0 });

  CHECK(universe.getChainsWithCrosslinker(2).size() == 12);
  CHECK(universe.getVertexDegree(
          universe.getIdxByAtomId(static_cast<long>(nrOfBeadsPerChain))) == 3);
  CHECK(universe.getVertexDegree(0) == 4);

  // check that the Hookean evaluator collapses free chains
  const std::shared_ptr<pcm::HookeanForceEvaluator>
    simpleSpringForceEvaluator =
      std::make_shared<pcm::HookeanForceEvaluator>();
  pcm::MEHPSimulator forceRelaxerSimpleSpring =
    pcm::MEHPSimulator(universe,
                       {},
                       2,
                       pcm::EntanglementMode::NONE,
                       simpleSpringForceEvaluator,
                       false,
                       true,
                       false,
                       false);

  CHECK(forceRelaxerSimpleSpring.getNrOfSprings() == 12);

  const pcm::MEHPNonlinearSolverConfig nlConfig;
  REQUIRE_NOTHROW(forceRelaxerSimpleSpring.run(nlConfig));
  REQUIRE(forceRelaxerSimpleSpring.getNrOfIterations() > 0);
  CHECK(forceRelaxerSimpleSpring.getNrOfActiveSprings() == 0);
  CHECK(forceRelaxerSimpleSpring.getAverageSpringLength() >= 0.0);
  CHECK(forceRelaxerSimpleSpring.getAverageSpringLength() <= 1e-6);
}

TEST_CASE("Fully active chains are fully active",
          "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"Fully active chains are fully active\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  const std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Relaxation 3D case")
  {
    // perfect diamond network = fully connected =>
    // maximum is at perfect crystal structure -> must be all active.
    const std::string inputFile =
      suspectedPath +
      "/structure/"
      "3d-diamond-lattice_10x10x10_a_3_d_0.85_v_0.V-fixed.structure.out";
    if (std::filesystem::exists(inputFile)) {
      std::cout << "Reading file " << inputFile << std::endl;
      universeSeq.initializeFromDataSequence({ { inputFile } });
      const pe::Universe universe = universeSeq.atIndex(0);
      std::cout << "Read file " << inputFile << std::endl;

      pcm::MEHPSimulator forceRelaxer =
        pcm::MEHPSimulator(universe,
                           {},
                           2,
                           pcm::EntanglementMode::NONE,
                           nullptr,
                           false,
                           true,
                           false,
                           false);
      const pcm::MEHPNonlinearSolverConfig nlConfig;
      REQUIRE_NOTHROW(forceRelaxer.run(nlConfig));
      REQUIRE(forceRelaxer.getNrOfIterations() > 0);
      CHECK(forceRelaxer.getExitReason() == pcm::ExitReason::F_TOLERANCE);
      CHECK(forceRelaxer.getNrOfActiveSprings() ==
            forceRelaxer.getNrOfSprings());
    }
  }
}

TEST_CASE("MEHP Bugfixes and Regression Tests",
          "[analysis][MEHPSimulator][bugfixes]")
{
  SECTION("1. Test langevin_inv symmetry and restorative force")
  {
    const std::vector<double> testValues = { 0.0, 0.01, 0.1, 0.5, 0.8, 0.99,
                                             1.5 };
    for (const double val : testValues) {
      const double pos = pcm::langevin_inv(val);
      const double neg = pcm::langevin_inv(-val);
      CHECK(neg == Catch::Approx(-pos));
      if (val > 0.0) {
        CHECK(pos > 0.0);
        CHECK(neg < 0.0);
      }
    }
  }

  SECTION(
    "2. Test evaluateNetworkEnergyAndGradient and csch at r=0: returns 0 "
    "and no NaN")
  {
    pcm::MEHPNetwork net;
    net.nrOfLinks = 2;
    net.nrOfNodes = 2;
    net.nrOfStrands = 1;
    net.nrOfSprings = 1;
    net.springIndexA = Eigen::ArrayXi(1);
    net.springIndexA << 0;
    net.springIndexB = Eigen::ArrayXi(1);
    net.springIndexB << 1;
    net.springContourLength = Eigen::VectorXd::Constant(1, 5.0);
    net.meanSpringContourLength = 5.0;
    net.meanStrandContourLength = 5.0;

    const pcm::NonGaussianSpringForceEvaluator evaluator(1.0, 1.0);

    const Eigen::VectorXd springDistances = Eigen::VectorXd::Zero(3);
    std::array<double, 6> grad = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    const double energy = evaluator.evaluateNetworkEnergyAndGradient(
      net, springDistances, false, grad.data());
    CHECK(energy == 0.0);
    for (size_t i = 0; i < grad.size(); ++i) {
      CHECK_FALSE(std::isnan(grad[i]));
      CHECK(grad[i] == 0.0);
    }

    const Eigen::Vector3d s = Eigen::Vector3d::Zero();
    const double stress = evaluator.evaluateStressContribution(s, 0, 0, 5.0);
    CHECK_FALSE(std::isnan(stress));
    CHECK(stress == 0.0);
  }

  SECTION("3. Test NonGaussian evaluateStressContribution incorporates "
          "non-linear tension")
  {
    const pcm::NonGaussianSpringForceEvaluator nonGaussianEval(1.0, 1.0);
    const pcm::HookeanForceEvaluator gaussianEval(1.0);

    // Highly stretched spring: r = 9.0 (close to contour length N*l = 10.0)
    const Eigen::Vector3d s(9.0, 0.0, 0.0);
    const double nonGaussianStress =
      nonGaussianEval.evaluateStressContribution(s, 0, 0, 10.0);
    const double gaussianStress =
      gaussianEval.evaluateStressContribution(s, 0, 0, 10.0);

    // Non-linear Langevin tension diverges near contour length, so
    // nonGaussianStress >> gaussianStress
    CHECK_FALSE(std::isnan(nonGaussianStress));
    CHECK(nonGaussianStress > gaussianStress * 2.0);
  }

  SECTION("4. Test MEHPSimulator move and copy semantics")
  {
    STATIC_CHECK(std::is_copy_constructible_v<pcm::MEHPSimulator>);
    STATIC_CHECK(std::is_copy_assignable_v<pcm::MEHPSimulator>);
    STATIC_CHECK(std::is_move_constructible_v<pcm::MEHPSimulator>);
    STATIC_CHECK(std::is_move_assignable_v<pcm::MEHPSimulator>);

    const pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
    pcm::MEHPSimulator r1(universe);
    pcm::MEHPSimulator rCopy(r1);
    CHECK(rCopy.getNrOfIterations() == 0);
    pcm::MEHPSimulator r2 = std::move(r1);
    CHECK(r2.getNrOfIterations() == 0);
  }
}
