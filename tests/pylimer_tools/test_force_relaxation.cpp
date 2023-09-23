#include "../../src/pylimer_tools_cpp/calc/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/calc/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <vector>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::calc::mehp;

void
testGradient(pcm::MEHPForceEvaluator* forceEvaluator)
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
  pcm::MEHPForceRelaxation forceRelaxer2 =
    pcm::MEHPForceRelaxation(universe, 2);
  // assemble network
  pcm::Network net;
  Eigen::VectorXd coordinates = Eigen::VectorXd(12);
  coordinates << 0.0, 0.0, 0.0, // 1
    1.0, 0.0, 0.0,              // 2
    2.0, 0.0, 0.0,              // 3
    3.0, 0.0, 0.0;              // 4
  net.coordinates = coordinates;
  Eigen::ArrayXi springCoordinateIndexA = Eigen::ArrayXi::Zero(12);
  net.springCoordinateIndexA = springCoordinateIndexA;
  Eigen::ArrayXi springCoordinateIndexB = Eigen::ArrayXi::Zero(12);
  net.springCoordinateIndexB = springCoordinateIndexB;
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
  net.nrOfNodes = 4;
  net.nrOfSprings = 4;
  net.vol = universe.getVolume();
  net.L[0] = 4.0;
  net.L[1] = 4.0;
  net.L[2] = 4.0;
  net.nrOfLoops = 1;
  net.springsContourLength = Eigen::VectorXd::Constant(4, 1.0);
  net.meanSpringContourLength = 1.0;
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
    (coordinatesSpringEndA - coordinatesSpringEndB);
  for (size_t i = 0; i < 12; ++i) {
    REQUIRE(coordinatesSpringEndA[i] == expectedCoordsA[i]);
    REQUIRE(coordinatesSpringEndB[i] == expectedCoordsB[i]);
    REQUIRE(springDistances[i] == expectedCoordsA[i] - expectedCoordsB[i]);
  }
  // setup gradient and coordinates
  const double h = 1.e-5;
  double grad[12];
  double x[12];
  for (int i = 0; i < 12; ++i) {
    grad[i] = 0.0;
    x[i] = 0.0;
  }
  forceEvaluator->setIs2D(false);
  forceEvaluator->setNetwork(net);
  forceEvaluator->prepareForEvaluations();
  // actual computation to test gradient
  for (size_t i = 0; i < 12; ++i) {
    // std::cout << "MEHP Gradient Test coordinate " << i << std::endl;
    // evaluate gradient
    double f = forceEvaluator->evaluateForceSetGradient(12, x, grad, nullptr);
    if (i % 3 != 0) {
      // in x and y direction, we expect no spring distance -> 0 gradient
      CHECK(grad[i] == Catch::Approx(0.0));
    } else {
      // test finite difference vs. gradient
      x[i] = -h;
      double fm =
        forceEvaluator->evaluateForceSetGradient(12, x, nullptr, nullptr);
      x[i] = h;
      double fp =
        forceEvaluator->evaluateForceSetGradient(12, x, nullptr, nullptr);
      x[i] = 0.0; // reset
      // std::cout << i << " " << fm << " " << fp << " " << f << std::endl;
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
          "[analysis][MEHPForceRelaxation][SimpleSpringMEHPForceEvaluator]")
{
  SECTION("Test SimpleSpringMEHPForceEvaluator force gradient")
  {
    pcm::SimpleSpringMEHPForceEvaluator forceEvaluatorInstance =
      pcm::SimpleSpringMEHPForceEvaluator(1.0);
    testGradient(&forceEvaluatorInstance);
  }

  SECTION("Test NonGaussianSpringForceEvaluator force gradient")
  {
    pcm::NonGaussianSpringForceEvaluator forceEvaluatorInstance =
      pcm::NonGaussianSpringForceEvaluator(1.0, 2.0, 1.0);
    testGradient(&forceEvaluatorInstance);
  }
}

TEST_CASE(
  "MEHP Force Relaxation2 runs",
  "[analysis][MEHPForceRelaxation][SimpleSpringMEHPForceEvaluator][long]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("MEHP Force Relaxation2 3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);

      // BENCHMARK_ADVANCED("MEHP LD_MMA " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceRelaxation forceRelaxer3 =
      //     pcm::MEHPForceRelaxation(universe2, 2);
      //   meter.measure([&forceRelaxer3] {
      //     forceRelaxer3.runForceRelaxation("LD_MMA");
      //     return forceRelaxer3.getNrOfIterations();
      //   });
      // };
      // BENCHMARK_ADVANCED("MEHP LD_LBFGS " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceRelaxation forceRelaxer3 =
      //     pcm::MEHPForceRelaxation(universe2, 2);
      //   meter.measure([&forceRelaxer3] {
      //     forceRelaxer3.runForceRelaxation("LD_LBFGS");
      //     return forceRelaxer3.getNrOfIterations();
      //   });
      // };

      double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      pcm::MEHPForceRelaxation forceRelaxer2 =
        pcm::MEHPForceRelaxation(universe2, 2);
      REQUIRE(forceRelaxer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE(forceRelaxer2.getNrOfIterations() == 0);
      REQUIRE(forceRelaxer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
      CHECK(forceRelaxer2.getVolume() ==
            Catch::Approx(97.383096 * 97.383096 * 97.383096));
      // initial system values
      CHECK(forceRelaxer2.getPressure() ==
            Catch::Approx(0.39911682390778536 / 79.));
      CHECK(forceRelaxer2.getForce() * 79 == Catch::Approx(552894.1903005145));
      CHECK(forceRelaxer2.getAverageContourLength() == Catch::Approx(79.0));
      Eigen::VectorXd contourLengths = forceRelaxer2.getSpringContourLength();
      for (int i = 0; i < contourLengths.size(); ++i) {
        CHECK(contourLengths[i] ==
              Catch::Approx(forceRelaxer2.getAverageContourLength()));
      }
      CHECK(forceRelaxer2.getResidualNorm() * 79 ==
            Catch::Approx(1457.465048151));
      REQUIRE_NOTHROW(forceRelaxer2.runForceRelaxation());
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
      double Nb = 80.; // nr of bonds per strand
      double conversionFactor =
        3. * kb * T / (slope * beadMass * (Nb / 79.)); // J/sigma^2
      CHECK((conversionFactor / 79.) / (sigmaToM * sigmaToM) ==
            Catch::Approx(0.000245543));
      double nu =
        nrOfChains / (forceRelaxer2.getVolume() * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      CHECK(forceRelaxer2.getPressure() ==
            Catch::Approx(0.153806 / 79.)); // LJ Units [?]
      CHECK(forceRelaxer2.getPressure() * conversionFactor /
              (sigmaToM * sigmaToM * sigmaToM) ==
            Catch::Approx(61308.3)); // shear modulus from the pressure, MPa
      double nrOfChainCorrection =
        (forceRelaxer2.getDefaultNrOfChains() / nrOfChains);
      double expectedNb2 = slope * Nb * beadMass;
      double nb2Correction =
        (forceRelaxer2.getDefaultR0Square() / (expectedNb2));
      double gammaCorrectionFactor = nrOfChainCorrection * nb2Correction;
      CHECK(
        forceRelaxer2.getGammaFactor() * nrOfChainCorrection *
          forceRelaxer2.getDefaultR0Square() ==
        Catch::Approx(42.6132)); // as from conversion-less Mathematica script
      CHECK(forceRelaxer2.getGammaFactor() * gammaCorrectionFactor * kb * T *
              nu ==
            Catch::Approx(61308.3)); // ANT shear modulus, Pa
      CHECK(forceRelaxer2.getGammaFactor() * gammaCorrectionFactor ==
            Catch::Approx(0.319446)); // "correct" gamma factor
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
      { { suspectedPath + "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
                          "333_2d_t_7500001.structure.out" } });
    REQUIRE(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPForceRelaxation forceRelaxer =
      pcm::MEHPForceRelaxation(universe, 2, true);
    REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::UNSET);
    REQUIRE(forceRelaxer.getNrOfIterations() == 0);
    REQUIRE(forceRelaxer.getVolume() == Catch::Approx(universe.getVolume()));
    forceRelaxer.runForceRelaxation("LD_LBFGS", 5);
    REQUIRE(forceRelaxer.getNrOfNodes() != universe.getNrOfAtoms());
    REQUIRE(forceRelaxer.getNrOfIterations() <= 5);
    REQUIRE(forceRelaxer.getNrOfIterations() >= 1);
    REQUIRE(universe.getAtomsOfType(2).size() == 7200);
    REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::MAX_STEPS);
    // run again, this time fully
    pcm::MEHPForceRelaxation forceRelaxer2 =
      pcm::MEHPForceRelaxation(universe, 2, true);
    forceRelaxer2.runForceRelaxation("LD_LBFGS", 10000, 1e-5, 1e-18);
    REQUIRE(forceRelaxer2.getNrOfIterations() > 5);
    CHECK(forceRelaxer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceRelaxer2.getGammaFactor(25, forceRelaxer2.getNrOfSprings()) ==
          Catch::Approx(1. / 3.).epsilon(0.001));
    auto stressTensor = forceRelaxer2.getStressTensor();
    CHECK(forceRelaxer2.getPressure() ==
          Catch::Approx(
            (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.)
            .epsilon(0.02));
    CHECK(forceRelaxer2.getResidualNorm() > 0.0);
    CHECK(forceRelaxer2.getForce() > 0.0);
    // TODO: find better, more accurate tests here
    CHECK(forceRelaxer2.getNrOfActiveNodes() > 1);
    CHECK(forceRelaxer2.getNrOfActiveSprings() > 1);
    CHECK(forceRelaxer2.getAverageSpringLength() > 1.0);
    CHECK(forceRelaxer2.getEffectiveFunctionalityOfAtoms().size() ==
          forceRelaxer2.getNrOfNodes());
    CHECK(forceRelaxer2.getSolubleWeightFraction() > 0.);

    pe::Universe universe3 = forceRelaxer2.getCrosslinkerVerse();
    CHECK(universe3.getNrOfAtoms() == forceRelaxer2.getNrOfNodes());
    CHECK(universe3.getNrOfBonds() == forceRelaxer2.getNrOfSprings());
    CHECK(universe3.getAtomsOfType(2).size() == universe3.getNrOfAtoms());

    // try out different algorithms
    std::vector<std::string> algorithms = { "LD_MMA",
                                            // "LD_TNEWTON_PRECOND_RESTART",
                                            // "GD_STOGO",
                                            "LD_SLSQP",
                                            "GN_DIRECT" };

    // for (std::string algorithm : algorithms) {
    //   pcm::MEHPForceRelaxation forceRelaxerN =
    //     pcm::MEHPForceRelaxation(universe, 2);
    //   std::cout << "Testing algorithm " << algorithm << std::endl;
    //   auto start = std::chrono::high_resolution_clock::now();
    //   forceRelaxerN.runForceRelaxation(true, 15, algorithm.c_str(), 10000);
    //   auto stop = std::chrono::high_resolution_clock::now();
    //   CHECK(forceRelaxerN.getGammaEq() ==
    //         Catch::Approx(forceRelaxer2.getGammaEq()));
    //   CHECK(forceRelaxerN.getFinalPressure() ==
    //         Catch::Approx(forceRelaxer2.getFinalPressure()));
    //   auto duration = duration_cast<std::chrono::microseconds>(stop - start);
    //   std::cout << "Took: " << duration.count() << std::endl;
    // }
  }

  SECTION("Simple melt case")
  {
    std::string dataFile = suspectedPath + "lammps_data_file.out";
    universeSeq.initializeFromDataSequence({ { dataFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPForceRelaxation forceRelaxer =
      pcm::MEHPForceRelaxation(universe, 2, true);
    // not a proper network -> 0 springs
    CHECK(forceRelaxer.getNrOfSprings() == 0);
    CHECK(forceRelaxer.getSolubleWeightFraction() == 1.);
  }
}

TEST_CASE(
  "MEHP Force Relaxation2 runs with non-gaussian force evaluator",
  "[analysis][MEHPForceRelaxation][NonGaussianSpringForceEvaluator][long]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("3D case")
  {
    // TODO: correct values (& forces?)
    std::string largeInputFile =
      suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);

      // BENCHMARK_ADVANCED("MEHP LD_MMA " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceRelaxation forceRelaxer3 =
      //     pcm::MEHPForceRelaxation(universe2, 2);
      //   meter.measure([&forceRelaxer3] {
      //     forceRelaxer3.runForceRelaxation("LD_MMA");
      //     return forceRelaxer3.getNrOfIterations();
      //   });
      // };
      // BENCHMARK_ADVANCED("MEHP LD_LBFGS " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceRelaxation forceRelaxer3 =
      //     pcm::MEHPForceRelaxation(universe2, 2);
      //   meter.measure([&forceRelaxer3] {
      //     forceRelaxer3.runForceRelaxation("LD_LBFGS");
      //     return forceRelaxer3.getNrOfIterations();
      //   });
      // };

      double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      pcm::NonGaussianSpringForceEvaluator nonGaussianForceEvaluator =
        pcm::NonGaussianSpringForceEvaluator(1.0, 78., 0.98);
      pcm::MEHPForceRelaxation forceRelaxer2 = pcm::MEHPForceRelaxation(
        universe2, 2, false, &nonGaussianForceEvaluator);
      REQUIRE(forceRelaxer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE(forceRelaxer2.getNrOfIterations() == 0);
      REQUIRE(forceRelaxer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
      CHECK(forceRelaxer2.getVolume() ==
            Catch::Approx(97.383096 * 97.383096 * 97.383096));
      // initial system values
      CHECK(forceRelaxer2.getForce() ==
            Catch::Approx(22103.6441026747).epsilon(1e-5));
      CHECK(forceRelaxer2.getResidualNorm() ==
            Catch::Approx(58.7568113327).epsilon(1e-5));
      CHECK(forceRelaxer2.getPressure() *
              forceRelaxer2.getNetwork().meanSpringContourLength ==
            Catch::Approx(0.39911682390778536));
      CHECK(forceRelaxer2.getAverageContourLength() == Catch::Approx(79.0));
      Eigen::VectorXd contourLengths = forceRelaxer2.getSpringContourLength();
      for (int i = 0; i < contourLengths.size(); ++i) {
        CHECK(contourLengths[i] ==
              Catch::Approx(forceRelaxer2.getAverageContourLength()));
      }
      REQUIRE_NOTHROW(forceRelaxer2.runForceRelaxation("LD_MMA"));
      // As long as gradient is unclear: gradient free methods, e.g.:
      // "LN_SBPLX", "LN_BOBYQA", "LN_NELDERMEAD",
      // "LN_COBYLA", "LN_NEWUOA_BOUND"
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
      double nu =
        nrOfChains / (forceRelaxer2.getVolume() * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      auto stressTensor = forceRelaxer2.getStressTensor();
      CHECK(
        forceRelaxer2.getPressure() ==
        Catch::Approx(
          (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.)
          .epsilon(0.02));
      CHECK(forceRelaxer2.getPressure() * 79. ==
            Catch::Approx(0.1538073308)); // LJ Units [?]
      CHECK(forceRelaxer2.getPressure() * conversionFactor /
              (sigmaToM * sigmaToM * sigmaToM * 79.) ==
            Catch::Approx(61308.9809826224 * 80. /
                          79.)); // shear modulus from the pressure, MPa
      double nrOfChainCorrection =
        (forceRelaxer2.getDefaultNrOfChains() / nrOfChains);
      double expectedNb2 = slope * Nb * beadMass;
      double nb2Correction =
        (forceRelaxer2.getDefaultR0Square() / (expectedNb2));
      double gammaCorrectionFactor = nrOfChainCorrection * nb2Correction;
      CHECK(forceRelaxer2.getGammaFactor() * nrOfChainCorrection *
              forceRelaxer2.getDefaultR0Square() ==
            Catch::Approx(42.613678259));
      CHECK(forceRelaxer2.getGammaFactor() * gammaCorrectionFactor * kb * T *
              nu ==
            Catch::Approx(61308.9809826224)); // ANT shear modulus, Pa
      CHECK(forceRelaxer2.getGammaFactor() * gammaCorrectionFactor ==
            Catch::Approx(0.3194493682)); // "correct" gamma factor
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

TEST_CASE(
  "MEHP Force Relaxation2 runs with Langevin force evaluator and non-network",
  "[analysis][MEHPForceRelaxation][NonGaussianSpringForceEvaluator]")
{
  return;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      pcm::NonGaussianSpringForceEvaluator nonGaussianForceEvaluator =
        pcm::NonGaussianSpringForceEvaluator(1.0, 79, 0.98);
      pcm::MEHPForceRelaxation forceRelaxer2 = pcm::MEHPForceRelaxation(
        universe2, 2, false, &nonGaussianForceEvaluator);
      REQUIRE(forceRelaxer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE_NOTHROW(forceRelaxer2.runForceRelaxation("LD_MMA"));
      CHECK(forceRelaxer2.getNrOfIterations() > 1);
    }
  }
}

TEST_CASE("Inverse Langevin test",
          "[analysis][MEHPForceRelaxation][NonGaussianSpringForceEvaluator]")
{
  // simple test to make sure the inverse langevin approximation is fine
  CHECK(pcm::langevin_inv(0.01) == Catch::Approx(0.0300018).epsilon(0.02));
  CHECK(pcm::langevin_inv(0.1) == Catch::Approx(0.301817).epsilon(0.02));
  CHECK(pcm::langevin_inv(0.25) == Catch::Approx(0.779897).epsilon(0.02));
  CHECK(pcm::langevin_inv(0.5) == Catch::Approx(1.79676).epsilon(0.02));
  CHECK(pcm::langevin_inv(0.9999) == Catch::Approx(10000).epsilon(0.02));
  CHECK(pcm::langevin_inv(1. / 3.) == Catch::Approx(1.07456).epsilon(0.02));
  CHECK(pcm::langevin_inv(1.03078 / 3.) ==
        Catch::Approx(1.11306).epsilon(0.02));
}

TEST_CASE("Manual NonGaussianSpringForceEvaluator gradient test",
          "[analysis][MEHPForceRelaxation][NonGaussianSpringForceEvaluator]")
{
  // return;
  // setup network
  pcm::Network net;
  net.nrOfSprings = 1;
  net.nrOfNodes = 2;
  net.springIndexA = Eigen::VectorXi(1);
  net.springIndexA << 0;
  net.springIndexB = Eigen::VectorXi(1);
  net.springIndexB << 1;
  net.springsContourLength = Eigen::VectorXd::Constant(1, 3.0);
  net.meanSpringContourLength = 3.0;
  // setup evaluator
  pcm::NonGaussianSpringForceEvaluator forceEvaluatorInstance =
    pcm::NonGaussianSpringForceEvaluator(1.0, 3.0, 1.0);
  forceEvaluatorInstance.setIs2D(false);
  forceEvaluatorInstance.setNetwork(net);
  forceEvaluatorInstance.prepareForEvaluations();
  // setup other variables
  Eigen::VectorXd springDistances = Eigen::VectorXd::Zero(3 * net.nrOfSprings);
  Eigen::VectorXd u = Eigen::VectorXd::Zero(3 * net.nrOfNodes);
  double* r = new double[3 * net.nrOfNodes];
  // actually check values
  // first, the zero positions
  CHECK(forceEvaluatorInstance.evaluateForceSetGradient(
          net.nrOfNodes * 3, springDistances, u, r) == 0.0);
  for (int i = 0; i < net.nrOfNodes * 3; ++i) {
    CHECK(r[i] == 0.0);
  }
  // then, some values as compared to what is obtained from a Mathematica script
  springDistances[0] = 1.0;
  CHECK(forceEvaluatorInstance.evaluateForceSetGradient(
          net.nrOfNodes * 3, springDistances, u, r) ==
        Catch::Approx(0.517942).epsilon(0.02));
  CHECK(r[0] == Catch::Approx(1.07456).epsilon(0.02));
  CHECK(r[1] == 0.0);
  CHECK(r[2] == 0.0);
  CHECK(r[3] == Catch::Approx(-r[0]));
  CHECK(r[4] == Catch::Approx(-r[1]));
  CHECK(r[5] == Catch::Approx(-r[2]));

  // and again some others
  springDistances[1] = -1.0;
  double rDist = std::sqrt(2.0);
  CHECK(forceEvaluatorInstance.evaluateForceSetGradient(
          net.nrOfNodes * 3, springDistances, u, r) ==
        Catch::Approx(1.07797).epsilon(0.02));
  CHECK(r[0] == Catch::Approx(1.6542 * 1.0 / rDist).epsilon(0.02));
  CHECK(r[1] == Catch::Approx(-1.6542 * 1.0 / rDist).epsilon(0.02));
  CHECK(r[2] == 0.0);
  CHECK(r[3] == Catch::Approx(-r[0]));
  CHECK(r[4] == Catch::Approx(-r[1]));
  CHECK(r[5] == Catch::Approx(-r[2]));
  // and a final one
  springDistances[2] = 0.25;
  CHECK(forceEvaluatorInstance.evaluateForceSetGradient(
          net.nrOfNodes * 3, springDistances, u, r) ==
        Catch::Approx(1.11463).epsilon(0.02));
  double OneOverRDist = 1.0 / std::sqrt(2.0 + 0.25 * 0.25);
  CHECK(r[0] == Catch::Approx(1.68968 * OneOverRDist).epsilon(0.02));
  CHECK(r[1] == Catch::Approx(-1.68968 * OneOverRDist).epsilon(0.02));
  CHECK(r[2] == Catch::Approx(1.68968 * 0.25 * OneOverRDist).epsilon(0.02));
  CHECK(r[3] == Catch::Approx(-r[0]));
  CHECK(r[4] == Catch::Approx(-r[1]));
  CHECK(r[5] == Catch::Approx(-r[2]));

  // cleanup
  delete[] (r);
}

TEST_CASE("Free chains collapse",
          "[analysis][MEHPForceRelaxation][NonGaussianSpringForceEvaluator]["
          "SimpleSpringMEHPForceEvaluator]")
{
  size_t nrOfBeads = 30;
  size_t nrOfBeadsPerChain = 3;
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
  double offset = 10.0;
  for (int i = 0; i < nrOfBeads; ++i) {
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i);
    atomTypes.push_back(i % nrOfBeadsPerChain == 0 ? 2 : 1);
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

  // now, check for every force evaluator, that the maximum entropy is when all
  // these beads overlap first, the gaussian spring one
  pcm::SimpleSpringMEHPForceEvaluator simpleSpringForceEvaluator =
    pcm::SimpleSpringMEHPForceEvaluator();
  pcm::MEHPForceRelaxation forceRelaxerSimpleSpring =
    pcm::MEHPForceRelaxation(universe, 2, false, &simpleSpringForceEvaluator);
  REQUIRE_NOTHROW(forceRelaxerSimpleSpring.runForceRelaxation());
  REQUIRE(forceRelaxerSimpleSpring.getNrOfIterations() > 0);
  CHECK(forceRelaxerSimpleSpring.getNrOfActiveSprings() == 0);
  // CHECK(forceRelaxerSimpleSpring.getAverageSpringLength() ==
  //       Catch::Approx(0.0));
  CHECK(forceRelaxerSimpleSpring.getAverageSpringLength() >= 0.0);
  CHECK(forceRelaxerSimpleSpring.getAverageSpringLength() <= 1e-6);

  // then, the non-gaussian one
  pcm::NonGaussianSpringForceEvaluator langevinForceEvaluator =
    pcm::NonGaussianSpringForceEvaluator(1.0, nrOfBeadsPerChain * 2, 1.0);
  pcm::MEHPForceRelaxation forceRelaxerLangevin =
    pcm::MEHPForceRelaxation(universe, 2, false, &langevinForceEvaluator);
  pe::Universe resultingUniverse0 = forceRelaxerLangevin.getCrosslinkerVerse();
  // auto distances0 = resultingUniverse0.computeBondLengths();
  // for (auto i : distances0) {
  //   std::cout << i << std::endl;
  // }
  // std::cout << forceRelaxerLangevin.getNrOfIterations() << ", "
  //           << forceRelaxerLangevin.getForce() << ", "
  //           << forceRelaxerLangevin.getResidualNorm() << std::endl;

  forceRelaxerLangevin.runForceRelaxation(
    "LD_MMA", 500000, 1e-15, 1e-19); // "LD_SLSQP", "LD_MMA"
  CHECK(forceRelaxerLangevin.getNrOfActiveSprings() == 0);
  // CHECK(forceRelaxerLangevin.getAverageSpringLength() == Catch::Approx(0.0));
  CHECK(forceRelaxerLangevin.getAverageSpringLength() >= 0.0);
  CHECK(forceRelaxerLangevin.getAverageSpringLength() <= 1e-6);
  REQUIRE(forceRelaxerLangevin.getNrOfIterations() > 0);
  CHECK(forceRelaxerLangevin.getExitReason() == pcm::ExitReason::X_TOLERANCE);

  pe::Universe resultingUniverse = forceRelaxerLangevin.getCrosslinkerVerse();
  // auto distances = resultingUniverse.computeBondLengths();
  // for (auto i : distances) {
  //   std::cout << i << std::endl;
  // }
  // auto residuals = forceRelaxerLangevin.getResiduals();
  // for (auto i : residuals) {
  //   std::cout << i << " ";
  // }
  // std::cout << std::endl;
  // std::cout << forceRelaxerLangevin.getNrOfIterations() << ", "
  //           << forceRelaxerLangevin.getForce() << ", "
  //           << forceRelaxerLangevin.getResidualNorm() << std::endl;
  REQUIRE(forceRelaxerLangevin.getForce() >=
          forceRelaxerSimpleSpring.getForce());
  REQUIRE(forceRelaxerLangevin.getResidualNorm() >=
          forceRelaxerSimpleSpring.getResidualNorm());

  std::array<std::array<double, 3>, 3> stressTensorSimpleSpring =
    forceRelaxerSimpleSpring.getStressTensor();
  std::array<std::array<double, 3>, 3> stressTensorLangevin =
    forceRelaxerLangevin.getStressTensor();
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      CHECK(stressTensorLangevin[i][j] + 1e-5 ==
            Catch::Approx(stressTensorSimpleSpring[i][j] + 1e-5).margin(5e-7));
    }
  }
}
