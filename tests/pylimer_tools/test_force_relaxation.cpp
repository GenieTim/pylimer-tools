#include "../../src/pylimer_tools_cpp/calc/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
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

TEST_CASE("MEHP Force Relaxation2 computes correct gradients",
          "[analysis][MEHPForceRelaxation]")
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
  net.averageSpringLength = 1.0;
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
  // actual computation to test gradient
  for (size_t i = 0; i < 12; ++i) {
    // std::cout << "MEHP Gradient Test coordinate " << i << std::endl;
    // evaluate gradient
    double f = forceRelaxer2.evaluateForceSetGradient(
      &net, 1.0, false, 12, x, grad, NULL);
    if (i % 3 != 0) {
      // in x and y direction, we expect no spring distance -> 0 gradient
      REQUIRE(grad[i] == Catch::Approx(0.0));
    } else {
      // test finite difference vs. gradient
      x[i] = -h;
      double fm = forceRelaxer2.evaluateForceSetGradient(
        &net, 1.0, false, 12, x, NULL, NULL);
      x[i] = h;
      double fp = forceRelaxer2.evaluateForceSetGradient(
        &net, 1.0, false, 12, x, NULL, NULL);
      x[i] = 0.0; // reset
      // std::cout << i << " " << fm << " " << fp << " " << f << std::endl;
      // require gradient to be similar to finite difference
      REQUIRE(grad[i] == Catch::Approx((1.0 / (2.0 * h)) * (fp - fm)));
    }
  }
}

TEST_CASE("MEHP Force Relaxation2 runs", "[analysis][MEHPForceRelaxation]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
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
      CHECK(forceRelaxer2.getPressure() == Catch::Approx(0.39911682390778536));
      CHECK(forceRelaxer2.getForce() == Catch::Approx(552894.1903005145));
      CHECK(forceRelaxer2.getResidualNorm() ==
            Catch::Approx(2.1242043665818353e6));
      REQUIRE_NOTHROW(forceRelaxer2.runForceRelaxation());
      CHECK(forceRelaxer2.getNrOfSprings() == 8142);
      CHECK(forceRelaxer2.getNrOfIterations() > 1);

      // conversion factors
      double kb = 1.381e-23; // Boltzmann, J/K
      double T = 300.;       // Temperature, K
      double sigmaToNm = 0.616;
      double sigmaToM = sigmaToNm * 1.e-9;
      double slope = 0.00393 / (sigmaToNm * sigmaToNm); // sigma^2/(g/mol)
      double beadMass = 161.;                           // g/mol
      double Nb = 80.; // nr of beads per strand
      double conversionFactor =
        3. * kb * T / (slope * beadMass * Nb); // J/sigma^2
      CHECK(conversionFactor / (sigmaToM * sigmaToM) ==
            Catch::Approx(0.000245543));
      double nu =
        nrOfChains / (forceRelaxer2.getVolume() * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      CHECK(forceRelaxer2.getPressure() ==
            Catch::Approx(0.153806)); // LJ Units [?]
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
    forceRelaxer.runForceRelaxation("LD_MMA", 5);
    REQUIRE(forceRelaxer.getNrOfNodes() != universe.getNrOfAtoms());
    REQUIRE(forceRelaxer.getNrOfIterations() <= 5);
    REQUIRE(forceRelaxer.getNrOfIterations() >= 1);
    REQUIRE(universe.getAtomsOfType(2).size() == 7200);
    REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::MAX_STEPS);
    // run again, this time fully
    pcm::MEHPForceRelaxation forceRelaxer2 =
      pcm::MEHPForceRelaxation(universe, 2, true);
    forceRelaxer2.runForceRelaxation("LD_MMA", 10000, 1e-5, 1e-18);
    REQUIRE(forceRelaxer2.getNrOfIterations() > 5);
    CHECK(forceRelaxer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceRelaxer2.getGammaFactor(25, forceRelaxer2.getNrOfSprings()) ==
          Catch::Approx(1. / 3.).epsilon(0.001));
    auto stressTensor = forceRelaxer2.getStressTensor();
    CHECK(forceRelaxer2.getPressure() ==
          (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.);
    CHECK(forceRelaxer2.getResidualNorm() > 0.0);
    CHECK(forceRelaxer2.getForce() > 0.0);
    // TODO: find better, more accurate tests here
    CHECK(forceRelaxer2.getNrOfActiveNodes() > 1);
    CHECK(forceRelaxer2.getNrOfActiveSprings() > 1);
    CHECK(forceRelaxer2.getAverageSpringLength() > 1.0);

    pe::Universe universe3 = forceRelaxer2.getCrosslinkerVerse();
    CHECK(universe3.getNrOfAtoms() == forceRelaxer2.getNrOfNodes());
    CHECK(universe3.getNrOfBonds() == forceRelaxer2.getNrOfSprings());
    CHECK(universe3.getAtomsOfType(2).size() == universe3.getNrOfAtoms());

    // try out different algorithms
    std::vector<std::string> algorithms = { "LD_LBFGS",
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
}
