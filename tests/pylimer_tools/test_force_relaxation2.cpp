#include "../../src/pylimer_tools_cpp/calc/MEHPForceRelaxation2.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <map>
#include <vector>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::calc::mehp;

TEST_CASE("MEHP Force Relaxation2 runs", "[analysis][MEHPForceRelaxation2]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_0.90015_pdms_1e4_a_78_bs_t_807536.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      pcm::MEHPForceRelaxation2 forceRelaxer2 =
        pcm::MEHPForceRelaxation2(universe2, 2);
      REQUIRE(forceRelaxer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE(forceRelaxer2.getNrOfIterations() == 0);
      REQUIRE(forceRelaxer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
      REQUIRE_NOTHROW(forceRelaxer2.runForceRelaxation());
      REQUIRE(forceRelaxer2.getNrOfIterations() > 1);
      REQUIRE(forceRelaxer2.getFinalPressure() == Catch::Approx(2.49685e8));
      double kb = 1.381e-23; // Boltzmann, J/K
      double T = 300.;       // Temperature, K
      double sigmaToNm = 0.62;
      double sigmaToM = sigmaToNm * 1.e-9;
      double slope = sigmaToNm * sigmaToNm * 0.00393; // sigma^2/(g/mol)
      double beadMass = 161.;                         // g/mol
      double Nb = 80;                                 // nr of beads per strand
      double conversionFactor = 3 * kb * T / (1e-18 * slope * beadMass * Nb);
      REQUIRE(conversionFactor ==
              Catch::Approx(0.000245543 / (sigmaToNm * sigmaToNm)));
      REQUIRE(forceRelaxer2.getGammaEq() * forceRelaxer2.getNrOfSprings() /
                (forceRelaxer2.getVolume() * sigmaToM * sigmaToM * sigmaToM) ==
              Catch::Approx(61308.3));
      REQUIRE(forceRelaxer2.getGammaEq() / (slope * Nb * beadMass) ==
              Catch::Approx(0.319446));
      REQUIRE(forceRelaxer2.getFinalPressure() / conversionFactor ==
              Catch::Approx(61308.3));
      REQUIRE(forceRelaxer2.getExitReason() == pcm::ExitReason::TOLERANCE);
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
    pcm::MEHPForceRelaxation2 forceRelaxer =
      pcm::MEHPForceRelaxation2(universe, 2);
    REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::UNSET);
    REQUIRE(forceRelaxer.getNrOfIterations() == 0);
    REQUIRE(forceRelaxer.getVolume() == Catch::Approx(universe.getVolume()));
    // setting the steps so low does not change much,
    // as the internal line-search steps are not counted
    // -> not a fast experiment anymore
    // forceRelaxer.runForceRelaxation(5, 1e-5, 15, true);
    // REQUIRE(forceRelaxer.getNrOfIterations() <= 5);
    // REQUIRE(forceRelaxer.getNrOfIterations() >= 1);
    // // REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::MAX_STEPS);
    // REQUIRE(forceRelaxer.getNrOfNodes() != universe.getNrOfAtoms());
    // REQUIRE(universe.getAtomsOfType(2).size() == 7200);
    // REQUIRE(forceRelaxer.getGammaEq() == Catch::Approx(1.
    // / 3.).epsilon(0.01));
  }
}
