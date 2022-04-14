#include "../../src/pylimer_tools_cpp/calc/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <map>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::calc::mehp;

TEST_CASE("MEHP Force Relaxation runs", "[analysis][MEHPForceRelaxation]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));
  universeSeq.initializeFromDataSequence(
    { { suspectedPath + "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
                        "333_2d_t_7500001.structure.out" } });
  REQUIRE(universeSeq.getLength() == 1);
  pe::Universe universe = universeSeq.atIndex(0);
  pcm::MEHPForceRelaxation forceRelaxer = pcm::MEHPForceRelaxation(universe, 2);
  REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::UNSET);
  REQUIRE(forceRelaxer.getVolume() == Catch::Approx(universe.getVolume()));
  forceRelaxer.runForceRelaxation(2, 5, 1e-6, 25, true, 0.075);
  REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::MAX_STEPS);
  REQUIRE(forceRelaxer.getNrOfNodes() != universe.getNrOfAtoms());

  // forceRelaxer.runForceRelaxation(2, 500000, 1e-6, 25, true, 0.075);
  // REQUIRE(forceRelaxer.getGammaEq() == Catch::Approx(1. / 3.).epsilon(0.01));
  // REQUIRE(forceRelaxer.getVolume() == Catch::Approx(universe.getVolume()));
  // REQUIRE(forceRelaxer.getExitReason() == pcm::ExitReason::TOLERANCE);
}
