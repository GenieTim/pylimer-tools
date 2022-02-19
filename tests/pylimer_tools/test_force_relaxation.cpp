#include "../../src/pylimer_tools_cpp/calc/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <map>
#include <vector>
extern "C" {
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::calc::mehp;

TEST_CASE("MEHP Force Relaxation runs", "[analysis][MEHPForceRelaxation]") {
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));
  universeSeq.initializeFromDataSequence(
      {{suspectedPath + "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0.666_2d_t_5000001.structure.out"}});
  REQUIRE(universeSeq.getLength() == 1);
  pe::Universe universe = universeSeq.atIndex(0);
  pcm::MEHPForceRelaxation forceRelaxer = pcm::MEHPForceRelaxation(universe, 2);
  forceRelaxer.runForceRelaxation(2, 250000, 1e-6);
  REQUIRE(forceRelaxer.getVolume() == Catch::Approx(universe.getVolume()));
  REQUIRE(forceRelaxer.getGammaEq() == Catch::Approx(1./3.));
}
