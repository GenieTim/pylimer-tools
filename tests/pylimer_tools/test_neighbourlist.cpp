#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/NeighbourList.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;

TEST_CASE("NeighbourList works as intended", "[entity][NeighbourList]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  universeSeq.initializeFromDataSequence(
    { { suspectedPath + "lammps_data_file_small.out" } });
  REQUIRE(universeSeq.getLength() == 1);
  REQUIRE(universeSeq.atIndex(0).getNrOfAtoms() == 12);
  REQUIRE(universeSeq.atIndex(0).getNrOfBonds() == 5);

  pe::Universe universe = universeSeq.atIndex(0);

  pe::NeighbourList neighbourList =
    pe::NeighbourList(universe.getAtoms(), universe.getBox(), 3.0);

  std::vector<pe::Atom> neighbours =
    neighbourList.getAtomsCloseTo(universe.getAtom(10000));
  REQUIRE(neighbours.size() > 0);

  neighbours = neighbourList.getAtomsCloseTo(universe.getAtom(10000), 0.0);
  REQUIRE(neighbours.size() == 0);
}
