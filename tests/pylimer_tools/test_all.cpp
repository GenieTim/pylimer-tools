#include <catch2/catch_test_macros.hpp>
#include "../../src/_pylimer_tools/entities/Universe.h"
#include "../../src/_pylimer_tools/entities/Molecule.h"
#include "../../src/_pylimer_tools/entities/Atom.h"
#include "../../src/_pylimer_tools/entities/Box.h"
#include <iostream>
#include <igraph/igraph.h>

namespace pe = pylimer_tools::entities;

TEST_CASE("Universe can be used", "[Universe]")
{
  REQUIRE(1 == 1);

  pe::Universe universe = pe::Universe(1.0, 1.0, 1.0);
  REQUIRE(universe.getVolume() == 1.0);

  SECTION("resizing smaller changes volume")
  {
    universe.setBox(pe::Box(0.0, 1.0, 2.0));
    REQUIRE(universe.getVolume() == 0.0);
  }

  SECTION("resizing bigger changes volume")
  {
    universe.setBoxLengths(1.0, 1.0, 1.0);
    REQUIRE(universe.getVolume() == 1.0);
    universe.setBox(pe::Box(2.0, 1.0, 2.0));
    REQUIRE(universe.getVolume() == 4.0);
  }

  SECTION("atoms can be added")
  {
    universe.addAtoms(2, {{0, 1}}, {{1, 1}}, {{0.0, 1.0}}, {{0.0, 1.0}}, {{0.0, 1.0}}, {{0, 0}}, {{0, 0}}, {{0, 0}});
    REQUIRE(universe.getNrOfAtoms() == 2);
    universe.addAtoms(2, {{3, 4}}, {{1, 1}}, {{0.0, 1.0}}, {{0.0, 1.0}}, {{0.0, 1.0}}, {{0, 0}}, {{0, 0}}, {{0, 0}});
    REQUIRE(universe.getNrOfAtoms() == 4);
    CHECK(universe.getAtomsWithType(1).size() == 4);
    CHECK(universe.getNrOfBonds() == 0);

    SECTION("molecules are found")
    {
      REQUIRE(universe.getNrOfAtoms() == 4);
      REQUIRE(universe.getMolecules(2).size() == 4);
      universe.addBonds(2, {{0, 3}}, {{3, 4}});
      REQUIRE(universe.getNrOfAtoms() == 4);
      REQUIRE(universe.getNrOfBonds() == 2);
      REQUIRE(universe.getMolecules(1).size() == 0);
      std::vector<pe::Molecule> molecules = universe.getMolecules(2);
      REQUIRE(molecules.size() == 2);
    }
  }
}
