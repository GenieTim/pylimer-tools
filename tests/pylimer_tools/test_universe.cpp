#include <catch2/catch_test_macros.hpp>
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include <iostream>
extern "C" {
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;

TEST_CASE("Universe can be used", "[entity][Universe]")
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
      // works twice: make sure there is no removal of anything happening
      REQUIRE(universe.getMolecules(2).size() == 4);
      universe.addBonds(2, {{0, 3}}, {{3, 4}});
      REQUIRE(universe.getNrOfAtoms() == 4);
      REQUIRE(universe.getNrOfBonds() == 2);
      REQUIRE(universe.getMolecules(1).size() == 0);
      std::vector<pe::Molecule> molecules = universe.getMolecules(2);
      REQUIRE(molecules.size() == 2);
      // works twice: make sure there is no removal of anything happening
      REQUIRE(universe.getMolecules(2).size() == 2);
    }
  }

  SECTION("Molecules with crosslinkers are found")
  {
    /**
    # The system looks like this (in terms of bonds, not 3D placement):
    # 1-2-3-*6
    # |      |
    # *7-5---|
    # 8
    #
    # *4
    */
    universe.addAtoms(8,
                      {{1, 2, 3, 4, 5, 6, 7, 8}},    // id
                      {{1, 1, 1, 2, 1, 2, 2, 1}},    // type
                      {{1.25, 2, 3, 2, 1, 1, 1, 2}}, // x
                      {{1, 1, 1, 2, 3, 1, 1, 2}},    // y
                      {{1, 1, 1, 1, 1, 2, 3, 2}},    // z
                      {{1, 1, 1, 1, 1, 1, 1, 1}},    // nx
                      {{1, 1, 1, 1, 1, 1, 1, 1}},    // ny
                      {{1, 1, 1, 1, 1, 1, 1, 1}}     // nz
    );
    universe.addBonds(7, {{1, 3, 5, 1, 5, 3, 7}}, {{2, 2, 6, 7, 7, 6, 8}});
    // get atoms with type returns
    REQUIRE(universe.getAtomsWithType(2).size() == 3);
    REQUIRE(universe.getAtomsWithType(1).size() == 5);
    REQUIRE(universe.getAtomsWithType(0).size() == 0);
    // get atoms with type returns atoms with properties
    REQUIRE(universe.getMolecules(2)[0].getAtomsWithType(1).size() == 3);
    REQUIRE(universe.getMolecules(2)[1].getAtomsWithType(1)[0].getId() == 5);
    REQUIRE(universe.getMolecules(2)[2].getAtomsWithType(1)[0].getType() == 1);
    // get atoms with crosslinkers returns
    auto chains = universe.getChainsWithCrosslinker(2);
    REQUIRE(chains.size() == 3);
    REQUIRE(chains[0].getAtoms()[0].getId() == 1);
    REQUIRE(chains[0].getAtoms()[0].getX() == 1.25);
    REQUIRE(chains[0].getNrOfAtoms() == 5);
    REQUIRE(chains[0].getKey() == "1-2-3-6-7");
    REQUIRE(chains[0].getAtomsWithType(2).size() == 2);
    //
    auto functionalityPerType = universe.determineFunctionalityPerType();
    REQUIRE(functionalityPerType[1] == 2);
    REQUIRE(functionalityPerType[2] == 3);
    REQUIRE(functionalityPerType.size() == 2);
  }

  SECTION("Molecule Types are determined correctly")
  {
    /**
    # The system looks like this (in terms of bonds, not 3D placement):
    # 1-2-3
    #
    # *7-*6-5
    */
    universe.addAtoms(6, {{1, 2, 3, 5, 6, 7}}, // id
                      {{1, 1, 1, 1, 2, 2}},    // type
                      {{2, 2, 2, 2, 2, 2}},    // x
                      {{1, 1, 1, 1, 1, 1}},    // y
                      {{1, 1, 1, 1, 1, 1}},    // z
                      {{1, 1, 1, 1, 1, 1}},    // nx
                      {{1, 1, 1, 1, 1, 1}},    // ny
                      {{1, 1, 1, 1, 1, 1}}     // nz
    );
    universe.addBonds(4, {{1, 2, 5, 6}}, {{2, 3, 6, 7}});
    REQUIRE(universe.getNrOfBonds() == 4);
    REQUIRE(universe.getMolecules(2).size() == 2);
    REQUIRE(universe.getChainsWithCrosslinker(2).size() == 2);
    auto chains = universe.getChainsWithCrosslinker(2);
    REQUIRE(chains.size() == 2);
    REQUIRE(chains[0].getNrOfAtoms() == 3);
    REQUIRE(chains[1].getNrOfAtoms() == 2);
  }
}
