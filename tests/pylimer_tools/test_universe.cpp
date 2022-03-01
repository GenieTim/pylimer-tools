#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <map>
#include <vector>
extern "C" {
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;

void outputLoops(std::map<int, std::vector<std::vector<pe::Atom>>>  loops) {
  for (const auto &myPair : loops) {
    std::cout << myPair.first << std::endl;
    for (const std::vector<pe::Atom> &as : myPair.second) {
      for (const pe::Atom &a : as) {
        std::cout << "\t" << a.getId() << std::endl;
      }
    }
  }
}

TEST_CASE("Universe can be used", "[entity][Universe]") {
  REQUIRE(1 == 1);

  pe::Universe universe = pe::Universe(1.0, 1.0, 1.0);
  REQUIRE(universe.getVolume() == 1.0);

  SECTION("resizing smaller changes volume") {
    universe.setBox(pe::Box(0.0, 1.0, 2.0));
    REQUIRE(universe.getVolume() == 0.0);
  }

  SECTION("resizing bigger changes volume") {
    universe.setBoxLengths(1.0, 1.0, 1.0);
    REQUIRE(universe.getVolume() == 1.0);
    universe.setBox(pe::Box(2.0, 1.0, 2.0));
    REQUIRE(universe.getVolume() == 4.0);
  }

  SECTION("atoms can be added") {
    universe.addAtoms(2, {{0, 1}}, {{1, 1}}, {{0.0, 1.0}}, {{0.0, 1.0}},
                      {{0.0, 1.0}}, {{0, 0}}, {{0, 0}}, {{0, 0}});
    REQUIRE(universe.getNrOfAtoms() == 2);
    universe.addAtoms(2, {{3, 4}}, {{1, 1}}, {{0.0, 1.0}}, {{0.0, 1.0}},
                      {{0.0, 1.0}}, {{0, 0}}, {{0, 0}}, {{0, 0}});
    REQUIRE(universe.getNrOfAtoms() == 4);
    CHECK(universe.getAtomsWithType(1).size() == 4);
    CHECK(universe.getNrOfBonds() == 0);

    SECTION("molecules are found") {
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

  SECTION("Molecules with crosslinkers are found") {
    /**
    # The system looks like this (in terms of bonds, not 3D placement):
    # 1-2-3-*6
    # |      |
    # *7-5---|
    # 8
    #
    # *4
    */
    universe.addAtoms(8, {{1, 2, 3, 4, 5, 6, 7, 8}},    // id
                      {{1, 1, 1, 2, 1, 2, 2, 1}},       // type
                      {{1.25, 2, 3, 1.01, 2, 4, 1, 1}}, // x
                      {{1, 1, 1, 4.01, 2, 1, 2, 3}},    // y
                      {{1, 1, 1, 1.01, 1, 1, 1, 1}},    // z
                      {{1, 1, 1, 1, 1, 1, 1, 1}},       // nx
                      {{1, 1, 1, 1, 1, 1, 1, 1}},       // ny
                      {{1, 1, 1, 1, 1, 1, 1, 1}}        // nz
    );
    universe.addBonds(7, {{1, 2, 3, 6, 5, 7, 7}}, {{2, 3, 6, 5, 7, 1, 8}},
                      {{1, 1, 1, 1, 1, 1, 11}});
    SECTION("get bonds returns") {
      auto edges = universe.getEdges();
      REQUIRE(edges["bond_from"][0] == 1 - 1);
      REQUIRE(edges["bond_to"][0] == 2 - 1);
      auto bonds = universe.getBonds();
      REQUIRE(bonds.at("bond_from")[0] == 1);
      REQUIRE(bonds.at("bond_to")[0] == 2);
      REQUIRE(bonds["bond_from"].size() == bonds["bond_to"].size());
      REQUIRE(bonds["bond_from"].size() == edges["bond_to"].size());
      REQUIRE(bonds["bond_from"].size() == edges["bond_from"].size());
      REQUIRE(bonds["bond_from"].size() == 7);
      REQUIRE(bonds["bond_from"][3] ==
              universe.getAtomIdByIdx(edges["bond_from"][3]));
      REQUIRE(bonds["bond_to"][3] ==
              universe.getAtomIdByIdx(edges["bond_to"][3]));
      REQUIRE(universe.getIdxByAtomId(bonds["bond_from"][2]) ==
              edges["bond_from"][2]);
      REQUIRE(universe.getIdxByAtomId(bonds["bond_to"][2]) ==
              edges["bond_to"][2]);
      // REQUIRE(bonds["bond_type"][5] == 1);
      // REQUIRE(bonds["bond_type"][6] == 11);
      // get atoms with type returns
      REQUIRE(universe.getAtomsWithType(2).size() == 3);
      REQUIRE(universe.getAtomsWithType(1).size() == 5);
      REQUIRE(universe.getAtomsWithType(0).size() == 0);
    }
    SECTION("get atoms returns") {
      // get atoms with type returns atoms with properties
      std::vector<pe::Molecule> molecules = universe.getMolecules(2);
      REQUIRE(molecules[0].getAtomsWithType(1).size() == 3);
      REQUIRE(molecules[1].getAtomsWithType(1)[0].getId() == 5);
      REQUIRE(molecules[2].getAtomsWithType(1)[0].getType() == 1);
      // get molecules allows to fetch atoms with degree
      REQUIRE(molecules[0].getAtomsOfDegree(2).size() == 1);
      REQUIRE(molecules[0].getAtomsOfDegree(1).size() == 2);
      REQUIRE(molecules[0].getAtomsOfDegree(0).size() == 0);
    }
    SECTION("get atoms with crosslinkers returns") {
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
      //
      REQUIRE(chains[0].getAtomsOfDegree(2).size() == 3);
      REQUIRE(chains[0].getAtomsOfDegree(1).size() == 2);
      REQUIRE(chains[0].getAtomsOfDegree(0).size() == 0);
    }
    SECTION("Loops are found") {
      std::map<int, std::vector<std::vector<pe::Atom>>> loops =
          universe.findLoops(2, -1, true);
      REQUIRE(loops.contains(2));
      REQUIRE(loops.size() == 1);
    }
    SECTION("Loops are found in reduced universe") {
      universe.addAtoms(2, {{9, 10}}, {{1, 1}}, {{0.0, 0.0}}, {{1.0, 0.0}},
                        {{1.0, 0.0}}, {{0, 1}}, {{0, 1}}, {{0, 1}});
      universe.addBonds(2, {{4, 9}}, {{9, 4}}, {{1, 1}}, false, false);
      pe::Universe reducedUniverse = universe.getNetworkOfCrosslinker(2);
      std::map<int, std::vector<std::vector<pe::Atom>>> loops =
          reducedUniverse.findLoops(2, -1);
      REQUIRE(loops.size() == 2);
    }
    SECTION("Infinite Strands are found") {
      universe.setBoxLengths(4.0, 4.0, 2.0);
      REQUIRE(universe.hasInfiniteStrand(2, -1) == false);
      universe.addBonds(2, {{8, 4}}, {{4, 1}});
      std::map<int, std::vector<std::vector<pe::Atom>>> loops =
          universe.findLoops(2, -1);
      REQUIRE(loops.size() == 3);
      REQUIRE(universe.hasInfiniteStrand(2, -1) == true);
    }
    SECTION("Reduction to Cross-linker-verse works") {
      pe::Universe reducedUniverse = universe.getNetworkOfCrosslinker(2);
      REQUIRE(reducedUniverse.getNrOfAtoms() == 3);
      REQUIRE(reducedUniverse.getAtomsWithType(2).size() == 3);

      auto bonds = reducedUniverse.getEdges();

      REQUIRE(reducedUniverse.getNrOfBonds() == 2);
      REQUIRE(bonds["bond_from"][0] == 1);
      REQUIRE(bonds["bond_to"][0] == 2);
      REQUIRE(bonds["bond_from"][1] == 1);
      REQUIRE(bonds["bond_to"][1] == 2);

      auto atom = reducedUniverse.getAtom(4);
      REQUIRE(atom.getX() == 1.01);
      REQUIRE(atom.getY() == 4.01);
      REQUIRE(atom.getZ() == 1.01);
    }
  }

  SECTION("Molecule Types are determined correctly") {
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
