#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
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

TEST_CASE("Molecules work as intended", "[entity][Molecule]")
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
  std::vector<pe::Molecule> molecules = universe.getMolecules(0);

  SECTION("Special constructors work")
  {
    pe::Molecule molecule1 = molecules[0];
    pe::Molecule molecule2 = molecules[1];
    molecule2 = molecule1;
    REQUIRE(molecule2.getNrOfAtoms() == molecule1.getNrOfAtoms());
    REQUIRE(molecule2.getKey() == molecule1.getKey());
    REQUIRE(molecule2.getNrOfBonds() == molecule1.getNrOfBonds());
  }

  SECTION("Molecules can calculate")
  {
    REQUIRE(universe.getNrOfBonds() == 5);
    REQUIRE(universe.validate());
    REQUIRE(universe.getVolume() == Catch::Approx(64.3659 * 64.3659 * 64.3659));
    REQUIRE(molecules.size() == 7);

    // thankfully deterministic thanks to order of vertex ids
    std::vector<int> expectedLengths = { { 3, 2, 1, 2, 1, 2, 1 } };
    std::vector<double> expectedEndToEndDistances = {
      { 6.7651120738, 0.0324414569, 0.0, 0.0213056989, 0.0, 6.721720793, 0.0 }
    };

    int iteration = 0;
    for (pe::Molecule molecule : molecules) {
      REQUIRE(molecule.getLength() == expectedLengths[iteration]);
      REQUIRE(molecule.getBox()->getLx() ==
              Catch::Approx(32.182950030000001 * 2));
      REQUIRE(molecule.getBox()->getLy() ==
              Catch::Approx(32.182950030000001 * 2));
      REQUIRE(molecule.getBox()->getLz() ==
              Catch::Approx(32.182950030000001 * 2));
      REQUIRE(molecule.computeEndToEndDistance() ==
              Catch::Approx(expectedEndToEndDistances[iteration]));
      ++iteration;
    }

    pe::Molecule molecule1 = molecules[0];
    std::vector<double> bondLengths = molecule1.computeBondLengths();
    REQUIRE(bondLengths[0] == Catch::Approx(0.021513));
    REQUIRE(bondLengths[1] == Catch::Approx(6.74369));
    REQUIRE(molecule1.getNrOfAtoms() == 3);
    REQUIRE(molecule1.computeWeight() == Catch::Approx(3.0));
    REQUIRE(molecule1.getType() == pe::MoleculeType::UNDEFINED);
    REQUIRE(molecule1.computeRadiusOfGyration() == Catch::Approx(0.0));
    REQUIRE_THROWS(molecule1.getIdxByAtomId(19234121));

    std::vector<pe::Atom> atomsInLine = molecule1.getAtomsLinedUp();
    REQUIRE(atomsInLine.size() == molecule1.getNrOfAtoms());
    REQUIRE(molecule1.getAtomsOfDegree(2).size() == 1);
    REQUIRE(molecule1.getAtomsOfDegree(1).size() == 2);

    pe::Atom firstAtom = molecule1.getAtomByVertexIdx(0);
    REQUIRE(firstAtom.getId() == molecule1.getAtomIdByIdx(0));
    REQUIRE(0 == molecule1.getIdxByAtomId(firstAtom.getId()));
  }
}
