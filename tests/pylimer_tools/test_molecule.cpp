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
    REQUIRE(molecule1.computeTotalMass() == Catch::Approx(3.0));
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

TEST_CASE("Molecules compute radius of gyration", "[entity][Molecule]")
{

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  /**
   * The system looks like this (in terms of bonds, not 3D placement):
   *
   * 1-2
   * | |
   * 4-3
   *
   * 5-6
   * | |
   * 8-7
   */
  universe.setBox(pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0));
  universe.addAtoms(8,
                    { { 1, 2, 3, 4, 5, 6, 7, 8 } },   // id
                    { { 2, 1, 1, 1, 2, 1, 1, 1 } },   // type
                    { { 1, 2, 3, 4, 9, -10, 9, 8 } }, // x
                    { { 1, 2, 3, 4, 9, -10, 9, 8 } }, // y
                    { { 1, 2, 3, 4, 9, -10, 9, 8 } }, // z
                    { { 1, 1, 1, 1, 0, 1, 1, 1 } },   // nx
                    { { 1, 1, 1, 1, 0, 1, 1, 1 } },   // ny
                    { { 1, 1, 1, 1, 0, 1, 1, 1 } }    // nz
  );
  universe.addBonds(8,
                    { { 1, 2, 3, 4, 5, 6, 7, 8 } },
                    { { 2, 3, 4, 1, 6, 7, 8, 5 } },
                    { { 1, 1, 1, 1, 1, 1, 1, 1 } },
                    false,
                    false);

  std::vector<pe::Molecule> molecules = universe.getMolecules();
  CHECK(molecules.size() == 2);

  SECTION("R_g can be computed in both ways")
  {
    // first method
    CHECK(molecules[0].computeRadiusOfGyration() == Catch::Approx(15. / 4.));
    // second method
    CHECK(molecules[0].computeRadiusOfGyrationWithDerivedImageFlags(2) ==
          Catch::Approx(15. / 4.));

    // first method
    CHECK(molecules[1].computeRadiusOfGyration() == Catch::Approx(15. / 4.));
    // second method
    CHECK(molecules[1].computeRadiusOfGyrationWithDerivedImageFlags(2) ==
          Catch::Approx(15. / 4.));
  }
}
