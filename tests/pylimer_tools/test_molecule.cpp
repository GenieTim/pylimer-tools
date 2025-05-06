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
  std::cout << "Running test \"Molecules work as intended\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  REQUIRE(std::filesystem::exists(suspectedPath));

  universeSeq.initializeFromDataSequence(
    { { suspectedPath + "/lammps_data_file_small.out" } });
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
    std::vector<std::string> expectedKeys = { { "10000-20000-90000",
                                                "30000-40000",
                                                "50360",
                                                "50000-60000",
                                                "50908",
                                                "70000-80000",
                                                "74363" } };
    std::vector<double> expectedEndToEndDistances = {
      { 6.7651120738, 0.0324414569, 0.0, 0.0213056989, 0.0, 6.721720793, 0.0 }
    };

    CHECK(molecules.size() == expectedLengths.size());
    CHECK(molecules[0].getAtomsOfDegree(2).size() == 1);
    CHECK(molecules[0].getAtomsOfDegree(1).size() == 2);

    int iteration = 0;
    for (pe::Molecule molecule : molecules) {
      CHECK(molecule.getKey() == expectedKeys[iteration]);
      REQUIRE(molecule.getLength() == expectedLengths[iteration]);
      REQUIRE(molecule.getBox().getLx() ==
              Catch::Approx(32.182950030000001 * 2));
      REQUIRE(molecule.getBox().getLy() ==
              Catch::Approx(32.182950030000001 * 2));
      REQUIRE(molecule.getBox().getLz() ==
              Catch::Approx(32.182950030000001 * 2));
      CHECK(molecule.computeEndToEndDistance() ==
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
    // REQUIRE(molecule1.computeRadiusOfGyration() == Catch::Approx(0.0));
    REQUIRE_THROWS(molecule1.getIdxByAtomId(19234121));

    std::vector<pe::Atom> atomsInLine = molecule1.getAtomsLinedUp();
    REQUIRE(atomsInLine.size() == molecule1.getNrOfAtoms());
    REQUIRE(molecule1.getAtomsOfDegree(2).size() == 1);
    REQUIRE(molecule1.getAtomsOfDegree(1).size() == 2);
    REQUIRE(atomsInLine[0].getId() == 20000);
    REQUIRE(atomsInLine[2].getId() == 90000);

    Eigen::Vector3d boxDistance = molecule1.getOverallBondSum(2);
    // TODO: check
    CHECK(boxDistance[0] == Catch::Approx(2.106517));
    CHECK(boxDistance[1] == Catch::Approx(5.01925855));
    CHECK(boxDistance[2] == Catch::Approx(4.01701022));

    pe::Atom firstAtom = molecule1.getAtomByVertexIdx(0);
    REQUIRE(firstAtom.getId() == molecule1.getAtomIdByIdx(0));
    REQUIRE(0 == molecule1.getIdxByAtomId(firstAtom.getId()));
  }
}

TEST_CASE("Molecules sum the bonds correctly", "[entity][Molecule]")
{
  std::cout << "Running test \"Molecules sum the bonds correctly\""
            << std::endl;
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);

  /**
   * This Universe contains one single chain of 5 atoms,
   *   1 - 2 - 3 - 4 - 5
   */
  universe.addAtoms({ { 1, 2, 3, 4, 5 } },        // id
                    { { 2, 1, 1, 1, 2 } },        // type
                    { { 5., 7., 9., 11., 13. } }, // x
                    { { 0., 0., 0., 0., 0 } },    // y
                    { { 0., 0., 0., 0., 0 } },    // z
                    { { 0, 0, 0, 0, 0 } },
                    { { 0, 0, 0, 0, 0 } },
                    { { 0, 0, 0, 0, 0 } });
  universe.addBonds({ { 1, 2, 3, 4 } }, { { 2, 3, 4, 5 } });

  std::vector<pe::Molecule> molecules = universe.getMolecules();
  REQUIRE(molecules.size() == 1);

  Eigen::Vector3d overallSum = molecules[0].getOverallBondSum(2);
  CHECK(overallSum[0] == 13. - 5.);
  CHECK(overallSum[1] == 0.);
  CHECK(overallSum[2] == 0.);

  Eigen::Vector3d partialSum = molecules[0].getOverallBondSumFromTo(1, 2);
  CHECK(partialSum[0] == 2.);
  CHECK(partialSum[1] == 0.);
  CHECK(partialSum[2] == 0.);
  CHECK(molecules[0].getNrOfBondsFromTo(1, 2) == 1);

  // reverse order of indices
  CHECK_THROWS(molecules[0].getOverallBondSumFromTo(2, 1));
  CHECK_THROWS(molecules[0].getNrOfBondsFromTo(2, 1));
  partialSum = molecules[0].getOverallBondSumFromTo(2, 1, 2, false);
  CHECK(partialSum[0] == -2.);
  CHECK(partialSum[1] == 0.);
  CHECK(partialSum[2] == 0.);

  CHECK(molecules[0].getNrOfBondsFromTo(1, 2, 2, false) == 1);
  CHECK(molecules[0].getNrOfBondsFromTo(2, 4, 2, true) == 2);
  CHECK(molecules[0].getNrOfBondsFromTo(2, 4, 2, false) == 2);
  CHECK(molecules[0].getNrOfBondsFromTo(4, 2, 2, false) == 2);
  CHECK(molecules[0].getNrOfBondsFromTo(1, 5) == molecules[0].getNrOfBonds());

  CHECK(molecules[0].getPathLength(0, 1) == 1);
  CHECK(molecules[0].getPathLength(1, 3, 1) == 0);
  CHECK(molecules[0].getPathLength(1, 3) == 2);
  CHECK(molecules[0].getPathLength(3, 1) == 2);
  CHECK(molecules[0].getPathLength(0, 4) == molecules[0].getNrOfBonds());
}

TEST_CASE("Molecules compute radius of gyration", "[entity][Molecule]")
{
  std::cout << "Running test \"Molecules compute radius of gyration\""
            << std::endl;

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
  const pe::Box box = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
  universe.setBox(box);
  universe.addAtoms({ { 1, 2, 3, 4, 5, 6, 7, 8 } },     // id
                    { { 2, 1, 1, 1, 2, 1, 1, 1 } },     // type
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // x
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // y
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // z
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } },     // nx
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } },     // ny
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } }      // nz
  );

  std::map<int, double> masses;
  masses[1] = 1.0;
  masses[2] = 1.0;
  universe.addBonds(8,
                    { { 1, 2, 3, 4, 5, 6, 7, 8 } },
                    { { 2, 3, 4, 1, 6, 7, 8, 5 } },
                    { { 1, 1, 1, 1, 1, 1, 1, 1 } },
                    false,
                    false);

  std::vector<pe::Molecule> molecules = universe.getMolecules();
  SECTION("Coordinates are computed appropriately")
  {
    igraph_vector_int_t all_vertices;
    igraph_vector_int_init(&all_vertices, universe.getNrOfAtoms());
    for (int i = 1; i <= universe.getNrOfAtoms(); ++i) {
      igraph_vector_int_set(&all_vertices, i - 1, universe.getIdxByAtomId(i));
    }
    CHECK(igraph_vector_int_size(&all_vertices) == universe.getNrOfAtoms());
    Eigen::VectorXd coordinates =
      universe.getUnwrappedVertexCoordinates(all_vertices, box);
    CHECK(coordinates.size() == universe.getNrOfAtoms() * 3);
    igraph_vector_int_destroy(&all_vertices);

    for (size_t i = 0; i < 4; ++i) {
      pe::Atom atomI = universe.getAtom(5 + i);
      CHECK(atomI.getUnwrappedX(box) ==
            Catch::Approx(29. + static_cast<double>(i)));
      CHECK(atomI.getUnwrappedY(box) ==
            Catch::Approx(29. + static_cast<double>(i)));
      CHECK(atomI.getUnwrappedZ(box) ==
            Catch::Approx(29. + static_cast<double>(i)));
      CHECK(coordinates[(4 + i) * 3] ==
            Catch::Approx(29. + static_cast<double>(i)));
      CHECK(coordinates[(4 + i) * 3 + 1] ==
            Catch::Approx(29. + static_cast<double>(i)));
      CHECK(coordinates[(4 + i) * 3 + 2] ==
            Catch::Approx(29. + static_cast<double>(i)));
    }
  }
  CHECK_THROWS(molecules[0].computeRadiusOfGyration());
  CHECK_THROWS(molecules[0].computeRadiusOfGyrationWithDerivedImageFlags());
  universe.setMasses(masses);
  molecules = universe.getMolecules();
  CHECK(molecules.size() == 2);

  SECTION("Coordinates are derived correctly")
  {
    std::vector<igraph_integer_t> vertices;
    for (int i = 0; i < molecules[1].getNrOfAtoms(); i++) {
      vertices.push_back(i);
    }
    Eigen::VectorXd assumedCoordinates =
      Eigen::VectorXd::Zero(molecules[1].getNrOfAtoms() * 3);
    assumedCoordinates =
      molecules[1].getAssumedVertexCoordinates<Eigen::VectorXd>(
        assumedCoordinates, box, vertices);
    for (int i = 0; i < 4; ++i) {
      CHECK(assumedCoordinates[(i) * 3] ==
            Catch::Approx(9. + static_cast<double>(i)));
      CHECK(assumedCoordinates[(i) * 3 + 1] ==
            Catch::Approx(9. + static_cast<double>(i)));
      CHECK(assumedCoordinates[(i) * 3 + 2] ==
            Catch::Approx(9. + static_cast<double>(i)));
    }
  }

  SECTION("R_g can be computed in both ways")
  {
    // first method
    CHECK(molecules[0].computeRadiusOfGyration() == Catch::Approx(15. / 4.));
    // second method
    CHECK(molecules[0].computeRadiusOfGyrationWithDerivedImageFlags() ==
          Catch::Approx(15. / 4.));

    // first method
    CHECK(molecules[1].computeRadiusOfGyration() == Catch::Approx(15. / 4.));
    // second method
    CHECK(molecules[1].computeRadiusOfGyrationWithDerivedImageFlags() ==
          Catch::Approx(15. / 4.));
  }
}

TEST_CASE("Strand Ends are found", "[entity][Molecule]")
{
  std::cout << "Running test \"Strand Ends are found\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath =
    "../pylimer_tools/fixtures/structure/crosslinked_M10000_N39_p_0.9.out";
  if (std::filesystem::exists(suspectedPath)) {
    universeSeq.initializeFromDataSequence({ { suspectedPath } });
    REQUIRE(universeSeq.getLength() == 1);

    pe::Universe universe = universeSeq.atIndex(0);
    std::vector<pe::Molecule> molecules = universe.getChainsWithCrosslinker(2);

    for (pe::Molecule chain : molecules) {
      std::vector<pe::Atom> chainEnds = chain.getChainEnds(2, true);
      CHECK(chainEnds.size() == 2);
      if (chain.getType() == pe::MoleculeType::PRIMARY_LOOP) {
        CHECK(chainEnds[0] == chainEnds[1]);
        CHECK(chain.getChainEnds(2, false).size() == 1);
      }
    }
  } else {
    std::cout << "Skipping test: file " << suspectedPath << " not found."
              << std::endl;
  }
}

TEST_CASE("Molecule equality works", "[entity][Molecule]")
{
  std::cout << "Running test \"Molecule equality works\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string suspectedPath = "../pylimer_tools/fixtures/structure/";

  std::string inputFile =
    suspectedPath + "square_lattice_2x2_a_5.2d.structure.out";
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);

  std::vector<pe::Molecule> chains = universe.getChainsWithCrosslinker(2);
  pe::Molecule chain0 = chains[0];

  CHECK(chain0 == chains[0]);
  CHECK(chain0 != chains[1]);

  CHECK(chain0.containsAtom(chain0.getAtomByVertexIdx(1)));
}
