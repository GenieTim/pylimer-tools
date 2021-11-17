#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
extern "C" {
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;

TEST_CASE("Molecules work as intended", "[entity][Molecule]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  universeSeq.initializeFromDataSequence({{suspectedPath + "lammps_data_file_small.out"}});
  REQUIRE(universeSeq.getLength() == 1);
  REQUIRE(universeSeq.atIndex(0).getNrOfAtoms() == 12);
  REQUIRE(universeSeq.atIndex(0).getNrOfBonds() == 5);

  pe::Universe universe = universeSeq.atIndex(0);

  SECTION("Molecules can calculate")
  {
    std::vector<pe::Molecule> molecules = universe.getMolecules(0);
    REQUIRE(universe.getNrOfBonds() == 5);
    REQUIRE(universe.validate());
    REQUIRE(molecules.size() == 7);

    // thankfully deterministic thanks to order of vertex ids
    std::vector<int> expectedLengths = {{3, 2, 1, 2, 1, 2, 1}};
    std::vector<double> expectedEndToEndDistances = {{6.7651120738, 0.0324414569, 0.0, 0.0213056989, 0.0, 6.721720793, 0.0}};

    int iteration = 0;
    for (pe::Molecule molecule : molecules)
    {
      REQUIRE(molecule.getLength() == expectedLengths[iteration]);
      REQUIRE(molecule.computeEndToEndDistance() == Catch::Approx(expectedEndToEndDistances[iteration]));
      ++iteration;
    }

    pe::Molecule molecule1 = molecules[0];
    std::vector<double> bondLengths = molecule1.computeBondLengths();
    REQUIRE(bondLengths[0] == Catch::Approx(0.021513));
    REQUIRE(bondLengths[1] == Catch::Approx(6.74369));
  }
}
