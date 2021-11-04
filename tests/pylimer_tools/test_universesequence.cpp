#include <catch2/catch_test_macros.hpp>
#include "../../src/_pylimer_tools/entities/UniverseSequence.h"
#include "../../src/_pylimer_tools/entities/Universe.h"
#include "../../src/_pylimer_tools/entities/Molecule.h"
#include "../../src/_pylimer_tools/entities/Atom.h"
#include "../../src/_pylimer_tools/entities/Box.h"
#include <iostream>
#include <string>
#include <filesystem>
#include <igraph/igraph.h>

namespace pe = pylimer_tools::entities;

TEST_CASE("UniverseSequence can be used", "[entity][UniverseSequence]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  SECTION("Reading from dump file works")
  {
    universeSeq.initializeFromDumpFile(suspectedPath + "lammps_data_file_small.out", suspectedPath + "lammps_dump_small.lammpstrj");
    REQUIRE(universeSeq.getLength() == 1);
    REQUIRE(universeSeq.atIndex(0).getNrOfAtoms() == 12);
    REQUIRE(universeSeq.atIndex(0).getNrOfBonds() == 5);
  }

  SECTION("Reading from data files works")
  {
    universeSeq.initializeFromDataSequence({{suspectedPath + "lammps_data_file.out"}});
    REQUIRE(universeSeq.getLength() == 1);
    REQUIRE(universeSeq.atIndex(0).getNrOfAtoms() == 2594);
    REQUIRE(universeSeq.atIndex(0).getNrOfBonds() == 1050);
  }
}
