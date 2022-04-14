#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <string>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;

TEST_CASE("UniverseSequence can be used", "[entity][UniverseSequence]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  REQUIRE(std::filesystem::exists(suspectedPath));

  SECTION("Reading from dump file works")
  {
    universeSeq.initializeFromDumpFile(
      suspectedPath + "lammps_data_file_small.out",
      suspectedPath + "lammps_dump_small.lammpstrj");
    REQUIRE(universeSeq.getLength() == 1);
    REQUIRE_THROWS(universeSeq.atIndex(1));
    REQUIRE(universeSeq.atIndex(0).getNrOfAtoms() == 12);
    REQUIRE(universeSeq.atIndex(0).getNrOfBonds() == 5);
    REQUIRE(universeSeq.atIndex(0).getTimestep() == 70764);
    auto universes = universeSeq.getAll();
    REQUIRE(universes.size() == 1);
  }

  SECTION("Reading from data files works")
  {
    universeSeq.initializeFromDataSequence(
      { { suspectedPath + "lammps_data_file.out" } });
    REQUIRE(universeSeq.getLength() == 1);
    REQUIRE(universeSeq.atIndex(0).getNrOfAtoms() == 3000);
    REQUIRE(universeSeq.atIndex(0).getNrOfBonds() == 2900);
  }

  SECTION("Reading large files is sensibly fast")
  {
    universeSeq.initializeFromDumpFile(suspectedPath + "big_dump_file_data.out",
                                       suspectedPath +
                                         "big_dump_file.lammpstrj");
    // can be queried "randomly"
    pe::Universe universe = universeSeq.atIndex(10);
    REQUIRE(universe.getNrOfAtoms() == 32);
    REQUIRE(universeSeq.getLength() == 74322);
    // and can be queried again
    universeSeq.forgetAtIndex(10);
    pe::Universe universeAgain = universeSeq.atIndex(10);
    REQUIRE(universeAgain.getNrOfAtoms() == 32);
    REQUIRE(universeAgain.getAtom(1).getX() == universeAgain.getAtom(1).getX());
    // and the last one
    std::cout << "Requesting last index" << std::endl;
    REQUIRE_THROWS(universeSeq.atIndex(74322));
    std::cout << "Requesting last existing index" << std::endl;
    pe::Universe thirdUniverse = universeSeq.atIndex(74321);
    REQUIRE(thirdUniverse.getNrOfAtoms() == 32);
    REQUIRE(universeAgain.getAtom(1).getX() != thirdUniverse.getAtom(1).getX());
    // and back again
    std::cout << "Requesting fourth existing index" << std::endl;
    pe::Universe fourthUniverse = universeSeq.atIndex(9);
    REQUIRE(fourthUniverse.getNrOfAtoms() == 32);
    REQUIRE(fourthUniverse.getAtom(1).getX() !=
            thirdUniverse.getAtom(1).getX());
  }
}
