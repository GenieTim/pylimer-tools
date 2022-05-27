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

  SECTION("Empty dump files throw")
  {
    universeSeq.initializeFromDumpFile(suspectedPath +
                                         "lammps_data_file_small.out",
                                       suspectedPath + "empty_file.txt");
    REQUIRE_THROWS(universeSeq.next());
  }

  SECTION("Broken box is detected")
  {
    universeSeq.initializeFromDumpFile(
      suspectedPath + "lammps_data_file_small.out",
      suspectedPath + "lammps_dump_small_broken_box.lammpstrj");
    REQUIRE_THROWS(universeSeq.next());
  }

  SECTION("Missing box is augmented")
  {
    universeSeq.initializeFromDumpFile(
      suspectedPath + "lammps_data_file_small.out",
      suspectedPath + "lammps_dump_small_no_box.lammpstrj");
    pe::Universe universe = universeSeq.atIndex(0);
    REQUIRE(universe.getBox().getLx() ==
            3.2182950030000001e+01 + 3.2182950030000001e+01);
  }

  SECTION("Unwrapped atoms are read too")
  {
    universeSeq.initializeFromDumpFile(
      suspectedPath + "lammps_data_file_small.out",
      suspectedPath + "lammps_dump_small_unwrapped.lammpstrj");
    pe::Universe universe = universeSeq.atIndex(0);
    REQUIRE(universe.getNrOfAtoms() == 12);
  }

  SECTION("Angles are read, too")
  {
    universeSeq.initializeFromDataSequence(
      { { suspectedPath + "lammps_data_file_small_wangles.out" } });
    REQUIRE(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.next();
    REQUIRE(universe.getNrOfAngles() == 1);
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
    REQUIRE_THROWS(universeSeq.atIndex(74322));
    pe::Universe thirdUniverse = universeSeq.atIndex(74321);
    REQUIRE(thirdUniverse.getNrOfAtoms() == 32);
    REQUIRE(universeAgain.getAtom(1).getX() != thirdUniverse.getAtom(1).getX());
    // and back again
    pe::Universe fourthUniverse = universeSeq.atIndex(9);
    REQUIRE(fourthUniverse.getNrOfAtoms() == 32);
    REQUIRE(fourthUniverse.getAtom(1).getX() !=
            thirdUniverse.getAtom(1).getX());
  }

  // SECTION("Some non-commited large file reading works")
  // {
  //   universeSeq.initializeFromDumpFile(
  //     suspectedPath +
  //     "uncrosslinked_MD_melt_M10000_N79_equil_50M_bs_w_extra_"
  //                     "chain_near_min.out",
  //     suspectedPath + "melt_fene_N_78_rev.lammpstrj");
  //   REQUIRE(universeSeq.getLength() > 8638005/10);
  // }
}
