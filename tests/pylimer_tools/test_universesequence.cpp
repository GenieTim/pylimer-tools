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
extern "C"
{
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;

TEST_CASE("UniverseSequence can be used", "[entity][UniverseSequence]")
{
  std::cout << "Running test \"UniverseSequence can be used\"" << std::endl;
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
    REQUIRE(universeSeq.atIndex(0).getMasses()[1] == 1.0);
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

  SECTION("MSD computation works")
  {
    universeSeq.initializeFromDumpFile(
      suspectedPath + "lammps_data_file_small.out",
      suspectedPath + "lammps_dump_small_3step.lammpstrj");
    REQUIRE(universeSeq.getLength() == 3);
    std::vector<long int> atomIds = { 10000, 20000, 30000 };
    std::unordered_map<long int, double> msdForAtoms =
      universeSeq.computeMsdForAtoms(atomIds, 1, true);
    REQUIRE(msdForAtoms[1] == Catch::Approx(0.2364648387));
    // std::cout << "MSD computation works." << std::endl;
  }

  SECTION("Ree computation works")
  {
    universeSeq.initializeFromDumpFile(
      suspectedPath + "lammps_data_file_small.out",
      suspectedPath + "lammps_dump_small_3step.lammpstrj");
    REQUIRE(universeSeq.getLength() == 3);
    std::vector<long int> atomIdsFrom = { 10000, 20000, 30000 };
    std::vector<long int> atomIdsTo = { 10000, 20000, 30000 };
    std::unordered_map<long int, double> ree =
      universeSeq.computeDistanceAutocorrelationFromToAtoms(
        atomIdsFrom, atomIdsTo, 1, true);
    for (size_t dt = 0; dt < 2; ++dt) {
      CHECK(ree[dt + 1] == 0.0);
    }
    atomIdsFrom = { 70000 };
    atomIdsTo = { 80000 };
    ree = universeSeq.computeDistanceAutocorrelationFromToAtoms(
      atomIdsFrom, atomIdsTo, 1, true);
    pe::Universe universe = universeSeq.atIndex(0);
    pe::Box box = universe.getBox();
    Eigen::Vector3d coords1 =
      universe.getAtom(70000).getUnwrappedCoordinates(box);
    Eigen::Vector3d coords2 =
      universe.getAtom(80000).getUnwrappedCoordinates(box);
    for (size_t dt = 0; dt < 2; ++dt) {
      CHECK(ree[dt + 1] == Catch::Approx((coords2 - coords1).squaredNorm()));
    }
    // std::cout << "R_ee computation works." << std::endl;
    atomIdsFrom = { 70000, 80000, 74363, 70000 };
    atomIdsTo = { 80000, 74363, 70000, 90000 };
    for (size_t universe_i = 0; universe_i < universeSeq.getLength();
         universe_i++) {
      // make sure the atoms chosen above have the same coordinates in all
      // universes
      for (size_t atomI = 0; atomI < atomIdsFrom.size(); atomI++) {
        Eigen::Vector3d coordsFrom = universeSeq.atIndex(universe_i)
                                       .getAtom(atomIdsFrom[atomI])
                                       .getUnwrappedCoordinates(box);
        Eigen::Vector3d coordsTo = universeSeq.atIndex(universe_i)
                                     .getAtom(atomIdsTo[atomI])
                                     .getUnwrappedCoordinates(box);
        CHECK(coordsFrom.isApprox(
          universe.getAtom(atomIdsFrom[atomI]).getUnwrappedCoordinates(box)));
        CHECK(coordsTo.isApprox(
          universe.getAtom(atomIdsTo[atomI]).getUnwrappedCoordinates(box)));
      }
    }
    ree = universeSeq.computeDistanceAutocorrelationFromToAtoms(
      atomIdsFrom, atomIdsTo, 1, true);
    for (size_t dt = 0; dt < 2; ++dt) {
      double distMean = 0.0;
      for (size_t i = 0; i < atomIdsTo.size(); ++i) {
        coords1 = universe.getAtom(atomIdsFrom[i]).getUnwrappedCoordinates(box);
        coords2 = universe.getAtom(atomIdsTo[i]).getUnwrappedCoordinates(box);
        distMean += (coords2 - coords1).squaredNorm() /
                    static_cast<double>(atomIdsTo.size());
      }
      CHECK(distMean == Catch::Approx(ree[dt + 1]));
    }

    std::cout << "R_ee computation works." << std::endl;
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

  // SECTION("Some non-commited large file reading works")
  // {
  //   std::string dataFilePath =
  //     "/Volumes/drHobbies43/doctorate-hobbies-4/simulations/"
  //     "atomistic-simulations/structure/"
  //     "atomistic_melt_pdms_1000_a_43_v_1.structure.out";
  //   std::vector<std::string> dataFileSeq = { dataFilePath };
  //   std::vector<pu::AtomStyle> atomStyle = { pu::AtomStyle::FULL };
  //   universeSeq.setDataFileAtomStyle(atomStyle);
  //   universeSeq.initializeFromDataSequence(dataFileSeq);

  //   REQUIRE_NOTHROW(universeSeq.atIndex(0));
  // }
};
