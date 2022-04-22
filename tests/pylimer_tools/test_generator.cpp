#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/utils/DataFileWriter.h"
#include "../../src/pylimer_tools_cpp/utils/MCUniverseGenerator.h"
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <map>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;

TEST_CASE("Certain configurations do not lead to memory corruption")
{
  // the following parameters have led to a `double free or corruption` error?!?
  int nrOfCrosslinkers = static_cast<int>(5e4 * 2 * 0.7 / 7);
  double sideLength = std::cbrt((10 * 5e4 * nrOfCrosslinkers) / 0.85);
  pu::MCUniverseGenerator generator =
    pu::MCUniverseGenerator(sideLength, sideLength, sideLength);
  REQUIRE_NOTHROW(generator.setSeed(68419));
  REQUIRE_NOTHROW(generator.setBeadDistance(0.965));
}

TEST_CASE("Universe can be generated", "[generator][MCUniverseGenerator]")
{
  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(100, 2);
  REQUIRE_THROWS(generator.addCrosslinkers(100, 3));
  generator.addSolventChains(100, 16, 3);
  generator.addAndLinkStrands((4 / 2) * 100, 16, 0.8);

  pe::Universe universe = generator.getUniverse();
  std::map<int, double> weights;
  weights[1] = 1.0;
  weights[2] = 1.0;
  weights[3] = 1.0;
  universe.setMasses(weights);
  REQUIRE(universe.getVolume() == 10.0 * 10.0 * 10.0);

  auto angles = universe.detectAngles();
  universe.addAngles(
    angles["angle_from"], angles["angle_via"], angles["angle_to"]);
  REQUIRE(universe.getNrOfAngles() > 0);

  SECTION("Nrs of chains is correct")
  {
    REQUIRE(universe.getAtomsWithType(2).size() == 100);
    REQUIRE(universe.getAtomsWithType(1).size() == (4 / 2) * 100 * 16);
    REQUIRE(universe.getMolecules(2).size() == (4 / 2) * 100 + 100);
  }

  SECTION("Universe is generated deterministically")
  {
    pu::MCUniverseGenerator generator2 =
      pu::MCUniverseGenerator(10.0, 10.0, 10.0);
    generator2.setSeed(8804);
    generator2.setBeadDistance(0.964);
    generator2.addCrosslinkers(100, 2);
    generator2.addSolventChains(100, 16, 3);
    generator2.addAndLinkStrands((4 / 2) * 100, 16, 0.8);

    pe::Universe universe2 = generator2.getUniverse();

    REQUIRE(universe.getNrOfAtoms() == universe2.getNrOfAtoms());
    REQUIRE(universe.getNrOfBonds() == universe2.getNrOfBonds());
    REQUIRE(universe.getAtom(3) == universe2.getAtom(3));
  }

  SECTION("Errors are thrown")
  {
    // only one type of cross-linker supported yet
    REQUIRE_THROWS(generator.addCrosslinkers(100, 3));
    // nr of strands and strand lengths must be same:
    REQUIRE_THROWS(generator.addAndLinkStrands(3, { { 10, 100 } }, 0.1, 4, 1));
    // strands may only be added once:
    REQUIRE_THROWS(generator.addAndLinkStrands((4 / 2) * 100, 16, 0.8));
    // not enough strands to reach conversion:
    REQUIRE_THROWS(generator.addAndLinkStrands(2, 10, 1.0, 4, 1));
  }

  // SECTION("Universe can be written and read again") {
  //   pu::DataFileWriter writer = pu::DataFileWriter(universe);
  //   writer.configIncludeAngles(true);
  //   std::string file = "tmp_data_file_with_mc.structure.out";
  //   writer.writeToFile(file);
  //   pe::UniverseSequence seq = pe::UniverseSequence();
  //   seq.initializeFromDataSequence({{file}});
  //   pe::Universe readUniverse = seq.atIndex(0);

  //   REQUIRE(universe.getNrOfAtoms() == readUniverse.getNrOfAtoms());
  //   REQUIRE(universe.getNrOfBonds() == readUniverse.getNrOfBonds());
  //   REQUIRE(universe.getNrOfAngles() == readUniverse.getNrOfAngles());
  // }
}

TEST_CASE("Large Universe can be generated", "[generator][MCUniverseGenerator]")
{

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(1200, 2);

  pe::Universe universe = generator.getUniverse();
  REQUIRE(universe.getVolume() == 10.0 * 10.0 * 10.0);
  REQUIRE(universe.getAtomsWithType(2).size() == 1200);
}
