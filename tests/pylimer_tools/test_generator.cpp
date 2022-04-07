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

TEST_CASE("Universe can be generated", "[generator][MCUniverseGenerator]")
{
  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(100, 2);
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
