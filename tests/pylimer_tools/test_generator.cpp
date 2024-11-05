#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/io/DataFileWriter.h"
#include "../../src/pylimer_tools_cpp/utils/MCUniverseGenerator.h"
#include <catch2/catch_approx.hpp>
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

TEST_CASE("Certain configurations do not lead to memory corruption",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"Certain configurations do not lead to memory "
               "corruption\""
            << std::endl;
  // the following parameters have led to a `double free or corruption` error?!?
  int nrOfCrosslinkers = static_cast<int>(5e4 * 2 * 0.7 / 7);
  double sideLength = std::cbrt((10 * 5e4 * nrOfCrosslinkers) / 0.85);
  pu::MCUniverseGenerator generator =
    pu::MCUniverseGenerator(sideLength, sideLength, sideLength);
  REQUIRE_NOTHROW(generator.setSeed(68419));
  REQUIRE_NOTHROW(generator.setBeadDistance(0.965));

  pe::Universe universe = generator.getUniverse();
  REQUIRE(universe.getNrOfAtoms() == 0);
}

TEST_CASE("Universe can be generated", "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"Universe can be generated\"" << std::endl;
  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(100, 4, 2);
  generator.addSolventChains(100, 16, 3);
  generator.addAndLinkStrandsToConversion((4 / 2) * 100, 16, 0.8);
  generator.configNrOfMCSteps(100);

  pe::Universe universe = generator.getUniverse();
  std::map<int, double> weights;
  weights[1] = 1.0;
  weights[2] = 1.0;
  weights[3] = 1.0;
  universe.setMasses(weights);
  REQUIRE(universe.getVolume() == 10.0 * 10.0 * 10.0);

  auto angles = universe.detectAngles();
  std::vector<int> angleTypes;
  angleTypes.reserve(angles["angle_from"].size());
  for (size_t i = 0; i < angles["angle_from"].size(); i++) {
    angleTypes.push_back(1);
  }
  universe.addAngles(
    angles["angle_from"], angles["angle_via"], angles["angle_to"], angleTypes);
  REQUIRE(universe.getNrOfAngles() > 0);

  SECTION("Nrs of chains is correct")
  {
    REQUIRE(universe.getAtomsOfType(2).size() == 100);
    REQUIRE(universe.getAtomsOfType(1).size() == (4 / 2) * 100 * 16);
    REQUIRE(universe.getMolecules(2).size() == (4 / 2) * 100 + 100);
  }

  SECTION("Universe is generated deterministically")
  {
    pu::MCUniverseGenerator generator2 =
      pu::MCUniverseGenerator(10.0, 10.0, 10.0);
    generator2.setSeed(8804);
    generator2.setBeadDistance(0.964);
    generator2.setMeanSquaredBeadDistance(0.964 * 0.964);
    generator2.addCrosslinkers(100, 4, 2);
    generator2.addSolventChains(100, 16, 3);
    generator2.addAndLinkStrandsToConversion((4 / 2) * 100, 16, 0.8);

    pe::Universe universe2 = generator2.getUniverse();

    REQUIRE(universe.getNrOfAtoms() == universe2.getNrOfAtoms());
    REQUIRE(universe.getNrOfBonds() == universe2.getNrOfBonds());
    REQUIRE(universe.getAtom(3) == universe2.getAtom(3));

    auto molecules = universe.getMolecules(2);
    CHECK(molecules.size() == (4 / 2) * 100 + 100);
    for (const auto& molecule : molecules) {
      CHECK(molecule.getNrOfAtoms() == 16);
    }

    auto bondLengths = universe.computeBondLengths();
    for (const double bondLength : bondLengths) {
      CHECK(bondLength > 0.0);
      CHECK(bondLength < 3.0);
    }
  }

  SECTION("Errors are thrown")
  {
    // nr of strands and strand lengths must be same:
    REQUIRE_THROWS(
      generator.addAndLinkStrandsToConversion(3, { { 10, 100 } }, 0.1, 1));
    // not enough strands to reach conversion:
    REQUIRE_THROWS(generator.addAndLinkStrandsToConversion(2, 10, 1.0, 1));
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
  std::cout << "Running test \"Large Universe can be generated\"" << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(1200, 2);

  pe::Universe universe = generator.getUniverse();
  REQUIRE(universe.getVolume() == 10.0 * 10.0 * 10.0);
  REQUIRE(universe.getAtomsOfType(2).size() == 1200);
}

TEST_CASE("MCUniverseGenerator knows about <b> vs. <b^2>",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"MCUniverseGenerator knows about <b> vs. <b^2>\""
            << std::endl;
  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);

  generator.setBeadDistance(2.);
  CHECK(generator.getConfiguredBeadDistance() == Catch::Approx(2.));
  CHECK(generator.getConfiguredMeanSquaredBeadDistance() ==
        Catch::Approx((2. * 2.) * (3. * M_PI / 8.)));

  generator.setMeanSquaredBeadDistance(2.0);
  CHECK(generator.getConfiguredMeanSquaredBeadDistance() == Catch::Approx(2.));
  CHECK(generator.getConfiguredBeadDistance() ==
        Catch::Approx(std::sqrt(2. / ((3. * M_PI) / 8.))));

  generator.setBeadDistance(0.5, false);
  CHECK(generator.getConfiguredMeanSquaredBeadDistance() == Catch::Approx(2.));
  CHECK(generator.getConfiguredBeadDistance() == Catch::Approx(0.5));

  generator.setMeanSquaredBeadDistance(0.25, false);
  CHECK(generator.getConfiguredMeanSquaredBeadDistance() ==
        Catch::Approx(0.25));
  CHECK(generator.getConfiguredBeadDistance() == Catch::Approx(0.5));

  CHECK_THROWS(generator.setMeanSquaredBeadDistance(-1.));
  CHECK_THROWS(generator.setBeadDistance(-1.));
}

TEST_CASE("MCUniverseGenerator can generate without primary loops",
          "[generator][MCUniverseGenerator]")
{
  std::cout
    << "Running test \"MCUniverseGenerator can generate without primary loops\""
    << std::endl;
  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(100, 4, 2);
  generator.configPrimaryLoopProbability(0.);
  generator.addAndLinkStrandsToConversion(200, 10, 0.9, 1, 1.);

  pe::Universe universe = generator.getUniverse();
  CHECK(universe.getAtomsOfType(2).size() == 100);
  CHECK(universe.getMolecules(2).size() == 200);

  // verify that no primary loops are present
  // method 1: detection of molecules
  auto chains = universe.getChainsWithCrosslinker(2);
  for (const auto& chain : chains) {
    CHECK(chain.getType() != pe::MoleculeType::PRIMARY_LOOP);
  }

  // method 2: detection of loops
  auto loops = universe.countLoopLengths(14);
  CHECK(loops.size() == 0);
}

TEST_CASE("Universe can cross-link up to w_sol",
          "[generator][MCUniverseGenerator][long]")
{
  std::cout << "Running test \"Universe can cross-link up to w_sol\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(100, 4, 2);

  generator.addAndLinkStrandsToSolubleFraction(200, 19, 0.9, 1, 1.);

  pe::Universe universe = generator.getUniverse();
  REQUIRE(universe.getVolume() == 10.0 * 10.0 * 10.0);
  REQUIRE(universe.getAtomsOfType(2).size() == 100);
  REQUIRE(universe.getAtomsOfType(1).size() == 200 * 19);

  auto clusters = universe.getClusters();
  CHECK(clusters.size() < 20);
}
