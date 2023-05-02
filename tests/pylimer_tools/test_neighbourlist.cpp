#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/EigenNeighbourList.h"
#include "../../src/pylimer_tools_cpp/entities/NeighbourList.h"
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

TEST_CASE("NeighbourList works as intended", "[entity][NeighbourList]")
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

  pe::NeighbourList neighbourList =
    pe::NeighbourList(universe.getAtoms(), universe.getBox(), 3.0);

  std::vector<pe::Atom> neighbours =
    neighbourList.getAtomsCloseTo(universe.getAtom(10000), 0.000001);
  REQUIRE(neighbours.size() == 0);

  neighbours =
    neighbourList.getAtomsCloseTo(universe.getAtom(10000), 1.0, 0.9999);
  REQUIRE(neighbours.size() == 0);

  neighbours = neighbourList.getAtomsCloseTo(universe.getAtom(20000));
  REQUIRE(neighbours.size() > 0);

  REQUIRE_THROWS(
    neighbourList.getAtomsCloseTo(universe.getAtom(10000), 2.0, 3.0));

  // remove
  neighbourList.removeAtom(universe.getAtom(10000));
  // make sure we cannot query the remove atom
  REQUIRE_THROWS(neighbourList.getAtomsCloseTo(universe.getAtom(10000)));
  // nor find it in other neighbour lists
  std::vector<pe::Atom> neighbours2 =
    neighbourList.getAtomsCloseTo(universe.getAtom(20000));
  REQUIRE(neighbours2.size() == neighbours.size() - 1);
}

TEST_CASE("Manually accurate NeighbourList", "[entity][NeighbourList]")
{
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  const pe::Box box = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
  universe.setBox(box);
  universe.addAtoms(8,
                    { { 1, 2, 3, 4, 5, 6, 7, 8 } },     // id
                    { { 2, 1, 1, 1, 2, 1, 1, 1 } },     // type
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // x
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // y
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // z
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } },     // nx
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } },     // ny
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } }      // nz
  );

  pe::NeighbourList neighbourList =
    pe::NeighbourList(universe.getAtoms(), universe.getBox(), 3.0);

  std::vector<pe::Atom> neighbours =
    neighbourList.getAtomsCloseTo(universe.getAtom(2), 1.0);
  REQUIRE(neighbours.size() == 0);
  neighbours = neighbourList.getAtomsCloseTo(universe.getAtom(2), 2.0);
  REQUIRE(neighbours.size() == 2);
  neighbours = neighbourList.getAtomsCloseTo(universe.getAtom(2), 2.0, 1.8);
  REQUIRE(neighbours.size() == 0);
  neighbours = neighbourList.getAtomsCloseTo(universe.getAtom(1), 2.0);
  REQUIRE(neighbours.size() == 1);
}

TEST_CASE("Manually accurate EigenNeighbourList",
          "[entity][EigenNeighbourList]")
{
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  const pe::Box box = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
  universe.setBox(box);
  universe.addAtoms(8,
                    { { 1, 2, 3, 4, 5, 6, 7, 8 } },     // id
                    { { 2, 1, 1, 1, 2, 1, 1, 1 } },     // type
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // x
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // y
                    { { 1, 2, 3, 4, 9, -10, -9, -8 } }, // z
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } },     // nx
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } },     // ny
                    { { 0, 0, 0, 0, 1, 2, 2, 2 } }      // nz
  );

  Eigen::VectorXd coordinates = universe.getUnwrappedVertexCoordinates(&box);
  pe::EigenNeighbourList neighbourList =
    pe::EigenNeighbourList(coordinates, universe.getBox(), 3.0);

  Eigen::ArrayXi neighbours = neighbourList.getIndicesCloseToCoordinates(
    coordinates.segment(3 * universe.getIdxByAtomId(2), 3), 1.0);
  REQUIRE(neighbours.size() == 0);
  neighbours = neighbourList.getIndicesCloseToCoordinates(
    coordinates.segment(3 * universe.getIdxByAtomId(2), 3), 2.0);
  REQUIRE(neighbours.size() == 2);
  // neighbours = neighbourList.getIndicesCloseToCoordinates(
  //   coordinates.segment(3 * universe.getIdxByAtomId(2), 3), 2.0, 1.8);
  // REQUIRE(neighbours.size() == 0);
  neighbours = neighbourList.getIndicesCloseToCoordinates(
    coordinates.segment(3 * universe.getIdxByAtomId(1), 3), 2.0);
  REQUIRE(neighbours.size() == 1);
}
