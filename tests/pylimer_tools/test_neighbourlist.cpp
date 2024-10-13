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
  std::cout << "Running test \"NeighbourList works as intended\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";
  CHECK(std::filesystem::exists(suspectedPath));

  universeSeq.initializeFromDataSequence(
    { { suspectedPath + "lammps_data_file_small.out" } });
  CHECK(universeSeq.getLength() == 1);
  CHECK(universeSeq.atIndex(0).getNrOfAtoms() == 12);
  CHECK(universeSeq.atIndex(0).getNrOfBonds() == 5);

  pe::Universe universe = universeSeq.atIndex(0);

  pe::NeighbourList neighbourList =
    pe::NeighbourList(universe.getAtoms(), universe.getBox(), 3.0);

  pe::Atom testAtom10000 = universe.getAtom(10000);
  std::vector<pe::Atom> neighbours =
    neighbourList.getAtomsCloseTo(testAtom10000, 0.000001);
  CHECK(neighbours.size() == 0);

  neighbours = neighbourList.getAtomsCloseTo(testAtom10000, 1.0, 0.9999);
  CHECK(neighbours.size() == 0);

  pe::Atom testAtom20000 = universe.getAtom(20000);
  neighbours = neighbourList.getAtomsCloseTo(testAtom20000);
  CHECK(neighbours.size() > 0);

  CHECK_THROWS(neighbourList.getAtomsCloseTo(testAtom10000, 2.0, 3.0, true));

  // remove
  neighbourList.removeAtom(testAtom10000);
  // make sure we cannot query the remove atom
  CHECK_THROWS(neighbourList.getAtomsCloseTo(testAtom10000, 3.0, 0.0, false, true));
  // nor find it in other neighbour lists
  std::vector<pe::Atom> neighbours2 =
    neighbourList.getAtomsCloseTo(testAtom20000);
  CHECK(neighbours2.size() == neighbours.size() - 1);
}

TEST_CASE("Manually accurate NeighbourList", "[entity][NeighbourList]")
{
  std::cout << "Running test \"Manually accurate NeighbourList\"" << std::endl;
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
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

  pe::NeighbourList neighbourList =
    pe::NeighbourList(universe.getAtoms(), universe.getBox(), 3.0);

  pe::Atom atom2 = universe.getAtom(2);
  std::vector<pe::Atom> neighbours = neighbourList.getAtomsCloseTo(atom2, 1.0);
  CHECK(neighbours.size() == 0);
  neighbours = neighbourList.getAtomsCloseTo(atom2, 2.0);
  CHECK(neighbours.size() == 2);
  neighbours = neighbourList.getAtomsCloseTo(atom2, 2.0, 1.8);
  CHECK(neighbours.size() == 0);
  pe::Atom atom3 = universe.getAtom(3);
  neighbourList.removeAtom(atom3);
  neighbours = neighbourList.getAtomsCloseTo(atom2, 2.0);
  CHECK(neighbours.size() == 1);
}

TEST_CASE("Random coordinates EigenNeighbourList",
          "[entity][EigenNeighbourList][long]")
{
  std::cout << "Running test \"Random coordinates EigenNeighbourList\""
            << std::endl;
  int numAtoms = 8300;
  Eigen::VectorXd coordinates = Eigen::VectorXd::Random(numAtoms * 3) * 100.;
  pe::Box box = pe::Box(14, 14, 14);
  pe::EigenNeighbourList neighbourList =
    pe::EigenNeighbourList(coordinates, box, 1.0, 1.0);
  double cutoff = 1.0;
  // pre-allocate the neighbor indices array
  Eigen::ArrayXi neighbors = Eigen::ArrayXi(static_cast<int>(
    numAtoms * (std::ceil((3.1 * cutoff) * (3.1 * cutoff) * (3.1 * cutoff)) /
                box.getVolume())));

  for (int i = 0; i < numAtoms; ++i) {
    int numNeighbors = neighbourList.getIndicesCloseToCoordinates(
      neighbors, coordinates.segment(3 * i, 3), cutoff);

    std::vector<size_t> relevantNeighbors;
    std::vector<size_t> relevantPairs;
    for (size_t neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
      const size_t j = neighbors[neigh_idx];
      if (j <= i) {
        continue;
      }
      Eigen::Vector3d pairdistance =
        coordinates.segment(3 * i, 3) - coordinates.segment(3 * j, 3);
      box.handlePBC(pairdistance);
      const double rNorm = pairdistance.norm();
      if (rNorm >= cutoff || rNorm < 1e-12) {
        continue;
      }

      relevantNeighbors.push_back(j);
    }
    for (size_t j = i + 1; j < numAtoms; ++j) {
      Eigen::Vector3d pairdistance =
        coordinates.segment(3 * i, 3) - coordinates.segment(3 * j, 3);
      box.handlePBC(pairdistance);
      const double rNorm = pairdistance.norm();
      if (rNorm >= cutoff || rNorm < 1e-12) {
        continue;
      }

      relevantPairs.push_back(j);
    }

    CHECK(relevantPairs.size() == relevantNeighbors.size());
    std::sort(relevantPairs.begin(), relevantPairs.end());
    std::sort(relevantNeighbors.begin(), relevantNeighbors.end());

    CHECK(relevantNeighbors == relevantPairs);
  }
  CHECK(true);
}

TEST_CASE("Manually accurate EigenNeighbourList",
          "[entity][EigenNeighbourList]")
{
  std::cout << "Running test \"Manually accurate EigenNeighbourList\""
            << std::endl;
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
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

  Eigen::VectorXd coordinates = universe.getUnwrappedVertexCoordinates(box);
  pe::EigenNeighbourList neighbourList =
    pe::EigenNeighbourList(coordinates, universe.getBox(), 2.0, 1.0);

  Eigen::ArrayXi neighbours = neighbourList.getIndicesCloseToCoordinates(
    coordinates.segment(3 * universe.getIdxByAtomId(2), 3), 1.0);
  CHECK(neighbours.size() == 3);
  neighbours = neighbourList.getIndicesCloseToCoordinates(
    coordinates.segment(3 * universe.getIdxByAtomId(2), 3), 2.0);
  CHECK(neighbours.size() == 4);
  // neighbours = neighbourList.getIndicesCloseToCoordinates(
  //   coordinates.segment(3 * universe.getIdxByAtomId(2), 3), 2.0, 1.8);
  // CHECK(neighbours.size() == 0);
  neighbours = neighbourList.getIndicesCloseToCoordinates(
    coordinates.segment(3 * universe.getIdxByAtomId(1), 3), 2.0);
  CHECK(neighbours.size() == 3);

  std::vector<std::vector<pe::bucket_idx_t>> buckets =
    neighbourList.getNeighbourBuckets();
  Eigen::ArrayXi resizableResults = Eigen::ArrayXi(10);
  std::vector<pe::coordinate_idx_t> indices;
  Eigen::VectorXi bucketSizes = neighbourList.getNeighbourBucketSizes();
  for (size_t bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    CHECK(buckets[bucketIndex].size() == bucketSizes[bucketIndex]);
    Eigen::Vector3d centralCoordinates =
      neighbourList.getCentralCoordinatesOfBucket(bucketIndex);
    indices = neighbourList.getCombinedBucketIndicesForCoordinates(
      centralCoordinates, 0.01);
    CHECK(indices.size() == 1);
    CHECK(indices[0] == bucketIndex);
    // small enough cut-off to only query this one bucket
    Eigen::ArrayXi atomsInThisBucket =
      neighbourList.getIndicesCloseToCoordinates(centralCoordinates, 0.01);
    CHECK(atomsInThisBucket.size() == bucketSizes[bucketIndex]);
    int numA = neighbourList.getIndicesCloseToCoordinates(
      resizableResults, centralCoordinates, 0.01);
    CHECK(numA == bucketSizes[bucketIndex]);
    for (size_t i = 0; i < numA; ++i) {
      CHECK(atomsInThisBucket[i] == resizableResults[i]);
    }
  }
  CHECK(buckets.size() == 10 * 10 * 10);
  Eigen::Vector3d zeroVec = Eigen::Vector3d(0., 0., 0.);
  indices = neighbourList.getCombinedBucketIndicesForCoordinates(zeroVec, 1.9);
  CHECK(indices.size() == 8);
  indices = neighbourList.getCombinedBucketIndicesForCoordinates(zeroVec, 0.2);
  Eigen::Vector3d oneVec = Eigen::Vector3d(1., 1., 1.);
  CHECK(indices.size() == 8);
  indices = neighbourList.getCombinedBucketIndicesForCoordinates(oneVec, 2.0);
  CHECK(indices.size() == 9 * 3);
  indices = neighbourList.getCombinedBucketIndicesForCoordinates(oneVec, 0.2);
  CHECK(indices.size() == 1);
  CHECK(true);
}
