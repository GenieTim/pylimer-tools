#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/EigenNeighbourList.h"
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
  CHECK(neighbourList == neighbourList);
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
    for (Eigen::Index neigh_idx = 0;
         neigh_idx < static_cast<Eigen::Index>(numNeighbors);
         ++neigh_idx) {
      const int j = neighbors[neigh_idx];
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

      relevantNeighbors.push_back(static_cast<size_t>(j));
    }
    for (int j = i + 1; j < numAtoms; ++j) {
      Eigen::Vector3d pairdistance =
        coordinates.segment(3 * i, 3) - coordinates.segment(3 * j, 3);
      box.handlePBC(pairdistance);
      const double rNorm = pairdistance.norm();
      if (rNorm >= cutoff || rNorm < 1e-12) {
        continue;
      }

      relevantPairs.push_back(static_cast<size_t>(j));
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
    CHECK(
      buckets[bucketIndex].size() ==
      static_cast<size_t>(bucketSizes[static_cast<Eigen::Index>(bucketIndex)]));
    Eigen::Vector3d centralCoordinates =
      neighbourList.getCentralCoordinatesOfBucket(bucketIndex);
    indices = neighbourList.getCombinedBucketIndicesForCoordinates(
      centralCoordinates, 0.01);
    CHECK(indices.size() == 1);
    CHECK(static_cast<size_t>(indices[0]) == bucketIndex);
    // small enough cut-off to only query this one bucket
    Eigen::ArrayXi atomsInThisBucket =
      neighbourList.getIndicesCloseToCoordinates(centralCoordinates, 0.01);
    CHECK(atomsInThisBucket.size() ==
          bucketSizes[static_cast<Eigen::Index>(bucketIndex)]);
    int numA = neighbourList.getIndicesCloseToCoordinates(
      resizableResults, centralCoordinates, 0.01);
    CHECK(numA == bucketSizes[static_cast<Eigen::Index>(bucketIndex)]);
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numA); ++i) {
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

TEST_CASE("EigenNeighbourList Works as intended",
          "[entity][EigenNeighbourList]")
{
  std::cout << "Running test \"EigenNeighbourList Works as intended\""
            << std::endl;
  Eigen::VectorXd coordinates = Eigen::VectorXd::Zero(3 * 20);
  // one diagonal line of coordinates, breaking outside the box boundaries
  for (int i = 0; i < 20; ++i) {
    coordinates.segment(3 * i, 3) = Eigen::Vector3d(i, i, i);
  }
  pe::Box box = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);

  double distanceBetweenAtoms = std::sqrt(3.);

  pe::EigenNeighbourList neighbourList =
    pe::EigenNeighbourList(coordinates, box, 1.1 * distanceBetweenAtoms);

  Eigen::ArrayXi neighbours;
  for (int i = 0; i < 20; ++i) {
    INFO("Checking neighbours within default distance of atom " << i);
    neighbourList.getIndicesCloseToCoordinates(neighbours,
                                               coordinates.segment(3 * i, 3),
                                               1.1 * distanceBetweenAtoms,
                                               true);

    CHECK(neighbours.size() >= 3);
    // expect one of the neighbours to be itself, one before and one after
    bool foundItself = false;
    bool foundBefore = false;
    bool foundAfter = false;
    for (Eigen::Index j = 0; j < neighbours.size(); ++j) {
      foundItself = (neighbours[j] == i) || foundItself;
      foundBefore = (neighbours[j] == ((20 + (i - 1)) % 20)) || foundBefore;
      foundAfter = (neighbours[j] == ((i + 1) % 20)) || foundAfter;
    }
    CHECK(foundItself);
    CHECK(foundBefore);
    CHECK(foundAfter);

    INFO("Checking neighbours within new distance of atom " << i);
    neighbourList.getIndicesCloseToCoordinates(neighbours,
                                               coordinates.segment(3 * i, 3),
                                               1.3 * distanceBetweenAtoms,
                                               false);

    CHECK(neighbours.size() >= 3);
    // expect one of the neighbours to be itself, one before and one after
    foundItself = false;
    foundBefore = false;
    foundAfter = false;
    for (Eigen::Index j = 0; j < neighbours.size(); ++j) {
      foundItself = (neighbours[j] == i) || foundItself;
      foundBefore = (neighbours[j] == ((20 + (i - 1)) % 20)) || foundBefore;
      foundAfter = (neighbours[j] == ((i + 1) % 20)) || foundAfter;
    }
    CHECK(foundItself);
    CHECK(foundBefore);
    CHECK(foundAfter);
  }
}

#ifdef CEREALIZABLE
TEST_CASE("EigenNeighbourList can be serialized and deserialized",
          "[EigenNeighbourList][entity][serialization]")
{
  pe::Box box = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.);
  Eigen::VectorXd coordinates = Eigen::VectorXd::Random(3 * 10);
  pe::EigenNeighbourList neighbourList =
    pe::EigenNeighbourList(coordinates, box, M_PI);
  std::string neighbourListStr =
    pylimer_tools::utils::serializeToString(neighbourList);
  pe::EigenNeighbourList deserializedNeighbourList;
  pylimer_tools::utils::deserializeFromString(deserializedNeighbourList,
                                              neighbourListStr);
  CHECK(neighbourList == deserializedNeighbourList);

  Eigen::VectorXd coordinates2 = Eigen::VectorXd::Random(3 * 10);
  neighbourList.resetCoordinates(coordinates2);
  CHECK(neighbourList != deserializedNeighbourList);
}
#endif

TEST_CASE("EigenNeighbourList diagnostic validateWhyNotIncluded works",
          "[EigenNeighbourList][entity]")
{
  pe::Box box = pe::Box(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
  Eigen::VectorXd coordinates = Eigen::VectorXd::Zero(6);
  coordinates << 1.0, 1.0, 1.0, 8.0, 8.0, 8.0;
  pe::EigenNeighbourList neighbourList =
    pe::EigenNeighbourList(coordinates, box, 2.0);

  // Default cutoff, target outside and inside
  CHECK_NOTHROW(neighbourList.validateWhyNotIncluded(
    Eigen::Vector3d(1.0, 1.0, 1.0), Eigen::Vector3d(8.0, 8.0, 8.0), -1.0));
  CHECK_NOTHROW(neighbourList.validateWhyNotIncluded(
    Eigen::Vector3d(1.0, 1.0, 1.0), Eigen::Vector3d(1.5, 1.5, 1.5), 2.0));

  // Custom cutoff, target outside and inside
  CHECK_NOTHROW(neighbourList.validateWhyNotIncluded(
    Eigen::Vector3d(1.0, 1.0, 1.0), Eigen::Vector3d(8.0, 8.0, 8.0), 1.0));
  CHECK_NOTHROW(neighbourList.validateWhyNotIncluded(
    Eigen::Vector3d(1.0, 1.0, 1.0), Eigen::Vector3d(2.0, 2.0, 2.0), 1.5));
}

TEST_CASE("EigenNeighbourList handles small boxes and non-unit scaling factors "
          "without duplicates or misses",
          "[EigenNeighbourList][entity]")
{
  const pe::Box smallBox = pe::Box(0.0, 2.5, 0.0, 2.5, 0.0, 2.5);
  const int numAtoms = 40;
  Eigen::VectorXd coordinates = Eigen::VectorXd::Zero(3 * numAtoms);
  for (int i = 0; i < numAtoms; ++i) {
    coordinates[3 * i + 0] = std::fmod(0.37 * i, 2.5);
    coordinates[3 * i + 1] = std::fmod(0.61 * i, 2.5);
    coordinates[3 * i + 2] = std::fmod(0.83 * i, 2.5);
  }

  for (const double scalingFactor : { 1.0, 1.5, 3.0 }) {
    DYNAMIC_SECTION("scalingFactor = " << scalingFactor)
    {
      const double cutoff = 1.1;
      const pe::EigenNeighbourList constList(
        coordinates, smallBox, cutoff, scalingFactor);
      CHECK(constList.checkIfCoordinatesAreCurrent(coordinates));
      CHECK(constList.getNumBinnedCoordinates() == numAtoms);

      Eigen::VectorXd perturbed = coordinates;
      perturbed[0] += 1.25;
      CHECK_FALSE(constList.checkIfCoordinatesAreCurrent(perturbed));

      Eigen::ArrayXi neighbours(4);
      for (int i = 0; i < numAtoms; ++i) {
        const Eigen::Vector3d query = coordinates.segment<3>(3 * i);
        const int nFound = constList.getIndicesCloseToCoordinates(
          neighbours, query, cutoff, true);

        // Must never return duplicate atom indices on small periodic grids
        std::vector<int> returned(neighbours.data(),
                                  neighbours.data() + nFound);
        std::sort(returned.begin(), returned.end());
        CHECK(std::adjacent_find(returned.begin(), returned.end()) ==
              returned.end());

        // Must recall every true PBC neighbour within cutoff
        for (int j = 0; j < numAtoms; ++j) {
          Eigen::Vector3d diff = query - coordinates.segment<3>(3 * j);
          smallBox.handlePBC(diff);
          if (diff.norm() < cutoff) {
            CHECK(std::binary_search(returned.begin(), returned.end(), j));
          }
        }
      }
    }
  }

  SECTION("Single-bucket and empty neighbour list preserve preallocated buffer")
  {
    const pe::EigenNeighbourList emptyList(
      Eigen::VectorXd(0), pe::Box(1.0, 1.0, 1.0), 2.0, 1.0);
    Eigen::ArrayXi preallocated(16);
    const int count = emptyList.getIndicesCloseToCoordinates(
      preallocated, Eigen::Vector3d::Zero(), 2.0, true);
    CHECK(count == 0);
    CHECK(preallocated.size() == 16);
  }
}
