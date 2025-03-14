#include "../../src/pylimer_tools_cpp/topo/EntanglementDetector.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <iostream>

namespace pe = pylimer_tools::entities;
namespace pt = pylimer_tools::topo::entanglement_detection;

TEST_CASE("Entanglement Detector can ignore soluble chains",
          "[analysis][EntanglementDetector]")
{
  size_t nrOfBeads = 100;
  size_t nrOfBeadsPerChain = 5;
  pe::Universe universe = pe::Universe(15.0, 15.0, 15.0);
  std::vector<double> xPositions;
  xPositions.reserve(nrOfBeads);
  std::vector<double> yPositions;
  yPositions.reserve(nrOfBeads);
  std::vector<double> zPositions;
  zPositions.reserve(nrOfBeads);
  std::vector<long int> atomIds;
  atomIds.reserve(nrOfBeads);
  std::vector<int> atomTypes;
  atomTypes.reserve(nrOfBeads);
  std::vector<int> zeroInts;
  zeroInts.reserve(nrOfBeads);
  std::vector<long int> bondFrom;
  bondFrom.reserve(nrOfBeads - 1);
  std::vector<long int> bondTo;
  bondTo.reserve(nrOfBeads - 1);
  double offset = 10.0;
  for (int i = 0; i < nrOfBeads; ++i) {
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i + 1);
    atomTypes.push_back(i % nrOfBeadsPerChain == 0 ? 2 : 1);
    zeroInts.push_back(0);
    if (i > 0 && !(i % nrOfBeadsPerChain == 0)) {
      bondFrom.push_back(i);
      bondTo.push_back(i + 1);
    }
  }
  universe.addAtoms(atomIds,
                    atomTypes,
                    xPositions,
                    yPositions,
                    zPositions,
                    zeroInts,
                    zeroInts,
                    zeroInts);
  universe.addBonds(bondFrom, bondTo);
  CHECK(universe.getNrOfAtoms() == nrOfBeads);

  SECTION("Ignoring soluble & dangling chains")
  {
    pt::AtomPairEntanglements es = pt::randomlyFindEntanglements(
      universe, 10, 10.0, 0., 0, 0, "seed123987", 2, true, true);
    CHECK(es.pairsOfAtoms.size() == 0);
  }

  SECTION("Allowing soluble & dangling chains")
  {
    pt::AtomPairEntanglements es = pt::randomlyFindEntanglements(
      universe, 10, 10.0, 0., 0, 0, "seed123987", 2, true, false);
    CHECK(es.pairsOfAtoms.size() == 10);
  }
}

TEST_CASE("Entanglement Detector throws", "[analysis][EntanglementDetector]")
{
  pe::Universe universe = pe::Universe(15.0, 15.0, 15.0);
  CHECK_THROWS(pt::randomlyFindEntanglements(universe, 10, 1., 2.));
}
