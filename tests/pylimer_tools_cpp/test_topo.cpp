#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/topo/EntanglementDetector.h"
#include "../../src/pylimer_tools_cpp/utils/RandomWalker.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <iostream>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
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

TEST_CASE("Entanglement Detector respects distance cut-off",
          "[analysis][EntanglementDetector][long]")
{
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::vector<std::string> files = {
    "3d-diamond-lattice_10x10x10_a_3_d_0.85_imperfect.structure.out",
    "crosslinked_M10000_N39_p_0.9.out",
    "crosslinked_p_0.98_melt_100_a_38_50_xlinks_v_22.structure.out-equilibration_do_crosslink.structure.out",
    "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out"
  };

  for (const std::string& file : files) {
    pe::UniverseSequence universeSeq = pe::UniverseSequence();
    REQUIRE(universeSeq.getLength() == 0);
    std::string inputFile = suspectedPath + "/structure/" + file;
    if (!std::filesystem::exists(inputFile)) {
      std::cerr << "File not found: " << inputFile << std::endl;
      continue;
    }

    std::cout << "Processing file: " << file << std::endl;

    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);

    for (const bool filtered : {
      true,
      false
    }) {
      std::cout << "Filtering: " << (filtered ? "true" : "false") << std::endl;
      for (const double lowerCutoff : {
             0.,
             1.,
           }) {
        std::cout << "\tLower cutoff: " << lowerCutoff << std::endl;
        for (const double upperCutoff : {
               2.,
               5.,
             }) {
          std::cout << "\t\tUpper cutoff: " << upperCutoff << std::endl;
          for (const double sameStrandCutoff : { -1., 2., 5. }) {
            INFO("Processing file: " << inputFile << " with lowerCutoff: "
                                     << lowerCutoff
                                     << ", upperCutoff: " << upperCutoff
                                     << ", sameStrandCutoff: "
                                     << sameStrandCutoff
                                     << ", filtered: " << filtered);
            std::cout << "\t\t\tSame strand cutoff: " << sameStrandCutoff;
            auto start_ref = std::chrono::high_resolution_clock::now();

            pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
              entanglements = pylimer_tools::topo::entanglement_detection::
                randomlyFindEntanglementsV2(universe,
                                            0.1 * universe.getNrOfAtoms(),
                                            upperCutoff,
                                            lowerCutoff,
                                            0,
                                            sameStrandCutoff,
                                            "",
                                            2,
                                            true,
                                            filtered);

            auto end_ref = std::chrono::high_resolution_clock::now();
            auto duration_ref =
              std::chrono::duration_cast<std::chrono::microseconds>(end_ref -
                                                                    start_ref);
            std::cout << "\tEntanglements v1: "
                      << std::duration_to_string(duration_ref) << " "
                      << std::endl;

            CHECK(entanglements.pairsOfAtoms.size() >=
                  0.05 * universe.getNrOfAtoms());
            CHECK(entanglements.pairsOfAtoms.size() <=
                  0.11 * universe.getNrOfAtoms());

            for (auto& pair : entanglements.pairsOfAtoms) {
              pe::Atom a1 = universe.getAtom(pair.first);
              pe::Atom a2 = universe.getAtom(pair.second);
              CHECK(a1.distanceTo(a2, universe.getBox()) <= upperCutoff);
              CHECK(a1.distanceTo(a2, universe.getBox()) >= lowerCutoff);
            }
          }
        }
      }
    }
  }
}
