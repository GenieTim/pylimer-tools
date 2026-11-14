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

  pylimer_tools::sim::mehp::CrosslinkerUniverse crosslinkerUniverse;
  pe::Box box(15.0, 15.0, 15.0);
  CHECK_THROWS(pt::randomlyFindEntanglementsAnalytical(
    box, crosslinkerUniverse, 10, 1.0, 2.0));
  CHECK_THROWS(pt::randomlyFindEntanglementsAnalytical(
    box, crosslinkerUniverse, 5, 2.0, 0.0, 10));
}

TEST_CASE("Entanglement Detector enforces sameStrandCutoff on rings and "
          "crosslinker-adjacent beads",
          "[analysis][EntanglementDetector]")
{
  SECTION("Isolated cyclic ring of 6 beads")
  {
    pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
    std::vector<long int> ids = { 1, 2, 3, 4, 5, 6 };
    std::vector<int> types(6, 1);
    std::vector<double> coords(6, 1.0);
    std::vector<int> zeros(6, 0);
    universe.addAtoms(ids, types, coords, coords, coords, zeros, zeros, zeros);
    universe.addBonds({ 1, 2, 3, 4, 5, 6 }, { 2, 3, 4, 5, 6, 1 });

    // Max distance on a 6-ring is 3 bonds
    pt::AtomPairEntanglements es3 = pt::randomlyFindEntanglements(
      universe, 1, 5.0, 0.0, 1, 3.0, "ring_seed", 2, true, false);
    REQUIRE(es3.pairsOfAtoms.size() == 1);
    const size_t v1 = universe.getIdxByAtomId(es3.pairsOfAtoms[0].first);
    const size_t v2 = universe.getIdxByAtomId(es3.pairsOfAtoms[0].second);
    CHECK(universe.getPathLength(static_cast<long int>(v1),
                                 static_cast<long int>(v2)) >= 3);

    // No pair has distance >= 4 on a 6-ring
    pt::AtomPairEntanglements es4 = pt::randomlyFindEntanglements(
      universe, 1, 5.0, 0.0, 0, 4.0, "ring_seed", 2, true, false);
    CHECK(es4.pairsOfAtoms.empty());
  }

  SECTION("Crosslinker-adjacent beads across a primary loop / junction")
  {
    // Atom 1 is a degree-4 crosslinker with two 2-bead primary loops:
    // 1 - 2 - 3 - 1 and 1 - 4 - 5 - 1
    // Every bead (2, 3, 4, 5) has distToEnd == 1, so distToEnd[u] +
    // distToEnd[v] == 2 < 3, forcing the bounded BFS check for
    // sameStrandCutoff=3.
    pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
    std::vector<long int> ids = { 1, 2, 3, 4, 5 };
    std::vector<int> types = { 2, 1, 1, 1, 1 };
    std::vector<double> coords(5, 1.0);
    std::vector<int> zeros(5, 0);
    universe.addAtoms(ids, types, coords, coords, coords, zeros, zeros, zeros);
    universe.addBonds({ 1, 2, 3, 1, 4, 5 }, { 2, 3, 1, 4, 5, 1 });

    // All pairs among {2, 3, 4, 5} have path length <= 2 through atom 1
    pt::AtomPairEntanglements es3 = pt::randomlyFindEntanglements(
      universe, 1, 5.0, 0.0, 0, 3.0, "loop_seed", 2, true, false);
    CHECK(es3.pairsOfAtoms.empty());

    // With sameStrandCutoff = 2.0, pairs across different loops (e.g. 2 and 4)
    // have path length 2 >= 2.0 and can pair, while same-loop pairs (2 and 3)
    // have path length 1 < 2.0 and cannot.
    pt::AtomPairEntanglements es2 = pt::randomlyFindEntanglements(
      universe, 1, 5.0, 0.0, 1, 2.0, "loop_seed", 2, true, false);
    REQUIRE(es2.pairsOfAtoms.size() == 1);
    for (const std::pair<size_t, size_t>& pair : es2.pairsOfAtoms) {
      const long int v1 =
        static_cast<long int>(universe.getIdxByAtomId(pair.first));
      const long int v2 =
        static_cast<long int>(universe.getIdxByAtomId(pair.second));
      CHECK(universe.getPathLength(v1, v2) >= 2);
    }
  }
}

TEST_CASE("Analytical Entanglement Detector samples valid entanglements on "
          "CrosslinkerUniverse",
          "[analysis][EntanglementDetector][analytical]")
{
  pe::Box box(20.0, 20.0, 20.0);
  pylimer_tools::sim::mehp::CrosslinkerUniverse cu;

  // Create 8 crosslinkers in a cube and connect them with active strands,
  // plus some dangling and soluble strands to test topological leaf peeling
  cu.xlinkTypes = {
    2, 2, 2, 2, 2, 2, 2, 2, 2
  }; // 8 active + 1 leaf crosslinker
  cu.xlinkX = { 5., 15., 5., 15., 5., 15., 5., 15., 10. };
  cu.xlinkY = { 5., 5., 15., 15., 5., 5., 15., 15., 10. };
  cu.xlinkZ = { 5., 5., 5., 5., 15., 15., 15., 15., 10. };
  cu.xlinkChainId = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  cu.strandsOfXlink.resize(9);

  auto addStrand = [&](long int from, long int to, int nBeads) {
    size_t sIdx = cu.strandFrom.size();
    cu.strandFrom.push_back(from);
    cu.strandTo.push_back(to);
    cu.beadsInStrand.push_back(nBeads);
    cu.strandBeadType.push_back(1);
    cu.beadDistanceInStrand.push_back(0.96);
    cu.meanSquaredBeadDistanceInStrand.push_back(1.0);
    if (from >= 0) {
      cu.strandsOfXlink[from].push_back(sIdx);
    }
    if (to >= 0 && to != from) {
      cu.strandsOfXlink[to].push_back(sIdx);
    }
  };

  // 12 edges of the cube (each crosslinker 0..7 has degree 3 -> active!)
  addStrand(0, 1, 25);
  addStrand(2, 3, 25);
  addStrand(4, 5, 25);
  addStrand(6, 7, 25);
  addStrand(0, 2, 25);
  addStrand(1, 3, 25);
  addStrand(4, 6, 25);
  addStrand(5, 7, 25);
  addStrand(0, 4, 25);
  addStrand(1, 5, 25);
  addStrand(2, 6, 25);
  addStrand(3, 7, 25);

  // Strand 12: bridge to leaf crosslinker 8 (degree 1 -> peeled away by leaf
  // peeling!)
  addStrand(0, 8, 25);
  // Strand 13: dangling strand attached to crosslinker 1
  addStrand(1, -1, 25);
  // Strand 14: free soluble strand
  addStrand(-1, -1, 25);

  SECTION("Respects topological filtering, cutoffs, and determinism")
  {
    const double lowerCutoff = 0.5;
    const double upperCutoff = 2.5;
    const double sameStrandCutoff = 4.0;

    pt::StrandPairEntanglements es1 =
      pt::randomlyFindEntanglementsAnalytical(box,
                                              cu,
                                              40,
                                              upperCutoff,
                                              lowerCutoff,
                                              10,
                                              sameStrandCutoff,
                                              "analytical_seed_42",
                                              true);
    pt::StrandPairEntanglements es2 =
      pt::randomlyFindEntanglementsAnalytical(box,
                                              cu,
                                              40,
                                              upperCutoff,
                                              lowerCutoff,
                                              10,
                                              sameStrandCutoff,
                                              "analytical_seed_42",
                                              true);

    REQUIRE(es1.size() >= 10);
    CHECK(es1.size() <= 40);
    // Determinism check
    REQUIRE(es1.size() == es2.size());
    CHECK(es1.strandIdx1 == es2.strandIdx1);
    CHECK(es1.strandIdx2 == es2.strandIdx2);
    CHECK(es1.beadIdx1 == es2.beadIdx1);
    CHECK(es1.beadIdx2 == es2.beadIdx2);
    CHECK(es1.x.isApprox(es2.x));
    CHECK(es1.y.isApprox(es2.y));
    CHECK(es1.z.isApprox(es2.z));

    for (size_t e = 0; e < es1.size(); ++e) {
      // Strands 12, 13, 14 are topologically inactive and must be filtered out
      CHECK(es1.strandIdx1[e] < 12);
      CHECK(es1.strandIdx2[e] < 12);
      CHECK(es1.alpha1[e] > 0.0);
      CHECK(es1.alpha1[e] < 1.0);
      CHECK(es1.alpha2[e] > 0.0);
      CHECK(es1.alpha2[e] < 1.0);
      if (es1.strandIdx1[e] == es1.strandIdx2[e]) {
        const double contourDiff =
          std::abs(static_cast<double>(es1.beadIdx1[e]) -
                   static_cast<double>(es1.beadIdx2[e]));
        CHECK(contourDiff >= sameStrandCutoff);
      }

      Eigen::Vector3d p1 = es1.position1.segment<3>(3 * e);
      Eigen::Vector3d p2 = es1.position2.segment<3>(3 * e);
      Eigen::Vector3d diff = p2 - p1;
      box.handlePBC(diff);
      CHECK(diff.norm() >= lowerCutoff - 1e-9);
      CHECK(diff.norm() <= upperCutoff + 1e-9);
    }

    // Conversion to AtomPairEntanglements
    pt::AtomPairEntanglements atomPairs = es1.toAtomPairEntanglements(cu);
    CHECK(atomPairs.pairsOfAtoms.size() == es1.size());
    for (const auto& pair : atomPairs.pairsOfAtoms) {
      CHECK(pair.first != pair.second);
      CHECK(atomPairs.pairOfAtom[pair.first - 1] >= 0);
      CHECK(atomPairs.pairOfAtom[pair.second - 1] >= 0);
      CHECK(atomPairs.pairOfAtom[pair.first - 1] ==
            atomPairs.pairOfAtom[pair.second - 1]);
    }
  }

  SECTION(
    "Purely soluble/dangling network yields 0 entanglements when filtered")
  {
    pylimer_tools::sim::mehp::CrosslinkerUniverse solubleCu;
    solubleCu.xlinkTypes = { 2, 2 };
    solubleCu.xlinkX = { 5., 7. };
    solubleCu.xlinkY = { 5., 7. };
    solubleCu.xlinkZ = { 5., 7. };
    solubleCu.xlinkChainId = { -1, -1 };
    solubleCu.strandsOfXlink = { { 0, 1 }, { 0 } };
    solubleCu.strandFrom = { 0, 0, -1 };
    solubleCu.strandTo = { 1, -1, -1 };
    solubleCu.beadsInStrand = { 20, 20, 20 };
    solubleCu.strandBeadType = { 1, 1, 1 };
    solubleCu.beadDistanceInStrand = { 1., 1., 1. };
    solubleCu.meanSquaredBeadDistanceInStrand = { 1., 1., 1. };

    pt::StrandPairEntanglements filtered =
      pt::randomlyFindEntanglementsAnalytical(
        box, solubleCu, 15, 4.0, 0.0, 0, 2.0, "seed_sol", true);
    CHECK(filtered.size() == 0);

    pt::StrandPairEntanglements unfiltered =
      pt::randomlyFindEntanglementsAnalytical(
        box, solubleCu, 15, 4.0, 0.0, 5, 2.0, "seed_sol", false);
    CHECK(unfiltered.size() >= 5);
  }
}

TEST_CASE("Universe toCrosslinkerUniverse and topological leaf peeling",
          "[analysis][Universe][CrosslinkerUniverse][leaf_peeling]")
{
  pe::Universe universe(30.0, 30.0, 30.0);
  std::vector<double> xPos;
  std::vector<double> yPos;
  std::vector<double> zPos;
  std::vector<long int> atomIds;
  std::vector<int> atomTypes;
  std::vector<int> zeroInts;
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;

  auto addAtom = [&](double x, double y, double z, int type) -> long int {
    long int id = static_cast<long int>(atomIds.size() + 1);
    atomIds.push_back(id);
    atomTypes.push_back(type);
    xPos.push_back(x);
    yPos.push_back(y);
    zPos.push_back(z);
    zeroInts.push_back(0);
    return id;
  };

  auto addStrand = [&](long int idFrom,
                       long int idTo,
                       int nBeads,
                       const Eigen::Vector3d& pFrom,
                       const Eigen::Vector3d& pTo) {
    long int prevId = idFrom;
    for (int b = 1; b <= nBeads; ++b) {
      const double frac =
        static_cast<double>(b) / static_cast<double>(nBeads + 1);
      Eigen::Vector3d p = pFrom + frac * (pTo - pFrom);
      long int bId = addAtom(p.x(), p.y(), p.z(), 1);
      if (prevId > 0) {
        bondFrom.push_back(prevId);
        bondTo.push_back(bId);
      }
      prevId = bId;
    }
    if (idTo > 0 && prevId > 0) {
      bondFrom.push_back(prevId);
      bondTo.push_back(idTo);
    }
  };

  // 4 active crosslinkers (c0..c3) forming a tetrahedron (degree 3 each)
  Eigen::Vector3d p0(5., 5., 5.);
  Eigen::Vector3d p1(15., 5., 5.);
  Eigen::Vector3d p2(5., 15., 5.);
  Eigen::Vector3d p3(5., 5., 15.);

  long int c0 = addAtom(p0.x(), p0.y(), p0.z(), 2);
  long int c1 = addAtom(p1.x(), p1.y(), p1.z(), 2);
  long int c2 = addAtom(p2.x(), p2.y(), p2.z(), 2);
  long int c3 = addAtom(p3.x(), p3.y(), p3.z(), 2);

  // 6 active strands connecting c0..c3 (each with 3 intermediate beads)
  addStrand(c0, c1, 3, p0, p1);
  addStrand(c0, c2, 3, p0, p2);
  addStrand(c0, c3, 3, p0, p3);
  addStrand(c1, c2, 3, p1, p2);
  addStrand(c1, c3, 3, p1, p3);
  addStrand(c2, c3, 3, p2, p3);

  // 1 leaf crosslinker c4 (degree 1) connected only to c1
  Eigen::Vector3d p4(25., 5., 5.);
  long int c4 = addAtom(p4.x(), p4.y(), p4.z(), 2);
  addStrand(c1, c4, 3, p1, p4);

  // 1 dangling strand connected to c2 (free end)
  Eigen::Vector3d pDangEnd(5., 25., 5.);
  addStrand(c2, -1, 3, p2, pDangEnd);

  // 1 soluble chain (both ends free)
  Eigen::Vector3d pSolStart(20., 20., 20.);
  Eigen::Vector3d pSolEnd(25., 25., 25.);
  addStrand(-1, -1, 4, pSolStart, pSolEnd);

  universe.addAtoms(
    atomIds, atomTypes, xPos, yPos, zPos, zeroInts, zeroInts, zeroInts);
  universe.addBonds(bondFrom, bondTo);

  SECTION("toCrosslinkerUniverse converts accurately")
  {
    pylimer_tools::sim::mehp::CrosslinkerUniverse cu =
      universe.toCrosslinkerUniverse(2);
    CHECK(cu.xlinkTypes.size() == 5); // c0, c1, c2, c3, c4
    CHECK(cu.strandFrom.size() >=
          8); // 6 core + 1 leaf + 1 dangling (+ soluble)
  }

  SECTION(
    "peelTopologicalLeaves and identifyActiveAtoms correctly identify 2-core")
  {
    std::vector<bool> activeAtoms = universe.identifyActiveAtoms(2);
    REQUIRE(activeAtoms.size() == universe.getNrOfAtoms());

    // c0..c3 are active
    CHECK(activeAtoms[universe.getIdxByAtomId(c0)] == true);
    CHECK(activeAtoms[universe.getIdxByAtomId(c1)] == true);
    CHECK(activeAtoms[universe.getIdxByAtomId(c2)] == true);
    CHECK(activeAtoms[universe.getIdxByAtomId(c3)] == true);

    // c4 is a leaf crosslinker (degree 1) -> must be inactive
    CHECK(activeAtoms[universe.getIdxByAtomId(c4)] == false);

    // Soluble chain atoms must be inactive
    long int solEndId = atomIds.back();
    CHECK(activeAtoms[universe.getIdxByAtomId(solEndId)] == false);
  }

  SECTION("EntanglementFilteringMode behavior in randomlyFindEntanglements")
  {
    std::vector<bool> activeAtoms = universe.identifyActiveAtoms(2);

    // PRUNE_UPFRONT: all sampled atoms must be active
    pt::AtomPairEntanglements upfront = pt::randomlyFindEntanglements(
      universe,
      10,
      8.0,
      0.0,
      0,
      1.0,
      "seed_upfront",
      2,
      true,
      pt::EntanglementFilteringMode::PRUNE_UPFRONT);

    for (const std::pair<size_t, size_t>& pair : upfront.pairsOfAtoms) {
      size_t v1 = universe.getIdxByAtomId(pair.first);
      size_t v2 = universe.getIdxByAtomId(pair.second);
      CHECK(activeAtoms[v1] == true);
      CHECK(activeAtoms[v2] == true);
    }

    CHECK(pt::parseEntanglementFilteringMode("NONE") ==
          pt::EntanglementFilteringMode::NONE);
    CHECK(pt::parseEntanglementFilteringMode("Omit_Inactive") ==
          pt::EntanglementFilteringMode::OMIT_INACTIVE);
    CHECK(pt::parseEntanglementFilteringMode("PRUNE_UPFRONT") ==
          pt::EntanglementFilteringMode::PRUNE_UPFRONT);
    CHECK_THROWS_AS(pt::parseEntanglementFilteringMode("invalid_mode"),
                    std::invalid_argument);

    // NONE: allows inactive/dangling atoms to be sampled
    pt::AtomPairEntanglements none =
      pt::randomlyFindEntanglements(universe,
                                    15,
                                    8.0,
                                    0.0,
                                    0,
                                    1.0,
                                    "seed_none",
                                    2,
                                    true,
                                    pt::EntanglementFilteringMode::NONE);
    CHECK(none.pairsOfAtoms.size() > 0);
  }
}

TEST_CASE("Entanglement Detector respects distance cut-off",
          "[analysis][EntanglementDetector][long]")
{
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::vector<std::string> files = {
    "3d-diamond-lattice_10x10x10_a_3_d_0.85_imperfect.structure.out",
    "crosslinked_M10000_N39_p_0.9.out",
    "crosslinked_p_0.98_melt_100_a_38_50_xlinks_v_22.structure.out-"
    "equilibration_do_crosslink.structure.out",
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

    for (const bool filtered : { true, false }) {
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
            INFO("Processing file: "
                 << inputFile << " with lowerCutoff: " << lowerCutoff
                 << ", upperCutoff: " << upperCutoff << ", sameStrandCutoff: "
                 << sameStrandCutoff << ", filtered: " << filtered);
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
