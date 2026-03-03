#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/io/DataFileWriter.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance2.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/utils/MCUniverseGenerator.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include <map>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pc = pylimer_tools::calc;

/**
 * Utility functions to validate the resulting distributions
 */

double
normalCDF(double x) // Phi(-∞, x) aka N(x)
{
  return std::erfc(-x / std::sqrt(2)) / 2;
}

/**
 * @brief Performs the Anderson-Darling test to check if observations follow a
 * normal distribution
 *
 * This function implements the Anderson-Darling test to determine whether a
 * given set of observations follows a normal distribution with specified mean
 * and variance. The test is performed at a 10% significance level.
 *
 * Sources:
 * - https://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm
 * - https://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test
 *
 * @param observations Vector of observed values to be tested
 * @param expectedMean The expected mean of the normal distribution
 * @param expectedVariance The expected variance of the normal distribution
 * @return bool Returns true if the observations follow the specified normal
 * distribution (null hypothesis cannot be rejected at 10% significance level),
 *              false otherwise
 */
bool
andersonDarlingNormalDistributionTest(Eigen::VectorXd observations,
                                      const double expectedMean,
                                      const double expectedVariance)
{
  std::ranges::sort(observations);
  const double sampleMean = observations.mean();
  const double sampleVariance =
    (observations.array() - sampleMean).square().sum() /
    (observations.size() - 1);
  const double n = observations.size();

  // if (std::abs(sampleMean) < 1e-4 || std::abs(expectedMean) < 1e-4) {
  //   CHECK_THAT(sampleMean, Catch::Matchers::WithinAbs(expectedMean, 10. /
  //   n));
  // } else {
  //   CHECK_THAT(sampleMean, Catch::Matchers::WithinRel(expectedMean, 0.1));
  // }
  // CHECK_THAT(sampleVariance, Catch::Matchers::WithinRel(expectedVariance,
  // 0.1)); Eigen::VectorXd normalizedObservations =
  //   (observations.array() - sampleMean).matrix() / std::sqrt(sampleVariance);
  Eigen::VectorXd normalizedObservations =
    (observations.array() - expectedMean).matrix() /
    std::sqrt(expectedVariance);

  double A2 = -n;
  for (double i = 0; i < n; ++i) {
    A2 -= ((2. * (i + 1) - 1.) / n) *
          (std::log(
             normalCDF(normalizedObservations(static_cast<Eigen::Index>(i)))) +
           std::log(1. - normalCDF(normalizedObservations(
                           static_cast<Eigen::Index>(n - (i + 1))))));
  }

  // Case 0, 10 % significance level,
  // according to Marsaglia & Marsaglia
  return A2 < 1.933;
}

/**
 * Actual test cases
 */
TEST_CASE("Certain configurations do not lead to memory corruption",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"Certain configurations do not lead to memory "
               "corruption\""
            << std::endl;
  // the following parameters have led to a `double free or corruption` error?!?
  constexpr int nrOfCrosslinkers = static_cast<int>(5e4 * 2 * 0.7 / 7);
  const double sideLength = std::cbrt((10 * 5e4 * nrOfCrosslinkers) / 0.85);
  pu::MCUniverseGenerator generator =
    pu::MCUniverseGenerator(sideLength, sideLength, sideLength);
  REQUIRE_NOTHROW(generator.setSeed(68419));
  REQUIRE_NOTHROW(generator.setBeadDistance(0.965));

  pe::Universe universe = generator.getUniverse();
  REQUIRE(universe.getNrOfAtoms() == 0);

  generator.addCrosslinkers(nrOfCrosslinkers);
  REQUIRE_THROWS(generator.linkStrandsToConversion(0.2));
  universe = generator.getUniverse();
  REQUIRE(universe.getNrOfAtoms() == nrOfCrosslinkers);
}

TEST_CASE("Universe can be generated", "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"Universe can be generated\"" << std::endl;
  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(100, 4, 2);
  generator.addSolventChains(100, 16, 3);
  generator.addStrands((4 / 2) * 100, 16);
  generator.configNrOfMCSteps(10);
  generator.useLinearMaxDistance(1.);
  generator.linkStrandsToConversion(0.8);

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
    generator2.addCrosslinkers(100, 4, 2);
    generator2.addSolventChains(100, 16, 3);
    generator2.addStrands((4 / 2) * 100, 16);
    generator2.configNrOfMCSteps(10);
    generator.useLinearMaxDistance(1.);
    generator2.linkStrandsToConversion(0.8);

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
      CHECK(bondLength < 3.5);
    }
  }

  SECTION("Errors are thrown")
  {
    // nr of strands and strand lengths must be same:
    REQUIRE_THROWS(generator.addStrands(3, { { 10, 100 } }, 1));
    // not enough strands to reach conversion:
    generator.addStrands(2, 10, 1);
    REQUIRE_THROWS(generator.linkStrandsToConversion(2.0));
    REQUIRE_THROWS(generator.addRandomlyFunctionalizedStrands(
      2, { 10, 10 }, 7, 3, 2, 2, true));
    CHECK_THROWS(generator.addRandomlyFunctionalizedStrands(0, {}, 2));
    CHECK_THROWS(generator.addRandomlyFunctionalizedStrands(1, {}, 2));
    CHECK_THROWS(generator.addMonofunctionalStrands(1, {}, -1));
    CHECK_THROWS(generator.addCrosslinkersAt(Eigen::VectorXd::Zero(2)));
    CHECK_THROWS(generator.addCrosslinkersAt(Eigen::VectorXd::Zero(6), -1));
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

  const pe::Universe universe = generator.getUniverse();
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
  CHECK_THROWS(generator.setMeanSquaredBeadDistance(
    std::numeric_limits<double>::infinity()));
  CHECK_THROWS(
    generator.setBeadDistance(std::numeric_limits<double>::infinity()));
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
  generator.configSecondaryLoopProbability(0.5);
  generator.addStrands(200, 10, 1);
  CHECK(generator.getCurrentCrosslinkerConversion() == 0.0);
  CHECK(generator.getCurrentStrandsConversion() == 0.0);

  SECTION("Without max distance")
  {
    generator.linkStrandsToConversion(0.9, 1.);
    CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
               Catch::Matchers::WithinAbs(0.9, 1e-3));
    CHECK_THAT(generator.getCurrentStrandsConversion(),
               Catch::Matchers::WithinAbs(0.9, 1e-3));

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

  SECTION("With max distance")
  {
    generator.useLinearMaxDistance(3.);
    generator.linkStrandsToConversion(0.9, 1.);

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
}

TEST_CASE("Universe can crosslink up to w_sol",
          "[generator][MCUniverseGenerator][MEHPForceRelaxation]["
          "MEHPForceBalance][MEHPForceBalance2][long]")
{
  std::cout << "Running test \"Universe can crosslink up to w_sol\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(100, 4, 2);

  generator.addStrands(200, 19, 1);
  generator.linkStrandsToSolubleFraction(0.1, 1.);

  pe::Universe universe = generator.getUniverse();
  std::map<int, double> weights;
  weights[1] = 1.0;
  weights[2] = 1.0;
  universe.setMasses(weights);
  REQUIRE(universe.getVolume() == 10.0 * 10.0 * 10.0);
  REQUIRE(universe.getAtomsOfType(2).size() == 100);
  REQUIRE(universe.getAtomsOfType(1).size() == 200 * 19);

  size_t nFreeStrands = 0;
  for (const pe::Molecule& mol : universe.getChainsWithCrosslinker(2)) {
    nFreeStrands += mol.getType() == pe::MoleculeType::FREE_CHAIN;
  }

  auto clusters = universe.getClusters();
  CHECK(clusters.size() < 20);

  pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
    pylimer_tools::sim::mehp::MEHPForceRelaxation(universe);

  while (forceRelaxer.suggestsRerun()) {
    forceRelaxer.runForceRelaxation();
  }

  CHECK_THAT(
    0.1,
    Catch::Matchers::WithinAbs(forceRelaxer.getSolubleWeightFraction(), 0.05));

  pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer2 =
    generator.getForceRelaxation();

  pylimer_tools::sim::mehp::MEHPForceBalance forceBalance =
    generator.getForceBalance();
  forceBalance.configAssumeBoxLargeEnough(true);

  // compare structures before running optimization procedures
  CHECK(forceBalance.getNrOfSprings() == forceRelaxer2.getNrOfSprings());
  CHECK(forceBalance.getNrOfSprings() ==
        forceRelaxer.getNrOfSprings() + nFreeStrands);
  Eigen::VectorXd fbVecs = forceBalance.getCurrentPartialSpringDistances();
  Eigen::VectorXd frVecs = forceRelaxer2.getCurrentSpringDistances();
  REQUIRE(fbVecs.size() == 3 * forceBalance.getNrOfSprings());
  CHECK(frVecs.isApprox(fbVecs));
  CHECK_THAT(
    forceBalance.getNrOfActiveSprings(1e-2),
    Catch::Matchers::WithinAbs(forceRelaxer2.getNrOfActiveSprings(1e-2), 1));
  CHECK_THAT(
    forceBalance.getGammaFactor(1.),
    Catch::Matchers::WithinRel(forceRelaxer2.getGammaFactor(1.), 0.01));

  // run force relaxation
  while (forceRelaxer2.suggestsRerun()) {
    forceRelaxer2.runForceRelaxation();
  }

  CHECK_THAT(
    forceRelaxer2.countActiveClusteredAtoms(),
    Catch::Matchers::WithinAbs(forceRelaxer.countActiveClusteredAtoms(), 1));

  // also run force balance
  forceBalance.runForceRelaxation();

  CHECK(forceBalance.validateNetwork());
  CHECK(forceBalance.getNrOfSprings() ==
        forceRelaxer.getNrOfSprings() + nFreeStrands);
  CHECK_THAT(forceBalance.getGammaFactor(1.),
             Catch::Matchers::WithinRel(forceRelaxer.getGammaFactor(1.), 0.05));
  CHECK_THAT(
    forceBalance.getNrOfActiveSprings(1e-2),
    Catch::Matchers::WithinAbs(forceRelaxer.getNrOfActiveSprings(1e-2), 1));

  // and finally force balance 2
  pylimer_tools::sim::mehp::MEHPForceBalance2 forceBalance2 =
    pylimer_tools::sim::mehp::MEHPForceBalance2(universe);
  forceBalance2.runForceRelaxation();
  CHECK_THAT(forceBalance2.getSolubleWeightFraction(),
             Catch::Matchers::WithinAbs(0.1, 0.05));
  double activeClusteredAtomsFromUniverse =
    forceBalance2.countActiveClusteredAtoms();
  forceBalance2.configAssumeNetworkIsComplete(true);
  CHECK_THAT(forceBalance2.getSolubleWeightFraction(),
             Catch::Matchers::WithinAbs(0.1, 0.05));
  double activeClusteredAtomsFromNetwork =
    forceBalance2.countActiveClusteredAtoms();
  CHECK_THAT(activeClusteredAtomsFromUniverse,
             Catch::Matchers::WithinRel(activeClusteredAtomsFromNetwork, 0.01));

  pylimer_tools::sim::mehp::MEHPForceBalance2 forceBalance2Generator =
    generator.getForceBalance2();
  forceBalance2Generator.runForceRelaxation();
  CHECK_THAT(forceBalance2Generator.getSolubleWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalance2.getSolubleWeightFraction(), 0.005));
  double activeClusteredAtomsFromGenerator =
    forceBalance2Generator.countActiveClusteredAtoms();
  CHECK_THAT(
    activeClusteredAtomsFromNetwork,
    Catch::Matchers::WithinRel(activeClusteredAtomsFromGenerator, 0.01));

  Eigen::ArrayXi clusterIndicesPerLink =
    Eigen::ArrayXi::Constant(forceBalance2Generator.getNrOfLinks(), -1);
  Eigen::ArrayXi clusterIndicesPerStrand =
    Eigen::ArrayXi::Constant(forceBalance2Generator.getNrOfStrands(), -1);
  size_t minIndexWithoutCluster = 0;
  size_t currentClusterIdx = 0;
  pylimer_tools::sim::mehp::ForceBalance2Network net2 =
    forceBalance2Generator.getNetwork();
  while ((clusterIndicesPerLink < 0).any()) {
    while (clusterIndicesPerLink[minIndexWithoutCluster] >= 0) {
      minIndexWithoutCluster++;
    }
    // find clusters
    clusterIndicesPerLink[minIndexWithoutCluster] = currentClusterIdx;
    bool didChange = true;
    while (didChange) {
      didChange = false;
      for (Eigen::Index linkIdx = 0; linkIdx < net2.nrOfLinks; ++linkIdx) {
        if (clusterIndicesPerLink[linkIdx] != currentClusterIdx) {
          continue;
        }

        for (const int strandIdx : net2.strandIndicesOfLink[linkIdx]) {
          clusterIndicesPerStrand[strandIdx] = currentClusterIdx;
          for (const int strandsLinkIdx : net2.linkIndicesOfStrand[strandIdx]) {
            if (clusterIndicesPerLink[strandsLinkIdx] < 0) {
              didChange = true;
              clusterIndicesPerLink[strandsLinkIdx] = currentClusterIdx;
            }
          }
        }
      }
    }
    currentClusterIdx += 1;
  }

  CHECK(clusterIndicesPerLink.maxCoeff() + 1 == universe.getClusters().size());
  CHECK(clusterIndicesPerLink.maxCoeff() + 1 == currentClusterIdx);
  // output cluster sizes
  std::map<int, int> clusterSizes;
  for (const int clusterIdx : clusterIndicesPerLink) {
    if (!clusterSizes.contains(clusterIdx)) {
      clusterSizes[clusterIdx] = 0;
    }
    clusterSizes[clusterIdx] += 1;
  }
  for (int strandIdx = 0; strandIdx < net2.nrOfStrands; ++strandIdx) {
    double strandNAtoms = 0;
    for (size_t springIdx : net2.springIndicesOfStrand[strandIdx]) {
      strandNAtoms += net2.springContourLength[springIdx] - 1;
    }
    clusterSizes[clusterIndicesPerStrand[strandIdx]] += strandNAtoms;
  }

  std::cout << "Cluster sizes: ";
  for (const auto& [clusterIdx, size] : clusterSizes) {
    std::cout << clusterIdx << ": " << size << ", ";
  }
  std::cout << std::endl;
}

TEST_CASE("Generator can return force balance and relaxation",
          "[generator][MCUniverseGenerator][long]")
{
  std::cout
    << "Running test \"Generator can return force balance and relaxation\""
    << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.addCrosslinkers(100, 4, 2);

  generator.addStrands(200, 19, 1);
  generator.linkStrandsToConversion(0.925, 1.);

  pe::Universe universe = generator.getUniverse();
  REQUIRE(universe.getVolume() == 10.0 * 10.0 * 10.0);
  REQUIRE(universe.getAtomsOfType(2).size() == 100);
  REQUIRE(universe.getAtomsOfType(1).size() == 200 * 19);

  pylimer_tools::sim::mehp::ForceBalanceNetwork fbNet =
    generator.convertToForceBalanceNetwork();
  pylimer_tools::sim::mehp::Network frNet =
    generator.convertToForceRelaxationNetwork();

  REQUIRE(fbNet.nrOfSprings == frNet.nrOfSprings);
  REQUIRE(fbNet.coordinates.isApprox(frNet.coordinates));
  REQUIRE(fbNet.springIndexA.isApprox(frNet.springIndexA));
  REQUIRE(fbNet.springIndexB.isApprox(frNet.springIndexB));
  REQUIRE(fbNet.springPartIndexA.isApprox(frNet.springIndexA));
  REQUIRE(fbNet.springPartIndexB.isApprox(frNet.springIndexB));
  REQUIRE(
    fbNet.springPartCoordinateIndexA.isApprox(frNet.springCoordinateIndexA));
  REQUIRE(
    fbNet.springPartCoordinateIndexB.isApprox(frNet.springCoordinateIndexB));
  REQUIRE(fbNet.springPartBoxOffset.isApprox(frNet.springBoxOffset));
  REQUIRE(fbNet.vol == frNet.vol);
  for (size_t dir = 0; dir < 3; dir++) {
    REQUIRE(fbNet.L[0] == frNet.L[0]);
  }

  pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer2 =
    generator.getForceRelaxation();
  forceRelaxer2.configAssumeBoxLargeEnough(true);

  pylimer_tools::sim::mehp::MEHPForceBalance forceBalance =
    generator.getForceBalance();
  forceBalance.configAssumeBoxLargeEnough(true);
  // compare structures before running optimization procedures
  CHECK(forceBalance.getNrOfSprings() == forceRelaxer2.getNrOfSprings());
  Eigen::VectorXd fbVecs = forceBalance.getCurrentPartialSpringDistances();
  Eigen::VectorXd fbVecs2 = forceBalance.getCurrentSpringDistances();
  REQUIRE(fbVecs.isApprox(fbVecs2));
  Eigen::VectorXd frVecs = forceRelaxer2.getCurrentSpringDistances();
  Eigen::VectorXd frVecs2 = forceRelaxer2.getSpringDistances();
  REQUIRE(frVecs.isApprox(frVecs2));
  REQUIRE(fbVecs.size() == 3 * forceBalance.getNrOfSprings());
  REQUIRE(frVecs.size() == 3 * forceBalance.getNrOfSprings());
  CHECK(frVecs.isApprox(fbVecs));
  CHECK_THAT(
    forceBalance.getNrOfActiveSprings(1e-2),
    Catch::Matchers::WithinAbs(forceRelaxer2.getNrOfActiveSprings(1e-2), 1));
  CHECK_THAT(
    forceBalance.getGammaFactor(1.),
    Catch::Matchers::WithinRel(forceRelaxer2.getGammaFactor(1.), 0.01));

  // run force relaxation
  while (forceRelaxer2.suggestsRerun()) {
    forceRelaxer2.runForceRelaxation();
  }

  // also run force balance
  forceBalance.runForceRelaxation();

  CHECK(forceBalance.validateNetwork());
  CHECK(forceBalance.getNrOfSprings() == forceRelaxer2.getNrOfSprings());
  CHECK_THAT(
    forceBalance.getGammaFactor(1.),
    Catch::Matchers::WithinRel(forceRelaxer2.getGammaFactor(1.), 0.025));
  CHECK_THAT(
    forceBalance.getNrOfActiveSprings(1e-2),
    Catch::Matchers::WithinAbs(forceRelaxer2.getNrOfActiveSprings(1e-2), 1));
}

TEST_CASE("MCUniverseGenerator can remove w_sol",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"MCUniverseGenerator can remove w_sol\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.75);
  generator.addCrosslinkers(400, 4, 2);
  generator.configNrOfMCSteps(0);

  generator.addStrands(800, 10, 1);
  generator.linkStrandsToConversion(0.85, 1.);

  pe::Universe universeBeforeRemoval = generator.getUniverse();
  CHECK(universeBeforeRemoval.getAtomsOfType(2).size() == 400);
  CHECK(universeBeforeRemoval.getAtomsOfType(1).size() == 800 * 10);

  SECTION("Without rescaling")
  {
    generator.removeSolubleFraction(false);

    pe::Universe universeAfterRemoval = generator.getUniverse();
    CHECK(universeAfterRemoval.getVolume() ==
          universeBeforeRemoval.getVolume());

    CHECK(universeAfterRemoval.getAtomsOfType(2).size() < 400);
    CHECK(universeAfterRemoval.getAtomsOfType(1).size() < 800 * 10);
  }

  SECTION("With rescaling")
  {
    generator.removeSolubleFraction(true);

    pe::Universe universeAfterRemoval = generator.getUniverse();
    CHECK(universeAfterRemoval.getVolume() < universeBeforeRemoval.getVolume());
    double densityAfter =
      universeAfterRemoval.getNrOfAtoms() / universeAfterRemoval.getVolume();
    double densityBefore =
      universeBeforeRemoval.getNrOfAtoms() / universeBeforeRemoval.getVolume();
    CHECK(densityAfter == Catch::Approx(densityBefore));

    CHECK(universeAfterRemoval.getAtomsOfType(2).size() < 400);
    CHECK(universeAfterRemoval.getAtomsOfType(1).size() < 800 * 10);
  }
}

TEST_CASE("MCUniverseGenerator uses correct force relaxation network",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"MCUniverseGenerator uses correct force "
               "relaxation network\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.75);
  generator.addCrosslinkers(400, 4, 2);
  generator.configNrOfMCSteps(0);

  generator.addStrands(800, 10, 1);
  generator.useZScoreMaxDistance(3.29, 1.);
  generator.linkStrandsToConversion(0.925, 1.);

  pe::Universe universe = generator.getUniverse();
  CHECK(universe.getNrOfAtoms() == 400 + 800 * 10);
  CHECK(generator.getCurrentNrOfAtoms() == universe.getNrOfAtoms());
  CHECK(generator.getCurrentNrOfBonds() == universe.getNrOfBonds());

  size_t nFreeStrands = 0;
  for (const pe::Molecule& mol : universe.getChainsWithCrosslinker(2)) {
    nFreeStrands += mol.getType() == pe::MoleculeType::FREE_CHAIN;
  }

  pylimer_tools::sim::mehp::MEHPForceRelaxation relaxer =
    pylimer_tools::sim::mehp::MEHPForceRelaxation(
      universe, 2, false, nullptr, 1.0, false, false);
  relaxer.configAssumeBoxLargeEnough(true);

  pylimer_tools::sim::mehp::MEHPForceRelaxation relaxerFromGenerator =
    pylimer_tools::sim::mehp::MEHPForceRelaxation(
      generator.convertToForceRelaxationNetwork());
  relaxerFromGenerator.configAssumeBoxLargeEnough(true);

  // comparison only after running force relaxation,
  // since the box offsets, single nodes and stuff is relevant before that
  // CHECK(relaxer.countActiveClusteredAtoms() ==
  //       Catch::Approx(relaxerFromGenerator.countActiveClusteredAtoms()));
  // CHECK(relaxer.getNrOfNodes() == relaxerFromGenerator.getNrOfNodes());
  CHECK(relaxer.getNrOfSprings() + nFreeStrands ==
        relaxerFromGenerator.getNrOfSprings());

  relaxer.runForceRelaxation();
  relaxerFromGenerator.runForceRelaxation();

  CHECK(relaxer.countActiveClusteredAtoms() ==
        Catch::Approx(relaxerFromGenerator.countActiveClusteredAtoms()));
  // cannot compare, since each strand end becomes a node in the from-universe
  // generation
  // CHECK(relaxer.getNrOfNodes() == relaxerFromGenerator.getNrOfNodes());
  CHECK(relaxer.getNrOfSprings() + nFreeStrands ==
        relaxerFromGenerator.getNrOfSprings());
}

TEST_CASE("MUniverseGenerator can generate with crosslink chains",
          "[generator][MCUniverseGenerator]")
{
  std::cout
    << "Running test \"MUniverseGenerator can generate with crosslink chains\""
    << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.75);
  generator.configNrOfMCSteps(0);

  generator.addCrosslinkStrands(
    10, { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 }, 2, 2, 1);

  pe::Universe universe = generator.getUniverse();
  CHECK(generator.getCurrentNrOfAtoms() == universe.getNrOfAtoms());
  CHECK(generator.getCurrentNrOfBonds() == universe.getNrOfBonds());
  CHECK(universe.getAtomsOfType(2).size() == 2 * 10);
  CHECK(universe.getAtomsOfType(1).size() == 10 * 10);

  std::vector<pe::Molecule> chains = universe.getChainsWithCrosslinker(2);
  CHECK(chains.size() == 10);
  for (pe::Molecule& chain : chains) {
    CHECK(chain.getNrOfAtoms() == 12);
  }

  generator.configPrimaryLoopProbability(0.);
  generator.configSecondaryLoopProbability(0.);
  generator.addStrands(
    20, pylimer_tools::utils::initializeWithValue(20, 10), 3);
  generator.linkStrandsToConversion(0.5);
  universe = generator.getUniverse();

  CHECK(universe.getAtomsOfType(2).size() == 2 * 10);
  CHECK(universe.getAtomsOfType(1).size() == 10 * 10);
  CHECK(universe.getAtomsOfType(3).size() == 20 * 10);
}

TEST_CASE(
  "MCUniverseGenerator can generate with randomly functionalized chains",
  "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"MCUniverseGenerator can generate with randomly "
               "functionalized chains\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.75);
  generator.configNrOfMCSteps(0);

  SECTION("Without functionalization")
  {
    generator.addRandomlyFunctionalizedStrands(
      10, { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 }, 0., 2, 1, true);

    pe::Universe universe = generator.getUniverse();
    CHECK(universe.getAtomsOfType(2).size() == 0);
    CHECK(universe.getAtomsOfType(1).size() == 10 * 10);

    std::vector<pe::Molecule> chains = universe.getChainsWithCrosslinker(2);
    CHECK(chains.size() == 10);
    for (pe::Molecule& chain : chains) {
      CHECK(chain.getNrOfAtoms() == 10);
    }
  }

  SECTION("With functionalization")
  {
    std::vector<int> chainLengths = pu::initializeWithValue(1000, 100);
    generator.addRandomlyFunctionalizedStrands(
      1000, chainLengths, 0.2, 3, 2, 1, true);

    pe::Universe universe = generator.getUniverse();
    double nCrosslinks = static_cast<double>(universe.getAtomsOfType(2).size());
    CHECK_THAT(nCrosslinks, Catch::Matchers::WithinRel(1000 * 100 * 0.2, 0.01));
    CHECK(universe.getNrOfAtoms() == 1000 * 100);
    std::vector<pe::Molecule> chains = universe.getChainsWithCrosslinker(2);
    CHECK(chains.size() > 5000);
  }

  SECTION("With functionalization probability > 1")
  {
    std::vector<int> chainLengths = pu::initializeWithValue(1000, 100);
    generator.addRandomlyFunctionalizedStrands(
      1000, chainLengths, 3.2, 1, 2, 1, true);
    CHECK_THAT(generator.getCurrentNrOfAvailableCrosslinkSites(),
               Catch::Matchers::WithinRel(1000 * 100 * 3.2, 0.01));

    pe::Universe universe = generator.getUniverse();
    double nCrosslinks = static_cast<double>(universe.getAtomsOfType(2).size());
    CHECK_THAT(nCrosslinks, Catch::Matchers::WithinRel(1000 * 100, 0.05));
    CHECK(universe.getNrOfAtoms() == 1000 * 100);
    std::vector<pe::Molecule> chains = universe.getChainsWithCrosslinker(2);
    CHECK(chains.size() > 5000);
    chains = universe.getChainsWithCrosslinker(6);
    CHECK(chains.size() == 1000);
  }
}

TEST_CASE("MCUniverseGenerator can generate with monofunctional chains",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"MCUniverseGenerator can generate with "
               "monofunctional chains\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(457564875e2);
  generator.setBeadDistance(0.75);
  generator.configNrOfMCSteps(0);

  generator.addMonofunctionalStrands(100, 100, 5);
  CHECK(generator.getCurrentNrOfAtoms() == 100 * 100);
  generator.addCrosslinkers(25, 4, 2, false);
  CHECK(generator.getCurrentCrosslinkerConversion() == 0.);
  CHECK(generator.getCurrentNrOfAtoms() == 100 * 100 + 25);
  CHECK_NOTHROW(generator.linkStrandsToConversion(1.));
  CHECK_THROWS(generator.linkStrandsToConversion(0.));
  pe::Universe universe = generator.getUniverse();
  std::vector<pe::MoleculeType> moleculeTypes =
    universe.identifyObviouslyDanglingAtoms(true);
  CHECK(std::ranges::all_of(moleculeTypes, [](pe::MoleculeType type) {
    return type == pe::MoleculeType::FREE_CHAIN;
  }));

  generator.relaxCrosslinks();

  pylimer_tools::sim::mehp::MEHPForceBalance2 fb2 =
    generator.getForceBalance2();
  fb2.runForceRelaxation();
  CHECK(fb2.getNrOfActiveSprings() == 0);
  CHECK(fb2.getNrOfActiveStrands() == 0);

  pylimer_tools::sim::mehp::MEHPForceBalance fb = generator.getForceBalance();
  fb.runForceRelaxation();
  CHECK(fb.getNrOfActiveSprings() == 0);
}

TEST_CASE("Randomly functionalized chains collapse",
          "[generator][MCUniverseGenerator][long]")
{
  std::cout << "Running test \"Randomly functionalized chains collapse\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(457564875e2);
  generator.setBeadDistance(0.75);
  generator.configNrOfMCSteps(0);

  std::vector<int> chainLengths = pu::initializeWithValue(20, 100);
  generator.addRandomlyFunctionalizedStrands(
    20, chainLengths, 7.2, 1, 2, 1, true);

  CHECK_THAT(generator.getCurrentNrOfAvailableCrosslinkSites(),
             Catch::Matchers::WithinRel(20 * 100 * 7.2, 0.01));
  CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
             Catch::Matchers::WithinAbs(0.0, 1e-10));

  pe::Universe universe = generator.getUniverse();
  double nCrosslinks = static_cast<double>(universe.getAtomsOfType(2).size());
  CHECK_THAT(nCrosslinks, Catch::Matchers::WithinRel(20 * 100, 0.01));
  CHECK(universe.getNrOfAtoms() == 20 * 100);

  std::vector<pe::Molecule> chains = universe.getMolecules(6);
  CHECK(chains.size() == 20);
  for (pe::Molecule& chain : chains) {
    CHECK(chain.getNrOfAtoms() == 100);
    CHECK_THAT(chain.getAtomsOfType(2).size(),
               Catch::Matchers::WithinAbs(100., 2));
  }

  // first with "incorrect" spring type
  pylimer_tools::sim::mehp::MEHPForceBalance forceBalanceEasy =
    pylimer_tools::sim::mehp::MEHPForceBalance(universe, 6);
  CHECK(forceBalanceEasy.getNrOfSprings() == 0);

  forceBalanceEasy.runForceRelaxation();

  CHECK_THAT(forceBalanceEasy.getSolubleWeightFraction(),
             Catch::Matchers::WithinAbs(1.0, 0.001));

  // then, with many short springs
  pylimer_tools::sim::mehp::MEHPForceBalance forceBalance =
    pylimer_tools::sim::mehp::MEHPForceBalance(universe);

  pylimer_tools::sim::OutputConfiguration outputConfig;
  outputConfig.intValues = { pylimer_tools::sim::ComputedIntValues::STEP };
  outputConfig.doubleValues = {
    pylimer_tools::sim::ComputedDoubleValues::RESIDUAL,
    pylimer_tools::sim::ComputedDoubleValues::GAMMA
  };
  outputConfig.outputEvery = 100;
  std::vector<pylimer_tools::sim::OutputConfiguration> outputConfigs = {
    outputConfig
  };
  forceBalance.configStepOutput(outputConfigs);

  SECTION("With big box")
  {
    forceBalance.configAssumeBoxLargeEnough(true);
    CHECK(forceBalance.getResidual() > 0.);
    forceBalance.runForceRelaxation(50000, 1e-12);

    CHECK_THAT(forceBalance.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK(forceBalance.getNrOfNodes() <= 20 * 100);
  }

  SECTION("With small box")
  {
    forceBalance.configAssumeBoxLargeEnough(false);
    CHECK(forceBalance.getResidual() > 0.);
    forceBalance.runForceRelaxation(50000, 1e-12);

    CHECK_THAT(forceBalance.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK(forceBalance.getNrOfNodes() <= 20 * 100);
  }
}

TEST_CASE(
  "Single randomly functionalized chain can be prevented from building loops",
  "[generator][MCUniverseGenerator][long]")
{

  std::cout << "Running test \"Single randomly functionalized chain can be "
               "prevented from building loops\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.75);
  generator.configNrOfMCSteps(0);

  SECTION("No primary loops")
  {
    generator.addRandomlyFunctionalizedStrands(1, { 100 }, 0.8, 4, 2, 1, true);
    double nCrosslinkSites =
      static_cast<double>(generator.getCurrentNrOfAvailableCrosslinkSites());
    CHECK_THAT(nCrosslinkSites,
               Catch::Matchers::WithinRel(100 * 0.8 * 4, 0.05));
    CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
               Catch::Matchers::WithinAbs(0.0, 1e-10));

    generator.addStrands(5, { 10, 10, 10, 10, 10 }, 1);
    generator.configPrimaryLoopProbability(0.0);
    CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
               Catch::Matchers::WithinAbs(0.0, 1e-10));
    SECTION("To conversion")
    {
      generator.linkStrandsToConversion(9. / nCrosslinkSites);
      // we get all five dangling strands, not more
      CHECK(generator.getCurrentNrOfAvailableCrosslinkSites() ==
            nCrosslinkSites - 5);
      CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
                 Catch::Matchers::WithinRel(5. / nCrosslinkSites, 0.05));

      std::vector<pe::Molecule> chains =
        generator.getUniverse().getChainsWithCrosslinker(6);
      CHECK(chains.size() == 11);
      for (pe::Molecule& chain : chains) {
        CHECK(chain.getType() != pe::MoleculeType::PRIMARY_LOOP);
      }
    }

    SECTION("To soluble fraction")
    {
      generator.linkStrandsToSolubleFraction(0.);
      // we get all five dangling strands, not more
      CHECK(generator.getCurrentNrOfAvailableCrosslinkSites() ==
            nCrosslinkSites - 5);
      CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
                 Catch::Matchers::WithinRel(5. / nCrosslinkSites, 0.05));

      std::vector<pe::Molecule> chains =
        generator.getUniverse().getChainsWithCrosslinker(6);
      CHECK(chains.size() == 11);
      for (pe::Molecule& chain : chains) {
        CHECK(chain.getType() != pe::MoleculeType::PRIMARY_LOOP);
      }
    }
  }

  SECTION("No primary nor secondary loops")
  {
    generator.addRandomlyFunctionalizedStrands(
      2, { 100, 100 }, 0.8, 4, 2, 1, true);
    double nCrosslinkSites =
      static_cast<double>(generator.getCurrentNrOfAvailableCrosslinkSites());
    CHECK_THAT(nCrosslinkSites,
               Catch::Matchers::WithinRel(100 * 0.8 * 4 * 2, 0.05));
    CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
               Catch::Matchers::WithinAbs(0.0, 1e-10));

    // generator.addStrands(2, { 10, 10 }, 1);
    generator.addStrands(5, { 10, 10, 10, 10, 10 }, 1);
    generator.configPrimaryLoopProbability(0.0);
    generator.configSecondaryLoopProbability(0.0);
    CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
               Catch::Matchers::WithinAbs(0.0, 1e-10));

    SECTION("To soluble fraction")
    {
      generator.linkStrandsToSolubleFraction(0.);
      // we get all five dangling strands and one connecting strand,
      // afterwards, things would be primary/secondary loops
      CHECK(generator.getCurrentNrOfAvailableCrosslinkSites() ==
            nCrosslinkSites - 6);
      CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
                 Catch::Matchers::WithinRel(6. / (nCrosslinkSites), 0.05));

      std::vector<pe::Molecule> chains =
        generator.getUniverse().getChainsWithCrosslinker(6);
      for (pe::Molecule& chain : chains) {
        CHECK(chain.getType() != pe::MoleculeType::PRIMARY_LOOP);
      }
    }

    SECTION("To conversion")
    {
      generator.linkStrandsToConversion(9.1 / nCrosslinkSites);
      // we get all five dangling strands and one connecting strand,
      // afterwards, things would be primary/secondary loops
      CHECK(generator.getCurrentNrOfAvailableCrosslinkSites() ==
            nCrosslinkSites - 6);
      CHECK_THAT(generator.getCurrentCrosslinkerConversion(),
                 Catch::Matchers::WithinRel(6. / (nCrosslinkSites), 0.05));

      std::vector<pe::Molecule> chains =
        generator.getUniverse().getChainsWithCrosslinker(6);
      for (pe::Molecule& chain : chains) {
        CHECK(chain.getType() != pe::MoleculeType::PRIMARY_LOOP);
      }
    }
  }
}

TEST_CASE("Regularly spaced functionalized strands can be generated",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"Regularly spaced functionalized strands can be "
               "generated\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(8804);
  generator.setBeadDistance(0.964);
  generator.configNrOfMCSteps(10);

  SECTION("Basic functionality")
  {
    // Add strands with crosslinks every 5 beads, starting at offset 2
    std::vector<int> strandLengths = { 20, 25, 30 };
    generator.addRegularlySpacedFunctionalizedStrands(
      3, strandLengths, 5, 2, 4, 2, 1, true);

    pe::Universe universe = generator.getUniverse();

    // Check that we have the correct number of crosslinks
    // For strand 0 (20 beads):
    // crosslinks at positions 2, 7, 12, 17 = 4 crosslinks, 5 chains
    // For strand 1 (25 beads):
    // crosslinks at positions 2, 7, 12, 17, 22 = 5 crosslinks, 6 chains
    // For strand 2 (30 beads):
    // crosslinks at positions 2, 7, 12, 17, 22, 27 = 6 crosslinks, 7 chains
    // Total: 4 + 5 + 6 = 15 crosslinks
    CHECK(universe.getAtomsOfType(2).size() == 15);

    CHECK(universe.getMolecules(3).size() == 3); // only three long chains
    CHECK(universe.getMolecules(2).size() ==
          5 + 6 + 7); // but more short chains

    // Total atoms should be 20 + 25 + 30 = 75 (each bead is either crosslink or
    // strand atom)
    CHECK(universe.getNrOfAtoms() == 75);

    CHECK(generator.getCurrentNrOfAvailableCrosslinkSites() == 4 * 15);
  }

  SECTION("Zero offset")
  {
    // Test with zero offset - crosslinks at positions 0, 3, 6, 9, ...
    std::vector<int> strandLengths = { 10 };
    generator.addRegularlySpacedFunctionalizedStrands(
      1, strandLengths, 3, 0, 5, 2, 1, true);

    pe::Universe universe = generator.getUniverse();

    // For a 10-bead strand with spacing 3 and offset 0: crosslinks at 0, 3, 6,
    // 9 = 4 crosslinks
    CHECK(universe.getAtomsOfType(2).size() == 4);
    CHECK(universe.getNrOfAtoms() == 10);
    CHECK(universe.getMolecules(3).size() == 1);   // only one long chain
    REQUIRE(universe.getMolecules(2).size() == 3); // three short chains

    std::vector<pe::Molecule> chains = universe.getChainsWithCrosslinker(2);
    CHECK(chains[0].getNrOfAtoms() == 4);
    CHECK(chains[1].getNrOfAtoms() == 4);
    CHECK(chains[2].getNrOfAtoms() == 4);

    CHECK(generator.getCurrentNrOfAvailableCrosslinkSites() == 5 * 4);
  }

  SECTION("Large spacing")
  {
    // Test with large spacing - should have few crosslinks
    std::vector<int> strandLengths = { 10 };
    generator.addRegularlySpacedFunctionalizedStrands(
      1, strandLengths, 20, 0, 4, 2, 1, true);

    pe::Universe universe = generator.getUniverse();

    // For a 10-bead strand with spacing 20: only crosslink at position 0 = 1
    // crosslink
    CHECK(universe.getAtomsOfType(2).size() == 1);
    CHECK(universe.getNrOfAtoms() == 10);
  }

  SECTION("Large offset")
  {
    // Test with offset beyond strand length - should have no crosslinks
    std::vector<int> strandLengths = { 10 };
    generator.addRegularlySpacedFunctionalizedStrands(
      1, strandLengths, 3, 15, 4, 2, 1, true);

    pe::Universe universe = generator.getUniverse();

    // No crosslinks should be created since offset is beyond strand length
    CHECK(universe.getAtomsOfType(2).size() == 0);
    CHECK(universe.getNrOfAtoms() == 10);
  }

  SECTION("Error handling")
  {
    std::vector<int> strandLengths = { 10, 10 };

    // Test negative spacing
    CHECK_THROWS(generator.addRegularlySpacedFunctionalizedStrands(
      2, strandLengths, -1, 0, 4, 2, 1, true));

    // Test negative offset
    CHECK_THROWS(generator.addRegularlySpacedFunctionalizedStrands(
      2, strandLengths, 3, -1, 4, 2, 1, true));

    // Test zero spacing
    CHECK_THROWS(generator.addRegularlySpacedFunctionalizedStrands(
      2, strandLengths, 0, 0, 4, 2, 1, true));

    // Test inconsistent sizes
    CHECK_THROWS(generator.addRegularlySpacedFunctionalizedStrands(
      3, strandLengths, 3, 0, 4, 2, 1, true));

    // Test zero strands
    CHECK_THROWS(generator.addRegularlySpacedFunctionalizedStrands(
      0, {}, 3, 0, 4, 2, 1, true));
  }
}

TEST_CASE(
  "Universe generator with randomly functionalized chains reach correct w_sol",
  "[generator][MCUniverseGenerator][long]")
{
  std::cout << "Running test \"Universe generator with randomly functionalized "
               "chains reach correct w_sol\""
            << std::endl;

  pu::MCUniverseGenerator generator =
    pu::MCUniverseGenerator(35.375493, 35.375493, 35.375493);
  generator.setSeed(8804);
  generator.setMeanSquaredBeadDistance(1.107008);
  generator.configNrOfMCSteps(0);

  std::vector<int> chainLengths = pu::initializeWithValue(202, 210);
  generator.addRandomlyFunctionalizedStrands(
    202, chainLengths, 0.023809523, 1, 2, 1, true);
  generator.addStrands(212, 107, 3);
  generator.addMonofunctionalStrands(586, 57, 4);

  generator.useZScoreMaxDistance(3.29, 1.107008);
  generator.linkStrandsToSolubleFraction(0.31);

  pe::Universe universe = generator.getUniverse();

  size_t nFreeStrands = 0;
  for (const pe::Molecule& mol : universe.getChainsWithCrosslinker(2)) {
    nFreeStrands += mol.getType() == pe::MoleculeType::FREE_CHAIN;
  }

  CHECK(generator.getCurrentNrOfAtoms() == universe.getNrOfAtoms());
  CHECK(generator.getCurrentNrOfBonds() == universe.getNrOfBonds());

  CHECK(universe.getAtomsOfType(3).size() == 107 * 212);
  CHECK(universe.getAtomsOfType(4).size() == 57 * 586);

  pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
    pylimer_tools::sim::mehp::MEHPForceRelaxation(
      universe, 2, false, nullptr, 1.0, false, false);
  forceRelaxer.configAssumeBoxLargeEnough(true);

  pylimer_tools::sim::mehp::MEHPForceRelaxation relaxerFromGenerator =
    pylimer_tools::sim::mehp::MEHPForceRelaxation(
      generator.convertToForceRelaxationNetwork());
  relaxerFromGenerator.configAssumeBoxLargeEnough(true);

  // comparison only after running force relaxation,
  // since the box offsets, single nodes and stuff is relevant before that
  // CHECK(relaxer.countActiveClusteredAtoms() ==
  //       Catch::Approx(relaxerFromGenerator.countActiveClusteredAtoms()));
  // CHECK(forceRelaxer.getNrOfNodes() == relaxerFromGenerator.getNrOfNodes());
  CHECK(forceRelaxer.getNrOfSprings() + nFreeStrands ==
        relaxerFromGenerator.getNrOfSprings());

  while (forceRelaxer.suggestsRerun()) {
    forceRelaxer.runForceRelaxation("LD_MMA", 5000, 1e-11, 1e-8);
  }
  while (relaxerFromGenerator.suggestsRerun()) {
    relaxerFromGenerator.runForceRelaxation("LD_MMA", 5000, 1e-11, 1e-8);
  }

  CHECK_THAT(
    0.31,
    Catch::Matchers::WithinRel(forceRelaxer.getSolubleWeightFraction(), 0.05));
  CHECK_THAT(0.31,
             Catch::Matchers::WithinRel(
               relaxerFromGenerator.getSolubleWeightFraction(), 0.05));

  CHECK(forceRelaxer.countActiveClusteredAtoms() ==
        Catch::Approx(relaxerFromGenerator.countActiveClusteredAtoms()));

  // then, check that soluble fraction is removed correctly
  generator.removeSolubleFraction(true);
  pe::Universe universe2 = generator.getUniverse();

  CHECK(universe2.getNrOfAtoms() < universe.getNrOfAtoms());
  CHECK(generator.getCurrentNrOfAtoms() == universe2.getNrOfAtoms());

  pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer2 =
    pylimer_tools::sim::mehp::MEHPForceRelaxation(universe2);

  while (forceRelaxer2.suggestsRerun()) {
    forceRelaxer2.runForceRelaxation();
  }

  CHECK_THAT(
    0.,
    Catch::Matchers::WithinAbs(forceRelaxer2.getSolubleWeightFraction(), 0.05));
}

TEST_CASE("Universe generator with randomly functionalized chains use "
          "appropriate bond lengths",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"Universe generator with randomly functionalized "
               "chains use appropriate bond lengths\""
            << std::endl;

  pu::MCUniverseGenerator generator =
    pu::MCUniverseGenerator(35.375493, 35.375493, 35.375493);
  generator.setSeed(8804);
  constexpr double meanSquaredB = 1.107008;
  generator.setMeanSquaredBeadDistance(meanSquaredB);
  generator.configNrOfMCSteps(0);

  const std::vector<int> chainLengths = pu::initializeWithValue(50, 50);

  SECTION("Randomly functionalized chains")
  {
    generator.addRandomlyFunctionalizedStrands(
      chainLengths.size(), chainLengths, 0.7, 8, 2, 1, true);
    REQUIRE_NOTHROW(generator.validateInternalState());

    const pe::Universe universe = generator.getUniverse();
    std::vector<double> bondLengths = universe.computeBondLengths();
    const double maxBondLength = *std::ranges::max_element(bondLengths);
    CHECK(maxBondLength < 5. * meanSquaredB);
  }

  SECTION("End-functionalized chains")
  {
    generator.addCrosslinkStrands(
      chainLengths.size(), chainLengths, 8, 2, 1, true);
    REQUIRE_NOTHROW(generator.validateInternalState());

    const pe::Universe universe = generator.getUniverse();
    std::vector<double> bondLengths = universe.computeBondLengths();
    const double maxBondLength = *std::ranges::max_element(bondLengths);
    CHECK(maxBondLength < 5. * meanSquaredB);
  }
}

TEST_CASE("Universe generator uses correct w_sol even for strange structures",
          "[generator][MCUniverseGenerator][long]")
{
  std::cout
    << "Running test \"Universe generator uses correct w_sol even for strange "
       "structures\""
    << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(
    76.21419834207877, 76.21419834207877, 76.21419834207877);
  generator.setSeed(8804);
  generator.setMeanSquaredBeadDistance(1.107008);
  generator.configNrOfMCSteps(0);

  const std::vector<int> randomFchainLengths =
    pu::initializeWithValue(2020, 210);
  generator.addRandomlyFunctionalizedStrands(randomFchainLengths.size(),
                                             randomFchainLengths,
                                             0.024285714285714285,
                                             1,
                                             2,
                                             1,
                                             true);

  const std::vector<int> chainLengths = pu::initializeWithValue(2121, 107);
  generator.addStrands(chainLengths.size(), chainLengths, 1);

  {
    pylimer_tools::sim::mehp::MEHPForceBalance2 fb2_initial =
      generator.getForceBalance2();
    fb2_initial.runForceRelaxation();
    CHECK_THAT(fb2_initial.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(1.0, 0.05));
  }

  generator.useZScoreMaxDistance(3., 1.107008);
  generator.linkStrandsToSolubleFraction(0.31);

  const pylimer_tools::sim::mehp::ForceBalanceNetwork net1 =
    generator.getForceBalance().getNetwork();

  const pe::Universe universe = generator.getUniverse();

  double nBondsInStrands = net1.springsContourLength.sum();
  double nStrands = net1.springsContourLength.size();
  CHECK(nStrands == net1.nrOfSprings);
  CHECK(nBondsInStrands - nStrands + net1.nrOfLinks == universe.getNrOfAtoms());
  CHECK(nBondsInStrands == universe.getNrOfBonds());
  CHECK(generator.getCurrentNrOfAtoms() == universe.getNrOfAtoms());
  CHECK(generator.getCurrentNrOfBonds() == universe.getNrOfBonds());

  {
    // SCOPE 1:
    INFO("Force Balance 2 based on Universe");
    pylimer_tools::sim::mehp::MEHPForceBalance2 forceBalance =
      pylimer_tools::sim::mehp::MEHPForceBalance2(universe);

    forceBalance.runForceRelaxation();
    CHECK(forceBalance.inferNrOfAtomsFromNetwork() == universe.getNrOfAtoms());

    CHECK_THAT(forceBalance.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.31, 0.05));
    const pylimer_tools::sim::mehp::ForceBalance2Network net2 =
      forceBalance.getNetwork();
    CHECK(net2.springContourLength.sum() - net2.nrOfSprings + net2.nrOfLinks ==
          universe.getNrOfAtoms());
    CHECK(forceBalance.inferNrOfAtomsFromNetwork(net2) ==
          universe.getNrOfAtoms());
    forceBalance.configAssumeNetworkIsComplete(true);
    CHECK_THAT(forceBalance.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.31, 0.05));
  }

  {
    // SCOPE 2:
    INFO("Force Balance 2 based on force balance 1 network from generator");

    pylimer_tools::sim::mehp::MEHPForceBalance2 forceBalance =
      pylimer_tools::sim::mehp::MEHPForceBalance2(
        net1, generator.getForceBalance().getSpringPartitions());
    forceBalance.configAssumeNetworkIsComplete(true);

    forceBalance.runForceRelaxation();
    CHECK_THAT(forceBalance.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.31, 0.05));
  }

  {
    // SCOPE 3:
    INFO("Force Balance 2 based on force relaxation network from generator");

    pylimer_tools::sim::mehp::MEHPForceBalance2 forceBalance =
      pylimer_tools::sim::mehp::MEHPForceBalance2(
        generator.getForceRelaxation().getNetwork());

    forceBalance.runForceRelaxation();
    CHECK_THAT(forceBalance.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.31, 0.05));
  }

  {
    // SCOPE 4:
    INFO("Force Relaxation based on network from generator");

    pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
      generator.getForceRelaxation();
    forceRelaxer.configAssumeBoxLargeEnough(false);
    forceRelaxer.configAssumeNetworkIsComplete(true);

    while (forceRelaxer.suggestsRerun()) {
      forceRelaxer.runForceRelaxation();
    }

    CHECK_THAT(forceRelaxer.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.31, 0.05));

    INFO("Now with assuming box large enough");
    forceRelaxer.configAssumeBoxLargeEnough(true);
    while (forceRelaxer.suggestsRerun()) {
      forceRelaxer.runForceRelaxation();
    }

    CHECK_THAT(forceRelaxer.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.31, 0.05));
  }

  {
    // SCOPE 5:
    INFO("Force Relaxation based on Universe");

    pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
      pylimer_tools::sim::mehp::MEHPForceRelaxation(universe, 2);
    forceRelaxer.configAssumeBoxLargeEnough(false);

    while (forceRelaxer.suggestsRerun()) {
      forceRelaxer.runForceRelaxation();
    }

    CHECK_THAT(forceRelaxer.getSolubleWeightFraction(),
               Catch::Matchers::WithinAbs(0.31, 0.05));
  }
}

TEST_CASE("Universe generator handles structures without crosslinks",
          "[generator][MCUniverseGenerator]")
{
  std::cout << "Running test \"Universe generator handles structures without "
               "crosslinks\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(
    76.21419834207877, 76.21419834207877, 76.21419834207877);
  generator.setSeed(8804);
  generator.setMeanSquaredBeadDistance(1.107008);
  generator.configNrOfMCSteps(0);

  const std::vector<int> chainLengths = pu::initializeWithValue(50, 50);
  generator.addStrands(chainLengths.size(), chainLengths, 1);

  CHECK_NOTHROW(generator.validateInternalState());
  CHECK_NOTHROW(generator.linkStrandsToSolubleFraction(1.));
  CHECK_THROWS(generator.linkStrandsToConversion(0.));
  CHECK_THROWS(generator.linkStrandsToConversion(1.));
}

TEST_CASE("Linear walk chain can be generated", "[topo][RandomWalker]")
{
  std::cout << "Running test \"Linear walk chain can be generated\""
            << std::endl;
  pe::Box box = pe::Box(10., 10., 10.);

  SECTION("With ends")
  {
    Eigen::VectorXd coordinates = pu::doLinearWalkChainFromTo(
      box, Eigen::Vector3d(0., 0., 0.), Eigen::Vector3d(5., 5., 5.), 10, true);
    CHECK(coordinates.size() == (10 + 2) * 3);
    CHECK(coordinates.segment(0, 3) == Eigen::Vector3d(0., 0., 0.));
    CHECK(coordinates.segment(3 * (10 + 1), 3) == Eigen::Vector3d(5., 5., 5.));
    Eigen::Vector3d prevDistance =
      coordinates.segment(0, 3) - coordinates.segment(3, 3);
    for (size_t i = 3; i < coordinates.size(); i += 3) {
      Eigen::Vector3d currentDistance =
        coordinates.segment(i, 3) - coordinates.segment(i - 3, 3);
      CHECK(currentDistance.norm() == Catch::Approx(prevDistance.norm()));
    }
  }

  SECTION("Without ends")
  {
    Eigen::VectorXd coordinates = pu::doLinearWalkChainFromTo(
      box, Eigen::Vector3d(0., 0., 0.), Eigen::Vector3d(5., 5., 5.), 10, false);
    CHECK(coordinates.size() == (10) * 3);
    Eigen::Vector3d prevDistance =
      coordinates.segment(0, 3) - coordinates.segment(3, 3);
    for (size_t i = 3; i < coordinates.size(); i += 3) {
      Eigen::Vector3d currentDistance =
        coordinates.segment(i, 3) - coordinates.segment(i - 3, 3);
      CHECK(currentDistance.norm() == Catch::Approx(prevDistance.norm()));
    }
  }
}

TEST_CASE("randomWalkChain generates valid chains", "[topo][RandomWalker]")
{
  std::cout << "Running test \"randomWalkChain generates valid chains\"";

  // Create a box
  pylimer_tools::entities::Box box(10.0, 10.0, 10.0);

  // Parameters for the random walk
  int chainLen = 30;
  double beadDistance = 1.0;
  double meanSquaredBeadDistance = (3. * M_PI / 8.) * beadDistance;
  std::string seed = "test_seed";

  Eigen::VectorXd coordinates = pylimer_tools::utils::doRandomWalkChain(
    chainLen, beadDistance, meanSquaredBeadDistance, seed);

  // Check the size of the returned vector (should be 3 * chainLen)
  REQUIRE(coordinates.size() == 3 * chainLen);

  // Check that bond lengths are reasonable
  Eigen::VectorXd bondLengths =
    coordinates.tail((chainLen - 1) * 3) - coordinates.head((chainLen - 1) * 3);
  box.handlePBC(bondLengths);
  CHECK(andersonDarlingNormalDistributionTest(bondLengths, 0., M_PI / 8.));
}

TEST_CASE("doRandomWalkChainFromToMC generates valid chains",
          "[topo][RandomWalker][mc]")
{
  std::cout
    << "Running test \"doRandomWalkChainFromToMC generates valid chains\""
    << std::endl;
  // Create a box
  pylimer_tools::entities::Box box(10.0, 10.0, 10.0);

  // Define start and end points
  Eigen::Vector3d from(1.0, 1.0, 1.0);
  Eigen::Vector3d to(1.0, 1.0, 1.0);

  // Parameters for the random walk
  int chainLen = 30;
  double beadDistance = 1.0;
  double meanSquaredBeadDistance = (3. * M_PI / 8.) * beadDistance;
  std::string seed = "let_s_use_this_seed";
  // Parameters for MC
  int numIterations = 500;

  // Generate the chain
  Eigen::VectorXd coordinatesNoMC =
    pylimer_tools::utils::doRandomWalkChainFromTo(
      box, from, to, chainLen, beadDistance, meanSquaredBeadDistance, seed);
  Eigen::VectorXd coordinates =
    pylimer_tools::utils::doRandomWalkChainFromToMC(box,
                                                    from,
                                                    to,
                                                    chainLen,
                                                    beadDistance,
                                                    meanSquaredBeadDistance,
                                                    seed,
                                                    numIterations);

  // Check the size of the returned vector (should be 3 * chainLen)
  REQUIRE(coordinatesNoMC.size() == 3 * chainLen);
  REQUIRE(coordinates.size() == 3 * chainLen);

  // Check that bond lengths are reasonable
  Eigen::VectorXd bondLengthsNoMC = coordinatesNoMC.tail((chainLen - 1) * 3) -
                                    coordinatesNoMC.head((chainLen - 1) * 3);
  box.handlePBC(bondLengthsNoMC);
  CHECK(andersonDarlingNormalDistributionTest(bondLengthsNoMC, 0., M_PI / 8.));
  Eigen::VectorXd bondLengths =
    coordinates.tail((chainLen - 1) * 3) - coordinates.head((chainLen - 1) * 3);
  box.handlePBC(bondLengths);
  CHECK(andersonDarlingNormalDistributionTest(bondLengths, 0., M_PI / 8.));
}

TEST_CASE("doRandomWalkChain generates very long chains",
          "[topo][RandomWalker][mc]")
{
  // Create a box
  pylimer_tools::entities::Box box(100.0, 100.0, 100.0);

  // Define start and end points
  Eigen::Vector3d from(1.0, 1.0, 1.0);

  // Parameters for the random walk
  constexpr int chainLen = 32;
  constexpr double beadDistance = 1.0;
  constexpr double meanSquaredBeadDistance = (3. * M_PI / 8.) * beadDistance;
  const std::string seed = "some_other_seed";
  int numIterations = 500;

  // Generate the chain
  Eigen::VectorXd coordinates = pylimer_tools::utils::doRandomWalkChain(
    chainLen, beadDistance, meanSquaredBeadDistance, seed);

  // Check the size of the returned vector (should be 3 * chainLen)
  REQUIRE(coordinates.size() == 3 * chainLen);

  // Check that bond lengths are reasonable
  const Eigen::VectorXd bondLengths =
    coordinates.tail((chainLen - 1) * 3) - coordinates.head((chainLen - 1) * 3);
  CHECK(andersonDarlingNormalDistributionTest(bondLengths, 0., M_PI / 8.));
}

/**
 * Tests for new tetra-PEG network generation methods
 */

TEST_CASE("addStarCrosslinkers generates star-like crosslinkers with "
          "pre-connected strands",
          "[generator][MCUniverseGenerator][tetra-peg]")
{
  std::cout << "Running test \"addStarCrosslinkers generates star-like "
               "crosslinkers with pre-connected strands\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
  generator.setSeed(12345);
  generator.setBeadDistance(0.964);
  generator.configNrOfMCSteps(0);

  SECTION("Basic star crosslinker generation")
  {
    // Add star crosslinkers with f=6 functionality and strand length 10
    generator.addStarCrosslinkers(50, 6, 10, 2, 1);

    pe::Universe universe = generator.getUniverse();

    // Check total number of atoms: 50 crosslinkers + 50*6*10 strand atoms
    CHECK(universe.getNrOfAtoms() == 50 + 50 * 6 * 10);
    CHECK(universe.getAtomsOfType(2).size() == 50);          // crosslinkers
    CHECK(universe.getAtomsOfType(1).size() == 50 * 6 * 10); // strand atoms

    // Check number of bonds: 50*6*9 bonds in strands + 50*6 crosslinker-strand
    // bonds
    CHECK(universe.getNrOfBonds() == 50 * 6 * 10);

    // Check that we have exactly 50*6 molecules (6 per star)
    auto molecules = universe.getMolecules(2);
    CHECK(molecules.size() == 50 * 6);

    // Check that each molecule has the correct number of atoms
    for (const auto& molecule : molecules) {
      CHECK(molecule.getNrOfAtoms() == 10); // 10 atoms per strand
    }

    // then, check that we have a cross-link for each strand
    auto chains = universe.getChainsWithCrosslinker(2);
    CHECK(chains.size() == 50 * 6);

    for (const auto& chain : chains) {
      CHECK(chain.getNrOfAtoms() == 11); // 10 atoms per strand + 1 crosslinker
    }
  }

  SECTION("Star crosslinkers with different functionality")
  {
    // Add star crosslinkers with f=8 functionality
    generator.addStarCrosslinkers(25, 8, 15, 3, 2);

    pe::Universe universe = generator.getUniverse();

    CHECK(universe.getNrOfAtoms() == 25 + 25 * 8 * 15);
    CHECK(universe.getAtomsOfType(3).size() == 25);          // crosslinkers
    CHECK(universe.getAtomsOfType(2).size() == 25 * 8 * 15); // strand atoms

    auto molecules = universe.getMolecules(3);
    CHECK(molecules.size() == 25 * 8);

    for (const auto& molecule : molecules) {
      CHECK(molecule.getNrOfAtoms() == 15); // 15 atoms per strand
    }
  }

  SECTION("Error handling for invalid parameters")
  {
    // Test error cases
    CHECK_THROWS(
      generator.addStarCrosslinkers(-1, 4, 10, 2, 1)); // negative number
    CHECK_THROWS(
      generator.addStarCrosslinkers(10, 1, 10, 2, 1)); // functionality < 2
    CHECK_THROWS(
      generator.addStarCrosslinkers(10, 4, 0, 2, 1)); // zero strand length
  }
}

TEST_CASE(
  "linkStrandToStrand connects strand ends based on distance probability",
  "[generator][MCUniverseGenerator][tetra-peg]")
{
  std::cout << "Running test \"linkStrandToStrand connects strand ends based "
               "on distance probability\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(15.0, 15.0, 15.0);
  generator.setSeed(54321);
  generator.setBeadDistance(0.964);
  generator.configNrOfMCSteps(0);

  SECTION("Basic strand end linking")
  {
    // First add 20 star crosslinkers to have free strand ends
    generator.addStarCrosslinkers(20, 4, 8, 2, 1);

    pe::Universe universeBefore = generator.getUniverse();
    size_t bondsBefore = universeBefore.getNrOfBonds();

    // Link strand ends - should create some new bonds
    bool linkMade = generator.linkStrandToStrand();

    pe::Universe universeAfter = generator.getUniverse();
    size_t bondsAfter = universeAfter.getNrOfBonds();

    // Check that a link was potentially made
    CHECK(bondsAfter == bondsBefore + 1);
    CHECK(universeAfter.getNrOfAtoms() ==
          universeBefore.getNrOfAtoms()); // no new atoms created
  }

  SECTION("Linking with custom cInfinity parameter")
  {
    // Add strands with free ends
    generator.addStarCrosslinkers(15, 6, 12, 2, 1);

    size_t bondsBefore = generator.getUniverse().getNrOfBonds();
    bool linkMade =
      generator.linkStrandToStrand(2.0); // custom cInfinity parameter
    size_t bondsAfter = generator.getUniverse().getNrOfBonds();

    CHECK(bondsAfter == bondsBefore + 1);
  }

  SECTION("No linking when no free strand ends exist")
  {
    // Add regular crosslinkers without free strand ends
    generator.addCrosslinkers(10, 4, 2);

    size_t bondsBefore = generator.getUniverse().getNrOfBonds();
    bool linkMade = generator.linkStrandToStrand();
    size_t bondsAfter = generator.getUniverse().getNrOfBonds();

    // Should not create any new bonds since there are no free strand ends
    CHECK_FALSE(linkMade);
    CHECK(bondsAfter == bondsBefore);
  }
}

TEST_CASE("linkStrandsToStrandsToConversion achieves target conversion through "
          "iterative linking",
          "[generator][MCUniverseGenerator][tetra-peg]")
{
  std::cout
    << "Running test \"linkStrandsToStrandsToConversion achieves target "
       "conversion through iterative linking\""
    << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(20.0, 20.0, 20.0);
  generator.setSeed(98765);
  generator.setBeadDistance(0.964);
  generator.configNrOfMCSteps(0);

  SECTION("Achieving target conversion")
  {
    // Add 30 star crosslinkers to create many free strand ends
    generator.addStarCrosslinkers(30, 6, 10, 2, 1);

    pe::Universe universeBefore = generator.getUniverse();
    size_t bondsBefore = universeBefore.getNrOfBonds();
    CHECK(bondsBefore == 30 * (6 * 10));

    // Check that the current conversion rate is accurate
    CHECK_THAT(generator.getCurrentStrandsConversion(),
               Catch::Matchers::WithinRel(0.5, 1e-4));
    // since this is the current conversion, no extra linking should occur
    generator.linkStrandsToStrandsToConversion(0.5);
    CHECK(generator.getUniverse().getNrOfBonds() == bondsBefore);

    // Target 50% conversion of free strand ends,
    // while 50% is already converted by being connected to the crosslinkers
    double targetConversion = 0.5 + 0.25;
    generator.linkStrandsToStrandsToConversion(targetConversion);

    pe::Universe universeAfter = generator.getUniverse();
    size_t bondsAfter = universeAfter.getNrOfBonds();

    double nStrands = 30 * 6;
    CHECK(universeBefore.getMolecules(2).size() == nStrands);

    // With 30 stars * 6 strands = 180 strands, total 360 strand ends
    // Starting conversion: 0.5 (180 ends connected to crosslinkers, 180 free)
    // Target conversion: 0.75 (270 ends connected, 90 free)
    // Need to connect 90 more ends, which means 45 strand-to-strand bonds
    CHECK(bondsAfter == bondsBefore + 45);
    CHECK(universeAfter.getMolecules(2).size() == nStrands * targetConversion);
    CHECK(generator.findFreeStrandEnds().size() ==
          (1. - targetConversion) * 30 * 6 * 2);
  }

  SECTION("High conversion rate")
  {
    // Test with higher conversion rate
    generator.addStarCrosslinkers(20, 4, 8, 2, 1);

    size_t bondsBefore = generator.getUniverse().getNrOfBonds();
    CHECK(bondsBefore == 20 * (4 * 8));

    // Target 80% conversion.
    // Again, 50% is already converted by being connected to the crosslinkers
    double targetConversion = 0.8 * 0.5 + 0.5;
    generator.linkStrandsToStrandsToConversion(targetConversion);

    size_t bondsAfter = generator.getUniverse().getNrOfBonds();

    CHECK_THAT(bondsAfter - bondsBefore,
               Catch::Matchers::WithinRel(20 * 4 * 0.8 / 2., 1e-4));
    CHECK_THAT(
      generator.findFreeStrandEnds().size(),
      Catch::Matchers::WithinRel((1. - targetConversion) * 20 * 4 * 2, 1e-4));
  }

  SECTION("Custom cInfinity parameter")
  {
    generator.addStarCrosslinkers(25, 4, 6, 2, 1);

    size_t bondsBefore = generator.getUniverse().getNrOfBonds();
    generator.linkStrandsToStrandsToConversion(0.6, 2.0); // custom cInfinity
    size_t bondsAfter = generator.getUniverse().getNrOfBonds();

    CHECK(bondsAfter >= bondsBefore);
  }

  SECTION("Error handling")
  {
    generator.addStarCrosslinkers(10, 4, 5, 2, 1);

    // Test invalid conversion values
    CHECK_THROWS(
      generator.linkStrandsToStrandsToConversion(-0.1)); // negative conversion
    CHECK_THROWS(
      generator.linkStrandsToStrandsToConversion(1.1)); // conversion > 1
  }
}

TEST_CASE("linkStrandsToStrandsToConversion also works with free strands",
          "[generator][MCUniverseGenerator][tetra-peg]")
{
  std::cout << "Running test \"linkStrandsToStrandsToConversion also works "
               "with free strands\""
            << std::endl;

  pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(80.0, 80.0, 80.0);
  generator.setSeed(12345);

  generator.setBeadDistance(1.0);

  generator.configNrOfMCSteps(0);

  generator.addStarCrosslinkers(150, 4, 20, 2, 1);

  generator.addStrands(50, 15, 1);

  CHECK_NOTHROW(generator.linkStrandsToStrandsToConversion(0.9));
}

#ifdef CEREALIZABLE
TEST_CASE(
  "MaxDistanceProvider implementations can be serialized and deserialized",
  "[generator][MaxDistanceProvider][serialization]")
{
  std::cout << "Running test \"MaxDistanceProvider implementations can be "
               "serialized and deserialized\""
            << std::endl;

  SECTION("LinearMaxDistanceProvider")
  {
    pu::LinearMaxDistanceProvider original = pu::LinearMaxDistanceProvider(2.5);
    CHECK(original.getMaxDistance(10.0) == Catch::Approx(25.0));

    std::string serialized = pu::serializeToString(original);
    pu::LinearMaxDistanceProvider deserialized;
    pu::deserializeFromString(deserialized, serialized);

    CHECK(deserialized.getMaxDistance(10.0) == Catch::Approx(25.0));
  }

  SECTION("ZScoreMaxDistanceProvider")
  {
    pu::ZScoreMaxDistanceProvider original =
      pu::ZScoreMaxDistanceProvider(3.0, 1.5);
    CHECK(original.getStdMultiplier() == Catch::Approx(3.0));
    CHECK(original.getInnerMultiplier() == Catch::Approx(1.5));
    CHECK(original.getMaxDistance(16.0) ==
          Catch::Approx(3.0 * std::sqrt(16.0 * 1.5)));

    std::string serialized = pu::serializeToString(original);
    pu::ZScoreMaxDistanceProvider deserialized;
    pu::deserializeFromString(deserialized, serialized);

    CHECK(deserialized.getMaxDistance(16.0) ==
          Catch::Approx(3.0 * std::sqrt(16.0 * 1.5)));
  }

  SECTION("NoMaxDistanceProvider")
  {
    pu::NoMaxDistanceProvider original = pu::NoMaxDistanceProvider();
    CHECK(original.getMaxDistance(10.0) == Catch::Approx(-1.0));

    std::string serialized = pu::serializeToString(original);
    pu::NoMaxDistanceProvider deserialized;
    pu::deserializeFromString(deserialized, serialized);

    CHECK(deserialized.getMaxDistance(10.0) == Catch::Approx(-1.0));
  }

  // SECTION("Polymorphic serialization")
  // {
  //   pu::ZScoreMaxDistanceProvider original =
  //   pu::ZScoreMaxDistanceProvider(2.0, 0.8);
  //   CHECK(original.getMaxDistance(25.0) == Catch::Approx(2.0 * std::sqrt(25.0
  //   * 0.8)));
  //
  //   std::string serialized = pu::serializeToString(original);
  //   std::unique_ptr<pu::MaxDistanceProvider> deserialized;
  //   pu::deserializeFromString(deserialized, serialized);
  //
  //   CHECK(deserialized->getMaxDistance(25.0) == Catch::Approx(2.0 *
  //   std::sqrt(25.0 * 0.8)));
  // }
}
#endif
