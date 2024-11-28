#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/io/DataFileWriter.h"
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
    generator.addStrands((4 / 2) * 100, 16);
    generator.configNrOfMCSteps(10);
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

    SECTION("Nrs of chains is correct") {
        REQUIRE(universe.getAtomsOfType(2).size() == 100);
        REQUIRE(universe.getAtomsOfType(1).size() == (4 / 2) * 100 * 16);
        REQUIRE(universe.getMolecules(2).size() == (4 / 2) * 100 + 100);
    }

    SECTION("Universe is generated deterministically") {
        pu::MCUniverseGenerator generator2 =
            pu::MCUniverseGenerator(10.0, 10.0, 10.0);
        generator2.setSeed(8804);
        generator2.setBeadDistance(0.964);
        generator2.addCrosslinkers(100, 4, 2);
        generator2.addSolventChains(100, 16, 3);
        generator2.addStrands((4 / 2) * 100, 16);
        generator2.configNrOfMCSteps(10);
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
            CHECK(bondLength < 3.2);
        }
    }

    SECTION("Errors are thrown") {
        // nr of strands and strand lengths must be same:
        REQUIRE_THROWS(generator.addStrands(3, { { 10, 100 } }, 1));
        // not enough strands to reach conversion:
        generator.addStrands(2, 10, 1);
        REQUIRE_THROWS(generator.linkStrandsToConversion(2.0));
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
    generator.addStrands(200, 10, 1);
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

TEST_CASE("Universe can cross-link up to w_sol",
          "[generator][MCUniverseGenerator]")
{
    std::cout << "Running test \"Universe can cross-link up to w_sol\""
              << std::endl;

    pu::MCUniverseGenerator generator = pu::MCUniverseGenerator(10.0, 10.0, 10.0);
    generator.setSeed(8804);
    generator.setBeadDistance(0.964);
    generator.addCrosslinkers(100, 4, 2);

    generator.addStrands(200, 19, 1);
    generator.linkStrandsToSolubleFraction(0.1, 1.);

    pe::Universe universe = generator.getUniverse();
    REQUIRE(universe.getVolume() == 10.0 * 10.0 * 10.0);
    REQUIRE(universe.getAtomsOfType(2).size() == 100);
    REQUIRE(universe.getAtomsOfType(1).size() == 200 * 19);

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

    SECTION("Without rescaling") {
        generator.removeSolubleFraction(false);

        pe::Universe universeAfterRemoval = generator.getUniverse();
        CHECK(universeAfterRemoval.getVolume() ==
                  universeBeforeRemoval.getVolume());

        CHECK(universeAfterRemoval.getAtomsOfType(2).size() < 400);
        CHECK(universeAfterRemoval.getAtomsOfType(1).size() < 800 * 10);
    }

    SECTION("With rescaling") {
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
    generator.linkStrandsToConversion(0.925, 1.);

    pe::Universe universe = generator.getUniverse();
    CHECK(universe.getNrOfAtoms() == 400 + 800 * 10);
    CHECK(generator.getCurrentNrOfAtoms() == universe.getNrOfAtoms());
    CHECK(generator.getCurrentNrOfBonds() == universe.getNrOfBonds());

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
    CHECK(relaxer.getNrOfSprings() == relaxerFromGenerator.getNrOfSprings());

    relaxer.runForceRelaxation();
    relaxerFromGenerator.runForceRelaxation();

    CHECK(relaxer.countActiveClusteredAtoms() ==
              Catch::Approx(relaxerFromGenerator.countActiveClusteredAtoms()));
    // cannot compare, since each strand end becomes a node in the from-universe
    // generation
    // CHECK(relaxer.getNrOfNodes() == relaxerFromGenerator.getNrOfNodes());
    CHECK(relaxer.getNrOfSprings() == relaxerFromGenerator.getNrOfSprings());
}

TEST_CASE("MUniverseGenerator can generate with cross-link chains",
          "[generator][MCUniverseGenerator]")
{
    std::cout
            << "Running test \"MUniverseGenerator can generate with cross-link chains\""
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
    for (auto chain : chains) {
        CHECK(chain.getNrOfAtoms() == 12);
    }
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

    SECTION("Without functionalization") {
        generator.addRandomlyFunctionalizedStrands(
            10, { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 }, 0., 2, 1, true);

        pe::Universe universe = generator.getUniverse();
        CHECK(universe.getAtomsOfType(2).size() == 0);
        CHECK(universe.getAtomsOfType(1).size() == 10 * 10);

        std::vector<pe::Molecule> chains = universe.getChainsWithCrosslinker(2);
        CHECK(chains.size() == 10);
        for (auto chain : chains) {
            CHECK(chain.getNrOfAtoms() == 10);
        }
    }

    SECTION("With functionalization") {
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

    generator.linkStrandsToSolubleFraction(0.31);

    pe::Universe universe = generator.getUniverse();

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
    CHECK(forceRelaxer.getNrOfSprings() == relaxerFromGenerator.getNrOfSprings());

    while (forceRelaxer.suggestsRerun()) {
        forceRelaxer.runForceRelaxation("LD_MMA", 5000, 1e-11, 1e-8);
    }
    while (relaxerFromGenerator.suggestsRerun()) {
        relaxerFromGenerator.runForceRelaxation("LD_MMA", 5000, 1e-11, 1e-8);
    }

    CHECK_THAT(
        0.31,
        Catch::Matchers::WithinRel(forceRelaxer.getSolubleWeightFraction(), 0.05));

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
