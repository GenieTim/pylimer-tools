#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <iostream>
#include <map>
#include <vector>

extern "C" {
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;

// void
// outputLoops(std::map<int, std::vector<std::vector<pe::Atom>>> loops)
// {
//   for (const auto& myPair : loops) {
//     std::cout << myPair.first << std::endl;
//     for (const std::vector<pe::Atom>& as : myPair.second) {
//       for (const pe::Atom& a : as) {
//         std::cout << "\t" << a.getId() << std::endl;
//       }
//     }
//   }
// }

template<typename T>
std::vector<T>
getVectorWithOne(T val) {
    std::vector<T> vec;
    vec.push_back(val);
    return vec;
}

TEST_CASE("Universe can be created", "[entity][Universe]") {
    std::cout << "Running test \"Universe can be created\"" << std::endl;

    pe::Universe universe(1.0, 1.0, 1.0);

    SECTION("atoms can be added") {
        universe.addAtoms({{0, 1}},
                          {{1, 1}},
                          {{0.0, 1.0}},
                          {{0.0, 1.0}},
                          {{0.0, 1.0}},
                          {{0, 0}},
                          {{0, 0}},
                          {{0, 0}});
        CHECK(universe.getNrOfAtoms() == 2);
        universe.addAtoms({{3, 4}},
                          {{1, 1}},
                          {{0.0, 1.0}},
                          {{0.0, 1.0}},
                          {{0.0, 1.0}},
                          {{0, 0}},
                          {{0, 0}},
                          {{0, 0}});
        CHECK(universe.getNrOfAtoms() == 4);
        CHECK(universe.getAtomsOfType(1).size() == 4);
        CHECK(universe.getNrOfBonds() == 0);
        REQUIRE_THROWS(universe.getIdxByAtomId(1000));

        SECTION("Atom coordinates can be rescaled") {
            CHECK(universe.getVolume() == 1.0);
            CHECK(universe.getAtom(0).getX() == 0.0);
            CHECK(universe.getAtom(1).getX() == 1.0);
            universe.setBox(pe::Box(2.0, 2.0, 2.0), true);
            CHECK(universe.getAtom(0).getX() == 0.0);
            CHECK(universe.getAtom(0).getY() == 0.0);
            CHECK(universe.getAtom(0).getZ() == 0.0);
            CHECK(universe.getAtom(1).getX() == 2.0);
            CHECK(universe.getAtom(1).getY() == 2.0);
            CHECK(universe.getAtom(1).getZ() == 2.0);
            universe.setBox(pe::Box(1.0, 2.0, 1.0, 2.0, 1.0, 2.0), true);
            CHECK(universe.getAtom(0).getX() == 1.0);
            CHECK(universe.getAtom(0).getY() == 1.0);
            CHECK(universe.getAtom(0).getZ() == 1.0);
            CHECK(universe.getAtom(1).getX() == 2.0);
            CHECK(universe.getAtom(1).getY() == 2.0);
            CHECK(universe.getAtom(1).getZ() == 2.0);
        }

        SECTION("Atoms can be removed") {
            universe.addBonds({{0, 1, 3, 4}}, {{1, 3, 4, 0}});
            pe::Universe universeCopy = pe::Universe(universe);
            CHECK_NOTHROW(universeCopy.getAtom(3));
            CHECK(universeCopy.getNrOfAtoms() == 4);
            CHECK(universeCopy.getNrOfBonds() == 4);
            universeCopy.removeAtoms({{1, 3}});
            CHECK(universeCopy.getAtom(4) == universe.getAtom(4));
            CHECK_THROWS(universeCopy.getAtom(3));
            CHECK(universeCopy.validate());
            CHECK(universeCopy.getNrOfAtoms() == 2);
            CHECK(universeCopy.getNrOfBonds() == 1);
        }

        SECTION("molecules are found") {
            CHECK(universe.getNrOfAtoms() == 4);
            CHECK(universe.getMolecules(2).size() == 4);
            // works twice: make sure there is no removal of anything happening
            CHECK(universe.getMolecules(2).size() == 4);
            universe.addBonds(2, {{0, 3}}, {{3, 4}});
            CHECK(universe.getNrOfAtoms() == 4);
            CHECK(universe.getNrOfBonds() == 2);
            CHECK(universe.getMolecules(1).size() == 0);
            std::vector<pe::Molecule> molecules = universe.getMolecules(2);
            CHECK(molecules.size() == 2);
            // works twice: make sure there is no removal of anything happening
            CHECK(universe.getMolecules(2).size() == 2);
        }

        SECTION("Special constructors work") {
            pe::Universe universe2 = universe;
            CHECK(universe2.getNrOfAtoms() == 4);
            CHECK(universe2.getMolecules(2).size() == 4);
        }
    }

    SECTION("Disallowed mutations are detected") {
        std::vector<int> threeZeros = {{0, 0, 0}};
        std::vector<long int> oneTwoThree = {{1, 2, 3}};
        std::vector<long int> threeLongZeros = {{0, 0, 0}};
        std::vector<double> threeDoubleZeros = {{0, 0, 0}};
        // atom not yet added
        REQUIRE_THROWS(universe.addBonds(
            3, oneTwoThree, threeLongZeros, threeZeros, false, true));
        // same, but ignore the error
        REQUIRE_NOTHROW(universe.addBonds(
            3, oneTwoThree, threeLongZeros, threeZeros, true, true));
        // all fine
        REQUIRE_NOTHROW(universe.addAtoms(oneTwoThree,
            threeZeros,
            threeDoubleZeros,
            threeDoubleZeros,
            threeDoubleZeros,
            threeZeros,
            threeZeros,
            threeZeros));
        // id already exists
        REQUIRE_THROWS(universe.addAtoms(threeLongZeros,
            threeZeros,
            threeDoubleZeros,
            threeDoubleZeros,
            threeDoubleZeros,
            threeZeros,
            threeZeros,
            threeZeros));
        std::vector<long int> fourLongZeros = {{0, 0, 0, 0}};
        // different lengths
        REQUIRE_THROWS(universe.addAtoms(fourLongZeros,
            threeZeros,
            threeDoubleZeros,
            threeDoubleZeros,
            threeDoubleZeros,
            threeZeros,
            threeZeros,
            threeZeros));
        // different lengths
        REQUIRE_THROWS(
            universe.addAngles(oneTwoThree, fourLongZeros, oneTwoThree, threeZeros));
        // different lengths
        REQUIRE_THROWS(
            universe.addBonds(3, oneTwoThree, fourLongZeros, threeZeros));
        REQUIRE_THROWS(universe.computeDxs(threeLongZeros, fourLongZeros));
    }
}

TEST_CASE("Universe can be resized", "[entity][Universe]") {
    std::cout << "Running test \"Universe can be resized\"" << std::endl;

    pe::Universe universe = pe::Universe(1.0, 1.0, 1.0);
    REQUIRE(universe.getVolume() == 1.0);

    SECTION("resizing smaller changes volume") {
        REQUIRE_THROWS(pe::Box(0.0, 1.0, 2.0));
        universe.setBox(pe::Box(0.25, 1.0, 2.0));
        REQUIRE(universe.getVolume() == 0.5);
        // check empty stuff
        REQUIRE(universe.getClusters().size() == 0);
        REQUIRE(universe.getMolecules(2).size() == 0);
        REQUIRE(universe.getChainsWithCrosslinker(2).size() == 0);
        REQUIRE(universe.determineEffectiveFunctionalityPerType().size() == 0);
        REQUIRE(universe.computeWeightFractions().size() == 0);
        REQUIRE(universe.computeMeanBondLength() == 0.0);
        REQUIRE(universe.computeBondLengths().size() == 0);
        REQUIRE(universe.computeMeanEndToEndDistance(2) == 0.0);
        REQUIRE(universe.computeMeanSquareEndToEndDistance(2) == 0.0);
        REQUIRE(universe.computeEndToEndDistances(2).size() == 0);
        REQUIRE(universe.findLoops(2).size() == 0);
    }

    SECTION("resizing bigger changes volume") {
        universe.setBoxLengths(1.0, 1.0, 1.0);
        REQUIRE(universe.getVolume() == 1.0);
        universe.setBox(pe::Box(2.0, 1.0, 2.0));
        REQUIRE(universe.getVolume() == 4.0);
        universe.setBox(pe::Box(-1.0, 1.0, 0.0, 1.0, 0.0, 2.0));
        REQUIRE(universe.getVolume() == 4.0);
    }
}

TEST_CASE("Universe can be used", "[entity][Universe]") {
    std::cout << "Running test \"Universe can be used\"" << std::endl;

    pe::Universe universe = pe::Universe(1.0, 1.0, 1.0);

    SECTION("Loop Entanglements are detected") {
        /**
         * The system looks like this (in terms of bonds, not 3D placement):
         *
         * 1-2
         * | |
         * 4-3
         *
         * 5-6
         * | |
         * 8-7
         */
        universe.setBox(pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0));
        universe.addAtoms({{1, 2, 3, 4, 5, 6, 7, 8}},
                          // id
                          {{1, 1, 1, 2, 1, 2, 2, 1}},
                          // type
                          {{-5, 5, 5, -5, 7, 1, 1, 7}},
                          // x
                          {{-5, -5, 5, 5, 0, 0, 0, 0}},
                          // y
                          {{1, 1, 1, 1, -5, -5, 5, 5}},
                          // z
                          {{1, 1, 1, 1, 1, 1, 1, 1}},
                          // nx
                          {{1, 1, 1, 1, 1, 1, 1, 1}},
                          // ny
                          {{1, 1, 1, 1, 1, 1, 1, 1}} // nz
        );
        universe.addBonds(8,
                          {{1, 2, 3, 4, 5, 6, 7, 8}},
                          {{2, 3, 4, 1, 6, 7, 8, 5}},
                          {{1, 1, 1, 1, 1, 1, 1, 1}},
                          false,
                          false);
        REQUIRE(universe.getNrOfBonds() == 8);
        const std::vector<long int> empty;
        CHECK(universe.findLoopEntanglements(empty, empty, empty, empty).size() ==
            0);
        CHECK(
            universe
            .findLoopEntanglements({ { 1, 2, 7 } }, { { 0, 3, 6 } }, empty, empty)
            .size() == 0);
        CHECK(
            universe
            .findLoopEntanglements({ { 0, 1, 2 } }, { { 5, 7, 6 } }, empty, empty)
            .size() == 2);
        CHECK(universe
            .findLoopEntanglements(
                { { 0, 1, 2, 3 } }, { { 4, 5, 6, 7 } }, empty, empty)
            .size() > 0);
        CHECK(universe
            .findLoopEntanglements(
                { { 4, 5, 6, 7 } }, { { 0, 1, 2, 3 } }, empty, empty)
            .size() > 0);
    }

    SECTION("Molecules with crossLinkers are found") {
        universe.setBox(pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0));
        /**
        # The system looks like this (in terms of bonds, not 3D placement):
        # 1-2-3-*6
        # |      |
        # *7-5---|
        # 8
        #
        # *4
        */
        universe.addAtoms({{1, 2, 3, 4, 5, 6, 7, 8}},
                          // id
                          {{1, 1, 1, 2, 1, 2, 2, 1}},
                          // type
                          {{1.25, 2, 3, 1.01, 2, 4, 1, 1}},
                          // x
                          {{1, 1, 1, 4.01, 2, 1, 2, 3}},
                          // y
                          {{1, 1, 1, 1.01, 1, 1, 1, 1}},
                          // z
                          {{1, 1, 1, 1, 1, 1, 1, 1}},
                          // nx
                          {{1, 1, 1, 1, 1, 1, 1, 1}},
                          // ny
                          {{1, 1, 1, 1, 1, 1, 1, 1}} // nz
        );
        universe.addBonds(7,
                          {{1, 2, 3, 6, 5, 7, 7}},
                          {{2, 3, 6, 5, 7, 1, 8}},
                          {{1, 1, 1, 1, 1, 1, 11}});

        SECTION("Graph copy can be accessed") {
            igraph_t graph = universe.getCopyOfGraph();
            CHECK(igraph_vcount(&graph) == universe.getNrOfAtoms());
            CHECK(igraph_ecount(&graph) == universe.getNrOfBonds());
            igraph_integer_t degree;
            igraph_degree_1(&graph, &degree, 0, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);
            CHECK(degree == 2);
        }

        SECTION("Atoms can be replaced") {
            REQUIRE_THROWS(universe.replaceAtom(2, universe.getAtom(3)));
            pe::Atom replacementAtom = pe::Atom(2, 2, 1.23, 1.23, 1.23, 0, 0, 0);
            universe.replaceAtom(2, replacementAtom);
            pe::Atom replacedAtom = universe.getAtom(2);
            REQUIRE(replacedAtom == replacementAtom);
        }

        SECTION("Atom types are counted") {
            std::map<int, int> atomCounts = universe.countAtomTypes();
            REQUIRE(atomCounts.size() == 2);
            REQUIRE(atomCounts[1] == 5);
            REQUIRE(atomCounts[2] == 3);
        }

        SECTION("Shortest paths are found") {
            std::vector<pe::Atom> path = universe.getShortestPath(0, 2);
            REQUIRE(path.size() == 3);
        }

        SECTION("Masses are persisted in session") {
            std::map<int, double> weightFractions = universe.computeWeightFractions();
            CHECK(!weightFractions.empty());
            CHECK(weightFractions[2] ==
                (3. / 8.)); // without masses, the default is 1.0 per type
            std::map<int, double> masses = universe.getMasses();
            REQUIRE(masses.size() == 0);
            masses[1] = 1.0;
            masses[2] = 2.0;
            universe.setMasses(masses);
            std::map<int, double> newMasses = universe.getMasses();
            CHECK(newMasses[1] == 1.0);
            CHECK(newMasses[2] == 2.0);
            weightFractions = universe.computeWeightFractions();
            CHECK(weightFractions[1] == (5.0 * 1.0) / (3.0 * 2.0 + 5.0 * 1.0));
            CHECK(weightFractions[2] == (2.0 * 3.0) / (3.0 * 2.0 + 5.0 * 1.0));
            CHECK(universe.computeTotalMass() == Catch::Approx(5 * 1.0 + 3 * 2.0));
            std::map<int, double> otherMasses(masses);
            otherMasses[2] = 0.0;
            CHECK(otherMasses[2] == 0.0);
            CHECK(universe.computeTotalMassWithMasses(otherMasses) ==
                Catch::Approx(5 * 1.0));
        }

        SECTION("Bonds can be removed") {
            int bondsBefore = universe.getNrOfBonds();
            REQUIRE(bondsBefore == 7);
            universe.writeGraphToFile("test-graph.gml");
            REQUIRE_THROWS(universe.removeBonds({ {} }, { { 1, 2, 3 } }));
            universe.removeBonds({{1, 2, 3}}, {{2, 3, 6}});
            CHECK(universe.getNrOfBonds() == bondsBefore - 3);
        }

        SECTION("Get bonds returns") {
            auto edges = universe.getEdges();
            CHECK(edges["edge_from"][0] == 1 - 1);
            CHECK(edges["edge_to"][0] == 2 - 1);
            auto bonds = universe.getBonds();
            CHECK(bonds.at("bond_from")[0] == 1);
            CHECK(bonds.at("bond_to")[0] == 2);
            CHECK(bonds["bond_from"].size() == bonds["bond_to"].size());
            CHECK(bonds["bond_from"].size() == edges["edge_to"].size());
            CHECK(bonds["bond_from"].size() == edges["edge_from"].size());
            CHECK(bonds["bond_from"].size() == 7);
            CHECK(bonds["bond_from"][3] ==
                universe.getAtomIdByIdx(edges["edge_from"][3]));
            CHECK(bonds["bond_to"][3] ==
                universe.getAtomIdByIdx(edges["edge_to"][3]));
            CHECK(universe.getIdxByAtomId(bonds["bond_from"][2]) ==
                edges["edge_from"][2]);
            CHECK(universe.getIdxByAtomId(bonds["bond_to"][2]) ==
                edges["edge_to"][2]);
            CHECK(bonds["bond_type"][5] == 1);
            CHECK(bonds["bond_type"][6] == 11);
            // get atoms with type returns
            CHECK(universe.getAtomsOfType(2).size() == 3);
            CHECK(universe.getAtomsOfType(1).size() == 5);
            CHECK(universe.getAtomsOfType(0).size() == 0);
        }

        SECTION("internal lengths are computed correctly") {
            REQUIRE(universe.computeMeanBondLength() == Catch::Approx(1.1452634834));
            auto chains = universe.getChainsWithCrosslinker(2);
            REQUIRE(chains.size() == 3);
            REQUIRE(chains[0].computeEndToEndDistance() == Catch::Approx(sqrt(10.0)));
            REQUIRE(
                universe.computeMeanEndToEndDistance(2) ==
                Catch::Approx((2.0 / 3.0) * sqrt(9.0 + 1.0) + (1.0 / 3.0) * sqrt(1.0)));
            REQUIRE(universe.computeDxs({ { 7, 6 } }, { { 6, 7 } }).size() == 2);
            std::vector<double> zeros;
            zeros.push_back(0.0);
            zeros.push_back(0.0);
            REQUIRE(universe.computeDzs({ { 7, 6 } }, { { 6, 7 } }) == zeros);
            std::vector<double> dys;
            dys.push_back(-1.0);
            dys.push_back(1.0);
            REQUIRE(universe.computeDys({ { 7, 6 } }, { { 6, 7 } }) == dys);
        }

        SECTION("get angles returns") {
            REQUIRE(universe.getAngles()["angle_from"].size() == 0);
            auto detectedAngles = universe.detectAngles();
            REQUIRE(detectedAngles["angle_from"].size() == 8);
            std::vector<int> angleTypes;
            angleTypes.reserve(detectedAngles["angle_from"].size());
            for (size_t i = 0; i < detectedAngles["angle_from"].size(); i++) {
                angleTypes.push_back(1);
            }
            universe.addAngles(detectedAngles["angle_from"],
                               detectedAngles["angle_via"],
                               detectedAngles["angle_to"],
                               angleTypes);
            CHECK(detectedAngles["angle_from"][0] == 2);
            CHECK(detectedAngles["angle_to"][0] == 7);
            CHECK(detectedAngles["angle_via"][0] == 1);
            REQUIRE(universe.getAngles()["angle_from"].size() == 8);

            // angles can be computed correctly
            std::vector<double> angles = universe.computeAngles();
            REQUIRE(angles.size() == 8);
            CHECK(Catch::Approx(angles[0]) == 1.815770741);
        }

        SECTION("get dihedral angles returns") {
            REQUIRE(universe.getDihedralAngles()["dihedral_angle_from"].size() == 0);
            auto detectedAngles = universe.detectDihedralAngles();
            REQUIRE(detectedAngles["dihedral_angle_from"].size() == 8);
            std::vector<int> angleTypes;
            angleTypes.reserve(detectedAngles["dihedral_angle_from"].size());
            for (size_t i = 0; i < detectedAngles["dihedral_angle_from"].size();
                 i++) {
                angleTypes.push_back(1);
            }
            universe.addDihedralAngles(detectedAngles["dihedral_angle_from"],
                                       detectedAngles["dihedral_angle_via1"],
                                       detectedAngles["dihedral_angle_via2"],
                                       detectedAngles["dihedral_angle_to"],
                                       angleTypes);
            REQUIRE(universe.getDihedralAngles()["dihedral_angle_from"].size() == 8);
        }

        SECTION("get atoms returns") {
            // get atoms with type returns atoms with properties
            std::vector<pe::Molecule> molecules = universe.getMolecules(2);
            REQUIRE(molecules[0].getAtomsOfType(1).size() == 3);
            REQUIRE(molecules[1].getAtomsOfType(1)[0].getId() == 5);
            REQUIRE(molecules[2].getAtomsOfType(1)[0].getType() == 1);
            // get molecules allows to fetch atoms with degree
            REQUIRE(molecules[0].getAtomsOfDegree(2).size() == 1);
            REQUIRE(molecules[0].getAtomsOfDegree(1).size() == 2);
            REQUIRE(molecules[0].getAtomsOfDegree(0).size() == 0);
            //
            REQUIRE_THROWS(universe.getAtomByVertexIdx(999));
        }

        SECTION("get atoms with crossLinkers returns") {
            // get atoms with crossLinkers returns
            auto chains = universe.getChainsWithCrosslinker(2);
            REQUIRE(chains.size() == 3);
            REQUIRE(chains[0].getAtoms()[0].getId() == 1);
            REQUIRE(chains[0].getAtoms()[0].getX() == 1.25);
            REQUIRE(chains[0].getNrOfAtoms() == 5);
            REQUIRE(chains[0].getKey() == "1-2-3-6-7");
            REQUIRE(chains[0].getAtomsOfType(2).size() == 2);
            //
            auto functionalityPerType = universe.determineFunctionalityPerType();
            REQUIRE(functionalityPerType[1] == 2);
            REQUIRE(functionalityPerType[2] == 3);
            REQUIRE(functionalityPerType.size() == 2);
            //
            REQUIRE(chains[0].getAtomsOfDegree(2).size() == 3);
            REQUIRE(chains[0].getAtomsOfDegree(1).size() == 2);
            REQUIRE(chains[0].getAtomsOfDegree(0).size() == 0);
        }

        SECTION("Loops are found") {
            std::map<int, std::vector<std::vector<pe::Atom> > > loops =
                    universe.findLoopsOfAtoms(2, -1, true);
            REQUIRE(pylimer_tools::utils::map_has_key(loops, 2));
            REQUIRE(loops.size() == 1);

            auto allBonds = universe.getBonds();
            std::vector<long int> bondsOf1;
            std::vector<long int> bondsOf2;
            for (int i = 0; i < allBonds["bond_from"].size(); ++i) {
                if (allBonds["bond_from"][i] == 1) {
                    bondsOf1.push_back(allBonds["bond_to"][i]);
                }
                if (allBonds["bond_from"][i] == 2) {
                    bondsOf2.push_back(allBonds["bond_to"][i]);
                }
            }
            REQUIRE(bondsOf1.size() == 2);
            REQUIRE_THAT(
                bondsOf1,
                Catch::Matchers::UnorderedEquals(std::vector<long int>{ 2, 7 }));
            REQUIRE(bondsOf2.size() == 1);
            REQUIRE(bondsOf2[0] == 3);
            REQUIRE(universe
                .getEdgeIdsFromTo(universe.getIdxByAtomId(1),
                    universe.getIdxByAtomId(2))
                .size() == 1);
            std::vector<pe::Atom> minimalLoop1 =
                    universe.findMinimalOrderLoopFrom(1, 2, -1, false);
            REQUIRE(minimalLoop1.size() == 6);
            // translate to ids for simpler matching
            std::vector<int> minimalLoopIds;
            minimalLoopIds.reserve(minimalLoop1.size());
            for (pe::Atom a: minimalLoop1) {
                minimalLoopIds.push_back(a.getId());
            }
            REQUIRE(*std::max_element(minimalLoopIds.begin(), minimalLoopIds.end()) ==
                7);
            REQUIRE(*std::min_element(minimalLoopIds.begin(), minimalLoopIds.end()) ==
                1);
            REQUIRE_THAT(minimalLoopIds, Catch::Matchers::VectorContains(2));
        }

        SECTION("Clusters are found") {
            std::vector<pe::Universe> clusters = universe.getClusters();
            REQUIRE(clusters.size() == 2);
            REQUIRE(clusters[1].getAtoms()[0].getId() == 4);
            REQUIRE_THROWS(clusters[0].getMolecules(-1)[0].getAtomsLinedUp());
        }

        SECTION("Loops are found in reduced universe") {
            /**
             * @brief increase system size, reduce to x-linkers, test
             *
             * The system looks like this (in terms of bonds, not 3D placement):
             *
             * 1-2-3-*6
             * |      |
             * *7-5---|
             * 8
             *
             * *4=9
             *
             * 10
             */
            universe.addAtoms({{9, 10}},
                              {{1, 1}},
                              {{0.0, 0.0}},
                              {{1.0, 0.0}},
                              {{1.0, 0.0}},
                              {{0, 1}},
                              {{0, 1}},
                              {{0, 1}});
            universe.addBonds(
                2, {{4, 9}}, {{9, 4}}, {{1, 1}}, false, false);
            REQUIRE(universe.getNrOfBonds() == 2 + 7);
            pe::Universe reducedUniverse = universe.getNetworkOfCrosslinker(2);
            REQUIRE(reducedUniverse.getNrOfAtoms() == 3);
            std::map<int, std::vector<std::vector<pe::Atom> > > loops =
                    reducedUniverse.findLoopsOfAtoms(2, -1);
            REQUIRE(loops.size() == 2);
            REQUIRE(reducedUniverse.getPropertyValues<int>("id").size() == 3);
            REQUIRE(reducedUniverse.getAtomTypes().size() == 3);
            for (const auto &myPair: loops) {
                for (const std::vector<pe::Atom> &loop: myPair.second) {
                    REQUIRE(loop.size() == myPair.first);
                }
            }
            REQUIRE(loops[1][0][0].getId() == 4);
            REQUIRE(loops[1].size() == 1);
            REQUIRE(loops[2].size() == 1);
            REQUIRE_THAT(loops[2][0],
                         Catch::Matchers::UnorderedEquals(std::vector<pe::Atom>{
                             reducedUniverse.getAtom(6), reducedUniverse.getAtom(7) }));
            std::map<std::string, std::vector<long int> > bonds =
                    reducedUniverse.getBonds();
            REQUIRE_THAT(
                bonds["bond_from"],
                Catch::Matchers::UnorderedEquals(std::vector<long int>{ 4, 6, 6 }));
            REQUIRE_THAT(
                bonds["bond_to"],
                Catch::Matchers::UnorderedEquals(std::vector<long int>{ 4, 7, 7 }));
            REQUIRE(reducedUniverse.getNrOfBonds() == 3);

            // second-order loop (in reduced universe)
            REQUIRE(reducedUniverse
                .getEdgeIdsFromTo(reducedUniverse.getIdxByAtomId(6),
                    reducedUniverse.getIdxByAtomId(7))
                .size() == 2);
            std::vector<pe::Atom> minimalLoop2 =
                    reducedUniverse.findMinimalOrderLoopFrom(6, 7, -1, false);
            REQUIRE(minimalLoop2.size() == 2);
            std::vector<pe::Atom> minimalLoop4 =
                    universe.findMinimalOrderLoopFrom(6, 5, -1, false);
            REQUIRE(minimalLoop4.size() == 6);
            // self-loop
            std::vector<pe::Atom> minimalLoop3 =
                    reducedUniverse.findMinimalOrderLoopFrom(4, 4, -1, false);
            REQUIRE(minimalLoop3.size() == 1);

            // second-order loop
            REQUIRE(universe
                .getEdgeIdsFromTo(universe.getIdxByAtomId(4),
                    universe.getIdxByAtomId(9))
                .size() == 2);
            std::vector<pe::Atom> minimalLoop1 =
                    universe.findMinimalOrderLoopFrom(4, 9, -1, false);
            REQUIRE(minimalLoop1.size() == 2);
            REQUIRE_THROWS(reducedUniverse.findMinimalOrderLoopFrom(4, 9, -1, false));

            std::vector<pe::Atom> minimalLoop5 =
                    reducedUniverse.findMinimalOrderLoopFrom(6, 7, -1, false);
            REQUIRE(minimalLoop5.size() == 2);
        }

        SECTION("Infinite Strands are found") {
            universe.setBoxLengths(4.0, 4.0, 2.0);
            REQUIRE(universe.hasInfiniteStrand(2, -1) == false);
            universe.addBonds(2, {{8, 4}}, {{4, 1}});
            std::map<int, std::vector<std::vector<pe::Atom> > > loops =
                    universe.findLoopsOfAtoms(2, -1);
            CHECK(loops.size() == 2);
            CHECK(universe.hasInfiniteStrand(2, -1) == true);
            CHECK(universe.getMeanStrandLength(2) ==
                Catch::Approx(((double)(3 + 1 + 1)) / 3.0));
        }

        SECTION("Reduction to Cross-linker-verse works") {
            pe::Universe reducedUniverse = universe.getNetworkOfCrosslinker(2);
            REQUIRE(reducedUniverse.getNrOfAtoms() == 3);
            REQUIRE(reducedUniverse.getAtomsOfType(2).size() == 3);

            auto edges = reducedUniverse.getEdges();

            REQUIRE(reducedUniverse.getNrOfBonds() == 2);
            REQUIRE(edges.at("edge_from").size() == 2);
            REQUIRE(edges.at("edge_to").size() == 2);
            REQUIRE(edges.at("edge_from")[0] == 1);
            REQUIRE(edges.at("edge_to")[0] == 2);
            REQUIRE(edges.at("edge_from")[1] == 1);
            REQUIRE(edges.at("edge_to")[1] == 2);

            auto atom = reducedUniverse.getAtom(4);
            REQUIRE(atom.getX() == 1.01);
            REQUIRE(atom.getY() == 4.01);
            REQUIRE(atom.getZ() == 1.01);

            SECTION("Copy constructors/assignment work") {
                pe::Universe newUniverse = pe::Universe(universe);
                newUniverse = reducedUniverse;
                REQUIRE(newUniverse.getNrOfBonds() == 2);
                REQUIRE(reducedUniverse.getNrOfAtoms() == 3);
                REQUIRE(newUniverse.getNrOfAtoms() == 3);
            }
        }

        SECTION("Computations work") {
            universe.setBoxLengths(10.0, 10.0, 10.0);

            // To understand where the following computations are coming from
            auto chains = universe.getChainsWithCrosslinker(2);
            REQUIRE(chains.size() == 3);
            REQUIRE(chains[0].getKey() == "1-2-3-6-7");
            REQUIRE(chains[1].getKey() == "5-6-7");
            REQUIRE(chains[2].getKey() == "7-8");

            auto endToEndDistances = universe.computeEndToEndDistances(2);
            CHECK(endToEndDistances.size() == 3);
            CHECK(endToEndDistances[0] == chains[0].computeEndToEndDistance());
            CHECK(endToEndDistances[0] == Catch::Approx(sqrt(3.0 * 3.0 + 1.0 * 1.0)));
            CHECK(endToEndDistances[1] == chains[1].computeEndToEndDistance());
            CHECK(endToEndDistances[1] == Catch::Approx(sqrt(3.0 * 3.0 + 1.0 * 1.0)));
            CHECK(endToEndDistances[2] == chains[2].computeEndToEndDistance());
            CHECK(endToEndDistances[2] == 1.0);
            CHECK(universe.computeMeanEndToEndDistance(2) ==
                Catch::Approx((2.0 * sqrt(3.0 * 3.0 + 1.0 * 1.0) + 1.0) / 3.0));
            CHECK(universe.computeMeanSquareEndToEndDistance(2) ==
                Catch::Approx((2.0 * (3.0 * 3.0 + 1.0 * 1.0) + 1.0) / 3.0));
        }
    }

    SECTION("Molecule Types are determined correctly") {
        /**
        # The system looks like this (in terms of bonds, not 3D placement):
        # 1-2-3
        #
        # *7-*6-5
        */
        universe.addAtoms({{1, 2, 3, 5, 6, 7}},
                          // id
                          {{1, 1, 1, 1, 2, 2}},
                          // type
                          {{3, 2, 2, 2, 2, 2}},
                          // x
                          {{1, 1, 1, 1, 1, 1}},
                          // y
                          {{1, 1, 1, 1, 1, 1}},
                          // z
                          {{1, 1, 1, 1, 1, 1}},
                          // nx
                          {{1, 1, 1, 1, 1, 1}},
                          // ny
                          {{1, 1, 1, 1, 1, 1}} // nz
        );
        universe.addBonds(4, {{1, 2, 5, 6}}, {{2, 3, 6, 7}});
        CHECK(universe.getNrOfBonds() == 4);
        CHECK(universe.getMolecules(2).size() == 2);
        CHECK(universe.getChainsWithCrosslinker(2).size() == 3);
        CHECK(universe.computeFunctionalityForVertex(universe.getIdxByAtomId(1)) ==
            1);
        CHECK(universe.computeFunctionalityForAtom(1) == 1);
        CHECK(universe.getAtomsConnectedTo(universe.getIdxByAtomId(1)).size() == 1);
        auto chains = universe.getChainsWithCrosslinker(2);
        CHECK(chains.size() == 3);
        CHECK(chains[0].getNrOfAtoms() == 3);
        CHECK(chains[0].getType() == pe::MoleculeType::FREE_CHAIN);
        CHECK(chains[1].getNrOfAtoms() == 2);
        CHECK(chains[1].getType() == pe::MoleculeType::DANGLING_CHAIN);
        CHECK(chains[2].getNrOfAtoms() == 2);
        // other checks
        CHECK(chains[0].getKey() == "1-2-3");
        CHECK(chains[1].getKey() == "5-6");
        CHECK(chains[2].getKey() == "6-7");
        universe.setBoxLengths(3.0, 3.0, 3.0);
        CHECK(universe.determineEffectiveFunctionalityPerType()[2] ==
            (1. + 2.) / 2.);
        CHECK(chains[0].computeEndToEndDistance() == 1.0);
        CHECK(chains[1].computeEndToEndDistance() == 0.0);
        CHECK(universe.computeMeanEndToEndDistance(2) == 1.0);
        CHECK(universe.computeMeanSquareEndToEndDistance(2) == 1.0);

        SECTION("Primary loops too") {
            // instead of rewriting above tests,
            // let's just add new atoms & bonds.
            // the resulting new structure will look like this (in terms of bonds, not
            // 3D placement):
            // 1-2-3
            //
            // 9-8
            // | |
            // ↳-*7-*6-5
            universe.addAtoms({{8, 9}},
                              // id
                              {{1, 1}},
                              // type
                              {{0, 1}},
                              // x
                              {{2, 3}},
                              // y
                              {{0, 0}},
                              // z
                              {{1, 1}},
                              // nx
                              {{1, 1}},
                              // ny
                              {{1, 1}} // nz
            );
            universe.addBonds(3, {{7, 8, 9}}, {{8, 9, 7}});
            CHECK(universe.getNrOfBonds() == 7);
            auto newChains = universe.getChainsWithCrosslinker(2);
            CHECK(newChains.size() == 4);
            CHECK(newChains[2].getType() == pe::MoleculeType::PRIMARY_LOOP);
            CHECK(newChains[2].getKey() == "7-8-9");
        }
    }

    SECTION("PolyDispersity Calculation") {
        SECTION("Literature system") {
            // as taken from https://www.pslc.ws/macrog/average.htm
            std::vector<int> nrOfMolecules = {
                {1, 3, 5, 8, 10, 13, 20, 13, 10, 8, 5, 3, 1}
            };
            std::vector<double> massOfMolecules = {
                {8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3., 2.5, 2}
            };

            std::map<int, double> weights;
            int currentId = 0;

            std::vector<double> oneZeroDouble = getVectorWithOne<double>(0.0);
            std::vector<int> oneZeroInt = getVectorWithOne(0);

            // first, add these molecules as single atoms to the universe
            for (size_t i = 0; i < nrOfMolecules.size(); ++i) {
                std::vector<int> type;
                type.push_back(i);
                // this is slow and could easily be optimized, but whatever
                for (size_t j = 0; j < nrOfMolecules[i]; ++j) {
                    std::vector<long int> atomId;
                    atomId.push_back(currentId);
                    universe.addAtoms(atomId,
                                      type,
                                      oneZeroDouble,
                                      oneZeroDouble,
                                      oneZeroDouble,
                                      oneZeroInt,
                                      oneZeroInt,
                                      oneZeroInt);
                    currentId += 1;
                }
                weights.emplace(i, massOfMolecules[i] * 100000.);
            }
            universe.setMasses(weights);
            // some asserts that everything has been correcly set by the test
            CHECK(universe.getMolecules(-1).size() == 100);
            CHECK(currentId == 100);
            // then, do the calculation
            CHECK(universe.computeTotalMass() == Catch::Approx(50000000));
            CHECK(universe.computeNumberAverageMolecularWeight(-1) ==
                Catch::Approx(500000));
            CHECK(universe.computeWeightAverageMolecularWeight(-1) ==
                Catch::Approx(531600));
            CHECK(universe.computePolydispersityIndex(-1) == Catch::Approx(1.0632));
        }

        SECTION("System to zero out") {
            // other system
            std::vector<double> oneZeroDouble = getVectorWithOne<double>(0.0);
            std::vector<int> oneZeroInt = getVectorWithOne(0);

            pe::Universe universe2 = pe::Universe(10., 10., 10.);
            std::map<int, double> weights2;
            weights2.emplace(1, 1.0);
            weights2.emplace(2, 3.0);
            universe2.setMasses(weights2);
            std::vector<int> type;
            type.push_back(1);
            long int currentId = 0;
            for (size_t i = 0; i < 10; ++i) {
                for (size_t j = 0; j < 10; ++j) {
                    currentId += 1;
                    universe2.addAtoms(getVectorWithOne<long int>(currentId),
                                       type,
                                       oneZeroDouble,
                                       oneZeroDouble,
                                       oneZeroDouble,
                                       oneZeroInt,
                                       oneZeroInt,
                                       oneZeroInt);
                    if (j > 0) {
                        universe2.addBonds(getVectorWithOne<long int>(currentId),
                                           getVectorWithOne<long int>(currentId - 1));
                    }
                }
                CHECK(universe2.computePolydispersityIndex(-1) == Catch::Approx(1.0));
            }
            // add cross-linkers
            SECTION("With cross-linkers") {
                std::vector<int> type2 = getVectorWithOne(2);
                for (size_t i = 0; i < 10; ++i) {
                    currentId += 1;
                    universe2.addAtoms(getVectorWithOne<long int>(currentId),
                                       type2,
                                       oneZeroDouble,
                                       oneZeroDouble,
                                       oneZeroDouble,
                                       oneZeroInt,
                                       oneZeroInt,
                                       oneZeroInt);
                }
                CHECK(universe2.computePolydispersityIndex(2) == Catch::Approx(1.0));
            }
        }
    }

    SECTION("Local Density Computation") {
        const pe::Box box = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
        universe.setBox(box);
        universe.addAtoms({{1, 2, 3, 4, 5, 6, 7, 8}},
                          // id
                          {{2, 1, 1, 1, 2, 1, 1, 1}},
                          // type
                          {{1., 2., 3., 4., 9., -10., -9., -8.}},
                          // x
                          {{1., 2., 3., 4., 9., -10., -9., -8.}},
                          // y
                          {{1., 2., 3., 4., 9., -10., -9., -8.}},
                          // z
                          {{0, 0, 0, 0, 1, 2, 2, 2}},
                          // nx
                          {{0, 0, 0, 0, 1, 2, 2, 2}},
                          // ny
                          {{0, 0, 0, 0, 1, 2, 2, 2}} // nz
        );

        std::vector<double> distances = {0., 1., 2., 3.};
        std::vector<size_t> result = universe.countAtomsInSkinDistance(distances);
        CHECK(result.size() == distances.size() - 1);
        CHECK(result.size() == 3);
        CHECK(result[0] == 0);
        // num neighbours between 1. and 2. (sqrt(3))
        CHECK(result[1] == 12);
        CHECK(result[2] == 0);
    }
}

TEST_CASE("Universe can be decomposed", "[Universe][entity]") {
    std::cout << "Running test \"Universe can be decomposed\"" << std::endl;
    size_t nrOfBeads = 30;
    size_t nrOfBeadsPerChain = 3;
    pe::Universe universe =
            pe::Universe(nrOfBeads * 10.0, nrOfBeads * 10.0, nrOfBeads * 10.0);
    std::vector<double> xPositions, yPositions, zPositions;
    std::vector<int> atomTypes, zeroInts;
    std::vector<long int> atomIds;
    std::vector<long int> bondFrom, bondTo;
    double offset = 10.0;
    for (int i = 0; i < nrOfBeads; ++i) {
        xPositions.push_back(i * 1.0 + offset);
        yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                             offset); // /!\ i needs to be int, not unsigned!
        zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
        atomIds.push_back(i);
        atomTypes.push_back(i % nrOfBeadsPerChain == 0 ? 2 : 1);
        zeroInts.push_back(0);
        if (i > 0) {
            bondFrom.push_back(i - 1);
            bondTo.push_back(i);
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
    REQUIRE(universe.getNrOfAtoms() == nrOfBeads);
    REQUIRE(universe.getNrOfBonds() == nrOfBeads - 1);
    //

    CHECK(universe.getChainsWithCrosslinker(2).size() == 10);

    SECTION("Also with secondary loop") {
        // add a secondary loop to check that it is returned as well
        universe.addBonds({0}, {static_cast<long>(nrOfBeadsPerChain)});
        CHECK(universe.getChainsWithCrosslinker(2).size() == 11);
    }

    SECTION("Also with primary loop") {
        // add a primary loop to check that it is returned as well
        universe.addBonds({0}, {0});
        CHECK(universe.getChainsWithCrosslinker(2).size() == 11);
    }

    SECTION("Also with more loops") {
        // add a primary loop to check that it is returned as well
        universe.addBonds({0, 0}, {static_cast<long>(nrOfBeadsPerChain), 0});
        CHECK(universe.getChainsWithCrosslinker(2).size() == 12);
    }
}

TEST_CASE("Coordinates work") {
    pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);

    SECTION("with molecules and lined up") {
        /**
         * @brief A grid of two rows, each one bead between the two crosslinkers
         *
         */
        universe.addAtoms(
            {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}},
            {{2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1}},
            {{0., 2.5, 5, 7.5, 0.1, 2.5, 5, 7.5, -0.1, 5., 0., 5.}},
            // x with slight (0.1) deviation, so we don't start perfect
            {{0.1, 0., -0.1, 0., 5., 5., 5., 5., 2.5, 2.5, 7.5, 7.5}},
            // y
            {{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}},
            // z
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
            {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}});
        universe.addBonds(
            {{1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7}},
            {{2, 9, 4, 11, 2, 4, 10, 12, 9, 11, 6, 8, 6, 8, 10, 12}});
        std::vector<pe::Molecule> molecules = universe.getChainsWithCrosslinker(2);

        for (const pe::Molecule &molecule: molecules) {
            CHECK(molecule.getLength() == 3);
            std::vector<pe::Atom> atoms = molecule.getAtomsLinedUp();
            std::vector<long int> alignedVertices = molecule.getVerticesLinedUp(2);
            CHECK(alignedVertices.size() == molecule.getLength());
            CHECK(atoms.size() == molecule.getLength());
            for (size_t i = 0; i < 3; ++i) {
                CHECK(molecule.getIdxByAtomId(atoms[i].getId()) == alignedVertices[i]);
            }
            CHECK(atoms[0].getType() == 2);
            CHECK(atoms[1].getType() == 1);
            CHECK(atoms[2].getType() == 2);
            Eigen::Vector3d tmpDist = molecule.getOverallBondSum(2);
            CHECK(tmpDist[2] == 0.0);
        }

        std::vector<pe::Atom> atomsLinedUp = molecules[0].getAtomsLinedUp();
        CHECK(atomsLinedUp[0].getId() == 1);
        CHECK(atomsLinedUp[1].getId() == 2);
        CHECK(atomsLinedUp[2].getId() == 3);
        Eigen::Vector3d dist = molecules[0].getOverallBondSum(2);
        CHECK(dist[0] == 5.);
        CHECK(dist[1] == -0.2);

        std::vector<long int> alignedVertices = molecules[0].getVerticesLinedUp(2);
        Eigen::VectorXd alignedCoordinates =
                Eigen::VectorXd::Zero(3 * alignedVertices.size());
        pe::Box box = universe.getBox();
        Eigen::VectorXd vertexCoordinates =
                molecules[0].getUnwrappedVertexCoordinates(alignedVertices, box);
        molecules[0].getAssumedVertexCoordinates(
            alignedCoordinates, box, alignedVertices);
        CHECK(vertexCoordinates[0] == 0.);
        CHECK(vertexCoordinates[1] == 0.1);
        CHECK(alignedCoordinates[0] == 0.);
        CHECK(alignedCoordinates[1] == 0.1);
        for (size_t i = 0; i < 3; ++i) {
            CHECK(alignedCoordinates[3 * i + 0] == atomsLinedUp[i].getX());
            CHECK(alignedCoordinates[3 * i + 1] == atomsLinedUp[i].getY());
            CHECK(alignedCoordinates[3 * i + 2] == atomsLinedUp[i].getZ());
            //
            CHECK(alignedCoordinates[3 * i + 0] == vertexCoordinates[3 * i + 0]);
            CHECK(alignedCoordinates[3 * i + 1] == vertexCoordinates[3 * i + 1]);
            CHECK(alignedCoordinates[3 * i + 2] == vertexCoordinates[3 * i + 2]);
        }
    }

    SECTION("Coordinates are fetched appropriately") {
        const pe::Box box = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
        universe.setBox(box);
        universe.addAtoms({{1, 2, 3, 4, 5, 6, 7, 8}},
                          // id
                          {{2, 1, 1, 1, 2, 1, 1, 1}},
                          // type
                          {{1, 2, 3, 4, 9, -10, -9, -8}},
                          // x
                          {{1, 2, 3, 4, 9, -10, -9, -8}},
                          // y
                          {{1, 2, 3, 4, 9, -10, -9, -8}},
                          // z
                          {{0, 0, 0, 0, 1, 2, 2, 2}},
                          // nx
                          {{0, 0, 0, 0, 1, 2, 2, 2}},
                          // ny
                          {{0, 0, 0, 0, 1, 2, 2, 2}} // nz
        );
        std::vector<long int> indices = {{1, 2, 3}};
        Eigen::VectorXd coordinates =
                universe.getUnwrappedVertexCoordinates(indices, box);
        CHECK(coordinates.size() == 9);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t dir = 0; dir < 3; ++dir) {
                CHECK(coordinates[i * 3 + dir] == i + 2);
            }
        }
    }
}

TEST_CASE("Large universe can be used", "[Universe][entity]") {
    pe::UniverseSequence universeSeq = pe::UniverseSequence();
    CHECK(universeSeq.getLength() == 0);
    std::string suspectedPath = "../pylimer_tools/fixtures/";
    universeSeq.initializeFromDataSequence(
        {{suspectedPath + "lammps_data_file.out"}});
    CHECK(universeSeq.getLength() == 1);

    pe::Universe universe = universeSeq.atIndex(0);

    SECTION("Edges can be interpolated") {
        std::map<std::string, std::vector<long int> > edges = universe.getEdges();
        std::vector<std::pair<size_t, size_t> > interpolatedEdges1 =
                universe.interpolateEdges(2, 1);
        CHECK(interpolatedEdges1.size() == edges["edge_from"].size());
        std::vector<std::pair<size_t, size_t> > interpolatedEdges2 =
                universe.interpolateEdges(2, 2.);
        CHECK(interpolatedEdges2.size() ==
            Catch::Approx(2. * edges["edge_from"].size()));
    }
}

TEST_CASE("Vertex coordinates are assumed for tree-like structures", "[Universe][entity]") {
    std::cout << "Running test \"Vertex coordinates are assumed for tree-like structures\"" << std::endl;

    pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
    std::vector<double> coords = {
        0.0, 1.0, 10.0, 21.0, 22.0, 3.0, 2.0, 3.0
    };
    std::vector<long int> ids = {
        0, 1, 2, 3, 4, 5, 6, 7
    };
    universe.addAtoms(
        ids, pylimer_tools::utils::initializeWithValue(ids.size(), 1),
        coords, coords, coords,
        pylimer_tools::utils::initializeWithValue(ids.size(), 0),
        pylimer_tools::utils::initializeWithValue(ids.size(), 0),
        pylimer_tools::utils::initializeWithValue(ids.size(), 0)
    );
    std::vector<long int> bondFrom = {0, 1, 2, 3, 3, 4, 6};
    std::vector<long int> bondTo = {1, 2, 3, 4, 6, 5, 7};
    universe.addBonds(
        bondFrom,
        bondTo
    );

    CHECK(universe.getNrOfBonds() == bondFrom.size());
    CHECK(universe.getNrOfAtoms() == ids.size());

    Eigen::VectorXd assumedCoordinates = universe.getAssumedVertexCoordinates(universe.getBox());
    // for each bond, check that the bond-length is shorter than the half box, even without PBC
    for (size_t i = 0; i < bondFrom.size(); ++i) {
        Eigen::Vector3d bondVector = (
            assumedCoordinates.segment(3 * bondTo[i], 3) - assumedCoordinates.segment(3 * bondFrom[i], 3)
        );
        CHECK(bondVector.norm() < 5.);
    }
}
