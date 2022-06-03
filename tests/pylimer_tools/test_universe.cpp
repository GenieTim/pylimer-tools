#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
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

TEST_CASE("Universe can be used", "[entity][Universe]")
{
  REQUIRE(1 == 1);

  pe::Universe universe = pe::Universe(1.0, 1.0, 1.0);
  REQUIRE(universe.getVolume() == 1.0);

  SECTION("resizing smaller changes volume")
  {
    universe.setBox(pe::Box(0.0, 1.0, 2.0));
    REQUIRE(universe.getVolume() == 0.0);
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

  SECTION("resizing bigger changes volume")
  {
    universe.setBoxLengths(1.0, 1.0, 1.0);
    REQUIRE(universe.getVolume() == 1.0);
    universe.setBox(pe::Box(2.0, 1.0, 2.0));
    REQUIRE(universe.getVolume() == 4.0);
    universe.setBox(pe::Box(-1.0, 1.0, 0.0, 1.0, 0.0, 2.0));
    REQUIRE(universe.getVolume() == 4.0);
  }

  SECTION("atoms can be added")
  {
    universe.addAtoms(2,
                      { { 0, 1 } },
                      { { 1, 1 } },
                      { { 0.0, 1.0 } },
                      { { 0.0, 1.0 } },
                      { { 0.0, 1.0 } },
                      { { 0, 0 } },
                      { { 0, 0 } },
                      { { 0, 0 } });
    REQUIRE(universe.getNrOfAtoms() == 2);
    universe.addAtoms(2,
                      { { 3, 4 } },
                      { { 1, 1 } },
                      { { 0.0, 1.0 } },
                      { { 0.0, 1.0 } },
                      { { 0.0, 1.0 } },
                      { { 0, 0 } },
                      { { 0, 0 } },
                      { { 0, 0 } });
    REQUIRE(universe.getNrOfAtoms() == 4);
    CHECK(universe.getAtomsOfType(1).size() == 4);
    CHECK(universe.getNrOfBonds() == 0);
    REQUIRE_THROWS(universe.getIdxByAtomId(1000));

    SECTION("Atoms can be removed")
    {
      universe.addBonds({ { 0, 1, 3, 4 } }, { { 1, 3, 4, 0 } });
      pe::Universe universeCopy = pe::Universe(universe);
      REQUIRE_NOTHROW(universeCopy.getAtom(3));
      REQUIRE(universeCopy.getNrOfAtoms() == 4);
      REQUIRE(universeCopy.getNrOfBonds() == 4);
      universeCopy.removeAtoms({ { 1, 3 } });
      REQUIRE(universeCopy.getAtom(4) == universe.getAtom(4));
      REQUIRE_THROWS(universeCopy.getAtom(3));
      REQUIRE(universeCopy.validate());
      REQUIRE(universeCopy.getNrOfAtoms() == 2);
      REQUIRE(universeCopy.getNrOfBonds() == 1);
    }

    SECTION("molecules are found")
    {
      REQUIRE(universe.getNrOfAtoms() == 4);
      REQUIRE(universe.getMolecules(2).size() == 4);
      // works twice: make sure there is no removal of anything happening
      REQUIRE(universe.getMolecules(2).size() == 4);
      universe.addBonds(2, { { 0, 3 } }, { { 3, 4 } });
      REQUIRE(universe.getNrOfAtoms() == 4);
      REQUIRE(universe.getNrOfBonds() == 2);
      REQUIRE(universe.getMolecules(1).size() == 0);
      std::vector<pe::Molecule> molecules = universe.getMolecules(2);
      REQUIRE(molecules.size() == 2);
      // works twice: make sure there is no removal of anything happening
      REQUIRE(universe.getMolecules(2).size() == 2);
    }

    SECTION("Special constructors work")
    {
      pe::Universe universe2 = universe;
      REQUIRE(universe2.getNrOfAtoms() == 4);
      REQUIRE(universe2.getMolecules(2).size() == 4);
    }
  }

  SECTION("Disallowed mutations are detected")
  {
    std::vector<int> threeZeros = { { 0, 0, 0 } };
    std::vector<long int> oneTwoThree = { { 1, 2, 3 } };
    std::vector<long int> threeLongZeros = { { 0, 0, 0 } };
    std::vector<double> threeDoubleZeros = { { 0, 0, 0 } };
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
    std::vector<long int> fourLongZeros = { { 0, 0, 0, 0 } };
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
    REQUIRE_THROWS(universe.addAngles(oneTwoThree, fourLongZeros, oneTwoThree));
    // different lengths
    REQUIRE_THROWS(
      universe.addBonds(3, oneTwoThree, fourLongZeros, threeZeros));
    REQUIRE_THROWS(universe.computeDxs(threeLongZeros, fourLongZeros));
  }

  SECTION("Molecules with crosslinkers are found")
  {
    /**
    # The system looks like this (in terms of bonds, not 3D placement):
    # 1-2-3-*6
    # |      |
    # *7-5---|
    # 8
    #
    # *4
    */
    universe.addAtoms(8,
                      { { 1, 2, 3, 4, 5, 6, 7, 8 } },       // id
                      { { 1, 1, 1, 2, 1, 2, 2, 1 } },       // type
                      { { 1.25, 2, 3, 1.01, 2, 4, 1, 1 } }, // x
                      { { 1, 1, 1, 4.01, 2, 1, 2, 3 } },    // y
                      { { 1, 1, 1, 1.01, 1, 1, 1, 1 } },    // z
                      { { 1, 1, 1, 1, 1, 1, 1, 1 } },       // nx
                      { { 1, 1, 1, 1, 1, 1, 1, 1 } },       // ny
                      { { 1, 1, 1, 1, 1, 1, 1, 1 } }        // nz
    );
    universe.addBonds(7,
                      { { 1, 2, 3, 6, 5, 7, 7 } },
                      { { 2, 3, 6, 5, 7, 1, 8 } },
                      { { 1, 1, 1, 1, 1, 1, 11 } });

    SECTION("Atom types are counted")
    {
      std::map<int, int> atomCounts = universe.countAtomTypes();
      REQUIRE(atomCounts.size() == 2);
      REQUIRE(atomCounts[1] == 5);
      REQUIRE(atomCounts[2] == 3);
    }

    SECTION("masses are persisted in session")
    {
      std::map<int, double> masses = universe.getMasses();
      REQUIRE(masses.size() == 0);
      masses[1] = 1.0;
      masses[2] = 2.0;
      universe.setMasses(masses);
      std::map<int, double> newMasses = universe.getMasses();
      REQUIRE(newMasses[1] == 1.0);
      REQUIRE(newMasses[2] == 2.0);
      std::map<int, double> weightFraction = universe.computeWeightFractions();
      REQUIRE(weightFraction[1] == (5.0 * 1.0) / (3.0 * 2.0 + 5.0 * 1.0));
      REQUIRE(weightFraction[2] == (2.0 * 3.0) / (3.0 * 2.0 + 5.0 * 1.0));
    }

    SECTION("get bonds returns")
    {
      auto edges = universe.getEdges();
      REQUIRE(edges["edge_from"][0] == 1 - 1);
      REQUIRE(edges["edge_to"][0] == 2 - 1);
      auto bonds = universe.getBonds();
      REQUIRE(bonds.at("bond_from")[0] == 1);
      REQUIRE(bonds.at("bond_to")[0] == 2);
      REQUIRE(bonds["bond_from"].size() == bonds["bond_to"].size());
      REQUIRE(bonds["bond_from"].size() == edges["edge_to"].size());
      REQUIRE(bonds["bond_from"].size() == edges["edge_from"].size());
      REQUIRE(bonds["bond_from"].size() == 7);
      REQUIRE(bonds["bond_from"][3] ==
              universe.getAtomIdByIdx(edges["edge_from"][3]));
      REQUIRE(bonds["bond_to"][3] ==
              universe.getAtomIdByIdx(edges["edge_to"][3]));
      REQUIRE(universe.getIdxByAtomId(bonds["bond_from"][2]) ==
              edges["edge_from"][2]);
      REQUIRE(universe.getIdxByAtomId(bonds["bond_to"][2]) ==
              edges["edge_to"][2]);
      // REQUIRE(bonds["bond_type"][5] == 1);
      // REQUIRE(bonds["bond_type"][6] == 11);
      // get atoms with type returns
      REQUIRE(universe.getAtomsOfType(2).size() == 3);
      REQUIRE(universe.getAtomsOfType(1).size() == 5);
      REQUIRE(universe.getAtomsOfType(0).size() == 0);
    }

    SECTION("internal lengths are computed correctly")
    {
      universe.setBoxLengths(10.0, 10.0, 10.0);
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

    SECTION("get angles returns")
    {
      REQUIRE(universe.getAngles()["angle_from"].size() == 0);
      auto detectedAngles = universe.detectAngles();
      universe.addAngles(detectedAngles["angle_from"],
                         detectedAngles["angle_to"],
                         detectedAngles["angle_via"]);
      REQUIRE(universe.getAngles()["angle_from"].size() == 6);
    }

    SECTION("get atoms returns")
    {
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

    SECTION("get atoms with crosslinkers returns")
    {
      // get atoms with crosslinkers returns
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

    SECTION("Loops are found")
    {
      std::map<int, std::vector<std::vector<pe::Atom>>> loops =
        universe.findLoops(2, -1, true);
      REQUIRE(loops.contains(2));
      REQUIRE(loops.size() == 1);
    }

    SECTION("Clusters are found")
    {
      std::vector<pe::Molecule> clusters = universe.getClusters();
      REQUIRE(clusters.size() == 2);
      REQUIRE(clusters[1].getAtoms()[0].getId() == 4);
      REQUIRE_THROWS(clusters[0].getAtomsLinedUp());
    }

    SECTION("Loops are found in reduced universe")
    {
      universe.addAtoms(2,
                        { { 9, 10 } },
                        { { 1, 1 } },
                        { { 0.0, 0.0 } },
                        { { 1.0, 0.0 } },
                        { { 1.0, 0.0 } },
                        { { 0, 1 } },
                        { { 0, 1 } },
                        { { 0, 1 } });
      universe.addBonds(
        2, { { 4, 9 } }, { { 9, 4 } }, { { 1, 1 } }, false, false);
      pe::Universe reducedUniverse = universe.getNetworkOfCrosslinker(2);
      REQUIRE(reducedUniverse.getNrOfAtoms() == 3);
      std::map<int, std::vector<std::vector<pe::Atom>>> loops =
        reducedUniverse.findLoops(2, -1);
      REQUIRE(loops.size() == 2);
      REQUIRE(reducedUniverse.getPropertyValues<int>("id").size() == 3);
      REQUIRE(reducedUniverse.getAtomTypes().size() == 3);
    }

    SECTION("Infinite Strands are found")
    {
      universe.setBoxLengths(4.0, 4.0, 2.0);
      REQUIRE(universe.hasInfiniteStrand(2, -1) == false);
      universe.addBonds(2, { { 8, 4 } }, { { 4, 1 } });
      std::map<int, std::vector<std::vector<pe::Atom>>> loops =
        universe.findLoops(2, -1);
      REQUIRE(loops.size() == 3);
      REQUIRE(universe.hasInfiniteStrand(2, -1) == true);
      REQUIRE(universe.getMeanStrandLength(2) ==
              Catch::Approx(((double)(3 + 1 + 1)) / 3.0));
    }

    SECTION("Reduction to Cross-linker-verse works")
    {
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

      SECTION("Copy constructors/assignment work")
      {
        pe::Universe newUniverse = pe::Universe(universe);
        newUniverse = reducedUniverse;
        REQUIRE(newUniverse.getNrOfBonds() == 2);
        REQUIRE(reducedUniverse.getNrOfAtoms() == 3);
        REQUIRE(newUniverse.getNrOfAtoms() == 3);
      }
    }

    SECTION("Computations work")
    {
      universe.setBoxLengths(10.0, 10.0, 10.0);

      // To understand where the following computations are coming from
      auto chains = universe.getChainsWithCrosslinker(2);
      REQUIRE(chains.size() == 3);
      REQUIRE(chains[0].getKey() == "1-2-3-6-7");
      REQUIRE(chains[1].getKey() == "5-6-7");
      REQUIRE(chains[2].getKey() == "7-8");

      auto endToEndDistances = universe.computeEndToEndDistances(2);
      REQUIRE(endToEndDistances.size() == 3);
      REQUIRE(endToEndDistances[0] == chains[0].computeEndToEndDistance());
      REQUIRE(endToEndDistances[0] ==
              Catch::Approx(sqrt(3.0 * 3.0 + 1.0 * 1.0)));
      REQUIRE(endToEndDistances[1] == chains[1].computeEndToEndDistance());
      REQUIRE(endToEndDistances[1] ==
              Catch::Approx(sqrt(3.0 * 3.0 + 1.0 * 1.0)));
      REQUIRE(endToEndDistances[2] == chains[2].computeEndToEndDistance());
      REQUIRE(endToEndDistances[2] == 1.0);
      REQUIRE(universe.computeMeanEndToEndDistance(2) ==
              Catch::Approx((2.0 * sqrt(3.0 * 3.0 + 1.0 * 1.0) + 1.0) / 3.0));
      REQUIRE(universe.computeMeanSquareEndToEndDistance(2) ==
              Catch::Approx((2.0 * (3.0 * 3.0 + 1.0 * 1.0) + 1.0) / 3.0));
    }
  }

  SECTION("Molecule Types are determined correctly")
  {
    /**
    # The system looks like this (in terms of bonds, not 3D placement):
    # 1-2-3
    #
    # *7-*6-5
    */
    universe.addAtoms(6,
                      { { 1, 2, 3, 5, 6, 7 } }, // id
                      { { 1, 1, 1, 1, 2, 2 } }, // type
                      { { 3, 2, 2, 2, 2, 2 } }, // x
                      { { 1, 1, 1, 1, 1, 1 } }, // y
                      { { 1, 1, 1, 1, 1, 1 } }, // z
                      { { 1, 1, 1, 1, 1, 1 } }, // nx
                      { { 1, 1, 1, 1, 1, 1 } }, // ny
                      { { 1, 1, 1, 1, 1, 1 } }  // nz
    );
    universe.addBonds(4, { { 1, 2, 5, 6 } }, { { 2, 3, 6, 7 } });
    REQUIRE(universe.getNrOfBonds() == 4);
    REQUIRE(universe.getMolecules(2).size() == 2);
    REQUIRE(universe.getChainsWithCrosslinker(2).size() == 2);
    REQUIRE(
      universe.computeFunctionalityForVertex(universe.getIdxByAtomId(1)) == 1);
    REQUIRE(universe.computeFunctionalityForAtom(1) == 1);
    REQUIRE(universe.getAtomsConnectedTo(universe.getIdxByAtomId(1)).size() ==
            1);
    auto chains = universe.getChainsWithCrosslinker(2);
    REQUIRE(chains.size() == 2);
    REQUIRE(chains[0].getNrOfAtoms() == 3);
    REQUIRE(chains[0].getType() == pe::MoleculeType::FREE_CHAIN);
    REQUIRE(chains[1].getNrOfAtoms() == 2);
    REQUIRE(chains[1].getType() == pe::MoleculeType::DANGLING_CHAIN);
    // other checks
    REQUIRE(chains[0].getKey() == "1-2-3");
    REQUIRE(chains[1].getKey() == "5-6");
    universe.setBoxLengths(3.0, 3.0, 3.0);
    REQUIRE(universe.determineEffectiveFunctionalityPerType()[2] ==
            (1. + 2.) / 2.);
    REQUIRE(chains[0].computeEndToEndDistance() == 1.0);
    REQUIRE(chains[1].computeEndToEndDistance() == 0.0);
    REQUIRE(universe.computeMeanEndToEndDistance(2) == 1.0);
    REQUIRE(universe.computeMeanSquareEndToEndDistance(2) == 1.0);

    SECTION("Primary loops too")
    {
      // instead of rewriting above tests,
      // let's just add new atoms & bonds.
      // the resulting new structure will look like this (in terms of bonds, not
      // 3D placement):
      // 1-2-3
      //
      // 9-8
      // | |
      // ↳-*7-*6-5
      universe.addAtoms(2,
                        { { 8, 9 } }, // id
                        { { 1, 1 } }, // type
                        { { 0, 1 } }, // x
                        { { 2, 3 } }, // y
                        { { 0, 0 } }, // z
                        { { 1, 1 } }, // nx
                        { { 1, 1 } }, // ny
                        { { 1, 1 } }  // nz
      );
      universe.addBonds(3, { { 7, 8, 9 } }, { { 8, 9, 7 } });
      REQUIRE(universe.getNrOfBonds() == 7);
      auto newChains = universe.getChainsWithCrosslinker(2);
      REQUIRE(newChains.size() == 3);
      REQUIRE(newChains[2].getType() == pe::MoleculeType::PRIMARY_LOOP);
      REQUIRE(newChains[2].getKey() == "7-8-9");
    }
  }

  SECTION("PolyDispersity Calculation")
  {
    // as taken from https://www.pslc.ws/macrog/average.htm
    std::vector<int> nrOfMolecules = {
      { 1, 3, 5, 8, 10, 13, 20, 13, 10, 8, 5, 3, 1 }
    };
    std::vector<double> massOfMolecules = {
      { 8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3., 2.5, 2 }
    };

    std::map<int, double> weights;
    int currentId = 0;

    std::vector<double> oneZeroDouble;
    oneZeroDouble.push_back(0.0);
    std::vector<int> oneZeroInt;
    oneZeroInt.push_back(0);

    // first, add these molecules as single atoms to the universe
    for (int i = 0; i < nrOfMolecules.size(); ++i) {
      std::vector<int> type;
      type.push_back(i);
      // this is slow and could easily be optimized, but whatever
      for (int j = 0; j < nrOfMolecules[i]; ++j) {
        std::vector<long int> atomId;
        atomId.push_back(currentId);
        universe.addAtoms(1,
                          atomId,
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
    REQUIRE(universe.getMolecules(-1).size() == 100);
    REQUIRE(currentId == 100);
    // then, do the calculation
    REQUIRE(universe.computeTotalMass() == Catch::Approx(50000000));
    REQUIRE(universe.computeNumberAverageMolecularWeight(-1) ==
            Catch::Approx(500000));
    REQUIRE(universe.computeWeightAverageMolecularWeight(-1) ==
            Catch::Approx(531600));
    REQUIRE(universe.computePolydispersityIndex(-1) == Catch::Approx(1.0632));
  }
}
