#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <string>
extern "C" {
#include <igraph/igraph.h>
}

#define CATCH_CONFIG_MAIN
namespace pe = pylimer_tools::entities;

TEST_CASE("Atoms can calculate distances", "[entity][Atom]") {
  pe::Atom atom1 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
  pe::Box unitBox = pe::Box(1.0, 1.0, 1.0);

  SECTION("Box does not change internal state") {
    REQUIRE(unitBox.getLowX() == 0.0);
    REQUIRE(unitBox.getLowY() == 0.0);
    REQUIRE(unitBox.getLowZ() == 0.0);
    REQUIRE(unitBox.getHighX() == 1.0);
    REQUIRE(unitBox.getHighY() == 1.0);
    REQUIRE(unitBox.getHighZ() == 1.0);
    REQUIRE(unitBox.getLx() == 1.0);
    REQUIRE(unitBox.getLy() == 1.0);
    REQUIRE(unitBox.getLz() == 1.0);
  }

  SECTION("Same box image distance") {
    pe::Atom atom1_2 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
    REQUIRE(atom1.distanceTo(atom1_2, &unitBox) == 0.0);
    REQUIRE(atom1_2.distanceTo(atom1, &unitBox) == 0.0);
  }

  SECTION("Different box image distance") {
    pe::Atom atom1_right = pe::Atom(0, 0, -1.0, 0.0, 0.0, 1, 0, 0);
    REQUIRE(atom1.distanceTo(atom1_right, &unitBox) == 0.0);
    pe::Atom atom1_below = pe::Atom(0, 0, 0.0, -1.0, 0.0, 0, 1, 0);
    REQUIRE(atom1.distanceTo(atom1_below, &unitBox) == 0.0);
    pe::Atom atom1_front = pe::Atom(0, 0, 0.0, 0.0, -1.0, 0, 0, 1);
    REQUIRE(atom1.distanceTo(atom1_below, &unitBox) == 0.0);
    pe::Atom atom1_lefttopback = pe::Atom(0, 0, 1.0, 1.0, 1.0, -1, -1, -1);
    REQUIRE(atom1_lefttopback.distanceTo(atom1, &unitBox) == 0.0);
  }
}
