#include <catch2/catch_test_macros.hpp>
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/Molecule.h"
#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include <iostream>
#include <string>
#include <filesystem>
#include <igraph/igraph.h>

namespace pe = pylimer_tools::entities;

TEST_CASE("Atoms can calculate distances", "[entity][Atom]")
{
  pe::Atom atom1 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
  pe::Box unitBox = pe::Box(1.0, 1.0, 1.0);

  SECTION("Same box image distance")
  {
    pe::Atom atom1_2 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
    REQUIRE(atom1.distanceTo(atom1_2, &unitBox) == 0.0);
    REQUIRE(atom1_2.distanceTo(atom1, &unitBox) == 0.0);
  }

  SECTION("Different box image distance")
  {
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
