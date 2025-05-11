#include "../../src/pylimer_tools_cpp/entities/Atom.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <iostream>
#include <string>

namespace pe = pylimer_tools::entities;

TEST_CASE("Atoms persist state", "[entity][Atom]")
{
  std::cout << "Running test \"Atoms persist state\"" << std::endl;
  pe::Atom atom1 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
  pe::Atom atom2 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
  REQUIRE(atom1 == atom2);

  pe::Atom atom3 = pe::Atom(1, 2, 1.0, 2.0, 3.0, 1, 2, 3);
  REQUIRE_FALSE(atom2 == atom3);
  REQUIRE(atom3.getId() == 1);
  REQUIRE(atom3.getType() == 2);
  REQUIRE(atom3.getX() == 1.0);
  REQUIRE(atom3.getY() == 2.0);
  REQUIRE(atom3.getZ() == 3.0);
  REQUIRE(atom3.getNX() == 1);
  REQUIRE(atom3.getNY() == 2);
  REQUIRE(atom3.getNZ() == 3);
  REQUIRE(atom3.getType() == 2);
}
