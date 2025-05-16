#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/utils/LammpsAtomStyle.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <iostream>
#include <string>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;

TEST_CASE("Atoms persist state", "[entity][Atom][header_tests]")
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


TEST_CASE("Atoms can calculate distances", "[entity][Atom][header_tests]")
{
  std::cout << "Running test \"Atoms can calculate distances\"" << std::endl;
  pe::Atom atom1 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
  pe::Box unitBox = pe::Box(1.0, 1.0, 1.0);

  SECTION("Box does not change internal state")
  {
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

  SECTION("Same box image distance")
  {
    pe::Atom atom1_2 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
    CHECK(atom1.distanceTo(atom1_2, unitBox) == 0.0);
    CHECK(atom1_2.distanceTo(atom1, unitBox) == 0.0);
  }

  SECTION("Different box image distance")
  {
    pe::Atom atom1_right = pe::Atom(0, 0, -1.0, 0.0, 0.0, 1, 0, 0);
    CHECK(atom1.distanceTo(atom1_right, unitBox) == 0.0);
    CHECK(atom1_right.distanceTo(atom1, unitBox) == 0.0);

    pe::Atom atom1_below = pe::Atom(0, 0, 0.0, -1.0, 0.0, 0, 1, 0);
    CHECK(atom1.distanceTo(atom1_below, unitBox) == 0.0);
    CHECK(atom1_below.distanceTo(atom1, unitBox) == 0.0);

    pe::Atom atom1_front = pe::Atom(0, 0, 0.0, 0.0, -1.0, 0, 0, 1);
    REQUIRE(atom1.distanceTo(atom1_front, unitBox) == 0.0);
    REQUIRE(atom1_front.distanceTo(atom1, unitBox) == 0.0);

    pe::Atom atom1_lefttopback = pe::Atom(0, 0, 1.0, 1.0, 1.0, -1, -1, -1);
    REQUIRE(atom1_lefttopback.distanceTo(atom1, unitBox) == 0.0);
  }

  SECTION("Move to the mean position in box")
  {
    pe::Atom atom1 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
    pe::Atom atom2 = pe::Atom(0, 0, 1.0, 1.0, 1.0, 0, 0, 0);
    Eigen::Vector3d meanPosition_12 = atom1.meanPositionWith(atom2, unitBox);
    CHECK(meanPosition_12[0] == 0.0);
    CHECK(meanPosition_12[1] == 0.0);
    CHECK(meanPosition_12[2] == 0.0);

    pe::Atom atom3 = pe::Atom(0, 0, 0.5, 0.5, 0.5, 0, 0, 0);
    Eigen::Vector3d meanPosition_13 = atom1.meanPositionWith(atom3, unitBox);
    CHECK(meanPosition_13[0] == 0.25);
    CHECK(meanPosition_13[1] == 0.25);
    CHECK(meanPosition_13[2] == 0.25);

    pe::Box unitNegBox = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
    Eigen::Vector3d meanPosition_13n =
      atom1.meanPositionWith(atom3, unitNegBox);
    CHECK(meanPosition_13n[0] == 0.25);
    CHECK(meanPosition_13n[1] == 0.25);
    CHECK(meanPosition_13n[2] == 0.25);

    pe::Atom negAtom4 = pe::Atom(0, 0, -2.0, 0.0, 0.0, 0, 0, 0);
    Eigen::Vector3d meanPosition_24 =
      atom2.meanPositionWith(negAtom4, unitNegBox);
    CHECK(meanPosition_24[0] == -0.5);
    CHECK(meanPosition_24[1] == 0.5);
    CHECK(meanPosition_24[2] == 0.5);

    pe::Atom negAtom5 = pe::Atom(0, 0, -2.0, -5.0, 1.0, 1, 1, 1);
    Eigen::Vector3d meanPosition_25 =
      atom2.meanPositionWithUnwrapped(negAtom5, unitNegBox);
    CHECK(meanPosition_25[0] == ((20. - 2.) + 1.) / 2.);
    CHECK(meanPosition_25[1] == ((20. - 5.) + 1.) / 2.);
    CHECK(meanPosition_25[2] == (((20. + 1.) + 1.) / 2.) - 20.);

    meanPosition_25 = atom2.meanPositionWithUnwrapped(negAtom5, unitNegBox);
    CHECK(meanPosition_25[0] == ((20. - 2.) + 1.) / 2.);
    CHECK(meanPosition_25[1] == ((20. - 5.) + 1.) / 2.);
    CHECK(meanPosition_25[2] == (((20. + 1.) + 1.) / 2.) - 20.);
  }
}

TEST_CASE("LAMMPS Atom styles resolve to correct name", "[entity][Atom][header_tests]")
{
  std::cout << "Running test \"LAMMPS Atom styles resolve to correct name\""
            << std::endl;

#define X(e, n)                                                                \
  CHECK(pu::getAtomStyleString(pu::AtomStyle::e) == std::string(n));
  LAMMPS_ATOM_STYLES
#undef X

  CHECK_THROWS(
    pu::getAtomStyleString(static_cast<pylimer_tools::utils::AtomStyle>(100)));
}

TEST_CASE("LAMMPS Atom styles resolve correctly from name", "[entity][Atom]")
{
  std::cout << "Running test \"LAMMPS Atom styles resolve correctly from name\""
            << std::endl;

#define X(e, n)                                                                \
  CHECK(pu::getAtomStyleFromString(std::string(n)) == pu::AtomStyle::e);
  LAMMPS_ATOM_STYLES
#undef X

  CHECK_THROWS(pu::getAtomStyleFromString("UnknownAtomStyle"));
}
