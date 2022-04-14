#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/utils/StringUtils.h"
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <string>

#define CATCH_CONFIG_MAIN
namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;

TEST_CASE("Atoms can calculate distances", "[entity][Atom]")
{
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

TEST_CASE("Atoms persist state", "[entity][Atom]")
{
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

TEST_CASE("CsvTokenizer works", "[utils][StringUtil]")
{
  pu::CsvTokenizer tk1 = pu::CsvTokenizer("test, test2");
  REQUIRE(tk1.getLength() == 2);
  REQUIRE(tk1.get<std::string>(1) == "test2");
  pu::CsvTokenizer tk2 = pu::CsvTokenizer("test 12  0.001", 3);
  REQUIRE(tk2.getLength() == 3);
  REQUIRE(tk2.get<std::string>(0) == "test");
  REQUIRE(tk2.get<unsigned int>(1) == 12);
  REQUIRE(tk2.get<int>(1) == 12);
  REQUIRE(tk2.get<unsigned long int>(1) == 12);
  REQUIRE(tk2.get<long int>(1) == 12);
  REQUIRE(tk2.get<float>(2) == 0.001);
  REQUIRE(tk2.get<double>(2) == 0.001);
  REQUIRE(tk2.get<long double>(2) == 0.001);
}

TEST_CASE("String utility functions work", "[utils][StringUtil]")
{
  REQUIRE(std::to_string("test") == "test");
  REQUIRE(pu::isUpper("TEST"));
  REQUIRE_FALSE(pu::isUpper("TeST"));
  REQUIRE(pu::contains("TEST", "ES"));
  REQUIRE_FALSE(pu::contains("TEST", "es"));
  REQUIRE(pu::ltrim("  test  ") == "test  ");
  REQUIRE(pu::rtrim("  test  ") == "  test");
  REQUIRE(pu::trim("  test ") == "test");
  REQUIRE(pu::rstrip("__test__", "_") == "");
  REQUIRE(pu::rstrip("test__", "_") == "test");
  // REQUIRE(pu::lstrip("__test__", "_") == "test__");
  REQUIRE(pu::startsWith("test", "te"));
  REQUIRE_FALSE(pu::startsWith("test", "es"));
  std::string testCommentString = "#  test  ";
  REQUIRE(pu::trimLineOmitComment("  test  ") == "  test");
  REQUIRE(pu::trimLineOmitComment(testCommentString) == "");
  REQUIRE(pu::trimLineOmitComment(testCommentString.c_str()) == "");
}
