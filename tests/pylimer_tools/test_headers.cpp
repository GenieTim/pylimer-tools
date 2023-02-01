#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/utils/GraphUtils.h"
#include "../../src/pylimer_tools_cpp/utils/StringUtils.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <string>
extern "C"
{
#include <igraph/igraph.h>
}

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
    REQUIRE(atom1_right.distanceTo(atom1, &unitBox) == 0.0);

    pe::Atom atom1_below = pe::Atom(0, 0, 0.0, -1.0, 0.0, 0, 1, 0);
    REQUIRE(atom1.distanceTo(atom1_below, &unitBox) == 0.0);
    REQUIRE(atom1_below.distanceTo(atom1, &unitBox) == 0.0);

    pe::Atom atom1_front = pe::Atom(0, 0, 0.0, 0.0, -1.0, 0, 0, 1);
    REQUIRE(atom1.distanceTo(atom1_front, &unitBox) == 0.0);
    REQUIRE(atom1_front.distanceTo(atom1, &unitBox) == 0.0);

    pe::Atom atom1_lefttopback = pe::Atom(0, 0, 1.0, 1.0, 1.0, -1, -1, -1);
    REQUIRE(atom1_lefttopback.distanceTo(atom1, &unitBox) == 0.0);
  }

  SECTION("Move to the mean position in box")
  {
    pe::Atom atom1 = pe::Atom(0, 0, 0.0, 0.0, 0.0, 0, 0, 0);
    pe::Atom atom2 = pe::Atom(0, 0, 1.0, 1.0, 1.0, 0, 0, 0);
    auto meanPosition_12 = atom1.meanPositionWith(atom2, &unitBox);
    REQUIRE(meanPosition_12[0] == 0.0);
    REQUIRE(meanPosition_12[1] == 0.0);
    REQUIRE(meanPosition_12[2] == 0.0);

    pe::Atom atom3 = pe::Atom(0, 0, 0.5, 0.5, 0.5, 0, 0, 0);
    auto meanPosition_13 = atom1.meanPositionWith(atom3, &unitBox);
    REQUIRE(meanPosition_13[0] == 0.25);
    REQUIRE(meanPosition_13[1] == 0.25);
    REQUIRE(meanPosition_13[2] == 0.25);
  }
}

TEST_CASE("Box can do PBC computations", "[entity][Box]")
{
  pe::Box testBox = pe::Box(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
  REQUIRE(testBox.getVolume() == Catch::Approx(10. * 10. * 10.));

  Eigen::Vector3d distances;
  distances << 10.2, 10.2, 10.2;
  REQUIRE_NOTHROW(testBox.handlePBC(distances));
  REQUIRE(distances[0] == Catch::Approx(0.2 - 5.));
  REQUIRE(distances[1] == Catch::Approx(0.2 - 5.));
  REQUIRE(distances[2] == Catch::Approx(0.2 - 5.));

  Eigen::VectorXd distances3(9);
  distances3 << 10.2, 10.2, 10.2, -0.2, -0.2, -0.2, 1.0, 1.0, 1.0;
  REQUIRE_NOTHROW(testBox.handlePBC(distances3));
  for (size_t i = 0; i < 3; ++i) {
    REQUIRE(distances3[i] == Catch::Approx(0.2 - 5.));
    REQUIRE(distances3[3 + i] == Catch::Approx(-.2 + 5.));
    REQUIRE(distances3[6 + i] == Catch::Approx(1.));
  }

  std::vector<double> distances2;
  REQUIRE_NOTHROW(testBox.handlePBC(distances2));
  distances2.reserve(3);
  distances2.push_back(10.2);
  distances2.push_back(10.2);
  distances2.push_back(10.2);
  REQUIRE_NOTHROW(testBox.handlePBC(distances2));
  REQUIRE(distances2[0] == Catch::Approx(0.2 - 5.));
  REQUIRE(distances2[1] == Catch::Approx(0.2 - 5.));
  REQUIRE(distances2[2] == Catch::Approx(0.2 - 5.));

  Eigen::Vector3d distances4;
  distances4 << 1e100, 0., 0.;
  REQUIRE_THROWS(testBox.handlePBC(distances4));

  Eigen::Vector3d distances5;
  distances5 << -1e100, 0., 0.;
  REQUIRE_THROWS(testBox.handlePBC(distances5));
}

TEST_CASE("Box can adjust coordinates", "[entity][Box]")
{
  pe::Box testBox = pe::Box(10.0, 10.0, 10.0);
  REQUIRE(testBox.getVolume() == Catch::Approx(10. * 10. * 10.));

  pe::Box testBox2 = pe::Box(20., 5., 5.);

  Eigen::Vector3d distances;
  distances << 1., 1., 1.;

  testBox.adjustCoordinatesTo(distances, testBox2);
  REQUIRE(distances[0] == Catch::Approx(2.));
  REQUIRE(distances[1] == Catch::Approx(0.5));
  REQUIRE(distances[2] == Catch::Approx(0.5));
}

TEST_CASE("Box works also after simple shear", "[entity][Box]")
{
  pe::Box testBox = pe::Box(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
  REQUIRE(testBox.getVolume() == Catch::Approx(10. * 10. * 10.));
  testBox.applySimpleShear(0.1, 0);
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
  REQUIRE(tk2.get<float>(2) == 0.001f);
  REQUIRE(tk2.get<double>(2) == static_cast<double>(0.001));
  REQUIRE(tk2.get<long double>(2) ==
          Catch::Approx(static_cast<long double>(0.001)));
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
  REQUIRE(pu::trimLineOmitComment("  test  ") == "test  ");
  REQUIRE(pu::trimLineOmitComment(testCommentString) == "");
  REQUIRE(pu::trimLineOmitComment(testCommentString.c_str()) == "");
}

TEST_CASE("Simple Cycles are found", "[utiles][GraphUtil][SimpleCycleFinder]")
{
  // turn on attribute handling
  igraph_set_attribute_table(&igraph_cattribute_table);

  SECTION("Star")
  {
    // setup simple test graph
    std::cout << "Testing star" << std::endl;
    igraph_t starGraph;
    igraph_star(&starGraph, 7, IGRAPH_STAR_UNDIRECTED, 1);
    pu::SimpleCycleFinder starCycleFinder = pu::SimpleCycleFinder(&starGraph);
    REQUIRE(starCycleFinder.findNext().empty());
    igraph_destroy(&starGraph);
  }

  // SECTION("Lattice")
  // {
  //   return;
  //   std::cout << "Testing lattice" << std::endl;
  //   igraph_t latticeGraph;
  //   std::vector<int> dimvector = { 1, 1 };
  //   igraph_vector_t dimvector_i;

  //   igraph_vector_init(&dimvector_i, 2);
  //   pu::StdVectorToIgraphVectorT(dimvector, &dimvector_i);
  //   igraph_lattice(&latticeGraph, &dimvector_i, 3, false, true, false);

  //   pu::SimpleCycleFinder latticeCycleFinder =
  //     pu::SimpleCycleFinder(&latticeGraph);
  //   CHECK(latticeCycleFinder.findAllSimpleCycles().size() == 1);
  //   igraph_vector_destroy(&dimvector_i);
  //   igraph_destroy(&latticeGraph);
  // }

  SECTION("Ring")
  {
    std::cout << "Testing ring" << std::endl;
    igraph_t ringGraph;
    igraph_ring(&ringGraph, 10, false, true, true);
    CHECK(igraph_vcount(&ringGraph) == 10);
    pu::SimpleCycleFinder ringCycleFinder = pu::SimpleCycleFinder(&ringGraph);
    REQUIRE(ringCycleFinder.findAllSimpleCycles().size() == 1);
    igraph_destroy(&ringGraph);
  }
}
