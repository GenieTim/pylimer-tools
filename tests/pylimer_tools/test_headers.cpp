#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/utils/GraphUtils.h"
#include "../../src/pylimer_tools_cpp/utils/StringUtils.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <iostream>
#include <string>
extern "C"
{
#include <igraph.h>
}

#define CATCH_CONFIG_MAIN
namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;

TEST_CASE("Atoms can calculate distances", "[entity][Atom]")
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

TEST_CASE("Box can do PBC computations", "[entity][Box]")
{
  std::cout << "Running test \"Box can do PBC computations\"" << std::endl;
  SECTION("Positive Box")
  {
    pe::Box testBox = pe::Box(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
    REQUIRE(testBox.getVolume() == Catch::Approx(10. * 10. * 10.));

    Eigen::Vector3d distances;
    distances << 10.2, 10.2, 10.2;
    REQUIRE_NOTHROW(testBox.handlePBC(distances));
    CHECK(distances[0] == Catch::Approx(10.2 - 10.));
    CHECK(distances[1] == Catch::Approx(10.2 - 10.));
    CHECK(distances[2] == Catch::Approx(10.2 - 10.));

    Eigen::VectorXd distances3(12);
    distances3 << 10.2, 10.2, 10.2, -10.2, -10.2, -10.2, 1.0, 1.0, 1.0, -4.,
      -4., -4.;
    REQUIRE_NOTHROW(testBox.handlePBC(distances3));
    for (size_t i = 0; i < 3; ++i) {
      CHECK(distances3[i] == Catch::Approx(10.2 - 10.));
      CHECK(distances3[3 + i] == Catch::Approx(-10.2 + 10.));
      CHECK(distances3[6 + i] == Catch::Approx(1.));
      CHECK(distances3[9 + i] == Catch::Approx(-4.));
    }

    std::vector<double> distances2;
    REQUIRE_NOTHROW(testBox.handlePBC(distances2));
    distances2.reserve(3);
    distances2.push_back(10.2);
    distances2.push_back(10.2);
    distances2.push_back(10.2);
    REQUIRE_NOTHROW(testBox.handlePBC(distances2));
    CHECK(distances2[0] == Catch::Approx(10.2 - 10.));
    CHECK(distances2[1] == Catch::Approx(10.2 - 10.));
    CHECK(distances2[2] == Catch::Approx(10.2 - 10.));

    // Eigen::Vector3d distances4;
    // distances4 << 1e100, 0., 0.;
    // REQUIRE_THROWS(testBox.handlePBC(distances4));

    // Eigen::Vector3d distances5;
    // distances5 << -1e100, 0., 0.;
    // REQUIRE_THROWS(testBox.handlePBC(distances5));

    Eigen::Vector3d distances6;
    distances6 << 5.0, 5.0, 5.0;
    REQUIRE_NOTHROW(testBox.handlePBC(distances6));
    CHECK(distances6[0] == Catch::Approx(5.));
    CHECK(distances6[1] == Catch::Approx(5.));
    CHECK(distances6[2] == Catch::Approx(5.));
    testBox.applySimpleShear(0.1, 2);
    REQUIRE_NOTHROW(testBox.handlePBC(distances6));
    CHECK(distances6[0] == Catch::Approx(5.));
    CHECK(distances6[1] == Catch::Approx(5.));
    CHECK(distances6[2] == Catch::Approx(5.));
    REQUIRE(testBox.getVolume() == Catch::Approx(10. * 10. * 10.));

    Eigen::Vector3d distances7;
    distances7 << 73.7435, 0.0657623, -5.26946;
    pe::Box testBox2 = pe::Box(75.8649, 75.8649, 75.8649);
    testBox2.applySimpleShear(-0.1, 2);
    REQUIRE_NOTHROW(testBox2.handlePBC(distances7));
    CHECK(distances7[0] == Catch::Approx(-2.12143).epsilon(0.001));
    CHECK_THAT(distances7[1], Catch::Matchers::WithinRel(0.0657623, 0.001));
    CHECK_THAT(distances7[2], Catch::Matchers::WithinRel(2.31703, 0.001));
  }

  SECTION("Mixed Box")
  {
    pe::Box testBox = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);
    REQUIRE(testBox.getVolume() == Catch::Approx(20. * 20. * 20.));

    Eigen::Vector3d distances;
    distances << 10.2, 10.2, 10.2;
    REQUIRE_NOTHROW(testBox.handlePBC(distances));
    CHECK(distances[0] == Catch::Approx(10.2 - 20.));
    CHECK(distances[1] == Catch::Approx(10.2 - 20.));
    CHECK(distances[2] == Catch::Approx(10.2 - 20.));
    REQUIRE_NOTHROW(testBox.minImageDistances(distances));
    CHECK(distances[0] == Catch::Approx(10.2 - 20.));
    CHECK(distances[1] == Catch::Approx(10.2 - 20.));
    CHECK(distances[2] == Catch::Approx(10.2 - 20.));

    Eigen::VectorXd distances3(12);
    distances3 << 10.2, 10.2, 10.2, -10.2, -10.2, -10.2, 1.0, 1.0, 1.0, -4.,
      -4., -4.;
    SECTION("PBC")
    {
      REQUIRE_NOTHROW(testBox.handlePBC(distances3));
      for (size_t i = 0; i < 3; ++i) {
        CHECK(distances3[i] == Catch::Approx(10.2 - 20.));
        CHECK(distances3[3 + i] == Catch::Approx(-10.2 + 20.));
        CHECK(distances3[6 + i] == Catch::Approx(1.));
        CHECK(distances3[9 + i] == Catch::Approx(-4.));
      }
    }

    SECTION("Box placement")
    {
      REQUIRE_NOTHROW(testBox.minImageDistances(distances3));
      for (size_t i = 0; i < 3; ++i) {
        CHECK(distances3[i] == Catch::Approx(10.2 - 20.));
        CHECK(distances3[3 + i] == Catch::Approx(-10.2 + 20.));
        CHECK(distances3[6 + i] == Catch::Approx(1.));
        CHECK(distances3[9 + i] == Catch::Approx(-4.));
      }
    }

    SECTION("PBC, smaller box")
    {
      pe::Box testBox2 = pe::Box(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0);
      REQUIRE_NOTHROW(testBox2.handlePBC(distances3));
      for (size_t i = 0; i < 3; ++i) {
        CHECK(distances3[i] == Catch::Approx(10.2 - 10.));
        CHECK(distances3[3 + i] == Catch::Approx(-10.2 + 10.));
        CHECK(distances3[6 + i] == Catch::Approx(1.));
        CHECK(distances3[9 + i] == Catch::Approx(-4.));
      }
    }
  }
}

TEST_CASE("Atoms compute distances PBC correctly", "[entity][Atoms][Box]")
{
  // related to test "Universe can be used" > "Local Density Computation"
  std::cout << "Running test \"Atoms compute distances PBC correctly\""
            << std::endl;
  pe::Box testBox = pe::Box(-10.0, 10.0, -10.0, 10.0, -10.0, 10.0);

  pe::Atom atom1 = pe::Atom(1, 1, 9., 9., 9., 1, 1, 1);
  pe::Atom atom2 = pe::Atom(1, 1, -10., -10., -10., 2, 2, 2);

  double distance = atom1.distanceTo(atom2, testBox);
  CHECK(distance < 2.0);
  CHECK(distance > 1.0);
}

TEST_CASE("Box can adjust coordinates", "[entity][Box]")
{
  std::cout << "Running test \"Box can adjust coordinates\"" << std::endl;
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
  std::cout << "Running test \"Box works also after simple shear\""
            << std::endl;
  pe::Box testBox = pe::Box(0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
  REQUIRE(testBox == testBox);
  REQUIRE(testBox.getVolume() == Catch::Approx(10. * 10. * 10.));
  testBox.applySimpleShear(0.1, 0);
  REQUIRE(testBox.getVolume() == Catch::Approx(10. * 10. * 10.));

  pe::Box undeformedBox = pe::Box(10.0, 10.0, 10.0);

  Eigen::VectorXd testCoords(12);
  testCoords << 10.2, 10.2, 10.2, -10.2, -10.2, -10.2, 1.0, 1.0, 1.0, -4., -4.,
    -4.;
  Eigen::VectorXd testCoordCopy = testCoords;
  Eigen::VectorXd testCoordCopy2 = testCoords;

  undeformedBox.adjustCoordinatesTo(testCoords, testBox);

  REQUIRE(testCoords[1] == Catch::Approx(10.2));
  REQUIRE(testCoords[0] == Catch::Approx(10.2 + 0.1 * 10.2));

  // deform back
  testBox.adjustCoordinatesTo(testCoords, undeformedBox);

  pe::Box deformedBox2 = pe::Box(10.0, 10.0, 10.0);
  pe::Box largeBox = pe::Box(40., 40., 40.);
  for (int dir = 0; dir < 3; ++dir) {
    deformedBox2.applySimpleShear(0.2, dir);
    // deform
    undeformedBox.adjustCoordinatesTo(testCoords, deformedBox2);
    CHECK(testCoords[(dir + 1) % 3] == Catch::Approx(10.2));
    CHECK(testCoords[dir] == Catch::Approx(10.2 + 0.2 * 10.2));
    // deform back
    deformedBox2.adjustCoordinatesTo(testCoords, undeformedBox);
    for (size_t i = 0; i < testCoords.size(); ++i) {
      CHECK(testCoords[i] == Catch::Approx(testCoordCopy[i]));
    }
    largeBox.applySimpleShear(0.0, dir);
    largeBox.handlePBC(testCoordCopy2);
    for (size_t i = 0; i < testCoordCopy.size(); ++i) {
      CHECK(testCoordCopy[i] == Catch::Approx(testCoordCopy2[i]));
    }
  }
}

TEST_CASE("Box's offset corresponds to PBC", "[entity][Box]")
{
  pe::Box box = pe::Box(10., 10., 10.);
  Eigen::Vector3d diff = Eigen::Vector3d::Constant(15.);
  CHECK(box.getOffset(diff).isApprox(Eigen::Vector3d::Constant(-20.)));
  Eigen::Vector3d diff2 = diff;
  box.handlePBC(diff2);
  CHECK(box.getOffset(diff2).isApprox(Eigen::Vector3d::Zero()));
  CHECK(diff2.isApprox(Eigen::Vector3d::Constant(-5.)));
  CHECK((diff + box.getOffset(diff)).isApprox(diff2));
  CHECK(box.isValidOffset(box.getOffset(diff)));
  CHECK_FALSE(box.isValidOffset(Eigen::Vector3d::Constant(3.2)));

  pe::Box box2 = pe::Box(97.383096, 97.383096, 97.383096);
  Eigen::Vector3d offset;
  offset << 97.3831, 2 * 97.3831, -3 * 97.3831;
  CHECK(box2.isValidOffset(offset, 1e-3));
}

TEST_CASE("Box throws", "[entity][Box]")
{
  std::cout << "Running test \"Box throws\"" << std::endl;
  REQUIRE_THROWS(pe::Box(0.0, -1.0, 0.0, -1.0, 0.0, -1.0));
}

TEST_CASE("Box can interpolate", "[entity][Box]")
{
  std::cout << "Running test \"Box can interpolate\"" << std::endl;
  pe::Box box = pe::Box(10., 10., 10.);
  pe::Box target = pe::Box(20., 5., 10.);
  CHECK(box.getVolume() == Catch::Approx(target.getVolume()));

  pe::Box result = box.interpolate(target, 0.5);
  CHECK(result.getLx() == Catch::Approx(15.));
  CHECK(result.getLy() == Catch::Approx(20. / 3.));
  CHECK(result.getLz() == Catch::Approx(10.));
  CHECK(result.getVolume() == Catch::Approx(box.getVolume()));
  CHECK(result.getShearDirection() == -1);

  result = box.interpolate(target, 1.0);
  CHECK((result.getL().isApprox(target.getL())));
  CHECK(result.getShearDirection() == -1);
  result = box.interpolate(target, 0.0);
  CHECK((result.getL().isApprox(box.getL())));
  CHECK(result.getShearDirection() == -1);

  Eigen::Array3d ls = target.getL();
  CHECK(ls[0] == Catch::Approx(20.));
  CHECK(ls[1] == Catch::Approx(5.));
  CHECK(ls[2] == Catch::Approx(10.));
  Eigen::ArrayXd lsLong = ls.replicate(3, 1);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t dir = 0; dir < 3; ++dir) {
      CHECK(lsLong[i * 3 + dir] == Catch::Approx(ls[dir]));
    }
  }
}

TEST_CASE("Box can compute bounding box", "[entity][Box]")
{
  std::cout << "Running test \"Box can compute bounding box\"" << std::endl;
  pe::Box box = pe::Box(10., 10., 10.);
  pe::Box boundingBox = box.getBoundingBox();
  CHECK(box == boundingBox);
  box.applySimpleShear(0.1, 0);
  boundingBox = box.getBoundingBox();
  CHECK(box.getLy() == boundingBox.getLy());
  CHECK(box.getLz() == boundingBox.getLz());
  CHECK(box.getLx() * 1.1 == Catch::Approx(boundingBox.getLx()));
}

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

TEST_CASE("CsvTokenizer works", "[utils][StringUtil]")
{
  std::cout << "Running test \"CsvTokenizer works\"" << std::endl;
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

  std::vector<std::string> split;
  split.reserve(2);
  int nSplit = pu::split(split, "test, test2", ", ");
  REQUIRE(nSplit == 2);
  CHECK(split[1] == "test2");
  nSplit = pu::split(split, "test\ttest2", "\t");
  REQUIRE(nSplit == 2);
}

TEST_CASE("String utility functions work", "[utils][StringUtil]")
{
  std::cout << "Running test \"String utility functions work\"" << std::endl;
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
