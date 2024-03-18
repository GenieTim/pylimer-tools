#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/utils/StringUtils.h"
#include "../../src/pylimer_tools_cpp/utils/VectorUtils.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>
extern "C"
{
#include <igraph/igraph.h>
}

#include <Eigen/Dense>

namespace pu = pylimer_tools::utils;

TEST_CASE("Eigen behaves as required", "[analysis][MEHPForceBalance][Eigen]")
{
  SECTION("Summation works with repeated indices")
  {
    Eigen::VectorXd testVec = Eigen::VectorXd::Zero(10);
    Eigen::ArrayXi testIdx = Eigen::ArrayXi::Zero(5);
    testIdx << 0, 0, 5, 5, 1;
    testVec(testIdx) += Eigen::VectorXd::Ones(5);
    REQUIRE(testVec[5] == Catch::Approx(2.));
    REQUIRE(testVec[0] == Catch::Approx(2.));
    REQUIRE(testVec[1] == Catch::Approx(1.));
    REQUIRE(testVec[2] == 0.0);
  }

  SECTION("Summation works with different repeated indices")
  {
    Eigen::VectorXd testVec = Eigen::VectorXd::Zero(10);
    Eigen::ArrayXi testIdx = Eigen::ArrayXi::Zero(5);
    testIdx << 0, 0, 5, 5, 1;
    Eigen::ArrayXi testIdx2 = Eigen::ArrayXi::Zero(5);
    testIdx2 << 5, 0, 5, 1, 5;
    Eigen::VectorXd seq = Eigen::VectorXd::Zero(7);
    seq << 0, 1, 2, 3, 4, 5, 6;
    Eigen::VectorXd resultsVec = Eigen::VectorXd::Zero(5);
    resultsVec = seq(testIdx) + seq(testIdx2);
    REQUIRE(resultsVec[0] == Catch::Approx(5.));
    REQUIRE(resultsVec[1] == Catch::Approx(0.));
    REQUIRE(resultsVec[2] == Catch::Approx(10.));
    REQUIRE(resultsVec[3] == Catch::Approx(6.));
    REQUIRE(resultsVec[4] == Catch::Approx(6.));
  }

  SECTION("Casting bool to double results in 1.0/0.0")
  {
    auto gen = std::bind(std::uniform_int_distribution<>(0, 1),
                         std::default_random_engine());
    Eigen::Array<bool, 1, 100> boolArray;
    for (int i = 0; i < 100; i++) {
      bool b = gen();
      boolArray[i] = b;
    }
    Eigen::ArrayXd castedBoolArray = boolArray.cast<double>();
    for (int i = 0; i < 100; i++) {
      if (boolArray[i]) {
        CHECK(castedBoolArray[i] == 1.0);
      } else {
        CHECK(castedBoolArray[i] + 1e-5 == 1e-5);
      }
    }
  }

  SECTION("Empty vectors are empty")
  {
    Eigen::VectorXd v = Eigen::VectorXd::Zero(0);
    CHECK(v.size() == 0);
    Eigen::ArrayXi a = Eigen::ArrayXi::Zero(0);
    CHECK(a.size() == 0);
  }
}

TEST_CASE("Vector Rows can be removed", "[Eigen]")
{
  SECTION("VectorXi")
  {
    Eigen::VectorXi testVec(11);
    testVec = Eigen::VectorXi::LinSpaced(11, 0, 10);
    CHECK(testVec[4] == 4);
    pu::removeRow(testVec, 2);
    CHECK(testVec[4] == 5);
    CHECK(testVec.size() == 11 - 1);
    pu::removeRow(testVec, 8);
    CHECK(testVec.size() == 11 - 2);
    CHECK(testVec[4] == 5);
    CHECK_THROWS(pu::removeRow(testVec, 10));

    pu::removeRows(testVec, 0, 2);
    CHECK(testVec.size() == 11 - 2 - 2);
    CHECK(testVec[4] == 7);
  }

  // same for VectorXd
  SECTION("VectorXd")
  {
    Eigen::VectorXd testVec(11);
    testVec = Eigen::VectorXd::LinSpaced(11, 0, 10);
    CHECK(testVec[4] == 4.);
    pu::removeRow(testVec, 2);
    CHECK(testVec[4] == 5.);
    CHECK(testVec.size() == 11 - 1);
    pu::removeRow(testVec, 8);
    CHECK(testVec.size() == 11 - 2);
    CHECK(testVec[4] == 5.);
    CHECK_THROWS(pu::removeRow(testVec, 10));

    pu::removeRows(testVec, 0, 2);
    CHECK(testVec.size() == 11 - 2 - 2);
    CHECK(testVec[4] == 7.);
  }

  // and for array xi
  SECTION("ArrayXi")
  {
    Eigen::ArrayXi testVec(11);
    testVec = Eigen::VectorXi::LinSpaced(11, 0, 10);
    CHECK(testVec[4] == 4);
    pu::removeRow(testVec, 2);
    CHECK(testVec[4] == 5);
    CHECK(testVec.size() == 11 - 1);
    pu::removeRow(testVec, 8);
    CHECK(testVec.size() == 11 - 2);
    CHECK(testVec[4] == 5);
    CHECK_THROWS(pu::removeRow(testVec, 10));

    pu::removeRows(testVec, 0, 2);
    CHECK(testVec.size() == 11 - 2 - 2);
    CHECK(testVec[4] == 7);
  }
}

TEST_CASE("Elements can be found and conditionally added", "[VectorUtils]")
{
  std::vector<int> testVec;
  testVec.push_back(1);
  testVec.push_back(10000);
  testVec.push_back(99);
  testVec.push_back(-6173800);

  CHECK(pu::last(testVec) == -6173800);

  for (int val : testVec) {
    CHECK(pu::contains(testVec, val));
  }
  CHECK_FALSE(pu::contains(testVec, 100));

  CHECK(testVec.size() == 4);
  pu::addIfNotContained(testVec, 1);
  CHECK(testVec.size() == 4);
  pu::addIfNotContained(testVec, 100);
  CHECK(testVec.size() == 5);
  CHECK(pu::contains(testVec, 100));
  CHECK(pu::last(testVec) == 100);
}

TEST_CASE("Elements are inserted to a sorted vector")
{
  std::vector<size_t> vec = { 1, 3, 5, 7, 9 };
  size_t size_before = vec.size();
  pylimer_tools::utils::addToSorted<size_t>(vec, 4);
  CHECK(vec.size() == size_before + 1);
  CHECK(vec[2] == 4);
  pylimer_tools::utils::addToSorted<size_t>(vec, 0);
  for (size_t i = 1; i < vec.size(); ++i) {
    CHECK(vec[i] > vec[i - 1]);
  }
}
