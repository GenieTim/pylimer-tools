#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/utils/GraphUtils.h"
#include "../../src/pylimer_tools_cpp/utils/StringUtils.h"
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
  SECTION("Summation works with same indices")
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
}

TEST_CASE("Vector Rows can be removed", "[Eigen]")
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
