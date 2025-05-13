#include "../../src/pylimer_tools_cpp/utils/ExtraEigenTypes.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <random>
#include <string>

TEST_CASE("Eigen behaves as required", "[analysis][MEHPForceBalance][Eigen]")
{
  std::cout << "Running test \"Eigen behaves as required\"" << std::endl;
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

  // SECTION("Vector addition/subtraction is as expected")
  // {
  //   Eigen::Vector3d coords = Eigen::Vector3d::Zero();
  //   coords << 15.0609, 1.663, 2.32802;
  //   Eigen::Vector3d other = Eigen::Vector3d::Zero();
  //   other << 12.2155, 1.36349, 6.70744;

  //   Eigen::Vector3d offset = Eigen::Vector3d::Zero();
  //   offset << -0., -0., -0.;

  //   Eigen::Vector3d diff = other - coords + offset;
  //   CHECK(diff[0] == Catch::Approx(-2.84534));
  //   CHECK(diff[1] == Catch::Approx(-0.299508));
  //   CHECK(diff[2] == Catch::Approx(4.37942));
  // }

  SECTION("Entries can be swapped")
  {
    Eigen::VectorXi testVec(11);
    testVec = Eigen::VectorXi::LinSpaced(11, 0, 10);
    CHECK(testVec[4] == 4);
    std::swap(testVec[3], testVec[4]);
    CHECK(testVec[3] == 4);
  }
}

TEST_CASE("Eigen Median computation", "[Eigen]")
{
  std::cout << "Running test \"Eigen Median computation\"" << std::endl;

  SECTION("Modifiable vector")
  {
    Eigen::VectorXd testVec(11);
    testVec = Eigen::VectorXd::LinSpaced(11, 0, 10);
    CHECK(Eigen::median(testVec) == 5.0);

    Eigen::VectorXd testVec2(10);
    testVec2 << 3, 4, 1, 1, 1, 3.5, 4, 1, 1, 1;
    CHECK(Eigen::median(testVec2) == 1.);

    Eigen::VectorXd testVec3(1);
    testVec3 << 1;
    CHECK(Eigen::median(testVec3) == 1);
  }

  SECTION("Constant vector")
  {
    const Eigen::Vector3d testVec = Eigen::Vector3d::Zero();
    CHECK(Eigen::median(testVec) == 0.);
  }
};

TEST_CASE("Self adjoint check", "[Eigen]")
{
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(2);
  Eigen::SparseMatrix<double> mat1 = Eigen::SparseMatrix<double>(5, 5);

  SECTION("Check passes for a self-adjoint matrix")
  {
    triplets.push_back(Eigen::Triplet<double>(0, 0, 1.0));
    triplets.push_back(Eigen::Triplet<double>(1, 1, 2.0));
    mat1.setFromTriplets(triplets.begin(), triplets.end());

    CHECK(Eigen::isSelfAdjoint(mat1));
  }

  SECTION("Check fails for a non-self-adjoint matrix")
  {
    triplets.push_back(Eigen::Triplet<double>(0, 0, 1.0));
    triplets.push_back(Eigen::Triplet<double>(1, 1, 2.0));
    triplets.push_back(Eigen::Triplet<double>(0, 1, 3.0));
    mat1.setFromTriplets(triplets.begin(), triplets.end());

    CHECK_FALSE(Eigen::isSelfAdjoint(mat1));
  }

  SECTION("Check fails for non-symmetric matrix")
  {
    const Eigen::SparseMatrix<double> mat2 = Eigen::SparseMatrix<double>(5, 6);
    CHECK_FALSE(Eigen::isSelfAdjoint(mat2));
  }
}
