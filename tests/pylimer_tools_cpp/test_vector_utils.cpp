#include "../../src/pylimer_tools_cpp/entities/Atom.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/utils/ExtraEigenTypes.h"
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

TEST_CASE("Vector Rows can be removed", "[Eigen]")
{
  std::cout << "Running test \"Vector Rows can be removed\"" << std::endl;
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

  SECTION("std::vector")
  {
    std::vector<int> testVec = { 1, 124, 12, 42, 41, 132, 12, 123, 5, 12, 412 };

    CHECK(testVec.size() == 11);
    CHECK(pu::contains(testVec, 1));
    CHECK(pu::contains(testVec, 12));
    CHECK_FALSE(pu::contains(testVec, -1));
    CHECK_NOTHROW(pu::removeIfContained(testVec, -1));
    CHECK(testVec.size() == 11);
    CHECK_NOTHROW(pu::removeIfContained(testVec, 12));
    CHECK_FALSE(pu::contains(testVec, 12));
    CHECK(pu::contains(testVec, 1));
    CHECK(testVec.size() == 8);
  }
}

TEST_CASE("First occurrence is found", "[VectorUtils]")
{
  std::cout << "Running test \"First occurrence is found\"" << std::endl;
  std::vector<int> testVec = { 1, 2, 3, 4, 4, 5, 5, 5, 6, 7, 8, 9, 10 };
  CHECK(pu::first_occuring_index(testVec, 5) == 5);
  CHECK(pu::first_occuring_index(testVec, 5, 6) == 5);
  CHECK(pu::first_occuring_index(testVec, 5, 7) == 5);
  CHECK(pu::first_occuring_index(testVec, 5, 2) == 5);
  CHECK(pu::first_occuring_index(testVec, 1, 0) == 0);
  CHECK(pu::first_occuring_index(testVec, 5, 8) == testVec.size());
  CHECK(pu::first_occuring_index(testVec, 11) == testVec.size());
}

TEST_CASE("Elements can be found and conditionally added", "[VectorUtils]")
{
  std::cout << "Running test \"Elements can be found and conditionally added\""
            << std::endl;
  std::vector<int> testVec;
  CHECK_THROWS(pu::last(testVec));
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

TEST_CASE("Elements are inserted to a sorted vector", "[VectorUtils]")
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

TEST_CASE("Row removal works", "[VectorUtils]")
{
  std::cout << "Running test \"Row removal works\"" << std::endl;

  SECTION("For 1D vector")
  {
    std::vector<size_t> vec = { 1, 7, 23, 4, 1, 4, 2, 5, 1,
                                2, 3, 4,  5, 6, 7, 8, 9, 10 };
    size_t sizeBefore = vec.size();

    std::vector<size_t> toRemove = std::vector<size_t>({ { 1, 2, 3, 2 } });
    pylimer_tools::utils::removeRows(vec, toRemove);
    CHECK(vec.size() == sizeBefore - 3);
    CHECK(vec[0] == 1);
    CHECK(vec[1] == 1);
    CHECK(vec[2] == 4);
  }

  SECTION("For 2D vector")
  {
    std::vector<std::vector<size_t>> vecOfVecs = {
      { 1, 2, 3, 4, 3 }, { 5, 6, 7 }, { 8, 9, 10 }, {}, { 11, 12, 13, 14, 15 }
    };
    size_t sizeBefore = vecOfVecs.size();
    std::vector<size_t> toRemove = std::vector<size_t>({ { 1, 2, 3, 2 } });
    pylimer_tools::utils::removeRows(vecOfVecs, toRemove);
    CHECK(vecOfVecs.size() == sizeBefore - 3);
    CHECK(vecOfVecs.size() == 2);
    CHECK(vecOfVecs[0][0] == 1);
    CHECK(vecOfVecs[0][1] == 2);
    CHECK(vecOfVecs[0][2] == 3);
    CHECK(vecOfVecs[0].size() == 5);
    CHECK(vecOfVecs[1][0] == 11);
    CHECK(vecOfVecs[1][1] == 12);
    CHECK(vecOfVecs[1].size() == 5);
  }

  SECTION("For Eigen vector")
  {
    Eigen::VectorXd eigenVec = Eigen::VectorXd::LinSpaced(10, 0, 90);
    CHECK(eigenVec.size() == 10);
    CHECK(eigenVec(0) == 0);
    CHECK(eigenVec(1) == 10);
    CHECK(eigenVec(2) == 20);
    std::vector<size_t> toRemove = std::vector<size_t>({ { 1, 2, 3, 2 } });
    pylimer_tools::utils::removeRows(eigenVec, toRemove);
    CHECK(eigenVec.size() == 10 - 3);
    CHECK(eigenVec(0) == 0);
    CHECK(eigenVec(1) == 40);
    CHECK(eigenVec(2) == 50);
  }
}

TEST_CASE("Index renumbering works", "[VectorUtils]")
{
  std::cout << "Running test \"Index renumbering works\"" << std::endl;
  std::vector<size_t> testVecWithIndices = {
    15, 1, 2, 3, 4, 5, 7, 9, 10, 13, 14
  };
  std::vector<size_t> removedIndices = { 0, 6, 8, 11, 12 };

  CHECK_THROWS(
    pylimer_tools::utils::getMappingForRenumbering(removedIndices, 16));

  std::ranges::sort(removedIndices, std::greater<size_t>());

  const std::vector<long int> mapping =
    pylimer_tools::utils::getMappingForRenumbering(removedIndices, 16);
  CHECK(mapping.size() == 16);

  CHECK_NOTHROW(
    pylimer_tools::utils::renumberWithMapping(testVecWithIndices, mapping));
  CHECK(testVecWithIndices[0] == 15 - removedIndices.size());
  CHECK(testVecWithIndices[1] == 0);
  CHECK(testVecWithIndices[2] == 1);
  CHECK(testVecWithIndices[3] == 2);
  CHECK(testVecWithIndices[4] == 3);
  CHECK(testVecWithIndices[5] == 4);
  CHECK(testVecWithIndices[6] == 5);
}

TEST_CASE("Append and Prepend works")
{
  std::cout << "Running test \"Append and Prepend works\"" << std::endl;
  std::vector<size_t> vec = { 1, 2, 3, 4, 5 };
  size_t sizeBefore = vec.size();

  SECTION("Basic, single value")
  {
    pylimer_tools::utils::prepend(vec, { 0 });
    CHECK(vec[0] == 0);
    CHECK(vec[1] == 1);
    pylimer_tools::utils::append(vec, { 6 });
    CHECK(vec[6] == 6);
    CHECK(vec.size() == 7);
  }

  SECTION("Prepend, multiple values")
  {
    pylimer_tools::utils::prepend(vec, { 0, 1, 2 });
    CHECK(vec[0] == 0);
    CHECK(vec[1] == 1);
    CHECK(vec[2] == 2);
    CHECK(vec[3] == 1);
    CHECK(vec.size() == 8);
  }

  SECTION("Append, multiple values")
  {
    pylimer_tools::utils::append(vec, { 3, 4, 5 });
    CHECK(vec[0] == 1);
    CHECK(vec[5] == 3);
    CHECK(vec[6] == 4);
    CHECK(vec[7] == 5);
    CHECK(vec.size() == 8);
  }

  SECTION("Inverse append")
  {
    pylimer_tools::utils::append_inverse(vec, { 10, 11, 12 });
    CHECK(vec[0] == 1);
    CHECK(vec[5] == 12);
    CHECK(vec[6] == 11);
    CHECK(vec.size() == 8);
  }

  SECTION("Inverse prepend")
  {
    pylimer_tools::utils::prepend_inverse(vec, { 7, 8, 9 });
    CHECK(vec[0] == 9);
    CHECK(vec[1] == 8);
    CHECK(vec[2] == 7);
    CHECK(vec[3] == 1);
    CHECK(vec.size() == 8);
  }

  SECTION("Append, with offset")
  {
    pylimer_tools::utils::append(vec, { 3, 4, 5 }, 1, 1);
    CHECK(vec[0] == 1);
    CHECK(vec.size() == sizeBefore + 1);
    CHECK(vec.back() == 4);
    pylimer_tools::utils::append_inverse(vec, { 6, 7, 8, 9, 10 }, 2, 3);
    CHECK(vec.size() == sizeBefore + 1);
    pylimer_tools::utils::append_inverse(vec, { 6, 7, 8, 9, 10 }, 2, 1);
    CHECK(vec.size() == sizeBefore + 3);
    CHECK(vec.back() == 7);
    CHECK(vec[vec.size() - 2] == 8);
  }

  SECTION("Prepend, with offset")
  {
    pylimer_tools::utils::prepend(vec, { 3, 4, 5 }, 1, 1);
    CHECK(vec[0] == 4);
    CHECK(vec[1] == 1);
    CHECK(vec.back() == 5);
    CHECK(vec.size() == sizeBefore + 1);
    pylimer_tools::utils::prepend_inverse(vec, { 6, 7, 8, 9, 10 }, 2, 3);
    CHECK(vec.size() == sizeBefore + 1);
    pylimer_tools::utils::prepend_inverse(vec, { 6, 7, 8, 9, 10 }, 2, 1);
    CHECK(vec.size() == sizeBefore + 3);
    CHECK(vec.back() == 5);
    CHECK(vec[0] == 8);
    CHECK(vec[1] == 7);
    CHECK(vec[2] == 4);
    CHECK(vec[3] == 1);
  }
}

TEST_CASE("Duplicates are removed", "[VectorUtils]")
{
  std::cout << "Running test \"Duplicates are removed\"" << std::endl;
  std::vector<size_t> vec = { 1, 7, 77, 7, 3, 3, 3, 4, 4, 2, 2, 4, 4, 5 };
  size_t sizeBefore = vec.size();

  pylimer_tools::utils::sort_remove_duplicates(vec);

  CHECK(vec.size() == sizeBefore - 7);
  CHECK(vec[0] == 1);
  CHECK(vec[1] == 2);
  CHECK(vec[2] == 3);
  CHECK(vec[3] == 4);
  CHECK(vec[4] == 5);
}

TEST_CASE("Maximum value is found", "[VectorUtils]")
{
  std::cout << "Running test \"Maximum value is found\"" << std::endl;
  std::vector<size_t> vec = { 1, 7, 77, 7, 3, 3, 3, 4, 4, 2, 2, 4, 4, 5 };
  size_t maxValue = pylimer_tools::utils::max_element<size_t>(vec, 0);
  CHECK(maxValue == 77);

  std::vector<double> emptyVec;
  CHECK(1. == pylimer_tools::utils::max_element<double>(emptyVec, 1.));

  std::vector<double> nearlyEmptyVec;
  nearlyEmptyVec.push_back(std::numeric_limits<double>::min());
  CHECK(std::numeric_limits<double>::min() ==
        pylimer_tools::utils::max_element<double>(nearlyEmptyVec, 1.));
}

TEST_CASE("Index of element is found", "[VectorUtils]")
{
  std::cout << "Running test \"Index of element is found\"" << std::endl;
  std::vector<size_t> vec = { 1, 7, 77, 7, 3, 3, 3, 4, 4, 2, 2, 4, 4, 5 };
  size_t index = pylimer_tools::utils::index_of<size_t>(vec, 7);
  CHECK(index == 1);
  index = pylimer_tools::utils::index_of<size_t>(vec, 3);
  CHECK(index == 4);

  CHECK_THROWS(pylimer_tools::utils::index_of<size_t>(vec, 99));

  // same for a different type of vector
  std::vector<double> vecD = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6 };
  index = pylimer_tools::utils::index_of<double>(vecD, 4.4);
  CHECK(index == 3);
  index = pylimer_tools::utils::index_of<double>(vecD, 6.6);
  CHECK(index == 5);
  CHECK_THROWS(pylimer_tools::utils::index_of<double>(vecD, 7.7));
}

TEST_CASE("Vector equality is checked", "[VectorUtils]")
{
  std::cout << "Running test \"Vector equality is checked\"" << std::endl;
  std::vector<size_t> vec1 = { 1, 2, 3, 4, 5 };
  std::vector<size_t> vec2 = { 1, 2, 3, 4, 5 };
  std::vector<size_t> vec3 = { 1, 2, 3, 4, 6 };
  CHECK(pylimer_tools::utils::equal(vec1, vec2));
  CHECK_FALSE(pylimer_tools::utils::equal(vec1, vec3));
  vec2.push_back(6);
  CHECK_FALSE(pylimer_tools::utils::equal(vec1, vec2));
  CHECK_FALSE(pylimer_tools::utils::equal(vec2, vec3));

  // similary for a different type of vector
  std::vector<double> vec4 = { 1.1, 2.2, 3.3, 4.4, 5.5 };
  std::vector<double> vec5 = { 1.1, 2.2, 3.3, 4.4, 5.5 };
  std::vector<double> vec6 = { 1.1, 2.2, 3.3, 4.4, 6.6 };
  CHECK(pylimer_tools::utils::equal(vec4, vec5));
  CHECK_FALSE(pylimer_tools::utils::equal(vec4, vec6));
}

TEST_CASE("Eigen and std::vector equality is checked", "[Eigen][VectorUtils]")
{
  std::cout << "Running test \"Eigen and std::vector equality is checked\""
            << std::endl;

  SECTION("Integer vectors")
  {
    std::vector<int> stdVec = { 1, 2, 3, 4, 5 };
    Eigen::VectorXi eigenVec(5);
    eigenVec << 1, 2, 3, 4, 5;

    CHECK(stdVec == eigenVec);
    CHECK(eigenVec == stdVec);

    stdVec[2] = 10;
    CHECK_FALSE(stdVec == eigenVec);
    CHECK_FALSE(eigenVec == stdVec);

    stdVec = { 1, 2, 3, 4, 5, 6 };
    CHECK_FALSE(stdVec == eigenVec);
    CHECK_FALSE(eigenVec == stdVec);
  }

  SECTION("Double vectors")
  {
    std::vector<double> stdVec = { 1.1, 2.2, 3.3, 4.4, 5.5 };
    Eigen::VectorXd eigenVec(5);
    eigenVec << 1.1, 2.2, 3.3, 4.4, 5.5;

    CHECK(stdVec == eigenVec);
    CHECK(eigenVec == stdVec);

    eigenVec(3) = 4.5;
    CHECK_FALSE(stdVec == eigenVec);
    CHECK_FALSE(eigenVec == stdVec);
  }

  SECTION("Array types")
  {
    std::vector<double> stdVec = { 1.1, 2.2, 3.3, 4.4, 5.5 };
    Eigen::ArrayXd eigenArray(5);
    eigenArray << 1.1, 2.2, 3.3, 4.4, 5.5;

    CHECK(stdVec == eigenArray);
    CHECK(eigenArray == stdVec);
  }

  SECTION("Empty vectors")
  {
    std::vector<int> emptyStdVec;
    Eigen::VectorXi emptyEigenVec(0);

    CHECK(emptyStdVec == emptyEigenVec);
    CHECK(emptyEigenVec == emptyStdVec);
  }

  SECTION("Approximate equality")
  {
    SECTION("Integer vectors")
    {
      std::vector<int> vec1 = { 1, 2, 3, 4, 5 };
      std::vector<int> vec2 = { 1, 2, 3, 4, 5 };
      std::vector<int> vec3 = { 1, 2, 3, 4, 6 };
      CHECK(pu::equal(vec1, vec2));
      CHECK_FALSE(pu::equal(vec1, vec3));
      CHECK(pu::vector_approx_equal(vec1, vec2));
      CHECK_FALSE(pu::vector_approx_equal(vec1, vec3));
      CHECK(pu::vector_approx_rel_equal(vec1, vec2));
      CHECK_FALSE(pu::vector_approx_rel_equal(vec1, vec3));
      vec2.push_back(6);
      CHECK_FALSE(pu::equal(vec1, vec2));
      CHECK_FALSE(pu::equal(vec2, vec3));
      CHECK_FALSE(pu::vector_approx_equal(vec1, vec2));
      CHECK_FALSE(pu::vector_approx_equal(vec2, vec3));
      CHECK_FALSE(pu::vector_approx_rel_equal(vec1, vec2));
      CHECK_FALSE(pu::vector_approx_rel_equal(vec2, vec3));
    }

    SECTION("Double Eigen vectors")
    {
      Eigen::VectorXd vec4 = Eigen::VectorXd::LinSpaced(6, 1, 6);
      Eigen::VectorXd vec5 = Eigen::VectorXd::LinSpaced(6, 1, 6);
      Eigen::VectorXd vec6 = Eigen::VectorXd::LinSpaced(7, 1, 7);
      CHECK(pu::vector_approx_equal(vec4, vec5));
      CHECK_FALSE(pu::vector_approx_equal(vec4, vec6));
      CHECK(pu::vector_approx_rel_equal(vec4, vec5));
      CHECK_FALSE(pu::vector_approx_rel_equal(vec4, vec6));
    }
  }
}

TEST_CASE("Segment-wise norm is calculated", "[VectorUtils]")
{
  std::cout << "Running test \"Segment-wise norm is calculated\"" << std::endl;
  Eigen::VectorXd vec1 = Eigen::VectorXd::LinSpaced(6, 1, 6);

  CHECK(pu::segmentwise_norm(vec1, 1) == vec1);
  const std::vector<double> res6 = pu::segmentwise_norm(vec1, 6);
  CHECK(res6.size() == 1);
  CHECK(res6[0] == std::sqrt(vec1.array().square().sum()));

  CHECK_THROWS(pu::segmentwise_norm(vec1, 0));
  CHECK_THROWS(pu::segmentwise_norm(vec1, 5));
}

TEST_CASE("Segment-wise norm maximum is calculated", "[VectorUtils]")
{
  std::cout << "Running test \"Segment-wise norm is calculated\"" << std::endl;
  Eigen::VectorXd vec1 = Eigen::VectorXd::LinSpaced(6, 1, 6);

  CHECK(pu::segmentwise_norm_max(vec1, 1) == 6);
  const double res6 = pu::segmentwise_norm_max(vec1, 6);
  CHECK(res6 == std::sqrt(vec1.array().square().sum()));

  CHECK_THROWS(pu::segmentwise_norm(vec1, 0));
  CHECK_THROWS(pu::segmentwise_norm(vec1, 5));
}

TEST_CASE("Finite component check is performed", "[VectorUtils]")
{
  std::cout << "Running test \"Finite component check is performed\""
            << std::endl;
  Eigen::VectorXd vec1 = Eigen::VectorXd::LinSpaced(6, 1, 6);
  Eigen::VectorXd vec2 = vec1.array() * std::numeric_limits<double>::infinity();

  CHECK(pu::all_components_finite(vec1));
  CHECK_FALSE(pu::all_components_finite(vec2));

  Eigen::VectorXd vec3 =
    vec1.array() * std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(pu::all_components_finite(vec3));

  Eigen::VectorXd vec4 =
    vec1.array() * std::numeric_limits<double>::signaling_NaN();
  CHECK_FALSE(pu::all_components_finite(vec4));
}
