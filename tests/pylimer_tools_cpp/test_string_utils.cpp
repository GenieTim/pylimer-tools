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
#include <igraph/igraph.h>
}

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;

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
  REQUIRE(pu::trimLineOmitComment("  ITEM: TIMESTEP # comment\r\n") ==
          "ITEM: TIMESTEP ");
  REQUIRE(pu::trimLineOmitComment("  ITEM: TIMESTEP\r\n") == "ITEM: TIMESTEP");
  REQUIRE(pu::trimLineOmitComment(testCommentString) == "");
  REQUIRE(pu::trimLineOmitComment(testCommentString.c_str()) == "");
  const std::string vecStr = std::to_string(Eigen::Vector3d(1.0, 2.0, 3.0));
  CHECK(vecStr.find("1.0") != std::string::npos);
  CHECK(vecStr.find("2.0") != std::string::npos);
  CHECK(vecStr.find("3.0") != std::string::npos);
  const auto dur = std::chrono::hours(1) + std::chrono::minutes(2) +
                   std::chrono::seconds(3) + std::chrono::milliseconds(45);
  CHECK(std::duration_to_string(dur) == "1:02:03.045");
}
