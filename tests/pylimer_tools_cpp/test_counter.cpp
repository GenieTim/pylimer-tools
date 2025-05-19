#include "../../src/pylimer_tools_cpp/utils/Counter.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <iostream>

namespace pu = pylimer_tools::utils;

TEST_CASE("Counter counts", "[Counter][IndexCounter][header_tests][utils]")
{
  SECTION("Strings are counted correctly")
  {
    pu::Counter counter = pu::Counter<std::string>(2);
    counter.increment("A");
    counter.increment("B");
    counter.increment("A");

    std::unordered_map<std::string, int> results = counter.asMap();
    REQUIRE(results["A"] == 2);
    REQUIRE(results["B"] == 1);
    CHECK(results["C"] == 0);
  }

  SECTION("Integers are counted correctly")
  {
    pu::IndexCounter counter = pu::IndexCounter(5);
    counter.increment(2);
    counter.increment(3);
    counter.increment(2);
    counter.increment(9);
    std::unordered_map<int, int> results = counter.asMap();
    REQUIRE(results[2] == 2);
    REQUIRE(results[3] == 1);
    REQUIRE(results[4] == 0);
    CHECK(results[9] == 1);
  }
}