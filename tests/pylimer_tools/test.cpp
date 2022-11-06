

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <unordered_map>
#include <vector>

TEST_CASE("TESTS ARE RUN", "[general]")
{
  REQUIRE(1 == 2 - 1);
}

TEST_CASE("std::map behaves as expected", "[general]") {
  std::unordered_map<size_t, size_t> testMap;
  std::vector<size_t> zeros = {0,0,0,0,0};
  CHECK(testMap.emplace(0, 100).second == true);
  for (size_t i = 0; i < zeros.size(); ++i) {
    CHECK(testMap.emplace(zeros[i], i).second == false);
  }
}
