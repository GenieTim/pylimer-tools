

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>

TEST_CASE("TESTS ARE RUN", "[general]")
{
  REQUIRE(1 == 2 - 1);
}
