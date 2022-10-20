#include "../../src/pylimer_tools_cpp/calc/TopologyCalc.h"
#include <catch2/catch_test_macros.hpp>

TEST_CASE("Segment Intersection is found", "[TopologyCalc]")
{
  REQUIRE(1 == 2 - 1);
  Eigen::Vector3d vertex1;
  vertex1 << -1.0, -1.0, 0.0;
  Eigen::Vector3d vertex2;
  vertex2 << 1.0, -1.0, 0.0;
  Eigen::Vector3d vertex3;
  vertex3 << 1.0, 1.0, 0.0;

  SECTION("Intersection is found")
  {
    Eigen::Vector3d rayOrigin;
    rayOrigin << 0.5, -0.5, -1.0;
    Eigen::Vector3d rayTarget;
    rayTarget << 0.5, -0.5, 1.0;

    Eigen::Vector3d intersectionPoint;

    REQUIRE(pylimer_tools::calc::segmentIntersectsTriangle(
      rayOrigin, rayTarget, vertex1, vertex2, vertex3, intersectionPoint));
  }

  SECTION("Intersection should not be found")
  {
    Eigen::Vector3d rayOrigin;
    rayOrigin << -0.5, 0.5, -1.0;
    Eigen::Vector3d rayTarget;
    rayTarget << -0.5, 0.5, 1.0;

    Eigen::Vector3d intersectionPoint;

    REQUIRE_FALSE(pylimer_tools::calc::segmentIntersectsTriangle(
      rayOrigin, rayTarget, vertex1, vertex2, vertex3, intersectionPoint));
  }

  SECTION("Intersection should not be found in too short segment")
  {
    Eigen::Vector3d rayOrigin;
    rayOrigin << 0.5, -0.5, -1.0;
    Eigen::Vector3d rayTarget;
    rayTarget << 0.5, -0.5, -2.0;

    Eigen::Vector3d intersectionPoint;

    REQUIRE_FALSE(pylimer_tools::calc::segmentIntersectsTriangle(
      rayOrigin, rayTarget, vertex1, vertex2, vertex3, intersectionPoint));
  }
}
