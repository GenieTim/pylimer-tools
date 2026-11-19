#include "../../src/pylimer_tools_cpp/calc/Correlator.h"
#include "../../src/pylimer_tools_cpp/entities/Box.h"
#include "../../src/pylimer_tools_cpp/topo/TopologyCalc.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <iostream>
#include <limits>

TEST_CASE("Radius intersection is found", "[TopologyCalc]")
{
  std::cout << "Running test \"Radius intersection is found\"" << std::endl;
  REQUIRE(1 == 2 - 1);
  Eigen::Vector3d center;
  center << 0.0, 0.0, -0.5;

  Eigen::Vector3d queryPoint;
  queryPoint << 0.0, 0.0, 0.5;

  SECTION("Intersection is found")
  {
    Eigen::Vector3d intersectionPoint =
      pylimer_tools::topo::sampleIntersectionPoint(
        center, 1., queryPoint, 1.0, 0.0);
    CHECK(intersectionPoint[2] == Catch::Approx(0.0));
  }

  SECTION("Intersection should not be found")
  {
    CHECK_THROWS(pylimer_tools::topo::sampleIntersectionPoint(
      center, 0.1, queryPoint, 0.1, 0.0));
  }
}

TEST_CASE("Segment Intersection is found", "[TopologyCalc]")
{
  std::cout << "Running test \"Segment Intersection is found\"" << std::endl;
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

    REQUIRE(pylimer_tools::topo::segmentIntersectsTriangle(
      rayOrigin, rayTarget, vertex1, vertex2, vertex3, intersectionPoint));
  }

  SECTION("Intersection should not be found")
  {
    Eigen::Vector3d rayOrigin;
    rayOrigin << -0.5, 0.5, -1.0;
    Eigen::Vector3d rayTarget;
    rayTarget << -0.5, 0.5, 1.0;

    Eigen::Vector3d intersectionPoint;

    REQUIRE_FALSE(pylimer_tools::topo::segmentIntersectsTriangle(
      rayOrigin, rayTarget, vertex1, vertex2, vertex3, intersectionPoint));
  }

  SECTION("Intersection should not be found in too short segment")
  {
    Eigen::Vector3d rayOrigin;
    rayOrigin << 0.5, -0.5, -1.0;
    Eigen::Vector3d rayTarget;
    rayTarget << 0.5, -0.5, -2.0;

    Eigen::Vector3d intersectionPoint;

    REQUIRE_FALSE(pylimer_tools::topo::segmentIntersectsTriangle(
      rayOrigin, rayTarget, vertex1, vertex2, vertex3, intersectionPoint));
  }

  // SECTION("Intersection found with PBC")
  // {
  //   Eigen::Vector3d rayOrigin;
  //   rayOrigin << 4.5, 3.5, 3.0;
  //   Eigen::Vector3d rayTarget;
  //   rayTarget << 4.5, 3.5, 5.0;

  //   Eigen::Vector3d intersectionPoint;

  //   pylimer_tools::entities::Box box =
  //     pylimer_tools::entities::Box(-2., 2., -2., 2., -2., 2.);

  //   REQUIRE(pylimer_tools::calc::segmentIntersectsTriangle(
  //     rayOrigin,
  //     rayTarget,
  //     vertex1,
  //     vertex2,
  //     vertex3,
  //     intersectionPoint,
  //     [&](Eigen::Vector3d vec) {
  //       box.handlePBC<Eigen::Vector3d>(vec);
  //       return vec;
  //     }));
  // }
}

TEST_CASE("Correlator validation and edge cases", "[calc][Correlator]")
{
  SECTION("Invalid configuration parameters are rejected")
  {
    CHECK_THROWS_AS(pylimer_tools::calc::Correlator(0, 16, 2),
                    std::invalid_argument);
    CHECK_THROWS_AS(pylimer_tools::calc::Correlator(4, 0, 2),
                    std::invalid_argument);
    CHECK_THROWS_AS(pylimer_tools::calc::Correlator(4, 16, 0),
                    std::invalid_argument);
    CHECK_THROWS_AS(pylimer_tools::calc::Correlator(4, 16, 1),
                    std::invalid_argument);
    CHECK_THROWS_AS(pylimer_tools::calc::Correlator(4, 2, 4),
                    std::invalid_argument);
    CHECK_THROWS_AS(pylimer_tools::calc::Correlator(100000, 100000, 2),
                    std::invalid_argument);
  }

  SECTION("Single-level correlator evaluates finite results")
  {
    pylimer_tools::calc::Correlator corr(1, 4, 1);
    for (int i = 1; i <= 6; ++i) {
      corr.add(static_cast<double>(i));
    }
    corr.evaluate(false);
    const Eigen::VectorXd times = corr.getTimes();
    const Eigen::VectorXd values = corr.getCorrelation();
    CHECK(corr.getNumSamples() == 6);
    CHECK(corr.getNumCorrelators() == 1);
    REQUIRE(times.size() == 4);
    REQUIRE(values.size() == 4);
    for (Eigen::Index i = 0; i < times.size(); ++i) {
      CHECK(std::isfinite(times(i)));
      CHECK(std::isfinite(values(i)));
      CHECK_THAT(times(i),
                 Catch::Matchers::WithinAbs(static_cast<double>(i), 1e-12));
    }

    corr.evaluate(true);
    const Eigen::VectorXd normValues = corr.getCorrelation();
    REQUIRE(normValues.size() == 4);
    for (Eigen::Index i = 0; i < normValues.size(); ++i) {
      CHECK(std::isfinite(normValues(i)));
    }
  }

  SECTION("setsize resets an already-populated correlator cleanly")
  {
    pylimer_tools::calc::Correlator corr(4, 8, 2);
    for (int i = 0; i < 10; ++i) {
      corr.add(static_cast<double>(i) + 1.0);
    }
    corr.evaluate(false);
    REQUIRE(corr.getNumSamples() == 10);
    REQUIRE(corr.npcorr > 0);

    corr.setsize(2, 4, 2);
    CHECK(corr.getNumSamples() == 0);
    CHECK(corr.getNumCorrelators() == 2);
    CHECK(corr.p == 4);
    CHECK(corr.m == 2);
    CHECK(corr.npcorr == 0);
    CHECK_THAT(corr.accval, Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK(corr.getTimes().size() == 0);
    CHECK(corr.getCorrelation().size() == 0);

    for (int i = 0; i < 4; ++i) {
      corr.add(2.0);
    }
    corr.evaluate(false);
    CHECK(corr.getNumSamples() == 4);
    const Eigen::VectorXd values = corr.getCorrelation();
    REQUIRE(values.size() == 4);
    for (Eigen::Index i = 0; i < values.size(); ++i) {
      CHECK_THAT(values(i), Catch::Matchers::WithinAbs(4.0, 1e-12));
    }
  }

  SECTION("Empty correlator evaluation does not produce NaN")
  {
    pylimer_tools::calc::Correlator corr(4, 8, 2);
    CHECK(corr.getNumSamples() == 0);
    CHECK_THAT(corr.getMean(), Catch::Matchers::WithinAbs(0.0, 1e-12));

    corr.evaluate(false);
    CHECK(corr.npcorr == 0);
    CHECK(corr.getTimes().size() == 0);
    CHECK(corr.getCorrelation().size() == 0);

    corr.evaluate(true);
    CHECK(corr.npcorr == 0);
  }

  SECTION("Non-finite sample inputs are rejected")
  {
    pylimer_tools::calc::Correlator corr(4, 8, 2);
    CHECK_THROWS_AS(corr.add(std::numeric_limits<double>::quiet_NaN()),
                    std::invalid_argument);
    CHECK_THROWS_AS(corr.add(std::numeric_limits<double>::infinity()),
                    std::invalid_argument);

    Eigen::VectorXd batch(3);
    batch << 1.0, std::numeric_limits<double>::quiet_NaN(), 2.0;
    CHECK_THROWS_AS(corr.addBatch(batch), std::invalid_argument);
    CHECK(corr.getNumSamples() == 0);
  }

  SECTION("Large negative values are handled accurately without sentinel "
          "collisions")
  {
    pylimer_tools::calc::Correlator corr(2, 4, 2);
    const double val = -5.0e10;
    for (int i = 0; i < 6; ++i) {
      corr.add(val);
    }
    corr.evaluate(false);
    REQUIRE(corr.getNumSamples() == 6);
    REQUIRE(corr.npcorr > 0);
    for (Eigen::Index i = 0; i < corr.getCorrelation().size(); ++i) {
      CHECK_THAT(corr.getCorrelation()(i),
                 Catch::Matchers::WithinRel(val * val, 1e-12));
    }
  }
}

TEST_CASE("Correlator multi-tau evaluation and helper API",
          "[calc][Correlator]")
{
  SECTION("Highest reached correlator level k == kmax is included")
  {
    pylimer_tools::calc::Correlator corr(2, 4, 2);
    const double c = 3.5;
    for (int i = 0; i < 8; ++i) {
      corr.add(c);
    }
    corr.evaluate(false);
    const Eigen::VectorXd times = corr.getTimes();
    const Eigen::VectorXd values = corr.getCorrelation();
    // Level 0: lags 0, 1, 2, 3. Level 1 (k == kmax == 1): lags 2*2=4, 3*2=6.
    REQUIRE(times.size() == 6);
    CHECK_THAT(times(0), Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(times(3), Catch::Matchers::WithinAbs(3.0, 1e-12));
    CHECK_THAT(times(4), Catch::Matchers::WithinAbs(4.0, 1e-12));
    CHECK_THAT(times(5), Catch::Matchers::WithinAbs(6.0, 1e-12));
    for (Eigen::Index i = 0; i < values.size(); ++i) {
      CHECK_THAT(values(i), Catch::Matchers::WithinAbs(c * c, 1e-12));
    }

    corr.evaluate(true);
    for (Eigen::Index i = 0; i < corr.getCorrelation().size(); ++i) {
      CHECK_THAT(corr.getCorrelation()(i),
                 Catch::Matchers::WithinAbs(0.0, 1e-11));
    }
  }

  SECTION("Cosine signal autocorrelation matches analytical 0.5*cos(omega*tau)")
  {
    pylimer_tools::calc::Correlator corr(6, 16, 2);
    const int N = 2048;
    const double omega = 0.25;
    Eigen::VectorXd signal(N);
    for (int i = 0; i < N; ++i) {
      signal(i) = std::cos(omega * static_cast<double>(i));
    }
    corr.addBatch(signal);
    corr.evaluate(true);

    CHECK(corr.getNumSamples() == static_cast<unsigned long int>(N));
    CHECK(corr.getNumCorrelators() == 6);
    CHECK_THAT(corr.getMean(), Catch::Matchers::WithinAbs(0.0, 5e-3));

    const Eigen::VectorXd times = corr.getTimes();
    const Eigen::VectorXd values = corr.getCorrelation();
    REQUIRE(times.size() > 16);
    for (Eigen::Index i = 0; i < 16; ++i) {
      const double expected = 0.5 * std::cos(omega * times(i));
      CHECK_THAT(values(i), Catch::Matchers::WithinAbs(expected, 1e-2));
    }

    corr.reset();
    CHECK(corr.getNumSamples() == 0);
    CHECK(corr.npcorr == 0);
  }

#ifdef CEREALIZABLE
  SECTION("Cereal serialization round-trip preserves Correlator state")
  {
    pylimer_tools::calc::Correlator corr(4, 8, 2);
    for (int i = 0; i < 20; ++i) {
      corr.add(static_cast<double>(i) * 0.5);
    }
    corr.evaluate(false);

    const std::string data = pylimer_tools::utils::serializeToString(corr);
    pylimer_tools::calc::Correlator restored =
      pylimer_tools::calc::Correlator::fromString(data);
    CHECK(restored.getNumSamples() == corr.getNumSamples());
    CHECK(restored.getTimes().isApprox(corr.getTimes()));
    CHECK(restored.getCorrelation().isApprox(corr.getCorrelation()));
  }
#endif
}
