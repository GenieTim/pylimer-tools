#ifndef RANDOM_WALKER_H
#define RANDOM_WALKER_H

#include "../entities/Box.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

namespace pylimer_tools {
namespace utils {

  struct Positions
  {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
  };

  /**
   * @brief Do a random walk of certain length from a certain starting point
   *
   * @param from the atom to start the random walk from
   * @param chainLen the number of atoms to add in between from and to
   */
  Positions doRandomWalkChain(int chainLen,
                              double beadDistance = 1.0,
                              std::string seed = "")
  {
    std::mt19937 rng;
    if (seed == "") {
      std::random_device rd;
      rng = std::mt19937(rd());
    } else {
      std::seed_seq seed2(seed.begin(), seed.end());
      rng = std::mt19937(seed2);
    }

    std::uniform_real_distribution<double> angleDistribution =
      std::uniform_real_distribution<double>(0, 2 * M_PI);

    std::vector<double> xs;
    xs.reserve(chainLen);
    std::vector<double> ys;
    ys.reserve(chainLen);
    std::vector<double> zs;
    zs.reserve(chainLen);

    double lastX = 0.0;
    double lastY = 0.0;
    double lastZ = 0.0;

    for (int i = 0; i < chainLen; ++i) {
      double alpha = angleDistribution(rng);
      double beta = angleDistribution(rng);

      // coordinate system conversion: confirmation e.g. in
      // https://math.stackexchange.com/a/1385150/738831 or
      // https://en.wikipedia.org/wiki/Spherical_coordinate_system
      xs.push_back(lastX + beadDistance * std::cos(beta) * std::sin(alpha));
      ys.push_back(lastY + beadDistance * std::sin(beta) * std::sin(alpha));
      zs.push_back(lastZ + beadDistance * std::cos(alpha));

#ifndef NDEBUG
      assert(
        APPROX_EQUAL(std::sqrt(SQUARE(xs[i] - lastX) + SQUARE(ys[i] - lastY) +
                               SQUARE(zs[i] - lastZ)),
                     beadDistance,
                     1e-10));
#endif

      lastX = xs[i];
      lastY = ys[i];
      lastZ = zs[i];
      assert(!isnan(lastX) && !isnan(lastY) && !isnan(lastZ));
    }

    Positions results;
    results.x = xs;
    results.y = ys;
    results.z = zs;
    return results;
  }

  /**
   * @brief Do a random walk of certain length to add a chain from one to
   * another atom
   *
   * @param from the atom to start the random walk from
   * @param to the atom to end the random walk at
   * @param chainLen the number of atoms to add in between from and to
   */
  Positions doRandomWalkChainFromTo(const pylimer_tools::entities::Box& box,
                                    std::array<double, 3> from,
                                    std::array<double, 3> to,
                                    int chainLen,
                                    double beadDistance = 1.0,
                                    std::string seed = "")
  {
    std::mt19937 rng;
    if (seed == "") {
      std::random_device rd;
      rng = std::mt19937(rd());
    } else {
      std::seed_seq seed2(seed.begin(), seed.end());
      rng = std::mt19937(seed2);
    }

    std::uniform_real_distribution<double> angleDistribution =
      std::uniform_real_distribution<double>(0, 2 * M_PI);

    std::vector<double> xs;
    xs.reserve(chainLen);
    std::vector<double> ys;
    ys.reserve(chainLen);
    std::vector<double> zs;
    zs.reserve(chainLen);

    Eigen::Vector3d fromVec(from[0], from[1], from[2]);
    Eigen::Vector3d toVec(to[0], to[1], to[2]);

    Eigen::Vector3d dist = toVec - fromVec;
    box.handlePBC(dist);

    double lastX = 0.;
    double lastY = 0.;
    double lastZ = 0.;
    for (size_t j = 0; j < chainLen; ++j) {
      double alpha = angleDistribution(rng);
      double beta = angleDistribution(rng);
      // coordinate system conversion: confirmation e.g. in
      // https://math.stackexchange.com/a/1385150/738831 or
      // https://en.wikipedia.org/wiki/Spherical_coordinate_system
      xs.push_back(lastX + beadDistance * std::cos(beta) * std::sin(alpha));
      ys.push_back(lastY + beadDistance * std::sin(beta) * std::sin(alpha));
      zs.push_back(lastZ + beadDistance * std::cos(alpha));

#ifndef NDEBUG
      assert(
        APPROX_EQUAL(std::sqrt(SQUARE(xs[j] - lastX) + SQUARE(ys[j] - lastY) +
                               SQUARE(zs[j] - lastZ)),
                     beadDistance,
                     1e-8));
#endif

      lastX = xs[j];
      lastY = ys[j];
      lastZ = zs[j];
      assert(!isnan(lastX) && !isnan(lastY) && !isnan(lastZ));
    }

    assert(xs.size() == ys.size() && xs.size() == zs.size() &&
           zs.size() == chainLen);
    // at this point, it was a "normal" random walk,
    // but now, we need to close the Brownian bridge
    // note that the bond size gets a bit destroyed
    // inspiration:
    // https://medium.com/@christopher.tabori/bridging-the-gap-an-introduction-to-brownian-bridge-simulations-10655b0baf02
    double bondLen = 0.;
    for (size_t i = 0; i < chainLen; ++i) {
      double pathFraction =
        static_cast<double>(i + 1) / static_cast<double>(chainLen + 1);
      Eigen::Vector3d deterministicPosition = fromVec + dist * pathFraction;
      xs[i] = xs[i] + deterministicPosition[0] - xs[chainLen] * pathFraction;
      ys[i] = ys[i] + deterministicPosition[1] - ys[chainLen] * pathFraction;
      zs[i] = zs[i] + deterministicPosition[2] - zs[chainLen] * pathFraction;

#ifndef NDEBUG
      if (i == 0) {
        bondLen = std::sqrt(SQUARE(xs[i] - from[0]) + SQUARE(ys[i] - from[1]) +
                            SQUARE(zs[i] - from[2]));
        assert(APPROX_REL_EQUAL(bondLen, beadDistance, 0.25));
      } else {
        assert(i > 0);
        bondLen =
          std::sqrt(SQUARE(xs[i] - xs[i - 1]) + SQUARE(ys[i] - ys[i - 1]) +
                    SQUARE(zs[i] - zs[i - 1]));
        assert(APPROX_REL_EQUAL(bondLen, beadDistance, 0.25));
      }
#endif
    }

#ifndef NDEBUG
    size_t ni = chainLen - 1;
    bondLen = std::sqrt(SQUARE(xs[ni] - from[0] + dist[0]) +
                        SQUARE(ys[ni] - from[1] + dist[1]) +
                        SQUARE(zs[ni] - from[2] + dist[2]));
    assert(APPROX_REL_EQUAL(bondLen, beadDistance, 0.25));
#endif

    Positions results;
    results.x = xs;
    results.y = ys;
    results.z = zs;
    return results;
  }

  /**
   * @brief Do a random walk of certain length to add a chain from one to
   * another atom
   *
   * @param from the atom to start the random walk from
   * @param to the atom to end the random walk at
   * @param chainLen the number of atoms to add in between from and to
   */
  Positions doRandomWalkChainFromToMC(const pylimer_tools::entities::Box& box,
                                      std::array<double, 3> from,
                                      std::array<double, 3> to,
                                      int chainLen,
                                      double beadDistance = 1.0,
                                      std::string seed = "",
                                      int numIterations = 1000)
  {
    std::mt19937 rng;
    if (seed == "") {
      std::random_device rd;
      rng = std::mt19937(rd());
    } else {
      std::seed_seq seed2(seed.begin(), seed.end());
      rng = std::mt19937(seed2);
    }

    std::uniform_real_distribution<double> probabilitySamplingDist =
      std::uniform_real_distribution<double>(0., 1.);
    std::uniform_real_distribution<double> positionSamplingDist =
      std::uniform_real_distribution<double>(-0.5, 0.5);

    Positions positions;
    positions.x.reserve(chainLen + 2);
    positions.y.reserve(chainLen + 2);
    positions.z.reserve(chainLen + 2);

    Eigen::Vector3d fromVec(from[0], from[1], from[2]);
    Eigen::Vector3d toVec(to[0], to[1], to[2]);

    Eigen::Vector3d dist = toVec - fromVec;
    box.handlePBC(dist);

    for (size_t j = 0; j < chainLen + 2; ++j) {
      double pathFraction =
        static_cast<double>(j) / static_cast<double>(chainLen + 1);
      Eigen::Vector3d deterministicPosition = fromVec + dist * pathFraction;
      positions.x.push_back(deterministicPosition[0]);
      positions.y.push_back(deterministicPosition[1]);
      positions.z.push_back(deterministicPosition[2]);
    }

    // after determining these deterministic positions,
    // first prepare for probability computations
    double normalisationFactor =
      std::pow(3. / (2. * M_PI * SQUARE(beadDistance)), 3. / 2.);
    double normalisationFactorInExponential = -3. / (2. * SQUARE(beadDistance));

    // do MC steps to update the positions
    int numLastStepsAccepted = 0;
    int iterations = 0;
    double stepSize = std::cbrt(beadDistance) * 0.05;
    do {
      iterations += 1;
      numLastStepsAccepted = 0;
      for (size_t i = 1; i <= chainLen; ++i) {
        double bondLen1 = (SQUARE(positions.x[i] - positions.x[i - 1]) +
                           SQUARE(positions.y[i] - positions.y[i - 1]) +
                           SQUARE(positions.z[i] - positions.z[i - 1]));
        double bondLen2 = (SQUARE(positions.x[i] - positions.x[i + 1]) +
                           SQUARE(positions.y[i] - positions.y[i + 1]) +
                           SQUARE(positions.z[i] - positions.z[i + 1]));

        double currentProbability =
          std::exp(normalisationFactorInExponential * (bondLen1 + bondLen2));

        double disX = positionSamplingDist(rng) * stepSize;
        double disY = positionSamplingDist(rng) * stepSize;
        double disZ = positionSamplingDist(rng) * stepSize;

        double newBondLen1 =
          (SQUARE(positions.x[i] + disX - positions.x[i - 1]) +
           SQUARE(positions.y[i] + disY - positions.y[i - 1]) +
           SQUARE(positions.z[i] + disZ - positions.z[i - 1]));
        double newBondLen2 =
          (SQUARE(positions.x[i] + disX - positions.x[i + 1]) +
           SQUARE(positions.y[i] + disY - positions.y[i + 1]) +
           SQUARE(positions.z[i] + disZ - positions.z[i + 1]));

        double newProbability =
          std::exp(normalisationFactorInExponential * (newBondLen1 + bondLen2));

        double alpha = newProbability / currentProbability;

        if (alpha >= probabilitySamplingDist(rng) && alpha >= 0.5) {
          // accept
          numLastStepsAccepted += 1;
          positions.x[i] += disX;
          positions.y[i] += disY;
          positions.z[i] += disZ;
        }
      }

      double acceptanceRatio = static_cast<double>(numLastStepsAccepted) /
                               static_cast<double>(chainLen);
      // target acceptance of 50%
      stepSize *= (1. + (acceptanceRatio - 0.5) / 10.);
    } while (iterations < numIterations && stepSize > 1e-5);

    // #ifndef NDEBUG
    //     // validate bond lengths after random walk
    for (size_t i = 1; i < chainLen + 1; ++i) {
      double bondLen = std::sqrt(SQUARE(positions.x[i] - positions.x[i - 1]) +
                                 SQUARE(positions.y[i] - positions.y[i - 1]) +
                                 SQUARE(positions.z[i] - positions.z[i - 1]));
      // assert(APPROX_REL_EQUAL(bondLen, beadDistance, 0.5));
      RUNTIME_EXP_IFN(bondLen < 10. * beadDistance,
                      "Bond length after random walk is too long: got " +
                        std::to_string(bondLen) + " at index " +
                        std::to_string(i) + ".");
    }
    // #endif

    Positions results;
    results.x = positions.x;
    results.y = positions.y;
    results.z = positions.z;

    // remove first and last again
    results.x.erase(results.x.begin());
    results.y.erase(results.y.begin());
    results.z.erase(results.z.begin());
    results.x.pop_back();
    results.y.pop_back();
    results.z.pop_back();

    // finally, when we are happy, return results
    return results;
  }

  /**
   * @brief Do a random walk of certain length to add a chain from one to
   * another atom
   *
   * @param from the atom to start the random walk from
   * @param to the atom to end the random walk at
   * @param chainLen the number of atoms to add in between from and to
   */
  Positions doLinearWalkChainFromTo(const pylimer_tools::entities::Box& box,
                                    std::array<double, 3> from,
                                    std::array<double, 3> to,
                                    int chainLen)
  {
    std::array<double, 3> dist = { to[0] - from[0],
                                   to[1] - from[1],
                                   to[2] - from[2] };

    std::array<double, 3> pbc_dist = box.minImageDistances(dist);

    std::vector<double> xs;
    xs.reserve(chainLen);
    std::vector<double> ys;
    ys.reserve(chainLen);
    std::vector<double> zs;
    zs.reserve(chainLen);

    for (size_t i = 0; i < chainLen; ++i) {
      double denominator = (i + 1) / static_cast<double>(chainLen + 1);
      std::vector<double> currentCoords = {
        from[0] + pbc_dist[0] * denominator,
        from[1] + pbc_dist[1] * denominator,
        from[2] + pbc_dist[2] * denominator,
      };
      std::vector<double> currentCoordsPbc =
        box.minImageDistances(currentCoords);
      xs.push_back(currentCoordsPbc[0]);
      ys.push_back(currentCoordsPbc[1]);
      zs.push_back(currentCoordsPbc[2]);
    }

    Positions results;
    results.x = xs;
    results.y = ys;
    results.z = zs;
    return results;
  }
}
}

#endif
