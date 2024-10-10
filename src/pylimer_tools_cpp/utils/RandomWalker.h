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

  /**
   * @brief Do a random walk of certain length from a certain starting point
   *
   * @param from the atom to start the random walk from
   * @param chainLen the number of atoms to add in between from and to
   */
  std::unordered_map<std::string, std::vector<double>> doRandomWalkChain(
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
      lastX = xs[i];
      ys.push_back(lastY + beadDistance * std::sin(beta) * std::sin(alpha));
      lastY = ys[i];
      zs.push_back(lastZ + beadDistance * std::cos(alpha));
      lastZ = zs[i];
      assert(!isnan(lastX) && !isnan(lastY) && !isnan(lastZ));
    }

    std::unordered_map<std::string, std::vector<double>> results;
    results.reserve(3);
    results.emplace("x", xs);
    results.emplace("y", ys);
    results.emplace("z", zs);
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
  std::unordered_map<std::string, std::vector<double>> doRandomWalkChainFromTo(
    const pylimer_tools::entities::Box& box,
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
    for (size_t i = 0; i < chainLen; ++i) {
      double alpha = angleDistribution(rng);
      double beta = angleDistribution(rng);
      // coordinate system conversion: confirmation e.g. in
      // https://math.stackexchange.com/a/1385150/738831 or
      // https://en.wikipedia.org/wiki/Spherical_coordinate_system
      xs.push_back(lastX + beadDistance * std::cos(beta) * std::sin(alpha));
      lastX = xs[i];
      ys.push_back(lastY + beadDistance * std::sin(beta) * std::sin(alpha));
      lastY = ys[i];
      zs.push_back(lastZ + beadDistance * std::cos(alpha));
      lastZ = zs[i];
      assert(!isnan(lastX) && !isnan(lastY) && !isnan(lastZ));
    }

    assert(xs.size() == ys.size() && xs.size() == zs.size() &&
           zs.size() == chainLen);
    // at this point, it was a "normal" random walk,
    // but now, we need to close the Brownian bridge
    // note that the bond size gets a bit destroyed
    // inspiration:
    // https://medium.com/@christopher.tabori/bridging-the-gap-an-introduction-to-brownian-bridge-simulations-10655b0baf02
    for (size_t i = 0; i < chainLen; ++i) {
      double pathFraction =
        static_cast<double>(i + 1) / static_cast<double>(chainLen + 1);
      Eigen::Vector3d deterministicPosition = fromVec + dist * pathFraction;
      xs[i] = xs[i] + deterministicPosition[0] - xs[chainLen] * pathFraction;
      ys[i] = ys[i] + deterministicPosition[1] - ys[chainLen] * pathFraction;
      zs[i] = zs[i] + deterministicPosition[2] - zs[chainLen] * pathFraction;
    }

    std::unordered_map<std::string, std::vector<double>> results;
    results.reserve(3);
    results.emplace("x", xs);
    results.emplace("y", ys);
    results.emplace("z", zs);
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
  std::unordered_map<std::string, std::vector<double>> doLinearWalkChainFromTo(
    const pylimer_tools::entities::Box& box,
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

    std::unordered_map<std::string, std::vector<double>> results;
    results.reserve(3);
    results.emplace("x", xs);
    results.emplace("y", ys);
    results.emplace("z", zs);
    return results;
  }
}
}

#endif
