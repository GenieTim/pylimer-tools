#ifndef RANDOM_WALKER_H
#define RANDOM_WALKER_H

#include "../entities/Box.h"
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

    double lastX = from[0];
    double lastY = from[1];
    double lastZ = from[2];

    std::array<double, 3> dist = { to[0] - from[0],
                                   to[1] - from[1],
                                   to[2] - from[2] };

    std::array<double, 3> target = box.minImageDistances(dist);

    for (int i = 0; i < chainLen; ++i) {
      std::array<double, 3> ds_uncorrected = { target[0] - lastX,
                                               target[1] - lastY,
                                               target[2] - lastZ };
      std::array<double, 3> ds = box.minImageDistances(ds);

      // for primary loops, dx, dy & dz are zero, initially.
      // therewith, alpha will be NaN
      double remainingDistance =
        std::sqrt(ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2]);
      // alpha = theta in Wikipedia
      double idealAlpha =
        std::acos(std::clamp(ds[2] / remainingDistance, -1.0, 1.0));
      // beta = phi in Wikipedia
      double idealBeta =
        (ds[0] == 0.0) ? (M_PI * 0.5) : (std::atan2(ds[1], ds[0]));
      double bondLenToUse = beadDistance;
      double idealWeight = 0.0;
      double bondsRemaining = ((chainLen - i) + 1);
      if (((remainingDistance) / (bondsRemaining)) > beadDistance) {
        // need to constrain, cannot use random alpha & beta
        // TODO: find some a bit more sophisticated probability adjustment (or
        // simply use constraints for probability)
        idealWeight = 1.0;
        bondLenToUse = remainingDistance / (bondsRemaining);
        if (bondLenToUse > 2 * beadDistance) {
          std::cout << "Using bond length: " << bondLenToUse << " for "
                    << bondsRemaining << " remaining bonds between " << lastX
                    << ", " << lastY << ", " << lastZ
                    << " to target: " << target[0] << ", " << target[1] << ", "
                    << target[2] << " with length " << remainingDistance
                    << " at i = " << i << " of " << chainLen << std::endl;
        }
        // std::min(((double)i) / ((double)chainLen) +
        //              (remainingDistance / (chainLen - i + 1)),
        //          1.0);
      }

      double alpha =
        (1.0 - idealWeight) * angleDistribution(rng) + idealWeight * idealAlpha;
      // happens e.g. for primary loops
      if (isnan(alpha)) {
        // std::cout << "Got nan for alpha with idealAlpha = " << idealAlpha
        //           << ", weight = " << idealWeight << ", dx = " << dx
        //           << ", dy = " << dy << ", dz = " << dz << " at i = " << i
        //           << std::endl;
        alpha = isnan(idealAlpha) ? angleDistribution(rng) : idealAlpha;
      };
      double beta =
        (1.0 - idealWeight) * angleDistribution(rng) + idealWeight * idealBeta;
      if (isnan(beta)) {
        // std::cout << "Got nan for beta with idealBeta = " << idealBeta
        //           << ", weight = " << idealWeight << ", dx = " << dx
        //           << ", dy = " << dy << ", dz = " << dz << " at i = " << i
        //           << std::endl;

        beta = isnan(idealBeta) ? angleDistribution(rng) : idealBeta;
      }
      // std::cout << "Using ideal weight " << idealWeight << " at " << i
      //           << ", remaining d: " << remainingDistance << " with alpha "
      //           << alpha << ", ideal " << idealAlpha << ", beta " << beta
      //           <<
      //           ", ideal " << idealBeta << std::endl;
      // coordinate system conversion: confirmation e.g. in
      // https://math.stackexchange.com/a/1385150/738831 or
      // https://en.wikipedia.org/wiki/Spherical_coordinate_system
      xs.push_back(lastX + bondLenToUse * std::cos(beta) * std::sin(alpha));
      lastX = xs[i];
      ys.push_back(lastY + bondLenToUse * std::sin(beta) * std::sin(alpha));
      lastY = ys[i];
      zs.push_back(lastZ + bondLenToUse * std::cos(alpha));
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
