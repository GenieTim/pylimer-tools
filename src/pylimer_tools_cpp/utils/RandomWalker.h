#pragma once

#include "../entities/Box.h"
#include "../sim/MCSimulator.h"
#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <random>
#include <string>
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
   * @param beadDistance
   * @param meanSquaredBeadDistance
   * @param rng
   */
  static Eigen::VectorXd doRandomWalkChain(const int chainLen,
                                           const double beadDistance,
                                           const double meanSquaredBeadDistance,
                                           std::mt19937& rng)
  {
    std::uniform_real_distribution<double> angleDistribution =
      std::uniform_real_distribution<double>(0, 2 * M_PI);
    std::normal_distribution<double> stepSizeDistribution =
      std::normal_distribution<double>(
        beadDistance,
        std::sqrt(meanSquaredBeadDistance - SQUARE(beadDistance)));

    Eigen::VectorXd coordinates = Eigen::VectorXd(3 * chainLen);

    Eigen::Vector3d lastPosition = Eigen::Vector3d::Zero();

    double stepSize = beadDistance;

    for (int i = 0; i < chainLen; ++i) {
      const double alpha = angleDistribution(rng);
      const double beta = angleDistribution(rng);

      // TODO: adjust step size according to mean squared bead distance
      stepSize = stepSizeDistribution(rng);

      // coordinate system conversion: confirmation e.g. in
      // https://math.stackexchange.com/a/1385150/738831 or
      // https://en.wikipedia.org/wiki/Spherical_coordinate_system
      Eigen::Vector3d displacement(stepSize * std::cos(beta) * std::sin(alpha),
                                   stepSize * std::sin(beta) * std::sin(alpha),
                                   stepSize * std::cos(alpha));
      coordinates.segment(3 * i, 3) = lastPosition + displacement;

#ifndef NDEBUG
      assert(APPROX_EQUAL(displacement.norm(), std::abs(stepSize), 1e-10));
#endif

      lastPosition = coordinates.segment(3 * i, 3);
      assert(!lastPosition.array().isNaN().any());
    }

    return coordinates;
  }

  static Eigen::VectorXd doRandomWalkChain(
    const int chainLen,
    const double beadDistance = 1.0,
    const double meanSquaredBeadDistance = 1.0,
    std::string seed = "")
  {
    std::mt19937 rng = [&]() {
      if (seed.empty()) {
        std::random_device rd;
        return std::mt19937(rd());
      } else {
        std::seed_seq seed2(seed.begin(), seed.end());
        return std::mt19937(seed2);
      }
    }();

    return doRandomWalkChain(
      chainLen, beadDistance, meanSquaredBeadDistance, rng);
  }

  /**
   * @brief Do a random walk of certain length to add a chain from one to
   * another atom
   *
   * @param box the simulation box for PBC
   * @param from the atom to start the random walk from
   * @param to the atom to end the random walk at
   * @param chainLen the number of atoms to add in between from and to
   * @param beadDistance the distance between beads
   * @param meanSquaredBeadDistance the mean squared distance between beads
   * @param rng the random number generator
   * @param includeEnds whether to include the first and the last positions
   */
  static Eigen::VectorXd doRandomWalkChainFromTo(
    const pylimer_tools::entities::Box& box,
    const Eigen::Vector3d& from,
    const Eigen::Vector3d& to,
    const int chainLen,
    const double beadDistance,
    const double meanSquaredBeadDistance,
    std::mt19937& rng,
    const bool includeEnds = false)
  {
    std::uniform_real_distribution<double> angleDistribution =
      std::uniform_real_distribution<double>(0, 2 * M_PI);

    Eigen::Vector3d dist = to - from;
    box.handlePBC(dist);

    // include the first and the last positions
    Eigen::VectorXd coordinates = doRandomWalkChain(
      chainLen + 2, beadDistance, meanSquaredBeadDistance, rng);
    coordinates += from.replicate(chainLen + 2, 1);

    assert(coordinates.size() == 3 * (chainLen + 2));

    // at this point, it was a "normal" random walk,
    // but now, we need to close the Brownian bridge
    // note that the bond size gets a bit destroyed
    // inspiration:
    // https://medium.com/@christopher.tabori/bridging-the-gap-an-introduction-to-brownian-bridge-simulations-10655b0baf02
    double bondLen = 0.;
    for (size_t i = 1; i < chainLen + 2; ++i) {
      double pathFraction =
        static_cast<double>(i) / static_cast<double>(chainLen + 1);
      Eigen::Vector3d deterministicPosition = (from + dist) * pathFraction;
      coordinates.segment(3 * i, 3) +=
        deterministicPosition -
        pathFraction * coordinates.segment(3 * (chainLen + 1), 3);
    }

    if (includeEnds) {
      return coordinates;
    }
    // omit first and last positions
    return coordinates.segment(3, coordinates.size() - 6);
  }

  static Eigen::VectorXd doRandomWalkChainFromTo(
    const pylimer_tools::entities::Box& box,
    const Eigen::Vector3d& from,
    const Eigen::Vector3d& to,
    const int chainLen,
    const double beadDistance = 1.0,
    const double meanSquaredBeadDistance = 1.0,
    std::string seed = "")
  {
    std::mt19937 rng = [&]() {
      if (seed.empty()) {
        std::random_device rd;
        return std::mt19937(rd());
      } else {
        std::seed_seq seed2(seed.begin(), seed.end());
        return std::mt19937(seed2);
      }
    }();
    return doRandomWalkChainFromTo(
      box, from, to, chainLen, beadDistance, meanSquaredBeadDistance, rng);
  }

  /**
   * @brief Do a random walk of certain length to add a chain from one to
   * another atom, with some MC steps to equilibrate the bond lengths
   *
   * @param box the box in which to perform the random walk
   * @param from the atom to start the random walk from
   * @param to the atom to end the random walk at
   * @param chainLen the number of atoms to add in between from and to
   * @param beadDistance
   * @param meanSquaredBeadDistance
   * @param rng the random number generator
   * @param numIterations
   */
  static Eigen::VectorXd doRandomWalkChainFromToMC(
    const pylimer_tools::entities::Box& box,
    const Eigen::Vector3d& from,
    const Eigen::Vector3d& to,
    const int chainLen,
    const double beadDistance,
    const double meanSquaredBeadDistance,
    std::mt19937& rng,
    const int numIterations)
  {
    // first, do a random walk already for initial guesses
    Eigen::VectorXd coordinates =
      doRandomWalkChainFromTo(box,
                              from,
                              to,
                              chainLen,
                              beadDistance,
                              meanSquaredBeadDistance,
                              rng,
                              true);

    // then, do some MC steps to equilibrate the bond lengths
    pylimer_tools::sim::equilibrateChainWithMC(
      coordinates, meanSquaredBeadDistance, rng, true, true, numIterations);

    // omit first and last positions
    return coordinates.segment(3, coordinates.size() - 6);
  }

  static Eigen::VectorXd doRandomWalkChainFromToMC(
    const pylimer_tools::entities::Box& box,
    const Eigen::Vector3d& from,
    const Eigen::Vector3d& to,
    const int chainLen,
    const double beadDistance = 1.0,
    const double meanSquaredBeadDistance = 1.0,
    std::string seed = "",
    const int numIterations = 1000)
  {
    std::mt19937 rng = [&]() {
      if (seed.empty()) {
        std::random_device rd;
        return std::mt19937(rd());
      } else {
        std::seed_seq seed2(seed.begin(), seed.end());
        return std::mt19937(seed2);
      }
    }();
    return doRandomWalkChainFromToMC(box,
                                     from,
                                     to,
                                     chainLen,
                                     beadDistance,
                                     meanSquaredBeadDistance,
                                     rng,
                                     numIterations);
  }

  /**
   * @brief Do a random walk of certain length to add a chain from one to
   * another atom
   *
   * @param box
   * @param from the atom to start the random walk from
   * @param to the atom to end the random walk at
   * @param chainLen the number of atoms to add in between from and to
   * @param includeEnds
   */
  static Eigen::VectorXd doLinearWalkChainFromTo(
    const pylimer_tools::entities::Box& box,
    const Eigen::Vector3d& from,
    const Eigen::Vector3d& to,
    const int chainLen,
    const bool includeEnds = false)
  {
    Eigen::Vector3d dist = to - from;
    box.handlePBC(dist);

    Eigen::VectorXd results = Eigen::VectorXd::Zero(3 * (chainLen + 2));

    for (Eigen::Index i = 0; i < chainLen + 2; ++i) {
      double denominator =
        static_cast<double>(i) / static_cast<double>(chainLen + 1);
      results.segment(3 * i, 3) = from + denominator * dist;
    }

    if (includeEnds) {
      return results;
    }
    // omit first and last positions
    return results.segment(3, results.size() - 6);
  }
}
}