
#pragma once

#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <numeric>
#include <random>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace pylimer_tools {
namespace sim {

  static void equilibrateChainWithMC(Eigen::VectorXd& coordinates,
                                     const double meanSquaredBeadDistance,
                                     std::mt19937& rng,
                                     const bool fixFirst = false,
                                     const bool fixLast = false,
                                     const size_t nSteps = 1000)
  {
    INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                       "Coordinates must have a multiple of 3 elements");
    int numLastStepsAccepted = 0;
    int iterations = 0;

    std::uniform_real_distribution<double> probabilitySamplingDist =
      std::uniform_real_distribution<double>(0., 1.);
    std::uniform_real_distribution<double> displacementSamplingDist =
      std::uniform_real_distribution<double>(-0.5, 0.5);

    const size_t nBeads = static_cast<size_t>(coordinates.size() / 3);
    const size_t nBeadsToMove = nBeads - (fixFirst ? 1 : 0) - (fixLast ? 1 : 0);
    if (nBeads <= 1 || nBeadsToMove == 0 || nSteps == 0) {
      return;
    }

    const double normalisationFactorInExponential =
      -3. / (2. * meanSquaredBeadDistance);
    double stepSize = std::cbrt(meanSquaredBeadDistance) * 0.5;
    double* const coords = coordinates.data();

    auto tryAccept = [&](double* const bead,
                         const double dx,
                         const double dy,
                         const double dz,
                         const double deltaBondLen2) -> void {
      // Metropolis acceptance criterion:
      // P_new / P_old = exp(normalisationFactorInExponential * deltaBondLen2).
      // Since normalisationFactorInExponential < 0, if deltaBondLen2 <= 0 the
      // move is accepted unconditionally without evaluating exp() or RNG.
      if (deltaBondLen2 <= 0.0 ||
          std::exp(normalisationFactorInExponential * deltaBondLen2) >
            probabilitySamplingDist(rng)) {
        bead[0] += dx;
        bead[1] += dy;
        bead[2] += dz;
        numLastStepsAccepted += 1;
      }
    };

    do {
      iterations += 1;
      numLastStepsAccepted = 0;

      if (!fixFirst) {
        const double dx = stepSize * displacementSamplingDist(rng);
        const double dy = stepSize * displacementSamplingDist(rng);
        const double dz = stepSize * displacementSamplingDist(rng);
        const double ux = coords[0] - coords[3];
        const double uy = coords[1] - coords[4];
        const double uz = coords[2] - coords[5];
        const double d2 = dx * dx + dy * dy + dz * dz;
        const double deltaBondLen2 = 2.0 * (ux * dx + uy * dy + uz * dz) + d2;
        tryAccept(coords, dx, dy, dz, deltaBondLen2);
      }

      for (size_t i = 1; i + 1 < nBeads; ++i) {
        double* const c = coords + 3 * i;
        const double dx = stepSize * displacementSamplingDist(rng);
        const double dy = stepSize * displacementSamplingDist(rng);
        const double dz = stepSize * displacementSamplingDist(rng);
        const double ux = 2.0 * c[0] - c[-3] - c[3];
        const double uy = 2.0 * c[1] - c[-2] - c[4];
        const double uz = 2.0 * c[2] - c[-1] - c[5];
        const double d2 = dx * dx + dy * dy + dz * dz;
        const double deltaBondLen2 = 2.0 * ((ux * dx + uy * dy + uz * dz) + d2);
        tryAccept(c, dx, dy, dz, deltaBondLen2);
      }

      if (!fixLast) {
        double* const c = coords + 3 * (nBeads - 1);
        const double dx = stepSize * displacementSamplingDist(rng);
        const double dy = stepSize * displacementSamplingDist(rng);
        const double dz = stepSize * displacementSamplingDist(rng);
        const double ux = c[0] - c[-3];
        const double uy = c[1] - c[-2];
        const double uz = c[2] - c[-1];
        const double d2 = dx * dx + dy * dy + dz * dz;
        const double deltaBondLen2 = 2.0 * (ux * dx + uy * dy + uz * dz) + d2;
        tryAccept(c, dx, dy, dz, deltaBondLen2);
      }

      const double acceptanceRatio = static_cast<double>(numLastStepsAccepted) /
                                     static_cast<double>(nBeadsToMove);
      // dynamic step to target acceptance of 50%
      stepSize *= (1. + (acceptanceRatio - 0.5) / 10.);
    } while (static_cast<size_t>(iterations) < nSteps && stepSize > 1e-5);

#ifndef NDEBUG
    // validate bond lengths
    for (size_t i = 1; i < nBeads; ++i) {
      double bondLen2 =
        (coordinates.segment(3 * i, 3) - coordinates.segment(3 * (i - 1), 3))
          .squaredNorm();
      RUNTIME_EXP_IFN(bondLen2 < 8. * meanSquaredBeadDistance,
                      "Invalid bond length after equilibration got " +
                        std::to_string(bondLen2) + " at index " +
                        std::to_string(i) + ".");
    }
#endif
  }

}
}
