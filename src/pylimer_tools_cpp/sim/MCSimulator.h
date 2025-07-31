
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

    size_t nBeads = coordinates.size() / 3;
    size_t nBeadsToMove = nBeads - (fixFirst ? 1 : 0) - (fixLast ? 1 : 0);

    double normalisationFactor =
      std::pow(3. / (2. * M_PI * meanSquaredBeadDistance), 3. / 2.);
    double normalisationFactorInExponential =
      -3. / (2. * meanSquaredBeadDistance);
    double stepSize = std::cbrt(meanSquaredBeadDistance) * 0.5;

    do {
      iterations += 1;
      numLastStepsAccepted = 0;

      // TODO: improve performance by doing more than one bead at a time
      for (Eigen::Index i = (fixFirst ? 1 : 0);
           i < (fixLast ? nBeads - 1 : nBeads);
           ++i) {
        double bondLen21 = i == 0 ? 0.
                                  : (coordinates.segment(3 * i, 3) -
                                     coordinates.segment(3 * (i - 1), 3))
                                      .squaredNorm();
        double bondLen22 = i == nBeads - 1
                             ? 0.
                             : (coordinates.segment(3 * i, 3) -
                                coordinates.segment(3 * (i + 1), 3))
                                 .squaredNorm();

        double currentProbability =
          std::exp(normalisationFactorInExponential * (bondLen21 + bondLen22));

        Eigen::Vector3d displacement =
          stepSize * Eigen::Vector3d(displacementSamplingDist(rng),
                                     displacementSamplingDist(rng),
                                     displacementSamplingDist(rng));

        double newBondLen21 =
          i == 0 ? 0.
                 : (coordinates.segment(3 * i, 3) + displacement -
                    coordinates.segment(3 * (i - 1), 3))
                     .squaredNorm();
        double newBondLen22 =
          i == nBeads - 1 ? 0.
                          : (coordinates.segment(3 * i, 3) + displacement -
                             coordinates.segment(3 * (i + 1), 3))
                              .squaredNorm();

        double newProbability = std::exp(normalisationFactorInExponential *
                                         (newBondLen21 + newBondLen22));

        if ((newProbability / currentProbability) >
            //((bondLen21 + bondLen22) - (newBondLen21 + newBondLen22)) >
            probabilitySamplingDist(rng)) {
          coordinates.segment(3 * i, 3) += displacement;
          numLastStepsAccepted += 1;
        }
      }

      double acceptanceRatio = static_cast<double>(numLastStepsAccepted) /
                               static_cast<double>(nBeadsToMove);
      // dynamic step to target acceptance of 50%
      stepSize *= (1. + (acceptanceRatio - 0.5) / 10.);
    } while (iterations < nSteps && stepSize > 1e-5);

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
