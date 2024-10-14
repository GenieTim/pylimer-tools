
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <numeric>
#include <random>

namespace pylimer_tools {
namespace sim {

    void equilibrateChainWithMC(Eigen::VectorXd& coordinates,
                                double beadDistance,
                                std::mt19937 rng,
                                bool fixFirst = false,
                                bool fixLast = false,
                                size_t nSteps = 1000)
    {
      INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                         "Coordinates must have a multiple of 3 elements");
      int numLastStepsAccepted = 0;
      int iterations = 0;
      double stepSize = beadDistance;

      std::uniform_real_distribution<double> probabilitySamplingDist =
        std::uniform_real_distribution<double>(0., 1.);

      size_t nBeads = coordinates.size() / 3;
      size_t nBeadsToMove = nBeads - (fixFirst ? 1 : 0) - (fixLast ? 1 : 0);

      std::vector<size_t> steppingPosIndices(nBeadsToMove);
      std::iota(
        steppingPosIndices.begin(), steppingPosIndices.end(), fixFirst ? 1 : 0);

      double normalisationFactor =
        std::pow(3. / (2. * M_PI * SQUARE(beadDistance)), 3. / 2.);
      double normalisationFactorInExponential =
        -3. / (2. * SQUARE(beadDistance));

      do {
        iterations += 1;
        numLastStepsAccepted = 0;

        if (!fixFirst) {
          double bondLen2 =
            (coordinates.segment(0, 3) - coordinates.segment(3, 3))
              .squaredNorm();

          double currentProbability =
            std::exp(normalisationFactorInExponential * bondLen2);

          Eigen::Vector3d displacement = stepSize * Eigen::Vector3d::Random();

          double newBondLen2 = (coordinates.segment(0, 3) + displacement -
                                coordinates.segment(3, 3))
                                 .squaredNorm();

          double newProbability =
            std::exp(normalisationFactorInExponential * newBondLen2);

          if ((newProbability / currentProbability) >
              probabilitySamplingDist(rng)) {
            coordinates.segment(0, 3) += displacement;
            numLastStepsAccepted += 1;
          }
        }

        for (size_t i = 1; i < nBeads - 1; ++i) {
          double bondLen21 = (coordinates.segment(3 * i, 3) -
                              coordinates.segment(3 * (i - 1), 3))
                               .squaredNorm();
          double bondLen22 = (coordinates.segment(3 * i, 3) -
                              coordinates.segment(3 * (i + 1), 3))
                               .squaredNorm();

          double currentProbability = std::exp(
            normalisationFactorInExponential * (bondLen21 + bondLen22));

          Eigen::Vector3d displacement = stepSize * Eigen::Vector3d::Random();

          double newBondLen21 = (coordinates.segment(3 * i, 3) + displacement -
                                 coordinates.segment(3 * (i - 1), 3))
                                  .squaredNorm();
          double newBondLen22 = (coordinates.segment(3 * i, 3) + displacement -
                                 coordinates.segment(3 * (i + 1), 3))
                                  .squaredNorm();

          double newProbability = std::exp(normalisationFactorInExponential *
                                           (newBondLen21 + newBondLen22));

          if ((newProbability / currentProbability) >
              probabilitySamplingDist(rng)) {
            coordinates.segment(3 * i, 3) += displacement;
            numLastStepsAccepted += 1;
          }
        }

        if (!fixLast) {
          double bondLen2 = (coordinates.segment(3 * (nBeads - 1), 3) -
                             coordinates.segment(3 * (nBeads - 2), 3))
                              .squaredNorm();

          double currentProbability =
            std::exp(normalisationFactorInExponential * bondLen2);

          Eigen::Vector3d displacement = stepSize * Eigen::Vector3d::Random();

          double newBondLen2 =
            (coordinates.segment(3 * (nBeads - 1), 3) + displacement -
             coordinates.segment(3 * (nBeads - 2), 3))
              .squaredNorm();

          double newProbability =
            std::exp(normalisationFactorInExponential * newBondLen2);

          if ((newProbability / currentProbability) >
              probabilitySamplingDist(rng)) {
            coordinates.segment(3 * (nBeads - 1), 3) += displacement;
            numLastStepsAccepted += 1;
          }
        }

        double acceptanceRatio = static_cast<double>(numLastStepsAccepted) /
                                 static_cast<double>(nBeadsToMove);
        // target acceptance of 50%
        stepSize *= (1. + (acceptanceRatio - 0.5) / 10.);
      } while (iterations < nSteps && stepSize > 1e-5);
    }

}
}
