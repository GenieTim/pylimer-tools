#ifndef DPD_SIMULATOR_H
#define DPD_SIMULATOR_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/EigenNeighbourList.h"
#include "../entities/Universe.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <nlopt.hpp>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {

  namespace dpd {
    typedef Eigen::Array<size_t, Eigen::Dynamic, 1> ArrayXst;
    typedef Eigen::Array<long int, Eigen::Dynamic, 1> ArrayXli;

    class DPDSimulator
    {
    public:
      DPDSimulator(const pylimer_tools::entities::Universe u,
                   const int crosslinkerType = 2,
                   const bool is2D = false,
                   const std::string seed = "");

      /**
       * @brief actually do run the simulation
       *
       */
      void runSimulation(const long int nSteps,
                         double dt = 0.06,
                         double lambda = 0.65,
                         bool withMC = false);

      /**
       * @brief Compute the force vector, and return the pressure
       *
       */
      double computeForces(
        Eigen::VectorXd& forces,
        Eigen::Matrix3d& stressTensor,
        const Eigen::VectorXd& coordinates,
        const Eigen::VectorXd& velocities,
        const double dt = 0.06,
        const double cutoff =
          1.0); // unfortunately not const because of the random nr generator

      /**
       * @brief Compute the temperature
       *
       * @param velocities
       * @return double
       */
      double computeTemperature(const Eigen::VectorXd& velocities) const;

    /**
     * @brief Get access to the current stress-tensor
     * 
     * @return Eigen::Matrix3d 
     */
      Eigen::Matrix3d getStressTensor() const;

      /**
       * @brief Randomly add new slip-springs
       */
      int createSlipSprings(const int num);

      int shiftSlipSprings(const double kbT = 1.);

      int relocateSlipSprings(const double kbT = 1.);

      void configLambda(const double l) { this->lambda = l; }

      void configSpringConstant(const double nk) { this->k = nk; }

      void configSlipspringLowCutoff(const double lowC)
      {
        INVALIDARG_EXP_IFN(
          lowC < this->highCutoff,
          "The low cutoff must be lower than the high cutoff.");
        this->lowCutoff = lowC;
      }

      void configSlipspringHighCutoff(const double highC)
      {
        INVALIDARG_EXP_IFN(
          this->lowCutoff < highC,
          "The low cutoff must be lower than the high cutoff.");
        this->highCutoff = highC;
      }

      void configSigma(const double sigma)
      {
        this->sigma = sigma;
        this->gamma = 0.5 * sigma * sigma;
      }

      void configA(const double newA) { this->A = newA; }

      void validateState();

    protected:
      void addSlipSprings(std::vector<size_t>& partnerA,
                          std::vector<size_t>& partnerB,
                          const int bondType = 9);

      bool attemptSlipSpringShift(const size_t springIdx,
                                  const size_t endToShift,
                                  const double kbT = 1.);

      void replaceSlipSpringPartner(const size_t springIdx,
                                    const size_t partnerBefore,
                                    const size_t partnerAfter);

    private:
      bool is2D;

      // configuration
      double lambda = 0.65;
      double k = 2.;
      double lowCutoff = 0.5;
      double highCutoff = 2.0;
      double A = 25.;
      double sigma = 3;
      double gamma = 0.5 * 3 * 3;

      // simulation state
      long int currentStep = 0;
      Eigen::VectorXd currentVelocitiesPlus;
      Eigen::VectorXd currentVelocities;
      Eigen::VectorXd currentForces;
      Eigen::Matrix3d currentStressTensor;

      // randomness
      std::mt19937 e2;
      std::uniform_real_distribution<double> uniform_rand_mean0std1;
      std::uniform_real_distribution<double> uniform_rand_between_0_1;

      // universe structure
      int numAtoms = 0;
      int numBonds = 0;
      int numSlipSprings = 0;
      pylimer_tools::entities::Box box;
      // atoms
      Eigen::VectorXd coordinates;
      Eigen::ArrayXi idxFunctionalities;
      std::vector<int> atomTypes;
      std::vector<long int> atomIds;
      std::vector<size_t> chainEndIndices;
      // bonds
      Eigen::ArrayXi bondPartnersA;
      Eigen::ArrayXi bondPartnersB;
      Eigen::ArrayXi bondTypes;
      // mapping between atoms and bonds
      std::vector<std::vector<size_t>> bondsOfIndex;

      pylimer_tools::entities::EigenNeighbourList neighbourlist;
    };
  }

}
}

#endif
