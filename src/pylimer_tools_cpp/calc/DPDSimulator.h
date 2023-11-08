#ifndef DPD_SIMULATOR_H
#define DPD_SIMULATOR_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/EigenNeighbourList.h"
#include "../entities/Universe.h"
#include "Correlator.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include "OutputSupportingSimulation.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
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

    class DPDSimulator : public pylimer_tools::calc::OutputSupportingSimulation
    {

    private:
      ////////////////////////////////////////////////////////////////
      // configuration
      double maxBondLen = 5.;
      bool is2D = false;
      bool shiftPossibilityEmpty = true;
      bool shiftOneAtATime = false;
      double lambda = 0.65;
      double k = 2.;
      double lowCutoff = 0.5;
      double highCutoff = 2.0;
      double A = 25.;
      double sigma = 3.;
      double gamma = 0.5 * 3. * 3.;
      long int nStepsMC = 500;
      long int nStepsDPD = 500;

      ////////////////////////////////////////////////////////////////
      // simulation state
      long int currentStep = 0;
      double currentTime = 0.;
      Eigen::VectorXd currentVelocitiesPlus;
      Eigen::VectorXd currentVelocities;
      Eigen::VectorXd currentForces;
      Eigen::Matrix3d currentStressTensor;

      ////////////////////////////////////////////////////////////////
      // randomness
      std::mt19937 e2;
      std::uniform_real_distribution<double> uniform_rand_mean0std1;
      std::uniform_real_distribution<double> uniform_rand_between_0_1;

      ////////////////////////////////////////////////////////////////
      // universe structure
      int numAtoms = 0;
      int numBonds = 0;
      int numSlipSprings = 0;
      pylimer_tools::entities::Box box;
      pylimer_tools::entities::Universe universe;

      // atoms
      Eigen::VectorXd coordinates;
      Eigen::ArrayXi idxFunctionalities;
      std::vector<int> atomTypes;
      std::vector<long int> atomIds;
      std::vector<size_t> chainEndIndices;
      // bonds
      Eigen::ArrayXi bondPartnerCoordinatesA;
      Eigen::ArrayXi bondPartnersA;
      Eigen::ArrayXi bondPartnerCoordinatesB;
      Eigen::ArrayXi bondPartnersB;
      Eigen::ArrayXi bondTypes;
      // mapping between atoms and bonds
      std::vector<std::vector<size_t>> bondsOfIndex;

      pylimer_tools::entities::EigenNeighbourList neighbourlist;

    public:
      DPDSimulator(const pylimer_tools::entities::Universe& u,
                   const int crosslinkerType = 2,
                   const bool is2D = false,
                   const std::string& seed = "");

      /**
       * @brief actually do run the simulation
       *
       */
      void runSimulation(const long int nSteps,
                         double dt,
                         bool withMC,
                         const std::function<bool()>& shouldInterrupt,
                         const std::function<void()>& cleanupInterrupt);

      void runSimulation(const long int nSteps,
                         double dt = 0.06,
                         bool withMC = false)
      {
        runSimulation(
          nSteps, dt, withMC, []() { return false; }, []() {});
      }

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
       * @brief Compute the temperature
       *
       * @param velocities
       * @return double
       */
      double getTemperature() const
      {
        return this->computeTemperature(this->currentVelocities);
      };

      /**
       * @brief Get access to the current stress-tensor
       *
       * @return Eigen::Matrix3d
       */
      Eigen::Matrix3d getStressTensor() const;

      double computeBondLength(int bondIdx) const
      {
        Eigen::Vector3d bondDistances =
          this->coordinates(
            this->bondPartnerCoordinatesA.segment(3 * bondIdx, 3)) -
          this->coordinates(
            this->bondPartnerCoordinatesB.segment(3 * bondIdx, 3));
        this->box.handlePBC(bondDistances);
        return bondDistances.norm();
      };

      long int getTimestep() const;

      ////////////////////////////////////////////////////////////////
      // MC Procedures
      /**
       * @brief Randomly add new slip-springs
       */
      int createSlipSprings(const int num, const int bondType);

      int shiftSlipSprings(const double kbT = 1.);

      int relocateSlipSprings(const double kbT = 1.);

      ////////////////////////////////////////////////////////////////
      // configuration
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

      void configSigma(const double newSigma)
      {
        this->sigma = newSigma;
        this->gamma = 0.5 * newSigma * newSigma;
      }

      void configA(const double newA) { this->A = newA; }

      void startMeasuringMSDForAtoms(const std::vector<size_t>& atomIds);

      void configNumStepsMC(long int steps = 500) { this->nStepsMC = steps; }

      void configNumStepsDPD(long int steps = 500) { this->nStepsDPD = steps; }

      void configShiftPossibilityEmpty(bool shiftPossibilityEmptyConfig = true)
      {
        this->shiftPossibilityEmpty = shiftPossibilityEmptyConfig;
      }

      void configShiftOneAtATime(bool shiftOne = false)
      {
        this->shiftOneAtATime = shiftOne;
      }

      ////////////////////////////////////////////////////////////////
      // results access & export

      pylimer_tools::entities::Universe getUniverse() const;

      ////////////////////////////////////////////////////////////////
      // validation
      void validateState();
      void validateNeighbourlist(double cutoff);

    protected:
      void addSlipSprings(std::vector<size_t>& partnerA,
                          std::vector<size_t>& partnerB,
                          const int bondType = 9);

      bool attemptSlipSpringShift(const size_t springIdx,
                                  const size_t endToShift,
                                  const double kbT = 1.);

      bool attemptSlipSpringShift(const size_t springIdx,
                                  const double kbT = 1.);

      void replaceSlipSpringPartner(const size_t springIdx,
                                    const size_t partnerBefore,
                                    const size_t partnerAfter);
    };
  }

}
}

#endif
