#ifndef DPD_SIMULATOR_H
#define DPD_SIMULATOR_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/EigenNeighbourList.h"
#include "../entities/Universe.h"
#include "../utils/ExtraEigenTypes.h"
#include "../utils/PerformanceTimer.h"
#include "Correlator.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include "OutputSupportingSimulation.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
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
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools {
namespace calc {

  namespace dpd {

    enum DPDPerformanceSections
    {
      TIME_STEPPING,
      FORCES,
      PAIR_FORCE,
      BOND_FORCE,
      OUTPUT,
      SHIFT,
      RELOCATION,
      NUM_PERFORMANCE_SECTIONS
    };

    static const std::array<std::string, 7> DPDPerformanceSectionNames = {
      "Time-Stepping", "Forces", "Pair-Forces", "Bond-Forces",
      "Output",        "Shift",  "Relocation"
    };

    class DPDSimulator : public pylimer_tools::calc::OutputSupportingSimulation
    {

    private:
      DPDSimulator(){}; // not exposed to users, only used by Cereal
      friend class cereal::access;

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
      double dt = 0.06;

      ////////////////////////////////////////////////////////////////
      // simulation state
      long int currentStep = 0;
      double currentTime = 0.;
      long int numShifts = 0;
      long int numRelocations = 0;
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
                         bool withMC,
                         const std::function<bool()>& shouldInterrupt,
                         const std::function<void()>& cleanupInterrupt);

      void runSimulation(const long int nSteps, bool withMC = false)
      {
        runSimulation(
          nSteps, withMC, []() { return false; }, []() {});
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
        pylimer_tools::utils::PerformanceTimer<
          DPDPerformanceSections::NUM_PERFORMANCE_SECTIONS>& timer,
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
      void configTimeStep(const double dt = 0.06) { this->dt = dt; }

      void configLambda(const double l) { this->lambda = l; }

      void configSpringConstant(const double nk) { this->k = nk; }

      double getSpringConstant() const { return this->k; }

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

      long int getNumStepsMC() { return this->nStepsMC; }

      void configNumStepsDPD(long int steps = 500) { this->nStepsDPD = steps; }

      long int getNumStepsDPD() { return this->nStepsDPD; }

      void configShiftPossibilityEmpty(bool shiftPossibilityEmptyConfig = true)
      {
        this->shiftPossibilityEmpty = shiftPossibilityEmptyConfig;
      }
      bool getShiftPossibilityEmpty() const
      {
        return this->shiftPossibilityEmpty;
      }

      void configShiftOneAtATime(bool shiftOne = false)
      {
        this->shiftOneAtATime = shiftOne;
      }

      bool getShiftOneAtATime() const { return this->shiftOneAtATime; }

      ////////////////////////////////////////////////////////////////
      // results access & export

      pylimer_tools::entities::Universe getUniverse() const;

      double getTimestep() override { return this->dt; }
      double getCurrentTime(double currentStep) override
      {
        return this->currentTime;
      }

      /**
       * @brief Get access to the current stress-tensor
       *
       * @return Eigen::Matrix3d
       */
      Eigen::Matrix3d getStressTensor() override
      {
        return this->currentStressTensor;
      }

      int getNumShifts() override { return this->numShifts; }

      int getNumRelocations() override { return this->numRelocations; }

      size_t getNumParticles() override { return this->numAtoms; }

      size_t getNumSlipSprings() const { return this->numSlipSprings; }

      double getVolume() override { return this->box.getVolume(); }

      Eigen::VectorXd getBondLengths() override
      {
        Eigen::VectorXd bondDistances =
          this->coordinates(this->bondPartnerCoordinatesA) -
          this->coordinates(this->bondPartnerCoordinatesB);
        this->box.handlePBC(bondDistances);

        Eigen::VectorXd bondLengths =
          Eigen::VectorXd::Zero(this->numBonds + this->numSlipSprings);
#pragma omp parallel for
        for (size_t i = 0; i < (this->numBonds + this->numSlipSprings); ++i) {
          double b = bondDistances.segment(3 * i, 3).norm();
          bondLengths[i] = b;
        }
        return bondLengths;
      }

      Eigen::VectorXd getCoordinates() override { return this->coordinates; }

      double getTemperature() override
      {
        return this->computeTemperature(this->currentVelocities);
      }

      ////////////////////////////////////////////////////////////////
      // validation
      void validateState();
      void validateNeighbourlist(double cutoff);

      ////////////////////////////////////////////////////////////////
      // serialization
      template<class Archive>
      void serialize(Archive& ar)
      {
        ar(cereal::virtual_base_class<OutputSupportingSimulation>(this),
           // configuration
           maxBondLen,
           is2D,
           shiftPossibilityEmpty,
           shiftOneAtATime,
           lambda,
           k,
           lowCutoff,
           highCutoff,
           A,
           sigma,
           gamma,
           nStepsDPD,
           nStepsMC,
           dt,
           // simulation state
           currentStep,
           currentTime,
           numShifts,
           numRelocations,
           currentVelocitiesPlus,
           currentVelocities,
           currentForces,
           currentStressTensor,
           // randomness
           e2,
           uniform_rand_mean0std1,
           uniform_rand_between_0_1,
           // universe structure
           numAtoms,
           numBonds,
           numSlipSprings,
           box,
           universe,
           // -> atoms
           coordinates,
           idxFunctionalities,
           atomTypes,
           atomIds,
           chainEndIndices,
           // -> bonds
           bondPartnerCoordinatesA,
           bondPartnerCoordinatesB,
           bondPartnersA,
           bondPartnersB,
           bondTypes,
           bondsOfIndex,
           // neighbourlist
           neighbourlist);
      }

      static DPDSimulator readRestartFile(std::string filename)
      {
        DPDSimulator res;
        pylimer_tools::utils::deserializeFromFile<DPDSimulator>(res, filename);
        return res;
      };

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

      void writeRestartFile(std::string& filename) override
      {
        pylimer_tools::utils::serializeToFile<DPDSimulator>(*this, filename);
      };
    };
  };
}
}

CEREAL_REGISTER_TYPE(pylimer_tools::calc::dpd::DPDSimulator);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
  pylimer_tools::calc::OutputSupportingSimulation,
  pylimer_tools::calc::dpd::DPDSimulator);

#endif
