#ifndef DPD_SIMULATOR_H
#define DPD_SIMULATOR_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/EigenNeighbourList.h"
#include "../entities/Universe.h"
#include "Correlator.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
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

    enum ComputedIntValues
    {
      STEP = 0,
      NUM_SHIFT = 1,
      NUM_RELOC = 2,
    };

    const std::array<std::string, 3> ComputedIntValuesNames = {
      "Step",
      "numShift",
      "numReloc",
    };

    enum ComputedDoubleValues
    {
      TIMESTEP = 0,
      TIME = 1,
      VOLUME = 2,
      PRESSURE = 3,
      TEMPERATURE = 4,
      STRESS_XX = 5,
      STRESS_YY = 6,
      STRESS_ZZ = 7,
      STRESS_XY = 8,
      STRESS_YZ = 9,
      STRESS_XZ = 10,
      STRESS_NXY = 11,
      STRESS_NYZ = 12,
      STRESS_NXZ = 13,
      MEAN_B = 14,
      MAX_B = 15,
      MSD = 16
    };

    const std::array<std::string, 17> ComputedDoubleValuesNames = {
      "TimeStep",
      "Time",
      "Volume",
      "Pressure",
      "Temperature",
      "Stress[0,0]",
      "Stress[1,1]",
      "Stress[2,2]",
      "Stress[0,1]",
      "Stress[1,2]",
      "Stress[0,2]",
      "Stress[0,0]-Stress[1,1]",
      "Stress[1,1]-Stress[2,2]",
      "Stress[0,0]-Stress[2,2]",
      "<b>",
      "max(b)",
      "MSD"
    };

    struct OutputConfiguration
    {
      std::vector<ComputedIntValues> intValues;
      std::vector<ComputedDoubleValues> doubleValues;
      std::string filename;
      int outputEvery;

      OutputConfiguration()
        : filename("")
        , outputEvery(10)
      {
      }
    };

    typedef Eigen::Array<size_t, Eigen::Dynamic, 1> ArrayXst;
    typedef Eigen::Array<long int, Eigen::Dynamic, 1> ArrayXli;

    class DPDSimulator
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
      std::vector<OutputConfiguration> outputConfigs;
      std::vector<OutputConfiguration> outputAverageConfigs;
      std::vector<OutputConfiguration> outputAutoCorrelationConfigs;
      std::vector<pylimer_tools::calc::Correlator> autocorrelators;
      std::vector<std::shared_ptr<std::ostream>> outputStreams;
      std::vector<int> outputFileStreams;

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

      ////////////////////////////////////////////////////////////////
      // computation state
      std::vector<Eigen::ArrayXi> msdMeasuredIndices;
      std::vector<Eigen::VectorXd> msdOrigins;
      std::vector<size_t> msdOriginTimesteps;

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

      void startMeasuringMSDForAtoms(const std::vector<size_t> atomIds);

      void configAutoCorrelatorOutput(std::vector<OutputConfiguration> vals,
                                      const unsigned int numcorrin = 32,
                                      const unsigned int pin = 16,
                                      const unsigned int min = 2)
      {
        int num_values_to_correlate = 0;
        for (size_t i = 0; i < vals.size(); ++i) {
          INVALIDARG_EXP_IFN(
            vals[i].intValues.size() == 0,
            "Correlation of integer values is not supported yet.");
          num_values_to_correlate += vals[i].doubleValues.size();
        }
        this->autocorrelators.clear();
        this->autocorrelators.reserve(num_values_to_correlate);
        for (size_t i = 0; i < num_values_to_correlate; ++i) {
          pylimer_tools::calc::Correlator correlator =
            pylimer_tools::calc::Correlator(numcorrin, pin, min);
          this->autocorrelators.push_back(correlator);
        }
        this->outputAutoCorrelationConfigs = vals;
      }

      void configAverageOutput(std::vector<OutputConfiguration> vals)
      {
        this->outputAverageConfigs = vals;
      }

      void configStepOutput(std::vector<OutputConfiguration> vals)
      {
        this->outputConfigs = vals;
      }

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
      int openFilesOutputHeader(std::vector<OutputConfiguration>& configs,
                                std::string prefix = "",
                                int streamIdx = 0);
      inline void doOutputValues(OutputConfiguration& oc,
                                 std::array<int, 3>& intvalues,
                                 std::array<double, 17>& doublevalues,
                                 std::string& outputBuffer,
                                 int streamIdx = 0)
      {
        assert(streamIdx <= this->outputStreams.size());
        for (ComputedIntValues val : oc.intValues) {
          switch (val) {
            default:
              outputBuffer += std::to_string(intvalues[val]) + "\t";
          }
        }
        for (ComputedDoubleValues val : oc.doubleValues) {
          switch (val) {
            case ComputedDoubleValues::MSD:
              // compute MSD
              for (size_t msdIdx = 0; msdIdx < msdMeasuredIndices.size();
                   ++msdIdx) {
                double result =
                  (this->msdOrigins[msdIdx] -
                   this->coordinates(this->msdMeasuredIndices[msdIdx]))
                    .squaredNorm() /
                  (static_cast<double>(this->msdMeasuredIndices[msdIdx].size() /
                                       3.));
                outputBuffer += std::to_string(result) + "\t";
              }
              break;
            default:
              outputBuffer += std::to_string(doublevalues[val]) + "\t";
          }
        }
        if (!outputBuffer.empty()) {
          outputBuffer.pop_back(); // remove last "\t"
          outputBuffer += "\n";
          // output the buffer, clear it
          (*(this->outputStreams[streamIdx])) << outputBuffer;
          outputBuffer.clear();
        }
      };

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
