#ifndef OUTPUT_SUPPORTING_SIM_H
#define OUTPUT_SUPPORTING_SIM_H

#include "../calc/Correlator.h"
#include "../utils/CerealUtils.h"
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#ifdef CEREALIZABLE
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

namespace pylimer_tools {
namespace sim {

#define NUM_COMPUTABLE_INT_VALUES 8

  enum ComputedIntValues
  {
    STEP = 0,
    NUM_SHIFT = 1,
    NUM_RELOC = 2,
    NUM_ATOMS = 3,
    NUM_EXTRA_ATOMS = 4,
    NUM_BONDS = 5,
    NUM_EXTRA_BONDS = 6,
    NUM_BONDS_TO_FORM = 7,
  };

  const std::array<std::string, NUM_COMPUTABLE_INT_VALUES>
    ComputedIntValuesNames = {
      "Step",          "numShift", "numReloc",      "numAtoms",
      "numExtraAtoms", "numBonds", "numExtraBonds", "numBondsToForm",
    };

#define NUM_COMPUTABLE_DOUBLE_VALUES 19

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
    GAMMA = 14,
    RESIDUAL = 15,
    MEAN_B = 16,
    MAX_B = 17,
    MSD = 18
  };

  const std::array<std::string, NUM_COMPUTABLE_DOUBLE_VALUES>
    ComputedDoubleValuesNames = { "TimeStep",
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
                                  "Gamma",
                                  "Residual",
                                  "<b>",
                                  "max(b)",
                                  "MSD" };

  struct OutputConfiguration
  {
    std::vector<ComputedIntValues> intValues;
    std::vector<ComputedDoubleValues> doubleValues;
    std::string filename = "";
    int outputEvery = 10;
    /**
     * @brief use every: for autocorrelation/averaging, how often to include
     * values
     */
    int useEvery = 1;
    /**
     * @brief Whether to append to the file or truncate it
     *
     * This does not need to be persisted, as when restarting, the output will
     * append anyway.
     */
    bool append = false;

    OutputConfiguration() {}

    template<class Archive>
    void serialize(Archive& ar)
    {
      ar(intValues, doubleValues, filename, outputEvery, useEvery);
    }
  };

  class OutputSupportingSimulation
  {
  protected:
    ////////////////////////////////////////////////////////////////
    // output configurations
    std::vector<OutputConfiguration> outputConfigs = {};
    std::vector<OutputConfiguration> outputAverageConfigs = {};
    std::vector<OutputConfiguration> outputAutoCorrelationConfigs = {};
    ////////////////////////////////////////////////////////////////
    // restart configurations
    int outputRestartEvery = 0;
    std::string restartOutputFile = "";
    ////////////////////////////////////////////////////////////////
    // output speedups
    std::array<int, NUM_COMPUTABLE_DOUBLE_VALUES> doubleValueRequiredEvery;
    std::array<int, NUM_COMPUTABLE_INT_VALUES> intValueRequiredEvery;
    int requireStressTensorEvery = 0;
    int requireBondLenEvery = 0;
    std::string outputBuffer = "";
    bool doAverage = false;
    ////////////////////////////////////////////////////////////////
    // output streams
    std::vector<std::shared_ptr<std::ostream>> outputStreams = {};
    std::vector<int> outputFileStreams = {};
    ////////////////////////////////////////////////////////////////
    // computation state
    std::vector<pylimer_tools::calc::Correlator> autocorrelators = {};
    std::vector<Eigen::ArrayXi> msdMeasuredIndices = {};
    std::vector<Eigen::VectorXd> msdOrigins = {};
    std::vector<size_t> msdOriginTimesteps = {};
    std::vector<double> runningAverages = {};

    /**
     * @brief Open all files required for output
     *
     */
    void prepareAllOutputs()
    {
      this->outputStreams.clear();
      this->outputFileStreams.clear();

      // output headers
      std::cout.flush();
      std::ios::sync_with_stdio(false);
      this->openFilesOutputHeader(this->outputConfigs);

      // prepare averages
      int numAverages = this->openFilesOutputHeader(this->outputAverageConfigs,
                                                    "# OutputStep\t",
                                                    this->outputConfigs.size());
      RUNTIME_EXP_IFN(runningAverages.size() == numAverages,
                      "The nr. of running averages is not consistent with the "
                      "number of output quantities.");

      // prepare autocorrelation
      this->openFilesOutputHeader(this->outputAutoCorrelationConfigs,
                                  "Step\t",
                                  this->outputConfigs.size() +
                                    this->outputAverageConfigs.size());
      std::string autocorrelationOutputBuffer;
      autocorrelationOutputBuffer.reserve(
        this->outputAutoCorrelationConfigs.size() * 50);
    }

    /**
     * @brief Close all files required for output
     *
     */
    void closeAllOutputs()
    {
      // finish up
      std::ios::sync_with_stdio(true);
      std::cout.flush();
      this->outputStreams.clear();
      this->outputFileStreams.clear();
    }

    int openFilesOutputHeader(const std::vector<OutputConfiguration>& configs,
                              const std::string& prefix = "",
                              int streamIdx = 0);

    inline bool requiresDEvaluation(const ComputedDoubleValues val,
                                    const long int currentStep) const
    {
      return (this->doubleValueRequiredEvery[val] != 0) &&
             ((currentStep % this->doubleValueRequiredEvery[val]) == 0);
    }

    inline bool requiresIEvaluation(const ComputedIntValues val,
                                    const long int currentStep) const
    {
      return (this->intValueRequiredEvery[val] != 0) &&
             ((currentStep % this->intValueRequiredEvery[val]) == 0);
    }

    /**
     * @brief Iterate all possible output configurations, handle them
     *
     * @param current_step
     */
    void handleOutput(const long int currentStep)
    {
      // "lazily" compute the values we need, others less lazily when they are
      // computationally inexpensive
      std::array<long int, NUM_COMPUTABLE_INT_VALUES> intvalues = {
        currentStep,
        this->requiresIEvaluation(NUM_SHIFT, currentStep) ? this->getNumShifts()
                                                          : 0,
        this->requiresIEvaluation(NUM_RELOC, currentStep)
          ? this->getNumRelocations()
          : 0,
        this->requiresIEvaluation(NUM_ATOMS, currentStep)
          ? static_cast<long int>(this->getNumAtoms())
          : 0,
        this->requiresIEvaluation(NUM_EXTRA_ATOMS, currentStep)
          ? static_cast<long int>(this->getNumExtraAtoms())
          : 0,
        this->requiresIEvaluation(NUM_BONDS, currentStep)
          ? static_cast<long int>(this->getNumBonds())
          : 0,
        this->requiresIEvaluation(NUM_EXTRA_BONDS, currentStep)
          ? static_cast<long int>(this->getNumExtraBonds())
          : 0,
        this->requiresIEvaluation(NUM_BONDS_TO_FORM, currentStep)
          ? this->getNumBondsToForm()
          : 0,
      };

      Eigen::Matrix3d stressTensor =
        ((this->requireStressTensorEvery > 0) &&
         (currentStep % this->requireStressTensorEvery) == 0)
          ? this->getStressTensor()
          : Eigen::Matrix3d::Zero();
      double pressure = stressTensor.trace() / 3.;
      // double kineticPressureTerm =
      //   requiresDEvaluation(PRESSURE, currentStep)
      //     ? ((getNumParticles() * this->getTemperature()) /
      //     this->getVolume()) : 0.0;
      Eigen::VectorXd bondLengths =
        (((this->requireBondLenEvery > 0)) &&
         ((currentStep % this->requireBondLenEvery) == 0))
          ? this->getBondLengths()
          : Eigen::Vector3d::Zero();
      if (bondLengths.size() == 0) {
        bondLengths = Eigen::Vector3d::Zero();
      }

      // assemble all computed values into an easy-to-access array
      std::array<double, NUM_COMPUTABLE_DOUBLE_VALUES> doublevalues = {
        this->getTimestep(),
        this->requiresDEvaluation(TIME, currentStep)
          ? this->getCurrentTime(currentStep)
          : 0.,
        this->requiresDEvaluation(VOLUME, currentStep) ? this->getVolume() : 0.,
        pressure, // + kineticPressureTerm,
        this->requiresDEvaluation(TEMPERATURE, currentStep)
          ? this->getTemperature()
          : 0.,
        stressTensor(0, 0),
        stressTensor(1, 1),
        stressTensor(2, 2),
        stressTensor(0, 1),
        stressTensor(1, 2),
        stressTensor(0, 2),
        stressTensor(0, 0) - stressTensor(1, 1),
        stressTensor(1, 1) - stressTensor(2, 2),
        stressTensor(0, 0) - stressTensor(2, 2),
        this->requiresDEvaluation(GAMMA, currentStep) ? this->getGamma() : 0.,
        this->requiresDEvaluation(RESIDUAL, currentStep) ? this->getResidual() : 0.,
        this->requiresDEvaluation(MEAN_B, currentStep) ? bondLengths.mean()
                                                       : 0.0,
        this->requiresDEvaluation(MAX_B, currentStep) ? bondLengths.maxCoeff()
                                                      : 0.0,
        0.
      };
      int streamIdx = 0;
      for (streamIdx = 0; streamIdx < this->outputConfigs.size(); ++streamIdx) {
        if (currentStep % this->outputConfigs[streamIdx].outputEvery == 0) {
          this->doOutputValues(
            this->outputConfigs[streamIdx], intvalues, doublevalues, streamIdx);
          outputBuffer.clear();
        }
      }

      // compute averages
      if (doAverage) {
        size_t msdIdx = 0;
        size_t averagesIdx = 0;
        for (const OutputConfiguration& oc : this->outputAverageConfigs) {
          size_t previousAverageIdx = averagesIdx;
          if ((currentStep % oc.useEvery) == 0) {
            double multiplier = (static_cast<double>(oc.useEvery) /
                                 static_cast<double>(oc.outputEvery));
            for (ComputedIntValues val : oc.intValues) {
              switch (val) {
                default:
                  runningAverages[averagesIdx] +=
                    static_cast<double>(intvalues[val]) * multiplier;
                  averagesIdx += 1;
                  break;
              }
            }
            for (ComputedDoubleValues val : oc.doubleValues) {
              switch (val) {
                case ComputedDoubleValues::MSD:
                  // compute MSD
                  for (msdIdx = 0; msdIdx < this->msdMeasuredIndices.size();
                       ++msdIdx) {
                    double result =
                      (this->msdOrigins[msdIdx] -
                       getCoordinates()(this->msdMeasuredIndices[msdIdx]))
                        .squaredNorm() /
                      (static_cast<double>(
                        this->msdMeasuredIndices[msdIdx].size() / 3.));
                    runningAverages[averagesIdx + msdIdx] +=
                      result * multiplier;
                  }
                  averagesIdx += msdIdx;
                  break;
                default:
                  runningAverages[averagesIdx] +=
                    doublevalues[val] * multiplier;
                  averagesIdx += 1;
                  break;
              }
            }
          }

          // check (and if, output) averages
          if (currentStep % oc.outputEvery == 0) {
            // output & start again
            outputBuffer += std::to_string(intvalues[ComputedIntValues::STEP]);
            for (size_t i = previousAverageIdx; i < averagesIdx; ++i) {
              outputBuffer += "\t" + std::to_string(runningAverages[i]);
              runningAverages[i] = 0.;
            }
            (*(this->outputStreams[streamIdx])) << outputBuffer << std::endl;
            outputBuffer.clear();
          }

          streamIdx += 1;
        }
      }

      // do autocorrelation
      size_t autocorrelator_idx = 0;
      for (const OutputConfiguration& oc : this->outputAutoCorrelationConfigs) {
        const size_t autocorrelator_idx_before = autocorrelator_idx;
        for (ComputedDoubleValues cv : oc.doubleValues) {
          assert(autocorrelator_idx < this->autocorrelators.size());
          RUNTIME_EXP_IFN(std::isfinite(doublevalues[cv]),
                          "Expect output quantities to be finite, found " +
                            std::to_string(doublevalues[cv]) +
                            " for property " + ComputedDoubleValuesNames[cv] +
                            ".");
          this->autocorrelators[autocorrelator_idx].add(doublevalues[cv]);
          autocorrelator_idx += 1;
        }
        if (currentStep % oc.outputEvery == 0) {
          outputBuffer += "# TimeStep " +
                          std::to_string(intvalues[ComputedIntValues::STEP]) +
                          "\n";
          this->autocorrelators[autocorrelator_idx_before].evaluate();
          const unsigned int npcorr =
            this->autocorrelators[autocorrelator_idx_before].npcorr;
          RUNTIME_EXP_IFN(npcorr > 0,
                          "Expected more than 0 correlator results.");
          for (int autocorr_idx_offset = 1;
               autocorr_idx_offset < oc.doubleValues.size();
               ++autocorr_idx_offset) {
            size_t idx = autocorrelator_idx_before + autocorr_idx_offset;
            this->autocorrelators[idx].evaluate();
            RUNTIME_EXP_IFN(this->autocorrelators[idx].npcorr == npcorr,
                            "Autocorrelation states are inconsistent.");
          }

          for (size_t output_idx = 0; output_idx < npcorr; output_idx += 1) {
            outputBuffer += std::to_string(
              this->autocorrelators[autocorrelator_idx_before].t[output_idx]);
            for (int autocorr_idx_offset = 0;
                 autocorr_idx_offset < oc.doubleValues.size();
                 ++autocorr_idx_offset) {
              size_t idx = autocorrelator_idx_before + autocorr_idx_offset;
              outputBuffer +=
                "\t" + std::to_string(this->autocorrelators[idx].f[output_idx]);
            }
            outputBuffer += "\n";
          }
          (*(this->outputStreams[streamIdx])) << outputBuffer << std::flush;
          streamIdx += 1;
          outputBuffer.clear();
        }

        streamIdx += 1;
      }

      // potentially write restart file
#ifdef CEREALIZABLE
      if (this->outputRestartEvery > 0 &&
          currentStep % this->outputRestartEvery == 0) {
        this->writeRestartFile(this->restartOutputFile);
      }
#endif

      if (currentStep % 50 == 0) {
        std::flush(std::cout);
      }
    }

    // static OutputSupportingSimulation readRestartFile(std::string filename)
    // {
    //   throw std::runtime_error(
    //     "Cannot read restart file on abstract base class.");
    // };

    /**
     * @brief Output the passed valus
     *
     * @param oc
     * @param intvalues
     * @param doublevalues
     * @param outputBuffer
     * @param coordinates
     * @param streamIdx
     */
    inline void doOutputValues(
      const OutputConfiguration& oc,
      const std::array<long int, NUM_COMPUTABLE_INT_VALUES>& intvalues,
      const std::array<double, NUM_COMPUTABLE_DOUBLE_VALUES>& doublevalues,
      int streamIdx = 0)
    {
      assert(streamIdx <= this->outputStreams.size());
      for (ComputedIntValues val : oc.intValues) {
        RUNTIME_EXP_IFN(std::isfinite(static_cast<double>(intvalues[val])),
                        "Expect output quantities to be finite, found " +
                          std::to_string(intvalues[val]) + " for property " +
                          ComputedIntValuesNames[val] + ".");
        switch (val) {
          default:
            outputBuffer += std::to_string(intvalues[val]) + "\t";
        }
      }
      for (ComputedDoubleValues val : oc.doubleValues) {
        RUNTIME_EXP_IFN(std::isfinite(doublevalues[val]),
                        "Expect output quantities to be finite, found " +
                          std::to_string(doublevalues[val]) + " for property " +
                          ComputedDoubleValuesNames[val] + ".");
        switch (val) {
          case ComputedDoubleValues::MSD:
            // compute MSD
            for (size_t msdIdx = 0; msdIdx < this->msdMeasuredIndices.size();
                 ++msdIdx) {
              assert(this->msdOrigins.size() > msdIdx);
              assert(this->getCoordinates().size() >
                     this->msdMeasuredIndices[msdIdx].maxCoeff());
              assert(this->msdOrigins[msdIdx].size() ==
                     this->msdMeasuredIndices[msdIdx].size());
              Eigen::ArrayXi nIndices = this->msdMeasuredIndices[msdIdx];
              assert(nIndices.size() > 0);
              assert(nIndices.size() == this->msdOrigins[msdIdx].size());
              Eigen::VectorXd relCoords = this->getCoordinates();
              assert(nIndices.minCoeff() >= 0 &&
                     nIndices.maxCoeff() < relCoords.size());
              double result =
                (this->msdOrigins[msdIdx] - relCoords(nIndices)).squaredNorm() /
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

    /**
     * @brief Remember how often a particular value is needed to be computed for
     * any of the outputs
     *
     * @param configs
     */
    void updateValuesRequiredEvery(
      const std::vector<OutputConfiguration>& configs)
    {
      for (OutputConfiguration c : configs) {
        for (ComputedDoubleValues v : c.doubleValues) {
          if (this->doubleValueRequiredEvery[v] == 0) {
            this->doubleValueRequiredEvery[v] = c.useEvery;
          } else {
            this->doubleValueRequiredEvery[v] =
              std::gcd(c.useEvery, this->doubleValueRequiredEvery[v]);
          }
        }
        for (ComputedIntValues i : c.intValues) {
          if (this->intValueRequiredEvery[i] == 0) {
            this->intValueRequiredEvery[i] = c.useEvery;
          } else {
            this->intValueRequiredEvery[i] =
              std::gcd(c.useEvery, this->intValueRequiredEvery[i]);
          }
        }
      }
      std::vector<ComputedDoubleValues> stressTensorRequiringValues = {
        STRESS_XX, STRESS_YY,  STRESS_ZZ,  STRESS_XY,  STRESS_YZ,
        STRESS_XZ, STRESS_NXY, STRESS_NYZ, STRESS_NXZ, PRESSURE
      };
      for (ComputedDoubleValues v : stressTensorRequiringValues) {
        if (requireStressTensorEvery == 0) {
          requireStressTensorEvery = this->doubleValueRequiredEvery[v];
        } else {
          requireStressTensorEvery = std::gcd(
            requireStressTensorEvery, this->doubleValueRequiredEvery[v]);
        }
      }
      // similarly for the bond length
      requireBondLenEvery =
        std::gcd(this->doubleValueRequiredEvery[ComputedDoubleValues::MAX_B],
                 this->doubleValueRequiredEvery[ComputedDoubleValues::MEAN_B]);
    }

  public:
    OutputSupportingSimulation()
    {
      this->doubleValueRequiredEvery.fill(0);
      this->intValueRequiredEvery.fill(0);
      this->outputBuffer.reserve(600);
    }

#ifdef CEREALIZABLE
    virtual void writeRestartFile(std::string& filename) = 0;
#endif

    void validateAndTruncateOutputFiles(
      const std::vector<OutputConfiguration>& vals) const
    {
      for (size_t i = 0; i < vals.size(); ++i) {
        if (vals[i].filename.size() > 0) {
          // empty the file
          std::ifstream file;
          file.open(
            vals[i].filename.c_str(),
            std::ifstream::out |
              (vals[i].append ? std::ifstream::app : std::ifstream::trunc));
          if (!file.is_open() || file.fail()) {
            file.close();
            throw std::invalid_argument("The file " + vals[i].filename +
                                        " could not be opened.");
          }
          file.close();
        }
      }
    }

    void configAutoCorrelatorOutput(std::vector<OutputConfiguration>& vals,
                                    const unsigned int numcorrin = 32,
                                    const unsigned int pin = 16,
                                    const unsigned int min = 2)
    {
      int num_values_to_correlate = 0;
      for (size_t i = 0; i < vals.size(); ++i) {
        INVALIDARG_EXP_IFN(
          vals[i].intValues.size() == 0,
          "Correlation of integer values is not supported yet.");
        INVALIDARG_EXP_IFN(vals[i].outputEvery >= vals[i].useEvery,
                           "Require useEvery to be smaller than output every");
        num_values_to_correlate += vals[i].doubleValues.size();
      }
      this->validateAndTruncateOutputFiles(vals);
      this->autocorrelators.clear();
      this->autocorrelators.reserve(num_values_to_correlate);
      this->updateValuesRequiredEvery(vals);
      for (size_t i = 0; i < num_values_to_correlate; ++i) {
        pylimer_tools::calc::Correlator correlator =
          pylimer_tools::calc::Correlator(numcorrin, pin, min);
        this->autocorrelators.push_back(correlator);
      }
      this->outputAutoCorrelationConfigs = vals;
    }

    void configAverageOutput(const std::vector<OutputConfiguration>& configs)
    {
      this->outputAverageConfigs = configs;
      this->updateValuesRequiredEvery(configs);

      int numAverages = 0;
      for (const OutputConfiguration& c : configs) {
        numAverages += c.doubleValues.size();
        numAverages += c.intValues.size();
        INVALIDARG_EXP_IFN(c.outputEvery >= c.useEvery,
                           "Require useEvery to be smaller than output every");
        INVALIDARG_EXP_IFN(c.outputEvery % c.useEvery == 0,
                           "Output every must be a multiple of useEvery");
      }

      this->validateAndTruncateOutputFiles(configs);

      this->runningAverages =
        pylimer_tools::utils::initializeWithValue<double>(numAverages, 0.);
      this->doAverage = numAverages > 0;
    }

    void configStepOutput(std::vector<OutputConfiguration>& vals)
    {
      for (size_t i = 0; i < vals.size(); ++i) {
        vals[i].useEvery = vals[i].outputEvery;
      }

      this->validateAndTruncateOutputFiles(vals);
      this->outputConfigs = vals;
      this->updateValuesRequiredEvery(vals);
    }

    void configRestartOutput(const std::string outputFile, int outputEvery)
    {
      this->outputRestartEvery = outputEvery;
      this->restartOutputFile = outputFile;
    }

#ifdef CEREALIZABLE
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version)
    {
      ar(
        // output configurations
        outputConfigs,
        outputAverageConfigs,
        outputAutoCorrelationConfigs);

      ar(
        // restart configurations - meta!
        outputRestartEvery,
        restartOutputFile,
        // output speedups – could also recompute instead ?!?
        doubleValueRequiredEvery,
        intValueRequiredEvery,
        requireStressTensorEvery,
        requireBondLenEvery,
        outputBuffer,
        doAverage,
        // output streams – here, it gets dangerous!
        outputFileStreams,
        // outputStreams,
        // computation state
        autocorrelators,
        msdMeasuredIndices,
        msdOrigins,
        msdOriginTimesteps,
        runningAverages);
    }
#endif

    virtual double getCurrentTime(double currentStep) = 0;
    virtual double getGamma() = 0;
    virtual double getResidual() = 0;
    virtual double getTemperature() = 0;
    virtual double getTimestep() = 0;
    virtual double getVolume() = 0;
    virtual Eigen::Matrix3d getStressTensor() = 0;
    virtual Eigen::VectorXd getBondLengths() = 0;
    virtual Eigen::VectorXd getCoordinates() = 0;
    virtual int getNumRelocations() = 0;
    virtual int getNumShifts() = 0;
    virtual long int getNumBondsToForm() = 0;
    virtual size_t getNumAtoms() = 0;
    virtual size_t getNumBonds() = 0;
    virtual size_t getNumExtraAtoms() = 0;
    virtual size_t getNumExtraBonds() = 0;
    virtual size_t getNumParticles() = 0;
  };
}
}

#ifdef CEREALIZABLE
CEREAL_CLASS_VERSION(pylimer_tools::sim::OutputSupportingSimulation, 2);
#endif
#endif
