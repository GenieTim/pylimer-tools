#ifndef MEHP_FORCE_RELAX2_H
#define MEHP_FORCE_RELAX2_H

#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlopt.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {

#define NUM_COMPUTABLE_INT_VALUES 3

  enum ComputedIntValues
  {
    STEP = 0,
    NUM_SHIFT = 1,
    NUM_RELOC = 2,
  };

  const std::array<std::string, NUM_COMPUTABLE_INT_VALUES>
    ComputedIntValuesNames = {
      "Step",
      "numShift",
      "numReloc",
    };

#define NUM_COMPUTABLE_DOUBLE_VALUES 17

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
                                  "<b>",
                                  "max(b)",
                                  "MSD" };

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

  class OutputSupportingSimulation
  {
  protected:
    ////////////////////////////////////////////////////////////////
    // output configurations
    std::vector<OutputConfiguration> outputConfigs;
    std::vector<OutputConfiguration> outputAverageConfigs;
    std::vector<OutputConfiguration> outputAutoCorrelationConfigs;
    ////////////////////////////////////////////////////////////////
    // output streams
    std::vector<std::shared_ptr<std::ostream>> outputStreams;
    std::vector<int> outputFileStreams;
    ////////////////////////////////////////////////////////////////
    // computation state
    std::vector<pylimer_tools::calc::Correlator> autocorrelators;
    std::vector<Eigen::ArrayXi> msdMeasuredIndices;
    std::vector<Eigen::VectorXd> msdOrigins;
    std::vector<size_t> msdOriginTimesteps;

    int openFilesOutputHeader(const std::vector<OutputConfiguration>& configs,
                              const std::string& prefix = "",
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

  public:
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

    void configAverageOutput(const std::vector<OutputConfiguration>& vals)
    {
      this->outputAverageConfigs = vals;
    }

    void configStepOutput(const std::vector<OutputConfiguration>& vals)
    {
      this->outputConfigs = vals;
    }
  }
}
}

#endif MEHP_FORCE_RELAX2_H