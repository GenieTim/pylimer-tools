#ifndef OUTPUT_SUPPORTING_SIM_H
#define OUTPUT_SUPPORTING_SIM_H

#include "../calc/Correlator.h"
#include "../utils/CerealUtils.h"
#include "../utils/VectorUtils.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#ifdef CEREALIZABLE
#include <cereal/types/polymorphic.hpp>
#endif

namespace pylimer_tools::sim {

#define COMPUTED_INT_VALUES                                                    \
  X(STEP, "Step")                                                              \
  X(NUM_SHIFT, "numShift")                                                     \
  X(NUM_RELOC, "numReloc")                                                     \
  X(NUM_ATOMS, "numAtoms")                                                     \
  X(NUM_EXTRA_ATOMS, "numExtraAtoms")                                          \
  X(NUM_BONDS, "numBonds")                                                     \
  X(NUM_EXTRA_BONDS, "numExtraBonds")                                          \
  X(NUM_BONDS_TO_FORM, "numBondsToForm")

#define NUM_COMPUTABLE_INT_VALUES 8

enum ComputedIntValues
{
#define X(name, str) name,
  COMPUTED_INT_VALUES
#undef X
};

const std::array<std::string, NUM_COMPUTABLE_INT_VALUES>
  ComputedIntValuesNames = {
#define X(name, str) str,
    COMPUTED_INT_VALUES
#undef X
  };

#define COMPUTED_DOUBLE_VALUES                                                 \
  X(TIMESTEP, "TimeStep")                                                      \
  X(TIME, "Time")                                                              \
  X(VOLUME, "Volume")                                                          \
  X(PRESSURE, "Pressure")                                                      \
  X(TEMPERATURE, "Temperature")                                                \
  X(STRESS_XX, "Stress[0,0]")                                                  \
  X(STRESS_YY, "Stress[1,1]")                                                  \
  X(STRESS_ZZ, "Stress[2,2]")                                                  \
  X(STRESS_XY, "Stress[0,1]")                                                  \
  X(STRESS_YZ, "Stress[1,2]")                                                  \
  X(STRESS_XZ, "Stress[0,2]")                                                  \
  X(STRESS_NXY, "Stress[0,0]-Stress[1,1]")                                     \
  X(STRESS_NYZ, "Stress[1,1]-Stress[2,2]")                                     \
  X(STRESS_NXZ, "Stress[0,0]-Stress[2,2]")                                     \
  X(GAMMA, "Gamma")                                                            \
  X(RESIDUAL, "Residual")                                                      \
  X(MEAN_B, "<b>")                                                             \
  X(MAX_B, "max(b)")                                                           \
  X(MSD, "MSD")

#define NUM_COMPUTABLE_DOUBLE_VALUES 19

enum ComputedDoubleValues
{
#define X(name, str) name,
  COMPUTED_DOUBLE_VALUES
#undef X
};

const std::array<std::string, NUM_COMPUTABLE_DOUBLE_VALUES>
  ComputedDoubleValuesNames = {
#define X(name, str) str,
    COMPUTED_DOUBLE_VALUES
#undef X
  };

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
  void prepareAllOutputs();

  /**
   * @brief Close all files required for output
   *
   */
  void closeAllOutputs();

  int openFilesOutputHeader(const std::vector<OutputConfiguration>& configs,
                            const std::string& prefix = "",
                            size_t streamIdx = 0);

  inline bool requiresDEvaluation(const ComputedDoubleValues val,
                                  const long int currentStep) const
  {
    const bool requiresEval =
      (this->doubleValueRequiredEvery[val] > 0) &&
      ((currentStep % this->doubleValueRequiredEvery[val]) == 0);
    // if (requiresEval) {
    //   std::cout << "Requiring evaluating " <<
    //   ComputedDoubleValuesNames[val]
    //             << " at step " << currentStep << std::endl;
    // }
    return requiresEval;
  }

  inline bool requiresIEvaluation(const ComputedIntValues val,
                                  const long int currentStep) const
  {
    return (this->intValueRequiredEvery[val] > 0) &&
           ((currentStep % this->intValueRequiredEvery[val]) == 0);
  }

  /**
   * @brief Iterate all possible output configurations, handle them
   *
   * @param currentStep
   */
  void handleOutput(const long int currentStep);

  // static OutputSupportingSimulation readRestartFile(std::string filename)
  // {
  //   throw std::runtime_error(
  //     "Cannot read restart file on abstract base class.");
  // };

  /**
   * @brief Output the passed values
   *
   * @param oc
   * @param intValues
   * @param doubleValues
   * @param streamIdx
   */
  void doOutputValues(
    const OutputConfiguration& oc,
    const std::array<long int, NUM_COMPUTABLE_INT_VALUES>& intValues,
    const std::array<double, NUM_COMPUTABLE_DOUBLE_VALUES>& doubleValues,
    const size_t streamIdx = 0);

  /**
   * @brief Remember how often a particular value is needed to be computed for
   * any of the outputs
   *
   * @param configs
   */
  void updateValuesRequiredEvery(
    const std::vector<OutputConfiguration>& configs);

public:
  virtual ~OutputSupportingSimulation() = default;

  OutputSupportingSimulation()
  {
    this->doubleValueRequiredEvery.fill(0);
    this->intValueRequiredEvery.fill(0);
    this->outputBuffer.reserve(600);
  }

#ifdef CEREALIZABLE
  virtual void writeRestartFile(std::string& filename) = 0;
#endif

  /**
   *
   * @param vals the output configurations whose output to validate and truncate
   */
  static void validateAndTruncateOutputFiles(
    const std::vector<OutputConfiguration>& vals);

  /**
   * Configure the simulation to do autocorrelation computation
   *
   * @param vals the autocorrelation output configurations
   * @param numcorrin a parameter for the autocorrelation computation
   * @param pin a parameter for the autocorrelation computation
   * @param min a parameter for the autocorrelation computation
   */
  void configAutoCorrelatorOutput(const std::vector<OutputConfiguration>& vals,
                                  const unsigned int numcorrin = 32,
                                  const unsigned int pin = 16,
                                  const unsigned int min = 2);

  /**
   *
   * @param configs the configuration for the output and computation of the
   * averages
   */
  void configAverageOutput(const std::vector<OutputConfiguration>& configs);

  /**
   *
   * @param vals the configuration for the output and computation per step
   */
  void configStepOutput(std::vector<OutputConfiguration>& vals);

  /**
   *
   * @param outputFile the file-path to which the restart file should be written
   * @param outputEvery how often to write the restart file
   */
  void configRestartOutput(const std::string& outputFile,
                           const int outputEvery);

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
  virtual double getTemperature() = 0;
  virtual double getGamma() = 0;
  virtual double getTimestep() = 0;
  virtual double getVolume() = 0;
  virtual double getResidual() = 0;
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

#ifdef CEREALIZABLE
CEREAL_CLASS_VERSION(pylimer_tools::sim::OutputSupportingSimulation, 2);
#endif
#endif
