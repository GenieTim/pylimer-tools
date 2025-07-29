#include "OutputSupportingSimulation.h"

namespace pylimer_tools::sim {

void
OutputSupportingSimulation::prepareAllOutputs()
{
  this->outputStreams.clear();
  this->outputFileStreams.clear();

  // output headers
  std::cout.flush();
  std::ios::sync_with_stdio(false);
  this->openFilesOutputHeader(this->outputConfigs);

  // prepare averages
  const int numAverages = this->openFilesOutputHeader(
    this->outputAverageConfigs, "# OutputStep\t", this->outputConfigs.size());
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
void
OutputSupportingSimulation::closeAllOutputs()
{
  // finish up
  std::ios::sync_with_stdio(true);
  std::cout.flush();
  this->outputStreams.clear();
  this->outputFileStreams.clear();
}
/**
 * @brief Open the specified files, and write the headers
 *
 * @param configs
 * @return int
 */
int
OutputSupportingSimulation::openFilesOutputHeader(
  const std::vector<OutputConfiguration>& configs,
  const std::string& prefix,
  size_t streamIdx)
{
  INVALIDARG_EXP_IFN(streamIdx == this->outputStreams.size(),
                     "The stream idx " + std::to_string(streamIdx) +
                       " hints at an invalid state.");
  size_t numComputes = 0;
  this->outputStreams.reserve(streamIdx + configs.size());
  std::string thisFileOutputBuffer = "";
  thisFileOutputBuffer.reserve(80 * 20);
  for (OutputConfiguration oc : configs) {
    if (oc.filename != "" && oc.filename != "stdio") {
      // always append, as the truncation happened already
      // (`this->validateAndTruncateOutputFiles`)
      this->outputStreams.push_back(std::make_shared<std::ofstream>(
        oc.filename, std::ios::out | std::ios::app));
      this->outputFileStreams.push_back(streamIdx);
    } else {
      this->outputStreams.push_back(
        std::shared_ptr<std::ostream>(&std::cout, [](void*) {}));
    }

    thisFileOutputBuffer = prefix;

    for (const ComputedIntValues val : oc.intValues) {
      switch (val) {
        default:
          numComputes += 1;
          thisFileOutputBuffer += ComputedIntValuesNames[val] + "\t";
      }
    }
    for (const ComputedDoubleValues val : oc.doubleValues) {
      switch (val) {
        case ComputedDoubleValues::MSD:
          for (size_t i = 0; i < this->msdOrigins.size(); ++i) {
            thisFileOutputBuffer +=
              "MSD" + std::to_string(i) + "_" +
              std::to_string(this->msdOriginTimesteps[i]) + "\t";
          }
          numComputes += this->msdOrigins.size();
          break;
        default:
          numComputes += 1;
          thisFileOutputBuffer += ComputedDoubleValuesNames[val] + "\t";
      }
    }

    if (!thisFileOutputBuffer.empty()) {
      thisFileOutputBuffer.pop_back(); // remove trailing tab
    }

    (*this->outputStreams[streamIdx]) << thisFileOutputBuffer << std::endl;
    streamIdx += 1;
    thisFileOutputBuffer.clear();
  }
  return numComputes;
}
void
OutputSupportingSimulation::handleOutput(const long int currentStep)
{
  // "lazily" compute the values we need, others less lazily when they are
  // computationally inexpensive
  const std::array<long int, NUM_COMPUTABLE_INT_VALUES> intvalues = {
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
  const double pressure = stressTensor.trace() / 3.;
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
  const std::array<double, NUM_COMPUTABLE_DOUBLE_VALUES> doublevalues = {
    { this->getTimestep(),
      this->requiresDEvaluation(ComputedDoubleValues::TIME, currentStep)
        ? this->getCurrentTime(currentStep)
        : 0.,
      this->requiresDEvaluation(ComputedDoubleValues::VOLUME, currentStep)
        ? this->getVolume()
        : 0.,
      pressure, // + kineticPressureTerm,
      this->requiresDEvaluation(ComputedDoubleValues::TEMPERATURE, currentStep)
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
      this->requiresDEvaluation(ComputedDoubleValues::GAMMA, currentStep)
        ? this->getGamma()
        : 0.,
      this->requiresDEvaluation(ComputedDoubleValues::RESIDUAL, currentStep)
        ? this->getResidual()
        : 0.,
      this->requiresDEvaluation(ComputedDoubleValues::MEAN_B, currentStep)
        ? bondLengths.mean()
        : 0.0,
      this->requiresDEvaluation(ComputedDoubleValues::MAX_B, currentStep)
        ? bondLengths.maxCoeff()
        : 0.0,
      0. }
  };
  size_t streamIdx = 0;
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
      const size_t previousAverageIdx = averagesIdx;
      if ((currentStep % oc.useEvery) == 0) {
        const double multiplier = (static_cast<double>(oc.useEvery) /
                                   static_cast<double>(oc.outputEvery));
        for (const ComputedIntValues val : oc.intValues) {
          switch (val) {
            default:
              runningAverages[averagesIdx] +=
                static_cast<double>(intvalues[val]) * multiplier;
              averagesIdx += 1;
              break;
          }
        }
        for (const ComputedDoubleValues val : oc.doubleValues) {
          switch (val) {
            case ComputedDoubleValues::MSD:
              // compute MSD
              for (msdIdx = 0; msdIdx < this->msdMeasuredIndices.size();
                   ++msdIdx) {
                const double result =
                  (this->msdOrigins[msdIdx] -
                   getCoordinates()(this->msdMeasuredIndices[msdIdx]))
                    .squaredNorm() /
                  (static_cast<double>(this->msdMeasuredIndices[msdIdx].size() /
                                       3.));
                runningAverages[averagesIdx + msdIdx] += result * multiplier;
              }
              averagesIdx += msdIdx;
              break;
            default:
              runningAverages[averagesIdx] += doublevalues[val] * multiplier;
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
    for (const ComputedDoubleValues cv : oc.doubleValues) {
      assert(autocorrelator_idx < this->autocorrelators.size());
      RUNTIME_EXP_IFN(std::isfinite(doublevalues[cv]),
                      "Expect output quantities to be finite, found " +
                        std::to_string(doublevalues[cv]) + " for property " +
                        ComputedDoubleValuesNames[cv] + ".");
      this->autocorrelators[autocorrelator_idx].add(doublevalues[cv]);
      autocorrelator_idx += 1;
    }
    if (currentStep % oc.outputEvery == 0) {
      outputBuffer += "# TimeStep " +
                      std::to_string(intvalues[ComputedIntValues::STEP]) + "\n";
      this->autocorrelators[autocorrelator_idx_before].evaluate();
      const unsigned int npcorr =
        this->autocorrelators[autocorrelator_idx_before].npcorr;
      RUNTIME_EXP_IFN(npcorr > 0, "Expected more than 0 correlator results.");
      for (size_t autocorr_idx_offset = 1;
           autocorr_idx_offset < oc.doubleValues.size();
           ++autocorr_idx_offset) {
        const size_t idx = autocorrelator_idx_before + autocorr_idx_offset;
        this->autocorrelators[idx].evaluate();
        RUNTIME_EXP_IFN(this->autocorrelators[idx].npcorr == npcorr,
                        "Autocorrelation states are inconsistent.");
      }

      for (size_t output_idx = 0; output_idx < static_cast<size_t>(npcorr);
           output_idx += 1) {
        outputBuffer += std::to_string(
          this->autocorrelators[autocorrelator_idx_before].t[output_idx]);
        for (size_t autocorr_idx_offset = 0;
             autocorr_idx_offset < oc.doubleValues.size();
             ++autocorr_idx_offset) {
          const size_t idx = autocorrelator_idx_before + autocorr_idx_offset;
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
void
OutputSupportingSimulation::doOutputValues(
  const OutputConfiguration& oc,
  const std::array<long int, 8>& intValues,
  const std::array<double, 19>& doubleValues,
  const size_t streamIdx)
{
  assert(streamIdx <= this->outputStreams.size());
  assert(doubleValues.size() == NUM_COMPUTABLE_DOUBLE_VALUES);
  for (const ComputedIntValues val : oc.intValues) {
    RUNTIME_EXP_IFN(std::isfinite(static_cast<double>(intValues[val])),
                    "Expect output quantities to be finite, found " +
                      std::to_string(intValues[val]) + " for property " +
                      ComputedIntValuesNames[val] + ".");
    switch (val) {
      default:
        outputBuffer += std::to_string(intValues[val]) + "\t";
    }
  }
  for (const ComputedDoubleValues val : oc.doubleValues) {
    RUNTIME_EXP_IFN(std::isfinite(doubleValues[val]),
                    "Expect output quantities to be finite, found " +
                      std::to_string(doubleValues[val]) + " for property " +
                      ComputedDoubleValuesNames[val] + ".");
    switch (val) {
      case ComputedDoubleValues::MSD:
        // compute MSD
        for (size_t msdIdx = 0; msdIdx < this->msdMeasuredIndices.size();
             ++msdIdx) {
          assert(this->msdOrigins.size() > msdIdx);
          assert(this->msdOrigins[msdIdx].size() ==
                 this->msdMeasuredIndices[msdIdx].size());
          Eigen::VectorXd relCoords = this->getCoordinates();
          const int maxCoeff = this->msdMeasuredIndices[msdIdx].maxCoeff();
          std::cout << "Size: " << relCoords.size() << ", max coeff "
                    << maxCoeff << std::endl;
          assert(relCoords.size() >
                 this->msdMeasuredIndices[msdIdx].maxCoeff());
          Eigen::ArrayXi nIndices = this->msdMeasuredIndices[msdIdx];
          assert(nIndices.size() > 0);
          assert(nIndices.size() == this->msdOrigins[msdIdx].size());
          assert(nIndices.minCoeff() >= 0 &&
                 nIndices.maxCoeff() < relCoords.size());
          const double result =
            (this->msdOrigins[msdIdx] - relCoords(nIndices)).squaredNorm() /
            (static_cast<double>(this->msdMeasuredIndices[msdIdx].size() / 3.));
          outputBuffer += std::to_string(result) + "\t";
        }
        break;
      default:
        outputBuffer += //"(" + std::to_string(val) + ")" +
          std::to_string(doubleValues[val]) + "\t";
    }
  }
  if (!outputBuffer.empty()) {
    outputBuffer.pop_back(); // remove last "\t"
    outputBuffer += "\n";
    // output the buffer, clear it
    (*(this->outputStreams[streamIdx])) << outputBuffer;
    outputBuffer.clear();
  }
}
void
OutputSupportingSimulation::updateValuesRequiredEvery(
  const std::vector<OutputConfiguration>& configs)
{
  for (OutputConfiguration c : configs) {
    for (const ComputedDoubleValues v : c.doubleValues) {
      if (this->doubleValueRequiredEvery[v] == 0) {
        this->doubleValueRequiredEvery[v] = c.useEvery;
      } else {
        this->doubleValueRequiredEvery[v] =
          std::gcd(c.useEvery, this->doubleValueRequiredEvery[v]);
      }
    }
    for (const ComputedIntValues i : c.intValues) {
      if (this->intValueRequiredEvery[i] == 0) {
        this->intValueRequiredEvery[i] = c.useEvery;
      } else {
        this->intValueRequiredEvery[i] =
          std::gcd(c.useEvery, this->intValueRequiredEvery[i]);
      }
    }
  }
  const std::vector<ComputedDoubleValues> stressTensorRequiringValues = {
    STRESS_XX, STRESS_YY,  STRESS_ZZ,  STRESS_XY,  STRESS_YZ,
    STRESS_XZ, STRESS_NXY, STRESS_NYZ, STRESS_NXZ, PRESSURE
  };
  for (const ComputedDoubleValues v : stressTensorRequiringValues) {
    if (this->requireStressTensorEvery == 0) {
      this->requireStressTensorEvery = this->doubleValueRequiredEvery[v];
    } else {
      this->requireStressTensorEvery = std::gcd(
        this->requireStressTensorEvery, this->doubleValueRequiredEvery[v]);
    }
  }
  // similarly for the bond length
  this->requireBondLenEvery =
    std::gcd(this->doubleValueRequiredEvery[ComputedDoubleValues::MAX_B],
             this->doubleValueRequiredEvery[ComputedDoubleValues::MEAN_B]);
}
void
OutputSupportingSimulation::validateAndTruncateOutputFiles(
  const std::vector<OutputConfiguration>& vals)
{
  for (size_t i = 0; i < vals.size(); ++i) {
    if (vals[i].filename.size() > 0) {
      // empty the file
      std::ifstream file;
      file.open(vals[i].filename.c_str(),
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
void
OutputSupportingSimulation::configAutoCorrelatorOutput(
  const std::vector<OutputConfiguration>& vals,
  const unsigned int numcorrin,
  const unsigned int pin,
  const unsigned int min)
{
  size_t num_values_to_correlate = 0;
  for (size_t i = 0; i < vals.size(); ++i) {
    INVALIDARG_EXP_IFN(vals[i].intValues.size() == 0,
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
void
OutputSupportingSimulation::configAverageOutput(
  const std::vector<OutputConfiguration>& configs)
{
  this->outputAverageConfigs = configs;
  this->updateValuesRequiredEvery(configs);

  size_t numAverages = 0;
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
void
OutputSupportingSimulation::configStepOutput(
  std::vector<OutputConfiguration>& vals)
{
  for (size_t i = 0; i < vals.size(); ++i) {
    vals[i].useEvery = vals[i].outputEvery;
  }

  this->validateAndTruncateOutputFiles(vals);
  this->outputConfigs = vals;
  this->updateValuesRequiredEvery(vals);
}
void
OutputSupportingSimulation::configRestartOutput(const std::string& outputFile,
                                                const int outputEvery)
{
  this->outputRestartEvery = outputEvery;
  this->restartOutputFile = outputFile;
};
}
