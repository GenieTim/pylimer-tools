#include "OutputSupportingSimulation.h"


namespace pylimer_tools {
namespace sim {

  /**
   * @brief Open the specified files, and write the headers
   *
   * @param configs
   * @return int
   */
  int OutputSupportingSimulation::openFilesOutputHeader(
    const std::vector<OutputConfiguration>& configs,
    const std::string& prefix,
    int streamIdx)
  {
    INVALIDARG_EXP_IFN(streamIdx == this->outputStreams.size(),
                       "The stream idx " + std::to_string(streamIdx) +
                         " hints at an invalid state.");
    int numComputes = 0;
    this->outputStreams.reserve(streamIdx + configs.size());
    std::string outputBuffer = "";
    outputBuffer.reserve(80 * 20);
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

      outputBuffer = prefix;

      for (ComputedIntValues val : oc.intValues) {
        switch (val) {
          default:
            numComputes += 1;
            outputBuffer += ComputedIntValuesNames[val] + "\t";
        }
      }
      for (ComputedDoubleValues val : oc.doubleValues) {
        switch (val) {
          case ComputedDoubleValues::MSD:
            for (size_t i = 0; i < this->msdOrigins.size(); ++i) {
              outputBuffer += "MSD" + std::to_string(i) + "_" +
                              std::to_string(this->msdOriginTimesteps[i]) +
                              "\t";
            }
            numComputes += this->msdOrigins.size();
            break;
          default:
            numComputes += 1;
            outputBuffer += ComputedDoubleValuesNames[val] + "\t";
        }
      }

      if (!outputBuffer.empty()) {
        outputBuffer.pop_back(); // remove trailing tab
      }

      (*this->outputStreams[streamIdx]) << outputBuffer << std::endl;
      streamIdx += 1;
      outputBuffer.clear();
    }
    return numComputes;
  };
}
}
