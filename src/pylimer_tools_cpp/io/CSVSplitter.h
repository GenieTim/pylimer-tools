
#include "../utils/StringUtils.h"
#include <filesystem>
#include <fstream> // std::ifstream
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {
  inline std::vector<std::string> splitCSV(const std::string& filePath,
                                           const std::string& delimiter)
  {
    if (!std::filesystem::exists(filePath)) {
      throw std::invalid_argument("File to read (" + filePath +
                                  ") does not exist.");
    }

    std::vector<std::string> results;

    std::string line;
    std::ifstream file;
    file.open(filePath);

    if (!file.is_open()) {
      throw std::invalid_argument("File to read (" + filePath +
                                  "): failed to open.");
    }

    std::ofstream outputFile;
    std::string currentTempFile = std::tmpnam(nullptr);
    outputFile.open(currentTempFile);
    int previousLen = -1;
    std::string shortenedLine = "";
    std::vector<std::string> tokenzierResults;
    tokenzierResults.reserve(16);

    while (std::getline(file, line)) {
      // omit comment
      shortenedLine = pylimer_tools::utils::rstrip(line, "#");
      shortenedLine = pylimer_tools::utils::trim(shortenedLine);
      int newLength =
        pylimer_tools::utils::split(tokenzierResults, shortenedLine, delimiter);
      if (newLength != previousLen || tokenzierResults[0] == "Step") {
        // todo: also check for headers with same length?!?
        previousLen = newLength;
        outputFile.close();
        currentTempFile = std::tmpnam(nullptr);
        // let's not care about the first one
        results.push_back(currentTempFile);
        outputFile.open(currentTempFile);
      }
      outputFile << shortenedLine << "\n";
    }

    outputFile.close();

    return results;
  }

  inline std::string mergeCSVFiles(const std::vector<std::string>& files)
  {
    std::set<std::string> columns;
    for (std::string file : files) {
      std::ifstream inFile;
      inFile.open(file);
      std::string line;
      std::getline(inFile, line);
      pylimer_tools::utils::CsvTokenizer tokenizer(line);
    }
    return "";
  }
}
}
