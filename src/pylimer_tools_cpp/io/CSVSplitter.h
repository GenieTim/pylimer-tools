
#include "../utils/StringUtils.h"
#include <filesystem>
#include <fstream> // std::ifstream
#include <iostream>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {
  std::vector<std::string> splitCSV(const std::string filePath)
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

    while (std::getline(file, line)) {
      shortenedLine = pylimer_tools::utils::trimLineOmitComment(line);
      pylimer_tools::utils::CsvTokenizer tokenizer(line);
      if (tokenizer.getLength() != previousLen) {
        // todo: also check for headers with same length?!?
        previousLen = tokenizer.getLength();
        outputFile.close();
        currentTempFile = std::tmpnam(nullptr);
        // let's not care about the first one
        results.push_back(currentTempFile);
        outputFile.open(currentTempFile);
      }
      outputFile << line;
    }

    outputFile.close();

    return results;
  }
}
}
