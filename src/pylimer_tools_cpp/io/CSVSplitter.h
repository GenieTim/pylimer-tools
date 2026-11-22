
#pragma once

#include "../utils/StringUtils.h"
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream> // std::ifstream
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {
  inline std::string makeTempCsvFilePath()
  {
    static std::atomic<size_t> counter{ 0 };
    const std::filesystem::path tempDir =
      std::filesystem::temp_directory_path();
    const long long timestamp =
      std::chrono::steady_clock::now().time_since_epoch().count();
    while (true) {
      const size_t id = counter.fetch_add(1, std::memory_order_relaxed);
      const std::filesystem::path candidate =
        tempDir / ("pylimer_csv_" + std::to_string(timestamp) + "_" +
                   std::to_string(id) + ".tmp");
      if (!std::filesystem::exists(candidate)) {
        return candidate.string();
      }
    }
  }

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
    std::string currentTempFile;
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
      if (newLength != previousLen ||
          (!tokenzierResults.empty() && tokenzierResults[0] == "Step")) {
        // todo: also check for headers with same length?!?
        previousLen = newLength;
        if (outputFile.is_open()) {
          outputFile.close();
        }
        currentTempFile = makeTempCsvFilePath();
        results.push_back(currentTempFile);
        outputFile.open(currentTempFile);
      }
      outputFile << shortenedLine << "\n";
    }

    if (outputFile.is_open()) {
      outputFile.close();
    }

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
