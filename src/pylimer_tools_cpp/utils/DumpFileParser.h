#ifndef DUMP_FILE_PARSER_H
#define DUMP_FILE_PARSER_H

#include "../entities/Universe.h"
#include "StringUtils.h"
#include "TYPES.h"
#include <algorithm>
#include <any>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream> // std::ifstream
#include <map>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {
// types
typedef std::map<std::string, std::vector<pylimer_tools::utils::CsvTokenizer>>
    data_item_t;

class DumpFileParser {
public:
  DumpFileParser(){};
  DumpFileParser(const std::string filePath);

  // rule of three:
  // 1. destructor (to destroy the graph)
  ~DumpFileParser();
  // 2. copy constructor
  DumpFileParser(const DumpFileParser &src);
  // 3. copy assignment operator
  DumpFileParser &operator=(DumpFileParser src);

  void read();
  void finish();
  void readNGroups(const INDEX_TYPE start, const int N);
  void readGroupByIdx(const INDEX_TYPE i);
  void forgetAt(const INDEX_TYPE index);

  template <typename OUT>
  std::vector<OUT> getValuesForAt(const INDEX_TYPE index,
                                  const std::string headerKey,
                                  const std::string &column);
  template <typename OUT>
  std::vector<OUT> getValuesForAt(const INDEX_TYPE index,
                                  const std::string headerKey,
                                  const INDEX_TYPE column);
  // the next two methods are specializations for easier py binding
  std::vector<std::string> getStringValuesForAt(const INDEX_TYPE index,
                                                const std::string headerKey,
                                                const std::string column) {
    return this->getValuesForAt<std::string>(index, headerKey, column);
  };
  std::vector<double> getNumericValuesForAt(const INDEX_TYPE index,
                                            const std::string headerKey,
                                            const std::string column) {
    return this->getValuesForAt<double>(index, headerKey, column);
  };
  int getLength() { return this->nrOfGroups; }
  bool hasKey(std::string headerKey) {
    if (this->data.size() == 0) {
      throw std::invalid_argument("Cannot check for header '" + headerKey + "' without reading a group first.");
    }
    if (this->getLength() > 0) {
      // assumption: each step has the same keys
      auto firstEl = this->data.begin();
      return firstEl->second.contains(headerKey);
    }
    return false;
  }
  bool keyHasColumn(std::string headerKey, std::string column) {
    const auto colItIdx =
        std::find(this->headerColMap.at(headerKey).begin(),
                  this->headerColMap.at(headerKey).end(), column);
    if (this->headerColMap.at(headerKey).end() == colItIdx) {
      return false;
    }
    return true;
  }
  bool keyHasDirectionalColumn(std::string headerKey, std::string dirPraefix,
                               std::string dirSuffix) {
    // std::cout << "Searching for " << headerKey << " " << dirPraefix <<
    // dirSuffix << " in " <<
    // pylimer_tools::utils::join(this->headerColMap.at(headerKey).begin(),
    // this->headerColMap.at(headerKey).end(), std::string(" ")) << std::endl;
    return this->keyHasColumn(headerKey, dirPraefix + "x" + dirSuffix) &&
           this->keyHasColumn(headerKey, dirPraefix + "y" + dirSuffix) &&
           this->keyHasColumn(headerKey, dirPraefix + "z" + dirSuffix);
  }

private:
  std::string cleanHeader(std::string header);

  template <typename OUT>
  inline std::vector<OUT> parseTypesInLine(std::string line) {
    std::vector<OUT> resultnumbers;
    pylimer_tools::utils::CsvTokenizer tokenizer(line);
    resultnumbers.reserve(tokenizer.getLength());
    for (size_t i = 0; i < tokenizer.getLength(); ++i) {
      resultnumbers.push_back(tokenizer.get<OUT>(i));
    }
    return resultnumbers;
  }

  //// data
  std::string filePath;
  std::string newGroupKey;
  std::string currentLine;
  std::ifstream file;
  int nrOfGroups;
  std::unordered_map<INDEX_TYPE, data_item_t> data;
  std::map<std::string, std::vector<std::string>> headerColMap;
  std::map<INDEX_TYPE, int> groupPosMap;
};

template <typename OUT>
std::vector<OUT> DumpFileParser::getValuesForAt(const INDEX_TYPE index,
                                                const std::string headerKey,
                                                const std::string &column) {
  // detect index of column
  INDEX_TYPE colIdx = 0;
  if (this->headerColMap.at(headerKey).size() > 1) {
    const auto colItIdx =
        std::find(this->headerColMap.at(headerKey).begin(),
                  this->headerColMap.at(headerKey).end(), column);
    if (this->headerColMap.at(headerKey).end() == colItIdx) {
      throw std::invalid_argument("Column '" + column +
                                  "' not found for header '" + headerKey + "'");
    }
    colIdx = colItIdx - this->headerColMap.at(headerKey).begin();
  }
  return this->getValuesForAt<OUT>(index, headerKey, colIdx);
}

template <typename OUT>
std::vector<OUT> DumpFileParser::getValuesForAt(const INDEX_TYPE index,
                                                const std::string headerKey,
                                                const INDEX_TYPE colIdx) {
  if (!this->data.contains(index)) {
    // std::cout << "Could not find index " << index << "in data yet" << std::endl;
    this->readGroupByIdx(index);
  }
  // std::cout << "Requested values for index " << index << ", key " << headerKey << " and column " << colIdx << std::endl; 
  data_item_t dataItem = this->data.at(index);
  //
  std::vector<pylimer_tools::utils::CsvTokenizer> relevantData =
      dataItem[headerKey];
  std::vector<OUT> results;

  for (pylimer_tools::utils::CsvTokenizer lineTok : relevantData) {
    results.push_back(lineTok.get<OUT>(colIdx));
  }

  return results;
}
} // namespace utils
} // namespace pylimer_tools

#endif
