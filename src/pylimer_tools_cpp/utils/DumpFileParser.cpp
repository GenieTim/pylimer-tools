#include "DumpFileParser.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"
#include "StringUtil.h"
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

/**
 * @brief Initialize the parser to read from a certain file path
 *
 * @param filePath
 */
void DumpFileParser::startReading(const std::string filePath) {

  if (!std::filesystem::exists(filePath)) {
    throw std::invalid_argument("File to read (" + filePath +
                                ") does not exist.");
  }

  std::string line;
  this->file.open(filePath);

  if (!this->file.is_open()) {
    throw std::invalid_argument("File to read file (" + filePath +
                                "): failed to open.");
  }

  // read everything until the first key
  while (getline(this->file, line)) {
    line = pylimer_tools::utils::trimLineOmitComment(line);
    // skip empty lines: break when not empty
    if (!line.empty()) {
      break;
    }
  }

  // Assemble CSV data for all keys
  this->newGroupKey = line; // new group key: key for a new timestep (group)
  this->currentLine = line; // current line
  this->groupPosMap[0] = this->file.tellg(); // record position of index to jump back at some point
};

// TODO: implement routine to read a group at any position
// e.g. by tellg() (https://www.cplusplus.com/reference/istream/istream/tellg/)
// together with seekg (https://www.cplusplus.com/reference/istream/istream/seekg/)
void DumpFileParser::readGroupByIdx(const int i) {
  // if (this->)
}


/**
 * @brief Read N timesteps
 *
 * @param N the nr of groups to read; a negative value results in all groups
 * being read.
 */
void DumpFileParser::readNGroups(const int N) {
  if (this->newGroupKey.empty()) {
    throw std::runtime_error("Starting reading first, specify the file");
  }

  int groupsRead = 0;
  std::string currentKey = this->cleanHeader(this->currentLine);
  int currentNrOfExpectedGroups = this->headerColMap[currentKey].size();
  data_item_t dataItem;
  std::string line = this->currentLine;

  while (getline(this->file, line)) {
    line = pylimer_tools::utils::trimLineOmitComment(line);
    // skip empty lines
    if (line.empty()) {
      continue;
    }
    // new header
    if (pylimer_tools::utils::startsWith(line, "ITEM:")) {
      currentKey = this->cleanHeader(line);
      currentNrOfExpectedGroups = this->headerColMap[currentKey].size();
    } else {
      dataItem[currentKey].push_back(
          pylimer_tools::utils::CsvTokenizer(line, currentNrOfExpectedGroups));
    }

    if (line == newGroupKey) {
      // new timestep
      groupsRead += 1;
      if ((N > 0 && groupsRead >= N)) {
        break;
      }
      this->data.push_back(dataItem);
      this->groupPosMap[this->data.size()] = this->file.tellg();
      dataItem = data_item_t();
    }
  }

  this->currentLine = line;
  // last timestep
  this->data.push_back(dataItem);
};

/**
 * @brief Forget the data at a certain index
 *
 * @param index
 */
void DumpFileParser::forgetAt(const int index) {
  this->data.erase(this->data.begin() + index);
};

/**
 * @brief Read a whole file
 *
 * @param filePath
 */
void DumpFileParser::read(const std::string filePath) {
  this->startReading(filePath);
  this->readNGroups(-1);
  this->finishedReading = true;
  this->finish();
}

void DumpFileParser::finish() { this->file.close(); }

std::string DumpFileParser::cleanHeader(std::string headerToClean) {
  // "ITEM: ".size() = 6
  headerToClean.erase(0, 6);
  pylimer_tools::utils::CsvTokenizer tokenizer(headerToClean);

  std::string newHeader = "";
  std::vector<std::string> columns;
  for (size_t i = 0; i < tokenizer.getLength(); ++i) {
    std::string beg = tokenizer.get<std::string>(i);
    if (isUpper(beg)) {
      newHeader.append(beg);
      newHeader.append(" ");
    } else {
      columns.push_back(beg);
    }
  }
  newHeader = pylimer_tools::utils::rtrim(newHeader);
  if (!this->headerColMap.contains(newHeader)) {
    this->headerColMap.insert_or_assign(newHeader, columns);
  }

  return newHeader;
}
} // namespace utils
} // namespace pylimer_tools
