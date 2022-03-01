#include "DumpFileParser.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"
#include "StringUtils.h"
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

// rule of three:
// 1. destructor (to destroy the graph)
DumpFileParser::~DumpFileParser() { this->finish(); };
// 2. copy constructor
DumpFileParser::DumpFileParser(const DumpFileParser &src) {
  this->currentLine = src.currentLine;
  this->newGroupKey = src.newGroupKey;
  this->nrOfGroups = src.nrOfGroups;
  this->data = src.data;
  this->headerColMap = src.headerColMap;
  this->groupPosMap = src.groupPosMap;
  this->filePath = src.filePath;
  // cannot clone ifstream
  this->file.clear();
  this->file.open(src.filePath);
};
// 3. copy assignment operator
DumpFileParser &DumpFileParser::operator=(DumpFileParser src) {
  std::swap(this->currentLine, src.currentLine);
  std::swap(this->newGroupKey, src.newGroupKey);
  std::swap(this->nrOfGroups, src.nrOfGroups);
  std::swap(this->data, src.data);
  std::swap(this->headerColMap, src.headerColMap);
  std::swap(this->groupPosMap, src.groupPosMap);
  // cannot clone ifstream
  if (src.file.is_open()) {
    src.file.close();
  }
  if (this->file.is_open()) {
    this->file.close();
  }
  this->file.clear();
  this->file.open(src.filePath);
  src.file.open(this->filePath);
  std::swap(this->filePath, src.filePath);
  return *this;
};

/**
 * @brief Initialize the parser to read from a certain file path
 *
 * @param filePath
 */
DumpFileParser::DumpFileParser(const std::string filePath) {
  this->filePath = filePath;
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
  this->groupPosMap.emplace(
      0, this->file
             .tellg()); // record position of index to jump back at some point

  // read the whole file, skipping all lines that are not the new group key
  // to record the positions
  int groupsFound = 0;

  while (getline(this->file, line)) {
    line = pylimer_tools::utils::trimLineOmitComment(line);
    // skip empty lines
    if (line.empty()) {
      continue;
    }

    if (line == this->newGroupKey) {
      // new timestep
      groupsFound += 1;
      this->groupPosMap.emplace(groupsFound, this->file.tellg());
    }
  }

  this->nrOfGroups = groupsFound + 1;
  // reset position to start of first group
  this->file.clear();
  this->file.seekg(this->groupPosMap.at(0));
}

/**
 * @brief Read a group by its index
 *
 * Useful for not-having to read all universes at once if only interested in
 * one. Position of groups is determined using tellg()
 * (https://www.cplusplus.com/reference/istream/istream/tellg/), whereas the
 * returned position is found again using seekg()
 * (https://www.cplusplus.com/reference/istream/istream/seekg/)
 *
 * @param i the index of the group to read
 */
void DumpFileParser::readGroupByIdx(const size_t i) { this->readNGroups(i, 1); }

/**
 * @brief Read N timesteps
 *
 * @param start the index to start at reading
 * @param N the nr of groups to read; a negative value results in all groups
 * being read.
 */
void DumpFileParser::readNGroups(const size_t start, const int N) {
  if (!this->file.is_open()) {
    throw std::runtime_error("Cannot read from closed file.");
  }

  if (start >= this->getLength() || ((int)start + N) > this->getLength()) {
    throw std::invalid_argument("Cannot read from outside the length of the "
                                "dump file. Tried to read from " +
                                std::to_string(start) + " to " +
                                std::to_string(N) + " for a file with " +
                                std::to_string(this->getLength()) +
                                " time-steps.");
  }

  if (this->file.eof()) {
    this->file.clear();
  }
  this->file.seekg(this->groupPosMap.at(start));

  size_t groupsRead = 0;
  std::string currentKey = this->cleanHeader(this->newGroupKey);
  int currentNrOfExpectedGroups = this->headerColMap.at(currentKey).size();
  data_item_t dataItem;
  std::string line = this->currentLine;
  // std::cout << "Starting to read at " << start << " for " << N << " with "
  //           << line << " and key " << currentKey << std::endl;
  int linesRead = 0;

  while (std::getline(this->file, line)) {
    linesRead += 1;
    // std::cout << "Read line: " << line << std::endl;
    line = pylimer_tools::utils::trimLineOmitComment(line);
    if (line.empty()) {
      continue;
    }
    // new header
    if (pylimer_tools::utils::startsWith(line, "ITEM:")) {
      currentKey = this->cleanHeader(line);
      currentNrOfExpectedGroups = this->headerColMap.at(currentKey).size();
    } else {
      // std::cout << "Appending data: " << currentNrOfExpectedGroups << " to "
      // << currentKey << std::endl;
      dataItem[currentKey].push_back(
          pylimer_tools::utils::CsvTokenizer(line, currentNrOfExpectedGroups));
    }

    if (line == this->newGroupKey) {
      // new timestep
      if ((N > 0 && (groupsRead + 1) >= N)) {
        break;
      }
      this->data.insert_or_assign(start + groupsRead, dataItem);
      groupsRead += 1;
      dataItem = data_item_t();
    }
  }

  if (linesRead > 1) {
    // last timestep
    this->data.insert_or_assign(start + groupsRead, dataItem);
    groupsRead += 1;
  }
  if (groupsRead != N) {
    throw std::runtime_error("Failed to read " + std::to_string(N) + ", read " +
                             std::to_string(groupsRead) + " groups on " +
                             std::to_string(linesRead) +
                             " lines. Stream error: " + std::strerror(errno));
    this->file.clear();
  }
  // std::cout << "Read " << groupsRead << " groups "
  //           << "(last " << start + groupsRead - 1 << ") of " <<
  //           this->getLength()
  //           << std::endl;
};

/**
 * @brief Forget the data at a certain index
 *
 * @param index
 */
void DumpFileParser::forgetAt(const size_t index) {
  if (this->data.contains(index)) {
    this->data.erase(index);
  }
};

/**
 * @brief Read a whole file
 *
 * @param filePath
 */
void DumpFileParser::read() {
  this->data.reserve(this->getLength());
  this->readNGroups(0, -1);
  this->finish();
}

void DumpFileParser::finish() {
  if (this->file.is_open()) {
    this->file.close();
  }
}

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
