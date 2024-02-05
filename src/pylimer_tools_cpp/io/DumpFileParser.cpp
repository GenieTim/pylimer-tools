#include "DumpFileParser.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
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
  DumpFileParser::~DumpFileParser()
  {
    this->finish();
  };
  // 2. copy constructor
  DumpFileParser::DumpFileParser(const DumpFileParser& src)
  {
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
  DumpFileParser& DumpFileParser::operator=(DumpFileParser src)
  {
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
  DumpFileParser::DumpFileParser(const std::string filePath)
  {
    if (!std::filesystem::exists(filePath)) {
      throw std::invalid_argument("File to read (" + filePath +
                                  ") does not exist.");
    }
    this->filePath = filePath;

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
      0,
      this->file
        .tellg()); // record position of index to jump back at some point

    // read the whole file, skipping all lines that are not the new group key
    // to record the positions
    int groupsFound = 0;
    size_t linesSinceLastIgnore = 0;
    while (getline(this->file, line)) {
      linesSinceLastIgnore += 1;
      // performance improvement: not skipping.
      // could be bad for certain files.
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

      // skip forward. Could be too far, in principle.
      if (linesSinceLastIgnore > 2) {
        // jump to next occurrence of the first character of the desired line
        this->file.ignore(std::numeric_limits<std::streamsize>::max(),
                          this->file.widen(this->newGroupKey[0]));
        this->file.unget(); // put the character back
        linesSinceLastIgnore = 0;
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
  void DumpFileParser::readGroupByIdx(const size_t i)
  {
    this->readNGroups(i, 1);
  }

  /**
   * @brief Get the nr of groups
   *
   * @return size_t
   */
  size_t DumpFileParser::getLength()
  {
    return this->nrOfGroups;
  }

  /**
   * @brief Check whether a header key exists
   *
   * @param headerKey
   * @return bool
   */
  bool DumpFileParser::hasKey(std::string headerKey)
  {
    if (this->data.size() == 0) {
      throw std::invalid_argument("Cannot check for header '" + headerKey +
                                  "' without reading a group first.");
    }
    return pylimer_tools::utils::map_has_key(this->headerColMap, headerKey);
  }

  /**
   * @brief Check whether a header key has a certain column
   *
   * @param headerKey
   * @param column
   * @return bool
   */
  bool DumpFileParser::keyHasColumn(std::string headerKey, std::string column)
  {
    const auto colItIdx = std::find(this->headerColMap.at(headerKey).begin(),
                                    this->headerColMap.at(headerKey).end(),
                                    column);
    if (this->headerColMap.at(headerKey).end() == colItIdx) {
      return false;
    }
    return true;
  }

  /**
   * @brief Check whether a header key has a certain column three times
   *
   * @param headerKey the key to check
   * @param dirPraefix the praefix in front of the "x", "y" and "z" of the
   * column
   * @param dirSuffix the suffix behind the "x", "y" and "z" of the column
   * @return bool
   */
  bool DumpFileParser::keyHasDirectionalColumn(std::string headerKey,
                                               std::string dirPraefix,
                                               std::string dirSuffix)
  {
    // std::cout << "Searching for " << headerKey << " " << dirPraefix <<
    // dirSuffix << " in " <<
    // pylimer_tools::utils::join(this->headerColMap.at(headerKey).begin(),
    // this->headerColMap.at(headerKey).end(), std::string(" ")) << std::endl;
    return this->keyHasColumn(headerKey, dirPraefix + "x" + dirSuffix) &&
           this->keyHasColumn(headerKey, dirPraefix + "y" + dirSuffix) &&
           this->keyHasColumn(headerKey, dirPraefix + "z" + dirSuffix);
  }

  /**
   * @brief Read N timesteps
   *
   * @param start the index to start at reading
   * @param N the nr of groups to read; a negative value results in all groups
   * being read.
   */
  void DumpFileParser::readNGroups(const size_t start, const int N)
  {
    if (!this->file.is_open()) {
      throw std::runtime_error("Cannot read from closed file.");
    }

    if (start >= this->getLength() ||
        (N != -1 && ((static_cast<int>(start)) + N) > this->getLength())) {
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
        // std::cout << "Appending data: " << currentNrOfExpectedGroups << " to
        // "
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
    if (groupsRead != N && N != -1) {
      throw std::runtime_error("Failed to read " + std::to_string(N) +
                               ", read " + std::to_string(groupsRead) +
                               " groups on " + std::to_string(linesRead) +
                               " lines. Stream error: " + std::strerror(errno));
      this->file.clear();
    }
    // std::cout << "Read " << groupsRead << " groups "
    //           << "(last " << start + groupsRead - 1 << ") of " <<
    //           this->getLength()
    //           << std::endl;
  };

  std::vector<long int> DumpFileParser::readTimeSteps()
  {
    return this->readDumpFileSections(ReadableDumpFileSections::TIMESTEP)
      .timesteps;
  };

  std::vector<pylimer_tools::entities::Box> DumpFileParser::readBoxes()
  {
    return this->readDumpFileSections(ReadableDumpFileSections::BOX).boxes;
  };

  constexpr unsigned int str2int(const char* str, int h = 0)
  {
    return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
  }
  constexpr unsigned int str2int(const std::string str, int h = 0)
  {
    return str2int(str.c_str(), h);
  }

  std::vector<std::vector<pylimer_tools::entities::Atom>>
  DumpFileParser::readAtoms()
  {
    return this->readDumpFileSections(ReadableDumpFileSections::ATOM).atoms;
  };

  ReadDumpFileSections DumpFileParser::readDumpFileSections(
    ReadableDumpFileSections sectionsToRead)
  {
    this->rewind();

    std::vector<std::vector<pylimer_tools::entities::Atom>> resultingAtoms;
    resultingAtoms.reserve(this->getLength());

    std::vector<pylimer_tools::entities::Box> resultingBoxes;
    resultingBoxes.reserve(this->getLength());

    std::vector<long int> resultingTimeSteps;
    resultingTimeSteps.reserve(this->getLength());

    std::string line = this->currentLine;
    std::string atomFormat = "";
    int numAtoms = 0;
    while (std::getline(this->file, line)) {
      if (pylimer_tools::utils::startsWith(line, "ITEM: TIMESTEP") &&
          (sectionsToRead & ReadableDumpFileSections::TIMESTEP)) {
        RUNTIME_EXP_IFN(std::getline(this->file, line),
                        "File ended before reading an indicated time-step.");
        resultingTimeSteps.push_back(std::stol(line));
      } else if (pylimer_tools::utils::startsWith(line, "ITEM: BOX BOUNDS") &&
                 (sectionsToRead & ReadableDumpFileSections::BOX)) {
        double loX, hiX, loY, hiY, loZ, hiZ;
        RUNTIME_EXP_IFN(std::getline(this->file, line),
                        "File ended before box could be read");
        // read the next 3 lines, the box
        RUNTIME_EXP_IFN(2 == sscanf(line.c_str(), "%le %le", &loX, &hiX),
                        "Could not read the expected two box coords");
        RUNTIME_EXP_IFN(
          std::getline(this->file, line),
          "File ended before reading box bounds of all three coordinates");
        RUNTIME_EXP_IFN(2 == sscanf(line.c_str(), "%le %le", &loY, &hiY),
                        "Could not read the expected two box coords");
        RUNTIME_EXP_IFN(
          std::getline(this->file, line),
          "File ended before reading box bounds of all three coordinates");
        RUNTIME_EXP_IFN(2 == sscanf(line.c_str(), "%le %le", &loZ, &hiZ),
                        "Could not read the expected two box coords");
        resultingBoxes.push_back(
          pylimer_tools::entities::Box(loX, hiX, loY, hiY, loZ, hiZ));
      } else if (pylimer_tools::utils::startsWith(line,
                                                  "ITEM: NUMBER OF ATOMS") &&
                 (sectionsToRead & ReadableDumpFileSections::ATOM)) {
        RUNTIME_EXP_IFN(std::getline(this->file, line),
                        "File ended before num atoms could be read");
        numAtoms = std::stoi(line);
      } else if (pylimer_tools::utils::startsWith(line, "ITEM: ATOMS") &&
                 (sectionsToRead & ReadableDumpFileSections::ATOM)) {
        atomFormat = pylimer_tools::utils::trimLineOmitComment(
          line.substr(std::string("ITEM: ATOMS ").length()));
        double x, y, z = 0.;
        int nx, ny, nz = 0;
        int id, type;
        std::vector<pylimer_tools::entities::Atom> localResults;
        localResults.reserve(numAtoms);
        for (size_t i = 0; i < numAtoms; ++i) {
          RUNTIME_EXP_IFN(std::getline(this->file, line),
                          "File ended before all atoms could be read");
          // special frequent case
          if (atomFormat == "id type x y z ix iy iz") {
            sscanf(line.c_str(),
                   "%d %d %le %le %le %d %d %d",
                   &id,
                   &type,
                   &x,
                   &y,
                   &z,
                   &nx,
                   &ny,
                   &nz);
          } else {
            std::vector<std::string> splitFormat;
            pylimer_tools::utils::split(splitFormat, atomFormat, " ");
            std::stringstream ss(line);
            for (const std::string& formatPart : splitFormat) {
              switch (str2int(formatPart)) {
                case str2int("id"):
                  ss >> id;
                case str2int("type"):
                  ss >> type;
                case str2int("x"):
                  ss >> x;
                case str2int("y"):
                  ss >> y;
                case str2int("z"):
                  ss >> z;
                case str2int("ix"):
                  ss >> nx;
                case str2int("iy"):
                  ss >> ny;
                case str2int("iz"):
                  ss >> nz;
                default:
                  throw std::runtime_error("Not implemented format part: " +
                                           formatPart);
              }
            }
          }
          localResults.push_back(
            pylimer_tools::entities::Atom(id, type, x, y, z, nx, ny, nz));
        }
        resultingAtoms.push_back(localResults);
      }
    }

    // go back to the beginning of the file
    this->file.clear();
    this->file.seekg(this->groupPosMap.at(0));

    // return the time-steps
    ReadDumpFileSections results;
    results.atoms = resultingAtoms;
    results.boxes = resultingBoxes;
    results.timesteps = resultingTimeSteps;
    return results;
  };

  void DumpFileParser::rewind()
  {
    if (!this->file.is_open()) {
      throw std::runtime_error("Cannot read from closed file.");
    }

    if (this->file.eof()) {
      this->file.clear();
    }
    this->file.seekg(0, std::ios::beg);
  }

  /**
   * @brief Forget the data at a certain index
   *
   * @param index
   */
  void DumpFileParser::forgetAt(const size_t index)
  {
    if (pylimer_tools::utils::map_has_key(this->data, index)) {
      this->data.erase(index);
    }
  };

  /**
   * @brief Read a whole file
   *
   * @param filePath
   */
  void DumpFileParser::read()
  {
    this->data.reserve(this->getLength());
    this->readNGroups(0, -1);
    this->finish();
  }

  void DumpFileParser::finish()
  {
    if (this->file.is_open()) {
      this->file.close();
    }
  }

  std::string DumpFileParser::cleanHeader(std::string headerToClean)
  {
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
    if (!pylimer_tools::utils::map_has_key(this->headerColMap, newHeader)) {
      this->headerColMap.insert_or_assign(newHeader, columns);
    }

    return newHeader;
  }
} // namespace utils
} // namespace pylimer_tools
