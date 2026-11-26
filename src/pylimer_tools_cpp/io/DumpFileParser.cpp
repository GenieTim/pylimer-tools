#include "DumpFileParser.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include <any>
#include <cassert>
#include <cstring>
#include <filesystem>
#include <fstream> // std::ifstream
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace pylimer_tools::utils {
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
  this->file.open(src.filePath, std::ios::binary);
};
// 3. copy assignment operator
DumpFileParser&
DumpFileParser::operator=(DumpFileParser src)
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
  this->file.open(src.filePath, std::ios::binary);
  src.file.open(this->filePath, std::ios::binary);
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

  this->openFile();
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
void
DumpFileParser::readGroupByIdx(const size_t i)
{
  this->readNGroups(i, 1);
}

/**
 * @brief Get the nr of groups
 *
 * @return size_t
 */
size_t
DumpFileParser::getLength() const
{
  return this->nrOfGroups;
}

/**
 * @brief Check whether a header key exists
 *
 * @param headerKey
 * @return bool
 */
bool
DumpFileParser::hasKey(std::string headerKey) const
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
bool
DumpFileParser::keyHasColumn(std::string headerKey, std::string column)
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
bool
DumpFileParser::keyHasDirectionalColumn(std::string headerKey,
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
void
DumpFileParser::readNGroups(const size_t start, const int N)
{
  if (!this->file.is_open()) {
    throw std::runtime_error("Cannot read from closed file.");
  }

  if (start >= this->getLength() ||
      (N != -1 &&
       (N < 0 || (start + static_cast<size_t>(N)) > this->getLength()))) {
    throw std::invalid_argument(
      "Cannot read from outside the length of the "
      "dump file. Tried to read from " +
      std::to_string(start) + " to " + std::to_string(N) + " for a file with " +
      std::to_string(this->getLength()) + " time-steps.");
  }

  if (this->file.eof()) {
    this->file.clear();
  }
  this->file.seekg(this->groupPosMap.at(start));

  size_t groupsRead = 0;
  std::string currentKey = this->cleanHeader(this->newGroupKey);
  size_t currentNrOfExpectedGroups = this->headerColMap.at(currentKey).size();
  data_item_t dataItem;
  std::string line = this->currentLine;
  // std::cout << "Starting to read at " << start << " for " << N << " with "
  //           << line << " and key " << currentKey << std::endl;
  int linesRead = 0;

  while (std::getline(this->file, line)) {
    linesRead += 1;
    // std::cout << "Read line: " << line << std::endl;
    line = pylimer_tools::utils::rtrim(
      pylimer_tools::utils::trimLineOmitComment(line));
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
      if (N > 0 && (groupsRead + 1) >= static_cast<size_t>(N)) {
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
  if (N != -1 && groupsRead != static_cast<size_t>(N)) {
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

std::vector<long int>
DumpFileParser::readTimeSteps()
{
  return this->readDumpFileSections(ReadableDumpFileSections::TIMESTEP)
    .timesteps;
};

std::vector<pylimer_tools::entities::Box>
DumpFileParser::readBoxes()
{
  return this->readDumpFileSections(ReadableDumpFileSections::BOX).boxes;
};

constexpr unsigned int
str2int(const char* str, const int h = 0)
{
  return !str[h]
           ? 5381
           : (str2int(str, h + 1) * 33) ^
               static_cast<unsigned int>(static_cast<unsigned char>(str[h]));
}
// constexpr
unsigned int
str2int(const std::string& str, const int h = 0)
{
  return str2int(str.c_str(), h);
}

std::vector<std::vector<pylimer_tools::entities::Atom>>
DumpFileParser::readAtoms()
{
  return this->readDumpFileSections(ReadableDumpFileSections::ATOM).atoms;
};

ReadDumpFileSections
DumpFileParser::readDumpFileSections(ReadableDumpFileSections sectionsToRead)
{
  this->rewind();

  std::vector<std::vector<pylimer_tools::entities::Atom>> resultingAtoms;
  resultingAtoms.reserve(this->getLength());

  std::vector<std::unordered_map<std::string, std::vector<double>>>
    additionalAtomData;
  if (sectionsToRead & ReadableDumpFileSections::EXTRA_ATOM) {
    additionalAtomData.reserve(this->getLength());
  }

  std::vector<pylimer_tools::entities::Box> resultingBoxes;
  resultingBoxes.reserve(this->getLength());

  std::vector<long int> resultingTimeSteps;
  resultingTimeSteps.reserve(this->getLength());

  std::string line = this->currentLine;
  std::string atomFormat = "";
  int numAtoms = 0;
  size_t sectionsRead = 0;
  while (std::getline(this->file, line)) {
    const std::string_view trimmedHeader =
      pylimer_tools::utils::trimLineOmitCommentView(line);
    if (trimmedHeader.starts_with("ITEM: TIMESTEP") &&
        (sectionsToRead & ReadableDumpFileSections::TIMESTEP)) {
      RUNTIME_EXP_IFN(std::getline(this->file, line),
                      "File ended before reading an indicated time-step.");
      resultingTimeSteps.push_back(std::stol(line));
    } else if (trimmedHeader.starts_with("ITEM: BOX BOUNDS") &&
               ((sectionsToRead & ReadableDumpFileSections::BOX) ||
                (sectionsToRead & ReadableDumpFileSections::ATOM))) {
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
    } else if (trimmedHeader.starts_with("ITEM: NUMBER OF ATOMS") &&
               (sectionsToRead & ReadableDumpFileSections::ATOM)) {
      RUNTIME_EXP_IFN(std::getline(this->file, line),
                      "File ended before num atoms could be read");
      numAtoms = std::stoi(line);
    } else if (trimmedHeader.starts_with("ITEM: ATOMS") &&
               (sectionsToRead & ReadableDumpFileSections::ATOM)) {
      constexpr std::string_view atomsPrefix = "ITEM: ATOMS";
      const std::string_view formatView =
        trimmedHeader.substr(atomsPrefix.size());
      const size_t firstChar =
        formatView.find_first_not_of(pylimer_tools::utils::WHITESPACE_VIEW);
      const size_t lastChar =
        formatView.find_last_not_of(pylimer_tools::utils::WHITESPACE_VIEW);
      atomFormat =
        (firstChar == std::string_view::npos)
          ? std::string()
          : std::string(formatView.substr(firstChar, lastChar - firstChar + 1));
      double x = 0., y = 0, z = 0.;
      int nx = 0, ny = 0, nz = 0;
      bool isUnwrappedX = false, isUnwrappedY = false, isUnwrappedZ = false;
      bool isScaledX = false, isScaledY = false, isScaledZ = false;
      int id = 0, type = 1;
      std::vector<pylimer_tools::entities::Atom> localResults;
      localResults.reserve(static_cast<size_t>(numAtoms));
      std::unordered_map<std::string, std::vector<double>> localExtraAtomData;
      for (int i = 0; i < numAtoms; ++i) {
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
#define CASE(idStr, targetVar)                                                 \
  case str2int(idStr):                                                         \
    ss >> targetVar;                                                           \
    /** always write to the extra data, for all properties */                  \
    if (sectionsToRead & ReadableDumpFileSections::EXTRA_ATOM) {               \
      localExtraAtomData[formatPart].push_back(                                \
        static_cast<double>(targetVar));                                       \
    }

            switch (str2int(formatPart)) {
              CASE("id", id)
              break;
              CASE("type", type)
              break;
              CASE("x", x)
              break;
              CASE("xu", x)
              isUnwrappedX = true;
              break;
              CASE("xs", x)
              isScaledX = true;
              break;
              CASE("xsu", x)
              isUnwrappedX = true;
              isScaledX = true;
              break;
              CASE("y", y)
              break;
              CASE("yu", y)
              isUnwrappedY = true;
              break;
              CASE("ys", y)
              isScaledY = true;
              break;
              CASE("ysu", y)
              isUnwrappedY = true;
              isScaledY = true;
              break;
              CASE("z", z)
              break;
              CASE("zu", z)
              isUnwrappedZ = true;
              break;
              CASE("zs", z)
              isScaledZ = true;
              break;
              CASE("zsu", z)
              isUnwrappedZ = true;
              isScaledZ = true;
              break;
              CASE("ix", nx)
              break;
              CASE("iy", ny)
              break;
              CASE("iz", nz)
              break;
              default:
                if (sectionsToRead & ReadableDumpFileSections::EXTRA_ATOM) {
                  double d;
                  ss >> d;
                  localExtraAtomData[formatPart].push_back(d);
                }
                // throw std::runtime_error("Not implemented format part: '" +
                //                          formatPart + "'");
            }

#undef CASE
          }
        }
        assert(resultingBoxes.size() > sectionsRead);
        localResults.push_back(pylimer_tools::entities::Atom(
          id,
          type,
          x * (isScaledX ? resultingBoxes[sectionsRead].getLx() : 1.),
          y * (isScaledY ? resultingBoxes[sectionsRead].getLy() : 1.),
          z * (isScaledZ ? resultingBoxes[sectionsRead].getLz() : 1.),
          isUnwrappedX ? 0 : nx,
          isUnwrappedY ? 0 : ny,
          isUnwrappedZ ? 0 : nz));
      }
      resultingAtoms.push_back(localResults);
      if (sectionsToRead & ReadableDumpFileSections::EXTRA_ATOM) {
        additionalAtomData.push_back(localExtraAtomData);
      }

      sectionsRead += 1;
      if ((sectionsToRead & ReadableDumpFileSections::TIMESTEP)) {
        assert(resultingTimeSteps.size() == sectionsRead);
      }
      if ((sectionsToRead & ReadableDumpFileSections::BOX)) {
        assert(resultingBoxes.size() == sectionsRead);
      }
      if ((sectionsToRead & ReadableDumpFileSections::ATOM)) {
        assert(resultingAtoms.size() == sectionsRead);
      }
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
  results.extraAtomsData = additionalAtomData;
  return results;
};

void
DumpFileParser::rewind()
{
  if (!this->file.is_open()) {
    this->openFile();
  }

  if (this->file.eof()) {
    this->file.clear();
  }
  this->file.seekg(0, std::ios::beg);
}

namespace {
  constexpr size_t DUMP_SCAN_BUFFER_SIZE = 256 * 1024;

  inline std::string_view trimDumpKeyView(const std::string_view line)
  {
    std::string_view trimmed =
      pylimer_tools::utils::trimLineOmitCommentView(line);
    const size_t end =
      trimmed.find_last_not_of(pylimer_tools::utils::WHITESPACE_VIEW);
    return (end == std::string_view::npos) ? std::string_view()
                                           : trimmed.substr(0, end + 1);
  }

  inline bool lineMatchesKey(const std::string_view line,
                             const std::string& key)
  {
    const size_t keyLen = key.size();
    const size_t lineLen = line.size();
    if (lineLen < keyLen) {
      return false;
    }
    if (lineLen == keyLen) {
      return std::memcmp(line.data(), key.data(), keyLen) == 0;
    }
    const char first = line[0];
    if (first != key[0] && first != ' ' && first != '\t' && first != '\r' &&
        first != '\n' && first != '\f' && first != '\v') {
      return false;
    }
    return trimDumpKeyView(line) == key;
  }
} // namespace

void
DumpFileParser::openFile()
{
  this->file.open(this->filePath, std::ios::binary);

  if (!this->file.is_open()) {
    throw std::invalid_argument("File to read ('" + this->filePath +
                                "'): failed to open.");
  }

  this->groupPosMap.clear();

  std::vector<char> buffer(DUMP_SCAN_BUFFER_SIZE);
  std::string carry;
  std::streamoff chunkStartOffset = 0;
  bool foundFirstKey = false;
  size_t groupsFound = 0;

  while (this->file) {
    this->file.read(buffer.data(),
                    static_cast<std::streamsize>(DUMP_SCAN_BUFFER_SIZE));
    const std::streamsize bytesRead = this->file.gcount();
    if (bytesRead <= 0) {
      break;
    }

    const char* ptr = buffer.data();
    const char* const end = buffer.data() + bytesRead;

    while (ptr < end) {
      const char* const newline = static_cast<const char*>(
        std::memchr(ptr, '\n', static_cast<size_t>(end - ptr)));
      if (newline == nullptr) {
        carry.append(ptr, static_cast<size_t>(end - ptr));
        break;
      }

      std::string_view lineView;
      if (carry.empty()) {
        lineView = std::string_view(ptr, static_cast<size_t>(newline - ptr));
      } else {
        carry.append(ptr, static_cast<size_t>(newline - ptr));
        lineView = std::string_view(carry);
      }

      const std::streampos nextLinePos = static_cast<std::streampos>(
        chunkStartOffset + (newline - buffer.data()) + 1);

      if (!foundFirstKey) {
        const std::string_view trimmed = trimDumpKeyView(lineView);
        if (!trimmed.empty()) {
          // new group key: key for a new timestep (group)
          this->newGroupKey.assign(trimmed.data(), trimmed.size());
          this->currentLine = this->newGroupKey;
          // record position of index to jump back at some point
          this->groupPosMap.emplace(0, nextLinePos);
          foundFirstKey = true;
        }
      } else if (lineMatchesKey(lineView, this->newGroupKey)) {
        // new timestep
        groupsFound += 1;
        this->groupPosMap.emplace(groupsFound, nextLinePos);
      }

      if (!carry.empty()) {
        carry.clear();
      }
      ptr = newline + 1;
    }

    chunkStartOffset += bytesRead;
  }

  if (!carry.empty()) {
    const std::streampos eofPos = static_cast<std::streampos>(chunkStartOffset);
    if (!foundFirstKey) {
      const std::string_view trimmed = trimDumpKeyView(carry);
      if (!trimmed.empty()) {
        this->newGroupKey.assign(trimmed.data(), trimmed.size());
        this->currentLine = this->newGroupKey;
        this->groupPosMap.emplace(0, eofPos);
        foundFirstKey = true;
      }
    } else if (lineMatchesKey(carry, this->newGroupKey)) {
      groupsFound += 1;
      this->groupPosMap.emplace(groupsFound, eofPos);
    }
  }

  if (!foundFirstKey) {
    this->newGroupKey.clear();
    this->currentLine.clear();
    this->groupPosMap.emplace(0, static_cast<std::streampos>(chunkStartOffset));
  }

  this->nrOfGroups = groupsFound + 1;
  // reset position to start of first group
  this->file.clear();
  this->file.seekg(this->groupPosMap.at(0));
}

/**
 * @brief Forget the data at a certain index
 *
 * @param index
 */
void
DumpFileParser::forgetAt(const size_t index)
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
void
DumpFileParser::read()
{
  this->data.reserve(this->getLength());
  this->readNGroups(0, -1);
  this->finish();
}

void
DumpFileParser::finish()
{
  if (this->file.is_open()) {
    this->file.close();
  }
}

std::string
DumpFileParser::cleanHeader(std::string headerToClean)
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
}
