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
  : filePath(src.filePath)
  , newGroupKey(src.newGroupKey)
  , currentLine(src.currentLine)
  , nrOfGroups(src.nrOfGroups)
  , data(src.data)
  , headerColMap(src.headerColMap)
  , groupPosMap(src.groupPosMap)
{
  if (!src.filePath.empty() && src.file.is_open()) {
    this->file.open(src.filePath, std::ios::binary);
    if (this->file.is_open() && !this->groupPosMap.empty()) {
      this->file.seekg(this->groupPosMap.at(0));
    }
  }
};
// 3. copy assignment operator
DumpFileParser&
DumpFileParser::operator=(DumpFileParser src)
{
  std::swap(this->filePath, src.filePath);
  std::swap(this->newGroupKey, src.newGroupKey);
  std::swap(this->currentLine, src.currentLine);
  std::swap(this->file, src.file);
  std::swap(this->nrOfGroups, src.nrOfGroups);
  std::swap(this->data, src.data);
  std::swap(this->headerColMap, src.headerColMap);
  std::swap(this->groupPosMap, src.groupPosMap);
  return *this;
};

/**
 * @brief Initialize the parser to read from a certain file path
 *
 * @param filePath
 */
DumpFileParser::DumpFileParser(const std::string& filePath)
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
DumpFileParser::hasKey(const std::string& headerKey) const
{
  if (this->data.empty()) {
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
DumpFileParser::keyHasColumn(const std::string& headerKey,
                             const std::string& column) const
{
  const std::vector<std::string>& cols = this->headerColMap.at(headerKey);
  return std::find(cols.begin(), cols.end(), column) != cols.end();
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
DumpFileParser::keyHasDirectionalColumn(const std::string& headerKey,
                                        const std::string& dirPraefix,
                                        const std::string& dirSuffix) const
{
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
      this->data.insert_or_assign(start + groupsRead, std::move(dataItem));
      groupsRead += 1;
      dataItem.clear();
    }
  }

  if (!dataItem.empty()) {
    // last timestep
    this->data.insert_or_assign(start + groupsRead, std::move(dataItem));
    groupsRead += 1;
  }
  if (N != -1 && groupsRead != static_cast<size_t>(N)) {
    throw std::runtime_error("Failed to read " + std::to_string(N) + ", read " +
                             std::to_string(groupsRead) + " groups on " +
                             std::to_string(linesRead) +
                             " lines. Stream error: " + std::strerror(errno));
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

namespace {
  enum class DumpAtomField
  {
    ID,
    TYPE,
    X,
    XU,
    XS,
    XSU,
    Y,
    YU,
    YS,
    YSU,
    Z,
    ZU,
    ZS,
    ZSU,
    IX,
    IY,
    IZ,
    EXTRA
  };

  struct DumpAtomColumn
  {
    DumpAtomField field;
    std::vector<double>* extraVec = nullptr;
  };

  DumpAtomField parseDumpAtomField(const std::string_view name)
  {
    if (name == "id") {
      return DumpAtomField::ID;
    }
    if (name == "type") {
      return DumpAtomField::TYPE;
    }
    if (name == "x") {
      return DumpAtomField::X;
    }
    if (name == "xu") {
      return DumpAtomField::XU;
    }
    if (name == "xs") {
      return DumpAtomField::XS;
    }
    if (name == "xsu") {
      return DumpAtomField::XSU;
    }
    if (name == "y") {
      return DumpAtomField::Y;
    }
    if (name == "yu") {
      return DumpAtomField::YU;
    }
    if (name == "ys") {
      return DumpAtomField::YS;
    }
    if (name == "ysu") {
      return DumpAtomField::YSU;
    }
    if (name == "z") {
      return DumpAtomField::Z;
    }
    if (name == "zu") {
      return DumpAtomField::ZU;
    }
    if (name == "zs") {
      return DumpAtomField::ZS;
    }
    if (name == "zsu") {
      return DumpAtomField::ZSU;
    }
    if (name == "ix") {
      return DumpAtomField::IX;
    }
    if (name == "iy") {
      return DumpAtomField::IY;
    }
    if (name == "iz") {
      return DumpAtomField::IZ;
    }
    return DumpAtomField::EXTRA;
  }
} // namespace

std::vector<std::vector<pylimer_tools::entities::Atom>>
DumpFileParser::readAtoms()
{
  return this->readDumpFileSections(ReadableDumpFileSections::ATOM).atoms;
};

ReadDumpFileSections
DumpFileParser::readDumpFileSections(ReadableDumpFileSections sectionsToRead)
{
  ReadDumpFileSections results;
  if (this->getLength() == 0) {
    return results;
  }

  this->rewind();

  const size_t expectedFrames = this->getLength();
  if (sectionsToRead & ReadableDumpFileSections::ATOM) {
    results.atoms.reserve(expectedFrames);
  }
  if (sectionsToRead & ReadableDumpFileSections::EXTRA_ATOM) {
    results.extraAtomsData.reserve(expectedFrames);
  }
  if ((sectionsToRead & ReadableDumpFileSections::BOX) ||
      (sectionsToRead & ReadableDumpFileSections::ATOM)) {
    results.boxes.reserve(expectedFrames);
  }
  if (sectionsToRead & ReadableDumpFileSections::TIMESTEP) {
    results.timesteps.reserve(expectedFrames);
  }

  std::string line = this->currentLine;
  int numAtoms = 0;
  size_t sectionsRead = 0;
  const bool readExtraAtom =
    static_cast<bool>(sectionsToRead & ReadableDumpFileSections::EXTRA_ATOM);
  while (std::getline(this->file, line)) {
    const std::string_view trimmedHeader =
      pylimer_tools::utils::trimLineOmitCommentView(line);
    if (trimmedHeader.starts_with("ITEM: TIMESTEP") &&
        (sectionsToRead & ReadableDumpFileSections::TIMESTEP)) {
      RUNTIME_EXP_IFN(std::getline(this->file, line),
                      "File ended before reading an indicated time-step.");
      const char* cursor = line.c_str();
      long int timestep = 0;
      RUNTIME_EXP_IFN(pylimer_tools::utils::parseNumericToken(cursor, timestep),
                      "Could not parse timestep from line '" + line + "'.");
      results.timesteps.push_back(timestep);
    } else if (trimmedHeader.starts_with("ITEM: BOX BOUNDS") &&
               ((sectionsToRead & ReadableDumpFileSections::BOX) ||
                (sectionsToRead & ReadableDumpFileSections::ATOM))) {
      double loX = 0.0, hiX = 0.0, loY = 0.0, hiY = 0.0, loZ = 0.0, hiZ = 0.0;
      RUNTIME_EXP_IFN(std::getline(this->file, line),
                      "File ended before box could be read");
      const char* pX = line.c_str();
      RUNTIME_EXP_IFN(pylimer_tools::utils::parseNumericToken(pX, loX) &&
                        pylimer_tools::utils::parseNumericToken(pX, hiX),
                      "Could not read the expected two box coords");
      RUNTIME_EXP_IFN(
        std::getline(this->file, line),
        "File ended before reading box bounds of all three coordinates");
      const char* pY = line.c_str();
      RUNTIME_EXP_IFN(pylimer_tools::utils::parseNumericToken(pY, loY) &&
                        pylimer_tools::utils::parseNumericToken(pY, hiY),
                      "Could not read the expected two box coords");
      RUNTIME_EXP_IFN(
        std::getline(this->file, line),
        "File ended before reading box bounds of all three coordinates");
      const char* pZ = line.c_str();
      RUNTIME_EXP_IFN(pylimer_tools::utils::parseNumericToken(pZ, loZ) &&
                        pylimer_tools::utils::parseNumericToken(pZ, hiZ),
                      "Could not read the expected two box coords");
      results.boxes.emplace_back(loX, hiX, loY, hiY, loZ, hiZ);
    } else if (trimmedHeader.starts_with("ITEM: NUMBER OF ATOMS") &&
               (sectionsToRead & ReadableDumpFileSections::ATOM)) {
      RUNTIME_EXP_IFN(std::getline(this->file, line),
                      "File ended before num atoms could be read");
      const char* cursor = line.c_str();
      RUNTIME_EXP_IFN(
        pylimer_tools::utils::parseNumericToken(cursor, numAtoms) &&
          numAtoms >= 0,
        "Could not parse valid number of atoms from line '" + line + "'.");
    } else if (trimmedHeader.starts_with("ITEM: ATOMS") &&
               (sectionsToRead & ReadableDumpFileSections::ATOM)) {
      constexpr std::string_view atomsPrefix = "ITEM: ATOMS";
      const std::string_view trimmedFormat =
        pylimer_tools::utils::trimWhitespaceView(
          trimmedHeader.substr(atomsPrefix.size()));

      bool isUnwrappedX = false, isUnwrappedY = false, isUnwrappedZ = false;
      bool isScaledX = false, isScaledY = false, isScaledZ = false;
      std::vector<pylimer_tools::entities::Atom> localResults;
      localResults.reserve(static_cast<size_t>(numAtoms));
      std::unordered_map<std::string, std::vector<double>> localExtraAtomData;

      const bool useDefaultFastPath =
        (trimmedFormat == "id type x y z ix iy iz") && !readExtraAtom;
      std::vector<DumpAtomColumn> columns;
      if (!useDefaultFastPath) {
        const std::vector<std::string_view> colViews =
          pylimer_tools::utils::splitWhitespaceView(trimmedFormat);
        columns.reserve(colViews.size());
        for (const std::string_view colView : colViews) {
          const DumpAtomField field = parseDumpAtomField(colView);
          switch (field) {
            case DumpAtomField::XU:
              isUnwrappedX = true;
              break;
            case DumpAtomField::XS:
              isScaledX = true;
              break;
            case DumpAtomField::XSU:
              isUnwrappedX = true;
              isScaledX = true;
              break;
            case DumpAtomField::YU:
              isUnwrappedY = true;
              break;
            case DumpAtomField::YS:
              isScaledY = true;
              break;
            case DumpAtomField::YSU:
              isUnwrappedY = true;
              isScaledY = true;
              break;
            case DumpAtomField::ZU:
              isUnwrappedZ = true;
              break;
            case DumpAtomField::ZS:
              isScaledZ = true;
              break;
            case DumpAtomField::ZSU:
              isUnwrappedZ = true;
              isScaledZ = true;
              break;
            default:
              break;
          }

          std::vector<double>* extraVec = nullptr;
          if (readExtraAtom) {
            extraVec = &localExtraAtomData[std::string(colView)];
            extraVec->reserve(static_cast<size_t>(numAtoms));
          }
          columns.push_back(DumpAtomColumn{ field, extraVec });
        }
      }

      const bool anyScaled = isScaledX || isScaledY || isScaledZ;
      const bool hasBox = results.boxes.size() > sectionsRead;
      RUNTIME_EXP_IFN(
        !anyScaled || hasBox,
        "Scaled atom coordinates require BOX BOUNDS in dump file.");
      const double scaleX =
        (isScaledX && hasBox) ? results.boxes[sectionsRead].getLx() : 1.0;
      const double scaleY =
        (isScaledY && hasBox) ? results.boxes[sectionsRead].getLy() : 1.0;
      const double scaleZ =
        (isScaledZ && hasBox) ? results.boxes[sectionsRead].getLz() : 1.0;
      const double offsetX =
        (isScaledX && hasBox) ? results.boxes[sectionsRead].getLowX() : 0.0;
      const double offsetY =
        (isScaledY && hasBox) ? results.boxes[sectionsRead].getLowY() : 0.0;
      const double offsetZ =
        (isScaledZ && hasBox) ? results.boxes[sectionsRead].getLowZ() : 0.0;

      for (int i = 0; i < numAtoms; ++i) {
        RUNTIME_EXP_IFN(std::getline(this->file, line),
                        "File ended before all atoms could be read");
        const char* cursor = line.c_str();
        long int id = 0;
        int type = 1;
        double x = 0.0, y = 0.0, z = 0.0;
        int nx = 0, ny = 0, nz = 0;
        bool ok = true;

        // special frequent case
        if (useDefaultFastPath) {
          ok = pylimer_tools::utils::parseNumericToken(cursor, id) &&
               pylimer_tools::utils::parseNumericToken(cursor, type) &&
               pylimer_tools::utils::parseNumericToken(cursor, x) &&
               pylimer_tools::utils::parseNumericToken(cursor, y) &&
               pylimer_tools::utils::parseNumericToken(cursor, z) &&
               pylimer_tools::utils::parseNumericToken(cursor, nx) &&
               pylimer_tools::utils::parseNumericToken(cursor, ny) &&
               pylimer_tools::utils::parseNumericToken(cursor, nz);
        } else {
          for (const DumpAtomColumn& col : columns) {
            double extraVal = 0.0;
            switch (col.field) {
              case DumpAtomField::ID:
                ok = pylimer_tools::utils::parseNumericToken(cursor, id);
                extraVal = static_cast<double>(id);
                break;
              case DumpAtomField::TYPE:
                ok = pylimer_tools::utils::parseNumericToken(cursor, type);
                extraVal = static_cast<double>(type);
                break;
              case DumpAtomField::X:
              case DumpAtomField::XU:
              case DumpAtomField::XS:
              case DumpAtomField::XSU:
                ok = pylimer_tools::utils::parseNumericToken(cursor, x);
                extraVal = x;
                break;
              case DumpAtomField::Y:
              case DumpAtomField::YU:
              case DumpAtomField::YS:
              case DumpAtomField::YSU:
                ok = pylimer_tools::utils::parseNumericToken(cursor, y);
                extraVal = y;
                break;
              case DumpAtomField::Z:
              case DumpAtomField::ZU:
              case DumpAtomField::ZS:
              case DumpAtomField::ZSU:
                ok = pylimer_tools::utils::parseNumericToken(cursor, z);
                extraVal = z;
                break;
              case DumpAtomField::IX:
                ok = pylimer_tools::utils::parseNumericToken(cursor, nx);
                extraVal = static_cast<double>(nx);
                break;
              case DumpAtomField::IY:
                ok = pylimer_tools::utils::parseNumericToken(cursor, ny);
                extraVal = static_cast<double>(ny);
                break;
              case DumpAtomField::IZ:
                ok = pylimer_tools::utils::parseNumericToken(cursor, nz);
                extraVal = static_cast<double>(nz);
                break;
              case DumpAtomField::EXTRA:
              default:
                ok = pylimer_tools::utils::parseNumericToken(cursor, extraVal);
                break;
            }
            if (!ok) {
              break;
            }
            if (col.extraVec != nullptr) {
              col.extraVec->push_back(extraVal);
            }
          }
        }
        RUNTIME_EXP_IFN(ok, "Failed to parse atom line '" + line + "'.");

        localResults.emplace_back(id,
                                  type,
                                  offsetX + x * scaleX,
                                  offsetY + y * scaleY,
                                  offsetZ + z * scaleZ,
                                  isUnwrappedX ? 0 : nx,
                                  isUnwrappedY ? 0 : ny,
                                  isUnwrappedZ ? 0 : nz);
      }
      results.atoms.push_back(std::move(localResults));
      if (readExtraAtom) {
        results.extraAtomsData.push_back(std::move(localExtraAtomData));
      }

      sectionsRead += 1;
      if ((sectionsToRead & ReadableDumpFileSections::TIMESTEP)) {
        assert(results.timesteps.size() == sectionsRead);
      }
      if ((sectionsToRead & ReadableDumpFileSections::BOX)) {
        assert(results.boxes.size() == sectionsRead);
      }
      if ((sectionsToRead & ReadableDumpFileSections::ATOM)) {
        assert(results.atoms.size() == sectionsRead);
      }
    }
  }

  // go back to the beginning of the file
  this->file.clear();
  this->file.seekg(this->groupPosMap.at(0));

  if (!(sectionsToRead & ReadableDumpFileSections::BOX)) {
    results.boxes.clear();
  }
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
    this->nrOfGroups = 0;
    this->groupPosMap.emplace(0, static_cast<std::streampos>(0));
  } else {
    this->nrOfGroups = groupsFound + 1;
  }
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
  if (this->getLength() == 0) {
    this->finish();
    return;
  }
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
  if (headerToClean.size() >= 6) {
    headerToClean.erase(0, 6);
  } else {
    headerToClean.clear();
  }
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
