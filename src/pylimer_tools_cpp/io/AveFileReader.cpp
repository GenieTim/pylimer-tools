
#include "./AveFileReader.h"

#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <cstring>
#include <filesystem>
#include <fstream> // std::ifstream

namespace pylimer_tools::utils {

std::ifstream
AveFileReader::openFile() const
{
  std::ifstream file(this->filePath);
  if (!file) {
    throw std::runtime_error("Failed open file " + this->filePath);
  }
  return file;
}

void
AveFileReader::readHeaderFromStream(std::ifstream& file,
                                    std::string& firstDataLine,
                                    bool& hasFirstDataLine)
{
  std::string line;
  std::string lastHeaderLine;
  int headerRowCount = 0;
  hasFirstDataLine = false;
  firstDataLine.clear();

  while (std::getline(file, line)) {
    const std::string_view trimmed =
      pylimer_tools::utils::trimWhitespaceView(line);
    if (trimmed.empty()) {
      continue;
    }
    if (trimmed.front() == '#') {
      headerRowCount += 1;
      lastHeaderLine.assign(trimmed.data(), trimmed.size());
    } else {
      firstDataLine = std::move(line);
      hasFirstDataLine = true;
      break;
    }
  }

  if (this->numHeaderRows < 0) {
    this->numHeaderRows = headerRowCount;
  }
  if (this->columnNames.empty() && !lastHeaderLine.empty()) {
    const std::vector<std::string_view> nameViews =
      pylimer_tools::utils::splitWhitespaceView(
        std::string_view(lastHeaderLine).substr(1));
    this->columnNames.reserve(nameViews.size());
    for (const std::string_view nameView : nameViews) {
      this->columnNames.emplace_back(nameView);
    }
  }
}

/**
 * @brief Count the number of lines with data
 *
 * @return int
 */
int
AveFileReader::getNrOfRows()
{
  if (this->numRows >= 0) {
    return this->numRows;
  }
  std::ifstream inStream = this->openFile();
  std::string line;
  bool hasFirstDataLine = false;
  this->readHeaderFromStream(inStream, line, hasFirstDataLine);

  int count = hasFirstDataLine ? 1 : 0;
  while (std::getline(inStream, line)) {
    if (!pylimer_tools::utils::trimLineOmitCommentView(line).empty()) {
      count += 1;
    }
  }

  this->numRows = count;
  return this->numRows;
}

/**
 * @brief Count the number of (comment/header) rows at the top of the file
 *
 * @return int
 */
int
AveFileReader::getNrOfHeaderRows()
{
  if (this->numHeaderRows >= 0) {
    return this->numHeaderRows;
  }

  std::ifstream file = this->openFile();
  std::string firstDataLine;
  bool hasFirstDataLine = false;
  this->readHeaderFromStream(file, firstDataLine, hasFirstDataLine);
  return this->numHeaderRows;
}

/**
 * @brief Find the columns from the header
 *
 * @return const std::vector<std::string>&
 */
const std::vector<std::string>&
AveFileReader::getColumnNames()
{
  if (!this->columnNames.empty()) {
    return this->columnNames;
  }

  std::ifstream file = this->openFile();
  std::string firstDataLine;
  bool hasFirstDataLine = false;
  this->readHeaderFromStream(file, firstDataLine, hasFirstDataLine);
  RUNTIME_EXP_IFN(!this->columnNames.empty(), "Detected header line is empty.");
  return this->columnNames;
}

/**
 * @brief Count the nr of columns in the file
 *
 * @return int
 */
int
AveFileReader::getNrOfColumns()
{
  return static_cast<int>(this->getColumnNames().size());
}

/**
 * @brief Read and cache the file's data if not already loaded
 */
void
AveFileReader::readData()
{
  if (this->dataIsRead) {
    return;
  }

  std::ifstream file = this->openFile();
  std::string line;
  bool hasFirstDataLine = false;
  this->readHeaderFromStream(file, line, hasFirstDataLine);
  RUNTIME_EXP_IFN(!this->columnNames.empty(), "Detected header line is empty.");

  const size_t nCols = this->columnNames.size();
  std::vector<std::vector<double>> results =
    pylimer_tools::utils::initializeWithValue(nCols, std::vector<double>());
  if (this->numRows > 0) {
    const size_t nRows = static_cast<size_t>(this->numRows);
    for (size_t i = 0; i < nCols; ++i) {
      results[i].reserve(nRows);
    }
  }

  size_t rowIdx = 0;
  const auto parseRow = [&](const std::string& rowLine) {
    const std::string_view cleaned =
      pylimer_tools::utils::trimLineOmitCommentView(rowLine);
    if (cleaned.empty()) {
      return;
    }
    const char* cursor = cleaned.data();
    for (size_t col = 0; col < nCols; ++col) {
      double val = 0.0;
      if (pylimer_tools::utils::parseNumericToken(cursor, val)) {
        results[col].push_back(val);
      } else {
        throw std::runtime_error("Failed to read col " + std::to_string(col) +
                                 " on row " + std::to_string(rowIdx) + ".");
      }
    }
    rowIdx += 1;
  };

  if (hasFirstDataLine) {
    parseRow(line);
  }
  while (std::getline(file, line)) {
    parseRow(line);
  }

  if (this->numRows < 0) {
    this->numRows = static_cast<int>(rowIdx);
  }
  this->data = std::move(results);
  this->dataIsRead = true;
}

/**
 * @brief Actually read the file's data
 *
 * @return const std::vector<std::vector<double>>&
 */
const std::vector<std::vector<double>>&
AveFileReader::getData()
{
  if (!this->dataIsRead) {
    this->readData();
  }
  return this->data;
}

namespace {
  std::vector<double> autocorrelateArray(
    const Eigen::Ref<const Eigen::ArrayXd>& colData,
    const Eigen::Index nRows,
    const std::vector<size_t>& dts)
  {
    for (size_t i = 0; i < dts.size(); ++i) {
      if (i > 0) {
        INVALIDARG_EXP_IFN(dts[i - 1] < dts[i],
                           "Invalid dts: They need to be sequential.");
      }
      INVALIDARG_EXP_IFN(nRows > 1 && dts[i] < static_cast<size_t>(nRows - 1),
                         "Invalid dts: got requested " +
                           std::to_string(dts[i]) + ", but only got " +
                           std::to_string(nRows) + " rows.");
    }
    RUNTIME_EXP_IFN(colData.size() == nRows, "Invalid row sizes");

    std::vector<double> results;
    results.reserve(dts.size());
    for (const size_t dt : dts) {
      const Eigen::Index dtIdx = static_cast<Eigen::Index>(dt);
      results.push_back((colData.segment(0, nRows - dtIdx) *
                         colData.segment(dtIdx, nRows - dtIdx))
                          .mean());
    }
    return results;
  }
} // namespace

/**
 * @brief Do autocorrelation for a given column and set of delta indices
 *
 * @param column
 * @param dts
 * @return std::vector<double>
 */
std::vector<double>
AveFileReader::autocorrelateColumn(const int column,
                                   const std::vector<size_t>& dts)
{
  INVALIDINDEX_EXP_IFN(column >= 0, "Invalid column");
  if (!this->dataIsRead) {
    this->readData();
  }
  INVALIDINDEX_EXP_IFN(column < this->getNrOfColumns(), "Invalid column");

  const Eigen::Index nRows = static_cast<Eigen::Index>(this->getNrOfRows());
  const size_t colIdx = static_cast<size_t>(column);
  const Eigen::Map<const Eigen::ArrayXd, Eigen::Unaligned> colData(
    this->data[colIdx].data(),
    static_cast<Eigen::Index>(this->data[colIdx].size()));

  return autocorrelateArray(colData, nRows, dts);
}

/**
 * @brief Do autocorrelation for the difference between two columns and a
 * given set of delta indices
 *
 * @param column1
 * @param column2
 * @param dts
 * @return std::vector<double>
 */
std::vector<double>
AveFileReader::autocorrelateColumnDifference(const int column1,
                                             const int column2,
                                             const std::vector<size_t>& dts)
{
  INVALIDINDEX_EXP_IFN(column1 >= 0 && column2 >= 0, "Invalid column");
  if (!this->dataIsRead) {
    this->readData();
  }
  INVALIDINDEX_EXP_IFN(column1 < this->getNrOfColumns() &&
                         column2 < this->getNrOfColumns(),
                       "Invalid column");

  const Eigen::Index nRows = static_cast<Eigen::Index>(this->getNrOfRows());
  const size_t col1Idx = static_cast<size_t>(column1);
  const size_t col2Idx = static_cast<size_t>(column2);
  const Eigen::ArrayXd colData =
    Eigen::Map<const Eigen::ArrayXd, Eigen::Unaligned>(
      this->data[col1Idx].data(),
      static_cast<Eigen::Index>(this->data[col1Idx].size())) -
    Eigen::Map<const Eigen::ArrayXd, Eigen::Unaligned>(
      this->data[col2Idx].data(),
      static_cast<Eigen::Index>(this->data[col2Idx].size()));

  return autocorrelateArray(colData, nRows, dts);
}

}
