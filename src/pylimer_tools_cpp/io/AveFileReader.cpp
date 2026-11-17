
#include "./AveFileReader.h"

#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <cstring>
#include <filesystem>
#include <fstream> // std::ifstream

namespace pylimer_tools::utils {

/**
 * @brief Count the number of lines with data
 *
 * @return int
 */
int
AveFileReader::getNrOfRows()
{
  if (this->numRows > 0) {
    return this->numRows;
  }
  std::ifstream in_stream(this->filePath);
  if (!in_stream) {
    throw std::runtime_error("Failed open file " + this->filePath);
  }

  int count = 0;
  std::string line;

  // skip all header lines
  while (getline(in_stream, line)) {
    if (line.at(0) != '#') {
      count += line.size() > 0;
      break;
    }
  }
  while (getline(in_stream, line)) {
    if (line.size() > 0) {
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
  if (this->numHeaderRows > 0) {
    return this->numHeaderRows;
  }

  std::ifstream file(this->filePath);
  if (!file) {
    throw std::runtime_error("Failed open file " + this->filePath);
  }

  std::string line;
  int count = 0;
  while (getline(file, line)) {
    if (line.at(0) == '#') {
      count += 1;
    } else {
      break;
    }
  }

  this->numHeaderRows = count;
  return count;
}

/**
 * @brief Find the columns from the header
 *
 * @return std::vector<std::string>
 */
std::vector<std::string>
AveFileReader::getColumnNames()
{
  if (this->columnNames.size() > 0) {
    return this->columnNames;
  }

  std::ifstream file(this->filePath);
  if (!file) {
    throw std::runtime_error("Failed open file " + this->filePath);
  }

  std::string line;
  std::string lastLine;
  while (getline(file, line)) {
    if (line.size() == 0) {
      break;
    }
    if (line.at(0) == '#') {
      lastLine = line;
    } else {
      break;
    }
  }

  RUNTIME_EXP_IFN(lastLine.size() > 0, "Detected header line is empty.");

  std::string headerLine =
    pylimer_tools::utils::trim(lastLine.substr(1, lastLine.size() - 1));
  this->columnNames = pylimer_tools::utils::split(headerLine, ' ');
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
  return this->getColumnNames().size();
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

  const size_t nRows = static_cast<size_t>(std::max(0, this->getNrOfRows()));
  const size_t nCols = static_cast<size_t>(std::max(0, this->getNrOfColumns()));
  std::vector<std::vector<double>> results =
    pylimer_tools::utils::initializeWithValue(nCols, std::vector<double>());
  for (size_t i = 0; i < nCols; ++i) {
    results[i].reserve(nCols);
  }

  std::ifstream file(this->filePath);
  if (!file) {
    throw std::runtime_error("Failed open file " + this->filePath);
  }

  std::string line;
  const size_t nHeaderRows =
    static_cast<size_t>(std::max(0, this->getNrOfHeaderRows()));
  for (size_t i = 0; i < nHeaderRows; ++i) {
    RUNTIME_EXP_IFN(std::getline(file, line),
                    "File ended before data was even reached");
  }

  for (size_t i = 0; i < nRows; ++i) {
    if (!std::getline(file, line)) {
      if (i + 1 < nRows) {
        throw std::runtime_error(
          "File ended before all rows could be read (reached row " +
          std::to_string(i) + " of " + std::to_string(nRows) + ")");
      } else {
        break;
      }
    }
    std::stringstream ss(line);
    for (size_t col = 0; col < nCols; ++col) {
      double val;
      if (ss >> val) {
        results[col].push_back(val);
      } else {
        throw std::runtime_error("Failed to read col " + std::to_string(col) +
                                 " on row " + std::to_string(i) + ".");
      }
    }
  }

  this->data = std::move(results);
  this->dataIsRead = true;
}

/**
 * @brief Actually read the file's data
 *
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>>
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
  INVALIDINDEX_EXP_IFN(column >= 0 && column < this->getNrOfColumns(),
                       "Invalid column");

  const Eigen::Index nRows = static_cast<Eigen::Index>(this->getNrOfRows());
  if (!this->dataIsRead) {
    this->readData();
  }
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
  INVALIDINDEX_EXP_IFN(column1 >= 0 && column1 < this->getNrOfColumns(),
                       "Invalid column");
  INVALIDINDEX_EXP_IFN(column2 >= 0 && column2 < this->getNrOfColumns(),
                       "Invalid column");

  const Eigen::Index nRows = static_cast<Eigen::Index>(this->getNrOfRows());
  if (!this->dataIsRead) {
    this->readData();
  }
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
