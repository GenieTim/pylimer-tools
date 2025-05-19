
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
 * @brief Actually read the file's data
 *
 * @return std::vector<std::vector<double>>
 */
std::vector<std::vector<double>>
AveFileReader::getData()
{
  if (this->data.size() > 0) {
    return this->data;
  }

  int nRows = this->getNrOfRows();
  int nCols = this->getNrOfColumns();
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
  for (size_t i = 0; i < this->getNrOfHeaderRows(); ++i) {
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

  this->data = results;
  return results;
}

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
  INVALIDINDEX_EXP_IFN(column < this->getNrOfColumns(), "Invalid column");

  int nRows = this->getNrOfRows();
  // validate dts
  for (size_t i = 1; i < dts.size(); ++i) {
    INVALIDARG_EXP_IFN(dts[i - 1] < dts[i],
                       "Invalid dts: They need to be sequential.");
    INVALIDARG_EXP_IFN(dts[i] < nRows - 1,
                       "Invalid dts: got requested " + std::to_string(dts[i]) +
                         ", but only got " + std::to_string(nRows) + " rows.");
  }

  this->getData();
  Eigen::ArrayXd colData = Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(
    this->data[column].data(), this->data[column].size());
  RUNTIME_EXP_IFN(colData.size() == nRows, "Invalid row sizes");

  std::vector<double> results;
  results.reserve(dts.size());
  for (size_t dt : dts) {
    results.push_back(
      (colData.segment(0, nRows - dt) * colData.segment(dt, nRows - dt))
        .mean());
  }

  return results;
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
  INVALIDINDEX_EXP_IFN(column1 < this->getNrOfColumns(), "Invalid column");
  INVALIDINDEX_EXP_IFN(column2 < this->getNrOfColumns(), "Invalid column");

  int nRows = this->getNrOfRows();
  // validate dts
  for (size_t i = 1; i < dts.size(); ++i) {
    INVALIDARG_EXP_IFN(dts[i - 1] < dts[i],
                       "Invalid dts: They need to be sequential.");
    INVALIDARG_EXP_IFN(dts[i] < nRows - 1,
                       "Invalid dts: got requested " + std::to_string(dts[i]) +
                         ", but only got " + std::to_string(nRows) + " rows.");
  }

  Eigen::ArrayXd colData =
    Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(
      this->getData()[column1].data(), this->getData()[column1].size()) -
    Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(
      this->getData()[column2].data(), this->getData()[column2].size());
  RUNTIME_EXP_IFN(colData.size() == nRows, "Invalid row sizes");

  std::vector<double> results;
  results.reserve(dts.size());
  for (size_t dt : dts) {
    results.push_back(
      (colData.segment(0, nRows - dt) * colData.segment(dt, nRows - dt))
        .mean());
  }

  return results;
}

}
