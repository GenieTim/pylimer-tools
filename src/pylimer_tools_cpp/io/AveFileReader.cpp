
#include "./AveFileReader.h"

#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <cstring>
#include <filesystem>
#include <fstream> // std::ifstream

namespace pylimer_tools {
namespace utils {

  /**
   * @brief Count the number of lines with data
   *
   * @return int
   */
  int AveFileReader::getNrOfRows()
  {
    if (this->numRows > 0) {
      return this->numRows;
    }
    std::ifstream in_stream(this->filePath);
    if (!in_stream) {
      throw std::runtime_error("Failed open file " + this->filePath);
    }

    auto count = std::count_if(std::istreambuf_iterator<char>{ in_stream },
                               {},
                               [](char c) { return c == '\n'; });
    this->numRows = count - this->getNrOfHeaderRows();
    return this->numRows;
  }

  /**
   * @brief Count the number of (comment/header) rows at the top of the file
   *
   * @return int
   */
  int AveFileReader::getNrOfHeaderRows()
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
  std::vector<std::string> AveFileReader::getColumnNames()
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
      if (line.at(0) == '#') {
        lastLine = line;
      } else {
        break;
      }
    }

    std::string headerLine =
      pylimer_tools::utils::trimLineOmitComment(lastLine);
    this->columnNames = pylimer_tools::utils::split(headerLine, ' ');
    return this->columnNames;
  }

  /**
   * @brief Count the nr of columns in the file
   *
   * @return int
   */
  int AveFileReader::getNrOfColumns()
  {
    return this->getColumnNames().size();
  }

  /**
   * @brief Actually read the file's data
   *
   * @return std::vector<std::vector<double>>
   */
  std::vector<std::vector<double>> AveFileReader::getData()
  {
    if (this->data.size() > 0) {
      return this->data;
    }

    int numRows = this->getNrOfRows();
    int numCols = this->getNrOfColumns();
    std::vector<std::vector<double>> results =
      pylimer_tools::utils::initializeWithValue(numCols, std::vector<double>());
    for (size_t i = 0; i < numCols; ++i) {
      results[i].reserve(numCols);
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

    for (size_t i = 0; i < numRows; ++i) {
      RUNTIME_EXP_IFN(std::getline(file, line),
                      "File ended before all rows could be read (reached row " +
                        std::to_string(i) + " of " + std::to_string(numRows) +
                        ")");
      std::stringstream ss(line);
      for (size_t col = 0; col < numCols; ++col) {
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

  std::vector<double> AveFileReader::autocorrelateColumn(
    int column,
    const std::vector<size_t>& dts)
  {
    INVALIDARG_EXP_IFN(column < this->getNrOfColumns(), "Invalid column");

    int numRows = this->getNrOfRows();
    // validate dts
    for (size_t i = 1; i < dts.size(); ++i) {
      INVALIDARG_EXP_IFN(dts[i - 1] < dts[i], "Invalid dts");
      INVALIDARG_EXP_IFN(dts[i] < numRows - 1, "Invalid dts");
    }

    Eigen::ArrayXd colData = Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(
      this->getData()[column].data(), this->getData()[column].size());
    RUNTIME_EXP_IFN(colData.size() == numRows, "Invalid row sizes");

    std::vector<double> results;
    results.reserve(dts.size());
    for (size_t dt : dts) {
      results.push_back(
        (colData.segment(0, numRows - dt) * colData.segment(dt, numRows))
          .mean());
    }

    return results;
  }

  std::vector<double> AveFileReader::autocorrelateColumnDifference(
    int column1,
    int column2,
    const std::vector<size_t>& dts)
  {
    INVALIDARG_EXP_IFN(column1 < this->getNrOfColumns(), "Invalid column");
    INVALIDARG_EXP_IFN(column2 < this->getNrOfColumns(), "Invalid column");

    int numRows = this->getNrOfRows();
    // validate dts
    for (size_t i = 1; i < dts.size(); ++i) {
      INVALIDARG_EXP_IFN(dts[i - 1] < dts[i], "Invalid dts");
      INVALIDARG_EXP_IFN(dts[i] < numRows - 1, "Invalid dts");
    }

    Eigen::ArrayXd colData =
      Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(
        this->getData()[column1].data(), this->getData()[column1].size()) -
      Eigen::Map<Eigen::ArrayXd, Eigen::Unaligned>(
        this->getData()[column2].data(), this->getData()[column2].size());
    RUNTIME_EXP_IFN(colData.size() == numRows, "Invalid row sizes");

    std::vector<double> results;
    results.reserve(dts.size());
    for (size_t dt : dts) {
      results.push_back(
        (colData.segment(0, numRows - dt) * colData.segment(dt, numRows))
          .mean());
    }

    return results;
  }

}
}
