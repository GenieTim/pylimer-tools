#include "DataFileParser2.h"
#include "StringUtils.h"
#include <algorithm>
#include <filesystem>
#include <fstream> // std::ifstream
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {

  void DataFileParser2::read(const std::string filePath)
  {
    if (!std::filesystem::exists(filePath)) {
      throw std::invalid_argument("File to read (" + filePath +
                                  ") does not exist.");
    }

    std::string line;
    std::ifstream file;
    file.open(filePath);

    if (!file.is_open()) {
      throw std::invalid_argument("File to read (" + filePath +
                                  "): failed to open.");
    }

    // read everything until "Masses"
    while (getline(file, line)) {
      line = pylimer_tools::utils::trimLineOmitComment(line);
      // skip empty lines
      if (line.empty()) {
        continue;
      }
      // read up until the masses
      if (line.find("Masses") != std::string::npos) {
        break;
      }
      // read the nr of data points to read afterwards
      this->readNs(line);
    }

    // reserve space
    // for atom data
    this->atomIds.reserve(this->nAtoms);
    this->atomTypes.reserve(this->nAtoms);
    this->moleculeIds.reserve(this->nAtoms);
    this->atomX.reserve(this->nAtoms);
    this->atomY.reserve(this->nAtoms);
    this->atomZ.reserve(this->nAtoms);
    this->atomNx.reserve(this->nAtoms);
    this->atomNy.reserve(this->nAtoms);
    this->atomNz.reserve(this->nAtoms);
    // and bond data
    this->bondIds.reserve(this->nBonds);
    this->bondTypes.reserve(this->nBonds);
    this->bondFrom.reserve(this->nBonds);
    this->bondTo.reserve(this->nBonds);

    // skip empty lines plus the line with "Masses"
    while (getline(file, line)) {
      line = pylimer_tools::utils::trimLineOmitComment(line);

      // skip empty lines
      if (!line.empty()) {
        break;
      }
    }

    // Then, read masses, up until the next section ("atoms")
    do {
      line = pylimer_tools::utils::trimLineOmitComment(line);

      // skip empty lines
      if (line.empty()) {
        continue;
      }
      // read masses until e.g. atoms section
      if (line.find("Atoms") != std::string::npos ||
          line.find("Coeffs") != std::string::npos) {
        break;
      }
      // read the mass...
      this->readMass(line);
    } while (getline(file, line));

    this->skipLinesToContains(line, file, "Atoms");
    // skip this line too
    if (!getline(file, line)) {
      throw std::runtime_error(
        "Data file ended too early. Not able to read any atoms.");
    }
    // then, skip empty lines
    this->skipEmptyLines(line, file);

    // Then, read atoms, up until the next section ("bonds")
    for (int i = 0; i < this->nAtoms; ++i) {
      this->readAtom(line);

      if (!getline(file, line)) {
        throw std::runtime_error(
          "Data file ended too early. Not enough atoms read.");
      }
    }

    // Then, read bonds
    this->skipLinesToContains(line, file, "Bonds");
    // skip this line too
    if (!getline(file, line)) {
      throw std::runtime_error(
        "Data file ended too early. Not able to read any bonds.");
    }
    // then, skip empty lines
    this->skipEmptyLines(line, file);

    for (int i = 0; i < this->nBonds; i++) {
      this->readBond(line);

      if (!getline(file, line) && i + 1 < this->nBonds) {
        throw std::runtime_error(
          "Data file ended too early. Not enough bonds read.");
      }
    }

    // Then, read angles
    if (this->nAngles > 0) {
      this->skipLinesToContains(line, file, "Angles");
      // skip this line too
      if (!getline(file, line)) {
        throw std::runtime_error(
          "Data file ended too early. Not able to read any angles.");
      }
      // then, skip empty lines
      this->skipEmptyLines(line, file);

      for (int i = 0; i < this->nAngles; i++) {
        this->readAngle(line);

        if (!getline(file, line) && i + 1 < this->nAngles) {
          throw std::runtime_error(
            "Data file ended too early. Not enough angles read.");
        }
      }
    }

    // we ignore dihedrals etc. for now.
    file.close();
  }

  void DataFileParser2::skipLinesToContains(std::string& line,
                                           std::ifstream& file,
                                           std::string upTo)
  {
    do {
      if (contains(line, upTo)) {
        break;
      }
    } while (getline(file, line));
  }

  void DataFileParser2::skipEmptyLines(std::string& line, std::ifstream& file)
  {
    do {
      line = pylimer_tools::utils::trimLineOmitComment(line);

      // skip until empty lines
      if (!line.empty()) {
        break;
      }
    } while (getline(file, line));
  }

  void DataFileParser2::readNs(const std::string line)
  {
    if (contains(line, "atoms")) {
      this->nAtoms = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "bonds")) {
      this->nBonds = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "angles")) {
      this->nAngles = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "atom types")) {
      this->nAtomTypes = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "bond types")) {
      this->nBondTypes = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "angle types")) {
      this->nAngleTypes = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "xlo xhi")) {
      std::vector<double> parsedL = this->parseTypesInLine<double>(line, 2);
      this->xHi = parsedL[1];
      this->xLo = parsedL[0];
    } else if (contains(line, "ylo yhi")) {
      std::vector<double> parsedL = this->parseTypesInLine<double>(line, 2);
      this->yHi = parsedL[1];
      this->yLo = parsedL[0];
    } else if (contains(line, "zlo zhi")) {
      std::vector<double> parsedL = this->parseTypesInLine<double>(line, 2);
      this->zHi = parsedL[1];
      this->zLo = parsedL[0];
    }
  }

  void DataFileParser2::readMass(const std::string line)
  {
    int iteration = 0;
    int key = 0;
    pylimer_tools::utils::CsvTokenizer tokenizer(line);
    if (tokenizer.getLength() != 2) {
      throw std::runtime_error(
        "Incorrect nr of fields tokenized when reading masses");
    }

    key = tokenizer.get<int>(0);
    // for now, we just override duplicate keys
    this->masses[key] = tokenizer.get<double>(1);
  }

  void DataFileParser2::readAtom(std::string line)
  {
    size_t atomId,nx, ny, nz;
    int  atomType, moleculeId;
    double x, y, z;
    int resFound = sscanf(line.c_str(),
                          "%zd %d %d %le %le %le %zd %zd %zd",
                          &atomId,
                          &moleculeId,
                          &atomType,
                          &x,
                          &y,
                          &z,
                          &nx,
                          &ny,
                          &nz);
    // pylimer_tools::utils::CsvTokenizer tokenizer(line, 9);

    this->atomIds.push_back(atomId);
    this->moleculeIds.push_back(moleculeId);
    this->atomTypes.push_back(atomType);
    this->atomX.push_back(x);
    this->atomY.push_back(y);
    this->atomZ.push_back(z);

    if (resFound > 5) {
      this->atomNx.push_back(nx);
      this->atomNy.push_back(ny);
      this->atomNz.push_back(nz);
    }
    // this->atomIds.push_back(tokenizer.get<int>(0));
    // this->moleculeIds.push_back(tokenizer.get<int>(1));
    // this->atomTypes.push_back(tokenizer.get<int>(2));
    // this->atomX.push_back(tokenizer.get<double>(3));
    // this->atomY.push_back(tokenizer.get<double>(4));
    // this->atomZ.push_back(tokenizer.get<double>(5));
    // // TODO: be more flexible towards
    // if (tokenizer.getLength() > 5) {
    //   this->atomNx.push_back(tokenizer.get<int>(6));
    //   this->atomNy.push_back(tokenizer.get<int>(7));
    //   this->atomNz.push_back(tokenizer.get<int>(8));
    // }
  }

  void DataFileParser2::readBond(std::string line)
  {
    size_t bondId, bondType, bondFrom, bondTo;
    // pylimer_tools::utils::CsvTokenizer tokenizer(line, 4);
    sscanf(line.c_str(), "%zu %zu %zu %zu", &bondId, &bondType, &bondFrom, &bondTo);
    this->bondIds.push_back(bondId);
    this->bondTypes.push_back(bondType);
    this->bondFrom.push_back(bondFrom);
    this->bondTo.push_back(bondTo);
    // this->bondIds.push_back(tokenizer.get<long int>(0));
    // this->bondTypes.push_back(tokenizer.get<int>(1));
    // this->bondFrom.push_back(tokenizer.get<long int>(2));
    // this->bondTo.push_back(tokenizer.get<long int>(3));
  }

  void DataFileParser2::readAngle(std::string line)
  {
    // pylimer_tools::utils::CsvTokenizer tokenizer(line, 5);
    size_t angleId, angleType, angleFrom, angleVia, angleTo;
    sscanf(line.c_str(),
           "%zu %zu %zu %zu %zu",
           &angleId,
           &angleType,
           &angleFrom,
           &angleVia,
           &angleTo);

    this->angleIds.push_back(angleId);
    this->angleTypes.push_back(angleType);
    this->angleFrom.push_back(angleFrom);
    this->angleVia.push_back(angleVia);
    this->angleTo.push_back(angleTo);
    // this->angleIds.push_back(tokenizer.get<long int>(0));
    // this->angleTypes.push_back(tokenizer.get<int>(1));
    // this->angleFrom.push_back(tokenizer.get<long int>(2));
    // this->angleVia.push_back(tokenizer.get<long int>(3));
    // this->angleTo.push_back(tokenizer.get<long int>(4));
  }
} // namespace utils
} // namespace pylimer_tools
