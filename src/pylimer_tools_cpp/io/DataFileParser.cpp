#include "DataFileParser.h"
#include "../utils/StringUtils.h"
#include <algorithm>
#include <filesystem>
#include <fstream> // std::ifstream
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {

  void DataFileParser::read(const std::string filePath,
                            const AtomStyle atomStyle,
                            const AtomStyle atomStyle2,
                            const AtomStyle atomStyle3)
  {
    if (!std::filesystem::exists(filePath)) {
      throw std::invalid_argument("Data file to read (" + filePath +
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
      switch (atomStyle) {
        case AtomStyle::ANGLE:
        case AtomStyle::BOND:
        case AtomStyle::MOLECULAR:
          this->readAtom(line);
          break;
        case AtomStyle::HYBRID:
          this->readAtomHybrid(line, atomStyle2, atomStyle3);
          break;
        case AtomStyle::CHARGE:
          this->readAtomCharge(line);
          break;
        case AtomStyle::FULL:
          this->readAtomFull(line);
        default:
          throw std::invalid_argument("This atom style is not supported yet.");
          break;
      }

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

    // then, read dihedral angles
    if (this->nDihedralAngles > 0) {
      this->skipLinesToContains(line, file, "Dihedrals");
      // skip this line too
      if (!getline(file, line)) {
        throw std::runtime_error(
          "Data file ended too early. Not able to read any dihedral angles.");
      }
      // then, skip empty lines
      this->skipEmptyLines(line, file);

      for (int i = 0; i < this->nDihedralAngles; i++) {
        this->readDihedralAngle(line);

        if (!getline(file, line) && i + 1 < this->nDihedralAngles) {
          throw std::runtime_error(
            "Data file ended too early. Not enough dihedral angles read.");
        }
      }
    }

    // we ignore dihedrals etc. for now.
    file.close();
  }

  void DataFileParser::skipLinesToContains(std::string& line,
                                           std::ifstream& file,
                                           std::string upTo)
  {
    do {
      if (contains(line, upTo)) {
        break;
      }
    } while (getline(file, line));
  }

  void DataFileParser::skipEmptyLines(std::string& line, std::ifstream& file)
  {
    do {
      line = pylimer_tools::utils::trimLineOmitComment(line);

      // skip until empty lines
      if (!line.empty()) {
        break;
      }
    } while (getline(file, line));
  }

  void DataFileParser::readNs(const std::string line)
  {
    if (contains(line, "atoms")) {
      this->nAtoms = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "bonds")) {
      this->nBonds = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "angles")) {
      this->nAngles = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "dihedrals")) {
      this->nDihedralAngles = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "atom types")) {
      this->nAtomTypes = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "bond types")) {
      this->nBondTypes = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "angle types")) {
      this->nAngleTypes = (this->parseTypesInLine<int>(line, 1))[0];
    } else if (contains(line, "dihedral types")) {
      this->nDihedralAngleTypes = (this->parseTypesInLine<int>(line, 1))[0];
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

  void DataFileParser::readMass(const std::string line)
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

  void DataFileParser::readAtomFull(std::string line)
  {
    size_t atomId, nx, ny, nz;
    int atomType, moleculeId;
    double charge;
    double x, y, z;
    int resFound = sscanf(line.c_str(),
                          "%zd %d %d %le %le %le %le %zd %zd %zd",
                          &atomId,
                          &moleculeId,
                          &atomType,
                          &charge,
                          &x,
                          &y,
                          &z,
                          &nx,
                          &ny,
                          &nz);

    this->atomIds.push_back(atomId);
    this->moleculeIds.push_back(moleculeId);
    this->atomTypes.push_back(atomType);
    this->atomX.push_back(x);
    this->atomY.push_back(y);
    this->atomZ.push_back(z);
    this->additionalAtomData["charge"].push_back(charge);

    if (resFound > 6) {
      this->atomNx.push_back(nx);
      this->atomNy.push_back(ny);
      this->atomNz.push_back(nz);
    }
  }

  void DataFileParser::readAtomCharge(std::string line)
  {
    size_t atomId, nx, ny, nz;
    int atomType, moleculeId;
    double charge;
    double x, y, z;
    int resFound = sscanf(line.c_str(),
                          "%zd %d %le %le %le %le %zd %zd %zd",
                          &atomId,
                          &atomType,
                          &charge,
                          &x,
                          &y,
                          &z,
                          &nx,
                          &ny,
                          &nz);

    this->atomIds.push_back(atomId);
    this->moleculeIds.push_back(moleculeId);
    this->atomTypes.push_back(atomType);
    this->atomX.push_back(x);
    this->atomY.push_back(y);
    this->atomZ.push_back(z);
    this->additionalAtomData["charge"].push_back(charge);

    if (resFound > 6) {
      this->atomNx.push_back(nx);
      this->atomNy.push_back(ny);
      this->atomNz.push_back(nz);
    }
  }

  void DataFileParser::readAtom(std::string line)
  {
    size_t atomId, nx, ny, nz;
    int atomType, moleculeId;
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

    this->atomIds.push_back(atomId);
    this->moleculeIds.push_back(moleculeId);
    this->atomTypes.push_back(atomType);
    this->atomX.push_back(x);
    this->atomY.push_back(y);
    this->atomZ.push_back(z);

    if (resFound > 6) {
      this->atomNx.push_back(nx);
      this->atomNy.push_back(ny);
      this->atomNz.push_back(nz);
    }
  }

  void DataFileParser::readAtomHybrid(std::string line,
                                      AtomStyle style1,
                                      AtomStyle style2)
  {
    if (style1 != AtomStyle::BOND && style2 != AtomStyle::EDPD) {
      throw std::runtime_error(
        "This combination is not implemented for hybrid atom style");
    }
    size_t atomId, nx, ny, nz;
    int atomType, moleculeId;
    double x, y, z;
    double edpdTemp, edpd;
    int resFound = sscanf(line.c_str(),
                          "%zd %d %le %le %le %d %le %le %zd %zd %zd",
                          &atomId,
                          &atomType,
                          &x,
                          &y,
                          &z,
                          &moleculeId,
                          &edpdTemp,
                          &edpd,
                          &nx,
                          &ny,
                          &nz);

    this->atomIds.push_back(atomId);
    this->moleculeIds.push_back(moleculeId);
    this->atomTypes.push_back(atomType);
    this->atomX.push_back(x);
    this->atomY.push_back(y);
    this->atomZ.push_back(z);
    this->additionalAtomData["edpd_temp"].push_back(edpdTemp);
    this->additionalAtomData["edpd"].push_back(edpd);

    if (resFound > 8) {
      this->atomNx.push_back(nx);
      this->atomNy.push_back(ny);
      this->atomNz.push_back(nz);
    }
  }

  void DataFileParser::readBond(std::string line)
  {
    size_t bondId, bondType, bondFrom, bondTo;
    sscanf(
      line.c_str(), "%zu %zu %zu %zu", &bondId, &bondType, &bondFrom, &bondTo);
    this->bondIds.push_back(bondId);
    this->bondTypes.push_back(bondType);
    this->bondFrom.push_back(bondFrom);
    this->bondTo.push_back(bondTo);
  }

  void DataFileParser::readAngle(std::string line)
  {
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
  }

  void DataFileParser::readDihedralAngle(std::string line)
  {
    size_t angleId, angleType, angleFrom, angleVia1, angleVia2, angleTo;
    sscanf(line.c_str(),
           "%zu %zu %zu %zu %zu %zu",
           &angleId,
           &angleType,
           &angleFrom,
           &angleVia1,
           &angleVia2,
           &angleTo);

    this->dihedralAngleIds.push_back(angleId);
    this->dihedralAngleTypes.push_back(angleType);
    this->dihedralAngleFrom.push_back(angleFrom);
    this->dihedralAngleVia1.push_back(angleVia1);
    this->dihedralAngleVia2.push_back(angleVia2);
    this->dihedralAngleTo.push_back(angleTo);
  }
} // namespace utils
} // namespace pylimer_tools
