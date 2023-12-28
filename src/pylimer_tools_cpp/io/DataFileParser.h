#ifndef DATA_FILE_PARSER_H
#define DATA_FILE_PARSER_H

#include "../utils/StringUtils.h"
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {

  enum AtomStyle
  {
    NONE,
    ANGLE,
    ATOMIC,
    BODY,
    BOND,
    CHARGE,
    DIELECTRIC,
    DIPOLE,
    DPD,
    EDPD,
    ELECTRON,
    ELLIPSOID,
    FULL,
    LINE,
    MDPD,
    MOLECULAR,
    PERI,
    SMD,
    SPH,
    SPHERE,
    BPM_SPHERE,
    SPIN,
    TDPD,
    TEMPLATE,
    TRI,
    WAVEPACKET,
    HYBRID
  };

  class DataFileParser
  {
  public:
    void read(const std::string filePath,
              const AtomStyle atomStyle = AtomStyle::ANGLE,
              const AtomStyle atomStyle2 = AtomStyle::NONE,
              const AtomStyle atomStyle3 = AtomStyle::NONE);

    // access atom data
    int getNrOfAtoms() { return this->nAtoms; }
    int getNrOfAtomTypes() { return this->nAtomTypes; }
    std::vector<long int> getAtomIds() { return this->atomIds; }
    std::vector<int> getMoleculeIds() { return this->moleculeIds; }
    std::vector<int> getAtomTypes() { return this->atomTypes; }
    std::vector<double> getAtomX() { return this->atomX; }
    std::vector<double> getAtomY() { return this->atomY; }
    std::vector<double> getAtomZ() { return this->atomZ; }
    std::vector<int> getAtomNx() { return this->atomNx; }
    std::vector<int> getAtomNy() { return this->atomNy; }
    std::vector<int> getAtomNz() { return this->atomNz; }
    std::map<int, double> getMasses() { return this->masses; }

    // access bond data
    int getNrOfBonds() { return this->nBonds; }
    int getNrOfBondTypes() { return this->nBondTypes; }
    std::vector<int> getBondTypes() { return this->bondTypes; }
    std::vector<long int> getBondFrom() { return this->bondFrom; }
    std::vector<long int> getBondTo() { return this->bondTo; }

    // access angle data
    int getNrOfAngles() { return this->nAngles; }
    int getNrOfAngleTypes() { return this->nAngleTypes; }
    std::vector<int> getAngleTypes() { return this->angleTypes; }
    std::vector<long int> getAngleFrom() { return this->angleFrom; }
    std::vector<long int> getAngleVia() { return this->angleVia; }
    std::vector<long int> getAngleTo() { return this->angleTo; }

    // access dihedral angle data
    int getNrOfDihedralAngles() { return this->nDihedralAngles; }
    int getNrOfDihedralAngleTypes() { return this->nDihedralAngleTypes; }
    std::vector<int> getDihedralAngleTypes()
    {
      return this->dihedralAngleTypes;
    }
    std::vector<long int> getDihedralAngleFrom()
    {
      return this->dihedralAngleFrom;
    }
    std::vector<long int> getDihedralAngleVia1()
    {
      return this->dihedralAngleVia1;
    }
    std::vector<long int> getDihedralAngleVia2()
    {
      return this->dihedralAngleVia2;
    }
    std::vector<long int> getDihedralAngleTo() { return this->dihedralAngleTo; }
    std::unordered_map<std::string, std::vector<double>> getAdditionalAtomData()
    {
      return this->additionalAtomData;
    }

    // get box info
    double getLowX() { return this->xLo; }
    double getHighX() { return this->xHi; }
    double getLx() { return this->xHi - this->xLo; }
    double getLowY() { return this->yLo; }
    double getHighY() { return this->yHi; }
    double getLy() { return this->yHi - this->yLo; }
    double getLowZ() { return this->zLo; }
    double getHighZ() { return this->zHi; }
    double getLz() { return this->zHi - this->zLo; }

  private:
    void readNs(const std::string &line);
    void readMass(const std::string &line);
    // different atom styles
    void readAtom(const std::string &line);
    void readAtomFull(const std::string &line);
    void readAtomCharge(const std::string& line);
    void readAtomHybrid(const std::string &line, AtomStyle style1, AtomStyle style2);
    // bonds, angles, etc.
    void readBond(const std::string &line);
    void readAngle(const std::string &line);
    void readDihedralAngle(const std::string &line);

    // utilities
    static void skipEmptyLines(std::string& line, std::ifstream& file);
    static void skipLinesToContains(std::string& line,
                                    std::ifstream& file,
                                    std::string upTo);

    template<typename OUT>
    inline std::vector<OUT> parseTypesInLine(const std::string &line,
                                             int nToRead)
    {
      std::vector<OUT> resultnumbers;
      pylimer_tools::utils::CsvTokenizer tokenizer(line, nToRead);
      resultnumbers.reserve(tokenizer.getLength());
      for (size_t i = 0; i < tokenizer.getLength(); ++i) {
        resultnumbers.push_back(tokenizer.get<OUT>(i));
      }
      return resultnumbers;
    }

    //// data
    // nr of data points to read
    int nAtoms; // number of atoms
    int nBonds;
    int nAngles;
    int nAtomTypes;
    int nBondTypes;
    int nAngleTypes;
    int nDihedralAngles;
    int nDihedralAngleTypes;

    // box sizes
    double xLo;
    double xHi;
    double yLo;
    double yHi;
    double zLo;
    double zHi;

    // actual dimensional values
    std::map<int, double> masses;
    std::unordered_map<std::string, std::vector<double>> additionalAtomData;
    std::vector<long int> atomIds;
    std::vector<int> moleculeIds;
    std::vector<int> atomTypes;
    std::vector<double> atomX;
    std::vector<double> atomY;
    std::vector<double> atomZ;
    std::vector<int> atomNx;
    std::vector<int> atomNy;
    std::vector<int> atomNz;

    // bonds
    std::vector<long int> bondIds;
    std::vector<int> bondTypes;
    std::vector<long int> bondFrom;
    std::vector<long int> bondTo;

    // angles
    std::vector<long int> angleIds;
    std::vector<int> angleTypes;
    std::vector<long int> angleFrom;
    std::vector<long int> angleVia;
    std::vector<long int> angleTo;

    // dihedrals
    std::vector<long int> dihedralAngleIds;
    std::vector<int> dihedralAngleTypes;
    std::vector<long int> dihedralAngleFrom;
    std::vector<long int> dihedralAngleVia1;
    std::vector<long int> dihedralAngleVia2;
    std::vector<long int> dihedralAngleTo;
  };
} // namespace utils
} // namespace pylimer_tools

#endif
