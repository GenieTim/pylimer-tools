#ifndef DATA_FILE_PARSER_H
#define DATA_FILE_PARSER_H

#include "../utils/LammpsAtomStyle.h"
#include "../utils/StringUtils.h"
#include <filesystem>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace pylimer_tools::utils {

class DataFileParser
{

public:
  void read(const std::string& filePath,
            AtomStyle atomStyle = AtomStyle::ANGLE,
            const AtomStyle atomStyle2 = AtomStyle::NONE,
            const AtomStyle atomStyle3 = AtomStyle::NONE);

  // access atom data
  int getNrOfAtoms() const { return this->nAtoms; }
  int getNrOfAtomTypes() const { return this->nAtomTypes; }
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
  int getNrOfBonds() const { return this->nBonds; }
  int getNrOfBondTypes() const { return this->nBondTypes; }
  std::vector<int> getBondTypes() { return this->bondTypes; }
  std::vector<long int> getBondFrom() { return this->bondFrom; }
  std::vector<long int> getBondTo() { return this->bondTo; }

  // access angle data
  int getNrOfAngles() const { return this->nAngles; }
  int getNrOfAngleTypes() const { return this->nAngleTypes; }
  std::vector<int> getAngleTypes() { return this->angleTypes; }
  std::vector<long int> getAngleFrom() { return this->angleFrom; }
  std::vector<long int> getAngleVia() { return this->angleVia; }
  std::vector<long int> getAngleTo() { return this->angleTo; }

  // access dihedral angle data
  int getNrOfDihedralAngles() const { return this->nDihedralAngles; }
  int getNrOfDihedralAngleTypes() const { return this->nDihedralAngleTypes; }
  std::vector<int> getDihedralAngleTypes() { return this->dihedralAngleTypes; }
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
  double getLowX() const { return this->xLo; }
  double getHighX() const { return this->xHi; }
  double getLx() const { return this->xHi - this->xLo; }
  double getLowY() const { return this->yLo; }
  double getHighY() const { return this->yHi; }
  double getLy() const { return this->yHi - this->yLo; }
  double getLowZ() const { return this->zLo; }
  double getHighZ() const { return this->zHi; }
  double getLz() const { return this->zHi - this->zLo; }

private:
  // Atom style parsing descriptors
  struct AtomFieldDescriptor
  {
    enum FieldType
    {
      // Basic atom identifiers
      ATOM_ID,
      MOLECULE_ID,
      ATOM_TYPE,

      // Position coordinates
      X,
      Y,
      Z,

      // Image flags
      NX,
      NY,
      NZ,

      // Physical properties
      ATOM_CHARGE,
      MASS,
      DENSITY,
      VOLUME,
      DIAMETER,

      // Dipole/magnetic properties
      MUX,
      MUY,
      MUZ,

      // Temperature and thermal properties
      THETA,
      RHO,
      EDPD_TEMP,
      ATOM_EDPD,
      EDPD_CV,
      ESPH,
      CV,
      ENERGY,

      // Flags and status
      BODYFLAG,
      ELLIPSOIDFLAG,
      LINEFLAG,
      TRIANGLEFLAG,
      STATUS,

      // Radii and geometric properties
      KRADIUS,
      CRADIUS,
      ERADIUS,

      // Initial positions
      X0,
      Y0,
      Z0,

      // Spin properties
      SPX,
      SPY,
      SPZ,
      SP,
      ESPIN,

      // Surface and dielectric properties
      AREA,
      ED,
      EM,
      EPSILON,
      CURVATURE,

      // Template properties
      TEMPLATE_INDEX,
      TEMPLATE_ATOM,

      // Wave packet properties
      ETAG,
      CS_RE,
      CS_IM,

      // Species concentrations (for tdpd)
      CC1,
      CC2 // Can be extended for more species as needed
    };
    std::vector<FieldType> fields;
    std::string format;
  };

  void readNs(const std::string& line);
  void readMass(const std::string& line);
  // different atom styles
  void readAtom(const std::string& line);
  void readAtomFull(const std::string& line);
  void readAtomCharge(const std::string& line);
  void readAtomDipole(const std::string& line);
  void readAtomDpd(const std::string& line);
  void readAtomMdpd(const std::string& line);
  void readAtomSphere(const std::string& line);
  void readAtomHybrid(const std::string& line,
                      AtomStyle style1,
                      AtomStyle style2);
  void readAtomAtomic(const std::string& line);
  void readAtomBody(const std::string& line);
  void readAtomBpmSphere(const std::string& line);
  void readAtomDielectric(const std::string& line);
  void readAtomEdpd(const std::string& line);
  void readAtomElectron(const std::string& line);
  void readAtomEllipsoid(const std::string& line);
  void readAtomLine(const std::string& line);
  void readAtomPeri(const std::string& line);
  void readAtomRheo(const std::string& line);
  void readAtomRheoThermal(const std::string& line);
  void readAtomSmd(const std::string& line);
  void readAtomSph(const std::string& line);
  void readAtomSpin(const std::string& line);
  void readAtomTemplate(const std::string& line);
  void readAtomTri(const std::string& line);
  void readAtomWavepacket(const std::string& line);

  template<typename... Args>
  void readAtomGeneric(const std::string& line,
                       const AtomFieldDescriptor& descriptor,
                       Args&... args);

  AtomFieldDescriptor getAtomStyleDescriptor(AtomStyle style) const;

  // bonds, angles, etc.
  void readBonds(std::ifstream& file, std::string& line);
  void readBond(const std::string& line);
  void readAngles(std::ifstream& file, std::string& line);
  void readAngle(const std::string& line);
  void readDihedralAngles(std::ifstream& file, std::string& line);
  void readDihedralAngle(const std::string& line);
  // additional atom data
  void readVelocities(std::ifstream& file, std::string& line);

  // utilities
  static void skipEmptyLines(std::string& line, std::ifstream& file);
  static void skipLinesToContains(std::string& line,
                                  std::ifstream& file,
                                  const std::string& upTo);
  static void skipLinesToContains(std::string& line,
                                  std::ifstream& file,
                                  const std::vector<std::string>& upToEitherOr);

  template<typename OUT>
  inline std::vector<OUT> parseTypesInLine(const std::string& line,
                                           const int nToRead)
  {
    std::vector<OUT> resultnumbers;
    pylimer_tools::utils::CsvTokenizer tokenizer(line,
                                                 static_cast<size_t>(nToRead));
    resultnumbers.reserve(tokenizer.getLength());
    for (size_t i = 0; i < tokenizer.getLength(); ++i) {
      resultnumbers.push_back(tokenizer.get<OUT>(i));
    }
    return resultnumbers;
  }

  //// data
  // nr of data points to read
  int nAtoms = 0; // number of atoms
  int nBonds = 0;
  int nAngles = 0;
  int nAtomTypes = 0;
  int nBondTypes = 0;
  int nAngleTypes = 0;
  int nDihedralAngles = 0;
  int nDihedralAngleTypes = 0;

  // box sizes
  double xLo = 0;
  double xHi = 0;
  double yLo = 0;
  double yHi = 0;
  double zLo = 0;
  double zHi = 0;

  // actual dimensional values
  std::map<int, double> masses;
  std::unordered_map<std::string, std::vector<double>> additionalAtomData;
  std::vector<long int> atomIds = {};
  std::vector<int> moleculeIds = {};
  std::vector<int> atomTypes = {};
  std::vector<double> atomX = {};
  std::vector<double> atomY = {};
  std::vector<double> atomZ = {};
  std::vector<int> atomNx = {};
  std::vector<int> atomNy = {};
  std::vector<int> atomNz = {};

  // bonds
  std::vector<long int> bondIds = {};
  std::vector<int> bondTypes = {};
  std::vector<long int> bondFrom = {};
  std::vector<long int> bondTo = {};

  // angles
  std::vector<long int> angleIds = {};
  std::vector<int> angleTypes = {};
  std::vector<long int> angleFrom = {};
  std::vector<long int> angleVia = {};
  std::vector<long int> angleTo = {};

  // dihedrals
  std::vector<long int> dihedralAngleIds = {};
  std::vector<int> dihedralAngleTypes = {};
  std::vector<long int> dihedralAngleFrom = {};
  std::vector<long int> dihedralAngleVia1 = {};
  std::vector<long int> dihedralAngleVia2 = {};
  std::vector<long int> dihedralAngleTo = {};
};
}

#endif
