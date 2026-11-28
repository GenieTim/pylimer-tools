#ifndef DATA_FILE_PARSER_H
#define DATA_FILE_PARSER_H

#include "../utils/LammpsAtomStyle.h"
#include "../utils/StringUtils.h"
#include <filesystem>
#include <fstream>
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
  const std::vector<long int>& getAtomIds() const { return this->atomIds; }
  const std::vector<int>& getMoleculeIds() const { return this->moleculeIds; }
  const std::vector<int>& getAtomTypes() const { return this->atomTypes; }
  const std::vector<double>& getAtomX() const { return this->atomX; }
  const std::vector<double>& getAtomY() const { return this->atomY; }
  const std::vector<double>& getAtomZ() const { return this->atomZ; }
  const std::vector<int>& getAtomNx() const { return this->atomNx; }
  const std::vector<int>& getAtomNy() const { return this->atomNy; }
  const std::vector<int>& getAtomNz() const { return this->atomNz; }
  const std::map<int, double>& getMasses() const { return this->masses; }

  // access bond data
  int getNrOfBonds() const { return this->nBonds; }
  int getNrOfBondTypes() const { return this->nBondTypes; }
  const std::vector<long int>& getBondIds() const { return this->bondIds; }
  const std::vector<int>& getBondTypes() const { return this->bondTypes; }
  const std::vector<long int>& getBondFrom() const { return this->bondFrom; }
  const std::vector<long int>& getBondTo() const { return this->bondTo; }

  // access angle data
  int getNrOfAngles() const { return this->nAngles; }
  int getNrOfAngleTypes() const { return this->nAngleTypes; }
  const std::vector<long int>& getAngleIds() const { return this->angleIds; }
  const std::vector<int>& getAngleTypes() const { return this->angleTypes; }
  const std::vector<long int>& getAngleFrom() const { return this->angleFrom; }
  const std::vector<long int>& getAngleVia() const { return this->angleVia; }
  const std::vector<long int>& getAngleTo() const { return this->angleTo; }

  // access dihedral angle data
  int getNrOfDihedralAngles() const { return this->nDihedralAngles; }
  int getNrOfDihedralAngleTypes() const { return this->nDihedralAngleTypes; }
  const std::vector<long int>& getDihedralAngleIds() const
  {
    return this->dihedralAngleIds;
  }
  const std::vector<int>& getDihedralAngleTypes() const
  {
    return this->dihedralAngleTypes;
  }
  const std::vector<long int>& getDihedralAngleFrom() const
  {
    return this->dihedralAngleFrom;
  }
  const std::vector<long int>& getDihedralAngleVia1() const
  {
    return this->dihedralAngleVia1;
  }
  const std::vector<long int>& getDihedralAngleVia2() const
  {
    return this->dihedralAngleVia2;
  }
  const std::vector<long int>& getDihedralAngleTo() const
  {
    return this->dihedralAngleTo;
  }
  const std::unordered_map<std::string, std::vector<double>>&
  getAdditionalAtomData() const
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
    bool hasMoleculeId = false;

    AtomFieldDescriptor(std::initializer_list<FieldType> fieldList)
      : fields(fieldList)
      , hasMoleculeId(std::find(fields.begin(),
                                fields.end(),
                                FieldType::MOLECULE_ID) != fields.end())
    {
    }
  };

  void reset();
  void readNs(std::string_view line);
  void readMass(std::string_view line);
  // atom styles
  void readAtomsSection(std::ifstream& file,
                        AtomStyle atomStyle,
                        AtomStyle atomStyle2,
                        AtomStyle atomStyle3);
  void readAtomHybrid(std::string_view line,
                      AtomStyle style1,
                      AtomStyle style2);

  static const AtomFieldDescriptor& getAtomStyleDescriptor(AtomStyle style);

  // bonds, angles, etc.
  template<typename ParseFn>
  void readSectionLines(std::ifstream& file,
                        int count,
                        const char* sectionName,
                        ParseFn&& parseLine);
  void readBonds(std::ifstream& file);
  void readBond(std::string_view line);
  void readAngles(std::ifstream& file);
  void readAngle(std::string_view line);
  void readDihedralAngles(std::ifstream& file);
  void readDihedralAngle(std::string_view line);
  // additional atom data
  void readVelocities(std::ifstream& file,
                      std::vector<long int>& pendingVelAtomIds,
                      std::vector<double>& pendingVx,
                      std::vector<double>& pendingVy,
                      std::vector<double>& pendingVz);
  void finalizePendingVelocities(const std::vector<long int>& pendingVelAtomIds,
                                 const std::vector<double>& pendingVx,
                                 const std::vector<double>& pendingVy,
                                 const std::vector<double>& pendingVz);
  static bool nextNonEmptyDataLine(std::ifstream& file, std::string& lineOut);

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
