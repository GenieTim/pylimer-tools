#include "DataFileParser.h"
#include "../utils/StringUtils.h"
#include "../utils/utilityMacros.h"
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream> // std::ifstream
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace pylimer_tools::utils {

namespace {
  bool isDataFileSectionHeader(const std::string_view trimmedNoComment)
  {
    return trimmedNoComment == "Masses" || trimmedNoComment == "Atoms" ||
           trimmedNoComment == "Bonds" || trimmedNoComment == "Angles" ||
           trimmedNoComment == "Dihedrals" || trimmedNoComment == "Impropers" ||
           trimmedNoComment == "Velocities" ||
           trimmedNoComment == "Ellipsoids" || trimmedNoComment == "Lines" ||
           trimmedNoComment == "Triangles" || trimmedNoComment == "Bodies" ||
           trimmedNoComment.ends_with("Coeffs") ||
           trimmedNoComment.ends_with("Labels");
  }

  bool hasRemainingTokens(const char* cursor)
  {
    while (*cursor == ' ' || *cursor == '\t' || *cursor == '\r' ||
           *cursor == '\n' || *cursor == '\f' || *cursor == '\v') {
      ++cursor;
    }
    return *cursor != '\0' && *cursor != '#';
  }
} // namespace

void
DataFileParser::reset()
{
  this->nAtoms = 0;
  this->nBonds = 0;
  this->nAngles = 0;
  this->nAtomTypes = 0;
  this->nBondTypes = 0;
  this->nAngleTypes = 0;
  this->nDihedralAngles = 0;
  this->nDihedralAngleTypes = 0;

  this->xLo = 0.0;
  this->xHi = 0.0;
  this->yLo = 0.0;
  this->yHi = 0.0;
  this->zLo = 0.0;
  this->zHi = 0.0;

  this->masses.clear();
  this->additionalAtomData.clear();
  this->atomIds.clear();
  this->moleculeIds.clear();
  this->atomTypes.clear();
  this->atomX.clear();
  this->atomY.clear();
  this->atomZ.clear();
  this->atomNx.clear();
  this->atomNy.clear();
  this->atomNz.clear();

  this->bondIds.clear();
  this->bondTypes.clear();
  this->bondFrom.clear();
  this->bondTo.clear();

  this->angleIds.clear();
  this->angleTypes.clear();
  this->angleFrom.clear();
  this->angleVia.clear();
  this->angleTo.clear();

  this->dihedralAngleIds.clear();
  this->dihedralAngleTypes.clear();
  this->dihedralAngleFrom.clear();
  this->dihedralAngleVia1.clear();
  this->dihedralAngleVia2.clear();
  this->dihedralAngleTo.clear();
}

bool
DataFileParser::nextNonEmptyDataLine(std::ifstream& file, std::string& lineOut)
{
  while (std::getline(file, lineOut)) {
    const size_t commentPos = lineOut.find('#');
    if (commentPos != std::string::npos) {
      lineOut.erase(commentPos);
    }
    const size_t endPos =
      lineOut.find_last_not_of(pylimer_tools::utils::WHITESPACE_VIEW);
    if (endPos == std::string::npos) {
      continue;
    }
    lineOut.erase(endPos + 1);
    return true;
  }
  return false;
}

void
DataFileParser::read(const std::string& filePath,
                     const AtomStyle atomStyle,
                     const AtomStyle atomStyle2,
                     const AtomStyle atomStyle3)
{
  if (!std::filesystem::exists(filePath)) {
    throw std::invalid_argument("Data file to read (" + filePath +
                                ") does not exist.");
  }

  std::ifstream file(filePath);
  if (!file.is_open()) {
    throw std::invalid_argument("File to read (" + filePath +
                                "): failed to open.");
  }

  this->reset();

  std::string rawLine;
  // First line of a LAMMPS data file is always a title/comment line
  if (!std::getline(file, rawLine)) {
    return;
  }

  bool hasPendingSection = false;
  while (std::getline(file, rawLine)) {
    const std::string_view cleaned = pylimer_tools::utils::trimWhitespaceView(
      pylimer_tools::utils::trimLineOmitCommentView(rawLine));
    if (cleaned.empty()) {
      continue;
    }
    if (isDataFileSectionHeader(cleaned)) {
      hasPendingSection = true;
      break;
    }
    this->readNs(cleaned);
  }

  RUNTIME_EXP_IFN(this->nAtoms >= 0 && this->nBonds >= 0 &&
                    this->nAngles >= 0 && this->nDihedralAngles >= 0,
                  "Negative entity count encountered in data file header.");

  // reserve space for atom, bond, angle, and dihedral data
  const size_t nAtomsSz = static_cast<size_t>(this->nAtoms);
  this->atomIds.reserve(nAtomsSz);
  this->atomTypes.reserve(nAtomsSz);
  this->moleculeIds.reserve(nAtomsSz);
  this->atomX.reserve(nAtomsSz);
  this->atomY.reserve(nAtomsSz);
  this->atomZ.reserve(nAtomsSz);
  this->atomNx.reserve(nAtomsSz);
  this->atomNy.reserve(nAtomsSz);
  this->atomNz.reserve(nAtomsSz);

  const size_t nBondsSz = static_cast<size_t>(this->nBonds);
  this->bondIds.reserve(nBondsSz);
  this->bondTypes.reserve(nBondsSz);
  this->bondFrom.reserve(nBondsSz);
  this->bondTo.reserve(nBondsSz);

  const size_t nAnglesSz = static_cast<size_t>(this->nAngles);
  this->angleIds.reserve(nAnglesSz);
  this->angleTypes.reserve(nAnglesSz);
  this->angleFrom.reserve(nAnglesSz);
  this->angleVia.reserve(nAnglesSz);
  this->angleTo.reserve(nAnglesSz);

  const size_t nDihedralAnglesSz = static_cast<size_t>(this->nDihedralAngles);
  this->dihedralAngleIds.reserve(nDihedralAnglesSz);
  this->dihedralAngleTypes.reserve(nDihedralAnglesSz);
  this->dihedralAngleFrom.reserve(nDihedralAnglesSz);
  this->dihedralAngleVia1.reserve(nDihedralAnglesSz);
  this->dihedralAngleVia2.reserve(nDihedralAnglesSz);
  this->dihedralAngleTo.reserve(nDihedralAnglesSz);

  std::vector<long int> pendingVelAtomIds;
  std::vector<double> pendingVx;
  std::vector<double> pendingVy;
  std::vector<double> pendingVz;

  while (hasPendingSection) {
    const std::string_view sectionName =
      pylimer_tools::utils::trimWhitespaceView(
        pylimer_tools::utils::trimLineOmitCommentView(rawLine));
    hasPendingSection = false;

    if (sectionName == "Masses") {
      while (std::getline(file, rawLine)) {
        const std::string_view cleaned =
          pylimer_tools::utils::trimWhitespaceView(
            pylimer_tools::utils::trimLineOmitCommentView(rawLine));
        if (cleaned.empty()) {
          continue;
        }
        if (isDataFileSectionHeader(cleaned)) {
          hasPendingSection = true;
          break;
        }
        this->readMass(cleaned);
      }
    } else if (sectionName == "Atoms") {
      AtomStyle effStyle1 = atomStyle;
      AtomStyle effStyle2 = atomStyle2;
      AtomStyle effStyle3 = atomStyle3;
      const size_t hashPos = rawLine.find('#');
      if (hashPos != std::string::npos) {
        const std::vector<std::string_view> styleTokens =
          pylimer_tools::utils::splitWhitespaceView(
            std::string_view(rawLine).substr(hashPos + 1));
        if (!styleTokens.empty()) {
          effStyle1 = pylimer_tools::utils::getAtomStyleFromString(
            std::string(styleTokens[0]));
          if (effStyle1 == AtomStyle::HYBRID) {
            if (styleTokens.size() >= 2) {
              effStyle2 = pylimer_tools::utils::getAtomStyleFromString(
                std::string(styleTokens[1]));
            }
            if (styleTokens.size() >= 3) {
              effStyle3 = pylimer_tools::utils::getAtomStyleFromString(
                std::string(styleTokens[2]));
            }
          }
        }
      }
      this->readAtomsSection(file, effStyle1, effStyle2, effStyle3);
    } else if (sectionName == "Bonds") {
      this->readBonds(file);
    } else if (sectionName == "Angles") {
      this->readAngles(file);
    } else if (sectionName == "Dihedrals") {
      this->readDihedralAngles(file);
    } else if (sectionName == "Velocities") {
      this->readVelocities(
        file, pendingVelAtomIds, pendingVx, pendingVy, pendingVz);
    }

    if (!hasPendingSection) {
      while (std::getline(file, rawLine)) {
        const std::string_view cleaned =
          pylimer_tools::utils::trimWhitespaceView(
            pylimer_tools::utils::trimLineOmitCommentView(rawLine));
        if (isDataFileSectionHeader(cleaned)) {
          hasPendingSection = true;
          break;
        }
      }
    }
  }

  RUNTIME_EXP_IFN(static_cast<int>(this->atomIds.size()) == this->nAtoms,
                  "Data file ended too early. Not able to read any atoms.");

  if (!pendingVelAtomIds.empty()) {
    this->finalizePendingVelocities(
      pendingVelAtomIds, pendingVx, pendingVy, pendingVz);
  }
}

void
DataFileParser::readNs(const std::string_view line)
{
  const auto parseOneInt = [&line]() -> int {
    const std::string lineStr(line);
    const char* cursor = lineStr.c_str();
    int val = 0;
    RUNTIME_EXP_IFN(pylimer_tools::utils::parseNumericToken(cursor, val),
                    "Failed to parse count in header line '" + lineStr + "'.");
    return val;
  };
  const auto parseTwoDoubles = [&line](double& lo, double& hi) {
    const std::string lineStr(line);
    const char* cursor = lineStr.c_str();
    RUNTIME_EXP_IFN(pylimer_tools::utils::parseNumericToken(cursor, lo) &&
                      pylimer_tools::utils::parseNumericToken(cursor, hi),
                    "Failed to parse box bounds in header line '" + lineStr +
                      "'.");
  };

  if (line.find("atoms") != std::string_view::npos) {
    this->nAtoms = parseOneInt();
  } else if (line.find("bonds") != std::string_view::npos) {
    this->nBonds = parseOneInt();
  } else if (line.find("angles") != std::string_view::npos) {
    this->nAngles = parseOneInt();
  } else if (line.find("dihedrals") != std::string_view::npos) {
    this->nDihedralAngles = parseOneInt();
  } else if (line.find("atom types") != std::string_view::npos) {
    this->nAtomTypes = parseOneInt();
  } else if (line.find("bond types") != std::string_view::npos) {
    this->nBondTypes = parseOneInt();
  } else if (line.find("angle types") != std::string_view::npos) {
    this->nAngleTypes = parseOneInt();
  } else if (line.find("dihedral types") != std::string_view::npos) {
    this->nDihedralAngleTypes = parseOneInt();
  } else if (line.find("xlo xhi") != std::string_view::npos) {
    parseTwoDoubles(this->xLo, this->xHi);
  } else if (line.find("ylo yhi") != std::string_view::npos) {
    parseTwoDoubles(this->yLo, this->yHi);
  } else if (line.find("zlo zhi") != std::string_view::npos) {
    parseTwoDoubles(this->zLo, this->zHi);
  }
}

void
DataFileParser::readMass(const std::string_view line)
{
  const std::string lineStr(line);
  const char* cursor = lineStr.c_str();
  int key = 0;
  double massVal = 0.0;
  if (!pylimer_tools::utils::parseNumericToken(cursor, key) ||
      !pylimer_tools::utils::parseNumericToken(cursor, massVal) ||
      hasRemainingTokens(cursor)) {
    throw std::runtime_error(
      "Incorrect nr of fields tokenized when reading masses");
  }
  this->masses[key] = massVal;
}

const DataFileParser::AtomFieldDescriptor&
DataFileParser::getAtomStyleDescriptor(const AtomStyle style)
{
  using FT = AtomFieldDescriptor::FieldType;

  switch (style) {
    case AtomStyle::ANGLE:
    case AtomStyle::BOND:
    case AtomStyle::MOLECULAR: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::MOLECULE_ID, FT::ATOM_TYPE, FT::X, FT::Y,
        FT::Z,       FT::NX,          FT::NY,        FT::NZ
      };
      return desc;
    }

    case AtomStyle::ATOMIC: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::X, FT::Y, FT::Z, FT::NX, FT::NY, FT::NZ
      };
      return desc;
    }

    case AtomStyle::BODY: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::BODYFLAG, FT::MASS, FT::X,
        FT::Y,       FT::Z,         FT::NX,       FT::NY,   FT::NZ
      };
      return desc;
    }

    case AtomStyle::BPM_SPHERE: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::MOLECULE_ID, FT::ATOM_TYPE, FT::DIAMETER,
        FT::DENSITY, FT::X,           FT::Y,         FT::Z,
        FT::NX,      FT::NY,          FT::NZ
      };
      return desc;
    }

    case AtomStyle::CHARGE: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::ATOM_CHARGE, FT::X, FT::Y,
        FT::Z,       FT::NX,        FT::NY,          FT::NZ
      };
      return desc;
    }

    case AtomStyle::DIELECTRIC: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::ATOM_CHARGE, FT::X,    FT::Y,  FT::Z,
        FT::MUX,     FT::MUY,       FT::MUZ,         FT::AREA, FT::ED, FT::EM,
        FT::EPSILON, FT::CURVATURE, FT::NX,          FT::NY,   FT::NZ
      };
      return desc;
    }

    case AtomStyle::DIPOLE: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::ATOM_CHARGE, FT::X,  FT::Y,  FT::Z,
        FT::MUX,     FT::MUY,       FT::MUZ,         FT::NX, FT::NY, FT::NZ
      };
      return desc;
    }

    case AtomStyle::DPD: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::THETA, FT::X, FT::Y,
        FT::Z,       FT::NX,        FT::NY,    FT::NZ
      };
      return desc;
    }

    case AtomStyle::EDPD: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::EDPD_TEMP, FT::EDPD_CV, FT::X,
        FT::Y,       FT::Z,         FT::NX,        FT::NY,      FT::NZ
      };
      return desc;
    }

    case AtomStyle::ELECTRON: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::ATOM_CHARGE, FT::ESPIN,
        FT::ERADIUS, FT::X,         FT::Y,           FT::Z,
        FT::NX,      FT::NY,        FT::NZ
      };
      return desc;
    }

    case AtomStyle::ELLIPSOID: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::ELLIPSOIDFLAG,
        FT::DENSITY, FT::X,         FT::Y,
        FT::Z,       FT::NX,        FT::NY,
        FT::NZ
      };
      return desc;
    }

    case AtomStyle::LINE: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::MOLECULE_ID, FT::ATOM_TYPE, FT::LINEFLAG,
        FT::DENSITY, FT::X,           FT::Y,         FT::Z,
        FT::NX,      FT::NY,          FT::NZ
      };
      return desc;
    }

    case AtomStyle::MDPD: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::RHO, FT::X, FT::Y,
        FT::Z,       FT::NX,        FT::NY,  FT::NZ
      };
      return desc;
    }

    case AtomStyle::PERI: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::VOLUME, FT::DENSITY, FT::X,
        FT::Y,       FT::Z,         FT::NX,     FT::NY,      FT::NZ
      };
      return desc;
    }

    case AtomStyle::RHEO: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::STATUS, FT::RHO, FT::X,
        FT::Y,       FT::Z,         FT::NX,     FT::NY,  FT::NZ
      };
      return desc;
    }

    case AtomStyle::RHEO_THERMAL: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::STATUS, FT::RHO, FT::ENERGY, FT::X,
        FT::Y,       FT::Z,         FT::NX,     FT::NY,  FT::NZ
      };
      return desc;
    }

    case AtomStyle::SMD: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::MOLECULE_ID, FT::VOLUME,
        FT::MASS,    FT::KRADIUS,   FT::CRADIUS,     FT::X0,
        FT::Y0,      FT::Z0,        FT::X,           FT::Y,
        FT::Z,       FT::NX,        FT::NY,          FT::NZ
      };
      return desc;
    }

    case AtomStyle::SPH: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::RHO, FT::ESPH, FT::CV, FT::X,
        FT::Y,       FT::Z,         FT::NX,  FT::NY,   FT::NZ
      };
      return desc;
    }

    case AtomStyle::SPHERE: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::DIAMETER, FT::DENSITY, FT::X,
        FT::Y,       FT::Z,         FT::NX,       FT::NY,      FT::NZ
      };
      return desc;
    }

    case AtomStyle::SPIN: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::X,  FT::Y,  FT::Z,  FT::SPX,
        FT::SPY,     FT::SPZ,       FT::SP, FT::NX, FT::NY, FT::NZ
      };
      return desc;
    }

    case AtomStyle::TDPD: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::X,  FT::Y,  FT::Z,
        FT::CC1,     FT::CC2,       FT::NX, FT::NY, FT::NZ
      };
      return desc;
    }

    case AtomStyle::TEMPLATE: {
      static const AtomFieldDescriptor desc = { FT::ATOM_ID,
                                                FT::ATOM_TYPE,
                                                FT::MOLECULE_ID,
                                                FT::TEMPLATE_INDEX,
                                                FT::TEMPLATE_ATOM,
                                                FT::X,
                                                FT::Y,
                                                FT::Z,
                                                FT::NX,
                                                FT::NY,
                                                FT::NZ };
      return desc;
    }

    case AtomStyle::TRI: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::MOLECULE_ID, FT::ATOM_TYPE, FT::TRIANGLEFLAG,
        FT::DENSITY, FT::X,           FT::Y,         FT::Z,
        FT::NX,      FT::NY,          FT::NZ
      };
      return desc;
    }

    case AtomStyle::WAVEPACKET: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::ATOM_TYPE, FT::ATOM_CHARGE, FT::ESPIN, FT::ERADIUS,
        FT::ETAG,    FT::CS_RE,     FT::CS_IM,       FT::X,     FT::Y,
        FT::Z,       FT::NX,        FT::NY,          FT::NZ
      };
      return desc;
    }

    case AtomStyle::FULL: {
      static const AtomFieldDescriptor desc = {
        FT::ATOM_ID, FT::MOLECULE_ID, FT::ATOM_TYPE, FT::ATOM_CHARGE, FT::X,
        FT::Y,       FT::Z,           FT::NX,        FT::NY,          FT::NZ
      };
      return desc;
    }

    default:
      throw std::invalid_argument("Unsupported atom style (" +
                                  getAtomStyleString(style) +
                                  ") for generic parsing");
  }
}

void
DataFileParser::readAtomsSection(std::ifstream& file,
                                 const AtomStyle atomStyle,
                                 const AtomStyle atomStyle2,
                                 const AtomStyle atomStyle3)
{
  if (this->nAtoms <= 0) {
    return;
  }

  std::string line;
  if (atomStyle == AtomStyle::HYBRID) {
    for (int i = 0; i < this->nAtoms; ++i) {
      if (!nextNonEmptyDataLine(file, line)) {
        throw std::runtime_error(
          "Data file ended too early. Not enough atoms read. Read " +
          std::to_string(i) + " of " + std::to_string(this->nAtoms) + ".");
      }
      this->readAtomHybrid(line, atomStyle2, atomStyle3);
    }
    return;
  }

  const AtomFieldDescriptor& descriptor = getAtomStyleDescriptor(atomStyle);
  using FT = AtomFieldDescriptor::FieldType;

  const auto getExtraFieldKey = [](const FT fieldType) -> const char* {
    switch (fieldType) {
      case FT::ATOM_ID:
      case FT::MOLECULE_ID:
      case FT::ATOM_TYPE:
      case FT::X:
      case FT::Y:
      case FT::Z:
      case FT::NX:
      case FT::NY:
      case FT::NZ:
        return nullptr;
      case FT::ATOM_CHARGE:
        return "charge";
      case FT::MUX:
        return "mux";
      case FT::MUY:
        return "muy";
      case FT::MUZ:
        return "muz";
      case FT::ATOM_EDPD:
        return "edpd";
      case FT::THETA:
        return "theta";
      case FT::RHO:
        return "rho";
      case FT::DIAMETER:
        return "diameter";
      case FT::DENSITY:
        return "density";
      case FT::EDPD_TEMP:
        return "edpd_temp";
      case FT::EDPD_CV:
        return "edpd_cv";
      case FT::ESPH:
        return "esph";
      case FT::CV:
        return "cv";
      case FT::ENERGY:
        return "energy";
      case FT::BODYFLAG:
        return "bodyflag";
      case FT::ELLIPSOIDFLAG:
        return "ellipsoidflag";
      case FT::LINEFLAG:
        return "lineflag";
      case FT::TRIANGLEFLAG:
        return "triangleflag";
      case FT::STATUS:
        return "status";
      case FT::KRADIUS:
        return "kradius";
      case FT::CRADIUS:
        return "cradius";
      case FT::ERADIUS:
        return "eradius";
      case FT::MASS:
        return "mass";
      case FT::VOLUME:
        return "volume";
      case FT::X0:
        return "x0";
      case FT::Y0:
        return "y0";
      case FT::Z0:
        return "z0";
      case FT::SPX:
        return "spx";
      case FT::SPY:
        return "spy";
      case FT::SPZ:
        return "spz";
      case FT::SP:
        return "sp";
      case FT::ESPIN:
        return "espin";
      case FT::AREA:
        return "area";
      case FT::ED:
        return "ed";
      case FT::EM:
        return "em";
      case FT::EPSILON:
        return "epsilon";
      case FT::CURVATURE:
        return "curvature";
      case FT::TEMPLATE_INDEX:
        return "template_index";
      case FT::TEMPLATE_ATOM:
        return "template_atom";
      case FT::ETAG:
        return "etag";
      case FT::CS_RE:
        return "cs_re";
      case FT::CS_IM:
        return "cs_im";
      case FT::CC1:
        return "cc1";
      case FT::CC2:
        return "cc2";
      default:
        throw std::invalid_argument(
          "Unknown field type in AtomFieldDescriptor: " +
          std::to_string(static_cast<int>(fieldType)));
    }
  };

  const size_t nAtomsSz = static_cast<size_t>(this->nAtoms);
  const size_t numFields = descriptor.fields.size();
  std::vector<std::vector<double>*> extraFieldVecs(numFields, nullptr);
  for (size_t f = 0; f < numFields; ++f) {
    const char* key = getExtraFieldKey(descriptor.fields[f]);
    if (key != nullptr) {
      std::vector<double>& vec = this->additionalAtomData[key];
      vec.reserve(nAtomsSz);
      extraFieldVecs[f] = &vec;
    }
  }

  const int minRequiredFields = static_cast<int>(numFields) - 3;
  const int fullRequiredFields = static_cast<int>(numFields);

  for (int i = 0; i < this->nAtoms; ++i) {
    if (!nextNonEmptyDataLine(file, line)) {
      throw std::runtime_error(
        "Data file ended too early. Not enough atoms read. Read " +
        std::to_string(i) + " of " + std::to_string(this->nAtoms) + ".");
    }

    const char* cursor = line.c_str();
    int resFound = 0;

    for (size_t f = 0; f < numFields; ++f) {
      if (extraFieldVecs[f] != nullptr) {
        double val = 0.0;
        if (!pylimer_tools::utils::parseNumericToken(cursor, val)) {
          break;
        }
        extraFieldVecs[f]->push_back(val);
        ++resFound;
        continue;
      }

      bool ok = false;
      switch (descriptor.fields[f]) {
        case FT::ATOM_ID: {
          long int val = 0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->atomIds.push_back(val);
          }
          break;
        }
        case FT::MOLECULE_ID: {
          int val = 0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->moleculeIds.push_back(val);
          }
          break;
        }
        case FT::ATOM_TYPE: {
          int val = 0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->atomTypes.push_back(val);
          }
          break;
        }
        case FT::X: {
          double val = 0.0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->atomX.push_back(val);
          }
          break;
        }
        case FT::Y: {
          double val = 0.0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->atomY.push_back(val);
          }
          break;
        }
        case FT::Z: {
          double val = 0.0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->atomZ.push_back(val);
          }
          break;
        }
        case FT::NX: {
          int val = 0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->atomNx.push_back(val);
          }
          break;
        }
        case FT::NY: {
          int val = 0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->atomNy.push_back(val);
          }
          break;
        }
        case FT::NZ: {
          int val = 0;
          ok = pylimer_tools::utils::parseNumericToken(cursor, val);
          if (ok) {
            this->atomNz.push_back(val);
          }
          break;
        }
        default:
          break;
      }
      if (!ok) {
        break;
      }
      ++resFound;
    }

    RUNTIME_EXP_IFN(
      (resFound == minRequiredFields || resFound == fullRequiredFields) &&
        !hasRemainingTokens(cursor),
      "Invalid atom format in line '" + line + "': parsed " +
        std::to_string(resFound) + " fields.");

    if (resFound == minRequiredFields) {
      this->atomNx.push_back(0);
      this->atomNy.push_back(0);
      this->atomNz.push_back(0);
    }

    // Set default molecule ID for styles that don't have it
    if (!descriptor.hasMoleculeId) {
      this->moleculeIds.push_back(0);
    }
  }
}

void
DataFileParser::readAtomHybrid(const std::string_view line,
                               const AtomStyle style1,
                               const AtomStyle style2)
{
  const bool validBondEdpd =
    (style1 == AtomStyle::BOND && style2 == AtomStyle::EDPD) ||
    (style1 == AtomStyle::EDPD && style2 == AtomStyle::BOND);
  if (!validBondEdpd) {
    throw std::runtime_error(
      "This combination is not implemented for hybrid atom style");
  }

  const std::string lineStr(line);
  const char* cursor = lineStr.c_str();
  long int atomId = 0;
  int atomType = 0;
  int moleculeId = 0;
  double x = 0.0, y = 0.0, z = 0.0;
  double edpdTemp = 0.0, edpd = 0.0;
  const bool baseOk =
    pylimer_tools::utils::parseNumericToken(cursor, atomId) &&
    pylimer_tools::utils::parseNumericToken(cursor, atomType) &&
    pylimer_tools::utils::parseNumericToken(cursor, x) &&
    pylimer_tools::utils::parseNumericToken(cursor, y) &&
    pylimer_tools::utils::parseNumericToken(cursor, z) &&
    pylimer_tools::utils::parseNumericToken(cursor, moleculeId) &&
    pylimer_tools::utils::parseNumericToken(cursor, edpdTemp) &&
    pylimer_tools::utils::parseNumericToken(cursor, edpd);
  RUNTIME_EXP_IFN(baseOk,
                  "Did not find enough data in line '" + lineStr + "'.");

  int nx = 0, ny = 0, nz = 0;
  if (hasRemainingTokens(cursor)) {
    const bool imgOk = pylimer_tools::utils::parseNumericToken(cursor, nx) &&
                       pylimer_tools::utils::parseNumericToken(cursor, ny) &&
                       pylimer_tools::utils::parseNumericToken(cursor, nz) &&
                       !hasRemainingTokens(cursor);
    RUNTIME_EXP_IFN(imgOk,
                    "Invalid trailing image flags in hybrid atom line '" +
                      lineStr + "'.");
  }

  this->atomIds.push_back(atomId);
  this->moleculeIds.push_back(moleculeId);
  this->atomTypes.push_back(atomType);
  this->atomX.push_back(x);
  this->atomY.push_back(y);
  this->atomZ.push_back(z);
  this->additionalAtomData["edpd_temp"].push_back(edpdTemp);
  this->additionalAtomData["edpd"].push_back(edpd);
  this->atomNx.push_back(nx);
  this->atomNy.push_back(ny);
  this->atomNz.push_back(nz);
}

template<typename ParseFn>
void
DataFileParser::readSectionLines(std::ifstream& file,
                                 const int count,
                                 const char* sectionName,
                                 ParseFn&& parseLine)
{
  if (count <= 0) {
    return;
  }

  std::string line;
  for (int i = 0; i < count; ++i) {
    if (!nextNonEmptyDataLine(file, line)) {
      if (i == 0) {
        throw std::runtime_error(
          std::string("Data file ended too early. Not able to read any ") +
          sectionName + ".");
      }
      throw std::runtime_error(
        std::string("Data file ended too early. Not enough ") + sectionName +
        " read. Read " + std::to_string(i) + " of " + std::to_string(count) +
        " expected " + sectionName + ".");
    }
    parseLine(line);
  }
}

void
DataFileParser::readBonds(std::ifstream& file)
{
  this->readSectionLines(file,
                         this->nBonds,
                         "bonds",
                         [this](const std::string& l) { this->readBond(l); });
}

void
DataFileParser::readBond(const std::string_view line)
{
  const std::string lineStr(line);
  const char* cursor = lineStr.c_str();
  long int bondId = 0, newBondFrom = 0, newBondTo = 0;
  int bondType = 0;
  const bool ok =
    pylimer_tools::utils::parseNumericToken(cursor, bondId) &&
    pylimer_tools::utils::parseNumericToken(cursor, bondType) &&
    pylimer_tools::utils::parseNumericToken(cursor, newBondFrom) &&
    pylimer_tools::utils::parseNumericToken(cursor, newBondTo);
  RUNTIME_EXP_IFN(ok, "Invalid bond line '" + lineStr + "'.");
  this->bondIds.push_back(bondId);
  this->bondTypes.push_back(bondType);
  this->bondFrom.push_back(newBondFrom);
  this->bondTo.push_back(newBondTo);
}

void
DataFileParser::readAngles(std::ifstream& file)
{
  this->readSectionLines(file,
                         this->nAngles,
                         "angles",
                         [this](const std::string& l) { this->readAngle(l); });
}

void
DataFileParser::readAngle(const std::string_view line)
{
  const std::string lineStr(line);
  const char* cursor = lineStr.c_str();
  long int newAngleId = 0, newAngleFrom = 0, newAngleVia = 0, newAngleTo = 0;
  int newAngleType = 0;
  const bool ok =
    pylimer_tools::utils::parseNumericToken(cursor, newAngleId) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleType) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleFrom) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleVia) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleTo);
  RUNTIME_EXP_IFN(ok, "Invalid angle line '" + lineStr + "'.");

  this->angleIds.push_back(newAngleId);
  this->angleTypes.push_back(newAngleType);
  this->angleFrom.push_back(newAngleFrom);
  this->angleVia.push_back(newAngleVia);
  this->angleTo.push_back(newAngleTo);
}

void
DataFileParser::readDihedralAngles(std::ifstream& file)
{
  this->readSectionLines(
    file,
    this->nDihedralAngles,
    "dihedral angles",
    [this](const std::string& l) { this->readDihedralAngle(l); });
}

void
DataFileParser::readDihedralAngle(const std::string_view line)
{
  const std::string lineStr(line);
  const char* cursor = lineStr.c_str();
  long int newAngleId = 0, newAngleFrom = 0, newAngleVia1 = 0, newAngleVia2 = 0,
           newAngleTo = 0;
  int newAngleType = 0;
  const bool ok =
    pylimer_tools::utils::parseNumericToken(cursor, newAngleId) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleType) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleFrom) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleVia1) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleVia2) &&
    pylimer_tools::utils::parseNumericToken(cursor, newAngleTo);
  RUNTIME_EXP_IFN(ok, "Invalid dihedral angle line '" + lineStr + "'.");

  this->dihedralAngleIds.push_back(newAngleId);
  this->dihedralAngleTypes.push_back(newAngleType);
  this->dihedralAngleFrom.push_back(newAngleFrom);
  this->dihedralAngleVia1.push_back(newAngleVia1);
  this->dihedralAngleVia2.push_back(newAngleVia2);
  this->dihedralAngleTo.push_back(newAngleTo);
}

void
DataFileParser::readVelocities(std::ifstream& file,
                               std::vector<long int>& pendingVelAtomIds,
                               std::vector<double>& pendingVx,
                               std::vector<double>& pendingVy,
                               std::vector<double>& pendingVz)
{
  if (this->nAtoms <= 0) {
    return;
  }

  const size_t nAtomsSz = static_cast<size_t>(this->nAtoms);
  const bool atomsAlreadyRead = (this->atomIds.size() == nAtomsSz);

  if (!atomsAlreadyRead) {
    pendingVelAtomIds.reserve(nAtomsSz);
    pendingVx.reserve(nAtomsSz);
    pendingVy.reserve(nAtomsSz);
    pendingVz.reserve(nAtomsSz);

    std::string line;
    for (size_t i = 0; i < nAtomsSz; ++i) {
      if (!nextNonEmptyDataLine(file, line)) {
        if (i == 0) {
          throw std::runtime_error("Data file ended too early. Not able to "
                                   "read any advertised velocities.");
        }
        throw std::runtime_error("Data file ended too early. Not enough "
                                 "velocities read. Read " +
                                 std::to_string(i) + " of " +
                                 std::to_string(this->nAtoms) +
                                 " expected velocities.");
      }
      const char* cursor = line.c_str();
      long int atomId = 0;
      double vx = 0.0, vy = 0.0, vz = 0.0;
      const bool ok = pylimer_tools::utils::parseNumericToken(cursor, atomId) &&
                      pylimer_tools::utils::parseNumericToken(cursor, vx) &&
                      pylimer_tools::utils::parseNumericToken(cursor, vy) &&
                      pylimer_tools::utils::parseNumericToken(cursor, vz);
      RUNTIME_EXP_IFN(ok, "Invalid velocity line '" + line + "'.");
      pendingVelAtomIds.push_back(atomId);
      pendingVx.push_back(vx);
      pendingVy.push_back(vy);
      pendingVz.push_back(vz);
    }
    return;
  }

  std::vector<double> orderedVx(nAtomsSz, 0.0);
  std::vector<double> orderedVy(nAtomsSz, 0.0);
  std::vector<double> orderedVz(nAtomsSz, 0.0);
  std::unordered_map<long int, size_t> atomIdToIdx;
  std::vector<uint8_t> seen;
  bool inOrder = true;

  std::string line;
  for (size_t i = 0; i < nAtomsSz; ++i) {
    if (!nextNonEmptyDataLine(file, line)) {
      if (i == 0) {
        throw std::runtime_error("Data file ended too early. Not able to read "
                                 "any advertised velocities.");
      }
      throw std::runtime_error("Data file ended too early. Not enough "
                               "velocities read. Read " +
                               std::to_string(i) + " of " +
                               std::to_string(this->nAtoms) +
                               " expected velocities.");
    }
    const char* cursor = line.c_str();
    long int atomId = 0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    const bool ok = pylimer_tools::utils::parseNumericToken(cursor, atomId) &&
                    pylimer_tools::utils::parseNumericToken(cursor, vx) &&
                    pylimer_tools::utils::parseNumericToken(cursor, vy) &&
                    pylimer_tools::utils::parseNumericToken(cursor, vz);
    RUNTIME_EXP_IFN(ok, "Invalid velocity line '" + line + "'.");

    if (inOrder && atomId == this->atomIds[i]) {
      orderedVx[i] = vx;
      orderedVy[i] = vy;
      orderedVz[i] = vz;
    } else {
      if (inOrder) {
        inOrder = false;
        atomIdToIdx.reserve(nAtomsSz);
        for (size_t k = 0; k < nAtomsSz; ++k) {
          atomIdToIdx.emplace(this->atomIds[k], k);
        }
        seen.assign(nAtomsSz, 0);
        for (size_t k = 0; k < i; ++k) {
          seen[k] = 1;
        }
      }
      const auto it = atomIdToIdx.find(atomId);
      RUNTIME_EXP_IFN(it != atomIdToIdx.end(),
                      "Velocity entry references unknown atom ID " +
                        std::to_string(atomId) + ".");
      const size_t targetIdx = it->second;
      orderedVx[targetIdx] = vx;
      orderedVy[targetIdx] = vy;
      orderedVz[targetIdx] = vz;
      seen[targetIdx] = 1;
    }
  }

  if (!inOrder) {
    for (size_t k = 0; k < nAtomsSz; ++k) {
      RUNTIME_EXP_IFN(seen[k] != 0,
                      "Missing velocity entry for atom ID " +
                        std::to_string(this->atomIds[k]) + ".");
    }
  }

  this->additionalAtomData["vx"] = std::move(orderedVx);
  this->additionalAtomData["vy"] = std::move(orderedVy);
  this->additionalAtomData["vz"] = std::move(orderedVz);
}

void
DataFileParser::finalizePendingVelocities(
  const std::vector<long int>& pendingVelAtomIds,
  const std::vector<double>& pendingVx,
  const std::vector<double>& pendingVy,
  const std::vector<double>& pendingVz)
{
  const size_t nAtomsSz = this->atomIds.size();
  RUNTIME_EXP_IFN(pendingVelAtomIds.size() == nAtomsSz,
                  "Mismatch between velocity count and atom count.");

  std::unordered_map<long int, size_t> atomIdToIdx;
  atomIdToIdx.reserve(nAtomsSz);
  for (size_t i = 0; i < nAtomsSz; ++i) {
    atomIdToIdx.emplace(pendingVelAtomIds[i], i);
  }

  std::vector<double> orderedVx(nAtomsSz, 0.0);
  std::vector<double> orderedVy(nAtomsSz, 0.0);
  std::vector<double> orderedVz(nAtomsSz, 0.0);
  for (size_t i = 0; i < nAtomsSz; ++i) {
    const auto it = atomIdToIdx.find(this->atomIds[i]);
    RUNTIME_EXP_IFN(it != atomIdToIdx.end(),
                    "Missing velocity entry for atom ID " +
                      std::to_string(this->atomIds[i]) + ".");
    const size_t realIdx = it->second;
    orderedVx[i] = pendingVx[realIdx];
    orderedVy[i] = pendingVy[realIdx];
    orderedVz[i] = pendingVz[realIdx];
  }

  this->additionalAtomData["vx"] = std::move(orderedVx);
  this->additionalAtomData["vy"] = std::move(orderedVy);
  this->additionalAtomData["vz"] = std::move(orderedVz);
}
} // namespace pylimer_tools::utils
