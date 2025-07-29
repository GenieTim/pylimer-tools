#include "DataFileParser.h"
#include "../utils/StringUtils.h"
#include "../utils/utilityMacros.h"
#include <filesystem>
#include <fstream> // std::ifstream
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace pylimer_tools::utils {

void
DataFileParser::read(const std::string& filePath,
                     AtomStyle atomStyle,
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
    // read the nr of data points to read afterward
    this->readNs(line);
  }

  // reserve space
  // for atom data
  this->atomIds.reserve(static_cast<size_t>(this->nAtoms));
  this->atomTypes.reserve(static_cast<size_t>(this->nAtoms));
  this->moleculeIds.reserve(static_cast<size_t>(this->nAtoms));
  this->atomX.reserve(static_cast<size_t>(this->nAtoms));
  this->atomY.reserve(static_cast<size_t>(this->nAtoms));
  this->atomZ.reserve(static_cast<size_t>(this->nAtoms));
  this->atomNx.reserve(static_cast<size_t>(this->nAtoms));
  this->atomNy.reserve(static_cast<size_t>(this->nAtoms));
  this->atomNz.reserve(static_cast<size_t>(this->nAtoms));
  // and bond data
  this->bondIds.reserve(static_cast<size_t>(this->nBonds));
  this->bondTypes.reserve(static_cast<size_t>(this->nBonds));
  this->bondFrom.reserve(static_cast<size_t>(this->nBonds));
  this->bondTo.reserve(static_cast<size_t>(this->nBonds));

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
  // detect atom style
  if (pylimer_tools::utils::contains(line, "#")) {
    std::string atomStyleString = pylimer_tools::utils::trim(
      pylimer_tools::utils::removeAllRegex(line, "Atoms[ ]+#", true));
    if (atomStyleString.size() > 2) {
      atomStyle = pylimer_tools::utils::getAtomStyleFromString(atomStyleString);
    }
  }

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
      case AtomStyle::ATOMIC:
        this->readAtomAtomic(line);
        break;
      case AtomStyle::BODY:
        this->readAtomBody(line);
        break;
      case AtomStyle::BPM_SPHERE:
        this->readAtomBpmSphere(line);
        break;
      case AtomStyle::CHARGE:
        this->readAtomCharge(line);
        break;
      case AtomStyle::DIELECTRIC:
        this->readAtomDielectric(line);
        break;
      case AtomStyle::DIPOLE:
        this->readAtomDipole(line);
        break;
      case AtomStyle::DPD:
        this->readAtomDpd(line);
        break;
      case AtomStyle::EDPD:
        this->readAtomEdpd(line);
        break;
      case AtomStyle::ELECTRON:
        this->readAtomElectron(line);
        break;
      case AtomStyle::ELLIPSOID:
        this->readAtomEllipsoid(line);
        break;
      case AtomStyle::FULL:
        this->readAtomFull(line);
        break;
      case AtomStyle::LINE:
        this->readAtomLine(line);
        break;
      case AtomStyle::MDPD:
        this->readAtomMdpd(line);
        break;
      case AtomStyle::PERI:
        this->readAtomPeri(line);
        break;
      case AtomStyle::RHEO:
        this->readAtomRheo(line);
        break;
      case AtomStyle::RHEO_THERMAL:
        this->readAtomRheoThermal(line);
        break;
      case AtomStyle::SMD:
        this->readAtomSmd(line);
        break;
      case AtomStyle::SPH:
        this->readAtomSph(line);
        break;
      case AtomStyle::SPHERE:
        this->readAtomSphere(line);
        break;
      case AtomStyle::SPIN:
        this->readAtomSpin(line);
        break;
      case AtomStyle::TEMPLATE:
        this->readAtomTemplate(line);
        break;
      case AtomStyle::TRI:
        this->readAtomTri(line);
        break;
      case AtomStyle::WAVEPACKET:
        this->readAtomWavepacket(line);
        break;
      case AtomStyle::HYBRID:
        this->readAtomHybrid(line, atomStyle2, atomStyle3);
        break;
      default:
        throw std::invalid_argument(
          "This atom style (" +
          pylimer_tools::utils::getAtomStyleString(atomStyle) +
          ") is not supported yet.");
    }

    if (!getline(file, line) && i < this->nAtoms - 1) {
      throw std::runtime_error(
        "Data file ended too early. Not enough atoms read. Read " +
        std::to_string(i) + " of " + std::to_string(this->nAtoms) + ".");
    }
  }

  // read the rest of the file
  while (file.peek() != EOF) {
    this->skipLinesToContains(
      line, file, { "Bonds", "Angles", "Dihedrals", "Velocities" });

    if (file.peek() == EOF) {
      break;
    }

    if (contains(line, "Bonds")) {
      this->readBonds(file, line);
    } else if (contains(line, "Angles")) {
      this->readAngles(file, line);
    } else if (contains(line, "Dihedrals")) {
      this->readDihedralAngles(file, line);
    } else if (contains(line, "Velocities")) {
      this->readVelocities(file, line);
    }
  }

  // we ignore dihedrals etc. for now.
  file.close();
}

void
DataFileParser::skipLinesToContains(std::string& line,
                                    std::ifstream& file,
                                    const std::string& upTo)
{
  do {
    if (contains(line, upTo)) {
      break;
    }
  } while (getline(file, line));
}

void
DataFileParser::skipLinesToContains(
  std::string& line,
  std::ifstream& file,
  const std::vector<std::string>& upToEitherOr)
{
  do {
    for (const std::string& upTo : upToEitherOr) {
      if (contains(line, upTo)) {
        return;
      }
    }
  } while (getline(file, line));
};

void
DataFileParser::skipEmptyLines(std::string& line, std::ifstream& file)
{
  do {
    line = pylimer_tools::utils::trimLineOmitComment(line);

    // skip until empty lines
    if (!line.empty()) {
      break;
    }
  } while (getline(file, line));
}

void
DataFileParser::readNs(const std::string& line)
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

void
DataFileParser::readMass(const std::string& line)
{
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

DataFileParser::AtomFieldDescriptor
DataFileParser::getAtomStyleDescriptor(AtomStyle style) const
{
  using FT = AtomFieldDescriptor::FieldType;

  switch (style) {
    case AtomStyle::ANGLE:
    case AtomStyle::BOND:
    case AtomStyle::MOLECULAR:
      return { { FT::ATOM_ID,
                 FT::MOLECULE_ID,
                 FT::ATOM_TYPE,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %le %le %le %lu %lu %lu" };

    case AtomStyle::ATOMIC:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %le %le %le %ld %ld %ld" };

    case AtomStyle::BODY:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::BODYFLAG,
                 FT::MASS,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %le %le %le %le %ld %ld %ld" };

    case AtomStyle::BPM_SPHERE:
      return { { FT::ATOM_ID,
                 FT::MOLECULE_ID,
                 FT::ATOM_TYPE,
                 FT::DIAMETER,
                 FT::DENSITY,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::CHARGE:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::ATOM_CHARGE,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%zd %d %le %le %le %le %ld %ld %ld" };

    case AtomStyle::DIELECTRIC:
      return {
        { FT::ATOM_ID,
          FT::ATOM_TYPE,
          FT::ATOM_CHARGE,
          FT::X,
          FT::Y,
          FT::Z,
          FT::MUX,
          FT::MUY,
          FT::MUZ,
          FT::AREA,
          FT::ED,
          FT::EM,
          FT::EPSILON,
          FT::CURVATURE,
          FT::NX,
          FT::NY,
          FT::NZ },
        "%lu %d %le %le %le %le %le %le %le %le %le %le %le %le %ld %ld %ld"
      };

    case AtomStyle::DIPOLE:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::ATOM_CHARGE,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::MUX,
                 FT::MUY,
                 FT::MUZ,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%zd %d %le %le %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::DPD:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::THETA,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%zd %d %le %le %le %le %ld %ld %ld" };

    case AtomStyle::EDPD:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::EDPD_TEMP,
                 FT::EDPD_CV,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::ELECTRON:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::ATOM_CHARGE,
                 FT::ESPIN,
                 FT::ERADIUS,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %le %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::ELLIPSOID:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::ELLIPSOIDFLAG,
                 FT::DENSITY,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %le %le %le %le %ld %ld %ld" };

    case AtomStyle::LINE:
      return { { FT::ATOM_ID,
                 FT::MOLECULE_ID,
                 FT::ATOM_TYPE,
                 FT::LINEFLAG,
                 FT::DENSITY,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %d %le %le %le %le %ld %ld %ld" };

    case AtomStyle::MDPD:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::RHO,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%zd %d %le %le %le %le %ld %ld %ld" };

    case AtomStyle::PERI:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::VOLUME,
                 FT::DENSITY,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::RHEO:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::STATUS,
                 FT::RHO,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %le %le %le %le %ld %ld %ld" };

    case AtomStyle::RHEO_THERMAL:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::STATUS,
                 FT::RHO,
                 FT::ENERGY,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::SMD:
      return {
        { FT::ATOM_ID,
          FT::ATOM_TYPE,
          FT::MOLECULE_ID,
          FT::VOLUME,
          FT::MASS,
          FT::KRADIUS,
          FT::CRADIUS,
          FT::X0,
          FT::Y0,
          FT::Z0,
          FT::X,
          FT::Y,
          FT::Z,
          FT::NX,
          FT::NY,
          FT::NZ },
        "%lu %d %d %le %le %le %le %le %le %le %le %le %le %ld %ld %ld"
      };

    case AtomStyle::SPH:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::RHO,
                 FT::ESPH,
                 FT::CV,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %le %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::SPHERE:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::DIAMETER,
                 FT::DENSITY,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%zd %d %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::SPIN:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::SPX,
                 FT::SPY,
                 FT::SPZ,
                 FT::SP,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %le %le %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::TDPD:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::CC1,
                 FT::CC2,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::TEMPLATE:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::MOLECULE_ID,
                 FT::TEMPLATE_INDEX,
                 FT::TEMPLATE_ATOM,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %d %d %le %le %le %ld %ld %ld" };

    case AtomStyle::TRI:
      return { { FT::ATOM_ID,
                 FT::MOLECULE_ID,
                 FT::ATOM_TYPE,
                 FT::TRIANGLEFLAG,
                 FT::DENSITY,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %d %d %le %le %le %le %ld %ld %ld" };

    case AtomStyle::WAVEPACKET:
      return { { FT::ATOM_ID,
                 FT::ATOM_TYPE,
                 FT::ATOM_CHARGE,
                 FT::ESPIN,
                 FT::ERADIUS,
                 FT::ETAG,
                 FT::CS_RE,
                 FT::CS_IM,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%lu %d %le %le %le %d %le %le %le %le %le %ld %ld %ld" };

    case AtomStyle::FULL:
      return { { FT::ATOM_ID,
                 FT::MOLECULE_ID,
                 FT::ATOM_TYPE,
                 FT::ATOM_CHARGE,
                 FT::X,
                 FT::Y,
                 FT::Z,
                 FT::NX,
                 FT::NY,
                 FT::NZ },
               "%zd %d %d %le %le %le %le %zd %zd %zd" };

    default:
      throw std::invalid_argument("Unsupported atom style for generic parsing");
  }
}

template<typename... Args>
void
DataFileParser::readAtomGeneric(const std::string& line,
                                const AtomFieldDescriptor& descriptor,
                                Args&... args)
{
  const int resFound =
    sscanf(line.c_str(), descriptor.format.c_str(), &args...);

  RUNTIME_EXP_IFN(resFound >= static_cast<int>(descriptor.fields.size()) - 3,
                  "Did not find enough data in line '" + line + "': only " +
                    std::to_string(resFound) + ".");

  // Map parsed values to storage vectors based on field types
  size_t argIdx = 0;
  auto processArg = [&](auto& arg) {
    if (argIdx < descriptor.fields.size()) {
      using FT = AtomFieldDescriptor::FieldType;
      switch (descriptor.fields[argIdx]) {
        case FT::ATOM_ID:
          this->atomIds.push_back(static_cast<long int>(arg));
          break;
        case FT::MOLECULE_ID:
          this->moleculeIds.push_back(static_cast<int>(arg));
          break;
        case FT::ATOM_TYPE:
          this->atomTypes.push_back(static_cast<int>(arg));
          break;
        case FT::X:
          this->atomX.push_back(arg);
          break;
        case FT::Y:
          this->atomY.push_back(arg);
          break;
        case FT::Z:
          this->atomZ.push_back(arg);
          break;
        case FT::NX:
          if (resFound >= static_cast<int>(descriptor.fields.size()) - 3)
            this->atomNx.push_back(static_cast<int>(arg));
          break;
        case FT::NY:
          if (resFound >= static_cast<int>(descriptor.fields.size()) - 2)
            this->atomNy.push_back(static_cast<int>(arg));
          break;
        case FT::NZ:
          if (resFound >= static_cast<int>(descriptor.fields.size()) - 1)
            this->atomNz.push_back(static_cast<int>(arg));
          break;
        case FT::ATOM_CHARGE:
          this->additionalAtomData["charge"].push_back(arg);
          break;
        case FT::MUX:
          this->additionalAtomData["mux"].push_back(arg);
          break;
        case FT::MUY:
          this->additionalAtomData["muy"].push_back(arg);
          break;
        case FT::MUZ:
          this->additionalAtomData["muz"].push_back(arg);
          break;
        case FT::ATOM_EDPD:
          this->additionalAtomData["edpd"].push_back(arg);
          break;
        case FT::THETA:
          this->additionalAtomData["theta"].push_back(arg);
          break;
        case FT::RHO:
          this->additionalAtomData["rho"].push_back(arg);
          break;
        case FT::DIAMETER:
          this->additionalAtomData["diameter"].push_back(arg);
          break;
        case FT::DENSITY:
          this->additionalAtomData["density"].push_back(arg);
          break;
        case FT::EDPD_TEMP:
          this->additionalAtomData["edpd_temp"].push_back(arg);
          break;
        case FT::EDPD_CV:
          this->additionalAtomData["edpd_cv"].push_back(arg);
          break;
        case FT::ESPH:
          this->additionalAtomData["esph"].push_back(arg);
          break;
        case FT::CV:
          this->additionalAtomData["cv"].push_back(arg);
          break;
        case FT::ENERGY:
          this->additionalAtomData["energy"].push_back(arg);
          break;
        case FT::BODYFLAG:
          this->additionalAtomData["bodyflag"].push_back(arg);
          break;
        case FT::ELLIPSOIDFLAG:
          this->additionalAtomData["ellipsoidflag"].push_back(arg);
          break;
        case FT::LINEFLAG:
          this->additionalAtomData["lineflag"].push_back(arg);
          break;
        case FT::TRIANGLEFLAG:
          this->additionalAtomData["triangleflag"].push_back(arg);
          break;
        case FT::STATUS:
          this->additionalAtomData["status"].push_back(arg);
          break;
        case FT::KRADIUS:
          this->additionalAtomData["kradius"].push_back(arg);
          break;
        case FT::CRADIUS:
          this->additionalAtomData["cradius"].push_back(arg);
          break;
        case FT::ERADIUS:
          this->additionalAtomData["eradius"].push_back(arg);
          break;
        case FT::MASS:
          this->additionalAtomData["mass"].push_back(arg);
          break;
        case FT::VOLUME:
          this->additionalAtomData["volume"].push_back(arg);
          break;
        case FT::X0:
          this->additionalAtomData["x0"].push_back(arg);
          break;
        case FT::Y0:
          this->additionalAtomData["y0"].push_back(arg);
          break;
        case FT::Z0:
          this->additionalAtomData["z0"].push_back(arg);
          break;
        case FT::SPX:
          this->additionalAtomData["spx"].push_back(arg);
          break;
        case FT::SPY:
          this->additionalAtomData["spy"].push_back(arg);
          break;
        case FT::SPZ:
          this->additionalAtomData["spz"].push_back(arg);
          break;
        case FT::SP:
          this->additionalAtomData["sp"].push_back(arg);
          break;
        case FT::ESPIN:
          this->additionalAtomData["espin"].push_back(arg);
          break;
        case FT::AREA:
          this->additionalAtomData["area"].push_back(arg);
          break;
        case FT::ED:
          this->additionalAtomData["ed"].push_back(arg);
          break;
        case FT::EM:
          this->additionalAtomData["em"].push_back(arg);
          break;
        case FT::EPSILON:
          this->additionalAtomData["epsilon"].push_back(arg);
          break;
        case FT::CURVATURE:
          this->additionalAtomData["curvature"].push_back(arg);
          break;
        case FT::TEMPLATE_INDEX:
          this->additionalAtomData["template_index"].push_back(arg);
          break;
        case FT::TEMPLATE_ATOM:
          this->additionalAtomData["template_atom"].push_back(arg);
          break;
        case FT::ETAG:
          this->additionalAtomData["etag"].push_back(arg);
          break;
        case FT::CS_RE:
          this->additionalAtomData["cs_re"].push_back(arg);
          break;
        case FT::CS_IM:
          this->additionalAtomData["cs_im"].push_back(arg);
          break;
        case FT::CC1:
          this->additionalAtomData["cc1"].push_back(arg);
          break;
        case FT::CC2:
          this->additionalAtomData["cc2"].push_back(arg);
          break;
        default:
          throw std::invalid_argument(
            "Unknown field type in AtomFieldDescriptor: " +
            std::to_string(static_cast<int>(descriptor.fields[argIdx])));
          break;
      }
    }
    argIdx++;
  };

  (processArg(args), ...);

  // Set default molecule ID for styles that don't have it
  if (std::find(descriptor.fields.begin(),
                descriptor.fields.end(),
                AtomFieldDescriptor::FieldType::MOLECULE_ID) ==
      descriptor.fields.end()) {
    this->moleculeIds.push_back(0);
  }
}

void
DataFileParser::readAtom(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::MOLECULAR);
  size_t atomId, nx, ny, nz;
  int atomType, moleculeId;
  double x, y, z;

  readAtomGeneric(
    line, descriptor, atomId, moleculeId, atomType, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomFull(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::FULL);
  size_t atomId, nx, ny, nz;
  int atomType, moleculeId;
  double charge, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  moleculeId,
                  atomType,
                  charge,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomCharge(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::CHARGE);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double charge, x, y, z;

  readAtomGeneric(
    line, descriptor, atomId, atomType, charge, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomDipole(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::DIPOLE);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double charge, x, y, z, dipoleX, dipoleY, dipoleZ;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  charge,
                  x,
                  y,
                  z,
                  dipoleX,
                  dipoleY,
                  dipoleZ,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomDpd(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::DPD);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double x, y, z, theta;

  readAtomGeneric(
    line, descriptor, atomId, atomType, theta, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomMdpd(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::MDPD);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double x, y, z, rho;

  readAtomGeneric(line, descriptor, atomId, atomType, rho, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomSphere(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::SPHERE);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double x, y, z, diameter, density;

  readAtomGeneric(
    line, descriptor, atomId, atomType, diameter, density, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomAtomic(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::ATOMIC);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double x, y, z;

  readAtomGeneric(line, descriptor, atomId, atomType, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomBody(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::BODY);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, bodyflag;
  double mass, x, y, z;

  readAtomGeneric(
    line, descriptor, atomId, atomType, bodyflag, mass, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomBpmSphere(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::BPM_SPHERE);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, moleculeId;
  double diameter, density, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  moleculeId,
                  atomType,
                  diameter,
                  density,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomDielectric(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::DIELECTRIC);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double charge, x, y, z, mux, muy, muz, area, ed, em, epsilon, curvature;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  charge,
                  x,
                  y,
                  z,
                  mux,
                  muy,
                  muz,
                  area,
                  ed,
                  em,
                  epsilon,
                  curvature,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomEdpd(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::EDPD);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double edpd_temp, edpd_cv, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  edpd_temp,
                  edpd_cv,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomElectron(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::ELECTRON);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double charge, espin, eradius, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  charge,
                  espin,
                  eradius,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomEllipsoid(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::ELLIPSOID);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, ellipsoidflag;
  double density, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  ellipsoidflag,
                  density,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomLine(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::LINE);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, moleculeId, lineflag;
  double density, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  moleculeId,
                  atomType,
                  lineflag,
                  density,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomPeri(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::PERI);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double volume, density, x, y, z;

  readAtomGeneric(
    line, descriptor, atomId, atomType, volume, density, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomRheo(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::RHEO);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, status;
  double rho, x, y, z;

  readAtomGeneric(
    line, descriptor, atomId, atomType, status, rho, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomRheoThermal(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::RHEO_THERMAL);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, status;
  double rho, energy, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  status,
                  rho,
                  energy,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomSmd(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::SMD);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, moleculeId;
  double volume, mass, kradius, cradius, x0, y0, z0, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  moleculeId,
                  volume,
                  mass,
                  kradius,
                  cradius,
                  x0,
                  y0,
                  z0,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomSph(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::SPH);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double rho, esph, cv, x, y, z;

  readAtomGeneric(
    line, descriptor, atomId, atomType, rho, esph, cv, x, y, z, nx, ny, nz);
}

void
DataFileParser::readAtomSpin(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::SPIN);
  size_t atomId;
  long int nx, ny, nz;
  int atomType;
  double x, y, z, spx, spy, spz, sp;

  readAtomGeneric(
    line, descriptor, atomId, atomType, x, y, z, spx, spy, spz, sp, nx, ny, nz);
}

void
DataFileParser::readAtomTemplate(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::TEMPLATE);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, moleculeId, templateIndex, templateAtom;
  double x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  moleculeId,
                  templateIndex,
                  templateAtom,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomTri(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::TRI);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, moleculeId, triangleflag;
  double density, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  moleculeId,
                  atomType,
                  triangleflag,
                  density,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

void
DataFileParser::readAtomWavepacket(const std::string& line)
{
  auto descriptor = getAtomStyleDescriptor(AtomStyle::WAVEPACKET);
  size_t atomId;
  long int nx, ny, nz;
  int atomType, etag;
  double charge, espin, eradius, cs_re, cs_im, x, y, z;

  readAtomGeneric(line,
                  descriptor,
                  atomId,
                  atomType,
                  charge,
                  espin,
                  eradius,
                  etag,
                  cs_re,
                  cs_im,
                  x,
                  y,
                  z,
                  nx,
                  ny,
                  nz);
}

// Keep readAtomHybrid as-is since it's truly unique
void
DataFileParser::readAtomHybrid(const std::string& line,
                               const AtomStyle style1,
                               const AtomStyle style2)
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

  RUNTIME_EXP_IFN(resFound >= 8,
                  "Did not find enough data in line '" + line + "': only " +
                    std::to_string(resFound) + ".");

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

void
DataFileParser::readBonds(std::ifstream& file, std::string& line)
{
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
        "Data file ended too early. Not enough bonds read. Read " +
        std::to_string(i + 1) + " of " + std::to_string(this->nBonds) +
        " expected bonds.");
    }
  }
}

void
DataFileParser::readBond(const std::string& line)
{
  size_t bondId, bondType, newBondFrom, newBondTo;
  sscanf(line.c_str(),
         "%zu %zu %zu %zu",
         &bondId,
         &bondType,
         &newBondFrom,
         &newBondTo);
  this->bondIds.push_back(bondId);
  this->bondTypes.push_back(bondType);
  this->bondFrom.push_back(newBondFrom);
  this->bondTo.push_back(newBondTo);
}

void
DataFileParser::readAngles(std::ifstream& file, std::string& line)
{
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
        "Data file ended too early. Not enough angles read. Read " +
        std::to_string(i + 1) + " of " + std::to_string(this->nAngles) +
        " expected angles.");
    }
  }
}

void
DataFileParser::readAngle(const std::string& line)
{
  size_t newAngleId, newAngleType, newAngleFrom, newAngleVia, newAngleTo;
  sscanf(line.c_str(),
         "%zu %zu %zu %zu %zu",
         &newAngleId,
         &newAngleType,
         &newAngleFrom,
         &newAngleVia,
         &newAngleTo);

  this->angleIds.push_back(newAngleId);
  this->angleTypes.push_back(newAngleType);
  this->angleFrom.push_back(newAngleFrom);
  this->angleVia.push_back(newAngleVia);
  this->angleTo.push_back(newAngleTo);
}

void
DataFileParser::readDihedralAngles(std::ifstream& file, std::string& line)
{
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
      throw std::runtime_error("Data file ended too early. Not enough "
                               "dihedral angles read. Read " +
                               std::to_string(i + 1) + " of " +
                               std::to_string(this->nDihedralAngles) +
                               " expected angles.");
    }
  }
}

void
DataFileParser::readDihedralAngle(const std::string& line)
{
  size_t newAngleId, newAngleType, newAngleFrom, newAngleVia1, newAngleVia2,
    newAngleTo;
  sscanf(line.c_str(),
         "%zu %zu %zu %zu %zu %zu",
         &newAngleId,
         &newAngleType,
         &newAngleFrom,
         &newAngleVia1,
         &newAngleVia2,
         &newAngleTo);

  this->dihedralAngleIds.push_back(newAngleId);
  this->dihedralAngleTypes.push_back(newAngleType);
  this->dihedralAngleFrom.push_back(newAngleFrom);
  this->dihedralAngleVia1.push_back(newAngleVia1);
  this->dihedralAngleVia2.push_back(newAngleVia2);
  this->dihedralAngleTo.push_back(newAngleTo);
}

void
DataFileParser::readVelocities(std::ifstream& file, std::string& line)
{
  // skip this line too
  if (!getline(file, line)) {
    throw std::runtime_error("Data file ended too early. Not able to read "
                             "any advertised velocities.");
  }
  // then, skip empty lines
  this->skipEmptyLines(line, file);
  std::unordered_map<size_t, size_t> atomIdToIdx;
  std::vector<double> unorderedVx, unorderedVy, unorderedVz;
  unorderedVx.reserve(static_cast<size_t>(this->nAtoms));
  unorderedVy.reserve(static_cast<size_t>(this->nAtoms));
  unorderedVz.reserve(static_cast<size_t>(this->nAtoms));
  // start reading, assuming as many as atoms
  for (size_t i = 0; i < static_cast<size_t>(this->nAtoms); ++i) {
    size_t atomId;
    double vx, vy, vz;
    sscanf(line.c_str(), "%zu %le %le %le", &atomId, &vx, &vy, &vz);
    unorderedVx.push_back(vx);
    unorderedVy.push_back(vy);
    unorderedVz.push_back(vz);
    atomIdToIdx[atomId] = i;

    if (!getline(file, line) && i + 1 < this->nAtoms) {
      throw std::runtime_error("Data file ended too early. Not enough "
                               "velocities read. Read " +
                               std::to_string(i + 1) + " of " +
                               std::to_string(this->nAtoms) +
                               " expected velocities.");
    }
  }
  assert(unorderedVx.size() == this->nAtoms);

  // re-order for the other stuff
  std::vector<double> orderedVx, orderedVy, orderedVz;
  orderedVx.reserve(static_cast<size_t>(this->nAtoms));
  orderedVy.reserve(static_cast<size_t>(this->nAtoms));
  orderedVz.reserve(static_cast<size_t>(this->nAtoms));

  for (size_t i = 0; i < static_cast<size_t>(this->nAtoms); ++i) {
    size_t realIdx = atomIdToIdx.at(static_cast<size_t>(this->atomIds[i]));
    orderedVx.push_back(unorderedVx[realIdx]);
    orderedVy.push_back(unorderedVy[realIdx]);
    orderedVz.push_back(unorderedVz[realIdx]);
  }

  this->additionalAtomData["vx"] = orderedVx;
  this->additionalAtomData["vy"] = orderedVy;
  this->additionalAtomData["vz"] = orderedVz;
}
}
