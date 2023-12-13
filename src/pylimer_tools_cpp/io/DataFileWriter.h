#ifndef DATA_FILE_WRITER_H
#define DATA_FILE_WRITER_H

#include "../entities/Atom.h"
#include "../entities/Universe.h"
#include "../utils/StringUtils.h"
#include <algorithm>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {

  class DataFileWriter
  {

  public:
    DataFileWriter(const pylimer_tools::entities::Universe u)
      : universe(u)
    {
      // this->universe = u;
    }
    void setUniverseToWrite(const pylimer_tools::entities::Universe u)
    {
      this->universe = u;
    }
    void configIncludeAngles(const bool includeAngles)
    {
      this->includeAngles = includeAngles;
    }
    void configIncludeDihedralAngles(const bool includeDihedralAngles)
    {
      this->includeDihedralAngles = includeDihedralAngles;
    }
    void configMoveIntoBox(const bool doMoveIntoBox = true)
    {
      this->moveIntoBox = doMoveIntoBox;
    }
    void configAttemptImageReset(const bool doImageReset = true)
    {
      this->attemptImageReset = doImageReset;
    }
    void setCustomAtomFormat(const std::string atomFormat)
    {
      this->customAtomFormat = atomFormat;
      this->customAtomFormatAdditionalProperties.clear();
      // find custom properties
      std::regex property_regex("(\\$[a-zA-Z0-9]+)");
      auto words_begin = std::sregex_iterator(
        atomFormat.begin(), atomFormat.end(), property_regex);
      auto words_end = std::sregex_iterator();
      for (std::sregex_iterator i = words_begin; i != words_end; ++i) {
        std::smatch match = *i;
        std::string match_str = match.str();
        match_str.erase(std::remove(match_str.begin(), match_str.end(), '$'),
                        match_str.end());
        if (match_str != "atomId" && match_str != "moleculeId" &&
            match_str != "atomType" && match_str != "nx" && match_str != "ny" &&
            match_str != "nz" && match_str != "x" && match_str != "y" &&
            match_str != "z") {
          this->customAtomFormatAdditionalProperties.push_back(match_str);
        }
      }
    }
    void configMoleculeIdxForSwap(const bool includeSwap)
    {
      this->moleculeIdxSwappable = includeSwap;
    }
    void configCrosslinkerType(const int crosslinkerType)
    {
      this->crosslinkerType = crosslinkerType;
    }
    void configReindexAtoms(const bool reindex = true)
    {
      this->reindexAtoms = reindex;
    }
    void writeToFile(const std::string filePath)
    {
      std::ofstream file;
      auto t = std::time(nullptr);
      auto tm = *std::localtime(&t);
      std::vector<int> allAtomTypes = this->universe.getAtomTypes();
      int nrOfAtomTypes =
        pylimer_tools::utils::max_element<int>(allAtomTypes, 1);

      std::map<std::string, std::vector<long int>> bonds =
        this->universe.getBonds();
      std::map<std::string, std::vector<long int>> angles =
        this->universe.getAngles();
      std::map<std::string, std::vector<long int>> dihedral_angles =
        this->universe.getDihedralAngles();

      long int nrOfAngleTypes =
        pylimer_tools::utils::max_element<long int>(angles["angle_type"], 0);
      if (nrOfAngleTypes < 1) {
        nrOfAngleTypes = 1;
      }
      if (!this->includeAngles) {
        nrOfAngleTypes = 0;
      }

      long int nrOfDihedralAngleTypes =
        pylimer_tools::utils::max_element<long int>(
          dihedral_angles["dihedral_angle_type"], 0);
      if (nrOfDihedralAngleTypes < 1) {
        nrOfDihedralAngleTypes = 1;
      }
      if (!this->includeDihedralAngles ||
          dihedral_angles["dihedral_angle_type"].size() == 0) {
        nrOfDihedralAngleTypes = 0;
      }

      long int nrOfBondTypes =
        pylimer_tools::utils::max_element<long int>(bonds["bond_type"], 0);
      if (nrOfBondTypes < 1) {
        nrOfBondTypes = 1;
      }
      if (bonds["bond_from"].size() == 0) {
        nrOfBondTypes = 0;
      }

      file.open(filePath, std::ios::out | std::ios::trunc);
      if (!file.is_open()) {
        throw std::invalid_argument("Failed to open '" + filePath +
                                    "' for writing.");
      }
      file << std::setprecision(std::numeric_limits<double>::digits10 + 1);

      // write header
      file << "LAMMPS file generated using pylimer_tools at "
           << std::put_time(&tm, "%Y/%m/%d %H-%M-%S") << ".\n\n";
      file << "\t " << this->universe.getNrOfAtoms() << " atoms\n";
      file << "\t " << this->universe.getNrOfBonds() << " bonds\n";
      file << "\t "
           << (this->includeAngles ? this->universe.getNrOfAngles() : 0)
           << " angles\n";
      file << "\t "
           << (this->includeDihedralAngles
                 ? this->universe.getNrOfDihedralAngles()
                 : 0)
           << " dihedrals\n";
      file << "\t " << 0 << " impropers\n";
      file << "\n";
      file << "\t " << nrOfAtomTypes << " atom types\n";
      file << "\t " << nrOfBondTypes << " bond types\n";
      file << "\t " << nrOfAngleTypes << " angle types\n";
      file << "\t " << nrOfDihedralAngleTypes << " dihedral types\n";
      file << "\t " << 0 << " improper types\n";
      file << "\n";
      file << "\t " << this->universe.getBox().getLowX() << " "
           << this->universe.getBox().getHighX() << " xlo xhi\n";
      file << "\t " << this->universe.getBox().getLowY() << " "
           << this->universe.getBox().getHighY() << " ylo yhi\n";
      file << "\t " << this->universe.getBox().getLowZ() << " "
           << this->universe.getBox().getHighZ() << " zlo zhi\n";
      file << "\n";

      // write masses
      file << "Masses\n\n";
      std::map<int, double> masses = this->universe.getMasses();
      for (const auto& massPair : masses) {
        file << "\t" << massPair.first << " " << massPair.second << "\n";
      }
      file << "\n";

      // write atoms
      this->writeAtoms(file);

      // write bonds
      file << "Bonds\n\n";
      for (size_t i = 0; i < this->universe.getNrOfBonds(); ++i) {
        long int bondType = bonds.at("bond_type")[i];
        if (bondType == -1) {
          bondType = 1;
        }
        file << "\t" << i << "\t" << bondType << "\t"
             << (this->oldNewAtomIdMap.at(bonds.at("bond_from")[i])) << "\t"
             << (this->oldNewAtomIdMap.at(bonds.at("bond_to")[i])) << "\n";
      }
      file << "\n";

      // write angles
      if (this->includeAngles && this->universe.getNrOfAngles() > 0) {
        file << "Angles\n\n";
        for (size_t i = 0; i < this->universe.getNrOfAngles(); ++i) {
          file << "\t" << i << "\t" << angles["angle_type"][i] << "\t"
               << (this->oldNewAtomIdMap[angles["angle_from"][i]]) << "\t"
               << (this->oldNewAtomIdMap[angles["angle_via"][i]]) << "\t"
               << (this->oldNewAtomIdMap[angles["angle_to"][i]]) << "\n";
        }
        file << "\n";
      }

      // write dihedarl angles
      if (this->includeDihedralAngles &&
          this->universe.getNrOfDihedralAngles() > 0) {
        file << "Dihedrals\n\n";
        for (size_t i = 0; i < this->universe.getNrOfDihedralAngles(); ++i) {
          file
            << "\t" << i << "\t" << dihedral_angles["dihedral_angle_type"][i]
            << "\t"
            << (this
                  ->oldNewAtomIdMap[dihedral_angles["dihedral_angle_from"][i]])
            << "\t"
            << (this
                  ->oldNewAtomIdMap[dihedral_angles["dihedral_angle_via1"][i]])
            << "\t"
            << (this
                  ->oldNewAtomIdMap[dihedral_angles["dihedral_angle_via2"][i]])
            << "\t"
            << (this->oldNewAtomIdMap[dihedral_angles["dihedral_angle_to"][i]])
            << "\n";
        }
        file << "\n";
      }

      file.close();
    };

  private:
    // properties
    pylimer_tools::entities::Universe universe;
    std::unordered_map<long int, int> oldNewAtomIdMap;
    bool includeAngles = true;
    bool includeDihedralAngles = true;
    bool moleculeIdxSwappable = false;
    int crosslinkerType = 2;
    bool reindexAtoms = false;
    bool moveIntoBox = false;
    bool attemptImageReset = false;
    std::string customAtomFormat = "";
    std::vector<std::string> customAtomFormatAdditionalProperties;
    // functions
    // TODO: move the following to the box
    int getImageFlagForCoordinate(double coord,
                                  double boxLo,
                                  double boxHi) const
    {
      assert(boxHi > boxLo);
      int imageFlag = 0;
      double L = (boxHi + boxLo);
      while (coord > boxHi) {
        coord -= L;
        imageFlag += 1;
      }
      while (coord < boxLo) {
        coord += L;
        imageFlag -= 1;
      }
      return imageFlag;
    }
    double conditionallyMoveCoordinateIntoBox(double coord,
                                              double boxLo,
                                              double boxHi) const
    {
      assert(boxHi > boxLo);
      if (this->moveIntoBox == false) {
        return coord;
      }
      double boxL = (boxHi - boxLo);
      while (coord > boxHi && coord > boxLo) {
        coord -= boxL;
      }
      while (coord < boxLo && coord < boxHi) {
        coord += boxL;
      }
      return coord;
    }
    void writeAtom(std::ofstream& file,
                   pylimer_tools::entities::Atom atom,
                   int moleculeIdx,
                   int nAtomsOutput)
    {
      long int atomId = this->reindexAtoms ? nAtomsOutput : atom.getId();
      const pylimer_tools::entities::Box box = this->universe.getBox();
      this->oldNewAtomIdMap[atom.getId()] = atomId;
      int nx = this->attemptImageReset
                 ? this->getImageFlagForCoordinate(
                     atom.getUnwrappedX(&box), box.getLowX(), box.getHighX())
                 : atom.getNX();
      int ny = this->attemptImageReset
                 ? this->getImageFlagForCoordinate(
                     atom.getUnwrappedY(&box), box.getLowY(), box.getHighY())
                 : atom.getNY();
      int nz = this->attemptImageReset
                 ? this->getImageFlagForCoordinate(
                     atom.getUnwrappedZ(&box), box.getLowZ(), box.getHighZ())
                 : atom.getNZ();
      double x = this->conditionallyMoveCoordinateIntoBox(
        atom.getUnwrappedX(&box), box.getLowX(), box.getHighX());
      double y = this->conditionallyMoveCoordinateIntoBox(
        atom.getUnwrappedY(&box), box.getLowY(), box.getHighY());
      double z = this->conditionallyMoveCoordinateIntoBox(
        atom.getUnwrappedZ(&box), box.getLowZ(), box.getHighZ());
      if (this->customAtomFormat.size() < 2) {
        file << "\t" << atomId << "\t" << moleculeIdx << "\t" << atom.getType()
             << "\t" << x << "\t" << y << "\t" << z << "\t" << nx << "\t" << ny
             << "\t" << nz << "\n";
      } else {
        std::string outputStr = this->customAtomFormat;
        outputStr = std::regex_replace(
          outputStr, std::regex("\\$atomId"), std::to_string(atomId));
        outputStr = std::regex_replace(
          outputStr, std::regex("\\$moleculeId"), std::to_string(moleculeIdx));
        outputStr = std::regex_replace(
          outputStr, std::regex("\\$atomType"), std::to_string(atom.getType()));
        outputStr = std::regex_replace(
          outputStr, std::regex("\\$nx"), std::to_string(nx));
        outputStr = std::regex_replace(
          outputStr, std::regex("\\$ny"), std::to_string(ny));
        outputStr = std::regex_replace(
          outputStr, std::regex("\\$nz"), std::to_string(nz));
        outputStr =
          std::regex_replace(outputStr, std::regex("\\$x"), std::to_string(x));
        outputStr =
          std::regex_replace(outputStr, std::regex("\\$y"), std::to_string(y));
        outputStr =
          std::regex_replace(outputStr, std::regex("\\$z"), std::to_string(z));
        for (std::string additionalProperty :
             this->customAtomFormatAdditionalProperties) {
          outputStr = std::regex_replace(
            outputStr,
            std::regex("\\$" + additionalProperty),
            std::to_string(this->universe.getPropertyValue<double>(
              additionalProperty.c_str(),
              this->universe.getIdxByAtomId(atom.getId()))));
        }

        file << outputStr << "\n";
      }
    }
    void writeAtoms(std::ofstream& file)
    {
      file << "Atoms\n\n";

      this->oldNewAtomIdMap.reserve(this->universe.getNrOfAtoms());
      int nAtomsOutput = 0;

      // to support molecule idxs, we need to adjust the order of atoms output
      // first, we output the crosslinker beads
      std::vector<pylimer_tools::entities::Atom> crosslinkers =
        this->universe.getAtomsOfType(this->crosslinkerType);
      for (pylimer_tools::entities::Atom crosslinker : crosslinkers) {
        nAtomsOutput += 1;
        this->writeAtom(file, crosslinker, 0, nAtomsOutput);
      }

      // then, we can output all others
      int nMoleculesOutput = 0;
      std::vector<pylimer_tools::entities::Molecule> molecules =
        this->universe.getMolecules(this->crosslinkerType);
      for (pylimer_tools::entities::Molecule molecule : molecules) {
        // image flag reset attempt might not be the best yet
        // could try to use ->getAssumedVertexCoordinates() for more reset
        // options
        std::vector<pylimer_tools::entities::Atom> atoms =
          (this->moleculeIdxSwappable || this->attemptImageReset)
            ? molecule.getAtomsLinedUp(this->crosslinkerType, true)
            : molecule.getAtoms();
        nMoleculesOutput += 1;

        for (size_t i = 0; i < atoms.size(); ++i) {
          pylimer_tools::entities::Atom atom = atoms[i];
          nAtomsOutput += 1;
          int ip1 = i + 1;
          int swappableMoleculeIdx =
            (i >= (atoms.size() * 0.5)) ? (atoms.size() - i) : ip1;
          int moleculeIdx = this->moleculeIdxSwappable ? swappableMoleculeIdx
                                                       : nMoleculesOutput;

          this->writeAtom(file, atom, moleculeIdx, nAtomsOutput);
        }
      }

      file << "\n";
    }
  };
} // namespace utils
} // namespace pylimer_tools

#endif
