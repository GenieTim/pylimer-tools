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
      int uniqueAtomTypes = std::max(this->universe.countAtomTypes().size(),
                                     this->universe.getMasses().size());

      file.open(filePath);
      file << std::setprecision(std::numeric_limits<double>::digits10 + 1);

      // write header
      file << "LAMMPS file generated using pylimer_tools at "
           << std::put_time(&tm, "%Y/%m/%d %H-%M-%S") << ".\n\n";
      file << "\t " << this->universe.getNrOfAtoms() << " atoms\n";
      file << "\t " << this->universe.getNrOfBonds() << " bonds\n";
      file << "\t "
           << (this->includeAngles ? this->universe.getNrOfAngles() : 0)
           << " angles\n";
      file << "\t " << 0 << " dihedrals\n";
      file << "\t " << 0 << " impropers\n";
      file << "\n";
      file << "\t " << uniqueAtomTypes << " atom types\n";
      file << "\t " << 1 << " bond types\n"; // TODO: fix bond types overall
      file << "\t " << (this->includeAngles ? 1 : 0) << " angle types\n";
      file << "\t " << 0 << " dihedral types\n";
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
      std::map<std::string, std::vector<long int>> bonds =
        this->universe.getBonds();
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
        std::map<std::string, std::vector<long int>> angles =
          this->universe.getAngles();
        for (size_t i = 0; i < this->universe.getNrOfAngles(); ++i) {
          int angleType = 1; // TODO: support angle types?
          file << "\t" << i << "\t" << angleType << "\t"
               << (this->oldNewAtomIdMap[angles["angle_from"][i]]) << "\t"
               << (this->oldNewAtomIdMap[angles["angle_via"][i]]) << "\t"
               << (this->oldNewAtomIdMap[angles["angle_to"][i]]) << "\n";
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
    bool moleculeIdxSwappable = false;
    int crosslinkerType = 2;
    bool reindexAtoms = false;
    bool moveIntoBox = false;
    bool attemptImageReset = false;
    std::string customAtomFormat = "";
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
      int nx =
        this->attemptImageReset
          ? this->getImageFlagForCoordinate(atom.getUnwrappedX(&box),
                                            box.getLowX(),
                                            box.getHighX())
          : atom.getNX();
      int ny =
        this->attemptImageReset
          ? this->getImageFlagForCoordinate(atom.getUnwrappedY(&box),
                                            box.getLowY(),
                                            box.getHighY())
          : atom.getNY();
      int nz =
        this->attemptImageReset
          ? this->getImageFlagForCoordinate(atom.getUnwrappedZ(&box),
                                            box.getLowZ(),
                                            box.getHighZ())
          : atom.getNZ();
      double x = this->conditionallyMoveCoordinateIntoBox(
        atom.getUnwrappedX(&box),
        box.getLowX(),
        box.getHighX());
      double y = this->conditionallyMoveCoordinateIntoBox(
        atom.getUnwrappedY(&box),
        box.getLowY(),
        box.getHighY());
      double z = this->conditionallyMoveCoordinateIntoBox(
        atom.getUnwrappedZ(&box),
        box.getLowZ(),
        box.getHighZ());
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
        std::vector<pylimer_tools::entities::Atom> atoms =
          (this->moleculeIdxSwappable || this->attemptImageReset)
            ? molecule.getAtomsLinedUp()
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
