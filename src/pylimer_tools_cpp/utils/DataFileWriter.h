#ifndef DATA_FILE_WRITER_H
#define DATA_FILE_WRITER_H

#include "../entities/Atom.h"
#include "../entities/Universe.h"
#include "StringUtils.h"
#include <algorithm>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace pylimer_tools {
namespace utils {

class DataFileWriter {

public:
  DataFileWriter(const pylimer_tools::entities::Universe u) : universe(u) {
    // this->universe = u;
  }
  void setUniverseToWrite(const pylimer_tools::entities::Universe u) {
    this->universe = u;
  };
  void configIncludeAngles(const bool includeAngles) {
    this->includeAngles = includeAngles;
  }
  void configMoleculeIdxForSwap(const bool includeSwap) {
    this->moleculeIdxSwappable = includeSwap;
  }
  void configCrosslinkerType(const int crosslinkerType) {
    this->crosslinkerType = crosslinkerType;
  }
  void configReindexAtoms(const bool reindex) { this->reindexAtoms = reindex; }
  void writeToFile(const std::string filePath) {
    std::ofstream file;
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    int uniqueAtomTypes = std::max(this->universe.countAtomTypes().size(),
                                   this->universe.getMasses().size());

    file.open(filePath);

    // write header
    file << "LAMMPS file generated using pylimer_tools at "
         << std::put_time(&tm, "%Y/%m/%d %H-%M-%S") << ".\n\n";
    file << "\t " << this->universe.getNrOfAtoms() << " atoms\n";
    file << "\t " << this->universe.getNrOfBonds() << " bonds\n";
    file << "\t " << (this->includeAngles ? this->universe.getNrOfAngles() : 0)
         << " angles\n";
    file << "\t " << 0 << " dihedrals\n";
    file << "\t " << 0 << " impropers\n";
    file << "\n";
    file << "\t " << uniqueAtomTypes << " atom types\n";
    file << "\t " << 1 << " bond types\n"; // TODO: fix bond types overall
    file << "\t " << 1 << " angle types\n";
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
    for (const auto &massPair : masses) {
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
  // functions
  double moveCoordinateIntoBox(double coord, double boxLo, double boxHi) const {
    double boxL = (boxHi - boxLo);
    while (coord > boxHi && coord > boxLo) {
      coord -= boxL;
    }
    while (coord < boxLo && coord < boxHi) {
      coord += boxL;
    }
    return coord;
  }
  void writeAtom(std::ofstream &file, pylimer_tools::entities::Atom atom,
                 int moleculeIdx, int nAtomsOutput) {
    long int atomId = this->reindexAtoms ? nAtomsOutput : atom.getId();
    this->oldNewAtomIdMap[atom.getId()] = atomId;
    file << "\t" << atomId << "\t" << moleculeIdx << "\t" << atom.getType()
         << "\t"
         << this->moveCoordinateIntoBox(atom.getX(),
                                        this->universe.getBox().getLowX(),
                                        this->universe.getBox().getHighX())
         << "\t"
         << this->moveCoordinateIntoBox(atom.getY(),
                                        this->universe.getBox().getLowY(),
                                        this->universe.getBox().getHighY())
         << "\t"
         << this->moveCoordinateIntoBox(atom.getZ(),
                                        this->universe.getBox().getLowZ(),
                                        this->universe.getBox().getHighZ())
         << "\t" << atom.getNX() << "\t" << atom.getNY() << "\t" << atom.getNZ()
         << "\n";
  }
  void writeAtoms(std::ofstream &file) {
    file << "Atoms\n\n";

    this->oldNewAtomIdMap.reserve(this->universe.getNrOfAtoms());
    int nAtomsOutput = 0;

    // to support molecule idxs, we need to adjust the order of atoms output
    // first, we output the crosslinker beads
    std::vector<pylimer_tools::entities::Atom> crosslinkers =
        this->universe.getAtomsWithType(this->crosslinkerType);
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
          this->moleculeIdxSwappable ? molecule.getAtomsLinedUp()
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
