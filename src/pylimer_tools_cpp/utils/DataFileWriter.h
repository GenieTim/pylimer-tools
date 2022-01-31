#ifndef DATA_FILE_WRITER_H
#define DATA_FILE_WRITER_H

#include "../entities/Atom.h"
#include "../entities/Universe.h"
#include "StringUtils.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
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
    int uniqueAtomTypes = std::max(this->universe.countAtomTypes().size(),
                                   this->universe.getMasses().size());

    file.open(filePath);

    // write header
    file << "LAMMPS file generated using pylimer_tools.\n\n";
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
    file << "\t " << 0 << " " << this->universe.getBox().getLx()
         << " xlo xhi\n";
    file << "\t " << 0 << " " << this->universe.getBox().getLy()
         << " ylo yhi\n";
    file << "\t " << 0 << " " << this->universe.getBox().getLz()
         << " zlo zhi\n";
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
    for (int i = 0; i < this->universe.getNrOfBonds(); ++i) {
      long int bondType = bonds["bond_type"][i];
      if (bondType == -1) {
        bondType = 1;
      }
      file << "\t" << i << "\t" << bondType << "\t"
           << (this->universe.getAtomIdByIdx(bonds["bond_from"][i])) << "\t"
           << (this->universe.getAtomIdByIdx(bonds["bond_to"][i])) << "\n";
    }
    file << "\n";

    // write angles
    if (this->includeAngles && this->universe.getNrOfAngles() > 0) {
      file << "Angles\n\n";
      std::map<std::string, std::vector<long int>> angles =
          this->universe.getAngles();
      for (int i = 0; i < this->universe.getNrOfAngles(); ++i) {
        int angleType = 1; // TODO: support angle types?
        file << "\t" << i << "\t" << angleType << "\t"
             << (angles["angle_from"][i]) << "\t" << (angles["angle_via"][i])
             << "\t" << (angles["angle_to"][i]) << "\n";
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
  void writeAtom(std::ofstream &file, pylimer_tools::entities::Atom atom,
                 int moleculeIdx, int nAtomsOutput) {
    long int atomId = this->reindexAtoms ? nAtomsOutput : atom.getId();
    file << "\t" << atom.getId() << "\t" << moleculeIdx << "\t"
         << atom.getType() << "\t" << atom.getX() << "\t" << atom.getY() << "\t"
         << atom.getZ() << "\t" << atom.getNX() << "\t" << atom.getNY() << "\t"
         << atom.getNZ() << "\n";
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
      for (int i = 0; i < atoms.size(); ++i) {
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
