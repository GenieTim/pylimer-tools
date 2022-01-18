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
    this->moleculeIdxSubsequent = !includeSwap;
  }
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
    file << "\t " << 1 << " bond types\n";
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
    file << "Atoms\n\n";
    // TODO: support molecule idxs
    int moleculeIdx = 0;
    for (int i = 0; i < this->universe.getNrOfAtoms(); ++i) {
      pylimer_tools::entities::Atom atom = this->universe.getAtomByVertexIdx(i);
      file << "\t" << atom.getId() << "\t" << moleculeIdx << "\t"
           << atom.getType() << "\t" << atom.getX() << "\t" << atom.getY()
           << "\t" << atom.getZ() << "\t" << atom.getNX() << "\t"
           << atom.getNY() << "\t" << atom.getNZ() << "\n";
    }
    file << "\n";

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
  bool includeAngles = false;
  bool moleculeIdxSubsequent = false;
};
} // namespace utils
} // namespace pylimer_tools

#endif
