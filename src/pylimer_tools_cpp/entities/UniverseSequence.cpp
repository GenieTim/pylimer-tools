#include "UniverseSequence.h"
#include "../utils/DataFileParser.h"
#include "../utils/DumpFileParser.h"
#include "../utils/VectorUtils.h"
#include "Universe.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools {
namespace entities {
// TODO: connectivity (& graphs) could be stored only once if they stay the same
// over a sequence.
void UniverseSequence::initializeFromDumpFile(
    const std::string initialStructureDataFile, const std::string dumpFile) {
  this->modeDataFiles = false;
  this->reset();

  this->dataFileParser = pylimer_tools::utils::DataFileParser();
  dataFileParser.read(initialStructureDataFile);

  this->dumpFileParser = pylimer_tools::utils::DumpFileParser(dumpFile);
  // this->dumpFileParser.read();

  size_t nrOfTimesteps = dumpFileParser.getLength();
  this->universeCache.reserve(nrOfTimesteps);

  this->length = dumpFileParser.getLength();
};

void UniverseSequence::initializeFromDataSequence(
    const std::vector<std::string> dataFiles) {
  this->modeDataFiles = true;
  this->reset();
  this->dataFiles = dataFiles;
  this->length = dataFiles.size();
  this->universeCache.reserve(this->length);
};

Universe UniverseSequence::next() { return this->atIndex(this->index++); }

Universe UniverseSequence::atIndex(INDEX_TYPE index) {
  if (index >= this->length) {
    throw std::invalid_argument("Index larger than nr. of universes.");
  }
  if (!this->universeCache.contains(index)) {
    this->universeCache.emplace(index, this->modeDataFiles
                                           ? this->readDataFileAtIndex(index)
                                           : this->readDumpFileAtIndex(index));
  }
  return this->universeCache.at(index);
}

Universe UniverseSequence::readDataFileAtIndex(INDEX_TYPE index) {
  return this->readDataFile(this->dataFiles[index]);
}

Universe UniverseSequence::readDumpFileAtIndex(INDEX_TYPE index) {
  this->dumpFileParser.readGroupByIdx(index);
  Universe newUniverse = Universe(0.0, 0.0, 0.0);
  newUniverse.setTimestep(
      this->dumpFileParser.getValuesForAt<long int>(index, "TIMESTEP", 0)[0]);
  if (this->dumpFileParser.hasKey("BOX BOUNDS")) {
    std::vector<double> lo =
        this->dumpFileParser.getValuesForAt<double>(index, "BOX BOUNDS", 0);
    std::vector<double> hi =
        this->dumpFileParser.getValuesForAt<double>(index, "BOX BOUNDS", 1);
    newUniverse.setBoxLengths(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
  } else {
    newUniverse.setBoxLengths(this->dataFileParser.getLx(),
                              this->dataFileParser.getLy(),
                              this->dataFileParser.getLz());
  }

  std::string positionSuffix = "";
  double xMultiplier = 1.0;
  double yMultiplier = 1.0;
  double zMultiplier = 1.0;
  bool isUnwrapped = false;
  if (!this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "")) {
    if (this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "u")) {
      isUnwrapped = true;
      positionSuffix = "u";
    } else {
      xMultiplier = newUniverse.getBox().getLx();
      yMultiplier = newUniverse.getBox().getLy();
      zMultiplier = newUniverse.getBox().getLz();
      if (this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "su")) {
        positionSuffix = "su";
        isUnwrapped = true;
      } else if (this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "",
                                                              "s")) {
        positionSuffix = "s";
      } else {
        throw std::runtime_error("Did not find neither positional atom "
                                 "fields in atom data of dump file.");
      }
    }
  }

  std::vector<double> positionsX = this->dumpFileParser.getValuesForAt<double>(
      index, "ATOMS", "x" + positionSuffix);
  std::vector<double> positionsY = this->dumpFileParser.getValuesForAt<double>(
      index, "ATOMS", "y" + positionSuffix);
  std::vector<double> positionsZ = this->dumpFileParser.getValuesForAt<double>(
      index, "ATOMS", "z" + positionSuffix);
  if (xMultiplier != 1.0 && yMultiplier != 1.0 && zMultiplier != 1.0) {
    for (size_t i = 0; i < positionsZ.size(); ++i) {
      positionsX[i] *= xMultiplier;
      positionsY[i] *= yMultiplier;
      positionsZ[i] *= zMultiplier;
    }
  }

  std::vector<int> nx;
  std::vector<int> ny;
  std::vector<int> nz;

  int nAtoms = 0;
  if (this->dumpFileParser.hasKey("NUMBER OF ATOMS")) {
    std::vector<int> nAtomVec =
        this->dumpFileParser.getValuesForAt<int>(index, "NUMBER OF ATOMS", 0);
    if (nAtomVec.size() > 0) {
      nAtoms = nAtomVec[0];
    }
  }

  if (nAtoms == 0) {
    nAtoms = this->dataFileParser.getNrOfAtoms();
  }

  if (this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "i", "")) {
    nx = this->dumpFileParser.getValuesForAt<int>(index, "ATOMS", "ix");
    ny = this->dumpFileParser.getValuesForAt<int>(index, "ATOMS", "iy");
    nz = this->dumpFileParser.getValuesForAt<int>(index, "ATOMS", "iz");
  } else {
    nx = pylimer_tools::utils::initializeWithValue(
        nAtoms, 0); // this->dataFileParser.getAtomNx();
    ny = pylimer_tools::utils::initializeWithValue(
        nAtoms, 0); // this->dataFileParser.getAtomNy();
    nz = pylimer_tools::utils::initializeWithValue(
        nAtoms, 0); // this->dataFileParser.getAtomNz();
  }

  // read/parse/process atom ids
  std::vector<long int> atomIds;
  bool hasAtomIds = false;
  if (this->dumpFileParser.keyHasColumn("ATOMS", "id")) {
    atomIds =
        this->dumpFileParser.getValuesForAt<long int>(index, "ATOMS", "id");
    hasAtomIds = true;
  } else {
    atomIds.reserve(nAtoms);
    for (long int j = 0; j < nAtoms; ++j) {
      atomIds.push_back(j);
    }
  }

  // read/parse/process atom types
  std::vector<int> atomTypes;
  atomTypes.reserve(nAtoms);
  if (this->dumpFileParser.keyHasColumn("ATOMS", "type")) {
    atomTypes =
        this->dumpFileParser.getValuesForAt<int>(index, "ATOMS", "type");
  } else {
    if (hasAtomIds) {
      // infer from data file
      Universe dataFileUniverse =
          Universe(this->dataFileParser.getLx(), this->dataFileParser.getLy(),
                   this->dataFileParser.getLz());
      dataFileUniverse.addAtoms(
          this->dataFileParser.getNrOfAtoms(),
          this->dataFileParser.getAtomIds(),
          this->dataFileParser.getAtomTypes(), this->dataFileParser.getAtomX(),
          this->dataFileParser.getAtomY(), this->dataFileParser.getAtomZ(),
          this->dataFileParser.getAtomNx(), this->dataFileParser.getAtomNy(),
          this->dataFileParser.getAtomNz());
      for (long int j = 0; j < nAtoms; ++j) {
        atomTypes.push_back(dataFileUniverse.getAtom(atomIds[j]).getType());
      }
    } else {
      for (long int j = 0; j < nAtoms; ++j) {
        atomTypes.push_back(-1);
      }
    }
  }

  // some checking
  if (atomTypes.size() != nAtoms || atomIds.size() != nAtoms ||
      nx.size() != nAtoms || ny.size() != nAtoms || nz.size() != nAtoms ||
      positionsX.size() != nAtoms || positionsY.size() != nAtoms ||
      positionsZ.size() != nAtoms) {
    throw std::runtime_error("Failed to read timestep " +
                             std::to_string(newUniverse.getTimestep()) +
                             " due to different nr of atom properties");
  }

  newUniverse.addAtoms(nAtoms, atomIds, atomTypes, positionsX, positionsY,
                       positionsZ, nx, ny, nz);
  // ignore it if the bond atoms do not exist, as we want to be compatible for
  // dumps of only certain atom groups
  newUniverse.addBonds(this->dataFileParser.getNrOfBonds(),
                       this->dataFileParser.getBondFrom(),
                       this->dataFileParser.getBondTo(),
                       this->dataFileParser.getBondTypes(), true);
  newUniverse.setMasses(this->dataFileParser.getMasses());
  return newUniverse;
};

Universe UniverseSequence::readDataFile(const std::string filePath) {
  pylimer_tools::utils::DataFileParser fileParser =
      pylimer_tools::utils::DataFileParser();
  fileParser.read(filePath);
  Universe universe =
      Universe(fileParser.getLx(), fileParser.getLy(), fileParser.getLz());
  universe.addAtoms(fileParser.getNrOfAtoms(), fileParser.getAtomIds(),
                    fileParser.getAtomTypes(), fileParser.getAtomX(),
                    fileParser.getAtomY(), fileParser.getAtomZ(),
                    fileParser.getAtomNx(), fileParser.getAtomNy(),
                    fileParser.getAtomNz());
  universe.addBonds(fileParser.getNrOfBonds(), fileParser.getBondFrom(),
                    fileParser.getBondTo(), fileParser.getBondTypes());
  universe.setMasses(fileParser.getMasses());
  if (fileParser.getNrOfAngles() > 0) {
    universe.addAngles(fileParser.getAngleFrom(), fileParser.getAngleVia(),
                       fileParser.getAngleTo());
  }
  return universe;
}

std::vector<Universe> UniverseSequence::getAll() {
  std::vector<Universe> results;
  results.reserve(this->getLength());
  for (size_t i = 0; i < this->getLength(); ++i) {
    results.push_back(this->atIndex(i));
  }
  return results;
}

void UniverseSequence::forgetAtIndex(INDEX_TYPE index) {
  if (!this->modeDataFiles) {
    this->dumpFileParser.forgetAt(index);
  }
  this->universeCache.erase(index);
}

void UniverseSequence::resetIterator() { this->index = 0; }

int UniverseSequence::getLength() { return this->length; }

void UniverseSequence::reset() {
  this->universeCache.clear();
  this->dataFiles.clear();
  this->length = 0;
  this->resetIterator();
}
} // namespace entities
} // namespace pylimer_tools
