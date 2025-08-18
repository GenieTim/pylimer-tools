#include "UniverseSequence.h"
#include "../io/DataFileParser.h"
#include "../io/DumpFileParser.h"
#include "../utils/LammpsAtomStyle.h"
#include "../utils/VectorUtils.h"
#include "Universe.h"
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::entities {
// TODO: connectivity (& graphs) could be stored only once if they stay the
// same over a sequence.
void
UniverseSequence::initializeFromDumpFile(
  const std::string& initialStructureDataFile,
  const std::string& dumpFile)
{
  this->modeDataFiles = false;
  this->reset();

  this->dataFileParser = pylimer_tools::utils::DataFileParser();
  dataFileParser.read(initialStructureDataFile);

  this->dumpFileParser = pylimer_tools::utils::DumpFileParser(dumpFile);
  // this->dumpFileParser.read();

  size_t nrOfTimeSteps = dumpFileParser.getLength();
  this->universeCache.reserve(nrOfTimeSteps);

  this->length = dumpFileParser.getLength();
};

void
UniverseSequence::initializeFromDataSequence(
  const std::vector<std::string>& newDataFiles)
{
  this->modeDataFiles = true;
  this->reset();
  this->dataFiles = newDataFiles;
  this->length = newDataFiles.size();
  this->universeCache.reserve(this->length);
};

Universe
UniverseSequence::next()
{
  return this->atIndex(this->index++);
}

Universe
UniverseSequence::atIndex(size_t idx)
{
  INVALIDARG_EXP_IFN(idx < this->length,
                     "Index (" + std::to_string(idx) +
                       ") larger than nr. of universes (" +
                       std::to_string(this->length) + ").");
  if (!pylimer_tools::utils::map_has_key(this->universeCache, idx)) {
    this->universeCache.emplace(idx,
                                this->modeDataFiles
                                  ? this->readDataFileAtIndex(idx)
                                  : this->readDumpFileAtIndex(idx));
  }
  return this->universeCache.at(idx);
}

void
UniverseSequence::setDataFileAtomStyle(
  const std::vector<pylimer_tools::utils::AtomStyle>& newDataFileAtomStyle)
{
  INVALIDARG_EXP_IFN(newDataFileAtomStyle.size() <= 3,
                     "Expect at most 3 atom styles");
  this->dataFileAtomStyle = newDataFileAtomStyle;
}

Universe
UniverseSequence::readDataFileAtIndex(const size_t index) const
{
  return this->readDataFile(this->dataFiles[index]);
}

Universe
UniverseSequence::readDumpFileAtIndex(size_t index)
{
  // std::cout << "Reading dump file at idx " << index << std::endl;
  this->dumpFileParser.readGroupByIdx(index);
  Universe newUniverse = Universe(1.0, 1.0, 1.0);
  std::vector<long int> timeStepData =
    this->dumpFileParser.getValuesForAt<long int>(index, "TIMESTEP", 0);
  if (timeStepData.size() == 0) {
    throw std::runtime_error(
      "Universe with index " + std::to_string(index) +
      " does not have enough data on its timestep. Is the file defect?");
  }
  newUniverse.setTimestep(timeStepData[0]);
  if (this->dumpFileParser.hasKey("BOX BOUNDS")) {
    std::vector<double> lo =
      this->dumpFileParser.getValuesForAt<double>(index, "BOX BOUNDS", 0);
    std::vector<double> hi =
      this->dumpFileParser.getValuesForAt<double>(index, "BOX BOUNDS", 1);
    if (lo.size() < 3 || hi.size() < 3) {
      throw std::runtime_error(
        "Universe with index " + std::to_string(index) +
        " does not have enough data on its box size, " +
        std::to_string(lo.size()) + " and " + std::to_string(hi.size()) +
        " instead of at least 3 each. Is the file defect?");
    }
    pylimer_tools::entities::Box box =
      Box(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]);
    newUniverse.setBox(box);
  } else {
    newUniverse.setBoxLengths(this->dataFileParser.getLx(),
                              this->dataFileParser.getLy(),
                              this->dataFileParser.getLz());
  }

  std::string positionSuffix = "";
  bool isUnwrapped = false;
  bool isScaled = false;
  if (!this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "")) {
    if (this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "u")) {
      isUnwrapped = true;
      positionSuffix = "u";
    } else {
      if (this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "su")) {
        positionSuffix = "su";
        isUnwrapped = true;
        isScaled = true;
      } else if (this->dumpFileParser.keyHasDirectionalColumn(
                   "ATOMS", "", "s")) {
        positionSuffix = "s";
        isScaled = true;
      } else {
        throw std::runtime_error("Did not find neither positional atom "
                                 "fields in atom data of dump file.");
      }
    }
  }

  // Pre-allocate vectors to avoid stack issues
  std::vector<double> positionsX, positionsY, positionsZ;
  positionsX = this->dumpFileParser.getValuesForAt<double>(
    index, "ATOMS", "x" + positionSuffix);
  positionsY = this->dumpFileParser.getValuesForAt<double>(
    index, "ATOMS", "y" + positionSuffix);
  positionsZ = this->dumpFileParser.getValuesForAt<double>(
    index, "ATOMS", "z" + positionSuffix);
  if (positionsZ.size() != positionsY.size() ||
      positionsY.size() != positionsX.size()) {
    throw std::runtime_error(
      "Atom coordinates for universe with index " + std::to_string(index) +
      " do not have the same size (" + std::to_string(positionsX.size()) +
      ", " + std::to_string(positionsY.size()) + ", " +
      std::to_string(positionsZ.size()) + "). Is the file defect?");
  };
  if (isScaled) {
    for (size_t i = 0; i < positionsZ.size(); ++i) {
      positionsX[i] *= newUniverse.getBox().getLx();
      positionsY[i] *= newUniverse.getBox().getLx();
      positionsZ[i] *= newUniverse.getBox().getLx();
    }
  }

  std::vector<int> nx, ny, nz;

  int nAtoms = 0;
  if (this->dumpFileParser.hasKey("NUMBER OF ATOMS")) {
    std::vector<int> nAtomVec =
      this->dumpFileParser.getValuesForAt<int>(index, "NUMBER OF ATOMS", 0);
    if (nAtomVec.size() > 0) {
      nAtoms = nAtomVec[0];
    }
  } else {
    // std::cout << "Number of atoms not found in dumpfile" << std::endl;
    nAtoms = this->dataFileParser.getNrOfAtoms();
  }

  if (this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "i", "") &&
      !isUnwrapped) {
    nx = this->dumpFileParser.getValuesForAt<int>(index, "ATOMS", "ix");
    ny = this->dumpFileParser.getValuesForAt<int>(index, "ATOMS", "iy");
    nz = this->dumpFileParser.getValuesForAt<int>(index, "ATOMS", "iz");
  } else {
    nx = pylimer_tools::utils::initializeWithValue(nAtoms, 0);
    ny = pylimer_tools::utils::initializeWithValue(nAtoms, 0);
    nz = pylimer_tools::utils::initializeWithValue(nAtoms, 0);
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
      Universe dataFileUniverse = Universe(this->dataFileParser.getLx(),
                                           this->dataFileParser.getLy(),
                                           this->dataFileParser.getLz());
      dataFileUniverse.addAtoms(this->dataFileParser.getAtomIds(),
                                this->dataFileParser.getAtomTypes(),
                                this->dataFileParser.getAtomX(),
                                this->dataFileParser.getAtomY(),
                                this->dataFileParser.getAtomZ(),
                                this->dataFileParser.getAtomNx(),
                                this->dataFileParser.getAtomNy(),
                                this->dataFileParser.getAtomNz());
      for (long int j = 0; j < nAtoms; ++j) {
        // std::cout << "Inferring type from data file for " << j << " atom id
        // "
        // << atomIds[j] << std::endl;
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

  newUniverse.addAtoms(
    atomIds, atomTypes, positionsX, positionsY, positionsZ, nx, ny, nz);
  // ignore it if the bond atoms do not exist, as we want to be compatible for
  // dumps of only certain atom groups
  newUniverse.addBonds(this->dataFileParser.getNrOfBonds(),
                       this->dataFileParser.getBondFrom(),
                       this->dataFileParser.getBondTo(),
                       this->dataFileParser.getBondTypes(),
                       true,
                       false);
  newUniverse.setMasses(this->dataFileParser.getMasses());
  return newUniverse;
};

Universe
UniverseSequence::readDataFile(const std::string& filePath) const
{
  pylimer_tools::utils::AtomStyle style1 =
    this->dataFileAtomStyle.size() > 0 ? this->dataFileAtomStyle[0]
                                       : pylimer_tools::utils::AtomStyle::ANGLE;
  pylimer_tools::utils::AtomStyle style2 =
    this->dataFileAtomStyle.size() > 1 ? this->dataFileAtomStyle[1]
                                       : pylimer_tools::utils::AtomStyle::ANGLE;
  pylimer_tools::utils::AtomStyle style3 =
    this->dataFileAtomStyle.size() > 2 ? this->dataFileAtomStyle[2]
                                       : pylimer_tools::utils::AtomStyle::ANGLE;
  pylimer_tools::utils::DataFileParser fileParser =
    pylimer_tools::utils::DataFileParser();
  fileParser.read(filePath, style1, style2, style3);
  Universe universe = Universe(Box(fileParser.getLowX(),
                                   fileParser.getHighX(),
                                   fileParser.getLowY(),
                                   fileParser.getHighY(),
                                   fileParser.getLowZ(),
                                   fileParser.getHighZ()));
  universe.addAtoms(fileParser.getAtomIds(),
                    fileParser.getAtomTypes(),
                    fileParser.getAtomX(),
                    fileParser.getAtomY(),
                    fileParser.getAtomZ(),
                    fileParser.getAtomNx(),
                    fileParser.getAtomNy(),
                    fileParser.getAtomNz(),
                    fileParser.getAdditionalAtomData());
  universe.addBonds(fileParser.getNrOfBonds(),
                    fileParser.getBondFrom(),
                    fileParser.getBondTo(),
                    fileParser.getBondTypes(),
                    false,
                    false);
  universe.setMasses(fileParser.getMasses());
  if (fileParser.getNrOfAngles() > 0) {
    universe.addAngles(fileParser.getAngleFrom(),
                       fileParser.getAngleVia(),
                       fileParser.getAngleTo(),
                       fileParser.getAngleTypes());
  }
  if (fileParser.getNrOfDihedralAngles() > 0) {
    universe.addDihedralAngles(fileParser.getDihedralAngleFrom(),
                               fileParser.getDihedralAngleVia1(),
                               fileParser.getDihedralAngleVia2(),
                               fileParser.getDihedralAngleTo(),
                               fileParser.getDihedralAngleTypes());
  }

  return universe;
}

std::vector<Universe>
UniverseSequence::getAll()
{
  std::vector<Universe> results;
  results.reserve(this->getLength());
  for (size_t i = 0; i < this->getLength(); ++i) {
    results.push_back(this->atIndex(i));
  }
  return results;
}

void
UniverseSequence::forgetAtIndex(const size_t idx)
{
  if (!this->modeDataFiles) {
    this->dumpFileParser.forgetAt(idx);
  }
  if (pylimer_tools::utils::map_has_key(this->universeCache, idx)) {
    this->universeCache.erase(idx);
  }
}

std::vector<double>
UniverseSequence::computeDistanceFromToAtoms(
  const std::vector<long int>& atomIdsFrom,
  const std::vector<long int>& atomIdsTo,
  bool reduceMemory)
{
  if (this->modeDataFiles) {
    throw std::runtime_error("Datafiles R_ee not implemented yet.");
  }

  INVALIDARG_EXP_IFN(atomIdsFrom.size() == atomIdsTo.size(),
                     "Same size from and to is required.");

  pylimer_tools::utils::ReadDumpFileSections sections =
    this->dumpFileParser.readDumpFileSections(
      pylimer_tools::utils::ReadableDumpFileSections::TIMESTEP |
      pylimer_tools::utils::ReadableDumpFileSections::ATOM |
      pylimer_tools::utils::ReadableDumpFileSections::BOX);
  std::vector<long int> timeSteps = sections.timesteps;
  std::vector<Box> boxes = sections.boxes;
  std::vector<std::vector<Atom>> atoms = sections.atoms;

  RUNTIME_EXP_IFN(timeSteps.size() == boxes.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(boxes.size()) + " boxes.");
  RUNTIME_EXP_IFN(timeSteps.size() == atoms.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(atoms.size()) + " atoms.");
  RUNTIME_EXP_IFN(
    timeSteps.size() >= this->getLength(),
    "Dump file seems inconsistent: read " + std::to_string(timeSteps.size()) +
      " time-steps, but expected " + std::to_string(this->getLength()) + ".");
  if (timeSteps.size() > this->getLength()) {
    std::cerr << "WARNING: Dump file is either inconsistent or still being "
                 "written to. Continuing with "
              << std::to_string(this->getLength()) << " instead of "
              << std::to_string(timeSteps.size()) << " time-steps. "
              << std::endl;
  }

  // first, check that we start at the beginning
  size_t startingIndex = 0;
  for (size_t i = 1; i < this->getLength(); ++i) {
    if (timeSteps[i] < timeSteps[i - 1]) {
      startingIndex = i;
      std::cerr << "Correcting starting index due to time-step order to "
                << startingIndex << std::endl;
    }
  }

  // assemble all distances
  std::vector<double> results;
  results.reserve(atomIdsTo.size() * (this->getLength() - startingIndex + 1));
  for (size_t universe_idx = startingIndex; universe_idx < this->getLength();
       ++universe_idx) {
    std::unordered_map<long int, int> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[universe_idx].size());
    for (size_t j = 0; j < atoms[universe_idx].size(); ++j) {
      atomIdToAtomIndex[atoms[universe_idx][j].getId()] = j;
    }
    for (size_t j = 0; j < atomIdsFrom.size(); ++j) {
      Atom atomFrom = atoms[universe_idx][atomIdToAtomIndex.at(atomIdsFrom[j])];
      Atom atomTo = atoms[universe_idx][atomIdToAtomIndex.at(atomIdsTo[j])];
      results.push_back((atomTo.getUnwrappedCoordinates(boxes[universe_idx]) -
                         atomFrom.getUnwrappedCoordinates(boxes[universe_idx]))
                          .norm());
    }
  }

  return results;
}

std::vector<Eigen::Vector3d>
UniverseSequence::computeVectorFromToAtoms(
  const std::vector<long int>& atomIdsFrom,
  const std::vector<long int>& atomIdsTo,
  bool reduceMemory)
{
  if (this->modeDataFiles) {
    throw std::runtime_error("Datafiles R_ee not implemented yet.");
  }

  INVALIDARG_EXP_IFN(atomIdsFrom.size() == atomIdsTo.size(),
                     "Same size from and to is required.");

  pylimer_tools::utils::ReadDumpFileSections sections =
    this->dumpFileParser.readDumpFileSections(
      pylimer_tools::utils::ReadableDumpFileSections::TIMESTEP |
      pylimer_tools::utils::ReadableDumpFileSections::ATOM |
      pylimer_tools::utils::ReadableDumpFileSections::BOX);
  std::vector<long int> timeSteps = sections.timesteps;
  std::vector<Box> boxes = sections.boxes;
  std::vector<std::vector<Atom>> atoms = sections.atoms;

  RUNTIME_EXP_IFN(timeSteps.size() == boxes.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(boxes.size()) + " boxes.");
  RUNTIME_EXP_IFN(timeSteps.size() == atoms.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(atoms.size()) + " atoms.");
  RUNTIME_EXP_IFN(
    timeSteps.size() >= this->getLength(),
    "Dump file seems inconsistent: read " + std::to_string(timeSteps.size()) +
      " time-steps, but expected " + std::to_string(this->getLength()) + ".");
  if (timeSteps.size() > this->getLength()) {
    std::cerr << "WARNING: Dump file is either inconsistent or still being "
                 "written to. Continuing with "
              << std::to_string(this->getLength()) << " instead of "
              << std::to_string(timeSteps.size()) << " time-steps. "
              << std::endl;
  }

  // first, check that we start at the beginning
  size_t startingIndex = 0;
  for (size_t i = 1; i < this->getLength(); ++i) {
    if (timeSteps[i] < timeSteps[i - 1]) {
      startingIndex = i;
      std::cerr << "Correcting starting index due to time-step order to "
                << startingIndex << std::endl;
    }
  }

  // assemble all distances
  std::vector<Eigen::Vector3d> results;
  results.reserve(atomIdsTo.size() * (this->getLength() - startingIndex + 1));
  for (size_t universe_idx = startingIndex; universe_idx < this->getLength();
       ++universe_idx) {
    std::unordered_map<long int, int> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[universe_idx].size());
    for (size_t j = 0; j < atoms[universe_idx].size(); ++j) {
      atomIdToAtomIndex[atoms[universe_idx][j].getId()] = j;
    }
    for (size_t j = 0; j < atomIdsFrom.size(); ++j) {
      Atom atomFrom = atoms[universe_idx][atomIdToAtomIndex.at(atomIdsFrom[j])];
      Atom atomTo = atoms[universe_idx][atomIdToAtomIndex.at(atomIdsTo[j])];
      results.push_back(
        (atomTo.getUnwrappedCoordinates(boxes[universe_idx]) -
         atomFrom.getUnwrappedCoordinates(boxes[universe_idx])));
    }
  }

  return results;
}

std::unordered_map<long int, double>
UniverseSequence::computeDistanceAutocorrelationFromToAtoms(
  const std::vector<long int>& atomIdsFrom,
  const std::vector<long int>& atomIdsTo,
  int nrOfOrigins,
  bool reduceMemory)
{
  if (this->modeDataFiles) {
    throw std::runtime_error("Datafiles R_ee not implemented yet.");
  }

  INVALIDARG_EXP_IFN(atomIdsFrom.size() == atomIdsTo.size(),
                     "Same size from and to is required.");

  pylimer_tools::utils::ReadDumpFileSections sections =
    this->dumpFileParser.readDumpFileSections(
      pylimer_tools::utils::ReadableDumpFileSections::TIMESTEP |
      pylimer_tools::utils::ReadableDumpFileSections::ATOM |
      pylimer_tools::utils::ReadableDumpFileSections::BOX);
  std::vector<long int> timeSteps = sections.timesteps;
  std::vector<Box> boxes = sections.boxes;
  std::vector<std::vector<Atom>> atoms = sections.atoms;

  RUNTIME_EXP_IFN(timeSteps.size() == boxes.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(boxes.size()) + " boxes.");
  RUNTIME_EXP_IFN(timeSteps.size() == atoms.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(atoms.size()) + " atoms.");
  RUNTIME_EXP_IFN(
    timeSteps.size() >= this->getLength(),
    "Dump file seems inconsistent: read " + std::to_string(timeSteps.size()) +
      " time-steps, but expected " + std::to_string(this->getLength()) + ".");
  if (timeSteps.size() > this->getLength()) {
    std::cerr << "WARNING: Dump file is either inconsistent or still being "
                 "written to. Continuing with "
              << std::to_string(this->getLength()) << " instead of "
              << std::to_string(timeSteps.size()) << " time-steps. "
              << std::endl;
  }

  // first, check that we start at the beginning
  size_t startingIndex = 0;
  for (size_t i = 1; i < this->getLength(); ++i) {
    if (timeSteps[i] < timeSteps[i - 1]) {
      startingIndex = i;
      std::cerr << "Correcting starting index due to time-step order to "
                << startingIndex << std::endl;
    }
  }

  // assemble all coordinates
  std::vector<Eigen::VectorXd> endToEndVectors;
  endToEndVectors.reserve(this->getLength() - startingIndex);
  for (size_t i = startingIndex; i < this->getLength(); ++i) {
    std::unordered_map<long int, int> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[i].size());
    for (size_t j = 0; j < atoms[i].size(); ++j) {
      atomIdToAtomIndex[atoms[i][j].getId()] = j;
    }
    Eigen::VectorXd localCoordinatesFrom =
      Eigen::VectorXd::Zero(3 * atomIdsFrom.size());
    Eigen::VectorXd localCoordinatesTo =
      Eigen::VectorXd::Zero(3 * atomIdsTo.size());
    for (size_t j = 0; j < atomIdsFrom.size(); ++j) {
      Atom atomFrom = atoms[i][atomIdToAtomIndex.at(atomIdsFrom[j])];
      localCoordinatesFrom.segment(3 * j, 3) =
        atomFrom.getUnwrappedCoordinates(boxes[i]);
      Atom atomTo = atoms[i][atomIdToAtomIndex.at(atomIdsTo[j])];
      localCoordinatesTo.segment(3 * j, 3) =
        atomTo.getUnwrappedCoordinates(boxes[i]);
    }
    endToEndVectors.push_back(localCoordinatesTo - localCoordinatesFrom);
  }

  std::cout << "Assembled end-to-end vectors" << std::endl;

  std::unordered_map<long int, std::vector<double>> results;
  results.reserve(this->getLength() - startingIndex);
  // next, we actually start computations
  // this is a highly inefficient algorithm, but no idea how to do better
  // (except for omitting some data, skipping the graph, or other minor
  // optimizations)
  const int stepSize =
    std::max(1,
             static_cast<int>(
               std::floor((this->getLength() - startingIndex) / nrOfOrigins)));
  for (size_t parent_universe_idx = startingIndex;
       parent_universe_idx < this->getLength();
       parent_universe_idx += stepSize) {

    for (size_t universe_idx = parent_universe_idx;
         universe_idx < this->getLength();
         ++universe_idx) {
      long int delta_t =
        (timeSteps[universe_idx] - timeSteps[parent_universe_idx]);
      if (delta_t < 0) {
        std::cerr << "Encountered a negative delta time-step for universes "
                  << universe_idx << " (" << timeSteps[universe_idx] << ")"
                  << " and " << parent_universe_idx << " ("
                  << timeSteps[parent_universe_idx] << ") in file "
                  << this->dumpFileParser.getFile() << std::endl;
      }

      // verify
      assert(endToEndVectors[universe_idx - startingIndex].size() ==
             atomIdsFrom.size() * 3);
      assert(endToEndVectors[parent_universe_idx - startingIndex].size() ==
             atomIdsFrom.size() * 3);
      // compute the mean for all atoms between these two universes
      double localMean =
        (endToEndVectors[universe_idx - startingIndex].dot(
          endToEndVectors[parent_universe_idx - startingIndex])) /
        (static_cast<double>(atomIdsFrom.size()));
      results[delta_t].push_back(localMean);
    }
    std::cout << "Universe " << parent_universe_idx
              << " as basis has been handled." << std::endl;
  }

  // actually compute the mean
  std::unordered_map<long int, double> actual_means;
  actual_means.reserve(results.size());
  for (const auto& result_pair : results) {
    std::vector<double> sds = result_pair.second;
    if (sds.size() == 0) {
      continue;
    }
    actual_means[result_pair.first] =
      (Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sds.data(), sds.size()))
        .mean();
  }

  return actual_means;
}

// computations
/**
 * @brief Compute the mean square displacement for the given atoms
 *
 * @param atomIds
 * @param nrOfOrigins
 * @param reduceMemory
 * @return std::unordered_map<int, std::vector<double>>
 */
std::unordered_map<long int, double>
UniverseSequence::computeMsdForAtoms(const std::vector<long int>& atomIds,
                                     const int nrOfOrigins,
                                     const bool reduceMemory,
                                     const int max_tau)
{
  if (this->modeDataFiles) {
    return this->computeMsdForAtomsFromDataFiles(
      atomIds, nrOfOrigins, reduceMemory, max_tau);
  } else {
    return this->computeMsdForAtomsFromDumpFile(
      atomIds, nrOfOrigins, reduceMemory, max_tau);
  }
}

/**
 * @brief Compute the mean square displacement for the given atoms in the dump
 * file
 *
 * @param atomIds
 * @param x
 * @param y
 * @param z
 * @param nrOfOrigins
 * @param reduceMemory
 * @return std::unordered_map<long int, double>
 */
std::unordered_map<long int, double>
UniverseSequence::computeMsdForAtomProperties(
  const std::vector<long int>& atomIds,
  std::string x,
  std::string y,
  std::string z,
  int nrOfOrigins,
  bool reduceMemory,
  int max_tau)
{
  RUNTIME_EXP_IFN(!this->modeDataFiles,
                  "UniverseSequence::computeMsdForAtomsFromDumpFile only "
                  "works with dump files, duh.");

  pylimer_tools::utils::ReadDumpFileSections sections =
    this->dumpFileParser.readDumpFileSections(
      pylimer_tools::utils::ReadableDumpFileSections::TIMESTEP |
      pylimer_tools::utils::ReadableDumpFileSections::ATOM |
      pylimer_tools::utils::ReadableDumpFileSections::EXTRA_ATOM |
      pylimer_tools::utils::ReadableDumpFileSections::BOX);
  std::vector<long int> timeSteps = sections.timesteps;
  std::vector<Box> boxes = sections.boxes;
  std::vector<std::vector<Atom>> atoms = sections.atoms;
  std::vector<std::unordered_map<std::string, std::vector<double>>>
    extraAtomData = sections.extraAtomsData;

  RUNTIME_EXP_IFN(timeSteps.size() == boxes.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(boxes.size()) + " boxes.");
  RUNTIME_EXP_IFN(timeSteps.size() == atoms.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(atoms.size()) + " atoms.");
  RUNTIME_EXP_IFN(extraAtomData.size() == atoms.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(extraAtomData.size()) +
                    " extra atom sets, but " + std::to_string(atoms.size()) +
                    " atoms.");
  RUNTIME_EXP_IFN(
    timeSteps.size() >= this->getLength(),
    "Dump file seems inconsistent: read " + std::to_string(timeSteps.size()) +
      " time-steps, but expected " + std::to_string(this->getLength()) + ".");
  if (timeSteps.size() > this->getLength()) {
    std::cerr << "WARNING: Dump file is either inconsistent or still being "
                 "written to. Continuing with "
              << std::to_string(this->getLength()) << " instead of "
              << std::to_string(timeSteps.size()) << " time-steps. "
              << std::endl;
  }

  // first, check that we start at the beginning
  size_t startingIndex = 0;
  for (size_t i = 1; i < this->getLength(); ++i) {
    if (timeSteps[i] < timeSteps[i - 1]) {
      startingIndex = i;
      std::cerr << "Correcting starting index due to time-step order to "
                << startingIndex << std::endl;
    }
  }

  // assemble all coordinates
  std::vector<Eigen::VectorXd> coordinates;
  coordinates.reserve(this->getLength() - startingIndex);
  for (size_t i = startingIndex; i < this->getLength(); ++i) {
    std::unordered_map<long int, int> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[i].size());
    for (size_t j = 0; j < atoms[i].size(); ++j) {
      atomIdToAtomIndex[atoms[i][j].getId()] = j;
    }
    Eigen::VectorXd localCoordinates =
      Eigen::VectorXd::Zero(3 * atomIds.size());
    RUNTIME_EXP_IFN(extraAtomData[i].at(x).size() == atoms[i].size(),
                    "Wrong size");
    RUNTIME_EXP_IFN(extraAtomData[i].at(y).size() == atoms[i].size(),
                    "Wrong size");
    RUNTIME_EXP_IFN(extraAtomData[i].at(z).size() == atoms[i].size(),
                    "Wrong size");
    for (size_t j = 0; j < atomIds.size(); ++j) {
      size_t row = atomIdToAtomIndex.at(atomIds[j]);
      Eigen::Vector3d coords;
      coords << extraAtomData[i].at(x)[row], extraAtomData[i].at(y)[row],
        extraAtomData[i].at(z)[row];
      localCoordinates.segment(3 * j, 3) = coords;
    }
    coordinates.push_back(localCoordinates);
  }

  std::cout << "Assembled coordinates" << std::endl;

  std::unordered_map<long int, std::vector<double>> results;
  results.reserve(this->getLength() - startingIndex);
  // next, we actually start computations
  // this is a highly inefficient algorithm, but no idea how to do better
  // (except for omitting some data, skipping the graph, or other minor
  // optimizations)
  const int stepSize =
    std::max(1,
             static_cast<int>(
               std::floor((this->getLength() - startingIndex) / nrOfOrigins)));
  Eigen::VectorXd distance = Eigen::VectorXd::Zero(3 * atomIds.size());
  for (size_t parent_universe_idx = startingIndex;
       parent_universe_idx < this->getLength();
       parent_universe_idx += stepSize) {

    for (size_t universe_idx = parent_universe_idx + 1;
         universe_idx < this->getLength();
         ++universe_idx) {
      long int delta_t =
        (timeSteps[universe_idx] - timeSteps[parent_universe_idx]);
      if (delta_t < 0) {
        std::cerr << "Encountered a negative delta time-step for universes "
                  << universe_idx << " (" << timeSteps[universe_idx] << ")"
                  << " and " << parent_universe_idx << " ("
                  << timeSteps[parent_universe_idx] << ") in file "
                  << this->dumpFileParser.getFile() << std::endl;
      }

      // Skip if delta_t exceeds max_tau for better statistics
      if (delta_t > max_tau && max_tau >= 0) {
        continue;
      }

      distance = coordinates[universe_idx - startingIndex] -
                 coordinates[parent_universe_idx - startingIndex];
      double localMean = 0.0;
      double localMeanDenominator = 1. / static_cast<double>(atomIds.size());
      for (size_t atom_id = 0; atom_id < atomIds.size(); ++atom_id) {
        localMean +=
          distance.segment(3 * atom_id, 3).squaredNorm() * localMeanDenominator;
      }
      results[delta_t].push_back(localMean);
    }
    std::cout << "Universe " << parent_universe_idx
              << " as basis has been handled." << std::endl;
  }

  // actually compute the mean
  std::unordered_map<long int, double> actual_means;
  actual_means.reserve(results.size());
  for (const auto& result_pair : results) {
    std::vector<double> sds = result_pair.second;
    if (sds.size() == 0) {
      continue;
    }
    actual_means[result_pair.first] =
      (Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sds.data(), sds.size()))
        .mean();
  }

  return actual_means;
}

/**
 * @brief Compute the mean square displacement for the given atoms in the dump
 * file
 *
 * @param atomIds
 * @param nrOfOrigins
 * @param reduceMemory
 * @return std::unordered_map<long int, double>
 */
std::unordered_map<long int, double>
UniverseSequence::computeMsdForAtomsFromDumpFile(
  const std::vector<long int>& atomIds,
  int nrOfOrigins,
  bool reduceMemory,
  int max_tau)
{
  RUNTIME_EXP_IFN(!this->modeDataFiles,
                  "UniverseSequence::computeMsdForAtomsFromDumpFile only "
                  "works with dump files, duh.");

  pylimer_tools::utils::ReadDumpFileSections sections =
    this->dumpFileParser.readDumpFileSections(
      pylimer_tools::utils::ReadableDumpFileSections::TIMESTEP |
      pylimer_tools::utils::ReadableDumpFileSections::ATOM |
      pylimer_tools::utils::ReadableDumpFileSections::BOX);
  std::vector<long int> timeSteps = sections.timesteps;
  std::vector<Box> boxes = sections.boxes;
  std::vector<std::vector<Atom>> atoms = sections.atoms;

  RUNTIME_EXP_IFN(timeSteps.size() == boxes.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(boxes.size()) + " boxes.");
  RUNTIME_EXP_IFN(timeSteps.size() == atoms.size(),
                  "Dump file seems inconsistent: read " +
                    std::to_string(timeSteps.size()) + " time-steps, but " +
                    std::to_string(atoms.size()) + " atoms.");
  RUNTIME_EXP_IFN(
    timeSteps.size() >= this->getLength(),
    "Dump file seems inconsistent: read " + std::to_string(timeSteps.size()) +
      " time-steps, but expected " + std::to_string(this->getLength()) + ".");
  if (timeSteps.size() > this->getLength()) {
    std::cerr << "WARNING: Dump file is either inconsistent or still being "
                 "written to. Continuing with "
              << std::to_string(this->getLength()) << " instead of "
              << std::to_string(timeSteps.size()) << " time-steps. "
              << std::endl;
  }

  // first, check that we start at the beginning
  size_t startingIndex = 0;
  for (size_t i = 1; i < this->getLength(); ++i) {
    if (timeSteps[i] < timeSteps[i - 1]) {
      startingIndex = i;
      std::cerr << "Correcting starting index due to time-step order to "
                << startingIndex << std::endl;
    }
  }

  // assemble all coordinates
  std::vector<Eigen::VectorXd> coordinates;
  coordinates.reserve(this->getLength() - startingIndex);
  for (size_t i = startingIndex; i < this->getLength(); ++i) {
    std::unordered_map<long int, int> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[i].size());
    for (size_t j = 0; j < atoms[i].size(); ++j) {
      atomIdToAtomIndex[atoms[i][j].getId()] = j;
    }
    Eigen::VectorXd localCoordinates =
      Eigen::VectorXd::Zero(3 * atomIds.size());
    for (size_t j = 0; j < atomIds.size(); ++j) {
      Atom atom = atoms[i][atomIdToAtomIndex.at(atomIds[j])];
      Eigen::Vector3d coords;
      atom.getUnwrappedCoordinates<Eigen::Vector3d>(coords, boxes[i]);
      localCoordinates.segment(3 * j, 3) = coords;
    }
    coordinates.push_back(localCoordinates);
  }

  std::cout << "Assembled coordinates" << std::endl;

  std::unordered_map<long int, std::vector<double>> results;
  results.reserve(this->getLength() - startingIndex);
  // next, we actually start computations
  // this is a highly inefficient algorithm, but no idea how to do better
  // (except for omitting some data, skipping the graph, or other minor
  // optimizations)
  const int stepSize =
    std::max(1,
             static_cast<int>(
               std::floor((this->getLength() - startingIndex) / nrOfOrigins)));
  Eigen::VectorXd distance = Eigen::VectorXd::Zero(3 * atomIds.size());
  for (size_t parent_universe_idx = startingIndex;
       parent_universe_idx < this->getLength();
       parent_universe_idx += stepSize) {

    for (size_t universe_idx = parent_universe_idx + 1;
         universe_idx < this->getLength();
         ++universe_idx) {
      long int delta_t =
        (timeSteps[universe_idx] - timeSteps[parent_universe_idx]);
      if (delta_t < 0) {
        std::cerr << "Encountered a negative delta time-step for universes "
                  << universe_idx << " (" << timeSteps[universe_idx] << ")"
                  << " and " << parent_universe_idx << " ("
                  << timeSteps[parent_universe_idx] << ") in file "
                  << this->dumpFileParser.getFile() << std::endl;
      }

      // Skip if delta_t exceeds max_tau for better statistics
      if (delta_t > max_tau && max_tau >= 0) {
        continue;
      }

      distance = coordinates[universe_idx - startingIndex] -
                 coordinates[parent_universe_idx - startingIndex];
      double localMean = 0.0;
      double localMeanDenominator = 1. / static_cast<double>(atomIds.size());
      for (size_t atom_id = 0; atom_id < atomIds.size(); ++atom_id) {
        localMean +=
          distance.segment(3 * atom_id, 3).squaredNorm() * localMeanDenominator;
      }
      results[delta_t].push_back(localMean);
    }
    std::cout << "Universe " << parent_universe_idx
              << " as basis has been handled." << std::endl;
  }

  // actually compute the mean
  std::unordered_map<long int, double> actual_means;
  actual_means.reserve(results.size());
  for (const auto& result_pair : results) {
    std::vector<double> sds = result_pair.second;
    if (sds.size() == 0) {
      continue;
    }
    actual_means[result_pair.first] =
      (Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sds.data(), sds.size()))
        .mean();
  }

  return actual_means;
}

/**
 * @brief
 *
 * @param atomIds
 * @param nrOfOrigins
 * @param reduceMemory
 * @return std::unordered_map<long int, double>
 */
std::unordered_map<long int, double>
UniverseSequence::computeMsdForAtomsFromDataFiles(
  const std::vector<long int>& atomIds,
  int nrOfOrigins,
  bool reduceMemory,
  int max_tau)
{
  igraph_vector_int_t vertex_ids;
  igraph_vector_int_init(&vertex_ids, atomIds.size());
  pylimer_tools::entities::Universe initialUniverse = this->atIndex(0);
  for (long int atomId : atomIds) {
    igraph_vector_int_push_back(&vertex_ids,
                                initialUniverse.getIdxByAtomId(atomId));
  }
  pylimer_tools::entities::Box box = initialUniverse.getBox();
  std::unordered_map<int, std::vector<double>> results;
  results.reserve(this->getLength());
  std::vector<Eigen::VectorXd> coordinates;
  coordinates.reserve(this->getLength());
  std::vector<int> timeSteps;
  timeSteps.reserve(this->getLength());

  // initial loop to reserve everything
  for (size_t universe_idx = 0; universe_idx < this->getLength();
       ++universe_idx) {
    pylimer_tools::entities::Universe current_universe =
      this->atIndex(universe_idx);
    Eigen::VectorXd current_coordinates =
      current_universe.getUnwrappedVertexCoordinates(vertex_ids, box);
    coordinates.push_back(current_coordinates);
    timeSteps.push_back(current_universe.getTimestep());
    int delta_t =
      (current_universe.getTimestep() - initialUniverse.getTimestep());
    std::vector<double> resultsVec;
    resultsVec.reserve(atomIds.size() * this->getLength() *
                       (this->getLength() - universe_idx));
    results[delta_t] = resultsVec;
    RUNTIME_EXP_IFN(current_universe.getBox() == initialUniverse.getBox(),
                    "Boxes must be the same for all universes.");

    if (reduceMemory) {
      this->forgetAtIndex(universe_idx);
    }
  }

  // next, we actually start computations
  // this is a highly inefficient algorithm, but no idea how to do better
  // (except for omitting some data, skipping the graph, or other minor
  // optimizations)
  const int stepSize =
    std::max(1, static_cast<int>(std::floor(this->getLength() / nrOfOrigins)));
  Eigen::VectorXd distance = Eigen::VectorXd::Zero(3 * atomIds.size());
  for (size_t parent_universe_idx = 0; parent_universe_idx < this->getLength();
       parent_universe_idx += stepSize) {

    for (size_t universe_idx = parent_universe_idx + 1;
         universe_idx < this->getLength();
         ++universe_idx) {

      int delta_t = (timeSteps[universe_idx] - timeSteps[parent_universe_idx]);

      // Skip if delta_t exceeds max_tau for better statistics
      if (delta_t > max_tau && max_tau >= 0) {
        continue;
      }

      distance = coordinates[parent_universe_idx] - coordinates[universe_idx];
      for (size_t atom_id = 0; atom_id < atomIds.size(); ++atom_id) {
        results[delta_t].push_back(
          distance.segment(3 * atom_id, 3).squaredNorm());
      }
    }
    std::cout << "Universe " << parent_universe_idx
              << " as basis has been handled." << std::endl;
  }

  igraph_vector_int_destroy(&vertex_ids);
  // actually compute the mean
  std::unordered_map<long int, double> actual_means;
  actual_means.reserve(this->getLength());
  for (const auto& result_pair : results) {
    std::vector<double> sds = result_pair.second;
    if (sds.size() == 0) {
      continue;
    }
    actual_means[result_pair.first] =
      (Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sds.data(), sds.size()))
        .mean();
  }

  return actual_means;
}

// resets
void
UniverseSequence::resetIterator()
{
  this->index = 0;
}

size_t
UniverseSequence::getLength() const
{
  return this->length;
}

void
UniverseSequence::reset()
{
  this->universeCache.clear();
  this->dataFiles.clear();
  this->length = 0;
  this->resetIterator();
}
}
