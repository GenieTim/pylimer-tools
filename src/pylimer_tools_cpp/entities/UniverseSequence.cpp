#include "UniverseSequence.h"
#include "../io/DataFileParser.h"
#include "../io/DumpFileParser.h"
#include "../utils/LammpsAtomStyle.h"
#include "../utils/VectorUtils.h"
#include "Universe.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::entities {

namespace {
  size_t computeStepSize(const size_t nFrames, const int nrOfOrigins)
  {
    if (nrOfOrigins <= 0) {
      return 1;
    }
    const double rawStepSize = std::floor(static_cast<double>(nFrames) /
                                          static_cast<double>(nrOfOrigins));
    return std::max<size_t>(1, static_cast<size_t>(std::max(0.0, rawStepSize)));
  }

  std::unordered_map<long int, double> averageVectorMap(
    const std::unordered_map<long int, std::vector<double>>& results)
  {
    std::unordered_map<long int, double> actualMeans;
    actualMeans.reserve(results.size());
    for (const std::pair<const long int, std::vector<double>>& resultPair :
         results) {
      const std::vector<double>& values = resultPair.second;
      if (values.empty()) {
        continue;
      }
      actualMeans[resultPair.first] =
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(
          values.data(), static_cast<Eigen::Index>(values.size()))
          .mean();
    }
    return actualMeans;
  }

  std::unordered_map<long int, double> computePairwiseMsdMeans(
    const std::vector<Eigen::VectorXd>& coordinates,
    const std::vector<long int>& timeSteps,
    const size_t startingIndex,
    const size_t totalLength,
    const size_t numAtoms,
    const int nrOfOrigins,
    const int maxTau,
    const std::string& filePathForWarning = "")
  {
    if (totalLength <= startingIndex || numAtoms == 0) {
      return {};
    }
    std::unordered_map<long int, std::vector<double>> results;
    results.reserve(totalLength - startingIndex);

    const size_t stepSize =
      computeStepSize(totalLength - startingIndex, nrOfOrigins);
    const double invNumAtoms = 1.0 / static_cast<double>(numAtoms);
    for (size_t parentUniverseIdx = startingIndex;
         parentUniverseIdx < totalLength;
         parentUniverseIdx += stepSize) {
      for (size_t universeIdx = parentUniverseIdx + 1;
           universeIdx < totalLength;
           ++universeIdx) {
        const long int deltaT =
          timeSteps[universeIdx] - timeSteps[parentUniverseIdx];
        if (deltaT < 0 && !filePathForWarning.empty()) {
          std::cerr << "Encountered a negative delta time-step for universes "
                    << universeIdx << " (" << timeSteps[universeIdx] << ")"
                    << " and " << parentUniverseIdx << " ("
                    << timeSteps[parentUniverseIdx] << ") in file "
                    << filePathForWarning << std::endl;
        }

        // Skip if deltaT exceeds maxTau for better statistics
        if (deltaT > maxTau && maxTau >= 0) {
          continue;
        }

        const double localMean =
          (coordinates[universeIdx - startingIndex] -
           coordinates[parentUniverseIdx - startingIndex])
            .squaredNorm() *
          invNumAtoms;
        results[deltaT].push_back(localMean);
      }
      std::cout << "Universe " << parentUniverseIdx
                << " as basis has been handled." << std::endl;
    }

    return averageVectorMap(results);
  }
} // namespace
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
UniverseSequence::readDataFileAtIndex(const size_t universeIndex) const
{
  return this->readDataFile(this->dataFiles[universeIndex]);
}

Universe
UniverseSequence::readDumpFileAtIndex(size_t universeIndex)
{
  // std::cout << "Reading dump file at idx " << universeIndex << std::endl;
  this->dumpFileParser.readGroupByIdx(universeIndex);
  Universe newUniverse = Universe(1.0, 1.0, 1.0);
  std::vector<long int> timeStepData =
    this->dumpFileParser.getValuesForAt<long int>(universeIndex, "TIMESTEP", 0);
  if (timeStepData.size() == 0) {
    throw std::runtime_error(
      "Universe with index " + std::to_string(universeIndex) +
      " does not have enough data on its timestep. Is the file defect?");
  }
  newUniverse.setTimestep(timeStepData[0]);
  if (this->dumpFileParser.hasKey("BOX BOUNDS")) {
    std::vector<double> lo = this->dumpFileParser.getValuesForAt<double>(
      universeIndex, "BOX BOUNDS", 0);
    std::vector<double> hi = this->dumpFileParser.getValuesForAt<double>(
      universeIndex, "BOX BOUNDS", 1);
    if (lo.size() < 3 || hi.size() < 3) {
      throw std::runtime_error(
        "Universe with index " + std::to_string(universeIndex) +
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
    universeIndex, "ATOMS", "x" + positionSuffix);
  positionsY = this->dumpFileParser.getValuesForAt<double>(
    universeIndex, "ATOMS", "y" + positionSuffix);
  positionsZ = this->dumpFileParser.getValuesForAt<double>(
    universeIndex, "ATOMS", "z" + positionSuffix);
  if (positionsZ.size() != positionsY.size() ||
      positionsY.size() != positionsX.size()) {
    throw std::runtime_error(
      "Atom coordinates for universe with index " +
      std::to_string(universeIndex) + " do not have the same size (" +
      std::to_string(positionsX.size()) + ", " +
      std::to_string(positionsY.size()) + ", " +
      std::to_string(positionsZ.size()) + "). Is the file defect?");
  };
  if (isScaled) {
    const double lx = newUniverse.getBox().getLx();
    const double ly = newUniverse.getBox().getLy();
    const double lz = newUniverse.getBox().getLz();
    for (size_t i = 0; i < positionsZ.size(); ++i) {
      positionsX[i] *= lx;
      positionsY[i] *= ly;
      positionsZ[i] *= lz;
    }
  }

  std::vector<int> nx, ny, nz;

  int nAtoms = 0;
  if (this->dumpFileParser.hasKey("NUMBER OF ATOMS")) {
    std::vector<int> nAtomVec = this->dumpFileParser.getValuesForAt<int>(
      universeIndex, "NUMBER OF ATOMS", 0);
    if (nAtomVec.size() > 0) {
      nAtoms = nAtomVec[0];
    }
  } else {
    // std::cout << "Number of atoms not found in dumpfile" << std::endl;
    nAtoms = this->dataFileParser.getNrOfAtoms();
  }

  const size_t nAtomsUnsigned = static_cast<size_t>(nAtoms);
  if (this->dumpFileParser.keyHasDirectionalColumn("ATOMS", "i", "") &&
      !isUnwrapped) {
    nx = this->dumpFileParser.getValuesForAt<int>(universeIndex, "ATOMS", "ix");
    ny = this->dumpFileParser.getValuesForAt<int>(universeIndex, "ATOMS", "iy");
    nz = this->dumpFileParser.getValuesForAt<int>(universeIndex, "ATOMS", "iz");
  } else {
    nx = pylimer_tools::utils::initializeWithValue(nAtomsUnsigned, 0);
    ny = pylimer_tools::utils::initializeWithValue(nAtomsUnsigned, 0);
    nz = pylimer_tools::utils::initializeWithValue(nAtomsUnsigned, 0);
  }

  // read/parse/process atom ids
  std::vector<long int> atomIds;
  bool hasAtomIds = false;
  if (this->dumpFileParser.keyHasColumn("ATOMS", "id")) {
    atomIds = this->dumpFileParser.getValuesForAt<long int>(
      universeIndex, "ATOMS", "id");
    hasAtomIds = true;
  } else {
    atomIds.reserve(nAtomsUnsigned);
    for (long int j = 0; j < nAtoms; ++j) {
      atomIds.push_back(j);
    }
  }

  // read/parse/process atom types
  std::vector<int> atomTypes;
  atomTypes.reserve(nAtomsUnsigned);
  if (this->dumpFileParser.keyHasColumn("ATOMS", "type")) {
    atomTypes =
      this->dumpFileParser.getValuesForAt<int>(universeIndex, "ATOMS", "type");
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
      for (size_t j = 0; j < nAtomsUnsigned; ++j) {
        // std::cout << "Inferring type from data file for " << j << " atom id
        // "
        // << atomIds[j] << std::endl;
        atomTypes.push_back(
          dataFileUniverse.getAtom(static_cast<int>(atomIds[j])).getType());
      }
    } else {
      for (long int j = 0; j < nAtoms; ++j) {
        atomTypes.push_back(-1);
      }
    }
  }

  // some checking
  if (atomTypes.size() != nAtomsUnsigned || atomIds.size() != nAtomsUnsigned ||
      nx.size() != nAtomsUnsigned || ny.size() != nAtomsUnsigned ||
      nz.size() != nAtomsUnsigned || positionsX.size() != nAtomsUnsigned ||
      positionsY.size() != nAtomsUnsigned ||
      positionsZ.size() != nAtomsUnsigned) {
    throw std::runtime_error("Failed to read timestep " +
                             std::to_string(newUniverse.getTimestep()) +
                             " due to different nr of atom properties");
  }

  newUniverse.addAtoms(
    atomIds, atomTypes, positionsX, positionsY, positionsZ, nx, ny, nz);
  // ignore it if the bond atoms do not exist, as we want to be compatible for
  // dumps of only certain atom groups
  newUniverse.addBonds(static_cast<size_t>(this->dataFileParser.getNrOfBonds()),
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
  universe.addBonds(static_cast<size_t>(fileParser.getNrOfBonds()),
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
  size_t numFrames =
    this->getLength() > startingIndex ? (this->getLength() - startingIndex) : 0;
  size_t pairsPerFrame = atomIdsTo.size();
  std::vector<double> results(numFrames * pairsPerFrame);
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic) if (numFrames > 4)
#endif
  for (size_t f = 0; f < numFrames; ++f) {
    size_t universe_idx = startingIndex + f;
    std::unordered_map<long int, size_t> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[universe_idx].size());
    for (size_t j = 0; j < atoms[universe_idx].size(); ++j) {
      atomIdToAtomIndex[atoms[universe_idx][j].getId()] = j;
    }
    for (size_t j = 0; j < pairsPerFrame; ++j) {
      Atom atomFrom = atoms[universe_idx][atomIdToAtomIndex.at(atomIdsFrom[j])];
      Atom atomTo = atoms[universe_idx][atomIdToAtomIndex.at(atomIdsTo[j])];
      results[f * pairsPerFrame + j] =
        (atomTo.getUnwrappedCoordinates(boxes[universe_idx]) -
         atomFrom.getUnwrappedCoordinates(boxes[universe_idx]))
          .norm();
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
  size_t numFrames =
    this->getLength() > startingIndex ? (this->getLength() - startingIndex) : 0;
  size_t pairsPerFrame = atomIdsTo.size();
  std::vector<Eigen::Vector3d> results(numFrames * pairsPerFrame);
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic) if (numFrames > 4)
#endif
  for (size_t f = 0; f < numFrames; ++f) {
    size_t universe_idx = startingIndex + f;
    std::unordered_map<long int, size_t> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[universe_idx].size());
    for (size_t j = 0; j < atoms[universe_idx].size(); ++j) {
      atomIdToAtomIndex[atoms[universe_idx][j].getId()] = j;
    }
    for (size_t j = 0; j < pairsPerFrame; ++j) {
      Atom atomFrom = atoms[universe_idx][atomIdToAtomIndex.at(atomIdsFrom[j])];
      Atom atomTo = atoms[universe_idx][atomIdToAtomIndex.at(atomIdsTo[j])];
      results[f * pairsPerFrame + j] =
        (atomTo.getUnwrappedCoordinates(boxes[universe_idx]) -
         atomFrom.getUnwrappedCoordinates(boxes[universe_idx]));
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
    std::unordered_map<long int, size_t> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[i].size());
    for (size_t j = 0; j < atoms[i].size(); ++j) {
      atomIdToAtomIndex[atoms[i][j].getId()] = j;
    }
    Eigen::VectorXd localCoordinatesFrom =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * atomIdsFrom.size()));
    Eigen::VectorXd localCoordinatesTo =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * atomIdsTo.size()));
    for (size_t j = 0; j < atomIdsFrom.size(); ++j) {
      Atom atomFrom = atoms[i][atomIdToAtomIndex.at(atomIdsFrom[j])];
      localCoordinatesFrom.segment(static_cast<Eigen::Index>(3 * j), 3) =
        atomFrom.getUnwrappedCoordinates(boxes[i]);
      Atom atomTo = atoms[i][atomIdToAtomIndex.at(atomIdsTo[j])];
      localCoordinatesTo.segment(static_cast<Eigen::Index>(3 * j), 3) =
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
  const size_t stepSize =
    computeStepSize(this->getLength() - startingIndex, nrOfOrigins);
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
                  << this->dumpFileParser.getFilePath() << std::endl;
      }

      // verify
      assert(static_cast<size_t>(
               endToEndVectors[universe_idx - startingIndex].size()) ==
             atomIdsFrom.size() * 3);
      assert(static_cast<size_t>(
               endToEndVectors[parent_universe_idx - startingIndex].size()) ==
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

  return averageVectorMap(results);
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
  const std::vector<long int>& timeSteps = sections.timesteps;
  const std::vector<Box>& boxes = sections.boxes;
  const std::vector<std::vector<Atom>>& atoms = sections.atoms;
  const std::vector<std::unordered_map<std::string, std::vector<double>>>&
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
    std::unordered_map<long int, size_t> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[i].size());
    for (size_t j = 0; j < atoms[i].size(); ++j) {
      atomIdToAtomIndex[atoms[i][j].getId()] = j;
    }
    Eigen::VectorXd localCoordinates =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * atomIds.size()));
    const std::vector<double>& propX = extraAtomData[i].at(x);
    const std::vector<double>& propY = extraAtomData[i].at(y);
    const std::vector<double>& propZ = extraAtomData[i].at(z);
    RUNTIME_EXP_IFN(propX.size() == atoms[i].size(), "Wrong size");
    RUNTIME_EXP_IFN(propY.size() == atoms[i].size(), "Wrong size");
    RUNTIME_EXP_IFN(propZ.size() == atoms[i].size(), "Wrong size");
    for (size_t j = 0; j < atomIds.size(); ++j) {
      const size_t row = atomIdToAtomIndex.at(atomIds[j]);
      localCoordinates.segment<3>(static_cast<Eigen::Index>(3 * j)) =
        Eigen::Vector3d(propX[row], propY[row], propZ[row]);
    }
    coordinates.push_back(std::move(localCoordinates));
  }

  std::cout << "Assembled coordinates" << std::endl;

  return computePairwiseMsdMeans(coordinates,
                                 timeSteps,
                                 startingIndex,
                                 this->getLength(),
                                 atomIds.size(),
                                 nrOfOrigins,
                                 max_tau,
                                 this->dumpFileParser.getFilePath());
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
  const std::vector<long int>& timeSteps = sections.timesteps;
  const std::vector<Box>& boxes = sections.boxes;
  const std::vector<std::vector<Atom>>& atoms = sections.atoms;

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
    std::unordered_map<long int, size_t> atomIdToAtomIndex;
    atomIdToAtomIndex.reserve(atoms[i].size());
    for (size_t j = 0; j < atoms[i].size(); ++j) {
      atomIdToAtomIndex[atoms[i][j].getId()] = j;
    }
    Eigen::VectorXd localCoordinates =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * atomIds.size()));
    for (size_t j = 0; j < atomIds.size(); ++j) {
      const Atom& atom = atoms[i][atomIdToAtomIndex.at(atomIds[j])];
      Eigen::Vector3d coords;
      atom.getUnwrappedCoordinates<Eigen::Vector3d>(coords, boxes[i]);
      localCoordinates.segment<3>(static_cast<Eigen::Index>(3 * j)) = coords;
    }
    coordinates.push_back(std::move(localCoordinates));
  }

  std::cout << "Assembled coordinates" << std::endl;

  return computePairwiseMsdMeans(coordinates,
                                 timeSteps,
                                 startingIndex,
                                 this->getLength(),
                                 atomIds.size(),
                                 nrOfOrigins,
                                 max_tau,
                                 this->dumpFileParser.getFilePath());
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
  if (this->getLength() == 0 || atomIds.empty()) {
    return {};
  }

  igraph_vector_int_t vertexIds;
  igraph_vector_int_init(&vertexIds, static_cast<igraph_int_t>(atomIds.size()));
  const pylimer_tools::entities::Universe initialUniverse = this->atIndex(0);
  const pylimer_tools::entities::Box box = initialUniverse.getBox();
  std::vector<Eigen::VectorXd> coordinates;
  coordinates.reserve(this->getLength());
  std::vector<long int> timeSteps;
  timeSteps.reserve(this->getLength());

  try {
    for (size_t universeIdx = 0; universeIdx < this->getLength();
         ++universeIdx) {
      const pylimer_tools::entities::Universe currentUniverse =
        this->atIndex(universeIdx);
      RUNTIME_EXP_IFN(currentUniverse.getBox() == box,
                      "Boxes must be the same for all universes.");
      for (size_t i = 0; i < atomIds.size(); ++i) {
        igraph_vector_int_set(&vertexIds,
                              static_cast<igraph_int_t>(i),
                              currentUniverse.getIdxByAtomId(atomIds[i]));
      }
      coordinates.push_back(
        currentUniverse.getUnwrappedVertexCoordinates(vertexIds, box));
      timeSteps.push_back(currentUniverse.getTimestep());

      if (reduceMemory) {
        this->forgetAtIndex(universeIdx);
      }
    }
  } catch (...) {
    igraph_vector_int_destroy(&vertexIds);
    throw;
  }
  igraph_vector_int_destroy(&vertexIds);

  // LAMMPS data files do not store explicit timesteps (all default to 0);
  // fall back to sequential frame indices when all timesteps are identical.
  if (std::adjacent_find(timeSteps.begin(),
                         timeSteps.end(),
                         std::not_equal_to<long int>()) == timeSteps.end()) {
    for (size_t i = 0; i < timeSteps.size(); ++i) {
      timeSteps[i] = static_cast<long int>(i);
    }
  }

  return computePairwiseMsdMeans(coordinates,
                                 timeSteps,
                                 0,
                                 this->getLength(),
                                 atomIds.size(),
                                 nrOfOrigins,
                                 max_tau);
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
