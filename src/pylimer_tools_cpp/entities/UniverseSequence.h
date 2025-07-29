#ifndef UNIVERSE_SEQ_H
#define UNIVERSE_SEQ_H

#include "../io/DataFileParser.h"
#include "../io/DumpFileParser.h"
#include "Universe.h"
#include <string>
#include <unordered_map>
#include <vector>

namespace pylimer_tools::entities {
class UniverseSequence
{
public:
  void initializeFromDumpFile(const std::string& initialStructureFile,
                              const std::string& dumpFile);
  void initializeFromDataSequence(const std::vector<std::string>& dataFiles);
  Universe next();
  Universe atIndex(size_t index);
  void resetIterator();
  size_t getLength() const;
  void forgetAtIndex(size_t index);
  std::vector<Universe> getAll();

  void setDataFileAtomStyle(
    const std::vector<pylimer_tools::utils::AtomStyle>& dataFileAtomStyle);

  // computations
  std::unordered_map<long int, double> computeMsdForAtomProperties(
    const std::vector<long int>& atomIds,
    std::string x,
    std::string y,
    std::string z,
    int nrOfOrigins = 10,
    bool reduceMemory = false,
    int max_tau = -1);
  std::unordered_map<long int, double> computeMsdForAtoms(
    const std::vector<long int>& atomIds,
    int nrOfOrigins = 10,
    bool reduceMemory = false,
    int max_tau = -1);
  std::unordered_map<long int, double>
  computeDistanceAutocorrelationFromToAtoms(
    const std::vector<long int>& atomIdsFrom,
    const std::vector<long int>& atomIdsTo,
    int nrOfOrigins = 10,
    bool reduceMemory = false);
  std::vector<double> computeDistanceFromToAtoms(
    const std::vector<long int>& atomIdsFrom,
    const std::vector<long int>& atomIdsTo,
    bool reduceMemory = false);
  std::vector<Eigen::Vector3d> computeVectorFromToAtoms(
    const std::vector<long int>& atomIdsFrom,
    const std::vector<long int>& atomIdsTo,
    bool reduceMemory = false);

protected:
  size_t index = 0; // current index of the iterator
  size_t length = 0;
  bool isInitialized = false;
  bool modeDataFiles = false;
  std::unordered_map<size_t, Universe> universeCache;
  std::vector<std::string> dataFiles;
  pylimer_tools::utils::DataFileParser dataFileParser;
  pylimer_tools::utils::DumpFileParser dumpFileParser;
  std::vector<pylimer_tools::utils::AtomStyle> dataFileAtomStyle;

  void reset();
  Universe readDataFile(const std::string& filePath) const;
  Universe readDataFileAtIndex(const size_t index) const;
  Universe readDumpFileAtIndex(const size_t index);

  std::unordered_map<long int, double> computeMsdForAtomsFromDataFiles(
    const std::vector<long int>& atomIds,
    int nrOfOrigins = 10,
    bool reduceMemory = false,
    int max_tau = -1);
  std::unordered_map<long int, double> computeMsdForAtomsFromDumpFile(
    const std::vector<long int>& atomIds,
    int nrOfOrigins = 10,
    bool reduceMemory = false,
    int max_tau = -1);
};
}

#endif
