#ifndef UNIVERSE_SEQ_H
#define UNIVERSE_SEQ_H

#include <string>
#include <vector>
#include "Universe.h"
#include <unordered_map>
#include "../utils/DataFileParser.h"
#include "../utils/DumpFileParser.h"

namespace pylimer_tools
{
  namespace entities
  {
    class UniverseSequence
    {
    public:
      void initializeFromDumpFile(const std::string initialStructureFile, const std::string dumpFile);
      void initializeFromDataSequence(const std::vector<std::string> dataFiles);
      Universe next();
      Universe atIndex(size_t index);
      void resetIterator();
      size_t getLength() const;
      void forgetAtIndex(size_t index);
      std::vector<Universe> getAll();

    protected:
      size_t index = 0; // current index of the iterator
      size_t length = 0;
      bool isInitialized = false;
      bool modeDataFiles = false;
      std::unordered_map<size_t, Universe> universeCache;
      std::vector<std::string> dataFiles;
      pylimer_tools::utils::DataFileParser dataFileParser;
      pylimer_tools::utils::DumpFileParser dumpFileParser;

      void reset();
      Universe readDataFile(const std::string filePath);
      Universe readDataFileAtIndex(const size_t index);
      Universe readDumpFileAtIndex(const size_t index);
    };
  }
}

#endif
