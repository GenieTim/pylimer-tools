#ifndef UNIVERSE_SEQ_H
#define UNIVERSE_SEQ_H

#include <string>
#include <vector>
#include "Universe.h"
#include <unordered_map>
#include "../utils/TYPES.h"
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
      Universe atIndex(INDEX_TYPE index);
      void resetIterator();
      int getLength();
      void forgetAtIndex(INDEX_TYPE index);
      std::vector<Universe> getAll();

    protected:
      int index = 0; // current index of the iterator
      int length = 0;
      bool isInitialized = false;
      bool modeDataFiles = false;
      std::unordered_map<INDEX_TYPE, Universe> universeCache;
      std::vector<std::string> dataFiles;
      pylimer_tools::utils::DataFileParser dataFileParser;
      pylimer_tools::utils::DumpFileParser dumpFileParser;

      void reset();
      Universe readDataFile(const std::string filePath);
      Universe readDataFileAtIndex(const INDEX_TYPE index);
      Universe readDumpFileAtIndex(const INDEX_TYPE index);
    };
  }
}

#endif
