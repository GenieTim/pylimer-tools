#ifndef UNIVERSE_SEQ_H
#define UNIVERSE_SEQ_H

#include <string>
#include <vector>
#include "Universe.h"

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
      Universe atIndex(int index);
      void resetIterator();
      int getLength();

    protected:
      int index = 0; // current index of the iterator
      int length = 0;
      bool isInitialized = false;
      std::map<int, Universe> universeCache;
      std::vector<std::string> dataFiles;

      void reset();
      Universe readDataFile(const std::string filePath);
    };
  }
}

#endif
