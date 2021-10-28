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
    };
  }
}

#endif
