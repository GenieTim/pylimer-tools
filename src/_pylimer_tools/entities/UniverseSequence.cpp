#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Universe.h"
#include "UniverseSequence.h"
#include "DataFileParser.h"
#include "DumpFileParser.h"
#include <boost/algorithm/string.hpp>

namespace pylimer_tools
{
  namespace entities
  {
    // TODO: connectivity (& graphs) could be stored only once if they stay the same over a sequence.
    class UniverseSequence
    {
    public:
      void initializeFromDumpFile(const std::string initialStructureDataFile, const std::string dumpFile)
      {
        this->reset();

        pylimer_tools::utils::DataFileParser dataFileParser = pylimer_tools::utils::DataFileParser();
        dataFileParser.read(initialStructureDataFile);

        pylimer_tools::utils::DumpFileParser dumpFileParser = pylimer_tools::utils::DumpFileParser();
        dumpFileParser.read(dumpFile);

        for (int i = 0; i < dumpFileParser.getLength(); ++i)
        {
          Universe newUniverse = Universe(0, 0, 0);
          if (dumpFileParser.hasKey("BOX BOUNDS"))
          {
            std::vector<double> lo = dumpFileParser.getValuesForAt<double>(i, "BOX BOUNDS", 0);
            std::vector<double> hi = dumpFileParser.getValuesForAt<double>(i, "BOX BOUNDS", 1);
            newUniverse.setBoxLengths(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
          }
          else
          {
            newUniverse.setBoxLengths(dataFileParser.getLx(), dataFileParser.getLy(), dataFileParser.getLz());
          }

          std::string positionSuffix = "";
          double xMultiplier = 1.0;
          double yMultiplier = 1.0;
          double zMultiplier = 1.0;
          if (!dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "") && dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "su"))
          {
            positionSuffix = "su";
            xMultiplier = newUniverse.getBox().getLx();
            yMultiplier = newUniverse.getBox().getLy();
            zMultiplier = newUniverse.getBox().getLz();
          }
          else
          {
            throw std::runtime_error("Did not find neighter 'x' nor 'xsu' fields in atom data of dump file");
          }

          std::vector<double> positionsX = dumpFileParser.getValuesForAt<double>(i, "ATOMS", "x" + positionSuffix);
          std::vector<double> positionsY = dumpFileParser.getValuesForAt<double>(i, "ATOMS", "y" + positionSuffix);
          std::vector<double> positionsZ = dumpFileParser.getValuesForAt<double>(i, "ATOMS", "z" + positionSuffix);
          if (xMultiplier != 1.0)
          {
            for (size_t i = 0; i < positionsZ.size(); ++i)
            {
              positionsX[i] *= xMultiplier;
              positionsY[i] *= yMultiplier;
              positionsZ[i] *= zMultiplier;
            }
          }

          std::vector<int> nx;
          std::vector<int> ny;
          std::vector<int> nz;

          if (dumpFileParser.keyHasDirectionalColumn("ATOMS", "n", ""))
          {
            nx = dumpFileParser.getValuesForAt<int>(i, "ATOMS", "nx");
            ny = dumpFileParser.getValuesForAt<int>(i, "ATOMS", "ny");
            nz = dumpFileParser.getValuesForAt<int>(i, "ATOMS", "nz");
          }
          else
          {
            nx = dataFileParser.getAtomNx();
            ny = dataFileParser.getAtomNy();
            nz = dataFileParser.getAtomNz();
          }

          newUniverse.addAtoms(
              dataFileParser.getNrOfAtoms(),
              dumpFileParser.getValuesForAt<long int>(i, "ATOMS", "id"),
              dumpFileParser.getValuesForAt<int>(i, "ATOMS", "type"),
              positionsX, positionsY, positionsZ,
              nx, ny, nz);
          newUniverse.addBonds(dataFileParser.getNrOfBonds(), dataFileParser.getBondFrom(), dataFileParser.getBondTo());
          this->universeCache[i] = newUniverse;
        }

        this->length = dumpFileParser.getLength();
      };

      void initializeFromDataSequence(const std::vector<std::string> dataFiles)
      {
        this->reset();
        this->dataFiles = dataFiles;
        this->length = dataFiles.size();
      };

      Universe next()
      {
        return this->atIndex(this->index++);
      }

      Universe atIndex(int index)
      {
        if (this->universeCache.contains(index))
        {
          return this->universeCache[index];
        }
        if (index >= this->length)
        {
          throw std::invalid_argument("Index larger than nr. of universes.");
        }
        this->universeCache[index] = this->readDataFile(this->dataFiles[index]);
        return this->universeCache[index];
      }

      Universe readDataFile(const std::string filePath)
      {
        pylimer_tools::utils::DataFileParser fileParser = pylimer_tools::utils::DataFileParser();
        fileParser.read(filePath);
        Universe universe = Universe(fileParser.getLx(), fileParser.getLy(), fileParser.getLz());
        universe.addAtoms(fileParser.getNrOfAtoms(), fileParser.getAtomIds(), fileParser.getAtomTypes(), fileParser.getAtomX(), fileParser.getAtomY(), fileParser.getAtomZ(), fileParser.getAtomNx(), fileParser.getAtomNy(), fileParser.getAtomNz());
        universe.addBonds(fileParser.getNrOfBonds(), fileParser.getBondFrom(), fileParser.getBondTo());
        return universe;
      }

      void resetIterator()
      {
        this->index = 0;
      }

      int getLength()
      {
        return this->length;
      }

    private:
      int index; // current index of the iterator
      int length = 0;
      bool isInitialized = false;
      std::map<int, Universe> universeCache;
      std::vector<std::string> dataFiles;

      void reset()
      {
        this->universeCache.clear();
        this->dataFiles.clear();
        this->length = 0;
        this->resetIterator();
      }
    };
  }
}
