#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Universe.h"
#include "UniverseSequence.h"
#include "../utils/DataFileParser.h"
#include "../utils/DumpFileParser.h"
#include "../utils/VectorUtils.h"
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools
{
  namespace entities
  {
    // TODO: connectivity (& graphs) could be stored only once if they stay the same over a sequence.
    void UniverseSequence::initializeFromDumpFile(const std::string initialStructureDataFile, const std::string dumpFile)
    {
      this->reset();

      pylimer_tools::utils::DataFileParser dataFileParser = pylimer_tools::utils::DataFileParser();
      dataFileParser.read(initialStructureDataFile);

      pylimer_tools::utils::DumpFileParser dumpFileParser = pylimer_tools::utils::DumpFileParser();
      dumpFileParser.read(dumpFile);

      size_t nrOfTimesteps = dumpFileParser.getLength();

#pragma omp for
      for (size_t i = 0; i < nrOfTimesteps; ++i)
      {
        Universe newUniverse = Universe(0.0, 0.0, 0.0);
        newUniverse.setTimestep(dumpFileParser.getValuesForAt<long int>(i, "TIMESTEP", 0)[0]);
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
        bool isUnwrapped = false;
        if (!dumpFileParser.keyHasDirectionalColumn("ATOMS", "", ""))
        {
          if (dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "u"))
          {
            isUnwrapped = true;
            positionSuffix = "u";
          }
          else
          {
            xMultiplier = newUniverse.getBox().getLx();
            yMultiplier = newUniverse.getBox().getLy();
            zMultiplier = newUniverse.getBox().getLz();
            if (dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "su"))
            {
              positionSuffix = "su";
              isUnwrapped = true;
            }
            else if (dumpFileParser.keyHasDirectionalColumn("ATOMS", "", "s"))
            {
              positionSuffix = "s";
            }
            else
            {
              throw std::runtime_error("Did not find neither positional atom fields in atom data of dump file " + dumpFile + ".");
            }
          }
        }

        std::vector<double> positionsX = dumpFileParser.getValuesForAt<double>(i, "ATOMS", "x" + positionSuffix);
        std::vector<double> positionsY = dumpFileParser.getValuesForAt<double>(i, "ATOMS", "y" + positionSuffix);
        std::vector<double> positionsZ = dumpFileParser.getValuesForAt<double>(i, "ATOMS", "z" + positionSuffix);
        if (xMultiplier != 1.0 && yMultiplier != 1.0 && zMultiplier != 1.0)
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

        int nAtoms = 0;
        if (dumpFileParser.hasKey("NUMBER OF ATOMS"))
        {
          std::vector<int> nAtomVec = dumpFileParser.getValuesForAt<int>(i, "NUMBER OF ATOMS", 0);
          if (nAtomVec.size() > 0)
          {
            nAtoms = nAtomVec[0];
          }
        }

        if (nAtoms == 0)
        {
          nAtoms = dataFileParser.getNrOfAtoms();
        }

        if (dumpFileParser.keyHasDirectionalColumn("ATOMS", "i", ""))
        {
          nx = dumpFileParser.getValuesForAt<int>(i, "ATOMS", "ix");
          ny = dumpFileParser.getValuesForAt<int>(i, "ATOMS", "iy");
          nz = dumpFileParser.getValuesForAt<int>(i, "ATOMS", "iz");
        }
        else
        {
          nx = pylimer_tools::utils::initializeWithValue(nAtoms, 0); //dataFileParser.getAtomNx();
          ny = pylimer_tools::utils::initializeWithValue(nAtoms, 0); //dataFileParser.getAtomNy();
          nz = pylimer_tools::utils::initializeWithValue(nAtoms, 0); //dataFileParser.getAtomNz();
        }

        newUniverse.addAtoms(
            nAtoms,
            dumpFileParser.getValuesForAt<long int>(i, "ATOMS", "id"),
            dumpFileParser.getValuesForAt<int>(i, "ATOMS", "type"),
            positionsX, positionsY, positionsZ,
            nx, ny, nz);
        newUniverse.addBonds(dataFileParser.getNrOfBonds(), dataFileParser.getBondFrom(), dataFileParser.getBondTo());
        newUniverse.setMasses(dataFileParser.getMasses());
        this->universeCache.insert_or_assign(i, newUniverse);
      }

      this->length = dumpFileParser.getLength();
    };

    void UniverseSequence::initializeFromDataSequence(const std::vector<std::string> dataFiles)
    {
      this->reset();
      this->dataFiles = dataFiles;
      this->length = dataFiles.size();
    };

    Universe UniverseSequence::next()
    {
      return this->atIndex(this->index++);
    }

    Universe UniverseSequence::atIndex(int index)
    {
      if (this->universeCache.contains(index))
      {
        return this->universeCache.at(index);
      }
      if (index >= this->length)
      {
        throw std::invalid_argument("Index larger than nr. of universes.");
      }
      this->universeCache.insert_or_assign(index, this->readDataFile(this->dataFiles[index]));
      return this->universeCache.at(index);
    }

    Universe UniverseSequence::readDataFile(const std::string filePath)
    {
      pylimer_tools::utils::DataFileParser fileParser = pylimer_tools::utils::DataFileParser();
      fileParser.read(filePath);
      Universe universe = Universe(fileParser.getLx(), fileParser.getLy(), fileParser.getLz());
      universe.addAtoms(fileParser.getNrOfAtoms(), fileParser.getAtomIds(), fileParser.getAtomTypes(), fileParser.getAtomX(), fileParser.getAtomY(), fileParser.getAtomZ(), fileParser.getAtomNx(), fileParser.getAtomNy(), fileParser.getAtomNz());
      universe.addBonds(fileParser.getNrOfBonds(), fileParser.getBondFrom(), fileParser.getBondTo());
      universe.setMasses(fileParser.getMasses());
      return universe;
    }

    std::vector<Universe> UniverseSequence::getAll()
    {
      std::vector<Universe> results;
      results.reserve(this->getLength());
      for (size_t i = 0; i < this->getLength(); ++i)
      {
        results.push_back(this->atIndex(i));
      }
      return results;
    }

    void UniverseSequence::resetIterator()
    {
      this->index = 0;
    }

    int UniverseSequence::getLength()
    {
      return this->length;
    }

    void UniverseSequence::reset()
    {
      this->universeCache.clear();
      this->dataFiles.clear();
      this->length = 0;
      this->resetIterator();
    }
  }
}
