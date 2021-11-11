#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "StringUtil.h"
#include <map>
#include <filesystem>
#include "DataFileParser.h"

namespace pylimer_tools
{
  namespace utils
  {

    void DataFileParser::read(const std::string filePath)
    {
      if (!std::filesystem::exists(filePath))
      {
        throw std::invalid_argument("File to read (" + filePath + ") does not exist.");
      }
      char *cline = NULL;

      size_t len = 0;

      FILE *fp = fopen(filePath.c_str(), "r");
      if (fp == NULL)
      {
        throw std::runtime_error("Failed to open data file to read.");
      }

      // read everything until "Masses"
      while ((getline(&cline, &len, fp)) != -1)
      {
        std::string line = pylimer_tools::utils::trimLineOmitComment(cline);
        // skip empty lines
        if (line.empty())
        {
          continue;
        }
        // read up until the masses
        if (line.find("Masses") != std::string::npos)
        {
          break;
        }
        // read the nr of data points to read afterwards
        this->readNs(line);
      }

      // reserve space
      // for atom data
      this->atomIds.reserve(this->nAtoms);
      this->atomTypes.reserve(this->nAtoms);
      this->moleculeIds.reserve(this->nAtoms);
      this->atomX.reserve(this->nAtoms);
      this->atomY.reserve(this->nAtoms);
      this->atomZ.reserve(this->nAtoms);
      this->atomNx.reserve(this->nAtoms);
      this->atomNy.reserve(this->nAtoms);
      this->atomNz.reserve(this->nAtoms);
      // and bond data
      this->bondIds.reserve(this->nBonds);
      this->bondTypes.reserve(this->nBonds);
      this->bondFrom.reserve(this->nBonds);
      this->bondTo.reserve(this->nBonds);

      // skip empty lines plus the line with "Masses"
      while ((getline(&cline, &len, fp)) != -1)
      {
        std::string line = pylimer_tools::utils::trimLineOmitComment(cline);

        // skip empty lines
        if (!line.empty())
        {
          break;
        }
      }

      // Then, read masses, up until the next section ("atoms")
      do
      {
        std::string line = pylimer_tools::utils::trimLineOmitComment(cline);

        // break at empty lines
        if (line.empty())
        {
          break;
        }
        // read masses until e.g. atoms section
        if (line.find("Atoms") != std::string::npos)
        {
          break;
        }
        // read the mass...
        this->readMass(line);
      } while ((getline(&cline, &len, fp)) != -1);

      this->skipLinesToContains(cline, &len, fp, "Atoms");
      // skip this line too
      if ((getline(&cline, &len, fp)) == -1)
      {
        throw std::runtime_error("Data file ended too early. Not able to read any atoms.");
      }
      // then, skip empty lines
      this->skipEmptyLines(cline, &len, fp);

      // Then, read atoms, up until the next section ("bonds")
      for (int i = 0; i < this->nAtoms; ++i)
      {
        this->readAtom(std::string(cline));

        if ((getline(&cline, &len, fp)) == -1)
        {
          throw std::runtime_error("Data file ended too early. Not enough atoms read.");
        }
      }

      // Then, read bonds
      this->skipLinesToContains(cline, &len, fp, "Bonds");
      // skip this line too
      if ((getline(&cline, &len, fp)) == -1)
      {
        throw std::runtime_error("Data file ended too early. Not able to read any bonds.");
      }
      // then, skip empty lines
      this->skipEmptyLines(cline, &len, fp);

      for (int i = 0; i < this->nBonds; i++)
      {
        this->readBond(std::string(cline));

        if ((getline(&cline, &len, fp)) == -1 && i + 1 < this->nBonds)
        {
          throw std::runtime_error("Data file ended too early. Not enough bonds read.");
        }
      }

      // we ignore angles etc. for now.

      fclose(fp);
      if (cline)
      {
        free(cline);
      }
    }

    void DataFileParser::skipLinesToContains(char *cline, size_t *len, FILE *fp, std::string upTo)
    {
      do
      {
        std::string line(cline);

        if (contains(&line, upTo))
        {
          break;
        }
      } while ((getline(&cline, len, fp)) != -1);
    }

    void DataFileParser::skipEmptyLines(char *cline, size_t *len, FILE *fp)
    {
      do
      {
        std::string line = pylimer_tools::utils::trimLineOmitComment(cline);

        // skip until empty lines
        if (!line.empty())
        {
          break;
        }
      } while ((getline(&cline, len, fp)) != -1);
    }

    void DataFileParser::readNs(const std::string line)
    {
      if (contains(&line, "atoms"))
      {
        this->nAtoms = (this->parseTypesInLine<int>(line, 1))[0];
      }
      else if (contains(&line, "bonds"))
      {
        this->nBonds = (this->parseTypesInLine<int>(line, 1))[0];
      }
      else if (contains(&line, "atom types"))
      {
        this->nAtomTypes = (this->parseTypesInLine<int>(line, 1))[0];
      }
      else if (contains(&line, "bond types"))
      {
        this->nBondTypes = (this->parseTypesInLine<int>(line, 1))[0];
      }
      else if (contains(&line, "xlo xhi"))
      {
        std::vector<double> parsedL = this->parseTypesInLine<double>(line, 2);
        this->Lx = parsedL[1] - parsedL[0];
      }
      else if (contains(&line, "ylo yhi"))
      {
        std::vector<double> parsedL = this->parseTypesInLine<double>(line, 2);
        this->Ly = parsedL[1] - parsedL[0];
      }
      else if (contains(&line, "zlo zhi"))
      {
        std::vector<double> parsedL = this->parseTypesInLine<double>(line, 2);
        this->Lz = parsedL[1] - parsedL[0];
      }
    }

    void DataFileParser::readMass(const std::string line)
    {
      int iteration = 0;
      int key = 0;
      pylimer_tools::utils::CsvTokenizer tokenizer(line);
      if (tokenizer.getLength() != 2)
      {
        throw std::runtime_error("Incorrect nr of fields tokenized when reading masses");
      }

      key = tokenizer.get<int>(0);
      // for now, we just override duplicate keys
      this->masses[key] = tokenizer.get<double>(1);
    }

    void DataFileParser::readAtom(std::string line)
    {
      pylimer_tools::utils::CsvTokenizer tokenizer(line);

      this->atomIds.push_back(tokenizer.get<int>(0));
      this->moleculeIds.push_back(tokenizer.get<int>(1));
      this->atomTypes.push_back(tokenizer.get<int>(2));
      this->atomX.push_back(tokenizer.get<double>(3));
      this->atomY.push_back(tokenizer.get<double>(4));
      this->atomZ.push_back(tokenizer.get<double>(5));
      // TODO: be more flexible towards
      if (tokenizer.getLength() > 5)
      {
        this->atomNx.push_back(tokenizer.get<int>(6));
        this->atomNy.push_back(tokenizer.get<int>(7));
        this->atomNz.push_back(tokenizer.get<int>(8));
      }
    }

    void DataFileParser::readBond(std::string line)
    {
      pylimer_tools::utils::CsvTokenizer tokenizer(line);
      this->bondIds.push_back(tokenizer.get<long int>(0));
      this->bondTypes.push_back(tokenizer.get<int>(1));
      this->bondFrom.push_back(tokenizer.get<long int>(2));
      this->bondTo.push_back(tokenizer.get<long int>(3));
    }
  }
}
