#ifndef DATA_FILE_PARSER_H
#define DATA_FILE_PARSER_H

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "StringUtil.h"
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <map>
#include <filesystem>

namespace pylimer_tools
{
  namespace utils
  {

    class DataFileParser
    {
    public:
      void read(const std::string filePath);

      // access atom data
      int getNrOfAtoms() { return this->nAtoms; }
      int getNrOfAtomTypes() { return this->nAtomTypes; }
      std::vector<long int> getAtomIds() { return this->atomIds; }
      std::vector<int> getMoleculeIds() { return this->moleculeIds; }
      std::vector<int> getAtomTypes() { return this->atomTypes; }
      std::vector<double> getAtomX() { return this->atomX; }
      std::vector<double> getAtomY() { return this->atomY; }
      std::vector<double> getAtomZ() { return this->atomZ; }
      std::vector<int> getAtomNx() { return this->atomNx; }
      std::vector<int> getAtomNy() { return this->atomNy; }
      std::vector<int> getAtomNz() { return this->atomNz; }
      std::map<int, double> getMasses() { return this->masses; }

      // access bond data
      int getNrOfBonds() { return this->nBonds; }
      int getNrOfBondTypes() { return this->nBondTypes; }
      std::vector<int> getBondTypes() { return this->bondTypes; }
      std::vector<long int> getBondFrom() { return this->bondFrom; }
      std::vector<long int> getBondTo() { return this->bondTo; }

      // get box info
      double getLx() { return this->Lx; }
      double getLy() { return this->Ly; }
      double getLz() { return this->Lz; }

    private:
      bool shortenLineToSkip(std::string *line);
      bool shortenLineToSkip(char *line);
      void readNs(const std::string line);
      void readMass(const std::string line);
      void readAtom(std::string line);
      void readBond(std::string line);
      void skipEmptyLines(char *cline, size_t *len, FILE *fp);
      void skipLinesToContains(char *cline, size_t *len, FILE *fp, std::string upTo);

      template <typename OUT>
      inline std::vector<OUT> parseTypeInLine(const std::string line, int nToRead)
      {
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        boost::char_separator<char> sep{" ,;\t\n"};
        std::vector<OUT> resultnumbers;
        tokenizer tok{line, sep};
        int iteration = 0;
        for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
        {
          // we are only interested in the first value.
          // also, we cannot cast the strings that follow
          resultnumbers.push_back(boost::lexical_cast<OUT>(*it));
          ++iteration;
          if (iteration == nToRead)
          {
            break;
          }
        }
        return resultnumbers;
      }

      //// data
      // nr of data points to read
      int nAtoms; // number of atoms
      int nBonds;
      int nAtomTypes;
      int nBondTypes;

      // box sizes
      double Lx;
      double Ly;
      double Lz;

      // actual dimensional values
      std::map<int, double> masses;
      std::vector<long int> atomIds;
      std::vector<int> moleculeIds;
      std::vector<int> atomTypes;
      std::vector<double> atomX;
      std::vector<double> atomY;
      std::vector<double> atomZ;
      std::vector<int> atomNx;
      std::vector<int> atomNy;
      std::vector<int> atomNz;
      std::vector<long int> bondIds;
      std::vector<int> bondTypes;
      std::vector<long int> bondFrom;
      std::vector<long int> bondTo;
    };

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
        std::string line(cline);
        // skip empty lines
        if (this->shortenLineToSkip(&line))
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
        std::string line(cline);

        // skip empty lines
        if (!this->shortenLineToSkip(&line))
        {
          break;
        }
      }

      // Then, read masses, up until the next section ("atoms")
      do
      {
        std::string line(cline);

        // break at empty lines
        if (this->shortenLineToSkip(&line))
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
        std::string line(cline);

        // skip until empty lines
        if (!this->shortenLineToSkip(&line))
        {
          break;
        }
      } while ((getline(&cline, len, fp)) != -1);
    }

    bool DataFileParser::shortenLineToSkip(char *line)
    {
      std::string tempString = std::string(line);
      return this->shortenLineToSkip(&tempString);
    }

    bool DataFileParser::shortenLineToSkip(std::string *line)
    {
      boost::trim_left(*line);
      // trim comments
      if (contains(line, "#"))
      {
        std::vector<std::string> split;
        boost::split(split, *line, boost::is_any_of("#"));
        line = &split[0];
      }
      return line->empty();
    }

    void DataFileParser::readNs(const std::string line)
    {
      if (contains(&line, "atoms"))
      {
        this->nAtoms = (this->parseTypeInLine<int>(line, 1))[0];
      }
      else if (contains(&line, "bonds"))
      {
        this->nBonds = (this->parseTypeInLine<int>(line, 1))[0];
      }
      else if (contains(&line, "atom types"))
      {
        this->nAtomTypes = (this->parseTypeInLine<int>(line, 1))[0];
      }
      else if (contains(&line, "bond types"))
      {
        this->nBondTypes = (this->parseTypeInLine<int>(line, 1))[0];
      }
      else if (contains(&line, "xlo xhi"))
      {
        std::vector<double> parsedL = this->parseTypeInLine<double>(line, 2);
        this->Lx = parsedL[1] - parsedL[0];
      }
      else if (contains(&line, "ylo yhi"))
      {
        std::vector<double> parsedL = this->parseTypeInLine<double>(line, 2);
        this->Ly = parsedL[1] - parsedL[0];
      }
      else if (contains(&line, "zlo zhi"))
      {
        std::vector<double> parsedL = this->parseTypeInLine<double>(line, 2);
        this->Lz = parsedL[1] - parsedL[0];
      }
    }

    void DataFileParser::readMass(const std::string line)
    {
      int iteration = 0;
      int key = 0;
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      boost::char_separator<char> sep{" ,;\t\n"};
      tokenizer tok{line, sep};
      for (tokenizer::iterator beg = tok.begin(); beg != tok.end(); ++beg)
      {
        if (iteration == 0)
        {
          key = boost::lexical_cast<int>(*beg);
        }
        else if (iteration == 1)
        {
          // we just override duplicate keys
          this->masses[key] = boost::lexical_cast<double>(*beg);
        }
        else
        {
          throw std::runtime_error("Too many fields tokenized when reading masses");
        }
        iteration++;
      }

      if (iteration != 2)
      {
        throw std::runtime_error("Incorrect nr of fields tokenized when reading masses");
      }
    }

    void DataFileParser::readAtom(std::string line)
    {
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      boost::char_separator<char> sep{" ,;\t\n"};
      tokenizer tok{line, sep};
      tokenizer::iterator it1, it2 = tok.begin();
      it1 = it2;
      this->atomIds.push_back(boost::lexical_cast<int>(*it2));
      std::advance(it2, 1);
      this->moleculeIds.push_back(boost::lexical_cast<int>(*it2));
      std::advance(it2, 1);
      this->atomTypes.push_back(boost::lexical_cast<int>(*it2));
      std::advance(it2, 1);
      this->atomX.push_back(boost::lexical_cast<double>(*it2));
      std::advance(it2, 1);
      this->atomY.push_back(boost::lexical_cast<double>(*it2));
      std::advance(it2, 1);
      this->atomZ.push_back(boost::lexical_cast<double>(*it2));
      std::advance(it2, 1);
      this->atomNx.push_back(boost::lexical_cast<int>(*it2));
      std::advance(it2, 1);
      this->atomNy.push_back(boost::lexical_cast<int>(*it2));
      std::advance(it2, 1);
      this->atomNz.push_back(boost::lexical_cast<int>(*it2));
    }

    void DataFileParser::readBond(std::string line)
    {
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      boost::char_separator<char> sep{" ,;\t\n"};
      tokenizer tok{line, sep};
      tokenizer::iterator it1, it2 = tok.begin();
      it1 = it2;
      this->bondIds.push_back(boost::lexical_cast<long int>(*it2));
      std::advance(it2, 1);
      this->bondTypes.push_back(boost::lexical_cast<int>(*it2));
      std::advance(it2, 1);
      this->bondFrom.push_back(boost::lexical_cast<long int>(*it2));
      std::advance(it2, 1);
      this->bondTo.push_back(boost::lexical_cast<long int>(*it2));
    }
  }
}

#endif
