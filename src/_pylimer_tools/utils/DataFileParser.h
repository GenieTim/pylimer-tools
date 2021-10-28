#ifndef DATA_FILE_PARSER_H
#define DATA_FILE_PARSER_H

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Universe.h"
#include "UniverseSequence.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <map>

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
      void readNs(const std::string line);
      void readMass(const std::string line);
      void readAtom(const char *line);
      void readBond(const char *line);
      void skipEmptyLines(char *cline, size_t *len, FILE *fp);
      void skipLinesTo(char *cline, size_t *len, FILE *fp, std::string upTo);

      template <typename IN>
      inline std::vector<IN> parseTypeInLine(std::string line)
      {
        std::vector<IN> resultnumbers;
        boost::tokenizer<> tok(line);
        std::transform(tok.begin(), tok.end(), std::back_inserter(resultnumbers),
                       &boost::lexical_cast<IN, std::string>);
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
      char *cline = NULL;
      char *eof;

      size_t len = 0;

      FILE *fp = fopen(filePath.c_str(), "r");
      if (fp == NULL)
        throw std::runtime_error("Failed to open data file to read.");

      // skip 1st line of file
      eof = fgets(cline, 256, fp);
      if (eof == nullptr)
      {
        throw std::runtime_error("Unexpected end of molecule file");
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
        //
        if (boost::algorithm::contains(line, "Masses"))
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

      // Then, read masses, up until the next section ("atoms")
      while ((getline(&cline, &len, fp)) != -1)
      {
        std::string line(cline);

        // skip empty lines
        if (this->shortenLineToSkip(&line))
        {
          continue;
        }
        // read masses until atoms
        if (boost::algorithm::contains(line, "Atoms"))
        {
          break;
        }
        // read the nr of data points to read afterwards
        this->readMass(line);
      }

      this->skipEmptyLines(cline, &len, fp);

      // Then, read atoms, up until the next section ("bonds")
      for (int i = 0; i < this->nAtoms; ++i)
      {
        if ((getline(&cline, &len, fp)) == -1)
        {
          throw std::runtime_error("Data file ended too early. Not enough atoms read.");
        }

        this->readAtom(cline);
      }

      // Then, read bonds
      this->skipLinesTo(cline, &len, fp, "Bonds");
      this->skipEmptyLines(cline, &len, fp);

      for (int i = 0; i < this->nBonds; i++)
      {
        if ((getline(&cline, &len, fp)) == -1)
        {
          throw std::runtime_error("Data file ended too early. Not enough atoms read.");
        }

        this->readBond(cline);
      }

      // we ignore angles etc. for now.

      fclose(fp);
      if (cline)
      {
        free(cline);
      }
    }

    void DataFileParser::skipLinesTo(char *cline, size_t *len, FILE *fp, std::string upTo)
    {
      while ((getline(&cline, len, fp)) != -1)
      {
        std::string line(cline);

        if (upTo.compare(line))
        {
          break;
        }
      }
    }

    void DataFileParser::skipEmptyLines(char *cline, size_t *len, FILE *fp)
    {
      while ((getline(&cline, len, fp)) != -1)
      {
        std::string line(cline);

        // skip empty lines
        if (!this->shortenLineToSkip(&line))
        {
          break;
        }
      }
    }

    bool DataFileParser::shortenLineToSkip(std::string *line)
    {
      boost::trim_left(line);
      // trim comments
      if (boost::algorithm::contains(line, "#"))
      {
        std::vector<std::string> split;
        boost::split(split, line, boost::is_any_of("#"));
        line = &split[0];
      }
      return line->compare("");
    }

    void DataFileParser::readNs(const std::string line)
    {
      if (boost::algorithm::contains(line, "atoms"))
      {
        this->nAtoms = (this->parseTypeInLine<int>(line))[0];
      }
      else if (boost::algorithm::contains(line, "bonds"))
      {
        this->nBonds = (this->parseTypeInLine<int>(line))[0];
      }
      else if (boost::algorithm::contains(line, "atom types"))
      {
        this->nAtomTypes = (this->parseTypeInLine<int>(line))[0];
      }
      else if (boost::algorithm::contains(line, "bond types"))
      {
        this->nBondTypes = (this->parseTypeInLine<int>(line))[0];
      }
      else if (boost::algorithm::contains(line, "xlo xhi"))
      {
        std::vector<double> parsedL = this->parseTypeInLine<double>(line);
        this->Lx = parsedL[1] - parsedL[0];
      }
      else if (boost::algorithm::contains(line, "ylo yhi"))
      {
        std::vector<double> parsedL = this->parseTypeInLine<double>(line);
        this->Ly = parsedL[1] - parsedL[0];
      }
      else if (boost::algorithm::contains(line, "zlo zhi"))
      {
        std::vector<double> parsedL = this->parseTypeInLine<double>(line);
        this->Lz = parsedL[1] - parsedL[0];
      }
    }

    void DataFileParser::readMass(const std::string line)
    {
      boost::tokenizer<> tok(line);
      int iteration = 0;
      int key = 0;
      for (boost::tokenizer<>::iterator beg = tok.begin(); beg != tok.end(); ++beg)
      {
        if (iteration == 0)
        {
          key = boost::lexical_cast<int>(beg);
        }
        else if (iteration == 1)
        {
          // we just override duplicate keys
          this->masses[key] = boost::lexical_cast<double>(beg);
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

    void DataFileParser::readAtom(const char *line)
    {
      boost::tokenizer<> tok(line);
      boost::tokenizer<>::iterator it1, it2 = tok.begin();
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
      this->atomNx.push_back(boost::lexical_cast<double>(*it2));
      std::advance(it2, 1);
      this->atomNy.push_back(boost::lexical_cast<double>(*it2));
      std::advance(it2, 1);
      this->atomNz.push_back(boost::lexical_cast<double>(*it2));
    }

    void DataFileParser::readBond(const char *line)
    {
      boost::tokenizer<> tok(line);
      boost::tokenizer<>::iterator it1, it2 = tok.begin();
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
