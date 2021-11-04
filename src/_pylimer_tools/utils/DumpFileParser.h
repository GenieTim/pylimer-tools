#ifndef DUMP_FILE_PARSER_H
#define DUMP_FILE_PARSER_H

#include <string>
#include <vector>
#include <algorithm>
#include "StringUtil.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <map>
#include <cctype>
#include <cstring>
#include <filesystem>

namespace pylimer_tools
{
  namespace utils
  {

    class DumpFileParser
    {
    public:
      void read(const std::string filePath);
      template <typename OUT>
      std::vector<OUT> getValuesForAt(const int index, const std::string &headerKey, const std::string &column);
      template <typename OUT>
      std::vector<OUT> getValuesForAt(const int index, const std::string &headerKey, const int column);
      int getLength() { return this->data.size(); }
      bool hasKey(std::string headerKey)
      {
        if (this->getLength() > 0)
        {
          // assumption: each step has the same keys
          return this->data[0].contains(headerKey);
        }
        return false;
      }
      bool keyHasColumn(std::string headerKey, std::string column)
      {
        const auto colItIdx = std::find(this->headerColMap[headerKey].begin(), this->headerColMap[headerKey].end(), column);
        if (this->headerColMap[headerKey].end() == colItIdx)
        {
          return false;
        }
        return true;
      }
      bool keyHasDirectionalColumn(std::string headerKey, std::string dirPraefix, std::string dirSuffix)
      {
        return this->keyHasColumn(headerKey, dirPraefix + "x" + dirSuffix) && this->keyHasColumn(headerKey, dirPraefix + "y" + dirSuffix) && this->keyHasColumn(headerKey, dirPraefix + "z" + dirSuffix);
      }

    private:
      bool shortenLineToSkip(std::string *line);
      void skipEmptyLines(char *cline, size_t *len, FILE *fp);
      std::string cleanHeader(std::string header);

      template <typename IN>
      inline std::vector<IN> parseTypesInLine(std::string line)
      {
        std::vector<IN> resultnumbers;
        boost::tokenizer<> tok(line);
        std::transform(tok.begin(), tok.end(), std::back_inserter(resultnumbers),
                       &boost::lexical_cast<IN, std::string>);
        return resultnumbers;
      }

      //// data
      typedef std::map<std::string, std::string> data_item_t;
      std::vector<data_item_t> data;
      std::map<std::string, std::vector<std::string> > headerColMap;
    };

    void DumpFileParser::read(const std::string filePath)
    {
      if (!std::filesystem::exists(filePath))
      {
        throw std::invalid_argument("File to read (" + filePath + ") does not exist.");
      }
      char *cline = NULL;

      size_t len = 0;

      FILE *fp = fopen(filePath.c_str(), "r");
      if (fp == NULL)
        throw std::runtime_error("Failed to open data file to read.");

      // read everything until the first key
      while ((getline(&cline, &len, fp)) != -1)
      {
        std::string line(cline);
        // skip empty lines: break when not empty
        if (!this->shortenLineToSkip(&line))
        {
          break;
        }
      }

      // Assemble CSV data for all keys
      std::string newGroupKey(cline); // new group key: key for a new timestep (group)
      std::string currentKey = this->cleanHeader(std::string(cline));
      data_item_t dataItem;

      while ((getline(&cline, &len, fp)) != -1)
      {
        std::string line(cline);
        // skip empty lines
        if (this->shortenLineToSkip(&line))
        {
          continue;
        }
        // new header
        if (boost::algorithm::starts_with(line, "ITEM:"))
        {
          currentKey = this->cleanHeader(line);
        }
        else
        {
          if (!dataItem.contains(currentKey))
          {
            dataItem.insert_or_assign(currentKey, "");
          }
          dataItem.insert_or_assign(currentKey, dataItem.at(currentKey) + line);
        }

        if (line == newGroupKey)
        {
          // new timestep
          this->data.push_back(dataItem);
          dataItem = data_item_t();
        }
      }
      // last timestep
      this->data.push_back(dataItem);

      fclose(fp);
      if (cline)
      {
        free(cline);
      }
    }

    std::string DumpFileParser::cleanHeader(std::string headerToClean)
    {
      boost::algorithm::replace_first(headerToClean, "ITEM: ", "");
      boost::tokenizer<> tok(headerToClean);
      std::string newHeader = "";
      std::vector<std::string> columns;
      for (boost::tokenizer<>::iterator beg = tok.begin(); beg != tok.end(); ++beg)
      {
        if (isUpper(*beg))
        {
          newHeader.append(*beg);
          newHeader.append(" ");
        }
        else
        {
          columns.push_back(*beg);
        }
      }
      boost::algorithm::trim_right(newHeader);
      if (!this->headerColMap.contains(newHeader))
      {
        this->headerColMap.insert_or_assign(newHeader, columns);
      }

      return newHeader;
    }

    template <typename OUT>
    std::vector<OUT> DumpFileParser::getValuesForAt(const int index, const std::string &headerKey, const std::string &column)
    {
      // detect index of column
      const auto colItIdx = std::find(this->headerColMap[headerKey].begin(), this->headerColMap[headerKey].end(), column);
      if (this->headerColMap[headerKey].end() == colItIdx)
      {
        throw std::invalid_argument("Column '" + column + "' not found for header '" + headerKey + "'");
      }
      const int colIdx = colItIdx - this->headerColMap[headerKey].begin();
      return this->getValuesForAt<OUT>(index, headerKey, colIdx);
    };

    template <typename OUT>
    std::vector<OUT> DumpFileParser::getValuesForAt(const int index, const std::string &headerKey, const int colIdx)
    {
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      // this is clearly not the fastest way to access more than one column,
      // as the tokenizing happens again for each column.
      data_item_t dataItem = this->data[index];
      //
      std::string relevantData = dataItem[headerKey];
      std::istringstream f(relevantData);
      std::string line;
      std::vector<OUT> result;
      int lineNr = 0;

      boost::char_separator<char> sep{" ,;\t\n"};

      // all other lines
      while (std::getline(f, line))
      {
        tokenizer tok{line, sep};
        // boost::tokenizer<>::iterator it1, it2 = tok.begin();it1 = it2;
        // std::advance(it2, colIdx); cannot be used as boost gives segfault and I have not found a way to check bounds in another way
        int iteration = 0;
        bool found = false;
        for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
        {
          if (iteration == colIdx)
          {
            found = true;
            result.push_back(boost::lexical_cast<OUT>(*it));
            break;
          }
          iteration++;
        }
        if (!found)
        {
          throw std::runtime_error("Column index " + std::to_string(colIdx) + " out of range to get values for key " + headerKey + " in line " + std::to_string(lineNr) + ". Values are '" + line + "'.");
        }
        ++lineNr;
      }

      return result;
    };

    bool DumpFileParser::shortenLineToSkip(std::string *line)
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
  }
}

#endif
