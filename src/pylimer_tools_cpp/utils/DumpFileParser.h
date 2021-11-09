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
#include <any>
#include <cstring>
#include <filesystem>

namespace pylimer_tools
{
  namespace utils
  {
    // types
    typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    typedef std::map<std::string, std::vector<pylimer_tools::utils::CsvTokenizer>> data_item_t;

    class DumpFileParser
    {
    public:
      void read(const std::string filePath);
      template <typename OUT>
      std::vector<OUT> getValuesForAt(const int index, const std::string &headerKey, const std::string &column);
      template <typename OUT>
      std::vector<OUT> getValuesForAt(const int index, const std::string &headerKey, const int column);
      template <typename... Ts>
      std::vector<std::variant<Ts...>> parseRow(std::string row);
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
      template <typename T>
      void pushBackParsedValue(tokenizer::iterator &it, std::vector<std::any> &target);
      template <typename T, typename... restTs>
      void pushBackParsedValues(tokenizer::iterator &it, std::vector<std::any> &target);
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
      std::vector<data_item_t> data;
      std::map<std::string, std::vector<std::string>> headerColMap;
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
      {
        throw std::runtime_error("Failed to open data file to read.");
      }

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
          dataItem[currentKey].push_back(pylimer_tools::utils::CsvTokenizer(line));
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
      data_item_t dataItem = this->data[index];
      //
      std::vector<pylimer_tools::utils::CsvTokenizer> relevantData = dataItem[headerKey];
      std::vector<OUT> results;

      for (pylimer_tools::utils::CsvTokenizer lineTok : relevantData)
      {
        int iteration = 0;
        results.push_back(lineTok.get<OUT>(colIdx));
      }

      return results;
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
