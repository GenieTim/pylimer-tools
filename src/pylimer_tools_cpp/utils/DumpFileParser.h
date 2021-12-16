#ifndef DUMP_FILE_PARSER_H
#define DUMP_FILE_PARSER_H

#include <string>
#include <vector>
#include <algorithm>
#include "StringUtil.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"
#include <map>
#include <cctype>
#include <any>
#include <cstring>
#include <filesystem>
#include <fstream> // std::ifstream

namespace pylimer_tools
{
  namespace utils
  {
    // types
    typedef std::map<std::string, std::vector<pylimer_tools::utils::CsvTokenizer>> data_item_t;

    class DumpFileParser
    {
    public:
      void read(const std::string filePath);
      void startReading(const std::string filePath);
      void finish();
      void readNGroups(const int N);
      void forgetAt(const int index);

      bool isFinishedReading() const { return this->finishedReading; }
      template <typename OUT>
      std::vector<OUT> getValuesForAt(const int index, const std::string headerKey, const std::string &column);
      template <typename OUT>
      std::vector<OUT> getValuesForAt(const int index, const std::string headerKey, const int column);
      // the next two methods are specializations for easier py binding
      std::vector<std::string> getStringValuesForAt(const int index, const std::string headerKey, const std::string column)
      {
        return this->getValuesForAt<std::string>(index, headerKey, column);
      };
      std::vector<double> getNumericValuesForAt(const int index, const std::string headerKey, const std::string column)
      {
        return this->getValuesForAt<double>(index, headerKey, column);
      };
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
        // std::cout << "Searching for " << headerKey << " " << dirPraefix << dirSuffix << " in " << pylimer_tools::utils::join(this->headerColMap[headerKey].begin(), this->headerColMap[headerKey].end(), std::string(" ")) << std::endl;
        return this->keyHasColumn(headerKey, dirPraefix + "x" + dirSuffix) && this->keyHasColumn(headerKey, dirPraefix + "y" + dirSuffix) && this->keyHasColumn(headerKey, dirPraefix + "z" + dirSuffix);
      }

    private:
      std::string cleanHeader(std::string header);

      template <typename OUT>
      inline std::vector<OUT> parseTypesInLine(std::string line)
      {
        std::vector<OUT> resultnumbers;
        pylimer_tools::utils::CsvTokenizer tokenizer(line);
        resultnumbers.reserve(tokenizer.getLength());
        for (size_t i = 0; i < tokenizer.getLength(); ++i)
        {
          resultnumbers.push_back(tokenizer.get<OUT>(i));
        }
        return resultnumbers;
      }

      //// data
      std::string newGroupKey;
      std::string currentLine;
      std::ifstream file;
      std::vector<data_item_t> data;
      std::map<std::string, std::vector<std::string>> headerColMap;
      bool finishedReading = false;
    };

    template <typename OUT>
    std::vector<OUT> DumpFileParser::getValuesForAt(const int index, const std::string headerKey, const std::string &column)
    {
      // detect index of column
      int colIdx = 0;
      if (this->headerColMap[headerKey].size() > 1)
      {
        const auto colItIdx = std::find(this->headerColMap[headerKey].begin(), this->headerColMap[headerKey].end(), column);
        if (this->headerColMap[headerKey].end() == colItIdx)
        {
          throw std::invalid_argument("Column '" + column + "' not found for header '" + headerKey + "'");
        }
        colIdx = colItIdx - this->headerColMap[headerKey].begin();
      }
      return this->getValuesForAt<OUT>(index, headerKey, colIdx);
    }

    template <typename OUT>
    std::vector<OUT> DumpFileParser::getValuesForAt(const int index, const std::string headerKey, const int colIdx)
    {
      data_item_t dataItem = this->data[index];
      //
      std::vector<pylimer_tools::utils::CsvTokenizer> relevantData = dataItem[headerKey];
      std::vector<OUT> results;

      for (pylimer_tools::utils::CsvTokenizer lineTok : relevantData)
      {
        results.push_back(lineTok.get<OUT>(colIdx));
      }

      return results;
    }
  }
}

#endif
