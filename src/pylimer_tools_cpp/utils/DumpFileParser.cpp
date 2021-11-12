#include <string>
#include <vector>
#include <algorithm>
#include "StringUtil.h"
#include "DumpFileParser.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"
#include <map>
#include <cctype>
#include <any>
#include <cstring>
#include <filesystem>
#include <fstream>      // std::ifstream

namespace pylimer_tools
{
  namespace utils
  {
    // types
    typedef std::map<std::string, std::vector<pylimer_tools::utils::CsvTokenizer>> data_item_t;

    void DumpFileParser::read(const std::string filePath)
    {
      if (!std::filesystem::exists(filePath))
      {
        throw std::invalid_argument("File to read (" + filePath + ") does not exist.");
      }

      std::string line;
      std::ifstream file;
      file.open(filePath);

      if (!file.is_open()) {
        throw std::invalid_argument("File to read file (" + filePath + "): failed to open.")
      }

      // read everything until the first key
      while (getline(file, line))
      {
        line = pylimer_tools::utils::trimLineOmitComment(line);
        // skip empty lines: break when not empty
        if (!line.empty())
        {
          break;
        }
      }

      // Assemble CSV data for all keys
      std::string newGroupKey(line); // new group key: key for a new timestep (group)
      std::string currentKey = this->cleanHeader(line);
      data_item_t dataItem;

      while ((getline(&cline, &len, fp)) != -1)
      {
        line = pylimer_tools::utils::trimLineOmitComment(line);
        // skip empty lines
        if (line.empty())
        {
          continue;
        }
        // new header
        if (pylimer_tools::utils::startsWith(line, "ITEM:"))
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

      file.close();
    }

    std::string DumpFileParser::cleanHeader(std::string headerToClean)
    {
      // "ITEM: ".size() = 6
      headerToClean.erase(0, 6);
      pylimer_tools::utils::CsvTokenizer tokenizer(headerToClean);

      std::string newHeader = "";
      std::vector<std::string> columns;
      for (size_t i = 0; i < tokenizer.getLength(); ++i)
      {
        std::string beg = tokenizer.get<std::string>(i);
        if (isUpper(beg))
        {
          newHeader.append(beg);
          newHeader.append(" ");
        }
        else
        {
          columns.push_back(beg);
        }
      }
      newHeader = pylimer_tools::utils::rtrim(newHeader);
      if (!this->headerColMap.contains(newHeader))
      {
        this->headerColMap.insert_or_assign(newHeader, columns);
      }

      return newHeader;
    }
  }
}
