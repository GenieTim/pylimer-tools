#ifndef DUMP_FILE_PARSER_H
#define DUMP_FILE_PARSER_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include <any>
#include <fstream> // std::ifstream
#include <map>
#include <string>
#include <vector>

namespace pylimer_tools::utils {
// types
typedef std::map<std::string, std::vector<pylimer_tools::utils::CsvTokenizer>>
  data_item_t;

enum ReadableDumpFileSections : uint32_t
{
  TIMESTEP = (1 << 0),
  BOX = (1 << 1),
  ATOM = (1 << 2),
  EXTRA_ATOM = (1 << 3)
};

MAKE_FLAGS_ENUM(ReadableDumpFileSections, uint32_t);

struct ReadDumpFileSections
{
  std::vector<long int> timesteps;
  std::vector<pylimer_tools::entities::Box> boxes;
  std::vector<std::vector<pylimer_tools::entities::Atom>> atoms;
  std::vector<std::unordered_map<std::string, std::vector<double>>>
    extraAtomsData;
};

class DumpFileParser
{
public:
  DumpFileParser() {};
  DumpFileParser(const std::string filePath);

  // rule of three:
  // 1. destructor (to destroy the graph)
  ~DumpFileParser();
  // 2. copy constructor
  DumpFileParser(const DumpFileParser& src);
  // 3. copy assignment operator
  DumpFileParser& operator=(DumpFileParser src);

  void read();
  void finish();
  void readNGroups(const size_t start, const int N);
  void readGroupByIdx(const size_t i);
  void forgetAt(const size_t index);

  template<typename OUT>
  std::vector<OUT> getValuesForAt(const size_t index,
                                  const std::string headerKey,
                                  const std::string& column);
  template<typename OUT>
  std::vector<OUT> getValuesForAt(const size_t index,
                                  const std::string headerKey,
                                  const size_t column);
  // the next two methods are specializations for easier py binding
  std::vector<std::string> getStringValuesForAt(const size_t index,
                                                const std::string headerKey,
                                                const std::string column)
  {
    return this->getValuesForAt<std::string>(index, headerKey, column);
  };
  std::vector<double> getNumericValuesForAt(const size_t index,
                                            const std::string headerKey,
                                            const std::string column)
  {
    return this->getValuesForAt<double>(index, headerKey, column);
  };
  size_t getLength() const;
  bool hasKey(std::string headerKey) const;
  bool keyHasColumn(std::string headerKey, std::string column);
  bool keyHasDirectionalColumn(std::string headerKey,
                               std::string dirPraefix,
                               std::string dirSuffix);

  std::string getFile() const { return this->filePath; }

  std::vector<long int> readTimeSteps();
  std::vector<pylimer_tools::entities::Box> readBoxes();
  std::vector<std::vector<pylimer_tools::entities::Atom>> readAtoms();
  ReadDumpFileSections readDumpFileSections(
    ReadableDumpFileSections sectionsToRead);

private:
  std::string cleanHeader(std::string header);
  void rewind();
  void openFile();

  template<typename OUT>
  inline std::vector<OUT> parseTypesInLine(const std::string line)
  {
    std::vector<OUT> resultnumbers;
    pylimer_tools::utils::CsvTokenizer tokenizer(line);
    resultnumbers.reserve(tokenizer.getLength());
    for (size_t i = 0; i < tokenizer.getLength(); ++i) {
      resultnumbers.push_back(tokenizer.get<OUT>(i));
    }
    return resultnumbers;
  }

  //// data
  std::string filePath = "";
  std::string newGroupKey = "";
  std::string currentLine = "";
  std::ifstream file;
  size_t nrOfGroups;
  std::unordered_map<size_t, data_item_t> data;
  std::map<std::string, std::vector<std::string>> headerColMap;
  std::map<size_t, std::streampos> groupPosMap;
};

template<typename OUT>
std::vector<OUT>
DumpFileParser::getValuesForAt(const size_t index,
                               const std::string headerKey,
                               const std::string& column)
{
  // detect index of column
  size_t colIdx = 0;
  if (this->headerColMap.at(headerKey).size() > 1) {
    const auto colItIdx = std::find(this->headerColMap.at(headerKey).begin(),
                                    this->headerColMap.at(headerKey).end(),
                                    column);
    if (this->headerColMap.at(headerKey).end() == colItIdx) {
      throw std::invalid_argument("Column '" + column +
                                  "' not found for header '" + headerKey + "'");
    }
    colIdx = colItIdx - this->headerColMap.at(headerKey).begin();
  }
  return this->getValuesForAt<OUT>(index, headerKey, colIdx);
}

template<typename OUT>
std::vector<OUT>
DumpFileParser::getValuesForAt(const size_t index,
                               const std::string headerKey,
                               const size_t colIdx)
{
  if (!pylimer_tools::utils::map_has_key(this->data, index)) {
    // std::cout << "Could not find index " << index << "in data yet" <<
    // std::endl;
    this->readGroupByIdx(index);
  }
  // std::cout << "Requested values for index " << index << ", key " <<
  // headerKey << " and column " << colIdx << std::endl;

  data_item_t dataItem = this->data.at(index);
  //
  if (!pylimer_tools::utils::map_has_key(dataItem, headerKey)) {
    throw std::invalid_argument(headerKey + " is not a key in dataItem.");
  }
  std::vector<pylimer_tools::utils::CsvTokenizer> relevantData =
    dataItem.at(headerKey);
  std::vector<OUT> results;
  results.reserve(relevantData.size());

  for (const pylimer_tools::utils::CsvTokenizer& lineTok : relevantData) {
    results.push_back(lineTok.get<OUT>(colIdx));
  }

  return results;
}
}

#endif
