#ifndef AVE_FILE_READER
#define AVE_FILE_READER

#include <fstream>
#include <string>
#include <vector>

namespace pylimer_tools::utils {

class AveFileReader
{
public:
  explicit AveFileReader(const std::string& file)
    : filePath(file) {};

  const std::vector<std::string>& getColumnNames();
  int getNrOfRows();
  int getNrOfColumns();
  const std::vector<std::vector<double>>& getData();

  std::vector<double> autocorrelateColumn(int column,
                                          const std::vector<size_t>& dts);
  std::vector<double> autocorrelateColumnDifference(
    int column1,
    int column2,
    const std::vector<size_t>& dts);

  const std::string& getFilePath() const { return this->filePath; }

private:
  std::string filePath;
  // "cache"
  int numRows = -1;
  int numHeaderRows = -1;
  std::vector<std::string> columnNames;
  std::vector<std::vector<double>> data;
  bool dataIsRead = false;

  std::ifstream openFile() const;
  void readHeaderFromStream(std::ifstream& file,
                            std::string& firstDataLine,
                            bool& hasFirstDataLine);
  int getNrOfHeaderRows();
  void readData();
};

}

#endif //
