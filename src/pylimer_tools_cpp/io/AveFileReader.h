#ifndef AVE_FILE_READER
#define AVE_FILE_READER

#include <string>
#include <vector>

namespace pylimer_tools::utils {

class AveFileReader
{
public:
  AveFileReader(const std::string& file)
    : filePath(file) {};

  std::vector<std::string> getColumnNames();
  int getNrOfRows();
  int getNrOfColumns();
  std::vector<std::vector<double>> getData();

  std::vector<double> autocorrelateColumn(int column,
                                          const std::vector<size_t>& dts);
  std::vector<double> autocorrelateColumnDifference(
    int column1,
    int column2,
    const std::vector<size_t>& dts);

private:
  std::string filePath;
  // "cache"
  int numRows = -1;
  int numHeaderRows = -1;
  std::vector<std::string> columnNames;
  std::vector<std::vector<double>> data;

  int getNrOfHeaderRows();
};

}

#endif //
