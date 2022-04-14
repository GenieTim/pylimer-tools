#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}
#include "VectorUtils.h"

namespace pylimer_tools {
namespace utils {

  template<typename IN>
  static bool graphHasVertexWithProperty(igraph_t* graph,
                                         std::string propertyName,
                                         IN propertyValue)
  {
    igraph_vector_t results;
    igraph_vector_init(&results, 1);
    if (igraph_cattribute_VANV(
          graph, propertyName.c_str(), igraph_vss_all(), &results)) {
      throw std::runtime_error("Failed to query property " + propertyName);
    };
    std::vector<IN> resultsV;
    igraphVectorTToStdVector<IN>(&results, resultsV);
    for (IN result : resultsV) {
      if (result == propertyValue) {
        igraph_vector_destroy(&results);
        return true;
      }
    }
    igraph_vector_destroy(&results);
    return false;
  }
} // namespace utils
} // namespace pylimer_tools

#endif
