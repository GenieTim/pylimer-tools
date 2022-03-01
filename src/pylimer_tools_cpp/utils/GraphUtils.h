#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
extern "C" {
#include <igraph/igraph.h>
}
#include "VectorUtils.h"

namespace pylimer_tools {
namespace utils {

static std::vector<long int> getVerticesWithDegree(const igraph_t *graph,
                                                   std::vector<int> ofDegrees) {
  int graphSize = igraph_vcount(graph);
  igraph_vector_t degrees;
  if (igraph_vector_init(&degrees, graphSize)) {
    throw std::runtime_error("Failed to instantiate result vector.");
  }
  igraph_vs_t allVertexIds;
  igraph_vs_all(&allVertexIds);
  // complexity: O(|v|*d)
  if (igraph_degree(graph, &degrees, allVertexIds, IGRAPH_ALL, false)) {
    throw std::runtime_error("Failed to determine degree of vertices");
  }

  // NOTE: this is to omit the assumption, that the returned degree is
  // sequential for vertex 0, ..., |V|
  std::vector<long int> toSelect;
  igraph_vit_t vit;
  igraph_vit_create(graph, allVertexIds, &vit);
  while (!IGRAPH_VIT_END(vit)) {
    long int vertexId = static_cast<long int>(IGRAPH_VIT_GET(vit));
    int currentDegree = igraph_vector_e(&degrees, vertexId);
    for (int degree : ofDegrees) {
      if (currentDegree == degree) {
        toSelect.push_back(vertexId);
        break;
      }
    }
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vector_destroy(&degrees);
  igraph_vit_destroy(&vit);
  igraph_vs_destroy(&allVertexIds);

  return toSelect;
}

static std::vector<long int> getVerticesWithDegree(const igraph_t *graph,
                                                   int degree) {
  return getVerticesWithDegree(graph, std::vector<int>{degree});
}

static igraph_vs_t getVerticesWithDegreeSelector(const igraph_t *graph,
                                                 int degree) {
  // NOTE: this is to omit the assumption, that the returned degree is
  // sequential for vertex 0, ..., |V|
  std::vector<long int> toSelect = getVerticesWithDegree(graph, degree);

  igraph_vs_t result;
  igraph_vector_t toSelectVec;
  igraph_vector_init(&toSelectVec, toSelect.size());
  pylimer_tools::utils::StdVectorToIgraphVectorT(toSelect, &toSelectVec);
  igraph_vs_vector(&result, &toSelectVec);

  return result;
}

template <typename IN>
static bool graphHasVertexWithProperty(igraph_t *graph,
                                       std::string propertyName,
                                       IN propertyValue) {
  igraph_vector_t results;
  igraph_vector_init(&results, 1);
  if (igraph_cattribute_VANV(graph, propertyName.c_str(), igraph_vss_all(),
                             &results)) {
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
