#pragma once

// #include <iostream>
#include <string>
#include <vector>

extern "C"
{
#include <igraph/igraph.h>
}

#include "./VectorUtils.h"

static inline igraph_integer_t
castToIgraphInt(const igraph_real_t c)
{
  return std::lround(c);
}

template<typename INTT>
static inline INTT
igraphRealToInt(const igraph_real_t c)
{
  return static_cast<INTT>(std::lround(c));
}

namespace pylimer_tools::utils {
template<typename IN>
static bool
graphHasVertexWithProperty(igraph_t* graph,
                           const std::string propertyName,
                           IN propertyValue)
{
  INVALIDARG_EXP_IFN(propertyValue != 0,
                     "Check for default property value not supported.");
  igraph_vector_t results;
  igraph_vector_init(&results, 1);
  if (!igraph_cattribute_has_attr(
        graph, IGRAPH_ATTRIBUTE_VERTEX, propertyName.c_str())) {
    return false;
  }
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

[[maybe_unused]] static void
copyVertexProperties(const igraph_t* sourceGraph,
                     const igraph_integer_t sourceVertex,
                     igraph_t* targetGraph,
                     const igraph_integer_t targetVertex,
                     const std::vector<std::string>& propertyNames)
{
  for (const std::string& propertyName : propertyNames) {
    igraph_cattribute_VAN_set(
      targetGraph,
      propertyName.c_str(),
      targetVertex,
      igraph_cattribute_VAN(sourceGraph, propertyName.c_str(), sourceVertex));
  }
}

[[maybe_unused]] static void
copyEdgeProperties(const igraph_t* sourceGraph,
                   const igraph_integer_t sourceEdge,
                   igraph_t* targetGraph,
                   const igraph_integer_t targetEdge,
                   const std::vector<std::string>& propertyNames)
{
  for (const std::string& propertyName : propertyNames) {
    igraph_cattribute_EAN_set(
      targetGraph,
      propertyName.c_str(),
      targetEdge,
      igraph_cattribute_EAN(sourceGraph, propertyName.c_str(), sourceEdge));
  }
}
}