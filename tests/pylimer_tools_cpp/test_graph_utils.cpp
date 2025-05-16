#include "../../src/pylimer_tools_cpp/utils/GraphUtils.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

extern "C"
{
#include <igraph/igraph.h>
}
#include <iostream>

namespace pu = pylimer_tools::utils;

TEST_CASE("igraph real casting works", "[GraphUtils][utils][header_tests]")
{
  SECTION("castToIgraphInt")
  {
    CHECK(castToIgraphInt(0.9) == 1);
    CHECK(castToIgraphInt(1.1) == 1);
    CHECK(castToIgraphInt(5.0) == 5);
  }

  SECTION("igraphRealToInt")
  {
    CHECK(igraphRealToInt<long int>(0.9) == 1);
    CHECK(igraphRealToInt<size_t>(1.1) == 1);
    CHECK(igraphRealToInt<int>(5.0) == 5);
  }
}

TEST_CASE("Graph Vertex property is detected",
          "[GraphUtils][utils][header_tests]")
{
  if (!igraph_has_attribute_table()) {
    igraph_set_attribute_table(&igraph_cattribute_table);
  }
  igraph_t graph;
  igraph_empty(&graph, 5, IGRAPH_UNDIRECTED);
  igraph_add_edge(&graph, 0, 1);
  igraph_cattribute_VAN_set(&graph, "test_attr", 0, 1.0);

  CHECK_FALSE(
    pu::graphHasVertexWithProperty<igraph_real_t>(&graph, "test_attr", 2.0));
  CHECK(
    pu::graphHasVertexWithProperty<igraph_real_t>(&graph, "test_attr", 1.0));
  CHECK_FALSE(pu::graphHasVertexWithProperty<size_t>(&graph, "test_attr", 2));
  CHECK(pu::graphHasVertexWithProperty<size_t>(&graph, "test_attr", 1));

  igraph_destroy(&graph);
}

TEST_CASE("Graph properties are copied correctly",
          "[GraphUtils][utils][header_tests]")
{
  if (!igraph_has_attribute_table()) {
    igraph_set_attribute_table(&igraph_cattribute_table);
  }
  igraph_t graph;
  igraph_empty(&graph, 5, IGRAPH_UNDIRECTED);
  igraph_add_edge(&graph, 0, 1);
  igraph_add_edge(&graph, 0, 2);
  igraph_t copyGraph;
  igraph_copy(&copyGraph, &graph);

  // add the attributes after the copy, so we can copy them manually
  igraph_cattribute_VAN_set(&graph, "test_attr", 0, 1.0);
  igraph_cattribute_VAN_set(&graph, "test_attr2", 0, 2.0);
  igraph_cattribute_EAN_set(&graph, "test_edge_attr", 0, 3.0);

  REQUIRE_FALSE(pu::graphHasVertexWithProperty(&copyGraph, "test_attr", 1.0));

  pylimer_tools::utils::copyVertexProperties(
    &graph, 0, &copyGraph, 0, { "test_attr", "test_attr2" });
  pu::copyEdgeProperties(&graph, 0, &copyGraph, 0, { "test_edge_attr" });

  CHECK(pu::graphHasVertexWithProperty(&copyGraph, "test_attr", 1.0));
}