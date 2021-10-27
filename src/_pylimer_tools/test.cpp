
#include <igraph/igraph.h>
#include <vector>
#include "utils/vector_utils.h"
#include <iostream>

int main(int argc, char **argv)
{
  std::cout << "Running" << std::endl;

  igraph_t graph;
  igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
  igraph_add_vertices(&graph, 10, 0);
  std::cout << "Added Vertices" << std::endl;
  igraph_vs_t allVertexIds;
  igraph_vs_all(&allVertexIds);
  igraph_vector_t degrees;
  igraph_vector_init(&degrees, 3);
  std::cout << "Checking degrees" << std::endl;
  igraph_degree(&graph, &degrees, allVertexIds, IGRAPH_ALL, false);
  std::cout << "Got degrees" << std::endl;
  igraph_vit_t vit;
  igraph_vit_create(&graph, allVertexIds, &vit);
  std::cout << "Created View" << std::endl;

  while (!IGRAPH_VIT_END(vit))
  {
    long int vertexId = (long int)IGRAPH_VIT_GET(vit);
    std::cout << "Vertex " << vertexId << " has degree " << igraph_vector_e(&degrees, vertexId) << std::endl;
    IGRAPH_VIT_NEXT(vit);
  }

  igraph_vit_destroy(&vit);
  std::cout << "Assembling vertices to remove" << std::endl;
  std::vector<long int> indices;
  //{ 0, 2, 3, 5 };
  indices.push_back(0);
  indices.push_back(2);
  indices.push_back(3);
  indices.push_back(5);
  igraph_vector_t verticesToRemove;
  igraph_vector_init(&verticesToRemove, 3);
  pylimer_tools::utils::StdVectorToIgraphVectorT(indices, &verticesToRemove);
  std::cout << "Starting to remove vertices" << std::endl;

  igraph_delete_vertices(&graph, igraph_vss_vector(&verticesToRemove));
  std::cout << "Removed vertices" << std::endl;

  igraph_vit_t vit3;
  igraph_vit_create(&graph, allVertexIds, &vit3);

  while (!IGRAPH_VIT_END(vit3))
  {
    long int vertexId = (long int)IGRAPH_VIT_GET(vit3);
    std::cout << "Vertex " << vertexId << std::endl;
    IGRAPH_VIT_NEXT(vit3);
  }

  std::cout << "Adding vertices again" << std::endl;

  igraph_vector_t degrees2;
  igraph_vector_init(&degrees2, 3);
  igraph_degree(&graph, &degrees2, allVertexIds, IGRAPH_ALL, false);
  igraph_add_vertices(&graph, 2, 0);
  igraph_vit_t vit2;
  igraph_vit_create(&graph, allVertexIds, &vit2);

  while (!IGRAPH_VIT_END(vit2))
  {
    long int vertexId = (long int)IGRAPH_VIT_GET(vit2);
    std::cout << "Vertex " << vertexId << " has degree " << (int)igraph_vector_e(&degrees2, vertexId) << std::endl;
    IGRAPH_VIT_NEXT(vit2);
  }
}
