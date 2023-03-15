#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#include <algorithm>
// #include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}
#include "./VectorUtils.h"

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

  class SimpleCycleFinder
  {
  public:
    SimpleCycleFinder(const igraph_t* ingraph, int minDegree = 2)
    {
      igraph_copy(&this->graph, ingraph);
      this->findStartingIndices(minDegree);
    };
    SimpleCycleFinder(const igraph_t* ingraph,
                      long int currentStartIndex,
                      long int currentNeighbourIndex,
                      std::set<long int> loopsFound,
                      std::vector<long int> startingIndices)
    {
      igraph_copy(&this->graph, ingraph);
      this->loopsFound = loopsFound;
      this->currentStartIndex = currentStartIndex;
      this->currentNeighbourIndex = currentNeighbourIndex;
      this->startingIndices = startingIndices;
    };

    // rule of three
    // 1. destructor (to destroy the graph)
    ~SimpleCycleFinder()
    {
      // in addition to basic fields being deleted, we need to clean up the
      // graph as is done in parent
      igraph_destroy(&this->graph);
    };
    // 2. copy constructor
    SimpleCycleFinder(const SimpleCycleFinder& other)
      : SimpleCycleFinder(&other.graph,
                          other.currentStartIndex,
                          other.currentNeighbourIndex,
                          other.loopsFound,
                          other.startingIndices){};
    // 3. copy assignment operator
    SimpleCycleFinder& operator=(SimpleCycleFinder other)
    {
      std::swap(this->graph, other.graph);
      std::swap(this->currentStartIndex, other.currentStartIndex);
      std::swap(this->currentNeighbourIndex, other.currentNeighbourIndex);
      std::swap(this->startingIndices, other.startingIndices);
      std::swap(this->loopsFound, other.loopsFound);
      return *this;
    };

    // actual computation stuff
    // find (next) simple cycles in the graph
    /**
     * @brief Find the next simple cycle in the graph
     *
     * Time Complexity: ca. O()
     *
     * @return std::vector<long int> the vertex indices of the found cycle.
     * Empty if no more cycle found.
     */
    std::vector<long int> findNext()
    {
      std::vector<long int> emptyResults;
      if (!this->hasNext()) {
        return emptyResults;
      }
      // the actual search.
      // naïve implementation: for a certain starting point, find the shortest
      // path to its neighbour, without actually using the direct connection.
      while (this->currentStartIndex < this->startingIndices.size()) {
        long int localStartingIndex =
          this->startingIndices[this->currentStartIndex];
        // find neighbourhood of the starting vertex
        igraph_vector_int_t neighbours;
        igraph_vector_int_init(&neighbours, 0);
        igraph_neighbors(
          &graph,
          &neighbours,
          localStartingIndex,
          IGRAPH_ALL); // TODO: do we want this unidirectional instead?
        // O(d)
        // loop through neighbours
        while (this->currentNeighbourIndex <
               igraph_vector_int_size(&neighbours)) {
          long int neighbourVertexId =
            igraph_vector_int_get(&neighbours, this->currentNeighbourIndex);
          this->currentNeighbourIndex += 1;
          // remove the edges to the neighbour
          std::vector<long int> edges =
            this->getEdgeIdsFromTo(localStartingIndex, neighbourVertexId);
          for (long int edgeId : edges) {
            igraph_delete_edges(&this->graph, igraph_ess_1(edgeId));
          }
          // find the loops
          igraph_vector_int_t verticesOfLoop;
          igraph_vector_int_init(&verticesOfLoop, 0);
          // LEMMA (TODO): it is sufficient to get only one shortest path,
          // all the others will be found from different starting points onwards
          igraph_get_shortest_path(&this->graph,
                                   &verticesOfLoop,
                                   NULL,
                                   localStartingIndex,
                                   neighbourVertexId,
                                   IGRAPH_ALL); // O(|V|+|E|)
          // re-add the edges
          for (long int edgeId : edges) {
            igraph_add_edge(&this->graph, currentStartIndex, neighbourVertexId);
          }
          // check if the loop found is one
          if (igraph_vector_int_size(&verticesOfLoop) > 0) {
            // assemble it and its hash
            long int loopHash = 0;
            std::vector<long int> results;
            results.reserve(igraph_vector_int_size(&verticesOfLoop));
            for (size_t i = 0; i < igraph_vector_int_size(&verticesOfLoop);
                 ++i) {
              long int loopVertexId = igraph_vector_int_get(&verticesOfLoop, i);
              loopHash =
                loopHash xor loopVertexId; // TODO: check whether xor is an
                                           // appropriate hash function
              results.push_back(loopVertexId);
            }
            // and return it if it is a new one
            if (!pylimer_tools::utils::map_has_key(this->loopsFound,loopHash)) {
              igraph_vector_int_destroy(&neighbours);
              igraph_vector_int_destroy(&verticesOfLoop);
              this->loopsFound.insert(loopHash);
              return results;
            }
          }
          igraph_vector_int_destroy(&verticesOfLoop);
        } // end loop through neighbours
        igraph_vector_int_destroy(&neighbours);
        this->currentStartIndex += 1;
        this->currentNeighbourIndex = 0;
      }
      return emptyResults;
    }

    /**
     * @brief Check whether more cycles could be found
     *
     * @return true if yes, more cycles might be found
     * @return false if no, certainly all cycles have been found
     */
    bool hasNext() const
    {
      return this->currentStartIndex < this->startingIndices.size();
    }

    /**
     * @brief Reset the internal iterator state.
     *
     */
    void reset()
    {
      this->currentStartIndex = 0;
      this->loopsFound.clear();
    }

    /**
     * @brief Find all simple cycles. 
     * NOTE: this might use more memory than you
     * might like.
     *
     * @return std::vector<std::vector<long int>>
     */
    std::vector<std::vector<long int>> findAllSimpleCycles()
    {
      std::vector<std::vector<long int>> results;
      results.reserve(
        this->startingIndices.size()); // an estimate: one cycle per junction
      while (this->hasNext()) {
        std::vector<long int> loop = this->findNext();
        if (!loop.empty()) {
          results.push_back(loop);
        }
      }
      results.shrink_to_fit();
      return results;
    }

  protected:
    /**
     * @brief Find the vertex indices of all junctions (more than {minDegree}
     * edges)
     *
     * @param minDegree The minimum degree that the starting vertices should
     * have.
     */
    void findStartingIndices(int minDegree = 2)
    {
      this->startingIndices.clear();
      this->startingIndices.reserve(
        igraph_vcount(&this->graph)); // this is the maximum. Maybe some other
                                      // number would be more efficient.
      igraph_vector_int_t degrees;
      igraph_vector_int_init(&degrees, igraph_vcount(&this->graph));
      igraph_degree(&this->graph,
                    &degrees,
                    igraph_vss_all(),
                    IGRAPH_ALL,
                    false); // TODO: want to count self-loops?
      for (size_t i = 0; i < igraph_vcount(&this->graph); ++i) {
        if (igraph_vector_int_get(&degrees, i) > minDegree) {
          this->startingIndices.push_back(i);
        }
      }
      igraph_vector_int_destroy(&degrees);
      // additionally, we need to respect e.g. free cycles
      if (minDegree >= 2) {
        // therfore, let's search for free clusters that do not fulfill the
        // requirement, and add random vertices from them to the starting
        // indices
        // to find the index of the current graph in the clusters, use the
        // attributes
        igraph_vector_t ids;
        igraph_vector_init(&ids, igraph_vcount(&this->graph));
        for (size_t i = 0; i < igraph_vcount(&this->graph); ++i) {
          igraph_vector_set(
            &ids, i, i + 1); // TODO: this assumes the vertices are 0...N
        }
        igraph_cattribute_VAN_setv(&this->graph, "id", &ids);
        // do the decomposition
        igraph_graph_list_t components;
        igraph_graph_list_init(&components, 0);
        if (igraph_decompose(&graph, &components, IGRAPH_WEAK, -1, 0)) {
          throw std::runtime_error("Failed to decompose graph.");
        }
        size_t NComponents = igraph_graph_list_size(&components);
        for (size_t i = 0; i < NComponents; ++i) {
          igraph_t* g = igraph_graph_list_get_ptr(&components, i);
          igraph_integer_t maxDegree;
          igraph_maxdegree(
            g, &maxDegree, igraph_vss_all(), IGRAPH_ALL, false); // loops again?
          if (maxDegree <= minDegree && maxDegree > 1) {
            this->startingIndices.push_back(
              igraph_cattribute_VAN(&this->graph, "id", 1));
          }
          igraph_destroy(g);
        }
        igraph_vector_destroy(&ids);
        igraph_graph_list_destroy(&components);
      }

      this->startingIndices.shrink_to_fit();
    }

    /**
     * @brief Get the edge ids of the edges between two vertices
     *
     * Main useage: check whether two vertices are connected twice
     *
     * @param vertexId1 id of the first vertex
     * @param vertexId2 id of the second vertex
     * @return std::vector<long int> the edge ids
     */
    std::vector<long int> getEdgeIdsFromTo(const long int vertexId1,
                                           const long int vertexId2) const
    {
      igraph_es_t edgeSelector;
      // igraph_es_pairs_small(
      //   &edgeSelector, IGRAPH_UNDIRECTED, vertexId1, vertexId2, -1);
      // // igraph_es_fromto(
      // //   &edgeSelector, igraph_vss_1(vertexId1), igraph_vss_1(vertexId2));
      // // igraph_es_pairs(
      // //   &edgeSelector, igraph_vss_1(vertexId1), igraph_vss_1(vertexId2));
      igraph_es_incident(&edgeSelector, vertexId1, IGRAPH_ALL);
      igraph_eit_t iterator;
      igraph_eit_create(&this->graph, edgeSelector, &iterator);
      std::vector<long int> results;
      results.reserve(IGRAPH_EIT_SIZE(iterator));
      while (!IGRAPH_EIT_END(iterator)) {
        long int edgeId = static_cast<long int>(IGRAPH_EIT_GET(iterator));

        igraph_integer_t vertex1OfEdge;
        igraph_integer_t vertex2OfEdge;
        igraph_edge(&this->graph, edgeId, &vertex1OfEdge, &vertex2OfEdge);

        if ((vertex1OfEdge == vertexId1 && vertex2OfEdge == vertexId2) ||
            (vertex1OfEdge == vertexId2 && vertex2OfEdge == vertexId1)) {
          results.push_back(edgeId);
        }
        IGRAPH_EIT_NEXT(iterator);
      }

      igraph_eit_destroy(&iterator);
      igraph_es_destroy(&edgeSelector);

      return results;
    }

    igraph_t graph;
    long int currentStartIndex = 0;
    long int currentNeighbourIndex = 0;
    std::set<long int> loopsFound;
    std::vector<long int> startingIndices;
  };
} // namespace utils
} // namespace pylimer_tools

#endif
