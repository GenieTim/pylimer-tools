#include "MEHPForceBalance2.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../utils/VectorUtils.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlopt.h>
#include <nlopt.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {

    //----------------------------------------------------------------
    // Structural Query
    //----------------------------------------------------------------
    /**
     * @brief Given a vertex id and a rail edge, returns the other two edges
     * that are not part of the rail
     */
    std::vector<size_t> MEHPForceBalance2::getOffRailConnectedEdgeIds(
      size_t vertexId,
      size_t railEdge)
    {
      INVALIDARG_EXP_IFN(igraph_cattribute_VAN(&this->graph,
                                               "type",
                                               vertexId) != this->splipLinkType,
                         "Can only search for rail around slip-links");

      // fetch the edges involved
      size_t otherRailEdge = this->getOtherRailEdgeId(vertexId, railEdge);
      igraph_es_t selector;
      igraph_es_incident(&selector, vertexId, IGRAPH_ALL);
      igraph_eit_t iterator;
      igraph_eit_create(&this->graph, selector, &iterator);
      std::vector<size_t> results;
      results.reserve(IGRAPH_EIT_SIZE(iterator) - 2);
      while (!IGRAPH_EIT_END(iterator)) {
        if ((IGRAPH_EIT_GET(iterator) != railEdge) &&
            (IGRAPH_EIT_GET(iterator) != otherRailEdge)) {
          results.push_back(IGRAPH_EIT_GET(iterator));
        }
        IGRAPH_EIT_NEXT(iterator);
      }

      igraph_es_destroy(&selector);
      igraph_eit_destroy(&iterator);
      assert(results.size() == 2 || results.size() == 0);

      return results;
    };

    /**
     * @brief Given a vertex and a connected edge, returns the edge in the
     * opposite direction
     *
     */
    size_t MEHPForceBalance2::getOtherRailEdgeId(size_t vertexId,
                                                 size_t railEdge)
    {
      INVALIDARG_EXP_IFN(igraph_cattribute_VAN(&this->graph,
                                               "type",
                                               vertexId) != this->splipLinkType,
                         "Can only search for rail around slip-links");
      // fetch the edges involved
      igraph_es_t selector;
      igraph_es_incident(&selector, vertexId, IGRAPH_ALL);
      igraph_eit_t iterator;
      igraph_eit_create(&this->graph, selector, &iterator);
      igraph_vector_int_t edgeIds;
      igraph_vector_int_init(&edgeIds, IGRAPH_EIT_SIZE(iterator));
      igraph_integer_t rail_idx = -1;
      size_t current_idx = 0;
      while (!IGRAPH_EIT_END(iterator)) {
        igraph_vector_int_push_back(&edgeIds, IGRAPH_EIT_GET(iterator));
        if (IGRAPH_EIT_GET(iterator) == railEdge) {
          rail_idx = current_idx;
        }
        IGRAPH_EIT_NEXT(iterator);
        current_idx += 1;
      }
      assert(rail_idx >= 0);

      igraph_vector_t parents;
      igraph_vector_init(&parents, igraph_vector_int_size(&edgeIds));
      igraph_cattribute_EANV(&this->graph, "parent_edge", selector, &parents);

      igraph_integer_t railParent = igraph_vector_get(&parents, rail_idx);

      // iff the link is a cross-link, or the
      // slip-link is not involved twice with the same strand (rail),
      // the following method is sufficient
      std::vector<size_t> results;
      for (size_t i = 0; i < igraph_vector_size(&parents); ++i) {
        if ((static_cast<igraph_integer_t>(igraph_vector_get(&parents, i)) ==
             railParent) &&
            (i != rail_idx)) {
          results.push_back(igraph_vector_int_get(&edgeIds, i));
        }
      }

      igraph_es_destroy(&selector);
      igraph_eit_destroy(&iterator);
      igraph_vector_destroy(&parents);

      assert(results.size() > 0);
      if (results.size() == 1) {
        return results[0];
      }

      assert((results.size() - 1) % 2 == 0);

      // otherwise, we have to go the complicated way:
      // figure out, where the rail is currently going through.
      igraph_vector_t parentEdges;
      igraph_vector_init(&parentEdges, igraph_ecount(&this->graph));
      igraph_cattribute_EANV(&this->graph,
                             "parent_edge",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &parentEdges);
      igraph_vector_int_t allEdgesOfParent;
      igraph_vector_int_init(&allEdgesOfParent, 0);
      igraph_vector_int_reserve(&allEdgesOfParent, 4);
      for (size_t i = 0; i < igraph_ecount(&this->graph); ++i) {
        if (static_cast<igraph_integer_t>(igraph_vector_get(&parentEdges, i)) ==
            railParent) {
          igraph_vector_int_push_back(&allEdgesOfParent, i);
          igraph_cattribute_EAN_set(&this->graph, "prev_edge_id", i, i);
          igraph_integer_t from, to;
          igraph_edge(&this->graph, i, &from, &to);
          igraph_cattribute_VAN_set(&this->graph, "prev_vertex_id", from, from);
          igraph_cattribute_VAN_set(&this->graph, "prev_vertex_id", to, to);
        }
      }

      igraph_t subgraph;
      igraph_empty(&subgraph, 0, IGRAPH_UNDIRECTED);
      igraph_subgraph_from_edges(
        &this->graph, &subgraph, igraph_ess_vector(&allEdgesOfParent), true);

      igraph_vector_int_t edgesOnPath;
      igraph_vector_int_init(&edgesOnPath, igraph_ecount(&subgraph));
      igraph_vector_int_t verticesOnPath;
      igraph_vector_int_init(&verticesOnPath, igraph_vcount(&subgraph));
      igraph_eulerian_path(&subgraph, &edgesOnPath, &verticesOnPath);
      assert(igraph_vector_int_size(&edgesOnPath) == igraph_ecount(&subgraph));

      size_t result = 0;
      bool foundResult = false;
      for (size_t i = 0; i < igraph_vector_int_size(&edgesOnPath); ++i) {
        if (igraph_cattribute_EAN(&subgraph, "prev_edge_id", i) == railEdge) {
          // yay, found the rail.
          // now: on either side of this is the
          // link we talk about. We are interested in the edge on the other
          // side!
          igraph_integer_t from, to;
          igraph_edge(&subgraph, i - 1, &from, &to);
          if ((igraph_cattribute_VAN(&subgraph, "prev_vertex_id", from) ==
               vertexId) ||
              (igraph_cattribute_VAN(&subgraph, "prev_vertex_id", to) ==
               vertexId)) {
            // i-1 is on-rail
            result = igraph_cattribute_EAN(&subgraph, "prev_edge_id", i - 1);
            foundResult = true;
          } else {
            igraph_edge(&subgraph, i + 1, &from, &to);
            assert((igraph_cattribute_VAN(&subgraph, "prev_vertex_id", from) ==
                    vertexId) ||
                   (igraph_cattribute_VAN(&subgraph, "prev_vertex_id", to) ==
                    vertexId));
            // i+1 is on-rail
            result = igraph_cattribute_EAN(&subgraph, "prev_edge_id", i + 1);
            foundResult = true;
          }
        }
      }

      igraph_vector_destroy(&parentEdges);
      igraph_vector_int_destroy(&allEdgesOfParent);
      igraph_vector_int_destroy(&edgesOnPath);
      igraph_vector_int_destroy(&verticesOnPath);

      assert(foundResult);
      return result;
    }

    //----------------------------------------------------------------
    // Structural Adjustments
    //----------------------------------------------------------------

    /**
     * @brief Remove a certain, 2-functional link from the structures
     *
     * @param net
     * @param displacements
     * @param linkIdx
     */
    void MEHPForceBalance2::remove2fLink(const size_t linkIdx)
    {
      if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
        this->updateGraph();
      }

      std::vector<size_t> neighbours = this->getNeighbourLinkIndices(linkIdx);
      if (neighbours.size() == 0) {
        igraph_delete_vertices(&this->graph, igraph_vss_1(linkIdx));
        this->net.isUpToDate = false;
        return;
      }
      RUNTIME_EXP_IFN(neighbours.size() == 2,
                      "Expect f = 2 to have 2 neighbours");

      igraph_add_edge(&this->graph, neighbours[0], neighbours[1]);
      // verify our assumption of the new edge id
      size_t newEdgeId = igraph_ecount(&this->graph) - 1;
      igraph_integer_t from, to;
      igraph_edge(&this->graph, newEdgeId, &from, &to);
      RUNTIME_EXP_IFN((from == neighbours[0] && to == neighbours[1]),
                      "Assumption made on edges is incorrect");

      // fetch the edges involved
      igraph_es_t selector;
      igraph_es_incident(&selector, linkIdx, IGRAPH_ALL);
      igraph_eit_t iterator;
      igraph_eit_create(&this->graph, selector, &iterator);
      igraph_vector_int_t edgeIds;
      igraph_vector_int_init(&edgeIds, 0);
      igraph_vector_int_reserve(&edgeIds, IGRAPH_EIT_SIZE(iterator));
      RUNTIME_EXP_IFN(IGRAPH_EIT_SIZE(iterator) == 2,
                      "Expected f = 2 to have 2 edges");
      while (!IGRAPH_EIT_END(iterator)) {
        igraph_vector_int_push_back(&edgeIds, IGRAPH_EIT_GET(iterator));
        IGRAPH_EIT_NEXT(iterator);
      }
      igraph_es_destroy(&selector);
      igraph_eit_destroy(&iterator);

      size_t removedEdge1 = std::min(igraph_vector_int_get(&edgeIds, 0),
                                     igraph_vector_int_get(&edgeIds, 1));
      size_t removedEdge2 = std::max(igraph_vector_int_get(&edgeIds, 0),
                                     igraph_vector_int_get(&edgeIds, 1));

      igraph_cattribute_EAN_set(
        &this->graph,
        "parent_edge",
        newEdgeId,
        igraph_cattribute_EAN(&this->graph, "parent_edge", removedEdge1));

      // merge edge properties
      if (igraph_cattribute_EAN(&this->graph, "parent_edge", removedEdge1) !=
          igraph_cattribute_EAN(&this->graph, "parent_edge", removedEdge2)) {
        // two different "parent" springs -> need to recalculate some stuff
        // probably only if link is cross-link
        this->combineParentSprings(
          igraph_cattribute_EAN(&this->graph, "parent_edge", removedEdge1),
          igraph_cattribute_EAN(&this->graph, "parent_edge", removedEdge2));
      } else {
        // the same "parent" spring
        igraph_cattribute_EAN_set(
          &this->graph,
          "partition_fraction",
          newEdgeId,
          igraph_cattribute_EAN(
            &this->graph, "partition_fraction", removedEdge1) +
            igraph_cattribute_EAN(
              &this->graph, "partition_fraction", removedEdge2));
      }

      // "new" spawning edge will be between the two neighbours
      igraph_integer_t from, to;
      igraph_edge(&this->graph, removedEdge1, &from, &to);
      double removedEdgeBoxPrefix = (from == linkIdx) ? 1. : -1.;
      igraph_edge(&this->graph, removedEdge2, &from, &to);
      double keptEdgeBondPrefix = (from == linkIdx) ? -1. : 1.;

      for (std::string dir : { "x", "y", "z" }) {
        igraph_cattribute_EAN_set(
          &this->graph,
          ("bond_box_" + dir).c_str(),
          newEdgeId,
          // TODO: rethink the prefixes
          igraph_cattribute_EAN(
            &this->graph, ("bond_box_" + dir).c_str(), removedEdge2) *
              keptEdgeBondPrefix +
            igraph_cattribute_EAN(
              &this->graph, ("bond_box_" + dir).c_str(), removedEdge1) *
              removedEdgeBoxPrefix);
      }

      igraph_delete_edges(&this->graph,
                          igraph_ess_1(std::max(removedEdge1, removedEdge2)));
      igraph_delete_edges(&this->graph,
                          igraph_ess_1(std::min(removedEdge1, removedEdge2)));

      this->net.isUpToDate = false;
    }

    /**
     * @brief marks a certain "parent" spring as non-existing
     */
    void MEHPForceBalance2::combineParentSprings(size_t springIdxBefore,
                                                 size_t springIdxNow)
    {
      INVALIDARG_EXP_IFN(springIdxBefore != springIdxNow,
                         "Will only replace higher with lower spring idx");
      if (springIdxBefore < springIdxNow) {
        std::swap(springIdxBefore, springIdxNow);
      }
      if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
        this->updateGraph();
      }

      igraph_vector_t parentEdges;
      igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
      igraph_cattribute_EANV(&this->graph,
                             "parent_edge",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &parentEdges);
      igraph_vector_t partitionFraction;
      igraph_vector_init(&partitionFraction, this->net.nrOfPartialSprings);
      igraph_cattribute_EANV(&this->graph,
                             "partition_fraction",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &partitionFraction);

      double scalingFactorRemoved =
        (this->net.springsContourLength[springIdxBefore]) /
        (this->net.springsContourLength[springIdxBefore] +
         this->net.springsContourLength[springIdxNow]);
      double scalingFactorNow =
        (this->net.springsContourLength[springIdxNow]) /
        (this->net.springsContourLength[springIdxBefore] +
         this->net.springsContourLength[springIdxNow]);

      for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        // update the new value of the parent edge
        if (static_cast<size_t>(igraph_vector_get(&parentEdges, i)) ==
            springIdxBefore) {
          // ...if the edge was merged,
          igraph_vector_set(&parentEdges, i, springIdxNow);
          igraph_vector_set(&partitionFraction,
                            i,
                            igraph_vector_get(&partitionFraction, i) *
                              scalingFactorRemoved);
        } else if (igraph_vector_get(&parentEdges, i) > springIdxBefore) {
          igraph_vector_set(
            &parentEdges, i, igraph_vector_get(&parentEdges, i) - 1);
        } else if (igraph_vector_get(&parentEdges, i) == springIdxNow) {
          // ...the edge was merged into
          // update partition fraction
          igraph_vector_set(&partitionFraction,
                            i,
                            igraph_vector_get(&partitionFraction, i) *
                              scalingFactorNow);
        }
      }

      igraph_cattribute_EAN_setv(&this->graph, "parent_edge", &parentEdges);
      igraph_cattribute_EAN_setv(
        &this->graph, "partition_fraction", &partitionFraction);

      igraph_vector_destroy(&parentEdges);
      igraph_vector_destroy(&partitionFraction);

      this->net.springsContourLength(springIdxNow) +=
        this->net.springsContourLength(springIdxBefore);
      pylimer_tools::utils::removeRow(this->net.springsContourLength,
                                      springIdxBefore);

      this->net.isUpToDate = false;
    }
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
