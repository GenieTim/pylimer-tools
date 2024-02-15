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
    // Simulation / Optimization Procedures
    //----------------------------------------------------------------

    /**
     * FORCE RELAXATION
     */
    void MEHPForceBalance2::runForceRelaxation(
      BalanceRunMode mode,
      double damping,
      long int maxNrOfSteps, // default: 10000
      double xtol,
      const double initialResidualToUse,
      const StructureSimplificationMode simplificationMode,
      const double inactiveRemovalCutoff,
      const int outputFrequency,
      bool doInnerIterations,
      LinkSwappingMode allowSlipLinksToPassEachOther,
      const int swappingFrequency,
      const double oneOverSpringPartitionUpperLimit,
      const int nrOfCrosslinkSwapsAllowedPerSliplink)
    {
      // INVALIDARG_EXP_IFN(
      //   shouldRemoveInactiveCrosslinks == false &&
      //     remove2functionalCrosslinkers == true,
      //   "Removing 2-functional cross-links only makes sense when inactive "
      //   "cross-links may be removed too, during the procedure.");
      this->simulationHasRun = true;

      double removalTolerance =
        (inactiveRemovalCutoff > 0.0)
          ? inactiveRemovalCutoff
          : (0.25 *
             std::pow(this->net.vol / this->universe.getNrOfAtoms(), 1. / 3.));

      /* array allocation */
      std::vector<Eigen::ArrayXi> independentVertexSets;
      double maxDistanceMoved = 0.0;
      size_t indexOfMaxDistanceMoved = 0;
      // default = all
      std::vector<Eigen::ArrayXi> independentVertexsSpringSets;
      if (mode == BalanceRunMode::EIGEN_HEURISTIC) {
        std::tie(independentVertexSets, independentVertexsSpringSets) =
          getHeuristicallyIndependentCoordinateSets(this->net);
      } else if (mode == BalanceRunMode::EIGEN_RANDOM) {
        independentVertexSets = this->getRandomCoordinateSets(this->net);
      } else if (mode == BalanceRunMode::EIGEN_STRANDS) {
        independentVertexSets = { this->net.springPartCoordinateIndexA,
                                  this->net.springPartCoordinateIndexB };
      } else if (mode == BalanceRunMode::EIGEN_ALL) {
        independentVertexSets = { Eigen::ArrayXi::LinSpaced(
          3 * this->net.nrOfLinks, 0, 3 * this->net.nrOfLinks - 1) };
      }
      // this->getIndependentCoordinateSets(this->net);
      // { this->net.springPartCoordinateIndexA,
      // this->net.springPartCoordinateIndexB };
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(this->net,
                                             this->currentSpringPartitionsVec,
                                             oneOverSpringPartitionUpperLimit);
      const double initialResidual = (initialResidualToUse > 0.)
                                       ? initialResidualToUse
                                       : this->getDisplacementResidualNormFor(
                                           this->net, oneOverSpringPartitions);
      const double minN = this->net.springsContourLength.minCoeff();
      std::cout << "Starting force balance procedure "
                << "with " << initialResidual
                << " as initial residual, got requested "
                << initialResidualToUse
                // "with " << independentVertexSets.size() << "vertex sets."
                << std::endl;
      std::cout << "Swapping mode is " << allowSlipLinksToPassEachOther
                << " while simplification mode is " << simplificationMode
                << std::endl;
      std::cout << "Using oneOverSpringPartitionUpperLimit = "
                << oneOverSpringPartitionUpperLimit;
      double currentResidual = 0.0;
      double intermediateResidual = 0.0;
      size_t iterationsDone = 0;

      this->prepareAllOutputs();

      // actual loop
      do {
        if (allowSlipLinksToPassEachOther != LinkSwappingMode::NO_SWAPPING) {
          if (swappingFrequency > 0 &&
              (iterationsDone % swappingFrequency) == 0) {
            if (allowSlipLinksToPassEachOther ==
                LinkSwappingMode::SLIPLINKS_ONLY) {
              this->swapSlipLinks(this->net,
                                  this->currentSpringPartitionsVec,
                                  oneOverSpringPartitionUpperLimit);
            } else if (allowSlipLinksToPassEachOther == LinkSwappingMode::ALL) {
              this->swapSlipLinksInclXlinks(this->net,
                                            this->currentSpringPartitionsVec,
                                            oneOverSpringPartitionUpperLimit);
            } else if (allowSlipLinksToPassEachOther ==
                       LinkSwappingMode::ALL_MC) {
              this->moveSlipLinksToTheirBestBranch(
                this->net,
                this->currentSpringPartitionsVec,
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                false);
            } else if (allowSlipLinksToPassEachOther ==
                       LinkSwappingMode::ALL_MC_CYCLE) {
              this->moveSlipLinksToTheirBestBranch(
                this->net,
                this->currentSpringPartitionsVec,
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                true);
            } else {
              throw std::invalid_argument(
                "This swapping mode is currently not supported.");
            }
            oneOverSpringPartitions = this->assembleOneOverSpringPartition(
              this->net,
              this->currentSpringPartitionsVec,
              oneOverSpringPartitionUpperLimit);
          }

          if (!this->net.isUpToDate) {
            this->convertFromGraph();
          }
        }
        maxDistanceMoved = 0.0;
        currentResidual = 0.0;

        // place slip-link
        for (size_t link_idx = this->net.nrOfNodes;
             link_idx < this->net.nrOfLinks;
             ++link_idx) {
          // std::cout << "Handling " << link_idx << " of " << net.nrOfNodes
          //           << " / " << net.nrOfLinks << std::endl;

          // std::cout << "Still handling " << link_idx << " of " <<
          // net.nrOfNodes
          //           << " / " << net.nrOfLinks << std::endl;
          int innerIterationsDone = 0;
          do {
            double r2 =
              this->updateSpringPartition(this->net,
                                          this->currentSpringPartitionsVec,
                                          oneOverSpringPartitions,
                                          link_idx,
                                          oneOverSpringPartitionUpperLimit,
                                          allowSlipLinksToPassEachOther);
            double displacementDone =
              this->displaceToMeanPosition(this->net,
                                           this->currentSpringPartitionsVec,
                                           link_idx,
                                           oneOverSpringPartitionUpperLimit);
            innerIterationsDone += 1;
          } while (doInnerIterations && innerIterationsDone < 50);
        }

        intermediateResidual = this->getDisplacementResidualNormFor(
          this->net, oneOverSpringPartitions);

        if (mode == BalanceRunMode::EIGEN_RANDOM) {
          independentVertexSets = getRandomCoordinateSets(this->net);
        }
        // place cross-links
        if (mode == BalanceRunMode::ITERATIVE) {
          for (size_t link_idx = 0; link_idx < this->net.nrOfNodes;
               ++link_idx) {
            double distanceMoved =
              this->displaceToMeanPosition(this->net,
                                           this->currentSpringPartitionsVec,
                                           link_idx,
                                           oneOverSpringPartitionUpperLimit);
            if (distanceMoved > maxDistanceMoved) {
              maxDistanceMoved = distanceMoved;
              indexOfMaxDistanceMoved = link_idx;
            }
            // maxDistanceMoved = std::max(maxDistanceMoved, distanceMoved);
          }
        } else {
          for (size_t i = 0; i < independentVertexSets.size(); ++i) {
            Eigen::ArrayXi vertexSet = independentVertexSets[i];
            if (independentVertexsSpringSets.size() ==
                independentVertexSets.size()) {
              Eigen::ArrayXi independentVertexsSpringSet =
                independentVertexsSpringSets[i];
              maxDistanceMoved = std::max(
                maxDistanceMoved,
                this->displaceLinksToMeanPosition(this->net,
                                                  oneOverSpringPartitions,
                                                  independentVertexsSpringSet,
                                                  vertexSet,
                                                  damping));
            } else {
              maxDistanceMoved = std::max(
                maxDistanceMoved,
                this->displaceLinksToMeanPosition(
                  this->net, oneOverSpringPartitions, vertexSet, damping));
            }
          }
        }

        currentResidual = this->getDisplacementResidualNormFor(
          this->net, oneOverSpringPartitions);
        iterationsDone += 1;
        if (iterationsDone % 10 == 0) {
          if (simplificationMode ==
                StructureSimplificationMode::INACTIVE_ONLY ||
              simplificationMode == StructureSimplificationMode::ALL_TIM) {
            // std::cout << "Removing inactive cross-links" << std::endl;
            // default tolerance: 0.25*atom's cube length
            size_t nRemoved = this->removeInactiveCrosslinks(
              this->currentSpringPartitionsVec, removalTolerance);
            this->net.meanSpringContourLength =
              this->net.springsContourLength.size() > 0
                ? this->net.springsContourLength.mean()
                : 0.;
            if (nRemoved > 0) {
              std::cout << "Removed " << nRemoved << " inactive springs. "
                        << std::endl;
            }
            // this->validateNetwork(this->net,
            // this->currentDisplacements, this->currentSpringPartitionsVec);
          }
          if (simplificationMode == StructureSimplificationMode::X2F_ONLY ||
              simplificationMode == StructureSimplificationMode::ALL_TIM) {
            // std::cout << "Removing 2-f cross-links" << std::endl;
            size_t nRemoved = this->removeTwofunctionalCrosslinks();
            this->net.meanSpringContourLength =
              this->net.springsContourLength.size() > 0
                ? this->net.springsContourLength.mean()
                : 0.;
            if (nRemoved > 0) {
              std::cout << "Removed " << nRemoved
                        << " cross-linkers with f = 2. " << std::endl;
            }
            // this->validateNetwork(this->net,
            // this->currentDisplacements, this->currentSpringPartitionsVec);
          }
          if (simplificationMode == StructureSimplificationMode::ALL_ANDREI) {
            std::cout << "Removing cross-links and springs, Andrei's way"
                      << std::endl;
            this->doRemovalAndreisWay(this->currentSpringPartitionsVec,
                                      removalTolerance);
          }
          if (simplificationMode !=
              StructureSimplificationMode::NO_SIMPLIFICATION) {
            this->cleanupPrimaryLoopsInStructure();
            if (!this->net.isUpToDate) {
              this->convertFromGraph();
            }
            this->validateNetwork(this->net, this->currentSpringPartitionsVec);
          }
        }
        if (outputFrequency > 0 && iterationsDone % outputFrequency == 0) {
          this->handleOutput(iterationsDone);
        }
      } while (currentResidual / initialResidual > xtol &&
               iterationsDone < maxNrOfSteps && this->net.nrOfSprings > 0);

      // query solution & exit reason
      this->exitReason = (iterationsDone == maxNrOfSteps)
                           ? ExitReason::MAX_STEPS
                           : ExitReason::X_TOLERANCE;
      this->nrOfStepsDone += iterationsDone;
      std::cout << iterationsDone << " steps done. "
                << "Last max distance moved: " << maxDistanceMoved << std::endl;

      this->validateNetwork();
      this->currentSpringDistances = this->evaluateSpringDistances();
      this->currentPartialSpringDistances =
        this->evaluatePartialSpringDistances();
    }

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
        if ((castToIgraphInt(igraph_vector_get(&parents, i)) == railParent) &&
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
        if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) == railParent) {
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
     * @brief Remove a certain link from the structures, removing all
     * connections
     */
    void MEHPForceBalance2::removeLink(const size_t linkIdx)
    {
      if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
        this->updateGraph();
      }

      igraph_delete_vertices(&this->graph, igraph_vss_1(linkIdx));
      this->net.isUpToDate = false;
    };

    /**
     * @brief Remove a certain, 2-functional link from the structures, combining
     * the two strands
     *
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

      igraph_delete_vertices(&this->graph, igraph_vss_1(linkIdx));
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

    //----------------------------------------------------------------
    // Concrete structure modification
    //----------------------------------------------------------------

    /**
     * @brief Remove cross-links which do not have any springs with a certain
     * minimum length
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param tolerance
     */
    size_t MEHPForceBalance2::removeInactiveCrosslinks(
      Eigen::VectorXd& springPartitions,
      double tolerance)
    {
      size_t numRemoved = 0;
      size_t numRemovedInIteration = 0;
      do {
        this->removeInactiveParentEdges(tolerance);
        numRemovedInIteration = this->removeSubfunctionalVertices();
        numRemoved += numRemovedInIteration;
      } while (numRemovedInIteration > 0);

      return numRemoved;
    };

    /**
     * @brief Remove cross-linkers, springs and associated slip-links with the
     * scheme suggested by Andrei
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param tolerance
     * @return size_t
     */
    size_t MEHPForceBalance2::doRemovalAndreisWay(
      Eigen::VectorXd& springPartitions,
      double tolerance)
    {
      size_t numRemovedTotal = this->removeSubfunctionalVertices(2);
      // then, replace f = 2
      this->removeTwofunctionalCrosslinks();
      // and remove all springs that are inactive
      size_t numSpringsRemoved = this->removeInactiveParentEdges(tolerance);

      if (numSpringsRemoved > 0) {
        numRemovedTotal +=
          this->doRemovalAndreisWay(springPartitions, tolerance);
      }
      return numRemovedTotal;
    };

  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
