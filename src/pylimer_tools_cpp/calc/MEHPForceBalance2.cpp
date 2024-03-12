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
    // MARK: Simulation / Optimization Procedures
    //----------------------------------------------------------------

    /**
     * FORCE RELAXATION
     */
    void MEHPForceBalance2::runForceRelaxation(
      long int maxNrOfSteps, // default: 10000
      double xtol,
      const double initialResidualToUse,
      const StructureSimplificationMode simplificationMode,
      const double inactiveRemovalCutoff,
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
          : (0.25 * std::pow(this->universe.getVolume() /
                               this->universe.getNrOfAtoms(),
                             1. / 3.));

      /* array allocation */ double maxDistanceMoved = 0.0;
      size_t indexOfMaxDistanceMoved = 0;
      const double initialResidual =
        (initialResidualToUse > 0.)
          ? initialResidualToUse
          : this->getDisplacementResidualNorm(oneOverSpringPartitionUpperLimit);
      const double minN = this->getNetwork().springsContourLength.minCoeff();
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
                << oneOverSpringPartitionUpperLimit << std::endl;
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
              this->swapSlipLinks(oneOverSpringPartitionUpperLimit);
            } else if (allowSlipLinksToPassEachOther == LinkSwappingMode::ALL) {
              this->swapSlipLinksInclXlinks(oneOverSpringPartitionUpperLimit);
            } else if (allowSlipLinksToPassEachOther ==
                       LinkSwappingMode::ALL_MC) {
              this->moveSlipLinksToTheirBestBranch(
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                false);
            } else if (allowSlipLinksToPassEachOther ==
                       LinkSwappingMode::ALL_MC_CYCLE) {
              this->moveSlipLinksToTheirBestBranch(
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                true);
            } else {
              throw std::invalid_argument(
                "This swapping mode is currently not supported.");
            }
          }
        }
        maxDistanceMoved = 0.0;
        currentResidual = 0.0;

        // place slip-link
        for (size_t link_idx = 0; link_idx < igraph_vcount(&this->graph);
             ++link_idx) {
          if (castToIgraphInt(igraph_cattribute_VAN(
                &this->graph, "type", link_idx)) != this->slipLinkType) {
            continue;
          }
          // std::cout << "Handling " << link_idx << " of " << net.nrOfNodes
          //           << " / " << net.nrOfLinks << std::endl;

          // std::cout << "Still handling " << link_idx << " of " <<
          // net.nrOfNodes
          //           << " / " << net.nrOfLinks << std::endl;
          int innerIterationsDone = 0;
          do {
            double r2 =
              this->updateSpringPartition(link_idx,
                                          oneOverSpringPartitionUpperLimit,
                                          allowSlipLinksToPassEachOther);
            double displacementDone = this->displaceToMeanPosition(
              link_idx, oneOverSpringPartitionUpperLimit);
            innerIterationsDone += 1;
          } while (doInnerIterations && innerIterationsDone < 50);
        }

        intermediateResidual =
          this->getDisplacementResidualNorm(oneOverSpringPartitionUpperLimit);

        // place cross-links
        for (size_t link_idx = 0; link_idx < igraph_vcount(&this->graph);
             ++link_idx) {
          if (castToIgraphInt(igraph_cattribute_VAN(
                &this->graph, "type", link_idx)) != this->crosslinkerType) {
            continue;
          }
          double distanceMoved = this->displaceToMeanPosition(
            link_idx, oneOverSpringPartitionUpperLimit);
          if (distanceMoved > maxDistanceMoved) {
            maxDistanceMoved = distanceMoved;
            indexOfMaxDistanceMoved = link_idx;
          }
          // maxDistanceMoved = std::max(maxDistanceMoved, distanceMoved);
        }

        currentResidual =
          this->getDisplacementResidualNorm(oneOverSpringPartitionUpperLimit);
        iterationsDone += 1;
        if (iterationsDone % 10 == 0) {
          if (simplificationMode ==
                StructureSimplificationMode::INACTIVE_ONLY ||
              simplificationMode == StructureSimplificationMode::ALL_TIM) {
            // std::cout << "Removing inactive cross-links" << std::endl;
            // default tolerance: 0.25*atom's cube length
            size_t nRemoved = this->removeInactiveCrosslinks(removalTolerance);
            if (nRemoved > 0) {
              this->net.isUpToDate = false;
              // std::cout << "Removed " << nRemoved << " inactive springs. "
              //           << std::endl;
            }
          }
          if (simplificationMode == StructureSimplificationMode::X2F_ONLY ||
              simplificationMode == StructureSimplificationMode::ALL_TIM) {
            // std::cout << "Removing 2-f cross-links" << std::endl;
            size_t nRemoved = this->removeTwofunctionalLinks();
            if (nRemoved > 0) {
              this->net.isUpToDate = false;
              // std::cout << "Removed " << nRemoved
              //           << " cross-linkers with f = 2. " << std::endl;
            }
          }
          if (simplificationMode == StructureSimplificationMode::ALL_ANDREI) {
            // std::cout << "Removing cross-links and springs, Andrei's way"
            //           << std::endl;
            size_t nRemoved = this->doRemovalAndreisWay(removalTolerance);
            if (nRemoved > 0) {
              this->net.isUpToDate = false;
            }
          }
          if (simplificationMode !=
              StructureSimplificationMode::NO_SIMPLIFICATION) {
            this->cleanupPrimaryLoopsInStructure();
#ifndef NDEBUG
            this->validateNetwork(this->getNetwork(),
                                  this->currentSpringPartitionsVec);
#endif
          }
        }
        this->handleOutput(iterationsDone);

      } while ((currentResidual / initialResidual > xtol) &&
               (iterationsDone < maxNrOfSteps) &&
               (igraph_ecount(&this->graph) > 0));

      // finish up
      this->closeAllOutputs();

      // query solution & exit reason
      this->exitReason = (iterationsDone >= maxNrOfSteps)
                           ? ExitReason::MAX_STEPS
                           : ExitReason::X_TOLERANCE;
      if (igraph_ecount(&this->graph) == 0) {
        this->exitReason = ExitReason::NO_STEPS_POSSIBLE;
      }
      this->nrOfStepsDone += iterationsDone;
      std::cout << iterationsDone << " steps done. "
                << "Last max distance moved: " << maxDistanceMoved << ". "
                << "Current residual: " << currentResidual << ". "
                << "Initial residual: " << initialResidual << ". " << std::endl;

      assert(this->getNetwork().isUpToDate);
      this->validateNetwork();
      this->currentSpringDistances = this->evaluateSpringDistances();
      this->currentPartialSpringDistances =
        this->evaluatePartialSpringDistances();
    }

    /**
     * @brief Displace one link to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @param linkIdx the idx of the link to displace
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance2::displaceToMeanPosition(
      const size_t linkIdx,
      const double oneOverSpringPartitionUpperLimit)
    {
      INVALIDARG_EXP_IFN(
        linkIdx < igraph_vcount(&this->graph),
        "Invalid link-idx: cannot displace link " + std::to_string(linkIdx) +
          " to mean position, as only " +
          std::to_string(igraph_vcount(&this->graph)) + " links are present.");

      Eigen::Vector3d objectiveDisplacement = Eigen::Vector3d::Zero();
      double objectiveDisplacementContributors = 0.0;

      // fetch the edges involved
      igraph_vector_int_t edgesOfLink;
      igraph_vector_int_init(&edgesOfLink, 2);
      igraph_incident(&this->graph, &edgesOfLink, linkIdx, IGRAPH_ALL);

      for (size_t i = 0; i < igraph_vector_int_size(&edgesOfLink); ++i) {
        igraph_integer_t edgeId = igraph_vector_int_get(&edgesOfLink, i);
        Eigen::Vector3d partialDistance =
          this->evaluatePartialSpringDistanceFrom(edgeId, linkIdx, this->is2D);
        double contourLengthFraction =
          igraph_cattribute_EAN(&this->graph, "partition_fraction", edgeId);
        const double N =
          igraph_cattribute_EAN(&this->graph, "contour_length", edgeId);
        double oneOverContourLengthFraction = 1.0 / (N * contourLengthFraction);
        if (std::isfinite(oneOverContourLengthFraction)) {
          objectiveDisplacement +=
            (partialDistance)*oneOverContourLengthFraction; // /
          // totalDistance.array());
          objectiveDisplacementContributors += oneOverContourLengthFraction;
        }
      }

      igraph_vector_int_destroy(&edgesOfLink);

      // take mean for displacement
      // prevent NaN from division by zero
      Eigen::Vector3d coordsBefore = this->getVertexCoordinates(linkIdx);
      double denominator = (objectiveDisplacementContributors == 0.0
                              ? 0.0
                              : 1. / objectiveDisplacementContributors);
      this->setVertexCoordinates(
        linkIdx, coordsBefore + objectiveDisplacement * denominator);

      double dist = (objectiveDisplacement * denominator).squaredNorm();
      this->net.isUpToDate = false;
      // if (dist > 0.1) {
      //   std::cout << "Moving " << linkIdx << " for " << dist
      //             << " with displacements " << u.segment(3 * linkIdx, 3)[0]
      //             << ", " << u.segment(3 * linkIdx, 3)[1] << ", "
      //             << u.segment(3 * linkIdx, 3)[2] << std::endl;
      //   std::cout << "For objective displacements " <<
      //   objectiveDisplacement[0]
      //             << ", " << objectiveDisplacement[1] << ", "
      //             << objectiveDisplacement[2] << ", for "
      //             << objectiveDisplacementContributors << "." << std::endl;
      // }
      // igraph_cattribute_GAB_set(&this->graph, "is_up_to_date", false);
      return dist;
    }

    //----------------------------------------------------------------
    // MARK: Structural Query
    //----------------------------------------------------------------

    /**
     * @brief List the edges and vertices of one spring, in order
     *
     * @param springIdx
     * @param vertices
     * @param edges
     */
    void MEHPForceBalance2::findEdgesAndVerticesOfSpring(
      size_t springIdx,
      igraph_vector_int_t* vertices,
      igraph_vector_int_t* edges)
    {
      if (this->net.isUpToDate) {
        if (vertices != nullptr) {
          igraph_vector_int_clear(vertices);
          for (size_t i = 0;
               i < this->net.linkIndicesOfSprings[springIdx].size();
               ++i) {
            igraph_vector_int_push_back(
              vertices, this->net.linkIndicesOfSprings[springIdx][i]);
          }
        }
        if (edges != nullptr) {
          igraph_vector_int_clear(edges);
          for (size_t i = 0;
               i < this->net.localToGlobalSpringIndex[springIdx].size();
               ++i) {
            igraph_vector_int_push_back(
              edges, this->net.localToGlobalSpringIndex[springIdx][i]);
          }
        }
        return;
      }
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
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
        if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) == springIdx) {
          igraph_vector_int_push_back(&allEdgesOfParent, i);
          igraph_cattribute_EAN_set(&this->graph, "prev_edge_id", i, i);
          igraph_integer_t from, to;
          igraph_edge(&this->graph, i, &from, &to);
          igraph_cattribute_VAN_set(&this->graph, "prev_vertex_id", from, from);
          igraph_cattribute_VAN_set(&this->graph, "prev_vertex_id", to, to);
        }
      }

      this->findEdgesAndVerticesOfSpring(&allEdgesOfParent, vertices, edges);
      if (vertices != nullptr && edges != nullptr) {
        assert(igraph_vector_int_size(vertices) ==
               igraph_vector_int_size(edges) + 1);
      }

      igraph_vector_destroy(&parentEdges);
      igraph_vector_int_destroy(&allEdgesOfParent);
    }

    /**
     * @brief List the edges and vertices of one spring, in order
     *
     * @param springIdx
     * @param vertices
     * @param edges
     */
    void MEHPForceBalance2::findEdgesAndVerticesOfSpring(
      igraph_vector_int_t* unorderedEdges,
      igraph_vector_int_t* vertices,
      igraph_vector_int_t* edges)
    {
      assert(igraph_vector_int_size(unorderedEdges) > 0);
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
      // figure out, where the rail is currently going through.
      igraph_t subgraph;
      igraph_subgraph_from_edges(
        &this->graph,
        &subgraph,
        igraph_ess_vector(unorderedEdges),
        false); // keep vertices in order to keep vertex ids

      igraph_vector_int_t edgesOnPath;
      igraph_vector_int_init(&edgesOnPath, 0);
      igraph_vector_int_reserve(&edgesOnPath, igraph_ecount(&subgraph));
      igraph_vector_int_t verticesOnPath;
      igraph_vector_int_init(&verticesOnPath, 0);
      igraph_vector_int_reserve(&verticesOnPath, igraph_ecount(&subgraph) + 1);
      RUNTIME_EXP_IFN(
        igraph_eulerian_path(&subgraph, &edgesOnPath, &verticesOnPath) == 0,
        "Could not find expected Eulerian path.");
      assert(igraph_vector_int_size(&edgesOnPath) == igraph_ecount(&subgraph));
      assert(igraph_vector_int_size(unorderedEdges) ==
             igraph_ecount(&subgraph));

      if (edges != nullptr) {
        igraph_vector_int_clear(edges);
        for (size_t i = 0; i < igraph_vector_int_size(&edgesOnPath); ++i) {
          igraph_vector_int_push_back(
            edges,
            castToIgraphInt(
              igraph_cattribute_EAN(&subgraph,
                                    "prev_edge_id",
                                    igraph_vector_int_get(&edgesOnPath, i))));
        }
        assert(igraph_vector_int_size(&edgesOnPath) ==
               igraph_vector_int_size(edges));
      }
      if (vertices != nullptr) {
        igraph_vector_int_clear(vertices);
        for (size_t i = 0; i < igraph_vector_int_size(&verticesOnPath); ++i) {
          igraph_vector_int_push_back(
            vertices, igraph_vector_int_get(&verticesOnPath, i));
        }
        assert(igraph_vector_int_size(&verticesOnPath) ==
               igraph_vector_int_size(vertices));
      }
      igraph_integer_t vertex0Id = igraph_vector_int_get(&verticesOnPath, 0);
      igraph_integer_t vertexEndId = igraph_vector_int_get(
        &verticesOnPath, igraph_vector_int_size(&verticesOnPath) - 1);
      //  make sure the ordering is consistent between different calls to this
      //  function
      if (vertex0Id > vertexEndId) {
        if (vertices != nullptr) {
          igraph_vector_int_reverse(vertices);
        }
        if (edges != nullptr) {
          igraph_vector_int_reverse(edges);
        }
      }

      if (vertices != nullptr && edges != nullptr) {
        assert(igraph_vector_int_size(vertices) >
               igraph_vector_int_size(edges));
      }

      igraph_destroy(&subgraph);
      igraph_vector_int_destroy(&edgesOnPath);
      igraph_vector_int_destroy(&verticesOnPath);
    }

    /**
     * @brief Find the index of a partial spring, given the fractn of the
     * total spring to traverse
     *
     * @param springIdx
     * @param alpha
     * @param fractionTillThere the sum of the fractions before reaching the
     * searched index
     * @return igraph_integer_t
     */
    igraph_integer_t MEHPForceBalance2::findPartialSpringByFraction(
      size_t springIdx,
      double alpha,
      double& fractionTillThere)
    {
      igraph_vector_int_t edges;
      igraph_vector_int_init(&edges, 0);

      this->findEdgesAndVerticesOfSpring(springIdx, nullptr, &edges);
      assert(igraph_vector_int_size(&edges) > 0);

      double currentAlpha = 0.0;
      igraph_integer_t edgeId = igraph_vector_int_get(&edges, 0);
      for (size_t i = 0; i < igraph_vector_int_size(&edges); ++i) {
        edgeId = igraph_vector_int_get(&edges, i);
        igraph_real_t thisEdgeFraction =
          igraph_cattribute_EAN(&this->graph, "partition_fraction", edgeId);
        assert(APPROX_WITHIN(thisEdgeFraction, 0.0, 1.0, 1e-5));
        currentAlpha += thisEdgeFraction;

        if (currentAlpha > alpha) {
          fractionTillThere = currentAlpha - thisEdgeFraction;
          break;
        }
      }
      assert(APPROX_WITHIN(currentAlpha, 0.0, 1.0, 1e-5));

      igraph_vector_int_destroy(&edges);

      return edgeId;
    }

    /**
     * @brief Given a vertex id and a rail edge, returns the other two edges
     * that are not part of the rail
     */
    std::vector<igraph_integer_t> MEHPForceBalance2::getOffRailConnectedEdgeIds(
      igraph_integer_t vertexId,
      igraph_integer_t railEdgeId,
      igraph_integer_t avoidEdgeId)
    {
      INVALIDARG_EXP_IFN(igraph_cattribute_VAN(&this->graph,
                                               "type",
                                               vertexId) == this->slipLinkType,
                         "Can only search for rail around slip-links");

      // fetch the edges involved
      igraph_integer_t otherRailEdge =
        this->getOtherRailEdgeId(vertexId, railEdgeId, avoidEdgeId);
      igraph_es_t selector;
      igraph_es_incident(&selector, vertexId, IGRAPH_ALL);
      igraph_eit_t iterator;
      igraph_eit_create(&this->graph, selector, &iterator);
      std::vector<igraph_integer_t> results;
      results.reserve(IGRAPH_EIT_SIZE(iterator) - 2);
      while (!IGRAPH_EIT_END(iterator)) {
        if ((IGRAPH_EIT_GET(iterator) != railEdgeId) &&
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
    igraph_integer_t MEHPForceBalance2::getOtherRailEdgeId(
      igraph_integer_t vertexId,
      igraph_integer_t railEdgeId,
      igraph_integer_t avoidEdgeId)
    {
      INVALIDARG_EXP_IFN(igraph_cattribute_VAN(&this->graph,
                                               "type",
                                               vertexId) == this->slipLinkType,
                         "Can only search for rail around slip-links");
      igraph_integer_t fromRail, toRail;
      igraph_edge(&this->graph, railEdgeId, &fromRail, &toRail);
      INVALIDARG_EXP_IFN(fromRail == vertexId || toRail == vertexId,
                         "Vertex must be part of the rail");

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
        igraph_vector_int_set(&edgeIds, current_idx, IGRAPH_EIT_GET(iterator));
        if (IGRAPH_EIT_GET(iterator) == railEdgeId) {
          rail_idx = current_idx;
        }
        IGRAPH_EIT_NEXT(iterator);
        current_idx += 1;
      }
      assert(rail_idx >= 0);
      assert(current_idx == igraph_vector_int_size(&edgeIds));

      igraph_vector_t parents;
      igraph_vector_init(&parents, igraph_vector_int_size(&edgeIds));
      igraph_cattribute_EANV(&this->graph, "parent_edge", selector, &parents);

      igraph_integer_t railParent =
        castToIgraphInt(igraph_vector_get(&parents, rail_idx));
      assert(railParent ==
             igraph_cattribute_EAN(&this->graph, "parent_edge", railEdgeId));

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
      igraph_vector_int_destroy(&edgeIds);

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
      igraph_subgraph_from_edges(&this->graph,
                                 &subgraph,
                                 igraph_ess_vector(&allEdgesOfParent),
                                 false); // keep vertices for ids

      igraph_vector_int_t edgesOnPath;
      igraph_vector_int_init(&edgesOnPath, igraph_ecount(&subgraph));
      igraph_vector_int_t verticesOnPath;
      igraph_vector_int_init(&verticesOnPath, igraph_ecount(&subgraph) + 1);
      igraph_eulerian_path(&subgraph, &edgesOnPath, &verticesOnPath);
      assert(igraph_vector_int_size(&edgesOnPath) == igraph_ecount(&subgraph));
      if (igraph_vector_int_get(&verticesOnPath, 0) >
          igraph_vector_int_get(&verticesOnPath,
                                igraph_vector_int_size(&verticesOnPath) - 1)) {
        // make sure we are always looking at the chain in the same direction
        igraph_vector_int_reverse(&verticesOnPath);
        igraph_vector_int_reverse(&edgesOnPath);
      }
      this->printVerticesOnPath(&verticesOnPath);

      size_t result = 0;
      bool foundResult = false;
      for (size_t pathStepIdx = 0;
           pathStepIdx < igraph_vector_int_size(&edgesOnPath);
           ++pathStepIdx) {
        igraph_integer_t thisEdgeId =
          igraph_vector_int_get(&edgesOnPath, pathStepIdx);
        if (castToIgraphInt(igraph_cattribute_EAN(
              &subgraph, "prev_edge_id", thisEdgeId)) == railEdgeId) {
          igraph_integer_t from, to;
          igraph_edge(&subgraph, thisEdgeId, &from, &to);
          igraph_integer_t fromV = castToIgraphInt(
            igraph_cattribute_VAN(&subgraph, "prev_vertex_id", from));
          igraph_integer_t toV = castToIgraphInt(
            igraph_cattribute_VAN(&subgraph, "prev_vertex_id", to));
          assert(from == fromRail && to == toRail);
          igraph_integer_t prevEdgeId = -1;
          igraph_integer_t nextEdgeId = -1;
          // yay, found the rail.
          // now: on either side of this is the link we talk about.
          // We are interested in the edge on the other
          // side!
          igraph_integer_t prevfrom = -1, prevto = -1;
          igraph_integer_t nextfrom = -1, nextto = -1;
          if (pathStepIdx > 0) {
            prevEdgeId = igraph_vector_int_get(&edgesOnPath, pathStepIdx - 1);
            igraph_edge(&subgraph, prevEdgeId, &prevfrom, &prevto);
          }
          if (pathStepIdx < igraph_ecount(&subgraph) - 1) {
            nextEdgeId = igraph_vector_int_get(&edgesOnPath, pathStepIdx + 1);
            igraph_edge(&subgraph, nextEdgeId, &nextfrom, &nextto);
          }
          if (((prevfrom == vertexId) || (prevto == vertexId)) &&
              ((nextfrom == vertexId) || (nextto == vertexId))) {
            std::cerr << "Handling special case for partial spring "
                      << railEdgeId << std::endl;
            igraph_integer_t otherVertexId =
              vertexId == fromRail ? toRail : fromRail;
            igraph_integer_t prevPrevEdgeId = castToIgraphInt(
              igraph_cattribute_EAN(&subgraph, "prev_edge_id", prevEdgeId));
            igraph_integer_t prevNextEdgeId = castToIgraphInt(
              igraph_cattribute_EAN(&subgraph, "prev_edge_id", nextEdgeId));
            // both the previous and next edges are on the rail.
            // this can happen if we have a secondary loop.
            // momentary resolution (works for many examples, at least: choose
            // the less familiar one)
            if (((prevfrom == otherVertexId) || (prevto == otherVertexId)) ||
                (prevPrevEdgeId == avoidEdgeId)) {
              result = prevNextEdgeId;
              foundResult = true;
            } else {
              result = prevPrevEdgeId;
              foundResult = true;
            }
            break;
          } else if ((prevfrom == vertexId) || (prevto == vertexId)) {
            // i-1 is on-rail
            result = castToIgraphInt(
              igraph_cattribute_EAN(&subgraph, "prev_edge_id", prevEdgeId));
            foundResult = true;
            break;
          } else {
            assert((nextfrom == vertexId) || (nextto == vertexId));
            // i+1 is on-rail
            result = castToIgraphInt(
              igraph_cattribute_EAN(&subgraph, "prev_edge_id", nextEdgeId));
            foundResult = true;
            break;
          }
        }
      }

      igraph_destroy(&subgraph);
      igraph_vector_destroy(&parentEdges);
      igraph_vector_int_destroy(&allEdgesOfParent);
      igraph_vector_int_destroy(&edgesOnPath);
      igraph_vector_int_destroy(&verticesOnPath);

      RUNTIME_EXP_IFN(foundResult, "Did not find other edge.");
      return result;
    }

    //----------------------------------------------------------------
    // MARK: Structural Adjustments
    //----------------------------------------------------------------

    /**
     * @brief Combine two partial springs to be only one
     *
     * @param edge1Id
     * @param edge2Id
     */
    void MEHPForceBalance2::combinePartialSprings(
      const igraph_integer_t edge1Id,
      const igraph_integer_t edge2Id,
      const igraph_integer_t centralLink)
    {
      if (edge1Id > edge2Id) {
        return this->combinePartialSprings(edge2Id, edge1Id, centralLink);
      }
      assert(edge1Id < edge2Id);
      igraph_integer_t from1, to1;
      igraph_edge(&this->graph, edge1Id, &from1, &to1);
      igraph_integer_t from2, to2;
      igraph_edge(&this->graph, edge2Id, &from2, &to2);

      assert(from1 == centralLink || to1 == centralLink);
      assert(from2 == centralLink || to2 == centralLink);

      igraph_integer_t parent1 = castToIgraphInt(
        igraph_cattribute_EAN(&this->graph, "parent_edge", edge1Id));
      igraph_integer_t parent2 = castToIgraphInt(
        igraph_cattribute_EAN(&this->graph, "parent_edge", edge2Id));

#ifndef NDEBUG
      this->validateIgraphSpring(parent1);
      this->validateIgraphSpring(parent2);
#endif

      igraph_integer_t newTo = to2 == centralLink ? from2 : to2;
      igraph_integer_t newFrom = to1 == centralLink ? from1 : to1;

      size_t newEdgeId = this->createEdge(
        newFrom,
        newTo,
        std::min(parent1, parent2),
        igraph_cattribute_EAN(&this->graph, "partition_fraction", edge1Id) +
          igraph_cattribute_EAN(&this->graph, "partition_fraction", edge2Id),
        igraph_cattribute_EAN(&this->graph, "contour_length", edge1Id),
        this->getBondBoxOffsetForEdgeFrom(edge2Id, centralLink) +
          this->getBondBoxOffsetForEdgeTo(edge1Id, centralLink));
      // #ifndef NDEBUG
      //       std::cout << "Combining " << edge1Id << " (parent " << parent1
      //                 << ") with " << edge2Id << " (" << parent2 << ") to "
      //                 << newEdgeId << std::endl;
      // #endif

      if (parent1 == parent2) {
        this->debugParentEdge(parent1);
      }

      // need this iff edge2 is the last remaining of its parent
      double contourLengthBefore =
        igraph_cattribute_EAN(&this->graph, "contour_length", edge2Id);

      // delete the old edges – order matters!
      igraph_delete_edges(&this->graph, igraph_ess_1(edge2Id));
      igraph_delete_edges(&this->graph, igraph_ess_1(edge1Id));

      // merge edge properties
      if (parent1 != parent2) {
        // std::cout << "Merging " << parent1 << " & " << parent2 << std::endl;
        assert(castToIgraphInt(igraph_cattribute_VAN(
                 &this->graph, "type", centralLink)) == this->crosslinkerType);
        // two different "parent" springs -> need to recalculate some stuff
        // probably only if link is cross-link
        this->combineParentSprings(std::max(parent1, parent2),
                                   std::min(parent1, parent2),
                                   contourLengthBefore);
      }

#ifndef NDEBUG
      this->validateIgraphSpring(parent1);
#endif
      this->net.isUpToDate = false;
    }

    /**
     * @brief Remove a slip-link and combine the two edges corresponding to
     * the rail
     *
     * @param vertexId
     * @param railEdgeId
     */
    void MEHPForceBalance2::unlinkSlipLinkFromRail(
      const igraph_integer_t vertexId,
      const igraph_integer_t railEdgeId)
    {
      igraph_integer_t otherRailEdge =
        this->getOtherRailEdgeId(vertexId, railEdgeId);

      this->combinePartialSprings(railEdgeId, otherRailEdge, vertexId);
    }

    /**
     * @brief Inserts the given slip-link into a partial spring
     *
     * @param vertexId the slip-link to insert into the spring
     * @param railEdgeId the spring to be halfed
     */
    void MEHPForceBalance2::insertSlipLinkIntoRail(
      const igraph_integer_t vertexId,
      const igraph_integer_t railEdgeId)
    {
      igraph_integer_t from, to;
      igraph_edge(&this->graph, railEdgeId, &from, &to);

      igraph_integer_t newEdge1 = this->createEdge(
        vertexId,
        from,
        igraph_cattribute_EAN(&this->graph, "parent_edge", railEdgeId),
        0.5 *
          igraph_cattribute_EAN(&this->graph, "partition_fraction", railEdgeId),
        igraph_cattribute_EAN(&this->graph, "contour_length", railEdgeId),
        Eigen::Vector3d::Zero());
      igraph_integer_t newEdge2 = this->createEdge(
        vertexId,
        to,
        igraph_cattribute_EAN(&this->graph, "parent_edge", railEdgeId),
        0.5 *
          igraph_cattribute_EAN(&this->graph, "partition_fraction", railEdgeId),
        igraph_cattribute_EAN(&this->graph, "contour_length", railEdgeId),
        this->getBondBoxOffsetForEdge(railEdgeId));

      igraph_delete_edges(&this->graph, igraph_ess_1(railEdgeId));
      this->net.isUpToDate = false;
    }

    void MEHPForceBalance2::moveSlipLinkFromRailToRail(
      const igraph_integer_t vertexId,
      const igraph_integer_t sourceRailEdgeId,
      const igraph_integer_t targetRailEdgeId)
    {
      // the difficulty in this implementation is that
      // the edge ids change

      igraph_integer_t otherRailEdge =
        this->getOtherRailEdgeId(vertexId, sourceRailEdgeId);

      // replace the target edge with two new edges, via the vertex
      igraph_integer_t from, to;
      igraph_edge(&this->graph, targetRailEdgeId, &from, &to);

      igraph_integer_t newEdge1 = this->createEdge(
        vertexId,
        from,
        igraph_cattribute_EAN(&this->graph, "parent_edge", targetRailEdgeId),
        0.5 * igraph_cattribute_EAN(
                &this->graph, "partition_fraction", targetRailEdgeId),
        igraph_cattribute_EAN(&this->graph, "contour_length", targetRailEdgeId),
        Eigen::Vector3d::Zero());
      igraph_integer_t newEdge2 = this->createEdge(
        vertexId,
        to,
        igraph_cattribute_EAN(&this->graph, "parent_edge", targetRailEdgeId),
        0.5 * igraph_cattribute_EAN(
                &this->graph, "partition_fraction", targetRailEdgeId),
        igraph_cattribute_EAN(&this->graph, "contour_length", targetRailEdgeId),
        this->getBondBoxOffsetForEdge(targetRailEdgeId));

      // delete the old one
      igraph_delete_edges(&this->graph, igraph_ess_1(targetRailEdgeId));

      // and re-weld the source
      this->combinePartialSprings(
        sourceRailEdgeId - (targetRailEdgeId < sourceRailEdgeId ? 1 : 0),
        otherRailEdge - (targetRailEdgeId < otherRailEdge ? 1 : 0),
        vertexId);

      this->net.isUpToDate = false;
    }

    /**
     * @brief Remove a spring (and all its parts, incl. slip-links) from the
     * structures
     *
     * NOTE: this may result in slip-links with f = 2.
     * They are not automatically removed in order to preserve vertex
     * iterations.
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance2::removeParentSpring(const size_t springIdx)
    {
      if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
        this->updateGraph();
      }

      igraph_vector_t parentEdges;
      igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
      igraph_cattribute_EANV(&this->graph,
                             "parent_edge",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &parentEdges);

      std::vector<igraph_integer_t> edgeIdsToRemove;
      edgeIdsToRemove.reserve(4);
      for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        if (igraph_vector_get(&parentEdges, i) == springIdx) {
          edgeIdsToRemove.push_back(i);
        }
      }

      igraph_vector_destroy(&parentEdges);
      // actually remove all edges
      this->removePartialSprings(edgeIdsToRemove);
      // as the above may lead to slip-links with f = 2,
      // we want to remove those as well
      // this->removeTwofunctionalLinks();
    };

    /**
     * @brief Remove a set of edges from the graph
     *
     * @param edgeIdsToRemove
     */
    void MEHPForceBalance2::removePartialSprings(
      std::vector<igraph_integer_t>& edgeIdsToRemove)
    {
      igraph_vector_int_t edgeIdsToRemoveVec;
      igraph_vector_int_init(&edgeIdsToRemoveVec, edgeIdsToRemove.size());
      pylimer_tools::utils::StdVectorToIgraphVectorT(edgeIdsToRemove,
                                                     &edgeIdsToRemoveVec);
      igraph_delete_edges(&this->graph, igraph_ess_vector(&edgeIdsToRemoveVec));
      igraph_vector_int_destroy(&edgeIdsToRemoveVec);
      // remove vertices that "got lost"
      // this->removeOrphanedVertices();

      this->net.isUpToDate = false;
    }

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
     * Invalidates both vertex and edge iterators
     */
    void MEHPForceBalance2::remove2fLink(const size_t linkIdx)
    {
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

      // fetch the edges involved
      igraph_vector_int_t edgesOfLink;
      igraph_vector_int_init(&edgesOfLink, 2);
      igraph_incident(&this->graph, &edgesOfLink, linkIdx, IGRAPH_ALL);

      if (igraph_vector_int_size(&edgesOfLink) == 0) {
        REQUIRE_IGRAPH_SUCCESS(
          igraph_delete_vertices(&this->graph, igraph_vss_1(linkIdx)));
        igraph_vector_int_destroy(&edgesOfLink);
        this->net.isUpToDate = false;
        return;
      }

      if (igraph_vector_int_size(&edgesOfLink) == 1) {
        igraph_integer_t from, to;
        igraph_edge(
          &this->graph, igraph_vector_int_get(&edgesOfLink, 0), &from, &to);
        if (from != to) {
          throw std::runtime_error("Only self-loops are acceptable 'f = 1'");
        }

        REQUIRE_IGRAPH_SUCCESS(
          igraph_delete_vertices(&this->graph, igraph_vss_1(linkIdx)));
        igraph_vector_int_destroy(&edgesOfLink);
        this->net.isUpToDate = false;
        return;
      }

      if (igraph_vector_int_size(&edgesOfLink) != 2) {
        throw std::runtime_error(
          "Expect f = 2 to have 2 edges, got " +
          std::to_string(igraph_vector_int_size(&edgesOfLink)) + ".");
      }

      size_t removedEdge1 = std::min(igraph_vector_int_get(&edgesOfLink, 0),
                                     igraph_vector_int_get(&edgesOfLink, 1));
      size_t removedEdge2 = std::max(igraph_vector_int_get(&edgesOfLink, 0),
                                     igraph_vector_int_get(&edgesOfLink, 1));

      igraph_vector_int_destroy(&edgesOfLink);

      this->combinePartialSprings(removedEdge1, removedEdge2, linkIdx);

      igraph_delete_vertices(&this->graph, igraph_vss_1(linkIdx));
      this->net.isUpToDate = false;
    }

    /**
     * @brief Remove a certain, 3-functional link from the structures,
     * combining the two strands
     *
     * Invalidates both vertex and edge iterators
     */
    void MEHPForceBalance2::remove3fLink(const size_t linkIdx)
    {
      // before combining the two strands, need to know which ones.
      // easiest case: the
      igraph_es_t selector;
      igraph_es_incident(&selector, linkIdx, IGRAPH_ALL);
      igraph_eit_t iterator;
      igraph_eit_create(&this->graph, selector, &iterator);
      std::vector<size_t> parents;
      std::vector<igraph_integer_t> edgeIds;
      RUNTIME_EXP_IFN(IGRAPH_EIT_SIZE(iterator) == 3, "Expect f = 3");
      parents.reserve(IGRAPH_EIT_SIZE(iterator));
      edgeIds.reserve(IGRAPH_EIT_SIZE(iterator));
      while (!IGRAPH_EIT_END(iterator)) {
        parents.push_back(igraph_cattribute_EAN(
          &this->graph, "parent_edge", IGRAPH_EIT_GET(iterator)));
        edgeIds.push_back(IGRAPH_EIT_GET(iterator));
        IGRAPH_EIT_NEXT(iterator);
      }

      igraph_es_destroy(&selector);
      igraph_eit_destroy(&iterator);

      igraph_integer_t edgeToRemove;
      if (std::all_of(parents.begin(), parents.end(), [&](size_t i) {
            return i == parents[0];
          })) {
        // all equal -> this spring will be removed as a whole
        // CAUTION: any iteration will be invalidated
        this->removeParentSpring(parents[0]);
        return;
      }
      // not all equal -> want to re-connect the two that belong to the same
      // parent after removing the link
      if (parents[0] == parents[1]) {
        edgeToRemove = edgeIds[2];
      } else if (parents[0] == parents[2]) {
        edgeToRemove = edgeIds[1];
      } else {
        assert(parents[1] == parents[2]);
        edgeToRemove = edgeIds[0];
      }
      igraph_delete_edges(&this->graph, igraph_ess_1(edgeToRemove));
      this->remove2fLink(linkIdx);
    };

    /**
     * @brief marks a certain "parent" spring as non-existing
     */
    void MEHPForceBalance2::combineParentSprings(size_t springIdxBefore,
                                                 size_t springIdxNow,
                                                 double contourLengthBefore)
    {
      if (springIdxBefore == springIdxNow) {
        throw std::invalid_argument(
          "Will only replace higher with lower spring idx");
      }
      assert(springIdxBefore > springIdxNow);
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

      // #ifndef NDEBUG
      //       std::cout << "Combining parent " << springIdxBefore << " into "
      //                 << springIdxNow << std::endl;
      // #endif

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
      igraph_vector_t contourLengths;
      igraph_vector_init(&contourLengths, this->net.nrOfPartialSprings);
      igraph_cattribute_EANV(&this->graph,
                             "contour_length",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &contourLengths);

      // we also need to normalize the spring partition,
      // since this function may be called when the sum is incorrect
      double springContourLengthNow = -1.;
      double springParitionBefore = 0.;
      double springContourLengthBefore = -1.;
      double springPartitionNow = 0.;
      for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) ==
            springIdxBefore) {
          springContourLengthBefore = igraph_vector_get(&contourLengths, i);
          springParitionBefore += igraph_vector_get(&partitionFraction, i);
        }
        if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) ==
            springIdxNow) {
          springContourLengthNow = igraph_vector_get(&contourLengths, i);
          springPartitionNow += igraph_vector_get(&partitionFraction, i);
        }
      }
      assert(springContourLengthNow > 0);
      bool otherEndExist = true;
      if (springContourLengthBefore < 0.) {
        otherEndExist = false;
        assert(contourLengthBefore > 0.);
        springContourLengthBefore = contourLengthBefore;
      }
      assert(springContourLengthBefore > 0);
      double springContourLengthNew =
        springContourLengthNow + springContourLengthBefore;
      double scalingFactorRemoved =
        (springContourLengthBefore) /
        (springContourLengthNew * springParitionBefore);
      double scalingFactorNow =
        otherEndExist ? ((springContourLengthNow) /
                         (springContourLengthNew * springPartitionNow))
                      : (1. / springPartitionNow);

      for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        // update the new value of the parent edge
        if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) ==
            springIdxBefore) {
          // ...if the edge was merged,
          igraph_vector_set(&parentEdges, i, springIdxNow);
          igraph_vector_set(&partitionFraction,
                            i,
                            igraph_vector_get(&partitionFraction, i) *
                              scalingFactorRemoved);
          igraph_vector_set(&contourLengths, i, springContourLengthNew);
        } else if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) >
                   springIdxBefore) {
          // if the edge id has to change
          igraph_vector_set(
            &parentEdges, i, igraph_vector_get(&parentEdges, i) - 1);
        } else if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) ==
                   springIdxNow) {
          // ...the edge was merged into
          // update partition fraction
          igraph_vector_set(&partitionFraction,
                            i,
                            igraph_vector_get(&partitionFraction, i) *
                              scalingFactorNow);
          igraph_vector_set(&contourLengths, i, springContourLengthNew);
        }
      }

      igraph_cattribute_EAN_setv(&this->graph, "parent_edge", &parentEdges);
      igraph_cattribute_EAN_setv(
        &this->graph, "contour_length", &contourLengths);
      igraph_cattribute_EAN_setv(
        &this->graph, "partition_fraction", &partitionFraction);

      igraph_vector_destroy(&parentEdges);
      igraph_vector_destroy(&contourLengths);
      igraph_vector_destroy(&partitionFraction);

#ifndef NDEBUG
      this->validateIgraphSpring(springIdxNow);
#endif

      this->net.isUpToDate = false;
    }

    //----------------------------------------------------------------
    // MARK: Concrete structure modification
    //----------------------------------------------------------------

    /**
     * @brief Replace the two springs traversinga a two-functional cross-links
     * with a single spring
     *
     * @param net
     * @param displacements
     * @param springPartitions
     */
    size_t MEHPForceBalance2::removeTwofunctionalLinks()
    {
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

      size_t numRemoved = 0;
      long int vcount = igraph_vcount(&this->graph);
      for (long int i = vcount - 1; i >= 0; --i) {
        igraph_integer_t degree;
        // count self-loops; they should be removed in another way
        igraph_degree_1(&this->graph, &degree, i, IGRAPH_ALL, true);
        if (degree == 2) {
          // remove this link
          numRemoved += 1;
          this->remove2fLink(i);
        }
      }

      return numRemoved;
    };

    /**
     * @brief Remove double listed springs from cross-links (if they have
     * length 0)
     *
     * @param net
     * @return size_t the nr of removed edges
     */
    size_t MEHPForceBalance2::cleanupPrimaryLoopsInStructure()
    {
      igraph_vector_int_t allEdges;
      igraph_vector_int_init(&allEdges, igraph_ecount(&this->graph) * 2);
      if (igraph_edges(
            &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges)) {
        throw std::runtime_error("Failed to get all edges");
      }

      std::vector<igraph_integer_t> edgesToRemove;
      for (igraph_integer_t i = 0; i < igraph_ecount(&this->graph); ++i) {
        if (igraph_vector_int_get(&allEdges, 2 * i + 0) ==
            igraph_vector_int_get(&allEdges, 2 * i + 1)) {
          // check if this primary loop has length 0 -> remove
          if (this->getBondBoxOffsetForEdge(i).norm() <
              1e-3 * this->universe.getBox().getL().minCoeff()) {
            edgesToRemove.push_back(i);
          }
        }
      }

      this->removePartialSprings(edgesToRemove);
      if (edgesToRemove.size() > 0) {
        this->renumberParentSprings();
        // since we remove partial springs without adjusting any partitions,
        // we need to do that now
        this->renormalizePartitions();
      }
      return edgesToRemove.size();
    };

    /**
     * @brief Remove all "parent" springs that have no active "children"
     *
     * @param tolerance the acceptance tolerance, partial springs longer than
     * this are active
     */
    size_t MEHPForceBalance2::removeInactiveParentEdges(double tolerance)
    {
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

      igraph_vector_t parentEdges;
      igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
      igraph_cattribute_EANV(&this->graph,
                             "parent_edge",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &parentEdges);

      std::unordered_map<igraph_integer_t, bool> isRemovalCandidate;
      std::unordered_map<igraph_integer_t, std::vector<size_t>>
        parentsPartialSprings;
      for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        igraph_integer_t parentEdgeId =
          castToIgraphInt(igraph_vector_get(&parentEdges, i));
        isRemovalCandidate[parentEdgeId] = true;
        std::vector<size_t> vec;
        parentsPartialSprings[parentEdgeId] = vec;
      }
      for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        igraph_integer_t parentEdgeId =
          castToIgraphInt(igraph_vector_get(&parentEdges, i));
        if (isRemovalCandidate.at(parentEdgeId)) {
          isRemovalCandidate[parentEdgeId] =
            this->computeEdgeDistance(i).squaredNorm() <= tolerance;
        }
        parentsPartialSprings[parentEdgeId].push_back(i);
      }

      igraph_vector_destroy(&parentEdges);

      std::vector<igraph_integer_t> edgeIdsToRemove;
      for (auto& [key, value] : isRemovalCandidate) {
        if (value) {
          edgeIdsToRemove.insert(edgeIdsToRemove.end(),
                                 parentsPartialSprings[key].begin(),
                                 parentsPartialSprings[key].end());
        }
      }

      this->removePartialSprings(edgeIdsToRemove);
      size_t numRemoved = edgeIdsToRemove.size();

      if (numRemoved > 0) {
        this->renumberParentSprings();
      }

      return numRemoved;
    }

    /**
     * @brief remove all vertices that don't have any connections
     *
     * @return size_t the nr of vertices removed
     */
    size_t MEHPForceBalance2::removeOrphanedVertices()
    {
      if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
        this->updateGraph();
      }

      igraph_vector_int_t degrees;
      igraph_vector_int_init(&degrees, igraph_vcount(&this->graph));
      igraph_degree(&this->graph, &degrees, igraph_vss_all(), IGRAPH_ALL, true);

      igraph_vector_int_t vertexIds;
      igraph_vector_int_init(&vertexIds, 0);

      for (size_t i = 0; i < igraph_vector_int_size(&degrees); ++i) {
        if (igraph_vector_int_get(&degrees, i) == 0) {
          igraph_vector_int_push_back(&vertexIds, i);
        }
      }
      igraph_vector_int_destroy(&degrees);

      size_t numVertexIds = igraph_vector_int_size(&vertexIds);

      if (numVertexIds == 0) {
        igraph_vector_int_destroy(&vertexIds);
        return 0;
      }

      igraph_delete_vertices(&this->graph, igraph_vss_vector(&vertexIds));
      igraph_vector_int_destroy(&vertexIds);

      this->net.isUpToDate = false;
      return numVertexIds;
    }

    /**
     * @brief Remove chains with links with f = 1.
     *
     * NOTE: this function may result in even more chains with links with f
     * = 1.
     *
     * @return size_t
     */
    size_t MEHPForceBalance2::removeDanglingChains()
    {
      size_t numRemoved = 0;
      for (long int i = igraph_vcount(&this->graph) - 1; i >= 0; --i) {
        if (this->getVertexDegree(i) == 1) {
          igraph_vector_int_t edgesOfVertex;
          igraph_vector_int_init(&edgesOfVertex, 1);
          igraph_incident(&this->graph, &edgesOfVertex, i, IGRAPH_ALL);
          assert(igraph_vector_int_size(&edgesOfVertex) == 1);

          this->removeParentSpring(castToIgraphInt(
            igraph_cattribute_EAN(&this->graph,
                                  "parent_edge",
                                  igraph_vector_int_get(&edgesOfVertex, 0))));
          numRemoved += 1;
        }
      }

      // since some springs have been removed, we have to re-index the others
      if (numRemoved > 0) {
        this->renumberParentSprings();
      }

      return numRemoved;
    }

    /**
     * @brief Remove all vertices (incl. edges!) with a functionality < 3 for
     * cross-links, < 4 for slip-links
     *
     * @param minCrosslinkFunctionalityToBeKept
     * @return size_t the number of removed vertices
     */
    size_t MEHPForceBalance2::removeSubfunctionalVertices()
    {
      size_t numRemovedTotal = 0;
      size_t numRemovedInIteration = 0;
      size_t primaryLoopsRemovedTotal = this->cleanupPrimaryLoopsInStructure();
      // std::cout << "Removed " << primaryLoopsRemovedTotal << " primary loops"
      //           << std::endl;
      size_t primaryLoopsRemovedInIteration = 0;
      do {
        this->removeDanglingChains();
        numRemovedInIteration = this->removeOrphanedVertices();

#ifndef VERBOSE_DEBUG
        this->validateIgraphSprings();
#endif

        // std::cout << "Removed " << numRemovedInIteration
        //           << " vertices with degree < 2" << std::endl;

        primaryLoopsRemovedInIteration = this->cleanupPrimaryLoopsInStructure();
#ifndef VERBOSE_DEBUG
        this->validateIgraphSprings();
#endif
        // std::cout << "Removed " << primaryLoopsRemovedInIteration
        //           << " primary loops" << std::endl;

        igraph_vector_t types;
        igraph_vector_init(&types, 0);
        igraph_cattribute_VANV(&this->graph, "type", igraph_vss_all(), &types);

        // f = 1 and below have been removed
        // -> cleanup remaining f = 2 and f = 3
        if (numRemovedInIteration > 0 || primaryLoopsRemovedInIteration > 0) {
          for (long int i = igraph_vcount(&this->graph) - 1; i >= 0; --i) {
            igraph_integer_t degree = this->getVertexDegree(i);
            if (degree == 2) {
              this->remove2fLink(i);
              numRemovedInIteration += 1;
            } else if (degree == 3 && castToIgraphInt(igraph_vector_get(
                                        &types, i)) == this->slipLinkType) {
              std::cerr << "Found slip-link with f = 3, unexpected!"
                        << std::endl;
              this->remove3fLink(i);
              numRemovedInIteration += 1;
            }
          }
        }

#ifndef VERBOSE_DEBUG
        this->validateIgraphSprings();
#endif
        // std::cout << "Removed " << numRemovedInIteration
        //           << " vertices with degree == 2 or slip-Links degree == 3"
        //           << std::endl;

        igraph_vector_destroy(&types);
        numRemovedTotal += numRemovedInIteration;
        primaryLoopsRemovedInIteration +=
          this->cleanupPrimaryLoopsInStructure();
        primaryLoopsRemovedTotal += primaryLoopsRemovedInIteration;
#ifndef VERBOSE_DEBUG
        this->validateIgraphSprings();
#endif
        // std::cout << "Removed " << primaryLoopsRemovedInIteration
        //           << " primary loops" << std::endl;
      } while (numRemovedInIteration > 0 || primaryLoopsRemovedInIteration > 0);

      if (numRemovedTotal > 0 || primaryLoopsRemovedTotal > 0) {
        this->renumberParentSprings();
        this->net.isUpToDate = false;
      }
#ifndef NDEBUG
      this->validateIgraphSprings();
#endif
      return numRemovedTotal;
    }

    /**
     * @brief Remove chains that have two otherwise not connected ends
     *
     * @return size_t the nr of chains removed
     */
    size_t MEHPForceBalance2::removeFreeChains()
    {
      igraph_vector_int_t degrees;
      igraph_vector_int_init(&degrees, igraph_vcount(&this->graph));
      igraph_degree(&this->graph, &degrees, igraph_vss_all(), IGRAPH_ALL, true);

      igraph_vector_int_t allEdges;
      igraph_vector_int_init(&allEdges, igraph_ecount(&this->graph) * 2);
      if (igraph_edges(
            &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges)) {
        throw std::runtime_error("Failed to get all edges");
      }

      igraph_vector_int_t verticesToDelete;
      igraph_vector_int_init(&verticesToDelete, 0);
      for (size_t i = 0; i < igraph_ecount(&this->graph); ++i) {
        igraph_integer_t from = igraph_vector_int_get(&allEdges, 2 * i + 0);
        igraph_integer_t to = igraph_vector_int_get(&allEdges, 2 * i + 1);
        if (igraph_vector_int_get(&degrees, from) == 1 &&
            igraph_vector_int_get(&degrees, to) == 1) {
          igraph_vector_int_push_back(&verticesToDelete, from);
          igraph_vector_int_push_back(&verticesToDelete, to);
        }
      }

      igraph_delete_vertices(&this->graph,
                             igraph_vss_vector(&verticesToDelete));

      size_t numVerticesDeleted = igraph_vector_int_size(&verticesToDelete);
      assert(numVerticesDeleted % 2 == 0);

      igraph_vector_int_destroy(&verticesToDelete);
      igraph_vector_int_destroy(&allEdges);
      igraph_vector_int_destroy(&degrees);

      if (numVerticesDeleted > 0) {
        this->renumberParentSprings();
        this->net.isUpToDate = false;
      }

      return numVerticesDeleted / 2;
    }

    /**
     * @brief When springs have been removed, it is possible that the
     * numbering is not sequential anymore. This function fixes that.
     *
     */
    void MEHPForceBalance2::renumberParentSprings()
    {
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

      igraph_vector_t parentEdges;
      igraph_vector_init(&parentEdges, igraph_ecount(&this->graph));
      igraph_cattribute_EANV(&this->graph,
                             "parent_edge",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &parentEdges);

      // find unique parent edge ids
      std::unordered_set<igraph_integer_t> parentEdgeIds;
      for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        igraph_integer_t parentEdgeId =
          castToIgraphInt(igraph_vector_get(&parentEdges, i));
        parentEdgeIds.insert(parentEdgeId);
      }

      // convert to vector, sort
      std::vector<igraph_integer_t> parentEdgeIdsVec;
      std::copy(parentEdgeIds.begin(),
                parentEdgeIds.end(),
                std::back_inserter(parentEdgeIdsVec));
      std::sort(parentEdgeIdsVec.begin(), parentEdgeIdsVec.end());

      // iterate to find gaps
      for (long int i = parentEdgeIdsVec.size() - 2; i >= 0; --i) {
        // gap found if difference > 0
        int diffForEdgeId = parentEdgeIdsVec[i + 1] - (parentEdgeIdsVec[i] + 1);
        assert(diffForEdgeId >= 0);
        if (diffForEdgeId > 0) {
          for (size_t j = 0; j < igraph_vector_size(&parentEdges); ++j) {
            igraph_integer_t parentEdgeId =
              castToIgraphInt(igraph_vector_get(&parentEdges, j));
            if (parentEdgeId >= parentEdgeIdsVec[i + 1]) {
              igraph_vector_set(
                &parentEdges,
                j,
                static_cast<igraph_real_t>(parentEdgeId - diffForEdgeId));
            }
          }
        }
      }
      // finally, close the last gap: if the smallest parent id is not 0
      for (size_t j = 0; j < igraph_vector_size(&parentEdges); ++j) {
        int lastGapDiff = parentEdgeIdsVec[0];
        igraph_integer_t parentEdgeId =
          castToIgraphInt(igraph_vector_get(&parentEdges, j));
        igraph_vector_set(
          &parentEdges,
          j,
          static_cast<igraph_real_t>(parentEdgeId - lastGapDiff));
      }

      igraph_cattribute_EAN_setv(&this->graph, "parent_edge", &parentEdges);

      // validate the new numbering
      std::unordered_set<igraph_integer_t> parentEdgeIdsNow;
      igraph_integer_t maxParentEdgeId = 0;
      for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        igraph_integer_t parentEdgeId =
          castToIgraphInt(igraph_vector_get(&parentEdges, i));
        parentEdgeIdsNow.insert(parentEdgeId);
        maxParentEdgeId = std::max(parentEdgeId, maxParentEdgeId);
      }
      assert(maxParentEdgeId == parentEdgeIdsNow.size() - 1 ||
             parentEdgeIdsNow.size() == 0);

      igraph_vector_destroy(&parentEdges);
    };

    /**
     * @brief When spring has been removed, it is possible that the
     * sum of the partitions don't add up to one anymore. This function fixes
     * that.
     *
     */
    void MEHPForceBalance2::renormalizePartitions()
    {
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
      igraph_vector_t parentEdges;
      igraph_vector_init(&parentEdges, igraph_ecount(&this->graph));
      igraph_cattribute_EANV(&this->graph,
                             "parent_edge",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &parentEdges);

      igraph_vector_t partitions;
      igraph_vector_init(&partitions, igraph_ecount(&this->graph));
      igraph_cattribute_EANV(&this->graph,
                             "partition_fraction",
                             igraph_ess_all(IGRAPH_EDGEORDER_ID),
                             &partitions);

      std::unordered_map<size_t, double> sumPerParent;
      for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        size_t parentEdge = castToIgraphInt(igraph_vector_get(&parentEdges, i));
        sumPerParent[parentEdge] = 0.;
      }
      for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        size_t parentEdge = castToIgraphInt(igraph_vector_get(&parentEdges, i));
        sumPerParent[parentEdge] += igraph_vector_get(&partitions, i);
      }
      for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
        size_t parentEdge = castToIgraphInt(igraph_vector_get(&parentEdges, i));
        igraph_vector_set(&partitions,
                          i,
                          igraph_vector_get(&partitions, i) /
                            (sumPerParent[parentEdge]));
      }

      igraph_cattribute_EAN_setv(
        &this->graph, "partition_fraction", &partitions);
      igraph_vector_destroy(&partitions);
      igraph_vector_destroy(&parentEdges);
    };

    /**
     * @brief Remove cross-links which do not have any springs with a certain
     * minimum length
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param tolerance
     */
    size_t MEHPForceBalance2::removeInactiveCrosslinks(double tolerance)
    {
      size_t numBefore = this->getNrOfSlipLinks();

      size_t numRemovedInIteration = 0;
      do {
        this->removeInactiveParentEdges(tolerance);
        numRemovedInIteration = this->removeSubfunctionalVertices();
      } while (numRemovedInIteration > 0);

      size_t numRemoved = numBefore - this->getNrOfSlipLinks();

      if (numRemoved > 0) {
        this->renumberParentSprings();
      }

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
    size_t MEHPForceBalance2::doRemovalAndreisWay(double tolerance)
    {
      size_t numRemovedTotal = this->removeSubfunctionalVertices();
      // and remove all springs that are inactive
      size_t numSpringsRemoved = this->removeInactiveParentEdges(tolerance);

      if (numSpringsRemoved > 0) {
        numRemovedTotal += this->doRemovalAndreisWay(tolerance);
      }

      if (numRemovedTotal > 0) {
        this->renumberParentSprings();
      }

      return numRemovedTotal;
    };

    /**
     * @brief Add slip-links to the structure
     *
     * @param strandIdx1 for each slip-link, the index of the first strand
     * @param strandIdx2 for each slip-link, the index of the second strand
     * @param x the x coordinate of the slip-link
     * @param y the y coordinate of the slip-link
     * @param z the z coordinate of the slip-link
     * @param alpha1 the fraction at which to insert the slip-link into strand 1
     * @param alpha2 the fraction at which to insert the slip-link into strand 2
     * @param loops
     * @param loopsOfSliplinks
     * @param clampAlpha
     */
    void MEHPForceBalance2::addSlipLinks(
      const std::vector<size_t>& strandIdx1,
      const std::vector<size_t>& strandIdx2,
      const std::vector<double>& x,
      const std::vector<double>& y,
      const std::vector<double>& z,
      const std::vector<double>& alpha1,
      const std::vector<double>& alpha2,
      std::vector<std::vector<size_t>> loops,
      std::vector<std::vector<size_t>> loopsOfSliplinks,
      bool clampAlpha)
    {
      size_t additionalLen = strandIdx1.size();
      if (additionalLen == 0) {
        return;
      }
      INVALIDARG_EXP_IFN((loopsOfSliplinks.size() == 0 && loops.size() == 0),
                         "Loops are not yet supported.");
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
      // validate inputs
      size_t currentNrOfPartialSprings = igraph_ecount(&this->graph);
      if (additionalLen != x.size() || additionalLen != y.size() ||
          additionalLen != z.size()) {
        throw std::invalid_argument("x, y and z must have the same dimensions");
      }
      if (additionalLen != strandIdx2.size() ||
          additionalLen != alpha1.size() || additionalLen != alpha2.size()) {
        throw std::invalid_argument(
          "Strand indices and alpha estimates must have the same length");
      }
      for (size_t i = 0; i < additionalLen; ++i) {
        INVALIDARG_EXP_IFN(
          strandIdx1[i] < this->net.nrOfSprings,
          "Invalid spring index " + std::to_string(strandIdx1[i]) +
            ", expected below " + std::to_string(this->net.nrOfSprings) + ".");
        INVALIDARG_EXP_IFN(
          strandIdx2[i] < this->net.nrOfSprings,
          "Invalid spring index " + std::to_string(strandIdx2[i]) +
            ", expected below " + std::to_string(this->net.nrOfSprings) + ".");
        INVALIDARG_EXP_IFN(APPROX_WITHIN(alpha1[i], 0.0, 1.0, 1e-12),
                           "Expected alpha within [0, 1], got " +
                             std::to_string(alpha1[i]) + ".");
        INVALIDARG_EXP_IFN(APPROX_WITHIN(alpha2[i], 0.0, 1.0, 1e-12),
                           "Expected alpha within [0, 1], got " +
                             std::to_string(alpha2[i]) + ".");
      }

      // then, loop the slip-links to add
      igraph_integer_t numEdgesBefore = igraph_ecount(&this->graph);
      igraph_integer_t numVerticesBefore = igraph_vcount(&this->graph);
      igraph_add_vertices(&this->graph, additionalLen, nullptr);
      for (size_t i = 0; i < additionalLen; ++i) {
        igraph_cattribute_VAN_set(
          &this->graph, "x", numVerticesBefore + i, x[i]);
        igraph_cattribute_VAN_set(
          &this->graph, "y", numVerticesBefore + i, y[i]);
        igraph_cattribute_VAN_set(
          &this->graph, "z", numVerticesBefore + i, z[i]);
        igraph_cattribute_VAN_set(
          &this->graph, "type", numVerticesBefore + i, this->slipLinkType);
        igraph_cattribute_VAN_set(
          &this->graph, "num_link_swaps", numVerticesBefore + i, 0);
        // check that our assumption regarding the vertex id is correct
        igraph_integer_t degree;
        igraph_degree_1(
          &this->graph, &degree, numVerticesBefore + i, IGRAPH_ALL, true);
        assert(degree == 0);
      }

      this->net.isUpToDate = false;

      for (size_t i = 0; i < additionalLen; ++i) {
        // repeat the insertion twice, for both strands
        for (size_t v = 0; v < 2; ++v) {
          size_t strandIdx = v == 0 ? strandIdx1[i] : strandIdx2[i];
          double alpha = v == 0 ? alpha1[i] : alpha2[i];
          this->validateIgraphSpring(strandIdx);
          double prevAlphaSum = 0.;
          igraph_integer_t edgeToReplace =
            this->findPartialSpringByFraction(strandIdx, alpha, prevAlphaSum);
          double remainingAlpha = alpha - prevAlphaSum;
          assert(igraph_cattribute_EAN(
                   &this->graph, "parent_edge", edgeToReplace) == strandIdx);
          igraph_integer_t from, to;
          igraph_edge(&this->graph, edgeToReplace, &from, &to);
          // add the four new edges
          double currentPartition = igraph_cattribute_EAN(
            &this->graph, "partition_fraction", edgeToReplace);
          assert(prevAlphaSum + currentPartition > alpha);
          double currentContourLength = igraph_cattribute_EAN(
            &this->graph, "contour_length", edgeToReplace);
          assert(APPROX_WITHIN(currentPartition, 0.0, 1.0, 1e-5));
          // for the box offsets, we can simply give them to one of the two new
          // edges – all the rest should be handled by the optimisation,
          // actually.
          assert(APPROX_WITHIN(remainingAlpha, 0.0, 1.0, 1e-5));
          this->createEdge(from,
                           numVerticesBefore + i,
                           strandIdx,
                           remainingAlpha,
                           currentContourLength,
                           this->getBondBoxOffsetForEdge(edgeToReplace));
          assert(
            APPROX_WITHIN(currentPartition - remainingAlpha, 0.0, 1.0, 1e-5));
          igraph_integer_t newEdgeId =
            this->createEdge(to,
                             numVerticesBefore + i,
                             strandIdx,
                             currentPartition - remainingAlpha,
                             currentContourLength,
                             Eigen::Vector3d::Zero());
          // check assumptions
          assert(igraph_ecount(&this->graph) == newEdgeId + 1);
          // mark the old ones to be removed
          igraph_delete_edges(&this->graph, igraph_ess_1(edgeToReplace));
          assert(igraph_ecount(&this->graph) == newEdgeId);
          this->validateIgraphSpring(strandIdx);
        }
      }

      // some more sanity checks
      size_t numNewEdges = additionalLen * 2;

      size_t numEdgesNow = igraph_ecount(&this->graph);
      assert(numEdgesNow == numNewEdges + numEdgesBefore);

      this->convertFromGraph();
    }

    //----------------------------------------------------------------
    // MARK: Structural analysis
    //----------------------------------------------------------------

    /**
     * @brief Get the Effective Functionality Of each node
     *
     * Returns the number of active springs connected to each atom, atomId
     * used as index
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return std::unordered_map<long int, int>
     */
    std::unordered_map<long int, int>
    MEHPForceBalance2::getEffectiveFunctionalityOfAtoms(double tolerance)
    {
      RUNTIME_EXP_IFN(
        this->getNetwork().isUpToDate,
        "Network is not up to date, cannot be used for computations yet.");
      std::unordered_map<long int, int> results;
      results.reserve(this->net.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        this->getNrOfActiveSpringsConnected(tolerance);
      for (size_t i = 0; i < this->net.nrOfNodes; i++) {
        results.emplace(this->net.oldAtomIds[i], nrOfActiveSpringsConnected[i]);
      }
      return results;
    }

    /**
     * @brief Get the Nr Of Active Springs connected to each node
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return Eigen::VectorXi
     */
    Eigen::VectorXi MEHPForceBalance2::getNrOfActiveSpringsConnected(
      double tolerance)
    {
      RUNTIME_EXP_IFN(
        this->getNetwork().isUpToDate,
        "Network is not up to date, cannot be used for computations yet.");
      Eigen::VectorXi nrOfActiveSpringsConnected =
        Eigen::VectorXi::Zero(this->net.nrOfNodes);
      Eigen::ArrayXb springIsActive =
        this->findActiveSprings(this->currentSpringDistances, tolerance);
      for (size_t i = 0; i < this->net.nrOfSprings; i++) {
        if (springIsActive[i] == true) /* active spring */
        {
          int a = this->net.springIndexA[i];
          int b = this->net.springIndexB[i];
          ++(nrOfActiveSpringsConnected[a]);
          ++(nrOfActiveSpringsConnected[b]);
        }
      }
      return nrOfActiveSpringsConnected;
    }

    /**
     * @brief Get the Nr Of Active Springs connected to each node
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return Eigen::VectorXi
     */
    Eigen::VectorXi MEHPForceBalance2::getNrOfActivePartialSpringsConnected(
      double tolerance)
    {
      RUNTIME_EXP_IFN(
        this->getNetwork().isUpToDate,
        "Network is not up to date, cannot be used for computations yet.");
      Eigen::VectorXi nrOfActivePartialSpringsConnected =
        Eigen::VectorXi::Zero(this->net.nrOfNodes);
      Eigen::ArrayXb springIsActive =
        this->findActiveSprings(this->currentPartialSpringDistances, tolerance);
      RUNTIME_EXP_IFN(
        springIsActive.size() == this->net.nrOfPartialSprings,
        "Expect findActiveSprings to return an "
        "appropriately sized result. Got only " +
          std::to_string(springIsActive.size()) + " entries for " +
          std::to_string(this->net.nrOfPartialSprings) + " partial springs.");
      for (size_t i = 0; i < this->net.nrOfPartialSprings; ++i) {
        if (springIsActive[i] == true) /* active spring */
        {
          int a = this->net.springPartIndexA[i];
          int b = this->net.springPartIndexB[i];
          if (!this->net.linkIsSliplink[a]) {
            ++(nrOfActivePartialSpringsConnected[a]);
          }

          if (!this->net.linkIsSliplink[b]) {
            ++(nrOfActivePartialSpringsConnected[b]);
          }
        }
      }
      return nrOfActivePartialSpringsConnected;
    }

    //----------------------------------------------------------------
    // MARK: Computations
    //----------------------------------------------------------------

    /**
     * @brief Compute the force acting on a slip-link
     *
     * TODO: use "global" partial distances
     *
     * @param linkIdx
     * @param net
     * @param u
     * @param springPartitions
     * @param kappa0
     * @param minCutoff
     * @return Eigen::Vector3d
     */
    Eigen::Matrix3d MEHPForceBalance2::evaluateForceOnLink(
      const size_t linkIdx,
      Eigen::VectorXi& debugNrSpringsVisited,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

      Eigen::Matrix3d force = Eigen::Matrix3d::Zero();

      igraph_vector_int_t edgesOfSlipLink;
      igraph_vector_int_init(&edgesOfSlipLink, 4);
      igraph_incident(&this->graph, &edgesOfSlipLink, linkIdx, IGRAPH_ALL);

      for (size_t i = 0; i < igraph_vector_int_size(&edgesOfSlipLink); ++i) {
        igraph_integer_t edgeId = igraph_vector_int_get(&edgesOfSlipLink, i);
        Eigen::Vector3d vec =
          this->evaluatePartialSpringDistanceFrom(edgeId, linkIdx, this->is2D);
        debugNrSpringsVisited[edgeId] += 1;
        const double N =
          igraph_cattribute_EAN(&this->graph, "contour_length", edgeId);
        double denominator =
          1. / (N * igraph_cattribute_EAN(
                      &this->graph, "partition_fraction", edgeId));
        if (oneOverSpringPartitionUpperLimit > 0 ||
            !std::isfinite(denominator)) {
          denominator =
            CLAMP_ONE_OVER_SPRINGPARTITION(net.partialSpringIsPartial[edgeId],
                                           denominator,
                                           N,
                                           oneOverSpringPartitionUpperLimit);
        }

        force += kappa0 * denominator * (vec * vec.transpose());
      }

      igraph_vector_int_destroy(&edgesOfSlipLink);

      return force;
    };

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3>
    MEHPForceBalance2::evaluateStressTensorLinkBased(
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit,
      const bool xlinksOnly) const
    {
      Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

      double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

      Eigen::VectorXi debugNrSpringsVisited =
        Eigen::VectorXi::Zero(net.nrOfPartialSprings);

      size_t nrOfLinksToInspect = xlinksOnly ? net.nrOfNodes : net.nrOfLinks;
      for (size_t linkIdx = 0; linkIdx < nrOfLinksToInspect; ++linkIdx) {
        Eigen::Matrix3d force =
          this->evaluateForceOnLink(linkIdx,
                                    debugNrSpringsVisited,
                                    kappa0,
                                    oneOverSpringPartitionUpperLimit);
        /* spring contribution to the overall stress tensor */
        RUNTIME_EXP_IFN(std::isfinite(force.squaredNorm()),
                        "Got non-finite force contribution to stress tensor: " +
                          std::to_string(force.squaredNorm()) + " at link " +
                          std::to_string(linkIdx) + "!");
        stress += force;
      }

      std::array<std::array<double, 3>, 3> stressA;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          stressA[i][j] = halfOverVolume * stress(i, j);
        }
      }

      if (!xlinksOnly) {
        RUNTIME_EXP_IFN(
          debugNrSpringsVisited.sum() == 2 * net.nrOfPartialSprings,
          "Every spring must be visited twice, got min " +
            std::to_string(debugNrSpringsVisited.minCoeff()) + " and max " +
            std::to_string(debugNrSpringsVisited.maxCoeff()) + ". Sum is " +
            std::to_string(debugNrSpringsVisited.sum()) + " instead of " +
            std::to_string(2 * net.nrOfPartialSprings) + ".");
        RUNTIME_EXP_IFN(
          (debugNrSpringsVisited.array() == 2).all(),
          "Every spring must be visited twice, got min " +
            std::to_string(debugNrSpringsVisited.minCoeff()) + " and max " +
            std::to_string(debugNrSpringsVisited.maxCoeff()) + ".");
      }

      return stressA;
    }

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3>
    MEHPForceBalance2::evaluateStressTensor(
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {
      std::array<std::array<double, 3>, 3> stress;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          stress[i][j] = 0.0;
        }
      }

      double oneOverVolume = 1. / (net.L[0] * net.L[1] * net.L[2]);

      for (igraph_integer_t i = 0; i < igraph_ecount(&this->graph); ++i) {
        Eigen::Vector3d distance = this->computeEdgeDistance(i);
        double denominator =
          this->getEdgeDenominator(i, oneOverSpringPartitionUpperLimit);

        /* spring contribution to the overall stress tensor */
        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            double contribution =
              distance[j] * distance[k] * kappa0 * denominator;
            RUNTIME_EXP_IFN(
              std::isfinite(contribution),
              "Got non-finite contribution to stress tensor: " +
                std::to_string(contribution) + " at coordinates " +
                std::to_string(k) + ", " + std::to_string(j) +
                " for partial spring " + std::to_string(i) +
                " from distances " + std::to_string(distance[j]) + ", " +
                std::to_string(distance[k]) + " and denominator " +
                std::to_string(denominator) + ".");
            // if (std::isfinite(denominator) && std::isfinite(contribution))
            // {
            stress[j][k] += contribution;
            // }
          }
        }
      }

      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          stress[i][j] *= oneOverVolume;
          RUNTIME_EXP_IFN(std::isfinite(stress[i][j]),
                          "Got non-finite stress tensor component: " +
                            std::to_string(stress[i][j]) + " at coordinates " +
                            std::to_string(i) + ", " + std::to_string(j) +
                            " from denominator " +
                            std::to_string(oneOverVolume) + ".");
        }
      }

      return stress;
    }

    /**
     * @brief Sum the partial spring distances for the total distance of each
     * spring
     *
     * @param net
     * @param is2D
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance2::evaluateSpringDistances()
    {
      RUNTIME_EXP_IFN(
        this->getNetwork().isUpToDate,
        "Network is not up to date, cannot be used for computations yet.");

      Eigen::VectorXd partialSpringDistances =
        this->evaluatePartialSpringDistances();
      Eigen::VectorXd result = Eigen::VectorXd::Zero(net.nrOfSprings * 3);
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        result.segment(3 * net.partialToFullSpringIndex(i), 3) +=
          partialSpringDistances.segment(3 * i, 3);
      }

      return result;
    }

    /**
     * @brief Compute the distance of each partial spring
     *
     * @param net
     * @param is2D
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance2::evaluatePartialSpringDistances() const
    {
      Eigen::VectorXd springDistances =
        Eigen::VectorXd::Zero(3 * igraph_ecount(&this->graph));

      for (igraph_integer_t i = 0; i < igraph_ecount(&this->graph); ++i) {
        springDistances.segment(3 * i, 3) = this->computeEdgeDistance(i);
      }

      return springDistances;
    }

    /**
     * @brief Get the average spring length at the current step
     *
     * @return double
     */
    double MEHPForceBalance2::getAverageSpringLength()
    {
      RUNTIME_EXP_IFN(
        this->getNetwork().isUpToDate,
        "Network is not up to date, cannot be used for computations yet.");

      double r2 = 0.0;
      // TODO: this is sketchy
      for (int i = 0; i < this->net.nrOfSprings; i++) {
        double r2local = 0.0;
        for (int j = 0; j < 3; ++j) {
          r2local += this->currentSpringDistances[i * 3 + j] *
                     this->currentSpringDistances[i * 3 + j];
        }
        r2 += sqrt(r2local);
      }
      return r2 / this->net.nrOfSprings;
    }

    /**
     * @brief Compute the displacement residual
     *
     * @return double
     */
    double MEHPForceBalance2::getDisplacementResidualNorm(
      double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::VectorXd overallForces =
        Eigen::VectorXd::Zero(3 * igraph_vcount(&this->graph));
      for (igraph_integer_t i = 0; i < igraph_ecount(&this->graph); ++i) {
        Eigen::Vector3d distance = this->computeEdgeDistance(i);
        double denominator =
          this->getEdgeDenominator(i, oneOverSpringPartitionUpperLimit);
        igraph_integer_t from, to;
        igraph_edge(&this->graph, i, &from, &to);
        overallForces.segment(3 * from, 3) += denominator * distance;
        overallForces.segment(3 * to, 3) -= denominator * distance;
      }

      return overallForces.squaredNorm();
    }

    //----------------------------------------------------------------
    // MARK: Data access
    //----------------------------------------------------------------

    /**
     * @brief Get the Ids Of active Nodes
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @param minimumNrOfActiveConnections the number of active springs
     * required for this node to qualify as active
     * @return std::vector<long int> the atom ids
     */
    std::vector<long int> MEHPForceBalance2::getIdsOfActiveNodes(
      double tolerance,
      int minimumNrOfActiveConnections,
      int maximumNrOfActiveConnections,
      bool usePartial)
    {
      RUNTIME_EXP_IFN(
        this->getNetwork().isUpToDate,
        "Network is not up to date, cannot be used for computations yet.");

      std::vector<long int> results;
      results.reserve(this->net.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        usePartial ? this->getNrOfActiveSpringsConnected(tolerance)
                   : this->getNrOfActivePartialSpringsConnected(tolerance);
      for (size_t i = 0; i < this->net.nrOfNodes; i++) {
        if (nrOfActiveSpringsConnected[i] >= minimumNrOfActiveConnections &&
            (maximumNrOfActiveConnections < 0 ||
             maximumNrOfActiveConnections >= nrOfActiveSpringsConnected[i])) {
          results.push_back(this->net.oldAtomIds[i]);
        }
      }

      return results;
    }

    /**
     * @brief Convert the net to a Universe back
     *
     * @return pylimer_tools::entities::Universe
     */
    pylimer_tools::entities::Universe MEHPForceBalance2::getCrosslinkerVerse()
    {
      RUNTIME_EXP_IFN(
        this->getNetwork().isUpToDate,
        "Network is not up to date, cannot be used for computations yet.");

      // convert nodes & springs back to a universe
      pylimer_tools::entities::Universe xlinkUniverse =
        pylimer_tools::entities::Universe(this->universe.getBox());
      std::vector<long int> ids;
      std::vector<int> types = pylimer_tools::utils::initializeWithValue(
        this->net.nrOfNodes, crosslinkerType);
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      std::vector<int> zeros =
        pylimer_tools::utils::initializeWithValue(this->net.nrOfNodes, 0);
      ids.reserve(this->net.nrOfNodes);
      x.reserve(this->net.nrOfNodes);
      y.reserve(this->net.nrOfNodes);
      z.reserve(this->net.nrOfNodes);
      for (int i = 0; i < this->net.nrOfNodes; ++i) {
        x.push_back(this->net.coordinates[3 * i + 0]);
        y.push_back(this->net.coordinates[3 * i + 1]);
        z.push_back(this->net.coordinates[3 * i + 2]);
        ids.push_back(this->net.oldAtomIds[i]);
        // override type, since the types may be different from crosslinkerType
        // if converted with dangling chains
        types[i] = this->net.linkIsSliplink[i] ? this->slipLinkType
                                               : this->crosslinkerType;
      }
      xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
      std::vector<long int> bondFrom;
      std::vector<long int> bondTo;
      bondFrom.reserve(this->net.nrOfSprings);
      bondTo.reserve(this->net.nrOfSprings);
      for (int i = 0; i < this->net.nrOfSprings; ++i) {
        bondFrom.push_back(this->net.oldAtomIds[this->net.springIndexA[i]]);
        bondTo.push_back(this->net.oldAtomIds[this->net.springIndexB[i]]);
      }
      xlinkUniverse.addBonds(
        bondFrom.size(),
        bondFrom,
        bondTo,
        pylimer_tools::utils::initializeWithValue(bondFrom.size(), 1),
        false,
        false); // disable simplify to keep the self-loops etc.
      return xlinkUniverse;
    }

    //----------------------------------------------------------------
    // MARK: Network Validation
    //----------------------------------------------------------------

    bool MEHPForceBalance2::validateNetwork(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& springPartitions)
    {
      RUNTIME_EXP_IFN(net.isUpToDate,
                      "Network is not up to date. Validation is fruitless.");
      // std::cout << "Validating network..." << std::endl;
      /**
       * First, test dimensions
       */
      RUNTIME_EXP_IFN(!std::isinf(net.L[0]) && !std::isnan(net.L[0]),
                      "Box direction x must be scalar");
      RUNTIME_EXP_IFN(!std::isinf(net.L[1]) && !std::isnan(net.L[1]),
                      "Box direction y must be scalar");
      RUNTIME_EXP_IFN(!std::isinf(net.L[2]) && !std::isnan(net.L[2]),
                      "Box direction z must be scalar");
      RUNTIME_EXP_IFN(net.coordinates.size() == net.nrOfLinks * 3,
                      "Invalid size of coordinates");
      RUNTIME_EXP_IFN(net.localToGlobalSpringIndex.size() == net.nrOfSprings,
                      "Invalid size of connectivity map, got " +
                        std::to_string(net.localToGlobalSpringIndex.size()) +
                        " for " + std::to_string(net.nrOfSprings) +
                        " springs.");
      RUNTIME_EXP_IFN(net.springsContourLength.size() == net.nrOfSprings,
                      "Invalid size of contour lengths, got " +
                        std::to_string(net.springsContourLength.size()) +
                        " for " + std::to_string(net.nrOfSprings) +
                        " springs.");
      RUNTIME_EXP_IFN(net.springIndicesOfLinks.size() == net.nrOfLinks,
                      "Invalid size of spring indices of links, got " +
                        std::to_string(net.linkIndicesOfSprings.size()) +
                        " for " + std::to_string(net.nrOfSprings) +
                        " springs.");
      RUNTIME_EXP_IFN(net.linkIndicesOfSprings.size() == net.nrOfSprings,
                      "Invalid size of link indices of springs, got " +
                        std::to_string(net.linkIndicesOfSprings.size()) +
                        " for " + std::to_string(net.nrOfSprings) +
                        " springs.");
      RUNTIME_EXP_IFN(net.linkIsSliplink.size() == net.nrOfLinks,
                      "Invalid size of link is sliplink");
      RUNTIME_EXP_IFN(net.nrOfCrosslinkSwapsEndured.size() ==
                        net.nrOfLinks - net.nrOfNodes,
                      "Invalid size of link is sliplink");
      RUNTIME_EXP_IFN(
        net.linkIsSliplink.count() == (net.nrOfLinks - net.nrOfNodes),
        "Nr of nodes plus nr of slp-links should give the total nr of links. "
        "Got " +
          std::to_string(net.linkIsSliplink.count()) +
          " links marked as slip-link, but " + std::to_string(net.nrOfLinks) +
          " total links, of which " + std::to_string(net.nrOfNodes) +
          " are cross-links.");
      RUNTIME_EXP_IFN(net.oldAtomIds.size() == net.nrOfNodes,
                      "Invalid size of old atom ids");
      RUNTIME_EXP_IFN(net.springCoordinateIndexA.size() == net.nrOfSprings * 3,
                      "Invalid size of springCoordinateIndexA");
      RUNTIME_EXP_IFN(net.springCoordinateIndexB.size() == net.nrOfSprings * 3,
                      "Invalid size of springCoordinateIndexB");
      RUNTIME_EXP_IFN(net.springPartCoordinateIndexA.size() ==
                        net.nrOfPartialSprings * 3,
                      "Invalid size of springPartCoordinateIndexA");
      RUNTIME_EXP_IFN(net.springPartCoordinateIndexB.size() ==
                        net.nrOfPartialSprings * 3,
                      "Invalid size of springPartCoordinateIndexB");
      RUNTIME_EXP_IFN(net.springIndexA.size() == net.nrOfSprings,
                      "Invalid size of springIndexA");
      RUNTIME_EXP_IFN(net.springIndexB.size() == net.nrOfSprings,
                      "Invalid size of springIndexB");
      RUNTIME_EXP_IFN(net.springPartIndexA.size() == net.nrOfPartialSprings,
                      "Invalid size of springPartIndexA");
      RUNTIME_EXP_IFN(net.springPartBoxOffset.size() ==
                        net.nrOfPartialSprings * 3,
                      "Invalid size of springPartBoxOffset");
      RUNTIME_EXP_IFN(net.springPartIndexB.size() == net.nrOfPartialSprings,
                      "Invalid size of springPartIndexB");
      RUNTIME_EXP_IFN(net.springIsActive.size() == net.nrOfSprings,
                      "Invalid size of springIsActive");
      RUNTIME_EXP_IFN(springPartitions.size() == net.nrOfPartialSprings,
                      "Invalid size of spring partitions, got " +
                        std::to_string(springPartitions.size()) + " for " +
                        std::to_string(net.nrOfPartialSprings) +
                        " partial springs.");
      RUNTIME_EXP_IFN(net.partialSpringIsPartial.size() ==
                        net.nrOfPartialSprings,
                      "Invalid size of partialSpringIsPartial");
      RUNTIME_EXP_IFN(
        APPROX_EQUAL(springPartitions.sum(), net.nrOfSprings, 1e-3),
        "Spring partitions should sum to 1 per spring, got " +
          std::to_string(springPartitions.sum()) + " for " +
          std::to_string(net.nrOfSprings) + " springs.");
      RUNTIME_EXP_IFN(
        net.partialToFullSpringIndex.size() == net.nrOfPartialSprings,
        "Every partial spring must be able to map to the full spring.");

      /**
       * Test maximum values
       */
      if (net.nrOfSprings > 0) {
        RUNTIME_EXP_IFN(
          net.partialToFullSpringIndex.maxCoeff() < net.nrOfSprings,
          "Partial spring must map to full spring, which must have "
          "a lower index.");
        RUNTIME_EXP_IFN(net.springPartCoordinateIndexA.maxCoeff() <
                          3 * net.nrOfLinks,
                        "Part coordinates must map to coordinates.");
        RUNTIME_EXP_IFN(net.springPartCoordinateIndexB.maxCoeff() <
                          3 * net.nrOfLinks,
                        "Part coordinates must map to coordinates.");
        RUNTIME_EXP_IFN(net.springPartIndexA.maxCoeff() < net.nrOfLinks,
                        "Part indices must map to links.");
        RUNTIME_EXP_IFN(net.springPartIndexB.maxCoeff() < net.nrOfLinks,
                        "Part indices must map to links.");
        RUNTIME_EXP_IFN(net.springIndexA.maxCoeff() < net.nrOfNodes,
                        "Full springs must consist of cross-links only.");
        RUNTIME_EXP_IFN(net.springIndexB.maxCoeff() < net.nrOfNodes,
                        "Full springs must consist of cross-links only.");
      }

      /**
       * Test spring partition assumptions
       */
      for (size_t i = 0; i < springPartitions.size(); i++) {
        RUNTIME_EXP_IFN(APPROX_WITHIN(springPartitions[i], 0.0, 1.0, 1e-9),
                        "Spring partitions must be between 0. & 1., got " +
                          std::to_string(springPartitions[i]) +
                          " at i = " + std::to_string(i) + ".");
      }

      /**
       * Test reversibility of link <-> spring mapping
       */
      for (size_t link_idx = 0; link_idx < net.nrOfLinks; ++link_idx) {
        RUNTIME_EXP_IFN(
          net.linkIsSliplink[link_idx] == (link_idx >= net.nrOfNodes),
          "Expected slip-links to come sequentially after cross-links.");
        std::vector<size_t> thisLinksSprings =
          net.springIndicesOfLinks[link_idx];
        for (size_t spring_idx : thisLinksSprings) {
          std::vector<size_t> thisSpringsLinks =
            net.linkIndicesOfSprings[spring_idx];
          RUNTIME_EXP_IFN(std::find(thisSpringsLinks.begin(),
                                    thisSpringsLinks.end(),
                                    link_idx) != thisSpringsLinks.end(),
                          "Spring must have a connection to the link, too. Did "
                          "not find link " +
                            std::to_string(link_idx) + " in spring " +
                            std::to_string(spring_idx) + ".");
        }
      }

      /**
       * Test the assumptions on slip-links
       */
      for (size_t slipLinkIdx = net.nrOfNodes; slipLinkIdx < net.nrOfLinks;
           ++slipLinkIdx) {
        RUNTIME_EXP_IFN(
          net.springIndicesOfLinks[slipLinkIdx].size() == 2 ||
            net.springIndicesOfLinks[slipLinkIdx].size() == 1,
          "Expect each slip-link to be involved in exactly one or two "
          "springs, "
          "got " +
            std::to_string(net.springIndicesOfLinks[slipLinkIdx].size()) + ".");
        RUNTIME_EXP_IFN(net.linkIsSliplink[slipLinkIdx],
                        "Expected slip-links to know what they are.");
      }

      /**
       * Test the validitiy of springs and their mapping
       */
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        RUNTIME_EXP_IFN(net.linkIndicesOfSprings[i].size() >= 2,
                        "Each spring requires at least two links, got " +
                          std::to_string(net.linkIndicesOfSprings[i].size()) +
                          " at i = " + std::to_string(i) + ".");
        RUNTIME_EXP_IFN(
          net.localToGlobalSpringIndex[i].size() ==
            net.linkIndicesOfSprings[i].size() - 1,
          "Require a global index for each local one, got " +
            std::to_string(net.localToGlobalSpringIndex[i].size()) +
            " != " + std::to_string(net.linkIndicesOfSprings[i].size() - 1) +
            " for spring " + std::to_string(i) + ".");
        for (size_t partialIdx = 0;
             partialIdx < net.localToGlobalSpringIndex[i].size();
             ++partialIdx) {
          size_t partialSpringIdx = net.localToGlobalSpringIndex[i][partialIdx];
          size_t partner0 = net.linkIndicesOfSprings[i][partialIdx];
          size_t partner1 = net.linkIndicesOfSprings[i][partialIdx + 1];
          RUNTIME_EXP_IFN(
            ((net.springPartIndexA[partialSpringIdx] == partner0 &&
              net.springPartIndexB[partialSpringIdx] == partner1) ||
             (net.springPartIndexB[partialSpringIdx] == partner0 &&
              net.springPartIndexA[partialSpringIdx] == partner1)),
            "Expect linkIndicesOfSprings and localToGlobalSpringIndex ordering "
            "to correspond. Got partner0 = " +
              std::to_string(partner0) + ", partner1 = " +
              std::to_string(partner1) + " vs. springs part indices " +
              std::to_string(net.springPartIndexA[partialSpringIdx]) + " and " +
              std::to_string(net.springPartIndexB[partialSpringIdx]) +
              " in spring " + std::to_string(i) + " (partial: " +
              std::to_string(partialSpringIdx) + ") with global indices " +
              pylimer_tools::utils::join(
                net.localToGlobalSpringIndex[i].begin(),
                net.localToGlobalSpringIndex[i].end(),
                std::string(", ")) +
              " and links " +
              pylimer_tools::utils::join(net.linkIndicesOfSprings[i].begin(),
                                         net.linkIndicesOfSprings[i].end(),
                                         std::string(", ")) +
              ".");
        }
        // the following is not guraranteed anymore with the removal of links
        // while running RUNTIME_EXP_IFN(
        //   net.linkIndicesOfSprings[i][0] <=
        //     net.linkIndicesOfSprings[i][net.linkIndicesOfSprings[i].size()
        //     - 1],
        //   "Springs must have increasing end-point indices");
        std::vector<size_t> links = net.linkIndicesOfSprings[i];
        for (size_t j = 0; j < links.size(); ++j) {
          size_t link_idx = links[j];
          RUNTIME_EXP_IFN(net.linkIsSliplink[link_idx] ==
                            ((j != 0) && (j != (links.size() - 1))),
                          "Cross-links must be first and last in a spring, "
                          "slip-links in-between. Found discrepancy at " +
                            std::to_string(j) + "/" +
                            std::to_string(links.size()) + " in spring " +
                            std::to_string(i) + ".")
          std::vector<size_t> thisLinksSprings =
            net.springIndicesOfLinks[link_idx];
          RUNTIME_EXP_IFN(
            std::find(thisLinksSprings.begin(), thisLinksSprings.end(), i) !=
              thisLinksSprings.end(),
            "Link must have a connection to the spring, too. Did not find "
            "spring " +
              std::to_string(i) + " in link " + std::to_string(link_idx) + ".");
        }
        // also check the sum of the partials
        std::vector<size_t> globalSpringIndices =
          net.localToGlobalSpringIndex[i];
        double sum = 0.0;
        for (size_t globalIdx : globalSpringIndices) {
          sum += springPartitions[globalIdx];
        }
        RUNTIME_EXP_IFN(
          APPROX_EQUAL(sum, 1.0, 1e-10),
          "Spring partitions of one spring must sum to one, got " +
            std::to_string(sum) + " for spring " + std::to_string(i) + ".");
      }

      /**
       * Test the validity of partial springs and their mapping
       */
      for (size_t i = 0; i < net.nrOfPartialSprings; i++) {
        size_t fullIdx = net.partialToFullSpringIndex[i];
        size_t partialEndA = net.springPartIndexA[i];
        size_t partialEndB = net.springPartIndexB[i];
        RUNTIME_EXP_IFN(partialEndA < net.nrOfLinks,
                        "Cannot have a spring (" + std::to_string(i) +
                          ") part larger " + std::to_string(partialEndA) +
                          " than the nr of links (" +
                          std::to_string(net.nrOfLinks) + ").")
        RUNTIME_EXP_IFN(partialEndB < net.nrOfLinks,
                        "Cannot have a spring (" + std::to_string(i) +
                          ") part larger " + std::to_string(partialEndB) +
                          " than the nr of links (" +
                          std::to_string(net.nrOfLinks) + ").")
        RUNTIME_EXP_IFN(
          (net.linkIsSliplink[partialEndA] ||
           net.linkIsSliplink[partialEndB]) == net.partialSpringIsPartial[i],
          "Springs involving slip-links must be marked partial. Spring " +
            std::to_string(i) + " is marked: " +
            std::to_string(net.partialSpringIsPartial[i]) + ".");
        RUNTIME_EXP_IFN(
          (net.linkIndicesOfSprings[net.partialToFullSpringIndex[i]].size() >
           2) == net.partialSpringIsPartial[i],
          "Springs involving slip-links must be marked partial. Spring " +
            std::to_string(i) + " is marked: " +
            std::to_string(net.partialSpringIsPartial[i]) + ".");
        if (!net.linkIsSliplink[partialEndA]) {
          RUNTIME_EXP_IFN(
            net.springIndexA[fullIdx] == partialEndA ||
              net.springIndexB[fullIdx] == partialEndA,
            "Expect mapping of springs to work: " +
              std::to_string(partialEndA) +
              " is a cross-link, yet not part of the two ends of spring " +
              std::to_string(fullIdx) + ", where we have " +
              std::to_string(net.springIndexA[fullIdx]) + " and " +
              std::to_string(net.springIndexB[fullIdx]) + ".");
        }
        if (!net.linkIsSliplink[partialEndB]) {
          RUNTIME_EXP_IFN(
            net.springIndexA[fullIdx] == partialEndB ||
              net.springIndexB[fullIdx] == partialEndB,
            "Expect mapping of springs to work: " +
              std::to_string(partialEndB) +
              " is a cross-link, yet not part of the two ends of spring " +
              std::to_string(fullIdx) + ", where we have " +
              std::to_string(net.springIndexA[fullIdx]) + " and " +
              std::to_string(net.springIndexB[fullIdx]) + ".");
        }
        RUNTIME_EXP_IFN(net.springPartCoordinateIndexA[3 * i] % 3 == 0,
                        "Expected spring part coordinates to be sequentially "
                        "built from spring parts.");
        RUNTIME_EXP_IFN(net.springPartCoordinateIndexB[3 * i] % 3 == 0,
                        "Expected spring part coordinates to be sequentially "
                        "built from spring parts.");
        for (size_t dir = 0; dir < 3; ++dir) {
          RUNTIME_EXP_IFN(
            net.springPartCoordinateIndexA[3 * i + dir] ==
              3 * partialEndA + dir,
            "Spring part index and coordinate index must match. Got " +
              std::to_string(net.springPartCoordinateIndexA[3 * i + dir]) +
              " but expected " + std::to_string(3 * partialEndA + dir) +
              " with dir = " + std::to_string(dir) + ".");
          RUNTIME_EXP_IFN(
            net.springPartCoordinateIndexB[3 * i + dir] ==
              3 * partialEndB + dir,
            "Spring part index and coordinate index must match.Got " +
              std::to_string(net.springPartCoordinateIndexB[3 * i + dir]) +
              " but expected " + std::to_string(3 * partialEndB + dir) +
              " with dir = " + std::to_string(dir) + ".");
        }
      }

      /**
       * Check that we do not have any nan or inf values in our vectors
       */
      for (size_t coordI = 0; coordI < net.coordinates.size(); coordI++) {
        RUNTIME_EXP_IFN(std::isfinite(net.coordinates[coordI]),
                        "Coordinate component " + std::to_string(coordI) +
                          " must be finite, got " +
                          std::to_string(net.coordinates[coordI]) + ".");
      }
      for (size_t dir = 0; dir < 3; ++dir) {
        RUNTIME_EXP_IFN(std::isfinite(net.L[dir]),
                        "Expected box size to be finite, got " +
                          std::to_string(net.L[dir]) + " in dir " +
                          std::to_string(dir) + ".");
        RUNTIME_EXP_IFN(net.L[dir] > 0.0,
                        "Expected box size to be positive, got " +
                          std::to_string(net.L[dir]) + " in dir " +
                          std::to_string(dir) + ".");
        RUNTIME_EXP_IFN(
          APPROX_EQUAL(net.boxHalfs[dir], 0.5 * net.L[dir], 1e-12),
          "Expected box half to be half of box length");
      }

      /**
       * Validate additional loop-specific data that might not apply
       */
      if (net.loopsOfSliplink.size() > 0) {
        RUNTIME_EXP_IFN(net.loops.size() > 0, "Inconsistent use of loops.");
        RUNTIME_EXP_IFN(
          net.loopsOfSliplink.size() == (net.nrOfLinks - net.nrOfNodes),
          "Each slip-link must have associated list of loops, or none.");
        for (const std::vector<size_t>& loopsOfSliplink : net.loopsOfSliplink) {
          RUNTIME_EXP_IFN(
            loopsOfSliplink.size() <= 2,
            "Cannot have a slip-link attributed to more than two loops.");
          for (size_t loopIdx : loopsOfSliplink) {
            RUNTIME_EXP_IFN(loopIdx < net.loops.size(),
                            "Loop index out of range.");
          }
        }
        for (const std::vector<size_t>& loop : net.loops) {
          for (size_t i : loop) {
            RUNTIME_EXP_IFN(i >= net.nrOfSprings,
                            "Loop's spring index out of range.");
          }
        }
      }

      // std::cout << "Validation passed." << std::endl;
      return true;
    }

  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
