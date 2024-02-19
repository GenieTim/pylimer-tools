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
        igraph_cattribute_GAB_set(&this->graph, "is_up_to_date", false);

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
        igraph_cattribute_GAB_set(&this->graph, "is_up_to_date", false);

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
            size_t nRemoved = this->removeTwofunctionalLinks();
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
      this->currentSpringDistances =
        this->evaluateSpringDistances(net, this->is2D);
      this->currentPartialSpringDistances =
        this->evaluatePartialSpringDistances(net, this->is2D);
    }

    /**
     * @brief Displace one link to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @param linkIdx the idx of the link to displace
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance2::displaceToMeanPosition(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& springPartitions,
      const size_t linkIdx,
      const double oneOverSpringPartitionUpperLimit) const
    {
      std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
      // Eigen::Vector3d currentDisplacement = u.segment(3 * linkIdx, 3);
      Eigen::Vector3d objectiveDisplacement =
        Eigen::Vector3d::Zero(); // = remainingDisplacement.array();
      double objectiveDisplacementContributors = 0.0;
      bool cautionPrimaryLoop = false;
      std::vector<size_t> handledSprings;
      for (size_t spring_index = 0; spring_index < springIndices.size();
           ++spring_index) {
        if (cautionPrimaryLoop) {
          if (std::find(handledSprings.begin(),
                        handledSprings.end(),
                        springIndices[spring_index]) != handledSprings.end()) {
            std::cout << "primary loop detection necessary" << std::endl;
            continue;
          }
        }
        // compute partial distances & total distance of this spring
        std::vector<size_t> springsPartners =
          net.linkIndicesOfSprings[springIndices[spring_index]];

        if (springsPartners[0] == linkIdx &&
            springsPartners[springsPartners.size() - 1] == linkIdx) {
          // this is a primary loop of some kind. Make sure not to go over it
          // twice.
          cautionPrimaryLoop = true;
          // skip unentangled loop contribution
          if (springsPartners.size() == 2) {
            continue;
          }
        }

        for (size_t partner_idx = 0; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          if (springsPartners[partner_idx] == linkIdx ||
              springsPartners[partner_idx + 1] == linkIdx) {
            size_t globalSpringIndex =
              net.localToGlobalSpringIndex[(springIndices[spring_index])]
                                          [partner_idx];
            Eigen::Vector3d partialDistance =
              this->evaluatePartialSpringDistanceFrom(
                net, globalSpringIndex, linkIdx, this->is2D);
            // add to displacement
            double contourLengthFraction = springPartitions[globalSpringIndex];
            // std::cout << "Contribution from " <<
            // springsPartners[partner_idx]
            //           << " to " << springsPartners[partner_idx + 1]
            //           << " with l = " << contourLengthFraction << " and N =
            //           "
            //           <<
            //           net.springsContourLength[springIndices[spring_index]]
            //           << ", partial distance " << partialDistance[0] << ",
            //           "
            //           << partialDistance[1] << ", " << partialDistance[2]
            //           << std::endl;
            const double N =
              net.springsContourLength[springIndices[spring_index]];
            double oneOverContourLengthFraction =
              1.0 / (N * contourLengthFraction);
            // if (!std::isfinite(oneOverContourLengthFraction)) {
            //   oneOverContourLengthFraction =
            //     1.0 / (1e-12 *
            //            net.springsContourLength[springIndices[spring_index]]);
            // }
            // if (oneOverSpringPartitionUpperLimit > 0.0) {
            //   oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
            //     net.partialSpringIsPartial[globalSpringIndex],
            //     oneOverContourLengthFraction,
            //     N,
            //     oneOverSpringPartitionUpperLimit);
            // }
            if (std::isfinite(oneOverContourLengthFraction)) {
              objectiveDisplacement +=
                (partialDistance)*oneOverContourLengthFraction; // /
              // totalDistance.array());
              objectiveDisplacementContributors += oneOverContourLengthFraction;
            }
            // else {
            //   objectiveDisplacement = 1e9 * (partialDistance);
            //   objectiveDisplacementContributors += 1e9;
            // }

            if (cautionPrimaryLoop) {
              handledSprings.push_back(springIndices[spring_index]);
            }
          }
        }
      }
      // take mean for displacement
      // prevent NaN from division by zero
      Eigen::Vector3d coordsBefore = net.coordinates.segment(3 * linkIdx, 3);
      net.coordinates.segment(3 * linkIdx, 3) +=
        objectiveDisplacement / (objectiveDisplacementContributors == 0.0
                                   ? 1.0
                                   : objectiveDisplacementContributors);

      double dist = (net.coordinates.segment(3 * linkIdx, 3) - coordsBefore)
                      .segment(3 * linkIdx, 3)
                      .squaredNorm();
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

    /**
     * @brief Iterate the links and displace them to their mean position
     *
     * @param net
     * @param springPartitions0
     * @param damping
     * @return double
     */
    double MEHPForceBalance2::displaceLinksToMeanPosition(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& springPartitions0,
      double damping) const
    {
      Eigen::ArrayXi mask =
        Eigen::ArrayXi::LinSpaced(3 * net.nrOfLinks, 0, 3 * net.nrOfLinks - 1);
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(net, springPartitions0);

      return this->displaceLinksToMeanPosition(
        net, oneOverSpringPartitions, mask, damping);
    }

    /**
     * @brief Displace one link to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance2::displaceLinksToMeanPosition(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& oneOverSpringPartitions,
      const Eigen::ArrayXi& resultingCoordinateIndexMask,
      const double damping) const
    {
      Eigen::ArrayXi involvedSpringPartCoordinateIndexMask =
        Eigen::ArrayXi::LinSpaced(
          3 * net.nrOfPartialSprings, 0, 3 * net.nrOfPartialSprings - 1);
      return this->displaceLinksToMeanPosition(
        net,
        oneOverSpringPartitions,
        involvedSpringPartCoordinateIndexMask,
        resultingCoordinateIndexMask,
        damping);
    }

    /**
     * @brief Displace one link to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @param u the current displacements, wherein the resulting coordinates
     * shall be stored
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance2::displaceLinksToMeanPosition(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& oneOverSpringPartitions,
      const Eigen::ArrayXi& involvedSpringPartCoordinateIndexMask,
      const Eigen::ArrayXi& resultingCoordinateIndexMask,
      const double damping) const
    {
      assert(this->net.isUpToDate);
      INVALIDARG_EXP_IFN(oneOverSpringPartitions.size() ==
                           net.springPartCoordinateIndexB.size(),
                         "Spring partitions must have the size of the nr of "
                         "spring coordinates");
      INVALIDARG_EXP_IFN(resultingCoordinateIndexMask.size() % 3 == 0,
                         "Mask is expected to mask the coordinates");
      INVALIDARG_EXP_IFN(involvedSpringPartCoordinateIndexMask.size() % 3 == 0,
                         "Mask is expected to mask the coordinates");

      Eigen::ArrayXi relevantSpringPartCoordinateIndexA =
        net.springPartCoordinateIndexA(involvedSpringPartCoordinateIndexMask);
      Eigen::ArrayXi relevantSpringPartCoordinateIndexB =
        net.springPartCoordinateIndexB(involvedSpringPartCoordinateIndexMask);

      // assert(relevantSpringPartCoordinateIndexB.size() ==
      //        relevantSpringPartCoordinateIndexA.size());
      // assert(relevantSpringPartCoordinateIndexA.size() ==
      //        involvedSpringPartCoordinateIndexMask.size());
      // assert(resultingCoordinateIndexMask.maxCoeff() <
      //        net.coordinates.size());
      // assert(involvedSpringPartCoordinateIndexMask.maxCoeff() <
      //        net.springCoordinateIndexA.size());

      // TODO: we could save some time and space by directly adjusting the
      // coordinates
      Eigen::VectorXd relevantPartialDistancesA =
        (net.coordinates(relevantSpringPartCoordinateIndexB) -
         net.coordinates(relevantSpringPartCoordinateIndexA)) +
        net.springPartBoxOffset;
      // TODO: implement box not large enough case
      if (this->assumeBoxLargeEnough) {
        this->universe.getBox().handlePBC(relevantPartialDistancesA);
      }
      // NOTE: we have many zeros too much, here, actually.
      Eigen::ArrayXd oneOverSumOfSpringPartials = Eigen::ArrayXd::Zero(
        3 * net.nrOfLinks); // 0.0 for equal = primary loop, 1.0 otherwise
      oneOverSumOfSpringPartials(relevantSpringPartCoordinateIndexA) +=
        oneOverSpringPartitions(involvedSpringPartCoordinateIndexMask).array();
      oneOverSumOfSpringPartials(relevantSpringPartCoordinateIndexB) +=
        oneOverSpringPartitions(involvedSpringPartCoordinateIndexMask).array();
      // prevent NaN values from dividing by zero afterwards
      oneOverSumOfSpringPartials = (oneOverSumOfSpringPartials <= 1e-12)
                                     .select(1.0, oneOverSumOfSpringPartials);

      Eigen::VectorXd partialDistancesOverSpringPartitions =
        (relevantPartialDistancesA.array() *
         oneOverSpringPartitions(involvedSpringPartCoordinateIndexMask).array())
          .matrix();
      // NOTE: we have many zeros too much, here, actually.
      Eigen::VectorXd objectiveDisplacements =
        Eigen::VectorXd::Zero(3 * net.nrOfLinks);
      objectiveDisplacements(relevantSpringPartCoordinateIndexA) +=
        partialDistancesOverSpringPartitions;
      objectiveDisplacements(relevantSpringPartCoordinateIndexB) -=
        partialDistancesOverSpringPartitions;
      // ...and take the average
      objectiveDisplacements(resultingCoordinateIndexMask) =
        (objectiveDisplacements(resultingCoordinateIndexMask).array() /
         oneOverSumOfSpringPartials(resultingCoordinateIndexMask))
          .matrix();

      // reset for 2D systems
      if (this->is2D) {
        objectiveDisplacements(Eigen::seq(2, net.nrOfLinks, 3)) =
          Eigen::VectorXd::Zero(net.nrOfLinks);
      }

      // find the actual (max) displacement we did
      double maxDiff = 0.;
      Eigen::VectorXd objectivesToSet =
        objectiveDisplacements(resultingCoordinateIndexMask);
      for (size_t i = 0; i < resultingCoordinateIndexMask.size() / 3; i++) {
        maxDiff =
          std::max(maxDiff, objectivesToSet.segment(3 * i, 3).squaredNorm());
      }
      // double maxDiff =
      // (objectiveDisplacements(mask)).cwiseAbs2().maxCoeff();
      net.coordinates(resultingCoordinateIndexMask) +=
        damping * objectivesToSet;

      // igraph_cattribute_GAB_set(&this->graph, "is_up_to_date", false);

      return maxDiff;
    };

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
      assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
      // figure out, where the rail is currently going through.
      igraph_t subgraph;
      igraph_empty(&subgraph, 0, IGRAPH_UNDIRECTED);
      igraph_subgraph_from_edges(
        &this->graph, &subgraph, igraph_ess_vector(unorderedEdges), true);

      igraph_vector_int_t edgesOnPath;
      igraph_vector_int_init(&edgesOnPath, igraph_ecount(&subgraph));
      igraph_vector_int_t verticesOnPath;
      igraph_vector_int_init(&verticesOnPath, igraph_vcount(&subgraph));
      igraph_eulerian_path(&subgraph, &edgesOnPath, &verticesOnPath);
      assert(igraph_vector_int_size(&edgesOnPath) == igraph_ecount(&subgraph));

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
      }
      if (vertices != nullptr) {
        igraph_vector_int_clear(vertices);
        for (size_t i = 0; i < igraph_vector_int_size(&verticesOnPath); ++i) {
          igraph_vector_int_push_back(
            vertices,
            castToIgraphInt(igraph_cattribute_EAN(
              &subgraph,
              "prev_vertex_id",
              igraph_vector_int_get(&verticesOnPath, i))));
        }
      }
      igraph_integer_t vertex0Id = castToIgraphInt(
        igraph_cattribute_EAN(&subgraph,
                              "prev_vertex_id",
                              igraph_vector_int_get(&verticesOnPath, 0)));
      igraph_integer_t vertexEndId = castToIgraphInt(igraph_cattribute_EAN(
        &subgraph,
        "prev_vertex_id",
        igraph_vector_int_get(&verticesOnPath,
                              igraph_vector_int_size(&verticesOnPath) - 1)));
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

      igraph_destroy(&subgraph);
      igraph_vector_int_destroy(&edgesOnPath);
      igraph_vector_int_destroy(&verticesOnPath);
    }

    /**
     * @brief Given a vertex id and a rail edge, returns the other two edges
     * that are not part of the rail
     */
    std::vector<igraph_integer_t> MEHPForceBalance2::getOffRailConnectedEdgeIds(
      igraph_integer_t vertexId,
      igraph_integer_t railEdge)
    {
      INVALIDARG_EXP_IFN(igraph_cattribute_VAN(&this->graph,
                                               "type",
                                               vertexId) != this->splipLinkType,
                         "Can only search for rail around slip-links");

      // fetch the edges involved
      igraph_integer_t otherRailEdge = this->getOtherRailEdgeId(vertexId, railEdge);
      igraph_es_t selector;
      igraph_es_incident(&selector, vertexId, IGRAPH_ALL);
      igraph_eit_t iterator;
      igraph_eit_create(&this->graph, selector, &iterator);
      std::vector<igraph_integer_t> results;
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
    igraph_integer_t MEHPForceBalance2::getOtherRailEdgeId(igraph_integer_t vertexId,
                                                 igraph_integer_t railEdge)
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

      igraph_destroy(&subgraph);
      igraph_vector_destroy(&parentEdges);
      igraph_vector_int_destroy(&allEdgesOfParent);
      igraph_vector_int_destroy(&edgesOnPath);
      igraph_vector_int_destroy(&verticesOnPath);

      assert(foundResult);
      return result;
    }

    //----------------------------------------------------------------
    // MARK: Structural Adjustments
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

      size_t newEdgeId = igraph_ecount(&this->graph);
      igraph_add_edge(&this->graph, neighbours[0], neighbours[1]);
      // verify our assumption of the new edge id
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
     * @brief Remove a certain, 3-functional link from the structures,
     * combining the two strands
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
    // MARK: Concrete structure modification
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
      size_t numRemovedTotal = this->removeSubfunctionalVertices();
      // and remove all springs that are inactive
      size_t numSpringsRemoved = this->removeInactiveParentEdges(tolerance);

      if (numSpringsRemoved > 0) {
        numRemovedTotal +=
          this->doRemovalAndreisWay(springPartitions, tolerance);
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
      if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
        this->updateGraph();
      }
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
          &this->graph, "type", numVerticesBefore + i, this->splipLinkType);
        // check that our assumption regarding the vertex id is correct
        igraph_integer_t degree;
        igraph_degree_1(
          &this->graph, &degree, numVerticesBefore + i, IGRAPH_ALL, true);
        assert(degree == 0);
      }

      igraph_vector_int_t edges_that_have_been_replaced;
      igraph_vector_int_init(&edges_that_have_been_replaced, 0);
      igraph_vector_int_reserve(&edges_that_have_been_replaced,
                                2 * additionalLen);
      for (size_t i = 0; i < additionalLen; ++i) {
        igraph_integer_t spring1ToReplace =
          this->findPartialSpringByFraction(strandIdx1[i], alpha1[i]);
        igraph_integer_t spring2ToReplace =
          this->findPartialSpringByFraction(strandIdx2[i], alpha2[i]);
        igraph_integer_t from1, to1, from2, to2;
        igraph_edge(&this->graph, spring1ToReplace, &from1, &to1);
        igraph_edge(&this->graph, spring2ToReplace, &from2, &to2);
        // add the four new edges
        igraph_integer_t newEdgeId = igraph_ecount(&this->graph);
        igraph_add_edge(&this->graph, from1, numVerticesBefore + i);
        double currentPartition1 = igraph_cattribute_EAN(
          &this->graph, "partition_fraction", spring1ToReplace);
        igraph_cattribute_EAN_set(&this->graph,
                                  "partition_fraction",
                                  newEdgeId,
                                  alpha1[i] * currentPartition1);
        igraph_cattribute_EAN_set(
          &this->graph,
          "parent_edge",
          newEdgeId,
          igraph_cattribute_EAN(&this->graph, "parent_edge", spring1ToReplace));
        // for the box offsets, we can simply give them to one of the two new
        // edges – all the rest should be handled by the optimisation, actually.
        for (std::string dir : { "x", "y", "z" }) {
          igraph_cattribute_EAN_set(
            &this->graph,
            ("bond_box_" + dir).c_str(),
            newEdgeId,
            igraph_cattribute_EAN(
              &this->graph, ("bond_box_" + dir).c_str(), spring1ToReplace));
        }
        newEdgeId += 1;
        igraph_add_edge(&this->graph, to1, numVerticesBefore + i);
        igraph_cattribute_EAN_set(&this->graph,
                                  "partition_fraction",
                                  newEdgeId,
                                  (1. - alpha1[i]) * currentPartition1);
        igraph_cattribute_EAN_set(
          &this->graph,
          "parent_edge",
          newEdgeId,
          igraph_cattribute_EAN(&this->graph, "parent_edge", spring1ToReplace));
        // ... and the other gets 0s as bond box offset.
        for (std::string dir : { "x", "y", "z" }) {
          igraph_cattribute_EAN_set(
            &this->graph, ("bond_box_" + dir).c_str(), newEdgeId, 0.0);
        }
        // second (third and fourth new) edge
        double currentPartition2 = igraph_cattribute_EAN(
          &this->graph, "partition_fraction", spring2ToReplace);
        newEdgeId += 1;
        igraph_add_edge(&this->graph, from2, numVerticesBefore + i);
        igraph_cattribute_EAN_set(&this->graph,
                                  "partition_fraction",
                                  newEdgeId,
                                  alpha2[i] * currentPartition2);
        igraph_cattribute_EAN_set(
          &this->graph,
          "parent_edge",
          newEdgeId,
          igraph_cattribute_EAN(&this->graph, "parent_edge", spring2ToReplace));
        for (std::string dir : { "x", "y", "z" }) {
          igraph_cattribute_EAN_set(
            &this->graph,
            ("bond_box_" + dir).c_str(),
            newEdgeId,
            igraph_cattribute_EAN(
              &this->graph, ("bond_box_" + dir).c_str(), spring2ToReplace));
        }
        newEdgeId += 1;
        igraph_add_edge(&this->graph, to2, numVerticesBefore + i);
        igraph_cattribute_EAN_set(&this->graph,
                                  "partition_fraction",
                                  newEdgeId,
                                  (1. - alpha2[i]) * currentPartition2);
        igraph_cattribute_EAN_set(
          &this->graph,
          "parent_edge",
          newEdgeId,
          igraph_cattribute_EAN(&this->graph, "parent_edge", spring2ToReplace));
        for (std::string dir : { "x", "y", "z" }) {
          igraph_cattribute_EAN_set(
            &this->graph, ("bond_box_" + dir).c_str(), newEdgeId, 0.0);
        }

        // mark the old ones to be removed
        igraph_vector_int_push_back(&edges_that_have_been_replaced,
                                    spring1ToReplace);
        igraph_vector_int_push_back(&edges_that_have_been_replaced,
                                    spring2ToReplace);
      }

      // delete the old edges, which have been "split up"
      igraph_delete_edges(&this->graph,
                          igraph_ess_vector(&edges_that_have_been_replaced));
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
    MEHPForceBalance2::getEffectiveFunctionalityOfAtoms(double tolerance) const
    {
      assert(this->net.isUpToDate);
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
      double tolerance) const
    {
      assert(this->net.isUpToDate);
      Eigen::VectorXi nrOfActiveSpringsConnected =
        Eigen::VectorXi::Zero(this->net.nrOfNodes);
      ArrayXb springIsActive =
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
      double tolerance) const
    {
      assert(this->net.isUpToDate);
      Eigen::VectorXi nrOfActivePartialSpringsConnected =
        Eigen::VectorXi::Zero(this->net.nrOfNodes);
      ArrayXb springIsActive =
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

    /**
     * @brief Estimate some sets of random vertices
     *
     * This is particularly useful (used) for parallelising the displacements.
     * Note that the returned Eigen::ArrayXi will contain the indices to the
     * coordinates rather than the actual vertex indices.
     *
     * @param net
     * @return std::vector<Eigen::ArrayXi>
     */
    std::vector<Eigen::ArrayXi> MEHPForceBalance2::getRandomCoordinateSets(
      const ForceBalanceNetwork& net) const
    {
      assert(this->net.isUpToDate);
      std::vector<Eigen::ArrayXi> results;
      results.reserve(2);
      Eigen::ArrayXf randomFloats = Eigen::ArrayXf::Random(net.nrOfLinks);
      size_t nrOfPositiveFloats = (randomFloats > 0).count();
      // TODO: implement something better
      Eigen::ArrayXi links1 = Eigen::ArrayXi(nrOfPositiveFloats * 3);
      size_t links1Idx = 0;
      Eigen::ArrayXi links2 =
        Eigen::ArrayXi((net.nrOfLinks - nrOfPositiveFloats) * 3);
      size_t links2Idx = 0;
      for (size_t i = 0; i < net.nrOfLinks; ++i) {
        if (randomFloats[i] > 0) {
          links1[links1Idx * 3] = 3 * i;
          links1[links1Idx * 3 + 1] = 3 * i + 1;
          links1[links1Idx * 3 + 2] = 3 * i + 2;
          links1Idx += 1;
        } else {
          links2[links2Idx * 3] = 3 * i;
          links2[links2Idx * 3 + 1] = 3 * i + 1;
          links2[links2Idx * 3 + 2] = 3 * i + 2;
          links2Idx += 1;
        }
      }
      results.push_back(links1);
      results.push_back(links2);
      return results;
    }

    /**
     * @brief Estimate some sets of independent vertices
     *
     * This is useful (used) for parallelising the displacements.
     * Note that the returned Eigen::ArrayXi will contain the indices to the
     * coordinates rather than the actual vertex indices.
     *
     * Time complexity: ca. O(|v||e|)
     *
     * @param net
     * @return std::vector<Eigen::ArrayXi>
     */
    std::pair<std::vector<Eigen::ArrayXi>, std::vector<Eigen::ArrayXi>>
    MEHPForceBalance2::getHeuristicallyIndependentCoordinateSets(
      const ForceBalanceNetwork& net) const
    {
      assert(this->net.isUpToDate);
      // std::vector<Eigen::ArrayXi> results = { net.springCoordinateIndexA,
      //                                         net.springCoordinateIndexB };
      std::vector<Eigen::ArrayXi> resultingCoordinateIndexMask;
      std::vector<Eigen::ArrayXi> involvedSpringPartCoordinateIndexMask;
      // global block list: block all indices that are added to a result already
      ArrayXb globalBlocked = ArrayXb::Constant(net.nrOfLinks, false);
      size_t remainingLinks = net.nrOfLinks;
      size_t globalStartingIdx = 0;
      while (remainingLinks > 0) {
        ArrayXb localBlocked = globalBlocked;
        std::vector<size_t> localIndexList;
        std::vector<size_t> localSpringIndexList;

        size_t localStartingIndex = globalStartingIdx;
        while (localStartingIndex < net.nrOfLinks) {
          while (localBlocked[localStartingIndex]) {
            localStartingIndex += 1;
            if (!(localStartingIndex < net.nrOfLinks)) {
              goto while2exit;
            }
          }
          // add this link to the current results list
          localIndexList.push_back(localStartingIndex);
          remainingLinks -= 1;
          localBlocked[localStartingIndex] = true;
          globalBlocked[localStartingIndex] = true;
          // block the neighbours
          std::vector<size_t> connections =
            net.springIndicesOfLinks[localStartingIndex];
          for (size_t springIdx : connections) {
            localSpringIndexList.push_back(springIdx);
            std::vector<size_t> springPartners =
              net.linkIndicesOfSprings[springIdx];
            for (size_t i = 0; i < springPartners.size(); i++) {
              if (springPartners[i] == localStartingIndex) {
                if (i > 0) {
                  localBlocked[springPartners[i - 1]] = true;
                }
                if (i < springPartners.size() - 1) {
                  localBlocked[springPartners[i + 1]] = true;
                }
                break;
              }
            }
          }
        }
      while2exit:

        while (globalStartingIdx < net.nrOfLinks &&
               globalBlocked[globalStartingIdx]) {
          globalStartingIdx += 1;
        }

        // translate the localIndexList to the results
        Eigen::ArrayXi localRes =
          Eigen::ArrayXi::Zero(3 * localIndexList.size());
        for (int i = 0; i < localIndexList.size(); i++) {
          localRes.segment(3 * i, 3) << 3 * localIndexList[i],
            3 * localIndexList[i] + 1, 3 * localIndexList[i] + 2;
        }
        resultingCoordinateIndexMask.push_back(localRes);
        Eigen::ArrayXi localSpringRes =
          Eigen::ArrayXi::Zero(3 * localSpringIndexList.size());
        for (int i = 0; i < localSpringIndexList.size(); i++) {
          localSpringRes.segment(3 * i, 3) << 3 * localSpringIndexList[i],
            3 * localSpringIndexList[i] + 1, 3 * localSpringIndexList[i] + 2;
        }
        involvedSpringPartCoordinateIndexMask.push_back(localSpringRes);
      }

      return std::make_pair(resultingCoordinateIndexMask,
                            involvedSpringPartCoordinateIndexMask);
    }

    //----------------------------------------------------------------
    // MARK: Computations
    //----------------------------------------------------------------

    /**
     * @brief Compute the force acting on a cross-link
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
    Eigen::Matrix3d MEHPForceBalance2::evaluateForceOnCrossLink(
      const size_t linkIdx,
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& springPartitions,
      Eigen::VectorXi& debugNrSpringsVisited,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {
      assert(net.isUpToDate);
      INVALIDARG_EXP_IFN(!net.linkIsSliplink[linkIdx],
                         "This link is not a cross-link");
      std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
      Eigen::Matrix3d force = Eigen::Matrix3d::Zero();
      for (size_t springIndex : springIndices) {
        std::vector<size_t> springsPartners =
          net.linkIndicesOfSprings[springIndex];
        // special case for primary loops
        std::vector<size_t> partnerIndices = { 0, springsPartners.size() - 1 };
        for (size_t partner_idx : partnerIndices) {
          if (springsPartners[partner_idx] != linkIdx) {
            continue;
          }
          size_t springGlobalIdx;
          if (partner_idx == 0) {
            springGlobalIdx = net.localToGlobalSpringIndex[springIndex][0];
          } else {
            RUNTIME_EXP_IFN(springsPartners[springsPartners.size() - 1] ==
                              linkIdx,
                            "Unexpected state. Cross-link must be either first "
                            "or last in springs' contributors.");
            springGlobalIdx = pylimer_tools::utils::last<size_t>(
              net.localToGlobalSpringIndex[springIndex]);
          }
          Eigen::Vector3d distance = this->evaluatePartialSpringDistanceFrom(
            net, springGlobalIdx, linkIdx, this->is2D);
          debugNrSpringsVisited[springGlobalIdx] += 1;
          const double N = net.springsContourLength[springIndex];
          double denominator = 1. / (springPartitions[springGlobalIdx] * N);
          if (oneOverSpringPartitionUpperLimit > 0 ||
              !std::isfinite(denominator)) {
            denominator = CLAMP_ONE_OVER_SPRINGPARTITION(
              net.partialSpringIsPartial[springGlobalIdx],
              denominator,
              N,
              oneOverSpringPartitionUpperLimit);
          }

          for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
              force(i, j) += kappa0 * denominator * distance[i] * distance[j];
            }
          }
        }
      }
      return force;
    };

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
    Eigen::Matrix3d MEHPForceBalance2::evaluateForceOnSlipLink(
      const size_t linkIdx,
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& springPartitions,
      Eigen::VectorXi& debugNrSpringsVisited,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {
      assert(net.isUpToDate);
      INVALIDARG_EXP_IFN(net.linkIsSliplink[linkIdx],
                         "This link is not a slip-link");

      std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
      Eigen::Matrix3d force = Eigen::Matrix3d::Zero();

      for (size_t springIndex : springIndices) {
        std::vector<size_t> springsPartners =
          net.linkIndicesOfSprings[springIndex];
        for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          // TODO: think what happens when multiple on same strand
          if (springsPartners[partner_idx] == linkIdx) {
            size_t backSpringGlobalIdx =
              net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
            size_t forwardSpringGlobalIdx =
              net.localToGlobalSpringIndex[springIndex][partner_idx];
            Eigen::Vector3d vecBack = this->evaluatePartialSpringDistanceFrom(
              net, backSpringGlobalIdx, linkIdx, this->is2D);
            Eigen::Vector3d vecForward =
              this->evaluatePartialSpringDistanceFrom(
                net, forwardSpringGlobalIdx, linkIdx, this->is2D);

            debugNrSpringsVisited[backSpringGlobalIdx] += 1;
            debugNrSpringsVisited[forwardSpringGlobalIdx] += 1;

            const double N = net.springsContourLength[springIndex];
            double denominatorBack =
              1. / (springPartitions[backSpringGlobalIdx] * N);
            if (oneOverSpringPartitionUpperLimit > 0 ||
                !std::isfinite(denominatorBack)) {
              denominatorBack = CLAMP_ONE_OVER_SPRINGPARTITION(
                true, denominatorBack, N, oneOverSpringPartitionUpperLimit);
            }

            double denominatorForward =
              1 / (springPartitions[forwardSpringGlobalIdx] *
                   net.springsContourLength[springIndex]);
            if (oneOverSpringPartitionUpperLimit > 0 ||
                !std::isfinite(denominatorForward)) {
              denominatorForward = CLAMP_ONE_OVER_SPRINGPARTITION(
                true, denominatorForward, N, oneOverSpringPartitionUpperLimit);
            }
            for (size_t i = 0; i < 3; ++i) {
              for (size_t j = 0; j < 3; ++j) {
                force(i, j) +=
                  kappa0 * denominatorBack * vecBack[i] * vecBack[j];
                force(i, j) +=
                  kappa0 * denominatorForward * vecForward[i] * vecForward[j];
              }
            }
          }
        }
      }
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
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& springPartitions,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit,
      const bool xlinksOnly) const
    {
      Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();
      INVALIDARG_EXP_IFN(
        springPartitions.size() == net.springPartIndexA.size(),
        "Spring partitions must have the size of partial springs.");

      double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

      Eigen::VectorXi debugNrSpringsVisited =
        Eigen::VectorXi::Zero(net.nrOfPartialSprings);

      size_t nrOfLinksToInspect = xlinksOnly ? net.nrOfNodes : net.nrOfLinks;
      for (size_t linkIdx = 0; linkIdx < nrOfLinksToInspect; ++linkIdx) {
        Eigen::Matrix3d force =
          net.linkIsSliplink[linkIdx]
            ? this->evaluateForceOnSlipLink(linkIdx,
                                            net,
                                            springPartitions,
                                            debugNrSpringsVisited,
                                            kappa0,
                                            oneOverSpringPartitionUpperLimit)
            : this->evaluateForceOnCrossLink(linkIdx,
                                             net,
                                             springPartitions,
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
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& springPartitions,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {
      assert(net.isUpToDate);
      std::array<std::array<double, 3>, 3> stress;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          stress[i][j] = 0.0;
        }
      }
      INVALIDARG_EXP_IFN(
        springPartitions.size() == net.springPartIndexA.size(),
        "Spring partitions must have the size of partial springs.");

      double oneOverVolume = 1. / (net.L[0] * net.L[1] * net.L[2]);

      Eigen::VectorXd displacedCoords = net.coordinates;
      Eigen::VectorXd relevantPartialDistancesA =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;
      if (this->assumeBoxLargeEnough) {
        this->universe.getBox().handlePBC(relevantPartialDistancesA);
      }

      for (size_t partialSpringIdx = 0;
           partialSpringIdx < net.nrOfPartialSprings;
           ++partialSpringIdx) {
        // Eigen::VectorXd force =
        //   net.linkIsSliplink[i]
        //     ? this->evaluateForceOnSlipLink(
        //         i, net, u, springPartitions, kappa0, minCutoff)
        //     : this->evaluateForceOnCrossLink(
        //         i, net, u, springPartitions, kappa0, minCutoff);
        Eigen::Vector3d distance =
          relevantPartialDistancesA.segment(3 * partialSpringIdx, 3);
        size_t totalSpringIndex =
          net.partialToFullSpringIndex[partialSpringIdx];
        const double N = net.springsContourLength[totalSpringIndex];
        double denominator = 1. / (springPartitions[partialSpringIdx] * N);
        if (oneOverSpringPartitionUpperLimit > 0. ||
            !std::isfinite(denominator)) {
          denominator = CLAMP_ONE_OVER_SPRINGPARTITION(
            net.partialSpringIsPartial[partialSpringIdx],
            denominator,
            N,
            oneOverSpringPartitionUpperLimit);
        }
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
                " for partial spring " + std::to_string(partialSpringIdx) +
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
    Eigen::VectorXd MEHPForceBalance2::evaluateSpringDistances(
      const ForceBalanceNetwork& net,
      const bool is2D) const
    {
      assert(net.isUpToDate);
      Eigen::VectorXd partialSpringDistances =
        this->evaluatePartialSpringDistances(net, is2D);
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
    Eigen::VectorXd MEHPForceBalance2::evaluatePartialSpringDistances(
      const ForceBalanceNetwork& net,
      const bool is2D) const
    {
      assert(net.isUpToDate);

      Eigen::VectorXd springDistances =
        (net.coordinates(net.springPartCoordinateIndexB) -
         net.coordinates(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;
      assert(springDistances.size() == net.nrOfPartialSprings * 3);
      if (this->assumeBoxLargeEnough) {
        this->universe.getBox().handlePBC(springDistances);
      }
      // reset for 2D systems
      if (is2D && net.nrOfPartialSprings > 0) {
        // springDistances(Eigen::seq(2, net.nrOfSprings, 3)) =
        //   Eigen::VectorXd::Zero(net.nrOfSprings);
        for (size_t i = 2; i < 3 * net.nrOfSprings; i += 3) {
          springDistances[i] = 0.0;
        }
      }

      return springDistances;
    }

    /**
     * @brief Get the average spring length at the current step
     *
     * @return double
     */
    double MEHPForceBalance2::getAverageSpringLength() const
    {
      assert(net.isUpToDate);
      double r2 = 0.0;
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
     * @brief Translate the spring partition vector to its 3*size
     *
     * @param net
     * @param springPartitions0
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance2::assembleOneOverSpringPartition(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& springPartitions0,
      const double oneOverSpringPartitionUpperLimit) const
    {
      INVALIDARG_EXP_IFN(
        springPartitions0.size() == net.nrOfPartialSprings,
        "Spring partitions must have the size of the nr of springs");
      Eigen::VectorXd oneOverSpringPartitions =
        Eigen::VectorXd(3 * net.nrOfPartialSprings);

      Eigen::ArrayXd primaryLoopCorrectionMultiplier =
        (net.springPartIndexA != net.springPartIndexB)
          .cast<double>(); // 0.0 for equal = primary loop, 1.0 otherwise

      for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
        const double N =
          net.springsContourLength[net.partialToFullSpringIndex[i]];
        double valueToSet = (springPartitions0[i] > 0.0) // 1e-18 //
                              ? 1.0 / (springPartitions0[i] * N)
                              : 0.0;
        valueToSet =
          CLAMP_ONE_OVER_SPRINGPARTITION(net.partialSpringIsPartial[i],
                                         valueToSet,
                                         N,
                                         oneOverSpringPartitionUpperLimit);

        // if (springPartitions0[i] < 1e-9) {
        //   std::cout << "Got close call for partial spring " << i <<
        //   std::endl;
        // }
        oneOverSpringPartitions.segment(3 * i, 3) = Eigen::Vector3d::Constant(
          valueToSet * primaryLoopCorrectionMultiplier[i]);
      }

      return oneOverSpringPartitions;
    }

    /**
     * @brief Compute the displacement residual
     *
     * @param net
     * @param oneOverSpringPartitions
     * @return double
     */
    double MEHPForceBalance2::getDisplacementResidualNormFor(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& oneOverSpringPartitions) const
    {
      Eigen::VectorXd relevantPartialDistancesA =
        (net.coordinates(net.springPartCoordinateIndexB) -
         net.coordinates(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;
      for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
        for (size_t dir = 0; dir < 3; ++dir) {
          if (std::abs(relevantPartialDistancesA[3 * i + dir]) >
              50. * net.L[dir]) {
            std::cerr
              << "WARNING: Spring " << i << " between "
              << net.springPartIndexA[i] << " and " << net.springPartIndexB[i]
              << " has a length of " << relevantPartialDistancesA[3 * i + dir]
              << " in dir " << dir << " from "
              << net.coordinates[net.springPartCoordinateIndexB[3 * i + dir]]
              << " minus "
              << net.coordinates[net.springPartCoordinateIndexA[3 * i + dir]]
              << " meaning it will probably fail in PBC." << std::endl;
          }
        }
      }
      if (this->assumeBoxLargeEnough) {
        this->universe.getBox().handlePBC(relevantPartialDistancesA);
      }
      Eigen::VectorXd partialDistancesOverSpringPartitions =
        (relevantPartialDistancesA.array() * oneOverSpringPartitions.array())
          .matrix();

      // return partialDistancesOverSpringPartitions.squaredNorm();

      Eigen::VectorXd overallForces = Eigen::VectorXd::Zero(3 * net.nrOfLinks);
      overallForces(net.springPartCoordinateIndexB) +=
        partialDistancesOverSpringPartitions;
      overallForces(net.springPartCoordinateIndexA) -=
        partialDistancesOverSpringPartitions;

      // Eigen::Index maxRow, maxCol;
      // double max = overallForces.maxCoeff(&maxRow, &maxCol);
      // Eigen::Index minRow, minCol;
      // double min = overallForces.minCoeff(&minRow, &minCol);
      // std::cout << "Got min = " << min << " at " << minRow << ", " <<
      // minCol
      //           << std::endl;
      // std::cout << "Got max = " << max << " at " << maxRow << ", " <<
      // maxCol
      //           << std::endl;

      return overallForces.squaredNorm();

      // Eigen::Vector3d sumOfResiduals = Eigen::Vector3d::Zero();
      // for (int i = 0; i < net.nrOfPartialSprings; ++i) {
      //   sumOfResiduals += partialDistancesOverSpringPartitions.segment(3 *
      //   i, 3).cwiseAbs2();
      // }
      // return sumOfResiduals.squaredNorm();
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
      bool usePartial) const
    {
      assert(net.isUpToDate);
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
      const
    {
      assert(this->net.isUpToDate);
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
        types[i] = this->universe.getPropertyValue<int>(
          "type", this->universe.getIdxByAtomId(this->net.oldAtomIds[i]));
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
      const Eigen::VectorXd& springPartitions) const
    {
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
        "Nr of nodes plus nr of slp-links should give the total nr of links");
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
        for (std::vector<size_t> loopsOfSliplink : net.loopsOfSliplink) {
          RUNTIME_EXP_IFN(
            loopsOfSliplink.size() <= 2,
            "Cannot have a slip-link attributed to more than two loops.");
          for (size_t loopIdx : loopsOfSliplink) {
            RUNTIME_EXP_IFN(loopIdx < net.loops.size(),
                            "Loop index out of range.");
          }
        }
        for (std::vector<size_t> loop : net.loops) {
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
