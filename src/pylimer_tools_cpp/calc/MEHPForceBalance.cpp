#include "MEHPForceBalance.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlopt.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {
/**
 * @brief a macro for doing the clamping in the routines using kappa,
 * to prevent deivision by zero issues / multiplications by infinity
 */
#define CLAMP_ONE_OVER_SPRINGPARTITION(                                        \
  val, N, oneOverSpringPartitionUpperLimit)                                    \
  std::clamp(val,                                                              \
             oneOverSpringPartitionUpperLimit > 0.                             \
               ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)              \
               : 0.0,                                                          \
             oneOverSpringPartitionUpperLimit > 0.                             \
               ? oneOverSpringPartitionUpperLimit                              \
               : N);

    /**
     * FORCE RELAXATION
     */
    void MEHPForceBalance::runForceRelaxation(
      BalanceRunMode mode,
      double damping,
      long int maxNrOfSteps, // default: 10000
      double xtol,
      long int innerMaxNrOfSteps,
      double innerAlphaTol,
      const double oneOverSpringPartitionUpperLimit,
      const int maxFlag,
      const bool allowRemovalOfSliplinks,
      const bool allowMoveOfSliplinks)
    {
      this->simulationHasRun = true;

      ForceBalanceNetwork net = this->initialConfig;
      const int M = this->universe.getMolecules(crosslinkerType).size();
      const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
      const double bM = this->universe.computeMeanBondLength();
      const int f =
        this->universe.determineFunctionalityPerType()[crosslinkerType];
      bool is2D = this->is2D;

      /* array allocation */
      Eigen::VectorXd u = this->currentDisplacements;
      Eigen::VectorXd springPartitions = this->currentSpringPartitionsVec;
      std::vector<Eigen::ArrayXi> independentVertexSets;
      double maxDistanceMoved = 0.0;
      size_t indexOfMaxDistanceMoved = 0;
      size_t iterationsDone = 0;
      size_t totalInnerIterationsDone = 0;
      // default = all
      std::vector<Eigen::ArrayXi> independentVertexsSpringSets;
      if (mode == BalanceRunMode::EIGEN_HEURISTIC) {
        std::tie(independentVertexSets, independentVertexsSpringSets) =
          getHeuristicallyIndependentCoordinateSets(net);
      } else if (mode == BalanceRunMode::EIGEN_RANDOM) {
        independentVertexSets = this->getRandomCoordinateSets(net);
      } else if (mode == BalanceRunMode::EIGEN_STRANDS) {
        independentVertexSets = { net.springPartCoordinateIndexA,
                                  net.springPartCoordinateIndexB };
      } else if (mode == BalanceRunMode::EIGEN_ALL) {
        independentVertexSets = { Eigen::ArrayXi::LinSpaced(
          3 * net.nrOfLinks, 0, 3 * net.nrOfLinks - 1) };
      }
      // this->getIndependentCoordinateSets(net);
      // { net.springPartCoordinateIndexA, net.springPartCoordinateIndexB };
      std::cout << "Starting force balance procedure "
                // "with " << independentVertexSets.size() << "vertex sets."
                << std::endl;
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(
          net, springPartitions, oneOverSpringPartitionUpperLimit);
      double initialResidual =
        this->getDisplacementResidualNormFor(net, u, oneOverSpringPartitions);
      double currentResidual = 0.0;
      double intermediateResidual = 0.0;
      do {
        maxDistanceMoved = 0.0;
        currentResidual = 0.0;
        // place slip-link
        for (size_t link_idx = net.nrOfNodes; link_idx < net.nrOfLinks;
             ++link_idx) {
          std::vector<size_t> relevantPartitionIndices =
            this->getSpringpartitionIndicesOfSliplink(link_idx);
          assert(relevantPartitionIndices.size() == 4);
          size_t innerIterationsDone = 0;
          double displacementDone = 0.0;
          double rOverr0 = 0.0;
          double r2 = 0.0;
          double r02 = this->computePartitionUpdateZeroResidual(
            link_idx, u, springPartitions, oneOverSpringPartitionUpperLimit);

          bool allAtEnd = false;
          int flags = 0;
          do {
            r2 = this->updateSpringPartition(net,
                                             u,
                                             springPartitions,
                                             link_idx,
                                             oneOverSpringPartitionUpperLimit);
            rOverr0 = r2 / r02;
            displacementDone =
              this->displaceToMeanPosition(net,
                                           u,
                                           springPartitions,
                                           link_idx,
                                           oneOverSpringPartitionUpperLimit);
            innerIterationsDone += 1;
            flags += (springPartitions(relevantPartitionIndices).array() <
                      1 / net.meanSpringContourLength)
                       .count();
          } while (innerIterationsDone < innerMaxNrOfSteps &&
                   rOverr0 > innerAlphaTol && std::isfinite(rOverr0) &&
                   flags < maxFlag);
          totalInnerIterationsDone += innerIterationsDone;
        }
        oneOverSpringPartitions = this->assembleOneOverSpringPartition(
          net, springPartitions, oneOverSpringPartitionUpperLimit);
        intermediateResidual = this->getDisplacementResidualNormFor(
          net, u, oneOverSpringPartitions);

        if (mode == BalanceRunMode::EIGEN_RANDOM) {
          independentVertexSets = getRandomCoordinateSets(net);
        }
        // place cross-links
        if (mode == BalanceRunMode::ITERATIVE) {
          for (size_t link_idx = 0; link_idx < net.nrOfNodes; ++link_idx) {
            double distanceMoved =
              this->displaceToMeanPosition(net,
                                           u,
                                           springPartitions,
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
                this->displaceLinksToMeanPosition(net,
                                                  u,
                                                  oneOverSpringPartitions,
                                                  independentVertexsSpringSet,
                                                  vertexSet,
                                                  damping));
            } else {
              maxDistanceMoved = std::max(
                maxDistanceMoved,
                this->displaceLinksToMeanPosition(
                  net, u, oneOverSpringPartitions, vertexSet, damping));
            }
          }
        }
        oneOverSpringPartitions = this->assembleOneOverSpringPartition(
          net, springPartitions, oneOverSpringPartitionUpperLimit);
        currentResidual = this->getDisplacementResidualNormFor(
          net, u, oneOverSpringPartitions);
        iterationsDone += 1;
        if (iterationsDone % 50 == 0) {
          std::cout << "Iteration " << iterationsDone << " " << maxDistanceMoved
                    << " by " << indexOfMaxDistanceMoved
                    << ". Residual: " << currentResidual
                    << " from: " << initialResidual << " via "
                    << intermediateResidual << "\n";
          std::array<std::array<double, 3>, 3> stressTensor =
            this->evaluateStressTensor(
              net, u, springPartitions, 1.0, oneOverSpringPartitionUpperLimit);
          std::cout << "To stress tensor diagonal: " << stressTensor[0][0]
                    << ", " << stressTensor[1][1] << ", " << stressTensor[2][2]
                    << " with ";
          std::cout << "Total inner: " << totalInnerIterationsDone << "\n";
        }
      } while (currentResidual / initialResidual > xtol &&
               iterationsDone < maxNrOfSteps);

      // query solution & exit reason
      assert(u.size() == 3 * net.nrOfLinks);
      this->currentDisplacements = u;
      this->currentSpringPartitionsVec = springPartitions;
      this->currentSpringDistances =
        this->evaluateSpringDistances(net, this->currentDisplacements, is2D);
      this->currentPartialSpringDistances =
        this->evaluatePartialSpringDistances(
          net, this->currentDisplacements, is2D);

      this->exitReason = iterationsDone == maxNrOfSteps
                           ? ExitReason::MAX_STEPS
                           : ExitReason::X_TOLERANCE;
      this->nrOfStepsDone += iterationsDone;
      std::cout << iterationsDone << " steps done, " << totalInnerIterationsDone
                << " inner iterations. Last max distance moved: "
                << maxDistanceMoved << std::endl;
      this->validateNetwork(net);
    }

    double MEHPForceBalance::displaceLinksToMeanPosition(
      const ForceBalanceNetwork &net,
      Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions0,
      double damping) const
    {

      Eigen::ArrayXi mask = Eigen::ArrayXi::LinSpaced(
        3 * net.nrOfLinks, 0, 3 * net.nrOfLinks - 1);
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(net, springPartitions0);

      return this->displaceLinksToMeanPosition(
        net, u, oneOverSpringPartitions, mask, damping);
    }

    /**
     * @brief Displace one link to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @param u the current displacements, wherein the resulting coordinates
     * shall be stored
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance::displaceLinksToMeanPosition(
      const ForceBalanceNetwork &net,
      Eigen::VectorXd& u,
      const Eigen::VectorXd& oneOverSpringPartitions,
      const Eigen::ArrayXi& resultingCoordinateIndexMask,
      const double damping) const
    {
      Eigen::ArrayXi involvedSpringPartCoordinateIndexMask =
        Eigen::ArrayXi::LinSpaced(
          3 * net.nrOfPartialSprings, 0, 3 * net.nrOfPartialSprings - 1);
      return this->displaceLinksToMeanPosition(
        net,
        u,
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
    double MEHPForceBalance::displaceLinksToMeanPosition(
      const ForceBalanceNetwork &net,
      Eigen::VectorXd& u,
      const Eigen::VectorXd& oneOverSpringPartitions,
      const Eigen::ArrayXi& involvedSpringPartCoordinateIndexMask,
      const Eigen::ArrayXi& resultingCoordinateIndexMask,
      const double damping) const
    {
      INVALIDARG_EXP_IFN(
        u.size() == net.coordinates.size(),
        "Coordinates and displacements must have the same size");
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
      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd relevantPartialDistancesA =
        (displacedCoords(relevantSpringPartCoordinateIndexB) -
         displacedCoords(relevantSpringPartCoordinateIndexA));
      this->handlePBC(net, relevantPartialDistancesA);

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
      double maxDiff = 0;
      Eigen::VectorXd objectivesToSet =
        objectiveDisplacements(resultingCoordinateIndexMask);
      for (int i = 0; i < resultingCoordinateIndexMask.size() / 3; i++) {
        maxDiff =
          std::max(maxDiff, objectivesToSet.segment(3 * i, 3).squaredNorm());
      }
      // double maxDiff =
      // (objectiveDisplacements(mask)).cwiseAbs2().maxCoeff();
      u(resultingCoordinateIndexMask) += damping * objectivesToSet;

      return maxDiff;
    };

    double MEHPForceBalance::getDisplacementResidualNorm(double cutoff) const
    {
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(
          this->initialConfig, this->currentSpringPartitionsVec, cutoff);
      Eigen::VectorXd displacements = this->currentDisplacements;
      return this->getDisplacementResidualNormFor(
        this->initialConfig, displacements, oneOverSpringPartitions);
    }

    double MEHPForceBalance::getDisplacementResidualNormFor(
      const ForceBalanceNetwork &net,
      Eigen::VectorXd& u,
      const Eigen::VectorXd& oneOverSpringPartitions) const
    {
      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd relevantPartialDistancesA =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA));
      this->handlePBC(net, relevantPartialDistancesA);
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
    std::vector<Eigen::ArrayXi> MEHPForceBalance::getRandomCoordinateSets(
      const ForceBalanceNetwork &net) const
    {
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
    MEHPForceBalance::getHeuristicallyIndependentCoordinateSets(
      const ForceBalanceNetwork &net) const
    {
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

    /**
     * @brief Build a graph of the current configuration and find all the sets
     * of independent vertices
     *
     * This is particularly useful (used) for parallelising the displacements.
     * Note that the returned Eigen::ArrayXi will contain the indices to the
     * coordinates rather than the actual vertex indices.
     *
     * @param net
     * @return std::vector<Eigen::ArrayXi>
     */
    std::vector<Eigen::ArrayXi> MEHPForceBalance::getIndependentCoordinateSets(
      const ForceBalanceNetwork &net) const
    {
      std::vector<Eigen::ArrayXi> results;
      // build graph
      igraph_t graph;
      igraph_empty(&graph, net.nrOfLinks, IGRAPH_UNDIRECTED);
      igraph_vector_int_t edges;
      igraph_vector_int_init(&edges, net.nrOfPartialSprings * 2);
      for (int i = 0; i < net.nrOfPartialSprings; ++i) {
        // igraph_vector_int_push_back(&edges, net.springPartIndexA[i]);
        // igraph_vector_int_push_back(&edges, net.springPartIndexB[i]);
        igraph_vector_int_set(&edges, 2 * i, net.springPartIndexA[i]);
        igraph_vector_int_set(&edges, 2 * i + 1, net.springPartIndexB[i]);
      }
      assert(igraph_vector_int_size(&edges) == net.nrOfPartialSprings * 2);
      igraph_add_edges(&graph, &edges, nullptr);
      igraph_vector_int_destroy(&edges);

      // find dependencies in graph
      igraph_vector_int_list_t dependencies;
      igraph_vector_int_list_init(&dependencies, 0);
      igraph_independent_vertex_sets(&graph, &dependencies, -1, -1);

      // assemble to results
      results.reserve(igraph_vector_int_list_size(&dependencies));
      for (size_t i = 0; i < igraph_vector_int_list_size(&dependencies); ++i) {
        igraph_vector_int_t* depsI =
          igraph_vector_int_list_get_ptr(&dependencies, i);
        Eigen::ArrayXi result = Eigen::ArrayXi(igraph_vector_int_size(depsI));
        for (size_t j = 0; j < igraph_vector_int_size(depsI); ++j) {
          igraph_integer_t resI = igraph_vector_int_get(depsI, j);
          for (size_t dir = 0; dir < 3; ++dir) {
            result[3 * j + dir] = 3 * resI + dir;
          }
        }
        results.push_back(result);
      }
      igraph_vector_int_list_destroy(&dependencies);

      return results;
    };

    /**
     * @brief Translate the spring partition vector to its 3*size
     *
     * @param net
     * @param springPartitions0
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance::assembleOneOverSpringPartition(
      const ForceBalanceNetwork &net,
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
          net.springsContourLength[net.partialToFullSpringIndex.at(i)];
        double valueToSet = springPartitions0[i] > 0.0 // 1e-18 //
                              ? 1.0 / (springPartitions0[i] * N)
                              : 0.0;
        valueToSet = CLAMP_ONE_OVER_SPRINGPARTITION(
          valueToSet, N, oneOverSpringPartitionUpperLimit);

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
     * @brief Assemble all indices of partial springs for a particular slip-link
     *
     * @param linkIdx
     * @return std::vector<size_t>
     */
    std::vector<size_t> MEHPForceBalance::getSpringpartitionIndicesOfSliplink(
      const size_t linkIdx) const
    {
      INVALIDARG_EXP_IFN(this->initialConfig.linkIsSliplink[linkIdx],
                         "Link must be slip-link");
      std::vector<size_t> springIndices =
        this->initialConfig.springIndicesOfLinks[linkIdx];
      std::vector<size_t> results;
      results.reserve(4);
      for (size_t springIndex : springIndices) {
        std::vector<size_t> springsPartners =
          this->initialConfig.linkIndicesOfSprings[springIndex];
        for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          if (springsPartners[partner_idx] == linkIdx) {
            size_t currentSpringGlobalIdx =
              this->initialConfig.localToGlobalSpringIndex.at(
                springIndex)[partner_idx - 1];
            size_t neighbourSpringGlobalIdx =
              this->initialConfig.localToGlobalSpringIndex.at(
                springIndex)[partner_idx];
            results.push_back(currentSpringGlobalIdx);
            results.push_back(neighbourSpringGlobalIdx);
          }
        }
      }
      return results;
    }

    void MEHPForceBalance::moveRemoveSlipLinks(
      const bool move,
      const bool remove,
      ForceBalanceNetwork &net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions,
      double tolerance)
    {
      std::vector<size_t> linksToRemove;
      for (size_t linkIdx = net.nrOfNodes; linkIdx < net.nrOfLinks;
           ++linkIdx) {
        std::vector<size_t> relevantPartitionIndices =
          this->getSpringpartitionIndicesOfSliplink(linkIdx);
        for (size_t partitonIdx : relevantPartitionIndices) {
          if (springPartitions[partitonIdx] < tolerance) {
            // possible move/remove.
            // decide!
            // have to find out with what we coincide
            size_t p1 = net.springPartIndexA[partitonIdx];
            size_t p2 = net.springPartIndexB[partitonIdx];
            assert(p1 == linkIdx || p2 == linkIdx);
            size_t relevantPartner = p1 == linkIdx ? p2 : p1;
            if (net.linkIsSliplink[relevantPartner]) {
              // NOTE: you see here, we do not let slip-links pass each other
              // yet
              continue;
            }
            // check functionality of cross-link
            // TODO: implement
          }
        }
      }
    }

    /**
     * @brief Updates the partition/parametrisation of a spring around one link
     *
     */
    double MEHPForceBalance::updateSpringPartition(
      const ForceBalanceNetwork &net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions, /* gives the parametrisation of N */
      const size_t linkIdx,
      double oneOverSpringPartitionUpperLimit) const
    {
      INVALIDARG_EXP_IFN(linkIdx < net.springIndicesOfLinks.size(),
                         "Link to update needs to be in the list");
      INVALIDARG_EXP_IFN(net.linkIsSliplink[linkIdx],
                         "Only slip-links may slip along a spring");
      std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
      double residualNorm = 0.0;
      for (size_t springIndex : springIndices) {
        std::vector<size_t> springsPartners =
          net.linkIndicesOfSprings[springIndex];
        for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          if (springsPartners[partner_idx] == linkIdx) {
            // found position of this link in this spring
            // want to find the ideal value for
            // net.springPartitions[springIndex][partner_idx-1]
            Eigen::Vector3d vecBack =
              (MEHPForceBalance::evaluateDistanceBetween(
                net,
                u,
                springsPartners[partner_idx],
                springsPartners[partner_idx - 1],
                this->is2D));
            double distanceBack = (vecBack.squaredNorm());
            Eigen::Vector3d vecForward =
              (MEHPForceBalance::evaluateDistanceBetween(
                net,
                u,
                springsPartners[partner_idx + 1],
                springsPartners[partner_idx],
                this->is2D));
            double distanceForward = vecForward.squaredNorm();
            double idealValue =
              1. / (1. + sqrt(distanceForward / distanceBack));
            if (distanceBack <= 0.0) {
              idealValue = 0.0; // TODO: really?
            }
            size_t currentSpringGlobalIdx =
              net.localToGlobalSpringIndex.at(springIndex)[partner_idx - 1];
            size_t neighbourSpringGlobalIdx =
              net.localToGlobalSpringIndex.at(springIndex)[partner_idx];
            double currentS = springPartitions[currentSpringGlobalIdx];
            double nextS = springPartitions[neighbourSpringGlobalIdx];
            double newS = idealValue * (nextS + currentS);
            double complementaryS = (1. - idealValue) * (nextS + currentS);
            double localResidualNorm = 0.0;
            if (complementaryS > 0.) {
              localResidualNorm +=
                distanceForward / (complementaryS * complementaryS);
            }
            if (newS > 0.) {
              localResidualNorm += distanceBack / (newS * newS);
            }

            RUNTIME_EXP_IFN(
              APPROX_WITHIN(newS + complementaryS, 0., 1., 1e-10),
              "Require newS + complementaryS to be within 0, 1, got " +
                std::to_string(newS + complementaryS) + " from " +
                std::to_string(newS) + " and " +
                std::to_string(complementaryS) +
                " with ideal = " + std::to_string(idealValue) + " of " +
                std::to_string(nextS + currentS) + " for link " +
                std::to_string(linkIdx) + ".");
            RUNTIME_EXP_IFN(
              APPROX_EQUAL(nextS + currentS, newS + complementaryS, 1e-10),
              "Require nextS + currentS == newS + complementaryS, got " +
                std::to_string(nextS + currentS) + " vs. " +
                std::to_string(newS + complementaryS) + " from " +
                std::to_string(nextS) + " and " + std::to_string(currentS) +
                ", " + std::to_string(newS) + " and " +
                std::to_string(complementaryS) + ".");
            RUNTIME_EXP_IFN(APPROX_WITHIN(nextS + currentS, 0., 1., 1e-10),
                            "Require nextS + currentS to be within 0, 1, got " +
                              std::to_string(nextS + currentS) + " from " +
                              std::to_string(nextS) + " and " +
                              std::to_string(currentS) + ".");

            // (complementaryS > residualNormSTolerance &&
            //  newS > residualNormSTolerance)
            //   ? ( -
            //      distanceBack / (newS * newS))
            //   : 0.0;
            // if (!(APPROX_EQUAL(newS, currentS, 0.2))) {
            //   std::cout << "Updating " << linkIdx << " to " << newS << " and
            //   "
            //             << complementaryS << " with global springs "
            //             << currentSpringGlobalIdx << " and "
            //             << neighbourSpringGlobalIdx << " from " << currentS
            //             << ", " << nextS << std::endl;
            // }
            // std::cout << "Contribution to " << linkIdx
            //           << " from global springs " << currentSpringGlobalIdx
            //           << " (" << springsPartners[partner_idx - 1] << ") "
            //           << vecBack[0] << ", " << vecBack[1] << ", " <<
            //           vecBack[2]
            //           << " and " << neighbourSpringGlobalIdx << " ("
            //           << springsPartners[partner_idx + 1] << ") "
            //           << vecForward[0] << ", " << vecForward[1] << ", "
            //           << vecForward[2] << "; "
            //           << " with " << currentS << ", " << nextS << std::endl;
            // std::cout << "Distances are " << distanceForward << ", "
            //           << distanceBack << " to get ideal value " << idealValue
            //           << " for " << (nextS) << " , " << currentS <<
            //           std::endl;
            residualNorm += localResidualNorm * localResidualNorm;
            springPartitions[currentSpringGlobalIdx] = newS;
            springPartitions[neighbourSpringGlobalIdx] = complementaryS;
          }
        }
      }
      return residualNorm;
    }

    /**
     * @brief Displace one link to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @param u the current displacements, wherein the resulting coordinates
     * shall be stored
     * @param linkIdx the idx of the link to displace
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance::displaceToMeanPosition(
      const ForceBalanceNetwork &net,
      Eigen::VectorXd& u,
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
            size_t globalSpringIndex = net.localToGlobalSpringIndex.at(
              springIndices[spring_index])[partner_idx];
            Eigen::Vector3d partialDistance;
            // add partial distance to the total distance
            if (springsPartners[partner_idx] == linkIdx) {
              partialDistance = MEHPForceBalance::evaluateDistanceBetween(
                net,
                u,
                springsPartners[partner_idx + 1],
                springsPartners[partner_idx],
                this->is2D);
            } else {
              assert(springsPartners[partner_idx + 1] == linkIdx);
              partialDistance = MEHPForceBalance::evaluateDistanceBetween(
                net,
                u,
                springsPartners[partner_idx],
                springsPartners[partner_idx + 1],
                this->is2D);
            }
            // add to displacement
            double contourLengthFraction = springPartitions[globalSpringIndex];
            // std::cout << "Contribution from " << springsPartners[partner_idx]
            //           << " to " << springsPartners[partner_idx + 1]
            //           << " with l = " << contourLengthFraction << " and N = "
            //           <<
            //           net.springsContourLength[springIndices[spring_index]]
            //           << ", partial distance " << partialDistance[0] << ", "
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
            if (oneOverSpringPartitionUpperLimit > 0.0) {
              oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
                oneOverContourLengthFraction,
                N,
                oneOverSpringPartitionUpperLimit);
            }
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
      u.segment(3 * linkIdx, 3) +=
        objectiveDisplacement / (objectiveDisplacementContributors == 0.0
                                   ? 1.0
                                   : objectiveDisplacementContributors);

      double dist = objectiveDisplacement.squaredNorm();
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
      return dist;
    }

    Eigen::Matrix3d MEHPForceBalance::evaluateForceOnSlipLink(
      const size_t linkIdx,
      const ForceBalanceNetwork &net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      Eigen::VectorXi& debugNrSpringsVisited,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {
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
            Eigen::Vector3d vecBack =
              (MEHPForceBalance::evaluateDistanceBetween(
                net,
                u,
                springsPartners[partner_idx - 1],
                springsPartners[partner_idx],
                this->is2D));
            double distanceBack = (vecBack.squaredNorm());
            Eigen::Vector3d vecForward =
              (MEHPForceBalance::evaluateDistanceBetween(
                net,
                u,
                springsPartners[partner_idx + 1],
                springsPartners[partner_idx],
                this->is2D));
            double distanceForward = vecForward.squaredNorm();
            size_t backSpringGlobalIdx =
              net.localToGlobalSpringIndex.at(springIndex)[partner_idx - 1];
            size_t forwardSpringGlobalIdx =
              net.localToGlobalSpringIndex.at(springIndex)[partner_idx];

            debugNrSpringsVisited[backSpringGlobalIdx] += 1;
            debugNrSpringsVisited[forwardSpringGlobalIdx] += 1;

            const double N = net.springsContourLength[springIndex];
            double denominatorBack =
              1. / (springPartitions[backSpringGlobalIdx] * N);
            if (oneOverSpringPartitionUpperLimit > 0 ||
                !std::isfinite(denominatorBack)) {
              denominatorBack = CLAMP_ONE_OVER_SPRINGPARTITION(
                denominatorBack, N, oneOverSpringPartitionUpperLimit);
            }

            double denominatorForward =
              1 / (springPartitions[forwardSpringGlobalIdx] *
                   net.springsContourLength[springIndex]);
            if (oneOverSpringPartitionUpperLimit > 0 ||
                !std::isfinite(denominatorForward)) {
              denominatorForward = CLAMP_ONE_OVER_SPRINGPARTITION(
                denominatorForward, N, oneOverSpringPartitionUpperLimit);
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
    }

    Eigen::Matrix3d MEHPForceBalance::evaluateForceOnCrossLink(
      const size_t linkIdx,
      const ForceBalanceNetwork &net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      Eigen::VectorXi& debugNrSpringsVisited,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {
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
          Eigen::Vector3d distance;
          size_t springGlobalIdx;
          if (partner_idx == 0) {
            distance = (MEHPForceBalance::evaluateDistanceBetween(
              net, u, springsPartners[1], springsPartners[0], this->is2D));
            springGlobalIdx = net.localToGlobalSpringIndex.at(springIndex)[0];
          } else {
            RUNTIME_EXP_IFN(springsPartners[springsPartners.size() - 1] ==
                              linkIdx,
                            "Unexpected state. Cross-link must be either first "
                            "or last in springs' contributors.");

            distance = (MEHPForceBalance::evaluateDistanceBetween(
              net,
              u,
              springsPartners[springsPartners.size() - 2],
              springsPartners[springsPartners.size() - 1],
              this->is2D));
            springGlobalIdx =
              net.localToGlobalSpringIndex.at(springIndex)
                [net.localToGlobalSpringIndex.at(springIndex).size() - 1];
          }
          debugNrSpringsVisited[springGlobalIdx] += 1;
          const double N = net.springsContourLength[springIndex];
          double denominator = 1. / (springPartitions[springGlobalIdx] * N);
          if (oneOverSpringPartitionUpperLimit > 0 ||
              !std::isfinite(denominator)) {
            denominator = CLAMP_ONE_OVER_SPRINGPARTITION(
              denominator, N, oneOverSpringPartitionUpperLimit);
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

    Eigen::Vector3d MEHPForceBalance::evaluateDistanceBetween(
      const ForceBalanceNetwork &net,
      const Eigen::VectorXd& u,
      const size_t linkIndexA,
      const size_t linkIndexB,
      const bool is2D) const
    {
      Eigen::Vector3d distances = net.coordinates.segment(3 * linkIndexA, 3) +
                                  u.segment(3 * linkIndexA, 3) -
                                  (net.coordinates.segment(3 * linkIndexB, 3) +
                                   u.segment(3 * linkIndexB, 3));

      // Possibly improvable PBC
      this->handlePBC<Eigen::Vector3d>(net, distances);

      if (is2D) {
        distances[2] = 0.0;
      }

      return distances;
    }

    Eigen::VectorXd MEHPForceBalance::evaluateSpringDistances(
      const ForceBalanceNetwork &net,
      const Eigen::VectorXd& u,
      const bool is2D) const
    {
      // first, the distances
      assert(u.size() == net.coordinates.size());

      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd springDistances =
        (displacedCoords(net.springCoordinateIndexB) -
         displacedCoords(net.springCoordinateIndexA));
      assert(springDistances.size() == net.nrOfSprings * 3);
      this->handlePBC(net, springDistances);

      // reset for 2D systems
      if (this->is2D && net.nrOfSprings > 0) {
        // springDistances(Eigen::seq(2, net.nrOfSprings, 3)) =
        //   Eigen::VectorXd::Zero(net.nrOfSprings);
        for (size_t i = 2; i < 3 * net.nrOfSprings; i += 3) {
          springDistances[i] = 0.0;
        }
      }

      return springDistances;
    }

    Eigen::VectorXd MEHPForceBalance::evaluatePartialSpringDistances(
      const ForceBalanceNetwork &net,
      const Eigen::VectorXd& u,
      const bool is2D) const
    {
      // first, the distances
      assert(u.size() == net.coordinates.size());

      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd partialDistances =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA));
      this->handlePBC(net, partialDistances);

      // reset for 2D systems
      if (this->is2D) {
        // partialDistances(Eigen::seq(2, net.nrOfPartialSprings, 3)) =
        //   Eigen::VectorXd::Zero(net.nrOfPartialSprings);
        for (size_t i = 2; i < 3 * net.nrOfPartialSprings; i += 3) {
          partialDistances[i] = 0.0;
        }
      }

      return partialDistances;
    }

    /**
     * FORCE RELAXATION DATA ACCESS
     */
    pylimer_tools::entities::Universe MEHPForceBalance::getCrosslinkerVerse(
      int newCrosslinkerType) const
    {
      // convert nodes & springs back to a universe
      pylimer_tools::entities::Universe xlinkUniverse =
        pylimer_tools::entities::Universe(this->universe.getBox());
      std::vector<long int> ids;
      std::vector<int> types = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, crosslinkerType);
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      std::vector<int> zeros = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, 0);
      ids.reserve(this->initialConfig.nrOfNodes);
      x.reserve(this->initialConfig.nrOfNodes);
      y.reserve(this->initialConfig.nrOfNodes);
      z.reserve(this->initialConfig.nrOfNodes);
      for (int i = 0; i < this->initialConfig.nrOfNodes; ++i) {
        x.push_back(this->initialConfig.coordinates[3 * i + 0] +
                    this->currentDisplacements[3 * i + 0]);
        y.push_back(this->initialConfig.coordinates[3 * i + 1] +
                    this->currentDisplacements[3 * i + 1]);
        z.push_back(this->initialConfig.coordinates[3 * i + 2] +
                    this->currentDisplacements[3 * i + 2]);
        ids.push_back(this->initialConfig.oldAtomIds[i]);
      }
      xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
      std::vector<long int> bondFrom;
      std::vector<long int> bondTo;
      bondFrom.reserve(this->initialConfig.nrOfSprings);
      bondTo.reserve(this->initialConfig.nrOfSprings);
      for (int i = 0; i < this->initialConfig.nrOfSprings; ++i) {
        bondFrom.push_back(
          this->initialConfig.oldAtomIds[this->initialConfig.springIndexA[i]]);
        bondTo.push_back(
          this->initialConfig.oldAtomIds[this->initialConfig.springIndexB[i]]);
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

    void MEHPForceBalance::addSlipLinks(const std::vector<size_t>& strandIdx1,
                                        const std::vector<size_t>& strandIdx2,
                                        const std::vector<double>& x,
                                        const std::vector<double>& y,
                                        const std::vector<double>& z,
                                        const std::vector<double>& alpha1,
                                        const std::vector<double>& alpha2,
                                        bool clampAlpha)
    {
      size_t additionalLen = strandIdx1.size();
      if (additionalLen == 0) {
        return;
      }
      // validate inputs
      size_t currentNrOfLinks = this->initialConfig.nrOfLinks;
      size_t currentNrOfPartialSprings = this->initialConfig.nrOfPartialSprings;
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
          strandIdx1[i] < this->initialConfig.nrOfSprings,
          "Invalid spring index " + std::to_string(strandIdx1[i]) +
            ", expected below " +
            std::to_string(this->initialConfig.nrOfSprings) + ".");
        INVALIDARG_EXP_IFN(
          strandIdx2[i] < this->initialConfig.nrOfSprings,
          "Invalid spring index " + std::to_string(strandIdx2[i]) +
            ", expected below " +
            std::to_string(this->initialConfig.nrOfSprings) + ".");
        INVALIDARG_EXP_IFN(APPROX_WITHIN(alpha1[i], 0.0, 1.0, 1e-12),
                           "Expected alpha within [0, 1], got " +
                             std::to_string(alpha1[i]) + ".");
        INVALIDARG_EXP_IFN(APPROX_WITHIN(alpha2[i], 0.0, 1.0, 1e-12),
                           "Expected alpha within [0, 1], got " +
                             std::to_string(alpha2[i]) + ".");
      }
      // actually start adding them
      this->initialConfig.nrOfLinks += additionalLen;
      // but first, indicate the resize
      this->initialConfig.springIndicesOfLinks.reserve(
        this->initialConfig.nrOfLinks);
      this->currentDisplacements.conservativeResize(
        3 * this->initialConfig.nrOfLinks);
      this->currentSpringPartitionsVec.conservativeResize(
        currentNrOfPartialSprings + 2 * additionalLen);
      this->initialConfig.springPartCoordinateIndexA.conservativeResize(
        3 * (currentNrOfPartialSprings + 2 * additionalLen));
      this->initialConfig.springPartCoordinateIndexB.conservativeResize(
        3 * (currentNrOfPartialSprings + 2 * additionalLen));
      this->initialConfig.springPartIndexA.conservativeResize(
        currentNrOfPartialSprings + 2 * additionalLen);
      this->initialConfig.springPartIndexB.conservativeResize(
        currentNrOfPartialSprings + 2 * additionalLen);
      this->initialConfig.linkIsSliplink.conservativeResize(
        this->initialConfig.nrOfLinks);
      this->initialConfig.coordinates.conservativeResize(
        3 * this->initialConfig.nrOfLinks);
      size_t partialSpringsAdded = 0;
      // then, loop the slip-links to add
      for (size_t i = 0; i < additionalLen; ++i) {
        // add the info that is straight-forward to add
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i] = x[i];
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i + 1] =
          y[i];
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i + 2] =
          z[i];
        this->initialConfig.linkIsSliplink[currentNrOfLinks + i] = true;
        std::vector<size_t> springIndices{ strandIdx1[i], strandIdx2[i] };
        std::vector<size_t> springIndicesOfLink =
          strandIdx1[i] == strandIdx2[i] ? std::vector<size_t>{ strandIdx1[i] }
                                         : springIndices;
        this->initialConfig.springIndicesOfLinks.push_back(springIndicesOfLink);
        // add to the springs
        int springIndexIndex = 0;
        for (size_t springIndex : springIndices) {
          std::vector<size_t> springParticipants =
            this->initialConfig.linkIndicesOfSprings[springIndex];
          double alpha = (springIndexIndex == 0 ? alpha1[i] : alpha2[i]);
          INVALIDARG_EXP_IFN(alpha >= 0.0 && alpha <= 1.0,
                             "alpha must be between 0 and 1, got " +
                               std::to_string(alpha) + ".");
          if (clampAlpha) {
            alpha = std::clamp(
              alpha,
              1 / (this->initialConfig.springsContourLength[springIndex]),
              1 -
                (1 / (this->initialConfig.springsContourLength[springIndex])));
          }
          // detect the position in the spring
          std::vector<double> partitionsStrand;
          partitionsStrand.reserve(springParticipants.size() - 1);
          for (size_t i = 0; i < springParticipants.size() - 1; ++i) {
            partitionsStrand.push_back(
              this->currentSpringPartitionsVec[this->initialConfig
                                                 .localToGlobalSpringIndex.at(
                                                   springIndex)[i]]);
          }

          bool wasAdded = false;
          size_t targetIndexInSpring = 0;
          double cumulativePartition = 0.0;
          for (size_t p_idx = 0; p_idx < partitionsStrand.size(); ++p_idx) {
            cumulativePartition += partitionsStrand[p_idx];
            if (cumulativePartition > alpha) {
              targetIndexInSpring = p_idx;
              if (p_idx > 0) {
                alpha = alpha - (cumulativePartition - partitionsStrand[p_idx]);
              }
              wasAdded = true;
              break;
            }
          }
          if (!wasAdded) {
            targetIndexInSpring = springParticipants.size() - 2;
            if (partitionsStrand.size() > 0) {
              alpha = alpha - (cumulativePartition -
                               partitionsStrand[partitionsStrand.size() - 1]);
            }
          }

          RUNTIME_EXP_IFN(APPROX_WITHIN(alpha, 0.0, 1.0, 1e-12),
                          "alpha must be between 0 and 1, got " +
                            std::to_string(alpha) + ".");

          // have to adjust the existing springs, too!
          size_t springPartner1 = springParticipants[targetIndexInSpring];
          size_t springPartner2 = springParticipants[targetIndexInSpring + 1];
          size_t newNodeIdx = currentNrOfLinks + i;

          // update connectivity
          size_t lastSpringIndex =
            this->initialConfig.localToGlobalSpringIndex.at(
              springIndex)[targetIndexInSpring];
          size_t newSpringIndex =
            currentNrOfPartialSprings + partialSpringsAdded;

          this->initialConfig.localToGlobalSpringIndex.at(springIndex)
            .insert(this->initialConfig.localToGlobalSpringIndex.at(springIndex)
                        .begin() +
                      targetIndexInSpring + 1,
                    newSpringIndex);
          this->initialConfig.partialToFullSpringIndex.emplace(newSpringIndex,
                                                               springIndex);

          // adjust also the coordinates
          this->currentDisplacements.segment(3 * newNodeIdx, 3) =
            Eigen::Vector3d::Zero();
          if (this->initialConfig.springPartIndexA[lastSpringIndex] ==
              springPartner1) {
            this->initialConfig.springPartIndexB[lastSpringIndex] = newNodeIdx;
            for (size_t offset = 0; offset < 3; ++offset) {
              this->initialConfig
                .springPartCoordinateIndexB[3 * lastSpringIndex + offset] =
                3 * newNodeIdx + offset;
            }
          } else {
            assert(this->initialConfig.springPartIndexA[lastSpringIndex] ==
                   springPartner2);
            this->initialConfig.springPartIndexA[lastSpringIndex] = newNodeIdx;
            for (size_t offset = 0; offset < 3; ++offset) {
              this->initialConfig
                .springPartCoordinateIndexA[3 * lastSpringIndex + offset] =
                3 * newNodeIdx + offset;
            }
          }
          // add the new one
          this->initialConfig.springPartIndexA[newSpringIndex] = newNodeIdx;
          this->initialConfig.springPartIndexB[newSpringIndex] = springPartner2;
          for (size_t offset = 0; offset < 3; ++offset) {
            this->initialConfig
              .springPartCoordinateIndexA[3 * newSpringIndex + offset] =
              3 * newNodeIdx + offset;
            this->initialConfig
              .springPartCoordinateIndexB[3 * newSpringIndex + offset] =
              3 * springPartner2 + offset;
          }

          this->currentSpringPartitionsVec[newSpringIndex] =
            this->currentSpringPartitionsVec[lastSpringIndex] - alpha;
          RUNTIME_EXP_IFN(
            APPROX_WITHIN(this->currentSpringPartitionsVec[newSpringIndex],
                          0.0,
                          1.0,
                          1e-12),
            "Spring partition must be between 0 and 1, got " +
              std::to_string(this->currentSpringPartitionsVec[newSpringIndex]) +
              ".");
          this->currentSpringPartitionsVec[lastSpringIndex] = alpha;
          RUNTIME_EXP_IFN(
            APPROX_WITHIN(this->currentSpringPartitionsVec[lastSpringIndex],
                          0.0,
                          1.0,
                          1e-12),
            "Spring partition must be between 0 and 1, got " +
              std::to_string(this->currentSpringPartitionsVec[newSpringIndex]) +
              ".");

          this->initialConfig.linkIndicesOfSprings[springIndex].insert(
            this->initialConfig.linkIndicesOfSprings[springIndex].begin() +
              targetIndexInSpring +
              1, // + 1 to compensate for the first cross-link
            newNodeIdx);

          partialSpringsAdded += 1;
          springIndexIndex += 1;
        }
      }
      this->initialConfig.nrOfPartialSprings += partialSpringsAdded;

      size_t nrOfPartitionedSprings = 0;
      for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
        if (this->initialConfig.linkIndicesOfSprings[i].size() > 2) {
          nrOfPartitionedSprings += 1;
        }
      }
      this->initialConfig.nrOfSpringsWithPartition = nrOfPartitionedSprings;

      // do we really want to?
      this->validateNetwork();
      assert(partialSpringsAdded == 2 * additionalLen);
    };

    /**
     * @brief Get the Average Spring Length at the current step
     *
     * @return double
     */
    double MEHPForceBalance::getAverageSpringLength() const
    {
      double r2 = 0.0;
      for (int i = 0; i < this->initialConfig.nrOfSprings; i++) {
        double r2local = 0.0;
        for (int j = 0; j < 3; ++j) {
          r2local += this->currentSpringDistances[i * 3 + j] *
                     this->currentSpringDistances[i * 3 + j];
        }
        r2 += sqrt(r2local);
      }
      return r2 / this->initialConfig.nrOfSprings;
    }

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3>
    MEHPForceBalance::evaluateStressTensorLinkBased(
      const ForceBalanceNetwork &net,
      const Eigen::VectorXd& u,
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
                                            u,
                                            springPartitions,
                                            debugNrSpringsVisited,
                                            kappa0,
                                            oneOverSpringPartitionUpperLimit)
            : this->evaluateForceOnCrossLink(linkIdx,
                                             net,
                                             u,
                                             springPartitions,
                                             debugNrSpringsVisited,
                                             kappa0,
                                             oneOverSpringPartitionUpperLimit);
        /* spring contribution to the overall stress tensor */
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
    std::array<std::array<double, 3>, 3> MEHPForceBalance::evaluateStressTensor(
      const ForceBalanceNetwork &net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {
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

      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd relevantPartialDistancesA =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA));
      this->handlePBC(net, relevantPartialDistancesA);

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
          net.partialToFullSpringIndex.at(partialSpringIdx);
          const double N = net.springsContourLength[totalSpringIndex];
        double denominator = 1 / (springPartitions[partialSpringIdx] * N);
        if (oneOverSpringPartitionUpperLimit > 0. ||
            !std::isfinite(denominator)) {
          denominator = CLAMP_ONE_OVER_SPRINGPARTITION(
            denominator, N, oneOverSpringPartitionUpperLimit);
        }
        /* spring contribution to the overall stress tensor */
        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            double contribution =
              distance[j] * distance[k] * kappa0 * denominator;
            // if (std::isfinite(denominator) && std::isfinite(contribution)) {
            stress[j][k] += contribution;
            // }
          }
        }
      }

      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          stress[i][j] *= oneOverVolume;
        }
      }

      return stress;
    }

    std::array<std::array<double, 3>, 3> MEHPForceBalance::getStressTensor(
      const double oneOverSpringPartitionUpperLimit) const
    {
      return this->evaluateStressTensor(this->initialConfig,
                                        this->currentDisplacements,
                                        this->currentSpringPartitionsVec,
                                        1.0,
                                        oneOverSpringPartitionUpperLimit);
    }

    std::array<std::array<double, 3>, 3>
    MEHPForceBalance::getStressTensorLinkBased(
      const double oneOverSpringPartitionUpperLimit,
      const bool xlinksOnly) const
    {
      return this->evaluateStressTensorLinkBased(
        this->initialConfig,
        this->currentDisplacements,
        this->currentSpringPartitionsVec,
        1.0,
        oneOverSpringPartitionUpperLimit,
        xlinksOnly);
    }

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
    MEHPForceBalance::getEffectiveFunctionalityOfAtoms(double tolerance) const
    {
      std::unordered_map<long int, int> results;
      results.reserve(this->initialConfig.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        this->getNrOfActiveSpringsConnected(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
        results.emplace(this->initialConfig.oldAtomIds[i],
                        nrOfActiveSpringsConnected[i]);
      }
      return results;
    }

    /**
     * @brief Get the Ids Of active Nodes
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @param minimumNrOfActiveConnections the number of active springs
     * required for this node to qualify as active
     * @return std::vector<long int> the atom ids
     */
    std::vector<long int> MEHPForceBalance::getIdsOfActiveNodes(
      double tolerance,
      int minimumNrOfActiveConnections,
      int maximumNrOfActiveConnections,
      bool usePartial) const
    {
      std::vector<long int> results;
      results.reserve(this->initialConfig.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        usePartial ? this->getNrOfActiveSpringsConnected(tolerance)
                   : this->getNrOfActivePartialSpringsConnected(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
        if (nrOfActiveSpringsConnected[i] >= minimumNrOfActiveConnections &&
            (maximumNrOfActiveConnections < 0 ||
             maximumNrOfActiveConnections >= nrOfActiveSpringsConnected[i])) {
          results.push_back(this->initialConfig.oldAtomIds[i]);
        }
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
    Eigen::VectorXi MEHPForceBalance::getNrOfActiveSpringsConnected(
      double tolerance) const
    {
      Eigen::VectorXi nrOfActiveSpringsConnected =
        Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
      ArrayXb springIsActive =
        this->findActiveSprings(this->currentSpringDistances, tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfSprings; i++) {
        if (springIsActive[i] == true) /* active spring */
        {
          int a = this->initialConfig.springIndexA[i];
          int b = this->initialConfig.springIndexB[i];
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
    Eigen::VectorXi MEHPForceBalance::getNrOfActivePartialSpringsConnected(
      double tolerance) const
    {
      Eigen::VectorXi nrOfActivePartialSpringsConnected =
        Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
      ArrayXb springIsActive =
        this->findActiveSprings(this->currentPartialSpringDistances, tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfPartialSprings; i++) {
        if (springIsActive[i] == true) /* active spring */
        {
          int a = this->initialConfig.springPartIndexA[i];
          int b = this->initialConfig.springPartIndexB[i];
          if (!this->initialConfig.linkIsSliplink[a]) {
            ++(nrOfActivePartialSpringsConnected[a]);
          }

          if (!this->initialConfig.linkIsSliplink[b]) {
            ++(nrOfActivePartialSpringsConnected[b]);
          }
        }
      }
      return nrOfActivePartialSpringsConnected;
    }

    /**
     * @brief Get the Gamma Factor at the current step
     *
     * @param r02 the melt <R_0^2>, for phantom = Nb^2
     * @param nrOfChains the nr of chains to average over (can be different
     * from the nr of springs thanks to omitted free chains or primary loops)
     * @return double
     */
    double MEHPForceBalance::getGammaFactor(double r02, int nrOfChains) const
    {
      if (r02 < 0) {
        r02 = this->defaultR0Squared;
      }
      if (nrOfChains < 1) {
        nrOfChains = this->defaultNrOfChains;
      }

      return this->evaluateGammaFactor(
        this->currentSpringDistances, r02, nrOfChains);
    }

    /**
     * @brief Convert the universe to a network
     *
     * @param net the target network
     * @param crosslinkerType the atom type of the crosslinker
     * @return true
     * @return false
     */
    bool MEHPForceBalance::ConvertNetwork(ForceBalanceNetwork &net,
                                          const int crosslinkerType,
                                          bool remove2functionalCrosslinkers)
    {
      std::vector<pylimer_tools::entities::Atom> xlinkers =
        this->universe.getAtomsOfType(crosslinkerType);

      if (remove2functionalCrosslinkers) {
        for (pylimer_tools::entities::Atom xlinker : xlinkers) {
          // change type of cross-linkers with a degree <= 2 to "normal",
          // non-cross-link beads
          size_t vertexId = this->universe.getIdxByAtomId(xlinker.getId());
          if (this->universe.computeFunctionalityForVertex(vertexId) <= 2) {
            this->universe.setPropertyValue(
              vertexId, "type", crosslinkerType - 1);
          }
        }
        xlinkers = this->universe.getAtomsOfType(crosslinkerType);
      }

      size_t nrOfXlinks = xlinkers.size();

      std::vector<pylimer_tools::entities::Molecule> crosslinkerChains =
        this->universe.getChainsWithCrosslinker(crosslinkerType);

      // need to include all but dangling and free chains in order to
      // model entanglement
      size_t nrOfSprings = 0;
      for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
        std::vector<pylimer_tools::entities::Atom> endAtoms =
          crosslinkerChains[i].getAtomsOfType(crosslinkerType);
        if (crosslinkerChains[i].getType() ==
              pylimer_tools::entities::MoleculeType::NETWORK_STRAND ||
            crosslinkerChains[i].getType() ==
              pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
          nrOfSprings += 1;
        }
      }

      // crosslinkerUniverse.simplify();
      pylimer_tools::entities::Box box = this->universe.getBox();
      net.L[0] = box.getLx();
      net.L[1] = box.getLy();
      net.L[2] = box.getLz();
      net.boxHalfs[0] = 0.5 * net.L[0];
      net.boxHalfs[1] = 0.5 * net.L[1];
      net.boxHalfs[2] = 0.5 * net.L[2];
      net.nrOfNodes = nrOfXlinks;
      net.nrOfLinks = nrOfXlinks;
      net.nrOfSprings = nrOfSprings;
      net.nrOfPartialSprings = nrOfSprings;
      net.nrOfSpringsWithPartition = 0;
      net.coordinates = Eigen::VectorXd::Zero(3 * net.nrOfLinks);
      net.oldAtomIds = Eigen::ArrayXi::Zero(net.nrOfLinks);
      net.linkIsSliplink = ArrayXb::Constant(net.nrOfLinks, false);
      net.springIndicesOfLinks.reserve(net.nrOfLinks);
      net.partialToFullSpringIndex.reserve(net.nrOfLinks);
      for (size_t i = 0; i < net.nrOfLinks; ++i) {
        net.springIndicesOfLinks.push_back(std::vector<size_t>());
      }
      net.linkIndicesOfSprings.reserve(net.nrOfSprings);
      this->currentSpringPartitionsVec =
        Eigen::VectorXd::Ones(net.nrOfSprings);
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        net.linkIndicesOfSprings.push_back(std::vector<size_t>());
      }
      net.springIndexA = Eigen::ArrayXi::Zero(net.nrOfSprings);
      net.springIndexB = Eigen::ArrayXi::Zero(net.nrOfSprings);
      net.springCoordinateIndexA = Eigen::ArrayXi::Zero(3 * net.nrOfSprings);
      net.springCoordinateIndexB = Eigen::ArrayXi::Zero(3 * net.nrOfSprings);
      net.springIsActive = ArrayXb::Constant(net.nrOfSprings, false);
      net.springsContourLength = Eigen::VectorXd::Zero(net.nrOfSprings);
      net.oldAtomIdToSpringIndex.reserve(this->universe.getNrOfAtoms());
      net.springToMoleculeIds.reserve(nrOfSprings);

      // convert beads
      std::map<int, int> atomIdToNode;
      for (size_t i = 0; i < xlinkers.size(); ++i) {
        pylimer_tools::entities::Atom atom = xlinkers[i];
        atomIdToNode[atom.getId()] = i;
        net.oldAtomIds[i] = atom.getId();
        net.coordinates[3 * i + 0] = atom.getX();
        net.coordinates[3 * i + 1] = atom.getY();
        net.coordinates[3 * i + 2] = atom.getZ();
      }

      // convert springs
      size_t spring_idx = 0;
      // net.connectivityToSpringIndex.reserve(nrOfSprings);
      for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
        std::vector<pylimer_tools::entities::Atom> xlinkersOfChain =
          crosslinkerChains[i].getAtomsOfType(crosslinkerType);
        long int nodeIdxFrom;
        long int nodeIdxTo;
        bool addChain = false;
        if (crosslinkerChains[i].getType() ==
            pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
          assert(xlinkersOfChain.size() == 2);
          nodeIdxFrom = atomIdToNode.at(xlinkersOfChain[0].getId());
          nodeIdxTo = atomIdToNode.at(xlinkersOfChain[1].getId());
          if (nodeIdxFrom > nodeIdxTo) {
            std::swap(nodeIdxFrom, nodeIdxTo);
          }
          addChain = true;

          net.springsContourLength[spring_idx] =
            crosslinkerChains[i].getNrOfAtoms() - 1; // TODO: -2?
        } else if (crosslinkerChains[i].getType() ==
                   pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
          assert(xlinkersOfChain.size() == 1 ||
                 (xlinkersOfChain.size() == 2 &&
                  xlinkersOfChain[0].getId() == xlinkersOfChain[1].getId()));

          nodeIdxFrom = atomIdToNode.at(xlinkersOfChain[0].getId());
          nodeIdxTo = nodeIdxFrom;
          addChain = true;

          net.springsContourLength[spring_idx] =
            crosslinkerChains[i].getNrOfAtoms(); // TODO: -1?
        }

        if (addChain) {
          net.springToMoleculeIds.push_back(i);
          std::vector<pylimer_tools::entities::Atom> allChainAtoms =
            crosslinkerChains[i].getAtoms();
          for (pylimer_tools::entities::Atom a : allChainAtoms) {
            net.oldAtomIdToSpringIndex[a.getId()] = spring_idx;
          }

          net.springIndicesOfLinks[nodeIdxFrom].push_back(spring_idx);
          if (nodeIdxFrom != nodeIdxTo) {
            net.springIndicesOfLinks[nodeIdxTo].push_back(spring_idx);
          }

          net.linkIndicesOfSprings[spring_idx].push_back(nodeIdxFrom);
          net.linkIndicesOfSprings[spring_idx].push_back(nodeIdxTo);

          net.springIndexA[spring_idx] = nodeIdxFrom;
          net.springIndexB[spring_idx] = nodeIdxTo;
          for (size_t j = 0; j < 3; j++) {
            net.springCoordinateIndexA[3 * spring_idx + j] =
              nodeIdxFrom * 3 + j;
            net.springCoordinateIndexB[3 * spring_idx + j] = nodeIdxTo * 3 + j;
          }

          std::vector<size_t> zeroMap;
          zeroMap.push_back(spring_idx);
          net.localToGlobalSpringIndex.emplace(spring_idx, zeroMap);
          net.partialToFullSpringIndex.emplace(spring_idx, spring_idx);

          spring_idx += 1;
        }
      }

      net.springPartCoordinateIndexA = net.springCoordinateIndexA;
      net.springPartCoordinateIndexB = net.springCoordinateIndexB;
      net.springPartIndexA = net.springIndexA;
      net.springPartIndexB = net.springIndexB;

      // box volume
      net.vol = net.L[0] * net.L[1] * net.L[2];
      if (net.springsContourLength.size() > 0) {
        net.meanSpringContourLength = net.springsContourLength.mean();
      } else {
        net.meanSpringContourLength = 0.0;
      }

      return spring_idx == net.nrOfSprings;
    };

    bool MEHPForceBalance::validateNetwork(const ForceBalanceNetwork &net)
    {
      RUNTIME_EXP_IFN(!std::isinf(net.L[0]) && !std::isnan(net.L[0]),
                      "Box direction x must be scalar");
      RUNTIME_EXP_IFN(!std::isinf(net.L[1]) && !std::isnan(net.L[1]),
                      "Box direction y must be scalar");
      RUNTIME_EXP_IFN(!std::isinf(net.L[2]) && !std::isnan(net.L[2]),
                      "Box direction z must be scalar");
      RUNTIME_EXP_IFN(net.coordinates.size() == net.nrOfLinks * 3,
                      "Invalid size of coordinates");
      RUNTIME_EXP_IFN(this->currentDisplacements.size() == net.nrOfLinks * 3,
                      "Invalid size of current displacement");
      RUNTIME_EXP_IFN(net.localToGlobalSpringIndex.size() == net.nrOfSprings,
                      "Invalid size of connectivity map");
      RUNTIME_EXP_IFN(net.springsContourLength.size() == net.nrOfSprings,
                      "Invalid size of contour lengths");
      RUNTIME_EXP_IFN(net.springIndicesOfLinks.size() == net.nrOfLinks,
                      "Invalid size of spring indices of links");
      RUNTIME_EXP_IFN(net.linkIndicesOfSprings.size() == net.nrOfSprings,
                      "Invalid size of link indices of springs");
      RUNTIME_EXP_IFN(net.linkIsSliplink.size() == net.nrOfLinks,
                      "Invalid size of link is sliplink");
      RUNTIME_EXP_IFN(
        net.linkIsSliplink.count() == net.nrOfLinks - net.nrOfNodes,
        "Nr of nodes plus nr of slp-links should give the total nr of links");
      RUNTIME_EXP_IFN(net.oldAtomIds.size() == net.nrOfNodes,
                      "Invalid size of old atom ids");
      RUNTIME_EXP_IFN(net.springCoordinateIndexA.size() ==
                        net.nrOfSprings * 3,
                      "Invalid size of springCoordinateIndexA");
      RUNTIME_EXP_IFN(net.springCoordinateIndexB.size() ==
                        net.nrOfSprings * 3,
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
      RUNTIME_EXP_IFN(net.springPartIndexB.size() == net.nrOfPartialSprings,
                      "Invalid size of springPartIndexB");
      RUNTIME_EXP_IFN(net.springIsActive.size() == net.nrOfSprings,
                      "Invalid size of springIsActive");
      RUNTIME_EXP_IFN(this->currentSpringPartitionsVec.size() ==
                        net.nrOfPartialSprings,
                      "Invalid size of currentSpringPartitionsVec");
      RUNTIME_EXP_IFN(
        APPROX_EQUAL(
          this->currentSpringPartitionsVec.sum(), net.nrOfSprings, 1e-3),
        "Spring partitions should sum to 1 per spring, got " +
          std::to_string(this->currentSpringPartitionsVec.sum()) + " for " +
          std::to_string(net.nrOfSprings) + " springs.");
      RUNTIME_EXP_IFN(
        net.partialToFullSpringIndex.size() == net.nrOfPartialSprings,
        "Every partial spring must be able to map to the full spring.");
      for (size_t i = 0; i < this->currentSpringPartitionsVec.size(); i++) {
        RUNTIME_EXP_IFN(
          APPROX_WITHIN(this->currentSpringPartitionsVec[i], 0.0, 1.0, 1e-12),
          "Spring partitions must be between 0 & 1, got " +
            std::to_string(this->currentSpringPartitionsVec[i]) +
            " at i = " + std::to_string(i) + ".");
      }
      for (size_t link_idx = 0; link_idx < net.nrOfLinks; ++link_idx) {
        std::vector<size_t> thisLinksSprings =
          net.springIndicesOfLinks[link_idx];
        for (size_t spring_idx : thisLinksSprings) {
          std::vector<size_t> thisSpringsLinks =
            net.linkIndicesOfSprings[spring_idx];
          RUNTIME_EXP_IFN(std::find(thisSpringsLinks.begin(),
                                    thisSpringsLinks.end(),
                                    link_idx) != thisSpringsLinks.end(),
                          "Spring must have a connection to the link, too.");
        }
      }
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        RUNTIME_EXP_IFN(net.linkIndicesOfSprings[i].size() >= 2,
                        "Each spring requires at least two links, got " +
                          std::to_string(net.linkIndicesOfSprings[i].size()) +
                          " at i = " + std::to_string(i) + ".");
        RUNTIME_EXP_IFN(net.localToGlobalSpringIndex.at(i).size() ==
                          net.linkIndicesOfSprings[i].size() - 1,
                        "Require a global index for each local one");
        RUNTIME_EXP_IFN(
          net.linkIndicesOfSprings[i][0] <=
            net.linkIndicesOfSprings[i]
                                     [net.linkIndicesOfSprings[i].size() - 1],
          "Springs must have increasing end-point indices");
        std::vector<size_t> links = net.linkIndicesOfSprings[i];
        for (size_t link_idx : links) {
          std::vector<size_t> thisLinksSprings =
            net.springIndicesOfLinks[link_idx];
          RUNTIME_EXP_IFN(std::find(thisLinksSprings.begin(),
                                    thisLinksSprings.end(),
                                    i) != thisLinksSprings.end(),
                          "Link must have a connection to the spring, too.");
        }
        // also check the sum of the partials
        std::vector<size_t> globalSpringIndices =
          net.localToGlobalSpringIndex.at(i);
        double sum = 0.0;
        for (size_t globalIdx : globalSpringIndices) {
          sum += this->currentSpringPartitionsVec[globalIdx];
        }
        RUNTIME_EXP_IFN(
          APPROX_EQUAL(sum, 1.0, 1e-10),
          "Spring partitions of one spring must sum to one, got " +
            std::to_string(sum) + ".");
      }
      for (size_t i = 0; i < net.nrOfPartialSprings; i++) {
        size_t fullIdx = net.partialToFullSpringIndex.at(i);
        size_t partialEndA = net.springPartIndexA[i];
        size_t partialEndB = net.springPartIndexB[i];
        if (!net.linkIsSliplink[partialEndA]) {
          RUNTIME_EXP_IFN(net.springIndexA[fullIdx] == partialEndA ||
                            net.springIndexB[fullIdx] == partialEndA,
                          "Expect mapping of springs to work");
        }
        if (!net.linkIsSliplink[partialEndB]) {
          RUNTIME_EXP_IFN(net.springIndexA[fullIdx] == partialEndB ||
                            net.springIndexB[fullIdx] == partialEndB,
                          "Expect mapping of springs to work");
        }
        for (size_t dir = 0; dir < 3; ++dir) {
          RUNTIME_EXP_IFN(net.springPartCoordinateIndexA[3 * i + dir] ==
                            3 * partialEndA + dir,
                          "Spring part index and coordinate index must match.");
          RUNTIME_EXP_IFN(net.springPartCoordinateIndexB[3 * i + dir] ==
                            3 * partialEndB + dir,
                          "Spring part index and coordinate index must match.");
        }
      }
      return true;
    }
  }
}
}
