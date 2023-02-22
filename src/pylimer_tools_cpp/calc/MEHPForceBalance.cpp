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
  isPartialSpring, val, N, oneOverSpringPartitionUpperLimit)                   \
  (!isPartialSpring)                                                           \
    ? val                                                                      \
    : std::clamp(val,                                                          \
                 (oneOverSpringPartitionUpperLimit > 0.)                       \
                   ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)          \
                   : 0.0,                                                      \
                 (oneOverSpringPartitionUpperLimit > 0.)                       \
                   ? oneOverSpringPartitionUpperLimit                          \
                   : N);

    /**
     * FORCE RELAXATION
     */
    void MEHPForceBalance::runForceRelaxation(
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

      ForceBalanceNetwork net = this->initialConfig;
      double removalTolerance =
        (inactiveRemovalCutoff > 0.0)
          ? inactiveRemovalCutoff
          : (0.25 * std::pow(net.vol / this->universe.getNrOfAtoms(), 1. / 3.));

      /* array allocation */
      Eigen::VectorXd u = this->currentDisplacements;
      Eigen::VectorXd springPartitions = this->currentSpringPartitionsVec;
      std::vector<Eigen::ArrayXi> independentVertexSets;
      double maxDistanceMoved = 0.0;
      size_t indexOfMaxDistanceMoved = 0;
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
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(
          net, springPartitions, oneOverSpringPartitionUpperLimit);
      const double initialResidual = (initialResidualToUse > 0.)
                                       ? initialResidualToUse
                                       : this->getDisplacementResidualNormFor(
                                           net, u, oneOverSpringPartitions);
      const double minN = net.springsContourLength.minCoeff();
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
      // actual loop
      do {
        if (allowSlipLinksToPassEachOther != LinkSwappingMode::NO_SWAPPING) {
          if (swappingFrequency > 0 &&
              (iterationsDone % swappingFrequency) == 0) {
            if (allowSlipLinksToPassEachOther ==
                LinkSwappingMode::SLIPLINKS_ONLY) {
              this->swapSlipLinks(
                net, u, springPartitions, oneOverSpringPartitionUpperLimit);
            } else if (allowSlipLinksToPassEachOther == LinkSwappingMode::ALL) {
              this->swapSlipLinksInclXlinks(
                net, u, springPartitions, oneOverSpringPartitionUpperLimit);
            } else if (allowSlipLinksToPassEachOther ==
                       LinkSwappingMode::ALL_MC) {
              this->moveSlipLinksToTheirBestBranch(
                net,
                u,
                springPartitions,
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                false);
            } else if (allowSlipLinksToPassEachOther ==
                       LinkSwappingMode::ALL_MC_CYCLE) {
              this->moveSlipLinksToTheirBestBranch(
                net,
                u,
                springPartitions,
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                true);
            } else {
              throw std::invalid_argument(
                "This swapping mode is currently not supported.");
            }
            oneOverSpringPartitions = this->assembleOneOverSpringPartition(
              net, springPartitions, oneOverSpringPartitionUpperLimit);
          }
        }
        maxDistanceMoved = 0.0;
        currentResidual = 0.0;

        // place slip-link
        for (size_t link_idx = net.nrOfNodes; link_idx < net.nrOfLinks;
             ++link_idx) {
          // std::cout << "Handling " << link_idx << " of " << net.nrOfNodes
          //           << " / " << net.nrOfLinks << std::endl;

          // std::cout << "Still handling " << link_idx << " of " <<
          // net.nrOfNodes
          //           << " / " << net.nrOfLinks << std::endl;
          int innerIterationsDone = 0;
          do {
            double r2 =
              this->updateSpringPartition(net,
                                          u,
                                          springPartitions,
                                          oneOverSpringPartitions,
                                          link_idx,
                                          oneOverSpringPartitionUpperLimit,
                                          allowSlipLinksToPassEachOther);
            double displacementDone =
              this->displaceToMeanPosition(net,
                                           u,
                                           springPartitions,
                                           link_idx,
                                           oneOverSpringPartitionUpperLimit);
            innerIterationsDone += 1;
          } while (doInnerIterations && innerIterationsDone < 50);
        }

        intermediateResidual =
          this->getDisplacementResidualNormFor(net, u, oneOverSpringPartitions);

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

        currentResidual =
          this->getDisplacementResidualNormFor(net, u, oneOverSpringPartitions);
        iterationsDone += 1;
        if (iterationsDone % 10 == 0) {
          if (simplificationMode ==
                StructureSimplificationMode::INACTIVE_ONLY ||
              simplificationMode == StructureSimplificationMode::ALL_TIM) {
            // std::cout << "Removing inactive cross-links" << std::endl;
            // default tolerance: 0.25*atom's cube length
            size_t nRemoved = this->removeInactiveCrosslinks(
              net, u, springPartitions, removalTolerance);
            net.meanSpringContourLength = net.springsContourLength.mean();
            if (nRemoved > 0) {
              std::cout << "Removed " << nRemoved << " inactive springs. "
                        << std::endl;
            }
            // this->validateNetwork(net, u, springPartitions);
          }
          if (simplificationMode == StructureSimplificationMode::X2F_ONLY ||
              simplificationMode == StructureSimplificationMode::ALL_TIM) {
            // std::cout << "Removing 2-f cross-links" << std::endl;
            size_t nRemoved =
              this->removeTwofunctionalCrosslinks(net, u, springPartitions);
            net.meanSpringContourLength = net.springsContourLength.mean();
            if (nRemoved > 0) {
              std::cout << "Removed " << nRemoved
                        << " cross-linkers with f = 2. " << std::endl;
            }
            // this->validateNetwork(net, u, springPartitions);
          }
          if (simplificationMode == StructureSimplificationMode::ALL_ANDREI) {
            std::cout << "Removing cross-links and springs, Andrei's way"
                      << std::endl;
            this->doRemovalAndreisWay(
              net, u, springPartitions, removalTolerance);
          }
          if (simplificationMode !=
              StructureSimplificationMode::NO_SIMPLIFICATION) {
            this->cleanupPrimaryLoopsInStructure(net);
            this->validateNetwork(net, u, springPartitions);
          }
        }
        if (outputFrequency > 0 && iterationsDone % outputFrequency == 0) {
          std::cout << "Iteration " << iterationsDone << " " << maxDistanceMoved
                    << " by " << indexOfMaxDistanceMoved
                    << ". Residual: " << currentResidual * minN * minN << " ("
                    << (currentResidual / initialResidual) << ") "
                    << " from: " << initialResidual * minN * minN << " via "
                    << intermediateResidual * minN * minN << " ("
                    << (intermediateResidual / initialResidual) << ") "
                    << "\n";
          std::array<std::array<double, 3>, 3> stressTensor =
            this->evaluateStressTensor(
              net, u, springPartitions, 1.0, oneOverSpringPartitionUpperLimit);
          std::cout << "To stress tensor diagonal: " << stressTensor[0][0]
                    << ", " << stressTensor[1][1] << ", " << stressTensor[2][2]
                    << std::endl;
        }
      } while (currentResidual / initialResidual > xtol &&
               iterationsDone < maxNrOfSteps);

      // query solution & exit reason
      this->exitReason = (iterationsDone == maxNrOfSteps)
                           ? ExitReason::MAX_STEPS
                           : ExitReason::X_TOLERANCE;
      this->nrOfStepsDone += iterationsDone;
      std::cout << iterationsDone << " steps done. "
                << "Last max distance moved: " << maxDistanceMoved << std::endl;

      assert(u.size() == 3 * net.nrOfLinks);
      this->initialConfig = net;
      this->currentDisplacements = u;
      this->currentSpringPartitionsVec = springPartitions;
      this->validateNetwork();
      this->currentSpringDistances = this->evaluateSpringDistances(
        net, this->currentDisplacements, this->is2D);
      this->currentPartialSpringDistances =
        this->evaluatePartialSpringDistances(
          net, this->currentDisplacements, is2D);
    }

    double MEHPForceBalance::displaceLinksToMeanPosition(
      const ForceBalanceNetwork& net,
      Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions0,
      double damping) const
    {

      Eigen::ArrayXi mask =
        Eigen::ArrayXi::LinSpaced(3 * net.nrOfLinks, 0, 3 * net.nrOfLinks - 1);
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
      const ForceBalanceNetwork& net,
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
      const ForceBalanceNetwork& net,
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
      this->box.handlePBC(relevantPartialDistancesA);

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
      u(resultingCoordinateIndexMask) += damping * objectivesToSet;

      return maxDiff;
    };

    double MEHPForceBalance::getDisplacementResidualNorm(
      double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(this->initialConfig,
                                             this->currentSpringPartitionsVec,
                                             oneOverSpringPartitionUpperLimit);
      Eigen::VectorXd displacements = this->currentDisplacements;
      return this->getDisplacementResidualNormFor(
        this->initialConfig, displacements, oneOverSpringPartitions);
    }

    double MEHPForceBalance::getDisplacementResidualNormFor(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& oneOverSpringPartitions) const
    {
      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd relevantPartialDistancesA =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA));
      for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
        for (size_t dir = 0; dir < 3; ++dir) {
          if (std::abs(relevantPartialDistancesA[3 * i + dir]) >
              50. * net.L[dir]) {
            std::cerr
              << "WARNING: Spring " << i << " between "
              << net.springPartIndexA[i] << " and " << net.springPartIndexB[i]
              << " has a length of " << relevantPartialDistancesA[3 * i + dir]
              << " in dir " << dir << " from "
              << displacedCoords[net.springPartCoordinateIndexB[3 * i + dir]]
              << " minus "
              << displacedCoords[net.springPartCoordinateIndexA[3 * i + dir]]
              << " meaning it will probably fail in PBC." << std::endl;
          }
        }
      }
      this->box.handlePBC(relevantPartialDistancesA);
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
      const ForceBalanceNetwork& net) const
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
      const ForceBalanceNetwork& net) const
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
      const ForceBalanceNetwork& net) const
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
     * @brief Assemble all indices of partial springs for a particular slip-link
     *
     * @param linkIdx
     * @return void
     */
    void MEHPForceBalance::setSpringpartitionIndicesOfSliplink(
      std::vector<size_t>& results,
      const ForceBalanceNetwork& net,
      const size_t linkIdx) const
    {
      INVALIDARG_EXP_IFN(
        linkIdx < net.nrOfLinks,
        "Cannot set spring partition of index higher than nr. of links.");
      INVALIDARG_EXP_IFN(net.linkIsSliplink[linkIdx], "Link must be slip-link");
      std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
      size_t indexIndex = 0;
      while (results.size() < 4) {
        results.push_back(0);
      }
      for (size_t springIndex : springIndices) {
        std::vector<size_t> springsPartners =
          net.linkIndicesOfSprings[springIndex];
        for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          if (springsPartners[partner_idx] == linkIdx) {
            RUNTIME_EXP_IFN(
              indexIndex < 4,
              "Expect spring partitions indices of link not to exceed 4.");
            size_t currentSpringGlobalIdx =
              net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
            size_t neighbourSpringGlobalIdx =
              net.localToGlobalSpringIndex[springIndex][partner_idx];
            results[indexIndex] = currentSpringGlobalIdx;
            indexIndex++;
            results[indexIndex] = neighbourSpringGlobalIdx;
            indexIndex++;
          }
        }
      }
      assert(indexIndex == 4);
    }

    /**
     * @brief Remove double listed springs from cross-links
     *
     * @param net
     */
    void MEHPForceBalance::cleanupPrimaryLoopsInStructure(
      ForceBalanceNetwork& net)
    {
      for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
        // remove duplicate mentions of the same spring index
        std::sort(net.springIndicesOfLinks[linkIdx].begin(),
                  net.springIndicesOfLinks[linkIdx].end());
        auto last = std::unique(net.springIndicesOfLinks[linkIdx].begin(),
                                net.springIndicesOfLinks[linkIdx].end());
        net.springIndicesOfLinks[linkIdx].erase(
          last, net.springIndicesOfLinks[linkIdx].end());
      }
    }

    /**
     * @brief Remove cross-links which do not have any springs with a certain
     * minimum length
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param tolerance
     */
    size_t MEHPForceBalance::removeInactiveCrosslinks(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions,
      double tolerance) const
    {
      size_t numRemoved = 0;
      // this->validateNetwork(net, displacements, springPartitions);
      // first, we remove all inactive springs
      for (long int springIdx = net.nrOfSprings - 1; springIdx >= 0;
           --springIdx) {
        std::vector<size_t> involvedPartialSprings =
          net.localToGlobalSpringIndex[springIdx];
        bool isActive = false;
        for (size_t partialSpringIdx : involvedPartialSprings) {
          RUNTIME_EXP_IFN(
            net.coordinates.size() == displacements.size(),
            "Expected coordinates and displacements to have same size, got " +
              std::to_string(net.coordinates.size()) + " and " +
              std::to_string(displacements.size()) + ".");
          // assert(net.coordinates.size() == displacements.size());
          Eigen::Vector3d distance =
            (net.coordinates.segment(3 * net.springPartIndexA[partialSpringIdx],
                                     3) +
             displacements.segment(3 * net.springPartIndexA[partialSpringIdx],
                                   3)) -
            (net.coordinates.segment(3 * net.springPartIndexB[partialSpringIdx],
                                     3) +
             displacements.segment(3 * net.springPartIndexB[partialSpringIdx],
                                   3));
          this->box.handlePBC(distance);
          if (distance.squaredNorm() > tolerance) {
            isActive = true;
            break;
          }
        }
        if (!isActive) {
          // remove this spring
          std::cout << "Removing spring " << springIdx << std::endl;
          this->validateNetwork(net, displacements, springPartitions);
          this->removeSpring(net, displacements, springPartitions, springIdx);
          this->validateNetwork(net, displacements, springPartitions);
          numRemoved += 1;
        }
      }

      // then, we remove all cross-links that are 0- or 1-functional
      for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
           --crosslinkIdx) {
        if (net.springIndicesOfLinks[crosslinkIdx].size() == 0 // f = 0
        ) {
          // std::cout << "Removing x-link " << crosslinkIdx << std::endl;
          this->removeLink(net, displacements, crosslinkIdx);
          this->validateNetwork(net, displacements, springPartitions);
        }

        if ( // or f = 1, NOT primary loop
          net.springIndicesOfLinks[crosslinkIdx].size() == 1 &&
          XOR(
            net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx][0]]
                                    [0] == crosslinkIdx,
            pylimer_tools::utils::last(
              net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx]
                                                               [0]]) ==
              crosslinkIdx)) {
          // need to first remove the spring
          this->removeSpring(net,
                             displacements,
                             springPartitions,
                             net.springIndicesOfLinks[crosslinkIdx][0]);
          numRemoved += 1;
          // to then remove the link
          this->removeLink(net, displacements, crosslinkIdx);
          this->validateNetwork(net, displacements, springPartitions);
        }
      }

      this->validateNetwork(net, displacements, springPartitions);

      return numRemoved;
    }

    /**
     * @brief Add slip-links to this system
     *
     * @param sliplinkDensity
     * @param cutoff
     */
    size_t MEHPForceBalance::randomlyAddSliplinks(
      const size_t nrOfSliplinksToSample,
      const double cutoff,
      const size_t minimumNrOfSliplinks,
      const double sameStrandCutoff,
      const bool excludeCrosslinks,
      const int seed)
    {
      INVALIDARG_EXP_IFN(nrOfSliplinksToSample > minimumNrOfSliplinks,
                         "Maximum nr. should be larger than minimum, got " +
                           std::to_string(nrOfSliplinksToSample) + " and " +
                           std::to_string(minimumNrOfSliplinks) + ".");
      INVALIDARG_EXP_IFN(cutoff > 0.0,
                         "Expected a cutoff > 0.0, got " +
                           std::to_string(cutoff) + ".");
      // RUNTIME_EXP_IFN(this->initialConfig.nrOfLinks ==
      //                   this->initialConfig.nrOfNodes,
      //                 "Slip-links are only added randomly when no other "
      //                 "slip-links are in place yet.");
      INVALIDARG_EXP_IFN(minimumNrOfSliplinks <
                           this->universe.getNrOfAtoms() / 2,
                         "Minimum number of slip-links must be less than the "
                         "possible number of slip-links to place.");
      INVALIDARG_EXP_IFN(nrOfSliplinksToSample <
                           this->universe.getNrOfAtoms() / 2,
                         "Number of slip-links to place must be less than "
                         "the possible number of slip-links to place.");
      // query all the cross-linker chains we actually use to place
      // slip-links on
      std::vector<pylimer_tools::entities::Molecule> crosslinkerChains =
        this->universe.getChainsWithCrosslinker(crosslinkerType);
      std::vector<size_t> usableSpringIdxs;
      usableSpringIdxs.reserve(crosslinkerChains.size());
      // and also query all the corresponding atoms we use to place slip-links
      // on
      size_t nrOfEligibleAtoms = 0;
      std::vector<pylimer_tools::entities::Atom> eligibleAtoms;
      eligibleAtoms.reserve(this->universe.getNrOfAtoms());
      std::vector<bool> vertexIdxIsEligible =
        pylimer_tools::utils::initializeWithValue<bool>(
          this->universe.getNrOfAtoms(), false);
      std::unordered_map<size_t, size_t> atomToStrand;
      atomToStrand.reserve(this->universe.getNrOfAtoms());
      std::unordered_map<size_t, size_t> atomIdxInStrand;
      atomIdxInStrand.reserve(this->universe.getNrOfAtoms());
      size_t springId = 0;
      for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
        pylimer_tools::entities::Molecule chain = crosslinkerChains[i];
        RUNTIME_EXP_IFN(chain.getType() !=
                          pylimer_tools::entities::MoleculeType::UNDEFINED,
                        "Couldn't determine molecule type.");
        if (chain.getType() ==
              pylimer_tools::entities::MoleculeType::PRIMARY_LOOP ||
            chain.getType() ==
              pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
          assert(i == this->initialConfig.springToMoleculeIds[springId]);
          // TODO: also check that this is not a higher order dangling strand
          usableSpringIdxs.push_back(i);
          nrOfEligibleAtoms +=
            crosslinkerChains[i].getNrOfAtoms() -
            crosslinkerChains[i].getAtomsOfType(this->crosslinkerType).size();
          std::vector<pylimer_tools::entities::Atom> atoms =
            crosslinkerChains[i].getAtomsLinedUp(this->crosslinkerType);
          for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
            pylimer_tools::entities::Atom atom = atoms[atomIdx];
            if (atom.getType() != this->crosslinkerType) {
              eligibleAtoms.push_back(atom);
              vertexIdxIsEligible[this->universe.getIdxByAtomId(atom.getId())] =
                true;
              atomToStrand.emplace(atom.getId(), springId);
              atomIdxInStrand.emplace(atom.getId(), atomIdx);
            }
          }
          springId += 1;
        }
      }
      // build neighbourlist
      std::vector<pylimer_tools::entities::Atom> atomsForNeighbourList =
        this->universe.getAtoms();
      if (excludeCrosslinks) {
        // TODO: check whether it is faster to just only query the other ones
        atomsForNeighbourList.erase(
          std::remove_if(atomsForNeighbourList.begin(),
                         atomsForNeighbourList.end(),
                         [&](pylimer_tools::entities::Atom a) -> bool {
                           return a.getType() == this->crosslinkerType;
                         }));
      }
      pylimer_tools::entities::NeighbourList neighbourList =
        pylimer_tools::entities::NeighbourList(
          atomsForNeighbourList, this->universe.getBox(), cutoff);

      std::random_device rd{};
      std::mt19937 rng = std::mt19937(seed > 0 ? seed : rd());
      // build list of random samples
      // this way is more performant than
      // sampling integers and checking whether they have been sampled already
      std::vector<size_t> toSampleFrom;
      std::vector<bool> isMasked;
      toSampleFrom.reserve(this->universe.getNrOfAtoms());
      isMasked.reserve(this->universe.getNrOfAtoms());
      for (size_t i = 0; i < this->universe.getNrOfAtoms(); ++i) {
        toSampleFrom.push_back(i);
        isMasked[i] = (excludeCrosslinks &&
                       this->universe.getAtomByVertexIdx(i).getType() ==
                         this->crosslinkerType);
      }
      std::shuffle(toSampleFrom.begin(), toSampleFrom.end(), rng);

      // the resulting vectors to fill
      std::vector<double> slipLinkXs;
      slipLinkXs.reserve(nrOfSliplinksToSample);
      std::vector<double> slipLinkYs;
      slipLinkYs.reserve(nrOfSliplinksToSample);
      std::vector<double> slipLinkZs;
      slipLinkZs.reserve(nrOfSliplinksToSample);
      std::vector<size_t> slipLinkStrandA;
      slipLinkStrandA.reserve(nrOfSliplinksToSample);
      std::vector<size_t> slipLinkStrandB;
      slipLinkStrandB.reserve(nrOfSliplinksToSample);
      std::vector<double> slipLinkStrandAlpha;
      slipLinkStrandAlpha.reserve(nrOfSliplinksToSample);
      std::vector<double> slipLinkStrandBeta;
      slipLinkStrandBeta.reserve(nrOfSliplinksToSample);

      pylimer_tools::entities::Box box = this->universe.getBox();

      size_t nrOfSlipLinksPlaced = 0;
      size_t nrOfAttempts = 0;
      size_t sampleIdx = 0;
      // the actual sampling loop
      while (nrOfSlipLinksPlaced < minimumNrOfSliplinks ||
             nrOfAttempts < nrOfSliplinksToSample) {
        // first, randomly sample an atom
        while (isMasked[toSampleFrom[sampleIdx]] &&
               sampleIdx < toSampleFrom.size()) {
          sampleIdx += 1;
        }
        if (sampleIdx >= toSampleFrom.size()) {
          // this is a path that should barely ever be reached
          std::cerr << "Sample index exceeds number of samples." << std::endl;
          break;
        }
        size_t sampledVertexId = toSampleFrom[sampleIdx];
        nrOfAttempts += 1;
        sampleIdx += 1;
        pylimer_tools::entities::Atom a1 =
          this->universe.getAtomByVertexIdx(sampledVertexId);
        isMasked[sampledVertexId] = true;
        // then, find neighbouring atoms (but not from the same strand?!)
        std::vector<pylimer_tools::entities::Atom> neighbours =
          neighbourList.getAtomsCloseTo(a1);
        neighbourList.removeAtom(a1);
        // filter the neighbours to include only those from other strands
        // NOTE: this skews the whole thing a bit
        neighbours.erase(
          std::remove_if(
            neighbours.begin(),
            neighbours.end(),
            [&](pylimer_tools::entities::Atom a) -> bool {
              return (atomToStrand[a.getId()] ==
                        atomToStrand[a1.getId()] // do not use "at", because not
                                                 // all atoms in the neighbours
                                                 // have been assigned a strand
                      && std::abs(static_cast<double>(
                           atomIdxInStrand[a.getId()] -
                           atomIdxInStrand[a1.getId()])) < sameStrandCutoff);
            }),
          neighbours.end());
        if (neighbours.size() == 0) {
          std::cerr << "Not enough neighbours found." << std::endl;
          continue;
        }
        // then, randomly select one of them
        pylimer_tools::entities::Atom a2 = neighbours[0];
        if (neighbours.size() > 1) {
          size_t randomA2Idx =
            std::uniform_int_distribution<size_t>{ 0,
                                                   neighbours.size() - 1 }(rng);
          a2 = neighbours[randomA2Idx];
        }
        // finally, remove them from the neighbour lists so that they are not
        // sampled more than once
        size_t sampledVertexId2 = this->universe.getIdxByAtomId(a2.getId());
        neighbourList.removeAtom(a2);
        isMasked[sampledVertexId2] = true;
        // it is actually quite a lot of expensive stuff done until we get to
        // this check but only this way we have the balance of removing atoms
        // to sample them only once
        if (!vertexIdxIsEligible[sampledVertexId] ||
            !vertexIdxIsEligible[sampledVertexId2]) {
          // std::cout << "Sampled vertices are not eligible" << std::endl;
          continue;
        }
        // take the mean and their index etc. to add as slip-link
        std::array<double, 3> meanPositions = a1.meanPositionWith(a2, &box);
        slipLinkXs.push_back(meanPositions[0]);
        slipLinkYs.push_back(meanPositions[1]);
        slipLinkZs.push_back(meanPositions[2]);
        slipLinkStrandA.push_back(atomToStrand.at(a1.getId()));
        slipLinkStrandB.push_back(atomToStrand.at(a2.getId()));
        slipLinkStrandAlpha.push_back(
          static_cast<double>(atomIdxInStrand.at(a1.getId())) /
          (static_cast<double>(
            this->initialConfig
              .springsContourLength[atomToStrand.at(a1.getId())])));
        slipLinkStrandBeta.push_back(
          static_cast<double>(atomIdxInStrand.at(a2.getId())) /
          (static_cast<double>(
            this->initialConfig
              .springsContourLength[atomToStrand.at(a2.getId())])));
        nrOfSlipLinksPlaced += 1;
      }

      RUNTIME_EXP_IFN(
        slipLinkStrandA.size() == slipLinkStrandB.size() &&
          slipLinkXs.size() == slipLinkYs.size() &&
          slipLinkZs.size() == slipLinkXs.size(),
        "Expect all slip-link relevant properties to have the same length");
      // with the data assembled, we can actually add them to the structure
      // and stuff
      this->addSlipLinks(slipLinkStrandA,
                         slipLinkStrandB,
                         slipLinkXs,
                         slipLinkYs,
                         slipLinkZs,
                         slipLinkStrandAlpha,
                         slipLinkStrandBeta,
                         false);
      return nrOfSlipLinksPlaced;
    }

    /**
     * @brief Find cycles (loops) and set slip-links based on those
     *
     * @param maxLoopLength
     * @return size_t
     */
    size_t MEHPForceBalance::addSliplinksBasedOnCycles(const int maxLoopLength)
    {
      std::vector<std::vector<long int>> loops =
        this->universe.findLoops(this->crosslinkerType, maxLoopLength, false);
      std::cout << "Detected " << loops.size() << " loops." << std::endl;
      std::vector<std::vector<size_t>> reducedLoops;
      reducedLoops.reserve(loops.size());
      for (std::vector<long int> loop : loops) {
        std::set<size_t> reducedLoop;
        for (size_t i = 0; i < loop.size(); ++i) {
          if (this->universe.getPropertyValue<long int>("type", loop[i]) !=
              this->crosslinkerType) {
            // set takes care of duplicates, yet this is not efficient at all.
            reducedLoop.insert(
              this->initialConfig.oldAtomIdToSpringIndex
                [this->universe.getPropertyValue<long int>("id", loop[i])]);
          }
        }
        std::vector<size_t> reducedLoopVec(reducedLoop.begin(),
                                           reducedLoop.end());
        reducedLoops.push_back(reducedLoopVec);
      }

      size_t estimateOfNrOfSliplinks =
        loops.size() * loops.size() *
        this->initialConfig.meanSpringContourLength;

      // fetch some data to later estimate alpha & beta
      std::vector<pylimer_tools::entities::Molecule> crosslinkerChains =
        this->universe.getChainsWithCrosslinker(crosslinkerType);
      std::unordered_map<size_t, size_t> atomToStrand;
      atomToStrand.reserve(this->universe.getNrOfAtoms());
      std::unordered_map<size_t, size_t> atomIdxInStrand;
      atomIdxInStrand.reserve(this->universe.getNrOfAtoms());
      size_t springId = 0;
      for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
        pylimer_tools::entities::Molecule chain = crosslinkerChains[i];
        std::vector<pylimer_tools::entities::Atom> atoms =
          crosslinkerChains[i].getAtomsLinedUp(this->crosslinkerType);
        for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
          pylimer_tools::entities::Atom atom = atoms[atomIdx];
          if (atom.getType() != this->crosslinkerType) {
            atomToStrand.emplace(atom.getId(), springId);
            atomIdxInStrand.emplace(atom.getId(), atomIdx);
          }
        }
      }

      // the resulting vectors to fill
      std::vector<double> slipLinkXs;
      slipLinkXs.reserve(estimateOfNrOfSliplinks);
      std::vector<double> slipLinkYs;
      slipLinkYs.reserve(estimateOfNrOfSliplinks);
      std::vector<double> slipLinkZs;
      slipLinkZs.reserve(estimateOfNrOfSliplinks);
      std::vector<size_t> slipLinkStrandA;
      slipLinkStrandA.reserve(estimateOfNrOfSliplinks);
      std::vector<size_t> slipLinkStrandB;
      slipLinkStrandB.reserve(estimateOfNrOfSliplinks);
      std::vector<double> slipLinkStrandAlpha;
      slipLinkStrandAlpha.reserve(estimateOfNrOfSliplinks);
      std::vector<double> slipLinkStrandBeta;
      slipLinkStrandBeta.reserve(estimateOfNrOfSliplinks);
      std::vector<std::vector<size_t>> slipLinksLoops;
      slipLinksLoops.reserve(estimateOfNrOfSliplinks);

      std::cout << "Searching for intersections..." << std::endl;
      for (size_t i = 0; i < loops.size(); ++i) {
        // TODO: instead of the N^2 loop, might want to try some filter at least
        // also, we ignore all self-entanglements of one loop with itself.
        for (size_t j = i + 1; j < loops.size(); ++j) {
          std::vector<pylimer_tools::entities::LoopIntersectionInfo>
            intersections =
              this->universe.findLoopEntanglements(loops[i], loops[j]);
          for (pylimer_tools::entities::LoopIntersectionInfo intersection :
               intersections) {
            // TODO: this is yet the most naïve way to add these.
            // ideally, we would also check the back-and-forth, etc.
            slipLinkXs.push_back(intersection.intersectionPoint[0]);
            slipLinkYs.push_back(intersection.intersectionPoint[1]);
            slipLinkZs.push_back(intersection.intersectionPoint[2]);
            // TODO: decide on the atoms to use as a reference
            slipLinkStrandA.push_back(
              this->initialConfig
                .oldAtomIdToSpringIndex[intersection.involvedAtoms[0].getId()]);
            slipLinkStrandB.push_back(
              this->initialConfig
                .oldAtomIdToSpringIndex[intersection.involvedAtoms[3].getId()]);
            slipLinkStrandAlpha.push_back(
              static_cast<double>(
                atomIdxInStrand[intersection.involvedAtoms[0].getId()]) /
              this->initialConfig.meanSpringContourLength);
            slipLinkStrandBeta.push_back(
              static_cast<double>(
                atomIdxInStrand[intersection.involvedAtoms[3].getId()]) /
              this->initialConfig.meanSpringContourLength);
            std::vector<size_t> localLoops = { i, j };
            slipLinksLoops.push_back(localLoops);
          }
        }
      }
      std::cout << "Found " << slipLinksLoops.size() << " intersections." << std::endl;
      RUNTIME_EXP_IFN(
        slipLinkStrandA.size() == slipLinkStrandB.size() &&
          slipLinkXs.size() == slipLinkYs.size() &&
          slipLinkZs.size() == slipLinkXs.size(),
        "Expect all slip-link relevant properties to have the same length");
      // with the data assembled, we can actually add them to the structure
      // and stuff
      this->addSlipLinks(slipLinkStrandA,
                         slipLinkStrandB,
                         slipLinkXs,
                         slipLinkYs,
                         slipLinkZs,
                         slipLinkStrandAlpha,
                         slipLinkStrandBeta,
                         reducedLoops,
                         slipLinksLoops,
                         false);
      return slipLinkStrandA.size();
    }

    /**
     * @brief Remove a spring (and all its parts, incl. slip-links) from the
     * structures
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance::removeSpring(ForceBalanceNetwork& net,
                                        Eigen::VectorXd& displacements,
                                        Eigen::VectorXd& springPartitions,
                                        const size_t springIdx) const
    {
      // std::cout << "Starting to remove spring " << springIdx << std::endl;
      INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                         "Can only remove springs, not partial springs.");
      std::vector<size_t> affectedLinks = net.linkIndicesOfSprings[springIdx];
      std::vector<size_t> uniqueAffectedLinks =
        net.linkIndicesOfSprings[springIdx];
      std::sort(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end());
      uniqueAffectedLinks.erase(
        std::unique(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end()),
        uniqueAffectedLinks.end());

      // remove the link to the link, höhö
      for (size_t affectedLinkIdx : uniqueAffectedLinks) {
        bool found = false;
        RUNTIME_EXP_IFN(
          std::find(net.springIndicesOfLinks[affectedLinkIdx].begin(),
                    net.springIndicesOfLinks[affectedLinkIdx].end(),
                    springIdx) !=
            net.springIndicesOfLinks[affectedLinkIdx].end(),
          "Link must have a connection to the spring, too. Did not find "
          "spring " +
            std::to_string(springIdx) + " in link " +
            std::to_string(affectedLinkIdx) + ", got " +
            pylimer_tools::utils::join(
              net.springIndicesOfLinks[affectedLinkIdx].begin(),
              net.springIndicesOfLinks[affectedLinkIdx].end(),
              std::string(", ")) +
            ".");
        if (net.linkIsSliplink[affectedLinkIdx]) {
          RUNTIME_EXP_IFN(
            net.springIndicesOfLinks[affectedLinkIdx].size() <= 2,
            "Expect slip-link to be associated with 2 springs only, got " +
              pylimer_tools::utils::join(
                net.springIndicesOfLinks[affectedLinkIdx].begin(),
                net.springIndicesOfLinks[affectedLinkIdx].end(),
                std::string(", ")) +
              ".");
        }

        for (int i = net.springIndicesOfLinks[affectedLinkIdx].size() - 1;
             i >= 0;
             --i) {
          if (net.springIndicesOfLinks[affectedLinkIdx][i] == springIdx) {
            net.springIndicesOfLinks[affectedLinkIdx].erase(
              net.springIndicesOfLinks[affectedLinkIdx].begin() + i);
            found = true;
          }
        }

        RUNTIME_EXP_IFN(found,
                        "Expected to find spring " + std::to_string(springIdx) +
                          " in link " + std::to_string(affectedLinkIdx) +
                          " but did not, got " +
                          pylimer_tools::utils::join(
                            net.springIndicesOfLinks[affectedLinkIdx].begin(),
                            net.springIndicesOfLinks[affectedLinkIdx].end(),
                            std::string(", ")) +
                          ".");
        if (net.linkIsSliplink[affectedLinkIdx]) {
          RUNTIME_EXP_IFN(net.springIndicesOfLinks[affectedLinkIdx].size() <= 1,
                          "Expect slip-link to be associated with 1 springs "
                          "only after removing one, got " +
                            pylimer_tools::utils::join(
                              net.springIndicesOfLinks[affectedLinkIdx].begin(),
                              net.springIndicesOfLinks[affectedLinkIdx].end(),
                              std::string(", ")) +
                            ".");
        }
      }

      std::vector<size_t> affectedPartialSprings =
        net.localToGlobalSpringIndex[springIdx];
      assert(affectedPartialSprings.size() > 0);
      net.nrOfSprings -= 1;
      net.nrOfPartialSprings -= affectedPartialSprings.size();

      // actually spring remove stuff
      net.localToGlobalSpringIndex.erase(net.localToGlobalSpringIndex.begin() +
                                         springIdx);
      net.springToMoleculeIds.erase(net.springToMoleculeIds.begin() +
                                    springIdx);
      net.linkIndicesOfSprings.erase(net.linkIndicesOfSprings.begin() +
                                     springIdx);
      pylimer_tools::utils::removeRow(net.springsContourLength, springIdx);
      pylimer_tools::utils::removeRow(net.springIndexA, springIdx);
      pylimer_tools::utils::removeRows(
        net.springCoordinateIndexA, springIdx * 3, 3);
      pylimer_tools::utils::removeRow(net.springIndexB, springIdx);
      pylimer_tools::utils::removeRows(
        net.springCoordinateIndexB, springIdx * 3, 3);
      pylimer_tools::utils::removeRow(net.springIsActive, springIdx);

      // need to remove descending
      std::sort(affectedPartialSprings.begin(),
                affectedPartialSprings.end(),
                std::greater<size_t>());
      for (size_t partialSpringIdx : affectedPartialSprings) {
        pylimer_tools::utils::removeRow(net.springPartIndexA, partialSpringIdx);
        pylimer_tools::utils::removeRows(
          net.springPartCoordinateIndexA, 3 * partialSpringIdx, 3);
        pylimer_tools::utils::removeRow(net.springPartIndexB, partialSpringIdx);
        pylimer_tools::utils::removeRows(
          net.springPartCoordinateIndexB, 3 * partialSpringIdx, 3);
        pylimer_tools::utils::removeRow(net.partialToFullSpringIndex,
                                        partialSpringIdx);
        pylimer_tools::utils::removeRow(net.partialSpringIsPartial,
                                        partialSpringIdx);
        pylimer_tools::utils::removeRow(springPartitions, partialSpringIdx);
      }
      assert(springPartitions.size() == net.nrOfPartialSprings);

      // renumber the remaining stuff
      // first, renumber the springs
      for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
        for (size_t i = 0; i < net.springIndicesOfLinks[linkIdx].size(); ++i) {
          assert(net.springIndicesOfLinks[linkIdx][i] != springIdx);
          if (net.springIndicesOfLinks[linkIdx][i] > springIdx) {
            net.springIndicesOfLinks[linkIdx][i] -= 1;
          }
        }
      }

      // then, renumber the loops
      for (size_t loopIdx = 0; loopIdx < net.loops.size(); ++loopIdx) {
        for (size_t i = 0; i < net.loops[loopIdx].size(); ++i) {
          if (net.loops[loopIdx][i] == springIdx) {
            net.loops[loopIdx].erase(net.loops[loopIdx].begin() + i);
          }
          if (net.loops[loopIdx][i] > springIdx) {
            net.loops[loopIdx][i] -= 1;
          }
        }
      }

      // then, update the partial springs
      for (size_t partialSpringIdx = 0;
           partialSpringIdx < net.nrOfPartialSprings;
           ++partialSpringIdx) {
        assert(net.partialToFullSpringIndex[partialSpringIdx] != springIdx);
        if (net.partialToFullSpringIndex[partialSpringIdx] > springIdx) {
          net.partialToFullSpringIndex[partialSpringIdx] -= 1;
        }
      }

      assert(net.localToGlobalSpringIndex.size() == net.nrOfSprings);
      for (size_t loopingSpringIdx = 0; loopingSpringIdx < net.nrOfSprings;
           ++loopingSpringIdx) {
        for (size_t i = 0;
             i < net.localToGlobalSpringIndex[loopingSpringIdx].size();
             ++i) {
          for (size_t partialSpringIdx : affectedPartialSprings) {
            if (net.localToGlobalSpringIndex[loopingSpringIdx][i] >
                partialSpringIdx) {
              net.localToGlobalSpringIndex[loopingSpringIdx][i] -= 1;
            }
          }
        }
      }

      //
      // remove the affected slip-links that are now only on one spring
      //
      std::vector<size_t> linksToRemove;
      for (size_t i = 1; i < affectedLinks.size() - 1; ++i) {
        linksToRemove.push_back(affectedLinks[i]);
      }

      // need to remove descending to remove need to renumber these as well
      std::sort(
        linksToRemove.begin(), linksToRemove.end(), std::greater<size_t>());
      linksToRemove.erase(
        std::unique(linksToRemove.begin(), linksToRemove.end()),
        linksToRemove.end());

      for (size_t outermostI = 0; outermostI < linksToRemove.size();
           ++outermostI) {
        size_t slipLinkIdx = linksToRemove[outermostI];
        assert(net.linkIsSliplink[slipLinkIdx]);
        // first, merge the two other partial springs
        std::vector<size_t> springsOfLink =
          net.springIndicesOfLinks[slipLinkIdx];
        RUNTIME_EXP_IFN(springsOfLink.size() <= 1,
                        "Expected slip-link " + std::to_string(slipLinkIdx) +
                          " to have only 1 remaining spring, got " +
                          std::to_string(springsOfLink.size()) + " due to " +
                          pylimer_tools::utils::join(springsOfLink.begin(),
                                                     springsOfLink.end(),
                                                     std::string(", ")) +
                          ".");
        std::vector<size_t> involvedPartialSprings;
        involvedPartialSprings.reserve(2);
        for (int springInLinkIdx = springsOfLink.size() - 1;
             springInLinkIdx >= 0;
             --springInLinkIdx) {
          for (size_t partialSpringIdx :
               net.localToGlobalSpringIndex[springsOfLink[springInLinkIdx]]) {
            if (net.springPartIndexA[partialSpringIdx] == slipLinkIdx ||
                net.springPartIndexB[partialSpringIdx] == slipLinkIdx) {
              involvedPartialSprings.push_back(partialSpringIdx);
            }
          }
        }
        RUNTIME_EXP_IFN(
          involvedPartialSprings.size() >= springsOfLink.size(),
          "Expected more or equal number of partial springs (" +
            std::to_string(involvedPartialSprings.size()) + "; " +
            pylimer_tools::utils::join(involvedPartialSprings.begin(),
                                       involvedPartialSprings.end(),
                                       std::string(", ")) +
            ") than springs (" + std::to_string(springsOfLink.size()) + "; " +
            pylimer_tools::utils::join(
              springsOfLink.begin(), springsOfLink.end(), std::string(", ")) +
            ").");
        // RUNTIME_EXP_IFN(springsOfLink.size() % 2 == 0, "Expected link to have
        // an even number of components, got " +
        // std::to_string(springsOfLink.size()) + ".");
        if (involvedPartialSprings.size() > 0) {
          RUNTIME_EXP_IFN(
            involvedPartialSprings.size() == 2,
            "Expected only 2 involved partial springs, got: " +
              pylimer_tools::utils::join(involvedPartialSprings.begin(),
                                         involvedPartialSprings.end(),
                                         std::string(", ")) +
              " for springs " +
              pylimer_tools::utils::join(
                springsOfLink.begin(), springsOfLink.end(), std::string(", ")) +
              " when removing spring " + std::to_string(springIdx) + ".");
          assert(involvedPartialSprings.size() == 2);
          size_t partialSpringToKeep = involvedPartialSprings[0];
          size_t springToKeepIdx =
            net.partialToFullSpringIndex[partialSpringToKeep];
          size_t partialSpringToRemove = involvedPartialSprings[1];
          assert(partialSpringToKeep != partialSpringToRemove);
          assert(net.partialToFullSpringIndex[partialSpringToKeep] ==
                 net.partialToFullSpringIndex[partialSpringToRemove]);
          // actually do the merge
          this->mergePartialSprings(net,
                                    springPartitions,
                                    partialSpringToRemove,
                                    partialSpringToKeep,
                                    slipLinkIdx);
        }

        assert(net.springIndicesOfLinks[slipLinkIdx].empty());

        // then, actually remove the slip-link
        // std::cout << "Removing link " << slipLinkIdx << std::endl;
        this->removeLink(net, displacements, slipLinkIdx);
      }
      // std::cout << "Removed spring " << springIdx << std::endl;
    }

    /**
     * @brief
     *
     * @param net
     * @param displacements
     * @param linkIdx
     */
    void MEHPForceBalance::removeLink(ForceBalanceNetwork& net,
                                      Eigen::VectorXd& displacements,
                                      const size_t linkIdx) const
    {
      INVALIDARG_EXP_IFN(net.springIndicesOfLinks[linkIdx].size() == 0,
                         "Please remove the springs before removing the link.");

      pylimer_tools::utils::removeRows(net.coordinates, linkIdx * 3, 3);
      pylimer_tools::utils::removeRows(displacements, linkIdx * 3, 3);

      if (!net.linkIsSliplink[linkIdx]) {
        net.nrOfNodes -= 1;
        net.oldAtomIdToSpringIndex.erase(net.oldAtomIds[linkIdx]);
        pylimer_tools::utils::removeRow(net.oldAtomIds, linkIdx);
      } else {
        pylimer_tools::utils::removeRow(net.nrOfCrosslinkSwapsEndured,
                                        linkIdx - net.nrOfNodes);
        net.loopsOfSliplink.erase(net.loopsOfSliplink.begin() + linkIdx);
      }
      net.nrOfLinks -= 1;
      pylimer_tools::utils::removeRow(net.linkIsSliplink, linkIdx);
      net.springIndicesOfLinks.erase(net.springIndicesOfLinks.begin() +
                                     linkIdx);

      // renumber the remaining links
      for (size_t i = 0; i < net.linkIndicesOfSprings.size(); ++i) {
        for (size_t j = 0; j < net.linkIndicesOfSprings[i].size(); ++j) {
          RUNTIME_EXP_IFN(
            net.linkIndicesOfSprings[i][j] != linkIdx,
            "Expected not to find link to remove " + std::to_string(linkIdx) +
              " in any spring, found in spring " + std::to_string(i) + ", " +
              pylimer_tools::utils::join(net.linkIndicesOfSprings[i].begin(),
                                         net.linkIndicesOfSprings[i].end(),
                                         std::string(", ")) +
              ".");
          if (net.linkIndicesOfSprings[i][j] > linkIdx) {
            net.linkIndicesOfSprings[i][j] -= 1;
          }
        }
      }
      //
      assert(net.springPartIndexA.size() == net.springPartIndexB.size());
      for (size_t i = 0; i < net.springPartIndexA.size(); ++i) {
        RUNTIME_EXP_IFN(
          net.springPartIndexA[i] != linkIdx,
          "Exected link " + std::to_string(linkIdx) +
            " to not be linked anywhere anymore, found in partial spring " +
            std::to_string(i) + ".");
        if (net.springPartIndexA[i] > linkIdx) {
          net.springPartIndexA[i] -= 1;
          net.springPartCoordinateIndexA[3 * i] -= 3;
          net.springPartCoordinateIndexA[3 * i + 1] -= 3;
          net.springPartCoordinateIndexA[3 * i + 2] -= 3;
        }
        RUNTIME_EXP_IFN(
          net.springPartIndexB[i] != linkIdx,
          "Exected link " + std::to_string(linkIdx) +
            " to not be linked anywhere anymore, found in partial spring " +
            std::to_string(i) + ".");
        if (net.springPartIndexB[i] > linkIdx) {
          net.springPartIndexB[i] -= 1;
          net.springPartCoordinateIndexB[3 * i] -= 3;
          net.springPartCoordinateIndexB[3 * i + 1] -= 3;
          net.springPartCoordinateIndexB[3 * i + 2] -= 3;
        }
      }
      //
      assert(net.springIndexA.size() == net.springIndexB.size());
      for (size_t i = 0; i < net.springIndexA.size(); ++i) {
        assert(net.springIndexA[i] != linkIdx);
        if (net.springIndexA[i] > linkIdx) {
          net.springIndexA[i] -= 1;
          net.springCoordinateIndexA[3 * i] -= 3;
          net.springCoordinateIndexA[3 * i + 1] -= 3;
          net.springCoordinateIndexA[3 * i + 2] -= 3;
        }
        assert(net.springIndexB[i] != linkIdx);
        if (net.springIndexB[i] > linkIdx) {
          net.springIndexB[i] -= 1;
          net.springCoordinateIndexB[3 * i] -= 3;
          net.springCoordinateIndexB[3 * i + 1] -= 3;
          net.springCoordinateIndexB[3 * i + 2] -= 3;
        }
      }
    }

    /**
     * @brief
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance::mergePartialSprings(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& springPartitions,
      const size_t removedSpringIdx,
      const size_t keptSpringIdx,
      const size_t linkToReduce,
      bool skipEigenResize) const
    {
      INVALIDARG_EXP_IFN(net.linkIsSliplink[linkToReduce],
                         "The link to reduce must be a slip-link");
      INVALIDARG_EXP_IFN(keptSpringIdx != removedSpringIdx,
                         "Cannot merge one spring with the same one.");
      INVALIDARG_EXP_IFN(
        net.partialToFullSpringIndex[keptSpringIdx] ==
          net.partialToFullSpringIndex[removedSpringIdx],
        "The partial springs must be part of the same spring to merge them.")
      INVALIDARG_EXP_IFN(
        net.springPartIndexA[removedSpringIdx] == linkToReduce ||
          net.springPartIndexB[removedSpringIdx] == linkToReduce,
        "Link to reduce must be part of the springs to partial spring that is "
        "removed, got " +
          std::to_string(net.springPartIndexA[removedSpringIdx]) + " and " +
          std::to_string(net.springPartIndexB[removedSpringIdx]) +
          " instead of " + std::to_string(linkToReduce) + ".");
      INVALIDARG_EXP_IFN(net.springPartIndexA[keptSpringIdx] == linkToReduce ||
                           net.springPartIndexB[keptSpringIdx] == linkToReduce,
                         "Link to reduce must be part of the partial spring "
                         "that is to be kept, got " +
                           std::to_string(net.springPartIndexA[keptSpringIdx]) +
                           " and " +
                           std::to_string(net.springPartIndexB[keptSpringIdx]) +
                           " instead of " + std::to_string(linkToReduce) + ".");
      size_t fullSpringIdx = net.partialToFullSpringIndex[keptSpringIdx];
      // start with removal
      net.nrOfPartialSprings -= 1;
      // tell the kept one their new end
      size_t newEnd = net.springPartIndexA[removedSpringIdx] == linkToReduce
                        ? net.springPartIndexB[removedSpringIdx]
                        : net.springPartIndexA[removedSpringIdx];
      if (net.springPartIndexA[keptSpringIdx] == linkToReduce) {
        net.springPartIndexA[keptSpringIdx] = newEnd;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexA[3 * keptSpringIdx + dir] =
            3 * newEnd + dir;
        }
      } else {
        net.springPartIndexB[keptSpringIdx] = newEnd;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexB[3 * keptSpringIdx + dir] =
            3 * newEnd + dir;
        }
      }
      // remove the spring from the link
      // NOTE: currently, we allow it not to be present,
      // as it might be removed earlier already
      // It is anyway the case, that this function does not necessarily
      // keep the network valid
      int found = 0;
      for (int i = net.springIndicesOfLinks[linkToReduce].size() - 1; i >= 0;
           --i) {
        if (net.springIndicesOfLinks[linkToReduce][i] == fullSpringIdx) {
          net.springIndicesOfLinks[linkToReduce].erase(
            net.springIndicesOfLinks[linkToReduce].begin() + i);
          found += 1;
        }
      }
      assert(found == 1 || found == 0);
      RUNTIME_EXP_IFN(
        net.localToGlobalSpringIndex[fullSpringIdx].size() ==
          net.linkIndicesOfSprings[fullSpringIdx].size() - 1,
        "Require a global index for each local one, got " +
          std::to_string(net.localToGlobalSpringIndex[fullSpringIdx].size()) +
          " != " +
          std::to_string(net.linkIndicesOfSprings[fullSpringIdx].size() - 1) +
          " for spring " + std::to_string(fullSpringIdx) + ".");
      found = 0;
      // tell the spring of the removed link
      int removed = 0;
      springPartitions[keptSpringIdx] += springPartitions[removedSpringIdx];
      for (int j = net.linkIndicesOfSprings[fullSpringIdx].size() - 1; j >= 0;
           --j) {
        if (net.linkIndicesOfSprings[fullSpringIdx][j] == linkToReduce) {
          if (removed == 0) {
            if ((j > 0 && net.localToGlobalSpringIndex[fullSpringIdx][j - 1] ==
                            removedSpringIdx) ||
                (j < net.localToGlobalSpringIndex[fullSpringIdx].size() &&
                 net.localToGlobalSpringIndex[fullSpringIdx][j] ==
                   removedSpringIdx)) {
              net.linkIndicesOfSprings[fullSpringIdx].erase(
                net.linkIndicesOfSprings[fullSpringIdx].begin() + j);
              removed += 1;
            }
          }
          if (found == 1) {
            // we are dealing with a double -> re-add
            net.springIndicesOfLinks[linkToReduce].push_back(fullSpringIdx);
          } else if (found > 1 && removed > 0) {
            break; // required for certain cases... dangerous, somewhat.
          }
          found += 1;
        }
      }
      assert(found >= 1 && removed == 1);
      found = 0;
      for (int j = net.localToGlobalSpringIndex[fullSpringIdx].size() - 1;
           j >= 0;
           --j) {
        if (net.localToGlobalSpringIndex[fullSpringIdx][j] ==
            removedSpringIdx) {
          net.localToGlobalSpringIndex[fullSpringIdx].erase(
            net.localToGlobalSpringIndex[fullSpringIdx].begin() + j);
          found += 1;
        }
      }
      assert(found == 1);
      RUNTIME_EXP_IFN(
        net.localToGlobalSpringIndex[fullSpringIdx].size() ==
          net.linkIndicesOfSprings[fullSpringIdx].size() - 1,
        "Require a global index for each local one, got " +
          std::to_string(net.localToGlobalSpringIndex[fullSpringIdx].size()) +
          " != " +
          std::to_string(net.linkIndicesOfSprings[fullSpringIdx].size() - 1) +
          " for spring " + std::to_string(fullSpringIdx) + ".");
      // recompute some values
      net.partialSpringIsPartial[keptSpringIdx] =
        net.linkIndicesOfSprings[fullSpringIdx].size() > 2;
      // actually remove the rows
      pylimer_tools::utils::removeRow(
        net.partialSpringIsPartial, removedSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.partialToFullSpringIndex, removedSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        springPartitions, removedSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.springPartIndexA, removedSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.springPartIndexB, removedSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRows(net.springPartCoordinateIndexA,
                                       3 * removedSpringIdx,
                                       3,
                                       skipEigenResize);
      pylimer_tools::utils::removeRows(net.springPartCoordinateIndexB,
                                       3 * removedSpringIdx,
                                       3,
                                       skipEigenResize);
      // renumber stuff
      for (size_t loopSpringIdx = 0;
           loopSpringIdx < net.localToGlobalSpringIndex.size();
           ++loopSpringIdx) {
        for (size_t i = 0;
             i < net.localToGlobalSpringIndex[loopSpringIdx].size();
             ++i) {
          if (net.localToGlobalSpringIndex[loopSpringIdx][i] >
              removedSpringIdx) {
            net.localToGlobalSpringIndex[loopSpringIdx][i] -= 1;
          }
        }
      }
    }

    /**
     * @brief
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance::mergeSprings(ForceBalanceNetwork& net,
                                        Eigen::VectorXd& springPartitions,
                                        const size_t removedSpringIdx,
                                        const size_t keptSpringIdx,
                                        const size_t linkToReduce) const
    {
      INVALIDARG_EXP_IFN(removedSpringIdx < net.nrOfSprings &&
                           keptSpringIdx < net.nrOfSprings,
                         "Only full springs can be merged.");
      INVALIDARG_EXP_IFN(!net.linkIsSliplink[linkToReduce],
                         "The link to reduce must be a cross-link");
      INVALIDARG_EXP_IFN(keptSpringIdx != removedSpringIdx,
                         "Cannot replace one spring with the same one.");
      net.nrOfSprings -= 1;
      net.nrOfPartialSprings -= 1;
      if (net.linkIndicesOfSprings[removedSpringIdx].size() > 2 &&
          net.linkIndicesOfSprings[keptSpringIdx].size() > 2) {
        net.nrOfSpringsWithPartition -= 1;
      }
      // handle links
      std::vector<size_t> removedSpringsLinks =
        net.linkIndicesOfSprings[removedSpringIdx];
      std::vector<size_t> keptSpringsLinks =
        net.linkIndicesOfSprings[keptSpringIdx];

      size_t removedPartialSpringIdx =
        (removedSpringsLinks[removedSpringsLinks.size() - 1] == linkToReduce)
          ? pylimer_tools::utils::last(
              net.localToGlobalSpringIndex[removedSpringIdx])
          : net.localToGlobalSpringIndex[removedSpringIdx][0];
      size_t remainingPartialSpringIdx =
        (keptSpringsLinks[keptSpringsLinks.size() - 1] == linkToReduce)
          ? pylimer_tools::utils::last(
              net.localToGlobalSpringIndex[keptSpringIdx])
          : net.localToGlobalSpringIndex[keptSpringIdx][0];

      net.linkIndicesOfSprings[keptSpringIdx].reserve(
        keptSpringsLinks.size() + removedSpringsLinks.size() - 2);
      net.localToGlobalSpringIndex[keptSpringIdx].reserve(
        keptSpringsLinks.size() + removedSpringsLinks.size() - 2);
      RUNTIME_EXP_IFN(net.localToGlobalSpringIndex[keptSpringIdx].size() ==
                        net.linkIndicesOfSprings[keptSpringIdx].size() - 1,
                      "Invalid sizes when merging springs");
      // tell the partial springs their new full spring
      for (size_t partialSpringIndex :
           net.localToGlobalSpringIndex[removedSpringIdx]) {
        net.partialToFullSpringIndex[partialSpringIndex] = keptSpringIdx;
      }
      // std::cout << "Kept spring is "
      //           << pylimer_tools::utils::join(keptSpringsLinks.begin(),
      //                                         keptSpringsLinks.end(),
      //                                         std::string(", "))
      //           << std::endl;
      // actually merge the springs
      if (keptSpringsLinks[keptSpringsLinks.size() - 1] == linkToReduce) {
        // add to end...
        if (removedSpringsLinks[removedSpringsLinks.size() - 1] ==
            linkToReduce) {
          // std::cout << "End end" << std::endl;
          // from end
          net.linkIndicesOfSprings[keptSpringIdx][keptSpringsLinks.size() - 1] =
            removedSpringsLinks[removedSpringsLinks.size() - 2];
          for (size_t i = 3; i <= removedSpringsLinks.size(); ++i) {
            net.linkIndicesOfSprings[keptSpringIdx].push_back(
              removedSpringsLinks[removedSpringsLinks.size() - i]);
          }
          for (int i =
                 net.localToGlobalSpringIndex[removedSpringIdx].size() - 2;
               i >= 0;
               --i) {
            net.localToGlobalSpringIndex[keptSpringIdx].push_back(
              net.localToGlobalSpringIndex[removedSpringIdx][i]);
          }
        } else {
          // from start
          // std::cout << "End start" << std::endl;
          RUNTIME_EXP_IFN(removedSpringsLinks[0] == linkToReduce,
                          "Things don't make sense anymore.");
          net.linkIndicesOfSprings[keptSpringIdx][keptSpringsLinks.size() - 1] =
            removedSpringsLinks[1];
          for (size_t i = 2; i < removedSpringsLinks.size(); ++i) {
            net.linkIndicesOfSprings[keptSpringIdx].push_back(
              removedSpringsLinks[i]);
          }
          for (size_t i = 1;
               i < net.localToGlobalSpringIndex[removedSpringIdx].size();
               ++i) {
            net.localToGlobalSpringIndex[keptSpringIdx].push_back(
              net.localToGlobalSpringIndex[removedSpringIdx][i]);
          }
        }
      } else {
        RUNTIME_EXP_IFN(keptSpringsLinks[0] == linkToReduce,
                        "How could this be?");
        // add to start...
        if (removedSpringsLinks[removedSpringsLinks.size() - 1] ==
            linkToReduce) {
          // std::cout << "Start end" << std::endl;
          // from end
          net.linkIndicesOfSprings[keptSpringIdx][0] =
            removedSpringsLinks[removedSpringsLinks.size() - 2];
          for (size_t i = 3; i <= removedSpringsLinks.size(); ++i) {
            net.linkIndicesOfSprings[keptSpringIdx].insert(
              net.linkIndicesOfSprings[keptSpringIdx].begin(),
              removedSpringsLinks[removedSpringsLinks.size() - i]);
          }
          for (int i =
                 net.localToGlobalSpringIndex[removedSpringIdx].size() - 2;
               i >= 0;
               --i) {
            net.localToGlobalSpringIndex[keptSpringIdx].insert(
              net.localToGlobalSpringIndex[keptSpringIdx].begin(),
              net.localToGlobalSpringIndex[removedSpringIdx][i]);
          }
        } else {
          // std::cout << "Start start" << std::endl;
          // from start
          RUNTIME_EXP_IFN(removedSpringsLinks[0] == linkToReduce,
                          "No way this expcetion is every shown, right?");
          net.linkIndicesOfSprings[keptSpringIdx][0] = removedSpringsLinks[1];
          // have to insert it reverse order
          // happens automatically if we always insert the next the start
          for (size_t i = 2; i < removedSpringsLinks.size(); ++i) {
            net.linkIndicesOfSprings[keptSpringIdx].insert(
              net.linkIndicesOfSprings[keptSpringIdx].begin(),
              removedSpringsLinks[i]);
          }
          // skip the first (removed) partial spring
          for (size_t i = 1;
               i < net.localToGlobalSpringIndex[removedSpringIdx].size();
               ++i) {
            net.localToGlobalSpringIndex[keptSpringIdx].insert(
              net.localToGlobalSpringIndex[keptSpringIdx].begin(),
              net.localToGlobalSpringIndex[removedSpringIdx][i]);
          }
        }
      }
      RUNTIME_EXP_IFN(
        std::find(net.linkIndicesOfSprings[keptSpringIdx].begin(),
                  net.linkIndicesOfSprings[keptSpringIdx].end(),
                  linkToReduce) ==
          net.linkIndicesOfSprings[keptSpringIdx].end(),
        "Link " + std::to_string(linkToReduce) +
          " to reduce should not be in the kept links anymore, found " +
          pylimer_tools::utils::join(
            net.linkIndicesOfSprings[keptSpringIdx].begin(),
            net.linkIndicesOfSprings[keptSpringIdx].end(),
            std::string(", ")) +
          ".");
      assert(net.localToGlobalSpringIndex[keptSpringIdx].size() ==
             net.linkIndicesOfSprings[keptSpringIdx].size() - 1);
      assert(net.linkIndicesOfSprings[keptSpringIdx].size() ==
             keptSpringsLinks.size() + removedSpringsLinks.size() - 2);

      // tell the links of their new spring index
      for (size_t linkOfRemovedSpring : removedSpringsLinks) {
        for (size_t i = 0;
             i < net.springIndicesOfLinks[linkOfRemovedSpring].size();
             ++i) {
          if (net.springIndicesOfLinks[linkOfRemovedSpring][i] ==
              removedSpringIdx) {
            net.springIndicesOfLinks[linkOfRemovedSpring][i] = keptSpringIdx;
          }
        }
      }

      for (int i = net.springIndicesOfLinks[linkToReduce].size() - 1; i >= 0;
           --i) {
        if (net.springIndicesOfLinks[linkToReduce][i] == removedSpringIdx ||
            net.springIndicesOfLinks[linkToReduce][i] == keptSpringIdx) {
          net.springIndicesOfLinks[linkToReduce].erase(
            net.springIndicesOfLinks[linkToReduce].begin() + i);
        }
      }
      net.linkIndicesOfSprings.erase(net.linkIndicesOfSprings.begin() +
                                     removedSpringIdx);
      // partial springs
      if (net.partialSpringIsPartial[removedPartialSpringIdx]) {
        net.partialSpringIsPartial[remainingPartialSpringIdx] = true;
      }
      pylimer_tools::utils::removeRow(net.partialSpringIsPartial,
                                      removedPartialSpringIdx);

      RUNTIME_EXP_IFN(
        net.springPartIndexA[removedPartialSpringIdx] == linkToReduce ||
          net.springPartIndexB[removedPartialSpringIdx] == linkToReduce,
        "");
      RUNTIME_EXP_IFN(
        net.springPartIndexA[remainingPartialSpringIdx] == linkToReduce ||
          net.springPartIndexB[remainingPartialSpringIdx] == linkToReduce,
        "");
      bool removedIsA =
        net.springPartIndexA[removedPartialSpringIdx] == linkToReduce;
      size_t otherEndOfRemovedSpring =
        removedIsA ? net.springPartIndexB[removedPartialSpringIdx]
                   : net.springPartIndexA[removedPartialSpringIdx];
      if (net.springPartIndexA[remainingPartialSpringIdx] == linkToReduce) {
        net.springPartIndexA[remainingPartialSpringIdx] =
          otherEndOfRemovedSpring;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexA[3 * remainingPartialSpringIdx + dir] =
            3 * otherEndOfRemovedSpring + dir;
        }
      } else {
        RUNTIME_EXP_IFN(
          net.springPartIndexB[remainingPartialSpringIdx] == linkToReduce, "");
        net.springPartIndexB[remainingPartialSpringIdx] =
          otherEndOfRemovedSpring;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexB[3 * remainingPartialSpringIdx + dir] =
            3 * otherEndOfRemovedSpring + dir;
        }
      }
      pylimer_tools::utils::removeRow(net.springPartIndexA,
                                      removedPartialSpringIdx);
      pylimer_tools::utils::removeRow(net.springPartIndexB,
                                      removedPartialSpringIdx);
      pylimer_tools::utils::removeRows(
        net.springPartCoordinateIndexA, 3 * removedPartialSpringIdx, 3);
      pylimer_tools::utils::removeRows(
        net.springPartCoordinateIndexB, 3 * removedPartialSpringIdx, 3);

      // spring indices & coordinates
      if (net.springIndexA[removedSpringIdx] == linkToReduce) {
        if (net.springIndexA[keptSpringIdx] == linkToReduce) {
          net.springIndexA[keptSpringIdx] = net.springIndexB[removedSpringIdx];
          net.springCoordinateIndexA.segment(3 * keptSpringIdx, 3) =
            net.springCoordinateIndexB.segment(3 * removedSpringIdx, 3);
        } else {
          assert(net.springIndexB[keptSpringIdx] == linkToReduce);
          net.springIndexB[keptSpringIdx] = net.springIndexB[removedSpringIdx];
          net.springCoordinateIndexB.segment(3 * keptSpringIdx, 3) =
            net.springCoordinateIndexB.segment(3 * removedSpringIdx, 3);
        }
      } else {
        assert(net.springIndexB[removedSpringIdx] == linkToReduce);
        if (net.springIndexA[keptSpringIdx] == linkToReduce) {
          net.springIndexA[keptSpringIdx] = net.springIndexA[removedSpringIdx];
          net.springCoordinateIndexA.segment(3 * keptSpringIdx, 3) =
            net.springCoordinateIndexA.segment(3 * removedSpringIdx, 3);
        } else {
          assert(net.springIndexB[keptSpringIdx] == linkToReduce);
          net.springIndexB[keptSpringIdx] = net.springIndexA[removedSpringIdx];
          net.springCoordinateIndexB.segment(3 * keptSpringIdx, 3) =
            net.springCoordinateIndexA.segment(3 * removedSpringIdx, 3);
        }
      }
      pylimer_tools::utils::removeRow(net.springIndexA, removedSpringIdx);
      pylimer_tools::utils::removeRow(net.springIndexB, removedSpringIdx);
      pylimer_tools::utils::removeRows(
        net.springCoordinateIndexA, 3 * removedSpringIdx, 3);
      pylimer_tools::utils::removeRows(
        net.springCoordinateIndexB, 3 * removedSpringIdx, 3);
      pylimer_tools::utils::removeRow(net.springIsActive, removedSpringIdx);
      net.springToMoleculeIds.erase(net.springToMoleculeIds.begin() +
                                    removedSpringIdx);
      net.oldAtomIdToSpringIndex.erase(net.oldAtomIds[linkToReduce]);

      pylimer_tools::utils::removeRow(net.partialToFullSpringIndex,
                                      removedPartialSpringIdx);
      net.localToGlobalSpringIndex.erase(net.localToGlobalSpringIndex.begin() +
                                         removedSpringIdx);
      // renumber the remaining springs
      for (size_t i = 0; i < net.springIndicesOfLinks.size(); ++i) {
        for (size_t j = 0; j < net.springIndicesOfLinks[i].size(); ++j) {
          RUNTIME_EXP_IFN(net.springIndicesOfLinks[i][j] != removedSpringIdx,
                          "");
          if (net.springIndicesOfLinks[i][j] > removedSpringIdx) {
            net.springIndicesOfLinks[i][j] -= 1;
          }
        }
      }

      // then, renumber the loops
      for (size_t loopIdx = 0; loopIdx < net.loops.size(); ++loopIdx) {
        for (size_t i = 0; i < net.loops[loopIdx].size(); ++i) {
          if (net.loops[loopIdx][i] == removedSpringIdx) {
            net.loops[loopIdx].erase(net.loops[loopIdx].begin() + i);
          }
          if (net.loops[loopIdx][i] > removedSpringIdx) {
            net.loops[loopIdx][i] -= 1;
          }
        }
      }

      // and the partial springs
      for (size_t i = 0; i < net.partialToFullSpringIndex.size(); ++i) {
        RUNTIME_EXP_IFN(net.partialToFullSpringIndex[i] != removedSpringIdx,
                        "");
        if (net.partialToFullSpringIndex[i] > removedSpringIdx) {
          net.partialToFullSpringIndex[i] -= 1;
        }
      }

      for (size_t i = 0; i < net.localToGlobalSpringIndex.size(); ++i) {
        for (size_t j = 0; j < net.localToGlobalSpringIndex[i].size(); ++j) {
          RUNTIME_EXP_IFN(
            net.localToGlobalSpringIndex[i][j] != removedPartialSpringIdx, "");
          if (net.localToGlobalSpringIndex[i][j] > removedPartialSpringIdx) {
            net.localToGlobalSpringIndex[i][j] -= 1;
          }
        }
      }

      // handle contour lengths
      double contourLengthBefore = net.springsContourLength[keptSpringIdx];
      net.springsContourLength[keptSpringIdx] +=
        net.springsContourLength[removedSpringIdx];
      pylimer_tools::utils::removeRow(net.springsContourLength,
                                      removedSpringIdx);
      RUNTIME_EXP_IFN(net.springsContourLength.size() == net.nrOfSprings, "");
      // and spring partitions
      springPartitions[remainingPartialSpringIdx] +=
        springPartitions[removedPartialSpringIdx];

      pylimer_tools::utils::removeRow(springPartitions,
                                      removedPartialSpringIdx);
      RUNTIME_EXP_IFN(springPartitions.size() == net.nrOfPartialSprings, "");
      size_t newKeptSpringIdx = (keptSpringIdx < removedSpringIdx)
                                  ? keptSpringIdx
                                  : (keptSpringIdx - 1);
      // addmittedly, this is possibly dangerous, as it could hide
      // other mistakes
      double newTotalForNormalization =
        springPartitions(net.localToGlobalSpringIndex[newKeptSpringIdx]).sum();
      for (size_t globalPartSpringIndex :
           net.localToGlobalSpringIndex[newKeptSpringIdx]) {
        springPartitions[globalPartSpringIndex] *=
          1. / newTotalForNormalization;
      }

      std::cout << "Removed springs around " << linkToReduce << " with spring "
                << removedSpringIdx << " and partial "
                << removedPartialSpringIdx << ", keeping " << keptSpringIdx
                << " and " << remainingPartialSpringIdx << std::endl;
      std::cout << "Spring partitions sum to " << springPartitions.sum()
                << " for " << net.nrOfSprings
                << " springs, contour length before was " << contourLengthBefore
                << " and is now " << net.springsContourLength[newKeptSpringIdx]
                << std::endl;
    }

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
    size_t MEHPForceBalance::doRemovalAndreisWay(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions,
      double tolerance) const
    {
      size_t numRemovedTotal = 0;
      size_t numRemovedInIteration = 0;
      do {
        numRemovedInIteration = 0;
        // do removal of f = 1
        // remove all cross-links that are 0- or 1-functional
        for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
             --crosslinkIdx) {
          if (net.springIndicesOfLinks[crosslinkIdx].size() == 0 // f = 0
          ) {
            // std::cout << "Removing x-link " << crosslinkIdx << std::endl;
            this->removeLink(net, displacements, crosslinkIdx);
            numRemovedInIteration += 1;
            // this->validateNetwork(net, displacements, springPartitions);
          }

          if ( // or f = 1, NOT primary loop
            net.springIndicesOfLinks[crosslinkIdx].size() == 1 &&
            XOR(
              net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx]
                                                               [0]][0] ==
                crosslinkIdx,
              pylimer_tools::utils::last(
                net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx]
                                                                 [0]]) ==
                crosslinkIdx)) {
            // need to first remove the spring
            this->removeSpring(net,
                               displacements,
                               springPartitions,
                               net.springIndicesOfLinks[crosslinkIdx][0]);
            // to then remove the link
            this->removeLink(net, displacements, crosslinkIdx);
            numRemovedInIteration += 1;
          }
        }
        numRemovedTotal += numRemovedInIteration;
      } while (numRemovedInIteration > 0);
      // then, replace f = 2
      this->removeTwofunctionalCrosslinks(net, displacements, springPartitions);
      // and remove all springs that are inactive
      size_t numSpringsRemoved = 0;
      for (long int springIdx = net.nrOfSprings - 1; springIdx >= 0;
           --springIdx) {
        Eigen::Vector3d distance =
          (net.coordinates.segment(3 * net.springIndexA[springIdx], 3) +
           displacements.segment(3 * net.springIndexA[springIdx], 3)) -
          (net.coordinates.segment(3 * net.springIndexB[springIdx], 3) +
           displacements.segment(3 * net.springIndexB[springIdx], 3));
        this->box.handlePBC(distance);
        if (distance.squaredNorm() < tolerance &&
            net.linkIndicesOfSprings[springIdx].size() <= 2) {
          // remove
          this->removeSpring(net, displacements, springPartitions, springIdx);
          numSpringsRemoved += 1;
        }
      }

      this->validateNetwork(net, displacements, springPartitions);

      if (numSpringsRemoved > 0) {
        numRemovedTotal += this->doRemovalAndreisWay(
          net, displacements, springPartitions, tolerance);
      }
      return numRemovedTotal;
    };

    /**
     * @brief Add a slip-link to a given partial spring
     *
     * @param net
     * @param springPartitions
     * @param partialSpringIdx
     * @param slipLinkIdx
     * @param alpha
     */
    size_t MEHPForceBalance::addSlipLinkToPartialSpring(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& springPartitions,
      const size_t partialSpringIdx,
      const size_t slipLinkIdx,
      const double oneOverSpringPartitionUpperLimit) const
    {
      INVALIDARG_EXP_IFN(
        !net.linkIsSliplink[net.springPartIndexA[partialSpringIdx]] ||
          !net.linkIsSliplink[net.springPartIndexB[partialSpringIdx]],
        "Require at least one part to be a cross-link.");
      net.nrOfPartialSprings += 1;
      const size_t newPartialSpringIdx = net.nrOfPartialSprings - 1;
      const size_t relevantSpring =
        net.partialToFullSpringIndex[partialSpringIdx];
      const double N = net.springsContourLength[relevantSpring];
      const double minAlpha =
        (oneOverSpringPartitionUpperLimit > 0.)
          ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
          : 1e-9;
      INVALIDARG_EXP_IFN(APPROX_WITHIN(minAlpha, 0.0, 1.0, 1e-12),
                         "minAlpha must be within 0. and 1.");
      RUNTIME_EXP_IFN(
        minAlpha * (net.localToGlobalSpringIndex[relevantSpring].size() + 1.) <=
          1.,
        "With this minimum alpha, the slip-link cannot be placed on this "
        "partial spring.");
      // std::cout << "Adding slip-link " << slipLinkIdx << " to spring "
      //           << relevantSpring << " (partial " << partialSpringIdx
      //           << ") with minAlpha = " << minAlpha << std::endl;
      // resize the structures
      springPartitions.conservativeResize(net.nrOfPartialSprings);
      assert(springPartitions.size() == net.nrOfPartialSprings);
      net.springPartIndexA.conservativeResize(net.nrOfPartialSprings);
      net.springPartIndexB.conservativeResize(net.nrOfPartialSprings);
      net.springPartCoordinateIndexA.conservativeResize(3 *
                                                        net.nrOfPartialSprings);
      net.springPartCoordinateIndexB.conservativeResize(3 *
                                                        net.nrOfPartialSprings);
      net.partialToFullSpringIndex.conservativeResize(net.nrOfPartialSprings);
      net.partialSpringIsPartial.conservativeResize(net.nrOfPartialSprings);
      // add the new info
      net.partialSpringIsPartial[partialSpringIdx] = true;
      net.partialSpringIsPartial[newPartialSpringIdx] = true;
      net.partialToFullSpringIndex[newPartialSpringIdx] = relevantSpring;

      pylimer_tools::utils::addIfNotContained(
        net.springIndicesOfLinks[slipLinkIdx], relevantSpring);

      size_t oldPartnerA = net.springPartIndexA[partialSpringIdx];
      size_t oldPartnerB = net.springPartIndexB[partialSpringIdx];

      // slightly change numbering to keep the numbering of
      // localToGlobalSpringIndex constant. I.e., we want the
      // `newPartialSpringIdx` to correspond to the spring with the cross-link
      size_t springPartToReplace;
      bool forward;
      if (net.localToGlobalSpringIndex[relevantSpring][0] == partialSpringIdx) {
        forward = true;
        if (net.linkIndicesOfSprings[relevantSpring][0] == oldPartnerA) {
          springPartToReplace = oldPartnerA;
          // std::cout << "Case 1a" << std::endl;
        } else {
          assert(net.linkIndicesOfSprings[relevantSpring][0] == oldPartnerB);
          springPartToReplace = oldPartnerB;
          // std::cout << "Case 1b" << std::endl;
        }
        net.linkIndicesOfSprings[relevantSpring].insert(
          net.linkIndicesOfSprings[relevantSpring].begin() + 1, slipLinkIdx);
        net.localToGlobalSpringIndex[relevantSpring].insert(
          net.localToGlobalSpringIndex[relevantSpring].begin(),
          newPartialSpringIdx);
      } else {
        forward = false;
        assert(pylimer_tools::utils::last(
                 net.localToGlobalSpringIndex[relevantSpring]) ==
               partialSpringIdx);
        if (pylimer_tools::utils::last(
              net.linkIndicesOfSprings[relevantSpring]) == oldPartnerA) {
          springPartToReplace = oldPartnerA;
          // std::cout << "Case 2a" << std::endl;
        } else {
          assert(pylimer_tools::utils::last(
                   net.linkIndicesOfSprings[relevantSpring]) == oldPartnerB);
          springPartToReplace = oldPartnerB;
          // std::cout << "Case 2b" << std::endl;
        }
        net.linkIndicesOfSprings[relevantSpring].insert(
          net.linkIndicesOfSprings[relevantSpring].begin() +
            (net.linkIndicesOfSprings[relevantSpring].size() - 1),
          slipLinkIdx);
        net.localToGlobalSpringIndex[relevantSpring].push_back(
          newPartialSpringIdx);
      }

      // rewire the springs
      if (net.springPartIndexB[partialSpringIdx] == springPartToReplace) {
        net.springPartIndexB[partialSpringIdx] = slipLinkIdx;
      } else {
        assert(net.springPartIndexA[partialSpringIdx] == springPartToReplace);
        net.springPartIndexA[partialSpringIdx] = slipLinkIdx;
      }
      net.springPartIndexA[newPartialSpringIdx] = slipLinkIdx;
      net.springPartIndexB[newPartialSpringIdx] = springPartToReplace;
      for (size_t dir = 0; dir < 3; ++dir) {
        net.springPartCoordinateIndexA[3 * partialSpringIdx + dir] =
          3 * net.springPartIndexA[partialSpringIdx] + dir;
        net.springPartCoordinateIndexB[3 * partialSpringIdx + dir] =
          3 * net.springPartIndexB[partialSpringIdx] + dir;

        net.springPartCoordinateIndexA[3 * newPartialSpringIdx + dir] =
          3 * net.springPartIndexA[newPartialSpringIdx] + dir;
        net.springPartCoordinateIndexB[3 * newPartialSpringIdx + dir] =
          3 * net.springPartIndexB[newPartialSpringIdx] + dir;

        net.springPartCoordinateIndexB[3 * newPartialSpringIdx + dir] =
          3 * springPartToReplace + dir;
      }

      // renormalize this spring
      // mostly by moving the next slip-link further
      springPartitions[newPartialSpringIdx] = minAlpha;
      double remainingNormalisationOffset = minAlpha;
      if (forward) {
        for (size_t globalPartSpringIndex :
             net.localToGlobalSpringIndex[relevantSpring]) {
          double currAlpha = springPartitions[globalPartSpringIndex];
          if (currAlpha > minAlpha) {
            springPartitions[globalPartSpringIndex] -=
              remainingNormalisationOffset;
            springPartitions[globalPartSpringIndex] =
              std::max(springPartitions[globalPartSpringIndex], minAlpha);
            remainingNormalisationOffset -=
              (currAlpha - springPartitions[globalPartSpringIndex]);
          }
          if (remainingNormalisationOffset <= 0.) {
            break;
          }
        }
      } else {
        for (int i = net.localToGlobalSpringIndex[relevantSpring].size() - 1;
             i >= 0;
             --i) {
          size_t globalPartSpringIndex =
            net.localToGlobalSpringIndex[relevantSpring][i];
          double currAlpha = springPartitions[globalPartSpringIndex];
          if (currAlpha > minAlpha) {
            springPartitions[globalPartSpringIndex] -=
              remainingNormalisationOffset;
            springPartitions[globalPartSpringIndex] =
              std::max(springPartitions[globalPartSpringIndex], minAlpha);
            remainingNormalisationOffset -=
              (currAlpha - springPartitions[globalPartSpringIndex]);
          }
          if (remainingNormalisationOffset <= 0.) {
            break;
          }
        }
      }
      // check that normalisation worked
      double newTotalForNormalization =
        springPartitions(net.localToGlobalSpringIndex[relevantSpring]).sum();
      RUNTIME_EXP_IFN(APPROX_EQUAL(newTotalForNormalization, 1.0, 1e-9), "");

      return newPartialSpringIdx;
    }

    /**
     * @brief Replace the two springs traversinga a two-functional cross-links
     * with a single spring
     *
     * @param net
     * @param displacements
     * @param springPartitions
     */
    size_t MEHPForceBalance::removeTwofunctionalCrosslinks(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions) const
    {
      size_t numRemoved = 0;
      for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
           --crosslinkIdx) {
        if (net.springIndicesOfLinks[crosslinkIdx].size() == 2) {
          std::vector<size_t> springsToMerge =
            net.springIndicesOfLinks[crosslinkIdx];
          assert(springsToMerge.size() == 2);

          // check that it's not a primary loop in any way:
          if (springsToMerge[0] != springsToMerge[1] &&
              (XOR(net.linkIndicesOfSprings[springsToMerge[0]][0] ==
                     crosslinkIdx,
                   pylimer_tools::utils::last(
                     net.linkIndicesOfSprings[springsToMerge[0]]) ==
                     crosslinkIdx)) &&
              (XOR(net.linkIndicesOfSprings[springsToMerge[1]][0] ==
                     crosslinkIdx,
                   pylimer_tools::utils::last(
                     net.linkIndicesOfSprings[springsToMerge[1]]) ==
                     crosslinkIdx))) {
            // std::cout << "Merging springs " << springsToMerge[0] << " and "
            //           << springsToMerge[1] << std::endl;
            // let's remove this
            // TODO: this is inefficient shit, so much data being moved
            this->mergeSprings(net,
                               springPartitions,
                               springsToMerge[0],
                               springsToMerge[1],
                               crosslinkIdx);

            // this->validateNetwork(net, displacements, springPartitions);
            // std::cout << "Removing link " << crosslinkIdx << std::endl;
            this->removeLink(net, displacements, crosslinkIdx);

            // std::cout << "Removed cross-link " << crosslinkIdx << std::endl;

            // this->validateNetwork(net, displacements, springPartitions);
            numRemoved += 1;
          }
          // else: TODO: decide
        }
      }
      this->validateNetwork(net, displacements, springPartitions);
      return numRemoved;
    }

    /**
     * @brief Updates the partition/parametrisation of a spring around one link
     *
     */
    double MEHPForceBalance::updateSpringPartition(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions, /* gives the parametrisation of N */
      Eigen::VectorXd&
        oneOverSpringPartitions, /* gives the parametrisation of N */
      const size_t linkIdx,
      double oneOverSpringPartitionUpperLimit,
      bool allowSlipLinksToPassEachOther) const
    {
      // std::cout << "Updating spring partition " << linkIdx << " of "
      //           << net.nrOfNodes << " / " << net.nrOfLinks << std::endl;

      INVALIDARG_EXP_IFN(linkIdx < net.springIndicesOfLinks.size(),
                         "Link to update needs to be in the list");
      INVALIDARG_EXP_IFN(net.linkIsSliplink[linkIdx],
                         "Only slip-links may slip along a spring, link " +
                           std::to_string(linkIdx) +
                           " is not one. Network has " +
                           std::to_string(net.nrOfNodes) + " cross- of " +
                           std::to_string(net.nrOfLinks) + " links.");
      std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
      assert(springIndices.size() == 1 || springIndices.size() == 2);
      if (springIndices.size() == 2 && springIndices[0] == springIndices[1]) {
        springIndices.pop_back();
      }
      double residualNorm = 0.0;
      int residualNormContributions = 0;
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
              net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
            size_t neighbourSpringGlobalIdx =
              net.localToGlobalSpringIndex[springIndex][partner_idx];
            double currentS = springPartitions[currentSpringGlobalIdx];
            double nextS = springPartitions[neighbourSpringGlobalIdx];
            const double N = net.springsContourLength[springIndex];
            const double l = (currentS + nextS);
            if (oneOverSpringPartitionUpperLimit > 0.) {
              // TODO: sketch theory why this should/not be necessary!!!
              const double limit =
                std::clamp(1. / (oneOverSpringPartitionUpperLimit *
                                 (nextS + currentS) * (N)),
                           0.,
                           1.);
              idealValue = std::clamp(idealValue, limit, 1. - limit);
              // double oneOverCurrent = 1. / (currentS * N);
              // double oneOverNext = 1. / (nextS * N);
              // double limitedOneOverCurrent =
              // CLAMP_ONE_OVER_SPRINGPARTITION(
              //   true, oneOverCurrent, N, oneOverSpringPartitionUpperLimit);
              // double limitedOneOverNext = CLAMP_ONE_OVER_SPRINGPARTITION(
              //   true, oneOverNext, N, oneOverSpringPartitionUpperLimit);
              // currentS = (1./N) * 1. / limitedOneOverCurrent;
              // nextS = (1./N) * 1. / limitedOneOverNext;
            }
            double newS = idealValue * l;
            double idealValueM1 = (1. - idealValue);
            double complementaryS = (1. - idealValue) * l;
            double localResidualNorm = 0.0;
            residualNormContributions += 2;
            if (idealValue > 0.0 && idealValue < 1.0) {
              double idealValue2 = idealValue * idealValue;
              double idealValueM12 = idealValueM1 * idealValueM1;
              // way too complicated expression to solve subtraction
              // truncation issues?
              localResidualNorm =
                (std::fma(distanceBack,
                          idealValueM12,
                          -1. * distanceForward * idealValue2)) /
                (idealValueM12 * idealValue2);
            }
            // if ((1. - idealValue) != 0. && l > 0.) {
            //   double term1 = -
            //     (distanceForward / ((1. - idealValue) * (1. -
            //     idealValue))); localResidualNorm += term1;
            //     std::cout.precision(std::numeric_limits<double>::max_digits10);
            //     std::cout << "localResidualNorm term 1: " << term1 <<
            //     std::endl;

            // } else {
            //   std::cout << "localResidualNorm Case 1: " << l << " "
            //             << idealValue << std::endl;
            // }
            // if (idealValue != 0. && l > 0.) {
            //   double term2 =
            //     (distanceBack / (idealValue * idealValue));
            //     localResidualNorm += term2;
            //     std::cout << "localResidualNorm term 2: " << term2 <<
            //     std::endl;

            // } else {
            //   std::cout << "localResidualNorm Case 2: " << l << " "
            //             << idealValue << std::endl;
            // }
            localResidualNorm /= (N * l);
            // std::cout
            //   << "localResidualNorm val 2: "
            //   << localResidualNorm
            //                  << std::endl;

            RUNTIME_EXP_IFN(
              APPROX_WITHIN(newS + complementaryS, 0., 1., 1e-9),
              "Require newS + complementaryS to be within 0, 1, got " +
                std::to_string(newS + complementaryS) + " from " +
                std::to_string(newS) + " and " +
                std::to_string(complementaryS) +
                " with ideal = " + std::to_string(idealValue) + " of " +
                std::to_string(nextS + currentS) + " for link " +
                std::to_string(linkIdx) + ". Diff: " +
                std::to_string(1. - (newS + complementaryS)) + ".");
            RUNTIME_EXP_IFN(
              APPROX_EQUAL(nextS + currentS, newS + complementaryS, 1e-9),
              "Require nextS + currentS == newS + complementaryS, got " +
                std::to_string(nextS + currentS) + " vs. " +
                std::to_string(newS + complementaryS) + " from " +
                std::to_string(nextS) + " and " + std::to_string(currentS) +
                ", " + std::to_string(newS) + " and " +
                std::to_string(complementaryS) + ". Diff: " +
                std::to_string((nextS + currentS) - (newS + complementaryS)) +
                ".");
            RUNTIME_EXP_IFN(
              APPROX_WITHIN(nextS + currentS, 0., 1., 1e-9),
              "Require nextS + currentS to be within 0, 1, got " +
                std::to_string(nextS + currentS) + " from " +
                std::to_string(nextS) + " and " + std::to_string(currentS) +
                ". Diff: " + std::to_string(1. - (nextS + currentS)) + ".");
            RUNTIME_EXP_IFN(
              nextS >= 0.0,
              "nextS must be >= 0., got " + std::to_string(nextS) + " from " +
                std::to_string(nextS) + " and " + std::to_string(currentS) +
                ", " + std::to_string(newS) + " and " +
                std::to_string(complementaryS) + ".");
            RUNTIME_EXP_IFN(complementaryS >= 0.0,
                            "complementaryS must be >= 0., got " +
                              std::to_string(complementaryS) + " from " +
                              std::to_string(nextS) + " and " +
                              std::to_string(currentS) + ", " +
                              std::to_string(newS) + " and " +
                              std::to_string(complementaryS) + ".");

            // (complementaryS > residualNormSTolerance &&
            //  newS > residualNormSTolerance)
            //   ? ( -
            //      distanceBack / (newS * newS))
            //   : 0.0;
            // if (!(APPROX_EQUAL(newS, currentS, 0.2))) {
            //   std::cout << "Updating " << linkIdx << " to " << newS << "
            //   and
            //   "
            //             << complementaryS << " with global springs "
            //             << currentSpringGlobalIdx << " and "
            //             << neighbourSpringGlobalIdx << " from " << currentS
            //             << ", " << nextS << std::endl;
            // }
            // std::cout
            //   << "Contribution to "
            //   << linkIdx
            //           << " from global springs " <<
            //           currentSpringGlobalIdx
            //           << " (" << springsPartners[partner_idx - 1] <<
            //           ") "
            //           << vecBack[0] << ", " << vecBack[1] << ", " <<
            //           vecBack[2]
            //           << " and " << neighbourSpringGlobalIdx << " ("
            //           << springsPartners[partner_idx + 1] << ") "
            //           << vecForward[0] << ", " << vecForward[1] << ",
            //           "
            //           << vecForward[2] << "; "
            //           << " with " << currentS << ", " << nextS <<
            //           std::endl;
            //        std::cout
            // << "Distances are " << distanceForward
            // << ", "
            //           << distanceBack << " to get ideal value " <<
            //           idealValue
            //           << " for " << (nextS) << " , " << currentS <<
            //           std::endl;
            residualNorm += localResidualNorm * localResidualNorm;
            springPartitions[currentSpringGlobalIdx] = newS;
            springPartitions[neighbourSpringGlobalIdx] = complementaryS;
            if (oneOverSpringPartitions.size() > 0) {
              double primaryCorrectionMultiplierC = static_cast<double>(
                net.springPartIndexA[currentSpringGlobalIdx] !=
                net.springPartIndexB[currentSpringGlobalIdx]);
              double oneOverCurrent =
                primaryCorrectionMultiplierC *
                CLAMP_ONE_OVER_SPRINGPARTITION(
                  net.partialSpringIsPartial[currentSpringGlobalIdx],
                  1.0 / (newS * N),
                  N,
                  oneOverSpringPartitionUpperLimit);
              oneOverSpringPartitions.segment(3 * currentSpringGlobalIdx, 3) =
                Eigen::Vector3d::Constant(oneOverCurrent);
              double primaryCorrectionMultiplierN = static_cast<double>(
                net.springPartIndexA[neighbourSpringGlobalIdx] !=
                net.springPartIndexB[neighbourSpringGlobalIdx]);
              double oneOverNeighbour =
                primaryCorrectionMultiplierN *
                CLAMP_ONE_OVER_SPRINGPARTITION(
                  net.partialSpringIsPartial[neighbourSpringGlobalIdx],
                  1.0 / (complementaryS * N),
                  N,
                  oneOverSpringPartitionUpperLimit);
              oneOverSpringPartitions.segment(3 * neighbourSpringGlobalIdx, 3) =
                Eigen::Vector3d::Constant(oneOverNeighbour);
            }
          }
        }
      }
      assert(residualNormContributions == 4);
      return residualNorm;
    }

    /**
     * @brief Loop all slip-links and move them if appropriate to other springs
     *
     * @param net
     * @param u
     * @param springPartitions
     * @param oneOverSpringPartitionUpperLimit
     */
    void MEHPForceBalance::moveSlipLinksToTheirBestBranch(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      const double oneOverSpringPartitionUpperLimit,
      const int nrOfCrosslinkSwapsAllowedPerSliplink,
      const bool respectLoops)
    {
      for (size_t sliplinkIdx = net.nrOfNodes; sliplinkIdx < net.nrOfLinks;
           ++sliplinkIdx) {
        // check this slip-link
        // std::cout << "Moving slip-link " << sliplinkIdx << " to its best
        // branch"
        //           << std::endl;
        this->moveSlipLinkToItsBestBranch(net,
                                          u,
                                          springPartitions,
                                          sliplinkIdx,
                                          oneOverSpringPartitionUpperLimit,
                                          nrOfCrosslinkSwapsAllowedPerSliplink,
                                          respectLoops);
        // this->validateNetwork(net, u, springPartitions);
      }
      this->validateNetwork(net, u, springPartitions);
    }

    /**
     * @brief Move a slip-link if appropriate to other springs
     *
     * @param net
     * @param u
     * @param springPartitions
     * @param oneOverSpringPartitionUpperLimit
     */
    void MEHPForceBalance::moveSlipLinkToItsBestBranch(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      size_t slipLinkIdx,
      const double oneOverSpringPartitionUpperLimit,
      const int nrOfCrosslinkSwapsAllowedPerSliplink,
      const bool respectLoops)
    {
      INVALIDARG_EXP_IFN(net.linkIsSliplink[slipLinkIdx],
                         "Passed slip-link must be one.");
      std::vector<size_t> associatedSprings =
        net.springIndicesOfLinks[slipLinkIdx];
      // skip slip-links that are with its own spring, for now.
      if (associatedSprings.size() <= 1) {
        return;
      }

      for (size_t springIdx : associatedSprings) {
        const double N = net.springsContourLength[springIdx];
        const double swappableCutoff =
          (oneOverSpringPartitionUpperLimit > 0.)
            ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
            : 1e-12;
        std::vector<size_t> linksOnSpring = net.linkIndicesOfSprings[springIdx];
        for (size_t linkI = 1; linkI < linksOnSpring.size() - 1; linkI++) {
          if (linksOnSpring[linkI] == slipLinkIdx) {
            // found index of this slip-link.
            double partitionBeforeIdx =
              net.localToGlobalSpringIndex[springIdx][linkI - 1];
            double partitionAfterIdx =
              net.localToGlobalSpringIndex[springIdx][linkI];
            double didSwap = false;
            // check whether swap is needed in either direction
            // swap if yes
            if (springPartitions[partitionBeforeIdx] <= swappableCutoff) {
              didSwap = this->swapSlipLinkReversibly(
                net,
                u,
                springPartitions,
                partitionBeforeIdx,
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                respectLoops);
            }
            if (springPartitions[partitionAfterIdx] <= swappableCutoff &&
                !didSwap) {
              didSwap = this->swapSlipLinkReversibly(
                net,
                u,
                springPartitions,
                partitionAfterIdx,
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                respectLoops);
            }
          }
        }
      }
    }

    /**
     * @brief
     *
     * @param net
     * @param u
     * @param springPartitions
     * @param partialSpringIdx
     * @param oneOverSpringPartitionUpperLimit
     * @return bool
     */
    bool MEHPForceBalance::swapSlipLinkReversibly(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      const size_t partialSpringIdx,
      const double oneOverSpringPartitionUpperLimit,
      const int nrOfCrosslinkSwapsAllowedPerSliplink,
      const bool respectLoops)
    {
      INVALIDARG_EXP_IFN(partialSpringIdx < net.nrOfPartialSprings,
                         "Partial spring index out of range: got " +
                           std::to_string(partialSpringIdx) + " for " +
                           std::to_string(net.nrOfPartialSprings) +
                           " partial springs.");
      size_t partnerA = net.springPartIndexA[partialSpringIdx];
      size_t partnerB = net.springPartIndexB[partialSpringIdx];
      INVALIDARG_EXP_IFN(net.linkIsSliplink[partnerA] ||
                           net.linkIsSliplink[partnerB],
                         "Cannot swap cross-link with cross-link.");
      if (partnerA == partnerB) {
        return false;
      }
      size_t fullSpringIdx = net.partialToFullSpringIndex[partialSpringIdx];
      // analyse spring
      bool involvesCrosslink =
        !(net.linkIsSliplink[partnerA] && net.linkIsSliplink[partnerB]);
      if (involvesCrosslink) {
        // first check if allowed.
        size_t indexOfSliplink =
          net.linkIsSliplink[partnerA] ? partnerA : partnerB;
        if ((nrOfCrosslinkSwapsAllowedPerSliplink < 0) ||
            (net.nrOfCrosslinkSwapsEndured[indexOfSliplink - net.nrOfNodes] <
             nrOfCrosslinkSwapsAllowedPerSliplink)) {
          bool didSwap = this->swapSlipLinkWithXlinkReversibly(
            net,
            u,
            springPartitions,
            partialSpringIdx,
            oneOverSpringPartitionUpperLimit,
            respectLoops);
          if (didSwap) {
            net.nrOfCrosslinkSwapsEndured[indexOfSliplink - net.nrOfNodes] += 1;
          }
          return didSwap;
        }
        return false;
      } else {
        return this->swapSlipLinksReversibly(net,
                                             u,
                                             springPartitions,
                                             partialSpringIdx,
                                             oneOverSpringPartitionUpperLimit);
      }
    }

    /**
     * @brief
     *
     * @param net
     * @param u
     * @param springPartitions
     * @param partialSpringIdx
     * @param oneOverSpringPartitionUpperLimit
     * @return bool
     */
    bool MEHPForceBalance::swapSlipLinkWithXlinkReversibly(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      const size_t partialSpringIdx,
      const double oneOverSpringPartitionUpperLimit,
      const bool respectLoops)
    {
      size_t partnerA = net.springPartIndexA[partialSpringIdx];
      size_t partnerB = net.springPartIndexB[partialSpringIdx];
      INVALIDARG_EXP_IFN(
        XOR(net.linkIsSliplink[partnerA], net.linkIsSliplink[partnerB]),
        "This method only swaps cross-link with slip-link.");
      if (partnerA == partnerB) {
        return false;
      }
      size_t fullSpringIdx = net.partialToFullSpringIndex[partialSpringIdx];
      // analyse spring
      size_t crosslinkIdx = net.linkIsSliplink[partnerA] ? partnerB : partnerA;
      size_t slipLinkIdx = net.linkIsSliplink[partnerB] ? partnerB : partnerA;
      // compute the residual
      // TODO: check if this is dangerous due to the differences being hidden in
      // the truncated digits
      std::vector<size_t> relevantNeighboursA =
        this->getNeighbourLinkIndices(net, partnerA);
      std::vector<size_t> relevantNeighboursB =
        this->getNeighbourLinkIndices(net, partnerB);
      // combine these
      std::vector<size_t> relevantNeighbours;
      relevantNeighbours.reserve(relevantNeighboursA.size() +
                                 relevantNeighboursB.size() +
                                 2); // preallocate memory
      std::vector<size_t> relevantPartialSprings;
      relevantPartialSprings.reserve(
        net.localToGlobalSpringIndex[fullSpringIdx].size());
      relevantNeighbours.insert(relevantNeighbours.end(),
                                relevantNeighboursA.begin(),
                                relevantNeighboursA.end());
      relevantNeighbours.insert(relevantNeighbours.end(),
                                relevantNeighboursB.begin(),
                                relevantNeighboursB.end());
      relevantNeighbours.push_back(partnerA);
      relevantNeighbours.push_back(partnerB);
      relevantPartialSprings.insert(
        relevantPartialSprings.end(),
        net.localToGlobalSpringIndex[fullSpringIdx].begin(),
        net.localToGlobalSpringIndex[fullSpringIdx].end());
      // for cross-links, we need to take all partners of all springs into
      // account
      for (size_t crosslinksSpringIdx :
           net.springIndicesOfLinks[crosslinkIdx]) {
        relevantNeighbours.insert(
          relevantNeighbours.end(),
          net.linkIndicesOfSprings[crosslinksSpringIdx].begin(),
          net.linkIndicesOfSprings[crosslinksSpringIdx].end());
        relevantPartialSprings.insert(
          relevantPartialSprings.end(),
          net.localToGlobalSpringIndex[crosslinksSpringIdx].begin(),
          net.localToGlobalSpringIndex[crosslinksSpringIdx].end());
      }
      std::vector<size_t> relevantNeighboursCoordIndices;
      relevantNeighboursCoordIndices.reserve(3 * relevantNeighbours.size());
      for (size_t relevantNeighbour : relevantNeighbours) {
        relevantNeighboursCoordIndices.push_back(relevantNeighbour * 3 + 0);
        relevantNeighboursCoordIndices.push_back(relevantNeighbour * 3 + 1);
        relevantNeighboursCoordIndices.push_back(relevantNeighbour * 3 + 2);
      }
      // maybe remove duplicates – might be unnecessary?
      // compute the residual before the deformation
      const double residualBefore =
        this
          ->evaluateStressTensorForLinks(relevantNeighbours,
                                         net,
                                         u,
                                         springPartitions,
                                         1.0,
                                         oneOverSpringPartitionUpperLimit)
          .diagonal()
          .squaredNorm();

      // remember the current positions and partitions
      Eigen::VectorXd displacementsBefore = u(relevantNeighboursCoordIndices);
      Eigen::VectorXd springPartitionsBefore =
        springPartitions(relevantPartialSprings);

      // do swap
      size_t newPartialSpringIdx = 0;
      // swap with cross-link
      std::vector<size_t> springsOfCrosslink =
        net.springIndicesOfLinks[crosslinkIdx];
      if (springsOfCrosslink.size() < 2) {
        return false;
      }
      newPartialSpringIdx =
        this->rotateSlipLinkAroundCrosslink(net,
                                            u,
                                            springPartitions,
                                            partialSpringIdx,
                                            oneOverSpringPartitionUpperLimit,
                                            respectLoops);
      if (newPartialSpringIdx < 0) {
        return false;
      }

      // relax the affected links
      for (size_t relaxSteps = 0; relaxSteps < 2; ++relaxSteps) {
        this->relaxationLight(
          net, springPartitions, u, partnerA, oneOverSpringPartitionUpperLimit);
        this->relaxationLight(
          net, springPartitions, u, partnerB, oneOverSpringPartitionUpperLimit);
      }

      // compute if the residual is lower now
      double residualAfter =
        this
          ->evaluateStressTensorForLinks(relevantNeighbours,
                                         net,
                                         u,
                                         springPartitions,
                                         1.0,
                                         oneOverSpringPartitionUpperLimit)
          .diagonal()
          .squaredNorm();

      if (residualAfter <= residualBefore) {
        return true;
      }

      // otherwise, swap back
      // rotate back to the first spring
      size_t rotations = 0;
      bool isBackToInitialSpring = false;
      while (residualAfter > residualBefore && !isBackToInitialSpring &&
             rotations < 5 && newPartialSpringIdx >= 0) {
        newPartialSpringIdx =
          this->rotateSlipLinkAroundCrosslink(net,
                                              u,
                                              springPartitions,
                                              newPartialSpringIdx,
                                              oneOverSpringPartitionUpperLimit,
                                              respectLoops);
        isBackToInitialSpring = pylimer_tools::utils::contains(
          net.springIndicesOfLinks[slipLinkIdx], fullSpringIdx);
        for (size_t relaxSteps = 0; relaxSteps < 2; ++relaxSteps) {
          // TODO: this is not good at all.
          for (size_t linkIdx : relevantNeighbours) {
            this->relaxationLight(net,
                                  springPartitions,
                                  u,
                                  linkIdx,
                                  oneOverSpringPartitionUpperLimit);
          }
        }
        // compute if the residual is lower now
        residualAfter =
          this
            ->evaluateStressTensorForLinks(relevantNeighbours,
                                           net,
                                           u,
                                           springPartitions,
                                           1.0,
                                           oneOverSpringPartitionUpperLimit)
            .diagonal()
            .squaredNorm();
        rotations += 1;
      }
      if (rotations >= 5) {
        std::cerr << "Could not rotate slip-link " << slipLinkIdx
                  << " back to initial spring. "
                  << "Initial spring was " << fullSpringIdx
                  << ", whereas current springs are "
                  << pylimer_tools::utils::join(
                       net.springIndicesOfLinks[slipLinkIdx].begin(),
                       net.springIndicesOfLinks[slipLinkIdx].end(),
                       std::string(", "))
                  << ". Cross-link is " << crosslinkIdx
                  << " which is associated with springs "
                  << pylimer_tools::utils::join(
                       net.springIndicesOfLinks[crosslinkIdx].begin(),
                       net.springIndicesOfLinks[crosslinkIdx].end(),
                       std::string(", "))
                  << std::endl;
      }

      // relax the affected links back
      // TODO: this is not nice, but currently required because the numbers
      // change
      for (size_t relaxSteps = 0; relaxSteps < 2; ++relaxSteps) {
        // TODO: this is not good at all.
        for (size_t linkIdx : relevantNeighbours) {
          this->relaxationLight(net,
                                springPartitions,
                                u,
                                linkIdx,
                                oneOverSpringPartitionUpperLimit);
        }
      }

      return !isBackToInitialSpring; //(residualBefore < residualAfter);
    }

    /**
     * @brief
     *
     * @param net
     * @param u
     * @param springPartitions
     * @param partialSpringIdx
     * @param oneOverSpringPartitionUpperLimit
     * @return bool
     */
    bool MEHPForceBalance::swapSlipLinksReversibly(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      const size_t partialSpringIdx,
      const double oneOverSpringPartitionUpperLimit)
    {
      size_t partnerA = net.springPartIndexA[partialSpringIdx];
      size_t partnerB = net.springPartIndexB[partialSpringIdx];
      INVALIDARG_EXP_IFN(net.linkIsSliplink[partnerA] &&
                           net.linkIsSliplink[partnerB],
                         "This method only swaps slip-links.");
      if (partnerA == partnerB) {
        return false;
      }
      size_t fullSpringIdx = net.partialToFullSpringIndex[partialSpringIdx];
      // compute the residual
      // TODO: check if this is dangerous due to the differences being hidden in
      // the truncated digits
      std::vector<size_t> relevantNeighboursA =
        this->getNeighbourLinkIndices(net, partnerA);
      std::vector<size_t> relevantNeighboursB =
        this->getNeighbourLinkIndices(net, partnerB);
      // combine these
      std::vector<size_t> relevantNeighbours;
      relevantNeighbours.reserve(relevantNeighboursA.size() +
                                 relevantNeighboursB.size() +
                                 2); // preallocate memory
      std::vector<size_t> relevantPartialSprings;
      relevantPartialSprings.reserve(
        net.localToGlobalSpringIndex[fullSpringIdx].size());
      relevantNeighbours.insert(relevantNeighbours.end(),
                                relevantNeighboursA.begin(),
                                relevantNeighboursA.end());
      relevantNeighbours.insert(relevantNeighbours.end(),
                                relevantNeighboursB.begin(),
                                relevantNeighboursB.end());
      relevantNeighbours.push_back(partnerA);
      relevantNeighbours.push_back(partnerB);
      relevantPartialSprings.insert(
        relevantPartialSprings.end(),
        net.localToGlobalSpringIndex[fullSpringIdx].begin(),
        net.localToGlobalSpringIndex[fullSpringIdx].end());

      std::vector<size_t> relevantNeighboursCoordIndices;
      relevantNeighboursCoordIndices.reserve(3 * relevantNeighbours.size());
      for (size_t relevantNeighbour : relevantNeighbours) {
        relevantNeighboursCoordIndices.push_back(relevantNeighbour * 3 + 0);
        relevantNeighboursCoordIndices.push_back(relevantNeighbour * 3 + 1);
        relevantNeighboursCoordIndices.push_back(relevantNeighbour * 3 + 2);
      }
      // maybe remove duplicates – might be unnecessary?
      // compute the residual before the deformation
      const double residualBefore =
        this
          ->evaluateStressTensorForLinks(relevantNeighbours,
                                         net,
                                         u,
                                         springPartitions,
                                         1.0,
                                         oneOverSpringPartitionUpperLimit)
          .diagonal()
          .squaredNorm();

      // remember the current positions and partitions
      Eigen::VectorXd displacementsBefore = u(relevantNeighboursCoordIndices);
      Eigen::VectorXd springPartitionsBefore =
        springPartitions(relevantPartialSprings);

      // do swap
      this->swapSlipLinks(net, partialSpringIdx);

      // relax the affected links
      for (size_t relaxSteps = 0; relaxSteps < 2; ++relaxSteps) {
        this->relaxationLight(
          net, springPartitions, u, partnerA, oneOverSpringPartitionUpperLimit);
        this->relaxationLight(
          net, springPartitions, u, partnerB, oneOverSpringPartitionUpperLimit);
      }

      // compute if the residual is lower now
      double residualAfter =
        this
          ->evaluateStressTensorForLinks(relevantNeighbours,
                                         net,
                                         u,
                                         springPartitions,
                                         1.0,
                                         oneOverSpringPartitionUpperLimit)
          .diagonal()
          .squaredNorm();

      if (residualAfter <= residualBefore) {
        return true;
      }

      // otherwise, swap back
      this->swapSlipLinks(net, partialSpringIdx);

      // relax the affected links back
      u(relevantNeighboursCoordIndices) = displacementsBefore;
      springPartitions(relevantPartialSprings) = springPartitionsBefore;

      return false;
    }

    /**
     * @brief Do one displacement step
     *
     * @param net
     * @param springPartitions
     * @param linkIdx
     * @param oneOverSpringPartitionUpperLimit
     */
    void MEHPForceBalance::relaxationLight(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& springPartitions,
      Eigen::VectorXd& oneOverSpringPartitions,
      Eigen::VectorXd& u,
      const size_t linkIdx,
      const double oneOverSpringPartitionUpperLimit)
    {
      if (net.linkIsSliplink[linkIdx]) {
        this->updateSpringPartition(net,
                                    u,
                                    springPartitions,
                                    oneOverSpringPartitions,
                                    linkIdx,
                                    oneOverSpringPartitionUpperLimit);
      }
      this->displaceToMeanPosition(
        net, u, springPartitions, linkIdx, oneOverSpringPartitionUpperLimit);
    }

    /**
     * @brief Loop all springs, swap slip-links on them if they are close
     * enough
     *
     * @param net
     * @param u
     * @param springPartitions
     * @param oneOverSpringPartitionUpperLimit
     */
    void MEHPForceBalance::swapSlipLinksInclXlinks(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      double oneOverSpringPartitionUpperLimit,
      bool respectLoops)
    {
      // this->validateNetwork(net, u, springPartitions);
      for (size_t springIdx = 0; springIdx < net.nrOfSprings; ++springIdx) {
        if (net.linkIndicesOfSprings[springIdx].size() <= 2) {
          // no need to handle springs without slip-links
          continue;
        }

        const double N = net.springsContourLength[springIdx];
        const double swappableCutoff =
          (oneOverSpringPartitionUpperLimit > 0.)
            ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
            : 1e-12;

        // loop the remaining partial springs
        // NOTE: it is slightly problematc, that e.g.
        // net.localToGlobalSpringIndex changes /!\  no idea ye how to easily
        // compensate that...
        for (int partialIdx =
               net.localToGlobalSpringIndex[springIdx].size() - 1;
             partialIdx >= 0;
             --partialIdx) {
          // check if they qualify for swapping
          if (springPartitions[net.localToGlobalSpringIndex[springIdx]
                                                           [partialIdx]] <=
              swappableCutoff) {
            size_t partialSpringIdx =
              net.localToGlobalSpringIndex[springIdx][partialIdx];
            if (net.springPartIndexA[partialSpringIdx] !=
                net.springPartIndexB[partialSpringIdx]) {
              // do the swap
              if (partialIdx == 0 ||
                  partialIdx ==
                    net.localToGlobalSpringIndex[springIdx].size() - 1) {
                // swap with x-link
                this->rotateSlipLinkAroundCrosslink(
                  net,
                  u,
                  springPartitions,
                  partialSpringIdx,
                  oneOverSpringPartitionUpperLimit,
                  respectLoops);
                // this->validateNetwork(net, u, springPartitions);
                // std::cout << "Finished moving link " << involvedSlipLink
                //           << " around cross-link " << involvedCrosslink
                //           << " from partial " << partialSpringIdx << " to "
                //           << targetPartialSpringIdx << std::endl;
              } else {
                this->swapSlipLinks(net, partialSpringIdx);
                // this->validateNetwork(net, u, springPartitions);
              }
            }
          }
        }
      }
      this->validateNetwork(net, u, springPartitions);
    }

    /**
     * @brief
     *
     * @param net
     * @param u
     * @param springPartitions
     * @param oneOverSpringPartitionUpperLimit
     */
    void MEHPForceBalance::swapSlipLinks(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      double oneOverSpringPartitionUpperLimit)
    {
      for (size_t springIdx = 0; springIdx < net.nrOfSprings; ++springIdx) {
        if (net.linkIndicesOfSprings[springIdx].size() <= 3) {
          // no need to handle springs with 1 or less slip-links
          continue;
        }

        const double N = net.springsContourLength[springIdx];
        const double swappableCutoff =
          (oneOverSpringPartitionUpperLimit > 0.)
            ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
            : 1e-12;

        // loop the remaining partial springs
        for (size_t partialIdx = 1;
             partialIdx < net.localToGlobalSpringIndex[springIdx].size() - 1;
             ++partialIdx) {
          // check if they qualify for swapping
          if (springPartitions[net.localToGlobalSpringIndex[springIdx]
                                                           [partialIdx]] <=
              swappableCutoff) {
            size_t partialSpringIdx =
              net.localToGlobalSpringIndex[springIdx][partialIdx];
            if (net.springPartIndexA[partialSpringIdx] !=
                net.springPartIndexB[partialSpringIdx]) {
              // do the swap
              this->swapSlipLinks(net, partialSpringIdx);
            }
          }
        }
      }
      this->validateNetwork(net, u, springPartitions);
    }

    /**
     * @brief Move a slip-link from one spring attached to a cross-link to
     * another spring attached to the same cross-link
     *
     * Returns the idx of the new partial spring that had to be introduced.
     * Returns a negative idx if the move was illegal or impossible.
     *
     * @param net
     * @param u
     * @param springPartitions
     * @param partialSpringIdx
     */
    long int MEHPForceBalance::rotateSlipLinkAroundCrosslink(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      const size_t partialSpringIdx,
      double oneOverSpringPartitionUpperLimit,
      const bool respectLoops)
    {
      INVALIDARG_EXP_IFN(net.springPartIndexA[partialSpringIdx] !=
                           net.springPartIndexB[partialSpringIdx],
                         "One of the two ends of the partial spring must be a "
                         "slip-link, one a cross-link");
      INVALIDARG_EXP_IFN(net.springPartIndexA[partialSpringIdx] !=
                           net.springPartIndexB[partialSpringIdx],
                         "Cannot rotate");
      // assemble required data
      size_t springIdx = net.partialToFullSpringIndex[partialSpringIdx];
      RUNTIME_EXP_IFN(
        net.localToGlobalSpringIndex[springIdx][0] == partialSpringIdx ||
          pylimer_tools::utils::last(net.localToGlobalSpringIndex[springIdx]) ==
            partialSpringIdx,
        "Partial spring assembly is not correct");
      const double N = net.springsContourLength[springIdx];
      const double swappableCutoff =
        (oneOverSpringPartitionUpperLimit > 0.)
          ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
          : 1e-9;
      const int maxNrOfSliplinksOnSpring = 1.0 / swappableCutoff;
      size_t partialIdx =
        net.localToGlobalSpringIndex[springIdx][0] == partialSpringIdx
          ? 0
          : net.localToGlobalSpringIndex[springIdx].size() - 1;
      // decide on the involved parties
      size_t otherInvolvedPartialSpring =
        partialIdx == 0
          ? net.localToGlobalSpringIndex[springIdx][1]
          : net.localToGlobalSpringIndex
              [springIdx][net.localToGlobalSpringIndex[springIdx].size() - 2];
      size_t involvedSlipLink =
        net.linkIndicesOfSprings
          [springIdx]
          [partialIdx == 0 ? 1
                           : net.linkIndicesOfSprings[springIdx].size() - 2];
      size_t involvedCrosslink =
        net.linkIndicesOfSprings
          [springIdx]
          [partialIdx == 0 ? 0
                           : net.linkIndicesOfSprings[springIdx].size() - 1];
      std::vector<size_t> possibleTargetPartialSprings;
      // find possible target partial springs
      std::vector<size_t> possibleTargetSprings =
        net.springIndicesOfLinks[involvedCrosslink];
      if (respectLoops) {
        // filter out the target springs that may not be a target based on the
        // involved loops
        possibleTargetPartialSprings.erase(
          std::remove_if(possibleTargetPartialSprings.begin(),
                         possibleTargetPartialSprings.end(),
                         [&net, involvedSlipLink](size_t springIdxToCheck) {
                           for (size_t loopIdx :
                                net.loopsOfSliplink[involvedSlipLink]) {
                             for (size_t springIdx : net.loops[loopIdx]) {
                               if (springIdx == springIdxToCheck) {
                                 return true;
                               }
                             }
                           }

                           return false;
                         }),
          possibleTargetPartialSprings.end());
      }
      if (possibleTargetSprings.size() <= 1) {
        // e.g. in the case of many loops :P
        // std::cerr << "Spring " << springIdx << "'s cross-link " <<
        // involvedCrosslink << " has too few attached springs to
        // reasonably make swaps." << std::endl;
        return -1;
      }
      possibleTargetPartialSprings.reserve(possibleTargetSprings.size());
      int currentPartialSpringTargetIdx = -1;
      for (size_t i = 0; i < possibleTargetSprings.size(); ++i) {
        assert(net.linkIndicesOfSprings[possibleTargetSprings[i]][0] ==
                 involvedCrosslink ||
               pylimer_tools::utils::last(
                 net.linkIndicesOfSprings[possibleTargetSprings[i]]) ==
                 involvedCrosslink);
        size_t currentPossibleTargetPartialSpringIdx =
          net.linkIndicesOfSprings[possibleTargetSprings[i]][0] ==
              involvedCrosslink
            ? net.localToGlobalSpringIndex[possibleTargetSprings[i]][0]
            : pylimer_tools::utils::last(
                net.localToGlobalSpringIndex[possibleTargetSprings[i]]);
        if (currentPossibleTargetPartialSpringIdx != partialSpringIdx &&
            net.localToGlobalSpringIndex[possibleTargetSprings[i]].size() <
              maxNrOfSliplinksOnSpring) { // let's not combine with the
                                          // to-be-removed partial spring
          possibleTargetPartialSprings.push_back(
            currentPossibleTargetPartialSpringIdx);
        }
        if (currentPossibleTargetPartialSpringIdx == partialSpringIdx ||
            currentPossibleTargetPartialSpringIdx ==
              otherInvolvedPartialSpring) {
          // we unfortunately cannot assert this due to primary loops
          // assert(currentPartialSpringTargetIdx == -1);
          currentPartialSpringTargetIdx = i;
        }
      }
      if (possibleTargetPartialSprings.size() == 0) {
        return -1;
      }
      // remove the slip-link from one branch of the x-link
      // but skip resizing the Eigen structures, since the additional rows are
      // still needed
      this->mergePartialSprings(net,
                                springPartitions,
                                partialSpringIdx,
                                otherInvolvedPartialSpring,
                                involvedSlipLink,
                                true);
      // this->validateNetwork(net, u, springPartitions);
      // ... and add it to another
      // assert(currentPartialSpringTargetIdx >= 0);
      size_t targetPartialSpringIdx =
        possibleTargetPartialSprings[(currentPartialSpringTargetIdx) %
                                     possibleTargetPartialSprings.size()];
      if (targetPartialSpringIdx > partialSpringIdx) {
        targetPartialSpringIdx -= 1;
      }
      // std::cout << "Handling moving link " << involvedSlipLink
      //           << " around cross-link " << involvedCrosslink
      //           << " from partial " << partialSpringIdx << " to "
      //           << targetPartialSpringIdx << std::endl;
      assert(net.springPartIndexA[targetPartialSpringIdx] ==
               involvedCrosslink ||
             net.springPartIndexB[targetPartialSpringIdx] == involvedCrosslink);
      size_t newPartialSpringIdx =
        this->addSlipLinkToPartialSpring(net,
                                         springPartitions,
                                         targetPartialSpringIdx,
                                         involvedSlipLink,
                                         oneOverSpringPartitionUpperLimit);

      if ((net.springPartIndexA[targetPartialSpringIdx] == involvedCrosslink &&
           net.springPartIndexB[targetPartialSpringIdx] == involvedSlipLink) ||
          (net.springPartIndexB[targetPartialSpringIdx] == involvedCrosslink &&
           net.springPartIndexA[targetPartialSpringIdx] == involvedSlipLink)) {
        return targetPartialSpringIdx;
      } else {
        RUNTIME_EXP_IFN(
          (net.springPartIndexA[newPartialSpringIdx] == involvedCrosslink &&
           net.springPartIndexB[newPartialSpringIdx] == involvedSlipLink) ||
            (net.springPartIndexB[newPartialSpringIdx] == involvedCrosslink &&
             net.springPartIndexA[newPartialSpringIdx] == involvedSlipLink),
          "Expected to find cross- and slip-link at either partial spring, but "
          "did not.");
        return newPartialSpringIdx;
      }
    }

    /**
     * @brief
     *
     * @param net
     * @param partialSpringIdx
     */
    void MEHPForceBalance::swapSlipLinks(ForceBalanceNetwork& net,
                                         const size_t partialSpringIdx)
    {
      const size_t linkIdx1 = net.springPartIndexA[partialSpringIdx];
      const size_t linkIdx2 = net.springPartIndexB[partialSpringIdx];
      INVALIDARG_EXP_IFN(linkIdx1 != linkIdx2, "Cannot swap link with itself.");
      // std::cout << "Swapping link " << linkIdx1 << " and " << linkIdx2
      //           << std::endl;
      INVALIDARG_EXP_IFN(
        net.linkIsSliplink[linkIdx1],
        "Only partial springs with only slip-links allow swapping.");
      INVALIDARG_EXP_IFN(
        net.linkIsSliplink[linkIdx2],
        "Only partial springs with only slip-links allow swapping.");
      const size_t springIdx = net.partialToFullSpringIndex[partialSpringIdx];
      // find the rest of the connectivity required for swapping
      long int otherPartialOfLinkIdx1 = -1;
      long int otherPartialOfLinkIdx2 = -1;
      long int firstPositionInSpring = -1;
      for (size_t inSpringIdx = 1;
           inSpringIdx < net.linkIndicesOfSprings[springIdx].size() - 1;
           ++inSpringIdx) {
        if (net.localToGlobalSpringIndex[springIdx][inSpringIdx] ==
            partialSpringIdx) {
          if (net.linkIndicesOfSprings[springIdx][inSpringIdx] == linkIdx1 &&
              net.linkIndicesOfSprings[springIdx][inSpringIdx + 1] ==
                linkIdx2) {
            RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 == -1,
                            "Expect to find sequence of links only once.");
            otherPartialOfLinkIdx1 =
              net.localToGlobalSpringIndex[springIdx][inSpringIdx - 1];
            otherPartialOfLinkIdx2 =
              net.localToGlobalSpringIndex[springIdx][inSpringIdx + 1];
            firstPositionInSpring = inSpringIdx;
          }
          // else
          if (net.linkIndicesOfSprings[springIdx][inSpringIdx] == linkIdx2 &&
              net.linkIndicesOfSprings[springIdx][inSpringIdx + 1] ==
                linkIdx1) {
            RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 == -1,
                            "Expect to find sequence of links only once.");
            otherPartialOfLinkIdx2 =
              net.localToGlobalSpringIndex[springIdx][inSpringIdx - 1];
            otherPartialOfLinkIdx1 =
              net.localToGlobalSpringIndex[springIdx][inSpringIdx + 1];
            firstPositionInSpring = inSpringIdx;
          }
        }
      }

      RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 >= 0,
                      "Did not find partial spring " +
                        std::to_string(partialSpringIdx) + " in spring " +
                        std::to_string(springIdx) + ".");
      RUNTIME_EXP_IFN(otherPartialOfLinkIdx2 >= 0,
                      "Did not find partial spring " +
                        std::to_string(partialSpringIdx) + " in spring " +
                        std::to_string(springIdx) + ".");
      RUNTIME_EXP_IFN(firstPositionInSpring >= 0,
                      "Did not find partial spring in spring.");
      RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 != otherPartialOfLinkIdx2,
                      "Required assumption not met.");
      RUNTIME_EXP_IFN(firstPositionInSpring <
                        net.linkIndicesOfSprings[springIdx].size() - 1,
                      "Required assumption not met.");
      // actually do the swapping
      // net.springPartIndexA[partialSpringIdx] = linkIdx2;
      // net.springPartIndexB[partialSpringIdx] = linkIdx1;
      if (net.springPartIndexA[otherPartialOfLinkIdx1] == linkIdx1) {
        net.springPartIndexA[otherPartialOfLinkIdx1] = linkIdx2;
        net.springPartCoordinateIndexA.segment(3 * otherPartialOfLinkIdx1, 3) =
          Eigen::ArrayXi::LinSpaced(3, 3 * linkIdx2, 3 * linkIdx2 + 2);
      } else {
        RUNTIME_EXP_IFN(net.springPartIndexB[otherPartialOfLinkIdx1] ==
                          linkIdx1,
                        "Required assumption apparently not met.");
        net.springPartIndexB[otherPartialOfLinkIdx1] = linkIdx2;
        net.springPartCoordinateIndexB.segment(3 * otherPartialOfLinkIdx1, 3) =
          Eigen::ArrayXi::LinSpaced(3, 3 * linkIdx2, 3 * linkIdx2 + 2);
      }
      if (net.springPartIndexA[otherPartialOfLinkIdx2] == linkIdx2) {
        net.springPartIndexA[otherPartialOfLinkIdx2] = linkIdx1;
        net.springPartCoordinateIndexA.segment(3 * otherPartialOfLinkIdx2, 3) =
          Eigen::ArrayXi::LinSpaced(3, 3 * linkIdx1, 3 * linkIdx1 + 2);
      } else {
        RUNTIME_EXP_IFN(net.springPartIndexB[otherPartialOfLinkIdx2] ==
                          linkIdx2,
                        "Required assumption apparently not met.");
        net.springPartIndexB[otherPartialOfLinkIdx2] = linkIdx1;
        net.springPartCoordinateIndexB.segment(3 * otherPartialOfLinkIdx2, 3) =
          Eigen::ArrayXi::LinSpaced(3, 3 * linkIdx1, 3 * linkIdx1 + 2);
      }

      // std::swap(net.linkIndicesOfSprings[springIdx][firstPositionInSpring],
      //           net.linkIndicesOfSprings[springIdx][firstPositionInSpring +
      //           1]);
      if (net.linkIndicesOfSprings[springIdx][firstPositionInSpring] ==
          linkIdx1) {
        net.linkIndicesOfSprings[springIdx][firstPositionInSpring] = linkIdx2;
        net.linkIndicesOfSprings[springIdx][firstPositionInSpring + 1] =
          linkIdx1;
      } else {
        RUNTIME_EXP_IFN(
          net.linkIndicesOfSprings[springIdx][firstPositionInSpring] ==
            linkIdx2,
          "Required assumption apparently not met.");
        net.linkIndicesOfSprings[springIdx][firstPositionInSpring] = linkIdx1;
        net.linkIndicesOfSprings[springIdx][firstPositionInSpring + 1] =
          linkIdx2;
      }
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
      const ForceBalanceNetwork& net,
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
      u.segment(3 * linkIdx, 3) +=
        objectiveDisplacement / (objectiveDisplacementContributors == 0.0
                                   ? 1.0
                                   : objectiveDisplacementContributors);

      double dist = u.segment(3 * linkIdx, 3).squaredNorm();
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
      const ForceBalanceNetwork& net,
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
              net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
            size_t forwardSpringGlobalIdx =
              net.localToGlobalSpringIndex[springIndex][partner_idx];

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
    }

    Eigen::Matrix3d MEHPForceBalance::evaluateForceOnCrossLink(
      const size_t linkIdx,
      const ForceBalanceNetwork& net,
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
            springGlobalIdx = net.localToGlobalSpringIndex[springIndex][0];
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
            springGlobalIdx = pylimer_tools::utils::last<size_t>(
              net.localToGlobalSpringIndex[springIndex]);
          }
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

    Eigen::Vector3d MEHPForceBalance::evaluateDistanceBetween(
      const ForceBalanceNetwork& net,
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
      this->box.handlePBC<Eigen::Vector3d>(distances);

      if (is2D) {
        distances[2] = 0.0;
      }

      return distances;
    }

    Eigen::VectorXd MEHPForceBalance::evaluateSpringDistances(
      const ForceBalanceNetwork& net,
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
      this->box.handlePBC(springDistances);

      // reset for 2D systems
      if (is2D && net.nrOfSprings > 0) {
        // springDistances(Eigen::seq(2, net.nrOfSprings, 3)) =
        //   Eigen::VectorXd::Zero(net.nrOfSprings);
        for (size_t i = 2; i < 3 * net.nrOfSprings; i += 3) {
          springDistances[i] = 0.0;
        }
      }

      return springDistances;
    }

    Eigen::VectorXd MEHPForceBalance::evaluatePartialSpringDistances(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const bool is2D) const
    {
      // first, the distances
      assert(u.size() == net.coordinates.size());

      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd partialDistances =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA));
      this->box.handlePBC(partialDistances);

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

    void MEHPForceBalance::addSlipLinks(
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
      INVALIDARG_EXP_IFN(
        loopsOfSliplinks.size() == 0 ||
          loopsOfSliplinks.size() == additionalLen,
        "You must provide either loops for all new slip-links, or non at all.");
      INVALIDARG_EXP_IFN(
        (loopsOfSliplinks.size() == 0 &&
         this->initialConfig.loopsOfSliplink.size() == 0) ||
          ((loopsOfSliplinks.size() > 0) &&
           this->initialConfig.loopsOfSliplink.size() ==
             (this->initialConfig.nrOfLinks - this->initialConfig.nrOfNodes)),
        "Cannot add slip-links with loops to structure without, or vice "
        "versa.");
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
      this->initialConfig.nrOfCrosslinkSwapsEndured.conservativeResize(
        this->initialConfig.nrOfLinks - this->initialConfig.nrOfNodes);
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
      this->initialConfig.partialToFullSpringIndex.conservativeResize(
        currentNrOfPartialSprings + 2 * additionalLen);
      this->initialConfig.partialSpringIsPartial.conservativeResize(
        currentNrOfPartialSprings + 2 * additionalLen);
      // handle loops if appropriate
      size_t previousNrOfLoops = this->initialConfig.loops.size();
      this->initialConfig.loops.reserve(previousNrOfLoops + loops.size());
      this->initialConfig.loops.insert(
        this->initialConfig.loops.end(), loops.begin(), loops.end());
      if (loopsOfSliplinks.size() > 0) {
        this->initialConfig.loopsOfSliplink.reserve(
          this->initialConfig.nrOfLinks - this->initialConfig.nrOfNodes);
        this->initialConfig.loopsOfSliplink.insert(
          this->initialConfig.loopsOfSliplink.end(),
          loopsOfSliplinks.begin(),
          loopsOfSliplinks.end());
        if (previousNrOfLoops > 0) {
          // adjust the numbering
          for (size_t i = 0; i < loopsOfSliplinks.size(); ++i) {
            for (size_t j = 0; j < loopsOfSliplinks[i].size(); ++j) {
              loopsOfSliplinks[i][j] += previousNrOfLoops;
            }
          }
        }
      }
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
        this->initialConfig
          .nrOfCrosslinkSwapsEndured[currentNrOfLinks + i -
                                     this->initialConfig.nrOfNodes] = 0;
        std::vector<size_t> springIndices{ strandIdx1[i], strandIdx2[i] };
        std::vector<size_t> springIndicesOfLink =
          (strandIdx1[i] == strandIdx2[i])
            ? std::vector<size_t>{ strandIdx1[i] }
            : springIndices;
        this->initialConfig.springIndicesOfLinks.push_back(springIndicesOfLink);
        // add to the springs
        int springIndexIndex = 0;
        for (size_t springIndex : springIndices) {
          std::vector<size_t> springParticipants =
            this->initialConfig.linkIndicesOfSprings[springIndex];
          double alpha = (springIndexIndex == 0) ? alpha1[i] : alpha2[i];
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
          for (size_t j = 0; j < springParticipants.size() - 1; ++j) {
            partitionsStrand.push_back(
              this->currentSpringPartitionsVec
                [this->initialConfig.localToGlobalSpringIndex[springIndex][j]]);
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
            this->initialConfig
              .localToGlobalSpringIndex[springIndex][targetIndexInSpring];
          size_t newSpringIndex =
            currentNrOfPartialSprings + partialSpringsAdded;

          this->initialConfig.partialSpringIsPartial[lastSpringIndex] = true;
          this->initialConfig.partialSpringIsPartial[newSpringIndex] = true;

          this->initialConfig.localToGlobalSpringIndex[springIndex].insert(
            this->initialConfig.localToGlobalSpringIndex[springIndex].begin() +
              targetIndexInSpring + 1,
            newSpringIndex);
          this->initialConfig.partialToFullSpringIndex[newSpringIndex] =
            (springIndex);

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
            APPROX_WITHIN(
              this->currentSpringPartitionsVec[newSpringIndex], 0.0, 1.0, 1e-9),
            "Spring partition must be between 0 and 1, got " +
              std::to_string(this->currentSpringPartitionsVec[newSpringIndex]) +
              ".");
          this->currentSpringPartitionsVec[lastSpringIndex] = alpha;
          RUNTIME_EXP_IFN(
            APPROX_WITHIN(this->currentSpringPartitionsVec[lastSpringIndex],
                          0.0,
                          1.0,
                          1e-9),
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
     * @param loopTol
     * @return std::array<std::array<double, 3>, 3>
     */
    Eigen::Matrix3d MEHPForceBalance::evaluateStressTensorForLinks(
      const std::vector<size_t> linkIndices,
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      const double kappa0,
      const double oneOverSpringPartitionUpperLimit) const
    {

      Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();
      INVALIDARG_EXP_IFN(
        springPartitions.size() == net.springPartIndexA.size(),
        "Spring partitions must have the size of partial springs.");

      double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

      Eigen::VectorXi debugNrSpringsVisited =
        Eigen::VectorXi::Zero(net.nrOfPartialSprings);

      for (size_t linkIdx : linkIndices) {
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
        RUNTIME_EXP_IFN(std::isfinite(force.squaredNorm()),
                        "Got non-finite force contribution to stress tensor: " +
                          std::to_string(force.squaredNorm()) + " at link " +
                          std::to_string(linkIdx) + "!");
        stress += force;
      }

      return halfOverVolume * stress;
    };

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3>
    MEHPForceBalance::evaluateStressTensorLinkBased(
      const ForceBalanceNetwork& net,
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
    std::array<std::array<double, 3>, 3> MEHPForceBalance::evaluateStressTensor(
      const ForceBalanceNetwork& net,
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
      this->box.handlePBC(relevantPartialDistancesA);

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
    bool MEHPForceBalance::ConvertNetwork(ForceBalanceNetwork& net,
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
      net.nrOfCrosslinkSwapsEndured = Eigen::ArrayXi::Zero(0);
      net.oldAtomIds = Eigen::ArrayXi::Zero(net.nrOfLinks);
      net.linkIsSliplink = ArrayXb::Constant(net.nrOfLinks, false);
      net.springIndicesOfLinks.reserve(net.nrOfLinks);
      net.partialToFullSpringIndex = Eigen::ArrayXi(net.nrOfPartialSprings);
      for (size_t i = 0; i < net.nrOfLinks; ++i) {
        net.springIndicesOfLinks.push_back(std::vector<size_t>());
      }
      net.linkIndicesOfSprings.reserve(net.nrOfSprings);
      this->currentSpringPartitionsVec = Eigen::VectorXd::Ones(net.nrOfSprings);
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
      net.partialSpringIsPartial = ArrayXb::Constant(net.nrOfSprings, false);

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

          pylimer_tools::utils::addIfNotContained(
            net.springIndicesOfLinks[nodeIdxFrom], spring_idx);
          if (nodeIdxFrom != nodeIdxTo) {
            pylimer_tools::utils::addIfNotContained(
              net.springIndicesOfLinks[nodeIdxTo], spring_idx);
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
          net.localToGlobalSpringIndex.push_back(zeroMap);
          net.partialToFullSpringIndex[spring_idx] = (spring_idx);

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

    bool MEHPForceBalance::validateNetwork(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
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
      RUNTIME_EXP_IFN(u.size() == net.nrOfLinks * 3,
                      "Invalid size of displacements");
      RUNTIME_EXP_IFN(u.size() == net.coordinates.size(),
                      "Invalid size of displacements or coordinates");
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
        RUNTIME_EXP_IFN(std::isfinite(u[coordI]),
                        "Displacement component " + std::to_string(coordI) +
                          " must be finite, got " + std::to_string(u[coordI]) +
                          ".");
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
  }
}
}
