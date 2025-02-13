#include "MEHPForceBalance.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
// #include "../utils/MemoryUtil.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef NDEBUG
#define DEBUG_REMOVAL
#endif

namespace pylimer_tools {
namespace sim {
  namespace mehp {
#ifndef CLAMP_ONE_OVER_SPRINGPARTITION
/**
 * @brief a macro for doing the clamping in the routines using kappa,
 * to prevent deivision by zero issues / multiplications by infinity
 */
#define CLAMP_ONE_OVER_SPRINGPARTITION(                                        \
  isPartialSpring, val, N, oneOverSpringPartitionUpperLimit)                   \
  ((!isPartialSpring)                                                          \
     ? val                                                                     \
     : std::clamp(val,                                                         \
                  (oneOverSpringPartitionUpperLimit > 0.)                      \
                    ? (1. / (N - 1. / oneOverSpringPartitionUpperLimit))       \
                    : (0.0),                                                   \
                  (oneOverSpringPartitionUpperLimit > 0.)                      \
                    ? (oneOverSpringPartitionUpperLimit)                       \
                    : (N)));
#endif

    /**
     * FORCE RELAXATION
     */
    void MEHPForceBalance::runForceRelaxation(
      long int maxNrOfSteps, // default: 10000
      double xtol,
      const double initialResidualToUse,
      const StructureSimplificationMode simplificationMode,
      const double inactiveRemovalCutoff,
      bool doInnerIterations,
      LinkSwappingMode allowSlipLinksToPassEachOther,
      const int swappingFrequency,
      const double oneOverSpringPartitionUpperLimit,
      const int nrOfCrosslinkSwapsAllowedPerSliplink,
      const bool disableSlipping,
      const std::function<bool()>& shouldInterrupt,
      const std::function<void()>& cleanupInterrupt)
    {
      this->validateNetwork();
      // INVALIDARG_EXP_IFN(
      //   shouldRemoveInactiveCrosslinks == false &&
      //     remove2functionalCrosslinkers == true,
      //   "Removing 2-functional crosslinkers only makes sense when inactive "
      //   "crosslinkers may be removed too, during the procedure.");
      this->simulationHasRun = true;

      INVALIDARG_EXP_IFN(
        inactiveRemovalCutoff > 0.0 ||
          simplificationMode == StructureSimplificationMode::NO_SIMPLIFICATION,
        "Removal cut-off must be positive when simplification is enabled.");

      /* array allocation */
      std::vector<Eigen::ArrayXi> independentVertexSets;
      double maxDistanceMoved = 0.0;
      size_t indexOfMaxDistanceMoved = 0;

      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(this->initialConfig,
                                             this->currentSpringPartitionsVec,
                                             oneOverSpringPartitionUpperLimit);
      const double initialResidual =
        (initialResidualToUse > 0.)
          ? initialResidualToUse
          : this->getDisplacementResidualNormFor(this->initialConfig,
                                                 this->currentDisplacements,
                                                 oneOverSpringPartitions);
      const double minN = this->initialConfig.springsContourLength.minCoeff();
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
      double currentResidual = initialResidual;
      double previousResidual = initialResidual;
      double intermediateResidual = 0.0;
      size_t iterationsDone = 0;
      size_t nRemoved = 0;

      this->prepareAllOutputs();

      // actual loop
      bool wasInterrupted = false;
      do {
        if (allowSlipLinksToPassEachOther != LinkSwappingMode::NO_SWAPPING) {
          if (swappingFrequency > 0 &&
              (iterationsDone % swappingFrequency) == 0) {
            if (allowSlipLinksToPassEachOther ==
                LinkSwappingMode::SLIPLINKS_ONLY) {
              this->swapSlipLinks(this->initialConfig,
                                  this->currentDisplacements,
                                  this->currentSpringPartitionsVec,
                                  oneOverSpringPartitionUpperLimit);
            } else if (allowSlipLinksToPassEachOther == LinkSwappingMode::ALL) {
              this->swapSlipLinksInclXlinks(this->initialConfig,
                                            this->currentDisplacements,
                                            this->currentSpringPartitionsVec,
                                            oneOverSpringPartitionUpperLimit);
            } else if (allowSlipLinksToPassEachOther ==
                         LinkSwappingMode::ALL_MC_TRY ||
                       allowSlipLinksToPassEachOther ==
                         LinkSwappingMode::ALL_MC) {
              this->moveSlipLinksToTheirBestBranch(
                this->initialConfig,
                this->currentDisplacements,
                this->currentSpringPartitionsVec,
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                false,
                allowSlipLinksToPassEachOther == LinkSwappingMode::ALL_MC_TRY);
            } else if (allowSlipLinksToPassEachOther ==
                         LinkSwappingMode::ALL_MC_TRY_CYCLE ||
                       allowSlipLinksToPassEachOther ==
                         LinkSwappingMode::ALL_MC_CYCLE) {
              this->moveSlipLinksToTheirBestBranch(
                this->initialConfig,
                this->currentDisplacements,
                this->currentSpringPartitionsVec,
                oneOverSpringPartitionUpperLimit,
                nrOfCrosslinkSwapsAllowedPerSliplink,
                true,
                allowSlipLinksToPassEachOther ==
                  LinkSwappingMode::ALL_MC_TRY_CYCLE);
            } else {
              throw std::invalid_argument(
                "This swapping mode is currently not supported.");
            }
            oneOverSpringPartitions = this->assembleOneOverSpringPartition(
              this->initialConfig,
              this->currentSpringPartitionsVec,
              oneOverSpringPartitionUpperLimit);
          }
        }
        maxDistanceMoved = 0.0;

        intermediateResidual =
          this->getDisplacementResidualNormFor(this->initialConfig,
                                               this->currentDisplacements,
                                               oneOverSpringPartitions);

        // place slip-link
        for (size_t link_idx = this->initialConfig.nrOfNodes;
             link_idx < this->initialConfig.nrOfLinks;
             ++link_idx) {
          assert(this->initialConfig.linkIsSliplink[link_idx]);
          // std::cout << "Handling " << link_idx << " of " << net.nrOfNodes
          //           << " / " << net.nrOfLinks << std::endl;

          // std::cout << "Still handling " << link_idx << " of " <<
          // net.nrOfNodes
          //           << " / " << net.nrOfLinks << std::endl;
          int innerIterationsDone = 0;
          do {
            if (!disableSlipping) {
              double r2 =
                this->updateSpringPartition(this->initialConfig,
                                            this->currentDisplacements,
                                            this->currentSpringPartitionsVec,
                                            oneOverSpringPartitions,
                                            link_idx,
                                            oneOverSpringPartitionUpperLimit,
                                            allowSlipLinksToPassEachOther);
            }
            double displacementDone =
              this->displaceToMeanPosition(this->initialConfig,
                                           this->currentDisplacements,
                                           this->currentSpringPartitionsVec,
                                           link_idx,
                                           oneOverSpringPartitionUpperLimit);
            maxDistanceMoved = std::max(maxDistanceMoved, displacementDone);
            innerIterationsDone += 1;
          } while (doInnerIterations && innerIterationsDone < 50);
        }

        oneOverSpringPartitions = this->assembleOneOverSpringPartition(
          this->initialConfig,
          this->currentSpringPartitionsVec,
          oneOverSpringPartitionUpperLimit);
        intermediateResidual =
          this->getDisplacementResidualNormFor(this->initialConfig,
                                               this->currentDisplacements,
                                               oneOverSpringPartitions);

        // place crosslinkers
        for (size_t link_idx = 0; link_idx < this->initialConfig.nrOfNodes;
             ++link_idx) {
          assert(!this->initialConfig.linkIsSliplink[link_idx]);
          double distanceMoved =
            this->displaceToMeanPosition(this->initialConfig,
                                         this->currentDisplacements,
                                         this->currentSpringPartitionsVec,
                                         link_idx,
                                         oneOverSpringPartitionUpperLimit);
          if (distanceMoved > maxDistanceMoved) {
            maxDistanceMoved = distanceMoved;
            indexOfMaxDistanceMoved = link_idx;
          }
        }

        currentResidual = this->getDisplacementResidualNormFor(
          this->initialConfig,
          this->currentDisplacements,
          this->currentSpringPartitionsVec,
          oneOverSpringPartitionUpperLimit);
        // if (previousResidual < currentResidual && nRemoved == 0) {
        //   throw std::runtime_error(
        //     "Residual is bigger (" + std::to_string(currentResidual) + " vs.
        //     " + std::to_string(previousResidual) +
        //     ") than in the previous iteration. This makes no "
        //     "sense and hints at a mistake.");
        // }
        previousResidual = currentResidual;
        iterationsDone += 1;
        if (iterationsDone % this->simplificationFrequency == 0) {
          this->breakTooLongSprings(this->initialConfig,
                                    this->currentDisplacements,
                                    this->currentSpringPartitionsVec);

          do {
#ifndef NDEBUG
            this->validateNetwork();
#endif
            nRemoved = 0;
            if (simplificationMode ==
                  StructureSimplificationMode::INACTIVE_ONLY ||
                simplificationMode == StructureSimplificationMode::ALL_TIM) {
#ifdef DEBUG_REMOVAL
              std::cout << "Checking and possibly removing inactive cross-links"
                        << std::endl;
#endif
              nRemoved +=
                this->removeInactiveCrosslinks(this->initialConfig,
                                               this->currentDisplacements,
                                               this->currentSpringPartitionsVec,
                                               inactiveRemovalCutoff);
              this->initialConfig.meanSpringContourLength =
                this->initialConfig.springsContourLength.size() > 0
                  ? this->initialConfig.springsContourLength.mean()
                  : 0.;
            }
            if (simplificationMode == StructureSimplificationMode::X2F_ONLY ||
                simplificationMode == StructureSimplificationMode::ALL_TIM) {
#ifdef DEBUG_REMOVAL
              std::cout
                << "Checking and possibly removing cross-links with f = 2"
                << std::endl;
#endif
              nRemoved += this->removeTwofunctionalCrosslinks(
                this->initialConfig,
                this->currentDisplacements,
                this->currentSpringPartitionsVec);
              this->initialConfig.meanSpringContourLength =
                this->initialConfig.springsContourLength.size() > 0
                  ? this->initialConfig.springsContourLength.mean()
                  : 0.;
            }
            if (simplificationMode == StructureSimplificationMode::ALL_ANDREI) {
#ifdef DEBUG_REMOVAL
              std::cout << "Checking and possibly removing cross-links and "
                           "springs, Andrei's way"
                        << std::endl;
#endif
              nRemoved +=
                this->doRemovalAndreisWay(this->initialConfig,
                                          this->currentDisplacements,
                                          this->currentSpringPartitionsVec,
                                          inactiveRemovalCutoff);
            }
            // cleanup some things
            if (simplificationMode !=
                StructureSimplificationMode::NO_SIMPLIFICATION) {
              this->validateNetwork(this->initialConfig,
                                    this->currentDisplacements,
                                    this->currentSpringPartitionsVec);
              oneOverSpringPartitions = this->assembleOneOverSpringPartition(
                this->initialConfig,
                this->currentSpringPartitionsVec,
                oneOverSpringPartitionUpperLimit);
            }
          } while (nRemoved > 0);
        }
        this->handleOutput(iterationsDone);

        if (shouldInterrupt()) {
          wasInterrupted = true;
          break;
        }
      } while (currentResidual / initialResidual > xtol &&
               iterationsDone < maxNrOfSteps &&
               this->initialConfig.nrOfSprings > 0);

      // finish up
      this->closeAllOutputs();

      // query solution & exit reason
      this->exitReason = (iterationsDone == maxNrOfSteps)
                           ? ExitReason::MAX_STEPS
                           : ExitReason::X_TOLERANCE;
      this->nrOfStepsDone += iterationsDone;
      std::cout << iterationsDone << " steps done. "
                << "Last max distance moved: " << maxDistanceMoved << ". "
                << "Current residual: " << currentResidual << ". "
                << "Initial residual: " << initialResidual << ". " << std::endl;

      assert(this->currentDisplacements.size() ==
             3 * this->initialConfig.nrOfLinks);
      this->validateNetwork();
      this->currentSpringVectors = this->evaluateSpringVectors(
        this->initialConfig, this->currentDisplacements);
      this->currentPartialSpringVectors = this->evaluatePartialSpringVectors(
        this->initialConfig, this->currentDisplacements);
      if (wasInterrupted) {
        this->exitReason = ExitReason::INTERRUPT;
        cleanupInterrupt();
      }
    }

    /**
     * @brief Compute the displacement residual norm for the current
     * configuration
     *
     * @param oneOverSpringPartitionUpperLimit
     * @return double
     */
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

    /**
     * @brief Compute the displacement residual norm for a specific
     * configuration
     *
     * @param net
     * @param u
     * @param oneOverSpringPartitions
     * @return double
     */
    double MEHPForceBalance::getDisplacementResidualNormFor(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      const double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(this->initialConfig,
                                             this->currentSpringPartitionsVec,
                                             oneOverSpringPartitionUpperLimit);
      double elasticFreeEnergy = 0.;

      Eigen::VectorXd forces = Eigen::VectorXd::Zero(3 * net.nrOfLinks);
      // for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
      //   if (net.springPartIndexA[i] == net.springPartIndexB[i]) {
      //     // primary loops cancel out anyway
      //     continue;
      //   }
      //   Eigen::Vector3d dist =
      //     this->evaluatePartialSpringDistance(net, u, i, this->is2D);
      //   const double contourLengthFraction = springPartitions[i];
      //   const double N =
      //     net.springsContourLength[net.partialToFullSpringIndex[i]];
      //   double oneOverContourLengthFraction =
      //   CLAMP_ONE_OVER_SPRINGPARTITION(
      //     net.partialSpringIsPartial[i],
      //     ((contourLengthFraction > 0.) ? 1.0 / (N * contourLengthFraction)
      //                                   : 0.),
      //     N,
      //     oneOverSpringPartitionUpperLimit);
      //   assert(APPROX_EQUAL(
      //     oneOverContourLengthFraction, oneOverSpringPartitions[3 * i],
      //     1e-9));
      //   elasticFreeEnergy +=
      //     dist.squaredNorm() * 0.5 * oneOverContourLengthFraction;
      //   forces.segment(3 * net.springPartIndexA[i], 3) +=
      //     dist * oneOverContourLengthFraction;
      //   forces.segment(3 * net.springPartIndexB[i], 3) -=
      //     dist * oneOverContourLengthFraction;
      // }
      Eigen::VectorXi debugNrSpringsVisited =
        Eigen::VectorXi::Zero(net.nrOfPartialSprings);
      for (size_t i = 0; i < net.nrOfLinks; ++i) {
        forces.segment(3 * i, 3) =
          this->evaluateForceOnLink(i,
                                    net,
                                    u,
                                    springPartitions,
                                    debugNrSpringsVisited,
                                    oneOverSpringPartitionUpperLimit);
      }
      assert((debugNrSpringsVisited.array() == 2).all());

      // return elasticFreeEnergy;
      return forces.squaredNorm();
    }

    /**
     * @brief Compute the displacement residual norm for a specific
     * configuration
     *
     * @param net
     * @param u
     * @param oneOverSpringPartitions
     * @return double
     */
    double MEHPForceBalance::getDisplacementResidualNormFor(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& oneOverSpringPartitions) const
    {
      Eigen::VectorXd displacedCoords = net.coordinates + u;

      Eigen::VectorXd relevantPartialDistances =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;

      if (this->assumeBoxLargeEnough) {
        this->box.handlePBC(relevantPartialDistances);
      }

      if (this->is2D) {
        for (size_t i = 2; i < relevantPartialDistances.size(); i += 3) {
          relevantPartialDistances[i] = 0.;
        }
      }

#ifndef NDEBUG
      for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
        Eigen::Vector3d dist = this->evaluatePartialSpringDistance(net, u, i);
        assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
          dist, relevantPartialDistances.segment(3 * i, 3), 1e-9));
      }
#endif

      assert(relevantPartialDistances.size() == oneOverSpringPartitions.size());
      Eigen::VectorXd partialDistancesOverSpringPartitions =
        (relevantPartialDistances.array() * oneOverSpringPartitions.array())
          .matrix();

      Eigen::VectorXd overallForces = Eigen::VectorXd::Zero(3 * net.nrOfLinks);
      overallForces(net.springPartCoordinateIndexB) -=
        partialDistancesOverSpringPartitions;
      overallForces(net.springPartCoordinateIndexA) +=
        partialDistancesOverSpringPartitions;

      return overallForces.squaredNorm();
    }

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
        double contourLengthFraction = springPartitions0[i];
        double valueToSet =
          CLAMP_ONE_OVER_SPRINGPARTITION(net.partialSpringIsPartial[i],
                                         1.0 / (N * contourLengthFraction),
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
     * @brief Assemble all indices of partial springs for a particular
     * slip-link
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
     * @brief Remove double listed springs from cross-linkers
     *
     * @param net
     */
    void MEHPForceBalance::removeDuplicateListedSpringsFromLinks(
      ForceBalanceNetwork& net) const
    {
      for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
        this->removeDuplicateListedSpringsFromLink(net, linkIdx);
      }

#ifndef NDEBUG
      this->validateNetwork();
#endif
    }

    void MEHPForceBalance::removeDuplicateListedSpringsFromLink(
      ForceBalanceNetwork& net,
      size_t linkIdx,
      bool allowOnEntanglement) const
    {
      // remove duplicate mentions of the same spring index
      std::sort(net.springIndicesOfLinks[linkIdx].begin(),
                net.springIndicesOfLinks[linkIdx].end());
      auto last = std::unique(net.springIndicesOfLinks[linkIdx].begin(),
                              net.springIndicesOfLinks[linkIdx].end());
      if (last != net.springIndicesOfLinks[linkIdx].end()) {
#ifdef DEBUG_REMOVAL
        std::cout << "Removed duplicate spring indices from link " << linkIdx
                  << std::endl;
#endif
        RUNTIME_EXP_IFN(
          linkIdx > net.nrOfNodes ||
            (net.oldAtomTypes[linkIdx] != this->entanglementType) ||
            allowOnEntanglement,
          "Require entanglement beads to not form primary loops.");
        net.springIndicesOfLinks[linkIdx].erase(
          last, net.springIndicesOfLinks[linkIdx].end());
      }
    }

    size_t MEHPForceBalance::removePrimaryLoops(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions) const
    {
      size_t numRemoved = 0;
      for (long int springIdx = net.nrOfSprings - 1; springIdx >= 0;
           --springIdx) {
        if (net.springIndexA[springIdx] == net.springIndexB[springIdx] &&
            net.localToGlobalSpringIndex[springIdx].size() == 1 &&
            net.springsType[springIdx] != this->entanglementType) {
          this->removeSpringFollowingEntanglementLinks(
            net, displacements, springPartitions, springIdx);
          springIdx = std::min<long int>(springIdx, net.nrOfSprings - 1);
        }
      }

#ifndef NDEBUG
      this->validateNetwork();
#endif
    }

    /**
     * @brief Remove crosslinkers which do not have any springs with a certain
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
      //        this->removePrimaryLoops(net, displacements, springPartitions);
      // this->validateNetwork(net, displacements, springPartitions);
      // first, we remove all inactive springs
      for (long int springIdx = net.nrOfSprings - 1; springIdx >= 0;
           --springIdx) {
        if (springIdx >= net.nrOfSprings) {
          springIdx = net.nrOfSprings - 1;
          continue;
        }
        if (net.springsType[springIdx] == this->entanglementType) {
          // let's not remove entanglement springs
          continue;
        }
        if (net.oldAtomTypes[net.springIndexA[springIdx]] ==
              this->entanglementType &&
            net.oldAtomTypes[net.springIndexB[springIdx]] ==
              this->entanglementType) {
          // this is a quasi-partial spring, will be handled differently
          continue;
        }
        std::vector<size_t> involvedPartialSprings =
          this->getAllPartialSpringIndicesAlong(net, springIdx);
        assert(involvedPartialSprings.size() > 0);
        bool isActive = false;
        for (size_t partialSpringIdx : involvedPartialSprings) {
          RUNTIME_EXP_IFN(
            net.coordinates.size() == displacements.size(),
            "Expected coordinates and displacements to have same size, got " +
              std::to_string(net.coordinates.size()) + " and " +
              std::to_string(displacements.size()) + ".");
          // assert(net.coordinates.size() == displacements.size());
          Eigen::Vector3d distance = this->evaluatePartialSpringDistance(
            net, displacements, partialSpringIdx);
          double contourLength =
            net.springsContourLength
              [net.partialToFullSpringIndex[partialSpringIdx]];
          double partition = springPartitions[partialSpringIdx];
          if (!this->distanceIsWithinTolerance(
                distance, tolerance, contourLength, partition)) {
            isActive = true;
            break;
          }
        }
        if (!isActive) {
// remove this spring
#ifdef DEBUG_REMOVAL
          std::cout << "Removing inactive spring " << springIdx
                    << " with all dependencies" << std::endl;
#endif
          this->removeSpringFollowingEntanglementLinks(
            net, displacements, springPartitions, springIdx);

#ifndef NDEBUG
          this->validateNetwork(net, displacements, springPartitions);
#endif
          numRemoved += 1;
        }
      }

      // then, we remove all crosslinkers that are 0- or 1-functional
      for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
           --crosslinkIdx) {
        assert(net.springIndicesOfLinks.size() > crosslinkIdx);
        if (net.springIndicesOfLinks[crosslinkIdx].size() == 0 // f = 0
        ) {
#ifdef DEBUG_REMOVAL
          std::cout << "Removing f = 0 x-link " << crosslinkIdx << std::endl;
#endif

          this->removeLink(net, displacements, crosslinkIdx);
          numRemoved += 1;
#ifndef NDEBUG
          this->validateNetwork(net, displacements, springPartitions);
#endif
        }

        else if ( // or f = 1, NOT primary loop
          (net.springIndicesOfLinks[crosslinkIdx].size() == 1) &&
          (XOR(
            net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx][0]]
                                    [0] == crosslinkIdx,
            pylimer_tools::utils::last(
              net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx]
                                                               [0]]) ==
              crosslinkIdx))) {
          assert(net.oldAtomTypes[crosslinkIdx] != this->entanglementType);
#ifdef DEBUG_REMOVAL
          std::cout << "Removing f = 1 x-link " << crosslinkIdx << std::endl;
#endif
          std::vector<size_t> entanglementLinksRemoved =
            this->getEntanglementLinkIndicesAlong(
              net, net.springIndicesOfLinks[crosslinkIdx][0]);
          // need to first remove the spring
          const long int previousId = net.oldAtomIds[crosslinkIdx];
          this->removeSpringFollowingEntanglementLinks(
            net,
            displacements,
            springPartitions,
            net.springIndicesOfLinks[crosslinkIdx][0]);
          // this also removed some intermediate entanglement links
          // -> crossLinkIdx is now at a different index
          size_t newCrosslinkIdx = crosslinkIdx;
          for (size_t linkIdx : entanglementLinksRemoved) {
            assert(linkIdx != crosslinkIdx);
            if (linkIdx < crosslinkIdx) {
              newCrosslinkIdx -= 1;
            }
          }
          assert(net.oldAtomIds[newCrosslinkIdx] == previousId);
          assert(net.oldAtomTypes[newCrosslinkIdx] != this->entanglementType);
          assert(net.springIndicesOfLinks[newCrosslinkIdx].size() == 0);
          numRemoved += 1;
          // to then remove the cross-link
          this->removeLink(net, displacements, newCrosslinkIdx);
#ifdef DEBUG_REMOVAL
          std::cout << "Effectively removed f = 1 x-link " << newCrosslinkIdx
                    << std::endl;
#endif
          crosslinkIdx = std::min<long int>(crosslinkIdx, net.nrOfNodes - 1);
          // => we should ever only have 2-functional entanglement links that
          // could be merged after this.

#ifndef NDEBUG
          this->validateNetwork(net, displacements, springPartitions);
#endif
        }
      }

#ifndef NDEBUG
      this->validateNetwork(net, displacements, springPartitions);
#endif

      return numRemoved;
    }

    /**
     * @brief Remove springs that exert a stress higher than
     * `this->springBreakingLength`
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @return size_t the number of springs broken
     */
    size_t MEHPForceBalance::breakTooLongSprings(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions) const
    {
      if (this->springBreakingLength <= 0.) {
        return 0;
      }

      size_t numBroken = 0;

      // iterate the springs, determine their distance, and determine if it
      // exceeds the breaking force
      for (long int partialSpringIdx = net.nrOfPartialSprings;
           partialSpringIdx >= 0;
           --partialSpringIdx) {
        if (partialSpringIdx >= net.nrOfPartialSprings) {
          partialSpringIdx = net.nrOfPartialSprings - 1;
        }
        double len = this->getWeightedPartialSpringLength(
          net, displacements, springPartitions, partialSpringIdx);
        if (len > this->springBreakingLength) {
          // break this spring
          numBroken += 1;
          this->breakPartialSpring(
            net, displacements, springPartitions, partialSpringIdx);
        }
      }

      return numBroken;
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
      std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
        this->universe.getChainsWithCrosslinker(crossLinkerType);
      bool danglingChainsAreKept =
        this->initialConfig.nrOfNodes >
        this->universe.getAtomsOfType(crossLinkerType).size();
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
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        pylimer_tools::entities::Molecule chain = crossLinkerChains[i];
        RUNTIME_EXP_IFN(chain.getType() !=
                          pylimer_tools::entities::MoleculeType::UNDEFINED,
                        "Couldn't determine molecule type.");
        if (chain.getType() ==
              pylimer_tools::entities::MoleculeType::PRIMARY_LOOP ||
            chain.getType() ==
              pylimer_tools::entities::MoleculeType::NETWORK_STRAND ||
            (chain.getType() ==
               pylimer_tools::entities::MoleculeType::DANGLING_CHAIN &&
             danglingChainsAreKept)) {
          assert(i == this->initialConfig.springToMoleculeIds[springId]);
          // TODO: also check that this is not a higher order dangling strand
          nrOfEligibleAtoms +=
            crossLinkerChains[i].getNrOfAtoms() -
            crossLinkerChains[i].getAtomsOfType(this->crossLinkerType).size();
          std::vector<pylimer_tools::entities::Atom> atoms =
            crossLinkerChains[i].getAtomsLinedUp(this->crossLinkerType);
          for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
            pylimer_tools::entities::Atom atom = atoms[atomIdx];
            if (atom.getType() != this->crossLinkerType) {
              eligibleAtoms.push_back(atom);
              vertexIdxIsEligible[this->universe.getIdxByAtomId(atom.getId())] =
                true;
              atomToStrand.emplace(atom.getId(), springId);
              atomIdxInStrand.emplace(atom.getId(), atomIdx);
              RUNTIME_EXP_IFN(this->initialConfig.oldAtomIdToSpringIndex.at(
                                atom.getId()) == springId,
                              "The spring numbering seems incorrect. Placing "
                              "slip-links will "
                              "lead to inappropriate placement.");
            }
          }
          RUNTIME_EXP_IFN(this->initialConfig.springToMoleculeIds[springId] ==
                            i,
                          "The spring numbering seems incorrect. Placing "
                          "slip-links will lead to inappropriate placement.");
          // #ifndef NDEBUG
          std::vector<pylimer_tools::entities::Atom> chainEnds =
            crossLinkerChains[i].getChainEnds(this->crossLinkerType, true);
          RUNTIME_EXP_IFN(
            this->initialConfig.oldAtomIds
                [this->initialConfig.linkIndicesOfSprings[springId][0]] ==
              chainEnds[0].getId(),
            "Atom ends are inconsistent when sampling.");
          RUNTIME_EXP_IFN(chainEnds[0].getId() == atoms[0].getId(),
                          "Atom ends are inconsistent when sampling.");
          // #endif
          springId += 1;
        }
      }
      // build neighbourlist
      std::vector<pylimer_tools::entities::Atom> atomsForNeighbourList =
        this->universe.getAtoms();
      std::vector<bool> isMasked = pylimer_tools::utils::initializeWithValue(
        this->universe.getNrOfAtoms(), false);
      if (excludeCrosslinks) {
        // TODO: check whether it is faster to just only query the other ones
        std::erase_if(atomsForNeighbourList,
                      [&](const pylimer_tools::entities::Atom& a) -> bool {
                        return a.getType() == this->crossLinkerType;
                      });
        for (size_t i = 0; i < this->universe.getNrOfAtoms(); ++i) {
          isMasked[i] = (this->universe.getAtomByVertexIdx(i).getType() ==
                         this->crossLinkerType);
        }
        for (pylimer_tools::entities::Atom a : atomsForNeighbourList) {
          RUNTIME_EXP_IFN(
            a.getType() != this->crossLinkerType,
            "Removing cross-linkers from neighborlist did not seem to work.");
        }
      }
      pylimer_tools::entities::NeighbourList neighbourList =
        pylimer_tools::entities::NeighbourList(
          atomsForNeighbourList, this->universe.getBox(), cutoff);

      std::random_device rd{};
      std::mt19937 rng = std::mt19937(seed > 0 ? seed : rd());
      // std::cout << "Initial sampling rng seed: " << rng << std::endl;
      // build list of random samples
      // this way is more performant than
      // sampling integers and checking whether they have been sampled already
      std::vector<size_t> toSampleFrom;
      toSampleFrom.reserve(this->universe.getNrOfAtoms());
      for (size_t i = 0; i < this->universe.getNrOfAtoms(); ++i) {
        toSampleFrom.push_back(i);
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
        while (isMasked[toSampleFrom[sampleIdx]]) {
          sampleIdx += 1;

          if (sampleIdx >= toSampleFrom.size()) {
            break;
          }
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
        RUNTIME_EXP_IFN(!excludeCrosslinks ||
                          a1.getType() != this->crossLinkerType,
                        "Sampled atom is cross-link, but may not be");
        isMasked[sampledVertexId] = true;
        // then, find neighbouring atoms (but not from the same strand?!)
        std::vector<pylimer_tools::entities::Atom> neighbours =
          neighbourList.getAtomsCloseTo(a1);
        neighbourList.removeAtom(a1,
                                 "After querying neighbours. Impossible case.");
        // filter the neighbours to include only those from other strands
        // NOTE: this skews the whole thing a bit
        neighbours.erase(
          std::remove_if(
            neighbours.begin(),
            neighbours.end(),
            [&](pylimer_tools::entities::Atom a) -> bool {
              return (
                (atomToStrand[a.getId()] ==
                   atomToStrand[a1.getId()] // do not use "at", because not
                 // all atoms in the neighbours
                 // have been assigned a strand
                 && (std::abs(static_cast<double>(
                       atomIdxInStrand[a.getId()] -
                       atomIdxInStrand[a1.getId()])) < sameStrandCutoff))
                // the following check should not be necessary?!?
                || isMasked[this->universe.getIdxByAtomId(a.getId())]);
            }),
          neighbours.end());
        if (neighbours.size() == 0) {
          std::cerr << "Not enough close neighbours found." << std::endl;
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
        RUNTIME_EXP_IFN(sampledVertexId2 != sampledVertexId,
                        "Second sample may not be equal to the first");
        RUNTIME_EXP_IFN(!excludeCrosslinks ||
                          a2.getType() != this->crossLinkerType,
                        "Sampled atom is cross-link, but may not be");
        RUNTIME_EXP_IFN(!isMasked[sampledVertexId2],
                        "Sampled vertex 2 is masked, but could not be");
        neighbourList.removeAtom(a2, "Removing sampled vertex 2");
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
        // std::cout << "OMerging atoms " << a1.getId() << " and " <<
        // a2.getId()
        // << " with distance " << a1.distanceTo(a2, universe.getBox()) <<
        // std::endl;
        Eigen::Vector3d meanPositions = a1.meanPositionWith(a2, box);
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
      // std::cout << "Detecting slip-links based on cycles. Base memory
      // useage:
      // "
      //           << getCurrentRSS() << ", peak " << getPeakRSS() <<
      //           std::endl;
      std::vector<std::vector<long int>> loopEdges;
      std::vector<std::vector<long int>> loops = this->universe.findLoops(
        this->crossLinkerType, maxLoopLength, false, &loopEdges);
      std::cout << "Detected " << loops.size() << " loops." << std::endl;
      // reduced loops = loops, but only the (new) spring indices
      std::vector<std::vector<size_t>> reducedLoops;
      reducedLoops.reserve(loops.size());
      // min & max coordinates of the loops in order to easily know
      // a priori whether it makes sense to compare two loops or not
      std::vector<std::array<double, 3>> loopMinCoords;
      loopMinCoords.reserve(loops.size());
      std::vector<std::array<double, 3>> loopMaxCoords;
      loopMaxCoords.reserve(loops.size());
      pylimer_tools::entities::Box box = this->universe.getBox();
      for (std::vector<long int> loop : loops) {
        // max & min coordinates
        // TODO: implement, such that the max & min work also with infinite
        // loops
        // assume each subsequent atom is bonded with a bond shorter
        // than half the bond
        Eigen::VectorXd alignedLoopCoordinates =
          this->universe.getUnwrappedVertexCoordinates(loop, box);

        std::array<double, 3> minCoords;
        std::array<double, 3> maxCoords;

        for (int i = 0; i < 3; ++i) {
          minCoords[i] =
            alignedLoopCoordinates(Eigen::seq(i, Eigen::last - 2 + i, 3))
              .minCoeff();
          maxCoords[i] =
            alignedLoopCoordinates(Eigen::seq(i, Eigen::last - 2 + i, 3))
              .maxCoeff();
          // want to see later that this direction is spanned fully
          // that's why we just push it so far that it is within the box
          if (maxCoords[i] - minCoords[i] > box.getL(i) ||
              (std::fabs(alignedLoopCoordinates[i] -
                         alignedLoopCoordinates[3 * (loop.size() - 1) + i]) >
               0.5 * box.getL(i))) {
            // special case in case the loop is periodic
            minCoords[i] = box.getLowL(i);
            maxCoords[i] = box.getHighL(i);
          } else {
            // adjust for box
            while (minCoords[i] > box.getHighL(i)) {
              minCoords[i] -= box.getL(i);
            }
            while (minCoords[i] < box.getLowL(i)) {
              minCoords[i] += box.getL(i);
            }
            while (maxCoords[i] > box.getHighL(i)) {
              maxCoords[i] -= box.getL(i);
            }
            while (maxCoords[i] < box.getLowL(i)) {
              maxCoords[i] += box.getL(i);
            }
            // at this point, it is not given anymore that minCoords[i] <
            // maxCoords[i], instead, we know that both are within the box
          }
        }

        loopMinCoords.push_back(minCoords);
        loopMaxCoords.push_back(maxCoords);

        // loop reduction to springs
        std::set<size_t> reducedLoop;
        for (size_t i = 0; i < loop.size(); ++i) {
          // reduced loop
          if (this->universe.getPropertyValue<long int>("type", loop[i]) !=
              this->crossLinkerType) {
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
      // std::cout << "After finding loops, memory useage: " <<
      // getCurrentRSS()
      //           << ", peak " << getPeakRSS() << std::endl;

      // fetch some data to later estimate alpha & beta
      std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
        this->universe.getChainsWithCrosslinker(crossLinkerType);
      std::unordered_map<size_t, size_t> atomToStrand;
      atomToStrand.reserve(this->universe.getNrOfAtoms());
      std::unordered_map<size_t, size_t> atomIdxInStrand;
      atomIdxInStrand.reserve(this->universe.getNrOfAtoms());
      size_t springId = 0;
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        pylimer_tools::entities::Molecule chain = crossLinkerChains[i];
        std::vector<pylimer_tools::entities::Atom> atoms =
          crossLinkerChains[i].getAtomsLinedUp(this->crossLinkerType);
        for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
          pylimer_tools::entities::Atom atom = atoms[atomIdx];
          if (atom.getType() != this->crossLinkerType) {
            atomToStrand.emplace(atom.getId(), springId);
            atomIdxInStrand.emplace(atom.getId(), atomIdx);
          }
        }
      }

      // std::cout << "After mapping chains again, memory useage: "
      //           << getCurrentRSS() << ", peak " << getPeakRSS() <<
      //           std::endl;

      // the resulting vectors to fill
      std::unordered_map<long int, long int> intersectionsOfEdges;
      typedef std::tuple<pylimer_tools::entities::LoopIntersectionInfo,
                         std::set<size_t>>
        intersection_loops_tuple;
      std::vector<intersection_loops_tuple> relevantIntersections;

      // as of
      // https://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way
      auto hash_integer_pair = [](long int a, long int b) -> long int {
        return a >= b ? a * a + a + b : a + b * b; // where a, b >= 0
      };
      std::cout << "Searching for intersections..." << std::endl;

      for (size_t i = 0; i < loops.size(); ++i) {
        // reserve enough space
        size_t estimateNrOfSliplinks =
          loops.size() * loops.size() *
          this->initialConfig.meanSpringContourLength;
        if (i > 0) {
          estimateNrOfSliplinks =
            std::max(estimateNrOfSliplinks,
                     (loops.size() / i) * relevantIntersections.size());
        }
        relevantIntersections.reserve(estimateNrOfSliplinks);

        const std::array<double, 3> loop_i_min = loopMinCoords[i];
        const std::array<double, 3> loop_i_max = loopMaxCoords[i];

        // TODO: instead of the N^2 loop, might want to try some filter at
        // least also, we ignore all self-entanglements of one loop with
        // itself.
        for (size_t j = i + 1; j < loops.size(); ++j) {
          // first check if the loops have any overlap in 3D, otherwise, we
          // can skip them anyway.

          std::array<double, 3> loop_j_min = loopMinCoords[j];
          std::array<double, 3> loop_j_max = loopMaxCoords[j];

          for (int dir = 0; dir < 3; ++dir) {
            if (
              // "normal" case
              ((loop_i_min[dir] <= loop_i_max[dir] &&
                loop_j_min[dir] <= loop_j_max[dir]) &&
               !(loop_j_min[dir] <= loop_i_max[dir] &&
                 loop_j_max[dir] >= loop_i_min[dir])) ||
              // periodic case in j
              ((loop_i_min[dir] <= loop_i_max[dir] &&
                loop_j_min[dir] > loop_j_max[dir]) &&
               (loop_i_min[dir] >= loop_j_max[dir] &&
                loop_i_max[dir] <= loop_j_min[dir])) ||
              // periodic case in i
              ((loop_i_min[dir] > loop_i_max[dir] &&
                loop_j_min[dir] <= loop_j_max[dir]) &&
               (loop_j_min[dir] >= loop_i_max[dir] &&
                loop_j_max[dir] <= loop_i_min[dir])) ||
              // both periodic -> both pass the boundary -> overlap anyway
              ((loop_i_min[dir] > loop_i_max[dir] &&
                loop_j_min[dir] > loop_j_max[dir]) &&
               false)) {
              //  -> no overlap in dir, skip this comparison
              goto loop_j_checked;
            }
          }

          // then, if we did not goto, we can actually find the entanglements
          {
            std::vector<pylimer_tools::entities::LoopIntersectionInfo>
              intersections = this->universe.findLoopEntanglements(
                loops[i], loops[j], loopEdges[i], loopEdges[j]);
            for (pylimer_tools::entities::LoopIntersectionInfo intersection :
                 intersections) {
              bool keepIntersection = true;

              // check if we want to keep it
              // TODO: check if we want to keep it based on the direction /
              // distance of the intersection

              // check if we want to keep it, based on the edges involved
              long int edgePairHash =
                hash_integer_pair(intersection.edge1, intersection.edge2);
              if (pylimer_tools::utils::map_has_key(intersectionsOfEdges,
                                                    edgePairHash)) {
                keepIntersection = false;
                long int involvedIntersection =
                  intersectionsOfEdges.at(edgePairHash);
                std::get<1>(relevantIntersections[involvedIntersection])
                  .insert(i);
                std::get<1>(relevantIntersections[involvedIntersection])
                  .insert(j);
              }

              if (keepIntersection) {
                std::set<size_t> localSet;
                localSet.insert(i);
                localSet.insert(j);
                intersection_loops_tuple relevantIntersection =
                  std::make_tuple(intersection, localSet);
                relevantIntersections.push_back(relevantIntersection);
                intersectionsOfEdges.insert_or_assign(
                  edgePairHash,
                  static_cast<long int>(relevantIntersections.size()) - 1);
              }
            }
          }

        loop_j_checked:;
        }

        // make space: cleanup the loop i
        loops[i] = std::vector<long int>();
        // std::cout
        //   << "After checking intersections of loop " << i << " ("
        //   << relevantIntersections.size()
        //   << " relevant intersections found in total yet), memory useage: "
        //   << getCurrentRSS() << ", peak " << getPeakRSS() << std::endl;
      }

      // Actually make them into slip-links...
      std::vector<double> slipLinkXs;
      slipLinkXs.reserve(relevantIntersections.size());
      std::vector<double> slipLinkYs;
      slipLinkYs.reserve(relevantIntersections.size());
      std::vector<double> slipLinkZs;
      slipLinkZs.reserve(relevantIntersections.size());
      std::vector<size_t> slipLinkStrandA;
      slipLinkStrandA.reserve(relevantIntersections.size());
      std::vector<size_t> slipLinkStrandB;
      slipLinkStrandB.reserve(relevantIntersections.size());
      std::vector<double> slipLinkStrandAlpha;
      slipLinkStrandAlpha.reserve(relevantIntersections.size());
      std::vector<double> slipLinkStrandBeta;
      slipLinkStrandBeta.reserve(relevantIntersections.size());
      std::vector<std::vector<size_t>> slipLinksLoops;
      slipLinksLoops.reserve(relevantIntersections.size());

      for (intersection_loops_tuple intersectionAndLoops :
           relevantIntersections) {
        pylimer_tools::entities::LoopIntersectionInfo intersection =
          std::get<0>(intersectionAndLoops);
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
        std::vector<size_t> localLoops(
          std::get<1>(intersectionAndLoops).begin(),
          std::get<1>(intersectionAndLoops).end());
        slipLinksLoops.push_back(localLoops);
      }

      std::cout << "Found " << slipLinksLoops.size() << " intersections."
                << std::endl;
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
#ifdef DEBUG_REMOVAL
      std::cout << "Starting to remove spring " << springIdx << std::endl;
#endif
      INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                         "Can only remove springs, not partial springs.");
      Eigen::VectorXd allTotalSpringDistancesBefore =
        this->evaluateSpringLengths(net, displacements, this->is2D);
      std::vector<size_t> affectedLinks = net.linkIndicesOfSprings[springIdx];
      std::vector<size_t> uniqueAffectedLinks =
        net.linkIndicesOfSprings[springIdx];
      std::sort(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end());
      uniqueAffectedLinks.erase(
        std::unique(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end()),
        uniqueAffectedLinks.end());

      // remove the link to the link, höhö
      for (size_t affectedLinkIdx : uniqueAffectedLinks) {
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

        size_t found =
          std::erase(net.springIndicesOfLinks[affectedLinkIdx], springIdx);

        RUNTIME_EXP_IFN(found > 0,
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
      pylimer_tools::utils::removeRow(net.springsType, springIdx);
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
        pylimer_tools::utils::removeRows(
          net.springPartBoxOffset, 3 * partialSpringIdx, 3);
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
      if (linksToRemove.size() > 2) {
        assert(linksToRemove[0] > linksToRemove[1]);
      }

      Eigen::ArrayXb springIsAffected =
        Eigen::ArrayXb::Constant(net.nrOfSprings, false);
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
        // RUNTIME_EXP_IFN(springsOfLink.size() % 2 == 0, "Expected link to
        // have an even number of components, got " +
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
          size_t partialSpringToKeep =
            std::min(involvedPartialSprings[0], involvedPartialSprings[1]);
          size_t partialSpringToRemove =
            std::max(involvedPartialSprings[0], involvedPartialSprings[1]);
          assert(partialSpringToKeep != partialSpringToRemove);
          assert(net.partialToFullSpringIndex[partialSpringToKeep] ==
                 net.partialToFullSpringIndex[partialSpringToRemove]);
          // actually do the merge
          this->mergePartialSprings(net,
                                    displacements,
                                    springPartitions,
                                    partialSpringToRemove,
                                    partialSpringToKeep,
                                    slipLinkIdx);
          // total distance changes -> cannot use for checking the total
          // distance
          springIsAffected[net.partialToFullSpringIndex[partialSpringToKeep]] =
            true;
        }

        assert(net.springIndicesOfLinks[slipLinkIdx].empty());

        // then, actually remove the slip-link
#ifdef DEBUG_REMOVAL
        std::cout << "Removing slip-link " << slipLinkIdx << std::endl;
#endif
        this->removeLink(net, displacements, slipLinkIdx);
      }
#ifdef DEBUG_REMOVAL
      std::cout << "Removed spring " << springIdx << std::endl;
#endif

      Eigen::VectorXd allTotalSpringDistancesAfter =
        this->evaluateSpringLengths(net, displacements, this->is2D);
      assert(allTotalSpringDistancesAfter.size() == net.nrOfSprings);
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        size_t correspondingOldIdx = i >= springIdx ? i + 1 : i;
        if (springIsAffected[i]) {
          // when slip-links are removed, the overall distance must reduce.
          RUNTIME_EXP_IFN(allTotalSpringDistancesBefore[correspondingOldIdx] +
                              1e-9 >=
                            allTotalSpringDistancesAfter[i],
                          "Expected that the total distances stay constant for "
                          "non-changed springs.");
        } else {
          RUNTIME_EXP_IFN(
            APPROX_EQUAL(allTotalSpringDistancesBefore[correspondingOldIdx],
                         allTotalSpringDistancesAfter[i],
                         1e-9),
            "Expected that the total distances stay constant for non-changed "
            "springs.");
        }
      }
    }

    /**
     * @brief Remove a spring, but also all springs that are connected to it
     * and are connected via entanglement links.
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param springIdx
     */
    void MEHPForceBalance::removeSpringFollowingEntanglementLinks(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions,
      const size_t springIdx) const
    {
      std::vector<size_t> springsToRemove =
        this->getAllFullSpringIndicesAlong(net, springIdx);
      std::vector<size_t> linksToRemove =
        this->getEntanglementLinkIndicesAlong(net, springIdx);
      assert(springsToRemove.size() > 0);
      std::sort(
        springsToRemove.begin(), springsToRemove.end(), std::greater<>());
      for (size_t springIdxToDelete : springsToRemove) {
        assert(net.springsType[springIdxToDelete] != this->entanglementType);
        this->removeSpring(
          net, displacements, springPartitions, springIdxToDelete);
      }
      std::sort(linksToRemove.begin(), linksToRemove.end(), std::greater<>());
      for (size_t linkIdxToDelete : linksToRemove) {
        assert(net.springIndicesOfLinks[linkIdxToDelete].size() <= 1);
        assert(net.oldAtomTypes[linkIdxToDelete] == this->entanglementType);
        if (net.springIndicesOfLinks[linkIdxToDelete].size() == 1) {
          assert(
            net.springsType[net.springIndicesOfLinks[linkIdxToDelete][0]] ==
            this->entanglementType);
#ifdef DEBUG_REMOVAL
          std::cout << "Removing additional spring between entanglement links "
                    << net.springIndicesOfLinks[linkIdxToDelete][0]
                    << std::endl;
#endif

          this->removeSpring(net,
                             displacements,
                             springPartitions,
                             net.springIndicesOfLinks[linkIdxToDelete][0]);
        }
        this->removeLink(net, displacements, linkIdxToDelete);
      }
    };

    /**
     * @brief break a spring, given its partial spring index
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param partialSpringIdx
     */
    void MEHPForceBalance::breakPartialSpring(
      ForceBalanceNetwork& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions,
      const size_t partialSpringIdx) const
    {
      this->removeSpringFollowingEntanglementLinks(
        net,
        displacements,
        springPartitions,
        net.partialToFullSpringIndex[partialSpringIdx]);
    };

    /**
     * @brief remove a link from the network
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
                         "The springs have to be removed or re-linked before "
                         "removing the link.");
#ifdef DEBUG_REMOVAL
      std::cout << "Removing link " << linkIdx << std::endl;
#endif

      pylimer_tools::utils::removeRows(net.coordinates, linkIdx * 3, 3);
      pylimer_tools::utils::removeRows(displacements, linkIdx * 3, 3);

      if (!net.linkIsSliplink[linkIdx]) {
        net.nrOfNodes -= 1;
        net.oldAtomIdToSpringIndex.erase(net.oldAtomIds[linkIdx]);
        pylimer_tools::utils::removeRow(net.oldAtomIds, linkIdx);
        pylimer_tools::utils::removeRow(net.oldAtomTypes, linkIdx);
      } else {
        pylimer_tools::utils::removeRow(net.nrOfCrosslinkSwapsEndured,
                                        linkIdx - net.nrOfNodes);
        if (net.loopsOfSliplink.size() > 0) {
          net.loopsOfSliplink.erase(net.loopsOfSliplink.begin() + linkIdx);
        }
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
     * @brief Combine two springs
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance::mergePartialSprings(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      const size_t removedPartialSpringIdx,
      const size_t keptPartialSpringIdx,
      const size_t linkToReduce,
      bool skipEigenResize) const
    {
      INVALIDARG_EXP_IFN(net.linkIsSliplink[linkToReduce],
                         "The link to reduce must be a slip-link");
      INVALIDARG_EXP_IFN(keptPartialSpringIdx != removedPartialSpringIdx,
                         "Cannot merge one spring with the same one.");
      INVALIDARG_EXP_IFN(
        net.partialToFullSpringIndex[keptPartialSpringIdx] ==
          net.partialToFullSpringIndex[removedPartialSpringIdx],
        "The partial springs must be part of the same spring to merge them.");
      INVALIDARG_EXP_IFN(
        (net.springPartIndexA[keptPartialSpringIdx] == linkToReduce &&
         net.springPartIndexB[removedPartialSpringIdx] == linkToReduce) ||
          (net.springPartIndexB[keptPartialSpringIdx] == linkToReduce &&
           net.springPartIndexA[removedPartialSpringIdx] == linkToReduce),
        "Link to reduce must be part of the partial springs "
        "that are to be removed and kept, got " +
          std::to_string(net.springPartIndexA[keptPartialSpringIdx]) + " and " +
          std::to_string(net.springPartIndexB[keptPartialSpringIdx]) +
          " in kept spring, and got " +
          std::to_string(net.springPartIndexA[removedPartialSpringIdx]) +
          " and " +
          std::to_string(net.springPartIndexB[removedPartialSpringIdx]) +
          "in removed spring, instead of " + std::to_string(linkToReduce) +
          ".");
#ifdef DEBUG_REMOVAL
      std::cout << "Merging partial springs " << removedPartialSpringIdx
                << " and " << keptPartialSpringIdx << " around " << linkToReduce
                << std::endl;
#endif

      Eigen::Vector3d distanceBefore =
        this->evaluatePartialSpringDistance(
          net, u, removedPartialSpringIdx, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, keptPartialSpringIdx, this->is2D, false);
      size_t fullSpringIdx = net.partialToFullSpringIndex[keptPartialSpringIdx];
      // start with removal
      net.nrOfPartialSprings -= 1;
      // tell the kept one their new end
      // NOTE: if is possible, if the removedPartialSpring is a primary loop,
      // that this procedure is ambiguous.
      bool removedIsA =
        net.springPartIndexA[removedPartialSpringIdx] == linkToReduce;
      size_t newEnd = removedIsA
                        ? net.springPartIndexB[removedPartialSpringIdx]
                        : net.springPartIndexA[removedPartialSpringIdx];
      if (removedIsA) {
        assert(net.springPartIndexB[keptPartialSpringIdx] == linkToReduce);
        net.springPartIndexB[keptPartialSpringIdx] = newEnd;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexB[3 * keptPartialSpringIdx + dir] =
            3 * newEnd + dir;
        }
      } else {
        assert(net.springPartIndexA[keptPartialSpringIdx] == linkToReduce);
        net.springPartIndexA[keptPartialSpringIdx] = newEnd;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexA[3 * keptPartialSpringIdx + dir] =
            3 * newEnd + dir;
        }
      }
      net.springPartBoxOffset.segment(3 * keptPartialSpringIdx, 3) +=
        net.springPartBoxOffset.segment(3 * removedPartialSpringIdx, 3);
      // remove the spring from the link
      // NOTE: currently, we allow it not to be present,
      // as it might be removed earlier already
      // It is anyway the case, that this function does not necessarily
      // keep the network valid
      int found =
        std::erase(net.springIndicesOfLinks[linkToReduce], fullSpringIdx);
      // TODO: check the origin of this assertion
      // assert(found == 1 || found == 0);
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
      springPartitions[keptPartialSpringIdx] +=
        springPartitions[removedPartialSpringIdx];
      for (int j = net.linkIndicesOfSprings[fullSpringIdx].size() - 1; j >= 0;
           --j) {
        if (net.linkIndicesOfSprings[fullSpringIdx][j] == linkToReduce) {
          if (removed == 0) {
            if ((j > 0 && net.localToGlobalSpringIndex[fullSpringIdx][j - 1] ==
                            removedPartialSpringIdx) ||
                (j < net.localToGlobalSpringIndex[fullSpringIdx].size() &&
                 net.localToGlobalSpringIndex[fullSpringIdx][j] ==
                   removedPartialSpringIdx)) {
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
      this->removeDuplicateListedSpringsFromLink(net, linkToReduce);
      assert(found >= 1 && removed == 1);
      found = 0;
      for (int j = net.localToGlobalSpringIndex[fullSpringIdx].size() - 1;
           j >= 0;
           --j) {
        if (net.localToGlobalSpringIndex[fullSpringIdx][j] ==
            removedPartialSpringIdx) {
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
      net.partialSpringIsPartial[keptPartialSpringIdx] =
        net.linkIndicesOfSprings[fullSpringIdx].size() > 2;
      // actually remove the rows
      pylimer_tools::utils::removeRow(
        net.partialSpringIsPartial, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.partialToFullSpringIndex, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        springPartitions, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.springPartIndexA, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.springPartIndexB, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRows(net.springPartCoordinateIndexA,
                                       3 * removedPartialSpringIdx,
                                       3,
                                       skipEigenResize);
      pylimer_tools::utils::removeRows(net.springPartCoordinateIndexB,
                                       3 * removedPartialSpringIdx,
                                       3,
                                       skipEigenResize);
      pylimer_tools::utils::removeRows(net.springPartBoxOffset,
                                       3 * removedPartialSpringIdx,
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
              removedPartialSpringIdx) {
            net.localToGlobalSpringIndex[loopSpringIdx][i] -= 1;
          }
        }
      }

      // validation
      size_t newSpringIdx =
        keptPartialSpringIdx +
        (keptPartialSpringIdx > removedPartialSpringIdx ? -1 : 0);
      Eigen::Vector3d newDistance = this->evaluatePartialSpringDistance(
        net, u, newSpringIdx, this->is2D, false);
      RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(
                        newDistance, distanceBefore, 1e-5),
                      "After merging two partial springs, the overall distance "
                      "is not consistent. Expected distance " +
                        std::to_string(distanceBefore) + ", but got " +
                        std::to_string(newDistance) + " for spring " +
                        std::to_string(newSpringIdx) + ".");
    }

    /**
     * @brief Combine two springs
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance::mergeSprings(ForceBalanceNetwork& net,
                                        const Eigen::VectorXd& u,
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
      INVALIDARG_EXP_IFN(net.springsType[keptSpringIdx] !=
                           this->entanglementType,
                         "Should not merge entanglement springs.");
      INVALIDARG_EXP_IFN(net.springsType[removedSpringIdx] !=
                           this->entanglementType,
                         "Should not merge entanglement springs.");
#ifndef NDEBUG
      if (net.oldAtomTypes[linkToReduce] != this->entanglementType) {
        std::vector<size_t> fullSprings1 =
          this->getAllFullSpringIndicesAlong(net, removedSpringIdx);
        std::vector<size_t> fullSprings2 =
          this->getAllFullSpringIndicesAlong(net, keptSpringIdx);
        std::sort(fullSprings1.begin(), fullSprings1.end());
        std::sort(fullSprings2.begin(), fullSprings2.end());
        INVALIDARG_EXP_IFN(
          !(pylimer_tools::utils::equal(fullSprings1, fullSprings2)),
          "Cannot merge such that a primary loop made of "
          "entanglements results.");
      }
#endif

      size_t fullSpringIdx = net.partialToFullSpringIndex[keptSpringIdx];
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

      RUNTIME_EXP_IFN(
        net.springPartIndexA[removedPartialSpringIdx] == linkToReduce ||
          net.springPartIndexB[removedPartialSpringIdx] == linkToReduce,
        "Did not detect correct partial springs");
      RUNTIME_EXP_IFN(
        net.springPartIndexA[remainingPartialSpringIdx] == linkToReduce ||
          net.springPartIndexB[remainingPartialSpringIdx] == linkToReduce,
        "Did not detect correct partial springs");

      Eigen::Vector3d distanceBefore = this->evaluatePartialSpringDistance(
        net, u, removedPartialSpringIdx, this->is2D, false);
      Eigen::Vector3d distanceBeforeRemainingSpring =
        this->evaluatePartialSpringDistance(
          net, u, remainingPartialSpringIdx, this->is2D, false);

      net.nrOfSprings -= 1;
      net.nrOfPartialSprings -= 1;
      if (net.linkIndicesOfSprings[removedSpringIdx].size() > 2 &&
          net.linkIndicesOfSprings[keptSpringIdx].size() > 2) {
        net.nrOfSpringsWithPartition -= 1;
      }

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
          // ...from end
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
          // invert the direction of these transferred partial springs
          for (size_t partialSpringIdxToInvert :
               net.localToGlobalSpringIndex[removedSpringIdx]) {
            if (partialSpringIdxToInvert == removedPartialSpringIdx) {
              continue;
            }
            std::swap(net.springPartIndexA[partialSpringIdxToInvert],
                      net.springPartIndexB[partialSpringIdxToInvert]);
            for (int dir = 0; dir < 3; ++dir) {
              std::swap(
                net.springPartCoordinateIndexA[3 * partialSpringIdxToInvert +
                                               dir],
                net.springPartCoordinateIndexB[3 * partialSpringIdxToInvert +
                                               dir]);
            }
            net.springPartBoxOffset.segment(3 * partialSpringIdxToInvert, 3) *=
              -1.;
          }
          distanceBefore -= distanceBeforeRemainingSpring;
          distanceBefore *= -1.;
        } else {
          // ...from start
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
          distanceBefore += distanceBeforeRemainingSpring;
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
          distanceBefore += distanceBeforeRemainingSpring;
        } else {
          // std::cout << "Start start" << std::endl;
          // from start
          RUNTIME_EXP_IFN(removedSpringsLinks[0] == linkToReduce,
                          "No way this exception is ever shown, right?");
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

          // invert the direction of these transferred partial springs
          for (size_t partialSpringIdxToInvert :
               net.localToGlobalSpringIndex[removedSpringIdx]) {
            if (partialSpringIdxToInvert == removedPartialSpringIdx) {
              continue;
            }
            std::swap(net.springPartIndexA[partialSpringIdxToInvert],
                      net.springPartIndexB[partialSpringIdxToInvert]);
            for (int dir = 0; dir < 3; ++dir) {
              std::swap(
                net.springPartCoordinateIndexA[3 * partialSpringIdxToInvert +
                                               dir],
                net.springPartCoordinateIndexB[3 * partialSpringIdxToInvert +
                                               dir]);
            }
            net.springPartBoxOffset.segment(3 * partialSpringIdxToInvert, 3) *=
              -1.;
          }
          distanceBefore -= distanceBeforeRemainingSpring;
          distanceBefore *= -1.;
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
        this->removeDuplicateListedSpringsFromLink(
          net, linkOfRemovedSpring, linkOfRemovedSpring == linkToReduce);
      }

      for (int i = net.springIndicesOfLinks[linkToReduce].size() - 1; i >= 0;
           --i) {
        if (net.springIndicesOfLinks[linkToReduce][i] == removedSpringIdx ||
            net.springIndicesOfLinks[linkToReduce][i] == keptSpringIdx) {
          net.springIndicesOfLinks[linkToReduce].erase(
            net.springIndicesOfLinks[linkToReduce].begin() + i);
        }
      }
      this->removeDuplicateListedSpringsFromLink(net, linkToReduce);

      net.linkIndicesOfSprings.erase(net.linkIndicesOfSprings.begin() +
                                     removedSpringIdx);
      // partial springs
      if (net.partialSpringIsPartial[removedPartialSpringIdx]) {
        net.partialSpringIsPartial[remainingPartialSpringIdx] = true;
      }
      pylimer_tools::utils::removeRow(net.partialSpringIsPartial,
                                      removedPartialSpringIdx);

      bool removedIsA =
        net.springPartIndexA[removedPartialSpringIdx] == linkToReduce;
      size_t otherEndOfRemovedSpring =
        removedIsA ? net.springPartIndexB[removedPartialSpringIdx]
                   : net.springPartIndexA[removedPartialSpringIdx];
      double offsetMultiplier = removedIsA ? -1. : 1.;
      if (net.springPartIndexA[remainingPartialSpringIdx] == linkToReduce) {
        net.springPartIndexA[remainingPartialSpringIdx] =
          otherEndOfRemovedSpring;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexA[3 * remainingPartialSpringIdx + dir] =
            3 * otherEndOfRemovedSpring + dir;
        };
      } else {
        RUNTIME_EXP_IFN(
          net.springPartIndexB[remainingPartialSpringIdx] == linkToReduce, "");
        net.springPartIndexB[remainingPartialSpringIdx] =
          otherEndOfRemovedSpring;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexB[3 * remainingPartialSpringIdx + dir] =
            3 * otherEndOfRemovedSpring + dir;
        }
        offsetMultiplier *= -1.;
      }
      net.springPartBoxOffset.segment(3 * remainingPartialSpringIdx, 3) +=
        offsetMultiplier *
        net.springPartBoxOffset.segment(3 * removedPartialSpringIdx, 3);
      pylimer_tools::utils::removeRow(net.springPartIndexA,
                                      removedPartialSpringIdx);
      pylimer_tools::utils::removeRow(net.springPartIndexB,
                                      removedPartialSpringIdx);
      pylimer_tools::utils::removeRows(
        net.springPartCoordinateIndexA, 3 * removedPartialSpringIdx, 3);
      pylimer_tools::utils::removeRows(
        net.springPartCoordinateIndexB, 3 * removedPartialSpringIdx, 3);
      pylimer_tools::utils::removeRows(
        net.springPartBoxOffset, 3 * removedPartialSpringIdx, 3);

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
                          "Removed spring found in spring indices of link " +
                            std::to_string(i) + ". Must not happen.");
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
      pylimer_tools::utils::removeRow(net.springsType, removedSpringIdx);
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

      // std::cout << "Removed springs around " << linkToReduce << " with
      // spring
      // "
      //           << removedSpringIdx << " and partial "
      //           << removedPartialSpringIdx << ", keeping " << keptSpringIdx
      //           << " and " << remainingPartialSpringIdx << std::endl;
      // std::cout << "Spring partitions sum to " << springPartitions.sum()
      //           << " for " << net.nrOfSprings
      //           << " springs, contour length before was " <<
      //           contourLengthBefore
      //           << " and is now " <<
      //           net.springsContourLength[newKeptSpringIdx]
      //           << std::endl;

      // validation
      size_t newPartialSpringIdx =
        remainingPartialSpringIdx +
        (remainingPartialSpringIdx > removedPartialSpringIdx ? -1 : 0);
      Eigen::Vector3d newDistance = this->evaluatePartialSpringDistance(
        net, u, newPartialSpringIdx, this->is2D, false);
      RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(
                        newDistance, distanceBefore, 1e-5),
                      "After merging two springs, the overall distance "
                      "is not consistent. Expected distance " +
                        std::to_string(distanceBefore) + ", but got " +
                        std::to_string(newDistance) + " for spring " +
                        std::to_string(newPartialSpringIdx) + ".");

#ifndef NDEBUG
      // check that the ordering is correct
      for (size_t i = 0;
           i < net.localToGlobalSpringIndex[newKeptSpringIdx].size();
           ++i) {
        size_t partialSpringIdx =
          net.localToGlobalSpringIndex[newKeptSpringIdx][i];
        size_t endA = net.springPartIndexA[partialSpringIdx];
        size_t endB = net.springPartIndexB[partialSpringIdx];
        std::vector<size_t> linkIndices =
          net.linkIndicesOfSprings[newKeptSpringIdx];
        assert(endA == linkIndices[i]);
        assert(endB == linkIndices[i + 1]);
      }
#endif
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
        // remove all crosslinkers that are 0- or 1-functional
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
        if (this->distanceIsWithinTolerance(distance,
                                            tolerance,
                                            net.springsContourLength[springIdx],
                                            springPartitions[springIdx]) &&
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
     * @param splitPartialSpringIdx
     * @param slipLinkIdx
     * @param alpha
     */
    size_t MEHPForceBalance::addSlipLinkToPartialSpring(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions,
      const size_t splitPartialSpringIdx,
      const size_t slipLinkIdx,
      const double oneOverSpringPartitionUpperLimit) const
    {
      INVALIDARG_EXP_IFN(
        !net.linkIsSliplink[net.springPartIndexA[splitPartialSpringIdx]] ||
          !net.linkIsSliplink[net.springPartIndexB[splitPartialSpringIdx]],
        "Require at least one part to be a cross-link.");
      const size_t newPartialSpringIdx = net.nrOfPartialSprings;
      net.nrOfPartialSprings += 1;
      const size_t relevantSpring =
        net.partialToFullSpringIndex[splitPartialSpringIdx];
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

      const size_t oldPartnerA = net.springPartIndexA[splitPartialSpringIdx];
      const size_t oldPartnerB = net.springPartIndexB[splitPartialSpringIdx];

      const Eigen::Vector3d distanceBefore =
        this->evaluatePartialSpringDistance(
          net, u, splitPartialSpringIdx, this->is2D, false);

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
      net.springPartBoxOffset.conservativeResize(3 * net.nrOfPartialSprings);
      net.partialToFullSpringIndex.conservativeResize(net.nrOfPartialSprings);
      net.partialSpringIsPartial.conservativeResize(net.nrOfPartialSprings);
      // add the new info
      net.partialSpringIsPartial[splitPartialSpringIdx] = true;
      net.partialSpringIsPartial[newPartialSpringIdx] = true;
      net.partialToFullSpringIndex[newPartialSpringIdx] = relevantSpring;

      pylimer_tools::utils::addIfNotContained(
        net.springIndicesOfLinks[slipLinkIdx], relevantSpring);

      // slightly change numbering to keep the numbering of
      // localToGlobalSpringIndex constant. I.e., we want the
      // `newPartialSpringIdx` to correspond to the spring with the cross-link
      const bool forward = !net.linkIsSliplink[oldPartnerA];
      if (forward) {
        assert(net.linkIndicesOfSprings[relevantSpring][0] == oldPartnerA);
        assert(net.localToGlobalSpringIndex[relevantSpring][0] ==
               splitPartialSpringIdx);
        // std::cout << "Case 1a" << std::endl;
        net.linkIndicesOfSprings[relevantSpring].insert(
          net.linkIndicesOfSprings[relevantSpring].begin() + 1, slipLinkIdx);
        net.localToGlobalSpringIndex[relevantSpring].insert(
          net.localToGlobalSpringIndex[relevantSpring].begin() + 1,
          newPartialSpringIdx);
      } else {
        assert(pylimer_tools::utils::last(
                 net.localToGlobalSpringIndex[relevantSpring]) ==
               splitPartialSpringIdx);
        assert(pylimer_tools::utils::last(
                 net.linkIndicesOfSprings[relevantSpring]) == oldPartnerB);
        // std::cout << "Case 2b" << std::endl;
        net.linkIndicesOfSprings[relevantSpring].insert(
          net.linkIndicesOfSprings[relevantSpring].begin() +
            (net.linkIndicesOfSprings[relevantSpring].size() - 1),
          slipLinkIdx);
        net.localToGlobalSpringIndex[relevantSpring].push_back(
          newPartialSpringIdx);
      }

      // rewire the springs
      net.springPartIndexB[splitPartialSpringIdx] = slipLinkIdx;
      net.springPartIndexA[newPartialSpringIdx] = slipLinkIdx;
      net.springPartIndexB[newPartialSpringIdx] = oldPartnerB;

      for (size_t dir = 0; dir < 3; ++dir) {
        net.springPartCoordinateIndexA[3 * splitPartialSpringIdx + dir] =
          3 * net.springPartIndexA[splitPartialSpringIdx] + dir;
        net.springPartCoordinateIndexB[3 * splitPartialSpringIdx + dir] =
          3 * net.springPartIndexB[splitPartialSpringIdx] + dir;

        net.springPartCoordinateIndexA[3 * newPartialSpringIdx + dir] =
          3 * net.springPartIndexA[newPartialSpringIdx] + dir;
        net.springPartCoordinateIndexB[3 * newPartialSpringIdx + dir] =
          3 * net.springPartIndexB[newPartialSpringIdx] + dir;
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

      // TODO: this is problematic: instead, the decision has to be made, how
      // the previous/existing offset should be split, etc.
      net.springPartBoxOffset.segment(3 * newPartialSpringIdx, 3) =
        Eigen::Vector3d::Zero();
      this->reAlignSlipLinkToImages(
        net, u, slipLinkIdx, splitPartialSpringIdx, newPartialSpringIdx);

      Eigen::Vector3d distanceAfter =
        this->evaluatePartialSpringDistance(
          net, u, splitPartialSpringIdx, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, newPartialSpringIdx, this->is2D, false);
      assert(pylimer_tools::utils::vector_approx_equal(distanceBefore,
                                                       distanceAfter));

      return newPartialSpringIdx;
    }

    /**
     * @brief Replace the two springs traversing a two-functional crosslinkers
     * with a single spring
     *
     * Also handles entanglement beads
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

          // special case: this is an entanglement bead
          if (net.oldAtomTypes[crosslinkIdx] == this->entanglementType) {
            RUNTIME_EXP_IFN(
              net.springsType[springsToMerge[0]] != this->entanglementType &&
                net.springsType[springsToMerge[1]] != this->entanglementType,
              "Got two-functional entanglement bead, expect it to be unlinked "
              "from its other entanglement bead.");
          } else {
            // two checks for two types of primary loops for the two types of
            // entanglements we could have this is the first
            std::vector<size_t> entanglementsAlong1 =
              this->getEntanglementLinkIndicesAlong(net, springsToMerge[0]);
            std::vector<size_t> entanglementsAlong2 =
              this->getEntanglementLinkIndicesAlong(net, springsToMerge[1]);

            std::sort(entanglementsAlong1.begin(), entanglementsAlong1.end());
            std::sort(entanglementsAlong2.begin(), entanglementsAlong2.end());
            // if equal, we don't merge, as that would result in a primary loop
            // with only entanglements
            if (entanglementsAlong1 == entanglementsAlong2 &&
                entanglementsAlong1.size() > 0) {
              continue;
            }
          }

          assert(springsToMerge.size() == 2);

          // second primary loop check for slip-link entanglements
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
#ifdef DEBUG_REMOVAL
            std::cout << "Merging springs " << springsToMerge[0] << " and "
                      << springsToMerge[1] << " around " << crosslinkIdx
                      << std::endl;
#endif

            // let's remove this
            // TODO: this is inefficient shit, so much data being moved
            this->mergeSprings(net,
                               displacements,
                               springPartitions,
                               springsToMerge[0],
                               springsToMerge[1],
                               crosslinkIdx);

            // this->validateNetwork(net, displacements, springPartitions);
            // std::cout << "Removing link " << crosslinkIdx << std::endl;
            this->removeLink(net, displacements, crosslinkIdx);

            // std::cout << "Removed cross-link " << crosslinkIdx << std::endl;

#ifndef NDEBUG
            this->validateNetwork(net, displacements, springPartitions);
#endif
            numRemoved += 1;
          }
          // else: TODO: decide
        }
      }
#ifndef NDEBUG
      this->validateNetwork(net, displacements, springPartitions);
#endif
      return numRemoved;
    }

    /**
     * @brief Updates the partition/parametrization of a spring around one
     * link
     *
     */
    double MEHPForceBalance::updateSpringPartition(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      Eigen::VectorXd& springPartitions, /* gives the parametrization of N */
      Eigen::VectorXd&
        oneOverSpringPartitions, /* gives the parametrization of N */
      const size_t linkIdx,
      double oneOverSpringPartitionUpperLimit,
      bool allowSlipLinksToPassEachOther) const
    {
      // std::cout << "Updating spring partition " << linkIdx << " of "
      //           << net.nrOfNodes << " / " << net.nrOfLinks << " with limit
      //           "
      //           << oneOverSpringPartitionUpperLimit << std::endl;

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
            size_t currentSpringGlobalIdx =
              net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
            size_t neighbourSpringGlobalIdx =
              net.localToGlobalSpringIndex[springIndex][partner_idx];
            // found position of this link in this spring
            // want to find the ideal value for
            // net.springPartitions[springIndex][partner_idx-1]
            // NOTE: the following is slightly problematic for primary loops!
            Eigen::Vector3d vecBack = this->evaluatePartialSpringDistanceFrom(
              net, u, currentSpringGlobalIdx, linkIdx);
            double distanceBack = (vecBack.squaredNorm());
            Eigen::Vector3d vecForward =
              this->evaluatePartialSpringDistanceFrom(
                net, u, neighbourSpringGlobalIdx, linkIdx);
            double distanceForward = vecForward.squaredNorm();
            double idealValue =
              1. / (1. + sqrt(distanceForward / distanceBack));
            if (distanceBack <= 0.0) {
              idealValue = 0.0; // TODO: really?
            }
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
              nextS >= -0.00000000000001,
              "nextS must be >= 0., got " + std::to_string(nextS) + " from " +
                std::to_string(nextS) + " and " + std::to_string(currentS) +
                ", " + std::to_string(newS) + " and " +
                std::to_string(complementaryS) + ".");
            RUNTIME_EXP_IFN(complementaryS >= -0.00000000000001,
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
                  (1.0 / (newS * N)),
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
     * @brief Loop all slip-links and move them if appropriate to other
     * springs
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
      const bool respectLoops,
      const bool moveAttempt)
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
                                          respectLoops,
                                          moveAttempt);
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
      const bool respectLoops,
      const bool moveAttempt)
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
        for (size_t linkI = 1;
             linkI < net.linkIndicesOfSprings[springIdx].size() - 1;
             ++linkI) {
          if (net.linkIndicesOfSprings[springIdx][linkI] == slipLinkIdx) {
            // found index of this slip-link.
            double partitionBeforeIdx =
              net.localToGlobalSpringIndex[springIdx][linkI - 1];
            double partitionAfterIdx =
              net.localToGlobalSpringIndex[springIdx][linkI];
            assert(this->isPartOfSpring(net, slipLinkIdx, partitionBeforeIdx));
            assert(this->isPartOfSpring(net, slipLinkIdx, partitionAfterIdx));
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
                respectLoops,
                moveAttempt);
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
                respectLoops,
                moveAttempt);
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
      const bool respectLoops,
      const bool moveAttempt) const
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
        (!net.linkIsSliplink[partnerA] || !net.linkIsSliplink[partnerB]);
      if (involvesCrosslink) {
        // first check if allowed.
        size_t slipLinkIdx = net.linkIsSliplink[partnerA] ? partnerA : partnerB;
        if ((nrOfCrosslinkSwapsAllowedPerSliplink < 0) ||
            (net.nrOfCrosslinkSwapsEndured[slipLinkIdx - net.nrOfNodes] <
             nrOfCrosslinkSwapsAllowedPerSliplink)) {
          bool didSwap = false;
          if (moveAttempt) {
            didSwap = this->swapSlipLinkWithXlinkReversibly(
              net,
              u,
              springPartitions,
              partialSpringIdx,
              oneOverSpringPartitionUpperLimit,
              respectLoops);
          } else {
            // check if energy is smaller
            size_t xlinkIdx =
              net.linkIsSliplink[partnerA] ? partnerB : partnerA;
            size_t otherRailPart = this->getOtherRailPartialSpringIdx(
              net, partialSpringIdx, slipLinkIdx);
            Eigen::Vector3d otherRailDistance =
              this->evaluatePartialSpringDistance(net, u, otherRailPart);
            Eigen::Vector3d thisRailDistance =
              this->evaluatePartialSpringDistance(net, u, partialSpringIdx);

            // make sure the vectors are in the direction of the cross-link
            if (slipLinkIdx == partnerA) {
              otherRailDistance *= -1.;
              thisRailDistance *= -1.;
            }
#ifndef NDEBUG
            this->validateNetwork();
#endif

            bool found = false;

            std::unordered_set<size_t> partialSpringIndices =
              this->getPartialSpringIndicesOfLink(net, xlinkIdx);
            for (size_t attemptedEdge : partialSpringIndices) {
              if (attemptedEdge == partialSpringIdx) {
                found = true;
                continue;
              }

              Eigen::Vector3d attemptSpringDistance =
                this->evaluatePartialSpringDistanceFrom(
                  net, u, attemptedEdge, xlinkIdx);

              // TODO: involve denominators
              // we only look at a quasi force,
              double forceEstimateBefore =
                // on the slip-link
                (-otherRailDistance + thisRailDistance).squaredNorm() +
                // and on the cross-link
                (attemptSpringDistance - thisRailDistance).squaredNorm();
              // and for the case that we did switch the slip-link onto the
              // attemptedEdge
              double forceEstimateAfter =
                // force on the slip-link (-this + this), cancel out
                (attemptSpringDistance).squaredNorm() +
                // and on the cross-link
                (-2 * thisRailDistance - otherRailDistance).squaredNorm();

              if (forceEstimateAfter < forceEstimateBefore) {
                long int newPartialSpringIdx = this->moveSlipLinkFromRailToRail(
                  net,
                  u,
                  springPartitions,
                  partialSpringIdx,
                  attemptedEdge,
                  oneOverSpringPartitionUpperLimit);
                didSwap = newPartialSpringIdx >= 0;
              }
              if (didSwap) {
                break;
              }
            }
            assert(found || didSwap);
          }
          if (didSwap) {
            net.nrOfCrosslinkSwapsEndured[slipLinkIdx - net.nrOfNodes] += 1;
          }
          return didSwap;
        }
        return false;
      } else {
        // does not involve cross-link
        // first, decide: do we attempt the move, or not?
        if (moveAttempt) {
          return this->swapSlipLinksReversibly(
            net,
            u,
            springPartitions,
            partialSpringIdx,
            oneOverSpringPartitionUpperLimit);
        } else {
          // check if energy is smaller
          size_t indexInSpring = pylimer_tools::utils::index_of(
            net.localToGlobalSpringIndex[fullSpringIdx], partialSpringIdx);
          assert(indexInSpring > 0 &&
                 indexInSpring <
                   net.localToGlobalSpringIndex[fullSpringIdx].size() - 1);
          size_t otherRailFrom =
            net.localToGlobalSpringIndex[fullSpringIdx][indexInSpring - 1];
          size_t otherRailTo =
            net.localToGlobalSpringIndex[fullSpringIdx][indexInSpring + 1];

          double thisPartialSpringDenominator =
            this->getDenominatorOfPartialSpring(
              net,
              springPartitions,
              partialSpringIdx,
              oneOverSpringPartitionUpperLimit);
          Eigen::Vector3d thisSpringDistance =
            this->evaluatePartialSpringDistance(net, u, partialSpringIdx);
          double otherRailFromDenominator = this->getDenominatorOfPartialSpring(
            net,
            springPartitions,
            otherRailFrom,
            oneOverSpringPartitionUpperLimit);
          Eigen::Vector3d otherRailFromSpringDistance =
            this->evaluatePartialSpringDistance(net, u, otherRailFrom);

          double otherRailToDenominator = this->getDenominatorOfPartialSpring(
            net,
            springPartitions,
            otherRailTo,
            oneOverSpringPartitionUpperLimit);
          Eigen::Vector3d otherRailToSpringDistance =
            this->evaluatePartialSpringDistance(net, u, otherRailTo);

          double forceEstimateBefore =
            // force estimate on `partnerA` (from)
            (-otherRailFromSpringDistance * otherRailFromDenominator +
             thisSpringDistance * thisPartialSpringDenominator)
              .squaredNorm() +
            // force estimate on `partnerB` (to)
            (-thisSpringDistance * thisPartialSpringDenominator +
             otherRailToSpringDistance * otherRailToDenominator)
              .squaredNorm();
          // force estimate if we do the mutation
          double forceEstimateAfter =
            // force estimate on `partnerA` (from)
            (thisSpringDistance * thisPartialSpringDenominator +
             (otherRailToSpringDistance + thisSpringDistance) *
               otherRailToDenominator)
              .squaredNorm() +
            // force estimate on `partnerB` (to)
            (-thisSpringDistance * thisPartialSpringDenominator -
             (thisSpringDistance + otherRailFromSpringDistance) *
               otherRailFromDenominator)
              .squaredNorm();

          if (forceEstimateAfter <= forceEstimateBefore) {
            this->swapSlipLinks(net, partialSpringIdx);
            return true;
          }

          return false;
        }
      }
#ifndef NDEBUG
      this->validateNetwork();
#endif
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
      const bool respectLoops) const
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
      // TODO: check if this is dangerous due to the differences being hidden
      // in the truncated digits
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
      // for crosslinkers, we need to take all partners of all springs into
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
      for (size_t relaxSteps = 0; relaxSteps < 3; ++relaxSteps) {
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
      const double oneOverSpringPartitionUpperLimit) const
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
      // TODO: check if this is dangerous due to the differences being hidden
      // in the truncated digits
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
      for (size_t relaxSteps = 0; relaxSteps < 3; ++relaxSteps) {
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
     * @brief Adjust the two spring's box offsets to work best with the
     * specified slip-link
     *
     * @param net the network to adjust
     * @param slipLinkIdx the slip-link around which to adjust the two springs
     * @param spring1 one of the two partial spring idx
     * @param spring2 the partial spring idx of the other spring
     */
    void MEHPForceBalance::reAlignSlipLinkToImages(
      ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const size_t slipLinkIdx,
      const size_t partialSpringIdx1,
      const size_t partialSpringIdx2) const
    {
      assert(net.springPartIndexB[partialSpringIdx1] == slipLinkIdx);
      assert(net.springPartIndexA[partialSpringIdx2] == slipLinkIdx);
      assert(net.linkIsSliplink[slipLinkIdx]);
      assert(net.partialToFullSpringIndex[partialSpringIdx1] ==
             net.partialToFullSpringIndex[partialSpringIdx2]);
      Eigen::Vector3d totalOffset =
        this->getPartialSpringBoxOffset(net, partialSpringIdx1) +
        this->getPartialSpringBoxOffset(net, partialSpringIdx2);
      Eigen::Vector3d totalDistanceBefore =
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx1, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx2, this->is2D, false);

      Eigen::Vector3d sourceCoords =
        net.coordinates.segment(3 * net.springPartIndexA[partialSpringIdx1],
                                3) +
        u.segment(3 * net.springPartIndexA[partialSpringIdx1], 3);
      Eigen::Vector3d targetCoords =
        net.coordinates.segment(3 * net.springPartIndexB[partialSpringIdx2],
                                3) +
        u.segment(3 * net.springPartIndexB[partialSpringIdx2], 3);
      Eigen::Vector3d viaCoords = net.coordinates.segment(3 * slipLinkIdx, 3) +
                                  u.segment(3 * slipLinkIdx, 3);

      // std::cout << net.springPartIndexA[partialSpringIdx1] << ": "
      //           << sourceCoords << " to "
      //           << net.springPartIndexB[partialSpringIdx2] << ": "
      //           << targetCoords << " via " << viaCoords << std::endl;

      double bestOffsetScore = -1.;
      Eigen::Vector3d bestOffset = Eigen::Vector3d::Zero();

      // ugly brute-force method to check all possible combinations (ideally,
      // more or less at least)
      Eigen::Array3i multiplicity1 =
        ((this->universe.getBox()
            .getOffset(viaCoords - sourceCoords)
            .array()
            .abs() +
          this->getPartialSpringBoxOffset(net, partialSpringIdx1)
            .array()
            .abs()) /
         this->universe.getBox().getL())
          .rint()
          .cast<int>()
          .abs();
      Eigen::Array3i multiplicity2 =
        ((this->universe.getBox().getOffset(targetCoords - viaCoords).array() +
          this->getPartialSpringBoxOffset(net, partialSpringIdx2)
            .array()
            .abs()) /
         this->universe.getBox().getL())
          .rint()
          .cast<int>()
          .abs();
      Eigen::Array3i multiplicity = multiplicity1 + multiplicity2;
      if (this->is2D) {
        multiplicity[2] = 0;
        sourceCoords[2] = 0.;
        targetCoords[2] = 0.;
        viaCoords[2] = 0.;
        totalOffset[2] = 0.;
      }
      for (int mx = std::min(0, -multiplicity[0]);
           mx <= std::max(0, multiplicity[0]);
           ++mx) {
        for (int my = std::min(0, -multiplicity[1]);
             my <= std::max(0, multiplicity[1]);
             ++my) {
          for (int mz = std::min(0, -multiplicity[2]);
               mz <= std::max(0, multiplicity[2]);
               ++mz) {
            Eigen::Vector3d currentOffset;
            currentOffset << mx * net.L[0], my * net.L[1], mz * net.L[2];

            Eigen::Vector3d vec1 = (viaCoords - sourceCoords) + currentOffset;
            Eigen::Vector3d vec2 =
              (targetCoords - viaCoords) + (totalOffset - currentOffset);
            assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
              (vec1 + vec2), totalDistanceBefore));

            double currentScore = vec1.squaredNorm() + vec2.squaredNorm();
            // std::cout << "Score: " << currentScore << " for offset "
            //           << currentOffset << std::endl;
            // std::cout << "vec 1: " << vec1 << std::endl;
            // std::cout << "vec 2: " << vec2 << std::endl;

            if (bestOffsetScore < 0 || bestOffsetScore > currentScore) {
              bestOffsetScore = currentScore;
              bestOffset = currentOffset;
            }
          }
        }
      }

      assert(bestOffsetScore >= 0.);
      net.springPartBoxOffset.segment(3 * partialSpringIdx1, 3) = bestOffset;
      net.springPartBoxOffset.segment(3 * partialSpringIdx2, 3) =
        totalOffset - bestOffset;
      Eigen::Vector3d totalDistanceNow =
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx1, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx2, this->is2D, false);
      assert(pylimer_tools::utils::vector_approx_equal(totalDistanceNow,
                                                       totalDistanceBefore));
    };

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
      const double oneOverSpringPartitionUpperLimit) const
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
#ifndef NDEBUG
      this->validateNetwork(net, u, springPartitions);
#endif

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
#ifndef NDEBUG
              try {
                this->validateNetwork(net, u, springPartitions);
              } catch (const std::runtime_error& e) {
                std::cerr << "Validation error: " << e.what() << std::endl;
                assert(false);
              }
#endif
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
      double oneOverSpringPartitionUpperLimit) const
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
      const bool respectLoops) const
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
      if (respectLoops && net.loopsOfSliplink.size() > 0) {
        // filter out the target springs that may not be a target based on the
        // involved loops
        std::erase_if(
          possibleTargetSprings,
          [&net, springIdx, involvedSlipLink](size_t springIdxToCheck) {
            // we want to be able to distinguish all associated loops
            // into being (a) the ones the slip-link is slipping on, or
            // (b) the one the second, currently non-slipping part is
            // associated with only the cases of (a) are allowed as possible
            // target springs
            bool allLoopsWithSpringIdxToCheckIncludespringIdx = false;
            for (size_t loopIdx : net.loopsOfSliplink[involvedSlipLink]) {
              bool loopContainsCheck = false;
              bool loopContainsOriginal = false;
              for (size_t lspringIdx : net.loops[loopIdx]) {
                if (lspringIdx == springIdxToCheck) {
                  loopContainsCheck = true;
                }
                if (lspringIdx == springIdx) {
                  loopContainsOriginal = true;
                }
              }
              if (loopContainsOriginal && !loopContainsCheck) {
                // YES, this is not allowed -> remove
                return true;
              }
            }
            // return allLoopsWithSpringIdxToCheckIncludespringIdx;
            return false;
          });
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

      size_t targetPartialSpringIdx =
        possibleTargetPartialSprings[(currentPartialSpringTargetIdx) %
                                     possibleTargetPartialSprings.size()];
      return this->moveSlipLinkFromRailToRail(net,
                                              u,
                                              springPartitions,
                                              partialSpringIdx,
                                              targetPartialSpringIdx,
                                              oneOverSpringPartitionUpperLimit);
    }

    /**
     * @brief
     *
     * @param net
     * @param partialSpringIdx
     */
    void MEHPForceBalance::swapSlipLinks(ForceBalanceNetwork& net,
                                         const size_t partialSpringIdx) const
    {
      const size_t linkIdx1 = net.springPartIndexA[partialSpringIdx];
      const size_t linkIdx2 = net.springPartIndexB[partialSpringIdx];
      INVALIDARG_EXP_IFN(linkIdx1 != linkIdx2,
                         "Cannot swap link with itself: got " +
                           std::to_string(linkIdx1) + " and " +
                           std::to_string(linkIdx2) + ".");
      INVALIDARG_EXP_IFN(
        net.linkIsSliplink[linkIdx1],
        "Only partial springs with only slip-links allow swapping.");
      INVALIDARG_EXP_IFN(
        net.linkIsSliplink[linkIdx2],
        "Only partial springs with only slip-links allow swapping.");
      const size_t springIdx = net.partialToFullSpringIndex[partialSpringIdx];

      long int firstPositionInSpring = pylimer_tools::utils::index_of(
        net.localToGlobalSpringIndex[springIdx], partialSpringIdx);
      assert(firstPositionInSpring > 0 &&
             firstPositionInSpring <
               net.localToGlobalSpringIndex[springIdx].size() - 1);

      // find the rest of the connectivity required for swapping
      long int otherPartialOfLinkIdx1 =
        net.localToGlobalSpringIndex[springIdx][firstPositionInSpring - 1];
      long int otherPartialOfLinkIdx2 =
        net.localToGlobalSpringIndex[springIdx][firstPositionInSpring + 1];
      assert(otherPartialOfLinkIdx1 != otherPartialOfLinkIdx2);

      const size_t unaffectedEnd1 =
        net.springPartIndexA[otherPartialOfLinkIdx1];
      const size_t unaffectedEnd2 =
        net.springPartIndexB[otherPartialOfLinkIdx2];

      // std::cout << "Swapping link " << linkIdx1 << " and " << linkIdx2
      //           << " on partial spring " << partialSpringIdx
      //           << ". Newly linked partial springs: " <<
      //           otherPartialOfLinkIdx1
      //           << " (" << unaffectedEnd1 << ") "
      //           << " and " << otherPartialOfLinkIdx2 << " (" <<
      //           unaffectedEnd2
      //           << ") " << std::endl;

      Eigen::VectorXd u = Eigen::VectorXd::Zero(net.coordinates.size());
      Eigen::Vector3d distanceBefore =
        this->evaluatePartialSpringDistance(
          net, u, otherPartialOfLinkIdx1, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, otherPartialOfLinkIdx2, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx, this->is2D, false);

      RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 >= 0,
                      "Did not find partial spring " +
                        std::to_string(partialSpringIdx) + " in spring " +
                        std::to_string(springIdx) + ".");
      RUNTIME_EXP_IFN(otherPartialOfLinkIdx2 >= 0,
                      "Did not find partial spring " +
                        std::to_string(partialSpringIdx) + " in spring " +
                        std::to_string(springIdx) + ".");
      RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 != otherPartialOfLinkIdx2,
                      "Required assumption not met.");
      RUNTIME_EXP_IFN(firstPositionInSpring <
                        net.linkIndicesOfSprings[springIdx].size() - 1,
                      "Required assumption not met.");

      // actually do the swapping
      assert(net.springPartIndexB[otherPartialOfLinkIdx1] == linkIdx1);
      assert(net.springPartIndexA[otherPartialOfLinkIdx2] == linkIdx2);
      // re-link
      net.springPartIndexB[otherPartialOfLinkIdx1] = linkIdx2;
      net.springPartIndexA[otherPartialOfLinkIdx2] = linkIdx1;
      net.springPartCoordinateIndexB.segment(3 * otherPartialOfLinkIdx1, 3) =
        Eigen::ArrayXi::LinSpaced(3, 3 * linkIdx2, 3 * linkIdx2 + 2);
      net.springPartCoordinateIndexA.segment(3 * otherPartialOfLinkIdx2, 3) =
        Eigen::ArrayXi::LinSpaced(3, 3 * linkIdx1, 3 * linkIdx1 + 2);

      // update box offset
      net.springPartBoxOffset.segment(3 * otherPartialOfLinkIdx1, 3) +=
        this->getPartialSpringBoxOffset(net, partialSpringIdx);
      net.springPartBoxOffset.segment(3 * otherPartialOfLinkIdx2, 3) +=
        this->getPartialSpringBoxOffset(net, partialSpringIdx);

      // actually change direction
      std::swap(net.springPartIndexA[partialSpringIdx],
                net.springPartIndexB[partialSpringIdx]);
      for (int dir = 0; dir < 3; ++dir) {
        std::swap(net.springPartCoordinateIndexA[3 * partialSpringIdx + dir],
                  net.springPartCoordinateIndexB[3 * partialSpringIdx + dir]);
      }
      net.springPartBoxOffset.segment(3 * partialSpringIdx, 3) *= -1.;

      assert(net.linkIndicesOfSprings[springIdx][firstPositionInSpring] ==
               linkIdx1 &&
             net.linkIndicesOfSprings[springIdx][firstPositionInSpring + 1] ==
               linkIdx2);
      net.linkIndicesOfSprings[springIdx][firstPositionInSpring] = linkIdx2;
      net.linkIndicesOfSprings[springIdx][firstPositionInSpring + 1] = linkIdx1;

      // finally, validate result
      Eigen::Vector3d distanceAfter =
        this->evaluatePartialSpringDistance(
          net, u, otherPartialOfLinkIdx1, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, otherPartialOfLinkIdx2, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx, this->is2D, false);

      assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
        distanceAfter, distanceBefore));
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
      Eigen::Vector3d forceBefore = this->getForceOn(
        net, u, springPartitions, linkIdx, oneOverSpringPartitionUpperLimit);
      // Eigen::Vector3d currentDisplacement = u.segment(3 * linkIdx, 3);
      Eigen::Vector3d objectiveDisplacement =
        Eigen::Vector3d::Zero(); // = remainingDisplacement.array();
      double objectiveDisplacementContributors = 0.0;
      bool cautionPrimaryLoop = false;

      std::unordered_set<size_t> partialSpringIndices =
        this->getPartialSpringIndicesOfLink(net, linkIdx);

      for (const size_t globalSpringIndex : partialSpringIndices) {
        assert(net.springPartIndexA[globalSpringIndex] == linkIdx ||
               net.springPartIndexB[globalSpringIndex] == linkIdx);
        if (net.springPartIndexA[globalSpringIndex] == linkIdx &&
            net.springPartIndexB[globalSpringIndex] == linkIdx) {
          // skip primary loops
          continue;
        }
        Eigen::Vector3d partialDistance =
          this->evaluatePartialSpringDistanceFrom(
            net, u, globalSpringIndex, linkIdx);
        // std::cout << "Partial distance of " << globalSpringIndex << " from
        // "
        //           << linkIdx << " to "
        //           << this->getOtherSpringIndex(net, globalSpringIndex,
        //           linkIdx)
        //           << ": " << partialDistance << std::endl;
        // add to displacement
        const double contourLengthFraction =
          springPartitions[globalSpringIndex];
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
        const double N = net.springsContourLength
                           [net.partialToFullSpringIndex[globalSpringIndex]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[globalSpringIndex],
          1.0 / (N * contourLengthFraction),
          N,
          oneOverSpringPartitionUpperLimit);
        // std::cout << "oneOverContourLengthFraction: "
        //           << oneOverContourLengthFraction << " from "
        //           << 1.0 / (N * contourLengthFraction) << " (" << N << ", "
        //           << contourLengthFraction << ")" << std::endl;
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
      }
      // take mean for displacement
      // prevent NaN from division by zero
      const double denominator = 1. / (objectiveDisplacementContributors == 0.0
                                         ? 1.0
                                         : objectiveDisplacementContributors);
      u.segment(3 * linkIdx, 3) += objectiveDisplacement * denominator;

      if (!this->assumeBoxLargeEnough) {
        Eigen::Vector3d forceAfter = this->getForceOn(
          net, u, springPartitions, linkIdx, oneOverSpringPartitionUpperLimit);

        // this is only true if we don't have "full" PBC
        assert((pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
          forceAfter, Eigen::Vector3d::Zero(), 0.01)));
        if (!pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
              forceBefore, Eigen::Vector3d::Zero(), 0.01)) {
          assert(forceBefore.squaredNorm() >= forceAfter.squaredNorm());
        }
      }

      double dist = (objectiveDisplacement * denominator).squaredNorm();
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

    /**
     * @brief Compute the stress tensor on one cross- or slip-link
     *
     * @param linkIdx
     * @param net
     * @param u
     * @param springPartitions
     * @param debugNrSpringsVisited
     * @param oneOverSpringPartitionUpperLimit
     * @return Eigen::Matrix3d
     */
    Eigen::Matrix3d MEHPForceBalance::evaluateStressOnLink(
      const size_t linkIdx,
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      Eigen::VectorXi& debugNrSpringsVisited,
      const double oneOverSpringPartitionUpperLimit) const
    {
      std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
      Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

      std::unordered_set<size_t> partialSpringIndices =
        this->getPartialSpringIndicesOfLink(net, linkIdx);

      for (const size_t globalSpringIndex : partialSpringIndices) {
        Eigen::Vector3d partialDistance =
          this->evaluatePartialSpringDistanceFrom(
            net, u, globalSpringIndex, linkIdx);
        const double contourLengthFraction =
          springPartitions[globalSpringIndex];
        const double N = net.springsContourLength
                           [net.partialToFullSpringIndex[globalSpringIndex]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[globalSpringIndex],
          1.0 / (N * contourLengthFraction),
          N,
          oneOverSpringPartitionUpperLimit);

        double multiplier = this->kappa * oneOverContourLengthFraction;

        stress += multiplier * partialDistance * partialDistance.transpose();
        debugNrSpringsVisited[globalSpringIndex] += 1;

        // also account for primary loops.
        // they may have non-zero length thanks to assuming the box is not
        // large enough...
        if (net.springPartIndexA[globalSpringIndex] ==
            net.springPartIndexB[globalSpringIndex]) {
          stress +=
            multiplier * (-partialDistance) * (-partialDistance).transpose();

          debugNrSpringsVisited[globalSpringIndex] += 1;
        }
      }

      return stress;
    }

    /**
     * @brief Compute the force acting on one cross- or slip-link
     *
     * @param linkIdx
     * @param net
     * @param u
     * @param springPartitions
     * @param debugNrSpringsVisited
     * @param oneOverSpringPartitionUpperLimit
     * @return Eigen::Vector3d
     */
    Eigen::Vector3d MEHPForceBalance::evaluateForceOnLink(
      const size_t linkIdx,
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      Eigen::VectorXi& debugNrSpringsVisited,
      const double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::Vector3d force = Eigen::Vector3d::Zero();

      std::unordered_set<size_t> partialSpringIndices =
        this->getPartialSpringIndicesOfLink(net, linkIdx);

      for (const size_t globalSpringIndex : partialSpringIndices) {
        // partial spring's force goes both ways -> is zero anyway
        // but, as it would not be included twice in the list,
        // we have to skip them
        if (net.springPartIndexA[globalSpringIndex] == linkIdx &&
            net.springPartIndexB[globalSpringIndex] == linkIdx) {
          if (debugNrSpringsVisited.size() > 0) {
            debugNrSpringsVisited[globalSpringIndex] += 2;
          }
          continue;
        }
        Eigen::Vector3d partialDistance =
          this->evaluatePartialSpringDistanceFrom(
            net, u, globalSpringIndex, linkIdx);
        const double contourLengthFraction =
          springPartitions[globalSpringIndex];
        const double N = net.springsContourLength
                           [net.partialToFullSpringIndex[globalSpringIndex]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[globalSpringIndex],
          1.0 / (N * contourLengthFraction),
          N,
          oneOverSpringPartitionUpperLimit);

        force += this->kappa * oneOverContourLengthFraction * partialDistance;
        if (debugNrSpringsVisited.size() > 0) {
          debugNrSpringsVisited[globalSpringIndex] += 1;
        }
      }

      return force;
    }

    /**
     * @brief Get a vector of all springs
     *
     * @param net
     * @param u
     * @param is2D
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance::evaluateSpringVectors(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const bool is2D,
      const bool assumeLarge) const
    {
      assert(u.size() == net.coordinates.size());

      Eigen::VectorXd springVectors =
        Eigen::VectorXd::Zero(net.nrOfSprings * 3);

      // rather than going via partial springs, we could,
      // in the case of a large enough box, directly use the spring's end
      // indices
      if (assumeLarge) {
        Eigen::VectorXd displacedCoords = (net.coordinates + u);
        springVectors = displacedCoords(net.springCoordinateIndexB) -
                        displacedCoords(net.springCoordinateIndexA);
        this->box.handlePBC(springVectors);
      } else {
        Eigen::VectorXd partialSpringVectors =
          this->evaluatePartialSpringVectors(net, u, is2D, assumeLarge);

        // CAREFUL: this assumes a certain direction
        for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
          springVectors.segment(3 * net.partialToFullSpringIndex[i], 3) +=
            partialSpringVectors.segment(3 * i, 3);
        }
      }

      // reset for 2D systems
      if (is2D && net.nrOfSprings > 0) {
        // springDistances(Eigen::seq(2, net.nrOfSprings, 3)) =
        //   Eigen::VectorXd::Zero(net.nrOfSprings);
        for (size_t i = 2; i < 3 * net.nrOfSprings; i += 3) {
          springVectors[i] = 0.0;
        }
      }

      return springVectors;
    }

    /**
     * @brief Evaluate the sum of the length of all partial springs per spring
     *
     * @param net
     * @param u
     * @param is2D
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance::evaluateSpringLengths(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const bool is2D) const
    {
      // first, the distances
      assert(u.size() == net.coordinates.size());

      Eigen::VectorXd partialSpringVectors =
        this->evaluatePartialSpringVectors(net, u);
      assert(partialSpringVectors.size() == net.nrOfPartialSprings * 3);

      Eigen::VectorXd springLengths = Eigen::VectorXd::Zero(net.nrOfSprings);
      for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
        springLengths[net.partialToFullSpringIndex[i]] +=
          partialSpringVectors.segment(3 * i, 3).norm();
      }

      return springLengths;
    }

    /**
     * @brief Count the number of intra-chain slip-links
     * i.e., slip-links that entangle a strand with itself
     *
     * @return int
     */
    int MEHPForceBalance::getNumIntraChainSlipLinks() const
    {
      int result = 0;
      for (size_t i = this->initialConfig.nrOfNodes;
           i < this->initialConfig.nrOfLinks;
           ++i) {
        if (this->initialConfig.springIndicesOfLinks[i].size() < 2) {
          result += 1;
        }
        if (this->initialConfig.springIndicesOfLinks[i].size() == 2 &&
            this->initialConfig.springIndicesOfLinks[i][0] ==
              this->initialConfig.springIndicesOfLinks[i][1]) {
          result += 1;
        }
      }

      return result;
    };

    /**
     * @brief Evaluate the vectors between the two ends of all partial springs
     *
     * @param net
     * @param u
     * @param is2D
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance::evaluatePartialSpringVectors(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const bool is2D,
      const bool assumeLarge) const
    {
      // first, the distances
      assert(u.size() == net.coordinates.size());

      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd partialDistances =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;
      if (assumeLarge) {
        this->box.handlePBC(partialDistances);
      }

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
     * FORCE BALANCE DATA ACCESS
     */
    /**
     * @brief Convert the current network back into a universe, consisting
     * only of crosslinkers
     */
    pylimer_tools::entities::Universe MEHPForceBalance::getCrosslinkerVerse()
      const
    {
      // convert nodes & springs back to a universe
      pylimer_tools::entities::Universe xlinkUniverse =
        pylimer_tools::entities::Universe(this->universe.getBox());
      std::vector<long int> ids;
      std::vector<int> types = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, crossLinkerType);
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
        // override type, since the types may be different from
        // crossLinkerType if converted with dangling chains
        types[i] = this->initialConfig.oldAtomTypes[i];
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

    /**
     * @brief Add slip-links to the current force-balance network
     *
     * @param strandIdx1
     * @param strandIdx2
     * @param x
     * @param y
     * @param z
     * @param alpha1
     * @param alpha2
     * @param loops
     * @param loopsOfSliplinks
     * @param clampAlpha
     */
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
      INVALIDARG_EXP_IFN(loopsOfSliplinks.size() == 0 ||
                           loopsOfSliplinks.size() == additionalLen,
                         "You must provide either loops for all new "
                         "slip-links, or none at all.");
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
      Eigen::VectorXd springVectorsBefore = this->evaluateSpringVectors(
        this->initialConfig, this->currentDisplacements, this->is2D, false);

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
      this->initialConfig.springPartBoxOffset.conservativeResize(
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
          size_t dividedPartialSpringIdx =
            this->initialConfig
              .localToGlobalSpringIndex[springIndex][targetIndexInSpring];
          size_t newPartialSpringIdx =
            currentNrOfPartialSprings + partialSpringsAdded;

          // for validation later
          Eigen::Vector3d distanceBefore =
            this->evaluatePartialSpringDistance(this->initialConfig,
                                                this->currentDisplacements,
                                                dividedPartialSpringIdx,
                                                this->is2D,
                                                false);

          this->initialConfig.partialSpringIsPartial[dividedPartialSpringIdx] =
            true;
          this->initialConfig.partialSpringIsPartial[newPartialSpringIdx] =
            true;

          this->initialConfig.localToGlobalSpringIndex[springIndex].insert(
            this->initialConfig.localToGlobalSpringIndex[springIndex].begin() +
              targetIndexInSpring + 1,
            newPartialSpringIdx);
          this->initialConfig.partialToFullSpringIndex[newPartialSpringIdx] =
            (springIndex);

          // adjust also the coordinates
          this->currentDisplacements.segment(3 * newNodeIdx, 3) =
            Eigen::Vector3d::Zero();
          assert(
            this->initialConfig.springPartIndexA[dividedPartialSpringIdx] ==
            springPartner1);
          this->initialConfig.springPartIndexB[dividedPartialSpringIdx] =
            newNodeIdx;
          for (size_t offset = 0; offset < 3; ++offset) {
            this->initialConfig
              .springPartCoordinateIndexB[3 * dividedPartialSpringIdx +
                                          offset] = 3 * newNodeIdx + offset;
          }

          // add the new one
          this->initialConfig.springPartIndexA[newPartialSpringIdx] =
            newNodeIdx;
          this->initialConfig.springPartIndexB[newPartialSpringIdx] =
            springPartner2;
          for (int dir = 0; dir < 3; ++dir) {
            this->initialConfig
              .springPartCoordinateIndexA[3 * newPartialSpringIdx + dir] =
              3 * newNodeIdx + dir;
            this->initialConfig
              .springPartCoordinateIndexB[3 * newPartialSpringIdx + dir] =
              3 * springPartner2 + dir;
          }

          // set box offsets
          this->initialConfig.springPartBoxOffset.segment(
            3 * newPartialSpringIdx, 3) = Eigen::Vector3d::Zero();
          this->reAlignSlipLinkToImages(this->initialConfig,
                                        this->currentDisplacements,
                                        newNodeIdx,
                                        dividedPartialSpringIdx,
                                        newPartialSpringIdx);

          this->currentSpringPartitionsVec[newPartialSpringIdx] =
            this->currentSpringPartitionsVec[dividedPartialSpringIdx] - alpha;
          RUNTIME_EXP_IFN(
            APPROX_WITHIN(this->currentSpringPartitionsVec[newPartialSpringIdx],
                          0.0,
                          1.0,
                          1e-9),
            "Spring partition must be between 0 and 1, got " +
              std::to_string(
                this->currentSpringPartitionsVec[newPartialSpringIdx]) +
              ".");
          this->currentSpringPartitionsVec[dividedPartialSpringIdx] = alpha;
          RUNTIME_EXP_IFN(
            APPROX_WITHIN(
              this->currentSpringPartitionsVec[dividedPartialSpringIdx],
              0.0,
              1.0,
              1e-9),
            "Spring partition must be between 0 and 1, got " +
              std::to_string(
                this->currentSpringPartitionsVec[newPartialSpringIdx]) +
              ".");

          this->initialConfig.linkIndicesOfSprings[springIndex].insert(
            this->initialConfig.linkIndicesOfSprings[springIndex].begin() +
              targetIndexInSpring +
              1, // + 1 to compensate for the first cross-link
            newNodeIdx);

          Eigen::Vector3d distanceAfter =
            this->evaluatePartialSpringDistance(this->initialConfig,
                                                this->currentDisplacements,
                                                dividedPartialSpringIdx,
                                                this->is2D,
                                                false) +
            this->evaluatePartialSpringDistance(this->initialConfig,
                                                this->currentDisplacements,
                                                newPartialSpringIdx,
                                                this->is2D,
                                                false);

          RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(
                            distanceAfter, distanceBefore),
                          "Expected that overall vector does not change upon "
                          "adding slip-springs");

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

#ifndef NDEBUG
      Eigen::VectorXd springVectorsAfter = this->evaluateSpringVectors(
        this->initialConfig, this->currentDisplacements, this->is2D, false);
      for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
        bool containsPrimaryLoop = false;
        for (size_t partialSpringIdx :
             this->initialConfig.localToGlobalSpringIndex[i]) {
          containsPrimaryLoop =
            containsPrimaryLoop ||
            (this->initialConfig.springPartIndexA[partialSpringIdx] ==
             this->initialConfig.springPartIndexB[partialSpringIdx]);
        }
        Eigen::Vector3d vecBefore = springVectorsBefore.segment(3 * i, 3);
        Eigen::Vector3d vecAfter = springVectorsAfter.segment(3 * i, 3);
        assert( // containsPrimaryLoop ||
          pylimer_tools::utils::vector_approx_equal(vecBefore, vecAfter));
      }
#endif

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
          r2local += this->currentSpringVectors[i * 3 + j] *
                     this->currentSpringVectors[i * 3 + j];
        }
        r2 += sqrt(r2local);
      }
      return r2 / this->initialConfig.nrOfSprings;
    }

    /**
     * @brief Get the denominator for a specified partial spring
     *
     * @param net
     * @param springPartitions
     * @param partialSpringIdx
     * @param oneOverSpringPartitionUpperLimit
     * @return double
     */
    double MEHPForceBalance::getDenominatorOfPartialSpring(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& springPartitions,
      const size_t partialSpringIdx,
      const double oneOverSpringPartitionUpperLimit) const
    {
      const double N =
        net
          .springsContourLength[net.partialToFullSpringIndex[partialSpringIdx]];
      const double fraction = springPartitions[partialSpringIdx];

      double denominator = 1. / (fraction * N);
      if (oneOverSpringPartitionUpperLimit > 0. ||
          !std::isfinite(denominator)) {
        denominator = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[partialSpringIdx],
          denominator,
          N,
          oneOverSpringPartitionUpperLimit);
      }

      assert(std::isfinite(denominator));
      return denominator;
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
          this->evaluateStressOnLink(linkIdx,
                                     net,
                                     u,
                                     springPartitions,
                                     debugNrSpringsVisited,
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
        Eigen::Matrix3d stressOnLink =
          this->evaluateStressOnLink(linkIdx,
                                     net,
                                     u,
                                     springPartitions,
                                     debugNrSpringsVisited,
                                     oneOverSpringPartitionUpperLimit);
        /* spring contribution to the overall stress tensor */
        RUNTIME_EXP_IFN(std::isfinite(stressOnLink.squaredNorm()),
                        "Got non-finite force contribution to stress tensor: " +
                          std::to_string(stressOnLink.squaredNorm()) +
                          " at link " + std::to_string(linkIdx) + "!");
        stress += stressOnLink;
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
         displacedCoords(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;

      if (this->assumeBoxLargeEnough) {
        this->box.handlePBC(relevantPartialDistancesA);
      }

      if (this->is2D) {
        for (size_t i = 2; i < relevantPartialDistancesA.size(); i += 3) {
          relevantPartialDistancesA[i] = 0.;
        }
      }

      for (size_t partialSpringIdx = 0;
           partialSpringIdx < net.nrOfPartialSprings;
           ++partialSpringIdx) {
        Eigen::Vector3d distance =
          relevantPartialDistancesA.segment(3 * partialSpringIdx, 3);
        const double contourLengthFraction = springPartitions[partialSpringIdx];
        const double N = net.springsContourLength
                           [net.partialToFullSpringIndex[partialSpringIdx]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[partialSpringIdx],
          1.0 / (N * contourLengthFraction),
          N,
          oneOverSpringPartitionUpperLimit);

        /* spring contribution to the overall stress tensor */
        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            double contribution = distance[j] * distance[k] * this->kappa *
                                  oneOverContourLengthFraction;
            RUNTIME_EXP_IFN(
              std::isfinite(contribution),
              "Got non-finite contribution to stress tensor: " +
                std::to_string(contribution) + " at coordinates " +
                std::to_string(k) + ", " + std::to_string(j) +
                " for partial spring " + std::to_string(partialSpringIdx) +
                " from distances " + std::to_string(distance[j]) + ", " +
                std::to_string(distance[k]) + " and denominator " +
                std::to_string(oneOverContourLengthFraction) + ".");
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

    Eigen::Matrix3d MEHPForceBalance::getStressTensor(
      const double oneOverSpringPartitionUpperLimit) const
    {
      std::array<std::array<double, 3>, 3> res =
        this->evaluateStressTensor(this->initialConfig,
                                   this->currentDisplacements,
                                   this->currentSpringPartitionsVec,
                                   oneOverSpringPartitionUpperLimit);

      // convert the array to an Eigen matrix
      Eigen::Matrix3d convertedRes = Eigen::Matrix3d::Zero();
      for (size_t i = 0; i < 3; ++i) {
        convertedRes.row(i) = Eigen::Vector3d::Map(res[i].data(), 3);
      }
      return convertedRes;
    }

    Eigen::Matrix3d MEHPForceBalance::getStressTensorLinkBased(
      const double oneOverSpringPartitionUpperLimit,
      const bool xlinksOnly) const
    {
      std::array<std::array<double, 3>, 3> res =
        this->evaluateStressTensorLinkBased(this->initialConfig,
                                            this->currentDisplacements,
                                            this->currentSpringPartitionsVec,
                                            oneOverSpringPartitionUpperLimit,
                                            xlinksOnly);
      Eigen::Matrix3d convertedRes = Eigen::Matrix3d::Zero();
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          convertedRes(i, j) = res[i][j];
        }
      }
      return convertedRes;
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
     * @brief Get the indices of active Nodes
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return std::vector<long int> the atom ids
     */
    std::vector<long int> MEHPForceBalance::getIndicesOfActiveNodes(
      const ForceBalanceNetwork* net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      double tolerance) const
    {
      std::vector<long int> results;
      results.reserve(net->nrOfNodes);

      // find all active springs
      Eigen::ArrayXb springIsActive =
        this->findActiveSprings(net, u, springPartitions, tolerance);

      for (size_t i = 0; i < net->nrOfNodes; i++) {
        if (net->oldAtomTypes[i] != this->entanglementType) {
          std::vector<size_t> springIndices =
            this->getInvolvedFullSpringIndices(*net, i);
          for (const size_t springIndex : springIndices) {
            if (springIsActive[springIndex]) {
              results.push_back(i);
              break;
            }
          }
        }
      }

      return results;
    };

    /**
     * @brief Get the atom ids of the active cross-links (not entanglement
     * beads/links)
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return std::vector<long int> the atom ids
     */
    std::vector<long int> MEHPForceBalance::getIdsOfActiveNodes(
      double tolerance) const
    {
      std::vector<long int> results;
      // find all active springs
      std::vector<long int> activeNodes =
        this->getIndicesOfActiveNodes(&this->initialConfig,
                                      this->currentDisplacements,
                                      this->currentSpringPartitionsVec,
                                      tolerance);

      results.reserve(activeNodes.size());

      for (long int nodeIdx : activeNodes) {
        results.push_back(this->initialConfig.oldAtomIds[nodeIdx]);
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
      Eigen::ArrayXb springIsActive = this->findActiveSprings(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfSprings; i++) {
        if (springIsActive[i]) { /* active spring */
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
      Eigen::ArrayXb partialSpringIsActive =
        this->findActivePartialSprings(tolerance);
      // translate this to the nodes
      for (size_t i = 0; i < this->initialConfig.nrOfPartialSprings; ++i) {
        if (partialSpringIsActive[i]) {
          /* active spring */
          // size_t a =
          //   this->initialConfig
          //     .springIndexA[this->initialConfig.partialToFullSpringIndex[i]];
          // size_t b =
          //   this->initialConfig
          //     .springIndexB[this->initialConfig.partialToFullSpringIndex[i]];
          size_t a = this->initialConfig.springPartIndexA[i];
          size_t b = this->initialConfig.springPartIndexB[i];
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
     * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
     * computed as phantom = N<b^2>.
     * @param nrOfChains the nr of chains to average over (can be different
     * from the nr of springs thanks to omitted free chains or primary loops)
     * @return double
     */
    double MEHPForceBalance::getGammaFactor(
      double b02,
      int nrOfChains,
      double oneOverSpringPartitionUpperLimit) const
    {
      if (b02 < 0) {
        b02 = this->defaultBondLength * this->defaultBondLength;
      }

      Eigen::VectorXd gammaFactors = this->getGammaFactors(b02);

      if (nrOfChains < 1) {
        return gammaFactors.mean();
      } else {
        return gammaFactors.sum() / static_cast<double>(nrOfChains);
      }
    }

    /**
     * @brief Get the per-(partial)-spring gamma factors
     *
     * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
     * computed as phantom = N<b^2>.
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance::getGammaFactors(
      double b02,
      double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::VectorXd springVectors = this->evaluatePartialSpringVectors(
        this->initialConfig, this->currentDisplacements);
      RUNTIME_EXP_IFN(this->currentSpringPartitionsVec.size() * 3 ==
                        springVectors.size(),
                      "Unexpected dimensions in springVectors.");

      Eigen::VectorXd gammaFactors(springVectors.size() / 3);
      const double commonDenominator = 1. / b02;
      for (size_t i = 0; i < springVectors.size() / 3; ++i) {
        double N = this->initialConfig.springsContourLength
                     [this->initialConfig.partialToFullSpringIndex[i]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          this->initialConfig.partialSpringIsPartial[i],
          1.0 / (N * this->currentSpringPartitionsVec(i)),
          N,
          oneOverSpringPartitionUpperLimit);
        gammaFactors[i] = springVectors.segment(3 * i, 3).squaredNorm() *
                          commonDenominator * oneOverContourLengthFraction;
        RUNTIME_EXP_IFN(
          std::isfinite(gammaFactors[i]),
          "Non-finite gamma factor for partial spring " + std::to_string(i) +
            ", computed from N = " + std::to_string(N) +
            ", oneOverSpringPartitionUpperLimit = " +
            std::to_string(oneOverSpringPartitionUpperLimit) +
            " and squared distance = " +
            std::to_string(springVectors.segment(3 * i, 3).squaredNorm()) +
            ".");
      }
      return gammaFactors;
    }

    /**
     * @brief Get the Weighted Partial Spring Length
     *
     * @return double
     */
    double MEHPForceBalance::getWeightedPartialSpringLength(
      const ForceBalanceNetwork& net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& springPartitions,
      size_t partialSpringIdx,
      double oneOverSpringPartitionUpperLimit) const
    {
      double N =
        net
          .springsContourLength[net.partialToFullSpringIndex[partialSpringIdx]];
      double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
        net.partialSpringIsPartial[partialSpringIdx],
        1.0 / (N * springPartitions(partialSpringIdx)),
        N,
        oneOverSpringPartitionUpperLimit);
      return this->evaluatePartialSpringDistance(net, u, partialSpringIdx)
               .norm() *
             oneOverContourLengthFraction;
    }

    /**
     * @brief Convert the universe to a network
     *
     * @param net the target network
     * @param crossLinkerType the atom type of the crossLinker
     * @return true
     * @return false
     */
    bool MEHPForceBalance::ConvertNetwork(ForceBalanceNetwork& net,
                                          const int crossLinkerType,
                                          bool remove2functionalCrosslinkers,
                                          bool removeDanglingChains)
    {
      if (remove2functionalCrosslinkers) {
        for (pylimer_tools::entities::Atom xlinker :
             this->universe.getAtomsOfType(crossLinkerType)) {
          // change type of cross-linkers with a degree <= 2 to "normal",
          // non-cross-link beads
          size_t vertexId = this->universe.getIdxByAtomId(xlinker.getId());
          if (this->universe.computeFunctionalityForVertex(vertexId) <= 2) {
            this->universe.setPropertyValue(
              vertexId, "type", crossLinkerType - 1);
          }
        }
      }

      std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
        this->universe.getChainsWithCrosslinker(crossLinkerType);

      // need to include all but dangling and free chains in order to
      // model entanglement
      size_t nrOfSprings = 0;
      std::vector<bool> useChain =
        pylimer_tools::utils::initializeWithValue<bool>(
          crossLinkerChains.size(), false);
      std::vector<long int> vertexIdToLinkIdx =
        pylimer_tools::utils::initializeWithValue<long int>(
          this->universe.getNrOfAtoms(), -1);
      size_t currentLinkIdx = 0;
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        RUNTIME_EXP_IFN(crossLinkerChains[i].getType() !=
                          pylimer_tools::entities::MoleculeType::UNDEFINED,
                        "Cross-linker chain's chain type could not be "
                        "detected. Cannot work like that.");
        if (crossLinkerChains[i].getType() ==
            pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
          assert(
            crossLinkerChains[i].getChainEnds(crossLinkerType, true).size() ==
            2);
          useChain[i] = true;
          nrOfSprings += 1;
        } else if (crossLinkerChains[i].getType() ==
                   pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
          // when omitting f=2 cross-links, it's possible that we end up with
          // "free" primary loops – let's not use those
          if (crossLinkerChains[i].getAtomsOfType(crossLinkerType).size() > 0) {
            useChain[i] = true;
            nrOfSprings += 1;
          }
        } else if (!removeDanglingChains &&
                   crossLinkerChains[i].getType() ==
                     pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) {
          std::vector<pylimer_tools::entities::Atom> endAtoms =
            crossLinkerChains[i].getAtomsOfDegree(1);
          RUNTIME_EXP_IFN(endAtoms.size() == 2,
                          "Expected a dangling chain to have two ends, got " +
                            std::to_string(endAtoms.size()) + ".");
          // cannot assert this without entanglement types being set
          // assert(XOR((endAtoms[0].getType() == crossLinkerType),
          //            endAtoms[1].getType() == crossLinkerType));

          useChain[i] = true;
          nrOfSprings += 1;
        }

        if (useChain[i]) {
          std::vector<pylimer_tools::entities::Atom> endAtoms =
            crossLinkerChains[i].getChainEnds(crossLinkerType);
          for (const pylimer_tools::entities::Atom& endAtom : endAtoms) {
            size_t atomVertexIdx =
              this->universe.getIdxByAtomId(endAtom.getId());
            if (vertexIdToLinkIdx[atomVertexIdx] == -1) {
              vertexIdToLinkIdx[atomVertexIdx] = currentLinkIdx;
              currentLinkIdx += 1;
            }
          }
        }
      }

      size_t nrOfXlinks = currentLinkIdx;

      // crossLinkerUniverse.simplify();
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
      net.oldAtomTypes = Eigen::ArrayXi::Zero(net.nrOfLinks);
      net.linkIsSliplink = Eigen::ArrayXb::Constant(net.nrOfLinks, false);
      net.springIndicesOfLinks.reserve(net.nrOfLinks);
      for (size_t i = 0; i < net.nrOfLinks; ++i) {
        net.springIndicesOfLinks.push_back(std::vector<size_t>());
      }
      net.linkIndicesOfSprings.reserve(net.nrOfSprings);
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        net.linkIndicesOfSprings.push_back(std::vector<size_t>());
      }
      net.partialToFullSpringIndex = Eigen::ArrayXi(net.nrOfPartialSprings);
      this->currentSpringPartitionsVec = Eigen::VectorXd::Ones(net.nrOfSprings);
      net.springIndexA = Eigen::ArrayXi::Zero(net.nrOfSprings);
      net.springIndexB = Eigen::ArrayXi::Zero(net.nrOfSprings);
      net.springCoordinateIndexA = Eigen::ArrayXi::Zero(3 * net.nrOfSprings);
      net.springCoordinateIndexB = Eigen::ArrayXi::Zero(3 * net.nrOfSprings);
      net.springPartBoxOffset = Eigen::VectorXd::Zero(3 * net.nrOfSprings);
      net.springIsActive = Eigen::ArrayXb::Constant(net.nrOfSprings, false);
      net.springsContourLength = Eigen::VectorXd::Zero(net.nrOfSprings);
      net.springsType = Eigen::ArrayXi::Zero(net.nrOfSprings);
      net.oldAtomIdToSpringIndex.reserve(this->universe.getNrOfAtoms());
      net.springToMoleculeIds.reserve(nrOfSprings);
      net.partialSpringIsPartial =
        Eigen::ArrayXb::Constant(net.nrOfSprings, false);

      // convert (cross-linker-)beads
      std::map<int, int> atomIdToNode;
      for (size_t i = 0; i < vertexIdToLinkIdx.size(); ++i) {
        if (vertexIdToLinkIdx[i] != -1) {
          size_t linkIdx = vertexIdToLinkIdx[i];
          pylimer_tools::entities::Atom atom =
            this->universe.getAtomByVertexIdx(i);
          atomIdToNode[atom.getId()] = linkIdx;
          net.oldAtomIds[linkIdx] = atom.getId();
          net.oldAtomTypes[linkIdx] = atom.getType();
          Eigen::Vector3d coords = atom.getCoordinates();
          this->universe.getBox().handlePBC(coords);
          net.coordinates.segment(3 * linkIdx, 3) = coords;
        }
      }

      // convert springs
      size_t spring_idx = 0;
      Eigen::Vector3d expectedDistance = Eigen::Vector3d::Zero();
      // net.connectivityToSpringIndex.reserve(nrOfSprings);
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        if (!useChain[i]) {
          continue;
        }
        std::vector<pylimer_tools::entities::Atom> xlinkersOfChain =
          crossLinkerChains[i].getAtomsOfType(crossLinkerType);
        std::vector<pylimer_tools::entities::Atom> chainEnds =
          crossLinkerChains[i].getChainEnds(crossLinkerType, true);
        RUNTIME_EXP_IFN(
          chainEnds.size() == 2,
          "Expected two chain ends when converting structure. Got " +
            std::to_string(chainEnds.size()) + ".");
        long int atomIdFrom = chainEnds[0].getId();
        long int atomIdTo = chainEnds[1].getId();
        bool addChain = false;
        if (crossLinkerChains[i].getType() ==
            pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
          addChain = true;

          // spring contour length = nr of bonds between two cross-linkers
          net.springsContourLength[spring_idx] =
            crossLinkerChains[i].getNrOfAtoms() - 1;
        } else if (crossLinkerChains[i].getType() ==
                   pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
          addChain = true;

          net.springsContourLength[spring_idx] =
            crossLinkerChains[i].getNrOfAtoms();
          if (xlinkersOfChain.size() == 2) {
            net.springsContourLength[spring_idx] =
              crossLinkerChains[i].getNrOfAtoms() - 1;
          }
        } else if (crossLinkerChains[i].getType() ==
                     pylimer_tools::entities::MoleculeType::DANGLING_CHAIN &&
                   !removeDanglingChains) {
          net.springsContourLength[spring_idx] =
            crossLinkerChains[i].getNrOfAtoms() - 1;
          addChain = true;
        }
        auto bondTypes = crossLinkerChains[i].getBonds()["bond_type"];
        net.springsType[spring_idx] = MEAN(bondTypes);
        assert(addChain);

        if (addChain) {
          long int nodeIdxFrom = atomIdToNode.at(atomIdFrom);
          long int nodeIdxTo = atomIdToNode.at(atomIdTo);
          // if (nodeIdxFrom > nodeIdxTo) {
          //   std::swap(nodeIdxFrom, nodeIdxTo);
          //   std::swap(atomIdFrom, atomIdTo);
          // }

          net.springToMoleculeIds.push_back(i);
          std::vector<pylimer_tools::entities::Atom> allChainAtoms =
            crossLinkerChains[i].getAtoms();
          for (const pylimer_tools::entities::Atom& a : allChainAtoms) {
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

          expectedDistance = crossLinkerChains[i].getOverallBondSumFromTo(
            atomIdFrom, atomIdTo, crossLinkerType, true);
          Eigen::Vector3d actualDistance =
            net.coordinates.segment(3 * nodeIdxTo, 3) -
            net.coordinates.segment(3 * nodeIdxFrom, 3);
          net.springPartBoxOffset.segment(3 * spring_idx, 3) =
            expectedDistance - actualDistance;
          assert(this->universe.getBox().isValidOffset(expectedDistance -
                                                       actualDistance));

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
      RUNTIME_EXP_IFN(net.springsType.size() == net.nrOfSprings,
                      "Invalid size of springs types, got " +
                        std::to_string(net.springsType.size()) + " for " +
                        std::to_string(net.nrOfSprings) + " springs.");
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
      RUNTIME_EXP_IFN(net.oldAtomTypes.size() == net.nrOfNodes,
                      "Invalid size of old atom types");
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
                        "Full springs must consist of crosslinkers only.");
        RUNTIME_EXP_IFN(net.springIndexB.maxCoeff() < net.nrOfNodes,
                        "Full springs must consist of crosslinkers only.");
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
          "Expected slip-links to come sequentially after crosslinkers.");
        std::vector<size_t> thisLinksSprings =
          net.springIndicesOfLinks[link_idx];
        std::sort(thisLinksSprings.begin(), thisLinksSprings.end());
        auto last =
          std::unique(thisLinksSprings.begin(), thisLinksSprings.end());
        RUNTIME_EXP_IFN(last == thisLinksSprings.end(),
                        "Expect each link to only have one back-link to the "
                        "springs, found back-links " +
                          pylimer_tools::utils::join(thisLinksSprings.begin(),
                                                     thisLinksSprings.end(),
                                                     std::string("_")) +
                          " for link " + std::to_string(link_idx) + ".");
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
      Eigen::ArrayXi nrOfMentions = Eigen::ArrayXi::Zero(net.nrOfLinks);
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        RUNTIME_EXP_IFN(net.linkIndicesOfSprings[i].size() >= 2,
                        "Each spring requires at least two links, got " +
                          std::to_string(net.linkIndicesOfSprings[i].size()) +
                          " at i = " + std::to_string(i) + ".");
        RUNTIME_EXP_IFN(net.springsContourLength[i] > 0,
                        "Unexpected spring contour length, got " +
                          std::to_string(net.springsContourLength[i]) +
                          " for spring " + std::to_string(i) + ".");
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
              net.springPartIndexB[partialSpringIdx] == partner1)),
            "Expect linkIndicesOfSprings and localToGlobalSpringIndex "
            "ordering to correspond. Got partner0 = " +
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
        // the following is not guaranteed anymore with the removal of links
        // while running RUNTIME_EXP_IFN(
        //   net.linkIndicesOfSprings[i][0] <=
        //     net.linkIndicesOfSprings[i][net.linkIndicesOfSprings[i].size()
        //     - 1],
        //   "Springs must have increasing end-point indices");
        std::vector<size_t> links = net.linkIndicesOfSprings[i];
        for (size_t j = 0; j < links.size(); ++j) {
          size_t link_idx = links[j];
          nrOfMentions[link_idx] += 1;
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
      for (size_t i = net.nrOfNodes; i < net.nrOfLinks; ++i) {
        RUNTIME_EXP_IFN(nrOfMentions[i] == 2,
                        "Expect each slip-link to be mentioned twice in the "
                        "links-of-springs mapping, but " +
                          std::to_string(i) + " was mentioned " +
                          std::to_string(nrOfMentions[i]) + " times.");
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
            net.springIndexA[fullIdx] == partialEndA,
            "Expect mapping of springs to work: " +
              std::to_string(partialEndA) +
              " is a cross-link, yet not part of the two ends of spring " +
              std::to_string(fullIdx) + ", where we have " +
              std::to_string(net.springIndexA[fullIdx]) + " and " +
              std::to_string(net.springIndexB[fullIdx]) + ".");
        }
        if (!net.linkIsSliplink[partialEndB]) {
          RUNTIME_EXP_IFN(
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
        for (int dir = 0; dir < 3; ++dir) {
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
            "Spring part index and coordinate index must match. Got " +
              std::to_string(net.springPartCoordinateIndexB[3 * i + dir]) +
              " but expected " + std::to_string(3 * partialEndB + dir) +
              " with dir = " + std::to_string(dir) + ".");
        }
      }

      /**
       * Check that we do not have any nan or inf values in our vectors
       */
      for (long int coordI = 0; coordI < net.coordinates.size(); coordI++) {
        RUNTIME_EXP_IFN(std::isfinite(net.coordinates[coordI]),
                        "Coordinate component " + std::to_string(coordI) +
                          " must be finite, got " +
                          std::to_string(net.coordinates[coordI]) + ".");
        RUNTIME_EXP_IFN(std::isfinite(u[coordI]),
                        "Displacement component " + std::to_string(coordI) +
                          " must be finite, got " + std::to_string(u[coordI]) +
                          ".");
      }
      for (int dir = 0; dir < 3; ++dir) {
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

      /**
       * Validate additional entanglement-atom specific data that might not
       * apply
       */
      for (size_t linkIdx = 0; linkIdx < net.nrOfNodes; ++linkIdx) {
        if (net.oldAtomTypes[linkIdx] == this->entanglementType) {
          RUNTIME_EXP_IFN(
            net.springIndicesOfLinks[linkIdx].size() <= 3,
            "Expect each entanglement atom to have up to three "
            "springs maximum, got " +
              std::to_string(net.springIndicesOfLinks[linkIdx].size()) +
              " for link " + std::to_string(linkIdx) + ".");
          size_t nEntanglementPartners = 0;
          size_t nShortEntanglementPartners = 0;
          for (size_t springIdx : net.springIndicesOfLinks[linkIdx]) {
            size_t partnerIdx = this->getOtherEnd(net, springIdx, linkIdx);
            if (partnerIdx <= net.nrOfNodes &&
                net.oldAtomTypes[partnerIdx] == this->entanglementType) {
              nEntanglementPartners += 1;
              nShortEntanglementPartners += static_cast<int>(
                net.springsType[springIdx] == this->entanglementType);
            }
          }
          if (net.springIndicesOfLinks[linkIdx].size() == 3) {
            RUNTIME_EXP_IFN(nShortEntanglementPartners == 1,
                            "Each entanglement atom must have exactly one "
                            "entanglement partner. Got " +
                              std::to_string(nShortEntanglementPartners) +
                              " out of " +
                              std::to_string(nEntanglementPartners) +
                              " for link + " + std::to_string(linkIdx) + ".");
          } else {
            RUNTIME_EXP_IFN(
              net.springIndicesOfLinks[linkIdx].size() == 2,
              "Expect each entanglement atom without "
              "entanglement partner to be twofunctional link, got "
              "functionality " +
                std::to_string(net.springIndicesOfLinks[linkIdx].size()) +
                " for link " + std::to_string(linkIdx) + " (originally " +
                std::to_string(net.oldAtomIds[linkIdx]) + ").");
            RUNTIME_EXP_IFN(nShortEntanglementPartners == 0,
                            "2-functional entanglement atoms are only allowed "
                            "if the entanglement springs has been removed.")
          }
        }
      }
      for (size_t springIdx = 0; springIdx < net.nrOfSprings; ++springIdx) {
        if (net.springsType[springIdx] == this->entanglementType) {
          RUNTIME_EXP_IFN(
            net.springsContourLength[springIdx] == 1,
            "Entanglement springs must have contour length 1. Got " +
              std::to_string(net.springsContourLength[springIdx]) +
              " for spring " + std::to_string(springIdx) + ".");
          RUNTIME_EXP_IFN(
            net.oldAtomTypes[net.springIndexA[springIdx]] ==
                this->entanglementType &&
              net.oldAtomTypes[net.springIndexB[springIdx]] ==
                this->entanglementType,
            "Expect each entanglement spring to connect two entanglement "
            "atoms. Got " +
              std::to_string(net.oldAtomTypes[net.springIndexA[springIdx]]) +
              " and " +
              std::to_string(net.oldAtomTypes[net.springIndexB[springIdx]]) +
              " for spring " + std::to_string(springIdx) + ".");
        }
      }

      // std::cout << "Validation passed." << std::endl;
      return true;
    }
  }
}
}
