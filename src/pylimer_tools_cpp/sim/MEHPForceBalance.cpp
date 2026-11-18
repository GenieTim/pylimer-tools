#include "MEHPForceBalance.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../topo/EntanglementDetector.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// #ifndef NDEBUG
// #define DEBUG_REMOVAL
// #endif

namespace pylimer_tools::sim::mehp {
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
void
MEHPForceBalance::runForceRelaxation(
  const long int maxNrOfSteps, // default: 10000
  const double xtol,
  const double initialResidualToUse,
  const StructureSimplificationMode simplificationMode,
  const double inactiveRemovalCutoff,
  const bool doInnerIterations,
  const LinkSwappingMode allowSlipLinksToPassEachOther,
  const int swappingFrequency,
  const double oneOverSpringPartitionUpperLimit,
  const int nrOfCrosslinkSwapsAllowedPerSliplink,
  const bool disableSlipping,
  const std::function<bool()>& shouldInterrupt,
  const std::function<void()>& cleanupInterrupt)
{
  RUNTIME_EXP_IFN(this->validateNetwork(), "Network validation failed.");
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

  if (this->getNrOfSprings() == 0) {
    return;
  }

  /* array allocation */
  std::vector<Eigen::ArrayXi> independentVertexSets;
  double maxDistanceMoved = 0.0;

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
  std::cout << "Starting force balance procedure "
            << "with " << initialResidual
            << " as initial residual, got requested "
            << initialResidualToUse
            // "with " << independentVertexSets.size() << "vertex sets."
            << std::endl;
  std::cout << "Swapping mode is " << allowSlipLinksToPassEachOther << " ("
            << LinkSwappingModeNames[allowSlipLinksToPassEachOther] << ")"
            << " while simplification mode is " << simplificationMode << " ("
            << StructureSimplificationModeNames[simplificationMode] << ")"
            << std::endl;
  std::cout << "Using oneOverSpringPartitionUpperLimit = "
            << oneOverSpringPartitionUpperLimit << std::endl;
  double currentResidual = initialResidual;
  size_t iterationsDone = 0;

  std::vector<std::vector<size_t>> partialSpringIndicesOfLinks(
    this->initialConfig.nrOfLinks);
  const std::function<void()> refreshPartialSpringIndices = [&]() {
    partialSpringIndicesOfLinks.resize(this->initialConfig.nrOfLinks);
    for (size_t i = 0; i < this->initialConfig.nrOfLinks; ++i) {
      partialSpringIndicesOfLinks[i] =
        this->getPartialSpringIndicesOfLink(this->initialConfig, i);
    }
  };
  refreshPartialSpringIndices();

  this->prepareAllOutputs();

  // actual loop
  bool wasInterrupted = false;
  bool iterateForDisplacements =
    !(disableSlipping ||
      (this->initialConfig.nrOfNodes == this->initialConfig.nrOfLinks));
  do {
    if (allowSlipLinksToPassEachOther != LinkSwappingMode::NO_SWAPPING) {
      if (swappingFrequency > 0 &&
          (iterationsDone % static_cast<size_t>(swappingFrequency)) == 0) {
        if (allowSlipLinksToPassEachOther == LinkSwappingMode::SLIPLINKS_ONLY) {
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
                   allowSlipLinksToPassEachOther == LinkSwappingMode::ALL_MC) {
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
        refreshPartialSpringIndices();
      }
    }
    maxDistanceMoved = 0.0;

    // place slip-link
    if (!disableSlipping) {
      for (size_t link_idx = this->initialConfig.nrOfNodes;
           link_idx < this->initialConfig.nrOfLinks;
           ++link_idx) {
        assert(this->initialConfig
                 .linkIsSliplink[static_cast<Eigen::Index>(link_idx)]);
        // std::cout << "Handling " << link_idx << " of " << net.nrOfNodes
        //           << " / " << net.nrOfLinks << std::endl;

        // std::cout << "Still handling " << link_idx << " of " <<
        // net.nrOfNodes
        //           << " / " << net.nrOfLinks << std::endl;
        int innerIterationsDone = 0;
        do {
          this->updateSpringPartition(this->initialConfig,
                                      this->currentDisplacements,
                                      this->currentSpringPartitionsVec,
                                      oneOverSpringPartitions,
                                      link_idx,
                                      oneOverSpringPartitionUpperLimit,
                                      allowSlipLinksToPassEachOther);

          double displacementDone =
            this->displaceToMeanPosition(this->initialConfig,
                                         this->currentDisplacements,
                                         this->currentSpringPartitionsVec,
                                         link_idx,
                                         oneOverSpringPartitionUpperLimit,
                                         partialSpringIndicesOfLinks[link_idx]);
          maxDistanceMoved = std::max(maxDistanceMoved, displacementDone);
          innerIterationsDone += 1;
        } while (doInnerIterations && innerIterationsDone < 50);
      }
      oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(this->initialConfig,
                                             this->currentSpringPartitionsVec,
                                             oneOverSpringPartitionUpperLimit);
    }

    // place crosslinkers
    if (iterateForDisplacements) {
      for (size_t link_idx = 0; link_idx < this->initialConfig.nrOfNodes;
           ++link_idx) {
        assert(!this->initialConfig
                  .linkIsSliplink[static_cast<Eigen::Index>(link_idx)]);
        const double distanceMoved =
          this->displaceToMeanPosition(this->initialConfig,
                                       this->currentDisplacements,
                                       this->currentSpringPartitionsVec,
                                       link_idx,
                                       oneOverSpringPartitionUpperLimit,
                                       partialSpringIndicesOfLinks[link_idx]);
        if (distanceMoved > maxDistanceMoved) {
          maxDistanceMoved = distanceMoved;
        }
      }
    } else {
      maxDistanceMoved =
        std::max(maxDistanceMoved,
                 this->displaceToMeanPosition(this->initialConfig,
                                              this->currentDisplacements,
                                              oneOverSpringPartitions.array()));
      RUNTIME_EXP_IFN(
        std::isfinite(maxDistanceMoved),
        "Something went wrong: The last distance moved is not finite but " +
          std::to_string(maxDistanceMoved) + " in iteration " +
          std::to_string(iterationsDone + 1) + ".");
    }

    currentResidual =
      this->getDisplacementResidualNormFor(this->initialConfig,
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

    iterationsDone += 1;
    if (iterationsDone % static_cast<size_t>(this->simplificationFrequency) ==
        0) {
      const size_t numBroken =
        this->breakTooLongSprings(this->initialConfig,
                                  this->currentDisplacements,
                                  this->currentSpringPartitionsVec);
      if (numBroken > 0) {
        refreshPartialSpringIndices();
      }
      size_t nRemoved = 0;
      size_t nRemovedThisLoop = 0;

      do {
#ifndef NDEBUG
        assert(this->validateNetwork());
#endif
        nRemovedThisLoop = 0;
        if (simplificationMode == StructureSimplificationMode::INACTIVE_ONLY ||
            simplificationMode ==
              StructureSimplificationMode::INACTIVE_THEN_X2F) {
#ifdef DEBUG_REMOVAL
          std::cout << "Checking and possibly removing inactive crosslinks"
                    << std::endl;
#endif
          nRemovedThisLoop +=
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
            simplificationMode ==
              StructureSimplificationMode::INACTIVE_THEN_X2F) {
#ifdef DEBUG_REMOVAL
          std::cout << "Checking and possibly removing crosslinks with f = 2"
                    << std::endl;
#endif
          nRemovedThisLoop += this->removeTwofunctionalCrosslinks(
            this->initialConfig,
            this->currentDisplacements,
            this->currentSpringPartitionsVec);
          this->initialConfig.meanSpringContourLength =
            this->initialConfig.springsContourLength.size() > 0
              ? this->initialConfig.springsContourLength.mean()
              : 0.;
        }
        if (simplificationMode ==
            StructureSimplificationMode::X1F_X2F_THEN_INACTIVE) {
#ifdef DEBUG_REMOVAL
          std::cout << "Checking and possibly removing crosslinks and "
                       "springs, Andrei's way"
                    << std::endl;
#endif
          nRemovedThisLoop +=
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

        nRemoved += nRemovedThisLoop;
      } while (nRemovedThisLoop > 0);
      // after removal, the residual changed, might even have increased
      // beyond initial
      // -> reset previous and current to prevent change to iterative
      // displacement
      if (nRemoved > 0) {
        oneOverSpringPartitions = this->assembleOneOverSpringPartition(
          this->initialConfig,
          this->currentSpringPartitionsVec,
          oneOverSpringPartitionUpperLimit);
        refreshPartialSpringIndices();
      }
    }
    this->handleOutput(static_cast<long int>(iterationsDone));

    if (shouldInterrupt()) {
      wasInterrupted = true;
      break;
    }
  } while (currentResidual / initialResidual > xtol &&
           static_cast<long int>(iterationsDone) < maxNrOfSteps &&
           this->initialConfig.nrOfSprings > 0);

  // finish up
  this->closeAllOutputs();

  // query solution & exit reason
  this->exitReason = (static_cast<long int>(iterationsDone) == maxNrOfSteps)
                       ? ExitReason::MAX_STEPS
                       : ExitReason::X_TOLERANCE;
  this->nrOfStepsDone += iterationsDone;
  std::cout << iterationsDone << " steps done. "
            << "Last max distance moved: " << maxDistanceMoved << ". "
            << "Current residual: " << currentResidual << ". "
            << "Initial residual: " << initialResidual << ". " << std::endl;

  assert(this->currentDisplacements.size() ==
         3 * static_cast<Eigen::Index>(this->initialConfig.nrOfLinks));
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
double
MEHPForceBalance::getDisplacementResidualNorm(
  const double oneOverSpringPartitionUpperLimit) const
{
  const Eigen::VectorXd oneOverSpringPartitions =
    this->assembleOneOverSpringPartition(this->initialConfig,
                                         this->currentSpringPartitionsVec,
                                         oneOverSpringPartitionUpperLimit);
  const Eigen::VectorXd displacements = this->currentDisplacements;
  return this->getDisplacementResidualNormFor(
    this->initialConfig, displacements, oneOverSpringPartitions);
}

/**
 * @brief Compute the displacement residual norm for a specific
 * configuration
 *
 * @param net
 * @param u
 * @param springPartitions
 * @param oneOverSpringPartitionUpperLimit
 * @return double
 */
double
MEHPForceBalance::getDisplacementResidualNormFor(
  const ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  const Eigen::VectorXd& springPartitions,
  const double oneOverSpringPartitionUpperLimit) const
{
  const Eigen::ArrayXd oneOverSpringPartitions =
    this
      ->assembleOneOverSpringPartition(
        net, springPartitions, oneOverSpringPartitionUpperLimit)
      .array();

  const Eigen::ArrayXd loopPartialSpringEliminator =
    (net.springPartCoordinateIndexA != net.springPartCoordinateIndexB)
      .cast<double>();
  Eigen::ArrayXd forces =
    Eigen::ArrayXd::Zero(3 * static_cast<Eigen::Index>(net.nrOfLinks));
  const Eigen::ArrayXd distances =
    this
      ->evaluatePartialSpringVectors(
        net, u, this->is2D, this->assumeBoxLargeEnough)
      .array();
  forces(net.springPartCoordinateIndexA) +=
    (this->kappa * oneOverSpringPartitions * distances *
     loopPartialSpringEliminator);
  forces(net.springPartCoordinateIndexB) -=
    (this->kappa * oneOverSpringPartitions * distances *
     loopPartialSpringEliminator);

  // #ifndef NDEBUG
  //       Eigen::ArrayXd forces2 = Eigen::ArrayXd::Zero(3 * net.nrOfLinks);
  //       Eigen::VectorXi debugNrSpringsVisited =
  //         Eigen::VectorXi::Zero(net.nrOfPartialSprings);
  //       for (size_t i = 0; i < net.nrOfLinks; ++i) {
  //         forces2.segment(3 * i, 3) =
  //           this
  //             ->evaluateForceOnLink(i,
  //                                   net,
  //                                   u,
  //                                   springPartitions,
  //                                   debugNrSpringsVisited,
  //                                   oneOverSpringPartitionUpperLimit)
  //             .array();
  //       }
  //       assert((debugNrSpringsVisited.array() == 2).all());
  //       double squareN1 = forces.matrix().squaredNorm();
  //       double squareN2 = forces2.matrix().squaredNorm();
  //       assert(APPROX_EQUAL(squareN1, squareN2, 1e-9));
  //       assert(
  //         pylimer_tools::utils::vector_approx_equal(forces, forces2,
  //         1e-9, true));
  // #endif

  return forces.matrix().squaredNorm();
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
double
MEHPForceBalance::getDisplacementResidualNormFor(
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
    for (Eigen::Index i = 2; i < relevantPartialDistances.size(); i += 3) {
      relevantPartialDistances[i] = 0.;
    }
  }

#ifndef NDEBUG
  for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
    Eigen::Vector3d dist = this->evaluatePartialSpringDistance(net, u, i);
    Eigen::Vector3d comparison =
      relevantPartialDistances.segment(3 * static_cast<Eigen::Index>(i), 3);
    assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
      dist, comparison, 1e-9));
  }
#endif

  assert(relevantPartialDistances.size() == oneOverSpringPartitions.size());
  const Eigen::VectorXd partialDistancesOverSpringPartitions =
    (relevantPartialDistances.array() * oneOverSpringPartitions.array())
      .matrix();

  Eigen::VectorXd overallForces =
    Eigen::VectorXd::Zero(3 * static_cast<Eigen::Index>(net.nrOfLinks));
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
Eigen::VectorXd
MEHPForceBalance::assembleOneOverSpringPartition(
  const ForceBalanceNetwork& net,
  const Eigen::VectorXd& springPartitions0,
  const double oneOverSpringPartitionUpperLimit)
{
  INVALIDARG_EXP_IFN(
    springPartitions0.size() ==
      static_cast<Eigen::Index>(net.nrOfPartialSprings),
    "Spring partitions must have the size of the nr of springs");
  Eigen::VectorXd oneOverSpringPartitions =
    Eigen::VectorXd(3 * net.nrOfPartialSprings);

  Eigen::ArrayXd primaryLoopCorrectionMultiplier =
    (net.springPartIndexA != net.springPartIndexB)
      .cast<double>(); // 0.0 for equal = primary loop, 1.0 otherwise

  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(net.nrOfPartialSprings);
       ++i) {
    const double N = net.springsContourLength[net.partialToFullSpringIndex[i]];
    const double contourLengthFraction = springPartitions0[i];
    const double valueToSet =
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
 * @param results
 * @param linkIdx
 * @return void
 */
void
MEHPForceBalance::setSpringpartitionIndicesOfSliplink(
  std::vector<size_t>& results,
  const ForceBalanceNetwork& net,
  const size_t linkIdx)
{
  INVALIDARG_EXP_IFN(
    linkIdx < net.nrOfLinks,
    "Cannot set spring partition of index higher than nr. of links.");
  INVALIDARG_EXP_IFN(net.linkIsSliplink[static_cast<Eigen::Index>(linkIdx)],
                     "Link must be slip-link");
  const std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
  size_t indexIndex = 0;
  while (results.size() < 4) {
    results.push_back(0);
  }
  for (const size_t springIndex : springIndices) {
    std::vector<size_t> springsPartners = net.linkIndicesOfSprings[springIndex];
    for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
         ++partner_idx) {
      if (springsPartners[partner_idx] == linkIdx) {
        RUNTIME_EXP_IFN(
          indexIndex < 4,
          "Expect spring partitions indices of link not to exceed 4.");
        const size_t currentSpringGlobalIdx =
          net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
        const size_t neighbourSpringGlobalIdx =
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
 * @brief Remove double listed springs from crosslinkers
 *
 * @param net
 */
void
MEHPForceBalance::removeDuplicateListedSpringsFromLinks(
  ForceBalanceNetwork& net) const
{
  for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
    this->removeDuplicateListedSpringsFromLink(net, linkIdx);
  }

#ifndef NDEBUG
  this->validateNetwork();
#endif
}

void
MEHPForceBalance::removeDuplicateListedSpringsFromLink(
  ForceBalanceNetwork& net,
  size_t linkIdx,
  const bool allowOnEntanglement) const
{
  INVALIDARG_EXP_IFN(linkIdx < net.nrOfLinks,
                     "Cannot remove duplicate spring indices of index "
                     "higher than nr. of links.");
  // remove duplicate mentions of the same spring index
  std::ranges::sort(net.springIndicesOfLinks[linkIdx]);
  const auto last =
    std::ranges::unique(net.springIndicesOfLinks[linkIdx]).begin();
  if (last != net.springIndicesOfLinks[linkIdx].end()) {
#ifdef DEBUG_REMOVAL
    std::cout << "Removed duplicate spring indices from link " << linkIdx
              << std::endl;
#endif
    if (!allowOnEntanglement) {
      const Eigen::Index linkIdxE = static_cast<Eigen::Index>(linkIdx);
      if (linkIdx < net.nrOfNodes) {
        RUNTIME_EXP_IFN(
          net.oldAtomTypes[linkIdxE] != this->entanglementType,
          "Require entanglement beads to not form primary loops. Link " +
            std::to_string(linkIdx) + " is slip-link " +
            std::to_string(net.linkIsSliplink[linkIdxE]) +
            " with old atom type" + std::to_string(net.oldAtomTypes[linkIdxE]) +
            ".");
      } else {
        RUNTIME_EXP_IFN(
          !net.linkIsSliplink[linkIdxE],
          "Require entanglement beads to not form primary loops. Link " +
            std::to_string(linkIdx) + " is slip-link " +
            std::to_string(net.linkIsSliplink[linkIdxE]) + ".");
      }
    }
    net.springIndicesOfLinks[linkIdx].erase(
      last, net.springIndicesOfLinks[linkIdx].end());
  }
}

size_t
MEHPForceBalance::removePrimaryLoops(ForceBalanceNetwork& net,
                                     Eigen::VectorXd& displacements,
                                     Eigen::VectorXd& springPartitions) const
{
  size_t numRemoved = 0;
  for (long int springIdx = static_cast<long int>(net.nrOfSprings) - 1;
       springIdx >= 0;
       --springIdx) {
    const size_t springIdxU = static_cast<size_t>(springIdx);
    if (net.springIndexA[springIdx] == net.springIndexB[springIdx] &&
        net.localToGlobalSpringIndex[springIdxU].size() == 1 &&
        net.springsType[springIdx] != this->entanglementType) {
      this->removeSpringFollowingEntanglementLinks(
        net, displacements, springPartitions, springIdxU);
      springIdx = std::min<long int>(
        springIdx, static_cast<long int>(net.nrOfSprings) - 1);
      numRemoved += 1;
    }
  }

#ifndef NDEBUG
  this->validateNetwork();
#endif
  return numRemoved;
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
size_t
MEHPForceBalance::removeInactiveCrosslinks(ForceBalanceNetwork& net,
                                           Eigen::VectorXd& displacements,
                                           Eigen::VectorXd& springPartitions,
                                           const double tolerance) const
{
  size_t numRemoved = 0;
  //        this->removePrimaryLoops(net, displacements, springPartitions);
  // this->validateNetwork(net, displacements, springPartitions);
  // first, we remove all inactive springs
  for (long int springIdx = static_cast<long int>(net.nrOfSprings) - 1;
       springIdx >= 0;
       --springIdx) {
    if (springIdx >= static_cast<long int>(net.nrOfSprings)) {
      springIdx = static_cast<long int>(net.nrOfSprings) - 1;
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
    const size_t springIdxU = static_cast<size_t>(springIdx);
    std::vector<size_t> involvedPartialSprings =
      this->getAllPartialSpringIndicesAlong(net, springIdxU);
    assert(involvedPartialSprings.size() > 0);
    bool isActive = false;
    for (const size_t partialSpringIdx : involvedPartialSprings) {
      RUNTIME_EXP_IFN(
        net.coordinates.size() == displacements.size(),
        "Expected coordinates and displacements to have same size, got " +
          std::to_string(net.coordinates.size()) + " and " +
          std::to_string(displacements.size()) + ".");
      // assert(net.coordinates.size() == displacements.size());
      Eigen::Vector3d distance = this->evaluatePartialSpringDistance(
        net, displacements, partialSpringIdx);
      const Eigen::Index partialSpringIdxE =
        static_cast<Eigen::Index>(partialSpringIdx);
      const double contourLength =
        net.springsContourLength
          [net.partialToFullSpringIndex[partialSpringIdxE]];
      const double partition = springPartitions[partialSpringIdxE];
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
        net, displacements, springPartitions, springIdxU);

#ifndef NDEBUG
      assert(this->validateNetwork(net, displacements, springPartitions));
#endif
      numRemoved += 1;
    }
  }

  // then, we remove all crosslinkers that are 0- or 1-functional
  for (long int crosslinkIdx = static_cast<long int>(net.nrOfNodes) - 1;
       crosslinkIdx >= 0;
       --crosslinkIdx) {
    const size_t crosslinkIdxU = static_cast<size_t>(crosslinkIdx);
    assert(net.springIndicesOfLinks.size() > crosslinkIdxU);
    if (net.springIndicesOfLinks[crosslinkIdxU].size() == 0 // f = 0
    ) {
#ifdef DEBUG_REMOVAL
      std::cout << "Removing f = 0 x-link " << crosslinkIdx << std::endl;
#endif

      this->removeLink(net, displacements, crosslinkIdxU);
      numRemoved += 1;
#ifndef NDEBUG
      this->validateNetwork(net, displacements, springPartitions);
#endif
    }

    else if ( // or f = 1, NOT primary loop
      (net.springIndicesOfLinks[crosslinkIdxU].size() == 1) &&
      (XOR(net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdxU][0]]
                                   [0] == crosslinkIdxU,
           pylimer_tools::utils::last(
             net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdxU]
                                                              [0]]) ==
             crosslinkIdxU))) {
      assert(net.oldAtomTypes[crosslinkIdx] != this->entanglementType);
#ifdef DEBUG_REMOVAL
      std::cout << "Removing f = 1 x-link " << crosslinkIdx << std::endl;
#endif
      std::vector<size_t> entanglementLinksRemoved =
        this->getEntanglementLinkIndicesAlong(
          net, net.springIndicesOfLinks[crosslinkIdxU][0]);
      // need to first remove the spring
      const long int previousId = net.oldAtomIds[crosslinkIdx];
      this->removeSpringFollowingEntanglementLinks(
        net,
        displacements,
        springPartitions,
        net.springIndicesOfLinks[crosslinkIdxU][0]);
      // this also removed some intermediate entanglement links
      // -> crossLinkIdx is now at a different index
      size_t newCrosslinkIdx = crosslinkIdxU;
      for (const size_t linkIdx : entanglementLinksRemoved) {
        assert(linkIdx != crosslinkIdxU);
        if (linkIdx < crosslinkIdxU) {
          newCrosslinkIdx -= 1;
        }
      }
      const Eigen::Index newCrosslinkIdxE =
        static_cast<Eigen::Index>(newCrosslinkIdx);
      assert(net.oldAtomIds[newCrosslinkIdxE] == previousId);
      assert(net.oldAtomTypes[newCrosslinkIdxE] != this->entanglementType);
      assert(net.springIndicesOfLinks[newCrosslinkIdx].size() == 0);
      numRemoved += 1;
      // to then remove the crosslink
      this->removeLink(net, displacements, newCrosslinkIdx);
#ifdef DEBUG_REMOVAL
      std::cout << "Effectively removed f = 1 x-link " << newCrosslinkIdx
                << std::endl;
#endif
      crosslinkIdx = std::min<long int>(
        crosslinkIdx, static_cast<long int>(net.nrOfNodes) - 1);
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
size_t
MEHPForceBalance::breakTooLongSprings(ForceBalanceNetwork& net,
                                      Eigen::VectorXd& displacements,
                                      Eigen::VectorXd& springPartitions) const
{
  if (this->springBreakingLength <= 0.) {
    return 0;
  }

  size_t numBroken = 0;

  // iterate the springs, determine their distance, and determine if it
  // exceeds the breaking force
  for (long int partialSpringIdx =
         static_cast<long int>(net.nrOfPartialSprings);
       partialSpringIdx >= 0;
       --partialSpringIdx) {
    if (partialSpringIdx >= static_cast<long int>(net.nrOfPartialSprings)) {
      partialSpringIdx = static_cast<long int>(net.nrOfPartialSprings) - 1;
    }
    if (partialSpringIdx < 0) {
      break;
    }
    const size_t partialSpringIdxU = static_cast<size_t>(partialSpringIdx);
    const double len = this->getWeightedPartialSpringLength(
      net, displacements, springPartitions, partialSpringIdxU);
    if (len > this->springBreakingLength) {
      // break this spring
      numBroken += 1;
      this->breakPartialSpring(
        net, displacements, springPartitions, partialSpringIdxU);
    }
  }

  return numBroken;
}

/**
 * @brief Add slip-links to this system
 *
 * @param nrOfSliplinksToSample
 * @param cutoff
 * @param minimumNrOfSliplinks
 * @param excludeCrosslinks
 * @param seed
 */
size_t
MEHPForceBalance::randomlyAddSliplinks(const size_t nrOfSliplinksToSample,
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
                     "Expected a cutoff > 0.0, got " + std::to_string(cutoff) +
                       ".");
  if (this->universe.getNrOfAtoms() == 0 &&
      this->initialConfig.nrOfSprings > 0) {
    CrosslinkerUniverse crosslinkerUniverse;
    const size_t nNodes = this->initialConfig.nrOfNodes;
    crosslinkerUniverse.xlinkX.resize(nNodes);
    crosslinkerUniverse.xlinkY.resize(nNodes);
    crosslinkerUniverse.xlinkZ.resize(nNodes);
    crosslinkerUniverse.xlinkTypes.resize(nNodes, this->crossLinkerType);
    for (size_t i = 0; i < nNodes; ++i) {
      const Eigen::Index iE = static_cast<Eigen::Index>(i);
      crosslinkerUniverse.xlinkX[i] =
        this->initialConfig.coordinates[3 * iE + 0];
      crosslinkerUniverse.xlinkY[i] =
        this->initialConfig.coordinates[3 * iE + 1];
      crosslinkerUniverse.xlinkZ[i] =
        this->initialConfig.coordinates[3 * iE + 2];
      if (i < static_cast<size_t>(this->initialConfig.oldAtomTypes.size())) {
        crosslinkerUniverse.xlinkTypes[i] =
          this->initialConfig.oldAtomTypes[iE];
      }
    }
    const size_t nSprings = this->initialConfig.nrOfSprings;
    // Nodes of degree one that are not crosslinkers are chain ends
    // (e.g., dangling or free ends, whose initial coordinates are mere
    // placeholders); model them as free (unanchored) random-walk ends
    // rather than as fixed anchors.
    std::vector<size_t> nodeDegree(nNodes, 0);
    for (size_t s = 0; s < nSprings; ++s) {
      const Eigen::Index sE = static_cast<Eigen::Index>(s);
      nodeDegree[static_cast<size_t>(this->initialConfig.springIndexA[sE])] +=
        1;
      nodeDegree[static_cast<size_t>(this->initialConfig.springIndexB[sE])] +=
        1;
    }
    const auto isFreeEnd = [&](const size_t node) {
      return nodeDegree[node] == 1 &&
             crosslinkerUniverse.xlinkTypes[node] != this->crossLinkerType;
    };
    crosslinkerUniverse.strandFrom.resize(nSprings);
    crosslinkerUniverse.strandTo.resize(nSprings);
    crosslinkerUniverse.beadsInStrand.resize(nSprings);
    const double b =
      (this->defaultBondLength > 0.0) ? this->defaultBondLength : 1.0;
    const double b2 = b * b;
    crosslinkerUniverse.beadDistanceInStrand.assign(nSprings, b);
    crosslinkerUniverse.meanSquaredBeadDistanceInStrand.assign(nSprings, b2);
    // springs whose anchored end is B, such that the contour fraction
    // sampled from the anchor has to be mirrored w.r.t. the spring direction
    std::vector<bool> springIsReversed(nSprings, false);
    for (size_t s = 0; s < nSprings; ++s) {
      const Eigen::Index sE = static_cast<Eigen::Index>(s);
      const size_t nodeA =
        static_cast<size_t>(this->initialConfig.springIndexA[sE]);
      const size_t nodeB =
        static_cast<size_t>(this->initialConfig.springIndexB[sE]);
      const bool aIsFree = isFreeEnd(nodeA);
      const bool bIsFree = isFreeEnd(nodeB);
      // the number of bonds between the two nodes
      const int nBonds = static_cast<int>(
        std::round(this->initialConfig.springsContourLength[sE]));
      // CrosslinkerUniverse strand ends are not beads themselves, unless
      // they are free ends: then, they are the first/last bead of the strand
      const int nFreeEnds =
        static_cast<int>(aIsFree) + static_cast<int>(bIsFree);
      crosslinkerUniverse.beadsInStrand[s] =
        std::max(0, nBonds - 1 + nFreeEnds);
      if (aIsFree && !bIsFree) {
        crosslinkerUniverse.strandFrom[s] = static_cast<long int>(nodeB);
        crosslinkerUniverse.strandTo[s] = -1;
        springIsReversed[s] = true;
      } else {
        crosslinkerUniverse.strandFrom[s] =
          aIsFree ? -1 : static_cast<long int>(nodeA);
        crosslinkerUniverse.strandTo[s] =
          bIsFree ? -1 : static_cast<long int>(nodeB);
      }
    }

    pylimer_tools::topo::entanglement_detection::StrandPairEntanglements
      entanglements = pylimer_tools::topo::entanglement_detection::
        randomlyFindEntanglementsOnDemand(this->box,
                                          crosslinkerUniverse,
                                          nrOfSliplinksToSample,
                                          cutoff,
                                          0.0,
                                          minimumNrOfSliplinks,
                                          sameStrandCutoff,
                                          seed > 0 ? std::to_string(seed) : "",
                                          false);
    for (size_t e = 0; e < entanglements.size(); ++e) {
      if (springIsReversed[entanglements.strandIdx1[e]]) {
        entanglements.alpha1[e] = 1.0 - entanglements.alpha1[e];
      }
      if (springIsReversed[entanglements.strandIdx2[e]]) {
        entanglements.alpha2[e] = 1.0 - entanglements.alpha2[e];
      }
    }

    this->addSlipLinks(entanglements);
    return entanglements.size();
  }
  // RUNTIME_EXP_IFN(this->initialConfig.nrOfLinks ==
  //                   this->initialConfig.nrOfNodes,
  //                 "Slip-links are only added randomly when no other "
  //                 "slip-links are in place yet.");
  INVALIDARG_EXP_IFN(minimumNrOfSliplinks < this->universe.getNrOfAtoms() / 2,
                     "Minimum number of slip-links must be less than the "
                     "possible number of slip-links to place.");
  INVALIDARG_EXP_IFN(nrOfSliplinksToSample < this->universe.getNrOfAtoms() / 2,
                     "Number of slip-links to place must be less than "
                     "the possible number of slip-links to place.");
  // query all the crosslinker chains we actually use to place
  // slip-links on
  std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
    this->universe.getChainsWithCrosslinker(crossLinkerType);
  bool danglingChainsAreKept =
    this->initialConfig.nrOfNodes >
    this->universe.getAtomsOfType(crossLinkerType).size();
  std::unordered_map<size_t, size_t> atomToStrand;
  atomToStrand.reserve(this->universe.getNrOfAtoms());
  std::unordered_map<size_t, size_t> atomIdxInStrand;
  atomIdxInStrand.reserve(this->universe.getNrOfAtoms());
  for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
    const pylimer_tools::entities::Molecule& chain = crossLinkerChains[i];
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
      std::vector<pylimer_tools::entities::Atom> atoms =
        chain.getAtomsLinedUp(this->crossLinkerType);
      for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
        const pylimer_tools::entities::Atom& atom = atoms[atomIdx];
        if (atom.getType() != this->crossLinkerType) {
          const auto springIt = this->initialConfig.oldAtomIdToSpringIndex.find(
            static_cast<size_t>(atom.getId()));
          if (springIt != this->initialConfig.oldAtomIdToSpringIndex.end()) {
            atomToStrand.emplace(atom.getId(), springIt->second);
            atomIdxInStrand.emplace(atom.getId(), atomIdx);
          }
        }
      }
    }
  }

  const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
        this->universe,
        nrOfSliplinksToSample,
        cutoff,
        0.0,
        minimumNrOfSliplinks,
        sameStrandCutoff,
        seed > 0 ? std::to_string(seed) : "",
        this->crossLinkerType,
        excludeCrosslinks,
        pylimer_tools::topo::entanglement_detection::EntanglementFilteringMode::
          NONE);

  std::vector<double> slipLinkXs;
  slipLinkXs.reserve(entanglements.pairsOfAtoms.size());
  std::vector<double> slipLinkYs;
  slipLinkYs.reserve(entanglements.pairsOfAtoms.size());
  std::vector<double> slipLinkZs;
  slipLinkZs.reserve(entanglements.pairsOfAtoms.size());
  std::vector<size_t> slipLinkStrandA;
  slipLinkStrandA.reserve(entanglements.pairsOfAtoms.size());
  std::vector<size_t> slipLinkStrandB;
  slipLinkStrandB.reserve(entanglements.pairsOfAtoms.size());
  std::vector<double> slipLinkStrandAlpha;
  slipLinkStrandAlpha.reserve(entanglements.pairsOfAtoms.size());
  std::vector<double> slipLinkStrandBeta;
  slipLinkStrandBeta.reserve(entanglements.pairsOfAtoms.size());

  size_t nrOfSlipLinksPlaced = 0;
  for (const std::pair<size_t, size_t>& pair : entanglements.pairsOfAtoms) {
    const size_t atomId1 = pair.first;
    const size_t atomId2 = pair.second;
    const auto it1 = atomToStrand.find(atomId1);
    const auto it2 = atomToStrand.find(atomId2);
    if (it1 == atomToStrand.end() || it2 == atomToStrand.end()) {
      continue;
    }
    if (it1->second == it2->second) {
      continue;
    }
    const pylimer_tools::entities::Atom a1 =
      this->universe.getAtom(static_cast<int>(atomId1));
    const pylimer_tools::entities::Atom a2 =
      this->universe.getAtom(static_cast<int>(atomId2));
    const Eigen::Vector3d meanPositions = a1.meanPositionWith(a2, this->box);

    slipLinkXs.push_back(meanPositions[0]);
    slipLinkYs.push_back(meanPositions[1]);
    slipLinkZs.push_back(meanPositions[2]);
    slipLinkStrandA.push_back(it1->second);
    slipLinkStrandB.push_back(it2->second);
    slipLinkStrandAlpha.push_back(
      static_cast<double>(atomIdxInStrand.at(atomId1)) /
      static_cast<double>(
        this->initialConfig
          .springsContourLength[static_cast<Eigen::Index>(it1->second)]));
    slipLinkStrandBeta.push_back(
      static_cast<double>(atomIdxInStrand.at(atomId2)) /
      static_cast<double>(
        this->initialConfig
          .springsContourLength[static_cast<Eigen::Index>(it2->second)]));
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
size_t
MEHPForceBalance::addSliplinksBasedOnCycles(const int maxLoopLength)
{
  // std::cout << "Detecting slip-links based on cycles. Base memory
  // useage:
  // "
  //           << getCurrentRSS() << ", peak " << getPeakRSS() <<
  //           std::endl;
  std::vector<std::vector<igraph_integer_t>> loopEdges;
  std::vector<std::vector<igraph_integer_t>> loops = this->universe.findLoops(
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
  pylimer_tools::entities::Box universeBox = this->box;
  for (std::vector<igraph_integer_t> loop : loops) {
    // max & min coordinates
    // TODO: implement, such that the max & min work also with infinite
    // loops
    // assume each subsequent atom is bonded with a bond shorter
    // than half the bond
    Eigen::VectorXd alignedLoopCoordinates =
      this->universe.getUnwrappedVertexCoordinates(loop, universeBox);

    std::array<double, 3> minCoords;
    std::array<double, 3> maxCoords;

    for (int i = 0; i < 3; ++i) {
      const size_t iU = static_cast<size_t>(i);
      minCoords[iU] = alignedLoopCoordinates(
                        Eigen::seq(i, Eigen::placeholders::last - 2 + i, 3))
                        .minCoeff();
      maxCoords[iU] = alignedLoopCoordinates(
                        Eigen::seq(i, Eigen::placeholders::last - 2 + i, 3))
                        .maxCoeff();
      // want to see later that this direction is spanned fully
      // that's why we just push it so far that it is within the box
      if (maxCoords[iU] - minCoords[iU] > box.getL(i) ||
          (std::fabs(alignedLoopCoordinates[i] -
                     alignedLoopCoordinates
                       [3 * (static_cast<Eigen::Index>(loop.size()) - 1) + i]) >
           0.5 * box.getL(i))) {
        // special case in case the loop is periodic
        minCoords[iU] = box.getLowL(i);
        maxCoords[iU] = box.getHighL(i);
      } else {
        // adjust for box
        while (minCoords[static_cast<size_t>(i)] > box.getHighL(i)) {
          minCoords[static_cast<size_t>(i)] -= box.getL(i);
        }
        while (minCoords[static_cast<size_t>(i)] < box.getLowL(i)) {
          minCoords[static_cast<size_t>(i)] += box.getL(i);
        }
        while (maxCoords[static_cast<size_t>(i)] > box.getHighL(i)) {
          maxCoords[static_cast<size_t>(i)] -= box.getL(i);
        }
        while (maxCoords[static_cast<size_t>(i)] < box.getLowL(i)) {
          maxCoords[static_cast<size_t>(i)] += box.getL(i);
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
          this->initialConfig.oldAtomIdToSpringIndex[static_cast<size_t>(
            this->universe.getPropertyValue<long int>("id", loop[i]))]);
      }
    }
    std::vector<size_t> reducedLoopVec(reducedLoop.begin(), reducedLoop.end());
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
  auto hash_integer_pair = [](const long int a, const long int b) -> long int {
    return a >= b ? a * a + a + b : a + b * b; // where a, b >= 0
  };
  std::cout << "Searching for intersections..." << std::endl;

  for (size_t i = 0; i < loops.size(); ++i) {
    // reserve enough space
    size_t estimateNrOfSliplinks =
      static_cast<size_t>(static_cast<double>(loops.size() * loops.size()) *
                          this->initialConfig.meanSpringContourLength);
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

      for (size_t dir = 0; dir < 3; ++dir) {
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
            const size_t involvedIntersection =
              static_cast<size_t>(intersectionsOfEdges.at(edgePairHash));
            std::get<1>(relevantIntersections[involvedIntersection]).insert(i);
            std::get<1>(relevantIntersections[involvedIntersection]).insert(j);
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
    loops[i] = std::vector<igraph_integer_t>();
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

  for (intersection_loops_tuple intersectionAndLoops : relevantIntersections) {
    pylimer_tools::entities::LoopIntersectionInfo intersection =
      std::get<0>(intersectionAndLoops);
    // TODO: this is yet the most naïve way to add these.
    // ideally, we would also check the back-and-forth, etc.
    slipLinkXs.push_back(intersection.intersectionPoint[0]);
    slipLinkYs.push_back(intersection.intersectionPoint[1]);
    slipLinkZs.push_back(intersection.intersectionPoint[2]);
    // TODO: decide on the atoms to use as a reference
    slipLinkStrandA.push_back(
      this->initialConfig.oldAtomIdToSpringIndex[static_cast<size_t>(
        intersection.involvedAtoms[0].getId())]);
    slipLinkStrandB.push_back(
      this->initialConfig.oldAtomIdToSpringIndex[static_cast<size_t>(
        intersection.involvedAtoms[3].getId())]);
    slipLinkStrandAlpha.push_back(
      static_cast<double>(atomIdxInStrand[static_cast<size_t>(
        intersection.involvedAtoms[0].getId())]) /
      this->initialConfig.meanSpringContourLength);
    slipLinkStrandBeta.push_back(
      static_cast<double>(atomIdxInStrand[static_cast<size_t>(
        intersection.involvedAtoms[3].getId())]) /
      this->initialConfig.meanSpringContourLength);
    std::vector<size_t> localLoops(std::get<1>(intersectionAndLoops).begin(),
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
void
MEHPForceBalance::removeSpring(ForceBalanceNetwork& net,
                               Eigen::VectorXd& displacements,
                               Eigen::VectorXd& springPartitions,
                               const size_t springIdx) const
{
#ifdef DEBUG_REMOVAL
  std::cout << "Starting to remove spring " << springIdx << std::endl;
#endif
  INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                     "Can only remove springs, not partial springs.");
#ifndef NDEBUG
  Eigen::VectorXd allTotalSpringDistancesBefore =
    this->evaluateSpringLengths(net, displacements, this->is2D);
#endif

  const std::vector<size_t> affectedLinks = net.linkIndicesOfSprings[springIdx];
  std::vector<size_t> uniqueAffectedLinks = net.linkIndicesOfSprings[springIdx];
  std::ranges::sort(uniqueAffectedLinks);
  uniqueAffectedLinks.erase(
    std::unique(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end()),
    uniqueAffectedLinks.end());

  // remove the link to the link, höhö
  for (const size_t affectedLinkIdx : uniqueAffectedLinks) {
    const Eigen::Index affectedLinkIdxE =
      static_cast<Eigen::Index>(affectedLinkIdx);
    RUNTIME_EXP_IFN(
      std::find(net.springIndicesOfLinks[affectedLinkIdx].begin(),
                net.springIndicesOfLinks[affectedLinkIdx].end(),
                springIdx) != net.springIndicesOfLinks[affectedLinkIdx].end(),
      "Link must have a connection to the spring, too. Did not find "
      "spring " +
        std::to_string(springIdx) + " in link " +
        std::to_string(affectedLinkIdx) + ", got " +
        pylimer_tools::utils::join(
          net.springIndicesOfLinks[affectedLinkIdx].begin(),
          net.springIndicesOfLinks[affectedLinkIdx].end(),
          std::string(", ")) +
        ".");
    if (net.linkIsSliplink[affectedLinkIdxE]) {
      RUNTIME_EXP_IFN(
        net.springIndicesOfLinks[affectedLinkIdx].size() <= 2,
        "Expect slip-link to be associated with 2 springs only, got " +
          pylimer_tools::utils::join(
            net.springIndicesOfLinks[affectedLinkIdx].begin(),
            net.springIndicesOfLinks[affectedLinkIdx].end(),
            std::string(", ")) +
          ".");
    }

    const size_t found =
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
    if (net.linkIsSliplink[affectedLinkIdxE]) {
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
  const std::ptrdiff_t springIdxDiff = static_cast<std::ptrdiff_t>(springIdx);
  net.localToGlobalSpringIndex.erase(net.localToGlobalSpringIndex.begin() +
                                     springIdxDiff);
  net.springToMoleculeIds.erase(net.springToMoleculeIds.begin() +
                                springIdxDiff);
  net.linkIndicesOfSprings.erase(net.linkIndicesOfSprings.begin() +
                                 springIdxDiff);
  const unsigned int springIdxU = static_cast<unsigned int>(springIdx);
  pylimer_tools::utils::removeRow(net.springsContourLength, springIdxU);
  pylimer_tools::utils::removeRow(net.springsType, springIdxU);
  pylimer_tools::utils::removeRow(net.springIndexA, springIdxU);
  pylimer_tools::utils::removeRows(
    net.springCoordinateIndexA, springIdxU * 3, 3);
  pylimer_tools::utils::removeRow(net.springIndexB, springIdxU);
  pylimer_tools::utils::removeRows(
    net.springCoordinateIndexB, springIdxU * 3, 3);
  pylimer_tools::utils::removeRow(net.springIsActive, springIdxU);

  // need to remove descending
  std::ranges::sort(affectedPartialSprings, std::greater<size_t>());
  for (const size_t partialSpringIdx : affectedPartialSprings) {
    const unsigned int partialSpringIdxU =
      static_cast<unsigned int>(partialSpringIdx);
    pylimer_tools::utils::removeRow(net.springPartIndexA, partialSpringIdxU);
    pylimer_tools::utils::removeRows(
      net.springPartCoordinateIndexA, 3 * partialSpringIdxU, 3);
    pylimer_tools::utils::removeRow(net.springPartIndexB, partialSpringIdxU);
    pylimer_tools::utils::removeRows(
      net.springPartCoordinateIndexB, 3 * partialSpringIdxU, 3);
    pylimer_tools::utils::removeRows(
      net.springPartBoxOffset, 3 * partialSpringIdxU, 3);
    pylimer_tools::utils::removeRow(net.partialToFullSpringIndex,
                                    partialSpringIdxU);
    pylimer_tools::utils::removeRow(net.partialSpringIsPartial,
                                    partialSpringIdxU);
    pylimer_tools::utils::removeRow(springPartitions, partialSpringIdxU);
  }
  assert(springPartitions.size() ==
         static_cast<Eigen::Index>(net.nrOfPartialSprings));

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
        net.loops[loopIdx].erase(net.loops[loopIdx].begin() +
                                 static_cast<std::ptrdiff_t>(i));
      }
      if (net.loops[loopIdx][i] > springIdx) {
        net.loops[loopIdx][i] -= 1;
      }
    }
  }

  // then, update the partial springs
  // decrease by one if larger then springIdx
  const int springIdxI = static_cast<int>(springIdx);
  assert((net.partialToFullSpringIndex != springIdxI).all());
  net.partialToFullSpringIndex -=
    (net.partialToFullSpringIndex > springIdxI).cast<int>();

  assert(net.localToGlobalSpringIndex.size() == net.nrOfSprings);
  for (size_t loopingSpringIdx = 0; loopingSpringIdx < net.nrOfSprings;
       ++loopingSpringIdx) {
    for (size_t i = 0;
         i < net.localToGlobalSpringIndex[loopingSpringIdx].size();
         ++i) {
      for (const size_t partialSpringIdx : affectedPartialSprings) {
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
  linksToRemove.reserve(affectedLinks.size() -
                        2); // keep the first and last links
  for (size_t i = 1; i < affectedLinks.size() - 1; ++i) {
    linksToRemove.push_back(affectedLinks[i]);
  }

  // need to remove descending to remove need to renumber these as well
  std::ranges::sort(linksToRemove, std::greater<size_t>());
  linksToRemove.erase(std::ranges::unique(linksToRemove).begin(),
                      linksToRemove.end());
  if (linksToRemove.size() > 2) {
    assert(linksToRemove[0] > linksToRemove[1]);
  }

  Eigen::ArrayXb springIsAffected =
    Eigen::ArrayXb::Constant(static_cast<Eigen::Index>(net.nrOfSprings), false);
  for (size_t outermostI = 0; outermostI < linksToRemove.size(); ++outermostI) {
    size_t slipLinkIdx = linksToRemove[outermostI];
    const int slipLinkIdxI = static_cast<int>(slipLinkIdx);
    assert(net.linkIsSliplink[static_cast<Eigen::Index>(slipLinkIdx)]);
    // first, merge the two other partial springs
    std::vector<size_t> springsOfLink = net.springIndicesOfLinks[slipLinkIdx];
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
    for (long int springInLinkIdx =
           static_cast<long int>(springsOfLink.size()) - 1;
         springInLinkIdx >= 0;
         --springInLinkIdx) {
      for (size_t partialSpringIdx :
           net.localToGlobalSpringIndex[springsOfLink[static_cast<size_t>(
             springInLinkIdx)]]) {
        const Eigen::Index partialSpringIdxE =
          static_cast<Eigen::Index>(partialSpringIdx);
        if (net.springPartIndexA[partialSpringIdxE] == slipLinkIdxI ||
            net.springPartIndexB[partialSpringIdxE] == slipLinkIdxI) {
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
      const size_t partialSpringToKeep =
        std::min(involvedPartialSprings[0], involvedPartialSprings[1]);
      const size_t partialSpringToRemove =
        std::max(involvedPartialSprings[0], involvedPartialSprings[1]);
      const Eigen::Index partialSpringToKeepE =
        static_cast<Eigen::Index>(partialSpringToKeep);
      const Eigen::Index partialSpringToRemoveE =
        static_cast<Eigen::Index>(partialSpringToRemove);
      assert(partialSpringToKeep != partialSpringToRemove);
      assert(net.partialToFullSpringIndex[partialSpringToKeepE] ==
             net.partialToFullSpringIndex[partialSpringToRemoveE]);
      // actually do the merge
      this->mergePartialSprings(net,
                                displacements,
                                springPartitions,
                                partialSpringToRemove,
                                partialSpringToKeep,
                                slipLinkIdx);
      // total distance changes -> cannot use for checking the total
      // distance
      springIsAffected[net.partialToFullSpringIndex[partialSpringToKeepE]] =
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
#ifndef NDEBUG
  Eigen::VectorXd allTotalSpringDistancesAfter =
    this->evaluateSpringLengths(net, displacements, this->is2D);
  assert(allTotalSpringDistancesAfter.size() ==
         static_cast<Eigen::Index>(net.nrOfSprings));
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(net.nrOfSprings);
       ++i) {
    const Eigen::Index correspondingOldIdx =
      i >= static_cast<Eigen::Index>(springIdx) ? i + 1 : i;
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
#endif
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
void
MEHPForceBalance::removeSpringFollowingEntanglementLinks(
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
  std::ranges::sort(springsToRemove, std::greater<>());
  for (const size_t springIdxToDelete : springsToRemove) {
    assert(net.springsType[static_cast<Eigen::Index>(springIdxToDelete)] !=
           this->entanglementType);
    this->removeSpring(net, displacements, springPartitions, springIdxToDelete);
  }
  std::sort(linksToRemove.begin(), linksToRemove.end(), std::greater<>());
  for (size_t linkIdxToDelete : linksToRemove) {
    assert(net.springIndicesOfLinks[linkIdxToDelete].size() <= 1);
    assert(net.oldAtomTypes[static_cast<Eigen::Index>(linkIdxToDelete)] ==
           this->entanglementType);
    if (net.springIndicesOfLinks[linkIdxToDelete].size() == 1) {
      assert(net.springsType[static_cast<Eigen::Index>(
               net.springIndicesOfLinks[linkIdxToDelete][0])] ==
             this->entanglementType);
#ifdef DEBUG_REMOVAL
      std::cout << "Removing additional spring between entanglement links "
                << net.springIndicesOfLinks[linkIdxToDelete][0] << std::endl;
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
void
MEHPForceBalance::breakPartialSpring(ForceBalanceNetwork& net,
                                     Eigen::VectorXd& displacements,
                                     Eigen::VectorXd& springPartitions,
                                     const size_t partialSpringIdx) const
{
  this->removeSpringFollowingEntanglementLinks(
    net,
    displacements,
    springPartitions,
    static_cast<size_t>(net.partialToFullSpringIndex[static_cast<Eigen::Index>(
      partialSpringIdx)]));
};

/**
 * @brief remove a link from the network
 *
 * @param net
 * @param displacements
 * @param linkIdx
 */
void
MEHPForceBalance::removeLink(ForceBalanceNetwork& net,
                             Eigen::VectorXd& displacements,
                             const size_t linkIdx)
{
  INVALIDARG_EXP_IFN(net.springIndicesOfLinks[linkIdx].size() == 0,
                     "The springs have to be removed or re-linked before "
                     "removing the link.");
#ifdef DEBUG_REMOVAL
  std::cout << "Removing link " << linkIdx << std::endl;
#endif

  const unsigned int linkIdxU = static_cast<unsigned int>(linkIdx);
  pylimer_tools::utils::removeRows(net.coordinates, linkIdxU * 3, 3);
  pylimer_tools::utils::removeRows(displacements, linkIdxU * 3, 3);

  const Eigen::Index linkIdxE = static_cast<Eigen::Index>(linkIdx);
  const std::ptrdiff_t linkIdxD = static_cast<std::ptrdiff_t>(linkIdx);
  if (!net.linkIsSliplink[linkIdxE]) {
    net.nrOfNodes -= 1;
    net.oldAtomIdToSpringIndex.erase(
      static_cast<size_t>(net.oldAtomIds[linkIdxE]));
    pylimer_tools::utils::removeRow(net.oldAtomIds, linkIdxU);
    pylimer_tools::utils::removeRow(net.oldAtomTypes, linkIdxU);
  } else {
    pylimer_tools::utils::removeRow(
      net.nrOfCrosslinkSwapsEndured,
      static_cast<unsigned int>(linkIdx - net.nrOfNodes));
    if (net.loopsOfSliplink.size() > 0) {
      net.loopsOfSliplink.erase(net.loopsOfSliplink.begin() + linkIdxD);
    }
  }
  net.nrOfLinks -= 1;
  pylimer_tools::utils::removeRow(net.linkIsSliplink, linkIdxU);
  net.springIndicesOfLinks.erase(net.springIndicesOfLinks.begin() + linkIdxD);

  // renumber the remaining links
  for (size_t i = 0; i < net.linkIndicesOfSprings.size(); ++i) {
    for (size_t j = 0; j < net.linkIndicesOfSprings[i].size(); ++j) {
#ifndef NDEBUG
      RUNTIME_EXP_IFN(
        net.linkIndicesOfSprings[i][j] != linkIdx,
        "Expected not to find link to remove " + std::to_string(linkIdx) +
          " in any spring, found in spring " + std::to_string(i) + ", " +
          pylimer_tools::utils::join(net.linkIndicesOfSprings[i].begin(),
                                     net.linkIndicesOfSprings[i].end(),
                                     std::string(", ")) +
          ".");
#endif
      if (net.linkIndicesOfSprings[i][j] > linkIdx) {
        net.linkIndicesOfSprings[i][j] -= 1;
      }
    }
  }
  //
  const int linkIdxI = static_cast<int>(linkIdx);
  assert(net.springPartIndexA.size() == net.springPartIndexB.size());
  for (Eigen::Index i = 0; i < net.springPartIndexA.size(); ++i) {
#ifndef NDEBUG
    RUNTIME_EXP_IFN(
      net.springPartIndexA[i] != linkIdxI,
      "Exected link " + std::to_string(linkIdx) +
        " to not be linked anywhere anymore, found in partial spring " +
        std::to_string(i) + ".");
#endif
    if (net.springPartIndexA[i] > linkIdxI) {
      net.springPartIndexA[i] -= 1;
      net.springPartCoordinateIndexA[3 * i] -= 3;
      net.springPartCoordinateIndexA[3 * i + 1] -= 3;
      net.springPartCoordinateIndexA[3 * i + 2] -= 3;
    }
#ifndef NDEBUG
    RUNTIME_EXP_IFN(
      net.springPartIndexB[i] != linkIdxI,
      "Exected link " + std::to_string(linkIdx) +
        " to not be linked anywhere anymore, found in partial spring " +
        std::to_string(i) + ".");
#endif
    if (net.springPartIndexB[i] > linkIdxI) {
      net.springPartIndexB[i] -= 1;
      net.springPartCoordinateIndexB[3 * i] -= 3;
      net.springPartCoordinateIndexB[3 * i + 1] -= 3;
      net.springPartCoordinateIndexB[3 * i + 2] -= 3;
    }
  }
  // some validation of current state
  assert(net.springIndexA.size() == net.springIndexB.size());
  assert((net.springIndexA != linkIdxI).all());
  assert((net.springIndexB != linkIdxI).all());
  // renumber springs
  net.springIndexA -= (net.springIndexA > linkIdxI).cast<int>();
  net.springCoordinateIndexA -=
    3 * (net.springCoordinateIndexA > 3 * linkIdxI).cast<int>();
  net.springIndexB -= (net.springIndexB > linkIdxI).cast<int>();
  net.springCoordinateIndexB -=
    3 * (net.springCoordinateIndexB > 3 * linkIdxI).cast<int>();
}

/**
 * @brief Combine two springs
 *
 * @param net
 * @param springPartitions
 */
void
MEHPForceBalance::mergePartialSprings(ForceBalanceNetwork& net,
                                      const Eigen::VectorXd& u,
                                      Eigen::VectorXd& springPartitions,
                                      const size_t removedPartialSpringIdx,
                                      const size_t keptPartialSpringIdx,
                                      const size_t linkToReduce,
                                      const bool skipEigenResize) const
{
  const Eigen::Index keptPartialSpringIdxE =
    static_cast<Eigen::Index>(keptPartialSpringIdx);
  const Eigen::Index removedPartialSpringIdxE =
    static_cast<Eigen::Index>(removedPartialSpringIdx);
  const int linkToReduceI = static_cast<int>(linkToReduce);
  INVALIDARG_EXP_IFN(
    net.linkIsSliplink[static_cast<Eigen::Index>(linkToReduce)],
    "The link to reduce must be a slip-link");
  INVALIDARG_EXP_IFN(keptPartialSpringIdx != removedPartialSpringIdx,
                     "Cannot merge one spring with the same one.");
  INVALIDARG_EXP_IFN(
    net.partialToFullSpringIndex[keptPartialSpringIdxE] ==
      net.partialToFullSpringIndex[removedPartialSpringIdxE],
    "The partial springs must be part of the same spring to merge them.");
  INVALIDARG_EXP_IFN(
    (net.springPartIndexA[keptPartialSpringIdxE] == linkToReduceI &&
     net.springPartIndexB[removedPartialSpringIdxE] == linkToReduceI) ||
      (net.springPartIndexB[keptPartialSpringIdxE] == linkToReduceI &&
       net.springPartIndexA[removedPartialSpringIdxE] == linkToReduceI),
    "Link to reduce must be part of the partial springs "
    "that are to be removed and kept, got " +
      std::to_string(net.springPartIndexA[keptPartialSpringIdxE]) + " and " +
      std::to_string(net.springPartIndexB[keptPartialSpringIdxE]) +
      " in kept spring, and got " +
      std::to_string(net.springPartIndexA[removedPartialSpringIdxE]) + " and " +
      std::to_string(net.springPartIndexB[removedPartialSpringIdxE]) +
      "in removed spring, instead of " + std::to_string(linkToReduce) + ".");
#ifdef DEBUG_REMOVAL
  std::cout << "Merging partial springs " << removedPartialSpringIdx << " and "
            << keptPartialSpringIdx << " around " << linkToReduce << std::endl;
#endif

  const Eigen::Vector3d distanceBefore =
    this->evaluatePartialSpringDistance(
      net, u, removedPartialSpringIdx, this->is2D, false) +
    this->evaluatePartialSpringDistance(
      net, u, keptPartialSpringIdx, this->is2D, false);
  const size_t fullSpringIdx =
    static_cast<size_t>(net.partialToFullSpringIndex[keptPartialSpringIdxE]);
  // start with removal
  net.nrOfPartialSprings -= 1;
  // tell the kept one their new end
  // NOTE: if is possible, if the removedPartialSpring is a primary loop,
  // that this procedure is ambiguous.
  const bool removedIsA =
    net.springPartIndexA[removedPartialSpringIdxE] == linkToReduceI;
  const size_t newEnd = static_cast<size_t>(
    removedIsA ? net.springPartIndexB[removedPartialSpringIdxE]
               : net.springPartIndexA[removedPartialSpringIdxE]);
  if (removedIsA) {
    assert(net.springPartIndexB[keptPartialSpringIdxE] == linkToReduceI);
    net.springPartIndexB[keptPartialSpringIdxE] = static_cast<int>(newEnd);
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      net.springPartCoordinateIndexB[3 * keptPartialSpringIdxE + dir] =
        static_cast<int>(3 * newEnd) + static_cast<int>(dir);
    }
  } else {
    assert(net.springPartIndexA[keptPartialSpringIdxE] == linkToReduceI);
    net.springPartIndexA[keptPartialSpringIdxE] = static_cast<int>(newEnd);
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      net.springPartCoordinateIndexA[3 * keptPartialSpringIdxE + dir] =
        static_cast<int>(3 * newEnd) + static_cast<int>(dir);
    }
  }
  net.springPartBoxOffset.segment(3 * keptPartialSpringIdxE, 3) +=
    net.springPartBoxOffset.segment(3 * removedPartialSpringIdxE, 3);
  // remove the spring from the link
  // NOTE: currently, we allow it not to be present,
  // as it might be removed earlier already
  // It is anyway the case, that this function does not necessarily
  // keep the network valid
  size_t found =
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
  springPartitions[keptPartialSpringIdxE] +=
    springPartitions[removedPartialSpringIdxE];
  for (long int j =
         static_cast<long int>(net.linkIndicesOfSprings[fullSpringIdx].size()) -
         1;
       j >= 0;
       --j) {
    const size_t jU = static_cast<size_t>(j);
    if (net.linkIndicesOfSprings[fullSpringIdx][jU] == linkToReduce) {
      if (removed == 0) {
        if ((j > 0 && net.localToGlobalSpringIndex[fullSpringIdx][jU - 1] ==
                        removedPartialSpringIdx) ||
            (jU < net.localToGlobalSpringIndex[fullSpringIdx].size() &&
             net.localToGlobalSpringIndex[fullSpringIdx][jU] ==
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
  for (long int j = static_cast<long int>(
                      net.localToGlobalSpringIndex[fullSpringIdx].size()) -
                    1;
       j >= 0;
       --j) {
    if (net.localToGlobalSpringIndex[fullSpringIdx][static_cast<size_t>(j)] ==
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
  net.partialSpringIsPartial[keptPartialSpringIdxE] =
    net.linkIndicesOfSprings[fullSpringIdx].size() > 2;
  // actually remove the rows
  const unsigned int removedPartialSpringIdxU =
    static_cast<unsigned int>(removedPartialSpringIdx);
  pylimer_tools::utils::removeRow(
    net.partialSpringIsPartial, removedPartialSpringIdxU, skipEigenResize);
  pylimer_tools::utils::removeRow(
    net.partialToFullSpringIndex, removedPartialSpringIdxU, skipEigenResize);
  pylimer_tools::utils::removeRow(
    springPartitions, removedPartialSpringIdxU, skipEigenResize);
  pylimer_tools::utils::removeRow(
    net.springPartIndexA, removedPartialSpringIdxU, skipEigenResize);
  pylimer_tools::utils::removeRow(
    net.springPartIndexB, removedPartialSpringIdxU, skipEigenResize);
  pylimer_tools::utils::removeRows(net.springPartCoordinateIndexA,
                                   3 * removedPartialSpringIdxU,
                                   3,
                                   skipEigenResize);
  pylimer_tools::utils::removeRows(net.springPartCoordinateIndexB,
                                   3 * removedPartialSpringIdxU,
                                   3,
                                   skipEigenResize);
  pylimer_tools::utils::removeRows(
    net.springPartBoxOffset, 3 * removedPartialSpringIdxU, 3, skipEigenResize);
  // renumber stuff
  for (size_t loopSpringIdx = 0;
       loopSpringIdx < net.localToGlobalSpringIndex.size();
       ++loopSpringIdx) {
    for (size_t i = 0; i < net.localToGlobalSpringIndex[loopSpringIdx].size();
         ++i) {
      if (net.localToGlobalSpringIndex[loopSpringIdx][i] >
          removedPartialSpringIdx) {
        net.localToGlobalSpringIndex[loopSpringIdx][i] -= 1;
      }
    }
  }

  // validation
  const size_t newSpringIdx =
    keptPartialSpringIdx -
    (keptPartialSpringIdx > removedPartialSpringIdx ? 1 : 0);
  const Eigen::Vector3d newDistance = this->evaluatePartialSpringDistance(
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
 * @param displacements
 * @param springPartitions
 * @param removedSpringIdx
 */
void
MEHPForceBalance::mergeSprings(ForceBalanceNetwork& net,
                               const Eigen::VectorXd& displacements,
                               Eigen::VectorXd& springPartitions,
                               const size_t removedSpringIdx,
                               const size_t keptSpringIdx,
                               const size_t linkToReduce) const
{
  const Eigen::Index linkToReduceE = static_cast<Eigen::Index>(linkToReduce);
  const Eigen::Index keptSpringIdxE = static_cast<Eigen::Index>(keptSpringIdx);
  const Eigen::Index removedSpringIdxE =
    static_cast<Eigen::Index>(removedSpringIdx);
  const int linkToReduceI = static_cast<int>(linkToReduce);
  INVALIDARG_EXP_IFN(removedSpringIdx < net.nrOfSprings &&
                       keptSpringIdx < net.nrOfSprings,
                     "Only full springs can be merged.");
  INVALIDARG_EXP_IFN(!net.linkIsSliplink[linkToReduceE],
                     "The link to reduce must be a crosslink");
  INVALIDARG_EXP_IFN(keptSpringIdx != removedSpringIdx,
                     "Cannot replace one spring with the same one.");
  INVALIDARG_EXP_IFN(net.springsType[keptSpringIdxE] != this->entanglementType,
                     "Should not merge entanglement springs.");
  INVALIDARG_EXP_IFN(net.springsType[removedSpringIdxE] !=
                       this->entanglementType,
                     "Should not merge entanglement springs.");
#ifndef NDEBUG
  if (net.oldAtomTypes[linkToReduceE] != this->entanglementType) {
    std::vector<size_t> fullSprings1 =
      this->getAllFullSpringIndicesAlong(net, removedSpringIdx);
    std::vector<size_t> fullSprings2 =
      this->getAllFullSpringIndicesAlong(net, keptSpringIdx);
    std::ranges::sort(fullSprings1);
    std::ranges::sort(fullSprings2);
    INVALIDARG_EXP_IFN(
      !(pylimer_tools::utils::equal(fullSprings1, fullSprings2)),
      "Cannot merge such that a primary loop made of "
      "entanglements results.");
  }
#endif

  // handle links
  const std::vector<size_t> removedSpringsLinks =
    net.linkIndicesOfSprings[removedSpringIdx];
  const std::vector<size_t> keptSpringsLinks =
    net.linkIndicesOfSprings[keptSpringIdx];

  const size_t removedPartialSpringIdx =
    (removedSpringsLinks[removedSpringsLinks.size() - 1] == linkToReduce)
      ? pylimer_tools::utils::last(
          net.localToGlobalSpringIndex[removedSpringIdx])
      : net.localToGlobalSpringIndex[removedSpringIdx][0];
  const size_t remainingPartialSpringIdx =
    (keptSpringsLinks[keptSpringsLinks.size() - 1] == linkToReduce)
      ? pylimer_tools::utils::last(net.localToGlobalSpringIndex[keptSpringIdx])
      : net.localToGlobalSpringIndex[keptSpringIdx][0];
  const Eigen::Index removedPartialSpringIdxE =
    static_cast<Eigen::Index>(removedPartialSpringIdx);
  const Eigen::Index remainingPartialSpringIdxE =
    static_cast<Eigen::Index>(remainingPartialSpringIdx);

  RUNTIME_EXP_IFN(
    net.springPartIndexA[removedPartialSpringIdxE] == linkToReduceI ||
      net.springPartIndexB[removedPartialSpringIdxE] == linkToReduceI,
    "Did not detect correct partial springs");
  RUNTIME_EXP_IFN(
    net.springPartIndexA[remainingPartialSpringIdxE] == linkToReduceI ||
      net.springPartIndexB[remainingPartialSpringIdxE] == linkToReduceI,
    "Did not detect correct partial springs");

  Eigen::Vector3d distanceBefore = this->evaluatePartialSpringDistance(
    net, displacements, removedPartialSpringIdx, this->is2D, false);
  const Eigen::Vector3d distanceBeforeRemainingSpring =
    this->evaluatePartialSpringDistance(
      net, displacements, remainingPartialSpringIdx, this->is2D, false);

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
  for (const size_t partialSpringIndex :
       net.localToGlobalSpringIndex[removedSpringIdx]) {
    net
      .partialToFullSpringIndex[static_cast<Eigen::Index>(partialSpringIndex)] =
      static_cast<int>(keptSpringIdx);
  }
  // std::cout << "Kept spring is "
  //           << pylimer_tools::utils::join(keptSpringsLinks.begin(),
  //                                         keptSpringsLinks.end(),
  //                                         std::string(", "))
  //           << std::endl;
  // actually merge the springs
  if (keptSpringsLinks[keptSpringsLinks.size() - 1] == linkToReduce) {
    // add to end...
    if (removedSpringsLinks[removedSpringsLinks.size() - 1] == linkToReduce) {
      // std::cout << "End end" << std::endl;
      // ...from end
      net.linkIndicesOfSprings[keptSpringIdx][keptSpringsLinks.size() - 1] =
        removedSpringsLinks[removedSpringsLinks.size() - 2];
      for (size_t i = 3; i <= removedSpringsLinks.size(); ++i) {
        net.linkIndicesOfSprings[keptSpringIdx].push_back(
          removedSpringsLinks[removedSpringsLinks.size() - i]);
      }
      for (long int i =
             static_cast<long int>(
               net.localToGlobalSpringIndex[removedSpringIdx].size()) -
             2;
           i >= 0;
           --i) {
        net.localToGlobalSpringIndex[keptSpringIdx].push_back(
          net.localToGlobalSpringIndex[removedSpringIdx]
                                      [static_cast<size_t>(i)]);
      }
      // invert the direction of these transferred partial springs
      for (const size_t partialSpringIdxToInvert :
           net.localToGlobalSpringIndex[removedSpringIdx]) {
        if (partialSpringIdxToInvert == removedPartialSpringIdx) {
          continue;
        }
        const Eigen::Index partialSpringIdxToInvertE =
          static_cast<Eigen::Index>(partialSpringIdxToInvert);
        std::swap(net.springPartIndexA[partialSpringIdxToInvertE],
                  net.springPartIndexB[partialSpringIdxToInvertE]);
        for (Eigen::Index dir = 0; dir < 3; ++dir) {
          std::swap(
            net.springPartCoordinateIndexA[3 * partialSpringIdxToInvertE + dir],
            net
              .springPartCoordinateIndexB[3 * partialSpringIdxToInvertE + dir]);
        }
        net.springPartBoxOffset.segment(3 * partialSpringIdxToInvertE, 3) *=
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
    RUNTIME_EXP_IFN(keptSpringsLinks[0] == linkToReduce, "How could this be?");
    // add to start...
    if (removedSpringsLinks[removedSpringsLinks.size() - 1] == linkToReduce) {
      // std::cout << "Start end" << std::endl;
      // from end
      net.linkIndicesOfSprings[keptSpringIdx][0] =
        removedSpringsLinks[removedSpringsLinks.size() - 2];
      for (size_t i = 3; i <= removedSpringsLinks.size(); ++i) {
        net.linkIndicesOfSprings[keptSpringIdx].insert(
          net.linkIndicesOfSprings[keptSpringIdx].begin(),
          removedSpringsLinks[removedSpringsLinks.size() - i]);
      }
      for (long int i =
             static_cast<long int>(
               net.localToGlobalSpringIndex[removedSpringIdx].size()) -
             2;
           i >= 0;
           --i) {
        net.localToGlobalSpringIndex[keptSpringIdx].insert(
          net.localToGlobalSpringIndex[keptSpringIdx].begin(),
          net.localToGlobalSpringIndex[removedSpringIdx]
                                      [static_cast<size_t>(i)]);
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
      for (const size_t partialSpringIdxToInvert :
           net.localToGlobalSpringIndex[removedSpringIdx]) {
        if (partialSpringIdxToInvert == removedPartialSpringIdx) {
          continue;
        }
        const Eigen::Index partialSpringIdxToInvertE =
          static_cast<Eigen::Index>(partialSpringIdxToInvert);
        std::swap(net.springPartIndexA[partialSpringIdxToInvertE],
                  net.springPartIndexB[partialSpringIdxToInvertE]);
        for (Eigen::Index dir = 0; dir < 3; ++dir) {
          std::swap(
            net.springPartCoordinateIndexA[3 * partialSpringIdxToInvertE + dir],
            net
              .springPartCoordinateIndexB[3 * partialSpringIdxToInvertE + dir]);
        }
        net.springPartBoxOffset.segment(3 * partialSpringIdxToInvertE, 3) *=
          -1.;
      }
      distanceBefore -= distanceBeforeRemainingSpring;
      distanceBefore *= -1.;
    }
  }
  RUNTIME_EXP_IFN(
    std::find(net.linkIndicesOfSprings[keptSpringIdx].begin(),
              net.linkIndicesOfSprings[keptSpringIdx].end(),
              linkToReduce) == net.linkIndicesOfSprings[keptSpringIdx].end(),
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
  for (const size_t linkOfRemovedSpring : removedSpringsLinks) {
    for (size_t i = 0; i < net.springIndicesOfLinks[linkOfRemovedSpring].size();
         ++i) {
      if (net.springIndicesOfLinks[linkOfRemovedSpring][i] ==
          removedSpringIdx) {
        net.springIndicesOfLinks[linkOfRemovedSpring][i] = keptSpringIdx;
      }
    }
    this->removeDuplicateListedSpringsFromLink(net, linkOfRemovedSpring, true);
  }

  for (long int i =
         static_cast<long int>(net.springIndicesOfLinks[linkToReduce].size()) -
         1;
       i >= 0;
       --i) {
    const size_t iU = static_cast<size_t>(i);
    if (net.springIndicesOfLinks[linkToReduce][iU] == removedSpringIdx ||
        net.springIndicesOfLinks[linkToReduce][iU] == keptSpringIdx) {
      net.springIndicesOfLinks[linkToReduce].erase(
        net.springIndicesOfLinks[linkToReduce].begin() + i);
    }
  }
  this->removeDuplicateListedSpringsFromLink(net, linkToReduce);

  net.linkIndicesOfSprings.erase(net.linkIndicesOfSprings.begin() +
                                 static_cast<std::ptrdiff_t>(removedSpringIdx));
  // partial springs
  if (net.partialSpringIsPartial[removedPartialSpringIdxE]) {
    net.partialSpringIsPartial[remainingPartialSpringIdxE] = true;
  }
  const unsigned int removedPartialSpringIdxU =
    static_cast<unsigned int>(removedPartialSpringIdx);
  pylimer_tools::utils::removeRow(net.partialSpringIsPartial,
                                  removedPartialSpringIdxU);

  const bool removedIsA =
    net.springPartIndexA[removedPartialSpringIdxE] == linkToReduceI;
  const size_t otherEndOfRemovedSpring = static_cast<size_t>(
    removedIsA ? net.springPartIndexB[removedPartialSpringIdxE]
               : net.springPartIndexA[removedPartialSpringIdxE]);
  double offsetMultiplier = removedIsA ? -1. : 1.;
  if (net.springPartIndexA[remainingPartialSpringIdxE] == linkToReduceI) {
    net.springPartIndexA[remainingPartialSpringIdxE] =
      static_cast<int>(otherEndOfRemovedSpring);
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      net.springPartCoordinateIndexA[3 * remainingPartialSpringIdxE + dir] =
        static_cast<int>(3 * otherEndOfRemovedSpring) + static_cast<int>(dir);
    };
  } else {
    RUNTIME_EXP_IFN(
      net.springPartIndexB[remainingPartialSpringIdxE] == linkToReduceI, "");
    net.springPartIndexB[remainingPartialSpringIdxE] =
      static_cast<int>(otherEndOfRemovedSpring);
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      net.springPartCoordinateIndexB[3 * remainingPartialSpringIdxE + dir] =
        static_cast<int>(3 * otherEndOfRemovedSpring) + static_cast<int>(dir);
    }
    offsetMultiplier *= -1.;
  }
  net.springPartBoxOffset.segment(3 * remainingPartialSpringIdxE, 3) +=
    offsetMultiplier *
    net.springPartBoxOffset.segment(3 * removedPartialSpringIdxE, 3);
  pylimer_tools::utils::removeRow(net.springPartIndexA,
                                  removedPartialSpringIdxU);
  pylimer_tools::utils::removeRow(net.springPartIndexB,
                                  removedPartialSpringIdxU);
  pylimer_tools::utils::removeRows(
    net.springPartCoordinateIndexA, 3 * removedPartialSpringIdxU, 3);
  pylimer_tools::utils::removeRows(
    net.springPartCoordinateIndexB, 3 * removedPartialSpringIdxU, 3);
  pylimer_tools::utils::removeRows(
    net.springPartBoxOffset, 3 * removedPartialSpringIdxU, 3);

  // spring indices & coordinates
  if (net.springIndexA[removedSpringIdxE] == linkToReduceI) {
    if (net.springIndexA[keptSpringIdxE] == linkToReduceI) {
      net.springIndexA[keptSpringIdxE] = net.springIndexB[removedSpringIdxE];
      net.springCoordinateIndexA.segment(3 * keptSpringIdxE, 3) =
        net.springCoordinateIndexB.segment(3 * removedSpringIdxE, 3);
    } else {
      assert(net.springIndexB[keptSpringIdxE] == linkToReduceI);
      net.springIndexB[keptSpringIdxE] = net.springIndexB[removedSpringIdxE];
      net.springCoordinateIndexB.segment(3 * keptSpringIdxE, 3) =
        net.springCoordinateIndexB.segment(3 * removedSpringIdxE, 3);
    }
  } else {
    assert(net.springIndexB[removedSpringIdxE] == linkToReduceI);
    if (net.springIndexA[keptSpringIdxE] == linkToReduceI) {
      net.springIndexA[keptSpringIdxE] = net.springIndexA[removedSpringIdxE];
      net.springCoordinateIndexA.segment(3 * keptSpringIdxE, 3) =
        net.springCoordinateIndexA.segment(3 * removedSpringIdxE, 3);
    } else {
      assert(net.springIndexB[keptSpringIdxE] == linkToReduceI);
      net.springIndexB[keptSpringIdxE] = net.springIndexA[removedSpringIdxE];
      net.springCoordinateIndexB.segment(3 * keptSpringIdxE, 3) =
        net.springCoordinateIndexA.segment(3 * removedSpringIdxE, 3);
    }
  }
  const unsigned int removedSpringIdxU =
    static_cast<unsigned int>(removedSpringIdx);
  pylimer_tools::utils::removeRow(net.springIndexA, removedSpringIdxU);
  pylimer_tools::utils::removeRow(net.springIndexB, removedSpringIdxU);
  pylimer_tools::utils::removeRows(
    net.springCoordinateIndexA, 3 * removedSpringIdxU, 3);
  pylimer_tools::utils::removeRows(
    net.springCoordinateIndexB, 3 * removedSpringIdxU, 3);
  pylimer_tools::utils::removeRow(net.springIsActive, removedSpringIdxU);
  net.springToMoleculeIds.erase(net.springToMoleculeIds.begin() +
                                static_cast<std::ptrdiff_t>(removedSpringIdx));
  net.oldAtomIdToSpringIndex.erase(
    static_cast<size_t>(net.oldAtomIds[linkToReduceE]));

  pylimer_tools::utils::removeRow(net.partialToFullSpringIndex,
                                  removedPartialSpringIdxU);
  net.localToGlobalSpringIndex.erase(
    net.localToGlobalSpringIndex.begin() +
    static_cast<std::ptrdiff_t>(removedSpringIdx));
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
        net.loops[loopIdx].erase(net.loops[loopIdx].begin() +
                                 static_cast<std::ptrdiff_t>(i));
      }
      if (net.loops[loopIdx][i] > removedSpringIdx) {
        net.loops[loopIdx][i] -= 1;
      }
    }
  }

  // and the partial springs
  for (Eigen::Index i = 0; i < net.partialToFullSpringIndex.size(); ++i) {
    RUNTIME_EXP_IFN(net.partialToFullSpringIndex[i] !=
                      static_cast<int>(removedSpringIdx),
                    "");
    if (net.partialToFullSpringIndex[i] > static_cast<int>(removedSpringIdx)) {
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
  net.springsContourLength[keptSpringIdxE] +=
    net.springsContourLength[removedSpringIdxE];
  pylimer_tools::utils::removeRow(net.springsContourLength, removedSpringIdxU);
  pylimer_tools::utils::removeRow(net.springsType, removedSpringIdxU);
  RUNTIME_EXP_IFN(net.springsContourLength.size() ==
                    static_cast<Eigen::Index>(net.nrOfSprings),
                  "");
  // and spring partitions
  springPartitions[remainingPartialSpringIdxE] +=
    springPartitions[removedPartialSpringIdxE];

  pylimer_tools::utils::removeRow(springPartitions, removedPartialSpringIdxU);
  RUNTIME_EXP_IFN(springPartitions.size() ==
                    static_cast<Eigen::Index>(net.nrOfPartialSprings),
                  "");
  const size_t newKeptSpringIdx =
    (keptSpringIdx < removedSpringIdx) ? keptSpringIdx : (keptSpringIdx - 1);
  // admittedly, this is possibly dangerous, as it could hide
  // other mistakes
  const double newTotalForNormalization =
    springPartitions(net.localToGlobalSpringIndex[newKeptSpringIdx]).sum();
  for (const size_t globalPartSpringIndex :
       net.localToGlobalSpringIndex[newKeptSpringIdx]) {
    springPartitions[static_cast<Eigen::Index>(globalPartSpringIndex)] *=
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
  const size_t newPartialSpringIdx =
    remainingPartialSpringIdx -
    (remainingPartialSpringIdx > removedPartialSpringIdx ? 1 : 0);
  const Eigen::Vector3d newDistance = this->evaluatePartialSpringDistance(
    net, displacements, newPartialSpringIdx, this->is2D, false);
  RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(
                    newDistance, distanceBefore, 1e-5),
                  "After merging two springs, the overall distance "
                  "is not consistent. Expected distance " +
                    std::to_string(distanceBefore) + ", but got " +
                    std::to_string(newDistance) + " for spring " +
                    std::to_string(newPartialSpringIdx) + ".");

#ifndef NDEBUG
  // check that the ordering is correct
  for (size_t i = 0; i < net.localToGlobalSpringIndex[newKeptSpringIdx].size();
       ++i) {
    const Eigen::Index partialSpringIdx = static_cast<Eigen::Index>(
      net.localToGlobalSpringIndex[newKeptSpringIdx][i]);
    const size_t endA =
      static_cast<size_t>(net.springPartIndexA[partialSpringIdx]);
    const size_t endB =
      static_cast<size_t>(net.springPartIndexB[partialSpringIdx]);
    std::vector<size_t> linkIndices =
      net.linkIndicesOfSprings[newKeptSpringIdx];
    assert(endA == linkIndices[i]);
    assert(endB == linkIndices[i + 1]);
  }
#endif
}

/**
 * @brief Remove crosslinkers, springs and associated slip-links with the
 * scheme suggested by Andrei
 *
 * @param net
 * @param displacements
 * @param springPartitions
 * @param tolerance
 * @return size_t
 */
size_t
MEHPForceBalance::doRemovalAndreisWay(ForceBalanceNetwork& net,
                                      Eigen::VectorXd& displacements,
                                      Eigen::VectorXd& springPartitions,
                                      const double tolerance) const
{
  size_t numRemovedTotal = 0;
  size_t numRemovedInIteration = 0;
  do {
    numRemovedInIteration = 0;
    // do removal of f = 1
    // remove all crosslinkers that are 0- or 1-functional
    for (long int crosslinkIdx = static_cast<long int>(net.nrOfNodes) - 1;
         crosslinkIdx >= 0;
         --crosslinkIdx) {
      const size_t crosslinkIdxU = static_cast<size_t>(crosslinkIdx);
      if (net.springIndicesOfLinks[crosslinkIdxU].size() == 0 // f = 0
      ) {
        // std::cout << "Removing x-link " << crosslinkIdx << std::endl;
        this->removeLink(net, displacements, crosslinkIdxU);
        numRemovedInIteration += 1;
        // this->validateNetwork(net, displacements, springPartitions);
      }

      if ( // or f = 1, NOT primary loop
        net.springIndicesOfLinks[crosslinkIdxU].size() == 1 &&
        XOR(net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdxU][0]]
                                    [0] == crosslinkIdxU,
            pylimer_tools::utils::last(
              net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdxU]
                                                               [0]]) ==
              crosslinkIdxU)) {
        // need to first remove the spring
        this->removeSpring(net,
                           displacements,
                           springPartitions,
                           net.springIndicesOfLinks[crosslinkIdxU][0]);
        // to then remove the link
        this->removeLink(net, displacements, crosslinkIdxU);
        numRemovedInIteration += 1;
      }
    }
    numRemovedTotal += numRemovedInIteration;
  } while (numRemovedInIteration > 0);
  // then, replace f = 2
  this->removeTwofunctionalCrosslinks(net, displacements, springPartitions);
  // and remove all springs that are inactive
  size_t numSpringsRemoved = 0;
  for (long int springIdx = static_cast<long int>(net.nrOfSprings) - 1;
       springIdx >= 0;
       --springIdx) {
    const size_t springIdxU = static_cast<size_t>(springIdx);
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
        net.linkIndicesOfSprings[springIdxU].size() <= 2) {
      // remove
      this->removeSpring(net, displacements, springPartitions, springIdxU);
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
 * @param oneOverSpringPartitionUpperLimit
 */
size_t
MEHPForceBalance::addSlipLinkToPartialSpring(
  ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,
  const size_t splitPartialSpringIdx,
  const size_t slipLinkIdx,
  const double oneOverSpringPartitionUpperLimit) const
{
  const Eigen::Index splitPartialSpringIdxE =
    static_cast<Eigen::Index>(splitPartialSpringIdx);
  INVALIDARG_EXP_IFN(
    !net.linkIsSliplink[net.springPartIndexA[splitPartialSpringIdxE]] ||
      !net.linkIsSliplink[net.springPartIndexB[splitPartialSpringIdxE]],
    "Require at least one part to be a crosslink.");
  const size_t newPartialSpringIdx = net.nrOfPartialSprings;
  const Eigen::Index newPartialSpringIdxE =
    static_cast<Eigen::Index>(newPartialSpringIdx);
  net.nrOfPartialSprings += 1;
  const size_t relevantSpring =
    static_cast<size_t>(net.partialToFullSpringIndex[splitPartialSpringIdxE]);
  const double N =
    net.springsContourLength[static_cast<Eigen::Index>(relevantSpring)];
  const double minAlpha = (oneOverSpringPartitionUpperLimit > 0.)
                            ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
                            : 1e-9;
  INVALIDARG_EXP_IFN(APPROX_WITHIN(minAlpha, 0.0, 1.0, 1e-12),
                     "minAlpha must be within 0. and 1.");
  RUNTIME_EXP_IFN(
    minAlpha * (static_cast<double>(
                  net.localToGlobalSpringIndex[relevantSpring].size()) +
                1.) <=
      1.,
    "With this minimum alpha, the slip-link cannot be placed on this "
    "partial spring.");

  const size_t oldPartnerA =
    static_cast<size_t>(net.springPartIndexA[splitPartialSpringIdxE]);
  const size_t oldPartnerB =
    static_cast<size_t>(net.springPartIndexB[splitPartialSpringIdxE]);

  const Eigen::Vector3d distanceBefore = this->evaluatePartialSpringDistance(
    net, u, splitPartialSpringIdx, this->is2D, false);

  // std::cout << "Adding slip-link " << slipLinkIdx << " to spring "
  //           << relevantSpring << " (partial " << partialSpringIdx
  //           << ") with minAlpha = " << minAlpha << std::endl;
  // resize the structures
  const Eigen::Index nrOfPartialSpringsE =
    static_cast<Eigen::Index>(net.nrOfPartialSprings);
  springPartitions.conservativeResize(nrOfPartialSpringsE);
  assert(springPartitions.size() == nrOfPartialSpringsE);
  net.springPartIndexA.conservativeResize(nrOfPartialSpringsE);
  net.springPartIndexB.conservativeResize(nrOfPartialSpringsE);
  net.springPartCoordinateIndexA.conservativeResize(3 * nrOfPartialSpringsE);
  net.springPartCoordinateIndexB.conservativeResize(3 * nrOfPartialSpringsE);
  net.springPartBoxOffset.conservativeResize(3 * nrOfPartialSpringsE);
  net.partialToFullSpringIndex.conservativeResize(nrOfPartialSpringsE);
  net.partialSpringIsPartial.conservativeResize(nrOfPartialSpringsE);
  // add the new info
  net.partialSpringIsPartial[splitPartialSpringIdxE] = true;
  net.partialSpringIsPartial[newPartialSpringIdxE] = true;
  net.partialToFullSpringIndex[newPartialSpringIdxE] =
    static_cast<int>(relevantSpring);

  pylimer_tools::utils::addIfNotContained(net.springIndicesOfLinks[slipLinkIdx],
                                          relevantSpring);

  // slightly change numbering to keep the numbering of
  // localToGlobalSpringIndex constant. I.e., we want the
  // `newPartialSpringIdx` to correspond to the spring with the crosslink
  const bool forward =
    !net.linkIsSliplink[static_cast<Eigen::Index>(oldPartnerA)];
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
        static_cast<std::ptrdiff_t>(
          net.linkIndicesOfSprings[relevantSpring].size() - 1),
      slipLinkIdx);
    net.localToGlobalSpringIndex[relevantSpring].push_back(newPartialSpringIdx);
  }

  // rewire the springs
  net.springPartIndexB[splitPartialSpringIdxE] = static_cast<int>(slipLinkIdx);
  net.springPartIndexA[newPartialSpringIdxE] = static_cast<int>(slipLinkIdx);
  net.springPartIndexB[newPartialSpringIdxE] = static_cast<int>(oldPartnerB);

  for (Eigen::Index dir = 0; dir < 3; ++dir) {
    net.springPartCoordinateIndexA[3 * splitPartialSpringIdxE + dir] =
      3 * net.springPartIndexA[splitPartialSpringIdxE] + static_cast<int>(dir);
    net.springPartCoordinateIndexB[3 * splitPartialSpringIdxE + dir] =
      3 * net.springPartIndexB[splitPartialSpringIdxE] + static_cast<int>(dir);

    net.springPartCoordinateIndexA[3 * newPartialSpringIdxE + dir] =
      3 * net.springPartIndexA[newPartialSpringIdxE] + static_cast<int>(dir);
    net.springPartCoordinateIndexB[3 * newPartialSpringIdxE + dir] =
      3 * net.springPartIndexB[newPartialSpringIdxE] + static_cast<int>(dir);
  }

  // renormalize this spring
  // mostly by moving the next slip-link further
  springPartitions[newPartialSpringIdxE] = minAlpha;
  double remainingNormalisationOffset = minAlpha;
  if (forward) {
    for (const size_t globalPartSpringIndex :
         net.localToGlobalSpringIndex[relevantSpring]) {
      const Eigen::Index globalPartSpringIndexE =
        static_cast<Eigen::Index>(globalPartSpringIndex);
      const double currAlpha = springPartitions[globalPartSpringIndexE];
      if (currAlpha > minAlpha) {
        springPartitions[globalPartSpringIndexE] -=
          remainingNormalisationOffset;
        springPartitions[globalPartSpringIndexE] =
          std::max(springPartitions[globalPartSpringIndexE], minAlpha);
        remainingNormalisationOffset -=
          (currAlpha - springPartitions[globalPartSpringIndexE]);
      }
      if (remainingNormalisationOffset <= 0.) {
        break;
      }
    }
  } else {
    for (long int i = static_cast<long int>(
                        net.localToGlobalSpringIndex[relevantSpring].size()) -
                      1;
         i >= 0;
         --i) {
      const Eigen::Index globalPartSpringIndexE = static_cast<Eigen::Index>(
        net.localToGlobalSpringIndex[relevantSpring][static_cast<size_t>(i)]);
      const double currAlpha = springPartitions[globalPartSpringIndexE];
      if (currAlpha > minAlpha) {
        springPartitions[globalPartSpringIndexE] -=
          remainingNormalisationOffset;
        springPartitions[globalPartSpringIndexE] =
          std::max(springPartitions[globalPartSpringIndexE], minAlpha);
        remainingNormalisationOffset -=
          (currAlpha - springPartitions[globalPartSpringIndexE]);
      }
      if (remainingNormalisationOffset <= 0.) {
        break;
      }
    }
  }
  // check that normalisation worked
  const double newTotalForNormalization =
    springPartitions(net.localToGlobalSpringIndex[relevantSpring]).sum();
  RUNTIME_EXP_IFN(APPROX_EQUAL(newTotalForNormalization, 1.0, 1e-9), "");

  // TODO: this is problematic: instead, the decision has to be made, how
  // the previous/existing offset should be split, etc.
  net.springPartBoxOffset.segment(3 * newPartialSpringIdxE, 3) =
    Eigen::Vector3d::Zero();
  this->reAlignSlipLinkToImages(
    net, u, slipLinkIdx, splitPartialSpringIdx, newPartialSpringIdx);

  const Eigen::Vector3d distanceAfter =
    this->evaluatePartialSpringDistance(
      net, u, splitPartialSpringIdx, this->is2D, false) +
    this->evaluatePartialSpringDistance(
      net, u, newPartialSpringIdx, this->is2D, false);
  assert(
    pylimer_tools::utils::vector_approx_equal(distanceBefore, distanceAfter));

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
size_t
MEHPForceBalance::removeTwofunctionalCrosslinks(
  ForceBalanceNetwork& net,
  Eigen::VectorXd& displacements,
  Eigen::VectorXd& springPartitions) const
{
  size_t numRemoved = 0;
  for (long int crosslinkIdx = static_cast<long int>(net.nrOfNodes) - 1;
       crosslinkIdx >= 0;
       --crosslinkIdx) {
    const size_t crosslinkIdxU = static_cast<size_t>(crosslinkIdx);
    if (net.springIndicesOfLinks[crosslinkIdxU].size() == 2) {
      std::vector<size_t> springsToMerge =
        net.springIndicesOfLinks[crosslinkIdxU];

      // special case: this is an entanglement bead
      if (net.oldAtomTypes[crosslinkIdx] == this->entanglementType) {
        RUNTIME_EXP_IFN(
          net.springsType[static_cast<Eigen::Index>(springsToMerge[0])] !=
              this->entanglementType &&
            net.springsType[static_cast<Eigen::Index>(springsToMerge[1])] !=
              this->entanglementType,
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
          (XOR(net.linkIndicesOfSprings[springsToMerge[0]][0] == crosslinkIdxU,
               pylimer_tools::utils::last(
                 net.linkIndicesOfSprings[springsToMerge[0]]) ==
                 crosslinkIdxU)) &&
          (XOR(net.linkIndicesOfSprings[springsToMerge[1]][0] == crosslinkIdxU,
               pylimer_tools::utils::last(
                 net.linkIndicesOfSprings[springsToMerge[1]]) ==
                 crosslinkIdxU))) {
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
                           crosslinkIdxU);

        // this->validateNetwork(net, displacements, springPartitions);
        // std::cout << "Removing link " << crosslinkIdx << std::endl;
        this->removeLink(net, displacements, crosslinkIdxU);

        // std::cout << "Removed crosslink " << crosslinkIdx << std::endl;

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
double
MEHPForceBalance::updateSpringPartition(
  const ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,        /* gives the parametrization of N */
  Eigen::VectorXd& oneOverSpringPartitions, /* gives the parametrization of N */
  const size_t linkIdx,
  const double oneOverSpringPartitionUpperLimit,
  bool allowSlipLinksToPassEachOther) const
{
  // std::cout << "Updating spring partition " << linkIdx << " of "
  //           << net.nrOfNodes << " / " << net.nrOfLinks << " with limit
  //           "
  //           << oneOverSpringPartitionUpperLimit << std::endl;

  INVALIDARG_EXP_IFN(linkIdx < net.springIndicesOfLinks.size(),
                     "Link to update needs to be in the list");
  INVALIDARG_EXP_IFN(net.linkIsSliplink[static_cast<Eigen::Index>(linkIdx)],
                     "Only slip-links may slip along a spring, link " +
                       std::to_string(linkIdx) + " is not one. Network has " +
                       std::to_string(net.nrOfNodes) + " cross- of " +
                       std::to_string(net.nrOfLinks) + " links.");
  const std::vector<size_t>& springIndices = net.springIndicesOfLinks[linkIdx];
  assert(springIndices.size() == 1 || springIndices.size() == 2);
  const size_t numSprings =
    (springIndices.size() == 2 && springIndices[0] == springIndices[1])
      ? 1
      : springIndices.size();
  double residualNorm = 0.0;
  int residualNormContributions = 0;
  for (size_t sIdx = 0; sIdx < numSprings; ++sIdx) {
    const size_t springIndex = springIndices[sIdx];
    const std::vector<size_t>& springsPartners =
      net.linkIndicesOfSprings[springIndex];
    for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
         ++partner_idx) {
      if (springsPartners[partner_idx] == linkIdx) {
        const size_t currentSpringGlobalIdx =
          net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
        const size_t neighbourSpringGlobalIdx =
          net.localToGlobalSpringIndex[springIndex][partner_idx];
        // found position of this link in this spring
        // want to find the ideal value for
        // net.springPartitions[springIndex][partner_idx-1]
        // NOTE: The following is slightly problematic for primary loops!
        Eigen::Vector3d vecBack = this->evaluatePartialSpringDistanceFrom(
          net, u, currentSpringGlobalIdx, linkIdx);
        const double distanceBack = (vecBack.squaredNorm());
        Eigen::Vector3d vecForward = this->evaluatePartialSpringDistanceFrom(
          net, u, neighbourSpringGlobalIdx, linkIdx);
        const double distanceForward = vecForward.squaredNorm();
        double idealValue = 1. / (1. + sqrt(distanceForward / distanceBack));
        if (distanceBack <= 0.0) {
          idealValue = 0.0; // TODO: really?
        }
        const Eigen::Index currentSpringGlobalIdxE =
          static_cast<Eigen::Index>(currentSpringGlobalIdx);
        const Eigen::Index neighbourSpringGlobalIdxE =
          static_cast<Eigen::Index>(neighbourSpringGlobalIdx);
        const double currentS = springPartitions[currentSpringGlobalIdxE];
        const double nextS = springPartitions[neighbourSpringGlobalIdxE];
        const double N =
          net.springsContourLength[static_cast<Eigen::Index>(springIndex)];
        const double l = (currentS + nextS);
        if (oneOverSpringPartitionUpperLimit > 0.) {
          // TODO: sketch theory why this should/not be necessary!!!
          const double limit = std::clamp(
            1. / (oneOverSpringPartitionUpperLimit * (nextS + currentS) * (N)),
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
        const double newS = idealValue * l;
        const double idealValueM1 = (1. - idealValue);
        const double complementaryS = (1. - idealValue) * l;
        double localResidualNorm = 0.0;
        residualNormContributions += 2;
        if (idealValue > 0.0 && idealValue < 1.0) {
          const double idealValue2 = idealValue * idealValue;
          const double idealValueM12 = idealValueM1 * idealValueM1;
          // way too complicated expression to solve subtraction
          // truncation issues?
          localResidualNorm = (std::fma(distanceBack,
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
            std::to_string(newS) + " and " + std::to_string(complementaryS) +
            " with ideal = " + std::to_string(idealValue) + " of " +
            std::to_string(nextS + currentS) + " for link " +
            std::to_string(linkIdx) +
            ". Diff: " + std::to_string(1. - (newS + complementaryS)) + ".");
        RUNTIME_EXP_IFN(
          APPROX_EQUAL(nextS + currentS, newS + complementaryS, 1e-9),
          "Require nextS + currentS == newS + complementaryS, got " +
            std::to_string(nextS + currentS) + " vs. " +
            std::to_string(newS + complementaryS) + " from " +
            std::to_string(nextS) + " and " + std::to_string(currentS) + ", " +
            std::to_string(newS) + " and " + std::to_string(complementaryS) +
            ". Diff: " +
            std::to_string((nextS + currentS) - (newS + complementaryS)) + ".");
        RUNTIME_EXP_IFN(
          APPROX_WITHIN(nextS + currentS, 0., 1., 1e-9),
          "Require nextS + currentS to be within 0, 1, got " +
            std::to_string(nextS + currentS) + " from " +
            std::to_string(nextS) + " and " + std::to_string(currentS) +
            ". Diff: " + std::to_string(1. - (nextS + currentS)) + ".");
        RUNTIME_EXP_IFN(nextS >= -0.00000000000001,
                        "nextS must be >= 0., got " + std::to_string(nextS) +
                          " from " + std::to_string(nextS) + " and " +
                          std::to_string(currentS) + ", " +
                          std::to_string(newS) + " and " +
                          std::to_string(complementaryS) + ".");
        RUNTIME_EXP_IFN(
          complementaryS >= -0.00000000000001,
          "complementaryS must be >= 0., got " +
            std::to_string(complementaryS) + " from " + std::to_string(nextS) +
            " and " + std::to_string(currentS) + ", " + std::to_string(newS) +
            " and " + std::to_string(complementaryS) + ".");

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
        springPartitions[currentSpringGlobalIdxE] = newS;
        springPartitions[neighbourSpringGlobalIdxE] = complementaryS;
        if (oneOverSpringPartitions.size() > 0) {
          const double primaryCorrectionMultiplierC =
            static_cast<double>(net.springPartIndexA[currentSpringGlobalIdxE] !=
                                net.springPartIndexB[currentSpringGlobalIdxE]);
          double oneOverCurrent =
            primaryCorrectionMultiplierC *
            CLAMP_ONE_OVER_SPRINGPARTITION(
              net.partialSpringIsPartial[currentSpringGlobalIdxE],
              (1.0 / (newS * N)),
              N,
              oneOverSpringPartitionUpperLimit);
          oneOverSpringPartitions.segment(3 * currentSpringGlobalIdxE, 3) =
            Eigen::Vector3d::Constant(oneOverCurrent);
          const double primaryCorrectionMultiplierN = static_cast<double>(
            net.springPartIndexA[neighbourSpringGlobalIdxE] !=
            net.springPartIndexB[neighbourSpringGlobalIdxE]);
          double oneOverNeighbour =
            primaryCorrectionMultiplierN *
            CLAMP_ONE_OVER_SPRINGPARTITION(
              net.partialSpringIsPartial[neighbourSpringGlobalIdxE],
              1.0 / (complementaryS * N),
              N,
              oneOverSpringPartitionUpperLimit);
          oneOverSpringPartitions.segment(3 * neighbourSpringGlobalIdxE, 3) =
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
void
MEHPForceBalance::moveSlipLinksToTheirBestBranch(
  ForceBalanceNetwork& net,
  Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,
  const double oneOverSpringPartitionUpperLimit,
  const int nrOfCrosslinkSwapsAllowedPerSliplink,
  const bool respectLoops,
  const bool moveAttempt) const
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
void
MEHPForceBalance::moveSlipLinkToItsBestBranch(
  ForceBalanceNetwork& net,
  Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,
  const size_t slipLinkIdx,
  const double oneOverSpringPartitionUpperLimit,
  const int nrOfCrosslinkSwapsAllowedPerSliplink,
  const bool respectLoops,
  const bool moveAttempt) const
{
  INVALIDARG_EXP_IFN(net.linkIsSliplink[static_cast<Eigen::Index>(slipLinkIdx)],
                     "Passed slip-link must be one.");
  const std::vector<size_t> associatedSprings =
    net.springIndicesOfLinks[slipLinkIdx];
  // skip slip-links that are with its own spring, for now.
  if (associatedSprings.size() <= 1) {
    return;
  }

  for (const size_t springIdx : associatedSprings) {
    const double N =
      net.springsContourLength[static_cast<Eigen::Index>(springIdx)];
    const double swappableCutoff =
      (oneOverSpringPartitionUpperLimit > 0.)
        ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
        : 1e-12;
    for (size_t linkI = 1;
         linkI < net.linkIndicesOfSprings[springIdx].size() - 1;
         ++linkI) {
      if (net.linkIndicesOfSprings[springIdx][linkI] == slipLinkIdx) {
        // found index of this slip-link.
        const size_t partitionBeforeIdx =
          net.localToGlobalSpringIndex[springIdx][linkI - 1];
        const size_t partitionAfterIdx =
          net.localToGlobalSpringIndex[springIdx][linkI];
        assert(this->isPartOfSpring(net, slipLinkIdx, partitionBeforeIdx));
        assert(this->isPartOfSpring(net, slipLinkIdx, partitionAfterIdx));
        bool didSwap = false;
        // check whether swap is needed in either direction
        // swap if yes
        if (springPartitions[static_cast<Eigen::Index>(partitionBeforeIdx)] <=
            swappableCutoff) {
          didSwap =
            this->swapSlipLinkReversibly(net,
                                         u,
                                         springPartitions,
                                         partitionBeforeIdx,
                                         oneOverSpringPartitionUpperLimit,
                                         nrOfCrosslinkSwapsAllowedPerSliplink,
                                         respectLoops,
                                         moveAttempt);
        }
        if (springPartitions[static_cast<Eigen::Index>(partitionAfterIdx)] <=
              swappableCutoff &&
            !didSwap) {
          didSwap =
            this->swapSlipLinkReversibly(net,
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
bool
MEHPForceBalance::swapSlipLinkReversibly(
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
  const Eigen::Index partialSpringIdxE =
    static_cast<Eigen::Index>(partialSpringIdx);
  const size_t partnerA =
    static_cast<size_t>(net.springPartIndexA[partialSpringIdxE]);
  const size_t partnerB =
    static_cast<size_t>(net.springPartIndexB[partialSpringIdxE]);
  const Eigen::Index partnerAE = static_cast<Eigen::Index>(partnerA);
  const Eigen::Index partnerBE = static_cast<Eigen::Index>(partnerB);
  INVALIDARG_EXP_IFN(net.linkIsSliplink[partnerAE] ||
                       net.linkIsSliplink[partnerBE],
                     "Cannot swap crosslink with crosslink.");
  if (partnerA == partnerB) {
    return false;
  }
  const size_t fullSpringIdx =
    static_cast<size_t>(net.partialToFullSpringIndex[partialSpringIdxE]);
  // analyse spring
  const bool involvesCrosslink =
    (!net.linkIsSliplink[partnerAE] || !net.linkIsSliplink[partnerBE]);
  if (involvesCrosslink) {
    // first check if allowed.
    const size_t slipLinkIdx =
      net.linkIsSliplink[partnerAE] ? partnerA : partnerB;
    if ((nrOfCrosslinkSwapsAllowedPerSliplink < 0) ||
        (net.nrOfCrosslinkSwapsEndured[static_cast<Eigen::Index>(
           slipLinkIdx - net.nrOfNodes)] <
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
        const size_t xlinkIdx =
          net.linkIsSliplink[partnerAE] ? partnerB : partnerA;
        const size_t otherRailPart = this->getOtherRailPartialSpringIdx(
          net, partialSpringIdx, slipLinkIdx);
        Eigen::Vector3d otherRailDistance =
          this->evaluatePartialSpringDistance(net, u, otherRailPart);
        Eigen::Vector3d thisRailDistance =
          this->evaluatePartialSpringDistance(net, u, partialSpringIdx);

        // make sure the vectors are in the direction of the crosslink
        if (slipLinkIdx == partnerA) {
          otherRailDistance *= -1.;
          thisRailDistance *= -1.;
        }
#ifndef NDEBUG
        this->validateNetwork();
#endif

        bool found = false;

        const std::vector<size_t> partialSpringIndices =
          this->getPartialSpringIndicesOfLink(net, xlinkIdx);
        for (const size_t attemptedEdge : partialSpringIndices) {
          if (attemptedEdge == partialSpringIdx) {
            found = true;
            continue;
          }

          Eigen::Vector3d attemptSpringDistance =
            this->evaluatePartialSpringDistanceFrom(
              net, u, attemptedEdge, xlinkIdx);

          // TODO: involve denominators
          // we only look at a quasi force,
          const double forceEstimateBefore =
            // on the slip-link
            (-otherRailDistance + thisRailDistance).squaredNorm() +
            // and on the crosslink
            (attemptSpringDistance - thisRailDistance).squaredNorm();
          // and for the case that we did switch the slip-link onto the
          // attemptedEdge
          const double forceEstimateAfter =
            // force on the slip-link (-this + this), cancel out
            (attemptSpringDistance).squaredNorm() +
            // and on the crosslink
            (-2 * thisRailDistance - otherRailDistance).squaredNorm();

          if (forceEstimateAfter < forceEstimateBefore) {
            const long int newPartialSpringIdx =
              this->moveSlipLinkFromRailToRail(
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
        net.nrOfCrosslinkSwapsEndured[static_cast<Eigen::Index>(
          slipLinkIdx - net.nrOfNodes)] += 1;
      }
      return didSwap;
    }
    return false;
  } else {
    // does not involve crosslink
    // first, decide: do we attempt the move, or not?
    if (moveAttempt) {
      return this->swapSlipLinksReversibly(net,
                                           u,
                                           springPartitions,
                                           partialSpringIdx,
                                           oneOverSpringPartitionUpperLimit);
    } else {
      // check if energy is smaller
      const size_t indexInSpring = pylimer_tools::utils::index_of(
        net.localToGlobalSpringIndex[fullSpringIdx], partialSpringIdx);
      assert(indexInSpring > 0 &&
             indexInSpring <
               net.localToGlobalSpringIndex[fullSpringIdx].size() - 1);
      const size_t otherRailFrom =
        net.localToGlobalSpringIndex[fullSpringIdx][indexInSpring - 1];
      const size_t otherRailTo =
        net.localToGlobalSpringIndex[fullSpringIdx][indexInSpring + 1];

      const double thisPartialSpringDenominator =
        this->getDenominatorOfPartialSpring(net,
                                            springPartitions,
                                            partialSpringIdx,
                                            oneOverSpringPartitionUpperLimit);
      const Eigen::Vector3d thisSpringDistance =
        this->evaluatePartialSpringDistance(net, u, partialSpringIdx);
      const double otherRailFromDenominator =
        this->getDenominatorOfPartialSpring(net,
                                            springPartitions,
                                            otherRailFrom,
                                            oneOverSpringPartitionUpperLimit);
      const Eigen::Vector3d otherRailFromSpringDistance =
        this->evaluatePartialSpringDistance(net, u, otherRailFrom);

      const double otherRailToDenominator = this->getDenominatorOfPartialSpring(
        net, springPartitions, otherRailTo, oneOverSpringPartitionUpperLimit);
      const Eigen::Vector3d otherRailToSpringDistance =
        this->evaluatePartialSpringDistance(net, u, otherRailTo);

      const double forceEstimateBefore =
        // force estimate on `partnerA` (from)
        (-otherRailFromSpringDistance * otherRailFromDenominator +
         thisSpringDistance * thisPartialSpringDenominator)
          .squaredNorm() +
        // force estimate on `partnerB` (to)
        (-thisSpringDistance * thisPartialSpringDenominator +
         otherRailToSpringDistance * otherRailToDenominator)
          .squaredNorm();
      // force estimate if we do the mutation
      const double forceEstimateAfter =
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
bool
MEHPForceBalance::swapSlipLinkWithXlinkReversibly(
  ForceBalanceNetwork& net,
  Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,
  const size_t partialSpringIdx,
  const double oneOverSpringPartitionUpperLimit,
  const bool respectLoops) const
{
  const Eigen::Index partialSpringIdxE =
    static_cast<Eigen::Index>(partialSpringIdx);
  const size_t partnerA =
    static_cast<size_t>(net.springPartIndexA[partialSpringIdxE]);
  const size_t partnerB =
    static_cast<size_t>(net.springPartIndexB[partialSpringIdxE]);
  const Eigen::Index partnerAE = static_cast<Eigen::Index>(partnerA);
  const Eigen::Index partnerBE = static_cast<Eigen::Index>(partnerB);
  INVALIDARG_EXP_IFN(
    XOR(net.linkIsSliplink[partnerAE], net.linkIsSliplink[partnerBE]),
    "This method only swaps crosslink with slip-link.");
  if (partnerA == partnerB) {
    return false;
  }
  const size_t fullSpringIdx =
    static_cast<size_t>(net.partialToFullSpringIndex[partialSpringIdxE]);
  // analyse spring
  const size_t crosslinkIdx =
    net.linkIsSliplink[partnerAE] ? partnerB : partnerA;
  const size_t slipLinkIdx =
    net.linkIsSliplink[partnerBE] ? partnerB : partnerA;
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
  for (const size_t crosslinksSpringIdx :
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
  for (const size_t relevantNeighbour : relevantNeighbours) {
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
  long int newPartialSpringIdx = 0;
  // swap with crosslink
  const std::vector<size_t> springsOfCrosslink =
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
    newPartialSpringIdx = this->rotateSlipLinkAroundCrosslink(
      net,
      u,
      springPartitions,
      static_cast<size_t>(newPartialSpringIdx),
      oneOverSpringPartitionUpperLimit,
      respectLoops);
    isBackToInitialSpring = pylimer_tools::utils::contains(
      net.springIndicesOfLinks[slipLinkIdx], fullSpringIdx);
    for (size_t relaxSteps = 0; relaxSteps < 2; ++relaxSteps) {
      // TODO: this is not good at all.
      for (const size_t linkIdx : relevantNeighbours) {
        this->relaxationLight(
          net, springPartitions, u, linkIdx, oneOverSpringPartitionUpperLimit);
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
              << ". Crosslink is " << crosslinkIdx
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
    for (const size_t linkIdx : relevantNeighbours) {
      this->relaxationLight(
        net, springPartitions, u, linkIdx, oneOverSpringPartitionUpperLimit);
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
bool
MEHPForceBalance::swapSlipLinksReversibly(
  ForceBalanceNetwork& net,
  Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,
  const size_t partialSpringIdx,
  const double oneOverSpringPartitionUpperLimit) const
{
  const Eigen::Index partialSpringIdxE =
    static_cast<Eigen::Index>(partialSpringIdx);
  const size_t partnerA =
    static_cast<size_t>(net.springPartIndexA[partialSpringIdxE]);
  const size_t partnerB =
    static_cast<size_t>(net.springPartIndexB[partialSpringIdxE]);
  INVALIDARG_EXP_IFN(net.linkIsSliplink[static_cast<Eigen::Index>(partnerA)] &&
                       net.linkIsSliplink[static_cast<Eigen::Index>(partnerB)],
                     "This method only swaps slip-links.");
  if (partnerA == partnerB) {
    return false;
  }
  const size_t fullSpringIdx =
    static_cast<size_t>(net.partialToFullSpringIndex[partialSpringIdxE]);
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
  for (const size_t relevantNeighbour : relevantNeighbours) {
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
  const Eigen::VectorXd displacementsBefore = u(relevantNeighboursCoordIndices);
  const Eigen::VectorXd springPartitionsBefore =
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
  const double residualAfter =
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
 * @param partialSpringIdx1 one of the two partial spring idx
 * @param partialSpringIdx2 the partial spring idx of the other spring
 */
void
MEHPForceBalance::reAlignSlipLinkToImages(ForceBalanceNetwork& net,
                                          const Eigen::VectorXd& u,
                                          const size_t slipLinkIdx,
                                          const size_t partialSpringIdx1,
                                          const size_t partialSpringIdx2) const
{
  const Eigen::Index partialSpringIdx1E =
    static_cast<Eigen::Index>(partialSpringIdx1);
  const Eigen::Index partialSpringIdx2E =
    static_cast<Eigen::Index>(partialSpringIdx2);
  const Eigen::Index slipLinkIdxE = static_cast<Eigen::Index>(slipLinkIdx);
  assert(net.springPartIndexB[partialSpringIdx1E] ==
         static_cast<int>(slipLinkIdx));
  assert(net.springPartIndexA[partialSpringIdx2E] ==
         static_cast<int>(slipLinkIdx));
  assert(net.linkIsSliplink[slipLinkIdxE]);
  assert(net.partialToFullSpringIndex[partialSpringIdx1E] ==
         net.partialToFullSpringIndex[partialSpringIdx2E]);
  Eigen::Vector3d totalOffset =
    this->getPartialSpringBoxOffset(net, partialSpringIdx1) +
    this->getPartialSpringBoxOffset(net, partialSpringIdx2);
  const Eigen::Vector3d totalDistanceBefore =
    this->evaluatePartialSpringDistance(
      net, u, partialSpringIdx1, this->is2D, false) +
    this->evaluatePartialSpringDistance(
      net, u, partialSpringIdx2, this->is2D, false);

  Eigen::Vector3d sourceCoords =
    net.coordinates.segment(3 * net.springPartIndexA[partialSpringIdx1E], 3) +
    u.segment(3 * net.springPartIndexA[partialSpringIdx1E], 3);
  Eigen::Vector3d targetCoords =
    net.coordinates.segment(3 * net.springPartIndexB[partialSpringIdx2E], 3) +
    u.segment(3 * net.springPartIndexB[partialSpringIdx2E], 3);
  Eigen::Vector3d viaCoords = net.coordinates.segment(3 * slipLinkIdxE, 3) +
                              u.segment(3 * slipLinkIdxE, 3);

  // std::cout << net.springPartIndexA[partialSpringIdx1] << ": "
  //           << sourceCoords << " to "
  //           << net.springPartIndexB[partialSpringIdx2] << ": "
  //           << targetCoords << " via " << viaCoords << std::endl;

  double bestOffsetScore = -1.;
  Eigen::Vector3d bestOffset = Eigen::Vector3d::Zero();

  // ugly brute-force method to check all possible combinations (ideally,
  // more or less at least)
  const Eigen::Array3i multiplicity1 =
    ((this->box.getOffset(viaCoords - sourceCoords).array().abs() +
      this->getPartialSpringBoxOffset(net, partialSpringIdx1).array().abs()) /
     this->box.getL())
      .rint()
      .cast<int>()
      .abs();
  const Eigen::Array3i multiplicity2 =
    ((this->box.getOffset(targetCoords - viaCoords).array() +
      this->getPartialSpringBoxOffset(net, partialSpringIdx2).array().abs()) /
     this->box.getL())
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

        const double currentScore = vec1.squaredNorm() + vec2.squaredNorm();
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
  net.springPartBoxOffset.segment(3 * partialSpringIdx1E, 3) = bestOffset;
  net.springPartBoxOffset.segment(3 * partialSpringIdx2E, 3) =
    totalOffset - bestOffset;
  const Eigen::Vector3d totalDistanceNow =
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
void
MEHPForceBalance::relaxationLight(
  ForceBalanceNetwork& net,
  Eigen::VectorXd& springPartitions,
  Eigen::VectorXd& oneOverSpringPartitions,
  Eigen::VectorXd& u,
  const size_t linkIdx,
  const double oneOverSpringPartitionUpperLimit) const
{
  if (net.linkIsSliplink[static_cast<Eigen::Index>(linkIdx)]) {
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
void
MEHPForceBalance::swapSlipLinksInclXlinks(
  ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,
  const double oneOverSpringPartitionUpperLimit,
  const bool respectLoops) const
{
#ifndef NDEBUG
  this->validateNetwork(net, u, springPartitions);
#endif

  for (size_t springIdx = 0; springIdx < net.nrOfSprings; ++springIdx) {
    if (net.linkIndicesOfSprings[springIdx].size() <= 2) {
      // no need to handle springs without slip-links
      continue;
    }

    const double N =
      net.springsContourLength[static_cast<Eigen::Index>(springIdx)];
    const double swappableCutoff =
      (oneOverSpringPartitionUpperLimit > 0.)
        ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
        : 1e-12;

    // loop the remaining partial springs
    // NOTE: it is slightly problematc, that e.g.
    // net.localToGlobalSpringIndex changes /!\  no idea ye how to easily
    // compensate that...
    for (long int partialIdx =
           static_cast<long int>(
             net.localToGlobalSpringIndex[springIdx].size()) -
           1;
         partialIdx >= 0;
         --partialIdx) {
      const size_t partialIdxU = static_cast<size_t>(partialIdx);
      // check if they qualify for swapping
      if (springPartitions[static_cast<Eigen::Index>(
            net.localToGlobalSpringIndex[springIdx][partialIdxU])] <=
          swappableCutoff) {
        const size_t partialSpringIdx =
          net.localToGlobalSpringIndex[springIdx][partialIdxU];
        const Eigen::Index partialSpringIdxE =
          static_cast<Eigen::Index>(partialSpringIdx);
        if (net.springPartIndexA[partialSpringIdxE] !=
            net.springPartIndexB[partialSpringIdxE]) {
          // do the swap
          if (partialIdxU == 0 ||
              partialIdxU ==
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
            //           << " around crosslink " << involvedCrosslink
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
void
MEHPForceBalance::swapSlipLinks(
  ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,
  const double oneOverSpringPartitionUpperLimit) const
{
  for (size_t springIdx = 0; springIdx < net.nrOfSprings; ++springIdx) {
    if (net.linkIndicesOfSprings[springIdx].size() <= 3) {
      // no need to handle springs with 1 or less slip-links
      continue;
    }

    const double N =
      net.springsContourLength[static_cast<Eigen::Index>(springIdx)];
    const double swappableCutoff =
      (oneOverSpringPartitionUpperLimit > 0.)
        ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
        : 1e-12;

    // loop the remaining partial springs
    for (size_t partialIdx = 1;
         partialIdx < net.localToGlobalSpringIndex[springIdx].size() - 1;
         ++partialIdx) {
      // check if they qualify for swapping
      if (springPartitions[static_cast<Eigen::Index>(
            net.localToGlobalSpringIndex[springIdx][partialIdx])] <=
          swappableCutoff) {
        const size_t partialSpringIdx =
          net.localToGlobalSpringIndex[springIdx][partialIdx];
        const Eigen::Index partialSpringIdxE =
          static_cast<Eigen::Index>(partialSpringIdx);
        if (net.springPartIndexA[partialSpringIdxE] !=
            net.springPartIndexB[partialSpringIdxE]) {
          // do the swap
          this->swapSlipLinks(net, partialSpringIdx);
        }
      }
    }
  }
  this->validateNetwork(net, u, springPartitions);
}

/**
 * @brief Move a slip-link from one spring attached to a crosslink to
 * another spring attached to the same crosslink
 *
 * Returns the idx of the new partial spring that had to be introduced.
 * Returns a negative idx if the move was illegal or impossible.
 *
 * @param net
 * @param u
 * @param springPartitions
 * @param partialSpringIdx
 */
long int
MEHPForceBalance::rotateSlipLinkAroundCrosslink(
  ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  Eigen::VectorXd& springPartitions,
  const size_t partialSpringIdx,
  const double oneOverSpringPartitionUpperLimit,
  const bool respectLoops) const
{
  const Eigen::Index partialSpringIdxE =
    static_cast<Eigen::Index>(partialSpringIdx);
  INVALIDARG_EXP_IFN(net.springPartIndexA[partialSpringIdxE] !=
                       net.springPartIndexB[partialSpringIdxE],
                     "One of the two ends of the partial spring must be a "
                     "slip-link, one a crosslink");
  INVALIDARG_EXP_IFN(net.springPartIndexA[partialSpringIdxE] !=
                       net.springPartIndexB[partialSpringIdxE],
                     "Cannot rotate");
  // assemble required data
  size_t springIdx =
    static_cast<size_t>(net.partialToFullSpringIndex[partialSpringIdxE]);
  RUNTIME_EXP_IFN(
    net.localToGlobalSpringIndex[springIdx][0] == partialSpringIdx ||
      pylimer_tools::utils::last(net.localToGlobalSpringIndex[springIdx]) ==
        partialSpringIdx,
    "Partial spring assembly is not correct");
  const double N =
    net.springsContourLength[static_cast<Eigen::Index>(springIdx)];
  const double swappableCutoff =
    (oneOverSpringPartitionUpperLimit > 0.)
      ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
      : 1e-9;
  const int maxNrOfSliplinksOnSpring = static_cast<int>(1.0 / swappableCutoff);
  const size_t partialIdx =
    net.localToGlobalSpringIndex[springIdx][0] == partialSpringIdx
      ? 0
      : net.localToGlobalSpringIndex[springIdx].size() - 1;
  // decide on the involved parties
  const size_t otherInvolvedPartialSpring =
    partialIdx == 0
      ? net.localToGlobalSpringIndex[springIdx][1]
      : net.localToGlobalSpringIndex
          [springIdx][net.localToGlobalSpringIndex[springIdx].size() - 2];
  size_t involvedSlipLink =
    net.linkIndicesOfSprings
      [springIdx]
      [partialIdx == 0 ? 1 : net.linkIndicesOfSprings[springIdx].size() - 2];
  const size_t involvedCrosslink =
    net.linkIndicesOfSprings
      [springIdx]
      [partialIdx == 0 ? 0 : net.linkIndicesOfSprings[springIdx].size() - 1];
  std::vector<size_t> possibleTargetPartialSprings;
  // find possible target partial springs
  std::vector<size_t> possibleTargetSprings =
    net.springIndicesOfLinks[involvedCrosslink];
  if (respectLoops && net.loopsOfSliplink.size() > 0) {
    // filter out the target springs that may not be a target based on the
    // involved loops
    std::erase_if(
      possibleTargetSprings,
      [&net, springIdx, involvedSlipLink](const size_t springIdxToCheck) {
        // we want to be able to distinguish all associated loops
        // into being (a) the ones the slip-link is slipping on, or
        // (b) the one the second, currently non-slipping part is
        // associated with only the cases of (a) are allowed as possible
        // target springs
        for (const size_t loopIdx : net.loopsOfSliplink[involvedSlipLink]) {
          bool loopContainsCheck = false;
          bool loopContainsOriginal = false;
          for (const size_t lspringIdx : net.loops[loopIdx]) {
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
    // std::cerr << "Spring " << springIdx << "'s crosslink " <<
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
      net.linkIndicesOfSprings[possibleTargetSprings[i]][0] == involvedCrosslink
        ? net.localToGlobalSpringIndex[possibleTargetSprings[i]][0]
        : pylimer_tools::utils::last(
            net.localToGlobalSpringIndex[possibleTargetSprings[i]]);
    if (currentPossibleTargetPartialSpringIdx != partialSpringIdx &&
        static_cast<Eigen::Index>(
          net.localToGlobalSpringIndex[possibleTargetSprings[i]].size()) <
          static_cast<Eigen::Index>(
            maxNrOfSliplinksOnSpring)) { // let's not combine with the
      // to-be-removed partial spring
      possibleTargetPartialSprings.push_back(
        currentPossibleTargetPartialSpringIdx);
    }
    if (currentPossibleTargetPartialSpringIdx == partialSpringIdx ||
        currentPossibleTargetPartialSpringIdx == otherInvolvedPartialSpring) {
      // we unfortunately cannot assert this due to primary loops
      // assert(currentPartialSpringTargetIdx == -1);
      currentPartialSpringTargetIdx = static_cast<int>(i);
    }
  }
  if (possibleTargetPartialSprings.size() == 0) {
    return -1;
  }

  const size_t targetPartialSpringIdx =
    possibleTargetPartialSprings[static_cast<size_t>(
                                   currentPartialSpringTargetIdx) %
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
void
MEHPForceBalance::swapSlipLinks(ForceBalanceNetwork& net,
                                const size_t partialSpringIdx) const
{
  const Eigen::Index partialSpringIdxE =
    static_cast<Eigen::Index>(partialSpringIdx);
  const size_t linkIdx1 =
    static_cast<size_t>(net.springPartIndexA[partialSpringIdxE]);
  const size_t linkIdx2 =
    static_cast<size_t>(net.springPartIndexB[partialSpringIdxE]);
  INVALIDARG_EXP_IFN(linkIdx1 != linkIdx2,
                     "Cannot swap link with itself: got " +
                       std::to_string(linkIdx1) + " and " +
                       std::to_string(linkIdx2) + ".");
  INVALIDARG_EXP_IFN(
    net.linkIsSliplink[static_cast<Eigen::Index>(linkIdx1)],
    "Only partial springs with only slip-links allow swapping.");
  INVALIDARG_EXP_IFN(
    net.linkIsSliplink[static_cast<Eigen::Index>(linkIdx2)],
    "Only partial springs with only slip-links allow swapping.");
  const size_t springIdx =
    static_cast<size_t>(net.partialToFullSpringIndex[partialSpringIdxE]);

  const long int firstPositionInSpring =
    static_cast<long int>(pylimer_tools::utils::index_of(
      net.localToGlobalSpringIndex[springIdx], partialSpringIdx));
  assert(
    firstPositionInSpring > 0 &&
    firstPositionInSpring <
      static_cast<long int>(net.localToGlobalSpringIndex[springIdx].size()) -
        1);

  // find the rest of the connectivity required for swapping
  const size_t firstPositionInSpringU =
    static_cast<size_t>(firstPositionInSpring);
  const size_t otherPartialOfLinkIdx1 =
    net.localToGlobalSpringIndex[springIdx][firstPositionInSpringU - 1];
  const size_t otherPartialOfLinkIdx2 =
    net.localToGlobalSpringIndex[springIdx][firstPositionInSpringU + 1];
  const Eigen::Index otherPartialOfLinkIdx1E =
    static_cast<Eigen::Index>(otherPartialOfLinkIdx1);
  const Eigen::Index otherPartialOfLinkIdx2E =
    static_cast<Eigen::Index>(otherPartialOfLinkIdx2);
  assert(otherPartialOfLinkIdx1 != otherPartialOfLinkIdx2);

  // std::cout << "Swapping link " << linkIdx1 << " and " << linkIdx2
  //           << " on partial spring " << partialSpringIdx
  //           << ". Newly linked partial springs: " <<
  //           otherPartialOfLinkIdx1
  //           << " (" << unaffectedEnd1 << ") "
  //           << " and " << otherPartialOfLinkIdx2 << " (" <<
  //           unaffectedEnd2
  //           << ") " << std::endl;

  const Eigen::VectorXd u = Eigen::VectorXd::Zero(net.coordinates.size());
  const Eigen::Vector3d distanceBefore =
    this->evaluatePartialSpringDistance(
      net, u, otherPartialOfLinkIdx1, this->is2D, false) +
    this->evaluatePartialSpringDistance(
      net, u, otherPartialOfLinkIdx2, this->is2D, false) +
    this->evaluatePartialSpringDistance(
      net, u, partialSpringIdx, this->is2D, false);

  RUNTIME_EXP_IFN(otherPartialOfLinkIdx1 != otherPartialOfLinkIdx2,
                  "Required assumption not met.");
  RUNTIME_EXP_IFN(
    firstPositionInSpring <
      static_cast<long int>(net.linkIndicesOfSprings[springIdx].size()) - 1,
    "Required assumption not met.");

  // actually do the swapping
  assert(net.springPartIndexB[otherPartialOfLinkIdx1E] ==
         static_cast<int>(linkIdx1));
  assert(net.springPartIndexA[otherPartialOfLinkIdx2E] ==
         static_cast<int>(linkIdx2));
  // re-link
  net.springPartIndexB[otherPartialOfLinkIdx1E] = static_cast<int>(linkIdx2);
  net.springPartIndexA[otherPartialOfLinkIdx2E] = static_cast<int>(linkIdx1);
  net.springPartCoordinateIndexB.segment(3 * otherPartialOfLinkIdx1E, 3) =
    Eigen::ArrayXi::LinSpaced(
      3, static_cast<int>(3 * linkIdx2), static_cast<int>(3 * linkIdx2 + 2));
  net.springPartCoordinateIndexA.segment(3 * otherPartialOfLinkIdx2E, 3) =
    Eigen::ArrayXi::LinSpaced(
      3, static_cast<int>(3 * linkIdx1), static_cast<int>(3 * linkIdx1 + 2));

  // update box offset
  net.springPartBoxOffset.segment(3 * otherPartialOfLinkIdx1E, 3) +=
    this->getPartialSpringBoxOffset(net, partialSpringIdx);
  net.springPartBoxOffset.segment(3 * otherPartialOfLinkIdx2E, 3) +=
    this->getPartialSpringBoxOffset(net, partialSpringIdx);

  // actually change direction
  std::swap(net.springPartIndexA[partialSpringIdxE],
            net.springPartIndexB[partialSpringIdxE]);
  for (Eigen::Index dir = 0; dir < 3; ++dir) {
    std::swap(net.springPartCoordinateIndexA[3 * partialSpringIdxE + dir],
              net.springPartCoordinateIndexB[3 * partialSpringIdxE + dir]);
  }
  net.springPartBoxOffset.segment(3 * partialSpringIdxE, 3) *= -1.;

  assert(net.linkIndicesOfSprings[springIdx][firstPositionInSpringU] ==
           linkIdx1 &&
         net.linkIndicesOfSprings[springIdx][firstPositionInSpringU + 1] ==
           linkIdx2);
  net.linkIndicesOfSprings[springIdx][firstPositionInSpringU] = linkIdx2;
  net.linkIndicesOfSprings[springIdx][firstPositionInSpringU + 1] = linkIdx1;

  // finally, validate result
  const Eigen::Vector3d distanceAfter =
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
 * @brief Displace all links to the mean of all connected neighbours
 *
 * @param net the force balance network
 * @param u the current displacements, wherein the resulting coordinates
 * shall be stored
 * @param oneOverSpringPartitions the inverse of the spring partitions
 * @return double, the distance (squared norm) displaced
 */
double
MEHPForceBalance::displaceToMeanPosition(
  const ForceBalanceNetwork& net,
  Eigen::VectorXd& u,
  const Eigen::ArrayXd& oneOverSpringPartitions) const
{
  assert(oneOverSpringPartitions.size() ==
         static_cast<Eigen::Index>(net.nrOfPartialSprings * 3));
  Eigen::ArrayXd objectiveDisplacement =
    Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
  const Eigen::ArrayXd partialSpringDistances =
    this
      ->evaluatePartialSpringVectors(
        net, u, this->is2D, this->assumeBoxLargeEnough)
      .array();
  objectiveDisplacement(net.springPartCoordinateIndexA) +=
    (oneOverSpringPartitions * partialSpringDistances);
  objectiveDisplacement(net.springPartCoordinateIndexB) -=
    (oneOverSpringPartitions * partialSpringDistances);

  Eigen::ArrayXd springPartWeightingFactor =
    Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(net.nrOfLinks * 3));
  const Eigen::ArrayXd loopPartialSpringEliminator =
    (net.springPartCoordinateIndexA != net.springPartCoordinateIndexB)
      .cast<double>();

  springPartWeightingFactor(net.springPartCoordinateIndexA) +=
    oneOverSpringPartitions * loopPartialSpringEliminator;
  springPartWeightingFactor(net.springPartCoordinateIndexB) +=
    oneOverSpringPartitions * loopPartialSpringEliminator;
  springPartWeightingFactor = springPartWeightingFactor.unaryExpr(
    [](const double v) { return v > 0. ? v : 1.0; });
  Eigen::ArrayXd remainingDisplacement =
    (objectiveDisplacement / springPartWeightingFactor);
#ifndef NDEBUG
  RUNTIME_EXP_IFN(
    pylimer_tools::utils::all_components_finite(remainingDisplacement),
    "Some displacements are not finite");
#endif
  // at this point, we have the ideal displacement if we were to do it
  // just one link at a time.
  // by doing all at once, as here, though, e.g. a pair of links would
  // oscillate back and forth to compensate for that:

  // NOTE: this stays mostly static, could be stored on the network
  Eigen::ArrayXd nSpringsPerLink =
    Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(net.nrOfLinks * 3));
  // add a one for every partial spring that's not a primary loop
  nSpringsPerLink(net.springPartCoordinateIndexA) +=
    loopPartialSpringEliminator;
  nSpringsPerLink(net.springPartCoordinateIndexB) +=
    loopPartialSpringEliminator;
  nSpringsPerLink =
    nSpringsPerLink.unaryExpr([](const double v) { return v > 0. ? v : 1.0; });
  // make sure there are no infinite back-and-forth
  // and actually displace
  Eigen::ArrayXd backForthDisplacement =
    Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(net.nrOfLinks * 3));
  backForthDisplacement(net.springPartCoordinateIndexA) +=
    loopPartialSpringEliminator *
    (remainingDisplacement(net.springPartCoordinateIndexB) /
     (nSpringsPerLink(net.springPartCoordinateIndexA) * 2.));
  backForthDisplacement(net.springPartCoordinateIndexB) +=
    loopPartialSpringEliminator *
    (remainingDisplacement(net.springPartCoordinateIndexA) /
     (nSpringsPerLink(net.springPartCoordinateIndexB) * 2.));
#ifndef NDEBUG
  RUNTIME_EXP_IFN(
    pylimer_tools::utils::all_components_finite(backForthDisplacement),
    "Some displacements are not finite");
#endif

  // actually displace
  const Eigen::VectorXd finalDisplacement =
    (remainingDisplacement + backForthDisplacement).matrix();
  RUNTIME_EXP_IFN(pylimer_tools::utils::all_components_finite<Eigen::VectorXd>(
                    finalDisplacement),
                  "Some displacements are not finite");
  // this->box.handlePBC(finalDisplacement);
  u += finalDisplacement;

  const double max_disp =
    pylimer_tools::utils::segmentwise_norm_max(finalDisplacement, 3);

  return max_disp;
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
double
MEHPForceBalance::displaceToMeanPosition(
  const ForceBalanceNetwork& net,
  Eigen::VectorXd& u,
  const Eigen::VectorXd& springPartitions,
  const size_t linkIdx,
  const double oneOverSpringPartitionUpperLimit) const
{
  return this->displaceToMeanPosition(
    net,
    u,
    springPartitions,
    linkIdx,
    oneOverSpringPartitionUpperLimit,
    this->getPartialSpringIndicesOfLink(net, linkIdx));
}

double
MEHPForceBalance::displaceToMeanPosition(
  const ForceBalanceNetwork& net,
  Eigen::VectorXd& u,
  const Eigen::VectorXd& springPartitions,
  const size_t linkIdx,
  const double oneOverSpringPartitionUpperLimit,
  const std::vector<size_t>& partialSpringIndices) const
{
#ifndef NDEBUG
  const Eigen::Vector3d forceBefore = this->getForceOn(
    net, u, springPartitions, linkIdx, oneOverSpringPartitionUpperLimit);
#endif
  // Eigen::Vector3d currentDisplacement = u.segment(3 * linkIdx, 3);
  Eigen::Vector3d objectiveDisplacement =
    Eigen::Vector3d::Zero(); // = remainingDisplacement.array();
  double objectiveDisplacementContributors = 0.0;

  const int linkIdxI = static_cast<int>(linkIdx);
  for (const size_t globalSpringIndex : partialSpringIndices) {
    const Eigen::Index globalSpringIndexE =
      static_cast<Eigen::Index>(globalSpringIndex);
    assert(net.springPartIndexA[globalSpringIndexE] == linkIdxI ||
           net.springPartIndexB[globalSpringIndexE] == linkIdxI);
    if (net.springPartIndexA[globalSpringIndexE] == linkIdxI &&
        net.springPartIndexB[globalSpringIndexE] == linkIdxI) {
      // skip primary loops
      continue;
    }
    Eigen::Vector3d partialDistance = this->evaluatePartialSpringDistanceFrom(
      net, u, globalSpringIndex, linkIdx);
    // std::cout << "Partial distance of " << globalSpringIndex << " from
    // "
    //           << linkIdx << " to "
    //           << this->getOtherSpringIndex(net, globalSpringIndex,
    //           linkIdx)
    //           << ": " << partialDistance << std::endl;
    // add to displacement
    const double contourLengthFraction = springPartitions[globalSpringIndexE];
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
      net
        .springsContourLength[net.partialToFullSpringIndex[globalSpringIndexE]];
    double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
      net.partialSpringIsPartial[globalSpringIndexE],
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
  u.segment(3 * static_cast<Eigen::Index>(linkIdx), 3) +=
    objectiveDisplacement * denominator;

#ifndef NDEBUG
  if (!this->assumeBoxLargeEnough) {
    const Eigen::Vector3d forceAfter = this->getForceOn(
      net, u, springPartitions, linkIdx, oneOverSpringPartitionUpperLimit);

    // this is only true if we don't have "full" PBC
    assert((pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
      forceAfter, Eigen::Vector3d::Zero(), 0.01)));
    if (!pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
          forceBefore, Eigen::Vector3d::Zero(), 0.01)) {
      assert(forceBefore.squaredNorm() >= forceAfter.squaredNorm());
    }
  }
#endif

  const double dist = (objectiveDisplacement * denominator).squaredNorm();
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
Eigen::Matrix3d
MEHPForceBalance::evaluateStressOnLink(
  const size_t linkIdx,
  const ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  const Eigen::VectorXd& springPartitions,
  Eigen::VectorXi& debugNrSpringsVisited,
  const double oneOverSpringPartitionUpperLimit) const
{
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

  const std::vector<size_t> partialSpringIndices =
    this->getPartialSpringIndicesOfLink(net, linkIdx);

  for (const size_t globalSpringIndex : partialSpringIndices) {
    const Eigen::Index globalSpringIndexE =
      static_cast<Eigen::Index>(globalSpringIndex);
    Eigen::Vector3d partialDistance = this->evaluatePartialSpringDistanceFrom(
      net, u, globalSpringIndex, linkIdx);
    const double contourLengthFraction = springPartitions[globalSpringIndexE];
    const double N =
      net
        .springsContourLength[net.partialToFullSpringIndex[globalSpringIndexE]];
    const double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
      net.partialSpringIsPartial[globalSpringIndexE],
      1.0 / (N * contourLengthFraction),
      N,
      oneOverSpringPartitionUpperLimit);

    double multiplier = this->kappa * oneOverContourLengthFraction;

    stress += multiplier * partialDistance * partialDistance.transpose();
    debugNrSpringsVisited[globalSpringIndexE] += 1;

    // also account for primary loops.
    // they may have non-zero length thanks to assuming the box is not
    // large enough...
    if (net.springPartIndexA[globalSpringIndexE] ==
        net.springPartIndexB[globalSpringIndexE]) {
      stress +=
        multiplier * (-partialDistance) * (-partialDistance).transpose();

      debugNrSpringsVisited[globalSpringIndexE] += 1;
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
Eigen::Vector3d
MEHPForceBalance::evaluateForceOnLink(
  const size_t linkIdx,
  const ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  const Eigen::VectorXd& springPartitions,
  Eigen::VectorXi& debugNrSpringsVisited,
  const double oneOverSpringPartitionUpperLimit) const
{
  Eigen::Vector3d force = Eigen::Vector3d::Zero();

  const std::vector<size_t> partialSpringIndices =
    this->getPartialSpringIndicesOfLink(net, linkIdx);

  const int linkIdxI = static_cast<int>(linkIdx);
  for (const size_t globalSpringIndex : partialSpringIndices) {
    const Eigen::Index globalSpringIndexE =
      static_cast<Eigen::Index>(globalSpringIndex);
    // partial spring's force goes both ways -> is zero anyway
    // but, as it would not be included twice in the list,
    // we have to skip them
    if (net.springPartIndexA[globalSpringIndexE] == linkIdxI &&
        net.springPartIndexB[globalSpringIndexE] == linkIdxI) {
      if (debugNrSpringsVisited.size() > 0) {
        debugNrSpringsVisited[globalSpringIndexE] += 2;
      }
      continue;
    }
    Eigen::Vector3d partialDistance = this->evaluatePartialSpringDistanceFrom(
      net, u, globalSpringIndex, linkIdx);
    const double contourLengthFraction = springPartitions[globalSpringIndexE];
    const double N =
      net
        .springsContourLength[net.partialToFullSpringIndex[globalSpringIndexE]];
    const double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
      net.partialSpringIsPartial[globalSpringIndexE],
      1.0 / (N * contourLengthFraction),
      N,
      oneOverSpringPartitionUpperLimit);

    force += this->kappa * oneOverContourLengthFraction * partialDistance;
    if (debugNrSpringsVisited.size() > 0) {
      debugNrSpringsVisited[globalSpringIndexE] += 1;
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
Eigen::VectorXd
MEHPForceBalance::evaluateSpringVectors(const ForceBalanceNetwork& net,
                                        const Eigen::VectorXd& u,
                                        const bool use2D,
                                        const bool assumeLarge) const
{
  assert(u.size() == net.coordinates.size());

  Eigen::VectorXd springVectors =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(net.nrOfSprings * 3));

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
      this->evaluatePartialSpringVectors(net, u, use2D, assumeLarge);

    // CAREFUL: this assumes a certain direction
    for (Eigen::Index i = 0;
         i < static_cast<Eigen::Index>(net.nrOfPartialSprings);
         ++i) {
      springVectors.segment(3 * net.partialToFullSpringIndex[i], 3) +=
        partialSpringVectors.segment(3 * i, 3);
    }
  }

  // reset for 2D systems
  if (use2D && net.nrOfSprings > 0) {
    // springDistances(Eigen::seq(2, net.nrOfSprings, 3)) =
    //   Eigen::VectorXd::Zero(net.nrOfSprings);
    for (Eigen::Index i = 2; i < static_cast<Eigen::Index>(3 * net.nrOfSprings);
         i += 3) {
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
Eigen::VectorXd
MEHPForceBalance::evaluateSpringLengths(const ForceBalanceNetwork& net,
                                        const Eigen::VectorXd& u,
                                        const bool use2D) const
{
  // first, the distances
  assert(u.size() == net.coordinates.size());

  Eigen::VectorXd partialSpringVectors =
    this->evaluatePartialSpringVectors(net, u, use2D, true);
  assert(partialSpringVectors.size() ==
         static_cast<Eigen::Index>(net.nrOfPartialSprings * 3));

  Eigen::VectorXd springLengths =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(net.nrOfSprings));
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(net.nrOfPartialSprings);
       ++i) {
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
int
MEHPForceBalance::getNumIntraChainSlipLinks() const
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
Eigen::VectorXd
MEHPForceBalance::evaluatePartialSpringVectors(const ForceBalanceNetwork& net,
                                               const Eigen::VectorXd& u,
                                               const bool use2D,
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
  if (use2D) {
    // partialDistances(Eigen::seq(2, net.nrOfPartialSprings, 3)) =
    //   Eigen::VectorXd::Zero(net.nrOfPartialSprings);
    for (Eigen::Index i = 2;
         i < static_cast<Eigen::Index>(3 * net.nrOfPartialSprings);
         i += 3) {
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
pylimer_tools::entities::Universe
MEHPForceBalance::getCrosslinkerVerse() const
{
  // convert nodes & springs back to a universe
  pylimer_tools::entities::Universe xlinkUniverse =
    pylimer_tools::entities::Universe(this->box);
  std::vector<long int> ids;
  std::vector<int> types = pylimer_tools::utils::initializeWithValue(
    this->initialConfig.nrOfNodes, crossLinkerType);
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  const std::vector<int> zeros =
    pylimer_tools::utils::initializeWithValue(this->initialConfig.nrOfNodes, 0);
  ids.reserve(this->initialConfig.nrOfNodes);
  x.reserve(this->initialConfig.nrOfNodes);
  y.reserve(this->initialConfig.nrOfNodes);
  z.reserve(this->initialConfig.nrOfNodes);
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->initialConfig.nrOfNodes);
       ++i) {
    x.push_back(this->initialConfig.coordinates[3 * i + 0] +
                this->currentDisplacements[3 * i + 0]);
    y.push_back(this->initialConfig.coordinates[3 * i + 1] +
                this->currentDisplacements[3 * i + 1]);
    z.push_back(this->initialConfig.coordinates[3 * i + 2] +
                this->currentDisplacements[3 * i + 2]);
    ids.push_back(this->initialConfig.oldAtomIds[i]);
    // override type, since the types may be different from
    // crossLinkerType if converted with dangling chains
    types[static_cast<size_t>(i)] = this->initialConfig.oldAtomTypes[i];
  }
  xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  bondFrom.reserve(this->initialConfig.nrOfSprings);
  bondTo.reserve(this->initialConfig.nrOfSprings);
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->initialConfig.nrOfSprings);
       ++i) {
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
void
MEHPForceBalance::addSlipLinksImpl(
  const std::vector<size_t>& strandIdx1,
  const std::vector<size_t>& strandIdx2,
  const Eigen::Ref<const Eigen::VectorXd>& x,
  const Eigen::Ref<const Eigen::VectorXd>& y,
  const Eigen::Ref<const Eigen::VectorXd>& z,
  const std::vector<double>& alpha1,
  const std::vector<double>& alpha2,
  std::vector<std::vector<size_t>> loops,
  std::vector<std::vector<size_t>> loopsOfSliplinks,
  const bool clampAlpha)
{
  const size_t additionalLen = strandIdx1.size();
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
  const size_t currentNrOfLinks = this->initialConfig.nrOfLinks;
  const size_t currentNrOfPartialSprings =
    this->initialConfig.nrOfPartialSprings;
  const Eigen::Index additionalLenE = static_cast<Eigen::Index>(additionalLen);
  if (additionalLenE != x.size() || additionalLenE != y.size() ||
      additionalLenE != z.size()) {
    throw std::invalid_argument("x, y and z must have the same dimensions");
  }
  if (additionalLen != strandIdx2.size() || additionalLen != alpha1.size() ||
      additionalLen != alpha2.size()) {
    throw std::invalid_argument(
      "Strand indices and alpha estimates must have the same length");
  }
  for (size_t i = 0; i < additionalLen; ++i) {
    INVALIDARG_EXP_IFN(strandIdx1[i] < this->initialConfig.nrOfSprings,
                       "Invalid spring index " + std::to_string(strandIdx1[i]) +
                         ", expected below " +
                         std::to_string(this->initialConfig.nrOfSprings) + ".");
    INVALIDARG_EXP_IFN(strandIdx2[i] < this->initialConfig.nrOfSprings,
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
  const Eigen::Index nrOfLinksE =
    static_cast<Eigen::Index>(this->initialConfig.nrOfLinks);
  const Eigen::Index newNrOfPartialSpringsE =
    static_cast<Eigen::Index>(currentNrOfPartialSprings + 2 * additionalLen);
  this->initialConfig.nrOfCrosslinkSwapsEndured.conservativeResize(
    static_cast<Eigen::Index>(this->initialConfig.nrOfLinks -
                              this->initialConfig.nrOfNodes));
  this->currentDisplacements.conservativeResize(3 * nrOfLinksE);
  this->currentSpringPartitionsVec.conservativeResize(newNrOfPartialSpringsE);
  this->initialConfig.springPartCoordinateIndexA.conservativeResize(
    3 * newNrOfPartialSpringsE);
  this->initialConfig.springPartCoordinateIndexB.conservativeResize(
    3 * newNrOfPartialSpringsE);
  this->initialConfig.springPartBoxOffset.conservativeResize(
    3 * newNrOfPartialSpringsE);
  this->initialConfig.springPartIndexA.conservativeResize(
    newNrOfPartialSpringsE);
  this->initialConfig.springPartIndexB.conservativeResize(
    newNrOfPartialSpringsE);
  this->initialConfig.linkIsSliplink.conservativeResize(nrOfLinksE);
  this->initialConfig.coordinates.conservativeResize(3 * nrOfLinksE);
  this->initialConfig.partialToFullSpringIndex.conservativeResize(
    newNrOfPartialSpringsE);
  this->initialConfig.partialSpringIsPartial.conservativeResize(
    newNrOfPartialSpringsE);
  // handle loops if appropriate
  const size_t previousNrOfLoops = this->initialConfig.loops.size();
  this->initialConfig.loops.reserve(previousNrOfLoops + loops.size());
  this->initialConfig.loops.insert(
    this->initialConfig.loops.end(), loops.begin(), loops.end());
  if (loopsOfSliplinks.size() > 0) {
    this->initialConfig.loopsOfSliplink.reserve(this->initialConfig.nrOfLinks -
                                                this->initialConfig.nrOfNodes);
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
    const Eigen::Index newLinkIdxE =
      static_cast<Eigen::Index>(currentNrOfLinks + i);
    const Eigen::Index iE = static_cast<Eigen::Index>(i);
    this->initialConfig.coordinates[3 * newLinkIdxE] = x[iE];
    this->initialConfig.coordinates[3 * newLinkIdxE + 1] = y[iE];
    this->initialConfig.coordinates[3 * newLinkIdxE + 2] = z[iE];
    this->initialConfig.linkIsSliplink[newLinkIdxE] = true;
    this->initialConfig.nrOfCrosslinkSwapsEndured
      [newLinkIdxE - static_cast<Eigen::Index>(this->initialConfig.nrOfNodes)] =
      0;
    std::vector<size_t> springIndices{ strandIdx1[i], strandIdx2[i] };
    std::vector<size_t> springIndicesOfLink =
      (strandIdx1[i] == strandIdx2[i]) ? std::vector<size_t>{ strandIdx1[i] }
                                       : springIndices;
    this->initialConfig.springIndicesOfLinks.push_back(springIndicesOfLink);
    // add to the springs
    int springIndexIndex = 0;
    for (const size_t springIndex : springIndices) {
      std::vector<size_t> springParticipants =
        this->initialConfig.linkIndicesOfSprings[springIndex];
      double alpha = (springIndexIndex == 0) ? alpha1[i] : alpha2[i];
      INVALIDARG_EXP_IFN(alpha >= 0.0 && alpha <= 1.0,
                         "alpha must be between 0 and 1, got " +
                           std::to_string(alpha) + ".");
      if (clampAlpha) {
        const Eigen::Index springIndexE =
          static_cast<Eigen::Index>(springIndex);
        alpha = std::clamp(
          alpha,
          1 / (this->initialConfig.springsContourLength[springIndexE]),
          1 - (1 / (this->initialConfig.springsContourLength[springIndexE])));
      }
      // detect the position in the spring
      std::vector<double> partitionsStrand;
      partitionsStrand.reserve(springParticipants.size() - 1);
      for (size_t j = 0; j < springParticipants.size() - 1; ++j) {
        partitionsStrand.push_back(
          this->currentSpringPartitionsVec[static_cast<Eigen::Index>(
            this->initialConfig.localToGlobalSpringIndex[springIndex][j])]);
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
      const size_t springPartner1 = springParticipants[targetIndexInSpring];
      const size_t springPartner2 = springParticipants[targetIndexInSpring + 1];
      size_t newNodeIdx = currentNrOfLinks + i;

      // update connectivity
      const size_t dividedPartialSpringIdx =
        this->initialConfig
          .localToGlobalSpringIndex[springIndex][targetIndexInSpring];
      const Eigen::Index dividedPartialSpringIdxE =
        static_cast<Eigen::Index>(dividedPartialSpringIdx);
      size_t newPartialSpringIdx =
        currentNrOfPartialSprings + partialSpringsAdded;
      const Eigen::Index newPartialSpringIdxE =
        static_cast<Eigen::Index>(newPartialSpringIdx);

      // for validation later
      Eigen::Vector3d distanceBefore =
        this->evaluatePartialSpringDistance(this->initialConfig,
                                            this->currentDisplacements,
                                            dividedPartialSpringIdx,
                                            this->is2D,
                                            false);

      this->initialConfig.partialSpringIsPartial[dividedPartialSpringIdxE] =
        true;
      this->initialConfig.partialSpringIsPartial[newPartialSpringIdxE] = true;

      this->initialConfig.localToGlobalSpringIndex[springIndex].insert(
        this->initialConfig.localToGlobalSpringIndex[springIndex].begin() +
          static_cast<std::ptrdiff_t>(targetIndexInSpring) + 1,
        newPartialSpringIdx);
      this->initialConfig.partialToFullSpringIndex[newPartialSpringIdxE] =
        static_cast<int>(springIndex);

      // adjust also the coordinates
      this->currentDisplacements.segment(
        3 * static_cast<Eigen::Index>(newNodeIdx), 3) = Eigen::Vector3d::Zero();
      assert(this->initialConfig.springPartIndexA[dividedPartialSpringIdxE] ==
             static_cast<int>(springPartner1));
      this->initialConfig.springPartIndexB[dividedPartialSpringIdxE] =
        static_cast<int>(newNodeIdx);
      for (Eigen::Index offset = 0; offset < 3; ++offset) {
        this->initialConfig
          .springPartCoordinateIndexB[3 * dividedPartialSpringIdxE + offset] =
          static_cast<int>(3 * newNodeIdx) + static_cast<int>(offset);
      }

      // add the new one
      this->initialConfig.springPartIndexA[newPartialSpringIdxE] =
        static_cast<int>(newNodeIdx);
      this->initialConfig.springPartIndexB[newPartialSpringIdxE] =
        static_cast<int>(springPartner2);
      for (Eigen::Index dir = 0; dir < 3; ++dir) {
        this->initialConfig
          .springPartCoordinateIndexA[3 * newPartialSpringIdxE + dir] =
          static_cast<int>(3 * newNodeIdx) + static_cast<int>(dir);
        this->initialConfig
          .springPartCoordinateIndexB[3 * newPartialSpringIdxE + dir] =
          static_cast<int>(3 * springPartner2) + static_cast<int>(dir);
      }

      // set box offsets
      this->initialConfig.springPartBoxOffset.segment(
        3 * newPartialSpringIdxE, 3) = Eigen::Vector3d::Zero();
      this->reAlignSlipLinkToImages(this->initialConfig,
                                    this->currentDisplacements,
                                    newNodeIdx,
                                    dividedPartialSpringIdx,
                                    newPartialSpringIdx);

      this->currentSpringPartitionsVec[newPartialSpringIdxE] =
        this->currentSpringPartitionsVec[dividedPartialSpringIdxE] - alpha;
      RUNTIME_EXP_IFN(
        APPROX_WITHIN(this->currentSpringPartitionsVec[newPartialSpringIdxE],
                      0.0,
                      1.0,
                      1e-9),
        "Spring partition must be between 0 and 1, got " +
          std::to_string(
            this->currentSpringPartitionsVec[newPartialSpringIdxE]) +
          ".");
      this->currentSpringPartitionsVec[dividedPartialSpringIdxE] = alpha;
      RUNTIME_EXP_IFN(
        APPROX_WITHIN(
          this->currentSpringPartitionsVec[dividedPartialSpringIdxE],
          0.0,
          1.0,
          1e-9),
        "Spring partition must be between 0 and 1, got " +
          std::to_string(
            this->currentSpringPartitionsVec[newPartialSpringIdxE]) +
          ".");

      this->initialConfig.linkIndicesOfSprings[springIndex].insert(
        this->initialConfig.linkIndicesOfSprings[springIndex].begin() +
          static_cast<std::ptrdiff_t>(targetIndexInSpring) +
          1, // + 1 to compensate for the first crosslink
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

      RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(distanceAfter,
                                                                distanceBefore),
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
    for (const size_t partialSpringIdx :
         this->initialConfig.localToGlobalSpringIndex[i]) {
      const Eigen::Index partialSpringIdxE =
        static_cast<Eigen::Index>(partialSpringIdx);
      containsPrimaryLoop =
        containsPrimaryLoop ||
        (this->initialConfig.springPartIndexA[partialSpringIdxE] ==
         this->initialConfig.springPartIndexB[partialSpringIdxE]);
    }
    Eigen::Vector3d vecBefore =
      springVectorsBefore.segment(3 * static_cast<Eigen::Index>(i), 3);
    Eigen::Vector3d vecAfter =
      springVectorsAfter.segment(3 * static_cast<Eigen::Index>(i), 3);
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
double
MEHPForceBalance::getAverageSpringLength() const
{
  double r2 = 0.0;
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->initialConfig.nrOfSprings);
       i++) {
    double r2local = 0.0;
    for (int j = 0; j < 3; ++j) {
      r2local += this->currentSpringVectors[i * 3 + j] *
                 this->currentSpringVectors[i * 3 + j];
    }
    r2 += sqrt(r2local);
  }
  return r2 / static_cast<double>(this->initialConfig.nrOfSprings);
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
double
MEHPForceBalance::getDenominatorOfPartialSpring(
  const ForceBalanceNetwork& net,
  const Eigen::VectorXd& springPartitions,
  const size_t partialSpringIdx,
  const double oneOverSpringPartitionUpperLimit)
{
  const Eigen::Index partialSpringIdxE =
    static_cast<Eigen::Index>(partialSpringIdx);
  const double N =
    net.springsContourLength[net.partialToFullSpringIndex[partialSpringIdxE]];
  const double fraction = springPartitions[partialSpringIdxE];

  double denominator = 1. / (fraction * N);
  if (oneOverSpringPartitionUpperLimit > 0. || !std::isfinite(denominator)) {
    denominator = CLAMP_ONE_OVER_SPRINGPARTITION(
      net.partialSpringIsPartial[partialSpringIdxE],
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
Eigen::Matrix3d
MEHPForceBalance::evaluateStressTensorForLinks(
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

  const double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXi debugNrSpringsVisited =
    Eigen::VectorXi::Zero(static_cast<Eigen::Index>(net.nrOfPartialSprings));

  for (const size_t linkIdx : linkIndices) {
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

  const double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXi debugNrSpringsVisited =
    Eigen::VectorXi::Zero(static_cast<Eigen::Index>(net.nrOfPartialSprings));

  const size_t nrOfLinksToInspect = xlinksOnly ? net.nrOfNodes : net.nrOfLinks;
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
                      std::to_string(stressOnLink.squaredNorm()) + " at link " +
                      std::to_string(linkIdx) + "!");
    stress += stressOnLink;
  }

  std::array<std::array<double, 3>, 3> stressA;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      stressA[i][j] = halfOverVolume * stress(static_cast<Eigen::Index>(i),
                                              static_cast<Eigen::Index>(j));
    }
  }

  if (!xlinksOnly) {
    RUNTIME_EXP_IFN(
      static_cast<Eigen::Index>(debugNrSpringsVisited.sum()) ==
        static_cast<Eigen::Index>(2 * net.nrOfPartialSprings),
      "Every spring must be visited twice, got min " +
        std::to_string(debugNrSpringsVisited.minCoeff()) + " and max " +
        std::to_string(debugNrSpringsVisited.maxCoeff()) + ". Sum is " +
        std::to_string(debugNrSpringsVisited.sum()) + " instead of " +
        std::to_string(2 * net.nrOfPartialSprings) + ".");
    RUNTIME_EXP_IFN((debugNrSpringsVisited.array() == 2).all(),
                    "Every spring must be visited twice, got min " +
                      std::to_string(debugNrSpringsVisited.minCoeff()) +
                      " and max " +
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
MEHPForceBalance::evaluateStressTensor(
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

  const double oneOverVolume = 1. / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXd displacedCoords = net.coordinates + u;
  Eigen::VectorXd relevantPartialDistancesA =
    (displacedCoords(net.springPartCoordinateIndexB) -
     displacedCoords(net.springPartCoordinateIndexA)) +
    net.springPartBoxOffset;

  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(relevantPartialDistancesA);
  }

  if (this->is2D) {
    for (Eigen::Index i = 2; i < relevantPartialDistancesA.size(); i += 3) {
      relevantPartialDistancesA[i] = 0.;
    }
  }

  for (Eigen::Index partialSpringIdx = 0;
       partialSpringIdx < static_cast<Eigen::Index>(net.nrOfPartialSprings);
       ++partialSpringIdx) {
    Eigen::Vector3d distance =
      relevantPartialDistancesA.segment(3 * partialSpringIdx, 3);
    const double contourLengthFraction = springPartitions[partialSpringIdx];
    const double N =
      net.springsContourLength[net.partialToFullSpringIndex[partialSpringIdx]];
    const double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
      net.partialSpringIsPartial[partialSpringIdx],
      1.0 / (N * contourLengthFraction),
      N,
      oneOverSpringPartitionUpperLimit);

    /* spring contribution to the overall stress tensor */
    for (Eigen::Index j = 0; j < 3; j++) {
      for (Eigen::Index k = 0; k < 3; k++) {
        const double contribution = distance[j] * distance[k] * this->kappa *
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
        stress[static_cast<size_t>(j)][static_cast<size_t>(k)] += contribution;
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
                        " from denominator " + std::to_string(oneOverVolume) +
                        ".");
    }
  }

  return stress;
}

Eigen::Matrix3d
MEHPForceBalance::getStressTensor(
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
    convertedRes.row(static_cast<Eigen::Index>(i)) =
      Eigen::Vector3d::Map(res[i].data(), 3);
  }
  return convertedRes;
}

Eigen::Matrix3d
MEHPForceBalance::getStressTensorLinkBased(
  const double oneOverSpringPartitionUpperLimit,
  const bool xlinksOnly) const
{
  const std::array<std::array<double, 3>, 3> res =
    this->evaluateStressTensorLinkBased(this->initialConfig,
                                        this->currentDisplacements,
                                        this->currentSpringPartitionsVec,
                                        oneOverSpringPartitionUpperLimit,
                                        xlinksOnly);
  Eigen::Matrix3d convertedRes = Eigen::Matrix3d::Zero();
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      convertedRes(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
        res[i][j];
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
MEHPForceBalance::getEffectiveFunctionalityOfAtoms(const double tolerance) const
{
  std::unordered_map<long int, int> results;
  results.reserve(this->initialConfig.nrOfNodes);

  Eigen::VectorXi nrOfActiveSpringsConnected =
    this->getNrOfActiveSpringsConnected(tolerance);
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->initialConfig.nrOfNodes);
       i++) {
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
std::vector<long int>
MEHPForceBalance::getIndicesOfActiveNodes(
  const ForceBalanceNetwork* net,
  const Eigen::VectorXd& u,
  const Eigen::VectorXd& springPartitions,
  const double tolerance) const
{
  std::vector<long int> results;
  results.reserve(net->nrOfNodes);

  // find all active springs
  Eigen::ArrayXb springIsActive =
    this->findActiveSprings(net, u, springPartitions, tolerance);

  for (size_t i = 0; i < net->nrOfNodes; i++) {
    if (net->oldAtomTypes[static_cast<Eigen::Index>(i)] !=
        this->entanglementType) {
      std::vector<size_t> springIndices =
        this->getInvolvedFullSpringIndices(*net, i);
      for (const size_t springIndex : springIndices) {
        if (springIsActive[static_cast<Eigen::Index>(springIndex)]) {
          results.push_back(static_cast<long int>(i));
          break;
        }
      }
    }
  }

  return results;
};

/**
 * @brief Get the atom ids of the active crosslinks (not entanglement
 * beads/links)
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return std::vector<long int> the atom ids
 */
std::vector<long int>
MEHPForceBalance::getIdsOfActiveNodes(const double tolerance) const
{
  std::vector<long int> results;
  // find all active springs
  const std::vector<long int> activeNodes =
    this->getIndicesOfActiveNodes(&this->initialConfig,
                                  this->currentDisplacements,
                                  this->currentSpringPartitionsVec,
                                  tolerance);

  results.reserve(activeNodes.size());

  for (const long int nodeIdx : activeNodes) {
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
Eigen::VectorXi
MEHPForceBalance::getNrOfActiveSpringsConnected(const double tolerance) const
{
  Eigen::VectorXi nrOfActiveSpringsConnected = Eigen::VectorXi::Zero(
    static_cast<Eigen::Index>(this->initialConfig.nrOfNodes));
  Eigen::ArrayXb springIsActive = this->findActiveSprings(tolerance);
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->initialConfig.nrOfSprings);
       i++) {
    if (springIsActive[i]) { /* active spring */
      const int a = this->initialConfig.springIndexA[i];
      const int b = this->initialConfig.springIndexB[i];
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
Eigen::VectorXi
MEHPForceBalance::getNrOfActivePartialSpringsConnected(
  const double tolerance) const
{
  Eigen::VectorXi nrOfActivePartialSpringsConnected = Eigen::VectorXi::Zero(
    static_cast<Eigen::Index>(this->initialConfig.nrOfNodes));
  Eigen::ArrayXb partialSpringIsActive =
    this->findActivePartialSprings(tolerance);
  // translate this to the nodes
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->initialConfig.nrOfPartialSprings);
       ++i) {
    if (partialSpringIsActive[i]) {
      /* active spring */
      // size_t a =
      //   this->initialConfig
      //     .springIndexA[this->initialConfig.partialToFullSpringIndex[i]];
      // size_t b =
      //   this->initialConfig
      //     .springIndexB[this->initialConfig.partialToFullSpringIndex[i]];
      const Eigen::Index a = this->initialConfig.springPartIndexA[i];
      const Eigen::Index b = this->initialConfig.springPartIndexB[i];
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
double
MEHPForceBalance::getGammaFactor(double b02,
                                 const int nrOfChains,
                                 double oneOverSpringPartitionUpperLimit) const
{
  if (b02 < 0) {
    b02 = this->defaultBondLength * this->defaultBondLength;
  }

  if (this->getNrOfSprings() == 0) {
    return 0.;
  }

  const Eigen::VectorXd gammaFactors = this->getGammaFactors(b02);

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
Eigen::VectorXd
MEHPForceBalance::getGammaFactors(
  const double b02,
  const double oneOverSpringPartitionUpperLimit) const
{
  Eigen::VectorXd springVectors = this->evaluatePartialSpringVectors(
    this->initialConfig, this->currentDisplacements);
  RUNTIME_EXP_IFN(this->currentSpringPartitionsVec.size() * 3 ==
                    springVectors.size(),
                  "Unexpected dimensions in springVectors.");

  Eigen::VectorXd gammaFactors(springVectors.size() / 3);
  const double commonDenominator = 1. / b02;
  for (Eigen::Index i = 0; i < springVectors.size() / 3; ++i) {
    const double N =
      this->initialConfig
        .springsContourLength[this->initialConfig.partialToFullSpringIndex[i]];
    const double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
      this->initialConfig.partialSpringIsPartial[i],
      1.0 / (N * this->currentSpringPartitionsVec(i)),
      N,
      oneOverSpringPartitionUpperLimit);
    gammaFactors[i] = springVectors.segment(3 * i, 3).squaredNorm() *
                      commonDenominator * oneOverContourLengthFraction;
    RUNTIME_EXP_IFN(
      std::isfinite(gammaFactors[i]),
      "Non-finite gamma factor for partial spring " + std::to_string(i) +
        ", computed from N = " + std::to_string(N) + ", b02 = " +
        std::to_string(b02) + ", oneOverSpringPartitionUpperLimit = " +
        std::to_string(oneOverSpringPartitionUpperLimit) +
        " and squared distance = " +
        std::to_string(springVectors.segment(3 * i, 3).squaredNorm()) + ".");
  }
  return gammaFactors;
}

/**
 * @brief Get the per-(partial)-spring gamma factors
 *
 * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
 * computed as phantom = N<b^2>.
 * @param dir the direction (0=x, 1=y, 2=z)
 * @return Eigen::VectorXd
 */
Eigen::VectorXd
MEHPForceBalance::getGammaFactorsInDir(
  const double b02,
  const int dir,
  const double oneOverSpringPartitionUpperLimit) const
{
  INVALIDARG_EXP_IFN(dir >= 0 && dir <= 2, "Invalid direction.");
  Eigen::VectorXd springVectors = this->evaluatePartialSpringVectors(
    this->initialConfig, this->currentDisplacements);
  RUNTIME_EXP_IFN(this->currentSpringPartitionsVec.size() * 3 ==
                    springVectors.size(),
                  "Unexpected dimensions in springVectors.");

  Eigen::VectorXd gammaFactors(springVectors.size() / 3);
  const double commonDenominator = 1. / b02;
  for (Eigen::Index i = 0; i < springVectors.size() / 3; ++i) {
    const double N =
      this->initialConfig
        .springsContourLength[this->initialConfig.partialToFullSpringIndex[i]];
    const double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
      this->initialConfig.partialSpringIsPartial[i],
      1.0 / (N * this->currentSpringPartitionsVec(i)),
      N,
      oneOverSpringPartitionUpperLimit);
    gammaFactors[i] =
      SQUARE(springVectors[3 * i + static_cast<Eigen::Index>(dir)]) *
      commonDenominator * oneOverContourLengthFraction;
    RUNTIME_EXP_IFN(
      std::isfinite(gammaFactors[i]),
      "Non-finite gamma factor for partial spring " + std::to_string(i) +
        ", computed from N = " + std::to_string(N) +
        ", oneOverSpringPartitionUpperLimit = " +
        std::to_string(oneOverSpringPartitionUpperLimit) +
        " and squared distance = " +
        std::to_string(
          SQUARE(springVectors[3 * i + static_cast<Eigen::Index>(dir)])) +
        " in dir " + std::to_string(dir) + ".");
  }
  return gammaFactors;
}

/**
 * @brief Get the Weighted Partial Spring Length
 *
 * @return double
 */
double
MEHPForceBalance::getWeightedPartialSpringLength(
  const ForceBalanceNetwork& net,
  const Eigen::VectorXd& u,
  const Eigen::VectorXd& springPartitions,
  const size_t partialSpringIdx,
  const double oneOverSpringPartitionUpperLimit) const
{
  const Eigen::Index partialSpringIdxE =
    static_cast<Eigen::Index>(partialSpringIdx);
  const double N =
    net.springsContourLength[net.partialToFullSpringIndex[partialSpringIdxE]];
  const double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
    net.partialSpringIsPartial[partialSpringIdxE],
    1.0 / (N * springPartitions(partialSpringIdxE)),
    N,
    oneOverSpringPartitionUpperLimit);
  return this->evaluatePartialSpringDistance(net, u, partialSpringIdx).norm() *
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
bool
MEHPForceBalance::ConvertNetwork(ForceBalanceNetwork& net,
                                 const int xlinkType,
                                 const bool remove2functionalCrosslinkers,
                                 const bool removeDanglingChains)
{
  const size_t nrOfAtoms = this->universe.getNrOfAtoms();
  const size_t nrOfBonds = this->universe.getNrOfBonds();

  const std::vector<int> vertexDegrees = this->universe.getVertexDegrees();
  std::vector<int> vertexTypes = this->universe.getPropertyValues<int>("type");
  const std::vector<long int> vertexAtomIds =
    this->universe.getPropertyValues<long int>("id");

  if (remove2functionalCrosslinkers) {
    for (size_t vIdx = 0; vIdx < nrOfAtoms; ++vIdx) {
      // change type of crosslinkers with a degree <= 2 to "normal",
      // non-crosslink beads
      if (vertexTypes[vIdx] == xlinkType && vertexDegrees[vIdx] <= 2) {
        vertexTypes[vIdx] = xlinkType - 1;
        this->universe.setPropertyValue(
          static_cast<igraph_int_t>(vIdx), "type", xlinkType - 1);
      }
    }
  }

  const igraph_t& graph = this->universe.getGraph();
  std::vector<long int> edgeTypes(nrOfBonds, -1);
  if (nrOfBonds > 0 &&
      igraph_cattribute_has_attr(&graph, IGRAPH_ATTRIBUTE_EDGE, "type")) {
    igraph_vector_t typesVec;
    REQUIRE_IGRAPH_SUCCESS(igraph_vector_init(&typesVec, 0));
    if (igraph_cattribute_EANV(
          &graph, "type", igraph_ess_all(IGRAPH_EDGEORDER_ID), &typesVec)) {
      igraph_vector_destroy(&typesVec);
      throw std::runtime_error("Failed to fetch edge type attribute");
    }
    assert(static_cast<size_t>(igraph_vector_size(&typesVec)) == nrOfBonds);
    for (size_t eIdx = 0; eIdx < nrOfBonds; ++eIdx) {
      edgeTypes[eIdx] = static_cast<long int>(std::llround(
        igraph_vector_get(&typesVec, static_cast<igraph_int_t>(eIdx))));
    }
    igraph_vector_destroy(&typesVec);
  }

  const Eigen::VectorXd vertexCoordinates =
    nrOfAtoms > 0 ? this->universe.getUnwrappedVertexCoordinates(this->box)
                  : Eigen::VectorXd::Zero(0);

  std::vector<bool> vertexIsJunction(nrOfAtoms, false);
  for (size_t vIdx = 0; vIdx < nrOfAtoms; ++vIdx) {
    if (vertexDegrees[vIdx] > 2 || vertexTypes[vIdx] == xlinkType) {
      vertexIsJunction[vIdx] = true;
    }
  }

  std::vector<bool> vertexVisited(nrOfAtoms, false);
  std::vector<bool> edgeVisited(nrOfBonds, false);
  std::vector<long int> vertexIdToLinkIdx(nrOfAtoms, -1);
  size_t currentLinkIdx = 0;

  struct DiscoveredStrand
  {
    size_t fromV;
    size_t toV;
    double contourLength;
    int springType;
    Eigen::Vector3d expectedDistance;
    std::vector<size_t> atomVertices;
  };
  std::vector<DiscoveredStrand> discoveredStrands;
  discoveredStrands.reserve(nrOfBonds / 2 + 1);

  igraph_adjlist_t adjacencyList;
  REQUIRE_IGRAPH_SUCCESS(igraph_adjlist_init(
    &graph, &adjacencyList, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
  igraph_inclist_t incidenceList;
  if (igraph_inclist_init(
        &graph, &incidenceList, IGRAPH_ALL, IGRAPH_LOOPS_TWICE)) {
    igraph_adjlist_destroy(&adjacencyList);
    throw std::runtime_error("Failed to initialize incidence list.");
  }

  auto addDiscoveredStrand = [&](const size_t endA,
                                 const size_t endB,
                                 const size_t nBonds,
                                 const long int bondTypeSum,
                                 Eigen::Vector3d expectedDistance,
                                 std::vector<size_t> atomVertices) {
    size_t fromV = endA;
    size_t toV = endB;
    if (vertexAtomIds[fromV] > vertexAtomIds[toV]) {
      std::swap(fromV, toV);
      expectedDistance = -expectedDistance;
    }

    if (vertexIdToLinkIdx[fromV] == -1) {
      vertexIdToLinkIdx[fromV] = static_cast<long int>(currentLinkIdx);
      currentLinkIdx += 1;
    }
    if (vertexIdToLinkIdx[toV] == -1) {
      vertexIdToLinkIdx[toV] = static_cast<long int>(currentLinkIdx);
      currentLinkIdx += 1;
    }

    discoveredStrands.push_back(
      DiscoveredStrand{ fromV,
                        toV,
                        static_cast<double>(nBonds),
                        static_cast<int>(static_cast<double>(bondTypeSum) /
                                         static_cast<double>(nBonds)),
                        expectedDistance,
                        std::move(atomVertices) });
  };

  // Helper to walk from a non-junction vertex along one incident edge until
  // reaching a junction, a degree-1 end, or returning to startBeadV.
  auto walkFromNonJunction = [&](const size_t startBeadV,
                                 const igraph_int_t dirIdx,
                                 size_t& outEndV,
                                 size_t& outNBonds,
                                 long int& outBondTypeSum,
                                 Eigen::Vector3d& outVec,
                                 std::vector<size_t>& outVertices) -> bool {
    const igraph_vector_int_t* startNeighbors =
      igraph_adjlist_get(&adjacencyList, static_cast<igraph_int_t>(startBeadV));
    const igraph_vector_int_t* startIncidentEdges =
      igraph_inclist_get(&incidenceList, static_cast<igraph_int_t>(startBeadV));
    const size_t firstEdge =
      static_cast<size_t>(igraph_vector_int_get(startIncidentEdges, dirIdx));
    if (edgeVisited[firstEdge]) {
      return false;
    }
    edgeVisited[firstEdge] = true;

    size_t currV =
      static_cast<size_t>(igraph_vector_int_get(startNeighbors, dirIdx));
    outNBonds = 1;
    outBondTypeSum = edgeTypes[firstEdge];
    outVec =
      vertexCoordinates.segment<3>(3 * static_cast<Eigen::Index>(currV)) -
      vertexCoordinates.segment<3>(3 * static_cast<Eigen::Index>(startBeadV));
    this->box.handlePBC(outVec);
    outVertices.push_back(currV);

    if (currV == startBeadV) {
      outEndV = currV;
      return true;
    }

    while (!vertexIsJunction[currV] && vertexDegrees[currV] == 2) {
      vertexVisited[currV] = true;
      const igraph_vector_int_t* neighbors =
        igraph_adjlist_get(&adjacencyList, static_cast<igraph_int_t>(currV));
      const igraph_vector_int_t* incidentEdges =
        igraph_inclist_get(&incidenceList, static_cast<igraph_int_t>(currV));
      assert(igraph_vector_int_size(neighbors) == 2);
      assert(igraph_vector_int_size(incidentEdges) == 2);

      const size_t e0 =
        static_cast<size_t>(igraph_vector_int_get(incidentEdges, 0));
      const igraph_int_t nextChoice = edgeVisited[e0] ? 1 : 0;
      const size_t nextEdge =
        static_cast<size_t>(igraph_vector_int_get(incidentEdges, nextChoice));
      assert(!edgeVisited[nextEdge]);
      edgeVisited[nextEdge] = true;

      const size_t nextV =
        static_cast<size_t>(igraph_vector_int_get(neighbors, nextChoice));
      Eigen::Vector3d stepVec =
        vertexCoordinates.segment<3>(3 * static_cast<Eigen::Index>(nextV)) -
        vertexCoordinates.segment<3>(3 * static_cast<Eigen::Index>(currV));
      this->box.handlePBC(stepVec);
      outVec += stepVec;
      outNBonds += 1;
      outBondTypeSum += edgeTypes[nextEdge];
      currV = nextV;
      outVertices.push_back(currV);

      if (currV == startBeadV) {
        break;
      }
    }

    if (!vertexIsJunction[currV]) {
      vertexVisited[currV] = true;
    }
    outEndV = currV;
    return true;
  };

  // Pass 1: Discover all strands containing >= 1 non-junction bead in order
  // of their minimum non-junction vertex index (matching igraph_decompose).
  for (size_t vIdx = 0; vIdx < nrOfAtoms; ++vIdx) {
    if (vertexIsJunction[vIdx] || vertexVisited[vIdx]) {
      continue;
    }
    vertexVisited[vIdx] = true;
    if (vertexDegrees[vIdx] == 0) {
      continue;
    }

    std::vector<size_t> atomVertices;
    atomVertices.push_back(vIdx);

    if (vertexDegrees[vIdx] == 1) {
      size_t endV = vIdx;
      size_t nBonds = 0;
      long int bondTypeSum = 0;
      Eigen::Vector3d vec = Eigen::Vector3d::Zero();
      const bool walked = walkFromNonJunction(
        vIdx, 0, endV, nBonds, bondTypeSum, vec, atomVertices);
      assert(walked);
      (void)walked;
      if (vertexIsJunction[endV] && !removeDanglingChains) {
        addDiscoveredStrand(
          vIdx, endV, nBonds, bondTypeSum, vec, std::move(atomVertices));
      }
      continue;
    }

    assert(vertexDegrees[vIdx] == 2);
    size_t end0 = vIdx;
    size_t nBonds0 = 0;
    long int bondTypeSum0 = 0;
    Eigen::Vector3d vec0 = Eigen::Vector3d::Zero();
    const bool walked0 = walkFromNonJunction(
      vIdx, 0, end0, nBonds0, bondTypeSum0, vec0, atomVertices);
    if (!walked0 || end0 == vIdx) {
      // Free primary loop of non-junction beads: skip
      continue;
    }

    size_t end1 = vIdx;
    size_t nBonds1 = 0;
    long int bondTypeSum1 = 0;
    Eigen::Vector3d vec1 = Eigen::Vector3d::Zero();
    const bool walked1 = walkFromNonJunction(
      vIdx, 1, end1, nBonds1, bondTypeSum1, vec1, atomVertices);
    assert(walked1);
    (void)walked1;

    const size_t totalBonds = nBonds0 + nBonds1;
    const long int totalBondTypeSum = bondTypeSum0 + bondTypeSum1;
    const Eigen::Vector3d expectedDistance = vec1 - vec0;

    bool useStrand = false;
    if (vertexIsJunction[end0] && vertexIsJunction[end1]) {
      if (end0 == end1) {
        // Primary loop: only use if the junction is of xlinkType
        if (vertexTypes[end0] == xlinkType) {
          useStrand = true;
        }
        // Avoid duplicate entry for the single junction in atomVertices
        atomVertices.pop_back();
      } else {
        // NETWORK_STRAND
        useStrand = true;
      }
    } else if ((vertexIsJunction[end0] || vertexIsJunction[end1]) &&
               !removeDanglingChains) {
      // DANGLING_CHAIN
      useStrand = true;
    }

    if (useStrand) {
      addDiscoveredStrand(end0,
                          end1,
                          totalBonds,
                          totalBondTypeSum,
                          expectedDistance,
                          std::move(atomVertices));
    }
  }

  // Pass 2: Direct bonds between two junctions and self-loops on junctions.
  for (size_t startV = 0; startV < nrOfAtoms; ++startV) {
    if (!vertexIsJunction[startV]) {
      continue;
    }
    const igraph_vector_int_t* startNeighbors =
      igraph_adjlist_get(&adjacencyList, static_cast<igraph_int_t>(startV));
    const igraph_vector_int_t* startIncidentEdges =
      igraph_inclist_get(&incidenceList, static_cast<igraph_int_t>(startV));
    const igraph_int_t startDeg = igraph_vector_int_size(startNeighbors);
    assert(startDeg == igraph_vector_int_size(startIncidentEdges));

    std::vector<size_t> selfLoopEdges;
    for (igraph_int_t i = 0; i < startDeg; ++i) {
      const size_t edgeIdx =
        static_cast<size_t>(igraph_vector_int_get(startIncidentEdges, i));
      if (edgeVisited[edgeIdx]) {
        continue;
      }
      edgeVisited[edgeIdx] = true;

      const size_t endV =
        static_cast<size_t>(igraph_vector_int_get(startNeighbors, i));
      assert(vertexIsJunction[endV]);
      if (endV == startV) {
        selfLoopEdges.push_back(edgeIdx);
        continue;
      }

      bool useStrand = false;
      if (vertexDegrees[startV] > 1 && vertexDegrees[endV] > 1) {
        useStrand = true;
      } else if (!removeDanglingChains &&
                 (vertexDegrees[startV] > 1 || vertexDegrees[endV] > 1)) {
        useStrand = true;
      }
      if (!useStrand) {
        continue;
      }

      Eigen::Vector3d expectedDistance =
        vertexCoordinates.segment<3>(3 * static_cast<Eigen::Index>(endV)) -
        vertexCoordinates.segment<3>(3 * static_cast<Eigen::Index>(startV));
      this->box.handlePBC(expectedDistance);
      addDiscoveredStrand(startV,
                          endV,
                          1,
                          edgeTypes[edgeIdx],
                          expectedDistance,
                          std::vector<size_t>{ startV, endV });
    }

    for (const size_t edgeIdx : selfLoopEdges) {
      if (vertexTypes[startV] == xlinkType) {
        addDiscoveredStrand(startV,
                            startV,
                            1,
                            edgeTypes[edgeIdx],
                            Eigen::Vector3d::Zero(),
                            std::vector<size_t>{ startV });
      }
    }
  }

  igraph_adjlist_destroy(&adjacencyList);
  igraph_inclist_destroy(&incidenceList);

  const size_t nrOfSprings = discoveredStrands.size();
  const size_t nrOfXlinks = currentLinkIdx;

  // crossLinkerUniverse.simplify();
  const pylimer_tools::entities::Box universeBox = this->box;
  net.L[0] = universeBox.getLx();
  net.L[1] = universeBox.getLy();
  net.L[2] = universeBox.getLz();
  net.boxHalfs[0] = 0.5 * net.L[0];
  net.boxHalfs[1] = 0.5 * net.L[1];
  net.boxHalfs[2] = 0.5 * net.L[2];
  net.nrOfNodes = nrOfXlinks;
  net.nrOfLinks = nrOfXlinks;
  net.nrOfSprings = nrOfSprings;
  net.nrOfPartialSprings = nrOfSprings;
  net.nrOfSpringsWithPartition = 0;
  const Eigen::Index nrOfLinksE = static_cast<Eigen::Index>(net.nrOfLinks);
  const Eigen::Index nrOfSpringsE = static_cast<Eigen::Index>(net.nrOfSprings);
  const Eigen::Index nrOfPartialSpringsE =
    static_cast<Eigen::Index>(net.nrOfPartialSprings);
  net.coordinates = Eigen::VectorXd::Zero(3 * nrOfLinksE);
  net.nrOfCrosslinkSwapsEndured = Eigen::ArrayXi::Zero(0);
  net.oldAtomIds = Eigen::ArrayXi::Zero(nrOfLinksE);
  net.oldAtomTypes = Eigen::ArrayXi::Zero(nrOfLinksE);
  net.linkIsSliplink = Eigen::ArrayXb::Constant(nrOfLinksE, false);
  net.springIndicesOfLinks.reserve(net.nrOfLinks);
  for (size_t i = 0; i < net.nrOfLinks; ++i) {
    net.springIndicesOfLinks.push_back(std::vector<size_t>());
  }
  net.linkIndicesOfSprings.reserve(net.nrOfSprings);
  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    net.linkIndicesOfSprings.push_back(std::vector<size_t>());
  }
  net.partialToFullSpringIndex = Eigen::ArrayXi(nrOfPartialSpringsE);
  this->currentSpringPartitionsVec = Eigen::VectorXd::Ones(nrOfSpringsE);
  net.springIndexA = Eigen::ArrayXi::Zero(nrOfSpringsE);
  net.springIndexB = Eigen::ArrayXi::Zero(nrOfSpringsE);
  net.springCoordinateIndexA = Eigen::ArrayXi::Zero(3 * nrOfSpringsE);
  net.springCoordinateIndexB = Eigen::ArrayXi::Zero(3 * nrOfSpringsE);
  net.springPartBoxOffset = Eigen::VectorXd::Zero(3 * nrOfSpringsE);
  net.springIsActive = Eigen::ArrayXb::Constant(nrOfSpringsE, false);
  net.springsContourLength = Eigen::VectorXd::Zero(nrOfSpringsE);
  net.springsType = Eigen::ArrayXi::Zero(nrOfSpringsE);
  net.oldAtomIdToSpringIndex.reserve(nrOfAtoms);
  net.springToMoleculeIds.reserve(nrOfSprings);
  net.partialSpringIsPartial = Eigen::ArrayXb::Constant(nrOfSpringsE, false);

  // convert (crosslinker-)beads
  for (size_t i = 0; i < vertexIdToLinkIdx.size(); ++i) {
    if (vertexIdToLinkIdx[i] != -1) {
      const size_t linkIdx = static_cast<size_t>(vertexIdToLinkIdx[i]);
      const Eigen::Index linkIdxE = static_cast<Eigen::Index>(linkIdx);
      net.oldAtomIds[linkIdxE] = static_cast<int>(vertexAtomIds[i]);
      net.oldAtomTypes[linkIdxE] = vertexTypes[i];
      Eigen::Vector3d coords =
        vertexCoordinates.segment<3>(3 * static_cast<Eigen::Index>(i));
      this->box.handlePBC(coords);
      net.coordinates.segment(3 * linkIdxE, 3) = coords;
    }
  }

  // convert springs
  size_t spring_idx = 0;
  for (size_t i = 0; i < nrOfSprings; ++i) {
    const DiscoveredStrand& strand = discoveredStrands[i];
    const Eigen::Index springIdxE = static_cast<Eigen::Index>(spring_idx);

    net.springsContourLength[springIdxE] = strand.contourLength;
    net.springsType[springIdxE] = strand.springType;

    const long int nodeIdxFrom = vertexIdToLinkIdx[strand.fromV];
    const long int nodeIdxTo = vertexIdToLinkIdx[strand.toV];
    assert(nodeIdxFrom >= 0 && nodeIdxTo >= 0);

    net.springToMoleculeIds.push_back(spring_idx);
    for (const size_t vIdx : strand.atomVertices) {
      net.oldAtomIdToSpringIndex[static_cast<size_t>(vertexAtomIds[vIdx])] =
        spring_idx;
    }

    pylimer_tools::utils::addIfNotContained(
      net.springIndicesOfLinks[static_cast<size_t>(nodeIdxFrom)], spring_idx);
    if (nodeIdxFrom != nodeIdxTo) {
      pylimer_tools::utils::addIfNotContained(
        net.springIndicesOfLinks[static_cast<size_t>(nodeIdxTo)], spring_idx);
    }

    net.linkIndicesOfSprings[spring_idx].push_back(
      static_cast<size_t>(nodeIdxFrom));
    net.linkIndicesOfSprings[spring_idx].push_back(
      static_cast<size_t>(nodeIdxTo));

    net.springIndexA[springIdxE] = static_cast<int>(nodeIdxFrom);
    net.springIndexB[springIdxE] = static_cast<int>(nodeIdxTo);
    for (Eigen::Index j = 0; j < 3; j++) {
      net.springCoordinateIndexA[3 * springIdxE + j] =
        static_cast<int>(nodeIdxFrom * 3 + j);
      net.springCoordinateIndexB[3 * springIdxE + j] =
        static_cast<int>(nodeIdxTo * 3 + j);
    }

    std::vector<size_t> zeroMap;
    zeroMap.push_back(spring_idx);
    net.localToGlobalSpringIndex.push_back(zeroMap);
    net.partialToFullSpringIndex[springIdxE] = static_cast<int>(spring_idx);

    const Eigen::Vector3d actualDistance =
      net.coordinates.segment(3 * static_cast<Eigen::Index>(nodeIdxTo), 3) -
      net.coordinates.segment(3 * static_cast<Eigen::Index>(nodeIdxFrom), 3);
    net.springPartBoxOffset.segment(3 * springIdxE, 3) =
      strand.expectedDistance - actualDistance;
    assert(this->box.isValidOffset(strand.expectedDistance - actualDistance));

    spring_idx += 1;
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

bool
MEHPForceBalance::validateNetwork(const ForceBalanceNetwork& net,
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
  RUNTIME_EXP_IFN(APPROX_EQUAL(net.L[0] * net.L[1] * net.L[2], net.vol, 1e-6),
                  "Inconsistent box volume and dimensions");
  RUNTIME_EXP_IFN(net.coordinates.size() ==
                    static_cast<Eigen::Index>(net.nrOfLinks * 3),
                  "Invalid size of coordinates");
  RUNTIME_EXP_IFN(u.size() == static_cast<Eigen::Index>(net.nrOfLinks * 3),
                  "Invalid size of displacements");
  RUNTIME_EXP_IFN(u.size() == net.coordinates.size(),
                  "Invalid size of displacements or coordinates");
  RUNTIME_EXP_IFN(net.localToGlobalSpringIndex.size() == net.nrOfSprings,
                  "Invalid size of connectivity map, got " +
                    std::to_string(net.localToGlobalSpringIndex.size()) +
                    " for " + std::to_string(net.nrOfSprings) + " springs.");
  RUNTIME_EXP_IFN(net.springsContourLength.size() ==
                    static_cast<Eigen::Index>(net.nrOfSprings),
                  "Invalid size of contour lengths, got " +
                    std::to_string(net.springsContourLength.size()) + " for " +
                    std::to_string(net.nrOfSprings) + " springs.");
  RUNTIME_EXP_IFN(net.springsType.size() ==
                    static_cast<Eigen::Index>(net.nrOfSprings),
                  "Invalid size of springs types, got " +
                    std::to_string(net.springsType.size()) + " for " +
                    std::to_string(net.nrOfSprings) + " springs.");
  RUNTIME_EXP_IFN(net.springIndicesOfLinks.size() == net.nrOfLinks,
                  "Invalid size of spring indices of links, got " +
                    std::to_string(net.linkIndicesOfSprings.size()) + " for " +
                    std::to_string(net.nrOfSprings) + " springs.");
  RUNTIME_EXP_IFN(net.linkIndicesOfSprings.size() == net.nrOfSprings,
                  "Invalid size of link indices of springs, got " +
                    std::to_string(net.linkIndicesOfSprings.size()) + " for " +
                    std::to_string(net.nrOfSprings) + " springs.");
  RUNTIME_EXP_IFN(net.linkIsSliplink.size() ==
                    static_cast<Eigen::Index>(net.nrOfLinks),
                  "Invalid size of link is sliplink");
  RUNTIME_EXP_IFN(net.nrOfCrosslinkSwapsEndured.size() ==
                    static_cast<Eigen::Index>(net.nrOfLinks - net.nrOfNodes),
                  "Invalid size of nr of crosslink swaps endured");
  RUNTIME_EXP_IFN(
    net.linkIsSliplink.count() ==
      static_cast<Eigen::Index>(net.nrOfLinks - net.nrOfNodes),
    "Nr of nodes plus nr of slp-links should give the total nr of links");
  RUNTIME_EXP_IFN(net.oldAtomIds.size() ==
                    static_cast<Eigen::Index>(net.nrOfNodes),
                  "Invalid size of old atom ids");
  RUNTIME_EXP_IFN(net.oldAtomTypes.size() ==
                    static_cast<Eigen::Index>(net.nrOfNodes),
                  "Invalid size of old atom types");
  RUNTIME_EXP_IFN(net.springCoordinateIndexA.size() ==
                    static_cast<Eigen::Index>(net.nrOfSprings * 3),
                  "Invalid size of springCoordinateIndexA");
  RUNTIME_EXP_IFN(net.springCoordinateIndexB.size() ==
                    static_cast<Eigen::Index>(net.nrOfSprings * 3),
                  "Invalid size of springCoordinateIndexB");
  RUNTIME_EXP_IFN(net.springPartCoordinateIndexA.size() ==
                    static_cast<Eigen::Index>(net.nrOfPartialSprings * 3),
                  "Invalid size of springPartCoordinateIndexA");
  RUNTIME_EXP_IFN(net.springPartCoordinateIndexB.size() ==
                    static_cast<Eigen::Index>(net.nrOfPartialSprings * 3),
                  "Invalid size of springPartCoordinateIndexB");
  RUNTIME_EXP_IFN(net.springIndexA.size() ==
                    static_cast<Eigen::Index>(net.nrOfSprings),
                  "Invalid size of springIndexA");
  RUNTIME_EXP_IFN(net.springIndexB.size() ==
                    static_cast<Eigen::Index>(net.nrOfSprings),
                  "Invalid size of springIndexB");
  RUNTIME_EXP_IFN(net.springPartIndexA.size() ==
                    static_cast<Eigen::Index>(net.nrOfPartialSprings),
                  "Invalid size of springPartIndexA");
  RUNTIME_EXP_IFN(net.springPartBoxOffset.size() ==
                    static_cast<Eigen::Index>(net.nrOfPartialSprings * 3),
                  "Invalid size of springPartBoxOffset");
  RUNTIME_EXP_IFN(net.springPartIndexB.size() ==
                    static_cast<Eigen::Index>(net.nrOfPartialSprings),
                  "Invalid size of springPartIndexB");
  RUNTIME_EXP_IFN(net.springIsActive.size() ==
                    static_cast<Eigen::Index>(net.nrOfSprings),
                  "Invalid size of springIsActive");
  RUNTIME_EXP_IFN(springPartitions.size() ==
                    static_cast<Eigen::Index>(net.nrOfPartialSprings),
                  "Invalid size of spring partitions, got " +
                    std::to_string(springPartitions.size()) + " for " +
                    std::to_string(net.nrOfPartialSprings) +
                    " partial springs.");
  RUNTIME_EXP_IFN(net.partialSpringIsPartial.size() ==
                    static_cast<Eigen::Index>(net.nrOfPartialSprings),
                  "Invalid size of partialSpringIsPartial");
  RUNTIME_EXP_IFN(APPROX_EQUAL(springPartitions.sum(),
                               static_cast<double>(net.nrOfSprings),
                               1e-3),
                  "Spring partitions should sum to 1 per spring, got " +
                    std::to_string(springPartitions.sum()) + " for " +
                    std::to_string(net.nrOfSprings) + " springs.");
  RUNTIME_EXP_IFN(
    net.partialToFullSpringIndex.size() ==
      static_cast<Eigen::Index>(net.nrOfPartialSprings),
    "Every partial spring must be able to map to the full spring.");

  /**
   * Test maximum values
   */
  if (net.nrOfSprings > 0) {
    RUNTIME_EXP_IFN(
      static_cast<Eigen::Index>(net.partialToFullSpringIndex.maxCoeff()) <
        static_cast<Eigen::Index>(net.nrOfSprings),
      "Partial spring must map to full spring, which must have "
      "a lower index.");
    RUNTIME_EXP_IFN(
      static_cast<Eigen::Index>(net.springPartCoordinateIndexA.maxCoeff()) <
        static_cast<Eigen::Index>(3 * net.nrOfLinks),
      "Part coordinates must map to coordinates.");
    RUNTIME_EXP_IFN(
      static_cast<Eigen::Index>(net.springPartCoordinateIndexB.maxCoeff()) <
        static_cast<Eigen::Index>(3 * net.nrOfLinks),
      "Part coordinates must map to coordinates.");
    RUNTIME_EXP_IFN(static_cast<Eigen::Index>(net.springPartIndexA.maxCoeff()) <
                      static_cast<Eigen::Index>(net.nrOfLinks),
                    "Part indices must map to links.");
    RUNTIME_EXP_IFN(static_cast<Eigen::Index>(net.springPartIndexB.maxCoeff()) <
                      static_cast<Eigen::Index>(net.nrOfLinks),
                    "Part indices must map to links.");
    RUNTIME_EXP_IFN(static_cast<Eigen::Index>(net.springIndexA.maxCoeff()) <
                      static_cast<Eigen::Index>(net.nrOfNodes),
                    "Full springs must consist of crosslinkers only.");
    RUNTIME_EXP_IFN(static_cast<Eigen::Index>(net.springIndexB.maxCoeff()) <
                      static_cast<Eigen::Index>(net.nrOfNodes),
                    "Full springs must consist of crosslinkers only.");
  }

  /**
   * Test spring partition assumptions
   */
  for (Eigen::Index i = 0; i < springPartitions.size(); i++) {
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
      net.linkIsSliplink[static_cast<Eigen::Index>(link_idx)] ==
        (link_idx >= net.nrOfNodes),
      "Expected slip-links to come sequentially after crosslinkers.");
    std::vector<size_t> thisLinksSprings = net.springIndicesOfLinks[link_idx];
    std::ranges::sort(thisLinksSprings);
    auto last = std::ranges::unique(thisLinksSprings).begin();
    RUNTIME_EXP_IFN(last == thisLinksSprings.end(),
                    "Expect each link to only have one back-link to the "
                    "springs, found back-links " +
                      pylimer_tools::utils::join(thisLinksSprings.begin(),
                                                 thisLinksSprings.end(),
                                                 std::string("_")) +
                      " for link " + std::to_string(link_idx) + ".");
    for (const size_t spring_idx : thisLinksSprings) {
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
    RUNTIME_EXP_IFN(net.linkIsSliplink[static_cast<Eigen::Index>(slipLinkIdx)],
                    "Expected slip-links to know what they are.");
  }

  /**
   * Test the validitiy of springs and their mapping
   */
  Eigen::ArrayXi nrOfMentions =
    Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(net.nrOfLinks));
  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    RUNTIME_EXP_IFN(net.linkIndicesOfSprings[i].size() >= 2,
                    "Each spring requires at least two links, got " +
                      std::to_string(net.linkIndicesOfSprings[i].size()) +
                      " at i = " + std::to_string(i) + ".");
    RUNTIME_EXP_IFN(
      net.springsContourLength[static_cast<Eigen::Index>(i)] > 0,
      "Unexpected spring contour length, got " +
        std::to_string(net.springsContourLength[static_cast<Eigen::Index>(i)]) +
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
      const Eigen::Index partialSpringIdx =
        static_cast<Eigen::Index>(net.localToGlobalSpringIndex[i][partialIdx]);
      const size_t partner0 = net.linkIndicesOfSprings[i][partialIdx];
      const size_t partner1 = net.linkIndicesOfSprings[i][partialIdx + 1];
      RUNTIME_EXP_IFN(
        ((net.springPartIndexA[partialSpringIdx] ==
            static_cast<int>(partner0) &&
          net.springPartIndexB[partialSpringIdx] ==
            static_cast<int>(partner1))),
        "Expect linkIndicesOfSprings and localToGlobalSpringIndex "
        "ordering to correspond. Got partner0 = " +
          std::to_string(partner0) + ", partner1 = " +
          std::to_string(partner1) + " vs. springs part indices " +
          std::to_string(net.springPartIndexA[partialSpringIdx]) + " and " +
          std::to_string(net.springPartIndexB[partialSpringIdx]) +
          " in spring " + std::to_string(i) + " (partial: " +
          std::to_string(partialSpringIdx) + ") with global indices " +
          pylimer_tools::utils::join(net.localToGlobalSpringIndex[i].begin(),
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
      const size_t link_idx = links[j];
      nrOfMentions[static_cast<Eigen::Index>(link_idx)] += 1;
      RUNTIME_EXP_IFN(net.linkIsSliplink[static_cast<Eigen::Index>(link_idx)] ==
                        ((j != 0) && (j != (links.size() - 1))),
                      "Crosslinks must be first and last in a spring, "
                      "slip-links in-between. Found discrepancy at " +
                        std::to_string(j) + "/" + std::to_string(links.size()) +
                        " in spring " + std::to_string(i) + ".")
      std::vector<size_t> thisLinksSprings = net.springIndicesOfLinks[link_idx];
      RUNTIME_EXP_IFN(
        std::find(thisLinksSprings.begin(), thisLinksSprings.end(), i) !=
          thisLinksSprings.end(),
        "Link must have a connection to the spring, too. Did not find "
        "spring " +
          std::to_string(i) + " in link " + std::to_string(link_idx) + ".");
    }
    // also check the sum of the partials
    std::vector<size_t> globalSpringIndices = net.localToGlobalSpringIndex[i];
    double sum = 0.0;
    for (const size_t globalIdx : globalSpringIndices) {
      sum += springPartitions[static_cast<Eigen::Index>(globalIdx)];
    }
    RUNTIME_EXP_IFN(APPROX_EQUAL(sum, 1.0, 1e-10),
                    "Spring partitions of one spring must sum to one, got " +
                      std::to_string(sum) + " for spring " + std::to_string(i) +
                      ".");
  }
  for (Eigen::Index i = static_cast<Eigen::Index>(net.nrOfNodes);
       i < static_cast<Eigen::Index>(net.nrOfLinks);
       ++i) {
    RUNTIME_EXP_IFN(nrOfMentions[i] == 2,
                    "Expect each slip-link to be mentioned twice in the "
                    "links-of-springs mapping, but " +
                      std::to_string(i) + " was mentioned " +
                      std::to_string(nrOfMentions[i]) + " times.");
  }

  /**
   * Test the validity of partial springs and their mapping
   */
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(net.nrOfPartialSprings);
       i++) {
    const Eigen::Index fullIdx =
      static_cast<Eigen::Index>(net.partialToFullSpringIndex[i]);
    const size_t partialEndA = static_cast<size_t>(net.springPartIndexA[i]);
    const size_t partialEndB = static_cast<size_t>(net.springPartIndexB[i]);
    const Eigen::Index partialEndAE = static_cast<Eigen::Index>(partialEndA);
    const Eigen::Index partialEndBE = static_cast<Eigen::Index>(partialEndB);
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
      (net.linkIsSliplink[partialEndAE] || net.linkIsSliplink[partialEndBE]) ==
        net.partialSpringIsPartial[i],
      "Springs involving slip-links must be marked partial. Spring " +
        std::to_string(i) +
        " is marked: " + std::to_string(net.partialSpringIsPartial[i]) + ".");
    RUNTIME_EXP_IFN(
      (net.linkIndicesOfSprings[static_cast<size_t>(fullIdx)].size() > 2) ==
        net.partialSpringIsPartial[i],
      "Springs involving slip-links must be marked partial. Spring " +
        std::to_string(i) +
        " is marked: " + std::to_string(net.partialSpringIsPartial[i]) + ".");
    if (!net.linkIsSliplink[partialEndAE]) {
      RUNTIME_EXP_IFN(
        net.springIndexA[fullIdx] == static_cast<int>(partialEndA),
        "Expect mapping of springs to work: " + std::to_string(partialEndA) +
          " is a crosslink, yet not part of the two ends of spring " +
          std::to_string(fullIdx) + ", where we have " +
          std::to_string(net.springIndexA[fullIdx]) + " and " +
          std::to_string(net.springIndexB[fullIdx]) + ".");
    }
    if (!net.linkIsSliplink[partialEndBE]) {
      RUNTIME_EXP_IFN(
        net.springIndexB[fullIdx] == static_cast<int>(partialEndB),
        "Expect mapping of springs to work: " + std::to_string(partialEndB) +
          " is a crosslink, yet not part of the two ends of spring " +
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
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      RUNTIME_EXP_IFN(
        static_cast<size_t>(net.springPartCoordinateIndexA[3 * i + dir]) ==
          3 * partialEndA + static_cast<size_t>(dir),
        "Spring part index and coordinate index must match. Got " +
          std::to_string(net.springPartCoordinateIndexA[3 * i + dir]) +
          " but expected " +
          std::to_string(3 * partialEndA + static_cast<size_t>(dir)) +
          " with dir = " + std::to_string(dir) + ".");
      RUNTIME_EXP_IFN(
        static_cast<size_t>(net.springPartCoordinateIndexB[3 * i + dir]) ==
          3 * partialEndB + static_cast<size_t>(dir),
        "Spring part index and coordinate index must match. Got " +
          std::to_string(net.springPartCoordinateIndexB[3 * i + dir]) +
          " but expected " +
          std::to_string(3 * partialEndB + static_cast<size_t>(dir)) +
          " with dir = " + std::to_string(dir) + ".");
    }
  }

  /**
   * Check that we do not have any nan or inf values in our vectors
   */
  for (Eigen::Index coordI = 0; coordI < net.coordinates.size(); coordI++) {
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
    RUNTIME_EXP_IFN(APPROX_EQUAL(net.boxHalfs[dir], 0.5 * net.L[dir], 1e-12),
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
      for (const size_t loopIdx : loopsOfSliplink) {
        RUNTIME_EXP_IFN(loopIdx < net.loops.size(), "Loop index out of range.");
      }
    }
    for (std::vector<size_t> loop : net.loops) {
      for (const size_t i : loop) {
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
    const Eigen::Index linkIdxE = static_cast<Eigen::Index>(linkIdx);
    if (net.oldAtomTypes[linkIdxE] == this->entanglementType) {
      RUNTIME_EXP_IFN(
        net.springIndicesOfLinks[linkIdx].size() <= 3,
        "Expect each entanglement atom to have up to three "
        "springs maximum, got " +
          std::to_string(net.springIndicesOfLinks[linkIdx].size()) +
          " for link " + std::to_string(linkIdx) + ".");
      size_t nEntanglementPartners = 0;
      size_t nShortEntanglementPartners = 0;
      for (const size_t springIdx : net.springIndicesOfLinks[linkIdx]) {
        const size_t partnerIdx = this->getOtherEnd(net, springIdx, linkIdx);
        if (partnerIdx <= net.nrOfNodes &&
            net.oldAtomTypes[static_cast<Eigen::Index>(partnerIdx)] ==
              this->entanglementType) {
          nEntanglementPartners += 1;
          nShortEntanglementPartners += static_cast<size_t>(
            net.springsType[static_cast<Eigen::Index>(springIdx)] ==
            this->entanglementType);
        }
      }
      if (net.springIndicesOfLinks[linkIdx].size() == 3) {
        RUNTIME_EXP_IFN(nShortEntanglementPartners == 1,
                        "Each entanglement atom must have exactly one "
                        "entanglement partner. Got " +
                          std::to_string(nShortEntanglementPartners) +
                          " out of " + std::to_string(nEntanglementPartners) +
                          " for link + " + std::to_string(linkIdx) + ".");
      } else {
        RUNTIME_EXP_IFN(
          net.springIndicesOfLinks[linkIdx].size() == 2,
          "Expect each entanglement atom without "
          "entanglement partner to be twofunctional link, got "
          "functionality " +
            std::to_string(net.springIndicesOfLinks[linkIdx].size()) +
            " for link " + std::to_string(linkIdx) + " (originally " +
            std::to_string(net.oldAtomIds[linkIdxE]) + ").");
        RUNTIME_EXP_IFN(nShortEntanglementPartners == 0,
                        "2-functional entanglement atoms are only allowed "
                        "if the entanglement springs has been removed.")
      }
    }
  }
  for (Eigen::Index springIdx = 0;
       springIdx < static_cast<Eigen::Index>(net.nrOfSprings);
       ++springIdx) {
    if (net.springsType[springIdx] == this->entanglementType) {
      RUNTIME_EXP_IFN(net.springsContourLength[springIdx] == 1,
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
