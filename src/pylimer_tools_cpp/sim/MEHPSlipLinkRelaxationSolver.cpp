#include "MEHPSlipLinkRelaxationSolver.h"
#include "../utils/utilityMacros.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace pylimer_tools::sim::mehp {

namespace {

  constexpr size_t kDefaultSimplificationFrequency = 10;
  constexpr int kMaxInnerIterations = 50;

  /**
   * @brief Compute the L2 residual norm of the net link forces for a
   * linear Hookean network with clamped partial-spring contour weights.
   */
  [[nodiscard]] double computeDisplacementResidualNorm(
    const MEHPNetwork& net,
    const Eigen::VectorXd& u,
    const double kappa,
    const bool is2D,
    const double oneOverSpringPartitionUpperLimit)
  {
    if (net.nrOfLinks == 0 || net.nrOfSprings == 0 || kappa == 0.0) {
      return 0.0;
    }

    const Eigen::VectorXd springVectors = net.evaluateSpringVectors(u, is2D);
    Eigen::VectorXd forces =
      Eigen::VectorXd::Zero(3 * static_cast<Eigen::Index>(net.nrOfLinks));

    for (size_t s = 0; s < net.nrOfSprings; ++s) {
      const Eigen::Index sE = static_cast<Eigen::Index>(s);
      const int idxA = net.springIndexA[sE];
      const int idxB = net.springIndexB[sE];
      if (idxA < 0 || idxB < 0 || idxA == idxB) {
        continue;
      }

      const int strandIdxInt = net.strandIndexOfSpring[sE];
      const bool isPartial =
        (strandIdxInt >= 0 && !net.springIsEntanglement[sE] &&
         net.springIndicesOfStrand[static_cast<size_t>(strandIdxInt)].size() >
           1);
      const double strandN =
        (strandIdxInt >= 0)
          ? net.strandContourLength[static_cast<Eigen::Index>(strandIdxInt)]
          : net.springContourLength[sE];
      const double invContour =
        MEHPNetwork::computeClampedInverseSpringContourLength(
          isPartial,
          net.springContourLength[sE],
          strandN,
          oneOverSpringPartitionUpperLimit);
      const double kOverN = kappa * invContour;
      const Eigen::Vector3d forceVec =
        kOverN * springVectors.segment<3>(3 * sE);

      const Eigen::Index a3 = 3 * static_cast<Eigen::Index>(idxA);
      const Eigen::Index b3 = 3 * static_cast<Eigen::Index>(idxB);
      forces.segment<3>(a3) += forceVec;
      forces.segment<3>(b3) -= forceVec;
    }

    return forces.norm();
  }

  /**
   * @brief Relax interior entanglement link partitions/positions and then
   * crosslink positions in Gauss-Seidel passes without assuming contiguous
   * index ranges for crosslinks vs. slip-links.
   */
  [[nodiscard]] double performSlidingRelaxationPass(
    MEHPNetwork& net,
    Eigen::VectorXd& u,
    const MEHPSlipLinkSolverConfig& cfg,
    const bool is2D)
  {
    double maxDistanceMoved = 0.0;

    for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
      if (!net.linkIsEntanglement[static_cast<Eigen::Index>(linkIdx)]) {
        continue;
      }
      int innerIterationsDone = 0;
      do {
        const double partitionRes = net.updateSpringPartitionAtLink(
          linkIdx, u, cfg.oneOverSpringPartitionUpperLimit, is2D);
        const double displacementDone = net.displaceToMeanPosition(
          u, linkIdx, is2D, cfg.oneOverSpringPartitionUpperLimit);
        maxDistanceMoved = std::max(maxDistanceMoved, displacementDone);
        innerIterationsDone += 1;
        if (partitionRes == 0.0 && displacementDone == 0.0) {
          break;
        }
      } while (cfg.doInnerIterations &&
               innerIterationsDone < kMaxInnerIterations);
    }

    for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
      if (net.linkIsEntanglement[static_cast<Eigen::Index>(linkIdx)]) {
        continue;
      }
      const double displacementDone = net.displaceToMeanPosition(
        u, linkIdx, is2D, cfg.oneOverSpringPartitionUpperLimit);
      maxDistanceMoved = std::max(maxDistanceMoved, displacementDone);
    }

    return maxDistanceMoved;
  }

} // namespace

SlipLinkRelaxationSolver::SlipLinkRelaxationSolver(
  const MEHPSlipLinkSolverConfig& solverConfig)
  : config(solverConfig)
{
}

MEHPSolveResult
SlipLinkRelaxationSolver::solve(MEHPNetwork& net,
                                Eigen::VectorXd& u,
                                const MEHPForceEvaluator& forceEvaluator,
                                const bool is2D,
                                const MEHPSolveCallbacks& callbacks) const
{
  // 1. Validate preconditions and configuration early
  INVALIDARG_EXP_IFN(
    forceEvaluator.isLinear(),
    "SlipLinkRelaxationSolver requires a linear Hookean force evaluator.");
  INVALIDARG_EXP_IFN(
    !(this->config.slipMode == SlipMode::SLIDING &&
      net.springIsEntanglement.size() > 0 && net.springIsEntanglement.any()),
    "SlipMode::SLIDING cannot be used when the network contains spring-mode "
    "entanglement coupling springs; use SlipMode::FIXED or "
    "EntanglementMode::MERGED_LINKS.");
  INVALIDARG_EXP_IFN(u.size() == static_cast<Eigen::Index>(3 * net.nrOfLinks),
                     "Displacement vector u size (" + std::to_string(u.size()) +
                       ") does not match 3 * net.nrOfLinks (" +
                       std::to_string(3 * net.nrOfLinks) + ").");
  INVALIDARG_EXP_IFN(
    this->config.simplificationMode !=
      StructureSimplificationMode::X1F_X2F_THEN_INACTIVE,
    "StructureSimplificationMode::X1F_X2F_THEN_INACTIVE is not supported; "
    "use INACTIVE_THEN_X2F instead.");
  INVALIDARG_EXP_IFN(
    std::isfinite(this->config.inactiveRemovalCutoff) &&
      (this->config.inactiveRemovalCutoff > 0.0 ||
       this->config.simplificationMode ==
         StructureSimplificationMode::NO_SIMPLIFICATION),
    "inactiveRemovalCutoff must be finite and positive when simplification is "
    "enabled.");
  INVALIDARG_EXP_IFN(this->config.maxNrOfSteps >= 0,
                     "maxNrOfSteps must be non-negative.");
  INVALIDARG_EXP_IFN(std::isfinite(this->config.xtol) &&
                       this->config.xtol >= 0.0,
                     "xtol must be finite and non-negative.");
  INVALIDARG_EXP_IFN(std::isfinite(this->config.initialResidualToUse),
                     "initialResidualToUse must be finite.");
  INVALIDARG_EXP_IFN(
    std::isfinite(this->config.oneOverSpringPartitionUpperLimit),
    "oneOverSpringPartitionUpperLimit must be finite.");
  INVALIDARG_EXP_IFN(!std::isnan(this->config.springBreakingLength),
                     "springBreakingLength must not be NaN.");

  const double kappa = forceEvaluator.getKappa();
  INVALIDARG_EXP_IFN(std::isfinite(kappa) && kappa >= 0.0,
                     "Force evaluator kappa must be finite and non-negative.");

  net.validate(u);

  MEHPSolveResult result;
  result.oneOverSpringPartitionUpperLimitUsed =
    this->config.oneOverSpringPartitionUpperLimit;

  // 2. Empty network fast path
  if (net.nrOfLinks == 0 || net.nrOfSprings == 0) {
    result.exitReason = ExitReason::X_TOLERANCE;
    result.energy = 0.0;
    result.residual = 0.0;
    return result;
  }

  if (is2D) {
    for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
      u[3 * static_cast<Eigen::Index>(linkIdx) + 2] = 0.0;
    }
  }

  // 3. Zero spring constant implies zero forces everywhere
  if (kappa == 0.0) {
    result.exitReason = ExitReason::X_TOLERANCE;
    result.residual = 0.0;
    result.energy = forceEvaluator.evaluateNetworkEnergyAndGradient(
      net, net.evaluateSpringVectors(u, is2D), is2D, nullptr);
    return result;
  }

  const double actualInitialResidual = computeDisplacementResidualNorm(
    net, u, kappa, is2D, this->config.oneOverSpringPartitionUpperLimit);
  double referenceResidual = (this->config.initialResidualToUse > 0.0)
                               ? this->config.initialResidualToUse
                               : actualInitialResidual;
  double currentResidual = actualInitialResidual;

  size_t iterationsDone = 0;
  bool wasInterrupted = false;
  int totalRelocations = 0;
  int totalShifts = 0;

  const bool hasSimplification =
    (this->config.simplificationMode !=
     StructureSimplificationMode::NO_SIMPLIFICATION) ||
    (this->config.springBreakingLength > 0.0);

  if (this->config.maxNrOfSteps > 0) {
    do {
      const bool canSlipAndSwap =
        (this->config.slipMode == SlipMode::SLIDING) &&
        (net.nrOfLinks > net.nrOfNodes);

      // Periodic slip-link swapping / crosslink branch-hopping pass
      if (canSlipAndSwap &&
          this->config.linkSwappingMode != LinkSwappingMode::NO_SWAPPING &&
          this->config.swappingFrequency > 0 &&
          (iterationsDone %
           static_cast<size_t>(this->config.swappingFrequency)) == 0) {
        const long int crosslinkSwapsBefore =
          static_cast<long int>(net.nrOfCrosslinkSwapsEndured.sum());
        const size_t totalSwapsInPass = net.performSlipLinkSwaps(
          u,
          this->config.linkSwappingMode,
          this->config.oneOverSpringPartitionUpperLimit,
          this->config.nrOfCrosslinkSwapsAllowedPerSlipLink,
          is2D);
        const long int crosslinkSwapsAfter =
          static_cast<long int>(net.nrOfCrosslinkSwapsEndured.sum());
        const long int relocationsInPass =
          std::max<long int>(0, crosslinkSwapsAfter - crosslinkSwapsBefore);
        const long int shiftsInPass = std::max<long int>(
          0, static_cast<long int>(totalSwapsInPass) - relocationsInPass);
        totalRelocations += static_cast<int>(relocationsInPass);
        totalShifts += static_cast<int>(shiftsInPass);
      }

      // Partition & coordinate update pass (re-check nrOfLinks > nrOfNodes)
      const bool useSlidingPass =
        (this->config.slipMode == SlipMode::SLIDING) &&
        (net.nrOfLinks > net.nrOfNodes);
      const double maxDistanceMoved =
        useSlidingPass
          ? performSlidingRelaxationPass(net, u, this->config, is2D)
          : net.displaceToMeanPosition(
              u, is2D, this->config.oneOverSpringPartitionUpperLimit);

      RUNTIME_EXP_IFN(std::isfinite(maxDistanceMoved),
                      "Non-finite displacement computed in iteration " +
                        std::to_string(iterationsDone + 1) + ".");

      currentResidual = computeDisplacementResidualNorm(
        net, u, kappa, is2D, this->config.oneOverSpringPartitionUpperLimit);
      RUNTIME_EXP_IFN(std::isfinite(currentResidual),
                      "Non-finite residual norm computed in iteration " +
                        std::to_string(iterationsDone + 1) + ".");

      iterationsDone += 1;

      const bool convergedBeforeSimplification =
        (currentResidual <= this->config.xtol * referenceResidual);
      if (hasSimplification &&
          ((iterationsDone % kDefaultSimplificationFrequency) == 0 ||
           convergedBeforeSimplification)) {
        const size_t nChanges = net.simplify(u,
                                             this->config.simplificationMode,
                                             this->config.inactiveRemovalCutoff,
                                             is2D,
                                             this->config.springBreakingLength);
        if (nChanges > 0) {
          currentResidual = computeDisplacementResidualNorm(
            net, u, kappa, is2D, this->config.oneOverSpringPartitionUpperLimit);
          RUNTIME_EXP_IFN(
            std::isfinite(currentResidual),
            "Non-finite residual norm after simplification in iteration " +
              std::to_string(iterationsDone) + ".");
          if (referenceResidual <= 0.0 && currentResidual > 0.0) {
            referenceResidual = currentResidual;
          }
        }
      }

      if (callbacks.onStep) {
        callbacks.onStep(iterationsDone);
      }
      if (callbacks.shouldInterrupt && callbacks.shouldInterrupt()) {
        wasInterrupted = true;
        break;
      }
    } while (currentResidual > this->config.xtol * referenceResidual &&
             static_cast<long int>(iterationsDone) <
               this->config.maxNrOfSteps &&
             net.nrOfSprings > 0);
  }

#ifndef NDEBUG
  net.validate(u);
#endif

  result.nrOfIterations = iterationsDone;
  result.residual = currentResidual;
  result.nrOfRelocations = totalRelocations;
  result.nrOfShifts = totalShifts;
  result.energy = (net.nrOfSprings > 0)
                    ? forceEvaluator.evaluateNetworkEnergyAndGradient(
                        net, net.evaluateSpringVectors(u, is2D), is2D, nullptr)
                    : 0.0;

  if (wasInterrupted) {
    result.exitReason = ExitReason::INTERRUPT;
    if (callbacks.cleanupInterrupt) {
      callbacks.cleanupInterrupt();
    }
  } else if (net.nrOfSprings == 0 ||
             currentResidual <= this->config.xtol * referenceResidual) {
    result.exitReason = ExitReason::X_TOLERANCE;
  } else {
    result.exitReason = ExitReason::MAX_STEPS;
  }

  return result;
}

} // namespace pylimer_tools::sim::mehp
