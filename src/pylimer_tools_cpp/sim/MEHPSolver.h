#ifndef MEHP_SOLVER_H
#define MEHP_SOLVER_H

#include "./MEHPForceEvaluator.h"
#include "./MEHPNetwork.h"
#include "./MEHPUtilityStructures.h"
#include <Eigen/Dense>
#include <cstddef>
#include <functional>

namespace pylimer_tools::sim::mehp {

/**
 * @brief Configuration for the sparse linear-system MEHP equilibrium solver.
 *
 * Precondition: requires a linear (Hookean) force evaluator and
 * `slipMode == SlipMode::FIXED`.
 */
struct MEHPLinearSolverConfig
{
  SLESolver solver = SLESolver::DEFAULT;
  StructureSimplificationMode simplificationMode =
    StructureSimplificationMode::NO_SIMPLIFICATION;
  double inactiveRemovalCutoff = 1e-6;
  double residualReduction = 1e-12;
  int maxIterations = 50000;
  double springBreakingLength = -1.0;
  double inPlaceRebuildThreshold = 0.15;
  int cgFallbackMaxIterations = 20000;
  double gradientDescentLearningRate = 0.01;
  SlipMode slipMode = SlipMode::FIXED;
};

/**
 * @brief Configuration for alternating coordinate/partition slip-link
 * relaxation.
 *
 * Precondition: requires a linear (Hookean) force evaluator. When
 * `slipMode == SlipMode::SLIDING`, spring-mode entanglement coupling springs
 * (`springIsEntanglement.any()`) are disallowed.
 */
struct MEHPSlipLinkSolverConfig
{
  long int maxNrOfSteps = 50000;
  double xtol = 1e-9;
  double initialResidualToUse = -1.0;
  StructureSimplificationMode simplificationMode =
    StructureSimplificationMode::NO_SIMPLIFICATION;
  double inactiveRemovalCutoff = 1e-3;
  bool doInnerIterations = false;
  LinkSwappingMode linkSwappingMode = LinkSwappingMode::NO_SWAPPING;
  int swappingFrequency = 10;
  double oneOverSpringPartitionUpperLimit = 1.0;
  int nrOfCrosslinkSwapsAllowedPerSlipLink = -1;
  SlipMode slipMode = SlipMode::SLIDING;
  double springBreakingLength = -1.0;
};

/**
 * @brief Configuration for general nonlinear optimization MEHP equilibrium
 * solving via NLopt.
 *
 * Precondition: requires `slipMode == SlipMode::FIXED`. Supports both Hookean
 * and arbitrary nonlinear constitutive models on phantom and fixed-entanglement
 * networks.
 */
struct MEHPNonlinearSolverConfig
{
  NLoptAlgorithm algorithm = NLoptAlgorithm::LD_MMA;
  long int maxNrOfSteps = 250000;
  double xtol = 1e-12;
  double ftol = 1e-9;
  StructureSimplificationMode simplificationMode =
    StructureSimplificationMode::NO_SIMPLIFICATION;
  double inactiveRemovalCutoff = 1e-6;
  double springBreakingLength = -1.0;
  SlipMode slipMode = SlipMode::FIXED;
};

/**
 * @brief Optional runtime callbacks passed to an MEHP solver invocation.
 */
struct MEHPSolveCallbacks
{
  std::function<bool()> shouldInterrupt = []() { return false; };
  std::function<void()> cleanupInterrupt = []() {};
  std::function<void(size_t)> onStep = nullptr;
};

/**
 * @brief Diagnostic summary returned by an @ref MEHPSolver execution.
 */
struct MEHPSolveResult
{
  ExitReason exitReason = ExitReason::UNSET;
  size_t nrOfIterations = 0;
  double residual = 0.0;
  double energy = 0.0;
  bool cgFallbackTriggered = false;
  SLESolver effectiveSLESolver = SLESolver::DEFAULT;
  int nrOfRelocations = 0;
  int nrOfShifts = 0;
  double oneOverSpringPartitionUpperLimitUsed = -1.0;

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(exitReason,
       nrOfIterations,
       residual,
       energy,
       cgFallbackTriggered,
       effectiveSLESolver,
       nrOfRelocations,
       nrOfShifts,
       oneOverSpringPartitionUpperLimitUsed);
  }
#endif
};

/**
 * @brief Abstract strategy interface for solving mechanical equilibrium on an
 * @ref MEHPNetwork under a constitutive @ref MEHPForceEvaluator.
 */
class MEHPSolver
{
protected:
  MEHPSolver(const MEHPSolver&) = default;
  MEHPSolver& operator=(const MEHPSolver&) = default;
  MEHPSolver(MEHPSolver&&) noexcept = default;
  MEHPSolver& operator=(MEHPSolver&&) noexcept = default;

public:
  MEHPSolver() = default;
  virtual ~MEHPSolver() = default;

  /**
   * @brief Solve for equilibrium link displacements (and optionally contour
   * partitions / topology simplifications) starting from @p u.
   *
   * @param net Canonical polymer network (may be simplified or have slip-link
   *        partitions updated in-place depending on solver strategy and config)
   * @param u Link displacement vector of size `3 * net.nrOfLinks`, updated
   *        in-place with the converged equilibrium displacements
   * @param forceEvaluator Stateless constitutive spring model
   * @param is2D If true, keep z-displacements and z-forces constrained to zero
   * @param callbacks Optional interrupt and per-step monitoring hooks
   * @return MEHPSolveResult Convergence metrics and exit reason
   * @throws std::invalid_argument on violated solver preconditions or invalid
   *         input dimensions
   */
  [[nodiscard]] virtual MEHPSolveResult solve(
    MEHPNetwork& net,
    Eigen::VectorXd& u,
    const MEHPForceEvaluator& forceEvaluator,
    bool is2D,
    const MEHPSolveCallbacks& callbacks = {}) const = 0;
};

} // namespace pylimer_tools::sim::mehp

#endif // MEHP_SOLVER_H
