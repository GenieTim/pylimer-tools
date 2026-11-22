#ifndef MEHP_NONLINEAR_OPTIMIZATION_SOLVER_H
#define MEHP_NONLINEAR_OPTIMIZATION_SOLVER_H

#include "./MEHPSolver.h"
#include <nlopt.hpp>

namespace pylimer_tools::sim::mehp {

/**
 * @brief Convert an @ref NLoptAlgorithm enumerator to its @c nlopt::algorithm
 * equivalent.
 *
 * @param alg Strongly-typed NLopt algorithm enumerator.
 * @return nlopt::algorithm Matching NLopt backend algorithm constant.
 * @throws std::invalid_argument if @p alg is not a supported enumerator.
 */
[[nodiscard]] nlopt::algorithm
nloptAlgorithmToNlopt(NLoptAlgorithm alg);

/**
 * @brief General nonlinear energy minimization solver for MEHP networks using
 * NLopt.
 *
 * Supports fixed entanglement links (`SlipMode::FIXED`) with any
 * @ref EntanglementMode (`NONE`, `MERGED_LINKS`, `SPRING`), arbitrary linear or
 * nonlinear C++/Python @ref MEHPForceEvaluator implementations, optional strand
 * breaking, and iterative topology simplification.
 */
class NonlinearOptimizationSolver final : public MEHPSolver
{
private:
  MEHPNonlinearSolverConfig config;

public:
  /**
   * @brief Construct a nonlinear optimization MEHP solver with the given
   * configuration.
   *
   * @param config NLopt algorithm, tolerances, step limit, and simplification
   *        configuration.
   */
  explicit NonlinearOptimizationSolver(
    const MEHPNonlinearSolverConfig& config = {});

  /**
   * @brief Return the current nonlinear solver configuration.
   *
   * @return Const reference to the stored @ref MEHPNonlinearSolverConfig.
   */
  [[nodiscard]] const MEHPNonlinearSolverConfig& getConfig() const;

  /**
   * @brief Replace the nonlinear solver configuration.
   *
   * @param newConfig New configuration to apply on subsequent @ref solve calls.
   */
  void setConfig(const MEHPNonlinearSolverConfig& newConfig);

  /**
   * @brief Solve for static energy minimum on an @ref MEHPNetwork using NLopt.
   *
   * @param net Canonical polymer network (modified in-place if strand breaking
   *        or topology simplification is enabled).
   * @param u Link displacement vector of size `3 * net.nrOfLinks`, updated
   *        in-place with converged equilibrium displacements.
   * @param forceEvaluator Constitutive spring evaluator (Hookean, non-Gaussian,
   *        or custom Python/C++ evaluator).
   * @param is2D If true, constrain z-displacements and z-gradients to zero.
   * @param callbacks Optional interrupt, cleanup, and per-evaluation callbacks.
   * @return MEHPSolveResult Summary containing exit reason, evaluation count,
   *         L2 gradient residual norm, and minimized elastic energy.
   * @throws std::invalid_argument on violated preconditions or invalid config.
   * @throws std::runtime_error on network validation failure or non-finite energy.
   */
  [[nodiscard]] MEHPSolveResult solve(
    MEHPNetwork& net,
    Eigen::VectorXd& u,
    const MEHPForceEvaluator& forceEvaluator,
    bool is2D,
    const MEHPSolveCallbacks& callbacks = {}) const override;
};

} // namespace pylimer_tools::sim::mehp

#endif // MEHP_NONLINEAR_OPTIMIZATION_SOLVER_H
