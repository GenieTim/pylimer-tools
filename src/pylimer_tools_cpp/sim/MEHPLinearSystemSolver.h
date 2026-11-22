#ifndef MEHP_LINEAR_SYSTEM_SOLVER_H
#define MEHP_LINEAR_SYSTEM_SOLVER_H

#include "./MEHPSolver.h"

namespace pylimer_tools::sim::mehp {

/**
 * @brief Equilibrium solver for linear (Hookean) MEHP networks via sparse
 * graph-Laplacian linear system assembly and direct/iterative/GPU SLE solvers.
 *
 * Supports fixed entanglement links (`SlipMode::FIXED`), optional strand
 * breaking with incremental in-place CSR matrix updates, and iterative network
 * topology simplification.
 */
class LinearSystemMEHPSolver final : public MEHPSolver
{
private:
  MEHPLinearSolverConfig config;

public:
  /**
   * @brief Construct a linear-system MEHP solver with the given configuration.
   *
   * @param config Solver, simplification, and fallback configuration.
   */
  explicit LinearSystemMEHPSolver(const MEHPLinearSolverConfig& config = {});

  /**
   * @brief Return the current linear solver configuration.
   *
   * @return Const reference to the stored @ref MEHPLinearSolverConfig.
   */
  [[nodiscard]] const MEHPLinearSolverConfig& getConfig() const;

  /**
   * @brief Replace the linear solver configuration.
   *
   * @param newConfig New configuration to apply on subsequent @ref solve calls.
   */
  void setConfig(const MEHPLinearSolverConfig& newConfig);

  /**
   * @brief Solve for static force equilibrium on a linear Hookean @ref MEHPNetwork.
   *
   * @param net Canonical polymer network (modified in-place if strand breaking
   *        or topology simplification is enabled).
   * @param u Link displacement vector of size `3 * net.nrOfLinks`, updated
   *        in-place with converged equilibrium displacements.
   * @param forceEvaluator Constitutive model; must satisfy `isLinear() == true`.
   * @param is2D If true, constrain z-displacements and z-forces to zero.
   * @param callbacks Optional interrupt, cleanup, and step callbacks.
   * @return MEHPSolveResult Summary containing exit reason, iteration count,
   *         squared force residual norm, elastic energy, and fallback status.
   * @throws std::invalid_argument on violated preconditions or invalid config.
   * @throws std::runtime_error on network validation failure or linear solver error.
   */
  [[nodiscard]] MEHPSolveResult solve(
    MEHPNetwork& net,
    Eigen::VectorXd& u,
    const MEHPForceEvaluator& forceEvaluator,
    bool is2D,
    const MEHPSolveCallbacks& callbacks = {}) const override;
};

} // namespace pylimer_tools::sim::mehp

#endif // MEHP_LINEAR_SYSTEM_SOLVER_H
