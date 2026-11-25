#ifndef MEHP_SLIP_LINK_RELAXATION_SOLVER_H
#define MEHP_SLIP_LINK_RELAXATION_SOLVER_H

#include "./MEHPSolver.h"

namespace pylimer_tools::sim::mehp {

/**
 * @brief Iterative alternating coordinate and contour-partition relaxation
 * solver for slip-link and linear Hookean polymer networks.
 *
 * Alternates between updating interior entanglement contour partitions
 * (when @ref SlipMode::SLIDING is active), moving links toward local weighted
 * equilibrium positions, performing optional slip-link swaps and crosslink
 * branch hops, and applying strand breaking or topological simplification.
 */
class SlipLinkRelaxationSolver final : public MEHPSolver
{
private:
  MEHPSlipLinkSolverConfig config = {};

public:
  /**
   * @brief Construct a slip-link relaxation solver with the given
   * configuration.
   *
   * @param config Solver configuration parameters.
   */
  explicit SlipLinkRelaxationSolver(
    const MEHPSlipLinkSolverConfig& config = {});

  /**
   * @brief Return the current solver configuration.
   *
   * @return Const reference to the stored @ref MEHPSlipLinkSolverConfig.
   */
  [[nodiscard]] const MEHPSlipLinkSolverConfig& getConfig() const
  {
    return this->config;
  }

  /**
   * @brief Replace the solver configuration.
   *
   * @param newConfig New configuration parameters.
   */
  void setConfig(const MEHPSlipLinkSolverConfig& newConfig)
  {
    this->config = newConfig;
  }

  /**
   * @brief Relax link displacements and (when sliding) slip-link contour
   * partitions to mechanical equilibrium on a linear Hookean network.
   *
   * @param net Canonical polymer network (modified in-place by partition
   *        updates, slip-link swaps, strand breaking, and/or simplification)
   * @param u Link displacement vector of size `3 * net.nrOfLinks`, updated
   *        in-place (and resized if simplification prunes zero-functional
   * links)
   * @param forceEvaluator Linear constitutive evaluator (`isLinear() == true`)
   * @param is2D If true, keep z-displacements and z-forces constrained to zero
   * @param callbacks Optional interrupt and per-step monitoring hooks
   * @return MEHPSolveResult Convergence metrics, exit reason, energy, and swap
   *         counters
   * @throws std::invalid_argument on violated preconditions or invalid config
   * @throws std::runtime_error on network validation failure or non-finite
   * state
   */
  [[nodiscard]] MEHPSolveResult solve(
    MEHPNetwork& net,
    Eigen::VectorXd& u,
    const MEHPForceEvaluator& forceEvaluator,
    bool is2D,
    const MEHPSolveCallbacks& callbacks = {}) const override;
};

} // namespace pylimer_tools::sim::mehp

#endif // MEHP_SLIP_LINK_RELAXATION_SOLVER_H
