#include "MEHPNonlinearOptimizationSolver.h"
#include "../utils/utilityMacros.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace pylimer_tools::sim::mehp {

namespace {

  [[nodiscard]] bool requiresFiniteBounds(const NLoptAlgorithm alg)
  {
    switch (alg) {
      case NLoptAlgorithm::LD_MMA:
      case NLoptAlgorithm::GN_DIRECT:
      case NLoptAlgorithm::GN_CRS2_LM:
      case NLoptAlgorithm::GD_STOGO:
        return true;
      default:
        return false;
    }
  }

  [[nodiscard]] ExitReason mapNloptResultToExitReason(
    const nlopt::result nloptRes,
    const bool interrupted)
  {
    if (interrupted) {
      return ExitReason::INTERRUPT;
    }
    switch (nloptRes) {
      case nlopt::FTOL_REACHED:
      case nlopt::STOPVAL_REACHED:
        return ExitReason::F_TOLERANCE;
      case nlopt::XTOL_REACHED:
      case nlopt::SUCCESS:
        return ExitReason::X_TOLERANCE;
      case nlopt::MAXEVAL_REACHED:
      case nlopt::MAXTIME_REACHED:
        return ExitReason::MAX_STEPS;
      case nlopt::ROUNDOFF_LIMITED:
        return ExitReason::NO_STEPS_POSSIBLE;
      case nlopt::FORCED_STOP:
        return ExitReason::OTHER;
      case nlopt::FAILURE:
      case nlopt::INVALID_ARGS:
      case nlopt::OUT_OF_MEMORY:
        return ExitReason::FAILURE;
      default:
        return ExitReason::OTHER;
    }
  }

  void zeroZDisplacementsIf2D(Eigen::VectorXd& u,
                              const size_t nrOfLinks,
                              const bool is2D)
  {
    if (!is2D || nrOfLinks == 0) {
      return;
    }
    u.reshaped(3, static_cast<Eigen::Index>(nrOfLinks)).row(2).setZero();
  }

  struct ObjectiveContext
  {
    const MEHPNetwork& net;
    const MEHPForceEvaluator& forceEvaluator;
    const MEHPSolveCallbacks& callbacks;
    nlopt::opt& opt;
    bool is2D;
    size_t totalEvalBase;
    bool interrupted;
    std::exception_ptr callbackException;
    size_t evalCount;
    Eigen::VectorXd uBuffer;

    ObjectiveContext(const MEHPNetwork& netRef,
                     const MEHPForceEvaluator& evaluatorRef,
                     const MEHPSolveCallbacks& callbacksRef,
                     nlopt::opt& optRef,
                     const bool is2DMode,
                     const size_t evalBase)
      : net(netRef)
      , forceEvaluator(evaluatorRef)
      , callbacks(callbacksRef)
      , opt(optRef)
      , is2D(is2DMode)
      , totalEvalBase(evalBase)
      , interrupted(false)
      , callbackException(nullptr)
      , evalCount(0)
      , uBuffer(static_cast<Eigen::Index>(3 * netRef.nrOfLinks))
    {
    }

    ObjectiveContext(const ObjectiveContext&) = delete;
    ObjectiveContext& operator=(const ObjectiveContext&) = delete;
  };

  double nloptObjectiveCallback(const unsigned n,
                                const double* const x,
                                double* const grad,
                                void* const fData)
  {
    ObjectiveContext* const ctx = static_cast<ObjectiveContext*>(fData);
    if (ctx->interrupted) {
      return HUGE_VAL;
    }

    try {
      if (ctx->callbacks.shouldInterrupt && ctx->callbacks.shouldInterrupt()) {
        ctx->interrupted = true;
        ctx->opt.force_stop();
        return HUGE_VAL;
      }

      ctx->uBuffer =
        Eigen::Map<const Eigen::VectorXd>(x, static_cast<Eigen::Index>(n));
      zeroZDisplacementsIf2D(ctx->uBuffer, ctx->net.nrOfLinks, ctx->is2D);

      const Eigen::VectorXd springVectors =
        ctx->net.evaluateSpringVectors(ctx->uBuffer, ctx->is2D);
      const double energy =
        ctx->forceEvaluator.evaluateNetworkEnergyAndGradient(
          ctx->net, springVectors, ctx->is2D, grad);

      RUNTIME_EXP_IFN(
        std::isfinite(energy),
        "Non-finite energy encountered during nonlinear optimization.");
      if (grad != nullptr) {
        Eigen::Map<Eigen::VectorXd> gradMap(grad, static_cast<Eigen::Index>(n));
        if (ctx->is2D && ctx->net.nrOfLinks > 0) {
          gradMap.reshaped(3, static_cast<Eigen::Index>(ctx->net.nrOfLinks))
            .row(2)
            .setZero();
        }
        RUNTIME_EXP_IFN(
          gradMap.allFinite(),
          "Non-finite gradient encountered during nonlinear optimization.");
      }

      ctx->evalCount += 1;
      if (ctx->callbacks.onStep) {
        ctx->callbacks.onStep(ctx->totalEvalBase + ctx->evalCount);
      }
      return energy;
    } catch (...) {
      ctx->callbackException = std::current_exception();
      ctx->opt.force_stop();
      return HUGE_VAL;
    }
  }

} // namespace

nlopt::algorithm
nloptAlgorithmToNlopt(const NLoptAlgorithm alg)
{
  switch (alg) {
    case NLoptAlgorithm::LD_MMA:
      return nlopt::LD_MMA;
    case NLoptAlgorithm::LD_LBFGS:
      return nlopt::LD_LBFGS;
    case NLoptAlgorithm::LD_SLSQP:
      return nlopt::LD_SLSQP;
    case NLoptAlgorithm::LD_TNEWTON:
      return nlopt::LD_TNEWTON;
    case NLoptAlgorithm::LD_TNEWTON_RESTART:
      return nlopt::LD_TNEWTON_RESTART;
    case NLoptAlgorithm::LD_TNEWTON_PRECOND:
      return nlopt::LD_TNEWTON_PRECOND;
    case NLoptAlgorithm::LD_TNEWTON_PRECOND_RESTART:
      return nlopt::LD_TNEWTON_PRECOND_RESTART;
    case NLoptAlgorithm::LD_VAR1:
      return nlopt::LD_VAR1;
    case NLoptAlgorithm::LD_VAR2:
      return nlopt::LD_VAR2;
    case NLoptAlgorithm::LN_SBPLX:
      return nlopt::LN_SBPLX;
    case NLoptAlgorithm::LN_BOBYQA:
      return nlopt::LN_BOBYQA;
    case NLoptAlgorithm::LN_NELDERMEAD:
      return nlopt::LN_NELDERMEAD;
    case NLoptAlgorithm::LN_COBYLA:
      return nlopt::LN_COBYLA;
    case NLoptAlgorithm::LN_NEWUOA:
      return nlopt::LN_NEWUOA;
    case NLoptAlgorithm::LN_NEWUOA_BOUND:
      return nlopt::LN_NEWUOA_BOUND;
    case NLoptAlgorithm::LN_PRAXIS:
      return nlopt::LN_PRAXIS;
    case NLoptAlgorithm::GN_DIRECT:
      return nlopt::GN_DIRECT;
    case NLoptAlgorithm::GN_CRS2_LM:
      return nlopt::GN_CRS2_LM;
    case NLoptAlgorithm::GD_STOGO:
      return nlopt::GD_STOGO;
  }
  throw std::invalid_argument(
    "Unsupported NLoptAlgorithm enumerator with integer value " +
    std::to_string(static_cast<int>(alg)) + ".");
}

NonlinearOptimizationSolver::NonlinearOptimizationSolver(
  const MEHPNonlinearSolverConfig& solverConfig)
  : config(solverConfig)
{
}

const MEHPNonlinearSolverConfig&
NonlinearOptimizationSolver::getConfig() const
{
  return this->config;
}

void
NonlinearOptimizationSolver::setConfig(
  const MEHPNonlinearSolverConfig& newConfig)
{
  this->config = newConfig;
}

MEHPSolveResult
NonlinearOptimizationSolver::solve(MEHPNetwork& net,
                                   Eigen::VectorXd& u,
                                   const MEHPForceEvaluator& forceEvaluator,
                                   const bool is2D,
                                   const MEHPSolveCallbacks& callbacks) const
{
  INVALIDARG_EXP_IFN(
    this->config.slipMode == SlipMode::FIXED,
    "NonlinearOptimizationSolver only supports SlipMode::FIXED.");
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * net.nrOfLinks),
    "Displacement vector u size must match 3 * net.nrOfLinks.");
  INVALIDARG_EXP_IFN(
    this->config.simplificationMode !=
      StructureSimplificationMode::X1F_X2F_THEN_INACTIVE,
    "StructureSimplificationMode::X1F_X2F_THEN_INACTIVE is not implemented; "
    "use INACTIVE_THEN_X2F instead.");
  INVALIDARG_EXP_IFN(
    (this->config.inactiveRemovalCutoff > 0.0 &&
     std::isfinite(this->config.inactiveRemovalCutoff)) ||
      (this->config.simplificationMode !=
         StructureSimplificationMode::INACTIVE_ONLY &&
       this->config.simplificationMode !=
         StructureSimplificationMode::INACTIVE_THEN_X2F),
    "Removal cut-off must be finite and positive when inactive link "
    "simplification is enabled.");
  INVALIDARG_EXP_IFN(this->config.maxNrOfSteps > 0,
                     "maxNrOfSteps must be positive.");
  INVALIDARG_EXP_IFN(std::isfinite(this->config.xtol) &&
                       this->config.xtol >= 0.0,
                     "xtol must be finite and non-negative.");
  INVALIDARG_EXP_IFN(std::isfinite(this->config.ftol) &&
                       this->config.ftol >= 0.0,
                     "ftol must be finite and non-negative.");

  const nlopt::algorithm nloptAlg =
    nloptAlgorithmToNlopt(this->config.algorithm);

  net.validate(u);

  MEHPSolveResult result;
  if (net.nrOfLinks == 0 || net.nrOfSprings == 0 || net.nrOfStrands == 0) {
    result.exitReason = ExitReason::X_TOLERANCE;
    result.energy = 0.0;
    result.residual = 0.0;
    return result;
  }

  if (this->config.simplificationMode !=
      StructureSimplificationMode::NO_SIMPLIFICATION) {
    net.assumeComplete = false;
  }

  const bool boundedDomain = requiresFiniteBounds(this->config.algorithm);
  const int maxEvalsPerPass = static_cast<int>(std::min<long int>(
    this->config.maxNrOfSteps, std::numeric_limits<int>::max()));

  size_t totalEvaluations = 0;
  bool wasInterrupted = false;

  do {
    if (net.nrOfLinks == 0 || net.nrOfSprings == 0 || net.nrOfStrands == 0) {
      result.exitReason = ExitReason::X_TOLERANCE;
      break;
    }

    if (callbacks.shouldInterrupt && callbacks.shouldInterrupt()) {
      wasInterrupted = true;
      break;
    }

    const size_t nDofs = 3 * net.nrOfLinks;
    std::vector<double> uVec(nDofs, 0.0);
    Eigen::Map<Eigen::VectorXd>(uVec.data(), static_cast<Eigen::Index>(nDofs)) =
      u;

    std::vector<double> lowerBounds(nDofs, -HUGE_VAL);
    std::vector<double> upperBounds(nDofs, HUGE_VAL);
    for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
      for (size_t dir = 0; dir < 3; ++dir) {
        const size_t idx = 3 * linkIdx + dir;
        if (is2D && dir == 2) {
          lowerBounds[idx] = 0.0;
          upperBounds[idx] = 0.0;
          uVec[idx] = 0.0;
        } else if (boundedDomain) {
          const double boundHalf =
            std::max(net.L[dir], std::abs(uVec[idx]) + net.L[dir]);
          lowerBounds[idx] = -boundHalf;
          upperBounds[idx] = boundHalf;
        }
      }
    }

    nlopt::opt opt(nloptAlg, static_cast<unsigned int>(nDofs));
    opt.set_exceptions_enabled(false);
    opt.set_lower_bounds(lowerBounds);
    opt.set_upper_bounds(upperBounds);
    opt.set_xtol_rel(this->config.xtol);
    opt.set_ftol_rel(this->config.ftol);
    opt.set_ftol_abs(0.0);
    opt.set_maxeval(maxEvalsPerPass);

    ObjectiveContext ctx(
      net, forceEvaluator, callbacks, opt, is2D, totalEvaluations);
    opt.set_min_objective(nloptObjectiveCallback, &ctx);

    double minEnergy = 0.0;
    const nlopt::result nloptRes = opt.optimize(uVec, minEnergy);

    if (ctx.callbackException != nullptr) {
      std::rethrow_exception(ctx.callbackException);
    }

    totalEvaluations += static_cast<size_t>(std::max(0, opt.get_numevals()));

    u = Eigen::Map<const Eigen::VectorXd>(uVec.data(),
                                          static_cast<Eigen::Index>(nDofs));
    zeroZDisplacementsIf2D(u, net.nrOfLinks, is2D);

    result.exitReason = mapNloptResultToExitReason(nloptRes, ctx.interrupted);
    if (ctx.interrupted) {
      wasInterrupted = true;
      break;
    }

    size_t nTopologyChanges = 0;
    if (this->config.simplificationMode !=
          StructureSimplificationMode::NO_SIMPLIFICATION ||
        this->config.springBreakingLength > 0.0) {
      nTopologyChanges = net.simplify(u,
                                      this->config.simplificationMode,
                                      this->config.inactiveRemovalCutoff,
                                      is2D,
                                      this->config.springBreakingLength);
      if (nTopologyChanges > 0) {
        net.validate(u);
      }
    }

    if (nTopologyChanges == 0) {
      break;
    }
  } while (net.nrOfStrands > 0 && net.nrOfLinks > 0 && net.nrOfSprings > 0);

  net.validate(u);
  result.nrOfIterations = totalEvaluations;

  if (wasInterrupted) {
    result.exitReason = ExitReason::INTERRUPT;
    if (callbacks.cleanupInterrupt) {
      callbacks.cleanupInterrupt();
    }
  } else if (result.exitReason == ExitReason::UNSET) {
    result.exitReason = ExitReason::X_TOLERANCE;
  }

  if (net.nrOfLinks == 0 || net.nrOfSprings == 0) {
    result.energy = 0.0;
    result.residual = 0.0;
  } else {
    Eigen::VectorXd gradVec(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    const Eigen::VectorXd springVectors = net.evaluateSpringVectors(u, is2D);
    result.energy = forceEvaluator.evaluateNetworkEnergyAndGradient(
      net, springVectors, is2D, gradVec.data());
    result.residual = gradVec.norm();
  }

  return result;
}

} // namespace pylimer_tools::sim::mehp
