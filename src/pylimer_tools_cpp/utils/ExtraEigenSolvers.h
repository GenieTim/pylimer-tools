#pragma once

#include <Eigen/Dense>

namespace Eigen {

/**
 * @brief Solves the linear system Ax = b using steepest descent method.
 *
 * This implementation uses optimal step size for each iteration, making it
 * particularly effective for symmetric positive definite (SPD) systems.
 * The method is equivalent to unpreconditioned conjugate gradient with
 * the search direction reset to the residual at each iteration.
 *
 * IMPORTANT: This method assumes A is symmetric positive definite (SPD).
 * For non-SPD systems, use gradientDescentBarzilaiBorwein instead.
 *
 * @tparam MatrixType The type of matrix A (must be SPD)
 * @param A The symmetric positive definite coefficient matrix
 * @param b The right-hand side vector
 * @param learningRate Unused (kept for API compatibility, use
 * gradientDescentBarzilaiBorwein for non-SPD)
 * @param tolerance Convergence tolerance based on relative residual norm
 * @param maxIterations Maximum number of iterations allowed
 * @param iteration Reference to counter storing iterations performed
 * @param initialX Initial solution vector (if empty, uses zero vector)
 * @param initialResidual Initial residual norm (if negative, computes from
 * initial solution)
 * @param iterationCallback Optional callback (return true to stop early)
 * @return The solution vector x approximating the solution to Ax = b
 */
template<typename MatrixType>
static inline Eigen::VectorXd
gradientDescent(const MatrixType& A,
                const Eigen::VectorXd& b,
                double learningRate,
                const double tolerance,
                const int maxIterations,
                int& iteration,
                const Eigen::VectorXd& initialX = Eigen::VectorXd(),
                const double initialResidual = -1.0,
                const std::function<bool(int, const Eigen::VectorXd&)>&
                  iterationCallback = nullptr)
{
  Eigen::VectorXd x =
    initialX.size() == b.size()
      ? initialX
      : Eigen::VectorXd::Zero(b.size()); // Initialize solution vector
  Eigen::VectorXd residual = b - A * x;  // Compute initial residual r = b - Ax
  Eigen::VectorXd Ar(b.size());

  const double initialNorm =
    initialResidual > 0.0 ? initialResidual : residual.norm();
  double residualNorm = residual.norm();

  while ((residualNorm / initialNorm) > tolerance &&
         iteration < maxIterations) {
    Ar.noalias() = A * residual;
    const double rDotR = residual.squaredNorm();
    const double stepSize =
      rDotR / residual.dot(Ar); // Optimal step size for SPD systems
    x += stepSize * residual;
    iteration += 1;

    // Recompute residual periodically to prevent floating point drift
    if (iteration % 50 == 0) {
      residual.noalias() = b - A * x;
    } else {
      // Efficient residual update: r_{k+1} = r_k - alpha * A * r_k
      residual -= stepSize * Ar;
    }
    residualNorm = residual.norm();

    // Call iteration callback if provided
    if (iterationCallback && iterationCallback(iteration, x)) {
      break; // Allow early termination via callback
    }
  }

  return x;
}

/**
 * @brief Solves the linear system Ax = b using the Barzilai-Borwein gradient
 * descent method.
 *
 * The Barzilai-Borwein (BB) method dynamically adjusts the step size to
 * accelerate convergence. Unlike the steepest descent method, BB works for
 * general (non-SPD) matrices. This implementation uses an adaptive step size
 * based on the history of iterates.
 *
 * Two BB step size formulas:
 * - BB1 (short): alpha = (s^T * y) / (y^T * y), where s = x_{k+1} - x_k, y =
 * g_{k+1} - g_k
 * - BB2 (long):  alpha = (s^T * s) / (s^T * y)
 *
 * @tparam MatrixType The type of matrix A (can be non-SPD)
 * @param A The coefficient matrix of the linear system
 * @param b The right-hand side vector
 * @param learningRate Initial learning rate (good default is 0.01-0.1)
 * @param tolerance Convergence tolerance based on relative residual norm
 * @param maxIterations Maximum number of iterations allowed
 * @param iteration Reference to counter storing iterations performed
 * @param shortAlpha If true, uses BB1 (short) step size, otherwise uses BB2
 * (long)
 * @param initialX Initial solution vector (if empty, uses zero vector)
 * @param initialResidual Initial residual norm (if negative, computes from
 * initial solution)
 * @param iterationCallback Optional callback (return true to stop early)
 * @return The solution vector x approximating the solution to Ax = b
 */
template<typename MatrixType>
static inline Eigen::VectorXd
gradientDescentBarzilaiBorwein(
  const MatrixType& A,
  const Eigen::VectorXd& b,
  const double learningRate,
  const double tolerance,
  const int maxIterations,
  int& iteration,
  const bool shortAlpha = false,
  const Eigen::VectorXd& initialX = Eigen::VectorXd(),
  const double initialResidual = -1.0,
  const std::function<bool(int, const Eigen::VectorXd&)>& iterationCallback =
    nullptr)
{
  Eigen::VectorXd x =
    initialX.size() == b.size()
      ? initialX
      : Eigen::VectorXd::Zero(b.size()); // Initialize solution vector
  Eigen::VectorXd gradient = A * x - b; // Compute initial gradient: ∇f = Ax - b
  double alpha = learningRate;          // Initial step size

  Eigen::VectorXd nextX(b.size());
  Eigen::VectorXd nextGradient(b.size());
  Eigen::VectorXd deltaGradient(b.size());

  const double initialNorm =
    initialResidual > 0.0 ? initialResidual : gradient.norm();
  double gradientNorm = gradient.norm();
  while ((gradientNorm / initialNorm) > tolerance &&
         iteration < maxIterations) {
    nextX.noalias() = x - alpha * gradient;
    nextGradient.noalias() = A * nextX - b;
    deltaGradient.noalias() = nextGradient - gradient;

    // Since deltaX = nextX - x = -alpha * gradient:
    // deltaX . deltaX = alpha^2 * gradient.squaredNorm()
    // deltaX . deltaGradient = -alpha * gradient.dot(deltaGradient)
    const double sTs = alpha * alpha * gradient.squaredNorm();
    const double sTy = -alpha * gradient.dot(deltaGradient);
    const double alphaLong = sTs / sTy;
    const double alphaShort = sTy / deltaGradient.squaredNorm();

    x.swap(nextX);
    gradient.swap(nextGradient);
    gradientNorm = gradient.norm();
    alpha = shortAlpha ? alphaShort : alphaLong;
    if (!std::isfinite(alpha)) {
      alpha = learningRate; // Reset alpha to a reasonable value
    }

    iteration++;

    // Call iteration callback if provided
    if (iterationCallback && iterationCallback(iteration, x)) {
      break; // Allow early termination via callback
    }
  }

  return x;
}

/**
 * @brief Solves Ax = b using Barzilai-Borwein with Heavy Ball momentum.
 *
 * This combines the adaptive step size of the Barzilai-Borwein method with
 * Heavy Ball momentum acceleration. The momentum term helps accelerate
 * convergence by incorporating information from previous iterations.
 *
 * Update formula (with momentum active):
 *   x_{k+1} = x_k - (1 - beta) * alpha_k * g_k + beta * (x_k - x_{k-1}),
 * where g_k = A x_k - b and beta is the momentum coefficient (typically
 * 0.5-0.9). The (1 - beta) scaling keeps the steady-state effective step equal
 * to alpha_k. Robustness features:
 *  - alpha_k follows the adaptive (selective long/short) Barzilai-Borwein rule
 *    and keeps the sign of the curvature s^T y, so both positive- and
 *    negative-definite A are supported (the sign is detected after the first
 *    step, which is then discarded from the momentum);
 *  - adaptive restart: momentum is dropped whenever it opposes -alpha_k g_k;
 *  - backtracking with momentum reset if the residual blows up (non-finite or
 *    more than 100x the best residual so far). Each backtracking trial counts
 *    as one iteration.
 *
 * @tparam MatrixType The type of matrix A
 * @param A The coefficient matrix (symmetric definite, either sign)
 * @param b The right-hand side vector
 * @param learningRate Initial (and fallback) step size magnitude (good default
 * is 0.01-0.1)
 * @param tolerance Convergence tolerance based on relative residual norm
 * @param maxIterations Maximum iterations allowed
 * @param iteration Reference to counter storing iterations performed
 * @param initialX Initial solution vector (if empty, uses zero vector)
 * @param initialResidual Initial residual norm (if negative, computes from
 * initial solution)
 * @param iterationCallback Optional callback (return true to stop early)
 * @param momentumCoeff Momentum coefficient beta (default 0.7, range [0,1))
 * @return The solution vector x approximating the solution to Ax = b
 */
template<typename MatrixType>
static inline Eigen::VectorXd
gradientDescentHeavyBallBarzilaiBorwein(
  const MatrixType& A,
  const Eigen::VectorXd& b,
  const double learningRate,
  const double tolerance,
  const int maxIterations,
  int& iteration,
  const Eigen::VectorXd& initialX = Eigen::VectorXd(),
  const double initialResidual = -1.0,
  const std::function<bool(int, const Eigen::VectorXd&)>& iterationCallback =
    nullptr,
  const double momentumCoeff = 0.7)
{
  Eigen::VectorXd x =
    initialX.size() == b.size()
      ? initialX
      : Eigen::VectorXd::Zero(b.size()); // Initialize solution vector
  Eigen::VectorXd gradient = A * x - b;  // Compute initial gradient
  double alpha = learningRate;           // Initial step size
  Eigen::VectorXd momentum = Eigen::VectorXd::Zero(b.size()); // Previous step
  const double beta =
    std::max(0.0, std::min(momentumCoeff, 0.99)); // Clamp momentum

  double gradientNorm = gradient.norm();
  const double initialNorm =
    initialResidual > 0.0 ? initialResidual : gradientNorm;
  if (initialNorm <= 0.0) {
    return x;
  }

  double bestGradientNorm = gradientNorm;
  Eigen::VectorXd step(b.size());
  Eigen::VectorXd nextX(b.size());
  Eigen::VectorXd nextGradient(b.size());
  Eigen::VectorXd deltaGradient(b.size());

  while ((gradientNorm / initialNorm) > tolerance &&
         iteration < maxIterations) {
    // Adaptive gradient restart: reset momentum if it opposes -alpha * gradient
    // (works for both positive-definite alpha > 0 and negative-definite alpha <
    // 0)
    const bool hasMomentum = beta > 0.0 && momentum.squaredNorm() > 0.0 &&
                             (-alpha * gradient.dot(momentum) >= 0.0);
    if (!hasMomentum) {
      momentum.setZero();
    }

    // Heavy Ball update: step = -(1 - beta)*alpha*g + beta*momentum.
    // Scaling the gradient term by (1 - beta) when momentum is active makes the
    // steady-state effective step size match alpha.
    if (hasMomentum) {
      step.noalias() = -(1.0 - beta) * alpha * gradient + beta * momentum;
    } else {
      step.noalias() = -alpha * gradient;
    }

    nextX.noalias() = x + step;
    nextGradient.noalias() = A * nextX - b;
    double nextGradientNorm = nextGradient.norm();
    const double currentGradientNorm = gradientNorm;

    // If the step caused a catastrophic blowup or non-finite residual (and not
    // the initial probing step where alpha sign might flip), reset momentum and
    // backtrack
    if (iteration > 0 && (!std::isfinite(nextGradientNorm) ||
                          nextGradientNorm > 100.0 * bestGradientNorm)) {
      momentum.setZero();
      double safeAlpha = hasMomentum ? alpha : 0.5 * alpha;
      for (int bt = 0; bt < 30 && iteration < maxIterations; ++bt) {
        if (bt > 0) {
          safeAlpha *= 0.5;
        }
        step.noalias() = -safeAlpha * gradient;
        nextX.noalias() = x + step;
        nextGradient.noalias() = A * nextX - b;
        nextGradientNorm = nextGradient.norm();
        iteration++;
        if (std::isfinite(nextGradientNorm) &&
            nextGradientNorm <= currentGradientNorm) {
          break;
        }
      }
      alpha = safeAlpha;
    }

    deltaGradient.noalias() = nextGradient - gradient;

    // Barzilai-Borwein step size (scale-invariant curvature check, preserving
    // sign of sTy for both positive- and negative-definite systems)
    const double sTy = step.dot(deltaGradient);
    const double sTs = step.squaredNorm();
    const double yTy = deltaGradient.squaredNorm();

    double nextAlpha = (sTy < 0.0) ? -learningRate : learningRate;
    if (sTs > 0.0 && yTy > 0.0 &&
        std::abs(sTy) > 1e-12 * std::sqrt(sTs * yTy)) {
      const double alphaLong = sTs / sTy;
      const double alphaShort = sTy / yTy;
      if (std::isfinite(alphaLong) && std::isfinite(alphaShort)) {
        // Selective ABB rule: ratio alphaShort / alphaLong is always in (0, 1]
        nextAlpha = (alphaShort / alphaLong < 0.5) ? alphaShort : alphaLong;
        const double clampedMag =
          std::max(1e-10, std::min(std::abs(nextAlpha), 1e6));
        nextAlpha = std::copysign(clampedMag, nextAlpha);
      }
      // If alpha flipped sign (e.g., on iteration 0 for negative-definite A),
      // reset momentum so the wrong-sign initial step is discarded
      if (nextAlpha * alpha < 0.0) {
        momentum.setZero();
      } else {
        momentum = step;
      }
    } else {
      // Degenerate curvature or zero step: reset momentum
      momentum.setZero();
    }

    if (std::isfinite(nextGradientNorm) &&
        nextGradientNorm < bestGradientNorm) {
      bestGradientNorm = nextGradientNorm;
    }

    // swap instead of copy-assign: next* are loop-local and not reused
    x.swap(nextX);
    gradient.swap(nextGradient);
    gradientNorm = nextGradientNorm;
    alpha = nextAlpha;
    iteration++;

    // Call iteration callback if provided
    if (iterationCallback && iterationCallback(iteration, x)) {
      break; // Allow early termination via callback
    }
  }

  return x;
}

}
