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
 * @param learningRate Unused (kept for API compatibility, use gradientDescentBarzilaiBorwein for non-SPD)
 * @param tolerance Convergence tolerance based on relative residual norm
 * @param maxIterations Maximum number of iterations allowed
 * @param iteration Reference to counter storing iterations performed
 * @param initialX Initial solution vector (if empty, uses zero vector)
 * @param initialResidual Initial residual norm (if negative, computes from initial solution)
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

  const double initialNorm =
    initialResidual > 0.0 ? initialResidual : residual.norm();
  
  while ((residual.norm() / initialNorm) > tolerance &&
         iteration < maxIterations) {
    Eigen::VectorXd Ar = A * residual;
    double stepSize = residual.dot(residual) / residual.dot(Ar);  // Optimal step size for SPD systems
    x = x + stepSize * residual;
    iteration += 1;

    // Recompute residual periodically to prevent floating point drift
    if (iteration % 50 == 0) {
      residual = b - A * x;
    } else {
      // Efficient residual update: r_{k+1} = r_k - alpha * A * r_k
      residual = residual - stepSize * Ar;
    }

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
 * - BB1 (short): alpha = (s^T * y) / (y^T * y), where s = x_{k+1} - x_k, y = g_{k+1} - g_k
 * - BB2 (long):  alpha = (s^T * s) / (s^T * y)
 *
 * @tparam MatrixType The type of matrix A (can be non-SPD)
 * @param A The coefficient matrix of the linear system
 * @param b The right-hand side vector
 * @param learningRate Initial learning rate (good default is 0.01-0.1)
 * @param tolerance Convergence tolerance based on relative residual norm
 * @param maxIterations Maximum number of iterations allowed
 * @param iteration Reference to counter storing iterations performed
 * @param shortAlpha If true, uses BB1 (short) step size, otherwise uses BB2 (long)
 * @param initialX Initial solution vector (if empty, uses zero vector)
 * @param initialResidual Initial residual norm (if negative, computes from initial solution)
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
  Eigen::VectorXd gradient = A * x - b;  // Compute initial gradient: ∇f = Ax - b
  double alpha = learningRate;           // Initial step size

  const double initialNorm =
    initialResidual > 0.0 ? initialResidual : gradient.squaredNorm();
  while ((gradient.squaredNorm() / initialNorm) > tolerance &&
         iteration < maxIterations) {
    Eigen::VectorXd nextX = x - alpha * gradient;
    Eigen::VectorXd nextGradient = A * nextX - b;

    Eigen::VectorXd deltaX = nextX - x;
    Eigen::VectorXd deltaGradient = nextGradient - gradient;

    const double alphaLong = (deltaX.dot(deltaX)) / (deltaX.dot(deltaGradient));
    const double alphaShort =
      (deltaX.dot(deltaGradient)) / (deltaGradient.dot(deltaGradient));

    x = nextX;
    gradient = nextGradient;
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
 * Update formula: x_{k+1} = x_k - alpha_k * g_k + beta * (x_k - x_{k-1})
 * where beta is the momentum coefficient (typically 0.5-0.9).
 *
 * @tparam MatrixType The type of matrix A
 * @param A The coefficient matrix
 * @param b The right-hand side vector
 * @param learningRate Initial learning rate (good default is 0.01-0.1)
 * @param tolerance Convergence tolerance based on relative residual norm
 * @param maxIterations Maximum iterations allowed
 * @param iteration Reference to counter storing iterations performed
 * @param initialX Initial solution vector (if empty, uses zero vector)
 * @param initialResidual Initial residual norm (if negative, computes from initial solution)
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
  const double beta = std::max(0.0, std::min(momentumCoeff, 0.99)); // Clamp momentum

  const double initialNorm =
    initialResidual > 0.0 ? initialResidual : gradient.norm();
  
  Eigen::VectorXd Ax = A * x;  // Cache A*x
  
  while ((gradient.norm() / initialNorm) > tolerance &&
         iteration < maxIterations) {
    // Heavy Ball update: step = -alpha*g + beta*momentum
    Eigen::VectorXd step = -alpha * gradient + beta * momentum;
    x = x + step;
    
    // Efficient gradient computation
    Eigen::VectorXd Astep = A * step;
    Eigen::VectorXd nextAx = Ax + Astep;
    Eigen::VectorXd nextGradient = nextAx - b;
    
    Eigen::VectorXd deltaGradient = nextGradient - gradient;

    // Barzilai-Borwein step size (using BB2 which often works better with momentum)
    double sTy = step.dot(deltaGradient);
    double sTs = step.squaredNorm();
    
    if (std::abs(sTy) > 1e-14 && std::isfinite(sTs / sTy)) {
      alpha = sTs / sTy;
      // Clamp to reasonable bounds
      alpha = std::max(1e-8, std::min(alpha, 1e6));
    } else {
      alpha = learningRate;
    }

    // Update for next iteration
    momentum = step;
    Ax = nextAx;
    gradient = nextGradient;
    iteration++;

    // Periodically recompute to prevent floating point drift
    if (iteration % 50 == 0) {
      Ax = A * x;
      gradient = Ax - b;
      // Optionally reset momentum on recomputation to avoid accumulated errors
      if (iteration % 200 == 0) {
        momentum.setZero();
      }
    }

    // Call iteration callback if provided
    if (iterationCallback && iterationCallback(iteration, x)) {
      break; // Allow early termination via callback
    }
  }

  return x;
}

}
