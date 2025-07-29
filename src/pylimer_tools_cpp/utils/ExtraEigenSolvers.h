#pragma once

#include <Eigen/Dense>

namespace Eigen {

// Function to perform Gradient Descent with residual-based stopping criterion
// as documented in Wikipedia
// TODO: Could benefit from a preconditioner
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
  Eigen::VectorXd gradient = b - A * x;  // Compute initial residual

  const double initialNorm =
    initialResidual > 0.0 ? initialResidual : gradient.squaredNorm();
  while ((gradient.squaredNorm() / initialNorm) > tolerance &&
         iteration < maxIterations) {
    Eigen::VectorXd Ar = A * gradient;
    double stepSize = gradient.dot(gradient) / gradient.dot(Ar);
    x = x + stepSize * gradient;
    iteration += 1;

    // reset to "correct" gradient after every 100 iterations
    // against floating point precision issues
    if (iteration % 100 == 0) {
      gradient = b - A * x;
    } else {
      // otherwise, stick with the more efficient gradient update
      gradient = gradient - stepSize * Ar;
    }

    // Call iteration callback if provided
    if (iterationCallback && iterationCallback(iteration, x)) {
      break; // Allow early termination via callback
    }
  }

  return x;
}

// Function to perform Gradient Descent with residual-based stopping criterion
// Barzilai-Borwein method to solve Ax = b
// TODO: Could benefit from a preconditioner
/**
 * @brief Solves the linear system Ax = b using the Barzilai-Borwein gradient
 * descent method.
 *
 * This implementation uses a residual-based stopping criterion and dynamically
 * adjusts the step size using the Barzilai-Borwein method, which can accelerate
 * convergence compared to standard gradient descent.
 *
 * @tparam MatrixType The type of matrix A (must support matrix-vector
 * multiplication)
 * @param A The coefficient matrix of the linear system
 * @param b The right-hand side vector of the linear system
 * @param learningRate Initial learning rate (good default is 0.01)
 * @param tolerance The convergence tolerance based on relative residual norm
 * @param maxIterations Maximum number of iterations allowed
 * @param iteration Reference to a counter that will store the number of
 * iterations performed
 * @param shortAlpha If true, uses the "short" step length formula, otherwise
 * uses the "long" formula
 * @param initialX Initial solution vector (if empty, uses zero vector)
 * @param initialResidual Initial residual norm (if negative, computes from
 * initial solution)
 * @param iterationCallback Optional callback function called each iteration
 * (return true to stop)
 * @return Eigen::VectorXd The solution vector x that approximates the solution
 * to Ax = b
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
  Eigen::VectorXd gradient = A * x - b;  // Compute initial residual
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
    nullptr)
{
  Eigen::VectorXd x =
    initialX.size() == b.size()
      ? initialX
      : Eigen::VectorXd::Zero(b.size()); // Initialize solution vector
  Eigen::VectorXd gradient = A * x - b;  // Compute initial residual
  double alpha = learningRate;           // Initial step size
  Eigen::VectorXd deltaX =
    Eigen::VectorXd::Zero(b.size()); // Initialize deltaX, the momentum

  const double initialNorm =
    initialResidual > 0.0 ? initialResidual : gradient.squaredNorm();
  while ((gradient.squaredNorm() / initialNorm) > tolerance &&
         iteration < maxIterations) {
    Eigen::VectorXd nextX = x - alpha * gradient + alpha * deltaX;
    Eigen::VectorXd nextGradient = A * nextX - b;

    deltaX = nextX - x;
    Eigen::VectorXd deltaGradient = nextGradient - gradient;

    const double alphaLong = (deltaX.dot(deltaX)) / (deltaX.dot(deltaGradient));
    const double alphaShort =
      (deltaX.dot(deltaGradient)) / (deltaGradient.dot(deltaGradient));

    x = nextX;
    gradient = nextGradient;

    if (std::isfinite(alphaShort) && std::isfinite(alphaLong)) {
      alpha = (alphaShort + alphaLong) / 2.0;
    } else if (std::isfinite(alphaShort)) {
      alpha = alphaShort;
    } else if (std::isfinite(alphaLong)) {
      alpha = alphaLong;
    } else {
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

}
