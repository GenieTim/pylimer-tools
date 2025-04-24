#pragma once

#include <Eigen/Dense>

namespace Eigen {
// Function to perform Gradient Descent with residual-based stopping criterion
// Barzilai-Borwein method to solve Ax = b
template<typename MatrixType>
static inline Eigen::VectorXd
gradientDescentBarzilaiBorwein(const MatrixType& A,
                               const Eigen::VectorXd& b,
                               double learningRate,
                               const double tolerance,
                               const int maxIterations,
                               int& iteration)
{
  int n = A.cols();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n); // Initialize solution vector
  Eigen::VectorXd gradient = A * x - b;         // Compute initial residual
  double alpha = 0.1;                           // Initial step size

  Eigen::VectorXd nextX = x - alpha * gradient;
  Eigen::VectorXd nextGradient = A * nextX - b;

  const double initialNorm = gradient.squaredNorm();
  while ((nextGradient.squaredNorm() / initialNorm) > tolerance &&
         iteration < maxIterations) {
    Eigen::VectorXd deltaX = nextX - x;
    Eigen::VectorXd deltaGradient = nextGradient - gradient;

    double alphaLong = (deltaX.dot(deltaX)) / (deltaX.dot(deltaGradient));
    double alphaShort =
      (deltaX.dot(deltaGradient)) / (deltaGradient.dot(deltaGradient));

    x = nextX;
    gradient = nextGradient;
    nextX = x - alphaLong * gradient;
    nextGradient = A * nextX - b;

    iteration++;
  }

  //  std::cout << "Converged in " << iteration << " iterations.\n";
  return x;
}

template<typename MatrixType>
static inline Eigen::VectorXd
gradientDescent(const MatrixType& A,
                const Eigen::VectorXd& b,
                double learningRate,
                const double tolerance,
                const int maxIterations,
                int& iteration)
{
  int n = A.cols();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n); // Initialize solution vector
  Eigen::VectorXd gradient = b - A * x;         // Compute initial residual

  const double initialNorm = gradient.squaredNorm();
  while ((gradient.squaredNorm() / initialNorm) > tolerance &&
         iteration < maxIterations) {
    Eigen::VectorXd Ar = A * gradient;
    double stepSize = gradient.dot(gradient) / gradient.dot(Ar);
    x = x + stepSize * gradient;
    iterations += 1;
    // reset to "correct" gradient after every 100 iterations
    // against floating point precision issues
    if (iteration % 100 == 0) {
      gradient = b - A * x;
    } else {
      // otherwise, stick with the more efficient gradient update
      gradient = gradient - stepSize * Ar;
    }
  }

  return x;
}
}
