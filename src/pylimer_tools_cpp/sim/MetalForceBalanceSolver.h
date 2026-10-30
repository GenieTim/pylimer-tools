#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace pylimer_tools::sim::mehp {

/**
 * @brief Checks if Metal GPU acceleration is available on this system.
 */
bool
isMetalAvailable();

/**
 * @brief Solves the sparse linear equation system Ax = b on Apple Silicon GPU
 * using the Barzilai-Borwein gradient descent method implemented with Metal
 * compute kernels.
 *
 * @param A Coefficient matrix (sparse)
 * @param b Right-hand side vector
 * @param learningRate Initial learning rate
 * @param tolerance Relative residual tolerance
 * @param maxIterations Maximum number of iterations
 * @param iterationsDone Output count of iterations performed
 * @param shortAlpha Whether to use BB1 (short) step size formula instead of BB2
 * (long)
 * @param initialX Optional initial guess vector
 * @return Solution vector x
 */
Eigen::VectorXd
solveSLEMetalBarzilaiBorwein(
  const Eigen::SparseMatrix<double>& A,
  const Eigen::VectorXd& b,
  double learningRate,
  double tolerance,
  int maxIterations,
  int& iterationsDone,
  bool shortAlpha = false,
  const Eigen::VectorXd& initialX = Eigen::VectorXd());

/**
 * @brief Solves the sparse linear equation system Ax = b on Apple Silicon GPU
 * using the Conjugate Gradient (CG) method implemented with Metal compute
 * kernels.
 *
 * @param A Symmetric positive (semi-)definite coefficient matrix
 * @param b Right-hand side vector
 * @param tolerance Relative residual tolerance
 * @param maxIterations Maximum number of iterations
 * @param iterationsDone Output count of iterations performed
 * @param initialX Optional initial guess vector
 * @return Solution vector x
 */
Eigen::VectorXd
solveSLEMetalCG(const Eigen::SparseMatrix<double>& A,
                const Eigen::VectorXd& b,
                double tolerance,
                int maxIterations,
                int& iterationsDone,
                const Eigen::VectorXd& initialX = Eigen::VectorXd());

} // namespace pylimer_tools::sim::mehp
