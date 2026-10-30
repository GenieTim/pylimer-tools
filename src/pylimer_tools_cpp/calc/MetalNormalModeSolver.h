#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace pylimer_tools::calc {

/**
 * @brief Checks if Metal GPU acceleration is available on this system.
 */
bool
isMetalAvailable();

/**
 * @brief Computes all eigenvalues and optionally eigenvectors of a real
 * symmetric matrix on Apple Silicon GPU using Metal compute kernels.
 *
 * @param mat The input symmetric matrix (sparse)
 * @param includeEigenvectors Whether to compute eigenvectors
 * @param outEvals Output vector for sorted eigenvalues
 * @param outEvecs Output matrix for eigenvectors
 */
void
computeAllEigenvaluesMetal(const Eigen::SparseMatrix<double>& mat,
                           const bool includeEigenvectors,
                           Eigen::VectorXd& outEvals,
                           Eigen::MatrixXd& outEvecs);

/**
 * @brief Computes Chebyshev moments for Kernel Polynomial Method (KPM)
 * directly on Apple Silicon GPU using Metal compute kernels for SpMV,
 * vector recurrence, and parallel dot product reductions.
 *
 * @param mat The assembled connectivity matrix
 * @param numMoments Number of Chebyshev moments
 * @param numRandomVectors Number of random vectors for stochastic trace
 * estimation
 * @param c Spectral shift parameter
 * @param w Spectral scale parameter
 * @return Damped Chebyshev moments vector
 */
Eigen::VectorXd
computeChebyshevMomentsKPMMetal(const Eigen::SparseMatrix<double>& mat,
                                size_t numMoments,
                                size_t numRandomVectors,
                                double c,
                                double w);

} // namespace pylimer_tools::calc
