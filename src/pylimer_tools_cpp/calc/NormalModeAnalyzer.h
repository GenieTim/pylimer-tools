#ifndef NORMAL_MODE_ANALYZER_H
#define NORMAL_MODE_ANALYZER_H

#include "../utils/CerealUtils.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifdef CEREALIZABLE
#include <cereal/access.hpp>
#endif
#include <string>
#include <vector>

namespace pylimer_tools::calc {

#define NORMAL_MODE_SOLVERS                                                    \
  X(AUTO, "Auto")                                                              \
  X(DENSE_LAPACK, "Dense LAPACK")                                              \
  X(DENSE_EIGEN, "Dense Eigen")                                                \
  X(SPARSE_SLICING, "Sparse Slicing")

/**
 * @brief Available eigenvalue solver algorithms for normal mode analysis.
 */
enum class NormalModeSolver
{
#define X(e, s) e,
  NORMAL_MODE_SOLVERS
#undef X
};

inline const std::string NormalModeSolverNames[] = {
#define X(e, s) s,
  NORMAL_MODE_SOLVERS
#undef X
};

inline const NormalModeSolver allNormalModeSolvers[] = {
#define X(e, s) NormalModeSolver::e,
  NORMAL_MODE_SOLVERS
#undef X
};

class NormalModeAnalyzer
{
public:
  NormalModeAnalyzer(const std::vector<size_t> springFrom,
                     const std::vector<size_t> springTo);

  void findSparseEigenvalues(const size_t nrOfEigenvalues,
                             const bool includeEigenvectors = false);

  /**
   * @brief Compute all eigenvalues using the specified solver.
   *
   * @param solver Eigenvalue solver method (AUTO, DENSE_LAPACK, DENSE_EIGEN, or
   * SPARSE_SLICING)
   * @param includeEigenvectors Whether to also compute and store eigenvectors
   * @param decomposeComponents Whether to decompose the graph into connected
   * components
   */
  void computeAllEigenvalues(
    const NormalModeSolver solver = NormalModeSolver::AUTO,
    const bool includeEigenvectors = false,
    const bool decomposeComponents = true);

  void computeAllEigenvalues(const bool includeEigenvectors)
  {
    computeAllEigenvalues(NormalModeSolver::AUTO, includeEigenvectors, true);
  }

  void computeAllEigenvaluesGPU(const bool includeEigenvectors = false);

  Eigen::VectorXd getEigenvalues() const;
  void setEigenvalues(Eigen::VectorXd e);

  Eigen::MatrixXd getEigenvectors() const;
  void setEigenvectors(Eigen::MatrixXd e);

  Eigen::SparseMatrix<double> getAssembledConnectivityMatrix() const;

  Eigen::ArrayXd evaluateStressAutocorrelation(const Eigen::ArrayXd& t) const;

  Eigen::ArrayXd evaluateStorageModulus(const Eigen::ArrayXd& omega) const;

  Eigen::ArrayXd evaluateLossModulus(const Eigen::ArrayXd& omega) const;

  size_t getNrOfSolubleClusters() const;

  // MARK: KPM (Kernel Polynomial Method) methods

  /**
   * @brief Compute Chebyshev moments using the Kernel Polynomial Method (KPM)
   *        with Rademacher random vector trace estimation and Jackson damping
   * kernel.
   *
   * @param numMoments Number of Chebyshev moments (default: 1024)
   * @param numRandomVectors Number of random vectors for stochastic trace
   * estimation (default: 20)
   * @return Vector of damped Chebyshev moments
   */
  Eigen::VectorXd computeChebyshevMomentsKPM(
    size_t numMoments = 1024,
    size_t numRandomVectors = 20) const;

  /**
   * @brief Compute continuous density of states (DOS) using KPM.
   *
   * @param lambdaGrid Eigenvalue points at which to evaluate DOS
   * @param numMoments Number of Chebyshev moments (default: 1024)
   * @param numRandomVectors Number of random vectors (default: 20)
   * @return Density of states values on lambdaGrid
   */
  Eigen::ArrayXd computeDensityOfStatesKPM(const Eigen::ArrayXd& lambdaGrid,
                                           size_t numMoments = 1024,
                                           size_t numRandomVectors = 20) const;

  /**
   * @brief Evaluate stress autocorrelation G(t) directly using KPM via
   * Gauss-Chebyshev quadrature.
   *
   * @param t Time points at which to evaluate stress autocorrelation
   * @param numMoments Number of Chebyshev moments (default: 1024)
   * @param numRandomVectors Number of random vectors (default: 20)
   * @param numQuadraturePoints Number of Gauss-Chebyshev quadrature points
   * (default: 2048)
   * @return Stress autocorrelation values
   */
  Eigen::ArrayXd evaluateStressAutocorrelationKPM(
    const Eigen::ArrayXd& t,
    size_t numMoments = 1024,
    size_t numRandomVectors = 20,
    size_t numQuadraturePoints = 2048) const;

  /**
   * @brief Evaluate storage modulus G'(omega) directly using KPM via
   * Gauss-Chebyshev quadrature.
   *
   * @param omega Angular frequencies to evaluate
   * @param numMoments Number of Chebyshev moments (default: 1024)
   * @param numRandomVectors Number of random vectors (default: 20)
   * @param numQuadraturePoints Number of Gauss-Chebyshev quadrature points
   * (default: 2048)
   * @return Storage modulus values
   */
  Eigen::ArrayXd evaluateStorageModulusKPM(
    const Eigen::ArrayXd& omega,
    size_t numMoments = 1024,
    size_t numRandomVectors = 20,
    size_t numQuadraturePoints = 2048) const;

  /**
   * @brief Evaluate loss modulus G''(omega) directly using KPM via
   * Gauss-Chebyshev quadrature.
   *
   * @param omega Angular frequencies to evaluate
   * @param numMoments Number of Chebyshev moments (default: 1024)
   * @param numRandomVectors Number of random vectors (default: 20)
   * @param numQuadraturePoints Number of Gauss-Chebyshev quadrature points
   * (default: 2048)
   * @return Loss modulus values
   */
  Eigen::ArrayXd evaluateLossModulusKPM(
    const Eigen::ArrayXd& omega,
    size_t numMoments = 1024,
    size_t numRandomVectors = 20,
    size_t numQuadraturePoints = 2048) const;

#ifdef CEREALIZABLE
  static NormalModeAnalyzer fromString(std::string in)
  {
    NormalModeAnalyzer n;
    pylimer_tools::utils::deserializeFromString(n, in);
    return n;
  }
#endif

  size_t getMatrixSize() const { return assembledConnectivityMatrix.rows(); }

protected:
  void requireEigenvaluesComputation() const;
  void requireEigenvectorsComputation() const;
  size_t countSolubleClusters() const;

private:
  Eigen::SparseMatrix<double> assembledConnectivityMatrix;
  // computation state
  bool isEigenvaluesComputed = false;
  bool isEigenvectorsComputed = false;
  size_t clusterCount = 0;
  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd eigenvectors;

  void computeAllEigenvaluesDense(const Eigen::SparseMatrix<double>& mat,
                                  const bool includeEigenvectors,
                                  const bool forceEigenSolver,
                                  Eigen::VectorXd& outEvals,
                                  Eigen::MatrixXd& outEvecs);

  void computeAllEigenvaluesSparseSlicing(
    const Eigen::SparseMatrix<double>& mat,
    const bool includeEigenvectors,
    Eigen::VectorXd& outEvals,
    Eigen::MatrixXd& outEvecs);

  std::pair<double, double> getSpectralScalingParameters(
    double eps = 0.01) const;

  void computeEffectiveSpectrumKPM(size_t numMoments,
                                   size_t numRandomVectors,
                                   size_t numQuadraturePoints,
                                   std::vector<double>& lambdaNodes,
                                   std::vector<double>& deltaN) const;

  // MARK: serialization
  NormalModeAnalyzer() = default;
#ifdef CEREALIZABLE
  friend class cereal::access;
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(isEigenvaluesComputed,
       isEigenvectorsComputed,
       clusterCount,
       eigenvalues,
       eigenvectors);
  }
#endif
};
}

#endif
