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
  X(SPARSE_SLICING, "Sparse Slicing")                                          \
  X(METAL, "Metal GPU")

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

#define STOCHASTIC_SPECTRAL_METHODS                                            \
  X(AUTO, "Auto")                                                              \
  X(SLQ, "SLQ")                                                                \
  X(HYBRID, "Hybrid")                                                          \
  X(KPM, "KPM")

/**
 * @brief Available stochastic / approximate spectral methods for normal mode
 *        analysis on large networks.
 */
enum class StochasticSpectralMethod
{
#define X(e, s) e,
  STOCHASTIC_SPECTRAL_METHODS
#undef X
};

inline const std::string StochasticSpectralMethodNames[] = {
#define X(e, s) s,
  STOCHASTIC_SPECTRAL_METHODS
#undef X
};

inline const StochasticSpectralMethod allStochasticSpectralMethods[] = {
#define X(e, s) StochasticSpectralMethod::e,
  STOCHASTIC_SPECTRAL_METHODS
#undef X
};

class NormalModeAnalyzer
{
public:
  NormalModeAnalyzer(const std::vector<size_t>& springFrom,
                     const std::vector<size_t>& springTo,
                     size_t nrOfNodes = 0);

  void findSparseEigenvalues(const size_t nrOfEigenvalues,
                             const bool includeEigenvectors = false);

  /**
   * @brief Compute all eigenvalues using the specified solver.
   *
   * @param solver Eigenvalue solver method (AUTO, DENSE_LAPACK, DENSE_EIGEN,
   * SPARSE_SLICING, or METAL)
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

  const Eigen::VectorXd& getEigenvalues() const;
  void setEigenvalues(const Eigen::VectorXd& e);

  const Eigen::MatrixXd& getEigenvectors() const;
  void setEigenvectors(const Eigen::MatrixXd& e);

  const Eigen::SparseMatrix<double>& getAssembledConnectivityMatrix() const;

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
   * @brief Compute Chebyshev moments using the Kernel Polynomial Method (KPM)
   *        accelerated with Apple Silicon Metal compute kernels.
   *
   * @param numMoments Number of Chebyshev moments (default: 1024)
   * @param numRandomVectors Number of random vectors (default: 20)
   * @param batchSize Batch size for multi-RHS block SpMM (default: 4)
   * @return Vector of damped Chebyshev moments
   */
  Eigen::VectorXd computeChebyshevMomentsKPMMetal(size_t numMoments = 1024,
                                                  size_t numRandomVectors = 20,
                                                  size_t batchSize = 4) const;

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

  // MARK: Approximate / Stochastic Spectral Methods (SLQ, Hybrid, KPM)

  /**
   * @brief Evaluate stress autocorrelation G(t) using an approximate spectral
   *        method (SLQ, Hybrid exact-low + SLQ, or KPM).
   *
   * @param t Time points at which to evaluate stress autocorrelation
   * @param method Spectral approximation method (AUTO, SLQ, HYBRID, or KPM)
   * @param numSteps Lanczos steps (SLQ/HYBRID) or Chebyshev moments (KPM)
   * @param numRandomVectors Number of Rademacher random vectors (default: 20)
   * @param numExactLowModes Number of exact lowest nonzero modes for HYBRID
   * (default: 100)
   * @return Stress autocorrelation values
   */
  Eigen::ArrayXd evaluateStressAutocorrelationApprox(
    const Eigen::ArrayXd& t,
    StochasticSpectralMethod method = StochasticSpectralMethod::AUTO,
    size_t numSteps = 200,
    size_t numRandomVectors = 20,
    size_t numExactLowModes = 100) const;

  /**
   * @brief Evaluate storage modulus G'(omega) using an approximate spectral
   *        method (SLQ, Hybrid exact-low + SLQ, or KPM).
   *
   * @param omega Angular frequencies to evaluate
   * @param method Spectral approximation method (AUTO, SLQ, HYBRID, or KPM)
   * @param numSteps Lanczos steps (SLQ/HYBRID) or Chebyshev moments (KPM)
   * @param numRandomVectors Number of Rademacher random vectors (default: 20)
   * @param numExactLowModes Number of exact lowest nonzero modes for HYBRID
   * (default: 100)
   * @return Storage modulus values
   */
  Eigen::ArrayXd evaluateStorageModulusApprox(
    const Eigen::ArrayXd& omega,
    StochasticSpectralMethod method = StochasticSpectralMethod::AUTO,
    size_t numSteps = 200,
    size_t numRandomVectors = 20,
    size_t numExactLowModes = 100) const;

  /**
   * @brief Evaluate loss modulus G''(omega) using an approximate spectral
   *        method (SLQ, Hybrid exact-low + SLQ, or KPM).
   *
   * @param omega Angular frequencies to evaluate
   * @param method Spectral approximation method (AUTO, SLQ, HYBRID, or KPM)
   * @param numSteps Lanczos steps (SLQ/HYBRID) or Chebyshev moments (KPM)
   * @param numRandomVectors Number of Rademacher random vectors (default: 20)
   * @param numExactLowModes Number of exact lowest nonzero modes for HYBRID
   * (default: 100)
   * @return Loss modulus values
   */
  Eigen::ArrayXd evaluateLossModulusApprox(
    const Eigen::ArrayXd& omega,
    StochasticSpectralMethod method = StochasticSpectralMethod::AUTO,
    size_t numSteps = 200,
    size_t numRandomVectors = 20,
    size_t numExactLowModes = 100) const;

#ifdef CEREALIZABLE
  static NormalModeAnalyzer fromString(std::string in)
  {
    NormalModeAnalyzer n;
    pylimer_tools::utils::deserializeFromString(n, in);
    return n;
  }
#endif

  size_t getMatrixSize() const
  {
    return static_cast<size_t>(assembledConnectivityMatrix.rows());
  }

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
  mutable size_t cachedKPMMoments_{ 0 };
  mutable size_t cachedKPMVectors_{ 0 };
  mutable Eigen::VectorXd cachedMu_{};

  mutable bool hasCachedApproxSpectrum_{ false };
  mutable StochasticSpectralMethod cachedApproxMethod_{
    StochasticSpectralMethod::AUTO
  };
  mutable size_t cachedApproxSteps_{ 0 };
  mutable size_t cachedApproxVectors_{ 0 };
  mutable size_t cachedApproxLowModes_{ 0 };
  mutable std::vector<double> cachedApproxLambdaNodes_{};
  mutable std::vector<double> cachedApproxDeltaN_{};

  void solveMatrix(const Eigen::SparseMatrix<double>& mat,
                   NormalModeSolver solver,
                   bool includeEigenvectors,
                   Eigen::VectorXd& outEvals,
                   Eigen::MatrixXd& outEvecs);

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

  void computeConnectedComponents(size_t& numComponents,
                                  std::vector<size_t>& labels,
                                  std::vector<size_t>& counts) const;

  void runDeflatedSLQ(size_t numSteps,
                      size_t numRandomVectors,
                      size_t numComponents,
                      const std::vector<size_t>& labels,
                      const std::vector<size_t>& counts,
                      const Eigen::MatrixXd* deflationBasis,
                      std::vector<double>& lambdaNodes,
                      std::vector<double>& deltaN) const;

  void computeEffectiveSpectrumHybrid(size_t numSteps,
                                      size_t numRandomVectors,
                                      size_t numExactLowModes,
                                      std::vector<double>& lambdaNodes,
                                      std::vector<double>& deltaN) const;

  void ensureApproxSpectrumCached(StochasticSpectralMethod method,
                                  size_t numSteps,
                                  size_t numRandomVectors,
                                  size_t numExactLowModes) const;

  static Eigen::ArrayXd evaluateStressAutocorrelationFromSpectrum(
    const Eigen::ArrayXd& t,
    const std::vector<double>& lambdaNodes,
    const std::vector<double>& deltaN);

  static Eigen::ArrayXd evaluateStorageModulusFromSpectrum(
    const Eigen::ArrayXd& omega,
    const std::vector<double>& lambdaNodes,
    const std::vector<double>& deltaN);

  static Eigen::ArrayXd evaluateLossModulusFromSpectrum(
    const Eigen::ArrayXd& omega,
    const std::vector<double>& lambdaNodes,
    const std::vector<double>& deltaN);

  static Eigen::VectorXd computeJacksonKernel(size_t numMoments);

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
       eigenvectors,
       assembledConnectivityMatrix);
  }
#endif
};
}

#endif
