#ifndef NORMAL_MODE_ANALYZER_H
#define NORMAL_MODE_ANALYZER_H

#include "../utils/CerealUtils.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifdef CEREALIZABLE
#include <cereal/access.hpp>
#endif
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace pylimer_tools::entities {
class Universe;
}

namespace pylimer_tools::calc {

#define NORMAL_MODE_SOLVERS                                                    \
  X(AUTO, "Auto")                                                              \
  X(DENSE_LAPACK, "Dense LAPACK")                                              \
  X(DENSE_EIGEN, "Dense Eigen")                                                \
  X(SPARSE_SLICING, "Sparse Slicing")                                          \
  X(METAL, "Metal GPU")                                                        \
  X(CUDA, "CUDA GPU")

/**
 * @brief Available exact eigenvalue solver algorithms for normal mode analysis.
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

#define COMPUTE_BACKENDS                                                       \
  X(AUTO, "Auto")                                                              \
  X(CPU, "CPU")                                                                \
  X(METAL, "Metal")                                                            \
  X(CUDA, "CUDA")

/**
 * @brief Hardware compute backend for spectral approximations.
 */
enum class ComputeBackend
{
#define X(e, s) e,
  COMPUTE_BACKENDS
#undef X
};

inline const std::string ComputeBackendNames[] = {
#define X(e, s) s,
  COMPUTE_BACKENDS
#undef X
};

/**
 * @brief Configuration for Stochastic Lanczos Quadrature (SLQ).
 */
struct SLQConfig
{
  size_t numSteps = 200;
  size_t numRandomVectors = 20;
  ComputeBackend backend = ComputeBackend::AUTO;
};

/**
 * @brief Configuration for Hybrid exact low-mode + deflated SLQ approximation.
 */
struct HybridConfig
{
  size_t numSteps = 200;
  size_t numRandomVectors = 20;
  size_t numExactLowModes = 100;
  ComputeBackend backend = ComputeBackend::AUTO;
};

/**
 * @brief Configuration for Kernel Polynomial Method (KPM) spectral
 *        approximation.
 */
struct KPMConfig
{
  size_t numMoments = 1024;
  size_t numRandomVectors = 20;
  size_t numQuadraturePoints = 2048;
  size_t batchSize = 4;
  ComputeBackend backend = ComputeBackend::AUTO;
};

/**
 * @brief Immutable value object representing a discrete or quadrature-weighted
 *        normal mode spectrum $\{(\lambda_i, w_i)\}$ and optional eigenvectors.
 *
 * Evaluates viscoelastic observables ($G(t)$, $G'(\omega)$, $G''(\omega)$)
 * uniformly across exact, sparse, SLQ, Hybrid, and KPM spectra.
 */
class NormalModeSpectrum
{
public:
  /**
   * @brief Construct a NormalModeSpectrum from eigenvalues, optional weights,
   *        and optional eigenvectors.
   *
   * @param eigenvalues Vector of eigenvalues $\lambda_i$
   * @param weights Optional quadrature weights $w_i$ (defaults to all 1.0)
   * @param eigenvectors Optional eigenvector matrix (columns correspond to
   * modes)
   * @param isComplete Whether this represents the full spectrum of the system
   */
  explicit NormalModeSpectrum(
    Eigen::VectorXd eigenvalues,
    std::optional<Eigen::VectorXd> weights = std::nullopt,
    std::optional<Eigen::MatrixXd> eigenvectors = std::nullopt,
    bool isComplete = true);

  const Eigen::VectorXd& getEigenvalues() const;
  const Eigen::VectorXd& getWeights() const;
  std::optional<Eigen::MatrixXd> getEigenvectors() const;
  bool isComplete() const;
  size_t getNumModes() const;
  size_t getNumZeroModes() const;

  /**
   * @brief Evaluate stress autocorrelation $G(t) = \sum_i w_i e^{-2\lambda_i
   * t}$.
   */
  Eigen::ArrayXd stressAutocorrelation(
    const Eigen::Ref<const Eigen::ArrayXd>& t) const;

  /**
   * @brief Evaluate storage modulus $G'(\omega)$.
   */
  Eigen::ArrayXd storageModulus(
    const Eigen::Ref<const Eigen::ArrayXd>& omega) const;

  /**
   * @brief Evaluate loss modulus $G''(\omega)$.
   */
  Eigen::ArrayXd lossModulus(
    const Eigen::Ref<const Eigen::ArrayXd>& omega) const;

  /**
   * @brief Evaluate both storage and loss moduli $(G'(\omega), G''(\omega))$
   *        in a single pass over the spectrum.
   */
  std::pair<Eigen::ArrayXd, Eigen::ArrayXd> moduli(
    const Eigen::Ref<const Eigen::ArrayXd>& omega) const;

#ifdef CEREALIZABLE
  static NormalModeSpectrum fromString(std::string in)
  {
    NormalModeSpectrum s;
    pylimer_tools::utils::deserializeFromString(s, in);
    return s;
  }
#endif

private:
  Eigen::VectorXd eigenvalues_;
  Eigen::VectorXd weights_;
  bool hasEigenvectors_ = false;
  Eigen::MatrixXd eigenvectors_;
  bool isComplete_ = true;

  std::pair<Eigen::ArrayXd, Eigen::ArrayXd> evaluateModuli(
    const Eigen::Ref<const Eigen::ArrayXd>& omega,
    bool computeStorage,
    bool computeLoss) const;

  NormalModeSpectrum() = default;
#ifdef CEREALIZABLE
  friend class cereal::access;
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(eigenvalues_, weights_, hasEigenvectors_, eigenvectors_, isComplete_);
  }
#endif
};

/**
 * @brief Stateless graph Laplacian builder and normal mode solver entry point.
 *
 * Holds the assembled connectivity (Laplacian) matrix and dispatches exact or
 * approximate spectral solvers that return immutable NormalModeSpectrum
 * objects.
 */
class NormalModeAnalyzer
{
public:
  /**
   * @brief Construct from spring/bond edge lists and optional spring weights.
   *
   * @param springFrom Starting node indices for springs/bonds
   * @param springTo Ending node indices for springs/bonds
   * @param nrOfNodes Optional total number of nodes (to include isolated
   * vertices)
   * @param weights Optional per-spring weights (defaults to unit weight 1.0)
   */
  NormalModeAnalyzer(
    const std::vector<size_t>& springFrom,
    const std::vector<size_t>& springTo,
    size_t nrOfNodes = 0,
    const std::optional<std::vector<double>>& weights = std::nullopt);

  /**
   * @brief Construct directly from a square sparse Laplacian / connectivity
   * matrix.
   */
  static NormalModeAnalyzer fromMatrix(
    const Eigen::SparseMatrix<double>& matrix);

  /**
   * @brief Construct from a Universe's atom-bond graph.
   */
  static NormalModeAnalyzer fromUniverse(
    const pylimer_tools::entities::Universe& universe);

  size_t getSize() const
  {
    return static_cast<size_t>(assembledConnectivityMatrix_.rows());
  }

  const Eigen::SparseMatrix<double>& getMatrix() const
  {
    return assembledConnectivityMatrix_;
  }

  /**
   * @brief Count connected components (soluble clusters / zero modes) directly
   *        from the connectivity matrix.
   */
  size_t getNumSolubleClusters() const;

  /**
   * @brief Compute exact eigenvalues (and optionally eigenvectors) and return
   *        an immutable NormalModeSpectrum.
   *
   * @param solver Exact solver algorithm (AUTO, DENSE_LAPACK, DENSE_EIGEN,
   *               SPARSE_SLICING, METAL, or CUDA). Ignored when numModes <
   *               matrix size.
   * @param computeEigenvectors Whether to compute and store eigenvectors
   * @param decomposeComponents Whether to decompose into connected components.
   *                            Ignored when numModes < matrix size.
   * @param numModes Optional number of smallest eigenvalues to compute via
   *                 Spectra sparse solver (if omitted or equal to matrix size,
   *                 computes full spectrum)
   */
  NormalModeSpectrum solveExact(
    NormalModeSolver solver = NormalModeSolver::AUTO,
    bool computeEigenvectors = false,
    bool decomposeComponents = true,
    std::optional<size_t> numModes = std::nullopt) const;

  /**
   * @brief Compute an approximate spectrum using automatic method selection
   *        (Hybrid for n <= 30000, SLQ otherwise).
   */
  NormalModeSpectrum solveApprox() const;

  /**
   * @brief Compute an approximate spectrum using Stochastic Lanczos Quadrature.
   */
  NormalModeSpectrum solveApprox(const SLQConfig& config) const;

  /**
   * @brief Compute an approximate spectrum using Hybrid exact low-modes + SLQ.
   */
  NormalModeSpectrum solveApprox(const HybridConfig& config) const;

  /**
   * @brief Compute an approximate spectrum using the Kernel Polynomial Method.
   */
  NormalModeSpectrum solveApprox(const KPMConfig& config) const;

  /**
   * @brief Compute Jackson-damped Chebyshev moments using KPM.
   */
  Eigen::VectorXd chebyshevMoments(const KPMConfig& config = KPMConfig{}) const;

  /**
   * @brief Compute continuous density of states (DOS) on lambdaGrid using KPM.
   */
  Eigen::ArrayXd densityOfStates(
    const Eigen::Ref<const Eigen::ArrayXd>& lambdaGrid,
    const KPMConfig& config = KPMConfig{}) const;

#ifdef CEREALIZABLE
  static NormalModeAnalyzer fromString(std::string in)
  {
    NormalModeAnalyzer n;
    pylimer_tools::utils::deserializeFromString(n, in);
    return n;
  }
#endif

private:
  Eigen::SparseMatrix<double> assembledConnectivityMatrix_;

  static void solveMatrix(const Eigen::SparseMatrix<double>& mat,
                          NormalModeSolver solver,
                          bool includeEigenvectors,
                          Eigen::VectorXd& outEvals,
                          Eigen::MatrixXd& outEvecs);

  static void computeAllEigenvaluesDense(const Eigen::SparseMatrix<double>& mat,
                                         bool includeEigenvectors,
                                         bool forceEigenSolver,
                                         Eigen::VectorXd& outEvals,
                                         Eigen::MatrixXd& outEvecs);

  static void computeAllEigenvaluesSparseSlicing(
    const Eigen::SparseMatrix<double>& mat,
    bool includeEigenvectors,
    Eigen::VectorXd& outEvals,
    Eigen::MatrixXd& outEvecs);

  static void computeAllEigenvaluesCUDA(const Eigen::SparseMatrix<double>& mat,
                                        bool includeEigenvectors,
                                        Eigen::VectorXd& outEvals,
                                        Eigen::MatrixXd& outEvecs);

  void findSparseEigenvaluesImpl(size_t nrOfEigenvalues,
                                 bool includeEigenvectors,
                                 Eigen::VectorXd& outEvals,
                                 Eigen::MatrixXd& outEvecs) const;

  std::pair<double, double> getSpectralScalingParameters(
    double eps = 0.01) const;

  Eigen::VectorXd computeChebyshevMomentsCPU(size_t numMoments,
                                             size_t numRandomVectors) const;

  void computeConnectedComponents(size_t& numComponents,
                                  std::vector<size_t>& labels,
                                  std::vector<size_t>& counts) const;

  Eigen::SparseMatrix<double> extractComponentSubmatrix(
    const std::vector<size_t>& compNodes,
    std::vector<Eigen::Index>& localIdxBuffer,
    std::vector<size_t>* outLocalDeg = nullptr,
    size_t* outNumInternalEdges = nullptr,
    bool* outAllUnitWeights = nullptr) const;

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

  static Eigen::VectorXd computeJacksonKernel(size_t numMoments);

  NormalModeAnalyzer() = default;
#ifdef CEREALIZABLE
  friend class cereal::access;
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(assembledConnectivityMatrix_);
  }
#endif
};

} // namespace pylimer_tools::calc

#endif
