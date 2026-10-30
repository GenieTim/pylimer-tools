#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

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
 * If b.size() == A.rows() * 3, automatically invokes the fused Multi-RHS 3D Block-SpMV
 * solver treating A as an N x N scalar network Laplacian.
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
 * @param projectZeroModes Whether to periodically project out center-of-mass
 * translational zero modes
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
  const Eigen::VectorXd& initialX = Eigen::VectorXd(),
  bool projectZeroModes = false);

/**
 * @brief Multi-RHS 3D Block Barzilai-Borwein solver using an N x N scalar Laplacian
 * and 3 RHS components (x, y, z) per node. Fused 3-RHS SpMM achieves 3x matrix memory
 * bandwidth reduction.
 */
Eigen::VectorXd
solveSLEMetalBarzilaiBorwein3D(
  const Eigen::SparseMatrix<double>& scalarA,
  const Eigen::VectorXd& b,
  double learningRate,
  double tolerance,
  int maxIterations,
  int& iterationsDone,
  bool shortAlpha = false,
  const Eigen::VectorXd& initialX = Eigen::VectorXd(),
  bool projectZeroModes = false);

/**
 * @brief Solves the sparse linear equation system Ax = b on Apple Silicon GPU
 * using the Conjugate Gradient (CG) method implemented with Metal compute
 * kernels and optional GPU-resident Jacobi diagonal preconditioning.
 *
 * If b.size() == A.rows() * 3, automatically invokes the fused Multi-RHS 3D Block-SpMV
 * solver treating A as an N x N scalar network Laplacian.
 *
 * @param A Symmetric positive (semi-)definite coefficient matrix
 * @param b Right-hand side vector
 * @param tolerance Relative residual tolerance
 * @param maxIterations Maximum number of iterations
 * @param iterationsDone Output count of iterations performed
 * @param initialX Optional initial guess vector
 * @param projectZeroModes Whether to periodically project out center-of-mass
 * translational zero modes
 * @param useJacobi Whether to apply GPU-resident Jacobi (diagonal scaling) preconditioning
 * @return Solution vector x
 */
Eigen::VectorXd
solveSLEMetalCG(const Eigen::SparseMatrix<double>& A,
                const Eigen::VectorXd& b,
                double tolerance,
                int maxIterations,
                int& iterationsDone,
                const Eigen::VectorXd& initialX = Eigen::VectorXd(),
                bool projectZeroModes = false,
                bool useJacobi = true);

/**
 * @brief Multi-RHS 3D Block Conjugate Gradient solver using an N x N scalar Laplacian
 * and 3 RHS components (x, y, z) per node with GPU-resident Jacobi preconditioning.
 * Fused 3-RHS SpMM achieves 3x matrix memory bandwidth reduction.
 */
Eigen::VectorXd
solveSLEMetalCG3D(const Eigen::SparseMatrix<double>& scalarA,
                  const Eigen::VectorXd& b,
                  double tolerance,
                  int maxIterations,
                  int& iterationsDone,
                  const Eigen::VectorXd& initialX = Eigen::VectorXd(),
                  bool projectZeroModes = false,
                  bool useJacobi = true);

/**
 * @brief Multiplies an N x N scalar matrix with a 3N dense vector (interpreted as N float3 vectors)
 * on GPU using the fused 3-RHS SpMM kernel.
 */
Eigen::VectorXd
spmvMetal3D(const Eigen::SparseMatrix<double>& scalarA, const Eigen::VectorXd& x);

/**
 * @brief Multiplies a sparse matrix with a dense vector on GPU.
 */
Eigen::VectorXd
spmvMetal(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& x);

/**
 * @brief Benchmark pure SpMV execution time (in milliseconds) for 3D Block-SpMV vs Full SpMV
 * over multiple repetitions directly on GPU.
 */
struct SpmvBenchmarkResult
{
  double time3D_ms = 0.0;
  double timeFull_ms = 0.0;
  double speedup = 0.0;
  size_t memory3D_bytes = 0;
  size_t memoryFull_bytes = 0;
  double bandwidthReductionRatio = 0.0;
};

SpmvBenchmarkResult
benchmarkSpMVMetal(const Eigen::SparseMatrix<double>& scalarA,
                   const Eigen::SparseMatrix<double>& fullA,
                   int numRepetitions = 100);

/**
 * @brief Data structure for GPU spring strain and inactivity evaluation.
 */
struct SpringRelaxationData
{
  const int* springIndexA = nullptr;
  const int* springIndexB = nullptr;
  const double* springBoxOffset = nullptr;     // 3 * numSprings doubles
  const double* springContourLength = nullptr;  // numSprings doubles
  size_t numSprings = 0;
  double breakingThreshold = 0.0;
  double inactiveCutoff = 0.0;
};

/**
 * @brief Evaluation result for GPU spring strain and inactivity check.
 */
struct SpringEvaluationResult
{
  uint32_t numBroken = 0;
  uint32_t numInactive = 0;
  uint32_t totalFlagged = 0;
  std::vector<uint32_t> springFlags; // bit 0 = broken, bit 1 = inactive
};

/**
 * @brief Evaluates spring strains and inactivity on GPU using the current GPU-resident
 * coordinates, avoiding CPU host-device synchronization when no springs are flagged.
 */
SpringEvaluationResult
evaluateSpringStrainAndInactivityMetal(
  size_t numNodes,
  const SpringRelaxationData& springData);

/**
 * @brief Retrieves the current GPU-resident coordinates from the Metal buffer pool.
 */
Eigen::VectorXd
getGpuResidentCoordinates(size_t numCoords);

} // namespace pylimer_tools::sim::mehp
