#include "./NormalModeAnalyzer.h"

#include "../entities/Universe.h"
#include "../utils/ExtraEigenTypes.h"
#include "../utils/GraphLaplacianUtils.h"
#include "../utils/utilityMacros.h"
#ifdef OPENMP_FOUND
#include <omp.h>
#endif
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <algorithm>
#include <cassert>
#include <cmath>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>
#include <vector>
#ifdef EIGEN_USE_LAPACKE
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#endif

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include <cusolverDn.h>
#endif

#ifdef USE_METAL
#include "MetalNormalModeSolver.h"
#endif

#define NORMAL_MODE_ZERO_TOLERANCE 1e-12

namespace pylimer_tools::calc {

// ============================================================================
// MARK: NormalModeSpectrum
// ============================================================================

NormalModeSpectrum::NormalModeSpectrum(
  Eigen::VectorXd eigenvalues,
  std::optional<Eigen::VectorXd> weights,
  std::optional<Eigen::MatrixXd> eigenvectors,
  const bool isComplete)
  : eigenvalues_(std::move(eigenvalues))
  , isComplete_(isComplete)
{
  for (Eigen::Index i = 0; i < this->eigenvalues_.size(); ++i) {
    if (std::abs(this->eigenvalues_[i]) < NORMAL_MODE_ZERO_TOLERANCE) {
      this->eigenvalues_[i] = 0.0;
    }
  }

  if (weights.has_value()) {
    INVALIDARG_EXP_IFN(weights->size() == this->eigenvalues_.size(),
                       "Weights size must match eigenvalues size.");
    this->weights_ = std::move(*weights);
  } else {
    this->weights_ = Eigen::VectorXd::Ones(this->eigenvalues_.size());
  }

  if (eigenvectors.has_value()) {
    INVALIDARG_EXP_IFN(
      eigenvectors->cols() == this->eigenvalues_.size(),
      "Eigenvector matrix column count must match eigenvalues size.");
    this->eigenvectors_ = std::move(*eigenvectors);
    this->hasEigenvectors_ = true;
  }
}

const Eigen::VectorXd&
NormalModeSpectrum::getEigenvalues() const
{
  return this->eigenvalues_;
}

const Eigen::VectorXd&
NormalModeSpectrum::getWeights() const
{
  return this->weights_;
}

std::optional<Eigen::MatrixXd>
NormalModeSpectrum::getEigenvectors() const
{
  if (!this->hasEigenvectors_) {
    return std::nullopt;
  }
  return this->eigenvectors_;
}

bool
NormalModeSpectrum::isComplete() const
{
  return this->isComplete_;
}

size_t
NormalModeSpectrum::getNumModes() const
{
  return static_cast<size_t>(this->eigenvalues_.size());
}

size_t
NormalModeSpectrum::getNumZeroModes() const
{
  return static_cast<size_t>((this->eigenvalues_.array() == 0.0).count());
}

Eigen::ArrayXd
NormalModeSpectrum::stressAutocorrelation(
  const Eigen::Ref<const Eigen::ArrayXd>& t) const
{
  Eigen::ArrayXd result = Eigen::ArrayXd::Zero(t.size());
  for (Eigen::Index i = 0; i < this->eigenvalues_.size(); ++i) {
    const double lambda = this->eigenvalues_[i];
    if (lambda <= 0.0 || !std::isfinite(lambda)) {
      continue;
    }
    const double w = this->weights_[i];
    if (w == 0.0 || !std::isfinite(w)) {
      continue;
    }
    result += w * (-2.0 * lambda * t).exp();
  }
  RUNTIME_EXP_IFN(result.allFinite(),
                  "Stress autocorrelation is not fully finite.");
  return result;
}

std::pair<Eigen::ArrayXd, Eigen::ArrayXd>
NormalModeSpectrum::evaluateModuli(
  const Eigen::Ref<const Eigen::ArrayXd>& omega,
  const bool computeStorage,
  const bool computeLoss) const
{
  const Eigen::Index nFreq = omega.size();
  Eigen::ArrayXd storage =
    computeStorage ? Eigen::ArrayXd::Zero(nFreq) : Eigen::ArrayXd();
  Eigen::ArrayXd loss =
    computeLoss ? Eigen::ArrayXd::Zero(nFreq) : Eigen::ArrayXd();

  Eigen::ArrayXd r(nFreq);
  Eigen::ArrayXd r2(nFreq);
  Eigen::ArrayXd invDenom(nFreq);

  for (Eigen::Index i = 0; i < this->eigenvalues_.size(); ++i) {
    const double lambda = this->eigenvalues_[i];
    if (lambda <= 0.0 || !std::isfinite(lambda)) {
      continue;
    }
    const double w = this->weights_[i];
    if (w == 0.0 || !std::isfinite(w)) {
      continue;
    }
    const double invTwoLambda = 0.5 / lambda;
    r = omega * invTwoLambda;
    r2 = r.square();
    invDenom = w / (1.0 + r2);
    if (computeStorage) {
      storage += r2 * invDenom;
    }
    if (computeLoss) {
      loss += r * invDenom;
    }
  }

  if (computeStorage) {
    RUNTIME_EXP_IFN(storage.allFinite(),
                    "Storage modulus is not fully finite.");
  }
  if (computeLoss) {
    RUNTIME_EXP_IFN(loss.allFinite(), "Loss modulus is not fully finite.");
  }
  return { std::move(storage), std::move(loss) };
}

Eigen::ArrayXd
NormalModeSpectrum::storageModulus(
  const Eigen::Ref<const Eigen::ArrayXd>& omega) const
{
  return this->evaluateModuli(omega, true, false).first;
}

Eigen::ArrayXd
NormalModeSpectrum::lossModulus(
  const Eigen::Ref<const Eigen::ArrayXd>& omega) const
{
  return this->evaluateModuli(omega, false, true).second;
}

std::pair<Eigen::ArrayXd, Eigen::ArrayXd>
NormalModeSpectrum::moduli(const Eigen::Ref<const Eigen::ArrayXd>& omega) const
{
  return this->evaluateModuli(omega, true, true);
}

// ============================================================================
// MARK: NormalModeAnalyzer Construction
// ============================================================================

NormalModeAnalyzer::NormalModeAnalyzer(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  const size_t nrOfNodes,
  const std::optional<std::vector<double>>& weights)
{
  if (weights.has_value()) {
    INVALIDARG_EXP_IFN(weights->size() == springFrom.size(),
                       "weights must have the same size as springFrom and "
                       "springTo");
  }
  static const std::vector<double> kEmptyWeights;
  this->assembledConnectivityMatrix_ =
    pylimer_tools::utils::assembleGraphLaplacian(
      springFrom,
      springTo,
      weights.has_value() ? *weights : kEmptyWeights,
      nrOfNodes);
}

NormalModeAnalyzer
NormalModeAnalyzer::fromMatrix(const Eigen::SparseMatrix<double>& matrix)
{
  INVALIDARG_EXP_IFN(matrix.rows() == matrix.cols(),
                     "Connectivity matrix must be square.");
  NormalModeAnalyzer analyzer;
  analyzer.assembledConnectivityMatrix_ = matrix;
  analyzer.assembledConnectivityMatrix_.makeCompressed();
  return analyzer;
}

NormalModeAnalyzer
NormalModeAnalyzer::fromUniverse(
  const pylimer_tools::entities::Universe& universe)
{
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  size_t nAtoms = 0;
  pylimer_tools::utils::extractUniverseEdges(
    universe, springFrom, springTo, nAtoms);
  return NormalModeAnalyzer(springFrom, springTo, nAtoms);
}

size_t
NormalModeAnalyzer::getNumSolubleClusters() const
{
  size_t numComponents = 0;
  std::vector<size_t> labels;
  std::vector<size_t> counts;
  this->computeConnectedComponents(numComponents, labels, counts);
  return numComponents;
}

// ============================================================================
// MARK: Exact Solvers
// ============================================================================

void
NormalModeAnalyzer::findSparseEigenvaluesImpl(const size_t nrOfEigenvalues,
                                              const bool includeEigenvectors,
                                              Eigen::VectorXd& outEvals,
                                              Eigen::MatrixXd& outEvecs) const
{
  INVALIDARG_EXP_IFN(nrOfEigenvalues > 0,
                     "nrOfEigenvalues must be greater than 0");
  const size_t n =
    static_cast<size_t>(this->assembledConnectivityMatrix_.rows());
  INVALIDARG_EXP_IFN(nrOfEigenvalues < n,
                     "nrOfEigenvalues must be less than matrix size");

  auto sortResults = [&](const Eigen::VectorXd& rawEvals,
                         const Eigen::MatrixXd& rawEvecs) {
    std::vector<Eigen::Index> order(static_cast<size_t>(rawEvals.size()));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(),
              order.end(),
              [&rawEvals](const Eigen::Index a, const Eigen::Index b) {
                return rawEvals[a] < rawEvals[b];
              });
    outEvals.resize(rawEvals.size());
    for (size_t i = 0; i < order.size(); ++i) {
      outEvals[static_cast<Eigen::Index>(i)] = rawEvals[order[i]];
    }
    if (includeEigenvectors) {
      outEvecs.resize(rawEvecs.rows(), rawEvecs.cols());
      for (size_t i = 0; i < order.size(); ++i) {
        outEvecs.col(static_cast<Eigen::Index>(i)) = rawEvecs.col(order[i]);
      }
    }
  };

  try {
    Spectra::SparseSymShiftSolve<double> op(this->assembledConnectivityMatrix_);
    size_t ncv = std::min(2 * nrOfEigenvalues + 1, n);
    if (ncv <= nrOfEigenvalues) {
      ncv = std::min(nrOfEigenvalues + 2, n);
    }
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(
      op,
      static_cast<Eigen::Index>(nrOfEigenvalues),
      static_cast<Eigen::Index>(ncv),
      -1e-5);
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestMagn);
    if (eigs.info() == Spectra::CompInfo::Successful) {
      sortResults(eigs.eigenvalues(),
                  includeEigenvectors ? eigs.eigenvectors()
                                      : Eigen::MatrixXd());
      return;
    }
  } catch (...) {
    // Fallback to standard solver below
  }

  Spectra::SparseGenMatProd<double> op(this->assembledConnectivityMatrix_);
  size_t ncv = std::min(2 * nrOfEigenvalues, n);
  if (ncv <= nrOfEigenvalues) {
    ncv = std::min(nrOfEigenvalues + 2, n);
  }
  Spectra::SymEigsSolver<Spectra::SparseGenMatProd<double>> eigs(
    op,
    static_cast<Eigen::Index>(nrOfEigenvalues),
    static_cast<Eigen::Index>(ncv));
  eigs.init();
  eigs.compute(Spectra::SortRule::SmallestAlge);

  RUNTIME_EXP_IFN(eigs.info() == Spectra::CompInfo::Successful,
                  "Spectra failed to compute sparse eigenvalues");
  sortResults(eigs.eigenvalues(),
              includeEigenvectors ? eigs.eigenvectors() : Eigen::MatrixXd());
}

void
NormalModeAnalyzer::computeAllEigenvaluesDense(
  const Eigen::SparseMatrix<double>& mat,
  const bool includeEigenvectors,
  const bool forceEigenSolver,
  Eigen::VectorXd& outEvals,
  Eigen::MatrixXd& outEvecs)
{
  const Eigen::Index nRows = mat.rows();
  if (nRows == 0) {
    outEvals = Eigen::VectorXd::Zero(0);
    outEvecs = Eigen::MatrixXd::Zero(0, 0);
    return;
  }

  Eigen::MatrixXd matDense = Eigen::MatrixXd(mat);

#ifndef EIGEN_USE_LAPACKE
  (void)forceEigenSolver;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(
    matDense,
    includeEigenvectors ? Eigen::DecompositionOptions::ComputeEigenvectors
                        : Eigen::DecompositionOptions::EigenvaluesOnly);
  outEvals = solver.eigenvalues();
  if (includeEigenvectors) {
    outEvecs = solver.eigenvectors();
  }
#else
  if (forceEigenSolver) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(
      matDense,
      includeEigenvectors ? Eigen::DecompositionOptions::ComputeEigenvectors
                          : Eigen::DecompositionOptions::EigenvaluesOnly);
    outEvals = solver.eigenvalues();
    if (includeEigenvectors) {
      outEvecs = solver.eigenvectors();
    }
    return;
  }

  Eigen::VectorXd eigenvaluesMemory = Eigen::VectorXd::Zero(nRows);
  const lapack_int size = static_cast<lapack_int>(nRows);

  if (includeEigenvectors) {
    const int info = LAPACKE_dsyevd(
      Eigen::MatrixXd::IsRowMajor ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR,
      'V',
      'U',
      size,
      matDense.data(),
      size,
      eigenvaluesMemory.data());
    RUNTIME_EXP_IFN(info == 0,
                    "Error in LAPACKE_dsyevd. Exit code " +
                      std::to_string(info) + ".");
    outEvals = eigenvaluesMemory;
    outEvecs = matDense;
  } else {
    Eigen::MatrixXd dummyZ = Eigen::MatrixXd::Zero(nRows, 1);
    const lapack_int ldz = 1;
    const lapack_int il = 0;
    const lapack_int iu = 0;
    const double abstol = 1e-10;
    lapack_int M = 0;
    Eigen::Vector<lapack_int, Eigen::Dynamic> support =
      Eigen::Vector<lapack_int, Eigen::Dynamic>::Zero(nRows * 2);

    const int info = LAPACKE_dsyevr(
      Eigen::MatrixXd::IsRowMajor ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR,
      'N',
      'A',
      'U',
      size,
      matDense.data(),
      size,
      0,
      0,
      il,
      iu,
      abstol,
      &M,
      eigenvaluesMemory.data(),
      dummyZ.data(),
      ldz,
      support.data());
    RUNTIME_EXP_IFN(info == 0,
                    "Error in LAPACKE_dsyevr. Exit code " +
                      std::to_string(info) + ".");
    outEvals = eigenvaluesMemory;
  }
#endif
}

void
NormalModeAnalyzer::computeAllEigenvaluesSparseSlicing(
  const Eigen::SparseMatrix<double>& mat,
  const bool includeEigenvectors,
  Eigen::VectorXd& outEvals,
  Eigen::MatrixXd& outEvecs)
{
  const Eigen::Index n = mat.rows();
  if (n == 0) {
    outEvals = Eigen::VectorXd::Zero(0);
    outEvecs = Eigen::MatrixXd::Zero(0, 0);
    return;
  }
  if (n == 1) {
    outEvals = Eigen::VectorXd::Zero(1);
    if (includeEigenvectors) {
      outEvecs = Eigen::MatrixXd::Identity(1, 1);
    }
    return;
  }

  if (n <= 30) {
    computeAllEigenvaluesDense(
      mat, includeEigenvectors, false, outEvals, outEvecs);
    return;
  }

  const double maxDiag = mat.diagonal().maxCoeff();
  const double bUpper = 2.0 * maxDiag + 0.1;
  const double aLower = -1e-7;

  Eigen::SparseMatrix<double> shifted = mat;
  for (Eigen::Index i = 0; i < n; ++i) {
    shifted.coeffRef(i, i) += 0.0;
  }
  shifted.makeCompressed();
  const Eigen::VectorXd origDiag = shifted.diagonal();
  std::vector<double*> diagValPtrs(static_cast<size_t>(n), nullptr);
  for (Eigen::Index col = 0; col < shifted.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(shifted, col); it;
         ++it) {
      if (it.row() == col) {
        diagValPtrs[static_cast<size_t>(col)] = &it.valueRef();
        break;
      }
    }
  }
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.analyzePattern(shifted);

  auto countLessThan = [&](const double sigma) -> size_t {
    const double perturbations[] = {
      0.0, 1.234567e-8, -2.345678e-8, 5.678901e-7
    };
    for (const double pert : perturbations) {
      const double currentSigma = sigma + pert;
      for (Eigen::Index i = 0; i < n; ++i) {
        *diagValPtrs[static_cast<size_t>(i)] = origDiag[i] - currentSigma;
      }

      ldlt.factorize(shifted);
      if (ldlt.info() == Eigen::Success) {
        const auto& D = ldlt.vectorD();
        return static_cast<size_t>((D.array() < 0.0).count());
      }
    }
    RUNTIME_EXP_IFN(
      false,
      "SimplicialLDLT factorization failed in countLessThan for shift " +
        std::to_string(sigma));
    return 0;
  };

  struct Slice
  {
    double a;
    double b;
    size_t countA;
    size_t countB;
  };

  const size_t totalCount = static_cast<size_t>(n);
  std::vector<Slice> stack;
  stack.push_back({ aLower, bUpper, 0, totalCount });

  std::vector<double> eigenvaluesList;
  eigenvaluesList.reserve(totalCount);

  const size_t targetSliceSize =
    std::min(static_cast<size_t>(50), static_cast<size_t>(n));

  while (!stack.empty()) {
    const Slice slice = stack.back();
    stack.pop_back();

    const size_t k = slice.countB - slice.countA;
    if (k == 0) {
      continue;
    }

    if ((slice.b - slice.a) < 1e-6) {
      const double mid = (slice.a + slice.b) / 2.0;
      for (size_t i = 0; i < k; ++i) {
        eigenvaluesList.push_back(mid);
      }
      continue;
    }

    if (k <= 3) {
      if (k == 1) {
        double low = slice.a;
        double high = slice.b;
        while ((high - low) > 1e-8) {
          const double mid = (low + high) / 2.0;
          const size_t c = countLessThan(mid);
          if (c <= slice.countA) {
            low = mid;
          } else {
            high = mid;
          }
        }
        eigenvaluesList.push_back((low + high) / 2.0);
      } else {
        const double mid = (slice.a + slice.b) / 2.0;
        const size_t countMid = countLessThan(mid);
        stack.push_back({ mid, slice.b, countMid, slice.countB });
        stack.push_back({ slice.a, mid, slice.countA, countMid });
      }
      continue;
    }

    if (k <= targetSliceSize) {
      bool solvedBySpectra = false;
      try {
        const double sigma = (slice.a + slice.b) / 2.0;
        Spectra::SparseSymShiftSolve<double> op(mat);
        const size_t nev = k;
        size_t ncv = std::min(2 * nev + 2, static_cast<size_t>(n));
        if (ncv <= nev) {
          ncv = std::min(nev + 2, static_cast<size_t>(n));
        }
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(
          op,
          static_cast<Eigen::Index>(nev),
          static_cast<Eigen::Index>(ncv),
          sigma);
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn);
        if (eigs.info() == Spectra::CompInfo::Successful) {
          const Eigen::VectorXd evals = eigs.eigenvalues();
          std::vector<double> inSlice;
          for (Eigen::Index i = 0; i < evals.size(); ++i) {
            if (evals[i] >= slice.a && evals[i] < slice.b) {
              inSlice.push_back(evals[i]);
            }
          }
          if (inSlice.size() == k) {
            for (const double val : inSlice) {
              eigenvaluesList.push_back(val);
            }
            solvedBySpectra = true;
          }
        }
      } catch (...) {
        solvedBySpectra = false;
      }

      if (solvedBySpectra) {
        continue;
      }
    }

    const double mid = (slice.a + slice.b) / 2.0;
    const size_t countMid = countLessThan(mid);
    stack.push_back({ mid, slice.b, countMid, slice.countB });
    stack.push_back({ slice.a, mid, slice.countA, countMid });
  }

  std::sort(eigenvaluesList.begin(), eigenvaluesList.end());
  outEvals = Eigen::Map<Eigen::VectorXd>(
    eigenvaluesList.data(), static_cast<Eigen::Index>(eigenvaluesList.size()));

  if (includeEigenvectors) {
    Eigen::VectorXd dummyEvals;
    computeAllEigenvaluesDense(mat, true, false, dummyEvals, outEvecs);
  }
}

void
NormalModeAnalyzer::computeAllEigenvaluesCUDA(
  const Eigen::SparseMatrix<double>& mat,
  const bool includeEigenvectors,
  Eigen::VectorXd& outEvals,
  Eigen::MatrixXd& outEvecs)
{
  const Eigen::Index nRows = mat.rows();
  if (nRows == 0) {
    outEvals = Eigen::VectorXd::Zero(0);
    if (includeEigenvectors) {
      outEvecs = Eigen::MatrixXd::Zero(0, 0);
    }
    return;
  }

#ifdef USE_CUDA
  Eigen::MatrixXd matDense = Eigen::MatrixXd(mat);

  cusolverDnHandle_t cusolverH = nullptr;
  double* d_A = nullptr;
  double* d_W = nullptr;
  int* d_info = nullptr;
  double* d_work = nullptr;

  struct CudaCleanup
  {
    cusolverDnHandle_t& h;
    double*& A;
    double*& W;
    int*& info;
    double*& work;
    ~CudaCleanup()
    {
      if (work)
        cudaFree(work);
      if (info)
        cudaFree(info);
      if (W)
        cudaFree(W);
      if (A)
        cudaFree(A);
      if (h)
        cusolverDnDestroy(h);
    }
  } cleanup{ cusolverH, d_A, d_W, d_info, d_work };

  cusolverStatus_t cusolverStatus = cusolverDnCreate(&cusolverH);
  RUNTIME_EXP_IFN(cusolverStatus == CUSOLVER_STATUS_SUCCESS,
                  "cuSOLVER initialization failed");

  cudaError_t cudaStat =
    cudaMalloc((void**)&d_A, sizeof(double) * nRows * nRows);
  RUNTIME_EXP_IFN(cudaStat == cudaSuccess, "cudaMalloc failed for d_A");
  cudaStat = cudaMalloc((void**)&d_W, sizeof(double) * nRows);
  RUNTIME_EXP_IFN(cudaStat == cudaSuccess, "cudaMalloc failed for d_W");
  cudaStat = cudaMalloc((void**)&d_info, sizeof(int));
  RUNTIME_EXP_IFN(cudaStat == cudaSuccess, "cudaMalloc failed for d_info");

  cudaStat = cudaMemcpy(d_A,
                        matDense.data(),
                        sizeof(double) * nRows * nRows,
                        cudaMemcpyHostToDevice);
  RUNTIME_EXP_IFN(cudaStat == cudaSuccess, "cudaMemcpy failed to device");

  const cusolverEigMode_t jobz =
    includeEigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR;
  int lwork = 0;
  cusolverStatus = cusolverDnDsyevd_bufferSize(
    cusolverH, jobz, CUBLAS_FILL_MODE_LOWER, nRows, d_A, nRows, d_W, &lwork);
  RUNTIME_EXP_IFN(cusolverStatus == CUSOLVER_STATUS_SUCCESS,
                  "cuSOLVER buffer size calculation failed");

  cudaStat = cudaMalloc((void**)&d_work, sizeof(double) * lwork);
  RUNTIME_EXP_IFN(cudaStat == cudaSuccess, "cudaMalloc failed for d_work");

  cusolverStatus = cusolverDnDsyevd(cusolverH,
                                    jobz,
                                    CUBLAS_FILL_MODE_LOWER,
                                    nRows,
                                    d_A,
                                    nRows,
                                    d_W,
                                    d_work,
                                    lwork,
                                    d_info);
  RUNTIME_EXP_IFN(cusolverStatus == CUSOLVER_STATUS_SUCCESS,
                  "cuSOLVER dsyevd failed");

  cudaStat = cudaDeviceSynchronize();
  RUNTIME_EXP_IFN(cudaStat == cudaSuccess, "cudaDeviceSynchronize failed");

  outEvals = Eigen::VectorXd::Zero(nRows);
  cudaStat = cudaMemcpy(
    outEvals.data(), d_W, sizeof(double) * nRows, cudaMemcpyDeviceToHost);
  RUNTIME_EXP_IFN(cudaStat == cudaSuccess,
                  "cudaMemcpy failed to host (eigenvalues)");

  if (includeEigenvectors) {
    outEvecs = Eigen::MatrixXd::Zero(nRows, nRows);
    cudaStat = cudaMemcpy(outEvecs.data(),
                          d_A,
                          sizeof(double) * nRows * nRows,
                          cudaMemcpyDeviceToHost);
    RUNTIME_EXP_IFN(cudaStat == cudaSuccess,
                    "cudaMemcpy failed to host (eigenvectors)");
  }

  int info = 0;
  cudaStat = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
  RUNTIME_EXP_IFN(
    info == 0, "cuSOLVER dsyevd info return non-zero: " + std::to_string(info));
#else
  computeAllEigenvaluesDense(
    mat, includeEigenvectors, false, outEvals, outEvecs);
#endif
}

void
NormalModeAnalyzer::solveMatrix(const Eigen::SparseMatrix<double>& mat,
                                const NormalModeSolver solver,
                                const bool includeEigenvectors,
                                Eigen::VectorXd& outEvals,
                                Eigen::MatrixXd& outEvecs)
{
  const Eigen::Index n = mat.rows();
  if (n == 0) {
    outEvals = Eigen::VectorXd::Zero(0);
    if (includeEigenvectors) {
      outEvecs = Eigen::MatrixXd::Zero(0, 0);
    }
    return;
  }
  if (solver == NormalModeSolver::METAL) {
#ifdef USE_METAL
    if (isMetalAvailable()) {
      computeAllEigenvaluesMetal(mat, includeEigenvectors, outEvals, outEvecs);
    } else {
      computeAllEigenvaluesDense(
        mat, includeEigenvectors, false, outEvals, outEvecs);
    }
#else
    computeAllEigenvaluesDense(
      mat, includeEigenvectors, false, outEvals, outEvecs);
#endif
  } else if (solver == NormalModeSolver::CUDA) {
    computeAllEigenvaluesCUDA(mat, includeEigenvectors, outEvals, outEvecs);
  } else if (solver == NormalModeSolver::SPARSE_SLICING ||
             (solver == NormalModeSolver::AUTO && !includeEigenvectors &&
              n > 5000)) {
    computeAllEigenvaluesSparseSlicing(
      mat, includeEigenvectors, outEvals, outEvecs);
  } else {
    computeAllEigenvaluesDense(mat,
                               includeEigenvectors,
                               (solver == NormalModeSolver::DENSE_EIGEN),
                               outEvals,
                               outEvecs);
  }
}

NormalModeSpectrum
NormalModeAnalyzer::solveExact(const NormalModeSolver solver,
                               const bool computeEigenvectors,
                               const bool decomposeComponents,
                               const std::optional<size_t> numModes) const
{
  RUNTIME_EXP_IFN(this->assembledConnectivityMatrix_.rows() ==
                    this->assembledConnectivityMatrix_.cols(),
                  "Expected square matrix");

  const Eigen::Index nRows = this->assembledConnectivityMatrix_.rows();
  if (numModes.has_value()) {
    INVALIDARG_EXP_IFN(*numModes > 0, "numModes must be greater than 0");
    INVALIDARG_EXP_IFN(*numModes <= static_cast<size_t>(nRows),
                       "numModes must not exceed matrix size");
    if (*numModes < static_cast<size_t>(nRows)) {
      Eigen::VectorXd evals;
      Eigen::MatrixXd evecs;
      this->findSparseEigenvaluesImpl(
        *numModes, computeEigenvectors, evals, evecs);
      return NormalModeSpectrum(
        std::move(evals),
        std::nullopt,
        computeEigenvectors ? std::optional<Eigen::MatrixXd>(std::move(evecs))
                            : std::nullopt,
        false);
    }
  }

  if (nRows == 0) {
    return NormalModeSpectrum(
      Eigen::VectorXd::Zero(0),
      std::nullopt,
      computeEigenvectors
        ? std::optional<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(0, 0))
        : std::nullopt,
      true);
  }

  if (!decomposeComponents) {
    Eigen::VectorXd evals;
    Eigen::MatrixXd evecs;
    solveMatrix(this->assembledConnectivityMatrix_,
                solver,
                computeEigenvectors,
                evals,
                evecs);
    return NormalModeSpectrum(
      std::move(evals),
      std::nullopt,
      computeEigenvectors ? std::optional<Eigen::MatrixXd>(std::move(evecs))
                          : std::nullopt,
      true);
  }

  // --- Component-by-Component Decomposition ---
  size_t numComponents = 0;
  std::vector<size_t> labels;
  std::vector<size_t> counts;
  this->computeConnectedComponents(numComponents, labels, counts);

  std::vector<std::vector<size_t>> componentNodes(numComponents);
  for (size_t c = 0; c < numComponents; ++c) {
    componentNodes[c].reserve(counts[c]);
  }
  for (size_t i = 0; i < static_cast<size_t>(nRows); ++i) {
    componentNodes[labels[i]].push_back(i);
  }

  std::vector<double> allEigenvalues;
  allEigenvalues.reserve(static_cast<size_t>(nRows));
  Eigen::MatrixXd finalEvecs;
  if (computeEigenvectors) {
    finalEvecs = Eigen::MatrixXd::Zero(nRows, nRows);
  }
  Eigen::Index colOffset = 0;
  std::vector<Eigen::Index> localIdxBuffer(static_cast<size_t>(nRows), -1);

  for (const std::vector<size_t>& comp : componentNodes) {
    const size_t M = comp.size();
    if (M == 1) {
      if (computeEigenvectors) {
        finalEvecs(static_cast<Eigen::Index>(comp[0]), colOffset) = 1.0;
        colOffset++;
      }
      allEigenvalues.push_back(0.0);
      continue;
    }

    std::vector<size_t> localDeg;
    size_t numInternalEdges = 0;
    bool allUnitWeights = true;
    const Eigen::SparseMatrix<double> subMat = extractComponentSubmatrix(
      comp,
      localIdxBuffer,
      computeEigenvectors ? nullptr : &localDeg,
      computeEigenvectors ? nullptr : &numInternalEdges,
      computeEigenvectors ? nullptr : &allUnitWeights);

    if (!computeEigenvectors) {
      if (M == 2 && allUnitWeights && numInternalEdges == 1) {
        allEigenvalues.push_back(0.0);
        allEigenvalues.push_back(2.0);
        continue;
      }

      size_t deg1Count = 0;
      size_t deg2Count = 0;
      for (const size_t d : localDeg) {
        if (d == 1) {
          deg1Count++;
        } else if (d == 2) {
          deg2Count++;
        }
      }

      if (allUnitWeights && numInternalEdges == M - 1 && deg1Count == 2 &&
          deg2Count == M - 2) {
        for (size_t p = 0; p < M; ++p) {
          const double val =
            4.0 * std::pow(std::sin((static_cast<double>(p) * M_PI) /
                                    (2.0 * static_cast<double>(M))),
                           2);
          allEigenvalues.push_back(val);
        }
        continue;
      }

      if (allUnitWeights && numInternalEdges == M && deg2Count == M) {
        for (size_t p = 0; p < M; ++p) {
          const double val =
            4.0 * std::pow(std::sin((static_cast<double>(p) * M_PI) /
                                    static_cast<double>(M)),
                           2);
          allEigenvalues.push_back(val);
        }
        continue;
      }
    }

    Eigen::VectorXd subEvals;
    Eigen::MatrixXd subEvecs;
    solveMatrix(subMat, solver, computeEigenvectors, subEvals, subEvecs);

    for (Eigen::Index j = 0; j < subEvals.size(); ++j) {
      allEigenvalues.push_back(subEvals[j]);
    }
    if (computeEigenvectors) {
      for (Eigen::Index j = 0; j < subEvals.size(); ++j) {
        for (Eigen::Index k = 0; k < subEvecs.rows(); ++k) {
          finalEvecs(static_cast<Eigen::Index>(comp[static_cast<size_t>(k)]),
                     colOffset + j) = subEvecs(k, j);
        }
      }
      colOffset += subEvals.size();
    }
  }

  if (computeEigenvectors) {
    std::vector<Eigen::Index> order(static_cast<size_t>(nRows));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(),
              order.end(),
              [&allEigenvalues](const Eigen::Index a, const Eigen::Index b) {
                return allEigenvalues[static_cast<size_t>(a)] <
                       allEigenvalues[static_cast<size_t>(b)];
              });
    Eigen::VectorXd finalEvals(nRows);
    for (Eigen::Index i = 0; i < nRows; ++i) {
      finalEvals[i] =
        allEigenvalues[static_cast<size_t>(order[static_cast<size_t>(i)])];
    }

    // In-place column permutation: new column i <- old column order[i]
    std::vector<bool> visited(static_cast<size_t>(nRows), false);
    Eigen::VectorXd tmpCol(nRows);
    for (Eigen::Index i = 0; i < nRows; ++i) {
      if (visited[static_cast<size_t>(i)] ||
          order[static_cast<size_t>(i)] == i) {
        visited[static_cast<size_t>(i)] = true;
        continue;
      }
      Eigen::Index cur = i;
      tmpCol = finalEvecs.col(cur);
      while (true) {
        visited[static_cast<size_t>(cur)] = true;
        const Eigen::Index nxt = order[static_cast<size_t>(cur)];
        if (nxt == i) {
          finalEvecs.col(cur) = tmpCol;
          break;
        }
        finalEvecs.col(cur) = finalEvecs.col(nxt);
        cur = nxt;
      }
    }
    return NormalModeSpectrum(
      std::move(finalEvals), std::nullopt, std::move(finalEvecs), true);
  }

  std::sort(allEigenvalues.begin(), allEigenvalues.end());
  Eigen::VectorXd finalEvals = Eigen::Map<Eigen::VectorXd>(
    allEigenvalues.data(), static_cast<Eigen::Index>(allEigenvalues.size()));
  return NormalModeSpectrum(
    std::move(finalEvals), std::nullopt, std::nullopt, true);
}

// ============================================================================
// MARK: KPM & Approximate Spectral Methods
// ============================================================================

std::pair<double, double>
NormalModeAnalyzer::getSpectralScalingParameters(const double eps) const
{
  const Eigen::Index n = this->assembledConnectivityMatrix_.rows();
  if (n == 0) {
    return { 0.0, 1.0 };
  }
  const double maxDiag =
    this->assembledConnectivityMatrix_.diagonal().maxCoeff();
  const double lambdaMax = 2.0 * maxDiag + 0.1;
  const double c = lambdaMax / 2.0;
  const double w = lambdaMax / (2.0 * (1.0 - eps));
  return { c, w };
}

Eigen::VectorXd
computeJacksonKernel(const size_t numMoments)
{
  const Eigen::Index numMomentsIdx = static_cast<Eigen::Index>(numMoments);
  Eigen::VectorXd kernel(numMomentsIdx);
  const double mD = static_cast<double>(numMoments);
  const double denom = std::tan(M_PI / (mD + 1.0));
  for (Eigen::Index m = 0; m < numMomentsIdx; ++m) {
    const double md = static_cast<double>(m);
    kernel[m] = ((mD - md + 1.0) * std::cos(M_PI * md / (mD + 1.0)) +
                 std::sin(M_PI * md / (mD + 1.0)) / denom) /
                (mD + 1.0);
  }
  return kernel;
}

Eigen::VectorXd
NormalModeAnalyzer::computeJacksonKernel(const size_t numMoments)
{
  return pylimer_tools::calc::computeJacksonKernel(numMoments);
}

Eigen::VectorXd
NormalModeAnalyzer::computeChebyshevMomentsCPU(
  const size_t numMoments,
  const size_t numRandomVectors) const
{
  const Eigen::Index n = this->assembledConnectivityMatrix_.rows();
  if (n == 0 || numMoments == 0) {
    return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numMoments));
  }

  const std::pair<double, double> scaling = getSpectralScalingParameters(0.01);
  const double c = scaling.first;
  const double w = scaling.second;

  Eigen::VectorXd avgMoments =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numMoments));

  std::mt19937_64 rng(42);
  std::uniform_int_distribution<int> coin(0, 1);

  const double scale = 1.0 / std::sqrt(static_cast<double>(n));

  std::vector<Eigen::VectorXd> randomVectors(numRandomVectors,
                                             Eigen::VectorXd(n));
  for (size_t r = 0; r < numRandomVectors; ++r) {
    for (Eigen::Index i = 0; i < n; ++i) {
      randomVectors[r][i] = coin(rng) ? scale : -scale;
    }
  }

#ifdef OPENMP_FOUND
#pragma omp parallel for reduction(+ : avgMoments) schedule(dynamic)
#endif
  for (size_t r = 0; r < numRandomVectors; ++r) {
    const Eigen::VectorXd& v = randomVectors[r];
    Eigen::VectorXd u0 = v;
    avgMoments[0] += v.dot(u0);

    if (numMoments > 1) {
      Eigen::VectorXd y = this->assembledConnectivityMatrix_ * u0;
      Eigen::VectorXd u1 = (y - c * u0) / w;
      avgMoments[1] += v.dot(u1);

      Eigen::VectorXd uPrev = u0;
      Eigen::VectorXd uCurr = u1;
      Eigen::VectorXd uNext(n);
      const Eigen::Index numMomentsIdx = static_cast<Eigen::Index>(numMoments);
      for (Eigen::Index m = 1; m < numMomentsIdx - 1; ++m) {
        y.noalias() = this->assembledConnectivityMatrix_ * uCurr;
        uNext = (2.0 / w) * (y - c * uCurr) - uPrev;
        avgMoments[m + 1] += v.dot(uNext);
        uPrev.swap(uCurr);
        uCurr.swap(uNext);
      }
    }
  }

  avgMoments /= static_cast<double>(numRandomVectors);
  avgMoments.array() *= computeJacksonKernel(numMoments).array();
  return avgMoments;
}

Eigen::VectorXd
NormalModeAnalyzer::chebyshevMoments(const KPMConfig& config) const
{
  INVALIDARG_EXP_IFN(config.numMoments > 0,
                     "numMoments must be greater than 0");
  INVALIDARG_EXP_IFN(config.numRandomVectors > 0,
                     "numRandomVectors must be greater than 0");
  INVALIDARG_EXP_IFN(config.batchSize > 0, "batchSize must be greater than 0");
  INVALIDARG_EXP_IFN(config.backend != ComputeBackend::CUDA,
                     "ComputeBackend::CUDA is not supported for KPM; use "
                     "AUTO, CPU, or METAL.");

  if (this->assembledConnectivityMatrix_.rows() == 0) {
    return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(config.numMoments));
  }

  const bool wantMetal = (config.backend == ComputeBackend::METAL ||
                          config.backend == ComputeBackend::AUTO);
  if (wantMetal) {
#ifdef USE_METAL
    if (isMetalAvailable()) {
      const std::pair<double, double> scaling =
        getSpectralScalingParameters(0.01);
      return pylimer_tools::calc::computeChebyshevMomentsKPMMetal(
        this->assembledConnectivityMatrix_,
        config.numMoments,
        config.numRandomVectors,
        scaling.first,
        scaling.second,
        config.batchSize);
    }
#endif
  }
  return computeChebyshevMomentsCPU(config.numMoments, config.numRandomVectors);
}

Eigen::ArrayXd
NormalModeAnalyzer::densityOfStates(
  const Eigen::Ref<const Eigen::ArrayXd>& lambdaGrid,
  const KPMConfig& config) const
{
  const Eigen::Index n = this->assembledConnectivityMatrix_.rows();
  if (n == 0) {
    return Eigen::ArrayXd::Zero(lambdaGrid.size());
  }

  const std::pair<double, double> scaling = getSpectralScalingParameters(0.01);
  const double c = scaling.first;
  const double w = scaling.second;
  const Eigen::VectorXd mu = chebyshevMoments(config);

  Eigen::ArrayXd dos = Eigen::ArrayXd::Zero(lambdaGrid.size());
  for (Eigen::Index i = 0; i < lambdaGrid.size(); ++i) {
    const double lam = lambdaGrid[i];
    const double x = (lam - c) / w;
    if (x <= -1.0 || x >= 1.0) {
      continue;
    }

    double tPrev = 1.0;
    double tCurr = x;
    double sum = mu[0];
    if (config.numMoments > 1) {
      sum += 2.0 * mu[1] * tCurr;
      const Eigen::Index numMomentsIdx =
        static_cast<Eigen::Index>(config.numMoments);
      for (Eigen::Index m = 1; m < numMomentsIdx - 1; ++m) {
        const double tNext = 2.0 * x * tCurr - tPrev;
        sum += 2.0 * mu[m + 1] * tNext;
        tPrev = tCurr;
        tCurr = tNext;
      }
    }
    const double factor = 1.0 / (w * M_PI * std::sqrt(1.0 - x * x));
    dos[i] = static_cast<double>(n) * factor * std::max(0.0, sum);
  }
  return dos;
}

Eigen::SparseMatrix<double>
NormalModeAnalyzer::extractComponentSubmatrix(
  const std::vector<size_t>& compNodes,
  std::vector<Eigen::Index>& localIdxBuffer,
  std::vector<size_t>* outLocalDeg,
  size_t* outNumInternalEdges,
  bool* outAllUnitWeights) const
{
  const size_t M = compNodes.size();
  for (size_t i = 0; i < M; ++i) {
    localIdxBuffer[compNodes[i]] = static_cast<Eigen::Index>(i);
  }

  if (outLocalDeg != nullptr) {
    outLocalDeg->assign(M, 0);
  }
  if (outNumInternalEdges != nullptr) {
    *outNumInternalEdges = 0;
  }
  if (outAllUnitWeights != nullptr) {
    *outAllUnitWeights = true;
  }

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(M * 4);
  for (const size_t u : compNodes) {
    const Eigen::Index lu = localIdxBuffer[u];
    for (Eigen::SparseMatrix<double>::InnerIterator it(
           this->assembledConnectivityMatrix_, static_cast<Eigen::Index>(u));
         it;
         ++it) {
      const size_t v = static_cast<size_t>(it.row());
      const Eigen::Index lv = localIdxBuffer[v];
      triplets.emplace_back(
        static_cast<int>(lu), static_cast<int>(lv), it.value());
      if (v > u) {
        if (outNumInternalEdges != nullptr) {
          (*outNumInternalEdges)++;
        }
        if (outLocalDeg != nullptr) {
          (*outLocalDeg)[static_cast<size_t>(lu)]++;
          (*outLocalDeg)[static_cast<size_t>(lv)]++;
        }
        if (outAllUnitWeights != nullptr &&
            std::abs(it.value() - (-1.0)) > 1e-9) {
          *outAllUnitWeights = false;
        }
      }
    }
  }

  Eigen::SparseMatrix<double> subMat(static_cast<Eigen::Index>(M),
                                     static_cast<Eigen::Index>(M));
  subMat.setFromTriplets(triplets.begin(), triplets.end());
  return subMat;
}

void
NormalModeAnalyzer::computeConnectedComponents(
  size_t& numComponents,
  std::vector<size_t>& labels,
  std::vector<size_t>& counts) const
{
  std::vector<std::vector<size_t>> components;
  pylimer_tools::utils::computeConnectedComponents(
    this->assembledConnectivityMatrix_, components, &labels);
  numComponents = components.size();
  counts.resize(numComponents);
  for (size_t c = 0; c < numComponents; ++c) {
    counts[c] = components[c].size();
  }
}

void
NormalModeAnalyzer::runDeflatedSLQ(const size_t numSteps,
                                   const size_t numRandomVectors,
                                   const size_t numComponents,
                                   const std::vector<size_t>& labels,
                                   const std::vector<size_t>& counts,
                                   const Eigen::MatrixXd* deflationBasis,
                                   std::vector<double>& lambdaNodes,
                                   std::vector<double>& deltaN) const
{
  const Eigen::Index n = this->assembledConnectivityMatrix_.rows();
  if (n == 0 || numSteps == 0 || numRandomVectors == 0) {
    return;
  }

  std::mt19937_64 rng(42);
  std::uniform_int_distribution<int> coin(0, 1);

  std::vector<Eigen::VectorXd> randomVectors(numRandomVectors,
                                             Eigen::VectorXd(n));
  for (size_t r = 0; r < numRandomVectors; ++r) {
    Eigen::VectorXd& v = randomVectors[r];
    for (Eigen::Index i = 0; i < n; ++i) {
      v[i] = coin(rng) ? 1.0 : -1.0;
    }
    std::vector<double> compSums(numComponents, 0.0);
    for (Eigen::Index i = 0; i < n; ++i) {
      compSums[labels[static_cast<size_t>(i)]] += v[i];
    }
    for (Eigen::Index i = 0; i < n; ++i) {
      const size_t c = labels[static_cast<size_t>(i)];
      v[i] -= compSums[c] / static_cast<double>(counts[c]);
    }
    if (deflationBasis != nullptr && deflationBasis->cols() > 0) {
      v.noalias() -= (*deflationBasis) * (deflationBasis->transpose() * v);
    }
  }

  const size_t maxSteps = std::min(numSteps, static_cast<size_t>(n));
  std::vector<std::vector<double>> perVecLambdas(numRandomVectors);
  std::vector<std::vector<double>> perVecWeights(numRandomVectors);

#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic)
#endif
  for (size_t r = 0; r < numRandomVectors; ++r) {
    const Eigen::VectorXd& v = randomVectors[r];
    const double normV = v.norm();
    if (normV < 1e-13) {
      continue;
    }

    Eigen::VectorXd q = v / normV;
    Eigen::VectorXd qPrev = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd w(n);
    std::vector<double> alpha;
    std::vector<double> beta;
    alpha.reserve(maxSteps);
    beta.reserve(maxSteps > 0 ? maxSteps - 1 : 0);

    for (size_t j = 0; j < maxSteps; ++j) {
      w.noalias() = this->assembledConnectivityMatrix_ * q;
      if (j > 0) {
        w -= beta.back() * qPrev;
      }
      const double a = q.dot(w);
      alpha.push_back(a);
      w -= a * q;
      if (j + 1 < maxSteps) {
        const double b = w.norm();
        if (b < 1e-13) {
          break;
        }
        beta.push_back(b);
        qPrev.swap(q);
        q = w / b;
      }
    }

    const Eigen::Index mActual = static_cast<Eigen::Index>(alpha.size());
    if (mActual == 0) {
      continue;
    }
    const double normV2 = normV * normV;
    const double invR = 1.0 / static_cast<double>(numRandomVectors);

    if (mActual == 1) {
      if (alpha[0] > 1e-11) {
        perVecLambdas[r].push_back(alpha[0]);
        perVecWeights[r].push_back(normV2 * invR);
      }
      continue;
    }

    const Eigen::VectorXd diag =
      Eigen::Map<const Eigen::VectorXd>(alpha.data(), mActual);
    const Eigen::VectorXd subdiag =
      Eigen::Map<const Eigen::VectorXd>(beta.data(), mActual - 1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> tridiagSolver;
    tridiagSolver.computeFromTridiagonal(
      diag, subdiag, Eigen::ComputeEigenvectors);
    if (tridiagSolver.info() != Eigen::Success) {
      continue;
    }

    const Eigen::VectorXd& evals = tridiagSolver.eigenvalues();
    const Eigen::MatrixXd& evecs = tridiagSolver.eigenvectors();
    perVecLambdas[r].reserve(static_cast<size_t>(mActual));
    perVecWeights[r].reserve(static_cast<size_t>(mActual));
    for (Eigen::Index j = 0; j < mActual; ++j) {
      const double lam = evals[j];
      if (lam > 1e-11) {
        const double q0j = evecs(0, j);
        perVecLambdas[r].push_back(lam);
        perVecWeights[r].push_back(normV2 * q0j * q0j * invR);
      }
    }
  }

  for (size_t r = 0; r < numRandomVectors; ++r) {
    lambdaNodes.insert(
      lambdaNodes.end(), perVecLambdas[r].begin(), perVecLambdas[r].end());
    deltaN.insert(
      deltaN.end(), perVecWeights[r].begin(), perVecWeights[r].end());
  }
}

void
NormalModeAnalyzer::computeEffectiveSpectrumHybrid(
  const size_t numSteps,
  const size_t numRandomVectors,
  const size_t numExactLowModes,
  std::vector<double>& lambdaNodes,
  std::vector<double>& deltaN) const
{
  const Eigen::Index n = this->assembledConnectivityMatrix_.rows();
  if (n == 0) {
    return;
  }

  size_t numComponents = 0;
  std::vector<size_t> labels;
  std::vector<size_t> counts;
  computeConnectedComponents(numComponents, labels, counts);

  std::vector<std::vector<size_t>> componentNodes(numComponents);
  for (size_t c = 0; c < numComponents; ++c) {
    componentNodes[c].reserve(counts[c]);
  }
  for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
    componentNodes[labels[i]].push_back(i);
  }

  std::vector<Eigen::VectorXd> deflationCols;
  std::vector<Eigen::Index> localIdxBuffer(static_cast<size_t>(n), -1);

  for (size_t c = 0; c < numComponents; ++c) {
    const std::vector<size_t>& nodes = componentNodes[c];
    const size_t compSize = nodes.size();
    if (compSize <= 1) {
      continue;
    }

    const Eigen::SparseMatrix<double> subMat =
      extractComponentSubmatrix(nodes, localIdxBuffer);

    auto appendDeflatedModes = [&](const Eigen::VectorXd& evals,
                                   const Eigen::MatrixXd& evecs) {
      for (Eigen::Index j = 0; j < evals.size(); ++j) {
        if (evals[j] > 1e-9) {
          lambdaNodes.push_back(evals[j]);
          deltaN.push_back(1.0);
          Eigen::VectorXd fullVec = Eigen::VectorXd::Zero(n);
          for (size_t i = 0; i < compSize; ++i) {
            fullVec[static_cast<Eigen::Index>(nodes[i])] =
              evecs(static_cast<Eigen::Index>(i), j);
          }
          deflationCols.push_back(std::move(fullVec));
        }
      }
    };

    if (compSize <= 64) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
        Eigen::MatrixXd(subMat), Eigen::ComputeEigenvectors);
      if (es.info() == Eigen::Success) {
        appendDeflatedModes(es.eigenvalues(), es.eigenvectors());
      }
    } else if (numExactLowModes > 0 && compSize >= 4) {
      const size_t kReq = std::min(numExactLowModes + 1, compSize - 2);
      const size_t ncv = std::min(2 * kReq + 1, compSize);
      try {
        Spectra::SparseSymShiftSolve<double> op(subMat);
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(
          op,
          static_cast<Eigen::Index>(kReq),
          static_cast<Eigen::Index>(ncv),
          -1e-5);
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn);
        if (eigs.info() == Spectra::CompInfo::Successful) {
          appendDeflatedModes(eigs.eigenvalues(), eigs.eigenvectors());
        }
      } catch (...) {
        // Fall back to SLQ on this component if shift-invert fails
      }
    }
  }

  Eigen::MatrixXd deflationBasis;
  if (!deflationCols.empty()) {
    deflationBasis.resize(n, static_cast<Eigen::Index>(deflationCols.size()));
    for (size_t j = 0; j < deflationCols.size(); ++j) {
      deflationBasis.col(static_cast<Eigen::Index>(j)) = deflationCols[j];
    }
  }

  runDeflatedSLQ(numSteps,
                 numRandomVectors,
                 numComponents,
                 labels,
                 counts,
                 deflationCols.empty() ? nullptr : &deflationBasis,
                 lambdaNodes,
                 deltaN);
}

NormalModeSpectrum
NormalModeAnalyzer::solveApprox() const
{
  const size_t n = this->getSize();
  if (n <= 30000) {
    return this->solveApprox(HybridConfig{});
  }
  return this->solveApprox(SLQConfig{});
}

NormalModeSpectrum
NormalModeAnalyzer::solveApprox(const SLQConfig& config) const
{
  INVALIDARG_EXP_IFN(config.numSteps > 0, "numSteps must be greater than 0");
  INVALIDARG_EXP_IFN(config.numRandomVectors > 0,
                     "numRandomVectors must be greater than 0");
  INVALIDARG_EXP_IFN(config.backend == ComputeBackend::AUTO ||
                       config.backend == ComputeBackend::CPU,
                     "Only ComputeBackend::AUTO and ComputeBackend::CPU are "
                     "supported for SLQ.");

  size_t numComponents = 0;
  std::vector<size_t> labels;
  std::vector<size_t> counts;
  computeConnectedComponents(numComponents, labels, counts);

  std::vector<double> lambdaNodes;
  std::vector<double> deltaN;
  runDeflatedSLQ(config.numSteps,
                 config.numRandomVectors,
                 numComponents,
                 labels,
                 counts,
                 nullptr,
                 lambdaNodes,
                 deltaN);

  Eigen::VectorXd evals = Eigen::Map<const Eigen::VectorXd>(
    lambdaNodes.data(), static_cast<Eigen::Index>(lambdaNodes.size()));
  Eigen::VectorXd weights = Eigen::Map<const Eigen::VectorXd>(
    deltaN.data(), static_cast<Eigen::Index>(deltaN.size()));
  return NormalModeSpectrum(
    std::move(evals), std::move(weights), std::nullopt, false);
}

NormalModeSpectrum
NormalModeAnalyzer::solveApprox(const HybridConfig& config) const
{
  INVALIDARG_EXP_IFN(config.numSteps > 0, "numSteps must be greater than 0");
  INVALIDARG_EXP_IFN(config.numRandomVectors > 0,
                     "numRandomVectors must be greater than 0");
  INVALIDARG_EXP_IFN(config.backend == ComputeBackend::AUTO ||
                       config.backend == ComputeBackend::CPU,
                     "Only ComputeBackend::AUTO and ComputeBackend::CPU are "
                     "supported for Hybrid.");

  std::vector<double> lambdaNodes;
  std::vector<double> deltaN;
  computeEffectiveSpectrumHybrid(config.numSteps,
                                 config.numRandomVectors,
                                 config.numExactLowModes,
                                 lambdaNodes,
                                 deltaN);

  Eigen::VectorXd evals = Eigen::Map<const Eigen::VectorXd>(
    lambdaNodes.data(), static_cast<Eigen::Index>(lambdaNodes.size()));
  Eigen::VectorXd weights = Eigen::Map<const Eigen::VectorXd>(
    deltaN.data(), static_cast<Eigen::Index>(deltaN.size()));
  return NormalModeSpectrum(
    std::move(evals), std::move(weights), std::nullopt, false);
}

NormalModeSpectrum
NormalModeAnalyzer::solveApprox(const KPMConfig& config) const
{
  INVALIDARG_EXP_IFN(config.numMoments > 0,
                     "numMoments must be greater than 0");
  INVALIDARG_EXP_IFN(config.numRandomVectors > 0,
                     "numRandomVectors must be greater than 0");
  INVALIDARG_EXP_IFN(config.numQuadraturePoints > 0,
                     "numQuadraturePoints must be greater than 0");

  const Eigen::Index n = this->assembledConnectivityMatrix_.rows();
  if (n == 0) {
    return NormalModeSpectrum(
      Eigen::VectorXd::Zero(0), Eigen::VectorXd::Zero(0), std::nullopt, false);
  }

  const std::pair<double, double> scaling = getSpectralScalingParameters(0.01);
  const double c = scaling.first;
  const double w = scaling.second;
  const Eigen::VectorXd mu = chebyshevMoments(config);

  const size_t qCount = config.numQuadraturePoints;
  std::vector<double> lambdaNodes;
  std::vector<double> deltaN;
  lambdaNodes.reserve(qCount);
  deltaN.reserve(qCount);

  const Eigen::Index numMomentsIdx =
    static_cast<Eigen::Index>(config.numMoments);
  for (size_t k = 1; k <= qCount; ++k) {
    const double thetaK =
      (static_cast<double>(k) - 0.5) * M_PI / static_cast<double>(qCount);
    const double xK = std::cos(thetaK);
    const double lamK = c + w * xK;
    if (lamK <= 1e-4) {
      continue;
    }

    double sum = mu[0];
    for (Eigen::Index m = 1; m < numMomentsIdx; ++m) {
      sum += 2.0 * mu[m] * std::cos(static_cast<double>(m) * thetaK);
    }
    const double dn = (static_cast<double>(n) / static_cast<double>(qCount)) *
                      std::max(0.0, sum);
    lambdaNodes.push_back(lamK);
    deltaN.push_back(dn);
  }

  Eigen::VectorXd evals = Eigen::Map<const Eigen::VectorXd>(
    lambdaNodes.data(), static_cast<Eigen::Index>(lambdaNodes.size()));
  Eigen::VectorXd weights = Eigen::Map<const Eigen::VectorXd>(
    deltaN.data(), static_cast<Eigen::Index>(deltaN.size()));
  return NormalModeSpectrum(
    std::move(evals), std::move(weights), std::nullopt, false);
}

} // namespace pylimer_tools::calc
