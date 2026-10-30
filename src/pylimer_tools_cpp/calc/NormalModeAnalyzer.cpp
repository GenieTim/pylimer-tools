#include "./NormalModeAnalyzer.h"

#include "../utils/utilityMacros.h"
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
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>
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
NormalModeAnalyzer::NormalModeAnalyzer(const std::vector<size_t> springFrom,
                                       const std::vector<size_t> springTo)
{
  INVALIDARG_EXP_IFN(springFrom.size() == springTo.size(),
                     "springFrom and springTo must have the same size");
  INVALIDARG_EXP_IFN(springFrom.size() > 0,
                     "springFrom and springTo cannot be empty");
  // find maximum index = nr. of cols/rows in the connectivity matrix
  size_t maxIdx =
    std::max(*std::max_element(springFrom.begin(), springFrom.end()),
             *std::max_element(springTo.begin(), springTo.end()));
  size_t nRows = maxIdx + 1;
  // assemble connectivity matrix
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(springFrom.size() * 2);

  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(nRows);
  for (size_t i = 0; i < springFrom.size(); i++) {
    // some sanity checks
    assert(springFrom[i] < nRows);
    assert(springTo[i] < nRows);
    // triplets will be summed up -> we can use the same indices multiple
    // times
    triplets.push_back(Eigen::Triplet<double>(
      static_cast<int>(springFrom[i]), static_cast<int>(springTo[i]), -1.0));
    triplets.push_back(Eigen::Triplet<double>(
      static_cast<int>(springTo[i]), static_cast<int>(springFrom[i]), -1.0));
    // it is a bit more efficient to sum up the diagonal elements manually
    diagonal[static_cast<Eigen::Index>(springFrom[i])] += 1.0;
    diagonal[static_cast<Eigen::Index>(springTo[i])] += 1.0;
  }
  // the summed up diagonal can be translated to triplets as well
  for (Eigen::Index i = 0; i < diagonal.size(); i++) {
    if (diagonal[i] != 0.) {
      triplets.push_back(Eigen::Triplet<double>(i, i, diagonal[i]));
    }
  }

  // finally, translate into the sparse matrix format
  this->assembledConnectivityMatrix = Eigen::SparseMatrix<double>(nRows, nRows);
  this->assembledConnectivityMatrix.setFromTriplets(triplets.begin(),
                                                    triplets.end());
};

void
NormalModeAnalyzer::findSparseEigenvalues(const size_t nrOfEigenvalues,
                                          const bool includeEigenvectors)
{
  INVALIDARG_EXP_IFN(nrOfEigenvalues > 0,
                     "nrOfEigenvalues must be greater than 0");
  size_t n = this->assembledConnectivityMatrix.rows();
  INVALIDARG_EXP_IFN(nrOfEigenvalues < n,
                     "nrOfEigenvalues must be less than matrix size");

  try {
    Spectra::SparseSymShiftSolve<double> op(this->assembledConnectivityMatrix);
    size_t ncv = std::min(2 * nrOfEigenvalues + 1, n);
    if (ncv <= nrOfEigenvalues) {
      ncv = std::min(nrOfEigenvalues + 2, n);
    }
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(
      op, nrOfEigenvalues, ncv, -1e-5);
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestMagn);
    if (eigs.info() == Spectra::CompInfo::Successful) {
      this->setEigenvalues(eigs.eigenvalues());
      if (includeEigenvectors) {
        this->setEigenvectors(eigs.eigenvectors());
      }
      return;
    }
  } catch (...) {
    // Fallback to standard solver below
  }

  // Fallback to standard solver
  Spectra::SparseGenMatProd<double> op =
    Spectra::SparseGenMatProd<double>(this->assembledConnectivityMatrix);
  size_t ncv = std::min(2 * nrOfEigenvalues, n);
  if (ncv <= nrOfEigenvalues) {
    ncv = std::min(nrOfEigenvalues + 2, n);
  }
  Spectra::SymEigsSolver<Spectra::SparseGenMatProd<double>> eigs(
    op, nrOfEigenvalues, ncv);
  eigs.init();
  eigs.compute(Spectra::SortRule::SmallestAlge);

  RUNTIME_EXP_IFN(eigs.info() == Spectra::CompInfo::Successful,
                  "Spectra failed to compute sparse eigenvalues");
  this->setEigenvalues(eigs.eigenvalues());
  if (includeEigenvectors) {
    this->setEigenvectors(eigs.eigenvectors());
  }
};

void
NormalModeAnalyzer::computeAllEigenvaluesDense(
  const Eigen::SparseMatrix<double>& mat,
  const bool includeEigenvectors,
  const bool forceEigenSolver,
  Eigen::VectorXd& outEvals,
  Eigen::MatrixXd& outEvecs)
{
  Eigen::Index nRows = mat.rows();
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
  lapack_int size = static_cast<lapack_int>(nRows);

  if (includeEigenvectors) {
    int info = LAPACKE_dsyevd(Eigen::MatrixXd::IsRowMajor ? LAPACK_ROW_MAJOR
                                                          : LAPACK_COL_MAJOR,
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
    lapack_int ldz = 1;
    lapack_int il = 0;
    lapack_int iu = 0;
    double abstol = 1e-10;
    lapack_int M = 0;
    Eigen::Vector<lapack_int, Eigen::Dynamic> support =
      Eigen::Vector<lapack_int, Eigen::Dynamic>::Zero(nRows * 2);

    int info = LAPACKE_dsyevr(Eigen::MatrixXd::IsRowMajor ? LAPACK_ROW_MAJOR
                                                          : LAPACK_COL_MAJOR,
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
                    "Error in LAPACK_dsyevr. Exit code " +
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
  Eigen::Index n = mat.rows();
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

  // If system is small (<= 30), dense solver is faster and lightweight
  if (n <= 30) {
    computeAllEigenvaluesDense(
      mat, includeEigenvectors, false, outEvals, outEvecs);
    return;
  }

  double maxDiag = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    maxDiag = std::max(maxDiag, mat.coeff(i, i));
  }
  double b_upper = 2.0 * maxDiag + 0.1;
  double a_lower = -1e-7;

  // Root counting via Sylvester's Law of Inertia
  auto countLessThan = [&](double sigma) -> size_t {
    const double perturbations[] = {
      0.0, 1.234567e-8, -2.345678e-8, 5.678901e-7
    };
    for (double pert : perturbations) {
      double current_sigma = sigma + pert;
      Eigen::SparseMatrix<double> shifted = mat;
      for (Eigen::Index i = 0; i < shifted.rows(); ++i) {
        shifted.coeffRef(i, i) -= current_sigma;
      }
      shifted.makeCompressed();

      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
      ldlt.compute(shifted);
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

  size_t totalCount = static_cast<size_t>(n);
  std::vector<Slice> stack;
  stack.push_back({ a_lower, b_upper, 0, totalCount });

  std::vector<double> eigenvaluesList;
  eigenvaluesList.reserve(totalCount);

  size_t targetSliceSize =
    std::min(static_cast<size_t>(50), static_cast<size_t>(n));

  while (!stack.empty()) {
    Slice slice = stack.back();
    stack.pop_back();

    size_t k = slice.countB - slice.countA;
    if (k == 0) {
      continue;
    }

    if ((slice.b - slice.a) < 1e-6) {
      double mid = (slice.a + slice.b) / 2.0;
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
          double mid = (low + high) / 2.0;
          size_t c = countLessThan(mid);
          if (c <= slice.countA) {
            low = mid;
          } else {
            high = mid;
          }
        }
        eigenvaluesList.push_back((low + high) / 2.0);
      } else {
        double mid = (slice.a + slice.b) / 2.0;
        size_t countMid = countLessThan(mid);
        stack.push_back({ mid, slice.b, countMid, slice.countB });
        stack.push_back({ slice.a, mid, slice.countA, countMid });
      }
      continue;
    }

    if (k <= targetSliceSize) {
      bool solvedBySpectra = false;
      try {
        double sigma = (slice.a + slice.b) / 2.0;
        Spectra::SparseSymShiftSolve<double> op(mat);
        size_t nev = k;
        size_t ncv = std::min(2 * nev + 2, static_cast<size_t>(n));
        if (ncv <= nev) {
          ncv = std::min(nev + 2, static_cast<size_t>(n));
        }
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(
          op, nev, ncv, sigma);
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn);
        if (eigs.info() == Spectra::CompInfo::Successful) {
          Eigen::VectorXd evals = eigs.eigenvalues();
          std::vector<double> in_slice;
          for (Eigen::Index i = 0; i < evals.size(); ++i) {
            if (evals[i] >= slice.a - 1e-7 && evals[i] < slice.b + 1e-7) {
              in_slice.push_back(evals[i]);
            }
          }
          if (in_slice.size() == k) {
            for (double val : in_slice) {
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

    // Split slice in half
    double mid = (slice.a + slice.b) / 2.0;
    size_t countMid = countLessThan(mid);
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
NormalModeAnalyzer::computeAllEigenvalues(const NormalModeSolver solver,
                                          const bool includeEigenvectors,
                                          const bool decomposeComponents)
{
  RUNTIME_EXP_IFN(this->assembledConnectivityMatrix.rows() ==
                    this->assembledConnectivityMatrix.cols(),
                  "Expected square matrix");

  Eigen::Index nRows = this->assembledConnectivityMatrix.rows();
  if (nRows == 0) {
    this->setEigenvalues(Eigen::VectorXd::Zero(0));
    if (includeEigenvectors) {
      this->setEigenvectors(Eigen::MatrixXd::Zero(0, 0));
    }
    return;
  }

  if (!decomposeComponents) {
    Eigen::VectorXd evals;
    Eigen::MatrixXd evecs;
    if (solver == NormalModeSolver::METAL) {
#ifdef USE_METAL
      computeAllEigenvaluesMetal(
        this->assembledConnectivityMatrix, includeEigenvectors, evals, evecs);
#else
      computeAllEigenvaluesDense(this->assembledConnectivityMatrix,
                                 includeEigenvectors,
                                 false,
                                 evals,
                                 evecs);
#endif
    } else if (solver == NormalModeSolver::SPARSE_SLICING ||
               (solver == NormalModeSolver::AUTO && nRows > 2000)) {
      computeAllEigenvaluesSparseSlicing(
        this->assembledConnectivityMatrix, includeEigenvectors, evals, evecs);
    } else {
      computeAllEigenvaluesDense(this->assembledConnectivityMatrix,
                                 includeEigenvectors,
                                 (solver == NormalModeSolver::DENSE_EIGEN),
                                 evals,
                                 evecs);
    }
    this->setEigenvalues(evals);
    if (includeEigenvectors) {
      this->setEigenvectors(evecs);
    }
    return;
  }

  // --- Component-by-Component Decomposition ---
  std::vector<size_t> parent(static_cast<size_t>(nRows));
  std::iota(parent.begin(), parent.end(), 0);
  auto find_root = [&](size_t i) {
    size_t root = i;
    while (root != parent[root])
      root = parent[root];
    size_t curr = i;
    while (curr != root) {
      size_t nxt = parent[curr];
      parent[curr] = root;
      curr = nxt;
    }
    return root;
  };
  auto unite = [&](size_t i, size_t j) {
    size_t ri = find_root(i);
    size_t rj = find_root(j);
    if (ri != rj)
      parent[ri] = rj;
  };

  for (Eigen::Index u = 0; u < nRows; ++u) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(
           this->assembledConnectivityMatrix, u);
         it;
         ++it) {
      Eigen::Index v = it.row();
      if (v > u) {
        unite(static_cast<size_t>(u), static_cast<size_t>(v));
      }
    }
  }

  std::unordered_map<size_t, std::vector<size_t>> componentsMap;
  for (size_t i = 0; i < static_cast<size_t>(nRows); ++i) {
    componentsMap[find_root(i)].push_back(i);
  }

  std::vector<double> allEigenvalues;
  allEigenvalues.reserve(static_cast<size_t>(nRows));

  for (auto& [root, comp] : componentsMap) {
    size_t M = comp.size();
    if (M == 1) {
      allEigenvalues.push_back(0.0);
      continue;
    }

    std::unordered_map<size_t, size_t> localIdx;
    localIdx.reserve(M);
    for (size_t i = 0; i < M; ++i) {
      localIdx[comp[i]] = i;
    }

    std::vector<size_t> localDeg(M, 0);
    size_t numInternalEdges = 0;
    std::vector<Eigen::Triplet<double>> subTriplets;
    subTriplets.reserve(M * 4);
    bool allUnitWeights = true;

    for (size_t u : comp) {
      size_t lu = localIdx[u];
      for (Eigen::SparseMatrix<double>::InnerIterator it(
             this->assembledConnectivityMatrix, static_cast<Eigen::Index>(u));
           it;
           ++it) {
        size_t v = static_cast<size_t>(it.row());
        size_t lv = localIdx[v];
        subTriplets.push_back(Eigen::Triplet<double>(
          static_cast<int>(lu), static_cast<int>(lv), it.value()));
        if (v > u) {
          numInternalEdges++;
          localDeg[lu]++;
          localDeg[lv]++;
          if (std::abs(it.value() - (-1.0)) > 1e-9) {
            allUnitWeights = false;
          }
        }
      }
    }

    if (M == 2 && allUnitWeights && numInternalEdges == 1) {
      allEigenvalues.push_back(0.0);
      allEigenvalues.push_back(2.0);
      continue;
    }

    // Check if linear path (Rouse chain)
    size_t deg1Count = 0;
    size_t deg2Count = 0;
    for (size_t d : localDeg) {
      if (d == 1)
        deg1Count++;
      else if (d == 2)
        deg2Count++;
    }

    if (allUnitWeights && numInternalEdges == M - 1 && deg1Count == 2 &&
        deg2Count == M - 2) {
      for (size_t p = 0; p < M; ++p) {
        double val = 4.0 * std::pow(std::sin((static_cast<double>(p) * M_PI) /
                                             (2.0 * static_cast<double>(M))),
                                    2);
        allEigenvalues.push_back(val);
      }
      continue;
    }

    // Check if simple cycle (ring)
    if (allUnitWeights && numInternalEdges == M && deg2Count == M) {
      for (size_t p = 0; p < M; ++p) {
        double val = 4.0 * std::pow(std::sin((static_cast<double>(p) * M_PI) /
                                             static_cast<double>(M)),
                                    2);
        allEigenvalues.push_back(val);
      }
      continue;
    }

    // General subcomponent
    Eigen::SparseMatrix<double> subMat(static_cast<Eigen::Index>(M),
                                       static_cast<Eigen::Index>(M));
    subMat.setFromTriplets(subTriplets.begin(), subTriplets.end());

    Eigen::VectorXd subEvals;
    Eigen::MatrixXd subEvecs;
    if (solver == NormalModeSolver::METAL) {
#ifdef USE_METAL
      computeAllEigenvaluesMetal(
        subMat, includeEigenvectors, subEvals, subEvecs);
#else
      computeAllEigenvaluesDense(
        subMat, includeEigenvectors, false, subEvals, subEvecs);
#endif
    } else if (solver == NormalModeSolver::SPARSE_SLICING ||
               (solver == NormalModeSolver::AUTO && M > 2000)) {
      computeAllEigenvaluesSparseSlicing(
        subMat, includeEigenvectors, subEvals, subEvecs);
    } else {
      computeAllEigenvaluesDense(subMat,
                                 includeEigenvectors,
                                 (solver == NormalModeSolver::DENSE_EIGEN),
                                 subEvals,
                                 subEvecs);
    }

    for (Eigen::Index i = 0; i < subEvals.size(); ++i) {
      allEigenvalues.push_back(subEvals[i]);
    }
  }

  std::sort(allEigenvalues.begin(), allEigenvalues.end());
  Eigen::VectorXd finalEvals = Eigen::Map<Eigen::VectorXd>(
    allEigenvalues.data(), static_cast<Eigen::Index>(allEigenvalues.size()));
  this->setEigenvalues(finalEvals);

  if (includeEigenvectors) {
    Eigen::VectorXd dummyEvals;
    Eigen::MatrixXd fullEvecs;
#ifdef USE_METAL
    if (solver == NormalModeSolver::METAL) {
      computeAllEigenvaluesMetal(
        this->assembledConnectivityMatrix, true, dummyEvals, fullEvecs);
    } else {
      computeAllEigenvaluesDense(this->assembledConnectivityMatrix,
                                 true,
                                 (solver == NormalModeSolver::DENSE_EIGEN),
                                 dummyEvals,
                                 fullEvecs);
    }
#else
    computeAllEigenvaluesDense(this->assembledConnectivityMatrix,
                               true,
                               (solver == NormalModeSolver::DENSE_EIGEN),
                               dummyEvals,
                               fullEvecs);
#endif
    this->setEigenvectors(fullEvecs);
  }
}

void
NormalModeAnalyzer::computeAllEigenvaluesGPU(const bool includeEigenvectors)
{
  RUNTIME_EXP_IFN(this->assembledConnectivityMatrix.rows() ==
                    this->assembledConnectivityMatrix.cols(),
                  "Expected square matrix");

  Eigen::Index nRows = this->assembledConnectivityMatrix.rows();
  if (nRows == 0) {
    return;
  }

#ifdef USE_CUDA
  Eigen::MatrixXd assembledConnectivityMatrixDense =
    Eigen::MatrixXd(this->assembledConnectivityMatrix);

  cusolverDnHandle_t cusolverH = nullptr;
  cusolverStatus_t cusolver_status = cusolverDnCreate(&cusolverH);
  RUNTIME_EXP_IFN(cusolver_status == CUSOLVER_STATUS_SUCCESS,
                  "cuSOLVER initialization failed");

  double* d_A = nullptr;
  double* d_W = nullptr;
  int* d_info = nullptr;

  cudaError_t cuda_stat =
    cudaMalloc((void**)&d_A, sizeof(double) * nRows * nRows);
  RUNTIME_EXP_IFN(cuda_stat == cudaSuccess, "cudaMalloc failed for d_A");
  cuda_stat = cudaMalloc((void**)&d_W, sizeof(double) * nRows);
  RUNTIME_EXP_IFN(cuda_stat == cudaSuccess, "cudaMalloc failed for d_W");
  cuda_stat = cudaMalloc((void**)&d_info, sizeof(int));
  RUNTIME_EXP_IFN(cuda_stat == cudaSuccess, "cudaMalloc failed for d_info");

  cuda_stat = cudaMemcpy(d_A,
                         assembledConnectivityMatrixDense.data(),
                         sizeof(double) * nRows * nRows,
                         cudaMemcpyHostToDevice);
  RUNTIME_EXP_IFN(cuda_stat == cudaSuccess, "cudaMemcpy failed to device");

  int lwork = 0;
  cusolver_status = cusolverDnDsyevd_bufferSize(
    cusolverH,
    CUSOLVER_EIG_MODE_VECTOR, // ALWAYS ask for vector size to be safe
    CUBLAS_FILL_MODE_LOWER,
    nRows,
    d_A,
    nRows,
    d_W,
    &lwork);
  RUNTIME_EXP_IFN(cusolver_status == CUSOLVER_STATUS_SUCCESS,
                  "cuSOLVER buffer size calculation failed");

  double* d_work = nullptr;
  cuda_stat = cudaMalloc((void**)&d_work, sizeof(double) * lwork);
  RUNTIME_EXP_IFN(cuda_stat == cudaSuccess, "cudaMalloc failed for d_work");

  cusolver_status = cusolverDnDsyevd(
    cusolverH,
    includeEigenvectors ? CUSOLVER_EIG_MODE_VECTOR : CUSOLVER_EIG_MODE_NOVECTOR,
    CUBLAS_FILL_MODE_LOWER,
    nRows,
    d_A,
    nRows,
    d_W,
    d_work,
    lwork,
    d_info);
  RUNTIME_EXP_IFN(cusolver_status == CUSOLVER_STATUS_SUCCESS,
                  "cuSOLVER dsyevd failed");

  cuda_stat = cudaDeviceSynchronize();
  RUNTIME_EXP_IFN(cuda_stat == cudaSuccess, "cudaDeviceSynchronize failed");

  Eigen::VectorXd eigenvaluesMemory = Eigen::VectorXd::Zero(nRows);
  Eigen::MatrixXd eigenvectorsMemory = includeEigenvectors
                                         ? Eigen::MatrixXd::Zero(nRows, nRows)
                                         : Eigen::MatrixXd::Zero(nRows, 3);

  cuda_stat = cudaMemcpy(eigenvaluesMemory.data(),
                         d_W,
                         sizeof(double) * nRows,
                         cudaMemcpyDeviceToHost);
  RUNTIME_EXP_IFN(cuda_stat == cudaSuccess,
                  "cudaMemcpy failed to host (eigenvalues)");

  if (includeEigenvectors) {
    cuda_stat = cudaMemcpy(eigenvectorsMemory.data(),
                           d_A,
                           sizeof(double) * nRows * nRows,
                           cudaMemcpyDeviceToHost);
    RUNTIME_EXP_IFN(cuda_stat == cudaSuccess,
                    "cudaMemcpy failed to host (eigenvectors)");
  }

  int info = 0;
  cuda_stat = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
  RUNTIME_EXP_IFN(
    info == 0, "cuSOLVER dsyevd info return non-zero: " + std::to_string(info));

  cudaFree(d_work);
  cudaFree(d_info);
  cudaFree(d_W);
  cudaFree(d_A);
  cusolverDnDestroy(cusolverH);

  this->setEigenvalues(eigenvaluesMemory);
  if (includeEigenvectors) {
    this->setEigenvectors(eigenvectorsMemory);
  }
#elif defined(USE_METAL)
  Eigen::VectorXd outEvals;
  Eigen::MatrixXd outEvecs;
  computeAllEigenvaluesMetal(
    this->assembledConnectivityMatrix, includeEigenvectors, outEvals, outEvecs);
  this->setEigenvalues(outEvals);
  if (includeEigenvectors) {
    this->setEigenvectors(outEvecs);
  }
#else
  std::cerr << "GPU computation requested but neither USE_CUDA nor USE_METAL "
               "was defined. Falling back to CPU..."
            << std::endl;
  this->computeAllEigenvalues(includeEigenvectors);
#endif
}

void
NormalModeAnalyzer::requireEigenvaluesComputation() const
{
  RUNTIME_EXP_IFN(this->isEigenvaluesComputed,
                  "Eigenvalues have not been computed yet");
}

Eigen::VectorXd
NormalModeAnalyzer::getEigenvalues() const
{
  this->requireEigenvaluesComputation();
  return this->eigenvalues;
}

void
NormalModeAnalyzer::setEigenvalues(const Eigen::VectorXd newEigenvalues)
{
  this->eigenvalues = newEigenvalues;
  for (Eigen::Index i = 0; i < this->eigenvalues.size(); ++i) {
    if (std::abs(this->eigenvalues[i]) < NORMAL_MODE_ZERO_TOLERANCE) {
      this->eigenvalues[i] = 0.0;
    }
  }
  this->isEigenvaluesComputed = true;
  this->clusterCount = this->countSolubleClusters();

  // TODO: possibly validate vector dimensions
}

void
NormalModeAnalyzer::requireEigenvectorsComputation() const
{
  RUNTIME_EXP_IFN(this->isEigenvectorsComputed,
                  "Eigenvectors have not been computed yet");
}

Eigen::MatrixXd
NormalModeAnalyzer::getEigenvectors() const
{
  this->requireEigenvectorsComputation();
  return this->eigenvectors;
}

void
NormalModeAnalyzer::setEigenvectors(const Eigen::MatrixXd newEigenvectors)
{
  this->eigenvectors = newEigenvectors;
  this->isEigenvectorsComputed = true;
  // TODO: possibly validate matrix dimensions
}

size_t
NormalModeAnalyzer::countSolubleClusters() const
{
  this->requireEigenvaluesComputation();
  return (this->eigenvalues.array() == 0.0).count();
}

size_t
NormalModeAnalyzer::getNrOfSolubleClusters() const
{
  return this->clusterCount;
}

Eigen::SparseMatrix<double>
NormalModeAnalyzer::getAssembledConnectivityMatrix() const
{
  return this->assembledConnectivityMatrix;
}

Eigen::ArrayXd
NormalModeAnalyzer::evaluateStressAutocorrelation(const Eigen::ArrayXd& t) const
{
  Eigen::ArrayXd result = Eigen::ArrayXd::Zero(t.size());
  for (Eigen::Index i = 0; i < this->eigenvalues.size(); ++i) {
    if (APPROX_EQUAL(this->eigenvalues[i], 0.0, NORMAL_MODE_ZERO_TOLERANCE)) {
      continue;
    }
    if (!std::isfinite(this->eigenvalues[i]) || this->eigenvalues[i] == 0.0) {
      continue; // skip non-finite eigenvalues
    }
    result += (-2. * this->eigenvalues[i] * t).exp();
    RUNTIME_EXP_IFN(result.allFinite(),
                    "Stress autocorrelation is not fully finite anymore "
                    "after adding Eigenvalue " +
                      std::to_string(i) + ": " +
                      std::to_string(this->eigenvalues[i]) + ".");
  }
  return result;
}

Eigen::ArrayXd
NormalModeAnalyzer::evaluateStorageModulus(const Eigen::ArrayXd& omega) const
{
  Eigen::ArrayXd result = Eigen::ArrayXd::Zero(omega.size());
  for (Eigen::Index i = 0; i < this->eigenvalues.size(); ++i) {
    if (APPROX_EQUAL(this->eigenvalues[i], 0.0, NORMAL_MODE_ZERO_TOLERANCE)) {
      continue;
    }
    if (!std::isfinite(this->eigenvalues[i]) || this->eigenvalues[i] == 0.0) {
      continue; // skip non-finite eigenvalues
    }
    result += (omega / (2. * this->eigenvalues[i])).square() /
              (1. + (omega / (2. * this->eigenvalues[i])).square());
    if (!result.allFinite()) {
      std::ostringstream oss;
      oss << "Storage modulus is not fully finite anymore after adding "
             "Eigenvalue "
          << i << ": " << this->eigenvalues[i] << ". Non-finite at indices: ";
      for (Eigen::Index idx = 0; idx < result.size(); ++idx) {
        if (!std::isfinite(result[idx])) {
          oss << idx << " ";
        }
      }
      RUNTIME_EXP_IFN(false, oss.str());
    }
  }
  return result;
}

Eigen::ArrayXd
NormalModeAnalyzer::evaluateLossModulus(const Eigen::ArrayXd& omega) const
{
  Eigen::ArrayXd result = Eigen::ArrayXd::Zero(omega.size());
  for (Eigen::Index i = 0; i < this->eigenvalues.size(); ++i) {
    if (APPROX_EQUAL(this->eigenvalues[i], 0.0, NORMAL_MODE_ZERO_TOLERANCE)) {
      continue;
    }
    if (!std::isfinite(this->eigenvalues[i]) || this->eigenvalues[i] == 0.0) {
      continue; // skip non-finite eigenvalues
    }
    result += (omega / (2. * this->eigenvalues[i])) /
              (1. + (omega / (2. * this->eigenvalues[i])).square());
    if (!result.allFinite()) {
      std::ostringstream oss;
      oss << "Loss modulus is not fully finite anymore after adding Eigenvalue "
          << i << ": " << this->eigenvalues[i] << ". Non-finite at indices: ";
      for (Eigen::Index idx = 0; idx < result.size(); ++idx) {
        if (!std::isfinite(result[idx])) {
          oss << idx << " ";
        }
      }
      RUNTIME_EXP_IFN(false, oss.str());
    }
  }
  return result;
}

std::pair<double, double>
NormalModeAnalyzer::getSpectralScalingParameters(double eps) const
{
  Eigen::Index n = this->assembledConnectivityMatrix.rows();
  if (n == 0) {
    return { 0.0, 1.0 };
  }
  double maxDiag = 0.0;
  for (Eigen::Index i = 0; i < n; ++i) {
    maxDiag = std::max(maxDiag, this->assembledConnectivityMatrix.coeff(i, i));
  }
  double lambdaMax = 2.0 * maxDiag + 0.1;
  double c = lambdaMax / 2.0;
  double w = lambdaMax / (2.0 * (1.0 - eps));
  return { c, w };
}

Eigen::VectorXd
NormalModeAnalyzer::computeChebyshevMomentsKPM(size_t numMoments,
                                               size_t numRandomVectors) const
{
  Eigen::Index n = this->assembledConnectivityMatrix.rows();
  if (n == 0 || numMoments == 0)
    return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numMoments));

  auto [c, w] = getSpectralScalingParameters(0.01);

  Eigen::VectorXd avgMoments =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numMoments));

  std::mt19937_64 rng(42);
  std::uniform_int_distribution<int> coin(0, 1);

  double scale = 1.0 / std::sqrt(static_cast<double>(n));

  // Preallocate work vectors to avoid heap allocation churn in the recurrence
  Eigen::VectorXd v(n);
  Eigen::VectorXd u0(n);
  Eigen::VectorXd u1(n);
  Eigen::VectorXd u_prev(n);
  Eigen::VectorXd u_curr(n);
  Eigen::VectorXd u_next(n);
  Eigen::VectorXd y(n);

  for (size_t r = 0; r < numRandomVectors; ++r) {
    for (Eigen::Index i = 0; i < n; ++i) {
      v[i] = coin(rng) ? scale : -scale;
    }

    u0 = v;
    avgMoments[0] += v.dot(u0);

    if (numMoments > 1) {
      y.noalias() = this->assembledConnectivityMatrix * u0;
      u1 = (y - c * u0) / w;
      avgMoments[1] += v.dot(u1);

      u_prev = u0;
      u_curr = u1;
      for (size_t m = 1; m < numMoments - 1; ++m) {
        y.noalias() = this->assembledConnectivityMatrix * u_curr;
        u_next = (2.0 / w) * (y - c * u_curr) - u_prev;
        avgMoments[static_cast<Eigen::Index>(m + 1)] += v.dot(u_next);
        u_prev = std::move(u_curr);
        u_curr = std::move(u_next);
      }
    }
  }

  avgMoments /= static_cast<double>(numRandomVectors);

  // Jackson damping kernel
  for (size_t m = 0; m < numMoments; ++m) {
    double m_d = static_cast<double>(m);
    double M_d = static_cast<double>(numMoments);
    double gm =
      ((M_d - m_d + 1.0) * std::cos(M_PI * m_d / (M_d + 1.0)) +
       std::sin(M_PI * m_d / (M_d + 1.0)) / std::tan(M_PI / (M_d + 1.0))) /
      (M_d + 1.0);
    avgMoments[static_cast<Eigen::Index>(m)] *= gm;
  }

  return avgMoments;
}

Eigen::VectorXd
NormalModeAnalyzer::computeChebyshevMomentsKPMMetal(
  size_t numMoments,
  size_t numRandomVectors) const
{
#ifdef USE_METAL
  auto [c, w] = getSpectralScalingParameters(0.01);
  return pylimer_tools::calc::computeChebyshevMomentsKPMMetal(
    this->assembledConnectivityMatrix, numMoments, numRandomVectors, c, w);
#else
  std::cerr << "Metal acceleration not enabled. Falling back to CPU KPM."
            << std::endl;
  return computeChebyshevMomentsKPM(numMoments, numRandomVectors);
#endif
}

void
NormalModeAnalyzer::computeEffectiveSpectrumKPM(
  size_t numMoments,
  size_t numRandomVectors,
  size_t numQuadraturePoints,
  std::vector<double>& lambdaNodes,
  std::vector<double>& deltaN) const
{
  Eigen::Index n = this->assembledConnectivityMatrix.rows();
  if (n == 0 || numMoments == 0 || numQuadraturePoints == 0)
    return;

  auto [c, w] = getSpectralScalingParameters(0.01);
  Eigen::VectorXd mu = computeChebyshevMomentsKPM(numMoments, numRandomVectors);

  size_t Q = numQuadraturePoints;
  lambdaNodes.reserve(Q);
  deltaN.reserve(Q);

  for (size_t k = 1; k <= Q; ++k) {
    double theta_k =
      (static_cast<double>(k) - 0.5) * M_PI / static_cast<double>(Q);
    double x_k = std::cos(theta_k);
    double lam_k = c + w * x_k;
    if (lam_k <= 1e-4)
      continue;

    double sum = mu[0];
    for (size_t m = 1; m < numMoments; ++m) {
      sum += 2.0 * mu[static_cast<Eigen::Index>(m)] *
             std::cos(static_cast<double>(m) * theta_k);
    }
    double dn =
      (static_cast<double>(n) / static_cast<double>(Q)) * std::max(0.0, sum);
    lambdaNodes.push_back(lam_k);
    deltaN.push_back(dn);
  }
}

Eigen::ArrayXd
NormalModeAnalyzer::computeDensityOfStatesKPM(const Eigen::ArrayXd& lambdaGrid,
                                              size_t numMoments,
                                              size_t numRandomVectors) const
{
  Eigen::Index n = this->assembledConnectivityMatrix.rows();
  if (n == 0)
    return Eigen::ArrayXd::Zero(lambdaGrid.size());

  auto [c, w] = getSpectralScalingParameters(0.01);
  Eigen::VectorXd mu = computeChebyshevMomentsKPM(numMoments, numRandomVectors);

  Eigen::ArrayXd dos = Eigen::ArrayXd::Zero(lambdaGrid.size());
  for (Eigen::Index i = 0; i < lambdaGrid.size(); ++i) {
    double lam = lambdaGrid[i];
    double x = (lam - c) / w;
    if (x <= -1.0 || x >= 1.0)
      continue;

    double t_prev = 1.0;
    double t_curr = x;
    double sum = mu[0];
    if (numMoments > 1) {
      sum += 2.0 * mu[1] * t_curr;
      for (size_t m = 1; m < numMoments - 1; ++m) {
        double t_next = 2.0 * x * t_curr - t_prev;
        sum += 2.0 * mu[static_cast<Eigen::Index>(m + 1)] * t_next;
        t_prev = t_curr;
        t_curr = t_next;
      }
    }
    double factor = 1.0 / (w * M_PI * std::sqrt(1.0 - x * x));
    dos[i] = static_cast<double>(n) * factor * std::max(0.0, sum);
  }
  return dos;
}

Eigen::ArrayXd
NormalModeAnalyzer::evaluateStressAutocorrelationKPM(
  const Eigen::ArrayXd& t,
  size_t numMoments,
  size_t numRandomVectors,
  size_t numQuadraturePoints) const
{
  Eigen::ArrayXd result = Eigen::ArrayXd::Zero(t.size());
  std::vector<double> lambdaNodes;
  std::vector<double> deltaN;
  computeEffectiveSpectrumKPM(
    numMoments, numRandomVectors, numQuadraturePoints, lambdaNodes, deltaN);

  for (size_t k = 0; k < lambdaNodes.size(); ++k) {
    result += deltaN[k] * (-2.0 * lambdaNodes[k] * t).exp();
  }
  return result;
}

Eigen::ArrayXd
NormalModeAnalyzer::evaluateStorageModulusKPM(const Eigen::ArrayXd& omega,
                                              size_t numMoments,
                                              size_t numRandomVectors,
                                              size_t numQuadraturePoints) const
{
  Eigen::ArrayXd result = Eigen::ArrayXd::Zero(omega.size());
  std::vector<double> lambdaNodes;
  std::vector<double> deltaN;
  computeEffectiveSpectrumKPM(
    numMoments, numRandomVectors, numQuadraturePoints, lambdaNodes, deltaN);

  for (size_t k = 0; k < lambdaNodes.size(); ++k) {
    Eigen::ArrayXd term = (omega / (2.0 * lambdaNodes[k])).square();
    result += deltaN[k] * (term / (1.0 + term));
  }
  return result;
}

Eigen::ArrayXd
NormalModeAnalyzer::evaluateLossModulusKPM(const Eigen::ArrayXd& omega,
                                           size_t numMoments,
                                           size_t numRandomVectors,
                                           size_t numQuadraturePoints) const
{
  Eigen::ArrayXd result = Eigen::ArrayXd::Zero(omega.size());
  std::vector<double> lambdaNodes;
  std::vector<double> deltaN;
  computeEffectiveSpectrumKPM(
    numMoments, numRandomVectors, numQuadraturePoints, lambdaNodes, deltaN);

  for (size_t k = 0; k < lambdaNodes.size(); ++k) {
    Eigen::ArrayXd term = (omega / (2.0 * lambdaNodes[k]));
    result += deltaN[k] * (term / (1.0 + term.square()));
  }
  return result;
}
}
