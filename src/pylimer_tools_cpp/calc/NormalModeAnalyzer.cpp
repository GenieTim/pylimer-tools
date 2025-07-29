#include "./NormalModeAnalyzer.h"

#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <algorithm>
#include <iostream>
#include <vector>
#ifdef EIGEN_USE_LAPACKE
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
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
  Spectra::SparseGenMatProd<double> op =
    Spectra::SparseGenMatProd<double>(this->assembledConnectivityMatrix);
  Spectra::SymEigsSolver<Spectra::SparseGenMatProd<double>> eigs(
    op, nrOfEigenvalues, 2 * nrOfEigenvalues);
  // Initialize and compute
  eigs.init();

  int nconv = eigs.compute(Spectra::SortRule::SmallestAlge);

  // Retrieve results
  if (eigs.info() == Spectra::CompInfo::Successful) {
    this->setEigenvalues(eigs.eigenvalues());
    if (includeEigenvectors) {
      this->setEigenvectors(eigs.eigenvectors());
    }
  }
};

void
NormalModeAnalyzer::computeAllEigenvalues(const bool includeEigenvectors)
{
  RUNTIME_EXP_IFN(this->assembledConnectivityMatrix.rows() ==
                    this->assembledConnectivityMatrix.cols(),
                  "Expected square matrix");

  Eigen::Index nRows = this->assembledConnectivityMatrix.rows();
  if (nRows == 0) {
    return;
  }

  Eigen::MatrixXd assembledConnectivityMatrixDense =
    Eigen::MatrixXd(this->assembledConnectivityMatrix);
  RUNTIME_EXP_IFN(assembledConnectivityMatrixDense.rows() == nRows,
                  "Expected square matrix also after conversion to dense.");
  RUNTIME_EXP_IFN(assembledConnectivityMatrixDense.cols() == nRows,
                  "Expected square matrix also after conversion to dense.");

#ifndef EIGEN_USE_LAPACKE
  std::cerr << "Eigen LAPACK is not available, using Eigen's "
               "SelfAdjointEigenSolver instead. Expect reduced performance."
            << std::endl;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver =
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(
      assembledConnectivityMatrixDense,
      includeEigenvectors ? Eigen::DecompositionOptions::ComputeEigenvectors
                          : Eigen::DecompositionOptions::EigenvaluesOnly);

  this->setEigenvalues(solver.eigenvalues());
  if (includeEigenvectors) {
    this->setEigenvectors(solver.eigenvectors());
  }
#else
  // see also/alternative: igraph_lapack_dsyevr()
  Eigen::VectorXd eigenvaluesMemory = Eigen::VectorXd::Zero(nRows);
  Eigen::MatrixXd eigenvectorsMemory = includeEigenvectors
                                         ? Eigen::MatrixXd::Zero(nRows, nRows)
                                         : Eigen::MatrixXd::Zero(nRows, 3);
  lapack_int ldz = static_cast<lapack_int>(eigenvectorsMemory.cols());
  lapack_int il = 0;
  lapack_int iu = 0;
  double abstol = 1e-10;
  lapack_int M = 0;
  lapack_int size = static_cast<lapack_int>(nRows);

  Eigen::Vector<lapack_int, Eigen::Dynamic> support =
    Eigen::Vector<lapack_int, Eigen::Dynamic>::Zero(nRows * 2);

  int info = LAPACKE_dsyevr(Eigen::MatrixXd::IsRowMajor
                              ? LAPACK_ROW_MAJOR
                              : LAPACK_COL_MAJOR,            // matrix_order
                            includeEigenvectors ? 'V' : 'N', // JOBZ
                            'A',                             // RANGE
                            'U',                             // UPLO
                            size,                            // N
                            assembledConnectivityMatrixDense.data(), // A
                            size,                                    // LDA
                            0,                                       // VL
                            0,                                       // VU
                            il,                                      // IL
                            iu,                                      // IU
                            abstol,                                  // ABSTOL
                            &M, // M, total number of eigenvalues found
                            eigenvaluesMemory.data(),  // W
                            eigenvectorsMemory.data(), // Z
                            ldz,                       // LDZ
                            support.data()             // ISUPPZ
  );

  RUNTIME_EXP_IFN(info == 0,
                  "Error in LAPACK_dsyevr. Exit code " + std::to_string(info) +
                    ".");
  this->setEigenvalues(eigenvaluesMemory);
  this->setEigenvectors(eigenvectorsMemory);
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
  this
    ->eigenvalues(this->eigenvalues.array() > -NORMAL_MODE_ZERO_TOLERANCE &&
                  this->eigenvalues.array() < NORMAL_MODE_ZERO_TOLERANCE)
    .setZero();
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
}
