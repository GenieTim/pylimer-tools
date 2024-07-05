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

namespace pylimer_tools {
namespace calc {
  NormalModeAnalyzer::NormalModeAnalyzer(const std::vector<size_t> springFrom,
                                         const std::vector<size_t> springTo)
  {
    INVALIDARG_EXP_IFN(springFrom.size() == springTo.size(),
                       "springFrom and springTo must have the same size");
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
      // triplets will be summed up -> we can use the same indices multiple
      // times
      triplets.push_back(
        Eigen::Triplet<double>(springFrom[i], springTo[i], -1.0));
      triplets.push_back(
        Eigen::Triplet<double>(springTo[i], springFrom[i], -1.0));
      // it is a bit more efficient to sum up the diagonal elements manually
      diagonal[springFrom[i]] += 1.0;
      diagonal[springTo[i]] += 1.0;
      // some sanity checks
      assert(springFrom[i] < nRows);
      assert(springTo[i] < nRows);
    }
    // the summed up diagonal can be translated to triplets as well
    for (size_t i = 0; i < diagonal.size(); i++) {
      if (diagonal[i] != 0.) {
        triplets.push_back(Eigen::Triplet<double>(i, i, diagonal[i]));
      }
    }

    // finally, translate into the sparse matrix format
    this->assembledConnectivityMatrix =
      Eigen::SparseMatrix<double>(nRows, nRows);
    this->assembledConnectivityMatrix.setFromTriplets(triplets.begin(),
                                                      triplets.end());
  };

  void NormalModeAnalyzer::findSparseEigenvalues(size_t nrOfEigenvalues)
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
    }
  };

  void NormalModeAnalyzer::computeAllEigenvalues(bool includeEigenvectors)
  {
    Eigen::MatrixXd assembledConnectivityMatrixDense =
      Eigen::MatrixXd(this->assembledConnectivityMatrix);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver =
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(
        assembledConnectivityMatrixDense,
        includeEigenvectors ? Eigen::DecompositionOptions::ComputeEigenvectors
                            : Eigen::DecompositionOptions::EigenvaluesOnly);

    this->setEigenvalues(solver.eigenvalues());
    if (includeEigenvectors) {
      this->setEigenvectors(solver.eigenvectors());
    }
  }

  void NormalModeAnalyzer::requireEigenvaluesComputation() const
  {
    RUNTIME_EXP_IFN(this->isEigenvaluesComputed,
                    "Eigenvalues have not been computed yet");
  }

  Eigen::VectorXd NormalModeAnalyzer::getEigenvalues() const
  {
    this->requireEigenvaluesComputation();
    return this->eigenvalues;
  }

  void NormalModeAnalyzer::setEigenvalues(Eigen::VectorXd eigenvalues)
  {
    this->eigenvalues = eigenvalues;
    this->isEigenvaluesComputed = true;
    this->clusterCount = this->countSolubleClusters();

    // TODO: possibly validate vector dimensions
  }

  void NormalModeAnalyzer::requireEigenvectorsComputation() const
  {
    RUNTIME_EXP_IFN(this->isEigenvectorsComputed,
                    "Eigenvectors have not been computed yet");
  }

  Eigen::MatrixXd NormalModeAnalyzer::getEigenvectors() const
  {
    this->requireEigenvectorsComputation();
    return this->eigenvectors;
  }

  void NormalModeAnalyzer::setEigenvectors(Eigen::MatrixXd eigenvectors)
  {
    this->eigenvectors = eigenvectors;
    this->isEigenvectorsComputed = true;
    // TODO: possibly validate matrix dimensions
  }

  size_t NormalModeAnalyzer::countSolubleClusters() const
  {
    this->requireEigenvaluesComputation();
    return (this->eigenvalues.array() == 0.0).count();
  }

  size_t NormalModeAnalyzer::getNrOfSolubleClusters() const
  {
    return this->clusterCount;
  }

  Eigen::ArrayXd NormalModeAnalyzer::evaluateStressAutocorrelation(
    const Eigen::ArrayXd& t) const
  {
    Eigen::ArrayXd result = Eigen::ArrayXd::Zero(t.size());
    for (size_t i = 0; i < this->eigenvalues.size() - this->clusterCount; ++i) {
      result += (-2. * this->eigenvalues[i] * t).exp();
    }
    return result;
  };

  Eigen::ArrayXd NormalModeAnalyzer::evaluateStorageModulus(
    const Eigen::ArrayXd& omega) const
  {
    Eigen::ArrayXd result = Eigen::ArrayXd::Zero(omega.size());
    for (size_t i = 0; i < this->eigenvalues.size() - this->clusterCount; ++i) {
      result += (omega / (2. * this->eigenvalues[i])).square() /
                (1. + (omega / (2. * this->eigenvalues[i])).square());
    }
    return result;
  };

  Eigen::ArrayXd NormalModeAnalyzer::evaluateLossModulus(
    const Eigen::ArrayXd& omega) const
  {
    Eigen::ArrayXd result = Eigen::ArrayXd::Zero(omega.size());
    for (size_t i = 0; i < this->eigenvalues.size() - this->clusterCount; ++i) {
      result += (omega / (2. * this->eigenvalues[i])) /
                (1. + (omega / (2. * this->eigenvalues[i])).square());
    }
    return result;
  };
}
}
