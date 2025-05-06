#ifndef NORMAL_MODE_ANALYZER_H
#define NORMAL_MODE_ANALYZER_H

#include "../utils/CerealUtils.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifdef CEREALIZABLE
#include <cereal/access.hpp>
#endif
#include <vector>

namespace pylimer_tools::calc {
class NormalModeAnalyzer
{
public:
  NormalModeAnalyzer(const std::vector<size_t> springFrom,
                     const std::vector<size_t> springTo);

  void findSparseEigenvalues(const size_t nrOfEigenvalues,
                             const bool includeEigenvectors = false);

  void computeAllEigenvalues(const bool includeEigenvectors = false);

  Eigen::VectorXd getEigenvalues() const;
  void setEigenvalues(Eigen::VectorXd e);

  Eigen::MatrixXd getEigenvectors() const;
  void setEigenvectors(Eigen::MatrixXd e);

  Eigen::SparseMatrix<double> getAssembledConnectivityMatrix() const;

  Eigen::ArrayXd evaluateStressAutocorrelation(const Eigen::ArrayXd& t) const;

  Eigen::ArrayXd evaluateStorageModulus(const Eigen::ArrayXd& omega) const;

  Eigen::ArrayXd evaluateLossModulus(const Eigen::ArrayXd& omega) const;

  size_t getNrOfSolubleClusters() const;

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
