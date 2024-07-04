#ifndef NORMAL_MODE_ANALYZER_H
#define NORMAL_MODE_ANALYZER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace pylimer_tools {
namespace calc {
  class NormalModeAnalyzer
  {
  public:
    NormalModeAnalyzer(const std::vector<size_t> springFrom,
                       const std::vector<size_t> springTo);

    void findSparseEigenvalues(size_t nrOfEigenvalues);

    void computeAllEigenvalues();

    Eigen::VectorXd getEigenvalues() const;
    void setEigenvalues(Eigen::VectorXd e);

    Eigen::MatrixXd getEigenvectors() const;
    void setEigenvectors(Eigen::MatrixXd e);

    Eigen::ArrayXd evaluateStressAutocorrelation(const Eigen::ArrayXd &t) const;

    Eigen::ArrayXd evaluateStorageModulus(const Eigen::ArrayXd &omega) const;

    Eigen::ArrayXd evaluateLossModulus(const Eigen::ArrayXd &omega) const;

    size_t getNrOfSolubleClusters() const;

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
  };
}
}

#endif
