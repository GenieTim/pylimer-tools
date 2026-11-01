#ifndef RESISTANCE_DISTANCE_ANALYSIS_H
#define RESISTANCE_DISTANCE_ANALYSIS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstddef>
#include <utility>
#include <vector>

namespace pylimer_tools::topo {

/**
 * @brief Analysis of topological resistance distances, Kirchhoff index,
 *        and random-walk commute times on graph Laplacians.
 */
class ResistanceDistanceAnalysis
{
public:
  /**
   * @brief Construct from edge lists (springFrom, springTo), optional edge
   * weights, and optional number of nodes.
   *
   * @param springFrom Source node indices of edges
   * @param springTo Target node indices of edges
   * @param weights Optional weights (spring constants) of edges (defaults
   * to 1.0)
   * @param numNodes Optional total number of nodes (deduced if 0)
   */
  ResistanceDistanceAnalysis(const std::vector<size_t>& springFrom,
                             const std::vector<size_t>& springTo,
                             const std::vector<double>& weights = {},
                             size_t numNodes = 0);

  /**
   * @brief Construct from pre-assembled sparse scalar graph Laplacian matrix L.
   *
   * @param L Symmetric positive semidefinite graph Laplacian matrix
   */
  explicit ResistanceDistanceAnalysis(const Eigen::SparseMatrix<double>& L);

  // Factory methods for canonical topologies
  static ResistanceDistanceAnalysis createLinearChain(size_t numNodes);
  static ResistanceDistanceAnalysis createRing(size_t numNodes);
  static ResistanceDistanceAnalysis createComplete(size_t numNodes);
  static ResistanceDistanceAnalysis createStar(size_t numArms,
                                               size_t armLength = 1);

  // Graph topology queries
  size_t getNrOfNodes() const { return numNodes_; }
  size_t getNrOfEdges() const { return numEdges_; }
  double getTotalEdgeWeight() const { return totalWeight_; }
  bool isConnected() const { return connectedComponents_.size() <= 1; }
  size_t getNrOfConnectedComponents() const
  {
    return connectedComponents_.size();
  }
  const std::vector<std::vector<size_t>>& getConnectedComponents() const
  {
    return connectedComponents_;
  }

  const Eigen::SparseMatrix<double>& getLaplacian() const { return L_; }
  const Eigen::SparseMatrix<double>& getPinnedLaplacian() const
  {
    ensureFactorized();
    return L_pinned_;
  }

  /**
   * @brief Compute effective resistance between a single pair of nodes (i, j).
   *
   * Solves L v = (e_i - e_j) with a reference node pinned to 0.
   *
   * @param i First node index
   * @param j Second node index
   * @return Effective resistance R_ij >= 0. Returns infinity if i and j belong
   *         to different connected components.
   */
  double computeEffectiveResistance(size_t i, size_t j) const;

  /**
   * @brief Compute effective resistance for a batch of B pairs simultaneously
   *        using a single Multi-RHS linear solve:
   *        L [v_1, ..., v_B] = [b_1, ..., b_B].
   *
   * @param pairs Vector of (i, j) node index pairs
   * @return Vector of effective resistances for each query pair
   */
  Eigen::VectorXd computeEffectiveResistanceBatched(
    const std::vector<std::pair<size_t, size_t>>& pairs) const;

  /**
   * @brief Compute all-pairs effective resistance matrix R of size N x N.
   *
   * Solves L X = I with pinned reference nodes, giving all pairwise resistances
   * R_ij = X_ii + X_jj - 2 * X_ij.
   *
   * @return Symmetric N x N matrix of effective resistances
   */
  Eigen::MatrixXd computeAllPairsEffectiveResistance() const;

  /**
   * @brief Compute total Kirchhoff index K(G) = sum_{i < j} R_ij.
   *
   * For connected graphs, K(G) = N * Tr(L^+) = N * sum_{k=2}^N (1 / lambda_k).
   * Returns infinity if the graph is disconnected.
   */
  double computeKirchhoffIndex() const;

  /**
   * @brief Compute Kirchhoff index from Laplacian non-zero eigenvalues for
   * verification.
   */
  double computeKirchhoffIndexSpectral() const;

  /**
   * @brief Compute Moore-Penrose pseudo-inverse of the Laplacian L^+ (size N x
   * N).
   *
   * Satisfies L L^+ L = L and L^+ 1 = 0.
   */
  Eigen::MatrixXd computeMoorePenrosePseudoInverse() const;

  /**
   * @brief Compute random-walk commute time between node i and node j:
   *        C_ij = 2 * M * R_ij, where 2 * M = Tr(L) is the total graph volume.
   */
  double computeCommuteTime(size_t i, size_t j) const;

  /**
   * @brief Compute all-pairs commute time matrix: C = 2 * M * R.
   */
  Eigen::MatrixXd computeAllPairsCommuteTime() const;

  /**
   * @brief Compute total effective resistance from each node to all other
   * nodes: S_i = sum_{j} R_ij.
   */
  Eigen::VectorXd computeTotalResistancePerNode() const;

  /**
   * @brief Compute resistance centrality for each node: (N - 1) / sum_{j != i}
   * R_ij. Higher values indicate nodes closer to all other nodes in resistance
   * distance.
   */
  Eigen::VectorXd computeResistanceCentrality() const;

private:
  size_t numNodes_{ 0 };
  size_t numEdges_{ 0 };
  double totalWeight_{ 0.0 };
  Eigen::SparseMatrix<double> L_{};
  Eigen::SparseMatrix<double> L_pinned_{};

  std::vector<std::vector<size_t>> connectedComponents_{};
  std::vector<size_t> nodeComponentId_{};
  std::vector<size_t> referenceNodePerComponent_{};

  mutable Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_{};
  mutable bool solverFactorized_{ false };

  void initializeFromEdges(const std::vector<size_t>& springFrom,
                           const std::vector<size_t>& springTo,
                           const std::vector<double>& weights,
                           size_t numNodes);
  void findComponentsAndSetupPinnedSystem();
  void ensureFactorized() const;
};

} // namespace pylimer_tools::topo

#endif // RESISTANCE_DISTANCE_ANALYSIS_H
