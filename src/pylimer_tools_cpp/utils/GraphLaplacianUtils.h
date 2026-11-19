#ifndef GRAPH_LAPLACIAN_UTILS_H
#define GRAPH_LAPLACIAN_UTILS_H

#include "../entities/Universe.h"
#include "./utilityMacros.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace pylimer_tools::utils {

/**
 * @brief Assemble a symmetric positive semidefinite sparse graph Laplacian
 *        matrix L from edge lists (springFrom, springTo) and optional weights.
 *
 * Self-loops (u == v) are skipped from L and the returned edge/weight totals,
 * though all specified weights are validated to be finite and positive.
 *
 * @param springFrom Source node indices of edges
 * @param springTo Target node indices of edges
 * @param weights Optional positive edge weights (empty defaults to 1.0)
 * @param numNodes Optional total number of nodes (deduced from edges if 0)
 * @param outNumNodes Optional output receiving resolved number of nodes
 * @param outNumEdges Optional output receiving number of non-self-loop edges
 * @param outTotalWeight Optional output receiving sum of non-self-loop weights
 * @return Compressed symmetric sparse graph Laplacian L
 */
inline Eigen::SparseMatrix<double>
assembleGraphLaplacian(const std::vector<size_t>& springFrom,
                       const std::vector<size_t>& springTo,
                       const std::vector<double>& weights = {},
                       const size_t numNodes = 0,
                       size_t* outNumNodes = nullptr,
                       size_t* outNumEdges = nullptr,
                       double* outTotalWeight = nullptr)
{
  INVALIDARG_EXP_IFN(springFrom.size() == springTo.size(),
                     "springFrom and springTo must have the same size");
  INVALIDARG_EXP_IFN(!springFrom.empty() || numNodes > 0,
                     "Graph cannot be completely empty");
  INVALIDARG_EXP_IFN(
    weights.empty() || weights.size() == springFrom.size(),
    "weights must be empty or have the same size as springFrom");

  for (size_t e = 0; e < weights.size(); ++e) {
    INVALIDARG_EXP_IFN(
      std::isfinite(weights[e]) && weights[e] > 0.0,
      "Edge weight/spring constant must be finite and positive");
  }

  size_t resolvedNumNodes = numNodes;
  if (!springFrom.empty()) {
    const size_t maxIdx =
      std::max(*std::max_element(springFrom.begin(), springFrom.end()),
               *std::max_element(springTo.begin(), springTo.end()));
    if (numNodes > 0) {
      INVALIDARG_EXP_IFN(
        numNodes > maxIdx,
        "numNodes must be greater than the maximum node index in springs (" +
          std::to_string(maxIdx) + "), but got " + std::to_string(numNodes));
    } else {
      INVALIDARG_EXP_IFN(maxIdx <
                           static_cast<size_t>(std::numeric_limits<int>::max()),
                         "Node index exceeds Eigen sparse index limit");
      resolvedNumNodes = maxIdx + 1;
    }
  }

  INVALIDARG_EXP_IFN(resolvedNumNodes <=
                       static_cast<size_t>(std::numeric_limits<int>::max()),
                     "Number of nodes exceeds Eigen sparse index limit");

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(springFrom.size() * 2 + resolvedNumNodes);
  Eigen::VectorXd diag =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(resolvedNumNodes));

  size_t edgeCount = 0;
  double totalWeight = 0.0;

  for (size_t e = 0; e < springFrom.size(); ++e) {
    const size_t u = springFrom[e];
    const size_t v = springTo[e];
    if (u == v) {
      continue;
    }
    const double w = weights.empty() ? 1.0 : weights[e];
    triplets.emplace_back(static_cast<int>(u), static_cast<int>(v), -w);
    triplets.emplace_back(static_cast<int>(v), static_cast<int>(u), -w);
    diag(static_cast<Eigen::Index>(u)) += w;
    diag(static_cast<Eigen::Index>(v)) += w;
    totalWeight += w;
    ++edgeCount;
  }

  for (Eigen::Index i = 0; i < diag.size(); ++i) {
    if (diag(i) != 0.0) {
      triplets.emplace_back(static_cast<int>(i), static_cast<int>(i), diag(i));
    }
  }

  Eigen::SparseMatrix<double> L(static_cast<Eigen::Index>(resolvedNumNodes),
                                static_cast<Eigen::Index>(resolvedNumNodes));
  L.setFromTriplets(triplets.begin(), triplets.end());
  L.makeCompressed();

  if (outNumNodes != nullptr) {
    *outNumNodes = resolvedNumNodes;
  }
  if (outNumEdges != nullptr) {
    *outNumEdges = edgeCount;
  }
  if (outTotalWeight != nullptr) {
    *outTotalWeight = totalWeight;
  }

  return L;
}

/**
 * @brief Validate that a sparse matrix L is a valid symmetric graph Laplacian
 *        with non-negative diagonals, non-positive off-diagonals, and
 *        zero column/row sums.
 *
 * @param L Sparse matrix to validate
 * @param outNumEdges Optional output receiving number of off-diagonal edges
 * @param outTotalWeight Optional output receiving 0.5 * Tr(L)
 */
inline void
validateGraphLaplacian(const Eigen::SparseMatrix<double>& L,
                       size_t* outNumEdges = nullptr,
                       double* outTotalWeight = nullptr)
{
  INVALIDARG_EXP_IFN(L.rows() == L.cols(), "Laplacian matrix L must be square");
  INVALIDARG_EXP_IFN(L.rows() > 0, "Laplacian matrix L cannot be empty");

  double diagSum = 0.0;
  size_t offDiagCount = 0;
  for (int k = 0; k < L.outerSize(); ++k) {
    double colSum = 0.0;
    double colDiag = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
      const double val = it.value();
      INVALIDARG_EXP_IFN(std::isfinite(val),
                         "Laplacian matrix L must contain only finite entries");
      colSum += val;
      if (it.row() == it.col()) {
        INVALIDARG_EXP_IFN(
          val >= 0.0,
          "Laplacian matrix L must have non-negative diagonal entries");
        diagSum += val;
        colDiag += val;
      }
    }
    const double posTol = 1e-12 * std::max(1.0, colDiag);
    const double edgeTol = 1e-14 * std::max(colDiag, 0.0);
    for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
      if (it.row() != it.col()) {
        const double val = it.value();
        INVALIDARG_EXP_IFN(
          val <= posTol,
          "Laplacian matrix L must have non-positive off-diagonal entries");
        if (it.row() < it.col() && val < 0.0 && std::abs(val) > edgeTol) {
          ++offDiagCount;
        }
      }
    }
    INVALIDARG_EXP_IFN(std::abs(colSum) <= 1e-9 * std::max(1.0, colDiag),
                       "Laplacian matrix L columns/rows must sum to zero");
  }

  const Eigen::SparseMatrix<double> diff =
    L - Eigen::SparseMatrix<double>(L.transpose());
  INVALIDARG_EXP_IFN(diff.norm() <= 1e-10 * std::max(1.0, L.norm()),
                     "Laplacian matrix L must be symmetric");

  if (outNumEdges != nullptr) {
    *outNumEdges = offDiagCount;
  }
  if (outTotalWeight != nullptr) {
    *outTotalWeight = 0.5 * diagSum;
  }
}

/**
 * @brief Compute connected components of a graph from its sparse Laplacian L
 *        using Union-Find directly over L's non-zeros.
 *
 * Components are ordered by first-seen vertex index (0..n-1), and vertices
 * within each component are in ascending order so comp[0] is the smallest node
 * index in that component.
 *
 * @param L Sparse graph Laplacian matrix
 * @param outComponents Output list of components (each a sorted list of node
 *                      indices)
 * @param outNodeComponentId Optional output mapping each node to its component
 *                           index
 * @param outComponentTotalWeight Optional output receiving total edge weight
 *                                (0.5 * sum of diagonal entries) per component
 */
inline void
computeConnectedComponents(
  const Eigen::SparseMatrix<double>& L,
  std::vector<std::vector<size_t>>& outComponents,
  std::vector<size_t>* outNodeComponentId = nullptr,
  std::vector<double>* outComponentTotalWeight = nullptr)
{
  const size_t n = static_cast<size_t>(L.rows());
  outComponents.clear();
  if (outNodeComponentId != nullptr) {
    outNodeComponentId->assign(n, 0);
  }
  if (outComponentTotalWeight != nullptr) {
    outComponentTotalWeight->clear();
  }
  if (n == 0) {
    return;
  }

  std::vector<size_t> parent(n);
  std::iota(parent.begin(), parent.end(), 0);
  auto findRoot = [&parent](size_t i) -> size_t {
    while (parent[i] != i) {
      parent[i] = parent[parent[i]];
      i = parent[i];
    }
    return i;
  };

  std::vector<double> diagWeights;
  if (outComponentTotalWeight != nullptr) {
    diagWeights.assign(n, 0.0);
  }

  for (int k = 0; k < L.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
      const Eigen::Index r = it.row();
      const Eigen::Index c = it.col();
      if (r == c) {
        if (outComponentTotalWeight != nullptr) {
          diagWeights[static_cast<size_t>(r)] += it.value();
        }
      } else if (it.value() != 0.0) {
        const size_t ru = findRoot(static_cast<size_t>(r));
        const size_t rv = findRoot(static_cast<size_t>(c));
        if (ru != rv) {
          parent[ru] = rv;
        }
      }
    }
  }

  std::unordered_map<size_t, size_t> rootToCompId;
  rootToCompId.reserve(16);
  for (size_t i = 0; i < n; ++i) {
    const size_t root = findRoot(i);
    const std::unordered_map<size_t, size_t>::const_iterator it =
      rootToCompId.find(root);
    size_t compId = 0;
    if (it == rootToCompId.end()) {
      compId = outComponents.size();
      rootToCompId[root] = compId;
      outComponents.emplace_back();
      if (outComponentTotalWeight != nullptr) {
        outComponentTotalWeight->push_back(0.0);
      }
    } else {
      compId = it->second;
    }
    outComponents[compId].push_back(i);
    if (outNodeComponentId != nullptr) {
      (*outNodeComponentId)[i] = compId;
    }
    if (outComponentTotalWeight != nullptr) {
      (*outComponentTotalWeight)[compId] += 0.5 * diagWeights[i];
    }
  }
}

/**
 * @brief Build edge lists for a linear chain of numNodes vertices.
 */
inline void
buildLinearChainEdges(const size_t numNodes,
                      std::vector<size_t>& outFrom,
                      std::vector<size_t>& outTo)
{
  INVALIDARG_EXP_IFN(numNodes > 0, "numNodes must be > 0");
  if (numNodes == 1) {
    outFrom.clear();
    outTo.clear();
    return;
  }
  outFrom.resize(numNodes - 1);
  outTo.resize(numNodes - 1);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    outFrom[i] = i;
    outTo[i] = i + 1;
  }
}

/**
 * @brief Build edge lists for a ring graph of numNodes vertices.
 */
inline void
buildRingEdges(const size_t numNodes,
               std::vector<size_t>& outFrom,
               std::vector<size_t>& outTo)
{
  INVALIDARG_EXP_IFN(numNodes >= 3, "Ring graph requires at least 3 nodes");
  outFrom.resize(numNodes);
  outTo.resize(numNodes);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    outFrom[i] = i;
    outTo[i] = i + 1;
  }
  outFrom[numNodes - 1] = numNodes - 1;
  outTo[numNodes - 1] = 0;
}

/**
 * @brief Build edge lists for a complete graph K_{numNodes}.
 */
inline void
buildCompleteGraphEdges(const size_t numNodes,
                        std::vector<size_t>& outFrom,
                        std::vector<size_t>& outTo)
{
  INVALIDARG_EXP_IFN(numNodes >= 2, "Complete graph requires at least 2 nodes");
  const size_t nEdges = (numNodes * (numNodes - 1)) / 2;
  outFrom.clear();
  outTo.clear();
  outFrom.reserve(nEdges);
  outTo.reserve(nEdges);
  for (size_t i = 0; i < numNodes; ++i) {
    for (size_t j = i + 1; j < numNodes; ++j) {
      outFrom.push_back(i);
      outTo.push_back(j);
    }
  }
}

/**
 * @brief Build edge lists for a star graph with numArms arms of armLength
 *        edges each.
 */
inline void
buildStarEdges(const size_t numArms,
               const size_t armLength,
               std::vector<size_t>& outFrom,
               std::vector<size_t>& outTo,
               size_t& outNumNodes)
{
  INVALIDARG_EXP_IFN(numArms >= 2, "Star graph requires at least 2 arms");
  INVALIDARG_EXP_IFN(armLength >= 1, "armLength must be >= 1");
  outNumNodes = 1 + numArms * armLength;
  outFrom.clear();
  outTo.clear();
  outFrom.reserve(numArms * armLength);
  outTo.reserve(numArms * armLength);

  size_t nextNode = 1;
  for (size_t a = 0; a < numArms; ++a) {
    outFrom.push_back(0);
    outTo.push_back(nextNode);
    for (size_t l = 1; l < armLength; ++l) {
      outFrom.push_back(nextNode + l - 1);
      outTo.push_back(nextNode + l);
    }
    nextNode += armLength;
  }
}

/**
 * @brief Extract vertex-indexed edge lists and node count from a Universe.
 */
inline void
extractUniverseEdges(const pylimer_tools::entities::Universe& universe,
                     std::vector<size_t>& outFrom,
                     std::vector<size_t>& outTo,
                     size_t& outNumNodes)
{
  outNumNodes = universe.getNrOfAtoms();
  const std::map<std::string, std::vector<long int>> edges =
    universe.getEdges();
  const std::vector<long int>& fromRaw = edges.at("edge_from");
  const std::vector<long int>& toRaw = edges.at("edge_to");

  outFrom.resize(fromRaw.size());
  outTo.resize(toRaw.size());
  for (size_t i = 0; i < fromRaw.size(); ++i) {
    outFrom[i] = static_cast<size_t>(fromRaw[i]);
    outTo[i] = static_cast<size_t>(toRaw[i]);
  }
}

} // namespace pylimer_tools::utils

#endif // GRAPH_LAPLACIAN_UTILS_H
