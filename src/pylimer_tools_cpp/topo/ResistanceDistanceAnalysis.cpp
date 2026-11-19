#include "ResistanceDistanceAnalysis.h"

#include "../entities/Universe.h"
#include "../utils/GraphLaplacianUtils.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::topo {

ResistanceDistanceAnalysis::ResistanceDistanceAnalysis(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  const std::vector<double>& weights,
  const size_t numNodes)
{
  initializeFromEdges(springFrom, springTo, weights, numNodes);
}

ResistanceDistanceAnalysis::ResistanceDistanceAnalysis(
  const Eigen::SparseMatrix<double>& L)
  : numNodes_(static_cast<size_t>(L.rows()))
  , numEdges_(0)
  , totalWeight_(0.0)
  , L_(L)
  , L_pinned_()
  , connectedComponents_()
  , nodeComponentId_()
  , referenceNodePerComponent_()
  , componentTotalWeight_()
  , solver_()
  , solverFactorized_(false)
{
  L_.makeCompressed();
  pylimer_tools::utils::validateGraphLaplacian(L_, &numEdges_, &totalWeight_);
  findComponentsAndSetupPinnedSystem();
}

ResistanceDistanceAnalysis::ResistanceDistanceAnalysis(
  const pylimer_tools::entities::Universe& universe)
{
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  size_t nNodes = 0;
  pylimer_tools::utils::extractUniverseEdges(
    universe, springFrom, springTo, nNodes);
  initializeFromEdges(springFrom, springTo, {}, nNodes);
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::fromMatrix(const Eigen::SparseMatrix<double>& L)
{
  return ResistanceDistanceAnalysis(L);
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::fromUniverse(
  const pylimer_tools::entities::Universe& universe)
{
  return ResistanceDistanceAnalysis(universe);
}

void
ResistanceDistanceAnalysis::initializeFromEdges(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  const std::vector<double>& weights,
  const size_t numNodes)
{
  L_ = pylimer_tools::utils::assembleGraphLaplacian(springFrom,
                                                    springTo,
                                                    weights,
                                                    numNodes,
                                                    &numNodes_,
                                                    &numEdges_,
                                                    &totalWeight_);
  findComponentsAndSetupPinnedSystem();
}

void
ResistanceDistanceAnalysis::findComponentsAndSetupPinnedSystem()
{
  pylimer_tools::utils::computeConnectedComponents(
    L_, connectedComponents_, &nodeComponentId_, &componentTotalWeight_);

  referenceNodePerComponent_.clear();
  referenceNodePerComponent_.reserve(connectedComponents_.size());
  for (const std::vector<size_t>& comp : connectedComponents_) {
    referenceNodePerComponent_.push_back(comp[0]);
  }

  std::vector<bool> isRefNode(numNodes_, false);
  for (const size_t r : referenceNodePerComponent_) {
    isRefNode[r] = true;
  }

  std::vector<Eigen::Triplet<double>> pinnedTriplets;
  pinnedTriplets.reserve(static_cast<size_t>(L_.nonZeros()) +
                         referenceNodePerComponent_.size());

  for (int k = 0; k < L_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L_, k); it; ++it) {
      const size_t r = static_cast<size_t>(it.row());
      const size_t c = static_cast<size_t>(it.col());
      if (isRefNode[r] || isRefNode[c]) {
        // Zero out off-diagonals on reference rows and columns
        if (r == c && isRefNode[r]) {
          // Reference node diagonal will be set to 1.0 below
          continue;
        }
      } else {
        pinnedTriplets.emplace_back(
          static_cast<int>(r), static_cast<int>(c), it.value());
      }
    }
  }

  // Ensure reference nodes always have diagonal 1.0 even if isolated
  for (const size_t r : referenceNodePerComponent_) {
    pinnedTriplets.emplace_back(static_cast<int>(r), static_cast<int>(r), 1.0);
  }

  L_pinned_.resize(static_cast<Eigen::Index>(numNodes_),
                   static_cast<Eigen::Index>(numNodes_));
  L_pinned_.setFromTriplets(pinnedTriplets.begin(), pinnedTriplets.end());
  L_pinned_.makeCompressed();

  solverFactorized_ = false;
}

void
ResistanceDistanceAnalysis::ensureFactorized() const
{
  if (!solverFactorized_) {
    solver_.compute(L_pinned_);
    RUNTIME_EXP_IFN(solver_.info() == Eigen::Success,
                    "SimplicialLDLT factorization of pinned Laplacian failed");
    solverFactorized_ = true;
  }
}

double
ResistanceDistanceAnalysis::computeEffectiveResistance(const size_t i,
                                                       const size_t j) const
{
  INVALIDINDEX_EXP_IFN(i < numNodes_ && j < numNodes_,
                       "Node index out of range");

  if (i == j) {
    return 0.0;
  }

  if (nodeComponentId_[i] != nodeComponentId_[j]) {
    return std::numeric_limits<double>::infinity();
  }

  ensureFactorized();
  const size_t comp = nodeComponentId_[i];
  const size_t r = referenceNodePerComponent_[comp];

  Eigen::VectorXd b =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numNodes_));
  if (i != r) {
    b(static_cast<Eigen::Index>(i)) = 1.0;
  }
  if (j != r) {
    b(static_cast<Eigen::Index>(j)) = -1.0;
  }

  const Eigen::VectorXd v = solver_.solve(b);
  const double vi = (i == r) ? 0.0 : v(static_cast<Eigen::Index>(i));
  const double vj = (j == r) ? 0.0 : v(static_cast<Eigen::Index>(j));

  return std::max(0.0, vi - vj);
}

Eigen::VectorXd
ResistanceDistanceAnalysis::computeEffectiveResistanceBatched(
  const std::vector<std::pair<size_t, size_t>>& pairs,
  const size_t chunkSize) const
{
  INVALIDARG_EXP_IFN(chunkSize > 0, "chunkSize must be > 0");
  const size_t B = pairs.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(B));
  if (B == 0) {
    return result;
  }

  std::vector<size_t> activeIndices;
  activeIndices.reserve(B);

  for (size_t b = 0; b < B; ++b) {
    const size_t i = pairs[b].first;
    const size_t j = pairs[b].second;
    INVALIDINDEX_EXP_IFN(i < numNodes_ && j < numNodes_,
                         "Node index out of range in batched resistance pairs");

    if (i == j) {
      result(static_cast<Eigen::Index>(b)) = 0.0;
    } else if (nodeComponentId_[i] != nodeComponentId_[j]) {
      result(static_cast<Eigen::Index>(b)) =
        std::numeric_limits<double>::infinity();
    } else {
      activeIndices.push_back(b);
    }
  }

  if (activeIndices.empty()) {
    return result;
  }

  ensureFactorized();

  const Eigen::Index N = static_cast<Eigen::Index>(numNodes_);
  const size_t numActive = activeIndices.size();
  const size_t numChunks = (numActive + chunkSize - 1) / chunkSize;

#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic) if (numChunks > 1)
#endif
  for (size_t chunkIdx = 0; chunkIdx < numChunks; ++chunkIdx) {
    const size_t start = chunkIdx * chunkSize;
    const size_t curChunk = std::min(chunkSize, numActive - start);
    Eigen::MatrixXd Brhs =
      Eigen::MatrixXd::Zero(N, static_cast<Eigen::Index>(curChunk));

    for (size_t c = 0; c < curChunk; ++c) {
      const size_t b = activeIndices[start + c];
      const size_t i = pairs[b].first;
      const size_t j = pairs[b].second;
      const size_t comp = nodeComponentId_[i];
      const size_t r = referenceNodePerComponent_[comp];
      if (i != r) {
        Brhs(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(c)) += 1.0;
      }
      if (j != r) {
        Brhs(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(c)) -= 1.0;
      }
    }

    const Eigen::MatrixXd V = solver_.solve(Brhs);

    for (size_t c = 0; c < curChunk; ++c) {
      const size_t b = activeIndices[start + c];
      const size_t i = pairs[b].first;
      const size_t j = pairs[b].second;
      const size_t comp = nodeComponentId_[i];
      const size_t r = referenceNodePerComponent_[comp];
      const double vi = (i == r) ? 0.0
                                 : V(static_cast<Eigen::Index>(i),
                                     static_cast<Eigen::Index>(c));
      const double vj = (j == r) ? 0.0
                                 : V(static_cast<Eigen::Index>(j),
                                     static_cast<Eigen::Index>(c));
      result(static_cast<Eigen::Index>(b)) = std::max(0.0, vi - vj);
    }
  }

  return result;
}

Eigen::MatrixXd
ResistanceDistanceAnalysis::solveIdentityChunked(const size_t chunkSize) const
{
  INVALIDARG_EXP_IFN(chunkSize > 0, "chunkSize must be > 0");
  ensureFactorized();
  if (numNodes_ == 0) {
    return Eigen::MatrixXd(0, 0);
  }

  const Eigen::Index N = static_cast<Eigen::Index>(numNodes_);
  Eigen::MatrixXd X(N, N);
  const Eigen::Index chunkSizeIdx = static_cast<Eigen::Index>(chunkSize);
  const Eigen::Index numChunks = (N + chunkSizeIdx - 1) / chunkSizeIdx;

  // Level-3 BLAS chunked forward/back-substitutions:
  // Solving in chunks of B columns converts triangular solves from memory-bound
  // dtrsv (Level-2 BLAS) into compute-bound dtrsm (Level-3 BLAS), exploiting
  // AMX/Accelerate
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic) if (N > 256 && numChunks > 1)
#endif
  for (Eigen::Index chunkIdx = 0; chunkIdx < numChunks; ++chunkIdx) {
    const Eigen::Index colStart = chunkIdx * chunkSizeIdx;
    const Eigen::Index curChunk = std::min(chunkSizeIdx, N - colStart);
    Eigen::MatrixXd I_chunk = Eigen::MatrixXd::Zero(N, curChunk);
    for (Eigen::Index c = 0; c < curChunk; ++c) {
      I_chunk(colStart + c, c) = 1.0;
    }
    X.middleCols(colStart, curChunk).noalias() = solver_.solve(I_chunk);
  }

  return X;
}

Eigen::MatrixXd
ResistanceDistanceAnalysis::computeAllPairsEffectiveResistance(
  const size_t chunkSize) const
{
  Eigen::MatrixXd X = solveIdentityChunked(chunkSize);

  // Zero out rows and columns corresponding to pinned reference nodes
  for (const size_t r : referenceNodePerComponent_) {
    X.row(static_cast<Eigen::Index>(r)).setZero();
    X.col(static_cast<Eigen::Index>(r)).setZero();
  }

  const Eigen::VectorXd diagX = X.diagonal();
  const Eigen::Index N = static_cast<Eigen::Index>(numNodes_);
  Eigen::MatrixXd R(N, N);
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(static) if (N > 128)
#endif
  for (Eigen::Index j = 0; j < N; ++j) {
    const size_t compJ = nodeComponentId_[static_cast<size_t>(j)];
    const double diagJ = diagX(j);
    for (Eigen::Index i = 0; i < N; ++i) {
      if (i == j) {
        R(i, j) = 0.0;
      } else if (compJ != nodeComponentId_[static_cast<size_t>(i)]) {
        R(i, j) = std::numeric_limits<double>::infinity();
      } else {
        const double res = diagX(i) + diagJ - 2.0 * X(i, j);
        R(i, j) = std::max(0.0, res);
      }
    }
  }

  return R;
}

double
ResistanceDistanceAnalysis::computeKirchhoffIndex(const size_t chunkSize) const
{
  if (!isConnected()) {
    return std::numeric_limits<double>::infinity();
  }
  if (numNodes_ <= 1) {
    return 0.0;
  }

  return 0.5 * computeTotalResistancePerNode(chunkSize).sum();
}

double
ResistanceDistanceAnalysis::computeKirchhoffIndexSpectral() const
{
  if (!isConnected()) {
    return std::numeric_limits<double>::infinity();
  }
  if (numNodes_ <= 1) {
    return 0.0;
  }

  const Eigen::MatrixXd denseL = Eigen::MatrixXd(L_);
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(denseL);
  RUNTIME_EXP_IFN(es.info() == Eigen::Success,
                  "Spectral decomposition of L failed");

  const Eigen::VectorXd eigenvalues = es.eigenvalues();
  double sumReciprocal = 0.0;
  for (Eigen::Index i = 1; i < eigenvalues.size(); ++i) {
    RUNTIME_EXP_IFN(
      eigenvalues(i) > 0.0,
      "Connected graph Laplacian has non-positive non-trivial eigenvalue");
    sumReciprocal += 1.0 / eigenvalues(i);
  }

  return static_cast<double>(numNodes_) * sumReciprocal;
}

Eigen::MatrixXd
ResistanceDistanceAnalysis::computeMoorePenrosePseudoInverse(
  const size_t chunkSize) const
{
  Eigen::MatrixXd X = solveIdentityChunked(chunkSize);

  for (const size_t r : referenceNodePerComponent_) {
    X.row(static_cast<Eigen::Index>(r)).setZero();
    X.col(static_cast<Eigen::Index>(r)).setZero();
  }

  if (connectedComponents_.size() == 1) {
    if (numNodes_ <= 1) {
      return Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(numNodes_),
                                   static_cast<Eigen::Index>(numNodes_));
    }
    // P * X * P equals double-centering rows and columns in O(N^2) time
    const Eigen::VectorXd rowMeans = X.rowwise().mean();
    const Eigen::RowVectorXd colMeans = X.colwise().mean();
    const double grandMean = rowMeans.mean();
    X.colwise() -= rowMeans;
    X.rowwise() -= colMeans;
    X.array() += grandMean;
    return X;
  }

  Eigen::MatrixXd Lplus = Eigen::MatrixXd::Zero(
    static_cast<Eigen::Index>(numNodes_), static_cast<Eigen::Index>(numNodes_));

  for (const std::vector<size_t>& comp : connectedComponents_) {
    const size_t nComp = comp.size();
    if (nComp <= 1) {
      continue;
    }

    const Eigen::Index nCompIdx = static_cast<Eigen::Index>(nComp);
    Eigen::VectorXd rowMeans = Eigen::VectorXd::Zero(nCompIdx);
    Eigen::VectorXd colMeans = Eigen::VectorXd::Zero(nCompIdx);

    for (size_t b = 0; b < nComp; ++b) {
      const Eigen::Index colIdx = static_cast<Eigen::Index>(comp[b]);
      double colSum = 0.0;
      for (size_t a = 0; a < nComp; ++a) {
        const double val = X(static_cast<Eigen::Index>(comp[a]), colIdx);
        rowMeans(static_cast<Eigen::Index>(a)) += val;
        colSum += val;
      }
      colMeans(static_cast<Eigen::Index>(b)) =
        colSum / static_cast<double>(nComp);
    }
    rowMeans /= static_cast<double>(nComp);
    const double grandMean = rowMeans.mean();

    for (size_t b = 0; b < nComp; ++b) {
      const Eigen::Index colIdx = static_cast<Eigen::Index>(comp[b]);
      const double cMean = colMeans(static_cast<Eigen::Index>(b));
      for (size_t a = 0; a < nComp; ++a) {
        const Eigen::Index rowIdx = static_cast<Eigen::Index>(comp[a]);
        Lplus(rowIdx, colIdx) = X(rowIdx, colIdx) -
                                rowMeans(static_cast<Eigen::Index>(a)) - cMean +
                                grandMean;
      }
    }
  }

  return Lplus;
}

double
ResistanceDistanceAnalysis::computeCommuteTime(const size_t i,
                                               const size_t j) const
{
  const double R = computeEffectiveResistance(i, j);
  if (std::isinf(R)) {
    return std::numeric_limits<double>::infinity();
  }
  return 2.0 * componentTotalWeight_[nodeComponentId_[i]] * R;
}

Eigen::MatrixXd
ResistanceDistanceAnalysis::computeAllPairsCommuteTime(
  const size_t chunkSize) const
{
  Eigen::MatrixXd C = computeAllPairsEffectiveResistance(chunkSize);
  const Eigen::Index N = C.rows();
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(static) if (N > 128)
#endif
  for (Eigen::Index j = 0; j < N; ++j) {
    const double twoM =
      2.0 * componentTotalWeight_[nodeComponentId_[static_cast<size_t>(j)]];
    for (Eigen::Index i = 0; i < N; ++i) {
      const double r = C(i, j);
      C(i, j) =
        std::isinf(r) ? std::numeric_limits<double>::infinity() : twoM * r;
    }
  }
  return C;
}

Eigen::VectorXd
ResistanceDistanceAnalysis::computeTotalResistancePerNode(
  const size_t chunkSize) const
{
  INVALIDARG_EXP_IFN(chunkSize > 0, "chunkSize must be > 0");
  if (numNodes_ <= 1) {
    return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numNodes_));
  }
  if (!isConnected()) {
    return Eigen::VectorXd::Constant(static_cast<Eigen::Index>(numNodes_),
                                     std::numeric_limits<double>::infinity());
  }

  ensureFactorized();

  const Eigen::Index N = static_cast<Eigen::Index>(numNodes_);
  const Eigen::Index chunkSizeIdx = static_cast<Eigen::Index>(chunkSize);
  const Eigen::Index numChunks = (N + chunkSizeIdx - 1) / chunkSizeIdx;
  const Eigen::Index rIdx =
    static_cast<Eigen::Index>(referenceNodePerComponent_[0]);

  // Compute rowSumsX = X * 1 via a single linear solve L_pinned * y = 1_pinned
  Eigen::VectorXd onesPinned = Eigen::VectorXd::Ones(N);
  onesPinned(rIdx) = 0.0;
  Eigen::VectorXd rowSumsX = solver_.solve(onesPinned);
  rowSumsX(rIdx) = 0.0;

  Eigen::VectorXd diagX = Eigen::VectorXd::Zero(N);
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic) if (N > 256 && numChunks > 1)
#endif
  for (Eigen::Index chunkIdx = 0; chunkIdx < numChunks; ++chunkIdx) {
    const Eigen::Index colStart = chunkIdx * chunkSizeIdx;
    const Eigen::Index curChunk = std::min(chunkSizeIdx, N - colStart);
    Eigen::MatrixXd I_chunk = Eigen::MatrixXd::Zero(N, curChunk);
    for (Eigen::Index c = 0; c < curChunk; ++c) {
      I_chunk(colStart + c, c) = 1.0;
    }

    const Eigen::MatrixXd X_chunk = solver_.solve(I_chunk);
    for (Eigen::Index c = 0; c < curChunk; ++c) {
      const Eigen::Index globalCol = colStart + c;
      diagX(globalCol) = (globalCol == rIdx) ? 0.0 : X_chunk(globalCol, c);
    }
  }

  const double traceX = diagX.sum();
  Eigen::VectorXd tot = Eigen::VectorXd::Zero(N);
  for (Eigen::Index i = 0; i < N; ++i) {
    tot(i) = std::max(0.0,
                      static_cast<double>(numNodes_) * diagX(i) + traceX -
                        2.0 * rowSumsX(i));
  }
  return tot;
}

Eigen::VectorXd
ResistanceDistanceAnalysis::computeResistanceCentrality(
  const size_t chunkSize) const
{
  const Eigen::VectorXd tot = computeTotalResistancePerNode(chunkSize);
  const Eigen::Index N = tot.size();
  Eigen::VectorXd centrality = Eigen::VectorXd::Zero(N);
  const double nMinusOne = static_cast<double>(N > 1 ? N - 1 : 1);

  for (Eigen::Index i = 0; i < N; ++i) {
    const double s = tot(i);
    if (s > 0.0 && !std::isinf(s)) {
      centrality(i) = nMinusOne / s;
    }
  }
  return centrality;
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::createLinearChain(const size_t numNodes)
{
  std::vector<size_t> from;
  std::vector<size_t> to;
  pylimer_tools::utils::buildLinearChainEdges(numNodes, from, to);
  return ResistanceDistanceAnalysis(from, to, {}, numNodes);
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::createRing(const size_t numNodes)
{
  std::vector<size_t> from;
  std::vector<size_t> to;
  pylimer_tools::utils::buildRingEdges(numNodes, from, to);
  return ResistanceDistanceAnalysis(from, to, {}, numNodes);
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::createComplete(const size_t numNodes)
{
  std::vector<size_t> from;
  std::vector<size_t> to;
  pylimer_tools::utils::buildCompleteGraphEdges(numNodes, from, to);
  return ResistanceDistanceAnalysis(from, to, {}, numNodes);
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::createStar(const size_t numArms,
                                       const size_t armLength)
{
  std::vector<size_t> from;
  std::vector<size_t> to;
  size_t numNodes = 0;
  pylimer_tools::utils::buildStarEdges(numArms, armLength, from, to, numNodes);
  return ResistanceDistanceAnalysis(from, to, {}, numNodes);
}

} // namespace pylimer_tools::topo
