#include "ResistanceDistanceAnalysis.h"

#include "../utils/utilityMacros.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_set>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::topo {

ResistanceDistanceAnalysis::ResistanceDistanceAnalysis(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  const std::vector<double>& weights,
  size_t numNodes)
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
  , solver_()
  , solverFactorized_(false)
{
  INVALIDARG_EXP_IFN(L.rows() == L.cols(), "Laplacian matrix L must be square");
  INVALIDARG_EXP_IFN(L.rows() > 0, "Laplacian matrix L cannot be empty");

  // Compute totalWeight_ = 0.5 * Tr(L) and count off-diagonal edges
  double diagSum = 0.0;
  size_t offDiagCount = 0;
  for (int k = 0; k < L_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L_, k); it; ++it) {
      if (it.row() == it.col()) {
        diagSum += it.value();
      } else if (it.row() < it.col()) {
        ++offDiagCount;
      }
    }
  }
  totalWeight_ = 0.5 * diagSum;
  numEdges_ = offDiagCount;

  findComponentsAndSetupPinnedSystem();
}

void
ResistanceDistanceAnalysis::initializeFromEdges(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  const std::vector<double>& weights,
  size_t numNodes)
{
  INVALIDARG_EXP_IFN(springFrom.size() == springTo.size(),
                     "springFrom and springTo must have the same size");
  INVALIDARG_EXP_IFN(springFrom.size() > 0 || numNodes > 0,
                     "Graph cannot be completely empty");
  INVALIDARG_EXP_IFN(
    weights.empty() || weights.size() == springFrom.size(),
    "weights must be empty or have the same size as springFrom");

  numEdges_ = springFrom.size();

  if (numNodes == 0) {
    size_t maxIdx = 0;
    for (size_t i = 0; i < springFrom.size(); ++i) {
      maxIdx = std::max(maxIdx, std::max(springFrom[i], springTo[i]));
    }
    numNodes_ = maxIdx + 1;
  } else {
    numNodes_ = numNodes;
    for (size_t i = 0; i < springFrom.size(); ++i) {
      INVALIDARG_EXP_IFN(springFrom[i] < numNodes_ && springTo[i] < numNodes_,
                         "Node indices in springFrom/springTo exceed numNodes");
    }
  }

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(springFrom.size() * 4);
  Eigen::VectorXd diag = Eigen::VectorXd::Zero(numNodes_);
  totalWeight_ = 0.0;

  for (size_t e = 0; e < springFrom.size(); ++e) {
    size_t u = springFrom[e];
    size_t v = springTo[e];
    if (u == v) {
      continue; // Skip self loops in Laplacian assembly
    }
    double w = weights.empty() ? 1.0 : weights[e];
    INVALIDARG_EXP_IFN(w > 0.0, "Edge weight/spring constant must be positive");

    triplets.emplace_back(static_cast<int>(u), static_cast<int>(v), -w);
    triplets.emplace_back(static_cast<int>(v), static_cast<int>(u), -w);
    diag(static_cast<Eigen::Index>(u)) += w;
    diag(static_cast<Eigen::Index>(v)) += w;
    totalWeight_ += w;
  }

  for (size_t i = 0; i < numNodes_; ++i) {
    if (diag(static_cast<Eigen::Index>(i)) != 0.0) {
      triplets.emplace_back(static_cast<int>(i),
                            static_cast<int>(i),
                            diag(static_cast<Eigen::Index>(i)));
    }
  }

  L_.resize(numNodes_, numNodes_);
  L_.setFromTriplets(triplets.begin(), triplets.end());

  findComponentsAndSetupPinnedSystem();
}

void
ResistanceDistanceAnalysis::findComponentsAndSetupPinnedSystem()
{
  std::vector<std::vector<size_t>> adj(numNodes_);
  for (int k = 0; k < L_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L_, k); it; ++it) {
      int r = it.row();
      int c = it.col();
      if (r != c && std::abs(it.value()) > 1e-14) {
        adj[static_cast<size_t>(r)].push_back(static_cast<size_t>(c));
      }
    }
  }

  std::vector<bool> visited(numNodes_, false);
  connectedComponents_.clear();
  nodeComponentId_.assign(numNodes_, 0);
  referenceNodePerComponent_.clear();

  for (size_t i = 0; i < numNodes_; ++i) {
    if (!visited[i]) {
      std::vector<size_t> comp;
      std::queue<size_t> q;
      visited[i] = true;
      q.push(i);
      size_t compId = connectedComponents_.size();

      while (!q.empty()) {
        size_t u = q.front();
        q.pop();
        comp.push_back(u);
        nodeComponentId_[u] = compId;
        for (size_t v : adj[u]) {
          if (!visited[v]) {
            visited[v] = true;
            q.push(v);
          }
        }
      }

      referenceNodePerComponent_.push_back(comp[0]);
      connectedComponents_.push_back(std::move(comp));
    }
  }

  std::vector<bool> isRefNode(numNodes_, false);
  for (size_t r : referenceNodePerComponent_) {
    isRefNode[r] = true;
  }

  std::vector<Eigen::Triplet<double>> pinnedTriplets;
  pinnedTriplets.reserve(L_.nonZeros() + referenceNodePerComponent_.size());

  for (int k = 0; k < L_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L_, k); it; ++it) {
      size_t r = static_cast<size_t>(it.row());
      size_t c = static_cast<size_t>(it.col());
      if (isRefNode[r] || isRefNode[c]) {
        // Zero out off-diagonals on reference rows and columns
        if (r == c && isRefNode[r]) {
          // Set reference diagonal entry to 1.0 to pin potential v_r = 0
          pinnedTriplets.emplace_back(
            static_cast<int>(r), static_cast<int>(c), 1.0);
        }
      } else {
        pinnedTriplets.emplace_back(
          static_cast<int>(r), static_cast<int>(c), it.value());
      }
    }
  }

  // Ensure reference nodes always have diagonal 1.0 even if isolated
  for (size_t r : referenceNodePerComponent_) {
    pinnedTriplets.emplace_back(static_cast<int>(r), static_cast<int>(r), 1.0);
  }

  L_pinned_.resize(numNodes_, numNodes_);
  L_pinned_.setFromTriplets(pinnedTriplets.begin(), pinnedTriplets.end());

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
ResistanceDistanceAnalysis::computeEffectiveResistance(size_t i, size_t j) const
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
  size_t comp = nodeComponentId_[i];
  size_t r = referenceNodePerComponent_[comp];

  Eigen::VectorXd b = Eigen::VectorXd::Zero(numNodes_);
  if (i != r) {
    b(static_cast<Eigen::Index>(i)) = 1.0;
  }
  if (j != r) {
    b(static_cast<Eigen::Index>(j)) = -1.0;
  }

  Eigen::VectorXd v = solver_.solve(b);
  double vi = (i == r) ? 0.0 : v(static_cast<Eigen::Index>(i));
  double vj = (j == r) ? 0.0 : v(static_cast<Eigen::Index>(j));

  return std::max(0.0, vi - vj);
}

Eigen::VectorXd
ResistanceDistanceAnalysis::computeEffectiveResistanceBatched(
  const std::vector<std::pair<size_t, size_t>>& pairs) const
{
  size_t B = pairs.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(B));
  if (B == 0) {
    return result;
  }

  ensureFactorized();

  Eigen::MatrixXd Brhs = Eigen::MatrixXd::Zero(numNodes_, B);
  std::vector<bool> isInfinite(B, false);

  for (size_t b = 0; b < B; ++b) {
    size_t i = pairs[b].first;
    size_t j = pairs[b].second;
    INVALIDINDEX_EXP_IFN(i < numNodes_ && j < numNodes_,
                         "Node index out of range in batched resistance pairs");

    if (i == j) {
      continue;
    }
    if (nodeComponentId_[i] != nodeComponentId_[j]) {
      isInfinite[b] = true;
      continue;
    }

    size_t comp = nodeComponentId_[i];
    size_t r = referenceNodePerComponent_[comp];
    if (i != r) {
      Brhs(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(b)) += 1.0;
    }
    if (j != r) {
      Brhs(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(b)) -= 1.0;
    }
  }

  // Single Multi-RHS linear solve for all B query pairs simultaneously
  Eigen::MatrixXd V = solver_.solve(Brhs);

  for (size_t b = 0; b < B; ++b) {
    size_t i = pairs[b].first;
    size_t j = pairs[b].second;
    if (i == j) {
      result(static_cast<Eigen::Index>(b)) = 0.0;
    } else if (isInfinite[b]) {
      result(static_cast<Eigen::Index>(b)) =
        std::numeric_limits<double>::infinity();
    } else {
      size_t comp = nodeComponentId_[i];
      size_t r = referenceNodePerComponent_[comp];
      double vi = (i == r) ? 0.0
                           : V(static_cast<Eigen::Index>(i),
                               static_cast<Eigen::Index>(b));
      double vj = (j == r) ? 0.0
                           : V(static_cast<Eigen::Index>(j),
                               static_cast<Eigen::Index>(b));
      result(static_cast<Eigen::Index>(b)) = std::max(0.0, vi - vj);
    }
  }

  return result;
}

Eigen::MatrixXd
ResistanceDistanceAnalysis::solveIdentityChunked(size_t chunkSize) const
{
  ensureFactorized();
  if (numNodes_ == 0) {
    return Eigen::MatrixXd(0, 0);
  }

  Eigen::Index N = static_cast<Eigen::Index>(numNodes_);
  Eigen::MatrixXd X(N, N);
  if (chunkSize == 0) {
    chunkSize = 128;
  }

  // Level-3 BLAS chunked forward/back-substitutions:
  // Solving in chunks of B columns converts triangular solves from memory-bound
  // dtrsv (Level-2 BLAS) into compute-bound dtrsm (Level-3 BLAS), exploiting
  // AMX/Accelerate
  for (size_t colStart = 0; colStart < numNodes_; colStart += chunkSize) {
    size_t curChunk = std::min(chunkSize, numNodes_ - colStart);
    Eigen::MatrixXd I_chunk =
      Eigen::MatrixXd::Zero(N, static_cast<Eigen::Index>(curChunk));
    for (size_t c = 0; c < curChunk; ++c) {
      I_chunk(static_cast<Eigen::Index>(colStart + c),
              static_cast<Eigen::Index>(c)) = 1.0;
    }
    X.middleCols(static_cast<Eigen::Index>(colStart),
                 static_cast<Eigen::Index>(curChunk)) = solver_.solve(I_chunk);
  }

  return X;
}

Eigen::MatrixXd
ResistanceDistanceAnalysis::computeAllPairsEffectiveResistance(
  size_t chunkSize) const
{
  Eigen::MatrixXd X = solveIdentityChunked(chunkSize);

  // Zero out rows and columns corresponding to pinned reference nodes
  for (size_t r : referenceNodePerComponent_) {
    X.row(static_cast<Eigen::Index>(r)).setZero();
    X.col(static_cast<Eigen::Index>(r)).setZero();
  }

  Eigen::Index N = static_cast<Eigen::Index>(numNodes_);
  Eigen::MatrixXd R(N, N);
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(static) if (numNodes_ > 128)
#endif
  for (size_t i = 0; i < numNodes_; ++i) {
    for (size_t j = 0; j < numNodes_; ++j) {
      if (i == j) {
        R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = 0.0;
      } else if (nodeComponentId_[i] != nodeComponentId_[j]) {
        R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          std::numeric_limits<double>::infinity();
      } else {
        double res =
          X(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) +
          X(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(j)) -
          2.0 * X(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
        R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          std::max(0.0, res);
      }
    }
  }

  return R;
}

double
ResistanceDistanceAnalysis::computeKirchhoffIndex() const
{
  if (!isConnected()) {
    return std::numeric_limits<double>::infinity();
  }
  if (numNodes_ <= 1) {
    return 0.0;
  }

  Eigen::MatrixXd R = computeAllPairsEffectiveResistance();
  double K = 0.0;
#ifdef OPENMP_FOUND
#pragma omp parallel for reduction(+ : K) schedule(static) if (numNodes_ > 128)
#endif
  for (size_t i = 0; i < numNodes_; ++i) {
    for (size_t j = i + 1; j < numNodes_; ++j) {
      K += R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
    }
  }
  return K;
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

  Eigen::MatrixXd denseL = Eigen::MatrixXd(L_);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(denseL);
  RUNTIME_EXP_IFN(es.info() == Eigen::Success,
                  "Spectral decomposition of L failed");

  Eigen::VectorXd eigenvalues = es.eigenvalues();
  double sumReciprocal = 0.0;
  for (int i = 0; i < eigenvalues.size(); ++i) {
    if (eigenvalues(i) > 1e-10) {
      sumReciprocal += 1.0 / eigenvalues(i);
    }
  }

  return static_cast<double>(numNodes_) * sumReciprocal;
}

Eigen::MatrixXd
ResistanceDistanceAnalysis::computeMoorePenrosePseudoInverse(
  size_t chunkSize) const
{
  Eigen::MatrixXd X = solveIdentityChunked(chunkSize);

  for (size_t r : referenceNodePerComponent_) {
    X.row(static_cast<Eigen::Index>(r)).setZero();
    X.col(static_cast<Eigen::Index>(r)).setZero();
  }

  Eigen::MatrixXd Lplus = Eigen::MatrixXd::Zero(numNodes_, numNodes_);

  for (const auto& comp : connectedComponents_) {
    size_t nComp = comp.size();
    if (nComp <= 1) {
      continue;
    }

    Eigen::MatrixXd Xc(nComp, nComp);
    for (size_t a = 0; a < nComp; ++a) {
      for (size_t b = 0; b < nComp; ++b) {
        Xc(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) =
          X(static_cast<Eigen::Index>(comp[a]),
            static_cast<Eigen::Index>(comp[b]));
      }
    }

    Eigen::MatrixXd Pc = Eigen::MatrixXd::Identity(nComp, nComp) -
                         (1.0 / static_cast<double>(nComp)) *
                           Eigen::MatrixXd::Constant(nComp, nComp, 1.0);

    Eigen::MatrixXd Lplus_c = Pc * Xc * Pc;

    for (size_t a = 0; a < nComp; ++a) {
      for (size_t b = 0; b < nComp; ++b) {
        Lplus(static_cast<Eigen::Index>(comp[a]),
              static_cast<Eigen::Index>(comp[b])) =
          Lplus_c(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b));
      }
    }
  }

  return Lplus;
}

double
ResistanceDistanceAnalysis::computeCommuteTime(size_t i, size_t j) const
{
  double R = computeEffectiveResistance(i, j);
  if (std::isinf(R)) {
    return std::numeric_limits<double>::infinity();
  }
  return 2.0 * totalWeight_ * R;
}

Eigen::MatrixXd
ResistanceDistanceAnalysis::computeAllPairsCommuteTime() const
{
  Eigen::MatrixXd R = computeAllPairsEffectiveResistance();
  double twoM = 2.0 * totalWeight_;
  Eigen::MatrixXd C(numNodes_, numNodes_);
  for (size_t i = 0; i < numNodes_; ++i) {
    for (size_t j = 0; j < numNodes_; ++j) {
      double r = R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
      C(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
        std::isinf(r) ? std::numeric_limits<double>::infinity() : twoM * r;
    }
  }
  return C;
}

Eigen::VectorXd
ResistanceDistanceAnalysis::computeTotalResistancePerNode() const
{
  Eigen::MatrixXd R = computeAllPairsEffectiveResistance();
  Eigen::VectorXd tot = Eigen::VectorXd::Zero(numNodes_);
  for (size_t i = 0; i < numNodes_; ++i) {
    double s = 0.0;
    for (size_t j = 0; j < numNodes_; ++j) {
      if (i != j) {
        s += R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
      }
    }
    tot(static_cast<Eigen::Index>(i)) = s;
  }
  return tot;
}

Eigen::VectorXd
ResistanceDistanceAnalysis::computeResistanceCentrality() const
{
  Eigen::VectorXd tot = computeTotalResistancePerNode();
  Eigen::VectorXd centrality = Eigen::VectorXd::Zero(numNodes_);
  double nMinusOne = static_cast<double>(numNodes_ > 1 ? numNodes_ - 1 : 1);

  for (size_t i = 0; i < numNodes_; ++i) {
    double s = tot(static_cast<Eigen::Index>(i));
    if (s > 0.0 && !std::isinf(s)) {
      centrality(static_cast<Eigen::Index>(i)) = nMinusOne / s;
    }
  }
  return centrality;
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::createLinearChain(size_t numNodes)
{
  INVALIDARG_EXP_IFN(numNodes > 0, "numNodes must be > 0");
  if (numNodes == 1) {
    return ResistanceDistanceAnalysis({}, {}, {}, 1);
  }
  std::vector<size_t> from(numNodes - 1);
  std::vector<size_t> to(numNodes - 1);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    from[i] = i;
    to[i] = i + 1;
  }
  return ResistanceDistanceAnalysis(from, to, {}, numNodes);
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::createRing(size_t numNodes)
{
  INVALIDARG_EXP_IFN(numNodes >= 3, "Ring graph requires at least 3 nodes");
  std::vector<size_t> from(numNodes);
  std::vector<size_t> to(numNodes);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    from[i] = i;
    to[i] = i + 1;
  }
  from[numNodes - 1] = numNodes - 1;
  to[numNodes - 1] = 0;
  return ResistanceDistanceAnalysis(from, to, {}, numNodes);
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::createComplete(size_t numNodes)
{
  INVALIDARG_EXP_IFN(numNodes >= 2, "Complete graph requires at least 2 nodes");
  std::vector<size_t> from;
  std::vector<size_t> to;
  size_t nEdges = (numNodes * (numNodes - 1)) / 2;
  from.reserve(nEdges);
  to.reserve(nEdges);
  for (size_t i = 0; i < numNodes; ++i) {
    for (size_t j = i + 1; j < numNodes; ++j) {
      from.push_back(i);
      to.push_back(j);
    }
  }
  return ResistanceDistanceAnalysis(from, to, {}, numNodes);
}

ResistanceDistanceAnalysis
ResistanceDistanceAnalysis::createStar(size_t numArms, size_t armLength)
{
  INVALIDARG_EXP_IFN(numArms >= 2, "Star graph requires at least 2 arms");
  INVALIDARG_EXP_IFN(armLength >= 1, "armLength must be >= 1");
  size_t numNodes = 1 + numArms * armLength;
  std::vector<size_t> from;
  std::vector<size_t> to;
  from.reserve(numArms * armLength);
  to.reserve(numArms * armLength);

  size_t nextNode = 1;
  for (size_t a = 0; a < numArms; ++a) {
    from.push_back(0);
    to.push_back(nextNode);
    for (size_t l = 1; l < armLength; ++l) {
      from.push_back(nextNode + l - 1);
      to.push_back(nextNode + l);
    }
    nextNode += armLength;
  }
  return ResistanceDistanceAnalysis(from, to, {}, numNodes);
}

} // namespace pylimer_tools::topo
