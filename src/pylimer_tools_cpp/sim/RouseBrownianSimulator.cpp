#include "RouseBrownianSimulator.h"

#include "../utils/utilityMacros.h"
#include <algorithm>
#include <cmath>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::sim {

RouseBrownianSimulator::RouseBrownianSimulator(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  double gamma,
  double k,
  double T,
  double dt,
  size_t numNodes)
  : gamma_(gamma)
  , k_(k)
  , T_(T)
  , dt_(dt)
{
  INVALIDARG_EXP_IFN(gamma_ > 0.0, "Friction gamma must be positive");
  INVALIDARG_EXP_IFN(k_ > 0.0, "Spring constant k must be positive");
  INVALIDARG_EXP_IFN(T_ >= 0.0, "Temperature T must be non-negative");
  INVALIDARG_EXP_IFN(dt_ > 0.0, "Time step dt must be positive");

  initializeFromEdges(springFrom, springTo, numNodes);
  positions_ = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(numNodes_), 3);
  rhsBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);
  cgSolutionBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);
  updateShiftedMatrixAndFactorization();
}

RouseBrownianSimulator::RouseBrownianSimulator(
  const Eigen::SparseMatrix<double>& L,
  double gamma,
  double k,
  double T,
  double dt)
  : numNodes_(L.rows())
  , gamma_(gamma)
  , k_(k)
  , T_(T)
  , dt_(dt)
  , L_(L)
{
  INVALIDARG_EXP_IFN(L.rows() == L.cols(), "Laplacian matrix L must be square");
  INVALIDARG_EXP_IFN(L.rows() > 0, "Laplacian matrix L cannot be empty");
  INVALIDARG_EXP_IFN(gamma_ > 0.0, "Friction gamma must be positive");
  INVALIDARG_EXP_IFN(k_ > 0.0, "Spring constant k must be positive");
  INVALIDARG_EXP_IFN(T_ >= 0.0, "Temperature T must be non-negative");
  INVALIDARG_EXP_IFN(dt_ > 0.0, "Time step dt must be positive");

  const Eigen::SparseMatrix<double> diff =
    L_ - Eigen::SparseMatrix<double>(L_.transpose());
  INVALIDARG_EXP_IFN(diff.norm() <= 1e-10 * std::max(1.0, L_.norm()),
                     "Laplacian matrix L must be symmetric");
  for (int col = 0; col < L_.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L_, col); it; ++it) {
      if (it.row() == it.col()) {
        INVALIDARG_EXP_IFN(
          it.value() >= 0.0,
          "Laplacian matrix L must have non-negative diagonal entries");
      }
    }
  }

  positions_ = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(numNodes_), 3);
  rhsBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);
  cgSolutionBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);
  updateShiftedMatrixAndFactorization();
}

void
RouseBrownianSimulator::initializeFromEdges(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  size_t numNodes)
{
  INVALIDARG_EXP_IFN(springFrom.size() == springTo.size(),
                     "springFrom and springTo must have the same size");
  INVALIDARG_EXP_IFN(springFrom.size() > 0 || numNodes > 0,
                     "Graph cannot be completely empty");

  if (numNodes == 0) {
    size_t maxIdx = 0;
    for (size_t i = 0; i < springFrom.size(); ++i) {
      maxIdx = std::max(maxIdx, std::max(springFrom[i], springTo[i]));
    }
    numNodes_ = static_cast<Eigen::Index>(maxIdx + 1);
  } else {
    numNodes_ = static_cast<Eigen::Index>(numNodes);
    for (size_t i = 0; i < springFrom.size(); ++i) {
      INVALIDARG_EXP_IFN(springFrom[i] < numNodes && springTo[i] < numNodes,
                         "Node indices in springFrom/springTo exceed numNodes");
    }
  }

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(springFrom.size() * 4);
  Eigen::VectorXd diag = Eigen::VectorXd::Zero(numNodes_);

  for (size_t e = 0; e < springFrom.size(); ++e) {
    const Eigen::Index u = static_cast<Eigen::Index>(springFrom[e]);
    const Eigen::Index v = static_cast<Eigen::Index>(springTo[e]);
    if (u == v) {
      continue;
    }
    triplets.emplace_back(static_cast<int>(u), static_cast<int>(v), -1.0);
    triplets.emplace_back(static_cast<int>(v), static_cast<int>(u), -1.0);
    diag(u) += 1.0;
    diag(v) += 1.0;
  }

  for (Eigen::Index i = 0; i < numNodes_; ++i) {
    if (diag(i) != 0.0) {
      triplets.emplace_back(static_cast<int>(i), static_cast<int>(i), diag(i));
    }
  }

  L_.resize(numNodes_, numNodes_);
  L_.setFromTriplets(triplets.begin(), triplets.end());
}

void
RouseBrownianSimulator::updateShiftedMatrixAndFactorization()
{
  // A = I + dt * (k / gamma) * L
  double alpha = dt_ * (k_ / gamma_);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<size_t>(L_.nonZeros() + numNodes_));

  for (int col = 0; col < L_.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L_, col); it; ++it) {
      triplets.emplace_back(it.row(), it.col(), alpha * it.value());
    }
  }

  for (Eigen::Index i = 0; i < numNodes_; ++i) {
    triplets.emplace_back(static_cast<int>(i), static_cast<int>(i), 1.0);
  }

  A_.resize(numNodes_, numNodes_);
  A_.setFromTriplets(triplets.begin(), triplets.end());

  solverFactorized_ = false;
}

void
RouseBrownianSimulator::ensureFactorized() const
{
  if (!solverFactorized_) {
    switch (solverType_) {
      case RouseSolverType::DIRECT_LDLT:
        ldltSolver_.compute(A_);
        RUNTIME_EXP_IFN(ldltSolver_.info() == Eigen::Success,
                        "SimplicialLDLT factorization of A failed");
        break;
      case RouseSolverType::DIRECT_LLT:
        lltSolver_.compute(A_);
        RUNTIME_EXP_IFN(lltSolver_.info() == Eigen::Success,
                        "SimplicialLLT factorization of A failed");
        break;
      case RouseSolverType::CONJUGATE_GRADIENT:
        cgSolver_.compute(A_);
        RUNTIME_EXP_IFN(cgSolver_.info() == Eigen::Success,
                        "ConjugateGradient setup of A failed");
        break;
      default:
        throw std::invalid_argument("Unknown RouseSolverType");
    }
    solverFactorized_ = true;
  }
}

void
RouseBrownianSimulator::solveShiftedRHS(const Eigen::MatrixXd& rhs,
                                        Eigen::MatrixXd& result)
{
  ensureFactorized();
  switch (solverType_) {
    case RouseSolverType::DIRECT_LDLT:
      result = ldltSolver_.solve(rhs);
      break;
    case RouseSolverType::DIRECT_LLT:
      result = lltSolver_.solve(rhs);
      break;
    case RouseSolverType::CONJUGATE_GRADIENT: {
      cgSolutionBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);
      for (int c = 0; c < 3; ++c) {
        cgSolutionBuffer_.col(c) =
          cgSolver_.solveWithGuess(rhs.col(c), positions_.col(c));
      }
      result.noalias() = cgSolutionBuffer_;
      break;
    }
    default:
      throw std::invalid_argument("Unknown RouseSolverType");
  }
}

void
RouseBrownianSimulator::setParameters(double gamma,
                                      double k,
                                      double T,
                                      double dt)
{
  INVALIDARG_EXP_IFN(gamma > 0.0, "Friction gamma must be positive");
  INVALIDARG_EXP_IFN(k > 0.0, "Spring constant k must be positive");
  INVALIDARG_EXP_IFN(T >= 0.0, "Temperature T must be non-negative");
  INVALIDARG_EXP_IFN(dt > 0.0, "Time step dt must be positive");

  bool needsUpdate = (gamma != gamma_ || k != k_ || dt != dt_);
  gamma_ = gamma;
  k_ = k;
  T_ = T;
  dt_ = dt;

  if (needsUpdate) {
    updateShiftedMatrixAndFactorization();
  }
}

void
RouseBrownianSimulator::setGamma(double gamma)
{
  setParameters(gamma, k_, T_, dt_);
}

void
RouseBrownianSimulator::setK(double k)
{
  setParameters(gamma_, k, T_, dt_);
}

void
RouseBrownianSimulator::setTemperature(double T)
{
  setParameters(gamma_, k_, T, dt_);
}

void
RouseBrownianSimulator::setDt(double dt)
{
  setParameters(gamma_, k_, T_, dt);
}

void
RouseBrownianSimulator::setKB(double kB)
{
  INVALIDARG_EXP_IFN(kB > 0.0, "Boltzmann constant kB must be positive");
  kB_ = kB;
}

void
RouseBrownianSimulator::setSolverType(RouseSolverType type)
{
  if (type != solverType_) {
    solverType_ = type;
    solverFactorized_ = false;
  }
}

void
RouseBrownianSimulator::setSeed(uint64_t seed)
{
  rng_.seed(seed);
}

void
RouseBrownianSimulator::setPositions(const Eigen::MatrixXd& positions)
{
  INVALIDARG_EXP_IFN(positions.rows() == numNodes_ && positions.cols() == 3,
                     "Positions matrix must have dimensions (numNodes, 3)");
  positions_ = positions;
}

void
RouseBrownianSimulator::setPositionsFlattened(const Eigen::VectorXd& r,
                                              bool interleaved)
{
  INVALIDARG_EXP_IFN(r.size() == numNodes_ * 3,
                     "Flattened positions vector must have size 3 * numNodes");

  if (interleaved) {
    for (Eigen::Index i = 0; i < numNodes_; ++i) {
      positions_(i, 0) = r(3 * i);
      positions_(i, 1) = r(3 * i + 1);
      positions_(i, 2) = r(3 * i + 2);
    }
  } else {
    positions_.col(0) = r.segment(0, numNodes_);
    positions_.col(1) = r.segment(numNodes_, numNodes_);
    positions_.col(2) = r.segment(2 * numNodes_, numNodes_);
  }
}

Eigen::VectorXd
RouseBrownianSimulator::getPositionsFlattened(bool interleaved) const
{
  Eigen::VectorXd r(numNodes_ * 3);
  if (interleaved) {
    for (Eigen::Index i = 0; i < numNodes_; ++i) {
      r(3 * i) = positions_(i, 0);
      r(3 * i + 1) = positions_(i, 1);
      r(3 * i + 2) = positions_(i, 2);
    }
  } else {
    r.segment(0, numNodes_) = positions_.col(0);
    r.segment(numNodes_, numNodes_) = positions_.col(1);
    r.segment(2 * numNodes_, numNodes_) = positions_.col(2);
  }
  return r;
}

void
RouseBrownianSimulator::initializeLinearPositions(double bondLength, int axis)
{
  INVALIDARG_EXP_IFN(axis >= 0 && axis < 3,
                     "axis must be 0 (x), 1 (y), or 2 (z)");
  positions_ = Eigen::MatrixXd::Zero(numNodes_, 3);
  for (Eigen::Index i = 0; i < numNodes_; ++i) {
    positions_(i, axis) = static_cast<double>(i) * bondLength;
  }
}

void
RouseBrownianSimulator::initializeRandomWalk(double bondLength, uint64_t seed)
{
  std::mt19937_64 gen(seed);
  std::normal_distribution<double> dist(0.0, 1.0);

  positions_ = Eigen::MatrixXd::Zero(numNodes_, 3);
  for (Eigen::Index i = 1; i < numNodes_; ++i) {
    Eigen::Vector3d step(dist(gen), dist(gen), dist(gen));
    step.normalize();
    positions_.row(i) = positions_.row(i - 1) + (bondLength * step).transpose();
  }
}

void
RouseBrownianSimulator::step(size_t numSteps)
{
  const double D = getDiffusionCoefficient();
  const double sigma = (T_ > 0.0) ? std::sqrt(2.0 * D * dt_) : 0.0;
  rhsBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);

  for (size_t s = 0; s < numSteps; ++s) {
    if (sigma > 0.0) {
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numNodes_); ++i) {
        rhsBuffer_(i, 0) = positions_(i, 0) + sigma * normalDist_(rng_);
        rhsBuffer_(i, 1) = positions_(i, 1) + sigma * normalDist_(rng_);
        rhsBuffer_(i, 2) = positions_(i, 2) + sigma * normalDist_(rng_);
      }
    } else {
      rhsBuffer_.noalias() = positions_;
    }
    solveShiftedRHS(rhsBuffer_, positions_);
  }
}

void
RouseBrownianSimulator::stepWithExternalForces(const Eigen::MatrixXd& forces,
                                               size_t numSteps)
{
  INVALIDARG_EXP_IFN(
    forces.rows() == numNodes_ && forces.cols() == 3,
    "External forces matrix must have dimensions (numNodes, 3)");

  const double D = getDiffusionCoefficient();
  const double sigma = (T_ > 0.0) ? std::sqrt(2.0 * D * dt_) : 0.0;
  const double forceCoeff = dt_ / gamma_;
  rhsBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);

  for (size_t s = 0; s < numSteps; ++s) {
    if (sigma > 0.0) {
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numNodes_); ++i) {
        rhsBuffer_(i, 0) = positions_(i, 0) + forceCoeff * forces(i, 0) +
                           sigma * normalDist_(rng_);
        rhsBuffer_(i, 1) = positions_(i, 1) + forceCoeff * forces(i, 1) +
                           sigma * normalDist_(rng_);
        rhsBuffer_(i, 2) = positions_(i, 2) + forceCoeff * forces(i, 2) +
                           sigma * normalDist_(rng_);
      }
    } else {
      rhsBuffer_.noalias() = positions_ + forceCoeff * forces;
    }
    solveShiftedRHS(rhsBuffer_, positions_);
  }
}

void
RouseBrownianSimulator::stepDeterministic(size_t numSteps)
{
  rhsBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);
  for (size_t s = 0; s < numSteps; ++s) {
    rhsBuffer_.noalias() = positions_;
    solveShiftedRHS(rhsBuffer_, positions_);
  }
}

void
RouseBrownianSimulator::stepDeterministicWithExternalForces(
  const Eigen::MatrixXd& forces,
  size_t numSteps)
{
  INVALIDARG_EXP_IFN(
    forces.rows() == numNodes_ && forces.cols() == 3,
    "External forces matrix must have dimensions (numNodes, 3)");

  const double forceCoeff = dt_ / gamma_;
  rhsBuffer_.resize(static_cast<Eigen::Index>(numNodes_), 3);
  for (size_t s = 0; s < numSteps; ++s) {
    rhsBuffer_.noalias() = positions_ + forceCoeff * forces;
    solveShiftedRHS(rhsBuffer_, positions_);
  }
}

Eigen::Vector3d
RouseBrownianSimulator::getCenterOfMass() const
{
  return positions_.colwise().mean().transpose();
}

Eigen::Vector3d
RouseBrownianSimulator::getEndToEndVector(size_t nodeA, size_t nodeB) const
{
  const Eigen::Index nodeAIdx = static_cast<Eigen::Index>(nodeA);
  const Eigen::Index nodeBIdx = (nodeB == static_cast<size_t>(-1))
                                  ? (numNodes_ - 1)
                                  : static_cast<Eigen::Index>(nodeB);
  INVALIDINDEX_EXP_IFN(nodeAIdx < numNodes_ && nodeBIdx < numNodes_,
                       "End-to-end node index out of range");

  return (positions_.row(nodeBIdx) - positions_.row(nodeAIdx)).transpose();
}

double
RouseBrownianSimulator::getEndToEndDistance(size_t nodeA, size_t nodeB) const
{
  return getEndToEndVector(nodeA, nodeB).norm();
}

double
RouseBrownianSimulator::getRadiusOfGyration() const
{
  if (numNodes_ <= 1) {
    return 0.0;
  }
  const Eigen::Vector3d cm = getCenterOfMass();
  double rgSq = 0.0;
  for (Eigen::Index i = 0; i < numNodes_; ++i) {
    rgSq += (positions_.row(i).transpose() - cm).squaredNorm();
  }
  return std::sqrt(rgSq / static_cast<double>(numNodes_));
}

double
RouseBrownianSimulator::getPotentialEnergy() const
{
  double u = 0.0;
  for (int c = 0; c < 3; ++c) {
    const Eigen::VectorXd Lr = L_ * positions_.col(c);
    u += positions_.col(c).dot(Lr);
  }
  return 0.5 * k_ * u;
}

RouseSimulationResult
RouseBrownianSimulator::run(size_t totalSteps,
                            size_t sampleInterval,
                            bool saveTrajectory,
                            size_t endToEndNodeA,
                            size_t endToEndNodeB)
{
  INVALIDARG_EXP_IFN(sampleInterval > 0, "sampleInterval must be > 0");
  if (endToEndNodeB == static_cast<size_t>(-1)) {
    endToEndNodeB = static_cast<size_t>(numNodes_ - 1);
  }

  const size_t numSamples = totalSteps / sampleInterval + 1;
  RouseSimulationResult result;
  result.times = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numSamples));
  result.msd = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numSamples));
  result.msdCenterOfMass =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numSamples));
  result.endToEndCorrelation =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numSamples));

  if (saveTrajectory) {
    result.trajectory.reserve(numSamples);
  }

  const Eigen::MatrixXd initialPos = positions_;
  const Eigen::Vector3d initialCM = getCenterOfMass();
  const Eigen::Vector3d initialR_ee =
    getEndToEndVector(endToEndNodeA, endToEndNodeB);
  const double initialR_eeSq = initialR_ee.squaredNorm();

  // Record t = 0
  result.times(0) = 0.0;
  result.msd(0) = 0.0;
  result.msdCenterOfMass(0) = 0.0;
  result.endToEndCorrelation(0) = (initialR_eeSq > 1e-12) ? 1.0 : 0.0;
  if (saveTrajectory) {
    result.trajectory.push_back(positions_);
  }

  size_t currentStep = 0;
  const Eigen::Index numSamplesIdx = static_cast<Eigen::Index>(numSamples);
  Eigen::Index sampleIdx = 1;
  while (currentStep < totalSteps) {
    const size_t stepsToTake =
      std::min(sampleInterval, totalSteps - currentStep);
    step(stepsToTake);
    currentStep += stepsToTake;

    if (currentStep % sampleInterval == 0 && sampleIdx < numSamplesIdx) {
      const double t = static_cast<double>(currentStep) * dt_;
      result.times(sampleIdx) = t;

      // Bead MSD: (1/N) * sum_i |r_i(t) - r_i(0)|^2
      double msdVal = 0.0;
      for (Eigen::Index i = 0; i < numNodes_; ++i) {
        msdVal += (positions_.row(i) - initialPos.row(i)).squaredNorm();
      }
      result.msd(sampleIdx) = msdVal / static_cast<double>(numNodes_);

      // Center of mass MSD: |r_cm(t) - r_cm(0)|^2
      const Eigen::Vector3d curCM = getCenterOfMass();
      result.msdCenterOfMass(sampleIdx) = (curCM - initialCM).squaredNorm();

      // End-to-end vector autocorrelation: <R_ee(t) . R_ee(0)> / |R_ee(0)|^2
      const Eigen::Vector3d curR_ee =
        getEndToEndVector(endToEndNodeA, endToEndNodeB);
      if (initialR_eeSq > 1e-12) {
        result.endToEndCorrelation(sampleIdx) =
          curR_ee.dot(initialR_ee) / initialR_eeSq;
      } else {
        result.endToEndCorrelation(sampleIdx) = curR_ee.dot(initialR_ee);
      }

      if (saveTrajectory) {
        result.trajectory.push_back(positions_);
      }
      ++sampleIdx;
    }
  }

  return result;
}

Eigen::VectorXd
RouseBrownianSimulator::computeMSD(
  const std::vector<Eigen::MatrixXd>& trajectory,
  bool timeAveraged)
{
  size_t M = trajectory.size();
  INVALIDARG_EXP_IFN(M > 0, "Trajectory cannot be empty");
  const Eigen::Index numFrames = static_cast<Eigen::Index>(M);
  const Eigen::Index N = trajectory[0].rows();

  Eigen::VectorXd msd = Eigen::VectorXd::Zero(numFrames);

  if (!timeAveraged) {
    const Eigen::MatrixXd& r0 = trajectory[0];
    for (size_t t = 0; t < M; ++t) {
      double sumSq = 0.0;
      const Eigen::MatrixXd& rt = trajectory[t];
      for (Eigen::Index i = 0; i < N; ++i) {
        sumSq += (rt.row(i) - r0.row(i)).squaredNorm();
      }
      msd(static_cast<Eigen::Index>(t)) = sumSq / static_cast<double>(N);
    }
  } else {
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic, 16) if (numFrames > 64)
#endif
    for (Eigen::Index tau = 0; tau < numFrames; ++tau) {
      double sumSq = 0.0;
      const size_t count = M - static_cast<size_t>(tau);
      for (size_t t = 0; t < count; ++t) {
        const Eigen::MatrixXd& rTau = trajectory[t + static_cast<size_t>(tau)];
        const Eigen::MatrixXd& rT = trajectory[t];
        for (Eigen::Index i = 0; i < N; ++i) {
          sumSq += (rTau.row(i) - rT.row(i)).squaredNorm();
        }
      }
      msd(tau) = sumSq / (static_cast<double>(count) * static_cast<double>(N));
    }
  }

  return msd;
}

Eigen::VectorXd
RouseBrownianSimulator::computeCenterOfMassMSD(
  const std::vector<Eigen::MatrixXd>& trajectory,
  bool timeAveraged)
{
  size_t M = trajectory.size();
  INVALIDARG_EXP_IFN(M > 0, "Trajectory cannot be empty");

  std::vector<Eigen::Vector3d> cms(M);
  for (size_t t = 0; t < M; ++t) {
    cms[t] = trajectory[t].colwise().mean().transpose();
  }

  Eigen::VectorXd msd = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(M));

  if (!timeAveraged) {
    for (size_t t = 0; t < M; ++t) {
      msd(static_cast<Eigen::Index>(t)) = (cms[t] - cms[0]).squaredNorm();
    }
  } else {
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic, 16) if (M > 64)
#endif
    for (size_t tau = 0; tau < M; ++tau) {
      double sumSq = 0.0;
      size_t count = M - tau;
      for (size_t t = 0; t < count; ++t) {
        sumSq += (cms[t + tau] - cms[t]).squaredNorm();
      }
      msd(static_cast<Eigen::Index>(tau)) = sumSq / static_cast<double>(count);
    }
  }

  return msd;
}

Eigen::VectorXd
RouseBrownianSimulator::computeEndToEndAutocorrelation(
  const std::vector<Eigen::MatrixXd>& trajectory,
  size_t nodeA,
  size_t nodeB,
  bool timeAveraged)
{
  size_t M = trajectory.size();
  INVALIDARG_EXP_IFN(M > 0, "Trajectory cannot be empty");
  size_t N = static_cast<size_t>(trajectory[0].rows());
  if (nodeB == static_cast<size_t>(-1)) {
    nodeB = N - 1;
  }
  INVALIDINDEX_EXP_IFN(nodeA < N && nodeB < N, "Node index out of range");

  const Eigen::Index nodeAIdx = static_cast<Eigen::Index>(nodeA);
  const Eigen::Index nodeBIdx = static_cast<Eigen::Index>(nodeB);
  std::vector<Eigen::Vector3d> r_ee(M);
  for (size_t t = 0; t < M; ++t) {
    r_ee[t] =
      (trajectory[t].row(nodeBIdx) - trajectory[t].row(nodeAIdx)).transpose();
  }

  Eigen::VectorXd corr = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(M));

  if (!timeAveraged) {
    double norm0Sq = r_ee[0].squaredNorm();
    double normFactor = (norm0Sq > 1e-12) ? norm0Sq : 1.0;
    for (size_t t = 0; t < M; ++t) {
      corr(static_cast<Eigen::Index>(t)) = r_ee[t].dot(r_ee[0]) / normFactor;
    }
  } else {
    double totalNormSq = 0.0;
    for (size_t t = 0; t < M; ++t) {
      totalNormSq += r_ee[t].squaredNorm();
    }
    double meanNormSq = totalNormSq / static_cast<double>(M);
    double normFactor = (meanNormSq > 1e-12) ? meanNormSq : 1.0;

#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic, 16) if (M > 64)
#endif
    for (size_t tau = 0; tau < M; ++tau) {
      double sumDot = 0.0;
      size_t count = M - tau;
      for (size_t t = 0; t < count; ++t) {
        sumDot += r_ee[t + tau].dot(r_ee[t]);
      }
      corr(static_cast<Eigen::Index>(tau)) =
        (sumDot / static_cast<double>(count)) / normFactor;
    }
  }

  return corr;
}

RouseBrownianSimulator
RouseBrownianSimulator::createLinearChain(size_t numNodes,
                                          double gamma,
                                          double k,
                                          double T,
                                          double dt)
{
  INVALIDARG_EXP_IFN(numNodes > 0, "numNodes must be > 0");
  if (numNodes == 1) {
    return RouseBrownianSimulator({}, {}, gamma, k, T, dt, 1);
  }
  std::vector<size_t> from(numNodes - 1);
  std::vector<size_t> to(numNodes - 1);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    from[i] = i;
    to[i] = i + 1;
  }
  return RouseBrownianSimulator(from, to, gamma, k, T, dt, numNodes);
}

RouseBrownianSimulator
RouseBrownianSimulator::createRing(size_t numNodes,
                                   double gamma,
                                   double k,
                                   double T,
                                   double dt)
{
  INVALIDARG_EXP_IFN(numNodes >= 3, "Ring requires at least 3 nodes");
  std::vector<size_t> from(numNodes);
  std::vector<size_t> to(numNodes);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    from[i] = i;
    to[i] = i + 1;
  }
  from[numNodes - 1] = numNodes - 1;
  to[numNodes - 1] = 0;
  return RouseBrownianSimulator(from, to, gamma, k, T, dt, numNodes);
}

RouseBrownianSimulator
RouseBrownianSimulator::createStar(size_t numArms,
                                   size_t armLength,
                                   double gamma,
                                   double k,
                                   double T,
                                   double dt)
{
  INVALIDARG_EXP_IFN(numArms >= 2, "Star requires at least 2 arms");
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
  return RouseBrownianSimulator(from, to, gamma, k, T, dt, numNodes);
}

} // namespace pylimer_tools::sim
