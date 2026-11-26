#include "RouseBrownianSimulator.h"

#include "../entities/Universe.h"
#include "../utils/GraphLaplacianUtils.h"
#include "../utils/utilityMacros.h"
#include <algorithm>
#include <cmath>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::sim {

namespace {

  void validateTrajectoryFrames(const std::vector<Eigen::MatrixXd>& trajectory)
  {
    INVALIDARG_EXP_IFN(!trajectory.empty(), "Trajectory cannot be empty");
    const Eigen::Index numNodes = trajectory[0].rows();
    INVALIDARG_EXP_IFN(
      numNodes > 0 && trajectory[0].cols() == 3,
      "Trajectory frames must have dimensions (numNodes, 3) with numNodes > 0");
    for (size_t t = 1; t < trajectory.size(); ++t) {
      INVALIDARG_EXP_IFN(
        trajectory[t].rows() == numNodes && trajectory[t].cols() == 3,
        "All trajectory frames must have uniform (numNodes, 3) dimensions");
    }
  }

} // namespace

void
RouseBrownianSimulator::validateSimulationParameters(const double gamma,
                                                     const double k,
                                                     const double T,
                                                     const double dt)
{
  INVALIDARG_EXP_IFN(std::isfinite(gamma) && gamma > 0.0,
                     "Friction gamma must be finite and positive");
  INVALIDARG_EXP_IFN(std::isfinite(k) && k > 0.0,
                     "Spring constant k must be finite and positive");
  INVALIDARG_EXP_IFN(std::isfinite(T) && T >= 0.0,
                     "Temperature T must be finite and non-negative");
  INVALIDARG_EXP_IFN(std::isfinite(dt) && dt > 0.0,
                     "Time step dt must be finite and positive");
}

RouseBrownianSimulator::RouseBrownianSimulator(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  const double gamma,
  const double k,
  const double T,
  const double dt,
  const size_t numNodes,
  const std::vector<double>& weights)
  : gamma_(gamma)
  , k_(k)
  , T_(T)
  , dt_(dt)
{
  validateSimulationParameters(gamma_, k_, T_, dt_);

  initializeFromEdges(springFrom, springTo, numNodes, weights);
  positions_ = Eigen::MatrixXd::Zero(numNodes_, 3);
  rhsBuffer_.resize(numNodes_, 3);
  cgSolutionBuffer_.resize(numNodes_, 3);
  updateShiftedMatrixAndFactorization();
}

RouseBrownianSimulator::RouseBrownianSimulator(
  const Eigen::SparseMatrix<double>& L,
  const double gamma,
  const double k,
  const double T,
  const double dt)
  : numNodes_(L.rows())
  , gamma_(gamma)
  , k_(k)
  , T_(T)
  , dt_(dt)
  , L_(L)
{
  validateSimulationParameters(gamma_, k_, T_, dt_);
  L_.makeCompressed();
  pylimer_tools::utils::validateGraphLaplacian(L_);

  positions_ = Eigen::MatrixXd::Zero(numNodes_, 3);
  rhsBuffer_.resize(numNodes_, 3);
  cgSolutionBuffer_.resize(numNodes_, 3);
  updateShiftedMatrixAndFactorization();
}

RouseBrownianSimulator::RouseBrownianSimulator(
  const pylimer_tools::entities::Universe& universe,
  const double gamma,
  const double k,
  const double T,
  const double dt)
  : gamma_(gamma)
  , k_(k)
  , T_(T)
  , dt_(dt)
{
  validateSimulationParameters(gamma_, k_, T_, dt_);

  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  size_t nNodes = 0;
  pylimer_tools::utils::extractUniverseEdges(
    universe, springFrom, springTo, nNodes);
  initializeFromEdges(springFrom, springTo, nNodes, {});

  positions_.resize(numNodes_, 3);
  const Eigen::VectorXd coords =
    (universe.vertexPropertyExists("nx") &&
     universe.vertexPropertyExists("ny") && universe.vertexPropertyExists("nz"))
      ? universe.getUnwrappedVertexCoordinates(universe.getBox())
      : universe.getVertexCoordinates();
  setPositionsFlattened(coords, true);

  rhsBuffer_.resize(numNodes_, 3);
  cgSolutionBuffer_.resize(numNodes_, 3);
  updateShiftedMatrixAndFactorization();
}

RouseBrownianSimulator
RouseBrownianSimulator::fromMatrix(const Eigen::SparseMatrix<double>& L,
                                   const double gamma,
                                   const double k,
                                   const double T,
                                   const double dt)
{
  return RouseBrownianSimulator(L, gamma, k, T, dt);
}

RouseBrownianSimulator
RouseBrownianSimulator::fromUniverse(
  const pylimer_tools::entities::Universe& universe,
  const double gamma,
  const double k,
  const double T,
  const double dt)
{
  return RouseBrownianSimulator(universe, gamma, k, T, dt);
}

void
RouseBrownianSimulator::initializeFromEdges(
  const std::vector<size_t>& springFrom,
  const std::vector<size_t>& springTo,
  const size_t numNodes,
  const std::vector<double>& weights)
{
  size_t resolvedNumNodes = 0;
  L_ = pylimer_tools::utils::assembleGraphLaplacian(
    springFrom, springTo, weights, numNodes, &resolvedNumNodes);
  numNodes_ = static_cast<Eigen::Index>(resolvedNumNodes);
}

void
RouseBrownianSimulator::updateShiftedMatrixAndFactorization()
{
  // A = I + dt * (k / gamma) * L
  const double alpha = dt_ * (k_ / gamma_);
  Eigen::SparseMatrix<double> identity(numNodes_, numNodes_);
  identity.setIdentity();
  A_ = identity + alpha * L_;
  A_.makeCompressed();

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
        cgSolver_.setTolerance(1e-10);
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
      cgSolutionBuffer_.resize(numNodes_, 3);
      for (int c = 0; c < 3; ++c) {
        cgSolutionBuffer_.col(c) =
          cgSolver_.solveWithGuess(rhs.col(c), positions_.col(c));
        RUNTIME_EXP_IFN(cgSolver_.info() == Eigen::Success,
                        "ConjugateGradient solve of A failed");
      }
      result.noalias() = cgSolutionBuffer_;
      break;
    }
    default:
      throw std::invalid_argument("Unknown RouseSolverType");
  }
}

void
RouseBrownianSimulator::setParameters(const double gamma,
                                      const double k,
                                      const double T,
                                      const double dt)
{
  validateSimulationParameters(gamma, k, T, dt);

  const bool needsUpdate = (gamma != gamma_ || k != k_ || dt != dt_);
  gamma_ = gamma;
  k_ = k;
  T_ = T;
  dt_ = dt;

  if (needsUpdate) {
    updateShiftedMatrixAndFactorization();
  }
}

void
RouseBrownianSimulator::setGamma(const double gamma)
{
  setParameters(gamma, k_, T_, dt_);
}

void
RouseBrownianSimulator::setK(const double k)
{
  setParameters(gamma_, k, T_, dt_);
}

void
RouseBrownianSimulator::setTemperature(const double T)
{
  setParameters(gamma_, k_, T, dt_);
}

void
RouseBrownianSimulator::setDt(const double dt)
{
  setParameters(gamma_, k_, T_, dt);
}

void
RouseBrownianSimulator::setKB(const double kB)
{
  INVALIDARG_EXP_IFN(std::isfinite(kB) && kB > 0.0,
                     "Boltzmann constant kB must be finite and positive");
  kB_ = kB;
}

void
RouseBrownianSimulator::setSolverType(const RouseSolverType type)
{
  if (type != solverType_) {
    solverType_ = type;
    solverFactorized_ = false;
  }
}

void
RouseBrownianSimulator::setSeed(const uint64_t seed)
{
  rng_.seed(seed);
  normalDist_.reset();
}

void
RouseBrownianSimulator::setPositions(const Eigen::MatrixXd& positions)
{
  INVALIDARG_EXP_IFN(positions.rows() == numNodes_ && positions.cols() == 3,
                     "Positions matrix must have dimensions (numNodes, 3)");
  INVALIDARG_EXP_IFN(positions.allFinite(),
                     "Positions matrix must contain only finite values");
  positions_ = positions;
}

void
RouseBrownianSimulator::setPositionsFlattened(const Eigen::VectorXd& r,
                                              const bool interleaved)
{
  INVALIDARG_EXP_IFN(r.size() == numNodes_ * 3,
                     "Flattened positions vector must have size 3 * numNodes");
  INVALIDARG_EXP_IFN(
    r.allFinite(),
    "Flattened positions vector must contain only finite values");

  if (interleaved) {
    positions_.noalias() = r.reshaped<Eigen::RowMajor>(numNodes_, 3);
  } else {
    positions_.col(0) = r.segment(0, numNodes_);
    positions_.col(1) = r.segment(numNodes_, numNodes_);
    positions_.col(2) = r.segment(2 * numNodes_, numNodes_);
  }
}

Eigen::VectorXd
RouseBrownianSimulator::getPositionsFlattened(const bool interleaved) const
{
  Eigen::VectorXd r(numNodes_ * 3);
  if (interleaved) {
    r.reshaped(3, numNodes_).noalias() = positions_.transpose();
  } else {
    r.segment(0, numNodes_) = positions_.col(0);
    r.segment(numNodes_, numNodes_) = positions_.col(1);
    r.segment(2 * numNodes_, numNodes_) = positions_.col(2);
  }
  return r;
}

void
RouseBrownianSimulator::initializeLinearPositions(const double bondLength,
                                                  const int axis)
{
  INVALIDARG_EXP_IFN(std::isfinite(bondLength) && bondLength >= 0.0,
                     "bondLength must be finite and non-negative");
  INVALIDARG_EXP_IFN(axis >= 0 && axis < 3,
                     "axis must be 0 (x), 1 (y), or 2 (z)");
  positions_ = Eigen::MatrixXd::Zero(numNodes_, 3);
  for (Eigen::Index i = 0; i < numNodes_; ++i) {
    positions_(i, axis) = static_cast<double>(i) * bondLength;
  }
}

void
RouseBrownianSimulator::initializeRandomWalk(const double bondLength,
                                             const uint64_t seed)
{
  INVALIDARG_EXP_IFN(std::isfinite(bondLength) && bondLength >= 0.0,
                     "bondLength must be finite and non-negative");
  std::mt19937_64 gen(seed);
  std::normal_distribution<double> dist(0.0, 1.0);

  positions_ = Eigen::MatrixXd::Zero(numNodes_, 3);
  for (Eigen::Index i = 1; i < numNodes_; ++i) {
    const double dx = dist(gen);
    const double dy = dist(gen);
    const double dz = dist(gen);
    Eigen::Vector3d step(dx, dy, dz);
    const double norm = step.norm();
    if (norm > 0.0) {
      step /= norm;
    } else {
      step = Eigen::Vector3d(1.0, 0.0, 0.0);
    }
    positions_.row(i) = positions_.row(i - 1) + (bondLength * step).transpose();
  }
}

void
RouseBrownianSimulator::step(const size_t numSteps)
{
  const double D = getDiffusionCoefficient();
  const double sigma = (T_ > 0.0) ? std::sqrt(2.0 * D * dt_) : 0.0;
  rhsBuffer_.resize(numNodes_, 3);

  for (size_t s = 0; s < numSteps; ++s) {
    if (sigma > 0.0) {
      for (Eigen::Index i = 0; i < numNodes_; ++i) {
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
                                               const size_t numSteps)
{
  INVALIDARG_EXP_IFN(
    forces.rows() == numNodes_ && forces.cols() == 3,
    "External forces matrix must have dimensions (numNodes, 3)");
  INVALIDARG_EXP_IFN(forces.allFinite(),
                     "External forces matrix must contain only finite values");

  const double D = getDiffusionCoefficient();
  const double sigma = (T_ > 0.0) ? std::sqrt(2.0 * D * dt_) : 0.0;
  const double forceCoeff = dt_ / gamma_;
  rhsBuffer_.resize(numNodes_, 3);

  for (size_t s = 0; s < numSteps; ++s) {
    if (sigma > 0.0) {
      for (Eigen::Index i = 0; i < numNodes_; ++i) {
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
RouseBrownianSimulator::stepDeterministic(const size_t numSteps)
{
  rhsBuffer_.resize(numNodes_, 3);
  for (size_t s = 0; s < numSteps; ++s) {
    rhsBuffer_.noalias() = positions_;
    solveShiftedRHS(rhsBuffer_, positions_);
  }
}

void
RouseBrownianSimulator::stepDeterministicWithExternalForces(
  const Eigen::MatrixXd& forces,
  const size_t numSteps)
{
  INVALIDARG_EXP_IFN(
    forces.rows() == numNodes_ && forces.cols() == 3,
    "External forces matrix must have dimensions (numNodes, 3)");
  INVALIDARG_EXP_IFN(forces.allFinite(),
                     "External forces matrix must contain only finite values");

  const double forceCoeff = dt_ / gamma_;
  rhsBuffer_.resize(numNodes_, 3);
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
RouseBrownianSimulator::getEndToEndVector(
  const size_t nodeA,
  const std::optional<size_t> nodeB) const
{
  const size_t numNodesSz = static_cast<size_t>(numNodes_);
  const size_t resolvedNodeB = nodeB.value_or(numNodesSz - 1);
  INVALIDINDEX_EXP_IFN(nodeA < numNodesSz && resolvedNodeB < numNodesSz,
                       "End-to-end node index out of range");

  return (positions_.row(static_cast<Eigen::Index>(resolvedNodeB)) -
          positions_.row(static_cast<Eigen::Index>(nodeA)))
    .transpose();
}

double
RouseBrownianSimulator::getEndToEndDistance(
  const size_t nodeA,
  const std::optional<size_t> nodeB) const
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
  const double rgSq = (positions_.rowwise() - cm.transpose()).squaredNorm();
  return std::sqrt(rgSq / static_cast<double>(numNodes_));
}

double
RouseBrownianSimulator::getPotentialEnergy() const
{
  return 0.5 * k_ * positions_.cwiseProduct(L_ * positions_).sum();
}

RouseSimulationResult
RouseBrownianSimulator::run(const size_t totalSteps,
                            const size_t sampleInterval,
                            const bool saveTrajectory,
                            const size_t endToEndNodeA,
                            const std::optional<size_t> endToEndNodeB)
{
  INVALIDARG_EXP_IFN(sampleInterval > 0, "sampleInterval must be > 0");
  const size_t resolvedNodeB =
    endToEndNodeB.value_or(static_cast<size_t>(numNodes_ - 1));

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
    getEndToEndVector(endToEndNodeA, resolvedNodeB);
  const double initialR_eeSq = initialR_ee.squaredNorm();
  const double normFactor = (initialR_eeSq > 1e-12) ? initialR_eeSq : 1.0;

  // Record t = 0
  result.times(0) = 0.0;
  result.msd(0) = 0.0;
  result.msdCenterOfMass(0) = 0.0;
  result.endToEndCorrelation(0) = initialR_eeSq / normFactor;
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
      result.msd(sampleIdx) = (positions_ - initialPos).squaredNorm() /
                              static_cast<double>(numNodes_);

      // Center of mass MSD: |r_cm(t) - r_cm(0)|^2
      const Eigen::Vector3d curCM = getCenterOfMass();
      result.msdCenterOfMass(sampleIdx) = (curCM - initialCM).squaredNorm();

      // End-to-end vector autocorrelation: <R_ee(t) . R_ee(0)> / |R_ee(0)|^2
      const Eigen::Vector3d curR_ee =
        getEndToEndVector(endToEndNodeA, resolvedNodeB);
      result.endToEndCorrelation(sampleIdx) =
        curR_ee.dot(initialR_ee) / normFactor;

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
  const bool timeAveraged)
{
  validateTrajectoryFrames(trajectory);
  const size_t M = trajectory.size();
  const Eigen::Index numFrames = static_cast<Eigen::Index>(M);
  const Eigen::Index N = trajectory[0].rows();

  Eigen::VectorXd msd = Eigen::VectorXd::Zero(numFrames);

  if (!timeAveraged) {
    const Eigen::MatrixXd& r0 = trajectory[0];
    for (size_t t = 0; t < M; ++t) {
      msd(static_cast<Eigen::Index>(t)) =
        (trajectory[t] - r0).squaredNorm() / static_cast<double>(N);
    }
  } else {
#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic, 16) if (numFrames > 64)
#endif
    for (Eigen::Index tau = 0; tau < numFrames; ++tau) {
      double sumSq = 0.0;
      const size_t count = M - static_cast<size_t>(tau);
      for (size_t t = 0; t < count; ++t) {
        sumSq += (trajectory[t + static_cast<size_t>(tau)] - trajectory[t])
                   .squaredNorm();
      }
      msd(tau) = sumSq / (static_cast<double>(count) * static_cast<double>(N));
    }
  }

  return msd;
}

Eigen::VectorXd
RouseBrownianSimulator::computeCenterOfMassMSD(
  const std::vector<Eigen::MatrixXd>& trajectory,
  const bool timeAveraged)
{
  validateTrajectoryFrames(trajectory);
  const size_t M = trajectory.size();

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
      const size_t count = M - tau;
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
  const size_t nodeA,
  const std::optional<size_t> nodeB,
  const bool timeAveraged)
{
  validateTrajectoryFrames(trajectory);
  const size_t M = trajectory.size();
  const size_t N = static_cast<size_t>(trajectory[0].rows());
  const size_t resolvedNodeB = nodeB.value_or(N - 1);
  INVALIDINDEX_EXP_IFN(nodeA < N && resolvedNodeB < N,
                       "Node index out of range");

  const Eigen::Index nodeAIdx = static_cast<Eigen::Index>(nodeA);
  const Eigen::Index nodeBIdx = static_cast<Eigen::Index>(resolvedNodeB);
  std::vector<Eigen::Vector3d> r_ee(M);
  for (size_t t = 0; t < M; ++t) {
    r_ee[t] =
      (trajectory[t].row(nodeBIdx) - trajectory[t].row(nodeAIdx)).transpose();
  }

  Eigen::VectorXd corr = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(M));

  if (!timeAveraged) {
    const double norm0Sq = r_ee[0].squaredNorm();
    const double normFactor = (norm0Sq > 1e-12) ? norm0Sq : 1.0;
    for (size_t t = 0; t < M; ++t) {
      corr(static_cast<Eigen::Index>(t)) = r_ee[t].dot(r_ee[0]) / normFactor;
    }
  } else {
    double totalNormSq = 0.0;
    for (size_t t = 0; t < M; ++t) {
      totalNormSq += r_ee[t].squaredNorm();
    }
    const double meanNormSq = totalNormSq / static_cast<double>(M);
    const double normFactor = (meanNormSq > 1e-12) ? meanNormSq : 1.0;

#ifdef OPENMP_FOUND
#pragma omp parallel for schedule(dynamic, 16) if (M > 64)
#endif
    for (size_t tau = 0; tau < M; ++tau) {
      double sumDot = 0.0;
      const size_t count = M - tau;
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
RouseBrownianSimulator::createLinearChain(const size_t numNodes,
                                          const double gamma,
                                          const double k,
                                          const double T,
                                          const double dt)
{
  std::vector<size_t> from;
  std::vector<size_t> to;
  pylimer_tools::utils::buildLinearChainEdges(numNodes, from, to);
  return RouseBrownianSimulator(from, to, gamma, k, T, dt, numNodes);
}

RouseBrownianSimulator
RouseBrownianSimulator::createRing(const size_t numNodes,
                                   const double gamma,
                                   const double k,
                                   const double T,
                                   const double dt)
{
  std::vector<size_t> from;
  std::vector<size_t> to;
  pylimer_tools::utils::buildRingEdges(numNodes, from, to);
  return RouseBrownianSimulator(from, to, gamma, k, T, dt, numNodes);
}

RouseBrownianSimulator
RouseBrownianSimulator::createComplete(const size_t numNodes,
                                       const double gamma,
                                       const double k,
                                       const double T,
                                       const double dt)
{
  std::vector<size_t> from;
  std::vector<size_t> to;
  pylimer_tools::utils::buildCompleteGraphEdges(numNodes, from, to);
  return RouseBrownianSimulator(from, to, gamma, k, T, dt, numNodes);
}

RouseBrownianSimulator
RouseBrownianSimulator::createStar(const size_t numArms,
                                   const size_t armLength,
                                   const double gamma,
                                   const double k,
                                   const double T,
                                   const double dt)
{
  std::vector<size_t> from;
  std::vector<size_t> to;
  size_t numNodes = 0;
  pylimer_tools::utils::buildStarEdges(numArms, armLength, from, to, numNodes);
  return RouseBrownianSimulator(from, to, gamma, k, T, dt, numNodes);
}

} // namespace pylimer_tools::sim
