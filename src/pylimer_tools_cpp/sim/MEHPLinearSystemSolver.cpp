#ifdef __APPLE__
// Avoid header symbol conflict between Eigen's lapacke.h and Apple Accelerate.h
#undef EIGEN_USE_LAPACKE
#endif
#include "MEHPLinearSystemSolver.h"
#ifdef USE_METAL
#include "MetalForceBalanceSolver.h"
#endif
#include "../utils/ExtraEigenSolvers.h"
#include "../utils/ExtraEigenTypes.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#ifdef __APPLE__
#include <Eigen/AccelerateSupport>
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace pylimer_tools::sim::mehp {

namespace {

  [[nodiscard]] SLESolver resolveEffectiveSLESolver(
    const SLESolver solverChoice)
  {
    if (solverChoice == SLESolver::DEFAULT || solverChoice == SLESolver::AUTO) {
      return SLESolver::CONJUGATE_GRADIENT;
    }
    return solverChoice;
  }

  void zeroZDisplacementsIf2D(Eigen::VectorXd& u,
                              const size_t nrOfLinks,
                              const bool is2D)
  {
    if (!is2D || nrOfLinks == 0) {
      return;
    }
    u.reshaped(3, static_cast<Eigen::Index>(nrOfLinks)).row(2).setZero();
  }

  [[nodiscard]] double computeResidualNorm(
    const MEHPNetwork& net,
    const Eigen::VectorXd& u,
    const MEHPForceEvaluator& forceEvaluator,
    const bool is2D,
    Eigen::VectorXd& gradScratch,
    double* const outEnergy = nullptr)
  {
    if (net.nrOfLinks == 0 || net.nrOfSprings == 0) {
      if (outEnergy != nullptr) {
        *outEnergy = 0.0;
      }
      return 0.0;
    }
    gradScratch.resize(3 * static_cast<Eigen::Index>(net.nrOfLinks));
    const Eigen::VectorXd springVectors = net.evaluateSpringVectors(u, is2D);
    const double energy = forceEvaluator.evaluateNetworkEnergyAndGradient(
      net, springVectors, is2D, gradScratch.data());
    if (outEnergy != nullptr) {
      *outEnergy = energy;
    }
    return gradScratch.norm();
  }

  void assembleLinearSystem(const MEHPNetwork& net,
                            const double kappa,
                            const bool is2D,
                            std::vector<Eigen::Triplet<double>>& triplets,
                            Eigen::VectorXd& diagonal,
                            Eigen::SparseMatrix<double>& sysMatrix,
                            Eigen::VectorXd& constants)
  {
    const Eigen::Index nLinks = static_cast<Eigen::Index>(net.nrOfLinks);
    triplets.clear();
    constants.resize(3 * nLinks);
    constants.setZero();
    diagonal.resize(nLinks);
    diagonal.setZero();

    for (size_t springIdx = 0; springIdx < net.nrOfSprings; ++springIdx) {
      const Eigen::Index sIdx = static_cast<Eigen::Index>(springIdx);
      const int idxAInt = net.springIndexA[sIdx];
      const int idxBInt = net.springIndexB[sIdx];
      if (idxAInt < 0 || idxBInt < 0 || idxAInt == idxBInt) {
        continue;
      }
      const Eigen::Index idxA = static_cast<Eigen::Index>(idxAInt);
      const Eigen::Index idxB = static_cast<Eigen::Index>(idxBInt);
      const double multiplier = kappa / net.springContourLength[sIdx];

      triplets.emplace_back(idxAInt, idxBInt, multiplier);
      triplets.emplace_back(idxBInt, idxAInt, multiplier);
      diagonal[idxA] -= multiplier;
      diagonal[idxB] -= multiplier;

      const Eigen::Vector3d offset = net.springBoxOffset.segment<3>(3 * sIdx);
      constants.segment<3>(3 * idxA) -= offset * multiplier;
      constants.segment<3>(3 * idxB) += offset * multiplier;
    }

    if (is2D && nLinks > 0) {
      constants.reshaped(3, nLinks).row(2).setZero();
    }

    for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
      const int linkIdxInt = static_cast<int>(linkIdx);
      triplets.emplace_back(
        linkIdxInt, linkIdxInt, diagonal[static_cast<Eigen::Index>(linkIdx)]);
    }

    sysMatrix.resize(nLinks, nLinks);
    if (triplets.empty()) {
      sysMatrix.setZero();
      return;
    }
    sysMatrix.setFromTriplets(triplets.begin(), triplets.end());
    sysMatrix.makeCompressed();
  }

  void applyInPlaceMatrixUpdates(
    Eigen::SparseMatrix<double>& sysMatrix,
    const std::vector<BrokenSpringInfo>& brokenUpdates)
  {
    for (const BrokenSpringInfo& bs : brokenUpdates) {
      const Eigen::Index idxA = static_cast<Eigen::Index>(bs.idxA);
      const Eigen::Index idxB = static_cast<Eigen::Index>(bs.idxB);
      sysMatrix.coeffRef(idxA, idxB) -= bs.multiplier;
      sysMatrix.coeffRef(idxB, idxA) -= bs.multiplier;
      sysMatrix.coeffRef(idxA, idxA) += bs.multiplier;
      sysMatrix.coeffRef(idxB, idxB) += bs.multiplier;
    }
  }

  [[nodiscard]] Eigen::MatrixX3d unpackVectorToMatrixX3(
    const Eigen::VectorXd& vec,
    const Eigen::Index nLinks)
  {
    if (nLinks == 0) {
      return Eigen::MatrixX3d(0, 3);
    }
    return vec.reshaped(3, nLinks).transpose();
  }

  [[nodiscard]] Eigen::VectorXd packMatrixX3ToCoordinates(
    Eigen::MatrixX3d& xMat,
    const Eigen::VectorXd& origCoords,
    const Eigen::Index nLinks,
    const bool is2D)
  {
    Eigen::VectorXd flat(3 * nLinks);
    if (nLinks == 0) {
      return flat;
    }
    flat.reshaped(3, nLinks) = xMat.transpose();
    if (is2D) {
      flat.reshaped(3, nLinks).row(2) = origCoords.reshaped(3, nLinks).row(2);
    }
    return flat;
  }

  [[noreturn]] void throwSolverFailure(const char* const action,
                                       const SLESolver effectiveSolverChoice,
                                       const Eigen::ComputationInfo info)
  {
    RUNTIME_EXP(std::string(action) + " with solver " +
                std::to_string(effectiveSolverChoice) + " (" +
                SLESolverNames[effectiveSolverChoice] +
                "). Solver info: " + std::to_string(info));
  }

  template<typename SolverType>
  [[nodiscard]] Eigen::VectorXd solveDirectScalar(
    SolverType& solver,
    const Eigen::SparseMatrix<double>& sysMatrix,
    const Eigen::VectorXd& constants,
    const Eigen::VectorXd& origCoords,
    const Eigen::Index nLinks,
    const bool is2D,
    const SLESolver effectiveSolverChoice,
    ExitReason& exitReason)
  {
    if (constants.isZero(1e-9)) {
      exitReason = ExitReason::X_TOLERANCE;
      return Eigen::VectorXd::Zero(constants.size());
    }

    solver.compute(sysMatrix);
    if (solver.info() != Eigen::Success) {
      throwSolverFailure("System matrix computation failed",
                         effectiveSolverChoice,
                         solver.info());
    }

    const Eigen::MatrixX3d rhs = unpackVectorToMatrixX3(constants, nLinks);
    Eigen::MatrixX3d sol(nLinks, 3);
    if (is2D) {
      sol.leftCols<2>().noalias() = solver.solve(rhs.leftCols<2>());
      if (solver.info() != Eigen::Success) {
        throwSolverFailure(
          "System could not be solved", effectiveSolverChoice, solver.info());
      }
    } else {
      sol = solver.solve(rhs);
      if (solver.info() != Eigen::Success) {
        throwSolverFailure(
          "System could not be solved", effectiveSolverChoice, solver.info());
      }
    }

    exitReason = ExitReason::X_TOLERANCE;
    return packMatrixX3ToCoordinates(sol, origCoords, nLinks, is2D);
  }

  template<typename SolverType>
  [[nodiscard]] Eigen::VectorXd solveIterativeScalar(
    SolverType& solver,
    const Eigen::SparseMatrix<double>& sysMatrix,
    const Eigen::VectorXd& constants,
    const Eigen::VectorXd& origCoords,
    const Eigen::Index nLinks,
    const bool is2D,
    const double residualReduction,
    const int maxIterations,
    const SLESolver effectiveSolverChoice,
    size_t& iterationsDone,
    ExitReason& exitReason)
  {
    if (constants.isZero(1e-9)) {
      exitReason = ExitReason::X_TOLERANCE;
      return Eigen::VectorXd::Zero(constants.size());
    }

    solver.setTolerance(residualReduction);
    solver.setMaxIterations(maxIterations);
    solver.compute(sysMatrix);
    if (solver.info() != Eigen::Success) {
      throwSolverFailure("System matrix computation failed",
                         effectiveSolverChoice,
                         solver.info());
    }

    const Eigen::MatrixX3d rhs = unpackVectorToMatrixX3(constants, nLinks);
    Eigen::MatrixX3d sol(nLinks, 3);
    const int numCols = is2D ? 2 : 3;
    bool anyNoConvergence = false;
    for (int col = 0; col < numCols; ++col) {
      sol.col(col) = solver.solve(rhs.col(col));
      if (solver.info() == Eigen::NoConvergence) {
        anyNoConvergence = true;
      } else if (solver.info() != Eigen::Success) {
        throwSolverFailure(
          "System could not be solved", effectiveSolverChoice, solver.info());
      }
      iterationsDone += static_cast<size_t>(solver.iterations());
    }

    exitReason =
      anyNoConvergence ? ExitReason::MAX_STEPS : ExitReason::X_TOLERANCE;
    return packMatrixX3ToCoordinates(sol, origCoords, nLinks, is2D);
  }

  template<typename SolveColumnFn>
  [[nodiscard]] Eigen::VectorXd solveGradientDescentScalar(
    SolveColumnFn&& solveColumn,
    const Eigen::VectorXd& constants,
    const Eigen::VectorXd& initialSolution,
    const Eigen::VectorXd& origCoords,
    const Eigen::Index nLinks,
    const bool is2D,
    const int maxIterations,
    size_t& iterationsDone,
    ExitReason& exitReason)
  {
    if (constants.isZero(1e-9)) {
      exitReason = ExitReason::X_TOLERANCE;
      return Eigen::VectorXd::Zero(constants.size());
    }

    const Eigen::MatrixX3d rhs = unpackVectorToMatrixX3(constants, nLinks);
    Eigen::MatrixX3d sol(nLinks, 3);
    const int numCols = is2D ? 2 : 3;
    int maxColIterations = 0;
    for (int col = 0; col < numCols; ++col) {
      const Eigen::VectorXd initCol =
        initialSolution.reshaped(3, nLinks).row(col).transpose();
      int colIterations = 0;
      sol.col(col) = solveColumn(rhs.col(col), initCol, colIterations);
      maxColIterations = std::max(maxColIterations, colIterations);
      iterationsDone += static_cast<size_t>(colIterations);
    }

    exitReason = (maxColIterations >= maxIterations) ? ExitReason::MAX_STEPS
                                                     : ExitReason::X_TOLERANCE;
    return packMatrixX3ToCoordinates(sol, origCoords, nLinks, is2D);
  }

  [[nodiscard]] Eigen::VectorXd solveConjugateGradientWithFallback(
    const Eigen::SparseMatrix<double>& sysMatrix,
    const Eigen::VectorXd& constants,
    const Eigen::VectorXd& initialSolution,
    const Eigen::VectorXd& origCoords,
    const Eigen::Index nLinks,
    const bool is2D,
    const double residualReduction,
    const int maxIterations,
    const int cgFallbackMaxIterations,
    SLESolver& effectiveSolverChoice,
    bool& cgFallbackTriggered,
    size_t& iterationsDone,
    ExitReason& exitReason)
  {
    if (constants.isZero(1e-9)) {
      exitReason = ExitReason::X_TOLERANCE;
      return Eigen::VectorXd::Zero(constants.size());
    }

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                             Eigen::Lower | Eigen::Upper,
                             Eigen::DiagonalPreconditioner<double>>
      cgSolver;
    const double targetTol = residualReduction;
    cgSolver.setTolerance(targetTol);
    const int cgMaxIters = (cgFallbackMaxIterations > 0)
                             ? std::min(maxIterations, cgFallbackMaxIterations)
                             : maxIterations;
    cgSolver.setMaxIterations(cgMaxIters);

    bool fallbackNeeded = false;
    cgSolver.compute(sysMatrix);
    if (cgSolver.info() != Eigen::Success) {
      if (cgFallbackMaxIterations > 0) {
        fallbackNeeded = true;
      } else {
        throwSolverFailure("System matrix computation failed",
                           effectiveSolverChoice,
                           cgSolver.info());
      }
    }

    const Eigen::MatrixX3d rhs = unpackVectorToMatrixX3(constants, nLinks);
    Eigen::MatrixX3d sol = unpackVectorToMatrixX3(initialSolution, nLinks);
    const int numCols = is2D ? 2 : 3;
    bool anyNoConvergence = false;

    if (!fallbackNeeded) {
      for (int col = 0; col < numCols; ++col) {
        const bool useGuess =
          sol.col(col).allFinite() && !sol.col(col).isZero();
        if (useGuess) {
          sol.col(col) = cgSolver.solveWithGuess(rhs.col(col), sol.col(col));
        } else {
          sol.col(col) = cgSolver.solve(rhs.col(col));
        }
        if (cgFallbackMaxIterations > 0 &&
            (cgSolver.info() != Eigen::Success ||
             (cgSolver.iterations() >= cgMaxIters &&
              cgSolver.error() > targetTol))) {
          iterationsDone += static_cast<size_t>(cgSolver.iterations());
          fallbackNeeded = true;
          break;
        }
        if (cgSolver.info() == Eigen::NoConvergence) {
          anyNoConvergence = true;
        } else if (cgSolver.info() != Eigen::Success) {
          throwSolverFailure("System could not be solved",
                             effectiveSolverChoice,
                             cgSolver.info());
        }
        iterationsDone += static_cast<size_t>(cgSolver.iterations());
      }
    }

    if (fallbackNeeded) {
      cgFallbackTriggered = true;
      effectiveSolverChoice = SLESolver::BICGSTAB;

      Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                      Eigen::DiagonalPreconditioner<double>>
        bicgstabSolver;
      bicgstabSolver.setTolerance(targetTol);
      bicgstabSolver.setMaxIterations(maxIterations);
      bicgstabSolver.compute(sysMatrix);
      if (bicgstabSolver.info() != Eigen::Success) {
        throwSolverFailure(
          "System matrix computation failed for BiCGSTAB fallback",
          effectiveSolverChoice,
          bicgstabSolver.info());
      }

      anyNoConvergence = false;
      for (int col = 0; col < numCols; ++col) {
        const bool useGuess =
          sol.col(col).allFinite() && !sol.col(col).isZero();
        if (useGuess) {
          const Eigen::VectorXd xGuess =
            bicgstabSolver.solveWithGuess(rhs.col(col), sol.col(col));
          if (bicgstabSolver.info() == Eigen::Success) {
            sol.col(col) = xGuess;
          } else {
            sol.col(col) = bicgstabSolver.solve(rhs.col(col));
          }
        } else {
          sol.col(col) = bicgstabSolver.solve(rhs.col(col));
        }

        if (bicgstabSolver.info() == Eigen::NoConvergence) {
          const double rhsNorm = rhs.col(col).norm();
          const double colResidualNorm =
            (sysMatrix * sol.col(col) - rhs.col(col)).norm();
          const bool convergedByResidual =
            colResidualNorm <= targetTol * std::max(rhsNorm, 1.0);
          if (!convergedByResidual) {
            anyNoConvergence = true;
          }
        } else if (bicgstabSolver.info() != Eigen::Success) {
          throwSolverFailure("System could not be solved by BiCGSTAB fallback",
                             effectiveSolverChoice,
                             bicgstabSolver.info());
        }
        iterationsDone += static_cast<size_t>(bicgstabSolver.iterations());
      }
    }

    exitReason =
      anyNoConvergence ? ExitReason::MAX_STEPS : ExitReason::X_TOLERANCE;
    return packMatrixX3ToCoordinates(sol, origCoords, nLinks, is2D);
  }

#ifdef USE_METAL
  [[nodiscard]] SpringRelaxationData buildMetalSpringRelaxationData(
    const MEHPNetwork& net,
    const MEHPLinearSolverConfig& cfg,
    const bool is2D)
  {
    SpringRelaxationData springData;
    springData.springIndexA = net.springIndexA.data();
    springData.springIndexB = net.springIndexB.data();
    springData.springBoxOffset = net.springBoxOffset.data();
    springData.springContourLength = net.springContourLength.data();
    springData.numSprings = net.nrOfSprings;
    springData.breakingThreshold = cfg.springBreakingLength;
    springData.inactiveCutoff =
      (cfg.simplificationMode == StructureSimplificationMode::INACTIVE_ONLY ||
       cfg.simplificationMode == StructureSimplificationMode::INACTIVE_THEN_X2F)
        ? cfg.inactiveRemovalCutoff
        : 0.0;
    springData.is2D = is2D;
    return springData;
  }
#endif

} // namespace

LinearSystemMEHPSolver::LinearSystemMEHPSolver(
  const MEHPLinearSolverConfig& solverConfig)
  : config(solverConfig)
{
}

const MEHPLinearSolverConfig&
LinearSystemMEHPSolver::getConfig() const
{
  return this->config;
}

void
LinearSystemMEHPSolver::setConfig(const MEHPLinearSolverConfig& newConfig)
{
  this->config = newConfig;
}

MEHPSolveResult
LinearSystemMEHPSolver::solve(MEHPNetwork& net,
                              Eigen::VectorXd& u,
                              const MEHPForceEvaluator& forceEvaluator,
                              const bool is2D,
                              const MEHPSolveCallbacks& callbacks) const
{
  INVALIDARG_EXP_IFN(
    forceEvaluator.isLinear(),
    "LinearSystemMEHPSolver requires a linear (Hookean) force evaluator.");
  INVALIDARG_EXP_IFN(this->config.slipMode == SlipMode::FIXED,
                     "LinearSystemMEHPSolver only supports SlipMode::FIXED.");
  INVALIDARG_EXP_IFN(
    u.size() == static_cast<Eigen::Index>(3 * net.nrOfLinks),
    "Displacement vector u size must match 3 * net.nrOfLinks.");
  INVALIDARG_EXP_IFN(
    this->config.simplificationMode !=
      StructureSimplificationMode::X1F_X2F_THEN_INACTIVE,
    "StructureSimplificationMode::X1F_X2F_THEN_INACTIVE is not implemented; "
    "use INACTIVE_THEN_X2F instead.");
  INVALIDARG_EXP_IFN(
    (this->config.inactiveRemovalCutoff > 0.0 &&
     std::isfinite(this->config.inactiveRemovalCutoff)) ||
      (this->config.simplificationMode !=
         StructureSimplificationMode::INACTIVE_ONLY &&
       this->config.simplificationMode !=
         StructureSimplificationMode::INACTIVE_THEN_X2F),
    "Removal cut-off must be finite and positive when inactive link "
    "simplification is enabled.");
  INVALIDARG_EXP_IFN(std::isfinite(this->config.residualReduction) &&
                       this->config.residualReduction >= 0.0,
                     "residualReduction must be finite and non-negative.");
  INVALIDARG_EXP_IFN(this->config.maxIterations > 0,
                     "maxIterations must be positive.");
  INVALIDARG_EXP_IFN(
    std::isfinite(this->config.gradientDescentLearningRate) &&
      this->config.gradientDescentLearningRate > 0.0,
    "gradientDescentLearningRate must be finite and positive.");

  const double kappa = forceEvaluator.getKappa();
  INVALIDARG_EXP_IFN(
    std::isfinite(kappa) && kappa >= 0.0,
    "Force evaluator spring constant kappa must be finite and non-negative.");

  net.validate(u);

  MEHPSolveResult result;
  SLESolver effectiveSolverChoice =
    resolveEffectiveSLESolver(this->config.solver);
  result.effectiveSLESolver = effectiveSolverChoice;

  if (net.nrOfLinks == 0 || net.nrOfSprings == 0 || net.nrOfStrands == 0) {
    result.exitReason = ExitReason::X_TOLERANCE;
    result.energy = 0.0;
    result.residual = 0.0;
    return result;
  }

  if (this->config.simplificationMode !=
      StructureSimplificationMode::NO_SIMPLIFICATION) {
    net.assumeComplete = false;
  }

  Eigen::VectorXd gradScratch;
  if (kappa == 0.0) {
    zeroZDisplacementsIf2D(u, net.nrOfLinks, is2D);
    result.exitReason = ExitReason::X_TOLERANCE;
    result.residual = computeResidualNorm(
      net, u, forceEvaluator, is2D, gradScratch, &result.energy);
    return result;
  }

#ifdef USE_METAL
  const double initialResidual =
    computeResidualNorm(net, u, forceEvaluator, is2D, gradScratch);
#endif

  Eigen::VectorXd initialSolution(static_cast<Eigen::Index>(3 * net.nrOfLinks));
  if (u.isZero()) {
    initialSolution.setZero();
  } else {
    initialSolution = net.coordinates + u;
  }

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(net.nrOfLinks + net.nrOfSprings * 2);
  Eigen::VectorXd constants(static_cast<Eigen::Index>(net.nrOfLinks * 3));
  Eigen::VectorXd diagonal(static_cast<Eigen::Index>(net.nrOfLinks));
  Eigen::SparseMatrix<double> sysMatrix;

  bool matrixNeedsRebuild = true;
  size_t totalSpringsAtLastRebuild = net.nrOfSprings;
  size_t cumulativeBrokenSinceRebuild = 0;
  bool wasInterrupted = false;
  size_t iterationsDone = 0;

#ifdef USE_METAL
  const bool isMetalSolver =
    (effectiveSolverChoice == SLESolver::METAL_CG ||
     effectiveSolverChoice == SLESolver::METAL_BARZILAI_BORWEIN);
#endif

  size_t nRemovedInIteration = 0;
  do {
    nRemovedInIteration = 0;

    if (matrixNeedsRebuild) {
      assembleLinearSystem(
        net, kappa, is2D, triplets, diagonal, sysMatrix, constants);
      if (triplets.empty()) {
        break;
      }
      matrixNeedsRebuild = false;
      totalSpringsAtLastRebuild = net.nrOfSprings;
      cumulativeBrokenSinceRebuild = 0;
    }

    if (sysMatrix.nonZeros() == 0) {
      break;
    }

#ifndef NDEBUG
    const bool isSelfAdjoint = Eigen::isSelfAdjoint(sysMatrix);
    assert(isSelfAdjoint);
#endif

    const Eigen::Index nLinks = static_cast<Eigen::Index>(net.nrOfLinks);
    Eigen::VectorXd finalCoordinates;
    bool interruptRequested = false;

#ifdef USE_METAL
    bool hasMetalSpringEval = false;
    uint32_t metalTotalFlagged = 0;
#endif

    switch (effectiveSolverChoice) {
      case SLESolver::CONJUGATE_GRADIENT:
      case SLESolver::CONJUGATE_GRADIENT_DIAGONALIZED: {
        finalCoordinates = solveConjugateGradientWithFallback(
          sysMatrix,
          constants,
          initialSolution,
          net.coordinates,
          nLinks,
          is2D,
          this->config.residualReduction,
          this->config.maxIterations,
          this->config.cgFallbackMaxIterations,
          effectiveSolverChoice,
          result.cgFallbackTriggered,
          iterationsDone,
          result.exitReason);
        break;
      }
      case SLESolver::CONJUGATE_GRADIENT_IDENTITY: {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower | Eigen::Upper,
                                 Eigen::IdentityPreconditioner>
          solver;
        finalCoordinates = solveIterativeScalar(solver,
                                                sysMatrix,
                                                constants,
                                                net.coordinates,
                                                nLinks,
                                                is2D,
                                                this->config.residualReduction,
                                                this->config.maxIterations,
                                                effectiveSolverChoice,
                                                iterationsDone,
                                                result.exitReason);
        break;
      }
      case SLESolver::CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY: {
        Eigen::ConjugateGradient<
          Eigen::SparseMatrix<double>,
          Eigen::Lower | Eigen::Upper,
          Eigen::IncompleteCholesky<double, Eigen::Lower | Eigen::Upper>>
          solver;
        finalCoordinates = solveIterativeScalar(solver,
                                                sysMatrix,
                                                constants,
                                                net.coordinates,
                                                nLinks,
                                                is2D,
                                                this->config.residualReduction,
                                                this->config.maxIterations,
                                                effectiveSolverChoice,
                                                iterationsDone,
                                                result.exitReason);
        break;
      }
      case SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT:
      case SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED: {
        Eigen::LeastSquaresConjugateGradient<
          Eigen::SparseMatrix<double>,
          Eigen::LeastSquareDiagonalPreconditioner<double>>
          solver;
        finalCoordinates = solveIterativeScalar(solver,
                                                sysMatrix,
                                                constants,
                                                net.coordinates,
                                                nLinks,
                                                is2D,
                                                this->config.residualReduction,
                                                this->config.maxIterations,
                                                effectiveSolverChoice,
                                                iterationsDone,
                                                result.exitReason);
        break;
      }
      case SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY: {
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>,
                                             Eigen::IdentityPreconditioner>
          solver;
        finalCoordinates = solveIterativeScalar(solver,
                                                sysMatrix,
                                                constants,
                                                net.coordinates,
                                                nLinks,
                                                is2D,
                                                this->config.residualReduction,
                                                this->config.maxIterations,
                                                effectiveSolverChoice,
                                                iterationsDone,
                                                result.exitReason);
        break;
      }
      case SLESolver::BICGSTAB:
      case SLESolver::BICGSTAB_DIAGONALIZED: {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                        Eigen::DiagonalPreconditioner<double>>
          solver;
        finalCoordinates = solveIterativeScalar(solver,
                                                sysMatrix,
                                                constants,
                                                net.coordinates,
                                                nLinks,
                                                is2D,
                                                this->config.residualReduction,
                                                this->config.maxIterations,
                                                effectiveSolverChoice,
                                                iterationsDone,
                                                result.exitReason);
        break;
      }
      case SLESolver::BICGSTAB_IDENTITY: {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                        Eigen::IdentityPreconditioner>
          solver;
        finalCoordinates = solveIterativeScalar(solver,
                                                sysMatrix,
                                                constants,
                                                net.coordinates,
                                                nLinks,
                                                is2D,
                                                this->config.residualReduction,
                                                this->config.maxIterations,
                                                effectiveSolverChoice,
                                                iterationsDone,
                                                result.exitReason);
        break;
      }
      case SLESolver::BICGSTAB_INCOMPLETE_LU: {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                        Eigen::IncompleteLUT<double>>
          solver;
        finalCoordinates = solveIterativeScalar(solver,
                                                sysMatrix,
                                                constants,
                                                net.coordinates,
                                                nLinks,
                                                is2D,
                                                this->config.residualReduction,
                                                this->config.maxIterations,
                                                effectiveSolverChoice,
                                                iterationsDone,
                                                result.exitReason);
        break;
      }
      case SLESolver::GRADIENT_DESCENT: {
        finalCoordinates = solveGradientDescentScalar(
          [&](const Eigen::VectorXd& bCol,
              const Eigen::VectorXd& initCol,
              int& colIterations) {
            return Eigen::gradientDescent(
              sysMatrix,
              bCol,
              this->config.gradientDescentLearningRate,
              this->config.residualReduction,
              this->config.maxIterations,
              colIterations,
              initCol,
              -1.0,
              [&callbacks, &interruptRequested](int, const Eigen::VectorXd&) {
                if (callbacks.shouldInterrupt && callbacks.shouldInterrupt()) {
                  interruptRequested = true;
                  return true;
                }
                return false;
              });
          },
          constants,
          initialSolution,
          net.coordinates,
          nLinks,
          is2D,
          this->config.maxIterations,
          iterationsDone,
          result.exitReason);
        break;
      }
      case SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT:
      case SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG: {
        const bool shortAlpha =
          (effectiveSolverChoice ==
           SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT);
        finalCoordinates = solveGradientDescentScalar(
          [&](const Eigen::VectorXd& bCol,
              const Eigen::VectorXd& initCol,
              int& colIterations) {
            return Eigen::gradientDescentBarzilaiBorwein(
              sysMatrix,
              bCol,
              this->config.gradientDescentLearningRate,
              this->config.residualReduction,
              this->config.maxIterations,
              colIterations,
              shortAlpha,
              initCol,
              -1.0,
              [&callbacks, &interruptRequested](int, const Eigen::VectorXd&) {
                if (callbacks.shouldInterrupt && callbacks.shouldInterrupt()) {
                  interruptRequested = true;
                  return true;
                }
                return false;
              });
          },
          constants,
          initialSolution,
          net.coordinates,
          nLinks,
          is2D,
          this->config.maxIterations,
          iterationsDone,
          result.exitReason);
        break;
      }
      case SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM: {
        finalCoordinates = solveGradientDescentScalar(
          [&](const Eigen::VectorXd& bCol,
              const Eigen::VectorXd& initCol,
              int& colIterations) {
            return Eigen::gradientDescentHeavyBallBarzilaiBorwein(
              sysMatrix,
              bCol,
              this->config.gradientDescentLearningRate,
              this->config.residualReduction,
              this->config.maxIterations,
              colIterations,
              initCol,
              -1.0,
              [&callbacks, &interruptRequested](int, const Eigen::VectorXd&) {
                if (callbacks.shouldInterrupt && callbacks.shouldInterrupt()) {
                  interruptRequested = true;
                  return true;
                }
                return false;
              });
          },
          constants,
          initialSolution,
          net.coordinates,
          nLinks,
          is2D,
          this->config.maxIterations,
          iterationsDone,
          result.exitReason);
        break;
      }
      case SLESolver::SIMPLICIAL_LLT: {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        finalCoordinates = solveDirectScalar(solver,
                                             sysMatrix,
                                             constants,
                                             net.coordinates,
                                             nLinks,
                                             is2D,
                                             effectiveSolverChoice,
                                             result.exitReason);
        break;
      }
      case SLESolver::SIMPLICIAL_LDLT: {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        finalCoordinates = solveDirectScalar(solver,
                                             sysMatrix,
                                             constants,
                                             net.coordinates,
                                             nLinks,
                                             is2D,
                                             effectiveSolverChoice,
                                             result.exitReason);
        break;
      }
      case SLESolver::SPARSE_LU: {
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        finalCoordinates = solveDirectScalar(solver,
                                             sysMatrix,
                                             constants,
                                             net.coordinates,
                                             nLinks,
                                             is2D,
                                             effectiveSolverChoice,
                                             result.exitReason);
        break;
      }
      case SLESolver::SPARSE_QR: {
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
          solver;
        finalCoordinates = solveDirectScalar(solver,
                                             sysMatrix,
                                             constants,
                                             net.coordinates,
                                             nLinks,
                                             is2D,
                                             effectiveSolverChoice,
                                             result.exitReason);
        break;
      }
      case SLESolver::ACCELERATE_LLT: {
#ifdef __APPLE__
        Eigen::AccelerateLLT<Eigen::SparseMatrix<double>> solver;
        finalCoordinates = solveDirectScalar(solver,
                                             sysMatrix,
                                             constants,
                                             net.coordinates,
                                             nLinks,
                                             is2D,
                                             effectiveSolverChoice,
                                             result.exitReason);
        break;
#else
        throw std::invalid_argument(
          "AccelerateLLT solver is only available on Apple platforms.");
#endif
      }
      case SLESolver::ACCELERATE_LDLT: {
#ifdef __APPLE__
        Eigen::AccelerateLDLT<Eigen::SparseMatrix<double>> solver;
        finalCoordinates = solveDirectScalar(solver,
                                             sysMatrix,
                                             constants,
                                             net.coordinates,
                                             nLinks,
                                             is2D,
                                             effectiveSolverChoice,
                                             result.exitReason);
        break;
#else
        throw std::invalid_argument(
          "AccelerateLDLT solver is only available on Apple platforms.");
#endif
      }
      case SLESolver::METAL_CG: {
#ifdef USE_METAL
        if (constants.isZero(1e-9)) {
          finalCoordinates = Eigen::VectorXd::Zero(constants.size());
          result.exitReason = ExitReason::X_TOLERANCE;
          break;
        }
        const SpringRelaxationData springData =
          buildMetalSpringRelaxationData(net, this->config, is2D);
        const GpuRelaxationStepResult stepResult =
          runGpuRelaxationStep(sysMatrix,
                               constants,
                               this->config.residualReduction,
                               this->config.maxIterations,
                               springData,
                               initialSolution,
                               /*projectZeroModes=*/true,
                               /*useJacobi=*/true,
                               /*pinnedNode=*/-1,
                               /*enableCpuFallback=*/true);

        iterationsDone += static_cast<size_t>(stepResult.solverIterations);
        result.exitReason =
          (stepResult.solverIterations >= this->config.maxIterations)
            ? ExitReason::MAX_STEPS
            : ExitReason::X_TOLERANCE;
        hasMetalSpringEval = true;
        metalTotalFlagged = stepResult.totalFlagged;
        finalCoordinates = getGpuResidentCoordinates(3 * net.nrOfLinks);
        break;
#else
        throw std::runtime_error(
          "Metal acceleration is not enabled in this build.");
#endif
      }
      case SLESolver::METAL_BARZILAI_BORWEIN: {
#ifdef USE_METAL
        if (constants.isZero(1e-9)) {
          finalCoordinates = Eigen::VectorXd::Zero(constants.size());
          result.exitReason = ExitReason::X_TOLERANCE;
          break;
        }
        int solverIterations = 0;
        finalCoordinates =
          solveSLEMetalBarzilaiBorwein(sysMatrix,
                                       constants,
                                       this->config.gradientDescentLearningRate,
                                       this->config.residualReduction,
                                       this->config.maxIterations,
                                       solverIterations,
                                       true,
                                       initialSolution,
                                       /*projectZeroModes=*/true);
        iterationsDone += static_cast<size_t>(solverIterations);
        result.exitReason = (solverIterations >= this->config.maxIterations)
                              ? ExitReason::MAX_STEPS
                              : ExitReason::X_TOLERANCE;
        break;
#else
        throw std::runtime_error(
          "Metal acceleration is not enabled in this build.");
#endif
      }
      default:
        throw std::invalid_argument("Unsupported SLESolver enumerator.");
    }

    u = finalCoordinates - net.coordinates;
    zeroZDisplacementsIf2D(u, net.nrOfLinks, is2D);

    iterationsDone += 1;

    if (interruptRequested) {
      wasInterrupted = true;
      if (callbacks.onStep) {
        callbacks.onStep(iterationsDone);
      }
      break;
    }

#ifdef USE_METAL
    const double currentResidual =
      computeResidualNorm(net, u, forceEvaluator, is2D, gradScratch);
    if (isMetalSolver) {
      if (!hasMetalSpringEval && !constants.isZero(1e-9)) {
        const SpringRelaxationData springData =
          buildMetalSpringRelaxationData(net, this->config, is2D);
        const SpringEvaluationResult evalResult =
          evaluateSpringStrainAndInactivityMetal(net.nrOfLinks, springData);
        hasMetalSpringEval = true;
        metalTotalFlagged = evalResult.totalFlagged;
      }

      if (hasMetalSpringEval && metalTotalFlagged == 0 &&
          this->config.simplificationMode ==
            StructureSimplificationMode::NO_SIMPLIFICATION &&
          currentResidual <= initialResidual * this->config.residualReduction) {
        if (callbacks.onStep) {
          callbacks.onStep(iterationsDone);
        }
        if (callbacks.shouldInterrupt && callbacks.shouldInterrupt()) {
          wasInterrupted = true;
        }
        break;
      }
    }
#endif

    if (this->config.simplificationMode ==
        StructureSimplificationMode::NO_SIMPLIFICATION) {
      if (this->config.springBreakingLength > 0.0) {
        size_t nBrokenStrands = 0;
        const std::vector<size_t> brokenSpringIndices =
          net.findTooLongStrandSprings(
            u, this->config.springBreakingLength, is2D, &nBrokenStrands);
        nRemovedInIteration += nBrokenStrands;
        if (!brokenSpringIndices.empty()) {
          const std::vector<BrokenSpringInfo> brokenUpdates =
            net.removeSpringsAndCollectMatrixUpdates(
              brokenSpringIndices, kappa, &constants, is2D);
          if (!brokenUpdates.empty()) {
            cumulativeBrokenSinceRebuild += brokenUpdates.size();
            const bool shouldRebuildDueToThreshold =
              (this->config.inPlaceRebuildThreshold >= 0.0 &&
               totalSpringsAtLastRebuild > 0 &&
               (static_cast<double>(cumulativeBrokenSinceRebuild) /
                static_cast<double>(totalSpringsAtLastRebuild)) >=
                 this->config.inPlaceRebuildThreshold);
            if (shouldRebuildDueToThreshold) {
              matrixNeedsRebuild = true;
            } else {
              applyInPlaceMatrixUpdates(sysMatrix, brokenUpdates);
              matrixNeedsRebuild = false;
            }
          }
        }
      }
    } else {
      nRemovedInIteration = net.simplify(u,
                                         this->config.simplificationMode,
                                         this->config.inactiveRemovalCutoff,
                                         is2D,
                                         this->config.springBreakingLength);
      if (nRemovedInIteration > 0) {
        net.validate(u);
        matrixNeedsRebuild = true;
      }
    }

    initialSolution = net.coordinates + u;

    if (callbacks.onStep) {
      callbacks.onStep(iterationsDone);
    }

    if (callbacks.shouldInterrupt && callbacks.shouldInterrupt()) {
      wasInterrupted = true;
      break;
    }
  } while (net.nrOfStrands > 0 && nRemovedInIteration > 0);

  net.validate(u);

  result.effectiveSLESolver = effectiveSolverChoice;
  result.nrOfIterations = iterationsDone;

  if (wasInterrupted) {
    result.exitReason = ExitReason::INTERRUPT;
    if (callbacks.cleanupInterrupt) {
      callbacks.cleanupInterrupt();
    }
  } else if (result.exitReason == ExitReason::UNSET) {
    result.exitReason = ExitReason::X_TOLERANCE;
  }

  result.residual = computeResidualNorm(
    net, u, forceEvaluator, is2D, gradScratch, &result.energy);

  return result;
}

} // namespace pylimer_tools::sim::mehp
