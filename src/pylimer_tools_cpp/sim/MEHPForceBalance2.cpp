#include "MEHPForceBalance2.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
// #include "../utils/MemoryUtil.h"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// #ifndef NDEBUG
// #define DEBUG_REMOVAL
// #endif

namespace pylimer_tools {
namespace sim {
  namespace mehp {
#ifndef CLAMP_ONE_OVER_SPRINGPARTITION
    /**
     * @brief a macro for doing the clamping in the routines using kappa,
     * to prevent deivision by zero issues / multiplications by infinity
     */
#define CLAMP_ONE_OVER_SPRINGPARTITION(                                        \
  isPartialSpring, val, N, oneOverSpringPartitionUpperLimit)                   \
  ((!isPartialSpring)                                                          \
     ? val                                                                     \
     : std::clamp(val,                                                         \
                  (oneOverSpringPartitionUpperLimit > 0.)                      \
                    ? (1. / (N - 1. / oneOverSpringPartitionUpperLimit))       \
                    : (0.0),                                                   \
                  (oneOverSpringPartitionUpperLimit > 0.)                      \
                    ? (oneOverSpringPartitionUpperLimit)                       \
                    : (N)));
#endif

    /**
     * FORCE RELAXATION
     */
    void MEHPForceBalance2::runForceRelaxation(
      const StructureSimplificationMode simplificationMode,
      const double inactiveRemovalCutoff,
      const double oneOverSpringPartitionUpperLimit,
      const SLESolver solverChoice,
      const std::function<bool()>& shouldInterrupt,
      const std::function<void()>& cleanupInterrupt)
    {
      assert(this->validateNetwork());
      // INVALIDARG_EXP_IFN(
      //   shouldRemoveInactiveCrosslinks == false &&
      //     remove2functionalCrosslinkers == true,
      //   "Removing 2-functional crosslinkers only makes sense when inactive "
      //   "crosslinkers may be removed too, during the procedure.");
      this->simulationHasRun = true;

      INVALIDARG_EXP_IFN(
        inactiveRemovalCutoff > 0.0 ||
          simplificationMode == StructureSimplificationMode::NO_SIMPLIFICATION,
        "Removal cut-off must be positive when simplification is enabled.");

      if (this->getNrOfSprings() == 0) {
        return;
      }

      /* array allocation */
      double maxDistanceMoved = 0.0;
      size_t indexOfMaxDistanceMoved = 0;

      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(this->initialConfig,

                                             oneOverSpringPartitionUpperLimit);
      const double initialResidual =
        this->getDisplacementResidualNormFor(this->initialConfig,
                                             this->currentDisplacements,
                                             oneOverSpringPartitions);
      const double minN = this->initialConfig.springsContourLength.minCoeff();
      std::cout << "Starting force balance procedure "
                << "with " << initialResidual << " as initial residual."
                << std::endl;
      std::cout << "Simplification mode is " << simplificationMode << std::endl;
      std::cout << "Using oneOverSpringPartitionUpperLimit = "
                << oneOverSpringPartitionUpperLimit << std::endl;
      double currentResidual = initialResidual;
      double previousResidual = initialResidual;
      double intermediateResidual = 0.0;
      size_t iterationsDone = 0;

      this->prepareAllOutputs();

      // actual loop
      bool wasInterrupted = false;
      size_t nRemovedInIteration;
      do {
        nRemovedInIteration = 0;
        std::vector<Eigen::Triplet<double>> triplets;
        // diagonal + the lower of the two components of each spring
        triplets.reserve(this->initialConfig.nrOfLinks * 3 +
                         this->initialConfig.nrOfPartialSprings * 3 * 2);
        Eigen::VectorXd constants =
          Eigen::VectorXd::Zero(this->initialConfig.nrOfLinks * 3);
        // it's a bit more efficient to sum the diagonal ourselves
        Eigen::VectorXd diagonal =
          Eigen::VectorXd::Zero(this->initialConfig.nrOfLinks * 3);

        for (size_t springIdx = 0;
             springIdx < this->initialConfig.nrOfPartialSprings;
             ++springIdx) {
          if (this->initialConfig.springPartIndexA[springIdx] ==
              this->initialConfig.springPartIndexB[springIdx]) {
            continue;
          }
          const double contourLengthFraction =
            this->currentSpringPartitionsVec[springIdx];
          const double N =
            this->initialConfig.springsContourLength
              [this->initialConfig.partialToFullSpringIndex[springIdx]];
          double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
            this->initialConfig.partialSpringIsPartial[springIdx],
            1.0 / (N * contourLengthFraction),
            N,
            oneOverSpringPartitionUpperLimit);
          double multiplier =
            this->kappa *
            oneOverContourLengthFraction; // oneOverSpringPartitions(springIdx
          // * 3);
          assert(this->initialConfig.springPartIndexA[springIdx] <
                 this->initialConfig.nrOfLinks);
          assert(this->initialConfig.springPartIndexB[springIdx] <
                 this->initialConfig.nrOfLinks);
          // triplets will be summed up -> we can use the same indices multiple
          // times
          for (size_t dir = 0; dir < 3; ++dir) {
            // store only the lower part
            // if (this->initialConfig.springPartIndexA[springIdx] <
            //     this->initialConfig.springPartIndexB[springIdx]) {
            triplets.push_back(Eigen::Triplet<double>(
              this->initialConfig.springPartIndexA[springIdx] * 3 + dir,
              this->initialConfig.springPartIndexB[springIdx] * 3 + dir,
              1. * multiplier));
            // } else {
            triplets.push_back(Eigen::Triplet<double>(
              this->initialConfig.springPartIndexB[springIdx] * 3 + dir,
              this->initialConfig.springPartIndexA[springIdx] * 3 + dir,
              1. * multiplier));
            // }
          }
          diagonal.segment(3 * this->initialConfig.springPartIndexA[springIdx],
                           3) -= Eigen::Vector3d::Constant(multiplier);
          diagonal.segment(3 * this->initialConfig.springPartIndexB[springIdx],
                           3) -= Eigen::Vector3d::Constant(multiplier);

          // the constants, b in Ax = b
          constants.segment(3 * this->initialConfig.springPartIndexA[springIdx],
                            3) -=
            this->initialConfig.springPartBoxOffset.segment(3 * springIdx, 3) *
            multiplier;
          constants.segment(3 * this->initialConfig.springPartIndexB[springIdx],
                            3) +=
            this->initialConfig.springPartBoxOffset.segment(3 * springIdx, 3) *
            multiplier;
        }

        for (size_t linkIdx = 0; linkIdx < this->initialConfig.nrOfLinks;
             ++linkIdx) {
          for (size_t dir = 0; dir < 3; ++dir) {
            triplets.push_back(
              Eigen::Triplet<double>(linkIdx * 3 + dir,
                                     linkIdx * 3 + dir,
                                     diagonal(linkIdx * 3 + dir)));
          }
        }

        Eigen::SparseMatrix<double> sysMatrix(
          this->initialConfig.nrOfLinks * 3, this->initialConfig.nrOfLinks * 3);
        sysMatrix.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::VectorXd finalCoordinates;

#define COMPUTE_SOLVE(solver)                                                  \
  solver.compute(sysMatrix);                                                   \
  RUNTIME_EXP_IFN(solver.info() == Eigen::Success,                             \
                  "System matrix computation failed. Solver info: " +          \
                    std::to_string(solver.info()));                            \
  finalCoordinates = solver.solve(constants);                                  \
  RUNTIME_EXP_IFN(solver.info() == Eigen::Success,                             \
                  "System could not be solved. Solver info: " +                \
                    std::to_string(solver.info()));

#define ADD_ITERATIONS(solver) iterationsDone += solver.iterations();

#define SOLVE_ITERATIVE(solver)                                                \
  COMPUTE_SOLVE(solver);                                                       \
  ADD_ITERATIONS(solver);                                                      \
  break;

#define SOLVE_DIRECT(solver)                                                   \
  COMPUTE_SOLVE(solver);                                                       \
  break;

        switch (solverChoice) {
          // iterative solvers
          case SLESolver::DEFAULT:
          case SLESolver::CONJUGATE_GRADIENT:
          case SLESolver::CONJUGATE_GRADIENT_DIAGONALIZED: {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                     Eigen::Lower,
                                     Eigen::DiagonalPreconditioner<double>>
              solver;
            SOLVE_ITERATIVE(solver);
          }
          case SLESolver::CONJUGATE_GRADIENT_IDENTITY: {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                     Eigen::Lower,
                                     Eigen::IdentityPreconditioner>
              solver;
            SOLVE_ITERATIVE(solver);
          }
          case SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT: {
            Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>
              solver;
            SOLVE_ITERATIVE(solver);
          }
          case SLESolver::BICGSTAB: {
            Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
            SOLVE_ITERATIVE(solver);
          }
          // direct solvers
          case SLESolver::SIMPLICIAL_LLT: {
            Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
            SOLVE_DIRECT(solver);
          }
          case SLESolver::SIMPLICIAL_DLT: {
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            SOLVE_DIRECT(solver);
          }
          case SLESolver::SPARSE_LU: {
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            SOLVE_DIRECT(solver);
          }
          case SLESolver::SPARSE_QR: {
            Eigen::SparseQR<Eigen::SparseMatrix<double>,
                            Eigen::COLAMDOrdering<int>>
              solver;
            SOLVE_DIRECT(solver);
          }
          default:
            throw std::runtime_error("This solver is not implemented.");
        }

        this->currentDisplacements =
          finalCoordinates - this->initialConfig.coordinates;

        currentResidual = this->getDisplacementResidualNormFor(
          this->initialConfig,
          this->currentDisplacements,

          oneOverSpringPartitionUpperLimit);

        previousResidual = currentResidual;
        iterationsDone += 1;
        if (iterationsDone % this->simplificationFrequency == 0) {
          this->breakTooLongSprings(this->initialConfig,
                                    this->currentDisplacements,
                                    this->currentSpringPartitionsVec);
          size_t nRemovedThisLoop = 0;

          do {
#ifndef NDEBUG
            this->validateNetwork();
#endif
            nRemovedThisLoop = 0;
            if (simplificationMode ==
                  StructureSimplificationMode::INACTIVE_ONLY ||
                simplificationMode == StructureSimplificationMode::ALL_TIM) {
#ifdef DEBUG_REMOVAL
              std::cout << "Checking and possibly removing inactive cross-links"
                        << std::endl;
#endif
              nRemovedThisLoop +=
                this->removeInactiveCrosslinks(this->initialConfig,
                                               this->currentDisplacements,

                                               inactiveRemovalCutoff);
            }
            if (simplificationMode == StructureSimplificationMode::X2F_ONLY ||
                simplificationMode == StructureSimplificationMode::ALL_TIM) {
#ifdef DEBUG_REMOVAL
              std::cout
                << "Checking and possibly removing cross-links with f = 2"
                << std::endl;
#endif
              nRemovedThisLoop += this->removeTwofunctionalCrosslinks(
                this->initialConfig,
                this->currentDisplacements,
                this->currentSpringPartitionsVec);
            }
            if (simplificationMode == StructureSimplificationMode::ALL_ANDREI) {
#ifdef DEBUG_REMOVAL
              std::cout << "Checking and possibly removing cross-links and "
                           "springs, Andrei's way"
                        << std::endl;
#endif
              nRemovedThisLoop +=
                this->doRemovalAndreisWay(this->initialConfig,
                                          this->currentDisplacements,

                                          inactiveRemovalCutoff);
            }

            // cleanup some things
            if (simplificationMode !=
                StructureSimplificationMode::NO_SIMPLIFICATION) {
              this->validateNetwork(this->initialConfig,
                                    this->currentDisplacements,
                                    this->currentSpringPartitionsVec);
              oneOverSpringPartitions = this->assembleOneOverSpringPartition(
                this->initialConfig,

                oneOverSpringPartitionUpperLimit);
            }

            nRemovedInIteration += nRemovedThisLoop;
          } while (nRemovedThisLoop > 0);
          // after removal, the residual changed, might even have increased
          // beyond initial
          // -> reset previous and current to prevent change to iterative
          // displacement
          if (nRemovedInIteration > 0) {
            previousResidual = initialResidual;

            oneOverSpringPartitions = this->assembleOneOverSpringPartition(
              this->initialConfig,

              oneOverSpringPartitionUpperLimit);
          }
        }
        this->handleOutput(iterationsDone);

        if (shouldInterrupt()) {
          wasInterrupted = true;
          break;
        }
      } while (this->initialConfig.nrOfSprings > 0 && nRemovedInIteration > 0);

      // finish up
      this->closeAllOutputs();

      // query solution & exit reason
      this->exitReason = ExitReason::X_TOLERANCE;
      this->nrOfStepsDone += iterationsDone;
      std::cout << iterationsDone << " steps done. "
                << "Last max distance moved: " << maxDistanceMoved << ". "
                << "Current residual: " << currentResidual << ". "
                << "Initial residual: " << initialResidual << ". " << std::endl;

      assert(this->currentDisplacements.size() ==
             3 * this->initialConfig.nrOfLinks);
      this->validateNetwork();
      this->currentSpringVectors = this->evaluateSpringVectors(
        this->initialConfig, this->currentDisplacements);
      this->currentPartialSpringVectors = this->evaluatePartialSpringVectors(
        this->initialConfig, this->currentDisplacements);
      if (wasInterrupted) {
        this->exitReason = ExitReason::INTERRUPT;
        cleanupInterrupt();
      }
    }

    /**
     * @brief Compute the displacement residual norm for the current
     * configuration
     *
     * @param oneOverSpringPartitionUpperLimit
     * @return double
     */
    double MEHPForceBalance2::getDisplacementResidualNorm(
      double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::VectorXd oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(this->initialConfig,

                                             oneOverSpringPartitionUpperLimit);
      Eigen::VectorXd displacements = this->currentDisplacements;
      return this->getDisplacementResidualNormFor(
        this->initialConfig, displacements, oneOverSpringPartitions);
    }

    /**
     * @brief Compute the displacement residual norm for a specific
     * configuration
     *
     * @param net
     * @param u
     * @param oneOverSpringPartitions
     * @return double
     */
    double MEHPForceBalance2::getDisplacementResidualNormFor(
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,

      const double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::ArrayXd oneOverSpringPartitions =
        this
          ->assembleOneOverSpringPartition(net,

                                           oneOverSpringPartitionUpperLimit)
          .array();

      Eigen::ArrayXd loopPartialSpringEliminator =
        (net.springPartCoordinateIndexA != net.springPartCoordinateIndexB)
          .cast<double>();
      Eigen::ArrayXd forces = Eigen::ArrayXd::Zero(3 * net.nrOfLinks);
      Eigen::ArrayXd distances =
        this
          ->evaluatePartialSpringVectors(
            net, u, this->is2D, this->assumeBoxLargeEnough)
          .array();
      forces(net.springPartCoordinateIndexA) +=
        (this->kappa * oneOverSpringPartitions * distances *
         loopPartialSpringEliminator);
      forces(net.springPartCoordinateIndexB) -=
        (this->kappa * oneOverSpringPartitions * distances *
         loopPartialSpringEliminator);

      // #ifndef NDEBUG
      //       Eigen::VectorXi debugNrSpringsVisited =
      //         Eigen::VectorXi::Zero(net.nrOfPartialSprings);
      //       for (size_t i = 0; i < net.nrOfLinks; ++i) {
      //         Eigen::Array3d forces2 =
      //           this
      //             ->evaluateForceOnLink(i,
      //                                   net,
      //                                   u,
      //
      //                                   debugNrSpringsVisited,
      //                                   oneOverSpringPartitionUpperLimit)
      //             .array();
      //         double squareN1 = forces.segment(3 * i,
      //         3).matrix().squaredNorm(); double squareN2 =
      //         forces2.matrix().squaredNorm(); if (i == 129) {
      //           this->debugAtomVicinity(net.oldAtomIds[i]);
      //         }
      //         assert(APPROX_EQUAL(squareN1, squareN2, 1e-9));
      //         // assert(pylimer_tools::utils::vector_approx_equal(
      //         //   forces.segment(3 * i, 3), forces2, 1e-9, true));
      //       }
      //       assert((debugNrSpringsVisited.array() == 2).all());
      // #endif

      return forces.matrix().squaredNorm();
    }

    /**
     * @brief Compute the displacement residual norm for a specific
     * configuration
     *
     * @param net
     * @param u
     * @param oneOverSpringPartitions
     * @return double
     */
    double MEHPForceBalance2::getDisplacementResidualNormFor(
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,
      const Eigen::VectorXd& oneOverSpringPartitions) const
    {
      Eigen::VectorXd displacedCoords = net.coordinates + u;

      Eigen::VectorXd relevantPartialDistances =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;

      if (this->assumeBoxLargeEnough) {
        this->box.handlePBC(relevantPartialDistances);
      }

      if (this->is2D) {
        for (size_t i = 2; i < relevantPartialDistances.size(); i += 3) {
          relevantPartialDistances[i] = 0.;
        }
      }

#ifndef NDEBUG
      for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
        Eigen::Vector3d dist = this->evaluatePartialSpringDistance(net, u, i);
        Eigen::Vector3d comparison = relevantPartialDistances.segment(3 * i, 3);
        assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
          dist, comparison, 1e-9));
      }
#endif

      assert(relevantPartialDistances.size() == oneOverSpringPartitions.size());
      Eigen::VectorXd partialDistancesOverSpringPartitions =
        (relevantPartialDistances.array() * oneOverSpringPartitions.array())
          .matrix();

      Eigen::VectorXd overallForces = Eigen::VectorXd::Zero(3 * net.nrOfLinks);
      overallForces(net.springPartCoordinateIndexB) -=
        partialDistancesOverSpringPartitions;
      overallForces(net.springPartCoordinateIndexA) +=
        partialDistancesOverSpringPartitions;

      return overallForces.squaredNorm();
    }

    /**
     * @brief Translate the spring partition vector to its 3*size
     *
     * @param net
     * @param springPartitions0
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance2::assembleOneOverSpringPartition(
      const ForceBalance2Network& net,
      0,
      const double oneOverSpringPartitionUpperLimit) const
    {
      INVALIDARG_EXP_IFN(
        springPartitions0.size() == net.nrOfPartialSprings,
        "Spring partitions must have the size of the nr of springs");
      Eigen::VectorXd oneOverSpringPartitions =
        Eigen::VectorXd(3 * net.nrOfPartialSprings);

      Eigen::ArrayXd primaryLoopCorrectionMultiplier =
        (net.springPartIndexA != net.springPartIndexB)
          .cast<double>(); // 0.0 for equal = primary loop, 1.0 otherwise

      for (size_t i = 0; i < net.nrOfPartialSprings; ++i) {
        const double N =
          net.springsContourLength[net.partialToFullSpringIndex[i]];
        double contourLengthFraction = springPartitions0[i];
        double valueToSet =
          CLAMP_ONE_OVER_SPRINGPARTITION(net.partialSpringIsPartial[i],
                                         1.0 / (N * contourLengthFraction),
                                         N,
                                         oneOverSpringPartitionUpperLimit);

        // if (springPartitions0[i] < 1e-9) {
        //   std::cout << "Got close call for partial spring " << i <<
        //   std::endl;
        // }
        oneOverSpringPartitions.segment(3 * i, 3) = Eigen::Vector3d::Constant(
          valueToSet * primaryLoopCorrectionMultiplier[i]);
      }

      return oneOverSpringPartitions;
    }

    /**
     * @brief Remove double listed springs from cross-linkers
     *
     * @param net
     */
    void MEHPForceBalance2::removeDuplicateListedSpringsFromLinks(
      ForceBalance2Network& net) const
    {
      for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
        this->removeDuplicateListedSpringsFromLink(net, linkIdx);
      }

#ifndef NDEBUG
      this->validateNetwork();
#endif
    }

    void MEHPForceBalance2::removeDuplicateListedSpringsFromLink(
      ForceBalance2Network& net,
      size_t linkIdx,
      bool allowOnEntanglement) const
    {
      INVALIDARG_EXP_IFN(linkIdx < net.nrOfLinks,
                         "Cannot remove duplicate spring indices of index "
                         "higher than nr. of links.");
      // remove duplicate mentions of the same spring index
      std::sort(net.springIndicesOfLinks[linkIdx].begin(),
                net.springIndicesOfLinks[linkIdx].end());
      auto last = std::unique(net.springIndicesOfLinks[linkIdx].begin(),
                              net.springIndicesOfLinks[linkIdx].end());
      if (last != net.springIndicesOfLinks[linkIdx].end()) {
#ifdef DEBUG_REMOVAL
        std::cout << "Removed duplicate spring indices from link " << linkIdx
                  << std::endl;
#endif
        if (!allowOnEntanglement) {
          if (linkIdx < net.nrOfNodes) {
            RUNTIME_EXP_IFN(
              net.oldAtomTypes[linkIdx] != this->entanglementType,
              "Require entanglement beads to not form primary loops. Link " +
                std::to_string(linkIdx) + " is slip-link " +
                std::to_string(net.linkIsEntanglement[linkIdx]) +
                " with old atom type" +
                std::to_string(net.oldAtomTypes[linkIdx]) + ".");
          } else {
            RUNTIME_EXP_IFN(
              !net.linkIsEntanglement[linkIdx],
              "Require entanglement beads to not form primary loops. Link " +
                std::to_string(linkIdx) + " is slip-link " +
                std::to_string(net.linkIsEntanglement[linkIdx]) + ".");
          }
        }
        net.springIndicesOfLinks[linkIdx].erase(
          last, net.springIndicesOfLinks[linkIdx].end());
      }
    }

    size_t MEHPForceBalance2::removePrimaryLoops(
      ForceBalance2Network& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions) const
    {
      size_t numRemoved = 0;
      for (long int springIdx = net.nrOfSprings - 1; springIdx >= 0;
           --springIdx) {
        if (net.springIndexA[springIdx] == net.springIndexB[springIdx] &&
            net.localToGlobalSpringIndex[springIdx].size() == 1 &&
            net.springsType[springIdx] != this->entanglementType) {
          this->removeSpringFollowingEntanglementLinks(net,
                                                       displacements,

                                                       springIdx);
          springIdx = std::min<long int>(springIdx, net.nrOfSprings - 1);
          numRemoved += 1;
        }
      }

#ifndef NDEBUG
      this->validateNetwork();
#endif
      return numRemoved;
    }

    /**
     * @brief Remove crosslinkers which do not have any springs with a certain
     * minimum length
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param tolerance
     */
    size_t MEHPForceBalance2::removeInactiveCrosslinks(
      ForceBalance2Network& net,
      Eigen::VectorXd& displacements,

      double tolerance) const
    {
      size_t numRemoved = 0;
      //        this->removePrimaryLoops(net, displacements, springPartitions);
      // this->validateNetwork(net, displacements, springPartitions);
      // first, we remove all inactive springs
      for (long int springIdx = net.nrOfSprings - 1; springIdx >= 0;
           --springIdx) {
        if (springIdx >= net.nrOfSprings) {
          springIdx = net.nrOfSprings - 1;
          continue;
        }
        if (net.springsType[springIdx] == this->entanglementType) {
          // let's not remove entanglement springs
          continue;
        }
        if (net.oldAtomTypes[net.springIndexA[springIdx]] ==
              this->entanglementType &&
            net.oldAtomTypes[net.springIndexB[springIdx]] ==
              this->entanglementType) {
          // this is a quasi-partial spring, will be handled differently
          continue;
        }
        std::vector<size_t> involvedPartialSprings =
          this->getAllPartialSpringIndicesAlong(net, springIdx);
        assert(involvedPartialSprings.size() > 0);
        bool isActive = false;
        for (size_t partialSpringIdx : involvedPartialSprings) {
          RUNTIME_EXP_IFN(
            net.coordinates.size() == displacements.size(),
            "Expected coordinates and displacements to have same size, got " +
              std::to_string(net.coordinates.size()) + " and " +
              std::to_string(displacements.size()) + ".");
          // assert(net.coordinates.size() == displacements.size());
          Eigen::Vector3d distance = this->evaluatePartialSpringDistance(
            net, displacements, partialSpringIdx);
          double contourLength =
            net.springsContourLength
              [net.partialToFullSpringIndex[partialSpringIdx]];
          double partition = springPartitions[partialSpringIdx];
          if (!this->distanceIsWithinTolerance(
                distance, tolerance, contourLength, partition)) {
            isActive = true;
            break;
          }
        }
        if (!isActive) {
          // remove this spring
#ifdef DEBUG_REMOVAL
          std::cout << "Removing inactive spring " << springIdx
                    << " with all dependencies" << std::endl;
#endif
          this->removeSpringFollowingEntanglementLinks(net,
                                                       displacements,

                                                       springIdx);

#ifndef NDEBUG
          this->validateNetwork(net, displacements, springPartitions);
#endif
          numRemoved += 1;
        }
      }

      // then, we remove all crosslinkers that are 0- or 1-functional
      for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
           --crosslinkIdx) {
        assert(net.springIndicesOfLinks.size() > crosslinkIdx);
        if (net.springIndicesOfLinks[crosslinkIdx].size() == 0 // f = 0
        ) {
#ifdef DEBUG_REMOVAL
          std::cout << "Removing f = 0 x-link " << crosslinkIdx << std::endl;
#endif

          this->removeLink(net, displacements, crosslinkIdx);
          numRemoved += 1;
#ifndef NDEBUG
          this->validateNetwork(net, displacements, springPartitions);
#endif
        } else if ( // or f = 1, NOT primary loop
          (net.springIndicesOfLinks[crosslinkIdx].size() == 1) &&
          (XOR(
            net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx][0]]
                                    [0] == crosslinkIdx,
            pylimer_tools::utils::last(
              net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx]
                                                               [0]]) ==
              crosslinkIdx))) {
          assert(net.oldAtomTypes[crosslinkIdx] != this->entanglementType);
#ifdef DEBUG_REMOVAL
          std::cout << "Removing f = 1 x-link " << crosslinkIdx << std::endl;
#endif
          std::vector<size_t> entanglementLinksRemoved =
            this->getEntanglementLinkIndicesAlong(
              net, net.springIndicesOfLinks[crosslinkIdx][0]);
          // need to first remove the spring
          const long int previousId = net.oldAtomIds[crosslinkIdx];
          this->removeSpringFollowingEntanglementLinks(
            net,
            displacements,

            net.springIndicesOfLinks[crosslinkIdx][0]);
          // this also removed some intermediate entanglement links
          // -> crossLinkIdx is now at a different index
          size_t newCrosslinkIdx = crosslinkIdx;
          for (size_t linkIdx : entanglementLinksRemoved) {
            assert(linkIdx != crosslinkIdx);
            if (linkIdx < crosslinkIdx) {
              newCrosslinkIdx -= 1;
            }
          }
          assert(net.oldAtomIds[newCrosslinkIdx] == previousId);
          assert(net.oldAtomTypes[newCrosslinkIdx] != this->entanglementType);
          assert(net.springIndicesOfLinks[newCrosslinkIdx].size() == 0);
          numRemoved += 1;
          // to then remove the cross-link
          this->removeLink(net, displacements, newCrosslinkIdx);
#ifdef DEBUG_REMOVAL
          std::cout << "Effectively removed f = 1 x-link " << newCrosslinkIdx
                    << std::endl;
#endif
          crosslinkIdx = std::min<long int>(crosslinkIdx, net.nrOfNodes - 1);
          // => we should ever only have 2-functional entanglement links that
          // could be merged after this.

#ifndef NDEBUG
          this->validateNetwork(net, displacements, springPartitions);
#endif
        }
      }

#ifndef NDEBUG
      this->validateNetwork(net, displacements, springPartitions);
#endif

      return numRemoved;
    }

    /**
     * @brief Remove springs that exert a stress higher than
     * `this->springBreakingLength`
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @return size_t the number of springs broken
     */
    size_t MEHPForceBalance2::breakTooLongSprings(
      ForceBalance2Network& net,
      Eigen::VectorXd& displacements,
      Eigen::VectorXd& springPartitions) const
    {
      if (this->springBreakingLength <= 0.) {
        return 0;
      }

      size_t numBroken = 0;

      // iterate the springs, determine their distance, and determine if it
      // exceeds the breaking force
      for (long int partialSpringIdx = net.nrOfPartialSprings;
           partialSpringIdx >= 0;
           --partialSpringIdx) {
        if (partialSpringIdx >= net.nrOfPartialSprings) {
          partialSpringIdx = net.nrOfPartialSprings - 1;
        }
        double len = this->getWeightedPartialSpringLength(net,
                                                          displacements,

                                                          partialSpringIdx);
        if (len > this->springBreakingLength) {
          // break this spring
          numBroken += 1;
          this->breakPartialSpring(net,
                                   displacements,

                                   partialSpringIdx);
        }
      }

      return numBroken;
    }

    /**
     * @brief Add slip-links to this system
     *
     * @param sliplinkDensity
     * @param cutoff
     */
    size_t MEHPForceBalance2::randomlyAddSliplinks(
      const size_t nrOfSliplinksToSample,
      const double cutoff,
      const size_t minimumNrOfSliplinks,
      const double sameStrandCutoff,
      const bool excludeCrosslinks,
      const int seed)
    {
      INVALIDARG_EXP_IFN(nrOfSliplinksToSample > minimumNrOfSliplinks,
                         "Maximum nr. should be larger than minimum, got " +
                           std::to_string(nrOfSliplinksToSample) + " and " +
                           std::to_string(minimumNrOfSliplinks) + ".");
      INVALIDARG_EXP_IFN(cutoff > 0.0,
                         "Expected a cutoff > 0.0, got " +
                           std::to_string(cutoff) + ".");
      // RUNTIME_EXP_IFN(this->initialConfig.nrOfLinks ==
      //                   this->initialConfig.nrOfNodes,
      //                 "Slip-links are only added randomly when no other "
      //                 "slip-links are in place yet.");
      INVALIDARG_EXP_IFN(minimumNrOfSliplinks <
                           this->universe.getNrOfAtoms() / 2,
                         "Minimum number of slip-links must be less than the "
                         "possible number of slip-links to place.");
      INVALIDARG_EXP_IFN(nrOfSliplinksToSample <
                           this->universe.getNrOfAtoms() / 2,
                         "Number of slip-links to place must be less than "
                         "the possible number of slip-links to place.");
      // query all the cross-linker chains we actually use to place
      // slip-links on
      std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
        this->universe.getChainsWithCrosslinker(crossLinkerType);
      bool danglingChainsAreKept =
        this->initialConfig.nrOfNodes >
        this->universe.getAtomsOfType(crossLinkerType).size();
      // and also query all the corresponding atoms we use to place slip-links
      // on
      size_t nrOfEligibleAtoms = 0;
      std::vector<pylimer_tools::entities::Atom> eligibleAtoms;
      eligibleAtoms.reserve(this->universe.getNrOfAtoms());
      std::vector<bool> vertexIdxIsEligible =
        pylimer_tools::utils::initializeWithValue<bool>(
          this->universe.getNrOfAtoms(), false);
      std::unordered_map<size_t, size_t> atomToStrand;
      atomToStrand.reserve(this->universe.getNrOfAtoms());
      std::unordered_map<size_t, size_t> atomIdxInStrand;
      atomIdxInStrand.reserve(this->universe.getNrOfAtoms());
      size_t springId = 0;
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        pylimer_tools::entities::Molecule chain = crossLinkerChains[i];
        RUNTIME_EXP_IFN(chain.getType() !=
                          pylimer_tools::entities::MoleculeType::UNDEFINED,
                        "Couldn't determine molecule type.");
        if (chain.getType() ==
              pylimer_tools::entities::MoleculeType::PRIMARY_LOOP ||
            chain.getType() ==
              pylimer_tools::entities::MoleculeType::NETWORK_STRAND ||
            (chain.getType() ==
               pylimer_tools::entities::MoleculeType::DANGLING_CHAIN &&
             danglingChainsAreKept)) {
          assert(i == this->initialConfig.springToMoleculeIds[springId]);
          // TODO: also check that this is not a higher order dangling strand
          nrOfEligibleAtoms +=
            crossLinkerChains[i].getNrOfAtoms() -
            crossLinkerChains[i].getAtomsOfType(this->crossLinkerType).size();
          std::vector<pylimer_tools::entities::Atom> atoms =
            crossLinkerChains[i].getAtomsLinedUp(this->crossLinkerType);
          for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
            pylimer_tools::entities::Atom atom = atoms[atomIdx];
            if (atom.getType() != this->crossLinkerType) {
              eligibleAtoms.push_back(atom);
              vertexIdxIsEligible[this->universe.getIdxByAtomId(atom.getId())] =
                true;
              atomToStrand.emplace(atom.getId(), springId);
              atomIdxInStrand.emplace(atom.getId(), atomIdx);
              RUNTIME_EXP_IFN(this->initialConfig.oldAtomIdToSpringIndex.at(
                                atom.getId()) == springId,
                              "The spring numbering seems incorrect. Placing "
                              "slip-links will "
                              "lead to inappropriate placement.");
            }
          }
          RUNTIME_EXP_IFN(this->initialConfig.springToMoleculeIds[springId] ==
                            i,
                          "The spring numbering seems incorrect. Placing "
                          "slip-links will lead to inappropriate placement.");
          // #ifndef NDEBUG
          std::vector<pylimer_tools::entities::Atom> chainEnds =
            crossLinkerChains[i].getChainEnds(this->crossLinkerType, true);
          RUNTIME_EXP_IFN(
            this->initialConfig.oldAtomIds
                [this->initialConfig.linkIndicesOfSprings[springId][0]] ==
              chainEnds[0].getId(),
            "Atom ends are inconsistent when sampling.");
          RUNTIME_EXP_IFN(chainEnds[0].getId() == atoms[0].getId(),
                          "Atom ends are inconsistent when sampling.");
          // #endif
          springId += 1;
        }
      }
      // build neighbourlist
      std::vector<pylimer_tools::entities::Atom> atomsForNeighbourList =
        this->universe.getAtoms();
      std::vector<bool> isMasked = pylimer_tools::utils::initializeWithValue(
        this->universe.getNrOfAtoms(), false);
      if (excludeCrosslinks) {
        // TODO: check whether it is faster to just only query the other ones
        std::erase_if(atomsForNeighbourList,
                      [&](const pylimer_tools::entities::Atom& a) -> bool {
                        return a.getType() == this->crossLinkerType;
                      });
        for (size_t i = 0; i < this->universe.getNrOfAtoms(); ++i) {
          isMasked[i] = (this->universe.getAtomByVertexIdx(i).getType() ==
                         this->crossLinkerType);
        }
        for (pylimer_tools::entities::Atom a : atomsForNeighbourList) {
          RUNTIME_EXP_IFN(
            a.getType() != this->crossLinkerType,
            "Removing cross-linkers from neighborlist did not seem to work.");
        }
      }
      pylimer_tools::entities::NeighbourList neighbourList =
        pylimer_tools::entities::NeighbourList(
          atomsForNeighbourList, this->box, cutoff);

      std::random_device rd{};
      std::mt19937 rng = std::mt19937(seed > 0 ? seed : rd());
      // std::cout << "Initial sampling rng seed: " << rng << std::endl;
      // build list of random samples
      // this way is more performant than
      // sampling integers and checking whether they have been sampled already
      std::vector<size_t> toSampleFrom;
      toSampleFrom.reserve(this->universe.getNrOfAtoms());
      for (size_t i = 0; i < this->universe.getNrOfAtoms(); ++i) {
        toSampleFrom.push_back(i);
      }
      std::shuffle(toSampleFrom.begin(), toSampleFrom.end(), rng);

      // the resulting vectors to fill
      std::vector<double> slipLinkXs;
      slipLinkXs.reserve(nrOfSliplinksToSample);
      std::vector<double> slipLinkYs;
      slipLinkYs.reserve(nrOfSliplinksToSample);
      std::vector<double> slipLinkZs;
      slipLinkZs.reserve(nrOfSliplinksToSample);
      std::vector<size_t> slipLinkStrandA;
      slipLinkStrandA.reserve(nrOfSliplinksToSample);
      std::vector<size_t> slipLinkStrandB;
      slipLinkStrandB.reserve(nrOfSliplinksToSample);
      std::vector<double> slipLinkStrandAlpha;
      slipLinkStrandAlpha.reserve(nrOfSliplinksToSample);
      std::vector<double> slipLinkStrandBeta;
      slipLinkStrandBeta.reserve(nrOfSliplinksToSample);

      pylimer_tools::entities::Box box = this->box;

      size_t nrOfSlipLinksPlaced = 0;
      size_t nrOfAttempts = 0;
      size_t sampleIdx = 0;
      // the actual sampling loop
      while (nrOfSlipLinksPlaced < minimumNrOfSliplinks ||
             nrOfAttempts < nrOfSliplinksToSample) {
        // first, randomly sample an atom
        while (isMasked[toSampleFrom[sampleIdx]]) {
          sampleIdx += 1;

          if (sampleIdx >= toSampleFrom.size()) {
            break;
          }
        }
        if (sampleIdx >= toSampleFrom.size()) {
          // this is a path that should barely ever be reached
          std::cerr << "Sample index exceeds number of samples." << std::endl;
          break;
        }
        size_t sampledVertexId = toSampleFrom[sampleIdx];
        nrOfAttempts += 1;
        sampleIdx += 1;
        pylimer_tools::entities::Atom a1 =
          this->universe.getAtomByVertexIdx(sampledVertexId);
        RUNTIME_EXP_IFN(!excludeCrosslinks ||
                          a1.getType() != this->crossLinkerType,
                        "Sampled atom is cross-link, but may not be");
        isMasked[sampledVertexId] = true;
        // then, find neighbouring atoms (but not from the same strand?!)
        std::vector<pylimer_tools::entities::Atom> neighbours =
          neighbourList.getAtomsCloseTo(a1);
        neighbourList.removeAtom(a1,
                                 "After querying neighbours. Impossible case.");
        // filter the neighbours to include only those from other strands
        // NOTE: this skews the whole thing a bit
        neighbours.erase(
          std::remove_if(
            neighbours.begin(),
            neighbours.end(),
            [&](pylimer_tools::entities::Atom a) -> bool {
              return (
                (atomToStrand[a.getId()] ==
                   atomToStrand[a1.getId()] // do not use "at", because not
                 // all atoms in the neighbours
                 // have been assigned a strand
                 && (std::abs(static_cast<double>(
                       atomIdxInStrand[a.getId()] -
                       atomIdxInStrand[a1.getId()])) < sameStrandCutoff))
                // the following check should not be necessary?!?
                || isMasked[this->universe.getIdxByAtomId(a.getId())]);
            }),
          neighbours.end());
        if (neighbours.size() == 0) {
          std::cerr << "Not enough close neighbours found." << std::endl;
          continue;
        }
        // then, randomly select one of them
        pylimer_tools::entities::Atom a2 = neighbours[0];
        if (neighbours.size() > 1) {
          size_t randomA2Idx =
            std::uniform_int_distribution<size_t>{ 0,
                                                   neighbours.size() - 1 }(rng);
          a2 = neighbours[randomA2Idx];
        }
        // finally, remove them from the neighbour lists so that they are not
        // sampled more than once
        size_t sampledVertexId2 = this->universe.getIdxByAtomId(a2.getId());
        RUNTIME_EXP_IFN(sampledVertexId2 != sampledVertexId,
                        "Second sample may not be equal to the first");
        RUNTIME_EXP_IFN(!excludeCrosslinks ||
                          a2.getType() != this->crossLinkerType,
                        "Sampled atom is cross-link, but may not be");
        RUNTIME_EXP_IFN(!isMasked[sampledVertexId2],
                        "Sampled vertex 2 is masked, but could not be");
        neighbourList.removeAtom(a2, "Removing sampled vertex 2");
        isMasked[sampledVertexId2] = true;
        // it is actually quite a lot of expensive stuff done until we get to
        // this check but only this way we have the balance of removing atoms
        // to sample them only once
        if (!vertexIdxIsEligible[sampledVertexId] ||
            !vertexIdxIsEligible[sampledVertexId2]) {
          // std::cout << "Sampled vertices are not eligible" << std::endl;
          continue;
        }
        // take the mean and their index etc. to add as slip-link
        // std::cout << "OMerging atoms " << a1.getId() << " and " <<
        // a2.getId()
        // << " with distance " << a1.distanceTo(a2, universe.getBox()) <<
        // std::endl;
        Eigen::Vector3d meanPositions = a1.meanPositionWith(a2, box);
        slipLinkXs.push_back(meanPositions[0]);
        slipLinkYs.push_back(meanPositions[1]);
        slipLinkZs.push_back(meanPositions[2]);
        slipLinkStrandA.push_back(atomToStrand.at(a1.getId()));
        slipLinkStrandB.push_back(atomToStrand.at(a2.getId()));
        slipLinkStrandAlpha.push_back(
          static_cast<double>(atomIdxInStrand.at(a1.getId())) /
          (static_cast<double>(
            this->initialConfig
              .springsContourLength[atomToStrand.at(a1.getId())])));
        slipLinkStrandBeta.push_back(
          static_cast<double>(atomIdxInStrand.at(a2.getId())) /
          (static_cast<double>(
            this->initialConfig
              .springsContourLength[atomToStrand.at(a2.getId())])));
        nrOfSlipLinksPlaced += 1;
      }

      RUNTIME_EXP_IFN(
        slipLinkStrandA.size() == slipLinkStrandB.size() &&
          slipLinkXs.size() == slipLinkYs.size() &&
          slipLinkZs.size() == slipLinkXs.size(),
        "Expect all slip-link relevant properties to have the same length");
      // with the data assembled, we can actually add them to the structure
      // and stuff
      this->addSlipLinks(slipLinkStrandA,
                         slipLinkStrandB,
                         slipLinkXs,
                         slipLinkYs,
                         slipLinkZs,
                         slipLinkStrandAlpha,
                         slipLinkStrandBeta,
                         false);
      return nrOfSlipLinksPlaced;
    }

    /**
     * @brief Find cycles (loops) and set slip-links based on those
     *
     * @param maxLoopLength
     * @return size_t
     */
    size_t MEHPForceBalance2::addSliplinksBasedOnCycles(const int maxLoopLength)
    {
      // std::cout << "Detecting slip-links based on cycles. Base memory
      // useage:
      // "
      //           << getCurrentRSS() << ", peak " << getPeakRSS() <<
      //           std::endl;
      std::vector<std::vector<long int>> loopEdges;
      std::vector<std::vector<long int>> loops = this->universe.findLoops(
        this->crossLinkerType, maxLoopLength, false, &loopEdges);
      std::cout << "Detected " << loops.size() << " loops." << std::endl;
      // reduced loops = loops, but only the (new) spring indices
      std::vector<std::vector<size_t>> reducedLoops;
      reducedLoops.reserve(loops.size());
      // min & max coordinates of the loops in order to easily know
      // a priori whether it makes sense to compare two loops or not
      std::vector<std::array<double, 3>> loopMinCoords;
      loopMinCoords.reserve(loops.size());
      std::vector<std::array<double, 3>> loopMaxCoords;
      loopMaxCoords.reserve(loops.size());
      pylimer_tools::entities::Box box = this->box;
      for (std::vector<long int> loop : loops) {
        // max & min coordinates
        // TODO: implement, such that the max & min work also with infinite
        // loops
        // assume each subsequent atom is bonded with a bond shorter
        // than half the bond
        Eigen::VectorXd alignedLoopCoordinates =
          this->universe.getUnwrappedVertexCoordinates(loop, box);

        std::array<double, 3> minCoords;
        std::array<double, 3> maxCoords;

        for (int i = 0; i < 3; ++i) {
          minCoords[i] =
            alignedLoopCoordinates(Eigen::seq(i, Eigen::last - 2 + i, 3))
              .minCoeff();
          maxCoords[i] =
            alignedLoopCoordinates(Eigen::seq(i, Eigen::last - 2 + i, 3))
              .maxCoeff();
          // want to see later that this direction is spanned fully
          // that's why we just push it so far that it is within the box
          if (maxCoords[i] - minCoords[i] > box.getL(i) ||
              (std::fabs(alignedLoopCoordinates[i] -
                         alignedLoopCoordinates[3 * (loop.size() - 1) + i]) >
               0.5 * box.getL(i))) {
            // special case in case the loop is periodic
            minCoords[i] = box.getLowL(i);
            maxCoords[i] = box.getHighL(i);
          } else {
            // adjust for box
            while (minCoords[i] > box.getHighL(i)) {
              minCoords[i] -= box.getL(i);
            }
            while (minCoords[i] < box.getLowL(i)) {
              minCoords[i] += box.getL(i);
            }
            while (maxCoords[i] > box.getHighL(i)) {
              maxCoords[i] -= box.getL(i);
            }
            while (maxCoords[i] < box.getLowL(i)) {
              maxCoords[i] += box.getL(i);
            }
            // at this point, it is not given anymore that minCoords[i] <
            // maxCoords[i], instead, we know that both are within the box
          }
        }

        loopMinCoords.push_back(minCoords);
        loopMaxCoords.push_back(maxCoords);

        // loop reduction to springs
        std::set<size_t> reducedLoop;
        for (size_t i = 0; i < loop.size(); ++i) {
          // reduced loop
          if (this->universe.getPropertyValue<long int>("type", loop[i]) !=
              this->crossLinkerType) {
            // set takes care of duplicates, yet this is not efficient at all.
            reducedLoop.insert(
              this->initialConfig.oldAtomIdToSpringIndex
                [this->universe.getPropertyValue<long int>("id", loop[i])]);
          }
        }
        std::vector<size_t> reducedLoopVec(reducedLoop.begin(),
                                           reducedLoop.end());
        reducedLoops.push_back(reducedLoopVec);
      }
      // std::cout << "After finding loops, memory useage: " <<
      // getCurrentRSS()
      //           << ", peak " << getPeakRSS() << std::endl;

      // fetch some data to later estimate alpha & beta
      std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
        this->universe.getChainsWithCrosslinker(crossLinkerType);
      std::unordered_map<size_t, size_t> atomToStrand;
      atomToStrand.reserve(this->universe.getNrOfAtoms());
      std::unordered_map<size_t, size_t> atomIdxInStrand;
      atomIdxInStrand.reserve(this->universe.getNrOfAtoms());
      size_t springId = 0;
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        pylimer_tools::entities::Molecule chain = crossLinkerChains[i];
        std::vector<pylimer_tools::entities::Atom> atoms =
          crossLinkerChains[i].getAtomsLinedUp(this->crossLinkerType);
        for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
          pylimer_tools::entities::Atom atom = atoms[atomIdx];
          if (atom.getType() != this->crossLinkerType) {
            atomToStrand.emplace(atom.getId(), springId);
            atomIdxInStrand.emplace(atom.getId(), atomIdx);
          }
        }
      }

      // std::cout << "After mapping chains again, memory useage: "
      //           << getCurrentRSS() << ", peak " << getPeakRSS() <<
      //           std::endl;

      // the resulting vectors to fill
      std::unordered_map<long int, long int> intersectionsOfEdges;
      typedef std::tuple<pylimer_tools::entities::LoopIntersectionInfo,
                         std::set<size_t>>
        intersection_loops_tuple;
      std::vector<intersection_loops_tuple> relevantIntersections;

      // as of
      // https://stackoverflow.com/questions/919612/mapping-two-integers-to-one-in-a-unique-and-deterministic-way
      auto hash_integer_pair = [](long int a, long int b) -> long int {
        return a >= b ? a * a + a + b : a + b * b; // where a, b >= 0
      };
      std::cout << "Searching for intersections..." << std::endl;

      for (size_t i = 0; i < loops.size(); ++i) {
        // reserve enough space
        size_t estimateNrOfSliplinks =
          loops.size() * loops.size() *
          this->initialConfig.meanSpringContourLength;
        if (i > 0) {
          estimateNrOfSliplinks =
            std::max(estimateNrOfSliplinks,
                     (loops.size() / i) * relevantIntersections.size());
        }
        relevantIntersections.reserve(estimateNrOfSliplinks);

        const std::array<double, 3> loop_i_min = loopMinCoords[i];
        const std::array<double, 3> loop_i_max = loopMaxCoords[i];

        // TODO: instead of the N^2 loop, might want to try some filter at
        // least also, we ignore all self-entanglements of one loop with
        // itself.
        for (size_t j = i + 1; j < loops.size(); ++j) {
          // first check if the loops have any overlap in 3D, otherwise, we
          // can skip them anyway.

          std::array<double, 3> loop_j_min = loopMinCoords[j];
          std::array<double, 3> loop_j_max = loopMaxCoords[j];

          for (int dir = 0; dir < 3; ++dir) {
            if (
              // "normal" case
              ((loop_i_min[dir] <= loop_i_max[dir] &&
                loop_j_min[dir] <= loop_j_max[dir]) &&
               !(loop_j_min[dir] <= loop_i_max[dir] &&
                 loop_j_max[dir] >= loop_i_min[dir])) ||
              // periodic case in j
              ((loop_i_min[dir] <= loop_i_max[dir] &&
                loop_j_min[dir] > loop_j_max[dir]) &&
               (loop_i_min[dir] >= loop_j_max[dir] &&
                loop_i_max[dir] <= loop_j_min[dir])) ||
              // periodic case in i
              ((loop_i_min[dir] > loop_i_max[dir] &&
                loop_j_min[dir] <= loop_j_max[dir]) &&
               (loop_j_min[dir] >= loop_i_max[dir] &&
                loop_j_max[dir] <= loop_i_min[dir])) ||
              // both periodic -> both pass the boundary -> overlap anyway
              ((loop_i_min[dir] > loop_i_max[dir] &&
                loop_j_min[dir] > loop_j_max[dir]) &&
               false)) {
              //  -> no overlap in dir, skip this comparison
              goto loop_j_checked;
            }
          }

          // then, if we did not goto, we can actually find the entanglements
          {
            std::vector<pylimer_tools::entities::LoopIntersectionInfo>
              intersections = this->universe.findLoopEntanglements(
                loops[i], loops[j], loopEdges[i], loopEdges[j]);
            for (pylimer_tools::entities::LoopIntersectionInfo intersection :
                 intersections) {
              bool keepIntersection = true;

              // check if we want to keep it
              // TODO: check if we want to keep it based on the direction /
              // distance of the intersection

              // check if we want to keep it, based on the edges involved
              long int edgePairHash =
                hash_integer_pair(intersection.edge1, intersection.edge2);
              if (pylimer_tools::utils::map_has_key(intersectionsOfEdges,
                                                    edgePairHash)) {
                keepIntersection = false;
                long int involvedIntersection =
                  intersectionsOfEdges.at(edgePairHash);
                std::get<1>(relevantIntersections[involvedIntersection])
                  .insert(i);
                std::get<1>(relevantIntersections[involvedIntersection])
                  .insert(j);
              }

              if (keepIntersection) {
                std::set<size_t> localSet;
                localSet.insert(i);
                localSet.insert(j);
                intersection_loops_tuple relevantIntersection =
                  std::make_tuple(intersection, localSet);
                relevantIntersections.push_back(relevantIntersection);
                intersectionsOfEdges.insert_or_assign(
                  edgePairHash,
                  static_cast<long int>(relevantIntersections.size()) - 1);
              }
            }
          }

        loop_j_checked:;
        }

        // make space: cleanup the loop i
        loops[i] = std::vector<long int>();
        // std::cout
        //   << "After checking intersections of loop " << i << " ("
        //   << relevantIntersections.size()
        //   << " relevant intersections found in total yet), memory useage: "
        //   << getCurrentRSS() << ", peak " << getPeakRSS() << std::endl;
      }

      // Actually make them into slip-links...
      std::vector<double> slipLinkXs;
      slipLinkXs.reserve(relevantIntersections.size());
      std::vector<double> slipLinkYs;
      slipLinkYs.reserve(relevantIntersections.size());
      std::vector<double> slipLinkZs;
      slipLinkZs.reserve(relevantIntersections.size());
      std::vector<size_t> slipLinkStrandA;
      slipLinkStrandA.reserve(relevantIntersections.size());
      std::vector<size_t> slipLinkStrandB;
      slipLinkStrandB.reserve(relevantIntersections.size());
      std::vector<double> slipLinkStrandAlpha;
      slipLinkStrandAlpha.reserve(relevantIntersections.size());
      std::vector<double> slipLinkStrandBeta;
      slipLinkStrandBeta.reserve(relevantIntersections.size());
      std::vector<std::vector<size_t>> slipLinksLoops;
      slipLinksLoops.reserve(relevantIntersections.size());

      for (intersection_loops_tuple intersectionAndLoops :
           relevantIntersections) {
        pylimer_tools::entities::LoopIntersectionInfo intersection =
          std::get<0>(intersectionAndLoops);
        // TODO: this is yet the most naïve way to add these.
        // ideally, we would also check the back-and-forth, etc.
        slipLinkXs.push_back(intersection.intersectionPoint[0]);
        slipLinkYs.push_back(intersection.intersectionPoint[1]);
        slipLinkZs.push_back(intersection.intersectionPoint[2]);
        // TODO: decide on the atoms to use as a reference
        slipLinkStrandA.push_back(
          this->initialConfig
            .oldAtomIdToSpringIndex[intersection.involvedAtoms[0].getId()]);
        slipLinkStrandB.push_back(
          this->initialConfig
            .oldAtomIdToSpringIndex[intersection.involvedAtoms[3].getId()]);
        slipLinkStrandAlpha.push_back(
          static_cast<double>(
            atomIdxInStrand[intersection.involvedAtoms[0].getId()]) /
          this->initialConfig.meanSpringContourLength);
        slipLinkStrandBeta.push_back(
          static_cast<double>(
            atomIdxInStrand[intersection.involvedAtoms[3].getId()]) /
          this->initialConfig.meanSpringContourLength);
        std::vector<size_t> localLoops(
          std::get<1>(intersectionAndLoops).begin(),
          std::get<1>(intersectionAndLoops).end());
        slipLinksLoops.push_back(localLoops);
      }

      std::cout << "Found " << slipLinksLoops.size() << " intersections."
                << std::endl;
      RUNTIME_EXP_IFN(
        slipLinkStrandA.size() == slipLinkStrandB.size() &&
          slipLinkXs.size() == slipLinkYs.size() &&
          slipLinkZs.size() == slipLinkXs.size(),
        "Expect all slip-link relevant properties to have the same length");
      // with the data assembled, we can actually add them to the structure
      // and stuff
      this->addSlipLinks(slipLinkStrandA,
                         slipLinkStrandB,
                         slipLinkXs,
                         slipLinkYs,
                         slipLinkZs,
                         slipLinkStrandAlpha,
                         slipLinkStrandBeta,
                         reducedLoops,
                         slipLinksLoops,
                         false);
      return slipLinkStrandA.size();
    }

    /**
     * @brief Remove a spring (and all its parts, incl. slip-links) from the
     * structures
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance2::removeSpring(ForceBalance2Network& net,
                                         Eigen::VectorXd& displacements,

                                         const size_t springIdx) const
    {
#ifdef DEBUG_REMOVAL
      std::cout << "Starting to remove spring " << springIdx << std::endl;
#endif
      INVALIDARG_EXP_IFN(springIdx < net.nrOfSprings,
                         "Can only remove springs, not partial springs.");
#ifndef NDEBUG
      Eigen::VectorXd allTotalSpringDistancesBefore =
        this->evaluateSpringLengths(net, displacements, this->is2D);
#endif

      std::vector<size_t> affectedLinks = net.linkIndicesOfSprings[springIdx];
      std::vector<size_t> uniqueAffectedLinks =
        net.linkIndicesOfSprings[springIdx];
      std::sort(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end());
      uniqueAffectedLinks.erase(
        std::unique(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end()),
        uniqueAffectedLinks.end());

      // remove the link to the link, höhö
      for (size_t affectedLinkIdx : uniqueAffectedLinks) {
        RUNTIME_EXP_IFN(
          std::find(net.springIndicesOfLinks[affectedLinkIdx].begin(),
                    net.springIndicesOfLinks[affectedLinkIdx].end(),
                    springIdx) !=
            net.springIndicesOfLinks[affectedLinkIdx].end(),
          "Link must have a connection to the spring, too. Did not find "
          "spring " +
            std::to_string(springIdx) + " in link " +
            std::to_string(affectedLinkIdx) + ", got " +
            pylimer_tools::utils::join(
              net.springIndicesOfLinks[affectedLinkIdx].begin(),
              net.springIndicesOfLinks[affectedLinkIdx].end(),
              std::string(", ")) +
            ".");
        if (net.linkIsEntanglement[affectedLinkIdx]) {
          RUNTIME_EXP_IFN(
            net.springIndicesOfLinks[affectedLinkIdx].size() <= 2,
            "Expect slip-link to be associated with 2 springs only, got " +
              pylimer_tools::utils::join(
                net.springIndicesOfLinks[affectedLinkIdx].begin(),
                net.springIndicesOfLinks[affectedLinkIdx].end(),
                std::string(", ")) +
              ".");
        }

        size_t found =
          std::erase(net.springIndicesOfLinks[affectedLinkIdx], springIdx);

        RUNTIME_EXP_IFN(found > 0,
                        "Expected to find spring " + std::to_string(springIdx) +
                          " in link " + std::to_string(affectedLinkIdx) +
                          " but did not, got " +
                          pylimer_tools::utils::join(
                            net.springIndicesOfLinks[affectedLinkIdx].begin(),
                            net.springIndicesOfLinks[affectedLinkIdx].end(),
                            std::string(", ")) +
                          ".");
        if (net.linkIsEntanglement[affectedLinkIdx]) {
          RUNTIME_EXP_IFN(net.springIndicesOfLinks[affectedLinkIdx].size() <= 1,
                          "Expect slip-link to be associated with 1 springs "
                          "only after removing one, got " +
                            pylimer_tools::utils::join(
                              net.springIndicesOfLinks[affectedLinkIdx].begin(),
                              net.springIndicesOfLinks[affectedLinkIdx].end(),
                              std::string(", ")) +
                            ".");
        }
      }

      std::vector<size_t> affectedPartialSprings =
        net.localToGlobalSpringIndex[springIdx];
      assert(affectedPartialSprings.size() > 0);
      net.nrOfSprings -= 1;
      net.nrOfPartialSprings -= affectedPartialSprings.size();

      // actually spring remove stuff
      net.localToGlobalSpringIndex.erase(net.localToGlobalSpringIndex.begin() +
                                         springIdx);
      net.springToMoleculeIds.erase(net.springToMoleculeIds.begin() +
                                    springIdx);
      net.linkIndicesOfSprings.erase(net.linkIndicesOfSprings.begin() +
                                     springIdx);
      pylimer_tools::utils::removeRow(net.springsContourLength, springIdx);
      pylimer_tools::utils::removeRow(net.springsType, springIdx);
      pylimer_tools::utils::removeRow(net.springIndexA, springIdx);
      pylimer_tools::utils::removeRows(
        net.springCoordinateIndexA, springIdx * 3, 3);
      pylimer_tools::utils::removeRow(net.springIndexB, springIdx);
      pylimer_tools::utils::removeRows(
        net.springCoordinateIndexB, springIdx * 3, 3);
      pylimer_tools::utils::removeRow(net.springIsActive, springIdx);

      // need to remove descending
      std::sort(affectedPartialSprings.begin(),
                affectedPartialSprings.end(),
                std::greater<size_t>());
      for (size_t partialSpringIdx : affectedPartialSprings) {
        pylimer_tools::utils::removeRow(net.springPartIndexA, partialSpringIdx);
        pylimer_tools::utils::removeRows(
          net.springPartCoordinateIndexA, 3 * partialSpringIdx, 3);
        pylimer_tools::utils::removeRow(net.springPartIndexB, partialSpringIdx);
        pylimer_tools::utils::removeRows(
          net.springPartCoordinateIndexB, 3 * partialSpringIdx, 3);
        pylimer_tools::utils::removeRows(
          net.springPartBoxOffset, 3 * partialSpringIdx, 3);
        pylimer_tools::utils::removeRow(net.partialToFullSpringIndex,
                                        partialSpringIdx);
        pylimer_tools::utils::removeRow(net.partialSpringIsPartial,
                                        partialSpringIdx);
        pylimer_tools::utils::removeRow(partialSpringIdx);
      }
      assert(springPartitions.size() == net.nrOfPartialSprings);

      // renumber the remaining stuff
      // first, renumber the springs
      for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
        for (size_t i = 0; i < net.springIndicesOfLinks[linkIdx].size(); ++i) {
          assert(net.springIndicesOfLinks[linkIdx][i] != springIdx);
          if (net.springIndicesOfLinks[linkIdx][i] > springIdx) {
            net.springIndicesOfLinks[linkIdx][i] -= 1;
          }
        }
      }

      // then, renumber the loops
      for (size_t loopIdx = 0; loopIdx < net.loops.size(); ++loopIdx) {
        for (size_t i = 0; i < net.loops[loopIdx].size(); ++i) {
          if (net.loops[loopIdx][i] == springIdx) {
            net.loops[loopIdx].erase(net.loops[loopIdx].begin() + i);
          }
          if (net.loops[loopIdx][i] > springIdx) {
            net.loops[loopIdx][i] -= 1;
          }
        }
      }

      // then, update the partial springs
      // decrease by one if larger then springIdx
      assert((net.partialToFullSpringIndex != springIdx).all());
      net.partialToFullSpringIndex -=
        (net.partialToFullSpringIndex > springIdx).cast<int>();

      assert(net.localToGlobalSpringIndex.size() == net.nrOfSprings);
      for (size_t loopingSpringIdx = 0; loopingSpringIdx < net.nrOfSprings;
           ++loopingSpringIdx) {
        for (size_t i = 0;
             i < net.localToGlobalSpringIndex[loopingSpringIdx].size();
             ++i) {
          for (size_t partialSpringIdx : affectedPartialSprings) {
            if (net.localToGlobalSpringIndex[loopingSpringIdx][i] >
                partialSpringIdx) {
              net.localToGlobalSpringIndex[loopingSpringIdx][i] -= 1;
            }
          }
        }
      }

      //
      // remove the affected slip-links that are now only on one spring
      //
      std::vector<size_t> linksToRemove;
      linksToRemove.reserve(affectedLinks.size() -
                            2); // keep the first and last links
      for (size_t i = 1; i < affectedLinks.size() - 1; ++i) {
        linksToRemove.push_back(affectedLinks[i]);
      }

      // need to remove descending to remove need to renumber these as well
      std::sort(
        linksToRemove.begin(), linksToRemove.end(), std::greater<size_t>());
      linksToRemove.erase(
        std::unique(linksToRemove.begin(), linksToRemove.end()),
        linksToRemove.end());
      if (linksToRemove.size() > 2) {
        assert(linksToRemove[0] > linksToRemove[1]);
      }

      Eigen::ArrayXb springIsAffected =
        Eigen::ArrayXb::Constant(net.nrOfSprings, false);
      for (size_t outermostI = 0; outermostI < linksToRemove.size();
           ++outermostI) {
        size_t slipLinkIdx = linksToRemove[outermostI];
        assert(net.linkIsEntanglement[slipLinkIdx]);
        // first, merge the two other partial springs
        std::vector<size_t> springsOfLink =
          net.springIndicesOfLinks[slipLinkIdx];
        RUNTIME_EXP_IFN(springsOfLink.size() <= 1,
                        "Expected slip-link " + std::to_string(slipLinkIdx) +
                          " to have only 1 remaining spring, got " +
                          std::to_string(springsOfLink.size()) + " due to " +
                          pylimer_tools::utils::join(springsOfLink.begin(),
                                                     springsOfLink.end(),
                                                     std::string(", ")) +
                          ".");
        std::vector<size_t> involvedPartialSprings;
        involvedPartialSprings.reserve(2);
        for (int springInLinkIdx = springsOfLink.size() - 1;
             springInLinkIdx >= 0;
             --springInLinkIdx) {
          for (size_t partialSpringIdx :
               net.localToGlobalSpringIndex[springsOfLink[springInLinkIdx]]) {
            if (net.springPartIndexA[partialSpringIdx] == slipLinkIdx ||
                net.springPartIndexB[partialSpringIdx] == slipLinkIdx) {
              involvedPartialSprings.push_back(partialSpringIdx);
            }
          }
        }
        RUNTIME_EXP_IFN(
          involvedPartialSprings.size() >= springsOfLink.size(),
          "Expected more or equal number of partial springs (" +
            std::to_string(involvedPartialSprings.size()) + "; " +
            pylimer_tools::utils::join(involvedPartialSprings.begin(),
                                       involvedPartialSprings.end(),
                                       std::string(", ")) +
            ") than springs (" + std::to_string(springsOfLink.size()) + "; " +
            pylimer_tools::utils::join(
              springsOfLink.begin(), springsOfLink.end(), std::string(", ")) +
            ").");
        // RUNTIME_EXP_IFN(springsOfLink.size() % 2 == 0, "Expected link to
        // have an even number of components, got " +
        // std::to_string(springsOfLink.size()) + ".");
        if (involvedPartialSprings.size() > 0) {
          RUNTIME_EXP_IFN(
            involvedPartialSprings.size() == 2,
            "Expected only 2 involved partial springs, got: " +
              pylimer_tools::utils::join(involvedPartialSprings.begin(),
                                         involvedPartialSprings.end(),
                                         std::string(", ")) +
              " for springs " +
              pylimer_tools::utils::join(
                springsOfLink.begin(), springsOfLink.end(), std::string(", ")) +
              " when removing spring " + std::to_string(springIdx) + ".");
          assert(involvedPartialSprings.size() == 2);
          size_t partialSpringToKeep =
            std::min(involvedPartialSprings[0], involvedPartialSprings[1]);
          size_t partialSpringToRemove =
            std::max(involvedPartialSprings[0], involvedPartialSprings[1]);
          assert(partialSpringToKeep != partialSpringToRemove);
          assert(net.partialToFullSpringIndex[partialSpringToKeep] ==
                 net.partialToFullSpringIndex[partialSpringToRemove]);
          // actually do the merge
          this->mergePartialSprings(net,
                                    displacements,

                                    partialSpringToRemove,
                                    partialSpringToKeep,
                                    slipLinkIdx);
          // total distance changes -> cannot use for checking the total
          // distance
          springIsAffected[net.partialToFullSpringIndex[partialSpringToKeep]] =
            true;
        }

        assert(net.springIndicesOfLinks[slipLinkIdx].empty());

        // then, actually remove the slip-link
#ifdef DEBUG_REMOVAL
        std::cout << "Removing slip-link " << slipLinkIdx << std::endl;
#endif
        this->removeLink(net, displacements, slipLinkIdx);
      }
#ifdef DEBUG_REMOVAL
      std::cout << "Removed spring " << springIdx << std::endl;
#endif
#ifndef NDEBUG
      Eigen::VectorXd allTotalSpringDistancesAfter =
        this->evaluateSpringLengths(net, displacements, this->is2D);
      assert(allTotalSpringDistancesAfter.size() == net.nrOfSprings);
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        size_t correspondingOldIdx = i >= springIdx ? i + 1 : i;
        if (springIsAffected[i]) {
          // when slip-links are removed, the overall distance must reduce.
          RUNTIME_EXP_IFN(allTotalSpringDistancesBefore[correspondingOldIdx] +
                              1e-9 >=
                            allTotalSpringDistancesAfter[i],
                          "Expected that the total distances stay constant for "
                          "non-changed springs.");
        } else {
          RUNTIME_EXP_IFN(
            APPROX_EQUAL(allTotalSpringDistancesBefore[correspondingOldIdx],
                         allTotalSpringDistancesAfter[i],
                         1e-9),
            "Expected that the total distances stay constant for non-changed "
            "springs.");
        }
      }
#endif
    }

    /**
     * @brief Remove a spring, but also all springs that are connected to it
     * and are connected via entanglement links.
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param springIdx
     */
    void MEHPForceBalance2::removeSpringFollowingEntanglementLinks(
      ForceBalance2Network& net,
      Eigen::VectorXd& displacements,

      const size_t springIdx) const
    {
      std::vector<size_t> springsToRemove = { springIdx };
      std::vector<size_t> linksToRemove = net.linkIndicesOfSprings[springIdx];
      assert(springsToRemove.size() > 0);
      std::sort(
        springsToRemove.begin(), springsToRemove.end(), std::greater<>());
      for (size_t springIdxToDelete : springsToRemove) {
        this->removeSpring(net,
                           displacements,

                           springIdxToDelete);
      }
      std::sort(linksToRemove.begin(), linksToRemove.end(), std::greater<>());
      for (size_t linkIdxToDelete : linksToRemove) {
        assert(net.springIndicesOfLinks[linkIdxToDelete].size() <= 1);
        if (net.springIndicesOfLinks[linkIdxToDelete].size() == 1) {
#ifdef DEBUG_REMOVAL
          std::cout << "Removing additional spring between entanglement links "
                    << net.springIndicesOfLinks[linkIdxToDelete][0]
                    << std::endl;
#endif

          this->removeSpring(
            net, displacements, net.springIndicesOfLinks[linkIdxToDelete][0]);
        }
        this->removeLink(net, displacements, linkIdxToDelete);
      }
    };

    /**
     * @brief break a spring, given its partial spring index
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param partialSpringIdx
     */
    void MEHPForceBalance2::breakPartialSpring(
      ForceBalance2Network& net,
      Eigen::VectorXd& displacements,

      const size_t partialSpringIdx) const
    {
      this->removeSpringFollowingEntanglementLinks(
        net, displacements, net.partialToFullSpringIndex[partialSpringIdx]);
    };

    /**
     * @brief remove a link from the network
     *
     * @param net
     * @param displacements
     * @param linkIdx
     */
    void MEHPForceBalance2::removeLink(ForceBalance2Network& net,
                                       Eigen::VectorXd& displacements,
                                       const size_t linkIdx) const
    {
      INVALIDARG_EXP_IFN(net.springIndicesOfLinks[linkIdx].size() == 0,
                         "The springs have to be removed or re-linked before "
                         "removing the link.");
#ifdef DEBUG_REMOVAL
      std::cout << "Removing link " << linkIdx << std::endl;
#endif

      pylimer_tools::utils::removeRows(net.coordinates, linkIdx * 3, 3);
      pylimer_tools::utils::removeRows(displacements, linkIdx * 3, 3);

      if (!net.linkIsEntanglement[linkIdx]) {
        net.nrOfNodes -= 1;
        pylimer_tools::utils::removeRow(net.oldAtomIds, linkIdx);
        pylimer_tools::utils::removeRow(net.oldAtomTypes, linkIdx);
      }
      net.nrOfLinks -= 1;
      pylimer_tools::utils::removeRow(net.linkIsEntanglement, linkIdx);
      net.springIndicesOfLinks.erase(net.springIndicesOfLinks.begin() +
                                     linkIdx);

      // renumber the remaining links
      for (size_t i = 0; i < net.linkIndicesOfSprings.size(); ++i) {
        for (size_t j = 0; j < net.linkIndicesOfSprings[i].size(); ++j) {
#ifndef NDEBUG
          RUNTIME_EXP_IFN(
            net.linkIndicesOfSprings[i][j] != linkIdx,
            "Expected not to find link to remove " + std::to_string(linkIdx) +
              " in any spring, found in spring " + std::to_string(i) + ", " +
              pylimer_tools::utils::join(net.linkIndicesOfSprings[i].begin(),
                                         net.linkIndicesOfSprings[i].end(),
                                         std::string(", ")) +
              ".");
#endif
          if (net.linkIndicesOfSprings[i][j] > linkIdx) {
            net.linkIndicesOfSprings[i][j] -= 1;
          }
        }
      }
      //
      assert(net.springPartIndexA.size() == net.springPartIndexB.size());
      for (size_t i = 0; i < net.springPartIndexA.size(); ++i) {
#ifndef NDEBUG
        RUNTIME_EXP_IFN(
          net.springPartIndexA[i] != linkIdx,
          "Exected link " + std::to_string(linkIdx) +
            " to not be linked anywhere anymore, found in partial spring " +
            std::to_string(i) + ".");
#endif
        if (net.springPartIndexA[i] > linkIdx) {
          net.springPartIndexA[i] -= 1;
          net.springPartCoordinateIndexA[3 * i] -= 3;
          net.springPartCoordinateIndexA[3 * i + 1] -= 3;
          net.springPartCoordinateIndexA[3 * i + 2] -= 3;
        }
#ifndef NDEBUG
        RUNTIME_EXP_IFN(
          net.springPartIndexB[i] != linkIdx,
          "Exected link " + std::to_string(linkIdx) +
            " to not be linked anywhere anymore, found in partial spring " +
            std::to_string(i) + ".");
#endif
        if (net.springPartIndexB[i] > linkIdx) {
          net.springPartIndexB[i] -= 1;
          net.springPartCoordinateIndexB[3 * i] -= 3;
          net.springPartCoordinateIndexB[3 * i + 1] -= 3;
          net.springPartCoordinateIndexB[3 * i + 2] -= 3;
        }
      }
    }

    /**
     * @brief Combine two springs
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance2::mergePartialSprings(
      ForceBalance2Network& net,
      const Eigen::VectorXd& u,
      const size_t removedPartialSpringIdx,
      const size_t keptPartialSpringIdx,
      const size_t linkToReduce,
      bool skipEigenResize) const
    {
      INVALIDARG_EXP_IFN(net.linkIsEntanglement[linkToReduce],
                         "The link to reduce must be a slip-link");
      INVALIDARG_EXP_IFN(keptPartialSpringIdx != removedPartialSpringIdx,
                         "Cannot merge one spring with the same one.");
      INVALIDARG_EXP_IFN(
        net.partialToFullSpringIndex[keptPartialSpringIdx] ==
          net.partialToFullSpringIndex[removedPartialSpringIdx],
        "The partial springs must be part of the same spring to merge them.");
      INVALIDARG_EXP_IFN(
        (net.springPartIndexA[keptPartialSpringIdx] == linkToReduce &&
         net.springPartIndexB[removedPartialSpringIdx] == linkToReduce) ||
          (net.springPartIndexB[keptPartialSpringIdx] == linkToReduce &&
           net.springPartIndexA[removedPartialSpringIdx] == linkToReduce),
        "Link to reduce must be part of the partial springs "
        "that are to be removed and kept, got " +
          std::to_string(net.springPartIndexA[keptPartialSpringIdx]) + " and " +
          std::to_string(net.springPartIndexB[keptPartialSpringIdx]) +
          " in kept spring, and got " +
          std::to_string(net.springPartIndexA[removedPartialSpringIdx]) +
          " and " +
          std::to_string(net.springPartIndexB[removedPartialSpringIdx]) +
          "in removed spring, instead of " + std::to_string(linkToReduce) +
          ".");
#ifdef DEBUG_REMOVAL
      std::cout << "Merging partial springs " << removedPartialSpringIdx
                << " and " << keptPartialSpringIdx << " around " << linkToReduce
                << std::endl;
#endif

      Eigen::Vector3d distanceBefore =
        this->evaluatePartialSpringDistance(
          net, u, removedPartialSpringIdx, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, keptPartialSpringIdx, this->is2D, false);
      size_t fullSpringIdx = net.partialToFullSpringIndex[keptPartialSpringIdx];
      // start with removal
      net.nrOfPartialSprings -= 1;
      // tell the kept one their new end
      // NOTE: if is possible, if the removedPartialSpring is a primary loop,
      // that this procedure is ambiguous.
      bool removedIsA =
        net.springPartIndexA[removedPartialSpringIdx] == linkToReduce;
      size_t newEnd = removedIsA
                        ? net.springPartIndexB[removedPartialSpringIdx]
                        : net.springPartIndexA[removedPartialSpringIdx];
      if (removedIsA) {
        assert(net.springPartIndexB[keptPartialSpringIdx] == linkToReduce);
        net.springPartIndexB[keptPartialSpringIdx] = newEnd;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexB[3 * keptPartialSpringIdx + dir] =
            3 * newEnd + dir;
        }
      } else {
        assert(net.springPartIndexA[keptPartialSpringIdx] == linkToReduce);
        net.springPartIndexA[keptPartialSpringIdx] = newEnd;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexA[3 * keptPartialSpringIdx + dir] =
            3 * newEnd + dir;
        }
      }
      net.springPartBoxOffset.segment(3 * keptPartialSpringIdx, 3) +=
        net.springPartBoxOffset.segment(3 * removedPartialSpringIdx, 3);
      // remove the spring from the link
      // NOTE: currently, we allow it not to be present,
      // as it might be removed earlier already
      // It is anyway the case, that this function does not necessarily
      // keep the network valid
      int found =
        std::erase(net.springIndicesOfLinks[linkToReduce], fullSpringIdx);
      // TODO: check the origin of this assertion
      // assert(found == 1 || found == 0);
      RUNTIME_EXP_IFN(
        net.localToGlobalSpringIndex[fullSpringIdx].size() ==
          net.linkIndicesOfSprings[fullSpringIdx].size() - 1,
        "Require a global index for each local one, got " +
          std::to_string(net.localToGlobalSpringIndex[fullSpringIdx].size()) +
          " != " +
          std::to_string(net.linkIndicesOfSprings[fullSpringIdx].size() - 1) +
          " for spring " + std::to_string(fullSpringIdx) + ".");
      found = 0;
      // tell the spring of the removed link
      int removed = 0;
      springPartitions[keptPartialSpringIdx] +=
        springPartitions[removedPartialSpringIdx];
      for (int j = net.linkIndicesOfSprings[fullSpringIdx].size() - 1; j >= 0;
           --j) {
        if (net.linkIndicesOfSprings[fullSpringIdx][j] == linkToReduce) {
          if (removed == 0) {
            if ((j > 0 && net.localToGlobalSpringIndex[fullSpringIdx][j - 1] ==
                            removedPartialSpringIdx) ||
                (j < net.localToGlobalSpringIndex[fullSpringIdx].size() &&
                 net.localToGlobalSpringIndex[fullSpringIdx][j] ==
                   removedPartialSpringIdx)) {
              net.linkIndicesOfSprings[fullSpringIdx].erase(
                net.linkIndicesOfSprings[fullSpringIdx].begin() + j);
              removed += 1;
            }
          }
          if (found == 1) {
            // we are dealing with a double -> re-add
            net.springIndicesOfLinks[linkToReduce].push_back(fullSpringIdx);
          } else if (found > 1 && removed > 0) {
            break; // required for certain cases... dangerous, somewhat.
          }
          found += 1;
        }
      }
      this->removeDuplicateListedSpringsFromLink(net, linkToReduce);
      assert(found >= 1 && removed == 1);
      found = 0;
      for (int j = net.localToGlobalSpringIndex[fullSpringIdx].size() - 1;
           j >= 0;
           --j) {
        if (net.localToGlobalSpringIndex[fullSpringIdx][j] ==
            removedPartialSpringIdx) {
          net.localToGlobalSpringIndex[fullSpringIdx].erase(
            net.localToGlobalSpringIndex[fullSpringIdx].begin() + j);
          found += 1;
        }
      }
      assert(found == 1);
      RUNTIME_EXP_IFN(
        net.localToGlobalSpringIndex[fullSpringIdx].size() ==
          net.linkIndicesOfSprings[fullSpringIdx].size() - 1,
        "Require a global index for each local one, got " +
          std::to_string(net.localToGlobalSpringIndex[fullSpringIdx].size()) +
          " != " +
          std::to_string(net.linkIndicesOfSprings[fullSpringIdx].size() - 1) +
          " for spring " + std::to_string(fullSpringIdx) + ".");
      // recompute some values
      net.partialSpringIsPartial[keptPartialSpringIdx] =
        net.linkIndicesOfSprings[fullSpringIdx].size() > 2;
      // actually remove the rows
      pylimer_tools::utils::removeRow(
        net.partialSpringIsPartial, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.partialToFullSpringIndex, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(

        removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.springPartIndexA, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRow(
        net.springPartIndexB, removedPartialSpringIdx, skipEigenResize);
      pylimer_tools::utils::removeRows(net.springPartCoordinateIndexA,
                                       3 * removedPartialSpringIdx,
                                       3,
                                       skipEigenResize);
      pylimer_tools::utils::removeRows(net.springPartCoordinateIndexB,
                                       3 * removedPartialSpringIdx,
                                       3,
                                       skipEigenResize);
      pylimer_tools::utils::removeRows(net.springPartBoxOffset,
                                       3 * removedPartialSpringIdx,
                                       3,
                                       skipEigenResize);
      // renumber stuff
      for (size_t loopSpringIdx = 0;
           loopSpringIdx < net.localToGlobalSpringIndex.size();
           ++loopSpringIdx) {
        for (size_t i = 0;
             i < net.localToGlobalSpringIndex[loopSpringIdx].size();
             ++i) {
          if (net.localToGlobalSpringIndex[loopSpringIdx][i] >
              removedPartialSpringIdx) {
            net.localToGlobalSpringIndex[loopSpringIdx][i] -= 1;
          }
        }
      }

      // validation
      size_t newSpringIdx =
        keptPartialSpringIdx +
        (keptPartialSpringIdx > removedPartialSpringIdx ? -1 : 0);
      Eigen::Vector3d newDistance = this->evaluatePartialSpringDistance(
        net, u, newSpringIdx, this->is2D, false);
      RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(
                        newDistance, distanceBefore, 1e-5),
                      "After merging two partial springs, the overall distance "
                      "is not consistent. Expected distance " +
                        std::to_string(distanceBefore) + ", but got " +
                        std::to_string(newDistance) + " for spring " +
                        std::to_string(newSpringIdx) + ".");
    }

    /**
     * @brief Combine two springs
     *
     * @param net
     * @param springPartitions
     */
    void MEHPForceBalance2::mergeSprings(ForceBalance2Network& net,
                                         const Eigen::VectorXd& u,
                                         const size_t removedSpringIdx,
                                         const size_t keptSpringIdx,
                                         const size_t linkToReduce) const
    {
      INVALIDARG_EXP_IFN(removedSpringIdx < net.nrOfSprings &&
                           keptSpringIdx < net.nrOfSprings,
                         "Only full springs can be merged.");
      INVALIDARG_EXP_IFN(!net.linkIsEntanglement[linkToReduce],
                         "The link to reduce must be a cross-link");
      INVALIDARG_EXP_IFN(keptSpringIdx != removedSpringIdx,
                         "Cannot replace one spring with the same one.");

      size_t fullSpringIdx = net.partialToFullSpringIndex[keptSpringIdx];
      // handle links
      std::vector<size_t> removedSpringsLinks =
        net.linkIndicesOfSprings[removedSpringIdx];
      std::vector<size_t> keptSpringsLinks =
        net.linkIndicesOfSprings[keptSpringIdx];

      size_t removedPartialSpringIdx =
        (removedSpringsLinks[removedSpringsLinks.size() - 1] == linkToReduce)
          ? pylimer_tools::utils::last(
              net.localToGlobalSpringIndex[removedSpringIdx])
          : net.localToGlobalSpringIndex[removedSpringIdx][0];
      size_t remainingPartialSpringIdx =
        (keptSpringsLinks[keptSpringsLinks.size() - 1] == linkToReduce)
          ? pylimer_tools::utils::last(
              net.localToGlobalSpringIndex[keptSpringIdx])
          : net.localToGlobalSpringIndex[keptSpringIdx][0];

      RUNTIME_EXP_IFN(
        net.springPartIndexA[removedPartialSpringIdx] == linkToReduce ||
          net.springPartIndexB[removedPartialSpringIdx] == linkToReduce,
        "Did not detect correct partial springs");
      RUNTIME_EXP_IFN(
        net.springPartIndexA[remainingPartialSpringIdx] == linkToReduce ||
          net.springPartIndexB[remainingPartialSpringIdx] == linkToReduce,
        "Did not detect correct partial springs");

      Eigen::Vector3d distanceBefore = this->evaluatePartialSpringDistance(
        net, u, removedPartialSpringIdx, this->is2D, false);
      Eigen::Vector3d distanceBeforeRemainingSpring =
        this->evaluatePartialSpringDistance(
          net, u, remainingPartialSpringIdx, this->is2D, false);

      net.nrOfSprings -= 1;
      net.nrOfPartialSprings -= 1;

      net.linkIndicesOfSprings[keptSpringIdx].reserve(
        keptSpringsLinks.size() + removedSpringsLinks.size() - 2);
      net.localToGlobalSpringIndex[keptSpringIdx].reserve(
        keptSpringsLinks.size() + removedSpringsLinks.size() - 2);
      RUNTIME_EXP_IFN(net.localToGlobalSpringIndex[keptSpringIdx].size() ==
                        net.linkIndicesOfSprings[keptSpringIdx].size() - 1,
                      "Invalid sizes when merging springs");
      // tell the partial springs their new full spring
      for (size_t partialSpringIndex :
           net.localToGlobalSpringIndex[removedSpringIdx]) {
        net.partialToFullSpringIndex[partialSpringIndex] = keptSpringIdx;
      }
      // std::cout << "Kept spring is "
      //           << pylimer_tools::utils::join(keptSpringsLinks.begin(),
      //                                         keptSpringsLinks.end(),
      //                                         std::string(", "))
      //           << std::endl;
      // actually merge the springs
      if (keptSpringsLinks[keptSpringsLinks.size() - 1] == linkToReduce) {
        // add to end...
        if (removedSpringsLinks[removedSpringsLinks.size() - 1] ==
            linkToReduce) {
          // std::cout << "End end" << std::endl;
          // ...from end
          net.linkIndicesOfSprings[keptSpringIdx][keptSpringsLinks.size() - 1] =
            removedSpringsLinks[removedSpringsLinks.size() - 2];
          for (size_t i = 3; i <= removedSpringsLinks.size(); ++i) {
            net.linkIndicesOfSprings[keptSpringIdx].push_back(
              removedSpringsLinks[removedSpringsLinks.size() - i]);
          }
          for (int i =
                 net.localToGlobalSpringIndex[removedSpringIdx].size() - 2;
               i >= 0;
               --i) {
            net.localToGlobalSpringIndex[keptSpringIdx].push_back(
              net.localToGlobalSpringIndex[removedSpringIdx][i]);
          }
          // invert the direction of these transferred partial springs
          for (size_t partialSpringIdxToInvert :
               net.localToGlobalSpringIndex[removedSpringIdx]) {
            if (partialSpringIdxToInvert == removedPartialSpringIdx) {
              continue;
            }
            std::swap(net.springPartIndexA[partialSpringIdxToInvert],
                      net.springPartIndexB[partialSpringIdxToInvert]);
            for (int dir = 0; dir < 3; ++dir) {
              std::swap(
                net.springPartCoordinateIndexA[3 * partialSpringIdxToInvert +
                                               dir],
                net.springPartCoordinateIndexB[3 * partialSpringIdxToInvert +
                                               dir]);
            }
            net.springPartBoxOffset.segment(3 * partialSpringIdxToInvert, 3) *=
              -1.;
          }
          distanceBefore -= distanceBeforeRemainingSpring;
          distanceBefore *= -1.;
        } else {
          // ...from start
          // std::cout << "End start" << std::endl;
          RUNTIME_EXP_IFN(removedSpringsLinks[0] == linkToReduce,
                          "Things don't make sense anymore.");
          net.linkIndicesOfSprings[keptSpringIdx][keptSpringsLinks.size() - 1] =
            removedSpringsLinks[1];
          for (size_t i = 2; i < removedSpringsLinks.size(); ++i) {
            net.linkIndicesOfSprings[keptSpringIdx].push_back(
              removedSpringsLinks[i]);
          }
          for (size_t i = 1;
               i < net.localToGlobalSpringIndex[removedSpringIdx].size();
               ++i) {
            net.localToGlobalSpringIndex[keptSpringIdx].push_back(
              net.localToGlobalSpringIndex[removedSpringIdx][i]);
          }
          distanceBefore += distanceBeforeRemainingSpring;
        }
      } else {
        RUNTIME_EXP_IFN(keptSpringsLinks[0] == linkToReduce,
                        "How could this be?");
        // add to start...
        if (removedSpringsLinks[removedSpringsLinks.size() - 1] ==
            linkToReduce) {
          // std::cout << "Start end" << std::endl;
          // from end
          net.linkIndicesOfSprings[keptSpringIdx][0] =
            removedSpringsLinks[removedSpringsLinks.size() - 2];
          for (size_t i = 3; i <= removedSpringsLinks.size(); ++i) {
            net.linkIndicesOfSprings[keptSpringIdx].insert(
              net.linkIndicesOfSprings[keptSpringIdx].begin(),
              removedSpringsLinks[removedSpringsLinks.size() - i]);
          }
          for (int i =
                 net.localToGlobalSpringIndex[removedSpringIdx].size() - 2;
               i >= 0;
               --i) {
            net.localToGlobalSpringIndex[keptSpringIdx].insert(
              net.localToGlobalSpringIndex[keptSpringIdx].begin(),
              net.localToGlobalSpringIndex[removedSpringIdx][i]);
          }
          distanceBefore += distanceBeforeRemainingSpring;
        } else {
          // std::cout << "Start start" << std::endl;
          // from start
          RUNTIME_EXP_IFN(removedSpringsLinks[0] == linkToReduce,
                          "No way this exception is ever shown, right?");
          net.linkIndicesOfSprings[keptSpringIdx][0] = removedSpringsLinks[1];
          // have to insert it reverse order
          // happens automatically if we always insert the next the start
          for (size_t i = 2; i < removedSpringsLinks.size(); ++i) {
            net.linkIndicesOfSprings[keptSpringIdx].insert(
              net.linkIndicesOfSprings[keptSpringIdx].begin(),
              removedSpringsLinks[i]);
          }
          // skip the first (removed) partial spring
          for (size_t i = 1;
               i < net.localToGlobalSpringIndex[removedSpringIdx].size();
               ++i) {
            net.localToGlobalSpringIndex[keptSpringIdx].insert(
              net.localToGlobalSpringIndex[keptSpringIdx].begin(),
              net.localToGlobalSpringIndex[removedSpringIdx][i]);
          }

          // invert the direction of these transferred partial springs
          for (size_t partialSpringIdxToInvert :
               net.localToGlobalSpringIndex[removedSpringIdx]) {
            if (partialSpringIdxToInvert == removedPartialSpringIdx) {
              continue;
            }
            std::swap(net.springPartIndexA[partialSpringIdxToInvert],
                      net.springPartIndexB[partialSpringIdxToInvert]);
            for (int dir = 0; dir < 3; ++dir) {
              std::swap(
                net.springPartCoordinateIndexA[3 * partialSpringIdxToInvert +
                                               dir],
                net.springPartCoordinateIndexB[3 * partialSpringIdxToInvert +
                                               dir]);
            }
            net.springPartBoxOffset.segment(3 * partialSpringIdxToInvert, 3) *=
              -1.;
          }
          distanceBefore -= distanceBeforeRemainingSpring;
          distanceBefore *= -1.;
        }
      }
      RUNTIME_EXP_IFN(
        std::find(net.linkIndicesOfSprings[keptSpringIdx].begin(),
                  net.linkIndicesOfSprings[keptSpringIdx].end(),
                  linkToReduce) ==
          net.linkIndicesOfSprings[keptSpringIdx].end(),
        "Link " + std::to_string(linkToReduce) +
          " to reduce should not be in the kept links anymore, found " +
          pylimer_tools::utils::join(
            net.linkIndicesOfSprings[keptSpringIdx].begin(),
            net.linkIndicesOfSprings[keptSpringIdx].end(),
            std::string(", ")) +
          ".");
      assert(net.localToGlobalSpringIndex[keptSpringIdx].size() ==
             net.linkIndicesOfSprings[keptSpringIdx].size() - 1);
      assert(net.linkIndicesOfSprings[keptSpringIdx].size() ==
             keptSpringsLinks.size() + removedSpringsLinks.size() - 2);

      // tell the links of their new spring index
      for (size_t linkOfRemovedSpring : removedSpringsLinks) {
        for (size_t i = 0;
             i < net.springIndicesOfLinks[linkOfRemovedSpring].size();
             ++i) {
          if (net.springIndicesOfLinks[linkOfRemovedSpring][i] ==
              removedSpringIdx) {
            net.springIndicesOfLinks[linkOfRemovedSpring][i] = keptSpringIdx;
          }
        }
        this->removeDuplicateListedSpringsFromLink(
          net, linkOfRemovedSpring, true);
      }

      for (int i = net.springIndicesOfLinks[linkToReduce].size() - 1; i >= 0;
           --i) {
        if (net.springIndicesOfLinks[linkToReduce][i] == removedSpringIdx ||
            net.springIndicesOfLinks[linkToReduce][i] == keptSpringIdx) {
          net.springIndicesOfLinks[linkToReduce].erase(
            net.springIndicesOfLinks[linkToReduce].begin() + i);
        }
      }
      this->removeDuplicateListedSpringsFromLink(net, linkToReduce);

      net.linkIndicesOfSprings.erase(net.linkIndicesOfSprings.begin() +
                                     removedSpringIdx);
      // partial springs
      if (net.partialSpringIsPartial[removedPartialSpringIdx]) {
        net.partialSpringIsPartial[remainingPartialSpringIdx] = true;
      }
      pylimer_tools::utils::removeRow(net.partialSpringIsPartial,
                                      removedPartialSpringIdx);

      bool removedIsA =
        net.springPartIndexA[removedPartialSpringIdx] == linkToReduce;
      size_t otherEndOfRemovedSpring =
        removedIsA ? net.springPartIndexB[removedPartialSpringIdx]
                   : net.springPartIndexA[removedPartialSpringIdx];
      double offsetMultiplier = removedIsA ? -1. : 1.;
      if (net.springPartIndexA[remainingPartialSpringIdx] == linkToReduce) {
        net.springPartIndexA[remainingPartialSpringIdx] =
          otherEndOfRemovedSpring;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexA[3 * remainingPartialSpringIdx + dir] =
            3 * otherEndOfRemovedSpring + dir;
        };
      } else {
        RUNTIME_EXP_IFN(
          net.springPartIndexB[remainingPartialSpringIdx] == linkToReduce, "");
        net.springPartIndexB[remainingPartialSpringIdx] =
          otherEndOfRemovedSpring;
        for (size_t dir = 0; dir < 3; ++dir) {
          net.springPartCoordinateIndexB[3 * remainingPartialSpringIdx + dir] =
            3 * otherEndOfRemovedSpring + dir;
        }
        offsetMultiplier *= -1.;
      }
      net.springPartBoxOffset.segment(3 * remainingPartialSpringIdx, 3) +=
        offsetMultiplier *
        net.springPartBoxOffset.segment(3 * removedPartialSpringIdx, 3);
      pylimer_tools::utils::removeRow(net.springPartIndexA,
                                      removedPartialSpringIdx);
      pylimer_tools::utils::removeRow(net.springPartIndexB,
                                      removedPartialSpringIdx);
      pylimer_tools::utils::removeRows(
        net.springPartCoordinateIndexA, 3 * removedPartialSpringIdx, 3);
      pylimer_tools::utils::removeRows(
        net.springPartCoordinateIndexB, 3 * removedPartialSpringIdx, 3);
      pylimer_tools::utils::removeRows(
        net.springPartBoxOffset, 3 * removedPartialSpringIdx, 3);

      // spring indices & coordinates
      if (net.springIndexA[removedSpringIdx] == linkToReduce) {
        if (net.springIndexA[keptSpringIdx] == linkToReduce) {
          net.springIndexA[keptSpringIdx] = net.springIndexB[removedSpringIdx];
          net.springCoordinateIndexA.segment(3 * keptSpringIdx, 3) =
            net.springCoordinateIndexB.segment(3 * removedSpringIdx, 3);
        } else {
          assert(net.springIndexB[keptSpringIdx] == linkToReduce);
          net.springIndexB[keptSpringIdx] = net.springIndexB[removedSpringIdx];
          net.springCoordinateIndexB.segment(3 * keptSpringIdx, 3) =
            net.springCoordinateIndexB.segment(3 * removedSpringIdx, 3);
        }
      } else {
        assert(net.springIndexB[removedSpringIdx] == linkToReduce);
        if (net.springIndexA[keptSpringIdx] == linkToReduce) {
          net.springIndexA[keptSpringIdx] = net.springIndexA[removedSpringIdx];
          net.springCoordinateIndexA.segment(3 * keptSpringIdx, 3) =
            net.springCoordinateIndexA.segment(3 * removedSpringIdx, 3);
        } else {
          assert(net.springIndexB[keptSpringIdx] == linkToReduce);
          net.springIndexB[keptSpringIdx] = net.springIndexA[removedSpringIdx];
          net.springCoordinateIndexB.segment(3 * keptSpringIdx, 3) =
            net.springCoordinateIndexA.segment(3 * removedSpringIdx, 3);
        }
      }
      pylimer_tools::utils::removeRow(net.springIndexA, removedSpringIdx);
      pylimer_tools::utils::removeRow(net.springIndexB, removedSpringIdx);
      pylimer_tools::utils::removeRows(
        net.springCoordinateIndexA, 3 * removedSpringIdx, 3);
      pylimer_tools::utils::removeRows(
        net.springCoordinateIndexB, 3 * removedSpringIdx, 3);
      pylimer_tools::utils::removeRow(net.springIsActive, removedSpringIdx);
      net.springToMoleculeIds.erase(net.springToMoleculeIds.begin() +
                                    removedSpringIdx);
      net.oldAtomIdToSpringIndex.erase(net.oldAtomIds[linkToReduce]);

      pylimer_tools::utils::removeRow(net.partialToFullSpringIndex,
                                      removedPartialSpringIdx);
      net.localToGlobalSpringIndex.erase(net.localToGlobalSpringIndex.begin() +
                                         removedSpringIdx);
      // renumber the remaining springs
      for (size_t i = 0; i < net.springIndicesOfLinks.size(); ++i) {
        for (size_t j = 0; j < net.springIndicesOfLinks[i].size(); ++j) {
          RUNTIME_EXP_IFN(net.springIndicesOfLinks[i][j] != removedSpringIdx,
                          "Removed spring found in spring indices of link " +
                            std::to_string(i) + ". Must not happen.");
          if (net.springIndicesOfLinks[i][j] > removedSpringIdx) {
            net.springIndicesOfLinks[i][j] -= 1;
          }
        }
      }

      // then, renumber the loops
      for (size_t loopIdx = 0; loopIdx < net.loops.size(); ++loopIdx) {
        for (size_t i = 0; i < net.loops[loopIdx].size(); ++i) {
          if (net.loops[loopIdx][i] == removedSpringIdx) {
            net.loops[loopIdx].erase(net.loops[loopIdx].begin() + i);
          }
          if (net.loops[loopIdx][i] > removedSpringIdx) {
            net.loops[loopIdx][i] -= 1;
          }
        }
      }

      // and the partial springs
      for (size_t i = 0; i < net.partialToFullSpringIndex.size(); ++i) {
        RUNTIME_EXP_IFN(net.partialToFullSpringIndex[i] != removedSpringIdx,
                        "");
        if (net.partialToFullSpringIndex[i] > removedSpringIdx) {
          net.partialToFullSpringIndex[i] -= 1;
        }
      }

      for (size_t i = 0; i < net.localToGlobalSpringIndex.size(); ++i) {
        for (size_t j = 0; j < net.localToGlobalSpringIndex[i].size(); ++j) {
          RUNTIME_EXP_IFN(
            net.localToGlobalSpringIndex[i][j] != removedPartialSpringIdx, "");
          if (net.localToGlobalSpringIndex[i][j] > removedPartialSpringIdx) {
            net.localToGlobalSpringIndex[i][j] -= 1;
          }
        }
      }

      // handle contour lengths
      double contourLengthBefore = net.springsContourLength[keptSpringIdx];
      net.springsContourLength[keptSpringIdx] +=
        net.springsContourLength[removedSpringIdx];
      pylimer_tools::utils::removeRow(net.springsContourLength,
                                      removedSpringIdx);
      pylimer_tools::utils::removeRow(net.springsType, removedSpringIdx);
      RUNTIME_EXP_IFN(net.springsContourLength.size() == net.nrOfSprings, "");
      // and spring partitions
      springPartitions[remainingPartialSpringIdx] +=
        springPartitions[removedPartialSpringIdx];

      pylimer_tools::utils::removeRow(removedPartialSpringIdx);
      RUNTIME_EXP_IFN(springPartitions.size() == net.nrOfPartialSprings, "");
      size_t newKeptSpringIdx = (keptSpringIdx < removedSpringIdx)
                                  ? keptSpringIdx
                                  : (keptSpringIdx - 1);
      // addmittedly, this is possibly dangerous, as it could hide
      // other mistakes
      double newTotalForNormalization =
        springPartitions(net.localToGlobalSpringIndex[newKeptSpringIdx]).sum();
      for (size_t globalPartSpringIndex :
           net.localToGlobalSpringIndex[newKeptSpringIdx]) {
        springPartitions[globalPartSpringIndex] *=
          1. / newTotalForNormalization;
      }

      // std::cout << "Removed springs around " << linkToReduce << " with
      // spring
      // "
      //           << removedSpringIdx << " and partial "
      //           << removedPartialSpringIdx << ", keeping " << keptSpringIdx
      //           << " and " << remainingPartialSpringIdx << std::endl;
      // std::cout << "Spring partitions sum to " << springPartitions.sum()
      //           << " for " << net.nrOfSprings
      //           << " springs, contour length before was " <<
      //           contourLengthBefore
      //           << " and is now " <<
      //           net.springsContourLength[newKeptSpringIdx]
      //           << std::endl;

      // validation
      size_t newPartialSpringIdx =
        remainingPartialSpringIdx +
        (remainingPartialSpringIdx > removedPartialSpringIdx ? -1 : 0);
      Eigen::Vector3d newDistance = this->evaluatePartialSpringDistance(
        net, u, newPartialSpringIdx, this->is2D, false);
      RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(
                        newDistance, distanceBefore, 1e-5),
                      "After merging two springs, the overall distance "
                      "is not consistent. Expected distance " +
                        std::to_string(distanceBefore) + ", but got " +
                        std::to_string(newDistance) + " for spring " +
                        std::to_string(newPartialSpringIdx) + ".");

#ifndef NDEBUG
      // check that the ordering is correct
      for (size_t i = 0;
           i < net.localToGlobalSpringIndex[newKeptSpringIdx].size();
           ++i) {
        size_t partialSpringIdx =
          net.localToGlobalSpringIndex[newKeptSpringIdx][i];
        size_t endA = net.springPartIndexA[partialSpringIdx];
        size_t endB = net.springPartIndexB[partialSpringIdx];
        std::vector<size_t> linkIndices =
          net.linkIndicesOfSprings[newKeptSpringIdx];
        assert(endA == linkIndices[i]);
        assert(endB == linkIndices[i + 1]);
      }
#endif
    }

    /**
     * @brief Remove cross-linkers, springs and associated slip-links with the
     * scheme suggested by Andrei
     *
     * @param net
     * @param displacements
     * @param springPartitions
     * @param tolerance
     * @return size_t
     */
    size_t MEHPForceBalance2::doRemovalAndreisWay(
      ForceBalance2Network& net,
      Eigen::VectorXd& displacements,
      double tolerance) const
    {
      size_t numRemovedTotal = 0;
      size_t numRemovedInIteration = 0;
      do {
        numRemovedInIteration = 0;
        // do removal of f = 1
        // remove all crosslinkers that are 0- or 1-functional
        for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
             --crosslinkIdx) {
          if (net.springIndicesOfLinks[crosslinkIdx].size() == 0 // f = 0
          ) {
            // std::cout << "Removing x-link " << crosslinkIdx << std::endl;
            this->removeLink(net, displacements, crosslinkIdx);
            numRemovedInIteration += 1;
            // this->validateNetwork(net, displacements, springPartitions);
          }

          if ( // or f = 1, NOT primary loop
            net.springIndicesOfLinks[crosslinkIdx].size() == 1 &&
            XOR(
              net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx]
                                                               [0]][0] ==
                crosslinkIdx,
              pylimer_tools::utils::last(
                net.linkIndicesOfSprings[net.springIndicesOfLinks[crosslinkIdx]
                                                                 [0]]) ==
                crosslinkIdx)) {
            // need to first remove the spring
            this->removeSpring(
              net, displacements, net.springIndicesOfLinks[crosslinkIdx][0]);
            // to then remove the link
            this->removeLink(net, displacements, crosslinkIdx);
            numRemovedInIteration += 1;
          }
        }
        numRemovedTotal += numRemovedInIteration;
      } while (numRemovedInIteration > 0);
      // then, replace f = 2
      this->removeTwofunctionalCrosslinks(net, displacements);
      // and remove all springs that are inactive
      size_t numSpringsRemoved = 0;
      for (long int springIdx = net.nrOfSprings - 1; springIdx >= 0;
           --springIdx) {
        size_t a = net.linkIndicesOfSprings[springIdx][0];
        size_t b =
          pylimer_tools::utils::last(net.linkIndicesOfSprings[springIdx]);
        Eigen::Vector3d distance =
          (net.coordinates.segment(3 * a, 3) +
           displacements.segment(3 * a, 3)) -
          (net.coordinates.segment(3 * b, 3) + displacements.segment(3 * b, 3));
        this->box.handlePBC(distance);
        if (this->distanceIsWithinTolerance(
              distance, tolerance, net.springsContourLength[springIdx] ^) &&
            net.linkIndicesOfSprings[springIdx].size() <= 2) {
          // remove
          this->removeSpring(net, displacements, springIdx);
          numSpringsRemoved += 1;
        }
      }

      this->validateNetwork(net, displacements);

      if (numSpringsRemoved > 0) {
        numRemovedTotal +=
          this->doRemovalAndreisWay(net, displacements, ^tolerance);
      }
      return numRemovedTotal;
    };

    /**
     * @brief Replace the two springs traversing a two-functional crosslinkers
     * with a single spring
     *
     * Also handles entanglement beads
     *
     * @param net
     * @param displacements
     * @param springPartitions
     */
    size_t MEHPForceBalance2::removeTwofunctionalCrosslinks(
      ForceBalance2Network& net,
      Eigen::VectorXd& displacements) const
    {
      size_t numRemoved = 0;
      for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
           --crosslinkIdx) {
        if (net.springIndicesOfLinks[crosslinkIdx].size() == 2) {
          std::vector<size_t> springsToMerge =
            net.springIndicesOfLinks[crosslinkIdx];

          assert(springsToMerge.size() == 2);

          // second primary loop check for slip-link entanglements
          // check that it's not a primary loop in any way:
          if (springsToMerge[0] != springsToMerge[1] &&
              (XOR(net.linkIndicesOfSprings[springsToMerge[0]][0] ==
                     crosslinkIdx,
                   pylimer_tools::utils::last(
                     net.linkIndicesOfSprings[springsToMerge[0]]) ==
                     crosslinkIdx)) &&
              (XOR(net.linkIndicesOfSprings[springsToMerge[1]][0] ==
                     crosslinkIdx,
                   pylimer_tools::utils::last(
                     net.linkIndicesOfSprings[springsToMerge[1]]) ==
                     crosslinkIdx))) {
#ifdef DEBUG_REMOVAL
            std::cout << "Merging springs " << springsToMerge[0] << " and "
                      << springsToMerge[1] << " around " << crosslinkIdx
                      << std::endl;
#endif

            // let's remove this
            // TODO: this is inefficient shit, so much data being moved
            this->mergeSprings(net,
                               displacements,

                               springsToMerge[0],
                               springsToMerge[1],
                               crosslinkIdx);

            // this->validateNetwork(net, displacements, springPartitions);
            // std::cout << "Removing link " << crosslinkIdx << std::endl;
            this->removeLink(net, displacements, crosslinkIdx);

            // std::cout << "Removed cross-link " << crosslinkIdx << std::endl;

#ifndef NDEBUG
            this->validateNetwork(net, displacements);
#endif
            numRemoved += 1;
          }
          // else: TODO: decide
        }
      }
#ifndef NDEBUG
      this->validateNetwork(net, displacements);
#endif
      return numRemoved;
    }

    /**
     * @brief Adjust the two spring's box offsets to work best with the
     * specified slip-link
     *
     * @param net the network to adjust
     * @param slipLinkIdx the slip-link around which to adjust the two springs
     * @param spring1 one of the two partial spring idx
     * @param spring2 the partial spring idx of the other spring
     */
    void MEHPForceBalance2::reAlignSlipLinkToImages(
      ForceBalance2Network& net,
      const Eigen::VectorXd& u,
      const size_t slipLinkIdx,
      const size_t partialSpringIdx1,
      const size_t partialSpringIdx2) const
    {
      assert(net.springPartIndexB[partialSpringIdx1] == slipLinkIdx);
      assert(net.springPartIndexA[partialSpringIdx2] == slipLinkIdx);
      assert(net.linkIsEntanglement[slipLinkIdx]);
      assert(net.partialToFullSpringIndex[partialSpringIdx1] ==
             net.partialToFullSpringIndex[partialSpringIdx2]);
      Eigen::Vector3d totalOffset =
        this->getPartialSpringBoxOffset(net, partialSpringIdx1) +
        this->getPartialSpringBoxOffset(net, partialSpringIdx2);
      Eigen::Vector3d totalDistanceBefore =
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx1, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx2, this->is2D, false);

      Eigen::Vector3d sourceCoords =
        net.coordinates.segment(3 * net.springPartIndexA[partialSpringIdx1],
                                3) +
        u.segment(3 * net.springPartIndexA[partialSpringIdx1], 3);
      Eigen::Vector3d targetCoords =
        net.coordinates.segment(3 * net.springPartIndexB[partialSpringIdx2],
                                3) +
        u.segment(3 * net.springPartIndexB[partialSpringIdx2], 3);
      Eigen::Vector3d viaCoords = net.coordinates.segment(3 * slipLinkIdx, 3) +
                                  u.segment(3 * slipLinkIdx, 3);

      // std::cout << net.springPartIndexA[partialSpringIdx1] << ": "
      //           << sourceCoords << " to "
      //           << net.springPartIndexB[partialSpringIdx2] << ": "
      //           << targetCoords << " via " << viaCoords << std::endl;

      double bestOffsetScore = -1.;
      Eigen::Vector3d bestOffset = Eigen::Vector3d::Zero();

      // ugly brute-force method to check all possible combinations (ideally,
      // more or less at least)
      Eigen::Array3i multiplicity1 =
        ((this->box.getOffset(viaCoords - sourceCoords).array().abs() +
          this->getPartialSpringBoxOffset(net, partialSpringIdx1)
            .array()
            .abs()) /
         this->box.getL())
          .rint()
          .cast<int>()
          .abs();
      Eigen::Array3i multiplicity2 =
        ((this->box.getOffset(targetCoords - viaCoords).array() +
          this->getPartialSpringBoxOffset(net, partialSpringIdx2)
            .array()
            .abs()) /
         this->box.getL())
          .rint()
          .cast<int>()
          .abs();
      Eigen::Array3i multiplicity = multiplicity1 + multiplicity2;
      if (this->is2D) {
        multiplicity[2] = 0;
        sourceCoords[2] = 0.;
        targetCoords[2] = 0.;
        viaCoords[2] = 0.;
        totalOffset[2] = 0.;
      }
      for (int mx = std::min(0, -multiplicity[0]);
           mx <= std::max(0, multiplicity[0]);
           ++mx) {
        for (int my = std::min(0, -multiplicity[1]);
             my <= std::max(0, multiplicity[1]);
             ++my) {
          for (int mz = std::min(0, -multiplicity[2]);
               mz <= std::max(0, multiplicity[2]);
               ++mz) {
            Eigen::Vector3d currentOffset;
            currentOffset << mx * net.L[0], my * net.L[1], mz * net.L[2];

            Eigen::Vector3d vec1 = (viaCoords - sourceCoords) + currentOffset;
            Eigen::Vector3d vec2 =
              (targetCoords - viaCoords) + (totalOffset - currentOffset);
            assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
              (vec1 + vec2), totalDistanceBefore));

            double currentScore = vec1.squaredNorm() + vec2.squaredNorm();
            // std::cout << "Score: " << currentScore << " for offset "
            //           << currentOffset << std::endl;
            // std::cout << "vec 1: " << vec1 << std::endl;
            // std::cout << "vec 2: " << vec2 << std::endl;

            if (bestOffsetScore < 0 || bestOffsetScore > currentScore) {
              bestOffsetScore = currentScore;
              bestOffset = currentOffset;
            }
          }
        }
      }

      assert(bestOffsetScore >= 0.);
      net.springPartBoxOffset.segment(3 * partialSpringIdx1, 3) = bestOffset;
      net.springPartBoxOffset.segment(3 * partialSpringIdx2, 3) =
        totalOffset - bestOffset;
      Eigen::Vector3d totalDistanceNow =
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx1, this->is2D, false) +
        this->evaluatePartialSpringDistance(
          net, u, partialSpringIdx2, this->is2D, false);
      assert(pylimer_tools::utils::vector_approx_equal(totalDistanceNow,
                                                       totalDistanceBefore));
    };

    /**
     * @brief Displace all links to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @param u the current displacements, wherein the resulting coordinates
     * shall be stored
     * @param linkIdx the idx of the link to displace
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance2::displaceToMeanPosition(
      const ForceBalance2Network& net,
      Eigen::VectorXd& u,
      const Eigen::ArrayXd& oneOverSpringPartitions) const
    {
      assert(oneOverSpringPartitions.size() == net.nrOfPartialSprings * 3);
      Eigen::ArrayXd objectiveDisplacement =
        Eigen::ArrayXd::Zero(3 * net.nrOfLinks);
      Eigen::ArrayXd partialSpringDistances =
        this
          ->evaluatePartialSpringVectors(
            net, u, this->is2D, this->assumeBoxLargeEnough)
          .array();
      objectiveDisplacement(net.springPartCoordinateIndexA) +=
        (oneOverSpringPartitions * partialSpringDistances);
      objectiveDisplacement(net.springPartCoordinateIndexB) -=
        (oneOverSpringPartitions * partialSpringDistances);

      Eigen::ArrayXd springPartWeightingFactor =
        Eigen::ArrayXd::Zero(net.nrOfLinks * 3);
      Eigen::ArrayXd loopPartialSpringEliminator =
        (net.springPartCoordinateIndexA != net.springPartCoordinateIndexB)
          .cast<double>();

      springPartWeightingFactor(net.springPartCoordinateIndexA) +=
        oneOverSpringPartitions * loopPartialSpringEliminator;
      springPartWeightingFactor(net.springPartCoordinateIndexB) +=
        oneOverSpringPartitions * loopPartialSpringEliminator;
      springPartWeightingFactor = springPartWeightingFactor.unaryExpr(
        [](double v) { return v > 0. ? v : 1.0; });
      Eigen::ArrayXd remainingDisplacement =
        (objectiveDisplacement / springPartWeightingFactor);
#ifndef NDEBUG
      RUNTIME_EXP_IFN(
        pylimer_tools::utils::all_components_finite(remainingDisplacement),
        "Some displacements are not finite");
#endif
      // at this point, we have the ideal displacement if we were to do it
      // just one link at a time.
      // by doing all at once, as here, though, e.g. a pair of links would
      // oscillate back and forth to compensate for that:

      // NOTE: this stays mostly static, could be stored on the network
      Eigen::ArrayXd nSpringsPerLink = Eigen::ArrayXd::Zero(net.nrOfLinks * 3);
      // add a one for every partial spring that's not a primary loop
      nSpringsPerLink(net.springPartCoordinateIndexA) +=
        loopPartialSpringEliminator;
      nSpringsPerLink(net.springPartCoordinateIndexB) +=
        loopPartialSpringEliminator;
      nSpringsPerLink =
        nSpringsPerLink.unaryExpr([](double v) { return v > 0. ? v : 1.0; });
      // make sure there are no infinite back-and-forth
      // and actually displace
      Eigen::ArrayXd backForthDisplacement =
        Eigen::ArrayXd::Zero(net.nrOfLinks * 3);
      backForthDisplacement(net.springPartCoordinateIndexA) +=
        loopPartialSpringEliminator *
        (remainingDisplacement(net.springPartCoordinateIndexB) /
         (nSpringsPerLink(net.springPartCoordinateIndexA) * 2.));
      backForthDisplacement(net.springPartCoordinateIndexB) +=
        loopPartialSpringEliminator *
        (remainingDisplacement(net.springPartCoordinateIndexA) /
         (nSpringsPerLink(net.springPartCoordinateIndexB) * 2.));
#ifndef NDEBUG
      RUNTIME_EXP_IFN(
        pylimer_tools::utils::all_components_finite(backForthDisplacement),
        "Some displacements are not finite");
#endif

      // actually displace
      Eigen::VectorXd finalDisplacement =
        (remainingDisplacement + backForthDisplacement).matrix();
      RUNTIME_EXP_IFN(
        pylimer_tools::utils::all_components_finite(finalDisplacement),
        "Some displacements are not finite");
      // this->box.handlePBC(finalDisplacement);
      u += finalDisplacement;

      double max_disp =
        pylimer_tools::utils::segmentwise_norm_max(finalDisplacement, 3);

      return max_disp;
    }

    /**
     * @brief Displace one link to the mean of all connected neighbours
     *
     * @param net the force balance network
     * @param u the current displacements, wherein the resulting coordinates
     * shall be stored
     * @param linkIdx the idx of the link to displace
     * @return double, the distance (squared norm) displaced
     */
    double MEHPForceBalance2::displaceToMeanPosition(
      const ForceBalance2Network& net,
      Eigen::VectorXd& u,

      const size_t linkIdx,
      const double oneOverSpringPartitionUpperLimit) const
    {
#ifndef NDEBUG
      Eigen::Vector3d forceBefore =
        this->getForceOn(net, u, linkIdx, oneOverSpringPartitionUpperLimit);
#endif
      // Eigen::Vector3d currentDisplacement = u.segment(3 * linkIdx, 3);
      Eigen::Vector3d objectiveDisplacement =
        Eigen::Vector3d::Zero(); // = remainingDisplacement.array();
      double objectiveDisplacementContributors = 0.0;
      bool cautionPrimaryLoop = false;

      std::unordered_set<size_t> partialSpringIndices =
        this->getPartialSpringIndicesOfLink(net, linkIdx);

      for (const size_t globalSpringIndex : partialSpringIndices) {
        assert(net.springPartIndexA[globalSpringIndex] == linkIdx ||
               net.springPartIndexB[globalSpringIndex] == linkIdx);
        if (net.springPartIndexA[globalSpringIndex] == linkIdx &&
            net.springPartIndexB[globalSpringIndex] == linkIdx) {
          // skip primary loops
          continue;
        }
        Eigen::Vector3d partialDistance =
          this->evaluatePartialSpringDistanceFrom(
            net, u, globalSpringIndex, linkIdx);
        const double N = net.springsContourLength
                           [net.partialToFullSpringIndex[globalSpringIndex]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[globalSpringIndex],
          1.0 / (N),
          N,
          oneOverSpringPartitionUpperLimit);

        if (std::isfinite(oneOverContourLengthFraction)) {
          objectiveDisplacement +=
            (partialDistance)*oneOverContourLengthFraction;
          objectiveDisplacementContributors += oneOverContourLengthFraction;
        }
      }
      // take mean for displacement
      // prevent NaN from division by zero
      const double denominator = 1. / (objectiveDisplacementContributors == 0.0
                                         ? 1.0
                                         : objectiveDisplacementContributors);
      u.segment(3 * linkIdx, 3) += objectiveDisplacement * denominator;

#ifndef NDEBUG
      if (!this->assumeBoxLargeEnough) {
        Eigen::Vector3d forceAfter =
          this->getForceOn(net, u, linkIdx, oneOverSpringPartitionUpperLimit);

        // this is only true if we don't have "full" PBC
        assert((pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
          forceAfter, Eigen::Vector3d::Zero(), 0.01)));
        if (!pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
              forceBefore, Eigen::Vector3d::Zero(), 0.01)) {
          assert(forceBefore.squaredNorm() >= forceAfter.squaredNorm());
        }
      }
#endif

      double dist = (objectiveDisplacement * denominator).squaredNorm();
      // if (dist > 0.1) {
      //   std::cout << "Moving " << linkIdx << " for " << dist
      //             << " with displacements " << u.segment(3 * linkIdx, 3)[0]
      //             << ", " << u.segment(3 * linkIdx, 3)[1] << ", "
      //             << u.segment(3 * linkIdx, 3)[2] << std::endl;
      //   std::cout << "For objective displacements " <<
      //   objectiveDisplacement[0]
      //             << ", " << objectiveDisplacement[1] << ", "
      //             << objectiveDisplacement[2] << ", for "
      //             << objectiveDisplacementContributors << "." << std::endl;
      // }
      return dist;
    }

    /**
     * @brief Compute the stress tensor on one cross- or slip-link
     *
     * @param linkIdx
     * @param net
     * @param u
     * @param springPartitions
     * @param debugNrSpringsVisited
     * @param oneOverSpringPartitionUpperLimit
     * @return Eigen::Matrix3d
     */
    Eigen::Matrix3d MEHPForceBalance2::evaluateStressOnLink(
      const size_t linkIdx,
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,

      Eigen::VectorXi& debugNrSpringsVisited,
      const double oneOverSpringPartitionUpperLimit) const
    {
      std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
      Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

      std::unordered_set<size_t> partialSpringIndices =
        this->getPartialSpringIndicesOfLink(net, linkIdx);

      for (const size_t globalSpringIndex : partialSpringIndices) {
        Eigen::Vector3d partialDistance =
          this->evaluatePartialSpringDistanceFrom(
            net, u, globalSpringIndex, linkIdx);
        const double N = net.springsContourLength
                           [net.partialToFullSpringIndex[globalSpringIndex]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[globalSpringIndex],
          1.0 / (N),
          N,
          oneOverSpringPartitionUpperLimit);

        double multiplier = this->kappa * oneOverContourLengthFraction;

        stress += multiplier * partialDistance * partialDistance.transpose();
        debugNrSpringsVisited[globalSpringIndex] += 1;

        // also account for primary loops.
        // they may have non-zero length thanks to assuming the box is not
        // large enough...
        if (net.springPartIndexA[globalSpringIndex] ==
            net.springPartIndexB[globalSpringIndex]) {
          stress +=
            multiplier * (-partialDistance) * (-partialDistance).transpose();

          debugNrSpringsVisited[globalSpringIndex] += 1;
        }
      }

      return stress;
    }

    /**
     * @brief Compute the force acting on one cross- or slip-link
     *
     * @param linkIdx
     * @param net
     * @param u
     * @param springPartitions
     * @param debugNrSpringsVisited
     * @param oneOverSpringPartitionUpperLimit
     * @return Eigen::Vector3d
     */
    Eigen::Vector3d MEHPForceBalance2::evaluateForceOnLink(
      const size_t linkIdx,
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,

      Eigen::VectorXi& debugNrSpringsVisited,
      const double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::Vector3d force = Eigen::Vector3d::Zero();

      std::unordered_set<size_t> partialSpringIndices =
        this->getPartialSpringIndicesOfLink(net, linkIdx);

      for (const size_t globalSpringIndex : partialSpringIndices) {
        // partial spring's force goes both ways -> is zero anyway
        // but, as it would not be included twice in the list,
        // we have to skip them
        if (net.springPartIndexA[globalSpringIndex] == linkIdx &&
            net.springPartIndexB[globalSpringIndex] == linkIdx) {
          if (debugNrSpringsVisited.size() > 0) {
            debugNrSpringsVisited[globalSpringIndex] += 2;
          }
          continue;
        }
        Eigen::Vector3d partialDistance =
          this->evaluatePartialSpringDistanceFrom(
            net, u, globalSpringIndex, linkIdx);
        const double N = net.springsContourLength
                           [net.partialToFullSpringIndex[globalSpringIndex]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[globalSpringIndex],
          1.0 / (N),
          N,
          oneOverSpringPartitionUpperLimit);

        force += this->kappa * oneOverContourLengthFraction * partialDistance;
        if (debugNrSpringsVisited.size() > 0) {
          debugNrSpringsVisited[globalSpringIndex] += 1;
        }
      }

      return force;
    }

    /**
     * @brief Count the number of intra-chain slip-links
     * i.e., slip-links that entangle a strand with itself
     *
     * @return int
     */
    int MEHPForceBalance2::getNumIntraChainSlipLinks() const
    {
      int result = 0;
      for (size_t i = this->initialConfig.nrOfNodes;
           i < this->initialConfig.nrOfLinks;
           ++i) {
        if (this->initialConfig.springIndicesOfLinks[i].size() < 2) {
          result += 1;
        }
        if (this->initialConfig.springIndicesOfLinks[i].size() == 2 &&
            this->initialConfig.springIndicesOfLinks[i][0] ==
              this->initialConfig.springIndicesOfLinks[i][1]) {
          result += 1;
        }
      }

      return result;
    };

    /**
     * @brief Evaluate the vectors between the two ends of all partial springs
     *
     * @param net
     * @param u
     * @param is2D
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance2::evaluatePartialSpringVectors(
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,
      const bool is2D,
      const bool assumeLarge) const
    {
      // first, the distances
      assert(u.size() == net.coordinates.size());

      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd partialDistances =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;

      if (assumeLarge) {
        this->box.handlePBC(partialDistances);
      }

      // reset for 2D systems
      if (is2D) {
        // partialDistances(Eigen::seq(2, net.nrOfPartialSprings, 3)) =
        //   Eigen::VectorXd::Zero(net.nrOfPartialSprings);
        for (size_t i = 2; i < 3 * net.nrOfPartialSprings; i += 3) {
          partialDistances[i] = 0.0;
        }
      }

      return partialDistances;
    }

    /**
     * FORCE BALANCE DATA ACCESS
     */
    /**
     * @brief Convert the current network back into a universe, consisting
     * only of crosslinkers
     */
    pylimer_tools::entities::Universe MEHPForceBalance2::getCrosslinkerVerse()
      const
    {
      // convert nodes & springs back to a universe
      pylimer_tools::entities::Universe xlinkUniverse =
        pylimer_tools::entities::Universe(this->box);
      std::vector<long int> ids;
      std::vector<int> types = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, crossLinkerType);
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      std::vector<int> zeros = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, 0);
      ids.reserve(this->initialConfig.nrOfNodes);
      x.reserve(this->initialConfig.nrOfNodes);
      y.reserve(this->initialConfig.nrOfNodes);
      z.reserve(this->initialConfig.nrOfNodes);
      for (int i = 0; i < this->initialConfig.nrOfNodes; ++i) {
        x.push_back(this->initialConfig.coordinates[3 * i + 0] +
                    this->currentDisplacements[3 * i + 0]);
        y.push_back(this->initialConfig.coordinates[3 * i + 1] +
                    this->currentDisplacements[3 * i + 1]);
        z.push_back(this->initialConfig.coordinates[3 * i + 2] +
                    this->currentDisplacements[3 * i + 2]);
        ids.push_back(this->initialConfig.oldAtomIds[i]);
        // override type, since the types may be different from
        // crossLinkerType if converted with dangling chains
        types[i] = this->initialConfig.oldAtomTypes[i];
      }
      xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
      std::vector<long int> bondFrom;
      std::vector<long int> bondTo;
      bondFrom.reserve(this->initialConfig.nrOfSprings);
      bondTo.reserve(this->initialConfig.nrOfSprings);
      for (int i = 0; i < this->initialConfig.nrOfSprings; ++i) {
        bondFrom.push_back(
          this->initialConfig
            .oldAtomIds[this->initialConfig.linkIndicesOfSprings[i][0]]);
        bondTo.push_back(
          this->initialConfig.oldAtomIds[pylimer_tools::utils::last(
            this->initialConfig.linkIndicesOfSprings[i])]);
      }
      xlinkUniverse.addBonds(
        bondFrom.size(),
        bondFrom,
        bondTo,
        pylimer_tools::utils::initializeWithValue(bondFrom.size(), 1),
        false,
        false); // disable simplify to keep the self-loops etc.
      return xlinkUniverse;
    }

    /**
     * @brief Get the Average Spring Length at the current step
     *
     * @return double
     */
    double MEHPForceBalance2::getAverageSpringLength() const
    {
      Eigen::VectorXd partialSpringVectors = this->evaluatePartialSpringVectors(
        this->initialConfig, this->currentDisplacements);

      Eigen::VectorXd springVectors =
        Eigen::VectorXd::Zero(3 * this->initialConfig.nrOfSprings);
      for (size_t i = 0; i < this->initialConfig.nrOfPartialSprings; ++i) {
        springVectors.segment(this->initialConfig.partialToFullSpringIndex[i],
                              3) += partialSpringVectors.segment(3 * i, 3);
      }

      return pylimer_tools::utils::segmentwise_norm_mean(springVectors, 3);
    }

    /**
     * @brief Get the denominator for a specified partial spring
     *
     * @param net
     * @param springPartitions
     * @param partialSpringIdx
     * @param oneOverSpringPartitionUpperLimit
     * @return double
     */
    double MEHPForceBalance2::getDenominatorOfPartialSpring(
      const ForceBalance2Network& net,

      const size_t partialSpringIdx,
      const double oneOverSpringPartitionUpperLimit) const
    {
      const double N =
        net
          .springsContourLength[net.partialToFullSpringIndex[partialSpringIdx]];

      double denominator = 1. / (N);
      if (oneOverSpringPartitionUpperLimit > 0. ||
          !std::isfinite(denominator)) {
        denominator = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[partialSpringIdx],
          denominator,
          N,
          oneOverSpringPartitionUpperLimit);
      }

      assert(std::isfinite(denominator));
      return denominator;
    }

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @param loopTol
     * @return std::array<std::array<double, 3>, 3>
     */
    Eigen::Matrix3d MEHPForceBalance2::evaluateStressTensorForLinks(
      const std::vector<size_t> linkIndices,
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,

      const double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

      double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

      Eigen::VectorXi debugNrSpringsVisited =
        Eigen::VectorXi::Zero(net.nrOfPartialSprings);

      for (size_t linkIdx : linkIndices) {
        Eigen::Matrix3d force =
          this->evaluateStressOnLink(linkIdx,
                                     net,
                                     u,
                                     debugNrSpringsVisited,
                                     oneOverSpringPartitionUpperLimit);
        /* spring contribution to the overall stress tensor */
        RUNTIME_EXP_IFN(std::isfinite(force.squaredNorm()),
                        "Got non-finite force contribution to stress tensor: " +
                          std::to_string(force.squaredNorm()) + " at link " +
                          std::to_string(linkIdx) + "!");
        stress += force;
      }

      return halfOverVolume * stress;
    };

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3>
    MEHPForceBalance2::evaluateStressTensorLinkBased(
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,
      const double oneOverSpringPartitionUpperLimit,
      const bool xlinksOnly) const
    {
      Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

      double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

      Eigen::VectorXi debugNrSpringsVisited =
        Eigen::VectorXi::Zero(net.nrOfPartialSprings);

      size_t nrOfLinksToInspect = xlinksOnly ? net.nrOfNodes : net.nrOfLinks;
      for (size_t linkIdx = 0; linkIdx < nrOfLinksToInspect; ++linkIdx) {
        Eigen::Matrix3d stressOnLink =
          this->evaluateStressOnLink(linkIdx,
                                     net,
                                     u,
                                     debugNrSpringsVisited,
                                     oneOverSpringPartitionUpperLimit);
        /* spring contribution to the overall stress tensor */
        RUNTIME_EXP_IFN(std::isfinite(stressOnLink.squaredNorm()),
                        "Got non-finite force contribution to stress tensor: " +
                          std::to_string(stressOnLink.squaredNorm()) +
                          " at link " + std::to_string(linkIdx) + "!");
        stress += stressOnLink;
      }

      std::array<std::array<double, 3>, 3> stressA;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          stressA[i][j] = halfOverVolume * stress(i, j);
        }
      }

      if (!xlinksOnly) {
        RUNTIME_EXP_IFN(
          debugNrSpringsVisited.sum() == 2 * net.nrOfPartialSprings,
          "Every spring must be visited twice, got min " +
            std::to_string(debugNrSpringsVisited.minCoeff()) + " and max " +
            std::to_string(debugNrSpringsVisited.maxCoeff()) + ". Sum is " +
            std::to_string(debugNrSpringsVisited.sum()) + " instead of " +
            std::to_string(2 * net.nrOfPartialSprings) + ".");
        RUNTIME_EXP_IFN(
          (debugNrSpringsVisited.array() == 2).all(),
          "Every spring must be visited twice, got min " +
            std::to_string(debugNrSpringsVisited.minCoeff()) + " and max " +
            std::to_string(debugNrSpringsVisited.maxCoeff()) + ".");
      }

      return stressA;
    }

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3>
    MEHPForceBalance2::evaluateStressTensor(
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,

      const double oneOverSpringPartitionUpperLimit) const
    {
      std::array<std::array<double, 3>, 3> stress;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          stress[i][j] = 0.0;
        }
      }

      double oneOverVolume = 1. / (net.L[0] * net.L[1] * net.L[2]);

      Eigen::VectorXd displacedCoords = net.coordinates + u;
      Eigen::VectorXd relevantPartialDistancesA =
        (displacedCoords(net.springPartCoordinateIndexB) -
         displacedCoords(net.springPartCoordinateIndexA)) +
        net.springPartBoxOffset;

      if (this->assumeBoxLargeEnough) {
        this->box.handlePBC(relevantPartialDistancesA);
      }

      if (this->is2D) {
        for (size_t i = 2; i < relevantPartialDistancesA.size(); i += 3) {
          relevantPartialDistancesA[i] = 0.;
        }
      }

      for (size_t partialSpringIdx = 0;
           partialSpringIdx < net.nrOfPartialSprings;
           ++partialSpringIdx) {
        Eigen::Vector3d distance =
          relevantPartialDistancesA.segment(3 * partialSpringIdx, 3);
        const double N = net.springsContourLength
                           [net.partialToFullSpringIndex[partialSpringIdx]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          net.partialSpringIsPartial[partialSpringIdx],
          1.0 / (N),
          N,
          oneOverSpringPartitionUpperLimit);

        /* spring contribution to the overall stress tensor */
        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            double contribution = distance[j] * distance[k] * this->kappa *
                                  oneOverContourLengthFraction;
            RUNTIME_EXP_IFN(
              std::isfinite(contribution),
              "Got non-finite contribution to stress tensor: " +
                std::to_string(contribution) + " at coordinates " +
                std::to_string(k) + ", " + std::to_string(j) +
                " for partial spring " + std::to_string(partialSpringIdx) +
                " from distances " + std::to_string(distance[j]) + ", " +
                std::to_string(distance[k]) + " and denominator " +
                std::to_string(oneOverContourLengthFraction) + ".");
            // if (std::isfinite(denominator) && std::isfinite(contribution))
            // {
            stress[j][k] += contribution;
            // }
          }
        }
      }

      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          stress[i][j] *= oneOverVolume;
          RUNTIME_EXP_IFN(std::isfinite(stress[i][j]),
                          "Got non-finite stress tensor component: " +
                            std::to_string(stress[i][j]) + " at coordinates " +
                            std::to_string(i) + ", " + std::to_string(j) +
                            " from denominator " +
                            std::to_string(oneOverVolume) + ".");
        }
      }

      return stress;
    }

    Eigen::Matrix3d MEHPForceBalance2::getStressTensor(
      const double oneOverSpringPartitionUpperLimit) const
    {
      std::array<std::array<double, 3>, 3> res =
        this->evaluateStressTensor(this->initialConfig,
                                   this->currentDisplacements,
                                   oneOverSpringPartitionUpperLimit);

      // convert the array to an Eigen matrix
      Eigen::Matrix3d convertedRes = Eigen::Matrix3d::Zero();
      for (size_t i = 0; i < 3; ++i) {
        convertedRes.row(i) = Eigen::Vector3d::Map(res[i].data(), 3);
      }
      return convertedRes;
    }

    Eigen::Matrix3d MEHPForceBalance2::getStressTensorLinkBased(
      const double oneOverSpringPartitionUpperLimit,
      const bool xlinksOnly) const
    {
      std::array<std::array<double, 3>, 3> res =
        this->evaluateStressTensorLinkBased(this->initialConfig,
                                            this->currentDisplacements,
                                            oneOverSpringPartitionUpperLimit,
                                            xlinksOnly);
      Eigen::Matrix3d convertedRes = Eigen::Matrix3d::Zero();
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          convertedRes(i, j) = res[i][j];
        }
      }
      return convertedRes;
    }

    /**
     * @brief Get the Effective Functionality Of each node
     *
     * Returns the number of active springs connected to each atom, atomId
     * used as index
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return std::unordered_map<long int, int>
     */
    std::unordered_map<long int, int>
    MEHPForceBalance2::getEffectiveFunctionalityOfAtoms(double tolerance) const
    {
      std::unordered_map<long int, int> results;
      results.reserve(this->initialConfig.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        this->getNrOfActiveSpringsConnected(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
        results.emplace(this->initialConfig.oldAtomIds[i],
                        nrOfActiveSpringsConnected[i]);
      }
      return results;
    }

    /**
     * @brief Get the indices of active Nodes
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return std::vector<long int> the atom ids
     */
    std::vector<long int> MEHPForceBalance2::getIndicesOfActiveNodes(
      const ForceBalance2Network* net,
      const Eigen::VectorXd& u,

      double tolerance) const
    {
      std::vector<long int> results;
      results.reserve(net->nrOfNodes);

      // find all active springs
      Eigen::ArrayXb springIsActive =
        this->findActiveSprings(net, u, tolerance);

      for (size_t i = 0; i < net->nrOfNodes; i++) {
        std::vector<size_t> springIndices = net->springIndicesOfLinks[i];
        for (const size_t springIndex : springIndices) {
          if (springIsActive[springIndex]) {
            results.push_back(i);
            break;
          }
        }
      }

      return results;
    };

    /**
     * @brief Get the atom ids of the active cross-links (not entanglement
     * beads/links)
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return std::vector<long int> the atom ids
     */
    std::vector<long int> MEHPForceBalance2::getIdsOfActiveNodes(
      double tolerance) const
    {
      std::vector<long int> results;
      // find all active springs
      std::vector<long int> activeNodes = this->getIndicesOfActiveNodes(
        &this->initialConfig, this->currentDisplacements, tolerance);

      results.reserve(activeNodes.size());

      for (long int nodeIdx : activeNodes) {
        results.push_back(this->initialConfig.oldAtomIds[nodeIdx]);
      }

      return results;
    }

    /**
     * @brief Get the Nr Of Active Springs connected to each node
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return Eigen::VectorXi
     */
    Eigen::VectorXi MEHPForceBalance2::getNrOfActiveSpringsConnected(
      double tolerance) const
    {
      Eigen::VectorXi nrOfActiveSpringsConnected =
        Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
      Eigen::ArrayXb springIsActive = this->findActiveSprings(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfSprings; i++) {
        if (springIsActive[i]) {
          /* active spring */
          int a = this->initialConfig.linkIndicesOfSprings[i][0];
          int b = pylimer_tools::utils::last(
            this->initialConfig.linkIndicesOfSprings[i]);
          ++(nrOfActiveSpringsConnected[a]);
          ++(nrOfActiveSpringsConnected[b]);
        }
      }
      return nrOfActiveSpringsConnected;
    }

    /**
     * @brief Get the Nr Of Active Springs connected to each node
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return Eigen::VectorXi
     */
    Eigen::VectorXi MEHPForceBalance2::getNrOfActivePartialSpringsConnected(
      double tolerance) const
    {
      Eigen::VectorXi nrOfActivePartialSpringsConnected =
        Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
      Eigen::ArrayXb partialSpringIsActive =
        this->findActivePartialSprings(tolerance);
      // translate this to the nodes
      for (size_t i = 0; i < this->initialConfig.nrOfPartialSprings; ++i) {
        if (partialSpringIsActive[i]) {
          /* active spring */
          // size_t a =
          //   this->initialConfig
          //     .springIndexA[this->initialConfig.partialToFullSpringIndex[i]];
          // size_t b =
          //   this->initialConfig
          //     .springIndexB[this->initialConfig.partialToFullSpringIndex[i]];
          size_t a = this->initialConfig.springPartIndexA[i];
          size_t b = this->initialConfig.springPartIndexB[i];
          if (!this->initialConfig.linkIsEntanglement[a]) {
            ++(nrOfActivePartialSpringsConnected[a]);
          }

          if (!this->initialConfig.linkIsEntanglement[b]) {
            ++(nrOfActivePartialSpringsConnected[b]);
          }
        }
      }
      return nrOfActivePartialSpringsConnected;
    }

    /**
     * @brief Get the Gamma Factor at the current step
     *
     * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
     * computed as phantom = N<b^2>.
     * @param nrOfChains the nr of chains to average over (can be different
     * from the nr of springs thanks to omitted free chains or primary loops)
     * @return double
     */
    double MEHPForceBalance2::getGammaFactor(
      double b02,
      int nrOfChains,
      double oneOverSpringPartitionUpperLimit) const
    {
      if (b02 < 0) {
        b02 = this->defaultBondLength * this->defaultBondLength;
      }

      Eigen::VectorXd gammaFactors = this->getGammaFactors(b02);

      if (nrOfChains < 1) {
        return gammaFactors.mean();
      } else {
        return gammaFactors.sum() / static_cast<double>(nrOfChains);
      }
    }

    /**
     * @brief Get the per-(partial)-spring gamma factors
     *
     * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
     * computed as phantom = N<b^2>.
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance2::getGammaFactors(
      double b02,
      double oneOverSpringPartitionUpperLimit) const
    {
      Eigen::VectorXd springVectors = this->evaluatePartialSpringVectors(
        this->initialConfig, this->currentDisplacements);

      Eigen::VectorXd gammaFactors(springVectors.size() / 3);
      const double commonDenominator = 1. / b02;
      for (size_t i = 0; i < springVectors.size() / 3; ++i) {
        double N = this->initialConfig.springsContourLength
                     [this->initialConfig.partialToFullSpringIndex[i]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          this->initialConfig.partialSpringIsPartial[i],
          1.0 / (N),
          N,
          oneOverSpringPartitionUpperLimit);
        gammaFactors[i] = springVectors.segment(3 * i, 3).squaredNorm() *
                          commonDenominator * oneOverContourLengthFraction;
        RUNTIME_EXP_IFN(
          std::isfinite(gammaFactors[i]),
          "Non-finite gamma factor for partial spring " + std::to_string(i) +
            ", computed from N = " + std::to_string(N) + ", b02 = " +
            std::to_string(b02) + ", oneOverSpringPartitionUpperLimit = " +
            std::to_string(oneOverSpringPartitionUpperLimit) +
            " and squared distance = " +
            std::to_string(springVectors.segment(3 * i, 3).squaredNorm()) +
            ".");
      }
      return gammaFactors;
    }

    /**
     * @brief Get the per-(partial)-spring gamma factors
     *
     * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
     * computed as phantom = N<b^2>.
     * @param dir the direction (0=x, 1=y, 2=z)
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceBalance2::getGammaFactorsInDir(
      double b02,
      int dir,
      double oneOverSpringPartitionUpperLimit) const
    {
      INVALIDARG_EXP_IFN(dir >= 0 && dir <= 2, "Invalid direction.");
      Eigen::VectorXd springVectors = this->evaluatePartialSpringVectors(
        this->initialConfig, this->currentDisplacements);

      Eigen::VectorXd gammaFactors(springVectors.size() / 3);
      const double commonDenominator = 1. / b02;
      for (size_t i = 0; i < springVectors.size() / 3; ++i) {
        double N = this->initialConfig.springsContourLength
                     [this->initialConfig.partialToFullSpringIndex[i]];
        double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
          this->initialConfig.partialSpringIsPartial[i],
          1.0 / (N),
          N,
          oneOverSpringPartitionUpperLimit);
        gammaFactors[i] = SQUARE(springVectors[3 * i + dir]) *
                          commonDenominator * oneOverContourLengthFraction;
        RUNTIME_EXP_IFN(std::isfinite(gammaFactors[i]),
                        "Non-finite gamma factor for partial spring " +
                          std::to_string(i) +
                          ", computed from N = " + std::to_string(N) +
                          ", oneOverSpringPartitionUpperLimit = " +
                          std::to_string(oneOverSpringPartitionUpperLimit) +
                          " and squared distance = " +
                          std::to_string(SQUARE(springVectors[3 * i + dir])) +
                          " in dir " + std::to_string(dir) + ".");
      }
      return gammaFactors;
    }

    /**
     * @brief Get the Weighted Partial Spring Length
     *
     * @return double
     */
    double MEHPForceBalance2::getWeightedPartialSpringLength(
      const ForceBalance2Network& net,
      const Eigen::VectorXd& u,

      size_t partialSpringIdx,
      double oneOverSpringPartitionUpperLimit) const
    {
      double N =
        net
          .springsContourLength[net.partialToFullSpringIndex[partialSpringIdx]];
      double oneOverContourLengthFraction = CLAMP_ONE_OVER_SPRINGPARTITION(
        net.partialSpringIsPartial[partialSpringIdx],
        1.0 / (N),
        N,
        oneOverSpringPartitionUpperLimit);
      return this->evaluatePartialSpringDistance(net, u, partialSpringIdx)
               .norm() *
             oneOverContourLengthFraction;
    }

    /**
     * @brief Convert the universe to a network
     *
     * @param net the target network
     * @param crossLinkerType the atom type of the crossLinker
     * @return true
     * @return false
     */
    bool MEHPForceBalance2::ConvertNetwork(ForceBalance2Network& net,
                                           const int crossLinkerType,
                                           bool remove2functionalCrosslinkers,
                                           bool removeDanglingChains)
    {
      if (remove2functionalCrosslinkers) {
        for (pylimer_tools::entities::Atom xlinker :
             this->universe.getAtomsOfType(crossLinkerType)) {
          // change type of cross-linkers with a degree <= 2 to "normal",
          // non-cross-link beads
          size_t vertexId = this->universe.getIdxByAtomId(xlinker.getId());
          if (this->universe.computeFunctionalityForVertex(vertexId) <= 2) {
            this->universe.setPropertyValue(
              vertexId, "type", crossLinkerType - 1);
          }
        }
      }

      std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
        this->universe.getChainsWithCrosslinker(crossLinkerType);

      // need to include all but dangling and free chains in order to
      // model entanglement
      size_t nrOfSprings = 0;
      std::vector<bool> useChain =
        pylimer_tools::utils::initializeWithValue<bool>(
          crossLinkerChains.size(), false);
      std::vector<long int> vertexIdToLinkIdx =
        pylimer_tools::utils::initializeWithValue<long int>(
          this->universe.getNrOfAtoms(), -1);
      size_t currentLinkIdx = 0;
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        RUNTIME_EXP_IFN(crossLinkerChains[i].getType() !=
                          pylimer_tools::entities::MoleculeType::UNDEFINED,
                        "Cross-linker chain's chain type could not be "
                        "detected. Cannot work like that.");
        if (crossLinkerChains[i].getType() ==
            pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
          assert(
            crossLinkerChains[i].getChainEnds(crossLinkerType, true).size() ==
            2);
          useChain[i] = true;
          nrOfSprings += 1;
        } else if (crossLinkerChains[i].getType() ==
                   pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
          // when omitting f=2 cross-links, it's possible that we end up with
          // "free" primary loops – let's not use those
          if (crossLinkerChains[i].getAtomsOfType(crossLinkerType).size() > 0) {
            useChain[i] = true;
            nrOfSprings += 1;
          }
        } else if (!removeDanglingChains &&
                   crossLinkerChains[i].getType() ==
                     pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) {
          std::vector<pylimer_tools::entities::Atom> endAtoms =
            crossLinkerChains[i].getAtomsOfDegree(1);
          RUNTIME_EXP_IFN(endAtoms.size() == 2,
                          "Expected a dangling chain to have two ends, got " +
                            std::to_string(endAtoms.size()) + ".");
          // cannot assert this without entanglement types being set
          // assert(XOR((endAtoms[0].getType() == crossLinkerType),
          //            endAtoms[1].getType() == crossLinkerType));

          useChain[i] = true;
          nrOfSprings += 1;
        }

        if (useChain[i]) {
          std::vector<pylimer_tools::entities::Atom> endAtoms =
            crossLinkerChains[i].getChainEnds(crossLinkerType);
          for (const pylimer_tools::entities::Atom& endAtom : endAtoms) {
            size_t atomVertexIdx =
              this->universe.getIdxByAtomId(endAtom.getId());
            if (vertexIdToLinkIdx[atomVertexIdx] == -1) {
              vertexIdToLinkIdx[atomVertexIdx] = currentLinkIdx;
              currentLinkIdx += 1;
            }
          }
        }
      }

      size_t nrOfXlinks = currentLinkIdx;

      // crossLinkerUniverse.simplify();
      pylimer_tools::entities::Box box = this->box;
      net.L[0] = box.getLx();
      net.L[1] = box.getLy();
      net.L[2] = box.getLz();
      net.boxHalfs[0] = 0.5 * net.L[0];
      net.boxHalfs[1] = 0.5 * net.L[1];
      net.boxHalfs[2] = 0.5 * net.L[2];
      net.nrOfNodes = nrOfXlinks;
      net.nrOfLinks = nrOfXlinks;
      net.nrOfSprings = nrOfSprings;
      net.nrOfPartialSprings = nrOfSprings;
      net.coordinates = Eigen::VectorXd::Zero(3 * net.nrOfLinks);
      net.oldAtomIds = Eigen::ArrayXi::Zero(net.nrOfLinks);
      net.oldAtomTypes = Eigen::ArrayXi::Zero(net.nrOfLinks);
      net.linkIsEntanglement = Eigen::ArrayXb::Constant(net.nrOfLinks, false);
      net.springIndicesOfLinks.reserve(net.nrOfLinks);
      for (size_t i = 0; i < net.nrOfLinks; ++i) {
        net.springIndicesOfLinks.push_back(std::vector<size_t>());
      }
      net.linkIndicesOfSprings.reserve(net.nrOfSprings);
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        net.linkIndicesOfSprings.push_back(std::vector<size_t>());
      }
      net.partialToFullSpringIndex = Eigen::ArrayXi(net.nrOfPartialSprings);
      net.springPartBoxOffset = Eigen::VectorXd::Zero(3 * net.nrOfSprings);
      net.springsContourLength = Eigen::VectorXd::Zero(net.nrOfSprings);
      net.springsType = Eigen::ArrayXi::Zero(net.nrOfSprings);
      net.partialSpringIsPartial =
        Eigen::ArrayXb::Constant(net.nrOfSprings, false);

      // convert (cross-linker-)beads
      std::map<int, int> atomIdToNode;
      for (size_t i = 0; i < vertexIdToLinkIdx.size(); ++i) {
        if (vertexIdToLinkIdx[i] != -1) {
          size_t linkIdx = vertexIdToLinkIdx[i];
          pylimer_tools::entities::Atom atom =
            this->universe.getAtomByVertexIdx(i);
          atomIdToNode[atom.getId()] = linkIdx;
          net.oldAtomIds[linkIdx] = atom.getId();
          net.oldAtomTypes[linkIdx] = atom.getType();
          Eigen::Vector3d coords = atom.getCoordinates();
          this->box.handlePBC(coords);
          net.coordinates.segment(3 * linkIdx, 3) = coords;
        }
      }

      // convert springs
      size_t spring_idx = 0;
      Eigen::Vector3d expectedDistance = Eigen::Vector3d::Zero();
      // net.connectivityToSpringIndex.reserve(nrOfSprings);
      for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
        if (!useChain[i]) {
          continue;
        }
        std::vector<pylimer_tools::entities::Atom> xlinkersOfChain =
          crossLinkerChains[i].getAtomsOfType(crossLinkerType);
        std::vector<pylimer_tools::entities::Atom> chainEnds =
          crossLinkerChains[i].getChainEnds(crossLinkerType, true);
        RUNTIME_EXP_IFN(
          chainEnds.size() == 2,
          "Expected two chain ends when converting structure. Got " +
            std::to_string(chainEnds.size()) + ".");
        long int atomIdFrom = chainEnds[0].getId();
        long int atomIdTo = chainEnds[1].getId();
        bool addChain = false;
        if (crossLinkerChains[i].getType() ==
            pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
          addChain = true;

          // spring contour length = nr of bonds between two cross-linkers
          net.springsContourLength[spring_idx] =
            crossLinkerChains[i].getNrOfAtoms() - 1;
        } else if (crossLinkerChains[i].getType() ==
                   pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
          addChain = true;

          net.springsContourLength[spring_idx] =
            crossLinkerChains[i].getNrOfAtoms();
          if (xlinkersOfChain.size() == 2) {
            net.springsContourLength[spring_idx] =
              crossLinkerChains[i].getNrOfAtoms() - 1;
          }
        } else if (crossLinkerChains[i].getType() ==
                     pylimer_tools::entities::MoleculeType::DANGLING_CHAIN &&
                   !removeDanglingChains) {
          net.springsContourLength[spring_idx] =
            crossLinkerChains[i].getNrOfAtoms() - 1;
          addChain = true;
        }
        auto bondTypes = crossLinkerChains[i].getBonds()["bond_type"];
        net.springsType[spring_idx] = MEAN(bondTypes);
        assert(addChain);

        if (addChain) {
          long int nodeIdxFrom = atomIdToNode.at(atomIdFrom);
          long int nodeIdxTo = atomIdToNode.at(atomIdTo);
          // if (nodeIdxFrom > nodeIdxTo) {
          //   std::swap(nodeIdxFrom, nodeIdxTo);
          //   std::swap(atomIdFrom, atomIdTo);
          // }

          std::vector<pylimer_tools::entities::Atom> allChainAtoms =
            crossLinkerChains[i].getAtoms();

          pylimer_tools::utils::addIfNotContained(
            net.springIndicesOfLinks[nodeIdxFrom], spring_idx);
          if (nodeIdxFrom != nodeIdxTo) {
            pylimer_tools::utils::addIfNotContained(
              net.springIndicesOfLinks[nodeIdxTo], spring_idx);
          }

          net.linkIndicesOfSprings[spring_idx].push_back(nodeIdxFrom);
          net.linkIndicesOfSprings[spring_idx].push_back(nodeIdxTo);

          std::vector<size_t> zeroMap;
          zeroMap.push_back(spring_idx);
          net.localToGlobalSpringIndex.push_back(zeroMap);
          net.partialToFullSpringIndex[spring_idx] = (spring_idx);

          expectedDistance = crossLinkerChains[i].getOverallBondSumFromTo(
            atomIdFrom, atomIdTo, crossLinkerType, true);
          Eigen::Vector3d actualDistance =
            net.coordinates.segment(3 * nodeIdxTo, 3) -
            net.coordinates.segment(3 * nodeIdxFrom, 3);
          net.springPartBoxOffset.segment(3 * spring_idx, 3) =
            expectedDistance - actualDistance;
          assert(this->box.isValidOffset(expectedDistance - actualDistance));

          spring_idx += 1;
        }
      }

      return spring_idx == net.nrOfSprings;
    };

    bool MEHPForceBalance2::validateNetwork(const ForceBalance2Network& net,
                                            const Eigen::VectorXd& u) const
    {
      // std::cout << "Validating network..." << std::endl;
      /**
       * First, test dimensions
       */
      RUNTIME_EXP_IFN(!std::isinf(net.L[0]) && !std::isnan(net.L[0]),
                      "Box direction x must be scalar");
      RUNTIME_EXP_IFN(!std::isinf(net.L[1]) && !std::isnan(net.L[1]),
                      "Box direction y must be scalar");
      RUNTIME_EXP_IFN(!std::isinf(net.L[2]) && !std::isnan(net.L[2]),
                      "Box direction z must be scalar");
      RUNTIME_EXP_IFN(net.coordinates.size() == net.nrOfLinks * 3,
                      "Invalid size of coordinates");
      RUNTIME_EXP_IFN(u.size() == net.nrOfLinks * 3,
                      "Invalid size of displacements");
      RUNTIME_EXP_IFN(u.size() == net.coordinates.size(),
                      "Invalid size of displacements or coordinates");
      RUNTIME_EXP_IFN(net.localToGlobalSpringIndex.size() == net.nrOfSprings,
                      "Invalid size of connectivity map, got " +
                        std::to_string(net.localToGlobalSpringIndex.size()) +
                        " for " + std::to_string(net.nrOfSprings) +
                        " springs.");
      RUNTIME_EXP_IFN(net.springsContourLength.size() == net.nrOfSprings,
                      "Invalid size of contour lengths, got " +
                        std::to_string(net.springsContourLength.size()) +
                        " for " + std::to_string(net.nrOfSprings) +
                        " springs.");
      RUNTIME_EXP_IFN(net.springsType.size() == net.nrOfSprings,
                      "Invalid size of springs types, got " +
                        std::to_string(net.springsType.size()) + " for " +
                        std::to_string(net.nrOfSprings) + " springs.");
      RUNTIME_EXP_IFN(net.springIndicesOfLinks.size() == net.nrOfLinks,
                      "Invalid size of spring indices of links, got " +
                        std::to_string(net.linkIndicesOfSprings.size()) +
                        " for " + std::to_string(net.nrOfSprings) +
                        " springs.");
      RUNTIME_EXP_IFN(net.linkIndicesOfSprings.size() == net.nrOfSprings,
                      "Invalid size of link indices of springs, got " +
                        std::to_string(net.linkIndicesOfSprings.size()) +
                        " for " + std::to_string(net.nrOfSprings) +
                        " springs.");
      RUNTIME_EXP_IFN(net.linkIsEntanglement.size() == net.nrOfLinks,
                      "Invalid size of link is sliplink");
      RUNTIME_EXP_IFN(
        net.linkIsEntanglement.count() == (net.nrOfLinks - net.nrOfNodes),
        "Nr of nodes plus nr of slp-links should give the total nr of links");
      RUNTIME_EXP_IFN(net.oldAtomIds.size() == net.nrOfNodes,
                      "Invalid size of old atom ids");
      RUNTIME_EXP_IFN(net.oldAtomTypes.size() == net.nrOfNodes,
                      "Invalid size of old atom types");
      RUNTIME_EXP_IFN(net.springPartCoordinateIndexA.size() ==
                        net.nrOfPartialSprings * 3,
                      "Invalid size of springPartCoordinateIndexA");
      RUNTIME_EXP_IFN(net.springPartCoordinateIndexB.size() ==
                        net.nrOfPartialSprings * 3,
                      "Invalid size of springPartCoordinateIndexB");
      RUNTIME_EXP_IFN(net.springPartIndexA.size() == net.nrOfPartialSprings,
                      "Invalid size of springPartIndexA");
      RUNTIME_EXP_IFN(net.springPartBoxOffset.size() ==
                        net.nrOfPartialSprings * 3,
                      "Invalid size of springPartBoxOffset");
      RUNTIME_EXP_IFN(net.springPartIndexB.size() == net.nrOfPartialSprings,
                      "Invalid size of springPartIndexB");
      RUNTIME_EXP_IFN(net.partialSpringIsPartial.size() ==
                        net.nrOfPartialSprings,
                      "Invalid size of partialSpringIsPartial");
      RUNTIME_EXP_IFN(
        net.partialToFullSpringIndex.size() == net.nrOfPartialSprings,
        "Every partial spring must be able to map to the full spring.");

      /**
       * Test maximum values
       */
      if (net.nrOfSprings > 0) {
        RUNTIME_EXP_IFN(
          net.partialToFullSpringIndex.maxCoeff() < net.nrOfSprings,
          "Partial spring must map to full spring, which must have "
          "a lower index.");
        RUNTIME_EXP_IFN(net.springPartCoordinateIndexA.maxCoeff() <
                          3 * net.nrOfLinks,
                        "Part coordinates must map to coordinates.");
        RUNTIME_EXP_IFN(net.springPartCoordinateIndexB.maxCoeff() <
                          3 * net.nrOfLinks,
                        "Part coordinates must map to coordinates.");
        RUNTIME_EXP_IFN(net.springPartIndexA.maxCoeff() < net.nrOfLinks,
                        "Part indices must map to links.");
        RUNTIME_EXP_IFN(net.springPartIndexB.maxCoeff() < net.nrOfLinks,
                        "Part indices must map to links.");
      }

      /**
       * Test reversibility of link <-> spring mapping
       */
      for (size_t link_idx = 0; link_idx < net.nrOfLinks; ++link_idx) {
        RUNTIME_EXP_IFN(
          net.linkIsEntanglement[link_idx] == (link_idx >= net.nrOfNodes),
          "Expected slip-links to come sequentially after crosslinkers.");
        std::vector<size_t> thisLinksSprings =
          net.springIndicesOfLinks[link_idx];
        std::sort(thisLinksSprings.begin(), thisLinksSprings.end());
        auto last =
          std::unique(thisLinksSprings.begin(), thisLinksSprings.end());
        RUNTIME_EXP_IFN(last == thisLinksSprings.end(),
                        "Expect each link to only have one back-link to the "
                        "springs, found back-links " +
                          pylimer_tools::utils::join(thisLinksSprings.begin(),
                                                     thisLinksSprings.end(),
                                                     std::string("_")) +
                          " for link " + std::to_string(link_idx) + ".");
        for (size_t spring_idx : thisLinksSprings) {
          std::vector<size_t> thisSpringsLinks =
            net.linkIndicesOfSprings[spring_idx];
          RUNTIME_EXP_IFN(std::find(thisSpringsLinks.begin(),
                                    thisSpringsLinks.end(),
                                    link_idx) != thisSpringsLinks.end(),
                          "Spring must have a connection to the link, too. Did "
                          "not find link " +
                            std::to_string(link_idx) + " in spring " +
                            std::to_string(spring_idx) + ".");
        }
      }

      /**
       * Test the assumptions on slip-links
       */
      for (size_t slipLinkIdx = net.nrOfNodes; slipLinkIdx < net.nrOfLinks;
           ++slipLinkIdx) {
        RUNTIME_EXP_IFN(
          net.springIndicesOfLinks[slipLinkIdx].size() == 2 ||
            net.springIndicesOfLinks[slipLinkIdx].size() == 1,
          "Expect each slip-link to be involved in exactly one or two "
          "springs, "
          "got " +
            std::to_string(net.springIndicesOfLinks[slipLinkIdx].size()) + ".");
        RUNTIME_EXP_IFN(net.linkIsEntanglement[slipLinkIdx],
                        "Expected slip-links to know what they are.");
      }

      /**
       * Test the validitiy of springs and their mapping
       */
      Eigen::ArrayXi nrOfMentions = Eigen::ArrayXi::Zero(net.nrOfLinks);
      for (size_t i = 0; i < net.nrOfSprings; ++i) {
        RUNTIME_EXP_IFN(net.linkIndicesOfSprings[i].size() >= 2,
                        "Each spring requires at least two links, got " +
                          std::to_string(net.linkIndicesOfSprings[i].size()) +
                          " at i = " + std::to_string(i) + ".");
        RUNTIME_EXP_IFN(net.springsContourLength[i] > 0,
                        "Unexpected spring contour length, got " +
                          std::to_string(net.springsContourLength[i]) +
                          " for spring " + std::to_string(i) + ".");
        RUNTIME_EXP_IFN(
          net.localToGlobalSpringIndex[i].size() ==
            net.linkIndicesOfSprings[i].size() - 1,
          "Require a global index for each local one, got " +
            std::to_string(net.localToGlobalSpringIndex[i].size()) +
            " != " + std::to_string(net.linkIndicesOfSprings[i].size() - 1) +
            " for spring " + std::to_string(i) + ".");
        for (size_t partialIdx = 0;
             partialIdx < net.localToGlobalSpringIndex[i].size();
             ++partialIdx) {
          size_t partialSpringIdx = net.localToGlobalSpringIndex[i][partialIdx];
          size_t partner0 = net.linkIndicesOfSprings[i][partialIdx];
          size_t partner1 = net.linkIndicesOfSprings[i][partialIdx + 1];
          RUNTIME_EXP_IFN(
            ((net.springPartIndexA[partialSpringIdx] == partner0 &&
              net.springPartIndexB[partialSpringIdx] == partner1)),
            "Expect linkIndicesOfSprings and localToGlobalSpringIndex "
            "ordering to correspond. Got partner0 = " +
              std::to_string(partner0) + ", partner1 = " +
              std::to_string(partner1) + " vs. springs part indices " +
              std::to_string(net.springPartIndexA[partialSpringIdx]) + " and " +
              std::to_string(net.springPartIndexB[partialSpringIdx]) +
              " in spring " + std::to_string(i) + " (partial: " +
              std::to_string(partialSpringIdx) + ") with global indices " +
              pylimer_tools::utils::join(
                net.localToGlobalSpringIndex[i].begin(),
                net.localToGlobalSpringIndex[i].end(),
                std::string(", ")) +
              " and links " +
              pylimer_tools::utils::join(net.linkIndicesOfSprings[i].begin(),
                                         net.linkIndicesOfSprings[i].end(),
                                         std::string(", ")) +
              ".");
        }
        // the following is not guaranteed anymore with the removal of links
        // while running RUNTIME_EXP_IFN(
        //   net.linkIndicesOfSprings[i][0] <=
        //     net.linkIndicesOfSprings[i][net.linkIndicesOfSprings[i].size()
        //     - 1],
        //   "Springs must have increasing end-point indices");
        std::vector<size_t> links = net.linkIndicesOfSprings[i];
        for (size_t j = 0; j < links.size(); ++j) {
          size_t link_idx = links[j];
          nrOfMentions[link_idx] += 1;
          RUNTIME_EXP_IFN(net.linkIsEntanglement[link_idx] ==
                            ((j != 0) && (j != (links.size() - 1))),
                          "Cross-links must be first and last in a spring, "
                          "slip-links in-between. Found discrepancy at " +
                            std::to_string(j) + "/" +
                            std::to_string(links.size()) + " in spring " +
                            std::to_string(i) + ".")
          std::vector<size_t> thisLinksSprings =
            net.springIndicesOfLinks[link_idx];
          RUNTIME_EXP_IFN(
            std::find(thisLinksSprings.begin(), thisLinksSprings.end(), i) !=
              thisLinksSprings.end(),
            "Link must have a connection to the spring, too. Did not find "
            "spring " +
              std::to_string(i) + " in link " + std::to_string(link_idx) + ".");
        }
      }
      for (size_t i = net.nrOfNodes; i < net.nrOfLinks; ++i) {
        RUNTIME_EXP_IFN(nrOfMentions[i] == 2,
                        "Expect each slip-link to be mentioned twice in the "
                        "links-of-springs mapping, but " +
                          std::to_string(i) + " was mentioned " +
                          std::to_string(nrOfMentions[i]) + " times.");
      }

      /**
       * Test the validity of partial springs and their mapping
       */
      for (size_t i = 0; i < net.nrOfPartialSprings; i++) {
        size_t fullIdx = net.partialToFullSpringIndex[i];
        size_t partialEndA = net.springPartIndexA[i];
        size_t partialEndB = net.springPartIndexB[i];
        RUNTIME_EXP_IFN(partialEndA < net.nrOfLinks,
                        "Cannot have a spring (" + std::to_string(i) +
                          ") part larger " + std::to_string(partialEndA) +
                          " than the nr of links (" +
                          std::to_string(net.nrOfLinks) + ").")
        RUNTIME_EXP_IFN(partialEndB < net.nrOfLinks,
                        "Cannot have a spring (" + std::to_string(i) +
                          ") part larger " + std::to_string(partialEndB) +
                          " than the nr of links (" +
                          std::to_string(net.nrOfLinks) + ").")
        RUNTIME_EXP_IFN(
          (net.linkIsEntanglement[partialEndA] ||
           net.linkIsEntanglement[partialEndB]) ==
            net.partialSpringIsPartial[i],
          "Springs involving slip-links must be marked partial. Spring " +
            std::to_string(i) + " is marked: " +
            std::to_string(net.partialSpringIsPartial[i]) + ".");
        RUNTIME_EXP_IFN(
          (net.linkIndicesOfSprings[net.partialToFullSpringIndex[i]].size() >
           2) == net.partialSpringIsPartial[i],
          "Springs involving slip-links must be marked partial. Spring " +
            std::to_string(i) + " is marked: " +
            std::to_string(net.partialSpringIsPartial[i]) + ".");
        RUNTIME_EXP_IFN(net.springPartCoordinateIndexA[3 * i] % 3 == 0,
                        "Expected spring part coordinates to be sequentially "
                        "built from spring parts.");
        RUNTIME_EXP_IFN(net.springPartCoordinateIndexB[3 * i] % 3 == 0,
                        "Expected spring part coordinates to be sequentially "
                        "built from spring parts.");
        for (int dir = 0; dir < 3; ++dir) {
          RUNTIME_EXP_IFN(
            net.springPartCoordinateIndexA[3 * i + dir] ==
              3 * partialEndA + dir,
            "Spring part index and coordinate index must match. Got " +
              std::to_string(net.springPartCoordinateIndexA[3 * i + dir]) +
              " but expected " + std::to_string(3 * partialEndA + dir) +
              " with dir = " + std::to_string(dir) + ".");
          RUNTIME_EXP_IFN(
            net.springPartCoordinateIndexB[3 * i + dir] ==
              3 * partialEndB + dir,
            "Spring part index and coordinate index must match. Got " +
              std::to_string(net.springPartCoordinateIndexB[3 * i + dir]) +
              " but expected " + std::to_string(3 * partialEndB + dir) +
              " with dir = " + std::to_string(dir) + ".");
        }
      }

      /**
       * Check that we do not have any nan or inf values in our vectors
       */
      for (long int coordI = 0; coordI < net.coordinates.size(); coordI++) {
        RUNTIME_EXP_IFN(std::isfinite(net.coordinates[coordI]),
                        "Coordinate component " + std::to_string(coordI) +
                          " must be finite, got " +
                          std::to_string(net.coordinates[coordI]) + ".");
        RUNTIME_EXP_IFN(std::isfinite(u[coordI]),
                        "Displacement component " + std::to_string(coordI) +
                          " must be finite, got " + std::to_string(u[coordI]) +
                          ".");
      }
      for (int dir = 0; dir < 3; ++dir) {
        RUNTIME_EXP_IFN(std::isfinite(net.L[dir]),
                        "Expected box size to be finite, got " +
                          std::to_string(net.L[dir]) + " in dir " +
                          std::to_string(dir) + ".");
        RUNTIME_EXP_IFN(net.L[dir] > 0.0,
                        "Expected box size to be positive, got " +
                          std::to_string(net.L[dir]) + " in dir " +
                          std::to_string(dir) + ".");
        RUNTIME_EXP_IFN(
          APPROX_EQUAL(net.boxHalfs[dir], 0.5 * net.L[dir], 1e-12),
          "Expected box half to be half of box length");
      }

      // std::cout << "Validation passed." << std::endl;
      return true;
    }
  }
}
}