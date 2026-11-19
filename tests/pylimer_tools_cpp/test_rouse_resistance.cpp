#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/sim/RouseBrownianSimulator.h"
#include "../../src/pylimer_tools_cpp/topo/ResistanceDistanceAnalysis.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

namespace ps = pylimer_tools::sim;
namespace pt = pylimer_tools::topo;

// ============================================================================
// PART 1: Topological & Resistance Distance Analysis Tests
// ============================================================================

TEST_CASE("ResistanceDistanceAnalysis linear chain properties",
          "[topo][ResistanceDistanceAnalysis]")
{
  const size_t N = 6;
  auto rda = pt::ResistanceDistanceAnalysis::createLinearChain(N);

  CHECK(rda.getNrOfNodes() == N);
  CHECK(rda.getNrOfEdges() == N - 1);
  CHECK(rda.isConnected());
  CHECK(rda.getNrOfConnectedComponents() == 1);

  SECTION("Adjacent nodes have unit effective resistance")
  {
    for (size_t i = 0; i < N - 1; ++i) {
      CHECK_THAT(rda.computeEffectiveResistance(i, i + 1),
                 Catch::Matchers::WithinAbs(1.0, 1e-9));
    }
  }

  SECTION("Series resistance between any pair is path distance |i - j|")
  {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        double expected = static_cast<double>(i > j ? i - j : j - i);
        CHECK_THAT(rda.computeEffectiveResistance(i, j),
                   Catch::Matchers::WithinAbs(expected, 1e-9));
      }
    }
  }

  SECTION("Kirchhoff index matches analytical formula N(N^2 - 1) / 6")
  {
    double expectedK = static_cast<double>(N * (N * N - 1)) / 6.0;
    CHECK_THAT(rda.computeKirchhoffIndex(),
               Catch::Matchers::WithinAbs(expectedK, 1e-8));
    CHECK_THAT(rda.computeKirchhoffIndexSpectral(),
               Catch::Matchers::WithinAbs(expectedK, 1e-8));
  }

  SECTION("Batched Multi-RHS solve matches pairwise calculations")
  {
    std::vector<std::pair<size_t, size_t>> queryPairs = { { 0, 1 }, { 0, 5 },
                                                          { 2, 4 }, { 3, 3 },
                                                          { 4, 1 }, { 5, 0 } };
    Eigen::VectorXd batchedR =
      rda.computeEffectiveResistanceBatched(queryPairs);

    REQUIRE(batchedR.size() == static_cast<Eigen::Index>(queryPairs.size()));
    for (size_t b = 0; b < queryPairs.size(); ++b) {
      double singleR = rda.computeEffectiveResistance(queryPairs[b].first,
                                                      queryPairs[b].second);
      CHECK_THAT(batchedR(static_cast<Eigen::Index>(b)),
                 Catch::Matchers::WithinAbs(singleR, 1e-9));
    }
  }

  SECTION("All-pairs matrix matches single pair queries")
  {
    Eigen::MatrixXd allR = rda.computeAllPairsEffectiveResistance();
    REQUIRE(allR.rows() == static_cast<Eigen::Index>(N));
    REQUIRE(allR.cols() == static_cast<Eigen::Index>(N));

    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        double singleR = rda.computeEffectiveResistance(i, j);
        CHECK_THAT(
          allR(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)),
          Catch::Matchers::WithinAbs(singleR, 1e-9));
      }
    }
  }

  SECTION("Commute times satisfy C_ij = 2 * M * R_ij")
  {
    double twoM = 2.0 * static_cast<double>(N - 1);
    Eigen::MatrixXd C = rda.computeAllPairsCommuteTime();

    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        double expectedC = twoM * static_cast<double>(i > j ? i - j : j - i);
        CHECK_THAT(
          C(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)),
          Catch::Matchers::WithinAbs(expectedC, 1e-8));
        CHECK_THAT(rda.computeCommuteTime(i, j),
                   Catch::Matchers::WithinAbs(expectedC, 1e-8));
      }
    }
  }

  SECTION("Resistance centrality is symmetric and highest in the center")
  {
    Eigen::VectorXd centrality = rda.computeResistanceCentrality();
    REQUIRE(centrality.size() == static_cast<Eigen::Index>(N));

    // Middle nodes (2, 3) must have higher centrality than end nodes (0, 5)
    CHECK(centrality(2) > centrality(0));
    CHECK(centrality(3) > centrality(5));
    // Symmetric nodes must have identical centrality
    CHECK_THAT(centrality(0), Catch::Matchers::WithinAbs(centrality(5), 1e-9));
    CHECK_THAT(centrality(1), Catch::Matchers::WithinAbs(centrality(4), 1e-9));
    CHECK_THAT(centrality(2), Catch::Matchers::WithinAbs(centrality(3), 1e-9));
  }
}

TEST_CASE("ResistanceDistanceAnalysis ring graph properties",
          "[topo][ResistanceDistanceAnalysis]")
{
  const size_t N = 8;
  auto rda = pt::ResistanceDistanceAnalysis::createRing(N);

  CHECK(rda.getNrOfNodes() == N);
  CHECK(rda.getNrOfEdges() == N);
  CHECK(rda.isConnected());

  SECTION("Effective resistance matches parallel path formula d*(N-d)/N")
  {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        size_t diff = (i > j) ? (i - j) : (j - i);
        size_t d = std::min(diff, N - diff);
        double expected =
          static_cast<double>(d * (N - d)) / static_cast<double>(N);
        CHECK_THAT(rda.computeEffectiveResistance(i, j),
                   Catch::Matchers::WithinAbs(expected, 1e-9));
      }
    }
  }

  SECTION("Kirchhoff index matches analytical formula (N^3 - N) / 12")
  {
    double expectedK = static_cast<double>(N * N * N - N) / 12.0;
    CHECK_THAT(rda.computeKirchhoffIndex(),
               Catch::Matchers::WithinAbs(expectedK, 1e-8));
    CHECK_THAT(rda.computeKirchhoffIndexSpectral(),
               Catch::Matchers::WithinAbs(expectedK, 1e-8));
  }

  SECTION("All nodes have identical resistance centrality due to symmetry")
  {
    Eigen::VectorXd centrality = rda.computeResistanceCentrality();
    for (size_t i = 1; i < N; ++i) {
      CHECK_THAT(centrality(static_cast<Eigen::Index>(i)),
                 Catch::Matchers::WithinAbs(centrality(0), 1e-9));
    }
  }
}

TEST_CASE("ResistanceDistanceAnalysis complete graph K_N",
          "[topo][ResistanceDistanceAnalysis]")
{
  const size_t N = 5;
  auto rda = pt::ResistanceDistanceAnalysis::createComplete(N);

  CHECK(rda.getNrOfNodes() == N);
  CHECK(rda.getNrOfEdges() == (N * (N - 1)) / 2);

  SECTION("Pairwise effective resistance is 2 / N for all distinct pairs")
  {
    double expectedR = 2.0 / static_cast<double>(N);
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        double expected = (i == j) ? 0.0 : expectedR;
        CHECK_THAT(rda.computeEffectiveResistance(i, j),
                   Catch::Matchers::WithinAbs(expected, 1e-9));
      }
    }
  }

  SECTION("Kirchhoff index is N - 1")
  {
    double expectedK = static_cast<double>(N - 1);
    CHECK_THAT(rda.computeKirchhoffIndex(),
               Catch::Matchers::WithinAbs(expectedK, 1e-8));
    CHECK_THAT(rda.computeKirchhoffIndexSpectral(),
               Catch::Matchers::WithinAbs(expectedK, 1e-8));
  }
}

TEST_CASE("ResistanceDistanceAnalysis star graph properties",
          "[topo][ResistanceDistanceAnalysis]")
{
  const size_t numArms = 4;
  auto rda = pt::ResistanceDistanceAnalysis::createStar(numArms, 1);
  const size_t N = 1 + numArms; // Center node 0 + 4 leaves

  CHECK(rda.getNrOfNodes() == N);

  SECTION("Center-to-leaf resistance is 1, leaf-to-leaf is 2")
  {
    for (size_t leaf = 1; leaf < N; ++leaf) {
      CHECK_THAT(rda.computeEffectiveResistance(0, leaf),
                 Catch::Matchers::WithinAbs(1.0, 1e-9));
    }
    for (size_t l1 = 1; l1 < N; ++l1) {
      for (size_t l2 = l1 + 1; l2 < N; ++l2) {
        CHECK_THAT(rda.computeEffectiveResistance(l1, l2),
                   Catch::Matchers::WithinAbs(2.0, 1e-9));
      }
    }
  }

  SECTION("Kirchhoff index matches (N - 1)^2")
  {
    double expectedK = static_cast<double>((N - 1) * (N - 1));
    CHECK_THAT(rda.computeKirchhoffIndex(),
               Catch::Matchers::WithinAbs(expectedK, 1e-8));
  }

  SECTION("Center has strictly higher centrality than leaves")
  {
    Eigen::VectorXd centrality = rda.computeResistanceCentrality();
    for (size_t leaf = 1; leaf < N; ++leaf) {
      CHECK(centrality(0) > centrality(static_cast<Eigen::Index>(leaf)));
    }
  }
}

TEST_CASE("ResistanceDistanceAnalysis weighted edges and parallel springs",
          "[topo][ResistanceDistanceAnalysis]")
{
  // Graph with two parallel edges between 0 and 1 with weights 1.0 and 3.0
  // Equivalent resistance = 1 / (1 + 3) = 0.25
  std::vector<size_t> from = { 0, 0 };
  std::vector<size_t> to = { 1, 1 };
  std::vector<double> weights = { 1.0, 3.0 };

  pt::ResistanceDistanceAnalysis rda(from, to, weights, 2);
  CHECK_THAT(rda.computeEffectiveResistance(0, 1),
             Catch::Matchers::WithinAbs(0.25, 1e-9));
}

TEST_CASE("ResistanceDistanceAnalysis disconnected graph handling",
          "[topo][ResistanceDistanceAnalysis]")
{
  // Two disconnected components: (0-1) and (2-3-4)
  std::vector<size_t> from = { 0, 2, 3 };
  std::vector<size_t> to = { 1, 3, 4 };
  pt::ResistanceDistanceAnalysis rda(from, to, {}, 5);

  CHECK_FALSE(rda.isConnected());
  CHECK(rda.getNrOfConnectedComponents() == 2);

  SECTION("Intra-component pairs have correct finite resistance")
  {
    CHECK_THAT(rda.computeEffectiveResistance(0, 1),
               Catch::Matchers::WithinAbs(1.0, 1e-9));
    CHECK_THAT(rda.computeEffectiveResistance(2, 4),
               Catch::Matchers::WithinAbs(2.0, 1e-9));
  }

  SECTION("Cross-component pairs have infinite resistance")
  {
    CHECK(std::isinf(rda.computeEffectiveResistance(0, 2)));
    CHECK(std::isinf(rda.computeEffectiveResistance(1, 4)));
    CHECK(std::isinf(rda.computeKirchhoffIndex()));
  }

  SECTION("Batched solve handles mixed finite and infinite queries")
  {
    std::vector<std::pair<size_t, size_t>> pairs = {
      { 0, 1 }, { 0, 3 }, { 2, 3 }, { 1, 2 }
    };
    Eigen::VectorXd res = rda.computeEffectiveResistanceBatched(pairs);

    CHECK_THAT(res(0), Catch::Matchers::WithinAbs(1.0, 1e-9));
    CHECK(std::isinf(res(1)));
    CHECK_THAT(res(2), Catch::Matchers::WithinAbs(1.0, 1e-9));
    CHECK(std::isinf(res(3)));
  }

  SECTION("Total resistance per node is infinite and centrality is zero")
  {
    Eigen::VectorXd tot = rda.computeTotalResistancePerNode();
    Eigen::VectorXd centrality = rda.computeResistanceCentrality();
    for (Eigen::Index i = 0; i < tot.size(); ++i) {
      CHECK(std::isinf(tot(i)));
      CHECK_THAT(centrality(i), Catch::Matchers::WithinAbs(0.0, 1e-12));
    }
  }
}

TEST_CASE("ResistanceDistanceAnalysis and RouseBrownianSimulator edge cases "
          "and input validation",
          "[topo][sim][ResistanceDistanceAnalysis][RouseBrownianSimulator]")
{
  SECTION("Self-loops are ignored in edge count and Laplacian assembly")
  {
    std::vector<size_t> from = { 0, 1, 1, 2 };
    std::vector<size_t> to = { 1, 1, 2, 2 };
    pt::ResistanceDistanceAnalysis rda(from, to, {}, 3);
    CHECK(rda.getNrOfNodes() == 3);
    CHECK(rda.getNrOfEdges() == 2);
    CHECK_THAT(rda.computeEffectiveResistance(0, 2),
               Catch::Matchers::WithinAbs(2.0, 1e-9));
  }

  SECTION("Empty graph is rejected in both classes")
  {
    CHECK_THROWS_AS(pt::ResistanceDistanceAnalysis({}, {}, {}, 0),
                    std::invalid_argument);
    CHECK_THROWS_AS(ps::RouseBrownianSimulator({}, {}, 1.0, 1.0, 1.0, 0.01, 0),
                    std::invalid_argument);
  }

  SECTION("Asymmetric and negative-diagonal Laplacians are rejected")
  {
    Eigen::SparseMatrix<double> asymL(2, 2);
    asymL.insert(0, 0) = 1.0;
    asymL.insert(0, 1) = -1.0;
    asymL.insert(1, 1) = 1.0;
    asymL.makeCompressed();

    CHECK_THROWS_AS(pt::ResistanceDistanceAnalysis(asymL),
                    std::invalid_argument);
    CHECK_THROWS_AS(ps::RouseBrownianSimulator(asymL), std::invalid_argument);

    Eigen::SparseMatrix<double> negDiagL(2, 2);
    negDiagL.insert(0, 0) = -1.0;
    negDiagL.insert(1, 1) = 1.0;
    negDiagL.makeCompressed();

    CHECK_THROWS_AS(pt::ResistanceDistanceAnalysis(negDiagL),
                    std::invalid_argument);
    CHECK_THROWS_AS(ps::RouseBrownianSimulator(negDiagL),
                    std::invalid_argument);
  }
}

TEST_CASE("ResistanceDistanceAnalysis Moore-Penrose pseudo-inverse properties",
          "[topo][ResistanceDistanceAnalysis]")
{
  const size_t N = 5;
  auto rda = pt::ResistanceDistanceAnalysis::createLinearChain(N);
  Eigen::MatrixXd Lplus = rda.computeMoorePenrosePseudoInverse();
  Eigen::MatrixXd L = Eigen::MatrixXd(rda.getLaplacian());

  SECTION("Pseudo-inverse projects out the all-ones null mode: L^+ * 1 = 0")
  {
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd Lplus_ones = Lplus * ones;
    CHECK_THAT(Lplus_ones.norm(), Catch::Matchers::WithinAbs(0.0, 1e-10));
  }

  SECTION("First Moore-Penrose condition: L * L^+ * L = L")
  {
    Eigen::MatrixXd cond1 = L * Lplus * L;
    CHECK((cond1 - L).norm() < 1e-9);
  }

  SECTION("Second Moore-Penrose condition: L^+ * L * L^+ = L^+")
  {
    Eigen::MatrixXd cond2 = Lplus * L * Lplus;
    CHECK((cond2 - Lplus).norm() < 1e-9);
  }

  SECTION("Resistance distance equals L^+_ii + L^+_jj - 2 * L^+_ij")
  {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        double expected =
          Lplus(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) +
          Lplus(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(j)) -
          2.0 *
            Lplus(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
        CHECK_THAT(rda.computeEffectiveResistance(i, j),
                   Catch::Matchers::WithinAbs(expected, 1e-9));
      }
    }
  }

  SECTION("Kirchhoff index equals N * Tr(L^+)")
  {
    double traceK = static_cast<double>(N) * Lplus.trace();
    CHECK_THAT(rda.computeKirchhoffIndex(),
               Catch::Matchers::WithinAbs(traceK, 1e-9));
  }
}

// ============================================================================
// PART 2: Polymer Brownian & Rouse Dynamics Tests
// ============================================================================

TEST_CASE("RouseBrownianSimulator deterministic relaxation",
          "[sim][RouseBrownianSimulator]")
{
  const size_t N = 8;
  double gamma = 2.0;
  double k = 1.5;
  double T = 0.0; // Deterministic (no thermal noise)
  double dt = 0.05;

  auto sim = ps::RouseBrownianSimulator::createLinearChain(N, gamma, k, T, dt);
  sim.initializeLinearPositions(1.0, 0); // Along x-axis: [0, 1, 2, ..., N-1]

  Eigen::Vector3d initialCM = sim.getCenterOfMass();
  double initialEnergy = sim.getPotentialEnergy();
  CHECK(initialEnergy > 0.0);

  SECTION("Center of mass remains invariant during deterministic relaxation")
  {
    sim.stepDeterministic(50);
    Eigen::Vector3d curCM = sim.getCenterOfMass();
    CHECK_THAT((curCM - initialCM).norm(),
               Catch::Matchers::WithinAbs(0.0, 1e-12));
  }

  SECTION("Potential energy decreases monotonically during relaxation")
  {
    double prevEnergy = initialEnergy;
    for (size_t s = 0; s < 20; ++s) {
      sim.stepDeterministic(1);
      double curEnergy = sim.getPotentialEnergy();
      CHECK(curEnergy <= prevEnergy + 1e-14);
      prevEnergy = curEnergy;
    }
  }

  SECTION("Normal mode relaxation rate matches exact theoretical eigenvalue")
  {
    // Find non-zero eigenvector of Laplacian L
    Eigen::MatrixXd denseL = Eigen::MatrixXd(sim.getLaplacian());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(denseL);
    Eigen::VectorXd eigvals = es.eigenvalues();
    Eigen::MatrixXd eigvecs = es.eigenvectors();

    // Mode p = 1 (first non-zero mode)
    int modeIdx = 1;
    double lambda1 = eigvals(modeIdx);
    REQUIRE(lambda1 > 1e-6);

    Eigen::VectorXd mode1 = eigvecs.col(modeIdx);
    Eigen::MatrixXd pos = Eigen::MatrixXd::Zero(N, 3);
    pos.col(0) = mode1;
    sim.setPositions(pos);

    // After 1 step: r_{t+dt} = (1 + dt * (k / gamma) * lambda1)^{-1} * r_t
    double theoreticalDecayFactor = 1.0 / (1.0 + dt * (k / gamma) * lambda1);
    sim.stepDeterministic(1);

    Eigen::MatrixXd relaxedPos = sim.getPositions();
    double observedDecayFactor = relaxedPos.col(0).norm() / mode1.norm();

    CHECK_THAT(observedDecayFactor,
               Catch::Matchers::WithinAbs(theoreticalDecayFactor, 1e-10));
  }

  SECTION(
    "Direct LDLT, Direct LLT, and Conjugate Gradient yield identical steps")
  {
    auto simLDLT =
      ps::RouseBrownianSimulator::createLinearChain(N, gamma, k, T, dt);
    auto simLLT =
      ps::RouseBrownianSimulator::createLinearChain(N, gamma, k, T, dt);
    auto simCG =
      ps::RouseBrownianSimulator::createLinearChain(N, gamma, k, T, dt);

    simLDLT.setSolverType(ps::RouseSolverType::DIRECT_LDLT);
    simLLT.setSolverType(ps::RouseSolverType::DIRECT_LLT);
    simCG.setSolverType(ps::RouseSolverType::CONJUGATE_GRADIENT);

    simLDLT.initializeLinearPositions(1.0, 0);
    simLLT.initializeLinearPositions(1.0, 0);
    simCG.initializeLinearPositions(1.0, 0);

    simLDLT.stepDeterministic(10);
    simLLT.stepDeterministic(10);
    simCG.stepDeterministic(10);

    CHECK((simLDLT.getPositions() - simLLT.getPositions()).norm() < 1e-12);
    CHECK((simLDLT.getPositions() - simCG.getPositions()).norm() < 1e-5);
  }
}

TEST_CASE("RouseBrownianSimulator center-of-mass diffusion",
          "[sim][RouseBrownianSimulator]")
{
  const size_t N = 10;
  double gamma = 1.0;
  double k = 2.0;
  double T = 1.0;
  double dt = 0.01;

  auto sim = ps::RouseBrownianSimulator::createLinearChain(N, gamma, k, T, dt);
  sim.initializeLinearPositions(1.0, 0);
  sim.setSeed(12345);

  size_t totalSteps = 500;
  size_t sampleInterval = 10;
  auto result = sim.run(totalSteps, sampleInterval, true);

  REQUIRE(result.times.size() ==
          static_cast<Eigen::Index>(totalSteps / sampleInterval + 1));
  REQUIRE(result.msdCenterOfMass.size() == result.times.size());

  SECTION("Center of mass MSD increases with time")
  {
    CHECK(result.msdCenterOfMass(result.times.size() - 1) >
          result.msdCenterOfMass(0));
  }

  SECTION(
    "Post-processing computeMSD and computeCenterOfMassMSD match trajectory")
  {
    Eigen::VectorXd msd =
      ps::RouseBrownianSimulator::computeMSD(result.trajectory, false);
    Eigen::VectorXd cmMsd = ps::RouseBrownianSimulator::computeCenterOfMassMSD(
      result.trajectory, false);

    CHECK_THAT(msd(0), Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK_THAT(cmMsd(0), Catch::Matchers::WithinAbs(0.0, 1e-12));

    for (int i = 0; i < result.times.size(); ++i) {
      CHECK_THAT(msd(i), Catch::Matchers::WithinAbs(result.msd(i), 1e-10));
      CHECK_THAT(cmMsd(i),
                 Catch::Matchers::WithinAbs(result.msdCenterOfMass(i), 1e-10));
    }
  }

  SECTION("Time-averaged end-to-end autocorrelation decays from 1.0")
  {
    Eigen::VectorXd corr =
      ps::RouseBrownianSimulator::computeEndToEndAutocorrelation(
        result.trajectory, 0, N - 1, true);
    CHECK_THAT(corr(0), Catch::Matchers::WithinAbs(1.0, 1e-10));
  }
}

TEST_CASE("RouseBrownianSimulator external force drift",
          "[sim][RouseBrownianSimulator]")
{
  const size_t N = 5;
  double gamma = 2.0;
  double k = 1.0;
  double T = 0.0; // Deterministic
  double dt = 0.02;

  auto sim = ps::RouseBrownianSimulator::createLinearChain(N, gamma, k, T, dt);
  sim.initializeLinearPositions(1.0, 0);

  // Apply a constant force along the x-axis on node 0: F = (10, 0, 0)
  Eigen::MatrixXd forces = Eigen::MatrixXd::Zero(N, 3);
  forces(0, 0) = 10.0;

  Eigen::Vector3d initialCM = sim.getCenterOfMass();
  size_t steps = 100;
  sim.stepDeterministicWithExternalForces(forces, steps);

  // Theoretical center of mass drift:
  // Total external force = 10.0. Total friction = N * gamma = 5 * 2 = 10.
  // Terminal drift velocity v_CM = F_tot / (N * gamma) = 10.0 / 10.0 = 1.0
  // Displacement in time t = steps * dt = 100 * 0.02 = 2.0 is Delta x = 2.0
  // * 1.0 = 2.0
  double totalTime = static_cast<double>(steps) * dt;
  double expectedDisplacement =
    (10.0 / (static_cast<double>(N) * gamma)) * totalTime;

  Eigen::Vector3d finalCM = sim.getCenterOfMass();
  double observedDisplacement = finalCM(0) - initialCM(0);

  CHECK_THAT(observedDisplacement,
             Catch::Matchers::WithinAbs(expectedDisplacement, 1e-8));
  // y and z should remain 0
  CHECK_THAT(finalCM(1) - initialCM(1), Catch::Matchers::WithinAbs(0.0, 1e-12));
  CHECK_THAT(finalCM(2) - initialCM(2), Catch::Matchers::WithinAbs(0.0, 1e-12));
}

TEST_CASE("RouseBrownianSimulator coordinates and parameter manipulation",
          "[sim][RouseBrownianSimulator]")
{
  const size_t N = 4;
  auto sim = ps::RouseBrownianSimulator::createLinearChain(N);

  SECTION("Flattened position conversion: interleaved vs non-interleaved")
  {
    Eigen::MatrixXd pos(N, 3);
    pos << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0;

    sim.setPositions(pos);

    // Interleaved
    Eigen::VectorXd flatInter = sim.getPositionsFlattened(true);
    REQUIRE(flatInter.size() == 12);
    CHECK(flatInter(0) == 1.0);
    CHECK(flatInter(1) == 2.0);
    CHECK(flatInter(2) == 3.0);
    CHECK(flatInter(11) == 12.0);

    // Reconstruct via setPositionsFlattened
    sim.initializeLinearPositions();
    sim.setPositionsFlattened(flatInter, true);
    CHECK((sim.getPositions() - pos).norm() < 1e-12);

    // Non-interleaved
    Eigen::VectorXd flatNonInter = sim.getPositionsFlattened(false);
    REQUIRE(flatNonInter.size() == 12);
    CHECK(flatNonInter(0) == 1.0);
    CHECK(flatNonInter(4) == 2.0);
    CHECK(flatNonInter(8) == 3.0);

    sim.initializeLinearPositions();
    sim.setPositionsFlattened(flatNonInter, false);
    CHECK((sim.getPositions() - pos).norm() < 1e-12);
  }

  SECTION("Dynamic parameter updates and validation")
  {
    CHECK_THROWS_AS(sim.setGamma(-1.0), std::invalid_argument);
    CHECK_THROWS_AS(sim.setK(0.0), std::invalid_argument);
    CHECK_THROWS_AS(sim.setDt(-0.1), std::invalid_argument);
    CHECK_THROWS_AS(sim.setTemperature(-1.0), std::invalid_argument);

    sim.setParameters(3.0, 4.0, 2.5, 0.05);
    CHECK(sim.getGamma() == 3.0);
    CHECK(sim.getK() == 4.0);
    CHECK(sim.getTemperature() == 2.5);
    CHECK(sim.getDt() == 0.05);

    // Diffusion coefficient D = kB * T / gamma = 1.0 * 2.5 / 3.0
    CHECK_THAT(sim.getDiffusionCoefficient(),
               Catch::Matchers::WithinAbs(2.5 / 3.0, 1e-10));
  }

  SECTION("Radius of gyration and end-to-end calculations")
  {
    sim.initializeLinearPositions(1.0, 0); // Beads at x = 0, 1, 2, 3
    CHECK_THAT(sim.getEndToEndDistance(0, 3),
               Catch::Matchers::WithinAbs(3.0, 1e-10));

    // CM is at (1.5, 0, 0)
    // Displacements from CM: -1.5, -0.5, 0.5, 1.5
    // Sum of squares = 2.25 + 0.25 + 0.25 + 2.25 = 5.0
    // Rg = sqrt(5.0 / 4) = sqrt(1.25)
    double expectedRg = std::sqrt(1.25);
    CHECK_THAT(sim.getRadiusOfGyration(),
               Catch::Matchers::WithinAbs(expectedRg, 1e-10));
  }
}

TEST_CASE("ResistanceDistanceAnalysis Level-3 BLAS Chunking for All-Pairs & "
          "Pseudo-Inverse",
          "[ResistanceDistanceAnalysis][topo][benchmark]")
{
  const size_t N = 300;
  std::vector<size_t> from, to;
  std::vector<double> weights;
  for (size_t i = 0; i < N - 1; ++i) {
    from.push_back(i);
    to.push_back(i + 1);
    weights.push_back(1.0);
  }
  from.push_back(N - 1);
  to.push_back(0);
  weights.push_back(1.0);

  // Cross-links
  for (size_t i = 0; i < N / 4; ++i) {
    from.push_back(i * 3);
    to.push_back((i * 3 + N / 2) % N);
    weights.push_back(0.5);
  }

  pt::ResistanceDistanceAnalysis rda(from, to, weights, N);

  SECTION(
    "Chunked solves produce numerically identical results across chunk sizes")
  {
    Eigen::MatrixXd R_chunk1 = rda.computeAllPairsEffectiveResistance(1);
    Eigen::MatrixXd R_chunk32 = rda.computeAllPairsEffectiveResistance(32);
    Eigen::MatrixXd R_chunk64 = rda.computeAllPairsEffectiveResistance(64);
    Eigen::MatrixXd R_chunk128 = rda.computeAllPairsEffectiveResistance(128);

    double diff32 = (R_chunk1 - R_chunk32).cwiseAbs().maxCoeff();
    double diff64 = (R_chunk1 - R_chunk64).cwiseAbs().maxCoeff();
    double diff128 = (R_chunk1 - R_chunk128).cwiseAbs().maxCoeff();

    CHECK(diff32 < 1e-12);
    CHECK(diff64 < 1e-12);
    CHECK(diff128 < 1e-12);

    Eigen::MatrixXd P_chunk1 = rda.computeMoorePenrosePseudoInverse(1);
    Eigen::MatrixXd P_chunk128 = rda.computeMoorePenrosePseudoInverse(128);
    double pDiff = (P_chunk1 - P_chunk128).cwiseAbs().maxCoeff();
    CHECK(pDiff < 1e-12);
  }
}

TEST_CASE(
  "ResistanceDistanceAnalysis Level-3 BLAS Chunking Benchmark (B=1 vs B=128)",
  "[ResistanceDistanceAnalysis][topo][benchmark][long]")
{
  const size_t N = 300;
  std::vector<size_t> from, to;
  std::vector<double> weights;
  for (size_t i = 0; i < N - 1; ++i) {
    from.push_back(i);
    to.push_back(i + 1);
    weights.push_back(1.0);
  }
  from.push_back(N - 1);
  to.push_back(0);
  weights.push_back(1.0);

  // Cross-links
  for (size_t i = 0; i < N / 4; ++i) {
    from.push_back(i * 3);
    to.push_back((i * 3 + N / 2) % N);
    weights.push_back(0.5);
  }

  pt::ResistanceDistanceAnalysis rda(from, to, weights, N);

  auto start1 = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd R1 = rda.computeAllPairsEffectiveResistance(1);
  auto end1 = std::chrono::high_resolution_clock::now();
  double time1_ms =
    std::chrono::duration<double, std::milli>(end1 - start1).count();

  auto start128 = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd R128 = rda.computeAllPairsEffectiveResistance(128);
  auto end128 = std::chrono::high_resolution_clock::now();
  double time128_ms =
    std::chrono::duration<double, std::milli>(end128 - start128).count();

  double speedup = time1_ms / time128_ms;
  std::cout << "\n=== All-Pairs Resistance Distance Level-3 BLAS Chunking "
               "Benchmark (N="
            << N << ") ===\n"
            << "Level-2 BLAS (chunk=1, dtrsv):   " << time1_ms << " ms\n"
            << "Level-3 BLAS (chunk=128, dtrsm): " << time128_ms << " ms\n"
            << "Speedup:                         " << speedup << "x\n"
            << "============================================================="
               "============\n";

  CHECK(speedup >= 1.0);
}

TEST_CASE("ResistanceDistanceAnalysis and RouseBrownianSimulator Universe, "
          "fromMatrix, and validation extensions",
          "[topo][sim][ResistanceDistanceAnalysis][RouseBrownianSimulator]")
{
  SECTION("Construction from Universe and fromMatrix")
  {
    pylimer_tools::entities::Universe u(10.0, 10.0, 10.0);
    u.addAtoms({ 1, 2, 3, 4 },
               { 1, 1, 1, 1 },
               { 0.0, 1.5, 3.0, 4.5 },
               { 0.0, 0.0, 0.0, 0.0 },
               { 0.0, 0.0, 0.0, 0.0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 },
               { 0, 0, 0, 0 });
    u.addBonds({ 1, 2, 3 }, { 2, 3, 4 });

    pt::ResistanceDistanceAnalysis rdaU(u);
    pt::ResistanceDistanceAnalysis rdaFromU =
      pt::ResistanceDistanceAnalysis::fromUniverse(u);
    CHECK(rdaU.getNrOfNodes() == 4);
    CHECK(rdaU.getNrOfEdges() == 3);
    CHECK_THAT(rdaU.computeEffectiveResistance(0, 3),
               Catch::Matchers::WithinAbs(3.0, 1e-9));
    CHECK_THAT(rdaFromU.computeKirchhoffIndex(),
               Catch::Matchers::WithinAbs(10.0, 1e-8));

    pt::ResistanceDistanceAnalysis rdaFromM =
      pt::ResistanceDistanceAnalysis::fromMatrix(rdaU.getLaplacian());
    CHECK_THAT(rdaFromM.computeEffectiveResistance(0, 3),
               Catch::Matchers::WithinAbs(3.0, 1e-9));

    ps::RouseBrownianSimulator simU(u, 1.0, 2.0, 0.5, 0.01);
    ps::RouseBrownianSimulator simFromU =
      ps::RouseBrownianSimulator::fromUniverse(u, 1.0, 2.0, 0.5, 0.01);
    CHECK(simU.getNrOfNodes() == 4);
    CHECK_THAT(simU.getEndToEndDistance(),
               Catch::Matchers::WithinAbs(4.5, 1e-12));
    CHECK((simU.getPositions() - simFromU.getPositions()).norm() < 1e-12);

    ps::RouseBrownianSimulator simFromM =
      ps::RouseBrownianSimulator::fromMatrix(simU.getLaplacian(), 1.0, 2.0);
    CHECK(simFromM.getNrOfNodes() == 4);
  }

  SECTION("Reject positive off-diagonals and non-zero row/column sums")
  {
    Eigen::SparseMatrix<double> adj(2, 2);
    adj.insert(0, 1) = 1.0;
    adj.insert(1, 0) = 1.0;
    adj.makeCompressed();
    CHECK_THROWS_AS(pt::ResistanceDistanceAnalysis(adj), std::invalid_argument);
    CHECK_THROWS_AS(ps::RouseBrownianSimulator(adj), std::invalid_argument);

    Eigen::SparseMatrix<double> shiftedL(2, 2);
    shiftedL.insert(0, 0) = 2.0;
    shiftedL.insert(0, 1) = -1.0;
    shiftedL.insert(1, 0) = -1.0;
    shiftedL.insert(1, 1) = 2.0;
    shiftedL.makeCompressed();
    CHECK_THROWS_AS(pt::ResistanceDistanceAnalysis(shiftedL),
                    std::invalid_argument);
    CHECK_THROWS_AS(ps::RouseBrownianSimulator(shiftedL),
                    std::invalid_argument);
  }

  SECTION("Explicit sparse zeros are not counted as edges")
  {
    Eigen::SparseMatrix<double> sparseWithZeros(3, 3);
    sparseWithZeros.insert(0, 0) = 1.0;
    sparseWithZeros.insert(0, 1) = -1.0;
    sparseWithZeros.insert(1, 0) = -1.0;
    sparseWithZeros.insert(1, 1) = 1.0;
    sparseWithZeros.insert(1, 2) = 0.0;
    sparseWithZeros.insert(2, 1) = 0.0;
    sparseWithZeros.insert(2, 2) = 0.0;
    sparseWithZeros.makeCompressed();

    pt::ResistanceDistanceAnalysis rda(sparseWithZeros);
    CHECK(rda.getNrOfNodes() == 3);
    CHECK(rda.getNrOfEdges() == 1);
    CHECK_FALSE(rda.isConnected());
    CHECK_THAT(rda.computeEffectiveResistance(0, 1),
               Catch::Matchers::WithinAbs(1.0, 1e-9));
    CHECK(std::isinf(rda.computeEffectiveResistance(0, 2)));
  }

  SECTION("Disconnected Moore-Penrose pseudo-inverse properties")
  {
    std::vector<size_t> from = { 0, 1, 3, 4 };
    std::vector<size_t> to = { 1, 2, 4, 3 };
    pt::ResistanceDistanceAnalysis rda(from, to, {}, 6);
    Eigen::MatrixXd Lplus = rda.computeMoorePenrosePseudoInverse(2);
    Eigen::MatrixXd L = Eigen::MatrixXd(rda.getLaplacian());

    CHECK((L * Lplus * L - L).norm() < 1e-9);
    CHECK((Lplus * L * Lplus - Lplus).norm() < 1e-9);
    CHECK((Lplus * Eigen::VectorXd::Ones(6)).norm() < 1e-10);
    CHECK_THAT(Lplus(0, 0) + Lplus(2, 2) - 2.0 * Lplus(0, 2),
               Catch::Matchers::WithinAbs(2.0, 1e-9));
    CHECK_THAT(Lplus(5, 5), Catch::Matchers::WithinAbs(0.0, 1e-12));
  }

  SECTION("Total resistance per node matches row sums across chunk sizes")
  {
    const pt::ResistanceDistanceAnalysis rda =
      pt::ResistanceDistanceAnalysis::createLinearChain(20);
    Eigen::MatrixXd allR = rda.computeAllPairsEffectiveResistance(5);
    Eigen::VectorXd expectedRowSums = allR.rowwise().sum();
    Eigen::VectorXd tot1 = rda.computeTotalResistancePerNode(1);
    Eigen::VectorXd tot7 = rda.computeTotalResistancePerNode(7);
    CHECK((tot1 - expectedRowSums).cwiseAbs().maxCoeff() < 1e-9);
    CHECK((tot7 - expectedRowSums).cwiseAbs().maxCoeff() < 1e-9);
  }

  SECTION("RouseBrownianSimulator createComplete, weighted springs, and "
          "trajectory validation")
  {
    ps::RouseBrownianSimulator simComplete =
      ps::RouseBrownianSimulator::createComplete(4, 1.0, 2.0);
    CHECK(simComplete.getNrOfNodes() == 4);
    simComplete.initializeLinearPositions(1.0, 0);
    // All 6 pairs (0,1),(0,2),(0,3),(1,2),(1,3),(2,3) have dx = 1,2,3,1,2,1
    // sum dx^2 = 1 + 4 + 9 + 1 + 4 + 1 = 20 -> U = 0.5 * 2.0 * 20 = 20.0
    CHECK_THAT(simComplete.getPotentialEnergy(),
               Catch::Matchers::WithinAbs(20.0, 1e-10));

    ps::RouseBrownianSimulator simWeighted(
      { 0 }, { 1 }, 1.0, 1.0, 0.0, 0.01, 2, { 3.0 });
    simWeighted.initializeLinearPositions(2.0, 0);
    // U = 0.5 * k * w * r^2 = 0.5 * 1.0 * 3.0 * 4.0 = 6.0
    CHECK_THAT(simWeighted.getPotentialEnergy(),
               Catch::Matchers::WithinAbs(6.0, 1e-10));

    std::vector<Eigen::MatrixXd> emptyTraj;
    CHECK_THROWS_AS(ps::RouseBrownianSimulator::computeMSD(emptyTraj),
                    std::invalid_argument);
    std::vector<Eigen::MatrixXd> badCols = { Eigen::MatrixXd::Zero(3, 2) };
    CHECK_THROWS_AS(ps::RouseBrownianSimulator::computeCenterOfMassMSD(badCols),
                    std::invalid_argument);
    std::vector<Eigen::MatrixXd> mismatchedFrames = {
      Eigen::MatrixXd::Zero(3, 3), Eigen::MatrixXd::Zero(4, 3)
    };
    CHECK_THROWS_AS(ps::RouseBrownianSimulator::computeEndToEndAutocorrelation(
                      mismatchedFrames),
                    std::invalid_argument);
  }

  SECTION("Disconnected graph commute times use per-component edge weights")
  {
    // Component 1: 3-node chain (0-1-2) with weight 1.0 -> M_1 = 2.0, 2*M_1 = 4
    // Component 2: 2-node edge (3-4) with weight 5.0 -> M_2 = 5.0, 2*M_2 = 10
    const std::vector<size_t> from = { 0, 1, 3 };
    const std::vector<size_t> to = { 1, 2, 4 };
    const std::vector<double> weights = { 1.0, 1.0, 5.0 };
    const pt::ResistanceDistanceAnalysis rda(from, to, weights, 5);

    REQUIRE_FALSE(rda.isConnected());
    REQUIRE(rda.getNrOfConnectedComponents() == 2);

    // Per-component commute times: C_ij = 2 * M_C * R_ij
    CHECK_THAT(rda.computeCommuteTime(0, 1),
               Catch::Matchers::WithinAbs(4.0, 1e-9));
    CHECK_THAT(rda.computeCommuteTime(0, 2),
               Catch::Matchers::WithinAbs(8.0, 1e-9));
    CHECK_THAT(rda.computeCommuteTime(1, 2),
               Catch::Matchers::WithinAbs(4.0, 1e-9));
    // R_34 = 1 / 5 = 0.2 -> C_34 = 2 * 5 * 0.2 = 2.0
    CHECK_THAT(rda.computeCommuteTime(3, 4),
               Catch::Matchers::WithinAbs(2.0, 1e-9));
    CHECK_THAT(rda.computeCommuteTime(0, 0),
               Catch::Matchers::WithinAbs(0.0, 1e-12));
    CHECK(std::isinf(rda.computeCommuteTime(0, 3)));

    const Eigen::MatrixXd allC = rda.computeAllPairsCommuteTime(2);
    REQUIRE(allC.rows() == 5);
    REQUIRE(allC.cols() == 5);
    CHECK_THAT(allC(0, 1), Catch::Matchers::WithinAbs(4.0, 1e-9));
    CHECK_THAT(allC(0, 2), Catch::Matchers::WithinAbs(8.0, 1e-9));
    CHECK_THAT(allC(1, 2), Catch::Matchers::WithinAbs(4.0, 1e-9));
    CHECK_THAT(allC(3, 4), Catch::Matchers::WithinAbs(2.0, 1e-9));
    CHECK_THAT(allC(4, 3), Catch::Matchers::WithinAbs(2.0, 1e-9));
    CHECK(std::isinf(allC(0, 3)));
    CHECK(std::isinf(allC(2, 4)));
  }

  SECTION("Small-weight Laplacian preserves connectivity and resistance scale")
  {
    const double w = 1e-15;
    const std::vector<size_t> from = { 0, 1, 2 };
    const std::vector<size_t> to = { 1, 2, 3 };
    const std::vector<double> weights = { w, w, w };
    const pt::ResistanceDistanceAnalysis rda(from, to, weights, 4);

    CHECK(rda.getNrOfNodes() == 4);
    CHECK(rda.getNrOfEdges() == 3);
    CHECK(rda.isConnected());
    CHECK(rda.getNrOfConnectedComponents() == 1);
    CHECK_THAT(rda.computeEffectiveResistance(0, 1),
               Catch::Matchers::WithinRel(1e15, 1e-9));
    CHECK_THAT(rda.computeEffectiveResistance(0, 3),
               Catch::Matchers::WithinRel(3e15, 1e-9));

    const double kirchhoff = rda.computeKirchhoffIndex();
    const double kirchhoffSpectral = rda.computeKirchhoffIndexSpectral();
    CHECK_THAT(kirchhoff, Catch::Matchers::WithinRel(10.0 * 1e15, 1e-9));
    CHECK_THAT(kirchhoffSpectral, Catch::Matchers::WithinRel(kirchhoff, 1e-9));
  }

  SECTION("Chunk size handling and validation in ResistanceDistanceAnalysis")
  {
    const pt::ResistanceDistanceAnalysis rda =
      pt::ResistanceDistanceAnalysis::createLinearChain(6);
    const std::vector<std::pair<size_t, size_t>> pairs = {
      { 0, 1 }, { 0, 5 }, { 1, 4 }, { 2, 5 }, { 3, 3 }
    };
    const Eigen::VectorXd defaultBatched =
      rda.computeEffectiveResistanceBatched(pairs);
    const Eigen::VectorXd chunkedBatched =
      rda.computeEffectiveResistanceBatched(pairs, 2);
    CHECK((chunkedBatched - defaultBatched).cwiseAbs().maxCoeff() < 1e-9);

    CHECK_THROWS_AS(rda.computeEffectiveResistanceBatched(pairs, 0),
                    std::invalid_argument);
    CHECK_THROWS_AS(rda.computeResistanceCentrality(0), std::invalid_argument);
    CHECK_THROWS_AS(rda.computeAllPairsCommuteTime(0), std::invalid_argument);
  }

  SECTION("RouseBrownianSimulator setSeed is reproducible after odd normal "
          "draws")
  {
    // 3 nodes * 3 dimensions = 9 normal draws per step (odd)
    ps::RouseBrownianSimulator simA =
      ps::RouseBrownianSimulator::createLinearChain(3, 1.0, 1.0, 1.0, 0.01);
    simA.initializeLinearPositions(1.0, 0);
    simA.setSeed(999);
    simA.step(1);

    simA.initializeLinearPositions(1.0, 0);
    simA.setSeed(123);

    ps::RouseBrownianSimulator simB =
      ps::RouseBrownianSimulator::createLinearChain(3, 1.0, 1.0, 1.0, 0.01);
    simB.initializeLinearPositions(1.0, 0);
    simB.setSeed(123);

    simA.step(5);
    simB.step(5);
    CHECK((simA.getPositions() - simB.getPositions()).norm() < 1e-12);
  }

  SECTION("RouseBrownianSimulator rejects non-finite inputs and invalid "
          "weights")
  {
    ps::RouseBrownianSimulator sim =
      ps::RouseBrownianSimulator::createLinearChain(3, 1.0, 1.0, 0.0, 0.01);
    const double nanVal = std::numeric_limits<double>::quiet_NaN();
    const double infVal = std::numeric_limits<double>::infinity();

    Eigen::MatrixXd badPos = Eigen::MatrixXd::Zero(3, 3);
    badPos(1, 2) = nanVal;
    CHECK_THROWS_AS(sim.setPositions(badPos), std::invalid_argument);
    badPos(1, 2) = infVal;
    CHECK_THROWS_AS(sim.setPositions(badPos), std::invalid_argument);

    Eigen::VectorXd badFlat = Eigen::VectorXd::Zero(9);
    badFlat(4) = nanVal;
    CHECK_THROWS_AS(sim.setPositionsFlattened(badFlat), std::invalid_argument);
    badFlat(4) = infVal;
    CHECK_THROWS_AS(sim.setPositionsFlattened(badFlat), std::invalid_argument);

    Eigen::MatrixXd badForces = Eigen::MatrixXd::Zero(3, 3);
    badForces(0, 0) = nanVal;
    CHECK_THROWS_AS(sim.stepWithExternalForces(badForces),
                    std::invalid_argument);
    CHECK_THROWS_AS(sim.stepDeterministicWithExternalForces(badForces),
                    std::invalid_argument);
    badForces(0, 0) = infVal;
    CHECK_THROWS_AS(sim.stepWithExternalForces(badForces),
                    std::invalid_argument);
    CHECK_THROWS_AS(sim.stepDeterministicWithExternalForces(badForces),
                    std::invalid_argument);

    CHECK_THROWS_AS(ps::RouseBrownianSimulator(
                      { 0 }, { 1 }, 1.0, 1.0, 0.0, 0.01, 2, { -1.0 }),
                    std::invalid_argument);
    CHECK_THROWS_AS(
      ps::RouseBrownianSimulator({ 0 }, { 1 }, 1.0, 1.0, 0.0, 0.01, 2, { 0.0 }),
      std::invalid_argument);
    CHECK_THROWS_AS(ps::RouseBrownianSimulator(
                      { 0 }, { 1 }, 1.0, 1.0, 0.0, 0.01, 2, { nanVal }),
                    std::invalid_argument);
    CHECK_THROWS_AS(ps::RouseBrownianSimulator(
                      { 0 }, { 1 }, 1.0, 1.0, 0.0, 0.01, 2, { 1.0, 2.0 }),
                    std::invalid_argument);
  }
}
