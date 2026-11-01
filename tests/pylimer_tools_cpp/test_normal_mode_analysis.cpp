#include "../../src/pylimer_tools_cpp/calc/NormalModeAnalyzer.h"

#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace pc = pylimer_tools::calc;

TEST_CASE("NormalModeAnalyzer does not crash", "[analysis][NormalModeAnalyzer]")
{
  // very small melt chain
  pc::NormalModeAnalyzer normalModeAnalyzer =
    pc::NormalModeAnalyzer({ 1, 2, 3 }, { 2, 3, 4 });

  REQUIRE_THROWS(normalModeAnalyzer.getEigenvalues());
  REQUIRE_THROWS(normalModeAnalyzer.getEigenvectors());
  REQUIRE_THROWS(normalModeAnalyzer.getNrOfSolubleClusters());
  Eigen::ArrayXd testGrid = Eigen::ArrayXd::LinSpaced(5, 0.1, 1.0);
  REQUIRE_THROWS(normalModeAnalyzer.evaluateStressAutocorrelation(testGrid));
  REQUIRE_THROWS(normalModeAnalyzer.evaluateStorageModulus(testGrid));
  REQUIRE_THROWS(normalModeAnalyzer.evaluateLossModulus(testGrid));
  REQUIRE_NOTHROW(normalModeAnalyzer.computeAllEigenvalues());
  Eigen::VectorXd eigenvalues = normalModeAnalyzer.getEigenvalues();
  CHECK(eigenvalues.size() == 5);
}

TEST_CASE("NormalModeAnalyzer constructor supports explicit nrOfNodes for "
          "isolated vertices",
          "[analysis][NormalModeAnalyzer]")
{
  // springFrom and springTo only reference indices 0, 1, 2. Max index is 2.
  // Explicitly specify nrOfNodes = 6.
  pc::NormalModeAnalyzer nma({ 0, 1 }, { 1, 2 }, 6);
  CHECK(nma.getMatrixSize() == 6);
  nma.computeAllEigenvalues();
  CHECK(nma.getEigenvalues().size() == 6);
  // Isolated nodes 3, 4, 5 and the connected component each give a zero
  // eigenvalue -> at least 4 zero eigenvalues
  CHECK(nma.getNrOfSolubleClusters() >= 4);
}

#ifdef CEREALIZABLE
TEST_CASE(
  "NormalModeAnalyzer serialization preserves assembledConnectivityMatrix",
  "[analysis][NormalModeAnalyzer][serialization]")
{
  pc::NormalModeAnalyzer nma({ 0, 1, 2 }, { 1, 2, 0 });
  nma.computeAllEigenvalues();
  std::string serialized = pylimer_tools::utils::serializeToString(nma);
  pc::NormalModeAnalyzer restored =
    pc::NormalModeAnalyzer::fromString(serialized);

  CHECK(restored.getMatrixSize() == 3);
  CHECK(restored.getAssembledConnectivityMatrix().nonZeros() ==
        nma.getAssembledConnectivityMatrix().nonZeros());
  CHECK(restored.getEigenvalues().isApprox(nma.getEigenvalues()));
}
#endif

TEST_CASE(
  "NormalModeAnalyzer component decomposition produces valid eigenvectors",
  "[analysis][NormalModeAnalyzer]")
{
  // Two disconnected components: triangle (0,1,2) and line (3,4)
  std::vector<size_t> from = { 0, 1, 2, 3 };
  std::vector<size_t> to = { 1, 2, 0, 4 };
  pc::NormalModeAnalyzer nma(from, to);
  nma.computeAllEigenvalues(pc::NormalModeSolver::DENSE_EIGEN, true, true);

  Eigen::VectorXd evals = nma.getEigenvalues();
  Eigen::MatrixXd evecs = nma.getEigenvectors();
  REQUIRE(evals.size() == 5);
  REQUIRE(evecs.rows() == 5);
  REQUIRE(evecs.cols() == 5);

  // Orthonormality check: V^T * V = I
  Eigen::MatrixXd identity = evecs.transpose() * evecs;
  CHECK(identity.isApprox(Eigen::MatrixXd::Identity(5, 5), 1e-10));

  // Eigenvector check: A * v_i = lambda_i * v_i
  Eigen::MatrixXd A = nma.getAssembledConnectivityMatrix();
  for (Eigen::Index i = 0; i < 5; ++i) {
    Eigen::VectorXd Av = A * evecs.col(i);
    Eigen::VectorXd lv = evals[i] * evecs.col(i);
    CHECK((Av - lv).norm() < 1e-8);
  }
}

TEST_CASE("NormalModeAnalyzer computes correct eigenvalues",
          "[analysis][NormalModeAnalyzer]")
{
  pc::NormalModeAnalyzer normalModeAnalyzer =
    pc::NormalModeAnalyzer({ 0, 1, 2 }, { 1, 2, 0 });

  Eigen::Matrix3d assembledMatrix =
    normalModeAnalyzer.getAssembledConnectivityMatrix();
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      CHECK(assembledMatrix(i, j) == (i == j ? 2 : -1));
    }
  }

  normalModeAnalyzer.computeAllEigenvalues(true);
  Eigen::VectorXd eigenvalues = normalModeAnalyzer.getEigenvalues();

  CHECK_THAT(eigenvalues[0], Catch::Matchers::WithinAbs(0., 1e-9));
  // could be 3. or 0.,
  // CHECK_THAT(eigenvalues[1], Catch::Matchers::WithinAbs(0., 1e-9));
  CHECK_THAT(eigenvalues[2], Catch::Matchers::WithinRel(3.));

  Eigen::MatrixXd eigenvectors = normalModeAnalyzer.getEigenvectors();
  // Verify eigenvectors are orthonormal
  Eigen::MatrixXd shouldBeIdentity = eigenvectors.transpose() * eigenvectors;
  CHECK(shouldBeIdentity.isApprox(Eigen::Matrix3d::Identity(), 1e-10));

  // Verify each column is a valid eigenvector of A by checking the Rayleigh
  // quotient matches one of the known eigenvalues (0 or 3)
  for (int i = 0; i < 3; ++i) {
    Eigen::VectorXd v = eigenvectors.col(i);
    double rayleighQuotient = v.dot(assembledMatrix * v);
    Eigen::VectorXd Av = assembledMatrix * v;
    // A*v should be parallel to v (i.e. A*v = rayleighQuotient * v)
    CHECK((Av - rayleighQuotient * v).norm() < 1e-9);
    // Rayleigh quotient should match one of the eigenvalues
    bool matchesAnEigenvalue = false;
    for (int j = 0; j < eigenvalues.size(); ++j) {
      if (std::abs(rayleighQuotient - eigenvalues[j]) < 1e-9) {
        matchesAnEigenvalue = true;
        break;
      }
    }
    CHECK(matchesAnEigenvalue);
  }

  Eigen::Array3d t;
  t << 0., 1., 10.;
  Eigen::ArrayXd stressCorrelation =
    normalModeAnalyzer.evaluateStressAutocorrelation(t);
  CHECK(stressCorrelation.size() == 3);
  CHECK(stressCorrelation[0] >= stressCorrelation[1]);
  CHECK(stressCorrelation[1] >= stressCorrelation[2]);
  Eigen::Array3d omega;
  omega << 0., 1., 2.;
  Eigen::ArrayXd lossModulus = normalModeAnalyzer.evaluateLossModulus(omega);
  CHECK(lossModulus.size() == 3);
  CHECK(lossModulus[0] <= lossModulus[1]);
  CHECK(lossModulus[1] <= lossModulus[2]);
  Eigen::ArrayXd storageModulus =
    normalModeAnalyzer.evaluateStorageModulus(omega);
  CHECK(storageModulus.size() == 3);
  CHECK(storageModulus[0] <= storageModulus[1]);
  CHECK(storageModulus[1] <= storageModulus[2]);
}

TEST_CASE("NormalModeAnalyzer can compute some eigenvalues",
          "[analysis][NormalModeAnalyzer]")
{
  pc::NormalModeAnalyzer normalModeAnalyzer =
    pc::NormalModeAnalyzer({ 0, 1, 2 }, { 1, 2, 0 });

  normalModeAnalyzer.findSparseEigenvalues(1, true);
  Eigen::VectorXd eigenvalues = normalModeAnalyzer.getEigenvalues();

  CHECK(eigenvalues.size() == 1);

  CHECK_THAT(eigenvalues[0], Catch::Matchers::WithinAbs(0., 1e-10));
}

TEST_CASE("NormalModeAnalyzer solvers give consistent eigenvalues",
          "[analysis][NormalModeAnalyzer]")
{
  // A test network with 40 nodes:
  // - a connected cycle/cluster
  // - a separate linear chain
  // - isolated beads
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;

  // 1. Cluster on nodes 0..24
  for (size_t i = 0; i < 25; ++i) {
    springFrom.push_back(i);
    springTo.push_back((i + 1) % 25);
    size_t cross = (i * 7 + 3) % 25;
    if (cross != i && cross != (i + 1) % 25) {
      springFrom.push_back(i);
      springTo.push_back(cross);
    }
  }

  // 2. Linear chain on nodes 25..34
  for (size_t i = 25; i < 34; ++i) {
    springFrom.push_back(i);
    springTo.push_back(i + 1);
  }

  // Nodes 35..39: isolated beads (implicit since maxIdx = 39 if we connect
  // 38-39 or add edge)
  springFrom.push_back(38);
  springTo.push_back(39);

  pc::NormalModeAnalyzer reference(springFrom, springTo);
  reference.computeAllEigenvalues(
    pc::NormalModeSolver::DENSE_LAPACK, false, false);
  Eigen::VectorXd expected = reference.getEigenvalues();

  for (auto solver : pc::allNormalModeSolvers) {
    DYNAMIC_SECTION("Solver: "
                    << pc::NormalModeSolverNames[static_cast<int>(solver)]
                    << " (without decomposition)")
    {
      pc::NormalModeAnalyzer nma(springFrom, springTo);
      nma.computeAllEigenvalues(solver, false, false);
      Eigen::VectorXd actual = nma.getEigenvalues();

      REQUIRE(actual.size() == expected.size());
      for (Eigen::Index i = 0; i < actual.size(); ++i) {
        CHECK_THAT(actual[i], Catch::Matchers::WithinAbs(expected[i], 1e-4));
      }
    }

    DYNAMIC_SECTION("Solver: "
                    << pc::NormalModeSolverNames[static_cast<int>(solver)]
                    << " (with decomposition)")
    {
      pc::NormalModeAnalyzer nma(springFrom, springTo);
      nma.computeAllEigenvalues(solver, false, true);
      Eigen::VectorXd actual = nma.getEigenvalues();

      REQUIRE(actual.size() == expected.size());
      for (Eigen::Index i = 0; i < actual.size(); ++i) {
        CHECK_THAT(actual[i], Catch::Matchers::WithinAbs(expected[i], 1e-4));
      }
    }
  }
}

TEST_CASE("NormalModeAnalyzer KPM computes continuous spectra and moduli",
          "[analysis][NormalModeAnalyzer]")
{
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  size_t n = 50;
  for (size_t i = 0; i < n; ++i) {
    springFrom.push_back(i);
    springTo.push_back((i + 1) % n);
    size_t cross = (i * 11 + 5) % n;
    if (cross != i) {
      springFrom.push_back(i);
      springTo.push_back(cross);
    }
  }

  pc::NormalModeAnalyzer nma(springFrom, springTo);
  nma.computeAllEigenvalues(pc::NormalModeSolver::AUTO);

  // 1. Direct Chebyshev moments test
  Eigen::VectorXd moments = nma.computeChebyshevMomentsKPM(256, 30);
  REQUIRE(moments.size() == 256);
  // mu[0] is the normalized trace of I / n = 1.0 (with Jackson damping g0
  // = 1.0)
  CHECK_THAT(moments[0], Catch::Matchers::WithinAbs(1.0, 0.15));

  Eigen::ArrayXd omega(5);
  omega << 0.1, 0.5, 1.0, 2.0, 5.0;

  Eigen::ArrayXd gPrimeExact = nma.evaluateStorageModulus(omega);
  Eigen::ArrayXd gDoublePrimeExact = nma.evaluateLossModulus(omega);

  Eigen::ArrayXd gPrimeKPM = nma.evaluateStorageModulusKPM(omega, 512, 30);
  Eigen::ArrayXd gDoublePrimeKPM = nma.evaluateLossModulusKPM(omega, 512, 30);

  CHECK(gPrimeKPM.size() == 5);
  CHECK(gDoublePrimeKPM.size() == 5);

  // KPM should track the general magnitude and monotonic trend
  CHECK(gPrimeKPM[0] <= gPrimeKPM[4]);
  CHECK(gDoublePrimeKPM[0] <= gDoublePrimeKPM[3]);

  // Moduli should be within comparable magnitude of exact solution
  for (Eigen::Index i = 1; i < 4; ++i) {
    CHECK(gPrimeKPM[i] > 0.0);
    CHECK(gDoublePrimeKPM[i] > 0.0);
    // Relative order of magnitude check
    CHECK(gPrimeKPM[i] > 0.3 * gPrimeExact[i]);
    CHECK(gDoublePrimeKPM[i] > 0.3 * gDoublePrimeExact[i]);
  }

  // 2. Stress autocorrelation test
  Eigen::ArrayXd t(3);
  t << 0.0, 0.5, 2.0;
  Eigen::ArrayXd gtKPM = nma.evaluateStressAutocorrelationKPM(t, 256, 20);
  CHECK(gtKPM.size() == 3);
  CHECK(gtKPM[0] >= gtKPM[1]);
  CHECK(gtKPM[1] >= gtKPM[2]);

  // 3. Density of states
  Eigen::ArrayXd grid(10);
  grid.setLinSpaced(10, 0.5, 5.0);
  Eigen::ArrayXd dos = nma.computeDensityOfStatesKPM(grid, 256, 10);
  CHECK(dos.size() == 10);
  for (Eigen::Index i = 0; i < dos.size(); ++i) {
    CHECK(dos[i] >= 0.0);
  }
}

TEST_CASE("NormalModeAnalyzer edge cases and multigraphs",
          "[analysis][NormalModeAnalyzer]")
{
  // Empty springs throw invalid_argument
  REQUIRE_THROWS_AS(pc::NormalModeAnalyzer({}, {}), std::invalid_argument);

  // Mismatched vector sizes throw invalid_argument
  REQUIRE_THROWS_AS(pc::NormalModeAnalyzer({ 0, 1 }, { 1 }),
                    std::invalid_argument);

  // Multigraph: 2 beads with 2 parallel springs
  // Diagonal is (2, 2), off-diagonal (-2, -2). Eigenvalues are 0.0 and 4.0.
  pc::NormalModeAnalyzer dimerParallel({ 0, 0 }, { 1, 1 });
  dimerParallel.computeAllEigenvalues(pc::NormalModeSolver::AUTO, false, true);
  Eigen::VectorXd dimerEvals = dimerParallel.getEigenvalues();
  REQUIRE(dimerEvals.size() == 2);
  CHECK_THAT(dimerEvals[0], Catch::Matchers::WithinAbs(0.0, 1e-9));
  CHECK_THAT(dimerEvals[1], Catch::Matchers::WithinAbs(4.0, 1e-9));

  // Decomposed component with >30 nodes to exercise SPARSE_SLICING inside
  // decomposed branch
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  size_t clusterSize = 35;
  for (size_t i = 0; i < clusterSize; ++i) {
    springFrom.push_back(i);
    springTo.push_back((i + 1) % clusterSize);
    size_t cross = (i * 7 + 3) % clusterSize;
    if (cross != i && cross != (i + 1) % clusterSize) {
      springFrom.push_back(i);
      springTo.push_back(cross);
    }
  }
  // Linear chain on 35..44
  for (size_t i = 35; i < 44; ++i) {
    springFrom.push_back(i);
    springTo.push_back(i + 1);
  }

  pc::NormalModeAnalyzer largeSystem(springFrom, springTo);
  largeSystem.computeAllEigenvalues(
    pc::NormalModeSolver::DENSE_LAPACK, false, false);
  Eigen::VectorXd refEvals = largeSystem.getEigenvalues();

  pc::NormalModeAnalyzer slicingSystem(springFrom, springTo);
  slicingSystem.computeAllEigenvalues(
    pc::NormalModeSolver::SPARSE_SLICING, false, true);
  Eigen::VectorXd slicingEvals = slicingSystem.getEigenvalues();

  REQUIRE(slicingEvals.size() == refEvals.size());
  for (Eigen::Index i = 0; i < slicingEvals.size(); ++i) {
    CHECK_THAT(slicingEvals[i], Catch::Matchers::WithinAbs(refEvals[i], 1e-4));
  }
}

TEST_CASE("NormalModeAnalyzer benchmark", "[!benchmark][NormalModeAnalyzer]")
{
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  size_t n = 1000;
  for (size_t i = 0; i < n; ++i) {
    springFrom.push_back(i);
    springTo.push_back((i + 1) % n);
    size_t cross = (i * 37 + 17) % n;
    if (cross != i) {
      springFrom.push_back(i);
      springTo.push_back(cross);
    }
  }

  BENCHMARK("computeAllEigenvalues (eigenvalues only)")
  {
    pc::NormalModeAnalyzer nma(springFrom, springTo);
    nma.computeAllEigenvalues(false);
    return nma.getEigenvalues().size();
  };

  BENCHMARK("computeAllEigenvalues (with eigenvectors)")
  {
    pc::NormalModeAnalyzer nma(springFrom, springTo);
    nma.computeAllEigenvalues(true);
    return nma.getEigenvalues().size();
  };
}

TEST_CASE("NormalModeAnalyzer Metal GPU solver computes accurate eigenvalues "
          "and eigenvectors",
          "[analysis][NormalModeAnalyzer][Metal]")
{
  // Test 1: Simple 3-bead ring
  {
    pc::NormalModeAnalyzer nmaCPU({ 0, 1, 2 }, { 1, 2, 0 });
    nmaCPU.computeAllEigenvalues(pc::NormalModeSolver::DENSE_LAPACK, true);

    pc::NormalModeAnalyzer nmaGPU({ 0, 1, 2 }, { 1, 2, 0 });
    nmaGPU.computeAllEigenvaluesGPU(true);

    REQUIRE(nmaGPU.getEigenvalues().size() == 3);
    for (Eigen::Index i = 0; i < 3; ++i) {
      CHECK_THAT(nmaGPU.getEigenvalues()[i],
                 Catch::Matchers::WithinAbs(nmaCPU.getEigenvalues()[i], 1e-4));
    }

    // Check eigenvector orthonormality: V^T * V = I
    Eigen::MatrixXd V = nmaGPU.getEigenvectors();
    Eigen::MatrixXd Id = V.transpose() * V;
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        CHECK_THAT(Id(r, c),
                   Catch::Matchers::WithinAbs(r == c ? 1.0 : 0.0, 1e-3));
      }
    }
  }

  // Test 2: 20-bead ring with crosslinks
  {
    std::vector<size_t> springFrom;
    std::vector<size_t> springTo;
    size_t n = 20;
    for (size_t i = 0; i < n; ++i) {
      springFrom.push_back(i);
      springTo.push_back((i + 1) % n);
      if (i % 5 == 0) {
        springFrom.push_back(i);
        springTo.push_back((i + 10) % n);
      }
    }

    pc::NormalModeAnalyzer nmaCPU(springFrom, springTo);
    nmaCPU.computeAllEigenvalues(pc::NormalModeSolver::DENSE_LAPACK, false);

    pc::NormalModeAnalyzer nmaGPU(springFrom, springTo);
    nmaGPU.computeAllEigenvalues(pc::NormalModeSolver::METAL, false);

    REQUIRE(nmaGPU.getEigenvalues().size() == static_cast<Eigen::Index>(n));
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(n); ++i) {
      CHECK_THAT(nmaGPU.getEigenvalues()[i],
                 Catch::Matchers::WithinAbs(nmaCPU.getEigenvalues()[i], 1e-3));
    }
  }
}

TEST_CASE("NormalModeAnalyzer Metal KPM computes accurate Chebyshev moments",
          "[analysis][NormalModeAnalyzer][Metal][KPM]")
{
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  size_t n = 100;
  for (size_t i = 0; i < n; ++i) {
    springFrom.push_back(i);
    springTo.push_back((i + 1) % n);
    size_t cross = (i * 7 + 3) % n;
    if (cross != i) {
      springFrom.push_back(i);
      springTo.push_back(cross);
    }
  }

  pc::NormalModeAnalyzer nma(springFrom, springTo);
  size_t numMoments = 128;
  size_t numRandomVectors = 10;

  Eigen::VectorXd cpuMoments =
    nma.computeChebyshevMomentsKPM(numMoments, numRandomVectors);
  Eigen::VectorXd gpuMomentsB4 =
    nma.computeChebyshevMomentsKPMMetal(numMoments, numRandomVectors, 4);
  Eigen::VectorXd gpuMomentsB8 =
    nma.computeChebyshevMomentsKPMMetal(numMoments, numRandomVectors, 8);
  Eigen::VectorXd gpuMomentsB1 =
    nma.computeChebyshevMomentsKPMMetal(numMoments, numRandomVectors, 1);

  REQUIRE(gpuMomentsB4.size() == static_cast<Eigen::Index>(numMoments));
  REQUIRE(gpuMomentsB8.size() == static_cast<Eigen::Index>(numMoments));
  REQUIRE(gpuMomentsB1.size() == static_cast<Eigen::Index>(numMoments));

  for (size_t m = 0; m < numMoments; ++m) {
    CHECK_THAT(gpuMomentsB4[static_cast<Eigen::Index>(m)],
               Catch::Matchers::WithinAbs(
                 cpuMoments[static_cast<Eigen::Index>(m)], 1e-3));
    CHECK_THAT(gpuMomentsB8[static_cast<Eigen::Index>(m)],
               Catch::Matchers::WithinAbs(
                 cpuMoments[static_cast<Eigen::Index>(m)], 1e-3));
    CHECK_THAT(gpuMomentsB1[static_cast<Eigen::Index>(m)],
               Catch::Matchers::WithinAbs(
                 cpuMoments[static_cast<Eigen::Index>(m)], 1e-3));
    // Verify agreement between batched and sequential implementations
    CHECK_THAT(gpuMomentsB4[static_cast<Eigen::Index>(m)],
               Catch::Matchers::WithinAbs(
                 gpuMomentsB1[static_cast<Eigen::Index>(m)], 1e-4));
    CHECK_THAT(gpuMomentsB8[static_cast<Eigen::Index>(m)],
               Catch::Matchers::WithinAbs(
                 gpuMomentsB1[static_cast<Eigen::Index>(m)], 1e-4));
  }
}

TEST_CASE("NormalModeAnalyzer Metal vs CPU KPM benchmark",
          "[!benchmark][NormalModeAnalyzer][Metal]")
{
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  size_t n = 2000;
  for (size_t i = 0; i < n; ++i) {
    springFrom.push_back(i);
    springTo.push_back((i + 1) % n);
    size_t cross = (i * 37 + 17) % n;
    if (cross != i) {
      springFrom.push_back(i);
      springTo.push_back(cross);
    }
  }

  pc::NormalModeAnalyzer nma(springFrom, springTo);

  BENCHMARK(
    "Metal GPU KPM Sequential (B=1, 2000 nodes, 256 moments, 16 vectors)")
  {
    return nma.computeChebyshevMomentsKPMMetal(256, 16, 1).size();
  };

  BENCHMARK("Metal GPU KPM Multi-RHS Block-SpMM (B=4, 2000 nodes, 256 moments, "
            "16 vectors)")
  {
    return nma.computeChebyshevMomentsKPMMetal(256, 16, 4).size();
  };

  BENCHMARK("Metal GPU KPM Multi-RHS Block-SpMM (B=8, 2000 nodes, 256 moments, "
            "16 vectors)")
  {
    return nma.computeChebyshevMomentsKPMMetal(256, 16, 8).size();
  };
}
