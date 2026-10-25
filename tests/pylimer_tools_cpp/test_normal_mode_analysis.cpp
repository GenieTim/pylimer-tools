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
  REQUIRE_NOTHROW(normalModeAnalyzer.computeAllEigenvalues());
  Eigen::VectorXd eigenvalues = normalModeAnalyzer.getEigenvalues();
  CHECK(eigenvalues.size() == 5);
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
