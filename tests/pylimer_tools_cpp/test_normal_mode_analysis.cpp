#include "../../src/pylimer_tools_cpp/calc/NormalModeAnalyzer.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#ifdef USE_METAL
#include "../../src/pylimer_tools_cpp/calc/MetalNormalModeSolver.h"
#endif

#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace pc = pylimer_tools::calc;
namespace pe = pylimer_tools::entities;

TEST_CASE("NormalModeAnalyzer does not crash", "[analysis][NormalModeAnalyzer]")
{
  pc::NormalModeAnalyzer normalModeAnalyzer({ 1, 2, 3 }, { 2, 3, 4 });

  CHECK(normalModeAnalyzer.getSize() == 5);
  CHECK(normalModeAnalyzer.getNumSolubleClusters() == 2);

  pc::NormalModeSpectrum spectrum = normalModeAnalyzer.solveExact();
  CHECK(spectrum.isComplete());
  CHECK(spectrum.getNumModes() == 5);
  CHECK(spectrum.getNumZeroModes() == 2);
  CHECK_FALSE(spectrum.getEigenvectors().has_value());

  Eigen::ArrayXd testGrid = Eigen::ArrayXd::LinSpaced(5, 0.1, 1.0);
  CHECK(spectrum.stressAutocorrelation(testGrid).allFinite());
  CHECK(spectrum.storageModulus(testGrid).allFinite());
  CHECK(spectrum.lossModulus(testGrid).allFinite());
}

TEST_CASE("NormalModeAnalyzer constructor supports explicit nrOfNodes for "
          "isolated vertices",
          "[analysis][NormalModeAnalyzer]")
{
  pc::NormalModeAnalyzer nma({ 0, 1 }, { 1, 2 }, 6);
  CHECK(nma.getSize() == 6);
  CHECK(nma.getNumSolubleClusters() == 4);

  pc::NormalModeSpectrum spectrum = nma.solveExact();
  CHECK(spectrum.getNumModes() == 6);
  CHECK(spectrum.getNumZeroModes() == 4);
}

TEST_CASE("NormalModeAnalyzer fromMatrix, fromUniverse, and weighted springs",
          "[analysis][NormalModeAnalyzer]")
{
  // Weighted springs: 2-node system with spring constant 2.5 -> evals {0, 5.0}
  pc::NormalModeAnalyzer weighted({ 0 }, { 1 }, 2, std::vector<double>{ 2.5 });
  pc::NormalModeSpectrum wSpec = weighted.solveExact();
  REQUIRE(wSpec.getNumModes() == 2);
  CHECK_THAT(wSpec.getEigenvalues()[0], Catch::Matchers::WithinAbs(0.0, 1e-9));
  CHECK_THAT(wSpec.getEigenvalues()[1], Catch::Matchers::WithinAbs(5.0, 1e-9));

  // fromMatrix round-trip
  pc::NormalModeAnalyzer fromMat =
    pc::NormalModeAnalyzer::fromMatrix(weighted.getMatrix());
  CHECK(fromMat.getSize() == 2);
  CHECK(fromMat.solveExact().getEigenvalues().isApprox(wSpec.getEigenvalues()));

  // Non-square matrix throws
  Eigen::SparseMatrix<double> rect(2, 3);
  REQUIRE_THROWS_AS(pc::NormalModeAnalyzer::fromMatrix(rect),
                    std::invalid_argument);

  // fromUniverse
  pe::Universe u(10.0, 10.0, 10.0);
  u.addAtoms({ 1, 2, 3 },
             { 1, 1, 1 },
             { 0., 1., 2. },
             { 0., 0., 0. },
             { 0., 0., 0. },
             { 0, 0, 0 },
             { 0, 0, 0 },
             { 0, 0, 0 });
  u.addBonds({ 1, 2 }, { 2, 3 });
  pc::NormalModeAnalyzer fromU = pc::NormalModeAnalyzer::fromUniverse(u);
  CHECK(fromU.getSize() == 3);
  CHECK(fromU.getNumSolubleClusters() == 1);
  pc::NormalModeSpectrum uSpec = fromU.solveExact();
  CHECK(uSpec.getNumModes() == 3);
  CHECK(uSpec.getNumZeroModes() == 1);
}

#ifdef CEREALIZABLE
TEST_CASE("NormalModeAnalyzer and NormalModeSpectrum serialization",
          "[analysis][NormalModeAnalyzer][serialization]")
{
  pc::NormalModeAnalyzer nma({ 0, 1, 2 }, { 1, 2, 0 });
  pc::NormalModeSpectrum spec =
    nma.solveExact(pc::NormalModeSolver::AUTO, true);

  std::string serializedNma = pylimer_tools::utils::serializeToString(nma);
  pc::NormalModeAnalyzer restoredNma =
    pc::NormalModeAnalyzer::fromString(serializedNma);
  CHECK(restoredNma.getSize() == 3);
  CHECK(restoredNma.getMatrix().nonZeros() == nma.getMatrix().nonZeros());

  std::string serializedSpec = pylimer_tools::utils::serializeToString(spec);
  pc::NormalModeSpectrum restoredSpec =
    pc::NormalModeSpectrum::fromString(serializedSpec);
  CHECK(restoredSpec.getEigenvalues().isApprox(spec.getEigenvalues()));
  CHECK(restoredSpec.getWeights().isApprox(spec.getWeights()));
  REQUIRE(restoredSpec.getEigenvectors().has_value());
  CHECK(restoredSpec.getEigenvectors()->isApprox(*spec.getEigenvectors()));
  CHECK(restoredSpec.isComplete() == spec.isComplete());
}
#endif

TEST_CASE(
  "NormalModeAnalyzer component decomposition produces valid eigenvectors",
  "[analysis][NormalModeAnalyzer]")
{
  std::vector<size_t> from = { 0, 1, 2, 3 };
  std::vector<size_t> to = { 1, 2, 0, 4 };
  pc::NormalModeAnalyzer nma(from, to);
  pc::NormalModeSpectrum spec =
    nma.solveExact(pc::NormalModeSolver::DENSE_EIGEN, true, true);

  const Eigen::VectorXd& evals = spec.getEigenvalues();
  REQUIRE(spec.getEigenvectors().has_value());
  const Eigen::MatrixXd evecs = *spec.getEigenvectors();
  REQUIRE(evals.size() == 5);
  REQUIRE(evecs.rows() == 5);
  REQUIRE(evecs.cols() == 5);

  Eigen::MatrixXd identity = evecs.transpose() * evecs;
  CHECK(identity.isApprox(Eigen::MatrixXd::Identity(5, 5), 1e-10));

  Eigen::MatrixXd A = nma.getMatrix();
  for (Eigen::Index i = 0; i < 5; ++i) {
    Eigen::VectorXd Av = A * evecs.col(i);
    Eigen::VectorXd lv = evals[i] * evecs.col(i);
    CHECK((Av - lv).norm() < 1e-8);
  }
}

TEST_CASE("NormalModeAnalyzer computes correct eigenvalues and observables",
          "[analysis][NormalModeAnalyzer]")
{
  pc::NormalModeAnalyzer normalModeAnalyzer({ 0, 1, 2 }, { 1, 2, 0 });

  Eigen::Matrix3d assembledMatrix = normalModeAnalyzer.getMatrix();
  for (Eigen::Index i = 0; i < 3; ++i) {
    for (Eigen::Index j = 0; j < 3; ++j) {
      CHECK(assembledMatrix(i, j) == (i == j ? 2 : -1));
    }
  }

  pc::NormalModeSpectrum spec =
    normalModeAnalyzer.solveExact(pc::NormalModeSolver::AUTO, true);
  const Eigen::VectorXd& eigenvalues = spec.getEigenvalues();

  CHECK_THAT(eigenvalues[0], Catch::Matchers::WithinAbs(0., 1e-9));
  CHECK_THAT(eigenvalues[2], Catch::Matchers::WithinRel(3.));

  REQUIRE(spec.getEigenvectors().has_value());
  const Eigen::MatrixXd eigenvectors = *spec.getEigenvectors();
  Eigen::MatrixXd shouldBeIdentity = eigenvectors.transpose() * eigenvectors;
  CHECK(shouldBeIdentity.isApprox(Eigen::Matrix3d::Identity(), 1e-10));

  for (int i = 0; i < 3; ++i) {
    Eigen::VectorXd v = eigenvectors.col(i);
    double rayleighQuotient = v.dot(assembledMatrix * v);
    Eigen::VectorXd Av = assembledMatrix * v;
    CHECK((Av - rayleighQuotient * v).norm() < 1e-9);
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
  Eigen::ArrayXd stressCorrelation = spec.stressAutocorrelation(t);
  CHECK(stressCorrelation.size() == 3);
  CHECK(stressCorrelation[0] >= stressCorrelation[1]);
  CHECK(stressCorrelation[1] >= stressCorrelation[2]);

  Eigen::Array3d omega;
  omega << 0., 1., 2.;
  Eigen::ArrayXd lossModulus = spec.lossModulus(omega);
  CHECK(lossModulus.size() == 3);
  CHECK(lossModulus[0] <= lossModulus[1]);
  CHECK(lossModulus[1] <= lossModulus[2]);
  Eigen::ArrayXd storageModulus = spec.storageModulus(omega);
  CHECK(storageModulus.size() == 3);
  CHECK(storageModulus[0] <= storageModulus[1]);
  CHECK(storageModulus[1] <= storageModulus[2]);

  auto [gpPair, gppPair] = spec.moduli(omega);
  CHECK(gpPair.isApprox(storageModulus));
  CHECK(gppPair.isApprox(lossModulus));
}

TEST_CASE("NormalModeAnalyzer can compute some eigenvalues (numModes)",
          "[analysis][NormalModeAnalyzer]")
{
  pc::NormalModeAnalyzer normalModeAnalyzer({ 0, 1, 2 }, { 1, 2, 0 });

  pc::NormalModeSpectrum spec =
    normalModeAnalyzer.solveExact(pc::NormalModeSolver::AUTO, true, true, 1);
  CHECK_FALSE(spec.isComplete());
  CHECK(spec.getEigenvalues().size() == 1);
  CHECK_THAT(spec.getEigenvalues()[0], Catch::Matchers::WithinAbs(0., 1e-10));

  pc::NormalModeAnalyzer chainNma({ 0, 1, 2, 3, 4, 5, 6, 7 },
                                  { 1, 2, 3, 4, 5, 6, 7, 8 });
  pc::NormalModeSpectrum sparseSpec =
    chainNma.solveExact(pc::NormalModeSolver::AUTO, true, true, 4);
  CHECK_FALSE(sparseSpec.isComplete());
  const Eigen::VectorXd& sparseEvals = sparseSpec.getEigenvalues();
  REQUIRE(sparseSpec.getEigenvectors().has_value());
  const Eigen::MatrixXd sparseEvecs = *sparseSpec.getEigenvectors();
  REQUIRE(sparseEvals.size() == 4);
  REQUIRE(sparseEvecs.rows() == 9);
  REQUIRE(sparseEvecs.cols() == 4);
  for (Eigen::Index i = 0; i + 1 < sparseEvals.size(); ++i) {
    CHECK(sparseEvals[i] <= sparseEvals[i + 1]);
  }
  const Eigen::MatrixXd A = chainNma.getMatrix();
  for (Eigen::Index i = 0; i < sparseEvals.size(); ++i) {
    CHECK(
      (A * sparseEvecs.col(i) - sparseEvals[i] * sparseEvecs.col(i)).norm() <
      1e-8);
  }
}

TEST_CASE("NormalModeAnalyzer ignores self-loops and validates spectrum inputs",
          "[analysis][NormalModeAnalyzer]")
{
  pc::NormalModeAnalyzer cleanChain({ 0, 1, 2 }, { 1, 2, 3 });
  pc::NormalModeAnalyzer loopChain({ 0, 0, 1, 2, 2 }, { 0, 1, 2, 2, 3 });

  CHECK(loopChain.getMatrix().nonZeros() == cleanChain.getMatrix().nonZeros());

  const Eigen::VectorXd expectedEvals =
    cleanChain.solveExact(pc::NormalModeSolver::DENSE_EIGEN, false, false)
      .getEigenvalues();

  CHECK(loopChain.solveExact(pc::NormalModeSolver::DENSE_EIGEN, false, false)
          .getEigenvalues()
          .isApprox(expectedEvals, 1e-9));
  CHECK(loopChain.solveExact(pc::NormalModeSolver::DENSE_EIGEN, false, true)
          .getEigenvalues()
          .isApprox(expectedEvals, 1e-9));

  // Validate NormalModeSpectrum dimension checks
  REQUIRE_THROWS_AS(
    pc::NormalModeSpectrum(Eigen::VectorXd::Zero(3), Eigen::VectorXd::Ones(2)),
    std::invalid_argument);
  REQUIRE_THROWS_AS(pc::NormalModeSpectrum(Eigen::VectorXd::Zero(3),
                                           std::nullopt,
                                           Eigen::MatrixXd::Zero(3, 4)),
                    std::invalid_argument);

  // Validate backend and batchSize checks
  pc::SLQConfig badSlq;
  badSlq.backend = pc::ComputeBackend::METAL;
  REQUIRE_THROWS_AS(loopChain.solveApprox(badSlq), std::invalid_argument);

  pc::HybridConfig badHybrid;
  badHybrid.backend = pc::ComputeBackend::CUDA;
  REQUIRE_THROWS_AS(loopChain.solveApprox(badHybrid), std::invalid_argument);

  pc::KPMConfig badKpmBackend;
  badKpmBackend.backend = pc::ComputeBackend::CUDA;
  REQUIRE_THROWS_AS(loopChain.chebyshevMoments(badKpmBackend),
                    std::invalid_argument);

  pc::KPMConfig badKpmBatch;
  badKpmBatch.batchSize = 0;
  REQUIRE_THROWS_AS(loopChain.chebyshevMoments(badKpmBatch),
                    std::invalid_argument);
}

TEST_CASE("NormalModeAnalyzer solvers give consistent eigenvalues",
          "[analysis][NormalModeAnalyzer]")
{
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;

  for (size_t i = 0; i < 25; ++i) {
    springFrom.push_back(i);
    springTo.push_back((i + 1) % 25);
    size_t cross = (i * 7 + 3) % 25;
    if (cross != i && cross != (i + 1) % 25) {
      springFrom.push_back(i);
      springTo.push_back(cross);
    }
  }

  for (size_t i = 25; i < 34; ++i) {
    springFrom.push_back(i);
    springTo.push_back(i + 1);
  }

  springFrom.push_back(38);
  springTo.push_back(39);

  pc::NormalModeAnalyzer nma(springFrom, springTo);
  Eigen::VectorXd expected =
    nma.solveExact(pc::NormalModeSolver::DENSE_LAPACK, false, false)
      .getEigenvalues();

  for (const pc::NormalModeSolver solver : pc::allNormalModeSolvers) {
    DYNAMIC_SECTION("Solver: "
                    << pc::NormalModeSolverNames[static_cast<int>(solver)]
                    << " (without decomposition)")
    {
      Eigen::VectorXd actual =
        nma.solveExact(solver, false, false).getEigenvalues();
      REQUIRE(actual.size() == expected.size());
      for (Eigen::Index i = 0; i < actual.size(); ++i) {
        CHECK_THAT(actual[i], Catch::Matchers::WithinAbs(expected[i], 1e-4));
      }
    }

    DYNAMIC_SECTION("Solver: "
                    << pc::NormalModeSolverNames[static_cast<int>(solver)]
                    << " (with decomposition)")
    {
      Eigen::VectorXd actual =
        nma.solveExact(solver, false, true).getEigenvalues();
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
  pc::NormalModeSpectrum exactSpec = nma.solveExact();

  pc::KPMConfig momentsCfg;
  momentsCfg.numMoments = 256;
  momentsCfg.numRandomVectors = 30;
  Eigen::VectorXd moments = nma.chebyshevMoments(momentsCfg);
  REQUIRE(moments.size() == 256);
  CHECK_THAT(moments[0], Catch::Matchers::WithinAbs(1.0, 0.15));

  Eigen::ArrayXd omega(5);
  omega << 0.1, 0.5, 1.0, 2.0, 5.0;

  auto [gPrimeExact, gDoublePrimeExact] = exactSpec.moduli(omega);

  pc::KPMConfig kpmCfg;
  kpmCfg.numMoments = 512;
  kpmCfg.numRandomVectors = 30;
  kpmCfg.numQuadraturePoints = 2048;
  pc::NormalModeSpectrum kpmSpec = nma.solveApprox(kpmCfg);
  CHECK_FALSE(kpmSpec.isComplete());

  auto [gPrimeKPM, gDoublePrimeKPM] = kpmSpec.moduli(omega);

  CHECK(gPrimeKPM.size() == 5);
  CHECK(gDoublePrimeKPM.size() == 5);
  CHECK(gPrimeKPM[0] <= gPrimeKPM[4]);
  CHECK(gDoublePrimeKPM[0] <= gDoublePrimeKPM[3]);

  for (Eigen::Index i = 1; i < 4; ++i) {
    CHECK(gPrimeKPM[i] > 0.0);
    CHECK(gDoublePrimeKPM[i] > 0.0);
    CHECK(gPrimeKPM[i] > 0.3 * gPrimeExact[i]);
    CHECK(gDoublePrimeKPM[i] > 0.3 * gDoublePrimeExact[i]);
  }

  Eigen::ArrayXd t(3);
  t << 0.0, 0.5, 2.0;
  Eigen::ArrayXd gtKPM = kpmSpec.stressAutocorrelation(t);
  CHECK(gtKPM.size() == 3);
  CHECK(gtKPM[0] >= gtKPM[1]);
  CHECK(gtKPM[1] >= gtKPM[2]);

  Eigen::ArrayXd grid(10);
  grid.setLinSpaced(10, 0.5, 5.0);
  pc::KPMConfig dosCfg;
  dosCfg.numMoments = 256;
  dosCfg.numRandomVectors = 10;
  Eigen::ArrayXd dos = nma.densityOfStates(grid, dosCfg);
  CHECK(dos.size() == 10);
  for (Eigen::Index i = 0; i < dos.size(); ++i) {
    CHECK(dos[i] >= 0.0);
  }
}

TEST_CASE("NormalModeAnalyzer SLQ and Hybrid approximate spectra match exact",
          "[analysis][NormalModeAnalyzer][SLQ]")
{
  std::vector<size_t> springFrom;
  std::vector<size_t> springTo;
  for (size_t i = 0; i < 119; ++i) {
    springFrom.push_back(i);
    springTo.push_back(i + 1);
  }
  for (size_t i = 0; i < 100; i += 10) {
    springFrom.push_back(i);
    springTo.push_back(i + 20);
  }
  for (size_t i = 120; i < 129; ++i) {
    springFrom.push_back(i);
    springTo.push_back(i + 1);
  }

  pc::NormalModeAnalyzer nma(springFrom, springTo, 132);
  pc::NormalModeSpectrum exactSpec = nma.solveExact();

  Eigen::ArrayXd omega(6);
  omega << 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0;
  auto [gpExact, gppExact] = exactSpec.moduli(omega);

  Eigen::ArrayXd t(3);
  t << 0.0, 0.5, 5.0;
  Eigen::ArrayXd gtExact = exactSpec.stressAutocorrelation(t);

  pc::SLQConfig slqCfg;
  slqCfg.numSteps = 120;
  slqCfg.numRandomVectors = 30;
  pc::NormalModeSpectrum slqSpec = nma.solveApprox(slqCfg);
  auto [gpSLQ, gppSLQ] = slqSpec.moduli(omega);
  Eigen::ArrayXd gtSLQ = slqSpec.stressAutocorrelation(t);
  REQUIRE(gpSLQ.size() == 6);
  REQUIRE(gppSLQ.size() == 6);
  REQUIRE(gtSLQ.size() == 3);
  CHECK(gtSLQ[0] >= gtSLQ[1]);
  CHECK(gtSLQ[1] >= gtSLQ[2]);
  CHECK_THAT(gtSLQ[0], Catch::Matchers::WithinRel(gtExact[0], 0.15));
  for (Eigen::Index i = 0; i < omega.size(); ++i) {
    CHECK_THAT(gpSLQ[i], Catch::Matchers::WithinRel(gpExact[i], 0.25));
    CHECK_THAT(gppSLQ[i], Catch::Matchers::WithinRel(gppExact[i], 0.25));
  }

  pc::HybridConfig hybridCfg;
  hybridCfg.numSteps = 80;
  hybridCfg.numRandomVectors = 30;
  hybridCfg.numExactLowModes = 30;
  pc::NormalModeSpectrum hybridSpec = nma.solveApprox(hybridCfg);
  auto [gpHybrid, gppHybrid] = hybridSpec.moduli(omega);
  for (Eigen::Index i = 0; i < omega.size(); ++i) {
    CHECK_THAT(gpHybrid[i], Catch::Matchers::WithinRel(gpExact[i], 0.15));
    CHECK_THAT(gppHybrid[i], Catch::Matchers::WithinRel(gppExact[i], 0.15));
  }

  pc::NormalModeSpectrum autoSpec = nma.solveApprox();
  pc::KPMConfig kpmCfg;
  kpmCfg.numMoments = 256;
  kpmCfg.numRandomVectors = 20;
  pc::NormalModeSpectrum kpmSpec = nma.solveApprox(kpmCfg);
  CHECK(autoSpec.storageModulus(omega).allFinite());
  CHECK(kpmSpec.storageModulus(omega).allFinite());
}

TEST_CASE("NormalModeAnalyzer edge cases and multigraphs",
          "[analysis][NormalModeAnalyzer]")
{
  REQUIRE_THROWS_AS(pc::NormalModeAnalyzer({}, {}), std::invalid_argument);
  REQUIRE_THROWS_AS(pc::NormalModeAnalyzer({ 0, 1 }, { 1 }),
                    std::invalid_argument);

  pc::NormalModeAnalyzer dimerParallel({ 0, 0 }, { 1, 1 });
  Eigen::VectorXd dimerEvals =
    dimerParallel.solveExact(pc::NormalModeSolver::AUTO, false, true)
      .getEigenvalues();
  REQUIRE(dimerEvals.size() == 2);
  CHECK_THAT(dimerEvals[0], Catch::Matchers::WithinAbs(0.0, 1e-9));
  CHECK_THAT(dimerEvals[1], Catch::Matchers::WithinAbs(4.0, 1e-9));

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
  for (size_t i = 35; i < 44; ++i) {
    springFrom.push_back(i);
    springTo.push_back(i + 1);
  }

  pc::NormalModeAnalyzer largeSystem(springFrom, springTo);
  Eigen::VectorXd refEvals =
    largeSystem.solveExact(pc::NormalModeSolver::DENSE_LAPACK, false, false)
      .getEigenvalues();
  Eigen::VectorXd slicingEvals =
    largeSystem.solveExact(pc::NormalModeSolver::SPARSE_SLICING, false, true)
      .getEigenvalues();

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

  pc::NormalModeAnalyzer nma(springFrom, springTo);

  BENCHMARK("solveExact (eigenvalues only)")
  {
    return nma.solveExact(pc::NormalModeSolver::AUTO, false).getNumModes();
  };

  BENCHMARK("solveExact (with eigenvectors)")
  {
    return nma.solveExact(pc::NormalModeSolver::AUTO, true).getNumModes();
  };
}

TEST_CASE("NormalModeAnalyzer Metal GPU solver computes accurate eigenvalues "
          "and eigenvectors",
          "[analysis][NormalModeAnalyzer][Metal]")
{
#ifdef USE_METAL
  if (!pc::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
#endif
  {
    pc::NormalModeAnalyzer nma({ 0, 1, 2 }, { 1, 2, 0 });
    pc::NormalModeSpectrum specCPU =
      nma.solveExact(pc::NormalModeSolver::DENSE_LAPACK, true);
    pc::NormalModeSpectrum specGPU =
      nma.solveExact(pc::NormalModeSolver::METAL, true, false);

    REQUIRE(specGPU.getNumModes() == 3);
    for (Eigen::Index i = 0; i < 3; ++i) {
      CHECK_THAT(specGPU.getEigenvalues()[i],
                 Catch::Matchers::WithinAbs(specCPU.getEigenvalues()[i], 1e-4));
    }

    REQUIRE(specGPU.getEigenvectors().has_value());
    const Eigen::MatrixXd V = *specGPU.getEigenvectors();
    Eigen::MatrixXd Id = V.transpose() * V;
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        CHECK_THAT(Id(r, c),
                   Catch::Matchers::WithinAbs(r == c ? 1.0 : 0.0, 1e-3));
      }
    }
  }

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

    pc::NormalModeAnalyzer nma(springFrom, springTo);
    pc::NormalModeSpectrum specCPU =
      nma.solveExact(pc::NormalModeSolver::DENSE_LAPACK, false);
    pc::NormalModeSpectrum specGPU =
      nma.solveExact(pc::NormalModeSolver::METAL, false);

    REQUIRE(specGPU.getNumModes() == n);
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(n); ++i) {
      CHECK_THAT(specGPU.getEigenvalues()[i],
                 Catch::Matchers::WithinAbs(specCPU.getEigenvalues()[i], 1e-3));
    }
  }
}

TEST_CASE("NormalModeAnalyzer Metal KPM computes accurate Chebyshev moments",
          "[analysis][NormalModeAnalyzer][Metal][KPM]")
{
#ifdef USE_METAL
  if (!pc::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
#endif
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
  const size_t numMoments = 128;
  const size_t numRandomVectors = 10;

  pc::KPMConfig cpuCfg{
    numMoments, numRandomVectors, 2048, 4, pc::ComputeBackend::CPU
  };
  pc::KPMConfig b4Cfg{
    numMoments, numRandomVectors, 2048, 4, pc::ComputeBackend::METAL
  };
  pc::KPMConfig b8Cfg{
    numMoments, numRandomVectors, 2048, 8, pc::ComputeBackend::METAL
  };
  pc::KPMConfig b1Cfg{
    numMoments, numRandomVectors, 2048, 1, pc::ComputeBackend::METAL
  };

  Eigen::VectorXd cpuMoments = nma.chebyshevMoments(cpuCfg);
  Eigen::VectorXd gpuMomentsB4 = nma.chebyshevMoments(b4Cfg);
  Eigen::VectorXd gpuMomentsB8 = nma.chebyshevMoments(b8Cfg);
  Eigen::VectorXd gpuMomentsB1 = nma.chebyshevMoments(b1Cfg);

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
#ifdef USE_METAL
  if (!pc::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
#endif
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
    pc::KPMConfig cfg{ 256, 16, 2048, 1, pc::ComputeBackend::METAL };
    return nma.chebyshevMoments(cfg).size();
  };

  BENCHMARK("Metal GPU KPM Multi-RHS Block-SpMM (B=4, 2000 nodes, 256 moments, "
            "16 vectors)")
  {
    pc::KPMConfig cfg{ 256, 16, 2048, 4, pc::ComputeBackend::METAL };
    return nma.chebyshevMoments(cfg).size();
  };

  BENCHMARK("Metal GPU KPM Multi-RHS Block-SpMM (B=8, 2000 nodes, 256 moments, "
            "16 vectors)")
  {
    pc::KPMConfig cfg{ 256, 16, 2048, 8, pc::ComputeBackend::METAL };
    return nma.chebyshevMoments(cfg).size();
  };
}
