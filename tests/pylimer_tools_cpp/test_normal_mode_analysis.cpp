#include "../../src/pylimer_tools_cpp/calc/NormalModeAnalyzer.h"

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
  Eigen::Matrix3d expectedEigenvectors;
  expectedEigenvectors << -0.57735026918962584, 0, 0.81649658092772615, //
    -0.57735026918962562, -0.70710678118654746, -0.40824829046386291,   //
    -0.57735026918962584, 0.70710678118654746, -0.40824829046386302;

  CHECK(eigenvectors.isApprox(expectedEigenvectors, 1e-10));

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