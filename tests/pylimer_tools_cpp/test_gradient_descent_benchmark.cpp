#include "../../src/pylimer_tools_cpp/utils/ExtraEigenSolvers.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <random>

// Helper function to create SPD (Symmetric Positive Definite) matrix
Eigen::MatrixXd
createSPDMatrix(int n, double conditionNumber = 100.0, unsigned seed = 42)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> dis(0.0, 1.0);

  // Generate random orthogonal matrix using QR decomposition
  Eigen::MatrixXd randomMat = Eigen::MatrixXd::NullaryExpr(
    n, n, [&]() { return dis(gen) - 0.5; });
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(randomMat);
  Eigen::MatrixXd Q = qr.householderQ();

  // Create diagonal matrix with controlled eigenvalues
  Eigen::VectorXd eigenvalues = Eigen::VectorXd::LinSpaced(n, 1.0, conditionNumber);
  Eigen::MatrixXd D = eigenvalues.asDiagonal();

  // A = Q * D * Q^T is SPD with known condition number
  return Q * D * Q.transpose();
}

// Helper function to create general (non-SPD) matrix for BB methods
Eigen::MatrixXd
createGeneralMatrix(int n, unsigned seed = 42)
{
  std::mt19937 gen(seed);
  std::normal_distribution<> dis(0.0, 1.0);

  Eigen::MatrixXd A = Eigen::MatrixXd::NullaryExpr(
    n, n, [&]() { return dis(gen); });
  
  // Make it diagonally dominant to ensure convergence
  for (int i = 0; i < n; ++i) {
    A(i, i) = std::abs(A(i, i)) + n * 2.0;
  }
  
  return A;
}

// Helper to validate solution accuracy
double
computeRelativeResidual(const Eigen::MatrixXd& A,
                        const Eigen::VectorXd& x,
                        const Eigen::VectorXd& b)
{
  Eigen::VectorXd residual = A * x - b;
  double residualNorm = residual.norm();
  double bNorm = b.norm();
  return residualNorm / bNorm;
}

TEST_CASE("Gradient Descent Solvers - Correctness Tests",
          "[GradientDescent][correctness]")
{
  std::cout << "Running test \"Gradient Descent Solvers - Correctness Tests\""
            << std::endl;

  const int n = 50;
  const double tolerance = 1e-6;
  const int maxIterations = 10000;

  SECTION("Steepest Descent (SPD) - Correctness")
  {
    Eigen::MatrixXd A = createSPDMatrix(n, 100.0);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    int iterations = 0;
    Eigen::VectorXd xSolved = Eigen::gradientDescent(
      A, b, 0.01, tolerance, maxIterations, iterations);

    double relResidual = computeRelativeResidual(A, xSolved, b);
    std::cout << "  Steepest Descent: iterations=" << iterations
              << ", rel_residual=" << relResidual << std::endl;

    REQUIRE(relResidual < 1e-5);
    REQUIRE(iterations < maxIterations);

    // Check solution accuracy
    Eigen::VectorXd error = xSolved - xTrue;
    REQUIRE(error.norm() / xTrue.norm() < 0.01);
  }

  SECTION("Barzilai-Borwein (BB2) - Correctness")
  {
    Eigen::MatrixXd A = createGeneralMatrix(n);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    int iterations = 0;
    Eigen::VectorXd xSolved = Eigen::gradientDescentBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterations, false);

    double relResidual = computeRelativeResidual(A, xSolved, b);
    std::cout << "  BB2 (long): iterations=" << iterations
              << ", rel_residual=" << relResidual << std::endl;

    REQUIRE(relResidual < 1e-5);
    REQUIRE(iterations < maxIterations);

    Eigen::VectorXd error = xSolved - xTrue;
    REQUIRE(error.norm() / xTrue.norm() < 0.01);
  }

  SECTION("Barzilai-Borwein (BB1) - Correctness")
  {
    Eigen::MatrixXd A = createGeneralMatrix(n);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    int iterations = 0;
    Eigen::VectorXd xSolved = Eigen::gradientDescentBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterations, true);

    double relResidual = computeRelativeResidual(A, xSolved, b);
    std::cout << "  BB1 (short): iterations=" << iterations
              << ", rel_residual=" << relResidual << std::endl;

    REQUIRE(relResidual < 1e-5);
    REQUIRE(iterations < maxIterations);

    Eigen::VectorXd error = xSolved - xTrue;
    REQUIRE(error.norm() / xTrue.norm() < 0.01);
  }

  SECTION("Heavy Ball + Barzilai-Borwein - Correctness")
  {
    Eigen::MatrixXd A = createGeneralMatrix(n);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    int iterations = 0;
    Eigen::VectorXd xSolved = Eigen::gradientDescentHeavyBallBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterations);

    double relResidual = computeRelativeResidual(A, xSolved, b);
    std::cout << "  Heavy Ball BB: iterations=" << iterations
              << ", rel_residual=" << relResidual << std::endl;

    REQUIRE(relResidual < 1e-5);
    REQUIRE(iterations < maxIterations);

    Eigen::VectorXd error = xSolved - xTrue;
    REQUIRE(error.norm() / xTrue.norm() < 0.01);
  }

  SECTION("Compare methods on well-conditioned SPD system")
  {
    Eigen::MatrixXd A = createSPDMatrix(n, 10.0); // well-conditioned
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    int iterSD = 0, iterBB1 = 0, iterBB2 = 0, iterHB = 0;

    Eigen::VectorXd xSD = Eigen::gradientDescent(
      A, b, 0.01, tolerance, maxIterations, iterSD);
    Eigen::VectorXd xBB1 = Eigen::gradientDescentBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterBB1, true);
    Eigen::VectorXd xBB2 = Eigen::gradientDescentBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterBB2, false);
    Eigen::VectorXd xHB = Eigen::gradientDescentHeavyBallBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterHB);

    std::cout << "  Well-conditioned SPD (cond=10):" << std::endl;
    std::cout << "    Steepest Descent: " << iterSD << " iterations" << std::endl;
    std::cout << "    BB1 (short):      " << iterBB1 << " iterations" << std::endl;
    std::cout << "    BB2 (long):       " << iterBB2 << " iterations" << std::endl;
    std::cout << "    Heavy Ball BB:    " << iterHB << " iterations" << std::endl;

    REQUIRE(computeRelativeResidual(A, xSD, b) < 1e-5);
    REQUIRE(computeRelativeResidual(A, xBB1, b) < 1e-5);
    REQUIRE(computeRelativeResidual(A, xBB2, b) < 1e-5);
    REQUIRE(computeRelativeResidual(A, xHB, b) < 1e-5);
  }

  SECTION("Compare methods on ill-conditioned SPD system")
  {
    Eigen::MatrixXd A = createSPDMatrix(n, 1000.0); // ill-conditioned
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    int iterSD = 0, iterBB1 = 0, iterBB2 = 0, iterHB = 0;

    Eigen::VectorXd xSD = Eigen::gradientDescent(
      A, b, 0.01, tolerance, maxIterations, iterSD);
    Eigen::VectorXd xBB1 = Eigen::gradientDescentBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterBB1, true);
    Eigen::VectorXd xBB2 = Eigen::gradientDescentBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterBB2, false);
    Eigen::VectorXd xHB = Eigen::gradientDescentHeavyBallBarzilaiBorwein(
      A, b, 0.01, tolerance, maxIterations, iterHB);

    std::cout << "  Ill-conditioned SPD (cond=1000):" << std::endl;
    std::cout << "    Steepest Descent: " << iterSD << " iterations" << std::endl;
    std::cout << "    BB1 (short):      " << iterBB1 << " iterations" << std::endl;
    std::cout << "    BB2 (long):       " << iterBB2 << " iterations" << std::endl;
    std::cout << "    Heavy Ball BB:    " << iterHB << " iterations" << std::endl;

    REQUIRE(computeRelativeResidual(A, xSD, b) < 1e-5);
    REQUIRE(computeRelativeResidual(A, xBB1, b) < 1e-5);
    REQUIRE(computeRelativeResidual(A, xBB2, b) < 1e-5);
    REQUIRE(computeRelativeResidual(A, xHB, b) < 1e-5);
  }
}

TEST_CASE("Gradient Descent Solvers - Small System Benchmarks",
          "[GradientDescent][benchmark][small]")
{
  std::cout << "Running test \"Gradient Descent Solvers - Small System Benchmarks\""
            << std::endl;

  const int n = 100;
  const double tolerance = 1e-8;
  const int maxIterations = 10000;

  SECTION("Small SPD system (n=100, cond=50)")
  {
    Eigen::MatrixXd A = createSPDMatrix(n, 50.0);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    BENCHMARK("Steepest Descent (SPD)")
    {
      int iterations = 0;
      return Eigen::gradientDescent(
        A, b, 0.01, tolerance, maxIterations, iterations);
    };

    BENCHMARK("BB1 (short)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, true);
    };

    BENCHMARK("BB2 (long)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, false);
    };

    BENCHMARK("Heavy Ball BB (momentum=0.7)")
    {
      int iterations = 0;
      return Eigen::gradientDescentHeavyBallBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, Eigen::VectorXd(), -1.0, nullptr, 0.7);
    };

    BENCHMARK("Eigen BiCGSTAB (reference)")
    {
      Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
      solver.setMaxIterations(maxIterations);
      solver.setTolerance(tolerance);
      solver.compute(A);
      return solver.solve(b);
    };

    BENCHMARK("Eigen ConjugateGradient (reference)")
    {
      Eigen::ConjugateGradient<Eigen::MatrixXd> solver;
      solver.setMaxIterations(maxIterations);
      solver.setTolerance(tolerance);
      solver.compute(A);
      return solver.solve(b);
    };
  }

  SECTION("Small general system (n=100)")
  {
    Eigen::MatrixXd A = createGeneralMatrix(n);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    BENCHMARK("BB1 (short)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, true);
    };

    BENCHMARK("BB2 (long)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, false);
    };

    BENCHMARK("Heavy Ball BB")
    {
      int iterations = 0;
      return Eigen::gradientDescentHeavyBallBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations);
    };

    BENCHMARK("Eigen BiCGSTAB (reference)")
    {
      Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
      solver.setMaxIterations(maxIterations);
      solver.setTolerance(tolerance);
      solver.compute(A);
      return solver.solve(b);
    };
  }
}

TEST_CASE("Gradient Descent Solvers - Medium System Benchmarks",
          "[GradientDescent][benchmark][medium]")
{
  std::cout << "Running test \"Gradient Descent Solvers - Medium System Benchmarks\""
            << std::endl;

  const int n = 500;
  const double tolerance = 1e-6;
  const int maxIterations = 10000;

  SECTION("Medium SPD system (n=500, cond=100)")
  {
    Eigen::MatrixXd A = createSPDMatrix(n, 100.0);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    BENCHMARK("Steepest Descent (SPD)")
    {
      int iterations = 0;
      return Eigen::gradientDescent(
        A, b, 0.01, tolerance, maxIterations, iterations);
    };

    BENCHMARK("BB1 (short)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, true);
    };

    BENCHMARK("BB2 (long)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, false);
    };

    BENCHMARK("Heavy Ball BB (momentum=0.7)")
    {
      int iterations = 0;
      return Eigen::gradientDescentHeavyBallBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, Eigen::VectorXd(), -1.0, nullptr, 0.7);
    };

    BENCHMARK("Eigen ConjugateGradient (reference)")
    {
      Eigen::ConjugateGradient<Eigen::MatrixXd> solver;
      solver.setMaxIterations(maxIterations);
      solver.setTolerance(tolerance);
      solver.compute(A);
      return solver.solve(b);
    };
  }
}

TEST_CASE("Gradient Descent Solvers - Large System Benchmarks",
          "[GradientDescent][benchmark][large][long]")
{
  std::cout << "Running test \"Gradient Descent Solvers - Large System Benchmarks\""
            << std::endl;

  const int n = 1000;
  const double tolerance = 1e-6;
  const int maxIterations = 10000;

  SECTION("Large SPD system (n=1000, cond=200)")
  {
    Eigen::MatrixXd A = createSPDMatrix(n, 200.0);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    BENCHMARK("Steepest Descent (SPD)")
    {
      int iterations = 0;
      return Eigen::gradientDescent(
        A, b, 0.01, tolerance, maxIterations, iterations);
    };

    BENCHMARK("BB2 (long)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, false);
    };

    BENCHMARK("Heavy Ball BB (momentum=0.7)")
    {
      int iterations = 0;
      return Eigen::gradientDescentHeavyBallBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, Eigen::VectorXd(), -1.0, nullptr, 0.7);
    };

    BENCHMARK("Eigen ConjugateGradient (reference)")
    {
      Eigen::ConjugateGradient<Eigen::MatrixXd> solver;
      solver.setMaxIterations(maxIterations);
      solver.setTolerance(tolerance);
      solver.compute(A);
      return solver.solve(b);
    };
  }

  SECTION("Large general system (n=1000)")
  {
    Eigen::MatrixXd A = createGeneralMatrix(n);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    BENCHMARK("BB2 (long)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, false);
    };

    BENCHMARK("Heavy Ball BB")
    {
      int iterations = 0;
      return Eigen::gradientDescentHeavyBallBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations);
    };

    BENCHMARK("Eigen BiCGSTAB (reference)")
    {
      Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
      solver.setMaxIterations(maxIterations);
      solver.setTolerance(tolerance);
      solver.compute(A);
      return solver.solve(b);
    };
  }
}

TEST_CASE("Gradient Descent Solvers - Sparse System Benchmarks",
          "[GradientDescent][benchmark][sparse]")
{
  std::cout << "Running test \"Gradient Descent Solvers - Sparse System Benchmarks\""
            << std::endl;

  const int n = 500;
  const double tolerance = 1e-6;
  const int maxIterations = 10000;

  SECTION("Sparse SPD system (n=500, 5% sparsity)")
  {
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Create sparse SPD matrix (tridiagonal with diagonal dominance)
    Eigen::SparseMatrix<double> Asparse(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    
    for (int i = 0; i < n; ++i) {
      triplets.push_back(Eigen::Triplet<double>(i, i, 4.0 + dis(gen)));
      if (i > 0) {
        double val = dis(gen);
        triplets.push_back(Eigen::Triplet<double>(i, i - 1, val));
        triplets.push_back(Eigen::Triplet<double>(i - 1, i, val));
      }
    }
    Asparse.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::MatrixXd A = Eigen::MatrixXd(Asparse);

    Eigen::VectorXd xTrue = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = A * xTrue;

    BENCHMARK("Steepest Descent (SPD)")
    {
      int iterations = 0;
      return Eigen::gradientDescent(
        A, b, 0.01, tolerance, maxIterations, iterations);
    };

    BENCHMARK("BB2 (long)")
    {
      int iterations = 0;
      return Eigen::gradientDescentBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations, false);
    };

    BENCHMARK("Heavy Ball BB")
    {
      int iterations = 0;
      return Eigen::gradientDescentHeavyBallBarzilaiBorwein(
        A, b, 0.01, tolerance, maxIterations, iterations);
    };

    BENCHMARK("Eigen Sparse ConjugateGradient (reference)")
    {
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
      solver.setMaxIterations(maxIterations);
      solver.setTolerance(tolerance);
      solver.compute(Asparse);
      return solver.solve(b);
    };

    BENCHMARK("Eigen Dense ConjugateGradient (reference)")
    {
      Eigen::ConjugateGradient<Eigen::MatrixXd> solver;
      solver.setMaxIterations(maxIterations);
      solver.setTolerance(tolerance);
      solver.compute(A);
      return solver.solve(b);
    };
  }
}

TEST_CASE("Gradient Descent Solvers - Convergence Analysis",
          "[GradientDescent][convergence]")
{
  std::cout << "Running test \"Gradient Descent Solvers - Convergence Analysis\""
            << std::endl;

  const int n = 100;
  const double tolerance = 1e-10;
  const int maxIterations = 5000;

  SECTION("Convergence rate comparison on well-conditioned problem")
  {
    Eigen::MatrixXd A = createSPDMatrix(n, 10.0);
    Eigen::VectorXd xTrue = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd b = A * xTrue;

    std::vector<double> residuals;
    
    // Test each method and track convergence
    auto testMethod = [&](auto solver, const std::string& name) {
      int iter = 0;
      residuals.clear();
      
      auto callback = [&](int iteration, const Eigen::VectorXd& x) {
        double relResidual = computeRelativeResidual(A, x, b);
        residuals.push_back(relResidual);
        return false; // Continue
      };
      
      Eigen::VectorXd xSolved = solver(A, b, callback, iter);
      
      std::cout << "  " << name << ": " << iter << " iterations, "
                << "final residual=" << residuals.back() << std::endl;
      
      return iter;
    };

    int iterSD = testMethod(
      [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, auto callback, int& iter) {
        return Eigen::gradientDescent(
          A, b, 0.01, tolerance, maxIterations, iter, Eigen::VectorXd(), -1.0, callback);
      },
      "Steepest Descent"
    );

    int iterBB2 = testMethod(
      [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, auto callback, int& iter) {
        return Eigen::gradientDescentBarzilaiBorwein(
          A, b, 0.01, tolerance, maxIterations, iter, false, Eigen::VectorXd(), -1.0, callback);
      },
      "BB2 (long)"
    );

    int iterHB = testMethod(
      [&](const Eigen::MatrixXd& A, const Eigen::VectorXd& b, auto callback, int& iter) {
        return Eigen::gradientDescentHeavyBallBarzilaiBorwein(
          A, b, 0.01, tolerance, maxIterations, iter, Eigen::VectorXd(), -1.0, callback, 0.7);
      },
      "Heavy Ball BB"
    );

    // All should converge
    REQUIRE(iterSD < maxIterations);
    REQUIRE(iterBB2 < maxIterations);
    REQUIRE(iterHB < maxIterations);
  }
}
