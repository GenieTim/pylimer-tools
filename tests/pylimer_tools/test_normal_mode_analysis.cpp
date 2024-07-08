#include "../../src/pylimer_tools_cpp/calc/NormalModeAnalyzer.h"


#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace pc = pylimer_tools::calc;

TEST_CASE("NormalModeAnalyzer does not crash", "[analysis][NormalModeAnalyzer]") {
  // very small melt chain
  pc::NormalModeAnalyzer normalModeAnalyzer = pc::NormalModeAnalyzer(
    {1, 2, 3},
    {2, 3, 4}
  );

  REQUIRE_THROWS(normalModeAnalyzer.getEigenvalues());
  REQUIRE_THROWS(normalModeAnalyzer.getEigenvectors());
  REQUIRE_NOTHROW(normalModeAnalyzer.computeAllEigenvalues());
  Eigen::VectorXd eigenvalues = normalModeAnalyzer.getEigenvalues();
  CHECK(eigenvalues.size() == 4);
}
