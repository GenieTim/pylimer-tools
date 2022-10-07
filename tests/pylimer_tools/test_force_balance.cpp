#include "../../src/pylimer_tools_cpp/calc/MEHPForceBalance.h"
#include "../../src/pylimer_tools_cpp/calc/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/calc/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <vector>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::calc::mehp;

TEST_CASE("MEHP Force Balance runs",
          "[analysis][MEHPForceBalance][SimpleSpringMEHPForceEvaluator]")
{
  // return;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);

      // BENCHMARK_ADVANCED("MEHP LD_MMA " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceBalance forceBalancer3 =
      //     pcm::MEHPForceBalance(universe2, 2);
      //   meter.measure([&forceBalancer3] {
      //     forceBalancer3.runForceRelaxation("LD_MMA");
      //     return forceBalancer3.getNrOfIterations();
      //   });
      // };
      // BENCHMARK_ADVANCED("MEHP LD_LBFGS " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceBalance forceBalancer3 =
      //     pcm::MEHPForceBalance(universe2, 2);
      //   meter.measure([&forceBalancer3] {
      //     forceBalancer3.runForceRelaxation("LD_LBFGS");
      //     return forceBalancer3.getNrOfIterations();
      //   });
      // };

      double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      pcm::MEHPForceBalance forceBalancer2 =
        pcm::MEHPForceBalance(universe2, 2);
      pcm::MEHPForceRelaxation forceRelaxer =
        pcm::MEHPForceRelaxation(universe2, 2);
      // the strands are different -> cannot compare the distances anymore
      // CHECK((forceBalancer2.getCurrentSpringDistances() -
      // forceRelaxer.getCurrentSpringDistances()).isMuchSmallerThan(1e-12));
      REQUIRE(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE(forceBalancer2.getNrOfIterations() == 0);
      REQUIRE(forceBalancer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
      CHECK(forceBalancer2.getVolume() ==
            Catch::Approx(97.383096 * 97.383096 * 97.383096));
      // initial system values
      CHECK(forceBalancer2.getPressure() ==
            Catch::Approx(forceRelaxer.getPressure()));
      CHECK(forceBalancer2.getPressure() == Catch::Approx(0.39911682390778536));
      REQUIRE_NOTHROW(forceBalancer2.runForceRelaxation());
      CHECK(forceBalancer2.getNrOfSprings() == 8142);
      CHECK(forceBalancer2.getNrOfIterations() > 1);

      // conversion factors
      double kb = 1.381e-23; // Boltzmann, J/K
      double T = 300.;       // Temperature, K
      double sigmaToNm = 0.616;
      double sigmaToM = sigmaToNm * 1.e-9;
      double slope = 0.00393 / (sigmaToNm * sigmaToNm); // sigma^2/(g/mol)
      double beadMass = 161.;                           // g/mol
      double Nb = 80.; // nr of beads per strand
      double conversionFactor =
        3. * kb * T / (slope * beadMass * Nb); // J/sigma^2
      CHECK(conversionFactor / (sigmaToM * sigmaToM) ==
            Catch::Approx(0.000245543));
      double nu =
        nrOfChains / (forceBalancer2.getVolume() * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      CHECK(forceBalancer2.getPressure() ==
            Catch::Approx(0.153806)); // LJ Units [?]
      CHECK(forceBalancer2.getPressure() * conversionFactor /
              (sigmaToM * sigmaToM * sigmaToM) ==
            Catch::Approx(61308.3)); // shear modulus from the pressure, MPa
      double nrOfChainCorrection =
        (forceBalancer2.getDefaultNrOfChains() / nrOfChains);
      double expectedNb2 = slope * Nb * beadMass;
      double nb2Correction =
        (forceBalancer2.getDefaultR0Square() / (expectedNb2));
      double gammaCorrectionFactor = nrOfChainCorrection * nb2Correction;
      CHECK(
        forceBalancer2.getGammaFactor() * nrOfChainCorrection *
          forceBalancer2.getDefaultR0Square() ==
        Catch::Approx(42.6132)); // as from conversion-less Mathematica script
      CHECK(forceBalancer2.getGammaFactor() * gammaCorrectionFactor * kb * T *
              nu ==
            Catch::Approx(61308.3)); // ANT shear modulus, Pa
      CHECK(forceBalancer2.getGammaFactor() * gammaCorrectionFactor ==
            Catch::Approx(0.319446)); // "correct" gamma factor
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      // TODO: find better, more accurate tests here
      CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
      CHECK(forceBalancer2.getNrOfActiveSprings() > 1);
      CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
      CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
            forceBalancer2.getNrOfNodes());

    } else {
      std::cout << "Skipping large file PDMS MEHP run" << std::endl;
      REQUIRE(true);
    }
  }

  SECTION("2D case")
  {
    REQUIRE(std::filesystem::exists(suspectedPath));
    universeSeq.initializeFromDataSequence(
      { { suspectedPath + "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
                          "333_2d_t_7500001.structure.out" } });
    REQUIRE(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPForceBalance forceBalancer =
      pcm::MEHPForceBalance(universe, 2, true);
    REQUIRE(forceBalancer.getExitReason() == pcm::ExitReason::UNSET);
    REQUIRE(forceBalancer.getNrOfIterations() == 0);
    REQUIRE(forceBalancer.getVolume() == Catch::Approx(universe.getVolume()));
    forceBalancer.runForceRelaxation(5);
    REQUIRE(forceBalancer.getNrOfNodes() != universe.getNrOfAtoms());
    REQUIRE(forceBalancer.getNrOfIterations() <= 5);
    REQUIRE(forceBalancer.getNrOfIterations() >= 1);
    REQUIRE(universe.getAtomsOfType(2).size() == 7200);
    REQUIRE(forceBalancer.getExitReason() == pcm::ExitReason::MAX_STEPS);
    // run again, this time fully
    pcm::MEHPForceBalance forceBalancer2 =
      pcm::MEHPForceBalance(universe, 2, true);
    forceBalancer2.runForceRelaxation(10000, 1e-10, 100, 1e-9);
    REQUIRE(forceBalancer2.getNrOfIterations() > 5);
    CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer2.getGammaFactor(25, forceBalancer2.getNrOfSprings()) ==
          Catch::Approx(1. / 3.).epsilon(0.001));
    auto stressTensor = forceBalancer2.getStressTensor();
    CHECK(forceBalancer2.getPressure() ==
          Catch::Approx(
            (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.)
            .epsilon(0.02));
    // TODO: find better, more accurate tests here
    CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
    CHECK(forceBalancer2.getNrOfActiveSprings() > 1);
    CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
    CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
          forceBalancer2.getNrOfNodes());

    pe::Universe universe3 = forceBalancer2.getCrosslinkerVerse();
    CHECK(universe3.getNrOfAtoms() == forceBalancer2.getNrOfNodes());
    CHECK(universe3.getNrOfBonds() == forceBalancer2.getNrOfSprings());
    CHECK(universe3.getAtomsOfType(2).size() == universe3.getNrOfAtoms());

    // try out different algorithms
    std::vector<std::string> algorithms = { "LD_MMA",
                                            // "LD_TNEWTON_PRECOND_RESTART",
                                            // "GD_STOGO",
                                            "LD_SLSQP",
                                            "GN_DIRECT" };

    // for (std::string algorithm : algorithms) {
    //   pcm::MEHPForceBalance forceBalancerN =
    //     pcm::MEHPForceBalance(universe, 2);
    //   std::cout << "Testing algorithm " << algorithm << std::endl;
    //   auto start = std::chrono::high_resolution_clock::now();
    //   forceBalancerN.runForceRelaxation(true, 15, algorithm.c_str(), 10000);
    //   auto stop = std::chrono::high_resolution_clock::now();
    //   CHECK(forceBalancerN.getGammaEq() ==
    //         Catch::Approx(forceBalancer2.getGammaEq()));
    //   CHECK(forceBalancerN.getFinalPressure() ==
    //         Catch::Approx(forceBalancer2.getFinalPressure()));
    //   auto duration = duration_cast<std::chrono::microseconds>(stop - start);
    //   std::cout << "Took: " << duration.count() << std::endl;
    // }
  }
}

TEST_CASE("MEHP Force Relaxation2 runs with non-gaussian force evaluator",
          "[analysis][MEHPForceBalance][NonGaussianSpringForceEvaluator]")
{
  // return;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("3D case")
  {
    // TODO: correct values (& forces?)
    std::string largeInputFile =
      suspectedPath + "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);

      // BENCHMARK_ADVANCED("MEHP LD_MMA " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceBalance forceBalancer3 =
      //     pcm::MEHPForceBalance(universe2, 2);
      //   meter.measure([&forceBalancer3] {
      //     forceBalancer3.runForceRelaxation("LD_MMA");
      //     return forceBalancer3.getNrOfIterations();
      //   });
      // };
      // BENCHMARK_ADVANCED("MEHP LD_LBFGS " + largeInputFile)
      // (Catch::Benchmark::Chronometer meter)
      // {
      //   pcm::MEHPForceBalance forceBalancer3 =
      //     pcm::MEHPForceBalance(universe2, 2);
      //   meter.measure([&forceBalancer3] {
      //     forceBalancer3.runForceRelaxation("LD_LBFGS");
      //     return forceBalancer3.getNrOfIterations();
      //   });
      // };

      double nrOfChains = 1.e4;
      CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
            Catch::Approx(nrOfChains));
      pcm::NonGaussianSpringForceEvaluator nonGaussianForceEvaluator =
        pcm::NonGaussianSpringForceEvaluator(1.0, 79, 0.98);
      pcm::MEHPForceBalance forceBalancer2 =
        pcm::MEHPForceBalance(universe2, 2, false, &nonGaussianForceEvaluator);
      REQUIRE(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE(forceBalancer2.getNrOfIterations() == 0);
      REQUIRE(forceBalancer2.getVolume() ==
              Catch::Approx(universe2.getVolume()));
      CHECK(forceBalancer2.getVolume() ==
            Catch::Approx(97.383096 * 97.383096 * 97.383096));
      // initial system values
      CHECK(forceBalancer2.getPressure() == Catch::Approx(0.39911682390778536));
      REQUIRE_NOTHROW(forceBalancer2.runForceRelaxation());
      // As long as gradient is unclear: gradient free methods, e.g.:
      // "LN_SBPLX", "LN_BOBYQA", "LN_NELDERMEAD",
      // "LN_COBYLA", "LN_NEWUOA_BOUND"
      CHECK(forceBalancer2.getNrOfSprings() == 8142);
      CHECK(forceBalancer2.getNrOfIterations() > 1);

      // conversion factors
      double kb = 1.381e-23; // Boltzmann, J/K
      double T = 300.;       // Temperature, K
      double sigmaToNm = 0.616;
      double sigmaToM = sigmaToNm * 1.e-9;
      double slope = 0.00393 / (sigmaToNm * sigmaToNm); // sigma^2/(g/mol)
      double beadMass = 161.;                           // g/mol
      double Nb = 80.; // nr of beads per strand
      double conversionFactor =
        3. * kb * T / (slope * beadMass * Nb); // J/sigma^2
      CHECK(conversionFactor / (sigmaToM * sigmaToM) ==
            Catch::Approx(0.000245543));
      double nu =
        nrOfChains / (forceBalancer2.getVolume() * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      auto stressTensor = forceBalancer2.getStressTensor();
      CHECK(
        forceBalancer2.getPressure() ==
        Catch::Approx(
          (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.)
          .epsilon(0.02));
      CHECK(forceBalancer2.getPressure() ==
            Catch::Approx(0.1538073308).epsilon(0.0001)); // LJ Units [?]
      CHECK(forceBalancer2.getPressure() * conversionFactor /
              (sigmaToM * sigmaToM * sigmaToM) ==
            Catch::Approx(61308.9809826224)
              .epsilon(0.0001)); // shear modulus from the pressure, MPa
      double nrOfChainCorrection =
        (forceBalancer2.getDefaultNrOfChains() / nrOfChains);
      double expectedNb2 = slope * Nb * beadMass;
      double nb2Correction =
        (forceBalancer2.getDefaultR0Square() / (expectedNb2));
      double gammaCorrectionFactor = nrOfChainCorrection * nb2Correction;
      CHECK(forceBalancer2.getGammaFactor() * nrOfChainCorrection *
              forceBalancer2.getDefaultR0Square() ==
            Catch::Approx(42.613678259).epsilon(0.0001));
      CHECK(forceBalancer2.getGammaFactor() * gammaCorrectionFactor * kb * T *
              nu ==
            Catch::Approx(61308.9809826224)
              .epsilon(0.0001)); // ANT shear modulus, Pa
      CHECK(
        forceBalancer2.getGammaFactor() * gammaCorrectionFactor ==
        Catch::Approx(0.3194493682).epsilon(0.0001)); // "correct" gamma factor
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      // TODO: find better, more accurate tests here
      CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
      CHECK(forceBalancer2.getNrOfActiveSprings() > 1);
      CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
      CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
            forceBalancer2.getNrOfNodes());
    } else {
      std::cout << "Skipping large file PDMS MEHP run" << std::endl;
      REQUIRE(true);
    }
  }
}

TEST_CASE(
  "MEHP Force Relaxation2 runs with Langevin force evaluator and non-network",
  "[analysis][MEHPForceBalance][NonGaussianSpringForceEvaluator]")
{
  // return;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      pcm::NonGaussianSpringForceEvaluator nonGaussianForceEvaluator =
        pcm::NonGaussianSpringForceEvaluator(1.0, 79, 0.98);
      pcm::MEHPForceBalance forceBalancer2 =
        pcm::MEHPForceBalance(universe2, 2, false, &nonGaussianForceEvaluator);
      REQUIRE(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      REQUIRE_NOTHROW(forceBalancer2.runForceRelaxation());
      CHECK(forceBalancer2.getNrOfIterations() > 1);
    }
  }
}
TEST_CASE("Free chains collapse",
          "[analysis][MEHPForceBalance][NonGaussianSpringForceEvaluator]["
          "SimpleSpringMEHPForceEvaluator]")
{
  size_t nrOfBeads = 30;
  size_t nrOfBeadsPerChain = 3;
  pe::Universe universe =
    pe::Universe(nrOfBeads * 10.0, nrOfBeads * 10.0, nrOfBeads * 10.0);
  std::vector<double> xPositions;
  xPositions.reserve(nrOfBeads);
  std::vector<double> yPositions;
  yPositions.reserve(nrOfBeads);
  std::vector<double> zPositions;
  zPositions.reserve(nrOfBeads);
  std::vector<long int> atomIds;
  atomIds.reserve(nrOfBeads);
  std::vector<int> atomTypes;
  atomTypes.reserve(nrOfBeads);
  std::vector<int> zeroInts;
  zeroInts.reserve(nrOfBeads);
  std::vector<long int> bondFrom;
  bondFrom.reserve(nrOfBeads - 1);
  std::vector<long int> bondTo;
  bondTo.reserve(nrOfBeads - 1);
  double offset = 10.0;
  for (int i = 0; i < nrOfBeads; ++i) {
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i);
    atomTypes.push_back(i % nrOfBeadsPerChain == 0 ? 2 : 1);
    zeroInts.push_back(0);
    if (i > 0) {
      bondFrom.push_back(i - 1);
      bondTo.push_back(i);
    }
  }
  universe.addAtoms(atomIds,
                    atomTypes,
                    xPositions,
                    yPositions,
                    zPositions,
                    zeroInts,
                    zeroInts,
                    zeroInts);
  universe.addBonds(bondFrom, bondTo);
  REQUIRE(universe.getNrOfAtoms() == nrOfBeads);
  REQUIRE(universe.getNrOfBonds() == nrOfBeads - 1);

  // now, check for every force evaluator, that the maximum entropy is when all
  // these beads overlap first, the gaussian spring one
  pcm::SimpleSpringMEHPForceEvaluator simpleSpringForceEvaluator =
    pcm::SimpleSpringMEHPForceEvaluator();
  pcm::MEHPForceBalance forceBalancerSimpleSpring =
    pcm::MEHPForceBalance(universe, 2, false, &simpleSpringForceEvaluator);
  REQUIRE_NOTHROW(forceBalancerSimpleSpring.runForceRelaxation());
  REQUIRE(forceBalancerSimpleSpring.getNrOfIterations() > 0);
  CHECK(forceBalancerSimpleSpring.getNrOfActiveSprings() == 0);
  // CHECK(forceBalancerSimpleSpring.getAverageSpringLength() ==
  //       Catch::Approx(0.0));
  CHECK(forceBalancerSimpleSpring.getAverageSpringLength() >= 0.0);
  CHECK(forceBalancerSimpleSpring.getAverageSpringLength() <= 3e-6);
  REQUIRE_NOTHROW(forceBalancerSimpleSpring.validateNetwork());

  // then, the non-gaussian one
  pcm::NonGaussianSpringForceEvaluator langevinForceEvaluator =
    pcm::NonGaussianSpringForceEvaluator(1.0, nrOfBeadsPerChain * 2, 1.0);
  pcm::MEHPForceBalance forceBalancerLangevin =
    pcm::MEHPForceBalance(universe, 2, false, &langevinForceEvaluator);
  pe::Universe resultingUniverse0 = forceBalancerLangevin.getCrosslinkerVerse();
  // auto distances0 = resultingUniverse0.computeBondLengths();
  // for (auto i : distances0) {
  //   std::cout << i << std::endl;
  // }
  // std::cout << forceBalancerLangevin.getNrOfIterations() << ", "
  //           << forceBalancerLangevin.getForce() << ", "
  //           << forceBalancerLangevin.getResidualNorm() << std::endl;

  forceBalancerLangevin.runForceRelaxation(500000,
                                           1e-15); // "LD_SLSQP", "LD_MMA"
  CHECK(forceBalancerLangevin.getNrOfActiveSprings() == 0);
  // CHECK(forceBalancerLangevin.getAverageSpringLength() ==
  // Catch::Approx(0.0));
  CHECK(forceBalancerLangevin.getAverageSpringLength() >= 0.0);
  CHECK(forceBalancerLangevin.getAverageSpringLength() <= 1e-6);
  REQUIRE(forceBalancerLangevin.getNrOfIterations() > 0);
  CHECK(forceBalancerLangevin.getExitReason() == pcm::ExitReason::X_TOLERANCE);

  pe::Universe resultingUniverse = forceBalancerLangevin.getCrosslinkerVerse();
  // auto distances = resultingUniverse.computeBondLengths();
  // for (auto i : distances) {
  //   std::cout << i << std::endl;
  // }
  // auto residuals = forceBalancerLangevin.getResiduals();
  // for (auto i : residuals) {
  //   std::cout << i << " ";
  // }
  // std::cout << std::endl;
  // std::cout << forceBalancerLangevin.getNrOfIterations() << ", "
  //           << forceBalancerLangevin.getForce() << ", "
  //           << forceBalancerLangevin.getResidualNorm() << std::endl;

  std::array<std::array<double, 3>, 3> stressTensorSimpleSpring =
    forceBalancerSimpleSpring.getStressTensor();
  std::array<std::array<double, 3>, 3> stressTensorLangevin =
    forceBalancerLangevin.getStressTensor();
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      CHECK(stressTensorLangevin[i][j] + 1e-5 ==
            Catch::Approx(stressTensorSimpleSpring[i][j] + 1e-5));
    }
  }
}
