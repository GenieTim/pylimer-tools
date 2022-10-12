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

TEST_CASE("Eigen behaves as required", "[analysis][MEHPForceBalance][Eigen]")
{
  Eigen::VectorXd testVec = Eigen::VectorXd::Zero(10);
  Eigen::ArrayXi testIdx = Eigen::ArrayXi::Zero(5);
  testIdx << 0, 0, 5, 5, 1;
  testVec(testIdx) += Eigen::VectorXd::Ones(5);
  REQUIRE(testVec[5] == Catch::Approx(2.));
  REQUIRE(testVec[0] == Catch::Approx(2.));
  REQUIRE(testVec[1] == Catch::Approx(1.));
  REQUIRE(testVec[2] == 0.0);
}

TEST_CASE("MEHP Force Balance runs", "[analysis][MEHPForceBalance]")
{
  // return;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("MEHP Force Balance 3D case")
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
      CHECK_NOTHROW(forceBalancer2.validateNetwork());
      CHECK(forceBalancer2.getNrOfSprings() == 8142);
      CHECK(forceBalancer2.getNrOfIterations() > 1);
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);

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

  SECTION("MEHP Force Balance 2D case")
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
    CHECK_NOTHROW(forceBalancer.validateNetwork());
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
    CHECK_NOTHROW(forceBalancer2.validateNetwork());
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

TEST_CASE("MEHP Force Balance handles slip-links",
          "[analysis][MEHPForceBalance]")
{
  // construct an example network
  pe::Universe universe = pe::Universe(1.0, 1.0, 1.0);
  /**
   * The system looks like this (in terms of bonds, not 3D placement):
   *
   * 1-5-2
   * |  /|
   * 8 9 6
   * |/  |
   * 4-7-3
   */
  universe.setBox(pe::Box(-15.0, 15.0, -15.0, 15.0, -15.0, 15.0));
  universe.addAtoms(9,
                    { { 1, 2, 3, 4, 5, 6, 7, 8, 9 } },     // id
                    { { 2, 2, 2, 2, 1, 1, 1, 1, 1 } },     // type
                    { { -5, 5, 5, -5, 0, 5, -5, -5, 0 } }, // x
                    { { 5, 5, -5, -5, 5, 0, 0, 0, 0 } },   // y
                    { { 2, 2, 2, 2, 2, 2, 2, 2, 2 } },     // z
                    { { 1, 1, 1, 1, 1, 1, 1, 1, 1 } },     // nx
                    { { 1, 1, 1, 1, 1, 1, 1, 1, 1 } },     // ny
                    { { 1, 1, 1, 1, 1, 1, 1, 1, 1 } }      // nz
  );
  universe.addBonds(10,
                    { { 1, 2, 2, 2, 3, 3, 7, 8, 1, 9 } },
                    { { 5, 5, 9, 6, 6, 7, 4, 4, 8, 4 } },
                    { { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 } },
                    false,
                    false);
  // now, construct the force balancer
  pcm::MEHPForceBalance forceBalancer2 =
    pcm::MEHPForceBalance(universe, 2, false);
  REQUIRE(forceBalancer2.getNrOfNodes() == forceBalancer2.getNrOfLinks());
  REQUIRE(forceBalancer2.getNrOfNodes() == 4);
  REQUIRE(forceBalancer2.getNrOfSprings() == 5);
  Eigen::VectorXd springPartitions0 = forceBalancer2.getSpringPartitions();
  pcm::ForceBalanceNetwork net0 = forceBalancer2.getNetwork();
  for (int i = 0; i < net0.nrOfPartialSprings; i++) {
    std::cout << net0.springPartIndexA[i] << ", " << net0.springPartIndexB[i]
              << ": ";
    std::cout << springPartitions0[i] << std::endl;
  }
  std::cout << std::endl;
  // and add slip-links
  REQUIRE_THROWS(forceBalancer2.addSlipLinks({ { 1, 2, 3 } },
                                             { { 1, 2, 3 } },
                                             { { 1., 2., 3. } },
                                             { { 1., 2., 3. } },
                                             { { 1., 2. } }));
  REQUIRE_THROWS(forceBalancer2.addSlipLinks({ { 1, 2, 3 } },
                                             { { 1, 2 } },
                                             { { 1., 2., 3. } },
                                             { { 1., 2., 3. } },
                                             { { 1., 2., 3.0 } }));
  // and actually add slip-links
  REQUIRE_NOTHROW(forceBalancer2.addSlipLinks({ { 0, 3 } },
                                              { { 1, 3 } },
                                              { { 4.2, 1.3 } },
                                              { { 3.9, 1.2 } },
                                              { { 1.3, 1.1 } }));
  REQUIRE(forceBalancer2.getNrOfNodes() + 2 == forceBalancer2.getNrOfLinks());
  // and run with them.
  // Expect the slip-link of between two strands to converge to the central atom
  // Expect the slip-link of two strands to stay at 0.5, 0.5.
  Eigen::VectorXd springPartitions = forceBalancer2.getSpringPartitions();
  Eigen::VectorXd displacements =
    Eigen::VectorXd::Zero(forceBalancer2.getNrOfLinks() * 3);
  pcm::ForceBalanceNetwork net = forceBalancer2.getNetwork();
  for (int i = 0; i < net.nrOfPartialSprings; i++) {
    std::cout << net.springPartIndexA[i] << ", " << net.springPartIndexB[i]
              << ": ";
    std::cout << springPartitions[i] << std::endl;
  }
  std::cout << std::endl;

  forceBalancer2.displaceToMeanPosition(
    &net, displacements, springPartitions, 4);
  CHECK(net.coordinates[4 * 3] + displacements[4 * 3] == Catch::Approx(2.5));
  CHECK(net.coordinates[4 * 3 + 1] + displacements[4 * 3 + 1] ==
        Catch::Approx(2.5));
  CHECK(net.coordinates[4 * 3 + 2] + displacements[4 * 3 + 2] ==
        Catch::Approx(2));

  // reset
  displacements.setZero();

  REQUIRE(net.springIndicesOfLinks.size() == net.nrOfLinks);
  for (int i = 4; i < 5; ++i) {
    forceBalancer2.displaceToMeanPosition(
      &net, displacements, springPartitions, 5);
    forceBalancer2.updateSpringPartition(
      &net, displacements, springPartitions, 5);
  for (int i = 0; i < net.nrOfPartialSprings; i++) {
    std::cout << net.springPartIndexA[i] << ", " << net.springPartIndexB[i]
              << ": ";
    std::cout << springPartitions[i] << std::endl;
  }
  std::cout << std::endl;
  }
  for (int i = 0; i < 125; ++i) {
    // do some random 125 steps with these two slip-links
    // NOTE: difficulty: finding out which node and spring it is actually
    // after the removal of strand atoms
    forceBalancer2.displaceToMeanPosition(
      &net, displacements, springPartitions, 4);
    // for (int dir = 0; dir < 3; ++dir) {
    //   std::cout << (net.coordinates[3 * 4 + dir] + displacements[3 * 4 +
    //   dir])
    //             << ", ";
    // }
    forceBalancer2.updateSpringPartition(
      &net, displacements, springPartitions, 4);
    // std::cout << std::endl;
    // std::cout << springPartitions[0][0] << ", " << springPartitions[1][0]
    //           << std::endl;
  }
  // assert expectations are met.
  // NOTE: difficulty: finding out which spring idx it actually is
  // for (int i = 0; i < net.nrOfSprings; ++i) {
  //   std::cout << "Spring " << i << " " << std::endl;
  //   for (int j = 0; j < net.linkIndicesOfSprings[i].size(); ++j) {
  //     std::cout << net.linkIndicesOfSprings[i][j] << " :";
  //     for (int dir = 0; dir < 3; ++dir) {
  //       std::cout
  //         << (net.coordinates[3 * net.linkIndicesOfSprings[i][j] + dir] +
  //             displacements[3 * net.linkIndicesOfSprings[i][j] + dir])
  //         << ", ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  // }
  for (int i = 0; i < net.nrOfPartialSprings; i++) {
    std::cout << net.springPartIndexA[i] << ", " << net.springPartIndexB[i]
              << ": ";
    std::cout << springPartitions[i] << std::endl;
  }
  std::cout << std::endl;
  CHECK(springPartitions[5]+1e-5 == Catch::Approx(0.0+1e-5).epsilon(1e-6)); // 4-1
  CHECK(springPartitions[6] == Catch::Approx(1.0).epsilon(1e-6)); // 4-2
  CHECK(springPartitions[0] == Catch::Approx(1.0).epsilon(1e-6)); // 4-0
  CHECK(springPartitions[1]+1e-5 == Catch::Approx(0.0+1e-5).epsilon(1e-6)); // 1-4
  // CHECK(springPartitions[8] == Catch::Approx(1.0).margin(1e-6)); // 5-3
  // CHECK(springPartitions[7] + 1e-5 == Catch::Approx(0.0 + 1e-5).margin(1e-6)); // 5-5
  // CHECK(springPartitions[3] == Catch::Approx(1.0).margin(1e-6)); // 5-0
}

TEST_CASE("MEHP Force Balance runs with non-network",
          "[analysis][MEHPForceBalance][NonGaussianSpringForceEvaluator]")
{
  // return;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  REQUIRE(universeSeq.getLength() == 0);
  std::string suspectedPath = "../pylimer_tools/fixtures/";

  SECTION("MEHP Force Balance 3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      pcm::MEHPForceBalance forceBalancer2 =
        pcm::MEHPForceBalance(universe2, 2, false);
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
  pcm::MEHPForceBalance forceBalancer =
    pcm::MEHPForceBalance(universe, 2, false);
  REQUIRE_NOTHROW(forceBalancer.runForceRelaxation(50000, 1e-12));
  REQUIRE(forceBalancer.getNrOfIterations() > 0);
  CHECK(forceBalancer.getNrOfActiveSprings() == 0);
  // CHECK(forceBalancer.getAverageSpringLength() ==
  //       Catch::Approx(0.0));
  CHECK(forceBalancer.getAverageSpringLength() >= 0.0);
  CHECK(forceBalancer.getAverageSpringLength() <= 3e-6);
  REQUIRE_NOTHROW(forceBalancer.validateNetwork());

  pe::Universe resultingUniverse = forceBalancer.getCrosslinkerVerse();
  // auto distances = resultingUniverse.computeBondLengths();
  // for (auto i : distances) {
  //   std::cout << i << std::endl;
  // }
  // auto residuals = forceBalancer.getResiduals();
  // for (auto i : residuals) {
  //   std::cout << i << " ";
  // }
  // std::cout << std::endl;
  // std::cout << forceBalancer.getNrOfIterations() << ", "
  //           << forceBalancer.getForce() << ", "
  //           << forceBalancer.getResidualNorm() << std::endl;

  std::array<std::array<double, 3>, 3> stressTensorSimpleSpring =
    forceBalancer.getStressTensor();
  std::array<std::array<double, 3>, 3> stressTensorLangevin =
    forceBalancer.getStressTensor();
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      CHECK(stressTensorLangevin[i][j] + 1e-5 ==
            Catch::Approx(stressTensorSimpleSpring[i][j] + 1e-5));
    }
  }
}
