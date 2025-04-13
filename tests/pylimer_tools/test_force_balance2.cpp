#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance2.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceRelaxation.h"
#include "../../src/pylimer_tools_cpp/topo/EntanglementDetector.h"
#include "../../src/pylimer_tools_cpp/utils/ExtraEigenTypes.h"
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include "../../src/pylimer_tools_cpp/sim/MEHPForceBalance.h"

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::sim::mehp;

#ifndef PYLIMER_TEST_FIXTURES_DIR
#define PYLIMER_TEST_FIXTURES_DIR "../pylimer_tools/tests/fixtures"
#endif

TEST_CASE("MEHP Force Balance2 runs", "[analysis][MEHPForceBalance2][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 runs\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Balance2 3D case")
  {
    std::string largeInputFile =
      suspectedPath +
      "/xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    REQUIRE(std::filesystem::exists(largeInputFile));
    std::cout << "Reading file " << largeInputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { largeInputFile } });
    pe::Universe universe2 = universeSeq.atIndex(0);
    std::cout << "Read file " << largeInputFile << std::endl;

    double nrOfChains = 1.e4;
    CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
          Catch::Approx(nrOfChains));
    pcm::MEHPForceBalance2 forceBalancer2 =
      pcm::MEHPForceBalance2(universe2, 2);
    size_t nStrandsFb = forceBalancer2.getNrOfStrands();
    size_t nSpringsFb = forceBalancer2.getNrOfSprings();

    // SECTION("Stress tensor computations are equivalent")
    // {
    Eigen::Matrix3d stressTensor1 = forceBalancer2.getStressTensor();
    Eigen::Matrix3d stressTensor2 = forceBalancer2.getStressTensorLinkBased();
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        CHECK(stressTensor1(i, j) == Catch::Approx(stressTensor2(i, j)));
      }
    }

    SECTION("Actual balance results in correct phantom results")
    {
      std::cout << "Doing phantom force balance" << std::endl;
      pcm::MEHPForceRelaxation forceRelaxer =
        pcm::MEHPForceRelaxation(universe2, 2);
      forceRelaxer.configAssumeBoxLargeEnough(false);

      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      CHECK(forceBalancer2.getNrOfIterations() == 0);
      CHECK(forceBalancer2.getVolume() == Catch::Approx(universe2.getVolume()));
      CHECK(forceBalancer2.getVolume() ==
            Catch::Approx(97.383096 * 97.383096 * 97.383096));
      // initial system values
      CHECK_THAT(forceBalancer2.getPressure(),
                 Catch::Matchers::WithinAbs(forceRelaxer.getPressure(), 1e-5));
      CHECK(forceBalancer2.getPressure() == Catch::Approx(0.0061105865));
      CHECK_THAT(
        forceRelaxer.getResidual(),
        Catch::Matchers::WithinAbs(forceBalancer2.getResidual(), 1e-5));

      // run force relaxation
      CHECK_NOTHROW(forceBalancer2.runForceRelaxation());
      CHECK_NOTHROW(forceBalancer2.validateNetwork());
      CHECK(forceBalancer2.getNrOfStrands() == nStrandsFb);
      CHECK(forceBalancer2.getNrOfSprings() == nSpringsFb);
      CHECK(forceBalancer2.getNrOfIterations() > 1);
      CHECK(forceBalancer2.getResidual() < 1e-5);
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
        (forceBalancer2.getNetwork().springContourLength.mean() / Nb) * 3. *
        kb * T / (slope * beadMass); // J/sigma^2
      CHECK(conversionFactor / (sigmaToM * sigmaToM * 79.) ==
            Catch::Approx(0.0002450018));
      double nu =
        nrOfChains / (forceBalancer2.getVolume() * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      CHECK(forceBalancer2.getPressure() ==
            Catch::Approx(0.153806 / 79.)); // LJ Units [?]
      CHECK(forceBalancer2.getPressure() * conversionFactor /
              (sigmaToM * sigmaToM * sigmaToM) ==
            Catch::Approx(61172.8878)); // shear modulus from the pressure, MPa
      double b02 = slope * beadMass;
      // (forceBalancer2.getDefaultR0Square() / (expectedNb2));
      CHECK(forceBalancer2.getGammaFactor(1.) ==
            Catch::Approx(forceBalancer2.getGammaFactors(1.).mean()));
      CHECK(
        forceBalancer2.getGammaFactors(b02).mean() ==
        Catch::Approx((forceBalancer2.getGammaFactorsInDir(b02, 0).mean() +
                       forceBalancer2.getGammaFactorsInDir(b02, 1).mean() +
                       forceBalancer2.getGammaFactorsInDir(b02, 2).mean())));
      CHECK(forceBalancer2.getGammaFactor(b02) ==
            Catch::Approx(forceBalancer2.getGammaFactors(b02).mean()));
      CHECK_THAT(
        forceBalancer2.getGammaFactor(b02) * kb * T * nu,
        Catch::Matchers::WithinRel(61308.3, 0.03)); // ANT shear modulus, Pa
      CHECK_THAT(forceBalancer2.getGammaFactor(b02),
                 Catch::Matchers::WithinRel(
                   forceBalancer2.getGammaFactor(b02, nrOfChains), 0.03));
      CHECK_THAT(
        forceBalancer2.getGammaFactor(b02, nrOfChains),
        Catch::Matchers::WithinRel(
          0.319446, 0.03)); // "correct" gamma factor, see Mathematica script
      // conversion-less Mathematica script: 42.6132
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      // TODO: find better, more accurate tests here
      CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
      CHECK(forceBalancer2.getNrOfActiveStrands() > 1);
      CHECK(forceBalancer2.getNrOfActiveStrands() <=
            (forceBalancer2.getNrOfActiveStrandsInDir(0) +
             forceBalancer2.getNrOfActiveStrandsInDir(1) +
             forceBalancer2.getNrOfActiveStrandsInDir(2)));
      CHECK(forceBalancer2.getAverageSpringLength() > 1.0);
      CHECK(forceBalancer2.getEffectiveFunctionalityOfAtoms().size() ==
            forceBalancer2.getNrOfNodes());
    }
  }

  SECTION("MEHP Force Balance2 2D case")
  {
    CHECK(std::filesystem::exists(suspectedPath));
    universeSeq.initializeFromDataSequence(
      { { suspectedPath +
          "/structure/equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0."
          "333_2d_t_7500001.structure.out" } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2(universe, 2, true);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::UNSET);
    CHECK(forceBalancer.getNrOfIterations() == 0);
    CHECK(forceBalancer.getVolume() == Catch::Approx(universe.getVolume()));
    forceBalancer.runForceRelaxation();
    CHECK(forceBalancer.getNrOfNodes() < universe.getNrOfAtoms());
    CHECK(forceBalancer.getNrOfIterations() >= 1);
    CHECK(universe.getAtomsOfType(2).size() == 7200);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(forceBalancer.getGammaFactor(1., forceBalancer.getNrOfStrands()) ==
          Catch::Approx(1. / 3.).epsilon(0.001));
    auto stressTensor = forceBalancer.getStressTensor();
    CHECK(forceBalancer.getPressure() ==
          Catch::Approx(
            (stressTensor(0, 0) + stressTensor(1, 1) + stressTensor(2, 2)) / 3.)
            .epsilon(0.02));
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    // TODO: find better, more accurate tests here
    CHECK(forceBalancer.getNrOfActiveNodes() > 1);
    CHECK(forceBalancer.getNrOfActiveStrands() > 1);
    CHECK(forceBalancer.getAverageSpringLength() > 1.0);
    CHECK(forceBalancer.getEffectiveFunctionalityOfAtoms().size() ==
          forceBalancer.getNrOfNodes());

    pe::Universe universe3 = forceBalancer.getCrosslinkerVerse();
    CHECK(universe3.getNrOfAtoms() == forceBalancer.getNrOfNodes());
    CHECK(universe3.getNrOfBonds() == forceBalancer.getNrOfStrands());
    CHECK(universe3.getAtomsOfType(2).size() == universe3.getNrOfAtoms());
  }
}

TEST_CASE(
  "MEHP Force Balance2 can randomly add entanglements ignoring crosslinkers",
  "[analysis][MEHPForceBalance2]")
{
  std::cout << "Running test \"MEHP Force Balance2 can randomly add slip-links "
               "ignoring crosslinkers\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  if (std::filesystem::exists(inputFile)) {
    CHECK(std::filesystem::exists(suspectedPath));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    CHECK(universeSeq.getLength() == 1);
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file. " << std::endl;
    pcm::MEHPForceBalance2 forceBalancer =
      pcm::MEHPForceBalance2(universe, 250, 2.0, 0.0, 100);
    size_t nrOfAddedLinks =
      forceBalancer.getNrOfLinks() - forceBalancer.getNrOfNodes();
    CHECK(nrOfAddedLinks >= 50);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

    pcm::ForceBalance2Network net = forceBalancer.getNetwork();
    Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
    CHECK(net.nrOfSprings > 0);
    size_t numRemoved =
      forceBalancer.unlinkBifunctionalLinks(net, displacements);
    CHECK_NOTHROW(forceBalancer.validateNetwork(net, displacements));
    CHECK(numRemoved > 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork(net, displacements));
    // remove all springs...
    numRemoved = forceBalancer.removeInactiveLinks(net, displacements, 1e5);
    CHECK(net.nrOfSprings == 0);
    CHECK(numRemoved > 0);
  }
}

TEST_CASE("MEHP Force Balance2 can randomly add and remove entanglements",
          "[analysis][MEHPForceBalance2]")
{
  std::cout << "Running test \"MEHP Force Balance2 can randomly add and remove "
               "entanglements\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  CHECK(std::filesystem::exists(suspectedPath));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  CHECK(universeSeq.getLength() == 1);
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file. " << std::endl;
  pcm::MEHPForceBalance2 forceBalancer =
    pcm::MEHPForceBalance2(universe, 1000, 2.0, 0.0, 100);
  size_t nrOfAddedLinks =
    forceBalancer.getNrOfLinks() - forceBalancer.getNrOfNodes();
  CHECK(nrOfAddedLinks >= 100);
  // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

  pcm::ForceBalance2Network net = forceBalancer.getNetwork();
  Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
  size_t numRemoved = forceBalancer.unlinkBifunctionalLinks(net, displacements);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
  CHECK(numRemoved > 0);

  // run a while to get inactive links
  forceBalancer.runForceRelaxation();
  net = forceBalancer.getNetwork();
  displacements = forceBalancer.getCurrentDisplacements();
  // due to the randomness, it _could_ be one day that actually all strands
  // are active. unlikely, but I can imagine it to be possible.
  size_t numInactiveRemoved =
    forceBalancer.removeInactiveLinks(net, displacements, 0.1);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
  CHECK(numInactiveRemoved > 0);
  CHECK_NOTHROW(forceBalancer.validateNetwork(net, displacements));

  ////////////////////////////////////////////////////////////////
  forceBalancer = pcm::MEHPForceBalance2(universe, 1000, 2.0, 0.0, 100);
  nrOfAddedLinks = forceBalancer.getNrOfLinks() - forceBalancer.getNrOfNodes();
  CHECK(nrOfAddedLinks >= 100);
  // std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;
  // check that all f = 2 have already been removed
  // they have not, since more f = 2 are produced by
  // numInactiveRemoved = forceBalancer.removeBifunctionalCrosslinks(
  //       net, displacements, partitions);
  // CHECK(numInactiveRemoved == 0);

  // run a while to get inactive links
  CHECK_NOTHROW(forceBalancer.runForceRelaxation());
  net = forceBalancer.getNetwork();
  displacements = forceBalancer.getCurrentDisplacements();
  // due to the randomness, it _could_ be one day that actually all strands
  // are active. unlikely, but I can imagine it to be possible.
  numInactiveRemoved =
    forceBalancer.removeInactiveLinks(net, displacements, 0.1);
  CHECK(numInactiveRemoved > 0);
  numInactiveRemoved =
    forceBalancer.unlinkBifunctionalLinks(net, displacements);
  CHECK(numInactiveRemoved > 0);
  CHECK_NOTHROW(forceBalancer.validateNetwork(net, displacements));
}

TEST_CASE(
  "MEHP Force Balance2 runs with non-network",
  "[analysis][MEHPForceBalance2][NonGaussianSpringForceEvaluator][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 runs with non-network\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Balance2 3D case")
  {
    std::string largeInputFile =
      suspectedPath + "xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      pcm::MEHPForceBalance2 forceBalancer2 =
        pcm::MEHPForceBalance2(universe2, 2, false);
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      CHECK_NOTHROW(forceBalancer2.runForceRelaxation());
      CHECK(forceBalancer2.getNrOfIterations() > 1);
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK_THAT(forceBalancer2.getSolubleWeightFraction(),
                 Catch::Matchers::WithinAbs(1., 0.01));
    }
  }
}

TEST_CASE("MEHP Force Balance2 Free chains collapse",
          "[analysis][MEHPForceBalance2]")
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
    atomIds.push_back(i + 1);
    atomTypes.push_back(i % nrOfBeadsPerChain == 0 ? 2 : 1);
    zeroInts.push_back(0);
    if (i > 0) {
      bondFrom.push_back(i);
      bondTo.push_back(i + 1);
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
  CHECK(universe.getNrOfAtoms() == nrOfBeads);
  CHECK(universe.getNrOfBonds() == nrOfBeads - 1);

  // now, check for every force evaluator, that the maximum entropy is when all
  // these beads overlap first, the gaussian spring one
  pcm::MEHPForceBalance2 forceBalancer =
    pcm::MEHPForceBalance2(universe, 2, false);
  CHECK(forceBalancer.getNrOfStrands() == forceBalancer.getNrOfSprings());
  CHECK(forceBalancer.getNrOfStrands() == nrOfBeads / nrOfBeadsPerChain);

  CHECK_NOTHROW(forceBalancer.runForceRelaxation());
  CHECK(forceBalancer.getNrOfIterations() > 0);
  CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  CHECK(forceBalancer.getNrOfActiveStrands() == 0);
  CHECK(forceBalancer.getActiveWeightFraction() == 0.);
  CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
             Catch::Matchers::WithinRel(1.));
  CHECK_THAT(forceBalancer.getDanglingWeightFraction(),
             Catch::Matchers::WithinAbs(0., 1e-9));
  CHECK(forceBalancer.getAverageSpringLength() >= 0.0);
  CHECK(forceBalancer.getAverageSpringLength() <= 3e-6);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
}

TEST_CASE("MEHP Force Balance2 Entanglement Beads Are Removed",
          "[analysis][MEHPForceBalance2]")
{
  // construct one long chain
  size_t nrOfBeads = 30;
  size_t nrOfBeadsPerChain = 3;
  pe::Universe universe =
    pe::Universe(nrOfBeads * 10.0, nrOfBeads * 10.0, nrOfBeads * 10.0);
  // compiler smart enough to reserve, can save some code
  std::vector<double> xPositions;
  std::vector<double> yPositions;
  std::vector<double> zPositions;
  std::vector<long int> atomIds;
  std::vector<int> atomTypes;
  std::vector<int> zeroInts;
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  std::vector<int> bondTypes;
  double offset = 10.0;
  for (int i = 0; i < nrOfBeads; ++i) {
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i + 1);
    atomTypes.push_back(i % nrOfBeadsPerChain == 0 ? 2 : 1);
    zeroInts.push_back(0);
    if (i > 0) {
      bondFrom.push_back(i);
      bondTo.push_back(i + 1);
      bondTypes.push_back(1);
    }
  }
  // actually construct the universe
  universe.addAtoms(atomIds,
                    atomTypes,
                    xPositions,
                    yPositions,
                    zPositions,
                    zeroInts,
                    zeroInts,
                    zeroInts);
  universe.addBonds(bondFrom, bondTo, bondTypes);
  CHECK(universe.getNrOfAtoms() == nrOfBeads);
  CHECK(universe.getNrOfBonds() == bondTo.size());

  // start with the MEHPForceBalance2
  SECTION("Entanglement links")
  {
    pcm::MEHPForceBalance2 forceBalancer = pcm::MEHPForceBalance2(
      universe, 3, 2.0, 0., 0, 0, "asdflseed", 2, false, false, false);
    CHECK(forceBalancer.getNrOfStrands() == forceBalancer.getNrOfSprings());
    CHECK(forceBalancer.getNrOfStrands() == nrOfBeads / nrOfBeadsPerChain + 3);

    CHECK_NOTHROW(forceBalancer.runForceRelaxation(
      pcm::StructureSimplificationMode::ALL_TIM));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getNrOfActiveStrands() == 0);
    CHECK(forceBalancer.getNrOfStrands() == 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }
  SECTION("Entanglement springs")
  {
    pcm::MEHPForceBalance2 forceBalancer = pcm::MEHPForceBalance2(
      universe, 3, 2.0, 0., 0, 0, "asdflseed", 2, false, false, true);
    CHECK(forceBalancer.getNrOfStrands() == forceBalancer.getNrOfSprings());
    CHECK(forceBalancer.getNrOfStrands() == nrOfBeads / nrOfBeadsPerChain + 3);

    CHECK_NOTHROW(forceBalancer.runForceRelaxation(
      pcm::StructureSimplificationMode::ALL_TIM));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getNrOfActiveStrands() == 0);
    CHECK(forceBalancer.getNrOfStrands() == 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }
}

TEST_CASE("MEHP Force Balance2 fully active chains are fully active",
          "[analysis][MEHPForceBalance2][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 fully active chains are "
               "fully active\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Balance2 3D case")
  {
    // perfect diamond network = fully connected =>
    // maximum is at perfect crystal structure -> must be all active.
    std::string inputFile =
      suspectedPath +
      "/structure/"
      "3d-diamond-lattice_10x10x10_a_3_d_0.85_v_0.V-fixed.structure.out";
    REQUIRE(std::filesystem::exists(inputFile));
    std::cout << "Reading file " << inputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);
    std::cout << "Read file " << inputFile << std::endl;

    std::vector<pe::Molecule> molecules = universe.getChainsWithCrosslinker(2);
    for (pe::Molecule mol : molecules) {
      CHECK(mol.getType() == pe::NETWORK_STRAND);
    }

    SECTION("Without entanglements")
    {
      pcm::MEHPForceBalance2 forceBalancer =
        pcm::MEHPForceBalance2(universe, 2, false);
      size_t initialNSprings = forceBalancer.getNrOfStrands();

      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      double initialResidual = forceBalancer.getDisplacementResidualNorm();
      CHECK(std::isfinite(initialResidual));
      REQUIRE_NOTHROW(forceBalancer.runForceRelaxation(
        pcm::StructureSimplificationMode::ALL_TIM));
      CHECK(forceBalancer.getNrOfIterations() > 0);
      CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      CHECK(forceBalancer.getNrOfStrands() == initialNSprings);
      CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                 Catch::Matchers::WithinRel(1.0));
      CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                 Catch::Matchers::WithinAbs(0.0, 1e-9));
      CHECK(initialResidual > forceBalancer.getDisplacementResidualNorm());
    }

    SECTION("With entanglement links")
    {
      pcm::MEHPForceBalance2 forceBalancer = pcm::MEHPForceBalance2(
        universe, 400, 2.0, 0.0, 100, 0.0, "a533d", 2, false, false, false);

      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      double initialResidual = forceBalancer.getDisplacementResidualNorm();
      CHECK(std::isfinite(initialResidual));
      REQUIRE_NOTHROW(forceBalancer.runForceRelaxation(
        pcm::StructureSimplificationMode::ALL_TIM));
      CHECK(forceBalancer.getNrOfIterations() > 0);
      CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                 Catch::Matchers::WithinRel(1.0));
      CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                 Catch::Matchers::WithinAbs(0.0, 1e-9));
      CHECK(initialResidual > forceBalancer.getDisplacementResidualNorm());
    }

    SECTION("With entanglement springs")
    {
      pcm::MEHPForceBalance2 forceBalancer = pcm::MEHPForceBalance2(
        universe, 400, 2.0, 0.0, 100, 0.0, "a533d", 2, false, false, true);

      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      double initialResidual = forceBalancer.getDisplacementResidualNorm();
      CHECK(std::isfinite(initialResidual));
      REQUIRE_NOTHROW(forceBalancer.runForceRelaxation(
        pcm::StructureSimplificationMode::ALL_TIM));
      CHECK(forceBalancer.getNrOfIterations() > 0);
      CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                 Catch::Matchers::WithinRel(1.0));
      CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                 Catch::Matchers::WithinAbs(0.0, 1e-9));
      CHECK(initialResidual > forceBalancer.getDisplacementResidualNorm());
    }
  }
}

TEST_CASE("MEHPForceBalance2 gives approx. same results for entanglement links "
          "& springs in diamond network",
          "[analysis][MEHPForceBalance2]")
{
  std::cout << "Running test \"MEHPForceBalance2 gives approx. same results "
               "for entanglement links & springs in diamond network\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  // perfect diamond network = fully connected =>
  // maximum is at perfect crystal structure -> must be all active.
  std::string inputFile =
    suspectedPath +
    "/structure/3d-diamond-lattice_5x5x5_a_3_d_0.85_imperfect.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  auto entanglements =
    pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
      universe, 250, 2.5, 0., 10, 3, "a533d");

  pcm::MEHPForceBalance2 forceBalanceLinks =
    pcm::MEHPForceBalance2(universe, entanglements, 2, false, false);

  pcm::MEHPForceBalance2 forceBalanceSprings =
    pcm::MEHPForceBalance2(universe, entanglements, 2, false, true);

  CHECK(forceBalanceLinks.getStressTensor().isApprox(
    forceBalanceSprings.getStressTensor()));

  double initialResidual = forceBalanceLinks.getDisplacementResidualNorm();
  CHECK(initialResidual ==
        Catch::Approx(forceBalanceSprings.getDisplacementResidualNorm()));

  forceBalanceLinks.runForceRelaxation();
  forceBalanceSprings.runForceRelaxation();

  CHECK(forceBalanceLinks.getStressTensor().isApprox(
    forceBalanceSprings.getStressTensor()));
  CHECK(forceBalanceLinks.getDisplacementResidualNorm() ==
        Catch::Approx(forceBalanceSprings.getDisplacementResidualNorm()));
  CHECK_THAT(
    forceBalanceLinks.getIdsOfActiveNodes(),
    Catch::Matchers::Equals(forceBalanceSprings.getIdsOfActiveNodes()));
  CHECK(forceBalanceLinks.getNrOfSprings() <
        forceBalanceSprings.getNrOfSprings());

  forceBalanceLinks.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM, 0.01);
  forceBalanceSprings.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM, 0.01);

  CHECK(forceBalanceLinks.getStressTensor().isApprox(
    forceBalanceSprings.getStressTensor()));
  CHECK(forceBalanceLinks.getDisplacementResidualNorm() ==
        Catch::Approx(forceBalanceSprings.getDisplacementResidualNorm()));
  CHECK(forceBalanceLinks.getNrOfSprings() ==
        forceBalanceSprings.getNrOfSprings());
  CHECK_THAT(
    forceBalanceLinks.getIdsOfActiveNodes(),
    Catch::Matchers::Equals(forceBalanceSprings.getIdsOfActiveNodes()));
  CHECK(forceBalanceLinks.getNrOfSprings() <
        forceBalanceSprings.getNrOfSprings());

  forceBalanceLinks.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM, 0.01);
  forceBalanceSprings.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM, 0.01);

  CHECK(forceBalanceLinks.getNrOfSprings() <
        forceBalanceSprings.getNrOfSprings());
  CHECK(forceBalanceLinks.getStressTensor().isApprox(
    forceBalanceSprings.getStressTensor()));
  CHECK(forceBalanceLinks.getDisplacementResidualNorm() ==
        Catch::Approx(forceBalanceSprings.getDisplacementResidualNorm()));
  CHECK_THAT(
    forceBalanceLinks.getIdsOfActiveNodes(),
    Catch::Matchers::Equals(forceBalanceSprings.getIdsOfActiveNodes()));
}

TEST_CASE("MEHPForceBalance2 gives approx. same results for entanglement links "
          "& springs in p = 1 network",
          "[analysis][MEHPForceBalance2][long]")
{
  std::cout << "Running test \"MEHPForceBalance2 gives approx. same results "
               "for entanglement links "
               "& springs in p = 1 network\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  // perfect diamond network = fully connected =>
  // maximum is at perfect crystal structure -> must be all active.
  std::string inputFile =
    suspectedPath +
    "/structure/"
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  auto entanglements =
    pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
      universe, 25, 2.5, 0., 5, 2, "*a533d");

  pcm::MEHPForceBalance2 forceBalanceLinks =
    pcm::MEHPForceBalance2(universe, entanglements, 2, false, false);

  pcm::MEHPForceBalance2 forceBalanceSprings =
    pcm::MEHPForceBalance2(universe, entanglements, 2, false, true);

  CHECK(forceBalanceLinks.getStressTensor().isApprox(
    forceBalanceSprings.getStressTensor()));

  double initialResidual = forceBalanceLinks.getDisplacementResidualNorm();
  CHECK(initialResidual ==
        Catch::Approx(forceBalanceSprings.getDisplacementResidualNorm()));

  forceBalanceLinks.runForceRelaxation();
  forceBalanceSprings.runForceRelaxation();

  CHECK(forceBalanceLinks.getStressTensor().isApprox(
    forceBalanceSprings.getStressTensor()));
  CHECK(forceBalanceLinks.getDisplacementResidualNorm() ==
        Catch::Approx(forceBalanceSprings.getDisplacementResidualNorm()));
  CHECK_THAT(
    forceBalanceLinks.getIdsOfActiveNodes(),
    Catch::Matchers::Equals(forceBalanceSprings.getIdsOfActiveNodes()));

  forceBalanceLinks.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM, 0.01);
  forceBalanceSprings.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM, 0.01);

  CHECK(forceBalanceLinks.getStressTensor().isApprox(
    forceBalanceSprings.getStressTensor()));
  CHECK(forceBalanceLinks.getDisplacementResidualNorm() ==
        Catch::Approx(forceBalanceSprings.getDisplacementResidualNorm()));
  CHECK_THAT(
    forceBalanceLinks.getIdsOfActiveNodes(),
    Catch::Matchers::Equals(forceBalanceSprings.getIdsOfActiveNodes()));
}

TEST_CASE("MEHP Force Balance2 does not collapse",
          "[analysis][MEHPForceBalance2]")
{
  std::cout << "Running test \"MEHP Force Balance2 does not collapse\""
            << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  /**
   * @brief A grid of two rows, each one bead between the two crosslinkers
   *
   */
  universe.addAtoms(
    { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } },
    { { 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1 } },
    { { 0., 2.5, 5, 7.5, 0.1, 2.5, 5, 7.5, -0.1, 5., 0., 5. } },
    // x with slight (0.1) deviation, so we don't start perfect
    { { 0.1, 0., -0.1, 0., 5., 5., 5., 5., 2.5, 2.5, 7.5, 7.5 } },
    // y
    { { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. } },
    // z
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } });
  universe.addBonds(
    { { 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7 } },
    { { 2, 9, 4, 11, 2, 4, 10, 12, 9, 11, 6, 8, 6, 8, 10, 12 } });

  pcm::MEHPForceBalance2 forceBalanceNew =
    pcm::MEHPForceBalance2(universe, 2, true);
  CHECK_NOTHROW(forceBalanceNew.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM));
  CHECK(forceBalanceNew.getNrOfIterations() > 0);
  CHECK(forceBalanceNew.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  CHECK(forceBalanceNew.getNrOfActiveStrands() ==
        forceBalanceNew.getNrOfStrands());
  // compare to what we expect
  CHECK(forceBalanceNew.getNrOfActiveStrands() == 8);

  CHECK(forceBalanceNew.getNrOfActiveStrands() <=
        (forceBalanceNew.getNrOfActiveStrandsInDir(0) +
         forceBalanceNew.getNrOfActiveStrandsInDir(1) +
         forceBalanceNew.getNrOfActiveStrandsInDir(2)));
  CHECK(forceBalanceNew.getNrOfActiveNodes() == 4);
  CHECK(forceBalanceNew.getAverageSpringLength() == Catch::Approx(5.0));
  // forceBalanceNew.setSpringContourLengths(
  //   Eigen::VectorXd::Constant(forceBalanceNew.getNrOfSprings(), 5.));
  // TODO: check this again
  CHECK_THAT(forceBalanceNew.getGammaFactor(1., 100),
             Catch::Matchers::WithinAbs(1.0, 1e-2));
  CHECK_THAT(forceBalanceNew.getPressure(),
             Catch::Matchers::WithinAbs(1. / 30., 1e-3));
  // outputNetwork(forceBalanceNew.getNetwork(),
  //               forceBalanceNew.getCurrentDisplacements(),
  //               forceBalanceNew.getSpringPartitions());
  CHECK_THAT(forceBalanceNew.getSolubleWeightFraction(),
             Catch::Matchers::WithinAbs(0.0, 1e-9));
};

TEST_CASE(
  "MEHP Force Balance2 correctly collapses melts even with random slip-links",
  "[analysis][MEHPForceBalance2]")
{
  std::cout << "Running test \"MEHP Force Balance2 correctly collapses melts "
               "even with random slip-links\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath + "/structure/melt_83_a_100.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPForceBalance2 forceBalancer =
    pcm::MEHPForceBalance2(universe, 1000, 2.0, 0.0, 100, 5, "my_seed_fb12");

  CHECK_NOTHROW(forceBalancer.validateNetwork());
  double initialResidual = forceBalancer.getDisplacementResidualNorm();
  CHECK(std::isfinite(initialResidual));
  CHECK(forceBalancer.getNumExtraAtoms() > 100);
  CHECK(forceBalancer.getNetwork().nrOfSprings >
        forceBalancer.getNetwork().nrOfStrands);
  CHECK_NOTHROW(forceBalancer.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM));
  CHECK(forceBalancer.getNrOfIterations() > 0);
  CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  CHECK(forceBalancer.getNrOfActiveStrands() == 0);
  CHECK(forceBalancer.getNrOfActiveStrandsInDir(0) == 0);
  CHECK(forceBalancer.getNrOfActiveStrandsInDir(1) == 0);
  CHECK(forceBalancer.getNrOfActiveStrandsInDir(2) == 0);
  CHECK(initialResidual > forceBalancer.getDisplacementResidualNorm());
}

TEST_CASE("Particular MEHP Force Balance2 Example",
          "[analysis][MEHPForceBalance2]")
{
  std::cout << "Running test \"Particular MEHP Force Balance2 Example\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath +
    "/structure/"
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure.out-"
    "equilibration_do_crosslink.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPForceBalance2 forceBalancer =
    pcm::MEHPForceBalance2(universe, 2);
  double pressBefore = forceBalancer.getPressure();
  CHECK_NOTHROW(forceBalancer.runForceRelaxation());
  CHECK(pressBefore > forceBalancer.getPressure());
}

TEST_CASE(
  "MEHPForceBalance2 Adding slip-links does not influence other springs",
  "[analysis][MEHPForceBalance2][long]")
{
  std::cout
    << "Running test \"Adding slip-links does not influence other springs\""
    << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath +
    "/structure/"
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out";
  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  std::cout << "Read file " << inputFile << std::endl;

  // generate the same slip-links twice,
  // once for each assumption
  pcm::MEHPForceBalance2 forceBalancerWithoutEntanglements =
    pcm::MEHPForceBalance2(universe);
  pcm::MEHPForceBalance2 forceBalancerOldSamplingSmall =
    pcm::MEHPForceBalance2(universe, 1000, 6.0, 0.0, 900);
  Eigen::VectorXd springVectorsWithoutEntanglements =
    forceBalancerWithoutEntanglements.evaluateStrandVectors(
      forceBalancerWithoutEntanglements.getNetwork(),
      forceBalancerWithoutEntanglements.getCurrentDisplacements());

  Eigen::VectorXd springVectorsWithEntanglements =
    forceBalancerOldSamplingSmall.evaluateStrandVectors(
      forceBalancerOldSamplingSmall.getNetwork(),
      forceBalancerOldSamplingSmall.getCurrentDisplacements());

  CHECK(
    springVectorsWithoutEntanglements.isApprox(springVectorsWithEntanglements));

  // and the same for sampling method 2
  pcm::MEHPForceBalance2 forceBalancer2 =
    pcm::MEHPForceBalance2(universe, 1000, 6.0, 0.0, 900, 3.0, "53467829");

  pcm::MEHPForceBalance2 forceBalancer2Without =
    pcm::MEHPForceBalance2(universe, 0, 6.0, 0.0, 0, 3.0, "53467829");

  CHECK(forceBalancer2.getNrOfSprings() >
        forceBalancer2Without.getNrOfSprings());

  CHECK(forceBalancer2
          .evaluateStrandVectors(forceBalancer2.getNetwork(),
                                 forceBalancer2.getCurrentDisplacements())
          .isApprox(forceBalancer2Without.evaluateStrandVectors(
            forceBalancer2Without.getNetwork(),
            forceBalancer2Without.getCurrentDisplacements())));
}

TEST_CASE(
  "MEHP Force Balance2 Entanglement Links are just like Entanglement Springs",
  "[analysis][MEHPForceBalance2][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 Entanglement Links are just "
               "like Entanglement Springs\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  // a structure with lots of dangling things that can and will be entangled,
  // yet the entanglements removed
  std::string inputFile =
    suspectedPath +
    "/structure/mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0."
    "0114_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  auto masses = universe.getMasses();
  std::cout << "Read file " << inputFile << std::endl;

  // sample entanglements
  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
        universe, 4368, 4., 0., 190, 0, "af1346lhkdsaöf123", 2, true);

  // initialize the two force balance
  pcm::MEHPForceBalance2 forceBalancerEntanglementSprings =
    pcm::MEHPForceBalance2(universe, entanglements, 2, false, true);

  pcm::MEHPForceBalance2 forceBalancerEntanglementLinks =
    pcm::MEHPForceBalance2(universe, entanglements, 2, false, false);

  REQUIRE(forceBalancerEntanglementSprings.getNrOfStrands() >
          forceBalancerEntanglementLinks.getNrOfStrands());
  REQUIRE(forceBalancerEntanglementLinks.getNrOfStrands() == 464);
  REQUIRE(forceBalancerEntanglementSprings.getNrOfActiveStrands() >
          forceBalancerEntanglementLinks.getNrOfActiveStrands());

  std::vector<long int> activeNodes0_1 =
    forceBalancerEntanglementSprings.getIdsOfActiveNodes();
  std::sort(activeNodes0_1.begin(), activeNodes0_1.end());
  std::vector<long int> activeNodes0_2 =
    forceBalancerEntanglementLinks.getIdsOfActiveNodes();
  std::sort(activeNodes0_2.begin(), activeNodes0_2.end());
  REQUIRE(activeNodes0_1 == activeNodes0_2);

  SECTION("With removal of inactive springs")
  {
    forceBalancerEntanglementSprings.runForceRelaxation(
      pcm::StructureSimplificationMode::ALL_TIM);
    forceBalancerEntanglementLinks.runForceRelaxation(
      pcm::StructureSimplificationMode::ALL_TIM);

    CHECK(forceBalancerEntanglementSprings.getNrOfStrands() >
          forceBalancerEntanglementLinks.getNrOfStrands());
    std::vector<long int> activeNodes1 =
      forceBalancerEntanglementSprings.getIdsOfActiveNodes();
    std::sort(activeNodes1.begin(), activeNodes1.end());
    std::vector<long int> activeNodes2 =
      forceBalancerEntanglementLinks.getIdsOfActiveNodes();
    std::sort(activeNodes2.begin(), activeNodes2.end());
    CHECK(activeNodes1 == activeNodes2);

    CHECK_THAT(
      forceBalancerEntanglementSprings.getSolubleWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getSolubleWeightFraction(), 0.005));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getActiveWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getActiveWeightFraction(), 0.005));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getDanglingWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getDanglingWeightFraction(), 0.005));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getStressTensor().trace(),
      Catch::Matchers::WithinRel(
        forceBalancerEntanglementLinks.getStressTensor().trace(), 0.1));
    CHECK(forceBalancerEntanglementSprings.getStressTensor().trace() <
          forceBalancerEntanglementLinks.getStressTensor().trace());
    CHECK(forceBalancerEntanglementSprings.getGammaFactors(1.).sum() <
          forceBalancerEntanglementLinks.getGammaFactors(1.).sum());
  }

  SECTION("Without simplification of the structure")
  {
    // need to disable slipping to test entanglement links
    forceBalancerEntanglementSprings.runForceRelaxation(
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION);
    forceBalancerEntanglementLinks.runForceRelaxation(
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION);

    CHECK(forceBalancerEntanglementSprings.getNrOfStrands() >
          forceBalancerEntanglementLinks.getNrOfStrands());
    std::vector<long int> activeNodes1 =
      forceBalancerEntanglementSprings.getIdsOfActiveNodes();
    std::sort(activeNodes1.begin(), activeNodes1.end());
    std::vector<long int> activeNodes2 =
      forceBalancerEntanglementLinks.getIdsOfActiveNodes();
    std::sort(activeNodes2.begin(), activeNodes2.end());
    CHECK(activeNodes1 == activeNodes2);

    CHECK_THAT(
      forceBalancerEntanglementSprings.getSolubleWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getSolubleWeightFraction(), 0.0025));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getActiveWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getActiveWeightFraction(), 0.005));
    CHECK_THAT(
      forceBalancerEntanglementSprings.getDanglingWeightFraction(),
      Catch::Matchers::WithinAbs(
        forceBalancerEntanglementLinks.getDanglingWeightFraction(), 0.005));
    CHECK(forceBalancerEntanglementSprings.getStressTensor().trace() <
          forceBalancerEntanglementLinks.getStressTensor().trace());
    CHECK(forceBalancerEntanglementSprings.getGammaFactors(1.).sum() <
          forceBalancerEntanglementLinks.getGammaFactors(1.).sum());
  }
}

TEST_CASE(
  "MEHP Force Balance2 adding Entanglements vs. Phantom behaves as expected",
  "[analysis][MEHPForceBalance2][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 adding Entanglements vs. "
               "Phantom behaves as expected\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  // a structure with lots of dangling things that can and will be entangled,
  // yet the entanglements removed
  std::string inputFile =
    suspectedPath +
    "/structure/mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0."
    "0114_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  auto masses = universe.getMasses();
  std::cout << "Read file " << inputFile << std::endl;

  // sample entanglements
  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
        universe, 4368, 4., 0., 190, 0, "af1346lhkdsaöf123", 2, true);

  pcm::MEHPForceBalance2 forceBalancerEntanglements =
    pcm::MEHPForceBalance2(universe, entanglements, 2, false);

  pcm::MEHPForceBalance2 forceBalancerPhantom =
    pcm::MEHPForceBalance2(universe, 2, false);

  double initialResidualE =
    forceBalancerEntanglements.getDisplacementResidualNorm();
  double initialResidualP = forceBalancerPhantom.getDisplacementResidualNorm();

  CHECK(forceBalancerEntanglements.getNrOfSprings() >
        forceBalancerPhantom.getNrOfSprings());
  CHECK(forceBalancerEntanglements.getNrOfActiveStrands() >
        forceBalancerPhantom.getNrOfActiveStrands());
  CHECK(forceBalancerEntanglements.getSolubleWeightFraction() <
        forceBalancerPhantom.getSolubleWeightFraction());

  forceBalancerEntanglements.runForceRelaxation();
  forceBalancerPhantom.runForceRelaxation();

  CHECK(forceBalancerEntanglements.getNrOfSprings() >
        forceBalancerPhantom.getNrOfSprings());
  CHECK(forceBalancerEntanglements.getNrOfActiveStrands() >
        forceBalancerPhantom.getNrOfActiveStrands());
  CHECK(forceBalancerEntanglements.getSolubleWeightFraction() <
        forceBalancerPhantom.getSolubleWeightFraction());

  forceBalancerEntanglements.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM);
  forceBalancerPhantom.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM);

  CHECK(forceBalancerEntanglements.getNrOfStrands() >
        forceBalancerPhantom.getNrOfStrands());
  CHECK(forceBalancerEntanglements.getNrOfActiveStrands() >
        forceBalancerPhantom.getNrOfActiveStrands());
  CHECK(forceBalancerEntanglements.getSolubleWeightFraction() <
        forceBalancerPhantom.getSolubleWeightFraction());
};

TEST_CASE("MEHPForceBalance2 phantom with and without removal is the same",
          "[analysis][MEHPForceBalance2]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  // a structure with lots of dangling things that can and will be entangled,
  // yet the entanglements removed
  std::string inputFile =
    suspectedPath +
    "/structure/mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0."
    "0114_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  auto masses = universe.getMasses();
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPForceBalance2 forceBalancerPhantom =
    pcm::MEHPForceBalance2(universe, 2, false);
  pcm::MEHPForceBalance2 forceBalancerPhantomRem =
    pcm::MEHPForceBalance2(universe, 2, false);

  CHECK(forceBalancerPhantom.getNrOfSprings() ==
        forceBalancerPhantomRem.getNrOfSprings());
  CHECK_THAT(
    forceBalancerPhantom.getResidual(),
    Catch::Matchers::WithinRel(forceBalancerPhantomRem.getResidual(), 0.001));
  CHECK_THAT(forceBalancerPhantom.getSolubleWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getSolubleWeightFraction(), 0.001));
  CHECK_THAT(
    forceBalancerPhantom.getGamma(),
    Catch::Matchers::WithinRel(forceBalancerPhantomRem.getGamma(), 0.001));
  CHECK_THAT(forceBalancerPhantom.getDanglingWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getDanglingWeightFraction(), 0.001));

  // run with and without simplification
  forceBalancerPhantom.runForceRelaxation(
    pcm::StructureSimplificationMode::NO_SIMPLIFICATION);
  forceBalancerPhantomRem.runForceRelaxation(
    pcm::StructureSimplificationMode::ALL_TIM);

  // compare results
  CHECK(forceBalancerPhantom.getNrOfSprings() <
        forceBalancerPhantomRem.getNrOfSprings());
  CHECK_THAT(
    forceBalancerPhantom.getResidual(),
    Catch::Matchers::WithinRel(forceBalancerPhantomRem.getResidual(), 0.001));
  CHECK_THAT(forceBalancerPhantom.getSolubleWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getSolubleWeightFraction(), 0.001));
  CHECK_THAT(
    forceBalancerPhantom.getGamma(),
    Catch::Matchers::WithinRel(forceBalancerPhantomRem.getGamma(), 0.001));
  CHECK_THAT(forceBalancerPhantom.getDanglingWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getDanglingWeightFraction(), 0.001));
}

TEST_CASE("Temporary force balance 2 test case",
          "[analysis][MEHPForceBalance2][long]")
{
  std::cout << "Running test \"Temporary force balance test case\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  // a structure with lots of dangling things that can and will be entangled,
  // yet the entanglements removed
  std::string inputFile =
    suspectedPath +
    "/tmp/mc_own-si_pdms_crosslinked_melt_2590_a_221_2410_monoa_271_"
    "r_1.44_wsol_0.278_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  auto masses = universe.getMasses();
  std::cout << "Read file " << inputFile << std::endl;

  // sample entanglements
  pcm::MEHPForceBalance2 forceBalancerEntanglements =
    pcm::MEHPForceBalance2(universe, 26832, 3.25, 0., 26000, 0, "", 2, false);

  std::cout << "Sampled entanglements " << std::endl;

  CHECK_NOTHROW(forceBalancerEntanglements.runForceRelaxation(
    pcm::StructureSimplificationMode::NO_SIMPLIFICATION));
};

TEST_CASE("MEHPFB2 Basic conversion test",
          "[analysis][MEHPForceBalance2][MEHPForceBalance]")
{
  std::cout << "Running test \"MEHPFB2 Basic conversion test\"" << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  std::vector<double> coordsX = {
    5.0,                  // 0
    7.,  9., 10., 1., 3., // 1-5: horizontal line
    5.,  5., -5., 5., 5., // 6-10: vertical line
    7.,  9., 10., 1., 3., // 11-15: diagonal line
  };
  std::vector<double> coordsY = {
    5.0,                   // 0
    5.,  5., 5.,  -5., 5., // 1-5: horizontal line in x-plane
    7.,  9., 10., 1.,  3., // 6-10: horizontal line,
    7.,  9., 10., 1.,  3., // 11-15: diagonal line
  };
  std::vector<double> coordsZ = {
    5.,                    // 0
    5., 15., 5.,  -5., 5., // 1-5: horizontal line in x-plane
    5., 5.,  5.,  15., 5., // 6-10: horizontal line in y-plane,
    7., 9.,  10., 1.,  3., // 11-15: diagonal line
  };
  REQUIRE(coordsX.size() == coordsY.size());
  REQUIRE(coordsZ.size() == coordsY.size());
  std::vector<long int> ids = {};
  ids.reserve(coordsX.size());
  for (long int i = 0; i < coordsX.size(); ++i) {
    ids.push_back(i);
  }
  std::vector<int> types =
    pylimer_tools::utils::initializeWithValue(ids.size(), 1);
  types[0] = 2;
  types[2] = 2;
  universe.addAtoms(ids,
                    types,
                    coordsX,
                    coordsY,
                    coordsZ,
                    pylimer_tools::utils::initializeWithValue(ids.size(), 0),
                    pylimer_tools::utils::initializeWithValue(ids.size(), 0),
                    pylimer_tools::utils::initializeWithValue(ids.size(), 0));
  // add the straight lines going through the box,
  // as well as add a primary loop at vertex 2, since otherwise,
  // there is only 1 bead and no displacement needed
  std::vector<long int> bondFrom = { 0, 1,  2, 3,  4,  5,  0,  6,  7, 8,
                                     9, 10, 0, 11, 12, 13, 14, 15, 2 };
  std::vector<long int> bondTo = { 1,  2, 3,  4,  5,  0,  6,  7, 8, 9,
                                   10, 0, 11, 12, 13, 14, 15, 0, 2 };
  universe.addBonds(bondFrom, bondTo);

  CHECK(universe.getNrOfBonds() == bondFrom.size());
  CHECK(universe.getNrOfAtoms() == ids.size());

  pylimer_tools::sim::mehp::MEHPForceBalance fb1 =
    pylimer_tools::sim::mehp::MEHPForceBalance(universe);
  fb1.configAssumeBoxLargeEnough(false);

  pylimer_tools::sim::mehp::MEHPForceBalance2 fb2 =
    pylimer_tools::sim::mehp::MEHPForceBalance2(universe);

  pylimer_tools::sim::mehp::ForceBalance2Network net = fb2.getNetwork();

  CHECK(net.nrOfLinks == 2);
  CHECK(net.nrOfSprings == 5);
  CHECK(net.nrOfStrands == 5);
  CHECK(net.nrOfNodes == 2);
  CHECK_FALSE(net.springBoxOffset.isZero());

  CHECK(fb1.getNrOfSprings() == 5);

  CHECK_THAT(fb1.getGammaFactors(1.0).mean(),
             Catch::Matchers::WithinAbs(fb2.getGammaFactors(1.0).mean(), 1e-6));
  CHECK_THAT(fb1.getResidual(),
             Catch::Matchers::WithinRel(fb2.getResidual(), 1e-6));

  fb1.runForceRelaxation();
  fb2.runForceRelaxation();

  CHECK_THAT(fb1.getResidual(),
             Catch::Matchers::WithinAbs(fb2.getResidual(), 1e-6));
  CHECK_THAT(fb1.getGammaFactors(1.0).mean(),
             Catch::Matchers::WithinAbs(fb2.getGammaFactors(1.0).mean(), 1e-6));

  CHECK(fb2.getSolubleWeightFraction() == 0.);
  CHECK(fb2.getNrOfActiveStrands() == 4); // the primary loop is collapsed
  CHECK(fb1.getSolubleWeightFraction() == 0.);
  CHECK(fb1.getNrOfActiveSprings() == 4); // the primary loop is collapsed

  // some other network validations
  REQUIRE(net.oldAtomIds.size() == 2);
  // NOTE: the following checks are a bit too much checking the implementation
  // rather than just the functionality.
  CHECK(net.oldAtomIds[0] == 0);
  CHECK(net.oldAtomIds[1] == 2);
  CHECK(net.oldAtomTypes[0] == 2);
  CHECK(net.oldAtomTypes[1] == 2);

  CHECK(net.coordinates.size() == 6);
  CHECK(net.coordinates.segment<3>(0, 3).isApprox(Eigen::Vector3d(5., 5., 5.)));
  CHECK(net.coordinates.segment<3>(3, 3).isApprox(Eigen::Vector3d(9., 5., 5.)));
  CHECK(net.springIndexA.size() == 5);
  // first horizontal line
  CHECK(net.springIndexA[0] == 0);
  CHECK(net.springIndexB[0] == 1);
  CHECK(net.springBoxOffset.segment<3>(0 * 3, 3).isZero());
  CHECK(net.springContourLength[0] == 2);
  // primary loop without offset
  CHECK(net.springIndexA[1] == 1);
  CHECK(net.springIndexB[1] == 1);
  CHECK(net.springBoxOffset.segment<3>(1 * 3, 3).isZero());
  CHECK(net.springContourLength[1] == 1);
  // closing first line
  CHECK(net.springIndexA[2] == 1);
  CHECK(net.springIndexB[2] == 0);
  CHECK(net.springBoxOffset.segment<3>(2 * 3, 3).isApprox(
    Eigen::Vector3d(10., 0., 0.)));
  CHECK(net.springContourLength[2] == 4);
  // second line
  CHECK(net.springIndexA[3] == 0);
  CHECK(net.springIndexB[3] == 0);
  CHECK(net.springBoxOffset.segment<3>(3 * 3, 3).isApprox(
    Eigen::Vector3d(0., 10., 0.)));
  CHECK(net.springContourLength[3] == 6);
  // third line
  CHECK(net.springIndexA[4] == 0);
  CHECK(net.springIndexB[4] == 0);
  CHECK(net.springBoxOffset.segment<3>(4 * 3, 3).isApprox(
    Eigen::Vector3d(10., 10., 10.)));
  CHECK(net.springContourLength[4] == 6);
}

TEST_CASE("All MEHP Force Balance2 vs. Force Relaxation Phantom Comparisons",
          "[analysis][MEHPForceBalance2][MEHPForceRelaxation][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 vs. Force Relaxation "
               "Phantom Comparisons\""
            << std::endl;

  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::vector<std::string> files = {
    "3d-diamond-lattice_10x10x10_a_3_d_0.85_imperfect.structure.out",
    "3d-diamond-lattice_10x10x10_a_3_d_0.85_v_0.V-fixed.structure.out",
    "3d-diamond-lattice_3x3x3_a_23_d_3_v_0.structure.out",
    "3d-diamond-lattice_5x5x5_a_3_d_0.85_imperfect.structure.out",
    "3d-diamond-lattice_5x5x5_a_3_d_0.85_v_0.V-fixed.structure.out",
    "crosslinked_M10000_N39_p_0.9.out",
    "crosslinked_p_0.98_melt_100_a_3_50_xlinks_v_14.converted.structure.out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_0.98_melt_100_a_38_50_xlinks_v_22.structure.out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out",
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure.out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_1_1_melt_100_a_3_50_xlinks_v_1.V-fixed.structure.out-finish_"
    "crosslinking.structure.out",
    "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0.333_2d_t_7500001."
    "structure.out",
    "mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0.0114_f_4_v_1."
    "structure.out",
    "melt_213_a_47_106_xlinks_v_1.structure.out",
    "melt_83_a_100.structure.out",
    "network_100_a_46.structure.out",
    "network_p_1_100_a_38_50_xlinks.structure.out",
    "square_lattice_2x2_a_5.2d.structure.out"
  };

  size_t nFilesFound = 0;
  for (const std::string& file : files) {
    pe::UniverseSequence universeSeq = pe::UniverseSequence();
    REQUIRE(universeSeq.getLength() == 0);
    std::string inputFile = suspectedPath + "/structure/" + file;
    if (!std::filesystem::exists(inputFile)) {
      std::cerr << "File not found: " << inputFile << std::endl;
      continue;
    }

    nFilesFound += 1;
    std::cout << "Processing file: " << file << std::endl;

    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);

    pcm::MEHPForceBalance2 forceBalance2 = pcm::MEHPForceBalance2(universe);

    pcm::MEHPForceRelaxation forceRelaxation =
      pcm::MEHPForceRelaxation(universe);
    forceRelaxation.configAssumeBoxLargeEnough(false);

    CHECK_THAT(forceBalance2.getGamma(),
               Catch::Matchers::WithinRel(forceRelaxation.getGamma(), 1e-3));
    CHECK_THAT(forceBalance2.getStressTensor().trace(),
               Catch::Matchers::WithinRel(
                 forceRelaxation.getStressTensor().trace(), 1e-3));
    CHECK_THAT(forceBalance2.getResidual(),
               Catch::Matchers::WithinRel(forceRelaxation.getResidual(), 1e-3));
    CHECK_THAT(forceBalance2.getSolubleWeightFraction(),
               Catch::Matchers::WithinRel(
                 forceRelaxation.getSolubleWeightFraction(), 1e-3));
    CHECK_THAT(forceBalance2.getDanglingWeightFraction(),
               Catch::Matchers::WithinRel(
                 forceRelaxation.getDanglingWeightFraction(), 1e-3));
    // CHECK_THAT(forceBalance2.getActiveWeightFraction(),
    // Catch::Matchers::WithinRel(forceRelaxation.getActiveWeightFraction(),
    // 1e-3));

    auto start_fb2 = std::chrono::high_resolution_clock::now();
    forceBalance2.runForceRelaxation();
    auto end_fb2 = std::chrono::high_resolution_clock::now();
    auto duration_fb2 = std::chrono::duration_cast<std::chrono::microseconds>(
      end_fb2 - start_fb2);
    std::cout << "Time of Force Balance 2: "
              << std::duration_to_string(duration_fb2) << " " << std::endl;

    auto start_fr = std::chrono::high_resolution_clock::now();
    while (forceRelaxation.suggestsRerun()) {
      forceRelaxation.runForceRelaxation();
    }
    auto end_fr = std::chrono::high_resolution_clock::now();
    auto duration_fr =
      std::chrono::duration_cast<std::chrono::microseconds>(end_fr - start_fr);
    std::cout << "Time of Force Relaxation: "
              << std::duration_to_string(duration_fr) << " " << std::endl;

    CHECK_THAT(forceBalance2.getGamma(),
               Catch::Matchers::WithinRel(forceRelaxation.getGamma(), 1e-3));
    CHECK_THAT(forceBalance2.getStressTensor().trace(),
               Catch::Matchers::WithinRel(
                 forceRelaxation.getStressTensor().trace(), 1e-3));
    CHECK_THAT(forceBalance2.getResidual(),
               Catch::Matchers::WithinRel(forceRelaxation.getResidual(), 1e-3));
    CHECK_THAT(forceBalance2.getSolubleWeightFraction(),
               Catch::Matchers::WithinRel(
                 forceRelaxation.getSolubleWeightFraction(), 1e-3));
    CHECK_THAT(forceBalance2.getDanglingWeightFraction(),
               Catch::Matchers::WithinRel(
                 forceRelaxation.getDanglingWeightFraction(), 1e-3));
    // CHECK_THAT(forceBalance2.getActiveWeightFraction(),
    // Catch::Matchers::WithinRel(forceRelaxation.getActiveWeightFraction(),
    // 1e-3));
  }

  REQUIRE(static_cast<double>(nFilesFound) >
          static_cast<double>(files.size()) * 0.75);
}

TEST_CASE("All MEHP Force Balance 1 vs. 2 Comparisons with Entanglements and "
          "Simplification",
          "[analysis][MEHPForceBalance2][MEHPForceBalance][long]")
{
  std::cout << "Running test \"All MEHP Force Balance 1 vs. 2 Comparisons with "
               "Entanglements and Simplification\""
            << std::endl;

  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::vector<std::string> files = {
    "3d-diamond-lattice_10x10x10_a_3_d_0.85_imperfect.structure.out",
    "3d-diamond-lattice_10x10x10_a_3_d_0.85_v_0.V-fixed.structure.out",
    "3d-diamond-lattice_3x3x3_a_23_d_3_v_0.structure.out",
    "3d-diamond-lattice_5x5x5_a_3_d_0.85_imperfect.structure.out",
    "3d-diamond-lattice_5x5x5_a_3_d_0.85_v_0.V-fixed.structure.out",
    "crosslinked_M10000_N39_p_0.9.out",
    "crosslinked_p_0.98_melt_100_a_3_50_xlinks_v_14.converted.structure.out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_0.98_melt_100_a_38_50_xlinks_v_22.structure.out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_0.99145_0.99145_melt_10000_a_3_5000_xlinks_v_1.V-fixed."
    "structure.out-equilibration_do_crosslink.structure.out",
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure.out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_1_1_melt_100_a_3_50_xlinks_v_1.V-fixed.structure.out-finish_"
    "crosslinking.structure.out",
    "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0.333_2d_t_7500001."
    "structure.out",
    "mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0.0114_f_4_v_1."
    "structure.out",
    "melt_213_a_47_106_xlinks_v_1.structure.out",
    "melt_83_a_100.structure.out",
    "network_100_a_46.structure.out",
    "network_p_1_100_a_38_50_xlinks.structure.out",
    "square_lattice_2x2_a_5.2d.structure.out"
  };

  size_t nFilesFound = 0;
  for (const std::string& file : files) {
    pe::UniverseSequence universeSeq = pe::UniverseSequence();
    REQUIRE(universeSeq.getLength() == 0);
    std::string inputFile = suspectedPath + "/structure/" + file;
    if (!std::filesystem::exists(inputFile)) {
      std::cerr << "File not found: " << inputFile << std::endl;
      continue;
    }

    nFilesFound += 1;
    std::cout << "Processing file: " << inputFile << std::endl;

    universeSeq.initializeFromDataSequence({ { inputFile } });
    pe::Universe universe = universeSeq.atIndex(0);

    // use the same entanglements for both force balances
    auto entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
        universe,
        0.01 * universe.getNrOfAtoms(),
        2.5,
        0.,
        0,
        2,
        "t35ts33d",
        2,
        true,
        true);

    pcm::MEHPForceBalance2 forceBalance2 =
      pcm::MEHPForceBalance2(universe, entanglements);

    pcm::MEHPForceBalance forceBalance =
      pcm::MEHPForceBalance::constructWithSlipLinks(universe, entanglements);
    forceBalance.configAssumeBoxLargeEnough(false);

    // validate initial values
    CHECK_THAT(forceBalance2.getGamma(),
               Catch::Matchers::WithinRel(forceBalance.getGamma(), 1e-3));
    CHECK_THAT(
      forceBalance2.getStressTensor().trace(),
      Catch::Matchers::WithinRel(forceBalance.getStressTensor().trace(), 1e-3));
    CHECK_THAT(forceBalance2.getResidual(),
               Catch::Matchers::WithinRel(forceBalance.getResidual(), 1e-3));
    CHECK_THAT(forceBalance2.getSolubleWeightFraction(),
               Catch::Matchers::WithinRel(
                 forceBalance.getSolubleWeightFraction(), 1e-3));
    CHECK_THAT(forceBalance2.getDanglingWeightFraction(),
               Catch::Matchers::WithinRel(
                 forceBalance.getDanglingWeightFraction(), 1e-3));
    CHECK_THAT(
      forceBalance2.getActiveWeightFraction(),
      Catch::Matchers::WithinRel(forceBalance.getActiveWeightFraction(), 1e-3));

    auto start_fr = std::chrono::high_resolution_clock::now();
    forceBalance.runForceRelaxation(
      5000, 1e-9, -1., pcm::StructureSimplificationMode::ALL_TIM);
    auto end_fr = std::chrono::high_resolution_clock::now();
    auto duration_fr =
      std::chrono::duration_cast<std::chrono::microseconds>(end_fr - start_fr);
    std::cout << "Time of Force Balance 1: "
              << std::duration_to_string(duration_fr) << " " << std::endl;

    auto start_fb2 = std::chrono::high_resolution_clock::now();
    forceBalance2.runForceRelaxation(pcm::StructureSimplificationMode::ALL_TIM);
    auto end_fb2 = std::chrono::high_resolution_clock::now();
    auto duration_fb2 = std::chrono::duration_cast<std::chrono::microseconds>(
      end_fb2 - start_fb2);
    std::cout << "Time of Force Balance 2: "
              << std::duration_to_string(duration_fb2) << " " << std::endl;

    if (forceBalance.getExitReason() != pcm::ExitReason::MAX_STEPS) {
      CHECK_THAT(forceBalance2.getGamma(),
                 Catch::Matchers::WithinRel(forceBalance.getGamma(), 1e-3));
      CHECK_THAT(forceBalance2.getStressTensor().trace(),
                 Catch::Matchers::WithinRel(
                   forceBalance.getStressTensor().trace(), 1e-3));
      CHECK_THAT(forceBalance2.getResidual(),
                 Catch::Matchers::WithinRel(forceBalance.getResidual(), 1e-3));
      CHECK_THAT(forceBalance2.getSolubleWeightFraction(),
                 Catch::Matchers::WithinRel(
                   forceBalance.getSolubleWeightFraction(), 1e-3));
      CHECK_THAT(forceBalance2.getDanglingWeightFraction(),
                 Catch::Matchers::WithinRel(
                   forceBalance.getDanglingWeightFraction(), 1e-3));
      CHECK_THAT(forceBalance2.getActiveWeightFraction(),
                 Catch::Matchers::WithinRel(
                   forceBalance.getActiveWeightFraction(), 1e-3));
    } else {
      std::cerr << "Force Balance 1 did not converge for file: " << inputFile
                << std::endl;
    }
  }

  REQUIRE(static_cast<double>(nFilesFound) >
          static_cast<double>(files.size()) * 0.75);
}