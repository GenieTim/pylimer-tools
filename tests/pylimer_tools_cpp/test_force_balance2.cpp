#include "../../src/pylimer_tools_cpp/entities/Universe.h"
#include "../../src/pylimer_tools_cpp/entities/UniverseSequence.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPForceEvaluator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPLinearSystemSolver.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPNetwork.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPNonlinearOptimizationSolver.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPSimulator.h"
#include "../../src/pylimer_tools_cpp/sim/MEHPSlipLinkRelaxationSolver.h"
#include "../../src/pylimer_tools_cpp/topo/EntanglementDetector.h"
#include "../../src/pylimer_tools_cpp/utils/ExtraEigenTypes.h"
#include "../../src/pylimer_tools_cpp/utils/MCUniverseGenerator.h"
#ifdef USE_METAL
#include "../../src/pylimer_tools_cpp/sim/MetalForceBalanceSolver.h"
#endif
#include "../../src/pylimer_tools_cpp/utils/ExtraEigenSolvers.h"
#include <algorithm>
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
#include <memory>
#include <numeric>
#include <random>
#include <ranges>
#include <stdexcept>
#include <vector>

#include <unistd.h>

namespace pe = pylimer_tools::entities;
namespace pu = pylimer_tools::utils;
namespace pcm = pylimer_tools::sim::mehp;

#ifndef PYLIMER_TEST_FIXTURES_DIR
#define PYLIMER_TEST_FIXTURES_DIR "../pylimer_tools/tests/fixtures"
#endif

static std::vector<long int>
getAtomIdsOfActiveNodes(const pcm::MEHPSimulator& sim)
{
  const pcm::MEHPNetwork& net = sim.getNetwork();
  std::vector<long int> atomIds;
  for (const int linkIdx : sim.getIndicesOfActiveNodes()) {
    const Eigen::Index idx = static_cast<Eigen::Index>(linkIdx);
    if (!net.linkIsEntanglement[idx]) {
      atomIds.push_back(net.oldAtomIds[idx]);
    }
  }
  return atomIds;
}

static size_t
countActiveStrandsInDir(const pcm::MEHPSimulator& sim, const int dir)
{
  return static_cast<size_t>(
    sim.getNetwork()
      .findActiveStrandsInDir(
        sim.getCurrentDisplacements(), dir, 1e-3, sim.getIs2D())
      .count());
}

static Eigen::ArrayXi
countActiveSpringsConnected(const pcm::MEHPSimulator& sim)
{
  const pcm::MEHPNetwork& net = sim.getNetwork();
  const Eigen::ArrayXb active =
    net.findActiveSprings(sim.getCurrentDisplacements(), 1e-3, sim.getIs2D());
  Eigen::ArrayXi counts =
    Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(net.nrOfLinks));
  for (Eigen::Index s = 0; s < static_cast<Eigen::Index>(net.nrOfSprings);
       ++s) {
    if (active[s]) {
      counts[net.springIndexA[s]] += 1;
      counts[net.springIndexB[s]] += 1;
    }
  }
  return counts;
}

static Eigen::ArrayXi
countActiveStrandsConnected(const pcm::MEHPSimulator& sim)
{
  const pcm::MEHPNetwork& net = sim.getNetwork();
  const Eigen::ArrayXb activeStrands =
    net.findActiveStrands(sim.getCurrentDisplacements(), 1e-3, sim.getIs2D());
  Eigen::ArrayXi counts =
    Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(net.nrOfLinks));
  for (size_t s = 0; s < net.nrOfStrands; ++s) {
    if (activeStrands[static_cast<Eigen::Index>(s)]) {
      counts[static_cast<Eigen::Index>(net.linkIndicesOfStrand[s].front())] +=
        1;
      counts[static_cast<Eigen::Index>(net.linkIndicesOfStrand[s].back())] += 1;
    }
  }
  return counts;
}

TEST_CASE("MEHP Force Balance2 runs", "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 runs\"" << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Balance2 3D case")
  {
    std::string largeInputFile =
      suspectedPath +
      "/structure/xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out";
    REQUIRE(std::filesystem::exists(largeInputFile));
    std::cout << "Reading file " << largeInputFile << std::endl;
    universeSeq.initializeFromDataSequence({ { largeInputFile } });
    pe::Universe universe2 = universeSeq.atIndex(0);
    std::cout << "Read file " << largeInputFile << std::endl;

    double nrOfChains = 1.e4;
    CHECK(static_cast<double>(universe2.getMolecules(2).size()) ==
          Catch::Approx(nrOfChains));
    pcm::MEHPSimulator forceBalancer2 =
      pcm::MEHPSimulator(universe2, {}, 2, pcm::EntanglementMode::MERGED_LINKS);
    size_t nStrandsFb = forceBalancer2.getNrOfStrands();
    CHECK(nStrandsFb <= nrOfChains);
    size_t nSpringsFb = forceBalancer2.getNrOfSprings();
    CHECK(nSpringsFb <= nrOfChains);

    Eigen::Matrix3d stressTensor1 = forceBalancer2.getStressTensor();
    std::vector<size_t> allLinks(forceBalancer2.getNrOfLinks());
    std::iota(allLinks.begin(), allLinks.end(), 0);
    Eigen::Matrix3d stressTensor2 =
      forceBalancer2.getNetwork().evaluateStressTensorForLinks(
        allLinks, forceBalancer2.getCurrentDisplacements(), 1.0, false, -1.0);
    for (Eigen::Index i = 0; i < 3; ++i) {
      for (Eigen::Index j = 0; j < 3; ++j) {
        CHECK(stressTensor1(i, j) == Catch::Approx(stressTensor2(i, j)));
      }
    }

    SECTION("Actual balance results in correct phantom results")
    {
      std::cout << "Doing phantom force balance" << std::endl;
      pcm::MEHPSimulator forceRelaxer =
        pcm::MEHPSimulator(universe2,
                           {},
                           2,
                           pcm::EntanglementMode::NONE,
                           nullptr,
                           false,
                           true,
                           false,
                           false);

      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      CHECK(forceBalancer2.getNrOfIterations() == 0);
      CHECK(forceBalancer2.getVolume() == Catch::Approx(universe2.getVolume()));
      CHECK(forceBalancer2.getVolume() ==
            Catch::Approx(97.383096 * 97.383096 * 97.383096));
      // initial system values
      CHECK_THAT(forceBalancer2.getPressure(),
                 Catch::Matchers::WithinRel(forceRelaxer.getPressure(), 3e-2));
      CHECK_THAT(forceBalancer2.getPressure(),
                 Catch::Matchers::WithinRel(0.0061105865, 3e-2));

      // run force relaxation
      CHECK_NOTHROW(forceBalancer2.run(pcm::MEHPLinearSolverConfig{}));
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
        (forceBalancer2.getNetwork().springContourLength.sum() / nrOfChains /
         Nb) *
        3. * kb * T / (slope * beadMass); // J/sigma^2
      CHECK_THAT(conversionFactor / (sigmaToM * sigmaToM * 79.),
                 Catch::Matchers::WithinRel(0.000245, 1e-3));
      double nu =
        nrOfChains / (forceBalancer2.getVolume() * sigmaToM * sigmaToM *
                      sigmaToM); // chain number density, m^-3
      CHECK(nu == Catch::Approx(4.63241e25));

      // final values
      CHECK(forceBalancer2.getPressure() ==
            Catch::Approx(0.153806 / 79.)); // LJ Units [?]
      CHECK_THAT(forceBalancer2.getPressure() * conversionFactor /
                   (sigmaToM * sigmaToM * sigmaToM),
                 Catch::Matchers::WithinRel(
                   61172.8878, 1e-3)); // shear modulus from the pressure, MPa
      double b02 = slope * beadMass;
      CHECK(forceBalancer2.getGammaFactor(1.) ==
            Catch::Approx(forceBalancer2.getGammaFactors(1.).mean()));
      CHECK(
        forceBalancer2.getGammaFactors(b02).mean() ==
        Catch::Approx((forceBalancer2.getGammaFactorsInDir(0, b02).mean() +
                       forceBalancer2.getGammaFactorsInDir(1, b02).mean() +
                       forceBalancer2.getGammaFactorsInDir(2, b02).mean())));
      CHECK(forceBalancer2.getGammaFactor(b02) ==
            Catch::Approx(forceBalancer2.getGammaFactors(b02).mean()));
      CHECK_THAT(
        forceBalancer2.getGammaFactors(b02).sum() / nrOfChains * kb * T * nu,
        Catch::Matchers::WithinRel(61308.3, 0.03)); // ANT shear modulus, Pa
      CHECK_THAT(
        forceBalancer2.getGammaFactor(b02, static_cast<int>(nrOfChains)),
        Catch::Matchers::WithinRel(
          forceBalancer2.getGammaFactor(b02, static_cast<int>(nrOfChains)),
          0.03));
      CHECK_THAT(
        forceBalancer2.getGammaFactor(b02, static_cast<int>(nrOfChains)),
        Catch::Matchers::WithinRel(
          0.319446, 0.03)); // "correct" gamma factor, see Mathematica script
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK(forceBalancer2.getNrOfActiveNodes() > 1);
      CHECK(forceBalancer2.getNrOfActiveStrands() > 1);
      CHECK(forceBalancer2.getNrOfActiveStrands() <=
            (countActiveStrandsInDir(forceBalancer2, 0) +
             countActiveStrandsInDir(forceBalancer2, 1) +
             countActiveStrandsInDir(forceBalancer2, 2)));
      CHECK(forceBalancer2.getNetwork()
              .evaluateStrandLengths(forceBalancer2.getCurrentDisplacements(),
                                     forceBalancer2.getIs2D())
              .mean() > 1.0);
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
    pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, true);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::UNSET);
    CHECK(forceBalancer.getNrOfIterations() == 0);
    CHECK(forceBalancer.getVolume() == Catch::Approx(universe.getVolume()));
    forceBalancer.run(pcm::MEHPLinearSolverConfig{});
    CHECK(forceBalancer.getNrOfNodes() < universe.getNrOfAtoms());
    CHECK(forceBalancer.getNrOfIterations() >= 1);
    CHECK(universe.getAtomsOfType(2).size() == 7200);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(forceBalancer.getGammaFactor(
            1., static_cast<int>(forceBalancer.getNrOfStrands())) ==
          Catch::Approx(1. / 3.).epsilon(0.001));
    auto stressTensor = forceBalancer.getStressTensor();
    CHECK(forceBalancer.getPressure() ==
          Catch::Approx(
            (stressTensor(0, 0) + stressTensor(1, 1) + stressTensor(2, 2)) / 3.)
            .epsilon(0.02));
    CHECK_NOTHROW(forceBalancer.validateNetwork());
    CHECK(forceBalancer.getNrOfActiveNodes() > 1);
    CHECK(forceBalancer.getNrOfActiveStrands() > 1);
    CHECK(forceBalancer.getNetwork()
            .evaluateStrandLengths(forceBalancer.getCurrentDisplacements(),
                                   forceBalancer.getIs2D())
            .mean() > 1.0);
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
  "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance2 can randomly add slip-links "
               "ignoring crosslinkers\""
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
  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe, 250, 2.0, 0.0, 100);
  size_t nrOfAddedLinks =
    forceBalancer.getNrOfLinks() - forceBalancer.getNrOfNodes();
  CHECK(nrOfAddedLinks >= 50);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
  std::cout << "Added " << nrOfAddedLinks << " slip-links" << std::endl;

  pcm::MEHPNetwork net = forceBalancer.getNetwork();
  Eigen::VectorXd displacements = forceBalancer.getCurrentDisplacements();
  CHECK(net.nrOfSprings > 0);
  CHECK_NOTHROW(net.validate(displacements));
  pe::Universe xlinkerVerse0 = forceBalancer.getCrosslinkerVerse();
  CHECK(xlinkerVerse0.getNrOfBonds() > 0);
  Eigen::VectorXd gammaFactorsX0 = forceBalancer.getGammaFactorsInDir(0, 1.);
  CHECK(gammaFactorsX0.size() ==
        static_cast<Eigen::Index>(forceBalancer.getNrOfSprings()));
  CHECK(gammaFactorsX0.minCoeff() >= 0);
  CHECK(gammaFactorsX0.maxCoeff() > 0);

  // remove all springs...
  size_t numRemoved = forceBalancer.removeInactiveLinks(1e5);
  net = forceBalancer.getNetwork();
  CHECK(forceBalancer.getNrOfActiveSprings() == 0);
  CHECK(net.nrOfSprings == 0);
  CHECK(numRemoved > 0);

  pe::Universe xlinkerVerse = forceBalancer.getCrosslinkerVerse();
  CHECK(xlinkerVerse.getNrOfBonds() == 0);

  std::unordered_map<long int, int> effF =
    forceBalancer.getEffectiveFunctionalityOfAtoms();
  for (const auto& pair : effF) {
    CHECK(pair.second == 0);
  }
  Eigen::VectorXd gammaFactorsX = forceBalancer.getGammaFactorsInDir(0, 1.);
  CHECK(gammaFactorsX.isZero());
  CHECK(countActiveSpringsConnected(forceBalancer).isZero());
}

TEST_CASE("MEHP Force Balance2 can randomly add and remove entanglements",
          "[analysis][MEHPSimulator]")
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
  pcm::MEHPSimulator forceBalancer =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe, 1000, 2.0, 0.0, 100);
  size_t nrOfAddedLinks =
    forceBalancer.getNrOfLinks() - forceBalancer.getNrOfNodes();
  CHECK(nrOfAddedLinks >= 100);

  size_t numRemoved = forceBalancer.unlinkBifunctionalLinks();
  CHECK_NOTHROW(forceBalancer.validateNetwork());
  // since we don't even consider f = 2 crosslinks,
  // none will be removed.
  CHECK(numRemoved == 0);

  // run a while to get inactive links
  forceBalancer.run(pcm::MEHPLinearSolverConfig{});
  size_t numInactiveRemoved = forceBalancer.removeInactiveLinks(0.1);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
  CHECK(numInactiveRemoved > 0);

  ////////////////////////////////////////////////////////////////
  forceBalancer = pcm::MEHPSimulator::constructWithRandomEntanglements(
    universe, 1000, 2.0, 0.0, 100);
  nrOfAddedLinks = forceBalancer.getNrOfLinks() - forceBalancer.getNrOfNodes();
  CHECK(nrOfAddedLinks >= 100);

  // run a while to get inactive links
  CHECK_NOTHROW(forceBalancer.run(pcm::MEHPLinearSolverConfig{}));
  numInactiveRemoved = forceBalancer.removeInactiveLinks(0.1);
  CHECK(numInactiveRemoved > 0);
  numInactiveRemoved = forceBalancer.unlinkBifunctionalLinks();
  // the `removeInactiveLinks` method already took care of these bifunctional
  // links
  CHECK(numInactiveRemoved == 0);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
}

TEST_CASE("MEHP Force Balance2 runs with non-network",
          "[analysis][MEHPSimulator][NonGaussianSpringForceEvaluator][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 runs with non-network\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  SECTION("MEHP Force Balance2 3D case")
  {
    std::string largeInputFile =
      suspectedPath +
      "structure/xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out";
    if (std::filesystem::exists(largeInputFile)) {
      universeSeq.initializeFromDataSequence({ { largeInputFile } });
      pe::Universe universe2 = universeSeq.atIndex(0);
      pcm::MEHPSimulator forceBalancer2 = pcm::MEHPSimulator(
        universe2, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::UNSET);
      CHECK_NOTHROW(forceBalancer2.run(pcm::MEHPLinearSolverConfig{}));
      CHECK(forceBalancer2.getNrOfIterations() > 1);
      CHECK(forceBalancer2.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK_THAT(forceBalancer2.getSolubleWeightFraction(),
                 Catch::Matchers::WithinAbs(1., 0.01));
    }
  }
}

TEST_CASE("MEHP Force Balance2 Free chains collapse",
          "[analysis][MEHPSimulator]")
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
  for (int i = 0; i < static_cast<int>(nrOfBeads); ++i) {
    const double offset = 10.0;
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i + 1);
    atomTypes.push_back(static_cast<size_t>(i) % nrOfBeadsPerChain == 0 ? 2
                                                                        : 1);
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
  pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
  CHECK(forceBalancer.getNrOfStrands() == forceBalancer.getNrOfSprings());
  CHECK(forceBalancer.getNrOfStrands() ==
        1); // we ignore crosslinks, so no "nrOfBeads / nrOfBeadsPerChain"

  CHECK_NOTHROW(forceBalancer.run(pcm::MEHPLinearSolverConfig{}));
  CHECK(forceBalancer.getNrOfIterations() > 0);
  CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  CHECK(forceBalancer.getNrOfActiveStrands() == 0);
  CHECK(forceBalancer.getActiveWeightFraction() == 0.);
  CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
             Catch::Matchers::WithinRel(1.));
  CHECK_THAT(forceBalancer.getDanglingWeightFraction(),
             Catch::Matchers::WithinAbs(0., 1e-9));
  double avgStrandLen =
    forceBalancer.getNetwork()
      .evaluateStrandLengths(forceBalancer.getCurrentDisplacements(), false)
      .mean();
  CHECK(avgStrandLen >= 0.0);
  CHECK(avgStrandLen <= 3e-6);
  CHECK_NOTHROW(forceBalancer.validateNetwork());
}

TEST_CASE("MEHP Force Balance2 Entanglement Beads Are Removed",
          "[analysis][MEHPSimulator]")
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
  for (int i = 0; i < static_cast<int>(nrOfBeads); ++i) {
    const double offset = 10.0;
    xPositions.push_back(i * 1.0 + offset);
    yPositions.push_back(0.1 * static_cast<double>(i % 4 - i % 3) +
                         offset); // /!\ i needs to be int, not unsigned!
    zPositions.push_back(0.1 * static_cast<double>(i % 5 - i % 7) + offset); //
    atomIds.push_back(i + 1);
    atomTypes.push_back(static_cast<size_t>(i) % nrOfBeadsPerChain == 0 ? 2
                                                                        : 1);
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

  pcm::MEHPLinearSolverConfig simpCfg;
  simpCfg.simplificationMode =
    pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;

  SECTION("Entanglement links")
  {
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        3,
        2.0,
        0.,
        0,
        0,
        "asdflseed",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        nullptr,
        false,
        false,
        false);

    CHECK((forceBalancer.getInitialCoordinates() +
           forceBalancer.getCurrentDisplacements())
            .isApprox(forceBalancer.getCoordinates()));

    CHECK(forceBalancer.getNrOfStrands() < forceBalancer.getNrOfSprings());
    CHECK(forceBalancer.getNrOfStrands() == 1); // since we ignore crosslinks
    CHECK(forceBalancer.getNrOfSprings() == 1 + 3 * 2);
    CHECK_NOTHROW(forceBalancer.run(simpCfg));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getNrOfActiveStrands() == 0);
    CHECK(forceBalancer.getNrOfStrands() == 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork());
  }
  SECTION("Entanglement springs")
  {
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        3,
        2.0,
        0.,
        0,
        0,
        "asdflseed",
        2,
        pcm::EntanglementMode::SPRING,
        true,
        false,
        nullptr,
        false,
        false,
        false);
    CHECK(forceBalancer.getNrOfStrands() < forceBalancer.getNrOfSprings());
    CHECK(forceBalancer.getNrOfStrands() ==
          1 + 3); // 1 for the strand, 3 per entanglement spring

    CHECK_NOTHROW(forceBalancer.run(simpCfg));
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getNrOfActiveStrands() == 0);
    CHECK(forceBalancer.getNrOfStrands() == 0);
    CHECK_NOTHROW(forceBalancer.validateNetwork());

    CHECK((forceBalancer.getInitialCoordinates() +
           forceBalancer.getCurrentDisplacements())
            .isApprox(forceBalancer.getCoordinates()));
  }
}

TEST_CASE("MEHP Force Balance2 fully active chains are fully active",
          "[analysis][MEHPSimulator][long]")
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

    pcm::MEHPLinearSolverConfig simpCfg;
    simpCfg.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;

    SECTION("Without entanglements")
    {
      pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(
        universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
      size_t initialNSprings = forceBalancer.getNrOfStrands();

      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      double initialResidual = forceBalancer.getResidualNorm();
      CHECK(std::isfinite(initialResidual));
      REQUIRE_NOTHROW(forceBalancer.run(simpCfg));
      CHECK(forceBalancer.getNrOfIterations() > 0);
      CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      CHECK(forceBalancer.getNrOfStrands() == initialNSprings);
      CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                 Catch::Matchers::WithinRel(1.0));
      CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                 Catch::Matchers::WithinAbs(0.0, 1e-9));
      CHECK(initialResidual > forceBalancer.getResidualNorm());
    }

    SECTION("With entanglement links")
    {
      pcm::MEHPSimulator forceBalancer =
        pcm::MEHPSimulator::constructWithRandomEntanglements(
          universe,
          400,
          2.0,
          0.0,
          100,
          0.0,
          "a533d",
          2,
          pcm::EntanglementMode::MERGED_LINKS,
          true,
          false,
          nullptr,
          false,
          false,
          false);

      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      double initialResidual = forceBalancer.getResidualNorm();
      CHECK(std::isfinite(initialResidual));
      REQUIRE_NOTHROW(forceBalancer.run(simpCfg));
      CHECK(forceBalancer.getNrOfIterations() > 0);
      CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                 Catch::Matchers::WithinRel(1.0, 1e-6));
      CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                 Catch::Matchers::WithinAbs(0.0, 1e-6));
      CHECK(initialResidual > forceBalancer.getResidualNorm());
    }

    SECTION("With entanglement springs")
    {
      pcm::MEHPSimulator forceBalancer =
        pcm::MEHPSimulator::constructWithRandomEntanglements(
          universe,
          400,
          2.0,
          0.0,
          100,
          0.0,
          "a533d",
          2,
          pcm::EntanglementMode::SPRING,
          true,
          false,
          nullptr,
          false,
          false,
          false);

      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      double initialResidual = forceBalancer.getResidualNorm();
      CHECK(std::isfinite(initialResidual));
      REQUIRE_NOTHROW(forceBalancer.run(simpCfg));
      CHECK(forceBalancer.getNrOfIterations() > 0);
      CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
      CHECK(forceBalancer.getNrOfActiveStrands() ==
            forceBalancer.getNrOfStrands());
      CHECK_THAT(forceBalancer.getActiveWeightFraction(),
                 Catch::Matchers::WithinRel(1.0, 1e-6));
      CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                 Catch::Matchers::WithinAbs(0.0, 1e-9));
      CHECK(initialResidual > forceBalancer.getResidualNorm());
    }
  }
}

TEST_CASE("MEHPForceBalance2 gives approx. same results for entanglement links "
          "& springs in diamond network",
          "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHPForceBalance2 gives approx. same results "
               "for entanglement links & springs in diamond network\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

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

  pcm::MEHPSimulator forceBalanceLinks = pcm::MEHPSimulator(
    universe, entanglements, 2, pcm::EntanglementMode::MERGED_LINKS);

  pcm::MEHPSimulator forceBalanceSprings = pcm::MEHPSimulator(
    universe, entanglements, 2, pcm::EntanglementMode::SPRING);

  CHECK_THAT(forceBalanceLinks.getStressTensor().trace(),
             Catch::Matchers::WithinRel(
               forceBalanceSprings.getStressTensor().trace(), 1e-1));
  CHECK_THAT(
    forceBalanceLinks.getResidualNorm(),
    Catch::Matchers::WithinRel(forceBalanceSprings.getResidualNorm(), 1e-1));

  forceBalanceLinks.run(pcm::MEHPLinearSolverConfig{});
  forceBalanceSprings.run(pcm::MEHPLinearSolverConfig{});

  CHECK_THAT(forceBalanceLinks.getStressTensor().trace(),
             Catch::Matchers::WithinRel(
               forceBalanceSprings.getStressTensor().trace(), 1e-1));
  CHECK_THAT(
    forceBalanceLinks.getResidualNorm(),
    Catch::Matchers::WithinAbs(forceBalanceSprings.getResidualNorm(), 1e-6));
  CHECK_THAT(forceBalanceLinks.getActiveWeightFraction(),
             Catch::Matchers::WithinAbs(
               forceBalanceSprings.getActiveWeightFraction(), 1e-3));
  CHECK(forceBalanceLinks.getNrOfSprings() <
        forceBalanceSprings.getNrOfSprings());

  pcm::MEHPLinearSolverConfig simpCfg;
  simpCfg.simplificationMode =
    pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
  simpCfg.inactiveRemovalCutoff = 0.01;
  forceBalanceLinks.run(simpCfg);
  forceBalanceSprings.run(simpCfg);

  CHECK_THAT(forceBalanceLinks.getStressTensor().trace(),
             Catch::Matchers::WithinRel(
               forceBalanceSprings.getStressTensor().trace(), 0.1));
  CHECK_THAT(
    forceBalanceLinks.getResidualNorm(),
    Catch::Matchers::WithinAbs(forceBalanceSprings.getResidualNorm(), 1e-3));
  CHECK(forceBalanceLinks.getNrOfStrands() <
        forceBalanceSprings.getNrOfStrands());
  CHECK_THAT(forceBalanceLinks.getActiveWeightFraction(),
             Catch::Matchers::WithinAbs(
               forceBalanceSprings.getActiveWeightFraction(), 1e-3));
  CHECK(forceBalanceLinks.getNrOfSprings() <
        forceBalanceSprings.getNrOfSprings());

  forceBalanceLinks.run(simpCfg);
  forceBalanceSprings.run(simpCfg);

  CHECK(forceBalanceLinks.getNrOfSprings() <
        forceBalanceSprings.getNrOfSprings());
  CHECK_THAT(forceBalanceLinks.getStressTensor().trace(),
             Catch::Matchers::WithinRel(
               forceBalanceSprings.getStressTensor().trace(), 1e-1));
  CHECK_THAT(
    forceBalanceLinks.getResidualNorm(),
    Catch::Matchers::WithinAbs(forceBalanceSprings.getResidualNorm(), 1e-6));
  CHECK_THAT(forceBalanceLinks.getActiveWeightFraction(),
             Catch::Matchers::WithinAbs(
               forceBalanceSprings.getActiveWeightFraction(), 1e-3));
}

TEST_CASE("MEHPForceBalance2 gives approx. same results for entanglement links "
          "& springs in p = 1 network",
          "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHPForceBalance2 gives approx. same results "
               "for entanglement links "
               "& springs in p = 1 network\""
            << std::endl;

  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

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

  pcm::MEHPSimulator forceBalanceLinks = pcm::MEHPSimulator(
    universe, entanglements, 2, pcm::EntanglementMode::MERGED_LINKS);

  pcm::MEHPSimulator forceBalanceSprings = pcm::MEHPSimulator(
    universe, entanglements, 2, pcm::EntanglementMode::SPRING);

  CHECK_THAT(forceBalanceLinks.getStressTensor().trace(),
             Catch::Matchers::WithinRel(
               forceBalanceSprings.getStressTensor().trace(), 1e-3));
  CHECK_THAT(
    forceBalanceLinks.getResidualNorm(),
    Catch::Matchers::WithinAbs(forceBalanceSprings.getResidualNorm(), 1e-2));

  forceBalanceLinks.run(pcm::MEHPLinearSolverConfig{});
  forceBalanceSprings.run(pcm::MEHPLinearSolverConfig{});

  CHECK_THAT(forceBalanceLinks.getStressTensor().trace(),
             Catch::Matchers::WithinRel(
               forceBalanceSprings.getStressTensor().trace(), 3e-3));
  CHECK_THAT(
    forceBalanceLinks.getResidualNorm(),
    Catch::Matchers::WithinAbs(forceBalanceSprings.getResidualNorm(), 1e-2));
  CHECK_THAT(forceBalanceLinks.getActiveWeightFraction(),
             Catch::Matchers::WithinAbs(
               forceBalanceSprings.getActiveWeightFraction(), 1e-3));

  pcm::MEHPLinearSolverConfig simpCfg;
  simpCfg.simplificationMode =
    pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
  simpCfg.inactiveRemovalCutoff = 0.01;
  forceBalanceLinks.run(simpCfg);
  forceBalanceSprings.run(simpCfg);

  CHECK_THAT(forceBalanceLinks.getStressTensor().trace(),
             Catch::Matchers::WithinRel(
               forceBalanceSprings.getStressTensor().trace(), 3e-3));
  CHECK_THAT(
    forceBalanceLinks.getResidualNorm(),
    Catch::Matchers::WithinAbs(forceBalanceSprings.getResidualNorm(), 1e-3));
  CHECK_THAT(forceBalanceLinks.getActiveWeightFraction(),
             Catch::Matchers::WithinAbs(
               forceBalanceSprings.getActiveWeightFraction(), 1e-3));
}

TEST_CASE("MEHP Force Balance2 does not collapse", "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance2 does not collapse\""
            << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  universe.addAtoms(
    { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } },
    { { 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1 } },
    { { 0.0, 2.5, 5.0, 7.5, 0.1, 2.5, 5.0, 7.5, -.1, 5.0, 0.0, 5.0 } },
    { { 0.1, 0.0, -.1, 0.0, 5.0, 5.0, 5.0, 5.0, 2.5, 2.5, 7.5, 7.5 } },
    { { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } });
  universe.addBonds({ { 1,
                        1,
                        1, //
                        1,
                        3,
                        3, //
                        3,
                        3,
                        5, //
                        5,
                        5,
                        5, //
                        7,
                        7,
                        7, //
                        7 } },
                    { { 2,
                        9,
                        4, //
                        11,
                        2,
                        4, //
                        10,
                        12,
                        9, //
                        11,
                        6,
                        8, //
                        6,
                        8,
                        10, //
                        12 } });

  pcm::MEHPSimulator forceBalanceNew = pcm::MEHPSimulator(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, true);
  pcm::MEHPLinearSolverConfig cfg;
  cfg.simplificationMode = pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
  CHECK_NOTHROW(forceBalanceNew.run(cfg));
  CHECK(forceBalanceNew.getNrOfIterations() > 0);
  CHECK(forceBalanceNew.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  CHECK(forceBalanceNew.getNrOfActiveStrands() ==
        forceBalanceNew.getNrOfStrands());
  CHECK(forceBalanceNew.getNrOfActiveStrands() == 8);
  CHECK((forceBalanceNew.getNetwork().springContourLength.array() == 2).all());
  CHECK((countActiveSpringsConnected(forceBalanceNew).array() == 4).all());
  CHECK((countActiveStrandsConnected(forceBalanceNew).array() == 4).all());

  CHECK(forceBalanceNew.getNrOfActiveStrands() <=
        (countActiveStrandsInDir(forceBalanceNew, 0) +
         countActiveStrandsInDir(forceBalanceNew, 1) +
         countActiveStrandsInDir(forceBalanceNew, 2)));
  CHECK(forceBalanceNew.getNrOfActiveNodes() == 4);
  CHECK(forceBalanceNew.getNetwork()
          .evaluateStrandLengths(forceBalanceNew.getCurrentDisplacements(),
                                 forceBalanceNew.getIs2D())
          .mean() == Catch::Approx(5.0));
  CHECK_THAT(forceBalanceNew.getGammaFactor(1., 100),
             Catch::Matchers::WithinAbs(1.0, 1e-2));
  CHECK_THAT(forceBalanceNew.getPressure(),
             Catch::Matchers::WithinAbs(1. / 30., 1e-3));
  CHECK_THAT(forceBalanceNew.getSolubleWeightFraction(),
             Catch::Matchers::WithinAbs(0.0, 1e-9));
};

TEST_CASE("MEHP Force Balance2 correctly collapses melts even with random "
          "entanglements",
          "[analysis][MEHPSimulator]")
{
  std::cout << "Running test \"MEHP Force Balance2 correctly collapses melts "
               "even with random entanglements\""
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

  pcm::MEHPLinearSolverConfig simpCfg;
  simpCfg.simplificationMode =
    pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;

  SECTION("Phantom")
  {
    pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(universe);
    pcm::MEHPNetwork net0 = forceBalancer.getNetwork();
    CHECK(net0.nrOfStrands == 83);
    CHECK(net0.nrOfNodes == 2 * 83);
    for (size_t i = 0; i < net0.nrOfLinks; ++i) {
      CHECK(!net0.linkIsEntanglement[static_cast<Eigen::Index>(i)]);
      CHECK(net0.strandIndicesOfLink[i].size() == 1);
    }

    CHECK_NOTHROW(forceBalancer.validateNetwork());
    double initialResidual = forceBalancer.getResidualNorm();
    CHECK(std::isfinite(initialResidual));
    CHECK(forceBalancer.getNumExtraAtoms() == 0);
    CHECK(forceBalancer.getNetwork().nrOfSprings ==
          forceBalancer.getNetwork().nrOfStrands);
    CHECK_NOTHROW(forceBalancer.run(simpCfg));
    pcm::MEHPNetwork net = forceBalancer.getNetwork();
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveStrands() == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 0) == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 1) == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 2) == 0);
    CHECK(net.nrOfSprings == 0);
    CHECK(net.nrOfStrands == 0);
    CHECK(initialResidual > forceBalancer.getResidualNorm());
  }

  SECTION("Entangled with links")
  {
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        1000,
        2.0,
        0.0,
        100,
        5,
        "my_seed_fb12",
        2,
        pcm::EntanglementMode::MERGED_LINKS);
    pcm::MEHPNetwork net0 = forceBalancer.getNetwork();
    CHECK(net0.nrOfStrands == 83);
    CHECK(net0.nrOfNodes == 2 * 83);
    for (size_t i = 0; i < net0.nrOfLinks; ++i) {
      if (!net0.linkIsEntanglement[static_cast<Eigen::Index>(i)]) {
        CHECK(net0.strandIndicesOfLink[i].size() == 1);
        bool any =
          net0.linkIndicesOfStrand[net0.strandIndicesOfLink[i][0]][0] == i ||
          net0.linkIndicesOfStrand[net0.strandIndicesOfLink[i][0]].back() == i;
        CHECK(any);
      }
    }

    CHECK_NOTHROW(forceBalancer.validateNetwork());
    double initialResidual = forceBalancer.getResidualNorm();
    CHECK(std::isfinite(initialResidual));
    CHECK(forceBalancer.getNumExtraAtoms() > 100);
    CHECK(forceBalancer.getNetwork().nrOfSprings >
          forceBalancer.getNetwork().nrOfStrands);
    CHECK_NOTHROW(forceBalancer.run(simpCfg));
    pcm::MEHPNetwork net = forceBalancer.getNetwork();
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveStrands() == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 0) == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 1) == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 2) == 0);
    CHECK(net.nrOfSprings == 0);
    CHECK(net.nrOfStrands == 0);
    CHECK(initialResidual > forceBalancer.getResidualNorm());
  }

  SECTION("Entangled with springs")
  {
    pcm::MEHPSimulator forceBalancer =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        1000,
        2.0,
        0.0,
        100,
        5,
        "my_seed_fb12",
        2,
        pcm::EntanglementMode::SPRING);
    pcm::MEHPNetwork net0 = forceBalancer.getNetwork();
    CHECK(net0.nrOfStrands ==
          83 + static_cast<size_t>(net0.springIsEntanglement.count()));
    CHECK(net0.nrOfNodes == 2 * 83);

    CHECK_NOTHROW(forceBalancer.validateNetwork());
    double initialResidual = forceBalancer.getResidualNorm();
    CHECK(std::isfinite(initialResidual));
    CHECK(forceBalancer.getNumExtraAtoms() > 100);
    CHECK(forceBalancer.getNetwork().nrOfSprings >
          forceBalancer.getNetwork().nrOfStrands);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = pcm::SLESolver::DEFAULT;
    cfg.residualReduction = 1e-10;
    cfg.maxIterations = 50000;
    CHECK_NOTHROW(forceBalancer.run(cfg));
    pcm::MEHPNetwork net = forceBalancer.getNetwork();
    CHECK(forceBalancer.getNrOfIterations() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getNrOfActiveStrands() == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 0) == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 1) == 0);
    CHECK(countActiveStrandsInDir(forceBalancer, 2) == 0);
    CHECK(net.nrOfSprings == 0);
    CHECK(net.nrOfStrands == 0);
    CHECK(initialResidual > forceBalancer.getResidualNorm());
  }
}

TEST_CASE("Particular MEHP Force Balance2 Example", "[analysis][MEHPSimulator]")
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

  pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(universe);
  double pressBefore = forceBalancer.getPressure();
  CHECK_NOTHROW(forceBalancer.run(pcm::MEHPLinearSolverConfig{}));
  CHECK(pressBefore > forceBalancer.getPressure());
}

TEST_CASE(
  "MEHPForceBalance2 Adding slip-links does not influence other springs",
  "[analysis][MEHPSimulator][long]")
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

  pcm::MEHPSimulator forceBalancerWithoutEntanglements =
    pcm::MEHPSimulator(universe);
  pcm::MEHPSimulator forceBalancerOldSamplingSmall =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe, 1000, 6.0, 0.0, 900);
  Eigen::VectorXd springVectorsWithoutEntanglements =
    forceBalancerWithoutEntanglements.getNetwork().evaluateStrandVectors(
      forceBalancerWithoutEntanglements.getCurrentDisplacements());

  Eigen::VectorXd springVectorsWithEntanglements =
    forceBalancerOldSamplingSmall.getNetwork().evaluateStrandVectors(
      forceBalancerOldSamplingSmall.getCurrentDisplacements());

  CHECK(
    springVectorsWithoutEntanglements.isApprox(springVectorsWithEntanglements));

  pcm::MEHPSimulator forceBalancer2 =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe, 1000, 6.0, 0.0, 900, 3.0, "53467829");

  pcm::MEHPSimulator forceBalancer2Without =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe, 0, 6.0, 0.0, 0, 3.0, "53467829");

  CHECK(forceBalancer2.getNrOfSprings() >
        forceBalancer2Without.getNrOfSprings());

  CHECK(forceBalancer2.getNetwork()
          .evaluateStrandVectors(forceBalancer2.getCurrentDisplacements())
          .isApprox(forceBalancer2Without.getNetwork().evaluateStrandVectors(
            forceBalancer2Without.getCurrentDisplacements())));
}

TEST_CASE("MEHPForceBalance2 Bugfix Validations 1A.1 1A.2 1A.4",
          "[MEHPSimulator][bugfixes]")
{
  pe::UniverseSequence universeSeq;
  std::string inputFile = std::string(PYLIMER_TEST_FIXTURES_DIR) +
                          "/structure/network_100_a_46.structure.out";
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);

  // System with slip-links
  pcm::MEHPSimulator fb = pcm::MEHPSimulator::constructWithRandomEntanglements(
    universe, 50, 4.0, 0.0, 50, 3.0, "12345");
  REQUIRE(fb.getNetwork().nrOfLinks > fb.getNetwork().nrOfNodes);

  // 1A.1: evaluateStrandVectors vs individual evaluateStrandVector
  Eigen::VectorXd allStrandVectors =
    fb.getNetwork().evaluateStrandVectors(fb.getCurrentDisplacements());
  REQUIRE(allStrandVectors.size() ==
          static_cast<Eigen::Index>(3 * fb.getNrOfStrands()));
  for (size_t s = 0; s < fb.getNrOfStrands(); ++s) {
    Eigen::Vector3d single =
      fb.getNetwork().evaluateStrandVector(fb.getCurrentDisplacements(), s);
    Eigen::Vector3d fromAll =
      allStrandVectors.segment(3 * static_cast<Eigen::Index>(s), 3);
    CHECK(single.isApprox(fromAll));
  }

  // 1A.2: findSpringsAndLinksClusteredToActive size
  auto clustered = fb.getNetwork().findSpringsAndLinksClusteredToActive(
    fb.getCurrentDisplacements());
  CHECK(static_cast<size_t>(clustered.second.size()) ==
        fb.getNetwork().nrOfLinks);

  // 1A.4: residual and stress scale consistently with kappa != 1.0
  const double resKappa1 = fb.getResidualNorm();
  const Eigen::Matrix3d stressKappa1 = fb.getStressTensor();
  fb.setForceEvaluator(std::make_shared<pcm::HookeanForceEvaluator>(2.5));
  CHECK_THAT(fb.getResidualNorm(),
             Catch::Matchers::WithinRel(2.5 * resKappa1, 1e-6));
  CHECK(fb.getStressTensor().isApprox(2.5 * stressKappa1, 1e-6));
}

TEST_CASE(
  "MEHP Force Balance2 Entanglement Links are just like Entanglement Springs",
  "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 Entanglement Links are just "
               "like Entanglement Springs\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

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
  pcm::MEHPSimulator forceBalancerEntanglementSprings = pcm::MEHPSimulator(
    universe, entanglements, 2, pcm::EntanglementMode::SPRING);

  pcm::MEHPSimulator forceBalancerEntanglementLinks = pcm::MEHPSimulator(
    universe, entanglements, 2, pcm::EntanglementMode::MERGED_LINKS);

  REQUIRE(forceBalancerEntanglementSprings.getNrOfStrands() >
          forceBalancerEntanglementLinks.getNrOfStrands());
  REQUIRE(forceBalancerEntanglementLinks.getNrOfStrands() == 464);
  REQUIRE(forceBalancerEntanglementSprings.getNrOfActiveSprings() >
          forceBalancerEntanglementLinks.getNrOfActiveSprings());

  std::vector<long int> activeNodes0_1 =
    getAtomIdsOfActiveNodes(forceBalancerEntanglementSprings);
  std::ranges::sort(activeNodes0_1);
  std::vector<long int> activeNodes0_2 =
    getAtomIdsOfActiveNodes(forceBalancerEntanglementLinks);
  std::ranges::sort(activeNodes0_2);
  REQUIRE(activeNodes0_1 == activeNodes0_2);

  SECTION("With removal of inactive springs")
  {
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
    forceBalancerEntanglementSprings.run(cfg);
    forceBalancerEntanglementLinks.run(cfg);

    CHECK(forceBalancerEntanglementSprings.getNrOfStrands() >
          forceBalancerEntanglementLinks.getNrOfStrands());
    std::vector<long int> activeNodes1 =
      getAtomIdsOfActiveNodes(forceBalancerEntanglementSprings);
    std::ranges::sort(activeNodes1);
    std::vector<long int> activeNodes2 =
      getAtomIdsOfActiveNodes(forceBalancerEntanglementLinks);
    std::ranges::sort(activeNodes2);
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
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    forceBalancerEntanglementSprings.run(cfg);
    forceBalancerEntanglementLinks.run(cfg);

    CHECK(forceBalancerEntanglementSprings.getNrOfStrands() >
          forceBalancerEntanglementLinks.getNrOfStrands());
    std::vector<long int> activeNodes1 =
      getAtomIdsOfActiveNodes(forceBalancerEntanglementSprings);
    std::ranges::sort(activeNodes1);
    std::vector<long int> activeNodes2 =
      getAtomIdsOfActiveNodes(forceBalancerEntanglementLinks);
    std::ranges::sort(activeNodes2);
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
  "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"MEHP Force Balance2 adding Entanglements vs. "
               "Phantom behaves as expected\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath +
    "/structure/mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0."
    "0114_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  auto masses = universe.getMasses();
  std::cout << "Read file " << inputFile << std::endl;

  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
        universe, 4368, 4., 0., 190, 0, "af1346lhkdsaöf123", 2, true);

  pcm::MEHPSimulator forceBalancerEntanglements =
    pcm::MEHPSimulator(universe,
                       entanglements,
                       2,
                       pcm::EntanglementMode::MERGED_LINKS,
                       nullptr,
                       false);

  pcm::MEHPSimulator forceBalancerPhantom = pcm::MEHPSimulator(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);

  CHECK(forceBalancerEntanglements.getNrOfSprings() >
        forceBalancerPhantom.getNrOfSprings());

  forceBalancerEntanglements.run(pcm::MEHPLinearSolverConfig{});
  forceBalancerPhantom.run(pcm::MEHPLinearSolverConfig{});

  CHECK(forceBalancerEntanglements.getNrOfSprings() >
        forceBalancerPhantom.getNrOfSprings());
  CHECK(forceBalancerEntanglements.getNrOfActiveStrands() >
        forceBalancerPhantom.getNrOfActiveStrands());
  CHECK(forceBalancerEntanglements.getSolubleWeightFraction() <=
        forceBalancerPhantom.getSolubleWeightFraction());
  CHECK(forceBalancerEntanglements.getActiveWeightFraction() >=
        forceBalancerPhantom.getActiveWeightFraction());

  pcm::MEHPLinearSolverConfig simpCfg;
  simpCfg.simplificationMode =
    pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
  forceBalancerEntanglements.run(simpCfg);
  forceBalancerPhantom.run(simpCfg);

  CHECK(forceBalancerEntanglements.getNrOfStrands() >
        forceBalancerPhantom.getNrOfStrands());
  CHECK(forceBalancerEntanglements.getNrOfActiveStrands() >
        forceBalancerPhantom.getNrOfActiveStrands());
  CHECK(forceBalancerEntanglements.getSolubleWeightFraction() <
        forceBalancerPhantom.getSolubleWeightFraction());
};

TEST_CASE("MEHPForceBalance2 phantom with and without removal is the same",
          "[analysis][MEHPSimulator]")
{
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath +
    "/structure/mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0."
    "0114_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  auto masses = universe.getMasses();
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPSimulator forceBalancerPhantom = pcm::MEHPSimulator(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
  pcm::MEHPSimulator forceBalancerPhantomRem = pcm::MEHPSimulator(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);

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
  pcm::MEHPLinearSolverConfig cfgNoSimp;
  cfgNoSimp.simplificationMode =
    pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
  cfgNoSimp.inactiveRemovalCutoff = 1e-9;
  cfgNoSimp.solver = pcm::SLESolver::CONJUGATE_GRADIENT;
  forceBalancerPhantom.run(cfgNoSimp);

  pcm::MEHPLinearSolverConfig cfgSimp = cfgNoSimp;
  cfgSimp.simplificationMode =
    pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
  forceBalancerPhantomRem.run(cfgSimp);

  // compare results
  CHECK(forceBalancerPhantom.getNrOfSprings() >
        forceBalancerPhantomRem.getNrOfSprings());
  CHECK_THAT(
    forceBalancerPhantom.getResidual(),
    Catch::Matchers::WithinAbs(forceBalancerPhantomRem.getResidual(), 1e-6));
  CHECK_THAT(forceBalancerPhantom.getSolubleWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getSolubleWeightFraction(), 1e-3));
  CHECK_THAT(forceBalancerPhantom.getActiveWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getActiveWeightFraction(), 1e-3));
  CHECK_THAT(forceBalancerPhantom.getDanglingWeightFraction(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getDanglingWeightFraction(), 1e-3));
  CHECK_THAT(forceBalancerPhantom.getGammaFactors(1.).sum(),
             Catch::Matchers::WithinRel(
               forceBalancerPhantomRem.getGammaFactors(1.).sum(), 1e-4));
}

TEST_CASE("Temporary force balance 2 test case",
          "[analysis][MEHPSimulator][long]")
{
  std::cout << "Running test \"Temporary force balance test case\""
            << std::endl;
  pe::UniverseSequence universeSeq = pe::UniverseSequence();
  CHECK(universeSeq.getLength() == 0);
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;

  std::string inputFile =
    suspectedPath +
    "/tmp/mc_own-si_pdms_crosslinked_melt_2590_a_221_2410_monoa_271_"
    "r_1.44_wsol_0.278_f_4_v_1.structure.out";

  std::cout << "Reading file " << inputFile << std::endl;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);
  auto masses = universe.getMasses();
  std::cout << "Read file " << inputFile << std::endl;

  pcm::MEHPSimulator forceBalancerEntanglements =
    pcm::MEHPSimulator::constructWithRandomEntanglements(
      universe,
      26832,
      3.25,
      0.,
      26000,
      0,
      "",
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      true,
      false,
      nullptr,
      false);

  std::cout << "Sampled entanglements " << std::endl;

  pcm::MEHPLinearSolverConfig cfg;
  cfg.simplificationMode = pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
  CHECK_NOTHROW(forceBalancerEntanglements.run(cfg));
};

TEST_CASE("MEHPFB2 Basic conversion test",
          "[analysis][MEHPSimulator][StructureConversion]")
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
  for (size_t i = 0; i < coordsX.size(); ++i) {
    ids.push_back(static_cast<long int>(i));
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
  std::vector<long int> bondFrom = { 0, 1,  2, 3,  4,  5,  0,  6,  7, 8,
                                     9, 10, 0, 11, 12, 13, 14, 15, 2 };
  std::vector<long int> bondTo = { 1,  2, 3,  4,  5,  0,  6,  7, 8, 9,
                                   10, 0, 11, 12, 13, 14, 15, 0, 2 };
  universe.addBonds(bondFrom, bondTo);

  CHECK(universe.getNrOfBonds() == bondFrom.size());
  CHECK(universe.getNrOfAtoms() == ids.size());

  pcm::MEHPSimulator fb1 = pcm::MEHPSimulator(universe);
  pcm::MEHPSimulator fb2 = pcm::MEHPSimulator(universe);

  pcm::MEHPNetwork net = fb2.getNetwork();

  CHECK(net.nrOfLinks == 2);
  CHECK(net.nrOfSprings == 5);
  CHECK(net.nrOfStrands == 5);
  CHECK(net.nrOfNodes == 2);
  CHECK_FALSE(net.springBoxOffset.isZero());
  CHECK_THROWS(net.getNeighbourLinkIndices(100));

  CHECK(fb1.getNrOfSprings() == 5);

  CHECK_THAT(fb1.getGammaFactors(1.0).mean(),
             Catch::Matchers::WithinAbs(fb2.getGammaFactors(1.0).mean(), 1e-6));
  CHECK_THAT(fb1.getResidual(),
             Catch::Matchers::WithinRel(fb2.getResidual(), 1e-6));

  pcm::MEHPSlipLinkSolverConfig slCfg;
  slCfg.maxNrOfSteps = 25000;
  slCfg.xtol = 1e-13;
  slCfg.simplificationMode =
    pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
  slCfg.inactiveRemovalCutoff = 1e-12;
  slCfg.doInnerIterations = true;
  slCfg.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
  slCfg.oneOverSpringPartitionUpperLimit = -1.0;
  slCfg.slipMode = pcm::SlipMode::FIXED;
  fb1.run(slCfg);
  fb2.run(pcm::MEHPLinearSolverConfig{});

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

TEST_CASE("MEHPFB2 Conversion test with PBC-breaking entanglements",
          "[analysis][MEHPSimulator][StructureConversion]")
{
  std::cout << "Running test \"MEHPFB2 Conversion test with PBC-breaking "
               "entanglements\""
            << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  std::vector<double> coords;
  std::vector<long int> ids;
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  for (size_t i = 0; i < 40; ++i) {
    coords.push_back(i * 1.);
    ids.push_back(static_cast<long int>(i));
    bondFrom.push_back(static_cast<long int>(i));
    bondTo.push_back(static_cast<long int>((i + 1) % 40));
  }
  std::vector<int> types =
    pylimer_tools::utils::initializeWithValue(ids.size(), 1);
  types[2] = 0;
  types[2] = 2;
  universe.addAtoms(ids,
                    types,
                    coords,
                    coords,
                    coords,
                    pylimer_tools::utils::initializeWithValue(ids.size(), 0),
                    pylimer_tools::utils::initializeWithValue(ids.size(), 0),
                    pylimer_tools::utils::initializeWithValue(ids.size(), 0));
  bondFrom.push_back(2);
  bondTo.push_back(0);
  bondFrom.push_back(0);
  bondTo.push_back(2);
  universe.addBonds(bondFrom, bondTo);

  CHECK(universe.getNrOfBonds() == bondFrom.size());
  CHECK(universe.getNrOfAtoms() == ids.size());

  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements;
  entanglements.pairOfAtom =
    pylimer_tools::utils::initializeWithValue<long int>(ids.size(), -1);

  SECTION("Phantom")
  {
    pcm::MEHPSimulator fb1_phantom =
      pcm::MEHPSimulator(universe, entanglements);
    pcm::MEHPSimulator fb2_phantom =
      pcm::MEHPSimulator(universe, entanglements);
    CHECK_THAT(fb1_phantom.getGammaFactors(1.0).mean(),
               Catch::Matchers::WithinAbs(
                 fb2_phantom.getGammaFactors(1.0).mean(), 1e-6));
    CHECK_THAT(fb1_phantom.getResidual(),
               Catch::Matchers::WithinRel(fb2_phantom.getResidual(), 1e-6));

    Eigen::Vector3d relevantStrandVector =
      fb2_phantom.getNetwork().evaluateStrandVector(
        fb2_phantom.getCurrentDisplacements(), 3);
    CHECK(relevantStrandVector.isApprox(
      Eigen::Vector3d(40. - 2., 40. - 2., 40. - 2.)));
  }

  SECTION("Entanglements through the box")
  {
    entanglements.pairsOfAtoms.push_back(std::make_pair(3, 24));
    entanglements.pairsOfAtoms.push_back(std::make_pair(5, 14));
    entanglements.pairOfAtom[3] = 0;
    entanglements.pairOfAtom[24] = 0;
    entanglements.pairOfAtom[5] = 1;
    entanglements.pairOfAtom[14] = 1;

    pcm::MEHPSimulator fb1 = pcm::MEHPSimulator(universe, entanglements);
    pcm::MEHPSimulator fb2 = pcm::MEHPSimulator(universe, entanglements);

    pcm::MEHPNetwork net2 = fb2.getNetwork();
    pcm::MEHPNetwork net1 = fb1.getNetwork();

    CHECK(fb2.getCurrentDisplacements().isZero());
    CHECK(fb1.getCurrentDisplacements().isZero());

    CHECK_FALSE(net1.springBoxOffset.isZero());
    CHECK_FALSE(net2.springBoxOffset.isZero());
    CHECK(net1.nrOfLinks == 4); // 2x crosslink (0, 2), 2x entanglement link
    CHECK(net1.nrOfNodes == 2);
    CHECK(net2.nrOfLinks == 4); // 2x junction-link (0, 2), 2x entanglement link
    CHECK(net2.nrOfNodes == 2);
    CHECK(net2.oldAtomIds[0] == 0);
    CHECK(net2.oldAtomIds[1] == 2);
    CHECK(net2.oldAtomIds[2] == -1);
    CHECK(net2.oldAtomIds[3] == -1);
    CHECK(net2.linkIsEntanglement[2]);
    CHECK(net2.linkIsEntanglement[3]);

    CHECK(net1.coordinates.segment(0, 3).isZero());
    CHECK(net2.coordinates.segment(0, 3).isZero());
    CHECK(net1.coordinates.segment(3 * 1, 3).isApprox(
      net2.coordinates.segment(3 * 1, 3)));
    REQUIRE(net1.nrOfLinks == net2.nrOfLinks);
    REQUIRE(net1.nrOfSprings == net2.nrOfSprings);
    CHECK(net2.nrOfSprings == 8);
    CHECK(net2.nrOfStrands == 4);

    CHECK_THAT(
      fb1.getResiduals().segment<3>(3 * 1).norm(),
      Catch::Matchers::WithinRel(fb2.getResiduals().segment<3>(3 * 1).norm()));

    Eigen::ArrayXi expectedA =
      Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(net2.nrOfSprings));
    expectedA << 0, 1, 1, 1, 2, 3, 3, 2;
    CHECK(net2.springIndexA.isApprox(expectedA));
    Eigen::ArrayXi expectedB =
      Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(net2.nrOfSprings));
    expectedB << 1, 0, 0, 2, 3, 3, 2, 0;
    CHECK(net2.springIndexB.isApprox(expectedB));
    Eigen::VectorXd springVectors2 =
      net2.evaluateSpringVectors(fb2.getCurrentDisplacements());
    for (const Eigen::Index springIdx : { 0, 1, 2 }) {
      CHECK(springVectors2.segment(springIdx * 3, 3)
              .cwiseAbs()
              .isApproxToConstant(2. - 0.));
    }
    CHECK(springVectors2.segment(5 * 3, 3).isApproxToConstant(10.));
    CHECK(springVectors2.segment(6 * 3, 3).isApproxToConstant(23.5 - 14.5));
    CHECK(springVectors2.segment(7 * 3, 3).isApproxToConstant(40. - 23.5));
    Eigen::Vector3d relevantStrandVector =
      net2.evaluateStrandVector(fb2.getCurrentDisplacements(), 3);
    CHECK(relevantStrandVector.isApprox(
      Eigen::Vector3d(40. - 2., 40. - 2., 40. - 2.)));

    CHECK_THAT(
      fb1.getGammaFactors(1.0).mean(),
      Catch::Matchers::WithinRel(fb2.getGammaFactors(1.0).mean(), 1e-5));
    CHECK_THAT(fb1.getResidual(),
               Catch::Matchers::WithinRel(fb2.getResidual(), 1e-5));

    pcm::MEHPSlipLinkSolverConfig slCfg;
    slCfg.maxNrOfSteps = 25000;
    slCfg.xtol = 1e-13;
    slCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    slCfg.inactiveRemovalCutoff = 1e-12;
    slCfg.doInnerIterations = true;
    slCfg.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
    slCfg.oneOverSpringPartitionUpperLimit = -1.0;
    slCfg.slipMode = pcm::SlipMode::FIXED;
    fb1.run(slCfg);
    fb2.run(pcm::MEHPLinearSolverConfig{});

    CHECK_THAT(fb1.getResidual(),
               Catch::Matchers::WithinAbs(fb2.getResidual(), 1e-6));
    Eigen::VectorXd g1 = fb1.getGammaFactors(1.0);
    Eigen::VectorXd g2 = fb2.getGammaFactors(1.0);
    CHECK_THAT(g1.mean(), Catch::Matchers::WithinRel(g2.mean(), 1e-3));

    Eigen::VectorXd finalPositions1 =
      net1.coordinates + fb1.getCurrentDisplacements();
    Eigen::VectorXd finalPositions2 =
      net2.coordinates + fb2.getCurrentDisplacements();
    Eigen::VectorXd springVectorsAfter1 =
      net1.evaluateSpringVectors(fb1.getCurrentDisplacements());
    Eigen::VectorXd springVectorsAfter2 =
      net2.evaluateSpringVectors(fb2.getCurrentDisplacements());
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(net2.nrOfLinks);
         ++i) {
      CHECK(finalPositions1.segment(3 * i, 3).isApproxToConstant(
        finalPositions1[3 * i]));
      CHECK(finalPositions2.segment(3 * i, 3).isApproxToConstant(
        finalPositions2[3 * i]));
      CHECK(springVectorsAfter1.segment(3 * i, 3).isApproxToConstant(
        springVectorsAfter1[3 * i]));
      CHECK(springVectorsAfter2.segment(3 * i, 3).isApproxToConstant(
        springVectorsAfter2[3 * i]));
    }

    CHECK(fb2.getSolubleWeightFraction() == 0.);
    CHECK(fb1.getSolubleWeightFraction() == 0.);

    CHECK(springVectorsAfter2.segment(0 * 3, 3).isApprox(
      -1. * springVectorsAfter2.segment(1 * 3, 3)));
    CHECK(springVectorsAfter2.segment(1 * 3, 3).isApprox(
      springVectorsAfter2.segment(2 * 3, 3)));

    // from Mathematica solution
    CHECK_THAT(finalPositions1[3 * 0] - finalPositions1[3 * 1],
               Catch::Matchers::WithinRel(-0.45977, 1e-3));
    CHECK_THAT(finalPositions1[3 * 1] - finalPositions1[3 * 2],
               Catch::Matchers::WithinRel(-1.14943, 1e-3));
    CHECK_THAT(finalPositions1[3 * 2] - finalPositions1[3 * 3],
               Catch::Matchers::WithinRel(-1.66667, 1e-3));

    CHECK_THAT(finalPositions2[3 * 0] - finalPositions2[3 * 1],
               Catch::Matchers::WithinRel(-0.45977, 1e-3));
    CHECK_THAT(finalPositions2[3 * 1] - finalPositions2[3 * 2],
               Catch::Matchers::WithinRel(-1.14943, 1e-3));
    CHECK_THAT(finalPositions2[3 * 2] - finalPositions2[3 * 3],
               Catch::Matchers::WithinRel(-1.66667, 1e-3));

    CHECK_NOTHROW(fb2.validateNetwork());
  }

  SECTION("Entanglements through the box, direct copy of network")
  {
    entanglements.pairsOfAtoms.push_back(std::make_pair(3, 24));
    entanglements.pairsOfAtoms.push_back(std::make_pair(5, 14));
    entanglements.pairOfAtom[3] = 0;
    entanglements.pairOfAtom[24] = 0;
    entanglements.pairOfAtom[5] = 1;
    entanglements.pairOfAtom[14] = 1;

    pcm::MEHPSimulator fb1 = pcm::MEHPSimulator(universe, entanglements);
    pcm::MEHPNetwork net1 = fb1.getNetwork();
    Eigen::VectorXd partitions = fb1.getSpringPartitions();

    pcm::MEHPSimulator fb2 = pcm::MEHPSimulator(net1);
    pcm::MEHPNetwork net2 = fb2.getNetwork();

    CHECK_THAT(
      fb1.getGammaFactors(1.0).mean(),
      Catch::Matchers::WithinRel(fb2.getGammaFactors(1.0).mean(), 1e-5));
    CHECK_THAT(fb1.getResidual(),
               Catch::Matchers::WithinRel(fb2.getResidual(), 1e-5));
    for (size_t i = 0; i < net2.nrOfLinks; ++i) {
      CHECK_THAT(
        fb1.getResiduals().segment<3>(3 * static_cast<Eigen::Index>(i)).norm(),
        Catch::Matchers::WithinRel(
          fb2.getResiduals()
            .segment<3>(3 * static_cast<Eigen::Index>(i))
            .norm(),
          1e-5));
    }
    for (size_t i = 0; i < net2.nrOfSprings; ++i) {
      Eigen::Vector3d springVec1 =
        net1.evaluateSpringVector(fb1.getCurrentDisplacements(), i);
      Eigen::Vector3d springVec2 =
        net2.evaluateSpringVector(fb2.getCurrentDisplacements(), i);
      for (Eigen::Index dir = 0; dir < 3; ++dir) {
        CHECK_THAT(springVec1[dir],
                   Catch::Matchers::WithinRel(springVec2[dir], 1e-5));
      }
    }

    pcm::MEHPSlipLinkSolverConfig slCfg;
    slCfg.maxNrOfSteps = 25000;
    slCfg.xtol = 1e-13;
    slCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    slCfg.inactiveRemovalCutoff = 1e-12;
    slCfg.doInnerIterations = true;
    slCfg.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
    slCfg.oneOverSpringPartitionUpperLimit = -1.0;
    slCfg.slipMode = pcm::SlipMode::FIXED;
    fb1.run(slCfg);
    fb2.run(pcm::MEHPLinearSolverConfig{});

    CHECK_THAT(
      fb1.getGammaFactors(1.0).mean(),
      Catch::Matchers::WithinRel(fb2.getGammaFactors(1.0).mean(), 1e-5));
    CHECK_THAT(fb1.getResidualNorm(),
               Catch::Matchers::WithinAbs(fb2.getResidual(), 1e-5));
    for (size_t i = 0; i < net2.nrOfLinks; ++i) {
      CHECK_THAT(
        fb1.getResiduals().segment<3>(3 * static_cast<Eigen::Index>(i)).norm(),
        Catch::Matchers::WithinAbs(
          fb2.getResiduals()
            .segment<3>(3 * static_cast<Eigen::Index>(i))
            .norm(),
          1e-5));
    }
    for (size_t i = 0; i < net2.nrOfSprings; ++i) {
      Eigen::Vector3d springVec1 =
        net1.evaluateSpringVector(fb1.getCurrentDisplacements(), i);
      Eigen::Vector3d springVec2 =
        net2.evaluateSpringVector(fb2.getCurrentDisplacements(), i);
      for (Eigen::Index dir = 0; dir < 3; ++dir) {
        CHECK_THAT(springVec1[dir],
                   Catch::Matchers::WithinRel(springVec2[dir], 1e-5));
      }
    }
    fb2.setCurrentDisplacements(fb1.getCurrentDisplacements());
    CHECK(
      fb2.getCurrentDisplacements().isApprox(fb1.getCurrentDisplacements()));
    CHECK(fb1.getSpringPartitions().isApprox(partitions));
    CHECK(fb2.getNetwork().coordinates.isApprox(fb1.getNetwork().coordinates));
    CHECK(fb2.getNetwork().springBoxOffset.isApprox(
      fb1.getNetwork().springBoxOffset));
    CHECK(
      net2.evaluateSpringVectors(fb2.getCurrentDisplacements())
        .isApprox(net1.evaluateSpringVectors(fb1.getCurrentDisplacements())));
    CHECK_THAT(
      fb1.getGammaFactors(1.0).mean(),
      Catch::Matchers::WithinRel(fb2.getGammaFactors(1.0).mean(), 1e-5));
    CHECK_NOTHROW(fb2.validateNetwork());
  }
}

TEST_CASE("MEHPFB2 Conversion test of small grid",
          "[analysis][MEHPSimulator][StructureConversion]")
{
  std::cout << "Running test \"MEHPFB2 Conversion test of small grid\""
            << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  std::vector<double> coordsX;
  std::vector<double> coordsY;
  std::vector<long int> ids;
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;

  for (size_t y = 0; y < 2; ++y) {
    for (size_t i = 0; i < 10; ++i) {
      coordsX.push_back(static_cast<double>(i));
      coordsY.push_back(static_cast<double>(y) * 5.);
      assert(ids.size() == y * 10 + i);
      ids.push_back(static_cast<long int>(y * 10 + i));
      bondFrom.push_back(static_cast<long int>(y * 10 + i));
      bondTo.push_back(static_cast<long int>(y * 10 + ((i + 1) % 10)));
    }
  }
  for (size_t x = 0; x < 2; ++x) {
    for (size_t i = 1; i < 10; ++i) {
      if (i != 5) {
        coordsX.push_back(static_cast<double>(x) * 5.);
        coordsY.push_back(static_cast<double>(i));
        ids.push_back(static_cast<long int>(ids.size()));
      }
    }
  }
  std::vector<size_t> additionalBonds = {
    0,  20, 20, 21, 21, 22, 22, 23, 23, 10,
    10, 24, 24, 25, 25, 26, 26, 27, 27, 0, // "lower" line
    5,  28, 28, 29, 29, 30, 30, 31, 31, 15,
    15, 32, 32, 33, 33, 34, 34, 35, 35, 5 // "upper" line
  };
  for (size_t i = 0; i < additionalBonds.size(); ++i) {
    if (i % 2 == 0) {
      bondFrom.push_back(static_cast<long int>(additionalBonds[i]));
    } else {
      bondTo.push_back(static_cast<long int>(additionalBonds[i]));
    }
  }
  std::vector<int> types =
    pylimer_tools::utils::initializeWithValue(ids.size(), 1);
  types[0] = 2;
  types[5] = 2;
  types[10] = 2;
  types[15] = 2;
  coordsX[0] += 0.1;
  coordsX[5] -= 0.1;
  coordsX[10] += 0.2;
  coordsX[15] -= 0.05;

  universe.addAtoms(
    ids,
    types,
    coordsX,
    coordsY,
    pylimer_tools::utils::initializeWithValue(coordsX.size(), 0.),
    pylimer_tools::utils::initializeWithValue(ids.size(), 0),
    pylimer_tools::utils::initializeWithValue(ids.size(), 0),
    pylimer_tools::utils::initializeWithValue(ids.size(), 0));

  universe.addBonds(bondFrom, bondTo);

  CHECK(universe.getNrOfBonds() == bondFrom.size());
  CHECK(universe.getNrOfAtoms() == ids.size());

  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements;
  entanglements.pairOfAtom =
    pylimer_tools::utils::initializeWithValue<long int>(ids.size(), -1);

  SECTION("Phantom")
  {
    pcm::MEHPSimulator fb1_phantom =
      pcm::MEHPSimulator(universe, entanglements);
    pcm::MEHPSimulator fb2_phantom =
      pcm::MEHPSimulator(universe, entanglements);
    CHECK_THAT(fb1_phantom.getGammaFactors(1.0).mean(),
               Catch::Matchers::WithinAbs(
                 fb2_phantom.getGammaFactors(1.0).mean(), 1e-6));
    CHECK_THAT(fb1_phantom.getResidual(),
               Catch::Matchers::WithinRel(fb2_phantom.getResidual(), 1e-6));

    pcm::MEHPSlipLinkSolverConfig slCfg;
    slCfg.slipMode = pcm::SlipMode::FIXED;
    fb1_phantom.run(slCfg);
    fb2_phantom.run(pcm::MEHPLinearSolverConfig{});

    CHECK_THAT(fb1_phantom.getResidual(), Catch::Matchers::WithinAbs(0., 1e-6));
    CHECK_THAT(fb2_phantom.getResidual(), Catch::Matchers::WithinAbs(0., 1e-6));
    CHECK_THAT(fb1_phantom.getGammaFactors(1.0).mean(),
               Catch::Matchers::WithinAbs(
                 fb2_phantom.getGammaFactors(1.0).mean(), 1e-6));

    pcm::MEHPNetwork net2 = fb2_phantom.getNetwork();
    pcm::MEHPNetwork net1 = fb1_phantom.getNetwork();

    CHECK(net1.nrOfLinks == 4);
    CHECK(net2.nrOfLinks == 4);
    CHECK(net1.nrOfNodes == 4);
    CHECK(net2.nrOfNodes == 4);
    CHECK(net2.getNeighbourLinkIndices(0).size() == 4);
    CHECK(net2.getNeighbourLinkIndices(1).size() == 4);
    CHECK(net2.getNeighbourLinkIndices(2).size() == 4);
    CHECK(net2.getNeighbourLinkIndices(3).size() == 4);
    CHECK_FALSE(net1.springBoxOffset.isZero());
    CHECK_FALSE(net2.springBoxOffset.isZero());
    CHECK(net2.nrOfStrands == 8);
    CHECK(net2.nrOfSprings == net2.nrOfStrands);
    CHECK(net1.nrOfSprings == 8);
    CHECK(fb1_phantom.getNrOfActiveSprings() == net1.nrOfSprings);
    CHECK(fb2_phantom.getNrOfActiveSprings() == net2.nrOfSprings);
    CHECK(fb2_phantom.getNrOfActiveStrands() == net2.nrOfStrands);
    CHECK(fb2_phantom.getNrOfActiveNodes() == net2.nrOfNodes);
    CHECK(net2.springContourLength.isApproxToConstant(5.));

    Eigen::VectorXd displacements = fb2_phantom.getCurrentDisplacements();
    CHECK(net2.breakTooLongStrands(displacements, 100.) == 0);
    CHECK(net2.breakTooLongStrands(displacements, 1e-9) == 8);
    CHECK(net2.nrOfSprings == 0);
    CHECK_NOTHROW(net2.validate(displacements));
    CHECK_NOTHROW(net2.displaceToMeanPosition(displacements, size_t{ 0 }));
    CHECK_NOTHROW(net2.displaceToMeanPosition(displacements));
  }

  SECTION("Entangled")
  {
    entanglements.pairsOfAtoms.push_back(std::make_pair(3, 13));
    entanglements.pairsOfAtoms.push_back(std::make_pair(16, 32));
    entanglements.pairOfAtom[3] = 0;
    entanglements.pairOfAtom[13] = 0;
    entanglements.pairOfAtom[16] = 1;
    entanglements.pairOfAtom[32] = 1;

    pcm::MEHPSimulator fb1 = pcm::MEHPSimulator(universe, entanglements);
    pcm::MEHPSimulator fb2 = pcm::MEHPSimulator(universe, entanglements);

    pcm::MEHPNetwork net2 = fb2.getNetwork();
    pcm::MEHPNetwork net1 = fb1.getNetwork();

    CHECK(fb2.getCurrentDisplacements().isZero());
    CHECK(fb1.getCurrentDisplacements().isZero());
    CHECK(net1.nrOfNodes == 4);
    CHECK(net2.nrOfNodes == 4);
    CHECK(fb2.getNumParticles() == fb2.getNumAtoms());
    CHECK(fb2.getNumAtoms() == 4);
    CHECK(fb2.getNumParticles() == net1.nrOfNodes);
    CHECK(fb2.getNrOfActiveSprings() == net2.nrOfSprings);
    CHECK(fb1.getNrOfActiveSprings() == net1.nrOfSprings);

    CHECK_FALSE(net1.springBoxOffset.isZero());
    CHECK_FALSE(net2.springBoxOffset.isZero());
    REQUIRE(net1.nrOfSprings == net2.nrOfSprings);
    CHECK(net2.nrOfSprings == 12);
    CHECK(net2.nrOfStrands == 8);

    for (size_t strandIdx = 0; strandIdx < net2.nrOfStrands; ++strandIdx) {
      Eigen::Vector3d strandVector =
        net2.evaluateStrandVector(fb2.getCurrentDisplacements(), strandIdx);
      CHECK_THAT(strandVector.squaredNorm(),
                 Catch::Matchers::WithinAbs(5 * 5. + 0 + 0, 3));
    }
    CHECK_THAT(net2.springBoxOffset.cwiseAbs().maxCoeff(),
               Catch::Matchers::WithinRel(10., 1e-6));
    CHECK_THAT(
      fb1.getGammaFactors(1.0).mean(),
      Catch::Matchers::WithinRel(fb2.getGammaFactors(1.0).mean(), 1e-2));
    CHECK_THAT(fb1.getResidualNorm(),
               Catch::Matchers::WithinRel(fb2.getResidual(), 1e-1));

    pcm::MEHPSlipLinkSolverConfig slCfg;
    slCfg.maxNrOfSteps = 25000;
    slCfg.xtol = 1e-13;
    slCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    slCfg.inactiveRemovalCutoff = 1e-12;
    slCfg.doInnerIterations = true;
    slCfg.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
    slCfg.oneOverSpringPartitionUpperLimit = -1.0;
    slCfg.slipMode = pcm::SlipMode::FIXED;
    fb1.run(slCfg);
    fb2.run(pcm::MEHPLinearSolverConfig{});

    CHECK_THAT(fb1.getResidualNorm(),
               Catch::Matchers::WithinAbs(fb2.getResidual(), 1e-6));
    Eigen::VectorXd g1 = fb1.getGammaFactors(1.0);
    Eigen::VectorXd g2 = fb2.getGammaFactors(1.0);
    CHECK_THAT(g1.mean(), Catch::Matchers::WithinRel(g2.mean(), 1e-3));
    CHECK(fb2.getNrOfActiveSprings() == net2.nrOfSprings);
    CHECK(fb1.getNrOfActiveSprings() == net1.nrOfSprings);

    CHECK(fb2.getSolubleWeightFraction() == 0.);
    CHECK(fb1.getSolubleWeightFraction() == 0.);
  }

  SECTION("Entangled, direct copy of network")
  {
    entanglements.pairsOfAtoms.push_back(std::make_pair(3, 13));
    entanglements.pairsOfAtoms.push_back(std::make_pair(16, 32));
    entanglements.pairOfAtom[3] = 0;
    entanglements.pairOfAtom[13] = 0;
    entanglements.pairOfAtom[16] = 1;
    entanglements.pairOfAtom[32] = 1;

    pcm::MEHPSimulator fb1 = pcm::MEHPSimulator(universe, entanglements);
    pcm::MEHPNetwork net1 = fb1.getNetwork();

    pcm::MEHPSimulator fb2 = pcm::MEHPSimulator(net1);

    CHECK(fb2.getResiduals().segment<3>(0).isApprox(
      fb1.getResiduals().segment<3>(0)));
    CHECK(fb2.getNetwork()
            .evaluateStressTensorForLinks(
              { 0 }, fb2.getCurrentDisplacements(), 1.0, false, -1.0)
            .isApprox(fb1.getNetwork().evaluateStressTensorForLinks(
              { 0 }, fb1.getCurrentDisplacements(), 1.0, false, -1.0)));
    std::vector<size_t> allLinks1(fb1.getNrOfLinks());
    std::iota(allLinks1.begin(), allLinks1.end(), 0);
    std::vector<size_t> allLinks2(fb2.getNrOfLinks());
    std::iota(allLinks2.begin(), allLinks2.end(), 0);
    CHECK(fb2.getNetwork()
            .evaluateStressTensorForLinks(
              allLinks2, fb2.getCurrentDisplacements(), 1.0, false, -1.0)
            .isApprox(fb1.getNetwork().evaluateStressTensorForLinks(
              allLinks1, fb1.getCurrentDisplacements(), 1.0, false, -1.0)));
    CHECK(fb2.getStressTensor().isApprox(fb1.getStressTensor()));
    CHECK_THAT(
      fb1.getGammaFactors(1.0).mean(),
      Catch::Matchers::WithinRel(fb2.getGammaFactors(1.0).mean(), 1e-3));
    CHECK_THAT(fb1.getResidualNorm(),
               Catch::Matchers::WithinRel(fb2.getResidual(), 1e-3));

    pcm::MEHPSlipLinkSolverConfig slCfg;
    slCfg.maxNrOfSteps = 25000;
    slCfg.xtol = 1e-13;
    slCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    slCfg.inactiveRemovalCutoff = 1e-12;
    slCfg.doInnerIterations = true;
    slCfg.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
    slCfg.oneOverSpringPartitionUpperLimit = -1.0;
    slCfg.slipMode = pcm::SlipMode::FIXED;
    fb1.run(slCfg);
    fb2.run(pcm::MEHPLinearSolverConfig{});

    CHECK_THAT(fb1.getResidualNorm(),
               Catch::Matchers::WithinAbs(fb2.getResidual(), 1e-6));
    Eigen::VectorXd g1 = fb1.getGammaFactors(1.0);
    Eigen::VectorXd g2 = fb2.getGammaFactors(1.0);
    CHECK_THAT(g1.mean(), Catch::Matchers::WithinRel(g2.mean(), 1e-3));
  }
}

TEST_CASE("MEHPFB2 Basic deformation test", "[analysis][MEHPSimulator][Box]")
{
  std::cout << "Running test \"MEHPFB2 Basic deformation test\"" << std::endl;

  const std::string suspectedPath = std::string(PYLIMER_TEST_FIXTURES_DIR);
  REQUIRE(std::filesystem::exists(suspectedPath));

  SECTION("Non-empty universe")
  {
    std::string inputFile =
      suspectedPath +
      "/structure/3d-diamond-lattice_3x3x3_a_23_d_3_v_0.structure.out";
    REQUIRE(std::filesystem::exists(inputFile));
    pe::UniverseSequence universeSequence = pe::UniverseSequence();
    REQUIRE(universeSequence.getLength() == 0);
    universeSequence.initializeFromDataSequence({ { inputFile } });
    REQUIRE(universeSequence.getLength() == 1);
    pe::Universe universe = universeSequence.atIndex(0);

    pcm::MEHPSimulator fb2 = pcm::MEHPSimulator(universe);
    double residualBeforeDeformation = fb2.getResidual();
    pe::Box prevBox = universe.getBox();
    pe::Box newBox = pe::Box(prevBox.getLowX() * 0.9,
                             prevBox.getHighX() * 1.1,
                             prevBox.getLowY(),
                             prevBox.getHighY(),
                             prevBox.getLowZ(),
                             prevBox.getHighZ());
    fb2.deformTo(newBox);
    double residualAfterDeformation = fb2.getResidual();
    CHECK(residualBeforeDeformation < residualAfterDeformation);
  }
}

TEST_CASE("MEHPFB2 Basic conversion test with entanglements",
          "[analysis][MEHPSimulator][StructureConversion]")
{
  std::cout
    << "Running test \"MEHPFB2 Basic conversion test with entanglements\""
    << std::endl;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);

  std::vector<long int> ids = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  std::vector<int> types =
    pylimer_tools::utils::initializeWithValue(ids.size(), 1);
  std::vector<double> coordsX = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9. };
  std::vector<int> ns =
    pylimer_tools::utils::initializeWithValue(ids.size(), 0);

  universe.addAtoms(ids, types, coordsX, coordsX, coordsX, ns, ns, ns);
  universe.addBonds({ 0, 1, 2, 3, 4, 5, 6, 7, 8 },
                    { 1, 2, 3, 4, 5, 6, 7, 8, 9 });

  pylimer_tools::topo::entanglement_detection::AtomPairEntanglements
    entanglements;
  entanglements.pairOfAtom =
    pylimer_tools::utils::initializeWithValue<long int>(ids.size(), -1);
  entanglements.pairsOfAtoms.push_back(std::make_pair(1, 5));
  entanglements.pairOfAtom[1] = 0;
  entanglements.pairOfAtom[5] = 0;
  entanglements.pairsOfAtoms.push_back(std::make_pair(2, 4));
  entanglements.pairOfAtom[2] = 1;
  entanglements.pairOfAtom[4] = 1;

  SECTION("Entanglements as links")
  {
    pcm::MEHPSimulator fb2 = pcm::MEHPSimulator(
      universe, entanglements, 2, pcm::EntanglementMode::MERGED_LINKS);
    CHECK_NOTHROW(fb2.validateNetwork());

    pcm::MEHPNetwork net = fb2.getNetwork();
    CHECK(net.nrOfLinks == 4);
    CHECK(net.nrOfNodes == 2);
    CHECK(net.nrOfStrands == 1);
    CHECK(net.nrOfSprings == 5);
    CHECK(net.oldAtomIds[0] == 0);
    CHECK(net.oldAtomIds[1] == -1);
    CHECK(net.oldAtomIds[2] == -1);
    CHECK(net.oldAtomIds[3] == 9);
    CHECK(net.linkIsEntanglement[1]);
    CHECK(net.linkIsEntanglement[2]);
    CHECK(fb2.getNumExtraAtoms() == entanglements.pairsOfAtoms.size());
  }

  SECTION("Entanglements as springs")
  {
    pcm::MEHPSimulator fb2 = pcm::MEHPSimulator(
      universe, entanglements, 2, pcm::EntanglementMode::SPRING);
    CHECK_NOTHROW(fb2.validateNetwork());

    pcm::MEHPNetwork net = fb2.getNetwork();
    CHECK(net.nrOfLinks == 6);
    CHECK(net.nrOfNodes == 2);
    CHECK(net.nrOfStrands == 3);
    CHECK(net.nrOfSprings == 7);
    CHECK(net.oldAtomIds[0] == 0);
    CHECK(net.oldAtomIds[1] == -1);
    CHECK(net.oldAtomIds[2] == -1);
    CHECK(net.oldAtomIds[3] == -1);
    CHECK(net.oldAtomIds[4] == -1);
    CHECK(net.oldAtomIds[5] == 9);
    CHECK(net.linkIsEntanglement[1]);
    CHECK(net.linkIsEntanglement[2]);
    CHECK(net.linkIsEntanglement[3]);
    CHECK(net.linkIsEntanglement[4]);
    CHECK(fb2.getNumExtraBonds() == entanglements.pairsOfAtoms.size());
    CHECK(fb2.getNumExtraAtoms() == entanglements.pairsOfAtoms.size() * 2);
  }
}

TEST_CASE("Selective MC Generated Structure Solves",
          "[MEHPSimulator][MCUniverseGenerator][long]")
{
  double scaleDownBy = 1. / 10.;

  pylimer_tools::utils::MCUniverseGenerator generator =
    pylimer_tools::utils::MCUniverseGenerator(
      21.118005279044613 / std::cbrt(scaleDownBy),
      21.118005279044613 / std::cbrt(scaleDownBy),
      21.118005279044613 / std::cbrt(scaleDownBy));
  generator.setMeanSquaredBeadDistance(1.0647);
  generator.addRandomlyFunctionalizedStrands(
    2,
    pylimer_tools::utils::initializeWithValue(
      2, static_cast<int>(217 / scaleDownBy)),
    7 * 0.9407504937458855,
    1,
    2,
    1);
  int nStrands = static_cast<int>(886 / scaleDownBy);
  generator.addStrands(static_cast<size_t>(nStrands),
                       pylimer_tools::utils::initializeWithValue(
                         static_cast<size_t>(nStrands), 318),
                       1);
  generator.useZScoreMaxDistance(2.);
  CHECK_NOTHROW(generator.linkStrandsToConversion(0.12327785585899244));

  pcm::MEHPSimulator fb2 = generator.getMEHPSimulator();
  CHECK_NOTHROW(fb2.validateNetwork());
  pcm::MEHPLinearSolverConfig cfg;
  cfg.simplificationMode = pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
  cfg.inactiveRemovalCutoff = 1e-5;
  cfg.solver = pcm::SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG;
  fb2.run(cfg);
  CHECK(fb2.getPressure() < 1e-3);
}

TEST_CASE("All MEHP Force Balance2 vs. Force Relaxation Phantom Comparisons",
          "[analysis][MEHPSimulator][long]")
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
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure."
    "out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_1_1_melt_100_a_3_50_xlinks_v_1.V-fixed.structure.out-"
    "finish_"
    "crosslinking.structure.out",
    "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0.333_2d_t_7500001."
    "structure.out",
    "mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0.0114_f_4_v_1."
    "structure.out",
    "melt_213_a_47_106_xlinks_v_1.structure.out",
    "melt_83_a_100.structure.out",
    "network_100_a_46.structure.out",
    "network_p_1_100_a_38_50_xlinks.structure.out",
    "square_lattice_2x2_a_5.2d.structure.out",
    "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out",
    "xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out"
  };

#define CHECK_THAT_OR_ZERO(expr, otherExpr)                                    \
  if (std::abs(expr) < 1e-5) {                                                 \
    CHECK(std::abs(otherExpr) < 2e-5);                                         \
  } else {                                                                     \
    CHECK_THAT(                                                                \
      expr,                                                                    \
      Catch::Matchers::WithinRel(otherExpr, 1e-2)); /** allow 1% deviation */  \
  }

  size_t nFilesFound = 0;
  for (const std::string& file : files) {
    INFO("Processing file: " << file);
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

    pcm::MEHPSimulator forceBalance2 = pcm::MEHPSimulator(universe);

    pcm::MEHPSimulator forceRelaxation =
      pcm::MEHPSimulator(universe,
                         {},
                         2,
                         pcm::EntanglementMode::NONE,
                         nullptr,
                         false,
                         true,
                         false,
                         false);

    auto start_fb2 = std::chrono::high_resolution_clock::now();
    forceBalance2.run(pcm::MEHPLinearSolverConfig{});
    auto end_fb2 = std::chrono::high_resolution_clock::now();
    auto duration_fb2 = std::chrono::duration_cast<std::chrono::microseconds>(
      end_fb2 - start_fb2);
    std::cout << "Time of Force Balance 2: "
              << std::duration_to_string(duration_fb2) << " " << std::endl;

    pcm::MEHPNonlinearSolverConfig nlCfg;
    nlCfg.algorithm = pcm::NLoptAlgorithm::LD_MMA;
    nlCfg.maxNrOfSteps = 50000;
    nlCfg.xtol = 1e-17;
    nlCfg.ftol = 1e-17;
    auto start_fr = std::chrono::high_resolution_clock::now();
    forceRelaxation.run(nlCfg);
    auto end_fr = std::chrono::high_resolution_clock::now();
    auto duration_fr =
      std::chrono::duration_cast<std::chrono::microseconds>(end_fr - start_fr);
    std::cout << "Time of Force Relaxation: "
              << std::duration_to_string(duration_fr) << " " << std::endl;

    CHECK_THAT_OR_ZERO(forceBalance2.getGammaFactors(1.).sum(),
                       forceRelaxation.getGammaFactors(1.).sum());
    CHECK_THAT_OR_ZERO(forceBalance2.getStressTensor().trace(),
                       forceRelaxation.getStressTensor().trace());
    CHECK_THAT_OR_ZERO(forceBalance2.getResidual(),
                       forceRelaxation.getResidual());
    CHECK_THAT_OR_ZERO(forceBalance2.getSolubleWeightFraction(),
                       forceRelaxation.getSolubleWeightFraction());
    CHECK_THAT_OR_ZERO(forceBalance2.getDanglingWeightFraction(),
                       forceRelaxation.getDanglingWeightFraction());
  }

#undef CHECK_THAT_OR_ZERO
  REQUIRE(static_cast<double>(nFilesFound) >
          static_cast<double>(files.size()) * 0.75);
}

TEST_CASE("All MEHP Force Balance 1 vs. 2 Comparisons with Entanglements and "
          "Simplification",
          "[analysis][MEHPSimulator][long]")
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
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure."
    "out-"
    "equilibration_do_crosslink.structure.out",
    "crosslinked_p_1_1_melt_100_a_3_50_xlinks_v_1.V-fixed.structure.out-"
    "finish_"
    "crosslinking.structure.out",
    "equil_phantom_hexa_lattice_60x60_25_bx_sqrtNbsqrt0.333_2d_t_7500001."
    "structure.out",
    "mc_own-si_pdms_crosslinked_melt_464_a_77_r_1.71_wsol_0.0114_f_4_v_1."
    "structure.out",
    "melt_213_a_47_106_xlinks_v_1.structure.out",
    "melt_83_a_100.structure.out",
    "network_100_a_46.structure.out",
    "network_p_1_100_a_38_50_xlinks.structure.out",
    "square_lattice_2x2_a_5.2d.structure.out",
    "xlinked_0.90005_pdms_1e4_a_78_bs_t_775036.structure.out",
    "xlinked_1e4_a_28_f_3_p_0.151515151515152.structure.out"
  };

  size_t nFilesFound = 0;
  for (const std::string& file : files) {
    INFO("Processing file: " << file);
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

    auto entanglements =
      pylimer_tools::topo::entanglement_detection::randomlyFindEntanglements(
        universe,
        static_cast<size_t>(0.005 *
                            static_cast<double>(universe.getNrOfAtoms())),
        2.5,
        0.,
        0,
        2,
        "t35ts33d",
        2,
        true,
        true);

    pcm::MEHPSimulator forceBalance2 = pcm::MEHPSimulator(
      universe, entanglements, 2, pcm::EntanglementMode::MERGED_LINKS);

    pcm::MEHPSimulator forceBalance = pcm::MEHPSimulator(
      universe, entanglements, 2, pcm::EntanglementMode::MERGED_LINKS);

#define CHECK_THAT_OR_ZERO(expr, otherExpr)                                    \
  if (std::abs(expr) < 1e-6) {                                                 \
    CHECK(std::abs(otherExpr) < 1e-6);                                         \
  } else {                                                                     \
    CHECK_THAT(expr, Catch::Matchers::WithinRel(otherExpr, 1e-3));             \
  }

    // validate initial values
    CHECK_THAT_OR_ZERO(forceBalance2.getGammaFactors(1.).sum(),
                       forceBalance.getGammaFactors(1.).sum());
    CHECK_THAT_OR_ZERO(forceBalance2.getStressTensor().trace(),
                       forceBalance.getStressTensor().trace());
    CHECK_THAT_OR_ZERO(forceBalance2.getResidual(), forceBalance.getResidual());

    pcm::MEHPSlipLinkSolverConfig slCfg;
    slCfg.maxNrOfSteps = 5000;
    slCfg.xtol = 1e-9 / std::max(forceBalance.getResidual(), 1.0);
    slCfg.initialResidualToUse = -1.;
    slCfg.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
    slCfg.inactiveRemovalCutoff = 1e-5;
    slCfg.doInnerIterations = false;
    slCfg.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
    slCfg.oneOverSpringPartitionUpperLimit = 1.0;
    slCfg.slipMode = pcm::SlipMode::FIXED;

    auto start_fr = std::chrono::high_resolution_clock::now();
    forceBalance.run(slCfg);
    auto end_fr = std::chrono::high_resolution_clock::now();
    auto duration_fr =
      std::chrono::duration_cast<std::chrono::microseconds>(end_fr - start_fr);
    std::cout << "Time of Force Balance 1: "
              << std::duration_to_string(duration_fr) << " " << std::endl;

    pcm::MEHPLinearSolverConfig linCfg;
    linCfg.simplificationMode =
      pcm::StructureSimplificationMode::INACTIVE_THEN_X2F;
    auto start_fb2 = std::chrono::high_resolution_clock::now();
    forceBalance2.run(linCfg);
    auto end_fb2 = std::chrono::high_resolution_clock::now();
    auto duration_fb2 = std::chrono::duration_cast<std::chrono::microseconds>(
      end_fb2 - start_fb2);
    std::cout << "Time of Force Balance 2: "
              << std::duration_to_string(duration_fb2) << " " << std::endl;

    if (forceBalance.getExitReason() != pcm::ExitReason::MAX_STEPS) {
      CHECK_THAT_OR_ZERO(forceBalance2.getGammaFactors(1.).sum(),
                         forceBalance.getGammaFactors(1.).sum());
      CHECK_THAT_OR_ZERO(forceBalance2.getStressTensor().trace(),
                         forceBalance.getStressTensor().trace());
      CHECK_THAT_OR_ZERO(forceBalance2.getResidual(),
                         forceBalance.getResidual());
      CHECK_THAT_OR_ZERO(forceBalance2.getSolubleWeightFraction(),
                         forceBalance.getSolubleWeightFraction());
      CHECK_THAT_OR_ZERO(forceBalance2.getDanglingWeightFraction(),
                         forceBalance.getDanglingWeightFraction());
      CHECK_THAT_OR_ZERO(forceBalance2.getActiveWeightFraction(),
                         forceBalance.getActiveWeightFraction());
    } else {
      std::cerr << "Force Balance 1 did not converge for file: " << inputFile
                << std::endl;
    }

    pcm::MEHPSimulator forceBalance21(forceBalance.getNetwork());
    forceBalance21.run(linCfg);
    CHECK_THAT_OR_ZERO(forceBalance21.getGammaFactors(1.).sum(),
                       forceBalance2.getGammaFactors(1.).sum());

    forceBalance21 = pcm::MEHPSimulator(forceBalance.getNetwork());
    pcm::MEHPLinearSolverConfig noSimpLinCfg;
    noSimpLinCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    forceBalance21.run(noSimpLinCfg);
    CHECK_THAT_OR_ZERO(forceBalance21.getGammaFactors(1.).sum(),
                       forceBalance2.getGammaFactors(1.).sum());

#undef CHECK_THAT_OR_ZERO
  }

  REQUIRE(static_cast<double>(nFilesFound) >
          static_cast<double>(files.size()) * 0.75);
}

TEST_CASE("All MEHPSimulator solvers solve simple melt systems",
          "[MEHPSimulator]")
{
  std::cout << "Running test \"All MEHPSimulator solvers solve simple melt "
               "systems\""
            << std::endl;

  size_t nrOfBeads = 30;
  size_t nrOfChains = 5;
  pe::Universe universe =
    pe::Universe(nrOfBeads * 10.0, nrOfBeads * 10.0, nrOfBeads * 10.0);
  for (size_t i = 0; i < nrOfChains; ++i) {
    long int beadOffset = static_cast<long int>(i * nrOfBeads);
    std::vector<long int> atomIds(nrOfBeads);
    std::iota(atomIds.begin(), atomIds.end(), 1 + beadOffset);
    std::vector<double> coords(nrOfBeads);
    std::iota(
      coords.begin(), coords.end(), 1. + static_cast<double>(beadOffset));
    universe.addAtoms(atomIds,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 1),
                      coords,
                      coords,
                      coords,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0));
    std::vector<long int> bondFrom(nrOfBeads - 1);
    std::vector<long int> bondTo(nrOfBeads - 1);
    std::iota(bondFrom.begin(), bondFrom.end(), 1 + beadOffset);
    std::iota(bondTo.begin(), bondTo.end(), 2 + beadOffset);
    universe.addBonds(bondFrom, bondTo);
  }
  CHECK(universe.getNrOfAtoms() == nrOfBeads * nrOfChains);
  CHECK(universe.getNrOfBonds() == (nrOfBeads - 1) * nrOfChains);
  universe.addBonds(
    { 1, 2, 1 }, { 1, 1, static_cast<long int>(nrOfChains * nrOfBeads - 1) });

  for (pcm::SLESolver solver : {
         pcm::SLESolver::DEFAULT,
         pcm::SLESolver::SIMPLICIAL_LLT,
         pcm::SLESolver::SIMPLICIAL_LDLT,
         pcm::SLESolver::SPARSE_LU,
         pcm::SLESolver::SPARSE_QR,
#ifdef __APPLE__
         pcm::SLESolver::ACCELERATE_LLT,
         pcm::SLESolver::ACCELERATE_LDLT,
#endif
         pcm::SLESolver::GRADIENT_DESCENT,
         pcm::SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT,
         pcm::SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG,
         pcm::SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM,
#ifdef USE_METAL
         pcm::SLESolver::METAL_BARZILAI_BORWEIN,
         pcm::SLESolver::METAL_CG,
#endif
       }) {
#ifdef USE_METAL
    if ((solver == pcm::SLESolver::METAL_BARZILAI_BORWEIN ||
         solver == pcm::SLESolver::METAL_CG) &&
        !pcm::isMetalAvailable()) {
      continue;
    }
#endif
    INFO("Testing solver: " << solver << " ("
                            << pcm::SLESolverNames[static_cast<int>(solver)]
                            << ")");
    pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    CHECK(forceBalancer.getNrOfSprings() == nrOfChains + 5);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = solver;
    CHECK_NOTHROW(forceBalancer.run(cfg));
    CHECK(forceBalancer.getNrOfSprings() == nrOfChains + 5);
    CHECK(forceBalancer.getNrOfActiveSprings() == 0);
  }
}

TEST_CASE("All MEHPSimulator solvers solve simple network systems",
          "[MEHPSimulator]")
{
  std::cout << "Running test \"All MEHPSimulator solvers solve simple network "
               "systems\""
            << std::endl;

  size_t nrOfBeads = 30;
  size_t nrOfChains = 5;
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  for (size_t i = 0; i < nrOfChains; ++i) {
    long int beadOffset = static_cast<long int>(i * nrOfBeads);
    std::vector<long int> atomIds(nrOfBeads);
    std::iota(atomIds.begin(), atomIds.end(), 1 + beadOffset);
    std::vector<double> coords(nrOfBeads);
    std::iota(
      coords.begin(), coords.end(), 1. + static_cast<double>(beadOffset));
    universe.addAtoms(atomIds,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 1),
                      coords,
                      coords,
                      coords,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0));
    std::vector<long int> bondFrom(nrOfBeads - 1);
    std::vector<long int> bondTo(nrOfBeads - 1);
    std::iota(bondFrom.begin(), bondFrom.end(), 1 + beadOffset);
    std::iota(bondTo.begin(), bondTo.end(), 2 + beadOffset);
    universe.addBonds(bondFrom, bondTo);
  }
  CHECK(universe.getNrOfAtoms() == nrOfBeads * nrOfChains);
  CHECK(universe.getNrOfBonds() == (nrOfBeads - 1) * nrOfChains);
  universe.addBonds(
    { 1, 2, 1, 4, 29 },
    { 1,
      1,
      static_cast<long int>(nrOfChains * nrOfBeads - 1),
      14,
      static_cast<long int>((nrOfChains - 1) * nrOfBeads + 1) });

  for (pcm::SLESolver solver : {
         pcm::SLESolver::DEFAULT,
         pcm::SLESolver::SPARSE_QR,
         pcm::SLESolver::CONJUGATE_GRADIENT,
         pcm::SLESolver::CONJUGATE_GRADIENT_IDENTITY,
         pcm::SLESolver::CONJUGATE_GRADIENT_DIAGONALIZED,
         pcm::SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT,
         pcm::SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED,
         pcm::SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY,
         pcm::SLESolver::BICGSTAB,
         pcm::SLESolver::BICGSTAB_DIAGONALIZED,
         pcm::SLESolver::BICGSTAB_IDENTITY,
         pcm::SLESolver::BICGSTAB_INCOMPLETE_LU,
         pcm::SLESolver::GRADIENT_DESCENT,
         pcm::SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT,
         pcm::SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG,
         pcm::SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM,
#ifdef USE_METAL
         pcm::SLESolver::METAL_BARZILAI_BORWEIN,
         pcm::SLESolver::METAL_CG,
#endif
       }) {
#ifdef USE_METAL
    if ((solver == pcm::SLESolver::METAL_BARZILAI_BORWEIN ||
         solver == pcm::SLESolver::METAL_CG) &&
        !pcm::isMetalAvailable()) {
      continue;
    }
#endif
    INFO("Testing solver: " << solver);
    pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    CHECK(forceBalancer.getNrOfSprings() == nrOfChains + 9);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.residualReduction = 1e-8;
    cfg.solver = solver;
    CHECK_NOTHROW(forceBalancer.run(cfg));
    CHECK(forceBalancer.getNrOfSprings() == nrOfChains + 9);
    CHECK(forceBalancer.getNrOfActiveSprings() > 0);
    CHECK(forceBalancer.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    CHECK(forceBalancer.getResidual() < 1e-5);
  }
}

TEST_CASE("All MEHP compute correct dangling and soluble fractions",
          "[MEHPSimulator]")
{
  std::cout << "Running test \"All MEHP compute correct dangling fractions\""
            << std::endl;

  size_t nSolventChains = 2;
  size_t nBeadsPerSolventChain = 5;

  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  universe.addAtoms(
    { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } },
    { { 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1 } },
    { { 0.0, 2.5, 5.0, 7.5, 0.1, 2.5, 5.0, 7.5, -.1, 5.0, 0.0, 5.0 } },
    { { 0.1, 0.0, -.1, 0.0, 5.0, 5.0, 5.0, 5.0, 2.5, 2.5, 7.5, 7.5 } },
    { { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
    { { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } });
  universe.addBonds({ { 1,
                        1,
                        1, //
                        1,
                        3,
                        3, //
                        3,
                        3,
                        5, //
                        5,
                        5,
                        5, //
                        7,
                        7,
                        7, //
                        7 } },
                    { { 2,
                        9,
                        4, //
                        11,
                        2,
                        4, //
                        10,
                        12,
                        9, //
                        11,
                        6,
                        8, //
                        6,
                        8,
                        10, //
                        12 } });
  CHECK(universe.getNrOfAtoms() == 12);
  size_t currentAtomId = 13;

  for (size_t solventChainI = 0; solventChainI < nSolventChains;
       ++solventChainI) {
    std::vector<long int> ids(nBeadsPerSolventChain);
    std::vector<double> coords(nBeadsPerSolventChain, 0.);
    std::vector<long int> bondFrom;
    std::vector<long int> bondTo;
    for (size_t beadI = 0; beadI < nBeadsPerSolventChain; ++beadI) {
      ids[beadI] = static_cast<long int>(currentAtomId);
      currentAtomId += 1;
      if (beadI > 0) {
        bondFrom.push_back(ids[beadI - 1]);
        bondTo.push_back(ids[beadI]);
      }
    }
    universe.addAtoms(ids,
                      pu::initializeWithValue(ids.size(), 3),
                      coords,
                      coords,
                      coords,
                      pu::initializeWithValue(ids.size(), 0),
                      pu::initializeWithValue(ids.size(), 0),
                      pu::initializeWithValue(ids.size(), 0));
    universe.addBonds(bondFrom, bondTo);
  }

  for (size_t nExtraChains :
       { size_t{ 0 }, size_t{ 1 }, size_t{ 3 }, size_t{ 5 } }) {
    for (size_t nDanglingAtoms : { size_t{ 1 },
                                   size_t{ 3 },
                                   size_t{ 5 },
                                   size_t{ 6 },
                                   size_t{ 9 },
                                   size_t{ 14 } }) {
      INFO("Testing with " << nExtraChains << " extra chains and "
                           << nDanglingAtoms << " dangling atoms each");
      pe::Universe u2 = pe::Universe(universe);
      std::vector<long int> bondFrom;
      std::vector<long int> bondTo;
      for (size_t extraChainI = 0; extraChainI < nExtraChains; ++extraChainI) {
        std::vector<long int> ids(nDanglingAtoms);
        long int attachedToId =
          static_cast<long int>((extraChainI * nDanglingAtoms) % 12 + 1);
        std::vector<double> coordX(
          nDanglingAtoms,
          universe.getAtom(attachedToId).getUnwrappedX(universe.getBox()));
        std::vector<double> coordY(
          nDanglingAtoms,
          universe.getAtom(attachedToId).getUnwrappedY(universe.getBox()));
        std::vector<double> coordZ(
          nDanglingAtoms,
          universe.getAtom(attachedToId).getUnwrappedZ(universe.getBox()));
        for (size_t danglingAtomI = 0; danglingAtomI < nDanglingAtoms;
             ++danglingAtomI) {
          if (danglingAtomI == 0) {
            bondFrom.push_back(attachedToId);
          } else {
            bondFrom.push_back(static_cast<long int>(currentAtomId - 1));
          }
          ids[danglingAtomI] = static_cast<long int>(currentAtomId);
          bondTo.push_back(static_cast<long int>(currentAtomId));
          currentAtomId += 1;
          coordX[danglingAtomI] += (std::rand() % 1000) / 1000.0 - 0.5;
          coordY[danglingAtomI] += (std::rand() % 1000) / 1000.0 - 0.5;
          coordZ[danglingAtomI] += (std::rand() % 1000) / 1000.0 - 0.5;
        }

        std::vector<int> atomTypes =
          pylimer_tools::utils::initializeWithValue(nDanglingAtoms, 1);
        atomTypes.back() = 2;
        u2.addAtoms(
          ids,
          atomTypes,
          coordX,
          coordY,
          coordZ,
          pylimer_tools::utils::initializeWithValue(nDanglingAtoms, 0),
          pylimer_tools::utils::initializeWithValue(nDanglingAtoms, 0),
          pylimer_tools::utils::initializeWithValue(nDanglingAtoms, 0));
      }
      u2.addBonds(bondFrom, bondTo);

      size_t totalNDanglingAtoms = nDanglingAtoms * nExtraChains;
      size_t totalNSolventAtoms = nSolventChains * nBeadsPerSolventChain;
      size_t totalNAtoms = 12 + totalNDanglingAtoms + totalNSolventAtoms;
      CHECK(u2.getNrOfAtoms() == totalNAtoms);

      double expectedDanglingFraction =
        static_cast<double>(totalNDanglingAtoms) /
        static_cast<double>(totalNAtoms);
      double expectedSolventFraction = static_cast<double>(totalNSolventAtoms) /
                                       static_cast<double>(totalNAtoms);

      pcm::MEHPSimulator forceBalancer2 = pcm::MEHPSimulator(u2);
      pcm::MEHPLinearSolverConfig linCfg;
      linCfg.simplificationMode =
        pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
      forceBalancer2.run(linCfg);
      CHECK_THAT(forceBalancer2.getDanglingWeightFraction(),
                 Catch::Matchers::WithinRel(expectedDanglingFraction));
      CHECK_THAT(forceBalancer2.getSolubleWeightFraction(),
                 Catch::Matchers::WithinRel(expectedSolventFraction));

      pcm::MEHPSimulator forceBalancer = pcm::MEHPSimulator(u2);
      pcm::MEHPSlipLinkSolverConfig slCfg;
      slCfg.maxNrOfSteps = 50000;
      slCfg.xtol = 1e-12;
      slCfg.initialResidualToUse = 0.;
      slCfg.simplificationMode =
        pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
      slCfg.inactiveRemovalCutoff = 0.;
      slCfg.doInnerIterations = false;
      slCfg.linkSwappingMode = pcm::LinkSwappingMode::NO_SWAPPING;
      slCfg.slipMode = pcm::SlipMode::FIXED;
      forceBalancer.run(slCfg);
      CHECK_THAT(forceBalancer.getDanglingWeightFraction(),
                 Catch::Matchers::WithinRel(expectedDanglingFraction));
      CHECK_THAT(forceBalancer.getSolubleWeightFraction(),
                 Catch::Matchers::WithinRel(expectedSolventFraction));

      pcm::MEHPSimulator forceRelaxation =
        pcm::MEHPSimulator(u2,
                           {},
                           2,
                           pcm::EntanglementMode::NONE,
                           nullptr,
                           false,
                           true,
                           false,
                           false);
      pcm::MEHPNonlinearSolverConfig nlCfg;
      nlCfg.simplificationMode =
        pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
      forceRelaxation.run(nlCfg);
      CHECK_THAT(forceRelaxation.getDanglingWeightFraction(),
                 Catch::Matchers::WithinRel(expectedDanglingFraction));
      CHECK_THAT(forceRelaxation.getSolubleWeightFraction(),
                 Catch::Matchers::WithinRel(expectedSolventFraction));
    }
  }
}

#ifdef USE_METAL
TEST_CASE("MEHPForceBalance2 Metal SLE Solvers Accuracy and Benchmark",
          "[MEHPSimulator][Metal]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "Running test \"MEHPForceBalance2 Metal SLE Solvers Accuracy "
               "and Benchmark\""
            << std::endl;

  const int N = 6000;
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(N * 7);

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  Eigen::VectorXd diag = Eigen::VectorXd::Zero(N);
  for (int i = 0; i < N; ++i) {
    int n1 = (i + 1) % N;
    int n2 = (i + 3) % N;
    int n3 = (i + 17) % N;
    triplets.push_back({ i, n1, 1.0 });
    triplets.push_back({ n1, i, 1.0 });
    triplets.push_back({ i, n2, 0.5 });
    triplets.push_back({ n2, i, 0.5 });
    triplets.push_back({ i, n3, 0.25 });
    triplets.push_back({ n3, i, 0.25 });
    diag[i] -= (1.0 + 0.5 + 0.25 + 0.05);
  }
  for (int i = 0; i < N; ++i) {
    triplets.push_back({ i, i, diag[i] });
  }

  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();

  Eigen::VectorXd x_exact(N);
  for (int i = 0; i < N; ++i)
    x_exact[i] = std::sin(0.01 * i);
  Eigen::VectorXd b = A * x_exact;

  int iters_cpu_cg = 0, iters_metal_cg = 0;
  int iters_cpu_bb = 0, iters_metal_bb = 0;

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::IdentityPreconditioner>
    cpu_cg;
  cpu_cg.setTolerance(1e-5);
  cpu_cg.setMaxIterations(200);
  cpu_cg.compute(A);
  Eigen::VectorXd x_cpu_cg = cpu_cg.solve(b);
  iters_cpu_cg = static_cast<int>(cpu_cg.iterations());

  Eigen::VectorXd x_metal_cg =
    pcm::solveSLEMetalCG(A, b, 1e-5, 200, iters_metal_cg);

  Eigen::VectorXd x_cpu_bb = Eigen::gradientDescentBarzilaiBorwein(
    A, b, 0.01, 1e-5, 500, iters_cpu_bb, true);

  Eigen::VectorXd x_metal_bb = pcm::solveSLEMetalBarzilaiBorwein(
    A, b, 0.01, 1e-5, 500, iters_metal_bb, true);

  double rel_err_metal_cg = (x_metal_cg - x_exact).norm() / x_exact.norm();
  double rel_err_metal_bb = (x_metal_bb - x_exact).norm() / x_exact.norm();
  double diff_cg = (x_metal_cg - x_cpu_cg).norm() / x_cpu_cg.norm();
  double diff_bb = (x_metal_bb - x_cpu_bb).norm() / x_cpu_bb.norm();

  std::cout << "Metal CG relative error to exact: " << rel_err_metal_cg
            << ", diff to CPU CG: " << diff_cg << std::endl;
  std::cout << "Metal BB relative error to exact: " << rel_err_metal_bb
            << ", diff to CPU BB: " << diff_bb << std::endl;
  std::cout << "Iterations: CPU CG=" << iters_cpu_cg
            << ", Metal CG=" << iters_metal_cg << ", CPU BB=" << iters_cpu_bb
            << ", Metal BB=" << iters_metal_bb << std::endl;

  CHECK(rel_err_metal_cg < 1e-3);
  CHECK(rel_err_metal_bb < 1e-3);
  CHECK(diff_cg < 1e-3);
  CHECK(diff_bb < 1e-3);

  BENCHMARK("CPU Conjugate Gradient (N=6000)")
  {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                             Eigen::Lower | Eigen::Upper,
                             Eigen::IdentityPreconditioner>
      solver;
    solver.setTolerance(1e-5);
    solver.setMaxIterations(200);
    solver.compute(A);
    return solver.solve(b);
  };

  BENCHMARK("Metal GPU Conjugate Gradient (N=6000)")
  {
    int iters = 0;
    return pcm::solveSLEMetalCG(A, b, 1e-5, 200, iters);
  };

  BENCHMARK("CPU Barzilai-Borwein (N=6000)")
  {
    int iters = 0;
    return Eigen::gradientDescentBarzilaiBorwein(
      A, b, 0.01, 1e-5, 500, iters, true);
  };

  BENCHMARK("Metal GPU Barzilai-Borwein (N=6000)")
  {
    int iters = 0;
    return pcm::solveSLEMetalBarzilaiBorwein(
      A, b, 0.01, 1e-5, 500, iters, true);
  };
}

TEST_CASE("MEHPForceBalance2 Metal vs CPU Network Relaxation Benchmark",
          "[!benchmark][MEHPSimulator][Metal]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  size_t nrOfBeads = 80;
  size_t nrOfChains = 25;
  pe::Universe universe = pe::Universe(20.0, 20.0, 20.0);
  for (size_t i = 0; i < nrOfChains; ++i) {
    long int beadOffset = static_cast<long int>(i * nrOfBeads);
    std::vector<long int> atomIds(nrOfBeads);
    std::iota(atomIds.begin(), atomIds.end(), 1 + beadOffset);
    std::vector<double> coords(nrOfBeads);
    std::iota(
      coords.begin(), coords.end(), 1. + static_cast<double>(beadOffset) * 0.5);
    universe.addAtoms(atomIds,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 1),
                      coords,
                      coords,
                      coords,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0));
    std::vector<long int> bondFrom(nrOfBeads - 1);
    std::vector<long int> bondTo(nrOfBeads - 1);
    std::iota(bondFrom.begin(), bondFrom.end(), 1 + beadOffset);
    std::iota(bondTo.begin(), bondTo.end(), 2 + beadOffset);
    universe.addBonds(bondFrom, bondTo);
  }
  std::vector<long int> crossFrom, crossTo;
  for (size_t i = 0; i < nrOfChains - 1; ++i) {
    crossFrom.push_back(static_cast<long int>(i * nrOfBeads + 10));
    crossTo.push_back(static_cast<long int>((i + 1) * nrOfBeads + 20));
    crossFrom.push_back(static_cast<long int>(i * nrOfBeads + 40));
    crossTo.push_back(
      static_cast<long int>(((i + 2) % nrOfChains) * nrOfBeads + 50));
  }
  universe.addBonds(crossFrom, crossTo);

  BENCHMARK("CPU Conjugate Gradient Network Relaxation")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-4;
    cfg.solver = pcm::SLESolver::CONJUGATE_GRADIENT;
    fb.run(cfg);
    return fb.getNrOfActiveSprings();
  };

  BENCHMARK("Metal GPU Conjugate Gradient Network Relaxation")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-4;
    cfg.solver = pcm::SLESolver::METAL_CG;
    fb.run(cfg);
    return fb.getNrOfActiveSprings();
  };

  BENCHMARK("CPU Barzilai-Borwein Network Relaxation")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-4;
    cfg.solver = pcm::SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT;
    fb.run(cfg);
    return fb.getNrOfActiveSprings();
  };

  BENCHMARK("Metal GPU Barzilai-Borwein Network Relaxation")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-4;
    cfg.solver = pcm::SLESolver::METAL_BARZILAI_BORWEIN;
    fb.run(cfg);
    return fb.getNrOfActiveSprings();
  };
}

TEST_CASE("MEHPForceBalance2 Metal Buffer Pooling and Reallocation",
          "[MEHPSimulator][Metal]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "Running test \"MEHPForceBalance2 Metal Buffer Pooling and "
               "Reallocation\""
            << std::endl;

  const int N1 = 3000;
  std::vector<Eigen::Triplet<double>> triplets1;
  triplets1.reserve(N1 * 3);
  for (int i = 0; i < N1; ++i) {
    triplets1.push_back({ i, i, 4.0 });
    triplets1.push_back({ i, (i + 1) % N1, -1.0 });
    triplets1.push_back({ (i + 1) % N1, i, -1.0 });
  }
  Eigen::SparseMatrix<double> A1(N1, N1);
  A1.setFromTriplets(triplets1.begin(), triplets1.end());
  A1.makeCompressed();
  Eigen::VectorXd b1 = Eigen::VectorXd::Ones(N1);

  int iters = 0;
  Eigen::VectorXd x1_cg = pcm::solveSLEMetalCG(A1, b1, 1e-4, 100, iters);
  CHECK(iters > 0);
  CHECK((A1 * x1_cg - b1).norm() / b1.norm() < 1e-3);

  const int N2 = 1500;
  std::vector<Eigen::Triplet<double>> triplets2;
  triplets2.reserve(N2 * 3);
  for (int i = 0; i < N2; ++i) {
    triplets2.push_back({ i, i, 4.0 });
    triplets2.push_back({ i, (i + 1) % N2, -1.0 });
    triplets2.push_back({ (i + 1) % N2, i, -1.0 });
  }
  Eigen::SparseMatrix<double> A2(N2, N2);
  A2.setFromTriplets(triplets2.begin(), triplets2.end());
  A2.makeCompressed();
  Eigen::VectorXd b2 = Eigen::VectorXd::Ones(N2);

  Eigen::VectorXd x2_cg = pcm::solveSLEMetalCG(A2, b2, 1e-4, 100, iters);
  CHECK(iters > 0);
  CHECK((A2 * x2_cg - b2).norm() / b2.norm() < 1e-3);

  const int N3 = 6000;
  std::vector<Eigen::Triplet<double>> triplets3;
  triplets3.reserve(N3 * 3);
  for (int i = 0; i < N3; ++i) {
    triplets3.push_back({ i, i, 4.0 });
    triplets3.push_back({ i, (i + 1) % N3, -1.0 });
    triplets3.push_back({ (i + 1) % N3, i, -1.0 });
  }
  Eigen::SparseMatrix<double> A3(N3, N3);
  A3.setFromTriplets(triplets3.begin(), triplets3.end());
  A3.makeCompressed();
  Eigen::VectorXd b3 = Eigen::VectorXd::Ones(N3);

  Eigen::VectorXd x3_cg = pcm::solveSLEMetalCG(A3, b3, 1e-4, 100, iters);
  CHECK(iters > 0);
  CHECK((A3 * x3_cg - b3).norm() / b3.norm() < 1e-3);

  Eigen::VectorXd x3_bb =
    pcm::solveSLEMetalBarzilaiBorwein(A3, b3, 0.1, 1e-4, 200, iters, true);
  CHECK(iters > 0);
  CHECK((A3 * x3_bb - b3).norm() / b3.norm() < 1e-3);
}

TEST_CASE("MEHPForceBalance2 Metal Center-of-Mass Zero-Mode Projection",
          "[MEHPSimulator][Metal]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "Running test \"MEHPForceBalance2 Metal Center-of-Mass "
               "Zero-Mode Projection\""
            << std::endl;

  const int S = 8;
  const int numNodes = S * S * S;
  const int N = numNodes * 3;
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(N * 7);

  auto idx = [](int x, int y, int z) {
    return (x % S) * S * S + (y % S) * S + (z % S);
  };

  for (int x = 0; x < S; ++x) {
    for (int y = 0; y < S; ++y) {
      for (int z = 0; z < S; ++z) {
        int u = idx(x, y, z);
        int neighbors[6] = { idx(x + 1, y, z), idx(x - 1 + S, y, z),
                             idx(x, y + 1, z), idx(x, y - 1 + S, z),
                             idx(x, y, z + 1), idx(x, y, z - 1 + S) };
        for (int d = 0; d < 3; ++d) {
          triplets.push_back({ 3 * u + d, 3 * u + d, 6.0 });
          for (int v : neighbors) {
            triplets.push_back({ 3 * u + d, 3 * v + d, -1.0 });
          }
        }
      }
    }
  }

  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();

  Eigen::VectorXd ones = Eigen::VectorXd::Zero(N);
  for (int i = 0; i < numNodes; ++i)
    ones[3 * i] = 1.0;
  CHECK((A * ones).norm() < 1e-12);

  Eigen::VectorXd x_exact = Eigen::VectorXd::Zero(N);
  for (int x = 0; x < S; ++x) {
    for (int y = 0; y < S; ++y) {
      for (int z = 0; z < S; ++z) {
        int u = idx(x, y, z);
        x_exact[3 * u + 0] = std::sin(2.0 * M_PI * x / S);
        x_exact[3 * u + 1] = std::sin(2.0 * M_PI * y / S);
        x_exact[3 * u + 2] = std::sin(2.0 * M_PI * z / S);
      }
    }
  }
  Eigen::VectorXd b = A * x_exact;

  Eigen::VectorXd initialX = Eigen::VectorXd::Zero(N);
  for (int i = 0; i < numNodes; ++i) {
    initialX[3 * i + 0] = 15.0;
    initialX[3 * i + 1] = -25.0;
    initialX[3 * i + 2] = 40.0;
  }

  int itersWithout = 0;
  Eigen::VectorXd x_no_proj = pcm::solveSLEMetalCG(
    A, b, 1e-4, 150, itersWithout, initialX, /*projectZeroModes=*/false);

  int itersWith = 0;
  Eigen::VectorXd x_with_proj = pcm::solveSLEMetalCG(
    A, b, 1e-4, 150, itersWith, initialX, /*projectZeroModes=*/true);

  double meanX_no = 0, meanY_no = 0, meanZ_no = 0;
  for (int i = 0; i < numNodes; ++i) {
    meanX_no += x_no_proj[3 * i + 0];
    meanY_no += x_no_proj[3 * i + 1];
    meanZ_no += x_no_proj[3 * i + 2];
  }
  meanX_no /= numNodes;
  meanY_no /= numNodes;
  meanZ_no /= numNodes;

  double meanX_with = 0, meanY_with = 0, meanZ_with = 0;
  for (int i = 0; i < numNodes; ++i) {
    meanX_with += x_with_proj[3 * i + 0];
    meanY_with += x_with_proj[3 * i + 1];
    meanZ_with += x_with_proj[3 * i + 2];
  }
  meanX_with /= numNodes;
  meanY_with /= numNodes;
  meanZ_with /= numNodes;

  std::cout << "COM without projection: (" << meanX_no << ", " << meanY_no
            << ", " << meanZ_no << ")" << std::endl;
  std::cout << "COM with projection: (" << meanX_with << ", " << meanY_with
            << ", " << meanZ_with << ")" << std::endl;

  CHECK(std::abs(meanX_no) > 1.0);
  CHECK(std::abs(meanY_no) > 1.0);
  CHECK(std::abs(meanZ_no) > 1.0);

  CHECK(std::abs(meanX_with) < 1e-4);
  CHECK(std::abs(meanY_with) < 1e-4);
  CHECK(std::abs(meanZ_with) < 1e-4);

  CHECK((A * x_with_proj - b).norm() / b.norm() < 1e-3);

  int itersBB = 0;
  Eigen::VectorXd x_bb_proj = pcm::solveSLEMetalBarzilaiBorwein(
    A, b, 0.05, 1e-4, 250, itersBB, true, initialX, /*projectZeroModes=*/true);
  double meanX_bb = 0, meanY_bb = 0, meanZ_bb = 0;
  for (int i = 0; i < numNodes; ++i) {
    meanX_bb += x_bb_proj[3 * i + 0];
    meanY_bb += x_bb_proj[3 * i + 1];
    meanZ_bb += x_bb_proj[3 * i + 2];
  }
  meanX_bb /= numNodes;
  meanY_bb /= numNodes;
  meanZ_bb /= numNodes;

  CHECK(std::abs(meanX_bb) < 1e-4);
  CHECK(std::abs(meanY_bb) < 1e-4);
  CHECK(std::abs(meanZ_bb) < 1e-4);
  CHECK((A * x_bb_proj - b).norm() / b.norm() < 1e-3);
}

TEST_CASE("MEHPForceBalance2 Metal SIMDgroup-Cooperative SpMV Degree Variation",
          "[MEHPSimulator][Metal]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "Running test \"MEHPForceBalance2 Metal SIMDgroup-Cooperative "
               "SpMV Degree Variation\""
            << std::endl;

  const int N = 9000;
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(N * 8);

  Eigen::VectorXd diag = Eigen::VectorXd::Zero(N);
  for (int i = 0; i < N; ++i) {
    int deg = (i % 11 == 0) ? 1 : ((i % 5 == 0) ? 6 : ((i % 2 == 0) ? 4 : 2));
    for (int d = 1; d <= deg; ++d) {
      int j = (i + d * 3) % N;
      double w = 0.5 / deg;
      triplets.push_back({ i, j, w });
      triplets.push_back({ j, i, w });
      diag[i] -= (w + 0.05);
      diag[j] -= (w + 0.05);
    }
  }
  for (int i = 0; i < N; ++i) {
    triplets.push_back({ i, i, diag[i] });
  }

  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();

  Eigen::VectorXd x_exact = Eigen::VectorXd::Random(N);
  Eigen::VectorXd b = A * x_exact;

  int iters_cg = 0, iters_bb = 0;
  Eigen::VectorXd x_cg = pcm::solveSLEMetalCG(A, b, 1e-4, 200, iters_cg);
  Eigen::VectorXd x_bb =
    pcm::solveSLEMetalBarzilaiBorwein(A, b, 0.05, 1e-4, 300, iters_bb, true);

  double rel_err_cg = (x_cg - x_exact).norm() / x_exact.norm();
  double rel_err_bb = (x_bb - x_exact).norm() / x_exact.norm();

  std::cout << "Degree variation CG rel err: " << rel_err_cg
            << ", iters: " << iters_cg << std::endl;
  std::cout << "Degree variation BB rel err: " << rel_err_bb
            << ", iters: " << iters_bb << std::endl;

  CHECK(rel_err_cg < 1e-3);
  CHECK(rel_err_bb < 1e-3);
}

TEST_CASE(
  "MEHPForceBalance2 Metal Multi-RHS 3D Block-SpMV Accuracy and Benchmark",
  "[MEHPSimulator][Metal]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "Running test \"MEHPForceBalance2 Metal Multi-RHS 3D Block-SpMV "
               "Accuracy and Benchmark\""
            << std::endl;

  const int numNodes = 2000;
  const int N = numNodes * 3;

  std::vector<Eigen::Triplet<double>> scalarTriplets;
  scalarTriplets.reserve(numNodes * 8);
  Eigen::VectorXd scalarDiag = Eigen::VectorXd::Zero(numNodes);

  for (int i = 0; i < numNodes; ++i) {
    int n1 = (i + 1) % numNodes;
    int n2 = (i + 3) % numNodes;
    int n3 = (i + 17) % numNodes;
    scalarTriplets.push_back({ i, n1, 1.0 });
    scalarTriplets.push_back({ n1, i, 1.0 });
    scalarTriplets.push_back({ i, n2, 0.5 });
    scalarTriplets.push_back({ n2, i, 0.5 });
    scalarTriplets.push_back({ i, n3, 0.25 });
    scalarTriplets.push_back({ n3, i, 0.25 });
    scalarDiag[i] -= (1.0 + 0.5 + 0.25);
    scalarDiag[n1] -= 1.0;
    scalarDiag[n2] -= 0.5;
    scalarDiag[n3] -= 0.25;
  }
  for (int i = 0; i < numNodes; ++i) {
    scalarDiag[i] -= 0.5;
    scalarTriplets.push_back({ i, i, scalarDiag[i] });
  }

  Eigen::SparseMatrix<double> scalarL(numNodes, numNodes);
  scalarL.setFromTriplets(scalarTriplets.begin(), scalarTriplets.end());
  scalarL.makeCompressed();

  std::vector<Eigen::Triplet<double>> fullTriplets;
  fullTriplets.reserve(numNodes * 8 * 3);
  for (int k = 0; k < scalarL.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(scalarL, k); it; ++it) {
      int r = static_cast<int>(it.row());
      int c = static_cast<int>(it.col());
      double v = it.value();
      for (int d = 0; d < 3; ++d) {
        fullTriplets.push_back({ 3 * r + d, 3 * c + d, v });
      }
    }
  }
  Eigen::SparseMatrix<double> fullA(N, N);
  fullA.setFromTriplets(fullTriplets.begin(), fullTriplets.end());
  fullA.makeCompressed();

  CHECK(scalarL.nonZeros() * 3 == fullA.nonZeros());
  std::cout << "Scalar Laplacian NNZ: " << scalarL.nonZeros()
            << ", Full 3Nx3N NNZ: " << fullA.nonZeros()
            << " (exact 3x reduction in matrix nonzeros and storage)"
            << std::endl;

  Eigen::VectorXd x_exact(N);
  for (int i = 0; i < N; ++i) {
    x_exact[i] = std::sin(0.01 * i) + 0.2 * std::cos(0.05 * i);
  }

  Eigen::VectorXd y_metal_3d = pcm::spmvMetal3D(scalarL, x_exact);
  Eigen::VectorXd y_metal_full = pcm::spmvMetal(fullA, x_exact);
  Eigen::VectorXd y_cpu = fullA * x_exact;

  double spmv_err_3d = (y_metal_3d - y_cpu).norm() / y_cpu.norm();
  double spmv_err_full = (y_metal_full - y_cpu).norm() / y_cpu.norm();
  double spmv_diff_3d_full = (y_metal_3d - y_metal_full).norm() / y_cpu.norm();

  std::cout << "Pure SpMV relative error to CPU: 3D Block-SpMV = "
            << spmv_err_3d << ", Full SpMV = " << spmv_err_full
            << ", 3D vs Full diff = " << spmv_diff_3d_full << std::endl;
  CHECK(spmv_err_3d < 1e-5);
  CHECK(spmv_err_full < 1e-5);
  CHECK(spmv_diff_3d_full < 1e-5);

  Eigen::VectorXd b = y_cpu;

  int iters_3d = 0;
  int iters_full = 0;

  Eigen::VectorXd x_3d = pcm::solveSLEMetalCG3D(
    scalarL, b, 1e-5, 300, iters_3d, Eigen::VectorXd(), false, false);

  Eigen::VectorXd x_full = pcm::solveSLEMetalCG(
    fullA, b, 1e-5, 300, iters_full, Eigen::VectorXd(), false, false);

  double rel_err_3d = (x_3d - x_exact).norm() / x_exact.norm();
  double rel_err_full = (x_full - x_exact).norm() / x_exact.norm();
  double diff_3d_full = (x_3d - x_full).norm() / x_full.norm();

  std::cout << "3D Block-SpMV rel err: " << rel_err_3d
            << ", iters: " << iters_3d << std::endl;
  std::cout << "Full 3Nx3N SpMV rel err: " << rel_err_full
            << ", iters: " << iters_full << std::endl;
  std::cout << "Diff between 3D Block and Full solves: " << diff_3d_full
            << std::endl;

  CHECK(rel_err_3d < 1e-3);
  CHECK(rel_err_full < 1e-3);
  CHECK(diff_3d_full < 1e-4);

  int iters_auto = 0;
  Eigen::VectorXd x_auto = pcm::solveSLEMetalCG(
    scalarL, b, 1e-5, 300, iters_auto, Eigen::VectorXd(), false, false);
  CHECK((x_auto - x_3d).norm() / x_3d.norm() < 1e-5);
  CHECK(iters_auto == iters_3d);

  int iters_bb_3d = 0;
  Eigen::VectorXd x_bb_3d = pcm::solveSLEMetalBarzilaiBorwein3D(
    scalarL, b, 0.01, 1e-5, 500, iters_bb_3d, true);
  double rel_err_bb_3d = (x_bb_3d - x_exact).norm() / x_exact.norm();
  std::cout << "3D BB rel err: " << rel_err_bb_3d << ", iters: " << iters_bb_3d
            << std::endl;
  CHECK(rel_err_bb_3d < 1e-3);

  const int benchNodes = 10000;
  std::vector<Eigen::Triplet<double>> benchScalarTriplets;
  benchScalarTriplets.reserve(benchNodes * 7);
  for (int i = 0; i < benchNodes; ++i) {
    int n1 = (i + 1) % benchNodes;
    int n2 = (i + 5) % benchNodes;
    int n3 = (i + 23) % benchNodes;
    benchScalarTriplets.push_back({ i, n1, 1.0 });
    benchScalarTriplets.push_back({ i, n2, 0.5 });
    benchScalarTriplets.push_back({ i, n3, 0.25 });
    benchScalarTriplets.push_back({ i, i, -2.5 });
  }
  Eigen::SparseMatrix<double> benchScalarL(benchNodes, benchNodes);
  benchScalarL.setFromTriplets(benchScalarTriplets.begin(),
                               benchScalarTriplets.end());
  benchScalarL.makeCompressed();

  std::vector<Eigen::Triplet<double>> benchFullTriplets;
  benchFullTriplets.reserve(benchNodes * 7 * 3);
  for (int k = 0; k < benchScalarL.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(benchScalarL, k); it;
         ++it) {
      int r = static_cast<int>(it.row());
      int c = static_cast<int>(it.col());
      double v = it.value();
      for (int d = 0; d < 3; ++d) {
        benchFullTriplets.push_back({ 3 * r + d, 3 * c + d, v });
      }
    }
  }
  Eigen::SparseMatrix<double> benchFullA(benchNodes * 3, benchNodes * 3);
  benchFullA.setFromTriplets(benchFullTriplets.begin(),
                             benchFullTriplets.end());
  benchFullA.makeCompressed();

  pcm::SpmvBenchmarkResult benchResult =
    pcm::benchmarkSpMVMetal(benchScalarL, benchFullA, 100);

  std::cout << "Pure SpMV Benchmark (100 repetitions, N=" << benchNodes * 3
            << "): "
            << "Full SpMV time = " << benchResult.timeFull_ms << " ms, "
            << "3D Block-SpMV time = " << benchResult.time3D_ms << " ms, "
            << "Speedup = " << benchResult.speedup << "x, "
            << "Bandwidth reduction = " << benchResult.bandwidthReductionRatio
            << "x "
            << "(" << benchResult.memory3D_bytes / 1024 << " KB vs "
            << benchResult.memoryFull_bytes / 1024 << " KB)" << std::endl;

  CHECK(benchResult.bandwidthReductionRatio >= 2.9);
  CHECK(benchResult.speedup >= 1.35);
}

TEST_CASE(
  "MEHPForceBalance2 Metal GPU-Resident Jacobi Preconditioner Convergence",
  "[MEHPSimulator][Metal]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "Running test \"MEHPForceBalance2 Metal GPU-Resident Jacobi "
               "Preconditioner Convergence\""
            << std::endl;

  const int numNodes = 3000;
  const int N = numNodes * 3;

  std::vector<Eigen::Triplet<double>> scalarTriplets;
  scalarTriplets.reserve(numNodes * 8);
  Eigen::VectorXd scalarDiag = Eigen::VectorXd::Zero(numNodes);

  for (int i = 0; i < numNodes; ++i) {
    int n1 = (i + 1) % numNodes;
    int n2 = (i + 7) % numNodes;
    double k1 = 0.2 + 3.0 * std::sin(0.1 * i) * std::sin(0.1 * i);
    double k2 = 0.1 + 2.0 * std::cos(0.05 * i) * std::cos(0.05 * i);
    scalarTriplets.push_back({ i, n1, k1 });
    scalarTriplets.push_back({ n1, i, k1 });
    scalarTriplets.push_back({ i, n2, k2 });
    scalarTriplets.push_back({ n2, i, k2 });
    scalarDiag[i] -= (k1 + k2);
    scalarDiag[n1] -= k1;
    scalarDiag[n2] -= k2;

    if (i % 12 == 0) {
      int x1 = (i + 31) % numNodes;
      int x2 = (i + 97) % numNodes;
      double kX = 4.0;
      scalarTriplets.push_back({ i, x1, kX });
      scalarTriplets.push_back({ x1, i, kX });
      scalarTriplets.push_back({ i, x2, kX });
      scalarTriplets.push_back({ x2, i, kX });
      scalarDiag[i] -= 2.0 * kX;
      scalarDiag[x1] -= kX;
      scalarDiag[x2] -= kX;
    }
  }

  for (int i = 0; i < numNodes; ++i) {
    scalarDiag[i] -= 0.05;
    scalarTriplets.push_back({ i, i, scalarDiag[i] });
  }

  Eigen::SparseMatrix<double> scalarL(numNodes, numNodes);
  scalarL.setFromTriplets(scalarTriplets.begin(), scalarTriplets.end());
  scalarL.makeCompressed();

  double minDiag = std::abs(scalarDiag.maxCoeff());
  double maxDiag = std::abs(scalarDiag.minCoeff());
  double condRatio = maxDiag / minDiag;
  std::cout << "Polymer network diagonal range: [" << minDiag << ", " << maxDiag
            << "], diagonal variation ratio = " << condRatio << std::endl;
  CHECK(condRatio > 4.0);

  Eigen::VectorXd x_exact(N);
  for (int i = 0; i < N; ++i) {
    x_exact[i] = std::sin(0.02 * i) + 0.5 * std::cos(0.01 * i);
  }

  Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
  for (int k = 0; k < scalarL.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(scalarL, k); it; ++it) {
      int r = static_cast<int>(it.row());
      int c = static_cast<int>(it.col());
      double v = it.value();
      b.segment(3 * r, 3) += v * x_exact.segment(3 * c, 3);
    }
  }

  int iters_unpreconditioned = 0;
  int iters_jacobi = 0;

  Eigen::VectorXd x_unprec = pcm::solveSLEMetalCG3D(scalarL,
                                                    b,
                                                    1e-5,
                                                    400,
                                                    iters_unpreconditioned,
                                                    Eigen::VectorXd(),
                                                    false,
                                                    false);

  Eigen::VectorXd x_jacobi = pcm::solveSLEMetalCG3D(
    scalarL, b, 1e-5, 400, iters_jacobi, Eigen::VectorXd(), false, true);

  double rel_err_unprec = (x_unprec - x_exact).norm() / x_exact.norm();
  double rel_err_jacobi = (x_jacobi - x_exact).norm() / x_exact.norm();

  double iterReduction =
    100.0 * (1.0 - static_cast<double>(iters_jacobi) /
                     static_cast<double>(iters_unpreconditioned));

  std::cout << "Unpreconditioned CG iters: " << iters_unpreconditioned
            << ", rel err: " << rel_err_unprec << std::endl;
  std::cout << "Jacobi Preconditioned CG iters: " << iters_jacobi
            << ", rel err: " << rel_err_jacobi << std::endl;
  std::cout << "Iteration reduction from GPU-Resident Jacobi: " << iterReduction
            << "%" << std::endl;

  CHECK(rel_err_jacobi < 1e-3);
  CHECK(rel_err_unprec < 1e-3);
  CHECK(iterReduction >= 18.0);
  CHECK(iters_jacobi < iters_unpreconditioned);
}

TEST_CASE("MEHPForceBalance2 Metal End-to-End GPU Relaxation",
          "[MEHPSimulator][Metal]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "Running test \"MEHPForceBalance2 Metal End-to-End GPU "
               "Relaxation\""
            << std::endl;

  {
    const size_t numNodes = 4;
    const size_t numSprings = 4;

    Eigen::VectorXd coords(12);
    coords << 0, 0, 0, 1, 0, 0, 4, 0, 0, 0.05, 0, 0;

    Eigen::SparseMatrix<double> dummyA(numNodes, numNodes);
    dummyA.insert(0, 0) = 1.0;
    dummyA.insert(1, 1) = 1.0;
    dummyA.insert(2, 2) = 1.0;
    dummyA.insert(3, 3) = 1.0;
    dummyA.makeCompressed();
    Eigen::VectorXd dummyB = coords;
    int dummyIters = 0;
    pcm::solveSLEMetalCG3D(
      dummyA, dummyB, 1e-5, 1, dummyIters, coords, false, false);

    std::vector<int> springA = { 0, 0, 0, 1 };
    std::vector<int> springB = { 1, 2, 3, 2 };
    std::vector<double> boxOffset = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector<double> contourLength = { 1.0, 1.0, 1.0, 1.0 };

    pcm::SpringRelaxationData springData;
    springData.springIndexA = springA.data();
    springData.springIndexB = springB.data();
    springData.springBoxOffset = boxOffset.data();
    springData.springContourLength = contourLength.data();
    springData.numSprings = numSprings;
    springData.breakingThreshold = 2.5;
    springData.inactiveCutoff = 0.1;

    pcm::SpringEvaluationResult evalRes =
      pcm::evaluateSpringStrainAndInactivityMetal(numNodes, springData);

    std::cout << "GPU Spring Evaluation: numBroken = " << evalRes.numBroken
              << ", numInactive = " << evalRes.numInactive
              << ", totalFlagged = " << evalRes.totalFlagged << std::endl;

    CHECK(evalRes.numBroken >= 2);
    CHECK(evalRes.numInactive >= 1);
    CHECK(evalRes.totalFlagged >= 3);
    CHECK((evalRes.springFlags[0] & 1) == 0);
    CHECK((evalRes.springFlags[0] & 2) == 0);
    CHECK((evalRes.springFlags[1] & 1) != 0);
    CHECK((evalRes.springFlags[2] & 2) != 0);
  }

  {
    size_t nrOfBeads = 60;
    size_t nrOfChains = 20;
    pe::Universe universe = pe::Universe(20.0, 20.0, 20.0);
    for (size_t i = 0; i < nrOfChains; ++i) {
      long int beadOffset = static_cast<long int>(i * nrOfBeads);
      std::vector<long int> atomIds(nrOfBeads);
      std::iota(atomIds.begin(), atomIds.end(), 1 + beadOffset);
      std::vector<double> coords(nrOfBeads);
      std::iota(coords.begin(),
                coords.end(),
                1. + static_cast<double>(beadOffset) * 0.5);
      universe.addAtoms(
        atomIds,
        pylimer_tools::utils::initializeWithValue(nrOfBeads, 1),
        coords,
        coords,
        coords,
        pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
        pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
        pylimer_tools::utils::initializeWithValue(nrOfBeads, 0));
      std::vector<long int> bondFrom(nrOfBeads - 1);
      std::vector<long int> bondTo(nrOfBeads - 1);
      std::iota(bondFrom.begin(), bondFrom.end(), 1 + beadOffset);
      std::iota(bondTo.begin(), bondTo.end(), 2 + beadOffset);
      universe.addBonds(bondFrom, bondTo);
    }
    std::vector<long int> crossFrom, crossTo;
    for (size_t i = 0; i < nrOfChains - 1; ++i) {
      crossFrom.push_back(static_cast<long int>(i * nrOfBeads + 10));
      crossTo.push_back(static_cast<long int>((i + 1) * nrOfBeads + 20));
      crossFrom.push_back(static_cast<long int>(i * nrOfBeads + 40));
      crossTo.push_back(
        static_cast<long int>(((i + 2) % nrOfChains) * nrOfBeads + 50));
    }
    universe.addBonds(crossFrom, crossTo);

    pcm::MEHPSimulator fb_cpu(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfgCpu;
    cfgCpu.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfgCpu.inactiveRemovalCutoff = 1e-4;
    cfgCpu.solver = pcm::SLESolver::CONJUGATE_GRADIENT;
    fb_cpu.run(cfgCpu);
    double res_cpu = fb_cpu.getResidualNorm();
    size_t active_cpu = fb_cpu.getNrOfActiveSprings();

    pcm::MEHPSimulator fb_metal(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfgMetal = cfgCpu;
    cfgMetal.solver = pcm::SLESolver::METAL_CG;
    fb_metal.run(cfgMetal);
    double res_metal = fb_metal.getResidualNorm();
    size_t active_metal = fb_metal.getNrOfActiveSprings();

    std::cout << "CPU CG Final Residual: " << res_cpu
              << ", Active Springs: " << active_cpu << std::endl;
    std::cout << "Metal CG Final Residual: " << res_metal
              << ", Active Springs: " << active_metal << std::endl;

    CHECK(active_metal > 0);
    CHECK(active_metal == active_cpu);
    CHECK(res_metal < 1e-3);
    CHECK_NOTHROW(fb_metal.validateNetwork());

    pcm::MEHPSimulator fb_bb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfgBb = cfgCpu;
    cfgBb.solver = pcm::SLESolver::METAL_BARZILAI_BORWEIN;
    fb_bb.run(cfgBb);
    double res_bb = fb_bb.getResidualNorm();
    size_t active_bb = fb_bb.getNrOfActiveSprings();

    std::cout << "Metal BB Final Residual: " << res_bb
              << ", Active Springs: " << active_bb << std::endl;

    CHECK(active_bb == active_cpu);
    CHECK(res_bb < 1e-3);
    CHECK_NOTHROW(fb_bb.validateNetwork());
  }

  SECTION("Metal SIMD Cooperative Width Benchmark on Polymer Networks")
  {
    std::cout << "\n=== Benchmarking Metal SIMD Cooperative Widths (coop2 vs "
                 "coop4 vs coop8) ==="
              << std::endl;
    const size_t numNodes = 10000;
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(numNodes * 3);
    for (size_t i = 0; i < numNodes - 1; ++i) {
      trips.emplace_back(i, i + 1, -1.0);
      trips.emplace_back(i + 1, i, -1.0);
    }
    for (size_t i = 0; i < numNodes / 20; ++i) {
      size_t a = (i * 17) % numNodes;
      size_t b = (i * 37 + 100) % numNodes;
      if (a != b) {
        trips.emplace_back(a, b, -1.0);
        trips.emplace_back(b, a, -1.0);
      }
    }
    Eigen::VectorXd diag =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numNodes));
    for (const auto& tr : trips) {
      diag[tr.row()] -= tr.value();
    }
    for (size_t i = 0; i < numNodes; ++i) {
      trips.emplace_back(i, i, diag[static_cast<Eigen::Index>(i)] + 0.1);
    }
    Eigen::SparseMatrix<double> scalarL(static_cast<Eigen::Index>(numNodes),
                                        static_cast<Eigen::Index>(numNodes));
    scalarL.setFromTriplets(trips.begin(), trips.end());
    scalarL.makeCompressed();

    auto bench = pcm::benchmarkCoopWidthsMetal(scalarL, 50);
    std::cout << "Coop2 time: " << bench.timeCoop2_ms << " ms" << std::endl;
    std::cout << "Coop4 time: " << bench.timeCoop4_ms << " ms" << std::endl;
    std::cout << "Coop8 time: " << bench.timeCoop8_ms << " ms" << std::endl;
    std::cout << "Speedup Coop2 vs Coop8: " << bench.speedupCoop2Vs8 << "x"
              << std::endl;
    std::cout << "Speedup Coop4 vs Coop8: " << bench.speedupCoop4Vs8 << "x"
              << std::endl;

    CHECK(bench.timeCoop2_ms > 0.0);
    CHECK(bench.timeCoop8_ms > 0.0);
  }

  SECTION("CPU Direct Solvers: N x N scalar + 3 RHS vs 3N x 3N benchmark")
  {
    std::cout << "\n=== Benchmarking CPU Direct Solvers (Scalar N x N + 3 RHS "
                 "vs 3N x 3N) ==="
              << std::endl;
    const size_t N = 3000;
    std::vector<Eigen::Triplet<double>> tripsScalar;
    tripsScalar.reserve(N * 3);
    for (size_t i = 0; i < N - 1; ++i) {
      tripsScalar.emplace_back(i, i + 1, -1.0);
      tripsScalar.emplace_back(i + 1, i, -1.0);
    }
    for (size_t i = 0; i < N / 20; ++i) {
      size_t a = (i * 17) % N;
      size_t b = (i * 37 + 50) % N;
      if (a != b) {
        tripsScalar.emplace_back(a, b, -1.0);
        tripsScalar.emplace_back(b, a, -1.0);
      }
    }
    Eigen::VectorXd diag = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(N));
    for (const auto& tr : tripsScalar)
      diag[tr.row()] -= tr.value();
    for (size_t i = 0; i < N; ++i)
      tripsScalar.emplace_back(i, i, diag[static_cast<Eigen::Index>(i)] + 1.0);

    Eigen::SparseMatrix<double> L(static_cast<Eigen::Index>(N),
                                  static_cast<Eigen::Index>(N));
    L.setFromTriplets(tripsScalar.begin(), tripsScalar.end());
    L.makeCompressed();

    std::vector<Eigen::Triplet<double>> trips3N;
    trips3N.reserve(tripsScalar.size() * 3);
    for (const auto& tr : tripsScalar) {
      for (int d = 0; d < 3; ++d) {
        trips3N.emplace_back(tr.row() * 3 + d, tr.col() * 3 + d, tr.value());
      }
    }
    Eigen::SparseMatrix<double> K3N(3 * static_cast<Eigen::Index>(N),
                                    3 * static_cast<Eigen::Index>(N));
    K3N.setFromTriplets(trips3N.begin(), trips3N.end());
    K3N.makeCompressed();

    Eigen::VectorXd b3N =
      Eigen::VectorXd::Random(3 * static_cast<Eigen::Index>(N));
    Eigen::MatrixX3d B(static_cast<Eigen::Index>(N), 3);
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(N); ++i) {
      B.row(i) = b3N.segment<3>(3 * i);
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver3N;
    solver3N.compute(K3N);
    Eigen::VectorXd sol3N = solver3N.solve(b3N);
    auto t1 = std::chrono::high_resolution_clock::now();
    double time3N_ms =
      std::chrono::duration<double, std::milli>(t1 - t0).count();

    auto t2 = std::chrono::high_resolution_clock::now();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverScalar;
    solverScalar.compute(L);
    Eigen::MatrixX3d solScalar = solverScalar.solve(B);
    auto t3 = std::chrono::high_resolution_clock::now();
    double timeScalar_ms =
      std::chrono::duration<double, std::milli>(t3 - t2).count();

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(N); ++i) {
      CHECK_THAT(sol3N[3 * i],
                 Catch::Matchers::WithinRel(solScalar(i, 0), 1e-10));
      CHECK_THAT(sol3N[3 * i + 1],
                 Catch::Matchers::WithinRel(solScalar(i, 1), 1e-10));
      CHECK_THAT(sol3N[3 * i + 2],
                 Catch::Matchers::WithinRel(solScalar(i, 2), 1e-10));
    }

    double speedup = time3N_ms / std::max(timeScalar_ms, 1e-6);
    std::cout << "3N x 3N LDLT solve time: " << time3N_ms << " ms" << std::endl;
    std::cout << "Scalar N x N + 3 RHS LDLT solve time: " << timeScalar_ms
              << " ms" << std::endl;
    std::cout << "Speedup of Scalar + 3 RHS vs 3N: " << speedup << "x"
              << std::endl;
    CHECK(speedup > 1.5);
  }
}

TEST_CASE("MEHPForceBalance2 Metal Single-Precision Solvers & Robust CPU "
          "Fallback on Ill-Conditioned Networks",
          "[MEHPSimulator][Metal][IllConditioned][Fallback]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "\n=== Testing Single-Precision GPU Solvers & CPU Fallback on "
               "Ill-Conditioned Networks ==="
            << std::endl;
  const size_t numNodes = 200;
  std::vector<Eigen::Triplet<double>> trips;
  trips.reserve(numNodes * 3);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    double stiffness = (i % 10 == 0) ? 1e4 : 1.0;
    trips.emplace_back(i, i + 1, -stiffness);
    trips.emplace_back(i + 1, i, -stiffness);
  }
  trips.emplace_back(numNodes - 2, numNodes - 1, -1e-3);
  trips.emplace_back(numNodes - 1, numNodes - 2, -1e-3);

  Eigen::VectorXd diag =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numNodes));
  for (const auto& tr : trips) {
    diag[tr.row()] -= tr.value();
  }
  for (size_t i = 0; i < numNodes; ++i) {
    trips.emplace_back(i, i, diag[static_cast<Eigen::Index>(i)] + 0.01);
  }
  Eigen::SparseMatrix<double> scalarL(static_cast<Eigen::Index>(numNodes),
                                      static_cast<Eigen::Index>(numNodes));
  scalarL.setFromTriplets(trips.begin(), trips.end());
  scalarL.makeCompressed();

  Eigen::VectorXd b =
    Eigen::VectorXd::Random(static_cast<Eigen::Index>(numNodes * 3));
  double bNorm = b.norm();
  const double targetTol = 1e-8;

  SECTION(
    "Float32 GPU without fallback on ill-conditioned network hits float floor")
  {
    int itersNoFallback = 0;
    Eigen::VectorXd xNoFallback = pcm::solveSLEMetalCG3D(scalarL,
                                                         b,
                                                         targetTol,
                                                         1000,
                                                         itersNoFallback,
                                                         Eigen::VectorXd(),
                                                         false,
                                                         true,
                                                         -1,
                                                         false);
    Eigen::VectorXd AxNo = pcm::spmvMetal3D(scalarL, xNoFallback);
    double relResNo = (b - AxNo).norm() / bNorm;
    std::cout << "GPU Float32 (no fallback) Rel Residual: " << relResNo
              << ", Iterations: " << itersNoFallback << std::endl;
    CHECK(relResNo > targetTol);
  }

  SECTION("Robust CPU Double-Precision Fallback achieves tight tolerance")
  {
    int itersFallback = 0;
    Eigen::VectorXd xFallback = pcm::solveSLEMetalCG3D(scalarL,
                                                       b,
                                                       targetTol,
                                                       1000,
                                                       itersFallback,
                                                       Eigen::VectorXd(),
                                                       false,
                                                       true,
                                                       -1,
                                                       true);

    Eigen::MatrixXd xMat(static_cast<Eigen::Index>(numNodes), 3),
      bMat(static_cast<Eigen::Index>(numNodes), 3);
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numNodes); ++i) {
      xMat(i, 0) = xFallback[3 * i + 0];
      xMat(i, 1) = xFallback[3 * i + 1];
      xMat(i, 2) = xFallback[3 * i + 2];
      bMat(i, 0) = b[3 * i + 0];
      bMat(i, 1) = b[3 * i + 1];
      bMat(i, 2) = b[3 * i + 2];
    }
    double relResFall = (bMat - scalarL * xMat).norm() / bMat.norm();
    std::cout
      << "GPU with CPU Double Fallback Rel Residual (double-precision): "
      << relResFall << ", Total Iterations: " << itersFallback << std::endl;
    CHECK(relResFall <= 1e-7);
  }
}

TEST_CASE(
  "MEHPForceBalance2 Metal Pinned Reference Node vs Center-of-Mass Projection",
  "[MEHPSimulator][Metal][PinnedNode][EliminationTree]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "\n=== Testing Elimination Tree / Pinned Reference Node vs "
               "Center-of-Mass Projection ==="
            << std::endl;
  const size_t numNodes = 400;
  std::vector<Eigen::Triplet<double>> trips;
  trips.reserve(numNodes * 4);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    trips.emplace_back(i, i + 1, -1.0);
    trips.emplace_back(i + 1, i, -1.0);
  }
  trips.emplace_back(numNodes - 1, 0, -1.0);
  trips.emplace_back(0, numNodes - 1, -1.0);
  for (size_t i = 0; i < numNodes / 8; ++i) {
    size_t u = (i * 7) % numNodes;
    size_t v = (i * 13 + 50) % numNodes;
    if (u != v) {
      trips.emplace_back(u, v, -0.5);
      trips.emplace_back(v, u, -0.5);
    }
  }

  Eigen::VectorXd diag =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numNodes));
  for (const auto& tr : trips) {
    diag[tr.row()] -= tr.value();
  }
  for (size_t i = 0; i < numNodes; ++i) {
    trips.emplace_back(i, i, diag[static_cast<Eigen::Index>(i)]);
  }
  Eigen::SparseMatrix<double> L_periodic(static_cast<Eigen::Index>(numNodes),
                                         static_cast<Eigen::Index>(numNodes));
  L_periodic.setFromTriplets(trips.begin(), trips.end());
  L_periodic.makeCompressed();

  Eigen::VectorXd b =
    Eigen::VectorXd::Random(static_cast<Eigen::Index>(numNodes * 3));
  for (int d = 0; d < 3; ++d) {
    double mean = 0.0;
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numNodes); ++i)
      mean += b[3 * i + d];
    mean /= static_cast<double>(numNodes);
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numNodes); ++i)
      b[3 * i + d] -= mean;
  }

  auto t0 = std::chrono::high_resolution_clock::now();
  int itersCom = 0;
  Eigen::VectorXd xCom = pcm::solveSLEMetalCG3D(
    L_periodic, b, 1e-4, 500, itersCom, Eigen::VectorXd(), true, true, -1);
  auto t1 = std::chrono::high_resolution_clock::now();
  double timeCom_ms =
    std::chrono::duration<double, std::milli>(t1 - t0).count();

  auto t2 = std::chrono::high_resolution_clock::now();
  int itersPinned = 0;
  Eigen::VectorXd xPinned = pcm::solveSLEMetalCG3DPinned(
    L_periodic, b, 1e-4, 500, itersPinned, Eigen::VectorXd(), 0, true);
  auto t3 = std::chrono::high_resolution_clock::now();
  double timePinned_ms =
    std::chrono::duration<double, std::milli>(t3 - t2).count();

  std::cout << "COM Projection Method:   Iterations = " << itersCom
            << ", Time = " << timeCom_ms << " ms\n"
            << "Pinned Reference Node:   Iterations = " << itersPinned
            << ", Time = " << timePinned_ms << " ms\n"
            << "Iteration Reduction:     "
            << (100.0 * (itersCom - itersPinned) / std::max(itersCom, 1))
            << "%\n"
            << "==============================================================="
               "==========\n";

  Eigen::MatrixXd xComMat(static_cast<Eigen::Index>(numNodes), 3),
    xPinnedMat(static_cast<Eigen::Index>(numNodes), 3),
    bMat(static_cast<Eigen::Index>(numNodes), 3);
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numNodes); ++i) {
    xComMat(i, 0) = xCom[3 * i + 0];
    xComMat(i, 1) = xCom[3 * i + 1];
    xComMat(i, 2) = xCom[3 * i + 2];
    xPinnedMat(i, 0) = xPinned[3 * i + 0];
    xPinnedMat(i, 1) = xPinned[3 * i + 1];
    xPinnedMat(i, 2) = xPinned[3 * i + 2];
    bMat(i, 0) = b[3 * i + 0];
    bMat(i, 1) = b[3 * i + 1];
    bMat(i, 2) = b[3 * i + 2];
  }
  double relResCom = (bMat - L_periodic * xComMat).norm() / bMat.norm();
  double relResPinned = (bMat - L_periodic * xPinnedMat).norm() / bMat.norm();
  std::cout << "COM Actual Rel Residual:    " << relResCom << "\n"
            << "Pinned Actual Rel Residual: " << relResPinned << "\n";

  CHECK(std::abs(xPinned[0]) < 1e-6);
  CHECK(std::abs(xPinned[1]) < 1e-6);
  CHECK(std::abs(xPinned[2]) < 1e-6);
  CHECK(relResPinned <= 1e-3);
  CHECK(timePinned_ms > 0.0);
}

TEST_CASE(
  "MEHPForceBalance2 Metal Hybrid SIMDgroup SpMM Degree-Aware Work-Queue",
  "[MEHPSimulator][Metal][Hybrid][SpMM]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout
    << "\n=== Testing Hybrid SIMDgroup SpMM: Degree-Aware Work-Queue ==="
    << std::endl;
  const size_t numNodes = 4000;
  std::vector<Eigen::Triplet<double>> trips;
  trips.reserve(numNodes * 3);
  for (size_t i = 0; i < numNodes - 1; ++i) {
    trips.emplace_back(i, i + 1, -1.0);
    trips.emplace_back(i + 1, i, -1.0);
  }
  for (size_t i = 0; i < numNodes / 15; ++i) {
    size_t a = (i * 29) % numNodes;
    size_t b = (i * 47 + 100) % numNodes;
    if (a != b) {
      trips.emplace_back(a, b, -1.0);
      trips.emplace_back(b, a, -1.0);
    }
  }

  Eigen::VectorXd diag =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numNodes));
  for (const auto& tr : trips) {
    diag[tr.row()] -= tr.value();
  }
  for (size_t i = 0; i < numNodes; ++i) {
    trips.emplace_back(i, i, diag[static_cast<Eigen::Index>(i)] + 0.1);
  }
  Eigen::SparseMatrix<double> scalarL(static_cast<Eigen::Index>(numNodes),
                                      static_cast<Eigen::Index>(numNodes));
  scalarL.setFromTriplets(trips.begin(), trips.end());
  scalarL.makeCompressed();

  Eigen::VectorXd x =
    Eigen::VectorXd::Random(static_cast<Eigen::Index>(numNodes * 3));

  Eigen::VectorXd yUniform = pcm::spmvMetal3D(scalarL, x);
  Eigen::VectorXd yHybrid = pcm::spmvMetal3DHybrid(scalarL, x, 8);
  double diff = (yUniform - yHybrid).cwiseAbs().maxCoeff();
  std::cout << "Hybrid vs Uniform SpMV Max Diff: " << diff << std::endl;
  CHECK(diff < 1e-4);

  auto bench = pcm::benchmarkHybridSpMVMetal(scalarL, 50);
  std::cout << "Polymer network beads: Low-degree (<4) = "
            << bench.numLowDegreeRows << " ("
            << (100.0 * bench.numLowDegreeRows / numNodes) << "%), "
            << "High-degree (>=4) = " << bench.numHighDegreeRows << "\n"
            << "Uniform Coop2 Time: " << bench.timeCoop2_ms << " ms\n"
            << "Uniform Coop8 Time: " << bench.timeCoop8_ms << " ms\n"
            << "Degree-Aware Hybrid Queue Time: " << bench.timeHybrid_ms
            << " ms\n"
            << "Speedup Hybrid vs Coop8: " << bench.speedupHybridVs8 << "x\n"
            << "==============================================================="
               "==========\n";

  CHECK(bench.numLowDegreeRows >=
        static_cast<uint32_t>(0.80 * static_cast<double>(numNodes)));
  CHECK(bench.timeHybrid_ms > 0.0);
}

TEST_CASE("MEHPForceBalance2 Metal Zero-Copy End-to-End Simulation Pipeline",
          "[MEHPSimulator][Metal][ZeroCopy]")
{
  if (!pcm::isMetalAvailable()) {
    SKIP("Metal GPU backend is not available on this device.");
  }
  std::cout << "\n=== Testing Zero-Copy End-to-End Simulation Pipeline ==="
            << std::endl;
  const size_t numNodes = 6;
  Eigen::SparseMatrix<double> scalarA(static_cast<Eigen::Index>(numNodes),
                                      static_cast<Eigen::Index>(numNodes));
  std::vector<Eigen::Triplet<double>> trips;
  for (size_t i = 0; i < numNodes - 1; ++i) {
    trips.emplace_back(i, i + 1, -1.0);
    trips.emplace_back(i + 1, i, -1.0);
  }
  for (size_t i = 0; i < numNodes; ++i) {
    trips.emplace_back(i, i, 3.0);
  }
  scalarA.setFromTriplets(trips.begin(), trips.end());
  scalarA.makeCompressed();

  Eigen::VectorXd b =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numNodes * 3));
  b[0] = 1.0;
  b[3 * static_cast<Eigen::Index>(numNodes - 1)] = -1.0;

  std::vector<int> sA = { 0, 1, 2, 3, 4 };
  std::vector<int> sB = { 1, 2, 3, 4, 5 };
  std::vector<double> boxOffset(15, 0.0);
  std::vector<double> contour(5, 1.0);

  pcm::SpringRelaxationData springData;
  springData.springIndexA = sA.data();
  springData.springIndexB = sB.data();
  springData.springBoxOffset = boxOffset.data();
  springData.springContourLength = contour.data();
  springData.numSprings = 5;
  springData.breakingThreshold = 100.0;
  springData.inactiveCutoff = 0.0;

  pcm::GpuRelaxationStepResult result =
    pcm::runGpuRelaxationStep(scalarA,
                              b,
                              1e-4,
                              100,
                              springData,
                              Eigen::VectorXd(),
                              false,
                              true,
                              -1,
                              true);

  CHECK(result.converged);
  CHECK(result.numBroken == 0);
  CHECK(result.totalFlagged == 0);

  Eigen::VectorXd coords = pcm::getGpuResidentCoordinates(numNodes * 3);
  CHECK(coords.size() == static_cast<Eigen::Index>(numNodes * 3));
  CHECK(coords.norm() > 0.0);
}

#endif

TEST_CASE("MEHPForceBalance2 in-place CSR matrix update on strand breaking",
          "[MEHPSimulator]")
{
  std::cout << "\n=== Testing In-Place CSR Matrix Update on Strand Breaking ==="
            << std::endl;
  size_t nrOfBeads = 20;
  size_t nrOfChains = 4;
  pe::Universe universe = pe::Universe(10.0, 10.0, 10.0);
  for (size_t i = 0; i < nrOfChains; ++i) {
    long int beadOffset = static_cast<long int>(i * nrOfBeads);
    std::vector<long int> atomIds(nrOfBeads);
    std::iota(atomIds.begin(), atomIds.end(), 1 + beadOffset);
    std::vector<double> coords(nrOfBeads);
    std::iota(
      coords.begin(), coords.end(), 1. + static_cast<double>(beadOffset));
    universe.addAtoms(atomIds,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 1),
                      coords,
                      coords,
                      coords,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0));
    std::vector<long int> bondFrom(nrOfBeads - 1);
    std::vector<long int> bondTo(nrOfBeads - 1);
    std::iota(bondFrom.begin(), bondFrom.end(), 1 + beadOffset);
    std::iota(bondTo.begin(), bondTo.end(), 2 + beadOffset);
    universe.addBonds(bondFrom, bondTo);
  }
  universe.addBonds({ 1, 21, 41 }, { 21, 41, 1 });

  pcm::MEHPSimulator fb(
    universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
  pcm::MEHPLinearSolverConfig cfg;
  cfg.simplificationMode = pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
  cfg.inactiveRemovalCutoff = 1e-6;
  cfg.solver = pcm::SLESolver::CONJUGATE_GRADIENT;
  cfg.springBreakingLength = 5.0;

  size_t initSprings = fb.getNrOfSprings();
  CHECK_NOTHROW(fb.run(cfg));
  CHECK_NOTHROW(fb.validateNetwork());
  std::cout << "Initial springs: " << initSprings
            << ", Final active springs: " << fb.getNrOfActiveSprings()
            << std::endl;
}

TEST_CASE("MEHPForceBalance2 In-Place CSR Updates with Threshold Rebuild",
          "[MEHPSimulator][CSR][threshold]")
{
  pe::Universe universe(30.0, 30.0, 30.0);
  const size_t nrOfChains = 3;
  const size_t nrOfBeads = 20;

  for (size_t chainIdx = 0; chainIdx < nrOfChains; ++chainIdx) {
    long int beadOffset = static_cast<long int>(chainIdx * nrOfBeads);
    std::vector<long int> atomIds(nrOfBeads);
    std::iota(atomIds.begin(), atomIds.end(), 1 + beadOffset);
    std::vector<double> coords(nrOfBeads);
    std::iota(
      coords.begin(), coords.end(), 1.0 + static_cast<double>(beadOffset));
    universe.addAtoms(atomIds,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 1),
                      coords,
                      coords,
                      coords,
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0),
                      pylimer_tools::utils::initializeWithValue(nrOfBeads, 0));
    std::vector<long int> bondFrom(nrOfBeads - 1);
    std::vector<long int> bondTo(nrOfBeads - 1);
    std::iota(bondFrom.begin(), bondFrom.end(), 1 + beadOffset);
    std::iota(bondTo.begin(), bondTo.end(), 2 + beadOffset);
    universe.addBonds(bondFrom, bondTo);
  }
  universe.addBonds(std::vector<long int>{ 1, 21, 41 },
                    std::vector<long int>{ 21, 41, 1 });

  SECTION("Threshold configuration getter and setter")
  {
    pcm::MEHPLinearSolverConfig cfg;
    CHECK(cfg.inPlaceRebuildThreshold == Catch::Approx(0.15));
    pcm::LinearSystemMEHPSolver solver(cfg);
    CHECK(solver.getConfig().inPlaceRebuildThreshold == Catch::Approx(0.15));
    cfg.inPlaceRebuildThreshold = 0.20;
    solver.setConfig(cfg);
    CHECK(solver.getConfig().inPlaceRebuildThreshold == Catch::Approx(0.20));
  }

  SECTION("In-place update with high threshold retains numerical zeros")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = pcm::SLESolver::SIMPLICIAL_LDLT;
    cfg.springBreakingLength = 5.0;
    cfg.inPlaceRebuildThreshold = 0.90;
    CHECK_NOTHROW(fb.run(cfg));
    CHECK_NOTHROW(fb.validateNetwork());
  }

  SECTION("Threshold trigger rebuilds CSR and re-factorizes cleanly")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = pcm::SLESolver::SIMPLICIAL_LDLT;
    cfg.springBreakingLength = 5.0;
    cfg.inPlaceRebuildThreshold = 0.05;
    CHECK_NOTHROW(fb.run(cfg));
    CHECK_NOTHROW(fb.validateNetwork());
  }
}

#ifdef CEREALIZABLE
TEST_CASE("MEHPForceBalance2 Cereal Serialization (1A.10)",
          "[bugfixes][MEHPSimulator][cereal]")
{
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));

  pe::UniverseSequence universeSeq;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);

  SECTION("Default Hookean evaluator round-trip")
  {
    pcm::MEHPSimulator fb =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe, 250, 2.0, 0.0, 100);
    REQUIRE(fb.getNrOfLinks() > fb.getNrOfNodes());
    fb.setForceEvaluator(std::make_shared<pcm::HookeanForceEvaluator>(2.5));
    pcm::MEHPLinearSolverConfig cfg;
    fb.run(cfg);

    std::string serialized = fb.serializeToString();
    CHECK(!serialized.empty());

    pcm::MEHPSimulator restored =
      pcm::MEHPSimulator::constructFromString(serialized);
    CHECK(restored.getNrOfNodes() == fb.getNrOfNodes());
    CHECK(restored.getNrOfLinks() == fb.getNrOfLinks());
    CHECK(restored.getNrOfStrands() == fb.getNrOfStrands());
    CHECK(restored.getNrOfSprings() == fb.getNrOfSprings());
    CHECK(restored.getCurrentDisplacements().size() ==
          fb.getCurrentDisplacements().size());
    CHECK((restored.getCurrentDisplacements() - fb.getCurrentDisplacements())
            .norm() < 1e-12);
    CHECK(restored.getForceEvaluator()->isLinear());
    CHECK(restored.getForceEvaluator()->getKappa() == Catch::Approx(2.5));
    CHECK(restored.getEnergy() == Catch::Approx(fb.getEnergy()));
    CHECK(restored.getStressTensor().isApprox(fb.getStressTensor(), 1e-12));
    CHECK(restored.getLastSolveResult().exitReason ==
          fb.getLastSolveResult().exitReason);
    CHECK_NOTHROW(restored.validateNetwork());
  }

  SECTION("NonGaussianSpringForceEvaluator round-trip with fixed entanglements")
  {
    pcm::MEHPSimulator fb =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        30,
        2.5,
        0.0,
        5,
        2.0,
        "cereal_ng_seed",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        std::make_shared<pcm::NonGaussianSpringForceEvaluator>(1.7, 0.9));
    pcm::MEHPNonlinearSolverConfig nlCfg;
    nlCfg.maxNrOfSteps = 200;
    fb.run(nlCfg);

    std::string serialized = fb.serializeToString();
    CHECK(!serialized.empty());

    pcm::MEHPSimulator restored =
      pcm::MEHPSimulator::constructFromString(serialized);
    CHECK_FALSE(restored.getForceEvaluator()->isLinear());
    CHECK(restored.getForceEvaluator()->getKappa() == Catch::Approx(1.7));
    CHECK(restored.getNrOfLinks() == fb.getNrOfLinks());
    CHECK(restored.getNrOfSprings() == fb.getNrOfSprings());
    CHECK(restored.getEnergy() == Catch::Approx(fb.getEnergy()));
    CHECK(restored.getStressTensor().isApprox(fb.getStressTensor(), 1e-12));
    CHECK(restored.getResidualNorm() == Catch::Approx(fb.getResidualNorm()));
    CHECK_NOTHROW(restored.validateNetwork());
  }
}
#endif

TEST_CASE("MEHPForceBalance2 2D Mode and Gamma Normalization (1A.5, 1A.8)",
          "[bugfixes][MEHPSimulator][2D]")
{
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::string inputFile =
    suspectedPath + "/structure/network_100_a_46.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));

  pe::UniverseSequence universeSeq;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);

  SECTION("2D simulation keeps z-displacements zero (1A.5)")
  {
    pcm::MEHPSimulator fb2D(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, true);
    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = pcm::SLESolver::DEFAULT;
    fb2D.run(cfg);
    const Eigen::VectorXd u = fb2D.getCurrentDisplacements();
    for (Eigen::Index link = 0;
         link < static_cast<Eigen::Index>(fb2D.getNrOfLinks());
         ++link) {
      CHECK(std::abs(u[3 * link + 2]) < 1e-12);
    }
  }

  SECTION("Gamma factor normalization with slip-links (1A.8)")
  {
    pcm::MEHPSimulator fbWithLinks =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe, 250, 2.0, 0.0, 100);
    REQUIRE(fbWithLinks.getNrOfSprings() >
            fbWithLinks.getNrOfStrands()); // partial springs > strands
    double gammaAuto = fbWithLinks.getGammaFactor(1.0, -1);
    Eigen::VectorXd gammaVec = fbWithLinks.getGammaFactors(1.0);
    double expectedGamma =
      gammaVec.sum() / static_cast<double>(fbWithLinks.getNrOfStrands());
    CHECK(gammaAuto == Catch::Approx(expectedGamma));
    CHECK(gammaAuto > gammaVec.mean());
  }
}

TEST_CASE("MEHPForceBalance2 Adaptive Solver Fallback (CG to BiCGSTAB)",
          "[MEHPSimulator][solver][fallback][CG][BiCGSTAB]")
{
  std::string suspectedPath = PYLIMER_TEST_FIXTURES_DIR;
  std::string inputFile =
    suspectedPath +
    "/structure/"
    "crosslinked_p_1_0.5_melt_100_a_158_100_xlinks_v_13.V-fixed.structure.out-"
    "equilibration_do_crosslink.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));

  pe::UniverseSequence universeSeq;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  pe::Universe universe = universeSeq.atIndex(0);

  SECTION("Configuration getter and setter")
  {
    pcm::MEHPLinearSolverConfig cfg;
    CHECK(cfg.cgFallbackMaxIterations == 20000);
    pcm::LinearSystemMEHPSolver solver(cfg);
    CHECK(solver.getConfig().cgFallbackMaxIterations == 20000);
    cfg.cgFallbackMaxIterations = 15000;
    solver.setConfig(cfg);
    CHECK(solver.getConfig().cgFallbackMaxIterations == 15000);
    cfg.cgFallbackMaxIterations = 0;
    solver.setConfig(cfg);
    CHECK(solver.getConfig().cgFallbackMaxIterations == 0);
  }

  SECTION("Normal system converges with CG without triggering fallback")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    REQUIRE(fb.getNrOfLinks() <= fb.getNrOfNodes());

    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = pcm::SLESolver::DEFAULT;
    cfg.residualReduction = 1e-12;
    cfg.maxIterations = 50000;
    fb.run(cfg);

    CHECK(!fb.getLastSolveResult().cgFallbackTriggered);
    CHECK(fb.getLastSolveResult().effectiveSLESolver ==
          pcm::SLESolver::CONJUGATE_GRADIENT);
    CHECK(fb.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  }

  SECTION("Default solver is ConjugateGradient even when links > nodes")
  {
    pcm::MEHPSimulator fb =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe, 50, 4.0, 0.0, 50, 3.0, "12345");
    REQUIRE(fb.getNetwork().nrOfLinks > fb.getNetwork().nrOfNodes);

    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = pcm::SLESolver::DEFAULT;
    cfg.residualReduction = 1e-12;
    cfg.maxIterations = 50000;
    fb.run(cfg);

    CHECK((fb.getLastSolveResult().effectiveSLESolver ==
             pcm::SLESolver::CONJUGATE_GRADIENT ||
           (fb.getLastSolveResult().cgFallbackTriggered &&
            fb.getLastSolveResult().effectiveSLESolver ==
              pcm::SLESolver::BICGSTAB)));
    CHECK(fb.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  }

  SECTION("Fallback triggers when CG reaches iteration limit on system")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    REQUIRE(fb.getNrOfLinks() <= fb.getNrOfNodes());

    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = pcm::SLESolver::DEFAULT;
    cfg.residualReduction = 1e-8;
    cfg.maxIterations = 50000;
    cfg.cgFallbackMaxIterations = 10;
    fb.run(cfg);

    CHECK(fb.getLastSolveResult().cgFallbackTriggered);
    CHECK(fb.getLastSolveResult().effectiveSLESolver ==
          pcm::SLESolver::BICGSTAB);
    CHECK(fb.getExitReason() == pcm::ExitReason::X_TOLERANCE);
  }

  SECTION("Fallback disabled when cgFallbackMaxIterations <= 0")
  {
    pcm::MEHPSimulator fb(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    REQUIRE(fb.getNrOfLinks() <= fb.getNrOfNodes());

    pcm::MEHPLinearSolverConfig cfg;
    cfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    cfg.inactiveRemovalCutoff = 1e-6;
    cfg.solver = pcm::SLESolver::DEFAULT;
    cfg.residualReduction = 1e-12;
    cfg.maxIterations = 50000;
    cfg.cgFallbackMaxIterations = 0;
    fb.run(cfg);

    CHECK(!fb.getLastSolveResult().cgFallbackTriggered);
    CHECK(fb.getLastSolveResult().effectiveSLESolver ==
          pcm::SLESolver::CONJUGATE_GRADIENT);
  }

  SECTION("Ill-conditioned network is rescued by BiCGSTAB instead of stalling")
  {
    pcm::MEHPSimulator fbStalled(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    REQUIRE(fbStalled.getNrOfLinks() <= fbStalled.getNrOfNodes());

    pcm::MEHPLinearSolverConfig stalledCfg;
    stalledCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    stalledCfg.inactiveRemovalCutoff = 1e-6;
    stalledCfg.solver = pcm::SLESolver::DEFAULT;
    stalledCfg.residualReduction = 1e-12;
    stalledCfg.maxIterations = 25;
    stalledCfg.cgFallbackMaxIterations = 0;
    fbStalled.run(stalledCfg);
    CHECK(!fbStalled.getLastSolveResult().cgFallbackTriggered);
    CHECK(fbStalled.getLastSolveResult().effectiveSLESolver ==
          pcm::SLESolver::CONJUGATE_GRADIENT);
    double stalledResidual = fbStalled.getResidualNorm();
    CHECK(stalledResidual > 1e-10);

    pcm::MEHPSimulator fbRescued(
      universe, {}, 2, pcm::EntanglementMode::MERGED_LINKS, nullptr, false);
    pcm::MEHPLinearSolverConfig rescuedCfg;
    rescuedCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;
    rescuedCfg.inactiveRemovalCutoff = 1e-6;
    rescuedCfg.solver = pcm::SLESolver::DEFAULT;
    rescuedCfg.residualReduction = 1e-8;
    rescuedCfg.maxIterations = 50000;
    rescuedCfg.cgFallbackMaxIterations = 25;
    fbRescued.run(rescuedCfg);
    CHECK(fbRescued.getLastSolveResult().cgFallbackTriggered);
    CHECK(fbRescued.getLastSolveResult().effectiveSLESolver ==
          pcm::SLESolver::BICGSTAB);
    CHECK(fbRescued.getExitReason() == pcm::ExitReason::X_TOLERANCE);
    double rescuedResidual = fbRescued.getResidualNorm();
    CHECK(rescuedResidual < 1e-6);
    CHECK(rescuedResidual * 1000.0 < stalledResidual);
  }
}

TEST_CASE("MEHPSimulator cross-solver parity on phantom and fixed-entanglement "
          "networks",
          "[analysis][MEHPSimulator]")
{
  const std::string inputFile = std::string(PYLIMER_TEST_FIXTURES_DIR) +
                                "/structure/network_100_a_46.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));

  pe::UniverseSequence universeSeq;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  const pe::Universe universe = universeSeq.atIndex(0);

  auto checkParity = [](const pcm::MEHPNetwork& baseNet) {
    pcm::MEHPSimulator simLinear(baseNet);
    pcm::MEHPSimulator simSlipFixed(baseNet);
    pcm::MEHPSimulator simNonlinear(baseNet);

    pcm::MEHPLinearSolverConfig linCfg;
    linCfg.solver = pcm::SLESolver::DEFAULT;
    linCfg.residualReduction = 1e-10;
    linCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;

    pcm::MEHPSlipLinkSolverConfig slCfg;
    slCfg.slipMode = pcm::SlipMode::FIXED;
    slCfg.maxNrOfSteps = 25000;
    slCfg.xtol = 1e-13;
    slCfg.oneOverSpringPartitionUpperLimit = -1.0;
    slCfg.doInnerIterations = true;
    slCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;

    pcm::MEHPNonlinearSolverConfig nlCfg;
    nlCfg.algorithm = pcm::NLoptAlgorithm::LD_LBFGS;
    nlCfg.maxNrOfSteps = 50000;
    nlCfg.xtol = 1e-14;
    nlCfg.ftol = 1e-14;
    nlCfg.simplificationMode =
      pcm::StructureSimplificationMode::NO_SIMPLIFICATION;

    const pcm::MEHPSolveResult resLin = simLinear.run(linCfg);
    const pcm::MEHPSolveResult resSl = simSlipFixed.run(slCfg);
    const pcm::MEHPSolveResult resNl = simNonlinear.run(nlCfg);

    CHECK(resLin.residual < 1e-6);
    CHECK(resSl.residual < 1e-4);
    CHECK(resNl.residual < 1e-3);

    CHECK_THAT(simSlipFixed.getEnergy(),
               Catch::Matchers::WithinRel(simLinear.getEnergy(), 1e-3));
    CHECK_THAT(simNonlinear.getEnergy(),
               Catch::Matchers::WithinRel(simLinear.getEnergy(), 1e-3));
    CHECK_THAT(simSlipFixed.getGammaFactor(1.0),
               Catch::Matchers::WithinRel(simLinear.getGammaFactor(1.0), 1e-3));
    CHECK_THAT(simNonlinear.getGammaFactor(1.0),
               Catch::Matchers::WithinRel(simLinear.getGammaFactor(1.0), 1e-3));
    CHECK_THAT(simSlipFixed.getPressure(),
               Catch::Matchers::WithinRel(simLinear.getPressure(), 1e-3));
    CHECK_THAT(simNonlinear.getPressure(),
               Catch::Matchers::WithinRel(simLinear.getPressure(), 1e-3));
    CHECK(simSlipFixed.getStressTensor().isApprox(simLinear.getStressTensor(),
                                                  1e-2));
    CHECK(simNonlinear.getStressTensor().isApprox(simLinear.getStressTensor(),
                                                  1e-2));
  };

  SECTION("Phantom network")
  {
    const pcm::MEHPNetwork net = pcm::MEHPNetwork::fromUniverse(
      universe, {}, 2, pcm::EntanglementMode::NONE);
    checkParity(net);
  }

  SECTION("Fixed merged-link entanglement network")
  {
    const pcm::MEHPNetwork net =
      pcm::MEHPNetwork::fromUniverseWithRandomEntanglements(
        universe,
        30,
        2.5,
        0.0,
        10,
        2.0,
        "parity_seed",
        2,
        pcm::EntanglementMode::MERGED_LINKS);
    REQUIRE(net.nrOfLinks > net.nrOfNodes);
    REQUIRE(net.nrOfSprings > net.nrOfStrands);
    checkParity(net);
  }
}

TEST_CASE("MEHPForceEvaluator stiffness and NonGaussian fixed-entanglement "
          "solve",
          "[analysis][MEHPSimulator][NonGaussianSpringForceEvaluator]")
{
  SECTION("Analytical spring stiffness matches central finite differences")
  {
    const pcm::HookeanForceEvaluator hookean(2.3);
    const pcm::NonGaussianSpringForceEvaluator langevin(1.7, 0.8);
    const Eigen::Vector3d r(1.2, -0.7, 0.4);
    const double N = 4.0;
    const double h = 1e-6;

    for (const pcm::MEHPForceEvaluator* evalPtr :
         { static_cast<const pcm::MEHPForceEvaluator*>(&hookean),
           static_cast<const pcm::MEHPForceEvaluator*>(&langevin) }) {
      const Eigen::Matrix3d K = evalPtr->evaluateSpringStiffness(r, N);
      Eigen::Matrix3d Kfd = Eigen::Matrix3d::Zero();
      for (Eigen::Index col = 0; col < 3; ++col) {
        Eigen::Vector3d dr = Eigen::Vector3d::Zero();
        dr[col] = h;
        const Eigen::Vector3d fPlus = evalPtr->evaluateSpringForce(r + dr, N);
        const Eigen::Vector3d fMinus = evalPtr->evaluateSpringForce(r - dr, N);
        Kfd.col(col) = (fPlus - fMinus) / (2.0 * h);
      }
      CHECK(K.isApprox(Kfd, 1e-6));
    }
  }

  SECTION("NonGaussian force evaluator on fixed merged-link entanglements")
  {
    const std::string inputFile = std::string(PYLIMER_TEST_FIXTURES_DIR) +
                                  "/structure/network_100_a_46.structure.out";
    REQUIRE(std::filesystem::exists(inputFile));

    pe::UniverseSequence universeSeq;
    universeSeq.initializeFromDataSequence({ { inputFile } });
    const pe::Universe universe = universeSeq.atIndex(0);

    const auto evaluator =
      std::make_shared<pcm::NonGaussianSpringForceEvaluator>(1.0, 1.0);
    CHECK_FALSE(evaluator->isLinear());

    pcm::MEHPSimulator sim =
      pcm::MEHPSimulator::constructWithRandomEntanglements(
        universe,
        25,
        2.5,
        0.0,
        5,
        2.0,
        "nongauss_ent",
        2,
        pcm::EntanglementMode::MERGED_LINKS,
        true,
        false,
        evaluator);
    REQUIRE(sim.getNrOfLinks() > sim.getNrOfNodes());
    REQUIRE(sim.getNrOfSprings() > sim.getNrOfStrands());

    const double initialEnergy = sim.getEnergy();
    const double initialResidual = sim.getResidualNorm();

    pcm::MEHPNonlinearSolverConfig nlCfg;
    nlCfg.algorithm = pcm::NLoptAlgorithm::LD_MMA;
    nlCfg.maxNrOfSteps = 20000;
    nlCfg.xtol = 1e-10;
    nlCfg.ftol = 1e-10;

    const pcm::MEHPSolveResult res = sim.run(nlCfg);
    CHECK(res.nrOfIterations > 0);
    CHECK(sim.getEnergy() <= initialEnergy);
    CHECK(sim.getResidualNorm() < initialResidual);
    CHECK_NOTHROW(sim.validateNetwork());
  }
}

TEST_CASE("MEHPSolver precondition validation throws std::invalid_argument",
          "[analysis][MEHPSimulator][exceptions]")
{
  const std::string inputFile = std::string(PYLIMER_TEST_FIXTURES_DIR) +
                                "/structure/network_100_a_46.structure.out";
  REQUIRE(std::filesystem::exists(inputFile));

  pe::UniverseSequence universeSeq;
  universeSeq.initializeFromDataSequence({ { inputFile } });
  const pe::Universe universe = universeSeq.atIndex(0);

  const pcm::HookeanForceEvaluator hookeanEval(1.0);
  const pcm::NonGaussianSpringForceEvaluator nonLinearEval(1.0, 1.0);

  pcm::MEHPNetwork mergedNet =
    pcm::MEHPNetwork::fromUniverseWithRandomEntanglements(
      universe,
      10,
      2.5,
      0.0,
      2,
      2.0,
      "precond_seed",
      2,
      pcm::EntanglementMode::MERGED_LINKS);
  pcm::MEHPNetwork springNet =
    pcm::MEHPNetwork::fromUniverseWithRandomEntanglements(
      universe,
      10,
      2.5,
      0.0,
      2,
      2.0,
      "precond_seed",
      2,
      pcm::EntanglementMode::SPRING);
  REQUIRE(springNet.springIsEntanglement.any());

  Eigen::VectorXd uMerged =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * mergedNet.nrOfLinks));
  Eigen::VectorXd uSpring =
    Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * springNet.nrOfLinks));

  SECTION("LinearSystemMEHPSolver preconditions")
  {
    pcm::LinearSystemMEHPSolver linDefault;
    CHECK_THROWS_AS(linDefault.solve(mergedNet, uMerged, nonLinearEval, false),
                    std::invalid_argument);

    pcm::MEHPLinearSolverConfig slidingCfg;
    slidingCfg.slipMode = pcm::SlipMode::SLIDING;
    pcm::LinearSystemMEHPSolver linSliding(slidingCfg);
    CHECK_THROWS_AS(linSliding.solve(mergedNet, uMerged, hookeanEval, false),
                    std::invalid_argument);
  }

  SECTION("SlipLinkRelaxationSolver preconditions")
  {
    pcm::SlipLinkRelaxationSolver slSolver;
    CHECK_THROWS_AS(slSolver.solve(mergedNet, uMerged, nonLinearEval, false),
                    std::invalid_argument);
    CHECK_THROWS_AS(slSolver.solve(springNet, uSpring, hookeanEval, false),
                    std::invalid_argument);
  }

  SECTION("NonlinearOptimizationSolver preconditions")
  {
    pcm::MEHPNonlinearSolverConfig slidingCfg;
    slidingCfg.slipMode = pcm::SlipMode::SLIDING;
    pcm::NonlinearOptimizationSolver nlSolver(slidingCfg);
    CHECK_THROWS_AS(nlSolver.solve(mergedNet, uMerged, nonLinearEval, false),
                    std::invalid_argument);
  }
}

TEST_CASE("MEHPNetwork merges periodic primary-loop and PBC-wrapping springs",
          "[analysis][MEHPNetwork][simplify]")
{
  const pe::Box box(10.0, 10.0, 10.0);

  SECTION("Periodic primary loop from unlinkBifunctionalLinks (keptIsA = true)")
  {
    pcm::CrosslinkerUniverse cu;
    cu.xlinkTypes = { 2, 2, 2 };
    cu.xlinkX = { 1.0, 5.0, 1.0 };
    cu.xlinkY = { 5.0, 5.0, 1.0 };
    cu.xlinkZ = { 5.0, 5.0, 5.0 };
    cu.strandFrom = { 0, 1, 0 };
    cu.strandTo = { 1, 0, 2 };
    cu.beadsInStrand = { 10, 10, 10 };

    pcm::MEHPNetwork net = pcm::MEHPNetwork::fromCrosslinkerUniverse(
      box,
      cu,
      {},
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      false,
      true,
      false,
      false);
    net.springBoxOffset.segment<3>(3 * 1) = Eigen::Vector3d(10.0, 0.0, 0.0);
    Eigen::VectorXd u =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    REQUIRE_NOTHROW(net.validate(u));

    size_t nRemoved = 0;
    REQUIRE_NOTHROW(nRemoved = net.unlinkBifunctionalLinks(u, false));
    CHECK(nRemoved == 1);
    REQUIRE_NOTHROW(net.validate(u));
    CHECK(net.nrOfSprings == 2);

    const Eigen::Vector3d loopVec = net.evaluateSpringVector(u, 0);
    CHECK_THAT(loopVec.norm(), Catch::Matchers::WithinAbs(10.0, 1e-9));
    CHECK_THAT(loopVec[0], Catch::Matchers::WithinAbs(10.0, 1e-9));
    CHECK_THAT(loopVec[1], Catch::Matchers::WithinAbs(0.0, 1e-9));
    CHECK_THAT(loopVec[2], Catch::Matchers::WithinAbs(0.0, 1e-9));
  }

  SECTION(
    "Periodic primary loop from unlinkBifunctionalLinks (keptIsA = false)")
  {
    pcm::CrosslinkerUniverse cu;
    cu.xlinkTypes = { 2, 2, 2 };
    cu.xlinkX = { 1.0, 5.0, 1.0 };
    cu.xlinkY = { 5.0, 5.0, 1.0 };
    cu.xlinkZ = { 5.0, 5.0, 5.0 };
    cu.strandFrom = { 1, 0, 0 };
    cu.strandTo = { 0, 1, 2 };
    cu.beadsInStrand = { 10, 10, 10 };

    pcm::MEHPNetwork net = pcm::MEHPNetwork::fromCrosslinkerUniverse(
      box,
      cu,
      {},
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      false,
      true,
      false,
      false);
    net.springBoxOffset.segment<3>(3 * 0) = Eigen::Vector3d(10.0, 0.0, 0.0);
    Eigen::VectorXd u =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    REQUIRE_NOTHROW(net.validate(u));

    size_t nRemoved = 0;
    REQUIRE_NOTHROW(nRemoved = net.unlinkBifunctionalLinks(u, false));
    CHECK(nRemoved == 1);
    REQUIRE_NOTHROW(net.validate(u));
    CHECK(net.nrOfSprings == 2);

    const Eigen::Vector3d loopVec = net.evaluateSpringVector(u, 0);
    CHECK_THAT(loopVec.norm(), Catch::Matchers::WithinAbs(10.0, 1e-9));
    CHECK_THAT(loopVec[0], Catch::Matchers::WithinAbs(10.0, 1e-9));
    CHECK_THAT(loopVec[1], Catch::Matchers::WithinAbs(0.0, 1e-9));
    CHECK_THAT(loopVec[2], Catch::Matchers::WithinAbs(0.0, 1e-9));
  }

  SECTION("PBC-wrapping spring merge with assumeBoxLargeEnough = true")
  {
    pcm::CrosslinkerUniverse cu;
    cu.xlinkTypes = { 2, 2, 2 };
    cu.xlinkX = { 1.0, 5.5, 0.0 };
    cu.xlinkY = { 5.0, 5.0, 5.0 };
    cu.xlinkZ = { 5.0, 5.0, 5.0 };
    cu.strandFrom = { 0, 1 };
    cu.strandTo = { 1, 2 };
    cu.beadsInStrand = { 10, 10 };

    pcm::MEHPNetwork net = pcm::MEHPNetwork::fromCrosslinkerUniverse(
      box,
      cu,
      {},
      2,
      pcm::EntanglementMode::MERGED_LINKS,
      false,
      true,
      false,
      true);
    Eigen::VectorXd u =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
    REQUIRE_NOTHROW(net.validate(u));

    bool merged = false;
    REQUIRE_NOTHROW(merged = net.mergeSpringsAtLink(1, u, true));
    CHECK(merged);
    REQUIRE_NOTHROW(net.validate(u));
    CHECK(net.nrOfSprings == 1);

    const Eigen::Vector3d mergedVec = net.evaluateSpringVector(u, 0);
    CHECK_THAT(mergedVec[0], Catch::Matchers::WithinAbs(-1.0, 1e-9));
    CHECK_THAT(mergedVec[1], Catch::Matchers::WithinAbs(0.0, 1e-9));
    CHECK_THAT(mergedVec[2], Catch::Matchers::WithinAbs(0.0, 1e-9));
  }
}
