#include "MEHPForceRelaxation.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <iostream>
#include <nlopt.hpp>
#include <random>
#include <string>
#include <vector>

namespace pylimer_tools::sim::mehp {

/**
 * FORCE RELAXATION
 */
void
MEHPForceRelaxation::runForceRelaxation(const char* algorithm,
                                        const long int maxNrOfSteps,
                                        // default: 10000
                                        const double xtol,
                                        const double ftol)
{
  this->simulationHasRun = true;
  RUNTIME_EXP_IFN(this->forceEvaluator != nullptr,
                  "Force evaluator is not set");
  this->forceEvaluator->setNetwork(this->forceRelaxationNetwork);
  this->forceEvaluator->setIs2D(this->is2D);
  this->forceEvaluator->prepareForEvaluations();
  double stress[3][3];

  for (size_t j = 0; j < 3; j++) {
    for (size_t k = 0; k < 3; k++) {
      stress[j][k] = 0.;
    }
  }

  const Network net = this->forceRelaxationNetwork;
  const bool is2D = this->is2D;

  /* array allocation */
  std::vector<double> u0 =
    pylimer_tools::utils::initializeWithValue(3 * net.nrOfNodes, 0.0);

  /* force relaxation */
  nlopt::opt opt(algorithm, 3 * net.nrOfNodes);

  const nlopt::func objectiveF = [](const unsigned n,
                                    const double* x,
                                    double* grad,
                                    void* f_data) -> double {
    const MEHPForceEvaluator* fEvaluator =
      static_cast<MEHPForceEvaluator*>(f_data);
    return fEvaluator->evaluateForceSetGradient(n, x, grad, f_data);
  };
  opt.set_min_objective(objectiveF, this->forceEvaluator);
  // set constraints to support more algorithms
  std::vector<double> upperBounds;
  upperBounds.reserve(3 * net.nrOfNodes);
  std::vector<double> lowerBounds;
  lowerBounds.reserve(3 * net.nrOfNodes);
  for (size_t i = 0; i < net.nrOfNodes; ++i) {
    for (size_t dir = 0; dir < 3; ++dir) {
      lowerBounds.push_back(-net.L[dir]);
      upperBounds.push_back(net.L[dir]);
    }
  }
  opt.set_upper_bounds(upperBounds);
  opt.set_lower_bounds(lowerBounds);
  // set exit conditions
  opt.set_xtol_rel(xtol);
  opt.set_ftol_rel(ftol);
  opt.set_ftol_abs(0.0);
  opt.set_maxeval(maxNrOfSteps);
  // opt.set_param("verbosity", 1.0);
  // start/set/run minimization
  double minf;
  nlopt::result res;
  std::exception_ptr nloptException = nullptr;
  try {
    res = opt.optimize(u0, minf);
  } catch (...) {
    nloptException = std::current_exception();
  }

  // query solution & exit reason
  assert(u0.size() == 3 * net.nrOfNodes);
  bool require_rerun = false;
  for (Eigen::Index i = 0; i < u0.size(); ++i) {
    this->forceRelaxationNetwork.coordinates[i] += u0[i];
    if (!(u0[i] < upperBounds[i] - this->suggestRerunEps) ||
        !(u0[i] > lowerBounds[i] + this->suggestRerunEps)) {
      require_rerun = true;
    }
  }
  this->forceEvaluator->setNetwork(this->forceRelaxationNetwork);
  this->simulationSuggestsRerun = require_rerun;
  this->currentSpringDistances =
    this->evaluateSpringDistances(&this->forceRelaxationNetwork, is2D);

  this->exitReason = ExitReason::OTHER;
  if (nloptException != nullptr) {
    this->exitReason = ExitReason::FAILURE;
    std::cout << "Nlopt exception: " << opt.get_errmsg() << std::endl;
  } else if (res == nlopt::result::FTOL_REACHED) {
    this->exitReason = ExitReason::F_TOLERANCE;
  } else if (res == nlopt::result::XTOL_REACHED) {
    this->exitReason = ExitReason::X_TOLERANCE;
  } else if (res == nlopt::result::MAXEVAL_REACHED) {
    this->exitReason = ExitReason::MAX_STEPS;
  }
  this->nrOfStepsDone += opt.get_numevals();
}

/**
 * @brief Evaluate spring distances for a given network
 *
 * @param net the network to evaluate
 * @param is2D whether to use 2D mode (ignore z-coordinates)
 * @return Eigen::VectorXd the spring distances
 */
Eigen::VectorXd
MEHPForceRelaxation::evaluateSpringDistances(const Network* net,
                                             const bool is2D)
{
  const Eigen::VectorXd u = Eigen::VectorXd::Zero(net->coordinates.size());
  return MEHPForceRelaxation::evaluateSpringDistances(net, u, is2D);
}

/**
 * @brief Evaluate spring distances for a given network with displacements
 *
 * @param net the network to evaluate
 * @param u the displacement vector
 * @param is2D whether to use 2D mode (ignore z-coordinates)
 * @return Eigen::VectorXd the spring distances
 */
Eigen::VectorXd
MEHPForceRelaxation::evaluateSpringDistances(const Network* net,
                                             const Eigen::VectorXd& u,
                                             const bool is2D)
{
  // this is unnecessary overhead :P
  const pylimer_tools::entities::Box box =
    pylimer_tools::entities::Box(net->L[0], net->L[1], net->L[2]);

  // first, the distances
  assert(u.size() == net->coordinates.size());
  Eigen::VectorXd actualCoordinates = net->coordinates + u;
  // It *could* be more efficient to index u instead of the coordinates
  Eigen::VectorXd springDistances =
    (actualCoordinates(net->springCoordinateIndexB) -
     actualCoordinates(net->springCoordinateIndexA)) +
    net->springBoxOffset;

  if (net->assumeBoxLargeEnough) {
    box.handlePBC(springDistances);
  }

  if (is2D) {
    // springDistances(Eigen::seq(2, Eigen::last, Eigen::fix<3>)) =
    //   Eigen::VectorXd::Zero(net->nrOfSprings / 3);
    for (size_t i = 2; i < 3 * net->nrOfSprings; i += 3) {
      springDistances[static_cast<Eigen::Index>(i)] = 0.0;
    }
  }
  assert(springDistances.size() == net->nrOfSprings * 3);

  return springDistances;
}

/**
 * FORCE RELAXATION DATA ACCESS
 */
/**
 * @brief Convert the current network back into a universe, consisting
 * only of crosslinkers
 *
 * @return pylimer_tools::entities::Universe the crosslinker universe
 */
pylimer_tools::entities::Universe
MEHPForceRelaxation::getCrosslinkerVerse() const
{
  // convert nodes & springs back to a universe
  pylimer_tools::entities::Universe xlinkUniverse =
    pylimer_tools::entities::Universe(this->universe.getBox());
  std::vector<long int> ids;
  std::vector<int> types = pylimer_tools::utils::initializeWithValue(
    this->forceRelaxationNetwork.nrOfNodes, crossLinkerType);
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  const std::vector<int> zeros = pylimer_tools::utils::initializeWithValue(
    this->forceRelaxationNetwork.nrOfNodes, 0);
  ids.reserve(this->forceRelaxationNetwork.nrOfNodes);
  x.reserve(this->forceRelaxationNetwork.nrOfNodes);
  y.reserve(this->forceRelaxationNetwork.nrOfNodes);
  z.reserve(this->forceRelaxationNetwork.nrOfNodes);
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->forceRelaxationNetwork.nrOfNodes);
       ++i) {
    x.push_back(this->forceRelaxationNetwork.coordinates[3 * i + 0]);
    y.push_back(this->forceRelaxationNetwork.coordinates[3 * i + 1]);
    z.push_back(this->forceRelaxationNetwork.coordinates[3 * i + 2]);
    ids.push_back(this->forceRelaxationNetwork.oldAtomIds[i]);
    // override type, since the types may be different from crossLinkerType
    // if converted with dangling chains
    types[i] = this->universe.getPropertyValue<int>(
      "type",
      this->universe.getIdxByAtomId(
        this->forceRelaxationNetwork.oldAtomIds[i]));
  }
  xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  bondFrom.reserve(this->forceRelaxationNetwork.nrOfSprings);
  bondTo.reserve(this->forceRelaxationNetwork.nrOfSprings);
  for (int i = 0; i < this->forceRelaxationNetwork.nrOfSprings; ++i) {
    bondFrom.push_back(
      this->forceRelaxationNetwork.oldAtomIds[static_cast<Eigen::Index>(
        this->forceRelaxationNetwork.springIndexA[i])]);
    bondTo.push_back(
      this->forceRelaxationNetwork.oldAtomIds[static_cast<Eigen::Index>(
        this->forceRelaxationNetwork.springIndexB[i])]);
  }
  xlinkUniverse.addBonds(
    bondFrom.size(),
    bondFrom,
    bondTo,
    pylimer_tools::utils::initializeWithValue(bondFrom.size(), 1),
    false,
    false); // disable simplify to keep the self-loops etc.
  return xlinkUniverse;
}

/**
 * @brief Get the Average Spring Length at the current step
 *
 * @return double
 */
double
MEHPForceRelaxation::getAverageSpringLength() const
{
  double r2 = 0.0;
  for (int i = 0; i < this->forceRelaxationNetwork.nrOfSprings; i++) {
    double r2local = 0.0;
    for (int j = 0; j < 3; ++j) {
      r2local +=
        this->currentSpringDistances[static_cast<Eigen::Index>(i * 3 + j)] *
        this->currentSpringDistances[static_cast<Eigen::Index>(i * 3 + j)];
    }
    r2 += sqrt(r2local);
  }
  return r2 / this->forceRelaxationNetwork.nrOfSprings;
}

Eigen::ArrayXb
MEHPForceRelaxation::findActiveSprings(const Network* net,
                                       const double tolerance) const
{
  Eigen::VectorXd springVectors =
    this->evaluateSpringDistances(net, this->is2D);

  Eigen::ArrayXb result =
    Eigen::ArrayXb::Constant(springVectors.size() / 3, false);
  for (size_t i = 0; i < springVectors.size() / 3; ++i) {
    result[i] = !this->distanceIsWithinTolerance(
      springVectors.segment(3 * i, 3), tolerance, net->springsContourLength[i]);
  }
  return result;
}

bool
MEHPForceRelaxation::ConvertNetwork(Network* net,
                                    const int crossLinkerType,
                                    bool remove2functionalCrosslinkers,
                                    bool removeDanglingChains)
{
  std::vector<pylimer_tools::entities::Atom> springEndAtoms =
    this->universe.getAtomsOfType(crossLinkerType);

  if (remove2functionalCrosslinkers) {
    for (pylimer_tools::entities::Atom xlinker : springEndAtoms) {
      // change type of crosslinkers with a degree <= 2 to "normal",
      // non-crosslink beads
      size_t vertexId = this->universe.getIdxByAtomId(xlinker.getId());
      if (this->universe.computeFunctionalityForVertex(vertexId) <= 2) {
        this->universe.setPropertyValue(vertexId, "type", crossLinkerType - 1);
      }
    }
    springEndAtoms = this->universe.getAtomsOfType(crossLinkerType);
  }

  std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
    this->universe.getChainsWithCrosslinker(crossLinkerType);
  net->moleculeIdxToSpring =
    Eigen::VectorXi::Constant(crossLinkerChains.size(), -1);

  // need to include all but dangling and free chains in order to
  // model entanglement
  size_t nrOfSprings = 0;

  std::vector<bool> vertexAdded = pylimer_tools::utils::initializeWithValue(
    this->universe.getNrOfAtoms(), false);
  for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
    std::vector<pylimer_tools::entities::Atom> endAtoms =
      crossLinkerChains[i].getChainEnds(crossLinkerType, true);
    for (pylimer_tools::entities::Atom endAtom : endAtoms) {
      size_t endAtomVertexId =
        static_cast<size_t>(this->universe.getIdxByAtomId(endAtom.getId()));
      if (endAtom.getType() != crossLinkerType &&
          !vertexAdded[endAtomVertexId]) {
        springEndAtoms.push_back(endAtom);
        vertexAdded[endAtomVertexId] = true;
      }
    }
    RUNTIME_EXP_IFN(crossLinkerChains[i].getType() !=
                      pylimer_tools::entities::MoleculeType::UNDEFINED,
                    "Crosslinker chain's chain type could not be "
                    "detected. Cannot work like that.");
    if (crossLinkerChains[i].getType() ==
          pylimer_tools::entities::MoleculeType::NETWORK_STRAND ||
        crossLinkerChains[i].getType() ==
          pylimer_tools::entities::MoleculeType::PRIMARY_LOOP ||
        (!removeDanglingChains &&
         crossLinkerChains[i].getType() ==
           pylimer_tools::entities::MoleculeType::DANGLING_CHAIN)) {
      net->moleculeIdxToSpring[i] = static_cast<int>(nrOfSprings);
      nrOfSprings += 1;
    }
  }

  size_t nrOfSpringEnds = springEndAtoms.size();

  // crossLinkerUniverse.simplify();
  pylimer_tools::entities::Box box = this->universe.getBox();
  net->L[0] = box.getLx();
  net->L[1] = box.getLy();
  net->L[2] = box.getLz();
  net->nrOfNodes = nrOfSpringEnds;
  net->nrOfSprings = nrOfSprings;
  net->coordinates = Eigen::VectorXd::Zero(3 * net->nrOfNodes);
  net->oldAtomIds = Eigen::ArrayXi::Zero(net->nrOfNodes);
  net->springIndexA = Eigen::ArrayXi::Zero(net->nrOfSprings);
  net->springIndexB = Eigen::ArrayXi::Zero(net->nrOfSprings);
  net->springCoordinateIndexA = Eigen::ArrayXi::Zero(3 * net->nrOfSprings);
  net->springBoxOffset = Eigen::VectorXd::Zero(3 * net->nrOfSprings);
  net->springCoordinateIndexB = Eigen::ArrayXi::Zero(3 * net->nrOfSprings);
  net->springIsActive = Eigen::ArrayXb::Constant(net->nrOfSprings, false);
  net->springsContourLength = Eigen::VectorXd::Zero(net->nrOfSprings);
  Eigen::VectorXd targetDistances = Eigen::VectorXd::Zero(3 * net->nrOfSprings);

  net->springIndicesOfLinks.reserve(net->nrOfNodes);
  for (size_t i = 0; i < net->nrOfNodes; ++i) {
    net->springIndicesOfLinks.push_back(std::vector<size_t>());
  }

  // convert beads
  std::map<long int, long int> atomIdToNode;
  for (long int i = 0; i < springEndAtoms.size(); ++i) {
    pylimer_tools::entities::Atom atom = springEndAtoms[i];
    atomIdToNode[atom.getId()] = i;
    net->oldAtomIds[static_cast<Eigen::Index>(i)] = atom.getId();
    net->coordinates[3 * static_cast<Eigen::Index>(i) + 0] = atom.getX();
    net->coordinates[3 * static_cast<Eigen::Index>(i) + 1] = atom.getY();
    net->coordinates[3 * static_cast<Eigen::Index>(i) + 2] = atom.getZ();
  }

  // convert springs
  size_t spring_idx = 0;
  for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
    if (crossLinkerChains[i].getNrOfBonds() == 0) {
      continue;
    }
    std::vector<pylimer_tools::entities::Atom> xlinkersOfChain =
      crossLinkerChains[i].getAtomsOfType(crossLinkerType);
    std::vector<pylimer_tools::entities::Atom> endsOfChain =
      crossLinkerChains[i].getChainEnds(crossLinkerType, true);
    assert(endsOfChain.size() == 2);
    long int nodeIdxFrom = atomIdToNode.at(endsOfChain[0].getId());
    long int nodeIdxTo = atomIdToNode.at(endsOfChain[1].getId());
    bool addChain = false;
    if (crossLinkerChains[i].getType() ==
        pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
      addChain = true;
      // spring contour length = nr of bonds between two crosslinkers
      net->springsContourLength[spring_idx] =
        crossLinkerChains[i].getNrOfBonds();
    } else if (crossLinkerChains[i].getType() ==
               pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
      addChain = true;

      net->springsContourLength[spring_idx] =
        crossLinkerChains[i].getNrOfBonds();
    } else if (crossLinkerChains[i].getType() ==
               pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) {
      if (!removeDanglingChains) {
        // to keep dangling chains, we convert the trailing atom to a
        // crosslink
        net->springsContourLength[spring_idx] =
          crossLinkerChains[i].getNrOfBonds();
        addChain = true;
      }
    }

    if (addChain) {
      RUNTIME_EXP_IFN(net->moleculeIdxToSpring[i] == spring_idx,
                      "Expected spring mapping to be consistent.");
      std::vector<pylimer_tools::entities::Atom> allChainAtoms =
        crossLinkerChains[i].getAtomsLinedUp();
      targetDistances.segment(3 * spring_idx, 3) =
        crossLinkerChains[i].getOverallBondSum(crossLinkerType);

      if (atomIdToNode.at(allChainAtoms[0].getId()) == nodeIdxTo) {
        targetDistances.segment(3 * spring_idx, 3) *= -1;
      } else {
        assert(atomIdToNode.at(allChainAtoms[0].getId()) == nodeIdxFrom);
      }

      pylimer_tools::utils::addIfNotContained(
        net->springIndicesOfLinks[nodeIdxFrom], spring_idx);
      if (nodeIdxFrom != nodeIdxTo) {
        pylimer_tools::utils::addIfNotContained(
          net->springIndicesOfLinks[nodeIdxTo], spring_idx);
      }

      net->springIndexA[spring_idx] = nodeIdxFrom;
      net->springIndexB[spring_idx] = nodeIdxTo;
      for (size_t j = 0; j < 3; j++) {
        net->springCoordinateIndexA[3 * spring_idx + j] =
          static_cast<int>(nodeIdxFrom * 3 + j);
        net->springCoordinateIndexB[3 * spring_idx + j] =
          static_cast<int>(nodeIdxTo * 3 + j);
      }

      spring_idx += 1;
    }
  }

  RUNTIME_EXP_IFN(spring_idx == net->nrOfSprings,
                  "Expected nr of springs not fulfilled.");

  // box volume
  net->vol = net->L[0] * net->L[1] * net->L[2];
  if (net->nrOfSprings > 0) {
    net->meanSpringContourLength = net->springsContourLength.mean();
  }

  net->springBoxOffset = ((net->coordinates(net->springCoordinateIndexA) -
                           net->coordinates(net->springCoordinateIndexB)) +
                          targetDistances);
  // net->springBoxOffset = this->universe.getBox().getOffset(
  //   net->coordinates(net->springCoordinateIndexB) -
  //   net->coordinates(net->springCoordinateIndexA));

  // check whether spring contour lengths are what we want them to be
  size_t numCrosslinkers = springEndAtoms.size();
  // this->universe.countPropertyValue<int>("type", crossLinkerType);
  size_t contourSum = net->springsContourLength.sum();
  assert(net->springsContourLength.size() == net->nrOfSprings);
  assert(numCrosslinkers == net->nrOfNodes);

  return true; // crossLinkerUniverse.getNrOfBonds() == net->nrOfSprings;
}

/**
 * @brief Compute the stress tensor from spring distances and volume
 *
 * @param springDistances the spring distance vectors
 * @param volume the volume of the system
 * @return std::array<std::array<double, 3>, 3> the stress tensor
 */
std::array<std::array<double, 3>, 3>
MEHPForceRelaxation::evaluateStressTensor(
  const Eigen::VectorXd& springDistances,
  const double volume) const
{
  std::array<std::array<double, 3>, 3> stress;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      stress[i][j] = 0.0;
    }
  }

  for (size_t i = 0; i < springDistances.size() / 3; ++i) {
    double s[3] = { springDistances[static_cast<Eigen::Index>(3 * i + 0)],
                    springDistances[static_cast<Eigen::Index>(3 * i + 1)],
                    springDistances[static_cast<Eigen::Index>(3 * i + 2)] };
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        const double contribution =
          this->forceEvaluator->evaluateStressContribution(s, j, k, i);
        stress[j][k] += contribution;
      }
    }
  }

  for (size_t j = 0; j < 3; j++) {
    for (size_t k = 0; k < 3; k++) {
      stress[j][k] /= volume;
    }
  }

  return stress;
}

/**
 * @brief Compute the stress tensor from network, displacements and tolerance
 *
 * @param net the network
 * @param u the displacements
 * @param loopTol the loop tolerance
 * @return std::array<std::array<double, 3>, 3> the stress tensor
 */
std::array<std::array<double, 3>, 3>
MEHPForceRelaxation::evaluateStressTensor(Network* net,
                                          const Eigen::VectorXd& u,
                                          const double loopTol) const
{
  const Eigen::VectorXd springDistances =
    this->evaluateSpringDistances(net, u, this->is2D);

  return this->evaluateStressTensor(springDistances, net->vol);
}

double
MEHPForceRelaxation::computeActiveWeightFraction(Network* net,
                                                 const double tolerance) const
{
  if (net->nrOfSprings < 1) {
    return 0.;
  }
  // find all active springs
  const Eigen::ArrayXb activeSprings = this->findActiveSprings(net, tolerance);
  if (activeSprings.count() == 0) {
    return 1. - this->computeSolubleWeightFraction(net, tolerance);
  }
  // as of now, the springsContourLength is equal to the number of bonds
  // from crosslink to crosslink. therefore, the number of atoms of each
  // of these springs is one less
  Eigen::ArrayXd allActiveAtomsPerChains =
    activeSprings.cast<double>() * (net->springsContourLength.array() -
                                    Eigen::ArrayXd::Ones(net->nrOfSprings));
  // finally, add the cross-links and normalize by the number of atoms.
  return static_cast<double>(allActiveAtomsPerChains.sum() +
                             this->getNrOfActiveNodes()) /
         static_cast<double>(this->universe.getNrOfAtoms());
}

double
MEHPForceRelaxation::countActiveClusteredAtoms(Network* net,
                                               const double tolerance) const
{
  if (net->nrOfSprings < 1) {
    return 0.;
  }

  // if we don't have a universe, we cannot compute the clusters
  // and we should fall back to counting what we have,
  // assuming that the network is complete, i.e., no dangling chains
  // had been omitted
  if (net->assumeComplete) {
    const std::pair<Eigen::ArrayXb, Eigen::ArrayXb> clusteredToActive =
      this->findClusteredToActive(net, tolerance);
    // find all active springs
    const Eigen::ArrayXb activeSprings = clusteredToActive.first;
    assert(activeSprings.size() == net->nrOfSprings);
    const Eigen::ArrayXb nodeIsActive = clusteredToActive.second;
    assert(nodeIsActive.size() == net->nrOfNodes);
    // as of now, the springsContourLength is equal to the number of bonds
    // from crosslink to crosslink. therefore, the number of atoms of each
    // of these springs is one less
    Eigen::ArrayXd allActiveAtomsPerChains =
      activeSprings.cast<double>() * (net->springsContourLength.array() -
                                      Eigen::ArrayXd::Ones(net->nrOfSprings));
    const double activeNodes = nodeIsActive.count();
    return ((allActiveAtomsPerChains).matrix().sum() + activeNodes);
  }

  INVALIDARG_EXP_IFN(
    this->universe.getNrOfAtoms() > net->nrOfNodes,
    "Cannot compute active clustered atoms if the universe does not "
    "correspond to the universe. ");

  // because our internal structure may not contain the full universe,
  // i.e., e.g., dangling chains have been omitted,
  // we need to query the clusters from the universe
  const std::vector<pylimer_tools::entities::Universe> clusters =
    this->universe.getClusters();
  std::vector<long int> vertexIdxToClusterIdx(this->universe.getNrOfAtoms());
  for (size_t i = 0; i < clusters.size(); ++i) {
    for (const pylimer_tools::entities::Atom& atom : clusters[i].getAtoms()) {
      vertexIdxToClusterIdx[this->universe.getIdxByAtomId(atom.getId())] = i;
    }
  }

  std::vector<bool> clusterIsActive(clusters.size(), false);

  // find active atoms
  const std::vector<int> activeNodeIndices =
    this->getIndicesOfActiveNodes(net, tolerance);

  for (const long int& nodeIdx : activeNodeIndices) {
    const long int universeAtomIdx =
      this->universe.getIdxByAtomId(net->oldAtomIds[nodeIdx]);
    clusterIsActive[vertexIdxToClusterIdx[universeAtomIdx]] = true;
  }

  double nClusteredAtoms = 0.;
  for (size_t i = 0; i < clusters.size(); ++i) {
    if (clusterIsActive[i]) {
      nClusteredAtoms += clusters[i].getNrOfAtoms();
    }
  }

  return nClusteredAtoms;
}

std::vector<int>
MEHPForceRelaxation::getIndicesOfActiveNodes(const Network* net,
                                             const double tolerance) const
{
  std::vector<int> results;
  results.reserve(net->nrOfNodes);

  // find all active springs
  Eigen::ArrayXb strandIsActive = this->findActiveSprings(net, tolerance);

  for (Eigen::Index i = 0; i < net->nrOfNodes; i++) {
    for (const int springIdx : net->springIndicesOfLinks[i]) {
      if (strandIsActive[springIdx]) {
        results.push_back(static_cast<int>(i));
        break; // no need to check further springs
      }
    }
  }

  return results;
}

double
MEHPForceRelaxation::computeDanglingWeightFraction(Network* net,
                                                   const double tolerance) const
{
  const double activeWeightFraction =
    this->computeActiveWeightFraction(net, tolerance);
  RUNTIME_EXP_IFN(APPROX_WITHIN(activeWeightFraction, 0., 1., 1e-6),
                  "Expect active weight fraction to be between 0 and 1, got " +
                    std::to_string(activeWeightFraction) + ".");
  const double solubleWeightFraction =
    this->computeSolubleWeightFraction(net, tolerance);
  RUNTIME_EXP_IFN(APPROX_WITHIN(solubleWeightFraction, 0., 1., 1e-6),
                  "Expect soluble weight fraction to be between 0 and 1, got " +
                    std::to_string(solubleWeightFraction) + ".");
  RUNTIME_EXP_IFN(
    APPROX_WITHIN(activeWeightFraction + solubleWeightFraction, 0., 1., 1e-6),
    "Expect active and soluble weight fraction to add up to maximum 1, "
    "got " +
      std::to_string(activeWeightFraction + solubleWeightFraction) + ".");

  // TODO: currently, the weight of the atoms is ignored
  return 1. - activeWeightFraction - solubleWeightFraction;
}

/**
 * @brief Get the Effective Functionality Of each node
 *
 * Returns the number of active springs connected to each atom, atomId
 * used as index
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return std::unordered_map<long int, int>
 */
std::unordered_map<long int, int>
MEHPForceRelaxation::getEffectiveFunctionalityOfAtoms(
  const double tolerance) const
{
  std::unordered_map<long int, int> results;
  results.reserve(this->forceRelaxationNetwork.nrOfNodes);

  Eigen::VectorXi nrOfActiveSpringsConnected =
    this->getNrOfActiveSpringsConnected(tolerance);
  for (size_t i = 0; i < this->forceRelaxationNetwork.nrOfNodes; i++) {
    results.emplace(
      this->forceRelaxationNetwork.oldAtomIds[static_cast<Eigen::Index>(i)],
      nrOfActiveSpringsConnected[static_cast<Eigen::Index>(i)]);
  }
  return results;
}

/**
 * @brief Get the Ids Of active Nodes
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @param minimumNrOfActiveConnections the number of active springs
 * required for this node to qualify as active
 * @return std::vector<long int> the atom ids
 */
std::vector<long int>
MEHPForceRelaxation::getIdsOfActiveNodes(
  const double tolerance,
  const int minimumNrOfActiveConnections,
  const int maximumNrOfActiveConnections) const
{
  std::vector<long int> results;
  results.reserve(this->forceRelaxationNetwork.nrOfNodes);

  Eigen::VectorXi nrOfActiveSpringsConnected =
    this->getNrOfActiveSpringsConnected(tolerance);
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->forceRelaxationNetwork.nrOfNodes);
       i++) {
    if (nrOfActiveSpringsConnected[i] >= minimumNrOfActiveConnections &&
        (maximumNrOfActiveConnections < 0 ||
         maximumNrOfActiveConnections >= nrOfActiveSpringsConnected[i])) {
      results.push_back(this->forceRelaxationNetwork.oldAtomIds[i]);
    }
  }

  return results;
}

/**
 * @brief Get the Nr Of Active Springs connected to each node
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return Eigen::VectorXi
 */
Eigen::VectorXi
MEHPForceRelaxation::getNrOfActiveSpringsConnected(const double tolerance) const
{
  Eigen::VectorXi nrOfActiveSpringsConnected =
    Eigen::VectorXi::Zero(this->forceRelaxationNetwork.nrOfNodes);
  Eigen::ArrayXb springIsActive =
    this->findActiveSprings(&this->forceRelaxationNetwork, tolerance);
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->forceRelaxationNetwork.nrOfSprings);
       i++) {
    if (springIsActive[i] == true) {
      /* active spring */
      const Eigen::Index a =
        static_cast<Eigen::Index>(this->forceRelaxationNetwork.springIndexA[i]);
      const Eigen::Index b =
        static_cast<Eigen::Index>(this->forceRelaxationNetwork.springIndexB[i]);
      ++(nrOfActiveSpringsConnected[a]);
      ++(nrOfActiveSpringsConnected[b]);
    }
  }
  return nrOfActiveSpringsConnected;
}

/**
 * @brief Get the residuals (gradient) at the current step
 *
 * @return Eigen::VectorXd
 */
Eigen::VectorXd
MEHPForceRelaxation::getResiduals() const
{
  double* r = new double[3 * this->forceRelaxationNetwork.nrOfNodes];
  for (size_t i = 0; i < this->forceRelaxationNetwork.nrOfNodes * 3; ++i) {
    r[i] = 0.0;
  }
  try {
    this->forceEvaluator->evaluateForceSetGradient(
      3 * this->forceRelaxationNetwork.nrOfNodes,
      this->currentSpringDistances,
      r);
  } catch (const std::exception& e) {
    delete[] (r);
    throw e;
  }

  Eigen::VectorXd results =
    Eigen::VectorXd::Zero(this->forceRelaxationNetwork.nrOfNodes * 3);
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(
                                 this->forceRelaxationNetwork.nrOfNodes * 3);
       ++i) {
    results[i] = r[i];
  }
  delete[] (r);
  return results;
}

/**
 * @brief Get the Residual Norm at the current step
 *
 * @return double
 */
double
MEHPForceRelaxation::getResidualNorm() const
{
  return this->getResiduals().norm();
}

/**
 * @brief Get the Force at the current step
 *
 * @return double
 */
double
MEHPForceRelaxation::getForce() const
{
  return this->forceEvaluator->evaluateForceSetGradient(
    3 * this->forceRelaxationNetwork.nrOfNodes,
    this->currentSpringDistances,
    nullptr);
}

/**
 * @brief Get the Gamma Factor at the current step
 *
 * @param b02 for the denominator, part of the melt <R_0^2> = b02 *
 * nrOfBondsInSpring
 * @param nrOfChains the nr of chains to average over (can be different
 * from the nr of springs thanks to omitted free chains or primary loops)
 * @return double
 */
double
MEHPForceRelaxation::getGammaFactor(double b02, int nrOfChains) const
{
  if (this->forceRelaxationNetwork.springsContourLength.size() == 0) {
    return 0.0;
  }

  if (b02 < 0) {
    b02 = this->defaultR0Squared /
          this->forceRelaxationNetwork.springsContourLength.mean();
  }
  if (nrOfChains < 1) {
    nrOfChains = this->currentSpringDistances.size() / 3.;
  }

  return this->evaluateGammaFactor(
    this->currentSpringDistances, b02, nrOfChains);
}

/**
 * @brief Get the Gamma Factors at the current step
 *
 * @param b2 the melt b^2 (to go to phantom's Nb^2 for <R_0^2>, using N as
 * the contour length per spring)
 * @return Eigen::VectorXd the gamma factors for each spring
 */
Eigen::VectorXd
MEHPForceRelaxation::getGammaFactors(double b2) const
{
  if (b2 < 0) {
    b2 = this->defaultR0Squared /
         this->forceRelaxationNetwork.springsContourLength.mean();
  }

  return this->evaluateGammaFactors(this->currentSpringDistances, b2);
}
}