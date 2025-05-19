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
MEHPForceRelaxation::runForceRelaxation(
  const char* algorithm,
  const long int maxNrOfSteps, // default: 10000
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
      lowerBounds.push_back(
        -net.L[dir]); //  * 0.5 -> lead to some few atoms not being where
      //  they should. Maybe one box is still not enough?!?
      upperBounds.push_back(net.L[dir]); //  * 0.5
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
  for (size_t i = 0; i < u0.size(); ++i) {
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

Eigen::VectorXd
MEHPForceRelaxation::evaluateSpringDistances(const Network* net,
                                             const bool is2D)
{
  const Eigen::VectorXd u = Eigen::VectorXd::Zero(net->coordinates.size());
  return MEHPForceRelaxation::evaluateSpringDistances(net, u, is2D);
}

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
      springDistances[i] = 0.0;
    }
  }
  assert(springDistances.size() == net->nrOfSprings * 3);

  return springDistances;
}

/**
 * FORCE RELAXATION DATA ACCESS
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
  for (int i = 0; i < this->forceRelaxationNetwork.nrOfNodes; ++i) {
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
      this->forceRelaxationNetwork
        .oldAtomIds[this->forceRelaxationNetwork.springIndexA[i]]);
    bondTo.push_back(
      this->forceRelaxationNetwork
        .oldAtomIds[this->forceRelaxationNetwork.springIndexB[i]]);
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
      r2local += this->currentSpringDistances[i * 3 + j] *
                 this->currentSpringDistances[i * 3 + j];
    }
    r2 += sqrt(r2local);
  }
  return r2 / this->forceRelaxationNetwork.nrOfSprings;
}

/**
 * @brief Compute the stress tensor
 *
 * @param net
 * @param u
 * @return std::array<std::array<double, 3>, 3>
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
    double s[3] = { springDistances[3 * i + 0],
                    springDistances[3 * i + 1],
                    springDistances[3 * i + 2] };
    /* spring contribution to the overall stress tensor */
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
 * @brief Compute the stress tensor
 *
 * @param net
 * @param u
 * @param loopTol
 * @return std::array<std::array<double, 3>, 3>
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
    results.emplace(this->forceRelaxationNetwork.oldAtomIds[i],
                    nrOfActiveSpringsConnected[i]);
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
  for (size_t i = 0; i < this->forceRelaxationNetwork.nrOfNodes; i++) {
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
  for (size_t i = 0; i < this->forceRelaxationNetwork.nrOfSprings; i++) {
    if (springIsActive[i] == true) { /* active spring */
      const int a = this->forceRelaxationNetwork.springIndexA[i];
      const int b = this->forceRelaxationNetwork.springIndexB[i];
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
  for (size_t i = 0; i < this->forceRelaxationNetwork.nrOfNodes * 3; ++i) {
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
