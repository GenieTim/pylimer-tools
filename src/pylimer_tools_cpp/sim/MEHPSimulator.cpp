#include "MEHPSimulator.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace pylimer_tools::sim::mehp {

void
MEHPSimulator::initializeAfterNetworkConstruction()
{
  if (this->forceEvaluator == nullptr) {
    this->forceEvaluator = std::make_shared<HookeanForceEvaluator>(1.0);
  }
  this->currentDisplacements = Eigen::VectorXd::Zero(
    static_cast<Eigen::Index>(3 * this->network.nrOfLinks));
  if (this->is2D) {
    for (size_t i = 0; i < this->network.nrOfLinks; ++i) {
      this->network.coordinates[3 * static_cast<Eigen::Index>(i) + 2] = 0.0;
    }
  }
  this->network.validate(this->currentDisplacements);
}

MEHPSimulator::MEHPSimulator(const MEHPNetwork& net,
                             std::shared_ptr<MEHPForceEvaluator> evaluator,
                             const bool is2DIn)
  : OutputSupportingSimulation()
  , network(net)
  , currentDisplacements()
  , forceEvaluator(std::move(evaluator))
  , is2D(is2DIn)
  , crossLinkerType(2)
  , entanglementType(3)
  , lastSolveResult()
  , simulationHasRun(false)
{
  this->initializeAfterNetworkConstruction();
}

MEHPSimulator::MEHPSimulator(
  const pylimer_tools::entities::Universe& universe,
  const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements&
    entanglements,
  const int crossLinkerTypeIn,
  const EntanglementMode entanglementMode,
  std::shared_ptr<MEHPForceEvaluator> evaluator,
  const bool is2DIn,
  const bool keep2FunctionalCrosslinkers,
  const bool removePrimaryLoops,
  const bool assumeBoxLargeEnough)
  : OutputSupportingSimulation()
  , network(MEHPNetwork::fromUniverse(universe,
                                      entanglements,
                                      crossLinkerTypeIn,
                                      entanglementMode,
                                      keep2FunctionalCrosslinkers,
                                      removePrimaryLoops,
                                      assumeBoxLargeEnough))
  , currentDisplacements()
  , forceEvaluator(std::move(evaluator))
  , is2D(is2DIn)
  , crossLinkerType(crossLinkerTypeIn)
  , entanglementType(3)
  , lastSolveResult()
  , simulationHasRun(false)
{
  this->initializeAfterNetworkConstruction();
}

MEHPSimulator::MEHPSimulator(
  const pylimer_tools::entities::Box& box,
  const CrosslinkerUniverse& cu,
  const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
    entanglements,
  const int crossLinkerTypeIn,
  const EntanglementMode entanglementMode,
  std::shared_ptr<MEHPForceEvaluator> evaluator,
  const bool is2DIn,
  const bool removeDanglingChains,
  const bool includeFreeChains,
  const bool removePrimaryLoops,
  const bool assumeBoxLargeEnough)
  : OutputSupportingSimulation()
  , network(MEHPNetwork::fromCrosslinkerUniverse(box,
                                                 cu,
                                                 entanglements,
                                                 crossLinkerTypeIn,
                                                 entanglementMode,
                                                 removeDanglingChains,
                                                 includeFreeChains,
                                                 removePrimaryLoops,
                                                 assumeBoxLargeEnough))
  , currentDisplacements()
  , forceEvaluator(std::move(evaluator))
  , is2D(is2DIn)
  , crossLinkerType(crossLinkerTypeIn)
  , entanglementType(3)
  , lastSolveResult()
  , simulationHasRun(false)
{
  this->initializeAfterNetworkConstruction();
}

MEHPSimulator
MEHPSimulator::constructWithRandomEntanglements(
  const pylimer_tools::entities::Universe& universe,
  const size_t nrOfEntanglementsToSample,
  const double upperSamplingCutoff,
  const double lowerSamplingCutoff,
  const size_t minNrOfEntanglements,
  const double sameStrandCutoff,
  const std::string& seed,
  const int crossLinkerType,
  const EntanglementMode entanglementMode,
  const bool withCrosslinkers,
  const bool filterDanglingAndSoluble,
  std::shared_ptr<MEHPForceEvaluator> evaluator,
  const bool is2D,
  const bool keep2FunctionalCrosslinkers,
  const bool removePrimaryLoops,
  const bool assumeBoxLargeEnough)
{
  const MEHPNetwork net = MEHPNetwork::fromUniverseWithRandomEntanglements(
    universe,
    nrOfEntanglementsToSample,
    upperSamplingCutoff,
    lowerSamplingCutoff,
    minNrOfEntanglements,
    sameStrandCutoff,
    seed,
    crossLinkerType,
    entanglementMode,
    withCrosslinkers,
    filterDanglingAndSoluble,
    keep2FunctionalCrosslinkers,
    removePrimaryLoops,
    assumeBoxLargeEnough);
  MEHPSimulator sim(net, std::move(evaluator), is2D);
  sim.crossLinkerType = crossLinkerType;
  return sim;
}

void
MEHPSimulator::setForceEvaluator(std::shared_ptr<MEHPForceEvaluator> evaluator)
{
  INVALIDARG_EXP_IFN(evaluator != nullptr,
                     "Force evaluator cannot be null.");
  this->forceEvaluator = std::move(evaluator);
}

void
MEHPSimulator::configMeanBondLength(const double meanBondLength)
{
  INVALIDARG_EXP_IFN(std::isfinite(meanBondLength) && meanBondLength >= 0.0,
                     "Mean bond length must be non-negative and finite.");
  this->network.defaultBondLength = meanBondLength;
  this->network.defaultR02 = meanBondLength * meanBondLength;
}

void
MEHPSimulator::configMeanSquaredBondLength(const double r02)
{
  INVALIDARG_EXP_IFN(std::isfinite(r02) && r02 >= 0.0,
                     "Mean squared bond length must be non-negative and finite.");
  this->network.defaultR02 = r02;
  this->network.defaultBondLength = std::sqrt(r02);
}

void
MEHPSimulator::configAssumeBoxLargeEnough(const bool assume)
{
  this->network.assumeBoxLargeEnough = assume;
}

void
MEHPSimulator::configAssumeNetworkIsComplete(const bool assume)
{
  this->network.assumeComplete = assume;
}

void
MEHPSimulator::configRemovePrimaryLoops(const bool remove)
{
  this->network.removePrimaryLoops = remove;
}

void
MEHPSimulator::setIs2D(const bool newIs2D)
{
  this->is2D = newIs2D;
  if (this->is2D) {
    for (size_t i = 0; i < this->network.nrOfLinks; ++i) {
      this->network.coordinates[3 * static_cast<Eigen::Index>(i) + 2] = 0.0;
      if (this->currentDisplacements.size() ==
          static_cast<Eigen::Index>(3 * this->network.nrOfLinks)) {
        this->currentDisplacements[3 * static_cast<Eigen::Index>(i) + 2] = 0.0;
      }
    }
  }
}

void
MEHPSimulator::deformTo(const pylimer_tools::entities::Box& targetBox)
{
  this->network.deformTo(targetBox, this->currentDisplacements);
}

size_t
MEHPSimulator::addSlipLinks(
  const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
    entanglements,
  const bool clampAlpha)
{
  return this->network.addSlipLinks(
    entanglements, &this->currentDisplacements, clampAlpha);
}

size_t
MEHPSimulator::addSlipLinks(const std::vector<size_t>& strandIdx1,
                            const std::vector<size_t>& strandIdx2,
                            const Eigen::VectorXd& x,
                            const Eigen::VectorXd& y,
                            const Eigen::VectorXd& z,
                            const std::vector<double>& alpha1,
                            const std::vector<double>& alpha2,
                            const bool clampAlpha)
{
  return this->network.addSlipLinks(strandIdx1,
                                    strandIdx2,
                                    x,
                                    y,
                                    z,
                                    alpha1,
                                    alpha2,
                                    &this->currentDisplacements,
                                    clampAlpha);
}

size_t
MEHPSimulator::addSlipLinks(const std::vector<size_t>& strandIdx1,
                            const std::vector<size_t>& strandIdx2,
                            const Eigen::VectorXd& x,
                            const Eigen::VectorXd& y,
                            const Eigen::VectorXd& z,
                            const bool clampAlpha)
{
  return this->network.addSlipLinks(strandIdx1,
                                    strandIdx2,
                                    x,
                                    y,
                                    z,
                                    &this->currentDisplacements,
                                    clampAlpha);
}

size_t
MEHPSimulator::randomlyAddSlipLinks(const size_t nrOfEntanglementsToSample,
                                    const double cutoff,
                                    const size_t minNrOfEntanglements,
                                    const double sameStrandCutoff,
                                    const std::string& seed,
                                    const bool clampAlpha)
{
  return this->network.randomlyAddSlipLinks(nrOfEntanglementsToSample,
                                            cutoff,
                                            minNrOfEntanglements,
                                            sameStrandCutoff,
                                            this->crossLinkerType,
                                            seed,
                                            &this->currentDisplacements,
                                            clampAlpha);
}

size_t
MEHPSimulator::breakTooLongStrands(const double breakingLength)
{
  return this->network.breakTooLongStrands(
    this->currentDisplacements, breakingLength, this->is2D);
}

size_t
MEHPSimulator::removeInactiveLinks(const double tolerance)
{
  return this->network.removeInactiveLinks(
    this->currentDisplacements, tolerance, this->is2D);
}

size_t
MEHPSimulator::unlinkBifunctionalLinks(const bool entanglementsOnly)
{
  return this->network.unlinkBifunctionalLinks(this->currentDisplacements,
                                               entanglementsOnly);
}

size_t
MEHPSimulator::simplifyNetwork(const StructureSimplificationMode mode,
                               const double tolerance,
                               const double breakingLength)
{
  return this->network.simplify(
    this->currentDisplacements, mode, tolerance, this->is2D, breakingLength);
}

MEHPSolveResult
MEHPSimulator::run(const MEHPSolver& solver,
                   const MEHPSolveCallbacks& callbacks)
{
  RUNTIME_EXP_IFN(this->forceEvaluator != nullptr,
                  "Cannot run simulation with null force evaluator.");

  MEHPSolveCallbacks wrappedCallbacks = callbacks;
  wrappedCallbacks.onStep = [this, &callbacks](const size_t step) {
    this->lastSolveResult.nrOfIterations = step;
    if (this->hasOutputs()) {
      this->handleOutput(static_cast<long int>(step));
    }
    if (callbacks.onStep) {
      callbacks.onStep(step);
    }
  };

  this->lastSolveResult = solver.solve(this->network,
                                       this->currentDisplacements,
                                       *this->forceEvaluator,
                                       this->is2D,
                                       wrappedCallbacks);
  this->simulationHasRun = true;
  if (this->hasOutputs()) {
    this->handleOutput(static_cast<long int>(this->lastSolveResult.nrOfIterations));
  }
  return this->lastSolveResult;
}

MEHPSolveResult
MEHPSimulator::run(const MEHPLinearSolverConfig& config,
                   const MEHPSolveCallbacks& callbacks)
{
  const LinearSystemMEHPSolver solver(config);
  return this->run(solver, callbacks);
}

MEHPSolveResult
MEHPSimulator::run(const MEHPSlipLinkSolverConfig& config,
                   const MEHPSolveCallbacks& callbacks)
{
  const SlipLinkRelaxationSolver solver(config);
  return this->run(solver, callbacks);
}

MEHPSolveResult
MEHPSimulator::run(const MEHPNonlinearSolverConfig& config,
                   const MEHPSolveCallbacks& callbacks)
{
  const NonlinearOptimizationSolver solver(config);
  return this->run(solver, callbacks);
}

MEHPSolveResult
MEHPSimulator::run()
{
  RUNTIME_EXP_IFN(this->forceEvaluator != nullptr,
                  "Cannot run simulation with null force evaluator.");
  if (this->forceEvaluator->isLinear()) {
    const MEHPLinearSolverConfig linearConfig{};
    return this->run(linearConfig);
  }
  const MEHPNonlinearSolverConfig nonlinearConfig{};
  return this->run(nonlinearConfig);
}

void
MEHPSimulator::setCurrentDisplacements(const Eigen::VectorXd& newDisplacements)
{
  INVALIDARG_EXP_IFN(
    newDisplacements.size() ==
      static_cast<Eigen::Index>(3 * this->network.nrOfLinks),
    "Displacement vector size must equal 3 * nrOfLinks.");
  this->currentDisplacements = newDisplacements;
  if (this->is2D) {
    for (size_t i = 0; i < this->network.nrOfLinks; ++i) {
      this->currentDisplacements[3 * static_cast<Eigen::Index>(i) + 2] = 0.0;
    }
  }
}

double
MEHPSimulator::getAverageSpringLength() const
{
  if (this->network.nrOfSprings == 0) {
    return 0.0;
  }
  return this->getSpringLengths().mean();
}

size_t
MEHPSimulator::getNrOfActiveNodes(const double tolerance) const
{
  return static_cast<size_t>(
    this->network.findActiveNodes(this->currentDisplacements, tolerance, this->is2D)
      .count());
}

size_t
MEHPSimulator::getNrOfActiveSprings(const double tolerance) const
{
  return static_cast<size_t>(
    this->network
      .findActiveSprings(this->currentDisplacements, tolerance, this->is2D)
      .count());
}

size_t
MEHPSimulator::getNrOfActiveStrands(const double tolerance) const
{
  return static_cast<size_t>(
    this->network
      .findActiveStrands(this->currentDisplacements, tolerance, this->is2D)
      .count());
}

double
MEHPSimulator::getEnergy() const
{
  RUNTIME_EXP_IFN(this->forceEvaluator != nullptr,
                  "Force evaluator cannot be null.");
  const Eigen::VectorXd springVectors = this->getSpringVectors();
  return this->forceEvaluator->evaluateNetworkEnergyAndGradient(
    this->network, springVectors, this->is2D, nullptr);
}

Eigen::VectorXd
MEHPSimulator::getResiduals() const
{
  RUNTIME_EXP_IFN(this->forceEvaluator != nullptr,
                  "Force evaluator cannot be null.");
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(
    static_cast<Eigen::Index>(3 * this->network.nrOfLinks));
  const Eigen::VectorXd springVectors = this->getSpringVectors();
  (void)this->forceEvaluator->evaluateNetworkEnergyAndGradient(
    this->network, springVectors, this->is2D, grad.data());
  return -grad;
}

double
MEHPSimulator::getResidualNorm() const
{
  if (this->network.nrOfLinks == 0 || this->network.nrOfSprings == 0) {
    return 0.0;
  }
  return this->getResiduals().norm();
}

Eigen::Matrix3d
MEHPSimulator::getStressTensor()
{
  RUNTIME_EXP_IFN(this->forceEvaluator != nullptr,
                  "Force evaluator cannot be null.");
  if (this->network.nrOfSprings == 0) {
    return Eigen::Matrix3d::Zero();
  }
  // When sliding slip-links were solved with a partition clamp limit, evaluate
  // the stress tensor consistently with the clamped contour partitions.
  if (this->lastSolveResult.oneOverSpringPartitionUpperLimitUsed > 0.0 &&
      this->forceEvaluator->isLinear()) {
    return this->network.evaluateStressTensor(
      this->currentDisplacements,
      this->forceEvaluator->getKappa(),
      this->is2D,
      this->lastSolveResult.oneOverSpringPartitionUpperLimitUsed);
  }
  const Eigen::VectorXd springVectors = this->getSpringVectors();
  return this->forceEvaluator->evaluateStressTensor(
    this->network, springVectors, this->network.vol, this->is2D);
}

double
MEHPSimulator::getPressure() const
{
  RUNTIME_EXP_IFN(this->forceEvaluator != nullptr,
                  "Force evaluator cannot be null.");
  if (this->network.nrOfSprings == 0) {
    return 0.0;
  }
  const Eigen::Matrix3d stress =
    (this->lastSolveResult.oneOverSpringPartitionUpperLimitUsed > 0.0 &&
     this->forceEvaluator->isLinear())
      ? this->network.evaluateStressTensor(
          this->currentDisplacements,
          this->forceEvaluator->getKappa(),
          this->is2D,
          this->lastSolveResult.oneOverSpringPartitionUpperLimitUsed)
      : this->forceEvaluator->evaluateStressTensor(
          this->network,
          this->getSpringVectors(),
          this->network.vol,
          this->is2D);
  return stress.trace() / 3.0;
}

Eigen::VectorXd
MEHPSimulator::getBondLengths()
{
  return this->getSpringLengths();
}

Eigen::VectorXd
MEHPSimulator::getCoordinates()
{
  return this->network.coordinates + this->currentDisplacements;
}

double
MEHPSimulator::getCurrentTime(const double /* currentStep */)
{
  return static_cast<double>(this->lastSolveResult.nrOfIterations);
}

double
MEHPSimulator::getTemperature()
{
  std::cerr << "Warning: Temperature is not a reasonable metric for this "
               "type of computation."
            << std::endl;
  return 0.0;
}

} // namespace pylimer_tools::sim::mehp
