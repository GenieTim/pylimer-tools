#ifndef MEHP_FORCE_RELAX2_H
#define MEHP_FORCE_RELAX2_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include "OutputSupportingSimulation.h"
#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <nlopt.hpp>
#include <string>
#include <vector>
#ifdef CEREALIZABLE
#include "../utils/CerealUtils.h"
#include <cereal/access.hpp>
#endif

namespace pylimer_tools::sim::mehp {
class MEHPForceRelaxation
  : public pylimer_tools::sim::OutputSupportingSimulation
{
private:
#ifdef CEREALIZABLE
  MEHPForceRelaxation() {}; // not exposed to users, only used by Cereal

  friend class cereal::access;
#endif

  // state
  pylimer_tools::entities::Universe universe;
  MEHPForceEvaluator* forceEvaluator;
  bool simulationHasRun = false;
  bool simulationSuggestsRerun = false;
  bool outputEndNodes = false;
  ExitReason exitReason = ExitReason::UNSET;
  int nrOfStepsDone = 0;
  Network forceRelaxationNetwork;
  Eigen::VectorXd currentSpringDistances;
  Eigen::VectorXd currentVelocities;
  Eigen::VectorXd currentVelocitiesPlus;
  Eigen::VectorXd currentForces;
  // config
  SimpleSpringMEHPForceEvaluator springForceEvaluator; // helper for memory time
  bool is2D = false;
  int defaultNrOfChains = 0;
  double defaultR0Squared = 0.0;
  double suggestRerunEps = 1e-3;
  int crossLinkerType = 2;
  double dt = 1;

public:
  MEHPForceRelaxation(const pylimer_tools::entities::Universe& u,
                      const int crossLinkerType = 2,
                      const bool is2D = false,
                      MEHPForceEvaluator* newForceEvaluator = nullptr,
                      const double kappa = 1.0,
                      const bool remove2functionalCrosslinkers = false,
                      const bool removeDanglingChains = false)
    : universe(u)
  {
    if (newForceEvaluator == nullptr) {
      this->springForceEvaluator = SimpleSpringMEHPForceEvaluator(kappa);
      newForceEvaluator = &this->springForceEvaluator;
    }
    this->crossLinkerType = crossLinkerType;
    // interpret network already to be able to give early results
    Network net;
    ConvertNetwork(&net,
                   crossLinkerType,
                   remove2functionalCrosslinkers,
                   removeDanglingChains);
    // this->defaultR0Squared =
    //   universe.computeMeanSquareEndToEndDistance(crossLinkerType);
    this->defaultNrOfChains = net.springsContourLength.size();
    if (this->defaultNrOfChains > 0) {
      this->defaultR0Squared = net.springsContourLength.mean() *
                               universe.computeMeanSquaredBondLength();
    }
    this->forceRelaxationNetwork = net;
    this->is2D = is2D;
    this->initializeDefaults();
    this->setForceEvaluator(newForceEvaluator);
  };

  MEHPForceRelaxation(const Network& net,
                      const bool is2D = false,
                      MEHPForceEvaluator* newForceEvaluator = nullptr,
                      const double kappa = 1.0)
  {
    this->forceRelaxationNetwork = net;
    this->is2D = is2D;

    if (newForceEvaluator == nullptr) {
      this->springForceEvaluator = SimpleSpringMEHPForceEvaluator(kappa);
      newForceEvaluator = &this->springForceEvaluator;
    }

    this->initializeDefaults();
    this->setForceEvaluator(newForceEvaluator);
  }

#ifdef CEREALIZABLE
  static MEHPForceRelaxation constructFromString(std::string s)
  {
    MEHPForceRelaxation res = MEHPForceRelaxation();
    pylimer_tools::utils::deserializeFromString(res, s);
    return res;
  }
#endif

  /**
   * @brief Actually do run the simulation
   *
   * @param algorithm
   * @param maxNrOfSteps
   * @param xtol
   * @param ftol
   */
  void runForceRelaxation(const char* algorithm = "LD_MMA",
                          long int maxNrOfSteps = 250000, // default: 10000
                          double xtol = 1e-12,
                          double ftol = 1e-9);

  /**
   * @brief Get the universe consisting of crosslinkers only
   *
   * @param newCrosslinkerType the type to give the crosslinkers
   * @return pylimer_tools::entities::Universe
   */
  pylimer_tools::entities::Universe getCrosslinkerVerse() const;

  /**
   * @brief Get the default R0 squared value
   *
   * @return double the default R_0^2
   */
  double getDefaultR0Square() const { return this->defaultR0Squared; }

  /**
   * @brief Get the volume of the simulation box
   *
   * @return double the volume
   */
  double getVolume() override { return this->forceRelaxationNetwork.vol; }

  /**
   * @brief Get the number of nodes in the network
   *
   * @return size_t the number of nodes
   */
  size_t getNrOfNodes() const { return this->forceRelaxationNetwork.nrOfNodes; }

  /**
   * @brief Get the number of springs in the network
   *
   * @return size_t the number of springs
   */
  size_t getNrOfSprings() const
  {
    return this->forceRelaxationNetwork.nrOfSprings;
  }

  /**
   * @brief Get the number of bonds (same as springs)
   *
   * @return size_t the number of bonds
   */
  size_t getNumBonds() override { return this->getNrOfSprings(); }

  /**
   * @brief Get the number of extra bonds
   *
   * @return size_t the number of extra bonds (always 0)
   */
  size_t getNumExtraBonds() override { return 0; }

  /**
   * @brief Get the number of bonds to form
   *
   * @return long int the number of bonds to form (always 0)
   */
  long int getNumBondsToForm() override { return 0; }

  /**
   * @brief Get the number of atoms (same as nodes)
   *
   * @return size_t the number of atoms
   */
  size_t getNumAtoms() override { return this->getNrOfNodes(); }

  /**
   * @brief Get the number of extra atoms
   *
   * @return size_t the number of extra atoms (always 0)
   */
  size_t getNumExtraAtoms() override { return 0; }

  /**
   * @brief Get the network structure
   *
   * @return Network the network structure
   */
  Network getNetwork() const { return this->forceRelaxationNetwork; }

  /**
   * @brief Configure assumption about box size
   *
   * @param assumption whether to assume the box is large enough
   */
  void configAssumeBoxLargeEnough(const bool assumption = true)
  {
    this->forceRelaxationNetwork.assumeBoxLargeEnough = assumption;
  }

  /**
   * @brief Configure whether the network is complete
   *
   * If true, the network is assumed to be complete,
   * i.e., no dangling or free chains have been omitted.
   *
   * @param assumeNetworkIsComplete whether to assume the network is complete
   */
  void configAssumeNetworkIsComplete(const bool assumeNetworkIsComplete = false)
  {
    this->forceRelaxationNetwork.assumeComplete = assumeNetworkIsComplete;
  }

  // MEHPForceEvaluator getForceEvaluator() const
  // {
  //   return *this->forceEvaluator;
  // }

  /**
   * @brief Set the force evaluator
   *
   * @param forceEvaluator pointer to the force evaluator to use
   */
  void setForceEvaluator(MEHPForceEvaluator* forceEvaluator)
  {
    this->forceEvaluator = forceEvaluator;
    this->forceEvaluator->setNetwork(this->forceRelaxationNetwork);
    this->forceEvaluator->setIs2D(this->is2D);
    this->forceEvaluator->prepareForEvaluations();
  }

  /**
   * @brief Get the Nr Of Active Nodes
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @param minimumNrOfActiveConnections the minimum number of active
   * connections required for a node to be considered active
   * @param maximumNrOfActiveConnections the maximum number of active
   * connections allowed (-1 for no limit)
   * @return int the number of active nodes
   */
  int getNrOfActiveNodes(const double tolerance = 1e-3,
                         const int minimumNrOfActiveConnections = 2,
                         const int maximumNrOfActiveConnections = -1) const
  {
    return this
      ->getIdsOfActiveNodes(
        tolerance, minimumNrOfActiveConnections, maximumNrOfActiveConnections)
      .size();
  }

  /**
   * @brief Get the Soluble Weight Fraction
   *
   * @param tolerance the tolerance for determining active springs
   * @return double the soluble weight fraction
   */
  double getSolubleWeightFraction(const double tolerance = 1e-3)
  {
    return this->computeSolubleWeightFraction(&this->forceRelaxationNetwork,
                                              tolerance);
  }

  /**
   * @brief Count the number of atoms that are in any way connected to an
   * active spring
   *
   * @param tolerance the tolerance for determining active springs
   * @return double the number of active clustered atoms
   */
  double countActiveClusteredAtoms(const double tolerance = 1e-3)
  {
    return this->countActiveClusteredAtoms(&this->forceRelaxationNetwork,
                                           tolerance);
  }

  /**
   * @brief Determine whether a node and spring, respectively, are in any
   * way connected to an active spring
   *
   * @param tolerance
   * @return std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
   */
  std::pair<Eigen::ArrayXb, Eigen::ArrayXb> findClusteredToActive(
    const double tolerance = 0.05) const
  {
    return this->findClusteredToActive(&this->forceRelaxationNetwork,
                                       tolerance);
  }

  /**
   * @brief Get the Dangling Weight Fraction
   *
   * @param tolerance the tolerance for determining active springs
   * @return double the dangling weight fraction
   */
  double getDanglingWeightFraction(const double tolerance = 1e-3)
  {
    return this->computeDanglingWeightFraction(&this->forceRelaxationNetwork,
                                               tolerance);
  }

  /**
   * @brief Get the crosslinker Chains that are active
   *
   * @param tolerance the tolerance for determining active springs
   * @return std::vector<pylimer_tools::entities::Molecule> vector of active
   * chains
   */
  std::vector<pylimer_tools::entities::Molecule> getActiveChains(
    const double tolerance = 1e-3) const
  {
    const std::vector<pylimer_tools::entities::Molecule> crossLinkerChains =
      this->universe.getChainsWithCrosslinker(crossLinkerType);
    std::vector<pylimer_tools::entities::Molecule> resultingChains;
    Eigen::ArrayXb springIsActive =
      this->findActiveSprings(&this->forceRelaxationNetwork, tolerance);
    for (size_t i = 0; i < crossLinkerChains.size(); ++i) {
      if (this->forceRelaxationNetwork.moleculeIdxToSpring[i] >= 0 &&
          springIsActive[this->forceRelaxationNetwork.moleculeIdxToSpring[i]]) {
        resultingChains.push_back(crossLinkerChains[i]);
      }
    }
    return resultingChains;
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
  std::unordered_map<long int, int> getEffectiveFunctionalityOfAtoms(
    double tolerance = 0.05) const;

  /**
   * @brief Get the Ids Of active Nodes
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @param minimumNrOfActiveConnections the number of active springs
   * required for this node to qualify as active
   * @return std::vector<long int> the atom ids
   */
  std::vector<long int> getIdsOfActiveNodes(
    double tolerance = 0.05,
    int minimumNrOfActiveConnections = 2,
    int maximumNrOfActiveConnections = -1) const;

  /**
   * @brief Get the Nr Of Active Springs connected to each node
   *
   * @param tolerance the tolerance: springs under a certain length are
   * considered inactive
   * @return Eigen::VectorXi
   */
  Eigen::VectorXi getNrOfActiveSpringsConnected(double tolerance = 0.05) const;

  /**
   * @brief Get the current spring distances
   *
   * @return Eigen::VectorXd the current spring distances
   */
  Eigen::VectorXd getCurrentSpringDistances() const
  {
    return this->currentSpringDistances;
  }

  /**
   * @brief Get the Nr Of Active Springs object
   *
   * @param tol the tolerance: springs under a certain length are considered
   * inactive
   * @return int the number of active springs
   */
  int getNrOfActiveSprings(const double tol = 0.1) const
  {
    return this->countNrOfActiveSprings(&this->forceRelaxationNetwork, tol);
  }

  /**
   * @brief Get the average contour length of all springs
   *
   * @return double the average contour length
   */
  double getAverageContourLength() const
  {
    return this->forceRelaxationNetwork.meanSpringContourLength;
  }

  /**
   * @brief Get the spring contour lengths
   *
   * @return Eigen::VectorXd the spring contour lengths
   */
  Eigen::VectorXd getSpringContourLength() const
  {
    return this->forceRelaxationNetwork.springsContourLength;
  }

  /**
   * @brief Get the Average Spring Length at the current step
   *
   * @return double the average spring length
   */
  double getAverageSpringLength() const;

  /**
   * @brief Get the Pressure
   *
   * @return double the pressure
   */
  double getPressure() const
  {
    return this->evaluatePressure(this->currentSpringDistances);
  }

  /**
   * @brief Get the Residual Norm at the current step
   *
   * @return double the residual norm
   */
  double getResidualNorm() const;

  /**
   * @brief Get the residual override
   *
   * @return double the residual
   */
  double getResidual() override
  {
    // std::cout << "ResidualNorm: " << this->getResidualNorm() <<
    // std::endl;
    return this->getResidualNorm();
  }

  /**
   * @brief Get the residuals (gradient) at the current step
   *
   * @return Eigen::VectorXd the residuals vector
   */
  Eigen::VectorXd getResiduals() const;

  /**
   * @brief Get the Force at the current step
   *
   * @return double the force
   */
  double getForce() const;

  /**
   * @brief Get the Gamma Factor at the current step
   *
   * @param b02 for the denominator, part of the melt <R_0^2> = b02 *
   * nrOfBondsInSpring
   * @param nrOfChains the nr of chains to average over (can be different
   * from the nr of springs thanks to omitted free chains or primary loops)
   * @return double
   */
  double getGammaFactor(double b02 = -1.0, int nrOfChains = -1) const;

  double getGamma() override { return this->getGammaFactor(-1.); }
  /**
   * @brief Get all the gamma factors for each spring
   *
   * @param b02 for the denominator, part of the melt <R_0^2> = b02 *
   * nrOfBondsInSpring
   * @return Eigen::VectorXd the gamma factors vector
   */
  Eigen::VectorXd getGammaFactors(double b2 = 1.) const;

  /**
   * @brief Get the number of iterations performed
   *
   * @return int the number of iterations
   */
  int getNrOfIterations() const { return this->nrOfStepsDone; }

  /**
   * @brief Get the exit reason of the simulation
   *
   * @return ExitReason the exit reason
   */
  ExitReason getExitReason() const { return this->exitReason; }

  /**
   * @brief Get the Spring Lengths
   *
   * @return Eigen::VectorXd the spring lengths vector
   */
  Eigen::VectorXd getSpringLengths() const
  {
    Eigen::VectorXd springDistances = this->getSpringDistances();
    Eigen::VectorXd springLengths = Eigen::VectorXd::Zero(
      static_cast<Eigen::Index>(this->forceRelaxationNetwork.nrOfSprings));
    for (size_t i = 0; i < this->forceRelaxationNetwork.nrOfSprings; ++i) {
      springLengths[static_cast<Eigen::Index>(i)] =
        springDistances.segment(3 * static_cast<Eigen::Index>(i), 3).norm();
    }
    return springLengths;
  }

  /**
   * @brief Get the Spring Distances
   *
   * @return Eigen::VectorXd the spring distances vector
   */
  Eigen::VectorXd getSpringDistances() const
  {
    return this->evaluateSpringDistances(&this->forceRelaxationNetwork,
                                         this->is2D);
  }

  /**
   * @brief Compute the spring lengths
   *
   * @param net the network to do the computation for
   * @param is2D whether to use 2D mode
   * @return Eigen::VectorXd the spring distances
   */
  static Eigen::VectorXd evaluateSpringDistances(const Network* net,
                                                 const bool is2D);

  /**
   * @brief Compute the spring lengths with displacement
   *
   * @param net the network to do the computation for
   * @param displacement the displacement vector
   * @param is2D whether to use 2D mode
   * @return Eigen::VectorXd the spring distances
   */
  static Eigen::VectorXd evaluateSpringDistances(
    const Network* net,
    const Eigen::VectorXd& displacement,
    const bool is2D);

  /**
   * @brief Return whether the simulation resulted in offsets close to the
   * limits
   *
   * @return true
   * @return false
   */
  bool suggestsRerun() const
  {
    return this->simulationSuggestsRerun || !this->simulationHasRun;
  }

  void configRerunEps(const double eps = 1e-3) { this->suggestRerunEps = eps; }

  /**
   * @brief Get the Stress Tensor
   *
   * @return Eigen::Matrix3d
   */
  Eigen::Matrix3d getStressTensor() override
  {
    const std::array<std::array<double, 3>, 3> stressTensor =
      this->evaluateStressTensor(this->currentSpringDistances,
                                 this->getVolume());
    Eigen::Matrix3d result = Eigen::Matrix3d::Zero();
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        result(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          stressTensor[i][j];
      }
    }
    return result;
  }

  double getTimestep() override { return this->dt; }

  double getCurrentTime(const double currentStep) override
  {
    return this->dt * currentStep;
  }

#ifdef CEREALIZABLE
  void writeRestartFile(std::string& filename) override
  {
    throw std::runtime_error("No restarts allowed here, yet");
  }
#endif

  int getNumShifts() override { return 0; }
  int getNumRelocations() override { return 0; }

  Eigen::VectorXd getBondLengths() override
  {
    Eigen::VectorXd lens =
      Eigen::VectorXd::Zero(this->currentSpringDistances.size() / 3);

    for (size_t i = 0; i < this->currentSpringDistances.size() / 3; ++i) {
      const double b = lens.segment(3 * static_cast<Eigen::Index>(i), 3).norm();
      lens[static_cast<Eigen::Index>(i)] = b;
    }

    return lens;
  }

  Eigen::VectorXd getCoordinates() override
  {
    return this->forceRelaxationNetwork.coordinates;
  }

  double getTemperature() override
  {
    std::cerr << "Warning: Temperature is not a reasonable metric for this "
                 "type of computation."
              << std::endl;
    return -1; // TODO: implement, since we have #runPhantomSteps?
  }

  size_t getNumParticles() override
  {
    return this->forceRelaxationNetwork.nrOfNodes;
  }

  /**
   * @brief Iterate all spring distances, mark active ones (length >
   * tolerance)
   *
   * @param tolerance
   * @return Eigen::ArrayXb
   */
  Eigen::ArrayXb findActiveSprings(const Network* net,
                                   const double tolerance = 1e-3) const;

protected:
  /**
   * @brief Convert the universe to a network
   *
   * @param net the target network
   * @param crossLinkerType the atom type of the crossLinker
   * @return true
   * @return false
   */
  bool ConvertNetwork(Network* net,
                      const int crossLinkerType = 2,
                      bool remove2functionalCrosslinkers = false,
                      bool removeDanglingChains = false);
  ;

  /**
   * @brief Compute the gamma factor from certain spring distances
   *
   * @param springDistances
   * @param b02 for the denominator, part of the melt <R_0^2> = b02 *
   * nrOfBondsInSpring
   * @param nrOfChains the nr of chains to average over (can be different
   * from the nr of springs thanks to omitted free chains or primary loops)
   * @return double
   */
  double evaluateGammaFactor(const Eigen::VectorXd& springDistances,
                             const double b02,
                             const int nrOfChains) const
  {
    return this->evaluateGammaFactors(springDistances, b02).sum() /
           static_cast<double>(nrOfChains);
  }

  Eigen::VectorXd evaluateGammaFactors(const Eigen::VectorXd& springDistances,
                                       const double b02) const
  {
    INVALIDARG_EXP_IFN(
      springDistances.size() ==
        this->forceRelaxationNetwork.springsContourLength.size() * 3,
      "Invalid sizes.");
    Eigen::VectorXd gammaFactors(springDistances.size() / 3);
    for (size_t i = 0; i < springDistances.size() / 3; ++i) {
      gammaFactors(i) =
        springDistances.segment(3 * i, 3).squaredNorm() /
        (this->forceRelaxationNetwork.springsContourLength(i) * b02);
    }
    return gammaFactors;
  }

  /**
   * @brief Evaluate the pressure of the network at specific spring
   * distances
   *
   * @param springDistances the spring distances
   * @return double
   */
  double evaluatePressure(const Eigen::VectorXd& springDistances) const
  {
    const auto stressTensor = this->evaluateStressTensor(
      springDistances, this->forceRelaxationNetwork.vol);
    return this->evaluatePressure(stressTensor);
  }

  /**
   * @brief Evaluate the pressure of the network at specific displacements
   *
   * @param net the network to evaluate the pressure for
   * @param u the displacements
   * @return double
   */
  double evaluatePressure(Network* net, const Eigen::VectorXd& u) const
  {
    const auto stressTensor = this->evaluateStressTensor(net, u, -1);
    return this->evaluatePressure(stressTensor);
  }

  /**
   * @brief Evaluate the pressure from the stress tensor
   *
   * @param stressTensor
   * @return double
   */
  static double evaluatePressure(
    const std::array<std::array<double, 3>, 3>& stressTensor)
  {
    return (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) / 3.0;
  }

  /**
   * @brief Compute the stress tensor
   *
   * @param net
   * @param u
   * @return std::array<std::array<double, 3>, 3>
   */
  std::array<std::array<double, 3>, 3> evaluateStressTensor(
    const Eigen::VectorXd& springDistances,
    const double volume) const;

  /**
   * @brief Compute the stress tensor
   *
   * @param net
   * @param u
   * @param loopTol
   * @return std::array<std::array<double, 3>, 3>
   */
  std::array<std::array<double, 3>, 3> evaluateStressTensor(
    Network* net,
    const Eigen::VectorXd& u,
    const double loopTol) const;

  /**
   * @brief Count how many of the springs are active (length > tolerance)
   *
   * @param springDistances
   * @param tolerance
   * @return int
   */
  int countNrOfActiveSprings(const Network* net,
                             const double tolerance = 1e-3) const
  {
    return (this->findActiveSprings(net, tolerance) == true).count();
  }

  /**
   * @brief Compute the elastically effective/active weight fraction of the
   * network
   *
   * @param net
   * @param tolerance
   * @return double
   */
  double computeActiveWeightFraction(Network* net,
                                     const double tolerance = 1e-3) const;

  /**
   * @brief Count the number of atoms that can be considered part of an
   * active cluster, i.e., are somehow connected to an active spring
   *
   * @param net
   * @param springDistances
   * @param tolerance
   * @return double
   */
  double countActiveClusteredAtoms(Network* net,
                                   const double tolerance = 1e-3) const;

  std::vector<int> getIndicesOfActiveNodes(const Network* net,
                                           const double tolerance = 1e-3) const;

  /**
   * @brief Find whether springs and nodes are in any way connected to an
   * active spring
   *
   * @param net the network that includes the connectivity
   * @param tolerance the tolerance for considering springs as active
   * @return std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
   */
  std::pair<Eigen::ArrayXb, Eigen::ArrayXb> findClusteredToActive(
    const Network* net,
    const double tolerance = 0.05) const
  {
    // find all active springs
    Eigen::ArrayXb activeSprings = this->findActiveSprings(net, tolerance);
    // then, iteratively walk along the springs to mark those as "active"
    // that are connected to active springs
    bool hadChanged = true;
    Eigen::ArrayXb nodeIsActive = Eigen::ArrayXb::Zero(net->nrOfNodes);
    while (hadChanged) {
      hadChanged = false;
      for (size_t i = 0; i < net->nrOfNodes; ++i) {
        if (nodeIsActive[i]) {
          continue; // already active
        }
        // check whether any of the springs connected to this node is active
        bool anyActive = false;
        for (const size_t spring_idx : net->springIndicesOfLinks[i]) {
          if (activeSprings[spring_idx]) {
            anyActive = true;
            break;
          }
        }
        // if so, mark this node as active as well as all the springs
        // connected to it
        if (anyActive) {
          hadChanged = true;
          nodeIsActive(i) = true;
          for (const size_t spring_idx : net->springIndicesOfLinks[i]) {
            activeSprings[spring_idx] = true;
          }
        }
      }
    }
    return std::make_pair(activeSprings, nodeIsActive);
  }

  /**
   * @brief Compute the weight fraction of non-active springs
   *
   * @param net
   * @param tolerance
   * @return double
   */
  double computeDanglingWeightFraction(Network* net,
                                       const double tolerance = 1e-3) const;

  /**
   * @brief Compute the weight fraction of springs connected to active
   * springs (any depth)
   *
   * @param net
   * @param tolerance
   * @return double
   */
  double computeSolubleWeightFraction(Network* net,
                                      const double tolerance = 1e-3) const
  {
    const double nActiveClusteredAtoms =
      this->countActiveClusteredAtoms(net, tolerance);

    double nAtoms = net->assumeComplete
                      ? static_cast<double>(net->nrOfNodes) +
                          net->springsContourLength.sum() - net->nrOfSprings
                      : static_cast<double>(this->universe.getNrOfAtoms());
    // finally, normalize by the number of atoms.
    // NOTE: currently, the weight of the atoms is ignored
    return 1. - (nActiveClusteredAtoms / nAtoms);
  }

  /**
   * @brief Decide whether a distance is within a given tolerance. This is
   * *the* distance criterion for determining whether a spring is active.
   *
   * @param dist
   * @param tolerance
   * @param contourLength

   * @return true
   * @return false
   */
  static bool distanceIsWithinTolerance(const Eigen::Vector3d& dist,
                                        const double tolerance = 1e-3,
                                        const double contourLength = 1.)
  {
    return dist.norm() <= (tolerance * std::max(contourLength, 1.));
  }

  void initializeDefaults()
  {
    this->currentForces =
      Eigen::VectorXd::Zero(this->forceRelaxationNetwork.coordinates.size());
    this->currentVelocities =
      Eigen::VectorXd::Zero(this->forceRelaxationNetwork.coordinates.size());
    this->currentVelocitiesPlus =
      Eigen::VectorXd::Zero(this->forceRelaxationNetwork.coordinates.size());
    this->currentSpringDistances =
      this->evaluateSpringDistances(&this->forceRelaxationNetwork, this->is2D);
  }
};
}
#endif
