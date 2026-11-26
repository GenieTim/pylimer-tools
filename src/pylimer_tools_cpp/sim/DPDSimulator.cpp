#include "./DPDSimulator.h"
#include "../utils/PerformanceTimer.h"
#include "../utils/StringUtils.h"

#include <fstream>
#include <iostream>
#include <random>
#include <unordered_set>

namespace pylimer_tools::sim::dpd {

DPDSimulator::DPDSimulator(const pylimer_tools::entities::Universe& u,
                           const int crossLinkerType,
                           const int slipspringBondType,
                           const bool is2D,
                           const std::string& seed,
                           const bool assumeBoxLargeEnough)
  : box(u.getBox())
  , universe(u)
  , neighbourlist(
      pylimer_tools::entities::EigenNeighbourList(Eigen::VectorXd(0),
                                                  this->box,
                                                  1.0))
{
  INVALIDARG_EXP_IFN(!is2D, "2D simulations are not supported yet.");
  this->is2D = is2D;
  this->assumeBoxLargeEnough = assumeBoxLargeEnough;
  this->reseedRandomness(seed);
  const double mean = 0.0;
  const double std = 1.0;
  const double a = mean - std::sqrt(3.) * std;
  const double b = mean + std::sqrt(3.) * std;
  this->uniform_rand_mean0std1 = std::uniform_real_distribution<double>(a, b);
  this->uniform_rand_between_0_1 =
    std::uniform_real_distribution<double>(0., 1.);

  // initialize the faster data structure
  this->box = u.getBox();
  this->coordinates = u.getUnwrappedVertexCoordinates(this->box);
  std::map<std::string, std::vector<long int>> edges = u.getEdges();
  this->bondPartnersA =
    Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(edges["edge_from"].size()));
  this->bondPartnersB =
    Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(edges["edge_to"].size()));
  this->bondTypes =
    Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(edges["edge_type"].size()));

  this->numBonds = 0;
  this->numSlipSprings = 0;
  // we have to "sort" the edges such that the slip-springs are at the end
  // so, first add "normal" bonds
  for (size_t i = 0; i < edges["edge_from"].size(); i++) {
    if (edges["edge_type"][i] != slipspringBondType) {
      this->bondPartnersA(this->numBonds) =
        static_cast<int>(edges["edge_from"][i]);
      this->bondPartnersB(this->numBonds) =
        static_cast<int>(edges["edge_to"][i]);
      this->bondTypes(this->numBonds) = static_cast<int>(edges["edge_type"][i]);
      this->numBonds += 1;
    }
  }
  // then, iterate again to add slip-springs
  for (size_t i = 0; i < edges["edge_from"].size(); i++) {
    if (edges["edge_type"][i] == slipspringBondType) {
      this->bondPartnersA(this->numBonds + this->numSlipSprings) =
        static_cast<int>(edges["edge_from"][i]);
      this->bondPartnersB(this->numBonds + this->numSlipSprings) =
        static_cast<int>(edges["edge_to"][i]);
      this->bondTypes(this->numBonds + this->numSlipSprings) =
        static_cast<int>(edges["edge_type"][i]);
      this->numSlipSprings += 1;
    }
  }

  this->neighbourlist =
    pylimer_tools::entities::EigenNeighbourList(coordinates, this->box, 1.0);
  this->numAtoms = static_cast<int>(this->coordinates.size() / 3);
  this->idxFunctionalities = Eigen::ArrayXi::Zero(this->numAtoms);
  this->atomTypes = u.getPropertyValues<int>("type");
  this->atomIds = u.getPropertyValues<long int>("id");
  this->maxBondLen = 0.45 * this->box.getL().maxCoeff();
  this->crossLinkerType = crossLinkerType;
  this->slipspringBondType = slipspringBondType;

  this->bondsOfIndex.reserve(static_cast<size_t>(this->numAtoms));
  for (size_t i = 0; i < static_cast<size_t>(this->numAtoms); ++i) {
    std::vector<size_t> bonds;
    bonds.reserve(4);
    this->bondsOfIndex.push_back(bonds);
  }
  this->bondPartnerCoordinatesA =
    Eigen::ArrayXi(3 * this->bondPartnersA.size());
  this->bondPartnerCoordinatesB =
    Eigen::ArrayXi(3 * this->bondPartnersB.size());
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(this->numBonds + this->numSlipSprings);
       ++i) {
    this->bondsOfIndex[static_cast<size_t>(this->bondPartnersA[i])].push_back(
      static_cast<size_t>(i));
    this->bondsOfIndex[static_cast<size_t>(this->bondPartnersB[i])].push_back(
      static_cast<size_t>(i));
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      this->bondPartnerCoordinatesA[i * 3 + dir] =
        this->bondPartnersA[i] * 3 + static_cast<int>(dir);
      this->bondPartnerCoordinatesB[i * 3 + dir] =
        this->bondPartnersB[i] * 3 + static_cast<int>(dir);
    }
  }
  this->isRelocationTarget = Eigen::ArrayXb::Constant(this->numAtoms, false);
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(this->numAtoms); ++i) {
    // if we already have slip-springs, the functionality must account for
    // that: The idxFunctionalities should not include the slip-springs
    const size_t atomIdx = static_cast<size_t>(i);
    this->idxFunctionalities[i] =
      this->numSlipSprings == 0
        ? static_cast<int>(this->bondsOfIndex[atomIdx].size())
        : std::accumulate(this->bondsOfIndex[atomIdx].begin(),
                          this->bondsOfIndex[atomIdx].end(),
                          0,
                          [&](const int val, const size_t bondIdx) {
                            return val + static_cast<int>(
                                           bondIdx <
                                           static_cast<size_t>(this->numBonds));
                          });
    if (this->idxFunctionalities[i] < 2) {
      this->chainEndIndices.push_back(atomIdx);
      this->isRelocationTarget[i] = true;
    }
  }

  this->resetBondOffsets();
  this->resetBondDuplicationPenalty();

  // simulation state
  this->currentVelocitiesPlus = Eigen::VectorXd::Zero(coordinates.size());
  this->currentVelocities = Eigen::VectorXd::Zero(coordinates.size());
  if (universe.vertexPropertyExists("vx") &&
      universe.vertexPropertyExists("vy") &&
      universe.vertexPropertyExists("vz")) {
    const std::vector<double> vx = universe.getPropertyValues<double>("vx");
    const std::vector<double> vy = universe.getPropertyValues<double>("vy");
    const std::vector<double> vz = universe.getPropertyValues<double>("vz");
    for (size_t i = 0; i < universe.getNrOfAtoms(); ++i) {
      const Eigen::Index coordBase = static_cast<Eigen::Index>(i * 3);
      this->currentVelocities(coordBase) = vx[i];
      this->currentVelocities(coordBase + 1) = vy[i];
      this->currentVelocities(coordBase + 2) = vz[i];
    }
  }

  this->currentForces = Eigen::VectorXd::Zero(coordinates.size());
  this->currentStressTensor = Eigen::Matrix3d::Zero();
  // with zero velocity, we don't have kinetic contributions
  // we want to compute forces such that the initial pressure etc. is
  // correct
  this->refreshCurrentState();
}

/**
 * @brief Run the main simulation loop with DPD and optionally Monte Carlo steps
 *
 * @param nSteps Number of simulation steps to run
 * @param withMC Whether to include Monte Carlo moves during simulation
 * @param shouldInterrupt Function to check if simulation should be interrupted
 * @param cleanupInterrupt Function to call when cleaning up after interruption
 */
void
DPDSimulator::runSimulation(const long int nSteps,
                            const bool withMC,
                            const std::function<bool()>& shouldInterrupt,
                            const std::function<void()>& cleanupInterrupt)
{
  bool wasInterrupted = false;

  const double halfDt = 0.5 * this->dt;
  const double lambdaDt = this->lambda * this->dt;
  double temperature = this->computeTemperature(this->currentVelocities);

  this->prepareAllOutputs();

  pylimer_tools::utils::PerformanceTimer timer =
    pylimer_tools::utils::PerformanceTimer<
      DPDPerformanceSections::NUM_PERFORMANCE_SECTIONS>();
  timer.registerSections(DPDPerformanceSectionNames);

  const pylimer_tools::entities::Box originalBox = this->box;

  // start iterating over the steps to take
  long int step = 0;
  for (; step < nSteps; ++step) {
    if (withMC && ((step % this->nStepsDPD) == 0)) {
      const double mcKbT =
        temperature > 0.0 ? temperature : this->getTargetKbT();
      this->numShifts = 0;
      timer.section(DPDPerformanceSections::RELOCATION);
      this->numRelocations = this->relocateSlipSprings(mcKbT);
      timer.section(DPDPerformanceSections::SHIFT);
      for (int i = 0; i < this->nStepsMC; ++i) {
        this->numShifts += this->shiftSlipSprings(mcKbT);
      }
    }
    // update coordinates & velocities
    timer.section(DPDPerformanceSections::TIME_STEPPING);
    this->currentVelocitiesPlus =
      this->currentVelocities + halfDt * this->currentForces;
    this->coordinates += dt * this->currentVelocitiesPlus;
    this->neighbourlist.resetCoordinates(this->coordinates);
    this->currentVelocities += lambdaDt * this->currentForces;

    // deformation
    if (this->doDeformation) {
      // adjust box
      pylimer_tools::entities::Box nextBox = originalBox.interpolate(
        this->deformationTargetBox,
        (1. + static_cast<double>(step)) / static_cast<double>(nSteps));
      this->deformBoxImmediately(nextBox);
    }

    // TODO: figure out, what the issue is, how the velocities
    // should be synchronized for the pressure
    // temperature = this->computeTemperature(velocities);

    // re-compute the forces with these updated coordinates & velocities
    timer.section(DPDPerformanceSections::FORCES);
    this->currentPressure = computeForces(this->currentForces,
                                          this->currentStressTensor,
                                          this->coordinates,
                                          this->currentVelocities,
                                          timer,
                                          this->dt,
                                          1.0);

    // correct the velocities
    timer.section(DPDPerformanceSections::TIME_STEPPING);
    this->currentVelocities =
      this->currentVelocitiesPlus + halfDt * this->currentForces;
    temperature = this->computeTemperature(this->currentVelocities);

    // kinetic term of the stress/pressure
    const double m_over_boxv = 1. / (this->box.getVolume());
    const auto velMat = this->currentVelocities.reshaped(
      3, static_cast<Eigen::Index>(this->numAtoms));
    this->currentStressTensor.noalias() +=
      m_over_boxv * (velMat * velMat.transpose());

    // bond formation
    if (this->bondsToForm > 0 && step % this->formBondsEvery == 0) {
      timer.section(DPDPerformanceSections::MODIFY);
      this->attemptBondFormation();
    }

    // output
    timer.section(DPDPerformanceSections::OUTPUT);
    this->currentStep += 1;
    this->currentTime += this->dt;
    this->handleOutput(this->currentStep);

    // allow Ctrl+C etc., without locking too much
    if (shouldInterrupt()) {
      wasInterrupted = true;
      break;
    }
  }

  // finish up
  this->closeAllOutputs();

  if (this->doDeformation) {
    this->deformBoxImmediately(this->deformationTargetBox);
    this->doDeformation = false;
  }

  if (wasInterrupted) {
    cleanupInterrupt();
  }

  timer.stop();
  timer.output();
}

/**
 * @brief re-calculate stress tensor & pressure
 *
 */
void
DPDSimulator::refreshCurrentState()
{
  pylimer_tools::utils::PerformanceTimer timer =
    pylimer_tools::utils::PerformanceTimer<
      DPDPerformanceSections::NUM_PERFORMANCE_SECTIONS>();
  this->currentPressure = computeForces(this->currentForces,
                                        this->currentStressTensor,
                                        this->coordinates,
                                        this->currentVelocities,
                                        timer,
                                        this->dt,
                                        1.0);
  // kinetic term of the stress/pressure
  const double m_over_boxv = 1. / (this->box.getVolume());
  const auto velMat = this->currentVelocities.reshaped(
    3, static_cast<Eigen::Index>(this->numAtoms));
  this->currentStressTensor.noalias() +=
    m_over_boxv * (velMat * velMat.transpose());
}

/**
 * @brief Set a new seed for the random generator
 *
 * @param seed
 */
void
DPDSimulator::reseedRandomness(const std::string& seed)
{
  this->currentSeedString = seed;
  // initialize the random number generator
  if (seed == "") {
    std::random_device rd;
    this->e2 = std::mt19937(rd());
  } else {
    std::seed_seq seed2(seed.begin(), seed.end());
    this->e2 = std::mt19937(seed2);
  }
  this->initThreadRNGs();
}

/**
 * @brief Sets the bond duplication penalty back to defaults.
 *
 * The bond duplication penalty is here, to enforce
 * "slip springs connecting two already bonded beads do not contribute in
 * the DPD steps"
 */
void
DPDSimulator::resetBondDuplicationPenalty()
{
  this->bondDuplicationPenalty =
    Eigen::ArrayXd::Constant(3 * (this->numBonds + this->numSlipSprings), 1.);
  // this->bondDuplicationPenalty.setConstant(1.);

  for (size_t i = 0; i < static_cast<size_t>(this->numAtoms); ++i) {
    this->resetBondDuplicationPenalty(i);
  }
}

void
DPDSimulator::resetBondDuplicationPenalty(const size_t atomIdx)
{
  std::unordered_set<size_t> partners;
  partners.reserve(this->bondsOfIndex[atomIdx].size());
  // here as well, we rely on the bonds of index being sorted
  for (const size_t bondIdx : this->bondsOfIndex[atomIdx]) {
    const Eigen::Index bondEigenIdx = static_cast<Eigen::Index>(bondIdx);
    const size_t atomPartnerIdx = static_cast<size_t>(
      this->bondPartnersA[bondEigenIdx] == static_cast<int>(atomIdx)
        ? this->bondPartnersB[bondEigenIdx]
        : this->bondPartnersA[bondEigenIdx]);
    if (bondIdx >= static_cast<size_t>(this->numBonds)) {
      // all others should stay 1 anyway
      this->bondDuplicationPenalty.segment(3 * bondEigenIdx, 3).setConstant(1.);
    }
    if (partners.contains(atomPartnerIdx)) {
      // "real" bonds always contribute -> check that this is a
      // slip-spring
      if (bondIdx >= static_cast<size_t>(this->numBonds)) {
        this->bondDuplicationPenalty.segment(3 * bondEigenIdx, 3).setZero();
      }
    } else {
      partners.insert(atomPartnerIdx);
    }
  }
}

/**
 * @brief Reset the offset required for the PBC bonds,
 *
 * This enables us having bonds longer than the box.
 * However, when calling this function, you reset that fact, i.e., only
 * call this if you are sure that you currently do not have bonds that
 * escape the box.
 *
 */
void
DPDSimulator::resetBondOffsets()
{
  this->bondBoxOffsets =
    this->box.getOffset(this->coordinates(this->bondPartnerCoordinatesB) -
                        this->coordinates(this->bondPartnerCoordinatesA));
}

void
DPDSimulator::resetBondOffset(const int bondIdx)
{
  this->bondBoxOffsets.segment(bondIdx * 3, 3) = this->box.getOffset(
    this->coordinates.segment(3 * this->bondPartnersB[bondIdx], 3) -
    this->coordinates.segment(3 * this->bondPartnersA[bondIdx], 3));
}

/**
 * @brief Try to form as many bonds as we can
 *
 */
void
DPDSimulator::attemptBondFormation()
{
  const double cutoff = this->bondFormationDistance;
  int bondsFormed = 0;
  // allocate possible neighbours
  Eigen::ArrayXi neighbors = Eigen::ArrayXi(static_cast<int>(
    this->numAtoms *
    (std::ceil((3.1 * cutoff) * (3.1 * cutoff) * (3.1 * cutoff)) /
     this->box.getVolume())));
  std::vector<int> possibleCandidates;
  // iterate atoms - we want to start from the crosslinkers
  for (size_t atom_idx = 0; atom_idx < static_cast<size_t>(this->numAtoms);
       ++atom_idx) {
    const Eigen::Index atomEigenIdx = static_cast<Eigen::Index>(atom_idx);
    possibleCandidates.clear();
    if (this->atomTypes[atom_idx] != this->atomTypeBondFormationFrom) {
      continue;
    }
    if (this->idxFunctionalities[atomEigenIdx] >=
        this->maxBondsPerType[this->atomTypes[atom_idx]]) {
      continue;
    }

    // find neighbours
    const int numNeighbors = this->neighbourlist.getIndicesCloseToCoordinates(
      neighbors, this->coordinates.segment(3 * atomEigenIdx, 3), cutoff);
    for (Eigen::Index neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
      // loop neighbours to find applicable partner
      const size_t j = static_cast<size_t>(neighbors[neigh_idx]);
      const Eigen::Index jEigen = static_cast<Eigen::Index>(j);
      Eigen::Vector3d vec = this->coordinates.segment(3 * atomEigenIdx, 3) -
                            this->coordinates.segment(3 * jEigen, 3);
      this->box.handlePBC(vec);
      const double r2 = (vec).norm();
      if (r2 <= cutoff && atom_idx != j &&
          (this->idxFunctionalities[jEigen] <
           this->maxBondsPerType[this->atomTypes[j]]) &&
          (this->atomTypeBondFormationTo == this->atomTypes[j])) {
        possibleCandidates.push_back(static_cast<int>(j));
      }
    }

    if (possibleCandidates.size() == 0) {
      continue;
    }

    // yay, we can form a bond
    // NOTE: currently, we build only 1
    std::ranges::shuffle(possibleCandidates, this->e2);
    this->addBond(
      static_cast<long int>(atom_idx),
      possibleCandidates[0],
      this->slipspringBondType == 3 ? 4 : 3); // TODO: decide on bond type
    this->bondsToForm -= 1;
    bondsFormed += 1;
    if (this->bondsToForm <= 0) {
      break;
    }
  }

  if (bondsFormed > 0) {
    this->validateState();
  }
}

/**
 * @brief Create a new bond between two nodes
 *
 * CAUTION: this is an expensive operation, involving resizing of Eigen
 * containers etc. Use sparsely.
 *
 * @param fromIdx
 * @param toIdx
 */
void
DPDSimulator::addBond(const long int fromIdx,
                      const long int toIdx,
                      const int bondType)
{
  // allocate space for the new bonds
  this->bondPartnerCoordinatesA.conservativeResize(
    this->bondPartnerCoordinatesA.size() + 3);
  this->bondPartnerCoordinatesB.conservativeResize(
    this->bondPartnerCoordinatesB.size() + 3);
  this->bondPartnersA.conservativeResize(this->bondPartnersA.size() + 1);
  this->bondPartnersB.conservativeResize(this->bondPartnersB.size() + 1);
  this->bondTypes.conservativeResize(this->bondTypes.size() + 1);
  this->bondBoxOffsets.conservativeResize(this->bondBoxOffsets.size() + 3);
  this->bondDuplicationPenalty.conservativeResize(
    this->bondDuplicationPenalty.size() + 3);
  assert(this->bondPartnersA.size() == this->bondPartnersB.size());
  assert(this->bondPartnersA.size() == this->bondTypes.size());
  assert(this->bondPartnersA.size() ==
         this->numBonds + this->numSlipSprings + 1);

  // move the data to make space (to keep slip-springs at the end of the
  // thing)
  const int newBondIdx = this->numBonds;
  for (int i = this->numBonds + this->numSlipSprings; i >= this->numBonds;
       --i) {
    this->bondPartnersA[i] = this->bondPartnersA[i - 1];
    this->bondPartnersB[i] = this->bondPartnersB[i - 1];
    this->bondTypes[i] = this->bondTypes[i - 1];

    for (int dir = 0; dir < 3; ++dir) {
      this->bondPartnerCoordinatesA[i * 3 + dir] =
        this->bondPartnerCoordinatesA[(i - 1) * 3 + dir];
      this->bondPartnerCoordinatesB[i * 3 + dir] =
        this->bondPartnerCoordinatesB[(i - 1) * 3 + dir];
      this->bondBoxOffsets[i * 3 + dir] =
        this->bondBoxOffsets[(i - 1) * 3 + dir];
      this->bondDuplicationPenalty[i * 3 + dir] =
        this->bondDuplicationPenalty[(i - 1) * 3 + dir];
    }
  }

  // also update all other references to the slip-springs
  for (size_t i = 0; i < static_cast<size_t>(this->numAtoms); ++i) {
    for (size_t j = 0; j < this->bondsOfIndex[i].size(); ++j) {
      if (this->bondsOfIndex[i][j] >= static_cast<size_t>(newBondIdx)) {
        this->bondsOfIndex[i][j] += 1;
      }
    }
  }

  // actually register the new bond
  this->bondTypes[newBondIdx] = bondType;
  this->bondPartnersA[newBondIdx] = static_cast<int>(fromIdx);
  this->bondPartnersB[newBondIdx] = static_cast<int>(toIdx);
  for (Eigen::Index i = 0; i < 3; ++i) {
    this->bondPartnerCoordinatesA[3 * newBondIdx + i] =
      3 * static_cast<int>(fromIdx) + static_cast<int>(i);
    this->bondPartnerCoordinatesB[3 * newBondIdx + i] =
      3 * static_cast<int>(toIdx) + static_cast<int>(i);
  }
  const size_t fromIdxU = static_cast<size_t>(fromIdx);
  const size_t toIdxU = static_cast<size_t>(toIdx);
  pylimer_tools::utils::addToSorted<size_t>(this->bondsOfIndex[fromIdxU],
                                            static_cast<size_t>(newBondIdx));
  pylimer_tools::utils::addToSorted<size_t>(this->bondsOfIndex[toIdxU],
                                            static_cast<size_t>(newBondIdx));

  this->resetBondOffset(newBondIdx);
  this->resetBondDuplicationPenalty(fromIdxU);
  this->resetBondDuplicationPenalty(toIdxU);

  // count the new bonds
  this->numBonds += 1;
  this->idxFunctionalities[fromIdx] += 1;
  this->idxFunctionalities[toIdx] += 1;
  // also adjust chain end indices
  if (!this->allowRelocationInNetwork) {
    this->isRelocationTarget[fromIdx] = (this->idxFunctionalities[fromIdx] < 2);
    this->isRelocationTarget[toIdx] = (this->idxFunctionalities[toIdx] < 2);
    // remove the atoms from chain ends, if they were in there
    if (!this->isRelocationTarget[fromIdx]) {
      pylimer_tools::utils::removeIfContained<size_t>(this->chainEndIndices,
                                                      fromIdxU);
    }
    if (!this->isRelocationTarget[toIdx]) {
      pylimer_tools::utils::removeIfContained<size_t>(this->chainEndIndices,
                                                      toIdxU);
    }
  }

#ifndef NDEBUG
  // std::cout << "Added bond " << newBondIdx << " between " << fromIdx
  //           << " and " << toIdx << ". Checking state..." << std::endl;
  this->validateState();
  // std::cout << "Added bond " << newBondIdx << ". State checked."
  //           << std::endl;
#endif
}

void
DPDSimulator::deformBoxImmediately(const pylimer_tools::entities::Box& newBox)
{
  const pylimer_tools::entities::Box previousBox = this->box;
  // rescale box offsets
  this->box = newBox;
  assert(newBox == this->box);
  const Eigen::Array3d scalingFactor = newBox.getL() / previousBox.getL();
  const Eigen::Index nBondsTotal =
    static_cast<Eigen::Index>(this->numBonds + this->numSlipSprings);
  const Eigen::Index nAtomsIdx = static_cast<Eigen::Index>(this->numAtoms);
  this->bondBoxOffsets.reshaped(3, nBondsTotal).array().colwise() *=
    scalingFactor;
  this->coordinates.reshaped(3, nAtomsIdx).array().colwise() *= scalingFactor;
  this->currentVelocities.reshaped(3, nAtomsIdx).array().colwise() *=
    scalingFactor;
  this->currentVelocitiesPlus.reshaped(3, nAtomsIdx).array().colwise() *=
    scalingFactor;
  this->currentForces.reshaped(3, nAtomsIdx).array().colwise() *= scalingFactor;

  // reset neighbour-list.
  // *maybe* it would be faster to implement a way to adjust the box, but
  // yeah...
  this->neighbourlist = pylimer_tools::entities::EigenNeighbourList(
    this->coordinates,
    this->box,
    1.0); //.resetCoordinates(this->coordinates);
}

/**
 * @brief Determine the temperature of the system
 *
 * @param velocities
 * @return double
 */
double
DPDSimulator::computeTemperature(const Eigen::VectorXd& velocities)
{
  if (velocities.size() == 0) {
    return 0.0;
  }
  // configuration
  constexpr double dim = 3;
  constexpr double kb = 1.0;
  constexpr double m = 1.0;

  const double KE = 0.5 * m * velocities.squaredNorm();

  return KE /
         ((dim / 2.) * (static_cast<double>(velocities.size()) / dim) * kb);
}

/**
 * @brief Compute the force vector, and return the pressure
 *
 * @param forces
 * @param stressTensor
 * @param coords
 * @param velocities
 * @param timer
 * @param dt
 * @param cutoff
 * @return double
 */
double
DPDSimulator::computeForces(
  Eigen::VectorXd& forces,
  Eigen::Matrix3d& stressTensor,
  const Eigen::VectorXd& coords,
  const Eigen::VectorXd& velocities,
  pylimer_tools::utils::PerformanceTimer<
    DPDPerformanceSections::NUM_PERFORMANCE_SECTIONS>& timer,
  const double timeStep,
  const double cutoff)
{
  // initialisation
  assert(coordinates.size() == velocities.size());
  assert(forces.size() == coords.size());
  double pressure = 0.0;
  forces.setZero();
  stressTensor.setZero();

  // actual computation
  // (attractive) bond forces
  // TODO: investigate whether we would be faster to just push the indices
  // of the bonds to ignore
  timer.section(DPDPerformanceSections::BOND_FORCE);
  Eigen::VectorXd bondDistances =
    ((coords(this->bondPartnerCoordinatesB) -
      coords(this->bondPartnerCoordinatesA) + this->bondBoxOffsets)
       .array() *
     this->bondDuplicationPenalty)
      .matrix();
  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(bondDistances);
#ifndef NDEBUG
    const double halfMinBoxL = 0.5 * this->box.getL().minCoeff();
    const double halfMinBoxLSq = halfMinBoxL * halfMinBoxL;
    for (Eigen::Index i = 0; i < this->bondPartnersA.size(); ++i) {
      const double bondLenSq = bondDistances.segment(3 * i, 3).squaredNorm();
      if (bondLenSq > halfMinBoxLSq) {
        std::cerr
          << "WARNING: Bond " << i << " has length " << std::sqrt(bondLenSq)
          << ", which violates the assumption that the box is large enough"
          << std::endl;
      }
    }
#endif
  } else {
#ifndef NDEBUG
    for (Eigen::Index i = 0; i < this->bondPartnersA.size(); ++i) {
      if (this->bondBoxOffsets.segment(3 * i, 3) !=
          this->box.getOffset(
            coords(this->bondPartnerCoordinatesB.segment(3 * i, 3)) -
            coords(this->bondPartnerCoordinatesA.segment(3 * i, 3)))) {
        std::cerr << "INFO: Bond " << i << " spans more than one image."
                  << std::endl;
      }
    }
#endif
  }
  // assert(bondDistances.minCoeff() > -this->box.getL().maxCoeff());
  // assert(bondDistances.maxCoeff() < this->box.getL().maxCoeff());
  forces(this->bondPartnerCoordinatesA) += this->k * bondDistances;
  forces(this->bondPartnerCoordinatesB) -= this->k * bondDistances;

  // we use our own parallelization -> disable the one by Eigen
  Eigen::setNbThreads(1);
  // attractive force -> reduces pressure in the system
  assert(bondDistances.allFinite());
  pressure -= this->k * bondDistances.squaredNorm();
  const auto bondDistancesMat =
    bondDistances.reshaped(3, this->bondPartnersA.size());
  stressTensor.noalias() -=
    this->k * (bondDistancesMat * bondDistancesMat.transpose());

  timer.section(DPDPerformanceSections::PAIR_FORCE);

  // actually loop the atoms
#pragma omp parallel num_threads(this->getNumOmpThreads())
  {
#ifdef OPENMP_FOUND
    const int tid = omp_get_thread_num();
    const int numThreads = omp_get_num_threads();
#else
    const int tid = 0;
    const int numThreads = 1;
#endif
    std::mt19937& localRng = (numThreads == 1 || this->threadRNGs.empty())
                               ? this->e2
                               : this->threadRNGs[static_cast<size_t>(tid)];
    std::uniform_real_distribution<double> localDist =
      this->uniform_rand_mean0std1;

    const double sigmaOverSqrtDt = this->sigma / std::sqrt(timeStep);
    // pre-allocate the neighbor indices array
    Eigen::ArrayXi neighbors = Eigen::ArrayXi(static_cast<int>(
      this->numAtoms *
      (std::ceil((3.1 * cutoff) * (3.1 * cutoff) * (3.1 * cutoff)) /
       this->box.getVolume())));
    Eigen::Vector3d pairdistance;
    Eigen::Vector3d pairdistanceNormed;
    Eigen::Vector3d velocitydiff;
    Eigen::Vector3d pairForce;
    Eigen::Vector3d zeroOneTwo;
    zeroOneTwo << 0, 1, 2;

    // need to fix the schedule as with higher i, the workload gets much
    // less
#pragma omp for reduction(+ : forces, stressTensor, pressure)                  \
  schedule(dynamic, 16)
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(this->numAtoms - 1);
         ++i) {
      const int numNeighbors = this->neighbourlist.getHigherIndicesWithinCutoff(
        neighbors, coords, static_cast<int>(i), cutoff);

      // pair forces
      // Eigen::ArrayXi neighbourCoordIndices =
      // (neighbors.head(numNeighbors) * 3)
      //            .replicate(1, 3)
      //            .transpose()
      //            .eval()
      //            .reshaped(3 * numNeighbors, 1) +
      //          zeroOneTwo.replicate(numNeighbors, 1);
      // Eigen::VectorXd pairDistances =
      //   coords.segment(3 * i, 3).replicate(numNeighbors, 1) -
      //   coords(neighbourCoordIndices);
      // this->box.handlePBC(pairDistances);
      // Eigen::VectorXd rNorms = pairDistances.reshaped(numNeighbors,
      // 3).rowwise().norm(); assert(rNorms.size() == numNeighbors);
      // Eigen::VectorXd pairDistancesNormed = pairDistances /
      // rNorms.replicate(1, 3).transpose().eval().reshaped(3*numNeighbors,
      // 1); Eigen::VectorXd velocityDiffs = velocities.segment(3 * i,
      // 3).replicate(numNeighbors, 1) - velocities(neighbourCoordIndices);

      // pair forces
      for (Eigen::Index neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
        const Eigen::Index j = static_cast<Eigen::Index>(neighbors[neigh_idx]);
        pairdistance = coords.segment<3>(3 * i) - coords.segment<3>(3 * j);
        this->box.handlePBC(pairdistance);
        // slight performance improvement, taking the square norm here
        const double rNorm2 = pairdistance.squaredNorm();
        if (rNorm2 == 0.) {
          std::cerr << "WARNING: zero distance between atoms " << i << " and "
                    << j << "." << std::endl;
          continue;
        }

        const double rNorm = std::sqrt(rNorm2);
        const double one_minus_rnorm = 1. - rNorm;
        const double one_minus_rnorm2 = (1. - rNorm) * (1. - rNorm);
        pairdistanceNormed = pairdistance / rNorm;

        // dissipative/drag force
        velocitydiff =
          velocities.segment<3>(3 * i) - velocities.segment<3>(3 * j);
        const double rij_dot_vij = pairdistanceNormed.dot(velocitydiff);
        const double gamma_weighted_rij_dot_vij =
          this->gamma * one_minus_rnorm2 * rij_dot_vij;

        // conservative repulsion force - dissipative/drag force
        // TODO: check if we get fma here
        double pairForceConst =
          this->A * one_minus_rnorm - gamma_weighted_rij_dot_vij;

        // random force
        const double constant_rnd_prefix = one_minus_rnorm * sigmaOverSqrtDt;
        const double random_val = localDist(localRng);

        pairForceConst += constant_rnd_prefix * random_val;

        // back to 3D/Eigen Vector
        pairForce = pairForceConst * pairdistanceNormed;

        // actually assign the new forces
        forces.segment<3>(3 * i) += pairForce;
        forces.segment<3>(3 * j) -= pairForce;

        // pressure update: repulsive force -> increases pressure
        pressure += pairForce.dot(pairdistance);
        assert(std::isfinite(pressure));
        stressTensor.noalias() += pairForce * pairdistance.transpose();
      }
    }
  }

  // reset Eigen threads
  // Eigen::setNbThreads(0);
  pressure /= (3. * this->box.getVolume());
  stressTensor /= this->box.getVolume();
#ifndef NDEBUG
  assert(APPROX_REL_EQUAL(stressTensor.trace() / 3., pressure, 1e-5));
#endif
  return pressure;
}

/**
 * @brief Register a set of atoms and time for being measured for msd
 *
 * @param atomIdsToMeasure
 */
void
DPDSimulator::startMeasuringMSDForAtoms(
  const std::vector<size_t>& atomIdsToMeasure)
{
  // Translate atom IDS to indices of the local structure
  Eigen::ArrayXi coordinateIndices =
    Eigen::ArrayXi(static_cast<Eigen::Index>(3 * atomIdsToMeasure.size()));
#pragma omp parallel for if (atomIdsToMeasure.size() > 500)
  for (Eigen::Index i = 0;
       i < static_cast<Eigen::Index>(atomIdsToMeasure.size());
       ++i) {
    const size_t atomId = atomIdsToMeasure[static_cast<size_t>(i)];
    const int index = static_cast<int>(
      this->universe.getIdxByAtomId(static_cast<long int>(atomId)));
    coordinateIndices[3 * i] = 3 * index;
    coordinateIndices[3 * i + 1] = 3 * index + 1;
    coordinateIndices[3 * i + 2] = 3 * index + 2;
  }

  // Remember to measure these relative to the current time-step
  msdMeasuredIndices.push_back(coordinateIndices);
  msdOrigins.push_back(this->coordinates(coordinateIndices));
  msdOriginTimesteps.push_back(static_cast<size_t>(this->currentStep));
}

bool
DPDSimulator::acceptMetropolisMove(const double deltaEnergy, const double kbT)
{
  if (deltaEnergy <= 0.0) {
    return true;
  }
  if (kbT <= 0.0) {
    return false;
  }
  const double factor = std::exp(-deltaEnergy / kbT);
  return this->uniform_rand_between_0_1(this->e2) < factor;
}

int
DPDSimulator::collectSlipSpringPartnerCandidates(
  const size_t centerAtom,
  Eigen::ArrayXi& neighbours,
  std::vector<size_t>& candidates)
{
  const Eigen::Index centerEigen = static_cast<Eigen::Index>(centerAtom);
  const double lowCutoffSq = this->lowCutoff * this->lowCutoff;
  const double highCutoffSq = this->highCutoff * this->highCutoff;
  int numCandidates = 0;
  const int numNeighs = this->neighbourlist.getIndicesCloseToCoordinates(
    neighbours,
    this->coordinates.segment(3 * centerEigen, 3),
    this->highCutoff);
  RUNTIME_EXP_IFN(numNeighs <= neighbours.size(),
                  "Neighbourlist does not act as it should.");
  for (Eigen::Index j = 0; j < numNeighs; ++j) {
    RUNTIME_EXP_IFN(neighbours[j] >= 0 && neighbours[j] < this->numAtoms,
                    "Neighbourlist seems inconsequent.");
    const size_t neighAtom = static_cast<size_t>(neighbours[j]);
    const Eigen::Index neighEigen = static_cast<Eigen::Index>(neighAtom);
    if (neighAtom == centerAtom ||
        this->atomTypes[neighAtom] == this->crossLinkerType ||
        this->idxFunctionalities[neighEigen] > 2) {
      continue;
    }
    Eigen::Vector3d distance = this->coordinates.segment(3 * centerEigen, 3) -
                               this->coordinates.segment(3 * neighEigen, 3);
    this->box.handlePBC(distance);
    const double distSq = distance.squaredNorm();
    if (distSq > lowCutoffSq && distSq <= highCutoffSq) {
      if (numCandidates < static_cast<int>(candidates.size())) {
        candidates[static_cast<size_t>(numCandidates)] = neighAtom;
      } else {
        candidates.push_back(neighAtom);
      }
      numCandidates += 1;
    }
  }
  return numCandidates;
}

/**
 * @brief Randomly add new slip-springs
 *
 * @param num
 * @param bondType
 * @return int
 */
int
DPDSimulator::createSlipSprings(const int num, const int bondType)
{
  INVALIDARG_EXP_IFN(num >= 0,
                     "Number of slip-springs to create must be non-negative.");
  INVALIDARG_EXP_IFN(bondType >= 0, "Bond type must be non-negative.");
  if (num == 0 || this->numAtoms == 0) {
    return 0;
  }
  const int effectiveBondType =
    bondType == 0 ? this->slipspringBondType : bondType;
  int createdLastIteration = 100;
  int totalCreated = 0;
  std::vector<size_t> candidates;
  candidates.reserve(5);

  std::vector<size_t> slipSpringFrom;
  slipSpringFrom.reserve(static_cast<size_t>(num));
  std::vector<size_t> slipSpringTo;
  slipSpringTo.reserve(static_cast<size_t>(num));

  // randomly permute the atoms to start the search with
  // we do this over just randomly selecting an index
  // in order to reduce the probability of finding the same start twice
  std::vector<size_t> sourceIds;
  sourceIds.reserve(static_cast<size_t>(this->numAtoms));
  for (size_t i = 0; i < static_cast<size_t>(this->numAtoms); ++i) {
    if (this->atomTypes[i] != this->crossLinkerType &&
        this->idxFunctionalities[static_cast<Eigen::Index>(i)] <= 2) {
      sourceIds.push_back(i);
    }
  }
  if (sourceIds.empty()) {
    return 0;
  }

  // search for neighbours that are eligible
  Eigen::ArrayXi neighbours = Eigen::ArrayXi(16);
  while (createdLastIteration > 0 && totalCreated < num) {
    createdLastIteration = 0;
    std::ranges::shuffle(sourceIds, this->e2);
    for (const size_t i : sourceIds) {
      const int numCandidates =
        this->collectSlipSpringPartnerCandidates(i, neighbours, candidates);
      if (numCandidates == 0) {
        continue;
      }
      std::uniform_int_distribution<int> dist(0, numCandidates - 1);
      const int candidateIndex = dist(this->e2);
      // found a candidate to create a slip-spring to
      slipSpringFrom.push_back(i);
      slipSpringTo.push_back(candidates[static_cast<size_t>(candidateIndex)]);
      totalCreated += 1;
      createdLastIteration += 1;
      if (totalCreated >= num) {
        break;
      }
    }
  }

  this->addSlipSprings(slipSpringFrom, slipSpringTo, effectiveBondType);
  return totalCreated;
}

/**
 * @brief Add slip-springs according to the specification
 *
 * @param partnerA
 * @param partnerB
 * @param bondType
 */
void
DPDSimulator::addSlipSprings(std::vector<size_t>& partnerA,
                             std::vector<size_t>& partnerB,
                             const int bondType)
{
  INVALIDARG_EXP_IFN(partnerA.size() == partnerB.size(),
                     "Require same size A & B");
  for (size_t i = 0; i < partnerA.size(); ++i) {
    INVALIDINDEX_EXP_IFN(partnerA[i] < static_cast<size_t>(this->numAtoms),
                         "Invalid partner id");
    INVALIDINDEX_EXP_IFN(partnerB[i] < static_cast<size_t>(this->numAtoms),
                         "Invalid partner id");
  }
  const Eigen::Index sizeBefore =
    static_cast<Eigen::Index>(this->numBonds + this->numSlipSprings);
  const Eigen::Index numNew = static_cast<Eigen::Index>(partnerA.size());
  this->bondPartnersA.conservativeResize(sizeBefore + numNew);
  this->bondPartnersB.conservativeResize(sizeBefore + numNew);
  this->bondPartnerCoordinatesA.conservativeResize(3 * (sizeBefore + numNew));
  this->bondPartnerCoordinatesB.conservativeResize(3 * (sizeBefore + numNew));
  this->bondTypes.conservativeResize(sizeBefore + numNew);
  this->bondBoxOffsets.conservativeResize(3 * (sizeBefore + numNew));
  this->bondDuplicationPenalty.conservativeResize(3 * (sizeBefore + numNew));

  this->bondPartnersA.segment(sizeBefore, numNew) =
    Eigen::Map<Eigen::ArrayXst, Eigen::Unaligned>(partnerA.data(), numNew)
      .cast<int>();
  this->bondPartnersB.segment(sizeBefore, numNew) =
    Eigen::Map<Eigen::ArrayXst, Eigen::Unaligned>(partnerB.data(), numNew)
      .cast<int>();
  this->bondTypes.segment(sizeBefore, numNew) = bondType;

  for (Eigen::Index i = sizeBefore; i < sizeBefore + numNew; ++i) {
    this->bondsOfIndex[static_cast<size_t>(this->bondPartnersA[i])].push_back(
      static_cast<size_t>(i));
    this->bondsOfIndex[static_cast<size_t>(this->bondPartnersB[i])].push_back(
      static_cast<size_t>(i));
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      this->bondPartnerCoordinatesA[i * 3 + dir] =
        this->bondPartnersA[i] * 3 + static_cast<int>(dir);
      this->bondPartnerCoordinatesB[i * 3 + dir] =
        this->bondPartnersB[i] * 3 + static_cast<int>(dir);
    }
    this->resetBondOffset(static_cast<int>(i));
    this->resetBondDuplicationPenalty(
      static_cast<size_t>(this->bondPartnersA[i]));
    this->resetBondDuplicationPenalty(
      static_cast<size_t>(this->bondPartnersB[i]));
  }

  this->numSlipSprings += static_cast<int>(partnerA.size());
}

/**
 * @brief Move slip-springs that are at ends to new ends
 *
 * @param kbT
 * @return int
 */
int
DPDSimulator::relocateSlipSprings(const double kbT)
{
  INVALIDARG_EXP_IFN(kbT >= 0.0, "Thermal energy kbT must be non-negative.");
  if (this->chainEndIndices.empty() || this->numSlipSprings == 0) {
    return 0;
  }
  int nAccept = 0;
  std::vector<size_t> candidates;
  Eigen::ArrayXi neighbours = Eigen::ArrayXi(12);
  std::uniform_int_distribution<size_t> chainendDist(
    0, this->chainEndIndices.size() - 1);
  for (size_t springIdx = static_cast<size_t>(this->numBonds);
       springIdx < static_cast<size_t>(this->numBonds + this->numSlipSprings);
       ++springIdx) {
    const Eigen::Index springEigenIdx = static_cast<Eigen::Index>(springIdx);
    const size_t partnerA =
      static_cast<size_t>(this->bondPartnersA[springEigenIdx]);
    const size_t partnerB =
      static_cast<size_t>(this->bondPartnersB[springEigenIdx]);
    if ((!(this->isRelocationTarget[static_cast<Eigen::Index>(partnerA)] ||
           this->isRelocationTarget[static_cast<Eigen::Index>(partnerB)])) &&
        // always relocate from crosslinkers away, e.g. in case the
        // crosslink has gained an additional bond and is no longer a
        // relocation target.
        // consequently, we are inconsistent for crosslinkers with f = 2
        !(this->atomTypes[partnerA] == this->crossLinkerType ||
          this->atomTypes[partnerB] == this->crossLinkerType)) {
      continue;
    }

    const size_t candidateIndex = this->chainEndIndices[chainendDist(this->e2)];
    const Eigen::Index candidateEigenIdx =
      static_cast<Eigen::Index>(candidateIndex);
    if (this->atomTypes[candidateIndex] == this->crossLinkerType ||
        this->idxFunctionalities[candidateEigenIdx] > 2) {
      continue;
    }

    const int numCandidates = this->collectSlipSpringPartnerCandidates(
      candidateIndex, neighbours, candidates);
    if (numCandidates == 0) {
      continue;
    }
    std::uniform_int_distribution<int> candidateDist(0, numCandidates - 1);
    const size_t candidatePartnerIndex =
      candidates[static_cast<size_t>(candidateDist(this->e2))];

    // compute the Metropolis criterion using the unwrapped PBC bond vector
    const Eigen::Vector3d bondDistanceNow =
      this->computeBondDistance(static_cast<int>(springIdx));
    const double bondEnergyNow = 0.5 * this->k * bondDistanceNow.squaredNorm();
    // for the new slip-spring, we cannot rely on bond box offsets:
    // it is a freshly formed minimum-image slip-spring
    Eigen::Vector3d bondDistanceNew =
      this->coordinates.segment(
        static_cast<Eigen::Index>(candidatePartnerIndex * 3), 3) -
      this->coordinates.segment(candidateEigenIdx * 3, 3);
    this->box.handlePBC(bondDistanceNew);
    const double bondEnergyNew = 0.5 * this->k * bondDistanceNew.squaredNorm();
    if (this->acceptMetropolisMove(bondEnergyNew - bondEnergyNow, kbT)) {
      this->setSlipSpringPartners(
        springIdx, candidateIndex, candidatePartnerIndex, nullptr);
      nAccept++;
      RUNTIME_EXP_IFN(
        this->computeBondLength(static_cast<int>(springIdx)) <=
          this->highCutoff + 1e-10,
        "By relocation, newly created bonds should not be too long.");
    }
  }
  return nAccept;
}

/**
 * @brief Perform Monte Carlo slip-spring shift moves along polymer chains
 *
 * @param kbT
 * @return int
 */
int
DPDSimulator::shiftSlipSprings(const double kbT)
{
  INVALIDARG_EXP_IFN(kbT >= 0.0, "Thermal energy kbT must be non-negative.");
  if (this->numSlipSprings == 0 || this->numAtoms == 0) {
    return 0;
  }
  int n_accept = 0;
  std::vector<size_t> candidates;
  std::vector<size_t> nonCrosslinkerIndices;
  Eigen::ArrayXi neighbours = Eigen::ArrayXi(12);

  for (size_t springIdx = static_cast<size_t>(this->numBonds);
       springIdx < static_cast<size_t>(this->numBonds + this->numSlipSprings);
       ++springIdx) {
    const Eigen::Index springEigenIdx = static_cast<Eigen::Index>(springIdx);
    const size_t prevPartnerA =
      static_cast<size_t>(this->bondPartnersA[springEigenIdx]);
    const size_t prevPartnerB =
      static_cast<size_t>(this->bondPartnersB[springEigenIdx]);
    const Eigen::Vector3d prevOffset =
      this->bondBoxOffsets.segment(3 * springEigenIdx, 3);

    int springAccepts = 0;
    if (this->shiftOneAtATime) {
      springAccepts += static_cast<int>(
        this->attemptSlipSpringShift(springIdx, prevPartnerA, kbT));
      if (this->bondPartnersA[springEigenIdx] !=
          this->bondPartnersB[springEigenIdx]) {
        const size_t currentPartnerB =
          static_cast<size_t>(this->bondPartnersB[springEigenIdx]);
        springAccepts += static_cast<int>(
          this->attemptSlipSpringShift(springIdx, currentPartnerB, kbT));
      }
    } else {
      springAccepts +=
        static_cast<int>(this->attemptSlipSpringShift(springIdx, kbT));
    }

    if (this->bondPartnersA[springEigenIdx] ==
        this->bondPartnersB[springEigenIdx]) {
      // Spring collapsed onto a single bead: relocate to a random eligible pair
      if (nonCrosslinkerIndices.empty()) {
        nonCrosslinkerIndices.reserve(static_cast<size_t>(this->numAtoms));
        for (size_t idx = 0; idx < static_cast<size_t>(this->numAtoms); ++idx) {
          if (this->atomTypes[idx] != this->crossLinkerType &&
              this->idxFunctionalities[static_cast<Eigen::Index>(idx)] <= 2) {
            nonCrosslinkerIndices.push_back(idx);
          }
        }
      }
      std::ranges::shuffle(nonCrosslinkerIndices, this->e2);
      bool relocated = false;
      for (const size_t firstPartner : nonCrosslinkerIndices) {
        const int numCandidates = this->collectSlipSpringPartnerCandidates(
          firstPartner, neighbours, candidates);
        if (numCandidates == 0) {
          continue;
        }
        std::uniform_int_distribution<int> candidateDist(0, numCandidates - 1);
        const size_t secondPartner =
          candidates[static_cast<size_t>(candidateDist(this->e2))];
        this->setSlipSpringPartners(
          springIdx, firstPartner, secondPartner, nullptr);
        RUNTIME_EXP_IFN(
          this->computeBondLength(static_cast<int>(springIdx)) <=
            this->highCutoff + 1e-10,
          "By shifting, newly created bonds should not be too long.");
        relocated = true;
        break;
      }
      if (!relocated) {
        // No eligible pair exists within (lowCutoff, highCutoff]: revert move
        // so the slip-spring never remains as an invalid self-loop.
        this->setSlipSpringPartners(
          springIdx, prevPartnerA, prevPartnerB, &prevOffset);
        springAccepts = 0;
      }
    }
    n_accept += springAccepts;
  }
  return n_accept;
}

bool
DPDSimulator::proposeRailStep(const size_t partner,
                              size_t& newPartner,
                              Eigen::Vector3d& deltaR,
                              Eigen::Vector3d& deltaOffset)
{
  const Eigen::Index partnerEigen = static_cast<Eigen::Index>(partner);
  int distrLimit = this->idxFunctionalities[partnerEigen] - 1;
  if (distrLimit < 0) {
    RUNTIME_EXP_IFN(
      this->atomTypes[partner] == this->crossLinkerType,
      "Only crosslinkers are allowed to be single beads. Found bead " +
        std::to_string(partner) + " to have type " +
        std::to_string(this->atomTypes[partner]) + " with functionality " +
        std::to_string(distrLimit) + ".");
    return false;
  }
  if (distrLimit == 0 && this->shiftPossibilityEmpty) {
    distrLimit += 1;
  }
  std::uniform_int_distribution<int> dist(0, distrLimit);
  const int randomIdx = dist(this->e2);
  if (randomIdx >= this->idxFunctionalities[partnerEigen]) {
    RUNTIME_EXP_IFN(this->shiftPossibilityEmpty, "Invalid state.");
    return false;
  }
  const size_t selectedRailBond =
    this->bondsOfIndex[partner][static_cast<size_t>(randomIdx)];
  const Eigen::Index railBondEigen =
    static_cast<Eigen::Index>(selectedRailBond);
  assert(this->bondPartnersA[railBondEigen] == static_cast<int>(partner) ||
         this->bondPartnersB[railBondEigen] == static_cast<int>(partner));
  const bool shiftEndIsFirstOnRailBond =
    this->bondPartnersA[railBondEigen] == static_cast<int>(partner);
  newPartner = static_cast<size_t>(shiftEndIsFirstOnRailBond
                                     ? this->bondPartnersB[railBondEigen]
                                     : this->bondPartnersA[railBondEigen]);
  deltaOffset = (shiftEndIsFirstOnRailBond ? 1.0 : -1.0) *
                this->bondBoxOffsets.segment(3 * railBondEigen, 3);
  deltaR =
    (this->coordinates.segment(static_cast<Eigen::Index>(newPartner) * 3, 3) -
     this->coordinates.segment(partnerEigen * 3, 3)) +
    deltaOffset;
  return true;
}

bool
DPDSimulator::attemptSlipSpringShift(const size_t springIdx, const double kbT)
{
  INVALIDARG_EXP_IFN(kbT >= 0.0, "Thermal energy kbT must be non-negative.");
  const Eigen::Index springEigenIdx = static_cast<Eigen::Index>(springIdx);
  const size_t partnerA =
    static_cast<size_t>(this->bondPartnersA[springEigenIdx]);
  const size_t partnerB =
    static_cast<size_t>(this->bondPartnersB[springEigenIdx]);

  size_t newPartnerA = partnerA;
  Eigen::Vector3d deltaRA = Eigen::Vector3d::Zero();
  Eigen::Vector3d deltaOffsetA = Eigen::Vector3d::Zero();
  if (!this->proposeRailStep(partnerA, newPartnerA, deltaRA, deltaOffsetA)) {
    return false;
  }

  size_t newPartnerB = partnerB;
  Eigen::Vector3d deltaRB = Eigen::Vector3d::Zero();
  Eigen::Vector3d deltaOffsetB = Eigen::Vector3d::Zero();
  if (!this->proposeRailStep(partnerB, newPartnerB, deltaRB, deltaOffsetB)) {
    return false;
  }

  if (this->idxFunctionalities[static_cast<Eigen::Index>(newPartnerA)] > 2 ||
      this->idxFunctionalities[static_cast<Eigen::Index>(newPartnerB)] > 2) {
    return false;
  }

  const Eigen::Vector3d bondDistanceNow =
    this->computeBondDistance(static_cast<int>(springIdx));
  const double bondEnergyNow = 0.5 * this->k * bondDistanceNow.squaredNorm();
  Eigen::Vector3d bondDistanceNew = bondDistanceNow - deltaRA + deltaRB;
  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(bondDistanceNew);
  }
  const double bondEnergyNew = 0.5 * this->k * bondDistanceNew.squaredNorm();
  const bool accept =
    this->acceptMetropolisMove(bondEnergyNew - bondEnergyNow, kbT);
  if (accept) {
    const Eigen::Vector3d newOffset =
      this->bondBoxOffsets.segment(3 * springEigenIdx, 3) - deltaOffsetA +
      deltaOffsetB;
    this->setSlipSpringPartners(
      springIdx, newPartnerA, newPartnerB, &newOffset);
#ifndef NDEBUG
    assert(
      (this->computeBondDistance(static_cast<int>(springIdx)) - bondDistanceNew)
        .cwiseAbs()
        .maxCoeff() < 1e-10);
#endif
  }
  return accept;
}

bool
DPDSimulator::attemptSlipSpringShift(const size_t springIdx,
                                     const size_t endToShift,
                                     const double kbT)
{
  INVALIDARG_EXP_IFN(kbT >= 0.0, "Thermal energy kbT must be non-negative.");
  const Eigen::Index springEigenIdx = static_cast<Eigen::Index>(springIdx);
  INVALIDARG_EXP_IFN(
    this->bondPartnersA[springEigenIdx] == static_cast<int>(endToShift) ||
      this->bondPartnersB[springEigenIdx] == static_cast<int>(endToShift),
    "This spring and its partners do not match, cannot attempt shift.");
  const bool shiftEndIsFirst =
    this->bondPartnersA[springEigenIdx] == static_cast<int>(endToShift);
  const size_t oldPartnerA =
    static_cast<size_t>(this->bondPartnersA[springEigenIdx]);
  const size_t oldPartnerB =
    static_cast<size_t>(this->bondPartnersB[springEigenIdx]);

  size_t replacementPartner = endToShift;
  Eigen::Vector3d deltaR = Eigen::Vector3d::Zero();
  Eigen::Vector3d deltaOffset = Eigen::Vector3d::Zero();
  if (!this->proposeRailStep(
        endToShift, replacementPartner, deltaR, deltaOffset)) {
    return false;
  }
  if (this->idxFunctionalities[static_cast<Eigen::Index>(replacementPartner)] >
      2) {
    return false;
  }

  const Eigen::Vector3d bondDistanceNow =
    this->computeBondDistance(static_cast<int>(springIdx));
  const double bondEnergyNow = 0.5 * this->k * bondDistanceNow.squaredNorm();
  Eigen::Vector3d bondDistanceNew =
    bondDistanceNow + (shiftEndIsFirst ? -1.0 : 1.0) * deltaR;
  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(bondDistanceNew);
  }
  const double bondEnergyNew = 0.5 * this->k * bondDistanceNew.squaredNorm();
  const bool accept =
    this->acceptMetropolisMove(bondEnergyNew - bondEnergyNow, kbT);
  if (accept) {
    const size_t newPartnerA =
      shiftEndIsFirst ? replacementPartner : oldPartnerA;
    const size_t newPartnerB =
      shiftEndIsFirst ? oldPartnerB : replacementPartner;
    const Eigen::Vector3d newOffset =
      this->bondBoxOffsets.segment(3 * springEigenIdx, 3) +
      (shiftEndIsFirst ? -1.0 : 1.0) * deltaOffset;
    this->setSlipSpringPartners(
      springIdx, newPartnerA, newPartnerB, &newOffset);
#ifndef NDEBUG
    assert(
      (this->computeBondDistance(static_cast<int>(springIdx)) - bondDistanceNew)
        .cwiseAbs()
        .maxCoeff() < 1e-10);
#endif
  }
  return accept;
}

void
DPDSimulator::setSlipSpringPartners(const size_t springIdx,
                                    const size_t newPartnerA,
                                    const size_t newPartnerB,
                                    const Eigen::Vector3d* const explicitOffset)
{
  const Eigen::Index springEigenIdx = static_cast<Eigen::Index>(springIdx);
  INVALIDINDEX_EXP_IFN(springIdx >= static_cast<size_t>(this->numBonds) &&
                         springEigenIdx < this->bondPartnersA.size(),
                       "Cannot replace on a non-existing slip-spring.");
  INVALIDINDEX_EXP_IFN(newPartnerA < static_cast<size_t>(this->numAtoms),
                       "Cannot replace with a non-existing bead.");
  INVALIDINDEX_EXP_IFN(newPartnerB < static_cast<size_t>(this->numAtoms),
                       "Cannot replace with a non-existing bead.");
  INVALIDARG_EXP_IFN(
    (this->atomTypes[newPartnerA] != this->crossLinkerType ||
     this->idxFunctionalities[static_cast<Eigen::Index>(newPartnerA)] <= 2) &&
      (this->atomTypes[newPartnerB] != this->crossLinkerType ||
       this->idxFunctionalities[static_cast<Eigen::Index>(newPartnerB)] <= 2),
    "Cannot form slip-springs with crosslinkers.");

  const size_t oldPartnerA =
    static_cast<size_t>(this->bondPartnersA[springEigenIdx]);
  const size_t oldPartnerB =
    static_cast<size_t>(this->bondPartnersB[springEigenIdx]);

  if (oldPartnerA == newPartnerA && oldPartnerB == newPartnerB) {
    if (explicitOffset != nullptr) {
      this->bondBoxOffsets.segment(3 * springEigenIdx, 3) = *explicitOffset;
    } else {
      this->resetBondOffset(static_cast<int>(springIdx));
    }
    return;
  }

  auto eraseOneFromPartner = [&](const size_t partnerIdx) {
    const Eigen::Index partnerEigen = static_cast<Eigen::Index>(partnerIdx);
    std::vector<size_t>& partnerBonds = this->bondsOfIndex[partnerIdx];
    for (size_t i = static_cast<size_t>(this->idxFunctionalities[partnerEigen]);
         i < partnerBonds.size();
         ++i) {
      if (partnerBonds[i] == springIdx) {
        partnerBonds.erase(partnerBonds.begin() +
                           static_cast<std::ptrdiff_t>(i));
        return;
      }
    }
    this->validateState();
    throw std::runtime_error("Invalid internal state: replacing slip-spring "
                             "partner, but did not find it internally.");
  };

  eraseOneFromPartner(oldPartnerA);
  eraseOneFromPartner(oldPartnerB);

  this->bondPartnersA[springEigenIdx] = static_cast<int>(newPartnerA);
  this->bondPartnersB[springEigenIdx] = static_cast<int>(newPartnerB);
  for (Eigen::Index dir = 0; dir < 3; ++dir) {
    this->bondPartnerCoordinatesA[3 * springEigenIdx + dir] =
      3 * static_cast<int>(newPartnerA) + static_cast<int>(dir);
    this->bondPartnerCoordinatesB[3 * springEigenIdx + dir] =
      3 * static_cast<int>(newPartnerB) + static_cast<int>(dir);
  }

  if (explicitOffset != nullptr) {
    this->bondBoxOffsets.segment(3 * springEigenIdx, 3) = *explicitOffset;
  } else {
    this->resetBondOffset(static_cast<int>(springIdx));
  }

  pylimer_tools::utils::addToSorted<size_t>(this->bondsOfIndex[newPartnerA],
                                            springIdx);
  pylimer_tools::utils::addToSorted<size_t>(this->bondsOfIndex[newPartnerB],
                                            springIdx);

  this->resetBondDuplicationPenalty(oldPartnerA);
  if (oldPartnerB != oldPartnerA) {
    this->resetBondDuplicationPenalty(oldPartnerB);
  }
  if (newPartnerA != oldPartnerA && newPartnerA != oldPartnerB) {
    this->resetBondDuplicationPenalty(newPartnerA);
  }
  if (newPartnerB != oldPartnerA && newPartnerB != oldPartnerB &&
      newPartnerB != newPartnerA) {
    this->resetBondDuplicationPenalty(newPartnerB);
  }
}

void
DPDSimulator::replaceSlipSpringPartner(const size_t springIdx,
                                       const size_t partnerBefore,
                                       const size_t partnerAfter)
{
  const Eigen::Index springEigenIdx = static_cast<Eigen::Index>(springIdx);
  INVALIDINDEX_EXP_IFN(springEigenIdx < this->bondPartnersA.size(),
                       "Cannot replace on a non-existing spring.");
  INVALIDARG_EXP_IFN(
    this->bondPartnersA[springEigenIdx] == static_cast<int>(partnerBefore) ||
      this->bondPartnersB[springEigenIdx] == static_cast<int>(partnerBefore),
    "This spring and its partners do not match.");
  if (partnerAfter == partnerBefore) {
    return;
  }
  if (this->bondPartnersA[springEigenIdx] == static_cast<int>(partnerBefore)) {
    this->setSlipSpringPartners(
      springIdx,
      partnerAfter,
      static_cast<size_t>(this->bondPartnersB[springEigenIdx]),
      nullptr);
  } else {
    this->setSlipSpringPartners(
      springIdx,
      static_cast<size_t>(this->bondPartnersA[springEigenIdx]),
      partnerAfter,
      nullptr);
  }
}

/**
 * @brief Convert the current structure to a Unvierse
 *
 * @return pylimer_tools::entities::Universe
 */
pylimer_tools::entities::Universe
DPDSimulator::getUniverse(const bool withSlipsprings) const
{
  pylimer_tools::entities::Universe result =
    pylimer_tools::entities::Universe(this->box);

  const size_t nAtomsU = static_cast<size_t>(this->numAtoms);
  std::vector<double> xs;
  xs.reserve(nAtomsU);
  std::vector<double> ys;
  ys.reserve(nAtomsU);
  std::vector<double> zs;
  zs.reserve(nAtomsU);
  std::vector<int> zeros;
  zeros.reserve(nAtomsU);

  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(this->numAtoms * 3);
       i += 3) {
    xs.push_back(this->coordinates[i]);
    ys.push_back(this->coordinates[i + 1]);
    zs.push_back(this->coordinates[i + 2]);
    zeros.push_back(0);
  }

  result.addAtoms(
    this->atomIds, this->atomTypes, xs, ys, zs, zeros, zeros, zeros);

  const size_t numBondsToAdd = static_cast<size_t>(
    withSlipsprings ? (this->numBonds + this->numSlipSprings) : this->numBonds);
  std::vector<long int> bondFrom;
  bondFrom.reserve(numBondsToAdd);
  std::vector<long int> bondTo;
  bondTo.reserve(numBondsToAdd);
  std::vector<int> newBondTypes;
  newBondTypes.reserve(numBondsToAdd);

  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numBondsToAdd); ++i) {
    bondTo.push_back(this->universe.getAtomIdByIdx(this->bondPartnersA[i]));
    bondFrom.push_back(this->universe.getAtomIdByIdx(this->bondPartnersB[i]));
    newBondTypes.push_back(this->bondTypes[i]);
  }

  result.addBonds(numBondsToAdd, bondFrom, bondTo, newBondTypes, false, false);

  return result;
}

void
DPDSimulator::validateNeighbourlist(const double cutoff)
{
  // this->neighbourlist.resetCoordinates(this->coordinates);
  RUNTIME_EXP_IFN(this->neighbourlist.getNumBinnedCoordinates() ==
                    this->coordinates.size() / 3,
                  "Not all coordinates are represented.");
  RUNTIME_EXP_IFN(
    this->neighbourlist.checkIfCoordinatesAreCurrent(this->coordinates),
    "Apparently, the neighbourlist's coordinates were not reset properly.");
  // pre-allocate the neighbor indices array
  Eigen::ArrayXi neighbors = Eigen::ArrayXi(static_cast<int>(
    this->numAtoms *
    (std::ceil((3.1 * cutoff) * (3.1 * cutoff) * (3.1 * cutoff)) /
     this->box.getVolume())));

  // actually loop the atoms
  for (size_t i = 0; i < static_cast<size_t>(this->numAtoms); ++i) {
    const Eigen::Index iEigen = static_cast<Eigen::Index>(i);
    const int numNeighbors = this->neighbourlist.getIndicesCloseToCoordinates(
      neighbors, this->coordinates.segment(3 * iEigen, 3), cutoff);
    Eigen::ArrayXi neighbors2 =
      this->neighbourlist.getIndicesCloseToCoordinates(
        this->coordinates.segment(3 * iEigen, 3), cutoff);
    RUNTIME_EXP_IFN(
      (neighbors.segment(0, numNeighbors) == neighbors2).all(),
      "Neighbors should be equal no matter the method, but apparently, "
      "they are not.");

    std::vector<size_t> relevantNeighbors;
    std::vector<size_t> relevantPairs;

    // neighbour-list
    for (Eigen::Index neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
      const size_t j = static_cast<size_t>(neighbors[neigh_idx]);
      if (j <= i) {
        continue;
      }
      const Eigen::Index jEigen = static_cast<Eigen::Index>(j);
      Eigen::Vector3d pairDistance = this->coordinates.segment(3 * iEigen, 3) -
                                     this->coordinates.segment(3 * jEigen, 3);
      this->box.handlePBC(pairDistance);
      const double rNorm = pairDistance.norm();
      if (rNorm >= cutoff || rNorm < 1e-12) {
        continue;
      }

      relevantNeighbors.push_back(j);
    }

    assert(relevantNeighbors.size() <= static_cast<size_t>(numNeighbors));

    // pairs
    for (size_t j = i + 1; j < static_cast<size_t>(this->numAtoms); ++j) {
      const Eigen::Index jEigen = static_cast<Eigen::Index>(j);
      Eigen::Vector3d pairDistance = this->coordinates.segment(3 * iEigen, 3) -
                                     this->coordinates.segment(3 * jEigen, 3);
      this->box.handlePBC(pairDistance);
      const double rNorm = pairDistance.norm();
      if (rNorm >= cutoff || rNorm < 1e-12) {
        continue;
      }

      relevantPairs.push_back(j);
      // bool found = false;
      // for (size_t k = 0; k < numNeighbors; ++k) {
      //   if (neighbors[k] == j) {
      //     found = true;
      //     break;
      //   }
      // }
      // RUNTIME_EXP_IFN(found,
      //                 "Did not find pair neighbour " + std::to_string(j)
      //                 +
      //                   " in list of neighbors of atom " +
      //                   std::to_string(i) + ".");
    }

    std::ranges::sort(relevantPairs);
    std::ranges::sort(relevantNeighbors);
    if (relevantPairs.size() > relevantNeighbors.size()) {
      std::cout << "Debugging neighbourlist. " << numNeighbors << std::endl;
      // debug why
      // find the difference
      std::vector<size_t> diff;
      std::ranges::set_difference(
        relevantPairs, relevantNeighbors, std::back_inserter(diff));
      assert(diff.size() >= (relevantPairs.size() - relevantNeighbors.size()));
      // figure out why not included
      for (const size_t diff_j : diff) {
        this->neighbourlist.validateWhyNotIncluded(
          this->coordinates.segment(3 * iEigen, 3),
          this->coordinates.segment(3 * static_cast<Eigen::Index>(diff_j), 3),
          cutoff);
      }
    }

    RUNTIME_EXP_IFN(
      relevantPairs.size() == relevantNeighbors.size(),
      "Pairs and neighbours resulted in different sized partners: " +
        std::to_string(relevantPairs.size()) + " vs. " +
        std::to_string(relevantNeighbors.size()) + " for atom at idx " +
        std::to_string(i) + ". Pair's neighbours are: " +
        pylimer_tools::utils::join(
          relevantPairs.begin(), relevantPairs.end(), std::string(", ")) +
        ". NeighbourList's neighbours are: " +
        pylimer_tools::utils::join(relevantNeighbors.begin(),
                                   relevantNeighbors.end(),
                                   std::string(", ")) +
        ".");

    RUNTIME_EXP_IFN(relevantNeighbors == relevantPairs,
                    "Pairs and neighbours are not equal.");
  }
}

/**
 * @brief Validate just what we are currently debugging
 *
 */
void
DPDSimulator::validateDebugState()
{
}

Eigen::VectorXd
DPDSimulator::getBondLengths()
{
  Eigen::VectorXd bondDistances =
    this->coordinates(this->bondPartnerCoordinatesB) -
    this->coordinates(this->bondPartnerCoordinatesA) + this->bondBoxOffsets;
  if (this->assumeBoxLargeEnough) {
    // this should not do anything anymore, here,
    // assuming the assumption holds.
    this->box.handlePBC(bondDistances);
  }

  return bondDistances
    .reshaped(3,
              static_cast<Eigen::Index>(this->numBonds + this->numSlipSprings))
    .colwise()
    .norm()
    .transpose();
}

Eigen::VectorXd
DPDSimulator::getCoordinates()
{
  assert(this->coordinates.size() == 3 * this->numAtoms);
  // std::cout << "DPDSim returning " << this->coordinates.size() << "
  // coordinates..."
  //           << std::endl;
  return this->coordinates;
}
Eigen::VectorXd
DPDSimulator::getVelocities() const
{
  return this->currentVelocities;
}

double
DPDSimulator::getTemperature()
{
  return this->computeTemperature(this->currentVelocities);
}

/**
 * @brief Make sure all the structures obey the expected form
 *
 */
void
DPDSimulator::validateState()
{
  // atoms
  RUNTIME_EXP_IFN(this->coordinates.size() == 3 * this->numAtoms,
                  "State violation: size of coordinates incorrect.");
  RUNTIME_EXP_IFN(this->idxFunctionalities.size() == this->numAtoms,
                  "State violation: size of functionalities incorrect.");
  RUNTIME_EXP_IFN(this->atomTypes.size() == static_cast<size_t>(this->numAtoms),
                  "State violation: size of atom types incorrect.");
  RUNTIME_EXP_IFN(this->atomIds.size() == static_cast<size_t>(this->numAtoms),
                  "State violation: size of atom ids incorrect.");
  RUNTIME_EXP_IFN(this->bondsOfIndex.size() ==
                    static_cast<size_t>(this->numAtoms),
                  "State violation: bonds of indices distributed incorrectly.");
  for (size_t i = 0; i < static_cast<size_t>(this->numAtoms); ++i) {
    const Eigen::Index iEigen = static_cast<Eigen::Index>(i);
    const int num_actual_bonds = std::accumulate(
      this->bondsOfIndex[i].begin(),
      this->bondsOfIndex[i].end(),
      0,
      [&](const int val, const size_t bondIdx) {
        return val +
               static_cast<int>(bondIdx < static_cast<size_t>(this->numBonds));
      });
    if (this->numSlipSprings == 0) {
      RUNTIME_EXP_IFN(num_actual_bonds ==
                        static_cast<int>(this->bondsOfIndex[i].size()),
                      "State violation: bonds of index " + std::to_string(i) +
                        " were not counted appropriately.");
    }
    RUNTIME_EXP_IFN(this->idxFunctionalities[iEigen] == num_actual_bonds,
                    "State violation: inconsistent idx functionalities");
    RUNTIME_EXP_IFN(
      this->idxFunctionalities[iEigen] >= 1 ||
        this->atomTypes[i] == this->crossLinkerType,
      "Only crosslinkers are allowed to be single beads. Found bead " +
        std::to_string(i) + " to have type " +
        std::to_string(this->atomTypes[i]) + " with functionality " +
        std::to_string(this->idxFunctionalities[iEigen]) + ".");
    if (this->idxFunctionalities[iEigen] == 1) {
      RUNTIME_EXP_IFN(
        this->isRelocationTarget[iEigen],
        "Expected functionalities of 1 to be relocation targets.");
    }
    for (size_t j = 0; j < this->bondsOfIndex[i].size(); ++j) {
      if (j < static_cast<size_t>(this->idxFunctionalities[iEigen])) {
        RUNTIME_EXP_IFN(this->bondsOfIndex[i][j] <
                          static_cast<size_t>(this->numBonds),
                        "Expect bonds to come before slip springs.");
      } else {
        RUNTIME_EXP_IFN(this->bondsOfIndex[i][j] >=
                          static_cast<size_t>(this->numBonds),
                        "Expect bonds to come before slip springs.");
      }
      if (j > 0) {
        RUNTIME_EXP_IFN(this->bondsOfIndex[i][j] > this->bondsOfIndex[i][j - 1],
                        "Require bonds of atom to be sorted.");
      }
    }
  }
  RUNTIME_EXP_IFN(this->isRelocationTarget.size() == this->numAtoms,
                  "State violation: size of relocation targets incorrect.");
  RUNTIME_EXP_IFN(this->isRelocationTarget.cast<int>().sum() ==
                    static_cast<Eigen::Index>(this->chainEndIndices.size()),
                  "Expect relocation targets to be equal to chain ends.");
  for (const size_t i : this->chainEndIndices) {
    const Eigen::Index iEigen = static_cast<Eigen::Index>(i);
    RUNTIME_EXP_IFN(this->isRelocationTarget[iEigen],
                    "Expect chain ends to be relocation targets.");
    RUNTIME_EXP_IFN(this->idxFunctionalities[iEigen] <= 1 ||
                      this->allowRelocationInNetwork,
                    "Expect chain ends to have a functionality of 1, have "
                    "relocation in network enabled");
  }
  // bonds
  RUNTIME_EXP_IFN(this->bondPartnersA.size() == this->bondPartnersB.size(),
                  "State violation: nr of bonds inconsistent.");
  RUNTIME_EXP_IFN(this->bondPartnerCoordinatesA.size() ==
                    3 * this->bondPartnersA.size(),
                  "State violation: nr of bonds inconsistent.");
  RUNTIME_EXP_IFN(this->bondPartnerCoordinatesB.size() ==
                    3 * this->bondPartnersB.size(),
                  "State violation: nr of bonds inconsistent.");
  RUNTIME_EXP_IFN(this->bondTypes.size() == this->bondPartnersA.size(),
                  "State violation: nr of bonds inconsistent.");
  RUNTIME_EXP_IFN(this->bondPartnersB.size() ==
                    this->numBonds + this->numSlipSprings,
                  "State violation: nr of bonds inconsistent.");
  RUNTIME_EXP_IFN(this->bondBoxOffsets.size() ==
                    (this->numBonds + this->numSlipSprings) * 3,
                  "State violation: nr of box bond offsets incorrect.");
  RUNTIME_EXP_IFN(
    this->bondDuplicationPenalty.size() ==
      (this->numBonds + this->numSlipSprings) * 3,
    "State violation: nr of bond duplication penalities incorrect: got " +
      std::to_string(this->bondDuplicationPenalty.size()) + " for " +
      std::to_string(this->numBonds + this->numSlipSprings) + "bonds.");

  // bond duplication penalty
  std::unordered_set<size_t> partners;
  for (size_t i = 0; i < static_cast<size_t>(this->numAtoms); ++i) {
    partners.clear();
    // here as well, we rely on the bonds of index being sorted
    for (const size_t bondIdx : this->bondsOfIndex[i]) {
      const Eigen::Index bondEigenIdx = static_cast<Eigen::Index>(bondIdx);
      const size_t atomPartnerIdx = static_cast<size_t>(
        this->bondPartnersA[bondEigenIdx] == static_cast<int>(i)
          ? this->bondPartnersB[bondEigenIdx]
          : this->bondPartnersA[bondEigenIdx]);
      bool expectZero = false;
      if (partners.contains(atomPartnerIdx)) {
        // "real" bonds always contribute -> check that this is a
        // slip-spring
        if (bondIdx >= static_cast<size_t>(this->numBonds)) {
          expectZero = true;
        }
      } else {
        partners.insert(atomPartnerIdx);
      }
      RUNTIME_EXP_IFN(
        this->bondDuplicationPenalty.segment(3 * bondEigenIdx, 3)
          .isApprox(Eigen::Array3d::Constant(expectZero ? 0. : 1.)),
        "Incorrect bondDuplicationPenalty: found bond " +
          std::to_string(bondIdx) + " of " + std::to_string(this->numBonds) +
          " (+ " + std::to_string(this->numSlipSprings) + ") to should have " +
          std::to_string(expectZero ? 0 : 1) + " but it did not (" +
          std::to_string(
            this->bondDuplicationPenalty.segment(3 * bondEigenIdx, 3).sum()) +
          ").");
    }
  }

  // internal structure
  RUNTIME_EXP_IFN((this->bondPartnersB.array() < this->numAtoms).all(),
                  "State violation: too large indices found (e.g. " +
                    std::to_string(this->bondPartnersB.maxCoeff()) + " for " +
                    std::to_string(this->numAtoms) + " atoms)");
  RUNTIME_EXP_IFN((this->bondPartnersA.array() < this->numAtoms).all(),
                  "State violation: too large indices found (e.g. " +
                    std::to_string(this->bondPartnersA.maxCoeff()) + " for " +
                    std::to_string(this->numAtoms) + " atoms)");
  RUNTIME_EXP_IFN(
    !(this->bondPartnersA.array() == this->bondPartnersB.array()).any(),
    "Bonds with oneself are not allowed.");
  // loop springs
  for (Eigen::Index i = 0; i < this->numBonds + this->numSlipSprings; ++i) {
    for (Eigen::Index dir = 0; dir < 3; ++dir) {
      RUNTIME_EXP_IFN(
        this->bondPartnerCoordinatesA[i * 3 + dir] ==
          this->bondPartnersA[i] * 3 + static_cast<int>(dir),
        "Bond partners not accurate: got partner coordinate idx " +
          std::to_string(this->bondPartnerCoordinatesA[i * 3 + dir]) +
          " but expected " +
          std::to_string(this->bondPartnersA[i] * 3 + static_cast<int>(dir)) +
          " at index " + std::to_string(i) + " and dir " + std::to_string(dir) +
          ".");
      RUNTIME_EXP_IFN(
        this->bondPartnerCoordinatesB[i * 3 + dir] ==
          this->bondPartnersB[i] * 3 + static_cast<int>(dir),
        "Bond partners not accurate: got partner coordinate idx " +
          std::to_string(this->bondPartnerCoordinatesB[i * 3 + dir]) +
          " but expected " +
          std::to_string(this->bondPartnersB[i] * 3 + static_cast<int>(dir)) +
          " at index " + std::to_string(i) + " and dir " + std::to_string(dir) +
          ".");
    }
    RUNTIME_EXP_IFN(
      pylimer_tools::utils::contains(
        this->bondsOfIndex[static_cast<size_t>(this->bondPartnersA[i])],
        static_cast<size_t>(i)),
      "Reverse-link is incorrect. Bond " + std::to_string(i) + " (of " +
        std::to_string(this->numBonds) +
        ") was not found in bonds of partner A, " +
        std::to_string(this->bondPartnersA[i]) + ".");
    RUNTIME_EXP_IFN(
      pylimer_tools::utils::contains(
        this->bondsOfIndex[static_cast<size_t>(this->bondPartnersB[i])],
        static_cast<size_t>(i)),
      "Reverse-link is incorrect. Bond " + std::to_string(i) + " (of " +
        std::to_string(this->numBonds) +
        ") was not found in bonds of partner B, " +
        std::to_string(this->bondPartnersA[i]) + ".");
    // the following check is incorrect, as crosslinkers may have
    // slip-springs for f < 2, and during crosslinking might have one while
    // being upgraded to f > 2
    // if (i >= this->numBonds) {
    //   RUNTIME_EXP_IFN(this->atomType[this->bondPartnersA[i]] !=
    //                     this->crossLinkerType,
    //                   "Expect slip-links to not involve crosslinkers");
    //   RUNTIME_EXP_IFN(this->atomType[this->bondPartnersB[i]] !=
    //                     this->crossLinkerType,
    //                   "Expect slip-links to not involve crosslinkers");
    // }
  }
}

}
