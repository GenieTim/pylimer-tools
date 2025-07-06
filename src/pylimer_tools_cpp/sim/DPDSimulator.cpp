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
                           const std::string& seed)
  : box(u.getBox())
  , universe(u)
  , neighbourlist(
      pylimer_tools::entities::EigenNeighbourList(Eigen::VectorXd(0),
                                                  this->box,
                                                  1.0))
{
  INVALIDARG_EXP_IFN(!is2D, "2D simulations are not supported yet.");
  this->is2D = is2D;
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
  this->bondPartnersA = Eigen::ArrayXi::Zero(edges["edge_from"].size());
  this->bondPartnersB = Eigen::ArrayXi::Zero(edges["edge_to"].size());
  this->bondTypes = Eigen::ArrayXi::Zero(edges["edge_type"].size());

  this->numBonds = 0;
  this->numSlipSprings = 0;
  // we have to "sort" the edges such that the slip-springs are at the end
  // so, first add "normal" bonds
  for (size_t i = 0; i < edges["edge_from"].size(); i++) {
    if (edges["edge_type"][i] != slipspringBondType) {
      this->bondPartnersA(this->numBonds) = edges["edge_from"][i];
      this->bondPartnersB(this->numBonds) = edges["edge_to"][i];
      this->bondTypes(this->numBonds) = edges["edge_type"][i];
      this->numBonds += 1;
    }
  }
  // then, iterate again to add slip-springs
  for (size_t i = 0; i < edges["edge_from"].size(); i++) {
    if (edges["edge_type"][i] == slipspringBondType) {
      this->bondPartnersA(this->numBonds + this->numSlipSprings) =
        edges["edge_from"][i];
      this->bondPartnersB(this->numBonds + this->numSlipSprings) =
        edges["edge_to"][i];
      this->bondTypes(this->numBonds + this->numSlipSprings) =
        edges["edge_type"][i];
      this->numSlipSprings += 1;
    }
  }

  this->neighbourlist =
    pylimer_tools::entities::EigenNeighbourList(coordinates, this->box, 1.0);
  this->numAtoms = this->coordinates.size() / 3;
  this->idxFunctionalities = Eigen::ArrayXi::Zero(this->numAtoms);
  this->atomTypes = u.getPropertyValues<int>("type");
  this->atomIds = u.getPropertyValues<long int>("id");
  this->maxBondLen = 0.45 * this->box.getL().maxCoeff();
  this->crossLinkerType = crossLinkerType;
  this->slipspringBondType = slipspringBondType;

  this->bondsOfIndex.reserve(this->numAtoms);
  for (size_t i = 0; i < this->numAtoms; ++i) {
    std::vector<size_t> bonds;
    bonds.reserve(4);
    this->bondsOfIndex.push_back(bonds);
  }
  this->bondPartnerCoordinatesA =
    Eigen::ArrayXi(3 * this->bondPartnersA.size());
  this->bondPartnerCoordinatesB =
    Eigen::ArrayXi(3 * this->bondPartnersB.size());
  for (size_t i = 0; i < this->numBonds + this->numSlipSprings; ++i) {
    this->bondsOfIndex[this->bondPartnersA[i]].push_back(i);
    this->bondsOfIndex[this->bondPartnersB[i]].push_back(i);
    for (int dir = 0; dir < 3; ++dir) {
      this->bondPartnerCoordinatesA[i * 3 + dir] =
        this->bondPartnersA[i] * 3 + dir;
      this->bondPartnerCoordinatesB[i * 3 + dir] =
        this->bondPartnersB[i] * 3 + dir;
    }
  }
  this->isRelocationTarget = Eigen::ArrayXb::Constant(this->numAtoms, false);
  for (size_t i = 0; i < this->numAtoms; ++i) {
    // if we already have slip-springs, the functionality must account for
    // that: The idxFunctionalities should not include the slip-springs
    this->idxFunctionalities[i] =
      this->numSlipSprings == 0
        ? this->bondsOfIndex[i].size()
        : std::accumulate(this->bondsOfIndex[i].begin(),
                          this->bondsOfIndex[i].end(),
                          0,
                          [&](const int val, const size_t bondIdx) {
                            return val +
                                   static_cast<int>(bondIdx < this->numBonds);
                          });
    if (this->idxFunctionalities[i] < 2) {
      this->chainEndIndices.push_back(i);
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
      this->currentVelocities(i * 3) = vx[i];
      this->currentVelocities(i * 3 + 1) = vy[i];
      this->currentVelocities(i * 3 + 2) = vz[i];
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
 * @brief
 *
 * @param nSteps
 * @param withMC
 * @param shouldInterrupt
 * @param cleanupInterrupt
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
      this->numShifts = 0;
      timer.section(DPDPerformanceSections::RELOCATION);
      this->numRelocations = this->relocateSlipSprings(1. * temperature);
      timer.section(DPDPerformanceSections::SHIFT);
      for (int i = 0; i < this->nStepsMC; ++i) {
        this->numShifts += this->shiftSlipSprings(1. * temperature);
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
    for (size_t i = 0; i < this->numAtoms; ++i) {
      this->currentStressTensor +=
        m_over_boxv * this->currentVelocities.segment(3 * i, 3) *
        this->currentVelocities.segment(3 * i, 3).transpose();
    }

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
  for (size_t i = 0; i < this->numAtoms; ++i) {
    this->currentStressTensor +=
      (m_over_boxv * this->currentVelocities.segment(3 * i, 3) *
       this->currentVelocities.segment(3 * i, 3).transpose());
  }
}

/**
 * @brief Set a new seed for the random generator
 *
 * @param seed
 */
void
DPDSimulator::reseedRandomness(const std::string& seed)
{
  // initialize the random number generator
  if (seed == "") {
    std::random_device rd;
    this->e2 = std::mt19937(rd());
  } else {
    std::seed_seq seed2(seed.begin(), seed.end());
    this->e2 = std::mt19937(seed2);
  }
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

  for (size_t i = 0; i < this->numAtoms - 1; ++i) {
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
    size_t atomPartnerIdx = this->bondPartnersA[bondIdx] == atomIdx
                              ? this->bondPartnersB[bondIdx]
                              : this->bondPartnersA[bondIdx];
    if (bondIdx >= this->numBonds) {
      // all others should stay 1 anyway
      this->bondDuplicationPenalty.segment(3 * bondIdx, 3).setConstant(1.);
    }
    if (partners.contains(atomPartnerIdx)) {
      // "real" bonds always contribute -> check that this is a
      // slip-spring
      if (bondIdx >= this->numBonds) {
        this->bondDuplicationPenalty.segment(3 * bondIdx, 3).setZero();
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
  for (size_t atom_idx = 0; atom_idx < this->numAtoms; ++atom_idx) {
    possibleCandidates.clear();
    if (this->atomTypes[atom_idx] != this->atomTypeBondFormationFrom) {
      continue;
    }
    if (this->idxFunctionalities[atom_idx] >=
        this->maxBondsPerType[this->atomTypes[atom_idx]]) {
      continue;
    }

    // find neighbours
    const int numNeighbors = this->neighbourlist.getIndicesCloseToCoordinates(
      neighbors, this->coordinates.segment(3 * atom_idx, 3), cutoff);
    for (size_t neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
      // loop neighbours to find applicable partner
      const size_t j = neighbors[neigh_idx];
      Eigen::Vector3d vec = this->coordinates.segment(3 * atom_idx, 3) -
                            this->coordinates.segment(3 * j, 3);
      this->box.handlePBC(vec);
      const double r2 = (vec).norm();
      if (r2 <= cutoff && atom_idx != j &&
          (this->idxFunctionalities[j] <
           this->maxBondsPerType[this->atomTypes[j]]) &&
          (this->atomTypeBondFormationTo == this->atomTypes[j])) {
        possibleCandidates.push_back(j);
      }
    }

    if (possibleCandidates.size() == 0) {
      continue;
    }

    // yay, we can form a bond
    // NOTE: currently, we build only 1
    std::ranges::shuffle(possibleCandidates, this->e2);
    this->addBond(
      atom_idx,
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
  for (size_t i = 0; i < this->numAtoms; ++i) {
    for (size_t j = 0; j < this->bondsOfIndex[i].size(); ++j) {
      if (this->bondsOfIndex[i][j] >= newBondIdx) {
        this->bondsOfIndex[i][j] += 1;
      }
    }
  }

  // actually register the new bond
  this->bondTypes[newBondIdx] = bondType;
  this->bondPartnersA[newBondIdx] = fromIdx;
  this->bondPartnersB[newBondIdx] = toIdx;
  for (size_t i = 0; i < 3; ++i) {
    this->bondPartnerCoordinatesA[3 * newBondIdx + i] = 3 * fromIdx + i;
    this->bondPartnerCoordinatesB[3 * newBondIdx + i] = 3 * toIdx + i;
  }
  pylimer_tools::utils::addToSorted<size_t>(this->bondsOfIndex[fromIdx],
                                            newBondIdx);
  pylimer_tools::utils::addToSorted<size_t>(this->bondsOfIndex[toIdx],
                                            newBondIdx);

  this->resetBondOffset(newBondIdx);
  this->resetBondDuplicationPenalty(fromIdx);
  this->resetBondDuplicationPenalty(toIdx);

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
                                                      fromIdx);
    }
    if (!this->isRelocationTarget[toIdx]) {
      pylimer_tools::utils::removeIfContained<size_t>(this->chainEndIndices,
                                                      toIdx);
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
  bool isDifferent = newBox == this->box;
  const pylimer_tools::entities::Box previousBox = this->box;
  // rescale box offsets
  this->box = newBox;
  assert(newBox == this->box);
  const Eigen::VectorXd newOffsets =
    (this->bondBoxOffsets.array() *
     (this->box.getL() / previousBox.getL())
       .replicate((this->numBonds + this->numSlipSprings), 1))
      .matrix();
  this->bondBoxOffsets = newOffsets;

  // rescale coordinates
  Eigen::VectorXd previousCoords = this->coordinates;

  Eigen::Array3d scalingFactor = newBox.getL() / previousBox.getL();

  // rescale per-atom quantities
  for (size_t i = 0; i < this->numAtoms; ++i) {
    for (size_t dir = 0; dir < 3; ++dir) {
      this->coordinates[3 * i + dir] *= scalingFactor[dir];
      this->currentVelocities[3 * i + dir] *= scalingFactor[dir];
      this->currentVelocitiesPlus[3 * i + dir] *= scalingFactor[dir];
      this->currentForces[3 * i + dir] *= scalingFactor[dir];
    }
  }

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
  // configuration
  constexpr double dim = 3;
  constexpr double kb = 1.0;
  constexpr double m = 1.0;

  const double KE = 0.5 * m * velocities.squaredNorm();

  return KE / ((dim / 2.) * (velocities.size() / dim) * kb);
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
  const double dt,
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
    for (size_t i = 0; i < this->bondPartnersA.size(); ++i) {
      if (bondDistances.segment(3 * i, 3).squaredNorm() >
          0.5 * this->box.getL().minCoeff()) {
        std::cerr
          << "WARNING: Bond " << i << " has length "
          << bondDistances.segment(3 * i, 3).squaredNorm()
          << ", which violates the assumption that the box is large enough"
          << std::endl;
      }
    }
#endif
  } else {
#ifndef NDEBUG
    for (size_t i = 0; i < this->bondPartnersA.size(); ++i) {
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
#pragma omp parallel for reduction(+ : stressTensor, pressure)                 \
  schedule(static, 16)
  for (size_t i = 0; i < this->bondPartnersA.size(); ++i) {
    // attractive force -> reduces pressure in the system
    pressure -= this->k * bondDistances.segment(3 * i, 3).squaredNorm();
    stressTensor -= this->k * bondDistances.segment(3 * i, 3) *
                    bondDistances.segment(3 * i, 3).transpose();
    assert(std::isfinite(bondDistances.segment(3 * i, 3).squaredNorm()));
  }

  timer.section(DPDPerformanceSections::PAIR_FORCE);

  // actually loop the atoms
#pragma omp parallel
  {
    const double sigmaOverSqrtDt = this->sigma / std::sqrt(dt);
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
    for (size_t i = 0; i < this->numAtoms - 1; ++i) {
      const int numNeighbors = this->neighbourlist.getHigherIndicesWithinCutoff(
        neighbors, coords, i, cutoff);

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
      for (size_t neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
        const size_t j = neighbors[neigh_idx];
        pairdistance = coords.segment(3 * i, 3) - coords.segment(3 * j, 3);
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
          velocities.segment(3 * i, 3) - velocities.segment(3 * j, 3);
        const double rij_dot_vij = pairdistanceNormed.dot(velocitydiff);
        const double gamma_weighted_rij_dot_vij =
          this->gamma * one_minus_rnorm2 * rij_dot_vij;

        // conservative repulsion force - dissipative/drag force
        // TODO: check if we get fma here
        double pairForceConst =
          this->A * one_minus_rnorm - gamma_weighted_rij_dot_vij;

        // random force
        const double constant_rnd_prefix = one_minus_rnorm * sigmaOverSqrtDt;
        const double random_val = this->uniform_rand_mean0std1(this->e2);

        pairForceConst += constant_rnd_prefix * random_val;

        // back to 3D/Eigen Vector
        pairForce = pairForceConst * pairdistanceNormed;

        // actually assign the new forces
        forces.segment(3 * i, 3) += pairForce;
        forces.segment(3 * j, 3) -= pairForce;

        // pressure update: repulsive force -> increases pressure
        pressure += pairForce.dot(pairdistance);
        assert(std::isfinite(pressure));
        stressTensor += pairForce * pairdistance.transpose();
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
    Eigen::ArrayXi(3 * atomIdsToMeasure.size());
#pragma omp parallel for
  for (size_t i = 0; i < atomIdsToMeasure.size(); ++i) {
    const size_t atomId = atomIdsToMeasure[i];
    const size_t index = this->universe.getIdxByAtomId(atomId);
    coordinateIndices[3 * i] = 3 * index;
    coordinateIndices[3 * i + 1] = 3 * index + 1;
    coordinateIndices[3 * i + 2] = 3 * index + 2;
  }

  // Remember to measure these relative to the current time-step
  msdMeasuredIndices.push_back(coordinateIndices);
  msdOrigins.push_back(this->coordinates(coordinateIndices));
  msdOriginTimesteps.push_back(this->currentStep);
}

/**
 * @brief Randomly add new slip-springs
 *
 * @param num
 * @param bondType
 * @return int
 */
int
DPDSimulator::createSlipSprings(const int num, int bondType)
{
  if (bondType == 0) {
    bondType = this->slipspringBondType;
  }
  int createdLastIteration = 100;
  int totalCreated = 0;
  std::vector<size_t> candidates;
  candidates.reserve(5);

  std::vector<size_t> slipSpringFrom;
  slipSpringFrom.reserve(num);
  std::vector<size_t> slipSpringTo;
  slipSpringTo.reserve(num);

  // randomly permute the atoms to start the search with
  // we do this over just randomly selecting an index
  // in order to reduce the probability of finding the same start twice
  std::vector<size_t> sourceIds;
  sourceIds.reserve(this->numAtoms);
  for (size_t i = 0; i < this->numAtoms; ++i) {
    if (this->atomTypes[i] != this->crossLinkerType) {
      sourceIds.push_back(i);
    }
  }

  // search for neighbours that are elibile
  Eigen::ArrayXi neighbours = Eigen::ArrayXi(16);
  while (createdLastIteration > 0 && totalCreated < num) {
    createdLastIteration = 0;
    std::ranges::shuffle(sourceIds, this->e2);
    for (size_t i : sourceIds) {
      int numCandidates = 0;
      // for each atom, search for possible partners
      const int numNeighs = this->neighbourlist.getIndicesCloseToCoordinates(
        neighbours, this->coordinates.segment(3 * i, 3), this->highCutoff);
      for (size_t j = 0; j < numNeighs; ++j) {
        if (this->atomTypes[neighbours[j]] == this->crossLinkerType) {
          continue;
        }
        Eigen::Vector3d distance =
          this->coordinates.segment(3 * i, 3) -
          this->coordinates.segment(3 * neighbours[j], 3);
        this->box.handlePBC(distance);
        if (distance.norm() > this->lowCutoff &&
            distance.norm() <= this->highCutoff) {
          if (numCandidates < candidates.size()) {
            candidates[numCandidates] = neighbours[j];
          } else {
            candidates.push_back(neighbours[j]);
          }
          numCandidates += 1;
        }
      }
      if (numCandidates == 0) {
        continue;
      }
      std::uniform_int_distribution<int> dist(0, numCandidates - 1);
      const int candidateIndex = dist(this->e2);
      // found a candidate to create a slip-spring to
      slipSpringFrom.push_back(i);
      slipSpringTo.push_back(candidates[candidateIndex]);
      totalCreated += 1;
      createdLastIteration += 1;
      if (totalCreated >= num) {
        break;
      }
    }
  }

  this->addSlipSprings(slipSpringFrom, slipSpringTo, this->slipspringBondType);
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
    INVALIDINDEX_EXP_IFN(partnerA[i] < this->numAtoms, "Invalid partner id");
    INVALIDINDEX_EXP_IFN(partnerB[i] < this->numAtoms, "Invalid partner id");
  }
  const size_t sizeBefore = this->numBonds + this->numSlipSprings;
  this->bondPartnersA.conservativeResize(sizeBefore + partnerA.size());
  this->bondPartnersB.conservativeResize(sizeBefore + partnerB.size());
  this->bondPartnerCoordinatesA.conservativeResize(
    3 * (sizeBefore + partnerA.size()));
  this->bondPartnerCoordinatesB.conservativeResize(
    3 * (sizeBefore + partnerB.size()));
  this->bondTypes.conservativeResize(sizeBefore + partnerB.size());
  this->bondBoxOffsets.conservativeResize(3 * (sizeBefore + partnerA.size()));
  this->bondDuplicationPenalty.conservativeResize(
    3 * (sizeBefore + partnerA.size()));

  this->bondPartnersA.segment(sizeBefore, partnerA.size()) =
    Eigen::Map<Eigen::ArrayXst, Eigen::Unaligned>(partnerA.data(),
                                                  partnerA.size())
      .cast<int>();
  this->bondPartnersB.segment(sizeBefore, partnerB.size()) =
    Eigen::Map<Eigen::ArrayXst, Eigen::Unaligned>(partnerB.data(),
                                                  partnerB.size())
      .cast<int>();
  this->bondTypes.segment(sizeBefore, partnerB.size()) = bondType;

  for (size_t i = sizeBefore; i < sizeBefore + partnerA.size(); ++i) {
    this->bondsOfIndex[this->bondPartnersA[i]].push_back(i);
    this->bondsOfIndex[this->bondPartnersB[i]].push_back(i);
    for (int dir = 0; dir < 3; ++dir) {
      this->bondPartnerCoordinatesA[i * 3 + dir] =
        this->bondPartnersA[i] * 3 + dir;
      this->bondPartnerCoordinatesB[i * 3 + dir] =
        this->bondPartnersB[i] * 3 + dir;
    }
    this->resetBondOffset(i);
    this->resetBondDuplicationPenalty(this->bondPartnersA[i]);
    this->resetBondDuplicationPenalty(this->bondPartnersB[i]);
  }

  this->numSlipSprings += partnerA.size();
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
  int nAccept = 0;
  std::vector<size_t> candidates;
  Eigen::ArrayXi neighbours = Eigen::ArrayXi(12);
  std::uniform_int_distribution<int> chainendDist(
    0, this->chainEndIndices.size() - 1);
  for (size_t springIdx = this->numBonds;
       springIdx < (this->numBonds + this->numSlipSprings);
       ++springIdx) {
    if ((!(this->isRelocationTarget[this->bondPartnersA[springIdx]] ||
           this->isRelocationTarget[this->bondPartnersB[springIdx]])) &&
        // always relocate from crosslinkers away, e.g. in case the
        // crosslink has gained an additional bond and is no longer a
        // relocation target.
        // consequently, we are inconsistent for crosslinkers with f = 2
        !(this->atomTypes[this->bondPartnersA[springIdx]] ==
            this->crossLinkerType ||
          this->atomTypes[this->bondPartnersB[springIdx]] ==
            this->crossLinkerType)) {
      continue;
    }
    // design & attempt move
    const size_t partnerA = this->bondPartnersA[springIdx];
    const size_t partnerB = this->bondPartnersB[springIdx];

    const int candidateIndex = this->chainEndIndices[chainendDist(this->e2)];

    // search for neighbours
    int numCandidates = 0;
    // for each atom, search for possible partners
    const int numNeighs = this->neighbourlist.getIndicesCloseToCoordinates(
      neighbours,
      this->coordinates.segment(3 * candidateIndex, 3),
      this->highCutoff);
    for (size_t j = 0; j < numNeighs; ++j) {
      // slip-springs to crosslinkers are not allowed
      if (this->atomTypes[neighbours[j]] == this->crossLinkerType) {
        continue;
      }
      Eigen::Vector3d distance =
        this->coordinates.segment(3 * candidateIndex, 3) -
        this->coordinates.segment(3 * neighbours[j], 3);
      this->box.handlePBC(distance);
      if (distance.norm() > this->lowCutoff &&
          distance.norm() <= this->highCutoff) {
        if (numCandidates >= candidates.size()) {
          candidates.push_back(neighbours[j]);
        } else {
          candidates[numCandidates] = neighbours[j];
        }
        numCandidates += 1;
      }
    }
    if (numCandidates == 0) {
      continue;
    }
    std::uniform_int_distribution<int> candidateDist(0, numCandidates - 1);
    const int candidatePartnerIndex = candidates[candidateDist(this->e2)];

    // compute the Metropolis criterion
    Eigen::Vector3d bondDistanceNow =
      (this->coordinates.segment(partnerA * 3, 3) -
       this->coordinates.segment(partnerB * 3, 3)) +
      this->bondBoxOffsets.segment(3 * springIdx, 3);
    if (this->assumeBoxLargeEnough) {
      this->box.handlePBC(bondDistanceNow);
    }
    const double bondEnergyNow = k * bondDistanceNow.squaredNorm();
    // for the new slip-spring, we cannot rely on bond box offsets.
    // it really is a new slip-spring
    Eigen::Vector3d bondDistanceNew =
      this->coordinates.segment(candidateIndex * 3, 3) -
      this->coordinates.segment(candidatePartnerIndex * 3, 3);
    this->box.handlePBC(bondDistanceNew);
    const double bondEnergyNew = k * bondDistanceNew.squaredNorm();
    const double deltaEnergy = bondEnergyNew - bondEnergyNow;
    bool accept = false;
    if (deltaEnergy < 0.0) {
      accept = true;
    } else {
      const double factor = std::exp(-deltaEnergy / kbT);
      if (this->uniform_rand_between_0_1(this->e2) < factor) {
        accept = true;
      }
    }
    if (accept) {
      this->replaceSlipSpringPartner(springIdx, partnerA, candidateIndex);
      this->replaceSlipSpringPartner(
        springIdx, partnerB, candidatePartnerIndex);
      nAccept++;
      RUNTIME_EXP_IFN(
        this->computeBondLength(springIdx) <= this->highCutoff,
        "By relocation, newly created bonds should not be too long.");
    }
  }
  return nAccept;
}

/**
 * @brief The first MC procedure:
 *
 * @param kbT
 * @return int
 */
int
DPDSimulator::shiftSlipSprings(const double kbT)
{
  int n_accept = 0;
  std::vector<size_t> candidates;
  Eigen::ArrayXi neighbours = Eigen::ArrayXi(12);

  std::uniform_int_distribution<int> uniformDistNatoms(0, this->numAtoms - 1);
  for (size_t springIdx = this->numBonds;
       springIdx < (this->numBonds + this->numSlipSprings);
       ++springIdx) {
    const size_t partnerA = this->bondPartnersA[springIdx];
    const size_t partnerB = this->bondPartnersB[springIdx];
    // attempt to shift the spring around partnerA
    if (this->shiftOneAtATime) {
      n_accept += this->attemptSlipSpringShift(springIdx, partnerA);
      n_accept += this->attemptSlipSpringShift(springIdx, partnerB);
    } else {
      n_accept += this->attemptSlipSpringShift(springIdx);
    }
    if (this->bondPartnersA[springIdx] == this->bondPartnersB[springIdx]) {
      // complete relocation of this bond
      int firstPartner = uniformDistNatoms(this->e2);
      size_t n_attempts = 0;
      while (this->atomTypes[firstPartner] == this->crossLinkerType) {
        firstPartner = uniformDistNatoms(this->e2);
        n_attempts++;
        RUNTIME_EXP_IFN(n_attempts < 1000,
                        "Too many times, a crosslink was chosen randomly.");
      }
      // search for neighbours
      int numCandidates = 0;
      // for this first chosen atom, search for possible partners
      const int numNeighs = this->neighbourlist.getIndicesCloseToCoordinates(
        neighbours,
        this->coordinates.segment(3 * firstPartner, 3),
        this->highCutoff);
      RUNTIME_EXP_IFN(numNeighs <= neighbours.size(),
                      "Neighbourlist does not act as it should.");
      for (size_t j = 0; j < numNeighs; ++j) {
        RUNTIME_EXP_IFN(neighbours[j] < this->numAtoms,
                        "Neighbourlist seems inconsequent.");
        if (this->idxFunctionalities[neighbours[j]] > 2) {
          // we won't relocate to crosslinkers
          continue;
        }
        Eigen::Vector3d distance =
          this->coordinates.segment(3 * firstPartner, 3) -
          this->coordinates.segment(3 * neighbours[j], 3);
        // here as well, we cannot rely on the bond box offsets,
        // as this is a completely new slip-spring (bond)
        this->box.handlePBC(distance);
        if (distance.norm() > this->lowCutoff &&
            distance.norm() <= this->highCutoff) {
          if (numCandidates >= candidates.size()) {
            candidates.push_back(neighbours[j]);
          } else {
            candidates[numCandidates] = neighbours[j];
          }
          numCandidates += 1;
        }
      }
      if (numCandidates == 0) {
        // not sure what to do in this case here...
        continue;
      }
      std::uniform_int_distribution<int> candidateDist(0, numCandidates - 1);
      const int secondPartner = candidates[candidateDist(this->e2)];
      // actually relocate both ends
      this->replaceSlipSpringPartner(
        springIdx, this->bondPartnersA[springIdx], firstPartner);
      this->replaceSlipSpringPartner(
        springIdx, this->bondPartnersB[springIdx], secondPartner);
      RUNTIME_EXP_IFN(
        this->computeBondLength(springIdx) <= this->highCutoff,
        "By shifting, newly created bonds should not be too long.");
    }
  }
  return n_accept;
}

/**
 *
 * @param springIdx
 * @param kbT
 * @return true
 * @return false
 */
bool
DPDSimulator::attemptSlipSpringShift(const size_t springIdx, const double kbT)
{
  const size_t partnerA = this->bondPartnersA[springIdx];
  const size_t partnerB = this->bondPartnersB[springIdx];
  size_t newPartnerA = partnerA;
  size_t newPartnerB = partnerB;
  // attempt to shift the spring around partnerA
  int distrLimitA = this->idxFunctionalities[partnerA] - 1;
  if (distrLimitA < 0) {
    RUNTIME_EXP_IFN(this->atomTypes[partnerA] == this->crossLinkerType,
                    "Only crosslinkers are allowed to be single beads.");
    return false;
  }
  if (distrLimitA == 0 && this->shiftPossibilityEmpty) {
    distrLimitA += 1;
  }
  std::uniform_int_distribution<int> dista(0, distrLimitA);
  size_t selectedRailBondA;
  bool shiftEndAIsFirstOnRailBond;
  const int randomIdxA = dista(this->e2);
  if (randomIdxA >= this->idxFunctionalities[partnerA]) {
    RUNTIME_EXP_IFN(this->shiftPossibilityEmpty, "Invalid state.");
    return false;
  } else {
    selectedRailBondA = this->bondsOfIndex[partnerA][randomIdxA];
    assert(this->bondPartnersA[selectedRailBondA] == partnerA ||
           this->bondPartnersB[selectedRailBondA] == partnerA);
    shiftEndAIsFirstOnRailBond =
      this->bondPartnersA[selectedRailBondA] == partnerA;
    newPartnerA = shiftEndAIsFirstOnRailBond
                    ? this->bondPartnersB[selectedRailBondA]
                    : this->bondPartnersA[selectedRailBondA];
  }

  // and around B
  int distrLimitB = this->idxFunctionalities[partnerB] - 1;
  if (distrLimitB < 0) {
    RUNTIME_EXP_IFN(this->atomTypes[partnerB] == this->crossLinkerType,
                    "Only crosslinkers are allowed to be single beads.");
    return false;
  }
  if (distrLimitB == 0 && this->shiftPossibilityEmpty) {
    distrLimitB += 1;
  }
  std::uniform_int_distribution<int> distb(0, distrLimitB);
  size_t selectedRailBondB;
  bool shiftEndBIsFirstOnRailBond;
  const int randomIdxB = distb(this->e2);
  if (randomIdxB >= this->idxFunctionalities[partnerB]) {
    RUNTIME_EXP_IFN(this->shiftPossibilityEmpty, "Invalid state.");
    return false;
  } else {
    selectedRailBondB = this->bondsOfIndex[partnerB][randomIdxB];
    assert(this->bondPartnersA[selectedRailBondB] == partnerB ||
           this->bondPartnersB[selectedRailBondB] == partnerB);
    shiftEndBIsFirstOnRailBond =
      this->bondPartnersA[selectedRailBondB] == partnerB;
    newPartnerB = shiftEndBIsFirstOnRailBond
                    ? this->bondPartnersB[selectedRailBondB]
                    : this->bondPartnersA[selectedRailBondB];
  }

  // if it's a crosslink, we don't allow shifting
  if (this->idxFunctionalities[newPartnerA] > 2 ||
      this->idxFunctionalities[newPartnerB] > 2) {
    return false;
  }

  // compute the Metropolis criterion
  Eigen::Vector3d bondDistanceNow =
    (this->coordinates.segment(partnerA * 3, 3) -
     this->coordinates.segment(partnerB * 3, 3)) +
    this->bondBoxOffsets.segment(3 * springIdx, 3);
  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(bondDistanceNow);
  }
  const double bondEnergyNow = this->k * bondDistanceNow.squaredNorm();
  const Eigen::Vector3d originalOffsets =
    this->bondBoxOffsets.segment(3 * springIdx, 3);
  const Eigen::Vector3d bondDistanceRailA =
    this->coordinates.segment(newPartnerA * 3, 3) -
    this->coordinates.segment(partnerA * 3, 3) +
    (shiftEndAIsFirstOnRailBond ? 1. : -1.) *
      this->bondBoxOffsets.segment(3 * selectedRailBondA, 3);
  const Eigen::Vector3d bondDistanceRailB =
    this->coordinates.segment(newPartnerB * 3, 3) -
    this->coordinates.segment(partnerB * 3, 3) +
    (shiftEndBIsFirstOnRailBond ? 1. : -1.) *
      this->bondBoxOffsets.segment(3 * selectedRailBondB, 3);
  Eigen::Vector3d bondDistanceNew =
    bondDistanceNow + bondDistanceRailA - bondDistanceRailB;
  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(bondDistanceNew);
  }
  const double bondEnergyNew = this->k * bondDistanceNew.squaredNorm();
  const double deltaEnergy = bondEnergyNew - bondEnergyNow;
  bool accept = false;
  if (deltaEnergy < 0.0) {
    accept = true;
  } else {
    const double factor = std::exp(-deltaEnergy / kbT);
    if (this->uniform_rand_between_0_1(this->e2) < factor) {
      accept = true;
    }
  }
  if (accept) {
    this->replaceSlipSpringPartner(springIdx, partnerA, newPartnerA);
    this->replaceSlipSpringPartner(springIdx, partnerB, newPartnerB);
#ifndef NDEBUG
    if (!((this->bondBoxOffsets.segment(3 * springIdx, 3) -
           (originalOffsets +
            ((shiftEndAIsFirstOnRailBond ? 1. : -1.) *
             this->bondBoxOffsets.segment(3 * selectedRailBondA, 3)) -
            ((shiftEndBIsFirstOnRailBond ? 1. : -1.) *
             this->bondBoxOffsets.segment(3 * selectedRailBondB, 3))))
            .maxCoeff() < 1e-10)) {
      std::cerr << "Expected bond offset is not fulfilled: got "
                << this->bondBoxOffsets.segment(3 * springIdx, 3)
                << " instead of "
                << originalOffsets +
                     ((shiftEndAIsFirstOnRailBond ? 1. : -1.) *
                      this->bondBoxOffsets.segment(3 * selectedRailBondA, 3)) -
                     ((shiftEndBIsFirstOnRailBond ? 1. : -1.) *
                      this->bondBoxOffsets.segment(3 * selectedRailBondB, 3))
                << std::endl;
    }
#endif
    this->bondBoxOffsets.segment(3 * springIdx, 3) =
      originalOffsets +
      ((shiftEndAIsFirstOnRailBond ? 1. : -1.) *
       this->bondBoxOffsets.segment(3 * selectedRailBondA, 3)) -
      ((shiftEndBIsFirstOnRailBond ? 1. : -1.) *
       this->bondBoxOffsets.segment(3 * selectedRailBondB, 3));
#ifndef NDEBUG
    if (!((this->computeBondDistance(springIdx).cwiseAbs() -
           bondDistanceNew.cwiseAbs())
            .maxCoeff() < 1e-10)) {
      std::cerr << "Expected bond distance " << bondDistanceNew << " based on "
                << bondDistanceNow << ", " << bondDistanceRailA << ", "
                << bondDistanceRailB << std::endl;
      std::cerr << "Instead, got " << this->computeBondDistance(springIdx)
                << std::endl;
      std::cerr << "Involved coordinates: "
                << this->coordinates.segment(partnerA * 3, 3) << ", "
                << this->coordinates.segment(partnerB * 3, 3) << ", and new "
                << this->coordinates.segment(newPartnerA * 3, 3) << ", "
                << this->coordinates.segment(newPartnerB * 3, 3) << " in box "
                << this->box.getL() << " with offsets " << originalOffsets
                << " on original slip-link between A & B, "
                << this->bondBoxOffsets.segment(3 * selectedRailBondA, 3)
                << " on rail A and "
                << this->bondBoxOffsets.segment(3 * selectedRailBondB, 3)
                << " on rail B. " << std::endl;
      std::cerr << "Currently stored offsets are "
                << this->bondBoxOffsets.segment(3 * springIdx, 3) << std::endl;
    }
#endif
    if (this->computeBondLength(springIdx) > this->maxBondLen) {
      std::cerr << "After shifting, managed to get too long bond with bond "
                   "energy before "
                << bondEnergyNow << " and now " << bondEnergyNew << " for bond "
                << springIdx << std::endl;
      std::cerr << "Bond length is " << this->computeBondLength(springIdx)
                << " (original: " << bondDistanceNow.norm()
                << ", now supposedly: " << bondDistanceNew.norm() << ")"
                << std::endl;
    }
  }
  return accept;
};

/**
 *
 * @param springIdx
 * @param endToShift
 * @param kbT
 * @return true
 * @return false
 */
bool
DPDSimulator::attemptSlipSpringShift(const size_t springIdx,
                                     const size_t endToShift,
                                     const double kbT)
{
  INVALIDARG_EXP_IFN(
    this->bondPartnersA[springIdx] == endToShift ||
      this->bondPartnersB[springIdx] == endToShift,
    "This spring and its partners do not match, cannot attempt shift.");
  const bool shiftEndIsFirst = this->bondPartnersA[springIdx] == endToShift;
  const size_t partnerA = shiftEndIsFirst ? this->bondPartnersA[springIdx]
                                          : this->bondPartnersB[springIdx];
  const size_t partnerB = shiftEndIsFirst ? this->bondPartnersB[springIdx]
                                          : this->bondPartnersA[springIdx];
  assert(partnerA == endToShift);
  // attempt to shift the spring around partnerA
  int distr_limit = this->idxFunctionalities[partnerA] - 1;
  if (distr_limit < 0) {
    RUNTIME_EXP_IFN(
      this->atomTypes[partnerA] == this->crossLinkerType,
      "Only crosslinkers are allowed to be single beads. Found bead " +
        std::to_string(partnerA) + " to have type " +
        std::to_string(this->atomTypes[partnerA]) + " with functionality " +
        std::to_string(distr_limit) + ".");
    return false;
  }
  RUNTIME_EXP_IFN(distr_limit >= 0,
                  "Invalid state: a slip-spring shift was attempted to "
                  "bead without enough bonds.");
  if (distr_limit == 0 && this->shiftPossibilityEmpty) {
    distr_limit += 1;
  }
  std::uniform_int_distribution<int> dist(0, distr_limit);
  const int random_idx = dist(this->e2);
  if (random_idx >= this->idxFunctionalities[partnerA]) {
    assert(this->shiftPossibilityEmpty);
    return false;
  }
  const size_t selectedRailBond = this->bondsOfIndex[partnerA][random_idx];
  const bool shiftEndIsFirstOnRailBond =
    this->bondPartnersA[selectedRailBond] == partnerA;
  const size_t replacementForA = shiftEndIsFirstOnRailBond
                                   ? this->bondPartnersB[selectedRailBond]
                                   : this->bondPartnersA[selectedRailBond];
  // if it's a crosslink, we don't allow shifting
  if (this->idxFunctionalities[replacementForA] > 2) {
    return false;
  }
  // compute the Metropolis criterion
  const Eigen::Vector3d originalOffsets =
    this->bondBoxOffsets.segment(3 * springIdx, 3);
  Eigen::Vector3d bondDistanceNow =
    (this->coordinates.segment(partnerA * 3, 3) -
     this->coordinates.segment(partnerB * 3, 3)) +
    ((shiftEndIsFirst ? 1. : -1.) *
     this->bondBoxOffsets.segment(3 * springIdx, 3));
  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(bondDistanceNow);
  }
  const double bondEnergyNow = this->k * bondDistanceNow.squaredNorm();
  const Eigen::Vector3d railDistance =
    (this->coordinates.segment(replacementForA * 3, 3) -
     this->coordinates.segment(partnerA * 3, 3)) +
    ((shiftEndIsFirstOnRailBond ? 1. : -1.) *
     this->bondBoxOffsets.segment(3 * selectedRailBond, 3));
  Eigen::Vector3d bondDistanceNew = bondDistanceNow + railDistance;
  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(bondDistanceNew);
  }
  const double bondEnergyNew = this->k * bondDistanceNew.squaredNorm();
  const double deltaEnergy = bondEnergyNew - bondEnergyNow;
  bool accept = false;
  if (deltaEnergy < 0.0) {
    accept = true;
  } else {
    const double factor = std::exp(-deltaEnergy / kbT);
    if (this->uniform_rand_between_0_1(this->e2) < factor) {
      accept = true;
    }
  }
  if (accept) {
    this->replaceSlipSpringPartner(springIdx, partnerA, replacementForA);
#ifndef NDEBUG
    if (!((this->bondBoxOffsets.segment(3 * springIdx, 3) -
           (originalOffsets +
            (shiftEndIsFirstOnRailBond ? 1. : -1.) *
              this->bondBoxOffsets.segment(3 * selectedRailBond, 3)))
            .maxCoeff() < 1e-10)) {
      std::cerr << "Expected bond offset is not fulfilled: got "
                << this->bondBoxOffsets.segment(3 * springIdx, 3)
                << " instead of "
                << (originalOffsets +
                    (shiftEndIsFirstOnRailBond ? 1. : -1.) *
                      this->bondBoxOffsets.segment(3 * selectedRailBond, 3))
                << std::endl;
    }
#endif
    this->bondBoxOffsets.segment(3 * springIdx, 3) =
      originalOffsets + (shiftEndIsFirstOnRailBond ? 1. : -1.) *
                          this->bondBoxOffsets.segment(3 * selectedRailBond, 3);
#ifndef NDEBUG
    if (!((this->computeBondDistance(springIdx).cwiseAbs() -
           bondDistanceNew.cwiseAbs())
            .maxCoeff() < 1e-10)) {
      std::cerr << "Expected bond distance " << bondDistanceNew << " based on "
                << bondDistanceNow << ", " << railDistance << "." << std::endl;
      std::cerr << "Instead, got " << this->computeBondDistance(springIdx)
                << std::endl;
      std::cerr << "Involved coordinates: "
                << this->coordinates.segment(partnerA * 3, 3) << ", "
                << this->coordinates.segment(partnerB * 3, 3) << ", and new "
                << this->coordinates.segment(replacementForA * 3, 3)
                << " in box " << this->box.getL() << " with offsets "
                << originalOffsets << " on original slip-link between A & B, "
                << this->bondBoxOffsets.segment(3 * selectedRailBond, 3)
                << " on rail A. " << std::endl;
      std::cerr << "Currently stored offsets are "
                << this->bondBoxOffsets.segment(3 * springIdx, 3) << std::endl;
    }
#endif
  }
  return accept;
};

/**
 * @brief
 *
 * @param springIdx
 * @param partnerBefore
 * @param partnerAfter
 */
void
DPDSimulator::replaceSlipSpringPartner(const size_t springIdx,
                                       const size_t partnerBefore,
                                       const size_t partnerAfter)
{
  if (partnerAfter == partnerBefore) {
    return;
  }
  INVALIDINDEX_EXP_IFN(springIdx < this->bondPartnersA.size(),
                       "Cannot replace on a non-existing spring.");
  INVALIDINDEX_EXP_IFN(partnerBefore < this->numAtoms,
                       "Cannot replace with a non-existing bead.");
  INVALIDINDEX_EXP_IFN(partnerAfter < this->numAtoms,
                       "Cannot replace with a non-existing bead.");
  INVALIDARG_EXP_IFN(this->bondPartnersA[springIdx] == partnerBefore ||
                       this->bondPartnersB[springIdx] == partnerBefore,
                     "This spring and its partners do not match.");
  INVALIDARG_EXP_IFN((this->atomTypes[partnerAfter] != this->crossLinkerType ||
                      this->idxFunctionalities[partnerAfter] <= 2),
                     "Cannot form slip-springs with crosslinkers.");
  size_t otherPartner;
  if (this->bondPartnersA[springIdx] == partnerBefore) {
    this->bondPartnersA[springIdx] = partnerAfter;
    for (int dir = 0; dir < 3; ++dir) {
      this->bondPartnerCoordinatesA[3 * springIdx + dir] =
        3 * partnerAfter + dir;
    }
    otherPartner = this->bondPartnersB[springIdx];
  } else {
    this->bondPartnersB[springIdx] = partnerAfter;
    for (int dir = 0; dir < 3; ++dir) {
      this->bondPartnerCoordinatesB[3 * springIdx + dir] =
        3 * partnerAfter + dir;
    }
    otherPartner = this->bondPartnersA[springIdx];
  }
  this->resetBondOffset(springIdx);
  // add to the bonds of the new bond partner
  pylimer_tools::utils::addToSorted<size_t>(this->bondsOfIndex[partnerAfter],
                                            springIdx);
  // remove from the bonds of the previous bond partner
  for (size_t i = this->idxFunctionalities[partnerBefore];
       i < this->bondsOfIndex[partnerBefore].size();
       ++i) {
    if (this->bondsOfIndex[partnerBefore][i] == springIdx) {
      // remove this
      this->bondsOfIndex[partnerBefore].erase(
        this->bondsOfIndex[partnerBefore].begin() + i);
      this->resetBondDuplicationPenalty(partnerAfter);
      this->resetBondDuplicationPenalty(partnerBefore);
      // this->resetBondDuplicationPenalty(otherPartner);
      return;
    }
  }
  this->validateState();
  throw std::runtime_error("Invalid internal state: replacing slip-spring "
                           "partner, but did not find it internally.");
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

  std::vector<double> xs;
  xs.reserve(this->numAtoms);
  std::vector<double> ys;
  ys.reserve(this->numAtoms);
  std::vector<double> zs;
  zs.reserve(this->numAtoms);
  std::vector<int> zeros;
  zeros.reserve(this->numAtoms);

  for (size_t i = 0; i < this->numAtoms * 3; i += 3) {
    xs.push_back(this->coordinates[i]);
    ys.push_back(this->coordinates[i + 1]);
    zs.push_back(this->coordinates[i + 2]);
    zeros.push_back(0);
  }

  result.addAtoms(
    this->atomIds, this->atomTypes, xs, ys, zs, zeros, zeros, zeros);

  const size_t numBondsToAdd =
    withSlipsprings ? (this->numBonds + this->numSlipSprings) : this->numBonds;
  std::vector<long int> bondFrom;
  bondFrom.reserve(numBondsToAdd);
  std::vector<long int> bondTo;
  bondTo.reserve(numBondsToAdd);
  std::vector<int> newBondTypes;
  newBondTypes.reserve(numBondsToAdd);

  for (size_t i = 0; i < numBondsToAdd; ++i) {
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
  for (size_t i = 0; i < this->numAtoms; ++i) {
    const int numNeighbors = this->neighbourlist.getIndicesCloseToCoordinates(
      neighbors, this->coordinates.segment(3 * i, 3), cutoff);
    Eigen::ArrayXi neighbors2 =
      this->neighbourlist.getIndicesCloseToCoordinates(
        this->coordinates.segment(3 * i, 3), cutoff);
    RUNTIME_EXP_IFN(
      (neighbors.segment(0, numNeighbors) == neighbors2).all(),
      "Neighbors should be equal no matter the method, but apparently, "
      "they are not.");

    std::vector<size_t> relevantNeighbors;
    std::vector<size_t> relevantPairs;

    // neighbour-list
    for (size_t neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
      const size_t j = neighbors[neigh_idx];
      if (j <= i) {
        continue;
      }
      Eigen::Vector3d pairDistance = this->coordinates.segment(3 * i, 3) -
                                     this->coordinates.segment(3 * j, 3);
      this->box.handlePBC(pairDistance);
      const double rNorm = pairDistance.norm();
      if (rNorm >= cutoff || rNorm < 1e-12) {
        continue;
      }

      relevantNeighbors.push_back(j);
    }

    assert(relevantNeighbors.size() <= numNeighbors);

    // pairs
    for (size_t j = i + 1; j < this->numAtoms; ++j) {
      Eigen::Vector3d pairDistance = this->coordinates.segment(3 * i, 3) -
                                     this->coordinates.segment(3 * j, 3);
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
          this->coordinates.segment(3 * i, 3),
          this->coordinates.segment(3 * diff_j, 3),
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

  Eigen::VectorXd bondLengths =
    Eigen::VectorXd::Zero(this->numBonds + this->numSlipSprings);
#pragma omp parallel for
  for (size_t i = 0; i < (this->numBonds + this->numSlipSprings); ++i) {
    const double b = bondDistances.segment(3 * i, 3).norm();
    bondLengths[i] = b;
  }
  return bondLengths;
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
  RUNTIME_EXP_IFN(this->atomTypes.size() == this->numAtoms,
                  "State violation: size of atom types incorrect.");
  RUNTIME_EXP_IFN(this->atomIds.size() == this->numAtoms,
                  "State violation: size of atom ids incorrect.");
  RUNTIME_EXP_IFN(this->bondsOfIndex.size() == this->numAtoms,
                  "State violation: bonds of indices distributed incorrectly.");
  for (size_t i = 0; i < this->numAtoms; ++i) {
    const int num_actual_bonds =
      std::accumulate(this->bondsOfIndex[i].begin(),
                      this->bondsOfIndex[i].end(),
                      0,
                      [&](const int val, const size_t bondIdx) {
                        return val + static_cast<int>(bondIdx < this->numBonds);
                      });
    if (this->numSlipSprings == 0) {
      RUNTIME_EXP_IFN(num_actual_bonds == this->bondsOfIndex[i].size(),
                      "State violation: bonds of index " + std::to_string(i) +
                        " were not counted appropriately.");
    }
    RUNTIME_EXP_IFN(this->idxFunctionalities[i] == num_actual_bonds,
                    "State violation: inconsistent idx functionalities");
    RUNTIME_EXP_IFN(
      this->idxFunctionalities[i] >= 1 ||
        this->atomTypes[i] == this->crossLinkerType,
      "Only crosslinkers are allowed to be single beads. Found bead " +
        std::to_string(i) + " to have type " +
        std::to_string(this->atomTypes[i]) + " with functionality " +
        std::to_string(this->idxFunctionalities[i]) + ".");
    if (this->idxFunctionalities[i] == 1) {
      RUNTIME_EXP_IFN(
        this->isRelocationTarget[i],
        "Expected functionalities of 1 to be relocation targets.");
    }
    for (size_t j = 0; j < this->bondsOfIndex[i].size(); ++j) {
      if (j < this->idxFunctionalities[i]) {
        RUNTIME_EXP_IFN(this->bondsOfIndex[i][j] < this->numBonds,
                        "Expect bonds to come before slip springs.");
      } else {
        RUNTIME_EXP_IFN(this->bondsOfIndex[i][j] >= this->numBonds,
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
                    this->chainEndIndices.size(),
                  "Expect relocation targets to be equal to chain ends.");
  for (const size_t i : this->chainEndIndices) {
    RUNTIME_EXP_IFN(this->isRelocationTarget[i],
                    "Expect chain ends to be relocation targets.");
    RUNTIME_EXP_IFN(this->idxFunctionalities[i] <= 1 ||
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
  for (size_t i = 0; i < this->numAtoms; ++i) {
    partners.clear();
    // here as well, we rely on the bonds of index being sorted
    for (const size_t bondIdx : this->bondsOfIndex[i]) {
      size_t atomPartnerIdx = this->bondPartnersA[bondIdx] == i
                                ? this->bondPartnersB[bondIdx]
                                : this->bondPartnersA[bondIdx];
      bool expectZero = false;
      if (partners.contains(atomPartnerIdx)) {
        // "real" bonds always contribute -> check that this is a
        // slip-spring
        if (bondIdx >= this->numBonds) {
          expectZero = true;
        }
      } else {
        partners.insert(atomPartnerIdx);
      }
      RUNTIME_EXP_IFN(
        this->bondDuplicationPenalty.segment(3 * bondIdx, 3)
          .isApprox(expectZero ? Eigen::Array3d::Zero()
                               : Eigen::Array3d::Constant(1.)),
        "Incorrect bondDuplicationPenalty: found bond " +
          std::to_string(bondIdx) + " of " + std::to_string(this->numBonds) +
          " (+ " + std::to_string(this->numSlipSprings) + ") to should have " +
          std::to_string(expectZero ? 0 : 1) + " but it did not (" +
          std::to_string(
            this->bondDuplicationPenalty.segment(3 * bondIdx, 3).sum()) +
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
  for (size_t i = 0; i < this->numBonds + this->numSlipSprings; ++i) {
    for (int dir = 0; dir < 3; ++dir) {
      RUNTIME_EXP_IFN(
        this->bondPartnerCoordinatesA[i * 3 + dir] ==
          this->bondPartnersA[i] * 3 + dir,
        "Bond partners not accurate: got partner coordinate idx " +
          std::to_string(this->bondPartnerCoordinatesA[i * 3 + dir]) +
          " but expected " + std::to_string(this->bondPartnersA[i] * 3 + dir) +
          " at index " + std::to_string(i) + " and dir " + std::to_string(dir) +
          ".");
      RUNTIME_EXP_IFN(
        this->bondPartnerCoordinatesB[i * 3 + dir] ==
          this->bondPartnersB[i] * 3 + dir,
        "Bond partners not accurate: got partner coordinate idx " +
          std::to_string(this->bondPartnerCoordinatesB[i * 3 + dir]) +
          " but expected " + std::to_string(this->bondPartnersB[i] * 3 + dir) +
          " at index " + std::to_string(i) + " and dir " + std::to_string(dir) +
          ".");
    }
    RUNTIME_EXP_IFN(pylimer_tools::utils::contains(
                      this->bondsOfIndex[this->bondPartnersA[i]], i),
                    "Reverse-link is incorrect. Bond " + std::to_string(i) +
                      " (of " + std::to_string(this->numBonds) +
                      ") was not found in bonds of partner A, " +
                      std::to_string(this->bondPartnersA[i]) + ".");
    RUNTIME_EXP_IFN(pylimer_tools::utils::contains(
                      this->bondsOfIndex[this->bondPartnersB[i]], i),
                    "Reverse-link is incorrect. Bond " + std::to_string(i) +
                      " (of " + std::to_string(this->numBonds) +
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
