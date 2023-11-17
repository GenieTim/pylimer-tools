#include "./DPDSimulator.h"
#include "../utils/PerformanceTimer.h"
#include "./Correlator.h"

#include <fstream>
#include <iostream>
#include <random>

namespace pylimer_tools {
namespace calc {
  namespace dpd {

    DPDSimulator::DPDSimulator(const pylimer_tools::entities::Universe& u,
                               const int crosslinkerType,
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
      // initialize the random number generator
      if (seed == "") {
        std::random_device rd;
        this->e2 = std::mt19937(rd());
      } else {
        std::seed_seq seed2(seed.begin(), seed.end());
        this->e2 = std::mt19937(seed2);
      }
      double mean = 0.0;
      double std = 1.0;
      double a = mean - std::sqrt(3.) * std;
      double b = mean + std::sqrt(3.) * std;
      this->uniform_rand_mean0std1 =
        std::uniform_real_distribution<double>(a, b);
      this->uniform_rand_between_0_1 =
        std::uniform_real_distribution<double>(0., 1.);

      // initialize the faster data structure
      this->box = u.getBox();
      this->coordinates = u.getUnwrappedVertexCoordinates(&this->box);
      std::map<std::string, std::vector<long int>> edges = u.getEdges();
      this->bondPartnersA =
        Eigen::Map<Eigen::ArrayXli, Eigen::Unaligned>(edges["edge_from"].data(),
                                                      edges["edge_from"].size())
          .cast<int>();
      this->bondPartnersB = Eigen::Map<Eigen::ArrayXli, Eigen::Unaligned>(
                              edges["edge_to"].data(), edges["edge_to"].size())
                              .cast<int>();
      this->bondTypes = Eigen::Map<Eigen::ArrayXli, Eigen::Unaligned>(
                          edges["edge_type"].data(), edges["edge_type"].size())
                          .cast<int>();

      this->neighbourlist = pylimer_tools::entities::EigenNeighbourList(
        coordinates, this->box, 1.0);
      this->numAtoms = this->coordinates.size() / 3;
      this->numBonds = this->bondPartnersA.size();
      this->idxFunctionalities = Eigen::ArrayXi::Zero(this->numAtoms);
      this->atomTypes = u.getPropertyValues<int>("type");
      this->atomIds = u.getPropertyValues<long int>("id");
      this->maxBondLen = 0.45 * this->box.getL().maxCoeff();

      this->bondsOfIndex.reserve(this->numAtoms);
      for (size_t i = 0; i < this->numAtoms; ++i) {
        std::vector<size_t> bonds;
        bonds.reserve(4);
        this->bondsOfIndex.push_back(bonds);
      }
      this->bondPartnerCoordinatesA = Eigen::ArrayXi(3 * this->numBonds);
      this->bondPartnerCoordinatesB = Eigen::ArrayXi(3 * this->numBonds);
      for (size_t i = 0; i < this->numBonds; ++i) {
        this->bondsOfIndex[this->bondPartnersA[i]].push_back(i);
        this->bondsOfIndex[this->bondPartnersB[i]].push_back(i);
        for (int dir = 0; dir < 3; ++dir) {
          this->bondPartnerCoordinatesA[i * 3 + dir] =
            this->bondPartnersA[i] * 3 + dir;
          this->bondPartnerCoordinatesB[i * 3 + dir] =
            this->bondPartnersB[i] * 3 + dir;
        }
      }
      for (size_t i = 0; i < this->numAtoms; ++i) {
        this->idxFunctionalities[i] = this->bondsOfIndex[i].size();
        if (this->idxFunctionalities[i] < 2 ||
            this->atomTypes[i] == crosslinkerType) {
          this->chainEndIndices.push_back(i);
        }
      }

      // simulation state
      this->currentVelocitiesPlus = Eigen::VectorXd::Zero(coordinates.size());
      this->currentVelocities = Eigen::VectorXd::Zero(coordinates.size());
      this->currentForces = Eigen::VectorXd::Zero(coordinates.size());
      this->currentStressTensor = Eigen::Matrix3d::Zero();
    }

    /**
     * @brief
     *
     * @param nSteps
     * @param dt
     * @param lambda
     * @param withMC
     */
    void DPDSimulator::runSimulation(
      const long int nSteps,
      bool withMC,
      const std::function<bool()>& shouldInterrupt,
      const std::function<void()>& cleanupInterrupt)
    {
      bool wasInterrupted = false;

      const double halfDt = 0.5 * this->dt;
      double temperature = this->computeTemperature(this->currentVelocities);

      this->prepareAllOutputs();

      pylimer_tools::utils::PerformanceTimer timer =
        pylimer_tools::utils::PerformanceTimer<
          DPDPerformanceSections::NUM_PERFORMANCE_SECTIONS>();
      timer.registerSections(DPDPerformanceSectionNames);

      // start iterating over the steps to do
      long int step = 0;
      for (; step < nSteps; step++) {
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
        this->currentVelocities += lambda * this->dt * this->currentForces;

        // TODO: figure out, what the issue is, how the velocities
        // should be synchronized for the pressure
        // temperature = this->computeTemperature(velocities);

        // re-compute the forces with these updated coordinates & velocities
        timer.section(DPDPerformanceSections::FORCES);
        double pressure = computeForces(this->currentForces,
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
        double kineticPressureTerm =
          ((this->numAtoms * temperature) / this->box.getVolume());
        const double m = 1.;
        for (size_t i = 0; i < this->numAtoms; ++i) {
          this->currentStressTensor -=
            m * this->currentVelocities.segment(3 * i, 3) *
            this->currentVelocities.segment(3 * i, 3).transpose();
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
      std::ios::sync_with_stdio(true);
      this->outputStreams.clear();
      this->outputFileStreams.clear();

      if (wasInterrupted) {
        cleanupInterrupt();
      }

      timer.stop();
      timer.output();
    }

    /**
     * @brief Determine the temperature of the system
     *
     * @param velocities
     * @return double
     */
    double DPDSimulator::computeTemperature(
      const Eigen::VectorXd& velocities) const
    {
      // configuration
      const double dim = 3;
      const double kb = 1.0;
      const double m = 1.0;

      double KE = 0.5 * m * velocities.squaredNorm();

      return KE / ((dim / 2.) * (velocities.size() / dim) * kb);
    }

    /**
     * @brief Compute the force vector, and return the pressure
     *
     * @param forces
     * @param coordinates
     * @param velocities
     * @param dt
     * @param cutoff
     * @param A
     * @param sigma
     * @param k
     * @return double
     */
    double DPDSimulator::computeForces(
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
      timer.section(DPDPerformanceSections::BOND_FORCE);
      Eigen::VectorXd bondDistances = coords(this->bondPartnerCoordinatesA) -
                                      coords(this->bondPartnerCoordinatesB);
      this->box.handlePBC(bondDistances);
      assert(bondDistances.minCoeff() > -this->box.getL().maxCoeff());
      assert(bondDistances.maxCoeff() < this->box.getL().maxCoeff());
      forces(this->bondPartnerCoordinatesA) -= this->k * bondDistances;
      forces(this->bondPartnerCoordinatesB) += this->k * bondDistances;

#pragma omp parallel for reduction(+ : stressTensor, pressure)
      for (size_t i = 0; i < this->bondPartnersA.size(); ++i) {
        // TODO: check sign
        pressure += -this->k * bondDistances.segment(3 * i, 3).squaredNorm();
        // pressure -= 0.5 * this->k * bondDistances.segment(3 * i,
        // 3).squaredNorm();
        // pressure += 0.5 * this->k * bondDistances.segment(3
        // * i, 3).squaredNorm();
        stressTensor += this->k * bondDistances.segment(3 * i, 3) *
                        bondDistances.segment(3 * i, 3).transpose();
      }

      timer.section(DPDPerformanceSections::PAIR_FORCE);

      // pre-allocate the neighbor indices array
      Eigen::ArrayXi neighbors = Eigen::ArrayXi(static_cast<int>(
        this->numAtoms *
        (std::ceil((3.1 * cutoff) * (3.1 * cutoff) * (3.1 * cutoff)) /
         this->box.getVolume())));

      // actually loop the atoms
      for (size_t i = 0; i < this->numAtoms; ++i) {
        int numNeighbors = this->neighbourlist.getIndicesCloseToCoordinates(
          neighbors, coords.segment(3 * i, 3), cutoff);

        // pair forces
        // for (size_t j = i + 1; j < this->numAtoms; ++j) {
        for (size_t neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
          const size_t j = neighbors[neigh_idx];
          if (j <= i) {
            continue;
          }
          Eigen::Vector3d pairdistance =
            coords.segment(3 * i, 3) - coords.segment(3 * j, 3);
          this->box.handlePBC(pairdistance);
          const double rNorm = pairdistance.norm();
          if (rNorm >= cutoff || rNorm < 1e-12) {
            continue;
          }

          const double one_minus_rnorm = 1. - rNorm;
          const double one_minus_rnorm2 = (1. - rNorm) * (1. - rNorm);
          Eigen::Vector3d pairdistanceNormed = pairdistance / rNorm;

          // conservative repulsion force
          double pairForceConst = this->A * one_minus_rnorm;
          // pairForce = this->A * one_minus_rnorm * pairdistanceNormed;

          // dissipative/drag force
          Eigen::Vector3d velocitydiff =
            velocities.segment(3 * i, 3) - velocities.segment(3 * j, 3);
          const double rij_dot_vij = pairdistanceNormed.dot(velocitydiff);
          const double gamma_weighted_rij_dot_vij =
            this->gamma * one_minus_rnorm2 * rij_dot_vij;

          pairForceConst -= gamma_weighted_rij_dot_vij;
          // pairForce += -gamma_weighted_rij_dot_vij * pairdistanceNormed;

          // random force
          const double constant_rnd_prefix =
            this->sigma * one_minus_rnorm / sqrt(dt);
          const double random_val = this->uniform_rand_mean0std1(this->e2);

          // pairForce += constant_rnd_prefix * random_val * pairdistanceNormed;
          pairForceConst += constant_rnd_prefix * random_val;
          Eigen::Vector3d pairForce = pairForceConst * pairdistanceNormed;

          // actually assign the new forces
          forces.segment(3 * i, 3) += pairForce;
          forces.segment(3 * j, 3) -= pairForce;

          // pressure update
          pressure += pairForce.dot(pairdistance);
          stressTensor -= pairForce * pairdistance.transpose();
        }
      }

      return pressure / (3. * this->box.getVolume());
    }

    /**
     * @brief Get access to the current stress-tensor
     *
     * @return Eigen::Matrix3d
     */
    Eigen::Matrix3d DPDSimulator::getStressTensor() const
    {
      return this->currentStressTensor;
    }

    /**
     * @brief Register a set of atoms and time for being measured for msd
     *
     * @param atomIds
     */
    void DPDSimulator::startMeasuringMSDForAtoms(
      const std::vector<size_t>& atomIdsToMeasure)
    {
      // Translate atom IDS to indices of the local structure
      Eigen::ArrayXi coordinateIndices =
        Eigen::ArrayXi(3 * atomIdsToMeasure.size());
#pragma omp parallel for
      for (size_t i = 0; i < atomIdsToMeasure.size(); ++i) {
        size_t atomId = atomIdsToMeasure[i];
        size_t index = this->universe.getIdxByAtomId(atomId);
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
     * @param loCutoff
     * @param hiCutoff
     * @return int
     */
    int DPDSimulator::createSlipSprings(const int num, const int bondType)
    {
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
        sourceIds.push_back(i);
      }

      // search for neighbours that are elibile
      Eigen::ArrayXi neighbours = Eigen::ArrayXi(16);
      while (createdLastIteration > 0 && totalCreated < num) {
        createdLastIteration = 0;
        std::shuffle(sourceIds.begin(), sourceIds.end(), this->e2);
        for (size_t i : sourceIds) {
          int numCandidates = 0;
          // for each atom, search for possible partners
          int numNeighs = this->neighbourlist.getIndicesCloseToCoordinates(
            neighbours, this->coordinates.segment(3 * i, 3), this->highCutoff);
          for (size_t j = 0; j < numNeighs; ++j) {
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
          int candidateIndex = dist(this->e2);
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

      this->addSlipSprings(slipSpringFrom, slipSpringTo, bondType);
      return totalCreated;
    }

    /**
     * @brief Add slip-springs according to the specification
     *
     * @param partnerA
     * @param partnerB
     * @param bondType
     */
    void DPDSimulator::addSlipSprings(std::vector<size_t>& partnerA,
                                      std::vector<size_t>& partnerB,
                                      const int bondType)
    {
      INVALIDARG_EXP_IFN(partnerA.size() == partnerB.size(),
                         "Require same size A & B");
      for (size_t i = 0; i < partnerA.size(); ++i) {
        INVALIDARG_EXP_IFN(partnerA[i] < this->numAtoms, "Invalid partner id");
        INVALIDARG_EXP_IFN(partnerB[i] < this->numAtoms, "Invalid partner id");
      }
      size_t sizeBefore = this->numBonds + this->numSlipSprings;
      this->bondPartnersA.conservativeResize(sizeBefore + partnerA.size());
      this->bondPartnersB.conservativeResize(sizeBefore + partnerB.size());
      this->bondPartnerCoordinatesA.conservativeResize(
        3 * (sizeBefore + partnerA.size()));
      this->bondPartnerCoordinatesB.conservativeResize(
        3 * (sizeBefore + partnerB.size()));
      this->bondTypes.conservativeResize(sizeBefore + partnerB.size());

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
      }

      this->numSlipSprings += partnerA.size();
    }

    /**
     * @brief
     *
     * @param kbT
     * @param k
     * @return int
     */
    int DPDSimulator::relocateSlipSprings(const double kbT)
    {
      int nAccept = 0;
      std::vector<size_t> candidates;
      Eigen::ArrayXi neighbours = Eigen::ArrayXi(12);
      std::uniform_int_distribution<int> chainendDist(
        0, this->chainEndIndices.size() - 1);
      for (size_t springIdx = this->numBonds;
           springIdx < (this->numBonds + this->numSlipSprings);
           ++springIdx) {
        if (!(this->idxFunctionalities[this->bondPartnersA[springIdx]] < 2 ||
              this->idxFunctionalities[this->bondPartnersB[springIdx]] < 2)) {
          continue;
        }
        // design & attempt move
        const size_t partnerA = this->bondPartnersA[springIdx];
        const size_t partnerB = this->bondPartnersB[springIdx];

        int candidateIndex = this->chainEndIndices[chainendDist(this->e2)];

        // search for neighbours
        int numCandidates = 0;
        // for each atom, search for possible partners
        int numNeighs = this->neighbourlist.getIndicesCloseToCoordinates(
          neighbours,
          this->coordinates.segment(3 * candidateIndex, 3),
          this->highCutoff);
        for (size_t j = 0; j < numNeighs; ++j) {
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
        int candidatePartnerIndex = candidates[candidateDist(this->e2)];

        // compute the Metropolis criterion
        Eigen::Vector3d bondDistanceNow =
          this->coordinates.segment(partnerA * 3, 3) -
          this->coordinates.segment(partnerB * 3, 3);
        this->box.handlePBC(bondDistanceNow);
        double bondEnergyNow = k * bondDistanceNow.squaredNorm();
        Eigen::Vector3d bondDistanceNew =
          this->coordinates.segment(candidateIndex * 3, 3) -
          this->coordinates.segment(candidatePartnerIndex * 3, 3);
        this->box.handlePBC(bondDistanceNew);
        double bondEnergyNew = k * bondDistanceNew.squaredNorm();
        double deltaEnergy = bondEnergyNew - bondEnergyNow;
        bool accept = false;
        if (deltaEnergy < 0.0) {
          accept = true;
        } else {
          double factor = std::exp(-deltaEnergy / kbT);
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
     * @param k
     * @return int
     */
    int DPDSimulator::shiftSlipSprings(const double kbT)
    {
      int n_accept = 0;
      std::vector<size_t> candidates;
      Eigen::ArrayXi neighbours = Eigen::ArrayXi(12);

      std::uniform_int_distribution<int> uniformDistNatoms(0,
                                                           this->numAtoms - 1);
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
          // search for neighbours
          int numCandidates = 0;
          // for this first chosen atom, search for possible partners
          int numNeighs = this->neighbourlist.getIndicesCloseToCoordinates(
            neighbours,
            this->coordinates.segment(3 * firstPartner, 3),
            this->highCutoff);
          RUNTIME_EXP_IFN(numNeighs <= neighbours.size(),
                          "Neighbourlist does not act as it should.");
          for (size_t j = 0; j < numNeighs; ++j) {
            RUNTIME_EXP_IFN(neighbours[j] < this->numAtoms,
                            "Neighbourlist seems inconsequent.");
            Eigen::Vector3d distance =
              this->coordinates.segment(3 * firstPartner, 3) -
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
            // not sure what to do in this case here...
            continue;
          }
          std::uniform_int_distribution<int> candidateDist(0,
                                                           numCandidates - 1);
          int secondPartner = candidates[candidateDist(this->e2)];
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
    bool DPDSimulator::attemptSlipSpringShift(const size_t springIdx,
                                              const double kbT)
    {
      const size_t partnerA = this->bondPartnersA[springIdx];
      const size_t partnerB = this->bondPartnersB[springIdx];
      size_t newPartnerA = partnerA;
      size_t newPartnerB = partnerB;
      // attempt to shift the spring around partnerA
      int distrLimitA = this->idxFunctionalities[partnerA] - 1;
      if (distrLimitA == 0 && this->shiftPossibilityEmpty) {
        distrLimitA += 1;
      }
      std::uniform_int_distribution<int> dista(0, distrLimitA);
      int randomIdxA = dista(this->e2);
      if (randomIdxA >= this->idxFunctionalities[partnerA]) {
        RUNTIME_EXP_IFN(this->shiftPossibilityEmpty, "Invalid state.");
      } else {
        const size_t selectedBondA = this->bondsOfIndex[partnerA][randomIdxA];
        newPartnerA = this->bondPartnersA[selectedBondA] == partnerA
                        ? this->bondPartnersB[selectedBondA]
                        : this->bondPartnersA[selectedBondA];
      }

      // and around B
      int distrLimitB = this->idxFunctionalities[partnerB] - 1;
      if (distrLimitB == 0 && this->shiftPossibilityEmpty) {
        distrLimitB += 1;
      }
      std::uniform_int_distribution<int> distb(0, distrLimitB);
      int randomIdxB = distb(this->e2);
      if (randomIdxB >= this->idxFunctionalities[partnerB]) {
        assert(this->shiftPossibilityEmpty);
      } else {
        const size_t selectedBondB = this->bondsOfIndex[partnerB][randomIdxB];
        newPartnerB = this->bondPartnersA[selectedBondB] == partnerB
                        ? this->bondPartnersB[selectedBondB]
                        : this->bondPartnersA[selectedBondB];
      }

      // compute the Metropolis criterion
      Eigen::Vector3d bondDistanceNow =
        (this->coordinates.segment(partnerA * 3, 3) -
         this->coordinates.segment(partnerB * 3, 3));
      this->box.handlePBC(bondDistanceNow);
      double bondEnergyNow = this->k * bondDistanceNow.squaredNorm();
      Eigen::Vector3d bondDistanceNew =
        (this->coordinates.segment(newPartnerA * 3, 3) -
         this->coordinates.segment(newPartnerB * 3, 3));
      this->box.handlePBC(bondDistanceNew);
      double bondEnergyNew = this->k * bondDistanceNew.squaredNorm();
      double deltaEnergy = bondEnergyNew - bondEnergyNow;
      bool accept = false;
      if (deltaEnergy < 0.0) {
        accept = true;
      } else {
        double factor = std::exp(-deltaEnergy / kbT);
        if (this->uniform_rand_between_0_1(this->e2) < factor) {
          accept = true;
        }
      }
      if (accept) {
        this->replaceSlipSpringPartner(springIdx, partnerA, newPartnerA);
        this->replaceSlipSpringPartner(springIdx, partnerB, newPartnerB);
        if (this->computeBondLength(springIdx) > this->maxBondLen) {
          std::cerr << "After shifting, managed to get too long bond with bond "
                       "energy before "
                    << bondEnergyNow << " and now " << bondEnergyNew
                    << " for bond " << springIdx << std::endl;
          std::cerr << "Bond length is " << this->computeBondLength(springIdx)
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
    bool DPDSimulator::attemptSlipSpringShift(const size_t springIdx,
                                              const size_t endToShift,
                                              const double kbT)
    {
      INVALIDARG_EXP_IFN(
        this->bondPartnersA[springIdx] == endToShift ||
          this->bondPartnersB[springIdx] == endToShift,
        "This spring and its partners do not match, cannot attempt shift.");
      const size_t partnerA = this->bondPartnersA[springIdx] == endToShift
                                ? this->bondPartnersA[springIdx]
                                : this->bondPartnersB[springIdx];
      const size_t partnerB = this->bondPartnersB[springIdx] == endToShift
                                ? this->bondPartnersA[springIdx]
                                : this->bondPartnersB[springIdx];
      assert(partnerA == endToShift);
      // attempt to shift the spring around partnerA
      int distr_limit = this->idxFunctionalities[partnerA] - 1;
      if (distr_limit == 0 && this->shiftPossibilityEmpty) {
        distr_limit += 1;
      }
      std::uniform_int_distribution<int> dist(0, distr_limit);
      int random_idx = dist(this->e2);
      if (random_idx >= this->idxFunctionalities[partnerA]) {
        assert(this->shiftPossibilityEmpty);
        return false;
      }
      const size_t selectedBond = this->bondsOfIndex[partnerA][random_idx];
      const size_t replacementForA =
        this->bondPartnersA[selectedBond] == partnerA
          ? this->bondPartnersB[selectedBond]
          : this->bondPartnersA[selectedBond];
      // compute the Metropolis criterion
      Eigen::Vector3d bondDistanceNow =
        (this->coordinates.segment(partnerA * 3, 3) -
         this->coordinates.segment(partnerB * 3, 3));
      this->box.handlePBC(bondDistanceNow);
      double bondEnergyNow = this->k * bondDistanceNow.squaredNorm();
      Eigen::Vector3d bondDistanceNew =
        (this->coordinates.segment(replacementForA * 3, 3) -
         this->coordinates.segment(partnerB * 3, 3));
      this->box.handlePBC(bondDistanceNew);
      double bondEnergyNew = this->k * bondDistanceNew.squaredNorm();
      double deltaEnergy = bondEnergyNew - bondEnergyNow;
      bool accept = false;
      if (deltaEnergy < 0.0) {
        accept = true;
      } else {
        double factor = std::exp(-deltaEnergy / kbT);
        if (this->uniform_rand_between_0_1(this->e2) < factor) {
          accept = true;
        }
      }
      if (accept) {
        this->replaceSlipSpringPartner(springIdx, partnerA, replacementForA);
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
    void DPDSimulator::replaceSlipSpringPartner(const size_t springIdx,
                                                const size_t partnerBefore,
                                                const size_t partnerAfter)
    {
      INVALIDARG_EXP_IFN(springIdx < this->bondPartnersA.size(),
                         "Cannot replace on a non-existing spring.");
      INVALIDARG_EXP_IFN(partnerBefore < this->numAtoms,
                         "Cannot replace with a non-existing bead.");
      INVALIDARG_EXP_IFN(partnerAfter < this->numAtoms,
                         "Cannot replace with a non-existing bead.");
      INVALIDARG_EXP_IFN(this->bondPartnersA[springIdx] == partnerBefore ||
                           this->bondPartnersB[springIdx] == partnerBefore,
                         "This spring and its partners do not match.");
      if (this->bondPartnersA[springIdx] == partnerBefore) {
        this->bondPartnersA[springIdx] = partnerAfter;
        for (int dir = 0; dir < 3; ++dir) {
          this->bondPartnerCoordinatesA[3 * springIdx + dir] =
            3 * partnerAfter + dir;
        }
      } else {
        this->bondPartnersB[springIdx] = partnerAfter;
        for (int dir = 0; dir < 3; ++dir) {
          this->bondPartnerCoordinatesB[3 * springIdx + dir] =
            3 * partnerAfter + dir;
        }
      }
      // add to the bonds of the new bond partner
      this->bondsOfIndex[partnerAfter].push_back(springIdx);
      // remove from the bonds of the previous bond partner
      for (size_t i = this->idxFunctionalities[partnerBefore];
           i < this->bondsOfIndex[partnerBefore].size();
           ++i) {
        if (this->bondsOfIndex[partnerBefore][i] == springIdx) {
          // remove this
          this->bondsOfIndex[partnerBefore].erase(
            this->bondsOfIndex[partnerBefore].begin() + i);
          return;
        }
      }
      throw std::runtime_error("Invalid internal state: replacing slip-spring "
                               "partner, but did not find it internally.");
    }

    /**
     * @brief Convert the current structure to a Unvierse
     *
     * @return pylimer_tools::entities::Universe
     */
    pylimer_tools::entities::Universe DPDSimulator::getUniverse() const
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

      std::vector<long int> bondFrom(this->bondPartnersA.data(),
                                     this->bondPartnersA.data() +
                                       this->bondPartnersA.size());
      std::vector<long int> bondTo(this->bondPartnersB.data(),
                                   this->bondPartnersB.data() +
                                     this->bondPartnersB.size());
      for (size_t i = 0; i < this->bondPartnersB.size(); ++i) {
        bondTo[i] = this->universe.getAtomIdByIdx(bondTo[i]);
        bondFrom[i] = this->universe.getAtomIdByIdx(bondFrom[i]);
      }
      std::vector<int> newBondTypes(this->bondTypes.data(),
                                    this->bondTypes.data() +
                                      this->bondTypes.size());
      result.addBonds(this->numBonds + this->numSlipSprings,
                      bondFrom,
                      bondTo,
                      newBondTypes,
                      false,
                      false);

      return result;
    }

    void DPDSimulator::validateNeighbourlist(double cutoff)
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
        int numNeighbors = this->neighbourlist.getIndicesCloseToCoordinates(
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

        // neigbhourlist
        for (size_t neigh_idx = 0; neigh_idx < numNeighbors; ++neigh_idx) {
          const size_t j = neighbors[neigh_idx];
          if (j <= i) {
            continue;
          }
          Eigen::Vector3d pairdistance = this->coordinates.segment(3 * i, 3) -
                                         this->coordinates.segment(3 * j, 3);
          this->box.handlePBC(pairdistance);
          const double rNorm = pairdistance.norm();
          if (rNorm >= cutoff || rNorm < 1e-12) {
            continue;
          }

          relevantNeighbors.push_back(j);
        }

        assert(relevantNeighbors.size() <= numNeighbors);

        // pairs
        for (size_t j = i + 1; j < this->numAtoms; ++j) {
          Eigen::Vector3d pairdistance = this->coordinates.segment(3 * i, 3) -
                                         this->coordinates.segment(3 * j, 3);
          this->box.handlePBC(pairdistance);
          const double rNorm = pairdistance.norm();
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

        std::sort(relevantPairs.begin(), relevantPairs.end());
        std::sort(relevantNeighbors.begin(), relevantNeighbors.end());
        if (relevantPairs.size() > relevantNeighbors.size()) {
          std::cout << "Debugging neighbourlist. " << numNeighbors << std::endl;
          // debug why
          // find the difference
          std::vector<size_t> diff;
          std::set_difference(relevantPairs.begin(),
                              relevantPairs.end(),
                              relevantNeighbors.begin(),
                              relevantNeighbors.end(),
                              std::back_inserter(diff));
          assert(diff.size() >=
                 (relevantPairs.size() - relevantNeighbors.size()));
          // figure out why not included
          for (size_t diff_j : diff) {
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
     * @brief Make sure all the structures obey the expected form
     *
     */
    void DPDSimulator::validateState()
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
      RUNTIME_EXP_IFN(
        this->bondsOfIndex.size() == this->numAtoms,
        "State violation: bonds of indices distributed incorrectly.");
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
      // internal structure
      RUNTIME_EXP_IFN((this->bondPartnersB.array() < this->numAtoms).all(),
                      "State violation: too large indices found (e.g. " +
                        std::to_string(this->bondPartnersB.maxCoeff()) +
                        " for " + std::to_string(this->numAtoms) + " atoms)");
      RUNTIME_EXP_IFN((this->bondPartnersA.array() < this->numAtoms).all(),
                      "State violation: too large indices found (e.g. " +
                        std::to_string(this->bondPartnersA.maxCoeff()) +
                        " for " + std::to_string(this->numAtoms) + " atoms)");
      for (size_t i = 0; i < this->numBonds + this->numSlipSprings; ++i) {
        for (int dir = 0; dir < 3; ++dir) {
          RUNTIME_EXP_IFN(
            this->bondPartnerCoordinatesA[i * 3 + dir] ==
              this->bondPartnersA[i] * 3 + dir,
            "Bond partners not accurate: got partner coordinate idx " +
              std::to_string(this->bondPartnerCoordinatesA[i * 3 + dir]) +
              " but expected " +
              std::to_string(this->bondPartnersA[i] * 3 + dir) + " at index " +
              std::to_string(i) + " and dir " + std::to_string(dir) + ".");
          RUNTIME_EXP_IFN(
            this->bondPartnerCoordinatesB[i * 3 + dir] ==
              this->bondPartnersB[i] * 3 + dir,
            "Bond partners not accurate: got partner coordinate idx " +
              std::to_string(this->bondPartnerCoordinatesB[i * 3 + dir]) +
              " but expected " +
              std::to_string(this->bondPartnersB[i] * 3 + dir) + " at index " +
              std::to_string(i) + " and dir " + std::to_string(dir) + ".");
        }
        RUNTIME_EXP_IFN(pylimer_tools::utils::contains(
                          this->bondsOfIndex[this->bondPartnersA[i]], i),
                        "Reverse-link is incorrect.");
        RUNTIME_EXP_IFN(pylimer_tools::utils::contains(
                          this->bondsOfIndex[this->bondPartnersB[i]], i),
                        "Reverse-link is incorrect.");
      }
    }

  }
}
}
