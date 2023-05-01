#include "./DPDSimulator.h"

namespace pylimer_tools {
namespace calc {
  namespace dpd {

    DPDSimulator::DPDSimulator(const pylimer_tools::entities::Universe u,
                               const int crosslinkerType,
                               const bool is2D,
                               const std::string seed)
      : box(u.getBox())
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
      double mean = 0.;
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
        Eigen::Map<ArrayXli, Eigen::Unaligned>(edges["edge_from"].data(),
                                               edges["edge_from"].size())
          .cast<int>();
      this->bondPartnersB = Eigen::Map<ArrayXli, Eigen::Unaligned>(
                              edges["edge_to"].data(), edges["edge_to"].size())
                              .cast<int>();
      this->bondTypes = Eigen::Map<ArrayXli, Eigen::Unaligned>(
                          edges["edge_type"].data(), edges["edge_type"].size())
                          .cast<int>();

      this->neighbourlist = pylimer_tools::entities::EigenNeighbourList(
        coordinates, this->box, 1.0);
      this->numAtoms = this->coordinates.size() / 3;
      this->numBonds = this->bondPartnersA.size();
      this->idxFunctionalities = Eigen::ArrayXi::Zero(this->numAtoms);
      this->atomTypes = u.getPropertyValues<int>("type");
      this->atomIds = u.getPropertyValues<long int>("id");

      for (size_t i = 0; i < this->numBonds; ++i) {
        std::vector<size_t> bonds;
        bonds.reserve(4);
        this->bondsOfIndex.push_back(bonds);
      }
      for (size_t i = 0; i < this->numBonds; ++i) {
        this->bondsOfIndex[this->bondPartnersA[i]].push_back(i);
        this->bondsOfIndex[this->bondPartnersB[i]].push_back(i);
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
    void DPDSimulator::runSimulation(const long int nSteps,
                                     double dt,
                                     double lambda,
                                     bool withMC)
    {
      Eigen::VectorXd velocitiesPlus = this->currentVelocitiesPlus;
      Eigen::VectorXd velocities = this->currentVelocities;
      Eigen::VectorXd forces = this->currentForces;
      Eigen::Matrix3d stressTensor = Eigen::Matrix3d::Zero();

      std::ios::sync_with_stdio(false);
      std::string outputBuffer = "";
      outputBuffer.reserve(80 * 20);
      std::cout << "Step\tTemperature\tPressure\t#Shift\t#Relocation\tStress[1,"
                   "1]\tStress[2,2]\tStress[3,3]\tStress[1,2]\tStress[1,3]"
                   "\tStress[2,3]\n";

      int numShifts = 0;
      int numRelocations = 0;

      const double halfDt = 0.5 * dt;
      double temperature = this->computeTemperature(velocities);
      for (long int step = 0; step < nSteps; step++) {
        if (withMC) {
          numShifts = this->shiftSlipSprings(1. * temperature);
          numRelocations = this->relocateSlipSprings(1. * temperature);
        }
        // update coordinates & velocities
        velocitiesPlus = velocities + halfDt * forces;
        this->coordinates += dt * velocitiesPlus;
        velocities += lambda * dt * forces;

        // re-compute the forces with these updated coordinates & velocities
        double pressure =
          computeForces(forces, stressTensor, coordinates, velocities, dt, 1.0);

        // correct the velocities
        velocities = velocitiesPlus + halfDt * forces;
        temperature = this->computeTemperature(velocities);

        // output
        outputBuffer = "";
        outputBuffer += std::to_string(step + this->currentStep) + "\t";
        outputBuffer += std::to_string(temperature) + "\t";
        outputBuffer += std::to_string(pressure) + "\t";
        outputBuffer += std::to_string(numShifts) + "\t";
        outputBuffer += std::to_string(numRelocations) + "\t";
        outputBuffer += std::to_string(stressTensor(0, 0)) + "\t";
        outputBuffer += std::to_string(stressTensor(1, 1)) + "\t";
        outputBuffer += std::to_string(stressTensor(2, 2)) + "\t";
        outputBuffer += std::to_string(stressTensor(0, 1)) + "\t";
        outputBuffer += std::to_string(stressTensor(0, 2)) + "\t";
        outputBuffer += std::to_string(stressTensor(1, 2)) + "\t";
        outputBuffer += "\n";
        std::cout << outputBuffer;
        if (step % 25 == 0) {
          std::flush(std::cout);
        }
      }

      // finish up
      std::ios::sync_with_stdio(true);
      this->currentForces = forces;
      this->currentVelocities = velocities;
      this->currentVelocitiesPlus = velocitiesPlus;
      this->currentStep += nSteps;
      this->currentStressTensor = stressTensor;
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
    double DPDSimulator::computeForces(Eigen::VectorXd& forces,
                                       Eigen::Matrix3d& stressTensor,
                                       const Eigen::VectorXd& coordinates,
                                       const Eigen::VectorXd& velocities,
                                       const double dt,
                                       const double cutoff)
    {
      // initialisation
      assert(coordinates.size() == velocities.size());
      assert(forces.size() == coordinates.size());
      double pressure = 0.0;
      forces = Eigen::VectorXd::Zero(coordinates.size());
      stressTensor = Eigen::Matrix3d::Zero();

      // actual computation
      // (attractive) bond forces
      Eigen::VectorXd bondDistances =
        coordinates(this->bondPartnersA) - coordinates(this->bondPartnersB);
      this->box.handlePBC(bondDistances);
      forces(this->bondPartnersA) -= this->k * bondDistances;
      forces(this->bondPartnersB) += this->k * bondDistances;

      Eigen::Vector3d pairdistance;
      Eigen::Vector3d pairdistanceNormed;
      Eigen::Vector3d velocitydiff;
      Eigen::Vector3d pairForce;

      for (size_t i = 0; i < coordinates.size() / 3; ++i) {
        Eigen::ArrayXi neighbors =
          this->neighbourlist.getIndicesCloseToCoordinates(
            coordinates.segment(3 * i, 3));

        // pair forces
        for (size_t neigh_idx = 0; neigh_idx < neighbors.size(); ++neigh_idx) {
          const size_t j = neighbors[neigh_idx];
          if (j <= i) {
            continue;
          }
          pairdistance =
            coordinates.segment(3 * i, 3) - coordinates.segment(3 * j, 3);
          const double rNorm = pairdistance.norm();
          if (rNorm >= cutoff || rNorm < 1e-12) {
            // this is a performance bottleneck, commonly solved by using
            // neighbour lists
            continue;
          }

          const double one_minus_rnorm = 1. - rNorm;
          const double one_minus_rnorm2 = (1. - rNorm) * (1. - rNorm);
          pairdistanceNormed = pairdistance / rNorm;

          // conservative repulsion force
          pairForce = this->A * one_minus_rnorm * pairdistanceNormed;

          // dissipative/drag force
          velocitydiff =
            velocities.segment(3 * i, 3) - velocities.segment(3 * j, 3);
          const double rij_dot_vij = pairdistanceNormed.dot(velocitydiff);
          const double gamma_weighted_rij_dot_vij =
            this->gamma * one_minus_rnorm2 * rij_dot_vij;

          pairForce += -gamma_weighted_rij_dot_vij * pairdistanceNormed;

          // random force
          const double constant_rnd_prefix =
            this->sigma * one_minus_rnorm / sqrt(dt);
          const double random_val = this->uniform_rand_mean0std1(this->e2);

          pairForce += constant_rnd_prefix * random_val * pairdistanceNormed;

          // actually assign the new forces
          forces.segment(3 * i, 3) += pairForce;
          forces.segment(3 * j, 3) -= pairForce;

          // pressure update
          pressure += forces.segment(3 * i, 3).dot(pairdistance);
          stressTensor += forces.segment(3 * i, 3) * pairdistance.transpose();
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
    };

    /**
     * @brief Randomly add new slip-springs
     *
     * @param num
     * @param loCutoff
     * @param hiCutoff
     * @return int
     */
    int DPDSimulator::createSlipSprings(const int num)
    {
      int createdLastIteration = 100;
      int totalCreated = 0;
      std::vector<size_t> candidates;

      std::vector<size_t> slipSpringFrom;
      std::vector<size_t> slipSpringTo;

      // randomly permute the atoms to start the search with
      std::vector<size_t> sourceIds;
      sourceIds.reserve(this->numAtoms);
      for (size_t i = 0; i < this->numAtoms; ++i) {
        sourceIds.push_back(i);
      }
      std::shuffle(sourceIds.begin(), sourceIds.end(), this->e2);

      // search for neighbours that are elibile
      while (createdLastIteration > 0 && totalCreated < num) {
        for (size_t i : sourceIds) {
          int numCandidates = 0;
          // for each atom, search for possible partners
          Eigen::ArrayXi pairs =
            this->neighbourlist.getIndicesCloseToCoordinates(
              this->coordinates.segment(3 * i, 3), this->highCutoff);
          for (size_t j = 0; j < pairs.size(); ++j) {
            Eigen::Vector3d distance =
              this->coordinates.segment(3 * i, 3) -
              this->coordinates.segment(3 * pairs[j], 3);
            if (distance.norm() > this->lowCutoff &&
                distance.norm() <= this->highCutoff) {
              candidates[numCandidates++] = j;
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
          if (totalCreated >= num) {
            break;
          }
        }
      }
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
      size_t sizeBefore = this->numBonds + this->numSlipSprings;
      this->bondPartnersA.conservativeResize(sizeBefore + partnerA.size());
      this->bondPartnersB.conservativeResize(sizeBefore + partnerB.size());
      this->bondTypes.conservativeResize(sizeBefore + partnerB.size());

      this->bondPartnersA.segment(sizeBefore, partnerA.size()) =
        Eigen::Map<ArrayXst, Eigen::Unaligned>(partnerA.data(), partnerA.size())
          .cast<int>();
      this->bondPartnersB.segment(sizeBefore, partnerB.size()) =
        Eigen::Map<ArrayXst, Eigen::Unaligned>(partnerB.data(), partnerB.size())
          .cast<int>();
      this->bondTypes.segment(sizeBefore, partnerB.size()) = bondType;
      for (size_t i = 0; i < partnerA.size(); ++i) {
        this->bondsOfIndex[this->bondPartnersA[sizeBefore + i]].push_back(
          sizeBefore + i);
        this->bondsOfIndex[this->bondPartnersB[sizeBefore + i]].push_back(
          sizeBefore + i);
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
      int nAccept;
      std::vector<size_t> candidates;
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
        Eigen::ArrayXi pairs = this->neighbourlist.getIndicesCloseToCoordinates(
          this->coordinates.segment(3 * candidateIndex, 3), this->highCutoff);
        for (size_t j = 0; j < pairs.size(); ++j) {
          Eigen::Vector3d distance =
            this->coordinates.segment(3 * candidateIndex, 3) -
            this->coordinates.segment(3 * pairs[j], 3);
          if (distance.norm() > this->lowCutoff &&
              distance.norm() <= this->highCutoff) {
            candidates[numCandidates++] = j;
          }
        }
        if (numCandidates == 0) {
          continue;
        }
        std::uniform_int_distribution<int> candidateDist(0, numCandidates - 1);
        int candidatePartnerIndex = candidates[candidateDist(this->e2)];

        // compute the Metropolis criterion
        double bondEnergyNow =
          -k * (this->coordinates.segment(partnerA * 3, 3) -
                this->coordinates.segment(partnerB * 3, 3))
                 .squaredNorm();
        double bondEnergyNew =
          -k * (this->coordinates.segment(candidateIndex * 3, 3) -
                this->coordinates.segment(candidatePartnerIndex * 3, 3))
                 .squaredNorm();
        double deltaEnergy = bondEnergyNew - bondEnergyNow;
        bool accept = false;
        if (deltaEnergy < 0.0) {
          accept = true;
        } else {
          double factor = std::exp(-deltaEnergy / kbT);
          if (this->uniform_rand_mean0std1(this->e2) * 0.5 + 1. < factor) {
            accept = true;
          }
        }
        if (accept) {
          this->replaceSlipSpringPartner(springIdx, partnerA, candidateIndex);
          this->replaceSlipSpringPartner(
            springIdx, partnerB, candidatePartnerIndex);
          nAccept++;
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
      for (size_t springIdx = this->numBonds;
           springIdx < (this->numBonds + this->numSlipSprings);
           ++springIdx) {
        const size_t partnerA = this->bondPartnersA[springIdx];
        const size_t partnerB = this->bondPartnersB[springIdx];
        // attempt to shift the spring around partnerA
        n_accept += this->attemptSlipSpringShift(springIdx, partnerA);
        n_accept += this->attemptSlipSpringShift(springIdx, partnerB);
        if (this->bondPartnersA[springIdx] == this->bondPartnersB[springIdx]) {
          // TODO:
        }
      }
      return n_accept;
    }

    /**
     * @brief
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
      const size_t partnerA = this->bondPartnersA[springIdx] == endToShift
                                ? this->bondPartnersA[springIdx]
                                : this->bondPartnersB[springIdx];
      const size_t partnerB = this->bondPartnersB[springIdx] == endToShift
                                ? this->bondPartnersA[springIdx]
                                : this->bondPartnersB[springIdx];
      assert(partnerA == endToShift);
      // attempt to shift the spring around partnerA
      std::uniform_int_distribution<int> dist(
        0, this->idxFunctionalities[partnerA] - 1);
      const size_t selectedBond = this->bondsOfIndex[partnerA][dist(this->e2)];
      const size_t replacementForA =
        this->bondPartnersA[selectedBond] == partnerA
          ? this->bondPartnersB[selectedBond]
          : this->bondPartnersA[selectedBond];
      // compute the Metropolis criterion
      double bondEnergyNow =
        -this->k * (this->coordinates.segment(partnerA * 3, 3) -
                    this->coordinates.segment(partnerB * 3, 3))
                     .squaredNorm();
      double bondEnergyNew =
        -this->k * (this->coordinates.segment(replacementForA * 3, 3) -
                    this->coordinates.segment(partnerB * 3, 3))
                     .squaredNorm();
      double deltaEnergy = bondEnergyNew - bondEnergyNow;
      bool accept = false;
      if (deltaEnergy < 0.0) {
        accept = true;
      } else {
        double factor = std::exp(-deltaEnergy / kbT);
        if (this->uniform_rand_mean0std1(this->e2) * 0.5 + 1. < factor) {
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
      assert(this->bondPartnersA[springIdx] == partnerBefore ||
             this->bondPartnersB[springIdx] == partnerBefore);
      if (this->bondPartnersA[springIdx] == partnerBefore) {
        this->bondPartnersA[springIdx] = partnerAfter;
      } else {
        this->bondPartnersB[springIdx] = partnerAfter;
      }
      this->bondsOfIndex[partnerAfter].push_back(springIdx);
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
      throw std::runtime_error("Invalid internal state.");
    }

    /**
     * @brief Make sure all the structures obey the expected form
     *
     */
    void DPDSimulator::validateState()
    {
      RUNTIME_EXP_IFN(this->coordinates.size() == 3 * this->numAtoms,
                      "State violation");
      RUNTIME_EXP_IFN(this->idxFunctionalities.size() == this->numAtoms,
                      "State violation");
      RUNTIME_EXP_IFN(this->atomTypes.size() == this->numAtoms,
                      "State violation");
      RUNTIME_EXP_IFN(this->atomIds.size() == this->numAtoms,
                      "State violation");
      RUNTIME_EXP_IFN(this->bondsOfIndex.size() == this->numAtoms,
                      "State violation");

      RUNTIME_EXP_IFN(this->bondPartnersA.size() == this->bondPartnersB.size(),
                      "State violation");
      RUNTIME_EXP_IFN(this->bondTypes.size() == this->bondPartnersA.size(),
                      "State violation");
      RUNTIME_EXP_IFN(this->bondPartnersB.size() ==
                        this->numBonds + this->numSlipSprings,
                      "State violation");
    }

  }
}
}
