#ifndef DPD_SIMULATOR_H
#define DPD_SIMULATOR_H

#include "../entities/Atom.h"
#include "../entities/EigenNeighbourList.h"
#include "../entities/Universe.h"
#include "../sim/OutputSupportingSimulation.h"
#include "../utils/ExtraEigenTypes.h"
#include "../utils/PerformanceTimer.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#ifdef CEREALIZABLE
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools::sim {

namespace dpd {

  enum DPDPerformanceSections
  {
    TIME_STEPPING,
    FORCES,
    PAIR_FORCE,
    BOND_FORCE,
    OUTPUT,
    SHIFT,
    RELOCATION,
    MODIFY,
    NUM_PERFORMANCE_SECTIONS
  };

  static const std::array<std::string, 8> DPDPerformanceSectionNames = {
    "Time-Stepping", "Forces", "Pair-Forces", "Bond-Forces",
    "Output",        "Shift",  "Relocation",  "Modify"
  };

  class DPDSimulator : public pylimer_tools::sim::OutputSupportingSimulation
  {

  private:
    DPDSimulator() {}; // not exposed to users, only used by Cereal

#ifdef CEREALIZABLE
    friend class cereal::access;
#endif

    ////////////////////////////////////////////////////////////////
    // configuration
    bool allowRelocationInNetwork = false;
    bool assumeBoxLargeEnough = false;
    bool doDeformation = false;
    bool is2D = false;
    bool shiftOneAtATime = false;
    bool shiftPossibilityEmpty = true;
    double A = 25.;
    double dt = 0.06;
    double gamma = 0.5 * 3. * 3.;
    double highCutoff = 2.0;
    double k = 2.;
    double lambda = 0.65;
    double lowCutoff = 0.5;
    double maxBondLen = 5.;
    double sigma = 3.;
    int crossLinkerType = 2;
    int slipspringBondType = 9;
    long int nStepsDPD = 500;
    long int nStepsMC = 500;

    ////////////////////////////////////////////////////////////////
    // simulation state
    long int currentStep = 0;
    double currentTime = 0.;
    long int numShifts = 0;
    long int numRelocations = 0;
    Eigen::VectorXd currentVelocitiesPlus;
    Eigen::VectorXd currentVelocities;
    Eigen::VectorXd currentForces;
    Eigen::Matrix3d currentStressTensor = Eigen::Matrix3d::Zero();
    double currentPressure = 0.;

    ////////////////////////////////////////////////////////////////
    // state of bond formation
    int bondsToForm = 0;
    std::unordered_map<int, int> maxBondsPerType;
    double bondFormationDistance = 1.1;
    int formBondsEvery = 50;
    int atomTypeBondFormationFrom = 1;
    int atomTypeBondFormationTo = 2;

    ////////////////////////////////////////////////////////////////
    // randomness
    std::mt19937 e2;
    std::uniform_real_distribution<double> uniform_rand_mean0std1;
    std::uniform_real_distribution<double> uniform_rand_between_0_1;

    ////////////////////////////////////////////////////////////////
    // universe structure
    int numAtoms = 0;
    int numBonds = 0;
    int numSlipSprings = 0;
    pylimer_tools::entities::Box box;
    pylimer_tools::entities::Box deformationTargetBox;
    pylimer_tools::entities::Universe universe;

    // atoms
    Eigen::VectorXd coordinates;
    Eigen::ArrayXi idxFunctionalities;
    Eigen::ArrayXb isRelocationTarget;
    std::vector<int> atomTypes;
    std::vector<long int> atomIds;
    std::vector<size_t> chainEndIndices;
    // bonds
    Eigen::ArrayXi bondPartnerCoordinatesA;
    Eigen::ArrayXi bondPartnersA;
    Eigen::ArrayXi bondPartnerCoordinatesB;
    Eigen::ArrayXi bondPartnersB;
    Eigen::ArrayXi bondTypes;
    Eigen::ArrayXd bondDuplicationPenalty;
    Eigen::VectorXd bondBoxOffsets;
    // mapping between atoms and bonds
    std::vector<std::vector<size_t>> bondsOfIndex;

    pylimer_tools::entities::EigenNeighbourList neighbourlist;

  public:
    /**
     * @brief Construct a DPD Simulator from a given Universe
     *
     * @param u The Universe containing the initial system configuration
     * @param crossLinkerType The atom type representing crosslinkers (default:
     * 2)
     * @param slipspringBondType The bond type for slip-springs (default: 9)
     * @param is2D Whether to run a 2D simulation (default: false, currently not
     * supported)
     * @param seed Random seed for reproducibility (default: empty string for
     * random seed)
     */
    DPDSimulator(const pylimer_tools::entities::Universe& u,
                 const int crossLinkerType = 2,
                 const int slipspringBondType = 9,
                 const bool is2D = false,
                 const std::string& seed = "");

    /**
     * @brief Run the main simulation loop with DPD and optionally Monte Carlo
     * steps
     *
     * @param nSteps Number of simulation steps to run
     * @param withMC Whether to include Monte Carlo moves during simulation
     * @param shouldInterrupt Function to check if simulation should be
     * interrupted
     * @param cleanupInterrupt Function to call when cleaning up after
     * interruption
     */
    void runSimulation(const long int nSteps,
                       bool withMC,
                       const std::function<bool()>& shouldInterrupt,
                       const std::function<void()>& cleanupInterrupt);

    /**
     * @brief Run the simulation with default interrupt handlers
     *
     * @param nSteps Number of simulation steps to run
     * @param withMC Whether to include Monte Carlo moves during simulation
     * (default: false)
     */
    void runSimulation(const long int nSteps, const bool withMC = false)
    {
      runSimulation(nSteps, withMC, []() { return false; }, []() {});
    }

    /**
     * @brief Re-calculate stress tensor and pressure from current system state
     *
     * This method updates the current stress tensor and pressure based on the
     * current coordinates and velocities, typically called after changing
     * force field parameters.
     */
    void refreshCurrentState();

    /**
     * @brief Set a new seed for the random number generator
     *
     * @param seed Random seed string (empty string for random initialization)
     */
    void reseedRandomness(const std::string& seed);

    /**
     * @brief Create a new bond between two nodes
     *
     * CAUTION: This is an expensive operation involving resizing of Eigen
     * containers. Use sparingly.
     *
     * @param fromIdx Index of the first node
     * @param toIdx Index of the second node
     * @param bondType Type identifier for the bond
     */
    void addBond(const long int fromIdx,
                 const long int toIdx,
                 const int bondType);

    /**
     * @brief Compute all forces and stress tensor for the current system
     * configuration
     *
     * @param forces Output vector to store computed forces (modified in-place)
     * @param stressTensor Output matrix to store computed stress tensor
     * (modified in-place)
     * @param coordinates Current particle coordinates
     * @param velocities Current particle velocities
     * @param timer Performance timer for profiling force computation
     * @param dt Time step for force computation (default: 0.06)
     * @param cutoff Cutoff distance for force interactions (default: 1.0)
     * @return The computed pressure scalar value
     */
    double computeForces(
      Eigen::VectorXd& forces,
      Eigen::Matrix3d& stressTensor,
      const Eigen::VectorXd& coordinates,
      const Eigen::VectorXd& velocities,
      pylimer_tools::utils::PerformanceTimer<
        DPDPerformanceSections::NUM_PERFORMANCE_SECTIONS>& timer,
      const double dt = 0.06,
      const double cutoff =
        1.0); // unfortunately not const because of the random nr generator

    /**
     * @brief Compute the temperature
     *
     * @param velocities
     * @return double
     */
    static double computeTemperature(const Eigen::VectorXd& velocities);

    inline Eigen::Vector3d computeBondDistance(const int bondIdx) const
    {
      Eigen::Vector3d bondDistances =
        this->coordinates(
          this->bondPartnerCoordinatesB.segment(3 * bondIdx, 3)) -
        this->coordinates(
          this->bondPartnerCoordinatesA.segment(3 * bondIdx, 3)) +
        this->bondBoxOffsets.segment(3 * bondIdx, 3);
      if (this->assumeBoxLargeEnough) {
        this->box.handlePBC(bondDistances);
      }
      return bondDistances;
    }

    /**
     * @brief Compute the length of one specific bond
     *
     * @param bondIdx Index of the bond to measure
     * @return The length of the specified bond
     */
    double computeBondLength(const int bondIdx) const
    {
      return this->computeBondDistance(bondIdx).norm();
    }

    ////////////////////////////////////////////////////////////////
    // MC Procedures
    /**
     * @brief Randomly add new slip-springs between neighboring nodes
     *
     * @param num Number of slip-springs to attempt to create
     * @param bondType Type identifier for the slip-spring bonds (default: 0)
     * @return Number of slip-springs actually created
     */
    int createSlipSprings(const int num, int bondType = 0);

    /**
     * @brief Perform Monte Carlo moves to shift slip-springs along polymer
     * chains
     *
     * @param kbT Thermal energy scale for accepting/rejecting moves
     * (default: 1.0)
     * @return Number of successful shift moves performed
     */
    int shiftSlipSprings(const double kbT = 1.);

    /**
     * @brief Perform Monte Carlo moves to relocate slip-springs to new
     * positions
     *
     * @param kbT Thermal energy scale for accepting/rejecting moves
     * (default: 1.0)
     * @return Number of successful relocation moves performed
     */
    int relocateSlipSprings(const double kbT = 1.);

    ////////////////////////////////////////////////////////////////
    // configuration

    /**
     * @brief Configure whether to assume the simulation box is large enough
     *
     * If bonds could get larger than half the box length, this must be kept
     * false (default). If true, periodic boundary condition handling for bonds
     * is disabled for performance.
     */
    void configAssumeBoxLargeEnough()
    {
      this->assumeBoxLargeEnough = true;
      // this->bondBoxOffsets.setZero();
    }

    /**
     * @brief Check if the simulator assumes the box is large enough
     *
     * @return true if PBC handling for bonds is disabled, false otherwise
     */
    bool assumesBoxLargeEnough() const { return this->assumeBoxLargeEnough; }

    /**
     * @brief Configure the simulation time step
     *
     * @param dt Time step size (default: 0.06)
     */
    void configTimeStep(const double dt = 0.06) { this->dt = dt; }

    /**
     * @brief Configure the lambda parameter for modified velocity Verlet
     * integration
     *
     * @param l Lambda parameter value (default: 0.65)
     */
    void configLambda(const double l) { this->lambda = l; }

    /**
     * @brief Configure the spring constant for bond interactions
     *
     * @param nk Spring constant value (default: 2.0)
     */
    void configSpringConstant(const double nk) { this->k = nk; }

    /**
     * @brief Get the current spring constant value
     *
     * @return Current spring constant
     */
    double getSpringConstant() const { return this->k; }

    /**
     * @brief Configure the lower cutoff distance for slip-spring creation
     *
     * @param lowC Lower cutoff distance (must be less than high cutoff)
     */
    void configSlipspringLowCutoff(const double lowC)
    {
      INVALIDARG_EXP_IFN(lowC < this->highCutoff,
                         "The low cutoff must be lower than the high cutoff.");
      this->lowCutoff = lowC;
    }

    /**
     * @brief Configure the higher cutoff distance for slip-spring creation
     *
     * @param highC Higher cutoff distance (must be greater than low cutoff)
     */
    void configSlipspringHighCutoff(const double highC)
    {
      INVALIDARG_EXP_IFN(this->lowCutoff < highC,
                         "The low cutoff must be lower than the high cutoff.");
      this->highCutoff = highC;
    }

    /**
     * @brief Configure the sigma parameter for DPD interactions
     *
     * This also automatically updates gamma = 0.5 * sigma^2
     *
     * @param newSigma Sigma parameter value for noise strength
     */
    void configSigma(const double newSigma)
    {
      this->sigma = newSigma;
      this->gamma = 0.5 * newSigma * newSigma;
    }

    /**
     * @brief Configure the A parameter for DPD repulsive interactions
     *
     * @param newA Repulsion strength parameter
     */
    void configA(const double newA) { this->A = newA; }

    void configAllowRelocationInNetwork(const bool allowReloc)
    {
      this->allowRelocationInNetwork = allowReloc;
      this->isRelocationTarget.setConstant(false);
      this->chainEndIndices.clear();
      for (size_t i = 0; i < this->numAtoms; ++i) {
        if (this->idxFunctionalities[i] < 2) {
          this->isRelocationTarget[i] = true;
          this->chainEndIndices.push_back(i);
        } else if (allowReloc && this->atomTypes[i] != this->crossLinkerType) {
          // check if this is a bead that's connected to a crosslinker
          if (this->idxFunctionalities[i] == 2) {
            // this "any of" here might be a bit overkill, given that we know
            // that there are only two cases
            bool oneOfTheBeadsIsXlink = std::ranges::any_of(
              this->bondsOfIndex[i], [&](const size_t bondIdx) {
                assert(this->bondPartnersA[bondIdx] == i ||
                       this->bondPartnersB[bondIdx] == i);
                return this->atomTypes[this->bondPartnersA[bondIdx]] ==
                         this->crossLinkerType ||
                       this->atomTypes[this->bondPartnersB[bondIdx]] ==
                         this->crossLinkerType;
              });
            this->isRelocationTarget[i] = true;
            this->chainEndIndices.push_back(i);
          }
        }
      }
    }

    void startMeasuringMSDForAtoms(const std::vector<size_t>& atomIds);

    void configNumStepsMC(const long int steps = 500)
    {
      this->nStepsMC = steps;
    }

    int getSlipSpringBondType() const { return this->slipspringBondType; }

    long int getNumStepsMC() const { return this->nStepsMC; }

    void configNumStepsDPD(const long int steps = 500)
    {
      this->nStepsDPD = steps;
    }

    long int getNumStepsDPD() const { return this->nStepsDPD; }

    void configShiftPossibilityEmpty(
      const bool shiftPossibilityEmptyConfig = true)
    {
      this->shiftPossibilityEmpty = shiftPossibilityEmptyConfig;
    }
    bool getShiftPossibilityEmpty() const
    {
      return this->shiftPossibilityEmpty;
    }

    void configShiftOneAtATime(const bool shiftOne = false)
    {
      this->shiftOneAtATime = shiftOne;
    }

    bool getShiftOneAtATime() const { return this->shiftOneAtATime; }

    void configBoxDeformation(const pylimer_tools::entities::Box& newBox)
    {
      this->deformationTargetBox = newBox;
      this->doDeformation = true;
    }

    void deformBoxImmediately(const pylimer_tools::entities::Box& newBox);

    void configBondFormation(
      const int numBondsToForm,
      const std::unordered_map<int, int>& numBondsPerType,
      const double bondFormationDist = 1.1,
      const int formBondEvery = 10,
      const int formFrom = 2,
      const int formTo = 1)
    {
      INVALIDARG_EXP_IFN(bondFormationDist > 0.0,
                         "Bond formation distance must be > 0, got " +
                           std::to_string(bondFormationDist) + ".");
      INVALIDARG_EXP_IFN(formBondEvery > 0,
                         "Bond formation iteration must be > 0, got " +
                           std::to_string(formBondEvery) + ".");
      const std::map<int, int> functionalities =
        this->universe.determineFunctionalityPerType();
      for (const auto& type : functionalities) {
        INVALIDARG_EXP_IFN(
          numBondsPerType.at(type.first) >= type.second,
          "Number of bonds per type must be bigger than their current "
          "functionality. Got issue with type " +
            std::to_string(type.first) +
            ": it currently has a functionality of " +
            std::to_string(type.second) + ", but " +
            std::to_string(numBondsPerType.at(type.first)) + " was requested.");
      }
      this->bondsToForm = numBondsToForm;
      this->maxBondsPerType = numBondsPerType;
      this->bondFormationDistance = bondFormationDist;
      this->formBondsEvery = formBondEvery;
      this->atomTypeBondFormationFrom = formFrom;
      this->atomTypeBondFormationTo = formTo;
    }

    int getNrOfBondsToForm() const { return this->bondsToForm; }

    /**
     * @brief Try to form as many bonds as we can
     *
     */
    void attemptBondFormation();

    /**
     * @brief Sets the bond duplication penalty back to defaults.
     *
     * The bond duplication penalty is here, to enforce
     * "slip springs connecting two already bonded beads do not contribute in
     * the DPD steps"
     */
    void resetBondDuplicationPenalty();

    void resetBondDuplicationPenalty(size_t atomIdx);

    ////////////////////////////////////////////////////////////////
    // results access & export

    pylimer_tools::entities::Universe getUniverse(
      bool withSlipsprings = true) const;

    double getTimestep() override { return this->dt; }
    double getCurrentTime(double currentStep) override
    {
      return this->currentTime;
    }

    int getCurrentTimestep() const { return this->currentStep; }

    /**
     * @brief Get access to the current stress tensor
     *
     * @return Current 3x3 stress tensor matrix
     */
    Eigen::Matrix3d getStressTensor() override
    {
      return this->currentStressTensor;
    }

    /**
     * @brief Get the number of successful slip-spring shift moves performed
     *
     * @return Number of shifts performed
     */
    int getNumShifts() override { return this->numShifts; }

    /**
     * @brief Get the number of successful slip-spring relocation moves
     * performed
     *
     * @return Number of relocations performed
     */
    int getNumRelocations() override { return this->numRelocations; }

    /**
     * @brief Get the total number of particles in the system
     *
     * @return Number of particles (same as getNumAtoms)
     */
    size_t getNumParticles() override { return this->numAtoms; }

    /**
     * @brief Get the current number of slip-springs
     *
     * @return Number of slip-springs in the system
     */
    size_t getNumSlipSprings() const { return this->numSlipSprings; }

    /**
     * @brief Get the number of regular bonds (excluding slip-springs)
     *
     * @return Number of regular bonds
     */
    size_t getNumBonds() override { return this->numBonds; }

    /**
     * @brief Get the number of extra bonds (slip-springs)
     *
     * @return Number of slip-springs
     */
    size_t getNumExtraBonds() override { return this->numSlipSprings; }

    /**
     * @brief Get the number of bonds to be formed during simulation
     *
     * @return Number of bonds to form
     */
    long int getNumBondsToForm() override { return this->bondsToForm; }

    /**
     * @brief Get the total number of atoms in the system
     *
     * @return Number of atoms
     */
    size_t getNumAtoms() override { return this->numAtoms; }

    /**
     * @brief Get the number of extra atoms (always 0 for DPD simulations)
     *
     * @return 0 (no extra atoms in DPD)
     */
    size_t getNumExtraAtoms() override { return 0; }

    /**
     * @brief Get the current simulation box volume
     *
     * @return Simulation box volume
     */
    double getVolume() override { return this->box.getVolume(); }

    /**
     * @brief Get gamma parameter (not used in DPD, returns -1)
     *
     * @return -1 (not applicable for DPD)
     */
    double getGamma() override { return -1.; }

    /**
     * @brief Get residual (not used in DPD, returns -1)
     *
     * @return -1 (not applicable for DPD)
     */
    double getResidual() override { return -1.; }

    /**
     * @brief Get the lengths of all bonds in the system
     *
     * @return Vector containing the length of each bond
     */
    Eigen::VectorXd getBondLengths() override;

    /**
     * @brief Get the current particle coordinates
     *
     * @return Vector of particle coordinates (x1,y1,z1,x2,y2,z2,...)
     */
    Eigen::VectorXd getCoordinates() override;

    /**
     * @brief Get the current particle velocities
     *
     * @return Vector of particle velocities (vx1,vy1,vz1,vx2,vy2,vz2,...)
     */
    Eigen::VectorXd getVelocities() const;

    /**
     * @brief Get the current system temperature
     *
     * @return Current temperature computed from kinetic energy
     */
    double getTemperature() override;

    ////////////////////////////////////////////////////////////////
    // validation

    /**
     * @brief Validate the current simulation state for debugging purposes
     *
     * Performs comprehensive checks on internal data structure consistency
     * and throws exceptions if issues are found.
     */
    void validateState();

    /**
     * @brief Static validation method for debugging (currently empty)
     */
    static void validateDebugState();

    /**
     * @brief Validate the neighbor list consistency for debugging purposes
     *
     * @param cutoff Cutoff distance to check against neighbor list
     */
    void validateNeighbourlist(double cutoff);

    /**
     * @brief Generate a random number from uniform distribution with mean 0,
     * std 1
     *
     * @return Random number from uniform distribution
     */
    double getUniformRandMean0Std1()
    {
      return this->uniform_rand_mean0std1(this->e2);
    }

    /**
     * @brief Generate a random number from uniform distribution between 0 and 1
     *
     * @return Random number between 0 and 1
     */
    double getUniformRandBetween0And1()
    {
      return this->uniform_rand_between_0_1(this->e2);
    }

    /**
     * @brief Get the number of OpenMP threads being used
     *
     * @return Number of OpenMP threads (1 if OpenMP not available)
     */
    int getNumOmpThreads() const
    {
#ifdef OPENMP_FOUND
      return omp_get_max_threads();
#endif
      return 1; // no OpenMP found, so only one thread
    }

    /**
     * @brief Set the number of OpenMP threads to use
     *
     * @param numThreads Number of threads (must be 1 if OpenMP not available)
     */
    void setNumOmpThreads(int numThreads)
    {
#ifdef OPENMP_FOUND
      omp_set_num_threads(numThreads);
#else
      RUNTIME_EXP_IFN(numThreads == 1,
                      "Cannot set number of OpenMP threads, OpenMP not found.");
#endif
    }

    ////////////////////////////////////////////////////////////////
    // serialization

#ifdef CEREALIZABLE
    template<class Archive>
    void serialize(Archive& ar, std::uint32_t const version)
    {
      ar(cereal::virtual_base_class<OutputSupportingSimulation>(this));
      // configuration
      ar(maxBondLen,
         is2D,
         shiftPossibilityEmpty,
         shiftOneAtATime,
         lambda,
         k,
         lowCutoff,
         highCutoff,
         A,
         sigma,
         gamma,
         nStepsDPD,
         nStepsMC,
         doDeformation,
         dt);
      if (version > 1) {
        ar(allowRelocationInNetwork, crossLinkerType);
      }
      if (version > 2) {
        ar(slipspringBondType);
      }
      if (version > 4) {
        ar(assumeBoxLargeEnough);
      }
      // state of bond formation
      if (version > 1) {
        ar(bondsToForm,
           maxBondsPerType,
           bondFormationDistance,
           formBondsEvery,
           atomTypeBondFormationFrom,
           atomTypeBondFormationTo);
      }
      // simulation state
      ar(currentStep,
         currentTime,
         numShifts,
         numRelocations,
         currentVelocitiesPlus,
         currentVelocities,
         currentForces,
         currentStressTensor);
      if (version > 5) {
        ar(bondDuplicationPenalty);
      } else {
        this->resetBondDuplicationPenalty();
      }
      // randomness
      ar(e2, uniform_rand_mean0std1, uniform_rand_between_0_1);
      // universe structure
      ar(numAtoms,
         numBonds,
         numSlipSprings,
         box,
         deformationTargetBox,
         universe,
         // -> atoms
         coordinates,
         idxFunctionalities,
         atomTypes,
         atomIds,
         chainEndIndices,
         // -> bonds
         bondPartnerCoordinatesA,
         bondPartnerCoordinatesB,
         bondPartnersA,
         bondPartnersB,
         bondTypes,
         bondsOfIndex);
      if (version > 1) {
        ar(isRelocationTarget);
      }
      if (version > 3) {
        ar(bondBoxOffsets);
      } else {
        this->resetBondOffsets();
      }
      // neighbourlist
      ar(neighbourlist);
    }

    static DPDSimulator readRestartFile(std::string filename)
    {
      DPDSimulator res;
      pylimer_tools::utils::deserializeFromFile<DPDSimulator>(res, filename);
      return res;
    };

    void writeRestartFile(std::string& filename) override
    {
      pylimer_tools::utils::serializeToFile<DPDSimulator>(*this, filename);
    };
#endif

  protected:
    void addSlipSprings(std::vector<size_t>& partnerA,
                        std::vector<size_t>& partnerB,
                        const int bondType = 9);

    bool attemptSlipSpringShift(const size_t springIdx,
                                const size_t endToShift,
                                const double kbT = 1.);

    bool attemptSlipSpringShift(const size_t springIdx, const double kbT = 1.);

    void replaceSlipSpringPartner(const size_t springIdx,
                                  const size_t partnerBefore,
                                  const size_t partnerAfter);

    /**
     * @brief Reset the offset required for the PBC bonds,
     *
     * This enables us having bonds longer than the box.
     * However, when calling this function, you reset that fact, i.e., only
     * call this if you are sure that you currently do not have bonds that
     * escape the box.
     *
     */
    void resetBondOffsets();

    void resetBondOffset(int bondIdx);
  };
};
}

#ifdef CEREALIZABLE
CEREAL_REGISTER_TYPE(pylimer_tools::sim::dpd::DPDSimulator);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
  pylimer_tools::sim::OutputSupportingSimulation,
  pylimer_tools::sim::dpd::DPDSimulator);
CEREAL_CLASS_VERSION(pylimer_tools::sim::dpd::DPDSimulator, 6);
#endif

#endif
