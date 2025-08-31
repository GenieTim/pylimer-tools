#pragma once

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/EigenNeighbourList.h"
#include "../entities/Universe.h"
#include "../sim/MEHPForceBalance.h"
#include "../sim/MEHPForceBalance2.h"
#include "../sim/MEHPForceRelaxation.h"
#include "../sim/MEHPUtilityStructures.h"
#include "../utils/BoolUtils.h"
#include "RandomWalker.h"
#include "StringUtils.h"
#include "VectorUtils.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stack>
#include <string>
#include <vector>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif
#ifdef CEREALIZABLE
#include "../utils/CerealUtils.h"
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

#include "utilityMacros.h"
#include <random>

namespace pylimer_tools::utils {

#define UNCONNECTED -1
#define EMPTY_BACKGROUND -2

enum BackTrackStatus
{
  STOP,
  TRACK_FORWARD,
  TRACK_BACKWARD,
};

struct CrosslinkerUniverse
{
  std::vector<int> xlinkTypes;
  std::vector<double> xlinkX;
  std::vector<double> xlinkY;
  std::vector<double> xlinkZ;
  std::vector<std::vector<long int>> strandsOfXlink;
  std::vector<int> xlinkChainId;

  std::vector<long int> strandFrom;
  std::vector<long int> strandTo;
  std::vector<int> beadsInStrand;
  std::vector<int> strandBeadType;
  std::vector<double> beadDistanceInStrand;
  std::vector<double> meanSquaredBeadDistanceInStrand;

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(xlinkTypes,
       xlinkX,
       xlinkY,
       xlinkZ,
       strandsOfXlink,
       xlinkChainId,
       strandFrom,
       strandTo,
       beadsInStrand,
       strandBeadType,
       beadDistanceInStrand,
       meanSquaredBeadDistanceInStrand);
  }
#endif
};

class MaxDistanceProvider
{
public:
  virtual ~MaxDistanceProvider() = default;

  virtual double getMaxDistance(double N) = 0;
  [[nodiscard]] virtual std::unique_ptr<MaxDistanceProvider> clone() const = 0;
};

class LinearMaxDistanceProvider : public MaxDistanceProvider
{
private:
  double multiplier;

public:
  explicit LinearMaxDistanceProvider(const double mult = 1.)
    : multiplier(mult)
  {
  }

  double getMaxDistance(const double N) override
  {
    return N * this->multiplier;
  }

  [[nodiscard]] std::unique_ptr<MaxDistanceProvider> clone() const override
  {
    return std::make_unique<LinearMaxDistanceProvider>(this->multiplier);
  }

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& archive)
  {
    archive(multiplier);
  }
#endif
};

class ZScoreMaxDistanceProvider : public MaxDistanceProvider
{
private:
  double innerMultiplier;
  double stdMultiplier;

public:
  explicit ZScoreMaxDistanceProvider(const double stdMult = 3.,
                                     const double innerMult = 1.)
    : innerMultiplier(innerMult)
    , stdMultiplier(stdMult)
  {
  }

  double getMaxDistance(const double N) override
  {
    return this->stdMultiplier * std::sqrt(N * this->innerMultiplier);
  }

  [[nodiscard]] std::unique_ptr<MaxDistanceProvider> clone() const override
  {
    return std::make_unique<ZScoreMaxDistanceProvider>(this->stdMultiplier,
                                                       this->innerMultiplier);
  }

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& archive)
  {
    archive(innerMultiplier, stdMultiplier);
  }
#endif
};

class NoMaxDistanceProvider : public MaxDistanceProvider
{
public:

  double getMaxDistance(const double N) override { return -1.; }

  [[nodiscard]] std::unique_ptr<MaxDistanceProvider> clone() const override
  {
    return std::make_unique<NoMaxDistanceProvider>();
  }

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& archive)
  {
    // No data to serialize – but needed for Cereal
    // archive();
  }
#endif
};
}
#ifdef CEREALIZABLE
CEREAL_REGISTER_TYPE(pylimer_tools::utils::ZScoreMaxDistanceProvider);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
  pylimer_tools::utils::MaxDistanceProvider,
  pylimer_tools::utils::ZScoreMaxDistanceProvider);

CEREAL_REGISTER_TYPE(pylimer_tools::utils::LinearMaxDistanceProvider);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
  pylimer_tools::utils::MaxDistanceProvider,
  pylimer_tools::utils::LinearMaxDistanceProvider);

CEREAL_REGISTER_TYPE(pylimer_tools::utils::NoMaxDistanceProvider);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
  pylimer_tools::utils::MaxDistanceProvider,
  pylimer_tools::utils::NoMaxDistanceProvider);
#endif

namespace pylimer_tools::utils {
class MCUniverseGenerator
{
public:
  MCUniverseGenerator(const double Lx = 10.,
                      const double Ly = 10.,
                      const double Lz = 10.)
  {
    std::random_device rd{};
    this->maxDistanceProvider.reset(new NoMaxDistanceProvider());
    this->rng = std::mt19937(rd());
    this->distX = std::uniform_real_distribution<double>(0.0, Lx);
    this->distY = std::uniform_real_distribution<double>(0.0, Ly);
    this->distZ = std::uniform_real_distribution<double>(0.0, Lz);
    this->box = pylimer_tools::entities::Box(Lx, Ly, Lz);
    this->setBeadDistance(0.965);
  }

  MCUniverseGenerator(const pylimer_tools::entities::Box& box)
    : MCUniverseGenerator(box.getLx(), box.getLy(), box.getLz())
  {
    this->box = box;
  }

  // copy constructor
  MCUniverseGenerator(const MCUniverseGenerator& other)
    : beadDistance(other.beadDistance)
    , meanSquaredBeadDistance(other.meanSquaredBeadDistance)
    , primaryLoopProbability(other.primaryLoopProbability)
    , secondaryLoopProbability(other.secondaryLoopProbability)
    , nMcSteps(other.nMcSteps)
    , rng(other.rng)
    , distX(other.distX)
    , distY(other.distY)
    , distZ(other.distZ)
    , simplifiedUniverse(other.simplifiedUniverse)
    , originalNrOfAvailableCrosslinkSites(
        other.originalNrOfAvailableCrosslinkSites)
    , nrOfAvailableCrosslinkSites(other.nrOfAvailableCrosslinkSites)
    , remainingCrossLinkerFunctionality(other.remainingCrossLinkerFunctionality)
    , xlinkNeighbourList(other.xlinkNeighbourList)
    , box(other.box)
  {
    // clone/copy the maxDistanceProvider
    if (other.maxDistanceProvider) {
      this->maxDistanceProvider = other.maxDistanceProvider->clone();
    } else {
      this->maxDistanceProvider = std::make_unique<NoMaxDistanceProvider>();
    }
  }

  void setSeed(const unsigned int seed) { this->rng.seed(seed); }

  void setBeadDistance(double newBeadDistance, bool updateMeanSquared = true);

  double getConfiguredBeadDistance() const { return this->beadDistance; }

  void setMeanSquaredBeadDistance(double newMeanSquaredBeadDistance,
                                  bool updateMean = true);

  double getConfiguredMeanSquaredBeadDistance() const
  {
    return this->meanSquaredBeadDistance;
  }

  void configNrOfMCSteps(const size_t newNrOfMCSteps)
  {
    this->nMcSteps = newNrOfMCSteps;
  }

  void configPrimaryLoopProbability(double newPrimaryLoopProbability);

  void configSecondaryLoopProbability(double newSecondaryLoopProbability);

  void disableMaxDistance();

  void useLinearMaxDistance(double newMultiplier);

  void useZScoreMaxDistance(double newStdMultiplier,
                            double innerMultiplier = 1.);

  void configMaxDistanceProvider(
    std::unique_ptr<MaxDistanceProvider> newMaxDistanceProvider);

  /**
   * @brief Get the Universe object after actually sampling strands' beads and
   * their positions.
   *
   * @return pylimer_tools::entities::Universe
   */
  pylimer_tools::entities::Universe getUniverse();

  /**
   * @brief Add free atoms with a specified type and functionality for
   * possible crosslinking later
   *
   * @param coordinates
   * @param crosslinkerFunctionality
   * @param crossLinkerAtomType
   */
  void addCrosslinkersAt(const Eigen::VectorXd& coordinates,
                         int crosslinkerFunctionality = 4,
                         int crossLinkerAtomType = 2);

  /**
   * @brief Add free atoms with a specified type and functionality for
   * possible crosslinking later
   *
   * @param nrOfCrosslinkers
   * @param crosslinkerFunctionality
   * @param crossLinkerAtomType
   * @param whiteNoise
   */
  void addCrosslinkers(int nrOfCrosslinkers,
                       int crosslinkerFunctionality = 4,
                       int crossLinkerAtomType = 2,
                       bool whiteNoise = true);

  /**
   * @brief Add strands which contain crosslinks at random positions
   *
   * @param nrOfStrands
   * @param beadsPerStrand
   * @param functionalizationProbability the probability of each bead being a
   * crosslink
   * @param crosslinkerFunctionality the functionality of the crosslinkers,
   * excl. connections to the strand
   * @param crosslinkerAtomType
   * @param strandAtomType
   * @param whiteNoise
   */
  void addRandomlyFunctionalizedStrands(int nrOfStrands,
                                        std::vector<int> beadsPerStrand,
                                        double functionalizationProbability,
                                        int crosslinkerFunctionality = 4,
                                        int crosslinkerAtomType = 2,
                                        int strandAtomType = 2,
                                        bool whiteNoise = true);

  /**
   * @brief Add strands which contain crosslinks at regularly spaced positions
   *
   * @param nrOfStrands
   * @param beadsPerStrand
   * @param spacingBetweenCrosslinks the number of beads between crosslinks
   * @param offsetToFirstCrosslink offset from start of strand to first
   * crosslink (0-based)
   * @param crosslinkerFunctionality the functionality of the crosslinkers,
   * excl. connections to the strand
   * @param crosslinkerAtomType
   * @param strandAtomType
   * @param whiteNoise
   */
  void addRegularlySpacedFunctionalizedStrands(int nrOfStrands,
                                               std::vector<int> beadsPerStrand,
                                               int spacingBetweenCrosslinks,
                                               int offsetToFirstCrosslink = 0,
                                               int crosslinkerFunctionality = 4,
                                               int crosslinkerAtomType = 2,
                                               int strandAtomType = 2,
                                               bool whiteNoise = true);

  /**
   * @brief add strands with crosslinks as ends
   *
   * @param nrOfCrosslinkStrands
   * @param crosslinkerFunctionality incl. the connection to the strand
   * @param crosslinkerAtomType
   * @param beadsPerStrand excl. crosslinks
   * @param strandAtomType
   * @param whiteNoise
   */
  void addCrosslinkStrands(int nrOfCrosslinkStrands,
                           std::vector<int> beadsPerStrand,
                           int crosslinkerFunctionality = 4,
                           int crosslinkerAtomType = 2,
                           int strandAtomType = 2,
                           bool whiteNoise = true);

  /**
   * @brief Randomly distribute free "chains"
   *
   * @param nrOfSolventChains
   * @param chainLength
   * @param solventAtomType
   */
  void addSolventChains(int nrOfSolventChains,
                        int chainLength,
                        int solventAtomType = 3,
                        bool whiteNoise = true);

  /**
   * @brief Add multiple monofunctional strands with specified bead types,
   * link them to crosslinks
   *
   * @param nrOfStrands
   * @param beadsPerChains
   * @param strandAtomType
   */
  void addMonofunctionalStrands(int nrOfStrands,
                                std::vector<int> beadsPerChains,
                                int strandAtomType = 1);

  /**
   * @brief Add strands which only connect on one end to crosslinks
   *
   * @param nrOfStrands
   * @param chainLength
   * @param strandAtomType
   */
  void addMonofunctionalStrands(int nrOfStrands,
                                int chainLength,
                                int strandAtomType = 1);

  /**
   * @brief Register one strand with specified info
   *
   * @param beadsOfChain the nr of beads in the strand
   * @param strandAtomType the atom type of the beads in the strand
   * @param connectionFrom the `strandFrom` property of the strand
   * @param connectionTo the `strandTo` property of the strand
   */
  void addStrand(const int beadsOfChain,
                 const int strandAtomType = 1,
                 const long int connectionFrom = UNCONNECTED,
                 const long int connectionTo = UNCONNECTED);

  void addStrands(const int nrOfStrands,
                  const std::vector<int> beadsPerChains,
                  const int strandAtomType = 1);

  void addStrands(const int nrOfStrands,
                  const int chainLength,
                  const int strandAtomType = 1)
  {
    const std::vector<int> chainLengths =
      pylimer_tools::utils::initializeWithValue<int>(nrOfStrands, chainLength);
    return this->addStrands(nrOfStrands, chainLengths, strandAtomType);
  }

  /**
   * @brief Link one strand end of one strand to a crosslink
   *
   * @param strandIdx
   */
  void linkStrand(const size_t strandIdx, const double cInfinity = 1.);

  void linkStrandTo(size_t strandIdx, size_t crosslinkIdx);

  /**
   * @brief Add strands, link them to the crosslinks, stop when the callback
   * says so
   *
   * @param linkingController
   * @param cInfinity `C_\infty` for `<R_ee^2>_0` from `N` and `<b^2>`
   */
  void linkStrandsCallback(
    std::function<BackTrackStatus(const MCUniverseGenerator&, long int)>
      linkingController,
    double cInfinity = 1.);

  /**
   * @brief Add strands in between the crosslinkers, link them as appropriate
   *
   * @param targetCrossLinkerConversion "p", the target conversion of the
   * crosslinkers
   * @param cInfinity `C_\infty` for `<R_ee^2>_0` from `N` and `<b^2>`
   */
  void linkStrandsToConversion(const double targetCrossLinkerConversion,
                               const double cInfinity = 1.);

  /**
   * @brief Add strands in between the crosslinkers, link them as appropriate
   *
   * @param targetSolubleFraction "w_sol", the target soluble fraction
   * @param cInfinity `C_\infty` for `<R_ee^2>_0` from `N` and `<b^2>`
   */
  void linkStrandsToSolubleFraction(double targetSolubleFraction,
                                    double cInfinity = 1.);

  /**
   * @brief Remove a strand from the current state
   *
   * @param strandIdx
   */
  void removeStrand(size_t strandIdx);

  /**
   * @brief Remove a crosslink from the current state
   *
   * @param crosslinkIdx
   */
  void removeCrosslink(size_t crosslinkIdx);

  /**
   * @brief Remove soluble parts from the current state
   *
   * @param rescale whether to rescale the box and coordinates, keeping the
   * density constant
   */
  void removeSolubleFraction(bool rescale = true);

  /**
   * @brief Convert the current simplified universe into a force relaxation
   * network
   *
   * @return pylimer_tools::sim::mehp::Network
   */
  pylimer_tools::sim::mehp::Network convertToForceRelaxationNetwork() const;

  pylimer_tools::sim::mehp::MEHPForceRelaxation getForceRelaxation() const;

  /**
   * @brief Convert the current simplified universe to a useable structure for
   * the force balance
   *
   * @return pylimer_tools::sim::mehp::ForceBalanceNetwork
   */
  pylimer_tools::sim::mehp::ForceBalanceNetwork convertToForceBalanceNetwork()
    const;

  pylimer_tools::sim::mehp::MEHPForceBalance getForceBalance() const;

  pylimer_tools::sim::mehp::MEHPForceBalance2 getForceBalance2() const;

  /**
   * @brief Run force relaxation to improve the statistics of the crosslinked
   * strands
   *
   */
  void relaxCrosslinks();

  /**
   * @brief Count how many atoms are currently in the system
   *
   * @return size_t
   */
  size_t getCurrentNrOfAtoms() const;

  /**
   * @brief Count how many bonds are currently in the system
   *
   * @return size_t
   */
  size_t getCurrentNrOfBonds() const;

  double getCurrentCrosslinkerConversion() const;

  double getCurrentStrandsConversion() const;

  size_t getCurrentNrOfAvailableCrosslinkSites() const
  {
    return this->nrOfAvailableCrosslinkSites;
  }

  /**
   * @brief Check if the internal state of the universe is consistent
   *
   */
  void validateInternalState() const;

  /**
   * @brief Add star-like crosslinkers with pre-connected strands
   *
   * @param nrOfStars number of star crosslinkers to add
   * @param functionality functionality of each star crosslinker (number of
   * strands)
   * @param beadsPerStrand number of beads in each strand
   * @param crosslinkerAtomType atom type for the crosslinker
   * @param strandAtomType atom type for the strand beads
   * @param whiteNoise whether to use white noise positioning
   */
  void addStarCrosslinkers(int nrOfStars,
                           int functionality,
                           int beadsPerStrand,
                           int crosslinkerAtomType = 2,
                           int strandAtomType = 1,
                           bool whiteNoise = true);

  /**
   * @brief Link two free strand ends together directly
   *
   * @param strandIdx index of the strand to link
   * @param candidateStrands list of candidate strands to link to
   * @param cInfinity statistical parameter for end-to-end distance calculation
   * @return the idx of the other chosen strand, or -1 if no suitable pair was
   * found
   */
  long int linkStrandToStrand(const size_t strandIdx,
                              const std::vector<size_t>& candidateStrands,
                              const double cInfinity = 1.);

  long int linkStrandToStrand(const size_t strandIdx,
                              const double cInfinity = 1.)
  {
    return this->linkStrandToStrand(
      strandIdx, this->findFreeStrandEnds(), cInfinity);
  };

  bool linkStrandToStrand(const double cInfinity = 1.)
  {
    std::vector<size_t> candidates = this->findFreeStrandEnds();

    if (candidates.empty()) {
      return false; // no free strand ends available
    }

    // shuffle candidates to sample a random start point
    std::ranges::shuffle(candidates, this->rng);

    return this->linkStrandToStrand(candidates[0], candidates, cInfinity) >= 0;
  };

  /**
   * @brief Link free strand ends to each other until target conversion is
   * reached
   *
   * @param targetStrandConversion target conversion of free strand ends
   * @param cInfinity statistical parameter for end-to-end distance calculation
   */
  void linkStrandsToStrandsToConversion(double targetStrandConversion,
                                        double cInfinity = 1.);

  /**
   * @brief Find all free strand ends (not connected to any crosslinker)
   *
   * @return std::vector<size_t> vector of strand indices
   */
  std::vector<size_t> findFreeStrandEnds() const;

  /**
   * @brief General implementation for adding functionalized strands with
   * customizable crosslink selection
   *
   * This method allows for flexible crosslink placement by accepting a custom
   * selector function. The selector function is called for each bead position
   * and returns whether that bead should be converted to a crosslink and what
   * functionality it should have.
   *
   * @param nrOfStrands number of strands to add
   * @param beadsPerStrand vector specifying beads per strand
   * @param crosslinkSelector function that determines if a bead should be
   * converted to crosslink. Takes (strandIndex, beadIndex, totalBeads) and
   * returns (shouldConvert, functionality)
   * @param defaultCrosslinkFunctionality default functionality for crosslinks
   * (used for adding crosslinker atoms)
   * @param crosslinkerAtomType atom type for crosslinkers
   * @param strandAtomType atom type for strand beads
   * @param whiteNoise whether to use white noise for positioning
   */
  void addFunctionalizedStrandsImpl(
    int nrOfStrands,
    const std::vector<int>& beadsPerStrand,
    std::function<std::pair<bool, int>(int, int, int)> crosslinkSelector,
    int defaultCrosslinkFunctionality,
    int crosslinkerAtomType,
    int strandAtomType,
    bool whiteNoise);

#ifdef CEREALIZABLE
  // This method lets cereal know which data members to serialize
  template<class Archive>
  void serialize(Archive& archive)
  {
    archive(beadDistance,
            meanSquaredBeadDistance,
            primaryLoopProbability,
            secondaryLoopProbability,
            maxDistanceProvider,
            nMcSteps,
            rng,
            distX,
            distY,
            distZ,
            simplifiedUniverse,
            originalNrOfAvailableCrosslinkSites,
            nrOfAvailableCrosslinkSites,
            remainingCrossLinkerFunctionality,
            xlinkNeighbourList,
            box); // serialize things by passing them to the archive
  }
#endif

private:
  /// settings
  double beadDistance;
  double meanSquaredBeadDistance;
  double primaryLoopProbability = 1.0;
  double secondaryLoopProbability = 1.0;
  std::unique_ptr<MaxDistanceProvider> maxDistanceProvider;
  size_t nMcSteps = 2000;
  /// random distributions
  std::mt19937 rng;
  std::uniform_real_distribution<double> distX;
  std::uniform_real_distribution<double> distY;
  std::uniform_real_distribution<double> distZ;

  /// state
  CrosslinkerUniverse simplifiedUniverse;
  long int originalNrOfAvailableCrosslinkSites = 0;
  long int nrOfAvailableCrosslinkSites = 0;
  std::vector<int> remainingCrossLinkerFunctionality;
  pylimer_tools::entities::EigenNeighbourList xlinkNeighbourList;
  pylimer_tools::entities::Box box;

  /**
   * @brief Get the Current random seed for reproducibility
   *
   * @return std::string
   */
  std::string getCurrentSeed() const
  {
    std::ostringstream oss;
    oss << this->rng;
    return oss.str();
  }

  /**
   * @brief Do a random walk of certain length to add a chain
   *
   * @param chainLen the number of additional atoms to add to the chain
   */
  Eigen::VectorXd sampleFreeChainCoordinates(int chainLen);

  /**
   * @brief Do a random walk of certain length starting somewhere to add a
   * chain
   *
   * @param idxFrom the starting crosslink of the dangling chain
   * @param chainLen the number of additional atoms to add to the chain
   */
  Eigen::VectorXd sampleDanglingChainCoordinates(size_t idxFrom, int chainLen);

  /**
   * @brief Do a random walk of certain length to add a chain from one to
   * another atom
   *
   * @param from the atom to start the random walk from
   * @param to the atom to end the random walk at
   * @param chainLen the number of atoms to add in between from and to
   */
  Eigen::VectorXd sampleStrandCoordinates(size_t from, size_t to, int chainLen);

  /**
   * @brief Tell a strand how it is connected to a crosslink.
   * Automatically decides whether it is `from` or `to`
   *
   * @param strandIdx
   * @param crosslinkIdx
   * @param ignoreInexistent whether to throw an error if the crosslink does
   * not exist yet, or not
   */
  void linkStrandToCrosslink(size_t strandIdx,
                             size_t crosslinkIdx,
                             bool ignoreInexistent = false);
  ;

  /**
   * @brief Add atoms (incl. type, id etc.) with given positions to the
   * universe
   *
   * @param nrOfAtomsToAdd the nr. of atoms to add to the universe
   * @param atomType the type of the atoms to add
   * @return std::vector<size_t> the ids of the inserted atoms
   */
  std::vector<size_t> addXlinkAtoms(int nrOfAtomsToAdd,
                                    int atomType,
                                    Eigen::VectorXd coordinates);

  /**
   * @brief Add atoms (incl. random positions, id etc.) to the universe
   *
   * @param nrOfAtomsToAdd the nr. of atoms to add to the universe
   * @param atomType the type of the atoms to add
   * @return std::vector<size_t> the ids of the inserted atoms
   */
  std::vector<size_t> addXlinkAtoms(int nrOfAtomsToAdd,
                                    int atomType,
                                    bool whiteNoise = true);

  /**
   * @brief Generate positions randomly
   *
   * @param nSamples the number of positions to generate
   * @return Positions
   */
  Eigen::VectorXd generateRandomPositions(const int nSamples,
                                          const bool whiteNoise = true)
  {
    if (whiteNoise) {
      return this->generateRandomWhitePositions(nSamples);
    } else {
      return this->generateRandomBluePositions(nSamples);
    }
  }

  /**
   * @brief Generate positions randomly (can lead to clustering)
   *
   * @param nSamples the number of positions to generate
   * @return Positions
   */
  Eigen::VectorXd generateRandomWhitePositions(int nSamples);

  /**
   * @brief Generate positions randomly according to blue noise type sampling
   *
   * @param nSamples the number of positions to generate
   * @return Positions
   */
  Eigen::VectorXd generateRandomBluePositions(int nSamples);

  /**
   * @brief find an Atom in a collection with an objective distance
   *
   * @param from the Atom to start the distance from
   * @param desiredR02 the distance to target if possible
   * @param maxDistance the maximum distance to accept random matches from
   * @return long int the index in possiblePartners that matches best,
   * negative if none found
   */
  long int findAppropriateLink(size_t from,
                               const double desiredR02,
                               const double maxDistance);

  double evaluatePartnerProbability(
    size_t from,
    size_t to,
    const double normalisationFactorInExponential,
    const double maxDistance,
    bool respectXlinkChains = true);

  double getIdealCutoff() const;

  /**
   * @brief Re-initialize the neighbour list of the crosslinks
   *
   */
  void resetNeighbourList();

  /**
   * @brief Update the coordinates stored in the neighbour list
   *
   */
  void updateNeighbourListCoordinates();

  ///// Utility functions
  /**
   * @brief compute the distance between two crosslinks, given by their
   * indices
   *
   * @param i
   * @param j
   * @return double
   */
  double distanceBetween(size_t i, size_t j) const;

  Eigen::Vector3d getVectorBetween(size_t i, size_t j) const;

  double getDistance(double x1,
                     double y1,
                     double z1,
                     double x2,
                     double y2,
                     double z2) const;

  Eigen::Vector3d sampleCoordinatesWithinNBeadDistance(int nBeads);

  /**
   * @brief Copy the crosslinker coordinates into an Eigen::VectorXd
   *
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getCrosslinkCoordinates() const;

  /**
   * @brief Combine two strands into one by merging them
   *
   * @param strandIdx1 index of first strand
   * @param strandIdx2 index of second strand
   */
  void combineStrands(size_t strandIdx1, size_t strandIdx2);

  /**
   * @brief Calculate end-to-end distance probability for strand combination
   *
   * @param strandIdx1 index of first strand
   * @param strandIdx2 index of second strand
   * @param cInfinity statistical parameter
   * @return probability weight for this combination
   */
  double calculateStrandCombinationProbability(size_t strandIdx1,
                                               size_t strandIdx2,
                                               double cInfinity);
};
} // namespace pylimer_tools::utils
