#ifndef MC_UNIVERSE_GENERATOR_H
#define MC_UNIVERSE_GENERATOR_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../sim/MEHPForceRelaxation.h"
#include "../sim/MEHPUtilityStructures.h"
#include "../utils/BoolUtils.h"
#include "RandomWalker.h"
#include "StringUtils.h"
#include "VectorUtils.h"
#include <Eigen/Dense>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>/* isnan, sqrt */
#include <stack>
#include <string>
#include <vector>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

#include <random>

namespace pylimer_tools {
namespace utils {

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

    std::vector<long int> strandFrom;
    std::vector<long int> strandTo;
    std::vector<int> beadsInStrand;
    std::vector<int> strandBeadType;
    std::vector<double> beadDistanceInStrand;
    std::vector<double> meanSquaredBeadDistanceInStrand;
  };

  class MCUniverseGenerator
  {
  public:
    MCUniverseGenerator(const double Lx, const double Ly, const double Lz)
    {
      std::random_device rd{};
      this->rng = std::mt19937(rd());
      this->distX = std::uniform_real_distribution<double>(0.0, Lx);
      this->distY = std::uniform_real_distribution<double>(0.0, Ly);
      this->distZ = std::uniform_real_distribution<double>(0.0, Lz);
      this->box = pylimer_tools::entities::Box(Lx, Ly, Lz);
      this->setBeadDistance(0.965);
    }

    void setSeed(unsigned int seed) { this->rng.seed(seed); }

    void setBeadDistance(double newBeadDistance, bool updateMeanSquared = true)
    {
      INVALIDARG_EXP_IFN(newBeadDistance > 0, "Invalid mean bead distance");
      INVALIDARG_EXP_IFN(
        !(std::isnan(newBeadDistance) || std::isinf(newBeadDistance)),
        "Invalid mean bead distance");
      this->beadDistance = newBeadDistance;
      if (updateMeanSquared) {
        this->meanSquaredBeadDistance =
          (3. / 8.) * M_PI * SQUARE(this->beadDistance);
        RUNTIME_EXP_IFN(this->meanSquaredBeadDistance > 0,
                        "Invalid mean squared bead distance");
      }
    }

    double getConfiguredBeadDistance() const { return this->beadDistance; }

    void setMeanSquaredBeadDistance(double newMeanSquaredBeadDistance,
                                    bool updateMean = true)
    {
      INVALIDARG_EXP_IFN(newMeanSquaredBeadDistance > 0,
                         "Invalid mean squared bead distance");
      INVALIDARG_EXP_IFN(!(std::isnan(newMeanSquaredBeadDistance) ||
                           std::isinf(newMeanSquaredBeadDistance)),
                         "Invalid mean squared bead distance");
      this->meanSquaredBeadDistance = newMeanSquaredBeadDistance;
      if (updateMean) {
        this->beadDistance =
          std::sqrt(this->meanSquaredBeadDistance / ((3. / 8.) * M_PI));
        RUNTIME_EXP_IFN(this->beadDistance > 0, "Invalid mean bead distance");
      }
    }

    double getConfiguredMeanSquaredBeadDistance() const
    {
      return this->meanSquaredBeadDistance;
    }

    void configNrOfMCSteps(size_t newNrOfMCSteps)
    {
      this->nMcSteps = newNrOfMCSteps;
    }

    void configPrimaryLoopProbability(double newPrimaryLoopProbability)
    {
      INVALIDARG_EXP_IFN(newPrimaryLoopProbability >= 0,
                         "Invalid primary loop formation probability");
      this->primaryLoopProbability = newPrimaryLoopProbability;
    }

    void configSecondaryLoopProbability(double newSecondaryLoopProbability)
    {
      INVALIDARG_EXP_IFN(newSecondaryLoopProbability >= 0,
                         "Invalid secondary loop formation probability");
      this->secondaryLoopProbability = newSecondaryLoopProbability;
    }

    /**
     * @brief Get the Universe object after actually sampling strands' beads and
     * their positions.
     *
     * @return pylimer_tools::entities::Universe
     */
    pylimer_tools::entities::Universe getUniverse()
    {
      pylimer_tools::entities::Universe universe =
        pylimer_tools::entities::Universe(this->box);
      size_t nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
      int nrOfAtoms =
        nCrosslinks +
        std::reduce(this->simplifiedUniverse.beadsInStrand.begin(),
                    this->simplifiedUniverse.beadsInStrand.end(),
                    0);
      std::vector<int> zeros = initializeWithValue(nrOfAtoms, 0);
      std::vector<double> xs;
      xs.reserve(nrOfAtoms);
      std::vector<double> ys;
      ys.reserve(nrOfAtoms);
      std::vector<double> zs;
      zs.reserve(nrOfAtoms);
      std::vector<long int> ids;
      ids.reserve(nrOfAtoms);
      std::vector<int> types;
      types.reserve(nrOfAtoms);

      std::vector<long int> bondsFrom = {};
      std::vector<long int> bondsTo = {};

      // Add the cross-linkers
      long int currentId = 1;
      for (size_t i = 0; i < nCrosslinks; ++i) {
        ids.push_back(currentId);
        xs.push_back(this->simplifiedUniverse.xlinkX[i]);
        ys.push_back(this->simplifiedUniverse.xlinkY[i]);
        zs.push_back(this->simplifiedUniverse.xlinkZ[i]);
        types.push_back(this->simplifiedUniverse.xlinkTypes[i]);
        currentId += 1;
      }

      // Sample the strands
      assert(this->simplifiedUniverse.strandFrom.size() ==
             this->simplifiedUniverse.strandTo.size());
      double prevBeadDistance = this->beadDistance;
      double prevMeanSquaredBeadDistance = this->meanSquaredBeadDistance;
      for (size_t strandI = 0;
           strandI < this->simplifiedUniverse.strandFrom.size();
           ++strandI) {
        this->beadDistance =
          this->simplifiedUniverse.beadDistanceInStrand[strandI];
        this->meanSquaredBeadDistance =
          this->simplifiedUniverse.meanSquaredBeadDistanceInStrand[strandI];
        // sample the bead coordinates, depending on the type of strand
        Eigen::VectorXd coordinates;
        long int strandEnd1 = this->simplifiedUniverse.strandFrom[strandI];
        long int strandEnd2 = this->simplifiedUniverse.strandTo[strandI];
        int nBeadsInStrand = this->simplifiedUniverse.beadsInStrand[strandI];
        if (strandEnd1 < 0) {
          RUNTIME_EXP_IFN(
            strandEnd2 < 0,
            "if first end is not associated, expected second to be as well");
          coordinates = this->sampleFreeChainCoordinates(nBeadsInStrand);
        } else if (strandEnd2 < 0) {
          coordinates =
            this->sampleDanglingChainCoordinates(strandEnd1, nBeadsInStrand);
          bondsFrom.push_back(strandEnd1 + 1);
          bondsTo.push_back(currentId);
        } else {
          RUNTIME_EXP_IFN(
            strandEnd1 != -1 && strandEnd2 != -1,
            "Expected both ends to have an associated cross-linker");
          coordinates = this->sampleStrandCoordinates(
            strandEnd1, strandEnd2, nBeadsInStrand);
          bondsFrom.push_back(strandEnd1 + 1);
          bondsTo.push_back(currentId);
          bondsFrom.push_back(strandEnd2 + 1);
          bondsTo.push_back(currentId + nBeadsInStrand - 1);
        }

        // actually add the new beads to our list of things to add
        RUNTIME_EXP_IFN(coordinates.size() == 3 * nBeadsInStrand,
                        "Inconsistent coordinate size");
        RUNTIME_EXP_IFN(!coordinates.array().isNaN().any(),
                        "Coordinates contain NaN in strand " +
                          std::to_string(strandI + 1));
        for (size_t i = 0; i < nBeadsInStrand; ++i) {
          ids.push_back(currentId);
          if (i > 0) {
            bondsFrom.push_back(currentId - 1);
            bondsTo.push_back(currentId);
          }
          types.push_back(this->simplifiedUniverse.strandBeadType[strandI]);
          xs.push_back(coordinates(3 * i));
          ys.push_back(coordinates(3 * i + 1));
          zs.push_back(coordinates(3 * i + 2));
          currentId += 1;
        }
      }

      this->beadDistance = prevBeadDistance;
      this->meanSquaredBeadDistance = prevMeanSquaredBeadDistance;

      universe.addAtoms(ids, types, xs, ys, zs, zeros, zeros, zeros);
      universe.addBonds(bondsFrom, bondsTo);
      return universe;
    };

    /**
     * @brief Add free atoms with a specified type and functionality for
     * possible cross-linking later
     *
     * @param nrOfCrosslinkers
     * @param crosslinkerFunctionality
     * @param crossLinkerAtomType
     * @param whiteNoise
     */
    void addCrosslinkers(int nrOfCrosslinkers,
                         int crosslinkerFunctionality = 4,
                         int crossLinkerAtomType = 2,
                         bool whiteNoise = true)
    {
      int nCrosslinkerBefore = this->remainingCrossLinkerFunctionality.size();

      this->addXlinkAtoms(nrOfCrosslinkers, crossLinkerAtomType, whiteNoise);

      this->remainingCrossLinkerFunctionality.reserve(nCrosslinkerBefore +
                                                      nrOfCrosslinkers);
      this->originalNrOfAvailableCrosslinkSites +=
        crosslinkerFunctionality * nrOfCrosslinkers;
      this->nrOfAvailableCrosslinkSites +=
        crosslinkerFunctionality * nrOfCrosslinkers;
      for (size_t i = 0; i < nrOfCrosslinkers; ++i) {
        this->remainingCrossLinkerFunctionality.push_back(
          crosslinkerFunctionality);
      }
    };

    /**
     * @brief Add strands which contain cross-links at random positions
     *
     * @param nrOfStrands
     * @param beadsPerStrand
     * @param functionalizationProbability the probability of each bead being a
     * cross-link
     * @param crosslinkerFunctionality the functionality of the cross-linkers,
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
                                          bool whiteNoise = true)
    {
      INVALIDARG_EXP_IFN(nrOfStrands > 0, "");
      INVALIDARG_EXP_IFN(
        functionalizationProbability >= 0. &&
          functionalizationProbability <= 1.,
        "Functionalisation probability must be between 0 and 1.");
      INVALIDARG_EXP_IFN(crosslinkerFunctionality >= 1,
                         "Cross-linker functionality must be at least 1. Use "
                         "normal chains instead for lower functionalities.");
      this->validateInternalState();

      size_t nCrosslinksBefore = this->remainingCrossLinkerFunctionality.size();
      size_t nStrandsBefore = this->simplifiedUniverse.strandFrom.size();

      std::uniform_real_distribution<double> randomDist(0., 1.);
      size_t nCrosslinks = 0;
      size_t nEffectiveStrands = 0;

      for (size_t strandI = 0; strandI < nrOfStrands; ++strandI) {
        std::vector<int> partialStrandLengths;
        std::vector<size_t> crosslinkIndicesInStrand;
        long int lastSampledBead = -1;
        for (long int i = 0; i < beadsPerStrand[strandI]; ++i) {
          if (randomDist(this->rng) < functionalizationProbability) {
            // yes, we want to replace this bead `i` with a cross-link
            size_t currentSpringIdx = nStrandsBefore + nEffectiveStrands;
            long int remainingLength =
              i - (lastSampledBead >= 0 ? (lastSampledBead + 1) : 0);
            if (i > 0) {
              this->addStrands(1, remainingLength, strandAtomType);
              nEffectiveStrands += 1;
              this->linkStrandToCrosslink(
                currentSpringIdx, nCrosslinksBefore + nCrosslinks, true);
            }

            if (lastSampledBead > 0) {
              this->linkStrandToCrosslink(
                currentSpringIdx - 1, nCrosslinksBefore + nCrosslinks, true);
            }

            nCrosslinks += 1;
            lastSampledBead = i;
          }
        }

        if (lastSampledBead < beadsPerStrand[strandI] - 1) {
          // add dangling spring
          size_t currentSpringIdx = nStrandsBefore + nEffectiveStrands;

          long int remainingLength =
            beadsPerStrand[strandI] -
            (lastSampledBead >= 0 ? (lastSampledBead + 1) : 0);
          this->addStrands(1, remainingLength, strandAtomType);
          nEffectiveStrands += 1;
          if (lastSampledBead > 0) {
            this->linkStrandToCrosslink(
              currentSpringIdx, nCrosslinksBefore + nCrosslinks, true);
          }
        }
      }

      // then, add the cross-link atoms
      this->addCrosslinkers(
        nCrosslinks, crosslinkerFunctionality, crosslinkerAtomType, whiteNoise);

      for (size_t newStrandIdx = nStrandsBefore;
           newStrandIdx < nStrandsBefore + nEffectiveStrands;
           ++newStrandIdx) {
        long int from = this->simplifiedUniverse.strandFrom[newStrandIdx];
        long int to = this->simplifiedUniverse.strandTo[newStrandIdx];
        if (from >= 0) {
          this->simplifiedUniverse.strandsOfXlink[from].push_back(newStrandIdx);
        }
        if (to >= 0) {
          this->simplifiedUniverse.strandsOfXlink[to].push_back(newStrandIdx);
        }

        if (from >= 0 && to >= 0) {
          // finally, make sure that the coordinates of two subsequent
          // cross-links match the required length of the respective strand

          std::normal_distribution<double> otherEndCoordinateDist =
            std::normal_distribution<double>(
              0.,
              std::sqrt(
                static_cast<double>(
                  this->simplifiedUniverse.beadsInStrand[newStrandIdx] + 1) *
                this->meanSquaredBeadDistance));

          this->simplifiedUniverse.xlinkX[to] =
            this->simplifiedUniverse.xlinkX[from] +
            otherEndCoordinateDist(this->rng);
          this->simplifiedUniverse.xlinkY[to] =
            this->simplifiedUniverse.xlinkY[from] +
            otherEndCoordinateDist(this->rng);
          this->simplifiedUniverse.xlinkZ[to] =
            this->simplifiedUniverse.xlinkZ[from] +
            otherEndCoordinateDist(this->rng);

          // validate the distances
          double distance = this->distanceBetween(from, to);
          RUNTIME_EXP_IFN(
            distance < this->simplifiedUniverse.beadsInStrand[newStrandIdx] *
                         this->beadDistance,
            "Distance adjustment does not seem to be correct.");
        }
      }

      this->validateInternalState();
    }

    /**
     * @brief add strands with cross-links as ends
     *
     * @param nrOfCrosslinkStrands
     * @param crosslinkerFunctionality incl. the connection to the strand
     * @param crosslinkerAtomType
     * @param beadsPerStrand excl. cross-links
     * @param strandAtomType
     * @param whiteNoise
     */
    void addCrosslinkStrands(int nrOfCrosslinkStrands,
                             std::vector<int> beadsPerStrand,
                             int crosslinkerFunctionality = 4,
                             int crosslinkerAtomType = 2,
                             int strandAtomType = 2,
                             bool whiteNoise = true)
    {
      INVALIDARG_EXP_IFN(nrOfCrosslinkStrands > 0, "");
      INVALIDARG_EXP_IFN(crosslinkerFunctionality >= 2,
                         "Cross-linker functionality must be at least 2. Use "
                         "solvent chains instead for lower functionalities.");
      this->validateInternalState();

      Eigen::VectorXd firstLinkCoordinates =
        this->generateRandomPositions(nrOfCrosslinkStrands);
      assert(firstLinkCoordinates.size() == 3 * nrOfCrosslinkStrands);
      firstLinkCoordinates.conservativeResize(2 * 3 * nrOfCrosslinkStrands);

      for (size_t i = 0; i < nrOfCrosslinkStrands; ++i) {
        std::normal_distribution<double> otherEndCoordinateDist =
          std::normal_distribution<double>(
            0.,
            std::sqrt(static_cast<double>(beadsPerStrand[i] + 1) *
                      this->meanSquaredBeadDistance));

        for (size_t dir = 0; dir < 3; ++dir) {
          firstLinkCoordinates(2 * 3 * i + dir) =
            firstLinkCoordinates(3 * i + dir) +
            otherEndCoordinateDist(this->rng);
        }
      }

      size_t nCrosslinksBefore = this->remainingCrossLinkerFunctionality.size();

      this->addXlinkAtoms(
        2 * nrOfCrosslinkStrands, crosslinkerAtomType, firstLinkCoordinates);

      size_t nStrandsBefore = this->simplifiedUniverse.strandFrom.size();

      this->addStrands(nrOfCrosslinkStrands, beadsPerStrand, strandAtomType);

      this->originalNrOfAvailableCrosslinkSites +=
        crosslinkerFunctionality * nrOfCrosslinkStrands * 2;
      this->nrOfAvailableCrosslinkSites +=
        (crosslinkerFunctionality - 1) * nrOfCrosslinkStrands * 2;

      for (size_t i = 0; i < nrOfCrosslinkStrands; ++i) {
        // actually link these strands
        this->simplifiedUniverse.strandFrom[nStrandsBefore + i] =
          nCrosslinksBefore + i;
        this->simplifiedUniverse.strandTo[nStrandsBefore + i] =
          nCrosslinksBefore + i + nrOfCrosslinkStrands;
        // remember remaining cross-link functionality
        this->remainingCrossLinkerFunctionality.push_back(
          crosslinkerFunctionality - 1);
        this->remainingCrossLinkerFunctionality.push_back(
          crosslinkerFunctionality - 1);
      }

      this->validateInternalState();
    }

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
                          bool whiteNoise = true)
    {
      for (size_t i = 0; i < nrOfSolventChains; ++i) {
        this->simplifiedUniverse.strandFrom.push_back(EMPTY_BACKGROUND);
        this->simplifiedUniverse.strandTo.push_back(EMPTY_BACKGROUND);
        this->simplifiedUniverse.strandBeadType.push_back(solventAtomType);
        this->simplifiedUniverse.beadsInStrand.push_back(chainLength);
        this->simplifiedUniverse.beadDistanceInStrand.push_back(
          this->beadDistance);
        this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.push_back(
          this->meanSquaredBeadDistance);
      }
    };

    /**
     * @brief Add multiple monofunctional strands with specified bead types,
     * link them to cross-links
     *
     * @param nrOfStrands
     * @param beadsPerChains
     * @param strandAtomType
     */
    void addMonofunctionalStrands(int nrOfStrands,
                                  std::vector<int> beadsPerChains,
                                  int strandAtomType = 1)
    {
      INVALIDARG_EXP_IFN(beadsPerChains.size() == nrOfStrands,
                         "Nr of strands (" + std::to_string(nrOfStrands) +
                           ") must be equal to the number "
                           "of chainLengths (" +
                           std::to_string(beadsPerChains.size()) +
                           ") provided.");

      for (size_t i = 0; i < nrOfStrands; ++i) {
        this->simplifiedUniverse.strandBeadType.push_back(strandAtomType);
        this->simplifiedUniverse.beadsInStrand.push_back(beadsPerChains[i]);
        this->simplifiedUniverse.beadDistanceInStrand.push_back(
          this->beadDistance);
        this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.push_back(
          this->meanSquaredBeadDistance);
        this->simplifiedUniverse.strandFrom.push_back(UNCONNECTED);
        this->simplifiedUniverse.strandTo.push_back(EMPTY_BACKGROUND);
      }
    }

    void addMonofunctionalStrands(int nrOfStrands,
                                  int chainLength,
                                  int strandAtomType = 1)
    {
      const std::vector<int> chainLengths =
        pylimer_tools::utils::initializeWithValue<int>(nrOfStrands,
                                                       chainLength);
      return this->addMonofunctionalStrands(
        nrOfStrands, chainLengths, strandAtomType);
    }

    void addStrands(const int nrOfStrands,
                    const std::vector<int> beadsPerChains,
                    const int strandAtomType = 1)
    {
      INVALIDARG_EXP_IFN(
        beadsPerChains.size() == nrOfStrands,
        "Inconsistent nr of strands and nr of beads per strand given.");

      for (size_t strandIdx = 0; strandIdx < nrOfStrands; ++strandIdx) {
        this->simplifiedUniverse.strandBeadType.push_back(strandAtomType);
        this->simplifiedUniverse.beadsInStrand.push_back(
          beadsPerChains[strandIdx]);
        this->simplifiedUniverse.beadDistanceInStrand.push_back(
          this->beadDistance);
        this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.push_back(
          this->meanSquaredBeadDistance);
        this->simplifiedUniverse.strandFrom.push_back(UNCONNECTED);
        this->simplifiedUniverse.strandTo.push_back(UNCONNECTED);
      }
    }

    void addStrands(int nrOfStrands, int chainLength, int strandAtomType = 1)
    {
      const std::vector<int> chainLengths =
        pylimer_tools::utils::initializeWithValue<int>(nrOfStrands,
                                                       chainLength);
      return this->addStrands(nrOfStrands, chainLengths, strandAtomType);
    };

    /**
     * @brief Add strands, link them to the cross-links, stop when the callback
     * says so
     *
     * @param stopLinking a callback to indicated whether to stop linking
     * @param cInfinity `C_\infty` for `<R_ee^2>_0` from `N` and `<b^2>`
     */
    void linkStrandsCallback(
      std::function<BackTrackStatus(const MCUniverseGenerator&, size_t)>
        linkingController,
      double cInfinity = 1.)
    {
      this->validateInternalState();

      // prepare sampling of partners
      long int nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
      int nrOfStrandsAdded = 0;

      double currentCrosslinkerConversion =
        1. - (static_cast<double>(this->nrOfAvailableCrosslinkSites) /
              static_cast<double>(this->originalNrOfAvailableCrosslinkSites));
      double conversionPerBond =
        1. / (static_cast<double>(this->originalNrOfAvailableCrosslinkSites));

      size_t nrOfStrands = this->simplifiedUniverse.strandFrom.size();
      std::vector<size_t> availableStrandEnds;
      availableStrandEnds.reserve(2 * nrOfStrands);
      for (size_t i = 0; i < nrOfStrands; ++i) {
        // each strand is available with two ends, assuming they are unconnected
        if (this->simplifiedUniverse.strandFrom[i] == UNCONNECTED) {
          availableStrandEnds.push_back(i);
        }
        if (this->simplifiedUniverse.strandTo[i] == UNCONNECTED) {
          availableStrandEnds.push_back(i);
        }
      }
      availableStrandEnds.shrink_to_fit();
      std::shuffle(
        availableStrandEnds.begin(), availableStrandEnds.end(), this->rng);

      std::vector<size_t> availableCrosslinkSites;
      availableCrosslinkSites.reserve(this->nrOfAvailableCrosslinkSites);
      for (size_t i = 0; i < nCrosslinks; ++i) {
        for (long int s = 0; s < this->remainingCrossLinkerFunctionality[i];
             ++s) {
          availableCrosslinkSites.push_back(i);
        }
      }
      std::shuffle(availableCrosslinkSites.begin(),
                   availableCrosslinkSites.end(),
                   this->rng);

      std::stack<size_t> removedCrosslinkSites;
      // removedCrosslinkSites.reserve(availableCrosslinkSites.size());
      const double timesNForR02 = this->meanSquaredBeadDistance * cInfinity;

      // link one strand at a time until we reach the target conversion
      for (size_t sampleIdx = 0; sampleIdx < availableStrandEnds.size();
           ++sampleIdx) {
        BackTrackStatus status =
          linkingController(*this, availableStrandEnds.size() - sampleIdx);
        if (status == BackTrackStatus::STOP) {
          break;
        } else if (status == BackTrackStatus::TRACK_FORWARD) {
          size_t strandIdx = availableStrandEnds[sampleIdx];
          RUNTIME_EXP_IFN(
            this->simplifiedUniverse.strandTo[strandIdx] < 0,
            "Expected second strand end to be free, got " +
              std::to_string(this->simplifiedUniverse.strandTo[strandIdx]) +
              " for strand " + std::to_string(strandIdx) + ".");

          if (this->simplifiedUniverse.strandFrom[strandIdx] >= 0) {
            RUNTIME_EXP_IFN(
              this->simplifiedUniverse.strandTo[strandIdx] == UNCONNECTED,
              "Expected second strand end to be free, got " +
                std::to_string(this->simplifiedUniverse.strandTo[strandIdx]) +
                " for strand " + std::to_string(strandIdx) + ".");
            // we don't have free cross-link choice
            // find one that follows the desired end-to-end distribution
            size_t partnerCrosslinker = this->findAppropriateLink(
              this->simplifiedUniverse.strandFrom[strandIdx],
              static_cast<double>(
                this->simplifiedUniverse.beadsInStrand[strandIdx] + 1) *
                timesNForR02,
              -1. // beadsPerChains[strandIdx] * this->beadDistance
            );

            this->simplifiedUniverse.strandTo[strandIdx] = partnerCrosslinker;
            this->remainingCrossLinkerFunctionality[partnerCrosslinker] -= 1;
            this->simplifiedUniverse.strandsOfXlink[partnerCrosslinker]
              .push_back(strandIdx);
            this->nrOfAvailableCrosslinkSites -= 1;
          } else {
            // otherwise, randomly choose a free cross-link
            long int crosslinkIdxIdx;
            do {
              // find the next "available" cross-linker
              crosslinkIdxIdx = availableCrosslinkSites.back();
              removedCrosslinkSites.push(crosslinkIdxIdx);
              availableCrosslinkSites.pop_back();
            } while (this->remainingCrossLinkerFunctionality[crosslinkIdxIdx] <
                       1 &&
                     availableCrosslinkSites.size() > 0);

            if (availableCrosslinkSites.size() == 0) {
              std::cerr << "No more cross-link sites available." << std::endl;
              break;
            }

            this->simplifiedUniverse.strandFrom[strandIdx] = crosslinkIdxIdx;
            this->remainingCrossLinkerFunctionality[crosslinkIdxIdx] -= 1;
            this->simplifiedUniverse.strandsOfXlink[crosslinkIdxIdx].push_back(
              strandIdx);
            this->nrOfAvailableCrosslinkSites -= 1;
          }
        } else {
          assert(status == BackTrackStatus::TRACK_BACKWARD);
          // track backward -> reset the last link done
          sampleIdx -= 1;
          long int strandIdx = availableStrandEnds[sampleIdx];
          // this strand has been assigned a partner last step, add it again
          long int linkedXlink = -1;
          if (this->simplifiedUniverse.strandTo[strandIdx] >= 0) {
            linkedXlink = this->simplifiedUniverse.strandTo[strandIdx];
            this->simplifiedUniverse.strandTo[strandIdx] = UNCONNECTED;
          } else {
            linkedXlink = this->simplifiedUniverse.strandFrom[strandIdx];
            this->simplifiedUniverse.strandFrom[strandIdx] = UNCONNECTED;
            size_t lastRemoved;
            // need to re-fill the available cross-link sites for the
            // first strand end
            do {
              lastRemoved = removedCrosslinkSites.top();
              removedCrosslinkSites.pop();
              availableCrosslinkSites.push_back(lastRemoved);
            } while (lastRemoved != linkedXlink);
          }

          this->remainingCrossLinkerFunctionality[linkedXlink] += 1;
          this->simplifiedUniverse.strandsOfXlink[linkedXlink].pop_back();
          this->nrOfAvailableCrosslinkSites += 1;

          // one step back more, since the iteration will iterate anyway
          sampleIdx -= 1;
        }
      }

      this->validateInternalState();
    }

    /**
     * @brief Add strands in between the cross-linkers, link them as appropriate
     *
     * @param targetCrossLinkerConversion "p", the target conversion of the
     * cross-linkers
     * @param cInfinity `C_\infty` for `<R_ee^2>_0` from `N` and `<b^2>`
     */
    void linkStrandsToConversion(const double targetCrossLinkerConversion,
                                 const double cInfinity = 1.)
    {
      const double conversionPerBond =
        (1.0) /
        (static_cast<double>(this->originalNrOfAvailableCrosslinkSites));
      const double currentCrosslinkerConversion =
        1.0 - (static_cast<double>(this->nrOfAvailableCrosslinkSites) /
               static_cast<double>(this->originalNrOfAvailableCrosslinkSites));

      INVALIDARG_EXP_IFN(
        targetCrossLinkerConversion >= currentCrosslinkerConversion &&
          targetCrossLinkerConversion <= 1.0,
        "Cross-linker conversion must be between " +
          std::to_string(currentCrosslinkerConversion) + " and 1, got " +
          std::to_string(targetCrossLinkerConversion) + ".");

      long int potentialNewBonds = 0;
      for (size_t i = 0; i < this->simplifiedUniverse.strandFrom.size(); ++i) {
        if (this->simplifiedUniverse.strandFrom[i] == UNCONNECTED) {
          potentialNewBonds += 1;
        }
        if (this->simplifiedUniverse.strandTo[i] == UNCONNECTED) {
          potentialNewBonds += 1;
        }
      }
      if (currentCrosslinkerConversion + potentialNewBonds * conversionPerBond <
          targetCrossLinkerConversion) {
        throw std::invalid_argument(
          "A cross-linker conversion of " +
          std::to_string(targetCrossLinkerConversion) +
          " is not reachable with this nr of strands.");
      }

      const long int targetNrOfAvailableCrosslinkSites = std::round(
        (1.0 - targetCrossLinkerConversion) *
        static_cast<double>(this->originalNrOfAvailableCrosslinkSites));
      const double timesNForR02 = this->meanSquaredBeadDistance * cInfinity;

      this->linkStrandsCallback(
        [targetNrOfAvailableCrosslinkSites](const MCUniverseGenerator& gen,
                                            size_t nStrandsRemaining) {
          if (gen.nrOfAvailableCrosslinkSites <=
              targetNrOfAvailableCrosslinkSites) {
            return BackTrackStatus::STOP;
          };
          return BackTrackStatus::TRACK_FORWARD;
        },
        cInfinity);
    }

    /**
     * @brief Add strands in between the cross-linkers, link them as appropriate
     *
     * @param targetSolubleFraction "w_sol", the target soluble fraction
     * @param cInfinity `C_\infty` for `<R_ee^2>_0` from `N` and `<b^2>`
     */
    void linkStrandsToSolubleFraction(double targetSolubleFraction,
                                      double cInfinity = 1.)
    {
      INVALIDARG_EXP_IFN(targetSolubleFraction >= 0. &&
                           targetSolubleFraction <= 1.,
                         "Soluble fraction must be between 0 and 1, got " +
                           std::to_string(targetSolubleFraction) + ".");

      size_t nAtomsTotal =
        std::reduce(this->simplifiedUniverse.beadsInStrand.begin(),
                    this->simplifiedUniverse.beadsInStrand.end(),
                    0) +
        this->simplifiedUniverse.xlinkTypes.size();

      BackTrackStatus status = BackTrackStatus::TRACK_FORWARD;
      long int nSteps = this->simplifiedUniverse.xlinkTypes.size();
      long int lastStep = 0;
      long int currentStep = 0;

      this->linkStrandsCallback(
        [targetSolubleFraction,
         nAtomsTotal,
         &status,
         &nSteps,
         &lastStep,
         &currentStep](const MCUniverseGenerator& gen,
                       size_t nStrandsRemaining) {
          currentStep += 1;
          nSteps = std::max<long int>(
            std::min<long int>(
              nSteps,
              static_cast<long int>(gen.simplifiedUniverse.xlinkTypes.size())),
            (long int)1);
          // make large step without force relaxation
          if ((currentStep < lastStep + nSteps) &&
              // but only if the last strand will not be reached with this large
              // step
              (nStrandsRemaining > 2)) {
            return status;
          }

          lastStep = currentStep;
          pylimer_tools::sim::mehp::Network frNet =
            gen.convertToForceRelaxationNetwork();

          // actually start force relaxation
          pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
            pylimer_tools::sim::mehp::MEHPForceRelaxation(frNet);
          forceRelaxer.configAssumeBoxLargeEnough(true);

          while (forceRelaxer.suggestsRerun()) {
            forceRelaxer.runForceRelaxation("LD_MMA", 5000, 1e-12, 1e-9);
          }

          // finally, calculate the soluble fraction
          double solubleFraction =
            1. - forceRelaxer.countActiveClusteredAtoms() /
                   static_cast<double>(nAtomsTotal);
          std::cout << "Got w_sol = " << solubleFraction << " at step "
                    << currentStep << " (+" << nSteps << ") with "
                    << nStrandsRemaining << " strands remaining." << std::endl;
          if (solubleFraction > targetSolubleFraction) {
            if (status == BackTrackStatus::TRACK_BACKWARD) {
              nSteps /= 2;
            }
            status = BackTrackStatus::TRACK_FORWARD;
            return status;
          } else if (solubleFraction == targetSolubleFraction) {
            return BackTrackStatus::STOP;
          } else {
            if (nSteps == 1) {
              return BackTrackStatus::STOP;
            }
            status = BackTrackStatus::TRACK_BACKWARD;
            nSteps /= 2;
            return BackTrackStatus::TRACK_BACKWARD;
          };
        },
        cInfinity);
    }

    /**
     * @brief Remove a strand from the current state
     *
     * @param strandIdx
     */
    void removeStrand(size_t strandIdx)
    {
      INVALIDARG_EXP_IFN(strandIdx < this->simplifiedUniverse.strandFrom.size(),
                         "Strand to be removed is out of range.");

      if (this->simplifiedUniverse.strandFrom[strandIdx] >= 0) {
        // remove this strand from the cross-linker
        size_t xlink = this->simplifiedUniverse.strandFrom[strandIdx];
        pylimer_tools::utils::removeIfContained<long int>(
          this->simplifiedUniverse.strandsOfXlink[xlink], strandIdx);
        this->remainingCrossLinkerFunctionality[xlink] += 1;
        this->nrOfAvailableCrosslinkSites += 1;
      }
      if (this->simplifiedUniverse.strandTo[strandIdx] >= 0) {
        // remove this strand from the cross-linker
        size_t xlink = this->simplifiedUniverse.strandTo[strandIdx];
        pylimer_tools::utils::removeIfContained<long int>(
          this->simplifiedUniverse.strandsOfXlink[xlink], strandIdx);
        this->remainingCrossLinkerFunctionality[xlink] += 1;
        this->nrOfAvailableCrosslinkSites += 1;
      }

      this->simplifiedUniverse.strandFrom.erase(
        this->simplifiedUniverse.strandFrom.begin() + strandIdx);
      this->simplifiedUniverse.strandTo.erase(
        this->simplifiedUniverse.strandTo.begin() + strandIdx);
      this->simplifiedUniverse.beadsInStrand.erase(
        this->simplifiedUniverse.beadsInStrand.begin() + strandIdx);
      this->simplifiedUniverse.strandBeadType.erase(
        this->simplifiedUniverse.strandBeadType.begin() + strandIdx);
      this->simplifiedUniverse.beadDistanceInStrand.erase(
        this->simplifiedUniverse.beadDistanceInStrand.begin() + strandIdx);
      this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.erase(
        this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.begin() +
        strandIdx);

      // renumber the link from cross-links to strands
      for (size_t i = 0; i < this->simplifiedUniverse.xlinkTypes.size(); ++i) {
        for (size_t j = 0;
             j < this->simplifiedUniverse.strandsOfXlink[i].size();
             ++j) {
          assert(this->simplifiedUniverse.strandsOfXlink[i][j] != strandIdx);
          if (this->simplifiedUniverse.strandsOfXlink[i][j] > strandIdx) {
            this->simplifiedUniverse.strandsOfXlink[i][j] -= 1;
          }
        }
      }
    }

    /**
     * @brief Remove a cross-link from the current state
     *
     * @param crosslinkIdx
     */
    void removeCrosslink(size_t crosslinkIdx)
    {
      INVALIDARG_EXP_IFN(crosslinkIdx <
                           this->simplifiedUniverse.xlinkTypes.size(),
                         "Crosslink to be removed is out of range.");

      // re-connect the strands to account for the removed cross-link
      for (size_t i = 0; i < this->simplifiedUniverse.strandFrom.size(); ++i) {
        INVALIDARG_EXP_IFN(
          this->simplifiedUniverse.strandFrom[i] != crosslinkIdx &&
            this->simplifiedUniverse.strandTo[i] != crosslinkIdx,
          "The to-be removed cross-link is still connected.");
        if (this->simplifiedUniverse.strandFrom[i] > crosslinkIdx) {
          this->simplifiedUniverse.strandFrom[i] -= 1;
        }
        if (this->simplifiedUniverse.strandTo[i] > crosslinkIdx) {
          this->simplifiedUniverse.strandTo[i] -= 1;
        }
      }

      // adjust sum of available crosslink sites
      size_t remainingFunctionality =
        this->remainingCrossLinkerFunctionality[crosslinkIdx];
      size_t originalFunctionality =
        remainingFunctionality +
        this->simplifiedUniverse.strandsOfXlink[crosslinkIdx].size();
      this->nrOfAvailableCrosslinkSites -= remainingFunctionality;
      this->originalNrOfAvailableCrosslinkSites -= originalFunctionality;

      this->simplifiedUniverse.xlinkTypes.erase(
        this->simplifiedUniverse.xlinkTypes.begin() + crosslinkIdx);
      this->simplifiedUniverse.xlinkX.erase(
        this->simplifiedUniverse.xlinkX.begin() + crosslinkIdx);
      this->simplifiedUniverse.xlinkY.erase(
        this->simplifiedUniverse.xlinkY.begin() + crosslinkIdx);
      this->simplifiedUniverse.xlinkZ.erase(
        this->simplifiedUniverse.xlinkZ.begin() + crosslinkIdx);
      this->simplifiedUniverse.strandsOfXlink.erase(
        this->simplifiedUniverse.strandsOfXlink.begin() + crosslinkIdx);

      this->remainingCrossLinkerFunctionality.erase(
        this->remainingCrossLinkerFunctionality.begin() + crosslinkIdx);
    }

    /**
     * @brief Remove soluble parts from the current state
     *
     * @param rescale whether to rescale the box and coordinates, keeping the
     * density constant
     */
    void removeSolubleFraction(bool rescale = true)
    {
      size_t nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
      size_t nAtomsTotal =
        std::reduce(this->simplifiedUniverse.beadsInStrand.begin(),
                    this->simplifiedUniverse.beadsInStrand.end(),
                    0) +
        nCrosslinks;

      // first, remove already the strands that we know will not be relevant
      for (long int i = this->simplifiedUniverse.strandFrom.size() - 1; i >= 0;
           --i) {
        if (this->simplifiedUniverse.strandFrom[i] < 0 &&
            this->simplifiedUniverse.strandTo[i] < 0) {
          this->removeStrand(i);
        }
      }

      pylimer_tools::sim::mehp::Network forceRelaxationNetwork =
        this->convertToForceRelaxationNetwork();
      pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
        pylimer_tools::sim::mehp::MEHPForceRelaxation(forceRelaxationNetwork);

      while (forceRelaxer.suggestsRerun()) {
        forceRelaxer.runForceRelaxation("LD_MMA", 5000, 1e-12, 1e-9);
      }

      auto activeComponents = forceRelaxer.findClusteredToActive();
      Eigen::ArrayXb activeSprings = activeComponents.first;
      Eigen::ArrayXb activeNodes = activeComponents.second;

      RUNTIME_EXP_IFN(activeNodes.count() > 0,
                      "No active nodes found during force relaxation.");

      // now that we know which springs and nodes are soluble,
      // we can remove them
      RUNTIME_EXP_IFN(activeSprings.size() ==
                        this->simplifiedUniverse.strandFrom.size(),
                      "Number of springs does not match the number of strands. "
                      "Therefore, the mapping would be incorrect.");

      for (long int i = activeSprings.size() - 1; i >= 0; --i) {
        if (!activeSprings[i]) {
          this->removeStrand(i);
        }
      }
      // then, remove the w_sol cross-links
      for (long int i = this->simplifiedUniverse.xlinkTypes.size() - 1; i >= 0;
           --i) {
        assert(forceRelaxationNetwork.oldAtomIds(i) == i);
        if (!activeNodes[i]) {
          this->removeCrosslink(i);
        }
      }

      // finally, rescale the box if requested
      if (rescale) {
        size_t newNrOfAtoms =
          std::reduce(this->simplifiedUniverse.beadsInStrand.begin(),
                      this->simplifiedUniverse.beadsInStrand.end(),
                      0) +
          this->simplifiedUniverse.xlinkTypes.size();
        double scalingFactor = std::cbrt(static_cast<double>(newNrOfAtoms) /
                                         static_cast<double>(nAtomsTotal));
        this->box =
          pylimer_tools::entities::Box(scalingFactor * this->box.getLx(),
                                       scalingFactor * this->box.getLy(),
                                       scalingFactor * this->box.getLz());
        for (size_t i = 0; i < this->simplifiedUniverse.xlinkTypes.size();
             ++i) {
          this->simplifiedUniverse.xlinkX[i] *= scalingFactor;
          this->simplifiedUniverse.xlinkY[i] *= scalingFactor;
          this->simplifiedUniverse.xlinkZ[i] *= scalingFactor;
        }
      }
    }

    /**
     * @brief Convert the current simplified universe into a force relaxation
     * network
     *
     * @return pylimer_tools::sim::mehp::Network
     */
    pylimer_tools::sim::mehp::Network convertToForceRelaxationNetwork() const
    {
      pylimer_tools::sim::mehp::Network forceRelaxationNetwork;
      forceRelaxationNetwork.L[0] = this->box.getLx();
      forceRelaxationNetwork.L[1] = this->box.getLy();
      forceRelaxationNetwork.L[2] = this->box.getLz();
      forceRelaxationNetwork.vol = this->box.getVolume();

      const size_t nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
      size_t nDanglingEnds = 0;
      std::vector<size_t> danglingEndPartners = {};
      for (size_t i = 0; i < this->simplifiedUniverse.strandFrom.size(); ++i) {
        if (this->simplifiedUniverse.strandFrom[i] >= 0 &&
            this->simplifiedUniverse.strandTo[i] < 0) {
          nDanglingEnds += 1;
          danglingEndPartners.push_back(this->simplifiedUniverse.strandFrom[i]);
        }
      }
      forceRelaxationNetwork.nrOfNodes = nCrosslinks + nDanglingEnds;
      forceRelaxationNetwork.oldAtomIds =
        Eigen::ArrayXi::Zero(forceRelaxationNetwork.nrOfNodes);
      forceRelaxationNetwork.coordinates =
        Eigen::VectorXd(forceRelaxationNetwork.nrOfNodes * 3);
      forceRelaxationNetwork.springIndicesOfLinks.reserve(nCrosslinks);
      for (size_t crosslinkIdx = 0; crosslinkIdx < nCrosslinks;
           ++crosslinkIdx) {
        forceRelaxationNetwork.coordinates(3 * crosslinkIdx + 0) =
          this->simplifiedUniverse.xlinkX[crosslinkIdx];
        forceRelaxationNetwork.coordinates(3 * crosslinkIdx + 1) =
          this->simplifiedUniverse.xlinkY[crosslinkIdx];
        forceRelaxationNetwork.coordinates(3 * crosslinkIdx + 2) =
          this->simplifiedUniverse.xlinkZ[crosslinkIdx];
        forceRelaxationNetwork.oldAtomIds(crosslinkIdx) = crosslinkIdx;
      }
      for (size_t danglingEndIdx = 0; danglingEndIdx < nDanglingEnds;
           ++danglingEndIdx) {
        size_t nodeIdx = nCrosslinks + danglingEndIdx;
        // collapse dangling ends already to connected cross-link
        forceRelaxationNetwork.coordinates(3 * nodeIdx + 0) =
          this->simplifiedUniverse.xlinkX[danglingEndPartners[danglingEndIdx]];
        forceRelaxationNetwork.coordinates(3 * nodeIdx + 1) =
          this->simplifiedUniverse.xlinkY[danglingEndPartners[danglingEndIdx]];
        forceRelaxationNetwork.coordinates(3 * nodeIdx + 2) =
          this->simplifiedUniverse.xlinkZ[danglingEndPartners[danglingEndIdx]];
      }
      for (size_t i = 0; i < forceRelaxationNetwork.nrOfNodes; ++i) {
        std::vector<size_t> empty = {};
        forceRelaxationNetwork.springIndicesOfLinks.push_back(empty);
      }

      size_t nSpringEstimate = this->simplifiedUniverse.strandFrom.size();
      forceRelaxationNetwork.springIndexA =
        Eigen::ArrayXi::Zero(nSpringEstimate);
      forceRelaxationNetwork.springIndexB =
        Eigen::ArrayXi::Zero(nSpringEstimate);
      forceRelaxationNetwork.springsContourLength =
        Eigen::VectorXd::Zero(nSpringEstimate);
      size_t newSpringIdx = 0;
      size_t handledDanglingIdx = 0;
      // we omit all free strands
      for (size_t springIdx = 0;
           springIdx < this->simplifiedUniverse.strandFrom.size();
           ++springIdx) {
        long int from = this->simplifiedUniverse.strandFrom[springIdx];
        long int to = this->simplifiedUniverse.strandTo[springIdx];
        if (from >= 0) {
          if (to < 0) {
            to = nCrosslinks + handledDanglingIdx;
            handledDanglingIdx += 1;
            forceRelaxationNetwork.springsContourLength(newSpringIdx) =
              this->simplifiedUniverse.beadsInStrand[springIdx];
          } else {
            forceRelaxationNetwork.springsContourLength(newSpringIdx) =
              this->simplifiedUniverse.beadsInStrand[springIdx] + 1;
          }
          forceRelaxationNetwork.springIndexA(newSpringIdx) = from;
          forceRelaxationNetwork.springIndexB(newSpringIdx) = to;
          forceRelaxationNetwork.springIndicesOfLinks[from].push_back(
            newSpringIdx);
          forceRelaxationNetwork.springIndicesOfLinks[to].push_back(
            newSpringIdx);
          newSpringIdx += 1;
        }
      }
      assert(handledDanglingIdx == nDanglingEnds);
      const size_t nActualSprings = newSpringIdx;
      forceRelaxationNetwork.springIndexA.conservativeResize(nActualSprings);
      forceRelaxationNetwork.springIndexB.conservativeResize(nActualSprings);
      forceRelaxationNetwork.springsContourLength.conservativeResize(
        nActualSprings);
      forceRelaxationNetwork.springBoxOffset =
        Eigen::VectorXd::Zero(3 * nActualSprings);
      forceRelaxationNetwork.nrOfSprings = nActualSprings;

      forceRelaxationNetwork.springCoordinateIndexA =
        Eigen::ArrayXi::Zero(3 * nActualSprings);
      forceRelaxationNetwork.springCoordinateIndexB =
        Eigen::ArrayXi::Zero(3 * nActualSprings);
      for (size_t i = 0; i < nActualSprings; ++i) {
        for (size_t dir = 0; dir < 3; ++dir) {
          forceRelaxationNetwork.springCoordinateIndexA(3 * i + dir) =
            forceRelaxationNetwork.springIndexA(i) * 3 + dir;
          forceRelaxationNetwork.springCoordinateIndexB(3 * i + dir) =
            forceRelaxationNetwork.springIndexB(i) * 3 + dir;
        }
      }

      forceRelaxationNetwork.meanSpringContourLength =
        forceRelaxationNetwork.springsContourLength.size() > 0
          ? forceRelaxationNetwork.springsContourLength.mean()
          : 1.;
      forceRelaxationNetwork.assumeBoxLargeEnough = true;

      return forceRelaxationNetwork;
    }

    /**
     * @brief Run force relaxation to improve the statistics of the cross-linked
     * strands
     *
     */
    void relaxCrosslinks()
    {
      // first, convert to a useable structure for the force relaxation
      pylimer_tools::sim::mehp::Network forceRelaxationNetwork =
        this->convertToForceRelaxationNetwork();

      // actually start force relaxation
      pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
        pylimer_tools::sim::mehp::MEHPForceRelaxation(forceRelaxationNetwork);
      forceRelaxer.configAssumeBoxLargeEnough(true);

      while (forceRelaxer.suggestsRerun()) {
        forceRelaxer.runForceRelaxation();
      }

      // copy results
      forceRelaxationNetwork = forceRelaxer.getNetwork();
      size_t nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
      RUNTIME_EXP_IFN(forceRelaxationNetwork.nrOfNodes == nCrosslinks,
                      "Expected force relaxation to preserve cross-links.");
      for (size_t i = 0; i < forceRelaxationNetwork.nrOfNodes; ++i) {
        this->simplifiedUniverse.xlinkX[i] =
          forceRelaxationNetwork.coordinates(3 * i + 0);
        this->simplifiedUniverse.xlinkY[i] =
          forceRelaxationNetwork.coordinates(3 * i + 1);
        this->simplifiedUniverse.xlinkZ[i] =
          forceRelaxationNetwork.coordinates(3 * i + 2);
      }
    }

  private:
    double beadDistance;
    double meanSquaredBeadDistance;
    double primaryLoopProbability = 1.0;
    double secondaryLoopProbability = 1.0;
    size_t nMcSteps = 2000;
    std::mt19937 rng;
    std::uniform_real_distribution<double> distX;
    std::uniform_real_distribution<double> distY;
    std::uniform_real_distribution<double> distZ;

    CrosslinkerUniverse simplifiedUniverse;
    long int originalNrOfAvailableCrosslinkSites = 0;
    long int nrOfAvailableCrosslinkSites = 0;
    std::vector<int> remainingCrossLinkerFunctionality;
    pylimer_tools::entities::Box box;

    /**
     * @brief Get the Current random seed for reproducibility
     *
     * @return std::string
     */
    std::string getCurrentSeed()
    {
      std::ostringstream oss;
      oss << this->rng;
      return oss.str();
    }

    /**
     * @brief Do a random walk of certain length to add a chain
     *
     * @param from the starting Atom of the chain
     * @param chainLen the number of additional atoms to add to the chain
     * @param atomType the atom type of the atoms in the chain
     */
    Eigen::VectorXd sampleFreeChainCoordinates(int chainLen)
    {
      Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
        chainLen, this->beadDistance, this->meanSquaredBeadDistance, this->rng);

      if (this->nMcSteps > 0) {
        pylimer_tools::sim::equilibrateChainWithMC(
          positions,
          this->meanSquaredBeadDistance,
          this->rng,
          true,
          false,
          this->nMcSteps);
      }

      Eigen::Vector3d from = Eigen::Vector3d(
        this->distX(this->rng), this->distY(this->rng), this->distZ(this->rng));

      positions += from.replicate(chainLen, 1);
      return positions;
    }

    /**
     * @brief Do a random walk of certain length starting somewhere to add a
     * chain
     *
     * @param idxFrom the starting cross-link of the dangling chain
     * @param chainLen the number of additional atoms to add to the chain
     */
    Eigen::VectorXd sampleDanglingChainCoordinates(size_t idxFrom, int chainLen)
    {
      Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
        chainLen, this->beadDistance, this->meanSquaredBeadDistance, this->rng);

      if (this->nMcSteps > 0) {
        pylimer_tools::sim::equilibrateChainWithMC(
          positions,
          this->meanSquaredBeadDistance,
          this->rng,
          true,
          false,
          this->nMcSteps);
      }

      Eigen::Vector3d from =
        Eigen::Vector3d(this->simplifiedUniverse.xlinkX[idxFrom],
                        this->simplifiedUniverse.xlinkY[idxFrom],
                        this->simplifiedUniverse.xlinkZ[idxFrom]);

      positions += from.replicate(chainLen, 1);
      return positions;
    }

    /**
     * @brief Do a random walk of certain length to add a chain from one to
     * another atom
     *
     * @param from the atom to start the random walk from
     * @param to the atom to end the random walk at
     * @param chainLen the number of atoms to add in between from and to
     */
    Eigen::VectorXd sampleStrandCoordinates(size_t from,
                                            size_t to,
                                            int chainLen)
    {
      // determine the positions
      Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChainFromTo(
        this->box,
        Eigen::Vector3d(this->simplifiedUniverse.xlinkX[from],
                        this->simplifiedUniverse.xlinkY[from],
                        this->simplifiedUniverse.xlinkZ[from]),
        Eigen::Vector3d(this->simplifiedUniverse.xlinkX[to],
                        this->simplifiedUniverse.xlinkY[to],
                        this->simplifiedUniverse.xlinkZ[to]),
        chainLen,
        this->beadDistance,
        this->meanSquaredBeadDistance,
        this->rng,
        true);

      if (this->nMcSteps > 0) {
        pylimer_tools::sim::equilibrateChainWithMC(
          positions,
          this->meanSquaredBeadDistance,
          this->rng,
          true,
          true,
          this->nMcSteps);
      }

      assert(positions.size() == (chainLen + 2) * 3);

      // omit first and last, which were required
      // for the MC stepping
      return positions.segment(3, positions.size() - 6);
    }

    /**
     * @brief Tell a strand how it is connected to a cross-link.
     * Automatically decides whether it is `from` or `to`
     *
     * @param strandIdx
     * @param crosslinkIdx
     * @param ignoreInexistent whether to throw an error if the crosslink does
     * not exist yet, or not
     */
    void linkStrandToCrosslink(size_t strandIdx,
                               size_t crosslinkIdx,
                               bool ignoreInexistent = false)
    {
      INVALIDARG_EXP_IFN(strandIdx < this->simplifiedUniverse.strandFrom.size(),
                         "The strand index is out of bounds.");
      INVALIDARG_EXP_IFN(this->simplifiedUniverse.strandFrom[strandIdx] < 0 ||
                           this->simplifiedUniverse.strandTo[strandIdx] < 0,
                         "Require one end to be free to be connected to");
      if (!ignoreInexistent) {
        INVALIDARG_EXP_IFN(
          crosslinkIdx < this->simplifiedUniverse.strandFrom.size(), "");
      }

      if (this->simplifiedUniverse.strandFrom[strandIdx] < 0) {
        this->simplifiedUniverse.strandFrom[strandIdx] = crosslinkIdx;
      } else {
        this->simplifiedUniverse.strandTo[strandIdx] = crosslinkIdx;
      }
    };

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
                                      Eigen::VectorXd coordinates)
    {
      INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                         "Coordinates must have a size multiple of 3");
      INVALIDARG_EXP_IFN(coordinates.size() / 3 == nrOfAtomsToAdd,
                         "Coordinates must match the promised number of atoms");
      size_t currentNrOfJunctions = this->simplifiedUniverse.xlinkX.size();
      this->simplifiedUniverse.xlinkTypes.reserve(currentNrOfJunctions +
                                                  nrOfAtomsToAdd);
      this->simplifiedUniverse.xlinkX.reserve(currentNrOfJunctions +
                                              nrOfAtomsToAdd);
      this->simplifiedUniverse.xlinkY.reserve(currentNrOfJunctions +
                                              nrOfAtomsToAdd);
      this->simplifiedUniverse.xlinkZ.reserve(currentNrOfJunctions +
                                              nrOfAtomsToAdd);

      std::vector<size_t> indicesAdded;
      indicesAdded.reserve(nrOfAtomsToAdd);

      for (size_t i = 0; i < nrOfAtomsToAdd; ++i) {
        this->simplifiedUniverse.xlinkTypes.push_back(atomType);
        indicesAdded.push_back(currentNrOfJunctions + i);

        this->simplifiedUniverse.xlinkX.push_back(coordinates(3 * i));
        this->simplifiedUniverse.xlinkY.push_back(coordinates(3 * i + 1));
        this->simplifiedUniverse.xlinkZ.push_back(coordinates(3 * i + 2));

        this->simplifiedUniverse.strandsOfXlink.push_back({});
      }

      return indicesAdded;
    }

    /**
     * @brief Add atoms (incl. random positions, id etc.) to the universe
     *
     * @param nrOfAtomsToAdd the nr. of atoms to add to the universe
     * @param atomType the type of the atoms to add
     * @return std::vector<size_t> the ids of the inserted atoms
     */
    std::vector<size_t> addXlinkAtoms(int nrOfAtomsToAdd,
                                      int atomType,
                                      bool whiteNoise = true)
    {
      Eigen::VectorXd randomPos =
        this->generateRandomPositions(nrOfAtomsToAdd, whiteNoise);
      return this->addXlinkAtoms(nrOfAtomsToAdd, atomType, randomPos);
    }

    /**
     * @brief Generate positions randomly
     *
     * @param nSamples the number of positions to generate
     * @return Positions
     */
    Eigen::VectorXd generateRandomPositions(int nSamples,
                                            bool whiteNoise = true)
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
    Eigen::VectorXd generateRandomWhitePositions(int nSamples)
    {
      Eigen::VectorXd coordinates = Eigen::VectorXd(3 * nSamples);

      for (size_t i = 0; i < nSamples; ++i) {
        coordinates(3 * i) = this->distX(this->rng);
        coordinates(3 * i + 1) = this->distY(this->rng);
        coordinates(3 * i + 2) = this->distZ(this->rng);
      }

      return coordinates;
    }

    /**
     * @brief Generate positions randomly according to blue noise type sampling
     *
     * @param nSamples the number of positions to generate
     * @return Positions
     */
    Eigen::VectorXd generateRandomBluePositions(int nSamples)
    {
      Eigen::VectorXd coordinates = Eigen::VectorXd(3 * nSamples);

      // blue noise
      // inspiration:
      // https://github.com/Atrix256/RandomCode/blob/master/Mitchell/Source.cpp
      for (size_t i = 0; i < nSamples; ++i) {
        size_t numCandidates =
          std::min((size_t)(i * 1 + 1),
                   (size_t)500); // decrease the multiplier to speed things up
        double bestDistance = 0.0;
        double bestCandidateX = 0.0;
        double bestCandidateY = 0.0;
        double bestCandidateZ = 0.0;

        for (size_t j = 0; j < numCandidates; ++j) {
          double x = this->distX(this->rng);
          double y = this->distY(this->rng);
          double z = this->distZ(this->rng);

          double minDistance = std::numeric_limits<double>::max();
          for (size_t k = 0; k < nSamples; ++k) {
            double dist = this->getDistance(x,
                                            y,
                                            z,
                                            coordinates(3 * k),
                                            coordinates(3 * k + 1),
                                            coordinates(3 * k + 2));
            if (dist < minDistance) {
              minDistance = dist;
            }
          }

          if (minDistance > bestDistance) {
            bestDistance = minDistance;

            bestCandidateX = x;
            bestCandidateY = y;
            bestCandidateZ = z;
          }
        }

        coordinates(3 * i) = bestCandidateX;
        coordinates(3 * i + 1) = bestCandidateY;
        coordinates(3 * i + 2) = bestCandidateZ;
      }

      return coordinates;
    }

    /**
     * @brief find an Atom in a collection with an objective distance
     *
     * TODO: check probabilities to not only select the very best.
     *
     * @param from the Atom to start the distance from
     * @param desiredR02 the distance to target if possible
     * @param maxDistance the maximum distance to accept random matches from
     * @return int the index in possiblePartners that matches best
     */
    size_t findAppropriateLink(size_t from,
                               const double desiredR02,
                               const double maxDistance)
    {
      assert(this->simplifiedUniverse.xlinkTypes.size() ==
             this->remainingCrossLinkerFunctionality.size());

      std::vector<size_t> suitableMatches;
      std::vector<double> matchWeights;
      double sumOfWeights = 0.0;
      const double normalisationFactorInExponential = -3. / (2. * desiredR02);
      size_t nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
      // TODO: use a neighbour list instead, maybe?
      for (int i = 0; i < nCrosslinks; ++i) {
        if (this->remainingCrossLinkerFunctionality[i] < 1) {
          continue;
        }
        size_t partner = i;
        Eigen::Vector3d dist = this->getVectorBetween(from, i);
        if (dist.norm() < maxDistance || maxDistance < 0.) {
          suitableMatches.push_back(partner);
          double thisWeight =
            static_cast<double>(this->remainingCrossLinkerFunctionality[i]) *
            std::exp(dist.squaredNorm() * normalisationFactorInExponential);
          if (partner == from) {
            thisWeight *= this->primaryLoopProbability;
          }
          if (this->secondaryLoopProbability != 1.) {
            // check whether this cross-link would lead to a secondary loop
            // => apply weight for every other already existing back-link
            for (size_t partnersSubStrand :
                 this->simplifiedUniverse.strandsOfXlink[partner]) {
              assert(this->simplifiedUniverse.strandFrom[partnersSubStrand] ==
                       partner ||
                     this->simplifiedUniverse.strandTo[partnersSubStrand] ==
                       partner);
              if (this->simplifiedUniverse.strandFrom[partnersSubStrand] ==
                    from ||
                  this->simplifiedUniverse.strandTo[partnersSubStrand] ==
                    from) {
                thisWeight *= this->secondaryLoopProbability;
              }
            }
          }
          matchWeights.push_back(thisWeight);
          sumOfWeights += thisWeight;
        }
      }

      if (suitableMatches.size() == 0) {
        throw std::runtime_error("No suitable partner found.");
      }

      std::discrete_distribution<int> weightDist(matchWeights.begin(),
                                                 matchWeights.end());
      return suitableMatches[weightDist(this->rng)];
    }

    ///// Utility functions
    /**
     * @brief compute the distance between two cross-links, given by their
     * indices
     *
     * @param i
     * @param j
     * @return double
     */
    double distanceBetween(size_t i, size_t j)
    {
      return this->getDistance(this->simplifiedUniverse.xlinkX[i],
                               this->simplifiedUniverse.xlinkY[i],
                               this->simplifiedUniverse.xlinkZ[i],
                               this->simplifiedUniverse.xlinkX[j],
                               this->simplifiedUniverse.xlinkY[j],
                               this->simplifiedUniverse.xlinkZ[j]);
    }

    Eigen::Vector3d getVectorBetween(size_t i, size_t j)
    {
      Eigen::Vector3d diff;
      diff << (this->simplifiedUniverse.xlinkX[j] -
               this->simplifiedUniverse.xlinkX[i]),
        (this->simplifiedUniverse.xlinkY[j] -
         this->simplifiedUniverse.xlinkY[i]),
        (this->simplifiedUniverse.xlinkZ[j] -
         this->simplifiedUniverse.xlinkZ[i]);
      this->box.handlePBC(diff);
      return diff;
    }

    double getDistance(double x1,
                       double y1,
                       double z1,
                       double x2,
                       double y2,
                       double z2)
    {
      Eigen::Vector3d diff;
      diff << x2 - x1, y2 - y1, z2 - z1;
      this->box.handlePBC(diff);
      return diff.norm();
    }

    void validateInternalState() const
    {
      RUNTIME_EXP_IFN(
        all_equal<size_t>(
          5,
          this->simplifiedUniverse.strandFrom.size(),
          this->simplifiedUniverse.strandTo.size(),
          this->simplifiedUniverse.beadsInStrand.size(),
          this->simplifiedUniverse.beadDistanceInStrand.size(),
          this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.size()),
        "Inconsistent sizes in simplified universe.");
      RUNTIME_EXP_IFN(
        all_equal<size_t>(5,
                          this->simplifiedUniverse.xlinkTypes.size(),
                          this->simplifiedUniverse.xlinkX.size(),
                          this->simplifiedUniverse.xlinkY.size(),
                          this->simplifiedUniverse.xlinkZ.size(),
                          this->remainingCrossLinkerFunctionality.size()),
        "Inconsistent sizes in simplified universe.");

      long int nCrosslinks = this->simplifiedUniverse.xlinkX.size();
      const long int nrOfAvailableSites =
        std::reduce(this->remainingCrossLinkerFunctionality.begin(),
                    this->remainingCrossLinkerFunctionality.end(),
                    0);

      RUNTIME_EXP_IFN(this->nrOfAvailableCrosslinkSites == nrOfAvailableSites,
                      "Inconsistent nr of cross-link sites.");

      for (size_t xlinkIdx = 0; xlinkIdx < nCrosslinks; ++xlinkIdx) {
        for (long int subStrandIdx :
             this->simplifiedUniverse.strandsOfXlink[xlinkIdx]) {
          RUNTIME_EXP_IFN(
            this->simplifiedUniverse.strandFrom[subStrandIdx] == xlinkIdx ||
              this->simplifiedUniverse.strandTo[subStrandIdx] == xlinkIdx,
            "Inconsistent links int list of cross-links <> strands.");
        }
      }
      for (size_t strandIdx;
           strandIdx < this->simplifiedUniverse.strandFrom.size();
           ++strandIdx) {
        RUNTIME_EXP_IFN(this->simplifiedUniverse.beadsInStrand[strandIdx] >= 0,
                        "Expected a positive number of beads per strand.");
      }
    }
  };
} // namespace utils
} // namespace pylimer_tools

#endif
