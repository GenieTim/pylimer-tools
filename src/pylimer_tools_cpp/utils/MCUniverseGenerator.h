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
#include <string>
#include <vector>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

#include <random>

namespace pylimer_tools {
namespace utils {

  struct CrosslinkerUniverse
  {
    std::vector<int> xlinkTypes;
    std::vector<double> xlinkX;
    std::vector<double> xlinkY;
    std::vector<double> xlinkZ;
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
      this->primaryLoopProbability = newPrimaryLoopProbability;
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

      std::vector<size_t> newCrosslinkerIdxs = this->addAtomsWithType(
        nrOfCrosslinkers, crossLinkerAtomType, whiteNoise);
      this->crossLinkerIdxs.insert(this->crossLinkerIdxs.end(),
                                   newCrosslinkerIdxs.begin(),
                                   newCrosslinkerIdxs.end());
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
        this->simplifiedUniverse.strandFrom.push_back(-1);
        this->simplifiedUniverse.strandTo.push_back(-1);
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
    void addAndLinkMonofunctionalStrands(int nrOfStrands,
                                         std::vector<int> beadsPerChains,
                                         int strandAtomType = 1)
    {
      INVALIDARG_EXP_IFN(beadsPerChains.size() == nrOfStrands,
                         "Nr of strands (" + std::to_string(nrOfStrands) +
                           ") must be equal to the number "
                           "of chainLengths (" +
                           std::to_string(beadsPerChains.size()) +
                           ") provided.");

      long int nCrosslinks = this->crossLinkerIdxs.size();
      long int nrOfAvailableSites =
        std::reduce(this->remainingCrossLinkerFunctionality.begin(),
                    this->remainingCrossLinkerFunctionality.end(),
                    0);

      RUNTIME_EXP_IFN(nrOfStrands <= nrOfAvailableSites,
                      "Not enough cross-link sites available to link all these "
                      "monofunctional strands.");

      std::vector<size_t> availableCrosslinkSites;
      availableCrosslinkSites.reserve(nrOfAvailableSites);
      for (size_t i = 0; i < nCrosslinks; ++i) {
        for (long int s = 0; s < this->remainingCrossLinkerFunctionality[i];
             ++s) {
          availableCrosslinkSites.push_back(this->crossLinkerIdxs[i]);
        }
      }
      std::shuffle(availableCrosslinkSites.begin(),
                   availableCrosslinkSites.end(),
                   this->rng);

      for (size_t i = 0; i < nrOfStrands; ++i) {
        this->simplifiedUniverse.strandBeadType.push_back(strandAtomType);
        this->simplifiedUniverse.beadsInStrand.push_back(beadsPerChains[i]);
        this->simplifiedUniverse.beadDistanceInStrand.push_back(
          this->beadDistance);
        this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.push_back(
          this->meanSquaredBeadDistance);
        this->simplifiedUniverse.strandFrom.push_back(
          availableCrosslinkSites[i]);
        this->remainingCrossLinkerFunctionality[availableCrosslinkSites[i]] -=
          1;
        this->nrOfAvailableCrosslinkSites -= 1;
        this->simplifiedUniverse.strandTo.push_back(-1);
      }
    }

    /**
     * @brief Add strands, link them to the cross-links, stop when the callback
     * says so
     *
     * @param nrOfStrands the number of strands to add
     * @param beadsPerChains the number of beads per strand `N`
     * @param stopLinking a callback to indicated whether to stop linking
     * @param strandAtomType the type of the atoms in the strands
     * @param cInfinity `C_\infty` for `<R_ee^2>_0` from `N` and `<b^2>`
     */
    void addAndLinkStrandsCallback(
      int nrOfStrands,
      std::vector<int> beadsPerChains,
      std::function<bool(const MCUniverseGenerator&)> stopLinking,
      int strandAtomType = 1,
      double cInfinity = 1.)
    {
      INVALIDARG_EXP_IFN(beadsPerChains.size() == nrOfStrands,
                         "Nr of strands (" + std::to_string(nrOfStrands) +
                           ") must be equal to the number "
                           "of chainLengths (" +
                           std::to_string(beadsPerChains.size()) +
                           ") provided.");

      RUNTIME_EXP_IFN(
        this->crossLinkerIdxs.size() ==
          this->remainingCrossLinkerFunctionality.size(),
        "Invalid internal state, the remaining cross-linker functionalities "
        "are inconsistent with the indices of the cross-links.");
      this->validateInternalState();

      // prepare sampling of partners
      long int nCrosslinks = this->crossLinkerIdxs.size();
      int nrOfStrandsAdded = 0;

      double currentCrosslinkerConversion =
        1. - (static_cast<double>(this->nrOfAvailableCrosslinkSites) /
              static_cast<double>(this->originalNrOfAvailableCrosslinkSites));
      double conversionPerBond =
        1. / (static_cast<double>(this->originalNrOfAvailableCrosslinkSites));

      size_t nStrandsBefore = this->simplifiedUniverse.strandFrom.size();

      std::vector<size_t> availableStrandEnds;
      availableStrandEnds.reserve(2 * nrOfStrands);
      for (size_t i = 0; i < nrOfStrands; ++i) {
        // each strand is available with two ends
        availableStrandEnds.push_back(i + nStrandsBefore);
        availableStrandEnds.push_back(i + nStrandsBefore);
      }
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

      for (size_t strandIdx = 0; strandIdx < nrOfStrands; ++strandIdx) {
        this->simplifiedUniverse.strandBeadType.push_back(strandAtomType);
        this->simplifiedUniverse.beadsInStrand.push_back(
          beadsPerChains[strandIdx]);
        this->simplifiedUniverse.beadDistanceInStrand.push_back(
          this->beadDistance);
        this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.push_back(
          this->meanSquaredBeadDistance);
        this->simplifiedUniverse.strandFrom.push_back(-1);
        this->simplifiedUniverse.strandTo.push_back(-1);
      }

      const double timesNForR02 = this->meanSquaredBeadDistance * cInfinity;

      // link one strand at a time until we reach the target conversion
      for (int sampleIdx = 0; sampleIdx < 2 * nrOfStrands; ++sampleIdx) {
        if (stopLinking(*this)) {
          break;
        }

        size_t strandIdx = availableStrandEnds[sampleIdx];

        if (this->simplifiedUniverse.strandFrom[strandIdx] != -1) {
          // we don't have free cross-link choice
          RUNTIME_EXP_IFN(this->simplifiedUniverse.strandTo[strandIdx] == -1,
                          "Expected second strand end to be free");

          size_t partnerCrosslinker = this->findAppropriateLink(
            this->simplifiedUniverse.strandFrom[strandIdx],
            static_cast<double>(beadsPerChains[strandIdx] + 1) * timesNForR02,
            -1. // beadsPerChains[strandIdx] * this->beadDistance
          );

          this->simplifiedUniverse.strandTo[strandIdx] = partnerCrosslinker;
          this->remainingCrossLinkerFunctionality[partnerCrosslinker] -= 1;
          this->nrOfAvailableCrosslinkSites -= 1;
        } else {
          // otherwise, randomly choose a free cross-link
          long int crosslinkIdxIdx = availableCrosslinkSites.back();
          availableCrosslinkSites.pop_back();
          while (this->remainingCrossLinkerFunctionality[crosslinkIdxIdx] < 1 &&
                 availableCrosslinkSites.size() > 0) {
            // find the next "available" cross-linker
            crosslinkIdxIdx = availableCrosslinkSites.back();
            availableCrosslinkSites.pop_back();
          }

          if (availableCrosslinkSites.size() == 0) {
            std::cerr << "No more cross-link sites available." << std::endl;
            break;
          }

          this->simplifiedUniverse.strandFrom[strandIdx] =
            this->crossLinkerIdxs[crosslinkIdxIdx];
          this->remainingCrossLinkerFunctionality[crosslinkIdxIdx] -= 1;
          this->nrOfAvailableCrosslinkSites -= 1;
        }
      }
    }

    /**
     * @brief Add strands in between the cross-linkers, link them as appropriate
     *
     * @param nrOfStrands the nr. of Strands to add
     * @param beadsPerChains the nr. of beads per strand (excl. cross-linkers)
     * @param targetCrossLinkerConversion "p", the target conversion of the
     * cross-linkers
     * @param strandAtomType the type of the strand atoms
     */
    void addAndLinkStrandsToConversion(const int nrOfStrands,
                                       const std::vector<int> beadsPerChains,
                                       const double targetCrossLinkerConversion,
                                       const int strandAtomType = 1,
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

      const int potentialNewBonds = 2 * nrOfStrands;
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

      this->addAndLinkStrandsCallback(
        nrOfStrands,
        beadsPerChains,
        [targetNrOfAvailableCrosslinkSites](const MCUniverseGenerator& gen) {
          return gen.nrOfAvailableCrosslinkSites <=
                 targetNrOfAvailableCrosslinkSites;
        },
        strandAtomType,
        cInfinity);
    }

    void addAndLinkStrandsToConversion(int nrOfStrands,
                                       int chainLength,
                                       double crossLinkerConversion,
                                       int strandAtomType = 1,
                                       double cInfinity = 1.)
    {
      const std::vector<int> chainLengths =
        pylimer_tools::utils::initializeWithValue<int>(nrOfStrands,
                                                       chainLength);
      return this->addAndLinkStrandsToConversion(nrOfStrands,
                                                 chainLengths,
                                                 crossLinkerConversion,
                                                 strandAtomType,
                                                 cInfinity);
    };

    /**
     * @brief Add strands in between the cross-linkers, link them as appropriate
     *
     * @param nrOfStrands the nr. of Strands to add
     * @param beadsPerChains the nr. of beads per strand (excl. cross-linkers)
     * @param targetSolubleFraction "w_sol", the target soluble fraction
     * @param strandAtomType the type of the strand atoms
     */
    void addAndLinkStrandsToSolubleFraction(int nrOfStrands,
                                            std::vector<int> beadsPerChains,
                                            double targetSolubleFraction,
                                            int strandAtomType = 1,
                                            double cInfinity = 1.)
    {
      INVALIDARG_EXP_IFN(targetSolubleFraction >= 0. &&
                           targetSolubleFraction <= 1.,
                         "Soluble fraction must be between 0 and 1, got " +
                           std::to_string(targetSolubleFraction) + ".");

      this->addAndLinkStrandsCallback(
        nrOfStrands,
        beadsPerChains,
        [targetSolubleFraction](const MCUniverseGenerator& gen) {
          pylimer_tools::sim::mehp::Network frNet =
            gen.convertToForceRelaxationNetwork();

          // this network only contains non-dangling and non-free chains
          // need to compute compensation for that.
          size_t nOmittedBeads = 0;
          for (size_t i = 0; i < gen.simplifiedUniverse.strandFrom.size();
               ++i) {
            if (gen.simplifiedUniverse.strandTo[i] == -1) {
              nOmittedBeads += gen.simplifiedUniverse.beadsInStrand[i];
            }
          }

          // actually start force relaxation
          pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
            pylimer_tools::sim::mehp::MEHPForceRelaxation(frNet);
          forceRelaxer.configAssumeBoxLargeEnough(true);

          while (forceRelaxer.suggestsRerun()) {
            forceRelaxer.runForceRelaxation();
          }

          // finally, calculate the soluble fraction
          return forceRelaxer.getSolubleWeightFraction() >=
                 targetSolubleFraction;
        },
        strandAtomType,
        cInfinity);
    }

    void addAndLinkStrandsToSolubleFraction(int nrOfStrands,
                                            int chainLength,
                                            double targetSolubleFraction,
                                            int strandAtomType = 1,
                                            double cInfinity = 1.)
    {
      std::vector<int> chainLengths =
        pylimer_tools::utils::initializeWithValue<int>(nrOfStrands,
                                                       chainLength);
      return this->addAndLinkStrandsToSolubleFraction(nrOfStrands,
                                                      chainLengths,
                                                      targetSolubleFraction,
                                                      strandAtomType,
                                                      cInfinity);
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

      forceRelaxationNetwork.nrOfNodes = this->crossLinkerIdxs.size();
      forceRelaxationNetwork.oldAtomIds =
        Eigen::ArrayXi::Zero(forceRelaxationNetwork.nrOfNodes);
      forceRelaxationNetwork.coordinates =
        Eigen::VectorXd(forceRelaxationNetwork.nrOfNodes * 3);
      size_t i = 0;
      for (size_t crosslinkIdx : this->crossLinkerIdxs) {
        forceRelaxationNetwork.coordinates(3 * i + 0) =
          this->simplifiedUniverse.xlinkX[crosslinkIdx];
        forceRelaxationNetwork.coordinates(3 * i + 1) =
          this->simplifiedUniverse.xlinkY[crosslinkIdx];
        forceRelaxationNetwork.coordinates(3 * i + 2) =
          this->simplifiedUniverse.xlinkZ[crosslinkIdx];

        forceRelaxationNetwork.springIndicesOfLinks.push_back({});

        i += 1;
      }

      size_t nSpringEstimate = this->simplifiedUniverse.strandFrom.size();
      forceRelaxationNetwork.springIndexA =
        Eigen::ArrayXi::Zero(nSpringEstimate);
      forceRelaxationNetwork.springIndexB =
        Eigen::ArrayXi::Zero(nSpringEstimate);
      forceRelaxationNetwork.springsContourLength =
        Eigen::VectorXd::Zero(nSpringEstimate);
      size_t nActualSprings = 0;
      // we omit all dangling and free strands
      for (size_t springIdx = 0;
           springIdx < this->simplifiedUniverse.strandFrom.size();
           ++springIdx) {
        if (this->simplifiedUniverse.strandTo[springIdx] != -1) {
          long int from = this->simplifiedUniverse.strandFrom[springIdx];
          long int to = this->simplifiedUniverse.strandTo[springIdx];
          forceRelaxationNetwork.springIndexA(springIdx) = from;
          forceRelaxationNetwork.springIndexB(springIdx) = to;
          forceRelaxationNetwork.springIndicesOfLinks[from].push_back(
            springIdx);
          forceRelaxationNetwork.springIndicesOfLinks[to].push_back(springIdx);
          forceRelaxationNetwork.springsContourLength(springIdx) =
            this->simplifiedUniverse.beadsInStrand[springIdx] + 1;
          nActualSprings += 1;
        }
      }
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
        forceRelaxationNetwork.springsContourLength.mean();
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
      RUNTIME_EXP_IFN(forceRelaxationNetwork.nrOfNodes ==
                        this->crossLinkerIdxs.size(),
                      "Expected force relaxation to preserve cross-links.");
      for (size_t i = 0; i < forceRelaxationNetwork.nrOfNodes; ++i) {
        this->simplifiedUniverse.xlinkX[this->crossLinkerIdxs[i]] =
          forceRelaxationNetwork.coordinates(3 * i + 0);
        this->simplifiedUniverse.xlinkY[this->crossLinkerIdxs[i]] =
          forceRelaxationNetwork.coordinates(3 * i + 1);
        this->simplifiedUniverse.xlinkZ[this->crossLinkerIdxs[i]] =
          forceRelaxationNetwork.coordinates(3 * i + 2);
      }
    }

  private:
    double beadDistance;
    double meanSquaredBeadDistance;
    double primaryLoopProbability = 1.0;
    size_t nMcSteps = 2000;
    std::mt19937 rng;
    std::uniform_real_distribution<double> distX;
    std::uniform_real_distribution<double> distY;
    std::uniform_real_distribution<double> distZ;

    CrosslinkerUniverse simplifiedUniverse;
    std::vector<size_t> crossLinkerIdxs;
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

      pylimer_tools::sim::equilibrateChainWithMC(positions,
                                                 this->meanSquaredBeadDistance,
                                                 this->rng,
                                                 true,
                                                 false,
                                                 this->nMcSteps);

      Eigen::Vector3d from = Eigen::Vector3d(
        this->distX(this->rng), this->distY(this->rng), this->distZ(this->rng));

      positions += from.replicate(chainLen, 1);
      return positions;
    }

    /**
     * @brief Do a random walk of certain length starting somewhere to add a
     * chain
     *
     * @param from the starting Atom of the chain
     * @param chainLen the number of additional atoms to add to the chain
     * @param atomType the atom type of the atoms in the chain
     */
    Eigen::VectorXd sampleDanglingChainCoordinates(size_t idxFrom, int chainLen)
    {
      Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
        chainLen, this->beadDistance, this->meanSquaredBeadDistance, this->rng);

      pylimer_tools::sim::equilibrateChainWithMC(positions,
                                                 this->meanSquaredBeadDistance,
                                                 this->rng,
                                                 true,
                                                 false,
                                                 this->nMcSteps);

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
      Eigen::VectorXd
        positions = // pylimer_tools::utils::doRandomWalkChainFromTo(
        pylimer_tools::utils::doRandomWalkChainFromToMC(
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
          this->nMcSteps);

      assert(positions.size() == chainLen * 3);

      return positions;
    }

    /**
     * @brief Add atoms (incl. type, id etc.) with given positions to the
     * universe
     *
     * @param nrOfAtomsToAdd the nr. of atoms to add to the universe
     * @param atomType the type of the atoms to add
     * @return std::vector<size_t> the ids of the inserted atoms
     */
    std::vector<size_t> addAtomsWithType(int nrOfAtomsToAdd,
                                         int atomType,
                                         Eigen::VectorXd coordinates)
    {
      INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                         "Coordinates must have a size multiple of 3");
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
    std::vector<size_t> addAtomsWithType(int nrOfAtomsToAdd,
                                         int atomType,
                                         bool whiteNoise = true)
    {
      Eigen::VectorXd randomPos =
        this->generateRandomPositions(nrOfAtomsToAdd, whiteNoise);
      return this->addAtomsWithType(nrOfAtomsToAdd, atomType, randomPos);
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
      assert(this->crossLinkerIdxs.size() ==
             this->remainingCrossLinkerFunctionality.size());

      std::vector<size_t> suitableMatches;
      std::vector<double> matchWeights;
      double sumOfWeights = 0.0;
      const double normalisationFactorInExponential = -3. / (2. * desiredR02);
      // TODO: use a neighbour list instead, maybe?
      for (int i = 0; i < this->crossLinkerIdxs.size(); ++i) {
        if (this->remainingCrossLinkerFunctionality[i] < 1) {
          continue;
        }
        size_t partner = this->crossLinkerIdxs[i];
        Eigen::Vector3d dist =
          this->getVectorBetween(from, this->crossLinkerIdxs[i]);
        if (dist.norm() < maxDistance || maxDistance < 0.) {
          suitableMatches.push_back(partner);
          double thisWeight =
            static_cast<double>(this->remainingCrossLinkerFunctionality[i]) *
            std::exp(dist.squaredNorm() * normalisationFactorInExponential);
          if (partner == from) {
            thisWeight *= this->primaryLoopProbability;
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

      long int nCrosslinks = this->crossLinkerIdxs.size();
      const long int nrOfAvailableSites =
        std::reduce(this->remainingCrossLinkerFunctionality.begin(),
                    this->remainingCrossLinkerFunctionality.end(),
                    0);

      RUNTIME_EXP_IFN(this->nrOfAvailableCrosslinkSites == nrOfAvailableSites,
                      "Inconsistent nr of cross-link sites.");
    }
  };
} // namespace utils
} // namespace pylimer_tools

#endif
