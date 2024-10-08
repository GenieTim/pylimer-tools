#ifndef MC_UNIVERSE_GENERATOR_H
#define MC_UNIVERSE_GENERATOR_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "RandomWalker.h"
#include "StringUtils.h"
#include "VectorUtils.h"
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

  struct Positions
  {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
  };

  struct SimplifiedUniverse
  {
    std::vector<long int> ids;
    std::vector<int> types;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<long int> bondsFrom;
    std::vector<long int> bondsTo;
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
      this->distSelect = std::uniform_real_distribution<double>(0.0, 1.0);
      this->box = pylimer_tools::entities::Box(Lx, Ly, Lz);
    }

    void setSeed(unsigned int seed) { this->rng.seed(seed); };

    void setBeadDistance(double newBeadDistance)
    {
      this->beadDistance = newBeadDistance;
    };

    pylimer_tools::entities::Universe getUniverse()
    {
      pylimer_tools::entities::Universe universe =
        pylimer_tools::entities::Universe(this->box);
      int nrOfAtoms = this->simplifiedUniverse.ids.size();
      std::vector<int> zeros = initializeWithValue(nrOfAtoms, 0);
      universe.addAtoms(this->simplifiedUniverse.ids,
                        this->simplifiedUniverse.types,
                        this->simplifiedUniverse.x,
                        this->simplifiedUniverse.y,
                        this->simplifiedUniverse.z,
                        zeros,
                        zeros,
                        zeros);
      universe.addBonds(this->simplifiedUniverse.bondsFrom,
                        this->simplifiedUniverse.bondsTo);
      return universe;
    };

    void addCrosslinkers(int nrOfCrosslinkers,
                         int crosslinkerFunctionality = 4,
                         int crossLinkerAtomType = 2)
    {
      int nCrosslinkerBefore = this->remainingCrossLinkerFunctionality.size();

      std::vector<size_t> newCrosslinkerIdxs =
        this->addAtomsWithType(nrOfCrosslinkers, crossLinkerAtomType);
      this->crossLinkerIdxs.insert(this->crossLinkerIdxs.end(),
                                   newCrosslinkerIdxs.begin(),
                                   newCrosslinkerIdxs.end());
      this->remainingCrossLinkerFunctionality.reserve(nCrosslinkerBefore +
                                                      nrOfCrosslinkers);
      for (size_t i = 0; i < nrOfCrosslinkers; ++i) {
        this->remainingCrossLinkerFunctionality.push_back(
          crosslinkerFunctionality);
      }
    };

    /**
     * @brief Randomly distribute
     *
     * @param nrOfSolventChains
     * @param chainLength
     * @param solventAtomType
     */
    void addSolventChains(int nrOfSolventChains,
                          int chainLength,
                          int solventAtomType = 3)
    {
      // std::cout << "Adding " << nrOfSolventChains << " atoms for solvent
      // chains."
      //           << std::endl;
      std::vector<size_t> startAtoms =
        this->addAtomsWithType(nrOfSolventChains, solventAtomType);

      for (size_t atomIdx : startAtoms) {
        // std::cout << "Adding solvent chain with length " << chainLength
        //           << std::endl;
        this->addRandomWalkChainFrom(atomIdx, chainLength - 1, solventAtomType);
      }
    };

    /**
     * @brief Add strands in between the cross-linkers, link them as appropriate
     *
     * @param nrOfStrands the nr. of Strands to add
     * @param beadsPerChains the nr. of beads per strand (excl. cross-linkers)
     * @param targetCrossLinkerConversion "p", the target conversion of the
     * cross-linkers
     * @param strandAtomType the type of the strand atoms
     */
    void addAndLinkStrands(int nrOfStrands,
                           std::vector<int> beadsPerChains,
                           double targetCrossLinkerConversion,
                           int strandAtomType = 1)
    {
      if (beadsPerChains.size() != nrOfStrands) {
        throw std::invalid_argument(
          "Nr of strands (" + std::to_string(nrOfStrands) +
          ") must be equal to the number "
          "of chainLengths (" +
          std::to_string(beadsPerChains.size()) + ") provided.");
      }

      INVALIDARG_EXP_IFN(
        targetCrossLinkerConversion >= 0.0 &&
          targetCrossLinkerConversion <= 1.0,
        "Cross-linker conversion must be between 0 and 1, got " +
          std::to_string(targetCrossLinkerConversion) + ".");

      assert(this->crossLinkerIdxs.size() ==
             this->remainingCrossLinkerFunctionality.size());

      // eager reserve of vectors
      // could be more elaborate, but with a certain probability,
      // this is better
      this->simplifiedUniverse.ids.reserve(this->simplifiedUniverse.ids.size() +
                                           nrOfStrands * beadsPerChains[0]);
      this->simplifiedUniverse.types.reserve(
        this->simplifiedUniverse.types.size() +
        nrOfStrands * beadsPerChains[0]);

      this->simplifiedUniverse.x.reserve(this->simplifiedUniverse.x.size() +
                                         nrOfStrands * beadsPerChains[0]);
      this->simplifiedUniverse.y.reserve(this->simplifiedUniverse.y.size() +
                                         nrOfStrands * beadsPerChains[0]);
      this->simplifiedUniverse.z.reserve(this->simplifiedUniverse.z.size() +
                                         nrOfStrands * beadsPerChains[0]);
      this->simplifiedUniverse.bondsFrom.reserve(
        this->simplifiedUniverse.bondsFrom.size() +
        nrOfStrands * beadsPerChains[0]);
      this->simplifiedUniverse.bondsTo.reserve(
        this->simplifiedUniverse.bondsTo.size() +
        nrOfStrands * beadsPerChains[0]);

      //
      long int nCrosslinks = this->crossLinkerIdxs.size();
      std::uniform_int_distribution<size_t> crosslinkerIdxIdxDist =
        std::uniform_int_distribution<size_t>(0,
                                              nCrosslinks - 1);
      int nrOfStrandsAdded = 0;
      long int nrOfAvailableSites =
        std::reduce(this->remainingCrossLinkerFunctionality.begin(),
                    this->remainingCrossLinkerFunctionality.end(),
                    0);

      double conversionPerBond = (1.0 - this->currentCrosslinkerConversion) /
                                 (static_cast<double>(nrOfAvailableSites));

      int potentialNewBonds = 2 * nrOfStrands;
      if (this->currentCrosslinkerConversion +
            potentialNewBonds * conversionPerBond <
          targetCrossLinkerConversion) {
        throw std::invalid_argument(
          "A cross-linker conversion of " +
          std::to_string(targetCrossLinkerConversion) +
          " is not reachable with this nr of strands.");
      }

      std::vector<size_t> availableStrandEnds;
      availableStrandEnds.reserve(2 * nrOfStrands);
      for (size_t i = 0; i < nrOfStrands; ++i) {
        availableStrandEnds.push_back(i);
        availableStrandEnds.push_back(i);
      }
      std::shuffle(
        availableStrandEnds.begin(), availableStrandEnds.end(), this->rng);

      std::vector<long int> strandEnd1 =
        pylimer_tools::utils::initializeWithValue<long int>(nrOfStrands, -1);
      std::vector<long int> strandEnd2 =
        pylimer_tools::utils::initializeWithValue<long int>(nrOfStrands, -1);

      // link one strand at a time until we reach the target conversion
      for (int sampleIdx = 0; sampleIdx < 2 * nrOfStrands; ++sampleIdx) {
        if (this->currentCrosslinkerConversion >= targetCrossLinkerConversion) {
          break;
        }

        size_t strandIdx = availableStrandEnds[sampleIdx];

        if (strandEnd1[strandIdx] != -1) {
          // we don't have free cross-link choice
          assert(strandEnd2[strandIdx] == -1);

          size_t partnerCrosslinker = this->findAppropriateLink(
            strandEnd1[strandIdx],
            beadsPerChains[strandIdx] *
              (this->beadDistance * this->beadDistance),
            -1. // beadsPerChains[strandIdx] * this->beadDistance
          );

          strandEnd2[strandIdx] = partnerCrosslinker;
          this->remainingCrossLinkerFunctionality[partnerCrosslinker] -= 1;
          this->currentCrosslinkerConversion += conversionPerBond;
        } else {
          // otherwise, randomly choose a free cross-link
          long int crosslinkIdxIdx = crosslinkerIdxIdxDist(this->rng);
          while (this->remainingCrossLinkerFunctionality[crosslinkIdxIdx] < 1) {
            // find the next "available" cross-linker
            // TODO: while in practice this is not a problem,
            // in theory there are ways to optimize this for larger systems
            crosslinkIdxIdx += 1;
            if (crosslinkIdxIdx >= nCrosslinks) {
              crosslinkIdxIdx = 0;
            }
          }
          strandEnd1[strandIdx] = this->crossLinkerIdxs[crosslinkIdxIdx];
          this->remainingCrossLinkerFunctionality[crosslinkIdxIdx] -= 1;
          this->currentCrosslinkerConversion += conversionPerBond;
        }
      }

      // now that we know which strands should be connected,
      // we can do the connection
      for (int strandIdx = 0; strandIdx < nrOfStrands; ++strandIdx) {
        if (strandEnd1[strandIdx] == -1) {
          assert(strandEnd2[strandIdx] == -1);
          // free chain
          this->addSolventChains(1, beadsPerChains[strandIdx], strandAtomType);
        } else if (strandEnd2[strandIdx] == -1) {
          // dangling strand
          this->addRandomWalkChainFrom(
            strandEnd1[strandIdx], beadsPerChains[strandIdx], strandAtomType);
        } else {
          assert(strandEnd1[strandIdx] != -1 && strandEnd2[strandIdx] != -1);
          this->addRandomWalkChainFromTo(strandEnd1[strandIdx],
                                         strandEnd2[strandIdx],
                                         beadsPerChains[strandIdx],
                                         strandAtomType);
        }
      }
    };

    void addAndLinkStrands(int nrOfStrands,
                           int chainLength,
                           double crossLinkerConversion,
                           int strandAtomType = 1)
    {
      std::vector<int> chainLengths;
      chainLengths.reserve(nrOfStrands);
      for (int i = 0; i < nrOfStrands; ++i) {
        chainLengths.push_back(chainLength);
      }
      return this->addAndLinkStrands(
        nrOfStrands, chainLengths, crossLinkerConversion, strandAtomType);
    };

  private:
    double beadDistance = 0.965;
    double currentCrosslinkerConversion = 0.0;
    int maximumAtomId = 1;
    std::mt19937 rng;
    std::uniform_real_distribution<double> distX;
    std::uniform_real_distribution<double> distY;
    std::uniform_real_distribution<double> distZ;
    std::uniform_real_distribution<double> distSelect;

    SimplifiedUniverse simplifiedUniverse;
    std::vector<size_t> crossLinkerIdxs;
    std::vector<int> remainingCrossLinkerFunctionality;
    pylimer_tools::entities::Box box;

    /**
     * @brief Do a random walk of certain length to add a chain
     *
     * @param from the starting Atom of the chain
     * @param chainLen the number of additional atoms to add to the chain
     * @param atomType the atom type of the atoms in the chain
     */
    void addRandomWalkChainFrom(size_t idxFrom, int chainLen, int atomType = 1)
    {
      // std::cout << "Doing random walk from" << std::endl;
      std::vector<double> xs;
      xs.reserve(chainLen);
      std::vector<double> ys;
      ys.reserve(chainLen);
      std::vector<double> zs;
      zs.reserve(chainLen);

      std::uniform_real_distribution<double> angleDistribution =
        std::uniform_real_distribution<double>(0, 2 * M_PI);

      double lastX = this->simplifiedUniverse.x[idxFrom];
      double lastY = this->simplifiedUniverse.y[idxFrom];
      double lastZ = this->simplifiedUniverse.z[idxFrom];

      for (int i = 0; i < chainLen; ++i) {
        const double alpha = angleDistribution(this->rng);
        const double beta = angleDistribution(this->rng);
        // coordinate system conversion: confirmation e.g. in
        // https://math.stackexchange.com/a/1385150/738831
        xs.push_back(lastX +
                     this->beadDistance * std::cos(beta) * std::sin(alpha));
        lastX = xs[i];
        ys.push_back(lastY +
                     this->beadDistance * std::cos(beta) * std::cos(alpha));
        lastY = ys[i];
        zs.push_back(lastZ + this->beadDistance * std::sin(beta));
        lastZ = zs[i];
      }

      Positions positions;
      positions.x = xs;
      positions.y = ys;
      positions.z = zs;
      std::vector<size_t> idxs =
        this->addAtomsWithType(chainLen, atomType, positions);
      // initalize some bond specific
      std::vector<size_t> bondsFrom;
      std::vector<size_t> bondsTo;
      bondsFrom.reserve(idxs.size());
      bondsTo.reserve(idxs.size());
      for (size_t idx : idxs) {
        bondsFrom.push_back(this->simplifiedUniverse.ids[idx]);
        bondsTo.push_back(this->simplifiedUniverse.ids[idx]);
      }
      std::vector<int> bondTypes = initializeWithValue(chainLen, 1);
      // make the first bond from the starting bead given
      bondsFrom.insert(bondsFrom.begin(),
                       this->simplifiedUniverse.ids[idxFrom]);
      bondsFrom.pop_back();
      // finally, add the bonds
      this->simplifiedUniverse.bondsFrom.reserve(
        this->simplifiedUniverse.bondsFrom.size() + bondsFrom.size());
      this->simplifiedUniverse.bondsFrom.insert(
        this->simplifiedUniverse.bondsFrom.end(),
        bondsFrom.begin(),
        bondsFrom.end());
      this->simplifiedUniverse.bondsTo.reserve(
        this->simplifiedUniverse.bondsTo.size() + bondsTo.size());
      this->simplifiedUniverse.bondsTo.insert(
        this->simplifiedUniverse.bondsTo.end(), bondsTo.begin(), bondsTo.end());
      // std::cout << "Done random walk from" << std::endl;
    }

    /**
     * @brief Do a random walk of certain length to add a chain from one to
     * another atom
     *
     * @param from the atom to start the random walk from
     * @param to the atom to end the random walk at
     * @param chainLen the number of atoms to add in between from and to
     * @param atomType the type of the atoms to add
     */
    void addRandomWalkChainFromTo(size_t from,
                                  size_t to,
                                  int chainLen,
                                  int atomType = 1)
    {
      double lastX = this->simplifiedUniverse.x[from];
      double lastY = this->simplifiedUniverse.y[from];
      double lastZ = this->simplifiedUniverse.z[from];

      // support crossing of boundary conditions: find nearest image as target
      // (accept image mismatches)
      double targetX =
        (this->box.getLx() < std::sqrt((double)chainLen) * this->beadDistance)
          ? this->simplifiedUniverse.x[to]
          : lastX + this->_getDeltaDistance(
                      this->simplifiedUniverse.x[to], lastX, this->box.getLx());
      double targetY =
        (this->box.getLy() < std::sqrt((double)chainLen) * this->beadDistance)
          ? this->simplifiedUniverse.y[to]
          : lastY + this->_getDeltaDistance(
                      this->simplifiedUniverse.y[to], lastY, this->box.getLy());
      double targetZ =
        (this->box.getLz() < std::sqrt((double)chainLen) * this->beadDistance)
          ? this->simplifiedUniverse.z[to]
          : lastZ + this->_getDeltaDistance(
                      this->simplifiedUniverse.z[to], lastZ, this->box.getLz());

      std::unordered_map<std::string, std::vector<double>> walk_results =
        pylimer_tools::utils::doRandomWalkChainFromTo(
          this->box,
          std::array<double, 3>{
            lastX,
            lastY,
            lastZ,
          },
          std::array<double, 3>{ targetX, targetY, targetZ },
          chainLen,
          this->beadDistance);

      // assemeble and add these atoms
      Positions positions;
      positions.x = walk_results["x"];
      positions.y = walk_results["y"];
      positions.z = walk_results["z"];
      std::vector<size_t> idxs =
        this->addAtomsWithType(chainLen, atomType, positions);
      // initalize some bond specific
      std::vector<size_t> bondsFrom;
      std::vector<size_t> bondsTo;
      bondsFrom.reserve(idxs.size());
      bondsTo.reserve(idxs.size());
      for (size_t idx : idxs) {
        bondsFrom.push_back(this->simplifiedUniverse.ids[idx]);
        bondsTo.push_back(this->simplifiedUniverse.ids[idx]);
      }
      std::vector<int> bondTypes = initializeWithValue(chainLen, 1);
      // make the first bond from the starting bead given
      bondsFrom.insert(bondsFrom.begin(), this->simplifiedUniverse.ids[from]);
      // and the last bond further to the end of the chain
      bondsTo.insert(bondsTo.end(), this->simplifiedUniverse.ids[to]);
      // finally, actually add the bonds
      // finally, add the bonds
      this->simplifiedUniverse.bondsFrom.reserve(
        this->simplifiedUniverse.bondsFrom.size() + bondsFrom.size());
      this->simplifiedUniverse.bondsFrom.insert(
        this->simplifiedUniverse.bondsFrom.end(),
        bondsFrom.begin(),
        bondsFrom.end());
      this->simplifiedUniverse.bondsTo.reserve(
        this->simplifiedUniverse.bondsTo.size() + bondsTo.size());
      this->simplifiedUniverse.bondsTo.insert(
        this->simplifiedUniverse.bondsTo.end(), bondsTo.begin(), bondsTo.end());
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
                                         Positions randomPos)
    {
      this->simplifiedUniverse.ids.reserve(this->simplifiedUniverse.ids.size() +
                                           nrOfAtomsToAdd);
      this->simplifiedUniverse.types.reserve(
        this->simplifiedUniverse.types.size() + nrOfAtomsToAdd);
      this->simplifiedUniverse.x.reserve(this->simplifiedUniverse.x.size() +
                                         nrOfAtomsToAdd);
      this->simplifiedUniverse.y.reserve(this->simplifiedUniverse.y.size() +
                                         nrOfAtomsToAdd);
      this->simplifiedUniverse.z.reserve(this->simplifiedUniverse.z.size() +
                                         nrOfAtomsToAdd);

      int baseId = this->maximumAtomId + 1;
      std::vector<size_t> indicesAdded;
      indicesAdded.reserve(nrOfAtomsToAdd);
      size_t indexToStartWith = this->simplifiedUniverse.ids.size();

      for (size_t i = 0; i < nrOfAtomsToAdd; ++i) {
        this->simplifiedUniverse.ids.push_back(i + baseId);
        this->simplifiedUniverse.types.push_back(atomType);
        indicesAdded.push_back(indexToStartWith + i);
      }

      this->simplifiedUniverse.x.insert(this->simplifiedUniverse.x.end(),
                                        randomPos.x.begin(),
                                        randomPos.x.end());
      this->simplifiedUniverse.y.insert(this->simplifiedUniverse.y.end(),
                                        randomPos.y.begin(),
                                        randomPos.y.end());
      this->simplifiedUniverse.z.insert(this->simplifiedUniverse.z.end(),
                                        randomPos.z.begin(),
                                        randomPos.z.end());

      this->maximumAtomId += nrOfAtomsToAdd + 1;
      return indicesAdded;
    }

    /**
     * @brief Add atoms (incl. random positions, id etc.) to the universe
     *
     * @param nrOfAtomsToAdd the nr. of atoms to add to the universe
     * @param atomType the type of the atoms to add
     * @return std::vector<size_t> the ids of the inserted atoms
     */
    std::vector<size_t> addAtomsWithType(int nrOfAtomsToAdd, int atomType)
    {
      Positions randomPos = this->generateRandomPositions(nrOfAtomsToAdd);
      return this->addAtomsWithType(nrOfAtomsToAdd, atomType, randomPos);
    }

    /**
     * @brief Generate positions randomly
     *
     * @param nSamples the number of positions to generate
     * @return Positions
     */
    Positions generateRandomPositions(int nSamples)
    {
      if (nSamples < 1000) {
        return this->generateRandomBluePositions(nSamples);
      } else {
        return this->generateRandomWhitePositions(nSamples);
      }
    }

    /**
     * @brief Generate positions randomly (can lead to clustering)
     *
     * @param nSamples the number of positions to generate
     * @return Positions
     */
    Positions generateRandomWhitePositions(int nSamples)
    {
      std::vector<double> xs;
      std::vector<double> ys;
      std::vector<double> zs;

      for (size_t i = 0; i < nSamples; ++i) {
        double x = this->distX(this->rng);
        double y = this->distY(this->rng);
        double z = this->distZ(this->rng);

        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(z);
      }

      Positions result;
      result.x = xs;
      result.y = ys;
      result.z = zs;
      return result;
    }

    /**
     * @brief Generate positions randomly according to blue noise type sampling
     *
     * @param nSamples the number of positions to generate
     * @return Positions
     */
    Positions generateRandomBluePositions(int nSamples)
    {

      std::vector<double> xs;
      std::vector<double> ys;
      std::vector<double> zs;

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
          for (size_t k = 0; k < xs.size(); ++k) {
            double dist = this->getDistance(x, y, z, xs[k], ys[k], zs[k]);
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

        xs.push_back(bestCandidateX);
        ys.push_back(bestCandidateY);
        zs.push_back(bestCandidateZ);
      }

      Positions result;
      result.x = xs;
      result.y = ys;
      result.z = zs;
      return result;
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
      // TODO: use a neighbour list instead, maybe?
      for (int i = 0; i < this->crossLinkerIdxs.size(); ++i) {
        if (this->remainingCrossLinkerFunctionality[i] < 1) {
          continue;
        }
        size_t partner = this->crossLinkerIdxs[i];
        double distBetween = this->distanceBetween(partner, from);
        if (distBetween < maxDistance || maxDistance < 0.) {
          suitableMatches.push_back(partner);
          double thisWeight =
            std::exp(-3. * (distBetween) * (distBetween) / (desiredR02));
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

    double distanceBetween(size_t i, size_t j)
    {
      return this->getDistance(this->simplifiedUniverse.x[i],
                               this->simplifiedUniverse.y[i],
                               this->simplifiedUniverse.z[i],
                               this->simplifiedUniverse.x[j],
                               this->simplifiedUniverse.y[j],
                               this->simplifiedUniverse.z[j]);
    }

    double getDistance(double x1,
                       double y1,
                       double z1,
                       double x2,
                       double y2,
                       double z2)
    {
      double dx = this->_getDeltaDistance(x1, x2, this->box.getLx());
      double dy = this->_getDeltaDistance(y1, y2, this->box.getLy());
      double dz = this->_getDeltaDistance(z1, z2, this->box.getLz());
      return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    double _getDeltaDistance(double c1, double c2, double boxL) const
    {
      double delta = c1 - c2;
      while (delta > 0.5 * boxL) {
        delta -= boxL;
      }
      while (delta < -0.5 * boxL) {
        delta += boxL;
      }

      return delta;
    }
  };
} // namespace utils
} // namespace pylimer_tools

#endif
