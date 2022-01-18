#ifndef MC_UNIVERSE_GENERATOR_H
#define MC_UNIVERSE_GENERATOR_H

#include "../entities/Atom.h"
#include "../entities/Universe.h"
#include "StringUtils.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include <vector>

#include <random>

namespace pylimer_tools {
namespace utils {

struct Positions {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
};

class MCUniverseGenerator {
public:
  MCUniverseGenerator(const double Lx, const double Ly, const double Lz)
      : universe(pylimer_tools::entities::Universe(Lx, Ly, Lz)) {
    std::random_device rd;
    this->rng = std::mt19937(rd());
    this->distX = std::uniform_real_distribution<double>(0.0, Lx);
    this->distY = std::uniform_real_distribution<double>(0.0, Ly);
    this->distZ = std::uniform_real_distribution<double>(0.0, Lz);
  }
  void setSeed(unsigned int seed) { this->rng.seed(seed); };
  void setBeadDistance(double beadDistance) {
    this->beadDistance = beadDistance;
  };
  pylimer_tools::entities::Universe getUniverse() { return this->universe; };

  void addCrosslinkers(int nrOfCrosslinkers, int crosslinkerAtomType = 2) {
    if (this->crosslinkerType != crosslinkerAtomType) {
      if (this->crosslinkerType != 0) {
        throw std::invalid_argument(
            "Crosslinkers must all have the same type.");
      }
      this->crosslinkerType = crosslinkerAtomType;
    }
    this->addAtomsWithType(nrOfCrosslinkers, crosslinkerAtomType);
  };

  /**
   * @brief Randomly distribute 
   * 
   * @param nrOfSolventChains 
   * @param chainLength 
   * @param solventAtomType 
   */
  void addSolventChains(int nrOfSolventChains, int chainLength,
                        int solventAtomType = 3) {
    std::vector<long int> startAtoms =
        this->addAtomsWithType(nrOfSolventChains, solventAtomType);

    for (long int atomId : startAtoms) {
      pylimer_tools::entities::Atom startAtom =
          this->universe.getAtomByVertexIdx(
              this->universe.getIdxByAtomId(atomId));
      this->addRandomWalkChainFrom(startAtom, chainLength - 1, solventAtomType);
    }
  };

  /**
   * @brief Add strands in between the cross-linkers, link them as appropriate
   * 
   * @param nrOfStrands the nr. of Strands to add
   * @param beadsPerChains the nr. of beads per strand (excl. cross-linkers)
   * @param crosslinkerConversion "p", the target conversion of the cross-linkers
   * @param crosslinkerFunctionality the functionality of the cross-linker beads to add
   * @param strandAtomType the type of the strand atoms
   */
  void addAndLinkStrands(int nrOfStrands, std::vector<int> beadsPerChains,
                         double crosslinkerConversion,
                         int crosslinkerFunctionality = 4,
                         int strandAtomType = 1) {
    if (beadsPerChains.size() != nrOfStrands) {
      throw std::invalid_argument("Nr of strands must be equal to the number "
                                  "of chainLengths provided.");
    }

    std::uniform_real_distribution<double> linkerProbabilityDist =
        std::uniform_real_distribution<double>(0.0, 1.0);

    std::vector<pylimer_tools::entities::Atom> crosslinkers =
        this->universe.getAtomsWithType(this->crosslinkerType);
    std::vector<int> availableCrosslinkerSites =
        initializeWithValue(crosslinkers.size(), crosslinkerFunctionality);
    int nrOfStrandsAdded = 0;
    double conversionPerBond =
        1 / (crosslinkerFunctionality * crosslinkers.size());
    double currentDegreeOfConversion =
        this->universe
            .determineEffectiveFunctionalityPerType()[this->crosslinkerType];

    if (currentDegreeOfConversion != 0.0) {
      throw std::runtime_error(
          "Not implemented yet. Strands may only be linked once.");
    }

    // process all cross-linkers
    for (int i = 0; i < crosslinkers.size(); ++i) {
      const int availableSitesForThisCrosslinker = availableCrosslinkerSites[i];
      for (int siteToHandle = availableSitesForThisCrosslinker;
           siteToHandle > 0; --siteToHandle) {
        if (availableCrosslinkerSites[i] == 0) {
          break;
        }
        if (currentDegreeOfConversion >= crosslinkerConversion) {
          break;
        }
        // decide what type of site this is
        bool isActiveSite =
            linkerProbabilityDist(this->rng) < crosslinkerConversion;
        if (isActiveSite) {
          bool isDanglingStrand =
              linkerProbabilityDist(this->rng) < crosslinkerConversion;
          if (isDanglingStrand) {
            // decision is: dangling
            this->addRandomWalkChainFrom(crosslinkers[i],
                                         beadsPerChains[nrOfStrandsAdded]);
            nrOfStrandsAdded += 1;
            currentDegreeOfConversion += conversionPerBond;
          } else {
            // decision is: connected / network
            // find close enough crosslinker to do the chain to
            int targetIdx = this->findAppropriateLink(
                crosslinkers[i], crosslinkers, availableCrosslinkerSites,
                std::sqrt(beadsPerChains[i]) * this->beadDistance);
            this->addRandomWalkChainFromTo(crosslinkers[i],
                                           crosslinkers[targetIdx],
                                           beadsPerChains[nrOfStrandsAdded]);
            availableCrosslinkerSites[targetIdx] -= 1;
            nrOfStrandsAdded += 1;
            currentDegreeOfConversion += 2 * conversionPerBond;
          }

          if (nrOfStrandsAdded > nrOfStrands) {
            throw std::invalid_argument("Not enough strands to satisfy the "
                                        "requested degree of convergence.");
          }
        }
        availableCrosslinkerSites[i] -= 1;
      }
    }

    // process the remaining (free) chains
    for (size_t i = nrOfStrandsAdded; i < nrOfStrands; ++i) {
      this->addSolventChains(1, beadsPerChains[i], strandAtomType);
    }
  };

  void addAndLinkStrands(int nrOfStrands, int chainLength,
                         double crosslinkerConversion) {
    std::vector<int> chainLengths;
    chainLengths.reserve(nrOfStrands);
    for (int i = 0; i < nrOfStrands; ++i) {
      chainLengths.push_back(chainLength);
    }
    return this->addAndLinkStrands(nrOfStrands, chainLengths,
                                   crosslinkerConversion);
  };

private:
  double beadDistance = 0.965;
  int maximumAtomId = 1;
  std::mt19937 rng;
  int crosslinkerType = 0;
  std::uniform_real_distribution<double> distX;
  std::uniform_real_distribution<double> distY;
  std::uniform_real_distribution<double> distZ;
  pylimer_tools::entities::Universe universe;

  /**
   * @brief Do a random walk of certain length to add a chain
   *
   * @param from the starting Atom of the chain
   * @param chainLen the number of additional atoms to add to the chain
   * @param atomType the atom type of the atoms in the chain
   */
  void addRandomWalkChainFrom(pylimer_tools::entities::Atom from, int chainLen,
                              int atomType = 1) {
    std::vector<double> xs;
    xs.reserve(chainLen);
    std::vector<double> ys;
    ys.reserve(chainLen);
    std::vector<double> zs;
    zs.reserve(chainLen);

    std::uniform_real_distribution<double> angleDistribution =
        std::uniform_real_distribution<double>(0, 2 * M_PI);

    double lastX = from.getX();
    double lastY = from.getY();
    double lastZ = from.getZ();

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
    std::vector<long int> ids =
        this->addAtomsWithType(chainLen, atomType, positions);
    // initalize some bond specific
    std::vector<long int> bondsFrom = ids;
    std::vector<long int> bondsTo = ids;
    std::vector<int> bondTypes = initializeWithValue(chainLen, 1);
    // make the first bond from the starting bead given
    bondsFrom.insert(bondsFrom.begin(), from.getId());
    bondsFrom.pop_back();
    // finally, add the bonds
    this->universe.addBonds(chainLen, bondsFrom, bondsTo, bondTypes);
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
  void addRandomWalkChainFromTo(pylimer_tools::entities::Atom from,
                                pylimer_tools::entities::Atom to, int chainLen,
                                int atomType = 1) {
    std::vector<double> xs;
    xs.reserve(chainLen);
    std::vector<double> ys;
    ys.reserve(chainLen);
    std::vector<double> zs;
    zs.reserve(chainLen);

    std::uniform_real_distribution<double> angleDistribution =
        std::uniform_real_distribution<double>(0, 2 * M_PI);

    double lastX = from.getX();
    double lastY = from.getY();
    double lastZ = from.getZ();

    // support crossing of boundary conditions: find nearest image as target
    // (accept image mismatches)
    double targetX =
        lastX - this->_getDeltaDistance(to.getX(), lastX,
                                        this->universe.getBox().getLx());
    double targetY =
        lastY - this->_getDeltaDistance(to.getY(), lastY,
                                        this->universe.getBox().getLy());
    double targetZ =
        lastZ - this->_getDeltaDistance(to.getZ(), lastZ,
                                        this->universe.getBox().getLz());

    for (int i = 0; i < chainLen; ++i) {
      double dx = targetX - lastX;
      double dy = targetY - lastY;
      double dz = targetZ - lastZ;
      double idealAlpha = std::acos(dz);
      double idealBeta =
          dx > 0 ? (std::atan2((dy), (dx)))
                 : (dx < 0 ? (std::atan2((dy), (dx)) + M_PI) : M_PI / 2);
      double idealWeight =
          std::min(((double)i) / ((double)chainLen) +
                       (this->getDistance(lastX, lastY, lastZ, targetX, targetY,
                                          targetZ) /
                        (chainLen - i + 1)),
                   1.0);

      // TODO: find some a bit more sophisticated probability adjustment
      const double alpha = (1 - idealWeight) * angleDistribution(this->rng) +
                           idealWeight * idealAlpha;
      const double beta = (1 - idealWeight) * angleDistribution(this->rng) +
                          idealWeight * idealBeta;
      // coordinate system conversion: confirmation e.g. in
      // https://math.stackexchange.com/a/1385150/738831 or
      // https://en.wikipedia.org/wiki/Spherical_coordinate_system
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
    std::vector<long int> ids =
        this->addAtomsWithType(chainLen, atomType, positions);
    // initalize some bond specific
    std::vector<long int> bondsFrom = ids;
    std::vector<long int> bondsTo = ids;
    std::vector<int> bondTypes = initializeWithValue(chainLen, 1);
    // make the first bond from the starting bead given
    bondsFrom.insert(bondsFrom.begin(), from.getId());
    // and the last bond further to the end of the chain
    bondsTo.insert(bondsTo.end(), to.getId());
    // finally, actually add the bonds
    this->universe.addBonds(chainLen, bondsFrom, bondsTo, bondTypes);
  }

  /**
   * @brief Add atoms (incl. type, id etc.) with given positions to the universe
   *
   * @param nrOfAtomsToAdd the nr. of atoms to add to the universe
   * @param atomType the type of the atoms to add
   * @return std::vector<long int> the ids of the inserted atoms
   */
  std::vector<long int> addAtomsWithType(int nrOfAtomsToAdd, int atomType,
                                         Positions randomPos) {
    std::vector<long int> ids;
    std::vector<int> types;
    std::vector<int> zeros = initializeWithValue(nrOfAtomsToAdd, 0);
    ids.reserve(nrOfAtomsToAdd);
    types.reserve(nrOfAtomsToAdd);
    int baseId = this->maximumAtomId + 1;

    for (int i = 0; i < nrOfAtomsToAdd; ++i) {
      ids.push_back(i + baseId);
      types.push_back(atomType);
    }

    this->universe.addAtoms(nrOfAtomsToAdd, ids, types, randomPos.x,
                            randomPos.y, randomPos.z, zeros, zeros, zeros);
    this->maximumAtomId += nrOfAtomsToAdd + 1;
    return ids;
  }

  /**
   * @brief Add atoms (incl. random positions, id etc.) to the universe
   *
   * @param nrOfAtomsToAdd the nr. of atoms to add to the universe
   * @param atomType the type of the atoms to add
   * @return std::vector<long int> the ids of the inserted atoms
   */
  std::vector<long int> addAtomsWithType(int nrOfAtomsToAdd, int atomType) {
    Positions randomPos = this->generateRandomPositions(nrOfAtomsToAdd);
    return this->addAtomsWithType(nrOfAtomsToAdd, atomType, randomPos);
  }

  /**
   * @brief Generate positions randomly
   *
   * @param nSamples the number of positions to generate
   * @return Positions
   */
  Positions generateRandomPositions(int nSamples) {

    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;

    // blue noise
    // inspiration:
    // https://github.com/Atrix256/RandomCode/blob/master/Mitchell/Source.cpp
    for (size_t i = 0; i < nSamples; ++i) {
      size_t numCandidates =
          i * 1 + 1; // decrease the multiplier to speed things up
      double bestDistance = 0.0;
      double bestCandidateX = 0.0;
      double bestCandidateY = 0.0;
      double bestCandidateZ = 0.0;

      for (size_t j = 0; j < numCandidates; ++j) {
        double x = this->distX(this->rng);
        double y = this->distY(this->rng);
        double z = this->distZ(this->rng);

        double minDistance = this->universe.getVolume();
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
   * TODO: add some probability to not only select the very best.
   *
   * @param from the Atom to start the distance from
   * @param possiblePartners the Atoms that could be the target
   * @param acceptableDistance the distance to target if possible
   * @return int the index in possiblePartners that matches best
   */
  int findAppropriateLink(
      pylimer_tools::entities::Atom from,
      std::vector<pylimer_tools::entities::Atom> possiblePartners,
      std::vector<int> availablePartnerSites, double acceptableDistance) {
    if (possiblePartners.size() == 0) {
      throw std::invalid_argument("Cannot find a partner in none.");
    }
    double bestDistance = std::numeric_limits<double>::max();
    pylimer_tools::entities::Box box = this->universe.getBox();
    int bestMatch = 0;
    for (int i = 0; i < possiblePartners.size(); ++i) {
      if (availablePartnerSites[i] < 1) {
        continue;
      }
      pylimer_tools::entities::Atom partner = possiblePartners[i];
      double currDist =
          std::abs(partner.distanceTo(from, &box) - acceptableDistance);
      if (currDist < bestDistance) {
        bestDistance = currDist;
        bestMatch = i;
      }
    }
    return bestMatch;
  }

  double getDistance(double x1, double y1, double z1, double x2, double y2,
                     double z2) {
    double dx =
        this->_getDeltaDistance(x1, x2, this->universe.getBox().getLx());
    double dy =
        this->_getDeltaDistance(y1, y2, this->universe.getBox().getLy());
    double dz =
        this->_getDeltaDistance(z1, z2, this->universe.getBox().getLz());
    return dx * dx + dy * dy + dz * dz;
  }

  double _getDeltaDistance(double c1, double c2, double boxL) const {
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
