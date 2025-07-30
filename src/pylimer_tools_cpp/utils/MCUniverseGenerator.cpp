#include "MCUniverseGenerator.h"

namespace pylimer_tools::utils {
void
MCUniverseGenerator::setBeadDistance(double newBeadDistance,
                                     bool updateMeanSquared)
{
  INVALIDARG_EXP_IFN(newBeadDistance > 0, "Invalid mean bead distance");
  INVALIDARG_EXP_IFN(std::isfinite(newBeadDistance),
                     "Invalid mean bead distance");
  this->beadDistance = newBeadDistance;
  if (updateMeanSquared) {
    this->meanSquaredBeadDistance =
      (3. / 8.) * M_PI * SQUARE(this->beadDistance);
    RUNTIME_EXP_IFN(this->meanSquaredBeadDistance > 0,
                    "Invalid mean squared bead distance");
  }
}

void
MCUniverseGenerator::setMeanSquaredBeadDistance(
  double newMeanSquaredBeadDistance,
  bool updateMean)
{
  INVALIDARG_EXP_IFN(newMeanSquaredBeadDistance > 0,
                     "Invalid mean squared bead distance");
  INVALIDARG_EXP_IFN(std::isfinite(newMeanSquaredBeadDistance),
                     "Invalid mean squared bead distance");
  this->meanSquaredBeadDistance = newMeanSquaredBeadDistance;
  if (updateMean) {
    this->beadDistance =
      std::sqrt(this->meanSquaredBeadDistance / ((3. / 8.) * M_PI));
    RUNTIME_EXP_IFN(this->beadDistance > 0, "Invalid mean bead distance");
  }
}

void
MCUniverseGenerator::configPrimaryLoopProbability(
  double newPrimaryLoopProbability)
{
  INVALIDARG_EXP_IFN(newPrimaryLoopProbability >= 0,
                     "Invalid primary loop formation probability");
  this->primaryLoopProbability = newPrimaryLoopProbability;
}

void
MCUniverseGenerator::configSecondaryLoopProbability(
  double newSecondaryLoopProbability)
{
  INVALIDARG_EXP_IFN(newSecondaryLoopProbability >= 0,
                     "Invalid secondary loop formation probability");
  this->secondaryLoopProbability = newSecondaryLoopProbability;
}

void
MCUniverseGenerator::disableMaxDistance()
{
  this->configMaxDistanceProvider(std::make_unique<NoMaxDistanceProvider>());
}

void
MCUniverseGenerator::useLinearMaxDistance(double newMultiplier)
{
  this->configMaxDistanceProvider(
    std::make_unique<LinearMaxDistanceProvider>(newMultiplier));
}

void
MCUniverseGenerator::useZScoreMaxDistance(double newStdMultiplier,
                                          double innerMultiplier)
{
  this->configMaxDistanceProvider(std::make_unique<ZScoreMaxDistanceProvider>(
    newStdMultiplier, innerMultiplier));
}

void
MCUniverseGenerator::configMaxDistanceProvider(
  std::unique_ptr<MaxDistanceProvider> newMaxDistanceProvider)
{
  this->maxDistanceProvider = std::move(newMaxDistanceProvider);
  this->resetNeighbourList();
}
entities::Universe
MCUniverseGenerator::getUniverse()
{
#ifndef NDEBUG
  this->validateInternalState();
#endif

  entities::Universe universe = entities::Universe(this->box);
  size_t nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
  long int nrOfAtoms = this->getCurrentNrOfAtoms();
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

  // Add the crosslinkers
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
  for (size_t strandI = 0; strandI < this->simplifiedUniverse.strandFrom.size();
       ++strandI) {
    this->beadDistance = this->simplifiedUniverse.beadDistanceInStrand[strandI];
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
      RUNTIME_EXP_IFN(nBeadsInStrand >= 1,
                      "Cannot generate free chain with 0 beads");
      coordinates = this->sampleFreeChainCoordinates(nBeadsInStrand);
    } else if (strandEnd2 < 0) {
      RUNTIME_EXP_IFN(nBeadsInStrand >= 1,
                      "Cannot generate dangling chain with 0 beads");
      coordinates =
        this->sampleDanglingChainCoordinates(strandEnd1, nBeadsInStrand);
      bondsFrom.push_back(strandEnd1 + 1);
      bondsTo.push_back(currentId);
    } else {
      RUNTIME_EXP_IFN(strandEnd1 != -1 && strandEnd2 != -1,
                      "Expected both ends to have an associated crosslinker");

      if (nBeadsInStrand > 0) {
        coordinates =
          this->sampleStrandCoordinates(strandEnd1, strandEnd2, nBeadsInStrand);
        bondsFrom.push_back(strandEnd1 + 1);
        bondsTo.push_back(currentId);
        bondsFrom.push_back(strandEnd2 + 1);
        bondsTo.push_back(currentId + nBeadsInStrand - 1);
      } else {
        RUNTIME_EXP_IFN(nBeadsInStrand == 0,
                        "Cannot generate chain with " +
                          std::to_string(nBeadsInStrand) + " beads.");
        bondsFrom.push_back(strandEnd1 + 1);
        bondsTo.push_back(strandEnd2 + 1);
      }
    }

    // actually add the new beads to our list of things to add
    if (nBeadsInStrand > 0) {
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
  }

  assert(xs.size() == ys.size() && xs.size() == zs.size() &&
         xs.size() == ids.size() && types.size() == nrOfAtoms);
  this->beadDistance = prevBeadDistance;
  this->meanSquaredBeadDistance = prevMeanSquaredBeadDistance;

  universe.addAtoms(ids, types, xs, ys, zs, zeros, zeros, zeros);
  universe.addBonds(bondsFrom, bondsTo);
  return universe;
}
void
MCUniverseGenerator::addCrosslinkersAt(Eigen::VectorXd coordinates,
                                       int crosslinkerFunctionality,
                                       int crossLinkerAtomType)
{
  INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                     "Length of coordinates must be a multiple of 3");
  INVALIDARG_EXP_IFN(crosslinkerFunctionality >= 0,
                     "Expecting positive crosslinker functionality, got " +
                       std::to_string(crosslinkerFunctionality) + ".");
  int nCrosslinkerBefore = this->remainingCrossLinkerFunctionality.size();

  size_t nrOfCrosslinkers = coordinates.size() / 3;

  this->addXlinkAtoms(nrOfCrosslinkers, crossLinkerAtomType, coordinates);

  this->remainingCrossLinkerFunctionality.reserve(nCrosslinkerBefore +
                                                  nrOfCrosslinkers);
  this->originalNrOfAvailableCrosslinkSites +=
    crosslinkerFunctionality * nrOfCrosslinkers;
  this->nrOfAvailableCrosslinkSites +=
    crosslinkerFunctionality * nrOfCrosslinkers;
  for (size_t i = 0; i < nrOfCrosslinkers; ++i) {
    this->remainingCrossLinkerFunctionality.push_back(crosslinkerFunctionality);
  }
  this->updateNeighbourListCoordinates();
#ifndef NDEBUG
  this->validateInternalState();
#endif
}
void
MCUniverseGenerator::addCrosslinkers(int nrOfCrosslinkers,
                                     int crosslinkerFunctionality,
                                     int crossLinkerAtomType,
                                     bool whiteNoise)
{
  Eigen::VectorXd randomPos =
    this->generateRandomPositions(nrOfCrosslinkers, whiteNoise);

  this->addCrosslinkersAt(
    randomPos, crosslinkerFunctionality, crossLinkerAtomType);
}
void
MCUniverseGenerator::addRandomlyFunctionalizedStrands(
  int nrOfStrands,
  std::vector<int> beadsPerStrand,
  double functionalizationProbability,
  int crosslinkerFunctionality,
  int crosslinkerAtomType,
  int strandAtomType,
  bool whiteNoise)
{
  INVALIDARG_EXP_IFN(nrOfStrands > 0, "Cannot add 0 or less strands");
  INVALIDARG_EXP_IFN(
    functionalizationProbability >= 0.,
    "Functionalization probability must be larger or equal to 0.");
  INVALIDARG_EXP_IFN(crosslinkerFunctionality >= 1,
                     "Crosslinker functionality must be at least 1. Use "
                     "normal chains instead for lower functionalities.");
  INVALIDARG_EXP_IFN(
    crosslinkerFunctionality == 1 || functionalizationProbability <= 1.,
    "Functionalization probability must be smaller than 1, or the "
    "crosslinker functionality must be 1.");
  INVALIDARG_EXP_IFN(nrOfStrands == beadsPerStrand.size(),
                     "Inconsistent sizes");
#ifndef NDEBUG
  this->validateInternalState();
#endif

  size_t nCrosslinksBefore = this->remainingCrossLinkerFunctionality.size();
  size_t nStrandsBefore = this->simplifiedUniverse.strandFrom.size();

  std::uniform_real_distribution<double> randomDist(0., 1.);
  size_t nCrosslinks = 0;
  size_t nEffectiveStrands = 0;

  std::vector<int> strandIdOfCrosslink;
  std::vector<int> newCrosslinkFunctionality;

  int nRepeat = 1;
  if (functionalizationProbability > 1.) {
    nRepeat = static_cast<int>(std::ceil(functionalizationProbability));
    functionalizationProbability /= static_cast<double>(nRepeat);
  }

  for (size_t strandI = 0; strandI < nrOfStrands; ++strandI) {
    std::vector<int> partialStrandLengths;
    std::vector<size_t> crosslinkIndicesInStrand;
    long int lastSampledBead = -1;
    for (long int i = 0; i < beadsPerStrand[strandI]; ++i) {
      bool convertThisBead = false;
      int functionality = 0;
      for (size_t j = 0; j < nRepeat; ++j) {
        if (randomDist(this->rng) < functionalizationProbability) {
          convertThisBead = true;
          functionality += crosslinkerFunctionality;
        }
      }

      if (convertThisBead) {
        // yes, we want to replace this bead `i` with a crosslink
        size_t currentSpringIdx = nStrandsBefore + nEffectiveStrands;
        long int lengthToPrevious =
          i - (lastSampledBead >= 0 ? (lastSampledBead + 1) : 0);
        if (i > 0) {
          this->addStrand(lengthToPrevious,
                          strandAtomType,
                          UNCONNECTED,
                          lastSampledBead < 0 ? EMPTY_BACKGROUND : UNCONNECTED);
          nEffectiveStrands += 1;
        }

        // order matters here for the displacements afterwards
        if (lastSampledBead >= 0) {
          this->linkStrandToCrosslink(
            currentSpringIdx, nCrosslinksBefore + nCrosslinks - 1, true);
        }

        if (i > 0) {
          this->linkStrandToCrosslink(
            currentSpringIdx, nCrosslinksBefore + nCrosslinks, true);
        }

        strandIdOfCrosslink.push_back(nCrosslinksBefore + strandI);
        newCrosslinkFunctionality.push_back(functionality);
        nCrosslinks += 1;
        lastSampledBead = i;
      }
    }

    if (lastSampledBead < beadsPerStrand[strandI] - 1) {
      // add dangling spring
      // possibly, this is the full spring if no crosslink was sampled
      size_t currentSpringIdx = nStrandsBefore + nEffectiveStrands;
      assert(this->simplifiedUniverse.strandTo.size() == currentSpringIdx);

      long int remainingLength =
        beadsPerStrand[strandI] -
        (lastSampledBead >= 0 ? (lastSampledBead + 1) : 0);
      this->addStrand(remainingLength,
                      strandAtomType,
                      lastSampledBead < 0 ? EMPTY_BACKGROUND : UNCONNECTED,
                      EMPTY_BACKGROUND);
      if (lastSampledBead >= 0) {
        this->linkStrandToCrosslink(
          currentSpringIdx, nCrosslinksBefore + nCrosslinks - 1, true);
      }

      nEffectiveStrands += 1;
    }
  }

  // then, add the crosslink atoms
  // distances need to be taken into account, are handled later
  RUNTIME_EXP_IFN(nCrosslinks == strandIdOfCrosslink.size(),
                  "Did not register the expected number of crosslinks.");
  RUNTIME_EXP_IFN(nCrosslinks == newCrosslinkFunctionality.size(),
                  "Did not register the expected number of crosslinks.");
  this->addCrosslinkers(
    nCrosslinks, crosslinkerFunctionality, crosslinkerAtomType, whiteNoise);
  for (size_t newXlinOffset = 0; newXlinOffset < nCrosslinks; ++newXlinOffset) {
    this->simplifiedUniverse.xlinkChainId[nCrosslinksBefore + newXlinOffset] =
      strandIdOfCrosslink[newXlinOffset];
    int functionalityDifferenceToDesired =
      newCrosslinkFunctionality[newXlinOffset] -
      this
        ->remainingCrossLinkerFunctionality[nCrosslinksBefore + newXlinOffset];
    this->originalNrOfAvailableCrosslinkSites +=
      functionalityDifferenceToDesired;
    this->nrOfAvailableCrosslinkSites += functionalityDifferenceToDesired;
    this
      ->remainingCrossLinkerFunctionality[nCrosslinksBefore + newXlinOffset] +=
      functionalityDifferenceToDesired;
  }

  // actually link the crosslinks to the crosslink strand
  for (size_t newStrandIdx = nStrandsBefore;
       newStrandIdx < nStrandsBefore + nEffectiveStrands;
       ++newStrandIdx) {
    long int from = this->simplifiedUniverse.strandFrom[newStrandIdx];
    long int to = this->simplifiedUniverse.strandTo[newStrandIdx];
    assert(from < static_cast<long int>(nCrosslinksBefore + nCrosslinks));
    assert(to < static_cast<long int>(nCrosslinksBefore + nCrosslinks));
    // do the linking that was omitted earlier
    if (from >= 0) {
      this->simplifiedUniverse.strandsOfXlink[from].push_back(newStrandIdx);
    }
    if (to >= 0) {
      this->simplifiedUniverse.strandsOfXlink[to].push_back(newStrandIdx);
    }

    if (from >= 0 && to >= 0) {
      // finally, make sure that the coordinates of two subsequent
      // crosslinks match the required length of the respective strand
      Eigen::VectorXd delta = this->sampleCoordinatesWithinNBeadDistance(
        this->simplifiedUniverse.beadsInStrand[newStrandIdx]);

      this->simplifiedUniverse.xlinkX[to] =
        this->simplifiedUniverse.xlinkX[from] + delta[0];
      this->simplifiedUniverse.xlinkY[to] =
        this->simplifiedUniverse.xlinkY[from] + delta[1];
      this->simplifiedUniverse.xlinkZ[to] =
        this->simplifiedUniverse.xlinkZ[from] + delta[2];

      // validate the distances
#ifndef NDEBUG
      double distance = this->distanceBetween(from, to);
      double maxDistance =
        (this->simplifiedUniverse.beadsInStrand[newStrandIdx] + 2) *
        this->beadDistance;
      if (maxDistance > 1.) {
        RUNTIME_EXP_IFN(
          distance < maxDistance * 2.,
          "Distance adjustment does not seem to be correct: got " +
            std::to_string(distance) + " for chain with " +
            std::to_string(
              this->simplifiedUniverse.beadsInStrand[newStrandIdx]) +
            " beads and " + std::to_string(this->beadDistance) + " distance.");
      }
#endif
    }
  }

  this->updateNeighbourListCoordinates();
  this->validateInternalState();
}
void
MCUniverseGenerator::addCrosslinkStrands(int nrOfCrosslinkStrands,
                                         std::vector<int> beadsPerStrand,
                                         int crosslinkerFunctionality,
                                         int crosslinkerAtomType,
                                         int strandAtomType,
                                         bool whiteNoise)
{
  INVALIDARG_EXP_IFN(nrOfCrosslinkStrands > 0, "");
  INVALIDARG_EXP_IFN(crosslinkerFunctionality >= 2,
                     "Crosslinker functionality must be at least 2. Use "
                     "solvent chains instead for lower functionalities.");
  INVALIDARG_EXP_IFN(nrOfCrosslinkStrands == beadsPerStrand.size(),
                     "Inconsistent sizes");
#ifndef NDEBUG
  this->validateInternalState();
#endif
  size_t nCrosslinksBefore = this->remainingCrossLinkerFunctionality.size();
  size_t nStrandsBefore = this->simplifiedUniverse.strandFrom.size();

  // start with adding the ends
  this->addCrosslinkers(2 * nrOfCrosslinkStrands,
                        crosslinkerFunctionality,
                        crosslinkerAtomType,
                        whiteNoise);

  // then, add the middle strands
  this->addStrands(nrOfCrosslinkStrands, beadsPerStrand, strandAtomType);

  // finally, do the linking and adjustment of the positions
  for (size_t i = 0; i < nrOfCrosslinkStrands; ++i) {
    size_t strandIdx = nStrandsBefore + i;
    size_t from = nCrosslinksBefore + i * 2;
    size_t to = nCrosslinksBefore + i * 2 + 1;
    // actually link these strands
    this->linkStrandToCrosslink(strandIdx, from, false);
    this->linkStrandToCrosslink(strandIdx, to, false);
    // make sure the positions make sense
    Eigen::VectorXd delta =
      this->sampleCoordinatesWithinNBeadDistance(beadsPerStrand[i]);

    this->simplifiedUniverse.xlinkX[to] =
      this->simplifiedUniverse.xlinkX[from] + delta[0];
    this->simplifiedUniverse.xlinkY[to] =
      this->simplifiedUniverse.xlinkY[from] + delta[1];
    this->simplifiedUniverse.xlinkZ[to] =
      this->simplifiedUniverse.xlinkZ[from] + delta[2];
    this->simplifiedUniverse.xlinkChainId[to] = from;
  }

  this->updateNeighbourListCoordinates();
  this->validateInternalState();
}
void
MCUniverseGenerator::addSolventChains(int nrOfSolventChains,
                                      int chainLength,
                                      int solventAtomType,
                                      bool whiteNoise)
{
  for (size_t i = 0; i < nrOfSolventChains; ++i) {
    this->simplifiedUniverse.strandFrom.push_back(EMPTY_BACKGROUND);
    this->simplifiedUniverse.strandTo.push_back(EMPTY_BACKGROUND);
    this->simplifiedUniverse.strandBeadType.push_back(solventAtomType);
    this->simplifiedUniverse.beadsInStrand.push_back(chainLength);
    this->simplifiedUniverse.beadDistanceInStrand.push_back(this->beadDistance);
    this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.push_back(
      this->meanSquaredBeadDistance);
  }
}
void
MCUniverseGenerator::addMonofunctionalStrands(int nrOfStrands,
                                              std::vector<int> beadsPerChains,
                                              int strandAtomType)
{
  INVALIDARG_EXP_IFN(beadsPerChains.size() == nrOfStrands,
                     "Nr of strands (" + std::to_string(nrOfStrands) +
                       ") must be equal to the number "
                       "of chainLengths (" +
                       std::to_string(beadsPerChains.size()) + ") provided.");

  for (size_t i = 0; i < nrOfStrands; ++i) {
    this->simplifiedUniverse.strandBeadType.push_back(strandAtomType);
    this->simplifiedUniverse.beadsInStrand.push_back(beadsPerChains[i]);
    this->simplifiedUniverse.beadDistanceInStrand.push_back(this->beadDistance);
    this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.push_back(
      this->meanSquaredBeadDistance);
    this->simplifiedUniverse.strandFrom.push_back(UNCONNECTED);
    this->simplifiedUniverse.strandTo.push_back(EMPTY_BACKGROUND);
  }

#ifndef NDEBUG
  this->validateInternalState();
#endif
}
void
MCUniverseGenerator::addMonofunctionalStrands(int nrOfStrands,
                                              int chainLength,
                                              int strandAtomType)
{
  const std::vector<int> chainLengths =
    pylimer_tools::utils::initializeWithValue<int>(nrOfStrands, chainLength);
  return this->addMonofunctionalStrands(
    nrOfStrands, chainLengths, strandAtomType);
}
void
MCUniverseGenerator::addStrand(const int beadsOfChain,
                               const int strandAtomType,
                               const long int connectionFrom,
                               const long int connectionTo)
{
  INVALIDARG_EXP_IFN(beadsOfChain >= 0, "Beads per chain must be positive.");
  INVALIDARG_EXP_IFN(
    connectionFrom < 0 && connectionTo < 0,
    "Only unconnected or dangling strands can be added with this method.");

  this->simplifiedUniverse.strandBeadType.push_back(strandAtomType);
  this->simplifiedUniverse.beadsInStrand.push_back(beadsOfChain);
  this->simplifiedUniverse.beadDistanceInStrand.push_back(this->beadDistance);
  this->simplifiedUniverse.meanSquaredBeadDistanceInStrand.push_back(
    this->meanSquaredBeadDistance);
  this->simplifiedUniverse.strandFrom.push_back(connectionFrom);
  this->simplifiedUniverse.strandTo.push_back(connectionTo);
}
void
MCUniverseGenerator::addStrands(const int nrOfStrands,
                                const std::vector<int> beadsPerChains,
                                const int strandAtomType)
{
  INVALIDARG_EXP_IFN(
    beadsPerChains.size() == nrOfStrands,
    "Inconsistent nr of strands and nr of beads per strand given.");

  for (size_t strandIdx = 0; strandIdx < nrOfStrands; ++strandIdx) {
    this->addStrand(beadsPerChains[strandIdx], strandAtomType);
  }
}
void
MCUniverseGenerator::linkStrand(const size_t strandIdx, const double cInfinity)
{
  // validation
  INVALIDARG_EXP_IFN(strandIdx >= 0 &&
                       strandIdx < this->simplifiedUniverse.strandFrom.size(),
                     "Strand index out of range.");
  INVALIDARG_EXP_IFN(
    this->simplifiedUniverse.strandTo[strandIdx] < 0,
    "Expected second strand end of strand " + std::to_string(strandIdx) +
      " to be free, got " +
      std::to_string(this->simplifiedUniverse.strandTo[strandIdx]) +
      " for strand " + std::to_string(strandIdx) + ".");
  RUNTIME_EXP_IFN(this->nrOfAvailableCrosslinkSites > 0,
                  "No available crosslink sites left.");

  // actual linking
  if (this->simplifiedUniverse.strandFrom[strandIdx] >= 0) {
    INVALIDARG_EXP_IFN(
      this->simplifiedUniverse.strandTo[strandIdx] == UNCONNECTED,
      "Expected second strand end of strand " + std::to_string(strandIdx) +
        " to be free, got " +
        std::to_string(this->simplifiedUniverse.strandTo[strandIdx]) +
        " for strand " + std::to_string(strandIdx) + ".");
    const double timesNForR02 =
      this->simplifiedUniverse.meanSquaredBeadDistanceInStrand[strandIdx] *
      cInfinity;
    // we don't have free crosslink choice
    // find one that follows the desired end-to-end distribution
    long int partnerCrosslinker = this->findAppropriateLink(
      this->simplifiedUniverse.strandFrom[strandIdx],
      static_cast<double>(this->simplifiedUniverse.beadsInStrand[strandIdx] +
                          1) *
        timesNForR02,
      this->maxDistanceProvider->getMaxDistance(static_cast<double>(
        this->simplifiedUniverse.beadsInStrand[strandIdx] + 1)));

    RUNTIME_EXP_IFN(partnerCrosslinker >= 0,
                    "No suitable crosslink partner found.");

    this->simplifiedUniverse.strandTo[strandIdx] = partnerCrosslinker;
    this->remainingCrossLinkerFunctionality[partnerCrosslinker] -= 1;
    this->simplifiedUniverse.strandsOfXlink[partnerCrosslinker].push_back(
      strandIdx);
    this->nrOfAvailableCrosslinkSites -= 1;
  } else {
    // sample random crosslink site
    std::discrete_distribution<long int> xlinkIdxDist(
      this->remainingCrossLinkerFunctionality.begin(),
      this->remainingCrossLinkerFunctionality.end());

    long int matchingCrosslink = xlinkIdxDist(this->rng);

    // link to this crosslink
    this->simplifiedUniverse.strandFrom[strandIdx] = matchingCrosslink;
    this->remainingCrossLinkerFunctionality[matchingCrosslink] -= 1;
    this->simplifiedUniverse.strandsOfXlink[matchingCrosslink].push_back(
      strandIdx);
    this->nrOfAvailableCrosslinkSites -= 1;
  }
}
void
MCUniverseGenerator::linkStrandsCallback(
  std::function<BackTrackStatus(const MCUniverseGenerator&, long int)>
    linkingController,
  double cInfinity)
{
#ifndef NDEBUG
  this->validateInternalState();
#endif

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
    for (long int s = 0; s < this->remainingCrossLinkerFunctionality[i]; ++s) {
      availableCrosslinkSites.push_back(i);
    }
  }
  std::shuffle(
    availableCrosslinkSites.begin(), availableCrosslinkSites.end(), this->rng);

  std::stack<size_t> removedCrosslinkSites;
  // removedCrosslinkSites.reserve(availableCrosslinkSites.size());

  // link one strand at a time until we reach the target conversion
  size_t maxStep =
    std::min(availableStrandEnds.size(), availableCrosslinkSites.size());
  for (size_t sampleIdx = 0; sampleIdx < maxStep; ++sampleIdx) {
    BackTrackStatus status = linkingController(
      *this, static_cast<long int>(maxStep) - static_cast<long int>(sampleIdx));
    if (status == BackTrackStatus::STOP) {
      break;
    } else if (status == BackTrackStatus::TRACK_FORWARD) {
      size_t strandIdx = availableStrandEnds[sampleIdx];
      const double timesNForR02 =
        this->simplifiedUniverse.meanSquaredBeadDistanceInStrand[strandIdx] *
        cInfinity;
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
        // we don't have free crosslink choice
        // find one that follows the desired end-to-end distribution
        long int partnerCrosslinker = this->findAppropriateLink(
          this->simplifiedUniverse.strandFrom[strandIdx],
          static_cast<double>(
            this->simplifiedUniverse.beadsInStrand[strandIdx] + 1) *
            timesNForR02,
          this->maxDistanceProvider->getMaxDistance(
            this->simplifiedUniverse.beadsInStrand[strandIdx] + 1));

        // it should be thought through again what to do in this case, where
        // there is no matching crosslink.
        // RUNTIME_EXP_IFN(partnerCrosslinker >= 0, "No suitable crosslink
        // partner found.");
        if (partnerCrosslinker < 0) {
          continue;
        }

        this->simplifiedUniverse.strandTo[strandIdx] = partnerCrosslinker;
        this->remainingCrossLinkerFunctionality[partnerCrosslinker] -= 1;
        this->simplifiedUniverse.strandsOfXlink[partnerCrosslinker].push_back(
          strandIdx);
        this->nrOfAvailableCrosslinkSites -= 1;
      } else {
        // otherwise, randomly choose a free crosslink
        long int crosslinkIdxIdx;
        do {
          // find the next "available" crosslink
          crosslinkIdxIdx = availableCrosslinkSites.back();
          removedCrosslinkSites.push(crosslinkIdxIdx);
          availableCrosslinkSites.pop_back();
        } while (this->remainingCrossLinkerFunctionality[crosslinkIdxIdx] < 1 &&
                 availableCrosslinkSites.size() > 0);

        if (availableCrosslinkSites.size() == 0) {
          std::cerr << "No more crosslink sites available." << std::endl;
          break;
        }

        // link to this crosslink
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
        // need to re-fill the available crosslink sites for the
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
void
MCUniverseGenerator::linkStrandsToConversion(
  const double targetCrossLinkerConversion,
  const double cInfinity)
{
  RUNTIME_EXP_IFN(this->originalNrOfAvailableCrosslinkSites > 0,
                  "No available crosslink sites.");

  const double conversionPerBond =
    (1.0) / (static_cast<double>(this->originalNrOfAvailableCrosslinkSites));
  const double currentCrosslinkerConversion =
    this->originalNrOfAvailableCrosslinkSites > 0
      ? (1.0 - (static_cast<double>(this->nrOfAvailableCrosslinkSites) /
                static_cast<double>(this->originalNrOfAvailableCrosslinkSites)))
      : 0.0;

  INVALIDARG_EXP_IFN(
    APPROX_WITHIN(
      targetCrossLinkerConversion, currentCrosslinkerConversion, 1., 1e-6),
    "Crosslinker conversion must be between " +
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
      "A crosslinker conversion of " +
      std::to_string(targetCrossLinkerConversion) + " is not reachable with " +
      std::to_string(potentialNewBonds) + " free strand ends and " +
      std::to_string(this->nrOfAvailableCrosslinkSites) + " of " +
      std::to_string(this->originalNrOfAvailableCrosslinkSites) +
      " available crosslink sites with current conversion of " +
      std::to_string(currentCrosslinkerConversion) +
      ". Maximum possible p is " +
      std::to_string(currentCrosslinkerConversion +
                     potentialNewBonds * conversionPerBond) +
      ".");
  }

  const long int targetNrOfAvailableCrosslinkSites =
    std::round((1.0 - targetCrossLinkerConversion) *
               static_cast<double>(this->originalNrOfAvailableCrosslinkSites));
  const double timesNForR02 = this->meanSquaredBeadDistance * cInfinity;

  this->linkStrandsCallback(
    [targetNrOfAvailableCrosslinkSites](const MCUniverseGenerator& gen,
                                        long int nStrandsRemaining) {
      if (gen.nrOfAvailableCrosslinkSites <=
          targetNrOfAvailableCrosslinkSites) {
        return BackTrackStatus::STOP;
      };
      return BackTrackStatus::TRACK_FORWARD;
    },
    cInfinity);
}
void
MCUniverseGenerator::linkStrandsToSolubleFraction(double targetSolubleFraction,
                                                  double cInfinity)
{
  INVALIDARG_EXP_IFN(APPROX_WITHIN(targetSolubleFraction, 0., 1., 1e-8),
                     "Soluble fraction must be between 0 and 1, got " +
                       std::to_string(targetSolubleFraction) + ".");

  size_t nAtomsTotal = this->getCurrentNrOfAtoms();

  BackTrackStatus status = BackTrackStatus::TRACK_FORWARD;
  long int nSteps = 0;
  long int lastStep = 0;
  long int currentStep = 0;

  this->linkStrandsCallback(
    [targetSolubleFraction,
     nAtomsTotal,
     &status,
     &nSteps,
     &lastStep,
     &currentStep](const MCUniverseGenerator& gen, long int nStrandsRemaining) {
      if (currentStep == 0) {
        // go all in, we want to jump to end to know when to stop
        nSteps = nStrandsRemaining;
      }
      currentStep += 1;
      nSteps = std::max<long int>(nSteps, (long int)1);
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
        forceRelaxer.runForceRelaxation("LD_MMA", 5000, 1e-11, 1e-8);
      }

      // finally, calculate the soluble fraction
      double solubleFraction = 1. - forceRelaxer.countActiveClusteredAtoms() /
                                      static_cast<double>(nAtomsTotal);
      std::cout << "Got w_sol = " << solubleFraction << " at step "
                << currentStep << " (+" << nSteps << ") with "
                << nStrandsRemaining << " strands remaining." << std::endl;
      if (solubleFraction > targetSolubleFraction) {
        nSteps /= 2;
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
void
MCUniverseGenerator::removeStrand(size_t strandIdx)
{
  INVALIDARG_EXP_IFN(strandIdx < this->simplifiedUniverse.strandFrom.size(),
                     "Strand to be removed is out of range.");

  if (this->simplifiedUniverse.strandFrom[strandIdx] >= 0) {
    // remove this strand from the crosslinker
    size_t xlink = this->simplifiedUniverse.strandFrom[strandIdx];
    pylimer_tools::utils::removeIfContained<long int>(
      this->simplifiedUniverse.strandsOfXlink[xlink], strandIdx);
    this->remainingCrossLinkerFunctionality[xlink] += 1;
    this->nrOfAvailableCrosslinkSites += 1;
  }
  if (this->simplifiedUniverse.strandTo[strandIdx] >= 0) {
    // remove this strand from the crosslinker
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

  // renumber the link from crosslinks to strands
  for (size_t i = 0; i < this->simplifiedUniverse.xlinkTypes.size(); ++i) {
    for (size_t j = 0; j < this->simplifiedUniverse.strandsOfXlink[i].size();
         ++j) {
      assert(this->simplifiedUniverse.strandsOfXlink[i][j] != strandIdx);
      if (this->simplifiedUniverse.strandsOfXlink[i][j] > strandIdx) {
        this->simplifiedUniverse.strandsOfXlink[i][j] -= 1;
      }
    }
  }
}
void
MCUniverseGenerator::removeCrosslink(size_t crosslinkIdx)
{
  INVALIDARG_EXP_IFN(crosslinkIdx < this->simplifiedUniverse.xlinkTypes.size(),
                     "Crosslink to be removed is out of range.");

  // re-connect the strands to account for the removed crosslink
  for (size_t i = 0; i < this->simplifiedUniverse.strandFrom.size(); ++i) {
    INVALIDARG_EXP_IFN(this->simplifiedUniverse.strandFrom[i] != crosslinkIdx &&
                         this->simplifiedUniverse.strandTo[i] != crosslinkIdx,
                       "The to-be removed crosslink is still connected.");
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
  this->simplifiedUniverse.xlinkChainId.erase(
    this->simplifiedUniverse.xlinkChainId.begin() + crosslinkIdx);
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
void
MCUniverseGenerator::removeSolubleFraction(bool rescale)
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

  Eigen::ArrayXb activeSprings =
    forceRelaxer.findActiveSprings(&forceRelaxationNetwork);
  Eigen::ArrayXb activeNodes =
    Eigen::ArrayXb::Zero(forceRelaxationNetwork.nrOfNodes);

  RUNTIME_EXP_IFN(activeSprings.size() ==
                    this->simplifiedUniverse.strandFrom.size(),
                  "Number of springs does not match the number of strands. "
                  "Therefore, the mapping would be incorrect.");

  // need to follow connections to mark the whole clusters as active
  bool hasChanged = true;
  while (hasChanged) {
    hasChanged = false;

    for (size_t i = 0; i < activeNodes.size(); ++i) {
      // had already marked/handled this node and its connections
      if (activeNodes[i]) {
        continue;
      }

      bool nodeIsConnectedToActive = false;
      // check if this node is connected to an active spring
      for (int springIdx : forceRelaxationNetwork.springIndicesOfLinks[i]) {
        if (activeSprings[springIdx]) {
          nodeIsConnectedToActive = true;
          break; // no need to check further springs
        }
      }
      if (nodeIsConnectedToActive) {
        activeNodes[i] = true;
        hasChanged = true;
        for (int springIdx : forceRelaxationNetwork.springIndicesOfLinks[i]) {
          activeSprings[springIdx] = true;
        }
      }
    }
  }

  // now that we know which springs and nodes are soluble,
  // we can remove them
  for (long int i = activeSprings.size() - 1; i >= 0; --i) {
    if (!activeSprings[i]) {
      this->removeStrand(i);
    }
  }
  // then, remove the w_sol crosslinks
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
    this->box = pylimer_tools::entities::Box(scalingFactor * this->box.getLx(),
                                             scalingFactor * this->box.getLy(),
                                             scalingFactor * this->box.getLz());
    for (size_t i = 0; i < this->simplifiedUniverse.xlinkTypes.size(); ++i) {
      this->simplifiedUniverse.xlinkX[i] *= scalingFactor;
      this->simplifiedUniverse.xlinkY[i] *= scalingFactor;
      this->simplifiedUniverse.xlinkZ[i] *= scalingFactor;
    }
  }
  this->updateNeighbourListCoordinates();
}
pylimer_tools::sim::mehp::Network
MCUniverseGenerator::convertToForceRelaxationNetwork() const
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
  // crosslinks first, in order to keep consistent numbering
  for (size_t crosslinkIdx = 0; crosslinkIdx < nCrosslinks; ++crosslinkIdx) {
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
    // collapse dangling ends already to connected crosslink
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
  forceRelaxationNetwork.springIndexA = Eigen::ArrayXi::Zero(nSpringEstimate);
  forceRelaxationNetwork.springIndexB = Eigen::ArrayXi::Zero(nSpringEstimate);
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
        // dangling strand
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
      forceRelaxationNetwork.springIndicesOfLinks[from].push_back(newSpringIdx);
      if (to != from) {
        forceRelaxationNetwork.springIndicesOfLinks[to].push_back(newSpringIdx);
      }
      newSpringIdx += 1;
    } else {
      assert(to < 0);
    }
  }
  assert(handledDanglingIdx == nDanglingEnds);
  const size_t nActualSprings = newSpringIdx;
  forceRelaxationNetwork.springIndexA.conservativeResize(nActualSprings);
  forceRelaxationNetwork.springIndexB.conservativeResize(nActualSprings);
  forceRelaxationNetwork.springsContourLength.conservativeResize(
    nActualSprings);
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

  // compute the box offset for the springs
  const Eigen::VectorXd actualDistance =
    forceRelaxationNetwork.coordinates(
      forceRelaxationNetwork.springCoordinateIndexB) -
    forceRelaxationNetwork.coordinates(
      forceRelaxationNetwork.springCoordinateIndexA);
  Eigen::VectorXd expectedDistance = actualDistance;
  this->box.handlePBC(expectedDistance);
  forceRelaxationNetwork.springBoxOffset = expectedDistance - actualDistance;

  forceRelaxationNetwork.meanSpringContourLength =
    forceRelaxationNetwork.springsContourLength.size() > 0
      ? forceRelaxationNetwork.springsContourLength.mean()
      : 1.;
  forceRelaxationNetwork.assumeBoxLargeEnough = true;

  return forceRelaxationNetwork;
}
pylimer_tools::sim::mehp::MEHPForceRelaxation
MCUniverseGenerator::getForceRelaxation() const
{
  // first, convert to a useable structure for the force relaxation
  pylimer_tools::sim::mehp::Network forceRelaxationNetwork =
    this->convertToForceRelaxationNetwork();

  // actually start force relaxation
  pylimer_tools::sim::mehp::MEHPForceRelaxation forceRelaxer =
    pylimer_tools::sim::mehp::MEHPForceRelaxation(forceRelaxationNetwork);

  return forceRelaxer;
}
pylimer_tools::sim::mehp::ForceBalanceNetwork
MCUniverseGenerator::convertToForceBalanceNetwork() const
{
  pylimer_tools::sim::mehp::ForceBalanceNetwork forceBalanceNetwork;
  forceBalanceNetwork.L[0] = this->box.getLx();
  forceBalanceNetwork.L[1] = this->box.getLy();
  forceBalanceNetwork.L[2] = this->box.getLz();
  for (size_t dir = 0; dir < 3; ++dir) {
    forceBalanceNetwork.boxHalfs[dir] = forceBalanceNetwork.L[dir] * 0.5;
  }
  forceBalanceNetwork.vol = this->box.getVolume();

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
  forceBalanceNetwork.nrOfNodes = nCrosslinks + nDanglingEnds;
  forceBalanceNetwork.nrOfLinks = forceBalanceNetwork.nrOfNodes;
  forceBalanceNetwork.oldAtomIds =
    Eigen::ArrayXi::Zero(forceBalanceNetwork.nrOfNodes);
  forceBalanceNetwork.coordinates =
    Eigen::VectorXd(forceBalanceNetwork.nrOfNodes * 3);
  forceBalanceNetwork.linkIsSliplink =
    Eigen::ArrayXb::Zero(forceBalanceNetwork.nrOfLinks);
  forceBalanceNetwork.nrOfCrosslinkSwapsEndured = Eigen::ArrayXi::Zero(0);
  forceBalanceNetwork.springIndicesOfLinks.reserve(nCrosslinks);
  forceBalanceNetwork.oldAtomTypes =
    Eigen::ArrayXi::Zero(forceBalanceNetwork.nrOfNodes);
  // crosslinks first, in order to keep consistent numbering
  for (size_t crosslinkIdx = 0; crosslinkIdx < nCrosslinks; ++crosslinkIdx) {
    forceBalanceNetwork.coordinates(3 * crosslinkIdx + 0) =
      this->simplifiedUniverse.xlinkX[crosslinkIdx];
    forceBalanceNetwork.coordinates(3 * crosslinkIdx + 1) =
      this->simplifiedUniverse.xlinkY[crosslinkIdx];
    forceBalanceNetwork.coordinates(3 * crosslinkIdx + 2) =
      this->simplifiedUniverse.xlinkZ[crosslinkIdx];
    forceBalanceNetwork.oldAtomIds(crosslinkIdx) = crosslinkIdx;
    forceBalanceNetwork.oldAtomTypes(crosslinkIdx) =
      this->simplifiedUniverse.xlinkTypes[crosslinkIdx];
  }
  for (size_t danglingEndIdx = 0; danglingEndIdx < nDanglingEnds;
       ++danglingEndIdx) {
    size_t nodeIdx = nCrosslinks + danglingEndIdx;
    // collapse dangling ends already to connected crosslink
    forceBalanceNetwork.coordinates(3 * nodeIdx + 0) =
      this->simplifiedUniverse.xlinkX[danglingEndPartners[danglingEndIdx]];
    forceBalanceNetwork.coordinates(3 * nodeIdx + 1) =
      this->simplifiedUniverse.xlinkY[danglingEndPartners[danglingEndIdx]];
    forceBalanceNetwork.coordinates(3 * nodeIdx + 2) =
      this->simplifiedUniverse.xlinkZ[danglingEndPartners[danglingEndIdx]];
  }
  for (size_t i = 0; i < forceBalanceNetwork.nrOfNodes; ++i) {
    std::vector<size_t> empty = {};
    forceBalanceNetwork.springIndicesOfLinks.push_back(empty);
  }

  size_t nSpringEstimate = this->simplifiedUniverse.strandFrom.size();
  forceBalanceNetwork.springIndexA = Eigen::ArrayXi::Zero(nSpringEstimate);
  forceBalanceNetwork.springIndexB = Eigen::ArrayXi::Zero(nSpringEstimate);
  forceBalanceNetwork.springsContourLength =
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
        // dangling strand
        to = nCrosslinks + handledDanglingIdx;
        handledDanglingIdx += 1;
        forceBalanceNetwork.springsContourLength(newSpringIdx) =
          this->simplifiedUniverse.beadsInStrand[springIdx];
      } else {
        forceBalanceNetwork.springsContourLength(newSpringIdx) =
          this->simplifiedUniverse.beadsInStrand[springIdx] + 1;
      }
      forceBalanceNetwork.springIndexA(newSpringIdx) = from;
      forceBalanceNetwork.springIndexB(newSpringIdx) = to;
      forceBalanceNetwork.springIndicesOfLinks[from].push_back(newSpringIdx);
      if (from != to) {
        forceBalanceNetwork.springIndicesOfLinks[to].push_back(newSpringIdx);
      }
      newSpringIdx += 1;
    } else {
      assert(to < 0);
    }
  }
  assert(handledDanglingIdx == nDanglingEnds);
  const size_t nActualSprings = newSpringIdx;
  forceBalanceNetwork.springIndexA.conservativeResize(nActualSprings);
  forceBalanceNetwork.springIndexB.conservativeResize(nActualSprings);
  forceBalanceNetwork.springsContourLength.conservativeResize(nActualSprings);
  forceBalanceNetwork.nrOfSprings = nActualSprings;
  forceBalanceNetwork.springsType = Eigen::ArrayXi::Ones(nActualSprings);
  forceBalanceNetwork.springIsActive = Eigen::ArrayXb::Ones(nActualSprings);
  forceBalanceNetwork.partialSpringIsPartial =
    Eigen::ArrayXb::Zero(nActualSprings);
  forceBalanceNetwork.nrOfPartialSprings = nActualSprings;
  forceBalanceNetwork.nrOfSpringsWithPartition = 0;

  forceBalanceNetwork.springCoordinateIndexA =
    Eigen::ArrayXi::Zero(3 * nActualSprings);
  forceBalanceNetwork.springCoordinateIndexB =
    Eigen::ArrayXi::Zero(3 * nActualSprings);
  forceBalanceNetwork.partialToFullSpringIndex =
    Eigen::ArrayXi::LinSpaced(nActualSprings, 0, nActualSprings - 1);
  for (size_t i = 0; i < nActualSprings; ++i) {
    for (size_t dir = 0; dir < 3; ++dir) {
      forceBalanceNetwork.springCoordinateIndexA(3 * i + dir) =
        forceBalanceNetwork.springIndexA(i) * 3 + dir;
      forceBalanceNetwork.springCoordinateIndexB(3 * i + dir) =
        forceBalanceNetwork.springIndexB(i) * 3 + dir;
    }
    std::vector<size_t> linkIndices =
      // forceBalanceNetwork.springIndexA(i) ==
      //     forceBalanceNetwork.springIndexB(i)
      //   ? { { forceBalanceNetwork.springIndexA(i) } }
      //   :
      { static_cast<size_t>(forceBalanceNetwork.springIndexA(i)),
        static_cast<size_t>(forceBalanceNetwork.springIndexB(i)) };
    forceBalanceNetwork.linkIndicesOfSprings.push_back(linkIndices);
    forceBalanceNetwork.localToGlobalSpringIndex.push_back({ i });
  }
  forceBalanceNetwork.springPartBoxOffset =
    this->box.getOffset(forceBalanceNetwork.coordinates(
                          forceBalanceNetwork.springCoordinateIndexB) -
                        forceBalanceNetwork.coordinates(
                          forceBalanceNetwork.springCoordinateIndexA));

  forceBalanceNetwork.springPartCoordinateIndexA =
    forceBalanceNetwork.springCoordinateIndexA;
  forceBalanceNetwork.springPartCoordinateIndexB =
    forceBalanceNetwork.springCoordinateIndexB;
  forceBalanceNetwork.springPartIndexA = forceBalanceNetwork.springIndexA;
  forceBalanceNetwork.springPartIndexB = forceBalanceNetwork.springIndexB;

  forceBalanceNetwork.meanSpringContourLength =
    forceBalanceNetwork.springsContourLength.size() > 0
      ? forceBalanceNetwork.springsContourLength.mean()
      : 1.;

  return forceBalanceNetwork;
}
pylimer_tools::sim::mehp::MEHPForceBalance
MCUniverseGenerator::getForceBalance() const
{
  const pylimer_tools::sim::mehp::ForceBalanceNetwork net =
    this->convertToForceBalanceNetwork();

  Eigen::VectorXd springPartitions = Eigen::VectorXd::Ones(net.nrOfSprings);
  pylimer_tools::sim::mehp::MEHPForceBalance balance =
    pylimer_tools::sim::mehp::MEHPForceBalance(net, springPartitions);
  assert(balance.validateNetwork());
  return balance;
}
pylimer_tools::sim::mehp::MEHPForceBalance2
MCUniverseGenerator::getForceBalance2() const
{
  const pylimer_tools::sim::mehp::ForceBalanceNetwork net =
    this->convertToForceBalanceNetwork();

  Eigen::VectorXd springPartitions = Eigen::VectorXd::Ones(net.nrOfSprings);
  pylimer_tools::sim::mehp::MEHPForceBalance2 balance =
    pylimer_tools::sim::mehp::MEHPForceBalance2(net, springPartitions);
  assert(balance.validateNetwork());
  return balance;
}
void
MCUniverseGenerator::relaxCrosslinks()
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
  // dangling strands lead to higher nr of nodes
  RUNTIME_EXP_IFN(forceRelaxationNetwork.nrOfNodes >= nCrosslinks,
                  "Expected force relaxation to preserve crosslinks.");
  for (size_t i = 0; i < nCrosslinks; ++i) {
    this->simplifiedUniverse.xlinkX[i] =
      forceRelaxationNetwork.coordinates(3 * i + 0);
    this->simplifiedUniverse.xlinkY[i] =
      forceRelaxationNetwork.coordinates(3 * i + 1);
    this->simplifiedUniverse.xlinkZ[i] =
      forceRelaxationNetwork.coordinates(3 * i + 2);
  }
  this->updateNeighbourListCoordinates();
}
size_t
MCUniverseGenerator::getCurrentNrOfAtoms() const
{
  size_t nAtomsTotal =
    std::reduce(this->simplifiedUniverse.beadsInStrand.begin(),
                this->simplifiedUniverse.beadsInStrand.end(),
                0) +
    this->simplifiedUniverse.xlinkTypes.size();
  return nAtomsTotal;
}
size_t
MCUniverseGenerator::getCurrentNrOfBonds() const
{
  size_t nBondsTotal = 0;
  for (size_t strandIdx = 0;
       strandIdx < this->simplifiedUniverse.strandFrom.size();
       ++strandIdx) {
    size_t nBeads = this->simplifiedUniverse.beadsInStrand[strandIdx];
    long int from = this->simplifiedUniverse.strandFrom[strandIdx];
    long int to = this->simplifiedUniverse.strandTo[strandIdx];
    if (nBeads > 0) {
      nBondsTotal += nBeads - 1;
      if (to >= 0) {
        nBondsTotal += 1;
      }
    }
    if (from >= 0) {
      nBondsTotal += 1;
    }
  }
  return nBondsTotal;
}
double
MCUniverseGenerator::getCurrentCrosslinkerConversion() const
{
  return 1. - (static_cast<double>(this->nrOfAvailableCrosslinkSites) /
               static_cast<double>(this->originalNrOfAvailableCrosslinkSites));
}
void
MCUniverseGenerator::validateInternalState() const
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
    all_equal<size_t>(6,
                      this->simplifiedUniverse.xlinkTypes.size(),
                      this->simplifiedUniverse.xlinkX.size(),
                      this->simplifiedUniverse.xlinkY.size(),
                      this->simplifiedUniverse.xlinkZ.size(),
                      this->remainingCrossLinkerFunctionality.size(),
                      this->simplifiedUniverse.xlinkChainId.size()),
    "Inconsistent sizes in simplified universe.");

  long int nCrosslinks = this->simplifiedUniverse.xlinkX.size();
  const long int nrOfAvailableSites =
    std::reduce(this->remainingCrossLinkerFunctionality.begin(),
                this->remainingCrossLinkerFunctionality.end(),
                0);

  RUNTIME_EXP_IFN(this->nrOfAvailableCrosslinkSites == nrOfAvailableSites,
                  "Inconsistent nr of crosslink sites.");

  for (size_t xlinkIdx = 0; xlinkIdx < nCrosslinks; ++xlinkIdx) {
    for (long int subStrandIdx :
         this->simplifiedUniverse.strandsOfXlink[xlinkIdx]) {
      RUNTIME_EXP_IFN(
        this->simplifiedUniverse.strandFrom[subStrandIdx] == xlinkIdx ||
          this->simplifiedUniverse.strandTo[subStrandIdx] == xlinkIdx,
        "Inconsistent links int list of crosslinks <> strands.");
    }
    if (xlinkIdx > 0) {
      RUNTIME_EXP_IFN(
        this->simplifiedUniverse.xlinkChainId[xlinkIdx] >=
          this->simplifiedUniverse.xlinkChainId[xlinkIdx - 1],
        "Expected the chain ID to increase with crosslink index.");
    }
  }
  for (size_t strandIdx = 0;
       strandIdx < this->simplifiedUniverse.strandFrom.size();
       ++strandIdx) {
    RUNTIME_EXP_IFN(this->simplifiedUniverse.beadsInStrand[strandIdx] >= 0,
                    "Expected a positive number of beads per strand.");
    RUNTIME_EXP_IFN(this->simplifiedUniverse.strandFrom[strandIdx] <
                      nCrosslinks,
                    "Expected a valid crosslink index as a strand origin.");
    RUNTIME_EXP_IFN(
      this->simplifiedUniverse.strandTo[strandIdx] < nCrosslinks,
      "Expected a valid crosslink index as a strand destination.");

    if (this->simplifiedUniverse.beadsInStrand[strandIdx] == 0) {
      RUNTIME_EXP_IFN(this->simplifiedUniverse.strandFrom[strandIdx] >= 0 &&
                        this->simplifiedUniverse.strandTo[strandIdx] >= 0,
                      "Expected a valid crosslink index as a strand "
                      "origin and destination for empty length strands.");
    }
  }
}
Eigen::VectorXd
MCUniverseGenerator::sampleFreeChainCoordinates(int chainLen)
{
  Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
    chainLen, this->beadDistance, this->meanSquaredBeadDistance, this->rng);

  if (this->nMcSteps > 0) {
    pylimer_tools::sim::equilibrateChainWithMC(positions,
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
Eigen::VectorXd
MCUniverseGenerator::sampleDanglingChainCoordinates(size_t idxFrom,
                                                    int chainLen)
{
  INVALIDARG_EXP_IFN(idxFrom < this->simplifiedUniverse.xlinkTypes.size(),
                     "Invalid index for dangling chain start");
  Eigen::VectorXd positions = pylimer_tools::utils::doRandomWalkChain(
    chainLen, this->beadDistance, this->meanSquaredBeadDistance, this->rng);

  if (this->nMcSteps > 0) {
    pylimer_tools::sim::equilibrateChainWithMC(positions,
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
Eigen::VectorXd
MCUniverseGenerator::sampleStrandCoordinates(size_t from,
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
    pylimer_tools::sim::equilibrateChainWithMC(positions,
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
void
MCUniverseGenerator::linkStrandToCrosslink(size_t strandIdx,
                                           size_t crosslinkIdx,
                                           bool ignoreInexistent)
{
  INVALIDARG_EXP_IFN(
    strandIdx < this->simplifiedUniverse.strandFrom.size(),
    "The strand index " + std::to_string(strandIdx) +
      " is out of bounds, only " +
      std::to_string(this->simplifiedUniverse.strandFrom.size()) +
      " strands have been registered.");
  INVALIDARG_EXP_IFN(
    this->simplifiedUniverse.strandFrom[strandIdx] == UNCONNECTED ||
      this->simplifiedUniverse.strandTo[strandIdx] == UNCONNECTED,
    "Require one end to be free to be connected to");
  if (!ignoreInexistent) {
    INVALIDARG_EXP_IFN(crosslinkIdx <
                         this->remainingCrossLinkerFunctionality.size(),
                       "The chosen crosslink " + std::to_string(crosslinkIdx) +
                         " does not exist.");
    INVALIDARG_EXP_IFN(this->remainingCrossLinkerFunctionality[crosslinkIdx] >
                         0,
                       "The crosslink " + std::to_string(crosslinkIdx) +
                         " is already fully linked.");
    this->remainingCrossLinkerFunctionality[crosslinkIdx] -= 1;
    this->nrOfAvailableCrosslinkSites -= 1;
    this->simplifiedUniverse.strandsOfXlink[crosslinkIdx].push_back(strandIdx);
  }

  if (this->simplifiedUniverse.strandFrom[strandIdx] == UNCONNECTED) {
    this->simplifiedUniverse.strandFrom[strandIdx] = crosslinkIdx;
  } else {
    INVALIDARG_EXP_IFN(
      this->simplifiedUniverse.strandFrom[strandIdx] >= 0,
      "Strands are required to have the first end be connected first.");
    this->simplifiedUniverse.strandTo[strandIdx] = crosslinkIdx;
    Eigen::VectorXd newCoordinates = this->getCrosslinkCoordinates();
  }
}
std::vector<size_t>
MCUniverseGenerator::addXlinkAtoms(int nrOfAtomsToAdd,
                                   int atomType,
                                   Eigen::VectorXd coordinates)
{
  INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                     "Coordinates must have a size multiple of 3");
  INVALIDARG_EXP_IFN(coordinates.size() / 3 == nrOfAtomsToAdd,
                     "Coordinates must match the promised number of atoms");
  size_t currentNrOfJunctions = this->simplifiedUniverse.xlinkX.size();
  this->simplifiedUniverse.xlinkChainId.reserve(currentNrOfJunctions +
                                                nrOfAtomsToAdd);
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

    this->simplifiedUniverse.xlinkChainId.push_back(currentNrOfJunctions + i);
    this->simplifiedUniverse.strandsOfXlink.push_back({});
  }

  return indicesAdded;
}
std::vector<size_t>
MCUniverseGenerator::addXlinkAtoms(int nrOfAtomsToAdd,
                                   int atomType,
                                   bool whiteNoise)
{
  Eigen::VectorXd randomPos =
    this->generateRandomPositions(nrOfAtomsToAdd, whiteNoise);
  return this->addXlinkAtoms(nrOfAtomsToAdd, atomType, randomPos);
}
Eigen::VectorXd
MCUniverseGenerator::generateRandomWhitePositions(int nSamples)
{
  Eigen::VectorXd coordinates = Eigen::VectorXd(3 * nSamples);

  for (size_t i = 0; i < nSamples; ++i) {
    coordinates(3 * i) = this->distX(this->rng);
    coordinates(3 * i + 1) = this->distY(this->rng);
    coordinates(3 * i + 2) = this->distZ(this->rng);
  }

  return coordinates;
}
Eigen::VectorXd
MCUniverseGenerator::generateRandomBluePositions(int nSamples)
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
long int
MCUniverseGenerator::findAppropriateLink(size_t from,
                                         const double desiredR02,
                                         const double maxDistance)
{
  RUNTIME_EXP_IFN(this->simplifiedUniverse.xlinkTypes.size() ==
                    this->remainingCrossLinkerFunctionality.size(),
                  "Invalid internal state, number of and remaining "
                  "crosslink functionalities mismatch");

  double sumOfWeights = 0.0;
  const double normalisationFactorInExponential = -3. / (2. * desiredR02);

  std::vector<size_t> suitableMatches;
  std::vector<double> matchWeights;
  if (maxDistance < 0.) {
    matchWeights.reserve(this->nrOfAvailableCrosslinkSites);
    suitableMatches.reserve(this->nrOfAvailableCrosslinkSites);
    size_t nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
    for (int i = 0; i < nCrosslinks; ++i) {
      double thisWeight = this->evaluatePartnerProbability(
        from, i, normalisationFactorInExponential, maxDistance);
      if (thisWeight <= 0.0) {
        continue;
      }
      suitableMatches.push_back(i);
      matchWeights.push_back(thisWeight);
      sumOfWeights += thisWeight;
    }
  } else {
    // else: use a neighbour list instead, to improve performance
    // we could possibly improve this further by pre-allocating the result
    // array
    Eigen::ArrayXi acceptableIndices =
      this->xlinkNeighbourList.getIndicesCloseToCoordinates(
        Eigen::Vector3d(this->simplifiedUniverse.xlinkX[from],
                        this->simplifiedUniverse.xlinkY[from],
                        this->simplifiedUniverse.xlinkZ[from]),
        maxDistance);
    suitableMatches.reserve(acceptableIndices.size());
    matchWeights.reserve(acceptableIndices.size());
    for (int i = 0; i < acceptableIndices.size(); ++i) {
      int crosslinkIdx = acceptableIndices(i);
      double thisWeight = this->evaluatePartnerProbability(
        from, crosslinkIdx, normalisationFactorInExponential, maxDistance);
      if (thisWeight <= 0.0) {
        continue;
      }
      suitableMatches.push_back(crosslinkIdx);
      matchWeights.push_back(thisWeight);
      sumOfWeights += thisWeight;
    }
  }

  if (suitableMatches.size() == 0) {
    if (maxDistance > 0.) {
      std::cerr << "No suitable partner found within max distance. "
                   "Increasing for this link."
                << std::endl;
      return this->findAppropriateLink(from, desiredR02, -1.);
    }
    return -1;
  }

  std::discrete_distribution<long int> weightDist(matchWeights.begin(),
                                                  matchWeights.end());
  return suitableMatches[weightDist(this->rng)];
}
double
MCUniverseGenerator::evaluatePartnerProbability(
  size_t from,
  size_t to,
  const double normalisationFactorInExponential,
  const double maxDistance,
  bool respectXlinkChains)
{
  if (this->remainingCrossLinkerFunctionality[to] < 1) {
    return -1.;
  }
  Eigen::Vector3d dist = this->getVectorBetween(from, to);
  double thisWeight = -1.;
  if (dist.norm() < maxDistance || maxDistance < 0.) {
    thisWeight =
      static_cast<double>(this->remainingCrossLinkerFunctionality[to]) *
      std::exp(dist.squaredNorm() * normalisationFactorInExponential);
    int strand1 = this->simplifiedUniverse.xlinkChainId[from];
    int strand2 = this->simplifiedUniverse.xlinkChainId[to];
    if (to == from || (respectXlinkChains && strand1 == strand2)) {
      thisWeight *= this->primaryLoopProbability;
    }
    if (this->secondaryLoopProbability != 1.) {
      // check whether this crosslink would lead to a secondary loop
      // => apply weight for every other already existing back-link

      // secondary loops in this case are any loops between two strands
      // we store the chain ids sequentially -> use to reduce the
      // iterations
      size_t xlinkIdxOnStrand1 =
        respectXlinkChains
          ? pylimer_tools::utils::first_occuring_index(
              this->simplifiedUniverse.xlinkChainId, strand1, from)
          : from;
      assert(xlinkIdxOnStrand1 < this->simplifiedUniverse.xlinkChainId.size());
      while (xlinkIdxOnStrand1 < this->simplifiedUniverse.xlinkChainId.size() &&
             this->simplifiedUniverse.xlinkChainId[xlinkIdxOnStrand1] ==
               strand1) {
        for (size_t partnersSubStrand :
             this->simplifiedUniverse.strandsOfXlink[xlinkIdxOnStrand1]) {
          assert(this->simplifiedUniverse.strandFrom[partnersSubStrand] ==
                   xlinkIdxOnStrand1 ||
                 this->simplifiedUniverse.strandTo[partnersSubStrand] ==
                   xlinkIdxOnStrand1);
          long int otherStrandEnd =
            this->simplifiedUniverse.strandFrom[partnersSubStrand] ==
                xlinkIdxOnStrand1
              ? this->simplifiedUniverse.strandTo[partnersSubStrand]
              : this->simplifiedUniverse.strandFrom[partnersSubStrand];
          bool wouldBeSecondaryLoop = false;
          if (otherStrandEnd >= 0) {
            // more complex logic to consider crosslink chains
            wouldBeSecondaryLoop =
              ((strand2 ==
                this->simplifiedUniverse.xlinkChainId[otherStrandEnd]));
          }
          if (wouldBeSecondaryLoop) {
            thisWeight *= this->secondaryLoopProbability;
          }
        }
        xlinkIdxOnStrand1 += 1;
        if (!respectXlinkChains) {
          break;
        }
      }
    }
  }
  return thisWeight;
}
double
MCUniverseGenerator::getIdealCutoff() const
{
  if (this->simplifiedUniverse.beadsInStrand.size() == 0) {
    return -1.;
  }
  return this->maxDistanceProvider->getMaxDistance(
    static_cast<double>(
      *std::max_element(std::begin(this->simplifiedUniverse.beadsInStrand),
                        std::end(this->simplifiedUniverse.beadsInStrand))) +
    1.);
}
void
MCUniverseGenerator::resetNeighbourList()
{
  double idealCutoff = this->getIdealCutoff();
  if (idealCutoff > 0.) {
    Eigen::VectorXd newCoordinates = this->getCrosslinkCoordinates();
    this->xlinkNeighbourList.initialize(newCoordinates, this->box, idealCutoff);
  }
}
void
MCUniverseGenerator::updateNeighbourListCoordinates()
{
  if (this->getIdealCutoff() > 0.) {
    Eigen::VectorXd newCoordinates = this->getCrosslinkCoordinates();
    this->xlinkNeighbourList.resetCoordinates(newCoordinates);
  }
}
double
MCUniverseGenerator::distanceBetween(size_t i, size_t j) const
{
  return this->getDistance(this->simplifiedUniverse.xlinkX[i],
                           this->simplifiedUniverse.xlinkY[i],
                           this->simplifiedUniverse.xlinkZ[i],
                           this->simplifiedUniverse.xlinkX[j],
                           this->simplifiedUniverse.xlinkY[j],
                           this->simplifiedUniverse.xlinkZ[j]);
}
Eigen::Vector3d
MCUniverseGenerator::getVectorBetween(size_t i, size_t j) const
{
  Eigen::Vector3d diff;
  diff << (this->simplifiedUniverse.xlinkX[j] -
           this->simplifiedUniverse.xlinkX[i]),
    (this->simplifiedUniverse.xlinkY[j] - this->simplifiedUniverse.xlinkY[i]),
    (this->simplifiedUniverse.xlinkZ[j] - this->simplifiedUniverse.xlinkZ[i]);
  this->box.handlePBC(diff);
  return diff;
}
double
MCUniverseGenerator::getDistance(double x1,
                                 double y1,
                                 double z1,
                                 double x2,
                                 double y2,
                                 double z2) const
{
  Eigen::Vector3d diff;
  diff << x2 - x1, y2 - y1, z2 - z1;
  this->box.handlePBC(diff);
  return diff.norm();
}
Eigen::Vector3d
MCUniverseGenerator::sampleCoordinatesWithinNBeadDistance(int nBeads)
{
  double stddev = std::sqrt(static_cast<double>(nBeads + 1) *
                            this->meanSquaredBeadDistance / 3.);
  std::normal_distribution<double> otherEndCoordinateDist =
    std::normal_distribution<double>(0., stddev);
  // sample the coordinates delta, return
  Eigen::Vector3d otherEndCoordinates;
  otherEndCoordinates << otherEndCoordinateDist(this->rng),
    otherEndCoordinateDist(this->rng), otherEndCoordinateDist(this->rng);
  return otherEndCoordinates;
}
Eigen::VectorXd
MCUniverseGenerator::getCrosslinkCoordinates() const
{
  size_t nCrosslinks = this->simplifiedUniverse.xlinkTypes.size();
  Eigen::VectorXd coordinates = Eigen::VectorXd(3 * nCrosslinks);
  for (size_t i = 0; i < nCrosslinks; ++i) {
    coordinates(3 * i) = this->simplifiedUniverse.xlinkX[i];
    coordinates(3 * i + 1) = this->simplifiedUniverse.xlinkY[i];
    coordinates(3 * i + 2) = this->simplifiedUniverse.xlinkZ[i];
  }
  return coordinates;
}
void
MCUniverseGenerator::linkStrandTo(size_t strandIdx, size_t crosslinkIdx)
{
  return this->linkStrandToCrosslink(strandIdx, crosslinkIdx, false);
}
}
