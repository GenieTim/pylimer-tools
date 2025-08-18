extern "C"
{
#include <igraph/igraph.h>
}
#include "../utils/VectorUtils.h"
#include "Atom.h"
#include "Box.h"
#include "NeighbourList.h"
#include <cmath>
#include <unordered_map>
#include <vector>

#include <algorithm>

namespace pylimer_tools {
namespace entities {

  NeighbourList::NeighbourList(const std::vector<Atom>& atomList,
                               const Box& boxGeometry,
                               const double cutoffDistance)
  {
    if (cutoffDistance <= 1e-3) {
      throw std::invalid_argument("Cutoff must be larger than zero");
    }

    this->cutoff = cutoffDistance;
    this->atoms = atomList;
    this->box = boxGeometry;

    this->nrOfBucketsX = std::max(
      static_cast<size_t>(1),
      static_cast<size_t>(std::floor(boxGeometry.getLx() / cutoffDistance)));
    this->nrOfBucketsY = std::max(
      static_cast<size_t>(1),
      static_cast<size_t>(std::floor(boxGeometry.getLy() / cutoffDistance)));
    this->nrOfBucketsZ = std::max(
      static_cast<size_t>(1),
      static_cast<size_t>(std::floor(boxGeometry.getLz() / cutoffDistance)));

    this->bucketWidthX =
      boxGeometry.getLx() / static_cast<double>(this->nrOfBucketsX);
    this->bucketWidthY =
      boxGeometry.getLy() / static_cast<double>(this->nrOfBucketsY);
    this->bucketWidthZ =
      boxGeometry.getLz() / static_cast<double>(this->nrOfBucketsZ);

    this->totalNrOfBuckets =
      this->nrOfBucketsX * this->nrOfBucketsY * this->nrOfBucketsZ;

    this->neighbourBuckets.reserve(this->totalNrOfBuckets);

    // prepare the buckets
    for (size_t bucketIndex = 0; bucketIndex < this->totalNrOfBuckets;
         ++bucketIndex) {
      std::vector<size_t> vectorToPlace = std::vector<size_t>();
      // reserve a sensible capacity as estimated
      vectorToPlace.reserve(atoms.size() / this->totalNrOfBuckets);
      this->neighbourBuckets.emplace(bucketIndex, vectorToPlace);
    }

    // fill neighbour buckets
    for (size_t i = 0; i < atoms.size(); ++i) {
      size_t bucketIndex =
        this->getBucketIndexForTriplet(this->getBucketIndicesForAtom(atoms[i]));
      this->neighbourBuckets[bucketIndex].push_back(i);
    }
  };

  std::vector<pylimer_tools::entities::Atom> NeighbourList::getAtomsCloseTo(
    const pylimer_tools::entities::Atom& atom)
  {
    return this->getAtomsCloseTo(atom, this->cutoff);
  }

  std::vector<pylimer_tools::entities::Atom> NeighbourList::getAtomsCloseTo(
    const pylimer_tools::entities::Atom& atom,
    double upperCutoff,
    const double lowerCutoff,
    const bool unwrapped,
    const bool expectSelf)
  {
    if (lowerCutoff > upperCutoff) {
      throw std::invalid_argument("Expected upper cutoff > lower cutoff, got " +
                                  std::to_string(upperCutoff) + " and " +
                                  std::to_string(lowerCutoff) + ".");
    }
    if (upperCutoff <= 0.) {
      upperCutoff = this->cutoff;
      // throw std::invalid_argument("Expected upper cutoff > 0, got " +
      //                             std::to_string(upperCutoff) + ".");
    }
    std::tuple<long int, long int, long int> bucketIndicesForAtom =
      this->getBucketIndicesForAtom(atom);
    size_t indexBasis = this->getBucketIndexForTriplet(bucketIndicesForAtom);
    bool foundSelf = !expectSelf;
    for (size_t atomIndex : this->neighbourBuckets[indexBasis]) {
      if (this->atoms[atomIndex].getId() == atom.getId()) {
        foundSelf = true;
      }
    }
    if (!foundSelf) {
      throw std::invalid_argument("The requested atom is not in the list");
    }

    std::vector<size_t> bucketIndices =
      this->getCombinedBucketIndicesForAtom(atom, upperCutoff);
    std::vector<pylimer_tools::entities::Atom> results =
      std::vector<pylimer_tools::entities::Atom>();
    // good estimate for nr of atoms to return
    results.reserve(bucketIndices.size() * this->atoms.size() /
                    this->totalNrOfBuckets);
    // actually loop the buckets, look for atoms that are close
    bool foundBasis = false;
    for (size_t bucketIndex : bucketIndices) {
      if (bucketIndex == indexBasis) {
        foundBasis = true;
      }
      std::vector<size_t> atomIndices = this->neighbourBuckets.at(bucketIndex);
      for (size_t atomIndex : atomIndices) {
        double distance =
          !unwrapped
            ? this->atoms[atomIndex].distanceTo(atom, this->box)
            : this->atoms[atomIndex].distanceToUnwrapped(atom, this->box);
        if (distance < upperCutoff && distance >= lowerCutoff &&
            this->atoms[atomIndex].getId() != atom.getId()) {
          results.push_back(this->atoms[atomIndex]);
        }
      }
    }
    if (!foundBasis) {
      throw std::runtime_error(
        "Did not find basis bucket. Something is wrong.");
    }
    results.shrink_to_fit();
    // return results
    return results;
  }

  void NeighbourList::removeAtom(const Atom& atom, const std::string debugHint)
  {
    if (this->idToAtomIdx.size() == 0) {
      this->idToAtomIdx.reserve(this->atoms.size());
      for (size_t i = 0; i < this->atoms.size(); ++i) {
        this->idToAtomIdx.emplace(this->atoms[i].getId(), i);
      }
    }
    // it is sufficient to remove this atom from this bucket
    size_t indexBasis =
      this->getBucketIndexForTriplet(this->getBucketIndicesForAtom(atom));
    // have to remove the element with value
    size_t valToRemove = this->idToAtomIdx.at(atom.getId());
    std::vector<size_t>::iterator position =
      std::find(this->neighbourBuckets.at(indexBasis).begin(),
                this->neighbourBuckets.at(indexBasis).end(),
                valToRemove);
    if (position != this->neighbourBuckets.at(indexBasis).end()) {
      this->neighbourBuckets.at(indexBasis).erase(position);
      // std::remove(this->neighbourBuckets.at(indexBasis).begin(),
      //                    this->neighbourBuckets.at(indexBasis).end(),
      //                    valToRemove),
      //        this->neighbourBuckets.at(indexBasis).end());
    } else {
      throw std::invalid_argument(
        "This atom with id " + std::to_string(atom.getId()) +
        " is not in a bucket, so that it cannot be removed." +
        " Debug-hint: " + debugHint + ".");
    }
  }

  /////////////////////////////////////////////////////////////
  // protected
  /////////////////////////////////////////////////////////////

  size_t NeighbourList::normalizeBucketIndex(const long int bucketIndex,
                                             const size_t nrOfBuckets) const
  {
    size_t result =
      bucketIndex - nrOfBuckets * std::floor(static_cast<double>(bucketIndex) /
                                             static_cast<double>(nrOfBuckets));
    assert(result >= 0 && result <= nrOfBuckets);
    return result;
    // while (bucketIndex < 0) {
    //   bucketIndex = bucketIndex + nrOfBuckets;
    // }
    // while (bucketIndex >= nrOfBuckets) {
    //   bucketIndex = bucketIndex - nrOfBuckets;
    // }
    // return static_cast<size_t>(bucketIndex);
  }

  size_t NeighbourList::getBucketIndexForTriplet(
    std::tuple<long int, long int, long int> ind) const
  {
    size_t bucketIndexX =
      this->normalizeBucketIndex(std::get<0>(ind), this->nrOfBucketsX);
    size_t bucketIndexY =
      this->normalizeBucketIndex(std::get<1>(ind), this->nrOfBucketsY);
    size_t bucketIndexZ =
      this->normalizeBucketIndex(std::get<2>(ind), this->nrOfBucketsZ);
    return bucketIndexX + bucketIndexY * this->nrOfBucketsX +
           bucketIndexZ * this->nrOfBucketsX * this->nrOfBucketsY;
  }

  std::tuple<long int, long int, long int>
  NeighbourList::getBucketIndicesForAtom(
    const pylimer_tools::entities::Atom& atom) const
  {
    return std::make_tuple(
      static_cast<long int>(
        std::floor(atom.getUnwrappedX(this->box) / this->bucketWidthX)),
      static_cast<long int>(
        std::floor(atom.getUnwrappedY(this->box) / this->bucketWidthY)),
      static_cast<long int>(
        std::floor(atom.getUnwrappedZ(this->box) / this->bucketWidthZ)));
  }

  std::vector<size_t> NeighbourList::getCombinedBucketIndicesForAtom(
    const pylimer_tools::entities::Atom& atom,
    const double newCutoff) const
  {
    std::vector<size_t> result = std::vector<size_t>();
    std::tuple<long int, long int, long int> indexBasis =
      this->getBucketIndicesForAtom(atom);

    int nrOfBucketsPerSide =
      (newCutoff <= this->cutoff) ? 3 : std::ceil(3 * newCutoff / this->cutoff);
    if (nrOfBucketsPerSide % 2 == 0) {
      nrOfBucketsPerSide += 1;
    }
    int nrOfBucketsPerQuarter = (nrOfBucketsPerSide - 1) / 2;

    result.reserve(nrOfBucketsPerQuarter * nrOfBucketsPerQuarter *
                   nrOfBucketsPerQuarter);
    for (int offsetX = -nrOfBucketsPerQuarter; offsetX <= nrOfBucketsPerQuarter;
         offsetX++) {
      for (int offsetY = -nrOfBucketsPerQuarter;
           offsetY <= nrOfBucketsPerQuarter;
           offsetY++) {
        for (int offsetZ = -nrOfBucketsPerQuarter;
             offsetZ <= nrOfBucketsPerQuarter;
             offsetZ++) {
          size_t newIndex = this->getBucketIndexForTriplet(
            std::make_tuple(std::get<0>(indexBasis) + offsetX,
                            std::get<1>(indexBasis) + offsetY,
                            std::get<2>(indexBasis) + offsetZ));
          result.push_back(newIndex);
          if (offsetX == 0 && offsetY == 0 && offsetZ == 0) {
            assert(newIndex == this->getBucketIndexForTriplet(indexBasis));
          }
        }
      }
    }

    return result;
  }
};
}
