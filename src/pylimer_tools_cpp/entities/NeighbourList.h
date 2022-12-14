#ifndef NEIGHBOURLIST_H
#define NEIGHBOURLIST_H

extern "C"
{
#include <igraph/igraph.h>
}
#include "../utils/VectorUtils.h"
#include "Atom.h"
#include "Box.h"
#include <Eigen/Dense>
#include <map>
#include <unordered_map>
#include <vector>

#include <algorithm>

namespace pylimer_tools {
namespace entities {

  class NeighbourList
  {
  public:
    NeighbourList(const std::vector<pylimer_tools::entities::Atom> atoms,
                  const pylimer_tools::entities::Box box,
                  double cutoff)
    {
      if (cutoff <= 1e-3) {
        throw std::invalid_argument("Cutoff must be larger than zero");
      }

      this->cutoff = cutoff;
      this->atoms = atoms;
      this->box = box;

      this->nrOfBucketsX =
        static_cast<size_t>(std::floor(box.getLx() / cutoff));
      this->nrOfBucketsY =
        static_cast<size_t>(std::floor(box.getLy() / cutoff));
      this->nrOfBucketsZ =
        static_cast<size_t>(std::floor(box.getLz() / cutoff));

      this->bucketWidthX =
        box.getLx() / static_cast<double>(this->nrOfBucketsX);
      this->bucketWidthY =
        box.getLy() / static_cast<double>(this->nrOfBucketsY);
      this->bucketWidthZ =
        box.getLz() / static_cast<double>(this->nrOfBucketsZ);

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
        size_t bucketIndex = this->getBucketIndexForTriplet(
          this->getBucketIndicesForAtom(atoms[i]));
        this->neighbourBuckets[bucketIndex].push_back(i);
      }
    };

    std::vector<pylimer_tools::entities::Atom> getAtomsCloseTo(
      pylimer_tools::entities::Atom atom)
    {
      return this->getAtomsCloseTo(atom, this->cutoff);
    }

    std::vector<pylimer_tools::entities::Atom> getAtomsCloseTo(
      pylimer_tools::entities::Atom atom,
      double newCutoff)
    {
      if (newCutoff < 0.) {
        newCutoff = this->cutoff;
      }

      size_t indexBasis =
        this->getBucketIndexForTriplet(this->getBucketIndicesForAtom(atom));
      bool foundSelf = false;
      for (size_t atomIndex : this->neighbourBuckets[indexBasis]) {
        if (this->atoms[atomIndex].getId() == atom.getId()) {
          foundSelf = true;
        }
      }
      if (!foundSelf) {
        throw std::invalid_argument("The requested atom is not in the list");
      }

      std::vector<size_t> bucketIndices =
        this->getCombinedBucketIndicesForAtom(atom, newCutoff);
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
        std::vector<size_t> atomIndices = this->neighbourBuckets[bucketIndex];
        for (size_t atomIndex : atomIndices) {
          if (this->atoms[atomIndex].distanceTo(atom, &this->box) < newCutoff &&
              this->atoms[atomIndex].getId() != atom.getId()) {
            results.push_back(this->atoms[atomIndex]);
          }
        }
      }
      if (!foundBasis) {
        throw std::runtime_error(
          "Did not find basis bucket. Something is wrong.");
      }
      // return results
      return results;
    }

    void removeAtom(Atom atom)
    {
      size_t indexBasis =
        this->getBucketIndexForTriplet(this->getBucketIndicesForAtom(atom));
      // it is sufficient to remove this atom from this bucket
      if (this->idToAtomIdx.size() == 0) {
        this->idToAtomIdx.reserve(this->atoms.size());
        for (size_t i = 0; i < this->atoms.size(); ++i) {
          this->idToAtomIdx.emplace(this->atoms[i].getId(), i);
        }
      }
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
        throw std::invalid_argument("This atom is not in a bucket anyway");
      }
    }

  protected:
    size_t normalizeBucketIndex(long int bucketIndex, size_t nrOfBuckets)
    {
      while (bucketIndex < 0) {
        bucketIndex = bucketIndex + nrOfBuckets;
      }
      while (bucketIndex >= nrOfBuckets) {
        bucketIndex = bucketIndex - nrOfBuckets;
      }
      return static_cast<size_t>(bucketIndex);
    }

    size_t getBucketIndexForTriplet(
      std::tuple<long int, long int, long int> ind)
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

    std::tuple<size_t, size_t, size_t> getBucketIndicesForAtom(
      const pylimer_tools::entities::Atom& atom)
    {
      return std::make_tuple(
        static_cast<size_t>(std::floor(atom.getX() / this->bucketWidthX)),
        static_cast<size_t>(std::floor(atom.getY() / this->bucketWidthY)),
        static_cast<size_t>(std::floor(atom.getZ() / this->bucketWidthZ)));
    }

    std::vector<size_t> getCombinedBucketIndicesForAtom(
      const pylimer_tools::entities::Atom& atom,
      double newCutoff)
    {
      std::vector<size_t> result = std::vector<size_t>();
      std::tuple<long int, long int, long int> indexBasis = this->getBucketIndicesForAtom(atom);
      
      int nrOfBucketsPerSide =
        (newCutoff <= this->cutoff) ? 3 : std::ceil(3 * newCutoff / this->cutoff);
      if (nrOfBucketsPerSide % 2 == 0) {
        nrOfBucketsPerSide += 1;
      }
      int nrOfBucketsPerQuarter = (nrOfBucketsPerSide - 1) / 2;

      result.reserve(nrOfBucketsPerQuarter * nrOfBucketsPerQuarter *
                     nrOfBucketsPerQuarter);
      for (int offsetX = -nrOfBucketsPerQuarter;
           offsetX <= nrOfBucketsPerQuarter;
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

  private:
    double bucketWidthX;
    double bucketWidthY;
    double bucketWidthZ;

    size_t nrOfBucketsX;
    size_t nrOfBucketsY;
    size_t nrOfBucketsZ;
    size_t totalNrOfBuckets;

    double cutoff;

    pylimer_tools::entities::Box box;

    std::unordered_map<size_t, std::vector<size_t>> neighbourBuckets;

    std::vector<pylimer_tools::entities::Atom> atoms;
    std::unordered_map<size_t, size_t> idToAtomIdx;
  };
};
}

#endif
