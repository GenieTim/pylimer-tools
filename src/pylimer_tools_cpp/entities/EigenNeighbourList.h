#ifndef EIGENNEIGHBOURLIST_H
#define EIGENNEIGHBOURLIST_H

extern "C"
{
#include <igraph/igraph.h>
}
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include "Box.h"
#include <Eigen/Dense>
#include <map>
#include <unordered_map>
#include <vector>

#include <algorithm>

namespace pylimer_tools {
namespace entities {

  class EigenNeighbourList
  {
  public:
    EigenNeighbourList(const Eigen::VectorXd& coordinates,
                       const Box box,
                       double cutoff)
    {
      INVALIDARG_EXP_IFN(cutoff > 1e-3, "Cutoff must be larger than zero");
      INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                         "Coordinates must be in three");

      this->cutoff = cutoff;
      this->box = box;

      this->nrOfBucketsX =
        static_cast<size_t>(std::floor(box.getLx() / cutoff));
      this->nrOfBucketsY =
        static_cast<size_t>(std::floor(box.getLy() / cutoff));
      this->nrOfBucketsZ =
        static_cast<size_t>(std::floor(box.getLz() / cutoff));

      this->bucketWidths = Eigen::Array3d();
      this->bucketWidths[0] =
        box.getLx() / static_cast<double>(this->nrOfBucketsX);
      this->bucketWidths[1] =
        box.getLy() / static_cast<double>(this->nrOfBucketsY);
      this->bucketWidths[2] =
        box.getLz() / static_cast<double>(this->nrOfBucketsZ);

      this->totalNrOfBuckets =
        this->nrOfBucketsX * this->nrOfBucketsY * this->nrOfBucketsZ;

      this->neighbourBuckets.reserve(this->totalNrOfBuckets);
      this->neighbourBucketNeighboursDefaultCutoff.reserve(
        this->totalNrOfBuckets);

      // prepare the buckets
      for (size_t bucketIndex = 0; bucketIndex < this->totalNrOfBuckets;
           ++bucketIndex) {
        std::vector<size_t> vectorToPlace = std::vector<size_t>();
        // reserve a sensible capacity as estimated
        vectorToPlace.reserve((coordinates.size() / 3) /
                              this->totalNrOfBuckets);
        this->neighbourBuckets.push_back(vectorToPlace);
      }
      Eigen::Array3i indexTriplet;
      for (size_t bucketIndexX = 0; bucketIndexX < this->nrOfBucketsX;
           ++bucketIndexX) {
        indexTriplet[0] = bucketIndexX;
        for (size_t bucketIndexY = 0; bucketIndexY < this->nrOfBucketsY;
             ++bucketIndexY) {
          indexTriplet[1] = bucketIndexY;
          for (size_t bucketIndexZ = 0; bucketIndexZ < this->nrOfBucketsZ;
               ++bucketIndexZ) {
            indexTriplet[2] = bucketIndexZ;
            this->neighbourBucketNeighboursDefaultCutoff.push_back(
              this->getCombinedBucketIndicesForBucket(indexTriplet, cutoff));
          }
        }
      }
      assert(this->neighbourBucketNeighboursDefaultCutoff.size() ==
             this->totalNrOfBuckets);

      // fill neighbour buckets
      this->neighbourBucketSizes = Eigen::ArrayXi::Zero(this->totalNrOfBuckets);
      for (size_t i = 0; i < (coordinates.size() / 3); ++i) {
        size_t bucketIndex = this->getBucketIndexForTriplet(
          this->getBucketTripletForCoordinates(coordinates.segment(3 * i, 3)));
        this->neighbourBuckets[bucketIndex].push_back(i);
        this->neighbourBucketSizes[bucketIndex]++;
      }
    };

    void resetCoordinates(Eigen::VectorXd& newCoordinates)
    {
      this->neighbourBucketSizes = Eigen::ArrayXi::Zero(this->totalNrOfBuckets);
      for (size_t i = 0; i < (newCoordinates.size() / 3); ++i) {
        size_t bucketIndex =
          this->getBucketIndexForTriplet(this->getBucketTripletForCoordinates(
            newCoordinates.segment(3 * i, 3)));
        if (this->neighbourBuckets[bucketIndex].size() <=
            this->neighbourBucketSizes[bucketIndex]) {
          this->neighbourBuckets[bucketIndex].push_back(i);
        } else {
          this->neighbourBuckets[bucketIndex]
                                [this->neighbourBucketSizes[bucketIndex]] = i;
        }
        this->neighbourBucketSizes[bucketIndex]++;
      }
    }

    Eigen::ArrayXi getIndicesCloseToCoordinates(
      Eigen::Vector3d coordinates) const
    {
      return this->getIndicesCloseToCoordinates(coordinates, this->cutoff);
    }

    /**
     * @brief Get the Indices of Coordinates Close To the Coordinates of another
     * Index
     *
     * This function must be O(1), otherwise, this whole neighbor list will be
     * useless
     *
     * NOTE: the resulting list will not be reduced, i.e., it will contain
     * indices that have a distance > upperCutoff. Additionally,
     * the requested coordinates will also be included!
     *
     * @param idx
     * @param upperCutoff
     * @return Eigen::ArrayXi
     */
    Eigen::ArrayXi getIndicesCloseToCoordinates(
      const Eigen::Vector3d& coordinates,
      const double upperCutoff) const
    {
#ifndef NDEBUG
      INVALIDARG_EXP_IFN(upperCutoff > 0.0,
                         "Expected upper cutoff > 0., got " +
                           std::to_string(upperCutoff) + ".");
#endif
      std::vector<size_t> bucketIndices;
      if (upperCutoff == this->cutoff) {
        bucketIndices = this->neighbourBucketNeighboursDefaultCutoff
                          [this->getBucketIndexForTriplet(
                            this->getBucketTripletForCoordinates(coordinates))];
      } else {
        bucketIndices = this->getCombinedBucketIndicesForCoordinates(
          coordinates, upperCutoff);
      }

      // first, count the number of results we will get
      int nResults = this->neighbourBucketSizes(bucketIndices).sum();

      Eigen::ArrayXi results = Eigen::ArrayXi(nResults);

      Eigen::Vector3d difference;
      int results_idx = 0;
      for (size_t bucketIndex : bucketIndices) {
        for (size_t indexInBucket = 0;
             indexInBucket < this->neighbourBucketSizes[bucketIndex];
             indexInBucket++) {
          size_t atomIndex = this->neighbourBuckets[bucketIndex][indexInBucket];
          results[results_idx++] = atomIndex;
          // difference = this->coordinates.segment(3 * idx, 3) -
          //                this->coordinates.segment(3 * atomIndex, 3)
          // this->box.handlePBC(difference);
          // double distance = difference.squaredNorm();
          // if (distance < upperCutoff2 && distance >= lowerCutoff2) {
          //   results.push_back(this->atoms[atomIndex]);
          // }
        }
      }

      // return results
      return results;
    }

    /**
     * @brief Get the Indices of Coordinates Close To the Coordinates of another
     * Index
     *
     * This function must be O(1), otherwise, this whole neighbor list will be
     * useless
     *
     * NOTE: the resulting list will not be reduced, i.e., it will contain
     * indices that have a distance > upperCutoff. Additionally,
     * the requested coordinates will also be included!
     * 
     * The results will be written to the argument `result`, which will be resized if needed.
     * It will not be downsized, only upsized; the returned int is the actual remaining size.
     */
    int getIndicesCloseToCoordinates(Eigen::ArrayXi& result,
                                     const Eigen::Vector3d coordinates,
                                     const double upperCutoff)
    {

#ifndef NDEBUG
      INVALIDARG_EXP_IFN(upperCutoff > 0.0,
                         "Expected upper cutoff > 0., got " +
                           std::to_string(upperCutoff) + ".");
#endif

      std::vector<size_t> bucketIndices;
      if (upperCutoff == this->cutoff) {
        bucketIndices = this->neighbourBucketNeighboursDefaultCutoff
                          [this->getBucketIndexForTriplet(
                            this->getBucketTripletForCoordinates(coordinates))];
      } else {
        bucketIndices = this->getCombinedBucketIndicesForCoordinates(
          coordinates, upperCutoff);
      }

      // first, count the number of results we will get
      int nResults = this->neighbourBucketSizes(bucketIndices).sum();

      if (nResults >= result.size()) {
        result.conservativeResize(nResults);
      }

      Eigen::Vector3d difference;
      int results_idx = 0;
      for (size_t bucketIndex : bucketIndices) {
        for (size_t indexInBucket = 0;
             indexInBucket < this->neighbourBucketSizes[bucketIndex];
             indexInBucket++) {
          size_t atomIndex = this->neighbourBuckets[bucketIndex][indexInBucket];
          result[results_idx++] = atomIndex;
          // difference = this->coordinates.segment(3 * idx, 3) -
          //                this->coordinates.segment(3 * atomIndex, 3)
          // this->box.handlePBC(difference);
          // double distance = difference.squaredNorm();
          // if (distance < upperCutoff2 && distance >= lowerCutoff2) {
          //   results.push_back(this->atoms[atomIndex]);
          // }
        }
      }

      // return results
      return results_idx;
    }

  protected:
    size_t normalizeBucketIndex(long int bucketIndex, size_t nrOfBuckets) const
    {
      while (bucketIndex < 0) {
        bucketIndex = bucketIndex + nrOfBuckets;
      }
      while (bucketIndex >= nrOfBuckets) {
        bucketIndex = bucketIndex - nrOfBuckets;
      }
      return static_cast<size_t>(bucketIndex);
    }

    size_t getBucketIndexForTriplet(Eigen::Array3i ind) const
    {
      size_t bucketIndexX =
        this->normalizeBucketIndex(ind[0], this->nrOfBucketsX);
      size_t bucketIndexY =
        this->normalizeBucketIndex(ind[0], this->nrOfBucketsY);
      size_t bucketIndexZ =
        this->normalizeBucketIndex(ind[0], this->nrOfBucketsZ);
      return bucketIndexX + bucketIndexY * this->nrOfBucketsX +
             bucketIndexZ * this->nrOfBucketsX * this->nrOfBucketsY;
    }

    Eigen::Array3i getBucketTripletForCoordinates(
      const Eigen::Vector3d& coordinates) const
    {
      return (coordinates.array() / this->bucketWidths).floor().cast<int>();
    }

    std::vector<size_t> getCombinedBucketIndicesForBucket(
      const Eigen::ArrayXi indexBasis,
      double newCutoff) const
    {
      std::vector<size_t> result = std::vector<size_t>();

      int nrOfBucketsPerSide = (newCutoff <= this->cutoff)
                                 ? 3
                                 : std::ceil(3 * newCutoff / this->cutoff);
      if (nrOfBucketsPerSide % 2 == 0) {
        nrOfBucketsPerSide += 1;
      }
      int nrOfBucketsPerQuarter = (nrOfBucketsPerSide - 1) / 2;

      result.reserve(nrOfBucketsPerQuarter * nrOfBucketsPerQuarter *
                     nrOfBucketsPerQuarter);
      Eigen::Array3i offset = Eigen::Array3i::Zero();
      for (int offsetX = -nrOfBucketsPerQuarter;
           offsetX <= nrOfBucketsPerQuarter;
           offsetX++) {
        offset[0] = offsetX;
        for (int offsetY = -nrOfBucketsPerQuarter;
             offsetY <= nrOfBucketsPerQuarter;
             offsetY++) {
          offset[1] = offsetY;
          for (int offsetZ = -nrOfBucketsPerQuarter;
               offsetZ <= nrOfBucketsPerQuarter;
               offsetZ++) {
            offset[2] = offsetZ;
            size_t newIndex =
              this->getBucketIndexForTriplet(indexBasis + offset);
            result.push_back(newIndex);
            if (offsetX == 0 && offsetY == 0 && offsetZ == 0) {
              assert(newIndex == this->getBucketIndexForTriplet(indexBasis));
            }
          }
        }
      }

      return result;
    }

    std::vector<size_t> getCombinedBucketIndicesForCoordinates(
      const Eigen::Vector3d& coordinates,
      double newCutoff) const
    {
      Eigen::Array3i indexBasis =
        this->getBucketTripletForCoordinates(coordinates);
      return this->getCombinedBucketIndicesForBucket(indexBasis, newCutoff);
    }

  private:
    Eigen::Array3d bucketWidths;

    size_t nrOfBucketsX;
    size_t nrOfBucketsY;
    size_t nrOfBucketsZ;
    size_t totalNrOfBuckets;

    double cutoff;

    pylimer_tools::entities::Box box;

    std::vector<std::vector<size_t>> neighbourBuckets;
    std::vector<std::vector<size_t>> neighbourBucketNeighboursDefaultCutoff;
    Eigen::VectorXi neighbourBucketSizes;
  };
};
}

#endif
