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
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include <algorithm>

namespace pylimer_tools {
namespace entities {

  typedef int bucket_idx_t;
  typedef int coordinate_idx_t;

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

      this->nrOfBuckets = (box.getL().array() / cutoff).floor().cast<int>();

      this->bucketWidths = box.getL() / this->nrOfBuckets.cast<double>();

      this->totalNrOfBuckets = this->nrOfBuckets.prod();

      this->neighbourBuckets.reserve(this->totalNrOfBuckets);
      this->neighbourBucketNeighboursDefaultCutoff.reserve(
        this->totalNrOfBuckets);
      // std::cout << "Preparing " << this->totalNrOfBuckets << " bins"
      //           << std::endl;

      // prepare the buckets
      for (bucket_idx_t bucketIndex = 0; bucketIndex < this->totalNrOfBuckets;
           ++bucketIndex) {
        std::vector<bucket_idx_t> vectorToPlace = std::vector<bucket_idx_t>();
        // reserve a sensible capacity as estimated
        vectorToPlace.reserve((coordinates.size() / 3) /
                              this->totalNrOfBuckets);
        this->neighbourBuckets.push_back(vectorToPlace);
        Eigen::Vector3d centralCoordinates =
          this->getCentralCoordinatesOfBucket(bucketIndex);
        this->neighbourBucketNeighboursDefaultCutoff.push_back(
          this->getCombinedBucketIndicesForCoordinates(centralCoordinates,
                                                       cutoff));
      }

      assert(this->neighbourBucketNeighboursDefaultCutoff.size() ==
             this->totalNrOfBuckets);
      // std::cout << "Filling buckets" << std::endl;

      // fill neighbour buckets
      this->neighbourBucketSizes = Eigen::ArrayXi::Zero(this->totalNrOfBuckets);
      for (size_t i = 0; i < (coordinates.size() / 3); ++i) {
        int bucketIndex = this->getBucketIndexForTriplet(
          this->getBucketTripletForCoordinates(coordinates.segment(3 * i, 3)));
        this->neighbourBuckets[bucketIndex].push_back(i);
        this->neighbourBucketSizes[bucketIndex]++;
      }
    };

    /**
     * @brief Re-bin with a new set of coordinates
     *
     * @param newCoordinates
     */
    void resetCoordinates(Eigen::VectorXd& newCoordinates)
    {
      // just override all the buckets.
      this->neighbourBucketSizes = Eigen::ArrayXi::Zero(this->totalNrOfBuckets);
      for (size_t i = 0; i < (newCoordinates.size() / 3); ++i) {
        int bucketIndex =
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

    bool checkIfCoordinatesAreCurrent(Eigen::VectorXd& newCoordinates)
    {
      Eigen::ArrayXi neighbourBucketSizes2 =
        Eigen::ArrayXi::Zero(this->totalNrOfBuckets);
      std::unordered_map<size_t, std::vector<size_t>> neighbourBuckets2;
      for (size_t i = 0; i < (newCoordinates.size() / 3); ++i) {
        int bucketIndex =
          this->getBucketIndexForTriplet(this->getBucketTripletForCoordinates(
            newCoordinates.segment(3 * i, 3)));
        if (neighbourBuckets2[bucketIndex].size() <=
            neighbourBucketSizes2[bucketIndex]) {
          neighbourBuckets2[bucketIndex].push_back(i);
        } else {
          neighbourBuckets2[bucketIndex][neighbourBucketSizes2[bucketIndex]] =
            i;
        }
        neighbourBucketSizes2[bucketIndex]++;
      }
      if (!(neighbourBucketSizes2 == this->neighbourBucketSizes).all()) {
        return false;
      }
      for (auto& it : neighbourBuckets2) {
        for (size_t i = 0; i < neighbourBuckets2[it.first].size(); ++i) {
          if (neighbourBuckets2[it.first][i] != this->neighbourBuckets[it.first][i]) {
            return false;
          }
        }
      }
      return true;
    }

    /**
     * @brief Get the Indices Close To Coordinates with the Default Cut-Off
     *
     * @param coordinates
     * @return Eigen::ArrayXi
     */
    Eigen::ArrayXi getIndicesCloseToCoordinates(
      Eigen::Vector3d coordinates) const
    {
      return this->getIndicesCloseToCoordinates(coordinates, this->cutoff);
    }

    long int getNumBinnedCoordinates() const
    {
      return this->neighbourBucketSizes.sum();
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
      std::vector<bucket_idx_t> bucketIndices;
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

      // Eigen::Vector3d difference;
      int results_idx = 0;
      for (bucket_idx_t bucketIndex : bucketIndices) {
        const Eigen::ArrayXi bucketContent =
          Eigen::Map<const Eigen::ArrayXi, Eigen::Unaligned>(
            this->neighbourBuckets[bucketIndex].data(),
            this->neighbourBucketSizes[bucketIndex]);
        results.segment(results_idx, this->neighbourBucketSizes[bucketIndex]) =
          bucketContent;
        results_idx += this->neighbourBucketSizes[bucketIndex];
      }

#ifndef NDEBUG
      assert(results_idx == nResults);
#endif
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
     * The results will be written to the argument `result`, which will be
     * resized if needed. It will not be downsized, only upsized; the returned
     * int is the actual remaining size.
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

      int results_idx = 0;
      if (upperCutoff == this->cutoff) {
        // first, count the number of results we will get
        int nResults =
          this
            ->neighbourBucketSizes(
              this->neighbourBucketNeighboursDefaultCutoff
                [this->getBucketIndexForTriplet(
                  this->getBucketTripletForCoordinates(coordinates))])
            .sum();

        if (result.size() < nResults) {
          // heuristic minimum
          result.conservativeResize(nResults);
        }
        for (bucket_idx_t bucketIndex :
             this->neighbourBucketNeighboursDefaultCutoff
               [this->getBucketIndexForTriplet(
                 this->getBucketTripletForCoordinates(coordinates))]) {
          for (int indexInBucket = 0;
               indexInBucket < this->neighbourBucketSizes[bucketIndex];
               indexInBucket++) {
            coordinate_idx_t atomIndex =
              this->neighbourBuckets[bucketIndex][indexInBucket];
            result[results_idx] = atomIndex;
            results_idx += 1;
          }
        }
      } else {
        // TODO: this is more or less identical to
        // EigenNeighbourList::getCombinedIndicesForCoordinates
        // with some minor additions here and there
        Eigen::Array3i indexBasis =
          this->getBucketTripletForCoordinates(coordinates);
        Eigen::Array3i maxIndices = this->getBucketTripletForCoordinates(
          coordinates + Eigen::Vector3d::Constant(upperCutoff));
        Eigen::Array3i minIndices = this->getBucketTripletForCoordinates(
          coordinates - Eigen::Vector3d::Constant(upperCutoff));

        if (result.size() < 12) {
          // heuristic minimum
          result.conservativeResize(12);
        }

        // use a set to avoid returning duplicates
        std::set<bucket_idx_t> buckets;
        Eigen::Array3i indexTriplet;
        // now, do permutations of all these
        for (int i = minIndices[0]; i <= maxIndices[0]; ++i) {
          indexTriplet[0] = i;
          for (int j = minIndices[1]; j <= maxIndices[1]; ++j) {
            indexTriplet[1] = j;
            for (int k = minIndices[2]; k <= maxIndices[2]; ++k) {
              indexTriplet[2] = k;
              bucket_idx_t bucketIndex =
                this->getBucketIndexForTriplet(indexTriplet);
              // found the bucket, insert its contents into the results if
              // desired
              if (buckets.insert(bucketIndex).second) {
                for (int indexInBucket = 0;
                     indexInBucket < this->neighbourBucketSizes[bucketIndex];
                     indexInBucket++) {
                  coordinate_idx_t atomIndex =
                    this->neighbourBuckets[bucketIndex][indexInBucket];
                  result[results_idx] = atomIndex;
                  results_idx += 1;

                  if (results_idx >= result.size()) {
                    result.conservativeResize(result.size() * 2);
                  }
                }
              }
            }
          }
        }
      }

      // return results
      return results_idx;
    }

    /**
     * @brief For debugging/test purposes: the actual buckets
     *
     * @return std::vector<std::vector<size_t>>
     */
    std::vector<std::vector<coordinate_idx_t>> getNeighbourBuckets()
    {
      return this->neighbourBuckets;
    }
    Eigen::VectorXi getNeighbourBucketSizes()
    {
      return this->neighbourBucketSizes;
    }
    Eigen::Vector3d getCentralCoordinatesOfBucket(int bucketIndex)
    {
      Eigen::Array3i coeffs = this->tripletFromIndex(bucketIndex);
      Eigen::Vector3d results = coeffs.cast<double>() * this->bucketWidths +
                                0.5 * this->bucketWidths + this->box.getLowL();

#ifndef NDEBUG
      Eigen::Array3i verifyCoeffs =
        this->getBucketTripletForCoordinates(results);
      this->normalizeTriplet(verifyCoeffs);
      int verifyBucketIdx = this->getBucketIndexForTriplet(coeffs);
      assert(verifyBucketIdx == bucketIndex);
      assert(verifyCoeffs.isApprox(coeffs));
#endif
      return results;
    }

    std::vector<bucket_idx_t> getCombinedBucketIndicesForCoordinates(
      const Eigen::Vector3d& coordinates,
      double newCutoff) const
    {
      Eigen::Array3i indexBasis =
        this->getBucketTripletForCoordinates(coordinates);
      Eigen::Array3i maxIndices = this->getBucketTripletForCoordinates(
        coordinates + Eigen::Vector3d::Constant(newCutoff));
      Eigen::Array3i minIndices = this->getBucketTripletForCoordinates(
        coordinates - Eigen::Vector3d::Constant(newCutoff));

      // use a set to avoid returning duplicates
      std::set<bucket_idx_t> result;
      Eigen::Array3i indexTriplet;
      // now, do permutations of all these
      for (int i = minIndices[0]; i <= maxIndices[0]; ++i) {
        indexTriplet[0] = i;
        for (int j = minIndices[1]; j <= maxIndices[1]; ++j) {
          indexTriplet[1] = j;
          for (int k = minIndices[2]; k <= maxIndices[2]; ++k) {
            indexTriplet[2] = k;
            result.insert(this->getBucketIndexForTriplet(indexTriplet));
          }
        }
      }

      std::vector<bucket_idx_t> res(result.begin(), result.end());
      return res;
    }

  protected:
    /**
     * @brief Do "PBC" with a bucket index in one direction
     *
     * @param bucketIndex
     * @param nrOfBuckets
     * @return size_t
     */
    bucket_idx_t normalizeBucketIndex(long int bucketIndex,
                                      size_t nrOfBuckets) const
    {
      bucketIndex %= static_cast<long int>(nrOfBuckets);
      bucketIndex += nrOfBuckets * static_cast<long int>(bucketIndex < 0);
      return static_cast<bucket_idx_t>(bucketIndex);
    }

    Eigen::Array3i normalizeTriplet(Eigen::Array3i& triplet)
    {
      triplet = (triplet - (triplet / this->nrOfBuckets) * this->nrOfBuckets);
      triplet += this->nrOfBuckets * (triplet < 0).cast<int>();
      return triplet;
    }

    Eigen::Array3i tripletFromIndex(bucket_idx_t index) const
    {
      bucket_idx_t bucketIndexZ =
        std::floor(index / (this->nrOfBuckets[0] * this->nrOfBuckets[1]));
      bucket_idx_t bucketIndexY = std::floor(
        (index - bucketIndexZ * (this->nrOfBuckets[0] * this->nrOfBuckets[1])) /
        (this->nrOfBuckets[0]));
      bucket_idx_t bucketIndexX =
        index - bucketIndexZ * (this->nrOfBuckets[0] * this->nrOfBuckets[1]) -
        bucketIndexY * this->nrOfBuckets[0];
      return Eigen::Array3i(bucketIndexX, bucketIndexY, bucketIndexZ);
    }

    bucket_idx_t getBucketIndexForTriplet(Eigen::Array3i ind) const
    {
      bucket_idx_t bucketIndexX =
        this->normalizeBucketIndex(ind[0], this->nrOfBuckets[0]);
      bucket_idx_t bucketIndexY =
        this->normalizeBucketIndex(ind[1], this->nrOfBuckets[1]);
      bucket_idx_t bucketIndexZ =
        this->normalizeBucketIndex(ind[2], this->nrOfBuckets[2]);
      return bucketIndexX + bucketIndexY * this->nrOfBuckets[0] +
             bucketIndexZ * this->nrOfBuckets[0] * this->nrOfBuckets[1];
    }

    Eigen::Array3i getBucketTripletForCoordinates(
      const Eigen::Vector3d& coordinates) const
    {
      return ((coordinates.array() + this->box.getLowL()) / this->bucketWidths)
        .floor()
        .cast<int>();
    }

  private:
    Eigen::Array3d bucketWidths;
    Eigen::Array3i nrOfBuckets;

    size_t totalNrOfBuckets;

    double cutoff;

    pylimer_tools::entities::Box box;

    std::vector<std::vector<coordinate_idx_t>> neighbourBuckets;
    std::vector<std::vector<bucket_idx_t>>
      neighbourBucketNeighboursDefaultCutoff;
    Eigen::VectorXi neighbourBucketSizes;
  };
};
};

#endif
