#include "EigenNeighbourList.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include "Box.h"
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

namespace pylimer_tools {
namespace entities {
  EigenNeighbourList::EigenNeighbourList(const Eigen::VectorXd& coordinates,
                                         const Box& boxGeometry,
                                         const double cutoffDistance,
                                         const double scalingParam)
  {
    this->initialize(coordinates, boxGeometry, cutoffDistance, scalingParam);
  }

  bool EigenNeighbourList::operator==(const EigenNeighbourList& other) const
  {
    // Compare basic properties
    if (!APPROX_EQUAL(cutoff, other.cutoff, 1e-9) ||
        !APPROX_EQUAL(scalingFactor, other.scalingFactor, 1e-9) ||
        !APPROX_EQUAL(actualCutoff, other.actualCutoff, 1e-9) ||
        totalNrOfBuckets != other.totalNrOfBuckets || box != other.box) {
      return false;
    }

    // Compare Eigen arrays
    if (!bucketWidths.isApprox(other.bucketWidths) ||
        (nrOfBuckets != other.nrOfBuckets).any() ||
        (neighbourBucketSizes.array() != other.neighbourBucketSizes.array())
          .any()) {
      return false;
    }

    // Compare bucket data
    if (neighbourBuckets.size() != other.neighbourBuckets.size()) {
      return false;
    }

    for (size_t i = 0; i < neighbourBuckets.size(); ++i) {
      if (neighbourBuckets[i].size() != other.neighbourBuckets[i].size()) {
        return false;
      }

      for (size_t j = 0; j < neighbourBuckets[i].size(); ++j) {
        if (neighbourBuckets[i][j] != other.neighbourBuckets[i][j]) {
          return false;
        }
      }
    }

    // Compare bucket neighbors
    if (neighbourBucketNeighboursDefaultCutoff.size() !=
        other.neighbourBucketNeighboursDefaultCutoff.size()) {
      return false;
    }

    for (size_t i = 0; i < neighbourBucketNeighboursDefaultCutoff.size(); ++i) {
      if (neighbourBucketNeighboursDefaultCutoff[i].size() !=
          other.neighbourBucketNeighboursDefaultCutoff[i].size()) {
        return false;
      }

      for (size_t j = 0; j < neighbourBucketNeighboursDefaultCutoff[i].size();
           ++j) {
        if (neighbourBucketNeighboursDefaultCutoff[i][j] !=
            other.neighbourBucketNeighboursDefaultCutoff[i][j]) {
          return false;
        }
      }
    }

    return true;
  }

  bool EigenNeighbourList::operator!=(const EigenNeighbourList& other) const
  {
    return !(*this == other);
  }

  void EigenNeighbourList::initialize(const Eigen::VectorXd& coordinates,
                                      const Box& newBox,
                                      const double newCutoff,
                                      const double newScalingFactor)
  {
    INVALIDARG_EXP_IFN(newCutoff > 1e-3, "Cutoff must be larger than zero");
    INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                       "Coordinates must be in three");

    this->actualCutoff = newCutoff / newScalingFactor;
    this->cutoff = newCutoff;
    this->scalingFactor = newScalingFactor;
    this->box = newBox;

    this->nrOfBuckets =
      (this->box.getL().array() / this->actualCutoff).floor().cast<long int>();
    Eigen::Array3li ones = Eigen::Array3li::Ones();
    this->nrOfBuckets = this->nrOfBuckets.max(ones);

    this->bucketWidths = this->box.getL() / this->nrOfBuckets.cast<double>();

    this->totalNrOfBuckets = this->nrOfBuckets.prod();

    this->neighbourBuckets.clear();
    this->neighbourBuckets.reserve(this->totalNrOfBuckets);
    this->neighbourBucketNeighboursDefaultCutoff.clear();
    this->neighbourBucketNeighboursDefaultCutoff.reserve(
      this->totalNrOfBuckets);
    // std::cout << "Preparing " << this->totalNrOfBuckets << " bins"
    //           << std::endl;

    // prepare the default buckets
    long int nB1x2 = (this->nrOfBuckets[0] * this->nrOfBuckets[1]);
    int currentX = this->nrOfBuckets[0] - 1;
    int currentY = nB1x2 - 1;
    int currentZ = this->totalNrOfBuckets - 1;
    for (bucket_idx_t bucketIndex = 0;
         bucketIndex < static_cast<bucket_idx_t>(this->totalNrOfBuckets);
         ++bucketIndex) {
      currentX += 1;
      if (bucketIndex % this->nrOfBuckets[0] == 0) {
        currentY += 1;
      }
      if (bucketIndex % nB1x2 == 0) {
        currentZ += 1;
      }

      assert(currentX >= 0 && currentY >= 0 && currentZ >= 0);
      std::vector<bucket_idx_t> vectorToPlace = std::vector<bucket_idx_t>();
      // reserve a sensible capacity as estimated
      vectorToPlace.reserve((coordinates.size() / 3) / this->totalNrOfBuckets);
      this->neighbourBuckets.push_back(vectorToPlace);
      // assemble the neighbouring buckets with the default cut-off
      std::vector<bucket_idx_t> neighbouringBuckets;
      int lowerLim = scalingFactor * -1;
      int upperLim = scalingFactor * 1;
      neighbouringBuckets.reserve(
        (upperLim - lowerLim) * (upperLim - lowerLim) * (upperLim - lowerLim));
      for (int dx = lowerLim; dx <= upperLim; ++dx) {
        for (int dy = lowerLim; dy <= upperLim; ++dy) {
          for (int dz = lowerLim; dz <= upperLim; ++dz) {
            neighbouringBuckets.push_back(
              (currentX + dx) % this->nrOfBuckets[0] +
              ((currentY + dy) * this->nrOfBuckets[0]) % nB1x2 +
              ((currentZ + dz) * nB1x2) % this->totalNrOfBuckets);
          }
        }
      }
      this->neighbourBucketNeighboursDefaultCutoff.push_back(
        neighbouringBuckets);

      // /
      // Eigen::Vector3d centralCoordinates =
      //   this->getCentralCoordinatesOfBucket(bucketIndex);
      // this->neighbourBucketNeighboursDefaultCutoff.push_back(
      //   this->getCombinedBucketIndicesForCoordinates(
      //     centralCoordinates, this->cutoff, true));
    }

    assert(this->neighbourBucketNeighboursDefaultCutoff.size() ==
           this->totalNrOfBuckets);
    assert(this->neighbourBuckets.size() == this->totalNrOfBuckets);
    // std::cout << "Filling buckets" << std::endl;

    // fill neighbour buckets
    this->neighbourBucketSizes = Eigen::ArrayXi::Zero(this->totalNrOfBuckets);
    for (Eigen::Index i = 0; i < (coordinates.size() / 3); ++i) {
      bucket_idx_t bucketIndex = this->getBucketIndexForTriplet(
        this->getBucketTripletForCoordinates(coordinates.segment(3 * i, 3)));
      assert(bucketIndex >= 0 && bucketIndex < this->totalNrOfBuckets);
      this->neighbourBuckets[bucketIndex].push_back(i);
      this->neighbourBucketSizes[bucketIndex]++;
    }
  };

  /**
   * @brief Re-bin with a new set of coordinates
   *
   * @param newCoordinates
   */
  void EigenNeighbourList::resetCoordinates(Eigen::VectorXd& newCoordinates)
  {
    // just override all the buckets.
    this->neighbourBucketSizes.setZero();
    for (Eigen::Index i = 0; i < (newCoordinates.size() / 3); ++i) {
      int bucketIndex = this->getBucketIndexForTriplet(
        this->getBucketTripletForCoordinates(newCoordinates.segment(3 * i, 3)));
      if (static_cast<int>(this->neighbourBuckets[bucketIndex].size()) <=
          this->neighbourBucketSizes[bucketIndex]) {
        this->neighbourBuckets[bucketIndex].push_back(i);
      } else {
        this->neighbourBuckets[bucketIndex]
                              [this->neighbourBucketSizes[bucketIndex]] = i;
      }
      this->neighbourBucketSizes[bucketIndex]++;
    }
  }

  bool EigenNeighbourList::checkIfCoordinatesAreCurrent(
    Eigen::VectorXd& newCoordinates)
  {
    Eigen::ArrayXi neighbourBucketSizes2 =
      Eigen::ArrayXi::Zero(this->totalNrOfBuckets);
    std::unordered_map<bucket_idx_t, std::vector<coordinate_idx_t>>
      neighbourBuckets2;
    for (size_t i = 0; i < (newCoordinates.size() / 3); ++i) {
      int bucketIndex = this->getBucketIndexForTriplet(
        this->getBucketTripletForCoordinates(newCoordinates.segment(3 * i, 3)));
      if (neighbourBuckets2[bucketIndex].size() <=
          neighbourBucketSizes2[bucketIndex]) {
        neighbourBuckets2[bucketIndex].push_back(i);
      } else {
        neighbourBuckets2[bucketIndex][neighbourBucketSizes2[bucketIndex]] = i;
      }
      neighbourBucketSizes2[bucketIndex]++;
    }
    if (!(neighbourBucketSizes2 == this->neighbourBucketSizes.array()).all()) {
      return false;
    }
    for (auto& it : neighbourBuckets2) {
      for (size_t i = 0; i < neighbourBuckets2[it.first].size(); ++i) {
        if (neighbourBuckets2[it.first][i] !=
            this->neighbourBuckets[it.first][i]) {
          return false;
        }
      }
    }
    return true;
  }

  void EigenNeighbourList::validateWhyNotIncluded(Eigen::Vector3d sourceCoords,
                                                  Eigen::Vector3d targetCoords,
                                                  double newCutoff) const
  {
    if (newCutoff <= 0.) {
      newCutoff = this->cutoff;
    }
    Eigen::Array3li baseTriplet =
      this->getBucketTripletForCoordinates(sourceCoords);
    int baseBucketIndex = this->getBucketIndexForTriplet(baseTriplet);

    Eigen::Array3li targetTriplet =
      this->getBucketTripletForCoordinates(targetCoords);
    int targetBucketIndex = this->getBucketIndexForTriplet(targetTriplet);

    bool foundSource = false;
    bool foundTarget = false;
    std::vector<bucket_idx_t> bucketIdxs;
    if (newCutoff == this->cutoff) {
      bucketIdxs =
        this->neighbourBucketNeighboursDefaultCutoff[baseBucketIndex];
      for (bucket_idx_t bucketIndex : bucketIdxs) {
        foundSource = foundSource || (baseBucketIndex == bucketIndex);
        foundTarget = foundTarget || (targetBucketIndex == bucketIndex);
      }
      if (!foundTarget || !foundSource) {
        std::cerr << "Used default neighbour cutoff's list" << std::endl;
      }
    } else {
      Eigen::Array3li maxIndices = this->getBucketTripletForCoordinates(
        sourceCoords + Eigen::Vector3d::Constant(newCutoff));
      Eigen::Array3li minIndices = this->getBucketTripletForCoordinates(
        sourceCoords - Eigen::Vector3d::Constant(newCutoff));
      Eigen::Array3li indexTriplet;
      // now, do permutations of all these
      for (long int i = minIndices[0]; i <= maxIndices[0]; ++i) {
        for (long int j = minIndices[1]; j <= maxIndices[1]; ++j) {
          for (long int k = minIndices[2]; k <= maxIndices[2]; ++k) {
            indexTriplet[0] = i;
            indexTriplet[1] = j;
            indexTriplet[2] = k;
            bucket_idx_t bucketIndex =
              this->getBucketIndexForTriplet(indexTriplet);
            bucketIdxs.push_back(bucketIndex);
            foundSource = foundSource || (baseBucketIndex == bucketIndex);
            foundTarget = foundTarget || (targetBucketIndex == bucketIndex);
          }
        }
      }

      if (!foundTarget || !foundSource) {
        std::cerr << "With source indices " << baseTriplet << ", got max "
                  << maxIndices << " and min " << minIndices
                  << " but not target " << targetTriplet << std::endl;
        this->normalizeTriplet(baseTriplet);
        this->normalizeTriplet(maxIndices);
        this->normalizeTriplet(minIndices);
        this->normalizeTriplet(targetTriplet);
        std::cerr << "Normalized, that's " << baseTriplet << ", got max "
                  << maxIndices << " and min " << minIndices << ", target "
                  << targetTriplet << std::endl;
      }
    }
    if (!foundSource) {
      std::cerr << "Couldn't find source bucket " << baseBucketIndex
                << " in list: "
                << pylimer_tools::utils::join(
                     bucketIdxs.begin(), bucketIdxs.end(), std::string(", "))
                << std::endl;
    }
    if (!foundTarget) {
      std::cerr << "Target bucket " << targetBucketIndex
                << " was not found in list: "
                << pylimer_tools::utils::join(
                     bucketIdxs.begin(), bucketIdxs.end(), std::string(", "))
                << std::endl;
    }
    if (foundSource && foundTarget) {
      std::cout << "Target coords should actually be included." << std::endl;
    }
  }

  /**
   * @brief Get the indices close to given coordinates with the default cut-off
   *
   * @param coordinates
   * @return Eigen::ArrayXi
   */
  Eigen::ArrayXi EigenNeighbourList::getIndicesCloseToCoordinates(
    Eigen::Vector3d coordinates,
    double newCutoff) const
  {
    if (newCutoff == -1.) {
      newCutoff = this->cutoff;
    }
    Eigen::ArrayXi result = Eigen::ArrayXi(12);
    int num =
      this->getIndicesCloseToCoordinates(result, coordinates, newCutoff);
    result.conservativeResize(num);
    return result;
  }

  /**
   * @brief Get the number of coordinates actually stored/binned
   *
   * @return long int
   */
  long int EigenNeighbourList::getNumBinnedCoordinates() const
  {
    return this->neighbourBucketSizes.sum();
  }

  /**
   * @brief Similar to getIndicesCloseToCoordinates, this function returns those
   * coordinates, but filtered for cutoff and higher index
   *
   * @param result the array to write the coordinates into
   * @param coordinates the coordinates underlying this neighbour list
   * @param source the source index. Only resulting atoms with higher index are
   * returned
   * @param cutoff the cut-off to use. Only atoms within this distance are
   * returned in results.
   * @return int the number of first (valid) entries in result.
   */
  int EigenNeighbourList::getHigherIndicesWithinCutoff(
    Eigen::ArrayXi& result,
    const Eigen::VectorXd& coordinates,
    const int source,
    const double newCutoff) const
  {
    int numNeighbours = this->getIndicesCloseToCoordinates(
      result, coordinates.segment(3 * source, 3), newCutoff);
    int actualNumNeighbours = 0;
    const double cutoff2 = newCutoff * newCutoff;
    for (size_t i = 0; i < numNeighbours; ++i) {
      if (result[i] > source) {
        Eigen::Vector3d dist = coordinates.segment(3 * source, 3) -
                               coordinates.segment(result[i] * 3, 3);
        this->box.handlePBC(dist);
        if (dist.squaredNorm() < cutoff2) {
          result[actualNumNeighbours] = result[i];
          actualNumNeighbours += 1;
        }
      }
    }
    return actualNumNeighbours;
  };

  /**
   * @brief Get the Indices of Coordinates Close To the Coordinates of another
   * Index
   *
   * This function must be O(1), otherwise, this whole neighbor list will be
   * useless
   *
   * NOTE: The resulting list will not be reduced, i.e., it will contain
   * indices that have a distance > upperCutoff. Additionally,
   * the requested coordinates will also be included!
   *
   * The results will be written to the argument `result`, which will be
   * resized if needed. It will not be downsized, only upsized; the returned
   * int is the actual remaining size.
   */
  int EigenNeighbourList::getIndicesCloseToCoordinates(
    Eigen::ArrayXi& result,
    const Eigen::Vector3d coordinates,
    const double upperCutoff,
    const bool expectDefault) const
  {
    if (this->totalNrOfBuckets < 2) {
      size_t nResults = this->neighbourBucketSizes.sum();
      result = Eigen::ArrayXi::LinSpaced(nResults, 0, nResults - 1);
      return result.size();
    }

#ifndef NDEBUG
    INVALIDARG_EXP_IFN(upperCutoff > 0.0,
                       "Expected upper cutoff > 0., got " +
                         std::to_string(upperCutoff) + ".");
#endif

    int results_idx = 0;
    if (APPROX_EQUAL(upperCutoff, this->cutoff, 1e-12)) {
      const bucket_idx_t coordinatesBucketIdx = this->getBucketIndexForTriplet(
        this->getBucketTripletForCoordinates(coordinates));
      // first, count the number of results we will get
      long int nResults =
        this
          ->neighbourBucketSizes(
            this->neighbourBucketNeighboursDefaultCutoff[coordinatesBucketIdx])
          .sum();

      if (result.size() < nResults) {
        // heuristic minimum
        result.conservativeResize(nResults);
      }
      for (bucket_idx_t bucketIndex :
           this->neighbourBucketNeighboursDefaultCutoff[coordinatesBucketIdx]) {
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
      INVALIDARG_EXP_IFN(!expectDefault,
                         "Expected default cutoff, but did not get it");
      // TODO: this is more or less identical to
      // EigenNeighbourList::getCombinedIndicesForCoordinates
      // with some minor additions here and there
      Eigen::Array3li maxIndices = this->getBucketTripletForCoordinates(
        coordinates + Eigen::Vector3d::Constant(upperCutoff));
      Eigen::Array3li minIndices = this->getBucketTripletForCoordinates(
        coordinates - Eigen::Vector3d::Constant(upperCutoff));

      if (result.size() < 12) {
        // heuristic minimum
        result.conservativeResize(12);
      }

      // use a bitset to avoid returning duplicates
      std::vector<bool> visitedBuckets =
        pylimer_tools::utils::initializeWithValue<bool>(this->totalNrOfBuckets,
                                                        false);
      Eigen::Array3li indexTriplet;
      // now, do permutations of all these
      assert(minIndices[0] <= maxIndices[0]);
      assert(minIndices[1] <= maxIndices[1]);
      assert(minIndices[2] <= maxIndices[2]);
      for (long int i = minIndices[0]; i <= maxIndices[0]; ++i) {
        for (long int j = minIndices[1]; j <= maxIndices[1]; ++j) {
          for (long int k = minIndices[2]; k <= maxIndices[2]; ++k) {
            indexTriplet[0] = i;
            indexTriplet[1] = j;
            indexTriplet[2] = k;
            //
            bucket_idx_t bucketIndex =
              this->getBucketIndexForTriplet(indexTriplet);
            // found the bucket, insert its contents into the results if
            // desired
            if (bucketIndex >= 0 && !visitedBuckets[bucketIndex]) {
              visitedBuckets[bucketIndex] = true;
              for (int indexInBucket = 0;
                   indexInBucket < this->neighbourBucketSizes[bucketIndex];
                   ++indexInBucket) {
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
   * @brief For debugging/test purposes: The actual buckets
   *
   * @return std::vector<std::vector<coordinate_idx_t>>
   */
  std::vector<std::vector<coordinate_idx_t>>
  EigenNeighbourList::getNeighbourBuckets() const
  {
    return this->neighbourBuckets;
  }
  Eigen::VectorXi EigenNeighbourList::getNeighbourBucketSizes() const
  {
    return this->neighbourBucketSizes;
  }
  Eigen::Vector3d EigenNeighbourList::getCentralCoordinatesOfBucket(
    const int bucketIndex) const
  {
    Eigen::Array3li coeffs = this->tripletFromIndex(bucketIndex);
    Eigen::Vector3d results = coeffs.cast<double>() * this->bucketWidths +
                              0.5 * this->bucketWidths + this->box.getLowL();

#ifndef NDEBUG
    Eigen::Array3li verifyCoeffs =
      this->getBucketTripletForCoordinates(results);
    this->normalizeTriplet(verifyCoeffs);
    int verifyBucketIdx = this->getBucketIndexForTriplet(coeffs);
    assert(verifyBucketIdx == bucketIndex);
    assert(verifyCoeffs.isApprox(coeffs));
#endif
    return results;
  }

  std::vector<bucket_idx_t>
  EigenNeighbourList::getCombinedBucketIndicesForCoordinates(
    const Eigen::Vector3d& coordinates,
    const double newCutoff,
    const bool sort) const
  {
    Eigen::Array3li indexBasis =
      this->getBucketTripletForCoordinates(coordinates);
    Eigen::Array3li maxIndices = this->getBucketTripletForCoordinates(
      coordinates + Eigen::Vector3d::Constant(newCutoff));
    Eigen::Array3li minIndices = this->getBucketTripletForCoordinates(
      coordinates - Eigen::Vector3d::Constant(newCutoff));

    // use a set to avoid returning duplicates
    std::set<bucket_idx_t> result;
    Eigen::Array3li indexTriplet;
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
    if (sort) {
      std::sort(res.begin(), res.end());
    }
    return res;
  }
  //////////////////////////////////////////////////////////
  // protected:
  //////////////////////////////////////////////////////////
  /**
   * @brief Do "PBC" with a bucket index in one direction
   *
   * @param bucketIndex
   * @param nrOfBuckets
   * @return size_t
   */
  bucket_idx_t EigenNeighbourList::normalizeBucketIndex(
    long int bucketIndex,
    const size_t newNrOfBuckets) const
  {
    bucketIndex %= static_cast<long int>(newNrOfBuckets);
    bucketIndex += newNrOfBuckets * static_cast<long int>(bucketIndex < 0);
    return static_cast<bucket_idx_t>(bucketIndex);
  }

  void EigenNeighbourList::normalizeTriplet(Eigen::Array3li& triplet) const
  {
    triplet = (triplet - (triplet / this->nrOfBuckets) * this->nrOfBuckets);
    triplet += this->nrOfBuckets * (triplet < 0).cast<long int>();
  }

  /**
   * @brief Decompose the one index into a triplet of indices
   *
   * @param index
   * @return Eigen::Array3li
   */
  Eigen::Array3li EigenNeighbourList::tripletFromIndex(
    const bucket_idx_t index) const
  {
    bucket_idx_t bucketIndexZ =
      std::floor(index / (this->nrOfBuckets[0] * this->nrOfBuckets[1]));
    bucket_idx_t bucketIndexY = std::floor(
      (index - bucketIndexZ * (this->nrOfBuckets[0] * this->nrOfBuckets[1])) /
      (this->nrOfBuckets[0]));
    bucket_idx_t bucketIndexX =
      index - bucketIndexZ * (this->nrOfBuckets[0] * this->nrOfBuckets[1]) -
      bucketIndexY * this->nrOfBuckets[0];
    return Eigen::Array3li(bucketIndexX, bucketIndexY, bucketIndexZ);
  }

  /**
   * @brief Get one index, normalized, given a triplet of indices
   *
   * @param ind
   * @return bucket_idx_t
   */
  bucket_idx_t EigenNeighbourList::getBucketIndexForTriplet(
    Eigen::Array3li ind) const
  {
    this->normalizeTriplet(ind);
    return ind[0] + ind[1] * this->nrOfBuckets[0] +
           ind[2] * this->nrOfBuckets[0] * this->nrOfBuckets[1];
    // bucket_idx_t bucketIndexX =
    //   this->normalizeBucketIndex(ind[0], this->nrOfBuckets[0]);
    // bucket_idx_t bucketIndexY =
    //   this->normalizeBucketIndex(ind[1], this->nrOfBuckets[1]);
    // bucket_idx_t bucketIndexZ =
    //   this->normalizeBucketIndex(ind[2], this->nrOfBuckets[2]);
    // return bucketIndexX + bucketIndexY * this->nrOfBuckets[0] +
    //        bucketIndexZ * this->nrOfBuckets[0] * this->nrOfBuckets[1];
  }

  Eigen::Array3li EigenNeighbourList::getBucketTripletForCoordinates(
    const Eigen::Vector3d& coordinates) const
  {
    return ((coordinates.array() + this->box.getLowL()) / this->bucketWidths)
      .floor()
      .cast<long int>();
  }
};
};
