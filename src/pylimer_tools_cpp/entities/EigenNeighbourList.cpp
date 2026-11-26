#include "EigenNeighbourList.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include "../utils/utilityMacros.h"
#include "Box.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

namespace pylimer_tools {
namespace entities {
  EigenNeighbourList::EigenNeighbourList(
    const Eigen::Ref<const Eigen::VectorXd>& coordinates,
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

    // Compare bucket data and bucket neighbors
    if (defaultCutoffStride != other.defaultCutoffStride ||
        flatBuckets != other.flatBuckets ||
        bucketOffsets != other.bucketOffsets ||
        neighbourBucketNeighboursDefaultCutoff !=
          other.neighbourBucketNeighboursDefaultCutoff) {
      return false;
    }

    return true;
  }

  bool EigenNeighbourList::operator!=(const EigenNeighbourList& other) const
  {
    return !(*this == other);
  }

  void EigenNeighbourList::initialize(
    const Eigen::Ref<const Eigen::VectorXd>& coordinates,
    const Box& newBox,
    const double newCutoff,
    const double newScalingFactor)
  {
    INVALIDARG_EXP_IFN(newCutoff > 1e-3, "Cutoff must be larger than zero");
    INVALIDARG_EXP_IFN(newScalingFactor > 0.0,
                       "Scaling factor must be larger than zero");
    INVALIDARG_EXP_IFN(coordinates.size() % 3 == 0,
                       "Coordinates must be in three");

    this->actualCutoff = newCutoff / newScalingFactor;
    this->cutoff = newCutoff;
    this->scalingFactor = newScalingFactor;
    this->box = newBox;

    this->nrOfBuckets =
      (this->box.getL().array() / this->actualCutoff).floor().cast<long int>();
    const Eigen::Array3li ones = Eigen::Array3li::Ones();
    this->nrOfBuckets = this->nrOfBuckets.max(ones);

    this->bucketWidths = this->box.getL() / this->nrOfBuckets.cast<double>();

    this->totalNrOfBuckets = static_cast<size_t>(this->nrOfBuckets.prod());

    // Prepare the deduplicated periodic bucket stencil for the default cutoff.
    // Along any axis where (2 * halfSpan + 1) >= nrOfBuckets[d], clamping the
    // axis span to nrOfBuckets[d] avoids duplicate bucket entries and negative
    // modulo indices while preserving a uniform stride per bucket.
    const long int halfSpan =
      static_cast<long int>(std::ceil(this->scalingFactor));
    const long int fullSpan = 2 * halfSpan + 1;
    const long int nB0 = this->nrOfBuckets[0];
    const long int nB1 = this->nrOfBuckets[1];
    const long int nB2 = this->nrOfBuckets[2];
    const long int spanX = std::min<long int>(nB0, fullSpan);
    const long int spanY = std::min<long int>(nB1, fullSpan);
    const long int spanZ = std::min<long int>(nB2, fullSpan);
    this->defaultCutoffStride = static_cast<size_t>(spanX * spanY * spanZ);
    this->neighbourBucketNeighboursDefaultCutoff.resize(
      this->totalNrOfBuckets * this->defaultCutoffStride);

    const long int nB0x1 = nB0 * nB1;
    const long int startDx = (fullSpan >= nB0) ? 0 : -halfSpan;
    const long int startDy = (fullSpan >= nB1) ? 0 : -halfSpan;
    const long int startDz = (fullSpan >= nB2) ? 0 : -halfSpan;

    for (long int cz = 0; cz < nB2; ++cz) {
      for (long int cy = 0; cy < nB1; ++cy) {
        for (long int cx = 0; cx < nB0; ++cx) {
          const size_t bucketIndex =
            static_cast<size_t>(cx + cy * nB0 + cz * nB0x1);
          size_t writeIdx = bucketIndex * this->defaultCutoffStride;
          for (long int kx = 0; kx < spanX; ++kx) {
            const long int nx =
              (fullSpan >= nB0)
                ? kx
                : this->normalizeBucketIndex(cx + startDx + kx,
                                             static_cast<size_t>(nB0));
            for (long int ky = 0; ky < spanY; ++ky) {
              const long int ny =
                (fullSpan >= nB1)
                  ? ky
                  : this->normalizeBucketIndex(cy + startDy + ky,
                                               static_cast<size_t>(nB1));
              const long int xyTerm = nx + ny * nB0;
              for (long int kz = 0; kz < spanZ; ++kz) {
                const long int nz =
                  (fullSpan >= nB2)
                    ? kz
                    : this->normalizeBucketIndex(cz + startDz + kz,
                                                 static_cast<size_t>(nB2));
                this->neighbourBucketNeighboursDefaultCutoff[writeIdx++] =
                  xyTerm + nz * nB0x1;
              }
            }
          }
        }
      }
    }

    assert(this->neighbourBucketNeighboursDefaultCutoff.size() ==
           this->totalNrOfBuckets * this->defaultCutoffStride);

    // Fill neighbour buckets using flat CSR storage
    this->neighbourBucketSizes =
      Eigen::VectorXi::Zero(static_cast<Eigen::Index>(this->totalNrOfBuckets));
    this->resetCoordinates(coordinates);
  };

  /**
   * @brief Re-bin with a new set of coordinates
   *
   * @param newCoordinates
   */
  void EigenNeighbourList::resetCoordinates(
    const Eigen::Ref<const Eigen::VectorXd>& newCoordinates)
  {
    INVALIDARG_EXP_IFN(newCoordinates.size() % 3 == 0,
                       "Coordinates must be in three");
    const Eigen::Index numAtoms = newCoordinates.size() / 3;
    this->atomBucketsScratch.resize(static_cast<size_t>(numAtoms));

    this->neighbourBucketSizes.setZero();
    for (Eigen::Index i = 0; i < numAtoms; ++i) {
      const bucket_idx_t bucketIndex = this->getBucketIndexForTriplet(
        this->getBucketTripletForCoordinates(newCoordinates.segment<3>(3 * i)));
      assert(bucketIndex >= 0 &&
             static_cast<size_t>(bucketIndex) < this->totalNrOfBuckets);
      this->atomBucketsScratch[static_cast<size_t>(i)] = bucketIndex;
      this->neighbourBucketSizes[bucketIndex]++;
    }

    this->bucketOffsets.resize(this->totalNrOfBuckets + 1);
    size_t runningOffset = 0;
    for (size_t b = 0; b < this->totalNrOfBuckets; ++b) {
      this->bucketOffsets[b] = runningOffset;
      runningOffset += static_cast<size_t>(
        this->neighbourBucketSizes[static_cast<Eigen::Index>(b)]);
    }
    this->bucketOffsets[this->totalNrOfBuckets] = runningOffset;

    this->flatBuckets.resize(static_cast<size_t>(numAtoms));
    this->neighbourBucketSizes.setZero();
    for (Eigen::Index i = 0; i < numAtoms; ++i) {
      const bucket_idx_t bucketIndex =
        this->atomBucketsScratch[static_cast<size_t>(i)];
      const size_t pos =
        this->bucketOffsets[static_cast<size_t>(bucketIndex)] +
        static_cast<size_t>(this->neighbourBucketSizes[bucketIndex]);
      this->flatBuckets[pos] = i;
      this->neighbourBucketSizes[bucketIndex]++;
    }
  }

  bool EigenNeighbourList::checkIfCoordinatesAreCurrent(
    const Eigen::Ref<const Eigen::VectorXd>& newCoordinates) const
  {
    if (newCoordinates.size() % 3 != 0) {
      return false;
    }
    const Eigen::Index numAtoms = newCoordinates.size() / 3;
    if (static_cast<size_t>(numAtoms) != this->flatBuckets.size()) {
      return false;
    }

    Eigen::ArrayXi seenInBucket =
      Eigen::ArrayXi::Zero(static_cast<Eigen::Index>(this->totalNrOfBuckets));
    for (Eigen::Index i = 0; i < numAtoms; ++i) {
      const bucket_idx_t bucketIndex = this->getBucketIndexForTriplet(
        this->getBucketTripletForCoordinates(newCoordinates.segment<3>(3 * i)));
      if (bucketIndex < 0 ||
          static_cast<size_t>(bucketIndex) >= this->totalNrOfBuckets) {
        return false;
      }
      const int count = seenInBucket[bucketIndex];
      if (count >= this->neighbourBucketSizes[bucketIndex]) {
        return false;
      }
      const size_t pos = this->bucketOffsets[static_cast<size_t>(bucketIndex)] +
                         static_cast<size_t>(count);
      if (this->flatBuckets[pos] != i) {
        return false;
      }
      seenInBucket[bucketIndex] = count + 1;
    }
    return true;
  }

  void EigenNeighbourList::validateWhyNotIncluded(
    const Eigen::Ref<const Eigen::Vector3d>& sourceCoords,
    const Eigen::Ref<const Eigen::Vector3d>& targetCoords,
    double newCutoff) const
  {
    if (newCutoff <= 0.) {
      newCutoff = this->cutoff;
    }
    Eigen::Array3li baseTriplet =
      this->getBucketTripletForCoordinates(sourceCoords);
    bucket_idx_t baseBucketIndex = this->getBucketIndexForTriplet(baseTriplet);

    Eigen::Array3li targetTriplet =
      this->getBucketTripletForCoordinates(targetCoords);
    bucket_idx_t targetBucketIndex =
      this->getBucketIndexForTriplet(targetTriplet);

    bool foundSource = false;
    bool foundTarget = false;
    std::vector<bucket_idx_t> bucketIdxs;
    if (APPROX_EQUAL(newCutoff, this->cutoff, 1e-12)) {
      const size_t startOffset =
        static_cast<size_t>(baseBucketIndex) * this->defaultCutoffStride;
      const size_t endOffset = startOffset + this->defaultCutoffStride;
      bucketIdxs.assign(this->neighbourBucketNeighboursDefaultCutoff.begin() +
                          static_cast<std::ptrdiff_t>(startOffset),
                        this->neighbourBucketNeighboursDefaultCutoff.begin() +
                          static_cast<std::ptrdiff_t>(endOffset));
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
    const Eigen::Ref<const Eigen::Vector3d>& coordinates,
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
    return static_cast<long int>(this->flatBuckets.size());
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
    const Eigen::Ref<const Eigen::VectorXd>& coordinates,
    const int source,
    const double newCutoff) const
  {
    const Eigen::Vector3d sourceCoords = coordinates.segment<3>(3 * source);
    int numNeighbours =
      this->getIndicesCloseToCoordinates(result, sourceCoords, newCutoff);
    int actualNumNeighbours = 0;
    const double cutoff2 = newCutoff * newCutoff;
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(numNeighbours);
         ++i) {
      if (result[i] > source) {
        Eigen::Vector3d dist =
          sourceCoords - coordinates.segment<3>(result[i] * 3);
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
    const Eigen::Ref<const Eigen::Vector3d>& coordinates,
    const double upperCutoff,
    const bool expectDefault) const
  {
    if (this->totalNrOfBuckets < 2) {
      const Eigen::Index nResults =
        static_cast<Eigen::Index>(this->flatBuckets.size());
      if (result.size() < nResults) {
        result.conservativeResize(nResults);
      }
      if (nResults > 0) {
        result.head(nResults) = Eigen::ArrayXi::LinSpaced(
          nResults, 0, static_cast<int>(nResults - 1));
      }
      return static_cast<int>(nResults);
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
      const size_t startOffset =
        static_cast<size_t>(coordinatesBucketIdx) * this->defaultCutoffStride;
      const size_t endOffset = startOffset + this->defaultCutoffStride;
      // first, count the number of results we will get
      long int nResults = 0;
      for (size_t offset = startOffset; offset < endOffset; ++offset) {
        nResults += this->neighbourBucketSizes
                      [this->neighbourBucketNeighboursDefaultCutoff[offset]];
      }

      if (result.size() < nResults) {
        // heuristic minimum
        result.conservativeResize(nResults);
      }
      for (size_t offset = startOffset; offset < endOffset; ++offset) {
        const bucket_idx_t bucketIndex =
          this->neighbourBucketNeighboursDefaultCutoff[offset];
        const size_t bucketStart =
          this->bucketOffsets[static_cast<size_t>(bucketIndex)];
        const int bucketSize = this->neighbourBucketSizes[bucketIndex];
        for (int indexInBucket = 0; indexInBucket < bucketSize;
             ++indexInBucket) {
          const coordinate_idx_t atomIndex =
            this->flatBuckets[bucketStart + static_cast<size_t>(indexInBucket)];
          result[results_idx] = static_cast<int>(atomIndex);
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
            if (bucketIndex >= 0 &&
                !visitedBuckets[static_cast<size_t>(bucketIndex)]) {
              visitedBuckets[static_cast<size_t>(bucketIndex)] = true;
              const size_t bucketStart =
                this->bucketOffsets[static_cast<size_t>(bucketIndex)];
              const int bucketSize = this->neighbourBucketSizes[bucketIndex];
              if (results_idx + bucketSize > result.size()) {
                result.conservativeResize(std::max<Eigen::Index>(
                  result.size() * 2,
                  static_cast<Eigen::Index>(results_idx + bucketSize)));
              }
              for (int indexInBucket = 0; indexInBucket < bucketSize;
                   ++indexInBucket) {
                const coordinate_idx_t atomIndex =
                  this->flatBuckets[bucketStart +
                                    static_cast<size_t>(indexInBucket)];
                result[results_idx] = static_cast<int>(atomIndex);
                results_idx += 1;
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
    std::vector<std::vector<coordinate_idx_t>> result(this->totalNrOfBuckets);
    for (size_t b = 0; b < this->totalNrOfBuckets; ++b) {
      const size_t start = this->bucketOffsets[b];
      const size_t end =
        start + static_cast<size_t>(
                  this->neighbourBucketSizes[static_cast<Eigen::Index>(b)]);
      result[b].assign(
        this->flatBuckets.begin() + static_cast<std::ptrdiff_t>(start),
        this->flatBuckets.begin() + static_cast<std::ptrdiff_t>(end));
    }
    return result;
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
    bucket_idx_t verifyBucketIdx = this->getBucketIndexForTriplet(coeffs);
    assert(verifyBucketIdx == bucketIndex);
    assert(verifyCoeffs.isApprox(coeffs));
#endif
    return results;
  }

  std::vector<bucket_idx_t>
  EigenNeighbourList::getCombinedBucketIndicesForCoordinates(
    const Eigen::Ref<const Eigen::Vector3d>& coordinates,
    const double newCutoff,
    const bool sort) const
  {
    Eigen::Array3li maxIndices = this->getBucketTripletForCoordinates(
      coordinates + Eigen::Vector3d::Constant(newCutoff));
    Eigen::Array3li minIndices = this->getBucketTripletForCoordinates(
      coordinates - Eigen::Vector3d::Constant(newCutoff));

    // use a set to avoid returning duplicates
    std::set<bucket_idx_t> result;
    Eigen::Array3li indexTriplet;
    // now, do permutations of all these
    for (long int i = minIndices[0]; i <= maxIndices[0]; ++i) {
      indexTriplet[0] = i;
      for (long int j = minIndices[1]; j <= maxIndices[1]; ++j) {
        indexTriplet[1] = j;
        for (long int k = minIndices[2]; k <= maxIndices[2]; ++k) {
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
    bucketIndex += static_cast<long int>(newNrOfBuckets) *
                   static_cast<long int>(bucketIndex < 0);
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
      index / (this->nrOfBuckets[0] * this->nrOfBuckets[1]);
    bucket_idx_t bucketIndexY =
      (index - bucketIndexZ * (this->nrOfBuckets[0] * this->nrOfBuckets[1])) /
      (this->nrOfBuckets[0]);
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
    const Eigen::Ref<const Eigen::Vector3d>& coordinates) const
  {
    return ((coordinates.array() + this->box.getLowL()) / this->bucketWidths)
      .floor()
      .cast<long int>();
  }
};
};
