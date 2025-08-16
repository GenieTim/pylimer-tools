#pragma once

#include "../utils/ExtraEigenTypes.h"
#include "Box.h"
#include <Eigen/Dense>
#include <vector>

namespace pylimer_tools {
namespace entities {

  typedef long int bucket_idx_t;
  typedef long int coordinate_idx_t;

  /**
   * @brief An implementation of a neighbour list using binning, using Eigen for
   * performance
   *
   */
  class EigenNeighbourList
  {
  public:
    EigenNeighbourList() {}
    EigenNeighbourList(const Eigen::VectorXd& coordinates,
                       const Box& boxGeometry,
                       double cutoffDistance,
                       double scalingParam = 1.);

    /**
     * @brief Equality operator to compare two EigenNeighbourList instances
     *
     * @param other The other EigenNeighbourList to compare with
     * @return bool True if the two instances are equal
     */
    bool operator==(const EigenNeighbourList& other) const;

    /**
     * @brief Inequality operator to compare two EigenNeighbourList instances
     *
     * @param other The other EigenNeighbourList to compare with
     * @return bool True if the two instances are not equal
     */
    bool operator!=(const EigenNeighbourList& other) const;

    void initialize(const Eigen::VectorXd& coordinates,
                    const Box& box,
                    double cutoff,
                    double scalingFactor = 1.);

    /**
     * @brief Re-bin with a new set of coordinates
     *
     * @param newCoordinates
     */
    void resetCoordinates(Eigen::VectorXd& newCoordinates);

    bool checkIfCoordinatesAreCurrent(Eigen::VectorXd& newCoordinates);

    void validateWhyNotIncluded(Eigen::Vector3d sourceCoords,
                                Eigen::Vector3d targetCoords,
                                double newCutoff = -1.0) const;

    /**
     * @brief Get the Indices Close To Coordinates with the Default Cut-Off
     *
     * @param coordinates
     * @return Eigen::ArrayXi
     */
    Eigen::ArrayXi getIndicesCloseToCoordinates(Eigen::Vector3d coordinates,
                                                double newCutoff = -1.) const;

    long int getNumBinnedCoordinates() const;

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
    int getIndicesCloseToCoordinates(Eigen::ArrayXi& result,
                                     const Eigen::Vector3d coordinates,
                                     const double upperCutoff,
                                     bool expectDefault = false) const;

    int getHigherIndicesWithinCutoff(Eigen::ArrayXi& result,
                                     const Eigen::VectorXd& coordinates,
                                     const int source,
                                     const double cutoff) const;

    /**
     * @brief For debugging/test purposes: The actual buckets
     *
     * @return std::vector<std::vector<coordinate_idx_t>>
     */
    std::vector<std::vector<coordinate_idx_t>> getNeighbourBuckets() const;
    Eigen::VectorXi getNeighbourBucketSizes() const;
    Eigen::Vector3d getCentralCoordinatesOfBucket(int bucketIndex) const;

    std::vector<bucket_idx_t> getCombinedBucketIndicesForCoordinates(
      const Eigen::Vector3d& coordinates,
      double newCutoff,
      bool sort = false) const;

#ifdef CEREALIZABLE
    template<class Archive>
    void serialize(Archive& ar)
    {
      ar(bucketWidths,
         nrOfBuckets,
         totalNrOfBuckets,
         cutoff,
         scalingFactor,
         actualCutoff,
         box,
         neighbourBuckets,
         neighbourBucketNeighboursDefaultCutoff,
         neighbourBucketSizes);
    }
#endif

  protected:
    /**
     * @brief Do "PBC" with a bucket index in one direction
     *
     * @param bucketIndex
     * @param nrOfBuckets
     * @return size_t
     */
    bucket_idx_t normalizeBucketIndex(long int bucketIndex,
                                      size_t nrOfBuckets) const;

    void normalizeTriplet(Eigen::Array3li& triplet) const;

    /**
     * @brief Decompose the one index into a triplet of indices
     *
     * @param index
     * @return Eigen::Array3li
     */
    Eigen::Array3li tripletFromIndex(bucket_idx_t index) const;

    /**
     * @brief Get one index, normalized, given a triplet of indices
     *
     * @param ind
     * @return bucket_idx_t
     */
    bucket_idx_t getBucketIndexForTriplet(Eigen::Array3li ind) const;

    Eigen::Array3li getBucketTripletForCoordinates(
      const Eigen::Vector3d& coordinates) const;

  private:
    // the width of the buckets in each dimension (double)
    Eigen::Array3d bucketWidths;
    // the number of buckets in each dimension (integer)
    Eigen::Array3li nrOfBuckets;

    // the total number of buckets (integer)
    size_t totalNrOfBuckets = 0;

    double cutoff = 1.0;
    double scalingFactor = 1.0;
    double actualCutoff = 1.0;

    pylimer_tools::entities::Box box;

    std::vector<std::vector<coordinate_idx_t>> neighbourBuckets;
    std::vector<std::vector<bucket_idx_t>>
      neighbourBucketNeighboursDefaultCutoff;
    Eigen::VectorXi neighbourBucketSizes;
  };
};
};
