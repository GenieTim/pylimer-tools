#ifndef EIGENNEIGHBOURLIST_H
#define EIGENNEIGHBOURLIST_H

extern "C"
{
#include <igraph/igraph.h>
}
#include "../utils/ExtraEigenTypes.h"
#include "Box.h"
#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <algorithm>

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
                       const Box& box,
                       double cutoff,
                       double scalingFactor = 2.);

    void initialize(const Eigen::VectorXd& coordinates,
                    const Box& box,
                    double cutoff,
                    double scalingFactor = 2.);

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
                                     const double upperCutoff,
                                     bool expectDefault = false) const;

    int getHigherIndicesWithinCutoff(Eigen::ArrayXi& result,
                                     const Eigen::VectorXd& coordinates,
                                     const int source,
                                     const double cutoff) const;

    /**
     * @brief For debugging/test purposes: the actual buckets
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

    template<class Archive>
    void serialize(Archive& ar)
    {
        ar(bucketWidths,
           nrOfBuckets,
           totalNrOfBuckets,
           cutoff,
           scalingFactor,
           box,
           neighbourBuckets,
           neighbourBucketNeighboursDefaultCutoff,
           neighbourBucketSizes);
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
                                      size_t nrOfBuckets) const;

    void normalizeTriplet(Eigen::Array3li& triplet) const;

    Eigen::Array3li tripletFromIndex(bucket_idx_t index) const;

    bucket_idx_t getBucketIndexForTriplet(Eigen::Array3li ind) const;

    Eigen::Array3li getBucketTripletForCoordinates(
        const Eigen::Vector3d& coordinates) const;

private:
    Eigen::Array3d bucketWidths;
    Eigen::Array3li nrOfBuckets;

    size_t totalNrOfBuckets;

    double cutoff;
    double scalingFactor;
    double actualCutoff;

    pylimer_tools::entities::Box box;

    std::vector<std::vector<coordinate_idx_t>> neighbourBuckets;
    std::vector<std::vector<bucket_idx_t>>
    neighbourBucketNeighboursDefaultCutoff;
    Eigen::VectorXi neighbourBucketSizes;
};
};
};

#endif
