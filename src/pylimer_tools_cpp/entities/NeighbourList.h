#pragma once

#include "Atom.h"
#include "Box.h"
#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

#include <algorithm>

namespace pylimer_tools {
namespace entities {

  /**
   * @brief A simple neighbour list implementation using binning
   *
   */
  class NeighbourList
  {
  public:
    NeighbourList(const std::vector<Atom>& atomList,
                  const Box& boxGeometry,
                  double cutoffDistance);

    std::vector<pylimer_tools::entities::Atom> getAtomsCloseTo(
      const pylimer_tools::entities::Atom& atom);

    std::vector<pylimer_tools::entities::Atom> getAtomsCloseTo(
      const pylimer_tools::entities::Atom& atom,
      double upperCutoff,
      double lowerCutoff = 0.0,
      bool unwrapped = false,
      bool expectSelf = false);

    void removeAtom(const Atom& atom, const std::string debugHint = "");

  protected:
    size_t normalizeBucketIndex(long int bucketIndex, size_t nrOfBuckets) const;

    size_t getBucketIndexForTriplet(
      std::tuple<long int, long int, long int> ind) const;

    std::tuple<long int, long int, long int> getBucketIndicesForAtom(
      const pylimer_tools::entities::Atom& atom) const;

    std::vector<size_t> getCombinedBucketIndicesForAtom(
      const pylimer_tools::entities::Atom& atom,
      double newCutoff) const;

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
