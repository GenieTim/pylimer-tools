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

#include <stdexcept>
#include <algorithm>

namespace pylimer_tools {
namespace entities {

  class NeighbourList
  {
  public:
    NeighbourList(const std::vector<Atom> &atoms, const Box &box, double cutoff);

    std::vector<pylimer_tools::entities::Atom> getAtomsCloseTo(
      pylimer_tools::entities::Atom &atom);

    std::vector<pylimer_tools::entities::Atom> getAtomsCloseTo(
      pylimer_tools::entities::Atom &atom,
      double upperCutoff,
      double lowerCutoff = 0.0,
      bool unwrapped = false);

    void removeAtom(Atom atom);

  protected:
    size_t normalizeBucketIndex(long int bucketIndex, size_t nrOfBuckets) const;

    size_t getBucketIndexForTriplet(
      std::tuple<long int, long int, long int> ind) const;

    std::tuple<size_t, size_t, size_t> getBucketIndicesForAtom(
      const pylimer_tools::entities::Atom& atom);

    std::vector<size_t> getCombinedBucketIndicesForAtom(
      const pylimer_tools::entities::Atom& atom,
      double newCutoff);

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
