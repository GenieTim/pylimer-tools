#ifndef MEHP_UTILITY_STRUCT_H
#define MEHP_UTILITY_STRUCT_H

#include <Eigen/Dense>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>
#include <set>
#include "../utils/utilityMacros.h"

namespace pylimer_tools {
namespace calc {
  namespace mehp {
    enum ExitReason
    {
      UNSET,
      F_TOLERANCE,
      X_TOLERANCE,
      MAX_STEPS,
      FAILURE,
      OTHER
    };

    typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
    // typedef Eigen::Array<Eigen::ArrayXi, Eigen::Dynamic, 1> ArrayXArrayXi;
    typedef std::vector<std::vector<size_t>> ArrayXArrayXi;
    typedef std::vector<std::set<size_t>> ArrayXArrayXiUnique;
    typedef std::vector<std::vector<double>> ArrayXArrayXd;

    // improved structures using Eigen
    struct Network
    {
      double L[3];                    /* box sizes */
      double vol;                     /* box volume */
      double meanSpringContourLength; /* mean N */
      long int nrOfNodes;             /* number of nodes */
      long int nrOfSprings;           /* number of springs */
      long int nrOfLoops;             /* loops */
      // coordinates & connectivity
      Eigen::VectorXd coordinates;
      Eigen::VectorXd springsContourLength; /* the N for each spring */
      Eigen::ArrayXi oldAtomIds;
      Eigen::ArrayXi springCoordinateIndexA;
      Eigen::ArrayXi springCoordinateIndexB;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
      // interesting properties
      ArrayXb springIsActive;
    };

    struct ForceBalanceNetwork
    {
      // TODO: some info is redundant.
      // adjust code to support one way of storing things only
      double L[3];                    /* box sizes */
      double boxHalfs[3];             /* half box sizes */
      double vol;                     /* box volume */
      double meanSpringContourLength; /* mean N */
      long int nrOfLinks; /* number of links, = nrOfNodes + nrOfSlipLinks */
      long int nrOfNodes; /* number of cross-links */
      long int nrOfSprings;
      long int nrOfPartialSprings;
      long int nrOfSpringsWithPartition;
      // coordinates & connectivity
      Eigen::VectorXd coordinates;
      Eigen::VectorXd springsContourLength; /* the N for each spring */
      ArrayXArrayXi springIndicesOfLinks;   // maps link -> springs
      ArrayXArrayXi linkIndicesOfSprings;   // maps spring -> links
      ArrayXb partialSpringIsPartial; // indicates whether a spring involves a slip-link
      // local to global: from the 2D structures to the 1D Eigen vector
      ArrayXArrayXi localToGlobalSpringIndex;
      // map the "local", partial, spring indices to the full-length springs
      Eigen::ArrayXi partialToFullSpringIndex;
      std::unordered_map<size_t, size_t> oldAtomIdToSpringIndex;

      ArrayXb linkIsSliplink;
      Eigen::ArrayXi springPartCoordinateIndexA;
      Eigen::ArrayXi springPartCoordinateIndexB;
      Eigen::ArrayXi springPartIndexA;
      Eigen::ArrayXi springPartIndexB;

      // old stuff used for conversion. Does not include slip-links
      Eigen::ArrayXi springCoordinateIndexA;
      Eigen::ArrayXi springCoordinateIndexB;
      Eigen::ArrayXi oldAtomIds;
      std::vector<size_t> springToMoleculeIds; // maps 
      ArrayXb springIsActive;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
    };
  } // mehp
} // calc
} // pylimer_tools

#endif
