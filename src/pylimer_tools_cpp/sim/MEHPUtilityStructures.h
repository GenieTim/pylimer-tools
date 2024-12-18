#ifndef MEHP_UTILITY_STRUCT_H
#define MEHP_UTILITY_STRUCT_H

#include "../utils/ExtraEigenTypes.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pylimer_tools {
namespace sim {
  namespace mehp {
    enum StructureSimplificationMode
    {
      NO_SIMPLIFICATION,
      X2F_ONLY,
      INACTIVE_ONLY,
      ALL_TIM,
      ALL_ANDREI
    };

    enum LinkSwappingMode
    {
      NO_SWAPPING,
      SLIPLINKS_ONLY,
      ALL,
      ALL_CYCLE,
      ALL_MC,
      ALL_MC_CYCLE,
      ALL_MC_TRY,
      ALL_MC_TRY_CYCLE,
    };

    enum ExitReason
    {
      UNSET,
      F_TOLERANCE,
      X_TOLERANCE,
      MAX_STEPS,
      NO_STEPS_POSSIBLE,
      FAILURE,
      INTERRUPT,
      OTHER
    };

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
      size_t nrOfNodes = 0;           /* number of nodes */
      size_t nrOfSprings = 0;         /* number of springs */
      size_t nrOfLoops = 0;           /* loops */
      // coordinates & connectivity
      Eigen::VectorXd coordinates;
      Eigen::VectorXd springsContourLength; /* the N for each spring */
      Eigen::ArrayXi oldAtomIds;
      Eigen::ArrayXi springCoordinateIndexA;
      Eigen::ArrayXi springCoordinateIndexB;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
      Eigen::VectorXd springBoxOffset;

      ArrayXArrayXi springIndicesOfLinks; // maps link -> springs
      // interesting properties
      Eigen::ArrayXb springIsActive;
      Eigen::ArrayXi moleculeIdxToSpring;

      // config
      bool assumeBoxLargeEnough = false;
    };

    struct ForceBalanceNetwork
    {
      // TODO: some info is redundant.
      // adjust code to support one way of storing things only
      double L[3];                          /* box sizes */
      double boxHalfs[3];                   /* half box sizes */
      double vol = 0.0;                     /* box volume */
      double meanSpringContourLength = 0.0; /* mean N */
      size_t nrOfLinks = 0; /* number of links, = nrOfNodes + nrOfSlipLinks */
      size_t nrOfNodes = 0; /* number of crosslinkers */
      size_t nrOfSprings = 0;
      size_t nrOfPartialSprings = 0;
      size_t nrOfSpringsWithPartition = 0;
      bool isUpToDate = true;
      // coordinates & connectivity
      Eigen::VectorXd coordinates;
      Eigen::VectorXd springsContourLength; /* the N for each spring */
      Eigen::ArrayXi springsType; // gives each spring a type. Needed for
                                  // entanglements modelled as springs
      ArrayXArrayXi springIndicesOfLinks;    // maps link -> springs
      ArrayXArrayXi linkIndicesOfSprings;    // maps spring -> links
      Eigen::ArrayXb partialSpringIsPartial; // indicates whether a spring
      // involves a slip-link
      // local to global: from the 2D structures to the 1D Eigen vector
      // equivalent to "partial spring indices of spring"
      ArrayXArrayXi localToGlobalSpringIndex;
      // map the "local", partial, spring indices to the full-length springs
      std::unordered_map<size_t, size_t> oldAtomIdToSpringIndex;

      Eigen::ArrayXb linkIsSliplink;
      Eigen::ArrayXi
        nrOfCrosslinkSwapsEndured; // count for slip-links how many crosslinkers
      // they swapped around

      // partial springs
      Eigen::ArrayXi springPartCoordinateIndexA;
      Eigen::ArrayXi springPartCoordinateIndexB;
      Eigen::ArrayXi springPartIndexA;
      Eigen::ArrayXi springPartIndexB;
      Eigen::VectorXd springPartBoxOffset;
      Eigen::ArrayXi partialToFullSpringIndex;

      // these may be empty, or not, depending on the method used
      // to determine the slip-links
      ArrayXArrayXi loops;           // each loop just records its spring idx
      ArrayXArrayXi loopsOfSliplink; // each slip-link has two loops, ideally

      // old stuff used for conversion. Does not include slip-links
      Eigen::ArrayXi springCoordinateIndexA;
      Eigen::ArrayXi springCoordinateIndexB;
      Eigen::ArrayXi oldAtomIds;
      Eigen::ArrayXi oldAtomTypes;
      std::vector<size_t> springToMoleculeIds; // maps
      Eigen::ArrayXb springIsActive;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
    };
  } // mehp
} // calc
} // pylimer_tools

#endif
