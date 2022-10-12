#ifndef MEHP_UTILITY_STRUCT_H
#define MEHP_UTILITY_STRUCT_H

#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <map>
#include <utility>

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
      double L[3];                    /* box sizes */
      double boxHalfs[3];             /* half box sizes */
      double vol;                     /* box volume */
      double meanSpringContourLength; /* mean N */
      long int nrOfLinks; /* number of links, = nrOfNodes + nrOfSlipLinks */
      long int nrOfNodes; /* number of cross-links */
      long int nrOfSprings;
      long int nrOfPartialSprings;
      // coordinates & connectivity
      Eigen::VectorXd coordinates;
      Eigen::VectorXd springsContourLength; /* the N for each spring */
      ArrayXArrayXi springIndicesOfLinks;   // maps link -> springs
      ArrayXArrayXi linkIndicesOfSprings;   // maps spring -> links
      // TODO: this will fail if we have two springs between the same links
      std::map<std::pair<size_t, size_t>, size_t> connectivityToSpringIndex;
      

      ArrayXb linkIsSliplink;
      Eigen::ArrayXi springPartCoordinateIndexA;
      Eigen::ArrayXi springPartCoordinateIndexB;
      Eigen::ArrayXi springPartIndexA;
      Eigen::ArrayXi springPartIndexB;
      // old stuff used for conversion. Does not include slip-links
      Eigen::ArrayXi springCoordinateIndexA;
      Eigen::ArrayXi springCoordinateIndexB;
      Eigen::ArrayXi oldAtomIds;
      ArrayXb springIsActive;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
    };
    // to string, without macro expansion
#define STRINGINFY(s) #s
    // to string, with macro expansion
#define XSTRINGINFY(s) STRINGINFY(s)

#define INVALIDARG_EXP_IFN(condition, message)                                 \
  if (!(condition)) {                                                          \
    throw std::invalid_argument(message "\nFailed condition: " #condition);    \
  }

#define RUNTIME_EXP_IFN(condition, message)                                    \
  if (!(condition)) {                                                          \
    throw std::runtime_error(message "\nFailed condition: " #condition);       \
  }
  } // mehp

#define APPROX_EQUAL(value1, value2, eps) \
  value1 + eps >= value2 && value1 - eps <= value2
} // calc
} // pylimer_tools

#endif
