#ifndef MEHP_UTILITY_STRUCT_H
#define MEHP_UTILITY_STRUCT_H

#include <Eigen/Dense>
#include <array>
#include <math.h> // fabs, log, copysign, fma

namespace pylimer_tools {
namespace calc {
  namespace mehp {

    typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

    // improved structures using Eigen
    struct Network
    {
      double L[3];                /* box sizes */
      double vol;                 /* box volume */
      long int nrOfNodes;         /* number of nodes */
      long int nrOfSprings;       /* number of springs */
      double averageSpringLength; /* average spring length */
      long int nrOfLoops;         /* loops */
      // coordinates & coonectivity
      Eigen::VectorXd coordinates;
      Eigen::ArrayXi oldAtomIds;
      Eigen::ArrayXi springCoordinateIndexA;
      Eigen::ArrayXi springCoordinateIndexB;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
      // interesting properties
      ArrayXb springIsActive;
    };
  }
}
}

#endif
