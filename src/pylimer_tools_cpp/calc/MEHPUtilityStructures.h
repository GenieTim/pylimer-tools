#ifndef MEHP_UTILITY_STRUCT_H
#define MEHP_UTILITY_STRUCT_H

#include <Eigen/Dense>
#include <array>
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {

    typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
    // typedef Eigen::Array<Eigen::ArrayXi, Eigen::Dynamic, 1> ArrayXArrayXi;
    typedef std::vector<std::vector<size_t>> ArrayXArrayXi;

    // improved structures using Eigen
    struct Network
    {
      double L[3];          /* box sizes */
      double vol;           /* box volume */
      long int nrOfNodes;   /* number of nodes */
      long int nrOfSprings; /* number of springs */
      long int nrOfLoops;   /* loops */
      // coordinates & connectivity
      Eigen::VectorXd coordinates;
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
      double L[3]; /* box sizes */
      double boxHalfs[3]; /* half box sizes */
      double vol;         /* box volume */
      long int nrOfLinks; /* number of links, = nrOfNodes + nrOfSlipLinks */
      long int nrOfNodes; /* number of nodes */
      long int nrOfSprings;
      // coordinates & connectivity
      Eigen::VectorXd coordinates;
      ArrayXArrayXi springIndicesOfLinks; // maps link -> springs
      ArrayXArrayXi linkIndicesOfSprings; // maps spring -> links
      ArrayXb linkIsSpliplink;
      // old stuff used for conversion. Does not include slip-links
      Eigen::ArrayXi oldAtomIds;
      Eigen::ArrayXi springCoordinateIndexA;
      Eigen::ArrayXi springCoordinateIndexB;
      Eigen::ArrayXi springIndexA;
      Eigen::ArrayXi springIndexB;
      ArrayXb springIsActive;
    };
  }
}
}

#endif
