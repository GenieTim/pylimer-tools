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
namespace sim::mehp {
#define STRUCTURE_SIMPLIFICATION_MODES                                         \
  X(NO_SIMPLIFICATION, "No Simplification")                                    \
  X(X2F_ONLY, "Two-function cross-links only")                                 \
  X(INACTIVE_ONLY, "Inactive links only")                                      \
  X(ALL_TIM, "All, à la Tim")                                                  \
  X(ALL_ANDREI, "All, à la Andrei")

  enum StructureSimplificationMode
  {
#define X(e, s) e,
    STRUCTURE_SIMPLIFICATION_MODES
#undef X
  };

  static std::string StructureSimplificationModeNames[] = {
#define X(e, s) s,
    STRUCTURE_SIMPLIFICATION_MODES
#undef X
  };

#define LINK_SWAPPING_MODES                                                    \
  X(NO_SWAPPING, "No Swapping")                                                \
  X(SLIPLINKS_ONLY, "Sliplinks only")                                          \
  X(ALL, "All")                                                                \
  X(ALL_CYCLE, "All, restrict to cycles")                                      \
  X(ALL_MC, "All MC")                                                          \
  X(ALL_MC_CYCLE, "All MC, restrict to cycles")                                \
  X(ALL_MC_TRY, "All MC, attempt the move")                                    \
  X(ALL_MC_TRY_CYCLE, "All MC, restrict to cycles, attempt the move")

  enum LinkSwappingMode
  {
#define X(e, s) e,
    LINK_SWAPPING_MODES
#undef X
  };

  static std::string LinkSwappingModeNames[] = {
#define X(e, s) s,
    LINK_SWAPPING_MODES
#undef X
  };

#define EXIT_REASONS                                                           \
  X(UNSET, "Unset")                                                            \
  X(F_TOLERANCE, "F Tolerance")                                                \
  X(X_TOLERANCE, "X Tolerance")                                                \
  X(MAX_STEPS, "Max Steps")                                                    \
  X(NO_STEPS_POSSIBLE, "No Steps Possible")                                    \
  X(FAILURE, "Failure")                                                        \
  X(INTERRUPT, "Interrupt")                                                    \
  X(OTHER, "Other")

  enum ExitReason
  {
#define X(e, s) e,
    EXIT_REASONS
#undef X
  };

  static std::string ExitReasonNames[] = {
#define X(e, s) s,
    EXIT_REASONS
#undef X
  };

#define SLE_SOLVERS                                                            \
  X(DEFAULT, "default")                                                        \
  /* direct methods */                                                         \
  X(SIMPLICIAL_LLT, "SimplicialLLT")                                           \
  X(SIMPLICIAL_LDLT, "SimplicialLDLT")                                         \
  X(SPARSE_LU, "SparseLU")                                                     \
  X(SPARSE_QR, "SparseQR")                                                     \
  /* iterative methods */                                                      \
  X(CONJUGATE_GRADIENT, "ConjugateGradient")                                   \
  X(CONJUGATE_GRADIENT_DIAGONALIZED,                                           \
    "ConjugateGradient, DiagonalPreconditioner")                               \
  X(CONJUGATE_GRADIENT_IDENTITY, "ConjugateGradient, IdentityPreconditioner")  \
  X(CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY,                                    \
    "ConjugateGradient, IncompleteCholeskyPreconditioner")                     \
  X(LEAST_SQUARES_CONJUGATE_GRADIENT, "LeastSquaresConjugateGradient")         \
  X(LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED,                             \
    "LeastSquaresConjugateGradient, DiagonalPreconditioner")                   \
  X(LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY,                                 \
    "LeastSquaresConjugateGradient, IdentityPreconditioner")                   \
  X(BICGSTAB, "BiCGSTAB")                                                      \
  X(BICGSTAB_DIAGONALIZED, "BiCGSTAB, DiagonalPreconditioner")                 \
  X(BICGSTAB_IDENTITY, "BiCGSTAB, IdentityPreconditioner")                     \
  X(BICGSTAB_INCOMPLETE_LU, "BiCGSTAB, IncompleteLUTPreconditioner")

  enum SLESolver
  {
#define X(e, s) e,
    SLE_SOLVERS
#undef X
  };

  static std::string SLESolverNames[] = {
#define X(e, s) s,
    SLE_SOLVERS
#undef X
  };

  static SLESolver allSLESolvers[] = {
#define X(e, s) SLESolver::e,
    SLE_SOLVERS
#undef X
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

  struct ForceBalance2Network
  {
    // TODO: some info is redundant.
    // adjust code to support one way of storing things only
    double L[3];          /* box sizes */
    double boxHalfs[3];   /* half box sizes */
    size_t nrOfLinks = 0; /* number of links, = nrOfNodes + nrOfSlipLinks */
    size_t nrOfNodes = 0; /* number of crosslinkers */
    size_t nrOfStrands = 0;
    size_t nrOfSprings = 0;
    // coordinates & connectivity
    Eigen::VectorXd coordinates;
    Eigen::VectorXd springContourLength; // the N for each spring
    Eigen::ArrayXb
      springIsEntanglement; // Needed for entanglements modelled as springs
    ArrayXArrayXi strandIndicesOfLink; // maps link -> strands
    ArrayXArrayXi linkIndicesOfStrand; // maps strands -> links
    ArrayXArrayXi springIndicesOfStrand;
    // map the "local", partial, spring indices to the full-length springs
    Eigen::ArrayXb linkIsEntanglement; // whether a link is a "cross-link" or an
    // "entanglement-link"

    // partial springs
    Eigen::ArrayXi springCoordinateIndexA;
    Eigen::ArrayXi springCoordinateIndexB;
    Eigen::ArrayXi springIndexA;
    Eigen::ArrayXi springIndexB;
    Eigen::VectorXd springBoxOffset;    // the "PBC" for each spring
    Eigen::ArrayXi strandIndexOfSpring; // the mapping from spring to strand

    Eigen::ArrayXi oldAtomIds;
    Eigen::ArrayXi oldAtomTypes;
  };
}
} // pylimer_tools

#endif
