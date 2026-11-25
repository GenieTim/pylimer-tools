#ifndef MEHP_UTILITY_STRUCT_H
#define MEHP_UTILITY_STRUCT_H

#include "../utils/ExtraEigenTypes.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <array>
#include <cctype>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pylimer_tools::sim::mehp {
#define STRUCTURE_SIMPLIFICATION_MODES                                         \
  X(NO_SIMPLIFICATION, "No Simplification")                                    \
  X(X2F_ONLY, "Removal of two-functional crosslinks only")                     \
  X(INACTIVE_ONLY, "Removal of inactive entanglement- and crosslinks only")    \
  X(INACTIVE_THEN_X2F,                                                         \
    "Removal of inactive, and then two-functional entanglement- and "          \
    "crosslinks, "                                                             \
    "one after the other")                                                     \
  X(X1F_X2F_THEN_INACTIVE,                                                     \
    "Removal of one- and twofunctional crosslinks, and then inactive "         \
    "entanglements and crosslinks. Deprecated, use "                           \
    "INACTIVE_THEN_X2F instead")

enum StructureSimplificationMode
{
#define X(e, s) e,
  STRUCTURE_SIMPLIFICATION_MODES
#undef X
};

inline const std::string StructureSimplificationModeNames[] = {
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

inline const std::string LinkSwappingModeNames[] = {
#define X(e, s) s,
  LINK_SWAPPING_MODES
#undef X
};

#define EXIT_REASONS                                                           \
  X(UNSET, "Unset")                                                            \
  X(F_TOLERANCE, "F (force) tolerance reached")                                \
  X(X_TOLERANCE, "X (displacement) tolerance reached")                         \
  X(MAX_STEPS, "Maximum number of steps reached")                              \
  X(NO_STEPS_POSSIBLE, "No (more) steps possible")                             \
  X(FAILURE, "Failure")                                                        \
  X(INTERRUPT, "Interrupt")                                                    \
  X(OTHER, "Other")

enum ExitReason
{
#define X(e, s) e,
  EXIT_REASONS
#undef X
};

inline const std::string ExitReasonNames[] = {
#define X(e, s) s,
  EXIT_REASONS
#undef X
};

#define SLE_SOLVERS                                                            \
  X(DEFAULT, "default")                                                        \
  X(AUTO, "auto")                                                              \
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
  X(BICGSTAB_INCOMPLETE_LU, "BiCGSTAB, IncompleteLUTPreconditioner")           \
  X(GRADIENT_DESCENT, "GradientDescent")                                       \
  X(GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT,                                   \
    "GradientDescent (Barzilai-Borwein method, short time-step)")              \
  X(GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG,                                    \
    "GradientDescent (Barzilai-Borwein method, long time-step)")               \
  X(GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM,                                \
    "GradientDescent (Barzilai-Borwein & heavy ball method, selective "        \
    "time-step)")                                                              \
  X(METAL_CG, "Metal Conjugate Gradient (GPU)")                                \
  X(METAL_BARZILAI_BORWEIN, "Metal Barzilai-Borwein (GPU)")

enum SLESolver
{
#define X(e, s) e,
  SLE_SOLVERS
#undef X
};

inline const std::string SLESolverNames[] = {
#define X(e, s) s,
  SLE_SOLVERS
#undef X
};

inline const SLESolver allSLESolvers[] = {
#define X(e, s) SLESolver::e,
  SLE_SOLVERS
#undef X
};

/**
 * @brief Representation of topological entanglements in an MEHP network.
 */
enum class EntanglementMode
{
  NONE,         ///< Phantom network (crosslinks only)
  MERGED_LINKS, ///< Entanglement links lie directly on strand contour (no extra
                ///< spring)
  SPRING        ///< Entangled beads coupled by compliant Hookean spring
                ///< (springIsEntanglement[s] == true)
};

/**
 * @brief Mobility of entanglement links along polymer strands.
 */
enum class SlipMode
{
  FIXED,  ///< Entanglement links remain fixed at initial strand positions
  SLIDING ///< Entanglement links can slide along polymer strands
};

/**
 * @brief Supported NLopt optimization algorithms for MEHP force relaxation.
 */
enum class NLoptAlgorithm
{
  LD_MMA,
  LD_LBFGS,
  LD_SLSQP,
  LD_TNEWTON,
  LD_TNEWTON_RESTART,
  LD_TNEWTON_PRECOND,
  LD_TNEWTON_PRECOND_RESTART,
  LD_VAR1,
  LD_VAR2,
  LN_SBPLX,
  LN_BOBYQA,
  LN_NELDERMEAD,
  LN_COBYLA,
  LN_NEWUOA,
  LN_NEWUOA_BOUND,
  LN_PRAXIS,
  GN_DIRECT,
  GN_CRS2_LM,
  GD_STOGO
};

namespace detail {
  inline constexpr std::array<std::pair<NLoptAlgorithm, std::string_view>, 19>
    kNLoptAlgorithmNames = { {
      { NLoptAlgorithm::LD_MMA, "LD_MMA" },
      { NLoptAlgorithm::LD_LBFGS, "LD_LBFGS" },
      { NLoptAlgorithm::LD_SLSQP, "LD_SLSQP" },
      { NLoptAlgorithm::LD_TNEWTON, "LD_TNEWTON" },
      { NLoptAlgorithm::LD_TNEWTON_RESTART, "LD_TNEWTON_RESTART" },
      { NLoptAlgorithm::LD_TNEWTON_PRECOND, "LD_TNEWTON_PRECOND" },
      { NLoptAlgorithm::LD_TNEWTON_PRECOND_RESTART,
        "LD_TNEWTON_PRECOND_RESTART" },
      { NLoptAlgorithm::LD_VAR1, "LD_VAR1" },
      { NLoptAlgorithm::LD_VAR2, "LD_VAR2" },
      { NLoptAlgorithm::LN_SBPLX, "LN_SBPLX" },
      { NLoptAlgorithm::LN_BOBYQA, "LN_BOBYQA" },
      { NLoptAlgorithm::LN_NELDERMEAD, "LN_NELDERMEAD" },
      { NLoptAlgorithm::LN_COBYLA, "LN_COBYLA" },
      { NLoptAlgorithm::LN_NEWUOA, "LN_NEWUOA" },
      { NLoptAlgorithm::LN_NEWUOA_BOUND, "LN_NEWUOA_BOUND" },
      { NLoptAlgorithm::LN_PRAXIS, "LN_PRAXIS" },
      { NLoptAlgorithm::GN_DIRECT, "GN_DIRECT" },
      { NLoptAlgorithm::GN_CRS2_LM, "GN_CRS2_LM" },
      { NLoptAlgorithm::GD_STOGO, "GD_STOGO" },
    } };
} // namespace detail

/**
 * @brief Convert an @ref NLoptAlgorithm enumerator to its canonical string
 * name.
 *
 * @param alg The NLopt algorithm enumerator.
 * @return Canonical string representation (e.g. "LD_LBFGS").
 * @throws std::invalid_argument if @p alg is not a valid enumerator.
 */
[[nodiscard]] inline std::string
nloptAlgorithmToString(const NLoptAlgorithm alg)
{
  for (const std::pair<NLoptAlgorithm, std::string_view>& entry :
       detail::kNLoptAlgorithmNames) {
    if (entry.first == alg) {
      return std::string(entry.second);
    }
  }
  throw std::invalid_argument(
    "Unknown NLoptAlgorithm enumerator with integer value " +
    std::to_string(static_cast<int>(alg)) + ".");
}

/**
 * @brief Parse an @ref NLoptAlgorithm enumerator from a string name.
 *
 * Accepts canonical names (case-insensitive) with an optional "NLOPT_" prefix.
 *
 * @param name Algorithm name string.
 * @return Matching @ref NLoptAlgorithm enumerator.
 * @throws std::invalid_argument if @p name does not match any supported
 * algorithm.
 */
[[nodiscard]] inline NLoptAlgorithm
nloptAlgorithmFromString(const std::string_view name)
{
  std::string normalized;
  normalized.reserve(name.size());
  for (const char ch : name) {
    normalized.push_back(
      static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
  }
  constexpr std::string_view nloptPrefix = "NLOPT_";
  std::string_view lookup = normalized;
  if (lookup.starts_with(nloptPrefix)) {
    lookup.remove_prefix(nloptPrefix.size());
  }
  for (const std::pair<NLoptAlgorithm, std::string_view>& entry :
       detail::kNLoptAlgorithmNames) {
    if (entry.second == lookup) {
      return entry.first;
    }
  }
  throw std::invalid_argument("Unknown NLopt algorithm: '" + std::string(name) +
                              "'.");
}

/**
 * @brief Metadata describing a spring removed during bond/strand breaking.
 */
struct BrokenSpringInfo
{
  size_t idxA = 0;
  size_t idxB = 0;
  double multiplier = 0.0;

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(idxA, idxB, multiplier);
  }
#endif
};

// typedef Eigen::Array<Eigen::ArrayXi, Eigen::Dynamic, 1> ArrayXArrayXi;
typedef std::vector<std::vector<size_t>> ArrayXArrayXi;
typedef std::vector<std::set<size_t>> ArrayXArrayXiUnique;
typedef std::vector<std::vector<double>> ArrayXArrayXd;

struct CrosslinkerUniverse
{
  std::vector<int> xlinkTypes;
  std::vector<double> xlinkX;
  std::vector<double> xlinkY;
  std::vector<double> xlinkZ;
  std::vector<std::vector<long int>> strandsOfXlink;
  std::vector<int> xlinkChainId;

  std::vector<long int> strandFrom;
  std::vector<long int> strandTo;
  std::vector<int> beadsInStrand;
  std::vector<int> strandBeadType;
  std::vector<double> beadDistanceInStrand;
  std::vector<double> meanSquaredBeadDistanceInStrand;
  std::vector<uint64_t> strandSeeds;

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(xlinkTypes,
       xlinkX,
       xlinkY,
       xlinkZ,
       strandsOfXlink,
       xlinkChainId,
       strandFrom,
       strandTo,
       beadsInStrand,
       strandBeadType,
       beadDistanceInStrand,
       meanSquaredBeadDistanceInStrand,
       strandSeeds);
  }
#endif
};
}

#endif
