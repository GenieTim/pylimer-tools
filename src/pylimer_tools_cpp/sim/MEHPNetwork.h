#ifndef MEHP_NETWORK_H
#define MEHP_NETWORK_H

#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../topo/EntanglementDetector.h"
#include "../utils/ExtraEigenTypes.h"
#include "../utils/utilityMacros.h"
#include "MEHPUtilityStructures.h"
#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef CEREALIZABLE
#include "../utils/CerealUtils.h"
#include <cereal/types/array.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#endif

namespace pylimer_tools::sim::mehp {

/**
 * @brief Canonical polymer network representation shared across MEHP force
 * balance, force relaxation, and slip-link solvers.
 *
 * Unifies crosslink-only phantom networks, merged-link entanglement networks,
 * spring-coupled entanglement networks, and sliding slip-link networks under a
 * single link/strand/spring layout. All links (crosslinks, chain ends, and
 * entanglement links) are indexed uniformly in `[0, nrOfLinks)`, with node
 * vs. entanglement identity distinguished strictly via `linkIsEntanglement`
 * (never via index threshold assumptions).
 */
struct MEHPNetwork
{
  // ===========================================================================
  // Box & dimensions
  // ===========================================================================
  std::array<double, 3> L = { 1., 1., 1. };
  std::array<double, 3> boxHalfs = { 0.5, 0.5, 0.5 };
  double vol = 1.0;
  pylimer_tools::entities::Box box = pylimer_tools::entities::Box(1., 1., 1.);

  // ===========================================================================
  // Counts & averages
  // ===========================================================================
  size_t nrOfLinks = 0; ///< Total links (nodes + entanglements)
  size_t nrOfNodes =
    0; ///< Count of non-entanglement links (!linkIsEntanglement)
  size_t nrOfStrands = 0; ///< Number of topological polymer/coupling strands
  size_t nrOfSprings = 0; ///< Number of sub-springs/springs
  double meanStrandContourLength = 0.0;
  double meanSpringContourLength = 0.0;

  // ===========================================================================
  // Topology flags
  // ===========================================================================
  bool assumeBoxLargeEnough = false;
  bool removePrimaryLoops = false;
  bool assumeComplete = false;

  // ===========================================================================
  // Link arrays (sized nrOfLinks or 3 * nrOfLinks)
  // ===========================================================================
  Eigen::VectorXd coordinates;
  Eigen::ArrayXb linkIsEntanglement;
  Eigen::ArrayXi oldAtomIds; ///< -1 for entanglement links
  Eigen::ArrayXi oldAtomTypes;
  ArrayXArrayXi strandIndicesOfLink;
  Eigen::ArrayXi
    nrOfCrosslinkSwapsEndured; ///< Always sized nrOfLinks, indexed by linkIdx

  // ===========================================================================
  // Strand arrays (sized nrOfStrands)
  // ===========================================================================
  Eigen::VectorXd strandContourLength; ///< Total N_s per strand; 1.0 for SPRING
                                       ///< coupling strands
  ArrayXArrayXi linkIndicesOfStrand;
  ArrayXArrayXi springIndicesOfStrand;

  // ===========================================================================
  // Spring arrays (sized nrOfSprings or 3 * nrOfSprings)
  // ===========================================================================
  Eigen::ArrayXi springIndexA;
  Eigen::ArrayXi springIndexB;
  Eigen::ArrayXi springCoordinateIndexA;
  Eigen::ArrayXi springCoordinateIndexB;
  Eigen::VectorXd springBoxOffset;
  Eigen::VectorXd springContourLength; ///< alpha_i * N_s
  Eigen::VectorXd springPartitions;    ///< alpha_i in (0, 1]
  Eigen::ArrayXi strandIndexOfSpring;
  Eigen::ArrayXb springIsEntanglement;

  // ===========================================================================
  // Lightweight cluster / Universe metadata (no stored Universe copy)
  // ===========================================================================
  size_t universeNrOfAtoms = 0;
  double defaultBondLength = 0.0;
  double defaultR02 = 1.0;
  std::vector<size_t> clusterSizes;
  std::unordered_map<long int, size_t> atomIdToClusterIdx;

  // ===========================================================================
  // Unified Network Builders
  // ===========================================================================

  /**
   * @brief Build a canonical MEHPNetwork from an atomistic Universe using DFS
   * strand contraction.
   *
   * @param universe Atomistic polymer Universe
   * @param entanglements Optional pre-sampled atom-pair entanglements
   * @param crossLinkerType Atom type identifying crosslinker beads
   * @param entanglementMode Entanglement representation mode
   * @param keep2FunctionalCrosslinkers Keep degree-2 atoms of type
   *        @p crossLinkerType as explicit network nodes
   * @param removePrimaryLoops Remove polymer primary loops upon construction
   * @param assumeBoxLargeEnough Whether springs can be evaluated using minimum
   *        image PBC without baked offsets
   * @return MEHPNetwork
   */
  [[nodiscard]] static MEHPNetwork fromUniverse(
    const pylimer_tools::entities::Universe& universe,
    const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements&
      entanglements = {},
    int crossLinkerType = 2,
    EntanglementMode entanglementMode = EntanglementMode::MERGED_LINKS,
    bool keep2FunctionalCrosslinkers = false,
    bool removePrimaryLoops = false,
    bool assumeBoxLargeEnough = false);

  /**
   * @brief Sample random entanglements from an atomistic Universe via
   * randomlyFindEntanglementsV2 and build a canonical MEHPNetwork.
   *
   * @param universe Atomistic polymer Universe
   * @param nrOfEntanglementsToSample Target number of entanglements to sample
   * @param upperSamplingCutoff Upper distance cutoff for pair sampling
   * @param lowerSamplingCutoff Lower distance cutoff for pair sampling
   * @param minNrOfEntanglements Minimum entanglements required if possible
   * @param sameStrandCutoff Minimum chemical distance for same-strand pairs
   * @param seed PRNG seed string (empty for random seed)
   * @param crossLinkerType Atom type identifying crosslinker beads
   * @param entanglementMode Entanglement representation mode
   * @param withCrosslinkers Whether the system contains crosslinkers of type
   *        @p crossLinkerType to exclude as entanglement targets (passed as
   *        ignoreCrosslinks)
   * @param filterDanglingAndSoluble Omit entanglements touching
   * dangling/soluble strands during sampling
   * @param keep2FunctionalCrosslinkers Keep degree-2 crosslinkers as nodes
   * @param removePrimaryLoops Remove polymer primary loops upon construction
   * @param assumeBoxLargeEnough Assume box large enough for minimum-image PBC
   * @return MEHPNetwork
   */
  [[nodiscard]] static MEHPNetwork fromUniverseWithRandomEntanglements(
    const pylimer_tools::entities::Universe& universe,
    size_t nrOfEntanglementsToSample,
    double upperSamplingCutoff,
    double lowerSamplingCutoff = 0.0,
    size_t minNrOfEntanglements = 0,
    double sameStrandCutoff = 0.0,
    const std::string& seed = "",
    int crossLinkerType = 2,
    EntanglementMode entanglementMode = EntanglementMode::MERGED_LINKS,
    bool withCrosslinkers = true,
    bool filterDanglingAndSoluble = false,
    bool keep2FunctionalCrosslinkers = false,
    bool removePrimaryLoops = false,
    bool assumeBoxLargeEnough = false);

  /**
   * @brief Build a canonical MEHPNetwork from a coarse-grained
   * CrosslinkerUniverse and optional StrandPairEntanglements.
   *
   * Preserves crosslink node indices `0..cu.xlinkTypes.size()-1` in
   * `0..nrOfLinks-1` with matching `oldAtomIds(i) == i`.
   *
   * @param box Simulation box
   * @param cu Crosslinker & strand representation
   * @param entanglements Optional pre-sampled strand-pair entanglements
   * @param crossLinkerType Atom type identifying crosslinkers
   * @param entanglementMode Entanglement representation mode
   * @param removeDanglingChains Omit dangling chains when true
   * @param includeFreeChains Include free/soluble chains when true
   * @param removePrimaryLoops Remove polymer primary loops upon construction
   * @param assumeBoxLargeEnough Assume box large enough for minimum-image PBC
   * @return MEHPNetwork
   */
  [[nodiscard]] static MEHPNetwork fromCrosslinkerUniverse(
    const pylimer_tools::entities::Box& box,
    const CrosslinkerUniverse& cu,
    const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
      entanglements = {},
    int crossLinkerType = 2,
    EntanglementMode entanglementMode = EntanglementMode::MERGED_LINKS,
    bool removeDanglingChains = false,
    bool includeFreeChains = true,
    bool removePrimaryLoops = false,
    bool assumeBoxLargeEnough = false);

  // ===========================================================================
  // Metadata Extraction & Averages
  // ===========================================================================

  /**
   * @brief Populate Universe atom count, bond-length statistics, and covalent
   * cluster membership metadata from an atomistic Universe.
   *
   * @param u Atomistic Universe
   */
  void extractUniverseMetadata(const pylimer_tools::entities::Universe& u);

  /**
   * @brief Infer Universe atom count and covalent cluster sizes solely from the
   * current network topology.
   */
  void extractSelfMetadata();

  /**
   * @brief Infer the number of atomistic beads represented by the network.
   *
   * @return double Estimated number of atoms
   */
  [[nodiscard]] double inferNrOfAtomsFromNetwork() const;

  /**
   * @brief Recompute `meanStrandContourLength` and `meanSpringContourLength`
   * from `strandContourLength` and `springContourLength`.
   */
  void recomputeMeanContourLengths();

  // ===========================================================================
  // Slip-Link Insertion & Kinematics
  // ===========================================================================

  /**
   * @brief Compute the clamped inverse spring contour length `(1 / N_part)` for
   * sliding slip-links and fixed sub-springs.
   *
   * @param isPartialSpring True if the spring is a partitioned segment on a
   *        multi-spring strand
   * @param springContourLen Contour length of this spring (`alpha_i * N_s`)
   * @param strandContourLen Total contour length of the parent strand (`N_s`)
   * @param oneOverSpringPartitionUpperLimit Maximum allowed `1 / N_part`
   *        (non-positive disables clamping)
   * @return double Clamped inverse spring contour length
   */
  [[nodiscard]] static double computeClampedInverseSpringContourLength(
    bool isPartialSpring,
    double springContourLen,
    double strandContourLen,
    double oneOverSpringPartitionUpperLimit = 1.0);

  /**
   * @brief Insert merged slip-links onto polymer strands from
   * StrandPairEntanglements.
   *
   * @param entanglements Entanglements specifying strand pairs, coordinates,
   *        and fractional positions
   * @param u Optional displacement vector resized alongside new links
   * @param clampAlpha Clamp extreme fractional positions cleanly away from 0/1
   * @return size_t Number of slip-links inserted
   */
  size_t addSlipLinks(
    const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
      entanglements,
    Eigen::VectorXd* u = nullptr,
    bool clampAlpha = false);

  /**
   * @brief Insert merged slip-links at explicit fractional contour positions.
   *
   * @param strandIdx1 First strand index per entanglement
   * @param strandIdx2 Second strand index per entanglement
   * @param x Initial x coordinates
   * @param y Initial y coordinates
   * @param z Initial z coordinates
   * @param alpha1 Fractional contour position along strand 1
   * @param alpha2 Fractional contour position along strand 2
   * @param u Optional displacement vector resized alongside new links
   * @param clampAlpha Clamp fractional positions cleanly away from boundaries
   * @return size_t Number of slip-links inserted
   */
  size_t addSlipLinks(const std::vector<size_t>& strandIdx1,
                      const std::vector<size_t>& strandIdx2,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& y,
                      const Eigen::VectorXd& z,
                      const std::vector<double>& alpha1,
                      const std::vector<double>& alpha2,
                      Eigen::VectorXd* u = nullptr,
                      bool clampAlpha = false);

  /**
   * @brief Insert merged slip-links at strand midpoints (alpha = 0.5).
   *
   * @param strandIdx1 First strand index per entanglement
   * @param strandIdx2 Second strand index per entanglement
   * @param x Initial x coordinates
   * @param y Initial y coordinates
   * @param z Initial z coordinates
   * @param u Optional displacement vector resized alongside new links
   * @param clampAlpha Clamp fractional positions cleanly away from boundaries
   * @return size_t Number of slip-links inserted
   */
  size_t addSlipLinks(const std::vector<size_t>& strandIdx1,
                      const std::vector<size_t>& strandIdx2,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& y,
                      const Eigen::VectorXd& z,
                      Eigen::VectorXd* u = nullptr,
                      bool clampAlpha = false);

  /**
   * @brief Sample random entanglements directly on the current network geometry
   * using a temporary CrosslinkerUniverse and insert them as slip-links.
   *
   * @param nrOfEntanglementsToSample Target number of slip-links to sample
   * @param cutoff Upper distance cutoff
   * @param minNrOfEntanglements Minimum entanglements required if possible
   * @param sameStrandCutoff Minimum bead separation along the same strand
   * @param crossLinkerType Atom type treated as fixed crosslinker endpoint
   * @param seed PRNG seed string
   * @param u Optional current displacements updated alongside added slip-links
   * @param clampAlpha Clamp sampled fractional positions cleanly away from 0/1
   * @return size_t Number of slip-links inserted
   */
  size_t randomlyAddSlipLinks(size_t nrOfEntanglementsToSample,
                              double cutoff = 2.0,
                              size_t minNrOfEntanglements = 0,
                              double sameStrandCutoff = 2.0,
                              int crossLinkerType = 2,
                              const std::string& seed = "",
                              Eigen::VectorXd* u = nullptr,
                              bool clampAlpha = false);

  /**
   * @brief Realign box offsets of two consecutive sub-springs meeting at an
   * interior slip-link to minimize squared spring lengths while preserving the
   * combined end-to-end vector.
   *
   * @param u Current displacements
   * @param slipLinkIdx Interior slip-link index
   * @param springIdx1 Incoming sub-spring ending at @p slipLinkIdx
   * @param springIdx2 Outgoing sub-spring starting at @p slipLinkIdx
   * @param is2D Treat system as 2D when computing squared distances
   */
  void realignSlipLinkToImages(const Eigen::VectorXd& u,
                               size_t slipLinkIdx,
                               size_t springIdx1,
                               size_t springIdx2,
                               bool is2D = false);

  /**
   * @brief Update `springPartitions` and synced `springContourLength` around an
   * interior slip-link to equalize tension along its strand(s).
   *
   * @param linkIdx Entanglement link index (`linkIsEntanglement[linkIdx] ==
   * true`)
   * @param u Current displacements
   * @param oneOverSpringPartitionUpperLimit Upper limit on `1 / N_part`
   * @param is2D Treat system as 2D
   * @return double Mean local contour-partition residual norm
   */
  double updateSpringPartitionAtLink(
    size_t linkIdx,
    const Eigen::VectorXd& u,
    double oneOverSpringPartitionUpperLimit = 1.0,
    bool is2D = false);

  /**
   * @brief Displace all links one Jacobi-like step toward the weighted mean
   * equilibrium position of connected springs.
   *
   * @param u Displacement vector updated in-place
   * @param is2D Treat system as 2D
   * @param oneOverSpringPartitionUpperLimit Upper limit for inverse partition
   * @return double Maximum 3D displacement norm among all links
   */
  double displaceToMeanPosition(
    Eigen::VectorXd& u,
    bool is2D = false,
    double oneOverSpringPartitionUpperLimit = 1.0) const;

  /**
   * @brief Displace a single link to the weighted equilibrium position of its
   * connected springs.
   *
   * @param u Displacement vector updated in-place
   * @param linkIdx Target link index
   * @param is2D Treat system as 2D
   * @param oneOverSpringPartitionUpperLimit Upper limit for inverse partition
   * @return double Norm of displacement step applied to @p linkIdx
   */
  double displaceToMeanPosition(
    Eigen::VectorXd& u,
    size_t linkIdx,
    bool is2D = false,
    double oneOverSpringPartitionUpperLimit = 1.0) const;

  /**
   * @brief Merge the two springs connected to a bifunctional link into one.
   *
   * Supports both bifunctional crosslinks (merging two strands into one while
   * keeping `springPartitions` normalized) and interior slip-links along a
   * single strand (coalescing two adjacent sub-springs). Does not assume
   * crosslinks or entanglements reside in contiguous index ranges.
   *
   * @param linkIdx Link index whose two springs are merged
   * @param u Displacement vector resized if @p removeIfZeroFunctional is true
   * @param removeIfZeroFunctional Also erase @p linkIdx if it becomes
   * 0-functional
   * @return true if two springs were merged; false if @p linkIdx is not
   *         eligible (e.g., functionality != 2, looping, or entanglement
   * spring)
   */
  bool mergeSpringsAtLink(size_t linkIdx,
                          Eigen::VectorXd& u,
                          bool removeIfZeroFunctional = true);

  /**
   * @brief Remove a single 0-functional link and renumber downstream indices.
   *
   * @param linkIdx Link index to remove (must have 0 attached strands)
   * @param u Optional displacement vector resized alongside @p coordinates
   */
  void removeLink(size_t linkIdx, Eigen::VectorXd* u = nullptr);

  /**
   * @brief Remove a single 0-functional link and resize displacement vector @p
   * u.
   *
   * @param linkIdx Link index to remove (must have 0 attached strands)
   * @param u Displacement vector resized alongside @p coordinates
   */
  void removeLink(size_t linkIdx, Eigen::VectorXd& u);

  /**
   * @brief Remove a single spring (and its parent strand if non-empty, along
   * with attached entanglement coupling springs) from the network.
   *
   * @param springIdx Spring index to remove
   * @param removeZeroFunctionalLinks Whether to prune newly 0-functional links
   * @param u Optional displacement vector resized if links are removed
   */
  void removeSpring(size_t springIdx,
                    bool removeZeroFunctionalLinks = true,
                    Eigen::VectorXd* u = nullptr);

  /**
   * @brief Swap the ordering of two adjacent interior slip-links sharing
   * sub-spring @p springIdx along a polymer strand.
   *
   * @param springIdx Interior sub-spring connecting two slip-links
   * @param u Current displacements (used for invariant verification)
   * @param is2D Treat system as 2D
   */
  void swapSlipLinksOnSpring(size_t springIdx,
                             const Eigen::VectorXd& u,
                             bool is2D = false);

  /**
   * @brief Move an interior slip-link attached to a crosslink along
   * @p sourceSpringIdx onto another sub-spring @p targetSpringIdx attached to
   * the same crosslink.
   *
   * @param sourceSpringIdx Sub-spring between crosslink and slip-link
   * @param targetSpringIdx Target sub-spring touching the same crosslink
   * @param u Current displacements
   * @param oneOverSpringPartitionUpperLimit Upper limit on `1 / N_part`
   * @param is2D Treat system as 2D
   * @return long int New spring index connecting the crosslink and moved
   *         slip-link, or -1 if insufficient contour space on target strand
   */
  long int moveSlipLinkFromRailToRail(
    size_t sourceSpringIdx,
    size_t targetSpringIdx,
    Eigen::VectorXd& u,
    double oneOverSpringPartitionUpperLimit = 1.0,
    bool is2D = false);

  /**
   * @brief Rotate an interior slip-link connected via @p springIdx around its
   * adjacent crosslink onto another attached polymer branch.
   *
   * @param springIdx Terminal sub-spring connecting slip-link and crosslink
   * @param u Current displacements
   * @param oneOverSpringPartitionUpperLimit Upper limit on `1 / N_part`
   * @param respectLoops Avoid moves that separate a primary loop strand
   * @param is2D Treat system as 2D
   * @return long int New spring index on success, or -1 if no valid move
   */
  long int rotateSlipLinkAroundCrosslink(
    size_t springIdx,
    Eigen::VectorXd& u,
    double oneOverSpringPartitionUpperLimit = 1.0,
    bool respectLoops = false,
    bool is2D = false);

  /**
   * @brief Attempt to move or swap a slip-link adjacent to sub-spring
   * @p springIdx according to @p mode.
   *
   * @param springIdx Candidate sub-spring index
   * @param u Current displacements (updated if MC acceptance relaxes neighbors)
   * @param mode Swapping policy
   * @param oneOverSpringPartitionUpperLimit Upper limit on `1 / N_part`
   * @param nrOfCrosslinkSwapsAllowedPerSliplink Maximum crosslink hops per
   *        slip-link (`< 0` for unlimited)
   * @param is2D Treat system as 2D
   * @return true if a swap or branch hop was accepted and performed
   */
  bool moveSlipLinkAlongStrand(size_t springIdx,
                               Eigen::VectorXd& u,
                               LinkSwappingMode mode,
                               double oneOverSpringPartitionUpperLimit = 1.0,
                               int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
                               bool is2D = false);

  /**
   * @brief Execute a full pass of slip-link swaps and/or crosslink hops across
   * the network according to @p mode.
   *
   * @param u Current displacements
   * @param mode Swapping policy
   * @param oneOverSpringPartitionUpperLimit Upper limit on `1 / N_part`
   * @param nrOfCrosslinkSwapsAllowedPerSliplink Maximum crosslink hops per
   *        slip-link (`< 0` for unlimited)
   * @param is2D Treat system as 2D
   * @return size_t Number of swaps/hops performed
   */
  size_t performSlipLinkSwaps(Eigen::VectorXd& u,
                              LinkSwappingMode mode,
                              double oneOverSpringPartitionUpperLimit = 1.0,
                              int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
                              bool is2D = false);

  // ===========================================================================
  // Simplification, Pruning, Deformation & Matrix Update Support
  // ===========================================================================

  /**
   * @brief Find spring indices belonging to strands exceeding
   * @p breakingLength (along with attached entanglement coupling springs).
   *
   * @param u Current displacements
   * @param breakingLength Maximum allowed weighted sub-spring length `|r|/N_i`
   * @param is2D Treat system as 2D
   * @param outNumBrokenStrands Optional output count of broken strands
   * @return std::vector<size_t> Unique spring indices sorted descending
   */
  [[nodiscard]] std::vector<size_t> findTooLongStrandSprings(
    const Eigen::VectorXd& u,
    double breakingLength,
    bool is2D = false,
    size_t* outNumBrokenStrands = nullptr) const;

  /**
   * @brief Break and remove all strands containing a spring whose weighted
   * length `|r_i| / N_i` exceeds @p breakingLength.
   *
   * Link indices are preserved (`removeZeroFunctionalLinks = false`) so that
   * the caller's `const Eigen::VectorXd& u` remains conformant.
   *
   * @param u Current displacements
   * @param breakingLength Breaking threshold on `|r_i| / N_i`
   * @param is2D Treat system as 2D
   * @return size_t Number of strands broken
   */
  size_t breakTooLongStrands(const Eigen::VectorXd& u,
                             double breakingLength,
                             bool is2D = false);

  /**
   * @brief Mark inactive springs to delete (including attached entanglement
   * coupling springs when a polymer strand is inactive).
   *
   * @param u Current displacements
   * @param tolerance Activity criterion tolerance
   * @param is2D Treat system as 2D
   * @param checkAll Also mark orphaned/dangling springs directly
   * @return Eigen::ArrayXb Boolean mask sized `nrOfSprings`
   */
  [[nodiscard]] Eigen::ArrayXb findInactiveSprings(const Eigen::VectorXd& u,
                                                   double tolerance = 1e-9,
                                                   bool is2D = false,
                                                   bool checkAll = false) const;

  /**
   * @brief Remove specified springs without renumbering link indices and
   * collect `BrokenSpringInfo` records for in-place CSR matrix updates.
   *
   * If @p constants is non-null (sized `3 * nrOfLinks`), its RHS force-balance
   * entries are also updated in-place to subtract removed box-offset terms.
   *
   * @param springIndicesToRemove Spring indices to remove
   * @param kappa Spring stiffness scale factor
   * @param constants Optional RHS vector updated in-place
   * @param is2D Treat system as 2D
   * @return std::vector<BrokenSpringInfo> Rank-1 matrix update descriptors
   */
  std::vector<BrokenSpringInfo> removeSpringsAndCollectMatrixUpdates(
    const std::vector<size_t>& springIndicesToRemove,
    double kappa = 1.0,
    Eigen::VectorXd* constants = nullptr,
    bool is2D = false);

  /**
   * @brief Remove springs (and any empty/deleted strands), optionally pruning
   * newly zero-functional links and resizing @p u.
   *
   * Accepts indices in arbitrary order and handles duplicates cleanly.
   *
   * @param springIndicesToRemove Spring indices to remove
   * @param removeZeroFunctionalLinks Whether to remove zero-functional links
   * @param u Optional displacement vector resized alongside @p coordinates
   */
  void removeSprings(const std::vector<size_t>& springIndicesToRemove,
                     bool removeZeroFunctionalLinks = true,
                     Eigen::VectorXd* u = nullptr);

  /**
   * @brief Remove inactive springs, unlink bifunctional links, and remove
   * zero-functional links.
   *
   * @param u Current displacements resized alongside network links
   * @param tolerance Activity tolerance
   * @param is2D Treat system as 2D
   * @return size_t Number of zero-functional links removed
   */
  size_t removeInactiveLinks(Eigen::VectorXd& u,
                             double tolerance = 1e-9,
                             bool is2D = false);

  /**
   * @brief Merge consecutive springs passing through bifunctional links.
   *
   * @param u Displacement vector resized when zero-functional links are pruned
   * @param entanglementsOnly If true, only unlink bifunctional entanglement
   *        links (`linkIsEntanglement[l] == true`)
   * @return size_t Number of bifunctional links unlinked
   */
  size_t unlinkBifunctionalLinks(Eigen::VectorXd& u,
                                 bool entanglementsOnly = false);

  /**
   * @brief Remove all links with no connected strands
   * (`strandIndicesOfLink[l].empty()`).
   *
   * @param u Displacement vector resized alongside @p coordinates
   * @return size_t Number of links removed
   */
  size_t removeZeroFunctionalLinks(Eigen::VectorXd& u);

  /**
   * @brief Overload accepting an optional displacement pointer.
   *
   * @param u Optional displacement vector pointer
   * @return size_t Number of links removed
   */
  size_t removeZeroFunctionalLinks(Eigen::VectorXd* u);

  /**
   * @brief Run iterative topology simplification (breaking, inactive removal,
   * and/or bifunctional unlinking) until convergence.
   *
   * @param u Current displacements updated in-place
   * @param mode Simplification strategy
   * @param tolerance Activity tolerance
   * @param is2D Treat system as 2D
   * @param breakingLength Optional spring breaking length (> 0 enables
   * breaking)
   * @return size_t Total strands broken plus springs and links removed
   */
  size_t simplify(Eigen::VectorXd& u,
                  StructureSimplificationMode mode,
                  double tolerance = 1e-9,
                  bool is2D = false,
                  double breakingLength = -1.0);

  /**
   * @brief Affinely deform coordinates, displacements @p u, and baked box
   * offsets to @p targetBox, updating `L`, `boxHalfs`, `vol`, and `box`.
   *
   * @param targetBox Target simulation box
   * @param u Current displacement vector updated in-place (ignored if empty)
   */
  void deformTo(const pylimer_tools::entities::Box& targetBox,
                Eigen::VectorXd& u);

  /**
   * @brief Affinely deform coordinates and baked box offsets to @p targetBox.
   *
   * @param targetBox Target simulation box
   */
  void deformTo(const pylimer_tools::entities::Box& targetBox);

  // ===========================================================================
  // Vector, Length & Stress Evaluation
  // ===========================================================================

  /**
   * @brief Evaluate all spring end-to-end vectors `(x_B + u_B) - (x_A + u_A) +
   * boxOffset`. Applies `box.handlePBC` when `assumeBoxLargeEnough` is true.
   *
   * @param u Current displacements (`3 * nrOfLinks`)
   * @param is2D Zero z components when true
   * @return Eigen::VectorXd Flat vector sized `3 * nrOfSprings`
   */
  [[nodiscard]] Eigen::VectorXd evaluateSpringVectors(const Eigen::VectorXd& u,
                                                      bool is2D = false) const;

  /**
   * @brief Evaluate the 3D vector for a single spring. Applies `box.handlePBC`
   * when `assumeBoxLargeEnough` is true.
   *
   * @param u Current displacements (`3 * nrOfLinks`)
   * @param springIdx Spring index
   * @param is2D Zero z component when true
   * @return Eigen::Vector3d Spring vector A -> B
   */
  [[nodiscard]] Eigen::Vector3d evaluateSpringVector(const Eigen::VectorXd& u,
                                                     size_t springIdx,
                                                     bool is2D = false) const;

  /**
   * @brief Evaluate the Euclidean lengths of all springs.
   *
   * @param u Current displacements (`3 * nrOfLinks`)
   * @param is2D Treat system as 2D
   * @return Eigen::VectorXd Spring lengths sized `nrOfSprings`
   */
  [[nodiscard]] Eigen::VectorXd evaluateSpringLengths(const Eigen::VectorXd& u,
                                                      bool is2D = false) const;

  /**
   * @brief Evaluate the vector of @p springIdx oriented starting at @p linkIdx.
   */
  [[nodiscard]] Eigen::Vector3d evaluateSpringVectorFrom(
    const Eigen::VectorXd& u,
    size_t springIdx,
    size_t linkIdx,
    bool is2D = false) const;

  /**
   * @brief Evaluate the vector of @p springIdx oriented ending at @p linkIdx.
   */
  [[nodiscard]] Eigen::Vector3d evaluateSpringVectorTo(const Eigen::VectorXd& u,
                                                       size_t springIdx,
                                                       size_t linkIdx,
                                                       bool is2D = false) const;

  /**
   * @brief Evaluate the contour-summed vector of each strand.
   */
  [[nodiscard]] Eigen::VectorXd evaluateStrandVectors(const Eigen::VectorXd& u,
                                                      bool is2D = false) const;

  /**
   * @brief Evaluate the end-to-end vector of a single strand.
   */
  [[nodiscard]] Eigen::Vector3d evaluateStrandVector(const Eigen::VectorXd& u,
                                                     size_t strandIdx,
                                                     bool is2D = false) const;

  /**
   * @brief Evaluate the summed sub-spring contour-path lengths per strand.
   */
  [[nodiscard]] Eigen::VectorXd evaluateStrandLengths(const Eigen::VectorXd& u,
                                                      bool is2D = false) const;

  /**
   * @brief Evaluate weighted spring length `|r_i| / N_i`.
   */
  [[nodiscard]] double evaluateWeightedSpringLength(const Eigen::VectorXd& u,
                                                    size_t springIdx,
                                                    bool is2D = false) const;

  /**
   * @brief Evaluate the 3x3 virial/spring stress tensor summed over all
   * springs.
   *
   * @param u Current displacements
   * @param kappa Spring constant prefactor
   * @param is2D Treat system as 2D
   * @param oneOverSpringPartitionUpperLimit Clamp inverse partial contour
   *        lengths when `> 0`; pass `<= 0` for raw `1 / springContourLength`
   * @return Eigen::Matrix3d Stress tensor
   */
  [[nodiscard]] Eigen::Matrix3d evaluateStressTensor(
    const Eigen::VectorXd& u,
    double kappa = 1.0,
    bool is2D = false,
    double oneOverSpringPartitionUpperLimit = -1.0) const;

  /**
   * @brief Evaluate the link-based virial stress tensor restricted to a set of
   * links (used for local reversible slip-link swaps).
   */
  [[nodiscard]] Eigen::Matrix3d evaluateStressTensorForLinks(
    const std::vector<size_t>& linkIndices,
    const Eigen::VectorXd& u,
    double kappa = 1.0,
    bool is2D = false,
    double oneOverSpringPartitionUpperLimit = 1.0) const;

  // ===========================================================================
  // Topology Queries
  // ===========================================================================

  [[nodiscard]] std::vector<size_t> getSpringIndicesOfLink(
    size_t linkIdx) const;
  [[nodiscard]] std::vector<size_t> getNeighbourLinkIndices(
    size_t linkIdx) const;
  [[nodiscard]] std::vector<size_t> getEntanglementSpringsAtStrand(
    size_t strandIdx) const;
  [[nodiscard]] size_t getOtherLinkOfSpring(size_t springIdx,
                                            size_t linkIdx) const;
  [[nodiscard]] bool isPartOfSpring(size_t linkIdx, size_t springIdx) const;
  [[nodiscard]] bool isLoopingSpring(size_t springIdx) const;
  [[nodiscard]] bool isLoopingStrand(size_t strandIdx) const;
  [[nodiscard]] bool isEntanglementStrand(size_t strandIdx) const;

  [[nodiscard]] static bool distanceIsWithinTolerance(
    const Eigen::Vector3d& dist,
    double tolerance = 1e-3,
    double contourLength = 1.0,
    double contourLengthFraction = 1.0);

  // ===========================================================================
  // Network Observables & Topological Analysis
  // ===========================================================================

  [[nodiscard]] Eigen::ArrayXb findActiveSprings(const Eigen::VectorXd& u,
                                                 double tolerance = 1e-3,
                                                 bool is2D = false) const;
  [[nodiscard]] Eigen::ArrayXb findActiveStrands(const Eigen::VectorXd& u,
                                                 double tolerance = 1e-3,
                                                 bool is2D = false) const;
  [[nodiscard]] Eigen::ArrayXb findActiveStrandsInDir(const Eigen::VectorXd& u,
                                                      int dir,
                                                      double tolerance = 1e-3,
                                                      bool is2D = false) const;
  [[nodiscard]] Eigen::ArrayXb findActiveLinks(const Eigen::VectorXd& u,
                                               double tolerance = 1e-3,
                                               bool is2D = false) const;
  [[nodiscard]] Eigen::ArrayXb findActiveNodes(const Eigen::VectorXd& u,
                                               double tolerance = 1e-3,
                                               bool is2D = false) const;
  [[nodiscard]] std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
  findSpringsAndLinksClusteredToActive(const Eigen::VectorXd& u,
                                       double tolerance = 1e-3,
                                       bool is2D = false) const;
  [[nodiscard]] std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
  findStrandsAndNodesClusteredToActive(const Eigen::VectorXd& u,
                                       double tolerance = 1e-3,
                                       bool is2D = false) const;
  [[nodiscard]] std::vector<int> getIndicesOfActiveNodes(
    const Eigen::VectorXd& u,
    double tolerance = 1e-3,
    bool is2D = false) const;

  [[nodiscard]] double computeActiveWeightFraction(const Eigen::VectorXd& u,
                                                   double tolerance = 1e-3,
                                                   bool is2D = false) const;
  [[nodiscard]] double computeSolubleWeightFraction(const Eigen::VectorXd& u,
                                                    double tolerance = 1e-3,
                                                    bool is2D = false) const;
  [[nodiscard]] double computeDanglingWeightFraction(const Eigen::VectorXd& u,
                                                     double tolerance = 1e-3,
                                                     bool is2D = false) const;
  [[nodiscard]] double countActiveClusteredAtoms(const Eigen::VectorXd& u,
                                                 double tolerance = 1e-3,
                                                 bool is2D = false) const;

  [[nodiscard]] std::unordered_map<long int, int>
  getEffectiveFunctionalityOfAtoms(const Eigen::VectorXd& u,
                                   double tolerance = 1e-3,
                                   bool is2D = false) const;

  [[nodiscard]] Eigen::VectorXd getGammaFactors(const Eigen::VectorXd& u,
                                                double b02 = -1.0,
                                                bool is2D = false) const;
  [[nodiscard]] Eigen::VectorXd getGammaFactorsInDir(const Eigen::VectorXd& u,
                                                     double b02,
                                                     int dir,
                                                     bool is2D = false) const;
  [[nodiscard]] double getGammaFactor(const Eigen::VectorXd& u,
                                      double b02 = -1.0,
                                      int nrOfChains = -1,
                                      bool is2D = false) const;

  /**
   * @brief Reconstruct a coarse-grained Universe of crosslinker nodes bonded by
   * polymer strands.
   *
   * Uses `linkIsEntanglement` strictly (never index range assumptions) and
   * ignores entanglement coupling strands.
   *
   * @param u Optional current displacements added to `coordinates`
   * @return pylimer_tools::entities::Universe
   */
  [[nodiscard]] pylimer_tools::entities::Universe getCrosslinkerVerse(
    const Eigen::VectorXd& u = Eigen::VectorXd()) const;

  /**
   * @brief Count interior entanglement links that appear more than once on a
   * single strand (intra-chain slip-links).
   *
   * @return size_t Number of intra-chain slip-link occurrences
   */
  [[nodiscard]] size_t getNumIntraChainSlipLinks() const;

  // ===========================================================================
  // Validation & Serialization
  // ===========================================================================

  /**
   * @brief Validate all network dimensions, index mappings, contour/partition
   * invariants, box dimensions, and finiteness.
   *
   * @param u Optional displacement vector (`3 * nrOfLinks` if non-empty)
   * @throws std::runtime_error on any structural or numerical inconsistency
   */
  void validate(const Eigen::VectorXd& u = Eigen::VectorXd()) const;

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(L,
       boxHalfs,
       vol,
       box,
       nrOfLinks,
       nrOfNodes,
       nrOfStrands,
       nrOfSprings,
       meanStrandContourLength,
       meanSpringContourLength,
       assumeBoxLargeEnough,
       removePrimaryLoops,
       assumeComplete,
       coordinates,
       linkIsEntanglement,
       oldAtomIds,
       oldAtomTypes,
       strandIndicesOfLink,
       nrOfCrosslinkSwapsEndured,
       strandContourLength,
       linkIndicesOfStrand,
       springIndicesOfStrand,
       springIndexA,
       springIndexB,
       springCoordinateIndexA,
       springCoordinateIndexB,
       springBoxOffset,
       springContourLength,
       springPartitions,
       strandIndexOfSpring,
       springIsEntanglement,
       universeNrOfAtoms,
       defaultBondLength,
       defaultR02,
       clusterSizes,
       atomIdToClusterIdx);
  }
#endif
};

} // namespace pylimer_tools::sim::mehp

#endif // MEHP_NETWORK_H
