#pragma once

#include "../entities/Box.h"
#include "../topo/EntanglementDetector.h"
#include "MEHPUtilityStructures.h"

namespace pylimer_tools::sim::mehp {

/**
 * @brief Build a Network (for MEHPForceRelaxation) from a CrosslinkerUniverse.
 *
 * @param box Simulation box
 * @param cu Canonical CrosslinkerUniverse representation
 * @param includeDanglingChains Whether to include dangling strands as springs
 * @param includeFreeChains Whether to include free/soluble strands as springs
 * @return Network
 */
Network
buildForceRelaxationNetwork(const pylimer_tools::entities::Box& box,
                            const CrosslinkerUniverse& cu,
                            bool includeDanglingChains = true,
                            bool includeFreeChains = true);

/**
 * @brief Build a ForceBalanceNetwork (for MEHPForceBalance) from a
 * CrosslinkerUniverse.
 *
 * @param box Simulation box
 * @param cu Canonical CrosslinkerUniverse representation
 * @param includeDanglingChains Whether to include dangling strands as springs
 * @param includeFreeChains Whether to include free/soluble strands as springs
 * @return ForceBalanceNetwork
 */
ForceBalanceNetwork
buildForceBalanceNetwork(const pylimer_tools::entities::Box& box,
                         const CrosslinkerUniverse& cu,
                         bool includeDanglingChains = true,
                         bool includeFreeChains = true);

/**
 * @brief Build a ForceBalance2Network (for MEHPForceBalance2) from a
 * CrosslinkerUniverse and optional StrandPairEntanglements.
 *
 * @param box Simulation box
 * @param cu Canonical CrosslinkerUniverse representation
 * @param entanglements Entanglements along strands of @p cu
 * @param crossLinkerType Atom type of crosslinkers (unused when cu already
 *                        carries xlinkTypes, kept for API symmetry)
 * @param entanglementsAsSprings Whether to model each entanglement as a spring
 *                               connecting two slip-links rather than a merged
 *                               node
 * @param includeDanglingChains Whether to include dangling strands
 * @param includeFreeChains Whether to include free/soluble strands
 * @return ForceBalance2Network
 */
ForceBalance2Network
buildForceBalance2Network(
  const pylimer_tools::entities::Box& box,
  const CrosslinkerUniverse& cu,
  const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
    entanglements = {},
  int crossLinkerType = 2,
  bool entanglementsAsSprings = false,
  bool includeDanglingChains = true,
  bool includeFreeChains = true);

} // namespace pylimer_tools::sim::mehp
