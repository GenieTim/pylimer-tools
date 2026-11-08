<a id="sphx-glr-auto-examples-force-balance"></a>

# Force Balance, Maximum Entropy Homogenization Procedures

pylimer-tools provides three implementations of the Force Balance:cite:p:bernhard_phantom_2025, Maximum Entropy Homogenization Procedure (MEHP):cite:p:gusev_numerical_2019 to reduce polymer networks to their minimum energy, maximum entropy homogenized state.
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- [`MEHPForceRelaxation`](../../api/pylimer_tools_cpp.MEHPForceRelaxation.html.md#pylimer_tools_cpp.MEHPForceRelaxation): An MEHP implementation that allows the use of non-linear force potentials,
- [`MEHPForceBalance`](../../api/pylimer_tools_cpp.MEHPForceBalance.html.md#pylimer_tools_cpp.MEHPForceBalance): A Force Balance implementation using Hookean springs that supports slip-links with optional dynamic contour slipping (or fixed links via `disable_slipping=True`),
- [`MEHPForceBalance2`](../../api/pylimer_tools_cpp.MEHPForceBalance2.html.md#pylimer_tools_cpp.MEHPForceBalance2): The faster implementation of the Force Balance procedure, representing entanglements as fixed 4-functional links (like tetrafunctional crosslinks) or compliant Hookean springs, without contour slipping.

## Sampling Entanglements

Entanglements can be sampled in two distinct ways:

1. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a [`Universe`](../../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe).
2. **Analytical Brownian-bridge sampling**: Directly on the crosslinkers and strand connectivity of a [`CrosslinkerUniverse`](../../api/pylimer_tools_cpp.CrosslinkerUniverse.html.md#pylimer_tools_cpp.CrosslinkerUniverse) without instantiating intermediate strand beads, reducing peak memory by up to 90%+ and supporting topological leaf peeling of dangling and soluble fractions.

See:

- [Comparing Slip-Link Sampling Methods: Discrete vs. Analytical](plot_entanglement_sampling_comparison.html.md) for a detailed benchmark and comparison of discrete vs. analytical sampling methods,
- plot_entanglement_link_vs_spring for the difference between modeling entanglements as rigid merged links versus compliant Hookean springs in Force Balance 2,
- plot_slipping_vs_fixed_links for a comparison of dynamic contour slipping versus fixed entanglement links in Force Balance 1, and
- plot_inactive_removal_and_sampling for the removal of inactive strands/entanglements in Force Balance 2 and the impact of pre-filtering dangling chains during entanglement sampling.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>
