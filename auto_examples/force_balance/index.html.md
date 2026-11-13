<a id="sphx-glr-auto-examples-force-balance"></a>

# Force Balance, Maximum Entropy Homogenization Procedures

pylimer-tools provides three implementations of the Force Balance [[BG25](../../acknowledgements.html.md#id12)], Maximum Entropy Homogenization Procedure (MEHP) [[Gus19](../../acknowledgements.html.md#id8)] to reduce polymer networks to their minimum energy, maximum entropy homogenized state.
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- [`MEHPForceRelaxation`](../../api/pylimer_tools_cpp.MEHPForceRelaxation.html.md#pylimer_tools_cpp.MEHPForceRelaxation): An MEHP implementation that allows the use of non-linear force potentials,
- [`MEHPForceBalance`](../../api/pylimer_tools_cpp.MEHPForceBalance.html.md#pylimer_tools_cpp.MEHPForceBalance): A Force Balance implementation using Hookean springs that supports slip-links with optional dynamic contour slipping (or fixed links via `disable_slipping=True`),
- [`MEHPForceBalance2`](../../api/pylimer_tools_cpp.MEHPForceBalance2.html.md#pylimer_tools_cpp.MEHPForceBalance2): The faster implementation of the Force Balance procedure, representing entanglements as fixed 4-functional links (like tetrafunctional crosslinks) or compliant Hookean springs, without contour slipping.

## Sampling Entanglements

Entanglements can be sampled in three ways:

1. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a [`Universe`](../../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe).
2. **On-demand bounding-box sampling**: Computes per-strand bounding boxes and regenerates candidate strand beads on demand, matching discrete bead-level fidelity with ~80–95% lower peak memory.
3. **Analytical Brownian-bridge sampling**: Directly on the crosslinkers and strand connectivity of a [`CrosslinkerUniverse`](../../api/pylimer_tools_cpp.CrosslinkerUniverse.html.md#pylimer_tools_cpp.CrosslinkerUniverse) without instantiating intermediate strand beads, reducing peak memory by ~80–95%.

All three methods support [`EntanglementFilteringMode`](../../api/pylimer_tools_cpp.EntanglementFilteringMode.html.md#pylimer_tools_cpp.EntanglementFilteringMode) for topological leaf peeling of dangling and soluble fractions.

See:

- [Entanglement Sampling Methods](plot_entanglement_sampling_comparison.html.md) for a detailed benchmark and comparison of discrete vs. analytical sampling methods,
- [Entanglement Links vs. Springs](plot_entanglement_link_vs_spring.html.md) for the difference between modeling entanglements as rigid merged links versus compliant Hookean springs in Force Balance 2,
- [Slipping vs. Fixed Links](plot_slipping_vs_fixed_links.html.md) for a comparison of dynamic contour slipping versus fixed entanglement links in Force Balance 1, and
- [Inactive Removal & Sampling](plot_inactive_removal_and_sampling.html.md) for the removal of inactive strands/entanglements in Force Balance 2 and the impact of pre-filtering dangling chains during entanglement sampling,
- [Polydispersity Study](plot_polydispersity_study.html.md) for the effect of precursor strand length polydispersity on equilibrium shear modulus and defect fractions, and
- [Sampling Cut-Off Effect](plot_sampling_cutoff_effect.html.md) for the effect of the slip-link sampling cut-off distance on the equilibrium shear modulus in Force Balance 2, investigating pre-strain distortions, link vs. spring representations, and interactions with strand length $N$, conversion $p$, and stoichiometric imbalance $r$.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example shows how to use the MEHPForceBalance2 class in a deformation experiment.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_fb_deformation_exp_thumb.png)

[Deformation Experiment](plot_fb_deformation_exp.html.md)

  <div class="sphx-glr-thumbnail-title">Deformation Experiment</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the Maximum Entropy Homogenization Procedure (MEHP) (gusev_numerical_2019) with slip-links to homogenize a polymer network and predict its equilibrium shear modulus.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_balance_thumb.png)

[Force Balance](plot_force_balance.html.md)

  <div class="sphx-glr-thumbnail-title">Force Balance</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="MEHPForceBalance2 is a modernized, high-performance implementation of the Maximum Entropy Homogenization Procedure (MEHP) (gusev_numerical_2019, bernhard_phantom_2025).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_balance_2_thumb.png)

[Force Balance 2](plot_force_balance_2.html.md)

  <div class="sphx-glr-thumbnail-title">Force Balance 2</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the &quot;original&quot; MEHP implementation (see gusev_numerical_2019) without slip-links or entanglements.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_relaxation_thumb.png)

[MEHP Force Relaxation](plot_force_relaxation.html.md)

  <div class="sphx-glr-thumbnail-title">MEHP Force Relaxation</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In polymer network elasticity simulations using the Maximum Entropy Homogenization Procedure (MEHP) (gusev_numerical_2019, bernhard_phantom_2025), topological entanglements impose constraints on polymer chain fluctuations.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_entanglement_link_vs_spring_thumb.png)

[Entanglement Links vs. Springs](plot_entanglement_link_vs_spring.html.md)

  <div class="sphx-glr-thumbnail-title">Entanglement Links vs. Springs</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In the first-generation Maximum Entropy Homogenization Procedure (MEHP) Force Balance (~pylimer_tools_cpp.MEHPForceBalance, bernhard_phantom_2025), topological entanglements can be modeled with dynamic mobility along the polymer chain contour.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_slipping_vs_fixed_links_thumb.png)

[Slipping vs. Fixed Links](plot_slipping_vs_fixed_links.html.md)

  <div class="sphx-glr-thumbnail-title">Slipping vs. Fixed Links</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In this example, we study the effect of polymer strand length polydispersity on the equilibrium shear modulus and topological network structure of end-linked polymer networks.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_polydispersity_study_thumb.png)

[Polydispersity Study](plot_polydispersity_study.html.md)

  <div class="sphx-glr-thumbnail-title">Polydispersity Study</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In polymer network simulations using the Maximum Entropy Homogenization Procedure (MEHP) and Force Balance (gusev_numerical_2019, bernhard_phantom_2025), topological entanglements are introduced by sampling candidate pairs of chain beads that lie within a spatial proximity cut-off:">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_sampling_cutoff_effect_thumb.png)

[Sampling Cut-Off Effect](plot_sampling_cutoff_effect.html.md)

  <div class="sphx-glr-thumbnail-title">Sampling Cut-Off Effect</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In polymer network modeling via the Maximum Entropy Homogenization Procedure (MEHP) and Force Balance (gusev_numerical_2019, bernhard_phantom_2025), networks generated below or near full conversion contain significant fractions of dangling chains, sol fractions, and topological loops.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_inactive_removal_and_sampling_thumb.png)

[Inactive Removal & Sampling](plot_inactive_removal_and_sampling.html.md)

  <div class="sphx-glr-thumbnail-title">Inactive Removal & Sampling</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In polymer network elasticity simulations using the Maximum Entropy Homogenization Procedure (MEHP) and Force Balance (gusev_numerical_2019, bernhard_phantom_2025), entanglements are represented as localized topological constraints (entanglement links or springs).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_entanglement_sampling_comparison_thumb.png)

[Entanglement Sampling Methods](plot_entanglement_sampling_comparison.html.md)

  <div class="sphx-glr-thumbnail-title">Entanglement Sampling Methods</div>
</div>
<!-- thumbnail-parent-div-close --></div>
