<a id="sphx-glr-auto-examples-force-balance"></a>

# Force Balance, Maximum Entropy Homogenization Procedures

pylimer-tools provides three implementations of the Force Balance:cite:p:bernhard_phantom_2025, Maximum Entropy Homogenization Procedure (MEHP):cite:p:gusev_numerical_2019 to reduce polymer networks to their minimum energy, maximum entropy homogenized state.
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- [`MEHPForceRelaxation`](../../api/pylimer_tools_cpp.MEHPForceRelaxation.md#pylimer_tools_cpp.MEHPForceRelaxation): An MEHP implementation that allows to use non-linear force potentials,
- [`MEHPForceBalance`](../../api/pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance): A Force Balance implementation, which uses just Hookean springs, but allows the use of slip-links to model entanglements,
- [`MEHPForceBalance2`](../../api/pylimer_tools_cpp.MEHPForceBalance2.md#pylimer_tools_cpp.MEHPForceBalance2): The faster implementation of the Force Balance procedure, without allowing for slip-links, but allowing the modelling of entanglements as static links, like tetrafunctional crosslinks.

## Sampling Entanglements

Entanglements (slip-links) can be sampled in two distinct ways:

1. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a [`Universe`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe).
2. **Analytical Brownian-bridge sampling**: Directly on the crosslinkers and strand connectivity of a [`CrosslinkerUniverse`](../../api/pylimer_tools_cpp.CrosslinkerUniverse.md#pylimer_tools_cpp.CrosslinkerUniverse) without instantiating intermediate strand beads, reducing peak memory by up to 90%+ and supporting topological leaf peeling of dangling and soluble fractions.

See [Comparing Slip-Link Sampling Methods: Discrete vs. Analytical](plot_entanglement_sampling_comparison.md) for a detailed benchmark and comparison of discrete vs. analytical sampling methods, and [Modeling Entanglements: Slip-Links vs. Slip-Springs](plot_sliplink_vs_slipspring.md) for the difference between modeling entanglements as merged slip-links versus compliant Hookean slip-springs.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="MEHPForceBalance2 is a faster implementation of the MEHP force balance method, which uses static links to model entanglements instead of slip-links. See bernhard_phantom_2025">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_balance_2_thumb.png)

[Force Balance 2](plot_force_balance_2.md)

  <div class="sphx-glr-thumbnail-title">Force Balance 2</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the MEHP implementation with slip-springs to homogenize a polymer network and predict its equilibrium shear modulus.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_balance_thumb.png)

[Force Balance](plot_force_balance.md)

  <div class="sphx-glr-thumbnail-title">Force Balance</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In this example, we study the effect of polydispersity on the shear modulus of end-linked polymer networks.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_polydispersity_study_thumb.png)

[Polydispersity Study](plot_polydispersity_study.md)

  <div class="sphx-glr-thumbnail-title">Polydispersity Study</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example shows how to use the MEHPForceBalance2 class in a deformation experiment.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_fb_deformation_exp_thumb.png)

[Deformation Experiment](plot_fb_deformation_exp.md)

  <div class="sphx-glr-thumbnail-title">Deformation Experiment</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In polymer network elasticity simulations using the Maximum Entropy Homogenization Procedure (MEHP) (gusev_numerical_2019, bernhard_phantom_2025), topological entanglements impose constraints on polymer chain fluctuations.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_sliplink_vs_slipspring_thumb.png)

[Modeling Entanglements: Slip-Links vs. Slip-Springs](plot_sliplink_vs_slipspring.md)

  <div class="sphx-glr-thumbnail-title">Modeling Entanglements: Slip-Links vs. Slip-Springs</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the &quot;original&quot; MEHP implementation (see gusev_numerical_2019) without slip-links or entanglements.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_relaxation_thumb.png)

[Maximum Entropy Homogenization Procedure (MEHP)](plot_force_relaxation.md)

  <div class="sphx-glr-thumbnail-title">Maximum Entropy Homogenization Procedure (MEHP)</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In polymer network elasticity simulations using the Maximum Entropy Homogenization Procedure (MEHP) and Force Balance (gusev_numerical_2019, bernhard_phantom_2025), entanglements are represented as slip-links or localized crosslink constraints.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_entanglement_sampling_comparison_thumb.png)

[Comparing Slip-Link Sampling Methods: Discrete vs. Analytical](plot_entanglement_sampling_comparison.md)

  <div class="sphx-glr-thumbnail-title">Comparing Slip-Link Sampling Methods: Discrete vs. Analytical</div>
</div>
<!-- thumbnail-parent-div-close --></div>
