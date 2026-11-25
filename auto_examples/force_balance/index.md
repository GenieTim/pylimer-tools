<a id="sphx-glr-auto-examples-force-balance"></a>

# Force Balance, Maximum Entropy Homogenization Procedures

pylimer-tools provides a unified implementation of the Force Balance [[BG25](../../acknowledgements.md#id12)], Maximum Entropy Homogenization Procedure (MEHP) [[Gus19](../../acknowledgements.md#id8)] to reduce polymer networks to their minimum energy, maximum entropy homogenized state.
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

The unified architecture coordinates a canonical [`MEHPNetwork`](../../api/pylimer_tools_cpp.MEHPNetwork.md#pylimer_tools_cpp.MEHPNetwork) (supporting phantom networks, fixed merged entanglement links via [`MERGED_LINKS`](../../api/pylimer_tools_cpp.EntanglementMode.md#pylimer_tools_cpp.EntanglementMode.MERGED_LINKS), and compliant coupling springs via [`SPRING`](../../api/pylimer_tools_cpp.EntanglementMode.md#pylimer_tools_cpp.EntanglementMode.SPRING)), a stateless [`MEHPForceEvaluator`](../../api/pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator) ([`HookeanForceEvaluator`](../../api/pylimer_tools_cpp.HookeanForceEvaluator.md#pylimer_tools_cpp.HookeanForceEvaluator), [`NonGaussianSpringForceEvaluator`](../../api/pylimer_tools_cpp.NonGaussianSpringForceEvaluator.md#pylimer_tools_cpp.NonGaussianSpringForceEvaluator), or custom Python subclasses), and interchangeable [`MEHPSolver`](../../api/pylimer_tools_cpp.MEHPSolver.md#pylimer_tools_cpp.MEHPSolver) strategies inside [`MEHPSimulator`](../../api/pylimer_tools_cpp.MEHPSimulator.md#pylimer_tools_cpp.MEHPSimulator):

- [`LinearSystemMEHPSolver`](../../api/pylimer_tools_cpp.LinearSystemMEHPSolver.md#pylimer_tools_cpp.LinearSystemMEHPSolver) ([`MEHPLinearSolverConfig`](../../api/pylimer_tools_cpp.MEHPLinearSolverConfig.md#pylimer_tools_cpp.MEHPLinearSolverConfig)): The fastest solver for linear Hookean networks with fixed entanglement links ([`FIXED`](../../api/pylimer_tools_cpp.SlipMode.md#pylimer_tools_cpp.SlipMode.FIXED)) or springs, assembling and solving a sparse graph-Laplacian linear system (including optional GPU Metal backends).
- [`SlipLinkRelaxationSolver`](../../api/pylimer_tools_cpp.SlipLinkRelaxationSolver.md#pylimer_tools_cpp.SlipLinkRelaxationSolver) ([`MEHPSlipLinkSolverConfig`](../../api/pylimer_tools_cpp.MEHPSlipLinkSolverConfig.md#pylimer_tools_cpp.MEHPSlipLinkSolverConfig)): An iterative alternating coordinate and contour-partition relaxation solver for Hookean networks that supports dynamic contour slipping ([`SLIDING`](../../api/pylimer_tools_cpp.SlipMode.md#pylimer_tools_cpp.SlipMode.SLIDING)) and slip-link swapping as well as fixed links ([`FIXED`](../../api/pylimer_tools_cpp.SlipMode.md#pylimer_tools_cpp.SlipMode.FIXED)).
- [`NonlinearOptimizationSolver`](../../api/pylimer_tools_cpp.NonlinearOptimizationSolver.md#pylimer_tools_cpp.NonlinearOptimizationSolver) ([`MEHPNonlinearSolverConfig`](../../api/pylimer_tools_cpp.MEHPNonlinearSolverConfig.md#pylimer_tools_cpp.MEHPNonlinearSolverConfig)): A general NLopt-based energy minimization solver supporting linear and non-linear constitutive force potentials (e.g., [`NonGaussianSpringForceEvaluator`](../../api/pylimer_tools_cpp.NonGaussianSpringForceEvaluator.md#pylimer_tools_cpp.NonGaussianSpringForceEvaluator) or custom Python [`MEHPForceEvaluator`](../../api/pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator) subclasses) on phantom and fixed-entanglement networks.

## Sampling Entanglements

Entanglements can be sampled in two ways:

1. **On-demand bounding-box sampling** (default): Computes per-strand bounding boxes and regenerates candidate strand beads on demand, matching discrete bead-level fidelity with ~80–95% lower peak memory.
2. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a [`Universe`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe).

Both methods support [`EntanglementFilteringMode`](../../api/pylimer_tools_cpp.EntanglementFilteringMode.md#pylimer_tools_cpp.EntanglementFilteringMode) for topological leaf peeling of dangling and soluble fractions.

See:

- [Entanglement Sampling Methods](plot_entanglement_sampling_comparison.md) for a detailed benchmark and comparison of discrete vs. on-demand sampling methods,
- [Entanglement Links vs. Springs](plot_entanglement_link_vs_spring.md) for the difference between modeling entanglements as rigid merged links ([`MERGED_LINKS`](../../api/pylimer_tools_cpp.EntanglementMode.md#pylimer_tools_cpp.EntanglementMode.MERGED_LINKS)) versus compliant Hookean springs ([`SPRING`](../../api/pylimer_tools_cpp.EntanglementMode.md#pylimer_tools_cpp.EntanglementMode.SPRING)),
- [Slipping vs. Fixed Entanglement Links](plot_slipping_vs_fixed_links.md) for a comparison of dynamic contour slipping ([`SLIDING`](../../api/pylimer_tools_cpp.SlipMode.md#pylimer_tools_cpp.SlipMode.SLIDING)) versus fixed entanglement links ([`FIXED`](../../api/pylimer_tools_cpp.SlipMode.md#pylimer_tools_cpp.SlipMode.FIXED)) using [`SlipLinkRelaxationSolver`](../../api/pylimer_tools_cpp.SlipLinkRelaxationSolver.md#pylimer_tools_cpp.SlipLinkRelaxationSolver),
- [Inactive Removal & Sampling](plot_inactive_removal_and_sampling.md) for the removal of inactive strands/entanglements via [`StructureSimplificationMode`](../../api/pylimer_tools_cpp.StructureSimplificationMode.md#pylimer_tools_cpp.StructureSimplificationMode) and the impact of pre-filtering dangling chains during entanglement sampling,
- [Polydispersity Study](plot_polydispersity_study.md) for the effect of precursor strand length polydispersity on equilibrium shear modulus and defect fractions,
- [Sampling Cut-Off Effect](plot_sampling_cutoff_effect.md) for the effect of the slip-link sampling cut-off distance on the equilibrium shear modulus, investigating pre-strain distortions, link vs. spring representations, and interactions with strand length $N$, conversion $p$, and stoichiometric imbalance $r$, and
- [Representative Volume Study for PDMS Networks](plot_representative_volume_pdms.md) for a representative volume element (RVE) study of realistic PDMS networks across chain length $N$, conversion $p$, stoichiometric imbalance $r$, and ensemble sample count.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="In MEHPSimulator with LinearSystemMEHPSolver (gusev_numerical_2019, bernhard_phantom_2025), users can choose between two representations for sampled entanglements via entanglement_mode:">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_entanglement_link_vs_spring_thumb.png)

[Entanglement Links vs. Springs](plot_entanglement_link_vs_spring.md)

  <div class="sphx-glr-thumbnail-title">Entanglement Links vs. Springs</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In this example, we study the effect of precursor strand length polydispersity on the equilibrium shear modulus and topological network structure of end-linked polymer networks (tsimouri_monte_2021, bernhard_phantom_2025).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_polydispersity_study_thumb.png)

[Polydispersity Study](plot_polydispersity_study.md)

  <div class="sphx-glr-thumbnail-title">Polydispersity Study</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In MEHP Force Balance (gusev_numerical_2019, bernhard_phantom_2025), entanglements are sampled between bead pairs within a spatial cut-off \\|\\mathbf{r}_i - \\mathbf{r}_j\\| \\le r_{\\text{cut}}.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_sampling_cutoff_effect_thumb.png)

[Sampling Cut-Off Effect](plot_sampling_cutoff_effect.md)

  <div class="sphx-glr-thumbnail-title">Sampling Cut-Off Effect</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In MEHPSimulator (gusev_numerical_2019, bernhard_phantom_2025) with SlipLinkRelaxationSolver (~pylimer_tools_cpp.MEHPSlipLinkSolverConfig), entanglement links can either remain frozen at their initial contour positions (~pylimer_tools_cpp.SlipMode.FIXED) or slide dynamically along the polymer strands to balance contour tension (~pylimer_tools_cpp.SlipMode.SLIDING).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_slipping_vs_fixed_links_thumb.png)

[Slipping vs. Fixed Entanglement Links](plot_slipping_vs_fixed_links.md)

  <div class="sphx-glr-thumbnail-title">Slipping vs. Fixed Entanglement Links</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Networks below full conversion contain dangling chains, sol fractions, and loops. In MEHPSimulator with fixed entanglements, links do not slip off chain ends, so two complementary controls prevent dangling arms from artificially stiffening the network:">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_inactive_removal_and_sampling_thumb.png)

[Inactive Removal & Sampling](plot_inactive_removal_and_sampling.md)

  <div class="sphx-glr-thumbnail-title">Inactive Removal & Sampling</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In polymer network elasticity simulations using the Maximum Entropy Homogenization Procedure (MEHP) and Force Balance (gusev_numerical_2019, bernhard_phantom_2025), entanglements are represented as localized topological constraints (entanglement links or springs).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_entanglement_sampling_comparison_thumb.png)

[Entanglement Sampling Methods](plot_entanglement_sampling_comparison.md)

  <div class="sphx-glr-thumbnail-title">Entanglement Sampling Methods</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="When predicting the equilibrium shear modulus G of polymer networks via Monte Carlo network generation (~pylimer_tools_cpp.MCUniverseGenerator), on-demand entanglement sampling (~pylimer_tools_cpp.MCUniverseGenerator.sample_entanglements), and MEHPSimulator (~pylimer_tools_cpp.LinearSystemMEHPSolver), the simulation box must be large enough to act as a Representative Volume Element (RVE)—avoiding both finite-size bias (e.g. periodic self-looping or insufficient topological sampling) and excessive sample-to-sample variance.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_representative_volume_pdms_thumb.png)

[Representative Volume Study for PDMS Networks](plot_representative_volume_pdms.md)

  <div class="sphx-glr-thumbnail-title">Representative Volume Study for PDMS Networks</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example shows how to use the MEHPSimulator class in a deformation experiment.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_fb_deformation_exp_thumb.png)

[Deformation Experiment](plot_fb_deformation_exp.md)

  <div class="sphx-glr-thumbnail-title">Deformation Experiment</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the Maximum Entropy Homogenization Procedure (MEHP) (gusev_numerical_2019) with slip-links to homogenize a polymer network and predict its equilibrium shear modulus.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_balance_thumb.png)

[Force Balance](plot_force_balance.md)

  <div class="sphx-glr-thumbnail-title">Force Balance</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the unified MEHPSimulator with interchangeable constitutive MEHPForceEvaluator models and MEHPSolver strategies (see gusev_numerical_2019) without slip-links or entanglements.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_relaxation_thumb.png)

[MEHP Force Relaxation](plot_force_relaxation.md)

  <div class="sphx-glr-thumbnail-title">MEHP Force Relaxation</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="MEHPSimulator paired with LinearSystemMEHPSolver provides a high-performance implementation of the Maximum Entropy Homogenization Procedure (MEHP) (gusev_numerical_2019, bernhard_phantom_2025).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_balance_2_thumb.png)

[Force Balance 2](plot_force_balance_2.md)

  <div class="sphx-glr-thumbnail-title">Force Balance 2</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Miller--Macosko recursive probability theory (MMT, macosko_new_1976, miller_new_1976a, miller_calculation_1979) provides closed-form mean-field predictions for the gelation threshold, soluble/dangling/active weight fractions, and equilibrium shear modulus of end-linked polymer networks (gusev_molecular_2022, tsimouri_comparison_2024). Because MMT treats the network as an infinite, loop-free Bethe tree, comparing its predictions against explicit 3D spatial networks generated with MCUniverseGenerator and evaluated with MEHPSimulator (gusev_numerical_2019, bernhard_phantom_2025, bernhard_pylimertools_2025) reveals both where mean-field theory excels and where spatial loops, dilution, and monofunctional dilution cause systematic deviations.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_mmt_vs_mehp_gelation_and_loops_thumb.png)

[Mean-Field Theory (MMT) vs. Spatial MC + MEHP Across Gelation & Loop Regimes](plot_mmt_vs_mehp_gelation_and_loops.md)

  <div class="sphx-glr-thumbnail-title">Mean-Field Theory (MMT) vs. Spatial MC + MEHP Across Gelation & Loop Regimes</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example investigates finite-chain extensibility and large-strain stiffening in end-linked PDMS networks using NonGaussianSpringForceEvaluator versus the Gaussian HookeanForceEvaluator inside the unified MEHPSimulator framework (gusev_numerical_2019, bernhard_phantom_2025, bernhard_pylimertools_2025).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_finite_extensibility_deformation_thumb.png)

[Finite Extensibility & Large-Strain Stiffening](plot_finite_extensibility_deformation.md)

  <div class="sphx-glr-thumbnail-title">Finite Extensibility & Large-Strain Stiffening</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Real polymer elastomers span more than an order of magnitude in entanglement plateau modulus (G_e = 0.25\\text{--}4.38\\,\\text{MPa}) and entanglement molar mass (M_e = \\rho R T / G_e \\approx 480\\text{--}13{,}300\\,\\text{g/mol}). Using the Everaers et al. polymer property database bundled in pylimer_tools.io.bead_spring_parameter_provider together with MCUniverseGenerator and MEHPSimulator (gusev_numerical_2019, bernhard_phantom_2025, bernhard_pylimertools_2025), this example investigates two fundamental questions for predictive network modeling:">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_multichemistry_and_coarse_graining_thumb.png)

[Multi-Chemistry Crossover & Coarse-Graining Invariance (GAUSSIAN vs. KUHN vs. KG_LJ)](plot_multichemistry_and_coarse_graining.md)

  <div class="sphx-glr-thumbnail-title">Multi-Chemistry Crossover & Coarse-Graining Invariance (GAUSSIAN vs. KUHN vs. KG_LJ)</div>
</div>
<!-- thumbnail-parent-div-close --></div>
