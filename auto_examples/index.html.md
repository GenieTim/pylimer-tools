# Examples

The following gallery contains example scripts demonstrating the usage of various features in pylimer-tools.
Each script is designed to illustrate a specific functionality or analysis method available in the library.

Please be aware that the system sizes and simulation trajectory times are not representative of real-world applications.
These examples are primarily for educational purposes and to showcase the capabilities and usage of pylimer-tools.

## Running the Examples

To run the examples, you can execute the scripts directly in your Python environment.
Make sure you have pylimer-tools installed and the necessary dependencies are met.
Some examples require specific input files, such as LAMMPS data or output files.
All these files are part of the pylimer-tools test suite.
You can either download them from the [GitHub repository](https://github.com/GenieTim/pylimer-tools),
where you find them in the folder tests/pylimer_tools/fixtures,
in which case you will have to adjust the file paths in the examples accordingly,
or you can clone the entire repository and switch to the directory of the example script you want to run.
Then, run the script directly from there.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

## Dissipative Particle Dynamics (DPD) Simulations

pylimer-tools provides a powerful framework for simulating coarse-grained polymer systems using Dissipative Particle Dynamics (DPD) with slip-springs.
This allows for realistic modeling of polymer behavior at a mesoscopic scale.
The implementation is based on Langeloth *et al.* [[LMBohmMullerPlathe13](../acknowledgements.html.md#id11)], Schneider *et al.* [[SFKarimiVarzanehMullerPlathe21](../acknowledgements.html.md#id10)] and Schneider and De Pablo [[SDP23](../acknowledgements.html.md#id9)].

Refer to the full documentation of the DPD simulator in the [`DPDSimulator`](../api/pylimer_tools_cpp.DPDSimulator.html.md#pylimer_tools_cpp.DPDSimulator) class for more details.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to set up and run a DPD simulation with slip-springs using pylimer-tools.">![](auto_examples/dpd_simulations/images/thumb/sphx_glr_plot_dpd_simulation_thumb.png)

[DPD Simulations](dpd_simulations/plot_dpd_simulation.html.md)

  <div class="sphx-glr-thumbnail-title">DPD Simulations</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Force Balance, Maximum Entropy Homogenization Procedures

pylimer-tools provides a unified implementation of the Force Balance [[BG25](../acknowledgements.html.md#id12)], Maximum Entropy Homogenization Procedure (MEHP) [[Gus19](../acknowledgements.html.md#id8)] to reduce polymer networks to their minimum energy, maximum entropy homogenized state.
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

The unified architecture coordinates a canonical [`MEHPNetwork`](../api/pylimer_tools_cpp.MEHPNetwork.html.md#pylimer_tools_cpp.MEHPNetwork) (supporting phantom networks, fixed merged entanglement links via [`MERGED_LINKS`](../api/pylimer_tools_cpp.EntanglementMode.html.md#pylimer_tools_cpp.EntanglementMode.MERGED_LINKS), and compliant coupling springs via [`SPRING`](../api/pylimer_tools_cpp.EntanglementMode.html.md#pylimer_tools_cpp.EntanglementMode.SPRING)), a stateless [`MEHPForceEvaluator`](../api/pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator) ([`HookeanForceEvaluator`](../api/pylimer_tools_cpp.HookeanForceEvaluator.html.md#pylimer_tools_cpp.HookeanForceEvaluator), [`NonGaussianSpringForceEvaluator`](../api/pylimer_tools_cpp.NonGaussianSpringForceEvaluator.html.md#pylimer_tools_cpp.NonGaussianSpringForceEvaluator), or custom Python subclasses), and interchangeable [`MEHPSolver`](../api/pylimer_tools_cpp.MEHPSolver.html.md#pylimer_tools_cpp.MEHPSolver) strategies inside [`MEHPSimulator`](../api/pylimer_tools_cpp.MEHPSimulator.html.md#pylimer_tools_cpp.MEHPSimulator):

- [`LinearSystemMEHPSolver`](../api/pylimer_tools_cpp.LinearSystemMEHPSolver.html.md#pylimer_tools_cpp.LinearSystemMEHPSolver) ([`MEHPLinearSolverConfig`](../api/pylimer_tools_cpp.MEHPLinearSolverConfig.html.md#pylimer_tools_cpp.MEHPLinearSolverConfig)): The fastest solver for linear Hookean networks with fixed entanglement links ([`FIXED`](../api/pylimer_tools_cpp.SlipMode.html.md#pylimer_tools_cpp.SlipMode.FIXED)) or springs, assembling and solving a sparse graph-Laplacian linear system (including optional GPU Metal backends).
- [`SlipLinkRelaxationSolver`](../api/pylimer_tools_cpp.SlipLinkRelaxationSolver.html.md#pylimer_tools_cpp.SlipLinkRelaxationSolver) ([`MEHPSlipLinkSolverConfig`](../api/pylimer_tools_cpp.MEHPSlipLinkSolverConfig.html.md#pylimer_tools_cpp.MEHPSlipLinkSolverConfig)): An iterative alternating coordinate and contour-partition relaxation solver for Hookean networks that supports dynamic contour slipping ([`SLIDING`](../api/pylimer_tools_cpp.SlipMode.html.md#pylimer_tools_cpp.SlipMode.SLIDING)) and slip-link swapping as well as fixed links ([`FIXED`](../api/pylimer_tools_cpp.SlipMode.html.md#pylimer_tools_cpp.SlipMode.FIXED)).
- [`NonlinearOptimizationSolver`](../api/pylimer_tools_cpp.NonlinearOptimizationSolver.html.md#pylimer_tools_cpp.NonlinearOptimizationSolver) ([`MEHPNonlinearSolverConfig`](../api/pylimer_tools_cpp.MEHPNonlinearSolverConfig.html.md#pylimer_tools_cpp.MEHPNonlinearSolverConfig)): A general NLopt-based energy minimization solver supporting linear and non-linear constitutive force potentials (e.g., [`NonGaussianSpringForceEvaluator`](../api/pylimer_tools_cpp.NonGaussianSpringForceEvaluator.html.md#pylimer_tools_cpp.NonGaussianSpringForceEvaluator) or custom Python [`MEHPForceEvaluator`](../api/pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator) subclasses) on phantom and fixed-entanglement networks.

### Sampling Entanglements

Entanglements can be sampled in two ways:

1. **On-demand bounding-box sampling** (default): Computes per-strand bounding boxes and regenerates candidate strand beads on demand, matching discrete bead-level fidelity with ~80–95% lower peak memory.
2. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a [`Universe`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe).

Both methods support [`EntanglementFilteringMode`](../api/pylimer_tools_cpp.EntanglementFilteringMode.html.md#pylimer_tools_cpp.EntanglementFilteringMode) for topological leaf peeling of dangling and soluble fractions.

See:

- [Entanglement Sampling Methods](force_balance/plot_entanglement_sampling_comparison.html.md) for a detailed benchmark and comparison of discrete vs. on-demand sampling methods,
- [Entanglement Links vs. Springs](force_balance/plot_entanglement_link_vs_spring.html.md) for the difference between modeling entanglements as rigid merged links ([`MERGED_LINKS`](../api/pylimer_tools_cpp.EntanglementMode.html.md#pylimer_tools_cpp.EntanglementMode.MERGED_LINKS)) versus compliant Hookean springs ([`SPRING`](../api/pylimer_tools_cpp.EntanglementMode.html.md#pylimer_tools_cpp.EntanglementMode.SPRING)),
- [Slipping vs. Fixed Entanglement Links](force_balance/plot_slipping_vs_fixed_links.html.md) for a comparison of dynamic contour slipping ([`SLIDING`](../api/pylimer_tools_cpp.SlipMode.html.md#pylimer_tools_cpp.SlipMode.SLIDING)) versus fixed entanglement links ([`FIXED`](../api/pylimer_tools_cpp.SlipMode.html.md#pylimer_tools_cpp.SlipMode.FIXED)) using [`SlipLinkRelaxationSolver`](../api/pylimer_tools_cpp.SlipLinkRelaxationSolver.html.md#pylimer_tools_cpp.SlipLinkRelaxationSolver),
- [Inactive Removal & Sampling](force_balance/plot_inactive_removal_and_sampling.html.md) for the removal of inactive strands/entanglements via [`StructureSimplificationMode`](../api/pylimer_tools_cpp.StructureSimplificationMode.html.md#pylimer_tools_cpp.StructureSimplificationMode) and the impact of pre-filtering dangling chains during entanglement sampling,
- [Polydispersity Study](force_balance/plot_polydispersity_study.html.md) for the effect of precursor strand length polydispersity on equilibrium shear modulus and defect fractions,
- [Sampling Cut-Off Effect](force_balance/plot_sampling_cutoff_effect.html.md) for the effect of the slip-link sampling cut-off distance on the equilibrium shear modulus, investigating pre-strain distortions, link vs. spring representations, and interactions with strand length $N$, conversion $p$, and stoichiometric imbalance $r$, and
- [Representative Volume Study for PDMS Networks](force_balance/plot_representative_volume_pdms.html.md) for a representative volume element (RVE) study of realistic PDMS networks across chain length $N$, conversion $p$, stoichiometric imbalance $r$, and ensemble sample count.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="In MEHPSimulator with LinearSystemMEHPSolver (gusev_numerical_2019, bernhard_phantom_2025), users can choose between two representations for sampled entanglements via entanglement_mode:">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_entanglement_link_vs_spring_thumb.png)

[Entanglement Links vs. Springs](force_balance/plot_entanglement_link_vs_spring.html.md)

  <div class="sphx-glr-thumbnail-title">Entanglement Links vs. Springs</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In this example, we study the effect of precursor strand length polydispersity on the equilibrium shear modulus and topological network structure of end-linked polymer networks (tsimouri_monte_2021, bernhard_phantom_2025).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_polydispersity_study_thumb.png)

[Polydispersity Study](force_balance/plot_polydispersity_study.html.md)

  <div class="sphx-glr-thumbnail-title">Polydispersity Study</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In MEHP Force Balance (gusev_numerical_2019, bernhard_phantom_2025), entanglements are sampled between bead pairs within a spatial cut-off \\|\\mathbf{r}_i - \\mathbf{r}_j\\| \\le r_{\\text{cut}}.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_sampling_cutoff_effect_thumb.png)

[Sampling Cut-Off Effect](force_balance/plot_sampling_cutoff_effect.html.md)

  <div class="sphx-glr-thumbnail-title">Sampling Cut-Off Effect</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In MEHPSimulator (gusev_numerical_2019, bernhard_phantom_2025) with SlipLinkRelaxationSolver (~pylimer_tools_cpp.MEHPSlipLinkSolverConfig), entanglement links can either remain frozen at their initial contour positions (~pylimer_tools_cpp.SlipMode.FIXED) or slide dynamically along the polymer strands to balance contour tension (~pylimer_tools_cpp.SlipMode.SLIDING).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_slipping_vs_fixed_links_thumb.png)

[Slipping vs. Fixed Entanglement Links](force_balance/plot_slipping_vs_fixed_links.html.md)

  <div class="sphx-glr-thumbnail-title">Slipping vs. Fixed Entanglement Links</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Networks below full conversion contain dangling chains, sol fractions, and loops. In MEHPSimulator with fixed entanglements, links do not slip off chain ends, so two complementary controls prevent dangling arms from artificially stiffening the network:">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_inactive_removal_and_sampling_thumb.png)

[Inactive Removal & Sampling](force_balance/plot_inactive_removal_and_sampling.html.md)

  <div class="sphx-glr-thumbnail-title">Inactive Removal & Sampling</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In polymer network elasticity simulations using the Maximum Entropy Homogenization Procedure (MEHP) and Force Balance (gusev_numerical_2019, bernhard_phantom_2025), entanglements are represented as localized topological constraints (entanglement links or springs).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_entanglement_sampling_comparison_thumb.png)

[Entanglement Sampling Methods](force_balance/plot_entanglement_sampling_comparison.html.md)

  <div class="sphx-glr-thumbnail-title">Entanglement Sampling Methods</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="When predicting the equilibrium shear modulus G of polymer networks via Monte Carlo network generation (~pylimer_tools_cpp.MCUniverseGenerator), on-demand entanglement sampling (~pylimer_tools_cpp.MCUniverseGenerator.sample_entanglements), and MEHPSimulator (~pylimer_tools_cpp.LinearSystemMEHPSolver), the simulation box must be large enough to act as a Representative Volume Element (RVE)—avoiding both finite-size bias (e.g. periodic self-looping or insufficient topological sampling) and excessive sample-to-sample variance.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_representative_volume_pdms_thumb.png)

[Representative Volume Study for PDMS Networks](force_balance/plot_representative_volume_pdms.html.md)

  <div class="sphx-glr-thumbnail-title">Representative Volume Study for PDMS Networks</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example shows how to use the MEHPSimulator class in a deformation experiment.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_fb_deformation_exp_thumb.png)

[Deformation Experiment](force_balance/plot_fb_deformation_exp.html.md)

  <div class="sphx-glr-thumbnail-title">Deformation Experiment</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the Maximum Entropy Homogenization Procedure (MEHP) (gusev_numerical_2019) with slip-links to homogenize a polymer network and predict its equilibrium shear modulus.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_balance_thumb.png)

[Force Balance](force_balance/plot_force_balance.html.md)

  <div class="sphx-glr-thumbnail-title">Force Balance</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the unified MEHPSimulator with interchangeable constitutive MEHPForceEvaluator models and MEHPSolver strategies (see gusev_numerical_2019) without slip-links or entanglements.">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_relaxation_thumb.png)

[MEHP Force Relaxation](force_balance/plot_force_relaxation.html.md)

  <div class="sphx-glr-thumbnail-title">MEHP Force Relaxation</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="MEHPSimulator paired with LinearSystemMEHPSolver provides a high-performance implementation of the Maximum Entropy Homogenization Procedure (MEHP) (gusev_numerical_2019, bernhard_phantom_2025).">![](auto_examples/force_balance/images/thumb/sphx_glr_plot_force_balance_2_thumb.png)

[Force Balance 2](force_balance/plot_force_balance_2.html.md)

  <div class="sphx-glr-thumbnail-title">Force Balance 2</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Network Generation

pylimer-tools includes a powerful Monte Carlo-based network generator for creating realistic polymer systems.
This tool is particularly useful for generating quick initial configurations for molecular dynamics simulations,
or for doing fast predictions of the viscoelasticity using the [`NormalModeAnalyzer`](../api/pylimer_tools_cpp.NormalModeAnalyzer.html.md#pylimer_tools_cpp.NormalModeAnalyzer) or the
[`MEHPSimulator`](../api/pylimer_tools_cpp.MEHPSimulator.html.md#pylimer_tools_cpp.MEHPSimulator) and its solver strategies.

For more options, examples and advanced usage, see the [test suite](https://github.com/GenieTim/pylimer-tools/tree/main/tests)
and the [src](../modules.html.md) API reference, particularly of the [`MCUniverseGenerator`](../api/pylimer_tools_cpp.MCUniverseGenerator.html.md#pylimer_tools_cpp.MCUniverseGenerator) class.

The network generation system allows you to:

- Create end-linked and vulcanized polymer networks with controllable topology
- Generate solvent chains and background molecules
- Control crosslinking degree and functionality
- Set realistic bead distances and box sizes

The [`MCUniverseGenerator`](../api/pylimer_tools_cpp.MCUniverseGenerator.html.md#pylimer_tools_cpp.MCUniverseGenerator) class is the main interface for network generation.
The initial method was introduced in Gusev [[Gus19](../acknowledgements.html.md#id8)],
the current implementation is more powerful, as needed for Bernhard and Gusev [[BG25](../acknowledgements.html.md#id12)].

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example shows how to save a generated polymer network to a LAMMPS data file.">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_save_network_thumb.png)

[Save Generated Network](network_generator/plot_save_network.html.md)

  <div class="sphx-glr-thumbnail-title">Save Generated Network</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the MCUniverseGenerator to create a polydisperse end-linked polymer network with varying strand lengths.">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_polydisperse_system_thumb.png)

[Polydisperse Network](network_generator/plot_polydisperse_system.html.md)

  <div class="sphx-glr-thumbnail-title">Polydisperse Network</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Here&#x27;s a minimal example of creating a crosslinked polymer network:">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_simple_crosslinked_network_thumb.png)

[Simple Crosslinked Network](network_generator/plot_simple_crosslinked_network.html.md)

  <div class="sphx-glr-thumbnail-title">Simple Crosslinked Network</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In the following, there are a few calls showing how to fine-tune the Monte Carlo network generation process:">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_monte_carlo_configuration_thumb.png)

[Monte Carlo Configuration](network_generator/plot_monte_carlo_configuration.html.md)

  <div class="sphx-glr-thumbnail-title">Monte Carlo Configuration</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The generated networks can include topological information.">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_add_angles_thumb.png)

[Adding Topology Information](network_generator/plot_add_angles.html.md)

  <div class="sphx-glr-thumbnail-title">Adding Topology Information</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In this example, we create a vulcanized network using the pylimer-tools library.">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_vulcanized_network_thumb.png)

[Generate Vulcanized Networks](network_generator/plot_vulcanized_network.html.md)

  <div class="sphx-glr-thumbnail-title">Generate Vulcanized Networks</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Create complex systems with multiple components:">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_multi_component_system_thumb.png)

[Multi-Component System](network_generator/plot_multi_component_system.html.md)

  <div class="sphx-glr-thumbnail-title">Multi-Component System</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the link_strands_callback method to create custom crosslinking procedures as an alternative to the built-in linking methods, link_strands_to_conversion and link_strands_to_soluble_fraction. The callback function allows users to control when the crosslinking process should stop based on any criteria they can implement.">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_custom_crosslinking_callback_thumb.png)

[Custom Crosslinking Callbacks](network_generator/plot_custom_crosslinking_callback.html.md)

  <div class="sphx-glr-thumbnail-title">Custom Crosslinking Callbacks</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to create a tetra-PEG (tetra-functional polyethylene glycol) network using the add_star_crosslinkers and link_strands_to_strands_to_conversion methods.">![](auto_examples/network_generator/images/thumb/sphx_glr_plot_tetra_peg_network_thumb.png)

[Tetra-PEG Network Generation](network_generator/plot_tetra_peg_network.html.md)

  <div class="sphx-glr-thumbnail-title">Tetra-PEG Network Generation</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Normal Mode Analysis

pylimer-tools provides a normal mode analysis (NMA) implementation to predict the loss and storage modulus of polymer networks
as published in Gusev and Bernhard [[GB24](../acknowledgements.html.md#id7)].
This is useful for understanding the viscoelastic properties of the network.
The corresponding class is [`NormalModeAnalyzer`](../api/pylimer_tools_cpp.NormalModeAnalyzer.html.md#pylimer_tools_cpp.NormalModeAnalyzer).

See [KPM vs. Full Eigenvalue Normal Mode Analysis for PDMS Networks](normal_mode_analysis/plot_kpm_vs_exact_pdms.html.md) for a large-scale comparison of the Kernel Polynomial Method (KPM)
against full eigenvalue decomposition on realistic PDMS networks generated via Monte Carlo with
on-demand entanglement sampling.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the NormalModeAnalyzer in pylimer-tools to predict the loss and storage modulus of polymer networks using exact eigenvalue solvers or approximate spectral methods (~pylimer_tools_cpp.KPMConfig, SLQConfig, and HybridConfig).">![](auto_examples/normal_mode_analysis/images/thumb/sphx_glr_plot_normal_mode_analysis_thumb.png)

[Normal Mode Analysis (NMA)](normal_mode_analysis/plot_normal_mode_analysis.html.md)

  <div class="sphx-glr-thumbnail-title">Normal Mode Analysis (NMA)</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Normal Mode Analysis (NMA) in NormalModeAnalyzer predicts the linear viscoelastic storage modulus G&#x27;(\\omega) and loss modulus G&#x27;&#x27;(\\omega) of entangled polymer networks from their bead-spring connectivity graph (gusev_molecular_2024). Computing the full eigenvalue spectrum scales as \\mathcal{O}(N_{\\text{atoms}}^3), which quickly becomes prohibitive for realistic representative volume elements (RVEs). The Kernel Polynomial Method (KPM) approximates the spectral density using Chebyshev moments and stochastic trace estimation in \\mathcal{O}(N_{\\text{atoms}}) time.">![](auto_examples/normal_mode_analysis/images/thumb/sphx_glr_plot_kpm_vs_exact_pdms_thumb.png)

[KPM vs. Full Eigenvalue Normal Mode Analysis for PDMS Networks](normal_mode_analysis/plot_kpm_vs_exact_pdms.html.md)

  <div class="sphx-glr-thumbnail-title">KPM vs. Full Eigenvalue Normal Mode Analysis for PDMS Networks</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## File I/O: Readers & Writers

pylimer-tools provides comprehensive support for reading and writing various LAMMPS file formats, enabling seamless integration with molecular dynamics workflows.

The library supports four main file I/O operations:

1. **Reading LAMMPS data files** - Complete system configurations
2. **Reading LAMMPS dump files** - Trajectory and snapshot data
3. **Reading LAMMPS output files** - Measurements and simulation results
4. **Writing LAMMPS data files** - Modified or generated systems

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="Read output from pylimer-tools&#x27; own simulators.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_pylimer_tools_output_reader_thumb.png)

[pylimer-tools Output Reader](readers_writers/plot_pylimer_tools_output_reader.html.md)

  <div class="sphx-glr-thumbnail-title">pylimer-tools Output Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Process histogram data from fix ave/hist commands.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_histogram_reader_thumb.png)

[Histogram Reader](readers_writers/plot_histogram_reader.html.md)

  <div class="sphx-glr-thumbnail-title">Histogram Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="All readers support automatic caching for improved performance.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_performance_caching_thumb.png)

[Performance and Caching](readers_writers/plot_performance_caching.html.md)

  <div class="sphx-glr-thumbnail-title">Performance and Caching</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read correlation functions from fix ave/correlate and fix ave/correlate/long commands.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_correlated_averages_reader_thumb.png)

[Correlated Averages Reader](readers_writers/plot_correlated_averages_reader.html.md)

  <div class="sphx-glr-thumbnail-title">Correlated Averages Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read LAMMPS log files containing thermodynamic output from the thermo command.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_log_thermo_file_reader_thumb.png)

[Log & Thermo File Reader](readers_writers/plot_log_thermo_file_reader.html.md)

  <div class="sphx-glr-thumbnail-title">Log & Thermo File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read output from LAMMPS fix ave/time commands.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_time_averaged_data_reader_thumb.png)

[Time-Averaged Data Reader](readers_writers/plot_time_averaged_data_reader.html.md)

  <div class="sphx-glr-thumbnail-title">Time-Averaged Data Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Convert LAMMPS units to SI or other unit systems:">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_unit_conversion_thumb.png)

[Unit Conversion](readers_writers/plot_unit_conversion.html.md)

  <div class="sphx-glr-thumbnail-title">Unit Conversion</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DumpFileReader processes LAMMPS trajectory files. Is has a practical interface to read frames from a trajectory file and access them as Universe objects in the UniverseSequence class.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_dump_file_reader_thumb.png)

[Dump File Reader](readers_writers/plot_dump_file_reader.html.md)

  <div class="sphx-glr-thumbnail-title">Dump File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Take a later entry from a LAMMPS trajectory and save it as a data file.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_convert_trajectory_to_frame_thumb.png)

[Extract Trajectory Frame](readers_writers/plot_convert_trajectory_to_frame.html.md)

  <div class="sphx-glr-thumbnail-title">Extract Trajectory Frame</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="LAMMPS bond swapping can significantly speed up equilibration for polymer networks. However, to preserve the strand lengths, you need to ensure that the molecule indices are set correctly. Here&#x27;s how you could re-set the molecule indices for bond swapping with constant chain lengths.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_convert_for_bond_swapping_thumb.png)

[Convert for Bond Swapping](readers_writers/plot_convert_for_bond_swapping.html.md)

  <div class="sphx-glr-thumbnail-title">Convert for Bond Swapping</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The UniverseSequence class provides automatic memory management.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_memory_efficient_dump_file_reading_thumb.png)

[Streaming Dump File Reader](readers_writers/plot_memory_efficient_dump_file_reading.html.md)

  <div class="sphx-glr-thumbnail-title">Streaming Dump File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DataFileParser handles LAMMPS data files containing complete system definitions. The UniverseSequence provides a convenient way to have these data files read into the Universe objects.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_data_file_reader_thumb.png)

[Data File Reader](readers_writers/plot_data_file_reader.html.md)

  <div class="sphx-glr-thumbnail-title">Data File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DataFileWriter generates LAMMPS-compatible data files.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_data_file_writer_thumb.png)

[Data File Writer](readers_writers/plot_data_file_writer.html.md)

  <div class="sphx-glr-thumbnail-title">Data File Writer</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Rouse & Brownian Dynamics Simulations

pylimer-tools provides an implicit time-stepping Rouse and Brownian dynamics
simulator for polymer chains and general Hookean spring networks via
[`RouseBrownianSimulator`](../api/pylimer_tools_cpp.RouseBrownianSimulator.html.md#pylimer_tools_cpp.RouseBrownianSimulator), returning dynamical observables
in [`RouseSimulationResult`](../api/pylimer_tools_cpp.RouseSimulationResult.html.md#pylimer_tools_cpp.RouseSimulationResult).

At each time step $\Delta t$, the positions $\mathbf{r}_t \in \mathbb{R}^{N \times 3}$
are advanced using the unconditionally stable backward-Euler discretization of the
overdamped Langevin equation on the scalar graph Laplacian $\mathbf{L}$:

$$
\left(\mathbf{I} + \Delta t \frac{k}{\gamma} \mathbf{L}\right) \mathbf{r}_{t+\Delta t}
= \mathbf{r}_t + \frac{\Delta t}{\gamma} \mathbf{F}_{\text{ext}} + \sqrt{2 D \Delta t}\,\boldsymbol{\xi}_t,
$$

where $k$ is the spring constant, $\gamma$ is the bead friction coefficient,
$D = k_B T / \gamma$ is the bare bead diffusion coefficient, and
$\boldsymbol{\xi}_t$ has i.i.d. standard normal entries. All three Cartesian
coordinates share a single sparse factorization of
$\mathbf{A} = \mathbf{I} + \Delta t \frac{k}{\gamma} \mathbf{L}$, configurable
through [`RouseSolverType`](../api/pylimer_tools_cpp.RouseSolverType.html.md#pylimer_tools_cpp.RouseSolverType) (direct sparse LDLT, direct sparse
LLT, or iterative Conjugate Gradient).

See:

- [Rouse & Brownian Dynamics Simulation](rouse_simulations/plot_rouse_brownian_simulation.html.md) for
  deterministic normal-mode relaxation, linear solver comparison, external-force
  drift and creep/recoil, and stochastic Rouse MSD and end-to-end relaxation dynamics
  (including multi-tau analysis with [`Correlator`](../api/pylimer_tools_cpp.Correlator.html.md#pylimer_tools_cpp.Correlator)).
- [Polymer Architectures & Resistance-Fluctuation Duality in Rouse Dynamics](rouse_simulations/plot_rouse_topologies_and_resistance.html.md) for
  dynamics across linear, ring, star, complete, and weighted topologies and the exact
  fluctuation-resistance duality connecting [`RouseBrownianSimulator`](../api/pylimer_tools_cpp.RouseBrownianSimulator.html.md#pylimer_tools_cpp.RouseBrownianSimulator)
  and [`ResistanceDistanceAnalysis`](../api/pylimer_tools_cpp.ResistanceDistanceAnalysis.html.md#pylimer_tools_cpp.ResistanceDistanceAnalysis).

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example bridges RouseBrownianSimulator and ResistanceDistanceAnalysis, illustrating how topological effective resistances govern the equilibrium conformational fluctuations of Gaussian polymer architectures and crosslinked networks.">![](auto_examples/rouse_simulations/images/thumb/sphx_glr_plot_rouse_topologies_and_resistance_thumb.png)

[Polymer Architectures & Resistance-Fluctuation Duality in Rouse Dynamics](rouse_simulations/plot_rouse_topologies_and_resistance.html.md)

  <div class="sphx-glr-thumbnail-title">Polymer Architectures & Resistance-Fluctuation Duality in Rouse Dynamics</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use RouseBrownianSimulator to simulate overdamped polymer chain dynamics with an implicit backward-Euler integrator:">![](auto_examples/rouse_simulations/images/thumb/sphx_glr_plot_rouse_brownian_simulation_thumb.png)

[Rouse & Brownian Dynamics Simulation](rouse_simulations/plot_rouse_brownian_simulation.html.md)

  <div class="sphx-glr-thumbnail-title">Rouse & Brownian Dynamics Simulation</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Structure Analysis

pylimer-tools offers various tools to manipulate and analyze the topology of polymer networks.

Here you can find examples of some of the functionality provided by the package.

### Decomposing Crosslinked Networks to Chains

There are a number of ways how a given [`Universe`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe) can be decomposed.

The most common use case is to decompose a crosslinked network into its constituent chains, which can then be used for further analysis or processing.

To find the chains, there are two main approaches:
- Removing the crosslinkers and analysing the remaining polymer strands using [`get_molecules()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.get_molecules)
- Decomposing the network to chains, but keeping the crosslinkers in place, possibly being repeated in all the chains it is associated with. This decomposition can be done using the [`get_chains_with_crosslinker()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.get_chains_with_crosslinker) method.

pylimer-tools also provides a more advanced method for decomposing networks:
- If you removed the junctions yourself, you can use [`get_clusters()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.get_clusters)

If you want even more control, you could assemble your own decomposition algorithm using the [`Universe`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe) class and its methods, such as [`get_atoms_connected_to()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.get_atoms_connected_to).

### Finding Loops

You can also find loops in the network using the [`find_loops()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.find_loops) method.

If you want to count the number of loops present in the network, as a function of their length, you can use the [`count_loop_lengths()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.count_loop_lengths) method.

Finally, to find the shortest loop a specific atom is involved in, you can use the [`find_minimal_order_loop_from()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.find_minimal_order_loop_from) method.

#### CAUTION
There are exponentially many paths between two crosslinkers of a network,
and you may run out of memory when using these functions, if your network is lattice-like.

### Topological Leaf Pruning & Active Core

To identify the topologically active 2-core of a network (pruning dangling branches,
1-core trees, and unconnected/free chains), you can use
[`peel_topological_leaves()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.peel_topological_leaves) or
[`identify_active_atoms()`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe.identify_active_atoms).

#### NOTE
Topological leaf pruning operates purely on crosslinker connectivity
(`strand_from` / `strand_to`) and ignores atom coordinates and periodic
boundary conditions (PBC). As a result, it cannot distinguish a local
(non-PBC, zero-tension) loop from a PBC-spanning (elastically active) loop.
(Note that in an entangled network, even a non-PBC primary loop can become
elastically active if an entanglement or slip-link is associated with it.)
See [Topological Leaf Pruning vs. PBC-Aware Force Balance](structure_analysis/plot_topological_leaf_pruning_vs_pbc.html.md)
for a side-by-side comparison with PBC- and/or entanglement-aware force-balance relaxation
([`MEHPSimulator`](../api/pylimer_tools_cpp.MEHPSimulator.html.md#pylimer_tools_cpp.MEHPSimulator)).

### Topological Resistance Distance & Commute Times

For cyclic and cross-linked networks where many parallel pathways connect pairs of nodes,
[`ResistanceDistanceAnalysis`](../api/pylimer_tools_cpp.ResistanceDistanceAnalysis.html.md#pylimer_tools_cpp.ResistanceDistanceAnalysis) computes effective resistances
$R_{ij}$, the Kirchhoff index $K(G)$, the Moore-Penrose pseudo-inverse
$\mathbf{L}^+$, random-walk commute times $C_{ij}$, and node resistance
centralities directly from a [`Universe`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe), sparse graph
Laplacian matrix, or weighted edge list.

See [Topological Resistance Distance & Commute Time Analysis](structure_analysis/plot_resistance_distance_analysis.html.md) for
analytical verifications, matrix diagnostics, and performance benchmarks on crosslinked networks.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to decompose a crosslinked polymer network into chains and analyze its structure using pylimer-tools.">![](auto_examples/structure_analysis/images/thumb/sphx_glr_plot_decompose_structure_thumb.png)

[Decompose/Split Structure](structure_analysis/plot_decompose_structure.html.md)

  <div class="sphx-glr-thumbnail-title">Decompose/Split Structure</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to find and analyze loops in polymer networks using pylimer-tools.">![](auto_examples/structure_analysis/images/thumb/sphx_glr_plot_loop_analysis_thumb.png)

[Loop Finding and Analysis](structure_analysis/plot_loop_analysis.html.md)

  <div class="sphx-glr-thumbnail-title">Loop Finding and Analysis</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to analyze the end-to-end distribution of polymer chains using pylimer-tools. It reads a structure file, computes the end-to-end distances, and plots the distribution of these distances.">![](auto_examples/structure_analysis/images/thumb/sphx_glr_plot_end_to_end_distribution_thumb.png)

[End-to-End Distribution](structure_analysis/plot_end_to_end_distribution.html.md)

  <div class="sphx-glr-thumbnail-title">End-to-End Distribution</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to analyze the structure of a crosslinked polymer network to get the stoichiometric imbalance r, the crosslinker conversion p or the functionality of the crosslinker f.">![](auto_examples/structure_analysis/images/thumb/sphx_glr_plot_structure_analysis_thumb.png)

[Extract Synthesis Parameters](structure_analysis/plot_structure_analysis.html.md)

  <div class="sphx-glr-thumbnail-title">Extract Synthesis Parameters</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="peel_topological_leaves and identify_active_atoms perform an iterative 2-core peeling strictly on crosslinker connectivity (\`\`strand_from\`\` / strand_to). Because this peeling does not inspect atom coordinates, box dimensions, or periodic boundary conditions (PBC), it cannot distinguish a loop that wraps across the PBC (elastically active) from a loop that closes locally inside the box (elastically inactive / solvent).">![](auto_examples/structure_analysis/images/thumb/sphx_glr_plot_topological_leaf_pruning_vs_pbc_thumb.png)

[Topological Leaf Pruning vs. PBC-Aware Force Balance](structure_analysis/plot_topological_leaf_pruning_vs_pbc.html.md)

  <div class="sphx-glr-thumbnail-title">Topological Leaf Pruning vs. PBC-Aware Force Balance</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use ResistanceDistanceAnalysis to analyze the connectivity, effective resistances, Kirchhoff index, random-walk commute times, and node resistance centrality of polymer architectures and crosslinked networks.">![](auto_examples/structure_analysis/images/thumb/sphx_glr_plot_resistance_distance_analysis_thumb.png)

[Topological Resistance Distance & Commute Time Analysis](structure_analysis/plot_resistance_distance_analysis.html.md)

  <div class="sphx-glr-thumbnail-title">Topological Resistance Distance & Commute Time Analysis</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Trajectory Analysis

pylimer-tools provides advanced tools for analyzing molecular dynamics trajectories, particularly for polymer systems.
This includes reading trajectory data and computing various properties.

For more information, refer to the documentation of the [`UniverseSequence`](../api/pylimer_tools_cpp.UniverseSequence.html.md#pylimer_tools_cpp.UniverseSequence) class and the [`Universe`](../api/pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe) class for a list of methods available for property computations.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to analyze a polymer trajectory using pylimer-tools, including mean square displacement (MSD) and radius of gyration calculation.">![](auto_examples/trajectory_analysis/images/thumb/sphx_glr_plot_trajectory_analysis_thumb.png)

[Trajectory Analysis](trajectory_analysis/plot_trajectory_analysis.html.md)

  <div class="sphx-glr-thumbnail-title">Trajectory Analysis</div>
</div>
<!-- thumbnail-parent-div-close --></div>

[`Download all examples in Python source code: auto_examples_python.zip`](_downloads/07fcc19ba03226cd3d83d4e40ec44385/auto_examples_python.zip)

[`Download all examples in Jupyter notebooks: auto_examples_jupyter.zip`](_downloads/6f1e7a639e0699d6164445b55e6c116d/auto_examples_jupyter.zip)

[Gallery generated by Sphinx-Gallery](https://sphinx-gallery.github.io)
