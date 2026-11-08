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
where you find them in in the folder tests/pylimer_tools/fixtures,
in which case you will have to adjust the file paths in the examples accordingly,
or you can clone the entire repository and switch to the directory of the example script you want to run.
Then, run the script directly from there.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

## Dissipative Particle Dynamics (DPD) Simulations

pylimer-tools provides a powerful framework for simulating coarse-grained polymer systems using Dissipative Particle Dynamics (DPD) with slip-springs.
This allows for realistic modeling of polymer behavior at a mesoscopic scale.
The implementation is based on Langeloth *et al.* [[LMBohmMullerPlathe13](../acknowledgements.md#id11)], Schneider *et al.* [[SFKarimiVarzanehMullerPlathe21](../acknowledgements.md#id10)] and Schneider and De Pablo [[SDP23](../acknowledgements.md#id9)].

Refer to the full documentation of the DPD simulator in the [`DPDSimulator`](../api/pylimer_tools_cpp.DPDSimulator.md#pylimer_tools_cpp.DPDSimulator) class for more details.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

## Force Balance, Maximum Entropy Homogenization Procedures

pylimer-tools provides three implementations of the Force Balance:cite:p:bernhard_phantom_2025, Maximum Entropy Homogenization Procedure (MEHP):cite:p:gusev_numerical_2019 to reduce polymer networks to their minimum energy, maximum entropy homogenized state.
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- [`MEHPForceRelaxation`](../api/pylimer_tools_cpp.MEHPForceRelaxation.md#pylimer_tools_cpp.MEHPForceRelaxation): An MEHP implementation that allows the use of non-linear force potentials,
- [`MEHPForceBalance`](../api/pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance): A Force Balance implementation using Hookean springs that supports slip-links with optional dynamic contour slipping (or fixed links via `disable_slipping=True`),
- [`MEHPForceBalance2`](../api/pylimer_tools_cpp.MEHPForceBalance2.md#pylimer_tools_cpp.MEHPForceBalance2): The faster implementation of the Force Balance procedure, representing entanglements as fixed 4-functional links (like tetrafunctional crosslinks) or compliant Hookean springs, without contour slipping.

### Sampling Entanglements

Entanglements can be sampled in two distinct ways:

1. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a [`Universe`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe).
2. **Analytical Brownian-bridge sampling**: Directly on the crosslinkers and strand connectivity of a [`CrosslinkerUniverse`](../api/pylimer_tools_cpp.CrosslinkerUniverse.md#pylimer_tools_cpp.CrosslinkerUniverse) without instantiating intermediate strand beads, reducing peak memory by up to 90%+ and supporting topological leaf peeling of dangling and soluble fractions.

See:

- plot_entanglement_sampling_comparison for a detailed benchmark and comparison of discrete vs. analytical sampling methods,
- plot_entanglement_link_vs_spring for the difference between modeling entanglements as rigid merged links versus compliant Hookean springs in Force Balance 2,
- plot_slipping_vs_fixed_links for a comparison of dynamic contour slipping versus fixed entanglement links in Force Balance 1, and
- plot_inactive_removal_and_sampling for the removal of inactive strands/entanglements in Force Balance 2 and the impact of pre-filtering dangling chains during entanglement sampling.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

## Network Generation

pylimer-tools includes a powerful Monte Carlo-based network generator for creating realistic polymer systems.
This tools is particularly useful for generating quick initial configurations for molecular dynamics simulations,
or for doing fast predictions of the viscoelasticity using the [`NormalModeAnalyzer`](../api/pylimer_tools_cpp.NormalModeAnalyzer.md#pylimer_tools_cpp.NormalModeAnalyzer) or the
[`pylimer_tools_cpp.MEHPForceBalance`](../api/pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance) and its variants.

For more options, examples and advanced usage, see the [test suite](https://github.com/GenieTim/pylimer-tools/tree/main/tests)
and the modules API reference, particularly of the [`MCUniverseGenerator`](../api/pylimer_tools_cpp.MCUniverseGenerator.md#pylimer_tools_cpp.MCUniverseGenerator) class.

The network generation system allows you to:

- Create end-linked and vulcanized polymer networks with controllable topology
- Generate solvent chains and background molecules
- Control crosslinking degree and functionality
- Set realistic bead distances and box sizes

The [`MCUniverseGenerator`](../api/pylimer_tools_cpp.MCUniverseGenerator.md#pylimer_tools_cpp.MCUniverseGenerator) class is the main interface for network generation.
The initial method was introduced in Gusev [[Gus19](../acknowledgements.md#id8)],
the current implementation is more powerful, as needed for Bernhard and Gusev [[BG25](../acknowledgements.md#id12)].

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

## Normal Mode Analysis

pylimer-tools provides a normal mode analysis (NMA) implementation to predict the loss and storage modulus of polymer networks
as published in Gusev and Bernhard [[GB24](../acknowledgements.md#id7)].
This is useful for understanding the viscoelastic properties of the network.
The corresponding class is [`NormalModeAnalyzer`](../api/pylimer_tools_cpp.NormalModeAnalyzer.md#pylimer_tools_cpp.NormalModeAnalyzer).

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

## File I/O: Readers & Writers

pylimer-tools provides comprehensive support for reading and writing various LAMMPS file formats, enabling seamless integration with molecular dynamics workflows.

The library supports four main file I/O operations:

1. **Reading LAMMPS data files** - Complete system configurations
2. **Reading LAMMPS dump files** - Trajectory and snapshot data
3. **Reading LAMMPS output files** - Measurements and simulation results
4. **Writing LAMMPS data files** - Modified or generated systems

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

## Structure Analysis

pylimer-tools offers various tools to manipulate and analyze the topology of polymer networks.

Here you can find examples of some of the functionality provided by the package.

### Decomposing Crosslinked Networks to Chains

There are a number of ways how a given [`Universe`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) can be decomposed.

The most common use case is to decompose a crosslinked network into its constituent chains, which can then be used for further analysis or processing.

To find the chains, there are two main approaches:
- Removing the crosslinkers and analysing the remaining polymer strands using [`get_molecules()`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_molecules)
- Decomposing the network to chains, but keeping the crosslinkers in place, possibly being repeated in all the chains it is associated with. This decomposition can be done using the [`get_chains_with_crosslinker()`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_chains_with_crosslinker) method.

pylimer-tools provides two more advanced methods for decomposing networks:
- If you removed the junctions yourself, you can use [`get_clusters()`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_clusters)

If you want even more control, you could assemble your own decomposition algorithm using the [`Universe`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) class and its methods, such as [`get_atoms_connected_to()`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_atoms_connected_to).

### Finding Loops

You can also find loops in the network using the [`find_loops()`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.find_loops) method.

If you want to count the number of loops present in the network, as a function of their length, you can use the [`count_loop_lengths()`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.count_loop_lengths) method.

Finally, to find the shortest loop a specific atom is involved in, you can use the [`find_minimal_order_loop_from()`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.find_minimal_order_loop_from) method.

#### CAUTION
There are exponentially many paths between two crosslinkers of a network,
and you may run out of memory when using these functions, if your network is lattice-like.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

## Trajectory Analysis

pylimer-tools provides advanced tools for analyzing molecular dynamics trajectories, particularly for polymer systems.
This includes reading trajectory data and computing various properties.

For more information, refer to the documentation of the [`UniverseSequence`](../api/pylimer_tools_cpp.UniverseSequence.md#pylimer_tools_cpp.UniverseSequence) class and the [`Universe`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) class for a list of methods available for property computations.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>

[`Download all examples in Python source code: auto_examples_python.zip`](_downloads/07fcc19ba03226cd3d83d4e40ec44385/auto_examples_python.zip)

[`Download all examples in Jupyter notebooks: auto_examples_jupyter.zip`](_downloads/6f1e7a639e0699d6164445b55e6c116d/auto_examples_jupyter.zip)

[Gallery generated by Sphinx-Gallery](https://sphinx-gallery.github.io)
