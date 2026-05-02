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
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to set up and run a DPD simulation with slip-springs using pylimer-tools.">  <div class="sphx-glr-thumbnail-title">Dissipative Particle Dynamics (DPD) Simulations</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Force Balance, Maximum Entropy Homogenization Procedures

pylimer-tools provides three implementations of the Force Balance:cite:p:bernhard_phantom_2025, Maximum Entropy Homogenization Procedure (MEHP):cite:p:gusev_numerical_2019 to reduce polymer networks to their minimum energy, maximum entropy homogenized state.
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- [`MEHPForceRelaxation`](../api/pylimer_tools_cpp.MEHPForceRelaxation.md#pylimer_tools_cpp.MEHPForceRelaxation): An MEHP implementation that allows to use non-linear force potentials,
- [`MEHPForceBalance`](../api/pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance): A Force Balance implementation, which uses just Hookean springs, but allows the use of slip-links to model entanglements,
- [`MEHPForceBalance2`](../api/pylimer_tools_cpp.MEHPForceBalance2.md#pylimer_tools_cpp.MEHPForceBalance2): The faster implementation of the Force Balance procedure, without allowing for slip-links, but allowing the modelling of entanglements as static links, like tetrafunctional crosslinks.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="MEHPForceBalance2 is a faster implementation of the MEHP force balance method, which uses static links to model entanglements instead of slip-links. See :citebernhard_phantom_2025">  <div class="sphx-glr-thumbnail-title">Force Balance 2</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the MEHP implementation with slip-springs to homogenize a polymer network and predict its equilibrium shear modulus.">  <div class="sphx-glr-thumbnail-title">Force Balance</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In this example, we study the effect of polydispersity on the shear modulus of end-linked polymer networks.">  <div class="sphx-glr-thumbnail-title">Polydispersity Study</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example shows how to use the MEHPForceBalance2 class in a deformation experiment.">  <div class="sphx-glr-thumbnail-title">Deformation Experiment</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the &quot;original&quot; MEHP implementation (see :citegusev_numerical_2019) without slip-links or entanglements.">  <div class="sphx-glr-thumbnail-title">Maximum Entropy Homogenization Procedure (MEHP)</div>
</div>
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
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example shows how to save a generated polymer network to a LAMMPS data file.">  <div class="sphx-glr-thumbnail-title">Save Generated Network</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Here&#x27;s a minimal example of creating a crosslinked polymer network:">  <div class="sphx-glr-thumbnail-title">Simple Crosslinked Network</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The generated networks can include topological information.">  <div class="sphx-glr-thumbnail-title">Adding Topology Information</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the MCUniverseGenerator to create a polydisperse end-linked polymer network with varying strand lengths.">  <div class="sphx-glr-thumbnail-title">Generate a Polydisperse Polymer Network</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Create complex systems with multiple components:">  <div class="sphx-glr-thumbnail-title">Multi-Component System</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In the following, there are a few calls showing how to fine-tune the Monte Carlo network generation process:">  <div class="sphx-glr-thumbnail-title">Monte Carlo Configuration</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In this example, we create a vulcanized network using the pylimer-tools library.">  <div class="sphx-glr-thumbnail-title">Generate Vulcanized Networks</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the link_strands_callback method to create custom crosslinking procedures as an alternative to the built-in linking methods, link_strands_to_conversion and link_strands_to_soluble_fraction. The callback function allows users to control when the crosslinking process should stop based on any criteria they can implement.">  <div class="sphx-glr-thumbnail-title">Custom Crosslinking with Callback Functions</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to create a tetra-PEG (tetra-functional polyethylene glycol) network using the add_star_crosslinkers and link_strands_to_strands_to_conversion methods.">  <div class="sphx-glr-thumbnail-title">Tetra-PEG Network Generation</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Normal Mode Analysis

pylimer-tools provides a normal mode analysis (NMA) implementation to predict the loss and storage modulus of polymer networks
as published in Gusev and Bernhard [[GB24](../acknowledgements.md#id7)].
This is useful for understanding the viscoelastic properties of the network.
The corresponding class is [`NormalModeAnalyzer`](../api/pylimer_tools_cpp.NormalModeAnalyzer.md#pylimer_tools_cpp.NormalModeAnalyzer).

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the NormalModeAnalyzer in pylimer-tools to predict the loss and storage modulus of polymer networks.">  <div class="sphx-glr-thumbnail-title">Normal Mode Analysis (NMA)</div>
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
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="Read output from pylimer-tools&#x27; own simulators.">  <div class="sphx-glr-thumbnail-title">pylimer-tools Output Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read correlation functions from fix ave/correlate and fix ave/correlate/long commands.">  <div class="sphx-glr-thumbnail-title">Correlated Averages Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Process histogram data from fix ave/hist commands.">  <div class="sphx-glr-thumbnail-title">Histogram Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="All readers support automatic caching for improved performance.">  <div class="sphx-glr-thumbnail-title">Performance and Caching</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Convert LAMMPS units to SI or other unit systems:">  <div class="sphx-glr-thumbnail-title">Unit Conversion</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read LAMMPS log files containing thermodynamic output from the thermo command.">  <div class="sphx-glr-thumbnail-title">Log & Thermo File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Take a later entry from a LAMMPS trajectory and save it as a data file.">  <div class="sphx-glr-thumbnail-title">Extract Trajectory Frame to Data File</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DumpFileReader processes LAMMPS trajectory files. Is has a practical interface to read frames from a trajectory file and access them as Universe objects in the UniverseSequence class.">  <div class="sphx-glr-thumbnail-title">Dump File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read output from LAMMPS fix ave/time commands.">  <div class="sphx-glr-thumbnail-title">Time-Averaged Data Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="LAMMPS bond swapping can significantly speed up equilibration for polymer networks. However, to preserve the strand lengths, you need to ensure that the molecule indices are set correctly. Here&#x27;s how you could re-set the molecule indices for bond swapping with constant chain lengths.">  <div class="sphx-glr-thumbnail-title">Convert for Bond Swapping</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DataFileParser handles LAMMPS data files containing complete system definitions. The UniverseSequence provides a convenient way to have these data files read into the Universe objects.">  <div class="sphx-glr-thumbnail-title">Data File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The UniverseSequence class provides automatic memory management.">  <div class="sphx-glr-thumbnail-title">Memory-Efficient Dump File Reading</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DataFileWriter generates LAMMPS-compatible data files.">  <div class="sphx-glr-thumbnail-title">Data File Writer</div>
</div>
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

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to decompose a crosslinked polymer network into chains and analyze its structure using pylimer-tools.">  <div class="sphx-glr-thumbnail-title">Decompose/Split Structure</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to find and analyze loops in polymer networks using pylimer-tools.">  <div class="sphx-glr-thumbnail-title">Loop Finding and Analysis</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to analyze the end-to-end distribution of polymer chains using pylimer-tools. It reads a structure file, computes the end-to-end distances, and plots the distribution of these distances.">  <div class="sphx-glr-thumbnail-title">End-to-End Distribution</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to analyze the structure of a crosslinked polymer network to get the stoichiometric imbalance r, the crosslinker conversion p or the functionality of the crosslinker f.">  <div class="sphx-glr-thumbnail-title">Extract Synthesis Parameters</div>
</div>
<!-- thumbnail-parent-div-close --></div>

## Trajectory Analysis

pylimer-tools provides advanced tools for analyzing molecular dynamics trajectories, particularly for polymer systems.
This includes reading trajectory data and computing various properties.

For more information, refer to the documentation of the [`UniverseSequence`](../api/pylimer_tools_cpp.UniverseSequence.md#pylimer_tools_cpp.UniverseSequence) class and the [`Universe`](../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) class for a list of methods available for property computations.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to analyze a polymer trajectory using pylimer-tools, including mean square displacement (MSD) and radius of gyration calculation.">  <div class="sphx-glr-thumbnail-title">Trajectory Analysis</div>
</div>
<!-- thumbnail-parent-div-close --></div>
