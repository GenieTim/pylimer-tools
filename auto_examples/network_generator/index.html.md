<a id="sphx-glr-auto-examples-network-generator"></a>

# Network Generation

pylimer-tools includes a powerful Monte Carlo-based network generator for creating realistic polymer systems.
This tools is particularly useful for generating quick initial configurations for molecular dynamics simulations,
or for doing fast predictions of the viscoelasticity using the [`NormalModeAnalyzer`](../../api/pylimer_tools_cpp.NormalModeAnalyzer.md#pylimer_tools_cpp.NormalModeAnalyzer) or the
[`pylimer_tools_cpp.MEHPForceBalance`](../../api/pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance) and its variants.

For more options, examples and advanced usage, see the [test suite](https://github.com/GenieTim/pylimer-tools/tree/main/tests)
and the modules API reference, particularly of the [`MCUniverseGenerator`](../../api/pylimer_tools_cpp.MCUniverseGenerator.md#pylimer_tools_cpp.MCUniverseGenerator) class.

The network generation system allows you to:

- Create end-linked and vulcanized polymer networks with controllable topology
- Generate solvent chains and background molecules
- Control crosslinking degree and functionality
- Set realistic bead distances and box sizes

The [`MCUniverseGenerator`](../../api/pylimer_tools_cpp.MCUniverseGenerator.md#pylimer_tools_cpp.MCUniverseGenerator) class is the main interface for network generation.
The initial method was introduced in Gusev [[Gus19](../../acknowledgements.md#id8)],
the current implementation is more powerful, as needed for Bernhard and Gusev [[BG25](../../acknowledgements.md#id12)].

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
