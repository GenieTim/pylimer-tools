<a id="sphx-glr-auto-examples-network-generator"></a>

# Network Generation

pylimer-tools includes a powerful Monte Carlo-based network generator for creating realistic polymer systems.
This tools is particularly useful for generating quick initial configurations for molecular dynamics simulations,
or for doing fast predictions of the viscoelasticity using the [`NormalModeAnalyzer`](../../api/pylimer_tools_cpp.NormalModeAnalyzer.html.md#pylimer_tools_cpp.NormalModeAnalyzer) or the
[`pylimer_tools_cpp.MEHPForceBalance`](../../api/pylimer_tools_cpp.MEHPForceBalance.html.md#pylimer_tools_cpp.MEHPForceBalance) and its variants.

For more options, examples and advanced usage, see the [test suite](https://github.com/GenieTim/pylimer-tools/tree/main/tests)
and the modules API reference, particularly of the [`MCUniverseGenerator`](../../api/pylimer_tools_cpp.MCUniverseGenerator.html.md#pylimer_tools_cpp.MCUniverseGenerator) class.

The network generation system allows you to:

- Create end-linked and vulcanized polymer networks with controllable topology
- Generate solvent chains and background molecules
- Control crosslinking degree and functionality
- Set realistic bead distances and box sizes

The [`MCUniverseGenerator`](../../api/pylimer_tools_cpp.MCUniverseGenerator.html.md#pylimer_tools_cpp.MCUniverseGenerator) class is the main interface for network generation.
The initial method was introduced in Gusev [[Gus19](../../acknowledgements.html.md#id8)],
the current implementation is more powerful, as needed for Bernhard and Gusev [[BG25](../../acknowledgements.html.md#id12)].

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open -->
<!-- thumbnail-parent-div-close --></div>
