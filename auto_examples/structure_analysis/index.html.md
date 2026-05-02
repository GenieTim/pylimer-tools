<a id="sphx-glr-auto-examples-structure-analysis"></a>

# Structure Analysis

pylimer-tools offers various tools to manipulate and analyze the topology of polymer networks.

Here you can find examples of some of the functionality provided by the package.

## Decomposing Crosslinked Networks to Chains

There are a number of ways how a given [`Universe`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) can be decomposed.

The most common use case is to decompose a crosslinked network into its constituent chains, which can then be used for further analysis or processing.

To find the chains, there are two main approaches:
- Removing the crosslinkers and analysing the remaining polymer strands using [`get_molecules()`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_molecules)
- Decomposing the network to chains, but keeping the crosslinkers in place, possibly being repeated in all the chains it is associated with. This decomposition can be done using the [`get_chains_with_crosslinker()`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_chains_with_crosslinker) method.

pylimer-tools provides two more advanced methods for decomposing networks:
- If you removed the junctions yourself, you can use [`get_clusters()`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_clusters)

If you want even more control, you could assemble your own decomposition algorithm using the [`Universe`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) class and its methods, such as [`get_atoms_connected_to()`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_atoms_connected_to).

## Finding Loops

You can also find loops in the network using the [`find_loops()`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.find_loops) method.

If you want to count the number of loops present in the network, as a function of their length, you can use the [`count_loop_lengths()`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.count_loop_lengths) method.

Finally, to find the shortest loop a specific atom is involved in, you can use the [`find_minimal_order_loop_from()`](../../api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.find_minimal_order_loop_from) method.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to decompose a crosslinked polymer network into chains and analyze its structure using pylimer-tools.">  <div class="sphx-glr-thumbnail-title">Decompose/Split Structure</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to find and analyze loops in polymer networks using pylimer-tools.">  <div class="sphx-glr-thumbnail-title">Loop Finding and Analysis</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to analyze the end-to-end distribution of polymer chains using pylimer-tools. It reads a structure file, computes the end-to-end distances, and plots the distribution of these distances.">  <div class="sphx-glr-thumbnail-title">End-to-End Distribution</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to analyze the structure of a crosslinked polymer network to get the stoichiometric imbalance r, the crosslinker conversion p or the functionality of the crosslinker f.">  <div class="sphx-glr-thumbnail-title">Extract Synthesis Parameters</div>
</div>
<!-- thumbnail-parent-div-close --></div>
