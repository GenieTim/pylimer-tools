# Assumptions

Quite a few assumptions were made in the making of this library,
some of which have been addressed and overcome, others still remaining.
They are all found in the documentation of the respective functions.

Here, we repeat them. If your system complies with these assumptions,
you have no reason to be worried.
On the other hand, if you encounter problems or get unexpected results,
you may want to check whether your system violates one of these assumptions.

1. Network strands and their atoms are at maximum bifunctional
2. There is only ever one atom type that is a crosslink

## Periodic Boundary Conditions

This is a heads-up regarding the handling of periodic boundary conditions (PBC).
In most of the simulators and MEHP implementations of pylimer-tools, there are two implementations of the PBC.

When working with molecular dynamics simulations, especially those involving polymer networks, it is crucial to ensure that the simulation box is large enough to accommodate the periodic boundary conditions.
If the box is too small, it can lead to artifacts in the simulation results, such as unphysical interactions between atoms that are supposed to be separated by the periodic boundaries, or bonds that are shorter than they should be and flip directions due to the minimum image convention.

Particularly the latter can happen in the Force Balance, Maximum Entropy Homogenization Procedure (MEHP), where springs represent entire polymer chains.
In those cases, it is the periodic boundary conditions that prevent chains from collapsing.
If the conventional implementations of PBC are used (such as the one provided in LAMMPS), an exemplaric chain that spans the whole box would collapse, rather than staying elastically effective.
To avoid these issues, pylimer-tools defaults to *not* assuming that the simulation box is large enough for minimum-image PBC (assume_box_large_enough = False). Instead, it tracks periodic box offsets per bond or strand so that bonds or chains longer than half the box length are handled without collapsing or flipping direction.

When initializing your simulation, leave assume_box_large_enough at its default (False) if you are not confident that all bonds will remain shorter than half the box length.
Note, however, that this default mode still assumes that the initial configuration of the structure has bonds that are all shorter than half the box length, since those are used to determine the initial offsets imposed by the PBC.

Only enable the large-box assumption if you are confident that your simulation box is sufficiently large for minimum-image PBC:

- Both [`DPDSimulator`](api/pylimer_tools_cpp.DPDSimulator.md#pylimer_tools_cpp.DPDSimulator) and [`MEHPSimulator`](api/pylimer_tools_cpp.MEHPSimulator.md#pylimer_tools_cpp.MEHPSimulator) accept an `assume_box_large_enough` keyword argument (default `False`) in their constructors.
- [`pylimer_tools_cpp.DPDSimulator.config_assume_box_large_enough()`](api/pylimer_tools_cpp.DPDSimulator.md#pylimer_tools_cpp.DPDSimulator.config_assume_box_large_enough) (also aliased as [`assume_box_large_enough()`](api/pylimer_tools_cpp.DPDSimulator.md#pylimer_tools_cpp.DPDSimulator.assume_box_large_enough)) and [`pylimer_tools_cpp.MEHPSimulator.config_assume_box_large_enough()`](api/pylimer_tools_cpp.MEHPSimulator.md#pylimer_tools_cpp.MEHPSimulator.config_assume_box_large_enough) accept a boolean flag (default `True` when called explicitly; pass `False` to re-enable tracked periodic box offsets).

#### NOTE
Purely topological utilities such as [`pylimer_tools_cpp.Universe.peel_topological_leaves()`](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.peel_topological_leaves),
[`pylimer_tools_cpp.Universe.identify_active_atoms()`](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.identify_active_atoms), and
[`pylimer_tools_cpp.peel_topological_leaves()`](api/pylimer_tools_cpp.peel_topological_leaves.md#pylimer_tools_cpp.peel_topological_leaves) do not use coordinates or PBC at all:
they retain all 2-core cycles (including local non-PBC loops, which can still become
elastically active when an entanglement or slip-link is associated with them) and prune all
crosslinker-less rings (even PBC-spanning ones). See
[Topological Leaf Pruning vs. PBC-Aware Force Balance](auto_examples/structure_analysis/plot_topological_leaf_pruning_vs_pbc.md) for an
example comparing topological leaf pruning with PBC-aware MEHP force balance.
