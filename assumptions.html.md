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
To avoid these issues, pylimer-tools provides a configuration option to assume that the simulation box is large enough to handle periodic boundary conditions correctly. This is particularly important when using the DPDSimulator or other simulation classes that rely on accurate PBC handling.

When initializing your simulation, make sure to set the assume_box_large_enough configuration option to False if you are not confident that your simulation box is sufficiently large.
This will enable the correct handling of periodic boundaries and hopefully prevent any artifacts related to PBC.
Note, however, that it assumes that the initial configuration of the structure has bonds that are all shorter than half the box length, since those are used to determine the initial offsets imposed by the PBC.

Relevant methods:
- [`pylimer_tools_cpp.DPDSimulator.assume_box_large_enough()`](api/pylimer_tools_cpp.DPDSimulator.md#pylimer_tools_cpp.DPDSimulator.assume_box_large_enough)
- [`pylimer_tools_cpp.MEHPForceRelaxation.assume_box_large_enough()`](api/pylimer_tools_cpp.MEHPForceRelaxation.md#pylimer_tools_cpp.MEHPForceRelaxation.assume_box_large_enough)
- [`pylimer_tools_cpp.MEHPForceBalance.config_assume_box_large_enough()`](api/pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance.config_assume_box_large_enough)
