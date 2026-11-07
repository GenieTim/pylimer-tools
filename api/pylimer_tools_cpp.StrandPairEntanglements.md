# StrandPairEntanglements

### *class* pylimer_tools_cpp.StrandPairEntanglements(self: [pylimer_tools_cpp.StrandPairEntanglements](#pylimer_tools_cpp.StrandPairEntanglements))

Bases: `pybind11_object`

A struct to store analytically sampled entanglements on strands (without requiring an explicit Universe of all beads).

Get an instance of this struct

### Attributes Summary

| [`alpha_1`](#pylimer_tools_cpp.StrandPairEntanglements.alpha_1)      | Contour fraction along the first strand in (0, 1).              |
|---------------------------------------------------------------|-----------------------------------------------------------------|
| [`alpha_2`](#pylimer_tools_cpp.StrandPairEntanglements.alpha_2)      | Contour fraction along the second strand in (0, 1).             |
| [`bead_idx_1`](#pylimer_tools_cpp.StrandPairEntanglements.bead_idx_1)   | Bead index along the first strand in [0, beads_in_strand - 1].  |
| [`bead_idx_2`](#pylimer_tools_cpp.StrandPairEntanglements.bead_idx_2)   | Bead index along the second strand in [0, beads_in_strand - 1]. |
| [`position_1`](#pylimer_tools_cpp.StrandPairEntanglements.position_1)   | Sampled 3D positions of the first entangled beads.              |
| [`position_2`](#pylimer_tools_cpp.StrandPairEntanglements.position_2)   | Sampled 3D positions of the second entangled beads.             |
| [`strand_idx_1`](#pylimer_tools_cpp.StrandPairEntanglements.strand_idx_1) | First strand index of each entanglement pair.                   |
| [`strand_idx_2`](#pylimer_tools_cpp.StrandPairEntanglements.strand_idx_2) | Second strand index of each entanglement pair.                  |
| [`x`](#pylimer_tools_cpp.StrandPairEntanglements.x)            | X-coordinates of the sampled entanglement midpoints.            |
| [`y`](#pylimer_tools_cpp.StrandPairEntanglements.y)            | Y-coordinates of the sampled entanglement midpoints.            |
| [`z`](#pylimer_tools_cpp.StrandPairEntanglements.z)            | Z-coordinates of the sampled entanglement midpoints.            |

### Methods Summary

| [`size`](#pylimer_tools_cpp.StrandPairEntanglements.size)(self)                            |                                                                                                                              |
|----------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------|
| [`to_atom_pair_entanglements`](#pylimer_tools_cpp.StrandPairEntanglements.to_atom_pair_entanglements)(self, ...) | Convert the strand pair entanglements to AtomPairEntanglements using the deterministic atom indexing of CrosslinkerUniverse. |

### Attributes Documentation

#### alpha_1

Contour fraction along the first strand in (0, 1).

#### alpha_2

Contour fraction along the second strand in (0, 1).

#### bead_idx_1

Bead index along the first strand in [0, beads_in_strand - 1].

#### bead_idx_2

Bead index along the second strand in [0, beads_in_strand - 1].

#### position_1

Sampled 3D positions of the first entangled beads.

#### position_2

Sampled 3D positions of the second entangled beads.

#### strand_idx_1

First strand index of each entanglement pair.

#### strand_idx_2

Second strand index of each entanglement pair.

#### x

X-coordinates of the sampled entanglement midpoints.

#### y

Y-coordinates of the sampled entanglement midpoints.

#### z

Z-coordinates of the sampled entanglement midpoints.

### Methods Documentation

#### size(self: [pylimer_tools_cpp.StrandPairEntanglements](#pylimer_tools_cpp.StrandPairEntanglements)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### to_atom_pair_entanglements(self: [pylimer_tools_cpp.StrandPairEntanglements](#pylimer_tools_cpp.StrandPairEntanglements), crosslinker_universe: [pylimer_tools_cpp.CrosslinkerUniverse](pylimer_tools_cpp.CrosslinkerUniverse.md#pylimer_tools_cpp.CrosslinkerUniverse)) → [pylimer_tools_cpp.AtomPairEntanglements](pylimer_tools_cpp.AtomPairEntanglements.md#pylimer_tools_cpp.AtomPairEntanglements)

Convert the strand pair entanglements to AtomPairEntanglements using the deterministic atom indexing of CrosslinkerUniverse.
