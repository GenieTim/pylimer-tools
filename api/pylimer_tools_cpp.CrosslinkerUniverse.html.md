# CrosslinkerUniverse

### *class* pylimer_tools_cpp.CrosslinkerUniverse(self: [pylimer_tools_cpp.CrosslinkerUniverse](#pylimer_tools_cpp.CrosslinkerUniverse))

Bases: `pybind11_object`

A struct representing the topological and spatial state of crosslinkers and strands
without explicitly instantiating all individual strand beads.

Get an instance of this struct

### Attributes Summary

| [`bead_distance_in_strand`](#pylimer_tools_cpp.CrosslinkerUniverse.bead_distance_in_strand)              |    |
|---------------------------------------------------------------------------------------|----|
| [`beads_in_strand`](#pylimer_tools_cpp.CrosslinkerUniverse.beads_in_strand)                      |    |
| [`mean_squared_bead_distance_in_strand`](#pylimer_tools_cpp.CrosslinkerUniverse.mean_squared_bead_distance_in_strand) |    |
| [`strand_bead_type`](#pylimer_tools_cpp.CrosslinkerUniverse.strand_bead_type)                     |    |
| [`strand_from`](#pylimer_tools_cpp.CrosslinkerUniverse.strand_from)                          |    |
| [`strand_seeds`](#pylimer_tools_cpp.CrosslinkerUniverse.strand_seeds)                         |    |
| [`strand_to`](#pylimer_tools_cpp.CrosslinkerUniverse.strand_to)                            |    |
| [`strands_of_xlink`](#pylimer_tools_cpp.CrosslinkerUniverse.strands_of_xlink)                     |    |
| [`xlink_chain_id`](#pylimer_tools_cpp.CrosslinkerUniverse.xlink_chain_id)                       |    |
| [`xlink_types`](#pylimer_tools_cpp.CrosslinkerUniverse.xlink_types)                          |    |
| [`xlink_x`](#pylimer_tools_cpp.CrosslinkerUniverse.xlink_x)                              |    |
| [`xlink_y`](#pylimer_tools_cpp.CrosslinkerUniverse.xlink_y)                              |    |
| [`xlink_z`](#pylimer_tools_cpp.CrosslinkerUniverse.xlink_z)                              |    |

### Methods Summary

| [`classify_strand`](#pylimer_tools_cpp.CrosslinkerUniverse.classify_strand)(self, strand_idx)                | Classify a strand into its topological StrandCategory.                                    |
|---------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| [`compute_bond_length_stats`](#pylimer_tools_cpp.CrosslinkerUniverse.compute_bond_length_stats)(self, box[, ...])      | Stream over all strands in O(1) strand memory to compute realized bond length statistics. |
| [`compute_end_to_end_stats`](#pylimer_tools_cpp.CrosslinkerUniverse.compute_end_to_end_stats)(self, box[, ...])       | Stream over strands matching category to compute realized end-to-end distance statistics. |
| [`compute_strand_bond_counts`](#pylimer_tools_cpp.CrosslinkerUniverse.compute_strand_bond_counts)(self[, category])     | Return the number of bonds for each strand matching category.                             |
| [`count_atom_types`](#pylimer_tools_cpp.CrosslinkerUniverse.count_atom_types)(self)                           | Count atoms per atom type in O(N_xlinks + N_strands).                                     |
| [`count_atoms_by_degree`](#pylimer_tools_cpp.CrosslinkerUniverse.count_atoms_by_degree)(self)                      | Count atoms per bond degree in O(N_xlinks + N_strands).                                   |
| [`determine_effective_functionality_per_type`](#pylimer_tools_cpp.CrosslinkerUniverse.determine_effective_functionality_per_type)(self) | Compute effective functionality (mean degree) per atom type in O(N_xlinks + N_strands).   |

### Attributes Documentation

#### bead_distance_in_strand

#### beads_in_strand

#### mean_squared_bead_distance_in_strand

#### strand_bead_type

#### strand_from

#### strand_seeds

#### strand_to

#### strands_of_xlink

#### xlink_chain_id

#### xlink_types

#### xlink_x

#### xlink_y

#### xlink_z

### Methods Documentation

#### classify_strand(self: [pylimer_tools_cpp.CrosslinkerUniverse](#pylimer_tools_cpp.CrosslinkerUniverse), strand_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [pylimer_tools_cpp.StrandCategory](pylimer_tools_cpp.StrandCategory.html.md#pylimer_tools_cpp.StrandCategory)

Classify a strand into its topological StrandCategory.

* **Parameters:**
  **strand_idx** – Index of the strand.
* **Returns:**
  StrandCategory of the strand.

#### compute_bond_length_stats(self: [pylimer_tools_cpp.CrosslinkerUniverse](#pylimer_tools_cpp.CrosslinkerUniverse), box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.html.md#pylimer_tools_cpp.Box), n_mc_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0) → [pylimer_tools_cpp.BondLengthStats](pylimer_tools_cpp.BondLengthStats.html.md#pylimer_tools_cpp.BondLengthStats)

Stream over all strands in O(1) strand memory to compute realized bond length statistics.

* **Parameters:**
  * **box** – Simulation box.
  * **n_mc_steps** – Number of MC equilibration steps per strand (default: 0).
* **Returns:**
  BondLengthStats object.

#### compute_end_to_end_stats(self: [pylimer_tools_cpp.CrosslinkerUniverse](#pylimer_tools_cpp.CrosslinkerUniverse), box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.html.md#pylimer_tools_cpp.Box), category: [pylimer_tools_cpp.StrandCategory](pylimer_tools_cpp.StrandCategory.html.md#pylimer_tools_cpp.StrandCategory) = StrandCategory.NETWORK_STRAND, n_mc_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0) → [pylimer_tools_cpp.EndToEndStats](pylimer_tools_cpp.EndToEndStats.html.md#pylimer_tools_cpp.EndToEndStats)

Stream over strands matching category to compute realized end-to-end distance statistics.

* **Parameters:**
  * **box** – Simulation box.
  * **category** – StrandCategory filter (default: NETWORK_STRAND).
  * **n_mc_steps** – Number of MC equilibration steps per strand (default: 0).
* **Returns:**
  EndToEndStats object.

#### compute_strand_bond_counts(self: [pylimer_tools_cpp.CrosslinkerUniverse](#pylimer_tools_cpp.CrosslinkerUniverse), category: [pylimer_tools_cpp.StrandCategory](pylimer_tools_cpp.StrandCategory.html.md#pylimer_tools_cpp.StrandCategory) = StrandCategory.ALL) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[int](https://docs.python.org/3/builtins/functions.html#int)]

Return the number of bonds for each strand matching category.

* **Parameters:**
  **category** – StrandCategory filter (default: ALL).
* **Returns:**
  List of bond counts per matching strand.

#### count_atom_types(self: [pylimer_tools_cpp.CrosslinkerUniverse](#pylimer_tools_cpp.CrosslinkerUniverse)) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[int](https://docs.python.org/3/builtins/functions.html#int), [int](https://docs.python.org/3/builtins/functions.html#int)]

Count atoms per atom type in O(N_xlinks + N_strands).

* **Returns:**
  Dict mapping atom type to count.

#### count_atoms_by_degree(self: [pylimer_tools_cpp.CrosslinkerUniverse](#pylimer_tools_cpp.CrosslinkerUniverse)) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[int](https://docs.python.org/3/builtins/functions.html#int), [int](https://docs.python.org/3/builtins/functions.html#int)]

Count atoms per bond degree in O(N_xlinks + N_strands).

* **Returns:**
  Dict mapping degree to atom count.

#### determine_effective_functionality_per_type(self: [pylimer_tools_cpp.CrosslinkerUniverse](#pylimer_tools_cpp.CrosslinkerUniverse)) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[int](https://docs.python.org/3/builtins/functions.html#int), [float](https://docs.python.org/3/builtins/functions.html#float)]

Compute effective functionality (mean degree) per atom type in O(N_xlinks + N_strands).

* **Returns:**
  Dict mapping atom type to mean degree.
