# Molecule

### *class* pylimer_tools_cpp.Molecule(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), arg0: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box), arg1: igraph_t, arg2: [pylimer_tools_cpp.MoleculeType](pylimer_tools_cpp.MoleculeType.md#pylimer_tools_cpp.MoleculeType), arg3: [collections.abc.Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)])

Bases: `pybind11_object`

An (ideally) connected series of atoms/beads.

Construct a molecule from a graph structure.

* **Parameters:**
  * **box** – The simulation box containing this molecule
  * **graph** – Pointer to the igraph structure representing connectivity
  * **type** – The type of molecule (see MoleculeType enum)
  * **mass_map** – Map of atom types to their masses

### Methods Summary

| [`compute_bond_lengths`](#pylimer_tools_cpp.Molecule.compute_bond_lengths)(self)                                                                 | Compute the length of each bond in the molecule, respecting periodic boundaries.                                                |
|--------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| [`compute_end_to_end_distance`](#pylimer_tools_cpp.Molecule.compute_end_to_end_distance)(self)                                                   | Compute the end-to-end distance ($R_{ee}$) of this molecule.                                                                    |
| [`compute_end_to_end_distance_with_derived_image_flags`](#pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags)(self) | Compute the end-to-end distance ($R_{ee}$) of this molecule, but ignoring the image flags attached to the atoms.                |
| [`compute_end_to_end_vector`](#pylimer_tools_cpp.Molecule.compute_end_to_end_vector)(self)                                                       | Compute the end-to-end vector ($\overrightarrow{R}_{ee}$) of this molecule.                                                     |
| [`compute_end_to_end_vector_with_derived_image_flags`](#pylimer_tools_cpp.Molecule.compute_end_to_end_vector_with_derived_image_flags)(self)     | Compute the end-to-end vector ($\overrightarrow{R}_{ee}$) of this molecule, but ignoring the image flags attached to the atoms. |
| [`compute_radius_of_gyration`](#pylimer_tools_cpp.Molecule.compute_radius_of_gyration)(self)                                                     | Compute the radius of gyration, $R_g^2$ of this molecule.                                                                       |
| [`compute_radius_of_gyration_with_derived_image_flags`](#pylimer_tools_cpp.Molecule.compute_radius_of_gyration_with_derived_image_flags)(self)   | Compute the radius of gyration, $R_g^2$ of this molecule, but ignoring the image flags attached to the atoms.                   |
| [`compute_total_length`](#pylimer_tools_cpp.Molecule.compute_total_length)(self)                                                                 | Compute the sum of the lengths of all bonds.                                                                                    |
| [`compute_total_mass`](#pylimer_tools_cpp.Molecule.compute_total_mass)(self)                                                                     | Compute the total mass of this molecule.                                                                                        |
| [`compute_total_vector`](#pylimer_tools_cpp.Molecule.compute_total_vector)(self[, ...])                                                          | Compute the sum of all bond vectors.                                                                                            |
| [`compute_vector_from_to`](#pylimer_tools_cpp.Molecule.compute_vector_from_to)(self, atom_id_from, ...)                                          | Compute the sum of all bond vectors between two specified atoms.                                                                |
| [`get_atom_by_id`](#pylimer_tools_cpp.Molecule.get_atom_by_id)(self, atom_id)                                                                    | Get an atom by its id.                                                                                                          |
| [`get_atom_by_vertex_idx`](#pylimer_tools_cpp.Molecule.get_atom_by_vertex_idx)(self, vertex_idx)                                                 | Get an atom for a specific vertex.                                                                                              |
| [`get_atom_id_by_vertex_idx`](#pylimer_tools_cpp.Molecule.get_atom_id_by_vertex_idx)(self, vertex_id)                                            | Get the ID of the atom by the vertex ID of the underlying graph.                                                                |
| [`get_atom_types`](#pylimer_tools_cpp.Molecule.get_atom_types)(self)                                                                             | Query all types (each one for each atom) ordered by atom vertex id.                                                             |
| [`get_atoms`](#pylimer_tools_cpp.Molecule.get_atoms)(self)                                                                                       | Return all atom objects enclosed in this molecule, ordered by vertex id.                                                        |
| [`get_atoms_by_degree`](#pylimer_tools_cpp.Molecule.get_atoms_by_degree)(self, degree)                                                           | Get the atoms that have the specified number of bonds.                                                                          |
| [`get_atoms_by_type`](#pylimer_tools_cpp.Molecule.get_atoms_by_type)(self, type)                                                                 | Get the atoms with the specified type.                                                                                          |
| [`get_atoms_connected_to`](#pylimer_tools_cpp.Molecule.get_atoms_connected_to)(self, atom)                                                       | Get the atoms connected to a specified atom.                                                                                    |
| [`get_atoms_connected_to_vertex`](#pylimer_tools_cpp.Molecule.get_atoms_connected_to_vertex)(self, vertex_idx)                                   | Get the atoms connected to a specified vertex id.                                                                               |
| [`get_atoms_lined_up`](#pylimer_tools_cpp.Molecule.get_atoms_lined_up)(self[, crosslinker_type, ...])                                            | Return all atom objects enclosed in this molecule based on the connectivity.                                                    |
| [`get_bonds`](#pylimer_tools_cpp.Molecule.get_bonds)(self)                                                                                       | Get all bonds.                                                                                                                  |
| [`get_edge_ids_from`](#pylimer_tools_cpp.Molecule.get_edge_ids_from)(self, vertex_id)                                                            | Get the edge IDs incident to a specific vertex.                                                                                 |
| [`get_edge_ids_from_to`](#pylimer_tools_cpp.Molecule.get_edge_ids_from_to)(self, vertex_id_from, ...)                                            | Get the edge IDs of the edges between two specific vertices.                                                                    |
| [`get_edges`](#pylimer_tools_cpp.Molecule.get_edges)(self)                                                                                       | Get all bonds.                                                                                                                  |
| [`get_key`](#pylimer_tools_cpp.Molecule.get_key)(self)                                                                                           | Get a unique identifier for this molecule.                                                                                      |
| [`get_nr_of_atoms`](#pylimer_tools_cpp.Molecule.get_nr_of_atoms)(self)                                                                           | Count and return the number of atoms associated with this molecule.                                                             |
| [`get_nr_of_bonds`](#pylimer_tools_cpp.Molecule.get_nr_of_bonds)(self)                                                                           | Count and return the number of bonds associated with this molecule.                                                             |
| [`get_nr_of_edges_from_to`](#pylimer_tools_cpp.Molecule.get_nr_of_edges_from_to)(self, ...[, max_length])                                        | Get the number of edges in the shortest path between two specific vertices.                                                     |
| [`get_strand_ends`](#pylimer_tools_cpp.Molecule.get_strand_ends)(self[, crosslinker_type, ...])                                                  | Get the ends of the given strand (= molecule).                                                                                  |
| [`get_strand_type`](#pylimer_tools_cpp.Molecule.get_strand_type)(self)                                                                           | Get the type of this molecule (see [`MoleculeType`](pylimer_tools_cpp.MoleculeType.md#pylimer_tools_cpp.MoleculeType) enum).    |
| [`get_vertex_idx_by_atom_id`](#pylimer_tools_cpp.Molecule.get_vertex_idx_by_atom_id)(self, atom_id)                                              | Get the vertex ID of the underlying graph for an atom with a specified ID.                                                      |

### Methods Documentation

#### compute_bond_lengths(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Compute the length of each bond in the molecule, respecting periodic boundaries.

* **Returns:**
  A vector of bond lengths

#### compute_end_to_end_distance(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the end-to-end distance ($R_{ee}$) of this molecule.

#### WARNING
Returns 0.0 if the molecule does not have two or more atoms.

* **Returns:**
  The end-to-end distance

#### compute_end_to_end_distance_with_derived_image_flags(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the end-to-end distance ($R_{ee}$) of this molecule,
but ignoring the image flags attached to the atoms.
This only works for Molecules that can be lined up with
[`get_atoms_lined_up()`](#pylimer_tools_cpp.Molecule.get_atoms_lined_up),
as it needs the atoms sorted such that the periodic box can still be respected somewhat.

#### WARNING
Returns 0.0 if the molecule does not have two or more atoms.
Requires bonds to be shorter than half the box length.

* **Returns:**
  The end-to-end distance with derived image flags

#### compute_end_to_end_vector(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the end-to-end vector ($\overrightarrow{R}_{ee}$) of this molecule.

#### WARNING
Returns 0.0 if the molecule does not have two or more atoms.

* **Returns:**
  The end-to-end vector

#### compute_end_to_end_vector_with_derived_image_flags(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the end-to-end vector ($\overrightarrow{R}_{ee}$) of this molecule,
but ignoring the image flags attached to the atoms.
This only works for Molecules that can be lined up with
[`get_atoms_lined_up()`](#pylimer_tools_cpp.Molecule.get_atoms_lined_up),
as it needs the atoms sorted such that the periodic box can still be respected somewhat.

#### WARNING
Returns 0.0 if the molecule does not have two or more atoms.
Requires bonds to be shorter than half the box length.

* **Returns:**
  The end-to-end vector with derived image flags

#### compute_radius_of_gyration(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the radius of gyration, $R_g^2$ of this molecule.

${R_g}^2 = \\frac{1}{M} \sum_i m_i (r_i - r_{cm})^2$,
where $M$ is the total mass of the molecule, $r_{cm}$
are the coordinates of the center of mass of the molecule and the
sum is over all contained atoms.

* **Returns:**
  The radius of gyration squared

#### compute_radius_of_gyration_with_derived_image_flags(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the radius of gyration, $R_g^2$ of this molecule,
but ignoring the image flags attached to the atoms.
This only works for Molecules that can be lined up with
[`get_atoms_lined_up()`](#pylimer_tools_cpp.Molecule.get_atoms_lined_up),
as it needs the atoms sorted such that the periodic box can still be respected somewhat.
In other words, this function computes the radius of gyration
assuming the distance between two lined-up beads
is smaller than half the periodic box in each direction.

See also: [`compute_radius_of_gyration()`](#pylimer_tools_cpp.Molecule.compute_radius_of_gyration).

* **Returns:**
  The radius of gyration squared with derived image flags

#### compute_total_length(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the sum of the lengths of all bonds.
In most cases, this is equal to the contour length.

* **Returns:**
  The total contour length of the molecule

#### compute_total_mass(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the total mass of this molecule.

* **Returns:**
  The total mass of all atoms in this molecule

#### compute_total_vector(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, close_loop: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the sum of all bond vectors.

* **Parameters:**
  * **crosslinker_type** – The type of crosslinker atoms
  * **close_loop** – Whether to close the loop for calculations
* **Returns:**
  The vector sum of all bonds

#### compute_vector_from_to(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), atom_id_from: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), atom_id_to: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, require_order: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the sum of all bond vectors between two specified atoms.

* **Parameters:**
  * **atom_id_from** – ID of the starting atom
  * **atom_id_to** – ID of the ending atom
  * **crosslinker_type** – The type of crosslinker atoms
  * **require_order** – Whether to require specific ordering
* **Returns:**
  The vector sum from start to end atom

#### get_atom_by_id(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), atom_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)

Get an atom by its id.

* **Parameters:**
  **atom_id** – The atom ID to search for
* **Returns:**
  The atom with the specified ID

#### get_atom_by_vertex_idx(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), vertex_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)

Get an atom for a specific vertex.

* **Parameters:**
  **vertex_idx** – The vertex index to query
* **Returns:**
  The atom at the specified vertex

#### get_atom_id_by_vertex_idx(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), vertex_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the ID of the atom by the vertex ID of the underlying graph.

* **Parameters:**
  **vertex_id** – The vertex index in the underlying graph
* **Returns:**
  The atom ID corresponding to the vertex

#### get_atom_types(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Query all types (each one for each atom) ordered by atom vertex id.

* **Returns:**
  A vector of atom types in vertex order

#### get_atoms(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Return all atom objects enclosed in this molecule, ordered by vertex id.

* **Returns:**
  List of atoms in vertex order

#### get_atoms_by_degree(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), degree: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get the atoms that have the specified number of bonds.

* **Parameters:**
  **degree** – The number of bonds (degree/functionality)
* **Returns:**
  List of atoms with the specified degree

#### get_atoms_by_type(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get the atoms with the specified type.

* **Parameters:**
  **type** – The atom type to search for
* **Returns:**
  List of atoms with the specified type

#### get_atoms_connected_to(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), atom: [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get the atoms connected to a specified atom.

Internally uses [`get_atoms_connected_to()`](#pylimer_tools_cpp.Molecule.get_atoms_connected_to).

* **Parameters:**
  **atom** – The atom to query connections for
* **Returns:**
  List of connected atoms

#### get_atoms_connected_to_vertex(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), vertex_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get the atoms connected to a specified vertex id.

* **Parameters:**
  **vertex_idx** – The vertex index to query
* **Returns:**
  List of connected atoms

#### get_atoms_lined_up(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, assumed_coordinates: [bool](https://docs.python.org/3/library/functions.html#bool) = False, close_loop: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Return all atom objects enclosed in this molecule based on the connectivity.

This method works only for lone chains, atoms and loops,
as it throws an error if the molecule does not allow such a “line-up”,
for example because of crosslinkers.

Use the crosslinker_type parameter to force the atoms in a primary loop
to start with the crosslink.

* **Parameters:**
  * **crosslinker_type** – The type of crosslinker atoms
  * **assumed_coordinates** – Whether to assume coordinates are valid
  * **close_loop** – Whether to close loops
* **Returns:**
  List of atoms in connected order

#### get_bonds(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]]

Get all bonds. Returns a dict with three properties: ‘bond_from’, ‘bond_to’ and ‘bond_type’.

* **Returns:**
  Dictionary with bond information

#### get_edge_ids_from(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), vertex_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the edge IDs incident to a specific vertex.

* **Parameters:**
  **vertex_id** – The vertex to query
* **Returns:**
  A vector of edge IDs connected to the vertex

#### get_edge_ids_from_to(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), vertex_id_from: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), vertex_id_to: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the edge IDs of the edges between two specific vertices.

* **Parameters:**
  * **vertex_id_from** – Starting vertex ID
  * **vertex_id_to** – Ending vertex ID
* **Returns:**
  Vector of edge IDs between the specified vertices

#### get_edges(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]]

Get all bonds. Returns a dict with three properties: ‘edge_from’, ‘edge_to’ and ‘edge_type’.
The order is not necessarily related to any structural property.

#### NOTE
The integer values returned refer to the vertex ids, not the atom ids.
Use `get_atom_id_by_idx()` to translate them to atom ids, or
[`get_bonds()`](#pylimer_tools_cpp.Molecule.get_bonds) to have that done for you.

* **Returns:**
  Dictionary with edge information

#### get_key(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [str](https://docs.python.org/3/library/stdtypes.html#str)

Get a unique identifier for this molecule.

* **Returns:**
  A unique string identifier for this molecule

#### get_nr_of_atoms(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [int](https://docs.python.org/3/library/functions.html#int)

Count and return the number of atoms associated with this molecule.

* **Returns:**
  The number of atoms in this molecule

#### get_nr_of_bonds(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [int](https://docs.python.org/3/library/functions.html#int)

Count and return the number of bonds associated with this molecule.

* **Returns:**
  The number of bonds in this molecule

#### get_nr_of_edges_from_to(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), vertex_id_from: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), vertex_id_to: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), max_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of edges in the shortest path between two specific vertices.

If max_length is provided and positive, it will only consider paths up to that length.

* **Parameters:**
  * **vertex_id_from** – Starting vertex ID
  * **vertex_id_to** – Ending vertex ID
  * **max_length** – Maximum path length to consider (-1 for no limit)
* **Returns:**
  Number of edges in the shortest path, or -1 if no path exists

#### get_strand_ends(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, close_loop: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get the ends of the given strand (= molecule).
In case of a primary loop, the crosslink is returned, if there is one.
Use the argument close_loop to decide, whether this should be returned once or twice.

#### NOTE
Currently only works for linear strands.

* **Parameters:**
  * **crosslinker_type** – The type of crosslinker atoms
  * **close_loop** – Whether to return the crosslinker twice for loops
* **Returns:**
  List of end atoms

#### get_strand_type(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule)) → [pylimer_tools_cpp.MoleculeType](pylimer_tools_cpp.MoleculeType.md#pylimer_tools_cpp.MoleculeType)

Get the type of this molecule (see [`MoleculeType`](pylimer_tools_cpp.MoleculeType.md#pylimer_tools_cpp.MoleculeType) enum).

#### NOTE
This type might be unset; currently, only
[`get_chains_with_crosslinker()`](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.get_chains_with_crosslinker) assigns them automatically.

#### get_vertex_idx_by_atom_id(self: [pylimer_tools_cpp.Molecule](#pylimer_tools_cpp.Molecule), atom_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the vertex ID of the underlying graph for an atom with a specified ID.

* **Parameters:**
  **atom_id** – The atom ID to look up
* **Returns:**
  The vertex index corresponding to the atom
