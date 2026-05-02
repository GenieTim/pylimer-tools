# Universe

### *class* pylimer_tools_cpp.Universe(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), Lx: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), Ly: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), Lz: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat))

Bases: `pybind11_object`

Represents a full Polymer Network structure, a collection of molecules.

This is the main class for representing molecular systems, containing
atoms, bonds, angles, and the simulation box.

Instantiate this Universe (Collection of Molecules) providing the box lengths.

* **Parameters:**
  * **Lx** – Box length in x direction
  * **Ly** – Box length in y direction
  * **Lz** – Box length in z direction

### Methods Summary

| [`add_angles`](#pylimer_tools_cpp.Universe.add_angles)(self, angles_from, angles_via, ...)                                   | Add angles to the Universe.                                                                                                                                                                                    |
|------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`add_atoms`](#pylimer_tools_cpp.Universe.add_atoms)(self, ids, types, x, y, z, nx, ny, nz)                                  | Add atoms to the Universe, vertices to the underlying graph.                                                                                                                                                   |
| [`add_bonds`](#pylimer_tools_cpp.Universe.add_bonds)(\*args, \*\*kwargs)                                                     | Overloaded function.                                                                                                                                                                                           |
| [`add_bonds_with_types`](#pylimer_tools_cpp.Universe.add_bonds_with_types)(self, bonds_from, ...)                            | Add bonds to the underlying atoms, edges to the underlying graph.                                                                                                                                              |
| [`add_dihedral_angles`](#pylimer_tools_cpp.Universe.add_dihedral_angles)(self, angles_from, ...)                             | Add dihedral angles to the Universe.                                                                                                                                                                           |
| [`compute_angles`](#pylimer_tools_cpp.Universe.compute_angles)(self)                                                         | Compute the angle of each angle in the molecule, respecting periodic boundaries.                                                                                                                               |
| [`compute_bond_lengths`](#pylimer_tools_cpp.Universe.compute_bond_lengths)(self)                                             | Compute the length of each bond in the molecule, respecting periodic boundaries.                                                                                                                               |
| [`compute_bond_vectors`](#pylimer_tools_cpp.Universe.compute_bond_vectors)(self)                                             | Compute the vectors of each bond in the molecule, respecting periodic boundaries.                                                                                                                              |
| [`compute_dxs`](#pylimer_tools_cpp.Universe.compute_dxs)(self, atom_ids_to, atom_ids_from)                                   | Compute the dx distance for certain bonds (length in x direction).                                                                                                                                             |
| [`compute_dys`](#pylimer_tools_cpp.Universe.compute_dys)(self, atom_ids_to, atom_ids_from)                                   | Compute the dy distance for certain bonds (length in y direction).                                                                                                                                             |
| [`compute_dzs`](#pylimer_tools_cpp.Universe.compute_dzs)(self, atom_ids_to, atom_ids_from)                                   | Compute the dz distance for certain bonds (length in z direction).                                                                                                                                             |
| [`compute_end_to_end_distances`](#pylimer_tools_cpp.Universe.compute_end_to_end_distances)(self, ...[, ...])                 | Compute the end-to-end distance of each strand in the network.                                                                                                                                                 |
| [`compute_mean_end_to_end_distance`](#pylimer_tools_cpp.Universe.compute_mean_end_to_end_distance)(self, ...)                | Compute the mean of the end-to-end distances of each strand in the network.                                                                                                                                    |
| [`compute_mean_squared_end_to_end_distance`](#pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance)(...)      | Compute the mean square of the end-to-end distances of each strand (incl.                                                                                                                                      |
| [`compute_mean_strand_length`](#pylimer_tools_cpp.Universe.compute_mean_strand_length)(self, ...)                            | Compute the mean number of beads per strand.                                                                                                                                                                   |
| [`compute_number_average_molecular_weight`](#pylimer_tools_cpp.Universe.compute_number_average_molecular_weight)(...)        | Compute the number average molecular weight.                                                                                                                                                                   |
| [`compute_polydispersity_index`](#pylimer_tools_cpp.Universe.compute_polydispersity_index)(self, ...)                        | Compute the polydispersity index: the weight average molecular weight over the number average molecular weight.                                                                                                |
| [`compute_temperature`](#pylimer_tools_cpp.Universe.compute_temperature)(self[, dimensions, k_b])                            | Use the velocities per atom to compute the temperature from the kinetic energy of the system.                                                                                                                  |
| [`compute_total_mass`](#pylimer_tools_cpp.Universe.compute_total_mass)(self)                                                 | Compute the total mass of this network/universe in whatever mass unit was used when [`set_masses()`](#pylimer_tools_cpp.Universe.set_masses) was called.                                                       |
| [`compute_weight_average_molecular_weight`](#pylimer_tools_cpp.Universe.compute_weight_average_molecular_weight)(...)        | Compute the weight average molecular weight.                                                                                                                                                                   |
| [`compute_weight_fractions`](#pylimer_tools_cpp.Universe.compute_weight_fractions)(self)                                     | Compute the weight fractions of each atom type in the network.                                                                                                                                                 |
| [`contract_vertices_along_bond_type`](#pylimer_tools_cpp.Universe.contract_vertices_along_bond_type)(self, ...)              | Merge vertices along a specific bond type.                                                                                                                                                                     |
| [`count_atom_types`](#pylimer_tools_cpp.Universe.count_atom_types)(self)                                                     | Count how often each atom type is present.                                                                                                                                                                     |
| [`count_atoms_in_skin_distance`](#pylimer_tools_cpp.Universe.count_atoms_in_skin_distance)(self, distances)                  | This is a function that may help you to compute the radial distribution function.                                                                                                                              |
| [`count_loop_lengths`](#pylimer_tools_cpp.Universe.count_loop_lengths)(self[, max_length])                                   | Find all loops (below a specific length) and count the number of atoms involved in them.                                                                                                                       |
| [`detect_angles`](#pylimer_tools_cpp.Universe.detect_angles)(self)                                                           | Detect angles in the network based on the current bonds.                                                                                                                                                       |
| [`detect_dihedral_angles`](#pylimer_tools_cpp.Universe.detect_dihedral_angles)(self)                                         | Detect dihedral angles in the network based on the current bonds.                                                                                                                                              |
| [`determine_effective_functionality_per_type`](#pylimer_tools_cpp.Universe.determine_effective_functionality_per_type)(self) | Find the average functionality of each atom type in the network.                                                                                                                                               |
| [`determine_functionality_per_type`](#pylimer_tools_cpp.Universe.determine_functionality_per_type)(self)                     | Find the maximum functionality of each atom type in the network.                                                                                                                                               |
| [`find_loops`](#pylimer_tools_cpp.Universe.find_loops)(self, crosslinker_type[, ...])                                        | Decompose the Universe into loops.                                                                                                                                                                             |
| [`find_minimal_order_loop_from`](#pylimer_tools_cpp.Universe.find_minimal_order_loop_from)(self, ...[, ...])                 | Find the loops in the network starting with one connection.                                                                                                                                                    |
| [`get_angles`](#pylimer_tools_cpp.Universe.get_angles)(self)                                                                 | Get all angles added to this network.                                                                                                                                                                          |
| [`get_atom`](#pylimer_tools_cpp.Universe.get_atom)(self, atom_id)                                                            | Find an atom by its ID.                                                                                                                                                                                        |
| [`get_atom_by_vertex_id`](#pylimer_tools_cpp.Universe.get_atom_by_vertex_id)(self, vertex_id)                                | Find an atom by the ID of the vertex of the underlying graph.                                                                                                                                                  |
| [`get_atom_id_by_vertex_idx`](#pylimer_tools_cpp.Universe.get_atom_id_by_vertex_idx)(self, vertex_id)                        | Get the ID of the atom by the vertex ID of the underlying graph.                                                                                                                                               |
| [`get_atom_types`](#pylimer_tools_cpp.Universe.get_atom_types)(self)                                                         | Get all types (each one for each atom) ordered by atom vertex ID.                                                                                                                                              |
| [`get_atoms`](#pylimer_tools_cpp.Universe.get_atoms)(self)                                                                   | Get all atoms in the universe.                                                                                                                                                                                 |
| [`get_atoms_by_degree`](#pylimer_tools_cpp.Universe.get_atoms_by_degree)(self, functionality)                                | Get the atoms that have the specified number of bonds.                                                                                                                                                         |
| [`get_atoms_by_type`](#pylimer_tools_cpp.Universe.get_atoms_by_type)(self, atom_type)                                        | Query all atoms by their type.                                                                                                                                                                                 |
| [`get_atoms_connected_to`](#pylimer_tools_cpp.Universe.get_atoms_connected_to)(self, atom)                                   | Get the atoms connected to a specified atom.                                                                                                                                                                   |
| [`get_atoms_connected_to_vertex`](#pylimer_tools_cpp.Universe.get_atoms_connected_to_vertex)(self, vertex_idx)               | Get the atoms connected to a specified vertex ID.                                                                                                                                                              |
| [`get_bonds`](#pylimer_tools_cpp.Universe.get_bonds)(self)                                                                   | Get all bonds.                                                                                                                                                                                                 |
| [`get_box`](#pylimer_tools_cpp.Universe.get_box)(self)                                                                       | Get the underlying bounding box object.                                                                                                                                                                        |
| [`get_chains_with_crosslinker`](#pylimer_tools_cpp.Universe.get_chains_with_crosslinker)(self, ...)                          | Decompose the Universe into strands (molecules, which could be either chains, or even lonely atoms), without omitting the crosslinkers (as in [`get_molecules()`](#pylimer_tools_cpp.Universe.get_molecules)). |
| [`get_clusters`](#pylimer_tools_cpp.Universe.get_clusters)(self)                                                             | Get the components of the universe that are not connected to each other.                                                                                                                                       |
| [`get_edge_ids_from`](#pylimer_tools_cpp.Universe.get_edge_ids_from)(self, vertex_id)                                        | Get the edge IDs incident to a specific vertex.                                                                                                                                                                |
| [`get_edge_ids_from_to`](#pylimer_tools_cpp.Universe.get_edge_ids_from_to)(self, vertex_id_from, ...)                        | Get the edge IDs of all edges between two specific vertices.                                                                                                                                                   |
| [`get_edges`](#pylimer_tools_cpp.Universe.get_edges)(self)                                                                   | Get all edges.                                                                                                                                                                                                 |
| [`get_masses`](#pylimer_tools_cpp.Universe.get_masses)(self)                                                                 | Get the mass of one atom per type.                                                                                                                                                                             |
| [`get_molecules`](#pylimer_tools_cpp.Universe.get_molecules)(self, atom_type_to_omit)                                        | Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.                                                                                                           |
| [`get_network_of_crosslinker`](#pylimer_tools_cpp.Universe.get_network_of_crosslinker)(self, ...)                            | Reduce the network to contain only crosslinkers, replacing all the strands with a single bond.                                                                                                                 |
| [`get_nr_of_angles`](#pylimer_tools_cpp.Universe.get_nr_of_angles)(self)                                                     | Query the number of angles that have been added to this universe.                                                                                                                                              |
| [`get_nr_of_atoms`](#pylimer_tools_cpp.Universe.get_nr_of_atoms)(self)                                                       | Query the number of atoms in this universe.                                                                                                                                                                    |
| [`get_nr_of_bonds`](#pylimer_tools_cpp.Universe.get_nr_of_bonds)(self)                                                       | Query the number of bonds associated with this universe.                                                                                                                                                       |
| [`get_nr_of_bonds_of_atom`](#pylimer_tools_cpp.Universe.get_nr_of_bonds_of_atom)(self, arg0)                                 | Count the number of immediate neighbors of an atom, specified by its ID.                                                                                                                                       |
| [`get_nr_of_bonds_of_vertex`](#pylimer_tools_cpp.Universe.get_nr_of_bonds_of_vertex)(self, arg0)                             | Count the number of immediate neighbors of an atom, specified by its vertex ID.                                                                                                                                |
| [`get_nr_of_dihedral_angles`](#pylimer_tools_cpp.Universe.get_nr_of_dihedral_angles)(self)                                   | Query the number of dihedral angles that have been added to this universe.                                                                                                                                     |
| [`get_nr_of_edges_from_to`](#pylimer_tools_cpp.Universe.get_nr_of_edges_from_to)(self, ...[, max_length])                    | Get the number of edges in the shortest path between two specific vertices.                                                                                                                                    |
| [`get_timestep`](#pylimer_tools_cpp.Universe.get_timestep)(self)                                                             | Query the timestep when this universe was captured.                                                                                                                                                            |
| [`get_vertex_degrees`](#pylimer_tools_cpp.Universe.get_vertex_degrees)(self)                                                 | Get the degree (functionality) of each vertex.                                                                                                                                                                 |
| [`get_vertex_idx_by_atom_id`](#pylimer_tools_cpp.Universe.get_vertex_idx_by_atom_id)(self, atom_id)                          | Get the vertex ID of the underlying graph for an atom with a specified ID.                                                                                                                                     |
| [`get_volume`](#pylimer_tools_cpp.Universe.get_volume)(self)                                                                 | Query the volume of the underlying bounding box.                                                                                                                                                               |
| [`has_atom_with_id`](#pylimer_tools_cpp.Universe.has_atom_with_id)(self, atom_id)                                            | Check whether this universe contains an atom with the specified ID.                                                                                                                                            |
| [`has_infinite_strand`](#pylimer_tools_cpp.Universe.has_infinite_strand)(self, arg0, arg1)                                   | Check whether there is a strand (with crosslinker) in the universe that loops through periodic images without coming back.                                                                                     |
| [`hash_angle_type`](#pylimer_tools_cpp.Universe.hash_angle_type)(self, angle_from, angle_via, ...)                           | Convert the three integers to one long number/hash.                                                                                                                                                            |
| [`hash_dihedral_angle_type`](#pylimer_tools_cpp.Universe.hash_dihedral_angle_type)(self, angle_from, ...)                    | Convert the four integers to one long number/hash.                                                                                                                                                             |
| [`interpolate_edges`](#pylimer_tools_cpp.Universe.interpolate_edges)(self, crosslinker_type, ...)                            | Get more or less edges than currently present, interpolating between junctions.                                                                                                                                |
| [`remove_all_angles`](#pylimer_tools_cpp.Universe.remove_all_angles)(self)                                                   | Remove all angles from the Universe.                                                                                                                                                                           |
| [`remove_all_dihedral_angles`](#pylimer_tools_cpp.Universe.remove_all_dihedral_angles)(self)                                 | Remove all dihedral angles from the Universe.                                                                                                                                                                  |
| [`remove_atoms`](#pylimer_tools_cpp.Universe.remove_atoms)(self, atom_ids)                                                   | Remove atoms and all associated bonds by their atom IDs.                                                                                                                                                       |
| [`remove_bonds`](#pylimer_tools_cpp.Universe.remove_bonds)(self, bonds_from, bonds_to)                                       | Remove bonds by their connected atom IDs.                                                                                                                                                                      |
| [`remove_bonds_by_type`](#pylimer_tools_cpp.Universe.remove_bonds_by_type)(self, bond_type)                                  | Remove bonds with a specific type.                                                                                                                                                                             |
| [`replace_atom`](#pylimer_tools_cpp.Universe.replace_atom)(self, atom_id, replacement_atom)                                  | Replace the properties of an atom with the properties of another given atom.                                                                                                                                   |
| [`replace_atom_type`](#pylimer_tools_cpp.Universe.replace_atom_type)(self, atom_id, new_type)                                | Replace the type of an atom with another type.                                                                                                                                                                 |
| [`resample_velocities`](#pylimer_tools_cpp.Universe.resample_velocities)(self, mean, variance[, ...])                        | Resample velocities for atoms in the universe.                                                                                                                                                                 |
| [`set_box`](#pylimer_tools_cpp.Universe.set_box)(self, box[, rescale_atoms])                                                 | Override the currently assigned box with the one specified.                                                                                                                                                    |
| [`set_box_lengths`](#pylimer_tools_cpp.Universe.set_box_lengths)(self, lx, ly, lz[, ...])                                    | Override the currently assigned box with one with the side lengths specified.                                                                                                                                  |
| [`set_mass`](#pylimer_tools_cpp.Universe.set_mass)(self, atom_type, mass)                                                    | Set the mass for a specific atom type.                                                                                                                                                                         |
| [`set_masses`](#pylimer_tools_cpp.Universe.set_masses)(self, mass_per_type)                                                  | Set the mass per type of atom.                                                                                                                                                                                 |
| [`set_timestep`](#pylimer_tools_cpp.Universe.set_timestep)(self, timestep)                                                   | Set the timestep when this Universe was captured.                                                                                                                                                              |
| [`set_vertex_property`](#pylimer_tools_cpp.Universe.set_vertex_property)(self, vertex_id, ...)                               | Set a specific property for a specific vertex.                                                                                                                                                                 |
| [`simplify`](#pylimer_tools_cpp.Universe.simplify)(self)                                                                     | Remove self links and double bonds.                                                                                                                                                                            |

### Methods Documentation

#### add_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), angles_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], angles_via: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], angles_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], angle_types: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [None](https://docs.python.org/3/library/constants.html#None)

Add angles to the Universe. No relation to the underlying graph,
just a method to preserve read & write capabilities.

* **Parameters:**
  * **angles_from** – Vector of atom IDs for angle start points
  * **angles_via** – Vector of atom IDs for angle middle points
  * **angles_to** – Vector of atom IDs for angle end points
  * **angle_types** – Vector of angle types

#### add_atoms(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), ids: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], types: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], x: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], y: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], z: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], nx: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], ny: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], nz: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [None](https://docs.python.org/3/library/constants.html#None)

Add atoms to the Universe, vertices to the underlying graph.

* **Parameters:**
  * **ids** – Vector of atom IDs
  * **types** – Vector of atom types
  * **x** – Vector of x coordinates
  * **y** – Vector of y coordinates
  * **z** – Vector of z coordinates
  * **nx** – Vector of periodic image flags in x direction
  * **ny** – Vector of periodic image flags in y direction
  * **nz** – Vector of periodic image flags in z direction

#### add_bonds(\*args, \*\*kwargs)

Overloaded function.

1. add_bonds(self: pylimer_tools_cpp.Universe, bonds_from: collections.abc.Sequence[typing.SupportsInt], bonds_to: collections.abc.Sequence[typing.SupportsInt]) -> None
   > Add bonds to the underlying atoms, edges to the underlying graph.
   > If the connected atoms are not found, the bonds are silently skipped.
   > * **param bonds_from:**
   >   Vector of atom IDs for bond start points
   > * **param bonds_to:**
   >   Vector of atom IDs for bond end points
2. add_bonds(self: pylimer_tools_cpp.Universe, nr_of_bonds: typing.SupportsInt, bonds_from: collections.abc.Sequence[typing.SupportsInt], bonds_to: collections.abc.Sequence[typing.SupportsInt], bond_types: collections.abc.Sequence[typing.SupportsInt], ignore_non_existent_atoms: bool = False, simplify_universe: bool = True) -> None
   > Add bonds to the underlying atoms, edges to the underlying graph.
   > * **param nr_of_bonds:**
   >   Number of bonds to add
   > * **param bonds_from:**
   >   Vector of atom IDs for bond start points
   > * **param bonds_to:**
   >   Vector of atom IDs for bond end points
   > * **param bond_types:**
   >   Vector of bond types
   > * **param ignore_non_existent_atoms:**
   >   Whether to skip bonds to non-existent atoms
   > * **param simplify_universe:**
   >   Whether to simplify the universe after adding bonds

#### add_bonds_with_types(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), bonds_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], bonds_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], bond_types: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [None](https://docs.python.org/3/library/constants.html#None)

Add bonds to the underlying atoms, edges to the underlying graph.
If the connected atoms are not found, the bonds are silently skipped.

* **Parameters:**
  * **bonds_from** – Vector of atom IDs for bond start points
  * **bonds_to** – Vector of atom IDs for bond end points
  * **bond_types** – Vector of bond types

#### add_dihedral_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), angles_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], angles_via1: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], angles_via2: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], angles_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], angle_types: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [None](https://docs.python.org/3/library/constants.html#None)

Add dihedral angles to the Universe. No relation to the underlying graph,
just a method to preserve read & write capabilities.

* **Parameters:**
  * **angles_from** – Vector of atom IDs for dihedral start points
  * **angles_via1** – Vector of atom IDs for first middle points
  * **angles_via2** – Vector of atom IDs for second middle points
  * **angles_to** – Vector of atom IDs for dihedral end points
  * **angle_types** – Vector of dihedral angle types

#### compute_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Compute the angle of each angle in the molecule, respecting periodic boundaries.

* **Returns:**
  A list of angle values in radians

#### compute_bond_lengths(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Compute the length of each bond in the molecule, respecting periodic boundaries.

* **Returns:**
  A list of bond lengths

#### compute_bond_vectors(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']]

Compute the vectors of each bond in the molecule, respecting periodic boundaries.

* **Returns:**
  A list of bond vectors

#### compute_dxs(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_ids_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], atom_ids_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Compute the dx distance for certain bonds (length in x direction).

* **Parameters:**
  * **atom_ids_to** – Vector of destination atom IDs
  * **atom_ids_from** – Vector of source atom IDs
* **Returns:**
  Vector of x-direction distances

#### compute_dys(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_ids_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], atom_ids_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Compute the dy distance for certain bonds (length in y direction).

* **Parameters:**
  * **atom_ids_to** – Vector of destination atom IDs
  * **atom_ids_from** – Vector of source atom IDs
* **Returns:**
  Vector of y-direction distances

#### compute_dzs(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_ids_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], atom_ids_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Compute the dz distance for certain bonds (length in z direction).

* **Parameters:**
  * **atom_ids_to** – Vector of destination atom IDs
  * **atom_ids_from** – Vector of source atom IDs
* **Returns:**
  Vector of z-direction distances

#### compute_end_to_end_distances(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), derive_image_flags: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Compute the end-to-end distance of each strand in the network.

#### NOTE
Internally, this uses either [`compute_end_to_end_distance()`](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule.compute_end_to_end_distance)
or [`compute_end_to_end_distance_with_derived_image_flags()`](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags),
depending on derive_image_flags.
Invalid strands (where said function returns 0.0 or -1.0) are ignored.

* **Parameters:**
  * **crosslinker_type** – The type of crosslinker atoms
  * **derive_image_flags** – Whether to derive image flags from connectivity
* **Returns:**
  List of end-to-end distances

#### compute_mean_end_to_end_distance(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), derive_image_flags: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the mean of the end-to-end distances of each strand in the network.

#### NOTE
Internally, this uses either [`compute_end_to_end_distance()`](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule.compute_end_to_end_distance)
or [`compute_end_to_end_distance_with_derived_image_flags()`](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags),
depending on derive_image_flags.
Invalid strands (where said function returns 0.0 or -1.0) are ignored.

* **Parameters:**
  * **crosslinker_type** – The type of crosslinker atoms
  * **derive_image_flags** – Whether to derive image flags from connectivity
* **Returns:**
  Mean end-to-end distance

#### compute_mean_squared_end_to_end_distance(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), only_those_with_two_crosslinkers: [bool](https://docs.python.org/3/library/functions.html#bool) = False, derive_image_flags: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the mean square of the end-to-end distances of each strand (incl. crosslinks) in the network.

#### NOTE
Internally, this uses either [`compute_end_to_end_distance()`](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule.compute_end_to_end_distance)
or [`compute_end_to_end_distance_with_derived_image_flags()`](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags),
depending on derive_image_flags.
Invalid strands (where said function returns 0.0 or -1.0) are ignored.

* **Parameters:**
  * **crosslinker_type** – The type of crosslinker atoms
  * **only_those_with_two_crosslinkers** – Whether to only consider strands with two crosslinkers
  * **derive_image_flags** – Whether to derive image flags from connectivity
* **Returns:**
  Mean squared end-to-end distance

#### compute_mean_strand_length(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the mean number of beads per strand.

* **Parameters:**
  **crosslinker_type** – The type of crosslinker atoms
* **Returns:**
  Mean strand length

#### compute_number_average_molecular_weight(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the number average molecular weight.

#### NOTE
Crosslinkers are ignored completely.

* **Parameters:**
  **crosslinker_type** – The type of crosslinker atoms
* **Returns:**
  Number average molecular weight

#### compute_polydispersity_index(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the polydispersity index:
the weight average molecular weight over the number average molecular weight.

#### compute_temperature(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), dimensions: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 3, k_b: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [float](https://docs.python.org/3/library/functions.html#float)

Use the velocities per atom to compute the temperature from the kinetic energy of the system.

* **Parameters:**
  * **dimensions** – Number of dimensions (typically 3)
  * **k_b** – Boltzmann constant value
* **Returns:**
  Computed temperature

#### compute_total_mass(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the total mass of this network/universe in whatever mass unit was used when
[`set_masses()`](#pylimer_tools_cpp.Universe.set_masses) was called.

* **Returns:**
  Total mass of the universe

#### compute_weight_average_molecular_weight(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight average molecular weight.

#### NOTE
Crosslinkers are ignored completely.

* **Parameters:**
  **crosslinker_type** – The type of crosslinker atoms
* **Returns:**
  Weight average molecular weight

#### compute_weight_fractions(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float)]

Compute the weight fractions of each atom type in the network.

If no masses are stored, assumes a mass of 1 for each atom.

If the total mass is 0., returns the total mass per atom type.

#### contract_vertices_along_bond_type(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), bond_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)

Merge vertices along a specific bond type.

May result in new self-loops; use [`simplify()`](#pylimer_tools_cpp.Universe.simplify) to remove them.

* **Parameters:**
  **bond_type** – The bond type to contract along

#### count_atom_types(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [int](https://docs.python.org/3/library/functions.html#int)]

Count how often each atom type is present.

* **Returns:**
  Dictionary mapping atom types to their counts

#### count_atoms_in_skin_distance(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), distances: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], unwrapped: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

This is a function that may help you to compute the radial distribution function.
It loops through all atoms and counts neighbors within specified distance ranges.

* **Parameters:**
  * **distances** – The edges of the bins for distance counting
  * **unwrapped** – Whether to measure the distance in unwrapped coordinates or as PBC-corrected distance
* **Returns:**
  Array of counts for each distance bin

#### count_loop_lengths(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), max_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [int](https://docs.python.org/3/library/functions.html#int)]

Find all loops (below a specific length) and count the number of atoms involved in them.
Returns the count, how many loops per length are found.

#### detect_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]]

Detect angles in the network based on the current bonds.
Return the result in the same format as [`get_angles()`](#pylimer_tools_cpp.Universe.get_angles),
but all angles that are detected in the network, rather than the ones already set.
Note that the angle types are determined by
[`hash_angle_type()`](#pylimer_tools_cpp.Universe.hash_angle_type),
which serves angle types that should be mapped by you back to smaller numbers,
before serving them again to [`add_angles()`](#pylimer_tools_cpp.Universe.add_angles),
if you want to have them written e.g. for LAMMPS.

* **Returns:**
  Dictionary with detected angle information

#### detect_dihedral_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]]

Detect dihedral angles in the network based on the current bonds.
Return the result in the same format as `get_dihedral_angles()`,
but all dihedral angles that are detected in the network, rather than the ones already set.
Note that the angle types are determined by
[`hash_dihedral_angle_type()`](#pylimer_tools_cpp.Universe.hash_dihedral_angle_type),
which serves angle types that should be mapped by you back to smaller numbers,
before serving them to [`add_dihedral_angles()`](#pylimer_tools_cpp.Universe.add_dihedral_angles),
if you want to have them written e.g. for LAMMPS.

* **Returns:**
  Dictionary with detected dihedral angle information

#### determine_effective_functionality_per_type(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float)]

Find the average functionality of each atom type in the network.

* **Returns:**
  Dictionary mapping atom types to average functionality

#### determine_functionality_per_type(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [int](https://docs.python.org/3/library/functions.html#int)]

Find the maximum functionality of each atom type in the network.

* **Returns:**
  Dictionary mapping atom types to maximum functionality

#### find_loops(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), max_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1, skip_self_loops: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [list](https://docs.python.org/3/library/stdtypes.html#list)[[list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]]]

Decompose the Universe into loops.
The primary index specifies the degree of the loop.

#### find_minimal_order_loop_from(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), loop_start: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), loop_step1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), max_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1, skip_self_loops: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Find the loops in the network starting with one connection.

#### WARNING
There are exponentially many paths between two crosslinkers of a network,
and you may run out of memory when using this function, if your Universe/Network is lattice-like.
You can use the max_length parameter to restrict the algorithm to only search for loops up to a certain length.
Use a negative value to find all loops and paths.

* **Parameters:**
  * **loop_start** – The atom ID to start the search from
  * **loop_step1** – The first step to take, i.e., the first bond to follow
  * **max_length** – The maximum length of the loop to find, or -1 for no limit
  * **skip_self_loops** – Whether to skip self-loops (i.e., loops that start and end at the same atom; only relevant if loop_start is equal to loop_step1)
* **Returns:**
  List of loops found

#### get_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]]

Get all angles added to this network.

Returns a dict with three properties: ‘angle_from’, ‘angle_via’ and ‘angle_to’.

#### NOTE
The integer values returned refer to the atom IDs, not the vertex IDs.
Use `get_idx_by_atom_id()` to translate them to vertex IDs.

* **Returns:**
  Dictionary with angle information

#### get_atom(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)

Find an atom by its ID.

* **Parameters:**
  **atom_id** – The ID of the atom to find
* **Returns:**
  The atom with the specified ID

#### get_atom_by_vertex_id(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), vertex_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)

Find an atom by the ID of the vertex of the underlying graph.

* **Parameters:**
  **vertex_id** – The vertex ID to query
* **Returns:**
  The atom at the specified vertex

#### get_atom_id_by_vertex_idx(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), vertex_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the ID of the atom by the vertex ID of the underlying graph.

* **Parameters:**
  **vertex_id** – The vertex index in the underlying graph
* **Returns:**
  The atom ID corresponding to the vertex

#### get_atom_types(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get all types (each one for each atom) ordered by atom vertex ID.

* **Returns:**
  Vector of atom types in vertex order

#### get_atoms(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get all atoms in the universe.

* **Returns:**
  List of all atoms

#### get_atoms_by_degree(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), functionality: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get the atoms that have the specified number of bonds.

#### get_atoms_by_type(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Query all atoms by their type.

#### get_atoms_connected_to(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom: [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get the atoms connected to a specified atom.

Internally uses [`get_atoms_connected_to()`](#pylimer_tools_cpp.Universe.get_atoms_connected_to).

* **Parameters:**
  **atom** – The atom to query connections for
* **Returns:**
  List of connected atoms

#### get_atoms_connected_to_vertex(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), vertex_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

Get the atoms connected to a specified vertex ID.

* **Parameters:**
  **vertex_idx** – The vertex index to query
* **Returns:**
  List of connected atoms

#### get_bonds(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]]

Get all bonds. Returns a dict with three properties: ‘bond_from’, ‘bond_to’ and ‘bond_type’.
The order is not necessarily related to any structural characteristic.

* **Returns:**
  Dictionary with bond information

#### get_box(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)

Get the underlying bounding box object.

* **Returns:**
  The simulation box

#### get_chains_with_crosslinker(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Molecule](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule)]

Decompose the Universe into strands (molecules, which could be either chains, or even lonely atoms), without omitting the crosslinkers
(as in [`get_molecules()`](#pylimer_tools_cpp.Universe.get_molecules)).
In turn, e.g. for a tetrafunctional crosslinker, it will be 4 times in the resulting molecules.

#### NOTE
Crosslinkers without bonds to non-crosslinkers are not returned
(i.e., single crosslinkers, are not counted as strands).

* **Parameters:**
  **crosslinker_type** – The type of crosslinker atoms
* **Returns:**
  List of chains including crosslinkers

#### get_clusters(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)]

Get the components of the universe that are not connected to each other.
Returns a list of [`Universe`](#pylimer_tools_cpp.Universe) objects.
Unconnected, free atoms/beads become their own [`Universe`](#pylimer_tools_cpp.Universe).

* **Returns:**
  List of disconnected Universe components

#### get_edge_ids_from(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), vertex_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the edge IDs incident to a specific vertex.

* **Parameters:**
  **vertex_id** – The ID of the vertex to query
* **Returns:**
  List of edge IDs connected to the vertex

#### get_edge_ids_from_to(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), vertex_id_from: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), vertex_id_to: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the edge IDs of all edges between two specific vertices.

* **Parameters:**
  * **vertex_id_from** – The starting vertex ID
  * **vertex_id_to** – The ending vertex ID
* **Returns:**
  List of edge IDs connecting the two vertices

#### get_edges(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]]

Get all edges. Returns a dict with three properties: ‘edge_from’, ‘edge_to’ and ‘edge_type’.
The order is not necessarily related to any structural characteristic.

#### NOTE
The integer values returned refer to the vertex IDs, not the atom IDs.
Use `get_atom_id_by_idx()` to translate them to atom IDs, or
[`get_bonds()`](#pylimer_tools_cpp.Universe.get_bonds) to have that done for you.

* **Returns:**
  Dictionary with edge information

#### get_masses(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float)]

Get the mass of one atom per type.

* **Returns:**
  Dictionary mapping atom types to masses

#### get_molecules(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_type_to_omit: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Molecule](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule)]

Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.

Reduces the Universe to a list of molecules.
Specify the crosslinker_type to an existing type ID,
then those atoms will be omitted, and this function returns chains instead.

* **Parameters:**
  **atom_type_to_omit** – The type of atom to omit from the universe to end up with the desired molecules (e.g., the type of the crosslinkers).
* **Returns:**
  List of molecules

#### get_network_of_crosslinker(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)

Reduce the network to contain only crosslinkers, replacing all the strands with a single bond.
Useful e.g. to reduce the memory usage and runtime of
[`find_loops()`](#pylimer_tools_cpp.Universe.find_loops) or
[`has_infinite_strand()`](#pylimer_tools_cpp.Universe.has_infinite_strand).

Further use [`simplify()`](#pylimer_tools_cpp.Universe.simplify) to remove primary loops.

* **Parameters:**
  **crosslinker_type** – The type of crosslinker atoms
* **Returns:**
  Reduced network containing only crosslinkers

#### get_nr_of_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [int](https://docs.python.org/3/library/functions.html#int)

Query the number of angles that have been added to this universe.

#### get_nr_of_atoms(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [int](https://docs.python.org/3/library/functions.html#int)

Query the number of atoms in this universe.

* **Returns:**
  Number of atoms

#### get_nr_of_bonds(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [int](https://docs.python.org/3/library/functions.html#int)

Query the number of bonds associated with this universe.

#### get_nr_of_bonds_of_atom(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [int](https://docs.python.org/3/library/functions.html#int)

Count the number of immediate neighbors of an atom, specified by its ID.

* **Parameters:**
  **atom_id** – The ID of the atom to query
* **Returns:**
  Number of bonds connected to the atom

#### get_nr_of_bonds_of_vertex(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [int](https://docs.python.org/3/library/functions.html#int)

Count the number of immediate neighbors of an atom, specified by its vertex ID.

* **Parameters:**
  **vertex_id** – The vertex ID of the atom to query
* **Returns:**
  Number of bonds connected to the vertex

#### get_nr_of_dihedral_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [int](https://docs.python.org/3/library/functions.html#int)

Query the number of dihedral angles that have been added to this universe.

#### get_nr_of_edges_from_to(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), vertex_id_from: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), vertex_id_to: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), max_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of edges in the shortest path between two specific vertices.

If max_length is provided and positive, it will only consider paths up to that length.

* **Parameters:**
  * **vertex_id_from** – The starting vertex ID
  * **vertex_id_to** – The ending vertex ID
  * **max_length** – Maximum path length to consider (default: -1 for no limit)
* **Returns:**
  Number of edges in shortest path, or -1 if no path exists

#### get_timestep(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [int](https://docs.python.org/3/library/functions.html#int)

Query the timestep when this universe was captured.

#### get_vertex_degrees(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the degree (functionality) of each vertex.

#### get_vertex_idx_by_atom_id(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the vertex ID of the underlying graph for an atom with a specified ID.

* **Parameters:**
  **atom_id** – The atom ID to look up
* **Returns:**
  The vertex index corresponding to the atom

#### get_volume(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [float](https://docs.python.org/3/library/functions.html#float)

Query the volume of the underlying bounding box.

* **Returns:**
  The volume of the simulation box

#### has_atom_with_id(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [bool](https://docs.python.org/3/library/functions.html#bool)

Check whether this universe contains an atom with the specified ID.

* **Parameters:**
  **atom_id** – The atom ID to check
* **Returns:**
  True if the atom exists in this universe, False otherwise

#### has_infinite_strand(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), arg1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [bool](https://docs.python.org/3/library/functions.html#bool)

Check whether there is a strand (with crosslinker) in the universe that loops through periodic images without coming back.

> #### WARNING
> There are exponentially many paths between two crosslinkers of a network,
> and you may run out of memory when using this function, if your Universe/Network is lattice-like.
* **Returns:**
  True if infinite strands are detected, False otherwise

#### hash_angle_type(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), angle_from: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), angle_via: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), angle_to: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [int](https://docs.python.org/3/library/functions.html#int)

Convert the three integers to one long number/hash.
Used internally for duplicate detection.

* **Parameters:**
  * **angle_from** – First atom ID in the angle
  * **angle_via** – Middle atom ID in the angle
  * **angle_to** – Last atom ID in the angle
* **Returns:**
  Hash value for the angle type

#### hash_dihedral_angle_type(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), angle_from: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), angle_via1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), angle_via2: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), angle_to: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [int](https://docs.python.org/3/library/functions.html#int)

Convert the four integers to one long number/hash.
Used internally for duplicate detection.

* **Parameters:**
  * **angle_from** – First atom ID in the dihedral
  * **angle_via1** – Second atom ID in the dihedral
  * **angle_via2** – Third atom ID in the dihedral
  * **angle_to** – Fourth atom ID in the dihedral
* **Returns:**
  Hash value for the dihedral angle type

#### interpolate_edges(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), interpolation_factor: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[int](https://docs.python.org/3/library/functions.html#int), [int](https://docs.python.org/3/library/functions.html#int)]]

Get more or less edges than currently present, interpolating between junctions.

* **Parameters:**
  * **crosslinker_type** – The type of crosslinker atoms
  * **interpolation_factor** – Factor for edge interpolation
* **Returns:**
  Interpolated edge structure

#### remove_all_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [None](https://docs.python.org/3/library/constants.html#None)

Remove all angles from the Universe.
This will not remove the atoms or bonds, just the angles.

#### WARNING
This will not remove dihedral angles, use [`remove_all_dihedral_angles()`](#pylimer_tools_cpp.Universe.remove_all_dihedral_angles) for that.

#### remove_all_dihedral_angles(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [None](https://docs.python.org/3/library/constants.html#None)

Remove all dihedral angles from the Universe.
This will not remove the atoms or bonds, just the stored dihedral angles.

#### remove_atoms(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_ids: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [None](https://docs.python.org/3/library/constants.html#None)

Remove atoms and all associated bonds by their atom IDs.

* **Parameters:**
  **atom_ids** – Vector of atom IDs to remove

#### remove_bonds(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), bonds_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], bonds_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [None](https://docs.python.org/3/library/constants.html#None)

Remove bonds by their connected atom IDs.

* **Parameters:**
  * **bonds_from** – Vector of starting atom IDs
  * **bonds_to** – Vector of ending atom IDs

#### remove_bonds_by_type(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), bond_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [None](https://docs.python.org/3/library/constants.html#None)

Remove bonds with a specific type.

* **Parameters:**
  **bond_type** – The bond type to remove

#### replace_atom(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), replacement_atom: [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)) → [None](https://docs.python.org/3/library/constants.html#None)

Replace the properties of an atom with the properties of another given atom.

* **Parameters:**
  * **atom_id** – ID of the atom to replace
  * **replacement_atom** – The atom with new properties

#### replace_atom_type(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), new_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [None](https://docs.python.org/3/library/constants.html#None)

Replace the type of an atom with another type.

* **Parameters:**
  * **atom_id** – ID of the atom to modify
  * **new_type** – The new atom type

#### resample_velocities(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), mean: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), variance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), seed: [str](https://docs.python.org/3/library/stdtypes.html#str) = '', is_2d: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Resample velocities for atoms in the universe.

* **Parameters:**
  * **mean** – Mean velocity
  * **variance** – Velocity variance
  * **seed** – Random seed for velocity generation
  * **is_2d** – Whether to omit sampling the z direction

#### set_box(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box), rescale_atoms: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Override the currently assigned box with the one specified.

* **Parameters:**
  * **box** – The new box to assign
  * **rescale_atoms** – Whether to rescale atom positions

#### set_box_lengths(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), lx: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), ly: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), lz: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), rescale_atoms: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Override the currently assigned box with one with the side lengths specified.

* **Parameters:**
  * **lx** – Length in x direction
  * **ly** – Length in y direction
  * **lz** – Length in z direction
  * **rescale_atoms** – Whether to rescale atom positions

#### set_mass(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), mass: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [None](https://docs.python.org/3/library/constants.html#None)

Set the mass for a specific atom type.

* **Parameters:**
  * **atom_type** – The atom type to set mass for
  * **mass** – The mass value to assign

#### set_masses(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), mass_per_type: [collections.abc.Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)]) → [None](https://docs.python.org/3/library/constants.html#None)

Set the mass per type of atom.

* **Parameters:**
  **mass_per_type** – Map of atom types to their masses

#### set_timestep(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), timestep: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [None](https://docs.python.org/3/library/constants.html#None)

Set the timestep when this Universe was captured.

* **Parameters:**
  **timestep** – The timestep value

#### set_vertex_property(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe), vertex_id: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), property_name: [str](https://docs.python.org/3/library/stdtypes.html#str), value: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [None](https://docs.python.org/3/library/constants.html#None)

Set a specific property for a specific vertex.

* **Parameters:**
  * **vertex_id** – The vertex ID to modify
  * **property_name** – Name of the property to set
  * **value** – Value to assign to the property

#### simplify(self: [pylimer_tools_cpp.Universe](#pylimer_tools_cpp.Universe)) → [None](https://docs.python.org/3/library/constants.html#None)

Remove self links and double bonds. This function is called
automatically after adding bonds.

This operation cleans up the graph structure by removing
redundant connections.
