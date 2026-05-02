# MEHPForceBalance

### *class* pylimer_tools_cpp.MEHPForceBalance(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, is_2d: [bool](https://docs.python.org/3/library/functions.html#bool) = False, remove_2functional_crosslinkers: [bool](https://docs.python.org/3/library/functions.html#bool) = False, remove_dangling_chains: [bool](https://docs.python.org/3/library/functions.html#bool) = False)

Bases: `pybind11_object`

A small simulation tool for quickly minimizing the force between the crosslinker beads.

This is the second implementation in the group of MEHP provided by this package.
The distinct feature here is the slip-links: a form of entanglement,
represented as an entanglement link, just like a four-functional crosslink,
but with the ability to slip along the two associated strands, therewith
adjusting the fraction of the contour length on both sides of the link.

Please cite Bernhard and Gusev [[BG25](../acknowledgements.md#id12)] if you use this method.

Instantiate the simulator for a certain universe.

* **Parameters:**
  * **universe** – The universe to simulate with
  * **crosslinker_type** – The atom type of the crosslinkers. Needed to reduce the network.
  * **is2D** – Whether to ignore the z direction.
  * **kappa** – The spring constant
  * **remove_2functionalCrosslinkers** – whether to keep or remove the 2-functional crosslinkers when setting up the network
  * **remove_dangling_chains** – whether to keep or remove obviously dangling chains when setting up the network

### Attributes Summary

| [`network`](#pylimer_tools_cpp.MEHPForceBalance.network)   |    |
|------------------------------------------------------------|----|

### Methods Summary

| [`add_sliplinks`](#pylimer_tools_cpp.MEHPForceBalance.add_sliplinks)(self, strand_idx_1, ...[, ...])                                    | Add new slip-links                                                                                                                                                                                                                                                                        |
|-----------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`add_sliplinks_based_on_cycles`](#pylimer_tools_cpp.MEHPForceBalance.add_sliplinks_based_on_cycles)(self[, ...])                       | Detect and add slip-links based on detected entanglements.                                                                                                                                                                                                                                |
| [`config_assume_box_large_enough`](#pylimer_tools_cpp.MEHPForceBalance.config_assume_box_large_enough)(self[, ...])                     | Configure whether to run PBC on the bonds or not.                                                                                                                                                                                                                                         |
| [`config_entanglement_type`](#pylimer_tools_cpp.MEHPForceBalance.config_entanglement_type)(self[, type])                                | To have certain crosslinks behave as entanglements in the removal process, you can specify here a type, that you have used in the universe to specify: - the type of entanglement atoms (expected with functionality f = 3), - and the entanglement-bonds between the entanglement atoms. |
| [`config_mean_bond_length`](#pylimer_tools_cpp.MEHPForceBalance.config_mean_bond_length)(self[, b])                                     | Configure the $b$ used e.g. for the topological Gamma-factor.                                                                                                                                                                                                                             |
| [`config_simplification_frequency`](#pylimer_tools_cpp.MEHPForceBalance.config_simplification_frequency)(self[, ...])                   | Configure every how many steps to simplify the structure.                                                                                                                                                                                                                                 |
| [`config_spring_breaking_distance`](#pylimer_tools_cpp.MEHPForceBalance.config_spring_breaking_distance)(self[, ...])                   | Configure the "force" (distance over contour length) at which the bonds break.                                                                                                                                                                                                            |
| [`config_spring_constant`](#pylimer_tools_cpp.MEHPForceBalance.config_spring_constant)(self[, kappa])                                   | Configure the spring constant used in the simulation.                                                                                                                                                                                                                                     |
| [`config_step_output`](#pylimer_tools_cpp.MEHPForceBalance.config_step_output)(self, output_configuration)                              | Set which values to log during the simulation.                                                                                                                                                                                                                                            |
| [`construct_with_random_sliplinks`](#pylimer_tools_cpp.MEHPForceBalance.construct_with_random_sliplinks)(universe, ...)                 | Instantiate this simulator with randomly chosen slip-links.                                                                                                                                                                                                                               |
| [`deform_to`](#pylimer_tools_cpp.MEHPForceBalance.deform_to)(self, new_box)                                                             | Perform a deformation of the system box to a different box.                                                                                                                                                                                                                               |
| [`evaluate_spring_distance`](#pylimer_tools_cpp.MEHPForceBalance.evaluate_spring_distance)(self, network, ...)                          | Evaluate the distance vector for a specific spring.                                                                                                                                                                                                                                       |
| [`evaluate_spring_distance_from`](#pylimer_tools_cpp.MEHPForceBalance.evaluate_spring_distance_from)(self, network, ...)                | Evaluate the spring distance from a specific link.                                                                                                                                                                                                                                        |
| [`evaluate_spring_distance_to`](#pylimer_tools_cpp.MEHPForceBalance.evaluate_spring_distance_to)(self, network, ...)                    | Evaluate the spring distance to a specific link.                                                                                                                                                                                                                                          |
| [`get_average_strand_length`](#pylimer_tools_cpp.MEHPForceBalance.get_average_strand_length)(self)                                      | Get the average length of the strands.                                                                                                                                                                                                                                                    |
| [`get_coordinates`](#pylimer_tools_cpp.MEHPForceBalance.get_coordinates)(self)                                                          | Get the current coordinates of the crosslinkers and entanglement links.                                                                                                                                                                                                                   |
| [`get_crosslinker_universe`](#pylimer_tools_cpp.MEHPForceBalance.get_crosslinker_universe)(self)                                        | Returns the universe [of crosslinkers] with the positions of the current state of the simulation.                                                                                                                                                                                         |
| [`get_current_spring_lengths`](#pylimer_tools_cpp.MEHPForceBalance.get_current_spring_lengths)(self)                                    | Get the spring distances.                                                                                                                                                                                                                                                                 |
| [`get_current_spring_vectors`](#pylimer_tools_cpp.MEHPForceBalance.get_current_spring_vectors)(self)                                    | Get the spring vectors.                                                                                                                                                                                                                                                                   |
| [`get_current_strand_vectors`](#pylimer_tools_cpp.MEHPForceBalance.get_current_strand_vectors)(self)                                    |                                                                                                                                                                                                                                                                                           |
| [`get_dangling_weight_fraction`](#pylimer_tools_cpp.MEHPForceBalance.get_dangling_weight_fraction)(self[, tolerance])                   | Compute the weight fraction of non-active strands                                                                                                                                                                                                                                         |
| [`get_default_mean_bond_length`](#pylimer_tools_cpp.MEHPForceBalance.get_default_mean_bond_length)(self)                                | Returns the value effectively used in [`get_gamma_factor()`](#pylimer_tools_cpp.MEHPForceBalance.get_gamma_factor) for $b$ in $\langle R_{0,\eta}^2 = N_{\eta} b^2\rangle$.                                                                                                               |
| [`get_displacement_residual_norm`](#pylimer_tools_cpp.MEHPForceBalance.get_displacement_residual_norm)(self[, ...])                     | Get the current link displacement residual norm.                                                                                                                                                                                                                                          |
| [`get_displacements`](#pylimer_tools_cpp.MEHPForceBalance.get_displacements)(self)                                                      | Get the current link displacements.                                                                                                                                                                                                                                                       |
| [`get_effective_functionality_of_atoms`](#pylimer_tools_cpp.MEHPForceBalance.get_effective_functionality_of_atoms)(self[, ...])         | Returns the number of active strands connected to each atom, atomId used as index                                                                                                                                                                                                         |
| [`get_exit_reason`](#pylimer_tools_cpp.MEHPForceBalance.get_exit_reason)(self)                                                          | Returns the reason for termination of the simulation                                                                                                                                                                                                                                      |
| [`get_force_magnitude_vector`](#pylimer_tools_cpp.MEHPForceBalance.get_force_magnitude_vector)(self, arg0)                              | Evaluate the norm of the force on each (slip- or cross-) link.                                                                                                                                                                                                                            |
| [`get_force_on`](#pylimer_tools_cpp.MEHPForceBalance.get_force_on)(self, link_idx[, ...])                                               | Evaluate the force on a particular (slip- or cross-) link.                                                                                                                                                                                                                                |
| [`get_gamma_factor`](#pylimer_tools_cpp.MEHPForceBalance.get_gamma_factor)(self[, b02, nr_of_chains, ...])                              | Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:                                                                                                                                                                                                                        |
| [`get_gamma_factors`](#pylimer_tools_cpp.MEHPForceBalance.get_gamma_factors)(self, b02[, ...])                                          | Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)                                                                                                                                                                   |
| [`get_gamma_factors_in_dir`](#pylimer_tools_cpp.MEHPForceBalance.get_gamma_factors_in_dir)(self, b02, direction)                        | Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)                                                                                                                                        |
| [`get_ids_of_active_nodes`](#pylimer_tools_cpp.MEHPForceBalance.get_ids_of_active_nodes)(self[, tolerance])                             | Get the atom ids of the nodes that are considered active.                                                                                                                                                                                                                                 |
| [`get_initial_coordinates`](#pylimer_tools_cpp.MEHPForceBalance.get_initial_coordinates)(self)                                          | Get the initial coordinates of the (remaining) crosslinkers and entanglement links.                                                                                                                                                                                                       |
| [`get_neighbour_link_indices`](#pylimer_tools_cpp.MEHPForceBalance.get_neighbour_link_indices)(network, link_idx)                       | Get the indices of neighboring links for a given link.                                                                                                                                                                                                                                    |
| [`get_nr_of_active_nodes`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_active_nodes)(self[, tolerance])                               | Get the number of active nodes (incl.                                                                                                                                                                                                                                                     |
| [`get_nr_of_active_springs`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_active_springs)(self[, tolerance])                           | Get the number of active springs remaining after running the simulation.                                                                                                                                                                                                                  |
| [`get_nr_of_active_strands`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_active_strands)(self[, tolerance])                           | Get the number of active strands remaining after running the simulation.                                                                                                                                                                                                                  |
| [`get_nr_of_active_strands_in_dir`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_active_strands_in_dir)(self, direction)               | Get the number of active strands remaining after running the simulation.                                                                                                                                                                                                                  |
| [`get_nr_of_atoms`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_atoms)(self)                                                          |                                                                                                                                                                                                                                                                                           |
| [`get_nr_of_bonds`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_bonds)(self)                                                          |                                                                                                                                                                                                                                                                                           |
| [`get_nr_of_extra_atoms`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_extra_atoms)(self)                                              |                                                                                                                                                                                                                                                                                           |
| [`get_nr_of_extra_bonds`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_extra_bonds)(self)                                              |                                                                                                                                                                                                                                                                                           |
| [`get_nr_of_intra_chain_sliplinks`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_intra_chain_sliplinks)(self)                          |                                                                                                                                                                                                                                                                                           |
| [`get_nr_of_iterations`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_iterations)(self)                                                | Returns the number of iterations used for force relaxation so far.                                                                                                                                                                                                                        |
| [`get_nr_of_nodes`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_nodes)(self)                                                          | Get the number of nodes (crosslinkers) considered in this simulation.                                                                                                                                                                                                                     |
| [`get_nr_of_strands`](#pylimer_tools_cpp.MEHPForceBalance.get_nr_of_strands)(self)                                                      | Get the number of strands considered in this simulation.                                                                                                                                                                                                                                  |
| [`get_overall_strand_lengths`](#pylimer_tools_cpp.MEHPForceBalance.get_overall_strand_lengths)(self)                                    | Get the sum of the lengths of the springs of each strand.                                                                                                                                                                                                                                 |
| [`get_pressure`](#pylimer_tools_cpp.MEHPForceBalance.get_pressure)(self)                                                                | Returns the pressure at the current state of the simulation.                                                                                                                                                                                                                              |
| [`get_soluble_weight_fraction`](#pylimer_tools_cpp.MEHPForceBalance.get_soluble_weight_fraction)(self[, tolerance])                     | Compute the weight fraction of strands connected to active strands (any depth).                                                                                                                                                                                                           |
| [`get_strand_partition_indices_of_sliplink`](#pylimer_tools_cpp.MEHPForceBalance.get_strand_partition_indices_of_sliplink)(...)         | Get the indices of strand partitions associated with a slip-link.                                                                                                                                                                                                                         |
| [`get_strand_partitions`](#pylimer_tools_cpp.MEHPForceBalance.get_strand_partitions)(self)                                              | Get the current strand partitions (the fraction of the contour length associated with each spring).                                                                                                                                                                                       |
| [`get_stress_on`](#pylimer_tools_cpp.MEHPForceBalance.get_stress_on)(self, link_idx[, ...])                                             | Evaluate the stress on a particular (slip- or cross-) link.                                                                                                                                                                                                                               |
| [`get_stress_tensor`](#pylimer_tools_cpp.MEHPForceBalance.get_stress_tensor)(self[, ...])                                               | Returns the stress tensor at the current state of the simulation.                                                                                                                                                                                                                         |
| [`get_stress_tensor_link_based`](#pylimer_tools_cpp.MEHPForceBalance.get_stress_tensor_link_based)(self[, ...])                         | Returns the stress tensor at the current state of the simulation.                                                                                                                                                                                                                         |
| [`get_weighted_spring_lengths`](#pylimer_tools_cpp.MEHPForceBalance.get_weighted_spring_lengths)(self[, ...])                           | Get the current spring lengths (norm of vector) divided by the strand partition times the contour length.                                                                                                                                                                                 |
| [`inspect_displacement_to_mean_position_update`](#pylimer_tools_cpp.MEHPForceBalance.inspect_displacement_to_mean_position_update)(...) | Helper method to debug and/or understand what happens to certain links when being displaced.                                                                                                                                                                                              |
| [`inspect_parametrisation_optimsation_for_link`](#pylimer_tools_cpp.MEHPForceBalance.inspect_parametrisation_optimsation_for_link)(...) | Helper method to debug and/or understand what happens to certain links when being displaced and their partition updated.                                                                                                                                                                  |
| [`inspect_strand_partition_update`](#pylimer_tools_cpp.MEHPForceBalance.inspect_strand_partition_update)(self, link_idx)                | Helper method to debug and/or understand what happens to certain links when the strand partition is being updated.                                                                                                                                                                        |
| [`move_sliplinks_to_their_best_branch`](#pylimer_tools_cpp.MEHPForceBalance.move_sliplinks_to_their_best_branch)(self, ...)             |                                                                                                                                                                                                                                                                                           |
| [`randomly_add_sliplinks`](#pylimer_tools_cpp.MEHPForceBalance.randomly_add_sliplinks)(self, ...[, cutoff, ...])                        | Randomly sample and add slip-links based on certain criteria.                                                                                                                                                                                                                             |
| [`run_force_relaxation`](#pylimer_tools_cpp.MEHPForceBalance.run_force_relaxation)(self[, ...])                                         | Run the simulation.                                                                                                                                                                                                                                                                       |
| [`set_displacements`](#pylimer_tools_cpp.MEHPForceBalance.set_displacements)(self, arg0)                                                | Set the current link displacements.                                                                                                                                                                                                                                                       |
| [`set_strand_contour_lengths`](#pylimer_tools_cpp.MEHPForceBalance.set_strand_contour_lengths)(self, arg0)                              | Set/overwrite the contour lengths.                                                                                                                                                                                                                                                        |
| [`set_strand_partitions`](#pylimer_tools_cpp.MEHPForceBalance.set_strand_partitions)(self, arg0)                                        | Set the current strand partitions.                                                                                                                                                                                                                                                        |
| [`swap_sliplinks_incl_xlinks`](#pylimer_tools_cpp.MEHPForceBalance.swap_sliplinks_incl_xlinks)(self, arg0, arg1, ...)                   |                                                                                                                                                                                                                                                                                           |
| [`validate_network`](#pylimer_tools_cpp.MEHPForceBalance.validate_network)(self)                                                        | Validates the internal structures.                                                                                                                                                                                                                                                        |

### Attributes Documentation

#### network

### Methods Documentation

#### add_sliplinks(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), strand_idx_1: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], strand_idx_2: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], x: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], y: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], z: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], alpha_1: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], alpha_2: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)], clamp_alpha: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Add new slip-links

#### add_sliplinks_based_on_cycles(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), maxLoopLength: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [int](https://docs.python.org/3/library/functions.html#int)

Detect and add slip-links based on detected entanglements.

#### WARNING
Does not work yet.

#### config_assume_box_large_enough(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), box_large_enough: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Configure whether to run PBC on the bonds or not.

* **Parameters:**
  **box_large_enough** – If True, assume the box is large enough and don’t apply PBC on bonds.
  If your bonds could get larger than half the box length,
  this must be kept False (default).

#### config_entanglement_type(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [None](https://docs.python.org/3/library/constants.html#None)

To have certain crosslinks behave as entanglements in the removal process,
you can specify here a type, that you have used in the universe to specify:
- the type of entanglement atoms (expected with functionality f = 3),
- and the entanglement-bonds between the entanglement atoms.

I.e., say you want to model some entanglements as non-slipping,
bonds between two strand beads resulting in f = 3 beads, for example,
you can call this method to have the “StructureSimplificationMode” also remove these atoms,
if they have a functionality of 2 or less while still being connected to its partner bead.

* **Parameters:**
  **type** – The atom type to treat as entanglement.

#### config_mean_bond_length(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), b: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Configure the $b$ used e.g. for the topological Gamma-factor.

* **Parameters:**
  **b** – The mean bond length to use.

#### config_simplification_frequency(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), frequency: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 10) → [None](https://docs.python.org/3/library/constants.html#None)

Configure every how many steps to simplify the structure.

* **Parameters:**
  **frequency** – The number of steps between simplification. Default: 10.

#### config_spring_breaking_distance(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), distance_over_contour_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = -1) → [None](https://docs.python.org/3/library/constants.html#None)

Configure the “force” (distance over contour length) at which the bonds break.
Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
Springs breaking will happen before the simplification procedure is run.
Negative values will disable spring breaking.
Default: -1.

* **Parameters:**
  **distance_over_contour_length** – The threshold for breaking springs.

#### config_spring_constant(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), kappa: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Configure the spring constant used in the simulation.

* **Parameters:**
  **kappa** – The spring constant.

#### config_step_output(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), output_configuration: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.OutputConfiguration](pylimer_tools_cpp.OutputConfiguration.md#pylimer_tools_cpp.OutputConfiguration)]) → [None](https://docs.python.org/3/library/constants.html#None)

Set which values to log during the simulation.

* **Parameters:**
  **output_configuration** – A list of OutputConfiguration structs
  specifying what values to log and how often

#### *static* construct_with_random_sliplinks(universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), nr_of_sliplinks_to_sample: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), upper_sampling_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.2, lower_sampling_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.0, minimum_nr_of_sliplinks: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0, same_strand_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 3, seed: [str](https://docs.python.org/3/library/stdtypes.html#str) = '', crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, is_2d: [bool](https://docs.python.org/3/library/functions.html#bool) = False, skip_dangling_soluble_entanglements: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)

Instantiate this simulator with randomly chosen slip-links.

#### deform_to(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), new_box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [None](https://docs.python.org/3/library/constants.html#None)

Perform a deformation of the system box to a different box.
All coordinates etc. will be scaled as needed.

#### evaluate_spring_distance(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), network: [pylimer_tools_cpp.SimplifiedBalanceNetwork](pylimer_tools_cpp.SimplifiedBalanceNetwork.md#pylimer_tools_cpp.SimplifiedBalanceNetwork), displacements: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], spring_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate the distance vector for a specific spring.

* **Parameters:**
  * **network** – The force balance network
  * **displacements** – Current displacement vector
  * **spring_idx** – Index of the spring to evaluate
* **Returns:**
  3D distance vector for the spring

#### evaluate_spring_distance_from(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), network: [pylimer_tools_cpp.SimplifiedBalanceNetwork](pylimer_tools_cpp.SimplifiedBalanceNetwork.md#pylimer_tools_cpp.SimplifiedBalanceNetwork), displacements: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], spring_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate the spring distance from a specific link.

* **Parameters:**
  * **network** – The force balance network
  * **displacements** – Current displacement vector
  * **spring_idx** – Index of the spring to evaluate
  * **link_idx** – Index of the starting link
* **Returns:**
  3D distance vector from the specified link

#### evaluate_spring_distance_to(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), network: [pylimer_tools_cpp.SimplifiedBalanceNetwork](pylimer_tools_cpp.SimplifiedBalanceNetwork.md#pylimer_tools_cpp.SimplifiedBalanceNetwork), displacements: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], spring_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate the spring distance to a specific link.

* **Parameters:**
  * **network** – The force balance network
  * **displacements** – Current displacement vector
  * **spring_idx** – Index of the spring to evaluate
  * **link_idx** – Index of the target link
* **Returns:**
  3D distance vector to the specified link

#### get_average_strand_length(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the average length of the strands. Note that in contrast to [`get_gamma_factor()`](#pylimer_tools_cpp.MEHPForceBalance.get_gamma_factor),
this value is normalized by the number of strands rather than the number of chains.

#### get_coordinates(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current coordinates of the crosslinkers and entanglement links.

#### get_crosslinker_universe(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

Returns the universe [of crosslinkers] with the positions of the current state of the simulation.

#### get_current_spring_lengths(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Get the spring distances.

#### get_current_spring_vectors(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the spring vectors.

#### get_current_strand_vectors(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_dangling_weight_fraction(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of non-active strands

Caution: ignores atom masses.

#### get_default_mean_bond_length(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [float](https://docs.python.org/3/library/functions.html#float)

Returns the value effectively used in [`get_gamma_factor()`](#pylimer_tools_cpp.MEHPForceBalance.get_gamma_factor) for
$b$ in $\langle R_{0,\eta}^2 = N_{\eta} b^2\rangle$.

#### get_displacement_residual_norm(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [float](https://docs.python.org/3/library/functions.html#float)

Get the current link displacement residual norm.

#### get_displacements(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current link displacements.

#### get_effective_functionality_of_atoms(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [int](https://docs.python.org/3/library/functions.html#int)]

Returns the number of active strands connected to each atom, atomId used as index

* **Parameters:**
  **tolerance** – strands under this length are considered inactive

#### get_exit_reason(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [pylimer_tools_cpp.ExitReason](pylimer_tools_cpp.ExitReason.md#pylimer_tools_cpp.ExitReason)

Returns the reason for termination of the simulation

#### get_force_magnitude_vector(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate the norm of the force on each (slip- or cross-) link.

#### get_force_on(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate the force on a particular (slip- or cross-) link.

#### get_gamma_factor(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = -1.0, nr_of_chains: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1, one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [float](https://docs.python.org/3/library/functions.html#float)

Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

$\Gamma = \langle\gamma_{\eta}\rangle$, with $\gamma_{\eta} = \\frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}$,
which you can use as $G_{\mathrm{ANT}} = \Gamma \nu k_B T$,
where $\eta$ is the index of a particular strand,
$R_{0}^2$ is the melt mean square end to end distance, in phantom systems $$= N_{\eta}*b^2$$
$N_{\eta}$ is the number of atoms in this strand $\eta$,
$b$ its mean square bond length,
$\nu$ the number density of network strands,
$T$ the temperature and
$k_B$ Boltzmann’s constant.

* **Parameters:**
  * **b02** – The melt $<b>_0^2$: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it’s the slope in a <R_0^2> vs. N plot, also sometimes labelled $C_\infinity b^2$.
  * **nr_of_chains** – The value to normalize the sum of square distances by. Usually (and default if $< 0$) the nr of springs.

#### get_gamma_factors(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)

#### get_gamma_factors_in_dir(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), direction: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)

* **Parameters:**
  * **b02** – The melt $<b>_0^2$: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it’s the slope in a <R_0^2> vs. N plot, also sometimes labelled $C_\infinity b^2$.
  * **direction** – The direction in which to compute the gamma factors (0: x, 1: y, 2: z)

#### get_ids_of_active_nodes(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the atom ids of the nodes that are considered active.
Only crosslink ids are returned (not e.g. entanglement links).

* **Parameters:**
  **tolerance** – strands under this length are considered inactive. A node is active if it has > 1 active strands.

#### get_initial_coordinates(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the initial coordinates of the (remaining) crosslinkers and entanglement links.

#### *static* get_neighbour_link_indices(network: [pylimer_tools_cpp.SimplifiedBalanceNetwork](pylimer_tools_cpp.SimplifiedBalanceNetwork.md#pylimer_tools_cpp.SimplifiedBalanceNetwork), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the indices of neighboring links for a given link.

* **Parameters:**
  * **network** – The force balance network
  * **link_idx** – Index of the link to find neighbors for
* **Returns:**
  Vector of connected link indices

#### get_nr_of_active_nodes(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of active nodes (incl. entanglement nodes [atoms with type = entanglementType, present in the universe when creating this simulator],
excl. entanglement links [the slip-links created internally when e.g. constructing the simulator with random slip-links]).

* **Parameters:**
  **tolerance** – strands under this length are considered inactive. A node is active if it has > 1 active strands.

#### get_nr_of_active_springs(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

> Get the number of active springs remaining after running the simulation.
* **Parameters:**
  **tolerance** – springs under this length are considered inactive

#### get_nr_of_active_strands(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

> Get the number of active strands remaining after running the simulation.
* **Parameters:**
  **tolerance** – strands under this length are considered inactive

#### get_nr_of_active_strands_in_dir(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), direction: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

> Get the number of active strands remaining after running the simulation.
* **Parameters:**
  * **direction** – The direction in which to compute the active strands (0: x, 1: y, 2: z)
  * **tolerance** – strands under this length are considered inactive

#### get_nr_of_atoms(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_bonds(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_extra_atoms(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_extra_bonds(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_intra_chain_sliplinks(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_iterations(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [int](https://docs.python.org/3/library/functions.html#int)

Returns the number of iterations used for force relaxation so far.

#### get_nr_of_nodes(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of nodes (crosslinkers) considered in this simulation.

#### get_nr_of_strands(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of strands considered in this simulation.

#### get_overall_strand_lengths(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Get the sum of the lengths of the springs of each strand.

#### get_pressure(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [float](https://docs.python.org/3/library/functions.html#float)

Returns the pressure at the current state of the simulation.

#### get_soluble_weight_fraction(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of strands connected to active
strands (any depth).

Caution: ignores atom masses.

#### get_strand_partition_indices_of_sliplink(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), network: [pylimer_tools_cpp.SimplifiedBalanceNetwork](pylimer_tools_cpp.SimplifiedBalanceNetwork.md#pylimer_tools_cpp.SimplifiedBalanceNetwork), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the indices of strand partitions associated with a slip-link.

* **Parameters:**
  * **network** – The force balance network
  * **link_idx** – Index of the slip-link
* **Returns:**
  Vector of strand partition indices

#### get_strand_partitions(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current strand partitions (the fraction of the contour length associated with each spring).

#### get_stress_on(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Evaluate the stress on a particular (slip- or cross-) link.

#### get_stress_tensor(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Returns the stress tensor at the current state of the simulation.

The units are in $[\text{units of }\kappa]/[\text{distance units}]$,
where the units of $\kappa$ should be $[\text{force}]/[\text{distance units}]^2$.
Make sure to multiply by $\kappa$ or configure it appropriately.

#### get_stress_tensor_link_based(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, xlinks_only: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Returns the stress tensor at the current state of the simulation.

#### get_weighted_spring_lengths(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current spring lengths (norm of vector) divided by the strand partition times the contour length.

#### inspect_displacement_to_mean_position_update(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Helper method to debug and/or understand what happens to certain links when being displaced.

#### inspect_parametrisation_optimsation_for_link(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), displacements: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], strand_partitions: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], max_nr_of_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 100, alpha_tol: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1e-09, min_nr_of_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1, one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]'], [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]'], [int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float)]

Helper method to debug and/or understand what happens to certain links
when being displaced and their partition updated.

#### inspect_strand_partition_update(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Helper method to debug and/or understand what happens to certain links
when the strand partition is being updated.

#### move_sliplinks_to_their_best_branch(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), arg0: [pylimer_tools_cpp.SimplifiedBalanceNetwork](pylimer_tools_cpp.SimplifiedBalanceNetwork.md#pylimer_tools_cpp.SimplifiedBalanceNetwork), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], arg3: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), arg4: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), arg5: [bool](https://docs.python.org/3/library/functions.html#bool), arg6: [bool](https://docs.python.org/3/library/functions.html#bool)) → [None](https://docs.python.org/3/library/constants.html#None)

#### randomly_add_sliplinks(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), nr_of_sliplinks_to_sample: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 2.0, minimum_nr_of_sliplinks: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0, same_strand_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 2, exclude_crosslinks: [bool](https://docs.python.org/3/library/functions.html#bool) = False, seed: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [int](https://docs.python.org/3/library/functions.html#int)

Randomly sample and add slip-links based on certain criteria.

#### run_force_relaxation(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), max_nr_of_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 250000, x_tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1e-12, initial_residual_norm: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = -1.0, simplification_mode: [pylimer_tools_cpp.StructureSimplificationMode](pylimer_tools_cpp.StructureSimplificationMode.md#pylimer_tools_cpp.StructureSimplificationMode) = StructureSimplificationMode.NO_SIMPLIFICATION, inactive_removal_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001, do_inner_iterations: [bool](https://docs.python.org/3/library/functions.html#bool) = False, allow_sliplinks_to_pass_each_other: [pylimer_tools_cpp.LinkSwappingMode](pylimer_tools_cpp.LinkSwappingMode.md#pylimer_tools_cpp.LinkSwappingMode) = LinkSwappingMode.NO_SWAPPING, swapping_frequency: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 10, one_over_strand_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, nr_of_crosslink_swaps_allowed_per_sliplink: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1, disable_slipping: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Run the simulation.
Note that the final state of the minimization is persisted and reused if you use this method again.
This is useful if you want to run a global optimization first and add a local one afterwards.
As a consequence though, you cannot simply benchmark only this method; you must include the setup.

* **Parameters:**
  * **max_nr_of_steps** – The maximum number of steps to do during the simulation.
  * **x_tolerance** – The tolerance of the displacements as an exit condition.
  * **initial_residual_norm** – The residual norm relative to which the relative tolerance is specified. Negative values mean, it will be replaced with the current one.
  * **simplification_mode** – How to simplify the structure during the minimization.
  * **inactive_removal_cutoff** – The tolerance in distance units of the partial spring length to count as active, relevant if simplification mode is specified to be something other than NO_SIMPLIFICATION.
  * **do_inner_iterations** – Whether to do inner iterations; usually, they are not helpful.
  * **allow_sliplinks_to_pass_each_other** – Whether slip-links can pass each other.
  * **swapping_frequency** – How often slip-links attempt to swap.
  * **one_over_strand_partition_upper_limit** – Super-secret parameter. Use 1, gradually increase (and then -1) if you want to publish.
  * **nr_of_crosslink_swaps_allowed_per_sliplink** – Use to steer whether slip-links can cross crosslinks when swapping is enabled.
  * **disable_slipping** – Whether slip-links should be prohibited from slipping.

#### set_displacements(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [None](https://docs.python.org/3/library/constants.html#None)

Set the current link displacements.

#### set_strand_contour_lengths(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [None](https://docs.python.org/3/library/constants.html#None)

Set/overwrite the contour lengths.

#### set_strand_partitions(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [None](https://docs.python.org/3/library/constants.html#None)

Set the current strand partitions.

#### swap_sliplinks_incl_xlinks(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance), arg0: [pylimer_tools_cpp.SimplifiedBalanceNetwork](pylimer_tools_cpp.SimplifiedBalanceNetwork.md#pylimer_tools_cpp.SimplifiedBalanceNetwork), arg1: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], arg2: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], arg3: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), arg4: [bool](https://docs.python.org/3/library/functions.html#bool)) → [None](https://docs.python.org/3/library/constants.html#None)

#### validate_network(self: [pylimer_tools_cpp.MEHPForceBalance](#pylimer_tools_cpp.MEHPForceBalance)) → [bool](https://docs.python.org/3/library/functions.html#bool)

Validates the internal structures.

Throws an error if something is not ok.
Otherwise, it returns true.

Can be used e.g. as `assert fb.validate_network()`.

* **Returns:**
  True if validation passes, raises exception otherwise
