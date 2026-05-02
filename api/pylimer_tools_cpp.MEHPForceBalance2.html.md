# MEHPForceBalance2

### *class* pylimer_tools_cpp.MEHPForceBalance2(\*args, \*\*kwargs)

Bases: `pybind11_object`

A small simulation tool for quickly minimizing the force between the crosslinker beads.

This is the third implementation of the MEHP.
It’s the fastest implementation by using a simple spring model only, disabling the non-linear
periodic boundary conditions, and instead builds a sparse linear system of equations that’s readily solved.
However, it allows entanglements to be represented as additional links or/and springs,
although without slipping along the chain.

Please cite Bernhard and Gusev [[BG25](../acknowledgements.md#id12)] if you use this method.

Overloaded function.

1. \_\_init_\_(self: pylimer_tools_cpp.MEHPForceBalance2, universe: pylimer_tools_cpp.Universe, crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False) -> None
   > Instantiate the simulator for a certain universe.
   > * **param universe:**
   >   The universe giving the basic connectivity to compute with
   > * **param crosslinker_type:**
   >   The atom type of the crosslinkers. Needed to reduce the network.
   > * **param is_2d:**
   >   Whether to ignore the z direction.
   > * **param kappa:**
   >   The spring constant
   > * **param remove_2functionalCrosslinkers:**
   >   whether to keep or remove the 2-functional crosslinkers when setting up the network
   > * **param remove_dangling_chains:**
   >   whether to keep or remove obviously dangling chains when setting up the network
2. \_\_init_\_(self: pylimer_tools_cpp.MEHPForceBalance2, universe: pylimer_tools_cpp.Universe, entanglements: pylimer_tools_cpp.AtomPairEntanglements, crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False, entanglements_as_springs: bool = False) -> None
   > Instantiate the simulator for a certain universe with the given entanglements.
   > * **param universe:**
   >   The universe giving the basic connectivity to compute with
   > * **param entanglements:**
   >   The entanglements to use in the computation
   > * **param crosslinker_type:**
   >   The atom type of the crosslinkers. Needed to reduce the network.
   > * **param is_2d:**
   >   Whether to ignore the z direction.
   > * **param entanglements_as_springs:**
   >   whether to use the entanglements as springs instead of links
3. \_\_init_\_(self: pylimer_tools_cpp.MEHPForceBalance2, universe: pylimer_tools_cpp.Universe, nr_of_entanglements_to_sample: typing.SupportsInt, upper_sampling_cutoff: typing.SupportsFloat = 1.2, lower_sampling_cutoff: typing.SupportsFloat = 0.0, minimum_nr_of_sliplinks: typing.SupportsInt = 0, same_strand_cutoff: typing.SupportsFloat = 3, seed: str = ‘’, crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False, skip_dangling_soluble_entanglements: bool = True, entanglements_as_springs: bool = True) -> None
   > Instantiate this simulator with randomly chosen slip-links.
   > * **param universe:**
   >   The universe containing the basic atoms and connectivity
   > * **param nr_of_entanglements_to_sample:**
   >   The number of entanglements to sample
   > * **param upper_cutoff:**
   >   maximum distance from one sampled bead to its partner
   > * **param lower_cutoff:**
   >   minimum distance from one sampled bead to its partner
   > * **param minimum_nr_of_entanglements:**
   >   The minimum number of entanglements that should be sampled
   > * **param same_strand_cutoff:**
   >   distance from one sampled bead to its pair within the same strand
   > * **param seed:**
   >   The seed for the random number generator
   > * **param cross_linker_type:**
   > * **param is_2d:**
   > * **param filter_entanglements:**
   > * **param entanglements_as_springs:**
   >   whether to model the entanglements as merged beads or beads with 1 spring in between

### Attributes Summary

| [`network`](#pylimer_tools_cpp.MEHPForceBalance2.network)   |    |
|-------------------------------------------------------------|----|

### Methods Summary

| [`config_assume_network_is_complete`](#pylimer_tools_cpp.MEHPForceBalance2.config_assume_network_is_complete)(self[, ...])               | Configure whether the network is assumed to be complete.                                                                                           |
|------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| [`config_mean_bond_length`](#pylimer_tools_cpp.MEHPForceBalance2.config_mean_bond_length)(self[, b])                                     | Configure the $b$ used e.g. for the topological Gamma-factor.                                                                                      |
| [`config_spring_breaking_distance`](#pylimer_tools_cpp.MEHPForceBalance2.config_spring_breaking_distance)(self[, ...])                   | Configure the "force" (distance over contour length) at which the bonds break.                                                                     |
| [`config_spring_constant`](#pylimer_tools_cpp.MEHPForceBalance2.config_spring_constant)(self[, kappa])                                   | Configure the spring constant used in the simulation.                                                                                              |
| [`config_step_output`](#pylimer_tools_cpp.MEHPForceBalance2.config_step_output)(self, values)                                            | Set which values to log.                                                                                                                           |
| [`deform_to`](#pylimer_tools_cpp.MEHPForceBalance2.deform_to)(self, new_box)                                                             | Perform a deformation of the system box to a different box.                                                                                        |
| [`evaluate_spring_vector`](#pylimer_tools_cpp.MEHPForceBalance2.evaluate_spring_vector)(self, network, ...)                              | Evaluate the distance vector for a specific spring.                                                                                                |
| [`evaluate_spring_vector_from`](#pylimer_tools_cpp.MEHPForceBalance2.evaluate_spring_vector_from)(self, network, ...)                    | Evaluate the spring vector in the direction from a specific link.                                                                                  |
| [`evaluate_spring_vector_to`](#pylimer_tools_cpp.MEHPForceBalance2.evaluate_spring_vector_to)(self, network, ...)                        | Evaluate the spring vector in the direction to a specific link.                                                                                    |
| [`get_average_spring_length`](#pylimer_tools_cpp.MEHPForceBalance2.get_average_spring_length)(self)                                      | Get the average length of the springs.                                                                                                             |
| [`get_coordinates`](#pylimer_tools_cpp.MEHPForceBalance2.get_coordinates)(self)                                                          | Get the current coordinates of the crosslinkers and entanglement links.                                                                            |
| [`get_crosslinker_universe`](#pylimer_tools_cpp.MEHPForceBalance2.get_crosslinker_universe)(self)                                        | Returns the universe [of crosslinkers] with the positions of the current state of the simulation.                                                  |
| [`get_current_spring_lengths`](#pylimer_tools_cpp.MEHPForceBalance2.get_current_spring_lengths)(self)                                    | Get the spring lengths (Euclidean distances of the spring vectors).                                                                                |
| [`get_current_spring_vectors`](#pylimer_tools_cpp.MEHPForceBalance2.get_current_spring_vectors)(self)                                    | Get the spring vectors.                                                                                                                            |
| [`get_dangling_weight_fraction`](#pylimer_tools_cpp.MEHPForceBalance2.get_dangling_weight_fraction)(self[, tolerance])                   | Compute the weight fraction of non-active springs                                                                                                  |
| [`get_default_mean_bond_length`](#pylimer_tools_cpp.MEHPForceBalance2.get_default_mean_bond_length)(self)                                | Returns the value effectively used in `getGammaFactor()` for $b$ in $\langle R_{0,\eta}^2 = N_{\eta} b^2\rangle$.                                  |
| [`get_displacement_residual_norm`](#pylimer_tools_cpp.MEHPForceBalance2.get_displacement_residual_norm)(self)                            | Get the current link displacement residual norm.                                                                                                   |
| [`get_displacements`](#pylimer_tools_cpp.MEHPForceBalance2.get_displacements)(self)                                                      | Get the current link displacements.                                                                                                                |
| [`get_effective_functionality_of_atoms`](#pylimer_tools_cpp.MEHPForceBalance2.get_effective_functionality_of_atoms)(self[, ...])         | Returns the number of active springs connected to each atom, atomId used as index                                                                  |
| [`get_exit_reason`](#pylimer_tools_cpp.MEHPForceBalance2.get_exit_reason)(self)                                                          | Returns the reason for termination of the simulation                                                                                               |
| [`get_force_magnitude_vector`](#pylimer_tools_cpp.MEHPForceBalance2.get_force_magnitude_vector)(self)                                    | Evaluate the norm of the force on each (slip- or cross-) link.                                                                                     |
| [`get_force_on`](#pylimer_tools_cpp.MEHPForceBalance2.get_force_on)(self, link_idx)                                                      | Evaluate the force on a particular (slip- or cross-) link.                                                                                         |
| [`get_gamma_factor`](#pylimer_tools_cpp.MEHPForceBalance2.get_gamma_factor)(self[, b02, nr_of_chains])                                   | Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:                                                                                 |
| [`get_gamma_factors`](#pylimer_tools_cpp.MEHPForceBalance2.get_gamma_factors)(self, b02)                                                 | Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)                            |
| [`get_gamma_factors_in_dir`](#pylimer_tools_cpp.MEHPForceBalance2.get_gamma_factors_in_dir)(self, b02, direction)                        | Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02) |
| [`get_ids_of_active_nodes`](#pylimer_tools_cpp.MEHPForceBalance2.get_ids_of_active_nodes)(self[, tolerance])                             | Get the atom ids of the nodes that are considered active.                                                                                          |
| [`get_initial_coordinates`](#pylimer_tools_cpp.MEHPForceBalance2.get_initial_coordinates)(self)                                          | Get the initial coordinates of the crosslinkers and entanglement links.                                                                            |
| [`get_neighbour_link_indices`](#pylimer_tools_cpp.MEHPForceBalance2.get_neighbour_link_indices)(network, link_idx)                       | Get the indices of neighboring links for a given link.                                                                                             |
| [`get_nr_of_active_nodes`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_active_nodes)(self[, tolerance])                               | Get the number of active nodes (incl.                                                                                                              |
| [`get_nr_of_active_spring`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_active_spring)(self[, tolerance])                             | Get the number of active springs remaining after running the simulation.                                                                           |
| [`get_nr_of_active_springs`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_active_springs)(self[, tolerance])                           | Get the number of active springs remaining after running the simulation.                                                                           |
| [`get_nr_of_active_springs_in_dir`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_active_springs_in_dir)(self, direction)               | Get the number of active springs remaining after running the simulation.                                                                           |
| [`get_nr_of_atoms`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_atoms)(self)                                                          |                                                                                                                                                    |
| [`get_nr_of_bonds`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_bonds)(self)                                                          |                                                                                                                                                    |
| [`get_nr_of_extra_atoms`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_extra_atoms)(self)                                              |                                                                                                                                                    |
| [`get_nr_of_extra_bonds`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_extra_bonds)(self)                                              |                                                                                                                                                    |
| [`get_nr_of_intra_chain_sliplinks`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_intra_chain_sliplinks)(self)                          |                                                                                                                                                    |
| [`get_nr_of_iterations`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_iterations)(self)                                                | Returns the number of iterations used for force relaxation so far.                                                                                 |
| [`get_nr_of_nodes`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_nodes)(self)                                                          | Get the number of nodes (crosslinkers) considered in this simulation.                                                                              |
| [`get_nr_of_springs`](#pylimer_tools_cpp.MEHPForceBalance2.get_nr_of_springs)(self)                                                      | Get the number of springs considered in this simulation.                                                                                           |
| [`get_overall_strand_lengths`](#pylimer_tools_cpp.MEHPForceBalance2.get_overall_strand_lengths)(self)                                    | Get the sum of the lengths of the springs of each strand.                                                                                          |
| [`get_pressure`](#pylimer_tools_cpp.MEHPForceBalance2.get_pressure)(self)                                                                | Returns the pressure at the current state of the simulation.                                                                                       |
| [`get_soluble_weight_fraction`](#pylimer_tools_cpp.MEHPForceBalance2.get_soluble_weight_fraction)(self[, tolerance])                     | Compute the weight fraction of springs connected to active springs (any depth).                                                                    |
| [`get_stress_on`](#pylimer_tools_cpp.MEHPForceBalance2.get_stress_on)(self, link_idx)                                                    | Evaluate the stress on a particular (slip- or cross-) link.                                                                                        |
| [`get_stress_tensor`](#pylimer_tools_cpp.MEHPForceBalance2.get_stress_tensor)(self)                                                      | Returns the stress tensor at the current state of the simulation.                                                                                  |
| [`get_stress_tensor_link_based`](#pylimer_tools_cpp.MEHPForceBalance2.get_stress_tensor_link_based)(self[, xlinks_only])                 | Returns the stress tensor at the current state of the simulation.                                                                                  |
| [`get_weighted_spring_lengths`](#pylimer_tools_cpp.MEHPForceBalance2.get_weighted_spring_lengths)(self)                                  | Get the current spring lengths (norm of vector) divided by the spring partition times the contour length.                                          |
| [`inspect_displacement_to_mean_position_update`](#pylimer_tools_cpp.MEHPForceBalance2.inspect_displacement_to_mean_position_update)(...) | Helper method to debug and/or understand what happens to certain links when being displaced.                                                       |
| [`run_force_relaxation`](#pylimer_tools_cpp.MEHPForceBalance2.run_force_relaxation)(self[, ...])                                         | Run the simulation.                                                                                                                                |
| [`set_displacements`](#pylimer_tools_cpp.MEHPForceBalance2.set_displacements)(self, arg0)                                                | Set the current link displacements.                                                                                                                |
| [`set_spring_contour_lengths`](#pylimer_tools_cpp.MEHPForceBalance2.set_spring_contour_lengths)(self, arg0)                              | Set/overwrite the contour lengths.                                                                                                                 |
| [`validate_network`](#pylimer_tools_cpp.MEHPForceBalance2.validate_network)(self)                                                        | Validates the internal structures.                                                                                                                 |

### Attributes Documentation

#### network

### Methods Documentation

#### config_assume_network_is_complete(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), assume_network_is_complete: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Configure whether the network is assumed to be complete.

This assumption means, that the universe instance is not queried for clusters when
computing the fraction of e.g. the soluble or dangling strands.
Do not set this to true if inactive
(dangling or free) strands have been removed from the network.

This is useful to reduce memory between MC generator and force balance,
since the universe representation with all the in-between beads does not need to be stored.
However, currently, the removal procedures are not loss-free.
Therefore, use this with care.

#### config_mean_bond_length(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), b: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Configure the $b$ used e.g. for the topological Gamma-factor.

* **Parameters:**
  **b** – The mean bond length to use.

#### config_spring_breaking_distance(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), distance_over_contour_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = -1) → [None](https://docs.python.org/3/library/constants.html#None)

Configure the “force” (distance over contour length) at which the bonds break.
Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
Springs breaking will happen before the simplification procedure is run.
Negative values will disable spring breaking.
Default: -1.

* **Parameters:**
  **distance_over_contour_length** – The threshold for breaking springs.

#### config_spring_constant(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), kappa: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Configure the spring constant used in the simulation.

* **Parameters:**
  **kappa** – The spring constant.

#### config_step_output(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), values: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.OutputConfiguration](pylimer_tools_cpp.OutputConfiguration.md#pylimer_tools_cpp.OutputConfiguration)]) → [None](https://docs.python.org/3/library/constants.html#None)

Set which values to log.

* **Parameters:**
  **values** – a list of OutputConfiguration structs

#### deform_to(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), new_box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [None](https://docs.python.org/3/library/constants.html#None)

Perform a deformation of the system box to a different box.
All coordinates etc. will be scaled as needed.

#### evaluate_spring_vector(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), network: [pylimer_tools_cpp.SimplifiedBalance2Network](pylimer_tools_cpp.SimplifiedBalance2Network.md#pylimer_tools_cpp.SimplifiedBalance2Network), displacements: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], spring_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate the distance vector for a specific spring.

* **Parameters:**
  * **network** – The force balance 2 network
  * **displacements** – Current displacement vector
  * **spring_idx** – Index of the spring to evaluate
* **Returns:**
  3D distance vector for the spring

#### evaluate_spring_vector_from(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), network: [pylimer_tools_cpp.SimplifiedBalance2Network](pylimer_tools_cpp.SimplifiedBalance2Network.md#pylimer_tools_cpp.SimplifiedBalance2Network), displacements: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], spring_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate the spring vector in the direction from a specific link.

* **Parameters:**
  * **network** – The force balance 2 network
  * **displacements** – Current displacement vector
  * **spring_idx** – Index of the spring to evaluate
  * **link_idx** – Index of the starting link
* **Returns:**
  3D distance vector from the specified link

#### evaluate_spring_vector_to(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), network: [pylimer_tools_cpp.SimplifiedBalance2Network](pylimer_tools_cpp.SimplifiedBalance2Network.md#pylimer_tools_cpp.SimplifiedBalance2Network), displacements: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], spring_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate the spring vector in the direction to a specific link.

* **Parameters:**
  * **network** – The force balance 2 network
  * **displacements** – Current displacement vector
  * **spring_idx** – Index of the spring to evaluate
  * **link_idx** – Index of the target link
* **Returns:**
  3D distance vector to the specified link

#### get_average_spring_length(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the average length of the springs. Note that in contrast to `getGammaFactor()`,
this value is normalized by the number of springs rather than the number of chains.

#### get_coordinates(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current coordinates of the crosslinkers and entanglement links.

#### get_crosslinker_universe(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

Returns the universe [of crosslinkers] with the positions of the current state of the simulation.

#### get_current_spring_lengths(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Get the spring lengths (Euclidean distances of the spring vectors).

#### get_current_spring_vectors(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the spring vectors.

#### get_dangling_weight_fraction(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of non-active springs

Caution: ignores atom masses.

#### get_default_mean_bond_length(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [float](https://docs.python.org/3/library/functions.html#float)

Returns the value effectively used in `getGammaFactor()` for
$b$ in $\langle R_{0,\eta}^2 = N_{\eta} b^2\rangle$.

#### get_displacement_residual_norm(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the current link displacement residual norm.

#### get_displacements(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current link displacements.

#### get_effective_functionality_of_atoms(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [int](https://docs.python.org/3/library/functions.html#int)]

Returns the number of active springs connected to each atom, atomId used as index

* **Parameters:**
  **tolerance** – springs under this length are considered inactive

#### get_exit_reason(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [pylimer_tools_cpp.ExitReason](pylimer_tools_cpp.ExitReason.md#pylimer_tools_cpp.ExitReason)

Returns the reason for termination of the simulation

#### get_force_magnitude_vector(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate the norm of the force on each (slip- or cross-) link.

#### get_force_on(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate the force on a particular (slip- or cross-) link.

#### get_gamma_factor(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = -1.0, nr_of_chains: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [float](https://docs.python.org/3/library/functions.html#float)

Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

$\Gamma = \langle\gamma_{\eta}\rangle$, with $\gamma_{\eta} = \\frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}$,
which you can use as $G_{\mathrm{ANT}} = \Gamma \nu k_B T$,
where $\eta$ is the index of a particular strand,
$R_{0}^2$ is the melt mean square end to end distance, in phantom systems $$= N_{\eta}*b^2$$
$N_{\eta}$ is the number of atoms in this strand $\eta$,
$b$ its mean square bond length,
$T$ the temperature and
$k_B$ Boltzmann’s constant.

* **Parameters:**
  * **b02** – The melt $<b>_0^2$: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it’s the slope in a <R_0^2> vs. N plot, also sometimes labelled $C_\infinity b^2$.
  * **nr_of_chains** – The value to normalize the sum of square distances by. Usually (and default if $< 0$) the nr of springs.

#### get_gamma_factors(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)

#### get_gamma_factors_in_dir(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), direction: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)

* **Parameters:**
  * **b02** – The melt $<b>_0^2$: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it’s the slope in a <R_0^2> vs. N plot, also sometimes labelled $C_\infinity b^2$.
  * **direction** – The direction in which to compute the gamma factors (0: x, 1: y, 2: z)

#### get_ids_of_active_nodes(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the atom ids of the nodes that are considered active.
Only crosslink ids are returned (not e.g. entanglement links).

* **Parameters:**
  **tolerance** – springs under this length are considered inactive. A node is active if it has > 1 active springs.

#### get_initial_coordinates(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the initial coordinates of the crosslinkers and entanglement links.

#### *static* get_neighbour_link_indices(network: [pylimer_tools_cpp.SimplifiedBalance2Network](pylimer_tools_cpp.SimplifiedBalance2Network.md#pylimer_tools_cpp.SimplifiedBalance2Network), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the indices of neighboring links for a given link.

* **Parameters:**
  * **network** – The force balance 2 network
  * **link_idx** – Index of the link to find neighbors for
* **Returns:**
  Vector of connected link indices

#### get_nr_of_active_nodes(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of active nodes (incl. entanglement nodes [atoms with type = entanglementType, present in the universe when creating this simulator],
excl. entanglement links [the slip-links created internally when e.g. constructing the simulator with random slip-links]).

* **Parameters:**
  **tolerance** – springs under this length are considered inactive. A node is active if it has > 1 active springs.

#### get_nr_of_active_spring(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

> Get the number of active springs remaining after running the simulation.
* **Parameters:**
  **tolerance** – springs under this length are considered inactive

#### get_nr_of_active_springs(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

> Get the number of active springs remaining after running the simulation.
* **Parameters:**
  **tolerance** – springs under this length are considered inactive

#### get_nr_of_active_springs_in_dir(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), direction: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

> Get the number of active springs remaining after running the simulation.
* **Parameters:**
  * **direction** – The direction in which to compute the active springs (0: x, 1: y, 2: z)
  * **tolerance** – springs under this length are considered inactive

#### get_nr_of_atoms(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_bonds(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_extra_atoms(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_extra_bonds(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_intra_chain_sliplinks(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [int](https://docs.python.org/3/library/functions.html#int)

#### get_nr_of_iterations(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [int](https://docs.python.org/3/library/functions.html#int)

Returns the number of iterations used for force relaxation so far.

#### get_nr_of_nodes(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of nodes (crosslinkers) considered in this simulation.

#### get_nr_of_springs(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of springs considered in this simulation.

* **Parameters:**
  **tolerance** – springs under this length are considered inactive

#### get_overall_strand_lengths(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Get the sum of the lengths of the springs of each strand.

#### get_pressure(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [float](https://docs.python.org/3/library/functions.html#float)

Returns the pressure at the current state of the simulation.

#### get_soluble_weight_fraction(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of springs connected to active
springs (any depth).

Caution: ignores atom masses.

#### get_stress_on(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Evaluate the stress on a particular (slip- or cross-) link.

#### get_stress_tensor(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Returns the stress tensor at the current state of the simulation.

The units are in $[\text{units of }\kappa]/[\text{distance units}]$,
where the units of $\kappa$ should be $[\text{force}]/[\text{distance units}]^2$.
Make sure to multiply by $\kappa$ or configure it appropriately.

#### get_stress_tensor_link_based(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), xlinks_only: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Returns the stress tensor at the current state of the simulation.

#### get_weighted_spring_lengths(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current spring lengths (norm of vector) divided by the spring partition times the contour length.

#### inspect_displacement_to_mean_position_update(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Helper method to debug and/or understand what happens to certain links when being displaced.

#### run_force_relaxation(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), simplification_mode: [pylimer_tools_cpp.StructureSimplificationMode](pylimer_tools_cpp.StructureSimplificationMode.md#pylimer_tools_cpp.StructureSimplificationMode) = StructureSimplificationMode.NO_SIMPLIFICATION, inactive_removal_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1e-06, sle_solver: [pylimer_tools_cpp.SLESolver](pylimer_tools_cpp.SLESolver.md#pylimer_tools_cpp.SLESolver) = SLESolver.DEFAULT, tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1e-15, max_iterations: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 10000) → [None](https://docs.python.org/3/library/constants.html#None)

Run the simulation.

* **Parameters:**
  * **simplification_mode** – How to simplify the structure during the minimization.
  * **inactive_removal_cutoff** – The tolerance in distance units of the partial spring length to count as active, relevant if simplification mode is specified to be something other than NO_SIMPLIFICATION.
  * **sle_solver** – The solver to use for the system of linear equations (SLE).
  * **tolerance** – The stopping condition/tolerance for the SLE solver if it’s an iterative one.
  * **max_iterations** – The maximum number of iterations to perform if the SLE solver is an iterative one.

#### set_displacements(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [None](https://docs.python.org/3/library/constants.html#None)

Set the current link displacements.

#### set_spring_contour_lengths(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2), arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [None](https://docs.python.org/3/library/constants.html#None)

Set/overwrite the contour lengths.

#### validate_network(self: [pylimer_tools_cpp.MEHPForceBalance2](#pylimer_tools_cpp.MEHPForceBalance2)) → [bool](https://docs.python.org/3/library/functions.html#bool)

Validates the internal structures.

Throws an error if something is not ok.
Otherwise, it returns true.

Can be used e.g. as `assert fb.validate_network()`.

* **Returns:**
  True if validation passes, raises exception otherwise
