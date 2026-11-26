# MEHPSimulator

### *class* pylimer_tools_cpp.MEHPSimulator(\*args, \*\*kwargs)

Bases: `pybind11_object`

Unified Mechanical Equilibrium of Heterogeneous Polymers (MEHP) simulator
coordinating a canonical MEHPNetwork, a stateless MEHPForceEvaluator,
and an interchangeable MEHPSolver strategy.

Overloaded function.

1. \_\_init_\_(self: pylimer_tools_cpp.MEHPSimulator, network: pylimer_tools_cpp.MEHPNetwork, force_evaluator: pylimer_tools_cpp.MEHPForceEvaluator = None, is_2d: bool = False) -> None

Instantiate MEHPSimulator from an existing MEHPNetwork.

1. \_\_init_\_(self: pylimer_tools_cpp.MEHPSimulator, universe: pylimer_tools_cpp.Universe, entanglements: pylimer_tools_cpp.AtomPairEntanglements = AtomPairEntanglements(), crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, entanglement_mode: pylimer_tools_cpp.EntanglementMode = EntanglementMode.MERGED_LINKS, force_evaluator: pylimer_tools_cpp.MEHPForceEvaluator = None, is_2d: bool = False, keep_2functional_crosslinkers: bool = False, remove_primary_loops: bool = False, assume_box_large_enough: bool = False) -> None

Instantiate MEHPSimulator from an atomistic Universe.

1. \_\_init_\_(self: pylimer_tools_cpp.MEHPSimulator, box: pylimer_tools_cpp.Box, crosslinker_universe: pylimer_tools_cpp.CrosslinkerUniverse, entanglements: pylimer_tools_cpp.StrandPairEntanglements = StrandPairEntanglements(), crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, entanglement_mode: pylimer_tools_cpp.EntanglementMode = EntanglementMode.MERGED_LINKS, force_evaluator: pylimer_tools_cpp.MEHPForceEvaluator = None, is_2d: bool = False, remove_dangling_chains: bool = False, include_free_chains: bool = True, remove_primary_loops: bool = False, assume_box_large_enough: bool = False) -> None

Instantiate MEHPSimulator from a Box and a CrosslinkerUniverse.

### Attributes Summary

| [`force_evaluator`](#pylimer_tools_cpp.MEHPSimulator.force_evaluator)   |    |
|--------------------------------------------------------------------|----|
| [`is_2d`](#pylimer_tools_cpp.MEHPSimulator.is_2d)             |    |
| [`last_solve_result`](#pylimer_tools_cpp.MEHPSimulator.last_solve_result) |    |
| [`network`](#pylimer_tools_cpp.MEHPSimulator.network)           |    |

### Methods Summary

| [`add_slip_links`](#pylimer_tools_cpp.MEHPSimulator.add_slip_links)(\*args, \*\*kwargs)                | Overloaded function.                                                                                 |
|----------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| [`add_sliplinks`](#pylimer_tools_cpp.MEHPSimulator.add_sliplinks)(\*args, \*\*kwargs)                 | Overloaded function.                                                                                 |
| [`break_too_long_strands`](#pylimer_tools_cpp.MEHPSimulator.break_too_long_strands)(self, breaking_length)     |                                                                                                      |
| [`config_assume_box_large_enough`](#pylimer_tools_cpp.MEHPSimulator.config_assume_box_large_enough)(self[, assume])    |                                                                                                      |
| [`config_assume_network_is_complete`](#pylimer_tools_cpp.MEHPSimulator.config_assume_network_is_complete)(self[, assume]) |                                                                                                      |
| [`config_auto_correlator_output`](#pylimer_tools_cpp.MEHPSimulator.config_auto_correlator_output)(self, values)       |                                                                                                      |
| [`config_average_output`](#pylimer_tools_cpp.MEHPSimulator.config_average_output)(self, values)               |                                                                                                      |
| [`config_mean_bond_length`](#pylimer_tools_cpp.MEHPSimulator.config_mean_bond_length)(self, mean_bond_length)   |                                                                                                      |
| [`config_mean_squared_bond_length`](#pylimer_tools_cpp.MEHPSimulator.config_mean_squared_bond_length)(self, r02)        |                                                                                                      |
| [`config_remove_primary_loops`](#pylimer_tools_cpp.MEHPSimulator.config_remove_primary_loops)(self[, remove])       |                                                                                                      |
| [`config_restart_output`](#pylimer_tools_cpp.MEHPSimulator.config_restart_output)(self, file[, output_every]) |                                                                                                      |
| [`config_step_output`](#pylimer_tools_cpp.MEHPSimulator.config_step_output)(self, output_configuration)    |                                                                                                      |
| [`construct_from_string`](#pylimer_tools_cpp.MEHPSimulator.construct_from_string)(s)                          |                                                                                                      |
| [`construct_with_random_entanglements`](#pylimer_tools_cpp.MEHPSimulator.construct_with_random_entanglements)(...[, ...])   | Construct MEHPSimulator from a Universe with randomly sampled entanglements.                         |
| [`count_active_clustered_atoms`](#pylimer_tools_cpp.MEHPSimulator.count_active_clustered_atoms)(self[, tolerance])   |                                                                                                      |
| [`deform_to`](#pylimer_tools_cpp.MEHPSimulator.deform_to)(self, target_box)                       |                                                                                                      |
| [`get_active_weight_fraction`](#pylimer_tools_cpp.MEHPSimulator.get_active_weight_fraction)(self[, tolerance])     |                                                                                                      |
| [`get_average_contour_length`](#pylimer_tools_cpp.MEHPSimulator.get_average_contour_length)(self)                  |                                                                                                      |
| [`get_average_spring_length`](#pylimer_tools_cpp.MEHPSimulator.get_average_spring_length)(self)                   |                                                                                                      |
| [`get_bond_lengths`](#pylimer_tools_cpp.MEHPSimulator.get_bond_lengths)(self)                            |                                                                                                      |
| [`get_coordinates`](#pylimer_tools_cpp.MEHPSimulator.get_coordinates)(self)                             |                                                                                                      |
| [`get_crosslinker_universe`](#pylimer_tools_cpp.MEHPSimulator.get_crosslinker_universe)(self)                    |                                                                                                      |
| [`get_crosslinker_verse`](#pylimer_tools_cpp.MEHPSimulator.get_crosslinker_verse)(self)                       |                                                                                                      |
| [`get_current_displacements`](#pylimer_tools_cpp.MEHPSimulator.get_current_displacements)(self)                   |                                                                                                      |
| [`get_current_time`](#pylimer_tools_cpp.MEHPSimulator.get_current_time)(self[, current_step])            |                                                                                                      |
| [`get_dangling_weight_fraction`](#pylimer_tools_cpp.MEHPSimulator.get_dangling_weight_fraction)(self[, tolerance])   |                                                                                                      |
| [`get_default_bond_length`](#pylimer_tools_cpp.MEHPSimulator.get_default_bond_length)(self)                     |                                                                                                      |
| [`get_default_r0_square`](#pylimer_tools_cpp.MEHPSimulator.get_default_r0_square)(self)                       |                                                                                                      |
| [`get_effective_functionality_of_atoms`](#pylimer_tools_cpp.MEHPSimulator.get_effective_functionality_of_atoms)(self[, ...]) |                                                                                                      |
| [`get_energy`](#pylimer_tools_cpp.MEHPSimulator.get_energy)(self)                                  |                                                                                                      |
| [`get_exit_reason`](#pylimer_tools_cpp.MEHPSimulator.get_exit_reason)(self)                             |                                                                                                      |
| [`get_force_evaluator`](#pylimer_tools_cpp.MEHPSimulator.get_force_evaluator)(self)                         |                                                                                                      |
| [`get_gamma`](#pylimer_tools_cpp.MEHPSimulator.get_gamma)(self)                                   |                                                                                                      |
| [`get_gamma_factor`](#pylimer_tools_cpp.MEHPSimulator.get_gamma_factor)(self[, b02, nr_of_chains])       | Compute the average topological Gamma factor `sum(gamma_factors) / denom`.                           |
| [`get_gamma_factors`](#pylimer_tools_cpp.MEHPSimulator.get_gamma_factors)(self[, b02])                    | Compute the per-spring topological Gamma factors $\Gamma_i = \mathbf{r}_i^2 / (N_i b_0^2)$.          |
| [`get_gamma_factors_in_dir`](#pylimer_tools_cpp.MEHPSimulator.get_gamma_factors_in_dir)(self, dir[, b02])        | Compute the per-spring directional Gamma factors $\Gamma_{i,\alpha} = r_{i,\alpha}^2 / (N_i b_0^2)$. |
| [`get_indices_of_active_nodes`](#pylimer_tools_cpp.MEHPSimulator.get_indices_of_active_nodes)(self[, tolerance])    |                                                                                                      |
| [`get_initial_coordinates`](#pylimer_tools_cpp.MEHPSimulator.get_initial_coordinates)(self)                     |                                                                                                      |
| [`get_is_2d`](#pylimer_tools_cpp.MEHPSimulator.get_is_2d)(self)                                   |                                                                                                      |
| [`get_last_solve_result`](#pylimer_tools_cpp.MEHPSimulator.get_last_solve_result)(self)                       |                                                                                                      |
| [`get_network`](#pylimer_tools_cpp.MEHPSimulator.get_network)(self)                                 |                                                                                                      |
| [`get_nr_of_active_nodes`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_active_nodes)(self[, tolerance])         |                                                                                                      |
| [`get_nr_of_active_springs`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_active_springs)(self[, tolerance])       |                                                                                                      |
| [`get_nr_of_active_strands`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_active_strands)(self[, tolerance])       |                                                                                                      |
| [`get_nr_of_atoms`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_atoms)(self)                             |                                                                                                      |
| [`get_nr_of_bonds`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_bonds)(self)                             |                                                                                                      |
| [`get_nr_of_extra_atoms`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_extra_atoms)(self)                       |                                                                                                      |
| [`get_nr_of_extra_bonds`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_extra_bonds)(self)                       |                                                                                                      |
| [`get_nr_of_intra_chain_sliplinks`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_intra_chain_sliplinks)(self)             |                                                                                                      |
| [`get_nr_of_iterations`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_iterations)(self)                        |                                                                                                      |
| [`get_nr_of_links`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_links)(self)                             |                                                                                                      |
| [`get_nr_of_nodes`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_nodes)(self)                             |                                                                                                      |
| [`get_nr_of_particles`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_particles)(self)                         |                                                                                                      |
| [`get_nr_of_springs`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_springs)(self)                           |                                                                                                      |
| [`get_nr_of_strands`](#pylimer_tools_cpp.MEHPSimulator.get_nr_of_strands)(self)                           |                                                                                                      |
| [`get_num_bonds_to_form`](#pylimer_tools_cpp.MEHPSimulator.get_num_bonds_to_form)(self)                       |                                                                                                      |
| [`get_num_intra_chain_slip_links`](#pylimer_tools_cpp.MEHPSimulator.get_num_intra_chain_slip_links)(self)              |                                                                                                      |
| [`get_num_intra_chain_sliplinks`](#pylimer_tools_cpp.MEHPSimulator.get_num_intra_chain_sliplinks)(self)               |                                                                                                      |
| [`get_num_relocations`](#pylimer_tools_cpp.MEHPSimulator.get_num_relocations)(self)                         |                                                                                                      |
| [`get_num_shifts`](#pylimer_tools_cpp.MEHPSimulator.get_num_shifts)(self)                              |                                                                                                      |
| [`get_pressure`](#pylimer_tools_cpp.MEHPSimulator.get_pressure)(self)                                |                                                                                                      |
| [`get_residual`](#pylimer_tools_cpp.MEHPSimulator.get_residual)(self)                                |                                                                                                      |
| [`get_residual_norm`](#pylimer_tools_cpp.MEHPSimulator.get_residual_norm)(self)                           |                                                                                                      |
| [`get_residuals`](#pylimer_tools_cpp.MEHPSimulator.get_residuals)(self)                               |                                                                                                      |
| [`get_simulation_has_run`](#pylimer_tools_cpp.MEHPSimulator.get_simulation_has_run)(self)                      |                                                                                                      |
| [`get_soluble_weight_fraction`](#pylimer_tools_cpp.MEHPSimulator.get_soluble_weight_fraction)(self[, tolerance])    |                                                                                                      |
| [`get_spring_contour_length`](#pylimer_tools_cpp.MEHPSimulator.get_spring_contour_length)(self)                   |                                                                                                      |
| [`get_spring_lengths`](#pylimer_tools_cpp.MEHPSimulator.get_spring_lengths)(self)                          |                                                                                                      |
| [`get_spring_partitions`](#pylimer_tools_cpp.MEHPSimulator.get_spring_partitions)(self)                       |                                                                                                      |
| [`get_spring_vectors`](#pylimer_tools_cpp.MEHPSimulator.get_spring_vectors)(self)                          |                                                                                                      |
| [`get_strand_contour_length`](#pylimer_tools_cpp.MEHPSimulator.get_strand_contour_length)(self)                   |                                                                                                      |
| [`get_stress_tensor`](#pylimer_tools_cpp.MEHPSimulator.get_stress_tensor)(self)                           |                                                                                                      |
| [`get_temperature`](#pylimer_tools_cpp.MEHPSimulator.get_temperature)(self)                             |                                                                                                      |
| [`get_timestep`](#pylimer_tools_cpp.MEHPSimulator.get_timestep)(self)                                |                                                                                                      |
| [`get_volume`](#pylimer_tools_cpp.MEHPSimulator.get_volume)(self)                                  |                                                                                                      |
| [`has_outputs`](#pylimer_tools_cpp.MEHPSimulator.has_outputs)(self)                                 |                                                                                                      |
| [`randomly_add_slip_links`](#pylimer_tools_cpp.MEHPSimulator.randomly_add_slip_links)(self, ...[, cutoff, ...]) |                                                                                                      |
| [`randomly_add_sliplinks`](#pylimer_tools_cpp.MEHPSimulator.randomly_add_sliplinks)(self, ...[, cutoff, ...])  |                                                                                                      |
| [`read_restart_file`](#pylimer_tools_cpp.MEHPSimulator.read_restart_file)(file)                           |                                                                                                      |
| [`remove_inactive_links`](#pylimer_tools_cpp.MEHPSimulator.remove_inactive_links)(self[, tolerance])          |                                                                                                      |
| [`run`](#pylimer_tools_cpp.MEHPSimulator.run)(\*args, \*\*kwargs)                           | Overloaded function.                                                                                 |
| [`serialize_to_string`](#pylimer_tools_cpp.MEHPSimulator.serialize_to_string)(self)                         |                                                                                                      |
| [`set_current_displacements`](#pylimer_tools_cpp.MEHPSimulator.set_current_displacements)(self, ...)              |                                                                                                      |
| [`set_force_evaluator`](#pylimer_tools_cpp.MEHPSimulator.set_force_evaluator)(self, force_evaluator)        |                                                                                                      |
| [`set_is_2d`](#pylimer_tools_cpp.MEHPSimulator.set_is_2d)(self, is_2d)                            |                                                                                                      |
| [`simplify_network`](#pylimer_tools_cpp.MEHPSimulator.simplify_network)(self, mode[, tolerance, ...])    |                                                                                                      |
| [`unlink_bifunctional_links`](#pylimer_tools_cpp.MEHPSimulator.unlink_bifunctional_links)(self[, ...])            |                                                                                                      |
| [`validate_network`](#pylimer_tools_cpp.MEHPSimulator.validate_network)(self)                            |                                                                                                      |
| [`write_restart_file`](#pylimer_tools_cpp.MEHPSimulator.write_restart_file)(self, file)                    |                                                                                                      |

### Attributes Documentation

#### force_evaluator

#### is_2d

#### last_solve_result

#### network

### Methods Documentation

#### add_slip_links(\*args, \*\*kwargs)

Overloaded function.

1. add_slip_links(self: pylimer_tools_cpp.MEHPSimulator, entanglements: pylimer_tools_cpp.StrandPairEntanglements, clamp_alpha: bool = False) -> int
2. add_slip_links(self: pylimer_tools_cpp.MEHPSimulator, strand_idx_1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_idx_2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], y: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], z: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], alpha_1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], alpha_2: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], clamp_alpha: bool = False) -> int
3. add_slip_links(self: pylimer_tools_cpp.MEHPSimulator, strand_idx_1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_idx_2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], y: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], z: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], clamp_alpha: bool = False) -> int

#### add_sliplinks(\*args, \*\*kwargs)

Overloaded function.

1. add_sliplinks(self: pylimer_tools_cpp.MEHPSimulator, entanglements: pylimer_tools_cpp.StrandPairEntanglements, clamp_alpha: bool = False) -> int
2. add_sliplinks(self: pylimer_tools_cpp.MEHPSimulator, strand_idx_1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_idx_2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], y: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], z: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], alpha_1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], alpha_2: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], clamp_alpha: bool = False) -> int
3. add_sliplinks(self: pylimer_tools_cpp.MEHPSimulator, strand_idx_1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_idx_2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], y: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], z: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, “[m, 1]”], clamp_alpha: bool = False) -> int

#### break_too_long_strands(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), breaking_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### config_assume_box_large_enough(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), assume: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### config_assume_network_is_complete(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), assume: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### config_auto_correlator_output(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), values: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.OutputConfiguration](pylimer_tools_cpp.OutputConfiguration.md#pylimer_tools_cpp.OutputConfiguration)], num_corr_in: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 32, p: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 16, m: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### config_average_output(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), values: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.OutputConfiguration](pylimer_tools_cpp.OutputConfiguration.md#pylimer_tools_cpp.OutputConfiguration)]) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### config_mean_bond_length(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), mean_bond_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### config_mean_squared_bond_length(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), r02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### config_remove_primary_loops(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), remove: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### config_restart_output(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), file: [str](https://docs.python.org/3/builtins/stdtypes.html#str), output_every: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 50000) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### config_step_output(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), output_configuration: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.OutputConfiguration](pylimer_tools_cpp.OutputConfiguration.md#pylimer_tools_cpp.OutputConfiguration)]) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### *static* construct_from_string(s: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)

#### *static* construct_with_random_entanglements(universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), nr_of_entanglements_to_sample: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), upper_sampling_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), lower_sampling_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.0, min_nr_of_entanglements: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, same_strand_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.0, seed: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = '', crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2, entanglement_mode: [pylimer_tools_cpp.EntanglementMode](pylimer_tools_cpp.EntanglementMode.md#pylimer_tools_cpp.EntanglementMode) = EntanglementMode.MERGED_LINKS, with_crosslinkers: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, filter_dangling_and_soluble: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, force_evaluator: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator) = None, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, keep_2functional_crosslinkers: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, remove_primary_loops: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, assume_box_large_enough: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)

Construct MEHPSimulator from a Universe with randomly sampled entanglements.

#### count_active_clustered_atoms(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### deform_to(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), target_box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### get_active_weight_fraction(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_average_contour_length(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_average_spring_length(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_bond_lengths(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_coordinates(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_crosslinker_universe(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

#### get_crosslinker_verse(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

#### get_current_displacements(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_current_time(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), current_step: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.0) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_dangling_weight_fraction(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_default_bond_length(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_default_r0_square(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_effective_functionality_of_atoms(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [dict](https://docs.python.org/3/builtins/stdtypes.html#dict)[[int](https://docs.python.org/3/builtins/functions.html#int), [int](https://docs.python.org/3/builtins/functions.html#int)]

#### get_energy(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_exit_reason(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [pylimer_tools_cpp.ExitReason](pylimer_tools_cpp.ExitReason.md#pylimer_tools_cpp.ExitReason)

#### get_force_evaluator(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator)

#### get_gamma(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_gamma_factor(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0, nr_of_chains: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the average topological Gamma factor `sum(gamma_factors) / denom`.

* **Parameters:**
  * **b02** – Mean-squared bond length $b_0^2$ (defaults to network bond length squared).
  * **nr_of_chains** – Denominator for averaging. When `< 1` (default `-1`), divides by
    `network.nr_of_strands` — note that in `EntanglementMode.SPRING`,
    `nr_of_strands` includes entanglement coupling springs, so comparing
    `MERGED_LINKS` vs. `SPRING` with the default `-1` is not apples-to-apples.
    Pass `1` or use `np.sum(get_gamma_factors(b02))` when computing shear modulus.

#### get_gamma_factors(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute the per-spring topological Gamma factors $\Gamma_i = \mathbf{r}_i^2 / (N_i b_0^2)$.

Summing this vector (`np.sum(sim.get_gamma_factors(b02))`) gives the unnormalized total
needed to compute the shear modulus $G = (k_B T / V) \sum_i \Gamma_i$.

* **Parameters:**
  **b02** – Mean-squared bond length $b_0^2$ (defaults to network bond length squared).

#### get_gamma_factors_in_dir(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), dir: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute the per-spring directional Gamma factors $\Gamma_{i,\alpha} = r_{i,\alpha}^2 / (N_i b_0^2)$.

* **Parameters:**
  * **dir** – Cartesian direction index (`0`, `1`, or `2`).
  * **b02** – Mean-squared bond length $b_0^2$ (defaults to network bond length squared).

#### get_indices_of_active_nodes(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[int](https://docs.python.org/3/builtins/functions.html#int)]

#### get_initial_coordinates(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_is_2d(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### get_last_solve_result(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [pylimer_tools_cpp.MEHPSolveResult](pylimer_tools_cpp.MEHPSolveResult.md#pylimer_tools_cpp.MEHPSolveResult)

#### get_network(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [pylimer_tools_cpp.MEHPNetwork](pylimer_tools_cpp.MEHPNetwork.md#pylimer_tools_cpp.MEHPNetwork)

#### get_nr_of_active_nodes(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_active_springs(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_active_strands(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_atoms(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_bonds(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_extra_atoms(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_extra_bonds(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_intra_chain_sliplinks(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_iterations(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_links(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_nodes(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_particles(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_springs(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_strands(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_num_bonds_to_form(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_num_intra_chain_slip_links(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_num_intra_chain_sliplinks(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_num_relocations(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_num_shifts(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_pressure(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_residual(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_residual_norm(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_residuals(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_simulation_has_run(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### get_soluble_weight_fraction(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_spring_contour_length(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_spring_lengths(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_spring_partitions(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_spring_vectors(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_strand_contour_length(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### get_stress_tensor(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

#### get_temperature(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_timestep(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_volume(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### has_outputs(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### randomly_add_slip_links(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), nr_of_entanglements_to_sample: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2.0, min_nr_of_entanglements: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, same_strand_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2.0, seed: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = '', clamp_alpha: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### randomly_add_sliplinks(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), nr_of_entanglements_to_sample: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2.0, min_nr_of_entanglements: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, same_strand_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2.0, seed: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = '', clamp_alpha: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### *static* read_restart_file(file: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)

#### remove_inactive_links(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-06) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### run(\*args, \*\*kwargs)

Overloaded function.

1. run(self: pylimer_tools_cpp.MEHPSimulator) -> pylimer_tools_cpp.MEHPSolveResult
2. run(self: pylimer_tools_cpp.MEHPSimulator, solver: pylimer_tools_cpp.MEHPSolver) -> pylimer_tools_cpp.MEHPSolveResult
3. run(self: pylimer_tools_cpp.MEHPSimulator, config: pylimer_tools_cpp.MEHPLinearSolverConfig) -> pylimer_tools_cpp.MEHPSolveResult
4. run(self: pylimer_tools_cpp.MEHPSimulator, config: pylimer_tools_cpp.MEHPSlipLinkSolverConfig) -> pylimer_tools_cpp.MEHPSolveResult
5. run(self: pylimer_tools_cpp.MEHPSimulator, config: pylimer_tools_cpp.MEHPNonlinearSolverConfig) -> pylimer_tools_cpp.MEHPSolveResult

#### serialize_to_string(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [bytes](https://docs.python.org/3/builtins/stdtypes.html#bytes)

#### set_current_displacements(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), new_displacements: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_force_evaluator(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), force_evaluator: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_is_2d(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### simplify_network(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), mode: [pylimer_tools_cpp.StructureSimplificationMode](pylimer_tools_cpp.StructureSimplificationMode.md#pylimer_tools_cpp.StructureSimplificationMode), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-06, breaking_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### unlink_bifunctional_links(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), entanglements_only: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### validate_network(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### write_restart_file(self: [pylimer_tools_cpp.MEHPSimulator](#pylimer_tools_cpp.MEHPSimulator), file: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [None](https://docs.python.org/3/builtins/constants.html#None)
