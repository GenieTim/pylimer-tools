# RouseBrownianSimulator

### *class* pylimer_tools_cpp.RouseBrownianSimulator(\*args, \*\*kwargs)

Bases: `pybind11_object`

Overloaded function.

1. \_\_init_\_(self: pylimer_tools_cpp.RouseBrownianSimulator, spring_from: collections.abc.Sequence[typing.SupportsInt], spring_to: collections.abc.Sequence[typing.SupportsInt], gamma: typing.SupportsFloat = 1.0, k: typing.SupportsFloat = 1.0, T: typing.SupportsFloat = 1.0, dt: typing.SupportsFloat = 0.01, num_nodes: typing.SupportsInt = 0) -> None
2. \_\_init_\_(self: pylimer_tools_cpp.RouseBrownianSimulator, L: scipy.sparse.csc_matrix[numpy.float64], gamma: typing.SupportsFloat = 1.0, k: typing.SupportsFloat = 1.0, T: typing.SupportsFloat = 1.0, dt: typing.SupportsFloat = 0.01) -> None

### Methods Summary

| [`compute_center_of_mass_msd`](#pylimer_tools_cpp.RouseBrownianSimulator.compute_center_of_mass_msd)(trajectory[, ...])     |    |
|----------------------------------------------------------------------------------------------------|----|
| [`compute_msd`](#pylimer_tools_cpp.RouseBrownianSimulator.compute_msd)(trajectory[, time_averaged])          |    |
| [`create_linear_chain`](#pylimer_tools_cpp.RouseBrownianSimulator.create_linear_chain)(num_nodes[, gamma, k, T, dt]) |    |
| [`create_ring`](#pylimer_tools_cpp.RouseBrownianSimulator.create_ring)(num_nodes[, gamma, k, T, dt])         |    |
| [`create_star`](#pylimer_tools_cpp.RouseBrownianSimulator.create_star)(num_arms[, arm_length, gamma, ...])   |    |
| [`get_center_of_mass`](#pylimer_tools_cpp.RouseBrownianSimulator.get_center_of_mass)(self)                          |    |
| [`get_diffusion_coefficient`](#pylimer_tools_cpp.RouseBrownianSimulator.get_diffusion_coefficient)(self)                   |    |
| [`get_dt`](#pylimer_tools_cpp.RouseBrownianSimulator.get_dt)(self)                                      |    |
| [`get_end_to_end_distance`](#pylimer_tools_cpp.RouseBrownianSimulator.get_end_to_end_distance)(self[, node_a, node_b])   |    |
| [`get_end_to_end_vector`](#pylimer_tools_cpp.RouseBrownianSimulator.get_end_to_end_vector)(self[, node_a, node_b])     |    |
| [`get_gamma`](#pylimer_tools_cpp.RouseBrownianSimulator.get_gamma)(self)                                   |    |
| [`get_k`](#pylimer_tools_cpp.RouseBrownianSimulator.get_k)(self)                                       |    |
| [`get_laplacian`](#pylimer_tools_cpp.RouseBrownianSimulator.get_laplacian)(self)                               |    |
| [`get_nr_of_nodes`](#pylimer_tools_cpp.RouseBrownianSimulator.get_nr_of_nodes)(self)                             |    |
| [`get_positions`](#pylimer_tools_cpp.RouseBrownianSimulator.get_positions)(self)                               |    |
| [`get_potential_energy`](#pylimer_tools_cpp.RouseBrownianSimulator.get_potential_energy)(self)                        |    |
| [`get_radius_of_gyration`](#pylimer_tools_cpp.RouseBrownianSimulator.get_radius_of_gyration)(self)                      |    |
| [`get_shifted_matrix`](#pylimer_tools_cpp.RouseBrownianSimulator.get_shifted_matrix)(self)                          |    |
| [`get_temperature`](#pylimer_tools_cpp.RouseBrownianSimulator.get_temperature)(self)                             |    |
| [`initialize_linear_positions`](#pylimer_tools_cpp.RouseBrownianSimulator.initialize_linear_positions)(self[, ...])          |    |
| [`initialize_random_walk`](#pylimer_tools_cpp.RouseBrownianSimulator.initialize_random_walk)(self[, bond_length, seed]) |    |
| [`run`](#pylimer_tools_cpp.RouseBrownianSimulator.run)(self, total_steps[, sample_interval, ...])    |    |
| [`set_dt`](#pylimer_tools_cpp.RouseBrownianSimulator.set_dt)(self, arg0)                                |    |
| [`set_gamma`](#pylimer_tools_cpp.RouseBrownianSimulator.set_gamma)(self, arg0)                             |    |
| [`set_k`](#pylimer_tools_cpp.RouseBrownianSimulator.set_k)(self, arg0)                                 |    |
| [`set_parameters`](#pylimer_tools_cpp.RouseBrownianSimulator.set_parameters)(self, arg0, arg1, arg2, arg3)      |    |
| [`set_positions`](#pylimer_tools_cpp.RouseBrownianSimulator.set_positions)(self, arg0)                         |    |
| [`set_seed`](#pylimer_tools_cpp.RouseBrownianSimulator.set_seed)(self, arg0)                              |    |
| [`set_solver_type`](#pylimer_tools_cpp.RouseBrownianSimulator.set_solver_type)(self, arg0)                       |    |
| [`set_temperature`](#pylimer_tools_cpp.RouseBrownianSimulator.set_temperature)(self, arg0)                       |    |
| [`step`](#pylimer_tools_cpp.RouseBrownianSimulator.step)(self[, num_steps])                           |    |
| [`step_deterministic`](#pylimer_tools_cpp.RouseBrownianSimulator.step_deterministic)(self[, num_steps])             |    |
| [`step_with_external_forces`](#pylimer_tools_cpp.RouseBrownianSimulator.step_with_external_forces)(self, forces[, ...])    |    |

### Methods Documentation

#### *static* compute_center_of_mass_msd(trajectory: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']], time_averaged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### *static* compute_msd(trajectory: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']], time_averaged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### *static* create_linear_chain(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

#### *static* create_ring(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

#### *static* create_star(num_arms: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), arm_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1, gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

#### get_center_of_mass(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

#### get_diffusion_coefficient(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_dt(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_end_to_end_distance(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), node_a: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0, node_b: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 18446744073709551615) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_end_to_end_vector(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), node_a: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0, node_b: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 18446744073709551615) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

#### get_gamma(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_k(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_laplacian(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]

#### get_nr_of_nodes(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_positions(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

#### get_potential_energy(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_radius_of_gyration(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### get_shifted_matrix(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]

#### get_temperature(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### initialize_linear_positions(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), bond_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, axis: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### initialize_random_walk(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), bond_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, seed: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 42) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### run(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), total_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), sample_interval: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1, save_trajectory: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, end_to_end_node_a: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0, end_to_end_node_b: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 18446744073709551615) → [pylimer_tools_cpp.RouseSimulationResult](pylimer_tools_cpp.RouseSimulationResult.md#pylimer_tools_cpp.RouseSimulationResult)

#### set_dt(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_gamma(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_k(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_parameters(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), arg1: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), arg2: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), arg3: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_positions(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), arg0: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_seed(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), arg0: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_solver_type(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), arg0: [pylimer_tools_cpp.RouseSolverType](pylimer_tools_cpp.RouseSolverType.md#pylimer_tools_cpp.RouseSolverType)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### set_temperature(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), arg0: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### step(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### step_deterministic(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### step_with_external_forces(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), forces: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]'], num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1) → [None](https://docs.python.org/3/builtins/constants.html#None)
