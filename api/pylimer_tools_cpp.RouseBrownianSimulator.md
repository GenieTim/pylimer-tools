# RouseBrownianSimulator

### *class* pylimer_tools_cpp.RouseBrownianSimulator(\*args, \*\*kwargs)

Bases: `pybind11_object`

Implicit time-stepping Rouse and Brownian dynamics simulator for polymer chains
and general spring networks.

Overloaded function.

1. \_\_init_\_(self: pylimer_tools_cpp.RouseBrownianSimulator, spring_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], spring_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01, num_nodes: typing.SupportsInt | typing.SupportsIndex = 0, weights: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] = []) -> None
   > Construct from edge lists (spring_from, spring_to), friction gamma,
   > spring constant k, temperature T, time step dt, optional num_nodes,
   > and optional edge weights.
   > * **param spring_from:**
   >   Source node indices of springs
   > * **param spring_to:**
   >   Target node indices of springs
   > * **param gamma:**
   >   Bead friction coefficient
   > * **param k:**
   >   Hookean spring constant
   > * **param T:**
   >   Temperature
   > * **param dt:**
   >   Integration time step
   > * **param num_nodes:**
   >   Total number of nodes (deduced if 0)
   > * **param weights:**
   >   Optional positive edge weights (defaults to 1.0)
2. \_\_init_\_(self: pylimer_tools_cpp.RouseBrownianSimulator, L: scipy.sparse.csc_matrix[numpy.float64], gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> None
   > Construct from a pre-assembled sparse scalar graph Laplacian matrix L.
   > * **param L:**
   >   Symmetric positive semidefinite graph Laplacian matrix
   > * **param gamma:**
   >   Bead friction coefficient
   > * **param k:**
   >   Hookean spring constant
   > * **param T:**
   >   Temperature
   > * **param dt:**
   >   Integration time step
3. \_\_init_\_(self: pylimer_tools_cpp.RouseBrownianSimulator, universe: pylimer_tools_cpp.Universe, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> None
   > Construct from a Universe’s atom-bond graph and initial coordinates.
   > * **param universe:**
   >   Universe instance
   > * **param gamma:**
   >   Bead friction coefficient
   > * **param k:**
   >   Hookean spring constant
   > * **param T:**
   >   Temperature
   > * **param dt:**
   >   Integration time step

### Methods Summary

| [`compute_center_of_mass_msd`](#pylimer_tools_cpp.RouseBrownianSimulator.compute_center_of_mass_msd)(trajectory[, ...])     | Compute center-of-mass mean-squared displacement from a trajectory of (N, 3) position matrices.                                             |
|----------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| [`compute_end_to_end_autocorrelation`](#pylimer_tools_cpp.RouseBrownianSimulator.compute_end_to_end_autocorrelation)(trajectory)    | Compute normalized end-to-end vector autocorrelation from a trajectory of (N, 3) position matrices (node_b=None defaults to N - 1).         |
| [`compute_msd`](#pylimer_tools_cpp.RouseBrownianSimulator.compute_msd)(trajectory[, time_averaged])          | Compute bead mean-squared displacement from a trajectory of (N, 3) position matrices.                                                       |
| [`create_complete`](#pylimer_tools_cpp.RouseBrownianSimulator.create_complete)(num_nodes[, gamma, k, T, dt])     | Create a RouseBrownianSimulator for a complete graph K_N of num_nodes beads.                                                                |
| [`create_linear_chain`](#pylimer_tools_cpp.RouseBrownianSimulator.create_linear_chain)(num_nodes[, gamma, k, T, dt]) | Create a RouseBrownianSimulator for a linear chain of num_nodes beads.                                                                      |
| [`create_ring`](#pylimer_tools_cpp.RouseBrownianSimulator.create_ring)(num_nodes[, gamma, k, T, dt])         | Create a RouseBrownianSimulator for a ring polymer of num_nodes beads.                                                                      |
| [`create_star`](#pylimer_tools_cpp.RouseBrownianSimulator.create_star)(num_arms[, arm_length, gamma, ...])   | Create a RouseBrownianSimulator for a star polymer with num_arms arms of arm_length bonds each.                                             |
| [`from_matrix`](#pylimer_tools_cpp.RouseBrownianSimulator.from_matrix)(L[, gamma, k, T, dt])                 | Create a RouseBrownianSimulator from a sparse graph Laplacian matrix L.                                                                     |
| [`from_universe`](#pylimer_tools_cpp.RouseBrownianSimulator.from_universe)(universe[, gamma, k, T, dt])        | Create a RouseBrownianSimulator from a Universe's atom-bond graph and coordinates.                                                          |
| [`get_center_of_mass`](#pylimer_tools_cpp.RouseBrownianSimulator.get_center_of_mass)(self)                          | Compute the 3D center of mass of all nodes.                                                                                                 |
| [`get_diffusion_coefficient`](#pylimer_tools_cpp.RouseBrownianSimulator.get_diffusion_coefficient)(self)                   | Get the bead diffusion coefficient D = kB \* T / gamma.                                                                                     |
| [`get_dt`](#pylimer_tools_cpp.RouseBrownianSimulator.get_dt)(self)                                      | Get the integration time step dt.                                                                                                           |
| [`get_end_to_end_distance`](#pylimer_tools_cpp.RouseBrownianSimulator.get_end_to_end_distance)(self[, node_a, node_b])   | Compute the Euclidean distance between node_a and node_b (None defaults to N - 1).                                                          |
| [`get_end_to_end_vector`](#pylimer_tools_cpp.RouseBrownianSimulator.get_end_to_end_vector)(self[, node_a, node_b])     | Compute the 3D vector from node_a to node_b (None defaults to N - 1).                                                                       |
| [`get_gamma`](#pylimer_tools_cpp.RouseBrownianSimulator.get_gamma)(self)                                   | Get the friction coefficient gamma.                                                                                                         |
| [`get_k`](#pylimer_tools_cpp.RouseBrownianSimulator.get_k)(self)                                       | Get the Hookean spring constant k.                                                                                                          |
| [`get_kb`](#pylimer_tools_cpp.RouseBrownianSimulator.get_kb)(self)                                      | Get the Boltzmann constant kB.                                                                                                              |
| [`get_laplacian`](#pylimer_tools_cpp.RouseBrownianSimulator.get_laplacian)(self)                               | Get the sparse graph Laplacian matrix L.                                                                                                    |
| [`get_nr_of_nodes`](#pylimer_tools_cpp.RouseBrownianSimulator.get_nr_of_nodes)(self)                             | Get the number of nodes (beads) in the simulator.                                                                                           |
| [`get_positions`](#pylimer_tools_cpp.RouseBrownianSimulator.get_positions)(self)                               | Get the (N, 3) bead positions matrix.                                                                                                       |
| [`get_positions_flattened`](#pylimer_tools_cpp.RouseBrownianSimulator.get_positions_flattened)(self[, interleaved])      | Get bead positions as a flattened 3N vector.                                                                                                |
| [`get_potential_energy`](#pylimer_tools_cpp.RouseBrownianSimulator.get_potential_energy)(self)                        | Compute the harmonic spring potential energy of the current configuration.                                                                  |
| [`get_radius_of_gyration`](#pylimer_tools_cpp.RouseBrownianSimulator.get_radius_of_gyration)(self)                      | Compute the radius of gyration of the current bead configuration.                                                                           |
| [`get_shifted_matrix`](#pylimer_tools_cpp.RouseBrownianSimulator.get_shifted_matrix)(self)                          | Get the shifted implicit stepping matrix A = I + dt \* (k / gamma) \* L.                                                                    |
| [`get_solver_type`](#pylimer_tools_cpp.RouseBrownianSimulator.get_solver_type)(self)                             | Get the active linear solver type.                                                                                                          |
| [`get_temperature`](#pylimer_tools_cpp.RouseBrownianSimulator.get_temperature)(self)                             | Get the temperature T.                                                                                                                      |
| [`initialize_linear_positions`](#pylimer_tools_cpp.RouseBrownianSimulator.initialize_linear_positions)(self[, ...])          | Initialize bead positions spaced along a coordinate axis.                                                                                   |
| [`initialize_random_walk`](#pylimer_tools_cpp.RouseBrownianSimulator.initialize_random_walk)(self[, bond_length, seed]) | Initialize bead positions as a 3D random walk.                                                                                              |
| [`run`](#pylimer_tools_cpp.RouseBrownianSimulator.run)(self, total_steps[, sample_interval, ...])    | Run the simulation for total_steps and record dynamical observables every sample_interval steps (end_to_end_node_b=None defaults to N - 1). |
| [`set_dt`](#pylimer_tools_cpp.RouseBrownianSimulator.set_dt)(self, dt)                                  | Set the integration time step dt.                                                                                                           |
| [`set_gamma`](#pylimer_tools_cpp.RouseBrownianSimulator.set_gamma)(self, gamma)                            | Set the friction coefficient gamma.                                                                                                         |
| [`set_k`](#pylimer_tools_cpp.RouseBrownianSimulator.set_k)(self, k)                                    | Set the Hookean spring constant k.                                                                                                          |
| [`set_kb`](#pylimer_tools_cpp.RouseBrownianSimulator.set_kb)(self, kb)                                  | Set the Boltzmann constant kB.                                                                                                              |
| [`set_parameters`](#pylimer_tools_cpp.RouseBrownianSimulator.set_parameters)(self, gamma, k, T, dt)             | Set gamma, k, T, and dt simultaneously.                                                                                                     |
| [`set_positions`](#pylimer_tools_cpp.RouseBrownianSimulator.set_positions)(self, positions)                    | Set the (N, 3) bead positions matrix.                                                                                                       |
| [`set_positions_flattened`](#pylimer_tools_cpp.RouseBrownianSimulator.set_positions_flattened)(self, positions[, ...])   | Set bead positions from a flattened 3N vector.                                                                                              |
| [`set_seed`](#pylimer_tools_cpp.RouseBrownianSimulator.set_seed)(self, seed)                              | Set the random number generator seed for thermal noise.                                                                                     |
| [`set_solver_type`](#pylimer_tools_cpp.RouseBrownianSimulator.set_solver_type)(self, solver_type)                | Set the linear solver algorithm used for implicit time-stepping.                                                                            |
| [`set_temperature`](#pylimer_tools_cpp.RouseBrownianSimulator.set_temperature)(self, T)                          | Set the temperature T.                                                                                                                      |
| [`step`](#pylimer_tools_cpp.RouseBrownianSimulator.step)(self[, num_steps])                           | Advance the simulation by num_steps implicit Brownian steps.                                                                                |
| [`step_deterministic`](#pylimer_tools_cpp.RouseBrownianSimulator.step_deterministic)(self[, num_steps])             | Advance the simulation by num_steps deterministic relaxation steps (without thermal noise).                                                 |
| [`step_deterministic_with_external_forces`](#pylimer_tools_cpp.RouseBrownianSimulator.step_deterministic_with_external_forces)(...)      | Advance the simulation by num_steps deterministic steps under external forces.                                                              |
| [`step_with_external_forces`](#pylimer_tools_cpp.RouseBrownianSimulator.step_with_external_forces)(self, forces[, ...])    | Advance the simulation by num_steps under external forces and thermal noise.                                                                |

### Methods Documentation

#### *static* compute_center_of_mass_msd(trajectory: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']], time_averaged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute center-of-mass mean-squared displacement from a trajectory of (N, 3) position matrices.

#### *static* compute_end_to_end_autocorrelation(trajectory: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']], node_a: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, node_b: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) | [None](https://docs.python.org/3/builtins/constants.html#None) = None, time_averaged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute normalized end-to-end vector autocorrelation from a trajectory of (N, 3) position matrices (node_b=None defaults to N - 1).

#### *static* compute_msd(trajectory: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']], time_averaged: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute bead mean-squared displacement from a trajectory of (N, 3) position matrices.

#### *static* create_complete(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

Create a RouseBrownianSimulator for a complete graph K_N of num_nodes beads.

* **Parameters:**
  * **num_nodes** – Number of nodes in the complete graph (at least 2)
  * **gamma** – Bead friction coefficient
  * **k** – Hookean spring constant
  * **T** – Temperature
  * **dt** – Integration time step
* **Returns:**
  RouseBrownianSimulator instance

#### *static* create_linear_chain(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

Create a RouseBrownianSimulator for a linear chain of num_nodes beads.

#### *static* create_ring(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

Create a RouseBrownianSimulator for a ring polymer of num_nodes beads.

#### *static* create_star(num_arms: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arm_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1, gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

Create a RouseBrownianSimulator for a star polymer with num_arms arms of arm_length bonds each.

#### *static* from_matrix(L: [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64], gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

Create a RouseBrownianSimulator from a sparse graph Laplacian matrix L.

* **Parameters:**
  * **L** – Symmetric positive semidefinite graph Laplacian matrix
  * **gamma** – Bead friction coefficient
  * **k** – Hookean spring constant
  * **T** – Temperature
  * **dt** – Integration time step
* **Returns:**
  RouseBrownianSimulator instance

#### *static* from_universe(universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.01) → [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)

Create a RouseBrownianSimulator from a Universe’s atom-bond graph and coordinates.

* **Parameters:**
  * **universe** – Universe instance
  * **gamma** – Bead friction coefficient
  * **k** – Hookean spring constant
  * **T** – Temperature
  * **dt** – Integration time step
* **Returns:**
  RouseBrownianSimulator instance

#### get_center_of_mass(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the 3D center of mass of all nodes.

#### get_diffusion_coefficient(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the bead diffusion coefficient D = kB \* T / gamma.

#### get_dt(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the integration time step dt.

#### get_end_to_end_distance(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), node_a: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, node_b: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) | [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the Euclidean distance between node_a and node_b (None defaults to N - 1).

#### get_end_to_end_vector(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), node_a: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, node_b: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) | [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the 3D vector from node_a to node_b (None defaults to N - 1).

#### get_gamma(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the friction coefficient gamma.

#### get_k(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the Hookean spring constant k.

#### get_kb(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the Boltzmann constant kB.

#### get_laplacian(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]

Get the sparse graph Laplacian matrix L.

#### get_nr_of_nodes(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of nodes (beads) in the simulator.

#### get_positions(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

Get the (N, 3) bead positions matrix.

#### get_positions_flattened(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), interleaved: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get bead positions as a flattened 3N vector.

#### get_potential_energy(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the harmonic spring potential energy of the current configuration.

#### get_radius_of_gyration(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the radius of gyration of the current bead configuration.

#### get_shifted_matrix(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]

Get the shifted implicit stepping matrix A = I + dt \* (k / gamma) \* L.

#### get_solver_type(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [pylimer_tools_cpp.RouseSolverType](pylimer_tools_cpp.RouseSolverType.md#pylimer_tools_cpp.RouseSolverType)

Get the active linear solver type.

#### get_temperature(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the temperature T.

#### initialize_linear_positions(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), bond_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, axis: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0) → [None](https://docs.python.org/3/builtins/constants.html#None)

Initialize bead positions spaced along a coordinate axis.

#### initialize_random_walk(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), bond_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, seed: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 42) → [None](https://docs.python.org/3/builtins/constants.html#None)

Initialize bead positions as a 3D random walk.

#### run(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), total_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), sample_interval: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1, save_trajectory: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, end_to_end_node_a: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, end_to_end_node_b: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) | [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [pylimer_tools_cpp.RouseSimulationResult](pylimer_tools_cpp.RouseSimulationResult.md#pylimer_tools_cpp.RouseSimulationResult)

Run the simulation for total_steps and record dynamical observables every sample_interval steps (end_to_end_node_b=None defaults to N - 1).

#### set_dt(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the integration time step dt.

#### set_gamma(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the friction coefficient gamma.

#### set_k(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the Hookean spring constant k.

#### set_kb(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), kb: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the Boltzmann constant kB.

#### set_parameters(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set gamma, k, T, and dt simultaneously.

#### set_positions(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), positions: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the (N, 3) bead positions matrix.

#### set_positions_flattened(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), positions: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], interleaved: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set bead positions from a flattened 3N vector.

#### set_seed(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), seed: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the random number generator seed for thermal noise.

#### set_solver_type(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), solver_type: [pylimer_tools_cpp.RouseSolverType](pylimer_tools_cpp.RouseSolverType.md#pylimer_tools_cpp.RouseSolverType)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the linear solver algorithm used for implicit time-stepping.

#### set_temperature(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), T: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the temperature T.

#### step(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1) → [None](https://docs.python.org/3/builtins/constants.html#None)

Advance the simulation by num_steps implicit Brownian steps.

#### step_deterministic(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1) → [None](https://docs.python.org/3/builtins/constants.html#None)

Advance the simulation by num_steps deterministic relaxation steps (without thermal noise).

#### step_deterministic_with_external_forces(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), forces: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]'], num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1) → [None](https://docs.python.org/3/builtins/constants.html#None)

Advance the simulation by num_steps deterministic steps under external forces.

#### step_with_external_forces(self: [pylimer_tools_cpp.RouseBrownianSimulator](#pylimer_tools_cpp.RouseBrownianSimulator), forces: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]'], num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1) → [None](https://docs.python.org/3/builtins/constants.html#None)

Advance the simulation by num_steps under external forces and thermal noise.
