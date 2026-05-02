# MEHPForceRelaxation

### *class* pylimer_tools_cpp.MEHPForceRelaxation(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, is_2d: [bool](https://docs.python.org/3/library/functions.html#bool) = False, force_evaluator: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator) = None, kappa: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, remove_2functional_crosslinkers: [bool](https://docs.python.org/3/library/functions.html#bool) = False, remove_dangling_chains: [bool](https://docs.python.org/3/library/functions.html#bool) = False)

Bases: `pybind11_object`

A small simulation tool for quickly minimizing the force between the crosslinker beads.

This is the first of three force relaxation methods available in this library.
The relevant feature of this implementation is the configurable spring potential.
Consequently, it offers a variety of configurable non-linear solvers using NLoptLib.

Please cite Gusev [[Gus19](../acknowledgements.md#id8)] if you use this method in your work.

Instantiate the simulator for a certain universe.

* **Parameters:**
  * **universe** – The universe to simulate with
  * **crosslinker_type** – The atom type of the crosslinkers. Needed to reduce the network.
  * **is2d** – Whether to ignore the z direction.
  * **force_evaluator** – The force evaluator to use
  * **kappa** – The spring constant
  * **remove_2functional_crosslinkers** – Whether to replace two-functional crosslinkers with a “normal” chain bead
  * **remove_dangling_chains** – Whether to remove dangling chains before running the simulation.
    **Caution**: Removing the dangling chains will result in incorrect results fo the computation of
    [`get_soluble_weight_fraction()`](#pylimer_tools_cpp.MEHPForceRelaxation.get_soluble_weight_fraction) and
    [`get_dangling_weight_fraction()`](#pylimer_tools_cpp.MEHPForceRelaxation.get_dangling_weight_fraction)

### Attributes Summary

| [`network`](#pylimer_tools_cpp.MEHPForceRelaxation.network)   | Get the network structure.   |
|---------------------------------------------------------------|------------------------------|

### Methods Summary

| [`assume_box_large_enough`](#pylimer_tools_cpp.MEHPForceRelaxation.assume_box_large_enough)(self[, box_large_enough])              | Configure whether to run PBC on the bonds or not.                                                                                                        |
|------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`config_rerun_epsilon`](#pylimer_tools_cpp.MEHPForceRelaxation.config_rerun_epsilon)(self[, epsilon])                             | Configure the offset from the lower and upper bounds for the simulation to suggest another run.                                                          |
| [`config_step_output`](#pylimer_tools_cpp.MEHPForceRelaxation.config_step_output)(\*args, \*\*kwargs)                              | Overloaded function.                                                                                                                                     |
| [`count_active_clustered_atoms`](#pylimer_tools_cpp.MEHPForceRelaxation.count_active_clustered_atoms)(self[, tolerance])           | Counts the active clustered atoms in the system.                                                                                                         |
| [`get_active_chains`](#pylimer_tools_cpp.MEHPForceRelaxation.get_active_chains)(self[, tolerance])                                 | Get the crosslinker chains that are active.                                                                                                              |
| [`get_average_contour_length`](#pylimer_tools_cpp.MEHPForceRelaxation.get_average_contour_length)(self)                            | Get the average contour length of all springs.                                                                                                           |
| [`get_average_spring_length`](#pylimer_tools_cpp.MEHPForceRelaxation.get_average_spring_length)(self)                              | Get the average length of the springs.                                                                                                                   |
| [`get_crosslinker_universe`](#pylimer_tools_cpp.MEHPForceRelaxation.get_crosslinker_universe)(self)                                | Returns the universe [of crosslinkers] with the positions of the current state of the simulation.                                                        |
| [`get_current_spring_distances`](#pylimer_tools_cpp.MEHPForceRelaxation.get_current_spring_distances)(self)                        | Get the current spring distances.                                                                                                                        |
| [`get_dangling_weight_fraction`](#pylimer_tools_cpp.MEHPForceRelaxation.get_dangling_weight_fraction)(self[, tolerance])           | Compute the weight fraction of non-active springs.                                                                                                       |
| [`get_default_r0_square`](#pylimer_tools_cpp.MEHPForceRelaxation.get_default_r0_square)(self)                                      | Returns the value effectively used in [`get_gamma_factor()`](#pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor) for $\langle R_{0,\eta}^2\rangle$. |
| [`get_effective_functionality_of_atoms`](#pylimer_tools_cpp.MEHPForceRelaxation.get_effective_functionality_of_atoms)(self[, ...]) | Returns the number of active springs connected to each atom, atomId used as index.                                                                       |
| [`get_exit_reason`](#pylimer_tools_cpp.MEHPForceRelaxation.get_exit_reason)(self)                                                  | Returns the reason for termination of the simulation.                                                                                                    |
| [`get_force`](#pylimer_tools_cpp.MEHPForceRelaxation.get_force)(self)                                                              | Returns the force at the current state of the simulation.                                                                                                |
| [`get_gamma_factor`](#pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor)(self[, b0_squared, ...])                             | Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:                                                                                       |
| [`get_gamma_factors`](#pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factors)(self[, b0_squared])                                | Computes the gamma factor for each spring as part of the ANT/MEHP formulism.                                                                             |
| [`get_ids_of_active_nodes`](#pylimer_tools_cpp.MEHPForceRelaxation.get_ids_of_active_nodes)(self[, tolerance, ...])                | Get the atom ids of the nodes that are considered active.                                                                                                |
| [`get_nr_of_active_nodes`](#pylimer_tools_cpp.MEHPForceRelaxation.get_nr_of_active_nodes)(self[, tolerance, ...])                  | Get the number of active nodes remaining after running the simulation.                                                                                   |
| [`get_nr_of_active_springs`](#pylimer_tools_cpp.MEHPForceRelaxation.get_nr_of_active_springs)(self[, tolerance])                   | Get the number of active springs remaining after running the simulation.                                                                                 |
| [`get_nr_of_active_springs_connected`](#pylimer_tools_cpp.MEHPForceRelaxation.get_nr_of_active_springs_connected)(self[, ...])     | Returns the number of active springs connected to each node.                                                                                             |
| [`get_nr_of_iterations`](#pylimer_tools_cpp.MEHPForceRelaxation.get_nr_of_iterations)(self)                                        | Returns the number of iterations used for force relaxation.                                                                                              |
| [`get_nr_of_nodes`](#pylimer_tools_cpp.MEHPForceRelaxation.get_nr_of_nodes)(self)                                                  | Get the number of nodes considered in this simulation.                                                                                                   |
| [`get_nr_of_springs`](#pylimer_tools_cpp.MEHPForceRelaxation.get_nr_of_springs)(self)                                              | Get the number of springs considered in this simulation.                                                                                                 |
| [`get_pressure`](#pylimer_tools_cpp.MEHPForceRelaxation.get_pressure)(self)                                                        | Returns the pressure at the current state of the simulation.                                                                                             |
| [`get_residual_norm`](#pylimer_tools_cpp.MEHPForceRelaxation.get_residual_norm)(self)                                              | Returns the residual norm at the current state of the simulation.                                                                                        |
| [`get_residuals`](#pylimer_tools_cpp.MEHPForceRelaxation.get_residuals)(self)                                                      | Returns the residuals at the current state of the simulation.                                                                                            |
| [`get_soluble_weight_fraction`](#pylimer_tools_cpp.MEHPForceRelaxation.get_soluble_weight_fraction)(self[, tolerance])             | Compute the weight fraction of springs connected to active springs (any depth).                                                                          |
| [`get_spring_contour_length`](#pylimer_tools_cpp.MEHPForceRelaxation.get_spring_contour_length)(self)                              | Get the spring contour lengths.                                                                                                                          |
| [`get_spring_distances`](#pylimer_tools_cpp.MEHPForceRelaxation.get_spring_distances)(self)                                        | Get the current coordinate differences for all the springs.                                                                                              |
| [`get_spring_lengths`](#pylimer_tools_cpp.MEHPForceRelaxation.get_spring_lengths)(self)                                            | Get the current lengths for all the springs.                                                                                                             |
| [`get_stress_tensor`](#pylimer_tools_cpp.MEHPForceRelaxation.get_stress_tensor)(self)                                              | Returns the stress tensor at the current state of the simulation.                                                                                        |
| [`requires_another_run`](#pylimer_tools_cpp.MEHPForceRelaxation.requires_another_run)(self)                                        | For performance reasons, the objective is only minimised within the distances of one box.                                                                |
| [`run_force_relaxation`](#pylimer_tools_cpp.MEHPForceRelaxation.run_force_relaxation)(self[, algorithm, ...])                      | Run the simulation.                                                                                                                                      |
| [`set_force_evaluator`](#pylimer_tools_cpp.MEHPForceRelaxation.set_force_evaluator)(self, force_evaluator)                         | Reset the currently used force evaluator.                                                                                                                |

### Attributes Documentation

#### network

Get the network structure.

* **Returns:**
  The network structure used in the simulation

### Methods Documentation

#### assume_box_large_enough(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), box_large_enough: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Configure whether to run PBC on the bonds or not.

* **Parameters:**
  **box_large_enough** – If True, assume the box is large enough and don’t apply PBC on bonds.
  If your bonds could get larger than half the box length,
  this must be kept False (default).

#### config_rerun_epsilon(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), epsilon: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [None](https://docs.python.org/3/library/constants.html#None)

Configure the offset from the lower and upper bounds for the simulation to suggest another run.

* **Parameters:**
  **epsilon** – The epsilon value to use for the rerun check
  (See: [`requires_another_run()`](#pylimer_tools_cpp.MEHPForceRelaxation.requires_another_run))

#### config_step_output(\*args, \*\*kwargs)

Overloaded function.

1. config_step_output(self: pylimer_tools_cpp.MEHPForceRelaxation, output_configuration: collections.abc.Sequence[pylimer_tools_cpp.OutputConfiguration]) -> None
   > Set which values to log during the simulation.
   > * **param output_configuration:**
   >   An OutputConfiguration struct or list of OutputConfiguration structs
   >   specifying what values to log and how often
2. config_step_output(self: pylimer_tools_cpp.MEHPForceRelaxation, output_configuration: collections.abc.Sequence[pylimer_tools_cpp.OutputConfiguration]) -> None
   > Set which values to log during the simulation.
   > * **param output_configuration:**
   >   An OutputConfiguration struct or list of OutputConfiguration structs
   >   specifying what values to log and how often

#### count_active_clustered_atoms(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [float](https://docs.python.org/3/library/functions.html#float)

Counts the active clustered atoms in the system.

* **Parameters:**
  **tolerance** – Springs under this length are considered inactive
* **Returns:**
  The number of active clustered atoms

#### get_active_chains(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Molecule](pylimer_tools_cpp.Molecule.md#pylimer_tools_cpp.Molecule)]

Get the crosslinker chains that are active.

* **Parameters:**
  **tolerance** – Springs under this length are considered inactive
* **Returns:**
  List of active crosslinker chains

#### get_average_contour_length(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the average contour length of all springs.

* **Returns:**
  The average contour length

#### get_average_spring_length(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the average length of the springs. Note that in contrast to [`get_gamma_factor()`](#pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor),
this value is normalized by the number of springs rather than the number of chains.

* **Returns:**
  The average spring length

#### get_crosslinker_universe(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

Returns the universe [of crosslinkers] with the positions of the current state of the simulation.

* **Returns:**
  A Universe object containing the crosslinkers with updated positions

#### get_current_spring_distances(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current spring distances.

* **Returns:**
  A vector containing the current spring distance vectors

#### get_dangling_weight_fraction(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of non-active springs.

Caution: ignores atom masses.

* **Parameters:**
  **tolerance** – Springs under this length are considered inactive
* **Returns:**
  The dangling weight fraction

#### get_default_r0_square(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [float](https://docs.python.org/3/library/functions.html#float)

Returns the value effectively used in [`get_gamma_factor()`](#pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor) for $\langle R_{0,\eta}^2\rangle$.

* **Returns:**
  The default R0 squared value used in gamma factor calculations

#### get_effective_functionality_of_atoms(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [int](https://docs.python.org/3/library/functions.html#int)]

Returns the number of active springs connected to each atom, atomId used as index.

* **Parameters:**
  **tolerance** – Springs under this length are considered inactive
* **Returns:**
  Vector with effective functionality for each atom (indexed by atom ID)

#### get_exit_reason(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [pylimer_tools_cpp.ExitReason](pylimer_tools_cpp.ExitReason.md#pylimer_tools_cpp.ExitReason)

Returns the reason for termination of the simulation.

* **Returns:**
  The exit reason enum value indicating why the simulation ended

#### get_force(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [float](https://docs.python.org/3/library/functions.html#float)

Returns the force at the current state of the simulation.

* **Returns:**
  The current force value

#### get_gamma_factor(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), b0_squared: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = -1.0, nr_of_chains: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [float](https://docs.python.org/3/library/functions.html#float)

Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

$\Gamma = \langle\gamma_{\eta}\rangle$, with $\gamma_{\eta} = \\frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}$,
which you can use as $G_{\mathrm{ANT}} = \Gamma \nu k_B T$,
where $\eta$ is the index of a particular strand,
$R_{0}^2$ is the melt mean square end to end distance, in phantom systems $= N_{\eta} b^2$$,
$N_{\eta}$ is the number of atoms in this strand $\eta$,
$b$ its mean square bond length,
$\nu$ the volume fraction,
$T$ the temperature and
$k_B$ Boltzmann’s constant.

* **Parameters:**
  * **b0_squared** – Part of the denominator in the equation of $\Gamma$.
    If $-1.0$ (default), the network is used for determination (which is not accurate), the system is assumed to be phantom.
    For real systems, the value could be determined by [`compute_mean_squared_end_to_end_distance()`](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance)
    on the melt system, with subsequent division by the nr of bonds in the chain.
  * **nr_of_chains** – The value to normalize the sum of square distances by. Usually (and default if $< 0$) the nr of chains.

#### get_gamma_factors(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), b0_squared: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = -1.0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Computes the gamma factor for each spring as part of the ANT/MEHP formulism.

$\gamma_{\eta} = \\frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}$, with (here)
$R_{0,\eta}^2 = N_\eta \cdot ` the parameter `b0_squared$.
You can obtain this parameter e.g. by doing melt simulations at different lengths,
it’s the slope you obtain.

* **Parameters:**
  **b0_squared** – Part of the denominator in the equation of $\Gamma$.
  If $-1.0$ (default), the network is used for determination (which is not accurate), the system is assumed to be phantom.
  For real systems, the value could be determined by [`compute_mean_squared_end_to_end_distance()`](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance)
  on the melt system, with subsequent division by the nr of bonds in the chain.

See also [`get_gamma_factor()`](#pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor) for the mean of these.

#### get_ids_of_active_nodes(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001, minimum_nr_of_active_connections: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, maximum_nr_of_active_connections: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get the atom ids of the nodes that are considered active.

* **Parameters:**
  * **tolerance** – Springs under this length are considered inactive. A node is active if it has > 2 active springs.
  * **minimum_nr_of_active_connections** – Minimum number of active connections required for a node to be considered active
  * **maximum_nr_of_active_connections** – Maximum number of active connections allowed for a node to be considered active
* **Returns:**
  List of atom IDs for active nodes

#### get_nr_of_active_nodes(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001, minimumNrOfActiveConnections: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, maximumNrOfActiveConnections: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of active nodes remaining after running the simulation.

* **Parameters:**
  * **tolerance** – Springs under this length are considered inactive
  * **minimumNrOfActiveConnections** – A node is active if it has equal or more than this number of active springs
  * **maximumNrOfActiveConnections** – A node is active if it has equal or less than this number of active springs.
    Use a value < 0 to indicate that there is no maximum number of active connections
* **Returns:**
  The number of active nodes

#### get_nr_of_active_springs(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of active springs remaining after running the simulation.

* **Parameters:**
  **tolerance** – Springs under this length are considered inactive
* **Returns:**
  The number of active springs

#### get_nr_of_active_springs_connected(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.05) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.int32], '[m, 1]']

Returns the number of active springs connected to each node.

* **Parameters:**
  **tolerance** – Springs under this length are considered inactive
* **Returns:**
  Vector with the number of active springs for each node

#### get_nr_of_iterations(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [int](https://docs.python.org/3/library/functions.html#int)

Returns the number of iterations used for force relaxation.

* **Returns:**
  The number of iterations performed during force relaxation

#### get_nr_of_nodes(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of nodes considered in this simulation.

* **Returns:**
  The number of nodes in the simulation

#### get_nr_of_springs(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of springs considered in this simulation.

* **Parameters:**
  **tolerance** – Springs under this length are considered inactive
* **Returns:**
  The number of springs in the simulation

#### get_pressure(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [float](https://docs.python.org/3/library/functions.html#float)

Returns the pressure at the current state of the simulation.

* **Returns:**
  The current pressure value

#### get_residual_norm(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [float](https://docs.python.org/3/library/functions.html#float)

Returns the residual norm at the current state of the simulation.

* **Returns:**
  The current residual norm value

#### get_residuals(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Returns the residuals at the current state of the simulation.

* **Returns:**
  The current residual vector

#### get_soluble_weight_fraction(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.001) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of springs connected to active
springs (any depth).

Caution: ignores atom masses.

* **Parameters:**
  **tolerance** – Springs under this length are considered inactive
* **Returns:**
  The soluble weight fraction

#### get_spring_contour_length(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the spring contour lengths.

* **Returns:**
  A vector containing the contour lengths of all springs

#### get_spring_distances(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current coordinate differences for all the springs.

* **Returns:**
  A vector of size 3\*nrOfSprings, with each x, y, z values of the springs

#### get_spring_lengths(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current lengths for all the springs.

* **Returns:**
  A vector of size nrOfSprings, with each the norm of the distances

#### get_stress_tensor(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Returns the stress tensor at the current state of the simulation.

* **Returns:**
  The current stress tensor matrix

#### requires_another_run(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation)) → [bool](https://docs.python.org/3/library/functions.html#bool)

For performance reasons, the objective is only minimised within the distances of one box.
This means, that there is a possibility, e.g. for a single strand longer than two boxes,
that it would not be globally minimised.

If the final displacement of one of the atoms is close
(1e-3, configurable via [`config_rerun_epsilon()`](#pylimer_tools_cpp.MEHPForceRelaxation.config_rerun_epsilon))
to the imposed min/max, after minimizing,
this method would return true.

* **Returns:**
  True if another run is suggested, False otherwise

#### run_force_relaxation(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), algorithm: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'LD_MMA', max_nr_of_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 250000, x_tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1e-12, f_tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1e-09) → [None](https://docs.python.org/3/library/constants.html#None)

Run the simulation.
Note that the final state of the minimization is persisted and reused if you use this method again.
This is useful if you want to run a global optimization first and add a local one afterwards.
As a consequence though, you cannot simply benchmark only this method; you must include the setup.

* **Parameters:**
  * **algorithm** – The algorithm to use for the force relaxation. Choices: see [NLopt Algorithms](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/)
  * **maxNrOfSteps** – The maximum number of steps to do during the simulation.
  * **xTolerance** – The tolerance of the displacements as an exit condition.
  * **fTolerance** – The tolerance of the force as an exit condition.
  * **is2d** – Specify true if you want to evaluate the force relation only in x and y direction.

#### set_force_evaluator(self: [pylimer_tools_cpp.MEHPForceRelaxation](#pylimer_tools_cpp.MEHPForceRelaxation), force_evaluator: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator)) → [None](https://docs.python.org/3/library/constants.html#None)

Reset the currently used force evaluator.

* **Parameters:**
  **force_evaluator** – The new force evaluator to use
