# NormalModeAnalyzer

### *class* pylimer_tools_cpp.NormalModeAnalyzer(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), spring_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)], spring_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)], nr_of_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, weights: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)] | [None](https://docs.python.org/3/builtins/constants.html#None) = None)

Bases: `pybind11_object`

Stateless graph Laplacian builder and normal mode solver.

Please cite Gusev and Bernhard [[GB24](../acknowledgements.md#id7)] if you use this method in your work.

Initialize the NormalModeAnalyzer from bond/spring edge lists.

* **Parameters:**
  * **spring_from** – Starting node indices for springs/bonds
  * **spring_to** – Ending node indices for springs/bonds
  * **nr_of_nodes** – Optional total number of nodes (to include isolated vertices)
  * **weights** – Optional per-spring weights (defaults to unit weight 1.0)

### Attributes Summary

| [`matrix`](#pylimer_tools_cpp.NormalModeAnalyzer.matrix)               | Assembled sparse connectivity matrix.                           |
|-----------------------------------------------------------------------|-----------------------------------------------------------------|
| [`num_soluble_clusters`](#pylimer_tools_cpp.NormalModeAnalyzer.num_soluble_clusters) | Number of connected components (soluble clusters / zero modes). |
| [`size`](#pylimer_tools_cpp.NormalModeAnalyzer.size)                 | Dimension of the connectivity matrix.                           |

### Methods Summary

| [`chebyshev_moments`](#pylimer_tools_cpp.NormalModeAnalyzer.chebyshev_moments)(self[, config])              | Compute Jackson-damped Chebyshev moments using KPM.                                  |
|-------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| [`density_of_states`](#pylimer_tools_cpp.NormalModeAnalyzer.density_of_states)(self, lambda_grid[, config]) | Compute continuous density of states (DOS) on lambda_grid using KPM.                 |
| [`from_matrix`](#pylimer_tools_cpp.NormalModeAnalyzer.from_matrix)(matrix)                            | Construct a NormalModeAnalyzer directly from a sparse Laplacian/connectivity matrix. |
| [`from_universe`](#pylimer_tools_cpp.NormalModeAnalyzer.from_universe)(universe)                        | Construct a NormalModeAnalyzer from a Universe's atom-bond graph.                    |
| [`solve_approx`](#pylimer_tools_cpp.NormalModeAnalyzer.solve_approx)(self[, config])                   | Compute an approximate / stochastic normal mode spectrum.                            |
| [`solve_exact`](#pylimer_tools_cpp.NormalModeAnalyzer.solve_exact)(self[, solver, ...])               | Compute exact eigenvalues (and optionally eigenvectors).                             |

### Attributes Documentation

#### matrix

Assembled sparse connectivity matrix.

#### num_soluble_clusters

Number of connected components (soluble clusters / zero modes).

#### size

Dimension of the connectivity matrix.

### Methods Documentation

#### chebyshev_moments(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), config: [pylimer_tools_cpp.KPMConfig](pylimer_tools_cpp.KPMConfig.md#pylimer_tools_cpp.KPMConfig) = KPMConfig()) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute Jackson-damped Chebyshev moments using KPM.

* **Parameters:**
  **config** – KPMConfig instance
* **Returns:**
  Vector of damped Chebyshev moments

#### density_of_states(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), lambda_grid: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], config: [pylimer_tools_cpp.KPMConfig](pylimer_tools_cpp.KPMConfig.md#pylimer_tools_cpp.KPMConfig) = KPMConfig()) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute continuous density of states (DOS) on lambda_grid using KPM.

* **Parameters:**
  * **lambda_grid** – Eigenvalue points at which to evaluate DOS
  * **config** – KPMConfig instance
* **Returns:**
  Density of states values on the grid

#### *static* from_matrix(matrix: [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]) → [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)

Construct a NormalModeAnalyzer directly from a sparse Laplacian/connectivity matrix.

* **Parameters:**
  **matrix** – Square sparse connectivity matrix
* **Returns:**
  NormalModeAnalyzer instance

#### *static* from_universe(universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)) → [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)

Construct a NormalModeAnalyzer from a Universe’s atom-bond graph.

* **Parameters:**
  **universe** – Universe instance
* **Returns:**
  NormalModeAnalyzer instance

#### solve_approx(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), config: [pylimer_tools_cpp.SLQConfig](pylimer_tools_cpp.SLQConfig.md#pylimer_tools_cpp.SLQConfig) | [pylimer_tools_cpp.HybridConfig](pylimer_tools_cpp.HybridConfig.md#pylimer_tools_cpp.HybridConfig) | [pylimer_tools_cpp.KPMConfig](pylimer_tools_cpp.KPMConfig.md#pylimer_tools_cpp.KPMConfig) | [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [pylimer_tools_cpp.NormalModeSpectrum](pylimer_tools_cpp.NormalModeSpectrum.md#pylimer_tools_cpp.NormalModeSpectrum)

Compute an approximate / stochastic normal mode spectrum.

* **Parameters:**
  **config** – SLQConfig, HybridConfig, KPMConfig, or None for AUTO
* **Returns:**
  NormalModeSpectrum instance

#### solve_exact(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), solver: [pylimer_tools_cpp.NormalModeSolver](pylimer_tools_cpp.NormalModeSolver.md#pylimer_tools_cpp.NormalModeSolver) = NormalModeSolver.AUTO, compute_eigenvectors: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, decompose_components: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, num_modes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) | [None](https://docs.python.org/3/builtins/constants.html#None) = None) → [pylimer_tools_cpp.NormalModeSpectrum](pylimer_tools_cpp.NormalModeSpectrum.md#pylimer_tools_cpp.NormalModeSpectrum)

Compute exact eigenvalues (and optionally eigenvectors).

* **Parameters:**
  * **solver** – Solver algorithm (NormalModeSolver enum; default: NormalModeSolver.AUTO). Ignored when `num_modes < size`.
  * **compute_eigenvectors** – Whether to also compute eigenvectors (default: False)
  * **decompose_components** – Whether to decompose into connected components (default: True). Ignored when `num_modes < size`.
  * **num_modes** – Optional number of smallest modes to compute via Spectra sparse solver
* **Returns:**
  NormalModeSpectrum instance
