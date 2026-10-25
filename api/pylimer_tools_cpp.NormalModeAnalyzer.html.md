# NormalModeAnalyzer

### *class* pylimer_tools_cpp.NormalModeAnalyzer(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), spring_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], spring_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)])

Bases: `pybind11_object`

Compute the normal modes and predict the loss/storage moduli.

Please cite Gusev and Bernhard [[GB24](../acknowledgements.html.md#id7)] if you use this method in your work.

Initialize the NormalModeAnalyzer with the bonds (edges).

Constructs the connectivity matrix from the given edges.

* **Parameters:**
  * **spring_from** – Vector of starting node indices for springs/bonds
  * **spring_to** – Vector of ending node indices for springs/bonds

### Methods Summary

| [`compute_chebyshev_moments_kpm`](#pylimer_tools_cpp.NormalModeAnalyzer.compute_chebyshev_moments_kpm)(self[, ...])       | Compute Chebyshev moments using the Kernel Polynomial Method (KPM).                   |
|---------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| [`compute_density_of_states_kpm`](#pylimer_tools_cpp.NormalModeAnalyzer.compute_density_of_states_kpm)(self, lambda_grid) | Compute continuous density of states (DOS) using KPM.                                 |
| [`evaluate_loss_modulus`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_loss_modulus)(self, omega)               | Evaluate the loss modulus $G''(\omega)$.                                              |
| [`evaluate_loss_modulus_kpm`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_loss_modulus_kpm)(self, omega[, ...])    | Evaluate loss modulus G'' directly using KPM without computing eigenvalues.           |
| [`evaluate_storage_modulus`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_storage_modulus)(self, omega)            | Evaluate the storage modulus $G'(\omega)$.                                            |
| [`evaluate_storage_modulus_kpm`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_storage_modulus_kpm)(self, omega[, ...]) | Evaluate storage modulus G' directly using KPM without computing eigenvalues.         |
| [`evaluate_stress_autocorrelation`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_stress_autocorrelation)(self, t)         | Evaluate stress autocorrelation $C(t)$.                                               |
| [`evaluate_stress_autocorrelation_kpm`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_stress_autocorrelation_kpm)(self, t)     | Evaluate stress autocorrelation directly using KPM without computing eigenvalues.     |
| [`find_all_eigenvalues`](#pylimer_tools_cpp.NormalModeAnalyzer.find_all_eigenvalues)(self[, solver, ...])        | Find all eigenvalues using the specified solver.                                      |
| [`find_sparse_eigenvalues`](#pylimer_tools_cpp.NormalModeAnalyzer.find_sparse_eigenvalues)(self, nr_of_eigenvalues) | Find the k smallest eigenvalues using a sparse solver.                                |
| [`get_eigenvalues`](#pylimer_tools_cpp.NormalModeAnalyzer.get_eigenvalues)(self)                            | Get the eigenvalues.                                                                  |
| [`get_eigenvectors`](#pylimer_tools_cpp.NormalModeAnalyzer.get_eigenvectors)(self)                           | Get eigenvectors.                                                                     |
| [`get_matrix`](#pylimer_tools_cpp.NormalModeAnalyzer.get_matrix)(self)                                 | Get the assembled connectivity matrix.                                                |
| [`get_matrix_size`](#pylimer_tools_cpp.NormalModeAnalyzer.get_matrix_size)(self)                            | Get the size of the matrix (the maximum number of eigenvalues that could be queried). |
| [`get_nr_of_soluble_clusters`](#pylimer_tools_cpp.NormalModeAnalyzer.get_nr_of_soluble_clusters)(self)                 | Get the number of soluble clusters (eigenvalues = 0).                                 |
| [`set_eigenvalues`](#pylimer_tools_cpp.NormalModeAnalyzer.set_eigenvalues)(self, eigenvalues)               | Set the eigenvalues, e.g. if you use an external solver.                              |
| [`set_eigenvectors`](#pylimer_tools_cpp.NormalModeAnalyzer.set_eigenvectors)(self, eigenvectors)             | Set eigenvectors, e.g. if you use an external solver.                                 |

### Methods Documentation

#### compute_chebyshev_moments_kpm(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), num_moments: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1024, num_random_vectors: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 20) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute Chebyshev moments using the Kernel Polynomial Method (KPM).

* **Parameters:**
  * **num_moments** – Number of moments to compute (default: 1024)
  * **num_random_vectors** – Number of random Rademacher vectors for trace estimation (default: 20)
* **Returns:**
  Vector of damped Chebyshev moments

#### compute_density_of_states_kpm(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), lambda_grid: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], num_moments: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1024, num_random_vectors: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 20) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute continuous density of states (DOS) using KPM.

* **Parameters:**
  * **lambda_grid** – Eigenvalue points at which to evaluate DOS
  * **num_moments** – Number of Chebyshev moments (default: 1024)
  * **num_random_vectors** – Number of random vectors (default: 20)
* **Returns:**
  Density of states values on the grid

#### evaluate_loss_modulus(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate the loss modulus $G''(\omega)$. Yet misses the conversion factor.

* **Parameters:**
  **omega** – Angular frequencies
* **Returns:**
  Loss modulus values

#### evaluate_loss_modulus_kpm(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], num_moments: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1024, num_random_vectors: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 20, num_quadrature_points: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2048) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate loss modulus G’’ directly using KPM without computing eigenvalues.

* **Parameters:**
  * **omega** – Angular frequencies to evaluate
  * **num_moments** – Number of Chebyshev moments (default: 1024)
  * **num_random_vectors** – Number of random vectors (default: 20)
  * **num_quadrature_points** – Number of Gauss-Chebyshev quadrature points (default: 2048)
* **Returns:**
  Loss modulus values

#### evaluate_storage_modulus(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate the storage modulus $G'(\omega)$. Yet misses the conversion factor.

* **Parameters:**
  **omega** – Angular frequencies
* **Returns:**
  Storage modulus values

#### evaluate_storage_modulus_kpm(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], num_moments: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1024, num_random_vectors: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 20, num_quadrature_points: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2048) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate storage modulus G’ directly using KPM without computing eigenvalues.

* **Parameters:**
  * **omega** – Angular frequencies to evaluate
  * **num_moments** – Number of Chebyshev moments (default: 1024)
  * **num_random_vectors** – Number of random vectors (default: 20)
  * **num_quadrature_points** – Number of Gauss-Chebyshev quadrature points (default: 2048)
* **Returns:**
  Storage modulus values

#### evaluate_stress_autocorrelation(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), t: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate stress autocorrelation $C(t)$.

* **Parameters:**
  **t** – The time at which to evaluate the stress autocorrelation
* **Returns:**
  Stress autocorrelation values

#### evaluate_stress_autocorrelation_kpm(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), t: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], num_moments: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1024, num_random_vectors: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 20, num_quadrature_points: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2048) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate stress autocorrelation directly using KPM without computing eigenvalues.

* **Parameters:**
  * **t** – Time points to evaluate
  * **num_moments** – Number of Chebyshev moments (default: 1024)
  * **num_random_vectors** – Number of random vectors (default: 20)
  * **num_quadrature_points** – Number of Gauss-Chebyshev quadrature points (default: 2048)
* **Returns:**
  Stress autocorrelation values

#### find_all_eigenvalues(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), solver: [object](https://docs.python.org/3/builtins/functions.html#object) = 'auto', compute_eigenvectors: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, decompose_components: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [None](https://docs.python.org/3/builtins/constants.html#None)

Find all eigenvalues using the specified solver.

* **Parameters:**
  * **solver** – Solver method to use (NormalModeSolver enum or string: ‘auto’, ‘dense_lapack’, ‘dense_eigen’, ‘sparse_slicing’; default: NormalModeSolver.AUTO)
  * **compute_eigenvectors** – Whether to also compute eigenvectors (default: False)
  * **decompose_components** – Whether to decompose the graph into connected components (default: True)

#### find_sparse_eigenvalues(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), nr_of_eigenvalues: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), compute_eigenvectors: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [None](https://docs.python.org/3/builtins/constants.html#None)

Find the k smallest eigenvalues using a sparse solver.

* **Parameters:**
  * **nr_of_eigenvalues** – Number of smallest eigenvalues to find
  * **compute_eigenvectors** – Whether to also compute eigenvectors (default: False)
* **Returns:**
  True if computation was successful

#### get_eigenvalues(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the eigenvalues.

* **Returns:**
  Vector of eigenvalues

#### get_eigenvectors(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

Get eigenvectors.

* **Returns:**
  Matrix of eigenvectors

#### get_matrix(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)) → [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]

Get the assembled connectivity matrix.

* **Returns:**
  The connectivity matrix

#### get_matrix_size(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the size of the matrix (the maximum number of eigenvalues that could be queried).

* **Returns:**
  Size of the connectivity matrix

#### get_nr_of_soluble_clusters(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of soluble clusters (eigenvalues = 0).

* **Returns:**
  Number of soluble clusters

#### set_eigenvalues(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), eigenvalues: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set the eigenvalues, e.g. if you use an external solver.

* **Parameters:**
  **eigenvalues** – Vector of eigenvalues to set

#### set_eigenvectors(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), eigenvectors: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set eigenvectors, e.g. if you use an external solver.

* **Parameters:**
  **eigenvectors** – Matrix of eigenvectors to set
