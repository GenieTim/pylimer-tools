# NormalModeAnalyzer

### *class* pylimer_tools_cpp.NormalModeAnalyzer(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), spring_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], spring_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)])

Bases: `pybind11_object`

Compute the normal modes and predict the loss/storage moduli.

Please cite Gusev and Bernhard [[GB24](../acknowledgements.md#id7)] if you use this method in your work.

Initialize the NormalModeAnalyzer with the bonds (edges).

Constructs the connectivity matrix from the given edges.

* **Parameters:**
  * **spring_from** – Vector of starting node indices for springs/bonds
  * **spring_to** – Vector of ending node indices for springs/bonds

### Methods Summary

| [`evaluate_loss_modulus`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_loss_modulus)(self, omega)                 | Evaluate the loss modulus $G''(\omega)$.                                              |
|---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| [`evaluate_storage_modulus`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_storage_modulus)(self, omega)           | Evaluate the storage modulus $G'(\omega)$.                                            |
| [`evaluate_stress_autocorrelation`](#pylimer_tools_cpp.NormalModeAnalyzer.evaluate_stress_autocorrelation)(self, t) | Evaluate stress autocorrelation $C(t)$.                                               |
| [`find_all_eigenvalues`](#pylimer_tools_cpp.NormalModeAnalyzer.find_all_eigenvalues)(self[, ...])                   | Find all eigenvalues using a dense solver.                                            |
| [`find_sparse_eigenvalues`](#pylimer_tools_cpp.NormalModeAnalyzer.find_sparse_eigenvalues)(self, nr_of_eigenvalues) | Find the k smallest eigenvalues using a sparse solver.                                |
| [`get_eigenvalues`](#pylimer_tools_cpp.NormalModeAnalyzer.get_eigenvalues)(self)                                    | Get the eigenvalues.                                                                  |
| [`get_eigenvectors`](#pylimer_tools_cpp.NormalModeAnalyzer.get_eigenvectors)(self)                                  | Get eigenvectors.                                                                     |
| [`get_matrix`](#pylimer_tools_cpp.NormalModeAnalyzer.get_matrix)(self)                                              | Get the assembled connectivity matrix.                                                |
| [`get_matrix_size`](#pylimer_tools_cpp.NormalModeAnalyzer.get_matrix_size)(self)                                    | Get the size of the matrix (the maximum number of eigenvalues that could be queried). |
| [`get_nr_of_soluble_clusters`](#pylimer_tools_cpp.NormalModeAnalyzer.get_nr_of_soluble_clusters)(self)              | Get the number of soluble clusters (eigenvalues = 0).                                 |
| [`set_eigenvalues`](#pylimer_tools_cpp.NormalModeAnalyzer.set_eigenvalues)(self, eigenvalues)                       | Set the eigenvalues, e.g. if you use an external solver.                              |
| [`set_eigenvectors`](#pylimer_tools_cpp.NormalModeAnalyzer.set_eigenvectors)(self, eigenvectors)                    | Set eigenvectors, e.g. if you use an external solver.                                 |

### Methods Documentation

#### evaluate_loss_modulus(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate the loss modulus $G''(\omega)$. Yet misses the conversion factor.

* **Parameters:**
  **omega** – Angular frequencies
* **Returns:**
  Loss modulus values

#### evaluate_storage_modulus(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate the storage modulus $G'(\omega)$. Yet misses the conversion factor.

* **Parameters:**
  **omega** – Angular frequencies
* **Returns:**
  Storage modulus values

#### evaluate_stress_autocorrelation(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), t: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate stress autocorrelation $C(t)$.

* **Parameters:**
  **t** – The time at which to evaluate the stress autocorrelation
* **Returns:**
  Stress autocorrelation values

#### find_all_eigenvalues(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), compute_eigenvectors: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Find all eigenvalues using a dense solver.

* **Parameters:**
  **compute_eigenvectors** – Whether to also compute eigenvectors (default: False)
* **Returns:**
  True if computation was successful

#### find_sparse_eigenvalues(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), nr_of_eigenvalues: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), compute_eigenvectors: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

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

#### get_matrix(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)) → scipy.sparse.csc_matrix[numpy.float64]

Get the assembled connectivity matrix.

* **Returns:**
  The connectivity matrix

#### get_matrix_size(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the size of the matrix (the maximum number of eigenvalues that could be queried).

* **Returns:**
  Size of the connectivity matrix

#### get_nr_of_soluble_clusters(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of soluble clusters (eigenvalues = 0).

* **Returns:**
  Number of soluble clusters

#### set_eigenvalues(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), eigenvalues: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [None](https://docs.python.org/3/library/constants.html#None)

Set the eigenvalues, e.g. if you use an external solver.

* **Parameters:**
  **eigenvalues** – Vector of eigenvalues to set

#### set_eigenvectors(self: [pylimer_tools_cpp.NormalModeAnalyzer](#pylimer_tools_cpp.NormalModeAnalyzer), eigenvectors: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]']) → [None](https://docs.python.org/3/library/constants.html#None)

Set eigenvectors, e.g. if you use an external solver.

* **Parameters:**
  **eigenvectors** – Matrix of eigenvectors to set
