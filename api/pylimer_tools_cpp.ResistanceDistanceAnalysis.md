# ResistanceDistanceAnalysis

### *class* pylimer_tools_cpp.ResistanceDistanceAnalysis(\*args, \*\*kwargs)

Bases: `pybind11_object`

Analysis of topological resistance distances, Kirchhoff index,
and random-walk commute times on graph Laplacians.

Overloaded function.

1. \_\_init_\_(self: pylimer_tools_cpp.ResistanceDistanceAnalysis, spring_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], spring_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], weights: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] = [], num_nodes: typing.SupportsInt | typing.SupportsIndex = 0) -> None
   > Construct from edge lists (spring_from, spring_to), optional edge weights,
   > and optional number of nodes.
   > * **param spring_from:**
   >   Source node indices of edges
   > * **param spring_to:**
   >   Target node indices of edges
   > * **param weights:**
   >   Optional weights (spring constants) of edges (defaults to 1.0)
   > * **param num_nodes:**
   >   Optional total number of nodes (deduced if 0)
2. \_\_init_\_(self: pylimer_tools_cpp.ResistanceDistanceAnalysis, L: scipy.sparse.csc_matrix[numpy.float64]) -> None
   > Construct from a pre-assembled sparse scalar graph Laplacian matrix L.
   > * **param L:**
   >   Symmetric positive semidefinite graph Laplacian matrix

### Methods Summary

| [`compute_all_pairs_commute_time`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_all_pairs_commute_time)(self)              | Compute the N x N matrix of all-pairs random-walk commute times.                              |
|----------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| [`compute_all_pairs_effective_resistance`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_all_pairs_effective_resistance)(self)      | Compute the N x N symmetric matrix of all-pairs effective resistances.                        |
| [`compute_commute_time`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_commute_time)(self, i, j)                  | Compute the random-walk commute time C_ij = 2 \* M \* R_ij between nodes i and j.             |
| [`compute_effective_resistance`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_effective_resistance)(self, i, j)          | Compute the effective resistance between nodes i and j.                                       |
| [`compute_effective_resistance_batched`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_effective_resistance_batched)(self, pairs) | Compute the effective resistance for a batch of node pairs using a single multi-RHS solve.    |
| [`compute_kirchhoff_index`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_kirchhoff_index)(self)                     | Compute the Kirchhoff index K(G) = sum_{i < j} R_ij without allocating the full N x N matrix. |
| [`compute_kirchhoff_index_spectral`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_kirchhoff_index_spectral)(self)            | Compute the Kirchhoff index from the non-zero eigenvalues of the Laplacian.                   |
| [`compute_moore_penrose_pseudo_inverse`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_moore_penrose_pseudo_inverse)(self[, ...]) | Compute the Moore-Penrose pseudo-inverse L^+ of the graph Laplacian.                          |
| [`compute_resistance_centrality`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_resistance_centrality)(self)               | Compute the resistance centrality (N - 1) / sum_{j != i} R_ij for each node.                  |
| [`compute_total_resistance_per_node`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_total_resistance_per_node)(self)           | Compute the sum of effective resistances from each node to all other nodes.                   |
| [`create_complete`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.create_complete)(num_nodes)                        | Create a ResistanceDistanceAnalysis instance for a complete graph K_N.                        |
| [`create_linear_chain`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.create_linear_chain)(num_nodes)                    | Create a ResistanceDistanceAnalysis instance for a linear chain of num_nodes nodes.           |
| [`create_ring`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.create_ring)(num_nodes)                            | Create a ResistanceDistanceAnalysis instance for a ring graph of num_nodes nodes.             |
| [`create_star`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.create_star)(num_arms[, arm_length])               | Create a ResistanceDistanceAnalysis instance for a star graph.                                |
| [`get_connected_components`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_connected_components)(self)                    | Get the list of node indices for each connected component.                                    |
| [`get_laplacian`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_laplacian)(self)                               | Get the sparse graph Laplacian matrix L.                                                      |
| [`get_nr_of_connected_components`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_nr_of_connected_components)(self)              | Get the number of connected components in the graph.                                          |
| [`get_nr_of_edges`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_nr_of_edges)(self)                             | Get the number of non-self-loop edges in the graph.                                           |
| [`get_nr_of_nodes`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_nr_of_nodes)(self)                             | Get the number of nodes in the graph.                                                         |
| [`get_pinned_laplacian`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_pinned_laplacian)(self)                        | Get the sparse pinned Laplacian matrix with reference nodes pinned to 1 on the diagonal.      |
| [`get_total_edge_weight`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_total_edge_weight)(self)                       | Get the total weight of all edges in the graph.                                               |
| [`is_connected`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.is_connected)(self)                                | Check whether the graph consists of a single connected component.                             |
| [`solve_identity_chunked`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.solve_identity_chunked)(self[, chunk_size])        | Solve L_pinned \* X = I in multi-column chunks.                                               |

### Methods Documentation

#### compute_all_pairs_commute_time(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

Compute the N x N matrix of all-pairs random-walk commute times.

* **Returns:**
  N x N matrix of commute times

#### compute_all_pairs_effective_resistance(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), chunk_size: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 128) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

Compute the N x N symmetric matrix of all-pairs effective resistances.

* **Parameters:**
  **chunk_size** – Number of columns per solve chunk (defaults to 128)
* **Returns:**
  N x N matrix of effective resistances

#### compute_commute_time(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), i: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), j: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the random-walk commute time C_ij = 2 \* M \* R_ij between nodes i and j.

* **Parameters:**
  * **i** – First node index
  * **j** – Second node index
* **Returns:**
  Commute time C_ij

#### compute_effective_resistance(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), i: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), j: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the effective resistance between nodes i and j.

* **Parameters:**
  * **i** – First node index
  * **j** – Second node index
* **Returns:**
  Effective resistance R_ij (or inf if in different components)

#### compute_effective_resistance_batched(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), pairs: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)]]) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute the effective resistance for a batch of node pairs using a single multi-RHS solve.

* **Parameters:**
  **pairs** – List of (i, j) node index pairs
* **Returns:**
  Vector of effective resistances for each pair

#### compute_kirchhoff_index(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the Kirchhoff index K(G) = sum_{i < j} R_ij without allocating the full N x N matrix.

* **Returns:**
  Kirchhoff index (or inf if disconnected)

#### compute_kirchhoff_index_spectral(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the Kirchhoff index from the non-zero eigenvalues of the Laplacian.

* **Returns:**
  Kirchhoff index (or inf if disconnected)

#### compute_moore_penrose_pseudo_inverse(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), chunk_size: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 128) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

Compute the Moore-Penrose pseudo-inverse L^+ of the graph Laplacian.

* **Parameters:**
  **chunk_size** – Number of columns per solve chunk (defaults to 128)
* **Returns:**
  N x N Moore-Penrose pseudo-inverse matrix

#### compute_resistance_centrality(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute the resistance centrality (N - 1) / sum_{j != i} R_ij for each node.

* **Returns:**
  Vector of resistance centralities

#### compute_total_resistance_per_node(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute the sum of effective resistances from each node to all other nodes.

* **Returns:**
  Vector of total resistance per node

#### *static* create_complete(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)

Create a ResistanceDistanceAnalysis instance for a complete graph K_N.

* **Parameters:**
  **num_nodes** – Number of nodes in the complete graph (at least 2)
* **Returns:**
  ResistanceDistanceAnalysis instance

#### *static* create_linear_chain(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)

Create a ResistanceDistanceAnalysis instance for a linear chain of num_nodes nodes.

* **Parameters:**
  **num_nodes** – Number of nodes in the linear chain
* **Returns:**
  ResistanceDistanceAnalysis instance

#### *static* create_ring(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)

Create a ResistanceDistanceAnalysis instance for a ring graph of num_nodes nodes.

* **Parameters:**
  **num_nodes** – Number of nodes in the ring (at least 3)
* **Returns:**
  ResistanceDistanceAnalysis instance

#### *static* create_star(num_arms: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arm_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1) → [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)

Create a ResistanceDistanceAnalysis instance for a star graph.

* **Parameters:**
  * **num_arms** – Number of arms (at least 2)
  * **arm_length** – Number of edges per arm (defaults to 1)
* **Returns:**
  ResistanceDistanceAnalysis instance

#### get_connected_components(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[list](https://docs.python.org/3/builtins/stdtypes.html#list)[[int](https://docs.python.org/3/builtins/functions.html#int)]]

Get the list of node indices for each connected component.

* **Returns:**
  List of node index lists per connected component

#### get_laplacian(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]

Get the sparse graph Laplacian matrix L.

* **Returns:**
  Sparse Laplacian matrix

#### get_nr_of_connected_components(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of connected components in the graph.

* **Returns:**
  Number of connected components

#### get_nr_of_edges(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of non-self-loop edges in the graph.

* **Returns:**
  Number of edges

#### get_nr_of_nodes(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of nodes in the graph.

* **Returns:**
  Number of nodes

#### get_pinned_laplacian(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]

Get the sparse pinned Laplacian matrix with reference nodes pinned to 1 on the diagonal.

* **Returns:**
  Sparse pinned Laplacian matrix

#### get_total_edge_weight(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the total weight of all edges in the graph.

* **Returns:**
  Total edge weight

#### is_connected(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

Check whether the graph consists of a single connected component.

* **Returns:**
  True if connected, False otherwise

#### solve_identity_chunked(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), chunk_size: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 128) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

Solve L_pinned \* X = I in multi-column chunks.

* **Parameters:**
  **chunk_size** – Number of columns per solve chunk (defaults to 128)
* **Returns:**
  N x N solution matrix X
