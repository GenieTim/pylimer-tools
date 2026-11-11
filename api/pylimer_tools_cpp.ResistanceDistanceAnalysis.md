# ResistanceDistanceAnalysis

### *class* pylimer_tools_cpp.ResistanceDistanceAnalysis(\*args, \*\*kwargs)

Bases: `pybind11_object`

Overloaded function.

1. \_\_init_\_(self: pylimer_tools_cpp.ResistanceDistanceAnalysis, spring_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], spring_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], weights: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] = [], num_nodes: typing.SupportsInt | typing.SupportsIndex = 0) -> None
2. \_\_init_\_(self: pylimer_tools_cpp.ResistanceDistanceAnalysis, L: scipy.sparse.csc_matrix[numpy.float64]) -> None

### Methods Summary

| [`compute_all_pairs_commute_time`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_all_pairs_commute_time)(self)              |    |
|----------------------------------------------------------------------------------------------------|----|
| [`compute_all_pairs_effective_resistance`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_all_pairs_effective_resistance)(self)      |    |
| [`compute_commute_time`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_commute_time)(self, i, j)                  |    |
| [`compute_effective_resistance`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_effective_resistance)(self, i, j)          |    |
| [`compute_effective_resistance_batched`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_effective_resistance_batched)(self, pairs) |    |
| [`compute_kirchhoff_index`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_kirchhoff_index)(self)                     |    |
| [`compute_kirchhoff_index_spectral`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_kirchhoff_index_spectral)(self)            |    |
| [`compute_moore_penrose_pseudo_inverse`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_moore_penrose_pseudo_inverse)(self[, ...]) |    |
| [`compute_resistance_centrality`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_resistance_centrality)(self)               |    |
| [`compute_total_resistance_per_node`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.compute_total_resistance_per_node)(self)           |    |
| [`create_complete`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.create_complete)(num_nodes)                        |    |
| [`create_linear_chain`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.create_linear_chain)(num_nodes)                    |    |
| [`create_ring`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.create_ring)(num_nodes)                            |    |
| [`create_star`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.create_star)(num_arms[, arm_length])               |    |
| [`get_connected_components`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_connected_components)(self)                    |    |
| [`get_laplacian`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_laplacian)(self)                               |    |
| [`get_nr_of_connected_components`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_nr_of_connected_components)(self)              |    |
| [`get_nr_of_edges`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_nr_of_edges)(self)                             |    |
| [`get_nr_of_nodes`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_nr_of_nodes)(self)                             |    |
| [`get_total_edge_weight`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.get_total_edge_weight)(self)                       |    |
| [`is_connected`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.is_connected)(self)                                |    |
| [`solve_identity_chunked`](#pylimer_tools_cpp.ResistanceDistanceAnalysis.solve_identity_chunked)(self[, chunk_size])        |    |

### Methods Documentation

#### compute_all_pairs_commute_time(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

#### compute_all_pairs_effective_resistance(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), chunk_size: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 128) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

#### compute_commute_time(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), i: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), j: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### compute_effective_resistance(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), i: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), j: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### compute_effective_resistance_batched(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), pairs: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)]]) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### compute_kirchhoff_index(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### compute_kirchhoff_index_spectral(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### compute_moore_penrose_pseudo_inverse(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), chunk_size: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 128) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']

#### compute_resistance_centrality(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### compute_total_resistance_per_node(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

#### *static* create_complete(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)

#### *static* create_linear_chain(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)

#### *static* create_ring(num_nodes: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)

#### *static* create_star(num_arms: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), arm_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1) → [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)

#### get_connected_components(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[list](https://docs.python.org/3/builtins/stdtypes.html#list)[[int](https://docs.python.org/3/builtins/functions.html#int)]]

#### get_laplacian(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [scipy.sparse.csc_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix)[numpy.float64]

#### get_nr_of_connected_components(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_edges(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_nr_of_nodes(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [int](https://docs.python.org/3/builtins/functions.html#int)

#### get_total_edge_weight(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [float](https://docs.python.org/3/builtins/functions.html#float)

#### is_connected(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

#### solve_identity_chunked(self: [pylimer_tools_cpp.ResistanceDistanceAnalysis](#pylimer_tools_cpp.ResistanceDistanceAnalysis), chunk_size: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 128) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, n]']
