# MEHPNetwork

### *class* pylimer_tools_cpp.MEHPNetwork(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork))

Bases: `pybind11_object`

Canonical polymer network representation shared across MEHP linear-system,
slip-link relaxation, and nonlinear optimization solvers.

Construct an empty MEHPNetwork.

### Attributes Summary

| [`L`](#pylimer_tools_cpp.MEHPNetwork.L)                             |    |
|--------------------------------------------------------------------------------|----|
| [`assume_box_large_enough`](#pylimer_tools_cpp.MEHPNetwork.assume_box_large_enough)       |    |
| [`assume_complete`](#pylimer_tools_cpp.MEHPNetwork.assume_complete)               |    |
| [`box`](#pylimer_tools_cpp.MEHPNetwork.box)                           |    |
| [`coordinates`](#pylimer_tools_cpp.MEHPNetwork.coordinates)                   |    |
| [`default_bond_length`](#pylimer_tools_cpp.MEHPNetwork.default_bond_length)           |    |
| [`default_r02`](#pylimer_tools_cpp.MEHPNetwork.default_r02)                   |    |
| [`link_indices_of_strand`](#pylimer_tools_cpp.MEHPNetwork.link_indices_of_strand)        |    |
| [`link_is_entanglement`](#pylimer_tools_cpp.MEHPNetwork.link_is_entanglement)          |    |
| [`mean_spring_contour_length`](#pylimer_tools_cpp.MEHPNetwork.mean_spring_contour_length)    |    |
| [`mean_strand_contour_length`](#pylimer_tools_cpp.MEHPNetwork.mean_strand_contour_length)    |    |
| [`nr_of_crosslink_swaps_endured`](#pylimer_tools_cpp.MEHPNetwork.nr_of_crosslink_swaps_endured) |    |
| [`nr_of_links`](#pylimer_tools_cpp.MEHPNetwork.nr_of_links)                   |    |
| [`nr_of_nodes`](#pylimer_tools_cpp.MEHPNetwork.nr_of_nodes)                   |    |
| [`nr_of_springs`](#pylimer_tools_cpp.MEHPNetwork.nr_of_springs)                 |    |
| [`nr_of_strands`](#pylimer_tools_cpp.MEHPNetwork.nr_of_strands)                 |    |
| [`old_atom_ids`](#pylimer_tools_cpp.MEHPNetwork.old_atom_ids)                  |    |
| [`old_atom_types`](#pylimer_tools_cpp.MEHPNetwork.old_atom_types)                |    |
| [`remove_primary_loops`](#pylimer_tools_cpp.MEHPNetwork.remove_primary_loops)          |    |
| [`spring_box_offset`](#pylimer_tools_cpp.MEHPNetwork.spring_box_offset)             |    |
| [`spring_contour_length`](#pylimer_tools_cpp.MEHPNetwork.spring_contour_length)         |    |
| [`spring_index_a`](#pylimer_tools_cpp.MEHPNetwork.spring_index_a)                |    |
| [`spring_index_b`](#pylimer_tools_cpp.MEHPNetwork.spring_index_b)                |    |
| [`spring_indices_of_strand`](#pylimer_tools_cpp.MEHPNetwork.spring_indices_of_strand)      |    |
| [`spring_is_entanglement`](#pylimer_tools_cpp.MEHPNetwork.spring_is_entanglement)        |    |
| [`spring_partitions`](#pylimer_tools_cpp.MEHPNetwork.spring_partitions)             |    |
| [`strand_contour_length`](#pylimer_tools_cpp.MEHPNetwork.strand_contour_length)         |    |
| [`strand_index_of_spring`](#pylimer_tools_cpp.MEHPNetwork.strand_index_of_spring)        |    |
| [`strand_indices_of_link`](#pylimer_tools_cpp.MEHPNetwork.strand_indices_of_link)        |    |
| [`vol`](#pylimer_tools_cpp.MEHPNetwork.vol)                           |    |

### Methods Summary

| [`compute_active_weight_fraction`](#pylimer_tools_cpp.MEHPNetwork.compute_active_weight_fraction)(self, u[, ...])   | Compute the weight fraction of active strands.                                                                |
|---------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| [`compute_dangling_weight_fraction`](#pylimer_tools_cpp.MEHPNetwork.compute_dangling_weight_fraction)(self, u[, ...]) | Compute the weight fraction of dangling strands.                                                              |
| [`compute_soluble_weight_fraction`](#pylimer_tools_cpp.MEHPNetwork.compute_soluble_weight_fraction)(self, u[, ...])  | Compute the weight fraction of soluble strands.                                                               |
| [`deform_to`](#pylimer_tools_cpp.MEHPNetwork.deform_to)(\*args, \*\*kwargs)                    | Overloaded function.                                                                                          |
| [`evaluate_spring_lengths`](#pylimer_tools_cpp.MEHPNetwork.evaluate_spring_lengths)(self, u[, is_2d])        | Evaluate Euclidean lengths of all springs.                                                                    |
| [`evaluate_spring_vectors`](#pylimer_tools_cpp.MEHPNetwork.evaluate_spring_vectors)(self, u[, is_2d])        | Evaluate all flattened spring end-to-end vectors (size 3 \* nr_of_springs).                                   |
| [`find_active_springs`](#pylimer_tools_cpp.MEHPNetwork.find_active_springs)(self, u[, tolerance, is_2d]) | Return boolean mask of active springs.                                                                        |
| [`find_active_strands`](#pylimer_tools_cpp.MEHPNetwork.find_active_strands)(self, u[, tolerance, is_2d]) | Return boolean mask of active strands.                                                                        |
| [`from_crosslinker_universe`](#pylimer_tools_cpp.MEHPNetwork.from_crosslinker_universe)(box, ...[, ...])       | Build a canonical MEHPNetwork from a coarse-grained CrosslinkerUniverse and optional StrandPairEntanglements. |
| [`from_universe`](#pylimer_tools_cpp.MEHPNetwork.from_universe)(universe[, entanglements, ...])    | Build a canonical MEHPNetwork from an atomistic Universe using DFS strand contraction.                        |
| [`from_universe_with_random_entanglements`](#pylimer_tools_cpp.MEHPNetwork.from_universe_with_random_entanglements)(...)     | Sample random entanglements from an atomistic Universe and build a canonical MEHPNetwork.                     |
| [`get_crosslinker_verse`](#pylimer_tools_cpp.MEHPNetwork.get_crosslinker_verse)(self[, u])                 | Reconstruct a coarse-grained Universe of crosslinker nodes bonded by polymer strands.                         |
| [`get_gamma_factor`](#pylimer_tools_cpp.MEHPNetwork.get_gamma_factor)(self, u[, b02, ...])            | Compute the topological Gamma factor.                                                                         |
| [`simplify`](#pylimer_tools_cpp.MEHPNetwork.simplify)(self, u, mode[, tolerance, is_2d, ...]) | Run iterative topology simplification until convergence.                                                      |
| [`validate`](#pylimer_tools_cpp.MEHPNetwork.validate)(self[, u])                              | Validate network dimensions, index mappings, and contour invariants.                                          |

### Attributes Documentation

#### L

#### assume_box_large_enough

#### assume_complete

#### box

#### coordinates

#### default_bond_length

#### default_r02

#### link_indices_of_strand

#### link_is_entanglement

#### mean_spring_contour_length

#### mean_strand_contour_length

#### nr_of_crosslink_swaps_endured

#### nr_of_links

#### nr_of_nodes

#### nr_of_springs

#### nr_of_strands

#### old_atom_ids

#### old_atom_types

#### remove_primary_loops

#### spring_box_offset

#### spring_contour_length

#### spring_index_a

#### spring_index_b

#### spring_indices_of_strand

#### spring_is_entanglement

#### spring_partitions

#### strand_contour_length

#### strand_index_of_spring

#### strand_indices_of_link

#### vol

### Methods Documentation

#### compute_active_weight_fraction(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the weight fraction of active strands.

#### compute_dangling_weight_fraction(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the weight fraction of dangling strands.

#### compute_soluble_weight_fraction(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the weight fraction of soluble strands.

#### deform_to(\*args, \*\*kwargs)

Overloaded function.

1. deform_to(self: pylimer_tools_cpp.MEHPNetwork, target_box: pylimer_tools_cpp.Box) -> None

Affinely deform network coordinates and baked box offsets to target_box.

1. deform_to(self: pylimer_tools_cpp.MEHPNetwork, target_box: pylimer_tools_cpp.Box, u: object) -> None

Affinely deform network coordinates, displacements u, and baked box offsets to target_box.

#### evaluate_spring_lengths(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate Euclidean lengths of all springs.

#### evaluate_spring_vectors(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate all flattened spring end-to-end vectors (size 3 \* nr_of_springs).

#### find_active_springs(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[[numpy.bool](https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.bool)], '[m, 1]']

Return boolean mask of active springs.

#### find_active_strands(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[[numpy.bool](https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.bool)], '[m, 1]']

Return boolean mask of active strands.

#### *static* from_crosslinker_universe(box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.html.md#pylimer_tools_cpp.Box), crosslinker_universe: [pylimer_tools_cpp.CrosslinkerUniverse](pylimer_tools_cpp.CrosslinkerUniverse.html.md#pylimer_tools_cpp.CrosslinkerUniverse), entanglements: [pylimer_tools_cpp.StrandPairEntanglements](pylimer_tools_cpp.StrandPairEntanglements.html.md#pylimer_tools_cpp.StrandPairEntanglements) = StrandPairEntanglements(), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2, entanglement_mode: [pylimer_tools_cpp.EntanglementMode](pylimer_tools_cpp.EntanglementMode.html.md#pylimer_tools_cpp.EntanglementMode) = EntanglementMode.MERGED_LINKS, remove_dangling_chains: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, include_free_chains: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, remove_primary_loops: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, assume_box_large_enough: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork)

Build a canonical MEHPNetwork from a coarse-grained CrosslinkerUniverse and optional StrandPairEntanglements.

* **Parameters:**
  * **box** – Simulation box.
  * **crosslinker_universe** – Crosslinker & strand representation.
  * **entanglements** – Optional pre-sampled strand-pair entanglements.
  * **crosslinker_type** – Atom type identifying crosslinkers.
  * **entanglement_mode** – Entanglement representation mode.
  * **remove_dangling_chains** – Omit dangling chains when True.
  * **include_free_chains** – Include free/soluble chains when True.
  * **remove_primary_loops** – Remove polymer primary loops upon construction.
  * **assume_box_large_enough** – Assume box large enough for minimum-image PBC.
* **Returns:**
  Constructed MEHPNetwork.

#### *static* from_universe(universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe), entanglements: [pylimer_tools_cpp.AtomPairEntanglements](pylimer_tools_cpp.AtomPairEntanglements.html.md#pylimer_tools_cpp.AtomPairEntanglements) = AtomPairEntanglements(), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2, entanglement_mode: [pylimer_tools_cpp.EntanglementMode](pylimer_tools_cpp.EntanglementMode.html.md#pylimer_tools_cpp.EntanglementMode) = EntanglementMode.MERGED_LINKS, keep_2functional_crosslinkers: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, remove_primary_loops: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, assume_box_large_enough: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork)

Build a canonical MEHPNetwork from an atomistic Universe using DFS strand contraction.

* **Parameters:**
  * **universe** – Atomistic polymer Universe.
  * **entanglements** – Optional pre-sampled atom-pair entanglements.
  * **crosslinker_type** – Atom type identifying crosslinker beads.
  * **entanglement_mode** – Entanglement representation mode.
  * **keep_2functional_crosslinkers** – Keep degree-2 crosslinkers as explicit network nodes.
  * **remove_primary_loops** – Remove polymer primary loops upon construction.
  * **assume_box_large_enough** – Whether springs can be evaluated with minimum-image PBC without baked offsets.
* **Returns:**
  Constructed MEHPNetwork.

#### *static* from_universe_with_random_entanglements(universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe), nr_of_entanglements_to_sample: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), upper_sampling_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), lower_sampling_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.0, min_nr_of_entanglements: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0, same_strand_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.0, seed: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = '', crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2, entanglement_mode: [pylimer_tools_cpp.EntanglementMode](pylimer_tools_cpp.EntanglementMode.html.md#pylimer_tools_cpp.EntanglementMode) = EntanglementMode.MERGED_LINKS, with_crosslinkers: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True, filter_dangling_and_soluble: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, keep_2functional_crosslinkers: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, remove_primary_loops: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, assume_box_large_enough: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork)

Sample random entanglements from an atomistic Universe and build a canonical MEHPNetwork.

* **Parameters:**
  * **universe** – Atomistic polymer Universe.
  * **nr_of_entanglements_to_sample** – Target number of entanglements to sample.
  * **upper_sampling_cutoff** – Upper distance cutoff for pair sampling.
  * **lower_sampling_cutoff** – Lower distance cutoff for pair sampling.
  * **min_nr_of_entanglements** – Minimum entanglements required if possible.
  * **same_strand_cutoff** – Minimum chemical distance for same-strand pairs.
  * **seed** – PRNG seed string (empty for random seed).
  * **crosslinker_type** – Atom type identifying crosslinker beads.
  * **entanglement_mode** – Entanglement representation mode.
  * **with_crosslinkers** – Exclude crosslinker beads as entanglement targets when True.
  * **filter_dangling_and_soluble** – Omit entanglements touching dangling/soluble strands during sampling.
  * **keep_2functional_crosslinkers** – Keep degree-2 crosslinkers as explicit nodes.
  * **remove_primary_loops** – Remove polymer primary loops upon construction.
  * **assume_box_large_enough** – Assume box large enough for minimum-image PBC.
* **Returns:**
  Constructed MEHPNetwork.

#### get_crosslinker_verse(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'] = array([], dtype=float64)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe)

Reconstruct a coarse-grained Universe of crosslinker nodes bonded by polymer strands.

#### get_gamma_factor(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], b02: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0, nr_of_chains: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [float](https://docs.python.org/3/builtins/functions.html#float)

Compute the topological Gamma factor.

#### simplify(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [object](https://docs.python.org/3/builtins/functions.html#object), mode: [pylimer_tools_cpp.StructureSimplificationMode](pylimer_tools_cpp.StructureSimplificationMode.html.md#pylimer_tools_cpp.StructureSimplificationMode), tolerance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-09, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, breaking_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0) → [int](https://docs.python.org/3/builtins/functions.html#int)

Run iterative topology simplification until convergence.

#### validate(self: [pylimer_tools_cpp.MEHPNetwork](#pylimer_tools_cpp.MEHPNetwork), u: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'] = array([], dtype=float64)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Validate network dimensions, index mappings, and contour invariants.
