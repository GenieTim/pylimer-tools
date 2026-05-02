# MCUniverseGenerator

### *class* pylimer_tools_cpp.MCUniverseGenerator(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), lx: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 10.0, ly: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 10.0, lz: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 10.0)

Bases: `pybind11_object`

A [`pylimer_tools_cpp.Universe`](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) generator using a Monte-Carlo procedure.
This generator creates phantom bead-spring polymer networks with crosslinkers and Gaussian chain and bond statistics.

#### NOTE
A note on the characteristic ratio $C_{\infty}$ and the bead distance $b$:
If you use parameters from `pylimer_tools.io.bead_spring_parameter_provider` with
type `pylimer_tools.io.bead_spring_parameter_provider.ParameterType.GAUSSIAN`,
there is no need to set $C_{\infty}$ in the MCUniverseGenerator methods,
as the bead distance $b$ you set via [`set_bead_distance()`](#pylimer_tools_cpp.MCUniverseGenerator.set_bead_distance)
from `pylimer_tools.io.bead_spring_parameter_provider.Parameters.get("<b>")()`
already incorporates the characteristic ratio.

This is the recommended way to use the MCUniverseGenerator, since using e.g.
Kuhn segment parameters conflicts with the generator producing Gaussian bond statistics.

Please cite Gusev and Bernhard [[GB24](../acknowledgements.md#id7)] and/or Bernhard and Gusev [[BG25](../acknowledgements.md#id12)] if you use this method in your work.

Initialize a new MCUniverseGenerator with specified box dimensions.

* **Parameters:**
  * **lx** – Box length in x-direction (default: 10.0).
  * **ly** – Box length in y-direction (default: 10.0).
  * **lz** – Box length in z-direction (default: 10.0).

### Methods Summary

| [`add_crosslinkers`](#pylimer_tools_cpp.MCUniverseGenerator.add_crosslinkers)(self, nr_of_crosslinkers[, ...])                           | Add crosslinkers at random positions.                                                                                               |
|------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
| [`add_crosslinkers_at`](#pylimer_tools_cpp.MCUniverseGenerator.add_crosslinkers_at)(self, coordinates[, ...])                            | Add crosslinkers at specific coordinates.                                                                                           |
| [`add_end_functionalized_strands`](#pylimer_tools_cpp.MCUniverseGenerator.add_end_functionalized_strands)(self, ...[, ...])              | Add strands with crosslinkers at the ends.                                                                                          |
| [`add_functionalized_strands`](#pylimer_tools_cpp.MCUniverseGenerator.add_functionalized_strands)(self, ...[, ...])                      | Add strands with custom crosslink placement using a selector function.                                                              |
| [`add_monofunctional_strands`](#pylimer_tools_cpp.MCUniverseGenerator.add_monofunctional_strands)(self, ...[, ...])                      | Add multiple monofunctional strands with specified bead types.                                                                      |
| [`add_randomly_functionalized_strands`](#pylimer_tools_cpp.MCUniverseGenerator.add_randomly_functionalized_strands)(self, ...)           | Add strands with randomly distributed crosslinkers in between.                                                                      |
| [`add_regularly_spaced_functionalized_strands`](#pylimer_tools_cpp.MCUniverseGenerator.add_regularly_spaced_functionalized_strands)(...) | Add strands with regularly spaced crosslinkers in between.                                                                          |
| [`add_solvent_chains`](#pylimer_tools_cpp.MCUniverseGenerator.add_solvent_chains)(self, ...[, ...])                                      | Randomly distribute additional, free chains.                                                                                        |
| [`add_star_crosslinkers`](#pylimer_tools_cpp.MCUniverseGenerator.add_star_crosslinkers)(self, nr_of_stars, ...)                          | Add star-like crosslinkers with pre-connected strands (useful e.g. for tetra-PEG networks).                                         |
| [`add_strands`](#pylimer_tools_cpp.MCUniverseGenerator.add_strands)(self, nr_of_strands, strand_lengths)                                 | Add strands.                                                                                                                        |
| [`config_nr_of_mc_steps`](#pylimer_tools_cpp.MCUniverseGenerator.config_nr_of_mc_steps)(self[, n_steps])                                 | Set the number of Monte-Carlo steps during bond length equilibration.                                                               |
| [`config_primary_loop_probability`](#pylimer_tools_cpp.MCUniverseGenerator.config_primary_loop_probability)(self[, ...])                 | Configure an additional weight reduction for the primary loop formation.                                                            |
| [`config_secondary_loop_probability`](#pylimer_tools_cpp.MCUniverseGenerator.config_secondary_loop_probability)(self[, ...])             | Configure an additional weight reduction for the secondary loop formation.                                                          |
| [`disable_max_distance`](#pylimer_tools_cpp.MCUniverseGenerator.disable_max_distance)(self)                                              | Disable the maximum distance provider to allow unlimited distance sampling.                                                         |
| [`get_current_crosslinker_conversion`](#pylimer_tools_cpp.MCUniverseGenerator.get_current_crosslinker_conversion)(self)                  | Get the current conversion of crosslinkers, i.e., the fraction of crosslinkers that have been linked to strands.                    |
| [`get_current_nr_of_atoms`](#pylimer_tools_cpp.MCUniverseGenerator.get_current_nr_of_atoms)(self)                                        | Get the current number of atoms that the universe would/will have.                                                                  |
| [`get_current_nr_of_bonds`](#pylimer_tools_cpp.MCUniverseGenerator.get_current_nr_of_bonds)(self)                                        | Get the current number of bonds that the universe would/will have.                                                                  |
| [`get_current_strand_conversion`](#pylimer_tools_cpp.MCUniverseGenerator.get_current_strand_conversion)(self)                            | Get the current conversion of strands, i.e., the fraction of strand ends that have been consumed.                                   |
| [`get_force_balance`](#pylimer_tools_cpp.MCUniverseGenerator.get_force_balance)(self)                                                    | Get an instance of the force balance procedure.                                                                                     |
| [`get_force_balance2`](#pylimer_tools_cpp.MCUniverseGenerator.get_force_balance2)(self)                                                  | Get an instance of the force balance procedure.                                                                                     |
| [`get_force_relaxation`](#pylimer_tools_cpp.MCUniverseGenerator.get_force_relaxation)(self)                                              | Get an instance of the force relaxation procedure.                                                                                  |
| [`get_mean_bead_distance`](#pylimer_tools_cpp.MCUniverseGenerator.get_mean_bead_distance)(self)                                          | Get the currently configured mean bead distance.                                                                                    |
| [`get_mean_squared_bead_distance`](#pylimer_tools_cpp.MCUniverseGenerator.get_mean_squared_bead_distance)(self)                          | Get the currently configured mean squared bead distance.                                                                            |
| [`get_nr_of_atoms`](#pylimer_tools_cpp.MCUniverseGenerator.get_nr_of_atoms)(self)                                                        | Get the current number of atoms that the universe would/will have.                                                                  |
| [`get_nr_of_bonds`](#pylimer_tools_cpp.MCUniverseGenerator.get_nr_of_bonds)(self)                                                        | Get the current number of bonds that the universe would/will have.                                                                  |
| [`get_universe`](#pylimer_tools_cpp.MCUniverseGenerator.get_universe)(self)                                                              | Fetch the current (or final) state of the universe.                                                                                 |
| [`link_strand`](#pylimer_tools_cpp.MCUniverseGenerator.link_strand)(self, strand_idx[, c_infinity])                                      | Link one particular strand to a crosslinker.                                                                                        |
| [`link_strand_to`](#pylimer_tools_cpp.MCUniverseGenerator.link_strand_to)(self, strand_idx, link_idx)                                    | Link a strand to a specific crosslinker.                                                                                            |
| [`link_strand_to_strand`](#pylimer_tools_cpp.MCUniverseGenerator.link_strand_to_strand)(self[, c_infinity])                              | Link two free strand ends together directly.                                                                                        |
| [`link_strands_callback`](#pylimer_tools_cpp.MCUniverseGenerator.link_strands_callback)(self, linking_controller)                        | Link strands to crosslinkers using a custom callback function to control when to stop.                                              |
| [`link_strands_to_conversion`](#pylimer_tools_cpp.MCUniverseGenerator.link_strands_to_conversion)(self, ...[, ...])                      | Actually link the previously added strands to the previously added crosslinkers, until a certain crosslink conversion is reached.   |
| [`link_strands_to_soluble_fraction`](#pylimer_tools_cpp.MCUniverseGenerator.link_strands_to_soluble_fraction)(self, ...)                 | Actually link the previously added strands to the previously added crosslinkers, until a certain soluble fraction is reached.       |
| [`link_strands_to_strands_to_conversion`](#pylimer_tools_cpp.MCUniverseGenerator.link_strands_to_strands_to_conversion)(self, ...)       | Link free strand ends to each other until target conversion is reached.                                                             |
| [`relax_crosslinks`](#pylimer_tools_cpp.MCUniverseGenerator.relax_crosslinks)(self)                                                      | Run force relaxation with the crosslinkers and their strands, to have the crosslinks in their statistically most probable position. |
| [`remove_soluble_fraction`](#pylimer_tools_cpp.MCUniverseGenerator.remove_soluble_fraction)(self[, rescale_box])                         | Remove soluble fraction (as determined by phantom force relaxation) of the strands and crosslinks.                                  |
| [`set_bead_distance`](#pylimer_tools_cpp.MCUniverseGenerator.set_bead_distance)(self, distance[, ...])                                   | Set the mean distance between beads when doing MC stepping.                                                                         |
| [`set_mean_squared_bead_distance`](#pylimer_tools_cpp.MCUniverseGenerator.set_mean_squared_bead_distance)(self, ...[, ...])              | Set the mean squared distance between beads.                                                                                        |
| [`set_seed`](#pylimer_tools_cpp.MCUniverseGenerator.set_seed)(self, seed)                                                                | Set the seed for the random generator.                                                                                              |
| [`use_linear_max_distance`](#pylimer_tools_cpp.MCUniverseGenerator.use_linear_max_distance)(self, multiplier)                            | Converts the $N$ to a maximum distance within which to sample.                                                                      |
| [`use_zscore_max_distance`](#pylimer_tools_cpp.MCUniverseGenerator.use_zscore_max_distance)(self, std_multiplier)                        | Converts the $N$ to a maximum distance within which to sample.                                                                      |
| [`validate`](#pylimer_tools_cpp.MCUniverseGenerator.validate)(self)                                                                      | Check whether the internal state of the generator is valid.                                                                         |

### Methods Documentation

#### add_crosslinkers(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_crosslinkers: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), crosslinker_functionality: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 4, crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, white_noise: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Add crosslinkers at random positions.

* **Parameters:**
  * **nr_of_crosslinkers** – Number of crosslinkers to add.
  * **crosslinker_functionality** – Functionality of the crosslinkers (default: 4).
  * **crosslinker_type** – Atom type of the crosslinkers (default: 2).
  * **white_noise** – Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).

#### add_crosslinkers_at(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), coordinates: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], crosslinker_functionality: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 4, crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2) → [None](https://docs.python.org/3/library/constants.html#None)

Add crosslinkers at specific coordinates.

* **Parameters:**
  * **coordinates** – Coordinates of the crosslinkers as a flat array [x1, y1, z1, x2, y2, z2, …].
  * **crosslinker_functionality** – Functionality of the crosslinkers (default: 4).
  * **crosslinker_type** – Atom type of the crosslinkers (default: 2).

#### add_end_functionalized_strands(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_strands: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), strand_length: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], crosslinker_functionality: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 4, crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, strand_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1, white_noise: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Add strands with crosslinkers at the ends.

* **Parameters:**
  * **nr_of_strands** – Number of strands to add.
  * **strand_length** – Length of each strand.
  * **crosslinker_functionality** – Functionality of the crosslinkers (default: 4).
  * **crosslinker_type** – Atom type of the crosslinkers (default: 2).
  * **strand_atom_type** – Atom type of the beads that are not at the ends (default: 1).
  * **white_noise** – Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).

#### add_functionalized_strands(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_strands: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), strand_length: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], crosslink_selector: [collections.abc.Callable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable)[[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[bool](https://docs.python.org/3/library/functions.html#bool), [int](https://docs.python.org/3/library/functions.html#int)]], default_crosslinker_functionality: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, strand_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1, white_noise: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Add strands with custom crosslink placement using a selector function.

This is a general implementation that allows for flexible crosslink placement patterns.
The crosslink_selector function is called for each bead position and should return
a tuple (should_convert, functionality) indicating whether that bead should be
converted to a crosslink and what functionality it should have.

Example usage:

```default
# Create a custom pattern with crosslinks every 5 beads starting at position 2
def my_selector(strand_index, bead_index, total_beads):
    if bead_index >= 2 and (bead_index - 2) % 5 == 0:
        return (True, 4)  # Convert to crosslink with functionality 4
    return (False, 0)     # Don't convert

generator.add_functionalized_strands(
    nr_of_strands=10,
    strand_length=[20] * 10,
    crosslink_selector=my_selector,
    default_crosslinker_functionality=4
)
```

* **Parameters:**
  * **nr_of_strands** – Number of strands to add.
  * **strand_length** – Length of each strand (list of integers).
  * **crosslink_selector** – Function that takes (strand_index, bead_index, total_beads) and returns (should_convert, functionality).
  * **default_crosslinker_functionality** – Default functionality for crosslinks (used for positioning).
  * **crosslinker_type** – Atom type of the crosslinkers (default: 2).
  * **strand_atom_type** – Atom type of the beads that stay (default: 1).
  * **white_noise** – Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).

#### add_monofunctional_strands(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_monofunctional_strands: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), monofunctional_strand_length: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], monofunctional_strand_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 4) → [None](https://docs.python.org/3/library/constants.html#None)

Add multiple monofunctional strands with specified bead types.

* **Parameters:**
  * **nr_of_monofunctional_strands** – Number of monofunctional strands to add.
  * **monofunctional_strand_length** – Vector specifying the length of each strand in beads.
  * **monofunctional_strand_atom_type** – Atom type for the strand beads (default: 4).

#### add_randomly_functionalized_strands(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_strands: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), strand_length: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], functionalization_probability: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), crosslinker_functionality: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 4, crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, strand_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1, white_noise: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Add strands with randomly distributed crosslinkers in between.

* **Parameters:**
  * **nr_of_strands** – Number of strands to add.
  * **strand_length** – Length of each strand.
  * **functionalization_probability** – Proportion of beads that are made crosslink.
  * **crosslinker_functionality** – Functionality of the crosslinkers (default: 4).
  * **crosslinker_type** – Atom type of the crosslinkers (default: 2).
  * **strand_atom_type** – Atom type of the beads that stay (default: 1).
  * **white_noise** – Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).

#### add_regularly_spaced_functionalized_strands(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_strands: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), strand_length: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], spacing_between_crosslinks: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), offset_to_first_crosslink: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0, crosslinker_functionality: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 4, crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, strand_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1, white_noise: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Add strands with regularly spaced crosslinkers in between.

* **Parameters:**
  * **nr_of_strands** – Number of strands to add.
  * **strand_length** – Length of each strand.
  * **spacing_between_crosslinks** – Number of beads between crosslinks.
  * **offset_to_first_crosslink** – Offset from start of strand to first crosslink (0-based, default: 0).
  * **crosslinker_functionality** – Functionality of the crosslinkers (default: 4).
  * **crosslinker_type** – Atom type of the crosslinkers (default: 2).
  * **strand_atom_type** – Atom type of the beads that stay (default: 1).
  * **white_noise** – Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).

#### add_solvent_chains(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_solvent_chains: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), solvent_chain_length: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), solvent_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 3, white_noise: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Randomly distribute additional, free chains.

* **Parameters:**
  * **nr_of_solvent_chains** – Number of solvent chains to add.
  * **solvent_chain_length** – Length of each solvent chain in beads.
  * **solvent_atom_type** – Atom type for solvent chain beads (default: 3).
  * **white_noise** – Whether to use white noise (true) or blue noise (false) for positioning (default: true).

#### add_star_crosslinkers(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_stars: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), functionality: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), beads_per_strand: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), crosslinker_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, strand_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1, white_noise: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Add star-like crosslinkers with pre-connected strands (useful e.g. for tetra-PEG networks).
Each star crosslinker will have the specified functionality with strands already attached.

* **Parameters:**
  * **nr_of_stars** – Number of star crosslinkers to add
  * **functionality** – Functionality of each star crosslinker (number of strands)
  * **beads_per_strand** – Number of beads in each strand
  * **crosslinker_atom_type** – Atom type for the crosslinker
  * **strand_atom_type** – Atom type for the strand beads
  * **white_noise** – Whether to use white noise positioning

#### add_strands(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), nr_of_strands: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), strand_lengths: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], strand_atom_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 1) → [None](https://docs.python.org/3/library/constants.html#None)

Add strands.
Adds them unconnected at first.
To link them to crosslinkers, use some of the `link_strand_` methods.

* **Parameters:**
  * **nr_of_strands** – Number of strands to add.
  * **strand_lengths** – A list of integers representing the number of beads of each of the strands.
  * **strand_atom_type** – Type of atoms for the strands (default: 1).

#### config_nr_of_mc_steps(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), n_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2000) → [None](https://docs.python.org/3/library/constants.html#None)

Set the number of Monte-Carlo steps during bond length equilibration.

* **Parameters:**
  **n_steps** – Number of MC steps to perform (default: 2000)

#### config_primary_loop_probability(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), probability: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Configure an additional weight reduction for the primary loop formation.

Defaults to 1., which means the general $P(\vec{R})$ is used without any bias.
This results in more primary loops for shorter chains than for longer ones.

Set to 0. to disable the formation of primary loops.

#### config_secondary_loop_probability(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), probability: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Configure an additional weight reduction for the secondary loop formation.

Defaults to 1., which means the general $P(\vec{R})$ is used without any bias.
This results in more secondary loops for shorter chains than for longer ones.

Set to 0. to disable the formation of secondary loops.

#### disable_max_distance(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [None](https://docs.python.org/3/library/constants.html#None)

Disable the maximum distance provider to allow unlimited distance sampling.
This may slow down performance for large systems but ensures complete sampling.

#### get_current_crosslinker_conversion(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the current conversion of crosslinkers, i.e., the fraction of crosslinkers
that have been linked to strands.

* **Returns:**
  Crosslinker conversion as a fraction between 0 and 1.

#### get_current_nr_of_atoms(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the current number of atoms that the universe would/will have.

* **Returns:**
  Number of atoms in the generated universe.

#### get_current_nr_of_bonds(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the current number of bonds that the universe would/will have.

* **Returns:**
  Number of bonds in the generated universe.

#### get_current_strand_conversion(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the current conversion of strands, i.e., the fraction of strand ends that have been consumed.

* **Returns:**
  Strand conversion as a fraction between 0 and 1.

#### get_force_balance(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [pylimer_tools_cpp.MEHPForceBalance](pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance)

Get an instance of the force balance procedure.
This is a useful shorthand e.g. to skip the sampling of beads within in the strands.

* **Returns:**
  Configured MEHPForceBalance instance.

#### get_force_balance2(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [pylimer_tools_cpp.MEHPForceBalance2](pylimer_tools_cpp.MEHPForceBalance2.md#pylimer_tools_cpp.MEHPForceBalance2)

Get an instance of the force balance procedure.
This is a useful shorthand e.g. to skip the sampling of beads within in the strands.

* **Returns:**
  Configured MEHPForceBalance2 instance.

#### get_force_relaxation(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [pylimer_tools_cpp.MEHPForceRelaxation](pylimer_tools_cpp.MEHPForceRelaxation.md#pylimer_tools_cpp.MEHPForceRelaxation)

Get an instance of the force relaxation procedure.
This is a useful shorthand e.g. to skip the sampling of beads within in the strands.

* **Returns:**
  Configured MEHPForceRelaxation instance.

#### get_mean_bead_distance(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the currently configured mean bead distance.

* **Returns:**
  The currently configured mean distance between beads

#### get_mean_squared_bead_distance(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the currently configured mean squared bead distance.

* **Returns:**
  The currently configured mean squared distance between beads

#### get_nr_of_atoms(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the current number of atoms that the universe would/will have.

* **Returns:**
  Number of atoms in the generated universe

#### get_nr_of_bonds(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the current number of bonds that the universe would/will have.

* **Returns:**
  Number of bonds in the generated universe

#### get_universe(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

Fetch the current (or final) state of the universe.

Use this method to actually (MC) place beads between the crosslinks and retrieve the generated structure.

* **Returns:**
  The generated Universe object containing all atoms, bonds, and their coordinates.

#### link_strand(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), strand_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), c_infinity: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Link one particular strand to a crosslinker.
This strand will have one bond made, to an appropriate crosslinker,
as chosen by the parameters associated with the strand.

* **Parameters:**
  * **strand_idx** – Index of the strand to be linked.
  * **c_infinity** – As needed for the end-to-end distribution, given by $\langle R^2\rangle_0 = C_{\infty} N b^2$ (default: 1.0).

#### link_strand_to(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), strand_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), link_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [None](https://docs.python.org/3/library/constants.html#None)

Link a strand to a specific crosslinker.
This assumes that you keep track of the order in which you added the crosslinkers and strands,
as those will determine the indices.

* **Parameters:**
  * **strand_idx** – Index of the strand to be linked.
  * **link_idx** – Index of the crosslinker to be linked.

#### link_strand_to_strand(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), c_infinity: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [bool](https://docs.python.org/3/library/functions.html#bool)

Link two free strand ends together directly.
This method finds free strand ends and combines them into a single strand based on
end-to-end distance probability distributions.

* **Parameters:**
  **c_infinity** – Statistical parameter for end-to-end distance calculation
* **Returns:**
  True if a successful link was made, False if no suitable pair was found

#### link_strands_callback(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), linking_controller: [collections.abc.Callable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Callable)[[[pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], [pylimer_tools_cpp.BackTrackStatus](pylimer_tools_cpp.BackTrackStatus.md#pylimer_tools_cpp.BackTrackStatus)], c_infinity: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Link strands to crosslinkers using a custom callback function to control when to stop.

The callback function receives the current MCUniverseGenerator state and the current step number,
and should return a BackTrackStatus value:
- STOP: Stop the linking process
- TRACK_FORWARD: Continue linking, make more bonds
- TRACK_BACKWARD: Track backward in the linking process, i.e., remove bonds again

* **Parameters:**
  * **linking_controller** – Callback function that controls the linking process.
    Function signature: (MCUniverseGenerator, int) -> BackTrackStatus
  * **c_infinity** – As needed for the end-to-end distribution, given by $\langle R^2\rangle_0 = C_{\infty} N b^2$.

#### link_strands_to_conversion(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), crosslinker_conversion: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), c_infinity: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Actually link the previously added strands to the previously added crosslinkers,
until a certain crosslink conversion is reached.

* **Parameters:**
  * **crosslinker_conversion** – Target conversion of crosslinkers (0: no connections to crosslinks; 1: all crosslinkers fully connected).
  * **c_infinity** – As needed for the end-to-end distribution, given by $\langle R^2\rangle_0 = C_{\infty} N b^2$.

#### link_strands_to_soluble_fraction(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), soluble_fraction: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), c_infinity: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Actually link the previously added strands to the previously added crosslinkers,
until a certain soluble fraction is reached.

* **Parameters:**
  * **soluble_fraction** – Target soluble fraction (0: all material is soluble; 1: no material is soluble).
  * **c_infinity** – As needed for the end-to-end distribution, given by $\langle R^2\rangle_0 = C_{\infty} N b^2$ (default: 1.0).

#### link_strands_to_strands_to_conversion(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), target_strand_conversion: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), c_infinity: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Link free strand ends to each other until target conversion is reached.
This method repeatedly calls [`link_strand_to_strand()`](#pylimer_tools_cpp.MCUniverseGenerator.link_strand_to_strand)
in a slightly more efficient manner than you could from Python,
until the specified conversion of free strand ends is achieved.

#### WARNING
The strands are merged during this process, i.e., from two strands results one strand.
Consequently, calling `get_current_strands_conversion()` after this method will not return
the requested target_strand_conversion, since that one is taken relative to the number of strands when calling this method.

* **Parameters:**
  * **target_strand_conversion** – Target conversion of free strand ends (0.0 to 1.0)
  * **c_infinity** – Statistical parameter for end-to-end distance calculation

#### relax_crosslinks(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [None](https://docs.python.org/3/library/constants.html#None)

Run force relaxation with the crosslinkers and their strands,
to have the crosslinks in their statistically most probable position.

#### remove_soluble_fraction(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), rescale_box: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Remove soluble fraction (as determined by phantom force relaxation) of the strands and crosslinks.

* **Parameters:**
  **rescale_box** – Whether to rescale the box dimensions to maintain constant density (default: true).

#### set_bead_distance(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), distance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), update_mean_squared: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Set the mean distance between beads when doing MC stepping.
Also used for the target crosslinker partner sampling.

NOTE: Mainly the mean squared bead distance is effectively used in the Monte-Carlo simulation.

* **Parameters:**
  * **distance** – Mean distance between beads.
  * **update_mean_squared** – Whether to update the mean squared distance as well, deduced from the assumed gaussian distribution in 3D (default: true).

#### set_mean_squared_bead_distance(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), mean_squared_distance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), update_mean: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Set the mean squared distance between beads.

* **Parameters:**
  * **mean_squared_distance** – Mean squared distance between beads
  * **update_mean** – Whether to update the mean bead distance as well, deduced from the assumed gaussian distribution in 3D (default: true)

#### set_seed(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), seed: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [None](https://docs.python.org/3/library/constants.html#None)

Set the seed for the random generator.

* **Parameters:**
  **seed** – Random seed value for reproducible results.

#### use_linear_max_distance(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), multiplier: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [None](https://docs.python.org/3/library/constants.html#None)

Converts the $N$ to a maximum distance within which to sample.
The distance will be calculated as $N \times \text{multiplier}$.
For example, commonly, a maximum distance is $N \times <b>$,
in which case the multiplier is $<b>$.

Useful only for performance improvements in large systems.

#### use_zscore_max_distance(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator), std_multiplier: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), in_sqrt_multiplier: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0) → [None](https://docs.python.org/3/library/constants.html#None)

Converts the $N$ to a maximum distance within which to sample.
The distance will be calculated as $\text{std_multiplier} \times \sqrt{N \times \text{in_sqrt_multiplier}}$.
Useful only for performance improvements in large systems.

* **Parameters:**
  * **std_multiplier** – The Z-Score, the multiplier of the standard deviation of the end-to-end distribution.
    E.g., 3. for 99.9994% of all conformations.
  * **in_sqrt_multiplier** – The multiplier with the $N$ in the square root. Probably $<b^2>$.

#### validate(self: [pylimer_tools_cpp.MCUniverseGenerator](#pylimer_tools_cpp.MCUniverseGenerator)) → [None](https://docs.python.org/3/library/constants.html#None)

Check whether the internal state of the generator is valid.
Throws errors if not.
Should in principle always be valid when called from Python; if not, there is a bug in the code.
