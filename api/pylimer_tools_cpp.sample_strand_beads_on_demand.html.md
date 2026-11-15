# sample_strand_beads_on_demand

### pylimer_tools_cpp.sample_strand_beads_on_demand(box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.html.md#pylimer_tools_cpp.Box), crosslinker_universe: [pylimer_tools_cpp.CrosslinkerUniverse](pylimer_tools_cpp.CrosslinkerUniverse.html.md#pylimer_tools_cpp.CrosslinkerUniverse), strand_idx: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), strand_seed: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), n_mc_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Deterministically realize the internal bead coordinates of a single strand
in a CrosslinkerUniverse using a per-strand PRNG seed.

* **Parameters:**
  * **box** – The simulation box.
  * **crosslinker_universe** – The CrosslinkerUniverse of crosslinkers and strands.
  * **strand_idx** – Index of the strand in crosslinker_universe.
  * **strand_seed** – Deterministic per-strand PRNG seed.
  * **n_mc_steps** – Number of Monte Carlo equilibration steps (default: 0).
* **Returns:**
  Vector of length 3 \* beads_in_strand[strand_idx] with bead coordinates.
