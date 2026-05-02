# randomly_sample_entanglements

### pylimer_tools_cpp.randomly_sample_entanglements(universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), nr_of_samples: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), upper_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), lower_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0, minimum_nr_of_samples: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0, same_strand_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 3.0, seed: [str](https://docs.python.org/3/library/stdtypes.html#str) = '', crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2, ignore_crosslinks: [bool](https://docs.python.org/3/library/functions.html#bool) = True, filter_dangling_and_soluble: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [pylimer_tools_cpp.AtomPairEntanglements](pylimer_tools_cpp.AtomPairEntanglements.md#pylimer_tools_cpp.AtomPairEntanglements)

Randomly find pairs of atoms that are close together and could be
entanglements

* **Parameters:**
  * **universe** – The universe of atoms from which to sample entanglements from.
  * **nr_of_samples** – The number of pairs of atoms to randomly sample.
  * **upper_cutoff** – The maximum distance between atoms for a pair to be considered a potential entanglement.
  * **lower_cutoff** – The minimum distance between atoms for a pair to be considered a potential entanglement.
  * **minimum_nr_of_samples** – The minimum number of entanglements to be found.
  * **same_strand_cutoff** – The maximum distance between atoms on the same strand for a pair to be considered a potential entanglement.
  * **seed** – A seed for the random number generator.
  * **crosslinker_type** – The type of crosslinker to consider when finding entanglements. Used for the splitting into strands.
  * **ignore_crosslinks** – Whether to ignore crosslinks when finding entanglements. Careful: if you don’t ignore them, the same-strand policy might not work correctly, since each crosslink should actually be associated with more than one strand.
  * **filter_dangling_and_soluble** – Whether to filter out dangling chains and soluble crosslinks when finding entanglements.
    This means, entanglements involving an obviously (1st order) dangling or soluble chain are
