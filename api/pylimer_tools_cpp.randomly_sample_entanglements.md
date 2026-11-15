# randomly_sample_entanglements

### pylimer_tools_cpp.randomly_sample_entanglements(\*args, \*\*kwargs)

Overloaded function.

1. randomly_sample_entanglements(universe: pylimer_tools_cpp.Universe, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = ‘’, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, ignore_crosslinks: bool = True, filtering_mode: pylimer_tools_cpp.EntanglementFilteringMode = EntanglementFilteringMode.NONE) -> pylimer_tools_cpp.AtomPairEntanglements
   > Randomly find pairs of atoms that are close together and could be entanglements.
   > * **param universe:**
   >   The universe of atoms from which to sample entanglements from.
   > * **param nr_of_samples:**
   >   The number of pairs of atoms to randomly sample.
   > * **param upper_cutoff:**
   >   The maximum distance between atoms for a pair to be considered a potential entanglement.
   > * **param lower_cutoff:**
   >   The minimum distance between atoms for a pair to be considered a potential entanglement.
   > * **param minimum_nr_of_samples:**
   >   The minimum number of entanglements to be found.
   > * **param same_strand_cutoff:**
   >   The minimum shortest path distance on the graph between atoms for an entanglement.
   > * **param seed:**
   >   A seed for the random number generator.
   > * **param crosslinker_type:**
   >   Atom type ID representing crosslinkers.
   > * **param ignore_crosslinks:**
   >   Whether to exclude crosslink atoms from forming entanglements.
   > * **param filtering_mode:**
   >   Filtering mode for entanglements (default: EntanglementFilteringMode.NONE).
2. randomly_sample_entanglements(box: pylimer_tools_cpp.Box, crosslinker_universe: pylimer_tools_cpp.CrosslinkerUniverse, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = ‘’, filtering_mode: pylimer_tools_cpp.EntanglementFilteringMode = EntanglementFilteringMode.OMIT_INACTIVE, n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> pylimer_tools_cpp.StrandPairEntanglements
3. randomly_sample_entanglements(crosslinker_universe: pylimer_tools_cpp.CrosslinkerUniverse, box: pylimer_tools_cpp.Box, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = ‘’, filtering_mode: pylimer_tools_cpp.EntanglementFilteringMode = EntanglementFilteringMode.OMIT_INACTIVE, n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> pylimer_tools_cpp.StrandPairEntanglements
