# randomly_sample_entanglements_analytical

### pylimer_tools_cpp.randomly_sample_entanglements_analytical(\*args, \*\*kwargs)

Overloaded function.

1. randomly_sample_entanglements_analytical(box: pylimer_tools_cpp.Box, crosslinker_universe: pylimer_tools_cpp.CrosslinkerUniverse, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = ‘’, filtering_mode: pylimer_tools_cpp.EntanglementFilteringMode = EntanglementFilteringMode.OMIT_INACTIVE) -> pylimer_tools_cpp.StrandPairEntanglements
   > Analytically and probabilistically sample pairs of strand beads that could form entanglements
   > using exact 3D Brownian bridge / random walk statistics and topological strand filtering,
   > without requiring all strand beads or a Universe to be instantiated.
   > * **param box:**
   >   The simulation box for periodic boundary conditions.
   > * **param crosslinker_universe:**
   >   The CrosslinkerUniverse of crosslinkers and strands.
   > * **param nr_of_samples:**
   >   The number of entanglement pairs to randomly sample.
   > * **param upper_cutoff:**
   >   The maximum distance between two beads to be considered an entanglement.
   > * **param lower_cutoff:**
   >   The minimum distance between two beads to be considered an entanglement.
   > * **param minimum_nr_of_samples:**
   >   The minimum number of entanglements to be found.
   > * **param same_strand_cutoff:**
   >   The minimum contour separation (in beads) between two beads on the same strand.
   > * **param seed:**
   >   A seed for the random number generator.
   > * **param filtering_mode:**
   >   Filtering mode for entanglements (default: EntanglementFilteringMode.OMIT_INACTIVE).
2. randomly_sample_entanglements_analytical(crosslinker_universe: pylimer_tools_cpp.CrosslinkerUniverse, box: pylimer_tools_cpp.Box, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = ‘’, filtering_mode: pylimer_tools_cpp.EntanglementFilteringMode = EntanglementFilteringMode.OMIT_INACTIVE) -> pylimer_tools_cpp.StrandPairEntanglements
