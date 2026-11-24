# peel_topological_leaves

### pylimer_tools_cpp.peel_topological_leaves(crosslinker_universe: [pylimer_tools_cpp.CrosslinkerUniverse](pylimer_tools_cpp.CrosslinkerUniverse.md#pylimer_tools_cpp.CrosslinkerUniverse)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[bool](https://docs.python.org/3/builtins/functions.html#bool)]

Peel 1-core dangling trees and subgraphs iteratively from the CrosslinkerUniverse topology.
Returns a boolean vector indicating for each strand whether it is ignored (pruned).

#### NOTE
This is a purely topological 2-core peeling on `strand_from`/`strand_to` and
does not account for coordinates or periodic boundary conditions (PBC):

- Free/dangling strands (including crosslinker-less solvent rings) are always pruned,
  even if they wrap across the PBC.
- Primary loops (`strand_from == strand_to >= 0`) and multi-crosslinker cycles are
  always kept in the 2-core, even if they do not span the PBC (note that a
  non-PBC primary loop can still become elastically active if an entanglement or
  slip-link is associated with it, whereas in a phantom network without
  entanglements it relaxes to zero tension).

To distinguish PBC-spanning (or entangled) active loops from non-PBC
zero-tension loops, use [`MEHPSimulator`](pylimer_tools_cpp.MEHPSimulator.md#pylimer_tools_cpp.MEHPSimulator) with
`assume_box_large_enough=False`.

* **Parameters:**
  **crosslinker_universe** – The CrosslinkerUniverse representation.
* **Returns:**
  Boolean vector of ignored/pruned strands.
