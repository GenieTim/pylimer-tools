# peel_topological_leaves

### pylimer_tools_cpp.peel_topological_leaves(crosslinker_universe: [pylimer_tools_cpp.CrosslinkerUniverse](pylimer_tools_cpp.CrosslinkerUniverse.md#pylimer_tools_cpp.CrosslinkerUniverse)) → [list](https://docs.python.org/3/builtins/stdtypes.html#list)[[bool](https://docs.python.org/3/builtins/functions.html#bool)]

Peel 1-core dangling trees and subgraphs iteratively from the CrosslinkerUniverse topology.
Returns a boolean vector indicating for each strand whether it is ignored (pruned).

* **Parameters:**
  **crosslinker_universe** – The CrosslinkerUniverse representation.
* **Returns:**
  Boolean vector of ignored/pruned strands.
