# LinearMaxDistanceProvider

### *class* pylimer_tools_cpp.LinearMaxDistanceProvider(self: [pylimer_tools_cpp.LinearMaxDistanceProvider](#pylimer_tools_cpp.LinearMaxDistanceProvider), max_distance_multiplier: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat))

Bases: [`MaxDistanceProvider`](pylimer_tools_cpp.MaxDistanceProvider.md#pylimer_tools_cpp.MaxDistanceProvider)

For MC generation, converts the $N$ to a maximum distance within which to sample.
The distance will be calculated as $N \times \text{max_distance_multiplier}$.
Useful only for performance improvements in large systems.

* **Parameters:**
  **max_distance_multiplier** – Multiplier for the maximum distance.

### Methods Summary

| [`get_max_distance`](#pylimer_tools_cpp.LinearMaxDistanceProvider.get_max_distance)(self, N)   | Get the maximum distance for a given N.   |
|------------------------------------------------------------------------------------------------|-------------------------------------------|

### Methods Documentation

#### get_max_distance(self: [pylimer_tools_cpp.LinearMaxDistanceProvider](#pylimer_tools_cpp.LinearMaxDistanceProvider), N: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the maximum distance for a given N.

* **Parameters:**
  **N** – Number of segments.
