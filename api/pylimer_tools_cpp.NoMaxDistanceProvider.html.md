# NoMaxDistanceProvider

### *class* pylimer_tools_cpp.NoMaxDistanceProvider(self: [pylimer_tools_cpp.NoMaxDistanceProvider](#pylimer_tools_cpp.NoMaxDistanceProvider))

Bases: [`MaxDistanceProvider`](pylimer_tools_cpp.MaxDistanceProvider.md#pylimer_tools_cpp.MaxDistanceProvider)

For MC generation, to disable the neighbour list usage.

### Methods Summary

| [`get_max_distance`](#pylimer_tools_cpp.NoMaxDistanceProvider.get_max_distance)(self, N)   | Get the maximum distance for a given N (always returns -1 to disable).   |
|--------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|

### Methods Documentation

#### get_max_distance(self: [pylimer_tools_cpp.NoMaxDistanceProvider](#pylimer_tools_cpp.NoMaxDistanceProvider), N: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the maximum distance for a given N (always returns -1 to disable).

* **Parameters:**
  **N** – Number of segments (ignored).
* **Returns:**
  Always returns -1 to disable maximum distance checks.
