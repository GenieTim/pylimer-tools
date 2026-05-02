# MaxDistanceProvider

### *class* pylimer_tools_cpp.MaxDistanceProvider

Bases: `pybind11_object`

A generic implementation of a class, that shall provide a maximum distance for the MC sampling.

### Methods Summary

| [`get_max_distance`](#pylimer_tools_cpp.MaxDistanceProvider.get_max_distance)(self, N)   | Get the maximum distance for a given N.   |
|------------------------------------------------------------------------------------------|-------------------------------------------|

### Methods Documentation

#### get_max_distance(self: [pylimer_tools_cpp.MaxDistanceProvider](#pylimer_tools_cpp.MaxDistanceProvider), N: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the maximum distance for a given N.

* **Parameters:**
  **N** – Number of segments.
* **Returns:**
  Maximum distance for sampling.
