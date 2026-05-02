# ZScoreMaxDistanceProvider

### *class* pylimer_tools_cpp.ZScoreMaxDistanceProvider(self: [pylimer_tools_cpp.ZScoreMaxDistanceProvider](#pylimer_tools_cpp.ZScoreMaxDistanceProvider), std_multiplier: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), in_sqrt_multiplier: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat))

Bases: [`MaxDistanceProvider`](pylimer_tools_cpp.MaxDistanceProvider.md#pylimer_tools_cpp.MaxDistanceProvider)

For MC generation, converts the $N$ to a maximum distance within which to sample.
The distance will be calculated as $\text{std_multiplier} \times \sqrt{N \times \text{in_sqrt_multiplier}}$.
Useful only for performance improvements in large systems.

* **Parameters:**
  * **std_multiplier** – The Z-Score, the multiplier of the standard deviation of the end-to-end distribution.
    E.g., 3.29 for 99.9% of all conformations.
  * **in_sqrt_multiplier** – The multiplier with the $N$ in the square root. Probably $<b^2>$.

### Methods Summary

| [`get_max_distance`](#pylimer_tools_cpp.ZScoreMaxDistanceProvider.get_max_distance)(self, N)   | Get the maximum distance for a given N.   |
|------------------------------------------------------------------------------------------------|-------------------------------------------|

### Methods Documentation

#### get_max_distance(self: [pylimer_tools_cpp.ZScoreMaxDistanceProvider](#pylimer_tools_cpp.ZScoreMaxDistanceProvider), N: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the maximum distance for a given N.

* **Parameters:**
  **N** – Number of segments.
