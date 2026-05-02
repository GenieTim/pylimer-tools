# do_random_walk

### pylimer_tools_cpp.do_random_walk(chain_len: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), bead_distance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, mean_squared_bead_distance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, seed: [str](https://docs.python.org/3/library/stdtypes.html#str) = '') → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Do a random walk, return the coordinates of each point visited.

* **Parameters:**
  * **chain_len** – Length of the chain to generate.
  * **bead_distance** – Mean distance between consecutive beads (default: 1.0).
  * **mean_squared_bead_distance** – Mean squared distance between consecutive beads (default: 1.0).
  * **seed** – Random seed for reproducibility (default: empty string for random seed).
* **Returns:**
  Coordinates of each point as a flat array (i.e., [x1, y1, z1, x2, y2, z2, …]).
