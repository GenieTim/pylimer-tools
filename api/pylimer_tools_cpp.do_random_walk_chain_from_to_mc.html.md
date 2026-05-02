# do_random_walk_chain_from_to_mc

### pylimer_tools_cpp.do_random_walk_chain_from_to_mc(box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box), from_coordinates: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[3, 1]'], to_coordinates: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[3, 1]'], chain_len: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), bead_distance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, mean_squared_bead_distance: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, seed: [str](https://docs.python.org/3/library/stdtypes.html#str) = '', n_iterations: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 10000) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Do a random walk from one point to another.
Then, relax the points in between using a Metropolis-Monte Carlo simulation.

* **Parameters:**
  * **box** – Simulation box for periodic boundary conditions.
  * **from_coordinates** – Starting coordinates as 3D vector.
  * **to_coordinates** – Target coordinates as 3D vector.
  * **chain_len** – Number of beads to place between start and end points.
  * **bead_distance** – Mean distance between consecutive beads (default: 1.0).
  * **mean_squared_bead_distance** – Mean squared distance between consecutive beads (default: 1.0).
  * **seed** – Random seed for reproducibility (default: empty string for no seed).
  * **n_iterations** – Number of Monte Carlo iterations for relaxation (default: 10000).
* **Returns:**
  Coordinates of the relaxed chain as a flat array (i.e., [x1, y1, z1, x2, y2, z2, …]).
