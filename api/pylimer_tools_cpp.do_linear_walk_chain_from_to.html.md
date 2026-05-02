# do_linear_walk_chain_from_to

### pylimer_tools_cpp.do_linear_walk_chain_from_to(box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box), from_coordinates: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[3, 1]'], to_coordinates: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[3, 1]'], chain_len: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), include_ends: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get coordinates linearly interpolated from one point to another (both exclusive).

* **Parameters:**
  * **box** – The box for doing PBC correction on the from/to.
  * **from_coordinates** – Coordinates of the start point.
  * **to_coordinates** – Coordinates of the end point.
  * **chain_len** – Number of coordinates to generate between the start and end-point.
  * **include_ends** – Whether to include the start and end points in the output (default: false).
    If yes, chain_len + 2 coordinates will be returned,
    where the first will be from_coordinates and the last will be to_coordinates.
