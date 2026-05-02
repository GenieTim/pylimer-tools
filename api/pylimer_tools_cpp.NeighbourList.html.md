# NeighbourList

### *class* pylimer_tools_cpp.NeighbourList(self: [pylimer_tools_cpp.NeighbourList](#pylimer_tools_cpp.NeighbourList), atoms: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)], box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box), cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat))

Bases: `pybind11_object`

Gives access to somewhat fast queries on the neighbourhood of atoms.

This class provides efficient spatial queries for finding atoms within
a specified distance of each other.

Instantiate a new neighbour list.

* **Parameters:**
  * **atoms** – Vector of atoms to include in the neighbour list
  * **box** – The simulation box
  * **cutoff** – Maximum distance for neighbour searches

### Methods Summary

| [`get_atoms_close_to`](#pylimer_tools_cpp.NeighbourList.get_atoms_close_to)(self, atom[, ...])   | List all atoms that are close to a given one.   |
|--------------------------------------------------------------------------------------------------|-------------------------------------------------|
| [`remove_atom`](#pylimer_tools_cpp.NeighbourList.remove_atom)(self, atom[, debug_hint])          | Remove an atom from this neighbour list.        |

### Methods Documentation

#### get_atoms_close_to(self: [pylimer_tools_cpp.NeighbourList](#pylimer_tools_cpp.NeighbourList), atom: [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom), upper_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1.0, lower_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 0.0, unwrapped: [bool](https://docs.python.org/3/library/functions.html#bool) = False, expect_self: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom)]

List all atoms that are close to a given one.

It is possible to request it within a new cutoff,
though the underlying neighbour list will not be regenerated.
For performance reasons, it is recommended to initialize a
new NeighbourList if you require a different cutoff, depending on your use case.

You can use a negative value for the upper_cutoff to use the cutoff used for
filling the neighbour list buckets.

* **Parameters:**
  * **atom** – The reference atom
  * **upper_cutoff** – Maximum distance for neighbours
  * **lower_cutoff** – Minimum distance for neighbours
  * **unwrapped** – Whether to use unwrapped coordinates
  * **expect_self** – Whether to expect the atom itself in results
* **Returns:**
  List of neighbouring atoms

#### remove_atom(self: [pylimer_tools_cpp.NeighbourList](#pylimer_tools_cpp.NeighbourList), atom: [pylimer_tools_cpp.Atom](pylimer_tools_cpp.Atom.md#pylimer_tools_cpp.Atom), debug_hint: [str](https://docs.python.org/3/library/stdtypes.html#str) = '') → [None](https://docs.python.org/3/library/constants.html#None)

Remove an atom from this neighbour list.
It will not show up when querying for neighbours,
but its neighbours cannot be queried either.

* **Parameters:**
  * **atom** – The atom to remove
  * **debug_hint** – Optional debug information
