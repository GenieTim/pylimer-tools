# Atom

### *class* pylimer_tools_cpp.Atom(\*args, \*\*kwargs)

Bases: `pybind11_object`

A single bead or atom.

Overloaded function.

1. \_\_init_\_(self: pylimer_tools_cpp.Atom, id: typing.SupportsInt, type: typing.SupportsInt, x: typing.SupportsFloat, y: typing.SupportsFloat, z: typing.SupportsFloat, nx: typing.SupportsInt, ny: typing.SupportsInt, nz: typing.SupportsInt) -> None
   > Construct this atom.
   > * **param id:**
   >   Unique identifier for the atom
   > * **param type:**
   >   Type classification of the atom
   > * **param x:**
   >   X coordinate position
   > * **param y:**
   >   Y coordinate position
   > * **param z:**
   >   Z coordinate position
   > * **param nx:**
   >   Periodic image flag in x direction
   > * **param ny:**
   >   Periodic image flag in y direction
   > * **param nz:**
   >   Periodic image flag in z direction
2. \_\_init_\_(self: pylimer_tools_cpp.Atom, properties: collections.abc.Mapping[str, typing.SupportsFloat]) -> None
   > Construct this atom from a properties dictionary.

   > The dictionary should contain at least the following keys:
   > - “id”: Unique identifier for the atom
   > - “type”: Type classification of the atom
   > - “x”, “y”, “z”: Coordinate positions
   > - “nx”, “ny”, “nz”: Periodic image flags

   > Any additional properties will be stored as extra data.
   > * **param properties:**
   >   Dictionary containing atom properties

### Methods Summary

| [`compute_vector_to`](#pylimer_tools_cpp.Atom.compute_vector_to)(self, to_atom, pbc_box)          | Compute the vector to another atom.                                                            |
|---------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|
| [`distance_to`](#pylimer_tools_cpp.Atom.distance_to)(self, to_atom, pbc_box)                      | Compute the distance to another atom.                                                          |
| [`distance_to_unwrapped`](#pylimer_tools_cpp.Atom.distance_to_unwrapped)(self, arg0, arg1)        | Compute the distance to another atom respecting the periodic image flags.                      |
| [`get_coordinates`](#pylimer_tools_cpp.Atom.get_coordinates)(self)                                | Get the coordinates of this atom as a vector.                                                  |
| [`get_extra_data`](#pylimer_tools_cpp.Atom.get_extra_data)(self)                                  | Get all extra data properties stored with this atom (e.g., charge, dipole, etc.).              |
| [`get_id`](#pylimer_tools_cpp.Atom.get_id)(self)                                                  | Get the ID of the atom.                                                                        |
| [`get_nx`](#pylimer_tools_cpp.Atom.get_nx)(self)                                                  | Get the box image that the atom is in in x direction (also known as ix or nx).                 |
| [`get_ny`](#pylimer_tools_cpp.Atom.get_ny)(self)                                                  | Get the box image that the atom is in in y direction (also known as iy or ny).                 |
| [`get_nz`](#pylimer_tools_cpp.Atom.get_nz)(self)                                                  | Get the box image that the atom is in in z direction (also known as iz or nz).                 |
| [`get_property`](#pylimer_tools_cpp.Atom.get_property)(self, property)                            | Get a specific property value from the extra data.                                             |
| [`get_type`](#pylimer_tools_cpp.Atom.get_type)(self)                                              | Get the type of the atom.                                                                      |
| [`get_unwrapped_coordinates`](#pylimer_tools_cpp.Atom.get_unwrapped_coordinates)(self, arg0)      | Get the unwrapped coordinates of this atom.                                                    |
| [`get_unwrapped_x`](#pylimer_tools_cpp.Atom.get_unwrapped_x)(self, box)                           | Get the unwrapped x coordinate of the atom.                                                    |
| [`get_unwrapped_y`](#pylimer_tools_cpp.Atom.get_unwrapped_y)(self, box)                           | Get the unwrapped y coordinate of the atom.                                                    |
| [`get_unwrapped_z`](#pylimer_tools_cpp.Atom.get_unwrapped_z)(self, box)                           | Get the unwrapped z coordinate of the atom.                                                    |
| [`get_x`](#pylimer_tools_cpp.Atom.get_x)(self)                                                    | Get the x coordinate of the atom.                                                              |
| [`get_y`](#pylimer_tools_cpp.Atom.get_y)(self)                                                    | Get the y coordinate of the atom.                                                              |
| [`get_z`](#pylimer_tools_cpp.Atom.get_z)(self)                                                    | Get the z coordinate of the atom.                                                              |
| [`mean_position_with`](#pylimer_tools_cpp.Atom.mean_position_with)(self, other_atom, pbc_box)     | Compute the mean position between this atom and another atom, considering periodic boundaries. |
| [`mean_position_with_unwrapped`](#pylimer_tools_cpp.Atom.mean_position_with_unwrapped)(self, ...) | Compute the mean position between this atom and another atom using unwrapped coordinates.      |
| [`vector_to_unwrapped`](#pylimer_tools_cpp.Atom.vector_to_unwrapped)(self, arg0, arg1)            | Compute the vector to another atom respecting the periodic image flags.                        |

### Methods Documentation

#### compute_vector_to(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), to_atom: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), pbc_box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the vector to another atom.

* **Parameters:**
  * **to_atom** – The target atom
  * **pbc_box** – The periodic boundary conditions box
* **Returns:**
  Vector pointing from this atom to the target atom

#### distance_to(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), to_atom: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), pbc_box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the distance to another atom.

* **Parameters:**
  * **to_atom** – The target atom
  * **pbc_box** – The periodic boundary conditions box
* **Returns:**
  Euclidean distance between the atoms

#### distance_to_unwrapped(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), arg0: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), arg1: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the distance to another atom respecting the periodic image flags.

* **Parameters:**
  **to_atom** – The target atom
* **Returns:**
  Unwrapped distance to the target atom

#### get_coordinates(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Get the coordinates of this atom as a vector.

* **Returns:**
  A vector containing the x, y, z coordinates

#### get_extra_data(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [float](https://docs.python.org/3/library/functions.html#float)]

Get all extra data properties stored with this atom (e.g., charge, dipole, etc.).

* **Returns:**
  Dictionary containing all extra properties

#### get_id(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the ID of the atom.

* **Returns:**
  The atom’s unique identifier

#### get_nx(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the box image that the atom is in in x direction (also known as ix or nx).

* **Returns:**
  The periodic image flag in x direction

#### get_ny(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the box image that the atom is in in y direction (also known as iy or ny).

* **Returns:**
  The periodic image flag in y direction

#### get_nz(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the box image that the atom is in in z direction (also known as iz or nz).

* **Returns:**
  The periodic image flag in z direction

#### get_property(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), property: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [float](https://docs.python.org/3/library/functions.html#float)

Get a specific property value from the extra data.

* **Parameters:**
  **property** – The name of the property to retrieve
* **Returns:**
  The value of the specified property
* **Raises:**
  std::out_of_range if the property doesn’t exist

#### get_type(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the type of the atom.

* **Returns:**
  The atom’s type classification

#### get_unwrapped_coordinates(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), arg0: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Get the unwrapped coordinates of this atom.

* **Parameters:**
  **box** – The simulation box to use for unwrapping
* **Returns:**
  A vector containing the unwrapped x, y, z coordinates

#### get_unwrapped_x(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the unwrapped x coordinate of the atom.

* **Parameters:**
  **box** – The simulation box to use for unwrapping
* **Returns:**
  The unwrapped x coordinate

#### get_unwrapped_y(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the unwrapped y coordinate of the atom.

* **Parameters:**
  **box** – The simulation box to use for unwrapping
* **Returns:**
  The unwrapped y coordinate

#### get_unwrapped_z(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the unwrapped z coordinate of the atom.

* **Parameters:**
  **box** – The simulation box to use for unwrapping
* **Returns:**
  The unwrapped z coordinate

#### get_x(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the x coordinate of the atom.

* **Returns:**
  The x coordinate

#### get_y(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the y coordinate of the atom.

* **Returns:**
  The y coordinate

#### get_z(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the z coordinate of the atom.

* **Returns:**
  The z coordinate

#### mean_position_with(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), other_atom: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), pbc_box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the mean position between this atom and another atom, considering periodic boundaries.

* **Parameters:**
  * **other_atom** – The other atom
  * **pbc_box** – The periodic boundary conditions box
* **Returns:**
  Vector representing the mean position

#### mean_position_with_unwrapped(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), other_atom: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), pbc_box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the mean position between this atom and another atom using unwrapped coordinates.

* **Parameters:**
  * **other_atom** – The other atom
  * **pbc_box** – The periodic boundary conditions box
* **Returns:**
  Vector representing the mean position (unwrapped)

#### vector_to_unwrapped(self: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), arg0: [pylimer_tools_cpp.Atom](#pylimer_tools_cpp.Atom), arg1: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.md#pylimer_tools_cpp.Box)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Compute the vector to another atom respecting the periodic image flags.

* **Parameters:**
  **to_atom** – The target atom
* **Returns:**
  Unwrapped vector to the target atom
