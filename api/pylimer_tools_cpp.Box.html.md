# Box

### *class* pylimer_tools_cpp.Box(\*args, \*\*kwargs)

Bases: `pybind11_object`

The box that the simulation is run in.

#### NOTE
Currently, only rectangular boxes are supported.

Overloaded function.

1. \_\_init_\_(self: pylimer_tools_cpp.Box, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat) -> None
2. \_\_init_\_(self: pylimer_tools_cpp.Box, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat, arg3: typing.SupportsFloat, arg4: typing.SupportsFloat, arg5: typing.SupportsFloat) -> None

### Methods Summary

| [`apply_pbc`](#pylimer_tools_cpp.Box.apply_pbc)(self, distances)                                | Apply periodic boundary conditions (PBC): adjust the specified distances to fit into this box.   |
|-------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------|
| [`apply_simple_shear`](#pylimer_tools_cpp.Box.apply_simple_shear)(self, shear_magnitude[, ...]) | Apply a simple shear to the box.                                                                 |
| [`get_bounding_box`](#pylimer_tools_cpp.Box.get_bounding_box)(self)                             | Get an orthogonal box that encloses this box.                                                    |
| [`get_high_x`](#pylimer_tools_cpp.Box.get_high_x)(self)                                         | Get the upper bound of the box in x direction.                                                   |
| [`get_high_y`](#pylimer_tools_cpp.Box.get_high_y)(self)                                         | Get the upper bound of the box in y direction.                                                   |
| [`get_high_z`](#pylimer_tools_cpp.Box.get_high_z)(self)                                         | Get the upper bound of the box in z direction.                                                   |
| [`get_l`](#pylimer_tools_cpp.Box.get_l)(self)                                                   | Get the three lengths of the box in an array/list.                                               |
| [`get_low_x`](#pylimer_tools_cpp.Box.get_low_x)(self)                                           | Get the lower bound of the box in x direction.                                                   |
| [`get_low_y`](#pylimer_tools_cpp.Box.get_low_y)(self)                                           | Get the lower bound of the box in y direction.                                                   |
| [`get_low_z`](#pylimer_tools_cpp.Box.get_low_z)(self)                                           | Get the lower bound of the box in z direction.                                                   |
| [`get_lx`](#pylimer_tools_cpp.Box.get_lx)(self)                                                 | Get the length of the box in x direction.                                                        |
| [`get_ly`](#pylimer_tools_cpp.Box.get_ly)(self)                                                 | Get the length of the box in y direction.                                                        |
| [`get_lz`](#pylimer_tools_cpp.Box.get_lz)(self)                                                 | Get the length of the box in z direction.                                                        |
| [`get_offset`](#pylimer_tools_cpp.Box.get_offset)(self, distances)                              | Compute the offset required to compensate for periodic boundary conditions.                      |
| [`get_volume`](#pylimer_tools_cpp.Box.get_volume)(self)                                         | Compute the volume of the box.                                                                   |
| [`is_valid_offset`](#pylimer_tools_cpp.Box.is_valid_offset)(self, potential_offset[, ...])      | Check whether the passed offset is a valid one in this box.                                      |

### Methods Documentation

#### apply_pbc(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box), distances: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Apply periodic boundary conditions (PBC): adjust the specified distances to fit into this box.

* **Parameters:**
  **distances** – The distances to adjust
* **Returns:**
  The adjusted distances

#### apply_simple_shear(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box), shear_magnitude: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), shear_direction: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 0) → [None](https://docs.python.org/3/library/constants.html#None)

Apply a simple shear to the box.

#### WARNING
Currently, this is not supported for all operations.

For shear magnitude, you specify the angle $\gamma$.

* **Parameters:**
  * **shear_magnitude** – The shear magnitude (angle $\gamma$)
  * **shear_direction** – Direction of shear: 0 for x, 1 for y, 2 for z.
    Use any other integer to disable shear.

#### get_bounding_box(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)

Get an orthogonal box that encloses this box.

For non-sheared boxes, the resulting box is identical to the current box.

* **Returns:**
  A new Box object representing the bounding box

#### get_high_x(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the upper bound of the box in x direction.

* **Returns:**
  The upper x-coordinate boundary

#### get_high_y(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the upper bound of the box in y direction.

* **Returns:**
  The upper y-coordinate boundary

#### get_high_z(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the upper bound of the box in z direction.

* **Returns:**
  The upper z-coordinate boundary

#### get_l(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Get the three lengths of the box in an array/list.

* **Returns:**
  Array containing [Lx, Ly, Lz] box dimensions

#### get_low_x(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the lower bound of the box in x direction.

* **Returns:**
  The lower x-coordinate boundary

#### get_low_y(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the lower bound of the box in y direction.

* **Returns:**
  The lower y-coordinate boundary

#### get_low_z(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the lower bound of the box in z direction.

* **Returns:**
  The lower z-coordinate boundary

#### get_lx(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the length of the box in x direction.

* **Returns:**
  The x-dimension length of the box

#### get_ly(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the length of the box in y direction.

* **Returns:**
  The y-dimension length of the box

#### get_lz(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the length of the box in z direction.

* **Returns:**
  The z-dimension length of the box

#### get_offset(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box), distances: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Compute the offset required to compensate for periodic boundary conditions.

Useful e.g. if you are using absolute coordinates for distances, but
still need an infinite network,
e.g., if the bonds need to be able to get longer than half the box.

* **Parameters:**
  **distances** – The distances to compute offset for
* **Returns:**
  The computed offset

#### get_volume(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the volume of the box.

$V = L_x \cdot L_y \cdot L_z$

* **Returns:**
  The volume of the box

#### is_valid_offset(self: [pylimer_tools_cpp.Box](#pylimer_tools_cpp.Box), potential_offset: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], abs_precision: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) = 1e-05) → [bool](https://docs.python.org/3/library/functions.html#bool)

Check whether the passed offset is a valid one in this box.

* **Parameters:**
  * **potential_offset** – The offset to validate
  * **abs_precision** – Absolute precision for the validation
* **Returns:**
  True if the offset is valid, False otherwise
