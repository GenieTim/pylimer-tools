# DumpFileReader

### *class* pylimer_tools_cpp.DumpFileReader(self: [pylimer_tools_cpp.DumpFileReader](#pylimer_tools_cpp.DumpFileReader), path_of_file_to_read: [str](https://docs.python.org/3/library/stdtypes.html#str))

Bases: `pybind11_object`

A reader for LAMMPS’s dump files.

Initialize the dump file reader.

* **Parameters:**
  **path_of_file_to_read** – Path to the dump file to read

### Methods Summary

| [`get_length`](#pylimer_tools_cpp.DumpFileReader.get_length)(self)                                              | Get the number of sections (time-steps) in the file.                                                     |
|-----------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| [`get_numeric_values_for_at`](#pylimer_tools_cpp.DumpFileReader.get_numeric_values_for_at)(self, rowIndex, ...) | Get numeric values for the section at index, the main header headerKey and the column columnIndex.       |
| [`get_string_values_for_at`](#pylimer_tools_cpp.DumpFileReader.get_string_values_for_at)(self, rowIndex, ...)   | Get string values for the section at index, the main header headerKey and the column columnIndex.        |
| [`has_key`](#pylimer_tools_cpp.DumpFileReader.has_key)(self, headerKey)                                         | Check whether the first section has the header specified.                                                |
| [`key_has_column`](#pylimer_tools_cpp.DumpFileReader.key_has_column)(self, headerKey, columnName)               | Check whether the header of the first section has the specified column.                                  |
| [`key_has_directional_column`](#pylimer_tools_cpp.DumpFileReader.key_has_directional_column)(self, header_key)  | Check whether the header of the first section has all the three columns {dir_prefix}{x|y|z}{dir_suffix}. |
| [`read`](#pylimer_tools_cpp.DumpFileReader.read)(self)                                                          | Read the whole file.                                                                                     |

### Methods Documentation

#### get_length(self: [pylimer_tools_cpp.DumpFileReader](#pylimer_tools_cpp.DumpFileReader)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of sections (time-steps) in the file.

* **Returns:**
  Number of timesteps/sections in the dump file

#### get_numeric_values_for_at(self: [pylimer_tools_cpp.DumpFileReader](#pylimer_tools_cpp.DumpFileReader), rowIndex: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), headerKey: [str](https://docs.python.org/3/library/stdtypes.html#str), columnIndex: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Get numeric values for the section at index, the main header headerKey and the column columnIndex.

* **Parameters:**
  * **rowIndex** – Index of the section/timestep to query
  * **headerKey** – Name of the header section to query
  * **columnIndex** – Index of the column within the header
* **Returns:**
  Vector of numeric values for the specified location

#### get_string_values_for_at(self: [pylimer_tools_cpp.DumpFileReader](#pylimer_tools_cpp.DumpFileReader), rowIndex: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), headerKey: [str](https://docs.python.org/3/library/stdtypes.html#str), columnIndex: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)]

Get string values for the section at index, the main header headerKey and the column columnIndex.

* **Parameters:**
  * **rowIndex** – Index of the section/timestep to query
  * **headerKey** – Name of the header section to query
  * **columnIndex** – Index of the column within the header
* **Returns:**
  Vector of string values for the specified location

#### has_key(self: [pylimer_tools_cpp.DumpFileReader](#pylimer_tools_cpp.DumpFileReader), headerKey: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [bool](https://docs.python.org/3/library/functions.html#bool)

Check whether the first section has the header specified.

* **Parameters:**
  **headerKey** – Name of the header to check for
* **Returns:**
  True if the header exists in the first section

#### key_has_column(self: [pylimer_tools_cpp.DumpFileReader](#pylimer_tools_cpp.DumpFileReader), headerKey: [str](https://docs.python.org/3/library/stdtypes.html#str), columnName: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [bool](https://docs.python.org/3/library/functions.html#bool)

Check whether the header of the first section has the specified column.

* **Parameters:**
  * **headerKey** – Name of the header section to check
  * **columnName** – Name of the column to look for
* **Returns:**
  True if the column exists in the specified header

#### key_has_directional_column(self: [pylimer_tools_cpp.DumpFileReader](#pylimer_tools_cpp.DumpFileReader), header_key: [str](https://docs.python.org/3/library/stdtypes.html#str), dir_prefix: [str](https://docs.python.org/3/library/stdtypes.html#str) = '', dir_suffix: [str](https://docs.python.org/3/library/stdtypes.html#str) = '') → [bool](https://docs.python.org/3/library/functions.html#bool)

Check whether the header of the first section has all the three columns {dir_prefix}{x|y|z}{dir_suffix}.

* **Parameters:**
  * **header_key** – Name of the header section to check
  * **dir_prefix** – Prefix for the directional columns (default: “”)
  * **dir_suffix** – Suffix for the directional columns (default: “”)
* **Returns:**
  True if all three directional columns (x, y, z) exist

#### read(self: [pylimer_tools_cpp.DumpFileReader](#pylimer_tools_cpp.DumpFileReader)) → [None](https://docs.python.org/3/library/constants.html#None)

Read the whole file.
