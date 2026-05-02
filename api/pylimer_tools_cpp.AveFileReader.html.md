# AveFileReader

### *class* pylimer_tools_cpp.AveFileReader(self: [pylimer_tools_cpp.AveFileReader](#pylimer_tools_cpp.AveFileReader), file_path: [str](https://docs.python.org/3/library/stdtypes.html#str))

Bases: `pybind11_object`

Alternative implementation of the data file reader implemented in
[`pylimer_tools.io.read_lammps_output_file.read_averages_file()`](../pylimer_tools.io.md#pylimer_tools.io.read_lammps_output_file.read_averages_file).

This implementation is better for certain use cases, worse for others.
In the end, only performance and memory usage are different.
For moderately sized and small files, we recommend to use the Python interface instead.

Initialize the AveFileReader with a file path.

* **Parameters:**
  **file_path** – Path to the averages file to read

### Methods Summary

| [`autocorrelate_column`](#pylimer_tools_cpp.AveFileReader.autocorrelate_column)(self, column_index, ...)         | Do autocorrelation on one particular column for a specified set of delta indices.                         |
|------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| [`autocorrelate_column_difference`](#pylimer_tools_cpp.AveFileReader.autocorrelate_column_difference)(self, ...) | Do autocorrelation on the difference between two particular columns for a specified set of delta indices. |
| [`get_column_names`](#pylimer_tools_cpp.AveFileReader.get_column_names)(self)                                    | Get the names of all columns in the file.                                                                 |
| [`get_data`](#pylimer_tools_cpp.AveFileReader.get_data)(self)                                                    | Get all data from the file.                                                                               |
| [`get_nr_of_columns`](#pylimer_tools_cpp.AveFileReader.get_nr_of_columns)(self)                                  | Get the number of columns in the file.                                                                    |
| [`get_nr_of_rows`](#pylimer_tools_cpp.AveFileReader.get_nr_of_rows)(self)                                        | Get the number of data rows in the file.                                                                  |

### Methods Documentation

#### autocorrelate_column(self: [pylimer_tools_cpp.AveFileReader](#pylimer_tools_cpp.AveFileReader), column_index: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), delta_indices: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Do autocorrelation on one particular column for a specified set of delta indices.

Assumes the data is equally spaced.

* **Parameters:**
  * **column_index** – Index of the column to autocorrelate
  * **delta_indices** – List of delta indices for the autocorrelation
* **Returns:**
  Autocorrelation values for the specified deltas

#### autocorrelate_column_difference(self: [pylimer_tools_cpp.AveFileReader](#pylimer_tools_cpp.AveFileReader), column_index1: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), column_index2: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), delta_indices: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)]) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Do autocorrelation on the difference between two particular columns for a specified set of delta indices.

Assumes the data is equally spaced.

* **Parameters:**
  * **column_index1** – Index of the first column
  * **column_index2** – Index of the second column
  * **delta_indices** – List of delta indices for the autocorrelation
* **Returns:**
  Autocorrelation values for the column differences at specified deltas

#### get_column_names(self: [pylimer_tools_cpp.AveFileReader](#pylimer_tools_cpp.AveFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)]

Get the names of all columns in the file.

* **Returns:**
  List of column names

#### get_data(self: [pylimer_tools_cpp.AveFileReader](#pylimer_tools_cpp.AveFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]]

Get all data from the file.

* **Returns:**
  2D array containing all the numerical data

#### get_nr_of_columns(self: [pylimer_tools_cpp.AveFileReader](#pylimer_tools_cpp.AveFileReader)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of columns in the file.

* **Returns:**
  Number of columns

#### get_nr_of_rows(self: [pylimer_tools_cpp.AveFileReader](#pylimer_tools_cpp.AveFileReader)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of data rows in the file.

* **Returns:**
  Number of rows
