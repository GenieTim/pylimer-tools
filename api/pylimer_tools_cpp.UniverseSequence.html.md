# UniverseSequence

### *class* pylimer_tools_cpp.UniverseSequence(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence))

Bases: `pybind11_object`

This class represents a sequence of Universes, with the Universe’s data
only being read on request. Dump files are read at once in order
to know how many timesteps/universes are available in total
(but the universes’ data is not read on first look through the file).
This, while it can lead to two (or more) reads of the whole file,
is a measure in order to enable low memory useage if needed (i.e. for large dump files).
Use Python’s iterator to have this UniverseSequence only ever retain one universe in memory.
Alternatively, use [`forget_at_index()`](#pylimer_tools_cpp.UniverseSequence.forget_at_index)
to have the UniverseSequence forget about already read universes.

Construct an empty UniverseSequence.

Use initialization methods to populate it with data.

### Methods Summary

| [`at_index`](#pylimer_tools_cpp.UniverseSequence.at_index)(self, index)                                                         | Get the Universe at the given index (as of in the sequence given by the dump file).                                                           |
|---------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------|
| [`compute_distance_autocorrelation_from_to`](#pylimer_tools_cpp.UniverseSequence.compute_distance_autocorrelation_from_to)(...) | Compute the autocorrelation of the dot product of the distance vector from certain to other atoms.                                            |
| [`compute_distance_from_to_atoms`](#pylimer_tools_cpp.UniverseSequence.compute_distance_from_to_atoms)(self, ...[, ...])        | Compute the root square norm of all the (unwrapped!) distances for the given pair of atoms.                                                   |
| [`compute_msd_for_atom_properties`](#pylimer_tools_cpp.UniverseSequence.compute_msd_for_atom_properties)(self, ...[, ...])      | Compute the mean square displacement for atoms using specified property names.                                                                |
| [`compute_msd_for_atoms`](#pylimer_tools_cpp.UniverseSequence.compute_msd_for_atoms)(self, atom_ids[, ...])                     | Compute the mean square displacement for atoms with the specified IDs.                                                                        |
| [`compute_vector_from_to_atoms`](#pylimer_tools_cpp.UniverseSequence.compute_vector_from_to_atoms)(self, ...[, ...])            | Compute the (unwrapped!) distances for the given pair of atoms.                                                                               |
| [`forget_at_index`](#pylimer_tools_cpp.UniverseSequence.forget_at_index)(self, index)                                           | Clear the memory of the Universe at the given index (as of in the sequence given by the dump file).                                           |
| [`get_all`](#pylimer_tools_cpp.UniverseSequence.get_all)(self)                                                                  | Get all universes initialized back in a list.                                                                                                 |
| [`get_length`](#pylimer_tools_cpp.UniverseSequence.get_length)(self)                                                            | Get the number of universes in this sequence.                                                                                                 |
| [`initialize_from_data_sequence`](#pylimer_tools_cpp.UniverseSequence.initialize_from_data_sequence)(self, data_files)          | Reset and initialize the Universes from an ordered list of Lammps data (`write_data`) files.                                                  |
| [`initialize_from_dump_file`](#pylimer_tools_cpp.UniverseSequence.initialize_from_dump_file)(self, ...)                         | Reset and initialize the Universes from a Lammps `dump` output.                                                                               |
| [`next`](#pylimer_tools_cpp.UniverseSequence.next)(self)                                                                        | Get the Universe that's next in the sequence.                                                                                                 |
| [`reset_iterator`](#pylimer_tools_cpp.UniverseSequence.reset_iterator)(self)                                                    | Reset the internal iterator, such that a subsequent call to [`next()`](#pylimer_tools_cpp.UniverseSequence.next) returns the first one again. |
| [`set_data_file_atom_style`](#pylimer_tools_cpp.UniverseSequence.set_data_file_atom_style)(self, atom_styles)                   | Set the format of the data files to be read.                                                                                                  |

### Methods Documentation

#### at_index(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), index: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

Get the Universe at the given index (as of in the sequence given
by the dump file).

* **Parameters:**
  **index** – The index of the universe to retrieve
* **Returns:**
  The Universe at the specified index

#### compute_distance_autocorrelation_from_to(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), atom_ids_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], atom_ids_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], nr_of_origins: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 25, reduce_memory: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float)]

Compute the autocorrelation of the dot product of the distance vector from certain to other atoms.

For example, this can be used to compute Eq. 4.51 from Masao Doi, Introduction to Polymer Physics, p. 74.

#### compute_distance_from_to_atoms(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), atom_ids_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], atom_ids_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], reduce_memory: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Compute the root square norm of all the (unwrapped!) distances for the given pair of atoms.

Can be used to somewhat faster compute e.g. all the end-to-end or bond distances.
Pay attention that the image flags are correct, otherwise, this data may not be useable.

#### compute_msd_for_atom_properties(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), atom_ids: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], x_property: [str](https://docs.python.org/3/library/stdtypes.html#str), y_property: [str](https://docs.python.org/3/library/stdtypes.html#str), z_property: [str](https://docs.python.org/3/library/stdtypes.html#str), nr_of_origins: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 25, reduce_memory: [bool](https://docs.python.org/3/library/functions.html#bool) = False, max_tau: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float)]

Compute the mean square displacement for atoms using specified property names.

* **Parameters:**
  * **atom_ids** – List of atom IDs for which to compute the MSD
  * **x_property** – Name of the x-coordinate property in the dump file (e.g., “x”, “xu”, “xs”)
  * **y_property** – Name of the y-coordinate property in the dump file (e.g., “y”, “yu”, “ys”)
  * **z_property** – Name of the z-coordinate property in the dump file (e.g., “z”, “zu”, “zs”)
  * **nr_of_origins** – Number of time origins to use for averaging. Higher values provide better statistics but increase computation time (default: 25)
  * **reduce_memory** – If True, reduces memory usage by forgetting universes after processing them (default: False)
  * **max_tau** – Maximum time lag (tau) to compute. If -1, computes for all possible tau values. For better statistics, consider setting this to approximately half the sequence length (default: -1)
* **Returns:**
  Dictionary mapping time lag (tau) to mean square displacement values

#### compute_msd_for_atoms(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), atom_ids: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], nr_of_origins: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 25, reduce_memory: [bool](https://docs.python.org/3/library/functions.html#bool) = False, max_tau: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = -1) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float)]

Compute the mean square displacement for atoms with the specified IDs.

* **Parameters:**
  * **atom_ids** – List of atom IDs for which to compute the MSD
  * **nr_of_origins** – Number of time origins to use for averaging. Higher values provide better statistics but increase computation time (default: 25)
  * **reduce_memory** – If True, reduces memory usage by forgetting universes after processing them (default: False)
  * **max_tau** – Maximum time lag (tau) to compute. If -1, computes for all possible tau values. For better statistics, consider setting this to approximately half the sequence length (default: -1)
* **Returns:**
  Dictionary mapping time lag (tau) to mean square displacement values

#### compute_vector_from_to_atoms(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), atom_ids_from: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], atom_ids_to: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)], reduce_memory: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']]

Compute the (unwrapped!) distances for the given pair of atoms.

Can be used to somewhat faster compute e.g. all the end-to-end or bond vectors.
Pay attention that the image flags are correct, otherwise, this data may not be useable.

#### forget_at_index(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), index: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [None](https://docs.python.org/3/library/constants.html#None)

Clear the memory of the Universe at the given index (as of in the
sequence given by the dump file).

#### get_all(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)]

Get all universes initialized back in a list.
For big dump files or lots of data files, this might lead to memory issues.
Use `__iter__()`
or [`at_index()`](#pylimer_tools_cpp.UniverseSequence.at_index)
and [`forget_at_index()`](#pylimer_tools_cpp.UniverseSequence.forget_at_index)
to craft a more memory-efficient retrieval mechanism.

* **Returns:**
  A list of all Universe objects in the sequence

#### get_length(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of universes in this sequence.

#### initialize_from_data_sequence(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), data_files: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[str](https://docs.python.org/3/library/stdtypes.html#str)]) → [None](https://docs.python.org/3/library/constants.html#None)

Reset and initialize the Universes from an ordered list of Lammps data (`write_data`) files.

#### initialize_from_dump_file(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), initial_data_file: [str](https://docs.python.org/3/library/stdtypes.html#str), dump_file: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [None](https://docs.python.org/3/library/constants.html#None)

Reset and initialize the Universes from a Lammps `dump` output.

#### NOTE
If you have not output the id of the atoms in the dump file, they will be assigned sequentially.
If you have not output the type of the atoms in the dump file, they will be set to -1 if they cannot be infered from the data file.

#### next(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence)) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

Get the Universe that’s next in the sequence.

* **Returns:**
  The next Universe in the sequence

#### reset_iterator(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence)) → [None](https://docs.python.org/3/library/constants.html#None)

Reset the internal iterator, such that a subsequent call to
[`next()`](#pylimer_tools_cpp.UniverseSequence.next) returns the first one again.

#### set_data_file_atom_style(self: [pylimer_tools_cpp.UniverseSequence](#pylimer_tools_cpp.UniverseSequence), atom_styles: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.AtomStyle](pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle)]) → [None](https://docs.python.org/3/library/constants.html#None)

Set the format of the data files to be read. See [`AtomStyle`](pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle).
