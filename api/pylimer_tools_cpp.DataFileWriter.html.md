# DataFileWriter

### *class* pylimer_tools_cpp.DataFileWriter(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe))

Bases: `pybind11_object`

A class to write a LAMMPS data file from a universe.

#### ATTENTION
The resulting file is not guaranteed to be a completely valid LAMMPS file.
In particular, this writer does not force you to set masses for all atom types,
and it does not have limits on the image flags.
You can use [`set_masses()`](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.set_masses) to set the masses for all atom types,
and either [`set_vertex_property()`](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.set_vertex_property) or
[`config_attempt_image_reset()`](#pylimer_tools_cpp.DataFileWriter.config_attempt_image_reset) and
[`config_move_into_box()`](#pylimer_tools_cpp.DataFileWriter.config_move_into_box) to ensure that the image flags are correct.

Initialize the writer with the universe to write.

* **Parameters:**
  **universe** – The universe to write to the data file

### Methods Summary

| [`config_atom_style`](#pylimer_tools_cpp.DataFileWriter.config_atom_style)(self[, atom_style])                    | Set the (LAMMPS) atom style to use for writing the atoms.                                                                                    |
|-------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| [`config_attempt_image_reset`](#pylimer_tools_cpp.DataFileWriter.config_attempt_image_reset)(self[, ...])         | Set whether to attempt to reset image flags so that output coordinates lie in the box.                                                       |
| [`config_crosslinker_type`](#pylimer_tools_cpp.DataFileWriter.config_crosslinker_type)(self[, crosslinker_type])  | Set which atom type represents crosslinkers.                                                                                                 |
| [`config_include_angles`](#pylimer_tools_cpp.DataFileWriter.config_include_angles)(self[, include_angles])        | Set whether to include the angles from the universe in the file or not.                                                                      |
| [`config_include_dihedral_angles`](#pylimer_tools_cpp.DataFileWriter.config_include_dihedral_angles)(self[, ...]) | Set whether to include the dihedral angles from the universe in the file or not.                                                             |
| [`config_include_velocities`](#pylimer_tools_cpp.DataFileWriter.config_include_velocities)(self[, ...])           | Set whether to include the velocities from the universe (if any) in the file or not.                                                         |
| [`config_molecule_idx_for_swap`](#pylimer_tools_cpp.DataFileWriter.config_molecule_idx_for_swap)(self[, ...])     | Swappable chains implies that their moleculeIdx in the LAMMPS data file is not identical per chain, but identical per position in the chain. |
| [`config_move_into_box`](#pylimer_tools_cpp.DataFileWriter.config_move_into_box)(self[, move_into_box])           | Set whether to change the output coordinates to lie in the box or not.                                                                       |
| [`config_reindex_atoms`](#pylimer_tools_cpp.DataFileWriter.config_reindex_atoms)(self[, reindex_atoms])           | Set whether to reindex the atoms or not.                                                                                                     |
| [`set_custom_atom_format`](#pylimer_tools_cpp.DataFileWriter.set_custom_atom_format)(self[, atom_format])         | Specify a custom format for the atom section.                                                                                                |
| [`set_universe_to_write`](#pylimer_tools_cpp.DataFileWriter.set_universe_to_write)(self, universe)                | Re-set the universe to write.                                                                                                                |
| [`write_to_file`](#pylimer_tools_cpp.DataFileWriter.write_to_file)(self, file)                                    | Actually do the writing to the disk.                                                                                                         |

### Methods Documentation

#### config_atom_style(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), atom_style: [pylimer_tools_cpp.AtomStyle](pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle) = AtomStyle.ANGLE) → [None](https://docs.python.org/3/library/constants.html#None)

Set the (LAMMPS) atom style to use for writing the atoms.

* **Parameters:**
  **atom_style** – The LAMMPS atom style to use (default: AtomStyle.ANGLE)

#### config_attempt_image_reset(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), attempt_image_reset: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Set whether to attempt to reset image flags so that output coordinates lie in the box.

* **Parameters:**
  **attempt_image_reset** – Whether to attempt image flag reset (default: False)

#### config_crosslinker_type(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) = 2) → [None](https://docs.python.org/3/library/constants.html#None)

Set which atom type represents crosslinkers.
Needed in case the moleculeIdx in the output file should have any meaning.
(e.g. with [`config_molecule_idx_for_swap()`](#pylimer_tools_cpp.DataFileWriter.config_molecule_idx_for_swap)).

* **Parameters:**
  **crosslinker_type** – The atom type representing crosslinkers (default: 2)

#### config_include_angles(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), include_angles: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Set whether to include the angles from the universe in the file or not.

* **Parameters:**
  **include_angles** – Whether to include angles (default: True)

#### config_include_dihedral_angles(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), include_dihedral_angles: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Set whether to include the dihedral angles from the universe in the file or not.

* **Parameters:**
  **include_dihedral_angles** – Whether to include dihedral angles (default: True)

#### config_include_velocities(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), include_velocities: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [None](https://docs.python.org/3/library/constants.html#None)

Set whether to include the velocities from the universe (if any) in the file or not.

* **Parameters:**
  **include_velocities** – Whether to include velocities (default: True)

#### config_molecule_idx_for_swap(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), enable_swappability: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Swappable chains implies that their moleculeIdx in the LAMMPS data file is not
identical per chain, but identical per position in the chain.
That’s how you can have bond swapping with constant chain length distribution.

* **Parameters:**
  **enable_swappability** – Whether to enable molecule index swappability (default: False)

#### config_move_into_box(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), move_into_box: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Set whether to change the output coordinates to lie in the box or not.

* **Parameters:**
  **move_into_box** – Whether to move coordinates into box (default: False)

#### config_reindex_atoms(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), reindex_atoms: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [None](https://docs.python.org/3/library/constants.html#None)

Set whether to reindex the atoms or not.
Re-indexing leads to atom IDs being in the range of 1 to the number of atoms.

* **Parameters:**
  **reindex_atoms** – Whether to reindex atoms (default: False)

#### set_custom_atom_format(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), atom_format: [str](https://docs.python.org/3/library/stdtypes.html#str) = '\\t$atomId\\t$moleculeId\\t$atomType\\t$x\\t$y\\t$z\\t$nx\\t$ny\\t$nz') → [None](https://docs.python.org/3/library/constants.html#None)

Specify a custom format for the atom section.

Placeholder options are:

- $atomId
- $moleculeId
- $atomType
- $x
- $y
- $z
- $nx
- $ny
- $nz

Additionally, you can use the keys used in
`set_property_value()`
as placeholders (as long as they are alphanumeric only; prefix in the format with ‘$’ as well).
Other placeholders are available if the universe was read from a LAMMPS data file with an
atom style with additional data.

This method is specifically useful if you need a different (or hybrid) atom style in LAMMPS.

Be sure to still call [`config_atom_style()`](#pylimer_tools_cpp.DataFileWriter.config_atom_style),
so that the file can be read correctly again.

* **Parameters:**
  **atom_format** – Custom format string for atoms (default: tab-separated standard format)

#### set_universe_to_write(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)) → [None](https://docs.python.org/3/library/constants.html#None)

Re-set the universe to write.

* **Parameters:**
  **universe** – The new universe to write to the data file

#### write_to_file(self: [pylimer_tools_cpp.DataFileWriter](#pylimer_tools_cpp.DataFileWriter), file: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [None](https://docs.python.org/3/library/constants.html#None)

Actually do the writing to the disk.

* **Parameters:**
  **file** – The path and file name to write to
