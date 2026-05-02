# DataFileReader

### *class* pylimer_tools_cpp.DataFileReader(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader))

Bases: `pybind11_object`

A reader for LAMMPS’s write_data files.

Initialize a new data file parser.

### Methods Summary

| [`get_atom_ids`](#pylimer_tools_cpp.DataFileReader.get_atom_ids)(self)                 | Get all atom IDs from the data file.                   |
|----------------------------------------------------------------------------------------|--------------------------------------------------------|
| [`get_atom_nx`](#pylimer_tools_cpp.DataFileReader.get_atom_nx)(self)                   | Get periodic image flags in x direction for all atoms. |
| [`get_atom_ny`](#pylimer_tools_cpp.DataFileReader.get_atom_ny)(self)                   | Get periodic image flags in y direction for all atoms. |
| [`get_atom_nz`](#pylimer_tools_cpp.DataFileReader.get_atom_nz)(self)                   | Get periodic image flags in z direction for all atoms. |
| [`get_atom_types`](#pylimer_tools_cpp.DataFileReader.get_atom_types)(self)             | Get all atom types from the data file.                 |
| [`get_atom_x`](#pylimer_tools_cpp.DataFileReader.get_atom_x)(self)                     | Get x coordinates of all atoms.                        |
| [`get_atom_y`](#pylimer_tools_cpp.DataFileReader.get_atom_y)(self)                     | Get y coordinates of all atoms.                        |
| [`get_atom_z`](#pylimer_tools_cpp.DataFileReader.get_atom_z)(self)                     | Get z coordinates of all atoms.                        |
| [`get_bond_from`](#pylimer_tools_cpp.DataFileReader.get_bond_from)(self)               | Get starting atom IDs for all bonds.                   |
| [`get_bond_to`](#pylimer_tools_cpp.DataFileReader.get_bond_to)(self)                   | Get ending atom IDs for all bonds.                     |
| [`get_bond_types`](#pylimer_tools_cpp.DataFileReader.get_bond_types)(self)             | Get all bond types from the data file.                 |
| [`get_high_x`](#pylimer_tools_cpp.DataFileReader.get_high_x)(self)                     | Get the upper bound of the box in x direction.         |
| [`get_high_y`](#pylimer_tools_cpp.DataFileReader.get_high_y)(self)                     | Get the upper bound of the box in y direction.         |
| [`get_high_z`](#pylimer_tools_cpp.DataFileReader.get_high_z)(self)                     | Get the upper bound of the box in z direction.         |
| [`get_low_x`](#pylimer_tools_cpp.DataFileReader.get_low_x)(self)                       | Get the lower bound of the box in x direction.         |
| [`get_low_y`](#pylimer_tools_cpp.DataFileReader.get_low_y)(self)                       | Get the lower bound of the box in y direction.         |
| [`get_low_z`](#pylimer_tools_cpp.DataFileReader.get_low_z)(self)                       | Get the lower bound of the box in z direction.         |
| [`get_lx`](#pylimer_tools_cpp.DataFileReader.get_lx)(self)                             | Get the box length in x direction.                     |
| [`get_ly`](#pylimer_tools_cpp.DataFileReader.get_ly)(self)                             | Get the box length in y direction.                     |
| [`get_lz`](#pylimer_tools_cpp.DataFileReader.get_lz)(self)                             | Get the box length in z direction.                     |
| [`get_masses`](#pylimer_tools_cpp.DataFileReader.get_masses)(self)                     | Get the mass values for each atom type.                |
| [`get_molecule_ids`](#pylimer_tools_cpp.DataFileReader.get_molecule_ids)(self)         | Get all molecule IDs from the data file.               |
| [`get_nr_of_atom_types`](#pylimer_tools_cpp.DataFileReader.get_nr_of_atom_types)(self) | Get the number of atom types in the data file.         |
| [`get_nr_of_atoms`](#pylimer_tools_cpp.DataFileReader.get_nr_of_atoms)(self)           | Get the number of atoms in the data file.              |
| [`get_nr_of_bond_types`](#pylimer_tools_cpp.DataFileReader.get_nr_of_bond_types)(self) | Get the number of bond types in the data file.         |
| [`get_nr_of_bonds`](#pylimer_tools_cpp.DataFileReader.get_nr_of_bonds)(self)           | Get the number of bonds in the data file.              |
| [`read`](#pylimer_tools_cpp.DataFileReader.read)(self, path_of_file_to_read[, ...])    | Actually read a LAMMPS's write_data file.              |

### Methods Documentation

#### get_atom_ids(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get all atom IDs from the data file.

* **Returns:**
  Vector of atom IDs

#### get_atom_nx(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get periodic image flags in x direction for all atoms.

* **Returns:**
  Vector of nx values (image flags)

#### get_atom_ny(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get periodic image flags in y direction for all atoms.

* **Returns:**
  Vector of ny values (image flags)

#### get_atom_nz(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get periodic image flags in z direction for all atoms.

* **Returns:**
  Vector of nz values (image flags)

#### get_atom_types(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get all atom types from the data file.

* **Returns:**
  Vector of atom types

#### get_atom_x(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Get x coordinates of all atoms.

* **Returns:**
  Vector of x coordinates

#### get_atom_y(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Get y coordinates of all atoms.

* **Returns:**
  Vector of y coordinates

#### get_atom_z(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]

Get z coordinates of all atoms.

* **Returns:**
  Vector of z coordinates

#### get_bond_from(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get starting atom IDs for all bonds.

* **Returns:**
  Vector of starting atom IDs for bonds

#### get_bond_to(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get ending atom IDs for all bonds.

* **Returns:**
  Vector of ending atom IDs for bonds

#### get_bond_types(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get all bond types from the data file.

* **Returns:**
  Vector of bond types

#### get_high_x(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the upper bound of the box in x direction.

* **Returns:**
  Upper x boundary

#### get_high_y(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the upper bound of the box in y direction.

* **Returns:**
  Upper y boundary

#### get_high_z(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the upper bound of the box in z direction.

* **Returns:**
  Upper z boundary

#### get_low_x(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the lower bound of the box in x direction.

* **Returns:**
  Lower x boundary

#### get_low_y(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the lower bound of the box in y direction.

* **Returns:**
  Lower y boundary

#### get_low_z(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the lower bound of the box in z direction.

* **Returns:**
  Lower z boundary

#### get_lx(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the box length in x direction.

* **Returns:**
  Box length in x direction

#### get_ly(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the box length in y direction.

* **Returns:**
  Box length in y direction

#### get_lz(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [float](https://docs.python.org/3/library/functions.html#float)

Get the box length in z direction.

* **Returns:**
  Box length in z direction

#### get_masses(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[int](https://docs.python.org/3/library/functions.html#int), [float](https://docs.python.org/3/library/functions.html#float)]

Get the mass values for each atom type.

* **Returns:**
  Map of atom types to their masses

#### get_molecule_ids(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Get all molecule IDs from the data file.

* **Returns:**
  Vector of molecule IDs

#### get_nr_of_atom_types(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of atom types in the data file.

* **Returns:**
  Number of atom types

#### get_nr_of_atoms(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of atoms in the data file.

* **Returns:**
  Number of atoms

#### get_nr_of_bond_types(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of bond types in the data file.

* **Returns:**
  Number of bond types

#### get_nr_of_bonds(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader)) → [int](https://docs.python.org/3/library/functions.html#int)

Get the number of bonds in the data file.

* **Returns:**
  Number of bonds

#### read(self: [pylimer_tools_cpp.DataFileReader](#pylimer_tools_cpp.DataFileReader), path_of_file_to_read: [str](https://docs.python.org/3/library/stdtypes.html#str), atom_style: [pylimer_tools_cpp.AtomStyle](pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle) = AtomStyle.ANGLE, atom_style_2: [pylimer_tools_cpp.AtomStyle](pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle) = AtomStyle.NONE, atom_style_3: [pylimer_tools_cpp.AtomStyle](pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle) = AtomStyle.NONE) → [None](https://docs.python.org/3/library/constants.html#None)

Actually read a LAMMPS’s write_data file.

* **Parameters:**
  * **path_of_file_to_read** – The path to the file to read
  * **atom_style** – The format of the “Atoms” section, see [https://docs.lammps.org/read_data.html](https://docs.lammps.org/read_data.html)
  * **atom_style2** – The format of the “Atoms” section if the previous parameter is equal to AtomStyle::HYBRID
  * **atom_style_3** – The format of the “Atoms” section if the second to last parameter is equal to AtomStyle::HYBRID
