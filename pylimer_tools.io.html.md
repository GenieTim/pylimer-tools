# pylimer_tools.io package

## Submodules

## pylimer_tools.io.extract_thermo_data module

### pylimer_tools.io.extract_thermo_data.detect_headers(file: [str](https://docs.python.org/3/library/stdtypes.html#str), max_nr_of_lines_to_read: [int](https://docs.python.org/3/library/functions.html#int) = 1500, use_cache: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)]

Read max_nr_of_lines_to_read lines from the given file and return all possible header lines.

Some assumptions are made regarding the columns, e.g., that 75% of them start with a character.

* **Parameters:**
  * **file** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The file to search for header lines
  * **max_nr_of_lines_to_read** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The number of lines to read in search for header lines.
    Use a negative number to read the whole file.
  * **use_cache** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Whether to read the result from cache or not.
    The cache is not read if the file changed meanwhile.
* **Returns:**
  List of detected header lines
* **Return type:**
  List[[str](https://docs.python.org/3/library/stdtypes.html#str)]

### pylimer_tools.io.extract_thermo_data.extract_thermo_params(file, header: [str](https://docs.python.org/3/library/stdtypes.html#str) | [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = 'Step Temp E_pair E_mol TotEng Press', texts_to_read: [int](https://docs.python.org/3/library/functions.html#int) = 50, min_line_len: [int](https://docs.python.org/3/library/functions.html#int) = 5, use_cache: [bool](https://docs.python.org/3/library/functions.html#bool) = True, lines_to_read_to_detect_header: [int](https://docs.python.org/3/library/functions.html#int) = 100000, lines_to_read_till_header: [float](https://docs.python.org/3/library/functions.html#float) = -1) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Extract the thermodynamic outputs produced for this simulation,
i.e., in LAMMPS, by the thermo command.

In particular, this function can handle log files,
handle sections with different columns,
and handles skipping over warnings as well as broken lines.

Note: The header parameter can be an array — make sure to pay attention
when reading a file with different header sections in them.

* **Parameters:**
  * **file** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The file path to the file to read from
  * **header** (*Union* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *,* *List* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *]* *,* *None* *]*) – The header of the CSV (where to start reading at).
    Can be a string, a list of strings, or None if you want to try the detection.
  * **texts_to_read** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The number of times to expect the header
  * **min_line_len** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The minimal length of a line to be accepted as data
  * **use_cache** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Whether to use cache or not (though it will be written anyway).
    The cache is not read if the file changed meanwhile.
  * **lines_to_read_to_detect_header** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The number of lines to read when trying to detect headers
  * **lines_to_read_till_header** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The number of lines that are acceptable to skip
    until a header should have been found.
    This is useful for (a) finding the header, and
    (b) exit early if you are unsure about the header(s)
* **Returns:**
  The thermodynamic parameters
* **Return type:**
  pd.DataFrame

### pylimer_tools.io.extract_thermo_data.get_thermo_cache_name_suffix(header: [str](https://docs.python.org/3/library/stdtypes.html#str) | [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = 'Step Temp E_pair E_mol TotEng Press', texts_to_read: [float](https://docs.python.org/3/library/functions.html#float) = 50, min_line_len: [float](https://docs.python.org/3/library/functions.html#float) = 5) → [str](https://docs.python.org/3/library/stdtypes.html#str)

Compose a cache file suffix in such a way, that it distinguishes different thermo reader parameters.

* **Parameters:**
  * **header** (*Union* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *,* *List* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *]* *,* *None* *]*) – The header of the CSV (where to start reading at)
  * **texts_to_read** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The number of times to expect the header
  * **min_line_len** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The minimal length of a line to be accepted as data
* **Returns:**
  A string to be used as cache file suffix
* **Return type:**
  [str](https://docs.python.org/3/library/stdtypes.html#str)

### pylimer_tools.io.extract_thermo_data.read_multi_section_separated_value_file(file: [str](https://docs.python.org/3/library/stdtypes.html#str), separator: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, use_cache: [bool](https://docs.python.org/3/library/functions.html#bool) = True, comment: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, skip_err: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Reads a file with multiple sections that have different headers throughout the file.

This function handles files with multiple data sections that may have different column structures.
It automatically detects the separator if not specified and combines all sections into a single DataFrame.

* **Parameters:**
  * **file** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the file to read
  * **separator** (*Union* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *,* *None* *]*) – Character used to separate values in the file (auto-detected if None)
  * **use_cache** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Whether to use cached results if available
  * **comment** (*Union* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *,* *None* *]*) – Character indicating the start of comments (e.g., “#”)
  * **skip_err** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Whether to skip errors when processing sections
* **Returns:**
  Combined DataFrame containing all data from the file
* **Return type:**
  pd.DataFrame

#### NOTE
Particularly useful for reading output files from the DPDSimulator or other
multi-section files where the structure may change between sections.

### pylimer_tools.io.extract_thermo_data.read_one_group(fp, header, min_line_len=4, additional_lines_skip=0, lines_to_read_till_header=1000.0) → [str](https://docs.python.org/3/library/stdtypes.html#str)

Read one group of csv lines from the file.

* **Parameters:**
  * **fp** (*file object*) – The file pointer to the file to read from
  * **header** ([*str*](https://docs.python.org/3/library/stdtypes.html#str) *or* [*list*](https://docs.python.org/3/library/stdtypes.html#list)) – The header of the CSV (where to start reading at)
  * **min_line_len** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The minimal length of a line to be accepted as data
  * **additional_lines_skip** ([*int*](https://docs.python.org/3/library/functions.html#int)) – Number of lines to skip after reading the header
  * **lines_to_read_till_header** ([*float*](https://docs.python.org/3/library/functions.html#float)) – Maximum number of lines to read until finding the header
* **Returns:**
  The filename of a temporary CSV file, or empty string if no data was read
* **Return type:**
  [str](https://docs.python.org/3/library/stdtypes.html#str)

## pylimer_tools.io.read_pylimer_tools_output_file module

This module provides a few functions to read output from pylimer_tools_cpp’s simulators.

### pylimer_tools.io.read_pylimer_tools_output_file.read_avg_file(filename: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Read an averages-output file from one of the simulators shipped with pylimer_tools.

This function parses the output file format used by pylimer_tools_cpp simulators,
handling multiple data sections and converting them to a pandas DataFrame.
The function also caches results to improve performance on subsequent reads.

* **Parameters:**
  **filename** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the averages file to read
* **Returns:**
  DataFrame containing the parsed averages data, grouped by OutputStep
* **Return type:**
  pd.DataFrame
* **Note:**
  The function automatically filters out lines containing “-nan” values,
  null characters, or fewer than 3 columns.
* **Note:**
  The returned DataFrame is grouped by OutputStep, keeping only the last
  entry for each step.

## pylimer_tools.io.read_lammps_output_file module

This module provides a few functions to read LAMMPS’ output files, including:

- log files (thermo output)
- dump files (focusing on the coordinates of atoms)
- data files (the LAMMPS structure)
- averaged data (from `fix ave/time...` or `fix ave/hist...`)
- correlation data (from `fix ave/correlate/...`)

### pylimer_tools.io.read_lammps_output_file.read_averages_file(filepath, use_cache: [bool](https://docs.python.org/3/library/functions.html#bool) = True, sep=' ') → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Read a file written by a fix ave/time command.

Uses pandas’ read_csv after detecting the columns.

Important assumption: The first 2 or 3 lines in the file are:
: - comment,
  - then one header indicating the columns,
  - and then either data or potentially a second header, if it is a sectioned file (e.g., from a fix ave/time … vector)

* **Parameters:**
  * **filepath** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the averages file
  * **use_cache** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Whether to use the cache to speed up reading & writing
  * **sep** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Delimiter used in the file (default is space)
* **Returns:**
  DataFrame containing the parsed average data
* **Return type:**
  pd.DataFrame
* **Raises:**
  [**FileNotFoundError**](https://docs.python.org/3/library/exceptions.html#FileNotFoundError) – If the averages file does not exist

### pylimer_tools.io.read_lammps_output_file.read_correlation_file(filepath, group_key='Timestep', use_cache: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Read a file written by a fix ave/correlate{/long} command.

* **Parameters:**
  * **filepath** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the correlation file
  * **group_key** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The key that denotes a new section
  * **use_cache** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Whether to use the cache to speed up reading & writing
* **Returns:**
  DataFrame containing the correlation data. Use the group_key with
  the DataFrame’s groupby() to restore the original sections.
* **Return type:**
  pd.DataFrame
* **Raises:**
  [**FileNotFoundError**](https://docs.python.org/3/library/exceptions.html#FileNotFoundError) – If the correlation file does not exist

### pylimer_tools.io.read_lammps_output_file.read_data_file(structure_file: [str](https://docs.python.org/3/library/stdtypes.html#str), atom_style: [List](https://docs.python.org/3/library/typing.html#typing.List)[[AtomStyle](api/pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle)] | [None](https://docs.python.org/3/library/constants.html#None) = None) → [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)

Read a file with LAMMPS’ data type of structure into a Universe.

* **Parameters:**
  * **structure_file** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the structure file
  * **atom_style** (*Union* *[**List* *[*[*AtomStyle*](api/pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle) *]* *,* *None* *]*) – The atom style(s) in the structure file (defaults to AtomStyle.Molecule if None)
* **Returns:**
  Universe object representing the molecular structure
* **Return type:**
  [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)
* **Raises:**
  [**FileNotFoundError**](https://docs.python.org/3/library/exceptions.html#FileNotFoundError) – If the structure file does not exist

### pylimer_tools.io.read_lammps_output_file.read_dump_file(data_file, dump_file, atom_style: [List](https://docs.python.org/3/library/typing.html#typing.List)[[AtomStyle](api/pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle)] | [None](https://docs.python.org/3/library/constants.html#None) = None) → [UniverseSequence](api/pylimer_tools_cpp.UniverseSequence.md#pylimer_tools_cpp.UniverseSequence)

Read a file with LAMMPS’ dump of snapshots of structures into a Universe.

* **Parameters:**
  * **data_file** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the LAMMPS data file containing structure information
  * **dump_file** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the LAMMPS dump file containing trajectory information
  * **atom_style** (*Union* *[**List* *[*[*AtomStyle*](api/pylimer_tools_cpp.AtomStyle.md#pylimer_tools_cpp.AtomStyle) *]* *,* *None* *]*) – The atom style(s) used in the data file
* **Returns:**
  Sequence of Universe objects representing the trajectory
* **Return type:**
  [UniverseSequence](api/pylimer_tools_cpp.UniverseSequence.md#pylimer_tools_cpp.UniverseSequence)

### pylimer_tools.io.read_lammps_output_file.read_histogram_file(filepath, use_cache: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Read a file written by fix ave/hist or similar.

This is a wrapper around read_sectioned_averages_file for histogram data.

* **Parameters:**
  * **filepath** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the histogram file
  * **use_cache** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Whether to use the cache to speed up reading & writing
* **Returns:**
  DataFrame containing the parsed histogram data
* **Return type:**
  pd.DataFrame
* **See:**
  [`read_sectioned_averages_file()`](#pylimer_tools.io.read_lammps_output_file.read_sectioned_averages_file)

### pylimer_tools.io.read_lammps_output_file.read_log_file(filepath, lines_to_read_to_detect_header=500000) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Read a LAMMPS’ log (thermo output) file.

* **Parameters:**
  * **filepath** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the LAMMPS log file
  * **lines_to_read_to_detect_header** ([*int*](https://docs.python.org/3/library/functions.html#int)) – Maximum number of lines to read when detecting the header
* **Returns:**
  DataFrame containing the parsed thermo data
* **Return type:**
  pd.DataFrame

### pylimer_tools.io.read_lammps_output_file.read_sectioned_averages_file(filepath, use_cache: [bool](https://docs.python.org/3/library/functions.html#bool) = True) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Read a file written by a fix ave/time command with multiple sections.

Use the section delimiter columns together with pandas’ groupby()
to restore the original sections.

* **Parameters:**
  * **filepath** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Path to the sectioned averages file
  * **use_cache** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Whether to use the cache to speed up reading & writing
* **Returns:**
  DataFrame containing the parsed sectioned data
* **Return type:**
  pd.DataFrame
* **Raises:**
  * [**FileNotFoundError**](https://docs.python.org/3/library/exceptions.html#FileNotFoundError) – If the file does not exist
  * [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If the file format is not recognized as a proper sectioned averages file

## pylimer_tools.io.unit_styles module

### *class* pylimer_tools.io.unit_styles.UnitStyle(unit_configuration: [dict](https://docs.python.org/3/library/stdtypes.html#dict), ureg: [UnitRegistry](https://pint.readthedocs.io/en/stable/api/base.html#pint.UnitRegistry))

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

UnitStyle: a collection of units of a particular LAMMPS unit style,
but in SI units
(i.e.: use this to convert your LAMMPS output data to SI units).

Example usage:

```python
unit_style_factory = UnitStyleFactory()
unit_style = unit_style_factory.get_unit_style(
    "lj", polymer="pdms", warning=False, accept_mol=True)

# multiply with the following factor to convert LJ stress to SI units,
# namely MPa in this example:
lj_stress_to_si_conversion_factor = (1.*unit_style.pressure).to("MPa").magnitude
```

Initialize a UnitStyle object.

* **Parameters:**
  * **unit_configuration** ([*dict*](https://docs.python.org/3/library/stdtypes.html#dict)) – Dictionary containing unit definitions
  * **ureg** (*UnitRegistry*) – Pint unit registry to use for unit conversions

#### \_\_getattr_\_(property: [str](https://docs.python.org/3/library/stdtypes.html#str))

Shorthand access for `get_base_unit_of()`.

* **Parameters:**
  **property** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The property name to get the unit for
* **Returns:**
  The unit object for the requested property
* **Return type:**
  [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity)

Example usage:

```python
units = get_unit_style("lj")
mass_with_units = mass_in_lj * units.mass
```

#### get_base_unit_of(property: [str](https://docs.python.org/3/library/stdtypes.html#str))

Returns the conversion factor from the unit style to SI units.

* **Parameters:**
  **property** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The property name to get the unit for (e.g., “mass”, “distance”)
* **Returns:**
  The unit object for the requested property
* **Return type:**
  [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity)

Example usage:

```python
units = get_unit_style("lj")
mass_in_si = mass_in_lj * units.get_base_unit_of("mass")
```

#### get_underlying_unit_registry()

Get the underlying Pint unit registry.

* **Returns:**
  The unit registry used by this UnitStyle
* **Return type:**
  UnitRegistry

### *class* pylimer_tools.io.unit_styles.UnitStyleFactory

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

This is a factory to get multiple instances of different
`UnitStyle`
using the same UnitRegistry, such that they are compatible.

Initialize the UnitStyleFactory with a new UnitRegistry.

#### get_available_polymers() → [list](https://docs.python.org/3/library/stdtypes.html#list)

List all available polymers for which we have lj unit conversions.

* **Returns:**
  List of polymer names
* **Return type:**
  [list](https://docs.python.org/3/library/stdtypes.html#list)

#### get_everares_et_al_data()

Load the Everaers et al. (2020) unit properties data.

* **Returns:**
  PolymerDataFrame containing polymer properties from Everaers et al.
* **Return type:**
  PolymerDataFrame

#### get_unit_registry()

Get the underlying unit registry.

* **Returns:**
  The unit registry used by this factory
* **Return type:**
  UnitRegistry

#### get_unit_style(unit_type: [str](https://docs.python.org/3/library/stdtypes.html#str), dimension: [int](https://docs.python.org/3/library/functions.html#int) = 3, \*\*kwargs) → [UnitStyle](#pylimer_tools.io.unit_styles.UnitStyle)

Get a UnitStyle instance corresponding to the unit system requested.

* **Parameters:**
  * **unit_type** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The unit type, e.g. “lj”, “nano”, “real”, “si”, etc.
  * **dimension** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The dimension of the box
  * **kwargs** ([*dict*](https://docs.python.org/3/library/stdtypes.html#dict)) – Additional arguments required for certain unit styles
* **Returns:**
  A UnitStyle object for the requested unit system
* **Return type:**
  [UnitStyle](#pylimer_tools.io.unit_styles.UnitStyle)
* **Raises:**
  * [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If required parameters are missing
  * [**NotImplementedError**](https://docs.python.org/3/library/exceptions.html#NotImplementedError) – If the requested unit type is not implemented

For LJ units, you must specify the polymer using the polymer parameter.

#### SEE ALSO
[https://docs.lammps.org/units.html](https://docs.lammps.org/units.html)

#### WARNING
Please check the source code of this function to see
whether the units you need are correctly implemented

## Module contents

This module contains various utility functions for working with
LAMMPS and pylimer_tools in- and output files.
