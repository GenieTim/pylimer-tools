File I/O: Readers & Writers
===========================

pylimer-tools provides comprehensive support for reading and writing various LAMMPS file formats, enabling seamless integration with molecular dynamics workflows.

Overview
--------

The library supports four main file I/O operations:

1. **Reading LAMMPS data files** - Complete system configurations
2. **Reading LAMMPS dump files** - Trajectory and snapshot data  
3. **Reading LAMMPS output files** - Measurements and simulation results
4. **Writing LAMMPS data files** - Modified or generated systems

Readers
-------

Data File Reader
~~~~~~~~~~~~~~~~

The :class:`~pylimer_tools_cpp.DataFileParser` handles LAMMPS data files containing complete system definitions.
The :class:`~pylimer_tools_cpp.UniverseSequence` provides a convenient way to have these data files read into the :class:`~pylimer_tools_cpp.Universe` objects.

.. code-block:: python

   from pylimer_tools_cpp import UniverseSequence
   
   # Load a single data file
   one_universe_seq = UniverseSequence()
   one_universe_seq.initialize_from_data_sequence(["system.structure.dat"])
   universe = one_universe_seq.at_index(0)
   
   # Load multiple data files as a sequence
   sequence = UniverseSequence()
   sequence.initialize_from_data_sequence([
       "config1.structure.dat",
       "config2.structure.dat", 
       "config3.structure.dat"
   ])
   for universe in sequence:
       print(f"Loaded universe with {universe.get_size()} atoms")

**Supported atom styles:**

- ``angle``
- ``bond``
- ``molecular``
- ``charge```
- ``full``
- ``hybrid`` of ``bond`` and ``edpd``

If your data files don't hint the atom style, you can specify it manually, 
see :func:`~pylimer_tools_cpp.DataFileParser.read()` and :func:`~pylimer_tools_cpp.UniverseSequence.set_data_file_atom_style()`.

Dump File Reader  
~~~~~~~~~~~~~~~~

The :class:`~pylimer_tools_cpp.DumpFileReader` processes LAMMPS trajectory files:

.. code-block:: python

   # Load trajectory data
   sequence = UniverseSequence()
   sequence.initialize_from_dump_sequence(initial_data_file="network.structure.dat", dump_file="trajectory.lammpstrj")
   
   print(f"Trajectory contains {sequence.get_length()} frames")
   
   # Access specific frames
   first_frame = sequence.at_index(0)
   last_frame = sequence.at_index(-1)  # Last frame

**Supported dump formats:**

- Standard LAMMPS dump format
- Custom dump formats with user-defined columns
- Compressed files (automatically detected)

Memory-Efficient Reading
~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~pylimer_tools_cpp.UniverseSequence` class provides automatic memory management:

.. code-block:: python

   # For large trajectories, frames are loaded on-demand
   sequence = UniverseSequence()
   sequence.initialize_from_dump_sequence(initial_data_file="network.structure.dat", dump_file="large_trajectory.lammpstrj")
   
   # Only the requested frame is loaded into memory
   for i in range(0, sequence.get_length(), 10):  # Every 10th frame
       universe = sequence.at_index(i)
       # Process universe
       # …
       # Then, optional: free memory for the frame 
       sequence.forget_at_index(i)

Output File Reader
~~~~~~~~~~~~~~~~~~

The Python module provides comprehensive support for reading various LAMMPS output file formats, 
including thermodynamic data, time-averaged quantities, histograms, and correlation functions.

Log & Thermo File Reader
^^^^^^^^^^^^^^^^^^^^^^^^

Read LAMMPS log files containing thermodynamic output from the ``thermo`` command:

.. code-block:: python

   from pylimer_tools.io.read_lammps_output_file import read_log_file
   
   # Read complete log file
   thermo_data = read_log_file("simulation.log")
   
   # Access thermodynamic properties
   print(f"Temperature range: {thermo_data['Temp'].min():.2f} - {thermo_data['Temp'].max():.2f}")
   print(f"Final energy: {thermo_data['TotEng'].iloc[-1]:.6f}")
   
   # Plot energy evolution
   import matplotlib.pyplot as plt
   plt.plot(thermo_data['Step'], thermo_data['TotEng'])
   plt.xlabel('Simulation Step')
   plt.ylabel('Total Energy')
   plt.show()

The log file reader automatically:

- Detects different thermodynamic output sections with varying columns
- Handles (ignores) warnings and broken lines in the log
- Skips non-numeric data
- Supports multi-run simulations with different ``thermo_style`` commands

Relevant documentation: :func:`~pylimer_tools.io.read_lammps_output_file.read_log_file` and :func:`~pylimer_tools.io.extract_thermo_data.extract_thermo_params`

Time-Averaged Data Reader
^^^^^^^^^^^^^^^^^^^^^^^^^

Read output from LAMMPS ``fix ave/time`` commands:

.. code-block:: python

   from pylimer_tools.io.read_lammps_output_file import read_averages_file
   
   # Read simple averages file
   averages = read_averages_file("pressure_averages.txt")
   
   # Read sectioned averages (e.g., from vector quantities from e.g. `fix ave/time ... vector`)
   sectioned_data = read_averages_file("vector_averages.txt")
   
   # Group by timestep for sectioned data
   for timestep, group in sectioned_data.groupby('Timestep'):
       print(f"Timestep {timestep}: {len(group)} data points")

Relevant documentation: :func:`~pylimer_tools.io.read_lammps_output_file.read_averages_file`
and :func:`~pylimer_tools.io.read_lammps_output_file.read_sectioned_averages_file`

**Supported averaging formats:**

- Scalar quantities (single values per output timestep)
- Vector quantities (multiple values per output timestep)
- Multi-section files with varying column structures

Histogram Reader
^^^^^^^^^^^^^^^^

Process histogram data from ``fix ave/hist`` commands:

.. code-block:: python

   from pylimer_tools.io.read_lammps_output_file import read_histogram_file
   
   # Read histogram data
   histogram_data = read_histogram_file("density_histogram.txt")
   
   # Plot histogram for specific timestep
   final_timestep = histogram_data['Timestep'].max()
   final_hist = histogram_data[histogram_data['Timestep'] == final_timestep]
   
   plt.bar(final_hist['Coord1'], final_hist['Count'])
   plt.xlabel('Density')
   plt.ylabel('Frequency')
   plt.title(f'Density Distribution at Step {final_timestep}')
   plt.show()

Relevant documentation: :func:`~pylimer_tools.io.read_lammps_output_file.read_histogram_file`

Correlated Averages Reader
^^^^^^^^^^^^^^^^^^^^^^^^^^

Read correlation functions from ``fix ave/correlate`` and ``fix ave/correlate/long`` commands:

.. code-block:: python

   from pylimer_tools.io.read_lammps_output_file import read_correlation_file
   
   # Read correlation data
   correlation_data = read_correlation_file("velocity_correlation.txt")
   
   # Access correlation data by timestep
   for timestep, group in correlation_data.groupby('Timestep'):
       plt.plot(group['SampleCorr'], group['correlation_value'], 
                label=f'Step {timestep}')
   
   plt.xlabel('Correlation Time')
   plt.ylabel('Correlation Value')
   plt.legend()
   plt.show()

Relevant documentation: :func:`~pylimer_tools.io.read_lammps_output_file.read_correlation_file`

pylimer-tools Output Reader
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read output from pylimer-tools' own simulators:

.. code-block:: python

   from pylimer_tools.io.read_pylimer_tools_output_file import read_avg_file
   
   # Read pylimer-tools averages file
   pylimer_data = read_avg_file("simulation_averages.txt")
   
   # Data is automatically grouped by OutputStep
   print(f"Simulation steps: {pylimer_data['OutputStep'].nunique()}")
   
   # Access specific measurements
   if 'Temperature' in pylimer_data.columns:
       plt.plot(pylimer_data['OutputStep'], pylimer_data['Temperature'])
       plt.xlabel('Output Step')
       plt.ylabel('Temperature')
       plt.show()

Relevant documentation: :func:`~pylimer_tools.io.read_pylimer_tools_output_file.read_avg_file`

Unit Conversion
^^^^^^^^^^^^^^^

Convert `LAMMPS units <https://docs.lammps.org/units.html>`_ to SI or other unit systems:

.. code-block:: python

   from pylimer_tools.io.unit_styles import UnitStyleFactory
   
   # Create unit style factory
   unit_factory = UnitStyleFactory()
   
   # Get unit conversion for specific LAMMPS unit style
   real_units = unit_factory.get_unit_style("real")
   metal_units = unit_factory.get_unit_style("metal")
   
   # For LJ units, specify the polymer
   lj_units = unit_factory.get_unit_style("lj", polymer="pdms")
   
   # Convert pressure from LAMMPS units to SI
   pressure_si = thermo_data['Press'] * real_units.pressure.to("Pa").magnitude
   
   # Convert temperature (already in Kelvin for 'real' units)
   temperature_si = thermo_data['Temp'] * real_units.temperature.to("K").magnitude
   
   print(f"Pressure in Pascal: {pressure_si}")

**Supported unit styles:**

- ``real`` - kcal/mol, Angstrom, fs
- ``metal`` - eV, Angstrom, ps  
- ``si`` - kg, m, s
- ``nano`` - attogram, nanometer, nanosecond
- ``lj`` - Lennard-Jones reduced units (requires polymer specification)

Performance and Caching
^^^^^^^^^^^^^^^^^^^^^^^

All readers support automatic caching for improved performance:

.. code-block:: python

   # Caching is enabled by default
   data1 = read_log_file("large_simulation.log")  # Reads from file
   data2 = read_log_file("large_simulation.log")  # Reads from cache (faster)
   
   # Disable caching if needed
   data3 = read_averages_file("averages.txt", use_cache=False)
   
   # Cache is automatically invalidated when files change
   # No manual cache management needed

Error Handling
^^^^^^^^^^^^^^

The readers handle various error conditions gracefully:

.. code-block:: python

   try:
       # Read potentially problematic file
       data = read_log_file("simulation.log")
   except FileNotFoundError:
       print("Log file not found")
   except RuntimeError as e:
       print(f"Error reading log file: {e}")
   
   # Check for empty results
   if data.empty:
       print("No data found in file")
   else:
       print(f"Successfully read {len(data)} data points")

**Common issues handled:**

- Missing or corrupted files
- Inconsistent column formats
- Mixed data types
- Incomplete simulations
- Unicode encoding issues
- Large file memory management

Writers
-------

Data File Writer
~~~~~~~~~~~~~~~~

The :class:`~pylimer_tools_cpp.DataFileWriter` generates LAMMPS-compatible data files:

Basic Usage
^^^^^^^^^^^

.. code-block:: python

   from pylimer_tools_cpp import DataFileWriter
   
   # Create writer for a universe
   writer = DataFileWriter(universe)
   
   # Basic output
   writer.write_to_file("output.structure.dat")

Advanced Configuration
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Configure output options
   writer = DataFileWriter(universe)
   
   # Include topology information  
   writer.config_include_angles(True)
   writer.config_include_dihedral_angles(True)
   
   # Control atom numbering
   writer.config_reindex_atoms(True)
   writer.config_crosslinker_type(2)
   
   # Handle periodic boundaries
   writer.config_move_into_box(True)
   writer.config_attempt_image_reset(True)
   
   # Include velocities
   writer.config_include_velocities(True)
   
   # Write configured output
   writer.write_to_file("configured_output.data")

Custom Atom Styles
^^^^^^^^^^^^^^^^^^

You can specify custom atom styles and formats:

.. code-block:: python

   from pylimer_tools_cpp import AtomStyle
   
   # Set predefined atom style
   writer.config_atom_style(AtomStyle.FULL)
   
   # Or define custom format
   custom_format = "$atomId $moleculeId $atomType $charge $x $y $z"
   writer.set_custom_atom_format(custom_format)

Note that currently, only the atom styles ``angle``, ``bond``, ``molecular``, and ``full`` are supported without you specifying a custom format.

Molecule Indexing Control
^^^^^^^^^^^^^^^^^^^^^^^^^

Control how molecules are numbered in the output:

.. code-block:: python

   # Standard molecular indexing: sets the molecule index for each atom 
   # based on the decomposition at the junctions
   writer.config_molecule_idx_for_swap(False)
   
   # Optimize for LAMMPS' bond/swap operations
   # Sets the molecule index for each atom based on the distance from the next junction
   writer.config_molecule_idx_for_swap(True)

Examples and Use Cases
----------------------

Converting Between Formats
~~~~~~~~~~~~~~~~~~~~~~~~~~

Take a later entry from a LAMMPS trajectory and save it as a data file:

.. code-block:: python

   # Load from dump file
   sequence = UniverseSequence()
   sequence.initialize_from_dump_file(
    initial_data_file="initial_complete_frame.data",
    dump_file="trajectory.lammpstrj"
   )
   
   # Get the last frame, for example
   final_universe = sequence.at_index(len(sequence) - 1)
   
   # Write as data file
   writer = DataFileWriter(final_universe)
   writer.config_include_angles(True)
   writer.write_to_file("final_configuration.data")

Preparing Systems for LAMMPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prepare a generated network for LAMMPS simulation:

.. code-block:: python

   from pylimer_tools_cpp import MCUniverseGenerator, DataFileWriter
   
   # Generate network (see network_generator documentation)
   generator = MCUniverseGenerator(50, 50, 50)
   # ... configure and generate network ...
   universe = generator.get_universe()
   
   # Set up for LAMMPS
   universe.set_masses({1: 1.0, 2: 2.0})  # Set masses per atom type
   
   # Detect and add angles
   angles = universe.detect_angles()
   universe.add_angles(
       angles["angle_from"],
       angles["angle_via"], 
       angles["angle_to"],
       angle_types=[1 for _ in range(len(angles["angle_from"]))]
   )
   
   # Write LAMMPS-ready file
   writer = DataFileWriter(universe)
   writer.config_include_angles(True)
   writer.config_atom_style(AtomStyle.ANGLE)
   writer.write_to_file("simulation_ready.data")

Faster Equilibration With Bond Swapping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LAMMPS bond swapping can significantly speed up equilibration for polymer networks.
However, to preserve the strand lengths, you need to ensure that the molecule indices are set correctly.
Here's how you could re-set the molecule indices for bond swapping with constant chain lengths:

.. code-block:: python

   import os
   from pylimer_tools_cpp import UniverseSequence, DataFileWriter
   
   # Load original system
   input_file = "original_system.data"
   sequence = UniverseSequence()
   sequence.initialize_from_data_sequence([input_file])
   universe = sequence.at_index(0)
   
   # Configure for optimal chain length distribution
   output_file = "for_swap_" + os.path.basename(input_file)
   
   writer = DataFileWriter(universe)
   writer.config_include_angles(True)
   writer.config_crosslinker_type(2) # Needed for decomposition to chains
   writer.config_molecule_idx_for_swap(True)  # Key setting
   
   writer.write_to_file(output_file)
   print(f"Optimized system written to: {output_file}")

Batch Processing
~~~~~~~~~~~~~~~~

Process multiple files in a workflow:

.. code-block:: python

   import glob
   from pathlib import Path
   
   # Find all data files
   input_files = glob.glob("input_*.structure.dat")
   
   for input_file in input_files:
       # Load system
       sequence = UniverseSequence()
       # Initialize from a single data file
       # Alternatively, you could initialize with the whole `input_files` list
       # and then iterate over the sequence, instead of iterating over the files
       sequence.initialize_from_data_sequence([input_file])
       universe = sequence.at_index(0)
       
       # Process system (example: add angles)
       if not universe.has_angles():
           angles = universe.detect_angles()
           universe.add_angles(
               angles["angle_from"],
               angles["angle_via"],
               angles["angle_to"],
               angle_types=[1 for _ in range(len(angles["angle_from"]))]
           )
       
       # Write processed system
       output_file = Path(input_file).stem + "_processed.structure.dat"
       writer = DataFileWriter(universe)
       writer.config_include_angles(True)
       writer.write_to_file(output_file)
       
       print(f"Processed: {input_file} -> {output_file}")

Performance Considerations
--------------------------

Large File Handling
~~~~~~~~~~~~~~~~~~~

For very large files:

.. code-block:: python

   # Use lazy loading for large trajectories
   sequence = UniverseSequence()
   sequence.initialize_from_dump_file(
     initial_data_file="initial_complete_frame.data",
     dump_file="very_large_trajectory.lammpstrj"
   )
   
   # Process in chunks to manage memory
   chunk_size = 100
   for start_idx in range(0, sequence.get_length(), chunk_size):
       end_idx = min(start_idx + chunk_size, sequence.get_length())
       
       for i in range(start_idx, end_idx):
           universe = sequence.at_index(i)
           # you can free the memory in the sequence if having this universe is enough
           sequence.forget_at_index(i)

           # Process frame
           pass

Write Optimization
~~~~~~~~~~~~~~~~~~

For better write performance:

.. code-block:: python

   # Pre-configure the writer once
   writer = DataFileWriter(universe)
   writer.config_include_angles(True)
   writer.config_reindex_atoms(True)
   
   # Write multiple related files efficiently
   for i, modified_universe in enumerate(universe_variations):
       writer.set_universe_to_write(modified_universe)
       writer.write_to_file(f"variation_{i}.data")

Error Handling
--------------

Don't forget to handle potential file I/O errors:

.. code-block:: python

   try:
       sequence = UniverseSequence()
       sequence.initialize_from_data_sequence(["system.data"])
   except FileNotFoundError:
       print("Error: Input file not found")
   except RuntimeError as e:
       print(f"Error reading file: {e}")
   
   try:
       writer = DataFileWriter(universe)
       writer.write_to_file("output.data")
   except PermissionError:
       print("Error: Cannot write to output file (permissions)")
   except OSError as e:
       print(f"Error writing file: {e}")
