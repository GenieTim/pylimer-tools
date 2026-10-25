# Getting Started Guide

This guide will help you get started with pylimer-tools, covering basic usage patterns and common workflows.

## Basic Usage

After [installing](installation.md) pylimer-tools, you can import the modules in Python:

```python
import pylimer_tools
import pylimer_tools_cpp
```

The main functionality is split between two modules:

- **pylimer_tools**: Pure Python utilities and high-level analysis functions
- **pylimer_tools_cpp**: High-performance C++ backend with Python bindings

## Core Concepts

Understanding the key data structures is essential for working with pylimer-tools:

### Universe

The [`Universe`](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) class represents a collection of atoms and their interactions at a specific time point.

```python
from pylimer_tools_cpp import Universe

# Create a universe (typically loaded from file)
universe = Universe()

# Basic properties
print(f"Number of atoms: {universe.get_size()}")
print(f"Simulation box volume: {universe.get_volume()}")
```

### UniverseSequence

The [`UniverseSequence`](api/pylimer_tools_cpp.UniverseSequence.md#pylimer_tools_cpp.UniverseSequence) class manages multiple universes (e.g., trajectory data) efficiently:

```python
from pylimer_tools_cpp import UniverseSequence

# Load trajectory data
sequence = UniverseSequence()
sequence.initialize_from_data_sequence(["frame1.structure.dat", "frame2.structure.dat"])

# Access individual frames
first_frame = sequence.at_index(0)
last_frame = sequence.at_index(sequence.get_length() - 1)
```

## Common Workflows

### 1. Loading and Analyzing LAMMPS Data

```python
import numpy as np
from pylimer_tools_cpp import UniverseSequence

# Load LAMMPS data file
sequence = UniverseSequence()
sequence.initialize_from_data_sequence(["your_system.data"])
universe = sequence.at_index(0)

# Basic system information
print(f"System contains {universe.get_size()} atoms")
print(f"Simulation box: {universe.get_box().get_x_length():.2f} x " +
      f"{universe.get_box().get_y_length():.2f} x " +
      f"{universe.get_box().get_z_length():.2f}")

# Analyze molecules
molecules = universe.get_molecules(atom_type_to_omit=2)  # Omit type 2 atoms
print(f"Number of molecules: {len(molecules)}")

# Calculate properties for each molecule
for i, molecule in enumerate(molecules[:5]):  # First 5 molecules
    bond_lengths = molecule.compute_bond_lengths()
    print(f"Molecule {i}: {len(bond_lengths)} bonds, "
          f"mean length = {bond_lengths.mean():.3f}")
```

### 2. Polymer Analysis

```python
# Calculate polymer-specific properties
radii_of_gyration = []
end_to_end_distances = []

for molecule in universe.get_molecules():
    if molecule.get_size() > 2:  # Skip single atoms/dimers
        rg = molecule.compute_radius_of_gyration()
        ree = molecule.compute_end_to_end_distance()

        radii_of_gyration.append(rg)
        end_to_end_distances.append(ree)

# Statistics
print(f"Mean radius of gyration: {np.mean(radii_of_gyration):.3f}")
print(f"Mean end-to-end distance: {np.mean(end_to_end_distances):.3f}")
```

### 3. Loading Trajectory Data

```python
# Load multiple frames from dump files
sequence = UniverseSequence()
sequence.initialize_from_dump_file(
   initial_data_file="initial_complete_frame.data",
   dump_file="trajectory.lammpstrj")

print(f"Trajectory contains {len(sequence)} frames")

# Analyze time evolution
times = []
total_energies = []

for universe in sequence:
    # Process each frame
    pass
```

### 5. Writing Modified Systems

```python
from pylimer_tools_cpp import DataFileWriter
from pylimer_tools.io.read_lammps_output_file import read_data_file

# Load system
universe = read_data_file("input.data")

# Create writer and configure output
writer = DataFileWriter(universe)
writer.config_include_angles(True)
writer.config_include_velocities(False)
writer.config_reindex_atoms(True)

# Write to file
writer.write_to_file("modified_system.data")
```

## Performance Tips

1. **Use C++ functions when possible**: The `pylimer_tools_cpp` module is optimized for performance
2. **Batch operations**: Process multiple properties in a single loop over molecules
3. **Memory efficiency**: For large systems, avoid storing all results in memory simultaneously
4. **Vectorized operations**: Use NumPy operations on the results from pylimer-tools functions

## Next Steps

- Dive into the [Assumptions](assumptions.md) and [Nomenclature](nomenclature.md) for understanding the underlying concepts and terminology used in pylimer-tools
- Explore the [File I/O: Readers & Writers](auto_examples/readers_writers/index.md) guide for detailed file I/O operations
- Check out [Network Generation](auto_examples/network_generator/index.md) for creating custom polymer networks
- See [Examples](auto_examples/index.md) for running simulations and computations with pylimer-tools
- Review the [src](modules.md) for complete API reference
- Look at the [test suite](https://github.com/GenieTim/pylimer-tools/tree/main/tests) for more examples

For help and support, visit the [GitHub repository](https://github.com/GenieTim/pylimer-tools) or check the [Assumptions](assumptions.md) and [Nomenclature](nomenclature.md) pages for important implementation details.
