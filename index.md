# pylimer-tools Documentation

[![Test Coverage of Python Code](https://github.com/GenieTim/pylimer-tools/blob/main/.github/coverage.svg?raw=true)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)[![C++ Code Coverage](https://github.com/GenieTim/pylimer-tools/blob/main/.github/cpp-coverage.svg?raw=true)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)[![Test Coverage](https://codecov.io/gh/GenieTim/pylimer-tools/branch/main/graph/badge.svg?token=5ZE1VSDXJQ)](https://codecov.io/gh/GenieTim/pylimer-tools)[![Tests](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml/badge.svg)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)[![PyPI version](https://badge.fury.io/py/pylimer-tools.svg)](https://badge.fury.io/py/pylimer-tools)[![Downloads](https://img.shields.io/pypi/dm/pylimer-tools.svg)](https://pypi.python.org/pypi/pylimer-tools/)

Welcome to pylimer-tools, a comprehensive Python library for working with bead-spring polymer systems and for analyzing LAMMPS molecular dynamics simulation output.
This library provides powerful tools for reading, processing, and analyzing polymer networks with a focus on performance and ease of use.

## 🚀 **Quick Start**

Install pylimer-tools with pip:

```bash
pip install pylimer-tools
```

Basic usage example:

```python
import numpy as np
from pylimer_tools_cpp import UniverseSequence

# Load a LAMMPS data file
filePath = "your_lammps_data_file.structure.out"
universeSequence = UniverseSequence()
universeSequence.initialize_from_data_sequence([filePath])
universe = universeSequence.at_index(0)

# Analyze the system
print(f"System size: {universe.get_size()}")
print(f"Volume: {universe.get_volume()} u³")
print(f"Mean bond length: {np.mean(universe.compute_bond_lengths())} u")
```

## 📚 **Documentation Sections**

## 🔧 **Key Features**

- **Polymer Analysis**: Calculate radius of gyration, end-to-end distances, molecular weights, bond and loop statistics
- **Network Generation**: Monte Carlo network generators for polymer systems with highly flexible parameters
- **Force Balance**: Reduce polymer networks to their minimum energy, maximum entropy homogenized state, predict the equilibrium shear modulus
- **Normal Mode Analysis**: Predict the loss and storage modulus of polymer networks
- **Dissipative Particle Dynamics**: Simulate coarse-grained polymer systems with slip-springs and DPD
- **High Performance**: C++ backend with Python bindings for optimal performance on large polymer structures
- **LAMMPS Integration**: Seamlessly read and write LAMMPS data and dump files with memory-efficient streaming
- **Comprehensive**: Support for angles, dihedrals, crosslinkers, velocities, charges, and custom polymer properties
- **Trajectory Analysis**: Advanced tools for analyzing molecular dynamics trajectories
- **Batch Processing**: Efficient tools for processing multiple files and time series data

## 💡 **Quick Navigation**

- **New to pylimer-tools?** Start with [Installation Guide](installation.md) and [Getting Started Guide](usage.md)
- **Working with files?** Check out [File I/O: Readers & Writers](auto_examples/readers_writers/index.md)
- **Generating networks?** See [Network Generation](auto_examples/network_generator/index.md)
- **Need examples?** Browse [Examples](auto_examples/index.md) for practical code samples
- **Understanding limitations?** Review [Assumptions](assumptions.md) and [Nomenclature](nomenclature.md)
- **Like papers?** See Bernhard *et al.* [[BSG25](acknowledgements.md#id25)] and [Acknowledgements](acknowledgements.md) for related publications

## 🔗 **Useful Links**

- [GitHub Repository](https://github.com/GenieTim/pylimer-tools)
- [PyPI Package](https://pypi.org/project/pylimer-tools/)
- [Issue Tracker](https://github.com/GenieTim/pylimer-tools/issues)
- [Tests](https://github.com/GenieTim/pylimer-tools/tree/main/tests)
- [Examples](https://github.com/GenieTim/pylimer-tools/tree/main/examples)

## 📖 **Indices and Tables**

* [Index](genindex.md)
* [Module Index](py-modindex.md)
* [Search Page](search.md)
* [References](acknowledgements.md#references)
