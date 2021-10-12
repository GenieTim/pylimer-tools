# Pylimer-Tools

[![Run Tests](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml/badge.svg)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)
[![Publish Documentation](https://github.com/GenieTim/pylimer-tools/actions/workflows/publish-documentation.yml/badge.svg)](https://github.com/GenieTim/pylimer-tools/actions/workflows/publish-documentation.yml)
[![PyPI version](https://badge.fury.io/py/pylimer-tools.svg)](https://badge.fury.io/py/pylimer-tools)
[![PyPI download month](https://img.shields.io/pypi/dm/pylimer-tools.svg)](https://pypi.python.org/pypi/pylimer-tools/)
[![PyPI license](https://img.shields.io/pypi/l/pylimer-tools.svg)](https://pypi.python.org/pypi/pylimer-tools/)

A collection of utility python functions for handling LAMMPS output and polymers in Python.

This toolbox provides means to read LAMMPS output: be it data, dump or thermodynamic data files. 
Additionally, it provides various methods to calculate with the read data, such as computing the 
radius of gyration, mean end to end distance, or simply splitting a polymer network back up into its chains.

## Installation

Use pip:

`pip install pylimer-tools`

## Usage

**NOTE**: currently, this release's API is _unstable_ and subject to change.

See the [documentation](https://genietim.github.io/pylimer-tools/index.html) (or in [markdown format](https://github.com/GenieTim/pylimer-tools/tree/main/docs/pylimer_tools)) for a current list of all available functions.

### Example

An example useage can be found in the [CLI application](./src/pylimer_tools/pylimer_tools.py) or in the following code snippet:

```python
import numpy as np

from pylimer_tools.entities.universum import Universum
from pylimer_tools.io.readLammpData import readLammpData

allData = readLammpData("some_lammps_output_file.dat")
universe = Universum(allData['atom_data'], allData['bond_data'], boxSizes=[
                      allData["Lx"], allData["Ly"], allData["Lz"]])
print("Size: {}. Volume: {} u^3".format(
    universe.getSize(), universe.getVolume()))
print("Mean bond length: {} u".format(
    np.mean([m.computeBondLengths().mean() for m in universe])))
print("Mean end to end distance: {} u".format(
    np.mean([m.computeEndToEndDistance() for m in universe])))
```
