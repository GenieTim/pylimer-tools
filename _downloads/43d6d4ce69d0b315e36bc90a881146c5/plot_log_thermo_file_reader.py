#!/usr/bin/env python
"""
Log & Thermo File Reader
========================

Read LAMMPS log files containing thermodynamic output from the ``thermo`` command.

The log file reader automatically:

- Detects different thermodynamic output sections with varying columns
- Handles (ignores) warnings and broken lines in the log
- Skips non-numeric data
- Supports multi-run simulations with different ``thermo_style`` commands

Relevant documentation: :func:`~pylimer_tools.io.read_lammps_output_file.read_log_file`
and :func:`~pylimer_tools.io.extract_thermo_data.extract_thermo_params`
"""

import os

import matplotlib.pyplot as plt

from pylimer_tools.io.read_lammps_output_file import read_log_file

# Read complete log file
log_file = os.path.join(
    os.getcwd(),
    "../..",
    "tests/pylimer_tools/fixtures/log.lammps",
)
thermo_data = read_log_file(log_file)

# Access thermodynamic properties
print(
    "Temperature range: {:.2f} - {:.2f}".format(
        thermo_data["Temp"].min(), thermo_data["Temp"].max()
    )
)
print("Final pressure: {:.6f}".format(thermo_data["Press"].iloc[-1]))

# Plot temperature evolution
plt.plot(thermo_data["Step"], thermo_data["Temp"])
plt.xlabel("Simulation Step")
plt.ylabel("Temperature [LJ-units]")
plt.show()
plt.ylabel("Temperature [LJ-units]")
plt.show()
