#!/usr/bin/env python
"""
pylimer-tools Output Reader
===========================

Read output from pylimer-tools' own simulators.

Relevant documentation: :func:`~pylimer_tools.io.read_pylimer_tools_output_file.read_avg_file`

"""

import os

import matplotlib.pyplot as plt

from pylimer_tools.io.read_pylimer_tools_output_file import read_avg_file

example_file = os.path.join(
    os.getcwd(),
    "..",
    "dpd_simulations/dpd_simulation_avg_output.txt",
)

# Read pylimer-tools averages file
pylimer_data = read_avg_file(example_file)

# Data is automatically grouped by OutputStep
print(f"Simulation steps: {pylimer_data['Step'].nunique()}")

# Access specific measurements
plt.plot(pylimer_data["Step"], pylimer_data["Temperature"])
plt.xlabel("Output Step")
plt.ylabel("Temperature")
plt.show()
