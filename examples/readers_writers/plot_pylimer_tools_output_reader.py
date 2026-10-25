#!/usr/bin/env python
"""
pylimer-tools Output Reader
===========================

Read output from pylimer-tools' own simulators.

Relevant documentation: :func:`~pylimer_tools.io.read_pylimer_tools_output_file.read_avg_file`

"""

import os
from pathlib import Path

import matplotlib.pyplot as plt

from pylimer_tools.io.read_pylimer_tools_output_file import read_avg_file

_examples_dir = (
    Path(__file__).resolve().parent.parent
    if "__file__" in globals()
    else (
        Path.cwd().parent / "examples"
        if Path.cwd().name == "docs"
        else Path("examples").resolve()
    )
)
example_file = str(
    _examples_dir /
    "dpd_simulations/dpd_simulation_avg_output.txt")

# Read pylimer-tools averages file
pylimer_data = read_avg_file(example_file)

# Data is automatically grouped by OutputStep
print(f"Simulation steps: {pylimer_data['Step'].nunique()}")

# Access specific measurements
plt.plot(pylimer_data["Step"], pylimer_data["Temperature"])
plt.xlabel("Output Step")
plt.ylabel("Temperature")
plt.show()
