#!/usr/bin/env python
"""
Correlated Averages Reader
==========================

Read correlation functions from ``fix ave/correlate`` and ``fix ave/correlate/long`` commands:
"""
import os

from pylimer_tools.io.read_lammps_output_file import read_correlation_file

correlation_file = os.path.join(
    os.getcwd(),
    "../..",
    "tests/pylimer_tools/fixtures/example_correlation_file.out.corr.txt",
)

# Read correlation data
correlation_data = read_correlation_file(correlation_file)

# Access correlation data by timestep
last_timestep = correlation_data["Timestep"].max()
last_correlation = correlation_data[correlation_data["Timestep"] == last_timestep]

last_correlation.plot(x="Time", y=["c_pA1*c_pA1", "c_pA2*c_pA2", "c_pA3*c_pA3"])
