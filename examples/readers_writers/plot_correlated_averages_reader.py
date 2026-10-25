#!/usr/bin/env python
"""
Correlated Averages Reader
==========================

Read correlation functions from ``fix ave/correlate`` and ``fix ave/correlate/long`` commands.
"""

import os
from pathlib import Path

from pylimer_tools.io.read_lammps_output_file import read_correlation_file

_fixtures_dir = (
    Path(__file__).resolve().parent.parent.parent /
    "tests/pylimer_tools/fixtures"
    if "__file__" in globals()
    else (
        Path.cwd().parent / "tests/pylimer_tools/fixtures"
        if Path.cwd().name == "docs"
        else Path("tests/pylimer_tools/fixtures").resolve()
    )
)
correlation_file = str(_fixtures_dir / "example_correlation_file.out.corr.txt")

# Read correlation data
correlation_data = read_correlation_file(correlation_file)

# Access correlation data by timestep
last_timestep = correlation_data["Timestep"].max()
last_correlation = correlation_data[correlation_data["Timestep"]
                                    == last_timestep]

_ = last_correlation.plot(
    x="Time", y=[
        "c_pA1*c_pA1", "c_pA2*c_pA2", "c_pA3*c_pA3"])
