#!/usr/bin/env python
"""
Correlated Averages Reader
==========================

Read correlation functions from ``fix ave/correlate`` and ``fix ave/correlate/long`` commands.
"""

import os
from pathlib import Path

from pylimer_tools.io.read_lammps_output_file import read_correlation_file

def _find_fixtures_dir() -> Path:
    for cur in ([Path(__file__).resolve().parent] if "__file__" in globals() else []) + [Path.cwd()]:
        for p in [cur] + list(cur.parents):
            candidate = p / "tests" / "pylimer_tools" / "fixtures"
            if candidate.is_dir():
                return candidate
    return Path("tests/pylimer_tools/fixtures").resolve()


_fixtures_dir = _find_fixtures_dir()
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
