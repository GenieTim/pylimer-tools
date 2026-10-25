#!/usr/bin/env python
"""
Performance and Caching
=======================

All readers support automatic caching for improved performance.
"""

import os
from pathlib import Path
import time

import matplotlib.pyplot as plt

from pylimer_tools.io.extract_thermo_data import extract_thermo_params


# Read complete log file
def _find_fixtures_dir() -> Path:
    for cur in (
        [Path(__file__).resolve().parent] if "__file__" in globals() else []
    ) + [Path.cwd()]:
        for p in [cur] + list(cur.parents):
            candidate = p / "tests" / "pylimer_tools" / "fixtures"
            if candidate.is_dir():
                return candidate
    return Path("tests/pylimer_tools/fixtures").resolve()


_fixtures_dir = _find_fixtures_dir()
log_file = str(_fixtures_dir / "log.lammps")

time_start = time.time()
# %%
# :func:`~pylimer_tools.io.extract_thermo_data.extract_thermo_params` is a more elaborate version of
# :func:`~pylimer_tools.io.read_lammps_output_file.read_log_file`
thermo_data = extract_thermo_params(log_file, use_cache=False)
time_end = time.time()

time_no_cache = time_end - time_start

time_start = time.time()
thermo_data = extract_thermo_params(log_file, use_cache=True)
time_end = time.time()

time_cache = time_end - time_start

plt.figure()
plt.bar(["Without Cache", "With Cache"], [time_no_cache, time_cache])
plt.ylabel("Time [s]")
plt.show()
