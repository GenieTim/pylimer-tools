#!/usr/bin/env python
"""
Performance and Caching
=======================

All readers support automatic caching for improved performance
"""

import os
import time

import matplotlib.pyplot as plt

from pylimer_tools.io.extract_thermo_data import extract_thermo_params

# Read complete log file
log_file = os.path.join(
    os.getcwd(),
    "../..",
    "tests/pylimer_tools/fixtures/log.lammps",
)

time_start = time.time()
# extract_thermo_params is a more elaborate version of
# pylimer_tools.io.read_lammps_output_file.read_log_file()
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
