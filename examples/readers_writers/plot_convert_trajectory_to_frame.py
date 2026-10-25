#!/usr/bin/env python
"""
Extract Trajectory Frame to Data File
=====================================

Take a later entry from a LAMMPS trajectory and save it as a data file.
"""

import os
from pathlib import Path

from pylimer_tools_cpp import DataFileWriter, UniverseSequence

if not os.path.exists("generated_structures"):
    os.makedirs("generated_structures")

# %%
# For large trajectories, frames are loaded on-demand
# Load a LAMMPS dump trajectory file (replace with your file)
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
file_path = str(_fixtures_dir)
sequence = UniverseSequence()
sequence.initialize_from_dump_file(
    initial_data_file=os.path.join(file_path, "big_dump_file_data.out"),
    dump_file=os.path.join(file_path, "big_dump_file.lammpstrj"),
)
# Get the last frame, for example
final_universe = sequence.at_index(len(sequence) - 1)

# Write as data file
writer = DataFileWriter(final_universe)
writer.config_include_angles(True)
writer.write_to_file("generated_structures/final_configuration.data")
