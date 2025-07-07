#!/usr/bin/env python
"""
Dump File Reader
================

The :class:`~pylimer_tools_cpp.DumpFileReader` processes LAMMPS trajectory files.
Is has a practical interface to read frames from a trajectory file
and access them as :class:`~pylimer_tools_cpp.Universe` objects
in the :class:`~pylimer_tools_cpp.UniverseSequence` class.
"""
import os

from pylimer_tools_cpp import DumpFileReader, UniverseSequence

# Load trajectory data
sequence = UniverseSequence()
file_path = os.path.join(
    os.getcwd(),
    "../..",
    "tests/pylimer_tools/fixtures/",
)
sequence.initialize_from_dump_file(
    initial_data_file=os.path.join(file_path, "lammps_data_file_small.out"),
    dump_file=os.path.join(file_path, "lammps_dump_small.lammpstrj"),
)

print(f"Trajectory contains {sequence.get_length()} frames")

# Access specific frames
first_frame = sequence.at_index(0)
last_frame = sequence.at_index(-1)  # Last frame

dump_file_reader = DumpFileReader(
    os.path.join(file_path, "lammps_dump_small_3step.lammpstrj")
)
print(f"Dump file contains {dump_file_reader.get_length()} frames")
