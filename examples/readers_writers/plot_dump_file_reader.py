#!/usr/bin/env python
"""
Dump File Reader
================

The :class:`~pylimer_tools_cpp.DumpFileReader` processes LAMMPS trajectory files:
"""
from pylimer_tools_cpp import UniverseSequence

# Load trajectory data
sequence = UniverseSequence()
sequence.initialize_from_dump_file(
    initial_data_file="network.structure.dat", dump_file="trajectory.lammpstrj"
)

print(f"Trajectory contains {sequence.get_length()} frames")

# Access specific frames
first_frame = sequence.at_index(0)
last_frame = sequence.at_index(-1)  # Last frame
