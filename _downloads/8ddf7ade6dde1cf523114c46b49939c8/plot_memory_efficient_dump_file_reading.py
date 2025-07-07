#!/usr/bin/env python
"""
Memory-Efficient Dump File Reading
================================

The :class:`~pylimer_tools_cpp.UniverseSequence` class provides automatic memory management:

"""

import os

import psutil

from pylimer_tools_cpp import UniverseSequence

# For large trajectories, frames are loaded on-demand
# Load a LAMMPS dump trajectory file (replace with your file)
file_path = os.path.join(
    os.getcwd(),
    "../..",
    "tests/pylimer_tools/fixtures/",
)
sequence = UniverseSequence()
sequence.initialize_from_dump_file(
    initial_data_file=os.path.join(file_path, "big_dump_file_data.out"),
    dump_file=os.path.join(file_path, "big_dump_file.lammpstrj"),
)

# Only the requested frame is loaded into memory
process = psutil.Process(os.getpid())
mem_before = process.memory_info().rss

for i in range(0, sequence.get_length(), 10):  # Every 10th frame
    universe = sequence.at_index(i)
    # Process universe
    # …
    print(f"Processing frame {i} with {universe.get_nr_of_atoms()} atoms")
    # Then, optional: free memory for the frame
    sequence.forget_at_index(i)

mem_after = process.memory_info().rss

dump_file_path = os.path.join(file_path, "big_dump_file.lammpstrj")
dump_file_size = os.path.getsize(dump_file_path)

print("\n--- Memory Usage Report ---")
print(f"Memory before processing: {mem_before / (1024*1024):.2f} MB")
print(f"Memory after processing:  {mem_after / (1024*1024):.2f} MB")
print(f"Memory used by script:    {(mem_after - mem_before) / (1024*1024):.2f} MB")
print(f"Dump file size:           {dump_file_size / (1024*1024):.2f} MB")
