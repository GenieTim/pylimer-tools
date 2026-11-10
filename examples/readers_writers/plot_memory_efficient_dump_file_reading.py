#!/usr/bin/env python
"""
Memory-Efficient Dump File Reading
==================================

The :class:`~pylimer_tools_cpp.UniverseSequence` class provides automatic memory management.
"""

import os
from pathlib import Path

import psutil

from pylimer_tools_cpp import UniverseSequence


# %%
# For large trajectories, frames are loaded on-demand
# Load a LAMMPS dump trajectory file (replace with your file)
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
file_path = str(_fixtures_dir)
sequence = UniverseSequence()
sequence.initialize_from_dump_file(
    initial_data_file=os.path.join(file_path, "big_dump_file_data.out"),
    dump_file=os.path.join(file_path, "big_dump_file.lammpstrj"),
)

# Only the requested frame is loaded into memory
process = psutil.Process(os.getpid())
mem_before = process.memory_info().rss

frame_indices = []
memory_samples = []

for i in range(0, sequence.get_length(), 10):  # Every 10th frame
    universe = sequence.at_index(i)
    frame_indices.append(i)
    memory_samples.append(process.memory_info().rss / (1024 * 1024))
    print(f"Processing frame {i} with {universe.get_nr_of_atoms()} atoms")
    # Then, optional: free memory for the frame
    sequence.forget_at_index(i)

mem_after = process.memory_info().rss

dump_file_path = os.path.join(file_path, "big_dump_file.lammpstrj")
dump_file_size = os.path.getsize(dump_file_path)

print("\n--- Memory Usage Report ---")
print(f"Memory before processing: {mem_before / (1024 * 1024):.2f} MB")
print(f"Memory after processing:  {mem_after / (1024 * 1024):.2f} MB")
print(
    f"Memory used by script:    {(mem_after - mem_before) / (1024 * 1024):.2f} MB")

# %%
# Memory Stability Visualization
# ------------------------------
# Demonstrating that on-demand loading with ``forget_at_index`` maintains
# a flat memory profile:

import matplotlib.pyplot as plt

plt.figure(figsize=(6.5, 4))
plt.plot(frame_indices, memory_samples, marker="o", color="#4C72B0", lw=2)
plt.xlabel("Trajectory Frame Index", fontsize=11)
plt.ylabel("Process RSS [MB]", fontsize=11)
plt.title("Memory Footprint with On-Demand Frame Loading", fontsize=12)
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plt.show()
print(f"Dump file size:           {dump_file_size / (1024 * 1024):.2f} MB")
