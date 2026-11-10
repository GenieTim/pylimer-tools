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
from pathlib import Path

from pylimer_tools_cpp import DumpFileReader, UniverseSequence

# Load trajectory data
sequence = UniverseSequence()


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
sequence.initialize_from_dump_file(
    initial_data_file=os.path.join(file_path, "lammps_data_file_small.out"),
    dump_file=os.path.join(file_path, "lammps_dump_small.lammpstrj"),
)

print(f"Trajectory contains {len(sequence)} frames")

# Access specific frames
first_frame = sequence.at_index(0)
last_frame = sequence.at_index(len(sequence) - 1)  # Last frame

dump_file_reader = DumpFileReader(
    os.path.join(file_path, "lammps_dump_small_3step.lammpstrj")
)
print(f"Dump file contains {dump_file_reader.get_length()} frames")

# %%
# Visualizing Box Dimensions Across Trajectory Frames
# ---------------------------------------------------
import matplotlib.pyplot as plt

frame_indices = list(range(len(sequence)))
box_volumes = [sequence.at_index(i).get_box().get_volume()
               for i in frame_indices]

plt.figure(figsize=(6, 4))
plt.plot(frame_indices, box_volumes, marker="o", color="#4C72B0", lw=2)
plt.xlabel("Frame Index", fontsize=11)
plt.ylabel("Box Volume [$\\sigma^3$]", fontsize=11)
plt.title("Simulation Box Volume Across Trajectory Frames", fontsize=12)
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plt.show()
