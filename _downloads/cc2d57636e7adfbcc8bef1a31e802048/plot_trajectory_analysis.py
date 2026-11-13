"""
Trajectory Analysis
===================

This example demonstrates how to analyze a polymer trajectory using pylimer-tools,
including mean square displacement (MSD) and radius of gyration calculation.
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pylimer_tools_cpp import UniverseSequence


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

# Compute the mean square displacement for all atoms
msd = sequence.compute_msd_for_atoms(
    atom_ids=[a.get_id() for a in sequence[0].get_atoms()],
    nr_of_origins=50,  # Number of origins to use for MSD calculation
)

# %%
# You can also compute other properties like radius of gyration, end-to-end distance, etc.
# by using the "lazy" iterator over the sequence
# This allows you to process large trajectories without loading everything
# into memory at once.

# Compute radius of gyration over time
rgs = []
timesteps = []
for universe in sequence:
    rgs.append(
        np.mean(
            [
                m.compute_radius_of_gyration()
                for m in universe.get_molecules(
                    -1
                )  # -1 because we don't have any crosslinks to filter out
            ]
        )
    )
    timesteps.append(universe.get_timestep())

# Plot 1: Mean Square Displacement (primary single-panel figure)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(
    list(msd.keys()),
    list(msd.values()),
    label="Mean Square Displacement",
)
ax.axhline(
    np.mean(list(msd.values())),
    color="black",
    linestyle="--",
    label="Mean",
)
ax.plot(
    list(msd.keys()),
    [np.mean(list(msd.values())[: i + 1]) for i in range(len(msd))],
    label="Cumulative Mean",
    color="red",
)
ax.set_xlabel(r"$\tau$")
ax.set_ylabel("MSD")
ax.set_title("Mean Square Displacement")
ax.set(xscale="log", yscale="log")
ax.legend()
plt.tight_layout()
plt.show()

# Plot 2: Radius of Gyration
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(timesteps, rgs, label="Radius of Gyration", color="orange")
ax.axhline(np.mean(rgs), color="black", linestyle="--", label="Mean")
ax.plot(
    timesteps,
    [np.mean(rgs[: i + 1]) for i in range(len(rgs))],
    label="Cumulative Mean",
    color="red",
)
ax.set_xlabel("Timestep")
ax.set_ylabel("Radius of Gyration")
ax.set_title("Radius of Gyration Over Time")
ax.set(xscale="log", yscale="log")
ax.legend()
plt.tight_layout()
plt.show()

# %%
# Please beware that the MSD is computed in an :math:`N^2` manner, which can be slow for large systems,
# many origins and long trajectories.
#
# If you want more origins or otherwise better performance,
# consider using the :py:mod:`tidynamics` package instead like so:
#
# .. code-block:: python
#
#    import tidynamics
#
#    # using sorted atom ids to make sure the order is consistent
#    # across all timesteps, which is important for MSD calculation
#    sorted_atom_ids = sorted([a.get_id() for a in sequence[0].get_atoms()])
#    msd_result = tidynamics.msd(
#          np.array(
#              [
#                  [
#                      coord
#                      for id in sorted_atom_ids
#                      for coord in u.get_atom_by_vertex_id(id).get_coordinates()
#                  ]
#                  for u in sequence
#              ]
#        )
#     )
#
