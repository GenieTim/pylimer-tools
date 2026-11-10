#!/usr/bin/env python
"""
Simple Crosslinked Network
==========================

Here's a minimal example of creating a crosslinked polymer network:
"""

import os
from pylimer_tools_cpp import DataFileWriter, MCUniverseGenerator

# Create generator for a 50x50x50 simulation box
generator = MCUniverseGenerator(50.0, 50.0, 50.0)

# Set random seed for reproducibility
generator.set_seed(12345)

# Set mean bead-to-bead distance
generator.set_bead_distance(1.0)

# Add crosslinking sites
generator.add_crosslinkers(100, 4)  # 100 crosslinkers of functionality 4

# Add polymer chains
generator.add_strands(
    nr_of_strands=200,  # 200 polymer chains
    strand_lengths=[20 for _ in range(200)],  # 20 beads per chain
)

# Crosslink polymer strands
generator.link_strands_to_conversion(
    crosslinker_conversion=0.925,  # 92.5% of crosslinker sites used
)

# Get the generated network
universe = generator.get_universe()

# Set masses, e.g. for LAMMPS
universe.set_masses({1: 1.0, 2: 2.0})

print(f"Generated network with {universe.get_nr_of_atoms()} atoms")

# Save the universe to a file
writer = DataFileWriter(universe)

if not os.path.exists("generated_structures"):
    os.makedirs("generated_structures")

writer.write_to_file("generated_structures/simple_crosslinked_network.data")

# %%
# Crosslinker Functionality Distribution
# --------------------------------------
# We analyze and visualize the effective functionality of the crosslinker
# nodes after linking:

from collections import Counter
import matplotlib.pyplot as plt

crosslinkers = universe.get_atoms_by_type(2)
degrees = universe.get_vertex_degrees()
func_counts = Counter(
    degrees[universe.get_vertex_idx_by_atom_id(a.get_id())] for a in crosslinkers
)

plt.figure(figsize=(6, 4))
labels = [f"{k}-functional" for k in sorted(func_counts.keys())]
values = [func_counts[k] for k in sorted(func_counts.keys())]
bars = plt.bar(
    labels, values, color="#4C72B0", edgecolor="black", linewidth=0.8, width=0.5
)

for bar in bars:
    yval = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2.0,
        yval + 1,
        f"{yval}",
        ha="center",
        va="bottom",
        fontsize=10,
    )

plt.xlabel("Effective Crosslinker Functionality", fontsize=11)
plt.ylabel("Number of Crosslinkers", fontsize=11)
plt.title("Crosslinker Functionality Distribution", fontsize=12)
plt.grid(True, axis="y", linestyle=":", alpha=0.6)
plt.tight_layout()
plt.show()
