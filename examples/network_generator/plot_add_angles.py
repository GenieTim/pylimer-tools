#!/usr/bin/env python
"""
Adding Topology Information
===========================

The generated networks can include topological information.
"""

# sphinx_gallery_thumbnail_path = '_static/thumbnails/network_generator/plot_add_angles.png'
import os
from pylimer_tools_cpp import DataFileWriter, MCUniverseGenerator

# Generate basic network
generator = MCUniverseGenerator(40, 40, 40)
generator.set_seed(54321)

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

# Get the generated universe
universe = generator.get_universe()

# Detect and add angles
angles = universe.detect_angles()
if len(angles["angle_from"]) > 0:
    universe.add_angles(
        angles["angle_from"],
        angles["angle_via"],
        angles["angle_to"],
        angle_types=[1 for _ in range(len(angles["angle_from"]))],
    )
    print(f"Added {len(angles['angle_from'])} angles to the system")

# Detect dihedral angles (if needed)
dihedrals = universe.detect_dihedral_angles()
if len(dihedrals["dihedral_angle_from"]) > 0:
    # Add dihedral information if your force field requires it
    universe.add_dihedral_angles(
        dihedrals["dihedral_angle_from"],
        dihedrals["dihedral_angle_via1"],
        dihedrals["dihedral_angle_via2"],
        dihedrals["dihedral_angle_to"],
        angle_types=[1 for _ in range(len(dihedrals["dihedral_angle_from"]))],
    )

print(
    f"Added {universe.get_nr_of_angles()} angles and {
        universe.get_nr_of_dihedral_angles()
    } dihedral angles."
)

# Save the universe to a file
writer = DataFileWriter(universe)
writer.config_include_angles(True)
writer.config_include_dihedral_angles(True)

if not os.path.exists("generated_structures"):
    os.makedirs("generated_structures")

writer.write_to_file("generated_structures/polymer_network_with_angles.data")

# %%
# Visualizing Topology Element Counts
# -----------------------------------
# We plot the counts of all generated topological elements:

import matplotlib.pyplot as plt

categories = ["Atoms", "Bonds", "Angles", "Dihedrals"]
counts = [
    universe.get_nr_of_atoms(),
    universe.get_nr_of_bonds(),
    universe.get_nr_of_angles(),
    universe.get_nr_of_dihedral_angles(),
]

plt.figure(figsize=(6.5, 4))
bars = plt.bar(
    categories,
    counts,
    color=["#4C72B0", "#55A868", "#DD8452", "#C44E52"],
    edgecolor="black",
    linewidth=0.8,
    width=0.5,
)
for bar in bars:
    yval = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2.0,
        yval + 50,
        f"{yval:,}",
        ha="center",
        va="bottom",
        fontsize=10,
    )

plt.ylabel("Count", fontsize=11)
plt.title("Network Topology Elements with Detected Angles", fontsize=12)
plt.grid(True, axis="y", linestyle=":", alpha=0.6)
plt.tight_layout()
plt.show()
