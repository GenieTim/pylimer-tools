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
