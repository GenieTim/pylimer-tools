#!/usr/bin/env python
"""
Generate Vulcanized Networks
============================

In this example, we create a vulcanized network using the `pylimer-tools` library.
"""

import matplotlib.pyplot as plt

from pylimer_tools_cpp import MCUniverseGenerator, randomly_sample_entanglements

# Create generator for a 50x50x50 simulation box
generator = MCUniverseGenerator(50.0, 50.0, 50.0)

# Set random seed for reproducibility
generator.set_seed(12345)

# Set mean bead-to-bead distance
generator.set_bead_distance(1.0)

# Add strands
generator.add_strands(
    nr_of_strands=200,  # 200 polymer chains
    strand_lengths=[100 for _ in range(200)],  # 100 beads per chain
)

universe = generator.get_universe()

# Sample vulcanization "crosslinks"
sampled_crosslinks = randomly_sample_entanglements(
    universe, nr_of_samples=200, upper_cutoff=2.5
)
n_sampled_crosslinks = len(sampled_crosslinks.pairs_of_atoms)

# Add sampled crosslinks to the universe
universe.add_bonds(
    n_sampled_crosslinks,
    [sampled_crosslinks.pairs_of_atoms[i][0]
        for i in range(n_sampled_crosslinks)],
    [sampled_crosslinks.pairs_of_atoms[i][1]
        for i in range(n_sampled_crosslinks)],
    [2 for _ in range(n_sampled_crosslinks)],  # Bond type 2
)

print("Generated vulcanized network with {} crosslinks".format(n_sampled_crosslinks))

# if you want, you can collapse the crosslinks to be a single atom
# as such:
universe.contract_vertices_along_bond_type(2)

# %%
# Refer to :doc:`auto_examples/readers_writers/index` for how you can save the universe to a file.
#
# For now, we will show that this produces a variety of different strand
# lengths:

strand_lengths = [m.get_nr_of_atoms()
                  for m in universe.get_chains_with_crosslinker(2)]

# Plot the distribution of strand lengths
plt.figure()
plt.hist(strand_lengths, bins=50, alpha=0.5, label="Strand Lengths")
plt.xlabel("Strand Length")
plt.ylabel("Frequency")
plt.title("Distribution of Strand Lengths in Vulcanized Network")
plt.show()
