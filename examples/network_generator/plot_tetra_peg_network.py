#!/usr/bin/env python
"""
Tetra-PEG Network Generation
=============================

This example demonstrates how to create a tetra-PEG (tetra-functional polyethylene glycol)
network using the new `~pylimer_tools_cpp.MCUniverseGenerator.add_star_crosslinkers`
and `~pylimer_tools_cpp.MCUniverseGenerator.link_strands_to_strands_to_conversion` methods.

Tetra-PEG networks are model systems commonly used in soft matter research where
4-functional star polymers are mixed and crosslinked to form a well-defined network structure.

The thumbnail was created using the `polymer-graph-sketcher <https://genietim.github.io/polymer-graph-sketcher/>`_ tool.
"""

# sphinx_gallery_thumbnail_path = '_static/thumbnails/network_generator/tetra-link.png'
import datetime
import os

from pylimer_tools_cpp import DataFileWriter, MCUniverseGenerator, MoleculeType

# System parameters
box_size = 80.0
seed = 12345
bead_distance = 1.0

# Network parameters - typical for tetra-PEG systems
n_star_crosslinkers = 150  # Number of 4-functional star polymers
star_functionality = 4  # Each star has 4 arms
arm_length = 20  # Beads per arm (PEG chain length)

# Additional free chains to simulate precursor chains before crosslinking
n_free_chains = 50
free_chain_length = 15

# Crosslinking parameters
strand_to_strand_conversion = 0.9  # Target conversion for strand-strand linking

print(f"Starting tetra-PEG network generation at {datetime.datetime.now()}")
print("System parameters:")
print(f"  Box size: {box_size}³")
print(f"  Star crosslinkers: {n_star_crosslinkers}")
print(f"  Star functionality: {star_functionality}")
print(f"  Arm length: {arm_length} beads")
print(f"  Free chains: {n_free_chains} of {free_chain_length} beads")
print(f"  Target strand-strand conversion: {strand_to_strand_conversion}")

# Initialize generator
generator = MCUniverseGenerator(box_size, box_size, box_size)
generator.set_seed(seed)
generator.set_bead_distance(bead_distance)

# Let's trust the Brownian Bridge and not run additional MC steps
generator.config_nr_of_mc_steps(0)

# Add star crosslinkers - these create pre-formed 4-functional stars
print("Adding star crosslinkers with pre-attached arms...")
generator.add_star_crosslinkers(
    nr_of_stars=n_star_crosslinkers,
    functionality=star_functionality,
    beads_per_strand=arm_length,
    crosslinker_atom_type=2,  # Star center atoms
    strand_atom_type=1,  # PEG chain atoms
    white_noise=True,
)

# Add some additional free chains to make the system more realistic
print("Adding additional free polymer chains...")
generator.add_strands(
    nr_of_strands=n_free_chains,
    strand_lengths=[free_chain_length for _ in range(n_free_chains)],
    strand_atom_type=1,  # Same type as PEG chains
)

print("System before strand-strand linking:")
print(f"  Total atoms: {generator.get_current_nr_of_atoms()}")
print(f"  Total bonds: {generator.get_current_nr_of_bonds()}")

# Now link free strand ends to each other to form the final network
print(
    f"Linking strand ends to achieve {
        strand_to_strand_conversion:.1%} conversion...")
generator.link_strands_to_strands_to_conversion(
    target_strand_conversion=strand_to_strand_conversion
)

print("System after strand-strand linking:")
print(f"  Total atoms: {generator.get_current_nr_of_atoms()}")
print(f"  Total bonds: {generator.get_current_nr_of_bonds()}")
print(
    f"  Final strand conversion: {
        generator.get_current_strand_conversion():.3f}")

# Optional: Relax crosslink positions for better equilibration
print("Relaxing crosslink positions...")
generator.relax_crosslinks()

# Generate the final universe with bead positions
print("Generating final universe with bead coordinates...")
universe = generator.get_universe()

# Set realistic masses (in atomic mass units)
# Type 1: PEG segments (CH2-CH2-O repeat unit ≈ 44 amu)
# Type 2: Star centers (could be larger, e.g., pentaerythritol-based ≈ 100 amu)
universe.set_masses({1: 44.0, 2: 100.0})

print(f"Network generation completed at {datetime.datetime.now()}")
print("Final network statistics:")
print(f"  Total atoms: {universe.get_nr_of_atoms()}")
print(f"  Total bonds: {universe.get_nr_of_bonds()}")
box = universe.get_box()
print(
    f"  Box dimensions: {
        box.get_lx():.1f} × {
            box.get_ly():.1f} × {
                box.get_lz():.1f}")

# Estimate network properties
total_star_arms = n_star_crosslinkers * star_functionality
total_free_chain_ends = (
    n_free_chains * 2  # Free chains have 2 ends
    + total_star_arms  # each star has one free end
)

print("Network analysis:")
print(f"  Star crosslinkers: {n_star_crosslinkers}")
print(f"  Total star arms: {total_star_arms}")
print(f"  Free chain ends (before linking): {total_free_chain_ends}")
print(
    "  Free chain ends (after linking): {}".format(
        sum(
            [
                (
                    1
                    if c.get_strand_type() == MoleculeType.DANGLING_CHAIN
                    else 2
                    if c.get_strand_type() == MoleculeType.FREE_CHAIN
                    else 0
                )
                for c in universe.get_chains_with_crosslinker(2)
            ]
        )
    )
)
print(
    "  Fully free chains (after linking): {}".format(
        sum(
            1
            for c in universe.get_chains_with_crosslinker(2)
            if c.get_strand_type() == MoleculeType.FREE_CHAIN
        )
    )
)


# Save the universe to a file
print("Saving network structure...")
writer = DataFileWriter(universe)

if not os.path.exists("generated_structures"):
    os.makedirs("generated_structures")

writer.write_to_file("generated_structures/tetra_peg_network.data")
print("Network saved to 'generated_structures/tetra_peg_network.data'")
