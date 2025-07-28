#!/usr/bin/env python
"""
Multi-Component System
======================

Create complex systems with multiple components:
"""

import datetime

from pylimer_tools_cpp import DataFileWriter, MCUniverseGenerator

# System parameters
box_size = 60.0
seed = 42
bead_distance = 1.0

# Network parameters
n_crosslinkers = 50

n_polymer_chains = 80
beads_per_chain = 25

n_solvent_chains = 200
beads_per_solvent = 5

n_monofunctional_chains = 50
beads_per_monofunctional_chain = 10

crosslinker_conversion = 0.85
crosslinker_functionality = 4

print(f"Starting network generation at {datetime.datetime.now()}")

# Initialize generator
generator = MCUniverseGenerator(box_size, box_size, box_size)
generator.set_seed(seed)
generator.set_bead_distance(bead_distance)

# Add components in order
print("Adding crosslinkers...")
generator.add_crosslinkers(
    n_crosslinkers,
    crosslinker_functionality=crosslinker_functionality,
    crosslinker_type=2,
)  # Type 2

print("Adding solvent chains...")
generator.add_solvent_chains(n_solvent_chains, beads_per_solvent, 3)  # Type 3

print("Adding and linking polymer strands...")
# Add bifunctional polymer chains
generator.add_strands(
    nr_of_strands=n_polymer_chains,
    strand_lengths=[beads_per_chain for _ in range(n_polymer_chains)],
    strand_atom_type=1,  # Type 1 for "normal" beads
)

generator.add_monofunctional_strands(
    nr_of_monofunctional_strands=n_monofunctional_chains,
    monofunctional_strand_length=[
        beads_per_monofunctional_chain for _ in range(n_monofunctional_chains)
    ],
    # Type 4 to recognize these as "monofunctional" strands later
    monofunctional_strand_atom_type=4,
)

# Crosslink polymer strands
generator.link_strands_to_conversion(
    crosslinker_conversion=crosslinker_conversion,  # 92.5% of crosslinker sites used
)

# Get final system
universe = generator.get_universe()
universe.set_masses({1: 1.0, 2: 2.0, 3: 0.5, 4: 1.0})  # Different masses

print(f"Network generation completed at {datetime.datetime.now()}")

# Save the universe to a file
writer = DataFileWriter(universe)
writer.write_to_file("generated_networks/multi_component_polymer_network.data")
