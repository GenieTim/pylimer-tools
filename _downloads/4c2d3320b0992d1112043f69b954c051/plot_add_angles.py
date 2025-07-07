#!/usr/bin/env python
"""
Adding Topology Information
===========================

The generated networks can include topological information
"""

from pylimer_tools_cpp import MCUniverseGenerator

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
    f"Added {universe.get_nr_of_angles()} angles and {universe.get_nr_of_dihedral_angles()} dihedral angles."
)
