#!/usr/bin/env python
"""
Save Generated Network
======================

This example shows how to save a generated polymer network to a LAMMPS data file.
"""

from pathlib import Path

from pylimer_tools_cpp import AtomStyle, DataFileWriter, MCUniverseGenerator

# Generate network
generator = MCUniverseGenerator(45, 45, 45)
# ... configure and generate ...
universe = generator.get_universe()

# Create output directory
output_dir = Path("generated_networks")
output_dir.mkdir(exist_ok=True)

# Write LAMMPS data file
writer = DataFileWriter(universe)
writer.config_include_angles(True)
writer.config_reindex_atoms(True)
writer.config_atom_style(AtomStyle.ANGLE)

output_file = output_dir / "polymer_network.data"
writer.write_to_file(str(output_file))

print(f"Network saved to: {output_file}")
