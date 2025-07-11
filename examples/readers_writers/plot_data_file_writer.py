#!/usr/bin/env python
"""
Data File Writer
================

The :class:`~pylimer_tools_cpp.DataFileWriter` generates LAMMPS-compatible data files:
"""

import os
import random

from pylimer_tools_cpp import AtomStyle, DataFileWriter, Universe

# create folder generated_structures/
if not os.path.exists("generated_structures"):
    os.makedirs("generated_structures")

# Generate a universe to write
universe = Universe(10.0, 10.0, 10.0)
universe.add_atoms(
    list(range(100)),
    [1 for _ in range(100)],
    [random.random() * 10.0 for _ in range(100)],
    [random.random() * 10.0 for _ in range(100)],
    [random.random() * 10.0 for _ in range(100)],
    [0 for _ in range(100)],
    [0 for _ in range(100)],
    [0 for _ in range(100)],
)
# universe.add_bonds(
#   …
# )

# Create writer for a universe
writer = DataFileWriter(universe)

# Basic output
writer.write_to_file("generated_structures/output.structure.dat")

"""
Advanced Configuration
----------------------

"""

# Configure output options
writer = DataFileWriter(universe)

# Include topology information
writer.config_include_angles(True)
writer.config_include_dihedral_angles(True)

# Control atom numbering
writer.config_reindex_atoms(True)
writer.config_crosslinker_type(2)

# Handle periodic boundaries
writer.config_move_into_box(True)
writer.config_attempt_image_reset(True)

# Include velocities
writer.config_include_velocities(True)

# Write configured output
writer.write_to_file("generated_structures/configured_output.data")

"""
Custom Atom Styles
------------------

You can specify custom atom styles and formats using the following format:

"""
# Set predefined atom style
writer.config_atom_style(AtomStyle.FULL)

# Or define custom format
custom_format = "$atomId $moleculeId $atomType $charge $x $y $z"
writer.set_custom_atom_format(custom_format)

"""
Note that currently, the atom styles `tdpd` and `hybrid` are not supported without you specifying this custom format.
"""

"""
Molecule Indexing Control
-------------------------

Control how molecules are numbered in the output
"""
# Standard molecular indexing: sets the molecule index for each atom
# based on the decomposition at the junctions
writer.config_molecule_idx_for_swap(False)

# Optimize for LAMMPS' bond/swap operations
# Sets the molecule index for each atom based on the distance from the next junction
writer.config_molecule_idx_for_swap(True)
