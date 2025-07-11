#!/usr/bin/env python
"""
Convert for Bond Swapping
=========================

LAMMPS bond swapping can significantly speed up equilibration for polymer networks.
However, to preserve the strand lengths, you need to ensure that the molecule indices are set correctly.
Here's how you could re-set the molecule indices for bond swapping with constant chain lengths.
"""

import os

from pylimer_tools_cpp import DataFileWriter, UniverseSequence

if not os.path.exists("generated_structures"):
    os.makedirs("generated_structures")

# Load original system
input_file = os.path.join(
    os.getcwd(),
    "../..",
    "tests/pylimer_tools/fixtures/structure/melt_83_a_100.structure.out",
)
sequence = UniverseSequence()
sequence.initialize_from_data_sequence([input_file])
universe = sequence.at_index(0)

# Configure for optimal chain length distribution
output_file = os.path.join(
    os.getcwd(), "generated_structures", "for_swap_" + os.path.basename(input_file)
)

writer = DataFileWriter(universe)
writer.config_include_angles(True)
writer.config_crosslinker_type(2)  # Needed for decomposition to chains
writer.config_molecule_idx_for_swap(True)  # Key setting

writer.write_to_file(output_file)
print(f"Optimized system capable of bond swapping written to: {output_file}")

# %%
# To convert it back to the molecule index relating to their current connectivity,
# simply repeat the process with ``config_molecule_idx_for_swap(False)``
