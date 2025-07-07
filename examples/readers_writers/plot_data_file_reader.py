#!/usr/bin/env python
"""
Data File Reader
================

The :class:`~pylimer_tools_cpp.DataFileParser` handles LAMMPS data files containing complete system definitions.
The :class:`~pylimer_tools_cpp.UniverseSequence` provides a convenient way to have these data files read
into the :class:`~pylimer_tools_cpp.Universe` objects.
"""

import os

from pylimer_tools_cpp import UniverseSequence

base_structure_path = os.path.join(
    os.getcwd(),
    "../..",
    "tests/pylimer_tools/fixtures/structure",
)

# Load a single data file
data_file = os.path.join(
    base_structure_path,
    "melt_83_a_100.structure.out",
)
one_universe_seq = UniverseSequence()
one_universe_seq.initialize_from_data_sequence([data_file])
universe = one_universe_seq.at_index(0)

# Load multiple data files as a sequence
sequence = UniverseSequence()
sequence.initialize_from_data_sequence(
    [
        os.path.join(
            base_structure_path,
            "3d-diamond-lattice_10x10x10_a_3_d_0.85_imperfect.structure.out",
        ),
        os.path.join(
            base_structure_path,
            "3d-diamond-lattice_10x10x10_a_3_d_0.85_v_0.V-fixed.structure.outt",
        ),
    ]
)
for universe in sequence:
    print(f"Loaded universe with {universe.get_nr_of_atoms()} atoms")

"""
**Supported atom styles:**

- ``angle``
- ``bond``
- ``molecular``
- ``charge```
- ``full``
- ``hybrid`` of ``bond`` and ``edpd``

If your data files don't hint the atom style, you can specify it manually, 
see :func:`~pylimer_tools_cpp.DataFileParser.read()` 
and :func:`~pylimer_tools_cpp.UniverseSequence.set_data_file_atom_style()`.
"""
