#!/usr/bin/env python
"""
Data File Reader
================

The :class:`~pylimer_tools_cpp.DataFileParser` handles LAMMPS data files containing complete system definitions.
The :class:`~pylimer_tools_cpp.UniverseSequence` provides a convenient way to have these data files read
into the :class:`~pylimer_tools_cpp.Universe` objects.
"""

import os
from pathlib import Path

from pylimer_tools_cpp import UniverseSequence


def _find_fixtures_dir() -> Path:
    for cur in (
        [Path(__file__).resolve().parent] if "__file__" in globals() else []
    ) + [Path.cwd()]:
        for p in [cur] + list(cur.parents):
            candidate = p / "tests" / "pylimer_tools" / "fixtures"
            if candidate.is_dir():
                return candidate
    return Path("tests/pylimer_tools/fixtures").resolve()


_fixtures_dir = _find_fixtures_dir()
base_structure_path = str(_fixtures_dir / "structure")

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
            "3d-diamond-lattice_10x10x10_a_3_d_0.85_v_0.V-fixed.structure.out",
        ),
    ]
)
for universe in sequence:
    print(f"Loaded universe with {universe.get_nr_of_atoms()} atoms")

# %%
# Please note that the ``hybrid`` and ``tdpd`` atom styles are not fully supported.
#
# If your data files don't hint the atom style in the line that starts with "Atom",
# you should specify it manually, see :meth:`~pylimer_tools_cpp.DataFileReader.read`
# and :func:`~pylimer_tools_cpp.UniverseSequence.set_data_file_atom_style()`.
