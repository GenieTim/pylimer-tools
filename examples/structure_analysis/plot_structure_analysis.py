#!/usr/bin/env python
"""
Extract Synthesis Parameters
============================

This example demonstrates how to analyze the structure of a crosslinked polymer network
to get the stoichiometric imbalance :math:`r`,
the crosslinker conversion :math:`p` or
the functionality of the crosslinker :math:`f`.
"""

import os
from pathlib import Path
from collections import Counter

import matplotlib.pyplot as plt

from pylimer_tools.analyse_networks import (
    compute_crosslinker_conversion,
    compute_stoichiometric_imbalance,
)
from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import Universe

# Replace with your crosslinker type
crosslinker_type = 2

# Load your network (replace path accordingly)
_fixtures_dir = (
    Path(__file__).resolve().parent.parent.parent /
    "tests/pylimer_tools/fixtures"
    if "__file__" in globals()
    else (
        Path.cwd().parent / "tests/pylimer_tools/fixtures"
        if Path.cwd().name == "docs"
        else Path("tests/pylimer_tools/fixtures").resolve()
    )
)
universe = read_data_file(
    str(_fixtures_dir / "structure/network_100_a_46.structure.out")
)
assert isinstance(universe, Universe)

# 1. Compute stoichiometric imbalance
r = compute_stoichiometric_imbalance(universe, crosslinker_type)
print(f"Stoichiometric imbalance (r): {r}")

# 2. Compute crosslinker conversion
p = compute_crosslinker_conversion(universe, crosslinker_type)
print(f"Crosslinker conversion (p): {p}")

# 3. Compute functionality of the crosslinker
f = universe.determine_functionality_per_type()[crosslinker_type]
print(f"Functionality of the crosslinker (f): {f}")

# Count functionalities
crosslinkers = universe.get_atoms_by_type(crosslinker_type)
functionalities = universe.get_vertex_degrees()
functionalities = [
    functionalities[universe.get_vertex_idx_by_atom_id(a.get_id())]
    for a in crosslinkers
]
functionality_counts = Counter(functionalities)

# Plot the functionality counts
plt.figure()
plt.bar(
    [f"{k} ({v})" for k, v in functionality_counts.items()],
    list(functionality_counts.values()),
)
plt.xlabel("Functionality")
plt.ylabel("Count")
plt.show()
