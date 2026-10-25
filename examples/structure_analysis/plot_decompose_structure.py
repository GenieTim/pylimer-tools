#!/usr/bin/env python
"""
Decompose/Split Structure
=========================

This example demonstrates how to decompose a crosslinked polymer network into chains
and analyze its structure using pylimer-tools.

The following approaches are shown:

- Removing crosslinkers and analyzing polymer strands
- Decomposing the network to chains while keeping crosslinkers
- Finding clusters after manual junction removal
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt

from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import Universe

# Replace with your crosslinker type
crosslinker_type = 2


# Load your network (replace path accordingly)
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
universe = read_data_file(
    str(_fixtures_dir / "structure/network_100_a_46.structure.out")
)
assert isinstance(universe, Universe)

# 1. Decompose by removing crosslinkers
chains = universe.get_molecules(crosslinker_type)
print(f"Number of chains (without crosslinkers): {len(chains)}")

# 2. Decompose keeping crosslinkers in all chains
chains_with_xl = universe.get_chains_with_crosslinker(crosslinker_type)
print(f"Number of chains (with crosslinkers): {len(chains_with_xl)}")

# 3. Find clusters (if you removed junctions yourself)
clusters = universe.get_clusters()
print(f"Number of clusters: {len(clusters)}")

# Plot the weights of the clusters, chains and molecules
plt.figure(figsize=(10, 6))
plt.bar(
    ["Clusters", "Chains (with XL)", "Chains (no XL)"],
    [len(clusters), len(chains_with_xl), len(chains)],
)
plt.xlabel("Type")
plt.ylabel("Count")
plt.show()
