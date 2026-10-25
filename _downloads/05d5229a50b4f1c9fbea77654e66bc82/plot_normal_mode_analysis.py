"""
Normal Mode Analysis (NMA)
==========================

This example demonstrates how to use the NormalModeAnalyzer in pylimer-tools
to predict the loss and storage modulus of polymer networks.

Reference: :cite:t:`gusev_molecular_2024`
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import NormalModeAnalyzer

# Load a LAMMPS data file (replace with your file)
def _find_fixtures_dir() -> Path:
    for cur in ([Path(__file__).resolve().parent] if "__file__" in globals() else []) + [Path.cwd()]:
        for p in [cur] + list(cur.parents):
            candidate = p / "tests" / "pylimer_tools" / "fixtures"
            if candidate.is_dir():
                return candidate
    return Path("tests/pylimer_tools/fixtures").resolve()


_fixtures_dir = _find_fixtures_dir()
filePath = str(_fixtures_dir / "structure/network_100_a_46.structure.out")
universe = read_data_file(filePath)
edges = universe.get_edges()

# Create the normal mode analyzer
nma = NormalModeAnalyzer(
    edges["edge_from"],
    edges["edge_to"],
)

# Perform the analysis
nma.find_all_eigenvalues(compute_eigenvectors=False)

# Get the storage and loss modulus at a given frequency
frequencies = np.logspace(-4, 3, num=int(1e3))  # Example frequency range
storage_modulus = nma.evaluate_storage_modulus(omega=frequencies)
loss_modulus = nma.evaluate_loss_modulus(omega=frequencies)

# %%
# .. todo:
#   Apply conversion factors where necessary
#
# Plot the results
df = pd.DataFrame(
    {
        "Frequency": frequencies,
        "Storage Modulus": storage_modulus,
        "Loss Modulus": loss_modulus,
    }
)
df.set_index("Frequency", inplace=True)
_ = df.plot(title="Storage and Loss Modulus vs Frequency", loglog=True)
