"""
Normal Mode Analysis (NMA)
==========================

This example demonstrates how to use the NormalModeAnalyzer in pylimer-tools
to predict the loss and storage modulus of polymer networks.

Reference: :cite:t:`gusev_molecular_2024`
"""

import os

import numpy as np
import pandas as pd

from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import NormalModeAnalyzer

# Load a LAMMPS data file (replace with your file)
filePath = os.path.join(
    os.getcwd(),
    "../..",
    "tests/pylimer_tools/fixtures/structure/network_100_a_46.structure.out",
)
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
# TODO: apply conversion factors if necessary
# plot the results
df = pd.DataFrame(
    {
        "Frequency": frequencies,
        "Storage Modulus": storage_modulus,
        "Loss Modulus": loss_modulus,
    }
)
df.set_index("Frequency", inplace=True)
_ = df.plot(title="Storage and Loss Modulus vs Frequency", loglog=True)
