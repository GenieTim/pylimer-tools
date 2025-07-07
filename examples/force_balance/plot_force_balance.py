"""
Force Balance
=============

This example demonstrates how to use the MEHP implementation with slip-springs
to homogenize a polymer network and predict its equilibrium shear modulus.
"""

import os

import numpy as np

from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import (
    MEHPForceBalance,
    Universe,
)

# Load your network (replace with your file)
universe = read_data_file(
    os.path.join(
        os.getcwd(),
        "../..",
        "tests/pylimer_tools/fixtures/structure/network_100_a_46.structure.out",
    )
)
assert isinstance(universe, Universe)

# 2. MEHPForceBalance: Hookean springs, allows slip-links
mehp_fb = MEHPForceBalance.construct_with_random_sliplinks(
    universe, nr_of_sliplinks_to_sample=200
)
mehp_fb.run_force_relaxation(max_nr_of_steps=1000)

print(
    "Final residual (should be close to 0):", mehp_fb.get_displacement_residual_norm()
)

# TODO: apply conversion factors
print(
    "Gamma factor (lacks the conversion factor): ",
    np.sum(mehp_fb.get_gamma_factors(b02=1.0)),
)
