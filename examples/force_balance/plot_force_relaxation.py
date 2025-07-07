"""
Maximum Entropy Homogenization Procedure (MEHP)
===============================================

This example demonstrates how to use the "original" MEHP implementation,
the one without slip-links or entanglements.

"""

import os

import numpy as np

from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import (
    MEHPForceRelaxation,
    NonGaussianSpringForceEvaluator,
    SimpleSpringMEHPForceEvaluator,
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

# 1. MEHPForceRelaxation with linear force potential
mehp_relax = MEHPForceRelaxation(universe)
force_evaluator = SimpleSpringMEHPForceEvaluator()
mehp_relax.set_force_evaluator(force_evaluator)
mehp_relax.run_force_relaxation()
print(
    "Gamma factor (lacks the conversion factor): ",
    np.sum(mehp_relax.get_gamma_factors(b0_squared=1.0)),
)

# 2. MEHPForceRelaxation with Langevin force potential
mehp_relax = MEHPForceRelaxation(universe)
force_evaluator = NonGaussianSpringForceEvaluator()
mehp_relax.set_force_evaluator(force_evaluator)
mehp_relax.run_force_relaxation()
print(
    "Gamma factor (lacks the conversion factor): ",
    np.sum(mehp_relax.get_gamma_factors(b0_squared=1.0)),
)
