"""
Maximum Entropy Homogenization Procedure (MEHP)
===============================================

This example demonstrates how to use the "original" MEHP implementation,
the one without slip-links or entanglements.

"""

import os

import numpy as np

from pylimer_tools.io.bead_spring_parameter_provider import get_parameters_for_polymer
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

# Prepare parameters for conversion factors
params = get_parameters_for_polymer("PDMS")
r02_slope = params.get("R02")
r02_slope_magnitude = r02_slope.to(params.get("distance_units") ** 2).magnitude
kbt = params.get("T") * params.get("kb")
gamma_conversion_factor = (
    (kbt / ((params.get("distance_units")) ** 3)).to("MPa").magnitude
)

# 1. MEHPForceRelaxation with linear force potential
mehp_relax = MEHPForceRelaxation(universe)
force_evaluator = SimpleSpringMEHPForceEvaluator()
mehp_relax.set_force_evaluator(force_evaluator)
while mehp_relax.requires_another_run():
    mehp_relax.run_force_relaxation()
shear_modulus = (
    gamma_conversion_factor
    * np.sum(mehp_relax.get_gamma_factors(r02_slope_magnitude))
    / universe.get_volume()
)
print(
    "Phantom Shear Modulus [MPa]: ",
    shear_modulus,
)

# 2. MEHPForceRelaxation with Langevin force potential
mehp_relax = MEHPForceRelaxation(universe)
force_evaluator = NonGaussianSpringForceEvaluator()
mehp_relax.set_force_evaluator(force_evaluator)
while mehp_relax.requires_another_run():
    mehp_relax.run_force_relaxation()
shear_modulus = (
    gamma_conversion_factor
    * np.sum(mehp_relax.get_gamma_factors(r02_slope_magnitude))
    / universe.get_volume()
)
print(
    "Phantom Shear Modulus [MPa]: ",
    shear_modulus,
)
