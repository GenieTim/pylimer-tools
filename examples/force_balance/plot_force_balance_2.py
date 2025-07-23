"""
Force Balance 2
===============

MEHPForceBalance2 is a faster implementation of the MEHP force balance method,
which uses static links to model entanglements instead of slip-links.
See :cite:t:`bernhard_phantom_2025`
"""

import os

import numpy as np

from pylimer_tools.io.bead_spring_parameter_provider import (
    Parameters,
    get_parameters_for_polymer,
)
from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import MEHPForceBalance2, Universe

# Get parameters for conversion factors
params = get_parameters_for_polymer("PDMS")
assert isinstance(params, Parameters)

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
mehp_fb = MEHPForceBalance2(
    universe,
    nr_of_entanglements_to_sample=int(
        params.get_entanglement_density() * universe.get_volume()
    ),
    upper_sampling_cutoff=params.get_sampling_cutoff(),
    entanglements_as_springs=False,
)
mehp_fb.run_force_relaxation()

print(
    "Final residual (should be close to 0):", mehp_fb.get_displacement_residual_norm()
)

# apply conversion factors, measure resulting shear modulus
r02_slope = params.get("R02")
r02_slope_magnitude = r02_slope.to(params.get("distance_units") ** 2).magnitude
kbt = params.get("T") * params.get("kb")
gamma_conversion_factor = (
    (kbt / ((params.get("distance_units")) ** 3)).to("MPa").magnitude
)

shear_modulus = (
    gamma_conversion_factor
    * np.sum(mehp_fb.get_gamma_factors(r02_slope_magnitude))
    / universe.get_volume()
)
print("Shear Modulus [MPa]: ", shear_modulus)
