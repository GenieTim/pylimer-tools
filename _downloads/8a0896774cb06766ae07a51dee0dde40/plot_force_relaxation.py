"""
MEHP Force Relaxation
=====================

This example demonstrates how to use
the "original" MEHP implementation (see :cite:t:`gusev_numerical_2019`)
without slip-links or entanglements.

"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pylimer_tools.io.bead_spring_parameter_provider import (
    ParameterType,
    get_parameters_for_polymer,
)
from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import (
    MEHPForceRelaxation,
    NonGaussianSpringForceEvaluator,
    SimpleSpringMEHPForceEvaluator,
    Universe,
    MEHPForceEvaluator,
)


# Example: Custom Force Evaluator (see also tests)
# This particular example is insofar a bad example, as it is already
# implemented faster (no Python<>C++ back & forth every time) in the
# SimpleSpringMEHPForceEvaluator.
class CustomForceEvaluator(MEHPForceEvaluator):
    """
    A custom force evaluator that implements a simple harmonic spring.
    """

    def __init__(self, kappa=1.0):
        super().__init__()
        self.kappa = kappa

    def evaluate_force_set_gradient(
            self, n, spring_distances, compute_gradient):
        network = self.network
        nr_springs = network.nr_of_springs
        force = 0.0
        for i in range(nr_springs):
            spring_vec = spring_distances[3 * i: 3 * i + 3]
            r_squared = np.sum(spring_vec**2)
            contour_length = network.spring_contour_length[i]
            force += r_squared / contour_length
        force *= 0.5 * self.kappa
        gradient = None
        if compute_gradient:
            gradient = [0.0] * n
            nrOfDim = 2 if self.is_2d else 3
            for j in range(nr_springs):
                a = network.spring_index_a[j]
                b = network.spring_index_b[j]
                contour_length = network.spring_contour_length[j]
                for dir_idx in range(nrOfDim):
                    spring_dist = spring_distances[3 * j + dir_idx]
                    grad_term = spring_dist * self.kappa / contour_length
                    gradient[3 * a + dir_idx] += grad_term
                    gradient[3 * b + dir_idx] -= grad_term
        return (force, gradient)

    def evaluate_stress_contribution(
            self, spring_distances, i, j, spring_index):
        network = self.network
        contour_length = network.spring_contour_length[spring_index]
        return self.kappa * spring_distances[i] * \
            spring_distances[j] / contour_length

    def prepare_for_evaluations(self):
        pass


# Load your network (replace with your file)
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

# Prepare parameters for conversion factors
params = get_parameters_for_polymer(
    "PDMS", parameter_type=ParameterType.GAUSSIAN)
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
shear_modulus_linear = (
    gamma_conversion_factor
    * np.sum(mehp_relax.get_gamma_factors(r02_slope_magnitude))
    / universe.get_volume()
)
print(f"Phantom Shear Modulus [MPa] (Linear): {shear_modulus_linear:.4f}")

# 2. MEHPForceRelaxation with Langevin force potential
mehp_relax = MEHPForceRelaxation(universe)
force_evaluator = NonGaussianSpringForceEvaluator()
mehp_relax.set_force_evaluator(force_evaluator)
while mehp_relax.requires_another_run():
    mehp_relax.run_force_relaxation()
shear_modulus_langevin = (
    gamma_conversion_factor
    * np.sum(mehp_relax.get_gamma_factors(r02_slope_magnitude))
    / universe.get_volume()
)
print(f"Phantom Shear Modulus [MPa] (Langevin): {shear_modulus_langevin:.4f}")

# 3. MEHPForceRelaxation with Custom Force Evaluator
custom_evaluator = CustomForceEvaluator(kappa=1.0)
mehp_relax = MEHPForceRelaxation(universe)
mehp_relax.set_force_evaluator(custom_evaluator)
while mehp_relax.requires_another_run():
    mehp_relax.run_force_relaxation()
shear_modulus_custom = (
    gamma_conversion_factor
    * np.sum(mehp_relax.get_gamma_factors(r02_slope_magnitude))
    / universe.get_volume()
)
print(
    f"Phantom Shear Modulus [MPa] (Custom Evaluator): {
        shear_modulus_custom:.4f}")

# %%
# Comparison of Force Evaluators
# ------------------------------
# We visualize the predicted equilibrium phantom shear modulus for each
# force evaluator.

labels = [
    "Simple Spring\n(Linear)",
    "Non-Gaussian\n(Langevin)",
    "Custom\n(Harmonic)"]
moduli = [shear_modulus_linear, shear_modulus_langevin, shear_modulus_custom]
colors = ["#4C72B0", "#55A868", "#DD8452"]

plt.figure(figsize=(7, 5))
bars = plt.bar(
    labels, moduli, color=colors, edgecolor="black", linewidth=0.8, width=0.5
)

for bar in bars:
    yval = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2.0,
        yval + 0.0001,
        f"{yval:.4f} MPa",
        ha="center",
        va="bottom",
        fontsize=10,
    )

plt.ylabel("Phantom Shear Modulus [MPa]", fontsize=11)
plt.title("Comparison of Force Evaluators in MEHP Force Relaxation", fontsize=12)
plt.grid(True, axis="y", linestyle=":", alpha=0.6)
plt.tight_layout()
plt.show()
