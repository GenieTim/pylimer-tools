#!/usr/bin/env python
"""
Unit Conversion
===============

Convert `LAMMPS units <https://docs.lammps.org/units.html>`_ to SI or other unit systems:


**Supported unit styles:**

- ``real`` - kcal/mol, Angstrom, fs
- ``metal`` - eV, Angstrom, ps
- ``si`` - kg, m, s
- ``nano`` - attogram, nanometer, nanosecond
- ``lj`` - Lennard-Jones reduced units (requires polymer specification)

"""

import matplotlib.pyplot as plt
from pylimer_tools.io.unit_styles import UnitStyleFactory

# Create unit style factory
unit_factory = UnitStyleFactory()

# Get unit conversion for specific LAMMPS unit style
real_units = unit_factory.get_unit_style("real")
metal_units = unit_factory.get_unit_style("metal")

# For LJ units, specify the polymer
lj_units = unit_factory.get_unit_style("lj", polymer="pdms")

# Convert pressure from LAMMPS units to SI
pressure_real_pa = (1 * real_units.pressure).to("Pa").magnitude  # type: ignore
pressure_lj_pa = (1 * lj_units.pressure).to("Pa").magnitude  # type: ignore

print(f"1 'real' pressure (atm) in Pascal: {pressure_real_pa:.2f} Pa")
print(f"1 'lj' pressure in Pascal: {pressure_lj_pa:.2f} Pa")

# %%
# Comparison Across Unit Styles
# -----------------------------
# We compare the SI equivalent of 1 pressure unit across all supported
# LAMMPS styles:

styles = ["real", "metal", "si", "nano", "lj"]
pressures_in_pa = []
for style in styles:
    u = (
        unit_factory.get_unit_style(style, polymer="pdms")
        if style == "lj"
        else unit_factory.get_unit_style(style)
    )
    pressures_in_pa.append((1 * u.pressure).to("Pa").magnitude)

plt.figure(figsize=(7, 4.5))
bars = plt.bar(
    styles,
    pressures_in_pa,
    color="#4C72B0",
    edgecolor="black",
    linewidth=0.8,
    width=0.5,
)
plt.yscale("log")
plt.xlabel("LAMMPS Unit Style", fontsize=11)
plt.ylabel("1 Unit of Pressure [Pa] (log scale)", fontsize=11)
plt.title(
    "Pressure Conversion to SI (Pa) Across LAMMPS Unit Styles",
    fontsize=12)
plt.grid(True, which="both", axis="y", linestyle=":", alpha=0.6)
plt.tight_layout()
plt.show()
