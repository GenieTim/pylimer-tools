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

from pylimer_tools.io.unit_styles import UnitStyleFactory

# Create unit style factory
unit_factory = UnitStyleFactory()

# Get unit conversion for specific LAMMPS unit style
real_units = unit_factory.get_unit_style("real")
metal_units = unit_factory.get_unit_style("metal")

# For LJ units, specify the polymer
lj_units = unit_factory.get_unit_style("lj", polymer="pdms")

# Convert pressure from LAMMPS units to SI
pressure_si = (1 * real_units.pressure).to("Pa").magnitude  # type: ignore

# Convert temperature (already in Kelvin for 'real' units)
temperature_si = (273.15 * real_units.temperature).to("K").magnitude  # type: ignore

print(f"LJ 1 Pressure in Pascal: {pressure_si}")
