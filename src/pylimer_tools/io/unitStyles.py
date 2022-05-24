import math
import os
import warnings

import pandas as pd
from pint import UnitRegistry


class UnitStyle(object):
    """
    UnitStyle: a collection of units of a particular lammps unit style,
    but in SI units 
    (i.e.: use this to convert your LAMMPS output data to SI units).
    """

    def __init__(self, unitConfiguration: dict, ureg: UnitRegistry):
        self.unitConfiguration = unitConfiguration
        self.underlyingUnitRegistry = ureg
        # add some auxiliary constants
        self.unitConfiguration['kb'] = 1.381e-23*ureg('joule/kelvin')
        if ('volume' not in self.unitConfiguration):
            self.unitConfiguration['volume'] = self.unitConfiguration['distance']**3

    def getUnderlyingUnitRegistry(self):
        return self.underlyingUnitRegistry

    def getBaseUnitOf(self, property: str):
        """
        Returns the conversion factor from the unit style to SI units.

        Example useage:

        .. code:: python

          units = getUnitStyle("lj")
          siMass = ljMass * units.getBaseUnitOf("mass")
        """
        translator = {
            "dynamic viscosity": "viscosity",
            "electric field": "electric_field",
        }
        property = property.lower()
        if (property in translator.keys()):
            property = translator[property]
        return self.unitConfiguration[property]

    def __getattr__(self, property: str):
        """
        Shorthand access for :func:`~pylimer_tools.io.unitStyles.UnitStyle.getBaseUnitOf`.

        Example useage:

        .. code:: python

          units = getUnitStyle("lj")
          siMass = ljMass * units.mass
        """
        if (property.lower() in self.unitConfiguration.keys()):
            return self.unitConfiguration[property.lower()]


class UnitStyleFactory(object):
    """
    This is a factory to get multiple instances of different 
    :obj:`~pylimer_tools.io.UnitStyle`
    using the same UnitRegistry, such that they are compatible.
    """

    def __init__(self):
        self.ureg = UnitRegistry()

    def getUnitRegistry(self):
        return self.ureg

    def getUnitStyle(self, unitType: str, dimension: int = 3, **kwargs) -> UnitStyle:
        """
        Get a UnitStyle instance corresponding to the unit system requested.

        See also:
          - https://docs.lammps.org/units.html

        Caution:
          - Please check the source code of this function to see 
            whether the units you need are correctly implemented

        Arguments:
          - unitType (str): The unit type, e.g. "lj", "nano", "real", "si", ...
          - dimension (int): The dimension of the box
          - **kwargs: additional arguments required force certain unit styles (e.g. lj)
        """
        ureg = self.ureg
        ELEMENTARY_CHARGE = (1.602176634e-19)*ureg.coulomb
        AVOGADRO_CONST = 6.02214076e23  # any/mol

        if (unitType == "lj"):
            if (("warning" not in kwargs or kwargs["warning"]) and ("polymer" in kwargs and not isinstance(kwargs["polymer"], dict))):
                warnings.warn(
                    "LJ unit styles are derived. Reference used: https://doi.org/10.1021/acs.macromol.9b02428")
            if ("polymer" not in kwargs):
                raise ValueError(
                    "LJ unit styles are derived. Please specify the polymer to use (`polymer=...`)")
            polymerData = kwargs["polymer"]
            if (isinstance(polymerData, str)):
                allPolymerData = pd.read_excel(
                    os.path.join(os.path.dirname(__file__), "..", "data/everaers_et_al_unit_properties.xlsx"))
                for row in allPolymerData.itertuples():
                    if (''.join(filter(str.isalnum, polymerData)).lower() == ''.join(filter(str.isalnum, row.name)).lower()):
                        polymerData = row
                        break
            if (not isinstance(polymerData, dict) and not isinstance(polymerData, tuple)):
                raise ValueError(
                    "No useable data for this polymer found to use for lj units. Check whether your useage is correct.")
            ureg.define("sigma = {} * nanometer".format(polymerData.sigma))
            ureg.define("eps = {}e-21 joule".format(polymerData.kB_Tref))
            # time is most difficult in lj — let's keep tau
            ureg.define('tau = 1 * tau')
            # NOTE: the formula in the LAMMPS documentaion contains \epsilon_0.
            # BUT: it does not add up in terms of units, so... the implementation here
            # *might* be wrong
            # epsZero = (8.8541878128e-12*ureg.farad/ureg.meter)
            return UnitStyle({
                'mass': polymerData.M_k*ureg('g/mol') if 'accept_mol' in kwargs else polymerData.M_k*ureg('g')/AVOGADRO_CONST,
                'distance': ureg.sigma,
                'time': ureg.tau,
                'energy': ureg.eps,
                'velocity': ureg.sigma / ureg.tau,
                'force': ureg.eps/(ureg.sigma),
                'torque': ureg.eps,
                'temperature': polymerData.T_ref*ureg.kelvin,
                'pressure': polymerData.kB_Tref_over_sigma_to_3*ureg("MPa") if hasattr(polymerData, "kB_Tref_over_sigma_to_3") else ureg.eps/(ureg.sigma**(3)),
                'viscosity': ureg.eps*ureg.tau/(ureg.sigma**(3)),
                # TODO: the use of elementary charge might not be correct, see above
                'charge': ELEMENTARY_CHARGE,
                'dipole': ELEMENTARY_CHARGE*ureg.sigma,
                'electric_field': ureg.eps/(ELEMENTARY_CHARGE*ureg.sigma),
                'density': polymerData.M_k*ureg('g/mol')/(ureg.sigma**(dimension)) if 'accept_mol' in kwargs else (polymerData.M_k/AVOGADRO_CONST)*ureg('g')/(ureg.sigma**(dimension)),
                'dt': 0.005*ureg.tau,
                'skin': 0.3*ureg.sigma
            }, ureg)
        elif (unitType == "real"):
            return UnitStyle({
                "mass": ureg('g/mol') if 'accept_mol' in kwargs else ureg('g')/AVOGADRO_CONST,
                "distance": ureg.angstrom,
                "time": ureg.femtosecond,
                "energy": ureg('kcal/mol') if 'accept_mol' in kwargs else ureg('kcal')/AVOGADRO_CONST,
                "velocity": ureg.angstrom/ureg.femtosecond,
                "force": ureg('kcal/(mol*angstrom)') if 'accept_mol' in kwargs else ureg('kcal')/AVOGADRO_CONST/ureg.angstrom,
                "torque": ureg('kcal/mol') if 'accept_mol' in kwargs else ureg('kcal')/AVOGADRO_CONST,
                "temperature": ureg.kelvin,
                "pressure": ureg.atmosphere,
                "viscosity": ureg.poise,
                "charge": ELEMENTARY_CHARGE,
                "dipole": ELEMENTARY_CHARGE*ureg.angstrom,
                "electric_field": ureg.volt/ureg.angstrom,
                "density": ureg.gram/(ureg.meter**(dimension)),
                "dt": 1.0*ureg.femtosecond,
                "skin": 2.0*ureg.angstrom
            }, ureg)
        elif(unitType == "si"):
            return UnitStyle({
                "mass": ureg.kilogram,
                "distance": ureg.meter,
                "time": ureg.second,
                "energy": ureg.joule,
                "velocity": ureg.meter/ureg.second,
                "force": ureg.newton,
                "torque": ureg.newton*ureg.meter,
                "temperature": ureg.kelvin,
                "pressure": ureg.pascal,
                "viscosity": ureg.pascal*ureg.second,
                "charge": ureg.coulomb,
                "dipole": ureg.coulomb*ureg.meter,
                "electric_field": ureg.volt/ureg.meter,
                "density": ureg.kilogram/(ureg.meter**(dimension)),
                "dt": 1e-8*ureg.second,
                "skin": 0.001*ureg.meter
            }, ureg)
        elif (unitType == "nano"):
            return UnitStyle({
                "mass": ureg.attogram,
                "distance": ureg.nanometer,
                "time": ureg.nanosecond,
                "energy": ureg.attogram*(ureg.nanometer**2)/(ureg.nanosecond**2),
                "velocity": ureg.nanometer/ureg.nanosecond,
                "force": ureg.attogram*ureg.nanometer/(ureg.nanosecond**2),
                "torque": ureg.attogram*(ureg.nanometer**2)/(ureg.nanosecond**2),
                "temperature": ureg.kelvin,
                "pressure": ureg.attogram/(ureg.nanometer*(ureg.nanosecond**2)),
                "viscosity": ureg.attogram/(ureg.nanometer*(ureg.nanosecond)),
                "charge": ELEMENTARY_CHARGE,
                "dipole": ELEMENTARY_CHARGE*ureg.nanometer,
                "electric_field": ureg.volt/ureg.nanometer,
                "density": ureg.attogram/(ureg.nanometer**(dimension)),
                "dt": 1e-8*ureg.second,
                "skin": 0.001*ureg.meter
            }, ureg)
        else:
            raise NotImplementedError(
                "Unity type '{}' is not implemented".format(unitType))
