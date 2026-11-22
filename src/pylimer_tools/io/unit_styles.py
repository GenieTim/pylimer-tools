import warnings
from typing import Any, Dict, Final, List

from pint import Quantity, UnitRegistry
from .polymer_data import (
    PolymerData,
    PolymerDataFrame,
    get_available_polymers,
    get_polymer_by_name,
    load_everaers_et_al_data,
)


class UnitStyle(object):
    """
    UnitStyle: a collection of units of a particular LAMMPS unit style,
    but in SI units
    (i.e.: use this to convert your LAMMPS output data to SI units).

    Example usage:

    .. code:: python

        unit_style_factory = UnitStyleFactory()
        unit_style = unit_style_factory.get_unit_style(
            "lj", polymer="pdms", warning=False, accept_mol=True)

        # multiply with the following factor to convert LJ stress to SI units,
        # namely MPa in this example:
        lj_stress_to_si_conversion_factor = (1.*unit_style.pressure).to("MPa").magnitude
    """

    def __init__(self, unit_configuration: Dict[str, Any], ureg: UnitRegistry) -> None:
        """
        Initialize a UnitStyle object.

        :param unit_configuration: Dictionary containing unit definitions
        :type unit_configuration: Dict[str, Any]
        :param ureg: Pint unit registry to use for unit conversions
        :type ureg: UnitRegistry
        """
        self.unit_configuration = unit_configuration
        self.underlying_unit_registry = ureg
        # add some auxiliary constants
        self.unit_configuration["kb"] = 1.381e-23 * ureg("joule/kelvin")
        if "volume" not in self.unit_configuration:
            self.unit_configuration["volume"] = self.unit_configuration["distance"] ** 3

    def get_underlying_unit_registry(self) -> UnitRegistry:
        """
        Get the underlying Pint unit registry.

        :return: The unit registry used by this UnitStyle
        :rtype: UnitRegistry
        """
        return self.underlying_unit_registry

    def get_base_unit_of(self, property: str) -> Quantity:
        """
        Returns the conversion factor from the unit style to SI units.

        :param property: The property name to get the unit for (e.g., "mass", "distance")
        :type property: str
        :return: The unit object for the requested property
        :rtype: pint.Quantity

        Example usage:

        .. code:: python

          units = get_unit_style("lj")
          mass_in_si = mass_in_lj * units.get_base_unit_of("mass")
        """
        translator = {
            "dynamic viscosity": "viscosity",
            "electric field": "electric_field",
        }
        property = property.lower()
        if property in translator.keys():
            property = translator[property]
        return self.unit_configuration[property]

    def __getattr__(self, property: str) -> Quantity:
        """
        Shorthand access for :func:`~pylimer_tools.io.unit_styles.UnitStyle.get_base_unit_of`.

        :param property: The property name to get the unit for
        :type property: str
        :return: The unit object for the requested property
        :rtype: pint.Quantity

        Example usage:

        .. code:: python

          units = get_unit_style("lj")
          mass_with_units = mass_in_lj * units.mass
        """
        if property.lower() in self.unit_configuration.keys():
            return self.unit_configuration[property.lower()]
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{property}'"
        )


class UnitStyleFactory(object):
    """
    This is a factory to get multiple instances of different
    :obj:`~pylimer_tools.io.UnitStyle`
    using the same UnitRegistry, such that they are compatible.
    """

    def __init__(self) -> None:
        """
        Initialize the UnitStyleFactory with a new UnitRegistry.
        """
        self.ureg = UnitRegistry()

    def get_unit_registry(self) -> UnitRegistry:
        """
        Get the underlying unit registry.

        :return: The unit registry used by this factory
        :rtype: UnitRegistry
        """
        return self.ureg

    def get_everares_et_al_data(self) -> PolymerDataFrame:
        """
        Load the Everaers et al. (2020) unit properties data.

        :return: PolymerDataFrame containing polymer properties from Everaers et al.
        :rtype: PolymerDataFrame
        """
        return load_everaers_et_al_data()

    def get_available_polymers(self) -> List[str]:
        """
        List all available polymers for which we have lj unit conversions.

        :return: List of polymer names
        :rtype: List[str]
        """
        return get_available_polymers()

    def get_unit_style(
        self, unit_type: str, dimension: int = 3, **kwargs: Any
    ) -> UnitStyle:
        """
        Get a UnitStyle instance corresponding to the unit system requested.

        :param unit_type: The unit type, e.g. "lj", "nano", "real", "si", etc.
        :type unit_type: str
        :param dimension: The dimension of the box
        :type dimension: int
        :param kwargs: Additional arguments required for certain unit styles
        :type kwargs: dict
        :return: A UnitStyle object for the requested unit system
        :rtype: UnitStyle
        :raises ValueError: If required parameters are missing
        :raises NotImplementedError: If the requested unit type is not implemented

        For LJ units, you must specify the polymer using the `polymer` parameter.

        See also:
            https://docs.lammps.org/units.html

        .. warning::
            Please check the source code of this function to see
            whether the units you need are correctly implemented
        """
        ureg = self.ureg
        elementary_charge: Final = (1.602176634e-19) * ureg.coulomb
        avogadro_constant: Final = 6.02214076e23  # any/mol

        accept_mol = "accept_mol" in kwargs and kwargs["accept_mol"]

        if unit_type == "lj":
            if ("warning" not in kwargs or kwargs["warning"]) and (
                "polymer" in kwargs and not isinstance(kwargs["polymer"], dict)
            ):
                warnings.warn(
                    "LJ unit styles are derived. Reference used: https://doi.org/10.1021/acs.macromol.9b02428"
                )
            if "polymer" not in kwargs:
                raise ValueError(
                    "LJ unit styles are derived. Please specify the polymer to use (`polymer=...`)"
                )
            polymer_data = kwargs["polymer"]
            if isinstance(polymer_data, str):
                polymer_data = get_polymer_by_name(polymer_data)
            elif isinstance(polymer_data, dict):
                polymer_data = PolymerData(polymer_data)
            if not isinstance(polymer_data, (tuple, PolymerData)):
                raise ValueError(
                    "No useable data for this polymer found to use for lj units. Check whether your usage is correct."
                )
            # follow derivation for more accurate results
            sigma_conversion = (
                0.1 * float(polymer_data.l_K) / (0.965 * float(polymer_data.Cb))
            )
            sigma = sigma_conversion * ureg.nanometer
            eps = float(polymer_data.kB_Tref) * 1e-21 * ureg.joule
            # time is most difficult in LJ — let's keep tau
            if "tau" not in ureg:
                ureg.define("tau = 1 * tau")
            # NOTE: The formula in the LAMMPS documentation contains \epsilon_0.
            # BUT: it does not add up in terms of units, so... the implementation here
            # *might* be wrong
            # epsZero = (8.8541878128e-12*ureg.farad/ureg.meter)
            return UnitStyle(
                {
                    "mass": (
                        polymer_data.Mb * ureg("g/mol")
                        if accept_mol
                        else polymer_data.Mb * ureg("g") / avogadro_constant
                    ),
                    "distance": sigma,
                    "time": ureg.tau,
                    "energy": eps,
                    "velocity": sigma / ureg.tau,
                    "force": eps / sigma,
                    "torque": eps,
                    "temperature": polymer_data.T_ref * ureg.kelvin,
                    "pressure": (
                        polymer_data.kB_Tref_over_sigma_to_3 * ureg("MPa")
                        if hasattr(polymer_data, "kB_Tref_over_sigma_to_3")
                        else eps / (sigma**3)
                    ),
                    "viscosity": eps * ureg.tau / (sigma**3),
                    # TODO: The use of elementary charge might not be correct, see
                    # above
                    "charge": elementary_charge,
                    "dipole": elementary_charge * sigma,
                    "electric_field": eps / (elementary_charge * sigma),
                    "density": (polymer_data.rho_bulk * ureg("g/(cm^3)")).to(
                        "g/(nanometer^3)"
                    ),
                    "dt": 0.005 * ureg.tau,
                    "skin": 0.3 * sigma,
                },
                ureg,
            )
        elif unit_type == "real":
            return UnitStyle(
                {
                    "mass": (
                        ureg("g/mol") if accept_mol else ureg("g") / avogadro_constant
                    ),
                    "distance": ureg.angstrom,
                    "time": ureg.femtosecond,
                    "energy": (
                        ureg("kcal/mol")
                        if accept_mol
                        else ureg("kcal") / avogadro_constant
                    ),
                    "velocity": ureg.angstrom / ureg.femtosecond,
                    "force": (
                        ureg("kcal/(mol*angstrom)")
                        if accept_mol
                        else ureg("kcal") / avogadro_constant / ureg.angstrom
                    ),
                    "torque": (
                        ureg("kcal/mol")
                        if accept_mol
                        else ureg("kcal") / avogadro_constant
                    ),
                    "temperature": ureg.kelvin,
                    "pressure": ureg.atmosphere,
                    "viscosity": ureg.poise,
                    "charge": elementary_charge,
                    "dipole": elementary_charge * ureg.angstrom,
                    "electric_field": ureg.volt / ureg.angstrom,
                    "density": ureg.gram / (ureg.centimeter ** (dimension)),
                    "dt": 1.0 * ureg.femtosecond,
                    "skin": 2.0 * ureg.angstrom,
                },
                ureg,
            )
        elif unit_type == "metal":
            return UnitStyle(
                {
                    "mass": (
                        ureg("g/mol") if accept_mol else ureg("g") / avogadro_constant
                    ),
                    "distance": ureg.angstrom,
                    "time": ureg.picosecond,
                    "energy": ureg("eV"),
                    "velocity": ureg.angstrom / ureg.picosecond,
                    "force": ureg("eV/angstrom"),
                    "torque": ureg("eV"),
                    "temperature": ureg.kelvin,
                    "pressure": ureg.bar,
                    "viscosity": ureg.poise,
                    "charge": elementary_charge,
                    "dipole": elementary_charge * ureg.angstrom,
                    "electric_field": ureg.volt / ureg.angstrom,
                    "density": ureg.gram / (ureg.centimeter ** (dimension)),
                    "dt": 0.001 * ureg.picosecond,
                    "skin": 2.0 * ureg.angstrom,
                },
                ureg,
            )
        elif unit_type == "si":
            return UnitStyle(
                {
                    "mass": ureg.kilogram,
                    "distance": ureg.meter,
                    "time": ureg.second,
                    "energy": ureg.joule,
                    "velocity": ureg.meter / ureg.second,
                    "force": ureg.newton,
                    "torque": ureg.newton * ureg.meter,
                    "temperature": ureg.kelvin,
                    "pressure": ureg.pascal,
                    "viscosity": ureg.pascal * ureg.second,
                    "charge": ureg.coulomb,
                    "dipole": ureg.coulomb * ureg.meter,
                    "electric_field": ureg.volt / ureg.meter,
                    "density": ureg.kilogram / (ureg.meter ** (dimension)),
                    "dt": 1e-8 * ureg.second,
                    "skin": 0.001 * ureg.meter,
                },
                ureg,
            )
        elif unit_type == "nano":
            return UnitStyle(
                {
                    "mass": ureg.attogram,
                    "distance": ureg.nanometer,
                    "time": ureg.nanosecond,
                    "energy": ureg.attogram
                    * (ureg.nanometer**2)
                    / (ureg.nanosecond**2),
                    "velocity": ureg.nanometer / ureg.nanosecond,
                    "force": ureg.attogram * ureg.nanometer / (ureg.nanosecond**2),
                    "torque": ureg.attogram
                    * (ureg.nanometer**2)
                    / (ureg.nanosecond**2),
                    "temperature": ureg.kelvin,
                    "pressure": ureg.attogram / (ureg.nanometer * (ureg.nanosecond**2)),
                    "viscosity": ureg.attogram / (ureg.nanometer * (ureg.nanosecond)),
                    "charge": elementary_charge,
                    "dipole": elementary_charge * ureg.nanometer,
                    "electric_field": ureg.volt / ureg.nanometer,
                    "density": ureg.attogram / (ureg.nanometer ** (dimension)),
                    "dt": 1e-8 * ureg.second,
                    "skin": 0.001 * ureg.meter,
                },
                ureg,
            )
        else:
            raise NotImplementedError(
                "Unit type '{}' is not implemented".format(unit_type)
            )
