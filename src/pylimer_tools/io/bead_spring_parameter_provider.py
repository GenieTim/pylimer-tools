import math
from enum import Enum
from typing import Optional, Union

import pint
from pint import Quantity, UnitRegistry

from pylimer_tools.io.unit_styles import UnitStyleFactory


class ParameterType(Enum):
    """
    Enum for parameter types used in the Parameters class.

    Attributes:
        GAUSSIAN: Gaussian chain model parameters for ideal bead-spring polymer chains
        KG_LJ: Kremer-Grest model parameters with Lennard-Jones units
        KUHN: Kuhn segment model parameters
    """

    GAUSSIAN = "Gaussian"
    KG_LJ = "Kremer-Grest, Lennard-Jones"
    KUHN = "Kuhn"


class Parameters:
    def __init__(self, data: dict, ureg: UnitRegistry, name: str = ""):
        """
        Initialize a Parameters object with polymer-related data.

        This constructor initializes a Parameters object with a dictionary of polymer-related
        quantities, a unit registry for unit conversions, and an optional name.

        Parameters:
        ----------
        data : dict
            A dictionary containing polymer-related quantities. Each value in the dictionary
            must be a Pint Quantity object. The dictionary must include the keys "Mw", "<b>",
            "rho" and "<b^2>" as these are required parameters.

        ureg : UnitRegistry
            A Pint UnitRegistry object used for unit conversions within the class.

        name : str, optional
            An optional name for this set of parameters. Defaults to an empty string.

        Raises:
        -------
        AssertionError
            If any value in the data dictionary is not a Quantity object, or if any of the
            required parameters ("Mw", "<b>", "<b^2>", "rho") are missing from the data dictionary.

        """
        for key, value in data.items():
            assert isinstance(value, Quantity), (
                f"Invalid value for parameter {key}: {value}"
            )

        # validate parameters
        required_params = ["Mw", "<b>", "<b^2>", "rho"]
        for param in required_params:
            assert param in data, f"Missing required parameter: {param}"

        assert data["Mw"].check(
            "[mass]/[substance]") or data["Mw"].check("[mass]")
        assert data["<b>"].check("[length]")
        assert data["<b^2>"].check("[length] ** 2")
        assert data["rho"].check("[mass]/[volume]") or data["rho"].check(
            "[substance]/[volume]"
        )

        if "R02" not in data:
            data["R02"] = data["<b^2>"]

        if "kb" not in data:
            data["kb"] = 1.380649e-23 * ureg.joule / ureg.kelvin

        self.data = data
        self.name = name
        self.ureg = ureg

    def get_unit_registry(self) -> UnitRegistry:
        return self.ureg

    def get_name(self) -> str:
        return self.name

    def get(self, name: str) -> Quantity:
        """
        Returns the Quantity associated with the given parameter.
        """
        if name in self.data:
            return self.data[(name)]
        synonyms = [
            ["Mw", "bead_mass", "bead_molar_mass"],
            ["rho", "density"],
            ["Ge", "g_e_1", "entanglement_modulus"],
            ["T", "temperature"],
        ]
        for synonym_group in synonyms:
            if name in synonym_group:
                for synonym in synonym_group:
                    if synonym in self.data:
                        return self.data[synonym]
        raise ValueError(f"Unknown parameter: {name}")

    def get_kappa(self) -> Quantity:
        return 3 * self.get("kb") * self.get("T") / self.get("R02")

    def get_sampling_cutoff(self) -> float:
        cutoff = Quantity(
            -0.45 * math.log(self.get("Ge").to("MPa").magnitude) + 1.97,
            "nm",
        )
        distance_units_to_nm = (
            (1 * self.get("distance_units")).to("nanometer").magnitude
        )
        return cutoff.magnitude / distance_units_to_nm

    def get_bead_density(self) -> float:
        """
        Returns the number of segments/beads per unit volume of this particular parameter set.
        """
        density = self.get("rho") / self.get("Mw")
        density *= self.get("distance_units") ** 3
        density = density.to_reduced_units()
        # if the density unit involves a mole, convert to number of beads
        if density.check("[substance]"):
            return density.to("mol").magnitude * 6.02214076e23
        else:
            # assert dimensionless density
            assert density.check("[dimensionless]")
            return density.magnitude

    def get_base_distance_units(self) -> pint.Quantity:
        return self.get("distance_units")

    def get_entanglement_density(
            self, g_e: Optional[pint.Quantity] = None) -> float:
        """
        Returns the number of entanglements per unit volume of this particular parameter set.

        Parameters:
        - g_e: the entanglement modulus
        """
        nm_to_actual_units = (
            (1 * self.ureg.nanometer).to(self.get("distance_units")).magnitude
        )

        if g_e is None:
            g_e = self.get("Ge")

        return (g_e / (self.get("kb") * self.get("T"))).to("nm^-3").magnitude / (
            nm_to_actual_units**3
        )

    def get_gamma_conversion_factor(self) -> Quantity:
        """
        Returns the conversion factor for shear modulus from the gamma factors
        used in the force balance method to MPa.

        Apply as:

        .. code-block:: python

            shear_modulus = gamma_conversion_factor * sum(gamma_factors) / universe.get_volume()

        where `gamma_factors` are the gamma factors calculated from the force balance method,
        using the appropriate `b02`, which can be obtained from the parameters as such:

        .. code-block:: python
            b02 = params.get("R02").to(params.get("distance_units") ** 2).magnitude
        """
        kbt = self.get("T") * self.get("kb")
        gamma_conversion_factor = (
            kbt / ((self.get("distance_units")) ** 3)).to("MPa")
        return gamma_conversion_factor  # type: ignore

    def get_fb_stress_conversion(self) -> float:
        return (self.get_kappa() / (1 * self.get("distance_units"))
                ).to("MPa").magnitude


def assemble_gaussian_parameters_from_kuhn(
    ureg: UnitRegistry,
    kuhn_length: Quantity,
    kuhn_mass: Quantity,
    density: Quantity,
    entanglement_modulus: Quantity,
    temperature: Union[Quantity, None] = None,
    name: str = "",
) -> Parameters:
    """
    Assembles a Parameters object from Kuhn length, mass, density, entanglement modulus, and an optional temperature.

    Parameters:
    - ureg (UnitRegistry): A Pint UnitRegistry instance for unit conversion.
    - kuhn_length (Quantity): The Kuhn length of the polymer chain.
    - kuhn_mass (Quantity): The Kuhn mass of the polymer chain.
    - density (Quantity): The density of the polymer solution.
    - entanglement_modulus (Quantity): The entanglement modulus of the polymer.
    - temperature (Quantity, optional): The temperature of the polymer solution. Defaults to 298 Kelvin.
    - name (str, optional): An optional name for this set of parameters. Defaults to an empty string.

    Returns:
    - Parameters: A Parameters object containing the calculated parameters.
    """
    alpha = 3 * math.pi / 8.0
    assert kuhn_length.check("[length]")
    assert kuhn_mass.check("[mass]/[substance]")
    assert density.check("[mass]/[volume]")
    assert entanglement_modulus.check("[pressure]")
    return Parameters(
        {
            "R02": alpha * ((kuhn_length / alpha) ** 2),
            "Mw": kuhn_mass / alpha,
            "<b>": kuhn_length / alpha,
            "<b^2>": alpha * ((kuhn_length / alpha) ** 2),
            "T": temperature if temperature is not None else 298 * ureg.kelvin,
            "kb": 1.380649e-23 * ureg.joule / ureg.kelvin,
            "distance_units": 1 * ureg.nanometer,
            "rho": density,
            "Ge": entanglement_modulus,
        },
        ureg,
        name,
    )


def get_parameters_for_polymer(
    polymer_name: str, parameter_type: ParameterType
) -> Parameters:
    """
    Returns a Parameters object for the specified polymer and parameter type.
    """
    if parameter_type == ParameterType.GAUSSIAN:
        return get_gaussian_parameters_for_polymer(polymer_name)
    elif parameter_type == ParameterType.KG_LJ:
        return get_kg_lj_parameters_for_polymer(polymer_name)
    elif parameter_type == ParameterType.KUHN:
        return get_kuhn_parameters_for_polymer(polymer_name)
    else:
        raise ValueError(f"Unsupported parameter type: {parameter_type}")


def get_kg_lj_parameters_for_polymer(polymer_name: str) -> Parameters:
    """
    Returns a Parameters object containing
    the parameters for a Kremer-Grest model with Lennard-Jones units
    for the specified polymer.
    """
    unit_style_factory = UnitStyleFactory()
    unit_style = unit_style_factory.get_unit_style("lj", polymer=polymer_name)
    ureg = unit_style.get_underlying_unit_registry()

    row = _get_relevant_everaers_row(polymer_name)

    bead_mass = row["Mb"] * ureg("g/mol")

    return Parameters(
        {
            "distance_units": 1 * unit_style.distance,
            "Mw": bead_mass,
            "<b>": 0.975 * ureg.sigma,
            "<b^2>": (
                row["R_to_2_over_M_c"] *
                ureg("angstrom^2 * mol / g") * bead_mass
            ).to(unit_style.distance**2),
            "rho": unit_style.density,
            "T": row["T_ref"] * ureg.kelvin,
            "kb": unit_style.kb,
            "Ge": row["G_e"] * ureg.megapascal,
        },
        ureg,
        name="kg-lj-{}".format(polymer_name),
    )


def get_kuhn_parameters_for_polymer(
    polymer_name: str, ureg: Union[UnitRegistry, None] = None
) -> Parameters:
    """
    Returns a Parameters object containing
    the parameters for a Kuhn segment model
    for the specified polymer.
    """
    if ureg is None:
        ureg = UnitRegistry()

    row = _get_relevant_everaers_row(polymer_name)

    density = row["rho_bulk"] * ureg("g/cm^3")  # kg/cm^3
    temperature = row["T_ref"] * ureg("K")  # Kelvin
    ge_1 = row["G_e"] * ureg("MPa")  # MPa
    kuhn_length = row["l_K"] * ureg("angstrom")  # °A
    kuhn_bead_mass = row["M_k"] * ureg("g/mol")  # kg/mol

    return Parameters(
        {
            "Mw": kuhn_bead_mass,
            "<b>": kuhn_length.to("nm"),
            "<b^2>": (
                row["R_to_2_over_M_c"] *
                ureg("angstrom^2 * mol / g") * kuhn_bead_mass
            ).to("nm^2"),
            "rho": density,
            "T": temperature,
            "Ge": ge_1,
            "distance_units": 1 * ureg.nanometer,
        },
        ureg,
        name="kuhn-si-{}".format(polymer_name),
    )


def get_gaussian_parameters_for_polymer(
    polymer_name: str, ureg: Union[UnitRegistry, None] = None
) -> Parameters:
    """
    Returns a Parameters object containing
    the parameters for a Gaussian chain model,
    (the ideal bead-spring polymer chain model)
    for the specified polymer.
    """
    if ureg is None:
        ureg = UnitRegistry()

    row = _get_relevant_everaers_row(polymer_name)

    density = row["rho_bulk"] * ureg("g/cm^3")  # kg/cm^3
    # rng.uniform(0.7, 1.5, 1)[0] * ureg("kg/cm^3")  # kg/cm^3
    temperature = row["T_ref"] * ureg("K")  # Kelvin
    ge_1 = row["G_e"] * ureg("MPa")  # MPa
    kuhn_length = row["l_K"] * ureg("angstrom")  # °A
    kuhn_bead_mass = row["M_k"] * ureg("g/mol")  # kg/mol

    return assemble_gaussian_parameters_from_kuhn(
        ureg,
        kuhn_length=kuhn_length.to("nm"),  # type: ignore
        kuhn_mass=kuhn_bead_mass,
        density=density,
        entanglement_modulus=ge_1,
        temperature=temperature,
        name="si-{}".format(polymer_name),
    )


def get_supported_polymer_names() -> list[str]:
    """
    Returns a list of supported polymer names based on the Everaers et al. data.
    """
    from .polymer_data import get_available_polymers

    # Convert to the format used by this module
    available_polymers = get_available_polymers()
    return [
        "".join(filter(str.isalnum, str(name))).lower() for name in available_polymers
    ]


def _get_relevant_everaers_row(polymer_name: str) -> dict:
    """
    Returns the relevant row from the Everaers et al. data for the specified polymer name.
    """
    from .polymer_data import get_polymer_by_name

    try:
        polymer_data = get_polymer_by_name(polymer_name)
        return polymer_data.to_dict()
    except ValueError:
        # Fallback to old method for backward compatibility
        unit_style_factory = UnitStyleFactory()
        everaers_data = unit_style_factory.get_everares_et_al_data()

        for _, row in everaers_data.iterrows():
            if (
                "".join(filter(str.isalnum, polymer_name)).lower()
                == "".join(filter(str.isalnum, str(row["name"]))).lower()
            ):
                return row.to_dict()

        raise ValueError(
            f"Polymer '{polymer_name}' not found in Everaers et al. data.")
