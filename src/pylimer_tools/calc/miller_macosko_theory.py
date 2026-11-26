"""
This module provides access to various computations introduced in the Miller-Macosko theory.
See :cite:t:`miller_new_1976a` and :cite:t:`macosko_new_1976`.

Additional references used in the development of this module include
:cite:t:`langley_elastically_1968`, :cite:t:`miller_average_1978`,
:cite:t:`valles_properties_1979`, :cite:t:`miller_calculation_1979`,
:cite:t:`venkataraman_critical_1989`, :cite:t:`patel_elastic_1992`,
:cite:t:`aoyama_nonlinear_2021a`, :cite:t:`gusev_molecular_2022`,
and :cite:t:`tsimouri_comparison_2024`.

Caution:
      Not all systems are supported yet.
      In particular, for most methods, only A_f and B_2 is supported.
      Also, the systems are mostly assumed to be end-linked and monodisperse.

"""

from __future__ import annotations

import math
import warnings
from typing import Any, Callable, List, Tuple, Union

import pint

from pylimer_tools.calc.structure_analysis import (
    _resolve_chain_functionality_per_type,
    compute_crosslinker_conversion,
    compute_fraction_of_bifunctional_reactive_sites,
    compute_stoichiometric_imbalance,
    count_strands_by_functionality,
    determine_chain_functionality_per_type,
)
from pylimer_tools.io.unit_styles import UnitStyle
from pylimer_tools_cpp import Universe

_KB_SI_J_PER_K: float = 1.380649e-23


def _compute_effective_junction_factor(f: int, alpha: float) -> float:
    """Compute MMT junction prefactor sum_{m=3}^f ((m - 2) / 2) * P(X_{m,f})."""
    gamma_mmt_sum = 0.0
    for m in range(3, f + 1):
        gamma_mmt_sum += (
            (m - 2) / 2
        ) * compute_probability_that_crosslink_is_effective(f, m, alpha)
        assert math.isfinite(gamma_mmt_sum), (
            "Non-finite gamma_mmt_sum computed for f = {}.".format(m)
        )
    return gamma_mmt_sum


def predict_shear_modulus(**kwargs) -> pint.Quantity:
    """
    Predict the shear modulus using MMT Analysis.

    Source:
      - :cite:t:`gusev_molecular_2022`

    :param kwargs: See :func:`~pylimer_tools.calc.miller_macosko_theory.compute_modulus_decomposition`
    :return: G: The predicted shear modulus

    ToDo:
      - Support more than one crosslinker type (as is supported by original formula)
    """
    g_mmt_phantom, g_mmt_entanglement, _, _ = compute_modulus_decomposition(**kwargs)
    return g_mmt_phantom + g_mmt_entanglement  # type: ignore


def _resolve_crosslinker_prefactor(
    crosslinker_number_density: Union[float, None],
    param: dict,
    crosslinker_type: int,
    func_name: str,
) -> float:
    """Resolve the crosslinker prefactor, warning on deprecated weight_fractions fallback."""
    if crosslinker_number_density is not None:
        return crosslinker_number_density
    warnings.warn(
        f"Calling {func_name} without 'crosslinker_number_density' is deprecated; "
        "the dimensionless weight_fractions[crosslinker_type] fallback will be "
        "removed in a future release.",
        DeprecationWarning,
        stacklevel=3,
    )
    weight_fractions = param["weight_fractions"]
    assert isinstance(weight_fractions, dict)
    return weight_fractions.get(crosslinker_type, 0.0)


def predict_number_density_of_junction_points(
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    b2: Union[float, None] = None,
    crosslinker_number_density: Union[float, None] = None,
) -> float:
    """
    Compute the number density of junction points using MMT.

    Source:
      - :cite:t:`aoyama_nonlinear_2021a` (see supporting information for formulae)

    :param network: The network to compute the number density for
    :param crosslinker_type: The atom type to use to split the molecules
    :param functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :param crosslinker_number_density: The initial number density of crosslinkers :math:`[A_f]_0`.
        If ``None``, falls back to ``weight_fractions[crosslinker_type]`` for backwards compatibility
        (deprecated; will be required or derived from ``network`` in a future release).
    :return: mu: The predicted number density of junction points
    """
    required_params = ["functionality_per_type", "p_f_a_out"]
    if crosslinker_number_density is None:
        required_params.append("weight_fractions")

    param = _compute_validate_parameters(
        {**locals()},
        required_params,
    )

    functionality_per_type, alpha = (
        param["functionality_per_type"],
        param["p_f_a_out"],
    )
    assert isinstance(functionality_per_type, dict) and alpha is not None

    prefactor = _resolve_crosslinker_prefactor(
        crosslinker_number_density,
        param,
        crosslinker_type,
        "predict_number_density_of_junction_points",
    )

    if functionality_per_type[crosslinker_type] == 3:
        return prefactor * (1 - alpha) ** 3
    elif functionality_per_type[crosslinker_type] == 4:
        return prefactor * (4 * alpha * (1 - alpha) ** 3 + (1 - alpha) ** 4)
    else:
        raise NotImplementedError(
            "Currently, only crosslinker functionalities of 3 and 4 are supported, {} given.".format(
                functionality_per_type[crosslinker_type]
            )
        )


def predict_number_density_of_network_strands(
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    r: Union[float, None] = None,
    p: Union[float, None] = None,
    b2: Union[float, None] = None,
    crosslinker_number_density: Union[float, None] = None,
) -> float:
    """
    Compute the number density of network strands using MMT.

    Source:
      - :cite:t:`aoyama_nonlinear_2021a` (see supporting information for formulae)

    :param network: The network to compute the number density for
    :param crosslinker_type: The atom type to use to split the molecules
    :param functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
    :param r: The stoichiometric imbalance
    :param p: The extent of reaction in terms of the crosslinkers
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :param crosslinker_number_density: The initial number density of crosslinkers :math:`[A_f]_0`.
        If ``None``, falls back to ``weight_fractions[crosslinker_type]`` for backwards compatibility
        (deprecated; will be required or derived from ``network`` in a future release).
    :return: nu: The predicted number density of network strands
    """
    required_params = ["functionality_per_type", "r", "p", "b2", "crosslinker_type"]
    if crosslinker_number_density is None:
        required_params.append("weight_fractions")

    param = _compute_validate_parameters(
        {**locals()},
        required_params,
    )

    functionality_per_type, r, p, b2, crosslinker_type = (
        param["functionality_per_type"],
        param["r"],
        param["p"],
        param["b2"],
        param["crosslinker_type"],
    )
    assert (
        isinstance(functionality_per_type, dict)
        and r is not None
        and p is not None
        and b2 is not None
        and isinstance(crosslinker_type, int)
    )

    prefactor = _resolve_crosslinker_prefactor(
        crosslinker_number_density,
        param,
        crosslinker_type,
        "predict_number_density_of_network_strands",
    )

    alpha, _ = compute_miller_macosko_probabilities(
        r=r,
        p=p,
        f=functionality_per_type[crosslinker_type],
        b2=b2,
    )

    if functionality_per_type[crosslinker_type] == 3:
        return (3 / 2) * prefactor * (1 - alpha) ** 3
    elif functionality_per_type[crosslinker_type] == 4:
        return prefactor * (6 * alpha * (1 - alpha) ** 3 + 2 * (1 - alpha) ** 4)
    else:
        raise NotImplementedError(
            "Currently, only junction functionalities of 3 and 4 are supported, {} given.".format(
                functionality_per_type[crosslinker_type]
            )
        )


def compute_weight_fraction_of_dangling_chains(
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    weight_fractions: Union[dict, None] = None,
    r: Union[float, None] = None,
    p: Union[float, None] = None,
    b2: Union[float, None] = None,
) -> float:
    """
    Compute the weight fraction of dangling (pendant) strands in infinite network

    Source:
      - Eq. 6.4 in :cite:t:`miller_calculation_1979`

    :param network: The network to compute the weight fraction for
    :param crosslinker_type: The atom type to use to split the molecules
    :param functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
    :param weight_fractions: a dictionary with the weight fraction of each type of atom
    :param r: The stoichiometric imbalance
    :param p: The extent of reaction in terms of the crosslinkers
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :return: weightFraction :math:`\\Phi_d = w_p`: weightDangling/weightTotal
    """
    if network is not None and network.get_nr_of_atoms() == 0:
        return 0

    param = _compute_validate_parameters(
        {**locals()},
        ["functionality_per_type", "weight_fractions", "p_f_a_out", "p_f_b_out"],
    )

    functionality_per_type, weight_fractions, alpha, beta = (
        param["functionality_per_type"],
        param["weight_fractions"],
        param["p_f_a_out"],
        param["p_f_b_out"],
    )
    assert (
        isinstance(functionality_per_type, dict)
        and isinstance(weight_fractions, dict)
        and alpha is not None
        and beta is not None
    )

    w_dangling = 0.0
    for atom_type, weight_fraction in weight_fractions.items():
        if atom_type == crosslinker_type:
            probabilities = compute_probability_that_crosslink_is_dangling(
                functionality_per_type[crosslinker_type], alpha
            )
            for i in range(functionality_per_type[crosslinker_type] - 1):
                probabilities += (
                    compute_probability_that_crosslink_with_degree_is_dangling(
                        functionality_of_monomer=functionality_per_type[
                            crosslinker_type
                        ],
                        degree_of_ineffectiveness=i,
                        p_f_a_out=alpha,
                    )
                    * (i / functionality_per_type[crosslinker_type])
                )
            w_dangling += probabilities * weight_fraction
        elif functionality_per_type[atom_type] == 2:
            w_dangling += (
                weight_fraction
                * compute_probability_that_bifunctional_monomer_is_dangling(beta)
            )
        elif functionality_per_type[atom_type] == 1:
            w_dangling += weight_fraction * (1.0 - beta)
        elif functionality_per_type[atom_type] == 0:
            # Solvent chains cannot be dangling
            w_dangling += 0.0
        else:
            raise NotImplementedError(
                "Currently, only monomeric, bifunctional, and junction functionalities are supported, {} given.".format(
                    functionality_per_type[atom_type]
                )
            )

    return w_dangling


def compute_weight_fraction_of_backbone(
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    weight_fractions: Union[dict, None] = None,
    r: Union[float, None] = None,
    p: Union[float, None] = None,
    b2: Union[float, None] = None,
) -> float:
    """
    Compute the weight fraction of the backbone (elastically effective) strands in an infinite network

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param functionality_per_type: a dictionary with key: atom type, and value: functionality atoms with this type.
    :param weight_fractions: a dictionary with the weight fraction of each type of atom
    :param r: The stoichiometric imbalance
    :param p: The extent of reaction in terms of the crosslinkers
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :return: :math:`\\Phi_{el} = w_e`: weight fraction of network backbone
    """
    if network is not None and network.get_nr_of_atoms() == 0:
        return 0

    param = _compute_validate_parameters(
        {**locals()},
        ["functionality_per_type", "weight_fractions", "p_f_a_out", "p_f_b_out"],
    )

    functionality_per_type, weight_fractions, alpha, beta = (
        param["functionality_per_type"],
        param["weight_fractions"],
        param["p_f_a_out"],
        param["p_f_b_out"],
    )
    assert isinstance(functionality_per_type, dict) and isinstance(
        weight_fractions, dict
    )

    w_elastic = 0.0
    for atom_type, weight_fraction in weight_fractions.items():
        if atom_type == crosslinker_type:
            probabilities = 0.0
            for i in range(2, functionality_per_type[crosslinker_type] + 1):
                probabilities += compute_probability_that_crosslink_is_effective(
                    functionality_of_monomer=functionality_per_type[crosslinker_type],
                    expected_degree_of_effect=i,
                    p_f_a_out=alpha,
                ) * (i / functionality_per_type[crosslinker_type])
            w_elastic += probabilities * weight_fraction
        elif functionality_per_type[atom_type] == 2:
            w_elastic += (
                weight_fraction
                * compute_probability_that_bifunctional_monomer_is_effective(beta)
            )
        elif functionality_per_type[atom_type] <= 1:
            # Monofunctional chains and solvent cannot form elastically effective strands
            w_elastic += 0.0
        else:
            raise NotImplementedError(
                "Currently, only monomeric, bifunctional, and junction functionalities are supported, {} given.".format(
                    functionality_per_type[atom_type]
                )
            )

    return w_elastic


def compute_weight_fraction_of_soluble_material(
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    weight_fractions: Union[dict, None] = None,
    r: Union[float, None] = None,
    p: Union[float, None] = None,
    b2: Union[float, None] = None,
) -> float:
    """
    Compute the weight fraction of soluble material by MMT.

    Source:
      - :cite:t:`patel_elastic_1992`
      - :cite:t:`aoyama_nonlinear_2021a`

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param weight_fractions: a dictionary with key: type, and value: weight fraction of type.
            Pass if you want to omit the network.
    :param functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
          See: :meth:`~pylimer_tools_cpp.Universe.determine_functionality_per_type`.
    :param r: The stoichiometric imbalance
    :param p: The extent of reaction in terms of the crosslinkers
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :return: :math:`W_{sol}` (float): The weight fraction of soluble material according to MMT.
    """
    if network is not None and network.get_nr_of_bonds() == 0:
        return 1.0

    param = _compute_validate_parameters(
        {**locals()},
        ["functionality_per_type", "weight_fractions", "p_f_a_out", "p_f_b_out"],
    )

    functionality_per_type, weight_fractions, alpha, beta = (
        param["functionality_per_type"],
        param["weight_fractions"],
        param["p_f_a_out"],
        param["p_f_b_out"],
    )
    assert (
        isinstance(functionality_per_type, dict)
        and isinstance(weight_fractions, dict)
        and alpha is not None
        and beta is not None
    )

    w_sol = 0
    for key in weight_fractions:
        coefficient = alpha if key == crosslinker_type else beta
        if key not in weight_fractions or math.isclose(
            weight_fractions[key], 0.0, abs_tol=1e-10
        ):
            continue
        w_sol += weight_fractions[key] * (
            math.pow(coefficient, functionality_per_type[key])
        )

    return w_sol


def compute_weight_fraction_of_soluble_material_from_weight_fractions(
    r: float, p: float, f: int, w_f: float, w_g: float, g: int = 2
) -> float:
    """
    Use MMT to compute the weight fraction of soluble material using

    .. math::
        W_{sol} = w_{A_f} P(F_A^{out})^f + w_{B_g} [rpP(F_A^{out})^{f-1}+1-rp]^g

    :param r: The stoichiometric imbalance
    :param p: The extent of reaction in terms of the crosslinkers
    :param f: The functionality of the the crosslinker
    :param w_f: The weight fraction of the crosslinkers
    :param w_g: The weight fraction of ordinary chains
    :param g: The functionality of the ordinary chains
    """
    alpha, _ = compute_miller_macosko_probabilities(r, p, f)
    return w_f * (alpha**f) + w_g * ((r * p * (alpha ** (f - 1)) + 1 - r * p) ** g)


def compute_miller_macosko_probabilities(
    r: float, p: float, f: int, b2: float = 1.0
) -> Tuple[float, float]:
    """
    Compute Macosko and Miller's probabilities :math:`P(F_A)` and :math:`P(F_B)`
    i.e., the probability that a randomly chosen A (crosslink) or B (strand-end),
    respectively, is the start of a finite chain.

    Sources:
        - :cite:t:`macosko_new_1976`
        - :cite:t:`miller_new_1976a`
        - :cite:t:`patel_elastic_1992` (with monofunctional chains, f = 4)
        - :cite:t:`urayama_damping_2004` (with monofunctional chains, f = 3)

    .. note::
        Currently, only systems with B_2, B_1 and A_f are supported.
        If you have a system with other functionality distributions,
        you would need to implement these formulas yourselves.

    :param r: The stoichiometric imbalance
    :param p: The extent of reaction in terms of the crosslinkers
    :param f: The functionality of the crosslinker
    :param b2: The fraction of reactive B sites belonging to bifunctional chains;
            defaults to 1.0 for no monofunctional chains (`0.0` for pure :math:`B_1`).
            Can be computed e.g. as :math:`b_2 = \\frac{2 \\cdot [B_2]}{[B_1] + 2 \\cdot [B_2]}`
    :return: alpha: :math:`P(F_A)`
    :return: beta: :math:`P(F_B)`
    """
    if b2 < 0.0 or b2 > 1.0 or not math.isfinite(b2):
        raise ValueError("b2 must be in [0, 1], got b2 = {}".format(b2))
    if r == 0 or p == 0 or f == 0:
        return 1.0, 1.0

    # Validate parameters and clamp tiny floating-point overshoot at p == p_max
    p = _validate_r_and_p(r, p, f)
    if b2 == 0.0:
        return 1.0, 1.0

    if not (1 / (f - 1) < b2 * (p**2) * r < 1):
        warnings.warn(
            "The resulting P(F_A) is probably unreliable, "
            + "as the detected root does not fulfill the required conditions."
        )

    # End validation.
    # actually do the calculations
    if f == 3:
        alpha = (1 - r * p * p * b2) / (r * p * p * b2)
    elif f == 4:
        alpha = ((1.0 / (r * p * p * b2)) - 3.0 / 4.0) ** (1.0 / 2.0) - (1.0 / 2.0)
    else:
        if not (f > 4):
            raise NotImplementedError(
                "A functionality of {} is not supported.".format(f)
            )
        if b2 * (p**2) * r <= 1.0 / (f - 1):
            alpha = 1.0
        else:
            try:
                from scipy import optimize
            except ImportError:
                raise ImportError(
                    "scipy is required for crosslinker functionality > 4. "
                    "Please install scipy or use f <= 4."
                )

            def fun_to_root_for_alpha(alpha_val: float) -> float:
                return (
                    r * b2 * p**2 * alpha_val ** (f - 1)
                    - alpha_val
                    - r * b2 * (p**2)
                    + 1
                )

            def fun_to_root_for_alpha_prime(alpha_val: float) -> float:
                return -1 + alpha_val ** (f - 2) * (-1 + f) * (p**2) * r * b2

            def fun_to_root_for_alpha_prime2(alpha_val: float) -> float:
                return alpha_val ** (f - 3) * (-2 + f) * (-1 + f) * (p**2) * r * b2

            alpha_sol = optimize.root_scalar(
                fun_to_root_for_alpha,
                bracket=(0, 1),
                method="halley",
                fprime=fun_to_root_for_alpha_prime,
                fprime2=fun_to_root_for_alpha_prime2,
                x0=0.5,
            )
            alpha = alpha_sol.root
    beta = r * p * alpha ** (f - 1) + 1 - r * p
    if alpha > 1 or alpha < 0:
        warnings.warn(
            "The resulting P(F_A) from r = {}, p = {}, b2 = {} for f = {} is probably unreliable, ".format(
                r, p, b2, f
            )
            + "as it will be clipped to [0,1] from {}".format(alpha)
        )
    if beta > 1 or beta < 0:
        warnings.warn(
            "The resulting P(F_B) from r = {}, p = {} , b2 = {}for f = {} is probably unreliable, ".format(
                r, p, b2, f
            )
            + "as it will be clipped to [0,1] from {}".format(beta)
        )
    # TODO: reconsider clipping.
    return max(0, min(alpha, 1)), max(0, min(beta, 1))


def compute_modulus_decomposition(
    network: Union[Universe, None] = None,
    ureg: Union[pint.UnitRegistry, None] = None,
    unit_style: Union[None, UnitStyle] = None,
    crosslinker_type: Union[int, None] = None,
    r: Union[float, None] = None,
    p: Union[float, None] = None,
    f: Union[int, None] = None,
    nu: Union[pint.Quantity, None] = None,
    temperature: Union[pint.Quantity, None] = None,
    functionality_per_type: Union[dict, None] = None,
    g_e_1: Union[pint.Quantity, None] = None,
    b2: Union[float, None] = None,
) -> Tuple[pint.Quantity, pint.Quantity, pint.Quantity, pint.Quantity]:
    """
    Compute four different estimates of the plateau modulus, using MMT, ANM and PNM.

    :param network: The polymer network to do the computation for
    :param ureg: The unit registry to use
    :param unit_style: The unit style to use to have the results in appropriate units
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param r: The stoichiometric imbalance. Optional if network is specified
    :param p: The extent of reaction. Optional if network is specified
    :param f: The functionality of the the crosslinker. Optional if network is specified
    :param nu: The strand number density (nr of reactive strands :math:`[B_2]_0 + [B_1]_0` per volume)
        (ideally with units). Optional if network is specified
    :param temperature: The temperature to compute the modulus at. Default: 298.15 K
        Optional, can be passed to improve performance
    :param functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
        Optional, can be passed to improve performance
    :param g_e_1: The melt entanglement modulus
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2. Defaults to 1.0 if ``network`` is not specified,
        otherwise inferred from ``network``.
    :return: G_MMT_phantom: The phantom contribution to the MMT modulus;
        see also :func:`pylimer_tools.calc.miller_macosko_theory.compute_junction_modulus`
    :return: G_MMT_entanglement: The entanglement contribution to the MMT modulus
    :return: g_anm: The ideal ANM upper-bound estimate (:math:`\\nu k_B T`)
    :return: g_pnm: The ideal PNM upper-bound estimate (:math:`(1 - 2/f) \\nu k_B T`)
    """
    if ureg is None:
        if unit_style is None:
            raise ValueError(
                "Unit style or unit registry must be specified to compute modulus."
            )
        ureg = unit_style.get_underlying_unit_registry()

    if b2 is None and network is None:
        b2 = 1.0

    param = _compute_validate_parameters(
        {**locals()},
        ["f", "nu", "p_f_a_out", "p_f_b_out", "p", "r", "b2"],
    )

    f, nu, alpha, beta, p, r, b2 = (
        param["f"],
        param["nu"],
        param["p_f_a_out"],
        param["p_f_b_out"],
        param["p"],
        param["r"],
        param["b2"],
    )
    assert (
        f is not None
        and nu is not None
        and alpha is not None
        and beta is not None
        and p is not None
        and r is not None
        and b2 is not None
    )

    if temperature is None:
        temperature = (273.15 + 25) * ureg.kelvin  # Temperature in Kelvin
    assert isinstance(temperature, ureg.Quantity) and temperature.check("[temperature]")
    if g_e_1 is None:
        g_e_1 = (
            8.3145  # gas constant, J/(mol*K)
            * temperature.to("kelvin").magnitude  # Temperature in Kelvin
            * 1e-6
            * 94.79281
        ) * ureg("MPa")
        # -> MPa, melt entanglement modulus of PDMS
    assert isinstance(g_e_1, ureg.Quantity) and g_e_1.check("[pressure]")

    # Ensure nu is a pint.Quantity.
    if not isinstance(nu, pint.Quantity):
        raise ValueError("nu must be a pint.Quantity.")

    # If given as molar concentration, convert to number density (molecules
    # per volume)
    if nu.check("[substance]/[volume]"):
        nu = nu.to("mol/meter**3").magnitude * ureg("1/meter**3") * 6.02214076e23

    if not nu.check("1/[volume]"):
        raise ValueError("nu must have dimensionality of number density (1/volume).")

    # Boltzmann constant
    kb = _KB_SI_J_PER_K * ureg.joule / ureg.kelvin
    # affine
    g_anm = nu * kb * temperature
    # phantom
    g_pnm = (1 - 2 / f) * nu * kb * temperature if f != 0 else 0.0
    # MMT:
    # Since nu = [B_2]_0 + [B_1]_0 and b2 = 2[B_2]_0 / (2[B_2]_0 + [B_1]_0),
    # total reactive B-site density is [B]_0 = 2 / (2 - b2) * nu,
    # hence [A_f]_0 = (r / f) * [B]_0 = (2 * r / (f * (2 - b2))) * nu.
    gamma_mmt_sum = _compute_effective_junction_factor(f, alpha)
    gamma_mmt = (2 * r / (f * (2 - b2))) * gamma_mmt_sum if f != 0 else 0.0
    g_mmt_phantom = gamma_mmt * nu * kb * temperature
    # fraction of elastically effective strands.
    g_mmt_entanglement = g_e_1 * compute_trapping_factor(beta=beta)
    # entanglement part. TODO : check adjustment with r
    return g_mmt_phantom, g_mmt_entanglement, g_anm, g_pnm


def compute_extracted_modulus(
    p: float,
    r: float,
    f: int,
    g_e_1: pint.Quantity,
    w_sol: float,
    xlink_concentration_0: pint.Quantity,
    ureg: pint.UnitRegistry,
    temperature: Union[pint.Quantity, None] = None,
    b2: float = 1.0,
) -> pint.Quantity:
    """
    Compute MMT's modulus, assuming the solvent is removed

    :param p: The crosslinker conversion
    :param r: The stoichiometric imbalance
    :param f: The functionality of the crosslinkers
    :param g_e_1: The melt entanglement modulus :math:`G_e(1) = k_B T \\epsilon_e`
    :param w_sol: The soluble fraction (to be removed)
    :param xlink_concentration_0: [A_f]_0, in 1/volume units
    :param ureg: The unit registry to use
    :param temperature: The temperatures; defaults to room temperature
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    """
    if temperature is None:
        if ureg is None:
            raise ValueError(
                "Unit registry must be initialized, or temperature specified."
            )
        temperature = (273.15 + 25) * ureg.kelvin

    junction_part = (1 - w_sol) ** (-1 / 3) * compute_junction_modulus(
        p=p,
        r=r,
        xlink_concentration_0=xlink_concentration_0,
        f=f,
        ureg=ureg,
        temperature=temperature,
        b2=b2,
    )
    entanglement_part = (1 - w_sol) ** (-2) * compute_entanglement_modulus(
        p=p,
        r=r,
        f=f,
        g_e_1=g_e_1,
        temperature=temperature,
        b2=b2,
    )
    return junction_part + entanglement_part


def compute_entanglement_modulus(
    g_e_1: pint.Quantity,
    temperature: Union[pint.Quantity, None] = None,
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    p: Union[float, None] = None,
    r: Union[float, None] = None,
    f: Union[int, None] = None,
    b2: Union[float, None] = None,
    beta: Union[float, None] = None,
) -> pint.Quantity:
    """
    Compute MMT's entanglement contribution to the equilibrium shear modulus, given by
    :math:`k_B T \\epsilon_e T_e`.

    :param g_e_1: The melt entanglement modulus :math:`G_e(1) = k_B T \\epsilon_e`
    :param temperature: Unused (kept for API compatibility; temperature is already included in ``g_e_1``)
    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param p: The crosslinker conversion
    :param r: The stoichiometric imbalance
    :param f: The functionality of the crosslinkers
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :param beta: :math:`P(F_b^{out})`, optional
    :returns: :math:`T_e \\cdot G_e(1)`
    """
    if beta is None:
        param = _compute_validate_parameters(
            {**locals()},
            ["p_f_b_out"],
        )
        beta = param["p_f_b_out"]
        assert beta is not None
    elif not (0 <= beta <= 1):
        raise ValueError("Invalid value for parameter 'beta' (got {}).".format(beta))

    return compute_trapping_factor(beta=beta) * g_e_1


def compute_junction_modulus(
    p: float,
    r: float,
    xlink_concentration_0: pint.Quantity,
    ureg: pint.UnitRegistry,
    f: int,
    temperature: Union[pint.Quantity, None] = None,
    b2: Union[float, None] = None,
) -> pint.Quantity:
    """
    Compute MMT's junction modulus, given by
    :math:`G_{junctions} = k_B T [A_f]_0 \\sum_{m=3}^{f} \\frac{m-2}{2} P(X_{m,f})`.

    :param p: The crosslinker conversion
    :param r: The stoichiometric imbalance
    :param xlink_concentration_0: [A_f]_0, in 1/volume units
    :param ureg: The unit registry to use
    :param f: The functionality of the crosslinkers
    :param temperature: The temperatures; defaults to room temperature (25 °C)
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :return: The junction modulus contribution
    """
    if b2 is None:
        b2 = 1.0
    param = _compute_validate_parameters(
        {**locals()},
        ["p_f_a_out"],
    )

    alpha = param["p_f_a_out"]
    assert alpha is not None
    if temperature is None:
        assert ureg is not None, "Unit registry must be initialized."
        temperature = (273.15 + 25) * ureg.kelvin
    gamma_mmt_sum = _compute_effective_junction_factor(f, alpha)

    kb = _KB_SI_J_PER_K * ureg.joule / ureg.kelvin
    return kb * temperature * xlink_concentration_0 * gamma_mmt_sum


def compute_trapping_factor(beta: float) -> float:
    """
    Compute the Langley trapping factor :math:`T_e`.

    Literature: :cite:t:`langley_elastically_1968`

    :param beta: :math:`P(F_b^{out})`, see :func:`~pylimer_tools.calc.miller_macosko_theory.compute_miller_macosko_probabilities()`
    :return: The Langley trapping factor
    """
    # for long B2s reacting with small A_fs
    return (1 - beta) ** 4
    # pel = ((1 / (p)) * (1 - alpha)) ** 2
    # return pel**2


def compute_probability_that_crosslink_is_effective(
    functionality_of_monomer: int, expected_degree_of_effect: int, p_f_a_out: float
) -> float:
    """
    Compute the probability that an Af, monomer will be an effective crosslink of exactly degree m

    :math:`P(X_m^f) = \\binom{f}{m} [P(F_A^{out})]^{f-m}[1-P(F_A^{out})]^m`

    Source:
        - Eq. 45 in :cite:t:`miller_new_1976a`

    :param functionality_of_monomer: f
    :param expected_degree_of_effect: m
    :param p_f_a_out: :math:`P(F_A^{out})`
    :return: The probability that a crosslink is effective
    """
    assert 0 <= p_f_a_out <= 1, "p_f_a_out must be between 0 and 1"
    f = functionality_of_monomer
    m = expected_degree_of_effect
    alpha = p_f_a_out
    return math.comb(f, m) * (alpha ** (f - m)) * ((1.0 - alpha) ** m)


def compute_probability_that_bifunctional_monomer_is_effective(
    p_f_b_out: float,
) -> float:
    """
    Consider a copolymerization of A_f with B_2.
    This function computes the probability that a random B_2 unit will be effective.

    :param p_f_b_out: :math:`P(F_B^{out})`
    :return: The probability that a bifunctional monomer is effective
    """
    assert 0 <= p_f_b_out <= 1, "p_f_b_out must be between 0 and 1"
    return (1 - p_f_b_out) ** 2


def compute_probability_that_crosslink_with_degree_is_dangling(
    functionality_of_monomer: int, degree_of_ineffectiveness: int, p_f_a_out: float
) -> float:
    """
    Consider a copolymerization of A_f with B_2.
    This function computes the probability that a random A_f unit will have i pendant arms.

    Source:
        - Eq. 6.3 in :cite:t:`miller_calculation_1979`

    :param functionality_of_monomer: f
    :param degree_of_ineffectiveness: i
    :param p_f_a_out: :math:`P(F_A^{out})`
    :return: The probability that a crosslink with degree is dangling
    """
    assert 0 <= p_f_a_out <= 1, "p_f_a_out must be between 0 and 1"
    f = functionality_of_monomer
    i = degree_of_ineffectiveness
    assert i <= f - 2, "degree_of_ineffectiveness must be less or equal to f-2"
    alpha = p_f_a_out
    # NOTE: verify that the last exponent is f - m, rather than f - 1 as in
    # the paper
    return math.comb(f, i) * (alpha ** (i)) * ((1.0 - alpha) ** (f - i))


def compute_probability_that_crosslink_is_dangling(
    functionality_of_monomer: int, p_f_a_out: float
) -> float:
    """
    Consider a copolymerization of A_f with B_2.
    This function computes the probability that a random A_f unit will be dangling (pendant).
    This is equal to the probability that only one of the arms is attached to the gel.

    Source:
        - Eq. 6.2 in :cite:t:`miller_calculation_1979`

    :param functionality_of_monomer: f
    :param p_f_a_out: :math:`P(F_A^{out})`
    :return: The probability that a crosslink is dangling
    """
    assert 0 <= p_f_a_out <= 1, "p_f_a_out must be between 0 and 1"
    f = functionality_of_monomer
    alpha = p_f_a_out
    return math.comb(f, 1) * (alpha ** (f - 1)) * (1.0 - alpha)


def compute_probability_that_bifunctional_monomer_is_dangling(
    p_f_b_out: float,
) -> float:
    """
    Consider a copolymerization of A_f with B_2.
    This function computes the probability that a random B_2 unit will be dangling.

    Source:
        - Eq. 6.1 in :cite:t:`miller_calculation_1979`

    :param p_f_b_out: :math:`P(F_B^{out})`
    :return: The probability that a bifunctional monomer is dangling
    """
    assert 0 <= p_f_b_out <= 1, "p_f_b_out must be between 0 and 1"
    return math.comb(2, 1) * (p_f_b_out) * (1.0 - p_f_b_out)


def predict_gelation_point(r: float, f: int, b2: float = 1.0) -> float:
    """
    Compute the gelation point :math:`p_{gel}` as theoretically predicted
    (gelation point = critical extent of reaction for gelation)

    Source:
      - :cite:t:`venkataraman_critical_1989`

    :param r: The stoichiometric imbalance of reactants (see: :func:`~pylimer_tools.calc.structure_analysis.compute_stoichiometric_imbalance`)
    :param f: functionality of the crosslinkers
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :return: p_gel: critical extent of reaction for gelation (``math.inf`` if gelation is impossible)
    """
    if r <= 0 or f <= 1 or b2 <= 0:
        return math.inf
    return math.sqrt(1 / (r * (f - 1) * b2))


def predict_maximum_p(r: float, f: int, b2: float = 1.0) -> Union[float, None]:
    """
    Compute the maximum crosslinker conversion possible given a stoichiometric imbalance.

    :param r: The stoichiometric imbalance of reactants (see: :func:`~pylimer_tools.calc.structure_analysis.compute_stoichiometric_imbalance`)
    :param f: functionality of the crosslinkers
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2. Since `r` already includes the number of active sites, this argument
        is not necessary.
    :return: p_max: The maximum crosslinker conversion possible
    """
    if r == 0 or f == 0:
        return None
    return min(1.0, 1.0 / r)


def predict_p_from_w_sol(
    w_sol: float,
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    weight_fractions: Union[dict, None] = None,
    r: Union[float, None] = None,
    b2: Union[float, None] = None,
) -> float:
    """
    Compute the extent of reaction based on the weight fraction of soluble material.

    :param w_sol: The weight fraction of soluble material
    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param functionality_per_type: a dictionary with key: type, and value: functionality of this atom type
    :param weight_fractions: a dictionary with the weight fraction of each type of atom
    :param r: The stoichiometric imbalance
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2
    :return: The extent of reaction p
    """
    param = _compute_validate_parameters(
        {**locals()},
        ["functionality_per_type", "weight_fractions", "r", "b2"],
    )

    functionality_per_type, weight_fractions, r, b2 = (
        param["functionality_per_type"],
        param["weight_fractions"],
        param["r"],
        param["b2"],
    )

    assert (
        isinstance(functionality_per_type, dict)
        and isinstance(weight_fractions, dict)
        and r is not None
        and b2 is not None
    )

    f_crosslinker = functionality_per_type[crosslinker_type]
    p_gel = predict_gelation_point(r=r, f=f_crosslinker, b2=b2)
    p_max = predict_maximum_p(r=r, f=f_crosslinker, b2=b2)
    if p_max is None or p_gel > p_max:
        raise ValueError(
            "Gelation is not possible for r = {}, f = {}, b2 = {} (p_gel = {}, p_max = {}).".format(
                r, f_crosslinker, b2, p_gel, p_max
            )
        )

    def compute_wsol(p_val: float) -> float:
        return compute_weight_fraction_of_soluble_material(
            network=None,
            crosslinker_type=crosslinker_type,
            functionality_per_type=functionality_per_type,
            weight_fractions=weight_fractions,
            r=r,
            p=p_val,
            b2=b2,
        )

    w_gel = compute_wsol(p_gel)
    w_max = compute_wsol(p_max)
    if math.isclose(w_sol, w_gel, abs_tol=1e-12):
        return p_gel
    if math.isclose(w_sol, w_max, abs_tol=1e-12):
        return p_max
    if w_sol < w_max or w_sol > w_gel:
        raise ValueError(
            "Target w_sol = {} is outside the achievable range [{}, {}] for p in [{}, {}].".format(
                w_sol, w_max, w_gel, p_gel, p_max
            )
        )

    try:
        from scipy import optimize
    except ImportError:
        raise ImportError(
            "scipy is required for predicting p from w_sol. "
            "Please install scipy to use this functionality."
        )

    res = optimize.root_scalar(
        lambda p_val: compute_wsol(p_val) - w_sol,
        bracket=(p_gel, p_max),
        method="brentq",
    )
    if not res.converged:
        warnings.warn("The p predicted from w_sol might be incorrect")
    return float(res.root)


def _validate_r_and_p(r: float, p: float, f: int) -> float:
    """
    Validate the parameters used in Miller-Macosko theory calculations.

    This function checks if the stoichiometric imbalance, crosslinker conversion,
    and crosslinker functionality are within valid ranges. It also verifies that
    the crosslinker conversion does not exceed the maximum possible value given
    the stoichiometric imbalance and functionality (allowing 1e-12 floating-point
    rounding tolerance and returning the clamped value).

    :param r: The stoichiometric imbalance (ratio of reactive groups)
    :param p: The crosslinker conversion (extent of reaction)
    :param f: The functionality of the crosslinker (number of reactive groups)
    :return: Validated (and tolerance-clamped) conversion ``p``
    :raises ValueError: If any parameter is outside its valid range or if p exceeds
                       the maximum possible conversion
    """
    if p < 0:
        raise ValueError(
            "The crosslinker conversion `p` must be positive, got {}".format(p)
        )
    if r < 0:
        raise ValueError(
            "The stoichiometric imbalance `r` must be positive, got {}".format(r)
        )
    if f < 2:
        raise ValueError(
            "The crosslinker functionality `f` must be >= 2, got {}".format(f)
        )
    p_max = predict_maximum_p(r=r, f=f)
    if p_max is None:
        raise ValueError(
            "For a system with r = {} and f = {}, maximum conversion could not be determined.".format(
                r, f
            )
        )
    if p > p_max:
        if math.isclose(p, p_max, rel_tol=1e-12, abs_tol=1e-12):
            return p_max
        raise ValueError(
            "For a system with r = {} and f = {}, p (in terms of crosslinkers) must be <= {}, {} given.".format(
                r, f, p_max, p
            )
        )
    return p


class _ParamValidatorAssembler:
    """
    A class to compute and validate one parameter.
    """

    def __init__(
        self,
        param_name: str,
        param_func: Callable,
        param_validator: Callable,
        dependencies: List[str],
    ):
        self.param_name = param_name
        self.param_func = param_func
        self.param_validator = param_validator
        self.dependencies = dependencies


"""
A list of parameter validators and assemblers used to compute and validate
parameters for Miller-Macosko theory calculations.

Each validator assembler contains:
- param_name: The name of the parameter to validate
- param_func: A function that computes the parameter value from other parameters
- param_validator: A function that validates the parameter value
- dependencies: A list of parameter names that the parameter depends on

This list is used by the _compute_validate_parameters function to automatically
compute missing parameters and validate existing ones
when performing Miller-Macosko theory calculations.
"""


def _infer_crosslinker_type(p: dict[str, Any]) -> int:
    """Infer crosslinker atom type from explicit functionality_per_type or network."""
    func_per_type = p.get("functionality_per_type")
    if isinstance(func_per_type, dict) and func_per_type:
        non_none = {k: v for k, v in func_per_type.items() if v is not None}
        if non_none:
            return max(non_none, key=lambda k: non_none[k])
    bead_func = p["network"].determine_functionality_per_type()
    if not bead_func:
        raise ValueError("Cannot determine crosslinker_type for an empty network.")
    return max(bead_func, key=lambda k: bead_func[k])


def _get_cached_strand_counts(p: dict[str, Any]) -> dict[int, int]:
    """Compute and cache strand functionality counts once per MMT evaluation."""
    cached = p.get("_strand_counts")
    if isinstance(cached, dict):
        return cached
    counts = count_strands_by_functionality(
        network=p["network"],
        crosslinker_type=p["crosslinker_type"],
        functionality_per_type=p["functionality_per_type"],
    )
    p["_strand_counts"] = counts
    return counts


def _compute_functionality_per_type(p: dict[str, Any]) -> dict[int, int]:
    """Infer complete functionality_per_type map from network and optional f."""
    inferred = determine_chain_functionality_per_type(
        p["network"], crosslinker_type=p["crosslinker_type"]
    )
    if p.get("f") is not None:
        inferred[p["crosslinker_type"]] = int(p["f"])
    return inferred


def _compute_r_from_cached(p: dict[str, Any]) -> float:
    """Compute stoichiometric imbalance sharing cached strand counts."""
    network: Universe = p["network"]
    if network.get_nr_of_atoms() == 0:
        return 0.0
    crosslinker_type: int = p["crosslinker_type"]
    atom_counts = network.count_atom_types()
    if crosslinker_type not in atom_counts:
        return 0.0
    strand_counts = _get_cached_strand_counts(p)
    other_formable = sum(func * cnt for func, cnt in strand_counts.items() if func > 0)
    if other_formable == 0:
        return math.inf
    crosslinker_formable = (
        atom_counts[crosslinker_type] * p["functionality_per_type"][crosslinker_type]
    )
    return crosslinker_formable / other_formable


def _compute_b2_from_cached(p: dict[str, Any]) -> float:
    """Compute b2 sharing cached strand counts."""
    if p["network"].get_nr_of_atoms() == 0:
        return 1.0
    strand_counts = _get_cached_strand_counts(p)
    n_mono = strand_counts.get(1, 0)
    n_bi = strand_counts.get(2, 0)
    total_sites = n_mono + 2 * n_bi
    if total_sites == 0:
        return 1.0
    return (2 * n_bi) / total_sites


def _compute_reactive_strand_number_density(p: dict[str, Any]) -> pint.Quantity:
    """Compute reactive strand number density ([B_2]_0 + [B_1]_0) / V, excluding solvent (B_0)."""
    network: Universe = p["network"]
    strand_counts = _get_cached_strand_counts(p)
    n_reactive_strands = sum(cnt for func, cnt in strand_counts.items() if func > 0)
    return n_reactive_strands / (
        network.get_volume() * p["unit_style"].get_base_unit_of("volume")
    )


_validators_assembler = [
    _ParamValidatorAssembler(
        "functionality_per_type",
        _compute_functionality_per_type,
        lambda x: isinstance(x, dict),
        ["network", "crosslinker_type"],
    ),
    _ParamValidatorAssembler(
        "weight_fractions",
        lambda p: p["network"].compute_weight_fractions(),
        lambda x: isinstance(x, dict),
        ["network"],
    ),
    _ParamValidatorAssembler(
        "crosslinker_type",
        _infer_crosslinker_type,
        lambda x: isinstance(x, int) and x >= 0,
        ["network"],
    ),
    _ParamValidatorAssembler(
        "f",
        lambda p: p["functionality_per_type"].get(p["crosslinker_type"], 0),
        lambda f: f >= 2 and math.isfinite(f),
        ["functionality_per_type", "crosslinker_type"],
    ),
    _ParamValidatorAssembler(
        "r",
        _compute_r_from_cached,
        lambda r: r > 0 and math.isfinite(r),
        ["network", "crosslinker_type", "functionality_per_type"],
    ),
    _ParamValidatorAssembler(
        "p",
        lambda p: compute_crosslinker_conversion(
            network=p["network"],
            crosslinker_type=p["crosslinker_type"],
            f=p.get("f"),
            functionality_per_type=p["functionality_per_type"],
        ),
        lambda p: 0 <= p <= 1,
        ["network", "crosslinker_type", "functionality_per_type"],
    ),
    _ParamValidatorAssembler(
        "nu",
        _compute_reactive_strand_number_density,
        lambda nu: nu.magnitude > 0 and math.isfinite(nu.magnitude),
        ["network", "crosslinker_type", "unit_style", "functionality_per_type"],
    ),
    _ParamValidatorAssembler(
        "b2",
        _compute_b2_from_cached,
        lambda b2: 0 <= b2 <= 1,
        ["network", "crosslinker_type", "functionality_per_type"],
    ),
    _ParamValidatorAssembler(
        "p_f_a_out",
        lambda p: compute_miller_macosko_probabilities(
            r=p["r"], p=p["p"], f=p["f"], b2=p["b2"]
        )[0],
        lambda alpha: 0 <= alpha <= 1,
        ["r", "p", "f", "b2"],
    ),
    _ParamValidatorAssembler(
        "p_f_b_out",
        lambda p: compute_miller_macosko_probabilities(
            r=p["r"], p=p["p"], f=p["f"], b2=p["b2"]
        )[1],
        lambda beta: 0 <= beta <= 1,
        ["r", "p", "f", "b2"],
    ),
    _ParamValidatorAssembler(
        "network",
        lambda p: p["network"],
        lambda x: isinstance(x, Universe),
        ["network"],
    ),
]
_validator_per_name = {v.param_name: v for v in _validators_assembler}


def _compute_validate_parameters(
    given_parameters: dict[str, Any], required_parameters: List[str]
) -> dict[str, Any]:
    """
    Compute missing parameters (e.g. from ``network``) and validate parameters.
    """
    given_parameters = dict(given_parameters)

    def _param_is_ready(param: str) -> bool:
        return param in given_parameters and given_parameters[param] is not None

    def _is_complete(dependencies: List[str]) -> bool:
        return all(_param_is_ready(param) for param in dependencies)

    def _can_be_computed(param: _ParamValidatorAssembler) -> bool:
        return all(_param_is_ready(dep) for dep in param.dependencies)

    def _validate(param_name: str) -> None:
        if not _validator_per_name[param_name].param_validator(
            given_parameters[param_name]
        ):
            raise ValueError(
                "Invalid value for parameter '{}' (got {}).".format(
                    param_name, given_parameters[param_name]
                )
            )

    # first, validate existing required parameters
    present = [d for d in required_parameters if _param_is_ready(d)]
    for p in present:
        _validate(p)

    # If crosslinker_type is missing, infer it from explicit functionality_per_type or network
    if not _param_is_ready("crosslinker_type"):
        if (
            _param_is_ready("functionality_per_type")
            and given_parameters["functionality_per_type"]
        ):
            given_parameters["crosslinker_type"] = _infer_crosslinker_type(
                given_parameters
            )
            _validate("crosslinker_type")
        elif _param_is_ready("network"):
            _validate("network")
            if given_parameters["network"].get_nr_of_atoms() > 0:
                given_parameters["crosslinker_type"] = _infer_crosslinker_type(
                    given_parameters
                )
                _validate("crosslinker_type")

    # If network and a (possibly partial) functionality_per_type were both provided,
    # complete missing atom types using chain functionality inference.
    if (
        _param_is_ready("network")
        and _param_is_ready("crosslinker_type")
        and _param_is_ready("functionality_per_type")
    ):
        _validate("network")
        _validate("crosslinker_type")
        _validate("functionality_per_type")
        resolved_func = _resolve_chain_functionality_per_type(
            given_parameters["network"],
            given_parameters["crosslinker_type"],
            given_parameters["functionality_per_type"],
        )
        if (
            _param_is_ready("f")
            and given_parameters["crosslinker_type"]
            not in given_parameters["functionality_per_type"]
        ):
            _validate("f")
            resolved_func[given_parameters["crosslinker_type"]] = int(
                given_parameters["f"]
            )
        given_parameters["functionality_per_type"] = resolved_func

    # Determine all parameters to compute and validate ready transitive dependencies
    to_compute = set(d for d in required_parameters if not _param_is_ready(d))
    validated_ready = set(present)
    found_last_iteration = True
    while found_last_iteration:
        found_last_iteration = False
        for p in list(to_compute):
            dependencies = _validator_per_name[p].dependencies
            for dep in dependencies:
                if not _param_is_ready(dep):
                    if dep not in to_compute:
                        to_compute.add(dep)
                        found_last_iteration = True
                elif dep in _validator_per_name and dep not in validated_ready:
                    _validate(dep)
                    validated_ready.add(dep)

    found_last_iteration = False
    while not _is_complete(required_parameters):
        # Instead of building a dependency graph,
        # we simply iterate over all validators and try to compute
        # the required parameters, as often as required.
        found_last_iteration = False

        for p in list(to_compute):
            if p not in to_compute:
                continue
            if _can_be_computed(_validator_per_name[p]):
                if p in ("p_f_a_out", "p_f_b_out"):
                    alpha, beta = compute_miller_macosko_probabilities(
                        r=given_parameters["r"],
                        p=given_parameters["p"],
                        f=given_parameters["f"],
                        b2=given_parameters["b2"],
                    )
                    if not _param_is_ready("p_f_a_out"):
                        given_parameters["p_f_a_out"] = alpha
                        to_compute.discard("p_f_a_out")
                        _validate("p_f_a_out")
                    if not _param_is_ready("p_f_b_out"):
                        given_parameters["p_f_b_out"] = beta
                        to_compute.discard("p_f_b_out")
                        _validate("p_f_b_out")
                else:
                    given_parameters[p] = _validator_per_name[p].param_func(
                        given_parameters
                    )
                    to_compute.remove(p)
                    _validate(p)
                found_last_iteration = True

        if not found_last_iteration:
            raise ValueError(
                "Missing required parameters: {}".format(", ".join(to_compute))
                + ".{}".format(
                    " Some of them may be computed as dependencies of others."
                    if len(to_compute) > 1
                    else ""
                )
            )

    return given_parameters
