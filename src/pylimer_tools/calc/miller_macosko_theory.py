import math
import warnings
from typing import Callable, List, Union

import numpy as np
import pint
import scipy.special
from scipy import optimize

from pylimer_tools.calc.structure_analysis import (
    compute_crosslinker_conversion,
    compute_fraction_of_bifunctional_reactive_sites,
    compute_stoichiometric_imbalance,
)
from pylimer_tools.io.unit_styles import UnitStyle
from pylimer_tools_cpp import Universe

"""
This module provides access to various computations introduced in the Miller-Macosko theory.

Caution:
    - not all systems are supported yet.
        In particular, for most methods, only A_f and B_2 is supported.
        Also, the systems are mostly assumed to be end-linked and monodisperse.

"""


def predict_shear_modulus(**kwargs):
    """
    Predict the shear modulus using MMT Analysis.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - see :func:`~pylimer_tools.calc.miller_macosko_theory.compute_modulus_decomposition`

    Returns:
      - G: the predicted shear modulus, or `None` if the universe is empty.

    ToDo:
      - Support more than one crosslinker type (as is supported by original formula)
    """
    g_mmt_phantom, g_mmt_entanglement, _, _ = compute_modulus_decomposition(**kwargs)
    return g_mmt_phantom + g_mmt_entanglement


def predict_number_density_of_junction_points(
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
) -> float:
    """
    Compute the number density of network strands using MMT

    Source:
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae)

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinker_type: the atom type to use to split the molecules
      - functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.

    Returns:
      - mu: The predicted number density of junction points
    """
    param = _compute_validate_parameters(
        {**locals()},
        ["functionality_per_type", "weight_fractions", "p_f_a_out"],
    )

    functionality_per_type, weight_fractions, alpha = (
        param["functionality_per_type"],
        param["weight_fractions"],
        param["p_f_a_out"],
    )
    if functionality_per_type[crosslinker_type] == 3:
        return weight_fractions[crosslinker_type] * (1 - alpha) ** 3
    elif functionality_per_type[crosslinker_type] == 4:
        return weight_fractions[crosslinker_type] * (
            4 * alpha * (1 - alpha) ** 3 + (1 - alpha) ** 4
        )
    else:
        raise NotImplementedError(
            "Currently, only cross-linker functionalities of 3 and 4 are supported, {} given.".format(
                functionality_per_type[crosslinker_type]
            )
        )


def predict_number_density_of_network_strands(
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    r: Union[float, None] = None,
    p: Union[float, None] = None,
) -> float:
    """
    Compute the number density of network strands using MMT

    Source:
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae)

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinker_type: the atom type to use to split the molecules
      - functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers

    Returns:
      - nu: The predicted number density of network strands
    """
    param = _compute_validate_parameters(
        {**locals()},
        ["functionality_per_type", "weight_fractions", "r", "p", "crosslinker_type"],
    )

    functionality_per_type, weight_fractions, r, p, crosslinker_type = (
        param["functionality_per_type"],
        param["weight_fractions"],
        param["r"],
        param["p"],
        param["crosslinker_type"],
    )

    if crosslinker_type not in weight_fractions:
        weight_fractions[crosslinker_type] = 0.0
    alpha, _ = compute_miller_macosko_probabilities(
        r=(
            r
            if r is not None
            else compute_stoichiometric_imbalance(
                network=network,
                crosslinker_type=crosslinker_type,
                functionality_per_type=functionality_per_type,
            )
        ),
        p=(
            p
            if p is not None
            else compute_crosslinker_conversion(
                network=network,
                crosslinker_type=crosslinker_type,
                f=functionality_per_type[crosslinker_type],
                functionality_per_type=functionality_per_type,
            )
        ),
        f=functionality_per_type[crosslinker_type],
    )

    if functionality_per_type[crosslinker_type] == 3:
        return (3 / 2) * weight_fractions[crosslinker_type] * (1 - alpha) ** 3
    elif functionality_per_type[crosslinker_type] == 4:
        return weight_fractions[crosslinker_type] * (
            6 * alpha * (1 - alpha) ** 3 + 2 * (1 - alpha) ** 4
        )
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
      - Eq. 6.4 in https://doi.org/10.1002/pen.760190409

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinker_type: the atom type to use to split the molecules
      - functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
      - weight_fractions: a dictionary with the weight fraction of each type of atom
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers
      - b2 (double, optional): the mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2

    Returns:
      - weightFraction :math:`$\\Phi_d = w_p$`: weightDangling/weightTotal
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
            # TODO: revise this, check if correct
            w_dangling += weight_fraction * (1.0 - beta)
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

    Arguments:
      - network: the polymer network to do the computation for
      - crosslinker_type: the type of the junctions/cross-linkers to select them in the network
      - functionality_per_type: a dictionary with key: atom type, and value: functionality atoms with this type.
      - weight_fractions: a dictionary with the weight fraction of each type of atom
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers
      - b2 (double, optional): the mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2

    Returns:
      - :math:`\\Phi_{el} = w_e`: weight fraction of network backbone
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
        else:
            w_elastic += (
                weight_fraction
                * compute_probability_that_bifunctional_monomer_is_effective(beta)
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
      - https://pubs.acs.org/doi/10.1021/ma00046a021
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737

    Arguments:
      - network: the polymer network to do the computation for
      - crosslinker_type: the type of the junctions/cross-linkers to select them in the network
      - weight_fractions (dict): a dictionary with key: type, and value: weight fraction of type.
            Pass if you want to omit the network.
      - functionality_per_type (dict): a dictionary with key: type, and value: functionality of this atom type.
          See: :func:`~pylimer_tools_cpp.Universe.determine_functionality_per_type`.
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers
      - b2 (double, optional): the mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2

    Returns:
      - :math:`W_{sol}` (float): the weight fraction of soluble material according to MMT.
      - weight_fractions (dict): a dictionary with key: type, and value: weight fraction of type
      - :math:`\\alpha` (float): Macosko & Miller's :math:`P(F_A)`
      - :math:`\\beta` (float): Macosko & Miller's :math:`P(F_B)`
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
):
    """
    Use MMT to compute the weight fraction of soluble material using

    :math:`W_{sol} = w_A_f P(F_A^{out})^f + w_B_g [rpP(F_A^{out})^{f-1}+1-rp]^g`

    Arguments:
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers
      - f: the functionality of the the crosslinker
      - w_f: the weight fraction of the crosslinkers
      - w_g: the weight fraction of ordinary chains
      - g: the functionality of the ordinary chains
    """
    alpha, _ = compute_miller_macosko_probabilities(r, p, f)
    return w_f * (alpha**f) + w_g * ((r * p * (alpha ** (f - 1)) + 1 - r * p) ** g)


def compute_miller_macosko_probabilities(r: float, p: float, f: int, b2: float = 1.0):
    """
    Compute Macosko and Miller's probabilities :math:`P(F_A)` and :math:`P(F_B)`
    i.e., the probability that a randomly chosen A (cross-link) or B (strand-end),
    respectively, is the start of a finite chain.

    Sources:
      - https://doi.org/10.1021/ma60050a003
      - https://doi.org/10.1021/ma60050a004
      - https://doi.org/10.1021/ma00046a021 (with monofunctional chains, f = 4)
      - https://doi.org/10.1021/cm0343507 (with monofunctional chains, f = 3)

    Note:
        Currently, only systems with B_2, B_1 and A_f are supported.

    Arguments:
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers
      - f: the functionality of the the crosslinker
      - b2: the fraction of bifunctional chains; defaults to 1.0 for no monofunctional chains.
            Can be computed e.g. as :math:`b_2 = \frac{2 \\cdot [B_2]}{[B_1] + 2 \\cdot [B_2]}`

    Returns:
      - alpha: :math:`P(F_A)`
      - beta: :math:`P(F_B)`
    """
    if r == 0 or p == 0 or f == 0:
        return 1.0, 1.0
    if b2 > 1 or b2 <= 0:
        raise ValueError("b2 must be in (0, 1), got b2 = {}".format(b2))

    # first, check a few things required by the formulae
    # since we want alpha, beta \in [0,1], given they are supposed to be
    # probabilities
    _validate_r_and_p(r, p, f)

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

        def fun_to_root_for_alpha(alpha):
            return r * b2 * p**2 * alpha ** (f - 1) - alpha - r * b2 * (p**2) + 1

        def fun_to_root_for_alpha_prime(alpha):
            return -1 + alpha ** (f - 2) * (-1 + f) * (p**2) * r * b2

        def fun_to_root_for_alpha_prime2(alpha):
            return alpha ** (f - 3) * (-2 + f) * (-1 + f) * (p**2) * r * b2

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
    return np.clip(alpha, 0, 1), np.clip(beta, 0, 1)


def compute_modulus_decomposition(
    network: Union[Universe, None] = None,
    ureg: pint.UnitRegistry = None,
    unit_style: Union[None, UnitStyle] = None,
    crosslinker_type: int = None,
    r: Union[float, None] = None,
    p: Union[float, None] = None,
    f: int = None,
    nu: Union[float, None] = None,
    temperature: pint.Quantity = None,
    functionality_per_type: Union[dict, None] = None,
    g_e_1: Union[float, None] = None,
    b2: float = 1.0,
):
    """
    Compute four different estimates of the plateau modulus, using MMT, ANM and PNM.

    Arguments:
      - network: the polymer network to do the computation for
      - unit_style: the unit style to use to have the results in appropriate units
      - crosslinker_type: the type of the junctions/cross-linkers to select them in the network
      - r: the stoichiometric imbalance. Optional if network is specified
      - p: the extent of reaction. Optional if network is specified
      - f: the functionality of the the crosslinker. Optional if network is specified
      - nu: the strand number density (nr of strands per volume) (ideally with units). Optional if network is specified
      - temperature: the temperature to compute the modulus at. Default: 298.15 K
          Optional, can be passed to improve performance
      - functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
          Optional, can be passed to improve performance
      - g_e_1: the melt entanglement modulus
      - b2 (double, optional): the mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2

    Returns:
      - G_MMT_phantom: the phantom contribution to the MMT modulus;
          see also :func:`pylimer_tools.calc.miller_macosko_theory.computeJunctionModulus`
      - G_MMT_entanglement: the entanglement contribution to the MMT modulus
      - g_anm: the ANM estimate of the modulus
      - g_pnm: the PNM estimate of the modulus
    """
    if ureg is None:
        if unit_style is None:
            raise ValueError(
                "Unit style or unit registry must be specified to compute modulus."
            )
        ureg = unit_style.get_underlying_unit_registry()

    param = _compute_validate_parameters(
        {**locals()},
        ["f", "nu", "p_f_a_out", "p_f_b_out", "p", "r"],
    )

    f, nu, alpha, beta, p, r = (
        param["f"],
        param["nu"],
        param["p_f_a_out"],
        param["p_f_b_out"],
        param["p"],
        param["r"],
    )

    if temperature is None:
        temperature = (273.15 + 25) * ureg.kelvin  # Temperature in Kelvin
    if g_e_1 is None:
        g_e_1 = (
            8.3145  # gas constant, J/(mol*K)
            * temperature.to("kelvin").magnitude  # Temperature in Kelvin
            * 1e-6
            * 94.79281
        ) * ureg("MPa")
        # -> MPa, melt entanglement modulus of PDMS

    # Boltzmann constant
    kb = 1.380649e-23 * ureg.joule / ureg.kelvin
    # affine
    g_anm = nu * kb * temperature
    # phantom
    g_pnm = (1 - 2 / f) * nu * kb * temperature if f != 0 else 0.0
    # MMT:
    gamma_mmt_sum = 0.0
    for m in range(3, f + 1):
        gamma_mmt_sum += (
            (m - 2) / 2
        ) * compute_probability_that_crosslink_is_effective(f, m, alpha)
    gamma_mmt = (2 * r * b2 / f) * gamma_mmt_sum if f != 0 else 0.0
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
    alpha: Union[float, None] = None,
    temperature: pint.Quantity = None,
    ureg: pint.UnitRegistry = None,
):
    """
    Compute MMT's modulus, assuming the solvent is removed

    Arguments:
        - p: the cross-linker conversion
        - r: the stoichiometric imbalance
        - f: the functionality of the crosslinkers
        - g_e_1: the melt entanglement modulus :math:`G_e(1) = k_B T \\epsilon_e`
        - xlink_concentration_0: [A_f]_0, in 1/volume units
        - alpha: :math:`P(F_a^{out})`, optional
        - temperature: the temperatures; defaults to room temperature
        - w_sol: the soluble fraction (to be removed)
    """
    if temperature is None:
        if ureg is None:
            raise ValueError(
                "Unit registry must be initialized, or temperature specified."
            )
        temperature = (273.15 + 25) * ureg.kelvin
    if alpha is None:
        alpha, _ = compute_miller_macosko_probabilities(r, p, f)

    junction_part = (1 - w_sol) ** (-1 / 3) * compute_junction_modulus(
        p=p,
        r=r,
        xlink_concentration_0=xlink_concentration_0,
        f=f,
        alpha=alpha,
        ureg=ureg,
        temperature=temperature,
    )
    entanglement_part = (1 - w_sol) ** (-2) * compute_entanglement_modulus(
        p=p,
        r=r,
        f=f,
        g_e_1=g_e_1,
        alpha=alpha,
        ureg=ureg,
        temperature=temperature,
    )
    return junction_part + entanglement_part


def compute_entanglement_modulus(
    g_e_1: pint.Quantity,
    temperature: pint.Quantity,
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    p: Union[float, None] = None,
    r: Union[float, None] = None,
    f: Union[int, None] = None,
    b2: Union[float, None] = None,
    beta: Union[float, None] = None,
):
    """
    Compute MMT's entanglement contribution to the equilibrium shear modulus, given by
    :math:`k_B T \\epsilon_e T_e`.

    Arguments:
      - g_e_1: the melt entanglement modulus :math:`G_e(1) = k_B T \\epsilon_e`
      - temperature: the temperatures; defaults to room temperature (25 °C)
      - network: the polymer network to do the computation for
      - crosslinker_type: the type of the junctions/cross-linkers to select them in the network
      - p: the cross-linker conversion
      - r: the stoichiometric imbalance
      - f: the functionality of the crosslinkers
      - beta: :math:`P(F_b^{out})`, optional
    """
    param = _compute_validate_parameters(
        {**locals()},
        ["p_f_b_out"],
    )

    return compute_trapping_factor(beta=param["p_f_b_out"]) * g_e_1


def compute_junction_modulus(
    p: float,
    r: float,
    xlink_concentration_0: pint.Quantity,
    ureg: pint.UnitRegistry,
    f: Union[int, None] = None,
    alpha: Union[float, None] = None,
    temperature: pint.Quantity = None,
):
    """
    Compute MMT's junction modulus, given by
    :math:`G_{junctions} = k_B T [A_f]_0 \\sum_{m=3}^{f} \frac{m-2}{2} P(X_{m,f})`.

    Arguments:
        - p: the cross-linker conversion
        - r: the stoichiometric imbalance
        - xlink_concentration_0: [A_f]_0, in 1/volume units
        - f: the functionality of the crosslinkers
        - alpha: :math:`P(F_a^{out})`, optional
        - temperature: the temperatures; defaults to room temperature (25 °C)
    """
    if temperature is None:
        temperature = (273.15 + 25) * ureg.kelvin
    if alpha is None:
        alpha, _ = compute_miller_macosko_probabilities(r, p, f)
    gamma_mmt_sum = 0.0
    for m in range(3, f + 1):
        gamma_mmt_sum += (
            (m - 2) / 2
        ) * compute_probability_that_crosslink_is_effective(f, m, alpha)

    kb = 1.380649e-23 * ureg.joule / ureg.kelvin
    return kb * temperature * xlink_concentration_0 * gamma_mmt_sum


def compute_trapping_factor(beta: float) -> float:
    """
    Compute the Langley trapping factor :math:`T_e`.

    Literature: https://doi.org/10.1021/ma60004a015

    Arguments:
        - beta: :math:`P(F_b^{out})`, see :func:`~pylimer_tools.calc.miller_macosko_theory.compute_mms_probabilities()`
        - p: the extent of reaction in terms of the crosslinkers.
    """
    # for long B2s reacting with small A_fs
    return (1 - beta) ** 4
    # pel = ((1 / (p)) * (1 - alpha)) ** 2
    # return pel**2


def compute_probability_that_crosslink_is_effective(
    functionality_of_monomer: int, expected_degree_of_effect: int, p_f_a_out: float
):
    """
    Compute the probability that an Af, monomer will be an effective cross-link of exactly degree m

    :math:`P(X_m^f) = \binom{f}{m} [P(F_A^{out})]^{f-m}[1-P(F_A^{out})]^m`

    Source:
        - Eq. 45 in Miller, Macosko 1976, A New Derivation of Post Gel Properties of Network

    Arguments:
        - functionality_of_monomer: f
        - expected_decree_of_effect: m
        - alpha: :math:`P(F_A^{out})`
    """
    assert 0 <= p_f_a_out <= 1, "p_f_a_out must be between 0 and 1"
    f = functionality_of_monomer
    m = expected_degree_of_effect
    alpha = p_f_a_out
    return scipy.special.binom(f, m) * (alpha ** (f - m)) * ((1.0 - alpha) ** m)


def compute_probability_that_bifunctional_monomer_is_effective(p_f_b_out: float):
    """
    Consider a copolymerization of A_f with B_2.
    This function computes the probability that a random B_2 unit will be effective.

    Arguments:
        - beta: :math:`P(F_B^{out})`
    """
    assert 0 <= p_f_b_out <= 1, "p_f_b_out must be between 0 and 1"
    return (1 - p_f_b_out) ** 2


def compute_probability_that_crosslink_with_degree_is_dangling(
    functionality_of_monomer: int, degree_of_ineffectiveness: int, p_f_a_out: float
):
    """
    Consider a copolymerization of A_f with B_2.
    This function computes the probability that a random A_f unit will have i pendant arms.

    Source:
        - Eq. 6.3 in https://doi.org/10.1002/pen.760190409

    Arguments:
        - functionality_of_monomer: f
        - degree_of_ineffectiveness: i
        - alpha: :math:`P(F_A^{out})`
    """
    assert 0 <= p_f_a_out <= 1, "p_f_a_out must be between 0 and 1"
    f = functionality_of_monomer
    i = degree_of_ineffectiveness
    assert i <= f - 2, "degree_of_ineffectiveness must be less or equal to f-2"
    alpha = p_f_a_out
    # NOTE: verify that the last exponent is f - m, rather than f - 1 as in
    # the paper
    return scipy.special.binom(f, i) * (alpha ** (i)) * ((1.0 - alpha) ** (f - i))


def compute_probability_that_crosslink_is_dangling(
    functionality_of_monomer: int, p_f_a_out: float
):
    """
    Consider a copolymerization of A_f with B_2.
    This function computes the probability that a random A_f unit will be dangling (pendant).
    This is equal to the probability that only one of the arms is attached to the gel.

    Source:
        - Eq. 6.2 in https://doi.org/10.1002/pen.760190409

    Arguments:
        - functionality_of_monomer: f
        - alpha: :math:`P(F_A^{out})`
    """
    assert 0 <= p_f_a_out <= 1, "p_f_a_out must be between 0 and 1"
    f = functionality_of_monomer
    alpha = p_f_a_out
    return scipy.special.binom(f, 1) * (alpha ** (f - 1)) * (1.0 - alpha)


def compute_probability_that_bifunctional_monomer_is_dangling(p_f_b_out: float):
    """
    Consider a copolymerization of A_f with B_2.
    This function computes the probability that a random B_2 unit will be dangling.

    Source:
        - Eq. 6.1 in https://doi.org/10.1002/pen.760190409

    Arguments:
        - beta: :math:`P(F_B^{out})`
    """
    assert 0 <= p_f_b_out <= 1, "p_f_b_out must be between 0 and 1"
    return scipy.special.binom(2, 1) * (p_f_b_out) * ((1.0 - p_f_b_out))


def predict_gelation_point(r: float, f: int, b2: int = 1) -> float:
    """
    Compute the gelation point :math:`p_{gel}` as theoretically predicted
    (gelation point = critical extent of reaction for gelation)

    Source:
      - https://www.sciencedirect.com/science/article/pii/003238618990253X

    Arguments:
      - r (double): the stoichiometric imbalance of reactants (see: #compute_stoichiometric_imbalance)
      - f (int): functionality of the crosslinkers
      - b2 (double, optional): the mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2

    Returns:
      - p_gel: critical extent of reaction for gelation
    """
    # if (r is None):
    #   r = calculateEffectiveCrosslinkerFunctionality(network, crosslinker_type, f)
    return math.sqrt(1 / (r * (f - 1) * b2))


def predict_maximum_p(r: float, f: int, b2: float = 1) -> float:
    """
    Compute the maximum cross-linker conversion possible given a stoichiometric inbalance.

    Arguments:
      - r (double): the stoichiometric imbalance of reactants (see: #compute_stoichiometric_imbalance)
      - f (int): functionality of the crosslinkers
      - b2 (double, optional): the mole fraction of reactive sites in B2 among all reactive sites
        in a mixture of B1 and B2

    Returns:
      - p_max: the maximum cross-linker conversion possible
    """
    # assume:
    n_chains = 1000
    # -> compute:
    n_xlinks = r * 2 * b2 * n_chains / f
    if n_xlinks == 0:
        return None
    max_possible_bonds = min(2 * n_chains, f * n_xlinks)
    p_max = max_possible_bonds / (n_xlinks * f)
    return p_max


def predict_p_from_w_sol(
    w_sol: float,
    network: Union[Universe, None] = None,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    weight_fractions: Union[dict, None] = None,
    r: Union[float, None] = None,
    b2: Union[float, None] = None,
):
    """
    Compute the extent of reaction based on the weight fraction of soluble material.

    Arguments:
      - w_sol: the weight fraction of soluble material
      - r: the stoichiometric imbalance
      - w_f: the weight fraction of crosslinkers with functionality f,
      - w_g: the weight fraction of precursor chains with functionality g,
      - f: the functionality of the crosslinkers
      - g: the functionality of the precursor chains
    """
    def compute_wsol(p):
        return compute_weight_fraction_of_soluble_material(
            network=network,
            crosslinker_type=crosslinker_type,
            functionality_per_type=functionality_per_type,
            weight_fractions=weight_fractions,
            r=r,
            p=p,
            b2=b2,
        )

    res = optimize.minimize_scalar(
        lambda p: abs(w_sol - compute_wsol(p)), bounds=[1e-5, 1.0 - 1e-5]
    )
    if not res.success:
        warnings.warn("The p predicted from w_sol might be incorrect")
    return res.x


def _validate_r_and_p(r: float, p: float, f: int):
    if p < 0:
        raise ValueError(
            "The cross-linker conversion `p` must be positive, got {}".format(p)
        )
    if r < 0:
        raise ValueError(
            "The stoichiometric imbalance `r` must be positive, got {}".format(r)
        )
    if f < 2:
        raise ValueError(
            "The cross-linker functionality `f` must be >= 2, got {}".format(f)
        )
    # assume:
    p_max = predict_maximum_p(r=r, f=f)
    if p > p_max:
        raise ValueError(
            "For a system with r = {} and f = {}, p (in terms of crosslinkers) must be < {}, {} given.".format(
                r, f, p_max, p
            )
        )


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


_validators_assembler = [
    _ParamValidatorAssembler(
        "functionality_per_type",
        lambda p: p["network"].determine_functionality_per_type(),
        lambda x: isinstance(x, dict),
        ["network"],
    ),
    _ParamValidatorAssembler(
        "weight_fractions",
        lambda p: p["network"].compute_weight_fractions(),
        lambda x: isinstance(x, dict),
        ["network"],
    ),
    _ParamValidatorAssembler(
        "crosslinker_type",
        lambda p: max(p["functionality_per_type"], key=p["functionality_per_type"].get),
        lambda x: isinstance(x, int) and x >= 0,
        ["functionality_per_type"],
    ),
    _ParamValidatorAssembler(
        "f",
        lambda p: p["functionality_per_type"].get(p["crosslinker_type"], 0),
        lambda f: f >= 2 and np.isfinite(f),
        ["functionality_per_type", "crosslinker_type"],
    ),
    _ParamValidatorAssembler(
        "r",
        lambda p: compute_stoichiometric_imbalance(
            network=p["network"],
            crosslinker_type=p["crosslinker_type"],
            functionality_per_type=p["functionality_per_type"],
        ),
        lambda r: r > 0 and np.isfinite(r),
        ["network", "crosslinker_type"],
    ),
    _ParamValidatorAssembler(
        "p",
        lambda p: compute_crosslinker_conversion(
            network=p["network"],
            crosslinker_type=p["crosslinker_type"],
            functionality_per_type=p["functionality_per_type"],
        ),
        lambda p: 0 <= p <= 1,
        ["network", "crosslinker_type"],
    ),
    _ParamValidatorAssembler(
        "nu",
        lambda p: (
            len(p["network"].get_molecules(p["crosslinker_type"]))
            / (p["network"].get_volume() * p["unit_style"].get_base_unit_of("volume"))
        ),
        lambda nu: nu.magnitude > 0 and np.isfinite(nu.magnitude),
        ["network", "crosslinker_type", "unit_style"],
    ),
    _ParamValidatorAssembler(
        "b2",
        lambda p: compute_fraction_of_bifunctional_reactive_sites(
            network=p["network"],
            crosslinker_type=p["crosslinker_type"],
            functionality_per_type=p["functionality_per_type"],
        ),
        lambda b2: 0 <= b2 <= 1,
        ["network", "crosslinker_type"],
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
    given_parameters: dict, required_parameters: List[str]
) -> dict:
    """
    Slightly overengineered function
    to compute missing parameters e.g. from the network,
    and validate parameters.
    """

    def _param_is_ready(param: str) -> bool:
        return param in given_parameters and not given_parameters[param] is None

    def _is_complete(dependencies: List[str]) -> bool:
        return all(_param_is_ready(param) for param in dependencies)

    def _can_be_computed(param: _ParamValidatorAssembler) -> bool:
        return all(_param_is_ready(dep) for dep in param.dependencies)

    def _validate(param_name: str):
        if not _validator_per_name[param_name].param_validator(given_parameters[p]):
            raise ValueError(
                "Invalid value for parameter '{}' (got {}).".format(
                    param_name, given_parameters[param_name]
                )
            )

    # first, validate existing parameters
    present = [d for d in required_parameters if _param_is_ready(d)]
    for p in present:
        _validate(p)

    # first, determine all parameters to compute
    to_compute = set([d for d in required_parameters if not _param_is_ready(d)])
    # add dependencies
    found_last_iteration = True
    while found_last_iteration:
        found_last_iteration = False
        for p in list(to_compute):
            dependencies = _validator_per_name[p].dependencies
            for dep in dependencies:
                if dep not in to_compute and not _param_is_ready(dep):
                    to_compute.add(dep)
                    found_last_iteration = True

    found_last_iteration = False
    while not _is_complete(required_parameters):
        # Instead of building a dependency graph,
        # we simply iterate over all validators and try to compute
        # the required parameters, as often as required.
        found_last_iteration = False

        for p in list(to_compute):
            if _can_be_computed(_validator_per_name[p]):
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
