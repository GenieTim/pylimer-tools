import math
import warnings
from typing import Union

import numpy as np
import pint
import scipy.special
from scipy import optimize

from pylimer_tools.calc.structure_analysis import (
    compute_crosslinker_conversion,
    compute_effective_crosslinker_functionality,
    compute_stoichiometric_imbalance, compute_weight_fractions,
    measure_weight_fraction_of_soluble_material)
from pylimer_tools.io.unit_styles import UnitStyle
from pylimer_tools_cpp import Universe

"""
This module provides access to various computations introduced in the Miller-Macosko theory.

Caution:
    - not all systems are supported yet.
        In particular, for most methods, only A_f and B_2 is supported.

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
    g_mmt_phantom, g_mmt_entanglement, _, _ = compute_modulus_decomposition(
        **kwargs)
    return g_mmt_phantom + g_mmt_entanglement


def predict_number_density_of_junction_points(
    network: Universe, crosslinker_type: int, functionality_per_type: dict = None
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
    if functionality_per_type is None or crosslinker_type not in functionality_per_type:
        functionality_per_type = network.determine_functionality_per_type()

    weight_fractions, alpha, _ = compute_weight_fractions_and_probabilities(
        network, crosslinker_type, functionality_per_type
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
    network: Universe,
    crosslinker_type: int = 2,
    functionality_per_type: dict = None,
    r: float = None,
    p: float = None,
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
    if (functionality_per_type is None) or (
        crosslinker_type not in functionality_per_type
    ):
        functionality_per_type = network.determine_functionality_per_type()

    if crosslinker_type not in functionality_per_type:
        raise ValueError("Could not determine cross-linker functionality")

    weight_fractions = network.compute_weight_fractions()
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
    network: Universe,
    crosslinker_type: int,
    functionality_per_type: dict = None,
    weight_fractions: dict = None,
    r: float = None,
    p: float = None,
) -> float:
    """
    Compute the weight fraction of dangling strands in infinite network

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinker_type: the atom type to use to split the molecules
      - functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
      - weight_fractions: a dictionary with the weight fraction of each type of atom
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers

    Returns:
      - weightFraction $\\Phi_d = 1 - \\Phi_{el} - w_{sol}$: weightDangling/weightTotal
    """

    # possible alternative?!:
    # 2*beta*(1-beta)

    return (
        1.0
        - compute_weight_fraction_of_backbone(
            network, crosslinker_type, functionality_per_type, weight_fractions, r, p
        )
        - compute_weight_fraction_of_soluble_material(
            network, crosslinker_type, functionality_per_type, weight_fractions, r, p
        )
    )


def compute_weight_fraction_of_backbone(
    network: Universe,
    crosslinker_type: int,
    functionality_per_type: dict = None,
    weight_fractions: dict = None,
    r: float = None,
    p: float = None,
) -> float:
    """
    Compute the weight fraction of the backbone strands in an infinite network

    Source:
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae)
      - https://pubs.acs.org/doi/10.1021/ma00046a021 (see appendix)

    Arguments:
      - network: the polymer network to do the computation for
      - crosslinker_type: the type of the junctions/cross-linkers to select them in the network
      - functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
      - weight_fractions: a dictionary with the weight fraction of each type of atom
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers

    Returns:
      - :math:`\\Phi_{el}`: weight fraction of network backbone
    """
    if network is not None and network.get_nr_of_atoms() == 0:
        return 0

    if functionality_per_type is None or crosslinker_type not in functionality_per_type:
        functionality_per_type = network.determine_functionality_per_type()

    weight_fractions, alpha, beta = compute_weight_fractions_and_probabilities(
        network, crosslinker_type, functionality_per_type, weight_fractions, r, p
    )
    w_sol = compute_weight_fraction_of_soluble_material(
        network, crosslinker_type, functionality_per_type, weight_fractions, r, p
    )
    if w_sol < 0 or w_sol > 1:
        warnings.warn(
            "The weight fraction w_sol predicted by MMT ({}) is outside accepted range. ".format(
                w_sol
            )
            + "Falling back to measurement."
        )
        w_sol = measure_weight_fraction_of_soluble_material(network)

    phi_el = 0
    w_a = weight_fractions[crosslinker_type] / \
        functionality_per_type[crosslinker_type]
    w_xl = weight_fractions[crosslinker_type]
    w_x2 = 1 - w_xl
    assert w_a <= 1 and w_a >= 0
    assert w_xl <= 1 and w_xl >= 0
    assert w_x2 <= 1 and w_x2 >= 0
    assert w_sol <= 1 and w_sol >= 0
    if functionality_per_type[crosslinker_type] == 3:
        phi_el = (
            (w_x2 * (1 - beta) ** 2)
            + (w_xl * ((1 - alpha) ** 3 + 3 * alpha * (1 - w_a) * ((1 - alpha) ** 2)))
        ) / (1 - w_sol)
    else:
        assert functionality_per_type[crosslinker_type] == 4
        phi_el = (
            (w_x2 * (1 - beta) ** 2)
            + (
                w_xl
                * (
                    ((1 - alpha) ** 4)
                    + 4 * alpha * (1 - w_a) * ((1 - alpha) ** 3)
                    + 6 * (alpha**2) * (1 - 2 * w_a) * (1 - alpha) ** 2
                )
            )
        ) / (1 - w_sol)

    return phi_el


def compute_weight_fraction_of_soluble_material(
    network: Universe,
    crosslinker_type: int,
    functionality_per_type: dict = None,
    weight_fractions: dict = None,
    r: float = None,
    p: float = None,
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

    Returns:
      - :math:`W_{sol}` (float): the weight fraction of soluble material according to MMT.
      - weight_fractions (dict): a dictionary with key: type, and value: weight fraction of type
      - :math:`\\alpha` (float): Macosko & Miller's :math:`P(F_A)`
      - :math:`\\beta` (float): Macosko & Miller's :math:`P(F_B)`
    """
    if network is not None and network.get_nr_of_bonds() == 0:
        return 1.0

    weight_fractions, alpha, beta = compute_weight_fractions_and_probabilities(
        network, crosslinker_type, functionality_per_type, weight_fractions, r, p
    )

    if functionality_per_type is not None and not np.all(
        [
            key in functionality_per_type
            for key in weight_fractions.keys()
            if weight_fractions[key] > 0.0
        ]
    ):
        warnings.warn(
            "functionality_per_type does not contain functionality for all types. Will be re-computed."
        )
        functionality_per_type = None

    if functionality_per_type is None:
        _require_network(network, "functionality_per_type")
        functionality_per_type = network.determine_functionality_per_type()

    w_sol = 0
    for key in weight_fractions:
        coefficient = alpha if key == crosslinker_type else beta
        if key not in weight_fractions or math.isclose(
            weight_fractions[key], 0, abs_tol=1e-10
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


def compute_weight_fractions_and_probabilities(
    network: Universe,
    crosslinker_type: int,
    functionality_per_type: dict = None,
    weight_fractions: dict = None,
    r: float = None,
    p: float = None,
):
    """
    Shortcut function filling all missing parameters,
    computing the weight fractions per type in the network,
    and the MMT probabilities :math:`P(F_a^{out})` and  :math:`P(F_b^{out})`

    Arguments:
      - network: the polymer network to do the computation for
      - crosslinker_type: the type of the junctions/cross-linkers to select them in the network
      - weight_fractions (dict): a dictionary with key: type, and value: weight fraction of type.
          Pass if you want to omit the network.
      - functionality_per_type (dict): a dictionary with key: type, and value: functionality of this atom type.
          See: :func:`~pylimer_tools_cpp.Universe.determine_functionality_per_type`.
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers
    """
    if functionality_per_type is None or crosslinker_type not in functionality_per_type:
        assert network is not None
        functionality_per_type = network.determine_functionality_per_type()
        if crosslinker_type not in functionality_per_type:
            raise ValueError(
                "The crosslinker type {} is not present in the network. Got types {}".format(
                    crosslinker_type,
                    ", ".join([str(t) for t in functionality_per_type.keys()]),
                )
            )

    if functionality_per_type[crosslinker_type] not in range(3, 7):
        raise NotImplementedError(
            "Currently, a crosslinker functionality of {} is not supported.".format(
                functionality_per_type[crosslinker_type]
            )
        )

    if weight_fractions is None:
        weight_fractions = compute_weight_fractions(network)
        assert math.isclose(
            sum(w for w in weight_fractions.values()), 1.0, abs_tol=1e-9
        )
        if crosslinker_type not in weight_fractions:
            weight_fractions[crosslinker_type] = 0.0

    for key in functionality_per_type:
        if (
            key != crosslinker_type
            and functionality_per_type[key] != 2
            and weight_fractions[key] > 0
        ):
            raise NotImplementedError(
                "Currently, only strand functionality of 2 is supported. {} given for type {}".format(
                    functionality_per_type[key], key
                )
            )

    if p is None:
        assert network is not None
        p = compute_crosslinker_conversion(
            network, crosslinker_type, functionality_per_type[crosslinker_type]
        )
        if p < 0 or p > 1:
            warnings.warn(
                "The p computed ({}) is outside the accepted range. ".format(p)
                + "Falling back to effective cross-linker functionality."
            )
            p = compute_effective_crosslinker_functionality(
                network, crosslinker_type)
    if p > 1 or p < 0:
        raise ValueError(
            "Detected p = {} for f = {}. Need p in (0, 1).".format(
                p, functionality_per_type[crosslinker_type]
            )
        )
    if r is None:
        assert network is not None
        r = compute_stoichiometric_imbalance(
            network, crosslinker_type, functionality_per_type=functionality_per_type
        )
    assert r >= 0

    alpha, beta = compute_miller_macosko_probabilities(
        r, p, functionality_per_type[crosslinker_type]
    )
    assert alpha <= 1 and alpha >= 0
    assert beta <= 1 and beta >= 0

    return weight_fractions, alpha, beta


def compute_miller_macosko_probabilities(r: float, p: float, f: int):
    """
    Compute Macosko and Miller's probabilities :math:`P(F_A)` and :math:`P(F_B)`
    i.e., the probability that a randomly chosen A (cross-link) or B (strand-end),
    respectively, is the start of a finite chain.

    Sources:
      - https://pubs.acs.org/doi/10.1021/ma60050a004
      - https://doi.org/10.1021/ma60050a003

    Note:
        Currently, only systems with B_2 and A_f are supported.

    Arguments:
      - r: the stoichiometric imbalance
      - p: the extent of reaction in terms of the crosslinkers
      - f: the functionality of the the crosslinker

    Returns:
      - alpha: :math:`P(F_A)`
      - beta: :math:`P(F_B)`
    """
    if r == 0 or p == 0 or f == 0:
        return 1.0, 1.0

    # first, check a few things required by the formulae
    # since we want alpha, beta \in [0,1], given they are supposed to be
    # probabilities
    validate_r_and_p(r, p, f)

    # actually do the calculations
    if f == 3:
        alpha = (1 - r * p * p) / (r * p * p)
        if not (1 / (p**2) < 2 * r and (1 / (p**2) > r)):
            warnings.warn(
                "The resulting P(F_A) is probably unreliable, "
                + "as the detected root does not fulfill the required conditions."
            )
    elif f == 4:
        alpha = ((1.0 / (r * p * p)) - 3.0 / 4.0) ** (1.0 / 2.0) - (1.0 / 2.0)
        if not (1 / (p**2) < 3 * r and (1 / (p**2) > r)):
            warnings.warn(
                "The resulting P(F_A) is probably unreliable, "
                + "as the detected root does not fulfill the required conditions."
            )
    else:

        def fun_to_root_for_alpha(alpha):
            return r * p**2 * alpha ** (f - 1) - alpha - r * (p**2) + 1

        def fun_to_root_for_alpha_prime(alpha):
            return -1 + alpha ** (f - 2) * (-1 + f) * (p**2) * r

        def fun_to_root_for_alpha_prime2(alpha):
            return alpha ** (f - 3) * (-2 + f) * (-1 + f) * (p**2) * r

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
            "The resulting P(F_A) from r = {}, p = {} for f = {} is probably unreliable, ".format(
                r, p, f
            )
            + "as it will be clipped to [0,1] from {}".format(alpha)
        )
    if beta > 1 or beta < 0:
        warnings.warn(
            "The resulting P(F_B) from r = {}, p = {} for f = {} is probably unreliable, ".format(
                r, p, f
            )
            + "as it will be clipped to [0,1] from {}".format(beta)
        )
    return np.clip(alpha, 0, 1), np.clip(beta, 0, 1)  # TODO: reconsider


def validate_r_and_p(r: float, p: float, f: int):
    if p < 0:
        raise ValueError(
            "The cross-linker conversion `p` must be positive, got {}".format(
                p)
        )
    if r < 0:
        raise ValueError(
            "The stoichiometric imbalance `r` must be positive, got {}".format(
                r)
        )
    if f < 2:
        raise ValueError(
            "The cross-linker functionality `f` must be >= 2, got {}".format(f)
        )
    # assume:
    n_chains = 1000
    # -> compute:
    n_xlinks = r * 2 * n_chains / f
    max_possible_bonds = min(2 * n_chains, f * n_xlinks)
    if n_xlinks == 0:
        return
    p_max = max_possible_bonds / (n_xlinks * f)
    if p > p_max:
        raise ValueError(
            "For a system with r = {} and f = {}, p (in terms of crosslinkers) must be < {}, {} given.".format(
                r, f, p_max, p
            )
        )


def compute_modulus_decomposition(
    network: Universe,
    unit_style: UnitStyle,
    crosslinker_type: int = None,
    r: float = None,
    p: float = None,
    f: int = None,
    nu: float = None,
    temperature: pint.Quantity = None,
    functionality_per_type: dict = None,
    g_e_1: float = None,
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

    Returns:
      - G_MMT_phantom: the phantom contribution to the MMT modulus;
          see also :func:`pylimer_tools.calc.miller_macosko_theory.computeJunctionModulus`
      - G_MMT_entanglement: the entanglement contribution to the MMT modulus
      - g_anm: the ANM estimate of the modulus
      - g_pnm: the PNM estimate of the modulus
    """
    if (crosslinker_type is None or network is None) and (
        r is None or f is None or p is None or nu is None
    ):
        raise ValueError(
            "Either the network and crosslinker_type or the required variables must be specified"
        )
    if r is None:
        r = compute_stoichiometric_imbalance(
            network,
            crosslinker_type=crosslinker_type,
            functionality_per_type=functionality_per_type,
        )
    if p is None:
        p = compute_crosslinker_conversion(
            network, crosslinker_type, functionality_per_type=functionality_per_type
        )
        if p < 0 or p > 1:
            warnings.warn(
                "The p computed ({}) is outside the accepted range. ".format(p)
                + "Falling back to effective cross-linker functionality."
            )
            p = compute_effective_crosslinker_functionality(
                network, crosslinker_type)
    if f is None:
        if functionality_per_type is None:
            functionality_per_type = network.determine_functionality_per_type()
        if crosslinker_type not in functionality_per_type:
            raise ValueError(
                "The cross-linker functionality could not be determined. "
                + "Please pass it explicitly (`f`, or in `functionality_per_type`)."
            )
        f = functionality_per_type[crosslinker_type]
    if nu is None:
        nu = len(network.get_molecules(crosslinker_type)) / (
            network.get_volume() * unit_style.get_base_unit_of("volume")
        )
    if temperature is None:
        temperature = (273.15 + 25) * unit_style.get_underlying_unit_registry()(
            "kelvin"
        )
    if g_e_1 is None:
        g_e_1 = (
            8.3145  # gas constant, J/(mol*K)
            * temperature.to("kelvin").magnitude  # Temperature in Kelvin
            * 1e-6
            * 94.79281
        ) * unit_style.get_underlying_unit_registry()("MPa")
        # -> MPa, melt entanglement modulus

    # affine
    g_anm = nu * unit_style.kB * temperature
    # phantom
    g_pnm = (1 - 2 / f) * nu * unit_style.kB * temperature if f != 0 else 0.0
    # MMT:
    alpha, beta = compute_miller_macosko_probabilities(r, p, f)
    gamma_mmt_sum = 0.0
    for m in range(3, f + 1):
        gamma_mmt_sum += ((m - 2) / 2) * compute_probability_that_monomer_is_effective(
            f, m, alpha
        )
    gamma_mmt = (2 * r / f) * gamma_mmt_sum if f != 0 else 0.0
    g_mmt_phantom = gamma_mmt * nu * unit_style.kB * temperature
    # fraction of elastically effective strands.
    g_mmt_entanglement = g_e_1 * \
        compute_trapping_factor(p=p, r=r, f=f, alpha=alpha)
    # entanglement part. TODO : check adjustment with r (and where the 0.22 is
    # coming from? Fabian' s fit!)
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
    unit_style: Union[None, UnitStyle] = None,
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
        - unit_style: the units used, needed for temperature if not defined
    """
    if temperature is None:
        temperature = (273.15 + 25) * unit_style.get_underlying_unit_registry()(
            "kelvin"
        )
    if alpha is None:
        alpha, _ = compute_miller_macosko_probabilities(r, p, f)

    junction_part = (1 - w_sol) ** (-1 / 3) * compute_junction_modulus(
        p=p,
        r=r,
        xlink_concentration_0=xlink_concentration_0,
        unit_style=unit_style,
        f=f,
        alpha=alpha,
        temperature=temperature,
    )
    entanglement_part = (1 - w_sol) ** (-2) * compute_entanglement_modulus(
        p=p,
        r=r,
        f=f,
        g_e_1=g_e_1,
        alpha=alpha,
        temperature=temperature,
        unit_style=unit_style,
    )
    return junction_part + entanglement_part


def compute_entanglement_modulus(
    p: float,
    r: float,
    f: int,
    g_e_1: pint.Quantity,
    alpha: Union[float, None] = None,
    temperature: pint.Quantity = None,
    unit_style: Union[None, UnitStyle] = None,
):
    """
    Compute MMT's entanglement contribution to the equilibrium shear modulus, given by
    :math:`k_B T \\epsilon_e T_e`.

    Arguments:
        - p: the cross-linker conversion
        - r: the stoichiometric imbalance
        - f: the functionality of the crosslinkers
        - g_e_1: the melt entanglement modulus :math:`G_e(1) = k_B T \\epsilon_e`
        - alpha: :math:`P(F_a^{out})`, optional
        - temperature: the temperatures; defaults to room temperature
        - unit_style: the units used, needed for temperature if not defined
    """
    if temperature is None:
        temperature = (273.15 + 25) * unit_style.get_underlying_unit_registry()(
            "kelvin"
        )
    if alpha is None:
        alpha, _ = compute_miller_macosko_probabilities(r, p, f)
    return compute_trapping_factor(p=p, r=r, f=f, alpha=alpha) * g_e_1


def compute_junction_modulus(
    p: float,
    r: float,
    xlink_concentration_0: pint.Quantity,
    unit_style: UnitStyle,
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
        - unit_style: the units used, for example for k_B
        - f: the functionality of the crosslinkers
        - alpha: :math:`P(F_a^{out})`, optional
        - temperature: the temperatures; defaults to room temperature
    """
    if temperature is None:
        temperature = (273.15 + 25) * unit_style.get_underlying_unit_registry()(
            "kelvin"
        )
    if alpha is None:
        alpha, _ = compute_miller_macosko_probabilities(r, p, f)
    gamma_mmt_sum = 0.0
    for m in range(3, f + 1):
        gamma_mmt_sum += ((m - 2) / 2) * compute_probability_that_monomer_is_effective(
            f, m, alpha
        )

    return unit_style.kB * temperature * xlink_concentration_0 * gamma_mmt_sum


def compute_trapping_factor(
    network: Universe = None,
    p: float = None,
    r: float = None,
    f: Union[int, None] = None,
    alpha: Union[float, None] = None,
    beta: Union[float, None] = None,
    crosslinker_type: int = 2,
) -> float:
    """
    Compute the Langley trapping factor :math:`T_e`.
    Not all parameters are required; they will be computed from the network if needed.
    Provide all parameters if you don't have the network.

    Literature: https://doi.org/10.1021/ma60004a015

    Arguments:
        - network: the network to compute the trapping factor for. 
        - p: the extent of reaction in terms of the crosslinkers.
        - r: the stoichiometric imbalance of reactants.
        - f: functionality of the crosslinkers. Only needed if alpha is None.
        - alpha: :math:`P(F_a^{out})`, see :func:`~pylimer_tools.calc.miller_macosko_theory.compute_mms_probabilities()`
        - beta: :math:`P(F_b^{out})`, see :func:`~pylimer_tools.calc.miller_macosko_theory.compute_mms_probabilities()`
    """
    if f is None:
        _require_network(network, "f")
        f = network.determine_functionality_per_type()[crosslinker_type]

    if p is None:
        _require_network(network, "p")
        p = compute_crosslinker_conversion(
            network, crosslinker_type=crosslinker_type, f=f)
    if r is None:
        _require_network(network, "r")
        r = compute_stoichiometric_imbalance(network, crosslinker_type=crosslinker_type, functionality_per_type={
            crosslinker_type: f
        })

    if alpha is None or beta is None:
        alpha, beta = compute_miller_macosko_probabilities(r, p, f)

    if p == 0 or r == 0:
        return 0.0

    # for long B2s reacting with small A_fs
    # return (1 - beta)**4
    pel = ((1 / (p)) * (1 - alpha)) ** 2
    return pel**2


def compute_probability_that_monomer_is_effective(
    functionality_of_monomer: int, expected_degree_of_effect: int, p_f_a_out: float
):
    """
    Compute the probability that an Af, monomer will be an effective cross-link of degree m

    :math:`P(X_m^f) = \binom{f}{m} [P(F_A^{out})]^{f-m}[1-P(F_A^{out})]^m`

    Source:
        - Eq. 45 in Miller, Macosko 1976, A New Derivation of Post Gel Properties of Network

    Arguments:
        - functionality_of_monomer: f
        - expected_decree_of_effect: m
        - alpha: :math:`P(F_A^{out})`
    """
    f = functionality_of_monomer
    m = expected_degree_of_effect
    alpha = p_f_a_out
    return scipy.special.binom(f, m) * (alpha ** (f - m)) * ((1.0 - alpha) ** m)


def predict_gelation_point(r: float, f: int, g: int = 2) -> float:
    """
    Compute the gelation point :math:`p_{gel}` as theoretically predicted
    (gelation point = critical extent of reaction for gelation)

    Source:
      - https://www.sciencedirect.com/science/article/pii/003238618990253X

    Arguments:
      - r (double): the stoichiometric imbalance of reactants (see: #compute_stoichiometric_imbalance)
      - f (int): functionality of the crosslinkers
      - g (int): functionality of the precursor polymer

    Returns:
      - p_gel: critical extent of reaction for gelation
    """
    # if (r is None):
    #   r = calculateEffectiveCrosslinkerFunctionality(network, crosslinker_type, f)
    return math.sqrt(1 / (r * (f - 1) * (g - 1)))


def predict_p_from_w_sol(
    w_sol: float, r: float, w_f: float, w_g: float, f: int, g: int = 2
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
        try:
            p_f_a_out, _ = compute_miller_macosko_probabilities(r, p, f)
        except ValueError:
            p_f_a_out = 1.0  # highest value -> this will not be the optimum
        return (
            w_f * p_f_a_out**f + w_g *
            (r * p * p_f_a_out ** (f - 1) + 1 - r * p) ** g
        )

    res = optimize.minimize_scalar(
        lambda p: abs(w_sol - compute_wsol(p)), bounds=[1e-3, 1.0 - 1e-3]
    )
    if not res.success:
        warnings.warn("The p predicted from w_sol might be incorrect")
    return res.x


def _require_network(network: Universe = None, instead_of: str = ""):
    if (network is None):
        if (instead_of == ""):
            raise ValueError("A network is required")
        else:
            raise ValueError(
                "If `{}` is not specified, a network is required".format(instead_of))
