"""
Module for calculating stoichiometric relationships in polymer networks.

In particular, it helps to work with the parameter `b2`, which is the mole fraction
of reactive sites in B2 among all reactive sites in a mixture of B1 and B2.
"""

import pint
from typing import Union, cast
import math


def compute_number_fractions(
    b2: float = 1.0,
) -> tuple[float, float]:
    """
    Compute the number fractions of bifunctional and monofunctional strands.

    :param b2: The mole fraction of reactive sites in B2 among all reactive sites in a mixture of B1 and B2
    :return: A tuple containing the number fractions (x_bifunctional, x_monofunctional)
    :rtype: tuple[float, float]
    """
    assert 0.0 <= b2 <= 1.0, "b2 must be between 0 and 1."

    n_total_strands = (
        1e3  # irrelevant, cancels out, but makes the calculation easier to read
    )
    """
    b2 = (2 * n_bifunctional) / (n_monofunctional + 2 * n_bifunctional)
    n_total = n_bifunctional + n_monofunctional
    """
    n_bifunctional = -b2 * n_total_strands / (-2 + b2)
    n_monofunctional = n_total_strands - n_bifunctional
    assert math.isclose(
        n_monofunctional,
        2 * (-n_total_strands + b2 * n_total_strands) / (-2 + b2),
        abs_tol=1e-6,
    )

    bifunctional_fraction = n_bifunctional / n_total_strands
    monofunctional_fraction = n_monofunctional / n_total_strands

    assert math.isclose(
        bifunctional_fraction + monofunctional_fraction,
        1.0,
        abs_tol=1e-6,
    )

    return bifunctional_fraction, monofunctional_fraction


def compute_strand_number_density(
    mw_bifunctional: pint.Quantity,
    density: pint.Quantity,
    mw_monofunctional: Union[pint.Quantity, None] = None,
    b2: float = 1.0,
) -> pint.Quantity:
    """
    Compute the strand number density of the polymer network.

    :param mw_bifunctional: The molecular weight of the bifunctional strands
    :param density: The density of the polymer network
    :param mw_monofunctional: The molecular weight of the monofunctional strands. None if not present
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites in a mixture of B1 and B2
    :return: The strand number density
    :rtype: pint.Quantity
    """
    if (
        mw_monofunctional is None
        or not math.isfinite(mw_monofunctional.magnitude)
        or math.isclose(mw_monofunctional.magnitude, 0.0)
    ):
        assert b2 == 1.0, "If no monofunctional strands are present, b2 must be 1.0."
        return density / mw_bifunctional

    bifunctional_fraction, monofunctional_fraction = compute_number_fractions(
        b2)

    return cast(
        pint.Quantity,
        (
            bifunctional_fraction * density / mw_bifunctional
            + monofunctional_fraction * density / mw_monofunctional
        ),
    )


def compute_weight_fractions(
    mw_bifunctional: pint.Quantity,
    mw_monofunctional: Union[pint.Quantity, None] = None,
    b2: float = 1.0,
) -> tuple[float, float]:
    """
    Compute the weight fractions of bifunctional and monofunctional strands.
    Makes some assumptions, e.g. the same density for both types of strands.

    :param mw_bifunctional: The molecular weight of the bifunctional strands
    :param mw_monofunctional: The molecular weight of the monofunctional strands. None if not present
    :param b2: The mole fraction of reactive sites in B2 among all reactive sites in a mixture of B1 and B2
    :return: A tuple containing the weight fractions (w_bifunctional, w_monofunctional)
    :rtype: tuple[float, float]
    """
    if mw_monofunctional is None or not math.isfinite(
            mw_monofunctional.magnitude):
        assert b2 == 1.0, "If no monofunctional strands are present, b2 must be 1.0."
        return 1.0, 0.0

    assert 0.0 <= b2 <= 1.0, "b2 must be between 0 and 1."

    bifunctional_fraction, monofunctional_fraction = compute_number_fractions(
        b2)

    mass_bifunctional = bifunctional_fraction * mw_bifunctional.magnitude
    mass_monofunctional = monofunctional_fraction * mw_monofunctional.magnitude
    total_mass = mass_bifunctional + mass_monofunctional

    w_bifunctional = mass_bifunctional / total_mass
    w_monofunctional = mass_monofunctional / total_mass

    return w_bifunctional, w_monofunctional
