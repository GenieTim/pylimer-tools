"""
This module is deprecated.
Use :func:`pylimer_tools_cpp.MEHPForceBalance` or :func:`pylimer_tools_cpp.MEHPForceRelaxation` instead.

In principle, it offers comparable functionality, but without the force minimization - the networks you pass
to these methods should be minimized first.
"""

# source: https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence, Union

import numpy as np

from pylimer_tools.calc.structure_analysis import compute_mean_end_to_end_distances
from pylimer_tools_cpp import MoleculeType

if TYPE_CHECKING:
    from pylimer_tools_cpp import Universe


def predict_shear_modulus(
    networks: Sequence[Universe],
    temperature: float = 1,
    k_boltzmann: float = 1,
    crosslinker_type: int = 2,
    total_mass=1,
) -> float:
    """
    Predict the shear modulus using ANT Analysis.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    :param networks: The polymer systems to predict the shear modulus for
    :type networks: Sequence[Universe]
    :param temperature: The temperature in your unit system
    :type temperature: float
    :param k_boltzmann: Boltzmann's constant in your unit system
    :type k_boltzmann: float
    :param crosslinker_type: The type of atoms to ignore (junctions, crosslinkers)
    :type crosslinker_type: int
    :param total_mass: The :math:`M` in the respective formula
    :type total_mass: float

    :return: The estimated shear modulus in pressure units
    :rtype: float
    """
    gamma = compute_topological_factor(networks, crosslinker_type, total_mass)
    nu = 0
    for network in networks:
        nu += (
            len(network.get_molecules(crosslinker_type))
            / (network.get_volume())
            / len(networks)
        )
    return gamma * nu * k_boltzmann * temperature


def compute_cycle_rank(
    networks: Union[Sequence[Universe], None] = None,
    nu: Union[float, None] = None,
    mu: Union[float, None] = None,
    abs_tol: float = 1,
    rel_tol: float = 1,
    crosslinker_type: int = 2,
) -> float:
    """
    Compute the cycle rank (:math:`\\chi`).
    Assumes the precursor-chains to be bifunctional.

    No need to provide all the parameters — either/or:

    * nu & mu
    * network, abs_tol, rel_tol, crosslinker_type

    :param networks: The networks to calculate the cycle rank for
    :type networks: Union[Sequence[Universe], None]
    :param nu: Number of elastically effective (active) strands per unit volume
    :type nu: Union[float, None]
    :param mu: Number density of the elastically effective crosslinks
    :type mu: Union[float, None]
    :param abs_tol: The absolute tolerance to categorize a chain as active (min. end-to-end distance).
                    Set to None to use only rel_tol
    :type abs_tol: float
    :param rel_tol: The relative tolerance to categorize a chain as active (0: all, 1: none)
    :type rel_tol: float
    :param crosslinker_type: The atom type of the crosslinkers/junctions
    :type crosslinker_type: int

    :return: The cycle rank (:math:`\\xi = \\nu_{eff} - \\mu_{eff}`) in units of 1/Volume
    :rtype: float
    """
    if nu is None:
        if crosslinker_type is None or networks is None:
            raise ValueError(
                "Argument missing: When not specifying nu, network and crosslinker_type need to be specified"
            )
        nu = compute_effective_nr_density_of_network(
            networks, abs_tol, rel_tol, crosslinker_type
        )
    if mu is None:
        if crosslinker_type is None or networks is None:
            raise ValueError(
                "Argument missing: When not specifying mu, network and crosslinker_type need to be specified"
            )
        mu = compute_effective_nr_density_of_junctions(
            networks, abs_tol, rel_tol, crosslinker_type
        )

    return nu - mu


def compute_effective_nr_density_of_network(
    networks: Sequence[Universe],
    abs_tol: float = 1,
    rel_tol: float = 1,
    crosslinker_type: int = 2,
) -> float:
    """
    Compute the effective number density :math:`\\nu_{eff}` of a network.
    Assumes the precursor-chains to be bifunctional.

    :math:`\\nu_{eff}` is the number of elastically effective (active) strands per unit volume,
    which are defined as the ones that can store elastic energy
    upon network deformation, resp. the effective number density of network strands.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    :param networks: The networks to compute :math:`\\nu_{eff}` for
    :type networks: Sequence[Universe]
    :param abs_tol: The absolute tolerance to categorize a chain as active (min. end-to-end distance).
                    Set to None to use only rel_tol
    :type abs_tol: float
    :param rel_tol: The relative tolerance to categorize a chain as active (0: all, 1: none)
    :type rel_tol: float
    :param crosslinker_type: The atom type of the crosslinkers/junctions
    :type crosslinker_type: int

    :return: The effective number density of network strands in units of 1/Volume
    :rtype: float
    """
    if len(networks) == 0:
        return 0

    # get the mean end to end distances
    r_taus = compute_mean_end_to_end_distances(networks, crosslinker_type)
    if len(r_taus) < 1:
        return 0.0
    r_taus = np.array(list(r_taus.values()))
    r_tau_max = np.max(r_taus)

    # process additional input parameters
    if abs_tol is None:
        abs_tol = r_tau_max

    # count how many effective strands there are
    num_effective = np.array(
        [r_tau > abs_tol or r_tau > rel_tol * r_tau_max for r_tau in r_taus]
    ).sum()
    mean_volume = compute_mean_universe_volume(networks)

    return num_effective / mean_volume


def compute_mean_universe_volume(
    networks: Sequence[Universe], accept_different_sizes: bool = False
) -> float:
    """
    Compute the mean volume of a list of universes.

    :param networks: A list of universes
    :type networks: Sequence[Universe]
    :param accept_different_sizes: Toggle whether to throw an error when the Universes have different numbers of atoms
    :type accept_different_sizes: bool

    :return: The mean volume of the universes
    :rtype: float

    :raises ValueError: If the input list is empty
    :raises NotImplementedError: If universes have different sizes and accept_different_sizes is False
    """
    if len(networks) < 1:
        raise ValueError("Must have at least one network")
    # compute the mean volume of the universes
    mean_volume = 0
    divisor = 1 / len(networks)
    network_size = networks[0].get_nr_of_atoms()
    for network in networks:
        if not accept_different_sizes and network.get_nr_of_atoms() != network_size:
            raise NotImplementedError(
                "Currently, only sequences of networks with the same size are supported"
                + " (got one with {} instead of {})".format(
                    network.get_nr_of_atoms(), network_size
                )
            )
        mean_volume += network.get_volume() * divisor
    return mean_volume


def compute_effective_nr_density_of_junctions(
    networks: Sequence[Universe],
    abs_tol: float = 0,
    rel_tol: float = 1,
    crosslinker_type: int = 2,
    min_num_effective_strands=2,
) -> float:
    """
    Compute the number density of the elastically effective crosslinks,
    defined as the ones that connect at least `min_num_effective_strands` elastically effective strands.
    Assumes the precursor-chains to be bifunctional.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    :param networks: The networks to compute :math:`\\mu_{eff}` for
    :type networks: Sequence[Universe]
    :param abs_tol: The absolute tolerance to categorize a chain as active (min. end-to-end distance).
                    Set to None to use only rel_tol
    :type abs_tol: float
    :param rel_tol: The relative tolerance to categorize a chain as active (0: all, 1: none)
    :type rel_tol: float
    :param crosslinker_type: The atom type of the crosslinkers/junctions
    :type crosslinker_type: int
    :param min_num_effective_strands: The number of elastically effective strands to qualify a junction as such
    :type min_num_effective_strands: int

    :return: The effective number density of junctions in units of 1/Volume
    :rtype: float
    """
    if len(networks) < 1:
        return 0.0
    if crosslinker_type is None:
        return 0.0

    mean_volume = compute_mean_universe_volume(networks)

    if min_num_effective_strands == 0:
        return len(networks[0].get_atoms_by_type(
            crosslinker_type)) / mean_volume

    # get the mean end to end distances
    r_taus = compute_mean_end_to_end_distances(networks, crosslinker_type)
    if len(r_taus) < 1:
        return 0.0
    r_tau_max = max(r_taus.values())

    # process additional input parameters
    if abs_tol is None:
        abs_tol = r_tau_max

    key_to_molecule = {}
    for molecule in list(networks)[
            0].get_chains_with_crosslinker(crosslinker_type):
        key_to_molecule[molecule.get_key()] = molecule

    # count how many active connections each junction has
    junction_activity = {}
    for key in r_taus:
        crosslinkers = key_to_molecule[key].get_atoms_by_type(crosslinker_type)
        assert len(crosslinkers) == 2
        is_active = r_taus[key] > abs_tol or r_taus[key] > rel_tol * r_tau_max
        if not (is_active):
            continue
        relevant_names = [crosslinkers[0].get_id(), crosslinkers[1].get_id()]
        for crosslinker_name in relevant_names:
            if crosslinker_name not in junction_activity:
                junction_activity[crosslinker_name] = 0
            junction_activity[crosslinker_name] += 1

    effective_junctions = np.array(
        [
            junction_activity[key] >= min_num_effective_strands
            for key in junction_activity
        ]
    )
    num_effective_junctions = effective_junctions.sum()
    return num_effective_junctions / mean_volume


def compute_topological_factor(
    networks: Sequence[Universe],
    crosslinker_type: int = 2,
    total_mass: float = 1,
    b: float | None = None,
) -> float:
    """
    Compute the topological factor of a polymer network.

    Assumptions:
      - the precursor-chains to be bifunctional
      - all Universes to have the same structure (with possibly differing positions)
      - crosslinkers do not count to the nr. of monomers in a strand

    Source:
      - eq. 16 in https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    :param networks: The networks to compute the topological factor for
    :type networks: Sequence[Universe]
    :param crosslinker_type: The type of atoms to ignore (junctions, crosslinkers)
    :type crosslinker_type: int
    :param total_mass: The :math:`M` in the respective formula
    :type total_mass: float
    :param b: The mean bond length. If None, it will be computed for each molecule in the first Universe
    :type b: float or None

    :return: The topological factor :math:`\\Gamma`
    :rtype: float
    """
    r_taus = compute_mean_end_to_end_distances(networks, crosslinker_type)

    # find the topological factor
    gamma_sum = 0
    network = networks[0]  # this is where the second assumption is made
    chains_to_process = network.get_chains_with_crosslinker(crosslinker_type)
    for molecule in chains_to_process:
        crosslinkers = molecule.get_atoms_by_type(crosslinker_type)
        if (
            len(crosslinkers) != 2
            or molecule.get_strand_type() == MoleculeType.PRIMARY_LOOP
            or molecule.get_strand_type() == MoleculeType.DANGLING_CHAIN
        ):
            # dangling, free chains and loops are irrelevant for our purposes
            continue
        if b is None:
            b = float(np.mean(molecule.compute_bond_lengths()))
        crosslinkers = [crosslinkers[0], crosslinkers[1]]
        # sort crosslinkers by name as a way to keep the vector directions
        # consistent between timesteps
        crosslinkers.sort(key=lambda a: a.get_id())
        key = molecule.get_key()
        gamma_sum += (
            r_taus[key] * r_taus[key] /
            ((molecule.get_nr_of_atoms() - 2) * b * b)
        )  # -2: remove crosslinkers again (assumption 3)

    return gamma_sum / total_mass
