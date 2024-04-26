# source: https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

from __future__ import annotations

from typing import Iterable

import numpy as np

from pylimer_tools.calc.structure_analysis import \
    compute_mean_end_to_end_distances
from pylimer_tools_cpp.pylimer_tools_cpp import MoleculeType, Universe


def predict_shear_modulus(networks: Iterable[Universe], temperature: float = 1, k_boltzmann: float = 1,
                          crossLinker_type: int = 2, total_mass=1) -> float:
    """
    Predict the shear modulus using ANT Analysis.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the polymer system to predict the shear modulus for
      - T: the temperature in your unit system
      - k_b: Boltzmann's constant in your unit system
      - crossLinker_type: the type of atoms to ignore (junctions, crossLinkers)
      - totalMass: the :math:`M` in the respective formula

    Returns:
      - shear modulus (float): the estimated shear modulus. Unit: [pressure]
    """
    gamma = compute_topological_factor(
        networks, crossLinker_type, total_mass)
    nu = 0
    for network in networks:
        nu += len(network.getMolecules(crossLinker_type)) / \
            (network.getVolume()) / len(networks)
    return gamma * nu * k_boltzmann * temperature


def compute_cycle_rank(networks: Iterable[Universe] = None, nu: int = None, mu: int = None,
                       abs_tol: float = 1, rel_tol: float = 1, crossLinker_type: int = 2) -> float:
    """
    Compute the cycle rank (:math:`\\chi`).
    Assumes the precursor-chains to be bifunctional.

    Arguments:
      - network: the network to calculate the cycle rank for
      - nu: number of elastically effective (active) strands per unit volume
      - mu: number density of the elastically effective crosslink
      - abs_tol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance)
            (None to use only rel_tol)
      - rel_tol (float): the relative tolerance to categorize a chain as active
            (0: all, 1: none (use only abs_tol))
      - crossLinker_type: the atom type of the crossLinkers/junctions

    No need to provide all the parameters — either/or:
    - nu & mu
    - network, abs_tol, rel_tol, crossLinker_type

    Returns:
      - cycleRank: the cycle rank ($\\xi = \\nu_{eff} - \\mu_{eff}$). Unit: [1/Volume]
    """
    if (nu is None):
        if (crossLinker_type is None or networks is None):
            raise ValueError(
                "Argument missing: When not specifying nu, network and crossLinker_type need to be specified")
        nu = compute_effective_nr_density_of_network(
            networks, abs_tol, rel_tol, crossLinker_type)
    if (mu is None):
        if (crossLinker_type is None or networks is None):
            raise ValueError(
                "Argument missing: When not specifying mu, network and crossLinker_type need to be specified")
        mu = compute_effective_nr_density_of_junctions(
            networks, abs_tol, rel_tol, crossLinker_type)

    return nu - mu


def compute_effective_nr_density_of_network(networks: Iterable[Universe], abs_tol: float = 1, rel_tol: float = 1,
                                            crossLinker_type: int = 2) -> float:
    """
    Compute the effective number density :math:`\\nu_{eff}` of a network.
    Assumes the precursor-chains to be bifunctional.

    :math:`\\nu_{eff}` is the number of elastically effective (active) strands per unit volume,
    which are defined as the ones that can store elastic energy
    upon network deformation, resp. the effective number density of network strands

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network (pylimer_tools.entities.Universe): the network to compute :math:`\\nu_{eff}` for
      - abs_tol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance)
            (None to use only rel_tol)
      - rel_tol (float): the relative tolerance to categorize a chain as active
            (0: all, 1: none (use only abs_tol))
      - crossLinker_type: the atom type of the cross-linkers/junctions

    Returns:
      - :math:`\\nu_{eff}` (float): the effective number density of network strands. Unit: [1/Volume]
    """
    if (len(networks) == 0):
        return None

    # get the mean end to end distances
    r_taus = compute_mean_end_to_end_distances(networks, crossLinker_type)
    if (len(r_taus) < 1):
        return 0.0
    r_taus = np.array(list(r_taus.values()))
    r_tau_max = np.max(r_taus)

    # process additional input parameters
    if (abs_tol is None):
        abs_tol = r_tau_max

    # count how many effective strands there are
    num_effective = np.array([r_tau > abs_tol or r_tau > rel_tol * r_tau_max
                             for r_tau in r_taus]).sum()
    mean_volume = compute_mean_universe_volume(networks)

    return num_effective / mean_volume


def compute_mean_universe_volume(networks: Iterable[Universe], accept_different_sizes: bool = False) -> float:
    """
    Compute the mean volume of a list of universes.

    Arguments:
      - networks: a list of universes
      - accept_different_sizes: toggle whether to throw an error when the Universe have different nr. of atoms

    Returns:
      - mean_volume (float): the mean volume of the universes
    """
    if (len(networks) < 1):
        raise ValueError('Must have at least one network')
    # compute the mean volume of the universes
    mean_volume = 0
    divisor = 1 / len(networks)
    network_size = networks[0].getNrOfAtoms()
    for network in networks:
        if (not accept_different_sizes and network.getNrOfAtoms() != network_size):
            raise NotImplementedError(
                "Currently, only sequences of networks with the same size are supported"
                + " (got one with {} instead of {})".format(
                    network.getNrOfAtoms(), network_size
                ))
        mean_volume += network.getVolume() * divisor
    return mean_volume


def compute_effective_nr_density_of_junctions(networks: Iterable[Universe], abs_tol: float = 0, rel_tol: float = 1,
                                              crossLinker_type: int = 2, min_num_effective_strands=2) -> float:
    """
    Compute the number density of the elastically effective crosslinks,
    defined as the ones that connect at least `min_num_effective_strands` elastically effective strands.
    Assumes the precursor-chains to be bifunctional.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network (pylimer_tools.entities.Universe): the network to compute :math:`\\nu_{eff}` for
      - abs_tol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance)
          (None to use only rel_tol)
      - rel_tol (float): the relative tolerance to categorize a chain as active
          (0: all, 1: none (use only abs_tol))
      - crossLinker_type: the atom type of the crossLinkers/junctions
      - min_num_effective_strands (int): the number of elastically effective strands to qualify a junction as such

    Returns:
      - :math:`\\mu_{eff}` (float): the effective number density of junctions. Unit: [1/Volume]
    """
    if (len(networks) < 1):
        return None
    if (crossLinker_type is None):
        return 0.0

    mean_volume = compute_mean_universe_volume(networks)

    if (min_num_effective_strands == 0):
        return len(networks[0].getAtomsOfType(crossLinker_type)) / mean_volume

    # get the mean end to end distances
    r_taus = compute_mean_end_to_end_distances(networks, crossLinker_type)
    if (len(r_taus) < 1):
        return 0.0
    r_tau_max = max(r_taus.values())

    # process additional input parameters
    if (abs_tol is None):
        abs_tol = r_tau_max

    key_to_molecule = {}
    for molecule in list(networks)[0].getChainsWithCrosslinker(crossLinker_type):
        key_to_molecule[molecule.getKey()] = molecule

    # count how many active connections each junction has
    junction_activity = {}
    for key in r_taus:
        crossLinkers = key_to_molecule[key].getAtomsOfType(crossLinker_type)
        assert (len(crossLinkers) == 2)
        is_active = r_taus[key] > abs_tol or r_taus[key] > rel_tol * r_tau_max
        if (not (is_active)):
            continue
        relevant_names = [crossLinkers[0].getId(), crossLinkers[1].getId()]
        for crossLinker_name in relevant_names:
            if (crossLinker_name not in junction_activity):
                junction_activity[crossLinker_name] = 0
            junction_activity[crossLinker_name] += 1

    effective_junctions = np.array(
        [junction_activity[key] >= min_num_effective_strands for key in junction_activity])
    num_effective_junctions = effective_junctions.sum()
    return num_effective_junctions / mean_volume


def compute_topological_factor(networks: Iterable[Universe], crossLinker_type: int = 2,
                               total_mass: float = 1, b: float | None = None) -> float:
    """
    Compute the topological factor of a polymer network.

    Assumptions:
      - the precursor-chains to be bifunctional
      - all Universes to have the same structure (with possibly differing positions)
      - crossLinkers do not count to the nr. of monomers in a strand

    Source:
      - eq. 16 in https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the network to compute the topological factor for
      - crossLinker_type: the type of atoms to ignore
      - total_mass: the :math:`M` in the respective formula
      - b: the mean bond length.
          If `None`, it will be computed for each molecule in the first Universe (Network).

    Returns:
      - the topological factor :math:`\\Gamma`
    """
    r_taus = compute_mean_end_to_end_distances(networks, crossLinker_type)

    # find the topological factor
    gamma_sum = 0
    network = networks[0]  # this is where the second assumption is made
    chains_to_process = network.getChainsWithCrosslinker(crossLinker_type)
    for molecule in chains_to_process:
        crossLinkers = molecule.getAtomsOfType(crossLinker_type)
        if (len(crossLinkers) != 2 or
            molecule.getType() == MoleculeType.PRIMARY_LOOP or
                molecule.getType() == MoleculeType.DANGLING_CHAIN):
            # dangling, free chains and loops are irrelevant for our purposes
            continue
        if (b is None):
            b = np.mean(molecule.computeBondLengths())
        crossLinkers = [crossLinkers[0], crossLinkers[1]]
        # sort crossLinkers by name as a way to keep the vector directions consistent between timesteps
        crossLinkers.sort(key=lambda a: a.getId())
        key = molecule.getKey()
        gamma_sum += r_taus[key] * r_taus[key] / \
            ((molecule.getLength() - 2) * b *
             b)  # -2: remove crossLinkers again (assumption 3)

    return gamma_sum / total_mass
