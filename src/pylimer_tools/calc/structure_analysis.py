from __future__ import annotations

import math
import warnings
from collections import Counter
from typing import Iterable, Tuple

import numpy as np

from pylimer_tools_cpp.pylimer_tools_cpp import MoleculeType, Universe


def compute_stoichiometric_imbalance(network: Universe, crossLinker_type: int = 2,
                                     functionality_per_type: dict = None, ignore_types: list = [],
                                     effective: bool = False) -> float:
    """
    Compute the stoichiometric imbalance
    ( nr. of bonds formable of crossLinker / (nr. of precursor chains * 2) )

    NOTE:
      if your system has a non-integer number of possible bonds (e.g. one site non-bonded),
      this will not be rounded/respected in any way.

    Arguments:
      - network: the polymer network to do the computation for
      - crossLinker_type: the type of the junctions/crossLinkers to select them in the network
      - functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
          If `None`: will use max functionality per type.
      - ignore_types: a list of integers, the types to ignore for the imbalance (e.g. solvent atom types)
      - effective: whether to use the effective functionality (if functionality_per_type is not passed) or the maximum

    Returns:
      - r (float): the stoichiometric imbalance
    """
    if (network.getNrOfAtoms() == 0):
        return 0.

    counts = Counter(network.getAtomTypes())

    if (functionality_per_type is not None and
            not np.all([k in functionality_per_type for k in counts if counts[k] > 0])):
        functionality_per_type = None
        warnings.warn("Not all atom types found in functionality_per_type, " +
                      "will ignore passed argument `functionality_per_type`.")

    if (functionality_per_type is None):
        functionality_per_type = network.determineEffectiveFunctionalityPerType(
        ) if effective else network.determineFunctionalityPerType()

    if (crossLinker_type not in counts):
        return 0.

    # TODO: use the data from the functionality_per_type to determine the
    # functionality per strand, maybe?
    strands = network.getMolecules(crossLinker_type)
    ignore_types.append(crossLinker_type)
    num_relevant_strands = len([m for m in strands if not np.all(
        [a.getType() in ignore_types for a in m.getAtoms()])])

    crossLinker_formable_bonds = counts[crossLinker_type] * \
        functionality_per_type[crossLinker_type]
    other_formable_bonds = num_relevant_strands * 2

    if (other_formable_bonds == 0):
        return math.inf

    # division by 2 is implicit
    return crossLinker_formable_bonds / (other_formable_bonds)


def compute_extent_of_reaction(network: Universe, crossLinker_type, functionality_per_type: dict = None) -> float:
    """
    Compute the extent of polymerization reaction
    (nr. of formed bonds in reaction / max. nr. of bonds formable)
    NOTE:
        - if your system has a non-integer number of possible bonds (e.g. one site non-bonded),
            this will not be rounded/respected in any way
        - if the system contains solvent or other molecules that should not be binding to
            cross-linkers, make sure to remove them before calling this function

    Arguments:
      - network: the polymer network to do the computation for
      - crossLinker_type: the atom type of crossLinker beads
      - functionality_per_type: a dictionary with key: type, and value: functionality of this atom type.
          If None: will use max functionality per type.

    Returns:
      - p (float): the extent of reaction
    """
    if (network.getNrOfAtoms() == 0):
        return 1

    if (functionality_per_type is not None and crossLinker_type not in functionality_per_type):
        functionality_per_type = None
        warnings.warn("Cross-linker type {} not found in passed functionality_per_type, ".format(crossLinker_type) +
                      "will ignore passed argument `functionality_per_type`.")

    if (functionality_per_type is None):
        functionality_per_type = network.determineFunctionalityPerType()

    num_strands = len(network.getMolecules(crossLinker_type))
    crosslinks = network.getAtomsOfType(crossLinker_type)
    num_crossLinkers = len(crosslinks)

    # assuming strand has functionality 2
    max_formable_bonds = min(num_strands * 2, num_crossLinkers *
                             functionality_per_type[crossLinker_type])

    if (max_formable_bonds == 0):
        return 1

    actually_formed_bonds = 0
    for crosslink in crosslinks:
        connected_to = network.getConnectedAtoms(crosslink)
        actually_formed_bonds += len(
            [a for a in connected_to if a.getType() != crossLinker_type])

    return actually_formed_bonds / (max_formable_bonds)


def compute_mean_end_to_end_distances(networks: Iterable[Universe], crossLinker_type: int = 2) -> dict:
    """
    Compute the mean end to end distance between each pair of (indirectly) connected crossLinker

    Arguments:
      - networks: the different configurations of the polymer network to do the computation for
      - crossLinker_type: the atom type to compute the in-between vectors for

    Returns:
      - endToEndDistances (dict): a dictionary with key: "{atom1.name}+{atom2.name}"
          and value: the norm of the mean difference vector
    """
    r_tau_vectors = compute_mean_end_to_end_vectors(networks, crossLinker_type)
    if (len(r_tau_vectors) < 1):
        return {}

    r_tau_vectors_array = np.array(list(r_tau_vectors.values()))
    r_taus = np.linalg.norm(r_tau_vectors_array, axis=1)

    return dict(zip(r_tau_vectors.keys(), r_taus))


def compute_mean_end_to_end_vectors(networks: Iterable[Universe], crossLinker_type: int = 2) -> dict:
    """
    Compute the mean end to end vectors between each pair of (indirectly) connected crossLinker

    Arguments:
      - networks: the different configurations of the polymer network to do the computation for
      - crossLinker_type: the atom type to compute the in-between vectors for

    Returns:
      - end_to_end_vectors (dict): a dictionary with key: "{atom1.name}+{atom2.name}"
          and value: their mean distance difference vector
    """
    if (len(networks) == 0):
        return {}
    end_to_end_vectors = {}
    key_counts = {}
    divider = 1 / len(networks)
    iteration = 0
    for network in networks:
        current_end_to_end_vectors = compute_end_to_end_vectors(
            network, crossLinker_type)
        # the mean calculation in this for loop
        # trades some memory for performance
        # there are still many performance and memory
        # improvements possible
        # (e.g. computing connectivity only once, storing it only once, ....)
        for key in current_end_to_end_vectors:
            if (key not in end_to_end_vectors):
                if (iteration > 0):
                    raise ValueError("Found molecule {} in network {}, but not in previous".format(
                        key, iteration
                    ))
                end_to_end_vectors[key] = [0, 0, 0]
                key_counts[key] = 0
            for i in range(3):
                end_to_end_vectors[key][i] += current_end_to_end_vectors[key][i] * divider
            key_counts[key] += 1
        iteration += 1
    if (len(key_counts) > 0 and not np.all([c == list(key_counts.values())[0] for c in key_counts.values()])):
        raise ValueError("The networks contain different molecules.")
    return end_to_end_vectors


def compute_end_to_end_vectors(network: Universe, crossLinker_type: int = 2) -> dict:
    """
    Compute the end to end vectors between each pair of (indirectly) connected crossLinker

    Arguments:
      - network: the polymer network to do the computation for
      - crossLinker_type: the atom type to compute the in-between vectors for

    Returns:
      - end_to_end_vectors (dict): a dictionary with key: "{atom1.name}+{atom2.name}"
          and value: their difference vector
    """
    # while we could do the decomposition again with explicit removal of irrelevant strand atoms,
    # this should not be any more expensive
    end_to_end_vectors = {}
    molecules = network.getChainsWithCrosslinker(crossLinker_type)
    for molecule in molecules:
        crossLinkers = molecule.getAtomsOfType(crossLinker_type)
        if (len(crossLinkers) != 2 or
            molecule.getType() == MoleculeType.PRIMARY_LOOP or
                molecule.getType() == MoleculeType.DANGLING_CHAIN):
            # dangling, free chains and loops are irrelevant for our purposes
            continue
        # igraph.VertexSeq is not sortable -> use a list
        crossLinkers = [crossLinkers[0], crossLinkers[1]]
        # sort crossLinkers by name as a way to keep the vector directions consistent between timesteps
        crossLinkers.sort(key=lambda a: a.getId())
        #
        end_to_end_vectors[molecule.getKey()] = crossLinkers[0].computeVectorTo(
            crossLinkers[1], network.getBox())

    return end_to_end_vectors


def compute_crossLinker_conversion(network: Universe, crossLinker_type: int = 2, f: int = None,
                                   functionality_per_type: dict = None) -> float:
    """
    Compute the extent of reaction of the crossLinkers
    (actual functionality divided by target functionality)

    Arguments:
      - network: the polymer network to do the computation for
      - crossLinker_type: the type of the junctions/crossLinkers to select them in the network
      - f: the functionality of the crossLinkers

    Returns:
      - r (float): the (mean) crossLinker conversion
    """
    if (f is None):
        if (functionality_per_type is None):
            functionality_per_type = network.determineFunctionalityPerType()
        if (crossLinker_type not in functionality_per_type):
            return 0.0
        f = functionality_per_type[crossLinker_type]

    if (f == 0.):
        warnings.warn("Crosslinker functionality = 0 is problematic.")

    return compute_effective_crossLinker_functionality(network, crossLinker_type) / f


def compute_effective_crossLinker_functionality(network: Universe, crossLinker_type: int = 2) -> float:
    """
    Compute the mean crossLinker functionality

    Arguments:
      - network: the polymer network to do the computation for
      - crossLinker_type: the type of the junctions/crossLinkers to select them in the network

    Returns:
      - f (float): the (mean) effective crossLinker functionality
    """
    junction_degrees = compute_effective_crossLinker_functionalities(
        network, crossLinker_type)
    return np.mean(junction_degrees) if len(junction_degrees) > 0 else 0.


def compute_effective_crossLinker_functionalities(network: Universe, crossLinker_type: int = 2) -> list[int]:
    """
    Compute the functionality of every crossLinker in the network

    Arguments:
      - network: the polymer network to do the computation for
      - crossLinker_type: the type of the junctions/crossLinkers to select them in the network

    Returns:
      - junctionDegrees (list[int]): the functionality of every crossLinker
    """
    if (network.getNrOfAtoms() == 0):
        return []
    junctions = network.getAtomsOfType(crossLinker_type)
    junction_ids = [v.getId() for v in junctions]
    junction_degrees = [network.getNrOfBondsOfAtom(id) for id in junction_ids]
    return junction_degrees


def compute_weight_fractions(network: Universe) -> dict:
    """
    Compute the weight fractions of each atom type in the network.

    Arguments:
      - network: the polymer network to do the computation for

    Returns:
      - :math:`\\vec{W_i}` (dict): using the type i as a key,
            this dict contains the weight fractions (:math:`\\frac{W_i}{W_{tot}}`)
    """
    return network.computeWeightFractions()


def measure_weight_fraction_of_backbone(network: Universe, crossLinker_type: int = 2):
    """
    Compute the weight fraction of network backbone in infinite network

    Arguments:
      - network: the network to compute the weight fraction for
      - crossLinker_type: the atom type to use to split the molecules

    Returns:
      - weightFraction (float): 1 - weightDangling/weightTotal,
    """
    if (network.getNrOfAtoms() < 1):
        return 0.0

    weight_fraction, _ = measure_weight_fraction_of_dangling_chains(
        network, crossLinker_type)
    return 1.0 - weight_fraction


def measure_weight_fraction_of_dangling_chains(network: Universe, crossLinker_type: int = 2) -> Tuple[float, float]:
    """
    Compute the weight fraction of dangling strands in infinite network

    NOTE:
        Currently, only primary dangling chains are taken into account.
        There are other methods that incorporate more.

    Arguments:
      - network: the network to compute the weight fraction for
      - crossLinker_type: the atom type to use to split the molecules

    Returns:
      - weightFraction: weightDangling/weightTotal,
      - numFraction: numDangling/numTotal
    """
    if (network.getNrOfAtoms() < 1):
        return 0.0, 0.0

    weights = network.getMasses()

    def get_weight_of_graph(graph):
        counts = Counter(graph.getAtomTypes())
        weight_total = 0
        for key in counts:
            weight_total += weights[key] * counts[key]
        return weight_total

    all_chains = network.getChainsWithCrosslinker(crossLinker_type)
    num_total = network.getNrOfAtoms()
    weight_total = get_weight_of_graph(network)

    num_dangling = 0
    weight_dangling = 0
    for chain in all_chains:
        if (chain.getType() == MoleculeType.DANGLING_CHAIN):
            num_dangling += chain.getNrOfAtoms()
            weight_dangling += get_weight_of_graph(chain)

    if (weight_total == 0):
        # warnings.warn("Total weight of network is = 0.")
        return 0.0, num_dangling / num_total

    return weight_dangling / weight_total, num_dangling / num_total


def measure_weight_fraction_of_soluble_material(network: Universe,
                                                rel_tol: float = 0.75, abs_tol: float = None) -> float:
    """
    Compute the weight fraction of soluble material by counting.

    Arguments:
      - network: the polymer network to do the computation for
      - rel_tol: the fraction of the maximum weight that counts as soluble. Ignored if abs_tol is specified
      - abs_tol: the weight from which on a component is not soluble anymore

    Returns:
      - :math:`W_{sol}` (float): the weight fraction of soluble material as counted.
            0. for an empty network

    """
    if (network.getNrOfAtoms() == 0):
        return 0.

    fractions = network.getClusters()
    weights = np.array([f.computeTotalMass() for f in fractions])
    total_weight = weights.sum()
    soluble_weight = 0
    for w in weights:
        if (abs_tol is not None):
            if (w < abs_tol):
                soluble_weight += w
        else:
            if (w < rel_tol * weights.max()):
                soluble_weight += w

    return soluble_weight / total_weight


def measure_lower_bound_weight_fraction_of_soluble_material(network: Universe, crossLinker_type: int = 2,
                                                            rel_tol: float = 0.75, abs_tol: float = None) -> float:
    """
    Compute a lower bound on the weight fraction of soluble material by counting.
    This works as:
        - only clusters, which do not contain loops and are smaller than the rel_tol of the biggest,
            are counted as soluble

    Arguments:
      - network: the polymer network to do the computation for
      - crossLinker_type: the type of the junctions/crossLinkers to select them in the network
      - rel_tol: the fraction of the maximum weight that counts as soluble. Ignored if abs_tol is specified
      - abs_tol: the weight from which on a component is not soluble anymore

    Returns:
      - :math:`W_{sol}` (float): the weight fraction of soluble material as counted.
            0. for an empty network

    """
    if (network.getNrOfAtoms() == 0):
        return 0.0

    def is_soluble_cluster(cluster):
        chains = cluster.getChainsWithCrosslinker(crossLinker_type)
        if (np.any([c.getType() == MoleculeType.PRIMARY_LOOP for c in chains])):
            return False
        loops = cluster.findLoops(crossLinker_type)
        return len(loops) == 0

    fractions = network.getClusters()
    weights = np.array([f.computeTotalMass() for f in fractions])
    total_weight = weights.sum()
    soluble_weight = 0
    for i in range(len(fractions)):
        w = weights[i]
        if (abs_tol is not None):
            if (w < abs_tol and is_soluble_cluster(fractions[i])):
                soluble_weight += w
        else:
            if (w < rel_tol * weights.max() and is_soluble_cluster(fractions[i])):
                soluble_weight += w

    return soluble_weight / total_weight
