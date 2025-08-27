"""
This module provides functions to compute various quantities related to polymer networks from their structure.
"""

from __future__ import annotations

import math
import statistics
import warnings
from collections import Counter
from typing import TYPE_CHECKING, Sequence, Tuple, Union

from pylimer_tools_cpp import MoleculeType

if TYPE_CHECKING:
    from pylimer_tools_cpp import Universe


def compute_stoichiometric_imbalance(
    network: Universe,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
    ignore_types: Union[list, None] = None,
    effective: bool = False,
) -> float:
    """
    Compute the stoichiometric imbalance
    ( nr. of bonds formable of crosslinker / (nr. of bonds formable of precursor chains) )

    r > 1 means an excess of crosslinkers, whereas r = 0 implies that there are not crosslinkers in the network.

    .. note::
        - if your system has a non-integer number of possible bonds (e.g. one site non-bonded),
          this will not be rounded/respected in any way.
        - Currently, only g = 2 chains are respected yet, and only one type of crosslinker is supported.
        - If you have e.g. solvent chains in the network, use `ignore_types` to not count them.

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param functionality_per_type: A dictionary with key: type, and value: functionality of this atom type.
        If `None`: will use max functionality per type.
    :param ignore_types: A list of integers, the types to ignore for the imbalance (e.g. solvent atom types)
    :param effective: Whether to use the effective functionality (if functionality_per_type is not passed)
        or the maximum functionality of the crosslinkers.
    :return: The stoichiometric imbalance. If the network is empty, or no crosslinkers are present 0. is returned.
    :rtype: float
    """
    if network.get_nr_of_atoms() == 0:
        return 0.0

    counts = Counter(network.get_atom_types())

    if functionality_per_type is not None and (
        crosslinker_type not in functionality_per_type
        or functionality_per_type[crosslinker_type] is None
    ):
        functionality_per_type = None
        warnings.warn(
            "Crosslink's atom type not found in functionality_per_type, "
            + "will ignore passed argument `functionality_per_type`."
        )

    if functionality_per_type is None:
        functionality_per_type = (
            network.determine_effective_functionality_per_type()
            if effective
            else network.determine_functionality_per_type()
        )

    if crosslinker_type not in counts:
        return 0.0

    # TODO: use the data from the functionality_per_type to determine the
    # functionality per strand, maybe?
    strands = network.get_molecules(crosslinker_type)
    if ignore_types is None:
        ignore_types = []
    ignore_types.append(crosslinker_type)
    num_relevant_strands = len(
        [
            m
            for m in strands
            if not all([a.get_type() in ignore_types for a in m.get_atoms()])
        ]
    )

    crosslinker_formable_bonds = (
        counts[crosslinker_type] * functionality_per_type[crosslinker_type]
    )
    other_formable_bonds = num_relevant_strands * 2

    if other_formable_bonds == 0:
        return math.inf

    # division by 2 is implicit
    return crosslinker_formable_bonds / (other_formable_bonds)


def compute_extent_of_reaction(
    network: Universe,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
) -> float:
    """
    Compute the extent of polymerization reaction
    (nr. of formed bonds in reaction / max. nr. of bonds formable)

    .. note::
        - if your system has a non-integer number of possible bonds (e.g. one site non-bonded),
          this will not be rounded/respected in any way
        - if the system contains solvent or other molecules that should not be binding to
          crosslinkers, make sure to remove them before calling this function
        - if you have multiple crosslinker types, be sure to replace them by one type only

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The atom type of crosslinker beads
    :param functionality_per_type: A dictionary with key: type, and value: functionality of this atom type.
        If None: will use max functionality per type.
    :return: The extent of reaction
    :rtype: float
    """
    if network.get_nr_of_atoms() == 0:
        return 1

    if len(network.get_atoms_by_type(crosslinker_type)) == 0:
        return 0

    if (functionality_per_type is not None) and (
        crosslinker_type not in functionality_per_type
    ):
        functionality_per_type = None
        warnings.warn(
            "Crosslinker type {} not found in passed functionality_per_type, ".format(
                crosslinker_type
            )
            + "will ignore passed argument `functionality_per_type`."
        )

    if functionality_per_type is None:
        functionality_per_type = network.determine_functionality_per_type()

    num_strands = len(network.get_molecules(crosslinker_type))
    crosslinks = network.get_atoms_by_type(crosslinker_type)
    num_crosslinkers = len(crosslinks)

    # assuming strand has functionality 2
    max_formable_bonds = min(
        num_strands * 2, num_crosslinkers *
        functionality_per_type[crosslinker_type]
    )

    if max_formable_bonds == 0:
        return 1

    actually_formed_bonds = 0
    for crosslink in crosslinks:
        connected_to = network.get_atoms_connected_to(atom=crosslink)
        actually_formed_bonds += len(
            [a for a in connected_to if a.get_type() != crosslinker_type]
        )

    return actually_formed_bonds / (max_formable_bonds)


def compute_fraction_of_bifunctional_reactive_sites(
    network: Universe,
    crosslinker_type: int = 2,
    functionality_per_type: Union[dict, None] = None,
) -> float:
    """
    Compute the mole fraction of reactive sites in B2
    among all reactive sites in a mixture of B1 and B2.
    Uses the network to detect what might be monofunctional chains,
    and then counts them and the bifunctional ones.

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The atom type of crosslinker beads
    :param functionality_per_type: A dictionary with key: type, and value: functionality of this atom type
    :return: The mole fraction of reactive sites in B2 among all reactive sites in a mixture of B1 and B2
    :rtype: float
    """
    if functionality_per_type is None:
        functionality_per_type = network.determine_functionality_per_type()

    monofunctional_types = [
        t for t, f in functionality_per_type.items() if f == 1]
    if len(monofunctional_types) == 0:
        """
        Assume the whole monofunctional chain has the same atom type
        """
        chains_with_crosslinks = network.get_chains_with_crosslinker(
            crosslinker_type=crosslinker_type
        )
        all_atom_types = set(network.get_atom_types())
        atom_type_has_only_monofunctional = {
            atype: True for atype in all_atom_types}
        for chain in chains_with_crosslinks:
            n_xlinks = len(chain.get_atoms_by_type(crosslinker_type))
            atom_types_in_chain = set(chain.get_atom_types())
            if n_xlinks > 1:
                for atype in atom_types_in_chain:
                    atom_type_has_only_monofunctional[atype] = False

        # the atom_type_has_only_monofunctional is now a dictionary that indicates
        # the atom types belonging to a monofunctional chain
        monofunctional_types = [
            t for t, i in atom_type_has_only_monofunctional.items() if i
        ]

    if len(monofunctional_types) == 0:
        return 1.0

    # then, count the number of chains where one end is of the monofunctional
    # type
    chains_with_crosslinks = network.get_chains_with_crosslinker(
        crosslinker_type=crosslinker_type
    )

    n_monofunctional_chains = sum(
        1
        for chain in chains_with_crosslinks
        if any(
            a.get_type() in monofunctional_types
            for a in chain.get_strand_ends(
                crosslinker_type=crosslinker_type, close_loop=True
            )
        )
        and chain.get_nr_of_atoms() > 1
    )
    n_bifunctional_chains = (
        len(network.get_molecules(crosslinker_type)) - n_monofunctional_chains
    )
    return (2 * n_bifunctional_chains) / (
        n_monofunctional_chains + 2 * n_bifunctional_chains
    )


def compute_mean_end_to_end_distances(
    networks: Sequence[Universe], crosslinker_type: int = 2
) -> dict:
    """
    Compute the mean end to end distance between each pair of (indirectly) connected crosslinker

    :param networks: The different configurations of the polymer network to do the computation for
    :type networks: Sequence[Universe]
    :param crosslinker_type: The atom type to compute the in-between vectors for
    :type crosslinker_type: int

    :return: a dictionary with key: "{atom1.name}+{atom2.name}" and value: The norm of the mean difference vector
    :rtype: dict
    """
    r_tau_vectors = compute_mean_end_to_end_vectors(networks, crosslinker_type)
    if len(r_tau_vectors) < 1:
        return {}

    def vector_norm(vector):
        return math.sqrt(sum(x * x for x in vector))

    r_taus = [vector_norm(vector) for vector in r_tau_vectors.values()]

    return dict(zip(r_tau_vectors.keys(), r_taus))


def compute_mean_end_to_end_vectors(
    networks: Sequence[Universe], crosslinker_type: int = 2
) -> dict:
    """
    Compute the mean end to end vectors between each pair of (indirectly) connected crosslinker

    :param networks: The different configurations of the polymer network to do the computation for
    :type networks: Sequence[Universe]
    :param crosslinker_type: The atom type to compute the in-between vectors for
    :type crosslinker_type: int

    :return: a dictionary with key: "{atom1.name}+{atom2.name}" and value: Their mean distance difference vector
    :rtype: dict
    """
    if len(networks) == 0:
        return {}
    end_to_end_vectors = {}
    key_counts = {}
    divider = 1 / len(networks)
    iteration = 0
    for network in networks:
        current_end_to_end_vectors = compute_end_to_end_vectors(
            network, crosslinker_type
        )
        # the mean calculation in this for loop
        # trades some memory for performance
        # there are still many performance and memory
        # improvements possible
        # (e.g. computing connectivity only once, storing it only once, ....)
        for key in current_end_to_end_vectors:
            if key not in end_to_end_vectors:
                if iteration > 0:
                    raise ValueError(
                        "Found molecule {} in network {}, but not in previous".format(
                            key, iteration
                        )
                    )
                end_to_end_vectors[key] = [0, 0, 0]
                key_counts[key] = 0
            for i in range(3):
                end_to_end_vectors[key][i] += (
                    current_end_to_end_vectors[key][i] * divider
                )
            key_counts[key] += 1
        iteration += 1
    if len(key_counts) > 0 and not all(
        c == list(key_counts.values())[0] for c in key_counts.values()
    ):
        raise ValueError("The networks contain different molecules.")
    return end_to_end_vectors


def compute_end_to_end_vectors(
        network: Universe, crosslinker_type: int = 2) -> dict:
    """
    Compute the end to end vectors between each pair of (indirectly) connected crosslinker

    :param network: The polymer network to do the computation for
    :type network: Universe
    :param crosslinker_type: The atom type to compute the in-between vectors for
    :type crosslinker_type: int

    :return: a dictionary with key: "{atom1.name}+{atom2.name}" and value: Their difference vector
    :rtype: dict
    """
    # while we could do the decomposition again with explicit removal of irrelevant strand atoms,
    # this should not be any more expensive
    end_to_end_vectors = {}
    molecules = network.get_chains_with_crosslinker(crosslinker_type)
    for molecule in molecules:
        crosslinkers = molecule.get_atoms_by_type(crosslinker_type)
        if (
            len(crosslinkers) != 2
            or molecule.get_strand_type() == MoleculeType.PRIMARY_LOOP
            or molecule.get_strand_type() == MoleculeType.DANGLING_CHAIN
        ):
            # dangling, free chains and loops are irrelevant for our purposes
            continue
        # igraph.VertexSeq is not sortable -> use a list
        crosslinkers = [crosslinkers[0], crosslinkers[1]]
        # sort crosslinkers by name as a way to keep the vector directions
        # consistent between timesteps
        crosslinkers.sort(key=lambda a: a.get_id())
        #
        end_to_end_vectors[molecule.get_key()] = crosslinkers[0].compute_vector_to(
            crosslinkers[1], network.get_box()
        )

    return end_to_end_vectors


def compute_crosslinker_conversion(
    network: Universe,
    crosslinker_type: int = 2,
    f: Union[int, None] = None,
    functionality_per_type: Union[dict, None] = None,
) -> float:
    """
    Compute the extent of reaction of the crosslinkers
    (actual functionality divided by target functionality)

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param f: The (target) functionality of the crosslinkers
    :param functionality_per_type: A dictionary with key: type, and value: (target) functionality of this atom type
    :return: The (mean) crosslinker conversion
    :rtype: float
    """
    if f is None:
        if functionality_per_type is None:
            functionality_per_type = network.determine_functionality_per_type()
        if crosslinker_type not in functionality_per_type:
            return 0.0
        f = functionality_per_type[crosslinker_type]

    if f is None or f <= 0.0 or not math.isfinite(f):
        raise ValueError(
            "Crosslinker functionality = {} is not reasonable.".format(f))

    return compute_effective_crosslinker_functionality(
        network, crosslinker_type) / f


def compute_effective_crosslinker_functionality(
    network: Universe, crosslinker_type: int = 2
) -> float:
    """
    Compute the mean crosslinker functionality

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :return: The (mean) effective crosslinker functionality
    :rtype: float
    """
    junction_degrees = compute_effective_crosslinker_functionalities(
        network, crosslinker_type
    )
    return statistics.mean(junction_degrees) if len(
        junction_degrees) > 0 else 0.0


def compute_effective_crosslinker_functionalities(
    network: Universe, crosslinker_type: int = 2
) -> list[int]:
    """
    Compute the functionality of every crosslinker in the network

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :return: The functionality of every crosslinker
    :rtype: list[int]
    """
    if network.get_nr_of_atoms() == 0:
        return []
    junctions = network.get_atoms_by_type(crosslinker_type)
    junction_ids = [v.get_id() for v in junctions]
    junction_degrees = [
        network.get_nr_of_bonds_of_atom(id) for id in junction_ids]
    return junction_degrees


def compute_weight_fractions(network: Universe) -> dict:
    """
    Compute the weight fractions of each atom type in the network.
    Kept for compatibility reasons; just call
    :func:`pylimer_tools_cpp.Universe.compute_weight_fractions()` instead.

    :param network: The polymer network to do the computation for
    :return: Using the type i as a key, this dict contains the weight fractions (:math:`\\frac{W_i}{W_{tot}}`)
    :rtype: dict
    """
    return network.compute_weight_fractions()


def measure_weight_fraction_of_backbone(
        network: Universe, crosslinker_type: int = 2):
    """
    Compute the weight fraction of network backbone in infinite network

    :param network: The network to compute the weight fraction for
    :param crosslinker_type: The atom type to use to split the molecules
    :return: 1 - weightDangling/weightTotal - weightSoluble/weightTotal
    :rtype: float

    See also:
      - :func:`pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_dangling_chains()`
      - :func:`pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_soluble_material()`
    """
    if network.get_nr_of_atoms() < 1:
        return 0.0

    weight_fraction_dangling, _ = measure_weight_fraction_of_dangling_chains(
        network, crosslinker_type
    )

    weight_fraction_soluble = measure_weight_fraction_of_soluble_material(
        network)

    return 1.0 - weight_fraction_dangling - weight_fraction_soluble


def measure_weight_fraction_of_dangling_chains(
    network: Universe, crosslinker_type: int = 2
) -> Tuple[float, float]:
    """
    Compute the weight fraction of dangling strands in infinite network

    .. note::
        Currently, only primary dangling chains are taken into account.
        There are other methods that incorporate more.

    :param network: The network to compute the weight fraction for
    :param crosslinker_type: The atom type to use to split the molecules
    :return: A tuple containing (weightDangling/weightTotal, numDangling/numTotal)
    :rtype: Tuple[float, float]
    """
    if network.get_nr_of_atoms() < 1:
        return 0.0, 0.0

    weights = network.get_masses()

    def get_weight_of_graph(graph):
        counts = Counter(graph.get_atom_types())
        weight_total = 0
        for key in counts:
            weight_total += weights[key] * counts[key]
        return weight_total

    all_chains = network.get_chains_with_crosslinker(crosslinker_type)
    num_total = network.get_nr_of_atoms()
    weight_total = get_weight_of_graph(network)

    num_dangling = 0
    weight_dangling = 0
    for chain in all_chains:
        if chain.get_strand_type() == MoleculeType.DANGLING_CHAIN:
            num_dangling += chain.get_nr_of_atoms()
            weight_dangling += get_weight_of_graph(chain)

    if weight_total == 0:
        # warnings.warn("Total weight of network is = 0.")
        return 0.0, num_dangling / num_total

    return weight_dangling / weight_total, num_dangling / num_total


def measure_weight_fraction_of_soluble_material(
    network: Universe, rel_tol: float = 0.5, abs_tol: Union[float, None] = None
) -> float:
    """
    Compute the weight fraction of soluble material by counting.
    Effectively, this method counts the weight of clusters
    that have a weight less than a certain fraction of the total weight.

    :param network: The polymer network to do the computation for
    :param rel_tol: The fraction of the maximum weight that counts as soluble. Ignored if abs_tol is specified
    :param abs_tol: The weight from which on a component is not soluble anymore
    :return: The weight fraction of soluble material as counted. 0 for an empty network
    :rtype: float
    """
    if network.get_nr_of_atoms() == 0:
        return 0.0

    fractions = network.get_clusters()
    weights = [f.compute_total_mass() for f in fractions]
    total_weight = sum(weights)
    soluble_weight = 0
    max_weight = max(weights) if weights else 0
    for w in weights:
        if abs_tol is not None:
            if w < abs_tol:
                soluble_weight += w
        else:
            if w < rel_tol * max_weight:
                soluble_weight += w

    if total_weight == 0:
        return 0.0

    return soluble_weight / total_weight


def measure_lower_bound_weight_fraction_of_soluble_material(
    network: Universe,
    crosslinker_type: int = 2,
    rel_tol: float = 0.75,
    abs_tol: Union[float, None] = None,
) -> float:
    """
    Compute a lower bound on the weight fraction of soluble material by counting.
    This is ones as such: only clusters, which do not contain loops and are smaller than the rel_tol of the biggest,
    are counted as soluble

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param rel_tol: The fraction of the maximum weight that counts as soluble. Ignored if abs_tol is specified
    :param abs_tol: The weight from which on a component is not soluble anymore
    :return: The weight fraction of soluble material as counted. 0 for an empty network
    :rtype: float
    """
    if network.get_nr_of_atoms() == 0:
        return 0.0

    def is_soluble_cluster(cluster):
        chains = cluster.get_chains_with_crosslinker(crosslinker_type)
        if any(c.get_strand_type() == MoleculeType.PRIMARY_LOOP for c in chains):
            return False
        loops = cluster.find_loops(crosslinker_type)
        return len(loops) == 0

    fractions = network.get_clusters()
    weights = [f.compute_total_mass() for f in fractions]
    total_weight = sum(weights)
    soluble_weight = 0
    max_weight = max(weights) if weights else 0
    for i in range(len(fractions)):
        w = weights[i]
        if abs_tol is not None:
            if w < abs_tol and is_soluble_cluster(fractions[i]):
                soluble_weight += w
        else:
            if w < rel_tol * max_weight and is_soluble_cluster(fractions[i]):
                soluble_weight += w

    return soluble_weight / total_weight
