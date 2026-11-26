"""
This module provides functions to compute various quantities related to polymer networks from their structure.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import math
import warnings
from typing import TYPE_CHECKING, Tuple, Union

from pylimer_tools_cpp import MoleculeType

if TYPE_CHECKING:
    from pylimer_tools_cpp import Molecule, Universe


def determine_chain_functionality_per_type(
    network: Universe,
    crosslinker_type: int = 2,
) -> dict[int, int]:
    """
    Determine the reactive chain functionality for each atom type in the network.

    Unlike :meth:`pylimer_tools_cpp.Universe.determine_functionality_per_type`, which
    returns the maximum vertex degree per atom type (and therefore reports ``2`` for
    internal beads of multi-bead monofunctional and solvent chains), this function
    infers non-crosslinker chain functionalities from the maximum number of bonds
    to crosslinker atoms observed across chains of each atom type:

    - ``crosslinker_type``: maximum vertex degree of crosslinker atoms (:math:`f`)
    - bifunctional (:math:`B_2`) chain types: ``2`` (up to 2 crosslinker bonds per chain)
    - monofunctional (:math:`B_1`) chain types: ``1`` (at most 1 crosslinker bond per chain)
    - solvent (:math:`B_0`) chain types: ``0`` (no crosslinker bonds)

    If no bonds between crosslinkers and non-crosslinker chains exist anywhere in the
    network (completely unreacted system), non-crosslinker atom types default to ``2``.

    .. note::
        - Each chain species (:math:`B_2`, :math:`B_1`, :math:`B_0`) is assumed to have
          distinct atom type(s) representing the chain's reactive functionality (rather
          than sharing interior backbone atom types across species or using distinct
          per-bead end-cap types).
        - Because functionalities are inferred from *observed* attachments in ``network``,
          partially reacted / low-conversion networks in which no chain or crosslinker of a
          given type has yet formed its maximum number of bonds will underestimate that
          type's target functionality; pass ``functionality_per_type`` explicitly in those
          cases.

    :param network: The polymer network to analyze
    :param crosslinker_type: The atom type of crosslinker beads
    :return: Dictionary mapping atom type to its chain/crosslinker functionality
    :rtype: dict[int, int]
    """
    bead_func: dict[int, int] = network.determine_functionality_per_type()
    if network.get_nr_of_atoms() == 0:
        return bead_func

    chains = network.get_chains_with_crosslinker(crosslinker_type=crosslinker_type)
    max_xlink_bonds_per_type: dict[int, int] = {}
    total_xlink_bonds = 0

    for chain in chains:
        xlink_atoms = chain.get_atoms_by_type(crosslinker_type)
        n_xlink_bonds = sum(
            1
            for a in xlink_atoms
            for nbr in chain.get_atoms_connected_to(a)
            if nbr.get_type() != crosslinker_type
        )
        total_xlink_bonds += n_xlink_bonds
        for atype in set(chain.get_atom_types()):
            if atype != crosslinker_type:
                if (
                    atype not in max_xlink_bonds_per_type
                    or n_xlink_bonds > max_xlink_bonds_per_type[atype]
                ):
                    max_xlink_bonds_per_type[atype] = n_xlink_bonds

    result: dict[int, int] = {}
    for atype, deg in bead_func.items():
        if atype == crosslinker_type:
            result[atype] = deg
        elif total_xlink_bonds == 0:
            result[atype] = 2
        else:
            result[atype] = max_xlink_bonds_per_type.get(atype, 0)

    return result


def _resolve_chain_functionality_per_type(
    network: Universe,
    crosslinker_type: int,
    functionality_per_type: Mapping[int, int | None] | None,
    warn_missing_crosslinker: bool = False,
) -> dict[int, int]:
    """
    Resolve a complete chain functionality map for all atom types in ``network``.

    Explicit entries in ``functionality_per_type`` take precedence; any missing
    atom types are inferred via :func:`determine_chain_functionality_per_type`.
    When ``functionality_per_type`` already covers every atom type present in
    ``network``, chain decomposition is skipped entirely.
    """
    if functionality_per_type is not None:
        if warn_missing_crosslinker and (
            crosslinker_type not in functionality_per_type
            or functionality_per_type[crosslinker_type] is None
        ):
            warnings.warn(
                f"Crosslinker type {crosslinker_type} not found in passed "
                "functionality_per_type; inferring its functionality from the network."
            )
        atom_types_in_net = network.count_atom_types()
        if all(
            t in functionality_per_type and functionality_per_type[t] is not None
            for t in atom_types_in_net
        ):
            return {
                int(atype): int(func)
                for atype, func in functionality_per_type.items()
                if func is not None
            }

    inferred = determine_chain_functionality_per_type(
        network, crosslinker_type=crosslinker_type
    )
    if functionality_per_type is None:
        return inferred
    resolved = dict(inferred)
    for atype, func in functionality_per_type.items():
        if func is not None:
            resolved[int(atype)] = int(func)
    return resolved


def _get_strand_functionality(
    strand: Molecule,
    functionality_per_type: Mapping[int, int],
    ignore_types_set: set[int],
) -> int:
    """
    Determine the number of reactive sites (functionality) of a single strand.

    Strands composed solely of ignored atom types contribute 0.
    """
    atom_types = set(strand.get_atom_types())
    if not atom_types or atom_types.issubset(ignore_types_set):
        return 0
    relevant_types = atom_types - ignore_types_set
    funcs = [
        int(functionality_per_type[t])
        for t in relevant_types
        if t in functionality_per_type and functionality_per_type[t] is not None
    ]
    if not funcs:
        return 2
    return min(funcs)


def count_strands_by_functionality(
    network: Universe,
    crosslinker_type: int = 2,
    functionality_per_type: Mapping[int, int | None] | None = None,
    ignore_types: Sequence[int] | None = None,
) -> dict[int, int]:
    """
    Count precursor strands (after splitting at ``crosslinker_type``) grouped by
    their reactive functionality (e.g. ``0`` for solvent, ``1`` for :math:`B_1`,
    ``2`` for :math:`B_2`).

    :param network: The polymer network to analyze
    :param crosslinker_type: The atom type of crosslinker beads
    :param functionality_per_type: Optional mapping from atom type to functionality
    :param ignore_types: Optional sequence of atom types to treat as non-reactive (functionality 0)
    :return: Dictionary mapping strand functionality to number of strands
    :rtype: dict[int, int]
    """
    if network.get_nr_of_atoms() == 0:
        return {}
    chain_func_per_type = _resolve_chain_functionality_per_type(
        network, crosslinker_type, functionality_per_type
    )
    ignore_types_set = set(ignore_types) if ignore_types is not None else set()
    ignore_types_set.discard(crosslinker_type)
    counts_by_func: dict[int, int] = {}
    for strand in network.get_molecules(crosslinker_type):
        func = _get_strand_functionality(strand, chain_func_per_type, ignore_types_set)
        counts_by_func[func] = counts_by_func.get(func, 0) + 1
    return counts_by_func


def compute_stoichiometric_imbalance(
    network: Universe,
    crosslinker_type: int = 2,
    functionality_per_type: Mapping[int, int | None] | None = None,
    ignore_types: Sequence[int] | None = None,
    effective: bool = False,
) -> float:
    """
    Compute the stoichiometric imbalance
    ( nr. of bonds formable of crosslinker / (nr. of bonds formable of precursor chains) )

    r > 1 means an excess of crosslinkers, whereas r = 0 implies that there are not crosslinkers in the network.

    .. note::
        - if your system has a non-integer number of possible bonds (e.g. one site non-bonded),
          this will not be rounded/respected in any way.
        - Only one type of crosslinker is supported.
        - Monofunctional (:math:`B_1`) and solvent (:math:`B_0`) chains are accounted for
          using ``functionality_per_type`` or inferred from crosslinker attachments;
          you can also use ``ignore_types`` to explicitly exclude atom types.

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The type of the junctions/crosslinkers to select them in the network
    :param functionality_per_type: A dictionary with key: type, and value: functionality of this atom type.
        If `None`: will infer chain functionalities from the network.
    :param ignore_types: A sequence of integers, the types to ignore for the imbalance (e.g. solvent atom types)
    :param effective: Whether to use the effective functionality (if crosslinker functionality is not passed)
        or the maximum functionality of the crosslinkers.
    :return: The stoichiometric imbalance. If the network is empty, or no crosslinkers are present 0. is returned.
    :rtype: float
    """
    if network.get_nr_of_atoms() == 0:
        return 0.0

    counts = network.count_atom_types()
    crosslinker_explicit = (
        functionality_per_type is not None
        and crosslinker_type in functionality_per_type
        and functionality_per_type[crosslinker_type] is not None
    )

    chain_func_per_type = _resolve_chain_functionality_per_type(
        network,
        crosslinker_type,
        functionality_per_type,
        warn_missing_crosslinker=True,
    )

    if crosslinker_type not in counts:
        return 0.0

    if not crosslinker_explicit and effective:
        crosslinker_func: float = (
            network.determine_effective_functionality_per_type().get(
                crosslinker_type, 0.0
            )
        )
    else:
        crosslinker_func = float(chain_func_per_type[crosslinker_type])

    strand_counts = count_strands_by_functionality(
        network, crosslinker_type, chain_func_per_type, ignore_types
    )
    crosslinker_formable_bonds = counts[crosslinker_type] * crosslinker_func
    other_formable_bonds = sum(
        func * cnt for func, cnt in strand_counts.items() if func > 0
    )

    if other_formable_bonds == 0:
        return math.inf

    return crosslinker_formable_bonds / other_formable_bonds


def compute_extent_of_reaction(
    network: Universe,
    crosslinker_type: int = 2,
    functionality_per_type: Mapping[int, int | None] | None = None,
    ignore_types: Sequence[int] | None = None,
) -> float:
    """
    Compute the extent of polymerization reaction
    (nr. of formed bonds in reaction / max. nr. of bonds formable)

    .. note::
        - if your system has a non-integer number of possible bonds (e.g. one site non-bonded),
          this will not be rounded/respected in any way
        - Monofunctional (:math:`B_1`) and solvent (:math:`B_0`) chains are accounted for
          using ``functionality_per_type`` or inferred from crosslinker attachments;
          you can also use ``ignore_types`` to explicitly exclude atom types.
        - if you have multiple crosslinker types, be sure to replace them by one type only

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The atom type of crosslinker beads
    :param functionality_per_type: A dictionary with key: type, and value: functionality of this atom type.
        If None: will infer chain functionalities from the network.
    :param ignore_types: A sequence of integers, the types to ignore (e.g. solvent atom types)
    :return: The extent of reaction
    :rtype: float
    """
    if network.get_nr_of_atoms() == 0:
        return 1.0

    crosslinks = network.get_atoms_by_type(crosslinker_type)
    num_crosslinkers = len(crosslinks)
    if num_crosslinkers == 0:
        return 0.0

    crosslinker_explicit = (
        functionality_per_type is not None
        and crosslinker_type in functionality_per_type
        and functionality_per_type[crosslinker_type] is not None
    )
    chain_func_per_type = _resolve_chain_functionality_per_type(
        network,
        crosslinker_type,
        functionality_per_type,
        warn_missing_crosslinker=True,
    )
    ignore_types_set = set(ignore_types) if ignore_types is not None else set()
    ignore_types_set.discard(crosslinker_type)

    strand_counts = count_strands_by_functionality(
        network, crosslinker_type, chain_func_per_type, ignore_types
    )
    other_formable_bonds = sum(
        func * cnt for func, cnt in strand_counts.items() if func > 0
    )

    actually_formed_bonds = 0
    for crosslink in crosslinks:
        connected_to = network.get_atoms_connected_to(atom=crosslink)
        actually_formed_bonds += sum(
            1
            for a in connected_to
            if a.get_type() != crosslinker_type and a.get_type() not in ignore_types_set
        )

    crosslinker_formable_bonds = (
        num_crosslinkers * chain_func_per_type[crosslinker_type]
    )
    max_formable_bonds = min(other_formable_bonds, crosslinker_formable_bonds)

    if max_formable_bonds == 0:
        if (
            not crosslinker_explicit
            and chain_func_per_type[crosslinker_type] == 0
            and other_formable_bonds > 0
            and actually_formed_bonds == 0
        ):
            return 0.0
        return 1.0

    return actually_formed_bonds / max_formable_bonds


def compute_fraction_of_bifunctional_reactive_sites(
    network: Universe,
    crosslinker_type: int = 2,
    functionality_per_type: Mapping[int, int | None] | None = None,
) -> float:
    """
    Compute the mole fraction of reactive sites in B2
    among all reactive sites in a mixture of B1 and B2.
    Uses the network to detect monofunctional and solvent chains,
    and then counts the monofunctional and bifunctional ones.

    :param network: The polymer network to do the computation for
    :param crosslinker_type: The atom type of crosslinker beads
    :param functionality_per_type: A dictionary with key: type, and value: functionality of this atom type
    :return: The mole fraction of reactive sites in B2 among all reactive sites in a mixture of B1 and B2
    :rtype: float
    """
    if network.get_nr_of_atoms() == 0:
        return 1.0

    strand_counts = count_strands_by_functionality(
        network, crosslinker_type, functionality_per_type
    )
    n_monofunctional_chains = strand_counts.get(1, 0)
    n_bifunctional_chains = strand_counts.get(2, 0)

    total_reactive_sites = n_monofunctional_chains + 2 * n_bifunctional_chains
    if total_reactive_sites == 0:
        return 1.0
    return (2 * n_bifunctional_chains) / total_reactive_sites


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
    if len(key_counts) > 0 and not all(c == len(networks) for c in key_counts.values()):
        raise ValueError("The networks contain different molecules.")
    return end_to_end_vectors


def compute_end_to_end_vectors(network: Universe, crosslinker_type: int = 2) -> dict:
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
    f: int | None = None,
    functionality_per_type: Mapping[int, int | None] | None = None,
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
        inferred_from_network = functionality_per_type is None
        if functionality_per_type is None:
            functionality_per_type = network.determine_functionality_per_type()
        if (
            crosslinker_type not in functionality_per_type
            or functionality_per_type[crosslinker_type] is None
        ):
            return 0.0
        f = functionality_per_type[crosslinker_type]
        if (
            inferred_from_network
            and f == 0
            and compute_effective_crosslinker_functionality(network, crosslinker_type)
            == 0.0
        ):
            return 0.0

    if f is None or f <= 0.0 or not math.isfinite(f):
        raise ValueError("Crosslinker functionality = {} is not reasonable.".format(f))

    return compute_effective_crosslinker_functionality(network, crosslinker_type) / f


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
    if network.get_nr_of_atoms() == 0:
        return 0.0
    return network.determine_effective_functionality_per_type().get(
        crosslinker_type, 0.0
    )


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
    junction_degrees = [network.get_nr_of_bonds_of_atom(id) for id in junction_ids]
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
    network: Universe, crosslinker_type: int = 2
) -> float:
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

    weight_fraction_soluble = measure_weight_fraction_of_soluble_material(network)

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
    has_all_masses = all(t in weights for t in network.count_atom_types())

    def get_weight_of_graph(graph: Union[Universe, Molecule]) -> float:
        if has_all_masses:
            return graph.compute_total_mass()
        return sum(weights.get(t, 1.0) for t in graph.get_atom_types())

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
