# source: https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

from __future__ import annotations

import warnings
from collections import Counter
from typing import Iterable, Tuple

import igraph
import numpy as np
from pylimer_tools.entities.molecule import Molecule
from pylimer_tools.entities.universum import Universum


def predictShearModulus(networks: Iterable[Universum], T: float = 1, k_B: float = 1, foreignAtomType=None, totalMass=1) -> float:
    """
    Predict the shear modulus using ANT Analysis.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the polymer system to predict the shear modulus for
      - T: the temperature in your unit system
      - k_b: Boltzmann's constant in your unit system
      - foreignAtomType: the type of atoms to ignore (junctions, crosslinkers)
      - totalMass: the $M$ in the respective formula
    
    Returns:
      - shear modulus (float): the estimated shear modulus. Unit: [pressure]
    """
    Gamma = calculateTopologicalFactor(networks, foreignAtomType, totalMass)
    nu = 0
    for network in networks:
      nu += len(network.getMolecules(foreignAtomType)) / (network.getVolume()) / len(networks)
    return Gamma*nu*k_B*T


def calculateCycleRank(networks: Iterable[Universum] = None, nu: int = None, mu: int = None, absTol: float = 1, relTol: float = 1, junctionType=None) -> float:
    """
    Compute the cycle rank ($\\chi$).
    Assumes the precursor-chains to be bifunctional.

    Arguments:
      - network: the network to calculate the cycle rank for
      - nu: number of elastically effective (active) strands per unit volume
      - mu: number density of the elastically effective crosslink
      - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol)
      - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol))
      - junctionType: the atom type of the crosslinkers/junctions

    No need to provide all the parameters — either/or:
    - nu & mu
    - network, absTol, relTol, junctionType

    Returns:
      - cycleRank: the cycle rank ($\\xi = \\nu_{eff} - \\mu_{eff}$). Unit: [1/Volume]
    """
    if (nu is None):
        if (junctionType is None or networks is None):
            raise ValueError(
                "Argument missing: When not specifiying nu, network and junctionType need to be specified")
        nu = calculateEffectiveNrDensityOfNetwork(
            networks, absTol, relTol, junctionType)
    if (mu is None):
        if (junctionType is None or networks is None):
            raise ValueError(
                "Argument missing: When not specifiying mu, network and junctionType need to be specified")
        mu = calculateEffectiveNrDensityOfJunctions(
            networks, absTol, relTol, junctionType)

    return nu - mu


def calculateEffectiveNrDensityOfNetwork(networks: Iterable[Universum], absTol: float = 1, relTol: float = 1, junctionType=None) -> float:
    """
    Compute the effective number density $\\nu_{eff}$ of a network.
    Assumes the precursor-chains to be bifunctional.

    $\\nu_{eff}$ is the number of elastically effective (active) strands per unit volume,
    which are defined as the ones that can store elastic energy
    upon network deformation, resp. the effective number density of network strands

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network (pylimer_tools.entities.Universum): the network to compute $\\nu_{eff}$ for
      - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol)
      - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol))
      - junctionType: the atom type of the crosslinkers/junctions

    Returns:
      - $\\nu_{eff}$ (float): the effective number density of network strands. Unit: [1/Volume]
    """
    if (len(networks) == 0):
        return None

    # get the mean end to end distances
    R_taus = computeMeanEndToEndDistances(networks, junctionType)
    if (len(R_taus) < 1):
        return 0.0
    R_taus = np.array(list(R_taus.values()))
    R_tau_max = np.max(R_taus)

    # process additional input parameters
    if (absTol is None):
        absTol = R_tau_max

    # count how many effective strands there are
    numEffective = np.array([R_tau > absTol or R_tau > relTol*R_tau_max
                             for R_tau in R_taus]).sum()
    meanVolume = calculateMeanUniverseVolume(networks)

    return numEffective / meanVolume


def calculateMeanUniverseVolume(networks: Iterable[Universum], acceptDifferentSizes: bool = False) -> float:
    """
    Compute the mean volume of a list of universes.

    Arguments:
      - networks: a list of universes
      - acceptDifferentSizes: toggle whether to throw an error when the Universe have different nr. of atoms

    Returns:
      - meanVolume (float): the mean volume of the universes
    """
    # compute the mean volume of the universes
    meanVolume = 0
    divisor = 1/len(networks)
    networkSize = networks[0].getSize()
    for network in networks:
        if (not acceptDifferentSizes and network.getSize() != networkSize):
            raise NotImplementedError(
                "Currently, only sequences of networks with the same size are supported (got one with {} instead of {})".format(network.getSize(), networkSize))
        meanVolume += network.getVolume()*divisor
    return meanVolume


def calculateEffectiveNrDensityOfJunctions(networks: Iterable[Universum], absTol: float = 0, relTol: float = 1, junctionType=None, minNumEffectiveStrands=2) -> float:
    """
    Compute the number density of the elastically effective crosslinks,
    defined as the ones that connect at least two elastically effective strands.
    Assumes the precursor-chains to be bifunctional.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network (pylimer_tools.entities.Universum): the network to compute $\\nu_{eff}$ for
      - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol)
      - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol))
      - junctionType: the atom type of the crosslinkers/junctions
      - minNumEffectiveStrands (int): the number of elastically effective strands to qualify a junction as such

    Returns:
      - $\\mu_{eff}$ (float): the effective number density of junctions. Unit: [1/Volume]
    """
    if (len(networks) < 1):
        return None
    if (junctionType is None):
        return 0.0

    meanVolume = calculateMeanUniverseVolume(networks)

    if (minNumEffectiveStrands == 0):
        return len(networks[0].getAtomsWithType(junctionType))/meanVolume

    # get the mean end to end distances
    R_taus = computeMeanEndToEndDistances(networks, junctionType)
    if (len(R_taus) < 1):
        return 0.0
    R_tau_max = max(R_taus.values())

    # process additional input parameters
    if (absTol is None):
        absTol = R_tau_max

    # count how many active connections each junction has
    junctionActivity = {}
    for key in R_taus:
        crosslinkerNames = key.split("+")
        assert(len(crosslinkerNames) == 3)
        isActive = R_taus[key] > absTol or R_taus[key] > relTol*R_tau_max
        if (not(isActive)):
            continue
        relevantNames = [crosslinkerNames[0], crosslinkerNames[1]]
        for crosslinkerName in relevantNames:
            if (crosslinkerName not in junctionActivity):
                junctionActivity[crosslinkerName] = 0
            junctionActivity[crosslinkerName] += 1

    effectiveJunctions = np.array(
        [junctionActivity[key] >= minNumEffectiveStrands for key in junctionActivity])
    numEffectiveJunctions = effectiveJunctions.sum()
    return numEffectiveJunctions/meanVolume


def calculateWeightFractionOfBackbone(network: Universum, crosslinkerType, weights=1):
    """
    Compute the weight fraction of network backbone in infinite network

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinkerType: the atom type to use to split the molecules
      - weights: either a dict with key: atomType and value: weight, or a scalar value if all atoms have the same weight

    Returns:
      - weightFraction (float): 1 - weightDangling/weightTotal,
    """
    weightFraction, _ = calculateWeightFractionOfDanglingChains(
        network, crosslinkerType, weights)
    return 1.0 - weightFraction


def calculateWeightFractionOfDanglingChains(network: Universum, crosslinkerType, weights=1) -> Tuple[float, float]:
    """
    Compute the weight fraction of dangling strands in infinite network

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinkerType: the atom type to use to split the molecules
      - weights: either a dict with key: atomType and value: weight, or a scalar value if all atoms have the same weight

    Returns:
      - weightFraction: weightDangling/weightTotal,
      - numFraction: numDangling/numTotal
    """
    if (network.getSize() < 1):
        return 0.0, 0.0

    def getWeightOfGraph(graph):
        counts = Counter(graph.vs["type"])
        weightTotal = 0
        for key in counts:
            if (type(weights) in (float, int)):
                weightTotal += weights*counts[key]
            else:
                weightTotal += weights[key]*counts[key]
        return weightTotal

    allChains = network.getChainsWithCrosslinker(crosslinkerType)
    numTotal = network.getSize()
    weightTotal = getWeightOfGraph(network.getUnderlyingGraph())

    numDangling = 0
    weightDangling = 0
    for chain in allChains:
        if (chain.getType() is Molecule.MoleculeType.DANGLING_CHAIN):
            numDangling += chain.getLength()
            weightDangling += getWeightOfGraph(chain.getUnderlyingGraph())

    if (weightTotal == 0):
        # warnings.warn("Weight total = 0")
        return 0.0, numDangling/numTotal

    return weightDangling/weightTotal, numDangling/numTotal


def computeMeanEndToEndDistances(networks: Iterable[Universum], crosslinkerType) -> dict:
    """
    Compute the mean end to end distance between each pair of (indirectly) connected crosslinker

    Arguments:
      - networks: the different configurations of the polymer network to do the computation for
      - crosslinkerType: the atom type to compute the in-between vectors for

    Returns:
      - endToEndDistances (dict): a dictionary with key: "{atom1.name}+{atom2.name}"
          and value: the norm of the mean difference vector
    """
    R_tau_vectors = computeMeanEndToEndVectors(networks, crosslinkerType)
    if (len(R_tau_vectors) < 1):
        return {}

    R_tau_vectors_array = np.array(list(R_tau_vectors.values()))
    # print(R_tau_vectors_array)
    R_taus = np.linalg.norm(R_tau_vectors_array, axis=1)

    return dict(zip(R_tau_vectors.keys(), R_taus))


def computeMeanEndToEndVectors(networks: Iterable[Universum], crosslinkerType) -> dict:
    """
    Compute the mean end to end vectors between each pair of (indirectly) connected crosslinker

    Arguments:
      - networks: the different configurations of the polymer network to do the computation for
      - crosslinkerType: the atom type to compute the in-between vectors for

    Returns:
      - endToEndVectors (dict): a dictionary with key: "{atom1.name}+{atom2.name}"
          and value: their mean difference vector
    """
    if (len(networks) == 0):
        return {}
    endToEndVectors = {}
    divider = 1/len(networks)
    for network in networks:
        currentEndToEndVectors = computeEndToEndVectors(
            network, crosslinkerType)
        # the mean calculation in this for loop
        # trades some memory for performance
        # there are still many performance and memory
        # improvements possible
        # (e.g. computing connectivity only once, storing it only once, ....)
        for key in currentEndToEndVectors:
            if (key not in endToEndVectors):
                endToEndVectors[key] = 0
            endToEndVectors[key] += currentEndToEndVectors[key]*divider
    return endToEndVectors


def computeEndToEndVectors(network: Universum, crosslinkerType) -> dict:
    """
    Compute the end to end vectors between each pair of (indirectly) connected crosslinker

    Arguments:
      - network: the polymer network to do the computation for
      - crosslinkerType: the atom type to compute the in-between vectors for

    Returns:
      - endToEndVectors (dict): a dictionary with key: "{atom1.name}+{atom2.name}"
          and value: their difference vector
    """
    # while we could do the decomposition again with explicit removal of irrelevant strand atoms,
    # this should not be any more expensive
    endToEndVectors = {}
    molecules = network.getChainsWithCrosslinker(crosslinkerType)
    for molecule in molecules:
        crosslinkers = molecule.getUnderlyingGraph().vs.select(type_eq=crosslinkerType)
        if (len(crosslinkers) != 2 or
            molecule.getType() == Molecule.MoleculeType.PRIMARY_LOOP or
                molecule.getType() == Molecule.MoleculeType.DANGLING_CHAIN):
            # dangling, free chains and loops are irrelevant for our purposes
            continue
        # igraph.VertexSeq is not sortable -> use a list
        crosslinkers = [crosslinkers[0], crosslinkers[1]]
        # sort crosslinkers by name as a way to keep the vector directions consistent between timesteps
        crosslinkers.sort(key=lambda a: a["name"])
        #
        key = _getKeyForMolecule(molecule, crosslinkers)
        endToEndVectors[key] = crosslinkers[0]["atom"].computeVectorTo(
            crosslinkers[1]["atom"])

    return endToEndVectors


def computeCrosslinkerConversion(network: Universum, junctionType, f: int) -> float:
    """
    Compute the extent of reaction of the crosslinkers
    (actual functionality divided by target functionality)

    Arguments:
      - network: the polymer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - f: the functionality of the crosslinkers

    Returns:
      - r (float): the (mean) crosslinker conversion
    """
    return calculateEffectiveCrosslinkerFunctionality(network, junctionType) / f


def calculateEffectiveCrosslinkerFunctionality(network: Universum, junctionType) -> float:
    """
    Compute the mean crosslinker functionality

    Arguments:
      - network: the polymer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network

    Returns:
      - f (float): the (mean) effective crosslinker functionality
    """
    junctionDegrees = calculateEffectiveCrosslinkerFunctionalities(
        network, junctionType)
    return np.mean(junctionDegrees)


def calculateEffectiveCrosslinkerFunctionalities(network: Universum, junctionType) -> list[int]:
    """
    Compute the functionality of every crosslinker in the network

    Arguments:
      - network: the polymer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network

    Returns:
      - junctionDegrees (list[int]): the functionality of every crosslinker
    """
    if (network.getSize() == 0):
        return []
    junctions = network.getUnderlyingGraph().vs.select(type_eq=junctionType)
    junctionIds = [v.index for v in junctions]
    junctionDegrees = network.getUnderlyingGraph().degree(
        junctionIds, mode="all", loops=False)
    return junctionDegrees


def calculateTopologicalFactor(networks: Iterable[Universum], foreignAtomType=None, totalMass=1, b=None) -> float:
    """
    Compute the topological factor of a polymer network.

    Assumptions: 
      - the precursor-chains to be bifunctional
      - all Universes to have the same structure (with possibly differing positions)
      - crosslinkers do not count to the nr. of monomers in a strand

    Source:
      - eq. 16 in https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the network to compute the topological factor for
      - foreignAtomType: the type of atoms to ignore
      - totalMass: the $M$ in the respective formula
      - b: the mean bond length. 
          If `None`, it will be computed for each molecule in the first Universum (Network).

    Returns:
      - the topological factor $\\Gamma$
    """
    R_taus = computeMeanEndToEndDistances(networks, foreignAtomType)

    # find the topological factor
    GammaSum = 0
    network = networks[0]  # this is where the second assumption is made
    chainsToProcess = network.getChainsWithCrosslinker(foreignAtomType)
    for molecule in chainsToProcess:
        crosslinkers = molecule.getUnderlyingGraph().vs.select(type_eq=foreignAtomType)
        if (len(crosslinkers) != 2 or
            molecule.getType() == Molecule.MoleculeType.PRIMARY_LOOP or
                molecule.getType() == Molecule.MoleculeType.DANGLING_CHAIN):
            # dangling, free chains and loops are irrelevant for our purposes
            continue
        if (b is None):
            b = molecule.computeBondLengths().mean()
        crosslinkers = [crosslinkers[0], crosslinkers[1]]
        # sort crosslinkers by name as a way to keep the vector directions consistent between timesteps
        crosslinkers.sort(key=lambda a: a["name"])
        key = _getKeyForMolecule(molecule, crosslinkers)
        GammaSum += R_taus[key]*R_taus[key] / \
            ((molecule.getLength()-2) * b *
             b)  # -2: remove crosslinkers again (assumption 3)

    return GammaSum / totalMass


def _getKeyForMolecule(molecule, crosslinkers):
    """
    Get a key to identify a molecule.
    The crosslinkers (ends of the molecule) are the first two components of the key, 
    whereas the names of all other atoms in the chain are used too to distinguish 
    e.g. two secondary loops
    """
    names = [a.name for a in molecule]
    names.sort()
    return "{}+{}+{}".format(crosslinkers[0]
                             ["name"], crosslinkers[1]["name"], "".join(names))
