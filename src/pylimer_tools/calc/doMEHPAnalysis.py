# source: https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

from collections import Counter

import numpy as np
from pylimer_tools.entities.molecule import Molecule
from pylimer_tools.entities.universum import Universum


def calculateCycleRank(network: Universum, nu: int = None, mu: int = None, absTol: float = 1, relTol: float = 1, junctionType=None):
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
      - cycleRank: the cycle rank ($\\chi = \\nu_{eff} - \\mu_{eff}$) 
    """
    if (nu is None):
        nu = calculateEffectiveNrDensityOfNetwork(
            network, absTol, relTol, junctionType)
    if (mu is None):
        mu = calculateEffectiveNrDensityOfJunctions(
            network, absTol, relTol, junctionType)

    return nu - mu


def calculateEffectiveNrDensityOfNetwork(network: Universum, absTol: float = 1, relTol: float = 1, junctionType=None):
    """
    Compute the effective number density $\\nu_{eff}$ of a network.
    Assumes the precursor-chains to be bifunctional.

    $\\nu_{eff}$ is the number of elastically effective (active) strands per unit volume, 
    which are defined as the ones that can store elastic energy 
    upon network deformation, resp. the effective number density of network strands

    Arguments:
      - network (pylimer_tools.entities.Universum): the network to compute $\\nu_{eff}$ for
      - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol)
      - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol))
      - junctionType: the atom type of the crosslinkers/junctions

    Returns:
      - $\\nu_{eff}$ (float): the effective number density of network strands
    """
    R_taus = []
    chainLengths = []
    for molecule in network.getMolecules(junctionType):
        R_tau = molecule.computeEndToEndDistance()
        R_taus.append(R_tau)
        chainLengths.append(molecule.getLength())

    chainLengths = np.asarray(chainLengths)
    R_taus = np.asarray(R_taus)
    R_tau_max = R_taus.max()
    if (absTol is None):
        absTol = R_tau_max
    numEffective = chainLengths[R_taus >
                                absTol or R_taus > relTol*R_tau_max].sum()

    return numEffective / network.getSize()


def calculateEffectiveNrDensityOfJunctions(network: Universum, absTol: float = 0, relTol: float = 0, junctionType=None, minNumEffectiveStrands=2):
    """
    Compute the number density of the elastically effective crosslinks, 
    defined as the ones that connect at least two elastically effective strands.
    Assumes the precursor-chains to be bifunctional.

    Arguments:
      - network (pylimer_tools.entities.Universum): the network to compute $\\nu_{eff}$ for
      - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol)
      - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol))
      - junctionType: the atom type of the crosslinkers/junctions
      - minNumEffectiveStrands (int): the number of elastically effective strands to qualify a junction as such

    Returns:
      - $\\mu_{eff}$ (float): the effective number density of junctions
    """
    if (junctionType is None):
        return 0.0
    numEffectiveJunctions = 0
    numIneffectiveJunctions = 0

    graph = network.getUnderlyingGraph()
    subgraphs = graph.decompose()
    for subgraph in subgraphs:
        # 1. find junctions
        junctions = subgraph.vs.select(type_eq=junctionType)
        for junction in junctions:
            # some junctions are easily ineffective
            if (junction.degree() < minNumEffectiveStrands):
                numIneffectiveJunctions += 1
            else:
                # 2. follow junction's strands to find those
                vertexNeighbors = junction.neighbors()
                for neighbor in vertexNeighbors:
                    # follow each neighbor
                    lastNeighbor = junction
                    currentNeighbor = neighbor
                    while(True):
                        # follow as long as we do not find a junction
                        nextNeighbors = currentNeighbor.neighbors()
                        # : this is an assumption that makes this function simpler. Assume: only junctions have more than 2 connections
                        assert(len(nextNeighbors) == 2)
                        nextKey = 1 if nextNeighbors[0] == lastNeighbor else 0
                        lastNeighbor = currentNeighbor
                        currentNeighbor = nextNeighbors[nextKey]
                        # found "end" of strand
                        if (currentNeighbor["type"] == junctionType):
                            # 3. decide if strand is elastically effective
                            Ree = currentNeighbor["atom"].computeDistanceTo(
                                junction["atom"])
                            if ((currentNeighbor != junction) and (Ree > absTol)):
                                numEffectiveJunctions += 1
                            else:
                                numIneffectiveJunctions += 1
                            break
                        # dangling chains are non-effective.
                        if (currentNeighbor.degree() < 2):
                            numIneffectiveJunctions += 1
                            break

    return numEffectiveJunctions/(numEffectiveJunctions+numIneffectiveJunctions)


def calculateWeightFractionOfBackbone(network: Universum, crosslinkerType):
    """
    Compute the weight fraction of network backbone in infinite network

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinkerType: the atom type to use to split the molecules
    """
    return 1 - calculateWeightFractionOfDanglingChains(network)


def calculateWeightFractionOfDanglingChains(network: Universum, crosslinkerType, weights=1):
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
    def getWeightOfGraph(graph):
        counts = Counter(graph["type"])
        weightTotal = 0
        for key in counts:
            weightTotal += weights[key]*counts[key]
        return weightTotal

    allChains = network.getChainsWithCrosslinker(crosslinkerType)
    numTotal = network.getSize()
    weightTotal = numTotal
    if (isinstance(weights, dict)):
        weightTotal = getWeightOfGraph(network.getUnderlyingGraph())
    if (type(weights) in (float, int)):
        weightTotal = numTotal*weights

    numDangling = 0
    weightDangling = 0
    for chain in allChains:
        if (chain.getType() is Molecule.MoleculeType.DANGLING_CHAIN):
            numDangling += chain.getLength()
            weightDangling += getWeightOfGraph(chain.getUnderlyingGraph())

    return weightDangling/weightTotal, numDangling/numTotal


def predictGelationPoint(r: float, f: int, g: int = 2) -> float:
    """
    Compute the gelation point $p_{gel}$ as theoretically predicted
    (gelation point = critical extent of reaction for gelation)

    Source:
      - https://www.sciencedirect.com/science/article/pii/003238618990253X

    Arguments:
      - r: the stoichiometric inbalance of reactants
      - f: functionality of the crosslinkers
      - g: functionality of the precursor polymer

    Returns:
      - p_gel: critical extent of reaction for gelation
    """
    # if (r is None):
    #   r = calculateEffectiveCrosslinkerFunctionality(network, junctionType, f)
    return 1/(r*(f-1)*(g-1))


def computeCrosslinkerConversion(network: Universum, junctionType, f: int) -> float:
    """
    Compute the extent of reaction of the crosslinkers 
    (actual functionality divided by target functionality)

    Arguments:
      - network: the poylmer network to do the computation for
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
      - network: the poylmer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network

    Returns:
      - f (float): the (mean) effective crosslinker functionality
    """
    junctionDegrees = calculateEffectiveCrosslinkerFunctionalities(
        network, junctionType)
    return np.mean(junctionDegrees)


def calculateEffectiveCrosslinkerFunctionalities(network: Universum, junctionType):
    """
    Compute the functionality of every crosslinker in the network

    Arguments:
      - network: the poylmer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network

    Returns:
      - junctionDegrees (list[int]): the functionality of every crosslinker
    """
    junctions = network.vs.select(type_eq=junctionType)
    junctionIds = [v.index for v in junctions]
    junctionDegrees = network.degree(junctionIds, mode="all", loops=False)
    return junctionDegrees


def calculateTopologicalFactor(network: Universum, foreignAtomType=None, totalMass=1):
    """
    Compute the topological factor of a polymer network.
    Assumes the precursor-chains to be bifunctional.

    Source:
      - eq. 16 in https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the network to compute the topological factor for
      - foreignAtomType: the type of atoms to ignore
      - totalMass: the $M$ in the respective formula

    Returns:
      - the topological factor $\\Gamma$
    """
    molecules = network.getMolecules(ignoreAtomType=foreignAtomType)
    Gamma_sum = 0
    for molecule in molecules:
        R_tau = molecule.computeEndToEndDistance()
        b = molecule.computeBondLengths().mean()
        Gamma_sum += R_tau*R_tau / (molecule.getLength() * b * b)

    return Gamma_sum / totalMass
