
import math
from collections import Counter

import igraph
import numpy as np
from pylimer_tools.entities.universum import Universum


def predictShearModulus(network: Universum, junctionType, weightPerType, strandLength: int = None, functionalityPerType=None, T: float = 1, k_B: float = 1, totalMass=1):
    """
    Predict the shear modulus using MMT Analysis.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the polymer system to predict the shear modulus for
      - junctionType: the type of atoms making up the junctions/crosslinkers
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - weightPerType: a dictionary with key: type, and value: weight per atom of this atom type. See: #computeWeightFractions
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - T: the temperature in your unit system
      - k_b: Boltzmann's constant in your unit system
      - totalMass: the $M$ in the respective formula

    Returns:
      - G: the predicted shear modulus, or `None` if the universe is empty.

    ToDo:
      - Support more than one crosslinker type (as is supported by original formula)
    """
    if (network.getSize() == 0):
        return None
    nu = len(network.getMolecules(junctionType)) / \
        network.getVolume()  # number of chains (network strands) per unit volume
    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType()
    p = computeExtentOfReaction(network, functionalityPerType)
    r = computeStoichiometricInbalance(
        network, junctionType, strandLength, functionalityPerType)
    f = functionalityPerType[junctionType]
    alpha, beta = computeMMsProbabilities(r, p, f)
    weightFractions = computeWeightFractions(network, weightPerType)
    return nu * k_B * T * (2*r/f) * (f - 2)/2 * weightFractions[junctionType] * alpha * (1-alpha)**f


def calculateWeightFractionOfDanglingChains(network: Universum, junctionType, weightPerType, strandLength: int = None, functionalityPerType: dict = None) -> float:
    """
    Compute the weight fraction of dangling strands in infinite network

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinkerType: the atom type to use to split the molecules
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - weights: either a dict with key: atomType and value: weight, or a scalar value if all atoms have the same weight

    Returns:
      - weightFraction $\\Phi_d = 1 - \\Phi_{el}$: weightDangling/weightTotal
    """
    return 1 - calculateWeightFractionOfBackbone(network, junctionType, weightPerType, strandLength, functionalityPerType)


def calculateWeightFractionOfBackbone(network: Universum, junctionType, weightPerType, strandLength: int = None, functionalityPerType: dict = None) -> float:
    """
    Compute the weight fraction of the backbone strands in an infinite network

    Source:
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae)

    Arguments:
      - network: the poylmer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - weightPerType: a dictionary with key: type, and value: weight per atom of this atom type. See: #computeWeightFractions
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 
          See: #computeExtentOfReaction

    Returns:
      - $\\Phi_{el}$: weight fraction of network backbone
    """
    if (network.getSize() == 0):
        return 0

    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType()

    W_sol, weightFractions, alpha, beta = computeWeightFractionOfSolubleMaterial(
        network, junctionType, weightPerType, strandLength, functionalityPerType)

    Phi_el = 0
    W_a = weightFractions[junctionType]/functionalityPerType[junctionType]
    W_xl = weightFractions[junctionType]
    W_x2 = 1-W_xl
    if (functionalityPerType[junctionType] == 3):
        Phi_el = ((W_x2*(1-beta)**2) +
                  (W_xl*((1-alpha)**3 + 3*alpha*(1-W_a)*((1-alpha)**2))))/(1-W_sol)
    else:
        assert(functionalityPerType[junctionType] == 4)
        Phi_el = ((W_x2*(1-beta)**2) +
                  (W_xl*(((1-alpha)**4) + 4*alpha*(1-W_a) * ((1-alpha)**3) +
                         6*(alpha**2)*(1-2*W_a)*(1-alpha)**2)))/(1-W_sol)

    return Phi_el


def computeWeightFractionOfSolubleMaterial(network: Universum, junctionType, weightPerType, strandLength: int = None, functionalityPerType: dict = None) -> float:
    """
    Compute the weight fraction of soluble material.

    Source:
      - https://pubs.acs.org/doi/10.1021/ma00046a021
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737

    Arguments:
      - network: the poylmer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - weightPerType: a dictionary with key: type, and value: weight per atom of this atom type. See: #computeWeightFractions
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 
          See: #computeExtentOfReaction

    Returns:
      - $W_{sol}$ (float): the weight fraction of soluble material according to MMT.
      - weightFractions (dict): a dictionary with key: type, and value: weight fraction of type
      - $\\alpha$ (float): Macosko & Miller's $P(F_A)$
      - $\\beta$ (float): Macosko & Miller's $P(F_B)$
    """
    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType()

    if (functionalityPerType[junctionType] not in [3, 4]):
        raise NotImplementedError(
            "Currently, only crosslinker functionality of 3 or 4 is supported. {} given.".format(functionalityPerType[junctionType]))

    for key in functionalityPerType:
        if (key != junctionType and functionalityPerType[key] != 2):
            raise NotImplementedError(
                "Currently, only strand functionality of 2 is supported. {} given for type {}".format(functionalityPerType[key], key))

    p = computeExtentOfReaction(network, functionalityPerType)
    r = computeStoichiometricInbalance(
        network, junctionType, strandLength, functionalityPerType)

    alpha, beta = computeMMsProbabilities(
        r, p, functionalityPerType[junctionType])
    weightFractions = computeWeightFractions(network, weightPerType)
    W_sol = 0
    for key in weightFractions:
        coeff = alpha if key == junctionType else beta
        W_sol += weightFractions[key]*coeff**functionalityPerType[key]

    return W_sol, weightFractions, alpha, beta


def computeMMsProbabilities(r, p, f):
    """
    Compute Macosko and Miller's probabilities $P(F_A)$ and $P(F_B)$

    Arguments:
      - r: the stoichiometric inbalance
      - p: the extent of reaction
      - f: the functionality of the the crosslinker

    Returns:
      - alpha: $P(F_A)$
      - beta: $P(F_B)$    
    """
    # first, check a few things required by the formulae
    # since we want alpha, beta \in [0,1], given they are supposed to be probabilities
    if (r > 1 or r < 0):
        raise ValueError(
            "A stoichiometric inbalance ouside of [0, 1] is not (yet) supported. Got {}".format(r))
    if (p < 1/math.sqrt(2) or p > 1):
        raise ValueError(
            "The extent of reaction has to be inside [1/sqrt(2), 1] for the result to be realistic. Got {}".format(p))
    if (r <= 1/(2*p*p)):
        raise ValueError(
            "The stoichiometric inbalance must be > 1/(2p^2) for the resulting alpha to be realisitic. Got p = {}, r = {}".format(p, r))

    # actually do the calculations
    if (f == 3):
        alpha = ((1 - r*p*p)/(r*p*p))
        beta = (r*p*alpha*alpha)
    else:
        assert(f == 4)
        alpha = (math.sqrt((1/(r*p*p)) - 0.5) - 0.5)
        beta = ((r*p*((math.sqrt((1/(r*p*p)) - 0.75) - 0.5)**3)) + 1 - r*p)
    return alpha, beta


def computeWeightFractions(network: Universum, weightPerType) -> dict:
    """
    Compute the weight fractions of each atom type in the network.

    Arguments:
      - network: the poylmer network to do the computation for
      - weightPerType: a dictionary with key: type, and value: weight per atom of this atom type. 

    Returns:
      - $\\vec{W_i}$ (dict): using the type i as a key, this dict contains the weight fractions ($\\frac{W_i}{W_{tot}}$)
    """
    if (network.getSize() == 0):
        return {}

    graph: igraph.Graph = network.getUnderlyingGraph()
    counts = Counter(graph.vs["type"])
    totalMass = 0
    partialMasses = {}
    for key in counts:
        totalMass += counts[key]*weightPerType[key]
        partialMasses[key] = counts[key]*weightPerType[key]

    if (totalMass == 0):
        return partialMasses

    weightFractions = {}
    for key in partialMasses:
        weightFractions[key] = partialMasses[key]/totalMass

    return weightFractions


def computeStoichiometricInbalance(network: Universum, junctionType, strandLength: int = None, functionalityPerType: dict = None) -> float:
    """
    Compute the stoichiometric inbalance
    ( nr. of bonds formable of crosslinker / nr. of formable bonds of precursor )

    NOTE: 
      - if your system has a non-integer number of possible bonds (e.g. one site unbonded),
          this will not be rounded/respected in any way. 

    Arguments:
      - network: the poylmer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - strandLength: the length of the network strands (in nr. of beads). 
          Used to infer the number of precursor strands. 
          If `None`: will use average length of each connected system when ignoring the crosslinkers.
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 
          If `None`: will use max functionality per type.

    Returns:
      - r (float): the stoichiometric inbalance
    """
    if (network.getSize() == 0):
        return 0

    graph: igraph.Graph = network.getUnderlyingGraph()
    counts = Counter(graph.vs["type"])

    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType(counts)

    if (strandLength is None):
        strands = network.getMolecules(junctionType)
        strandLength = np.mean([m.getLength() for m in strands])

    crosslinkerFormableBonds = counts[junctionType] * \
        functionalityPerType[junctionType]
    otherFormableBonds = 0
    for key in counts:
        if (key != junctionType):
            otherFormableBonds += counts[key]*functionalityPerType[key]

    # division by 2 is implicit
    return crosslinkerFormableBonds/(otherFormableBonds/strandLength)


def computeExtentOfReaction(network: Universum, functionalityPerType: dict = None) -> float:
    """
    Compute the extent of reaction
    (nr. of formed bonds in reaction / max. nr. of bonds formable)
    NOTE: if your system has a non-integer number of possible bonds (e.g. one site unbonded),
    this will not be rounded/respected in any way. 

    Arguments:
      - network: the poylmer network to do the computation for
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 
          If None: will use max functionality per type.

    Returns:
      - p (float): the extent of reaction
    """

    if (network.getSize() == 0):
        return 1

    graph: igraph.Graph = network.getUnderlyingGraph()
    counts = Counter(graph.vs["type"])
    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType(counts)

    maxFormableBonds = 0
    for key in counts:
        maxFormableBonds += functionalityPerType[key]*counts[key]

    if (maxFormableBonds == 0):
        return 1

    graph.simplify()
    # multiplication by 2 as each bond affects 2 possible bonds
    return graph.ecount()*2.0/(maxFormableBonds)


def predictGelationPoint(r: float, f: int, g: int = 2) -> float:
    """
    Compute the gelation point $p_{gel}$ as theoretically predicted
    (gelation point = critical extent of reaction for gelation)

    Source:
      - https://www.sciencedirect.com/science/article/pii/003238618990253X

    Arguments:
      - r (double): the stoichiometric inbalance of reactants (see: #computeStoichiometricInbalance)
      - f (int): functionality of the crosslinkers
      - g (int): functionality of the precursor polymer

    Returns:
      - p_gel: critical extent of reaction for gelation
    """
    # if (r is None):
    #   r = calculateEffectiveCrosslinkerFunctionality(network, junctionType, f)
    return 1/(r*(f-1)*(g-1))
