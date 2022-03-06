
import math
import warnings
from collections import Counter

import igraph
import numpy as np
from pylimer_tools_cpp import Universe
from scipy import optimize


def predictShearModulus(network: Universe, junctionType, strandLength: int = None, functionalityPerType=None, T: float = 1, k_B: float = 1):
    """
    Predict the shear modulus using MMT Analysis.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the polymer system to predict the shear modulus for
      - junctionType: the type of atoms making up the junctions/crosslinkers
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - T: the temperature in your unit system
      - k_b: Boltzmann's constant in your unit system
      - totalMass: the :math:`M` in the respective formula

    Returns:
      - G: the predicted shear modulus, or `None` if the universe is empty.

    ToDo:
      - Support more than one crosslinker type (as is supported by original formula)
    """
    if (network.getNrOfAtoms() == 0):
        return None
    nu = len(network.getMolecules(junctionType)) / \
        network.getVolume()  # number of chains (network strands) per unit volume
    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType()
    p = computeExtentOfReaction(network, junctionType, functionalityPerType)
    r = computeStoichiometricInbalance(
        network, junctionType, strandLength, functionalityPerType)
    f = functionalityPerType[junctionType]
    alpha, beta = computeMMsProbabilities(r, p, f)
    weightFractions = computeWeightFractions(network)
    return nu * k_B * T * (2*r/f) * (f - 2)/2 * weightFractions[junctionType] * alpha * (1-alpha)**f


def calculateWeightFractionOfDanglingChains(network: Universe, junctionType, strandLength: int = None, functionalityPerType: dict = None) -> float:
    """
    Compute the weight fraction of dangling strands in infinite network

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinkerType: the atom type to use to split the molecules
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance

    Returns:
      - weightFraction $\\Phi_d = 1 - \\Phi_{el}$: weightDangling/weightTotal
    """
    return 1 - calculateWeightFractionOfBackbone(network, junctionType, strandLength, functionalityPerType)


def calculateWeightFractionOfBackbone(network: Universe, junctionType, strandLength: int = None, functionalityPerType: dict = None) -> float:
    """
    Compute the weight fraction of the backbone strands in an infinite network

    Source:
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae)

    Arguments:
      - network: the poylmer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 
          See: #computeExtentOfReaction

    Returns:
      - :math:`\\Phi_{el}`: weight fraction of network backbone
    """
    if (network.getNrOfAtoms() == 0):
        return 0

    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType()

    W_sol, weightFractions, alpha, beta = computeWeightFractionOfSolubleMaterial(
        network, junctionType, strandLength, functionalityPerType)

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


def measureWeightFractioOfSolubleMaterial(network: Universe, relTol: float = 0.75, absTol: float = None) -> float:
    """
    Compute the weight fraction of soluble material by counting.

    Arguments:
      - network: the poylmer network to do the computation for
      - relTol: the fraction of the maximum weight that counts as soluble. Ignored if absTol is specified
      - absTol: the weight from which on a component is not soluble anymore

    Returns:
      - :math:`W_{sol}` (float): the weight fraction of soluble material as counted.

    """
    fractions = network.getClusters()
    weights = np.array([f.computeWeight() for f in fractions])
    totalWeight = weights.sum()
    solubleWeight = 0
    for w in weights:
        if (absTol is not None):
            if (w < absTol):
                solubleWeight += w
        else:
            if (w < relTol*weights.max()):
                solubleWeight += w

    return solubleWeight/totalWeight


def computeWeightFractionOfSolubleMaterial(network: Universe = None, junctionType: int = 2, strandLength: int = None, functionalityPerType: dict = None, weightFractions: dict = None, r: float = None, p: float = None) -> float:
    """
    Compute the weight fraction of soluble material by MMT.

    Source:
      - https://pubs.acs.org/doi/10.1021/ma00046a021
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737

    Arguments:
      - network: the poylmer network to do the computation for
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - weightFractions (dict): a dictionary with key: type, and value: weight fraction of type. Pass if you want to omit the network.
      - strandLength (int): the length of the network strands (in nr. of beads). 
          See: :func:`~pylimer_tools.calc.doMMTAnalysis.computeStoichiometricInbalance`.
      - functionalityPerType (dict): a dictionary with key: type, and value: functionality of this atom type. 
          See: :func:`~pylimer_tools.calc.doMMTAnalysis.computeExtentOfReaction`.

    Returns:
      - :math:`W_{sol}` (float): the weight fraction of soluble material according to MMT.
      - weightFractions (dict): a dictionary with key: type, and value: weight fraction of type
      - :math:`\\alpha` (float): Macosko & Miller's :math:`P(F_A)`
      - :math:`\\beta` (float): Macosko & Miller's :math:`P(F_B)`
    """
    if (functionalityPerType is None):
        assert(network is not None)
        functionalityPerType = network.determineFunctionalityPerType()

    if (functionalityPerType[junctionType] not in [3, 4]):
        raise NotImplementedError(
            "Currently, only crosslinker functionality of 3 or 4 is supported. {} given.".format(functionalityPerType[junctionType]))

    for key in functionalityPerType:
        if (key != junctionType and functionalityPerType[key] != 2):
            raise NotImplementedError(
                "Currently, only strand functionality of 2 is supported. {} given for type {}".format(functionalityPerType[key], key))

    if (weightFractions is None):
        weightFractions = computeWeightFractions(network)

    if (p is None):
        assert(network is not None)
        p = computeExtentOfReaction(
            network, junctionType, functionalityPerType)
    if (r is None):
        assert(network is not None)
        r = computeStoichiometricInbalance(
            network, junctionType, strandLength=strandLength, functionalityPerType=functionalityPerType)

    alpha, beta = computeMMsProbabilities(
        r, p, functionalityPerType[junctionType])
    W_sol = 0
    for key in weightFractions:
        coeff = alpha if key == junctionType else beta
        W_sol += weightFractions[key] * \
            (math.pow(coeff, functionalityPerType[key]))

    return W_sol, weightFractions, alpha, beta


def computeMMsProbabilities(r, p, f):
    """
    Compute Macosko and Miller's probabilities :math:`P(F_A)` and :math:`P(F_B)`

    Sources:
      - https://pubs.acs.org/doi/10.1021/ma60050a004
      - https://doi.org/10.1021/ma60050a003

    Arguments:
      - r: the stoichiometric inbalance
      - p: the extent of reaction
      - f: the functionality of the the crosslinker

    Returns:
      - alpha: :math:`P(F_A)`
      - beta: :math:`P(F_B)`    
    """
    # first, check a few things required by the formulae
    # since we want alpha, beta \in [0,1], given they are supposed to be probabilities
    # if (r > 1 or r < 0):
    #     raise ValueError(
    #         "A stoichiometric inbalance ouside of [0, 1] is not (yet) supported. Got {}".format(r))
    # if (p < 1/math.sqrt(2) or p > 1):
    #     raise ValueError(
    #         "The extent of reaction has to be inside [1/sqrt(2), 1] for the result to be realistic. Got {}".format(p))
    # if (r <= 1/(2*p*p) and f == 3):
    #     raise ValueError(
    #         "The stoichiometric inbalance must be > 1/(2p^2) for the resulting alpha to be realisitic. Got p = {}, r = {}".format(p, r))

    # actually do the calculations
    if (f == 3):
        alpha = ((1 - r*p*p)/(r*p*p))
        beta = (r*p*alpha**2)
    elif(f == 4):
        alpha = (math.sqrt((1/(r*p*p)) - 0.75) - 0.5)
        beta = ((r*p*(alpha**3)) + 1 - r*p)
    else:
        def funToRootForAlpha(alpha):
            return r*p**2*alpha**(f-1) - alpha - r*p ** 2 + 1
        alphaSol = optimize.root_scalar(
            funToRootForAlpha, bracket=(0, 1), method='brentq')
        alpha = alphaSol.root
        beta = ((r*p*alpha**(f-1)) + 1 - r*p)  # TODO: reconsider
    return alpha, beta


def computeWeightFractions(network: Universe) -> dict:
    """
    Compute the weight fractions of each atom type in the network.

    Arguments:
      - network: the poylmer network to do the computation for

    Returns:
      - :math:`\\vec{W_i}` (dict): using the type i as a key, this dict contains the weight fractions (:math:`\\frac{W_i}{W_{tot}}`)
    """
    if (network.getNrOfAtoms() == 0):
        return {}

    weightPerType = network.getMasses()
    counts = Counter(network.getAtomTypes())
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


def computeStoichiometricInbalance(network: Universe, junctionType, strandLength: int = None, functionalityPerType: dict = None) -> float:
    """
    Compute the stoichiometric inbalance
    ( nr. of bonds formable of crosslinker / nr. of formable bonds of precursor )

    NOTE: 
      if your system has a non-integer number of possible bonds (e.g. one site unbonded),
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
    if (network.getNrOfAtoms() == 0):
        return 0

    counts = Counter(network.getAtomTypes())

    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType()

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


def computeExtentOfReaction(network: Universe, crosslinkerType, functionalityPerType: dict = None, strandLength: float = None) -> float:
    """
    Compute the extent of polymerization reaction
    (nr. of formed bonds in reaction / max. nr. of bonds formable)
    NOTE: if your system has a non-integer number of possible bonds (e.g. one site unbonded),
    this will not be rounded/respected in any way. 

    Arguments:
      - network: the poylmer network to do the computation for
      - crosslinkerType: the atom type of crosslinker beads
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 
          If None: will use max functionality per type.
      - strandLength: the length of the network strands (in nr. of beads). 
          If None: will compute from network structure

    Returns:
      - p (float): the extent of reaction
    """

    if (network.getNrOfAtoms() == 0):
        return 1

    if (functionalityPerType is None):
        functionalityPerType = network.determineFunctionalityPerType()

    numStrands = len(network.getMolecules(crosslinkerType))
    numCrosslinkers = len(network.getAtomsWithType(crosslinkerType))

    # assuming strand has functionality 2
    maxFormableBonds = min(numStrands*2, numCrosslinkers *
                           functionalityPerType[crosslinkerType])

    if (maxFormableBonds == 0):
        return 1

    if (strandLength is None):
        strands = network.getMolecules(crosslinkerType)
        strandLength = np.mean([m.getLength() for m in strands])

    actuallyFormedBonds = (network.getNrOfBonds() -
                           (numStrands * (strandLength-1)))

    return actuallyFormedBonds/(maxFormableBonds)


def predictGelationPoint(r: float, f: int, g: int = 2) -> float:
    """
    Compute the gelation point :math:`p_{gel}` as theoretically predicted
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
    return math.sqrt(1/(r*(f-1)*(g-1)))
