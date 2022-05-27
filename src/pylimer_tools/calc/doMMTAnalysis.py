
import math
import warnings
from collections import Counter

import numpy as np
import pint
import scipy.special
from pylimer_tools.io.unitStyles import UnitStyle
from pylimer_tools_cpp import Universe
from scipy import optimize


def predictShearModulus(network: Universe, unitStyle: UnitStyle, junctionType: int = None, r: float = None, p: float = None, f: int = None, nu: float = None, T: pint.Quantity = None, strandLength: int = None, functionalityPerType: dict = None):
    """
    Predict the shear modulus using MMT Analysis.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the poylmer network to do the computation for
      - unitStyle: the unit style to use to have the results in appropriate units
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - r: the stoichiometric inbalance. Optional if network is specified
      - p: the extent of reaction. Optional if network is specified
      - f: the functionality of the the crosslinker. Optional if network is specified
      - nu: the strand number density (nr of strands per volume) (ideally with units). Optional if network is specified
      - T: the temperature to compute the modulus at. Default: 273.15 K
      - strandLength: the length of the network strands (in nr. of beads). Optional, can be passed to improve performance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. Optional, can be passed to improve performance

    Returns:
      - G: the predicted shear modulus, or `None` if the universe is empty.

    ToDo:
      - Support more than one crosslinker type (as is supported by original formula)
    """
    G_MMT_phantom, G_MMT_entanglement, _, _ = computeModulusDecomposition(
        network, unitStyle, junctionType, r, p, f, nu, T)
    return G_MMT_phantom + G_MMT_entanglement


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

    if (functionalityPerType is None or junctionType not in functionalityPerType):
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
    if (network.getNrOfAtoms() == 0):
        return None
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
    if (functionalityPerType is None or junctionType not in functionalityPerType):
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
        alpha = (((1./(r*p*p)) - 3./4.)**(1./2.) - (1./2.))
        beta = ((r*p*(alpha**3)) + 1 - r*p)
    else:
        def funToRootForAlpha(alpha):
            return r*p**2*alpha**(f-1) - alpha - r*p ** 2 + 1
        alphaSol = optimize.root_scalar(
            funToRootForAlpha, bracket=(0, 1), method='brentq')
        alpha = alphaSol.root
        beta = ((r*p*alpha**(f-1)) + 1 - r*p)  # TODO: reconsider
    return alpha, beta


def computeModulusDecomposition(network: Universe, unitStyle: UnitStyle, junctionType: int = None, r: float = None, p: float = None, f: int = None, nu: float = None, T: pint.Quantity = None, strandLength: int = None, functionalityPerType: dict = None):
    """
    Compute four different estimates of the plateau modulus, using MMT, ANM and PNM.

    Arguments:
      - network: the poylmer network to do the computation for
      - unitStyle: the unit style to use to have the results in appropriate units
      - junctionType: the type of the junctions/crosslinkers to select them in the network
      - r: the stoichiometric inbalance. Optional if network is specified
      - p: the extent of reaction. Optional if network is specified
      - f: the functionality of the the crosslinker. Optional if network is specified
      - nu: the strand number density (nr of strands per volume) (ideally with units). Optional if network is specified
      - T: the temperature to compute the modulus at. Default: 273.15 K
      - strandLength: the length of the network strands (in nr. of beads). Optional, can be passed to improve performance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. Optional, can be passed to improve performance

    Returns:
      - G_MMT_phantom: the phantom contribution to the MMT modulus
      - G_MMT_entanglement: the entanglement contribution to the MMT modulus
      - G_ANM: the ANM estimate of the modulus
      - G_PNM: the PNM estimate of the modulus

    """
    if (junctionType is None and (r is None or f is None or p is None or nu is None)):
        raise ValueError(
            "Either the junctionType or the required variables must be specified")
    if (r is None):
        r = computeStoichiometricInbalance(
            network, junctionType=junctionType, strandLength=strandLength, functionalityPerType=functionalityPerType)
    if (p is None):
        p = computeExtentOfReaction(
            network, junctionType, strandLength=strandLength)
    if (f is None):
        if (functionalityPerType is None):
            f = network.determineFunctionalityPerType()[junctionType]
        else:
            f = functionalityPerType[junctionType]
    if (nu is None):
        nu = len(network.getMolecules(junctionType)) / \
            (network.getVolume()*unitStyle.getBaseUnitOf('volume'))
    if (T is None):
        T = 273.15*unitStyle.getUnderlyingUnitRegistry()('kelvin')

    # affine
    G_ANM = nu*unitStyle.kB*T
    # phantom
    G_PNM = (1-2/f)*nu*unitStyle.kB*T
    # MMT:
    alpha, beta = computeMMsProbabilities(r, p, f)
    GammaMMTSum = 0.0
    for m in range(3, f+1):
        GammaMMTSum += (((m-2.)/2.)*scipy.special.binom(f, m)
                        * (alpha**(f-m))*((1.-alpha)**m))
    GammaMMT = (2*r/f) * GammaMMTSum
    G_MMT_phantom = GammaMMT*nu*unitStyle.kB*T
    # fraction of elastically effective strands. TODO : check adjustment with r
    pel = ((1/(p)) - alpha)**2
    G_MMT_entanglement = 0.22*unitStyle.getUnderlyingUnitRegistry()('MPa')*(pel**2)
    # entanglement part. TODO : check adjustment with r (and where the 0.22 is coming from? Fabian' s fit!)
    return G_MMT_phantom, G_MMT_entanglement, G_ANM, G_PNM


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


def computeStoichiometricInbalance(network: Universe, junctionType: int, strandLength: int = None, functionalityPerType: dict = None, ignoreTypes: list = [], effective: bool = False) -> float:
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
      - ignoreTypes: a list of integers, the types to ignore for the inbalance (e.g. solvent atom types)
      - effective: whether to use the effective functionality (if functionalityPerType is not passed) or the maximum

    Returns:
      - r (float): the stoichiometric inbalance
    """
    if (network.getNrOfAtoms() == 0):
        return 0

    counts = Counter(network.getAtomTypes())

    if (functionalityPerType is None or junctionType not in functionalityPerType):
        functionalityPerType = network.determineEffectiveFunctionalityPerType(
        ) if effective else network.determineFunctionalityPerType()

    if (junctionType not in counts):
        raise ValueError(
            "No junction with type {} seems to have been found in the network".format(junctionType))

    if (strandLength is None):
        strands = network.getMolecules(junctionType)
        strandLength = np.mean([m.getLength() for m in strands if not np.all(
            [a.getType() in ignoreTypes for a in m.getAtoms()])])

    crosslinkerFormableBonds = counts[junctionType] * \
        functionalityPerType[junctionType]
    otherFormableBonds = 0
    for key in counts:
        if (key in ignoreTypes):
            continue
        if (key not in functionalityPerType):
            raise ValueError(
                "Type {} must have an associated functionality".format(key))
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

    if (functionalityPerType is None or crosslinkerType not in functionalityPerType):
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
