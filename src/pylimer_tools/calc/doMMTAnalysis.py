
import math
import warnings
from collections import Counter

import numpy as np
import pint
import pylimer_tools.calc.doMEHPAnalysis as mehp
import scipy.special
from pylimer_tools.io.unitStyles import UnitStyle
from pylimer_tools_cpp import Universe
from scipy import optimize


def predictShearModulus(network: Universe, unitStyle: UnitStyle, crosslinkerType: int = None, r: float = None, p: float = None, f: int = None, nu: float = None, T: pint.Quantity = None, strandLength: int = None, functionalityPerType: dict = None):
    """
    Predict the shear modulus using MMT Analysis.

    Source:
      - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the poylmer network to do the computation for
      - unitStyle: the unit style to use to have the results in appropriate units
      - crosslinkerType: the type of the junctions/crosslinkers to select them in the network
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
        network, unitStyle, crosslinkerType, r, p, f, nu, T)
    return G_MMT_phantom + G_MMT_entanglement


def predictNumberDensityOfJunctionPoints(network: Universe, crosslinkerType: int, strandLength: int = None, functionalityPerType: dict = None) -> float:
    """
    Compute the number density of network strands using MMT

    Source:
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae)

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinkerType: the atom type to use to split the molecules
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 

    Returns:
      - mu: The predicted number density of junction points
    """
    if (functionalityPerType is None or crosslinkerType not in functionalityPerType):
        functionalityPerType = network.determineFunctionalityPerType()

    _, weightFractions, alpha, _ = computeWeightFractionOfSolubleMaterial(
        network, crosslinkerType, strandLength, functionalityPerType)

    if (functionalityPerType[crosslinkerType] == 3):
        return weightFractions[crosslinkerType]*(1-alpha)**3
    elif (functionalityPerType[crosslinkerType] == 4):
        return weightFractions[crosslinkerType]*(4*alpha*(1-alpha)**3 + (1-alpha)**4)
    else:
        raise NotImplementedError(
            "Currently, only cross-linker functionalities of 3 and 4 are supported")


def predictNumberDensityOfNetworkStrands(network: Universe, crosslinkerType: int, strandLength: int = None, functionalityPerType: dict = None, r: float = None, p: float = None) -> float:
    """
    Compute the number density of network strands using MMT

    Source:
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae)

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinkerType: the atom type to use to split the molecules
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 

    Returns:
      - nu: The predicted number density of network strands
    """
    if (functionalityPerType is None or crosslinkerType not in functionalityPerType):
        functionalityPerType = network.determineFunctionalityPerType()

    weightFractions = network.computeWeightFractions()
    alpha, _ = computeMMsProbabilities(r=r if r is not None else computeStoichiometricInbalance(network, crosslinkerType, strandLength, functionalityPerType),
                                       p=p if p is not None else mehp.computeCrosslinkerConversion(
                                           network, crosslinkerType, functionalityPerType[crosslinkerType]),
                                       f=functionalityPerType[crosslinkerType])

    if (functionalityPerType[crosslinkerType] == 3):
        return (3/2)*weightFractions[crosslinkerType]*(1-alpha)**3
    elif (functionalityPerType[crosslinkerType] == 4):
        return weightFractions[crosslinkerType]*(6*alpha*(1-alpha)**3 + 2*(1-alpha)**4)
    else:
        raise NotImplementedError(
            "Currently, only junction functionalities of 3 and 4 are supported")


def calculateWeightFractionOfDanglingChains(network: Universe, crosslinkerType: int, strandLength: int = None, functionalityPerType: dict = None, weightFractions: dict = None, r: float = None, p: float = None) -> float:
    """
    Compute the weight fraction of dangling strands in infinite network

    Arguments:
      - network: the network to compute the weight fraction for
      - crosslinkerType: the atom type to use to split the molecules
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 

    Returns:
      - weightFraction $\\Phi_d = 1 - \\Phi_{el}$: weightDangling/weightTotal
    """
    return 1 - calculateWeightFractionOfBackbone(network, crosslinkerType, strandLength, functionalityPerType, weightFractions, r, p)


def calculateWeightFractionOfBackbone(network: Universe, crosslinkerType: int, strandLength: int = None, functionalityPerType: dict = None, weightFractions: dict = None, r: float = None, p: float = None) -> float:
    """
    Compute the weight fraction of the backbone strands in an infinite network

    Source:
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae)
      - https://pubs.acs.org/doi/10.1021/ma00046a021 (see appendix)

    Arguments:
      - network: the poylmer network to do the computation for
      - crosslinkerType: the type of the junctions/crosslinkers to select them in the network
      - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. 

    Returns:
      - :math:`\\Phi_{el}`: weight fraction of network backbone
    """
    if (network is not None and network.getNrOfAtoms() == 0):
        return 0

    if (functionalityPerType is None or crosslinkerType not in functionalityPerType):
        functionalityPerType = network.determineFunctionalityPerType()

    W_sol, weightFractions, alpha, beta = computeWeightFractionOfSolubleMaterial(
        network, crosslinkerType, strandLength, functionalityPerType, weightFractions, r, p)
    if (W_sol < 0 or W_sol > 1):
        warnings.warn(
            "The weight fraction W_sol predicted by MMT ({}) is outside accepted range. Falling back to measurement.".format(W_sol))
        W_sol = measureWeightFractioOfSolubleMaterial(network)

    Phi_el = 0
    W_a = weightFractions[crosslinkerType] / \
        functionalityPerType[crosslinkerType]
    W_xl = weightFractions[crosslinkerType]
    W_x2 = 1-W_xl
    assert(W_a <= 1 and W_a >= 0)
    assert(W_xl <= 1 and W_xl >= 0)
    assert(W_x2 <= 1 and W_x2 >= 0)
    assert(W_sol <= 1 and W_sol >= 0)
    if (functionalityPerType[crosslinkerType] == 3):
        Phi_el = ((W_x2*(1-beta)**2) +
                  (W_xl*((1-alpha)**3 + 3*alpha*(1-W_a)*((1-alpha)**2))))/(1-W_sol)
    else:
        assert(functionalityPerType[crosslinkerType] == 4)
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
    weights = np.array([f.computeTotalMass() for f in fractions])
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


def computeWeightFractionOfSolubleMaterial(network: Universe, crosslinkerType: int, strandLength: int = None, functionalityPerType: dict = None, weightFractions: dict = None, r: float = None, p: float = None) -> float:
    """
    Compute the weight fraction of soluble material by MMT.

    Source:
      - https://pubs.acs.org/doi/10.1021/ma00046a021
      - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737

    Arguments:
      - network: the poylmer network to do the computation for
      - crosslinkerType: the type of the junctions/crosslinkers to select them in the network
      - weightFractions (dict): a dictionary with key: type, and value: weight fraction of type. Pass if you want to omit the network.
      - strandLength (int): the length of the network strands (in nr. of beads). 
          See: :func:`~pylimer_tools.calc.doMMTAnalysis.computeStoichiometricInbalance`.
      - functionalityPerType (dict): a dictionary with key: type, and value: functionality of this atom type. 
          See: :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.determineFunctionalityPerType`.

    Returns:
      - :math:`W_{sol}` (float): the weight fraction of soluble material according to MMT.
      - weightFractions (dict): a dictionary with key: type, and value: weight fraction of type
      - :math:`\\alpha` (float): Macosko & Miller's :math:`P(F_A)`
      - :math:`\\beta` (float): Macosko & Miller's :math:`P(F_B)`
    """
    if (functionalityPerType is None or crosslinkerType not in functionalityPerType):
        assert(network is not None)
        functionalityPerType = network.determineFunctionalityPerType()

    if (functionalityPerType[crosslinkerType] not in range(3, 7)):
        raise NotImplementedError(
            "Currently, a crosslinker functionality of {} is not supported.".format(functionalityPerType[crosslinkerType]))

    for key in functionalityPerType:
        if (key != crosslinkerType and functionalityPerType[key] != 2):
            raise NotImplementedError(
                "Currently, only strand functionality of 2 is supported. {} given for type {}".format(functionalityPerType[key], key))

    if (weightFractions is None):
        weightFractions = computeWeightFractions(network)

    if (p is None):
        assert(network is not None)
        p = mehp.computeCrosslinkerConversion(
            network, crosslinkerType, functionalityPerType[crosslinkerType])
        if (p < 0 or p > 1):
            warnings.warn(
                "The p computed ({}) is outside the accepted range. Falling back to effective cross-linker functionality.".format(p))
            p = mehp.calculateEffectiveCrosslinkerFunctionality(
                network, crosslinkerType)
    assert(p <= 1 and p >= 0)
    if (r is None):
        assert(network is not None)
        r = computeStoichiometricInbalance(
            network, crosslinkerType, strandLength=strandLength, functionalityPerType=functionalityPerType)
    assert(r >= 0)

    alpha, beta = computeMMsProbabilities(
        r, p, functionalityPerType[crosslinkerType])
    assert(alpha <= 1 and alpha >= 0)
    assert(beta <= 1 and beta >= 0)
    W_sol = 0
    for key in weightFractions:
        coeff = alpha if key == crosslinkerType else beta
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
    if (p > 1 or p < 0):
        raise ValueError(
            "An extent of reaction ouside of [0, 1] is not supported. Got {}".format(p))
    # if (p < 1/math.sqrt(2) or p > 1):
    #     raise ValueError(
    #         "The extent of reaction has to be inside [1/sqrt(2), 1] for the result to be realistic. Got {}".format(p))
    # if (r <= 1/(2*p*p) and f == 3):
    #     raise ValueError(
    #         "The stoichiometric inbalance must be > 1/(2p^2) for the resulting alpha to be realisitic. Got p = {}, r = {}".format(p, r))

    # actually do the calculations
    if (f == 3):
        alpha = ((1 - r*p*p)/(r*p*p))
        if (not (1/(p**2) < 2*r and (1/(p**2) > r))):
            warnings.warn(
                "The resulting P(F_A) is probably unreliable, as the detected root does not fulfill the required conditions.")
    elif(f == 4):
        alpha = (((1./(r*p*p)) - 3./4.)**(1./2.) - (1./2.))
        if (not (1/(p**2) < 3*r and (1/(p**2) > r))):
            warnings.warn(
                "The resulting P(F_A) is probably unreliable, as the detected root does not fulfill the required conditions.")
    else:
        def funToRootForAlpha(alpha):
            return r*p**2*alpha**(f-1) - alpha - r*(p ** 2) + 1

        def funToRootForAlphaPrime(alpha):
            return -1 + alpha**(f-2)*(-1+f)*(p**2)*r

        def funToRootForAlphaPrime2(alpha):
            return alpha**(f-3)*(-2+f)*(-1+f)*(p**2)*r

        alphaSol = optimize.root_scalar(
            funToRootForAlpha, bracket=(0, 1), method='halley', fprime=funToRootForAlphaPrime, fprime2=funToRootForAlphaPrime2, x0=0.5)
        alpha = alphaSol.root
    beta = r*p*alpha**(f-1) + 1 - r*p
    if (alpha > 1 or alpha < 0):
        warnings.warn(
            "The resulting P(F_A) is probably unreliable, as it will be clipped to [0,1] from {}".format(alpha))
    if (beta > 1 or beta < 0):
        warnings.warn(
            "The resulting P(F_B) is probably unreliable, as it will be clipped to [0,1] from {}".format(beta))
    return np.clip(alpha, 0, 1), np.clip(beta, 0, 1)  # TODO: reconsider


def computeModulusDecomposition(network: Universe, unitStyle: UnitStyle, crosslinkerType: int = None, r: float = None, p: float = None, f: int = None, nu: float = None, T: pint.Quantity = None, strandLength: int = None, functionalityPerType: dict = None, Ge1=0.22):
    """
    Compute four different estimates of the plateau modulus, using MMT, ANM and PNM.

    Arguments:
      - network: the polymer network to do the computation for
      - unitStyle: the unit style to use to have the results in appropriate units
      - crosslinkerType: the type of the junctions/crosslinkers to select them in the network
      - r: the stoichiometric inbalance. Optional if network is specified
      - p: the extent of reaction. Optional if network is specified
      - f: the functionality of the the crosslinker. Optional if network is specified
      - nu: the strand number density (nr of strands per volume) (ideally with units). Optional if network is specified
      - T: the temperature to compute the modulus at. Default: 298.15 K
      - strandLength: the length of the network strands (in nr. of beads). Optional, can be passed to improve performance
      - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type. Optional, can be passed to improve performance

    Returns:
      - G_MMT_phantom: the phantom contribution to the MMT modulus
      - G_MMT_entanglement: the entanglement contribution to the MMT modulus
      - G_ANM: the ANM estimate of the modulus
      - G_PNM: the PNM estimate of the modulus

    """
    if ((crosslinkerType is None or network is None) and (r is None or f is None or p is None or nu is None)):
        raise ValueError(
            "Either the network and crosslinkerType or the required variables must be specified")
    if (r is None):
        r = computeStoichiometricInbalance(
            network, crosslinkerType=crosslinkerType, strandLength=strandLength, functionalityPerType=functionalityPerType)
    if (p is None):
        p = mehp.computeCrosslinkerConversion(
            network, crosslinkerType, functionalityPerType=functionalityPerType)
        if (p < 0 or p > 1):
            warnings.warn(
                "The p computed ({}) is outside the accepted range. Falling back to effective cross-linker functionality.".format(p))
            p = mehp.calculateEffectiveCrosslinkerFunctionality(
                network, crosslinkerType)
    if (f is None):
        if (functionalityPerType is None):
            f = network.determineFunctionalityPerType()[crosslinkerType]
        else:
            f = functionalityPerType[crosslinkerType]
    if (nu is None):
        nu = len(network.getMolecules(crosslinkerType)) / \
            (network.getVolume()*unitStyle.getBaseUnitOf('volume'))
    if (T is None):
        T = (273.15+25)*unitStyle.getUnderlyingUnitRegistry()('kelvin')

    # affine
    G_ANM = nu*unitStyle.kB*T
    # phantom
    G_PNM = (1-2/f)*nu*unitStyle.kB*T
    # MMT:
    alpha, beta = computeMMsProbabilities(r, p, f)
    GammaMMTSum = 0.0
    for m in range(3, f+1):
        GammaMMTSum += (((m-2)/2) *
                        computeProbabilityThatMonomerIsEffective(f, m, alpha))
    GammaMMT = (2*r/f) * GammaMMTSum
    G_MMT_phantom = GammaMMT*nu*unitStyle.kB*T
    # fraction of elastically effective strands. TODO : check adjustment with r
    pel = ((1/(p)) - alpha/p)**2
    G_MMT_entanglement = Ge1*unitStyle.getUnderlyingUnitRegistry()('MPa')*(pel**2)
    # entanglement part. TODO : check adjustment with r (and where the 0.22 is coming from? Fabian' s fit!)
    return G_MMT_phantom, G_MMT_entanglement, G_ANM, G_PNM


def computeProbabilityThatMonomerIsEffective(functionalityOfMonomer: int, expectedDegreeOfEffect: int, PFaout: float):
    """
    Compute the probability that an Af, monomer will be an effective cross-link of degree m 

    Source:
        - Eq. 45 in Miller, Macosko 1976, A New Derivation of Post Gel Properties of Network

    Arguments:
        - functionalityOfMonomer: f
        - 
    """
    f = functionalityOfMonomer
    m = expectedDegreeOfEffect
    alpha = PFaout
    return scipy.special.binom(f, m) * (alpha**(f-m))*((1.-alpha)**m)


def computeWeightFractions(network: Universe) -> dict:
    """
    Compute the weight fractions of each atom type in the network.

    Arguments:
      - network: the poylmer network to do the computation for

    Returns:
      - :math:`\\vec{W_i}` (dict): using the type i as a key, this dict contains the weight fractions (:math:`\\frac{W_i}{W_{tot}}`)
    """
    return network.computeWeightFractions()


def computeStoichiometricInbalance(network: Universe, crosslinkerType: int, strandLength: int = None, functionalityPerType: dict = None, ignoreTypes: list = [], effective: bool = False) -> float:
    """
    Compute the stoichiometric inbalance
    ( nr. of bonds formable of crosslinker / nr. of formable bonds of precursor )

    NOTE: 
      if your system has a non-integer number of possible bonds (e.g. one site unbonded),
      this will not be rounded/respected in any way. 

    Arguments:
      - network: the poylmer network to do the computation for
      - crosslinkerType: the type of the junctions/crosslinkers to select them in the network
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

    if (functionalityPerType is None or crosslinkerType not in functionalityPerType):
        functionalityPerType = network.determineEffectiveFunctionalityPerType(
        ) if effective else network.determineFunctionalityPerType()

    if (crosslinkerType not in counts):
        raise ValueError(
            "No junction with type {} seems to have been found in the network".format(crosslinkerType))

    if (strandLength is None):
        strands = network.getMolecules(crosslinkerType)
        strandLength = np.mean([m.getLength() for m in strands if not np.all(
            [a.getType() in ignoreTypes for a in m.getAtoms()])])

    crosslinkerFormableBonds = counts[crosslinkerType] * \
        functionalityPerType[crosslinkerType]
    otherFormableBonds = 0
    for key in counts:
        if (key in ignoreTypes):
            continue
        if (key not in functionalityPerType):
            raise ValueError(
                "Type {} must have an associated functionality".format(key))
        if (key != crosslinkerType):
            otherFormableBonds += counts[key]*functionalityPerType[key]

    # division by 2 is implicit
    return crosslinkerFormableBonds/(otherFormableBonds/strandLength)


def computeExtentOfReaction(network: Universe, crosslinkerType, functionalityPerType: dict = None, strandLength: float = None) -> float:
    """
    Compute the extent of polymerization reaction
    (nr. of formed bonds in reaction / max. nr. of bonds formable)
    NOTE: 
        - if your system has a non-integer number of possible bonds (e.g. one site unbonded),
            this will not be rounded/respected in any way
        - if the system contains solvent or other molecules that should not be binding to 
            cross-linkers, make sure to remove them before calling this function
        - useage of mehp.computeCrosslinkerConversion is recommended instead for more consistent results (no dependency on strandLength)

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
    numCrosslinkers = len(network.getAtomsOfType(crosslinkerType))

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
    #   r = calculateEffectiveCrosslinkerFunctionality(network, crosslinkerType, f)
    return math.sqrt(1/(r*(f-1)*(g-1)))
