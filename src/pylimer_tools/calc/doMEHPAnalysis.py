# source: https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

import numpy as np
from pylimer_tools.entities.universum import Universum


def calculateCycleRank(network: Universum, nu: int = None, mu: int = None, absTol: float = 1, relTol: float = 1):
    """
    Compute the cycle rank (\(\chi\)) 

    Arguments:
      - network: the network to calculate the cycle rank for
      - nu: number of elastically effective (active) strands per unit volume
      - mu: number density of the elastically effective crosslink
    """
    if (nu is None):
        nu = calculateEffectiveNrDensityOfNetwork(network, absTol, relTol)
    if (mu is None):
        mu = calculateEffectiveNrDensityOfJunctions(network, absTol, relTol)

    return nu - mu


def calculateEffectiveNrDensityOfNetwork(network: Universum, absTol: float = 1, relTol: float = 1):
    """
    Compute the effective number density \(\nu_{eff}\) of a network

    \(\nu_{eff}\) is the number of elastically effective (active) strands per unit volume, 
    which are defined as the ones that can store elastic energy 
    upon network deformation, resp. the effective number density of network strands

    Arguments:
      - network (pylimer_tools.entities.Universum): the network to compute \(\nu_{eff}\) for
      - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol)
      - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol))
    """
    R_taus = []
    chainLengths = []
    for molecule in network.getMolecules():
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


def calculateEffectiveNrDensityOfJunctions(network: Universum, absTol: float = 0, relTol: float = 0):
    """
    
    """
    pass


def calculateTopologicalFactor(network: Universum, foreignAtomType=None, totalMass=1):
    """
    Compute the topological factor of a polymer network.

    Source:
      - eq. 16 in https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262

    Arguments:
      - network: the 
      - foreignAtomType: the type of atoms to ignore
      - totalMass: the \(M\) in the respective formula

    Returns:
      - the topological factor \(\Gamma\)
    """
    molecules = network.getMolecules(ignoreAtomType=foreignAtomType)
    Gamma_sum = 0
    for molecule in molecules:
        R_tau = molecule.computeEndToEndDistance()
        b = molecule.computeBondLengths().mean()
        Gamma_sum += R_tau*R_tau / (molecule.getLength() * b * b)

    return Gamma_sum / totalMass
