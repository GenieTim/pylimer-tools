import math

import numpy as np
import pandas as pd
from pandas.core.algorithms import isin


def normalizeData(atom):
    """ 
    Verify the atom to be a usable one

    Args:
        atom: the atom to check
    
    Returns: 
        atom (pd.Seroes): the resolved atom
    """
    if (isinstance(atom, pd.DataFrame)):
        return atom
    if (isinstance(atom, dict)):
        if ("atom" in atom.keys()):
            return atom['atom']
    return atom


def row_converter(row: pd.DataFrame, listy=None):
    """ 
    convert pandas row to a dictionary

    Args:
        row (pd.DataFrame|pd.Series): row as a tuple
        listy: a list of columns
    
    Returns: 
        pictionary (dictionary): the row's values
    """
    if (listy is None):
        if (isinstance(row, pd.Series)):
            listy = row.index
        else:
            listy = row.columns
    pictionary = {}
    for item in listy:
        pictionary[item] = row[item].values[0]
    return pictionary


def calculateDistanceBetweenAtoms(atomA, atomB):
    """
    Calculate the the distance between two atoms. 
    No translation between periodic images happens.

    Args: 
        atomA: the coordinates of atom 1
        atomB: the coordinates of atom 2

    Returns: 
        meanDistance: the norm of the connecting vector between the two coordinates
    """
    locaA = normalizeData(atomA)
    locaB = normalizeData(atomB)
    return math.sqrt((locaA["x"]-locaB["x"])**2+(locaA["y"]-locaB["y"])**2+(locaA["z"]-locaB["z"])**2)


def calculateNormalizedDistanceBetweenAtoms(atom1, atom2, boxLengths: list):
    """
    Calculate the the distance between two atoms. 

    Args: 
        atom1: the coordinates of atom 1
        atom2: the coordinates of atom 2
        boxLenghts: a list containing the box lengths (x, y, z) 

    Returns: 
        meanDistance: the norm of the connecting vector between the two coordinates
    """
    if (not isinstance(atom1, dict) and isinstance(atom1, pd.DataFrame)):
        atom1 = row_converter(atom1)
    if (not isinstance(atom2, dict) and isinstance(atom2, pd.DataFrame)):
        atom2 = row_converter(atom2)
    atom1 = normalizeData(atom1)
    atom2 = normalizeData(atom2)
    # first, unwrap this atom's coordinates
    dirsToCheck = ["x", "y", "z"]

    deltas = [None, None, None]

    for idx, dir in enumerate(dirsToCheck):
        delta = abs(atom1[dir]-atom2[dir])

        while(delta > boxLengths[idx]*0.5 and delta > 0):
            delta -= boxLengths[idx]

        deltas[idx] = delta

    return np.linalg.norm(deltas)
