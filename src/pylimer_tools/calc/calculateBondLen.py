from calc.calculateDistanceBetweenAtoms import calculateNormalizedDistanceBetweenAtoms
import warnings

import numpy as np
import pandas as pd
from pandas.core.algorithms import isin

from pylimer_tools.utils.getMolecules import getMolecules


def calculateMeanBondLen(coordsDf: pd.DataFrame, boxLengths: list):
    """
    Calculate the mean bond length 
    given the coordinates of the atoms in a pd.DataFrame 
    and the boxLengths in a list

    Assumes the bonds are by the coordsDf's ids, sequentially

    .. deprecated:: 0.0.1 
        This function is legacy compliant only. Use `calculateBondLen(...).mean()` instead.

    Arguments:
        - coordsDf: a dataframe containing the coordinates
        - boxLenghts: a list containing the box lengths (x, y, z) 

    Returns:
        - meanDistance: the mean of all distances
    """
    lastX, lastY, lastZ = [0, 0, 0]
    newX, newY, newZ = [0, 0, 0]
    distances = []
    minId = coordsDf["id"].min()
    maxId = coordsDf["id"].max()
    for fromId in range(minId, maxId):
        row = coordsDf.loc[coordsDf["id"] == fromId]

        def ilocIfNecessary(row, key):
            if (type(row[key]) in (int, str, bool, float, np.float64, np.float32, np.float16, np.float128)):
                return row[key]
            else:
                return row[key].iloc[0]

        if (fromId == minId):
            lastX = ilocIfNecessary(row, "xsu")*boxLengths[0]
            lastY = ilocIfNecessary(row, "ysu")*boxLengths[1]
            lastZ = ilocIfNecessary(row, "zsu")*boxLengths[2]
        else:
            lastX = ilocIfNecessary(row, "xsu")*boxLengths[0]
            lastY = ilocIfNecessary(row, "ysu")*boxLengths[1]
            lastZ = ilocIfNecessary(row, "zsu")*boxLengths[2]
            distance = np.linalg.norm([lastX-newX, lastY-newY, lastZ-newZ])
            distances.append(distance)
            lastX = newX
            lastY = newY
            lastZ = newZ
    distances = np.array(distances)
    return distances.mean()  # tolist()


def calculateBondLen(coordsDf: pd.DataFrame, bondsDf: pd.DataFrame, boxLengths: list, skipAtomType=None):
    """
    Calculate the bond lengths
    given the coordinates of the atoms in a pd.DataFrame 
    and the bonds of the atoms in a pd.DataFrame 
    and the boxLengths in a list

    Arguments:
        - coordsDf: a dataframe containing the coordinates
        - boxLenghts: a list containing the box lengths (x, y, z) 

    Returns:
        - a np.array of all bond lengths
    """
    Rs = []
    for bond in bondsDf.itertuples():
        atomTo = coordsDf[coordsDf['id'] == bond.to]
        atomFrom = coordsDf[coordsDf['id'] == bond.bondFrom]
        if (skipAtomType is not None):
            if (atomTo['type'] == skipAtomType or atomFrom['type'] == skipAtomType):
                continue
        Rdist = calculateNormalizedDistanceBetweenAtoms(
            atomTo, atomFrom, boxLengths)
        if (Rdist > 5):
            warnings.warn("Probably unrealistically long bond detected.")
        Rs.append(Rdist)

    return np.array(Rs)
