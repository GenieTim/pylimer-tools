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

    @deprecated This function is legacy compliant only

    Args: 
        coordsDf: a dataframe containing the coordinates
        boxLenghts: a list containing the box lengths (x, y, z) 

    Returns: 
        meanDistance: the mean of all distances
    """
    lastX, lastY, lastZ = [0, 0, 0]
    newX, newY, newZ = [0, 0, 0]
    distances = []
    minId = coordsDf["id"].min()
    maxId = coordsDf["id"].max()
    for fromId in range(minId, maxId):
        row = coordsDf.loc[coordsDf["id"] == fromId]
        if (lastX == 0 and lastY == 0 and lastZ == 0):
            lastX = row["xsu"].iloc[0]*boxLengths[0]
            lastY = row["ysu"].iloc[0]*boxLengths[1]
            lastZ = row["zsu"].iloc[0]*boxLengths[2]
        else:
            newX = row["xsu"].iloc[0]*boxLengths[0]
            newY = row["ysu"].iloc[0]*boxLengths[1]
            newZ = row["zsu"].iloc[0]*boxLengths[2]
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
    """
    Rs = []
    for bond in bondsDf.itertuples():
        atomTo = coordsDf[coordsDf['id'] == bond.to]
        atomFrom = coordsDf[coordsDf['id'] == bond.bondFrom]
        if (skipAtomType is not None):
            if (atomTo['type'] == skipAtomType or atomFrom['type'] == skipAtomType):
                continue
        Rdist = calculateNormalizedDistanceBetweenAtoms(atomTo, atomFrom, boxLengths)
        if (Rdist > 5):
            warnings.warn("Probably unrealistically long bond detected.")
        Rs.append(Rdist)

    return np.array(Rs)

