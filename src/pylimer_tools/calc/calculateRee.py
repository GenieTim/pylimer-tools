
import datetime
import hashlib
import os
import pathlib
import pickle
import warnings

import numpy as np

from pylimer_tools.calc.calculateDistanceBetweenAtoms import \
    calculateNormalizedDistanceBetweenAtoms
from pylimer_tools.utils.getMolecules import getFilteredMoleculesAndBonds
from pylimer_tools.io.readLammpData import readLammpData


def calculateRee2AvgForMolecules(molecules, atoms, bonds, boxLengths: list, expected_num_bonds: int = 2):
    """
    Calculate `<Ree^2>`

    Args: 
        molecules (List): a list of molecules as produced by pylimer_tools.utils.getMolecules
        atoms (pd.DataFrame): a collection of atom coordinates
        bonds (pd.DataFrame): a collection of bonds
        boxLenghts: a list containing the box lengths (x, y, z) 

    Returns: 
        retVal: `<Ree^2>`
        Ree: Ree for each molecule
    """
    Ree = []
    chainLen = []
    for chain in molecules:
        # first, find the outermost (non-connected), ends
        undefinedEnd = []
        chainLen.append(len(chain))
        # print("Ree len: {}".format(len(Ree)))
        for atom in chain:
            # print(atom)
            bondsTo = bonds[bonds['to'] == atom['atom']['id']]
            bondsFrom = bonds[bonds['bondFrom'] == atom['atom']['id']]
            # remove the atoms that are in both
            bondsTo = bondsTo[~(bondsTo['bondFrom'].isin(bondsFrom['to']))]
            # calculate length
            lenTo = len(bondsTo)
            lenFrom = len(bondsFrom)
            if (lenTo + lenFrom < expected_num_bonds or lenTo == 0 or lenFrom == 0):
                undefinedEnd.append(atom)
            elif (lenTo + lenFrom > expected_num_bonds):
                warnings.warn("Got atom with more than expected number of bonds: {} instead of {}".format(
                    lenTo+lenFrom, expected_num_bonds))
        if (len(undefinedEnd) != 2):
            warnings.warn(
                "Chain found with wrong ends: {} ends found".format(len(undefinedEnd)))
            continue
        else:
            # then, calculate the distance
            lRee = calculateNormalizedDistanceBetweenAtoms(
                undefinedEnd[0]['atom'], undefinedEnd[1]['atom'], boxLengths)
            # print("Calculated Ree: {}".format(lRee))
            Ree.append(lRee)

    retVal = np.square(Ree).mean()
    # print("Calculated {} Ree = {} ± {} for mean chain length of {} monomers".format(
        # len(Ree), retVal, np.square(Ree).std(), np.mean(chainLen)))
    if (np.mean(chainLen) < np.mean(Ree)):
        warnings.warn("This result is probably wrong/not realistic")
    return retVal, Ree


def calculateRee2Avg(atoms, bonds, boxDimensions: list, crosslinker_type: int = 2, expected_num_bonds: int = 2) -> float:
    """
    Calculate `<Ree^2>`

    Args: 
        atoms (pd.DataFrame): a collection of atom coordinates
        bonds (pd.DataFrame): a collection of bonds
        boxDimensions: a list containing the box lengths (x, y, z) 
        crosslinker_type: the type id of the crosslinker (in order to detect chains, filter those)
        expected_num_bonds: the expected number of bonds. Used for warnings.

    Returns: 
        retVal: `<Ree^2>`
    """
    molecules, atoms, bonds = getFilteredMoleculesAndBonds(
        atoms, bonds, boxDimensions, crosslinker_type)
    print("Got {} chains".format(len(molecules)))

    return calculateRee2AvgForMolecules(molecules, atoms, bonds, boxDimensions, expected_num_bonds)


def calculateRee2AvgForFile(file: str, fileType: str = "data", useCache: bool = True) -> float:
    """
    Calculate `<Ree^2>` for a specific file

    Args: 
        file: the path to the file to calculate Ree for
        fileType: the type of file to read (data or dump)
        useCache: wheter to use cache or not

    Returns: 
        retVal: `<Ree^2>`
    """
    # Custom Cache Stuff: keep for legacy purposes
    cacheFileName = os.path.dirname(
        __file__) + "/cache/" + hashlib.md5(file.encode()).hexdigest() + "-" + fileType + "-ree2-cache.pickle"

    if (os.path.isfile(cacheFileName) and useCache):
        mtimeCache = datetime.datetime.fromtimestamp(
            pathlib.Path(cacheFileName).stat().st_mtime)
        mtimeOrigin = datetime.datetime.fromtimestamp(
            pathlib.Path(file).stat().st_mtime)
        if (mtimeCache > mtimeOrigin):
            with open(cacheFileName, 'rb') as cacheFile:
                loadedData = pickle.load(cacheFile)
                if (isinstance(loadedData, float) or isinstance(loadedData, np.float64) or len(loadedData) < 2):
                    return loadedData
                else:
                    return np.square(loadedData).mean()
        else:
            print("Dump cache file is elder than dump. Reloading...")

    # new, revised code
    Ree = calculateReeForFile(file, fileType, useCache)
    return np.square(Ree).mean()


def calculateReeForFile(file: str, fileType: str = "data", useCache: bool = True) -> np.array:
    """
    Calculate Ree's for a specific file

    Args: 
        file: the path to the file to calculate Ree for
        fileType: the type of file to read (data or dump)
        useCache: wheter to use cache or not

    Returns: 
        retVal: <Ree^2>
    """
    # Custom Cache Stuff
    cacheFileName = os.path.dirname(
        __file__) + "/cache/" + hashlib.md5(file.encode()).hexdigest() + "-" + fileType + "-ree2-cache.pickle"

    if (os.path.isfile(cacheFileName) and useCache):
        mtimeCache = datetime.datetime.fromtimestamp(
            pathlib.Path(cacheFileName).stat().st_mtime)
        mtimeOrigin = datetime.datetime.fromtimestamp(
            pathlib.Path(file).stat().st_mtime)
        if (mtimeCache > mtimeOrigin):
            with open(cacheFileName, 'rb') as cacheFile:
                loadedData = pickle.load(cacheFile)
                if (isinstance(loadedData, float) or isinstance(loadedData, np.float64) or len(loadedData) < 2):
                    # nothing to do: the cache is legacy -> overwrite with new one
                    pass
                else:
                    return np.square(loadedData).mean()
        else:
            print("Dump cache file is elder than dump. Reloading...")

    if (fileType == "data"):
        fileData = readLammpData(file)
        retVal, Ree = calculateRee2Avg(fileData["atom_data"], fileData["bond_data"], [
            fileData["Lx"], fileData["Ly"], fileData["Lz"]])
        with open(cacheFileName, 'wb') as cacheFile:
            pickle.dump(Ree, cacheFile)
        return Ree
    else:
        raise Exception
