
import datetime
import hashlib
import os
import pathlib
import pickle
import warnings

import numpy as np

from pylimer_tools.utils.getMolecules import getFilteredMoleculesAndBonds
from pylimer_tools.io.readLammpData import readLammpData


def calculateRg(x, y, z):
    """
    Calculate `<R_g^2>`

    Args: 
        x (float): the x coordinate of the atom
        y (float): the y coordinate of the atom
        z (float): the z coordinate of the atom

    Returns: 
        retVal: `<R_g^2>`
    """
    if (not isinstance(x, np.ndarray)):
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
    Rxyz = np.matrix([(x - x.mean()), (y - y.mean()), (z - z.mean())])
    return sum(np.linalg.norm(Rxyz, 2, axis=0)**2)*(1/len(x))


def calculateRg2AvgForMolecules(molecules, atoms, bonds, boxLengths: list, expected_num_bonds: int = 2):
    """
    Calculate `<R_g^2>`

    Args: 
        molecules (list): list of molecules to calculate the R_g for
        atoms (pd.DataFrame): a collection of atom coordinates
        bonds (pd.DataFrame): a collection of bonds
        boxLengths: a list containing the box lengths (x, y, z) 
        expected_num_bonds: the expected number of bonds. Used for warnings.

    Returns: 
        retVal: `<R_g^2>`
    """
    Rg = []
    chainLen = []
    for chain in molecules:
        # first, find the coordinates of all atoms
        xs = []
        ys = []
        zs = []
        chainLen.append(len(chain))
        # print("Rg len: {}".format(len(Rg)))
        for atom in chain:
            # print(atom)
            # atomHere = atoms[atoms['id'] == atom['atom']['id']]
            # assert(len(atomHere) == 1)
            # atomHere = atomHere.iloc[0]
            atomHere = atom['atom']
            xs.append(atomHere['x'])
            ys.append(atomHere['y'])
            zs.append(atomHere['z'])
        if (len(xs)):
            # then, calculate the distance
            lRg = calculateRg(xs, ys, zs)
            # print("Calculated Rg: {}".format(lRg))
            Rg.append(lRg)

    retVal = np.mean(Rg)
    print("Calculated {} Rg = {} ± {} for mean chain length of {} monomers".format(
        len(Rg), retVal, np.std(Rg), np.mean(chainLen)))
    if (np.mean(chainLen) < retVal):
        warnings.warn("This result is probably wrong/not realistic")
    return retVal, Rg


def calculateRg2Avg(atoms, bonds, boxDimensions: list, crosslinker_type: int = 2, expected_num_bonds: int = 2) -> float:
    """
    Calculate `<R_g^2>`

    Args: 
        atoms (pd.DataFrame): a collection of atom coordinates
        bonds (pd.DataFrame): a collection of bonds
        boxDimensions: a list containing the box lengths (x, y, z) 
        crosslinker_type: the type id of the crosslinker (in order to detect chains, filter those)
        expected_num_bonds: the expected number of bonds. Used for warnings.

    Returns: 
        retVal: `<R_g^2>`
    """
    molecules, atoms, bonds = getFilteredMoleculesAndBonds(
        atoms, bonds, boxDimensions, crosslinker_type)
    print("Got {} chains".format(len(molecules)))

    return calculateRg2AvgForMolecules(molecules, atoms, bonds, boxDimensions, expected_num_bonds)


def calculateRg2AvgForFile(file: str, fileType: str = "data", useCache: bool = True) -> float:
    """
    Calculate `<Rg^2>` for a specific file

    Args: 
        file: the path to the file to calculate Ree for
        fileType: the type of file to read (data or dump)
        useCache: wheter to use cache or not

    Returns: 
        retVal: `<R_g^2>`
    """
    # Custom Cache Stuff: keep for legacy purposes
    cacheFileName = os.path.dirname(
        __file__) + "/cache/" + hashlib.md5(file.encode()).hexdigest() + "-" + fileType + "-Rg2-cache.pickle"

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
    Rg = calculateRgForFile(file, fileType, useCache)
    return np.square(Rg).mean()


def calculateRgForFile(file: str, fileType: str = "data", useCache: bool = True) -> np.array:
    """
    Calculate R_g for a specific file

    Args: 
        file: the path to the file to calculate Ree for
        fileType: the type of file to read (data or dump)
        useCache: wheter to use cache or not

    Returns: 
        retVal (list): R_gs
    """
    # Custom Cache Stuff
    cacheFileName = os.path.dirname(
        __file__) + "/cache/" + hashlib.md5(file.encode()).hexdigest() + "-" + fileType + "-Rg2-cache.pickle"

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
                    return loadedData
        else:
            print("Dump cache file is elder than dump. Reloading...")

    if (fileType == "data"):
        fileData = readLammpData(file)
        retVal, Rg = calculateRg2Avg(fileData["atom_data"], fileData["bond_data"], [
            fileData["Lx"], fileData["Ly"], fileData["Lz"]])
        with open(cacheFileName, 'wb') as cacheFile:
            pickle.dump(Rg, cacheFile)
        return Rg
    else:
        raise Exception
