import numpy as np
import pandas as pd

from pylimer_tools.utils.cacheUtility import doCache, loadCache

def findFloatsInLine(line: str) -> list:
    """
    Find as many floats as possible in a string (line)
    """
    return [float(s) for s in line.split() if (s[0].isdigit() or (len(s) > 0 and s[1].isdigit() and s[0] == "-"))]


def findIntsInLine(line: str) -> list:
    """
    Find as many integers as possible in a string (line)
    """
    return [int(s) for s in line.split() if s.isdigit()]


def readLammpData(file, useCache=True) -> dict:
    """
    Read a lammpstrj input data file

    Args: 
        file: the path to the file to read from
        useCache: wheter to use cache or not (respects modified date of file)
    """
    if (useCache):
        toReturn = loadCache(file, "-lammp-data-cache")
        if (toReturn is not None):
            toReturn['bond_data'].rename(
                columns={"from": "bondFrom"}, inplace=True)
            return toReturn

    out = {}

    with open(file, 'r') as fp:
        N_atoms = 0
        N_bonds = 0
        N_Atypes = 0
        N_Btypes = 0
        Lx_dat = [0, 0]
        Ly_dat = [0, 0]
        Lz_dat = [0, 0]

        line = fp.readline()
        while(line and "Masses" not in line):
            line = fp.readline()
            if ("LAMMPS" in line):
                continue
            if ("atoms" in line):
                N_atoms = findIntsInLine(line)[0]
            elif ("bonds" in line):
                N_bonds = findIntsInLine(line)[0]
            elif ("atom types" in line):
                N_Atypes = findIntsInLine(line)[0]
            elif("bond types" in line):
                N_Btypes = findIntsInLine(line)[0]
            elif("xlo xhi" in line):
                Lx_dat = findFloatsInLine(line)
            elif("ylo yhi" in line):
                Ly_dat = findFloatsInLine(line)
            elif("zlo zhi" in line):
                Lz_dat = findFloatsInLine(line)

        Lx = Lx_dat[1] - Lx_dat[0]
        Ly = Ly_dat[1] - Ly_dat[0]
        Lz = Lz_dat[1] - Lz_dat[0]

        masses = ""
        line = fp.readline()
        while(line and "Atoms" not in line):
            masses += line
            line = fp.readline()

        while(line and line.strip() == ""):
            line = fp.readline()

        atom_data = np.loadtxt(
            fp, dtype={
                'names': ("id", "molecule_id", "type", "x", "y", "z", "nx", "ny", "nz"),
                'formats': ("int", "int", "int", "float", "float", "float", "int", "int", "int")
            }, max_rows=N_atoms)

        while(line and "Bonds" not in line):
            line = fp.readline()

        while(line and line.strip() == ""):
            line = fp.readline()

        bond_data = np.loadtxt(fp, {
            'names': ("id", "type", "bondFrom", "to"),
            'formats': ("int", "int", "int", "int")
        }, max_rows=N_bonds)

        out = {
            "N_atoms": N_atoms,
            "N_Atypes": N_Atypes,
            "N_Btypes": N_Btypes,
            "N_bonds": N_bonds,
            "masses": masses,
            "Lx": Lx, "Ly": Ly, "Lz": Lz,
            "xlo": Lx_dat[0],
            "xhi": Lx_dat[1],
            "ylo": Ly_dat[0],
            "yhi": Ly_dat[1],
            "zlo": Lz_dat[0],
            "zhi": Lz_dat[1],
            "atom_data": pd.DataFrame(atom_data),
            "bond_data": pd.DataFrame(bond_data)
        }

    doCache(out, file, "-lammp-data-cache")

    return out
