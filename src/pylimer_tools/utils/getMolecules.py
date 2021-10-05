import pandas as pd


def getAndFilterMoleculesAndBonds(atom_data: pd.DataFrame, bond_data: pd.DataFrame, boxDimensions: list, crosslinker_type: int = None, recalculate_positions=True):
    """
    Find molecules from atom & bond data. Crosslinkers are omitted.

    Args:
        atom_data (pd.DataFrame): a collection of atom coordinates
        bond_data (pd.DataFrame): a collection of bonds
        boxDimensions: a list containing the box lengths (x, y, z) 
        crosslinker_type: the type of the crosslinker atoms to distinguish different molecules/chains

    Returns: 
        molecules (list): a list of lists of atoms representing molecules/chains
        atom_data (pd.DataFrame): the atoms, excluding crosslinkers
        bond_data: the bonds, excluding the ones of the crosslinker atoms
        crosslinker_atoms: the atoms of the crosslinkers
        crosslinker_bonds: the bonds of the crosslinker atoms
    """
    # first, need to split the data per molecule
    molecules = []
    atoms = {}
    crosslinker_atoms = []
    crosslinker_bonds = []
    assert(len(atom_data) > 1)
    if (crosslinker_type is not None):
        lenAtomsBefore = len(atom_data)
        lenBondsBefore = len(bond_data)
        # remove all crosslinkers from bond & atoms
        crosslinkers = atom_data[atom_data.type == crosslinker_type]
        crosslinker_atoms = crosslinkers
        crosslinker_bonds = bond_data[(
            bond_data["to"].isin(crosslinkers["id"]))]
        bond_data = bond_data[~(bond_data["to"].isin(crosslinkers["id"]))]
        crosslinker_bonds = crosslinker_bonds.merge(
            bond_data[(bond_data["bondFrom"].isin(crosslinkers["id"]))])
        bond_data = bond_data[~(bond_data["bondFrom"].isin(crosslinkers["id"]))]
        atom_data = atom_data[~(atom_data["id"].isin(crosslinkers['id']))]
        print("Removed crosslinkers. Got {} of {} atoms and {} of {} bonds left.".format(
            len(atom_data), lenAtomsBefore, len(bond_data), lenBondsBefore))

    # start with first atom, loop all
    atom_data = atom_data.sort_values(
        by="id", axis=0)
    for index, row in atom_data.iterrows():
        # find bonds "from" and "to" the current atom
        bondDf = bond_data
        bonds_from = bondDf.loc[bondDf['bondFrom'] == row['id']]
        bonds_to = bondDf.loc[bondDf['to'] == row['id']]
        atoms[row["id"]] = {
            "atom": row,
            "bonds_from": bonds_from,
            "bonds_to": bonds_to
        }

    atomKeysList = list(atoms.keys())
    for atomId in atomKeysList:
        if (atomId in atoms.keys()):
            molecule = []
            # use this image:
            startAtom = atoms[atomId]["atom"]
            periodic_image = [startAtom["nx"],
                              startAtom["ny"], startAtom["nz"]]
            # follow the chain in both directions
            atomsToFollow = [atoms.pop(atomId)]  # atoms[atomId]
            while(len(atomsToFollow) > 0):
                atomToFollow = atomsToFollow.pop(0)
                # first, unwrap this atom's coordinates
                if (recalculate_positions):
                    dirsToCheck = ["x", "y", "z"]
                    for idx, dir in enumerate(dirsToCheck):
                        dirStr = "n"+dir
                        if (atomToFollow["atom"][dirStr] != periodic_image[idx]):
                            atomToFollow["atom"][dir] -= (
                                periodic_image[idx] - atomToFollow["atom"][dirStr])*boxDimensions[idx]
                # continue following
                for toId in atomToFollow["bonds_from"]["to"]:
                    if (toId in atoms):
                        toAtom = atoms.pop(toId)
                        atomsToFollow.append(toAtom)
                for fromId in atomToFollow["bonds_to"]["bondFrom"]:
                    if (fromId in atoms):
                        toAtom = atoms.pop(fromId)
                        atomsToFollow.append(toAtom)
                molecule.append(atomToFollow)
            molecules.append(molecule)

    return molecules, atom_data, bond_data, crosslinker_atoms, crosslinker_bonds


def getFilteredMoleculesAndBonds(atom_data: pd.DataFrame, bond_data: pd.DataFrame, boxDimensions: list, crosslinker_type: int = None):
    """
    Find molecules from atom & bond data. Crosslinkers are omitted.

    Args:
        atom_data (pd.DataFrame): a collection of atom coordinates
        bond_data (pd.DataFrame): a collection of bonds
        boxDimensions: a list containing the box lengths (x, y, z) 
        crosslinker_type: the type of the crosslinker atoms to distinguish different molecules/chains

    Returns: 
        molecules (list): a list of lists of atoms representing molecules/chains
        atom_data (pd.DataFrame): the atoms, excluding crosslinkers
        bond_data: the bonds, excluding the ones of the crosslinker atoms
    """
    molecules, atom_data, bond_data, _, _ = getAndFilterMoleculesAndBonds(
        atom_data, bond_data, boxDimensions, crosslinker_type)
    return molecules, atom_data, bond_data


def getMolecules(atom_data: pd.DataFrame, bond_data: pd.DataFrame, boxDimensions: list, crosslinker_type: int = None):
    """
    Find molecules from atom & bond data. Crosslinkers are omitted.

    Args:
        atom_data (pd.DataFrame): a collection of atom coordinates
        bond_data (pd.DataFrame): a collection of bonds
        boxDimensions: a list containing the box lengths (x, y, z) 
        crosslinker_type: the type of the crosslinker atoms to distinguish different molecules/chains

    Returns: 
        molecules (list): a list of lists of atoms representing molecules/chains
    """
    molecules, _, _ = getFilteredMoleculesAndBonds(
        atom_data, bond_data, boxDimensions, crosslinker_type)
    return molecules
