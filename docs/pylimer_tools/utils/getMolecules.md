Module pylimer_tools.utils.getMolecules
=======================================

Functions
---------

    
`getAndFilterMoleculesAndBonds(atom_data: pandas.core.frame.DataFrame, bond_data: pandas.core.frame.DataFrame, boxDimensions: list, crosslinker_type: int = None, recalculate_positions=True)`
:   Find molecules from atom & bond data. Crosslinkers are omitted.
    
    Arguments:
        atom_data (pd.DataFrame): a collection of atom coordinates
        bond_data (pd.DataFrame): a collection of bonds
        boxDimensions: a list containing the box lengths (x, y, z) 
        crosslinker_type: the type of the crosslinker atoms to distinguish different molecules/chains
    
    Returns:
        - molecules (list): a list of lists of atoms representing molecules/chains
        - atom_data (pd.DataFrame): the atoms, excluding crosslinkers
        - bond_data: the bonds, excluding the ones of the crosslinker atoms
        - crosslinker_atoms: the atoms of the crosslinkers
        - crosslinker_bonds: the bonds of the crosslinker atoms

    
`getFilteredMoleculesAndBonds(atom_data: pandas.core.frame.DataFrame, bond_data: pandas.core.frame.DataFrame, boxDimensions: list, crosslinker_type: int = None)`
:   Find molecules from atom & bond data. Crosslinkers are omitted.
    
    Arguments:
        - atom_data (pd.DataFrame): a collection of atom coordinates
        - bond_data (pd.DataFrame): a collection of bonds
        - boxDimensions: a list containing the box lengths (x, y, z) 
        - crosslinker_type: the type of the crosslinker atoms to distinguish different molecules/chains
    
    Returns:
        - molecules (list): a list of lists of atoms representing molecules/chains
        - atom_data (pd.DataFrame): the atoms, excluding crosslinkers
        - bond_data: the bonds, excluding the ones of the crosslinker atoms

    
`getMolecules(atom_data: pandas.core.frame.DataFrame, bond_data: pandas.core.frame.DataFrame, boxDimensions: list, crosslinker_type: int = None)`
:   Find molecules from atom & bond data. Crosslinkers are omitted.
    
    Arguments:
        - atom_data (pd.DataFrame): a collection of atom coordinates
        - bond_data (pd.DataFrame): a collection of bonds
        - boxDimensions: a list containing the box lengths (x, y, z) 
        - crosslinker_type: the type of the crosslinker atoms to distinguish different molecules/chains
    
    Returns:
        - molecules (list): a list of lists of atoms representing molecules/chains