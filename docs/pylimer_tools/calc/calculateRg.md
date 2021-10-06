Module pylimer_tools.calc.calculateRg
=====================================

Functions
---------

    
`calculateRg(x, y, z)`
:   Calculate `<R_g^2>`
    
    Arguments:
        x (float): the x coordinate of the atom
        y (float): the y coordinate of the atom
        z (float): the z coordinate of the atom
    
    Returns:
        retVal: `<R_g^2>`

    
`calculateRg2Avg(atoms, bonds, boxDimensions: list, crosslinker_type: int = 2, expected_num_bonds: int = 2) ‑> float`
:   Calculate `<R_g^2>`
    
    Arguments:
        - atoms (pd.DataFrame): a collection of atom coordinates
        - bonds (pd.DataFrame): a collection of bonds
        - boxDimensions: a list containing the box lengths (x, y, z) 
        - crosslinker_type: the type id of the crosslinker (in order to detect chains, filter those)
        - expected_num_bonds: the expected number of bonds. Used for warnings.
    
    Returns:
        - retVal: `<R_g^2>`

    
`calculateRg2AvgForFile(file: str, fileType: str = 'data', useCache: bool = True) ‑> float`
:   Calculate `<Rg^2>` for a specific file
    
    Arguments:
        - file: the path to the file to calculate Ree for
        - fileType: the type of file to read (data or dump)
        - useCache: wheter to use cache or not
    
    Returns:
        - retVal: `<R_g^2>`

    
`calculateRg2AvgForMolecules(molecules, atoms, bonds, boxLengths: list, expected_num_bonds: int = 2)`
:   Calculate `<R_g^2>`
    
    Arguments:
        molecules (list): list of molecules to calculate the R_g for
        atoms (pd.DataFrame): a collection of atom coordinates
        bonds (pd.DataFrame): a collection of bonds
        boxLengths: a list containing the box lengths (x, y, z) 
        expected_num_bonds: the expected number of bonds. Used for warnings.
    
    Returns:
        retVal: `<R_g^2>`

    
`calculateRgForFile(file: str, fileType: str = 'data', useCache: bool = True) ‑> <built-in function array>`
:   Calculate R_g for a specific file
    
    Arguments:
        - file: the path to the file to calculate Ree for
        - fileType: the type of file to read (data or dump)
        - useCache: wheter to use cache or not
    
    Returns:
        - retVal (list): R_gs