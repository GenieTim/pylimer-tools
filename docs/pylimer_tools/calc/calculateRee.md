Module pylimer_tools.calc.calculateRee
======================================

Functions
---------

    
`calculateRee2Avg(atoms, bonds, boxDimensions: list, crosslinker_type: int = 2, expected_num_bonds: int = 2) ‑> float`
:   Calculate `<Ree^2>`
    
    Arguments:
        atoms (pd.DataFrame): a collection of atom coordinates
        bonds (pd.DataFrame): a collection of bonds
        boxDimensions: a list containing the box lengths (x, y, z) 
        crosslinker_type: the type id of the crosslinker (in order to detect chains, filter those)
        expected_num_bonds: the expected number of bonds. Used for warnings.
    
    Returns:
        retVal: `<Ree^2>`

    
`calculateRee2AvgForFile(file: str, fileType: str = 'data', useCache: bool = True) ‑> float`
:   Calculate `<Ree^2>` for a specific file
    
    Arguments:
        file: the path to the file to calculate Ree for
        fileType: the type of file to read (data or dump)
        useCache: wheter to use cache or not
    
    Returns:
        retVal: `<Ree^2>`

    
`calculateRee2AvgForMolecules(molecules, atoms, bonds, boxLengths: list, expected_num_bonds: int = 2)`
:   Calculate `<Ree^2>`
    
    Arguments:
        - molecules (List): a list of molecules as produced by pylimer_tools.utils.getMolecules
        - atoms (pd.DataFrame): a collection of atom coordinates
        - bonds (pd.DataFrame): a collection of bonds
        - boxLenghts: a list containing the box lengths (x, y, z) 
    
    Returns:
        - retVal: `<Ree^2>`
        - Ree: Ree for each molecule

    
`calculateReeForFile(file: str, fileType: str = 'data', useCache: bool = True) ‑> <built-in function array>`
:   Calculate Ree's for a specific file
    
    Arguments:
        file: the path to the file to calculate Ree for
        fileType: the type of file to read (data or dump)
        useCache: wheter to use cache or not
    
    Returns:
        retVal: <Ree^2>