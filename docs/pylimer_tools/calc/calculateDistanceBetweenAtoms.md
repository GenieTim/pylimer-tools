Module pylimer_tools.calc.calculateDistanceBetweenAtoms
=======================================================

Functions
---------

    
`calculateDistanceBetweenAtoms(atomA, atomB)`
:   Calculate the the distance between two atoms. 
    No translation between periodic images happens.
    
    Args: 
        atomA: the coordinates of atom 1
        atomB: the coordinates of atom 2
    
    Returns: 
        meanDistance: the norm of the connecting vector between the two coordinates

    
`calculateNormalizedDistanceBetweenAtoms(atom1, atom2, boxLengths: list)`
:   Calculate the the distance between two atoms. 
    
    Args: 
        atom1: the coordinates of atom 1
        atom2: the coordinates of atom 2
        boxLenghts: a list containing the box lengths (x, y, z) 
    
    Returns: 
        meanDistance: the norm of the connecting vector between the two coordinates

    
`normalizeData(atom)`
:   Verify the atom to be a usable one
    
    Args:
        atom: the atom to check
    
    Returns: 
        atom (pd.Seroes): the resolved atom

    
`row_converter(row: pandas.core.frame.DataFrame, listy=None)`
:   convert pandas row to a dictionary
    
    Args:
        row (pd.DataFrame|pd.Series): row as a tuple
        listy: a list of columns
    
    Returns: 
        pictionary (dictionary): the row's values