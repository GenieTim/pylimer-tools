Module pylimer_tools.calc.calculateBondLen
==========================================

Functions
---------

    
`calculateBondLen(coordsDf: pandas.core.frame.DataFrame, bondsDf: pandas.core.frame.DataFrame, boxLengths: list, skipAtomType=None)`
:   Calculate the bond lengths
    given the coordinates of the atoms in a pd.DataFrame 
    and the bonds of the atoms in a pd.DataFrame 
    and the boxLengths in a list

    
`calculateMeanBondLen(coordsDf: pandas.core.frame.DataFrame, boxLengths: list)`
:   Calculate the mean bond length 
    given the coordinates of the atoms in a pd.DataFrame 
    and the boxLengths in a list
    
    Assumes the bonds are by the coordsDf's ids, sequentially
    
    @deprecated This function is legacy compliant only
    
    Args: 
        coordsDf: a dataframe containing the coordinates
        boxLenghts: a list containing the box lengths (x, y, z) 
    
    Returns: 
        meanDistance: the mean of all distances