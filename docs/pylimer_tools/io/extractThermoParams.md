Module pylimer_tools.io.extractThermoParams
===========================================

Functions
---------

    
`extractThermoParams(file, header='Temp PotEng TotEng Press Volume c_3', textsToRead=5, minLineLen=5, useCache=True) ‑> pandas.core.frame.DataFrame`
:   Extract the thermodynamic outputs produced for this simulation.
    
    Note: the header parameter can be an array — make sure to pay attention
    when reading a file with different header sections in them
    
    Arguments:
        file: the file path to the file to read from
        header: the header of the CSV (where to start reading at)
        textsToRead: the number of times to expect the header
        minLineLen: the minimal length of a line to be accepted as data
        useCache: wheter to use cache or not
    
    Returns: 
        data (pd.DataFrame): the thermodynamic parameters

    
`readOneGroup(fp, header, minLineLen=4, additional_lines_skip=0) ‑> str`
:   Read one group of csv lines from the file
    
    Arguments:
        fp: the file pointer to the file to read from
        header: the header of the CSV (where to start reading at)
        minLineLen: the minimal length of a line to be accepted as data
        additional_lines_skip: number of lines to skip after reading the header