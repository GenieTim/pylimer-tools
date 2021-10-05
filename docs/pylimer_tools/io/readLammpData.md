Module pylimer_tools.io.readLammpData
=====================================

Functions
---------

    
`findFloatsInLine(line: str) ‑> list`
:   Find as many floats as possible in a string (line)

    
`findIntsInLine(line: str) ‑> list`
:   Find as many integers as possible in a string (line)

    
`readLammpData(file, useCache=True) ‑> dict`
:   Read a lammpstrj input data file
    
    Args: 
        file: the path to the file to read from
        useCache: wheter to use cache or not (respects modified date of file)