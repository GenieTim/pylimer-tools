Module pylimer_tools.io.readLammpDump
=====================================

Functions
---------

    
`readLammpDump(file: str, useCache: bool = True)`
:   Read a lammpstrj / LAMMPS dump output file to a dictonary
    
    Arguments:
        - file: the file path to read from
        - useCache: whether to use the cache or not. Cache does respect file modification time.