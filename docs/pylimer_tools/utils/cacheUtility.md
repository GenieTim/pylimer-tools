Module pylimer_tools.utils.cacheUtility
=======================================

Functions
---------

    
`doCache(obj, file: str, suffix: str)`
:   Store the object in the cache
    
    Arguments:
        - obj: the object to cache
        - file: the path to the file to save the object to
        - suffix: the file name's suffix

    
`getCacheFileName(file: str, suffix: str)`
:   Get the name and path of a cache file. Internal method.
    
    Arguments:
        - file: a cache name. Ideally the file that is read.
        - suffix: the file name's suffix
    
    Returns:
        - cacheFileName: the path to the cache file

    
`loadCache(file: str, suffix: str, disableWarnings: bool = False)`
:   Load an object from cache.
    
    Arguments:
        - file: a cache name. Ideally the file that is read, such that the filemtime of `file` can be used to check whether cache must be generated anew
        - suffix: the file name's suffix
        - disableWarnings: whether to disable warnings about missing possibilities to check for filemtime
    
    Returns:
        - cache: either the content of the cache, or None if the cache has to be loaded again / is non existant