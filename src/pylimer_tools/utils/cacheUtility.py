import datetime
import hashlib
import os
import pathlib
import pickle
import tempfile
import warnings


def doCache(obj, file: str, suffix: str):
    """
    Store the object in the cache

    Arguments:
        - obj: the object to cache
        - file: the path to the file to save the object to
        - suffix: the file name's suffix
    """
    cacheFileName = getCacheFileName(file, suffix)
    with open(cacheFileName, 'wb') as cacheFile:
        pickle.dump(obj, cacheFile)


def loadCache(file: str, suffix: str, disableWarnings: bool = False):
    """
    Load an object from cache.

    Arguments:
        - file: a cache name. Ideally the file that is read, such that the filemtime of `file` can be used to check whether cache must be generated anew
        - suffix: the file name's suffix
        - disableWarnings: whether to disable warnings about missing possibilities to check for filemtime

    Returns:
        - cache: either the content of the cache, or None if the cache has to be loaded again / is non existant
    """
    cacheFileName = getCacheFileName(file, suffix)
    if (os.path.isfile(cacheFileName)):
        if (not os.path.isfile(file)):
            if (not disableWarnings):
                warnings.warn(
                    'Cache called for non-existent file. Make sure the key is time-restricted')
            with open(cacheFileName, 'rb') as cacheFile:
                toReturn = pickle.load(cacheFile)
            return toReturn
        else:
            mtimeCache = datetime.datetime.fromtimestamp(
                pathlib.Path(cacheFileName).stat().st_mtime)
            mtimeOrigin = datetime.datetime.fromtimestamp(
                pathlib.Path(file).stat().st_mtime)
            if (mtimeCache > mtimeOrigin):
                with open(cacheFileName, 'rb') as cacheFile:
                    toReturn = pickle.load(cacheFile)
                return toReturn
            else:
                # print("Dump cache file is elder than dump. Reloading...")
                pass

    return None


def getCacheFileName(file: str, suffix: str):
    """
    Get the name and path of a cache file. Internal method.

    Arguments:
        - file: a cache name. Ideally the file that is read.
        - suffix: the file name's suffix

    Returns:
        - cacheFileName: the path to the cache file
    """
    cacheFileName = "{}/{}-{}.pickle".format(
        tempfile.gettempdir(),
        hashlib.md5(file.encode()).hexdigest(), suffix)
    return cacheFileName
