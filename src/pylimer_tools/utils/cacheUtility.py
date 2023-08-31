import datetime
import hashlib
import os
import pathlib
import pickle
import tempfile
import warnings
from typing import Iterable, List, Union

import numpy as np


def doCache(obj, file: str, suffix: str, tmp_dir: str = None):
    """
    Store the object in the cache

    Arguments:
        - obj: the thing to cache
        - file: the path to the file to save the object to
        - suffix: the file name's suffix
        - tmp_dir: the directory to store the cache in
    """
    cacheFileName = getCacheFileName(file, suffix, tmp_dir)
    with open(cacheFileName, 'wb') as cacheFile:
        pickle.dump(obj, cacheFile)


def loadCache(file: Union[str, List[str], None], suffix: str, disableWarnings: bool = False, tmp_dir: str = None, anyway: bool = False):
    """
    Load an object from cache.

    Arguments:
        - file: a cache name. Ideally the file that is read, such that the filemtime of `file` can be used to check whether cache must be generated anew
        - suffix: the file name's suffix
        - disableWarnings: whether to disable warnings about missing possibilities to check for filemtime
        - tmp_dir: the directory to load the cache from

    Returns:
        - cache: either the content of the cache, or None if the cache has to be loaded again / is non existant
    """
    if (file is None):
        file = ""
    if (not isinstance(file, list)):
        file = [file]
    cacheFileName = getCacheFileName(file, suffix, tmp_dir)
    if (os.path.isfile(cacheFileName)):
        if (not np.all([os.path.isfile(f) for f in file])):
            if (not disableWarnings):
                warnings.warn(
                    'Cache called for non-existent file. Make sure the key is time-restricted')
            with open(cacheFileName, 'rb') as cacheFile:
                toReturn = pickle.load(cacheFile)
            return toReturn
        else:
            mtimeCache = datetime.datetime.fromtimestamp(
                pathlib.Path(cacheFileName).stat().st_mtime)
            mtimesOrigin = [datetime.datetime.fromtimestamp(
                pathlib.Path(f).stat().st_mtime) for f in file]
            if (np.all(mtimeCache > np.array(mtimesOrigin)) or anyway):
                try:
                    with open(cacheFileName, 'rb') as cacheFile:
                        toReturn = pickle.load(cacheFile)
                    return toReturn
                except pickle.UnpicklingError as e:
                    warnings.warn(
                        "Unpickling of cache file {} failed: {}".format(file, e))
                except ModuleNotFoundError as e:
                    warnings.warn(
                        "Unpickling of cache file {} failed: {}".format(file, e))
            else:
                # print("Dump cache file is elder than dump. Reloading...")
                pass

    return None


def getCacheFileName(file: Union[str, List[str], None], suffix: str, tmp_dir: str = None):
    """
    Get the name and path of a cache file. Internal method.

    Arguments:
        - file: a cache name. Ideally the file that is read.
        - suffix: the file name's suffix
        - tmp_dir: the temporary directory

    Returns:
        - cacheFileName: the path to the cache file
    """
    if (isinstance(file, list)):
        file = "".join(file)
    if (file is None):
        file = ""
    cacheFileName = "{}/{}-{}.pickle".format(
        tempfile.gettempdir() if tmp_dir is None else tmp_dir,
        hashlib.md5(file.encode()).hexdigest(), suffix)
    return cacheFileName
