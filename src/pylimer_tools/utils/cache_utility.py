import datetime
import hashlib
import os
import pickle
import shutil
import tempfile
import warnings
from typing import List, Union


def do_cache(obj, file: str, suffix: str, tmp_dir: Union[str, None] = None):
    """
    Store the object in the cache.

    :param obj: The object to cache.
    :param file: A part of what's used for the cache's name. Ideally the file that is read,
        such that the filemtime of `file` can be used to check whether cache must be generated anew.
    :param suffix: The file name's suffix.
    :param tmp_dir: The directory to store the cache in.
    """
    cache_file_name = get_cache_file_name(file, suffix, tmp_dir)
    with open(cache_file_name, "wb") as cache_file:
        pickle.dump(obj, cache_file)


def load_cache(
    file: Union[str, List[str], None],
    suffix: str,
    disable_warnings: bool = False,
    tmp_dir: Union[str, None] = None,
    anyway: bool = False,
):
    """
    Load an object from cache, iff the cache is new enough.

    :param file: A part of what's used for the cache's name. Ideally the file that is read,
        such that the filemtime of `file` can be used to check whether cache must be generated anew.
    :param suffix: The file name's suffix.
    :param disable_warnings: Whether to disable warnings about missing possibilities to check for filemtime.
    :param tmp_dir: The directory to load the cache from.
    :param anyway: Whether to ignore the cache's modification time, and return the cached data anyway,
        as if it were current.
    :return: Either the content of the cache, or None if the cache has to be loaded again / is non existent.
    """
    if file is None:
        file = ""
    if not isinstance(file, list):
        file = [file]
    cache_file_name = get_cache_file_name(file, suffix, tmp_dir)
    old_cache_file_name = get_cache_file_name(file, suffix, tmp_dir, True)
    if os.path.isfile(old_cache_file_name) and not os.path.isfile(
            cache_file_name):
        shutil.copy2(old_cache_file_name, cache_file_name)
        os.remove(old_cache_file_name)
    if os.path.isfile(cache_file_name):
        if not all(os.path.isfile(f) for f in file):
            if not disable_warnings:
                warnings.warn(
                    "Cache called for non-existent file. Make sure the key is time-restricted"
                )
            with open(cache_file_name, "rb") as cache_file:
                to_return = pickle.load(cache_file)
            return to_return
        else:
            if anyway or is_current_cache(cache_file_name, file):
                try:
                    with open(cache_file_name, "rb") as cache_file:
                        to_return = pickle.load(cache_file)
                    return to_return
                except pickle.UnpicklingError as e:
                    warnings.warn(
                        "Unpickling of cache file {} failed: {}".format(
                            file, e)
                    )
                except ModuleNotFoundError as e:
                    warnings.warn(
                        "Unpickling of cache file {} failed: {}".format(
                            file, e)
                    )
            else:
                # print("Dump cache file is elder than dump. Reloading...")
                pass

    return None


def get_cache_file_name(
    file: Union[str, List[str], None],
    suffix: str,
    tmp_dir: Union[str, None] = None,
    old: bool = False,
):
    """
    Get the name and path of a cache file. Internal method.

    :param file: A part of what's used for the cache's name. Ideally the file that is read,
        such that the filemtime of `file` can be used to check whether cache must be generated anew.
    :param suffix: The file name's suffix.
    :param tmp_dir: The temporary directory.
    :param old: Whether to use the old file naming scheme.
    :return: The path to the cache file.
    """
    if isinstance(file, list):
        file = "".join(sorted(file) if not old else file)
    if file is None:
        file = ""
    cache_file_name = "{}/{}-{}.pickle".format(
        tempfile.gettempdir() if tmp_dir is None else tmp_dir,
        hashlib.md5(file.encode()).hexdigest(),
        suffix,
    )
    return cache_file_name


def is_current_cache(cache_file: str, dependencies: Union[str, List[str]]):
    """
    Determine whether the provided file is newer than all its dependencies.

    :param cache_file: The cache file that is required to be newer.
    :param dependencies: The list of files (or a single file path) that need to be older.
    :return: True if the file is newer than all its dependencies, False otherwise.
    """
    if not os.path.exists(cache_file):
        return False
    if not isinstance(dependencies, list):
        dependencies = [dependencies]
    mtime_cache = datetime.datetime.fromtimestamp(os.path.getmtime(cache_file))
    mtimes_origin = [
        datetime.datetime.fromtimestamp(os.path.getmtime(f)) for f in dependencies
    ]

    return all(mtime_cache > mt for mt in mtimes_origin)
