
import datetime
import hashlib
import os
import pathlib
import pickle
from io import StringIO

import pandas as pd
from pylimer_tools.utils.cacheUtility import doCache, loadCache


def readLammpDump(file: str, useCache: bool = True):
    """
    Read a lammpstrj / LAMMPS dump output file to a dictonary

    Args:
        file: the file path to read from
        useCache: whether to use the cache or not. Cache does respect file modification time.
    """
    if (useCache):
        cachedData = loadCache(file, "-lammp-dump-cache")
        if (cachedData is not None):
            return cachedData

    out = []
    currentKey = ""
    item = {}
    # read file: pretty slow currently.
    # Many speedup possibilities, e.g. reading multiple lines at once etc.
    with open(file, 'r') as fp:
        line = fp.readline()
        newGroupKey = line
        while line:
            # new header
            if (line.startswith("ITEM:")):
                currentKey = line.replace("ITEM: ", "", 1).strip()
            # elif (line[0].isalpha() and len(line.strip()) > 1):
            #     currentKey = line.strip()
            # new line
            else:
                if (currentKey not in item):
                    item[currentKey] = ""
                item[currentKey] += line + "\n"
            line = fp.readline()
            if (line == newGroupKey):
                # start new timestep
                out.append(item)
                item = {}
    out.append(item)
    # print("Read file. Parsing...")
    # parse items
    for i, item in enumerate(out):
        for key in item:
            # remove upper, which are indices but not header
            headerSplits = [s.strip()
                            for s in key.split(" ") if not s.isupper()]
            if (len(headerSplits) > 1):
                # make headers unique
                # TODO: generalize for potential other headers
                if (key.startswith("BOX BOUNDS")):
                    headerSplits = ["first", "second"]
                out[i][key] = pd.read_csv(
                    StringIO(item[key]), header=None, names=headerSplits, sep=" ", index_col=False)

    doCache(out, file, "-lammp-dump-cache")
    return out
