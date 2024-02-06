import os
import re
import warnings
from typing import Iterable, List, Union

import numpy as np
import pandas as pd

from pylimer_tools.io.extractThermoParams import extractThermoParams
from pylimer_tools.utils.cacheUtility import doCache, loadCache
from pylimer_tools_cpp.pylimer_tools_cpp import (AtomStyle, Universe,
                                                 UniverseSequence)


def readLogFile(filepath, lines_to_read_to_detect_header=500000) -> pd.DataFrame:
    return extractThermoParams(filepath, header=None, texts_to_read=500000, lines_to_read_to_detect_header=lines_to_read_to_detect_header)


def readDumpFile(dataFile, dumpFile, atom_style: Union[List[AtomStyle], None] = None) -> UniverseSequence:
    """
    Read a file with LAMMPS' dump of snapshots of structures into a Universe.
    """
    uS = UniverseSequence()
    if (atom_style is not None):
        uS.setDataFileAtomStyle(atom_style)
    uS.initializeFromDumpFile(dataFile, dumpFile)
    return uS


def readDataFile(structureFile: str, atom_style: Union[List[AtomStyle], None] = None) -> Universe:
    """
    Read a file with LAMMPS' data type of structure into a Universe.
    """
    uS = UniverseSequence()
    if (atom_style is not None):
        uS.setDataFileAtomStyle(atom_style)
    uS.initializeFromDataSequence([structureFile])
    universe = uS.atIndex(0)
    del uS
    return universe


def readAveragesFile(filepath, use_cache: bool = True, sep=" ") -> pd.DataFrame:
    """
    Read a file written by a `fix ave/time` command.

    Uses pandas' read_csv after detecting the columns.

    Important assumtion: the first 2 or 3 lines in the file are:
    - comment,
    - then one header indicating the columns,
    - and then either data or potentially a second header, if it is a sectioned file (e.g., from a `fix ave/time ... vector`)
    """
    assert (os.path.isfile(filepath))
    header_line = None
    with open(filepath, 'r') as f:
        line0 = f.readline()
        line1 = f.readline()
        line2 = f.readline()

        if (line2.startswith("#")):
            return readSectionedAveragesFile(filepath, use_cache=use_cache)

        header_line = line1 if line1.startswith("#") else line0
    header_line = header_line.removeprefix("#").strip()

    try:
        data = pd.read_csv(filepath, comment="#",
                           names=header_line.split(), sep=sep)
    except pd.errors.EmptyDataError as e:
        return pd.DataFrame()

    return data


def readSectionedAveragesFile(filepath, use_cache: bool = True) -> pd.DataFrame:
    """
    Read a file written by a `fix ave/time` command.

    Use the section delimeter columns together with pandas' groupby()
    to restore the original sections.
    """
    assert (os.path.isfile(filepath))

    cache_suffix = "sectionedavg-cache.pickle"
    cacheContent = loadCache(filepath, cache_suffix)

    if (cacheContent is not None and use_cache):
        return cacheContent

    data = {}
    with open(filepath, 'r') as f:
        f.readline()  # discard line 0
        line1 = f.readline()
        line2 = f.readline()

        if (not line2.startswith("#")):
            raise ValueError(
                "The file '{}' was not detected to be a proper sectioned averages file.".format(filepath))
            # return readSectionedAveragesFile(filepath)

        header_line1 = line1.removeprefix("#").strip()
        header_line2 = line2.removeprefix("#").strip()

        header_line1_split = header_line1.split()
        header_line2_split = header_line2.split()

        if (len(header_line1_split) == len(header_line2_split)):
            raise ValueError(
                "Cannot read this file, as we cannot distinguish between section header and main data")

        currentData = []
        currentKey = None
        for line in f:
            splitLine = line.split()
            if (currentKey == None):
                assert (len(splitLine) == len(header_line1.split()))
                currentKey = line
                continue
            if (len(splitLine) == len(header_line1_split)):
                data[currentKey] = currentData
                currentData = []
                currentKey = line
            else:
                assert (len(splitLine) == len(header_line2_split))
                currentData.append(splitLine)
        data[currentKey] = currentData

    # convert all the data to a dataframe
    dfs_to_concat = []

    if (header_line1_split is None):
        raise ValueError("Did not find a useable header line.")

    for key in data.keys():
        split_key = key.split()
        local_dataframe = pd.DataFrame(data[key], columns=header_line2_split)
        for i, col in enumerate(header_line1_split):
            local_dataframe[col] = split_key[i]
        dfs_to_concat.append(local_dataframe)

    df = pd.concat(dfs_to_concat, ignore_index=True)

    # convert all columns of DataFrame
    df = df.apply(pd.to_numeric, errors='ignore')
    doCache(df, filepath, cache_suffix)

    return df


def readHistogramFile(filepath, use_cache: bool = True) -> pd.DataFrame:
    """
    Read a file written by `fix ave/hist` or similar.

    See also:
        - `func:readSectionedAveragesFile`
    """
    return readSectionedAveragesFile(filepath, use_cache)


def readCorrelationFile(filepath, group_key="Timestep", use_cache: bool = True) -> pd.DataFrame:
    """
    Read a file written by a `fix ave/correlate{/long}` command.

    Parameters:
        - filepath: the path to the file to read
        - group_key: the key that denotes a new section

    Returns:
        - correlatedData: a DataFrame containing all the data of the file.
            Use the group_key with the DataFrame's groupby() to restore the original sections.
    """
    assert (os.path.isfile(filepath))

    cache_suffix = "{}-correlation-cache.pickle".format(
        group_key if isinstance(group_key, str) else "g")
    cacheContent = loadCache(filepath, cache_suffix)

    if (cacheContent is not None and use_cache):
        return cacheContent

    data = {}
    header_line = None
    with open(filepath, 'r') as f:
        currentData = []
        currentKey = None
        header_line = f.readline()
        if (header_line.startswith("#")):
            # in LAMMPS files, there is a title line that does not exist in our DPD output,
            # -> this line is needed for LAMMPS
            header_line = f.readline()
        cols = header_line.removeprefix("#").strip().split()
        normalLineLen = len(cols)
        lines_interpreted = 0

        def isGroupKey(line):
            # if (isinstance(group_key, list)):
            #     return np.any([x in line for x in group_key])
            # else:
            return group_key in line

        for line in f:
            if ((line.startswith("#") or len(line.strip()) == 0) and not isGroupKey(line)):
                if (lines_interpreted == 0):
                    header_line = line
                continue
            if (line == header_line):
                continue
            split = line.removeprefix("#").strip().split()
            if (len(split) == 2 or isGroupKey(line)):
                if (currentKey is not None and len(currentData) > 0):
                    data[currentKey] = currentData
                    currentData = []
                # new key
                currentKey = line
            elif (len(split) == normalLineLen or normalLineLen is None):
                # normalLineLen = len(split)
                currentData.append(split)
            else:
                raise ValueError(
                    "Did not expect {} splited values on line with content {} in correlation file {}".format(len(split), line, filepath))
            lines_interpreted += 1
        if (currentKey is not None and len(currentData) > 0):
            data[currentKey] = currentData

    cols.append(group_key)
    correlatedDataAssembled = []
    for key in data.keys():
        assert (group_key in str(key))
        compiledRegex = re.compile(r"{}:? ([\d]+)".format(group_key))
        results = compiledRegex.search(key)
        if (results is None):
            warnings.warn("Did not find {} with number in {} when reading {}".format(
                group_key, key, filepath))
        assert (results is not None)
        timestep = int(results.group(1))
        for row in data[key]:
            row.append(timestep)
            assert (len(row) == len(cols))
            correlatedDataAssembled.append(row)

    correlatedData = pd.DataFrame(
        correlatedDataAssembled, columns=cols)
    # convert all columns of DataFrame
    correlatedData = correlatedData.apply(pd.to_numeric, errors='ignore')
    doCache(correlatedData, filepath, cache_suffix)

    return correlatedData
