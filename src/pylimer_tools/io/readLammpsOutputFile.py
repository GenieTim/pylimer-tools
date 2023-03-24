import os
import re

import pandas as pd

from pylimer_tools.io.extractThermoParams import extractThermoParams
from pylimer_tools_cpp import Universe, UniverseSequence


def readLogFile(filepath) -> pd.DataFrame:
    return extractThermoParams(filepath, header=None, textsToRead=1e5)


def readDumpFile(filepath, structureFile) -> UniverseSequence:
    uS = UniverseSequence()
    uS.initializeFromDataSequence(filepath, structureFile)
    return uS


def readDataFile(structureFile: str) -> Universe:
    uS = UniverseSequence()
    uS.initializeFromDataSequence([structureFile])
    universe = uS.atIndex(0)
    del uS
    return universe


def readAveragesFile(filepath) -> pd.DataFrame:
    """
    Read a file written by a `fix ave/time` command.

    Uses pandas' read_csv after detecting the columns.

    Important assumtion: the first 2 or 3 lines in the file are:
    - comment,
    - then one header indicating the columns,
    - and then either data or potentially a second header, if it is a sectioned file (e.g., from a `fix ave/time ... vector`)
    """
    assert(os.path.isfile(filepath))
    header_line = None
    with open(filepath, 'r') as f:
        line0 = f.readline()
        line1 = f.readline()
        line2 = f.readline()

        if (line2.startswith("#")):
            return readSectionedAveragesFile(filepath)

        header_line = line1
    header_line = header_line.removesuffix("#").strip()

    data = pd.read_csv(filepath, comment="#",
                       names=header_line.split(), sep=" ")

    return data


def readSectionedAveragesFile(filepath) -> pd.DataFrame:
    """
    Read a file written by a `fix ave/time` command.

    Use the section delimeter columns together with pandas' groupby() 
    to restore the original sections.
    """
    assert(os.path.isfile(filepath))
    data = {}
    with open(filepath, 'r') as f:
        f.readline()  # discard line 0
        line1 = f.readline()
        line2 = f.readline()

        if (not line2.startswith("#")):
            raise ValueError(
                "This file was not detected to be a proper sectioned averages file.")
            # return readSectionedAveragesFile(filepath)

        header_line1 = line1.removesuffix("#").strip()
        header_line2 = line2.removesuffix("#").strip()

        header_line1_split = header_line1.split()
        header_line2_split = header_line2.split()

        if (len(header_line1_split) == 1 and len(header_line2_split)):
            raise ValueError(
                "Cannot read this file, as we cannot distinguish between section header and main data")

        currentData = []
        currentKey = None
        for line in f:
            if (currentKey == None):
                assert(len(line.split()) == len(header_line1.split()))
                currentKey = line
                continue
            splitLine = line.split()
            if (len(splitLine) == len(header_line1_split)):
                data[currentKey] = currentData
                currentData = []
                currentKey = line
            else:
                assert(len(splitLine) == len(header_line2_split))
                currentData.append(splitLine)
        data[currentKey] = currentData

    # convert all the data to a dataframe
    df = pd.DataFrame()
    for key in data.keys():
        split_key = key.split()
        local_dataframe = pd.DataFrame(data[key], columns=header_line2_split)
        for i, col in enumerate(header_line2_split):
            local_dataframe[col] = split_key[i]
        df = pd.concat([df, local_dataframe], as_index=False)

    # convert all columns of DataFrame
    df = df.apply(pd.to_numeric, errors='ignore')

    return df


def readHistogramFile(filepath) -> pd.DataFrame:
    """
    Read a file written by `fix ave/hist` or similar.

    See also:
        - `func:readSectionedAveragesFile`
    """
    return readSectionedAveragesFile(filepath)


def readCorrelationFile(filepath, group_key="Timestep") -> pd.DataFrame:
    """
    Read a file written by a `fix ave/correlate{/long}` command.

    Parameters:
        - filepath: the path to the file to read
        - group_key: the key that denotes a new section

    Returns:
        - correlatedData: a DataFrame containing all the data of the file.
            Use the group_key with the DataFrame's groupby() to restore the original sections.    
    """
    assert(os.path.isfile(filepath))
    data = {}
    header_line = None
    with open(filepath, 'r') as f:
        currentData = []
        currentKey = None
        normalLineLen = None
        line0 = f.readline()
        header_line = f.readline()
        lines_interpreted = 0
        for line in f:
            if ((line.startswith("#") or len(line.strip()) == 0) and group_key not in line):
                if (lines_interpreted == 0):
                    header_line = line
                continue
            split = line.split()
            if (len(split) == 2 or group_key in line):
                if (currentKey is not None):
                    data[currentKey] = currentData
                    currentData = []
                # new key
                currentKey = line
            elif (len(split) == normalLineLen or normalLineLen is None):
                normalLineLen = len(split)
                currentData.append([int(val) if (
                    "." not in val and "e" not in val) else float(val) for val in split])
            else:
                raise ValueError(
                    "Did not expect {} splited values on line with content {} in rdf file {}".format(len(split), line, filepath))
            lines_interpreted += 1
        data[currentKey] = currentData

    correlatedDataAssembled = []
    for key in data.keys():
        assert (group_key in key)
        compiledRegex = re.compile(r"{}: ([\d]+)".format(group_key))
        results = compiledRegex.search(key)
        assert(results is not None)
        timestep = int(results.group(1))
        for row in data[key]:
            row.append(timestep)
            correlatedDataAssembled.append(row)

    cols = header_line.split()
    cols.append(group_key)
    correlatedData = pd.DataFrame(
        correlatedDataAssembled, columns=cols)
    return correlatedData
