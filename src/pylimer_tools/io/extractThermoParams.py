
import csv
import datetime
import hashlib
import os
import pathlib
import pickle
from io import StringIO
from typing import Iterable

import pandas as pd
from pylimer_tools.utils.cacheUtility import doCache, loadCache
from pylimer_tools.utils.optimizeDf import optimize, reduce_mem_usage


# Helper functions
def readOneGroup(fp, header, minLineLen=4, additional_lines_skip=0) -> str:
    """
    Read one group of csv lines from the file

    Arguments:
        - fp: the file pointer to the file to read from
        - header: the header of the CSV (where to start reading at)
        - minLineLen: the minimal length of a line to be accepted as data
        - additional_lines_skip: number of lines to skip after reading the header


    Returns:
      A long CSV string
    """
    text = ""
    line = fp.readline()
    separator = ", "
    headerLen = None
    if (isinstance(header, str)):
        minLineLen = max(minLineLen, len(header.split()))
    else:
        minLineLen = max(minLineLen, min([len(h.split()) for h in header]))

    def checkSkipLine(line, header):
        return line and not line.startswith(header)

    def checkSkipLineHeaderList(line, header):
        if (not line):
            return False
        for headerL in header:
            if (line.startswith(headerL)):
                return False
        return True

    skipLineFun = checkSkipLineHeaderList if isinstance(
        header, list) else checkSkipLine
    # skip lines up until header (or file ending)
    while skipLineFun(line, header):
        line = fp.readline()
    # found header. Take next few lines:
    headerLen = len(line.split())
    if (not line):
        return ""
    else:
        text = (separator.join(line.split())).strip() + "\n"

    n_lines = 0
    while line and n_lines < additional_lines_skip:
        # skip ${additional_lines_skip} further
        line = fp.readline()
        # text += (', '.join(line.split())).strip() + "\n"
        n_lines += 1
    while line and not line.startswith("Loop time of"):
        line = fp.readline()
        if (len(line) < minLineLen or (len(line.split()) != headerLen) or (len(line) > 0 and (
                line.startswith("WARNING") or
                line[0].isalpha() or
                (line[0] == "-" and line[1] == "-") or
                (line[2].isalpha() or line[3].isalpha()) or
                (line[0] == "[") or
                ("src" in line) or
                ("fene" in line or ")" in line)  # from ":90)"
        ))):
            # skip line due to error, warning or similar
            continue
        text += (separator.join(line.split())).strip() + "\n"
        n_lines += 1
    return text


def getThermoCacheNameSuffix(header="Step Temp E_pair E_mol TotEng Press", textsToRead=5, minLineLen=5) -> str:
    """
    Compose a cache file suffix in such a way, that it distinguishes different thermo reader parameters

    Arguments:
        - header: the header of the CSV (where to start reading at)
        - textsToRead: the number of times to expect the header
        - minLineLen: the minimal length of a line to be accepted as data
    """
    if (isinstance(header, Iterable)):
        header = "{}{}".format("".join(header), len(header))

    return "{}{}{}-thermo-param-cache.pickle".format(header, textsToRead, minLineLen)


def extractThermoParams(file, header="Step Temp E_pair E_mol TotEng Press", textsToRead=5, minLineLen=5, useCache=True) -> pd.DataFrame:
    """
    Extract the thermodynamic outputs produced for this simulation.

    Note: the header parameter can be an array — make sure to pay attention
    when reading a file with different header sections in them

    Arguments:
        - file: the file path to the file to read from
        - header: the header of the CSV (where to start reading at)
        - textsToRead: the number of times to expect the header
        - minLineLen: the minimal length of a line to be accepted as data
        - useCache: wheter to use cache or not (though it will be written anyway)

    Returns:
        - data (pd.DataFrame): the thermodynamic parameters

    """
    df = None

    suffix = getThermoCacheNameSuffix(
        header, textsToRead, minLineLen)
    cacheContent = loadCache(file, suffix)

    if (cacheContent is not None and useCache):
        return cacheContent

    def stringToDf(text) -> pd.DataFrame:
        try:
            return pd.read_csv(StringIO(text), low_memory=False,
                               error_bad_lines=False, quoting=csv.QUOTE_NONE)
        except Exception as e:
            return pd.DataFrame()

    with open(file, 'r') as fp:
        text = readOneGroup(fp, header, minLineLen=minLineLen)
        textsRead = 1
        df = stringToDf(text)
        while(textsRead < textsToRead):
            text = readOneGroup(fp, header, minLineLen=minLineLen)
            textsRead += 1
            if (text != ""):
                newDf = stringToDf(text)
                if (not newDf.empty):
                    df = df.append(newDf)

    if (df is not None):
        # df.columns = df.columns.str.replace(' ', '')
        df.rename(columns=lambda x: x.strip(), inplace=True)
    else:
        df = pd.DataFrame()

    doCache(df, file, suffix)
    # print("Read {} rows for file {}".format(len(df), file))

    return df
