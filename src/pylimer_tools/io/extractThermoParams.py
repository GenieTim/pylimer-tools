
import csv
import datetime
import hashlib
import os
import pathlib
import pickle
from io import StringIO

import pandas as pd

from pylimer_tools.utils.optimizeDf import optimize, reduce_mem_usage


# Helper functions
def readOneGroup(fp, header, minLineLen=4, additional_lines_skip=0) -> str:
    """
    Read one group of csv lines from the file

    Arguments:
        fp: the file pointer to the file to read from
        header: the header of the CSV (where to start reading at)
        minLineLen: the minimal length of a line to be accepted as data
        additional_lines_skip: number of lines to skip after reading the header
    """
    text = ""
    line = fp.readline()
    separator = ", "
    headerLen = None
    if (isinstance(header, str)):
        minLineLen = max(minLineLen, len(header.split()))
        headerLen = len(header.split())
    else:
        minLineLen = max(minLineLen, min([len(header) for h in header]))

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
    if (headerLen is None):
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


def extractThermoParams(file, header="Temp PotEng TotEng Press Volume c_3", textsToRead=5, minLineLen=5, useCache=True) -> pd.DataFrame:
    """
    Extract the thermodynamic outputs produced for this simulation.

    Note: the header parameter can be an array — make sure to pay attention
    when reading a file with different header sections in them

    Arguments:
        file: the file path to the file to read from
        header: the header of the CSV (where to start reading at)
        textsToRead: the number of times to expect the header
        minLineLen: the minimal length of a line to be accepted as data
        useCache: wheter to use cache or not

    Returns: 
        data (pd.DataFrame): the thermodynamic parameters
    """
    df = None

    cacheFileName = os.path.dirname(
        __file__) + "/cache/" + hashlib.md5(file.encode()).hexdigest() + "-thermo-param-cache.pickle"

    if (os.path.isfile(cacheFileName) and useCache):
        mtimeCache = datetime.datetime.fromtimestamp(
            pathlib.Path(cacheFileName).stat().st_mtime)
        mtimeOrigin = datetime.datetime.fromtimestamp(
            pathlib.Path(file).stat().st_mtime)
        if (mtimeCache > mtimeOrigin):
            toReturn = None
            with open(cacheFileName, 'rb') as cacheFile:
                toReturn = pickle.load(cacheFile)
            if (toReturn is not None):
                print("Read {} rows for file {} from cache".format(
                    len(toReturn), file))
                return toReturn
        else:
            print("Dump cache file is elder than dump. Reloading...")

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

    if (not os.path.exists(os.path.dirname(cacheFileName))):
        os.makedirs(os.path.dirname(cacheFileName))

    with open(cacheFileName, 'wb') as cacheFile:
        pickle.dump(df, cacheFile)

    print("Read {} rows for file {}".format(len(df), file))

    return df
