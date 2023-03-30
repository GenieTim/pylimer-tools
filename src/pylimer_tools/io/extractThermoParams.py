
import csv
import hashlib
import os
import pathlib
import pickle
import tempfile
import warnings
from datetime import datetime
from io import StringIO
from typing import Iterable, List, Union

import numpy as np
import pandas as pd

from pylimer_tools.utils.cacheUtility import doCache, loadCache
from pylimer_tools.utils.optimizeDf import optimize, reduce_mem_usage


def detectHeaders(file: str, max_nr_of_lines_to_read: int = 1500) -> List[str]:
    """
    Read `max_nr_of_lines_to_read` lines from the given file and return all possible header lines.

    Some assumptions are made regarding the columns, e.g., that 75% of them start with a character.
    """
    lines_read = 0
    previous_line = None
    results = []
    with open(file, 'r') as f:
        for line in f:
            if (previous_line is not None and len(line.strip().split()) == len(previous_line.removeprefix("#").strip().split()) and np.sum([
                w[0].isalpha() for w in previous_line.split()
            ]) > 0.74*len(previous_line.split()) and np.sum([
                np.all([c.isnumeric() or c == "." or c == "+" or c == "-" or c == "e" for c in w.strip()]) for w in line.split()
            ]) > 0.5*len(line.split()) and "..." not in previous_line and len(previous_line.split()) > 2 and not np.any([
                previous_line.startswith(val) for val in [
                    "Memory usage per processor",
                    "Setting up Verlet run",
                    "Dangerous builds",
                    "<",
                    "Started at",
                    "Terminated at",
                    "Results reported at",
                    "WARNING"
                ]
            ])):
                results.append(previous_line.rstrip())
            previous_line = line
            lines_read += 1
            if (lines_read > max_nr_of_lines_to_read):
                break
    return results


def readOneGroup(fp, header, min_line_len=4, additional_lines_skip=0, lines_to_read_till_header=1e3) -> str:
    """
    Read one group of csv lines from the file

    Arguments:
        - fp: the file pointer to the file to read from
        - header: the header of the CSV (where to start reading at)
        - min_line_len: the minimal length of a line to be accepted as data
        - additional_lines_skip: number of lines to skip after reading the header


    Returns:
       The filename of a temporary CSV file
    """
    assert(isinstance(header, str) or (
        isinstance(header, list) and len(header) > 0))
    csvFileToWrite = "{}/{}_{}".format(
        tempfile.gettempdir(),
        hashlib.md5(datetime.now().strftime("%m.%d.%Y, %H:%M:%S.%f").encode()).hexdigest(), 'tmp_thermo_file.csv')
    n_lines = 0
    with open(csvFileToWrite, 'w') as output_csv:
        line = fp.readline()
        separator = ", "
        headerLen = None
        if (isinstance(header, str)):
            min_line_len = max(min_line_len, len(header.split()))
        else:
            min_line_len = max(min_line_len, min([len(h.split()) for h in header]))

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
        linesSkipped = 0
        while skipLineFun(line, header) and line.endswith("\n"):
            line = fp.readline()
            linesSkipped += 1
            if (linesSkipped > lines_to_read_till_header and lines_to_read_till_header > 0):
                raise RuntimeError(
                    "Skipped {} lines, not encountered any header yet.".format(linesSkipped))
        # found header. Take next few lines:
        headerLen = len(line.split())
        if (not line):
            return ""
        else:
            output_csv.write((separator.join(line.split())).strip() + "\n")

        n_lines = 0
        while line and n_lines < additional_lines_skip:
            # skip ${additional_lines_skip} further
            line = fp.readline()
            # text += (', '.join(line.split())).strip() + "\n"
            n_lines += 1
        while line and not line.startswith("Loop time of"):
            line = fp.readline()
            if (len(line) < min_line_len or (len(line.split()) != headerLen) or (len(line) > 0 and (
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
            output_csv.write((separator.join(line.split())).strip() + "\n")
            n_lines += 1
    return csvFileToWrite if n_lines > 0 else ""


def getThermoCacheNameSuffix(header: Union[str, List[str], None] = "Step Temp E_pair E_mol TotEng Press", texts_to_read:float=50, min_line_len:float=5) -> str:
    """
    Compose a cache file suffix in such a way, that it distinguishes different thermo reader parameters

    Arguments:
        - header: the header of the CSV (where to start reading at)
        - texts_to_read: the number of times to expect the header
        - min_line_len: the minimal length of a line to be accepted as data
    """
    if (isinstance(header, Iterable)):
        header = "{}{}".format("".join("".join(header).split()), len(header))

    # need to has header, as we could get a filename too long error otherwise. Addmittedly, still possible for certain inputs
    return "{}{}{}-thermo-param-cache.pickle".format(hashlib.md5(header.encode()).hexdigest() if header is not None else "", texts_to_read, min_line_len)


def extractThermoParams(file, header: Union[str, List[str], None] = "Step Temp E_pair E_mol TotEng Press", texts_to_read:float=50, min_line_len:float=5, use_cache:bool=True, lines_to_read_to_detect_header:float=1e5, lines_to_read_till_header:float=-1) -> pd.DataFrame:
    """
    Extract the thermodynamic outputs produced for this simulation.

    Note: the header parameter can be an array — make sure to pay attention
    when reading a file with different header sections in them

    Arguments:
        - file: the file path to the file to read from
        - header: the header of the CSV (where to start reading at). 
            Can be a string, a list of strings, or None if you want to try the detection.
        - texts_to_read: the number of times to expect the header
        - min_line_len: the minimal length of a line to be accepted as data
        - use_cache: wheter to use cache or not (though it will be written anyway)
        - lines_to_read_till_header: the number of lines that are acceptable to skip until a header should have been found.
            This is useful for (a) finding the header, and (b) exit early if you are unsure about the header(s)

    Returns:
        - data (pd.DataFrame): the thermodynamic parameters

    """
    df = None

    if (header is None):
        header = detectHeaders(
            file, max_nr_of_lines_to_read=lines_to_read_to_detect_header if lines_to_read_to_detect_header > 0 else 1500)

    suffix = getThermoCacheNameSuffix(
        header, texts_to_read, min_line_len)
    cacheContent = loadCache(file, suffix)

    if (cacheContent is not None and use_cache):
        return cacheContent

    def csvFileToDf(filePath) -> pd.DataFrame:
        try:
            tmpDf = pd.read_csv(filePath, low_memory=False,
                                on_bad_lines='skip', quoting=csv.QUOTE_NONE)
            try:
                os.remove(filePath)
            except Exception as e:
                pass
            return tmpDf
        except Exception as e:
            warnings.warn("Error reading temporary CSV thermo file '{}': {}".format(
                filePath, e), source=e)
            return pd.DataFrame()

    with open(file, 'r') as fp:
        tmpCsvFile = readOneGroup(
            fp, header, min_line_len=min_line_len, lines_to_read_till_header=lines_to_read_till_header)
        textsRead = 1
        if (tmpCsvFile == ""):
            df = pd.DataFrame()
        else:
            df = csvFileToDf(tmpCsvFile)
        while(textsRead < texts_to_read):
            tmpCsvFile = readOneGroup(
                fp, header, min_line_len=min_line_len, lines_to_read_till_header=lines_to_read_till_header)
            textsRead += 1
            if (tmpCsvFile != ""):
                newDf = csvFileToDf(tmpCsvFile)
                if (not newDf.empty):
                    df = pd.concat([df, newDf], ignore_index=True)
            else:
                break

    if (df is not None):
        # df.columns = df.columns.str.replace(' ', '')
        df.rename(columns=lambda x: x.strip(), inplace=True)
    else:
        df = pd.DataFrame()

    doCache(df, file, suffix)
    # print("Read {} rows for file {}".format(len(df), file))

    return df
