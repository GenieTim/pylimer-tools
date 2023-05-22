
import base64
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
from pylimer_tools_cpp import splitCSV


def _is_numeric_string(test: str):
    return np.all([c.isnumeric() or c == "." or c == "+" or c == "-" or c == "e" or c == "E" for c in test.strip()])

def detectHeaders(file: str, max_nr_of_lines_to_read: int = 1500, use_cache: bool = True) -> List[str]:
    """
    Read `max_nr_of_lines_to_read` lines from the given file and return all possible header lines.

    Some assumptions are made regarding the columns, e.g., that 75% of them start with a character.

    Arguments:
        - file: The file to search for header lines
        - max_nr_of_lines_to_read: The number of lines to read in search for header lines. 
            Use a negative number to read the whole file.
        - use_cache: Whether to read the result from cache or not. 
            The cache is not read if the file changed meanwhile.
    """
    suffix = str(max_nr_of_lines_to_read)
    cacheContent = loadCache(file, suffix)

    if (cacheContent is not None and use_cache):
        return cacheContent

    lines_read = 0
    previous_line = None
    results = []
    with open(file, 'r') as f:
        for line in f:
            if (previous_line is not None and len(line.strip().split()) == len(previous_line.removeprefix("#").strip().split()) and np.sum([
                w[0].isalpha() for w in previous_line.split()
            ]) > 0.74*len(previous_line.split()) and np.sum([
                _is_numeric_string(w) for w in line.split()
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
            if (lines_read > max_nr_of_lines_to_read and max_nr_of_lines_to_read > 0):
                break

    doCache(results, file, suffix)
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
    if (len(header) == 0):
        raise ValueError("header must have more than zero characters")
    assert(isinstance(header, str) or (
        isinstance(header, list) and len(header) > 0))
    csv_file_to_write = "{}/{}_{}".format(
        tempfile.gettempdir(),
        hashlib.md5(datetime.now().strftime("%m.%d.%Y, %H:%M:%S.%f").encode()).hexdigest(), 'tmp_thermo_file.csv')
    n_lines = 0
    with open(csv_file_to_write, 'w') as output_csv:
        line = fp.readline()
        separator = ", "
        headerLen = None
        if (isinstance(header, str)):
            min_line_len = max(min_line_len, len(header.split()))
        else:
            min_line_len = max(min_line_len, min(
                [len(h.split()) for h in header]))

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
    return csv_file_to_write if n_lines > 0 else ""


def getThermoCacheNameSuffix(header: Union[str, List[str], None] = "Step Temp E_pair E_mol TotEng Press", texts_to_read: float = 50, min_line_len: float = 5) -> str:
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


def extractThermoParams(file, header: Union[str, List[str], None] = "Step Temp E_pair E_mol TotEng Press", texts_to_read: float = 50, min_line_len: float = 5, use_cache: bool = True, lines_to_read_to_detect_header: float = 1e5, lines_to_read_till_header: float = -1) -> pd.DataFrame:
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
            The cache is not read if the file changed meanwhile.
        - lines_to_read_till_header: the number of lines that are acceptable to skip until a header should have been found.
            This is useful for (a) finding the header, and (b) exit early if you are unsure about the header(s)

    Returns:
        - data (pd.DataFrame): the thermodynamic parameters

    """
    df = None

    if (header is None):
        header = detectHeaders(
            file, max_nr_of_lines_to_read=lines_to_read_to_detect_header if lines_to_read_to_detect_header > 0 else 1500)
        if (len(header) == 0):
            raise RuntimeError(
                "Failed to find suitable header. Set a higher value of `lines_to_read_to_detect_header` if you insist that the file '{}' is appropriate.".format(file))

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
                warnings.warn("Could not remove file {}".format(filePath))
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
        tmpCsvFiles = []
        if (tmpCsvFile != ""):
            tmpCsvFiles.append(tmpCsvFile)
        while(textsRead < texts_to_read):
            tmpCsvFile = readOneGroup(
                fp, header, min_line_len=min_line_len, lines_to_read_till_header=lines_to_read_till_header)
            textsRead += 1
            if (tmpCsvFile != ""):
                tmpCsvFiles.append(tmpCsvFile)
            else:
                break
        if (len(tmpCsvFiles) == 1):
            df = csvFileToDf(tmpCsvFiles[0])
        elif (len(tmpCsvFiles) > 0):
            df = pd.concat([df for df in [csvFileToDf(f)
                           for f in tmpCsvFiles] if not df.empty], ignore_index=True)

    if (df is not None):
        # df.columns = df.columns.str.replace(' ', '')
        df.rename(columns=lambda x: x.strip(), inplace=True)
    else:
        df = pd.DataFrame()

    doCache(df, file, suffix)
    # print("Read {} rows for file {}".format(len(df), file))

    return df


def readMultiSectionSeparatedValueFile(file: str, separator: str = None, use_cache: bool = True, comment: str = None, skip_err: bool = False):
    """
    Reads a tsv-like file (could also be e.g. space separated, or csv if you use separator = ",")
    which contains multiple headers throughout the file.

    Particularly useful to read e.g. a file of output from the DPDSimulator.

    Parameters:
        - file: the path to the file to read
        - separator: the separator if not one of the defaults used/recognized by pandas' read_csv
        - use_cache: use to disable reading from cached results
        - comment: a character such as "#" to indicate the separator for where a comment starts
    """
    suffix = (base64.urlsafe_b64encode(comment.encode("utf-8")).decode('utf-8') if comment is not None else "") + \
        "mssv2-" + \
        base64.urlsafe_b64encode(
            separator.encode("utf-8")).decode('utf-8') if separator is not None else "-any"
    cacheContent = loadCache(file, suffix)

    if (cacheContent is not None and use_cache):
        return cacheContent

    print("Splitting CSV...")
    previous_len = -1

    tmp_csv_files = splitCSV(file, separator)
    print("CSV split to {} files... e.g. to {}, {} or {}".format(len(tmp_csv_files), tmp_csv_files[0], tmp_csv_files[1] if len(
        tmp_csv_files) > 1 else "", tmp_csv_files[2] if len(tmp_csv_files) > 2 else ""))

    if (len(tmp_csv_files) == 0):
        return pd.DataFrame()

    # determine the columns we want to have in the end
    all_headers = set()
    detected_dtypes = {}
    erronous_files = []
    for csv_file in tmp_csv_files:
        header_line = ""
        first_line = ""
        got_err = False
        with(open(csv_file, 'r')) as fp:
            try:
                header_line = next(fp)
                first_line = next(fp)
            except StopIteration as e:
                erronous_files.append(csv_file)
                got_err = True
        if (got_err):
            continue
        headers = header_line.strip().split(separator)
        for i, h in enumerate(headers):
            if (h not in all_headers):
                first_line_split = first_line.strip().split(separator)
                if (len(first_line_split) != len(headers)):
                    raise ValueError(
                        "Headers and first line do not match in nr of values", first_line, header_line)
                if (np.all([c.isdigit() or c == "-" for c in first_line_split[i]])):
                    detected_dtypes[h] = np.int64
                elif (np.all([c.isdigit() or c == "-" or c == "." or c == "e" or c == "E" for c in first_line_split[i]])):
                    detected_dtypes[h] = np.float64
                all_headers.add(h)
    all_headers = list(all_headers)
    csv_file_to_write = "{}/{}_{}".format(
        tempfile.gettempdir(),
        hashlib.md5(datetime.now().strftime("%m.%d.%Y, %H:%M:%S.%f").encode()).hexdigest(), 'tmp_mssv2_file.csv')

    print("{} Headers mapped...".format(len(all_headers)))

    # re-join the CSV files in one big file with all the columns
    # put NaN where we do not have a value for a column
    with open(csv_file_to_write, 'w') as outFile:
        outFile.write(separator.join(all_headers) + "\n")
        for csv_file in tmp_csv_files:
            if (csv_file in erronous_files):
                print("File {} skipped".format(csv_file))
                try:
                    os.remove(csv_file)
                except OSError as e:
                    warnings.warn("Could not remove file {}".format(csv_file))
                    pass
                continue
            with(open(csv_file, 'r')) as fp:
                header_line = next(fp)
                split_header = header_line.strip().split(separator)
                map_to_col = []
                n_found = 0
                for i, col in enumerate(all_headers):
                    if (col in split_header):
                        map_to_col.append(split_header.index(col))
                        n_found += 1
                    else:
                        map_to_col.append(-1)
                assert(n_found == len(split_header))
                for line in fp:
                    if (line == header_line or line.startswith("Step")):
                        continue
                    split_line = line.strip().split(separator)
                    str_to_write = separator.join(
                        [split_line[i] if i != -1 else "NaN" for i in map_to_col])
                    outFile.write(str_to_write + "\n")
            try:
                os.remove(csv_file)
            except OSError as e:
                warnings.warn("Could not remove file {}".format(csv_file))
                pass
            print("File {} handled".format(csv_file))
    # read the csv files again
    print("Reading final csv file {}".format(csv_file_to_write))
    try:
        df = pd.read_csv(
            csv_file_to_write, sep=separator, comment=comment, dtype=detected_dtypes, na_values=["NaN"])
    except pd.errors.EmptyDataError as e:
        warnings.warn("Data file {} turned out to be empty".format(file))
        return pd.DataFrame()
    doCache(df, file, suffix)
    try:
        os.remove(csv_file_to_write)
    except OSError as e:
        warnings.warn("Could not remove file {}".format(csv_file_to_write))
        pass
    # doCache(reduce_mem_usage(df), file, suffix)
    # print("Read {} rows for file {}".format(len(df), file))

    return df
