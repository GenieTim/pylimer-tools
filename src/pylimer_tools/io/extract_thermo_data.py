import base64
import csv
import hashlib
import os
import re
import tempfile
import warnings
from datetime import datetime
from typing import Iterable, List, Union

import numpy as np
import pandas as pd

from pylimer_tools.utils.cache_utility import do_cache, load_cache
from pylimer_tools_cpp import split_csv


def _is_numeric_string(test: str) -> bool:
    """
    Check if a string represents a numeric value.

    :param test: String to check
    :type test: str
    :return: True if the string represents a numeric value, False otherwise
    :rtype: bool
    """
    return bool(
        all(
            [
                c.isnumeric()
                or c == "."
                or c == "+"
                or c == "-"
                or c == "e"
                or c == "E"
                for c in test.strip()
            ]
        )
    )


def detect_headers(
    file: str, max_nr_of_lines_to_read: int = 1500, use_cache: bool = True
) -> List[str]:
    """
    Read `max_nr_of_lines_to_read` lines from the given file and return all possible header lines.

    Some assumptions are made regarding the columns, e.g., that 75% of them start with a character.

    :param file: The file to search for header lines
    :type file: str
    :param max_nr_of_lines_to_read: The number of lines to read in search for header lines.
                                   Use a negative number to read the whole file.
    :type max_nr_of_lines_to_read: int
    :param use_cache: Whether to read the result from cache or not.
                     The cache is not read if the file changed meanwhile.
    :type use_cache: bool
    :return: List of detected header lines
    :rtype: List[str]
    """
    suffix = str(max_nr_of_lines_to_read)
    cache_content = load_cache(file, suffix)

    if cache_content is not None and use_cache:
        return cache_content

    lines_read = 0
    previous_line = None
    results = []
    with open(file, "r") as f:
        for line in f:
            if (
                previous_line is not None
                and len(line.strip().split())
                == len(previous_line.removeprefix("#").strip().split())
                and sum([w[0].isalpha() for w in previous_line.split()])
                > 0.74 * len(previous_line.split())
                and sum([_is_numeric_string(w) for w in line.split()])
                > 0.5 * len(line.split())
                and "..." not in previous_line
                and len(previous_line.split()) > 2
                and not any(
                    [
                        previous_line.startswith(val)
                        for val in [
                            "Memory usage per processor",
                            "Setting up Verlet run",
                            "Dangerous builds",
                            "<",
                            "Started at",
                            "Terminated at",
                            "Results reported at",
                            "WARNING",
                        ]
                    ]
                )
            ):
                results.append(previous_line.rstrip())
            previous_line = line
            lines_read += 1
            if lines_read > max_nr_of_lines_to_read and max_nr_of_lines_to_read > 0:
                break

    do_cache(results, file, suffix)
    return results


def read_one_group(
    fp, header, min_line_len=4, additional_lines_skip=0, lines_to_read_till_header=1e3
) -> str:
    """
    Read one group of csv lines from the file.

    :param fp: The file pointer to the file to read from
    :type fp: file object
    :param header: The header of the CSV (where to start reading at)
    :type header: str or list
    :param min_line_len: The minimal length of a line to be accepted as data
    :type min_line_len: int
    :param additional_lines_skip: Number of lines to skip after reading the header
    :type additional_lines_skip: int
    :param lines_to_read_till_header: Maximum number of lines to read until finding the header
    :type lines_to_read_till_header: float
    :return: The filename of a temporary CSV file, or empty string if no data was read
    :rtype: str
    """
    if len(header) == 0:
        raise ValueError("header must have more than zero characters")
    assert isinstance(
        header,
        str) or (
        isinstance(
            header,
            list) and len(header) > 0)
    csv_file_to_write = "{}/{}_{}".format(
        tempfile.gettempdir(),
        hashlib.md5(
            datetime.now().strftime("%m.%d.%Y, %H:%M:%S.%f").encode()
        ).hexdigest(),
        "tmp_thermo_file.csv",
    )
    n_lines = 0
    with open(csv_file_to_write, "w") as output_csv:
        line = fp.readline()
        separator = ", "
        header_len = None
        if isinstance(header, str):
            min_line_len = max(min_line_len, len(header.split()))
        else:
            min_line_len = max(min_line_len,
                               min([len(h.split()) for h in header]))

        def check_skip_line(line, header):
            return line and not line.startswith(header)

        def check_skip_line_header_list(line, header):
            if not line:
                return False
            for header_line in header:
                if line.startswith(header_line):
                    return False
            return True

        skip_line_fun = (
            check_skip_line_header_list if isinstance(
                header, list) else check_skip_line
        )
        # skip lines up until header (or file ending)
        n_lines_skipped = 0
        while skip_line_fun(line, header) and line.endswith("\n"):
            line = fp.readline()
            n_lines_skipped += 1
            if (
                n_lines_skipped > lines_to_read_till_header
                and lines_to_read_till_header > 0
            ):
                raise RuntimeError(
                    "Skipped {} lines, not encountered any header yet.".format(
                        n_lines_skipped
                    )
                )
        # found header. Take next few lines:
        header_len = len(line.split())
        if not line:
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
            if (
                len(line) < min_line_len
                or (len(line.split()) != header_len)
                or (
                    len(line) > 0
                    and (
                        line.startswith("WARNING")
                        or line[0].isalpha()
                        or (line[0] == "-" and line[1] == "-")
                        or (line[2].isalpha() or line[3].isalpha())
                        or (line[0] == "[")
                        or ("src" in line)
                        or ("fene" in line or ")" in line)  # from ":90)"
                    )
                )
            ):
                # skip line due to error, warning or similar
                continue
            output_csv.write((separator.join(line.split())).strip() + "\n")
            n_lines += 1
    return csv_file_to_write if n_lines > 0 else ""


def get_thermo_cache_name_suffix(
    header: Union[str, List[str],
                  None] = "Step Temp E_pair E_mol TotEng Press",
    texts_to_read: float = 50,
    min_line_len: float = 5,
) -> str:
    """
    Compose a cache file suffix in such a way, that it distinguishes different thermo reader parameters.

    :param header: The header of the CSV (where to start reading at)
    :type header: Union[str, List[str], None]
    :param texts_to_read: The number of times to expect the header
    :type texts_to_read: float
    :param min_line_len: The minimal length of a line to be accepted as data
    :type min_line_len: float
    :return: A string to be used as cache file suffix
    :rtype: str
    """
    if isinstance(header, Iterable):
        header = "{}{}".format("".join("".join(header).split()), len(header))

    # need to has header, as we could get a filename too long error otherwise.
    # Admittedly, still possible for certain inputs
    return "{}{}{}-thermo-param-cache.pickle".format(
        hashlib.md5(header.encode()).hexdigest() if header is not None else "",
        texts_to_read,
        min_line_len,
    )


def extract_thermo_params(
    file,
    header: Union[str, List[str],
                  None] = "Step Temp E_pair E_mol TotEng Press",
    texts_to_read: int = 50,
    min_line_len: int = 5,
    use_cache: bool = True,
    lines_to_read_to_detect_header: int = int(1e5),
    lines_to_read_till_header: float = -1,
) -> pd.DataFrame:
    """
    Extract the thermodynamic outputs produced for this simulation,
    i.e., in LAMMPS, by the `thermo` command.

    In particular, this function can handle log files,
    handle sections with different columns,
    and handles skipping over warnings as well as broken lines.

    Note: The header parameter can be an array — make sure to pay attention
    when reading a file with different header sections in them.

    :param file: The file path to the file to read from
    :type file: str
    :param header: The header of the CSV (where to start reading at).
                  Can be a string, a list of strings, or None if you want to try the detection.
    :type header: Union[str, List[str], None]
    :param texts_to_read: The number of times to expect the header
    :type texts_to_read: int
    :param min_line_len: The minimal length of a line to be accepted as data
    :type min_line_len: int
    :param use_cache: Whether to use cache or not (though it will be written anyway).
                     The cache is not read if the file changed meanwhile.
    :type use_cache: bool
    :param lines_to_read_to_detect_header: The number of lines to read when trying to detect headers
    :type lines_to_read_to_detect_header: int
    :param lines_to_read_till_header: The number of lines that are acceptable to skip
                                     until a header should have been found.
                                     This is useful for (a) finding the header, and
                                     (b) exit early if you are unsure about the header(s)
    :type lines_to_read_till_header: float
    :return: The thermodynamic parameters
    :rtype: pd.DataFrame
    """
    df = None

    if header is None:
        header = detect_headers(
            file,
            max_nr_of_lines_to_read=(
                lines_to_read_to_detect_header
                if lines_to_read_to_detect_header > 0
                else 1500
            ),
        )
        if len(header) == 0:
            raise RuntimeError(
                "Failed to find suitable header. "
                + "Set a higher value of `lines_to_read_to_detect_header` if the file '{}' is appropriate.".format(
                    file
                )
            )

    suffix = get_thermo_cache_name_suffix(header, texts_to_read, min_line_len)
    cache_content = load_cache(file, suffix)

    if cache_content is not None and use_cache:
        return cache_content

    def csv_file_to_df(filepath) -> pd.DataFrame:
        try:
            tmp_df = pd.read_csv(
                filepath, low_memory=False, on_bad_lines="skip", quoting=csv.QUOTE_NONE
            )
            try:
                os.remove(filepath)
            except Exception as e:
                warnings.warn(
                    "Could not remove file {}: {}".format(
                        filepath, e))
                pass
            return tmp_df
        except Exception as e:
            warnings.warn(
                "Error reading temporary CSV thermo file '{}': {}".format(
                    filepath, e),
                source=e,
            )
            return pd.DataFrame()

    with open(file, "r") as fp:
        tmp_csv_file = read_one_group(
            fp,
            header,
            min_line_len=min_line_len,
            lines_to_read_till_header=lines_to_read_till_header,
        )
        n_texts_read = 1
        tmp_csv_files = []
        if tmp_csv_file != "":
            tmp_csv_files.append(tmp_csv_file)
        while n_texts_read < texts_to_read:
            tmp_csv_file = read_one_group(
                fp,
                header,
                min_line_len=min_line_len,
                lines_to_read_till_header=lines_to_read_till_header,
            )
            n_texts_read += 1
            if tmp_csv_file != "":
                tmp_csv_files.append(tmp_csv_file)
            else:
                break
        if len(tmp_csv_files) == 1:
            df = csv_file_to_df(tmp_csv_files[0])
        elif len(tmp_csv_files) > 0:
            df = pd.concat(
                [
                    df
                    for df in [csv_file_to_df(f) for f in tmp_csv_files]
                    if not df.empty
                ],
                ignore_index=True,
            )

    if df is not None:
        # df.columns = df.columns.str.replace(' ', '')
        df.rename(columns=lambda x: x.strip(), inplace=True)
    else:
        df = pd.DataFrame()

    do_cache(df, file, suffix)
    # print("Read {} rows for file {}".format(len(df), file))

    return df


def read_multi_section_separated_value_file(
    file: str,
    separator: Union[str, None] = None,
    use_cache: bool = True,
    comment: Union[str, None] = None,
    skip_err: bool = False,
) -> pd.DataFrame:
    """
    Reads a file with multiple sections that have different headers throughout the file.

    This function handles files with multiple data sections that may have different column structures.
    It automatically detects the separator if not specified and combines all sections into a single DataFrame.

    :param file: Path to the file to read
    :type file: str
    :param separator: Character used to separate values in the file (auto-detected if None)
    :type separator: Union[str, None]
    :param use_cache: Whether to use cached results if available
    :type use_cache: bool
    :param comment: Character indicating the start of comments (e.g., "#")
    :type comment: Union[str, None]
    :param skip_err: Whether to skip errors when processing sections
    :type skip_err: bool
    :return: Combined DataFrame containing all data from the file
    :rtype: pd.DataFrame

    .. note::
       Particularly useful for reading output files from the DPDSimulator or other
       multi-section files where the structure may change between sections.
    """
    suffix = (
        (
            base64.urlsafe_b64encode(comment.encode("utf-8")).decode("utf-8")
            if comment is not None
            else ""
        )
        + "mssv2-"
        + base64.urlsafe_b64encode(separator.encode("utf-8")).decode("utf-8")
        if separator is not None
        else "-any"
    )
    cache_content = load_cache(file, suffix)

    if cache_content is not None and use_cache:
        return cache_content

    if separator is None:
        # detect separator
        with open(file) as f:
            first_line = f.readline().strip("\n")
        possible_separators = [",", ";", " ", "\t"]
        best_sep = " "
        best_sep_count = 0
        for sep in possible_separators:
            if first_line.count(sep) > best_sep_count:
                best_sep_count = first_line.count(sep)
                best_sep = sep
        separator = best_sep

    print("Splitting CSV...")

    tmp_csv_files = split_csv(file, separator)
    print(
        "CSV split to {} files... e.g. to {}, {} or {}".format(
            len(tmp_csv_files),
            tmp_csv_files[0],
            tmp_csv_files[1] if len(tmp_csv_files) > 1 else "",
            tmp_csv_files[2] if len(tmp_csv_files) > 2 else "",
        )
    )

    if len(tmp_csv_files) == 0:
        return pd.DataFrame()

    # determine the columns we want to have in the end
    all_headers = set()
    detected_dtypes = {}
    erronous_files = []
    for csv_file in tmp_csv_files:
        header_line = ""
        first_line = ""
        got_err = False
        with open(csv_file, "r") as fp:
            try:
                header_line = next(fp)
                first_line = next(fp)
            except StopIteration:
                erronous_files.append(csv_file)
                got_err = True
        if got_err:
            continue
        headers = re.split("{}+".format(separator), header_line.strip())
        if sum([_is_numeric_string(h) for h in headers]) > 0.5 * len(headers):
            warnings.warn(
                "CSV file {} has header line {}, which does not seem to be a header.".format(
                    csv_file, header_line
                )
            )
        for i, h in enumerate(headers):
            if h not in all_headers:
                first_line_split = re.split(
                    "{}+".format(separator), first_line.strip())
                if len(first_line_split) != len(headers):
                    raise ValueError(
                        "Headers and first line do not match in nr of values",
                        first_line,
                        header_line,
                    )
                if all([c.isdigit() or c == "-" for c in first_line_split[i]]):
                    detected_dtypes[h] = np.int64
                elif all(
                    [
                        c.isdigit() or c == "-" or c == "." or c == "e" or c == "E"
                        for c in first_line_split[i]
                    ]
                ):
                    detected_dtypes[h] = np.float64
                all_headers.add(h)
    all_headers = list(all_headers)
    csv_file_to_write = "{}/{}_{}".format(
        tempfile.gettempdir(),
        hashlib.md5(
            datetime.now().strftime("%m.%d.%Y, %H:%M:%S.%f").encode()
        ).hexdigest(),
        "tmp_mssv2_file.csv",
    )

    print("{} Headers mapped...".format(len(all_headers)))

    # re-join the CSV files in one big file with all the columns
    # put NaN where we do not have a value for a column
    with open(csv_file_to_write, "w") as out_file:
        out_file.write(separator.join(all_headers) + "\n")
        for csv_file in tmp_csv_files:
            if csv_file in erronous_files:
                print("File {} skipped".format(csv_file))
                try:
                    os.remove(csv_file)
                except OSError as e:
                    warnings.warn(
                        "Could not remove file {}: {}".format(
                            csv_file, e))
                    pass
                continue
            with open(csv_file, "r") as fp:
                header_line = next(fp)
                split_header = re.split(
                    "{}+".format(separator), header_line.strip())
                map_to_col = []
                n_found = 0
                for i, col in enumerate(all_headers):
                    if col in split_header:
                        map_to_col.append(split_header.index(col))
                        n_found += 1
                    else:
                        map_to_col.append(-1)
                assert n_found == len(split_header)
                for line in fp:
                    if line == header_line or line.startswith("Step"):
                        continue
                    split_line = re.split(
                        "{}+".format(separator), line.strip())
                    str_to_write = separator.join(
                        [split_line[i] if i != -1 else "NaN" for i in map_to_col]
                    )
                    out_file.write(str_to_write + "\n")
            try:
                os.remove(csv_file)
            except OSError as e:
                warnings.warn(
                    "Could not remove file {}: {}".format(
                        csv_file, e))
                pass
            print("File {} handled".format(csv_file))
    # read the csv files again
    print("Reading final csv file {}".format(csv_file_to_write))
    try:
        df = pd.read_csv(
            csv_file_to_write,
            sep=separator + "+",
            comment=comment,
            dtype=detected_dtypes,
            na_values=["NaN"],
        )
    except pd.errors.EmptyDataError:
        warnings.warn("Data file '{}' turned out to be empty".format(file))
        return pd.DataFrame()
    do_cache(df, file, suffix)
    try:
        os.remove(csv_file_to_write)
    except OSError as e:
        warnings.warn(
            "Could not remove file {}: {}".format(
                csv_file_to_write, e))
        pass
    # doCache(reduce_mem_usage(df), file, suffix)
    # print("Read {} rows for file {}".format(len(df), file))

    return df
