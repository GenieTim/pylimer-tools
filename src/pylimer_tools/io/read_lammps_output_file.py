"""
This module provides a few functions to read LAMMPS' output files, including:

- log files (thermo output)
- dump files (focusing on the coordinates of atoms)
- data files (the LAMMPS structure)
- averaged data (from :code:`fix ave/time...` or :code:`fix ave/hist...`)
- correlation data (from :code:`fix ave/correlate/...`)

"""

import os
import re
import warnings
from typing import List, Union

import pandas as pd

from pylimer_tools.io.extract_thermo_data import extract_thermo_params
from pylimer_tools.utils.cache_utility import do_cache, load_cache
from pylimer_tools_cpp import AtomStyle, Universe, UniverseSequence


def read_log_file(
        filepath, lines_to_read_to_detect_header=500000) -> pd.DataFrame:
    """
    Read a LAMMPS' log (thermo output) file.

    :param filepath: Path to the LAMMPS log file
    :type filepath: str
    :param lines_to_read_to_detect_header: Maximum number of lines to read when detecting the header
    :type lines_to_read_to_detect_header: int
    :return: DataFrame containing the parsed thermo data
    :rtype: pd.DataFrame
    """
    return extract_thermo_params(
        filepath,
        header=None,
        texts_to_read=500000,
        lines_to_read_to_detect_header=lines_to_read_to_detect_header,
    )


def read_dump_file(
    data_file, dump_file, atom_style: Union[List[AtomStyle], None] = None
) -> UniverseSequence:
    """
    Read a file with LAMMPS' dump of snapshots of structures into a Universe.

    :param data_file: Path to the LAMMPS data file containing structure information
    :type data_file: str
    :param dump_file: Path to the LAMMPS dump file containing trajectory information
    :type dump_file: str
    :param atom_style: The atom style(s) used in the data file
    :type atom_style: Union[List[AtomStyle], None]
    :return: Sequence of Universe objects representing the trajectory
    :rtype: UniverseSequence
    """
    u_s = UniverseSequence()
    if atom_style is not None:
        u_s.set_data_file_atom_style(atom_style)
    u_s.initialize_from_dump_file(data_file, dump_file)
    return u_s


def read_data_file(
    structure_file: str, atom_style: Union[List[AtomStyle], None] = None
) -> Universe:
    """
    Read a file with LAMMPS' data type of structure into a Universe.

    :param structure_file: Path to the structure file
    :type structure_file: str
    :param atom_style: The atom style(s) in the structure file (defaults to AtomStyle.Molecule if None)
    :type atom_style: Union[List[AtomStyle], None]
    :return: Universe object representing the molecular structure
    :rtype: Universe
    :raises FileNotFoundError: If the structure file does not exist
    """
    if not (os.path.isfile(structure_file)):
        raise FileNotFoundError(
            f"Structure-file '{structure_file}' not found.")
    u_s = UniverseSequence()
    if atom_style is not None:
        u_s.set_data_file_atom_style(atom_style)
    u_s.initialize_from_data_sequence([structure_file])
    universe = u_s.at_index(0)
    del u_s
    return universe


def read_averages_file(filepath, use_cache: bool = True,
                       sep=" ") -> pd.DataFrame:
    """
    Read a file written by a `fix ave/time` command.

    Uses pandas' read_csv after detecting the columns.

    Important assumption: The first 2 or 3 lines in the file are:
        - comment,
        - then one header indicating the columns,
        - and then either data or potentially a second header, if it is a sectioned file (e.g., from a `fix ave/time ... vector`)

    :param filepath: Path to the averages file
    :type filepath: str
    :param use_cache: Whether to use the cache to speed up reading & writing
    :type use_cache: bool
    :param sep: Delimiter used in the file (default is space)
    :type sep: str
    :return: DataFrame containing the parsed average data
    :rtype: pd.DataFrame
    :raises FileNotFoundError: If the averages file does not exist
    """
    if not (os.path.isfile(filepath)):
        raise FileNotFoundError(f"Averages-file '{filepath}' not found.")
    header_line = None
    with open(filepath, "r") as f:
        line0 = f.readline()
        line1 = f.readline()
        line2 = f.readline()

        if line2.startswith("#"):
            return read_sectioned_averages_file(filepath, use_cache=use_cache)

        header_line = line1 if line1.startswith("#") else line0
    header_line = header_line.removeprefix("#").strip()

    try:
        data = pd.read_csv(
            filepath,
            comment="#",
            names=header_line.split(),
            sep=sep)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    return data


def read_sectioned_averages_file(
        filepath, use_cache: bool = True) -> pd.DataFrame:
    """
    Read a file written by a `fix ave/time` command with multiple sections.

    Use the section delimiter columns together with pandas' groupby()
    to restore the original sections.

    :param filepath: Path to the sectioned averages file
    :type filepath: str
    :param use_cache: Whether to use the cache to speed up reading & writing
    :type use_cache: bool
    :return: DataFrame containing the parsed sectioned data
    :rtype: pd.DataFrame
    :raises FileNotFoundError: If the file does not exist
    :raises ValueError: If the file format is not recognized as a proper sectioned averages file
    """
    if not (os.path.isfile(filepath)):
        raise FileNotFoundError(f"Averages-file '{filepath}' not found.")

    cache_suffix = "sectionedavg-cache.pickle"
    cache_content = load_cache(filepath, cache_suffix)

    if cache_content is not None and use_cache:
        return cache_content

    data = {}
    with open(filepath, "r") as f:
        f.readline()  # discard line 0
        line1 = f.readline()
        line2 = f.readline()

        if not line2.startswith("#"):
            raise ValueError(
                "The file '{}' was not detected to be a proper sectioned averages file.".format(
                    filepath
                )
            )
            # return readSectionedAveragesFile(filepath)

        header_line1 = line1.removeprefix("#").strip()
        header_line2 = line2.removeprefix("#").strip()

        header_line1_split = header_line1.split()
        header_line2_split = header_line2.split()

        if len(header_line1_split) == len(header_line2_split):
            raise ValueError(
                "Cannot read this file, as we cannot distinguish between section header and main data"
            )

        current_data = []
        current_key = None
        for line in f:
            split_line = line.split()
            if current_key is None:
                assert len(split_line) == len(header_line1.split())
                current_key = line
                continue
            if len(split_line) == len(header_line1_split):
                data[current_key] = current_data
                current_data = []
                current_key = line
            else:
                assert len(split_line) == len(header_line2_split)
                current_data.append(split_line)
        data[current_key] = current_data

    # convert all the data to a dataframe
    dfs_to_concat = []

    if header_line1_split is None:
        raise ValueError("Did not find a useable header line.")

    for key in data.keys():
        split_key = key.split()
        local_dataframe = pd.DataFrame(data[key], columns=header_line2_split)
        for i, col in enumerate(header_line1_split):
            local_dataframe[col] = split_key[i]
        dfs_to_concat.append(local_dataframe)

    df = pd.concat(dfs_to_concat, ignore_index=True)

    # convert all columns of DataFrame
    df = df.apply(pd.to_numeric, errors="ignore")
    do_cache(df, filepath, cache_suffix)

    return df


def read_histogram_file(filepath, use_cache: bool = True) -> pd.DataFrame:
    """
    Read a file written by `fix ave/hist` or similar.

    This is a wrapper around read_sectioned_averages_file for histogram data.

    :param filepath: Path to the histogram file
    :type filepath: str
    :param use_cache: Whether to use the cache to speed up reading & writing
    :type use_cache: bool
    :return: DataFrame containing the parsed histogram data
    :rtype: pd.DataFrame

    :see: :func:`~pylimer_tools.io.read_lammps_output_file.read_sectioned_averages_file`
    """
    return read_sectioned_averages_file(filepath, use_cache)


def read_correlation_file(
    filepath, group_key="Timestep", use_cache: bool = True
) -> pd.DataFrame:
    """
    Read a file written by a `fix ave/correlate{/long}` command.

    :param filepath: Path to the correlation file
    :type filepath: str
    :param group_key: The key that denotes a new section
    :type group_key: str
    :param use_cache: Whether to use the cache to speed up reading & writing
    :type use_cache: bool
    :return: DataFrame containing the correlation data. Use the group_key with
             the DataFrame's groupby() to restore the original sections.
    :rtype: pd.DataFrame
    :raises FileNotFoundError: If the correlation file does not exist
    """
    if not (os.path.isfile(filepath)):
        raise FileNotFoundError(f"Correlation-file '{filepath}' not found.")

    cache_suffix = "{}-correlation-cache.pickle".format(
        group_key if isinstance(group_key, str) else "g"
    )
    cache_content = load_cache(filepath, cache_suffix)

    if cache_content is not None and use_cache:
        return cache_content

    data = {}
    header_line = None
    with open(filepath, "r") as f:
        current_data = []
        current_key = None
        header_line = f.readline()
        if header_line.startswith("#"):
            # in LAMMPS files, there is a title line that does not exist in our DPD output,
            # -> this line is needed for LAMMPS
            header_line = f.readline()
        cols = header_line.removeprefix("#").strip().split()
        normal_line_len = len(cols)
        lines_interpreted = 0

        def is_group_key(line):
            # if (isinstance(group_key, list)):
            #     return np.any([x in line for x in group_key])
            # else:
            return group_key in line

        for line in f:
            if (line.startswith("#") or len(line.strip()) == 0) and not is_group_key(
                line
            ):
                if lines_interpreted == 0:
                    header_line = line
                continue
            if line == header_line:
                continue
            split = line.removeprefix("#").strip().split()
            if len(split) == 2 or is_group_key(line):
                if current_key is not None and len(current_data) > 0:
                    data[current_key] = current_data
                    current_data = []
                # new key
                current_key = line
            elif len(split) == normal_line_len or normal_line_len is None:
                # normal_line_len = len(split)
                current_data.append(split)
            else:
                raise ValueError(
                    "Did not expect {} splited values on line with content {} in correlation file {}".format(
                        len(split), line, filepath
                    )
                )
            lines_interpreted += 1
        if current_key is not None and len(current_data) > 0:
            data[current_key] = current_data

    cols.append(group_key)
    correlated_data_assembled = []
    for key in data.keys():
        assert group_key in str(key)
        compiled_regex = re.compile(r"{}:? ([\d]+)".format(group_key))
        results = compiled_regex.search(key)
        if results is None:
            warnings.warn(
                "Did not find {} with number in {} when reading {}".format(
                    group_key, key, filepath
                )
            )
        assert results is not None
        timestep = int(results.group(1))
        for row in data[key]:
            row.append(timestep)
            assert len(row) == len(cols)
            correlated_data_assembled.append(row)

    correlated_data = pd.DataFrame(correlated_data_assembled, columns=cols)
    # convert all columns of DataFrame
    correlated_data = correlated_data.apply(pd.to_numeric, errors="ignore")
    do_cache(correlated_data, filepath, cache_suffix)

    return correlated_data
