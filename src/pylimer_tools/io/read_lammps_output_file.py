import os
import re
import warnings
from typing import List, Union

import pandas as pd

from pylimer_tools.io.extract_thermo_data import extract_thermo_params
from pylimer_tools.utils.cache_utility import do_cache, load_cache
from pylimer_tools_cpp import AtomStyle, Universe, UniverseSequence


def read_log_file(filepath, lines_to_read_to_detect_header=500000) -> pd.DataFrame:
    return extract_thermo_params(filepath, header=None, texts_to_read=500000,
                                 lines_to_read_to_detect_header=lines_to_read_to_detect_header)


def read_dump_file(data_file, dump_file, atom_style: Union[List[AtomStyle], None] = None) -> UniverseSequence:
    """
    Read a file with LAMMPS' dump of snapshots of structures into a Universe.
    """
    u_s = UniverseSequence()
    if (atom_style is not None):
        u_s.set_data_file_atom_style(atom_style)
    u_s.initialize_from_dump_file(data_file, dump_file)
    return u_s


def read_data_file(structure_file: str, atom_style: Union[List[AtomStyle], None] = None) -> Universe:
    """
    Read a file with LAMMPS' data type of structure into a Universe.

    Arguments:
        - structure_file: the path to the structure file
        - atom_style: the atom style in the structure file. Defaults to AtomStyle.Molecule

    Returns:
        - universe (Universe): the Universe in the given structure file
    """
    u_s = UniverseSequence()
    if (atom_style is not None):
        u_s.set_data_file_atom_style(atom_style)
    u_s.initialize_from_data_sequence([structure_file])
    universe = u_s.atIndex(0)
    del u_s
    return universe


def read_averages_file(filepath, use_cache: bool = True, sep=" ") -> pd.DataFrame:
    """
    Read a file written by a `fix ave/time` command.

    Uses pandas' read_csv after detecting the columns.

    Important assumption: the first 2 or 3 lines in the file are:
    - comment,
    - then one header indicating the columns,
    - and then either data or potentially a second header, 
      if it is a sectioned file (e.g., from a `fix ave/time ... vector`)

    Arguments:
        - filepath: the path to the averages file
        - use_cache: whether to use the cache to speed up the reading & writing
        - sep: if, for some reason, the file uses a different delimiter than " "
    """
    assert (os.path.isfile(filepath))
    header_line = None
    with open(filepath, 'r') as f:
        line0 = f.readline()
        line1 = f.readline()
        line2 = f.readline()

        if (line2.startswith("#")):
            return read_sectioned_averages_file(filepath, use_cache=use_cache)

        header_line = line1 if line1.startswith("#") else line0
    header_line = header_line.removeprefix("#").strip()

    try:
        data = pd.read_csv(filepath, comment="#",
                           names=header_line.split(), sep=sep)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    return data


def read_sectioned_averages_file(filepath, use_cache: bool = True) -> pd.DataFrame:
    """
    Read a file written by a `fix ave/time` command.

    Use the section delimiter columns together with pandas' groupby()
    to restore the original sections.

    Arguments:
        - filepath: the path to the averages file
        - use_cache: whether to use the cache to speed up the reading & writing
    """
    assert (os.path.isfile(filepath))

    cache_suffix = "sectionedavg-cache.pickle"
    cache_content = load_cache(filepath, cache_suffix)

    if (cache_content is not None and use_cache):
        return cache_content

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

        current_data = []
        current_key = None
        for line in f:
            split_line = line.split()
            if (current_key is None):
                assert (len(split_line) == len(header_line1.split()))
                current_key = line
                continue
            if (len(split_line) == len(header_line1_split)):
                data[current_key] = current_data
                current_data = []
                current_key = line
            else:
                assert (len(split_line) == len(header_line2_split))
                current_data.append(split_line)
        data[current_key] = current_data

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
    do_cache(df, filepath, cache_suffix)

    return df


def read_histogram_file(filepath, use_cache: bool = True) -> pd.DataFrame:
    """
    Read a file written by `fix ave/hist` or similar.

    See also:
        - :func:`~pylimer_tools.io.read_lammps_output_file.readSectionedAveragesFile`
    """
    return read_sectioned_averages_file(filepath, use_cache)


def read_correlation_file(filepath, group_key="Timestep", use_cache: bool = True) -> pd.DataFrame:
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
    cache_content = load_cache(filepath, cache_suffix)

    if (cache_content is not None and use_cache):
        return cache_content

    data = {}
    header_line = None
    with open(filepath, 'r') as f:
        current_data = []
        current_key = None
        header_line = f.readline()
        if (header_line.startswith("#")):
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
            if ((line.startswith("#") or len(line.strip()) == 0) and not is_group_key(line)):
                if (lines_interpreted == 0):
                    header_line = line
                continue
            if (line == header_line):
                continue
            split = line.removeprefix("#").strip().split()
            if (len(split) == 2 or is_group_key(line)):
                if (current_key is not None and len(current_data) > 0):
                    data[current_key] = current_data
                    current_data = []
                # new key
                current_key = line
            elif (len(split) == normal_line_len or normal_line_len is None):
                # normal_line_len = len(split)
                current_data.append(split)
            else:
                raise ValueError(
                    "Did not expect {} splited values on line with content {} in correlation file {}".format(
                        len(split), line, filepath)
                )
            lines_interpreted += 1
        if (current_key is not None and len(current_data) > 0):
            data[current_key] = current_data

    cols.append(group_key)
    correlated_data_assembled = []
    for key in data.keys():
        assert (group_key in str(key))
        compiled_regex = re.compile(r"{}:? ([\d]+)".format(group_key))
        results = compiled_regex.search(key)
        if (results is None):
            warnings.warn("Did not find {} with number in {} when reading {}".format(
                group_key, key, filepath))
        assert (results is not None)
        timestep = int(results.group(1))
        for row in data[key]:
            row.append(timestep)
            assert (len(row) == len(cols))
            correlated_data_assembled.append(row)

    correlated_data = pd.DataFrame(
        correlated_data_assembled, columns=cols)
    # convert all columns of DataFrame
    correlated_data = correlated_data.apply(pd.to_numeric, errors='ignore')
    do_cache(correlated_data, filepath, cache_suffix)

    return correlated_data
