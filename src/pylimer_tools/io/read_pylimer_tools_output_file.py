
import pandas as pd

from pylimer_tools.utils.cache_utility import do_cache, load_cache


def read_avg_file(filename: str):
    """
    Read an averages-output file from one of the simulators shipped with pylimer_tools

    Args:
        - filename: The path to the file to read

    Returns:
        - results (pd.DataFrame): the content of the averages file

    """
    cache = load_cache(filename, "my-avg")
    if (cache is not None):
        return cache
    data_frames = []
    with open(filename, "r") as f:
        first_line_split = f.readline().removeprefix("#").strip().split()
        data = []
        for line in f:
            if ("-nan" in line or '\x00' in line or len(line.split()) < 3):
                continue
            stripped_line = line.removeprefix("#").strip()
            if (stripped_line.startswith(first_line_split[0])):
                data_frames.append(
                    pd.DataFrame(data, columns=first_line_split)
                )
                first_line_split = stripped_line.split()
                data = []
            elif (stripped_line != ""):
                data.append(stripped_line.split())
    if (not len(data) == 0):
        data_frames.append(pd.DataFrame(data, columns=first_line_split))
    df = pd.concat(data_frames, ignore_index=True)
    result = df.apply(pd.to_numeric, errors='ignore')
    result = result.groupby("OutputStep", as_index=False).last()
    assert (not result["OutputStep"].duplicated().any())
    do_cache(result, filename, "my-avg")
    return result
