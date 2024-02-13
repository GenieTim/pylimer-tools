
import pandas as pd
from pylimer_tools.utils.cacheUtility import doCache, loadCache


def read_avg_file(filename):
    cache = loadCache(filename, "my-avg")
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
    doCache(result, filename, "my-avg")
    return result
