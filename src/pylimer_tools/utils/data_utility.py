import warnings
import pandas as pd


def get_tail(data, percentage=0.2, min_n=25, max_percentage=0.5):
    """
    Extract the last few entries of a list

    Arguments:
        - data (list|pd.DataFrame|pd.Series): the list to extract the last few entries from
        - percentage: the percentage of entries to extract
        - min_n: the minimum number of entries to extract
        - max_percentage: the maximum percentage of entries to extract

    Returns:
        A (list|pd.DataFrame|pd.Series) with at maximum maxPercentage,
            at least minN entries (assuming the initial data is as large),
            but ideally `percentage` many percentage of the last entries.
    """
    assert (percentage <= 1)
    assert (max_percentage <= 1)
    tail_n = int(min(max(min(min_n, max_percentage * len(data)),
                         percentage * len(data)), len(data)))
    if (isinstance(data, pd.DataFrame) or isinstance(data, pd.Series)):
        return data.tail(tail_n)
    else:
        return data[-tail_n:]


def unify_data_stepsizes(data: pd.DataFrame, key: str, step_size: int = None,
                         max_expected_step_size: int = 100) -> pd.DataFrame:
    """
    Get a DataFrame where all data points have the same step between the values in column given by `key`
    NOTE: this function is rather unstable, as it has a few dirty assumptions, such as:
    - steps are modulo stepsize. Breaks e.g. with steps start with 1 and go up by step_size.
    - ideal step-size is max step difference. Breaks e.g. if there is one big gap

    Arguments:
        - data: the DataFrame to unify the step-size for
        - key: the column name indicating the column containing the step-nr.
        - max_expected_step_size: use to get a warning if the computed step-size is larger

    Returns:
        - data: a DataFrame with a consistent step-size
    """
    # lenBefore = len(data)
    if (step_size is None):
        step_size = data[key].sort_values().diff().max()
    if (step_size > max_expected_step_size):
        warnings.warn("Step size {} unexpectedly large, with max expected {}".format(
            step_size, max_expected_step_size))
    data = data[(data[key] % step_size) == 0]
    # print("Reduced from {} to {} data-points using step size of {}".format(lenBefore, len(data), step_size))
    return data
