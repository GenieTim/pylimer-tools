import warnings

import pandas as pd


def get_tail(data, percentage=0.2, min_n=25, max_percentage=0.5):
    """
    Extract the last few entries of a list

    :param data: The list, DataFrame, or Series to extract the last few entries from
    :type data: list or pd.DataFrame or pd.Series
    :param percentage: The percentage of entries to extract (default: 0.2)
    :type percentage: float
    :param min_n: The minimum number of entries to extract (default: 25)
    :type min_n: int
    :param max_percentage: The maximum percentage of entries to extract (default: 0.5)
    :type max_percentage: float
    :return: A subset of the input data containing the last entries according to the specified criteria
    :rtype: Same type as input data

    The function returns a subset with at maximum max_percentage,
    at least min_n entries (assuming the initial data is as large),
    but ideally `percentage` many percentage of the last entries.
    """
    assert percentage <= 1
    assert max_percentage <= 1
    tail_n = int(
        min(
            max(min(min_n, max_percentage * len(data)), percentage * len(data)),
            len(data),
        )
    )
    if isinstance(data, pd.DataFrame) or isinstance(data, pd.Series):
        return data.tail(tail_n)
    else:
        return data[-tail_n:]


def unify_data_stepsizes(
    data: pd.DataFrame,
    key: str,
    step_size: int = None,
    max_expected_step_size: int = 100,
) -> pd.DataFrame:
    """
    Get a DataFrame where all data points have the same step between the values in column given by `key`

    :param data: The DataFrame to unify the step-size for
    :type data: pd.DataFrame
    :param key: The column name indicating the column containing the step-nr
    :type key: str
    :param step_size: The step size to use for filtering (if None, computed automatically)
    :type step_size: int, optional
    :param max_expected_step_size: Used to get a warning if the computed step-size is larger
    :type max_expected_step_size: int, default=100
    :return: A DataFrame with a consistent step-size
    :rtype: pd.DataFrame

    NOTE: this function is rather unstable, as it has a few assumptions:
    - steps are modulo stepsize. Breaks e.g. with steps start with 1 and go up by step_size.
    - ideal step-size is max step difference. Breaks e.g. if there is one big gap
    """
    # lenBefore = len(data)
    if step_size is None:
        step_size = data[key].sort_values().diff().max()
    if step_size > max_expected_step_size:
        warnings.warn(
            "Step size {} unexpectedly large, with max expected {}".format(
                step_size, max_expected_step_size
            )
        )
    data = data[(data[key] % step_size) == 0]
    # print("Reduced from {} to {} data-points using step size of {}".format(lenBefore, len(data), step_size))
    return data
