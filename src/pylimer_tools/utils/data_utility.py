import warnings
from typing import Any, Optional, TypeVar

import pandas as pd

T = TypeVar("T")
DataFrameLike = TypeVar("DataFrameLike")


def apply_to_numeric_ignore_compat(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert DataFrame columns to numeric when possible.

    This preserves the historical behavior of ``errors='ignore'`` from
    ``pandas.to_numeric`` while remaining compatible with newer pandas versions.

    :param df: DataFrame whose columns should be converted where possible
    :type df: pd.DataFrame
    :return: DataFrame with numeric-like columns converted
    :rtype: pd.DataFrame
    """

    def _convert(series: pd.Series) -> pd.Series:
        try:
            return pd.to_numeric(series, errors="raise")
        except (ValueError, TypeError):
            return series

    return df.apply(_convert)


def get_tail(
    data: T,
    percentage: float = 0.2,
    min_n: int = 25,
    max_percentage: float = 0.5,
) -> T:
    """
    Extract the last few entries of a list, DataFrame, or Series.

    :param data: The list, DataFrame, or Series to extract the last few entries from
    :type data: list or pd.DataFrame or pd.Series or pl.DataFrame or pl.Series
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
            max(min(min_n, max_percentage * len(data)), percentage * len(data)),  # type: ignore[arg-type]
            len(data),  # type: ignore[arg-type]
        )
    )
    if hasattr(data, "tail") and callable(data.tail):
        return data.tail(tail_n)  # type: ignore[return-value]
    if tail_n == 0:
        return data[:0]  # type: ignore[index, return-value]
    return data[-tail_n:]  # type: ignore[index, return-value]


def unify_data_stepsizes(
    data: DataFrameLike,
    key: str,
    step_size: Optional[float] = None,
    max_expected_step_size: float = 100,
) -> DataFrameLike:
    """
    Get a DataFrame where all data points have the same step between the values in column given by `key`

    :param data: The DataFrame to unify the step-size for (pandas or polars DataFrame)
    :type data: pd.DataFrame or pl.DataFrame
    :param key: The column name indicating the column containing the step-nr
    :type key: str
    :param step_size: The step size to use for filtering (if None, computed automatically)
    :type step_size: Optional[float]
    :param max_expected_step_size: Used to get a warning if the computed step-size is larger
    :type max_expected_step_size: float, default=100
    :return: A DataFrame with a consistent step-size
    :rtype: Same type as input data

    NOTE: this function is rather unstable, as it has a few assumptions:
    - steps are modulo stepsize. Breaks e.g. with steps start with 1 and go up by step_size.
    - ideal step-size is max step difference. Breaks e.g. if there is one big gap
    """
    is_polars = hasattr(data, "filter") and not isinstance(data, pd.DataFrame)

    if step_size is None:
        if len(data) <= 1:  # type: ignore[arg-type]
            return data
        if is_polars:
            computed: Any = data[key].sort().diff().max()  # type: ignore[index]
        else:
            computed = data[key].sort_values().diff().max()  # type: ignore[index]
        if computed is None or pd.isna(computed):
            return data
        step_size = float(computed)

    if pd.isna(step_size) or step_size == 0:
        return data

    if step_size > max_expected_step_size:
        warnings.warn(
            "Step size {} unexpectedly large, with max expected {}".format(
                step_size, max_expected_step_size
            )
        )
    if is_polars:
        return data.filter((data[key] % step_size) == 0)  # type: ignore[attr-defined, index, return-value]
    return data[(data[key] % step_size) == 0]  # type: ignore[index, return-value]
