"""
Utility functions to reduce the memory usage of a pandas DataFrame.
Particularly useful when dealing with large datasets, e.g. output from long LAMMPS simulation runs.

Heavily inspired by the following sources:
- https://medium.com/bigdatarepublic/advanced-pandas-optimize-speed-and-memory-a654b53be6c2
- https://stackoverflow.com/questions/57531388/how-can-i-reduce-the-memory-of-a-pandas-dataframe
"""

import gc
from typing import List

import numpy as np
import pandas as pd


def reduce_mem_usage(
    df, obj_to_category=False, subset=None, inplace=True, print_stats=False
):
    """
    Iterate through all the columns of a dataframe and modify the data type to reduce memory usage.

    :param df: dataframe to reduce
    :type df: pd.DataFrame
    :param obj_to_category: convert non-datetime related objects to category dtype
    :type obj_to_category: bool
    :param subset: subset of columns to analyse
    :type subset: List[str] or None
    :param inplace: whether to modify the dataframe in place
    :type inplace: bool
    :param print_stats: whether to print memory usage statistics
    :type print_stats: bool

    :return: dataset with the column dtypes adjusted
    :rtype: pd.DataFrame
    """
    start_mem = df.memory_usage().sum() / 1024**2
    gc.collect()
    if print_stats:
        print("Memory usage of dataframe is {:.2f} MB".format(start_mem))

    if not inplace:
        df = df.copy()

    cols = subset if subset is not None else df.columns.tolist()

    for col in cols:
        col_type = df[col].dtype

        if (
            col_type != object
            and col_type.name != "category"
            and "datetime" not in col_type.name
        ):
            c_min = df[col].min()
            c_max = df[col].max()

            # test if column can be converted to an integer
            treat_as_int = str(col_type)[
                :3] == "int" or str(col_type)[
                :4] == "uint"

            if treat_as_int:
                if c_min > np.iinfo(
                        np.int8).min and c_max < np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif c_min > np.iinfo(np.uint8).min and c_max < np.iinfo(np.uint8).max:
                    df[col] = df[col].astype(np.uint8)
                elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif (
                    c_min > np.iinfo(
                        np.uint16).min and c_max < np.iinfo(
                        np.uint16).max
                ):
                    df[col] = df[col].astype(np.uint16)
                elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
                elif (
                    c_min > np.iinfo(
                        np.uint32).min and c_max < np.iinfo(
                        np.uint32).max
                ):
                    df[col] = df[col].astype(np.uint32)
                elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                    df[col] = df[col].astype(np.int64)
                elif (
                    c_min > np.iinfo(
                        np.uint64).min and c_max < np.iinfo(
                        np.uint64).max
                ):
                    df[col] = df[col].astype(np.uint64)
            else:
                if (
                    c_min > np.finfo(np.float16).min
                    and c_max < np.finfo(np.float16).max
                ):
                    df[col] = df[col].astype(np.float16)
                elif (
                    c_min > np.finfo(np.float32).min
                    and c_max < np.finfo(np.float32).max
                ):
                    df[col] = df[col].astype(np.float32)
                else:
                    df[col] = df[col].astype(np.float64)
        elif "datetime" not in col_type.name and obj_to_category:
            df[col] = df[col].astype("category")
    gc.collect()
    end_mem = df.memory_usage().sum() / 1024**2
    if print_stats:
        print("Memory usage after optimization is: {:.3f} MB".format(end_mem))
        print("Decreased by {:.1f}%".format(
            100 * (start_mem - end_mem) / start_mem))

    return df


def optimize_floats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Optimize the floating point type entries.

    :param df: dataframe to reduce
    :type df: pd.DataFrame

    :return: dataset with the column dtypes adjusted
    :rtype: pd.DataFrame
    """
    floats = df.select_dtypes(include=["float64"]).columns.tolist()
    df[floats] = df[floats].apply(pd.to_numeric, downcast="float")
    return df


def optimize_ints(df: pd.DataFrame) -> pd.DataFrame:
    """
    Optimize the integer point type entries.

    :param df: dataframe to reduce
    :type df: pd.DataFrame

    :return: dataset with the column dtypes adjusted
    :rtype: pd.DataFrame
    """
    ints = df.select_dtypes(include=["int64"]).columns.tolist()
    df[ints] = df[ints].apply(pd.to_numeric, downcast="integer")
    return df


def optimize_objects(df: pd.DataFrame,
                     datetime_features: List[str]) -> pd.DataFrame:
    """
    Optimize object type entries.

    :param df: dataframe to reduce
    :type df: pd.DataFrame
    :param datetime_features: list of column names that contain datetime data
    :type datetime_features: List[str]

    :return: dataset with the column dtypes adjusted
    :rtype: pd.DataFrame
    """
    for col in df.select_dtypes(include=["object"]):
        if col not in datetime_features:
            num_unique_values = len(df[col].unique())
            num_total_values = len(df[col])
            if float(num_unique_values) / num_total_values < 0.5:
                df[col] = df[col].astype("category")
        else:
            df[col] = pd.to_datetime(df[col])
    return df


def optimize(df: pd.DataFrame, datetime_features: List[str] = []):
    """
    Optimize all types of all columns in a dataframe.

    :param df: dataframe to reduce
    :type df: pd.DataFrame
    :param datetime_features: list of column names that contain datetime data
    :type datetime_features: List[str]

    :return: dataset with the column dtypes adjusted
    :rtype: pd.DataFrame
    """
    return optimize_floats(optimize_ints(
        optimize_objects(df, datetime_features)))
