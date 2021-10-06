Module pylimer_tools.utils.optimizeDf
=====================================

Functions
---------

    
`optimize(df: pandas.core.frame.DataFrame, datetime_features: List[str] = [])`
:   Optimize all types of all columns in a dataframe
    
    Arguments:
        - df (pd.DataFrame): dataframe to reduce
    
    Returns:
        - df (pd.DataFrame): dataset with the column dtypes adjusted

    
`optimize_floats(df: pandas.core.frame.DataFrame) ‑> pandas.core.frame.DataFrame`
:   Optimize the floating point type entries
    
    Arguments:
        df (pd.DataFrame): dataframe to reduce
    
    Returns:
        df (pd.DataFrame): dataset with the column dtypes adjusted

    
`optimize_ints(df: pandas.core.frame.DataFrame) ‑> pandas.core.frame.DataFrame`
:   Optimize the integer point type entries
    
    Arguments:
        df (pd.DataFrame): dataframe to reduce
    
    Returns:
        df (pd.DataFrame): dataset with the column dtypes adjusted

    
`optimize_objects(df: pandas.core.frame.DataFrame, datetime_features: List[str]) ‑> pandas.core.frame.DataFrame`
:   Optimize object type entries
    
    Arguments:
        - df (pd.DataFrame): dataframe to reduce
    
    Returns:
        - df (pd.DataFrame): dataset with the column dtypes adjusted

    
`reduce_mem_usage(df, obj_to_category=False, subset=None)`
:   Iterate through all the columns of a dataframe and modify the data type to reduce memory usage.
    
    Arguments:
        df (pd.DataFrame): dataframe to reduce
        obj_to_category (boolean): convert non-datetime related objects to category dtype
        subset (List): subset of columns to analyse
    
    Returns:
        df (pd.DataFrame): dataset with the column dtypes adjusted