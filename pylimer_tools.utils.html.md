# pylimer_tools.utils package

## Submodules

## pylimer_tools.utils.cache_utility module

### pylimer_tools.utils.cache_utility.do_cache(obj, file: [str](https://docs.python.org/3/library/stdtypes.html#str), suffix: [str](https://docs.python.org/3/library/stdtypes.html#str), tmp_dir: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Store the object in the cache.

* **Parameters:**
  * **obj** – The object to cache.
  * **file** – A part of what’s used for the cache’s name. Ideally the file that is read,
    such that the filemtime of file can be used to check whether cache must be generated anew.
  * **suffix** – The file name’s suffix.
  * **tmp_dir** – The directory to store the cache in.

### pylimer_tools.utils.cache_utility.get_cache_file_name(file: [str](https://docs.python.org/3/library/stdtypes.html#str) | [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None), suffix: [str](https://docs.python.org/3/library/stdtypes.html#str), tmp_dir: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, old: [bool](https://docs.python.org/3/library/functions.html#bool) = False)

Get the name and path of a cache file. Internal method.

* **Parameters:**
  * **file** – A part of what’s used for the cache’s name. Ideally the file that is read,
    such that the filemtime of file can be used to check whether cache must be generated anew.
  * **suffix** – The file name’s suffix.
  * **tmp_dir** – The temporary directory.
  * **old** – Whether to use the old file naming scheme.
* **Returns:**
  The path to the cache file.

### pylimer_tools.utils.cache_utility.is_current_cache(cache_file: [str](https://docs.python.org/3/library/stdtypes.html#str), dependencies: [str](https://docs.python.org/3/library/stdtypes.html#str) | [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)])

Determine whether the provided file is newer than all its dependencies.

* **Parameters:**
  * **cache_file** – The cache file that is required to be newer.
  * **dependencies** – The list of files (or a single file path) that need to be older.
* **Returns:**
  True if the file is newer than all its dependencies, False otherwise.

### pylimer_tools.utils.cache_utility.load_cache(file: [str](https://docs.python.org/3/library/stdtypes.html#str) | [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None), suffix: [str](https://docs.python.org/3/library/stdtypes.html#str), disable_warnings: [bool](https://docs.python.org/3/library/functions.html#bool) = False, tmp_dir: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, anyway: [bool](https://docs.python.org/3/library/functions.html#bool) = False)

Load an object from cache, iff the cache is new enough.

* **Parameters:**
  * **file** – A part of what’s used for the cache’s name. Ideally the file that is read,
    such that the filemtime of file can be used to check whether cache must be generated anew.
  * **suffix** – The file name’s suffix.
  * **disable_warnings** – Whether to disable warnings about missing possibilities to check for filemtime.
  * **tmp_dir** – The directory to load the cache from.
  * **anyway** – Whether to ignore the cache’s modification time, and return the cached data anyway,
    as if it were current.
* **Returns:**
  Either the content of the cache, or None if the cache has to be loaded again / is non existent.

## pylimer_tools.utils.data_utility module

### pylimer_tools.utils.data_utility.apply_to_numeric_ignore_compat(df: [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Convert DataFrame columns to numeric when possible.

This preserves the historical behavior of `errors='ignore'` from
`pandas.to_numeric` while remaining compatible with newer pandas versions.

* **Parameters:**
  **df** (*pd.DataFrame*) – DataFrame whose columns should be converted where possible
* **Returns:**
  DataFrame with numeric-like columns converted
* **Return type:**
  pd.DataFrame

### pylimer_tools.utils.data_utility.get_tail(data, percentage=0.2, min_n=25, max_percentage=0.5)

Extract the last few entries of a list

* **Parameters:**
  * **data** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *or* *pd.DataFrame* *or* *pd.Series*) – The list, DataFrame, or Series to extract the last few entries from
  * **percentage** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The percentage of entries to extract (default: 0.2)
  * **min_n** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The minimum number of entries to extract (default: 25)
  * **max_percentage** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The maximum percentage of entries to extract (default: 0.5)
* **Returns:**
  A subset of the input data containing the last entries according to the specified criteria
* **Return type:**
  Same type as input data

The function returns a subset with at maximum max_percentage,
at least min_n entries (assuming the initial data is as large),
but ideally percentage many percentage of the last entries.

### pylimer_tools.utils.data_utility.unify_data_stepsizes(data: [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame), key: [str](https://docs.python.org/3/library/stdtypes.html#str), step_size: [int](https://docs.python.org/3/library/functions.html#int) = None, max_expected_step_size: [int](https://docs.python.org/3/library/functions.html#int) = 100) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Get a DataFrame where all data points have the same step between the values in column given by key

* **Parameters:**
  * **data** (*pd.DataFrame*) – The DataFrame to unify the step-size for
  * **key** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – The column name indicating the column containing the step-nr
  * **step_size** ([*int*](https://docs.python.org/3/library/functions.html#int) *,* *optional*) – The step size to use for filtering (if None, computed automatically)
  * **max_expected_step_size** ([*int*](https://docs.python.org/3/library/functions.html#int) *,* *default=100*) – Used to get a warning if the computed step-size is larger
* **Returns:**
  A DataFrame with a consistent step-size
* **Return type:**
  pd.DataFrame

NOTE: this function is rather unstable, as it has a few assumptions:
- steps are modulo stepsize. Breaks e.g. with steps start with 1 and go up by step_size.
- ideal step-size is max step difference. Breaks e.g. if there is one big gap

## pylimer_tools.utils.optimize_dataframe module

Utility functions to reduce the memory usage of a pandas DataFrame.
Particularly useful when dealing with large datasets, e.g. output from long LAMMPS simulation runs.

Heavily inspired by the following sources:
- [https://medium.com/bigdatarepublic/advanced-pandas-optimize-speed-and-memory-a654b53be6c2](https://medium.com/bigdatarepublic/advanced-pandas-optimize-speed-and-memory-a654b53be6c2)
- [https://stackoverflow.com/questions/57531388/how-can-i-reduce-the-memory-of-a-pandas-dataframe](https://stackoverflow.com/questions/57531388/how-can-i-reduce-the-memory-of-a-pandas-dataframe)

### pylimer_tools.utils.optimize_dataframe.optimize(df: [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame), datetime_features: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] = [])

Optimize all types of all columns in a dataframe.

* **Parameters:**
  * **df** (*pd.DataFrame*) – dataframe to reduce
  * **datetime_features** (*List* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *]*) – list of column names that contain datetime data
* **Returns:**
  dataset with the column dtypes adjusted
* **Return type:**
  pd.DataFrame

### pylimer_tools.utils.optimize_dataframe.optimize_floats(df: [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Optimize the floating point type entries.

* **Parameters:**
  **df** (*pd.DataFrame*) – dataframe to reduce
* **Returns:**
  dataset with the column dtypes adjusted
* **Return type:**
  pd.DataFrame

### pylimer_tools.utils.optimize_dataframe.optimize_ints(df: [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Optimize the integer point type entries.

* **Parameters:**
  **df** (*pd.DataFrame*) – dataframe to reduce
* **Returns:**
  dataset with the column dtypes adjusted
* **Return type:**
  pd.DataFrame

### pylimer_tools.utils.optimize_dataframe.optimize_objects(df: [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame), datetime_features: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)]) → [DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html#pandas.DataFrame)

Optimize object type entries.

* **Parameters:**
  * **df** (*pd.DataFrame*) – dataframe to reduce
  * **datetime_features** (*List* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *]*) – list of column names that contain datetime data
* **Returns:**
  dataset with the column dtypes adjusted
* **Return type:**
  pd.DataFrame

### pylimer_tools.utils.optimize_dataframe.reduce_mem_usage(df, obj_to_category=False, subset=None, inplace=True, print_stats=False)

Iterate through all the columns of a dataframe and modify the data type to reduce memory usage.

* **Parameters:**
  * **df** (*pd.DataFrame*) – dataframe to reduce
  * **obj_to_category** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – convert non-datetime related objects to category dtype
  * **subset** (*List* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *] or* *None*) – subset of columns to analyse
  * **inplace** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – whether to modify the dataframe in place
  * **print_stats** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – whether to print memory usage statistics
* **Returns:**
  dataset with the column dtypes adjusted
* **Return type:**
  pd.DataFrame

## Module contents
