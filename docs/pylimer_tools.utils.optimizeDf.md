<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/utils/optimizeDf.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.utils.optimizeDf`





---

<a href="../src/pylimer_tools/utils/optimizeDf.py#L11"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `reduce_mem_usage`

```python
reduce_mem_usage(df, obj_to_category=False, subset=None)
```

Iterate through all the columns of a dataframe and modify the data type to reduce memory usage. 



**Arguments:**
 
 - <b>`df`</b> (pd.DataFrame):  dataframe to reduce 
 - <b>`obj_to_category`</b> (boolean):  convert non-datetime related objects to category dtype 
 - <b>`subset`</b> (List):  subset of columns to analyse 



**Returns:**
 
 - <b>`df`</b> (pd.DataFrame):  dataset with the column dtypes adjusted 


---

<a href="../src/pylimer_tools/utils/optimizeDf.py#L74"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `optimize_floats`

```python
optimize_floats(df: DataFrame) → DataFrame
```

Optimize the floating point type entries 



**Arguments:**
 
 - <b>`df`</b> (pd.DataFrame):  dataframe to reduce 



**Returns:**
 
 - <b>`df`</b> (pd.DataFrame):  dataset with the column dtypes adjusted 


---

<a href="../src/pylimer_tools/utils/optimizeDf.py#L89"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `optimize_ints`

```python
optimize_ints(df: DataFrame) → DataFrame
```

Optimize the integer point type entries 



**Arguments:**
 
 - <b>`df`</b> (pd.DataFrame):  dataframe to reduce 



**Returns:**
 
 - <b>`df`</b> (pd.DataFrame):  dataset with the column dtypes adjusted 


---

<a href="../src/pylimer_tools/utils/optimizeDf.py#L104"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `optimize_objects`

```python
optimize_objects(df: DataFrame, datetime_features: List[str]) → DataFrame
```

Optimize object type entries 



**Arguments:**
 
    - df (pd.DataFrame): dataframe to reduce 



**Returns:**
 
    - df (pd.DataFrame): dataset with the column dtypes adjusted 


---

<a href="../src/pylimer_tools/utils/optimizeDf.py#L125"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `optimize`

```python
optimize(df: DataFrame, datetime_features: List[str] = [])
```

Optimize all types of all columns in a dataframe 



**Arguments:**
 
    - df (pd.DataFrame): dataframe to reduce 



**Returns:**
 
    - df (pd.DataFrame): dataset with the column dtypes adjusted 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
