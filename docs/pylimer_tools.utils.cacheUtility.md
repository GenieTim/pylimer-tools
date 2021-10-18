<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/utils/cacheUtility.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.utils.cacheUtility`





---

<a href="../src/pylimer_tools/utils/cacheUtility.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `doCache`

```python
doCache(obj, file: str, suffix: str)
```

Store the object in the cache 



**Arguments:**
 
    - obj: the object to cache 
    - file: the path to the file to save the object to 
    - suffix: the file name's suffix 


---

<a href="../src/pylimer_tools/utils/cacheUtility.py#L24"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `loadCache`

```python
loadCache(file: str, suffix: str, disableWarnings: bool = False)
```

Load an object from cache. 



**Arguments:**
 
    - file: a cache name. Ideally the file that is read, such that the filemtime of `file` can be used to check whether cache must be generated anew 
    - suffix: the file name's suffix 
    - disableWarnings: whether to disable warnings about missing possibilities to check for filemtime 



**Returns:**
 
    - cache: either the content of the cache, or None if the cache has to be loaded again / is non existant 


---

<a href="../src/pylimer_tools/utils/cacheUtility.py#L61"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getCacheFileName`

```python
getCacheFileName(file: str, suffix: str)
```

Get the name and path of a cache file. Internal method. 



**Arguments:**
 
    - file: a cache name. Ideally the file that is read. 
    - suffix: the file name's suffix 



**Returns:**
 
    - cacheFileName: the path to the cache file 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
