<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/io/readLammpData.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.io.readLammpData`





---

<a href="../src/pylimer_tools/io/readLammpData.py#L7"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `findFloatsInLine`

```python
findFloatsInLine(line: str) → list
```

Find as many floats as possible in a string (line) 


---

<a href="../src/pylimer_tools/io/readLammpData.py#L14"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `findIntsInLine`

```python
findIntsInLine(line: str) → list
```

Find as many integers as possible in a string (line) 


---

<a href="../src/pylimer_tools/io/readLammpData.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `readLammpData`

```python
readLammpData(file, useCache=True) → dict
```

Read a lammpstrj input data file 



**Arguments:**
 
    - file: the path to the file to read from 
    - useCache: wheter to use cache or not (respects modified date of file) 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
