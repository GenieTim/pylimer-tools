<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/utils/getTail.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.utils.getTail`





---

<a href="../src/pylimer_tools/utils/getTail.py#L4"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getTail`

```python
getTail(data, percentage=0.2, minN=25, maxPercentage=0.5)
```

Extract the last few entries of a list 



**Arguments:**
 
    - data (list|pd.DataFrame|pd.Series): the list to extract the last few entries from 
    - percentage: the percentage of entries to extract 
    - minN: the minimum number of entries to extract 
    - maxPercentage: the maximum percentage of entries to extract 



**Returns:**
  A (list|pd.DataFrame|pd.Series) with at maximum maxPercentage,   at least minN entries (assuming the initial data is as large),   but ideally `percentage` many percentage of the last entries. 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
