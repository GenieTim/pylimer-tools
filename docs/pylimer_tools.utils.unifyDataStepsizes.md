<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/utils/unifyDataStepsizes.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.utils.unifyDataStepsizes`





---

<a href="../src/pylimer_tools/utils/unifyDataStepsizes.py#L5"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `unifyDataStepsizes`

```python
unifyDataStepsizes(
    data: DataFrame,
    key: str,
    stepSize: int = None,
    maxExpectedStepSize: int = 100
)
```

Get a DataFrame where all data points have the same step between the values in column given by `key` NOTE: this function is rather unstable, as it has a few dirty assumptions, such as: 
- steps are modulo stepsize. Breaks e.g. with steps start with 1 and go up by stepSize. 
- ideal step-size is max step difference. Breaks e.g. if there is one big gap 



**Arguments:**
 
    - data: the DataFrame to unify the step-size for 
    - key: the column name indicating the column containing the step-nr. 
    - maxExpectedStepSize: use to get a warning if the computed step-size is larger 



**Returns:**
 
    - data: a DataFrame with a consistent step-size 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
