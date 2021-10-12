<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/calc/calculateBondLen.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.calc.calculateBondLen`





---

<a href="../src/pylimer_tools/calc/calculateBondLen.py#L11"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateMeanBondLen`

```python
calculateMeanBondLen(coordsDf: DataFrame, boxLengths: list)
```

Calculate the mean bond length  given the coordinates of the atoms in a pd.DataFrame  and the boxLengths in a list 

Assumes the bonds are by the coordsDf's ids, sequentially 

.. deprecated:: 0.0.1   This function is legacy compliant only. Use `calculateBondLen(...).mean()` instead. 



**Arguments:**
 
    - coordsDf: a dataframe containing the coordinates 
    - boxLenghts: a list containing the box lengths (x, y, z)  



**Returns:**
 
    - meanDistance: the mean of all distances 


---

<a href="../src/pylimer_tools/calc/calculateBondLen.py#L60"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateBondLen`

```python
calculateBondLen(
    coordsDf: DataFrame,
    bondsDf: DataFrame,
    boxLengths: list,
    skipAtomType=None
)
```

Calculate the bond lengths given the coordinates of the atoms in a pd.DataFrame  and the bonds of the atoms in a pd.DataFrame  and the boxLengths in a list 



**Arguments:**
 
    - coordsDf: a dataframe containing the coordinates 
    - boxLenghts: a list containing the box lengths (x, y, z)  



**Returns:**
 
    - a np.array of all bond lengths 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
