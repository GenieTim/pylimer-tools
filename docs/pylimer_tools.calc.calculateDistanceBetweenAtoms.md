<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/calc/calculateDistanceBetweenAtoms.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.calc.calculateDistanceBetweenAtoms`





---

<a href="../src/pylimer_tools/calc/calculateDistanceBetweenAtoms.py#L8"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `normalizeData`

```python
normalizeData(atom)
```

Verify the atom to be a "usable" one 



**Arguments:**
 
    - atom: the atom to check 



**Returns:**
 
    - atom (pd.Seroes): the resolved atom 


---

<a href="../src/pylimer_tools/calc/calculateDistanceBetweenAtoms.py#L26"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `row_converter`

```python
row_converter(row: DataFrame, listy=None)
```

convert pandas row to a dictionary 



**Arguments:**
 
    - row (pd.DataFrame|pd.Series): row as a tuple 
    - listy: a list of columns 



**Returns:**
 
    - pictionary (dictionary): the row's values 


---

<a href="../src/pylimer_tools/calc/calculateDistanceBetweenAtoms.py#L48"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateDistanceBetweenAtoms`

```python
calculateDistanceBetweenAtoms(atomA, atomB)
```

Calculate the the distance between two atoms.  No translation between periodic images happens. 



**Arguments:**
 
    - atomA: the coordinates of atom 1 
    - atomB: the coordinates of atom 2 



**Returns:**
 
    - meanDistance: the norm of the connecting vector between the two coordinates 


---

<a href="../src/pylimer_tools/calc/calculateDistanceBetweenAtoms.py#L65"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateNormalizedDistanceBetweenAtoms`

```python
calculateNormalizedDistanceBetweenAtoms(atom1, atom2, boxLengths: list)
```

Calculate the the distance between two atoms.  



**Arguments:**
 
    - atom1: the coordinates of atom 1 
    - atom2: the coordinates of atom 2 
    - boxLenghts: a list containing the box lengths (x, y, z)  



**Returns:**
 
    - meanDistance: the norm of the connecting vector between the two coordinates 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
