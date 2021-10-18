<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/calc/calculateRg.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.calc.calculateRg`





---

<a href="../src/pylimer_tools/calc/calculateRg.py#L15"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateRg`

```python
calculateRg(x, y, z)
```

Calculate `<R_g^2>` 



**Arguments:**
 
 - <b>`x`</b> (float):  the x coordinate of the atom 
 - <b>`y`</b> (float):  the y coordinate of the atom 
 - <b>`z`</b> (float):  the z coordinate of the atom 



**Returns:**
 
 - <b>`retVal`</b>:  `<R_g^2>` 


---

<a href="../src/pylimer_tools/calc/calculateRg.py#L35"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateRg2AvgForMolecules`

```python
calculateRg2AvgForMolecules(
    molecules,
    atoms,
    bonds,
    boxLengths: list,
    expected_num_bonds: int = 2
)
```

Calculate `<R_g^2>` 



**Arguments:**
 
 - <b>`molecules`</b> (list):  list of molecules to calculate the R_g for 
 - <b>`atoms`</b> (pd.DataFrame):  a collection of atom coordinates 
 - <b>`bonds`</b> (pd.DataFrame):  a collection of bonds 
 - <b>`boxLengths`</b>:  a list containing the box lengths (x, y, z)  
 - <b>`expected_num_bonds`</b>:  the expected number of bonds. Used for warnings. 



**Returns:**
 
 - <b>`retVal`</b>:  `<R_g^2>` 


---

<a href="../src/pylimer_tools/calc/calculateRg.py#L81"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateRg2Avg`

```python
calculateRg2Avg(
    atoms,
    bonds,
    boxDimensions: list,
    crosslinker_type: int = 2,
    expected_num_bonds: int = 2
) → float
```

Calculate `<R_g^2>` 



**Arguments:**
 
    - atoms (pd.DataFrame): a collection of atom coordinates 
    - bonds (pd.DataFrame): a collection of bonds 
    - boxDimensions: a list containing the box lengths (x, y, z)  
    - crosslinker_type: the type id of the crosslinker (in order to detect chains, filter those) 
    - expected_num_bonds: the expected number of bonds. Used for warnings. 



**Returns:**
 
    - retVal: `<R_g^2>` 


---

<a href="../src/pylimer_tools/calc/calculateRg.py#L102"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateRg2AvgForFile`

```python
calculateRg2AvgForFile(
    file: str,
    fileType: str = 'data',
    useCache: bool = True
) → float
```

Calculate `<Rg^2>` for a specific file 



**Arguments:**
 
    - file: the path to the file to calculate Ree for 
    - fileType: the type of file to read (data or dump) 
    - useCache: wheter to use cache or not 



**Returns:**
 
    - retVal: `<R_g^2>` 


---

<a href="../src/pylimer_tools/calc/calculateRg.py#L138"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateRgForFile`

```python
calculateRgForFile(
    file: str,
    fileType: str = 'data',
    useCache: bool = True
) → <built-in function array>
```

Calculate R_g for a specific file 



**Arguments:**
 
    - file: the path to the file to calculate Ree for 
    - fileType: the type of file to read (data or dump) 
    - useCache: wheter to use cache or not 



**Returns:**
 
    - retVal (list): R_gs 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
