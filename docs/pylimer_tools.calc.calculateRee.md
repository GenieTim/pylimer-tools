<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/calc/calculateRee.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.calc.calculateRee`





---

<a href="../src/pylimer_tools/calc/calculateRee.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateRee2AvgForMolecules`

```python
calculateRee2AvgForMolecules(
    molecules,
    atoms,
    bonds,
    boxLengths: list,
    expected_num_bonds: int = 2
)
```

Calculate `<Ree^2>` 



**Arguments:**
 
    - molecules (List): a list of molecules as produced by pylimer_tools.utils.getMolecules 
    - atoms (pd.DataFrame): a collection of atom coordinates 
    - bonds (pd.DataFrame): a collection of bonds 
    - boxLenghts: a list containing the box lengths (x, y, z)  



**Returns:**
 
    - retVal: `<Ree^2>` 
    - Ree: Ree for each molecule 


---

<a href="../src/pylimer_tools/calc/calculateRee.py#L71"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateRee2Avg`

```python
calculateRee2Avg(
    atoms,
    bonds,
    boxDimensions: list,
    crosslinker_type: int = 2,
    expected_num_bonds: int = 2
) → float
```

Calculate `<Ree^2>` 



**Arguments:**
 
 - <b>`atoms`</b> (pd.DataFrame):  a collection of atom coordinates 
 - <b>`bonds`</b> (pd.DataFrame):  a collection of bonds 
 - <b>`boxDimensions`</b>:  a list containing the box lengths (x, y, z)  
 - <b>`crosslinker_type`</b>:  the type id of the crosslinker (in order to detect chains, filter those) 
 - <b>`expected_num_bonds`</b>:  the expected number of bonds. Used for warnings. 



**Returns:**
 
 - <b>`retVal`</b>:  `<Ree^2>` 


---

<a href="../src/pylimer_tools/calc/calculateRee.py#L92"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateRee2AvgForFile`

```python
calculateRee2AvgForFile(
    file: str,
    fileType: str = 'data',
    useCache: bool = True
) → float
```

Calculate `<Ree^2>` for a specific file 



**Arguments:**
 
 - <b>`file`</b>:  the path to the file to calculate Ree for 
 - <b>`fileType`</b>:  the type of file to read (data or dump) 
 - <b>`useCache`</b>:  wheter to use cache or not 



**Returns:**
 
 - <b>`retVal`</b>:  `<Ree^2>` 


---

<a href="../src/pylimer_tools/calc/calculateRee.py#L128"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateReeForFile`

```python
calculateReeForFile(
    file: str,
    fileType: str = 'data',
    useCache: bool = True
) → <built-in function array>
```

Calculate Ree's for a specific file 



**Arguments:**
 
 - <b>`file`</b>:  the path to the file to calculate Ree for 
 - <b>`fileType`</b>:  the type of file to read (data or dump) 
 - <b>`useCache`</b>:  wheter to use cache or not 



**Returns:**
 
 - <b>`retVal`</b>:  <Ree^2> 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
