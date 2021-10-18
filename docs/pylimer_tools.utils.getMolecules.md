<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/utils/getMolecules.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.utils.getMolecules`





---

<a href="../src/pylimer_tools/utils/getMolecules.py#L4"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getAndFilterMoleculesAndBonds`

```python
getAndFilterMoleculesAndBonds(
    atom_data: DataFrame,
    bond_data: DataFrame,
    boxDimensions: list,
    crosslinker_type: int = None,
    recalculate_positions=True
)
```

Find molecules from atom & bond data. Crosslinkers are omitted. 



**Arguments:**
 
    - atom_data (pd.DataFrame): a collection of atom coordinates 
    - bond_data (pd.DataFrame): a collection of bonds 
    - boxDimensions (list): a list containing the box lengths (x, y, z)  
    - crosslinker_type (int): the type of the crosslinker atoms to distinguish different molecules/chains 



**Returns:**
 
    - molecules (list): a list of lists of atoms representing molecules/chains 
    - atom_data (pd.DataFrame): the atoms, excluding crosslinkers 
    - bond_data: the bonds, excluding the ones of the crosslinker atoms 
    - crosslinker_atoms: the atoms of the crosslinkers 
    - crosslinker_bonds: the bonds of the crosslinker atoms 


---

<a href="../src/pylimer_tools/utils/getMolecules.py#L93"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getFilteredMoleculesAndBonds`

```python
getFilteredMoleculesAndBonds(
    atom_data: DataFrame,
    bond_data: DataFrame,
    boxDimensions: list,
    crosslinker_type: int = None
)
```

Find molecules from atom & bond data. Crosslinkers are omitted. 



**Arguments:**
 
    - atom_data (pd.DataFrame): a collection of atom coordinates 
    - bond_data (pd.DataFrame): a collection of bonds 
    - boxDimensions: a list containing the box lengths (x, y, z)  
    - crosslinker_type: the type of the crosslinker atoms to distinguish different molecules/chains 



**Returns:**
 
    - molecules (list): a list of lists of atoms representing molecules/chains 
    - atom_data (pd.DataFrame): the atoms, excluding crosslinkers 
    - bond_data: the bonds, excluding the ones of the crosslinker atoms 


---

<a href="../src/pylimer_tools/utils/getMolecules.py#L113"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `getMolecules`

```python
getMolecules(
    atom_data: DataFrame,
    bond_data: DataFrame,
    boxDimensions: list,
    crosslinker_type: int = None
)
```

Find molecules from atom & bond data. Crosslinkers are omitted. 



**Arguments:**
 
    - atom_data (pd.DataFrame): a collection of atom coordinates 
    - bond_data (pd.DataFrame): a collection of bonds 
    - boxDimensions: a list containing the box lengths (x, y, z)  
    - crosslinker_type: the type of the crosslinker atoms to distinguish different molecules/chains 



**Returns:**
 
    - molecules (list): a list of lists of atoms representing molecules/chains 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
