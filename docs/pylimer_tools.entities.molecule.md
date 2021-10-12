<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/entities/molecule.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.entities.molecule`






---

<a href="../src/pylimer_tools/entities/molecule.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Molecule`




<a href="../src/pylimer_tools/entities/molecule.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(molecule_graph, chainType=<MoleculeType.UNDEFINED: 0>)
```








---

<a href="../src/pylimer_tools/entities/molecule.py#L67"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `computeBondLengths`

```python
computeBondLengths()
```

Calculate the bond lengths 



**Returns:**
 
  - a np.array of all bond lengths 

---

<a href="../src/pylimer_tools/entities/molecule.py#L48"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `computeEndToEndDistance`

```python
computeEndToEndDistance()
```

Compute the end-to-end distance of this molecule/chain. 



**Returns:**
 
  - $R_ee$ (float): the end-to-end distance,   `None` if the molecule/chain does not have two distinct ends. 

---

<a href="../src/pylimer_tools/entities/molecule.py#L25"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `decomposeFurther`

```python
decomposeFurther(splitAtomType)
```

Split this molecule into smaller molecules by ignoring all atoms with a given type. 



**Arguments:**
 
  - splitAtomType: the type of the atom to omit 



**Returns:**
 
  - subMolecules (list): a list of Molecule objects         

---

<a href="../src/pylimer_tools/entities/molecule.py#L83"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getLength`

```python
getLength()
```

Query the length of the molecule/chain (the nr. of atoms) 



**Returns:**
 
  - len (int): the nr. of nodes (atoms) in this molecule 

---

<a href="../src/pylimer_tools/entities/molecule.py#L92"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getType`

```python
getType() → <enum 'MoleculeType'>
```

Query the type of the molecule/chain (loop, dangling, etc.) 



**Returns:**
 
  - moleculeType (MoleculeType): one of the enum's Molecule.MoleculeType... 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
