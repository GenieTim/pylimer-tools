<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/entities/atom.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.entities.atom`






---

<a href="../src/pylimer_tools/entities/atom.py#L7"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Atom`




<a href="../src/pylimer_tools/entities/atom.py#L11"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(data: 'Series', boxSizes: 'list', name: 'str' = None)
```

Instantiate the Atom. 



**Arguments:**
 
  - data: the data underlying the Atom 
  - boxSizes: the size of the box the atom is in. Used for periodic image computations. 




---

<a href="../src/pylimer_tools/entities/atom.py#L100"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `computeDistanceTo`

```python
computeDistanceTo(secondAtom: 'Atom') → float
```

Calculate the the distance between two atoms.  



**Arguments:**
 
    - secondAtom: the atom to compute the distance to 



**Returns:**
 
    - meanDistance: the norm of the connecting vector between the two coordinates 

---

<a href="../src/pylimer_tools/entities/atom.py#L85"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `computeVectorTo`

```python
computeVectorTo(secondAtom: 'Atom') → float
```

Calculate the the vector between two atoms.  



**Arguments:**
 
    - secondAtom: the atom to compute the distance to 



**Returns:**
 
    - difference (np.array): the connecting vector between the two coordinates 

---

<a href="../src/pylimer_tools/entities/atom.py#L46"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getDeltaX`

```python
getDeltaX(secondAtom: 'Atom') → float
```

Calculate the distance in the x dimension between this and another atom,  accounting for periodic displacements. 



**Arguments:**
 
  - secondAtom: the other atom to calculate the distance between 



**Returns:**
 
  - delta: the distance in the x direction 

---

<a href="../src/pylimer_tools/entities/atom.py#L59"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getDeltaY`

```python
getDeltaY(secondAtom: 'Atom') → float
```

Calculate the distance in the y dimension between this and another atom,  accounting for periodic displacements. 



**Arguments:**
 
  - secondAtom: the other atom to calculate the distance between 



**Returns:**
 
  - delta: the distance in the y direction 

---

<a href="../src/pylimer_tools/entities/atom.py#L72"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getDeltaZ`

```python
getDeltaZ(secondAtom: 'Atom') → float
```

Calculate the distance in the z dimension between this and another atom,  accounting for periodic displacements. 



**Arguments:**
 
  - secondAtom: the other atom to calculate the distance between 



**Returns:**
 
  - delta: the distance in the z direction 

---

<a href="../src/pylimer_tools/entities/atom.py#L113"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getUnderlyingData`

```python
getUnderlyingData() → Series
```

Auxilary method to get the pd.Series data associated with this atom 



**Returns:**
 
  - data (pd.Series): the data as given to this atom upon instantiation 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
