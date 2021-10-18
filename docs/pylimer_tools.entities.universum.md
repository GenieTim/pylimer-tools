<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/entities/universum.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.entities.universum`






---

<a href="../src/pylimer_tools/entities/universum.py#L18"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Universum`




<a href="../src/pylimer_tools/entities/universum.py#L22"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(boxSizes: 'list')
```

Instantiate this Universe (Collection of Molecules) 



**Arguments:**
 
  - boxSizes: a list containing the box lengths (x, y, z) 



**Returns:**
 
  - self (pylimer_tools.entities.Universum): the new Universum object 




---

<a href="../src/pylimer_tools/entities/universum.py#L35"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `addAtomBondData`

```python
addAtomBondData(atomData: 'DataFrame', bondData: 'DataFrame') → Universum
```

Add atoms and bonds to the underlying graph. 



**Arguments:**
 
  - atomData: the dataframe containing atoms with their positions, id, type etc. 
  - bondData: the dataframe containing two columns: one indicating where the bond originates, one where it goes. Direction irrelevant. 



**Returns:**
 
  - self: the Universum object for a fluent interface. 

---

<a href="../src/pylimer_tools/entities/universum.py#L159"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `determineFunctionalityPerType`

```python
determineFunctionalityPerType(typeCounts: 'Counter' = None) → dict
```

Find the maximum functionality of each atom type in the network 



**Arguments:**
 
  - typeCounts: the count of each type in the network. Optional to reduce duplicate counting costs. 



**Returns:**
 
  - functionalitites (dict): a dictionary with key: type, and value: functionality of this atom type.  

---

<a href="../src/pylimer_tools/entities/universum.py#L177"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getAtom`

```python
getAtom(atomId: 'int') → Atom
```

Find an atom by its ID 



**Arguments:**
 
  - atomId: the ID of the atom 



**Returns:**
 
  - atom (pylimer_tools.entities.Atom): the Atom object or None if it is not found 

---

<a href="../src/pylimer_tools/entities/universum.py#L196"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getAtomsWithType`

```python
getAtomsWithType(atomType) → list[Atom]
```

Find an atom by its type 



**Arguments:**
 
  - atomType: the type of the atom 



**Returns:**
 
  - atoms (list<pylimer_tools.entities.Atom>): the Atom objects or None if it is not found 

---

<a href="../src/pylimer_tools/entities/universum.py#L89"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getChainsWithCrosslinker`

```python
getChainsWithCrosslinker(crosslinkerType) → list[Molecule]
```

Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms, without omitting the crosslinkers. In turn, e.g. for a tetrafunctional crosslinker, it will be 4 times in the resulting molecules 



**Arguments:**
 
  - crosslinkerType: the atom type to use to split the molecules 



**Returns:**
 
  - molecules (list): a list of Molecule objects 

---

<a href="../src/pylimer_tools/entities/universum.py#L66"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getMolecules`

```python
getMolecules(ignoreAtomType=None) → list[Molecule]
```

Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms. 



**Arguments:**
 
  - ignoreAtomType: the atom type to ignore/omit from the molecules 



**Returns:**
 
  - molecules (list): a list of Molecule objects 

---

<a href="../src/pylimer_tools/entities/universum.py#L223"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getSize`

```python
getSize()
```

Get the number of atoms in this universe 



**Returns:**
 
  - nr (int): the number of atoms (nodes) 

---

<a href="../src/pylimer_tools/entities/universum.py#L214"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getVolume`

```python
getVolume()
```

Get this object's volume 



**Returns:**
 
  - volume (float): the volume of the box 

---

<a href="../src/pylimer_tools/entities/universum.py#L245"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `reset`

```python
reset() → Universum
```

Reset this Universe to be empty again. 



**Returns:**
 
  - self (pylimer_tools.entities.Universum): the Universum object for a fluent interface. 

---

<a href="../src/pylimer_tools/entities/universum.py#L232"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `setBoxSizes`

```python
setBoxSizes(boxSizes: 'list') → Universum
```

Re-set this Universe's size. 



**Arguments:**
 
  - boxSizes: a list containing the box lengths (x, y, z) 



**Returns:**
 
  - self (pylimer_tools.entities.Universum): the Universum object for a fluent interface. 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
