<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/entities/dumpFileUniverseSequence.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.entities.dumpFileUniverseSequence`






---

<a href="../src/pylimer_tools/entities/dumpFileUniverseSequence.py#L12"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `DumpFileUniverseSequence`




<a href="../src/pylimer_tools/entities/dumpFileUniverseSequence.py#L14"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(initialStructureFile, dumpFile)
```

Initialize this class. 



**Arguments:**
 
  - initialStructureFile: the file containing the initial structure,   file of type LAMMPS data  
  - dumpFile: the file containing the LAMMPS dump with the atom coordinates   at different time-steps. 




---

<a href="../src/pylimer_tools/entities/dumpFileUniverseSequence.py#L30"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `assembleUniverse`

```python
assembleUniverse(idx) → Universum
```

Create a Universe object for a specified time-step. 



**Arguments:**
 
    - idx: the index of the data 



**Returns:**
 
    - universe (Universum): the universe for that index 

---

<a href="../src/pylimer_tools/entities/dumpFileUniverseSequence.py#L63"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getUniverseAtIndex`

```python
getUniverseAtIndex(idx) → Universum
```

Get the Universe at the given index (as of in the sequence given by the dump file). 



**Arguments:**
 
    - idx: The index to get the universe at 



**Returns:**
 
    - universe (Universum): the Universe at the given index 

---

<a href="../src/pylimer_tools/entities/dumpFileUniverseSequence.py#L78"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `getUniverseAtTimestep`

```python
getUniverseAtTimestep(timestep) → Universum
```

Get the Universe at the given timestep. 



**Arguments:**
 
    - timestep: The timestep to get the universe at 



**Returns:**
 
    - universe (Universum): the Universe at the given timestep,   `None` if there is no entry with this timestep in the dump file. 

---

<a href="../src/pylimer_tools/entities/dumpFileUniverseSequence.py#L100"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `next`

```python
next()
```








---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
