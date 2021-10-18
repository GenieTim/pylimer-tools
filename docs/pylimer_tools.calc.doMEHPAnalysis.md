<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.calc.doMEHPAnalysis`





---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L15"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `predictShearModulus`

```python
predictShearModulus(
    networks: 'Iterable[Universum]',
    T: 'float' = 1,
    k_B: 'float' = 1,
    foreignAtomType=None,
    totalMass=1
) → float
```

Predict the shear modulus using ANT Analysis. 

Source: 
  - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262 



**Arguments:**
 
  - network: the polymer system to predict the shear modulus for 
  - T: the temperature in your unit system 
  - k_b: Boltzmann's constant in your unit system 
  - foreignAtomType: the type of atoms to ignore (junctions, crosslinkers) 
  - totalMass: the $M$ in the respective formula 



**Returns:**
 
  - shear modulus (float): the estimated shear modulus. Unit: [pressure] 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L39"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateCycleRank`

```python
calculateCycleRank(
    networks: 'Iterable[Universum]' = None,
    nu: 'int' = None,
    mu: 'int' = None,
    absTol: 'float' = 1,
    relTol: 'float' = 1,
    junctionType=None
) → float
```

Compute the cycle rank ($\chi$). Assumes the precursor-chains to be bifunctional. 



**Arguments:**
 
  - network: the network to calculate the cycle rank for 
  - nu: number of elastically effective (active) strands per unit volume 
  - mu: number density of the elastically effective crosslink 
  - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol) 
  - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol)) 
  - junctionType: the atom type of the crosslinkers/junctions 

No need to provide all the parameters — either/or: 
- nu & mu 
- network, absTol, relTol, junctionType 



**Returns:**
 
  - cycleRank: the cycle rank ($\xi = \nu_{eff} - \mu_{eff}$). Unit: [1/Volume] 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L75"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateEffectiveNrDensityOfNetwork`

```python
calculateEffectiveNrDensityOfNetwork(
    networks: 'Iterable[Universum]',
    absTol: 'float' = 1,
    relTol: 'float' = 1,
    junctionType=None
) → float
```

Compute the effective number density $\nu_{eff}$ of a network. Assumes the precursor-chains to be bifunctional. 

$\nu_{eff}$ is the number of elastically effective (active) strands per unit volume, which are defined as the ones that can store elastic energy upon network deformation, resp. the effective number density of network strands 

Source: 
  - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262 



**Arguments:**
 
  - network (pylimer_tools.entities.Universum): the network to compute $\nu_{eff}$ for 
  - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol) 
  - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol)) 
  - junctionType: the atom type of the crosslinkers/junctions 



**Returns:**
 
  - $\nu_{eff}$ (float): the effective number density of network strands. Unit: [1/Volume] 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L118"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateMeanUniverseVolume`

```python
calculateMeanUniverseVolume(
    networks: 'Iterable[Universum]',
    acceptDifferentSizes: 'bool' = False
) → float
```

Compute the mean volume of a list of universes. 



**Arguments:**
 
  - networks: a list of universes 
  - acceptDifferentSizes: toggle whether to throw an error when the Universe have different nr. of atoms 



**Returns:**
 
  - meanVolume (float): the mean volume of the universes 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L141"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateEffectiveNrDensityOfJunctions`

```python
calculateEffectiveNrDensityOfJunctions(
    networks: 'Iterable[Universum]',
    absTol: 'float' = 0,
    relTol: 'float' = 1,
    junctionType=None,
    minNumEffectiveStrands=2
) → float
```

Compute the number density of the elastically effective crosslinks, defined as the ones that connect at least two elastically effective strands. Assumes the precursor-chains to be bifunctional. 

Source: 
  - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262 



**Arguments:**
 
  - network (pylimer_tools.entities.Universum): the network to compute $\nu_{eff}$ for 
  - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol) 
  - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol)) 
  - junctionType: the atom type of the crosslinkers/junctions 
  - minNumEffectiveStrands (int): the number of elastically effective strands to qualify a junction as such 



**Returns:**
 
  - $\mu_{eff}$ (float): the effective number density of junctions. Unit: [1/Volume] 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L200"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateWeightFractionOfBackbone`

```python
calculateWeightFractionOfBackbone(
    network: 'Universum',
    crosslinkerType,
    weights=1
)
```

Compute the weight fraction of network backbone in infinite network 



**Arguments:**
 
  - network: the network to compute the weight fraction for 
  - crosslinkerType: the atom type to use to split the molecules 
  - weights: either a dict with key: atomType and value: weight, or a scalar value if all atoms have the same weight 



**Returns:**
 
  - weightFraction (float): 1 - weightDangling/weightTotal, 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L217"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateWeightFractionOfDanglingChains`

```python
calculateWeightFractionOfDanglingChains(
    network: 'Universum',
    crosslinkerType,
    weights=1
) → Tuple[float, float]
```

Compute the weight fraction of dangling strands in infinite network 



**Arguments:**
 
  - network: the network to compute the weight fraction for 
  - crosslinkerType: the atom type to use to split the molecules 
  - weights: either a dict with key: atomType and value: weight, or a scalar value if all atoms have the same weight 



**Returns:**
 
  - weightFraction: weightDangling/weightTotal, 
  - numFraction: numDangling/numTotal 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L261"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeMeanEndToEndDistances`

```python
computeMeanEndToEndDistances(
    networks: 'Iterable[Universum]',
    crosslinkerType
) → dict
```

Compute the mean end to end distance between each pair of (indirectly) connected crosslinker 



**Arguments:**
 
  - networks: the different configurations of the polymer network to do the computation for 
  - crosslinkerType: the atom type to compute the in-between vectors for 



**Returns:**
 
  - endToEndDistances (dict): a dictionary with key: "{atom1.name}+{atom2.name}" 
 - <b>`and value`</b>:  the norm of the mean difference vector 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L284"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeMeanEndToEndVectors`

```python
computeMeanEndToEndVectors(
    networks: 'Iterable[Universum]',
    crosslinkerType
) → dict
```

Compute the mean end to end vectors between each pair of (indirectly) connected crosslinker 



**Arguments:**
 
  - networks: the different configurations of the polymer network to do the computation for 
  - crosslinkerType: the atom type to compute the in-between vectors for 



**Returns:**
 
  - endToEndVectors (dict): a dictionary with key: "{atom1.name}+{atom2.name}" 
 - <b>`and value`</b>:  their mean difference vector 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L315"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeEndToEndVectors`

```python
computeEndToEndVectors(network: 'Universum', crosslinkerType) → dict
```

Compute the end to end vectors between each pair of (indirectly) connected crosslinker 



**Arguments:**
 
  - network: the polymer network to do the computation for 
  - crosslinkerType: the atom type to compute the in-between vectors for 



**Returns:**
 
  - endToEndVectors (dict): a dictionary with key: "{atom1.name}+{atom2.name}" 
 - <b>`and value`</b>:  their difference vector 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L350"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeCrosslinkerConversion`

```python
computeCrosslinkerConversion(
    network: 'Universum',
    junctionType,
    f: 'int'
) → float
```

Compute the extent of reaction of the crosslinkers (actual functionality divided by target functionality) 



**Arguments:**
 
  - network: the polymer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 
  - f: the functionality of the crosslinkers 



**Returns:**
 
  - r (float): the (mean) crosslinker conversion 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L366"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateEffectiveCrosslinkerFunctionality`

```python
calculateEffectiveCrosslinkerFunctionality(
    network: 'Universum',
    junctionType
) → float
```

Compute the mean crosslinker functionality 



**Arguments:**
 
  - network: the polymer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 



**Returns:**
 
  - f (float): the (mean) effective crosslinker functionality 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L382"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateEffectiveCrosslinkerFunctionalities`

```python
calculateEffectiveCrosslinkerFunctionalities(
    network: 'Universum',
    junctionType
) → list[int]
```

Compute the functionality of every crosslinker in the network 



**Arguments:**
 
  - network: the polymer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 



**Returns:**
 
  - junctionDegrees (list[int]): the functionality of every crosslinker 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L402"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateTopologicalFactor`

```python
calculateTopologicalFactor(
    networks: 'Iterable[Universum]',
    foreignAtomType=None,
    totalMass=1,
    b=None
) → float
```

Compute the topological factor of a polymer network. 

Assumptions:  
  - the precursor-chains to be bifunctional 
  - all Universes to have the same structure (with possibly differing positions) 
  - crosslinkers do not count to the nr. of monomers in a strand 

Source: 
  - eq. 16 in https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262 



**Arguments:**
 
  - network: the network to compute the topological factor for 
  - foreignAtomType: the type of atoms to ignore 
  - totalMass: the $M$ in the respective formula 
  - b: the mean bond length.   If `None`, it will be computed for each molecule in the first Universum (Network). 



**Returns:**
 
  - the topological factor $\Gamma$ 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
