<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.calc.doMEHPAnalysis`





---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateCycleRank`

```python
calculateCycleRank(
    network: Universum,
    nu: int = None,
    mu: int = None,
    absTol: float = 1,
    relTol: float = 1,
    junctionType=None
)
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
 
  - cycleRank: the cycle rank ($\chi = \nu_{eff} - \mu_{eff}$)  


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L40"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateEffectiveNrDensityOfNetwork`

```python
calculateEffectiveNrDensityOfNetwork(
    network: Universum,
    absTol: float = 1,
    relTol: float = 1,
    junctionType=None
)
```

Compute the effective number density $\nu_{eff}$ of a network. Assumes the precursor-chains to be bifunctional. 

$\nu_{eff}$ is the number of elastically effective (active) strands per unit volume,  which are defined as the ones that can store elastic energy  upon network deformation, resp. the effective number density of network strands 



**Arguments:**
 
  - network (pylimer_tools.entities.Universum): the network to compute $\nu_{eff}$ for 
  - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol) 
  - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol)) 
  - junctionType: the atom type of the crosslinkers/junctions 



**Returns:**
 
  - $\nu_{eff}$ (float): the effective number density of network strands 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L76"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateEffectiveNrDensityOfJunctions`

```python
calculateEffectiveNrDensityOfJunctions(
    network: Universum,
    absTol: float = 0,
    relTol: float = 0,
    junctionType=None,
    minNumEffectiveStrands=2
)
```

Compute the number density of the elastically effective crosslinks,  defined as the ones that connect at least two elastically effective strands. Assumes the precursor-chains to be bifunctional. 



**Arguments:**
 
  - network (pylimer_tools.entities.Universum): the network to compute $\nu_{eff}$ for 
  - absTol (float): the absolute tolerance to categorize a chain as active (min. end-to-end distance) (None to use only relTol) 
  - relTol (float): the relative tolerance to categorize a chain as active (0: all, 1: none (use only absTol)) 
  - junctionType: the atom type of the crosslinkers/junctions 
  - minNumEffectiveStrands (int): the number of elastically effective strands to qualify a junction as such 



**Returns:**
 
  - $\mu_{eff}$ (float): the effective number density of junctions 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L139"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateWeightFractionOfBackbone`

```python
calculateWeightFractionOfBackbone(network: Universum, crosslinkerType)
```

Compute the weight fraction of network backbone in infinite network 



**Arguments:**
 
  - network: the network to compute the weight fraction for 
  - crosslinkerType: the atom type to use to split the molecules 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L150"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateWeightFractionOfDanglingChains`

```python
calculateWeightFractionOfDanglingChains(
    network: Universum,
    crosslinkerType,
    weights=1
)
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

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L188"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `predictGelationPoint`

```python
predictGelationPoint(r: float, f: int, g: int = 2) → float
```

Compute the gelation point $p_{gel}$ as theoretically predicted (gelation point = critical extent of reaction for gelation) 

Source: 
  - https://www.sciencedirect.com/science/article/pii/003238618990253X 



**Arguments:**
 
  - r: the stoichiometric inbalance of reactants 
  - f: functionality of the crosslinkers 
  - g: functionality of the precursor polymer 



**Returns:**
 
  - p_gel: critical extent of reaction for gelation 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L209"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeCrosslinkerConversion`

```python
computeCrosslinkerConversion(network: Universum, junctionType, f: int) → float
```

Compute the extent of reaction of the crosslinkers  (actual functionality divided by target functionality) 



**Arguments:**
 
  - network: the poylmer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 
  - f: the functionality of the crosslinkers 



**Returns:**
 
  - r (float): the (mean) crosslinker conversion 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L225"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateEffectiveCrosslinkerFunctionality`

```python
calculateEffectiveCrosslinkerFunctionality(
    network: Universum,
    junctionType
) → float
```

Compute the mean crosslinker functionality 



**Arguments:**
 
  - network: the poylmer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 



**Returns:**
 
  - f (float): the (mean) effective crosslinker functionality 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L241"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateEffectiveCrosslinkerFunctionalities`

```python
calculateEffectiveCrosslinkerFunctionalities(network: Universum, junctionType)
```

Compute the functionality of every crosslinker in the network 



**Arguments:**
 
  - network: the poylmer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 



**Returns:**
 
  - junctionDegrees (list[int]): the functionality of every crosslinker 


---

<a href="../src/pylimer_tools/calc/doMEHPAnalysis.py#L258"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateTopologicalFactor`

```python
calculateTopologicalFactor(
    network: Universum,
    foreignAtomType=None,
    totalMass=1
)
```

Compute the topological factor of a polymer network. Assumes the precursor-chains to be bifunctional. 

Source: 
  - eq. 16 in https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262 



**Arguments:**
 
  - network: the network to compute the topological factor for 
  - foreignAtomType: the type of atoms to ignore 
  - totalMass: the $M$ in the respective formula 



**Returns:**
 
  - the topological factor $\Gamma$ 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
