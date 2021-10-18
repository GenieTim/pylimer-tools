<!-- markdownlint-disable -->

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pylimer_tools.calc.doMMTAnalysis`





---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `predictShearModulus`

```python
predictShearModulus(
    network: Universum,
    junctionType,
    weightPerType,
    strandLength: int = None,
    functionalityPerType=None,
    T: float = 1,
    k_B: float = 1,
    totalMass=1
)
```

Predict the shear modulus using MMT Analysis. 

Source: 
  - https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262 



**Arguments:**
 
  - network: the polymer system to predict the shear modulus for 
  - junctionType: the type of atoms making up the junctions/crosslinkers 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 
  - weightPerType: a dictionary with key: type, and value: weight per atom of this atom type. See: #computeWeightFractions 
  - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance 
  - T: the temperature in your unit system 
  - k_b: Boltzmann's constant in your unit system 
  - totalMass: the $M$ in the respective formula 



**Returns:**
 
  - G: the predicted shear modulus, or `None` if the universe is empty. 



**ToDo:**
 
  - Support more than one crosslinker type (as is supported by original formula) 


---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L48"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateWeightFractionOfDanglingChains`

```python
calculateWeightFractionOfDanglingChains(
    network: Universum,
    junctionType,
    weightPerType,
    strandLength: int = None,
    functionalityPerType: dict = None
) → float
```

Compute the weight fraction of dangling strands in infinite network 



**Arguments:**
 
  - network: the network to compute the weight fraction for 
  - crosslinkerType: the atom type to use to split the molecules 
  - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance 
  - weights: either a dict with key: atomType and value: weight, or a scalar value if all atoms have the same weight 



**Returns:**
 
  - weightFraction $\Phi_d = 1 - \Phi_{el}$: weightDangling/weightTotal 


---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L64"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `calculateWeightFractionOfBackbone`

```python
calculateWeightFractionOfBackbone(
    network: Universum,
    junctionType,
    weightPerType,
    strandLength: int = None,
    functionalityPerType: dict = None
) → float
```

Compute the weight fraction of the backbone strands in an infinite network 

Source: 
  - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 (see supporting information for formulae) 



**Arguments:**
 
  - network: the poylmer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 
  - weightPerType: a dictionary with key: type, and value: weight per atom of this atom type. See: #computeWeightFractions 
  - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance 
  - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type.  
 - <b>`See`</b>:  #computeExtentOfReaction 



**Returns:**
 
  - $\Phi_{el}$: weight fraction of network backbone 


---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L107"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeWeightFractionOfSolubleMaterial`

```python
computeWeightFractionOfSolubleMaterial(
    network: Universum,
    junctionType,
    weightPerType,
    strandLength: int = None,
    functionalityPerType: dict = None
) → float
```

Compute the weight fraction of soluble material. 

Source: 
  - https://pubs.acs.org/doi/10.1021/ma00046a021 
  - https://pubs.acs.org/doi/suppl/10.1021/acs.macromol.0c02737 



**Arguments:**
 
  - network: the poylmer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 
  - weightPerType: a dictionary with key: type, and value: weight per atom of this atom type. See: #computeWeightFractions 
  - strandLength: the length of the network strands (in nr. of beads). See: #computeStoichiometricInbalance 
  - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type.  
 - <b>`See`</b>:  #computeExtentOfReaction 



**Returns:**
 
  - $W_{sol}$ (float): the weight fraction of soluble material according to MMT. 
  - weightFractions (dict): a dictionary with key: type, and value: weight fraction of type 
  - $\alpha$ (float): Macosko & Miller's $P(F_A)$ 
  - $\beta$ (float): Macosko & Miller's $P(F_B)$ 


---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L156"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeMMsProbabilities`

```python
computeMMsProbabilities(r, p, f)
```

Compute Macosko and Miller's probabilities $P(F_A)$ and $P(F_B)$ 



**Arguments:**
 
  - r: the stoichiometric inbalance 
  - p: the extent of reaction 
  - f: the functionality of the the crosslinker 



**Returns:**
 
  - alpha: $P(F_A)$ 
  - beta: $P(F_B)$     


---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L192"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeWeightFractions`

```python
computeWeightFractions(network: Universum, weightPerType) → dict
```

Compute the weight fractions of each atom type in the network. 



**Arguments:**
 
  - network: the poylmer network to do the computation for 
  - weightPerType: a dictionary with key: type, and value: weight per atom of this atom type.  



**Returns:**
 
  - $\vec{W_i}$ (dict): using the type i as a key, this dict contains the weight fractions ($\frac{W_i}{W_{tot}}$) 


---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L224"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeStoichiometricInbalance`

```python
computeStoichiometricInbalance(
    network: Universum,
    junctionType,
    strandLength: int = None,
    functionalityPerType: dict = None
) → float
```

Compute the stoichiometric inbalance ( nr. of bonds formable of crosslinker / nr. of formable bonds of precursor ) 



**NOTE:**

> - if your system has a non-integer number of possible bonds (e.g. one site unbonded), this will not be rounded/respected in any way. 
>

**Arguments:**
 
  - network: the poylmer network to do the computation for 
  - junctionType: the type of the junctions/crosslinkers to select them in the network 
  - strandLength: the length of the network strands (in nr. of beads).   Used to infer the number of precursor strands.  
 - <b>`If `None``</b>:  will use average length of each connected system when ignoring the crosslinkers. 
  - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type.  
 - <b>`If `None``</b>:  will use max functionality per type. 



**Returns:**
 
  - r (float): the stoichiometric inbalance 


---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L269"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `computeExtentOfReaction`

```python
computeExtentOfReaction(
    network: Universum,
    functionalityPerType: dict = None
) → float
```

Compute the extent of reaction (nr. of formed bonds in reaction / max. nr. of bonds formable) NOTE: if your system has a non-integer number of possible bonds (e.g. one site unbonded), this will not be rounded/respected in any way.  



**Arguments:**
 
  - network: the poylmer network to do the computation for 
  - functionalityPerType: a dictionary with key: type, and value: functionality of this atom type.  
 - <b>`If None`</b>:  will use max functionality per type. 



**Returns:**
 
  - p (float): the extent of reaction 


---

<a href="../src/pylimer_tools/calc/doMMTAnalysis.py#L305"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `predictGelationPoint`

```python
predictGelationPoint(r: float, f: int, g: int = 2) → float
```

Compute the gelation point $p_{gel}$ as theoretically predicted (gelation point = critical extent of reaction for gelation) 

Source: 
  - https://www.sciencedirect.com/science/article/pii/003238618990253X 



**Arguments:**
 
  - r (double): the stoichiometric inbalance of reactants (see: #computeStoichiometricInbalance) 
  - f (int): functionality of the crosslinkers 
  - g (int): functionality of the precursor polymer 



**Returns:**
 
  - p_gel: critical extent of reaction for gelation 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
