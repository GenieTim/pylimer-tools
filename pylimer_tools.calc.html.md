# pylimer_tools.calc package

## Submodules

## pylimer_tools.calc.mehp_utilities module

This module is deprecated.
Use [`pylimer_tools_cpp.MEHPForceBalance()`](api/pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance) or [`pylimer_tools_cpp.MEHPForceRelaxation()`](api/pylimer_tools_cpp.MEHPForceRelaxation.md#pylimer_tools_cpp.MEHPForceRelaxation) instead.

In principle, it offers comparable functionality, but without the force minimization - the networks you pass
to these methods should be minimized first.

### pylimer_tools.calc.mehp_utilities.compute_cycle_rank(networks: Sequence[[Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)] | [None](https://docs.python.org/3/library/constants.html#None) = None, nu: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, mu: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, abs_tol: [float](https://docs.python.org/3/library/functions.html#float) = 1, rel_tol: [float](https://docs.python.org/3/library/functions.html#float) = 1, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the cycle rank ($\chi$).
Assumes the precursor-chains to be bifunctional.

No need to provide all the parameters — either/or:

* nu & mu
* network, abs_tol, rel_tol, crosslinker_type

* **Parameters:**
  * **networks** (*Union* *[**Sequence* *[*[*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) *]* *,* *None* *]*) – The networks to calculate the cycle rank for
  * **nu** (*Union* *[*[*float*](https://docs.python.org/3/library/functions.html#float) *,* *None* *]*) – Number of elastically effective (active) strands per unit volume
  * **mu** (*Union* *[*[*float*](https://docs.python.org/3/library/functions.html#float) *,* *None* *]*) – Number density of the elastically effective crosslinks
  * **abs_tol** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The absolute tolerance to categorize a chain as active (min. end-to-end distance).
    Set to None to use only rel_tol
  * **rel_tol** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The relative tolerance to categorize a chain as active (0: all, 1: none)
  * **crosslinker_type** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The atom type of the crosslinkers/junctions
* **Returns:**
  The cycle rank ($\xi = \nu_{eff} - \mu_{eff}$) in units of 1/Volume
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.mehp_utilities.compute_effective_nr_density_of_junctions(networks: Sequence[[Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)], abs_tol: [float](https://docs.python.org/3/library/functions.html#float) = 0, rel_tol: [float](https://docs.python.org/3/library/functions.html#float) = 1, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, min_num_effective_strands=2) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the number density of the elastically effective crosslinks,
defined as the ones that connect at least min_num_effective_strands elastically effective strands.
Assumes the precursor-chains to be bifunctional.

Source:
: - [https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262](https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262)

* **Parameters:**
  * **networks** (*Sequence* *[*[*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) *]*) – The networks to compute $\mu_{eff}$ for
  * **abs_tol** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The absolute tolerance to categorize a chain as active (min. end-to-end distance).
    Set to None to use only rel_tol
  * **rel_tol** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The relative tolerance to categorize a chain as active (0: all, 1: none)
  * **crosslinker_type** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The atom type of the crosslinkers/junctions
  * **min_num_effective_strands** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The number of elastically effective strands to qualify a junction as such
* **Returns:**
  The effective number density of junctions in units of 1/Volume
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.mehp_utilities.compute_effective_nr_density_of_network(networks: Sequence[[Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)], abs_tol: [float](https://docs.python.org/3/library/functions.html#float) = 1, rel_tol: [float](https://docs.python.org/3/library/functions.html#float) = 1, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the effective number density $\nu_{eff}$ of a network.
Assumes the precursor-chains to be bifunctional.

$\nu_{eff}$ is the number of elastically effective (active) strands per unit volume,
which are defined as the ones that can store elastic energy
upon network deformation, resp. the effective number density of network strands.

Source:
: - [https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262](https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262)

* **Parameters:**
  * **networks** (*Sequence* *[*[*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) *]*) – The networks to compute $\nu_{eff}$ for
  * **abs_tol** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The absolute tolerance to categorize a chain as active (min. end-to-end distance).
    Set to None to use only rel_tol
  * **rel_tol** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The relative tolerance to categorize a chain as active (0: all, 1: none)
  * **crosslinker_type** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The atom type of the crosslinkers/junctions
* **Returns:**
  The effective number density of network strands in units of 1/Volume
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.mehp_utilities.compute_mean_universe_volume(networks: Sequence[[Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)], accept_different_sizes: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the mean volume of a list of universes.

* **Parameters:**
  * **networks** (*Sequence* *[*[*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) *]*) – A list of universes
  * **accept_different_sizes** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Toggle whether to throw an error when the Universes have different numbers of atoms
* **Returns:**
  The mean volume of the universes
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)
* **Raises:**
  * [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If the input list is empty
  * [**NotImplementedError**](https://docs.python.org/3/library/exceptions.html#NotImplementedError) – If universes have different sizes and accept_different_sizes is False

### pylimer_tools.calc.mehp_utilities.compute_topological_factor(networks: Sequence[[Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)], crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, total_mass: [float](https://docs.python.org/3/library/functions.html#float) = 1, b: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the topological factor of a polymer network.

Assumptions:
: - the precursor-chains to be bifunctional
  - all Universes to have the same structure (with possibly differing positions)
  - crosslinkers do not count to the nr. of monomers in a strand

Source:
: - eq. 16 in [https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262](https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262)

* **Parameters:**
  * **networks** (*Sequence* *[*[*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) *]*) – The networks to compute the topological factor for
  * **crosslinker_type** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The type of atoms to ignore (junctions, crosslinkers)
  * **total_mass** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The $M$ in the respective formula
  * **b** ([*float*](https://docs.python.org/3/library/functions.html#float) *or* *None*) – The mean bond length. If None, it will be computed for each molecule in the first Universe
* **Returns:**
  The topological factor $\Gamma$
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.mehp_utilities.predict_shear_modulus(networks: Sequence[[Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)], temperature: [float](https://docs.python.org/3/library/functions.html#float) = 1, k_boltzmann: [float](https://docs.python.org/3/library/functions.html#float) = 1, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, total_mass=1) → [float](https://docs.python.org/3/library/functions.html#float)

Predict the shear modulus using ANT Analysis.

Source:
: - [https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262](https://pubs.acs.org/doi/10.1021/acs.macromol.9b00262)

* **Parameters:**
  * **networks** (*Sequence* *[*[*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) *]*) – The polymer systems to predict the shear modulus for
  * **temperature** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The temperature in your unit system
  * **k_boltzmann** ([*float*](https://docs.python.org/3/library/functions.html#float)) – Boltzmann’s constant in your unit system
  * **crosslinker_type** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The type of atoms to ignore (junctions, crosslinkers)
  * **total_mass** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The $M$ in the respective formula
* **Returns:**
  The estimated shear modulus in pressure units
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

## pylimer_tools.calc.miller_macosko_theory module

This module provides access to various computations introduced in the Miller-Macosko theory.
See Miller and Macosko [[MM76b](acknowledgements.md#id6)] and Macosko and Miller [[MM76a](acknowledgements.md#id5)].

Additional references used in the development of this module include
Langley [[Lan68](acknowledgements.md#id22)], Miller and Macosko [[MM78](acknowledgements.md#id15)],
Valles and Macosko [[VM79](acknowledgements.md#id17)], Miller *et al.* [[MVM79](acknowledgements.md#id16)],
Venkataraman *et al.* [[VCC+89](acknowledgements.md#id19)], Patel *et al.* [[PMC+92](acknowledgements.md#id14)],
Aoyama *et al.* [[AYU21](acknowledgements.md#id13)], Gusev and Schwarz [[GS22](acknowledgements.md#id21)],
and Tsimouri *et al.* [[TSBG24](acknowledgements.md#id23)].

### pylimer_tools.calc.miller_macosko_theory.compute_entanglement_modulus(g_e_1: [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity), temperature: [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity), network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) | [None](https://docs.python.org/3/library/constants.html#None) = None, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, p: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, r: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, f: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, b2: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, beta: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity)

Compute MMT’s entanglement contribution to the equilibrium shear modulus, given by
$k_B T \epsilon_e T_e$.

* **Parameters:**
  * **g_e_1** – The melt entanglement modulus $G_e(1) = k_B T \epsilon_e$
  * **temperature** – The temperatures; defaults to room temperature (25 °C)
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
  * **p** – The crosslinker conversion
  * **r** – The stoichiometric imbalance
  * **f** – The functionality of the crosslinkers
  * **beta** – $P(F_b^{out})$, optional
* **Returns:**
  $T_e \cdot G_e(1)$

### pylimer_tools.calc.miller_macosko_theory.compute_extracted_modulus(p: [float](https://docs.python.org/3/library/functions.html#float), r: [float](https://docs.python.org/3/library/functions.html#float), f: [int](https://docs.python.org/3/library/functions.html#int), g_e_1: [Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity), w_sol: [float](https://docs.python.org/3/library/functions.html#float), xlink_concentration_0: [Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity), ureg: [UnitRegistry](https://pint.readthedocs.io/en/stable/api/base.html#pint.UnitRegistry), temperature: [Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity) | [None](https://docs.python.org/3/library/constants.html#None) = None, b2: [float](https://docs.python.org/3/library/functions.html#float) = 1.0) → [Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity)

Compute MMT’s modulus, assuming the solvent is removed

* **Parameters:**
  * **p** – The crosslinker conversion
  * **r** – The stoichiometric imbalance
  * **f** – The functionality of the crosslinkers
  * **g_e_1** – The melt entanglement modulus $G_e(1) = k_B T \epsilon_e$
  * **xlink_concentration_0** – [A_f]_0, in 1/volume units
  * **ureg** – The unit registry to use
  * **alpha** – $P(F_a^{out})$, optional
  * **temperature** – The temperatures; defaults to room temperature
  * **w_sol** – The soluble fraction (to be removed)

### pylimer_tools.calc.miller_macosko_theory.compute_junction_modulus(p: [float](https://docs.python.org/3/library/functions.html#float), r: [float](https://docs.python.org/3/library/functions.html#float), xlink_concentration_0: [Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity), ureg: [UnitRegistry](https://pint.readthedocs.io/en/stable/api/base.html#pint.UnitRegistry), f: [int](https://docs.python.org/3/library/functions.html#int), temperature: [Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity) | [None](https://docs.python.org/3/library/constants.html#None) = None, b2: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity)

Compute MMT’s junction modulus, given by
$G_{junctions} = k_B T [A_f]_0 \sum_{m=3}^{f} \frac{m-2}{2} P(X_{m,f})$.

* **Parameters:**
  * **p** – The crosslinker conversion
  * **r** – The stoichiometric imbalance
  * **xlink_concentration_0** – [A_f]_0, in 1/volume units
  * **ureg** – The unit registry to use
  * **f** – The functionality of the crosslinkers
  * **alpha** – $P(F_a^{out})$, optional
  * **temperature** – The temperatures; defaults to room temperature (25 °C)
* **Returns:**
  The junction modulus contribution

### pylimer_tools.calc.miller_macosko_theory.compute_miller_macosko_probabilities(r: [float](https://docs.python.org/3/library/functions.html#float), p: [float](https://docs.python.org/3/library/functions.html#float), f: [int](https://docs.python.org/3/library/functions.html#int), b2: [float](https://docs.python.org/3/library/functions.html#float) = 1.0)

Compute Macosko and Miller’s probabilities $P(F_A)$ and $P(F_B)$
i.e., the probability that a randomly chosen A (crosslink) or B (strand-end),
respectively, is the start of a finite chain.

Sources:
: - Macosko and Miller [[MM76a](acknowledgements.md#id5)]
  - Miller and Macosko [[MM76b](acknowledgements.md#id6)]
  - Patel *et al.* [[PMC+92](acknowledgements.md#id14)] (with monofunctional chains, f = 4)
  - Urayama *et al.* [[UMTK04](acknowledgements.md#id18)] (with monofunctional chains, f = 3)

#### NOTE
Currently, only systems with B_2, B_1 and A_f are supported.
If you have a system with other functionality distributions,
you would need to implement these formulas yourselves.

* **Parameters:**
  * **r** – The stoichiometric imbalance
  * **p** – The extent of reaction in terms of the crosslinkers
  * **f** – The functionality of the the crosslinker
  * **b2** – The fraction of bifunctional chains; defaults to 1.0 for no monofunctional chains.
    Can be computed e.g. as $b_2 = \frac{2 \cdot [B_2]}{[B_1] + 2 \cdot [B_2]}$
* **Returns:**
  alpha: $P(F_A)$
* **Returns:**
  beta: $P(F_B)$

### pylimer_tools.calc.miller_macosko_theory.compute_modulus_decomposition(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) | [None](https://docs.python.org/3/library/constants.html#None) = None, ureg: [pint.UnitRegistry](https://pint.readthedocs.io/en/stable/api/base.html#pint.UnitRegistry) | [None](https://docs.python.org/3/library/constants.html#None) = None, unit_style: [None](https://docs.python.org/3/library/constants.html#None) | [UnitStyle](pylimer_tools.io.md#pylimer_tools.io.unit_styles.UnitStyle) = None, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, r: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, p: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, f: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, nu: [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity) | [None](https://docs.python.org/3/library/constants.html#None) = None, temperature: [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity) | [None](https://docs.python.org/3/library/constants.html#None) = None, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, g_e_1: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, b2: [float](https://docs.python.org/3/library/functions.html#float) = 1.0) → Tuple[[pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity), [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity), [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity), [pint.Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity)]

Compute four different estimates of the plateau modulus, using MMT, ANM and PNM.

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **ureg** – The unit registry to use
  * **unit_style** – The unit style to use to have the results in appropriate units
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
  * **r** – The stoichiometric imbalance. Optional if network is specified
  * **p** – The extent of reaction. Optional if network is specified
  * **f** – The functionality of the the crosslinker. Optional if network is specified
  * **nu** – The strand number density (nr of strands per volume) (ideally with units).
    Optional if network is specified
  * **temperature** – The temperature to compute the modulus at. Default: 298.15 K
    Optional, can be passed to improve performance
  * **functionality_per_type** – a dictionary with key: type, and value: functionality of this atom type.
    Optional, can be passed to improve performance
  * **g_e_1** – The melt entanglement modulus
  * **b2** – The mole fraction of reactive sites in B2 among all reactive sites
    in a mixture of B1 and B2
* **Returns:**
  G_MMT_phantom: The phantom contribution to the MMT modulus;
  see also [`pylimer_tools.calc.miller_macosko_theory.compute_junction_modulus()`](#pylimer_tools.calc.miller_macosko_theory.compute_junction_modulus)
* **Returns:**
  G_MMT_entanglement: The entanglement contribution to the MMT modulus
* **Returns:**
  g_anm: The ANM estimate of the modulus
* **Returns:**
  g_pnm: The PNM estimate of the modulus

### pylimer_tools.calc.miller_macosko_theory.compute_probability_that_bifunctional_monomer_is_dangling(p_f_b_out: [float](https://docs.python.org/3/library/functions.html#float)) → [float](https://docs.python.org/3/library/functions.html#float)

Consider a copolymerization of A_f with B_2.
This function computes the probability that a random B_2 unit will be dangling.

Source:
: - Eq. 6.1 in Miller *et al.* [[MVM79](acknowledgements.md#id16)]

* **Parameters:**
  **p_f_b_out** – $P(F_B^{out})$
* **Returns:**
  The probability that a bifunctional monomer is dangling

### pylimer_tools.calc.miller_macosko_theory.compute_probability_that_bifunctional_monomer_is_effective(p_f_b_out: [float](https://docs.python.org/3/library/functions.html#float)) → [float](https://docs.python.org/3/library/functions.html#float)

Consider a copolymerization of A_f with B_2.
This function computes the probability that a random B_2 unit will be effective.

* **Parameters:**
  **p_f_b_out** – $P(F_B^{out})$
* **Returns:**
  The probability that a bifunctional monomer is effective

### pylimer_tools.calc.miller_macosko_theory.compute_probability_that_crosslink_is_dangling(functionality_of_monomer: [int](https://docs.python.org/3/library/functions.html#int), p_f_a_out: [float](https://docs.python.org/3/library/functions.html#float)) → [float](https://docs.python.org/3/library/functions.html#float)

Consider a copolymerization of A_f with B_2.
This function computes the probability that a random A_f unit will be dangling (pendant).
This is equal to the probability that only one of the arms is attached to the gel.

Source:
: - Eq. 6.2 in Miller *et al.* [[MVM79](acknowledgements.md#id16)]

* **Parameters:**
  * **functionality_of_monomer** – f
  * **p_f_a_out** – $P(F_A^{out})$
* **Returns:**
  The probability that a crosslink is dangling

### pylimer_tools.calc.miller_macosko_theory.compute_probability_that_crosslink_is_effective(functionality_of_monomer: [int](https://docs.python.org/3/library/functions.html#int), expected_degree_of_effect: [int](https://docs.python.org/3/library/functions.html#int), p_f_a_out: [float](https://docs.python.org/3/library/functions.html#float)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the probability that an Af, monomer will be an effective crosslink of exactly degree m

$P(X_m^f) = \binom{f}{m} [P(F_A^{out})]^{f-m}[1-P(F_A^{out})]^m$

Source:
: - Eq. 45 in Miller and Macosko [[MM76b](acknowledgements.md#id6)]

* **Parameters:**
  * **functionality_of_monomer** – f
  * **expected_degree_of_effect** – m
  * **p_f_a_out** – $P(F_A^{out})$
* **Returns:**
  The probability that a crosslink is effective

### pylimer_tools.calc.miller_macosko_theory.compute_probability_that_crosslink_with_degree_is_dangling(functionality_of_monomer: [int](https://docs.python.org/3/library/functions.html#int), degree_of_ineffectiveness: [int](https://docs.python.org/3/library/functions.html#int), p_f_a_out: [float](https://docs.python.org/3/library/functions.html#float)) → [float](https://docs.python.org/3/library/functions.html#float)

Consider a copolymerization of A_f with B_2.
This function computes the probability that a random A_f unit will have i pendant arms.

Source:
: - Eq. 6.3 in Miller *et al.* [[MVM79](acknowledgements.md#id16)]

* **Parameters:**
  * **functionality_of_monomer** – f
  * **degree_of_ineffectiveness** – i
  * **p_f_a_out** – $P(F_A^{out})$
* **Returns:**
  The probability that a crosslink with degree is dangling

### pylimer_tools.calc.miller_macosko_theory.compute_trapping_factor(beta: [float](https://docs.python.org/3/library/functions.html#float)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the Langley trapping factor $T_e$.

Literature: Langley [[Lan68](acknowledgements.md#id22)]

* **Parameters:**
  **beta** – $P(F_b^{out})$, see [`compute_miller_macosko_probabilities()`](#pylimer_tools.calc.miller_macosko_theory.compute_miller_macosko_probabilities)
* **Returns:**
  The Langley trapping factor

### pylimer_tools.calc.miller_macosko_theory.compute_weight_fraction_of_backbone(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) | [None](https://docs.python.org/3/library/constants.html#None) = None, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, weight_fractions: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, r: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, p: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, b2: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of the backbone (elastically effective) strands in an infinite network

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
  * **functionality_per_type** – a dictionary with key: atom type, and value: functionality atoms with this type.
  * **weight_fractions** – a dictionary with the weight fraction of each type of atom
  * **r** – The stoichiometric imbalance
  * **p** – The extent of reaction in terms of the crosslinkers
  * **b2** – The mole fraction of reactive sites in B2 among all reactive sites
    in a mixture of B1 and B2
* **Returns:**
  $\Phi_{el} = w_e$: weight fraction of network backbone

### pylimer_tools.calc.miller_macosko_theory.compute_weight_fraction_of_dangling_chains(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) | [None](https://docs.python.org/3/library/constants.html#None) = None, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, weight_fractions: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, r: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, p: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, b2: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of dangling (pendant) strands in infinite network

Source:
: - Eq. 6.4 in Miller *et al.* [[MVM79](acknowledgements.md#id16)]

* **Parameters:**
  * **network** – The network to compute the weight fraction for
  * **crosslinker_type** – The atom type to use to split the molecules
  * **functionality_per_type** – a dictionary with key: type, and value: functionality of this atom type.
  * **weight_fractions** – a dictionary with the weight fraction of each type of atom
  * **r** – The stoichiometric imbalance
  * **p** – The extent of reaction in terms of the crosslinkers
  * **b2** – The mole fraction of reactive sites in B2 among all reactive sites
    in a mixture of B1 and B2
* **Returns:**
  weightFraction $$\Phi_d = w_p$$: weightDangling/weightTotal

### pylimer_tools.calc.miller_macosko_theory.compute_weight_fraction_of_soluble_material(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) | [None](https://docs.python.org/3/library/constants.html#None) = None, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, weight_fractions: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, r: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, p: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, b2: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of soluble material by MMT.

Source:
: - Patel *et al.* [[PMC+92](acknowledgements.md#id14)]
  - Aoyama *et al.* [[AYU21](acknowledgements.md#id13)]

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
  * **weight_fractions** – a dictionary with key: type, and value: weight fraction of type.
    Pass if you want to omit the network.
  * **functionality_per_type** – a dictionary with key: type, and value: functionality of this atom type.
    See: [`determine_functionality_per_type()`](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.determine_functionality_per_type).
  * **r** – The stoichiometric imbalance
  * **p** – The extent of reaction in terms of the crosslinkers
  * **b2** – The mole fraction of reactive sites in B2 among all reactive sites
    in a mixture of B1 and B2
* **Returns:**
  $W_{sol}$ (float): The weight fraction of soluble material according to MMT.
* **Returns:**
  weight_fractions (dict): a dictionary with key: type, and value: weight fraction of type
* **Returns:**
  $\alpha$ (float): Macosko & Miller’s $P(F_A)$
* **Returns:**
  $\beta$ (float): Macosko & Miller’s $P(F_B)$

### pylimer_tools.calc.miller_macosko_theory.compute_weight_fraction_of_soluble_material_from_weight_fractions(r: [float](https://docs.python.org/3/library/functions.html#float), p: [float](https://docs.python.org/3/library/functions.html#float), f: [int](https://docs.python.org/3/library/functions.html#int), w_f: [float](https://docs.python.org/3/library/functions.html#float), w_g: [float](https://docs.python.org/3/library/functions.html#float), g: [int](https://docs.python.org/3/library/functions.html#int) = 2)

Use MMT to compute the weight fraction of soluble material using

$$
`W_{sol} = w_A_f P(F_A^{out})^f + w_B_g [rpP(F_A^{out})^{f-1}+1-rp]^g`

$$

* **Parameters:**
  * **r** – The stoichiometric imbalance
  * **p** – The extent of reaction in terms of the crosslinkers
  * **f** – The functionality of the the crosslinker
  * **w_f** – The weight fraction of the crosslinkers
  * **w_g** – The weight fraction of ordinary chains
  * **g** – The functionality of the ordinary chains

### pylimer_tools.calc.miller_macosko_theory.predict_gelation_point(r: [float](https://docs.python.org/3/library/functions.html#float), f: [int](https://docs.python.org/3/library/functions.html#int), b2: [float](https://docs.python.org/3/library/functions.html#float) = 1) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the gelation point $p_{gel}$ as theoretically predicted
(gelation point = critical extent of reaction for gelation)

Source:
: - Venkataraman *et al.* [[VCC+89](acknowledgements.md#id19)]

* **Parameters:**
  * **r** – The stoichiometric imbalance of reactants (see: #compute_stoichiometric_imbalance)
  * **f** – functionality of the crosslinkers
  * **b2** – The mole fraction of reactive sites in B2 among all reactive sites
    in a mixture of B1 and B2
* **Returns:**
  p_gel: critical extent of reaction for gelation

### pylimer_tools.calc.miller_macosko_theory.predict_maximum_p(r: [float](https://docs.python.org/3/library/functions.html#float), f: [int](https://docs.python.org/3/library/functions.html#int), b2: [float](https://docs.python.org/3/library/functions.html#float) = 1) → [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)

Compute the maximum crosslinker conversion possible given a stoichiometric inbalance.

* **Parameters:**
  * **r** – The stoichiometric imbalance of reactants (see: #compute_stoichiometric_imbalance)
  * **f** – functionality of the crosslinkers
  * **b2** – The mole fraction of reactive sites in B2 among all reactive sites
    in a mixture of B1 and B2. Since r already includes the number of active sites, this argument
    is not necessary.
* **Returns:**
  p_max: The maximum crosslinker conversion possible

### pylimer_tools.calc.miller_macosko_theory.predict_number_density_of_junction_points(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) | [None](https://docs.python.org/3/library/constants.html#None) = None, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the number density of network strands using MMT

Source:
: - Aoyama *et al.* [[AYU21](acknowledgements.md#id13)] (see supporting information for formulae)

* **Parameters:**
  * **network** – The network to compute the weight fraction for
  * **crosslinker_type** – The atom type to use to split the molecules
  * **functionality_per_type** – a dictionary with key: type, and value: functionality of this atom type.
* **Returns:**
  mu: The predicted number density of junction points

### pylimer_tools.calc.miller_macosko_theory.predict_number_density_of_network_strands(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) | [None](https://docs.python.org/3/library/constants.html#None) = None, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, r: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, p: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the number density of network strands using MMT

Source:
: - Aoyama *et al.* [[AYU21](acknowledgements.md#id13)] (see supporting information for formulae)

* **Parameters:**
  * **network** – The network to compute the weight fraction for
  * **crosslinker_type** – The atom type to use to split the molecules
  * **functionality_per_type** – a dictionary with key: type, and value: functionality of this atom type.
  * **r** – The stoichiometric imbalance
  * **p** – The extent of reaction in terms of the crosslinkers
* **Returns:**
  nu: The predicted number density of network strands

### pylimer_tools.calc.miller_macosko_theory.predict_p_from_w_sol(w_sol: [float](https://docs.python.org/3/library/functions.html#float), network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) | [None](https://docs.python.org/3/library/constants.html#None) = None, crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, weight_fractions: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, r: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, b2: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Compute the extent of reaction based on the weight fraction of soluble material.

* **Parameters:**
  * **w_sol** – The weight fraction of soluble material
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
  * **functionality_per_type** – a dictionary with key: type, and value: functionality of this atom type
  * **weight_fractions** – a dictionary with the weight fraction of each type of atom
  * **r** – The stoichiometric imbalance
  * **b2** – The mole fraction of reactive sites in B2 among all reactive sites
    in a mixture of B1 and B2
* **Returns:**
  The extent of reaction p

### pylimer_tools.calc.miller_macosko_theory.predict_shear_modulus(\*\*kwargs) → [Quantity](https://pint.readthedocs.io/en/stable/api/base.html#pint.Quantity)

Predict the shear modulus using MMT Analysis.

Source:
: - Gusev and Schwarz [[GS22](acknowledgements.md#id21)]

* **Parameters:**
  **kwargs** – See [`compute_modulus_decomposition()`](#pylimer_tools.calc.miller_macosko_theory.compute_modulus_decomposition)
* **Returns:**
  G: The predicted shear modulus

## pylimer_tools.calc.structure_analysis module

This module provides functions to compute various quantities related to polymer networks from their structure.

### pylimer_tools.calc.structure_analysis.compute_crosslinker_conversion(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, f: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the extent of reaction of the crosslinkers
(actual functionality divided by target functionality)

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
  * **f** – The (target) functionality of the crosslinkers
  * **functionality_per_type** – A dictionary with key: type, and value: (target) functionality of this atom type
* **Returns:**
  The (mean) crosslinker conversion
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.structure_analysis.compute_effective_crosslinker_functionalities(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2) → [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

Compute the functionality of every crosslinker in the network

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
* **Returns:**
  The functionality of every crosslinker
* **Return type:**
  [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)]

### pylimer_tools.calc.structure_analysis.compute_effective_crosslinker_functionality(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the mean crosslinker functionality

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
* **Returns:**
  The (mean) effective crosslinker functionality
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.structure_analysis.compute_end_to_end_vectors(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

Compute the end to end vectors between each pair of (indirectly) connected crosslinker

* **Parameters:**
  * **network** ([*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)) – The polymer network to do the computation for
  * **crosslinker_type** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The atom type to compute the in-between vectors for
* **Returns:**
  a dictionary with key: “{atom1.name}+{atom2.name}” and value: Their difference vector
* **Return type:**
  [dict](https://docs.python.org/3/library/stdtypes.html#dict)

### pylimer_tools.calc.structure_analysis.compute_extent_of_reaction(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the extent of polymerization reaction
(nr. of formed bonds in reaction / max. nr. of bonds formable)

#### NOTE
- if your system has a non-integer number of possible bonds (e.g. one site non-bonded),
  this will not be rounded/respected in any way
- if the system contains solvent or other molecules that should not be binding to
  crosslinkers, make sure to remove them before calling this function
- if you have multiple crosslinker types, be sure to replace them by one type only

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The atom type of crosslinker beads
  * **functionality_per_type** – A dictionary with key: type, and value: functionality of this atom type.
    If None: will use max functionality per type.
* **Returns:**
  The extent of reaction
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.structure_analysis.compute_fraction_of_bifunctional_reactive_sites(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the mole fraction of reactive sites in B2
among all reactive sites in a mixture of B1 and B2.
Uses the network to detect what might be monofunctional chains,
and then counts them and the bifunctional ones.

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The atom type of crosslinker beads
  * **functionality_per_type** – A dictionary with key: type, and value: functionality of this atom type
* **Returns:**
  The mole fraction of reactive sites in B2 among all reactive sites in a mixture of B1 and B2
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.structure_analysis.compute_mean_end_to_end_distances(networks: Sequence[[Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)], crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

Compute the mean end to end distance between each pair of (indirectly) connected crosslinker

* **Parameters:**
  * **networks** (*Sequence* *[*[*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) *]*) – The different configurations of the polymer network to do the computation for
  * **crosslinker_type** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The atom type to compute the in-between vectors for
* **Returns:**
  a dictionary with key: “{atom1.name}+{atom2.name}” and value: The norm of the mean difference vector
* **Return type:**
  [dict](https://docs.python.org/3/library/stdtypes.html#dict)

### pylimer_tools.calc.structure_analysis.compute_mean_end_to_end_vectors(networks: Sequence[[Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)], crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

Compute the mean end to end vectors between each pair of (indirectly) connected crosslinker

* **Parameters:**
  * **networks** (*Sequence* *[*[*Universe*](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe) *]*) – The different configurations of the polymer network to do the computation for
  * **crosslinker_type** ([*int*](https://docs.python.org/3/library/functions.html#int)) – The atom type to compute the in-between vectors for
* **Returns:**
  a dictionary with key: “{atom1.name}+{atom2.name}” and value: Their mean distance difference vector
* **Return type:**
  [dict](https://docs.python.org/3/library/stdtypes.html#dict)

### pylimer_tools.calc.structure_analysis.compute_stoichiometric_imbalance(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, functionality_per_type: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None) = None, ignore_types: [list](https://docs.python.org/3/library/stdtypes.html#list) | [None](https://docs.python.org/3/library/constants.html#None) = None, effective: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the stoichiometric imbalance
( nr. of bonds formable of crosslinker / (nr. of bonds formable of precursor chains) )

r > 1 means an excess of crosslinkers, whereas r = 0 implies that there are not crosslinkers in the network.

#### NOTE
- if your system has a non-integer number of possible bonds (e.g. one site non-bonded),
  this will not be rounded/respected in any way.
- Currently, only g = 2 chains are respected yet, and only one type of crosslinker is supported.
- If you have e.g. solvent chains in the network, use ignore_types to not count them.

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
  * **functionality_per_type** – A dictionary with key: type, and value: functionality of this atom type.
    If None: will use max functionality per type.
  * **ignore_types** – A list of integers, the types to ignore for the imbalance (e.g. solvent atom types)
  * **effective** – Whether to use the effective functionality (if functionality_per_type is not passed)
    or the maximum functionality of the crosslinkers.
* **Returns:**
  The stoichiometric imbalance. If the network is empty, or no crosslinkers are present 0. is returned.
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.structure_analysis.compute_weight_fractions(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

Compute the weight fractions of each atom type in the network.
Kept for compatibility reasons; just call
[`pylimer_tools_cpp.Universe.compute_weight_fractions()`](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe.compute_weight_fractions) instead.

* **Parameters:**
  **network** – The polymer network to do the computation for
* **Returns:**
  Using the type i as a key, this dict contains the weight fractions ($\frac{W_i}{W_{tot}}$)
* **Return type:**
  [dict](https://docs.python.org/3/library/stdtypes.html#dict)

### pylimer_tools.calc.structure_analysis.measure_lower_bound_weight_fraction_of_soluble_material(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2, rel_tol: [float](https://docs.python.org/3/library/functions.html#float) = 0.75, abs_tol: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute a lower bound on the weight fraction of soluble material by counting.
This is ones as such: only clusters, which do not contain loops and are smaller than the rel_tol of the biggest,
are counted as soluble

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **crosslinker_type** – The type of the junctions/crosslinkers to select them in the network
  * **rel_tol** – The fraction of the maximum weight that counts as soluble. Ignored if abs_tol is specified
  * **abs_tol** – The weight from which on a component is not soluble anymore
* **Returns:**
  The weight fraction of soluble material as counted. 0 for an empty network
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_backbone(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2)

Compute the weight fraction of network backbone in infinite network

* **Parameters:**
  * **network** – The network to compute the weight fraction for
  * **crosslinker_type** – The atom type to use to split the molecules
* **Returns:**
  1 - weightDangling/weightTotal - weightSoluble/weightTotal
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

#### SEE ALSO
- [`pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_dangling_chains()`](#pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_dangling_chains)
- [`pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_soluble_material()`](#pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_soluble_material)

### pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_dangling_chains(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), crosslinker_type: [int](https://docs.python.org/3/library/functions.html#int) = 2) → Tuple[[float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float)]

Compute the weight fraction of dangling strands in infinite network

#### NOTE
Currently, only primary dangling chains are taken into account.
There are other methods that incorporate more.

* **Parameters:**
  * **network** – The network to compute the weight fraction for
  * **crosslinker_type** – The atom type to use to split the molecules
* **Returns:**
  A tuple containing (weightDangling/weightTotal, numDangling/numTotal)
* **Return type:**
  Tuple[[float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float)]

### pylimer_tools.calc.structure_analysis.measure_weight_fraction_of_soluble_material(network: [Universe](api/pylimer_tools_cpp.Universe.md#pylimer_tools_cpp.Universe), rel_tol: [float](https://docs.python.org/3/library/functions.html#float) = 0.5, abs_tol: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the weight fraction of soluble material by counting.
Effectively, this method counts the weight of clusters
that have a weight less than a certain fraction of the total weight.

* **Parameters:**
  * **network** – The polymer network to do the computation for
  * **rel_tol** – The fraction of the maximum weight that counts as soluble. Ignored if abs_tol is specified
  * **abs_tol** – The weight from which on a component is not soluble anymore
* **Returns:**
  The weight fraction of soluble material as counted. 0 for an empty network
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

## Module contents
