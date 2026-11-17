Force Balance, Maximum Entropy Homogenization Procedures
--------------------------------------------------------

pylimer-tools provides three implementations of the Force Balance :cite:p:`bernhard_phantom_2025`, Maximum Entropy Homogenization Procedure (MEHP) :cite:p:`gusev_numerical_2019` to reduce polymer networks to their minimum energy, maximum entropy homogenized state. 
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- :class:`~pylimer_tools_cpp.MEHPForceRelaxation`: An MEHP implementation that allows the use of non-linear force potentials,
- :class:`~pylimer_tools_cpp.MEHPForceBalance`: A Force Balance implementation using Hookean springs that supports slip-links with optional dynamic contour slipping (or fixed links via ``disable_slipping=True``),
- :class:`~pylimer_tools_cpp.MEHPForceBalance2`: The faster implementation of the Force Balance procedure, representing entanglements as fixed 4-functional links (like tetrafunctional crosslinks) or compliant Hookean springs, without contour slipping.

Sampling Entanglements
~~~~~~~~~~~~~~~~~~~~~~

Entanglements can be sampled in two ways:

1. **On-demand bounding-box sampling** (default): Computes per-strand bounding boxes and regenerates candidate strand beads on demand, matching discrete bead-level fidelity with ~80–95% lower peak memory.
2. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a :class:`~pylimer_tools_cpp.Universe`.

Both methods support :class:`~pylimer_tools_cpp.EntanglementFilteringMode` for topological leaf peeling of dangling and soluble fractions.

See:

- :doc:`plot_entanglement_sampling_comparison` for a detailed benchmark and comparison of discrete vs. on-demand sampling methods,
- :doc:`plot_entanglement_link_vs_spring` for the difference between modeling entanglements as rigid merged links versus compliant Hookean springs in Force Balance 2,
- :doc:`plot_slipping_vs_fixed_links` for a comparison of dynamic contour slipping versus fixed entanglement links in Force Balance 1, and
- :doc:`plot_inactive_removal_and_sampling` for the removal of inactive strands/entanglements in Force Balance 2 and the impact of pre-filtering dangling chains during entanglement sampling,
- :doc:`plot_polydispersity_study` for the effect of precursor strand length polydispersity on equilibrium shear modulus and defect fractions,
- :doc:`plot_sampling_cutoff_effect` for the effect of the slip-link sampling cut-off distance on the equilibrium shear modulus in Force Balance 2, investigating pre-strain distortions, link vs. spring representations, and interactions with strand length :math:`N`, conversion :math:`p`, and stoichiometric imbalance :math:`r`, and
- :doc:`plot_representative_volume_pdms` for a representative volume element (RVE) study of realistic PDMS networks across chain length :math:`N`, conversion :math:`p`, stoichiometric imbalance :math:`r`, and ensemble sample count.
