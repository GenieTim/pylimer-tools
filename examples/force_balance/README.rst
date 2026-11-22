Force Balance, Maximum Entropy Homogenization Procedures
--------------------------------------------------------

pylimer-tools provides a unified implementation of the Force Balance :cite:p:`bernhard_phantom_2025`, Maximum Entropy Homogenization Procedure (MEHP) :cite:p:`gusev_numerical_2019` to reduce polymer networks to their minimum energy, maximum entropy homogenized state. 
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

The unified architecture coordinates a canonical :class:`~pylimer_tools_cpp.MEHPNetwork` (supporting phantom networks, fixed merged entanglement links via :attr:`~pylimer_tools_cpp.EntanglementMode.MERGED_LINKS`, and compliant coupling springs via :attr:`~pylimer_tools_cpp.EntanglementMode.SPRING`), a stateless :class:`~pylimer_tools_cpp.MEHPForceEvaluator` (:class:`~pylimer_tools_cpp.HookeanForceEvaluator`, :class:`~pylimer_tools_cpp.NonGaussianSpringForceEvaluator`, or custom Python subclasses), and interchangeable :class:`~pylimer_tools_cpp.MEHPSolver` strategies inside :class:`~pylimer_tools_cpp.MEHPSimulator`:

- :class:`~pylimer_tools_cpp.LinearSystemMEHPSolver` (:class:`~pylimer_tools_cpp.MEHPLinearSolverConfig`): The fastest solver for linear Hookean networks with fixed entanglement links (:attr:`~pylimer_tools_cpp.SlipMode.FIXED`) or springs, assembling and solving a sparse graph-Laplacian linear system (including optional GPU Metal backends).
- :class:`~pylimer_tools_cpp.SlipLinkRelaxationSolver` (:class:`~pylimer_tools_cpp.MEHPSlipLinkSolverConfig`): An iterative alternating coordinate and contour-partition relaxation solver for Hookean networks that supports dynamic contour slipping (:attr:`~pylimer_tools_cpp.SlipMode.SLIDING`) and slip-link swapping as well as fixed links (:attr:`~pylimer_tools_cpp.SlipMode.FIXED`).
- :class:`~pylimer_tools_cpp.NonlinearOptimizationSolver` (:class:`~pylimer_tools_cpp.MEHPNonlinearSolverConfig`): A general NLopt-based energy minimization solver supporting linear and non-linear constitutive force potentials (e.g., :class:`~pylimer_tools_cpp.NonGaussianSpringForceEvaluator` or custom Python :class:`~pylimer_tools_cpp.MEHPForceEvaluator` subclasses) on phantom and fixed-entanglement networks.

Sampling Entanglements
~~~~~~~~~~~~~~~~~~~~~~

Entanglements can be sampled in two ways:

1. **On-demand bounding-box sampling** (default): Computes per-strand bounding boxes and regenerates candidate strand beads on demand, matching discrete bead-level fidelity with ~80–95% lower peak memory.
2. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a :class:`~pylimer_tools_cpp.Universe`.

Both methods support :class:`~pylimer_tools_cpp.EntanglementFilteringMode` for topological leaf peeling of dangling and soluble fractions.

See:

- :doc:`/auto_examples/force_balance/plot_entanglement_sampling_comparison` for a detailed benchmark and comparison of discrete vs. on-demand sampling methods,
- :doc:`/auto_examples/force_balance/plot_entanglement_link_vs_spring` for the difference between modeling entanglements as rigid merged links (:attr:`~pylimer_tools_cpp.EntanglementMode.MERGED_LINKS`) versus compliant Hookean springs (:attr:`~pylimer_tools_cpp.EntanglementMode.SPRING`),
- :doc:`/auto_examples/force_balance/plot_slipping_vs_fixed_links` for a comparison of dynamic contour slipping (:attr:`~pylimer_tools_cpp.SlipMode.SLIDING`) versus fixed entanglement links (:attr:`~pylimer_tools_cpp.SlipMode.FIXED`) using :class:`~pylimer_tools_cpp.SlipLinkRelaxationSolver`,
- :doc:`/auto_examples/force_balance/plot_inactive_removal_and_sampling` for the removal of inactive strands/entanglements via :class:`~pylimer_tools_cpp.StructureSimplificationMode` and the impact of pre-filtering dangling chains during entanglement sampling,
- :doc:`/auto_examples/force_balance/plot_polydispersity_study` for the effect of precursor strand length polydispersity on equilibrium shear modulus and defect fractions,
- :doc:`/auto_examples/force_balance/plot_sampling_cutoff_effect` for the effect of the slip-link sampling cut-off distance on the equilibrium shear modulus, investigating pre-strain distortions, link vs. spring representations, and interactions with strand length :math:`N`, conversion :math:`p`, and stoichiometric imbalance :math:`r`, and
- :doc:`/auto_examples/force_balance/plot_representative_volume_pdms` for a representative volume element (RVE) study of realistic PDMS networks across chain length :math:`N`, conversion :math:`p`, stoichiometric imbalance :math:`r`, and ensemble sample count.
