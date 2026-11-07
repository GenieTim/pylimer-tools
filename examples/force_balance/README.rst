Force Balance, Maximum Entropy Homogenization Procedures
--------------------------------------------------------

pylimer-tools provides three implementations of the Force Balance:cite:p:`bernhard_phantom_2025`, Maximum Entropy Homogenization Procedure (MEHP):cite:p:`gusev_numerical_2019` to reduce polymer networks to their minimum energy, maximum entropy homogenized state. 
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- :class:`~pylimer_tools_cpp.MEHPForceRelaxation`: An MEHP implementation that allows to use non-linear force potentials,
- :class:`~pylimer_tools_cpp.MEHPForceBalance`: A Force Balance implementation, which uses just Hookean springs, but allows the use of slip-links to model entanglements,
- :class:`~pylimer_tools_cpp.MEHPForceBalance2`: The faster implementation of the Force Balance procedure, without allowing for slip-links, but allowing the modelling of entanglements as static links, like tetrafunctional crosslinks.

Sampling Entanglements
~~~~~~~~~~~~~~~~~~~~~~

Entanglements (slip-links) can be sampled in two distinct ways:

1. **Discrete full-universe sampling**: Directly on explicit bead coordinates in a :class:`~pylimer_tools_cpp.Universe`.
2. **Analytical Brownian-bridge sampling**: Directly on the crosslinkers and strand connectivity of a :class:`~pylimer_tools_cpp.CrosslinkerUniverse` without instantiating intermediate strand beads, reducing peak memory by up to 90%+ and supporting topological leaf peeling of dangling and soluble fractions.

See :doc:`plot_entanglement_sampling_comparison` for a detailed benchmark and comparison of discrete vs. analytical sampling methods, and :doc:`plot_sliplink_vs_slipspring` for the difference between modeling entanglements as merged slip-links versus compliant Hookean slip-springs.
