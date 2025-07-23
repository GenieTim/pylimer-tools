Network Generation
------------------

pylimer-tools includes a powerful Monte Carlo-based network generator for creating realistic polymer systems. 
This tools is particularly useful for generating quick initial configurations for molecular dynamics simulations,
or for doing fast predictions of the viscoelasticity using the :class:`~pylimer_tools_cpp.NormalModeAnalyzer` or the 
:class:`pylimer_tools_cpp.MEHPForceBalance` and its variants.

For more options, examples and advanced usage, see the `test suite <https://github.com/GenieTim/pylimer-tools/tree/main/tests>`_ 
and the :doc:`modules` API reference, particularly of the :class:`~pylimer_tools_cpp.MCUniverseGenerator` class.

The network generation system allows you to:

- Create end-linked and vulcanized polymer networks with controllable topology
- Generate solvent chains and background molecules
- Control crosslinking degree and functionality
- Set realistic bead distances and box sizes

The :class:`~pylimer_tools_cpp.MCUniverseGenerator` class is the main interface for network generation.
The initial method was introduced in :cite:t:`gusev_numerical_2019`,
the current implementation is more powerful, as needed for :cite:t:`bernhard_phantom_2025`.
