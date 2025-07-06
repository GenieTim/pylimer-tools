Force Balance & Maximum Entropy Homogenization Procedures
=========================================================

pylimer-tools provides three implementations of the Maximum Entropy Homogenization Procedure (MEHP) to reduce polymer networks to their minimum energy, maximum entropy homogenized state. 
This is useful for predicting the equilibrium shear modulus of the network.
The three implementations are:
- :class:`~pylimer_tools_cpp.MEHPForceRelaxation`: An MEHP implementation that allows to use non-linear force potentials,
- :class:`~pylimer_tools_cpp.MEHPForceBalance`: The original Force Balance implementation, which uses just Hookean springs, but allows the use of slip-links to model entanglements,
- :class:`~pylimer_tools_cpp.MEHPForceBalance2`: The faster implementation of the Force Balance procedure, without allowing for slip-links, but allowing the modelling of entanglements as static links, like tetrafunctional cross-links.
