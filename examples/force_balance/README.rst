Force Balance, Maximum Entropy Homogenization Procedures
--------------------------------------------------------

pylimer-tools provides three implementations of the Force Balance:cite:p:`bernhard_phantom_2025`, Maximum Entropy Homogenization Procedure (MEHP):cite:p:`gusev_numerical_2019` to reduce polymer networks to their minimum energy, maximum entropy homogenized state. 
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- :class:`~pylimer_tools_cpp.MEHPForceRelaxation`: An MEHP implementation that allows to use non-linear force potentials,
- :class:`~pylimer_tools_cpp.MEHPForceBalance`: A Force Balance implementation, which uses just Hookean springs, but allows the use of slip-links to model entanglements,
- :class:`~pylimer_tools_cpp.MEHPForceBalance2`: The faster implementation of the Force Balance procedure, without allowing for slip-links, but allowing the modelling of entanglements as static links, like tetrafunctional crosslinks.
