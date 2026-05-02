<a id="sphx-glr-auto-examples-force-balance"></a>

# Force Balance, Maximum Entropy Homogenization Procedures

pylimer-tools provides three implementations of the Force Balance:cite:p:bernhard_phantom_2025, Maximum Entropy Homogenization Procedure (MEHP):cite:p:gusev_numerical_2019 to reduce polymer networks to their minimum energy, maximum entropy homogenized state.
This is useful for predicting the (phantom) equilibrium shear modulus of the network.

Here, we distinguish the following three different implementations:

- [`MEHPForceRelaxation`](../../api/pylimer_tools_cpp.MEHPForceRelaxation.md#pylimer_tools_cpp.MEHPForceRelaxation): An MEHP implementation that allows to use non-linear force potentials,
- [`MEHPForceBalance`](../../api/pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance): A Force Balance implementation, which uses just Hookean springs, but allows the use of slip-links to model entanglements,
- [`MEHPForceBalance2`](../../api/pylimer_tools_cpp.MEHPForceBalance2.md#pylimer_tools_cpp.MEHPForceBalance2): The faster implementation of the Force Balance procedure, without allowing for slip-links, but allowing the modelling of entanglements as static links, like tetrafunctional crosslinks.

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="MEHPForceBalance2 is a faster implementation of the MEHP force balance method, which uses static links to model entanglements instead of slip-links. See :citebernhard_phantom_2025">  <div class="sphx-glr-thumbnail-title">Force Balance 2</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the MEHP implementation with slip-springs to homogenize a polymer network and predict its equilibrium shear modulus.">  <div class="sphx-glr-thumbnail-title">Force Balance</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="In this example, we study the effect of polydispersity on the shear modulus of end-linked polymer networks.">  <div class="sphx-glr-thumbnail-title">Polydispersity Study</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example shows how to use the MEHPForceBalance2 class in a deformation experiment.">  <div class="sphx-glr-thumbnail-title">Deformation Experiment</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use the &quot;original&quot; MEHP implementation (see :citegusev_numerical_2019) without slip-links or entanglements.">  <div class="sphx-glr-thumbnail-title">Maximum Entropy Homogenization Procedure (MEHP)</div>
</div>
<!-- thumbnail-parent-div-close --></div>
