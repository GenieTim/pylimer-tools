<a id="sphx-glr-auto-examples-rouse-simulations"></a>

# Rouse & Brownian Dynamics Simulations

pylimer-tools provides an implicit time-stepping Rouse and Brownian dynamics
simulator for polymer chains and general Hookean spring networks via
[`RouseBrownianSimulator`](../../api/pylimer_tools_cpp.RouseBrownianSimulator.md#pylimer_tools_cpp.RouseBrownianSimulator), returning dynamical observables
in [`RouseSimulationResult`](../../api/pylimer_tools_cpp.RouseSimulationResult.md#pylimer_tools_cpp.RouseSimulationResult).

At each time step $\Delta t$, the positions $\mathbf{r}_t \in \mathbb{R}^{N \times 3}$
are advanced using the unconditionally stable backward-Euler discretization of the
overdamped Langevin equation on the scalar graph Laplacian $\mathbf{L}$:

$$
\left(\mathbf{I} + \Delta t \frac{k}{\gamma} \mathbf{L}\right) \mathbf{r}_{t+\Delta t}
= \mathbf{r}_t + \frac{\Delta t}{\gamma} \mathbf{F}_{\text{ext}} + \sqrt{2 D \Delta t}\,\boldsymbol{\xi}_t,
$$

where $k$ is the spring constant, $\gamma$ is the bead friction coefficient,
$D = k_B T / \gamma$ is the bare bead diffusion coefficient, and
$\boldsymbol{\xi}_t$ has i.i.d. standard normal entries. All three Cartesian
coordinates share a single sparse factorization of
$\mathbf{A} = \mathbf{I} + \Delta t \frac{k}{\gamma} \mathbf{L}$, configurable
through [`RouseSolverType`](../../api/pylimer_tools_cpp.RouseSolverType.md#pylimer_tools_cpp.RouseSolverType) (direct sparse LDLT, direct sparse
LLT, or iterative Conjugate Gradient).

See:

- [Rouse & Brownian Dynamics Simulation](plot_rouse_brownian_simulation.md) for
  deterministic normal-mode relaxation, linear solver comparison, external-force
  drift and creep/recoil, and stochastic Rouse MSD and end-to-end relaxation dynamics
  (including multi-tau analysis with [`Correlator`](../../api/pylimer_tools_cpp.Correlator.md#pylimer_tools_cpp.Correlator)).
- [Polymer Architectures & Resistance-Fluctuation Duality in Rouse Dynamics](plot_rouse_topologies_and_resistance.md) for
  dynamics across linear, ring, star, complete, and weighted topologies and the exact
  fluctuation-resistance duality connecting [`RouseBrownianSimulator`](../../api/pylimer_tools_cpp.RouseBrownianSimulator.md#pylimer_tools_cpp.RouseBrownianSimulator)
  and [`ResistanceDistanceAnalysis`](../../api/pylimer_tools_cpp.ResistanceDistanceAnalysis.md#pylimer_tools_cpp.ResistanceDistanceAnalysis).

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="This example bridges RouseBrownianSimulator and ResistanceDistanceAnalysis, illustrating how topological effective resistances govern the equilibrium conformational fluctuations of Gaussian polymer architectures and crosslinked networks.">![](auto_examples/rouse_simulations/images/thumb/sphx_glr_plot_rouse_topologies_and_resistance_thumb.png)

[Polymer Architectures & Resistance-Fluctuation Duality in Rouse Dynamics](plot_rouse_topologies_and_resistance.md)

  <div class="sphx-glr-thumbnail-title">Polymer Architectures & Resistance-Fluctuation Duality in Rouse Dynamics</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use RouseBrownianSimulator to simulate overdamped polymer chain dynamics with an implicit backward-Euler integrator:">![](auto_examples/rouse_simulations/images/thumb/sphx_glr_plot_rouse_brownian_simulation_thumb.png)

[Rouse & Brownian Dynamics Simulation](plot_rouse_brownian_simulation.md)

  <div class="sphx-glr-thumbnail-title">Rouse & Brownian Dynamics Simulation</div>
</div>
<!-- thumbnail-parent-div-close --></div>
