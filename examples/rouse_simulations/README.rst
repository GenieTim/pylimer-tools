Rouse & Brownian Dynamics Simulations
-------------------------------------

pylimer-tools provides an implicit time-stepping Rouse and Brownian dynamics
simulator for polymer chains and general Hookean spring networks via
:class:`~pylimer_tools_cpp.RouseBrownianSimulator`, returning dynamical observables
in :class:`~pylimer_tools_cpp.RouseSimulationResult`.

At each time step :math:`\Delta t`, the positions :math:`\mathbf{r}_t \in \mathbb{R}^{N \times 3}`
are advanced using the unconditionally stable backward-Euler discretization of the
overdamped Langevin equation on the scalar graph Laplacian :math:`\mathbf{L}`:

.. math::

   \left(\mathbf{I} + \Delta t \frac{k}{\gamma} \mathbf{L}\right) \mathbf{r}_{t+\Delta t}
   = \mathbf{r}_t + \frac{\Delta t}{\gamma} \mathbf{F}_{\text{ext}} + \sqrt{2 D \Delta t}\,\boldsymbol{\xi}_t,

where :math:`k` is the spring constant, :math:`\gamma` is the bead friction coefficient,
:math:`D = k_B T / \gamma` is the bare bead diffusion coefficient, and
:math:`\boldsymbol{\xi}_t` has i.i.d. standard normal entries. All three Cartesian
coordinates share a single sparse factorization of
:math:`\mathbf{A} = \mathbf{I} + \Delta t \frac{k}{\gamma} \mathbf{L}`, configurable
through :class:`~pylimer_tools_cpp.RouseSolverType` (direct sparse LDLT, direct sparse
LLT, or iterative Conjugate Gradient).

See:

- :doc:`/auto_examples/rouse_simulations/plot_rouse_brownian_simulation` for
  deterministic normal-mode relaxation, linear solver comparison, external-force
  drift and creep/recoil, and stochastic Rouse MSD and end-to-end relaxation dynamics
  (including multi-tau analysis with :class:`~pylimer_tools_cpp.Correlator`).
- :doc:`/auto_examples/rouse_simulations/plot_rouse_topologies_and_resistance` for
  dynamics across linear, ring, star, complete, and weighted topologies and the exact
  fluctuation-resistance duality connecting :class:`~pylimer_tools_cpp.RouseBrownianSimulator`
  and :class:`~pylimer_tools_cpp.ResistanceDistanceAnalysis`.
