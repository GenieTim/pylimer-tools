r"""
Inactive Removal & Sampling
===========================

Networks below full conversion contain **dangling chains**, **sol fractions**, and **loops**.
In :class:`~pylimer_tools_cpp.MEHPForceBalance2`, entanglements do not slip off chain ends,
so two complementary controls prevent dangling arms from artificially stiffening the network:

1. **Sampling-time topological filtering** (:class:`~pylimer_tools_cpp.EntanglementFilteringMode`):
   ``OMIT_INACTIVE`` (default) or ``PRUNE_UPFRONT`` avoids placing entanglements on dangling/soluble
   strands. Passing ``NONE`` places links everywhere, inflating :math:`G` by >300–4000% near
   percolation unless pruned by the solver.
2. **Solver-level simplification** (:class:`~pylimer_tools_cpp.StructureSimplificationMode`):
   ``INACTIVE_ONLY`` prunes zero-force dangling strands/entanglements; ``INACTIVE_THEN_X2F``
   also collapses 2-functional nodes, compressing the sparse linear system by up to 80%+.
"""

from pathlib import Path
import sys

if "__file__" in globals():
    _examples_dir = str(Path(__file__).resolve().parents[1])
    _fb_dir = str(Path(__file__).resolve().parent)
    for _p in (_examples_dir, _fb_dir):
        if _p not in sys.path:
            sys.path.insert(0, _p)

from _example_cache import get_example_cache_path, load_or_compute_cache  # noqa: E402
from _plotting import (  # noqa: E402
    plot_inactive_modulus_convergence,
    plot_inactive_relative_error,
    plot_inactive_solver_complexity,
    plot_inactive_topological_fractions,
)
from _study_runners import run_inactive_removal_and_sampling_simulation  # noqa: E402
import pylimer_tools_cpp as pt  # noqa: E402

# %%
# Comparing Filtering & Simplification Modes on a Defective Network
# -----------------------------------------------------------------
# At :math:`p = 0.80`, comparing ``EntanglementFilteringMode.NONE`` + ``NO_SIMPLIFICATION``
# against ``OMIT_INACTIVE`` + ``INACTIVE_THEN_X2F`` demonstrates how inactive pruning
# eliminates unphysical dead-end constraints and shrinks the linear system.

gen = pt.MCUniverseGenerator(18.0, 18.0, 18.0)
gen.set_seed(42)
gen.config_nr_of_mc_steps(0)
gen.add_strands(200, [50] * 200)
gen.add_crosslinkers(100, 4)
gen.link_strands_to_conversion(0.80)

es_unfilt = gen.sample_entanglements(
    250, upper_cutoff=2.6, seed="demo", filtering_mode=pt.EntanglementFilteringMode.NONE
)
fb_naive = gen.get_force_balance2(es_unfilt)
fb_naive.run_force_relaxation(pt.StructureSimplificationMode.NO_SIMPLIFICATION)

es_filt = gen.sample_entanglements(
    250,
    upper_cutoff=2.6,
    seed="demo",
    filtering_mode=pt.EntanglementFilteringMode.OMIT_INACTIVE,
)
fb_clean = gen.get_force_balance2(es_filt)
fb_clean.run_force_relaxation(pt.StructureSimplificationMode.INACTIVE_THEN_X2F)

print(
    f"Springs — naive: {fb_naive.get_nr_of_springs()}, "
    f"simplified: {fb_clean.get_nr_of_springs()}"
)

cache_path = get_example_cache_path(
    "inactive_removal_and_sampling_cache.npz", "force_balance"
)
data = load_or_compute_cache(
    cache_path,
    run_inactive_removal_and_sampling_simulation)

# %%
# Shear Modulus Convergence & The Dangling-Chain Trap
# ---------------------------------------------------
# Without filtering or simplification, pinned dangling chains inflate :math:`G` by hundreds to
# thousands of percent below :math:`p \approx 0.90`. Either pre-filtering or solver-level
# inactive removal restores the true physical percolation curve.

plot_inactive_modulus_convergence(data)
plot_inactive_relative_error(data)

# %%
# Topological Compression & Solver Complexity
# -------------------------------------------
# ``INACTIVE_THEN_X2F`` removes inactive springs and collapses 2-functional nodes, cutting
# active system size and solver iterations substantially near the gelation
# threshold.

plot_inactive_topological_fractions(data)
plot_inactive_solver_complexity(data)
