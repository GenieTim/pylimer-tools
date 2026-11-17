r"""
Slipping vs. Fixed Entanglement Links
=====================================

In :class:`~pylimer_tools_cpp.MEHPForceBalance` (:cite:t:`gusev_numerical_2019`, :cite:t:`bernhard_phantom_2025`),
entanglement links can either remain frozen at their initial contour positions (``disable_slipping=True``)
or slide dynamically along the polymer strands to balance contour tension (``disable_slipping=False``).

Dynamic contour slipping artificially softens the network modulus by relaxing affine tension
and increases solver iterations 2–3x. Consequently, :class:`~pylimer_tools_cpp.MEHPForceBalance2`
fixes entanglement positions by design.
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
    plot_slipping_vs_fixed_iterations,
    plot_slipping_vs_fixed_modulus,
)
from _study_runners import run_slipping_vs_fixed_simulation  # noqa: E402
import pylimer_tools_cpp as pt  # noqa: E402

# %%
# Toggling Dynamic Slipping in ``MEHPForceBalance``
# -------------------------------------------------
# On a generated network, ``disable_slipping=True`` freezes the sampled contour fractions,
# whereas ``disable_slipping=False`` updates them iteratively during force
# relaxation.

gen = pt.MCUniverseGenerator(18.0, 18.0, 18.0)
gen.set_seed(42)
gen.config_nr_of_mc_steps(0)
gen.add_strands(200, [50] * 200)
gen.add_crosslinkers(100, 4)
gen.link_strands_to_conversion(0.95)

es = gen.sample_entanglements(
    nr_of_entanglements_to_sample=300,
    upper_cutoff=2.6,
    seed="demo",
    filtering_mode=pt.EntanglementFilteringMode.OMIT_INACTIVE,
)

fb_fixed = gen.get_force_balance(es)
fb_fixed.run_force_relaxation(max_nr_of_steps=5000, disable_slipping=True)

fb_slip = gen.get_force_balance(es)
fb_slip.run_force_relaxation(max_nr_of_steps=5000, disable_slipping=False)

print(
    f"Iterations — fixed: {fb_fixed.get_nr_of_iterations()}, "
    f"slipping: {fb_slip.get_nr_of_iterations()}"
)

cache_path = get_example_cache_path("slipping_vs_fixed_cache.npz", "force_balance")
data = load_or_compute_cache(cache_path, run_slipping_vs_fixed_simulation)

# %%
# Modulus Reduction & Convergence Cost
# ------------------------------------
# Multi-seed ensemble results (95% confidence intervals) confirm that dynamic slipping reduces
# the predicted shear modulus by ~30–40% at high entanglement densities while requiring 2–3x more
# relaxation iterations.

plot_slipping_vs_fixed_modulus(data)
plot_slipping_vs_fixed_iterations(data)
