r"""
Entanglement Links vs. Springs
==============================

In :class:`~pylimer_tools_cpp.MEHPForceBalance2` (:cite:t:`gusev_numerical_2019`, :cite:t:`bernhard_phantom_2025`),
users can choose between two representations for sampled entanglements via ``entanglements_as_springs``:

1. **Entanglement Links** (``entanglements_as_springs=False``, default):
   The two entangled strand beads are merged into a single shared 4-functional junction node
   (:math:`\mathbf{r}_1 = \mathbf{r}_2`), introducing :math:`N_{\text{ent}}` extra nodes and 0 extra springs.
2. **Entanglement Springs** (``entanglements_as_springs=True``):
   The two beads remain distinct and are coupled by an additional Hookean spring, introducing
   :math:`2 N_{\text{ent}}` extra nodes and :math:`N_{\text{ent}}` extra springs.

In both cases, entanglements in ``MEHPForceBalance2`` **do not slip** along the chain contour
(see :doc:`plot_slipping_vs_fixed_links` for dynamic contour slipping in ``MEHPForceBalance``).
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
from _plotting import plot_link_vs_spring_dof, plot_link_vs_spring_modulus  # noqa: E402
from _study_runners import run_entanglement_link_vs_spring_simulation  # noqa: E402
import pylimer_tools_cpp as pt  # noqa: E402

# %%
# Comparing Links vs. Springs on a Generated Network
# --------------------------------------------------
# Passing the same sampled entanglements to :meth:`~pylimer_tools_cpp.MCUniverseGenerator.get_force_balance2`
# with ``entanglements_as_springs=False`` vs. ``True`` toggles between rigid merged junctions and
# compliant Hookean entanglement springs.

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

fb_links = gen.get_force_balance2(es, entanglements_as_springs=False)
fb_links.run_force_relaxation()

fb_springs = gen.get_force_balance2(es, entanglements_as_springs=True)
fb_springs.run_force_relaxation()

print(
    f"Extra atoms — links: {fb_links.get_nr_of_extra_atoms()}, "
    f"springs: {fb_springs.get_nr_of_extra_atoms()}"
)

cache_path = get_example_cache_path(
    "entanglement_link_vs_spring_cache.npz", "force_balance"
)
data = load_or_compute_cache(cache_path, run_entanglement_link_vs_spring_simulation)

# %%
# Equilibrium Shear Modulus & System Degrees of Freedom
# -----------------------------------------------------
# Multi-seed ensemble results (95% Student's t confidence intervals):
#
# - **Shear Modulus**: Merged entanglement links stiffen nearly linearly with entanglement density,
#   whereas entanglement springs saturate because the intermediate spring adds series compliance.
# - **System Size**: Entanglement springs add twice as many extra nodes to the linear system.

plot_link_vs_spring_modulus(data)
plot_link_vs_spring_dof(data)
