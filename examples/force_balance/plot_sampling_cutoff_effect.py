r"""
Sampling Cut-Off Effect
=======================

In MEHP Force Balance (:cite:t:`gusev_numerical_2019`, :cite:t:`bernhard_phantom_2025`),
entanglements are sampled between bead pairs within a spatial cut-off
:math:`\|\mathbf{r}_i - \mathbf{r}_j\| \le r_{\text{cut}}`.

Two distinct physical regimes emerge as :math:`r_{\text{cut}}` is varied:

1. **Contact Plateau (:math:`r_{\text{cut}} \lesssim 1.2 b`)**: Only beads in direct physical
   contact are paired, leaving :math:`G` insensitive to small changes in :math:`r_{\text{cut}}`.
2. **Pre-Strain Stiffening (:math:`r_{\text{cut}} \gtrsim 1.5 b`)**: Because pair counts in 3D
   scale as :math:`P(r) \propto r^2`, mean pair separations approach :math:`\bar{r} \approx \frac{3}{4} r_{\text{cut}}`.
   Merging distant beads into links introduces artificial pre-strain and inflates :math:`G`
   (buffered partially when using ``entanglements_as_springs=True``).
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
    plot_cutoff_modulus_vs_rcut,
    plot_cutoff_network_parameter_interaction,
    plot_cutoff_pair_distance_distribution,
    plot_cutoff_residual_vs_rcut,
    plot_cutoff_solver_iterations,
    plot_cutoff_strand_length_interaction,
)
from _study_runners import run_sampling_cutoff_simulation  # noqa: E402
import pylimer_tools_cpp as pt  # noqa: E402

# %%
# Varying ``upper_cutoff`` in Entanglement Sampling
# -------------------------------------------------
# Passing ``upper_cutoff`` to :meth:`~pylimer_tools_cpp.MCUniverseGenerator.sample_entanglements`
# compares rigid links (``EntanglementMode.MERGED_LINKS``) vs. compliant springs
# (``EntanglementMode.SPRING``).

gen = pt.MCUniverseGenerator(18.0, 18.0, 18.0)
gen.set_seed(42)
gen.config_nr_of_mc_steps(0)
gen.add_strands(200, [50] * 200)
gen.add_crosslinkers(100, 4)
gen.link_strands_to_conversion(0.95)

for r_cut in (1.2, 2.6):
    es = gen.sample_entanglements(
        nr_of_entanglements_to_sample=300,
        upper_cutoff=r_cut,
        seed=f"cut_{r_cut}",
    )
    fb = gen.get_mehp_simulator(es, entanglement_mode=pt.EntanglementMode.MERGED_LINKS)
    fb.run()
    print(f"r_cut={r_cut}: sampled={len(es)}, iterations={fb.get_nr_of_iterations()}")

cache_path = get_example_cache_path("sampling_cutoff_effect_cache.npz", "force_balance")
data = load_or_compute_cache(cache_path, run_sampling_cutoff_simulation)

# %%
# Modulus & Initial Force Residual vs. :math:`r_{\text{cut}}`
# -----------------------------------------------------------
# Increasing :math:`r_{\text{cut}}` inflates both the initial force residual (from merging
# separated beads) and the relaxed shear modulus :math:`G`, especially for
# rigid links.

plot_cutoff_modulus_vs_rcut(data)
plot_cutoff_residual_vs_rcut(data)

# %%
# Interactions with Strand Length :math:`N`, Conversion :math:`p`, and Stoichiometry :math:`r`
# ---------------------------------------------------------------------------------------------
# Longer strands (:math:`N = 70`) amplify entanglement sensitivity relative to crosslink
# elasticity, while incomplete conversion (:math:`p = 0.82`) or stoichiometric imbalance
# (:math:`r = 0.75`) reduce the active network modulus uniformly.

plot_cutoff_strand_length_interaction(data)
plot_cutoff_network_parameter_interaction(data)

# %%
# Sampled Pair Distance Distribution & Solver Iterations
# ------------------------------------------------------
# Sampled pair distances follow :math:`P(r) \propto r^2` up to :math:`r_{\text{cut}}`,
# and larger pre-strains require more relaxation iterations to converge.

plot_cutoff_pair_distance_distribution(data)
plot_cutoff_solver_iterations(data)
