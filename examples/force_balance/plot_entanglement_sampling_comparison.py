#!/usr/bin/env python
"""
Comparing Slip-Link Sampling Methods: Discrete vs. Analytical
=============================================================

In polymer network elasticity simulations using the Maximum Entropy Homogenization
Procedure (MEHP) and Force Balance (:cite:t:`gusev_numerical_2019`, :cite:t:`bernhard_phantom_2025`),
entanglements are represented as slip-links or localized crosslink constraints.

Historically, sampling entanglements required instantiating the complete bead-spring
:class:`~pylimer_tools_cpp.Universe` (:meth:`~pylimer_tools_cpp.MCUniverseGenerator.get_universe`)
with 3D coordinates for every single bead along every strand. While exact, this creates
an enormous memory bottleneck: networks with hundreds of thousands or millions of beads
require gigabytes of RAM just to store atom tables, neighbor lists, and graph structures
that are immediately discarded once the network is mapped into Hookean springs.

To eliminate this memory wall, ``pylimer-tools`` introduces an **analytical Brownian-bridge
sampling method** (:meth:`~pylimer_tools_cpp.MCUniverseGenerator.sample_entanglements` and
:func:`~pylimer_tools_cpp.randomly_sample_entanglements_analytical`). Operating directly on
the lightweight :class:`~pylimer_tools_cpp.CrosslinkerUniverse` without generating intermediate
beads, it evaluates conditional Gaussian distributions between crosslinker endpoints and
constructs :class:`~pylimer_tools_cpp.MEHPForceBalance` or :class:`~pylimer_tools_cpp.MEHPForceBalance2`
networks directly.

This example illustrates:

1. **Memory & Runtime Scaling**: The order-of-magnitude memory reduction achieved by bypassing full bead generation.
2. **Conformational Distributions**: Why analytical Brownian-bridge sampling faithfully replicates the discrete spatial statistics along strands.
3. **Network Topology & Modulus across Conversion**: How topological filtering of dangling/soluble fractions influences network properties and shear modulus.
4. **Practical Decision Guide**: When to choose which method, key trade-offs, and physical assumptions.
"""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

# %%
# Background: The Memory Wall in Slip-Link Sampling
# -------------------------------------------------
# In a polymer network, strands between crosslinkers fluctuate thermally.
# Under the affine network model and Gaussian chain approximation, the position
# :math:`\\mathbf{r}(\\alpha)` of a bead at contour fraction :math:`\\alpha \\in (0, 1)`
# on a strand connecting crosslinker anchors :math:`\\mathbf{R}_1` and :math:`\\mathbf{R}_2`
# follows a **Brownian bridge** conditional on its endpoints:
#
# .. math::
#
#    \\mathbf{r}(\\alpha) \\sim \\mathcal{N}\\left( (1-\\alpha)\\mathbf{R}_1 + \\alpha\\mathbf{R}_2,\\;
#    \\alpha(1-\\alpha) N b^2 / 3 \\,\\mathbf{I} \\right)
#
# Where:
#
# - :math:`N` is the number of bonds in the strand,
# - :math:`b` is the statistical segment length (Kuhn length or bead distance),
# - :math:`\\alpha = k / N` is the fractional contour position of bead :math:`k`.
#
# The traditional **Discrete Method** explicitly generates all :math:`N_{\\text{atoms}}` beads in
# 3D space and constructs spatial cell lists to find pairs within cutoff :math:`r_{\\text{cut}}`.
#
# The new **Analytical Method** skips bead coordinates entirely. It works directly on
# :class:`~pylimer_tools_cpp.CrosslinkerUniverse` (:math:`N_{\\text{xlinks}} \\ll N_{\\text{atoms}}`),
# selecting pairs of strands and sampling continuous contour positions :math:`\\alpha_1, \\alpha_2`
# from their conditional Gaussian distributions. Furthermore, it incorporates **topological leaf
# peeling** (:paramref:`filter_dangling_and_soluble=True`), guaranteeing that slip-links are
# placed only on elastically active, load-bearing strands.

# %%
# Loading Benchmark & Statistical Data
# ------------------------------------
# To comply with fast execution and CI limits while providing statistically
# significant comparisons on large systems (up to 200,000+ beads), precomputed
# benchmark data is loaded from the committed cache.
#
# The dataset can be recomputed from scratch at any time by passing ``--recompute``
# to this script or running ``data/generate_cache.py``.


def _get_cache_path() -> Path:
    base = (
        Path(__file__).resolve().parent
        if "__file__" in globals()
        else Path.cwd() / "examples" / "force_balance"
    )
    return base / "data" / "entanglement_sampling_comparison_cache.npz"


cache_path = _get_cache_path()

# Support CLI argument --recompute if run as main script
parser = argparse.ArgumentParser(
    description="Compare entanglement sampling methods")
parser.add_argument(
    "--recompute", action="store_true", help="Recompute benchmarks from scratch"
)
args, _ = parser.parse_known_args(
    sys.argv[1:] if "__file__" in globals() else [])

if args.recompute or not cache_path.exists():
    print(f"Recomputing benchmark data and saving to {cache_path}...")
    from data.generate_cache import main as recompute_cache  # type: ignore

    recompute_cache()

data = np.load(cache_path)

scaling_atoms = data["scaling_atoms"]
scaling_time_disc = data["scaling_time_discrete"]
scaling_rss_disc = data["scaling_rss_discrete"]
scaling_time_ana = data["scaling_time_analytical"]
scaling_rss_ana = data["scaling_rss_analytical"]

alpha_ana = data["alpha_analytical"]
dist_ana = data["dist_analytical"]
alpha_disc = data["alpha_discrete"]
dist_disc = data["dist_discrete"]

conversions = data["conversions"]
g_disc = data["g_discrete"]
w_act_disc = data["w_active_discrete"]
g_ana_filt = data["g_ana_filtered"]
w_act_ana_filt = data["w_active_ana_filtered"]
n_ent_ana_filt = data["n_ent_ana_filtered"]
g_ana_unfilt = data["g_ana_unfiltered"]
w_act_ana_unfilt = data["w_active_ana_unfiltered"]
n_ent_ana_unfilt = data["n_ent_ana_unfiltered"]
g_ana_springs = data["g_ana_springs"]
target_n_ent = int(data["target_n_ent"])

# %%
# Comparison 1: Memory Footprint & Runtime Scaling
# ------------------------------------------------
# As the network size increases, explicit bead generation rapidly inflates peak resident
# memory (RSS) due to bead coordinate arrays, neighbor search structures, and igraph graphs.
# In contrast, analytical sampling bypasses intermediate bead allocation, keeping memory
# footprint virtually flat and reducing peak RAM by 80% to over 90% on
# large networks.

fig, (ax_mem, ax_time) = plt.subplots(1, 2, figsize=(11, 4.5))

# Peak Memory
ax_mem.plot(
    scaling_atoms,
    scaling_rss_disc,
    "o-",
    color="#d95f02",
    linewidth=2,
    markersize=6,
    label="Discrete (Full Universe)",
)
ax_mem.plot(
    scaling_atoms,
    scaling_rss_ana,
    "s-",
    color="#1b9e77",
    linewidth=2,
    markersize=6,
    label="Analytical (Universe-Free)",
)

# Annotate memory savings on the largest system
max_atoms = scaling_atoms[-1]
mem_disc_last = scaling_rss_disc[-1]
mem_ana_last = scaling_rss_ana[-1]
pct_savings = (1.0 - mem_ana_last / mem_disc_last) * 100.0

ax_mem.annotate(
    f"{pct_savings:.0f}% Memory\nReduction",
    xy=(max_atoms, mem_ana_last),
    xytext=(max_atoms * 0.45, (mem_disc_last + mem_ana_last) / 2.0),
    arrowprops=dict(arrowstyle="->", color="#1b9e77", lw=1.5),
    fontsize=10,
    fontweight="bold",
    color="#1b9e77",
    bbox=dict(boxstyle="round,pad=0.3", fc="#e6f5f0", ec="#1b9e77", lw=1),
)

ax_mem.set_xscale("log")
ax_mem.set_yscale("log")
ax_mem.set_xlabel(r"Total Atoms in Network ($N_{\mathrm{atoms}}$)")
ax_mem.set_ylabel("Peak Resident Memory (RSS) [MB]")
ax_mem.set_title("Peak Memory Consumption")
ax_mem.grid(True, which="both", linestyle="--", alpha=0.5)
ax_mem.legend(frameon=True)

# Runtime
ax_time.plot(
    scaling_atoms,
    scaling_time_disc,
    "o-",
    color="#d95f02",
    linewidth=2,
    markersize=6,
    label="Discrete (Full Universe)",
)
ax_time.plot(
    scaling_atoms,
    scaling_time_ana,
    "s-",
    color="#1b9e77",
    linewidth=2,
    markersize=6,
    label="Analytical (Universe-Free)",
)

ax_time.set_xlabel(r"Total Atoms in Network ($N_{\mathrm{atoms}}$)")
ax_time.set_ylabel("Setup & Sampling Time [s]")
ax_time.set_title("Setup & Sampling Execution Time")
ax_time.grid(True, linestyle="--", alpha=0.5)
ax_time.legend(frameon=True)

plt.tight_layout()
plt.show()

# %%
# Comparison 2: Conformational Statistics (Contour & Spatial Distributions)
# -------------------------------------------------------------------------
# A natural question is: *Does the analytical Brownian-bridge method reproduce the
# same physical distribution of entanglements as explicit spatial bead sampling?*
#
# Below, we compare:
#
# 1. **Contour position** :math:`\\alpha \\in (0, 1)`: Both distributions exhibit the characteristic
#    bell curve peaked near :math:`\\alpha \\approx 0.5`. Beads near crosslinkers are spatially
#    pinned, whereas strand midpoints explore the largest conformational volume and are most
#    likely to encounter neighboring strands.
# 2. **Pair distance** :math:`r \\in [0, r_{\\text{cut}}]`: Both methods strictly adhere to the
# sampling cutoff :math:`r_{\\text{cut}}`, producing identical Euclidean
# distance distributions.

fig, (ax_alpha, ax_dist) = plt.subplots(1, 2, figsize=(11, 4.5))

# Contour fraction distribution
bins_alpha = np.linspace(0.0, 1.0, 35)
ax_alpha.hist(
    alpha_disc,
    bins=bins_alpha,
    density=True,
    alpha=0.5,
    color="#d95f02",
    label="Discrete Beads",
    edgecolor="white",
)
ax_alpha.hist(
    alpha_ana,
    bins=bins_alpha,
    density=True,
    alpha=0.5,
    color="#1b9e77",
    label="Analytical Bridge",
    edgecolor="white",
)

ax_alpha.set_xlabel(r"Contour Position $\alpha \in (0, 1)$ along Strand")
ax_alpha.set_ylabel("Probability Density")
ax_alpha.set_title("Entanglement Position along Strands")
ax_alpha.grid(True, linestyle="--", alpha=0.5)
ax_alpha.legend(frameon=True)

# Pair distance distribution
bins_dist = np.linspace(0.0, float(np.max(dist_disc)) * 1.02, 35)
ax_dist.hist(
    dist_disc,
    bins=bins_dist,
    density=True,
    alpha=0.5,
    color="#d95f02",
    label="Discrete Pairs",
    edgecolor="white",
)
ax_dist.hist(
    dist_ana,
    bins=bins_dist,
    density=True,
    alpha=0.5,
    color="#1b9e77",
    label="Analytical Pairs",
    edgecolor="white",
)

ax_dist.set_xlabel(r"Pair Separation Distance $r$ [code units]")
ax_dist.set_ylabel("Probability Density")
ax_dist.set_title(r"Entanglement Pair Separation ($r \leq r_{\mathrm{cut}}$)")
ax_dist.grid(True, linestyle="--", alpha=0.5)
ax_dist.legend(frameon=True)

plt.tight_layout()
plt.show()

# %%
# Comparison 3: Crosslinker Conversion & Topological Filtering
# ------------------------------------------------------------
# In real polymer networks, crosslinking reactions do not reach 100% conversion (:math:`p < 1.0`).
# This introduces topological defects: **dangling ends** (chains connected to only one crosslinker)
# and **soluble fractions** (unconnected clusters).
#
# When entanglements are sampled on defective networks, two strategies are possible:
#
# - **Unfiltered Sampling** (Discrete or Analytical with :paramref:`filter_dangling_and_soluble=False`):
#   Entanglements are placed blindly across all strands. Slip-links that land on dangling ends or
#   soluble molecules cannot bear macroscopic load and are mechanically inactive.
# - **Topologically Filtered Sampling** (Analytical with :paramref:`filter_dangling_and_soluble=True`, default):
#   Iterative leaf peeling prunes dangling chains and sol before sampling, placing all slip-links
#   exclusively onto elastically active strands.
#
# We also compare modeling slip-links as **merged crosslink nodes** (:paramref:`entanglements_as_springs=False`)
# versus **Hookean springs** (:paramref:`entanglements_as_springs=True`).

fig, (ax_mod, ax_frac) = plt.subplots(1, 2, figsize=(11, 4.5))

# Shear modulus vs conversion
ax_mod.plot(
    conversions,
    g_disc,
    "o--",
    color="#d95f02",
    linewidth=1.8,
    label="Discrete (Unfiltered)",
)
ax_mod.plot(
    conversions,
    g_ana_unfilt,
    "^:",
    color="#7570b3",
    linewidth=1.8,
    label="Analytical (Unfiltered)",
)
ax_mod.plot(
    conversions,
    g_ana_filt,
    "s-",
    color="#1b9e77",
    linewidth=2.2,
    label="Analytical (Filtered, Default)",
)
ax_mod.plot(
    conversions,
    g_ana_springs,
    "d-.",
    color="#e7298a",
    linewidth=1.8,
    label="Analytical (Filtered, Springs)",
)

ax_mod.set_xlabel("Crosslinker Conversion $p$")
ax_mod.set_ylabel("Shear Modulus $G$ [MPa]")
ax_mod.set_title("Equilibrium Shear Modulus vs. Conversion")
ax_mod.grid(True, linestyle="--", alpha=0.5)
ax_mod.legend(frameon=True, fontsize=9)

# Active entanglement fraction vs conversion
frac_active_filt = n_ent_ana_filt / target_n_ent
frac_active_unfilt = w_act_ana_unfilt  # proportion on active backbone

ax_frac.plot(
    conversions,
    w_act_disc,
    "o--",
    color="#d95f02",
    linewidth=1.8,
    label=r"Active Weight Fraction $w_{\mathrm{active}}$",
)
ax_frac.plot(
    conversions,
    frac_active_filt,
    "s-",
    color="#1b9e77",
    linewidth=2.2,
    label="Filtered Active Entanglement Fraction",
)

ax_frac.set_xlabel("Crosslinker Conversion $p$")
ax_frac.set_ylabel("Fraction")
ax_frac.set_title("Network Activity & Effective Entanglements")
ax_frac.grid(True, linestyle="--", alpha=0.5)
ax_frac.legend(frameon=True, fontsize=9)

plt.tight_layout()
plt.show()

# %%
# Method Comparison & Practical Guidance
# --------------------------------------
# The table below summarizes the key trade-offs between discrete and analytical sampling:
#
# .. list-table::
#    :widths: 22 38 40
#    :header-rows: 1
#
#    * - Criterion
#      - Discrete (Full Universe)
#      - Analytical (Universe-Free)
#    * - **Memory Footprint**
#      - High: :math:`O(N_{\\text{atoms}})`, reaches gigabytes for large networks.
#      - Minimal: :math:`O(N_{\\text{xlinks}} + N_{\\text{strands}})`, ~80–95% lower.
#    * - **Execution Speed**
#      - Slower: requires placing every bead and building 3D neighbor cell lists.
#      - Fast: samples conditional Gaussian distributions between crosslinkers directly.
#    * - **Strand Statistics**
#      - Exact bead coordinates as generated by the Monte Carlo universe generator.
#      - Assumes conditional Gaussian random walks (Brownian bridge) between anchors.
#    * - **Topological Filtering**
#      - Samples indiscriminately across all beads (unless manually cleaned).
#      - Built-in iterative leaf peeling (:paramref:`filter_dangling_and_soluble=True`).
#    * - **Ideal Use Case**
#      - Downstream MD/DPD simulations requiring full 3D atom configurations.
#      - MEHP Force Balance runs, high-throughput sweeps, and large networks.
#
# Drawbacks and Limitations to Keep in Mind
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# 1. **Gaussian Chain Assumption**:
#    The analytical method assumes chain segments follow ideal Gaussian Brownian bridges
#    between crosslinkers. For polymer melts under severe geometric confinement, highly
#    stretched chains near finite extensibility (Langevin regime), or semi-flexible polymers
#    with high persistence lengths, actual bead paths may deviate from ideal Brownian bridges.
#
# 2. **Effective Entanglement Density in Defective Networks**:
#    When setting a target number of entanglements :math:`N_{\\text{ent}}`:
#
#    - With :paramref:`filter_dangling_and_soluble=True` (default), all :math:`N_{\\text{ent}}`
#      links are concentrated onto the load-bearing gel fraction.
#    - With :paramref:`filter_dangling_and_soluble=False` or discrete sampling, some links
#      fall on dangling chains or sol complexes.
#    - If matching physical experiments where nominal entanglement density is measured relative
#      to total volume rather than the gel fraction, users can pass
#      ``filter_dangling_and_soluble=False`` or scale the target :math:`N_{\\text{ent}}` by
#      the active weight fraction.
