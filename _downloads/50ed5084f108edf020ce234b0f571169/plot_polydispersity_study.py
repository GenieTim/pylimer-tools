#!/usr/bin/env python
r"""
Polydispersity Study
====================

In this example, we study the effect of precursor strand length polydispersity
on the equilibrium shear modulus and topological network structure of end-linked
polymer networks (:cite:t:`tsimouri_monte_2021`, :cite:t:`bernhard_phantom_2025`).

Strand lengths are sampled from a log-normal distribution parameterized by the
polydispersity index :math:`D = M_w / M_n` at fixed mean strand length :math:`\langle N \rangle = 50`:

.. math::

   s = \sqrt{\ln(D)}, \quad m = \ln(\langle N \rangle) - \frac{1}{2} s^2

We evaluate both the total entangled modulus :math:`G_{\text{total}}` and the phantom
modulus :math:`G_{\text{phantom}}` (without entanglements), isolating the entanglement
contribution :math:`G_{\text{ent}} = G_{\text{total}} - G_{\text{phantom}}`.
"""

import math
from pathlib import Path
import random
import sys

if "__file__" in globals():
    _examples_dir = str(Path(__file__).resolve().parents[1])
    _fb_dir = str(Path(__file__).resolve().parent)
    for _p in (_examples_dir, _fb_dir):
        if _p not in sys.path:
            sys.path.insert(0, _p)

from _example_cache import get_example_cache_path, load_or_compute_cache  # noqa: E402
from _plotting import (  # noqa: E402
    plot_polydispersity_fractions,
    plot_polydispersity_modulus,
)
from _study_runners import run_polydispersity_simulation  # noqa: E402
import pylimer_tools_cpp as pt  # noqa: E402

# %%
# Generating a Polydisperse Network
# ---------------------------------
# Pass a per-strand list of lengths to :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_strands`
# and compare :meth:`~pylimer_tools_cpp.MCUniverseGenerator.get_force_balance2` with and
# without sampled entanglements.

rng = random.Random(42)
d, avg_len, n_strands = 1.5, 50, 200
s = math.sqrt(math.log(d))
m = math.log(avg_len) - 0.5 * s**2
strand_lengths = [
    max(1, int(round(rng.lognormvariate(m, s)))) for _ in range(n_strands)
]

gen = pt.MCUniverseGenerator(18.0, 18.0, 18.0)
gen.set_seed(42)
gen.config_nr_of_mc_steps(0)
gen.add_strands(nr_of_strands=n_strands, strand_lengths=strand_lengths)
gen.add_crosslinkers(int(0.5 * n_strands), 4)
gen.link_strands_to_conversion(crosslinker_conversion=0.925)

fb_phantom = gen.get_force_balance2()
fb_phantom.run_force_relaxation()

es = gen.sample_entanglements(300, upper_cutoff=2.6, seed="demo")
fb_total = gen.get_force_balance2(es)
fb_total.run_force_relaxation()

print(
    f"Active fraction: {fb_total.get_active_weight_fraction():.3f}, "
    f"Dangling fraction: {fb_total.get_dangling_weight_fraction():.3f}"
)

cache_path = get_example_cache_path(
    "polydispersity_study_cache.npz",
    "force_balance")
data = load_or_compute_cache(cache_path, run_polydispersity_simulation)

# %%
# Shear Modulus & Topological Fractions vs. Polydispersity
# --------------------------------------------------------
# Across moderate polydispersities (:math:`1.0 \le D \le 1.8`, :math:`N_{\text{strands}} = 10\,000`,
# 95% confidence intervals over 8 repeats), :math:`G_{\text{total}}` and :math:`w_{\text{active}}`
# remain remarkably stable because :math:`\langle N \rangle` is held constant.

plot_polydispersity_modulus(data)
plot_polydispersity_fractions(data)
