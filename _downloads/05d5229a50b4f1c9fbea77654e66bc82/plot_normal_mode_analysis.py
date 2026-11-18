"""
Normal Mode Analysis (NMA)
==========================

This example demonstrates how to use the NormalModeAnalyzer in pylimer-tools
to predict the loss and storage modulus of polymer networks using exact eigenvalue
solvers or approximate spectral methods (:class:`~pylimer_tools_cpp.KPMConfig`,
:class:`~pylimer_tools_cpp.SLQConfig`, and :class:`~pylimer_tools_cpp.HybridConfig`).

Reference: :cite:t:`gusev_molecular_2024`
"""

from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import (
    HybridConfig,
    KPMConfig,
    NormalModeAnalyzer,
    NormalModeSolver,
    SLQConfig,
)


# Load a LAMMPS data file (replace with your file)
def _find_fixtures_dir() -> Path:
    for cur in (
        [Path(__file__).resolve().parent] if "__file__" in globals() else []
    ) + [Path.cwd()]:
        for p in [cur] + list(cur.parents):
            candidate = p / "tests" / "pylimer_tools" / "fixtures"
            if candidate.is_dir():
                return candidate
    return Path("tests/pylimer_tools/fixtures").resolve()


_fixtures_dir = _find_fixtures_dir()
filePath = str(_fixtures_dir / "structure/network_100_a_46.structure.out")
universe = read_data_file(filePath)

# Create the normal mode analyzer
nma = NormalModeAnalyzer.from_universe(universe)

# 1. Exact analysis with connected component decomposition and
# dense/sparse solver
frequencies = np.logspace(-4, 3, num=int(1e3))

t0 = time.perf_counter()
exact_spectrum = nma.solve_exact(
    solver=NormalModeSolver.AUTO,
    compute_eigenvectors=False,
    decompose_components=True,
)
storage_modulus_exact, loss_modulus_exact = exact_spectrum.moduli(frequencies)
t_exact = time.perf_counter() - t0

# 2. Approximate spectral methods without computing full eigenvalues
# These methods scale linearly O(N) for large networks.
approx_configs: list[tuple[str, KPMConfig | SLQConfig | HybridConfig]] = [
    ("KPM", KPMConfig(num_moments=512, num_random_vectors=20)),
    ("SLQ", SLQConfig(num_steps=200, num_random_vectors=20)),
    (
        "Hybrid",
        HybridConfig(num_steps=200, num_random_vectors=20, num_exact_low_modes=100),
    ),
]

approx_results: dict[str, tuple[np.ndarray, np.ndarray]] = {}
summary_rows: list[dict[str, object]] = [
    {
        "Method": "Exact",
        "Time [s]": f"{t_exact:.3f}",
        "Speedup": "1.0x",
        "Mean Rel. Err. G' [%]": "0.00",
        "Max Rel. Err. G' [%]": "0.00",
        "Mean Rel. Err. G'' [%]": "0.00",
        "Max Rel. Err. G'' [%]": "0.00",
    }
]

for name, cfg in approx_configs:
    t0 = time.perf_counter()
    approx_spectrum = nma.solve_approx(cfg)
    gp_approx, gpp_approx = approx_spectrum.moduli(frequencies)
    dt = time.perf_counter() - t0
    approx_results[name] = (gp_approx, gpp_approx)

    rel_err_gp = (
        np.abs(gp_approx - storage_modulus_exact)
        / np.maximum(storage_modulus_exact, 1e-12)
        * 100.0
    )
    rel_err_gpp = (
        np.abs(gpp_approx - loss_modulus_exact)
        / np.maximum(loss_modulus_exact, 1e-12)
        * 100.0
    )
    summary_rows.append(
        {
            "Method": name,
            "Time [s]": f"{dt:.3f}",
            "Speedup": f"{t_exact / dt:.1f}x",
            "Mean Rel. Err. G' [%]": f"{np.mean(rel_err_gp):.2f}",
            "Max Rel. Err. G' [%]": f"{np.max(rel_err_gp):.2f}",
            "Mean Rel. Err. G'' [%]": f"{np.mean(rel_err_gpp):.2f}",
            "Max Rel. Err. G'' [%]": f"{np.max(rel_err_gpp):.2f}",
        }
    )

# %%
# Summary Table
# -------------
# Compare runtime, speedup, and relative error against exact diagonalization.

summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False))

# %%
# .. todo:
#   Apply conversion factors where necessary
#
# Pairwise Comparison Plots (Exact vs. Each Approximate Method)
# -------------------------------------------------------------
# Each panel compares Exact against a single approximate method so the curves
# remain easy to read.

fig, axes = plt.subplots(
    1,
    len(approx_configs),
    figsize=(15.0, 4.6),
    sharex=True,
    sharey=True,
)

for ax, (name, _) in zip(axes, approx_configs):
    gp_approx, gpp_approx = approx_results[name]
    ax.loglog(frequencies, storage_modulus_exact, "b-", label="G' (Exact)")
    ax.loglog(frequencies, loss_modulus_exact, "r-", label="G'' (Exact)")
    ax.loglog(frequencies, gp_approx, "b--", label=f"G' ({name})")
    ax.loglog(frequencies, gpp_approx, "r--", label=f"G'' ({name})")
    ax.set_title(f"Exact vs. {name}")
    ax.set_xlabel("Frequency [rad/s]")
    ax.grid(True, which="both", linestyle="--", alpha=0.45)
    ax.legend(frameon=True, fontsize=8.5, loc="lower right")

axes[0].set_ylabel("Modulus [a.u.]")
fig.suptitle("Normal Mode Analysis: Exact vs. Approximate Spectral Methods")
plt.tight_layout()
plt.show()
