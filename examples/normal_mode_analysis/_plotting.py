"""
Plotting utilities for Normal Mode Analysis (KPM vs. Exact) examples.
"""

from __future__ import annotations

from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np

_COLORS = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2")


def _plot_kpm_vs_exact_moduli_pair(
    freqs: np.ndarray,
    gp_ex: np.ndarray,
    gp_kpm: np.ndarray,
    gpp_ex: np.ndarray,
    gpp_kpm: np.ndarray,
    labels: Sequence[str],
    title_prefix: str,
) -> None:
    """Shared two-panel log-log comparison of G'(omega) and G''(omega) for Exact vs. KPM."""
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.2, 4.8), sharex=True, sharey=True)

    for idx, label in enumerate(labels):
        c = _COLORS[idx % len(_COLORS)]
        ax1.loglog(
            freqs, gp_ex[idx], "-", color=c, linewidth=2.0, label=f"{label} (Exact)"
        )
        ax1.loglog(
            freqs,
            gp_kpm[idx],
            "--",
            color=c,
            linewidth=2.0,
            alpha=0.85,
            label=f"{label} (KPM)",
        )

        ax2.loglog(
            freqs, gpp_ex[idx], "-", color=c, linewidth=2.0, label=f"{label} (Exact)"
        )
        ax2.loglog(
            freqs,
            gpp_kpm[idx],
            "--",
            color=c,
            linewidth=2.0,
            alpha=0.85,
            label=f"{label} (KPM)",
        )

    for ax, ylabel, subtitle in (
        (
            ax1,
            r"Storage Modulus $G'(\omega)$ [MPa]",
            f"{title_prefix}: Storage Modulus",
        ),
        (ax2,
         r"Loss Modulus $G''(\omega)$ [MPa]",
         f"{title_prefix}: Loss Modulus"),
    ):
        ax.set_xlabel(r"Reduced Frequency $\omega\tau_0$ [-]")
        ax.set_ylabel(ylabel)
        ax.set_title(subtitle)
        ax.set_ylim(bottom=1e-4)
        ax.grid(True, which="both", linestyle="--", alpha=0.45)

    ax1.legend(frameon=True, fontsize=7.8, ncol=2, loc="lower right")
    plt.tight_layout()
    plt.show()


def plot_kpm_vs_exact_chain_length(data: Any) -> None:
    """Plot Exact vs. KPM moduli across precursor chain lengths N."""
    freqs = data["freqs"]
    n_values = data["n_values"]
    n_repeats = int(data["n_repeats"])
    box_l = float(data["box_l_param_studies"])
    labels = [rf"$N = {int(n)}$" for n in n_values]
    _plot_kpm_vs_exact_moduli_pair(
        freqs=freqs,
        gp_ex=data["n_gp_ex_mean"],
        gp_kpm=data["n_gp_kpm_mean"],
        gpp_ex=data["n_gpp_ex_mean"],
        gpp_kpm=data["n_gpp_kpm_mean"],
        labels=labels,
        title_prefix=rf"Chain Length ($L={box_l:.0f}\,\mathrm{{nm}}, p=0.95, r=1.0, n={n_repeats}$)",
    )


def plot_kpm_vs_exact_conversion(data: Any) -> None:
    """Plot Exact vs. KPM moduli across cross-link conversions p."""
    freqs = data["freqs"]
    p_values = data["p_values"]
    n_repeats = int(data["n_repeats"])
    box_l = float(data["box_l_param_studies"])
    labels = [rf"$p = {float(p):.3g}$" for p in p_values]
    _plot_kpm_vs_exact_moduli_pair(
        freqs=freqs,
        gp_ex=data["p_gp_ex_mean"],
        gp_kpm=data["p_gp_kpm_mean"],
        gpp_ex=data["p_gpp_ex_mean"],
        gpp_kpm=data["p_gpp_kpm_mean"],
        labels=labels,
        title_prefix=rf"Conversion ($L={box_l:.0f}\,\mathrm{{nm}}, N=100, r=1.0, n={n_repeats}$)",
    )


def plot_kpm_vs_exact_stoichiometry(data: Any) -> None:
    """Plot Exact vs. KPM moduli across stoichiometric ratios r."""
    freqs = data["freqs"]
    r_values = data["r_values"]
    n_repeats = int(data["n_repeats"])
    box_l = float(data["box_l_param_studies"])
    labels = [rf"$r = {float(r):.2f}$" for r in r_values]
    _plot_kpm_vs_exact_moduli_pair(
        freqs=freqs,
        gp_ex=data["r_gp_ex_mean"],
        gp_kpm=data["r_gp_kpm_mean"],
        gpp_ex=data["r_gpp_ex_mean"],
        gpp_kpm=data["r_gpp_kpm_mean"],
        labels=labels,
        title_prefix=rf"Stoichiometry ($L={box_l:.0f}\,\mathrm{{nm}}, N=100, n={n_repeats}$)",
    )


def plot_kpm_vs_exact_rve_and_scaling(data: Any) -> None:
    """Plot RVE convergence of G'(omega), sample-count precision, and O(N) vs O(N^3) runtime."""
    from scipy import stats

    freqs = data["freqs"]
    box_lengths = data["rve_box_lengths"]
    gp_ex_mean = data["rve_gp_ex_mean"]
    gp_kpm_mean = data["rve_gp_kpm_mean"]
    # shape: (n_boxes, n_rve_repeats, n_freqs)
    gp_kpm_all = data["rve_gp_kpm_all"]
    n_atoms_mean = np.mean(data["rve_n_atoms"], axis=1)
    rve_t_ex = data["rve_t_ex"]
    has_exact = np.any(np.isfinite(rve_t_ex), axis=1)
    t_ex_mean = np.full(len(box_lengths), np.nan)
    t_ex_mean[has_exact] = np.nanmean(rve_t_ex[has_exact], axis=1)
    t_kpm_mean = np.mean(data["rve_t_kpm"], axis=1)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15.5, 4.6))

    # Panel 1: RVE convergence of G'(omega)
    for idx, box_l in enumerate(box_lengths):
        c = _COLORS[idx % len(_COLORS)]
        n_at = int(round(n_atoms_mean[idx]))
        ax1.loglog(
            freqs,
            gp_kpm_mean[idx],
            "--" if np.all(np.isfinite(gp_ex_mean[idx])) else "-",
            color=c,
            linewidth=1.8,
            label=rf"$L={box_l:.0f}\,\mathrm{{nm}}$ ($\sim{n_at:,}$ beads)",
        )
        if np.all(np.isfinite(gp_ex_mean[idx])):
            ax1.loglog(
                freqs,
                gp_ex_mean[idx],
                "-",
                color=c,
                linewidth=1.2,
                alpha=0.55)

    ax1.set_xlabel(r"Reduced Frequency $\omega\tau_0$ [-]")
    ax1.set_ylabel(r"Storage Modulus $G'(\omega)$ [MPa]")
    ax1.set_title("RVE Convergence (Solid: Exact, Dashed: KPM)")
    ax1.set_ylim(bottom=1e-4)
    ax1.grid(True, which="both", linestyle="--", alpha=0.45)
    ax1.legend(frameon=True, fontsize=7.5, loc="lower right")

    # Panel 2: Sample-count 95% CI half-width at omega ~ 1e-2
    # (rubbery/low-freq regime)
    freq_idx = int(np.argmin(np.abs(freqs - 1e-2)))
    n_rve_repeats = gp_kpm_all.shape[1]
    sample_counts = np.arange(2, n_rve_repeats + 1)
    selected_boxes = sorted({0, 2, 4, len(box_lengths) - 1} & set(range(len(box_lengths))))
    for box_idx in selected_boxes:
        c = _COLORS[box_idx % len(_COLORS)]
        box_l = box_lengths[box_idx]
        vals = gp_kpm_all[box_idx, :, freq_idx]
        mean_val = float(np.mean(vals))
        std_val = float(np.std(vals, ddof=1))
        rel_ci_pct = [
            100.0
            * float(stats.t.ppf(0.975, df=n - 1))
            * (std_val / np.sqrt(n))
            / mean_val
            for n in sample_counts
        ]
        ax2.plot(
            sample_counts,
            rel_ci_pct,
            "o-",
            color=c,
            markersize=4,
            linewidth=1.8,
            label=rf"$L={box_l:.0f}\,\mathrm{{nm}}$",
        )

    ax2.set_xlabel(r"Number of Independent Samples $n_{\mathrm{samples}}$")
    ax2.set_ylabel(r"Relative 95% CI Half-Width at $\omega\tau_0=10^{-2}$ [%]")
    ax2.set_title("Ensemble Precision vs. Sample Count")
    ax2.grid(True, linestyle="--", alpha=0.45)
    ax2.legend(frameon=True, fontsize=8.0)

    # Panel 3: Computational scaling O(N^3) vs O(N)
    valid_ex = np.isfinite(t_ex_mean)
    ax3.loglog(
        n_atoms_mean[valid_ex],
        t_ex_mean[valid_ex],
        "s-",
        color="#d62728",
        linewidth=2.0,
        markersize=6,
        label=r"Exact Eigenvalues ($\mathcal{O}(N^3)$)",
    )
    ax3.loglog(
        n_atoms_mean,
        t_kpm_mean,
        "o-",
        color="#1f77b4",
        linewidth=2.0,
        markersize=6,
        label=r"KPM ($\mathcal{O}(N)$)",
    )
    # Reference O(N^3) and O(N) slopes
    n_ref = np.array([n_atoms_mean[0], n_atoms_mean[-1]])
    ax3.loglog(
        n_ref,
        t_ex_mean[0] * (n_ref / n_atoms_mean[0]) ** 3,
        ":",
        color="#d62728",
        alpha=0.5,
        label=r"$\propto N^3$ guide",
    )
    ax3.loglog(
        n_ref,
        t_kpm_mean[0] * (n_ref / n_atoms_mean[0]),
        ":",
        color="#1f77b4",
        alpha=0.5,
        label=r"$\propto N$ guide",
    )

    ax3.set_xlabel(r"Number of Beads $N_{\mathrm{atoms}}$")
    ax3.set_ylabel("Wall-Clock Time per Network [s]")
    ax3.set_title("Computational Scaling: Exact vs. KPM")
    ax3.grid(True, which="both", linestyle="--", alpha=0.45)
    ax3.legend(frameon=True, fontsize=8.0)

    plt.tight_layout()
    plt.show()
