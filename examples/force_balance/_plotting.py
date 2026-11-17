"""
Matplotlib plotting helpers for the force_balance Sphinx-Gallery examples.

This module is prefixed with an underscore so Sphinx-Gallery ignores it when
building documentation, allowing each example section to pair concise
``pylimer_tools`` API usage directly with its resulting figure.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


def _format_k(x: float, _pos: Any) -> str:
    if x >= 1e6:
        return f"{x * 1e-6:.1f}".rstrip("0").rstrip(".") + "M"
    if x >= 1e3:
        return f"{int(x * 1e-3)}k"
    return str(int(x))


_K_FORMATTER = ticker.FuncFormatter(_format_k)


def _plot_with_iqr(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    iqr: Tuple[Optional[np.ndarray], Optional[np.ndarray]],
    fmt: str,
    color: str,
    label: str,
) -> None:
    q25, q75 = iqr
    if q25 is not None and q75 is not None:
        yerr = np.vstack([np.maximum(0.0, y - q25), np.maximum(0.0, q75 - y)])
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=fmt,
            color=color,
            linewidth=2,
            markersize=6,
            capsize=3,
            elinewidth=1.2,
            label=label,
        )
        ax.fill_between(x, q25, q75, color=color, alpha=0.14)
    else:
        ax.plot(x, y, fmt, color=color, linewidth=2, markersize=6, label=label)


def plot_link_vs_spring_modulus(data: Any) -> None:
    density_factors = data["density_factors"]
    n_repeats = int(data.get("n_repeats", 5))
    g_links_mean = data["g_links_mean"]
    g_links_ci = data["g_links_ci"]
    g_springs_mean = data["g_springs_mean"]
    g_springs_ci = data["g_springs_ci"]
    fig, ax_mod = plt.subplots(figsize=(7, 4.5))
    ax_mod.plot(
        density_factors,
        g_links_mean,
        "s-",
        color="#d95f02",
        linewidth=2,
        markersize=6,
        label=f"Entanglement Links (Merged Nodes, N={n_repeats})",
    )
    ax_mod.fill_between(
        density_factors,
        g_links_mean - g_links_ci,
        g_links_mean + g_links_ci,
        color="#d95f02",
        alpha=0.2,
        label="95% Confidence Interval",
    )
    ax_mod.errorbar(
        density_factors,
        g_links_mean,
        yerr=g_links_ci,
        fmt="none",
        ecolor="#d95f02",
        capsize=4,
        elinewidth=1.2,
    )

    ax_mod.plot(
        density_factors,
        g_springs_mean,
        "o-",
        color="#1b9e77",
        linewidth=2,
        markersize=6,
        label=f"Entanglement Springs (Hookean Springs, N={n_repeats})",
    )
    ax_mod.fill_between(
        density_factors,
        g_springs_mean - g_springs_ci,
        g_springs_mean + g_springs_ci,
        color="#1b9e77",
        alpha=0.2,
    )
    ax_mod.errorbar(
        density_factors,
        g_springs_mean,
        yerr=g_springs_ci,
        fmt="none",
        ecolor="#1b9e77",
        capsize=4,
        elinewidth=1.2,
    )

    ax_mod.set_xlabel(
        r"Entanglement Density Factor ($N_{\mathrm{ent}} / N_{\mathrm{nominal}}$)"
    )
    ax_mod.set_ylabel(r"Shear Modulus $G$ [Pa]")
    ax_mod.set_title("Equilibrium Shear Modulus vs. Entanglement Density")
    ax_mod.grid(True, linestyle="--", alpha=0.5)
    ax_mod.legend(frameon=True)
    plt.tight_layout()
    plt.show()


def plot_link_vs_spring_dof(data: Any) -> None:
    density_factors = data["density_factors"]
    extra_links_mean = data["extra_atoms_links_mean"]
    extra_springs_mean = data["extra_atoms_springs_mean"]
    fig, ax_dof = plt.subplots(figsize=(7, 4.5))
    ax_dof.plot(
        density_factors,
        extra_links_mean,
        "s-",
        color="#d95f02",
        linewidth=2,
        markersize=6,
        label=r"Entanglement Links: $N_{\mathrm{extra}} = N_{\mathrm{ent}}$",
    )
    ax_dof.plot(
        density_factors,
        extra_springs_mean,
        "o-",
        color="#1b9e77",
        linewidth=2,
        markersize=6,
        label=r"Entanglement Springs: $N_{\mathrm{extra}} = 2 N_{\mathrm{ent}}$",
    )

    ax_dof.set_xlabel(
        r"Entanglement Density Factor ($N_{\mathrm{ent}} / N_{\mathrm{nominal}}$)"
    )
    ax_dof.set_ylabel("Extra Nodes Added to SLE")
    ax_dof.set_title("System Size / Extra Degrees of Freedom")
    ax_dof.grid(True, linestyle="--", alpha=0.5)
    ax_dof.legend(frameon=True)

    plt.tight_layout()
    plt.show()


def plot_slipping_vs_fixed_modulus(data: Any) -> None:
    density_factors = data["density_factors"]
    n_repeats = int(data.get("n_repeats", 8))
    g_fixed_mean = data["g_fixed_mean"]
    g_fixed_ci = data["g_fixed_ci"]
    g_slip_mean = data["g_slip_mean"]
    g_slip_ci = data["g_slip_ci"]
    fig, ax_mod = plt.subplots(figsize=(7, 4.5))
    ax_mod.plot(
        density_factors,
        g_fixed_mean,
        "s-",
        color="#2b83ba",
        linewidth=2,
        markersize=6,
        label=f"Fixed Links (disable_slipping=True, N={n_repeats})",
    )
    ax_mod.fill_between(
        density_factors,
        g_fixed_mean - g_fixed_ci,
        g_fixed_mean + g_fixed_ci,
        color="#2b83ba",
        alpha=0.2,
        label="95% Confidence Interval",
    )
    ax_mod.errorbar(
        density_factors,
        g_fixed_mean,
        yerr=g_fixed_ci,
        fmt="none",
        ecolor="#2b83ba",
        capsize=4,
        elinewidth=1.2,
    )

    ax_mod.plot(
        density_factors,
        g_slip_mean,
        "o--",
        color="#d7191c",
        linewidth=2,
        markersize=6,
        label=f"Slipping Links (disable_slipping=False, N={n_repeats})",
    )
    ax_mod.fill_between(
        density_factors,
        g_slip_mean - g_slip_ci,
        g_slip_mean + g_slip_ci,
        color="#d7191c",
        alpha=0.2,
    )
    ax_mod.errorbar(
        density_factors,
        g_slip_mean,
        yerr=g_slip_ci,
        fmt="none",
        ecolor="#d7191c",
        capsize=4,
        elinewidth=1.2,
    )

    ax_mod.set_xlabel(
        r"Entanglement Density Factor ($N_{\mathrm{ent}} / N_{\mathrm{nominal}}$)"
    )
    ax_mod.set_ylabel(r"Shear Modulus $G$ [Pa]")
    ax_mod.set_title("Equilibrium Shear Modulus vs. Entanglement Density")
    ax_mod.grid(True, linestyle="--", alpha=0.5)
    ax_mod.legend(frameon=True)
    plt.tight_layout()
    plt.show()


def plot_slipping_vs_fixed_iterations(data: Any) -> None:
    density_factors = data["density_factors"]
    iter_fixed_mean = data["iter_fixed_mean"]
    iter_fixed_ci = data["iter_fixed_ci"]
    iter_slip_mean = data["iter_slip_mean"]
    iter_slip_ci = data["iter_slip_ci"]
    fig, ax_iter = plt.subplots(figsize=(7, 4.5))
    ax_iter.plot(
        density_factors,
        iter_fixed_mean,
        "s-",
        color="#2b83ba",
        linewidth=2,
        markersize=6,
        label=r"Fixed Links (disable_slipping=True)",
    )
    ax_iter.fill_between(
        density_factors,
        iter_fixed_mean - iter_fixed_ci,
        iter_fixed_mean + iter_fixed_ci,
        color="#2b83ba",
        alpha=0.2,
    )

    ax_iter.plot(
        density_factors,
        iter_slip_mean,
        "o--",
        color="#d7191c",
        linewidth=2,
        markersize=6,
        label=r"Slipping Links (disable_slipping=False)",
    )
    ax_iter.fill_between(
        density_factors,
        iter_slip_mean - iter_slip_ci,
        iter_slip_mean + iter_slip_ci,
        color="#d7191c",
        alpha=0.2,
    )

    ax_iter.set_xlabel(
        r"Entanglement Density Factor ($N_{\mathrm{ent}} / N_{\mathrm{nominal}}$)"
    )
    ax_iter.set_ylabel("Iterations to Convergence")
    ax_iter.set_title("Convergence Iterations vs. Entanglement Density")
    ax_iter.grid(True, linestyle="--", alpha=0.5)
    ax_iter.legend(frameon=True)

    plt.tight_layout()
    plt.show()


def plot_polydispersity_modulus(data: Any) -> None:
    pdi_values = data["pdi_values"]
    n_repeats = int(data.get("n_repeats", 8))
    n_strands = int(data.get("n_strands", 10000))
    g_total_mean = data["g_total_mean"]
    g_total_ci = data["g_total_ci"]
    g_phantom_mean = data["g_phantom_mean"]
    g_phantom_ci = data["g_phantom_ci"]
    g_ent_mean = data["g_ent_mean"]
    g_ent_ci = data["g_ent_ci"]
    fig, ax1 = plt.subplots(figsize=(7, 4.5))
    ax1.plot(
        pdi_values,
        g_total_mean,
        "o-",
        color="#1f77b4",
        lw=2,
        label=r"$G_{\mathrm{total}}$ (Entangled)",
    )
    ax1.fill_between(
        pdi_values,
        g_total_mean - g_total_ci,
        g_total_mean + g_total_ci,
        color="#1f77b4",
        alpha=0.25,
        label=r"95% CI ($G_{\mathrm{total}}$)",
    )

    ax1.plot(
        pdi_values,
        g_ent_mean,
        "d--",
        color="#2ca02c",
        lw=1.8,
        label=r"$G_{\mathrm{ent}}$ (Entanglements)",
    )
    ax1.fill_between(
        pdi_values,
        g_ent_mean - g_ent_ci,
        g_ent_mean + g_ent_ci,
        color="#2ca02c",
        alpha=0.2,
    )

    ax1.plot(
        pdi_values,
        g_phantom_mean,
        "s-.",
        color="#d62728",
        lw=1.8,
        label=r"$G_{\mathrm{phantom}}$ (Phantom)",
    )
    ax1.fill_between(
        pdi_values,
        g_phantom_mean - g_phantom_ci,
        g_phantom_mean + g_phantom_ci,
        color="#d62728",
        alpha=0.2,
    )

    ax1.set_xlabel("Polydispersity Index ($D = M_w / M_n$)", fontsize=11)
    ax1.set_ylabel("Shear Modulus [MPa]", fontsize=11)
    ax1.set_title(
        f"Elastic Shear Modulus vs. Polydispersity\n"
        f"($N_{{\\mathrm{{strands}}}} = {n_strands:,}$, $N_{{\\mathrm{{repeats}}}} = {
            n_repeats
        }$)",
        fontsize=12,
    )
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="best", framealpha=0.9)
    plt.tight_layout()
    plt.show()


def plot_polydispersity_fractions(data: Any) -> None:
    pdi_values = data["pdi_values"]
    w_active_mean = data["w_active_mean"]
    w_active_ci = data["w_active_ci"]
    w_dangling_mean = data["w_dangling_mean"]
    w_dangling_ci = data["w_dangling_ci"]
    fig, ax2 = plt.subplots(figsize=(7, 4.5))
    ax2.plot(
        pdi_values,
        w_active_mean,
        "o-",
        color="#4c72b0",
        lw=2,
        label=r"$w_{\mathrm{active}}$ (Elastically Active)",
    )
    ax2.fill_between(
        pdi_values,
        w_active_mean - w_active_ci,
        w_active_mean + w_active_ci,
        color="#4c72b0",
        alpha=0.25,
    )

    ax2.plot(
        pdi_values,
        w_dangling_mean,
        "s--",
        color="#c44e52",
        lw=1.8,
        label=r"$w_{\mathrm{dangling}}$ (Dangling Ends)",
    )
    ax2.fill_between(
        pdi_values,
        w_dangling_mean - w_dangling_ci,
        w_dangling_mean + w_dangling_ci,
        color="#c44e52",
        alpha=0.2,
    )

    ax2.set_xlabel("Polydispersity Index ($D = M_w / M_n$)", fontsize=11)
    ax2.set_ylabel("Weight Fraction", fontsize=11)
    ax2.set_title(
        "Network Topology vs. Polydispersity\n(Conversion $p = 92.5\\%$, $f = 4$)",
        fontsize=12,
    )
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="best", framealpha=0.9)

    plt.tight_layout()
    plt.show()


def plot_cutoff_modulus_vs_rcut(data: Any) -> None:
    cutoffs = data["cutoffs"]
    n_repeats = int(data.get("n_repeats", 8))
    b = float(data.get("b", 0.969))
    default_cutoff = float(data.get("default_cutoff", 2.594))
    g_l_mean = data["g_links_mean"]
    g_l_ci = data["g_links_ci"]
    g_s_mean = data["g_springs_mean"]
    g_s_ci = data["g_springs_ci"]
    fig, ax_mod = plt.subplots(figsize=(7, 4.5))
    ax_mod.plot(
        cutoffs,
        g_l_mean,
        "s-",
        color="#d95f02",
        linewidth=2,
        label=f"Entanglement Links (N={n_repeats})",
    )
    ax_mod.fill_between(
        cutoffs,
        g_l_mean - g_l_ci,
        g_l_mean + g_l_ci,
        color="#d95f02",
        alpha=0.2,
    )
    ax_mod.plot(
        cutoffs,
        g_s_mean,
        "o-",
        color="#1b9e77",
        linewidth=2,
        label="Entanglement Springs",
    )
    ax_mod.fill_between(
        cutoffs,
        g_s_mean - g_s_ci,
        g_s_mean + g_s_ci,
        color="#1b9e77",
        alpha=0.2,
    )
    ax_mod.axvline(
        b,
        color="#7570b3",
        linestyle="--",
        label=rf"Kuhn Length $b \approx {b:.2f}$",
    )
    ax_mod.axvline(
        default_cutoff,
        color="#e7298a",
        linestyle=":",
        linewidth=1.8,
        label=rf"Default $r_{{\mathrm{{cut}}}} \approx {default_cutoff:.2f}$",
    )
    ax_mod.set_xlabel(r"Sampling Cut-off $r_{\mathrm{cut}}$ [distance units]")
    ax_mod.set_ylabel(r"Shear Modulus $G$ [MPa]")
    ax_mod.set_title(
        r"(a) Modulus $G(r_{\mathrm{cut}})$: Links vs. Springs", fontweight="bold"
    )
    ax_mod.grid(True, linestyle="--", alpha=0.5)
    ax_mod.legend(fontsize=9, frameon=True)
    plt.tight_layout()
    plt.show()


def plot_cutoff_residual_vs_rcut(data: Any) -> None:
    cutoffs = data["cutoffs"]
    b = float(data.get("b", 0.969))
    default_cutoff = float(data.get("default_cutoff", 2.594))
    res_l_mean = data["res_links_mean"]
    res_l_ci = data["res_links_ci"]
    res_s_mean = data["res_springs_mean"]
    res_s_ci = data["res_springs_ci"]
    fig, ax_res = plt.subplots(figsize=(7, 4.5))
    ax_res.plot(
        cutoffs,
        res_l_mean,
        "s-",
        color="#d95f02",
        linewidth=2,
        label="Links Pre-Strain",
    )
    ax_res.fill_between(
        cutoffs,
        res_l_mean - res_l_ci,
        res_l_mean + res_l_ci,
        color="#d95f02",
        alpha=0.2,
    )
    ax_res.plot(
        cutoffs,
        res_s_mean,
        "o-",
        color="#1b9e77",
        linewidth=2,
        label="Springs Pre-Strain",
    )
    ax_res.fill_between(
        cutoffs,
        res_s_mean - res_s_ci,
        res_s_mean + res_s_ci,
        color="#1b9e77",
        alpha=0.2,
    )
    ax_res.axvline(b, color="#7570b3", linestyle="--")
    ax_res.axvline(default_cutoff, color="#e7298a", linestyle=":")
    ax_res.set_xlabel(r"Sampling Cut-off $r_{\mathrm{cut}}$ [distance units]")
    ax_res.set_ylabel(
        r"Initial Residual Norm $\|\mathbf{R}_{\mathrm{init}}\|$")
    ax_res.set_title(
        "(b) Topological Pre-Strain Distortion",
        fontweight="bold")
    ax_res.grid(True, linestyle="--", alpha=0.5)
    ax_res.legend(fontsize=9, frameon=True)
    plt.tight_layout()
    plt.show()


def plot_cutoff_strand_length_interaction(data: Any) -> None:
    cutoffs = data["cutoffs"]
    b = float(data.get("b", 0.969))
    g_l_mean = data["g_links_mean"]
    g_n30_mean = data["g_n30_mean"]
    g_n70_mean = data["g_n70_mean"]
    fig, ax_len = plt.subplots(figsize=(7, 4.5))
    # Normalize to value at cutoff = 1.2 (contact plateau)
    c12_idx = np.argmin(np.abs(cutoffs - 1.2))
    norm_n30 = g_n30_mean / g_n30_mean[c12_idx]
    norm_n70 = g_n70_mean / g_n70_mean[c12_idx]
    norm_n50 = g_l_mean / g_l_mean[c12_idx]

    ax_len.plot(
        cutoffs,
        norm_n30,
        "^-",
        color="#386cb0",
        linewidth=2,
        label=r"Short Strands ($N=30$, $G_c$-dominated)",
    )
    ax_len.plot(
        cutoffs,
        norm_n50,
        "s-",
        color="#d95f02",
        linewidth=2,
        label=r"Intermediate ($N=50$)",
    )
    ax_len.plot(
        cutoffs,
        norm_n70,
        "d-",
        color="#e41a1c",
        linewidth=2,
        label=r"Long Strands ($N=70$, $G_e$-dominated)",
    )
    ax_len.axhline(1.0, color="gray", linestyle="--", alpha=0.7)
    ax_len.axvline(b, color="#7570b3", linestyle="--")
    ax_len.set_xlabel(r"Sampling Cut-off $r_{\mathrm{cut}}$ [distance units]")
    ax_len.set_ylabel(r"Normalized Modulus $G / G(r_{\mathrm{cut}}=1.2)$")
    ax_len.set_title(
        r"(c) Relative Inflation vs. Strand Length $N$",
        fontweight="bold")
    ax_len.grid(True, linestyle="--", alpha=0.5)
    ax_len.legend(fontsize=9, frameon=True)
    plt.tight_layout()
    plt.show()


def plot_cutoff_network_parameter_interaction(data: Any) -> None:
    cutoffs = data["cutoffs"]
    g_l_mean = data["g_links_mean"]
    g_l_ci = data["g_links_ci"]
    g_p82_mean = data["g_p82_mean"]
    g_p82_ci = data["g_p82_ci"]
    g_r75_mean = data["g_r75_mean"]
    g_r75_ci = data["g_r75_ci"]
    fig, ax_net = plt.subplots(figsize=(7, 4.5))
    ax_net.plot(
        cutoffs,
        g_l_mean,
        "s-",
        color="#d95f02",
        linewidth=2,
        label=r"Ideal ($p=0.95, r=1.00$)",
    )
    ax_net.fill_between(
        cutoffs,
        g_l_mean - g_l_ci,
        g_l_mean + g_l_ci,
        color="#d95f02",
        alpha=0.15,
    )
    ax_net.plot(
        cutoffs,
        g_p82_mean,
        "v-",
        color="#984ea3",
        linewidth=2,
        label=r"Incomplete ($p=0.82, r=1.00$)",
    )
    ax_net.fill_between(
        cutoffs,
        g_p82_mean - g_p82_ci,
        g_p82_mean + g_p82_ci,
        color="#984ea3",
        alpha=0.15,
    )
    ax_net.plot(
        cutoffs,
        g_r75_mean,
        "x-",
        color="#ff7f00",
        linewidth=2,
        label=r"Imbalanced ($p=0.90, r=0.75$)",
    )
    ax_net.fill_between(
        cutoffs,
        g_r75_mean - g_r75_ci,
        g_r75_mean + g_r75_ci,
        color="#ff7f00",
        alpha=0.15,
    )
    ax_net.set_xlabel(r"Sampling Cut-off $r_{\mathrm{cut}}$ [distance units]")
    ax_net.set_ylabel(r"Shear Modulus $G$ [MPa]")
    ax_net.set_title(
        r"(d) Conversion $p$ & Stoichiometric Imbalance $r$", fontweight="bold"
    )
    ax_net.grid(True, linestyle="--", alpha=0.5)
    ax_net.legend(fontsize=9, frameon=True)
    plt.tight_layout()
    plt.show()


def plot_cutoff_pair_distance_distribution(data: Any) -> None:
    bin_centers = data["bin_centers"]
    hist_12 = data["hist_12"]
    hist_20 = data["hist_20"]
    hist_30 = data["hist_30"]
    fig, ax_dist = plt.subplots(figsize=(7, 4.5))
    ax_dist.plot(
        bin_centers,
        hist_12,
        "-",
        color="#1b9e77",
        linewidth=2,
        label=r"$r_{\mathrm{cut}} = 1.2$",
    )
    ax_dist.plot(
        bin_centers,
        hist_20,
        "-",
        color="#386cb0",
        linewidth=2,
        label=r"$r_{\mathrm{cut}} = 2.0$",
    )
    ax_dist.plot(
        bin_centers,
        hist_30,
        "-",
        color="#d95f02",
        linewidth=2,
        label=r"$r_{\mathrm{cut}} = 3.0$",
    )
    # Theoretical r^2 shell curve
    r_theory = np.linspace(0.01, 3.5, 100)
    ax_dist.plot(
        r_theory,
        0.45 * (r_theory**2) / (3.0**3 / 3.0),
        ":",
        color="black",
        label=r"Spherical Shell $\propto r^2$",
    )
    ax_dist.set_xlabel(r"Pair Separation Distance $r_{ij}$ [distance units]")
    ax_dist.set_ylabel(r"Probability Density $P(r)$")
    ax_dist.set_title(r"(e) Sampled Pair Separation $P(r)$", fontweight="bold")
    ax_dist.set_xlim(0, 3.5)
    ax_dist.grid(True, linestyle="--", alpha=0.5)
    ax_dist.legend(fontsize=9, frameon=True)
    plt.tight_layout()
    plt.show()


def plot_cutoff_solver_iterations(data: Any) -> None:
    cutoffs = data["cutoffs"]
    iter_l_mean = data["iter_links_mean"]
    iter_l_ci = data["iter_links_ci"]
    iter_s_mean = data["iter_springs_mean"]
    iter_s_ci = data["iter_springs_ci"]
    fig, ax_iter = plt.subplots(figsize=(7, 4.5))
    ax_iter.plot(
        cutoffs,
        iter_l_mean,
        "s-",
        color="#d95f02",
        linewidth=2,
        label="Links Iterations",
    )
    ax_iter.fill_between(
        cutoffs,
        iter_l_mean - iter_l_ci,
        iter_l_mean + iter_l_ci,
        color="#d95f02",
        alpha=0.15,
    )
    ax_iter.plot(
        cutoffs,
        iter_s_mean,
        "o-",
        color="#1b9e77",
        linewidth=2,
        label="Springs Iterations",
    )
    ax_iter.fill_between(
        cutoffs,
        iter_s_mean - iter_s_ci,
        iter_s_mean + iter_s_ci,
        color="#1b9e77",
        alpha=0.15,
    )
    ax_iter.set_xlabel(r"Sampling Cut-off $r_{\mathrm{cut}}$ [distance units]")
    ax_iter.set_ylabel("CG Iterations to Convergence")
    ax_iter.set_title("(f) Relaxation Solver Convergence", fontweight="bold")
    ax_iter.grid(True, linestyle="--", alpha=0.5)
    ax_iter.legend(fontsize=9, frameon=True)

    plt.tight_layout()
    plt.show()


def plot_inactive_modulus_convergence(data: Any) -> None:
    convs = data["conversions"]
    g_unf_nosimp_m = data["g_unf_nosimp_mean"]
    g_unf_nosimp_ci = data["g_unf_nosimp_ci"]
    g_unf_inact_m = data["g_unf_inact_mean"]
    g_unf_inact_ci = data["g_unf_inact_ci"]
    g_unf_x2f_m = data["g_unf_x2f_mean"]
    g_filt_nosimp_m = data["g_filt_nosimp_mean"]
    g_filt_nosimp_ci = data["g_filt_nosimp_ci"]
    g_filt_x2f_m = data["g_filt_x2f_mean"]
    g_filt_x2f_ci = data["g_filt_x2f_ci"]
    c_naive = "#d62728"
    c_filt_raw = "#1f77b4"
    c_inact = "#ff7f0e"
    c_truth = "#2ca02c"
    c_unf_x2f = "#9467bd"
    fig, ax1 = plt.subplots(figsize=(7, 5))
    ax1.plot(
        convs,
        g_unf_nosimp_m,
        "o--",
        color=c_naive,
        lw=2.2,
        ms=6,
        label="Unfiltered Sampling (NO_SIMPLIFICATION)",
    )
    ax1.fill_between(
        convs,
        g_unf_nosimp_m - g_unf_nosimp_ci,
        g_unf_nosimp_m + g_unf_nosimp_ci,
        color=c_naive,
        alpha=0.15,
    )

    ax1.plot(
        convs,
        g_filt_nosimp_m,
        "s:",
        color=c_filt_raw,
        lw=2.0,
        ms=6,
        label="Filtered Sampling (NO_SIMPLIFICATION)",
    )
    ax1.fill_between(
        convs,
        g_filt_nosimp_m - g_filt_nosimp_ci,
        g_filt_nosimp_m + g_filt_nosimp_ci,
        color=c_filt_raw,
        alpha=0.15,
    )

    ax1.plot(
        convs,
        g_unf_inact_m,
        "^--",
        color=c_inact,
        lw=1.8,
        ms=6,
        label="Unfiltered Sampling (INACTIVE_ONLY)",
    )
    ax1.fill_between(
        convs,
        g_unf_inact_m - g_unf_inact_ci,
        g_unf_inact_m + g_unf_inact_ci,
        color=c_inact,
        alpha=0.15,
    )

    ax1.plot(
        convs,
        g_filt_x2f_m,
        "d-",
        color=c_truth,
        lw=2.5,
        ms=7,
        label="Filtered Sampling (INACTIVE_THEN_X2F)",
    )
    ax1.fill_between(
        convs,
        g_filt_x2f_m - g_filt_x2f_ci,
        g_filt_x2f_m + g_filt_x2f_ci,
        color=c_truth,
        alpha=0.2,
    )

    ax1.plot(
        convs,
        g_unf_x2f_m,
        "x-.",
        color=c_unf_x2f,
        lw=1.8,
        ms=6,
        label="Unfiltered Sampling (INACTIVE_THEN_X2F)",
    )

    ax1.set_xlabel("Crosslinker Conversion $p$", fontsize=12)
    ax1.set_ylabel("Equilibrium Shear Modulus $G$ [kPa]", fontsize=12)
    ax1.set_title(
        "Shear Modulus across Network Conversion", fontsize=13, fontweight="bold"
    )
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper left", frameon=True, fontsize=9.5)
    plt.tight_layout()
    plt.show()


def plot_inactive_relative_error(data: Any) -> None:
    convs = data["conversions"]
    g_unf_nosimp_m = data["g_unf_nosimp_mean"]
    g_filt_nosimp_m = data["g_filt_nosimp_mean"]
    g_filt_x2f_m = data["g_filt_x2f_mean"]
    c_naive = "#d62728"
    c_filt_raw = "#1f77b4"
    fig, ax2 = plt.subplots(figsize=(7, 5))
    g_truth = g_filt_x2f_m
    # Avoid division by near-zero at lowest conversion if unpercolated
    inflation_unf = np.where(
        g_truth > 0.1, (g_unf_nosimp_m - g_truth) / g_truth * 100.0, 0.0
    )
    inflation_filt_raw = np.where(
        g_truth > 0.1, (g_filt_nosimp_m - g_truth) / g_truth * 100.0, 0.0
    )

    ax2.plot(
        convs,
        inflation_unf,
        "o--",
        color=c_naive,
        lw=2.2,
        ms=6,
        label="Unfiltered (NO_SIMPLIFICATION) Error",
    )
    ax2.plot(
        convs,
        inflation_filt_raw,
        "s:",
        color=c_filt_raw,
        lw=2.0,
        ms=6,
        label="Filtered (NO_SIMPLIFICATION) Discrepancy",
    )
    ax2.axhline(0, color="black", linestyle="-", lw=1.0, alpha=0.7)

    ax2.set_xlabel("Crosslinker Conversion $p$", fontsize=12)
    ax2.set_ylabel(
        "Artificial Modulus Inflation $\\Delta G / G_{\\text{true}}$ [%]",
        fontsize=12,
    )
    ax2.set_title(
        "Artificial Stiffening from Trapped Dangling Chains",
        fontsize=13,
        fontweight="bold",
    )
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", frameon=True, fontsize=10)

    plt.tight_layout()
    plt.show()


def plot_inactive_topological_fractions(data: Any) -> None:
    convs = data["conversions"]
    sp_unf_nosimp_m = data["sp_unf_nosimp_mean"]
    sp_unf_inact_m = data["sp_unf_inact_mean"]
    sp_unf_x2f_m = data["sp_unf_x2f_mean"]
    sp_filt_nosimp_m = data["sp_filt_nosimp_mean"]
    sp_filt_x2f_m = data["sp_filt_x2f_mean"]
    c_naive = "#d62728"
    c_filt_raw = "#1f77b4"
    c_inact = "#ff7f0e"
    c_truth = "#2ca02c"
    c_unf_x2f = "#9467bd"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel 1: Number of springs vs. conversion
    ax1.plot(
        convs,
        sp_unf_nosimp_m,
        "o--",
        color=c_naive,
        lw=2.0,
        ms=6,
        label="Initial Unfiltered Springs",
    )
    ax1.plot(
        convs,
        sp_filt_nosimp_m,
        "s:",
        color=c_filt_raw,
        lw=2.0,
        ms=6,
        label="Initial Filtered Springs",
    )
    ax1.plot(
        convs,
        sp_unf_inact_m,
        "^--",
        color=c_inact,
        lw=2.0,
        ms=6,
        label="Unfiltered + INACTIVE_ONLY",
    )
    ax1.plot(
        convs,
        sp_filt_x2f_m,
        "d-",
        color=c_truth,
        lw=2.5,
        ms=7,
        label="Filtered + INACTIVE_THEN_X2F",
    )
    ax1.plot(
        convs,
        sp_unf_x2f_m,
        "x-.",
        color=c_unf_x2f,
        lw=2.0,
        ms=6,
        label="Unfiltered + INACTIVE_THEN_X2F",
    )

    ax1.set_xlabel("Crosslinker Conversion $p$", fontsize=12)
    ax1.set_ylabel("Number of Springs $N_{\\text{springs}}$", fontsize=12)
    ax1.set_title(
        "Remaining Springs after Topological Pruning",
        fontsize=13,
        fontweight="bold",
    )
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="center left", frameon=True, fontsize=10)

    # Panel 2: Pruned percentage
    pruned_unf_x2f = (1.0 - sp_unf_x2f_m / sp_unf_nosimp_m) * 100.0
    pruned_filt_x2f = (1.0 - sp_filt_x2f_m / sp_filt_nosimp_m) * 100.0

    ax2.plot(
        convs,
        pruned_unf_x2f,
        "x-.",
        color=c_unf_x2f,
        lw=2.2,
        ms=7,
        label="Unfiltered Pruning (INACTIVE_THEN_X2F)",
    )
    ax2.plot(
        convs,
        pruned_filt_x2f,
        "d-",
        color=c_truth,
        lw=2.2,
        ms=7,
        label="Filtered Pruning (INACTIVE_THEN_X2F)",
    )

    ax2.set_xlabel("Crosslinker Conversion $p$", fontsize=12)
    ax2.set_ylabel("Springs Pruned / Collapsed [%]", fontsize=12)
    ax2.set_title(
        "Network Simplification Fraction vs. Conversion",
        fontsize=13,
        fontweight="bold",
    )
    ax2.yaxis.set_major_formatter(ticker.PercentFormatter())
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", frameon=True, fontsize=10.5)

    plt.tight_layout()
    plt.show()


def plot_inactive_solver_complexity(data: Any) -> None:
    convs = data["conversions"]
    links_unf_nosimp_m = data["links_unf_nosimp_mean"]
    links_unf_x2f_m = data["links_unf_x2f_mean"]
    links_filt_x2f_m = data["links_filt_x2f_mean"]
    iter_unf_nosimp_m = data["iter_unf_nosimp_mean"]
    iter_unf_inact_m = data["iter_unf_inact_mean"]
    iter_unf_x2f_m = data["iter_unf_x2f_mean"]
    iter_filt_nosimp_m = data["iter_filt_nosimp_mean"]
    iter_filt_x2f_m = data["iter_filt_x2f_mean"]
    c_naive = "#d62728"
    c_filt_raw = "#1f77b4"
    c_inact = "#ff7f0e"
    c_truth = "#2ca02c"
    c_unf_x2f = "#9467bd"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel 1: Degrees of freedom (3 * N_links)
    ax1.plot(
        convs,
        3 * links_unf_nosimp_m,
        "o--",
        color=c_naive,
        lw=2.0,
        ms=6,
        label="Unfiltered Unsimplified (3 $N_{\\text{links}}$)",
    )
    ax1.plot(
        convs,
        3 * links_filt_x2f_m,
        "d-",
        color=c_truth,
        lw=2.2,
        ms=7,
        label="Filtered + INACTIVE_THEN_X2F",
    )
    ax1.plot(
        convs,
        3 * links_unf_x2f_m,
        "x-.",
        color=c_unf_x2f,
        lw=2.0,
        ms=6,
        label="Unfiltered + INACTIVE_THEN_X2F",
    )

    ax1.set_xlabel("Crosslinker Conversion $p$", fontsize=12)
    ax1.set_ylabel(
        "Linear System Dimension ($3 N_{\\text{links}}$)",
        fontsize=12)
    ax1.set_title(
        "Degrees of Freedom in Sparse Force Balance",
        fontsize=13,
        fontweight="bold",
    )
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="upper left", frameon=True, fontsize=10)

    # Panel 2: Solver Iterations
    ax2.plot(
        convs,
        iter_unf_nosimp_m,
        "o--",
        color=c_naive,
        lw=2.0,
        ms=6,
        label="Unfiltered NO_SIMPLIFICATION",
    )
    ax2.plot(
        convs,
        iter_filt_nosimp_m,
        "s:",
        color=c_filt_raw,
        lw=2.0,
        ms=6,
        label="Filtered NO_SIMPLIFICATION",
    )
    ax2.plot(
        convs,
        iter_unf_inact_m,
        "^--",
        color=c_inact,
        lw=2.0,
        ms=6,
        label="Unfiltered INACTIVE_ONLY",
    )
    ax2.plot(
        convs,
        iter_filt_x2f_m,
        "d-",
        color=c_truth,
        lw=2.2,
        ms=7,
        label="Filtered INACTIVE_THEN_X2F",
    )
    ax2.plot(
        convs,
        iter_unf_x2f_m,
        "x-.",
        color=c_unf_x2f,
        lw=2.0,
        ms=6,
        label="Unfiltered INACTIVE_THEN_X2F",
    )

    ax2.set_xlabel("Crosslinker Conversion $p$", fontsize=12)
    ax2.set_ylabel("Relaxation Iterations", fontsize=12)
    ax2.set_title(
        "Solver Convergence & Pruning Iterations",
        fontsize=13,
        fontweight="bold",
    )
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right", frameon=True, fontsize=9.5)

    plt.tight_layout()
    plt.show()


def plot_sampling_memory_scaling(data: Any) -> None:
    """Plot peak RSS memory vs. network atom count for discrete and on-demand sampling."""
    scaling_atoms = data["scaling_atoms"]
    scaling_n_repeats = int(data.get("scaling_n_repeats", 1))
    scaling_rss_disc = data["scaling_rss_discrete"]
    scaling_rss_od = data["scaling_rss_on_demand"]
    disc_iqr = (
        data.get("scaling_rss_q25_discrete"),
        data.get("scaling_rss_q75_discrete"),
    )
    od_iqr = (
        data.get("scaling_rss_q25_on_demand"),
        data.get("scaling_rss_q75_on_demand"),
    )

    fig, ax_mem = plt.subplots(figsize=(7, 4.5))
    _plot_with_iqr(
        ax_mem,
        scaling_atoms,
        scaling_rss_disc,
        disc_iqr,
        "o-",
        "#d95f02",
        "Discrete (Full Universe)",
    )
    _plot_with_iqr(
        ax_mem,
        scaling_atoms,
        scaling_rss_od,
        od_iqr,
        "^-",
        "#7570b3",
        "On-Demand (Bounding Box)",
    )

    max_atoms = scaling_atoms[-1]
    pct_savings = (
        1.0 - float(scaling_rss_od[-1]) / float(scaling_rss_disc[-1])
    ) * 100.0
    ax_mem.annotate(
        f"{pct_savings:.0f}% Memory\nReduction",
        xy=(max_atoms, float(scaling_rss_od[-1])),
        xytext=(max_atoms * 0.55, 75),
        arrowprops=dict(arrowstyle="->", color="#7570b3", lw=1.5),
        fontsize=10,
        fontweight="bold",
        color="#7570b3",
        bbox=dict(boxstyle="round,pad=0.3", fc="#efedf5", ec="#7570b3", lw=1),
    )

    ax_mem.set_xscale("log")
    ax_mem.set_yscale("log")
    ax_mem.xaxis.set_major_formatter(_K_FORMATTER)
    ax_mem.set_xlabel(r"Total Atoms in Network ($N_{\mathrm{atoms}}$)")
    ax_mem.set_ylabel("Peak Resident Memory (RSS) [MB]")
    ax_mem.set_title(
        f"Peak Memory Consumption (median & IQR, N={scaling_n_repeats} seeds)"
    )
    ax_mem.grid(True, which="both", linestyle="--", alpha=0.5)
    ax_mem.legend(frameon=True)
    plt.tight_layout()
    plt.show()


def plot_sampling_runtime_scaling(data: Any) -> None:
    """Plot setup & sampling runtime vs. network atom count for discrete and on-demand sampling."""
    scaling_atoms = data["scaling_atoms"]
    scaling_n_repeats = int(data.get("scaling_n_repeats", 1))
    scaling_time_disc = data["scaling_time_discrete"]
    scaling_time_od = data["scaling_time_on_demand"]
    disc_iqr = (
        data.get("scaling_time_q25_discrete"),
        data.get("scaling_time_q75_discrete"),
    )
    od_iqr = (
        data.get("scaling_time_q25_on_demand"),
        data.get("scaling_time_q75_on_demand"),
    )

    fig, ax_time = plt.subplots(figsize=(7, 4.5))
    _plot_with_iqr(
        ax_time,
        scaling_atoms,
        scaling_time_disc,
        disc_iqr,
        "o-",
        "#d95f02",
        "Discrete (Full Universe)",
    )
    _plot_with_iqr(
        ax_time,
        scaling_atoms,
        scaling_time_od,
        od_iqr,
        "^-",
        "#7570b3",
        "On-Demand (Bounding Box)",
    )

    ax_time.set_xscale("log")
    ax_time.set_yscale("log")
    ax_time.xaxis.set_major_formatter(_K_FORMATTER)
    ax_time.set_xlabel(r"Total Atoms in Network ($N_{\mathrm{atoms}}$)")
    ax_time.set_ylabel("Sampling & Force-Balance Setup Time [s]")
    ax_time.set_title(
        f"Sampling & Force-Balance Setup Time (median & IQR, N={scaling_n_repeats} seeds)"
    )
    ax_time.grid(True, which="both", linestyle="--", alpha=0.5)
    ax_time.legend(frameon=True)
    plt.tight_layout()
    plt.show()


def plot_sampling_strand_length_scaling(data: Any) -> None:
    """Plot runtime and memory vs. strand length at fixed total atom count."""
    strlen_lens = data.get("strlen_strand_lens")
    if strlen_lens is None:
        return
    strlen_atoms = data.get("strlen_atoms")
    strlen_n_repeats = int(data.get("strlen_n_repeats", 1))

    fig, (ax_sl_t, ax_sl_m) = plt.subplots(1, 2, figsize=(11.5, 4.5))
    for mode, fmt, col, lbl in (
        ("discrete", "o-", "#d95f02", "Discrete (Full Universe)"),
        ("on_demand", "^-", "#7570b3", "On-Demand (Bounding Box)"),
    ):
        t_vals = data.get(f"strlen_time_{mode}")
        r_vals = data.get(f"strlen_rss_{mode}")
        if t_vals is None or r_vals is None:
            continue
        t_iqr = (
            data.get(f"strlen_time_q25_{mode}"),
            data.get(f"strlen_time_q75_{mode}"),
        )
        r_iqr = (
            data.get(
                f"strlen_rss_q25_{mode}"), data.get(
                f"strlen_rss_q75_{mode}"))
        _plot_with_iqr(ax_sl_t, strlen_lens, t_vals, t_iqr, fmt, col, lbl)
        _plot_with_iqr(ax_sl_m, strlen_lens, r_vals, r_iqr, fmt, col, lbl)

    n_atoms_label = (
        f"{int(np.median(strlen_atoms)):,}" if strlen_atoms is not None else "300,000"
    )
    for ax in (ax_sl_t, ax_sl_m):
        ax.set_xscale("log")
        ax.set_xticks(strlen_lens)
        ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        ax.set_xlabel(r"Strand Length ($N_{\mathrm{beads}}$)")
        ax.grid(True, which="both", linestyle="--", alpha=0.5)
        ax.legend(frameon=True)

    ax_sl_t.set_yscale("log")
    ax_sl_t.set_ylabel("Sampling & Force-Balance Setup Time [s]")
    ax_sl_t.set_title(
        f"Runtime vs. Strand Length ($N_{{\\mathrm{{atoms}}}} \\approx$ {n_atoms_label}, N={strlen_n_repeats})"
    )
    ax_sl_m.set_ylabel("Peak Resident Memory (RSS) [MB]")
    ax_sl_m.set_title(
        f"Memory vs. Strand Length ($N_{{\\mathrm{{atoms}}}} \\approx$ {n_atoms_label}, N={strlen_n_repeats})"
    )
    plt.tight_layout()
    plt.show()


def plot_sampling_phase_breakdown(data: Any) -> None:
    """Plot stacked per-phase runtime breakdown at the largest benchmarked network."""
    scaling_atoms = data["scaling_atoms"]
    t_samp_d = data.get("scaling_time_sample_discrete")
    if t_samp_d is None:
        return

    methods = ("discrete", "on_demand")
    labels = ["Discrete", "On-Demand"]
    t_univ = np.array(
        [float(data[f"scaling_time_universe_{m}"][-1]) for m in methods])
    t_samp = np.array(
        [float(data[f"scaling_time_sample_{m}"][-1]) for m in methods])
    t_build = np.array(
        [float(data[f"scaling_time_build_{m}"][-1]) for m in methods])

    fig, ax_ph = plt.subplots(figsize=(6.5, 4.5))
    x = np.arange(len(methods))
    w = 0.45
    ax_ph.bar(
        x, t_univ, w, label="Universe Realization (get_universe)", color="#e78ac3"
    )
    ax_ph.bar(
        x, t_samp, w, bottom=t_univ, label="Entanglement Sampling", color="#66c2a5"
    )
    ax_ph.bar(
        x,
        t_build,
        w,
        bottom=t_univ + t_samp,
        label="Force-Balance Build & Validation",
        color="#8da0cb",
    )
    for xi, u_v, s_v, b_v in zip(x, t_univ, t_samp, t_build):
        tot = u_v + s_v + b_v
        ax_ph.text(xi,
                   tot * 1.02,
                   f"{tot:.2f} s",
                   ha="center",
                   va="bottom",
                   fontsize=9)

    ax_ph.set_xticks(x)
    ax_ph.set_xticklabels(labels)
    ax_ph.set_ylabel("Median Elapsed Time [s]")
    ax_ph.set_title(
        f"Phase Breakdown at $N_{{\\mathrm{{atoms}}}} = {int(scaling_atoms[-1]):,}$"
    )
    ax_ph.set_ylim(0, float(np.max(t_univ + t_samp + t_build)) * 1.16)
    ax_ph.grid(True, axis="y", linestyle="--", alpha=0.5)
    ax_ph.legend(frameon=True)
    plt.tight_layout()
    plt.show()


def plot_sampling_distributions(data: Any) -> None:
    """Plot contour fraction alpha and pair distance r distributions for discrete vs. on-demand."""
    alpha_od = data["alpha_on_demand"]
    dist_od = data["dist_on_demand"]
    alpha_disc = data["alpha_discrete"]
    dist_disc = data["dist_discrete"]

    fig, (ax_alpha, ax_dist) = plt.subplots(1, 2, figsize=(11, 4.5))
    bins_a = np.linspace(0, 1, 26)
    ax_alpha.hist(
        alpha_disc,
        bins=bins_a,
        density=True,
        alpha=0.45,
        color="#d95f02",
        edgecolor="#a63603",
        label="Discrete Sampling",
    )
    ax_alpha.hist(
        alpha_od,
        bins=bins_a,
        density=True,
        histtype="step",
        linewidth=2.2,
        color="#7570b3",
        label="On-Demand Sampling",
    )
    ax_alpha.axhline(
        1.0, color="k", linestyle=":", linewidth=1.2, label="Uniform Theory (1.0)"
    )
    ax_alpha.set_xlabel(r"Strand Contour Position ($\alpha$)")
    ax_alpha.set_ylabel("Probability Density")
    ax_alpha.set_title(r"Distribution of Sampled Contour Positions $\alpha$")
    ax_alpha.set_xlim(0, 1)
    ax_alpha.set_ylim(0, 1.4)
    ax_alpha.grid(True, linestyle="--", alpha=0.4)
    ax_alpha.legend(frameon=True)

    r_max = float(max(np.max(dist_disc), np.max(dist_od)))
    bins_r = np.linspace(0, r_max * 1.02, 26)
    ax_dist.hist(
        dist_disc,
        bins=bins_r,
        density=True,
        alpha=0.45,
        color="#d95f02",
        edgecolor="#a63603",
        label="Discrete Sampling",
    )
    ax_dist.hist(
        dist_od,
        bins=bins_r,
        density=True,
        histtype="step",
        linewidth=2.2,
        color="#7570b3",
        label="On-Demand Sampling",
    )
    r_grid = np.linspace(0, r_max, 200)
    pdf_theory = 3.0 * (r_grid**2) / (r_max**3)
    ax_dist.plot(
        r_grid,
        pdf_theory,
        "k:",
        linewidth=1.5,
        label=r"3D Shell Theory ($3r^2 / r_{\mathrm{cut}}^3$)",
    )
    ax_dist.set_xlabel(r"Spatial Separation Distance ($r$)")
    ax_dist.set_ylabel("Probability Density")
    ax_dist.set_title(r"Distribution of Sampled Pair Distances $r$")
    ax_dist.grid(True, linestyle="--", alpha=0.4)
    ax_dist.legend(frameon=True)

    plt.tight_layout()
    plt.show()


def plot_sampling_conversion_sweep(data: Any) -> None:
    """Plot shear modulus, active entanglements, and active fractions across conversion."""
    conversions = data["conversions"]
    target_n_ent = int(data["target_n_ent"])
    n_repeats = int(data.get("n_repeats", 1))

    g_disc = data["g_discrete"]
    w_act_disc = data["w_active_discrete"]
    n_ent_disc = data.get("n_ent_discrete")
    g_od_filt = data["g_od_filtered"]
    w_act_od_filt = data["w_active_od_filtered"]
    n_ent_od_filt = data["n_ent_od_filtered"]
    g_od_unfilt = data["g_od_unfiltered"]
    n_ent_od_unfilt = data["n_ent_od_unfiltered"]

    g_disc_ci = data.get("g_discrete_ci")
    w_act_disc_ci = data.get("w_active_discrete_ci")
    n_ent_disc_ci = data.get("n_ent_discrete_ci")
    g_od_filt_ci = data.get("g_od_filtered_ci")
    w_act_od_filt_ci = data.get("w_active_od_filtered_ci")
    n_ent_od_filt_ci = data.get("n_ent_od_filtered_ci")
    g_od_unfilt_ci = data.get("g_od_unfiltered_ci")
    n_ent_od_unfilt_ci = data.get("n_ent_od_unfiltered_ci")

    fig, (ax_mod, ax_ent, ax_frac) = plt.subplots(1, 3, figsize=(15, 4.6))

    # Panel 1: Shear Modulus G
    for y, ci, fmt, col, lbl in (
        (
            g_disc * 1e3,
            g_disc_ci * 1e3 if g_disc_ci is not None else None,
            "o-",
            "#d95f02",
            f"Discrete (Links, N={n_repeats})",
        ),
        (
            g_od_filt * 1e3,
            g_od_filt_ci * 1e3 if g_od_filt_ci is not None else None,
            "^-",
            "#7570b3",
            "On-Demand Filtered (Links)",
        ),
        (
            g_od_unfilt * 1e3,
            g_od_unfilt_ci * 1e3 if g_od_unfilt_ci is not None else None,
            "x--",
            "#984ea3",
            "On-Demand Unfiltered (Links)",
        ),
    ):
        ax_mod.plot(
            conversions, y, fmt, color=col, linewidth=2, markersize=6, label=lbl
        )
        if ci is not None:
            ax_mod.fill_between(
                conversions,
                y - ci,
                y + ci,
                color=col,
                alpha=0.18)

    ax_mod.set_xlabel(r"Crosslinker Conversion $p$")
    ax_mod.set_ylabel("Equilibrium Shear Modulus $G$ [kPa]")
    ax_mod.set_title("Shear Modulus vs. Conversion")
    ax_mod.grid(True, linestyle="--", alpha=0.5)
    ax_mod.legend(frameon=True, fontsize=8)

    # Panel 2: Effective Active Entanglement Count
    for y, ci, fmt, col, lbl in (
        (n_ent_disc, n_ent_disc_ci, "o-", "#d95f02", "Discrete (Filtered)"),
        (n_ent_od_filt, n_ent_od_filt_ci, "^-", "#7570b3", "On-Demand (Filtered)"),
        (
            n_ent_od_unfilt,
            n_ent_od_unfilt_ci,
            "x--",
            "#984ea3",
            "On-Demand (Unfiltered + Simplified)",
        ),
    ):
        if y is None:
            continue
        ax_ent.plot(
            conversions, y, fmt, color=col, linewidth=2, markersize=6, label=lbl
        )
        if ci is not None:
            ax_ent.fill_between(
                conversions,
                y - ci,
                y + ci,
                color=col,
                alpha=0.18)

    ax_ent.axhline(
        target_n_ent,
        color="k",
        linestyle=":",
        linewidth=1.2,
        label=f"Target ({target_n_ent})",
    )
    ax_ent.set_xlabel(r"Crosslinker Conversion $p$")
    ax_ent.set_ylabel("Active Entanglement Links")
    ax_ent.set_title("Active Entanglements After Simplification")
    ax_ent.grid(True, linestyle="--", alpha=0.5)
    ax_ent.legend(frameon=True, fontsize=8)

    # Panel 3: Active Weight Fraction & Retention
    ax_frac.plot(
        conversions,
        w_act_disc,
        "o-",
        color="#d95f02",
        linewidth=2,
        markersize=6,
        label=r"Discrete $w_{\mathrm{active}}$",
    )
    if w_act_disc_ci is not None:
        ax_frac.fill_between(
            conversions,
            w_act_disc - w_act_disc_ci,
            w_act_disc + w_act_disc_ci,
            color="#d95f02",
            alpha=0.18,
        )
    ax_frac.plot(
        conversions,
        w_act_od_filt,
        "^-",
        color="#7570b3",
        linewidth=2,
        markersize=6,
        label=r"On-Demand $w_{\mathrm{active}}$",
    )
    if w_act_od_filt_ci is not None:
        ax_frac.fill_between(
            conversions,
            w_act_od_filt - w_act_od_filt_ci,
            w_act_od_filt + w_act_od_filt_ci,
            color="#7570b3",
            alpha=0.18,
        )
    ax_frac.plot(
        conversions,
        n_ent_od_filt / float(target_n_ent),
        "v--",
        color="#984ea3",
        linewidth=1.5,
        markersize=5,
        label=r"On-Demand $N_{\mathrm{active}} / N_{\mathrm{target}}$",
    )

    ax_frac.set_xlabel(r"Crosslinker Conversion $p$")
    ax_frac.set_ylabel("Fraction")
    ax_frac.set_title(r"Network Activity & Entanglement Retention")
    ax_frac.set_ylim(-0.05, 1.05)
    ax_frac.grid(True, linestyle="--", alpha=0.5)
    ax_frac.legend(frameon=True, fontsize=8)

    plt.tight_layout()
    plt.show()


_RVE_COLORS = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#a6761d",
]
_RVE_MARKERS = ["o", "s", "^", "D", "v", "P", "*"]


def _style_rve_box_axis(ax: plt.Axes, box_lengths: np.ndarray) -> None:
    """Apply standard log-scale and tick formatting for cubic box side length L."""
    ax.set_xscale("log")
    ax.set_xticks(box_lengths)
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.set_xlabel(r"Cubic Box Side Length $L$ [nm]")


def _style_rve_cv_axis(
    ax: plt.Axes,
    title: str,
    xlabel: str = r"Cubic Box Side Length $L$ [nm]",
    box_lengths: Optional[np.ndarray] = None,
    legend_fontsize: float = 8.5,
) -> None:
    """Apply reference threshold lines and log-log styling for CV scatter plots."""
    ax.axhline(5.0, color="gray", linestyle="--", linewidth=1.1, label="5% CV")
    ax.axhline(2.0, color="black", linestyle=":", linewidth=1.1, label="2% CV")
    if box_lengths is not None:
        _style_rve_box_axis(ax, box_lengths)
    else:
        ax.set_xscale("log")
        ax.set_xlabel(xlabel)
    ax.set_yscale("log")
    ax.set_ylabel(r"Sample-to-Sample CV ($\sigma / \langle G \rangle$) [%]")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle="--", alpha=0.45)
    ax.legend(frameon=True, fontsize=legend_fontsize)


def _plot_rve_normalized_sweep(
    box_lengths: np.ndarray,
    g_mean: np.ndarray,
    g_std: np.ndarray,
    g_ci: np.ndarray,
    conv_labels: list[str],
    cv_labels: list[str],
    conv_title: str,
    cv_title: str,
    ylim: Tuple[float, float],
    conv_legend_fontsize: float = 8.5,
) -> None:
    """Plot a two-panel normalized modulus convergence and CV scatter comparison."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    for i, (conv_lbl, cv_lbl) in enumerate(zip(conv_labels, cv_labels)):
        c = _RVE_COLORS[i % len(_RVE_COLORS)]
        m = _RVE_MARKERS[i % len(_RVE_MARKERS)]
        g_ref = g_mean[i, -1]
        ax1.errorbar(
            box_lengths,
            g_mean[i] / g_ref,
            yerr=g_ci[i] / g_ref,
            fmt=f"{m}-",
            color=c,
            linewidth=2,
            markersize=6,
            capsize=3.5,
            elinewidth=1.2,
            label=conv_lbl,
        )

        cv_pct = np.where(
            g_mean[i] > 1e-9,
            100.0 * g_std[i] / g_mean[i],
            np.nan,
        )
        ax2.plot(
            box_lengths,
            cv_pct,
            f"{m}-",
            color=c,
            linewidth=2,
            markersize=6,
            label=cv_lbl,
        )

    ax1.axhline(1.0, color="black", linestyle=":", linewidth=1.1)
    ax1.axhspan(0.98, 1.02, color="gray", alpha=0.12, label=r"$\pm 2\%$ band")
    _style_rve_box_axis(ax1, box_lengths)
    ax1.set_ylabel(r"Normalized Modulus $G(L) / G(L_{\max})$")
    ax1.set_title(conv_title)
    ax1.set_ylim(*ylim)
    ax1.grid(True, which="both", linestyle="--", alpha=0.45)
    ax1.legend(frameon=True, fontsize=conv_legend_fontsize)

    _style_rve_cv_axis(ax2, title=cv_title, box_lengths=box_lengths)

    plt.tight_layout()
    plt.show()


def plot_rve_chain_length_study(data: Any) -> None:
    """Plot shear modulus convergence and CV vs. box size across chain lengths N."""
    box_lengths = data["box_lengths"]
    n_values = data["n_values"]
    n_r0_nm = data["n_r0_nm"]
    n_repeats = int(data.get("n_repeats", 12))
    g_mean = data["n_g_total_mean"]
    g_std = data["n_g_total_std"]
    g_ci = data["n_g_total_ci"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    for i, (n_len, r0) in enumerate(zip(n_values, n_r0_nm)):
        c = _RVE_COLORS[i % len(_RVE_COLORS)]
        m = _RVE_MARKERS[i % len(_RVE_MARKERS)]
        ax1.errorbar(
            box_lengths,
            g_mean[i],
            yerr=g_ci[i],
            fmt=f"{m}-",
            color=c,
            linewidth=2,
            markersize=6,
            capsize=3.5,
            elinewidth=1.2,
            label=rf"$N = {n_len}$ ($R_0 \approx {r0:.1f}\,\mathrm{{nm}}$)",
        )
        ax1.axhline(
            g_mean[i, -1],
            color=c,
            linestyle=":",
            linewidth=1.0,
            alpha=0.65,
        )

        cv_pct = np.where(
            g_mean[i] > 1e-9,
            100.0 * g_std[i] / g_mean[i],
            np.nan,
        )
        ax2.plot(
            box_lengths / r0,
            cv_pct,
            f"{m}-",
            color=c,
            linewidth=2,
            markersize=6,
            label=rf"$N = {n_len}$",
        )

    _style_rve_box_axis(ax1, box_lengths)
    ax1.set_ylabel(r"Shear Modulus $G$ [MPa]")
    ax1.set_title(
        f"Modulus Convergence vs. Box Size ($p=0.95, r=1.0, n={n_repeats}$)"
    )
    ax1.grid(True, which="both", linestyle="--", alpha=0.45)
    ax1.legend(frameon=True, fontsize=8.5)

    _style_rve_cv_axis(
        ax2,
        title=r"Relative Scatter vs. Normalized Box Size $L / R_0$",
        xlabel=r"Normalized Box Size $L / R_0$",
    )

    plt.tight_layout()
    plt.show()


def plot_rve_conversion_study(data: Any) -> None:
    """Plot normalized modulus convergence and CV vs. box size across conversions p."""
    box_lengths = data["box_lengths"]
    p_values = data["p_values"]
    n_repeats = int(data.get("n_repeats", 12))
    g_mean = data["p_g_total_mean"]
    g_std = data["p_g_total_std"]
    g_ci = data["p_g_total_ci"]

    conv_labels = [
        rf"$p = {p_val:.3g}$ ($G_\infty = {g_mean[i, -1]:.3f}\,\mathrm{{MPa}}$)"
        for i, p_val in enumerate(p_values)
    ]
    cv_labels = [rf"$p = {p_val:.3g}$" for p_val in p_values]
    _plot_rve_normalized_sweep(
        box_lengths=box_lengths,
        g_mean=g_mean,
        g_std=g_std,
        g_ci=g_ci,
        conv_labels=conv_labels,
        cv_labels=cv_labels,
        conv_title=f"Modulus Convergence vs. Conversion ($N=100, r=1.0, n={n_repeats}$)",
        cv_title="Relative Scatter Across Conversions",
        ylim=(0.7, 1.45),
    )


def plot_rve_stoichiometry_study(data: Any) -> None:
    """Plot normalized modulus convergence and CV vs. box size across stoichiometry r."""
    box_lengths = data["box_lengths"]
    r_values = data["r_values"]
    r_conversions = data["r_conversions"]
    n_repeats = int(data.get("n_repeats", 12))
    g_mean = data["r_g_total_mean"]
    g_std = data["r_g_total_std"]
    g_ci = data["r_g_total_ci"]

    conv_labels = [
        rf"$r = {r_val:.2f}$ ($p={p_val:.2f}$, $G_\infty = {g_mean[i, -1]:.3f}\,\mathrm{{MPa}}$)"
        for i, (r_val, p_val) in enumerate(zip(r_values, r_conversions))
    ]
    cv_labels = [rf"$r = {r_val:.2f}$" for r_val in r_values]
    _plot_rve_normalized_sweep(
        box_lengths=box_lengths,
        g_mean=g_mean,
        g_std=g_std,
        g_ci=g_ci,
        conv_labels=conv_labels,
        cv_labels=cv_labels,
        conv_title=f"Modulus Convergence vs. Stoichiometry ($N=100, n={n_repeats}$)",
        cv_title="Relative Scatter Across Stoichiometric Ratios",
        ylim=(0.55, 1.55),
        conv_legend_fontsize=8.0,
    )


def plot_rve_sample_count_study(data: Any) -> None:
    """Plot 95% CI precision vs. sample count and finite-size bias vs. box size."""
    from scipy import stats

    box_lengths = data["box_lengths"]
    g_mean = data["s_g_total_mean"]
    g_std = data["s_g_total_std"]
    g_ci = data["s_g_total_ci"]
    n_samples_max = int(data.get("n_samples_study", 24))

    sample_counts = np.array(
        [2, 3, 4, 6, 8, 12, 16, 20, n_samples_max], dtype=int
    )
    sample_counts = sample_counts[sample_counts <= n_samples_max]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    for idx, box_l in enumerate(box_lengths):
        c = _RVE_COLORS[idx % len(_RVE_COLORS)]
        m = _RVE_MARKERS[idx % len(_RVE_MARKERS)]
        rel_ci_pct = [
            100.0
            * float(stats.t.ppf(0.975, df=n_s - 1))
            * g_std[idx]
            / (np.sqrt(n_s) * g_mean[idx])
            for n_s in sample_counts
        ]
        ax1.plot(
            sample_counts,
            rel_ci_pct,
            f"{m}-",
            color=c,
            linewidth=1.8,
            markersize=5,
            label=rf"$L = {box_l:.0f}\,\mathrm{{nm}}$",
        )

    ax1.axhline(5.0, color="gray", linestyle="--", linewidth=1.1, label="5% CI")
    ax1.axhline(2.0, color="black", linestyle=":", linewidth=1.1, label="2% CI")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"Number of Independent Samples $n_{\mathrm{samples}}$")
    ax1.set_ylabel(r"Relative 95% CI Half-Width [%]")
    ax1.set_title(
        r"Statistical Precision vs. Number of Samples ($N=100, p=0.925$)"
    )
    ax1.grid(True, which="both", linestyle="--", alpha=0.45)
    ax1.legend(frameon=True, fontsize=8, ncol=2)

    g_ref = g_mean[-1]
    rel_bias_pct = 100.0 * np.abs(g_mean - g_ref) / g_ref
    rel_ci_max_pct = 100.0 * g_ci / g_ref
    ax2.plot(
        box_lengths,
        rel_bias_pct,
        "o-",
        color="#d95f02",
        linewidth=2,
        markersize=6,
        label=r"Finite-size bias $|\langle G(L)\rangle - G(L_{\max})| / G(L_{\max})$",
    )
    ax2.plot(
        box_lengths,
        rel_ci_max_pct,
        "s--",
        color="#1b9e77",
        linewidth=2,
        markersize=6,
        label=rf"95% CI half-width ($n={n_samples_max}$)",
    )
    ax2.axhline(
        2.0, color="black", linestyle=":", linewidth=1.1, label="2% threshold"
    )
    _style_rve_box_axis(ax2, box_lengths)
    ax2.set_yscale("log")
    ax2.set_ylabel("Relative Error [%]")
    ax2.set_title("Systematic Finite-Size Bias vs. Statistical Uncertainty")
    ax2.grid(True, which="both", linestyle="--", alpha=0.45)
    ax2.legend(frameon=True, fontsize=8.5)

    plt.tight_layout()
    plt.show()
