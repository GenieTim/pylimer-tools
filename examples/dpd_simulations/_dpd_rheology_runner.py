"""
Simulation runner for multi-scale DPD slip-spring rheology vs. NMA and MEHP.

Provides parallelized, cache-friendly study functions comparing:
1. Equilibrium Green-Kubo shear stress autocorrelation G_GK(t) via multi-tau Correlator
2. Uniaxial step/ramp box deformation stress relaxation G_deform(t)
3. Static-graph NormalModeAnalyzer (NMA) and MEHPSimulator equilibrium moduli

across uncrosslinked Rouse melts, entangled slip-spring melts (with dynamic MC
shifting/relocation vs. frozen slip-springs), and crosslinked f=4 networks (without
slip-springs, with trapped slip-springs, and with crosslink relocation enabled).
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from enum import Enum
import os
from pathlib import Path
import sys
import tempfile
import time
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Ensure examples/ is importable for _example_cache
_EXAMPLES_DIR = str(Path(__file__).resolve().parents[1])
if _EXAMPLES_DIR not in sys.path:
    sys.path.insert(0, _EXAMPLES_DIR)

from _example_cache import (  # noqa: E402
    compute_ensemble_statistics,
    configure_process_threading,
    get_example_cache_path,
    load_or_compute_cache,
    resolve_max_workers,
    worker_process_initializer,
)
import pylimer_tools_cpp as pt  # noqa: E402


class TopologyKind(Enum):
    """Polymer system topology kind."""

    MELT = "melt"
    NETWORK = "network"


class SystemCondition(Enum):
    """Simulation condition identifiers for the DPD rheology study."""

    MELT_NO_SS = "melt_no_ss"
    MELT_SS = "melt_ss"
    MELT_SS_FROZEN = "melt_ss_frozen"
    NET_NO_SS = "net_no_ss"
    NET_SS_LOCKED = "net_ss_locked"
    NET_SS_RELOC = "net_ss_reloc"


@dataclass(frozen=True)
class _ConditionSpec:
    """Configuration for a single rheology simulation condition."""

    condition: SystemCondition
    topology: TopologyKind
    use_slipsprings: bool
    with_mc_production: bool
    allow_reloc_in_network: bool
    run_green_kubo: bool
    run_nma_mehp: bool


_CONDITION_SPECS: Tuple[_ConditionSpec, ...] = (
    _ConditionSpec(
        condition=SystemCondition.MELT_NO_SS,
        topology=TopologyKind.MELT,
        use_slipsprings=False,
        with_mc_production=False,
        allow_reloc_in_network=False,
        run_green_kubo=True,
        run_nma_mehp=True,
    ),
    _ConditionSpec(
        condition=SystemCondition.MELT_SS,
        topology=TopologyKind.MELT,
        use_slipsprings=True,
        with_mc_production=True,
        allow_reloc_in_network=False,
        run_green_kubo=True,
        run_nma_mehp=True,
    ),
    _ConditionSpec(
        condition=SystemCondition.MELT_SS_FROZEN,
        topology=TopologyKind.MELT,
        use_slipsprings=True,
        with_mc_production=False,
        allow_reloc_in_network=False,
        run_green_kubo=False,
        run_nma_mehp=False,
    ),
    _ConditionSpec(
        condition=SystemCondition.NET_NO_SS,
        topology=TopologyKind.NETWORK,
        use_slipsprings=False,
        with_mc_production=False,
        allow_reloc_in_network=False,
        run_green_kubo=True,
        run_nma_mehp=True,
    ),
    _ConditionSpec(
        condition=SystemCondition.NET_SS_LOCKED,
        topology=TopologyKind.NETWORK,
        use_slipsprings=True,
        with_mc_production=True,
        allow_reloc_in_network=False,
        run_green_kubo=True,
        run_nma_mehp=True,
    ),
    _ConditionSpec(
        condition=SystemCondition.NET_SS_RELOC,
        topology=TopologyKind.NETWORK,
        use_slipsprings=True,
        with_mc_production=True,
        allow_reloc_in_network=True,
        run_green_kubo=False,
        run_nma_mehp=True,
    ),
)

# Physical and simulation defaults (matching standard Langeloth/Schneider DPD parameters)
_RHO_BEAD: float = 3.0
_DPD_A: float = 25.0
_DPD_SIGMA: float = 3.0
_SPRING_K: float = 2.0
_KBT: float = 1.0
_DT: float = 0.05
_SS_FRACTION: float = 0.20

_N_CHAINS_MELT: int = 45
_N_MELT: int = 45
_N_XLINKS: int = 28
_N_STRAND_NET: int = 35
_CONVERSION_P: float = 0.96

_N_STEPS_EQUIL: int = 3000
_N_STEPS_GK: int = 14000
_N_STEPS_RAMP: int = 200
_N_STEPS_HOLD: int = 8000
_GK_OUTPUT_EVERY: int = 1
_DEF_OUTPUT_EVERY: int = 20
_LAMBDA_STRETCH: float = 1.4
_NMA_TAU0: float = 1.4


def _build_melt_universe(
    seed: int,
    rho: float = _RHO_BEAD,
    k_spring: float = _SPRING_K,
    n_chains: int = _N_CHAINS_MELT,
    chain_length: int = _N_MELT,
) -> Tuple[pt.Universe, float]:
    """Build an uncrosslinked Gaussian polymer melt at DPD bead density."""
    b2_ideal = 3.0 * _KBT / k_spring
    b_ideal = float(np.sqrt(b2_ideal))
    n_beads = n_chains * chain_length
    vol = n_beads / rho
    box_l = vol ** (1.0 / 3.0)

    gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(int(seed))
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(b_ideal)
    gen.set_mean_squared_bead_distance(b2_ideal)
    gen.add_solvent_chains(n_chains, chain_length, solvent_atom_type=1)
    return gen.get_universe(), vol


def _build_network_universe(
    seed: int,
    rho: float = _RHO_BEAD,
    k_spring: float = _SPRING_K,
    n_xlinks: int = _N_XLINKS,
    strand_length: int = _N_STRAND_NET,
    conversion: float = _CONVERSION_P,
) -> Tuple[pt.Universe, float]:
    """Build a stoichiometric f=4 crosslinked polymer network at DPD bead density."""
    b2_ideal = 3.0 * _KBT / k_spring
    b_ideal = float(np.sqrt(b2_ideal))
    n_strands = 2 * n_xlinks
    n_beads = n_xlinks + n_strands * strand_length
    vol = n_beads / rho
    box_l = vol ** (1.0 / 3.0)

    gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(int(seed))
    gen.config_nr_of_mc_steps(0)
    gen.config_primary_loop_probability(0.0)
    gen.set_bead_distance(b_ideal)
    gen.set_mean_squared_bead_distance(b2_ideal)
    gen.add_crosslinkers(n_xlinks, 4, crosslinker_type=2)
    gen.add_strands(n_strands, [strand_length] * n_strands, strand_atom_type=1)
    gen.link_strands_to_conversion(conversion)
    return gen.get_universe(), vol


def _create_and_equilibrate_dpd(
    universe: pt.Universe,
    seed: int,
    spec: _ConditionSpec,
) -> Tuple[pt.DPDSimulator, float]:
    """Instantiate DPDSimulator, add slip-springs if enabled, and equilibrate."""
    dpd = pt.DPDSimulator(
        universe,
        crosslinker_type=2,
        slipspring_bond_type=9,
        seed=f"dpd_{seed}_{spec.condition.value}",
    )
    dpd.assume_box_large_enough()
    dpd.config_a(_DPD_A)
    dpd.config_sigma(_DPD_SIGMA)
    dpd.config_spring_constant(_SPRING_K)
    dpd.config_allow_relocation_in_network(spec.allow_reloc_in_network)
    if spec.use_slipsprings:
        n_ss = int(round(_SS_FRACTION * universe.get_nr_of_atoms()))
        dpd.create_slip_springs(n_ss, 9)
    dpd.config_num_steps_dpd(50)
    dpd.config_num_steps_mc(50)
    # Always enable MC during initial slip-spring equilibration so slip-springs relax
    dpd.run_simulation(n_steps=_N_STEPS_EQUIL, dt=_DT, with_MC=spec.use_slipsprings)

    n_backbone_bonds = universe.get_nr_of_bonds()
    bond_lengths = dpd.get_bond_lengths()[:n_backbone_bonds]
    b2_measured = float(np.mean(bond_lengths**2))
    return dpd, b2_measured


def _compute_green_kubo_modulus(
    dpd: pt.DPDSimulator,
    vol: float,
    with_mc: bool,
    norm_correlator: bool,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """Run NVT production trajectory and compute Green-Kubo G(t) via Correlator.

    Uses the 5-component traceless symmetric stress tensor estimator (3 shear
    components plus 3 normal stress differences weighted by 1/6, divided by 5).
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = os.path.join(tmpdir, "gk_stress.txt")
        cfg = pt.OutputConfiguration()
        cfg.int_values = [pt.ComputedIntValues.STEP]
        cfg.double_values = [
            pt.ComputedDoubleValues.TEMPERATURE,
            pt.ComputedDoubleValues.PRESSURE,
            pt.ComputedDoubleValues.STRESS_XY,
            pt.ComputedDoubleValues.STRESS_XZ,
            pt.ComputedDoubleValues.STRESS_YZ,
            pt.ComputedDoubleValues.STRESS_NXY,
            pt.ComputedDoubleValues.STRESS_NXZ,
            pt.ComputedDoubleValues.STRESS_NYZ,
        ]
        cfg.output_every = _GK_OUTPUT_EVERY
        cfg.filename = out_path
        dpd.config_step_output([cfg])
        dpd.run_simulation(n_steps=_N_STEPS_GK, dt=_DT, with_MC=with_mc)
        # Clear step output configuration before subsequent runs
        dpd.config_step_output([])
        data = np.loadtxt(out_path, skiprows=1)

    mean_temp = float(np.mean(data[:, 1]))
    mean_press = float(np.mean(data[:, 2]))

    corrs: List[pt.Correlator] = []
    for col_idx in range(3, 9):
        corr = pt.Correlator(24, 16, 2)
        corr.add_batch(data[:, col_idx])
        corr.evaluate(norm=norm_correlator)
        corrs.append(corr)

    n_pts = int(min(c.npcorr for c in corrs))
    t_lags = np.asarray(corrs[0].raw_t[:n_pts], dtype=float) * (_GK_OUTPUT_EVERY * _DT)
    shear_sum = np.sum([c.correlation[:n_pts] for c in corrs[:3]], axis=0)
    normal_sum = np.sum([c.correlation[:n_pts] for c in corrs[3:]], axis=0)
    g_gk = (vol / (5.0 * _KBT)) * (shear_sum + normal_sum / 6.0)
    return t_lags, g_gk, mean_temp, mean_press


def _compute_uniaxial_deformation_modulus(
    dpd: pt.DPDSimulator,
    universe_ref: pt.Universe,
    with_mc: bool,
    t_bin_edges: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Apply volume-preserving uniaxial ramp + hold deformation and measure G_deform(t)."""
    box_orig = universe_ref.get_box()
    inv_sqrt_lam = 1.0 / np.sqrt(_LAMBDA_STRETCH)
    target_box = pt.Box(
        box_orig.get_low_x() * inv_sqrt_lam,
        box_orig.get_high_x() * inv_sqrt_lam,
        box_orig.get_low_y() * inv_sqrt_lam,
        box_orig.get_high_y() * inv_sqrt_lam,
        box_orig.get_low_z() * _LAMBDA_STRETCH,
        box_orig.get_high_z() * _LAMBDA_STRETCH,
    )
    strain_factor = _LAMBDA_STRETCH**2 - 1.0 / _LAMBDA_STRETCH

    with tempfile.TemporaryDirectory() as tmpdir:
        hold_path = os.path.join(tmpdir, "deform_hold.txt")
        cfg = pt.OutputConfiguration()
        cfg.int_values = [pt.ComputedIntValues.STEP]
        cfg.double_values = [
            pt.ComputedDoubleValues.STRESS_XX,
            pt.ComputedDoubleValues.STRESS_YY,
            pt.ComputedDoubleValues.STRESS_ZZ,
        ]
        cfg.output_every = _DEF_OUTPUT_EVERY
        cfg.filename = hold_path

        dpd.config_box_deformation(target_box)
        dpd.run_simulation(n_steps=_N_STEPS_RAMP, dt=_DT, with_MC=with_mc)

        step_start_hold = dpd.get_current_timestep()
        dpd.config_step_output([cfg])
        dpd.run_simulation(n_steps=_N_STEPS_HOLD, dt=_DT, with_MC=with_mc)
        dpd.config_step_output([])
        data = np.loadtxt(hold_path, skiprows=1)

    t_hold = (data[:, 0] - float(step_start_hold)) * _DT
    # DPD virial stress uses pressure-positive convention (sigma_DPD = -sigma_Cauchy)
    sigma_true = 0.5 * (data[:, 1] + data[:, 2]) - data[:, 3]
    g_raw = sigma_true / strain_factor

    # Log-bin the relaxation curve to suppress high-frequency thermal DPD noise
    n_bins = len(t_bin_edges) - 1
    g_binned = np.empty(n_bins, dtype=float)
    t_centers = np.empty(n_bins, dtype=float)
    for b_idx in range(n_bins):
        if b_idx < n_bins - 1:
            mask = (t_hold >= t_bin_edges[b_idx]) & (t_hold < t_bin_edges[b_idx + 1])
        else:
            mask = (t_hold >= t_bin_edges[b_idx]) & (t_hold <= t_bin_edges[b_idx + 1])
        if np.any(mask):
            g_binned[b_idx] = float(np.mean(g_raw[mask]))
            t_centers[b_idx] = float(np.mean(t_hold[mask]))
        else:
            t_mid = 0.5 * (t_bin_edges[b_idx] + t_bin_edges[b_idx + 1])
            nearest = int(np.argmin(np.abs(t_hold - t_mid)))
            g_binned[b_idx] = float(g_raw[nearest])
            t_centers[b_idx] = float(t_hold[nearest])

    late_mask = t_hold >= 250.0
    g_late = float(np.mean(g_raw[late_mask]))
    return t_centers, g_binned, g_late


def _evaluate_nma_and_mehp(
    dpd: pt.DPDSimulator,
    vol: float,
    b2_measured: float,
    with_slipsprings: bool,
    t_nma_grid: np.ndarray,
) -> Tuple[np.ndarray, float, float]:
    """Evaluate static-graph NormalModeAnalyzer G_NMA(t) and MEHPSimulator G_MEHP."""
    u_snap = dpd.get_universe(with_slipsprings=with_slipsprings)

    nma = pt.NormalModeAnalyzer.from_universe(u_snap)
    spec = nma.solve_exact(
        solver=pt.NormalModeSolver.AUTO,
        compute_eigenvectors=False,
        decompose_components=True,
    )
    g_nma_relax = (_KBT / vol) * np.asarray(
        spec.stress_autocorrelation(t_nma_grid / _NMA_TAU0), dtype=float
    )

    mehp = pt.MEHPSimulator(u_snap, crosslinker_type=2)
    mehp.run()
    g_mehp_meas = (_KBT / vol) * float(np.sum(mehp.get_gamma_factors(b2_measured)))
    b2_ideal = 3.0 * _KBT / _SPRING_K
    g_mehp_ideal = (_KBT / vol) * float(np.sum(mehp.get_gamma_factors(b2_ideal)))
    return g_nma_relax, g_mehp_meas, g_mehp_ideal


def _run_single_seed(
    seed_idx: int,
    seed: int,
    t_bin_edges: np.ndarray,
    t_nma_grid: np.ndarray,
) -> Dict[str, Any]:
    """Execute all rheology conditions for one independent random seed."""
    u_melt_base, vol_melt = _build_melt_universe(seed=seed)
    u_net_base, vol_net = _build_network_universe(seed=seed)

    result: Dict[str, Any] = {"seed_idx": seed_idx, "seed": seed}
    for spec in _CONDITION_SPECS:
        key = spec.condition.value
        u_base = u_melt_base if spec.topology == TopologyKind.MELT else u_net_base
        vol = vol_melt if spec.topology == TopologyKind.MELT else vol_net

        dpd, b2_meas = _create_and_equilibrate_dpd(u_base, seed, spec)
        result[f"{key}_b2"] = b2_meas

        if spec.run_nma_mehp:
            g_nma_relax, g_mehp_meas, g_mehp_ideal = _evaluate_nma_and_mehp(
                dpd=dpd,
                vol=vol,
                b2_measured=b2_meas,
                with_slipsprings=spec.use_slipsprings,
                t_nma_grid=t_nma_grid,
            )
            result[f"{key}_g_nma_relax"] = g_nma_relax
            result[f"{key}_g_mehp"] = g_mehp_meas
            result[f"{key}_g_mehp_ideal"] = g_mehp_ideal

        if spec.run_green_kubo:
            # For melts, norm=False avoids finite-sample mean-subtraction bias;
            # for crosslinked networks, norm=True removes quenched shear stress.
            norm_corr = spec.topology == TopologyKind.NETWORK
            t_gk, g_gk, temp, press = _compute_green_kubo_modulus(
                dpd=dpd,
                vol=vol,
                with_mc=spec.with_mc_production,
                norm_correlator=norm_corr,
            )
            result[f"{key}_t_gk"] = t_gk
            result[f"{key}_g_gk"] = g_gk
            result[f"{key}_temp"] = temp
            result[f"{key}_press"] = press

        t_def, g_def, g_late = _compute_uniaxial_deformation_modulus(
            dpd=dpd,
            universe_ref=u_base,
            with_mc=spec.with_mc_production,
            t_bin_edges=t_bin_edges,
        )
        result[f"{key}_t_def"] = t_def
        result[f"{key}_g_def"] = g_def
        result[f"{key}_g_late"] = g_late

    return result


def _store_ensemble_stats(
    out: Dict[str, Any],
    prefix: str,
    samples: np.ndarray,
) -> None:
    """Compute ensemble mean and 95% Student's t CI and store into flat dict."""
    stats = compute_ensemble_statistics(samples, axis=0, confidence=0.95)
    out[f"{prefix}_mean"] = stats["mean"]
    out[f"{prefix}_std"] = stats["std"]
    out[f"{prefix}_ci"] = stats["ci"]
    out[f"{prefix}_ci_lower"] = stats["ci_lower"]
    out[f"{prefix}_ci_upper"] = stats["ci_upper"]


def run_dpd_rheology_study(max_workers: Optional[int] = None) -> Dict[str, Any]:
    """Run the multi-seed DPD slip-spring vs. network rheology study.

    Parameters
    ----------
    max_workers : int, optional
        Maximum worker processes for parallel seed execution.

    Returns
    -------
    dict
        Flat dictionary of NumPy arrays and metadata suitable for ``np.savez_compressed``.
    """
    t0 = time.perf_counter()
    seeds: Tuple[int, ...] = (101, 102, 103, 104, 105, 106)
    n_seeds = len(seeds)

    effective_workers = min(resolve_max_workers(max_workers), n_seeds)
    configure_process_threading(effective_workers, verbose=True)

    t_bin_edges = np.geomspace(1.0, _N_STEPS_HOLD * _DT, 24)
    t_nma_grid = np.logspace(-1, 3, 80)

    print(
        f"Running DPD rheology study ({n_seeds} seeds, {len(_CONDITION_SPECS)} conditions) "
        f"using {effective_workers} worker(s)..."
    )

    seed_results: List[Optional[Dict[str, Any]]] = [None] * n_seeds
    if effective_workers == 1:
        for idx, seed in enumerate(seeds):
            seed_results[idx] = _run_single_seed(idx, seed, t_bin_edges, t_nma_grid)
    else:
        with ProcessPoolExecutor(
            max_workers=effective_workers,
            initializer=worker_process_initializer,
        ) as executor:
            futures = {
                executor.submit(
                    _run_single_seed, idx, seed, t_bin_edges, t_nma_grid
                ): idx
                for idx, seed in enumerate(seeds)
            }
            for fut in as_completed(futures):
                res = fut.result()
                seed_results[int(res["seed_idx"])] = res

    completed: List[Dict[str, Any]] = [r for r in seed_results if r is not None]
    assert len(completed) == n_seeds

    out: Dict[str, Any] = {
        "n_seeds": np.array(n_seeds, dtype=int),
        "seeds": np.array(seeds, dtype=int),
        "rho_bead": np.array(_RHO_BEAD, dtype=float),
        "dpd_a": np.array(_DPD_A, dtype=float),
        "dpd_sigma": np.array(_DPD_SIGMA, dtype=float),
        "spring_k": np.array(_SPRING_K, dtype=float),
        "kbt": np.array(_KBT, dtype=float),
        "dt": np.array(_DT, dtype=float),
        "ss_fraction": np.array(_SS_FRACTION, dtype=float),
        "n_chains_melt": np.array(_N_CHAINS_MELT, dtype=int),
        "n_melt": np.array(_N_MELT, dtype=int),
        "n_xlinks": np.array(_N_XLINKS, dtype=int),
        "n_strand_net": np.array(_N_STRAND_NET, dtype=int),
        "conversion_p": np.array(_CONVERSION_P, dtype=float),
        "lambda_stretch": np.array(_LAMBDA_STRETCH, dtype=float),
        "nma_tau0": np.array(_NMA_TAU0, dtype=float),
        "t_nma": t_nma_grid,
        "t_def": np.mean([r["melt_no_ss_t_def"] for r in completed], axis=0),
    }

    # Find common Green-Kubo lag grid length across all seeds and GK conditions
    min_npcorr = min(
        len(r[f"{spec.condition.value}_t_gk"])
        for r in completed
        for spec in _CONDITION_SPECS
        if spec.run_green_kubo
    )
    out["t_gk"] = completed[0]["melt_no_ss_t_gk"][:min_npcorr]

    for spec in _CONDITION_SPECS:
        key = spec.condition.value
        _store_ensemble_stats(
            out, f"{key}_b2", np.array([r[f"{key}_b2"] for r in completed])
        )
        _store_ensemble_stats(
            out, f"{key}_g_def", np.vstack([r[f"{key}_g_def"] for r in completed])
        )
        _store_ensemble_stats(
            out, f"{key}_g_late", np.array([r[f"{key}_g_late"] for r in completed])
        )

        if spec.run_green_kubo:
            gk_stack = np.vstack([r[f"{key}_g_gk"][:min_npcorr] for r in completed])
            _store_ensemble_stats(out, f"{key}_g_gk", gk_stack)
            _store_ensemble_stats(
                out, f"{key}_temp", np.array([r[f"{key}_temp"] for r in completed])
            )
            _store_ensemble_stats(
                out, f"{key}_press", np.array([r[f"{key}_press"] for r in completed])
            )

        if spec.run_nma_mehp:
            nma_relax_stack = np.vstack([r[f"{key}_g_nma_relax"] for r in completed])
            g_mehp_arr = np.array([r[f"{key}_g_mehp"] for r in completed])
            g_mehp_ideal_arr = np.array([r[f"{key}_g_mehp_ideal"] for r in completed])
            nma_total_stack = nma_relax_stack + g_mehp_arr[:, None]

            _store_ensemble_stats(out, f"{key}_g_nma_relax", nma_relax_stack)
            _store_ensemble_stats(out, f"{key}_g_nma_total", nma_total_stack)
            _store_ensemble_stats(out, f"{key}_g_mehp", g_mehp_arr)
            _store_ensemble_stats(out, f"{key}_g_mehp_ideal", g_mehp_ideal_arr)

    elapsed = time.perf_counter() - t0
    out["elapsed_seconds"] = np.array(elapsed, dtype=float)
    print(f"Completed DPD rheology study in {elapsed:.1f}s.")
    return out


if __name__ == "__main__":
    cache_path = get_example_cache_path(
        "dpd_slipspring_rheology_cache.npz", "dpd_simulations"
    )
    load_or_compute_cache(cache_path, run_dpd_rheology_study)
