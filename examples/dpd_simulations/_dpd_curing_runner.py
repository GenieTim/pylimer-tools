"""
Simulation runner for Study 4b: In-Situ DPD Dynamic Crosslinking vs. MCUniverseGenerator.

Provides parallelized, cache-friendly study functions comparing:
1. In-situ dynamic curing kinetics p(t) of a stoichiometric f=4 precursor melt
   in DPDSimulator (using config_bond_formation) against instantaneous
   MCUniverseGenerator networks at matched (N, p, r, rho, <b^2>).
2. Topological loop spectra (orders 1-4 via Universe.count_loop_lengths),
   spatial crosslink homogeneity (sub-box Fano factor and crosslink-crosslink
   radial distribution function g_xx(r)), and network strand end-to-end
   distance distributions P(R_ee).
3. Mechanical equilibrium shear moduli (MEHPSimulator with EntanglementMode.NONE
   and EntanglementMode.MERGED_LINKS) and active/dangling/soluble weight fractions.
4. OpenMP multi-threaded RNG non-determinism benchmark (OMP_NUM_THREADS=1 vs 2
   with identical initial configuration and random seed).

This module is prefixed with an underscore so Sphinx-Gallery ignores it when
building documentation.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from enum import Enum
import os
from pathlib import Path
import subprocess
import sys
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


class SynthesisMethod(Enum):
    """Network synthesis pathway compared in this study."""

    DPD_CURING = "dpd"
    MC_GENERATOR = "mc"


class MEHPEntanglementRegime(Enum):
    """MEHP entanglement evaluation regime."""

    PHANTOM = "phantom"
    ENTANGLED = "entangled"


@dataclass(frozen=True)
class _CuringStudyConfig:
    """Physical and numerical parameters for the DPD curing vs. MC study."""

    n_xlinks: int = 100
    f_xlink: int = 4
    strand_length: int = 20
    rho_bead: float = 3.0
    dpd_a: float = 25.0
    dpd_sigma: float = 3.0
    spring_k: float = 2.0
    kbt: float = 1.0
    dt: float = 0.05
    n_steps_equil: int = 2000
    bond_formation_dist: float = 0.85
    attempt_every: int = 25
    t_cure_max: float = 120.0
    target_conversions: Tuple[float, ...] = (0.55, 0.70, 0.80, 0.88, 0.94)
    entanglement_fraction: float = 0.05
    ent_upper_cutoff: float = 2.0
    ent_lower_cutoff: float = 0.5
    ent_same_strand_cutoff: float = 2.0
    n_subbox_grid: int = 3
    max_loop_order: int = 4
    n_rdf_bins: int = 25
    n_ree_bins: int = 20
    ree_hist_max: float = 14.0
    seeds: Tuple[int, ...] = (101, 102, 103, 104, 105, 106)

    @property
    def n_strands(self) -> int:
        """Number of bifunctional strands at stoichiometric balance r = 1.0."""
        return (self.f_xlink * self.n_xlinks) // 2

    @property
    def n_atoms(self) -> int:
        """Total number of beads (crosslinkers + strand beads)."""
        return self.n_xlinks + self.n_strands * self.strand_length

    @property
    def max_new_bonds(self) -> int:
        """Maximum number of crosslink-strand bonds at p = 1.0."""
        return self.f_xlink * self.n_xlinks

    @property
    def volume(self) -> float:
        """Simulation box volume at target bead density rho."""
        return float(self.n_atoms) / self.rho_bead

    @property
    def box_l(self) -> float:
        """Cubic simulation box edge length L."""
        return self.volume ** (1.0 / 3.0)


_DEFAULT_CONFIG: _CuringStudyConfig = _CuringStudyConfig()


def _build_precursor_universe(seed: int, cfg: _CuringStudyConfig) -> pt.Universe:
    """Create an uncrosslinked stoichiometric precursor melt of strands and crosslinkers."""
    b2_ideal = 3.0 * cfg.kbt / cfg.spring_k
    b_ideal = float(np.sqrt(b2_ideal))

    gen = pt.MCUniverseGenerator(cfg.box_l, cfg.box_l, cfg.box_l)
    gen.set_seed(int(seed))
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(b_ideal)
    gen.set_mean_squared_bead_distance(b2_ideal)
    gen.add_crosslinkers(cfg.n_xlinks, cfg.f_xlink, crosslinker_type=2)
    gen.add_strands(
        cfg.n_strands, [cfg.strand_length] * cfg.n_strands, strand_atom_type=1
    )
    return gen.get_universe()


def _build_matched_mc_universe(
    seed: int,
    conversion: float,
    b2_measured: float,
    cfg: _CuringStudyConfig,
) -> pt.Universe:
    """Generate a matched MCUniverseGenerator network at given conversion and <b^2>."""
    gen = pt.MCUniverseGenerator(cfg.box_l, cfg.box_l, cfg.box_l)
    gen.set_seed(int(seed))
    gen.config_nr_of_mc_steps(0)
    gen.config_primary_loop_probability(1.0)
    gen.config_secondary_loop_probability(1.0)
    gen.set_bead_distance(float(np.sqrt(b2_measured)))
    gen.set_mean_squared_bead_distance(b2_measured)
    gen.add_crosslinkers(cfg.n_xlinks, cfg.f_xlink, crosslinker_type=2)
    gen.add_strands(
        cfg.n_strands, [cfg.strand_length] * cfg.n_strands, strand_atom_type=1
    )
    gen.link_strands_to_conversion(float(conversion))
    return gen.get_universe()


def _compute_crosslink_spatial_metrics(
    universe: pt.Universe,
    rdf_edges: np.ndarray,
    n_grid: int,
) -> Tuple[float, np.ndarray, float]:
    """Compute crosslink sub-box Fano factor, g_xx(r), and mean nearest-neighbor distance."""
    x_coords = np.asarray(
        [
            atom.get_coordinates()
            for atom in universe.get_atoms()
            if atom.get_type() == 2
        ],
        dtype=float,
    )
    box = universe.get_box()
    lx, ly, lz = float(box.get_lx()), float(box.get_ly()), float(box.get_lz())
    vol = lx * ly * lz

    wrapped = np.column_stack(
        [
            np.mod(x_coords[:, 0] - box.get_low_x(), lx),
            np.mod(x_coords[:, 1] - box.get_low_y(), ly),
            np.mod(x_coords[:, 2] - box.get_low_z(), lz),
        ]
    )

    # Sub-box Fano factor: variance / mean of crosslink counts per cell
    ix = np.minimum((wrapped[:, 0] / lx * n_grid).astype(int), n_grid - 1)
    iy = np.minimum((wrapped[:, 1] / ly * n_grid).astype(int), n_grid - 1)
    iz = np.minimum((wrapped[:, 2] / lz * n_grid).astype(int), n_grid - 1)
    counts = np.zeros((n_grid, n_grid, n_grid), dtype=float)
    for a, b, c in zip(ix, iy, iz):
        counts[a, b, c] += 1.0
    fano = float(np.var(counts, ddof=1) / np.mean(counts))

    # Minimum-image crosslink-crosslink pair distances
    diff = wrapped[:, None, :] - wrapped[None, :, :]
    diff[:, :, 0] -= lx * np.round(diff[:, :, 0] / lx)
    diff[:, :, 1] -= ly * np.round(diff[:, :, 1] / ly)
    diff[:, :, 2] -= lz * np.round(diff[:, :, 2] / lz)
    dists = np.linalg.norm(diff, axis=-1)

    n_x = len(wrapped)
    triu = dists[np.triu_indices(n_x, k=1)]
    hist, _ = np.histogram(triu, bins=rdf_edges)
    shell_vol = (4.0 / 3.0) * np.pi * (rdf_edges[1:] ** 3 - rdf_edges[:-1] ** 3)
    rho_x = float(n_x) / vol
    g_xx = (2.0 * hist.astype(float)) / (float(n_x) * rho_x * shell_vol)

    # Mean nearest-neighbor crosslink distance
    np.fill_diagonal(dists, np.inf)
    mean_nn = float(np.mean(np.min(dists, axis=1)))
    return fano, g_xx, mean_nn


def _analyze_network_snapshot(
    universe: pt.Universe,
    b2_ref: float,
    ent_seed: str,
    cfg: _CuringStudyConfig,
    rdf_edges: np.ndarray,
    ree_edges: np.ndarray,
) -> Dict[str, Any]:
    """Evaluate topological, spatial, conformational, and MEHP metrics on a snapshot."""
    box = universe.get_box()
    vol = float(box.get_lx() * box.get_ly() * box.get_lz())

    # 1. Loop spectrum up to max_loop_order
    stride = cfg.strand_length + 1
    max_atoms = cfg.max_loop_order * stride + 2
    loop_counts_atoms = universe.count_loop_lengths(max_length=max_atoms)
    loop_counts = np.array(
        [
            int(loop_counts_atoms.get(order * stride, 0))
            for order in range(1, cfg.max_loop_order + 1)
        ],
        dtype=float,
    )

    # 2. Spatial crosslink homogeneity
    fano, g_xx, mean_nn = _compute_crosslink_spatial_metrics(
        universe, rdf_edges=rdf_edges, n_grid=cfg.n_subbox_grid
    )

    # 3. Network strand end-to-end distance distribution
    chains = universe.get_chains_with_crosslinker(2)
    ree_list: List[float] = [
        float(c.compute_end_to_end_distance_with_derived_image_flags())
        for c in chains
        if c.get_strand_type() == pt.MoleculeType.NETWORK_STRAND
    ]
    ree_arr = np.asarray(ree_list, dtype=float)
    if len(ree_arr) > 0:
        ree_rms = float(np.sqrt(np.mean(ree_arr**2)))
        ree_hist, _ = np.histogram(ree_arr, bins=ree_edges, density=True)
    else:
        ree_rms = float("nan")
        ree_hist = np.zeros(len(ree_edges) - 1, dtype=float)

    # 4. MEHP phantom (EntanglementMode.NONE)
    mehp_none = pt.MEHPSimulator(
        universe,
        crosslinker_type=2,
        entanglement_mode=pt.EntanglementMode.NONE,
    )
    mehp_none.config_mean_squared_bond_length(b2_ref)
    mehp_none.run()
    g_phantom = (cfg.kbt / vol) * float(np.sum(mehp_none.get_gamma_factors(b2_ref)))
    w_active = float(mehp_none.get_active_weight_fraction())
    w_dangling = float(mehp_none.get_dangling_weight_fraction())
    w_soluble = float(mehp_none.get_soluble_weight_fraction())

    # 5. MEHP entangled (EntanglementMode.MERGED_LINKS)
    n_ent = int(round(cfg.entanglement_fraction * universe.get_nr_of_atoms()))
    mehp_ent = pt.MEHPSimulator.construct_with_random_entanglements(
        universe=universe,
        nr_of_entanglements_to_sample=n_ent,
        upper_sampling_cutoff=cfg.ent_upper_cutoff,
        lower_sampling_cutoff=cfg.ent_lower_cutoff,
        same_strand_cutoff=cfg.ent_same_strand_cutoff,
        seed=ent_seed,
        crosslinker_type=2,
        entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
    )
    mehp_ent.config_mean_squared_bond_length(b2_ref)
    mehp_ent.run()
    g_entangled = (cfg.kbt / vol) * float(np.sum(mehp_ent.get_gamma_factors(b2_ref)))
    w_active_ent = float(mehp_ent.get_active_weight_fraction())

    return {
        "loops": loop_counts,
        "fano": fano,
        "g_xx": g_xx,
        "mean_nn": mean_nn,
        "ree_rms": ree_rms,
        "ree_hist": ree_hist,
        "g_phantom": g_phantom,
        "g_entangled": g_entangled,
        "w_active": w_active,
        "w_dangling": w_dangling,
        "w_soluble": w_soluble,
        "w_active_ent": w_active_ent,
    }


def _run_single_seed(
    seed_idx: int,
    seed: int,
    cfg: _CuringStudyConfig,
    t_cure_grid: np.ndarray,
    rdf_edges: np.ndarray,
    ree_edges: np.ndarray,
) -> Dict[str, Any]:
    """Run one complete DPD curing trajectory and matched MC comparison."""
    u0 = _build_precursor_universe(seed=seed, cfg=cfg)
    base_bonds = u0.get_nr_of_bonds()

    dpd = pt.DPDSimulator(
        u0,
        crosslinker_type=2,
        seed=f"dpd_curing_seed_{seed}",
    )
    dpd.assume_box_large_enough()
    dpd.config_a(cfg.dpd_a)
    dpd.config_sigma(cfg.dpd_sigma)
    dpd.config_spring_constant(cfg.spring_k)
    dpd.run_simulation(n_steps=cfg.n_steps_equil, dt=cfg.dt, with_MC=False)

    # Measure equilibrated DPD melt <b^2> and free-chain swelling C_inf_eff
    b2_meas = float(np.mean(dpd.get_bond_lengths()[:base_bonds] ** 2))
    u_melt = dpd.get_universe(with_slipsprings=False)
    ree_melt = np.asarray(u_melt.compute_end_to_end_distances(2), dtype=float)
    c_inf_eff = float(np.mean(ree_melt**2) / ((cfg.strand_length - 1) * b2_meas))

    # Staged dynamic curing to capture exact target conversions while recording p(t)
    t_step_start = dpd.get_current_timestep()
    kin_t: List[float] = [0.0]
    kin_p: List[float] = [0.0]
    dpd_snapshots: List[pt.Universe] = []
    actual_p_list: List[float] = []
    snapshot_t_list: List[float] = []

    max_cure_steps = int(round(cfg.t_cure_max / cfg.dt))
    for p_tgt in cfg.target_conversions:
        target_total = int(round(p_tgt * cfg.max_new_bonds))
        already = dpd.get_nr_of_bonds() - base_bonds
        remaining = max(0, target_total - already)
        dpd.config_bond_formation(
            num_bonds_to_form=remaining,
            max_bonds_per_atom_type={1: 2, 2: cfg.f_xlink},
            bond_formation_dist=cfg.bond_formation_dist,
            attempt_bond_formation_every=cfg.attempt_every,
            atom_type_form_from=2,
            atom_type_form_to=1,
        )
        while (dpd.get_nr_of_bonds() - base_bonds) < target_total and (
            dpd.get_current_timestep() - t_step_start
        ) < max_cure_steps * 2:
            elapsed_steps = dpd.get_current_timestep() - t_step_start
            step_inc = 25 if elapsed_steps < 400 else 50
            dpd.run_simulation(n_steps=step_inc, dt=cfg.dt, with_MC=False)
            formed = dpd.get_nr_of_bonds() - base_bonds
            kin_t.append((dpd.get_current_timestep() - t_step_start) * cfg.dt)
            kin_p.append(float(formed) / float(cfg.max_new_bonds))

        formed_now = dpd.get_nr_of_bonds() - base_bonds
        actual_p_list.append(float(formed_now) / float(cfg.max_new_bonds))
        snapshot_t_list.append((dpd.get_current_timestep() - t_step_start) * cfg.dt)
        dpd_snapshots.append(dpd.get_universe(with_slipsprings=False))

    # Continue curing up to t_cure_max so p(t) covers the full kinetics window
    already = dpd.get_nr_of_bonds() - base_bonds
    if already < cfg.max_new_bonds:
        dpd.config_bond_formation(
            num_bonds_to_form=cfg.max_new_bonds - already,
            max_bonds_per_atom_type={1: 2, 2: cfg.f_xlink},
            bond_formation_dist=cfg.bond_formation_dist,
            attempt_bond_formation_every=cfg.attempt_every,
            atom_type_form_from=2,
            atom_type_form_to=1,
        )
    while (dpd.get_current_timestep() - t_step_start) < max_cure_steps:
        dpd.run_simulation(n_steps=100, dt=cfg.dt, with_MC=False)
        formed = dpd.get_nr_of_bonds() - base_bonds
        kin_t.append((dpd.get_current_timestep() - t_step_start) * cfg.dt)
        kin_p.append(float(formed) / float(cfg.max_new_bonds))

    p_interp = np.interp(t_cure_grid, np.asarray(kin_t), np.asarray(kin_p))

    n_conv = len(cfg.target_conversions)
    result: Dict[str, Any] = {
        "seed_idx": seed_idx,
        "seed": seed,
        "b2_meas": b2_meas,
        "c_inf_eff": c_inf_eff,
        "p_kinetics": p_interp,
        "actual_p": np.asarray(actual_p_list, dtype=float),
        "snapshot_t": np.asarray(snapshot_t_list, dtype=float),
    }

    for method in (SynthesisMethod.DPD_CURING, SynthesisMethod.MC_GENERATOR):
        key = method.value
        loops_arr = np.zeros((n_conv, cfg.max_loop_order), dtype=float)
        fano_arr = np.zeros(n_conv, dtype=float)
        mean_nn_arr = np.zeros(n_conv, dtype=float)
        ree_rms_arr = np.zeros(n_conv, dtype=float)
        g_phantom_arr = np.zeros(n_conv, dtype=float)
        g_entangled_arr = np.zeros(n_conv, dtype=float)
        w_active_arr = np.zeros(n_conv, dtype=float)
        w_dangling_arr = np.zeros(n_conv, dtype=float)
        w_soluble_arr = np.zeros(n_conv, dtype=float)
        w_active_ent_arr = np.zeros(n_conv, dtype=float)
        g_xx_final = np.zeros(cfg.n_rdf_bins, dtype=float)
        ree_hist_final = np.zeros(cfg.n_ree_bins, dtype=float)

        for c_idx, p_actual in enumerate(actual_p_list):
            if method == SynthesisMethod.DPD_CURING:
                u_snap = dpd_snapshots[c_idx]
            else:
                u_snap = _build_matched_mc_universe(
                    seed=seed * 100 + c_idx,
                    conversion=p_actual,
                    b2_measured=b2_meas,
                    cfg=cfg,
                )

            metrics = _analyze_network_snapshot(
                universe=u_snap,
                b2_ref=b2_meas,
                ent_seed=f"ent_{key}_{seed}_{c_idx}",
                cfg=cfg,
                rdf_edges=rdf_edges,
                ree_edges=ree_edges,
            )
            loops_arr[c_idx] = metrics["loops"]
            fano_arr[c_idx] = metrics["fano"]
            mean_nn_arr[c_idx] = metrics["mean_nn"]
            ree_rms_arr[c_idx] = metrics["ree_rms"]
            g_phantom_arr[c_idx] = metrics["g_phantom"]
            g_entangled_arr[c_idx] = metrics["g_entangled"]
            w_active_arr[c_idx] = metrics["w_active"]
            w_dangling_arr[c_idx] = metrics["w_dangling"]
            w_soluble_arr[c_idx] = metrics["w_soluble"]
            w_active_ent_arr[c_idx] = metrics["w_active_ent"]
            if c_idx == n_conv - 1:
                g_xx_final = metrics["g_xx"]
                ree_hist_final = metrics["ree_hist"]

        result[f"{key}_loops"] = loops_arr
        result[f"{key}_fano"] = fano_arr
        result[f"{key}_mean_nn"] = mean_nn_arr
        result[f"{key}_ree_rms"] = ree_rms_arr
        result[f"{key}_g_phantom"] = g_phantom_arr
        result[f"{key}_g_entangled"] = g_entangled_arr
        result[f"{key}_w_active"] = w_active_arr
        result[f"{key}_w_dangling"] = w_dangling_arr
        result[f"{key}_w_soluble"] = w_soluble_arr
        result[f"{key}_w_active_ent"] = w_active_ent_arr
        result[f"{key}_g_xx_final"] = g_xx_final
        result[f"{key}_ree_hist_final"] = ree_hist_final

    return result


def _run_openmp_determinism_check() -> Dict[str, np.ndarray]:
    """Quantify DPD reproducibility under OMP_NUM_THREADS=1 vs. OMP_NUM_THREADS=2."""
    snippet = (
        "import json, numpy as np, pylimer_tools_cpp as pt\n"
        "gen = pt.MCUniverseGenerator(7.0, 7.0, 7.0)\n"
        "gen.set_seed(77)\n"
        "gen.config_nr_of_mc_steps(0)\n"
        "gen.set_bead_distance(float(np.sqrt(1.5)))\n"
        "gen.set_mean_squared_bead_distance(1.5)\n"
        "gen.add_crosslinkers(25, 4, crosslinker_type=2)\n"
        "gen.add_strands(50, [20] * 50, strand_atom_type=1)\n"
        "u0 = gen.get_universe()\n"
        "dpd = pt.DPDSimulator(u0, crosslinker_type=2, seed='omp_determinism_check')\n"
        "dpd.assume_box_large_enough()\n"
        "dpd.config_a(25.0)\n"
        "dpd.config_sigma(3.0)\n"
        "dpd.config_spring_constant(2.0)\n"
        "dpd.config_bond_formation(100, {1: 2, 2: 4}, 0.85, 25, 2, 1)\n"
        "dpd.run_simulation(n_steps=400, dt=0.05, with_MC=False)\n"
        "coords = np.asarray(dpd.get_coordinates(), dtype=float).tolist()\n"
        "payload = {'bonds': int(dpd.get_nr_of_bonds()), 'coords': coords}\n"
        "print('JSON_PAYLOAD:' + json.dumps(payload))\n"
    )

    import json

    def _run_one(n_threads: int) -> Tuple[int, np.ndarray]:
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(n_threads)
        proc = subprocess.run(
            [sys.executable, "-c", snippet],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        for line in proc.stdout.splitlines():
            if line.startswith("JSON_PAYLOAD:"):
                data = json.loads(line.split("JSON_PAYLOAD:", 1)[1])
                return int(data["bonds"]), np.asarray(data["coords"], dtype=float)
        raise RuntimeError("Failed to parse OpenMP determinism subprocess output.")

    b1a, c1a = _run_one(1)
    b1b, c1b = _run_one(1)
    b2a, c2a = _run_one(2)
    b2b, c2b = _run_one(2)

    rmsd_1t = float(np.sqrt(np.mean((c1a - c1b) ** 2)))
    rmsd_2t = float(np.sqrt(np.mean((c2a - c2b) ** 2)))
    rmsd_1t_vs_2t = float(np.sqrt(np.mean((c1a - c2a) ** 2)))

    return {
        "omp_1t_bonds": np.array([b1a, b1b], dtype=int),
        "omp_2t_bonds": np.array([b2a, b2b], dtype=int),
        "omp_1t_coord_rmsd": np.array(rmsd_1t, dtype=float),
        "omp_2t_coord_rmsd": np.array(rmsd_2t, dtype=float),
        "omp_1t_vs_2t_coord_rmsd": np.array(rmsd_1t_vs_2t, dtype=float),
    }


def _store_ensemble_stats(
    out: Dict[str, Any],
    prefix: str,
    samples: np.ndarray,
) -> None:
    """Compute ensemble mean, std, and 95% Student's t CI and store into flat dict."""
    stats = compute_ensemble_statistics(samples, axis=0, confidence=0.95)
    out[f"{prefix}_mean"] = stats["mean"]
    out[f"{prefix}_std"] = stats["std"]
    out[f"{prefix}_ci"] = stats["ci"]
    out[f"{prefix}_ci_lower"] = stats["ci_lower"]
    out[f"{prefix}_ci_upper"] = stats["ci_upper"]


def run_dpd_curing_study(max_workers: Optional[int] = None) -> Dict[str, Any]:
    """Run the multi-seed DPD dynamic crosslinking vs. MCUniverseGenerator study.

    Parameters
    ----------
    max_workers : int, optional
        Maximum number of parallel worker processes across seeds.

    Returns
    -------
    dict
        Flat dictionary of NumPy arrays and metadata suitable for ``np.savez_compressed``.
    """
    t0 = time.perf_counter()
    cfg = _DEFAULT_CONFIG
    n_seeds = len(cfg.seeds)
    effective_workers = min(resolve_max_workers(max_workers), n_seeds)
    configure_process_threading(effective_workers, verbose=True)

    t_cure_grid = np.linspace(0.0, cfg.t_cure_max, 61)
    rdf_edges = np.linspace(0.0, 0.5 * cfg.box_l, cfg.n_rdf_bins + 1)
    rdf_centers = 0.5 * (rdf_edges[:-1] + rdf_edges[1:])
    ree_edges = np.linspace(0.0, cfg.ree_hist_max, cfg.n_ree_bins + 1)
    ree_centers = 0.5 * (ree_edges[:-1] + ree_edges[1:])

    print(
        f"Running DPD dynamic curing vs. MC study ({n_seeds} seeds, "
        f"{cfg.n_strands} strands x N={cfg.strand_length}, L={cfg.box_l:.2f}) "
        f"using {effective_workers} worker(s)..."
    )

    seed_results: List[Optional[Dict[str, Any]]] = [None] * n_seeds
    if effective_workers == 1:
        for idx, seed in enumerate(cfg.seeds):
            seed_results[idx] = _run_single_seed(
                idx, seed, cfg, t_cure_grid, rdf_edges, ree_edges
            )
    else:
        with ProcessPoolExecutor(
            max_workers=effective_workers,
            initializer=worker_process_initializer,
        ) as executor:
            futures = {
                executor.submit(
                    _run_single_seed,
                    idx,
                    seed,
                    cfg,
                    t_cure_grid,
                    rdf_edges,
                    ree_edges,
                ): idx
                for idx, seed in enumerate(cfg.seeds)
            }
            for fut in as_completed(futures):
                res = fut.result()
                seed_results[int(res["seed_idx"])] = res

    completed: List[Dict[str, Any]] = [r for r in seed_results if r is not None]
    assert len(completed) == n_seeds

    print("Running OpenMP 1-thread vs. 2-thread determinism check...")
    omp_check = _run_openmp_determinism_check()

    out: Dict[str, Any] = {
        "n_seeds": np.array(n_seeds, dtype=int),
        "seeds": np.array(cfg.seeds, dtype=int),
        "n_xlinks": np.array(cfg.n_xlinks, dtype=int),
        "f_xlink": np.array(cfg.f_xlink, dtype=int),
        "n_strands": np.array(cfg.n_strands, dtype=int),
        "strand_length": np.array(cfg.strand_length, dtype=int),
        "rho_bead": np.array(cfg.rho_bead, dtype=float),
        "box_l": np.array(cfg.box_l, dtype=float),
        "volume": np.array(cfg.volume, dtype=float),
        "dpd_a": np.array(cfg.dpd_a, dtype=float),
        "dpd_sigma": np.array(cfg.dpd_sigma, dtype=float),
        "spring_k": np.array(cfg.spring_k, dtype=float),
        "dt": np.array(cfg.dt, dtype=float),
        "bond_formation_dist": np.array(cfg.bond_formation_dist, dtype=float),
        "attempt_every": np.array(cfg.attempt_every, dtype=int),
        "entanglement_fraction": np.array(cfg.entanglement_fraction, dtype=float),
        "target_conversions": np.array(cfg.target_conversions, dtype=float),
        "loop_orders": np.arange(1, cfg.max_loop_order + 1, dtype=int),
        "t_cure_grid": t_cure_grid,
        "rdf_centers": rdf_centers,
        "ree_centers": ree_centers,
    }
    out.update(omp_check)

    _store_ensemble_stats(
        out, "b2_meas", np.array([r["b2_meas"] for r in completed], dtype=float)
    )
    _store_ensemble_stats(
        out, "c_inf_eff", np.array([r["c_inf_eff"] for r in completed], dtype=float)
    )
    _store_ensemble_stats(
        out, "p_kinetics", np.vstack([r["p_kinetics"] for r in completed])
    )
    _store_ensemble_stats(
        out, "actual_p", np.vstack([r["actual_p"] for r in completed])
    )
    _store_ensemble_stats(
        out, "snapshot_t", np.vstack([r["snapshot_t"] for r in completed])
    )

    for method in (SynthesisMethod.DPD_CURING, SynthesisMethod.MC_GENERATOR):
        key = method.value
        for field in (
            "loops",
            "fano",
            "mean_nn",
            "ree_rms",
            "g_phantom",
            "g_entangled",
            "w_active",
            "w_dangling",
            "w_soluble",
            "w_active_ent",
            "g_xx_final",
            "ree_hist_final",
        ):
            stacked = np.stack([r[f"{key}_{field}"] for r in completed], axis=0)
            _store_ensemble_stats(out, f"{key}_{field}", stacked)

    elapsed = time.perf_counter() - t0
    out["elapsed_seconds"] = np.array(elapsed, dtype=float)
    print(f"Completed DPD curing vs. MC study in {elapsed:.1f}s.")
    return out


if __name__ == "__main__":
    cache_path = get_example_cache_path(
        "dpd_dynamic_curing_vs_mc_cache.npz", "dpd_simulations"
    )
    load_or_compute_cache(cache_path, run_dpd_curing_study)
