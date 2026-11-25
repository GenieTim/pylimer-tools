"""
Simulation runner for Study 5: Algorithmic & PBC Failure Modes
(``MaxDistanceProvider`` bias and ``assume_box_large_enough`` collapse).

This private module (prefixed with ``_`` so Sphinx-Gallery ignores it) runs
parallel multi-seed ensembles for:

1. **MC Linking Cutoff Sweep (``MaxDistanceProvider``)**:
   Evaluates ``use_zscore_max_distance(z, in_sqrt_multiplier=b2)``,
   ``use_linear_max_distance(m)``, and ``disable_max_distance()`` on melt PDMS
   networks (``ParameterType.GAUSSIAN``), measuring wall time, achieved
   conversion, global-search fallback count, end-to-end distance distribution
   :math:`P(R)`, primary/secondary loop fractions, and phantom + entangled
   MEHP shear moduli.
2. **Minimum-Image Collapse Sweep (``assume_box_large_enough``)**:
   Builds multi-image periodic Gaussian networks across box sizes
   :math:`L / \\sqrt{(N+1)b^2} \\in [0.8, 3.0]` and compares
   ``MEHPSimulator(..., assume_box_large_enough=False)`` against
   ``assume_box_large_enough=True`` on spring length distributions, flipped
   spring fractions, and shear modulus collapse.
3. **2D Visual Projection Demo**:
   Constructs a small periodic 2D network and stores relaxed node positions and
   spring vectors under both PBC modes to visualize periodic bond flipping and
   box-spanning collapse directly.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

if "__file__" in globals():
    _examples_dir = str(Path(__file__).resolve().parents[1])
    if _examples_dir not in sys.path:
        sys.path.insert(0, _examples_dir)

from _example_cache import (  # noqa: E402
    compute_ensemble_statistics,
    configure_process_threading,
    get_example_cache_path,
    load_or_compute_cache,
    resolve_max_workers,
    worker_process_initializer,
)
from pylimer_tools.io.bead_spring_parameter_provider import (  # noqa: E402
    ParameterType,
    Parameters,
    get_parameters_for_polymer,
)
import pylimer_tools_cpp as pt  # noqa: E402


class MaxDistanceMode(Enum):
    """Supported ``MaxDistanceProvider`` strategies in ``MCUniverseGenerator``."""

    ZSCORE = "zscore"
    LINEAR = "linear"
    DISABLED = "disabled"


class PBCMode(Enum):
    """Periodic boundary condition handling modes in ``MEHPSimulator``."""

    BAKED_OFFSETS = "baked_offsets"
    MINIMUM_IMAGE = "minimum_image"


@dataclass(frozen=True)
class CutoffConfig:
    """Single ``MaxDistanceProvider`` configuration in the MC cutoff sweep."""

    label: str
    mode: MaxDistanceMode
    value: float


@dataclass(frozen=True)
class PDMSParams:
    """Extracted scalar PDMS Gaussian bead-spring parameters in simulation units."""

    b: float
    b2: float
    r02: float
    bead_density: float
    entanglement_density: float
    sampling_cutoff: float
    gamma_to_mpa: float


# Study configuration constants
MC_N_BEADS: int = 30
MC_N_STRANDS: int = 26000
MC_TARGET_CONVERSION: float = 0.95
MC_SEEDS: Tuple[int, ...] = (101, 202, 303, 404, 505, 606)
MC_R_HIST_MAX: float = 3.5
MC_R_HIST_BINS: int = 70

CUTOFF_CONFIGS: Tuple[CutoffConfig, ...] = (
    CutoffConfig("z=1.0", MaxDistanceMode.ZSCORE, 1.0),
    CutoffConfig("z=1.5", MaxDistanceMode.ZSCORE, 1.5),
    CutoffConfig("z=2.0", MaxDistanceMode.ZSCORE, 2.0),
    CutoffConfig("z=2.5", MaxDistanceMode.ZSCORE, 2.5),
    CutoffConfig("z=3.0", MaxDistanceMode.ZSCORE, 3.0),
    CutoffConfig("z=3.5", MaxDistanceMode.ZSCORE, 3.5),
    CutoffConfig("z=4.0", MaxDistanceMode.ZSCORE, 4.0),
    CutoffConfig("lin=0.25", MaxDistanceMode.LINEAR, 0.25),
    CutoffConfig("lin=0.60", MaxDistanceMode.LINEAR, 0.60),
    CutoffConfig("disabled", MaxDistanceMode.DISABLED, -1.0),
)

PBC_N_BEADS: int = 50
PBC_TARGET_CONVERSION: float = 0.95
PBC_ALPHAS: Tuple[float, ...] = (0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0)
PBC_SEEDS: Tuple[int, ...] = tuple(2000 + i * 37 for i in range(16))
PBC_HIST_ALPHA_INDEX: int = 1  # alpha = 1.0
PBC_RL_HIST_MAX: float = 1.6
PBC_RL_HIST_BINS: int = 64


def extract_pdms_parameters() -> PDMSParams:
    """Return PDMS ``ParameterType.GAUSSIAN`` scalar parameters."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    assert isinstance(params, Parameters)
    dist_unit = params.get("distance_units")
    return PDMSParams(
        b=float(params.get("<b>").to(dist_unit).magnitude),
        b2=float(params.get("<b^2>").to(dist_unit**2).magnitude),
        r02=float(params.get("R02").to(dist_unit**2).magnitude),
        bead_density=float(params.get_bead_density()),
        entanglement_density=float(params.get_entanglement_density()),
        sampling_cutoff=float(params.get_sampling_cutoff()),
        gamma_to_mpa=float(params.get_gamma_conversion_factor().to("MPa").magnitude),
    )


def _configure_max_distance(
    gen: pt.MCUniverseGenerator,
    cfg: CutoffConfig,
    b2: float,
    n_bonds: int,
) -> float:
    """Apply a ``CutoffConfig`` to ``gen`` and return the physical cutoff distance."""
    if cfg.mode is MaxDistanceMode.ZSCORE:
        gen.use_zscore_max_distance(cfg.value, b2)
        return cfg.value * float(np.sqrt(n_bonds * b2))
    if cfg.mode is MaxDistanceMode.LINEAR:
        gen.use_linear_max_distance(cfg.value)
        return cfg.value * float(n_bonds)
    gen.disable_max_distance()
    return -1.0


def _compute_loop_fractions(net: pt.MEHPNetwork) -> Tuple[float, float]:
    """Return ``(primary_loop_fraction, secondary_loop_strand_fraction)``."""
    a = np.asarray(net.spring_index_a, dtype=np.int64)
    b_idx = np.asarray(net.spring_index_b, dtype=np.int64)
    n_tot = len(a)
    if n_tot == 0:
        return 0.0, 0.0
    n_prim = int(np.sum(a == b_idx))
    mask = a != b_idx
    if not np.any(mask):
        return n_prim / n_tot, 0.0
    pairs = np.stack(
        [np.minimum(a[mask], b_idx[mask]), np.maximum(a[mask], b_idx[mask])],
        axis=1,
    )
    _, counts = np.unique(pairs, axis=0, return_counts=True)
    n_sec_strands = int(np.sum(counts[counts >= 2]))
    return n_prim / n_tot, n_sec_strands / n_tot


def _run_single_mc_cutoff_seed(seed: int) -> Dict[str, Any]:
    """Evaluate all ``CUTOFF_CONFIGS`` for one MC random seed."""
    p = extract_pdms_parameters()
    n_str = MC_N_STRANDS
    n_xl = n_str // 2
    n_beads = MC_N_BEADS
    n_bonds = n_beads + 1
    r02_strand = n_bonds * p.b2
    r_rms = float(np.sqrt(r02_strand))

    vol = (n_str * n_beads + n_xl) / p.bead_density
    box_l = vol ** (1.0 / 3.0)
    n_ent = int(round(p.entanglement_density * vol))

    n_cfg = len(CUTOFF_CONFIGS)
    r_edges = np.linspace(0.0, MC_R_HIST_MAX, MC_R_HIST_BINS + 1)

    wall_time = np.zeros(n_cfg, dtype=float)
    conversion = np.zeros(n_cfg, dtype=float)
    r2_norm = np.zeros(n_cfg, dtype=float)
    r_max_norm = np.zeros(n_cfg, dtype=float)
    fallback_count = np.zeros(n_cfg, dtype=float)
    prim_frac = np.zeros(n_cfg, dtype=float)
    sec_frac = np.zeros(n_cfg, dtype=float)
    g_ph = np.zeros(n_cfg, dtype=float)
    g_ent = np.zeros(n_cfg, dtype=float)
    r_hist = np.zeros((n_cfg, MC_R_HIST_BINS), dtype=float)

    for i, cfg in enumerate(CUTOFF_CONFIGS):
        gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
        gen.set_seed(seed)
        gen.config_nr_of_mc_steps(0)
        gen.set_bead_distance(p.b)
        gen.set_mean_squared_bead_distance(p.b2)
        gen.add_crosslinkers(n_xl, 4)
        gen.add_strands(n_str, [n_beads] * n_str)
        r_cutoff = _configure_max_distance(gen, cfg, p.b2, n_bonds)

        t0 = time.perf_counter()
        gen.link_strands_to_conversion(MC_TARGET_CONVERSION)
        wall_time[i] = time.perf_counter() - t0
        conversion[i] = float(gen.get_current_crosslinker_conversion())

        fb_ph = gen.get_mehp_simulator(
            entanglement_mode=pt.EntanglementMode.NONE,
            crosslinker_type=2,
        )
        r_init = np.asarray(fb_ph.get_spring_lengths(), dtype=float)
        r_nz = r_init[r_init > 0.0]
        if len(r_nz) > 0:
            r2_norm[i] = float(np.mean(r_nz**2) / r02_strand)
            r_max_norm[i] = float(np.max(r_nz) / r_rms)
            hist_vals, _ = np.histogram(r_nz / r_rms, bins=r_edges, density=True)
            r_hist[i] = hist_vals
        if r_cutoff > 0.0:
            fallback_count[i] = float(np.sum(r_init > r_cutoff * (1.0 + 1e-9)))

        prim_f, sec_f = _compute_loop_fractions(fb_ph.network)
        prim_frac[i] = prim_f
        sec_frac[i] = sec_f

        fb_ph.run()
        g_ph[i] = p.gamma_to_mpa * float(np.sum(fb_ph.get_gamma_factors(p.r02))) / vol

        es = gen.sample_entanglements(
            n_ent,
            p.sampling_cutoff,
            seed=f"mc_cutoff_{cfg.label}_{seed}",
        )
        fb_ent = gen.get_mehp_simulator(
            es,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            crosslinker_type=2,
        )
        fb_ent.run()
        g_ent[i] = p.gamma_to_mpa * float(np.sum(fb_ent.get_gamma_factors(p.r02))) / vol

    return {
        "seed": seed,
        "wall_time": wall_time,
        "conversion": conversion,
        "r2_norm": r2_norm,
        "r_max_norm": r_max_norm,
        "fallback_count": fallback_count,
        "prim_frac": prim_frac,
        "sec_frac": sec_frac,
        "g_ph": g_ph,
        "g_ent": g_ent,
        "r_hist": r_hist,
    }


def build_pbc_aware_network(
    box_l: float,
    n_beads: int,
    n_strands: int,
    b2: float,
    p_target: float,
    seed: int,
    is_2d: bool = False,
) -> pt.Universe:
    """Build a melt-like periodic Gaussian network allowing multi-image wrapping.

    Strands are linked across all periodic images :math:`\\mathbf{n} \\in \\mathbb{Z}^d`
    within :math:`3.5\\,R_{\\text{rms}}` weighted by the exact Gaussian probability
    :math:`\\exp(-d\\,R^2 / (2(N+1)b^2))`. Interior beads are placed along the
    unwrapped strand path so every bead-to-bead bond satisfies
    :math:`|\\mathbf{b}_k| \\ll L/2`, enabling
    ``MEHPSimulator(..., assume_box_large_enough=False)`` to reconstruct exact
    baked periodic offsets :math:`\\mathbf{o}_s`.
    """
    rng = np.random.default_rng(seed)
    n_xl = n_strands // 2
    n_bonds = n_beads + 1
    r02_strand = n_bonds * b2
    r_rms = float(np.sqrt(r02_strand))

    xl_pos = rng.uniform(0.0, box_l, size=(n_xl, 3))
    if is_2d:
        xl_pos[:, 2] = 0.0
    rem_f = np.full(n_xl, 4, dtype=np.int64)

    n_img = max(1, int(np.ceil(3.5 * r_rms / box_l)))
    z_range = [0] if is_2d else list(range(-n_img, n_img + 1))
    shifts = (
        np.array(
            [
                [ix, iy, iz]
                for ix in range(-n_img, n_img + 1)
                for iy in range(-n_img, n_img + 1)
                for iz in z_range
            ],
            dtype=float,
        )
        * box_l
    )
    dim_factor = 1.0 if is_2d else 1.5  # d / 2

    n_links_target = int(round(p_target * 4 * n_xl))
    strand_from = np.full(n_strands, -1, dtype=np.int64)
    strand_to = np.full(n_strands, -1, dtype=np.int64)
    strand_shift = np.zeros((n_strands, 3), dtype=float)

    ends = np.repeat(np.arange(n_strands, dtype=np.int64), 2)
    rng.shuffle(ends)
    links_done = 0

    for s in ends:
        if links_done >= n_links_target:
            break
        if strand_from[s] < 0:
            avail = np.where(rem_f > 0)[0]
            if len(avail) == 0:
                break
            probs = rem_f[avail].astype(float)
            probs /= probs.sum()
            chosen = int(rng.choice(avail, p=probs))
            strand_from[s] = chosen
            rem_f[chosen] -= 1
            links_done += 1
        else:
            src = int(strand_from[s])
            avail = np.where(rem_f > 0)[0]
            if len(avail) == 0:
                break
            base_diff = xl_pos[avail] - xl_pos[src]
            all_diff = base_diff[:, None, :] + shifts[None, :, :]
            r2 = np.sum(all_diff**2, axis=2)
            weights = rem_f[avail, None].astype(float) * np.exp(
                -dim_factor * r2 / r02_strand
            )
            w_flat = weights.ravel()
            w_sum = float(w_flat.sum())
            if w_sum <= 0.0:
                continue
            idx = int(rng.choice(len(w_flat), p=w_flat / w_sum))
            a_i, m_i = divmod(idx, len(shifts))
            dst = int(avail[a_i])
            strand_to[s] = dst
            strand_shift[s] = shifts[m_i]
            rem_f[dst] -= 1
            links_done += 1

    ids: List[int] = []
    types: List[int] = []
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    b_from: List[int] = []
    b_to: List[int] = []

    for i in range(n_xl):
        ids.append(i + 1)
        types.append(2)
        xs.append(float(xl_pos[i, 0]))
        ys.append(float(xl_pos[i, 1]))
        zs.append(float(xl_pos[i, 2]))

    next_id = n_xl + 1
    for s in range(n_strands):
        u = int(strand_from[s])
        v = int(strand_to[s])
        if u < 0 or v < 0:
            continue
        p0 = xl_pos[u]
        p1 = xl_pos[v] + strand_shift[s]
        prev_id = u + 1
        for k in range(1, n_beads + 1):
            frac = k / float(n_bonds)
            pos = p0 + frac * (p1 - p0)
            pos[0] %= box_l
            pos[1] %= box_l
            pos[2] = 0.0 if is_2d else (pos[2] % box_l)
            ids.append(next_id)
            types.append(1)
            xs.append(float(pos[0]))
            ys.append(float(pos[1]))
            zs.append(float(pos[2]))
            b_from.append(prev_id)
            b_to.append(next_id)
            prev_id = next_id
            next_id += 1
        b_from.append(prev_id)
        b_to.append(v + 1)

    univ = pt.Universe(box_l, box_l, box_l)
    zeros = [0] * len(ids)
    univ.add_atoms(ids, types, xs, ys, zs, zeros, zeros, zeros)
    univ.add_bonds(b_from, b_to)
    return univ


def _run_pbc_alpha_sweep_for_seed(seed: int) -> Dict[str, Any]:
    """Evaluate all ``PBC_ALPHAS`` under both PBC modes for one random seed."""
    p = extract_pdms_parameters()
    n_beads = PBC_N_BEADS
    r_rms = float(np.sqrt((n_beads + 1) * p.b2))
    n_alpha = len(PBC_ALPHAS)

    g_baked = np.zeros(n_alpha, dtype=float)
    g_mic = np.zeros(n_alpha, dtype=float)
    g_ratio = np.zeros(n_alpha, dtype=float)
    frac_len_gt_half = np.zeros(n_alpha, dtype=float)
    frac_comp_gt_half = np.zeros(n_alpha, dtype=float)
    frac_flipped = np.zeros(n_alpha, dtype=float)

    rl_edges = np.linspace(0.0, PBC_RL_HIST_MAX, PBC_RL_HIST_BINS + 1)
    hist_baked = np.zeros(PBC_RL_HIST_BINS, dtype=float)
    hist_mic = np.zeros(PBC_RL_HIST_BINS, dtype=float)

    for idx, alpha in enumerate(PBC_ALPHAS):
        box_l = alpha * r_rms
        vol = box_l**3
        n_str_raw = int(round(p.bead_density * vol / (n_beads + 0.5)))
        n_str = max(12, n_str_raw + (n_str_raw % 2))

        univ = build_pbc_aware_network(
            box_l=box_l,
            n_beads=n_beads,
            n_strands=n_str,
            b2=p.b2,
            p_target=PBC_TARGET_CONVERSION,
            seed=seed + idx * 1009,
            is_2d=False,
        )

        sim_baked = pt.MEHPSimulator(
            univ,
            crosslinker_type=2,
            assume_box_large_enough=False,
        )
        sim_mic = pt.MEHPSimulator(
            univ,
            crosslinker_type=2,
            assume_box_large_enough=True,
        )
        sim_baked.run()
        sim_mic.run()

        gb = p.gamma_to_mpa * float(np.sum(sim_baked.get_gamma_factors(p.r02))) / vol
        gm = p.gamma_to_mpa * float(np.sum(sim_mic.get_gamma_factors(p.r02))) / vol
        g_baked[idx] = gb
        g_mic[idx] = gm
        g_ratio[idx] = gm / gb if gb > 1e-12 else np.nan

        r_b = np.asarray(sim_baked.get_spring_lengths(), dtype=float)
        r_m = np.asarray(sim_mic.get_spring_lengths(), dtype=float)
        v_b = np.asarray(sim_baked.get_spring_vectors(), dtype=float).reshape(-1, 3)
        v_m = np.asarray(sim_mic.get_spring_vectors(), dtype=float).reshape(-1, 3)

        if len(r_b) > 0:
            frac_len_gt_half[idx] = float(np.mean(r_b > 0.5 * box_l))
            frac_comp_gt_half[idx] = float(
                np.mean(np.max(np.abs(v_b), axis=1) > 0.5 * box_l)
            )
            frac_flipped[idx] = float(
                np.mean(np.linalg.norm(v_b - v_m, axis=1) > 0.25 * box_l)
            )

        if idx == PBC_HIST_ALPHA_INDEX and len(r_b) > 0:
            hb, _ = np.histogram(r_b / box_l, bins=rl_edges, density=True)
            hm, _ = np.histogram(r_m / box_l, bins=rl_edges, density=True)
            hist_baked = hb
            hist_mic = hm

    return {
        "seed": seed,
        "g_baked": g_baked,
        "g_mic": g_mic,
        "g_ratio": g_ratio,
        "frac_len_gt_half": frac_len_gt_half,
        "frac_comp_gt_half": frac_comp_gt_half,
        "frac_flipped": frac_flipped,
        "hist_baked": hist_baked,
        "hist_mic": hist_mic,
    }


def _build_2d_pbc_flip_demo() -> Dict[str, Any]:
    """Build and relax a small 2D periodic network under both PBC modes."""
    p = extract_pdms_parameters()
    n_beads = 30
    r_rms = float(np.sqrt((n_beads + 1) * p.b2))
    box_l = 1.05 * r_rms
    n_strands = 16

    univ = build_pbc_aware_network(
        box_l=box_l,
        n_beads=n_beads,
        n_strands=n_strands,
        b2=p.b2,
        p_target=1.0,
        seed=777,
        is_2d=True,
    )

    sim_baked = pt.MEHPSimulator(
        univ,
        crosslinker_type=2,
        is_2d=True,
        assume_box_large_enough=False,
    )
    sim_mic = pt.MEHPSimulator(
        univ,
        crosslinker_type=2,
        is_2d=True,
        assume_box_large_enough=True,
    )
    sim_baked.run()
    sim_mic.run()

    coords_baked = np.asarray(sim_baked.get_coordinates(), dtype=float).reshape(-1, 3)
    coords_mic = np.asarray(sim_mic.get_coordinates(), dtype=float).reshape(-1, 3)
    coords_baked_2d = np.mod(coords_baked[:, :2], box_l)
    coords_mic_2d = np.mod(coords_mic[:, :2], box_l)

    vec_baked_2d = np.asarray(sim_baked.get_spring_vectors(), dtype=float).reshape(
        -1, 3
    )[:, :2]
    vec_mic_2d = np.asarray(sim_mic.get_spring_vectors(), dtype=float).reshape(-1, 3)[
        :, :2
    ]

    net_b = sim_baked.network
    idx_a = np.asarray(net_b.spring_index_a, dtype=np.int64)
    idx_b = np.asarray(net_b.spring_index_b, dtype=np.int64)

    gamma_baked = float(np.sum(sim_baked.get_gamma_factors(p.r02)))
    gamma_mic = float(np.sum(sim_mic.get_gamma_factors(p.r02)))

    return {
        "demo2d_box_l": np.array(box_l, dtype=float),
        "demo2d_r_rms": np.array(r_rms, dtype=float),
        "demo2d_coords_baked": coords_baked_2d,
        "demo2d_coords_mic": coords_mic_2d,
        "demo2d_vec_baked": vec_baked_2d,
        "demo2d_vec_mic": vec_mic_2d,
        "demo2d_spring_a": idx_a,
        "demo2d_spring_b": idx_b,
        "demo2d_gamma_baked": np.array(gamma_baked, dtype=float),
        "demo2d_gamma_mic": np.array(gamma_mic, dtype=float),
    }


def _pack_stats(
    out: Dict[str, Any],
    prefix: str,
    samples: np.ndarray,
    axis: int = 0,
) -> None:
    """Compute ensemble statistics and store mean/std/ci into ``out``."""
    stats = compute_ensemble_statistics(samples, axis=axis)
    out[f"{prefix}_mean"] = stats["mean"]
    out[f"{prefix}_std"] = stats["std"]
    out[f"{prefix}_ci"] = stats["ci"]


def run_mc_pbc_artifacts_study(max_workers: Optional[int] = None) -> Dict[str, Any]:
    """Run Study 5 (MC ``MaxDistanceProvider`` sweep + ``assume_box_large_enough`` sweep).

    Parameters
    ----------
    max_workers : int, optional
        Maximum worker processes for parallel seed evaluation. Defaults to
        :func:`_example_cache.resolve_max_workers`.

    Returns
    -------
    dict of str to np.ndarray
        Aggregated ensemble statistics, histograms, and 2D demo geometry ready
        for ``np.savez_compressed``.
    """
    workers = resolve_max_workers(max_workers)
    configure_process_threading(workers=workers, verbose=True)

    p = extract_pdms_parameters()
    r_rms_mc = float(np.sqrt((MC_N_BEADS + 1) * p.b2))
    vol_mc = (MC_N_STRANDS * MC_N_BEADS + MC_N_STRANDS // 2) / p.bead_density
    box_l_mc = vol_mc ** (1.0 / 3.0)

    # Part 1: MC cutoff sweep across seeds
    mc_results: List[Dict[str, Any]] = []
    if workers == 1:
        for seed in MC_SEEDS:
            mc_results.append(_run_single_mc_cutoff_seed(seed))
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = {
                executor.submit(_run_single_mc_cutoff_seed, seed): seed
                for seed in MC_SEEDS
            }
            for fut in as_completed(futures):
                mc_results.append(fut.result())
    mc_results.sort(key=lambda r: int(r["seed"]))

    # Part 2: PBC assume_box_large_enough sweep across seeds
    pbc_results: List[Dict[str, Any]] = []
    if workers == 1:
        for seed in PBC_SEEDS:
            pbc_results.append(_run_pbc_alpha_sweep_for_seed(seed))
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = {
                executor.submit(_run_pbc_alpha_sweep_for_seed, seed): seed
                for seed in PBC_SEEDS
            }
            for fut in as_completed(futures):
                pbc_results.append(fut.result())
    pbc_results.sort(key=lambda r: int(r["seed"]))

    # Part 3: 2D visual demo
    demo2d = _build_2d_pbc_flip_demo()

    # Assemble output dictionary
    out: Dict[str, Any] = {}
    out["pdms_b"] = np.array(p.b, dtype=float)
    out["pdms_b2"] = np.array(p.b2, dtype=float)
    out["pdms_r02"] = np.array(p.r02, dtype=float)
    out["pdms_rho"] = np.array(p.bead_density, dtype=float)
    out["pdms_rho_e"] = np.array(p.entanglement_density, dtype=float)
    out["pdms_gamma_to_mpa"] = np.array(p.gamma_to_mpa, dtype=float)

    # Part 1 metadata & arrays
    out["mc_n_beads"] = np.array(MC_N_BEADS, dtype=np.int64)
    out["mc_n_strands"] = np.array(MC_N_STRANDS, dtype=np.int64)
    out["mc_target_conversion"] = np.array(MC_TARGET_CONVERSION, dtype=float)
    out["mc_box_l"] = np.array(box_l_mc, dtype=float)
    out["mc_r_rms"] = np.array(r_rms_mc, dtype=float)
    out["mc_n_seeds"] = np.array(len(MC_SEEDS), dtype=np.int64)
    out["mc_labels"] = np.array([c.label for c in CUTOFF_CONFIGS])
    out["mc_modes"] = np.array([c.mode.value for c in CUTOFF_CONFIGS])
    out["mc_values"] = np.array([c.value for c in CUTOFF_CONFIGS], dtype=float)
    out["mc_r_edges"] = np.linspace(0.0, MC_R_HIST_MAX, MC_R_HIST_BINS + 1)

    wall_times = np.stack([r["wall_time"] for r in mc_results], axis=0)
    ref_times = wall_times[:, -1:]  # last config is DISABLED
    speedups = ref_times / np.maximum(wall_times, 1e-9)

    g_ph_all = np.stack([r["g_ph"] for r in mc_results], axis=0)
    g_ent_all = np.stack([r["g_ent"] for r in mc_results], axis=0)
    g_ph_rel_err = np.abs(g_ph_all - g_ph_all[:, -1:]) / np.maximum(
        g_ph_all[:, -1:], 1e-12
    )
    g_ent_rel_err = np.abs(g_ent_all - g_ent_all[:, -1:]) / np.maximum(
        g_ent_all[:, -1:], 1e-12
    )

    _pack_stats(out, "mc_wall_time", wall_times)
    _pack_stats(out, "mc_speedup", speedups)
    _pack_stats(out, "mc_conversion", np.stack([r["conversion"] for r in mc_results]))
    _pack_stats(out, "mc_r2_norm", np.stack([r["r2_norm"] for r in mc_results]))
    _pack_stats(out, "mc_r_max_norm", np.stack([r["r_max_norm"] for r in mc_results]))
    _pack_stats(
        out, "mc_fallback_count", np.stack([r["fallback_count"] for r in mc_results])
    )
    _pack_stats(out, "mc_prim_frac", np.stack([r["prim_frac"] for r in mc_results]))
    _pack_stats(out, "mc_sec_frac", np.stack([r["sec_frac"] for r in mc_results]))
    _pack_stats(out, "mc_g_ph", g_ph_all)
    _pack_stats(out, "mc_g_ent", g_ent_all)
    _pack_stats(out, "mc_g_ph_rel_err", g_ph_rel_err)
    _pack_stats(out, "mc_g_ent_rel_err", g_ent_rel_err)
    out["mc_r_hist_mean"] = np.mean(
        np.stack([r["r_hist"] for r in mc_results], axis=0), axis=0
    )

    # Part 2 metadata & arrays
    out["pbc_n_beads"] = np.array(PBC_N_BEADS, dtype=np.int64)
    out["pbc_r_rms"] = np.array(np.sqrt((PBC_N_BEADS + 1) * p.b2), dtype=float)
    out["pbc_alphas"] = np.array(PBC_ALPHAS, dtype=float)
    out["pbc_n_seeds"] = np.array(len(PBC_SEEDS), dtype=np.int64)
    out["pbc_hist_alpha"] = np.array(PBC_ALPHAS[PBC_HIST_ALPHA_INDEX], dtype=float)
    out["pbc_rl_edges"] = np.linspace(0.0, PBC_RL_HIST_MAX, PBC_RL_HIST_BINS + 1)

    _pack_stats(out, "pbc_g_baked", np.stack([r["g_baked"] for r in pbc_results]))
    _pack_stats(out, "pbc_g_mic", np.stack([r["g_mic"] for r in pbc_results]))
    _pack_stats(out, "pbc_g_ratio", np.stack([r["g_ratio"] for r in pbc_results]))
    _pack_stats(
        out,
        "pbc_frac_len_gt_half",
        np.stack([r["frac_len_gt_half"] for r in pbc_results]),
    )
    _pack_stats(
        out,
        "pbc_frac_comp_gt_half",
        np.stack([r["frac_comp_gt_half"] for r in pbc_results]),
    )
    _pack_stats(
        out, "pbc_frac_flipped", np.stack([r["frac_flipped"] for r in pbc_results])
    )
    out["pbc_hist_baked_mean"] = np.mean(
        np.stack([r["hist_baked"] for r in pbc_results], axis=0), axis=0
    )
    out["pbc_hist_mic_mean"] = np.mean(
        np.stack([r["hist_mic"] for r in pbc_results], axis=0), axis=0
    )

    out.update(demo2d)
    return out


if __name__ == "__main__":
    cache_path = get_example_cache_path(
        "mc_cutoff_and_pbc_artifacts_cache.npz",
        example_folder="network_generator",
    )
    load_or_compute_cache(
        cache_path,
        run_mc_pbc_artifacts_study,
        recompute=True,
    )
