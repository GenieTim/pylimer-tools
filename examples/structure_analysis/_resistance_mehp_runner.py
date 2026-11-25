"""
Multi-seed simulation runner for Study 7:
Topological Effective Resistance vs. Spatial Stress Concentration in Defect-Rich Networks.

This module is prefixed with an underscore so Sphinx-Gallery ignores it when
building documentation.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from enum import Enum, IntEnum
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
from scipy import stats

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
from pylimer_tools.io.bead_spring_parameter_provider import (  # noqa: E402
    ParameterType,
    get_parameters_for_polymer,
)
import pylimer_tools_cpp as pt  # noqa: E402


class NetworkRegime(Enum):
    """Defect-rich and heterogeneous polymer network regimes."""

    NEAR_PERCOLATION = "near_percolation"
    BIMODAL_DEFECT = "bimodal_defect"


class StrandTopologyClass(IntEnum):
    """Topological and mechanical classification of network strands."""

    PRIMARY_LOOP = 0
    BRIDGE_CUT_EDGE = 1
    LOCAL_CYCLE = 2
    PBC_ACTIVE_BACKBONE = 3


class TopologicalDescriptor(Enum):
    """Per-strand topological predictors derived from the weighted Laplacian."""

    LEVERAGE = "leverage"
    STIFFNESS = "stiffness"
    MEAN_CENTRALITY = "mean_centrality"


class MechanicalObservable(Enum):
    """Per-strand mechanical load-bearing observables from MEHP force balance."""

    GAMMA_EQ = "gamma_eq"
    DELTA_GAMMA = "delta_gamma"
    FORCE_EQ = "force_eq"
    DELTA_FORCE = "delta_force"


@dataclass(frozen=True)
class _RegimeConfig:
    """Synthesis parameters for a single network regime."""

    regime: NetworkRegime
    label: str
    conversion: float
    strand_lengths: Tuple[int, ...]
    n_monofunctional: int
    monofunctional_length: int
    functionality: int


@dataclass(frozen=True)
class _PdmsConstants:
    """Numeric PDMS Gaussian bead-spring constants."""

    b: float
    b2_ms: float
    r02: float
    bead_density: float
    ent_density: float
    cutoff: float
    gamma_conv: float


SEEDS: Tuple[int, ...] = (101, 102, 103, 104, 105, 106, 107, 108)
REPRESENTATIVE_SEED: int = SEEDS[0]
DEFORMATION_LAMBDA: float = 1.5
ACTIVE_GAMMA_TOL: float = 1e-6
BRIDGE_LEVERAGE_TOL: float = 1e-6
LEVERAGE_BIN_EDGES: Tuple[float, ...] = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
ACTIVE_DEGREE_BINS: Tuple[int, ...] = (0, 2, 3, 4)

_REGIME_CONFIGS: Tuple[_RegimeConfig, ...] = (
    _RegimeConfig(
        regime=NetworkRegime.NEAR_PERCOLATION,
        label="Near-Percolation (p=0.68, N=50)",
        conversion=0.68,
        strand_lengths=(50,) * 1000,
        n_monofunctional=0,
        monofunctional_length=50,
        functionality=4,
    ),
    _RegimeConfig(
        regime=NetworkRegime.BIMODAL_DEFECT,
        label="Bimodal + Defects (p=0.92, N=15/90)",
        conversion=0.92,
        strand_lengths=(15,) * 500 + (90,) * 500,
        n_monofunctional=100,
        monofunctional_length=50,
        functionality=4,
    ),
)

_DESCRIPTORS: Tuple[TopologicalDescriptor, ...] = (
    TopologicalDescriptor.LEVERAGE,
    TopologicalDescriptor.STIFFNESS,
    TopologicalDescriptor.MEAN_CENTRALITY,
)

_OBSERVABLES: Tuple[MechanicalObservable, ...] = (
    MechanicalObservable.GAMMA_EQ,
    MechanicalObservable.DELTA_GAMMA,
    MechanicalObservable.FORCE_EQ,
    MechanicalObservable.DELTA_FORCE,
)


def _load_pdms_constants() -> _PdmsConstants:
    """Load PDMS Gaussian bead-spring constants from ParameterProvider."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    dist_unit = params.get("distance_units")
    return _PdmsConstants(
        b=float(params.get("<b>").to(dist_unit).magnitude),
        b2_ms=float(params.get("<b^2>").to(dist_unit**2).magnitude),
        r02=float(params.get("R02").to(dist_unit**2).magnitude),
        bead_density=float(params.get_bead_density()),
        ent_density=float(params.get_entanglement_density()),
        cutoff=float(params.get_sampling_cutoff()),
        gamma_conv=float(params.get_gamma_conversion_factor().to("MPa").magnitude),
    )


def _build_strand_laplacian(
    spring_a: np.ndarray,
    spring_b: np.ndarray,
    weights: np.ndarray,
    num_nodes: int,
) -> sp.csc_matrix:
    """Assemble the exact strand-level weighted graph Laplacian L = D - W.

    Primary loops (spring_a == spring_b) contribute zero to L and are skipped.
    """
    non_self = spring_a != spring_b
    a_ns = spring_a[non_self]
    b_ns = spring_b[non_self]
    w_ns = weights[non_self]
    rows = np.concatenate([a_ns, b_ns, a_ns, b_ns])
    cols = np.concatenate([a_ns, b_ns, b_ns, a_ns])
    vals = np.concatenate([w_ns, w_ns, -w_ns, -w_ns])
    return sp.csc_matrix((vals, (rows, cols)), shape=(num_nodes, num_nodes))


def _safe_correlation(
    x: np.ndarray,
    y: np.ndarray,
) -> Tuple[float, float]:
    """Compute (Spearman rho, Pearson r), returning (nan, nan) if degenerate."""
    if x.size < 3 or np.std(x) <= 0.0 or np.std(y) <= 0.0:
        return (float("nan"), float("nan"))
    rho = float(stats.spearmanr(x, y).statistic)
    r_val = float(stats.pearsonr(x, y).statistic)
    return (rho, r_val)


def _generate_universe(
    cfg: _RegimeConfig,
    seed: int,
    constants: _PdmsConstants,
) -> Tuple[pt.Universe, float]:
    """Synthesize a single MC PDMS network realization and return (universe, volume)."""
    n_strands = len(cfg.strand_lengths)
    n_xlinks = (2 * n_strands + cfg.n_monofunctional) // cfg.functionality
    total_beads = (
        sum(cfg.strand_lengths)
        + cfg.n_monofunctional * cfg.monofunctional_length
        + n_xlinks
    )
    volume = float(total_beads) / constants.bead_density
    box_l = volume ** (1.0 / 3.0)

    gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(seed)
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(constants.b)
    gen.set_mean_squared_bead_distance(constants.b2_ms)
    gen.add_crosslinkers(n_xlinks, cfg.functionality)
    gen.add_strands(n_strands, list(cfg.strand_lengths), strand_atom_type=1)
    if cfg.n_monofunctional > 0:
        gen.add_monofunctional_strands(
            cfg.n_monofunctional,
            [cfg.monofunctional_length] * cfg.n_monofunctional,
            monofunctional_strand_atom_type=4,
        )
    gen.link_strands_to_conversion(cfg.conversion)
    return gen.get_universe(), volume


def _classify_strands(
    spring_a: np.ndarray,
    spring_b: np.ndarray,
    leverage: np.ndarray,
    gamma_eq: np.ndarray,
) -> np.ndarray:
    """Partition strands into primary loops, bridges, local cycles, and PBC active."""
    classes = np.full(
        spring_a.shape[0], int(StrandTopologyClass.LOCAL_CYCLE), dtype=np.int32
    )
    is_primary = spring_a == spring_b
    is_bridge = (~is_primary) & (np.abs(leverage - 1.0) < BRIDGE_LEVERAGE_TOL)
    is_active = gamma_eq > ACTIVE_GAMMA_TOL

    classes[is_primary] = int(StrandTopologyClass.PRIMARY_LOOP)
    classes[is_bridge] = int(StrandTopologyClass.BRIDGE_CUT_EDGE)
    classes[(~is_primary) & (~is_bridge) & is_active] = int(
        StrandTopologyClass.PBC_ACTIVE_BACKBONE
    )
    # Any PBC-wrapping primary loop (if it occurs) is also active backbone
    classes[is_primary & is_active] = int(StrandTopologyClass.PBC_ACTIVE_BACKBONE)
    return classes


def _run_single_task(
    task: Tuple[str, int, _PdmsConstants],
) -> Dict[str, Any]:
    """Run topological resistance and MEHP force-balance analysis for one (regime, seed)."""
    regime_value, seed, constants = task
    regime = NetworkRegime(regime_value)
    cfg = next(c for c in _REGIME_CONFIGS if c.regime == regime)

    universe, volume = _generate_universe(cfg, seed, constants)
    b0 = float(np.sqrt(constants.r02))

    # 1. Phantom MEHP network (EntanglementMode.NONE) with 1-to-1 spring indexing
    net_ph = pt.MEHPNetwork.from_universe(
        universe,
        entanglement_mode=pt.EntanglementMode.NONE,
        crosslinker_type=2,
        assume_box_large_enough=False,
    )
    sim_ph = pt.MEHPSimulator(net_ph)
    solver_cfg_nosimp = pt.MEHPLinearSolverConfig()
    solver_cfg_nosimp.simplification_mode = (
        pt.StructureSimplificationMode.NO_SIMPLIFICATION
    )
    sim_ph.run(solver_cfg_nosimp)

    gamma_eq = np.asarray(sim_ph.get_gamma_factors(constants.r02), dtype=np.float64)
    lengths_eq = np.asarray(sim_ph.get_spring_lengths(), dtype=np.float64)
    n_s_all = np.asarray(net_ph.spring_contour_length, dtype=np.float64)
    # Dimensionless tension ||F_s|| / (3 k_B T / b_0) = ||Delta r_s|| / (N_s * b_0)
    force_eq = lengths_eq / (n_s_all * b0)
    w_act_ph = float(sim_ph.get_active_weight_fraction(ACTIVE_GAMMA_TOL))
    g_ph_mpa = float(constants.gamma_conv * np.sum(gamma_eq) / volume)

    # Uniaxial isochoric deformation (lambda = DEFORMATION_LAMBDA along x)
    box0 = net_ph.box
    inv_sqrt_lam = 1.0 / np.sqrt(DEFORMATION_LAMBDA)
    box_def = pt.Box(
        box0.get_lx() * DEFORMATION_LAMBDA,
        box0.get_ly() * inv_sqrt_lam,
        box0.get_lz() * inv_sqrt_lam,
    )
    sim_ph.deform_to(box_def)
    sim_ph.run(solver_cfg_nosimp)

    gamma_def = np.asarray(sim_ph.get_gamma_factors(constants.r02), dtype=np.float64)
    lengths_def = np.asarray(sim_ph.get_spring_lengths(), dtype=np.float64)
    force_def = lengths_def / (n_s_all * b0)
    delta_gamma = gamma_def - gamma_eq
    delta_force = force_def - force_eq

    # Quenched periodic boundary offsets per strand
    offsets_all = np.asarray(net_ph.spring_box_offset, dtype=np.float64).reshape(-1, 3)
    offset_norm_all = np.linalg.norm(offsets_all, axis=1)

    # 2. Exact strand-level weighted graph Laplacian on all nodes & giant component
    spring_a = np.asarray(net_ph.spring_index_a, dtype=np.int64)
    spring_b = np.asarray(net_ph.spring_index_b, dtype=np.int64)
    weights_all = 1.0 / n_s_all

    lap_full = _build_strand_laplacian(
        spring_a, spring_b, weights_all, int(net_ph.nr_of_links)
    )
    rda_full = pt.ResistanceDistanceAnalysis.from_matrix(lap_full)
    components = rda_full.get_connected_components()
    gcc_nodes = max(components, key=len)
    gcc_map = {int(old): idx for idx, old in enumerate(gcc_nodes)}

    # Classify all strands across all components
    pairs_all = [(int(u), int(v)) for u, v in zip(spring_a, spring_b)]
    r_ab_all = np.asarray(
        rda_full.compute_effective_resistance_batched(pairs_all), dtype=np.float64
    )
    leverage_all = weights_all * r_ab_all
    classes_all = _classify_strands(spring_a, spring_b, leverage_all, gamma_eq)

    # Restrict to GCC so R_ij and node centrality are finite and well-defined
    in_gcc = np.array(
        [
            (int(u) in gcc_map) and (int(v) in gcc_map)
            for u, v in zip(spring_a, spring_b)
        ],
        dtype=bool,
    )
    a_gcc = np.array([gcc_map[int(u)] for u in spring_a[in_gcc]], dtype=np.int64)
    b_gcc = np.array([gcc_map[int(v)] for v in spring_b[in_gcc]], dtype=np.int64)
    w_gcc = weights_all[in_gcc]
    n_s_gcc = n_s_all[in_gcc]
    gamma_eq_gcc = gamma_eq[in_gcc]
    gamma_def_gcc = gamma_def[in_gcc]
    dgamma_gcc = delta_gamma[in_gcc]
    force_eq_gcc = force_eq[in_gcc]
    force_def_gcc = force_def[in_gcc]
    dforce_gcc = delta_force[in_gcc]
    offset_norm_gcc = offset_norm_all[in_gcc]

    n_gcc = len(gcc_nodes)
    lap_gcc = _build_strand_laplacian(a_gcc, b_gcc, w_gcc, n_gcc)
    rda_gcc = pt.ResistanceDistanceAnalysis.from_matrix(lap_gcc)
    pairs_gcc = [(int(u), int(v)) for u, v in zip(a_gcc, b_gcc)]
    r_ab_gcc = np.asarray(
        rda_gcc.compute_effective_resistance_batched(pairs_gcc), dtype=np.float64
    )
    leverage_gcc = w_gcc * r_ab_gcc
    node_cent = np.asarray(rda_gcc.compute_resistance_centrality(), dtype=np.float64)
    mean_cent_gcc = 0.5 * (node_cent[a_gcc] + node_cent[b_gcc])
    classes_gcc = _classify_strands(a_gcc, b_gcc, leverage_gcc, gamma_eq_gcc)

    # Node topological degree and active connectivity in GCC
    node_topo_deg = np.zeros(n_gcc, dtype=np.int32)
    node_active_deg = np.zeros(n_gcc, dtype=np.int32)
    for u, v, g_val in zip(a_gcc, b_gcc, gamma_eq_gcc):
        if u == v:
            node_topo_deg[u] += 2
        else:
            node_topo_deg[u] += 1
            node_topo_deg[v] += 1
            if g_val > ACTIVE_GAMMA_TOL:
                node_active_deg[u] += 1
                node_active_deg[v] += 1

    # Class counts & fractions (GCC and full network)
    n_gcc_strands = int(classes_gcc.size)
    n_all_strands = int(classes_all.size)
    gcc_class_fracs = np.array(
        [
            float(np.sum(classes_gcc == int(cls))) / max(n_gcc_strands, 1)
            for cls in StrandTopologyClass
        ],
        dtype=np.float64,
    )
    all_class_counts = np.array(
        [int(np.sum(classes_all == int(cls))) for cls in StrandTopologyClass],
        dtype=np.float64,
    )
    bridge_mask = classes_gcc == int(StrandTopologyClass.BRIDGE_CUT_EDGE)
    bridge_max_gamma = (
        float(np.max(gamma_eq_gcc[bridge_mask])) if np.any(bridge_mask) else 0.0
    )

    # Correlations on active PBC-spanning backbone strands in GCC
    active_mask = classes_gcc == int(StrandTopologyClass.PBC_ACTIVE_BACKBONE)
    desc_arrays: Dict[TopologicalDescriptor, np.ndarray] = {
        TopologicalDescriptor.LEVERAGE: leverage_gcc[active_mask],
        TopologicalDescriptor.STIFFNESS: w_gcc[active_mask],
        TopologicalDescriptor.MEAN_CENTRALITY: mean_cent_gcc[active_mask],
    }
    obs_arrays: Dict[MechanicalObservable, np.ndarray] = {
        MechanicalObservable.GAMMA_EQ: gamma_eq_gcc[active_mask],
        MechanicalObservable.DELTA_GAMMA: dgamma_gcc[active_mask],
        MechanicalObservable.FORCE_EQ: force_eq_gcc[active_mask],
        MechanicalObservable.DELTA_FORCE: dforce_gcc[active_mask],
    }

    spearman_mat = np.full(
        (len(_DESCRIPTORS), len(_OBSERVABLES)), np.nan, dtype=np.float64
    )
    pearson_mat = np.full(
        (len(_DESCRIPTORS), len(_OBSERVABLES)), np.nan, dtype=np.float64
    )
    for i_d, desc in enumerate(_DESCRIPTORS):
        for j_o, obs in enumerate(_OBSERVABLES):
            rho, r_val = _safe_correlation(desc_arrays[desc], obs_arrays[obs])
            spearman_mat[i_d, j_o] = rho
            pearson_mat[i_d, j_o] = r_val

    # Fluctuation-suppression sum rule check: E[gamma_s] vs 1 - w_s R_AB
    mean_gamma_active = float(np.mean(gamma_eq_gcc[active_mask]))
    mean_one_minus_lev_active = float(np.mean(1.0 - leverage_gcc[active_mask]))
    mean_sum_active = float(
        np.mean(gamma_eq_gcc[active_mask] + leverage_gcc[active_mask])
    )

    # Binned mean gamma_eq vs leverage_gcc on active PBC-spanning strands
    n_bins = len(LEVERAGE_BIN_EDGES) - 1
    bin_gamma_mean = np.full(n_bins, np.nan, dtype=np.float64)
    bin_lev_mean = np.full(n_bins, np.nan, dtype=np.float64)
    bin_counts = np.zeros(n_bins, dtype=np.float64)
    for b_idx in range(n_bins):
        lo = LEVERAGE_BIN_EDGES[b_idx]
        hi = LEVERAGE_BIN_EDGES[b_idx + 1]
        in_bin = active_mask & (leverage_gcc >= lo) & (leverage_gcc < hi)
        cnt = int(np.sum(in_bin))
        bin_counts[b_idx] = float(cnt)
        if cnt >= 1:
            bin_gamma_mean[b_idx] = float(np.mean(gamma_eq_gcc[in_bin]))
            bin_lev_mean[b_idx] = float(np.mean(leverage_gcc[in_bin]))

    # Node centrality correlations and stratified means
    rho_cent_topo, _ = _safe_correlation(node_cent, node_topo_deg.astype(np.float64))
    rho_cent_active, _ = _safe_correlation(
        node_cent, node_active_deg.astype(np.float64)
    )
    cent_by_active_deg = np.full(len(ACTIVE_DEGREE_BINS), np.nan, dtype=np.float64)
    for d_idx, deg_val in enumerate(ACTIVE_DEGREE_BINS):
        mask_d = node_active_deg == deg_val
        if np.any(mask_d):
            cent_by_active_deg[d_idx] = float(np.mean(node_cent[mask_d]))

    # 3. Entanglement activation via EntanglementMode.MERGED_LINKS
    n_ent_target = int(constants.ent_density * volume)
    net_ent = pt.MEHPNetwork.from_universe_with_random_entanglements(
        universe,
        nr_of_entanglements_to_sample=n_ent_target,
        upper_sampling_cutoff=constants.cutoff,
        seed=f"ent_{regime.value}_{seed}",
        crosslinker_type=2,
        entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
        filter_dangling_and_soluble=True,
        assume_box_large_enough=False,
    )
    sim_ent = pt.MEHPSimulator(net_ent)
    solver_cfg_simp = pt.MEHPLinearSolverConfig()
    solver_cfg_simp.simplification_mode = (
        pt.StructureSimplificationMode.INACTIVE_THEN_X2F
    )
    sim_ent.run(solver_cfg_simp)

    gamma_ent = np.asarray(sim_ent.get_gamma_factors(constants.r02), dtype=np.float64)
    w_act_ent = float(sim_ent.get_active_weight_fraction(ACTIVE_GAMMA_TOL))
    g_ent_mpa = float(constants.gamma_conv * np.sum(gamma_ent) / volume)

    # Count primary loops activated by trapped entanglements
    net_ent_after = sim_ent.get_network()
    n_prim_ph_total = int(np.sum(classes_all == int(StrandTopologyClass.PRIMARY_LOOP)))
    n_prim_ent_active = 0
    for s_idx, link_seq in enumerate(net_ent_after.link_indices_of_strand):
        if len(link_seq) >= 2 and link_seq[0] == link_seq[-1]:
            sp_idxs = net_ent_after.spring_indices_of_strand[s_idx]
            if any(gamma_ent[si] > ACTIVE_GAMMA_TOL for si in sp_idxs):
                n_prim_ent_active += 1
    prim_activation_frac = (
        float(n_prim_ent_active) / float(n_prim_ph_total)
        if n_prim_ph_total > 0
        else float("nan")
    )

    result: Dict[str, Any] = {
        "regime": regime.value,
        "seed": seed,
        "n_all_strands": float(n_all_strands),
        "n_gcc_strands": float(n_gcc_strands),
        "n_gcc_nodes": float(n_gcc),
        "n_total_nodes": float(net_ph.nr_of_links),
        "gcc_class_fracs": gcc_class_fracs,
        "all_class_counts": all_class_counts,
        "bridge_max_gamma": bridge_max_gamma,
        "spearman_mat": spearman_mat,
        "pearson_mat": pearson_mat,
        "mean_gamma_active": mean_gamma_active,
        "mean_one_minus_lev_active": mean_one_minus_lev_active,
        "mean_sum_active": mean_sum_active,
        "bin_gamma_mean": bin_gamma_mean,
        "bin_lev_mean": bin_lev_mean,
        "bin_counts": bin_counts,
        "rho_cent_topo": rho_cent_topo,
        "rho_cent_active": rho_cent_active,
        "cent_by_active_deg": cent_by_active_deg,
        "w_act_ph": w_act_ph,
        "w_act_ent": w_act_ent,
        "g_ph_mpa": g_ph_mpa,
        "g_ent_mpa": g_ent_mpa,
        "prim_activation_frac": prim_activation_frac,
        "n_prim_ph_total": float(n_prim_ph_total),
        "n_prim_ent_active": float(n_prim_ent_active),
        "nonzero_offset_frac_active": float(
            np.mean(offset_norm_gcc[active_mask] > 1e-6)
        )
        if np.any(active_mask)
        else 0.0,
    }

    if seed == REPRESENTATIVE_SEED:
        result["rep_leverage"] = leverage_gcc
        result["rep_r_ab"] = r_ab_gcc
        result["rep_w_s"] = w_gcc
        result["rep_n_s"] = n_s_gcc
        result["rep_gamma_eq"] = gamma_eq_gcc
        result["rep_gamma_def"] = gamma_def_gcc
        result["rep_dgamma"] = dgamma_gcc
        result["rep_force_eq"] = force_eq_gcc
        result["rep_force_def"] = force_def_gcc
        result["rep_dforce"] = dforce_gcc
        result["rep_offset_norm"] = offset_norm_gcc
        result["rep_classes"] = classes_gcc
        result["rep_mean_cent"] = mean_cent_gcc
        result["rep_node_cent"] = node_cent
        result["rep_node_topo_deg"] = node_topo_deg
        result["rep_node_active_deg"] = node_active_deg

    return result


def _nan_ensemble_stats(
    arr: np.ndarray,
    axis: int = 0,
    confidence: float = 0.95,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute (mean, 95% Student-t CI half-width) ignoring NaNs along axis."""
    data = np.asarray(arr, dtype=np.float64)
    if not np.any(np.isnan(data)):
        st = compute_ensemble_statistics(data, axis=axis, confidence=confidence)
        return np.asarray(st["mean"], dtype=np.float64), np.asarray(
            st["ci"], dtype=np.float64
        )
    moved = np.moveaxis(data, axis, 0)
    out_shape = moved.shape[1:]
    flat = moved.reshape(moved.shape[0], -1)
    means = np.full(flat.shape[1], np.nan, dtype=np.float64)
    cis = np.full(flat.shape[1], np.nan, dtype=np.float64)
    for col in range(flat.shape[1]):
        valid = flat[:, col]
        valid = valid[~np.isnan(valid)]
        if valid.size == 0:
            continue
        means[col] = float(np.mean(valid))
        if valid.size > 1:
            std = float(np.std(valid, ddof=1))
            sem = std / np.sqrt(float(valid.size))
            t_crit = float(stats.t.ppf((1.0 + confidence) / 2.0, df=valid.size - 1))
            cis[col] = t_crit * sem
        else:
            cis[col] = 0.0
    return means.reshape(out_shape), cis.reshape(out_shape)


def run_resistance_mehp_study(
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Run the full 8-seed Resistance Distance vs. MEHP Stress study across both regimes.

    Parameters
    ----------
    max_workers : int, optional
        Number of worker processes for ``ProcessPoolExecutor``.

    Returns
    -------
    dict
        Dictionary of numpy arrays suitable for ``np.savez_compressed``.
    """
    constants = _load_pdms_constants()
    workers = resolve_max_workers(max_workers)
    configure_process_threading(workers=workers, verbose=True)

    tasks: List[Tuple[str, int, _PdmsConstants]] = [
        (cfg.regime.value, seed, constants) for cfg in _REGIME_CONFIGS for seed in SEEDS
    ]

    results_by_key: Dict[Tuple[str, int], Dict[str, Any]] = {}
    if workers <= 1:
        for task in tasks:
            res = _run_single_task(task)
            results_by_key[(res["regime"], int(res["seed"]))] = res
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as pool:
            futures = [pool.submit(_run_single_task, task) for task in tasks]
            for fut in as_completed(futures):
                res = fut.result()
                results_by_key[(res["regime"], int(res["seed"]))] = res

    cache: Dict[str, Any] = {
        "seeds": np.asarray(SEEDS, dtype=np.int32),
        "deformation_lambda": np.array(DEFORMATION_LAMBDA, dtype=np.float64),
        "leverage_bin_edges": np.asarray(LEVERAGE_BIN_EDGES, dtype=np.float64),
        "active_degree_bins": np.asarray(ACTIVE_DEGREE_BINS, dtype=np.int32),
        "regime_names": np.asarray([c.regime.value for c in _REGIME_CONFIGS]),
        "regime_labels": np.asarray([c.label for c in _REGIME_CONFIGS]),
        "descriptor_names": np.asarray([d.value for d in _DESCRIPTORS]),
        "observable_names": np.asarray([o.value for o in _OBSERVABLES]),
        "class_names": np.asarray([c.name for c in StrandTopologyClass]),
        "r02": np.array(constants.r02, dtype=np.float64),
        "gamma_conv_mpa": np.array(constants.gamma_conv, dtype=np.float64),
    }

    scalar_keys: Tuple[str, ...] = (
        "n_all_strands",
        "n_gcc_strands",
        "n_gcc_nodes",
        "n_total_nodes",
        "bridge_max_gamma",
        "mean_gamma_active",
        "mean_one_minus_lev_active",
        "mean_sum_active",
        "rho_cent_topo",
        "rho_cent_active",
        "w_act_ph",
        "w_act_ent",
        "g_ph_mpa",
        "g_ent_mpa",
        "prim_activation_frac",
        "n_prim_ph_total",
        "n_prim_ent_active",
        "nonzero_offset_frac_active",
    )
    array_keys: Tuple[str, ...] = (
        "gcc_class_fracs",
        "all_class_counts",
        "spearman_mat",
        "pearson_mat",
        "bin_gamma_mean",
        "bin_lev_mean",
        "bin_counts",
        "cent_by_active_deg",
    )

    for cfg in _REGIME_CONFIGS:
        prefix = cfg.regime.value
        seed_rows = [results_by_key[(prefix, s)] for s in SEEDS]

        for key in scalar_keys:
            vals = np.asarray([row[key] for row in seed_rows], dtype=np.float64)
            mean_v, ci_v = _nan_ensemble_stats(vals, axis=0)
            cache[f"{prefix}_{key}_seeds"] = vals
            cache[f"{prefix}_{key}_mean"] = mean_v
            cache[f"{prefix}_{key}_ci"] = ci_v

        for key in array_keys:
            stacked = np.stack([row[key] for row in seed_rows], axis=0)
            mean_a, ci_a = _nan_ensemble_stats(stacked, axis=0)
            cache[f"{prefix}_{key}_seeds"] = stacked
            cache[f"{prefix}_{key}_mean"] = mean_a
            cache[f"{prefix}_{key}_ci"] = ci_a

        rep = results_by_key[(prefix, REPRESENTATIVE_SEED)]
        for rep_key in (
            "rep_leverage",
            "rep_r_ab",
            "rep_w_s",
            "rep_n_s",
            "rep_gamma_eq",
            "rep_gamma_def",
            "rep_dgamma",
            "rep_force_eq",
            "rep_force_def",
            "rep_dforce",
            "rep_offset_norm",
            "rep_classes",
            "rep_mean_cent",
            "rep_node_cent",
            "rep_node_topo_deg",
            "rep_node_active_deg",
        ):
            cache[f"{prefix}_{rep_key}"] = np.asarray(rep[rep_key])

    return cache


def _parse_cli_workers(argv: Sequence[str]) -> Optional[int]:
    """Parse optional --workers / -j argument from CLI."""
    for i, arg in enumerate(argv):
        if arg in ("--workers", "-j") and i + 1 < len(argv):
            return int(argv[i + 1])
        if arg.startswith("--workers="):
            return int(arg.split("=", 1)[1])
    return None


if __name__ == "__main__":
    cache_path = get_example_cache_path(
        "resistance_vs_mehp_stress_cache.npz", example_folder="structure_analysis"
    )
    cli_workers = _parse_cli_workers(sys.argv[1:])
    load_or_compute_cache(
        cache_path,
        run_resistance_mehp_study,
        recompute="--recompute" in sys.argv,
        max_workers=cli_workers,
    )
