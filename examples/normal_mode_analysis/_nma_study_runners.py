"""
Simulation runners for Normal Mode Analysis (KPM vs. Exact) examples on PDMS networks.

Provides parallelized, cache-friendly study functions comparing the Kernel Polynomial
Method (KPM) against full eigenvalue decomposition across chain length N,
cross-link conversion p, stoichiometric imbalance r, and representative volume /
sample count.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from enum import Enum
import math
from pathlib import Path
import sys
import time
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

# Ensure examples/ is importable for _example_cache
_EXAMPLES_DIR = str(Path(__file__).resolve().parents[1])
if _EXAMPLES_DIR not in sys.path:
    sys.path.insert(0, _EXAMPLES_DIR)

from _example_cache import (  # noqa: E402
    compute_ensemble_statistics,
    configure_process_threading,
    resolve_max_workers,
    worker_process_initializer,
)
from pylimer_tools.io.bead_spring_parameter_provider import (  # noqa: E402
    ParameterType,
    Parameters,
    get_parameters_for_polymer,
)
import pylimer_tools_cpp as pt  # noqa: E402


class _NmaStudyKind(Enum):
    """Sub-study identifiers for parallel task dispatch."""

    CHAIN_LENGTH = "chain_length"
    CONVERSION = "conversion"
    STOICHIOMETRY = "stoichiometry"
    RVE = "rve"


def _build_pdms_nma(
    box_l: float,
    strand_len: int,
    conversion: float,
    stoichiometry: float,
    seed: int,
    b: float,
    b2: float,
    cutoff: float,
    ent_density: float,
    bead_density: float,
) -> Tuple[pt.NormalModeAnalyzer, int, int, float]:
    """Generate a PDMS network via MCUniverseGenerator with on-demand entanglements
    and construct a NormalModeAnalyzer."""
    vol = box_l**3
    n_strands = max(
        1, int(round(bead_density * vol / (strand_len + 0.5 * stoichiometry)))
    )
    n_xlinks = max(1, int(round(0.5 * stoichiometry * n_strands)))
    max_reachable_p = (2.0 * n_strands) / (4.0 * n_xlinks)
    effective_conversion = min(float(conversion), max_reachable_p)

    gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(int(seed))
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(b)
    gen.set_mean_squared_bead_distance(b2)
    gen.add_strands(n_strands, [int(strand_len)] * n_strands)
    gen.add_crosslinkers(n_xlinks, 4)
    gen.link_strands_to_conversion(effective_conversion)

    es = gen.sample_entanglements(
        nr_of_entanglements_to_sample=int(ent_density * vol),
        upper_cutoff=cutoff,
        seed=f"ent_{seed}",
        filtering_mode=pt.EntanglementFilteringMode.OMIT_INACTIVE,
    )
    atom_pairs = es.to_atom_pair_entanglements(gen.get_crosslinker_universe())

    univ = gen.get_universe()
    pairs = atom_pairs.pairs_of_atoms
    if pairs:
        univ.add_bonds([p[0] for p in pairs], [p[1] for p in pairs])

    edges = univ.get_edges()
    n_atoms = univ.get_nr_of_atoms()
    nma = pt.NormalModeAnalyzer(edges["edge_from"], edges["edge_to"], n_atoms)
    return nma, n_atoms, len(pairs), vol


def _evaluate_single_nma_task(
    kind: _NmaStudyKind,
    param_idx: int,
    rep_idx: int,
    seed: int,
    box_l: float,
    strand_len: int,
    conversion: float,
    stoichiometry: float,
    freqs: np.ndarray,
    compute_exact: bool,
    num_moments: int,
    num_random_vectors: int,
    num_quadrature_points: int,
    num_lanczos_steps: int,
    num_exact_low_modes: int,
    b: float,
    b2: float,
    cutoff: float,
    ent_density: float,
    bead_density: float,
    gamma_conv: float,
) -> Tuple[
    _NmaStudyKind,
    int,
    int,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
    float,
    float,
    float,
    int,
    int,
]:
    """Execute one MC PDMS network generation, entanglement sampling, and NMA evaluation."""
    nma, n_atoms, n_ent, vol = _build_pdms_nma(
        box_l=box_l,
        strand_len=strand_len,
        conversion=conversion,
        stoichiometry=stoichiometry,
        seed=seed,
        b=b,
        b2=b2,
        cutoff=cutoff,
        ent_density=ent_density,
        bead_density=bead_density,
    )
    scale = gamma_conv / vol

    t0_kpm = time.perf_counter()
    gp_kpm = (
        nma.evaluate_storage_modulus_kpm(
            freqs,
            num_moments=num_moments,
            num_random_vectors=num_random_vectors,
            num_quadrature_points=num_quadrature_points,
        )
        * scale
    )
    gpp_kpm = (
        nma.evaluate_loss_modulus_kpm(
            freqs,
            num_moments=num_moments,
            num_random_vectors=num_random_vectors,
            num_quadrature_points=num_quadrature_points,
        )
        * scale
    )
    t_kpm = time.perf_counter() - t0_kpm

    t0_slq = time.perf_counter()
    gp_slq = (
        nma.evaluate_storage_modulus_approx(
            freqs,
            method=pt.StochasticSpectralMethod.SLQ,
            num_steps=num_lanczos_steps,
            num_random_vectors=num_random_vectors,
        )
        * scale
    )
    gpp_slq = (
        nma.evaluate_loss_modulus_approx(
            freqs,
            method=pt.StochasticSpectralMethod.SLQ,
            num_steps=num_lanczos_steps,
            num_random_vectors=num_random_vectors,
        )
        * scale
    )
    t_slq = time.perf_counter() - t0_slq

    t0_hybrid = time.perf_counter()
    gp_hybrid = (
        nma.evaluate_storage_modulus_approx(
            freqs,
            method=pt.StochasticSpectralMethod.HYBRID,
            num_steps=num_lanczos_steps,
            num_random_vectors=num_random_vectors,
            num_exact_low_modes=num_exact_low_modes,
        )
        * scale
    )
    gpp_hybrid = (
        nma.evaluate_loss_modulus_approx(
            freqs,
            method=pt.StochasticSpectralMethod.HYBRID,
            num_steps=num_lanczos_steps,
            num_random_vectors=num_random_vectors,
            num_exact_low_modes=num_exact_low_modes,
        )
        * scale
    )
    t_hybrid = time.perf_counter() - t0_hybrid

    if compute_exact:
        t0_ex = time.perf_counter()
        nma.find_all_eigenvalues(
            solver=pt.NormalModeSolver.AUTO,
            compute_eigenvectors=False,
            decompose_components=True,
        )
        gp_ex = nma.evaluate_storage_modulus(freqs) * scale
        gpp_ex = nma.evaluate_loss_modulus(freqs) * scale
        t_ex = time.perf_counter() - t0_ex
    else:
        gp_ex = np.full_like(freqs, np.nan, dtype=float)
        gpp_ex = np.full_like(freqs, np.nan, dtype=float)
        t_ex = float("nan")

    return (
        kind,
        param_idx,
        rep_idx,
        gp_ex,
        gpp_ex,
        gp_kpm,
        gpp_kpm,
        gp_slq,
        gpp_slq,
        gp_hybrid,
        gpp_hybrid,
        t_ex,
        t_kpm,
        t_slq,
        t_hybrid,
        n_atoms,
        n_ent,
    )


def run_kpm_vs_exact_pdms_study(
    n_values: Sequence[int] = (50, 100, 200, 400),
    p_values: Sequence[float] = (0.75, 0.85, 0.925, 0.98),
    r_values: Sequence[float] = (0.65, 0.80, 1.00, 1.25),
    box_l_param_studies: float = 14.0,
    rve_box_lengths: Sequence[float] = (10.0, 12.0, 14.0, 16.0, 20.0, 28.0, 36.0),
    max_exact_box_l: float = 20.0,
    n_repeats: int = 8,
    n_rve_repeats: int = 12,
    num_moments: int = 1024,
    num_random_vectors: int = 20,
    num_quadrature_points: int = 2048,
    num_lanczos_steps: int = 200,
    num_exact_low_modes: int = 100,
    base_seed: int = 4200,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Run the KPM, SLQ, Hybrid vs. Exact normal mode analysis comparison for realistic PDMS networks.

    Evaluates:
      1. Chain length N sweep (at p=0.95, r=1.0, L=14 nm)
      2. Cross-link conversion p sweep (at N=100, r=1.0, L=14 nm)
      3. Stoichiometric imbalance r sweep (at N=100, p=0.95*min(1, 1/r), L=14 nm)
      4. Representative volume (L in [10, 36] nm) & sample count (n_rve_repeats=12)
         with exact eigenvalues computed up to L=20 nm and KPM/SLQ/Hybrid across all box sizes.
    """
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    assert isinstance(params, Parameters)
    dist_unit = params.get("distance_units")
    b = float(params.get("<b>").to(dist_unit).magnitude)
    b2 = float(params.get("<b^2>").to(dist_unit**2).magnitude)
    cutoff = float(params.get_sampling_cutoff())
    ent_density = float(params.get_entanglement_density())
    bead_density = float(params.get_bead_density())
    gamma_conv = float(params.get_gamma_conversion_factor().to("MPa").magnitude)

    freqs = np.logspace(-4, 2, 75)
    n_freqs = len(freqs)
    r_conversions = [0.95 * min(1.0, 1.0 / r_val) for r_val in r_values]

    # Preallocate arrays
    n_gp_ex = np.zeros((len(n_values), n_repeats, n_freqs))
    n_gpp_ex = np.zeros((len(n_values), n_repeats, n_freqs))
    n_gp_kpm = np.zeros((len(n_values), n_repeats, n_freqs))
    n_gpp_kpm = np.zeros((len(n_values), n_repeats, n_freqs))
    n_gp_slq = np.zeros((len(n_values), n_repeats, n_freqs))
    n_gpp_slq = np.zeros((len(n_values), n_repeats, n_freqs))
    n_gp_hybrid = np.zeros((len(n_values), n_repeats, n_freqs))
    n_gpp_hybrid = np.zeros((len(n_values), n_repeats, n_freqs))

    p_gp_ex = np.zeros((len(p_values), n_repeats, n_freqs))
    p_gpp_ex = np.zeros((len(p_values), n_repeats, n_freqs))
    p_gp_kpm = np.zeros((len(p_values), n_repeats, n_freqs))
    p_gpp_kpm = np.zeros((len(p_values), n_repeats, n_freqs))
    p_gp_slq = np.zeros((len(p_values), n_repeats, n_freqs))
    p_gpp_slq = np.zeros((len(p_values), n_repeats, n_freqs))
    p_gp_hybrid = np.zeros((len(p_values), n_repeats, n_freqs))
    p_gpp_hybrid = np.zeros((len(p_values), n_repeats, n_freqs))

    r_gp_ex = np.zeros((len(r_values), n_repeats, n_freqs))
    r_gpp_ex = np.zeros((len(r_values), n_repeats, n_freqs))
    r_gp_kpm = np.zeros((len(r_values), n_repeats, n_freqs))
    r_gpp_kpm = np.zeros((len(r_values), n_repeats, n_freqs))
    r_gp_slq = np.zeros((len(r_values), n_repeats, n_freqs))
    r_gpp_slq = np.zeros((len(r_values), n_repeats, n_freqs))
    r_gp_hybrid = np.zeros((len(r_values), n_repeats, n_freqs))
    r_gpp_hybrid = np.zeros((len(r_values), n_repeats, n_freqs))

    rve_gp_ex = np.full((len(rve_box_lengths), n_rve_repeats, n_freqs), np.nan)
    rve_gpp_ex = np.full((len(rve_box_lengths), n_rve_repeats, n_freqs), np.nan)
    rve_gp_kpm = np.zeros((len(rve_box_lengths), n_rve_repeats, n_freqs))
    rve_gpp_kpm = np.zeros((len(rve_box_lengths), n_rve_repeats, n_freqs))
    rve_gp_slq = np.zeros((len(rve_box_lengths), n_rve_repeats, n_freqs))
    rve_gpp_slq = np.zeros((len(rve_box_lengths), n_rve_repeats, n_freqs))
    rve_gp_hybrid = np.zeros((len(rve_box_lengths), n_rve_repeats, n_freqs))
    rve_gpp_hybrid = np.zeros((len(rve_box_lengths), n_rve_repeats, n_freqs))
    rve_t_ex = np.full((len(rve_box_lengths), n_rve_repeats), np.nan)
    rve_t_kpm = np.zeros((len(rve_box_lengths), n_rve_repeats))
    rve_t_slq = np.zeros((len(rve_box_lengths), n_rve_repeats))
    rve_t_hybrid = np.zeros((len(rve_box_lengths), n_rve_repeats))
    rve_n_atoms = np.zeros((len(rve_box_lengths), n_rve_repeats), dtype=int)
    rve_n_ent = np.zeros((len(rve_box_lengths), n_rve_repeats), dtype=int)

    tasks = []
    common_tail = (
        num_moments,
        num_random_vectors,
        num_quadrature_points,
        num_lanczos_steps,
        num_exact_low_modes,
        b,
        b2,
        cutoff,
        ent_density,
        bead_density,
        gamma_conv,
    )

    # Study 1: Chain length N
    for idx, n_len in enumerate(n_values):
        for rep in range(n_repeats):
            tasks.append(
                (
                    _NmaStudyKind.CHAIN_LENGTH,
                    idx,
                    rep,
                    base_seed + 10000 * idx + 313 * rep,
                    box_l_param_studies,
                    int(n_len),
                    0.95,
                    1.0,
                    freqs,
                    True,
                    *common_tail,
                )
            )

    # Study 2: Conversion p
    for idx, p_val in enumerate(p_values):
        for rep in range(n_repeats):
            tasks.append(
                (
                    _NmaStudyKind.CONVERSION,
                    idx,
                    rep,
                    base_seed + 100000 + 10000 * idx + 313 * rep,
                    box_l_param_studies,
                    100,
                    float(p_val),
                    1.0,
                    freqs,
                    True,
                    *common_tail,
                )
            )

    # Study 3: Stoichiometric imbalance r
    for idx, (r_val, conv_val) in enumerate(zip(r_values, r_conversions)):
        for rep in range(n_repeats):
            tasks.append(
                (
                    _NmaStudyKind.STOICHIOMETRY,
                    idx,
                    rep,
                    base_seed + 200000 + 10000 * idx + 313 * rep,
                    box_l_param_studies,
                    100,
                    float(conv_val),
                    float(r_val),
                    freqs,
                    True,
                    *common_tail,
                )
            )

    # Study 4: RVE & sample count
    for idx, box_l in enumerate(rve_box_lengths):
        do_exact = float(box_l) <= max_exact_box_l + 1e-9
        for rep in range(n_rve_repeats):
            tasks.append(
                (
                    _NmaStudyKind.RVE,
                    idx,
                    rep,
                    base_seed + 300000 + 10000 * idx + 313 * rep,
                    float(box_l),
                    100,
                    0.95,
                    1.0,
                    freqs,
                    do_exact,
                    *common_tail,
                )
            )

    workers = min(resolve_max_workers(max_workers), len(tasks))
    configure_process_threading(workers)
    total_tasks = len(tasks)

    def _store_result(
        res: Tuple[
            _NmaStudyKind,
            int,
            int,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            float,
            float,
            float,
            float,
            int,
            int,
        ],
    ) -> None:
        (
            kind,
            p_idx,
            rep_idx,
            gp_ex,
            gpp_ex,
            gp_kpm,
            gpp_kpm,
            gp_slq,
            gpp_slq,
            gp_hybrid,
            gpp_hybrid,
            t_ex,
            t_kpm,
            t_slq,
            t_hybrid,
            n_at,
            n_ent,
        ) = res
        if kind == _NmaStudyKind.CHAIN_LENGTH:
            n_gp_ex[p_idx, rep_idx] = gp_ex
            n_gpp_ex[p_idx, rep_idx] = gpp_ex
            n_gp_kpm[p_idx, rep_idx] = gp_kpm
            n_gpp_kpm[p_idx, rep_idx] = gpp_kpm
            n_gp_slq[p_idx, rep_idx] = gp_slq
            n_gpp_slq[p_idx, rep_idx] = gpp_slq
            n_gp_hybrid[p_idx, rep_idx] = gp_hybrid
            n_gpp_hybrid[p_idx, rep_idx] = gpp_hybrid
        elif kind == _NmaStudyKind.CONVERSION:
            p_gp_ex[p_idx, rep_idx] = gp_ex
            p_gpp_ex[p_idx, rep_idx] = gpp_ex
            p_gp_kpm[p_idx, rep_idx] = gp_kpm
            p_gpp_kpm[p_idx, rep_idx] = gpp_kpm
            p_gp_slq[p_idx, rep_idx] = gp_slq
            p_gpp_slq[p_idx, rep_idx] = gpp_slq
            p_gp_hybrid[p_idx, rep_idx] = gp_hybrid
            p_gpp_hybrid[p_idx, rep_idx] = gpp_hybrid
        elif kind == _NmaStudyKind.STOICHIOMETRY:
            r_gp_ex[p_idx, rep_idx] = gp_ex
            r_gpp_ex[p_idx, rep_idx] = gpp_ex
            r_gp_kpm[p_idx, rep_idx] = gp_kpm
            r_gpp_kpm[p_idx, rep_idx] = gpp_kpm
            r_gp_slq[p_idx, rep_idx] = gp_slq
            r_gpp_slq[p_idx, rep_idx] = gpp_slq
            r_gp_hybrid[p_idx, rep_idx] = gp_hybrid
            r_gpp_hybrid[p_idx, rep_idx] = gpp_hybrid
        elif kind == _NmaStudyKind.RVE:
            rve_gp_ex[p_idx, rep_idx] = gp_ex
            rve_gpp_ex[p_idx, rep_idx] = gpp_ex
            rve_gp_kpm[p_idx, rep_idx] = gp_kpm
            rve_gpp_kpm[p_idx, rep_idx] = gpp_kpm
            rve_gp_slq[p_idx, rep_idx] = gp_slq
            rve_gpp_slq[p_idx, rep_idx] = gpp_slq
            rve_gp_hybrid[p_idx, rep_idx] = gp_hybrid
            rve_gpp_hybrid[p_idx, rep_idx] = gpp_hybrid
            rve_t_ex[p_idx, rep_idx] = t_ex
            rve_t_kpm[p_idx, rep_idx] = t_kpm
            rve_t_slq[p_idx, rep_idx] = t_slq
            rve_t_hybrid[p_idx, rep_idx] = t_hybrid
            rve_n_atoms[p_idx, rep_idx] = n_at
            rve_n_ent[p_idx, rep_idx] = n_ent

    if workers > 1:
        print(
            f"Executing {total_tasks} PDMS NMA (KPM vs. Exact) tasks across {workers} workers..."
        )
        completed = 0
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = [executor.submit(_evaluate_single_nma_task, *t) for t in tasks]
            for fut in as_completed(futures):
                _store_result(fut.result())
                completed += 1
                if completed % 20 == 0 or completed == total_tasks:
                    print(
                        f"  [PDMS NMA Study] Completed {completed}/{total_tasks} tasks..."
                    )
    else:
        for idx, t in enumerate(tasks, start=1):
            _store_result(_evaluate_single_nma_task(*t))
            if idx % 20 == 0 or idx == total_tasks:
                print(f"  [PDMS NMA Study] Completed {idx}/{total_tasks} tasks...")

    stat_n_gp_ex = compute_ensemble_statistics(n_gp_ex, axis=1)
    stat_n_gpp_ex = compute_ensemble_statistics(n_gpp_ex, axis=1)
    stat_n_gp_kpm = compute_ensemble_statistics(n_gp_kpm, axis=1)
    stat_n_gpp_kpm = compute_ensemble_statistics(n_gpp_kpm, axis=1)
    stat_n_gp_slq = compute_ensemble_statistics(n_gp_slq, axis=1)
    stat_n_gpp_slq = compute_ensemble_statistics(n_gpp_slq, axis=1)
    stat_n_gp_hybrid = compute_ensemble_statistics(n_gp_hybrid, axis=1)
    stat_n_gpp_hybrid = compute_ensemble_statistics(n_gpp_hybrid, axis=1)

    stat_p_gp_ex = compute_ensemble_statistics(p_gp_ex, axis=1)
    stat_p_gpp_ex = compute_ensemble_statistics(p_gpp_ex, axis=1)
    stat_p_gp_kpm = compute_ensemble_statistics(p_gp_kpm, axis=1)
    stat_p_gpp_kpm = compute_ensemble_statistics(p_gpp_kpm, axis=1)
    stat_p_gp_slq = compute_ensemble_statistics(p_gp_slq, axis=1)
    stat_p_gpp_slq = compute_ensemble_statistics(p_gpp_slq, axis=1)
    stat_p_gp_hybrid = compute_ensemble_statistics(p_gp_hybrid, axis=1)
    stat_p_gpp_hybrid = compute_ensemble_statistics(p_gpp_hybrid, axis=1)

    stat_r_gp_ex = compute_ensemble_statistics(r_gp_ex, axis=1)
    stat_r_gpp_ex = compute_ensemble_statistics(r_gpp_ex, axis=1)
    stat_r_gp_kpm = compute_ensemble_statistics(r_gp_kpm, axis=1)
    stat_r_gpp_kpm = compute_ensemble_statistics(r_gpp_kpm, axis=1)
    stat_r_gp_slq = compute_ensemble_statistics(r_gp_slq, axis=1)
    stat_r_gpp_slq = compute_ensemble_statistics(r_gpp_slq, axis=1)
    stat_r_gp_hybrid = compute_ensemble_statistics(r_gp_hybrid, axis=1)
    stat_r_gpp_hybrid = compute_ensemble_statistics(r_gpp_hybrid, axis=1)

    stat_rve_gp_ex = compute_ensemble_statistics(rve_gp_ex, axis=1)
    stat_rve_gp_kpm = compute_ensemble_statistics(rve_gp_kpm, axis=1)
    stat_rve_gpp_kpm = compute_ensemble_statistics(rve_gpp_kpm, axis=1)
    stat_rve_gp_slq = compute_ensemble_statistics(rve_gp_slq, axis=1)
    stat_rve_gpp_slq = compute_ensemble_statistics(rve_gpp_slq, axis=1)
    stat_rve_gp_hybrid = compute_ensemble_statistics(rve_gp_hybrid, axis=1)
    stat_rve_gpp_hybrid = compute_ensemble_statistics(rve_gpp_hybrid, axis=1)

    n_r0_nm = np.array([math.sqrt(n_len * b2) for n_len in n_values])

    return {
        "freqs": freqs,
        "box_l_param_studies": np.array(box_l_param_studies),
        "n_repeats": np.array(n_repeats),
        "n_rve_repeats": np.array(n_rve_repeats),
        "num_moments": np.array(num_moments),
        "num_random_vectors": np.array(num_random_vectors),
        "num_quadrature_points": np.array(num_quadrature_points),
        "num_lanczos_steps": np.array(num_lanczos_steps),
        "num_exact_low_modes": np.array(num_exact_low_modes),
        # Study 1: N
        "n_values": np.array(n_values, dtype=int),
        "n_r0_nm": n_r0_nm,
        "n_gp_ex_mean": stat_n_gp_ex["mean"],
        "n_gp_ex_ci": stat_n_gp_ex["ci"],
        "n_gpp_ex_mean": stat_n_gpp_ex["mean"],
        "n_gp_kpm_mean": stat_n_gp_kpm["mean"],
        "n_gp_kpm_ci": stat_n_gp_kpm["ci"],
        "n_gpp_kpm_mean": stat_n_gpp_kpm["mean"],
        "n_gp_slq_mean": stat_n_gp_slq["mean"],
        "n_gp_slq_ci": stat_n_gp_slq["ci"],
        "n_gpp_slq_mean": stat_n_gpp_slq["mean"],
        "n_gp_hybrid_mean": stat_n_gp_hybrid["mean"],
        "n_gp_hybrid_ci": stat_n_gp_hybrid["ci"],
        "n_gpp_hybrid_mean": stat_n_gpp_hybrid["mean"],
        # Study 2: p
        "p_values": np.array(p_values, dtype=float),
        "p_gp_ex_mean": stat_p_gp_ex["mean"],
        "p_gp_ex_ci": stat_p_gp_ex["ci"],
        "p_gpp_ex_mean": stat_p_gpp_ex["mean"],
        "p_gp_kpm_mean": stat_p_gp_kpm["mean"],
        "p_gp_kpm_ci": stat_p_gp_kpm["ci"],
        "p_gpp_kpm_mean": stat_p_gpp_kpm["mean"],
        "p_gp_slq_mean": stat_p_gp_slq["mean"],
        "p_gp_slq_ci": stat_p_gp_slq["ci"],
        "p_gpp_slq_mean": stat_p_gpp_slq["mean"],
        "p_gp_hybrid_mean": stat_p_gp_hybrid["mean"],
        "p_gp_hybrid_ci": stat_p_gp_hybrid["ci"],
        "p_gpp_hybrid_mean": stat_p_gpp_hybrid["mean"],
        # Study 3: r
        "r_values": np.array(r_values, dtype=float),
        "r_conversions": np.array(r_conversions, dtype=float),
        "r_gp_ex_mean": stat_r_gp_ex["mean"],
        "r_gp_ex_ci": stat_r_gp_ex["ci"],
        "r_gpp_ex_mean": stat_r_gpp_ex["mean"],
        "r_gp_kpm_mean": stat_r_gp_kpm["mean"],
        "r_gp_kpm_ci": stat_r_gp_kpm["ci"],
        "r_gpp_kpm_mean": stat_r_gpp_kpm["mean"],
        "r_gp_slq_mean": stat_r_gp_slq["mean"],
        "r_gp_slq_ci": stat_r_gp_slq["ci"],
        "r_gpp_slq_mean": stat_r_gpp_slq["mean"],
        "r_gp_hybrid_mean": stat_r_gp_hybrid["mean"],
        "r_gp_hybrid_ci": stat_r_gp_hybrid["ci"],
        "r_gpp_hybrid_mean": stat_r_gpp_hybrid["mean"],
        # Study 4: RVE, samples & scaling
        "rve_box_lengths": np.array(rve_box_lengths, dtype=float),
        "rve_gp_ex_mean": stat_rve_gp_ex["mean"],
        "rve_gp_kpm_mean": stat_rve_gp_kpm["mean"],
        "rve_gp_kpm_std": stat_rve_gp_kpm["std"],
        "rve_gp_kpm_ci": stat_rve_gp_kpm["ci"],
        "rve_gpp_kpm_mean": stat_rve_gpp_kpm["mean"],
        "rve_gp_kpm_all": rve_gp_kpm,
        "rve_gp_slq_mean": stat_rve_gp_slq["mean"],
        "rve_gp_slq_ci": stat_rve_gp_slq["ci"],
        "rve_gpp_slq_mean": stat_rve_gpp_slq["mean"],
        "rve_gp_hybrid_mean": stat_rve_gp_hybrid["mean"],
        "rve_gp_hybrid_ci": stat_rve_gp_hybrid["ci"],
        "rve_gpp_hybrid_mean": stat_rve_gpp_hybrid["mean"],
        "rve_t_ex": rve_t_ex,
        "rve_t_kpm": rve_t_kpm,
        "rve_t_slq": rve_t_slq,
        "rve_t_hybrid": rve_t_hybrid,
        "rve_n_atoms": rve_n_atoms,
        "rve_n_ent": rve_n_ent,
    }
