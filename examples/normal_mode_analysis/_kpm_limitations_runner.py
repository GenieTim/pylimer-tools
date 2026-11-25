"""
Simulation runners for Study 3: When KPM Fails in Normal Mode Analysis.

Investigates:
  1. Very-low-frequency breakdown of KPM (and SLQ) vs. exact diagonalization and
     HybridConfig recovery on near-percolation (p >= p_c) and dangling-brush
     PDMS networks.
  2. Density-of-states (DOS) smearing near lambda = 0 via Jackson kernel damping
     and eigenvector localization (Inverse Participation Ratio & atom-type/cluster
     decomposition) of the softest nonzero modes.
  3. The physics boundary between NormalModeAnalyzer (Rouse graph Laplacian) and
     MEHPSimulator (3D force balance with PBC offsets) for uncrosslinked entangled
     melts vs. crosslinked networks and coordinate/PBC-offset invariance.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.sparse.csgraph import connected_components

_EXAMPLES_DIR = str(Path(__file__).resolve().parents[1])
if _EXAMPLES_DIR not in sys.path:
    sys.path.insert(0, _EXAMPLES_DIR)

from _example_cache import (  # noqa: E402
    compute_ensemble_statistics,
    configure_process_threading,
    resolve_max_workers,
    worker_process_initializer,
)
from pylimer_tools.calc.stoichiometric_relationships import (  # noqa: E402
    compute_number_fractions,
)
from pylimer_tools.io.bead_spring_parameter_provider import (  # noqa: E402
    ParameterType,
    Parameters,
    get_parameters_for_polymer,
)
import pylimer_tools_cpp as pt  # noqa: E402


class RegimeKind(Enum):
    """Challenging network regime identifiers."""

    NEAR_PERCOLATION = "near_percolation"
    DANGLING_BRUSH = "dangling_brush"


class SpectralMethod(Enum):
    """Spectral solver configurations evaluated in the frequency sweep."""

    EXACT = "exact"
    KPM_256 = "kpm_256"
    KPM_1024 = "kpm_1024"
    KPM_4096 = "kpm_4096"
    SLQ = "slq"
    HYBRID_50 = "hybrid_50"
    HYBRID_100 = "hybrid_100"


class BoundaryCase(Enum):
    """Network states for the NMA vs. MEHP physics-boundary comparison."""

    MELT_NO_ENT = "melt_no_ent"
    MELT_STATIC_ENT = "melt_static_ent"
    NET_NO_ENT = "net_no_ent"
    NET_STATIC_ENT = "net_static_ent"


@dataclass(frozen=True)
class PdmsParams:
    """Extracted Gaussian PDMS coarse-grained parameters."""

    b: float
    b2: float
    r02: float
    cutoff: float
    ent_density: float
    bead_density: float
    gamma_conv: float


def _load_pdms_params() -> PdmsParams:
    """Load PDMS ParameterType.GAUSSIAN unit-consistent parameters."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    assert isinstance(params, Parameters)
    dist_unit = params.get("distance_units")
    return PdmsParams(
        b=float(params.get("<b>").to(dist_unit).magnitude),
        b2=float(params.get("<b^2>").to(dist_unit**2).magnitude),
        r02=float(params.get("R02").to(dist_unit**2).magnitude),
        cutoff=float(params.get_sampling_cutoff()),
        ent_density=float(params.get_entanglement_density()),
        bead_density=float(params.get_bead_density()),
        gamma_conv=float(params.get_gamma_conversion_factor().to("MPa").magnitude),
    )


def _generate_challenging_network(
    regime: RegimeKind,
    box_l: float,
    seed: int,
    pdms: PdmsParams,
    near_perc_n: int = 50,
    near_perc_p: float = 0.65,
    dangle_bb_n: int = 35,
    dangle_n: int = 150,
    dangle_b2: float = 0.7,
    dangle_p: float = 0.92,
) -> Tuple[pt.MCUniverseGenerator, pt.Universe, pt.NormalModeAnalyzer, int]:
    """Generate a near-percolation or dangling-brush PDMS network with entanglements."""
    vol = box_l**3
    gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(int(seed))
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(pdms.b)
    gen.set_mean_squared_bead_distance(pdms.b2)

    if regime == RegimeKind.NEAR_PERCOLATION:
        n_strands = max(1, int(round(pdms.bead_density * vol / (near_perc_n + 0.5))))
        n_xlinks = max(1, int(round(0.5 * n_strands)))
        gen.add_strands(n_strands, [int(near_perc_n)] * n_strands, 1)
        gen.add_crosslinkers(n_xlinks, 4, 2)
        gen.link_strands_to_conversion(float(near_perc_p))
    elif regime == RegimeKind.DANGLING_BRUSH:
        x_bi, x_mono = compute_number_fractions(dangle_b2)
        beads_per_strand_mix = (
            x_bi * dangle_bb_n + x_mono * dangle_n + (2.0 * x_bi + x_mono) / 4.0
        )
        n_total = max(1, int(round(pdms.bead_density * vol / beads_per_strand_mix)))
        n_bi = max(1, int(round(x_bi * n_total)))
        n_mono = max(1, int(round(x_mono * n_total)))
        n_xlinks = max(1, int(round((2 * n_bi + n_mono) / 4.0)))
        gen.add_strands(n_bi, [int(dangle_bb_n)] * n_bi, 1)
        gen.add_monofunctional_strands(n_mono, [int(dangle_n)] * n_mono, 4)
        gen.add_crosslinkers(n_xlinks, 4, 2)
        gen.link_strands_to_conversion(float(dangle_p))
    else:
        raise ValueError(f"Unsupported RegimeKind: {regime}")

    es = gen.sample_entanglements(
        nr_of_entanglements_to_sample=int(pdms.ent_density * vol),
        upper_cutoff=pdms.cutoff,
        seed=f"ent_{regime.value}_{seed}",
        filtering_mode=pt.EntanglementFilteringMode.OMIT_INACTIVE,
    )
    atom_pairs = es.to_atom_pair_entanglements(gen.get_crosslinker_universe())
    univ = gen.get_universe()
    pairs = atom_pairs.pairs_of_atoms
    if pairs:
        univ.add_bonds([p[0] for p in pairs], [p[1] for p in pairs])

    nma = pt.NormalModeAnalyzer.from_universe(univ)
    return gen, univ, nma, len(pairs)


def _build_approx_config(
    method: SpectralMethod,
    num_random_vectors: int = 20,
    num_lanczos_steps: int = 200,
) -> pt.KPMConfig | pt.SLQConfig | pt.HybridConfig:
    """Map SpectralMethod enum to the corresponding C++ config object."""
    if method == SpectralMethod.KPM_256:
        return pt.KPMConfig(
            num_moments=256,
            num_random_vectors=num_random_vectors,
            num_quadrature_points=2048,
        )
    if method == SpectralMethod.KPM_1024:
        return pt.KPMConfig(
            num_moments=1024,
            num_random_vectors=num_random_vectors,
            num_quadrature_points=2048,
        )
    if method == SpectralMethod.KPM_4096:
        return pt.KPMConfig(
            num_moments=4096,
            num_random_vectors=num_random_vectors,
            num_quadrature_points=8192,
        )
    if method == SpectralMethod.SLQ:
        return pt.SLQConfig(
            num_steps=num_lanczos_steps,
            num_random_vectors=num_random_vectors,
        )
    if method == SpectralMethod.HYBRID_50:
        return pt.HybridConfig(
            num_steps=num_lanczos_steps,
            num_random_vectors=num_random_vectors,
            num_exact_low_modes=50,
        )
    if method == SpectralMethod.HYBRID_100:
        return pt.HybridConfig(
            num_steps=num_lanczos_steps,
            num_random_vectors=num_random_vectors,
            num_exact_low_modes=100,
        )
    raise ValueError(f"Method {method} is not an approximate config.")


def _evaluate_regime_task(
    regime: RegimeKind,
    rep_idx: int,
    seed: int,
    box_l: float,
    freqs: np.ndarray,
    methods: Sequence[SpectralMethod],
    pdms: PdmsParams,
) -> Tuple[
    RegimeKind,
    int,
    Dict[SpectralMethod, np.ndarray],
    Dict[SpectralMethod, np.ndarray],
    Dict[SpectralMethod, float],
    int,
    int,
    int,
    float,
    float,
    float,
]:
    """Evaluate exact and approximate NMA spectra on one network realization."""
    gen, univ, nma, n_ent = _generate_challenging_network(
        regime=regime,
        box_l=box_l,
        seed=seed,
        pdms=pdms,
    )
    vol = box_l**3
    scale = pdms.gamma_conv / vol
    n_atoms = univ.get_nr_of_atoms()
    n_zero = int(nma.num_soluble_clusters)

    fb = gen.get_mehp_simulator()
    fb.run()
    w_act = float(fb.get_active_weight_fraction())
    w_dang = float(fb.get_dangling_weight_fraction())
    w_sol = float(fb.get_soluble_weight_fraction())

    gp_map: Dict[SpectralMethod, np.ndarray] = {}
    gpp_map: Dict[SpectralMethod, np.ndarray] = {}
    time_map: Dict[SpectralMethod, float] = {}

    for method in methods:
        t0 = time.perf_counter()
        if method == SpectralMethod.EXACT:
            spec = nma.solve_exact(
                solver=pt.NormalModeSolver.DENSE_LAPACK,
                compute_eigenvectors=False,
                decompose_components=True,
            )
        else:
            spec = nma.solve_approx(_build_approx_config(method))
        gp, gpp = spec.moduli(freqs)
        time_map[method] = time.perf_counter() - t0
        gp_map[method] = gp * scale
        gpp_map[method] = gpp * scale

    return (
        regime,
        rep_idx,
        gp_map,
        gpp_map,
        time_map,
        n_atoms,
        n_ent,
        n_zero,
        w_act,
        w_dang,
        w_sol,
    )


def _evaluate_dos_and_localization(
    regime: RegimeKind,
    seed: int,
    box_l: float,
    lambda_grid: np.ndarray,
    bin_edges: np.ndarray,
    kpm_moments: Sequence[int],
    pdms: PdmsParams,
) -> Dict[str, np.ndarray]:
    """Compute KPM DOS vs. exact eigenvalue histogram and eigenvector localization."""
    _, univ, nma, n_ent = _generate_challenging_network(
        regime=regime,
        box_l=box_l,
        seed=seed,
        pdms=pdms,
    )
    spec = nma.solve_exact(
        solver=pt.NormalModeSolver.DENSE_LAPACK,
        compute_eigenvectors=True,
        decompose_components=True,
    )
    evals = np.asarray(spec.eigenvalues, dtype=float)
    evecs = spec.eigenvectors
    assert evecs is not None

    nz_mask = evals > 1e-9
    evals_nz = evals[nz_mask]
    evecs_nz = evecs[:, nz_mask]

    # Exact histogram normalized to modes per unit lambda (matches density_of_states)
    counts, _ = np.histogram(evals_nz, bins=bin_edges)
    exact_dos_hist = counts.astype(float) / np.diff(bin_edges)

    dos_kpm = np.zeros((len(kpm_moments), len(lambda_grid)), dtype=float)
    for idx, num_m in enumerate(kpm_moments):
        dos_kpm[idx] = nma.density_of_states(
            lambda_grid,
            pt.KPMConfig(num_moments=int(num_m), num_random_vectors=20),
        )

    # Eigenvector localization across all nonzero modes
    v2 = evecs_nz**2
    col_norm2 = np.sum(v2, axis=0)
    ipr = np.sum(v2**2, axis=0) / np.maximum(col_norm2**2, 1e-30)

    atom_types = np.asarray(univ.get_atom_types(), dtype=int)
    w_type1 = np.sum(v2[atom_types == 1], axis=0) / np.maximum(col_norm2, 1e-30)
    w_type2 = np.sum(v2[atom_types == 2], axis=0) / np.maximum(col_norm2, 1e-30)
    if np.any(atom_types == 4):
        w_type4 = np.sum(v2[atom_types == 4], axis=0) / np.maximum(col_norm2, 1e-30)
    else:
        w_type4 = np.zeros_like(w_type1)

    # Largest connected component vs. soluble clusters
    _, comp_labels = connected_components(nma.matrix, directed=False)
    comp_sizes = np.bincount(comp_labels)
    giant_comp_id = int(np.argmax(comp_sizes))
    in_giant = comp_labels == giant_comp_id
    w_giant = np.sum(v2[in_giant], axis=0) / np.maximum(col_norm2, 1e-30)
    w_soluble = 1.0 - w_giant

    return {
        "evals_nz": evals_nz,
        "exact_dos_hist": exact_dos_hist,
        "dos_kpm": dos_kpm,
        "ipr": ipr,
        "w_type1": w_type1,
        "w_type2": w_type2,
        "w_type4": w_type4,
        "w_giant": w_giant,
        "w_soluble": w_soluble,
        "n_atoms": np.array(univ.get_nr_of_atoms(), dtype=int),
        "n_ent": np.array(n_ent, dtype=int),
        "n_zero": np.array(int(nma.num_soluble_clusters), dtype=int),
        "lambda_max": np.array(float(np.max(evals)), dtype=float),
    }


def _evaluate_boundary_task(
    rep_idx: int,
    seed: int,
    box_l: float,
    strand_len: int,
    freqs: np.ndarray,
    pdms: PdmsParams,
) -> Tuple[
    int,
    Dict[BoundaryCase, np.ndarray],
    Dict[BoundaryCase, np.ndarray],
    Dict[BoundaryCase, float],
    Dict[BoundaryCase, float],
    Dict[BoundaryCase, int],
    Dict[BoundaryCase, float],
    int,
]:
    """Compare NMA and MEHP on uncrosslinked melt (p=0) vs. crosslinked network (p=0.95)."""
    vol = box_l**3
    scale = pdms.gamma_conv / vol
    n_strands = max(1, int(round(pdms.bead_density * vol / (strand_len + 0.5))))
    n_xlinks = max(1, int(round(0.5 * n_strands)))

    gp_nma: Dict[BoundaryCase, np.ndarray] = {}
    gpp_nma: Dict[BoundaryCase, np.ndarray] = {}
    g_mehp_gen: Dict[BoundaryCase, float] = {}
    g_mehp_static: Dict[BoundaryCase, float] = {}
    n_zero_map: Dict[BoundaryCase, int] = {}
    min_nz_map: Dict[BoundaryCase, float] = {}
    n_ent_sampled = 0

    for p_val, no_ent_case, ent_case in (
        (0.0, BoundaryCase.MELT_NO_ENT, BoundaryCase.MELT_STATIC_ENT),
        (0.95, BoundaryCase.NET_NO_ENT, BoundaryCase.NET_STATIC_ENT),
    ):
        gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
        gen.set_seed(int(seed))
        gen.config_nr_of_mc_steps(0)
        gen.set_bead_distance(pdms.b)
        gen.set_mean_squared_bead_distance(pdms.b2)
        gen.add_strands(n_strands, [int(strand_len)] * n_strands)
        gen.add_crosslinkers(n_xlinks, 4)
        if p_val > 0.0:
            gen.link_strands_to_conversion(p_val)

        # Without entanglements
        u_no = gen.get_universe()
        nma_no = pt.NormalModeAnalyzer.from_universe(u_no)
        spec_no = nma_no.solve_exact(
            solver=pt.NormalModeSolver.DENSE_LAPACK,
            compute_eigenvectors=False,
            decompose_components=True,
        )
        gp_no, gpp_no = spec_no.moduli(freqs)
        gp_nma[no_ent_case] = gp_no * scale
        gpp_nma[no_ent_case] = gpp_no * scale
        n_zero_map[no_ent_case] = int(nma_no.num_soluble_clusters)
        ev_no = np.asarray(spec_no.eigenvalues, dtype=float)
        nz_no = ev_no[ev_no > 1e-9]
        min_nz_map[no_ent_case] = float(nz_no[0]) if len(nz_no) > 0 else float("nan")

        fb_no = gen.get_mehp_simulator(entanglement_mode=pt.EntanglementMode.NONE)
        fb_no.run()
        g_no = float(scale * np.sum(fb_no.get_gamma_factors(pdms.r02)))
        g_mehp_gen[no_ent_case] = g_no
        g_mehp_static[no_ent_case] = g_no

        # With entanglements sampled at melt entanglement density
        es = gen.sample_entanglements(
            nr_of_entanglements_to_sample=int(pdms.ent_density * vol),
            upper_cutoff=pdms.cutoff,
            seed=f"bnd_{seed}_{p_val:.2f}",
            filtering_mode=pt.EntanglementFilteringMode.NONE,
        )
        pairs = es.to_atom_pair_entanglements(
            gen.get_crosslinker_universe()
        ).pairs_of_atoms
        n_ent_sampled = len(pairs)

        u_ent = gen.get_universe()
        if pairs:
            u_ent.add_bonds([pr[0] for pr in pairs], [pr[1] for pr in pairs])

        nma_ent = pt.NormalModeAnalyzer.from_universe(u_ent)
        spec_ent = nma_ent.solve_exact(
            solver=pt.NormalModeSolver.DENSE_LAPACK,
            compute_eigenvectors=False,
            decompose_components=True,
        )
        gp_ent, gpp_ent = spec_ent.moduli(freqs)
        gp_nma[ent_case] = gp_ent * scale
        gpp_nma[ent_case] = gpp_ent * scale
        n_zero_map[ent_case] = int(nma_ent.num_soluble_clusters)
        ev_ent = np.asarray(spec_ent.eigenvalues, dtype=float)
        nz_ent = ev_ent[ev_ent > 1e-9]
        min_nz_map[ent_case] = float(nz_ent[0]) if len(nz_ent) > 0 else float("nan")

        # MEHP via CrosslinkerUniverse (knows chemical crosslinks vs. melt)
        fb_gen_ent = gen.get_mehp_simulator(
            es, entanglement_mode=pt.EntanglementMode.MERGED_LINKS
        )
        fb_gen_ent.run()
        g_mehp_gen[ent_case] = float(
            scale * np.sum(fb_gen_ent.get_gamma_factors(pdms.r02))
        )

        # MEHP from universe where entanglements were turned into permanent bonds
        fb_static = pt.MEHPSimulator(u_ent, entanglement_mode=pt.EntanglementMode.NONE)
        fb_static.run()
        g_mehp_static[ent_case] = float(
            scale * np.sum(fb_static.get_gamma_factors(pdms.r02))
        )

    return (
        rep_idx,
        gp_nma,
        gpp_nma,
        g_mehp_gen,
        g_mehp_static,
        n_zero_map,
        min_nz_map,
        n_ent_sampled,
    )


def _evaluate_pbc_coordinate_blindness(
    box_l: float,
    strand_len: int,
    seed: int,
    pdms: PdmsParams,
) -> Dict[str, np.ndarray]:
    """Demonstrate NMA coordinate/PBC invariance vs. MEHP spatial sensitivity."""
    vol = box_l**3
    scale = pdms.gamma_conv / vol
    n_strands = max(1, int(round(pdms.bead_density * vol / (strand_len + 0.5))))
    n_xlinks = max(1, int(round(0.5 * n_strands)))

    gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(int(seed))
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(pdms.b)
    gen.set_mean_squared_bead_distance(pdms.b2)
    gen.add_strands(n_strands, [int(strand_len)] * n_strands)
    gen.add_crosslinkers(n_xlinks, 4)
    gen.link_strands_to_conversion(0.95)
    u_orig = gen.get_universe()

    # Identical connectivity graph collapsed to local coordinates (zero PBC crossings)
    atoms = u_orig.get_atoms()
    u_compact = pt.Universe(box_l, box_l, box_l)
    u_compact.add_atoms(
        [a.get_id() for a in atoms],
        [a.get_type() for a in atoms],
        [0.05 * a.get_x() for a in atoms],
        [0.05 * a.get_y() for a in atoms],
        [0.05 * a.get_z() for a in atoms],
        [0] * len(atoms),
        [0] * len(atoms),
        [0] * len(atoms),
    )
    bonds = u_orig.get_bonds()
    u_compact.add_bonds(bonds["bond_from"], bonds["bond_to"])

    ev_orig = (
        pt.NormalModeAnalyzer.from_universe(u_orig)
        .solve_exact(solver=pt.NormalModeSolver.DENSE_LAPACK)
        .eigenvalues
    )
    ev_comp = (
        pt.NormalModeAnalyzer.from_universe(u_compact)
        .solve_exact(solver=pt.NormalModeSolver.DENSE_LAPACK)
        .eigenvalues
    )
    max_eval_diff = float(np.max(np.abs(ev_orig - ev_comp)))

    fb_orig = pt.MEHPSimulator(u_orig)
    fb_orig.run()
    fb_comp = pt.MEHPSimulator(u_compact)
    fb_comp.run()

    return {
        "pbc_nma_max_eval_diff": np.array(max_eval_diff, dtype=float),
        "pbc_mehp_g_orig": np.array(
            float(scale * np.sum(fb_orig.get_gamma_factors(pdms.r02))), dtype=float
        ),
        "pbc_mehp_g_compact": np.array(
            float(scale * np.sum(fb_comp.get_gamma_factors(pdms.r02))), dtype=float
        ),
    }


def run_kpm_limitations_study(
    max_workers: Optional[int] = None,
    n_repeats: int = 8,
    n_boundary_repeats: int = 6,
    box_l: float = 14.0,
    box_l_boundary: float = 13.0,
    base_seed: int = 4200,
) -> Dict[str, Any]:
    """Run Study 3: KPM low-frequency failure, DOS smearing, localization & NMA/MEHP boundary."""
    pdms = _load_pdms_params()
    freqs = np.logspace(-5, 1, 75)
    n_freqs = len(freqs)
    methods: List[SpectralMethod] = [
        SpectralMethod.EXACT,
        SpectralMethod.KPM_256,
        SpectralMethod.KPM_1024,
        SpectralMethod.KPM_4096,
        SpectralMethod.SLQ,
        SpectralMethod.HYBRID_50,
        SpectralMethod.HYBRID_100,
    ]
    regimes: List[RegimeKind] = [
        RegimeKind.NEAR_PERCOLATION,
        RegimeKind.DANGLING_BRUSH,
    ]
    boundary_cases: List[BoundaryCase] = [
        BoundaryCase.MELT_NO_ENT,
        BoundaryCase.MELT_STATIC_ENT,
        BoundaryCase.NET_NO_ENT,
        BoundaryCase.NET_STATIC_ENT,
    ]

    # Allocate Part 1 arrays: shape (n_methods, n_repeats, n_freqs) per regime
    gp_store: Dict[RegimeKind, np.ndarray] = {
        reg: np.zeros((len(methods), n_repeats, n_freqs), dtype=float)
        for reg in regimes
    }
    gpp_store: Dict[RegimeKind, np.ndarray] = {
        reg: np.zeros((len(methods), n_repeats, n_freqs), dtype=float)
        for reg in regimes
    }
    t_store: Dict[RegimeKind, np.ndarray] = {
        reg: np.zeros((len(methods), n_repeats), dtype=float) for reg in regimes
    }
    meta_store: Dict[RegimeKind, Dict[str, np.ndarray]] = {
        reg: {
            "n_atoms": np.zeros(n_repeats, dtype=int),
            "n_ent": np.zeros(n_repeats, dtype=int),
            "n_zero": np.zeros(n_repeats, dtype=int),
            "w_act": np.zeros(n_repeats, dtype=float),
            "w_dang": np.zeros(n_repeats, dtype=float),
            "w_sol": np.zeros(n_repeats, dtype=float),
        }
        for reg in regimes
    }

    # Allocate Part 3 arrays
    bnd_gp = np.zeros((len(boundary_cases), n_boundary_repeats, n_freqs), dtype=float)
    bnd_gpp = np.zeros((len(boundary_cases), n_boundary_repeats, n_freqs), dtype=float)
    bnd_g_mehp_gen = np.zeros((len(boundary_cases), n_boundary_repeats), dtype=float)
    bnd_g_mehp_static = np.zeros((len(boundary_cases), n_boundary_repeats), dtype=float)
    bnd_n_zero = np.zeros((len(boundary_cases), n_boundary_repeats), dtype=int)
    bnd_min_nz = np.zeros((len(boundary_cases), n_boundary_repeats), dtype=float)
    bnd_n_ent = np.zeros(n_boundary_repeats, dtype=int)

    regime_tasks = []
    for reg_idx, reg in enumerate(regimes):
        for rep in range(n_repeats):
            seed = base_seed + 100000 * reg_idx + 313 * rep
            regime_tasks.append((reg, rep, seed, box_l, freqs, methods, pdms))

    boundary_tasks = []
    for rep in range(n_boundary_repeats):
        seed = base_seed + 300000 + 313 * rep
        boundary_tasks.append((rep, seed, box_l_boundary, 100, freqs, pdms))

    total_tasks = len(regime_tasks) + len(boundary_tasks)
    workers = min(resolve_max_workers(max_workers), total_tasks)
    configure_process_threading(workers)

    def _store_regime(res: Any) -> None:
        (
            reg,
            rep,
            gp_map,
            gpp_map,
            time_map,
            n_at,
            n_ent,
            n_zero,
            w_act,
            w_dang,
            w_sol,
        ) = res
        for m_idx, method in enumerate(methods):
            gp_store[reg][m_idx, rep] = gp_map[method]
            gpp_store[reg][m_idx, rep] = gpp_map[method]
            t_store[reg][m_idx, rep] = time_map[method]
        meta_store[reg]["n_atoms"][rep] = n_at
        meta_store[reg]["n_ent"][rep] = n_ent
        meta_store[reg]["n_zero"][rep] = n_zero
        meta_store[reg]["w_act"][rep] = w_act
        meta_store[reg]["w_dang"][rep] = w_dang
        meta_store[reg]["w_sol"][rep] = w_sol

    def _store_boundary(res: Any) -> None:
        (
            rep,
            gp_nma,
            gpp_nma,
            g_gen,
            g_stat,
            n_zero_map,
            min_nz_map,
            n_ent_val,
        ) = res
        for c_idx, bcase in enumerate(boundary_cases):
            bnd_gp[c_idx, rep] = gp_nma[bcase]
            bnd_gpp[c_idx, rep] = gpp_nma[bcase]
            bnd_g_mehp_gen[c_idx, rep] = g_gen[bcase]
            bnd_g_mehp_static[c_idx, rep] = g_stat[bcase]
            bnd_n_zero[c_idx, rep] = n_zero_map[bcase]
            bnd_min_nz[c_idx, rep] = min_nz_map[bcase]
        bnd_n_ent[rep] = n_ent_val

    if workers > 1:
        print(
            f"Executing {total_tasks} KPM limitations & boundary tasks across {workers} workers..."
        )
        completed = 0
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            reg_futs = [
                executor.submit(_evaluate_regime_task, *t) for t in regime_tasks
            ]
            bnd_futs = [
                executor.submit(_evaluate_boundary_task, *t) for t in boundary_tasks
            ]
            for fut in as_completed(reg_futs):
                _store_regime(fut.result())
                completed += 1
                print(
                    f"  [KPM Limitations] Completed {completed}/{total_tasks} tasks..."
                )
            for fut in as_completed(bnd_futs):
                _store_boundary(fut.result())
                completed += 1
                print(
                    f"  [KPM Limitations] Completed {completed}/{total_tasks} tasks..."
                )
    else:
        completed = 0
        for t in regime_tasks:
            _store_regime(_evaluate_regime_task(*t))
            completed += 1
            print(f"  [KPM Limitations] Completed {completed}/{total_tasks} tasks...")
        for t in boundary_tasks:
            _store_boundary(_evaluate_boundary_task(*t))
            completed += 1
            print(f"  [KPM Limitations] Completed {completed}/{total_tasks} tasks...")

    # Part 2: DOS smearing & eigenvector localization on representative percolating realizations
    lambda_grid = np.logspace(-6, 1.0, 200)
    bin_edges = np.logspace(-6, 1.0, 45)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
    kpm_moments = np.array([256, 1024, 4096], dtype=int)

    print("Evaluating DOS smearing & eigenvector localization (with eigenvectors)...")
    dos_np = _evaluate_dos_and_localization(
        RegimeKind.NEAR_PERCOLATION,
        seed=base_seed,
        box_l=box_l,
        lambda_grid=lambda_grid,
        bin_edges=bin_edges,
        kpm_moments=kpm_moments,
        pdms=pdms,
    )
    dos_db = _evaluate_dos_and_localization(
        RegimeKind.DANGLING_BRUSH,
        seed=base_seed + 100000,
        box_l=box_l,
        lambda_grid=lambda_grid,
        bin_edges=bin_edges,
        kpm_moments=kpm_moments,
        pdms=pdms,
    )
    pbc_info = _evaluate_pbc_coordinate_blindness(
        box_l=box_l_boundary, strand_len=60, seed=base_seed, pdms=pdms
    )

    out: Dict[str, Any] = {
        "freqs": freqs,
        "box_l": np.array(box_l, dtype=float),
        "box_l_boundary": np.array(box_l_boundary, dtype=float),
        "n_repeats": np.array(n_repeats, dtype=int),
        "n_boundary_repeats": np.array(n_boundary_repeats, dtype=int),
        "method_names": np.array([m.value for m in methods]),
        "boundary_case_names": np.array([c.value for c in boundary_cases]),
        "lambda_grid": lambda_grid,
        "dos_bin_edges": bin_edges,
        "dos_bin_centers": bin_centers,
        "dos_kpm_moments": kpm_moments,
    }

    for reg in regimes:
        prefix = "np" if reg == RegimeKind.NEAR_PERCOLATION else "db"
        gp_arr = gp_store[reg]
        gpp_arr = gpp_store[reg]
        stat_gp = compute_ensemble_statistics(gp_arr, axis=1)
        stat_gpp = compute_ensemble_statistics(gpp_arr, axis=1)

        # Compute per-realization relative errors vs Exact (method index 0), then average
        ex_gp = gp_arr[0:1, :, :]
        ex_gpp = gpp_arr[0:1, :, :]
        rel_err_gp = np.abs(gp_arr - ex_gp) / np.maximum(ex_gp, 1e-12) * 100.0
        rel_err_gpp = np.abs(gpp_arr - ex_gpp) / np.maximum(ex_gpp, 1e-12) * 100.0
        stat_err_gp = compute_ensemble_statistics(rel_err_gp, axis=1)
        stat_err_gpp = compute_ensemble_statistics(rel_err_gpp, axis=1)

        out[f"{prefix}_gp_mean"] = stat_gp["mean"]
        out[f"{prefix}_gp_ci"] = stat_gp["ci"]
        out[f"{prefix}_gpp_mean"] = stat_gpp["mean"]
        out[f"{prefix}_gpp_ci"] = stat_gpp["ci"]
        out[f"{prefix}_err_gp_mean"] = stat_err_gp["mean"]
        out[f"{prefix}_err_gp_ci"] = stat_err_gp["ci"]
        out[f"{prefix}_err_gpp_mean"] = stat_err_gpp["mean"]
        out[f"{prefix}_err_gpp_ci"] = stat_err_gpp["ci"]
        out[f"{prefix}_t_mean"] = np.mean(t_store[reg], axis=1)
        for m_key, m_arr in meta_store[reg].items():
            out[f"{prefix}_{m_key}"] = m_arr

    for prefix, dos_res in (("np_dos", dos_np), ("db_dos", dos_db)):
        for k, v in dos_res.items():
            out[f"{prefix}_{k}"] = v

    stat_bnd_gp = compute_ensemble_statistics(bnd_gp, axis=1)
    stat_bnd_gpp = compute_ensemble_statistics(bnd_gpp, axis=1)
    stat_bnd_g_gen = compute_ensemble_statistics(bnd_g_mehp_gen, axis=1)
    stat_bnd_g_stat = compute_ensemble_statistics(bnd_g_mehp_static, axis=1)

    out["bnd_gp_mean"] = stat_bnd_gp["mean"]
    out["bnd_gp_ci"] = stat_bnd_gp["ci"]
    out["bnd_gpp_mean"] = stat_bnd_gpp["mean"]
    out["bnd_gpp_ci"] = stat_bnd_gpp["ci"]
    out["bnd_g_mehp_gen_mean"] = stat_bnd_g_gen["mean"]
    out["bnd_g_mehp_gen_ci"] = stat_bnd_g_gen["ci"]
    out["bnd_g_mehp_static_mean"] = stat_bnd_g_stat["mean"]
    out["bnd_g_mehp_static_ci"] = stat_bnd_g_stat["ci"]
    out["bnd_n_zero_mean"] = np.mean(bnd_n_zero, axis=1)
    out["bnd_min_nz_mean"] = np.mean(bnd_min_nz, axis=1)
    out["bnd_n_ent_mean"] = np.mean(bnd_n_ent)
    out.update(pbc_info)

    return out
