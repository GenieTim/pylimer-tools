"""
Multi-seed simulation runner for Study 6:
Multi-Chemistry Crossover & Coarse-Graining Invariance (GAUSSIAN vs. KUHN vs. KG_LJ).

This module is prefixed with an underscore so Sphinx-Gallery ignores it when
building documentation.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from enum import Enum
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple
import warnings

import numpy as np

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
from pylimer_tools.calc import miller_macosko_theory as mmt  # noqa: E402
from pylimer_tools.io.bead_spring_parameter_provider import (  # noqa: E402
    Parameters,
    ParameterType,
    _get_relevant_everaers_row,
    get_parameters_for_polymer,
)
import pylimer_tools_cpp as pt  # noqa: E402

_GAS_CONSTANT_J_PER_MOL_K: float = 8.31446261815324
_AVOGADRO_NUMBER: float = 6.02214076e23
_CM3_TO_NM3: float = 1.0e21


class PolymerChemistry(Enum):
    """Contrasting polymer chemistries from the Everaers et al. dataset."""

    PE = "PE"
    CIS_PBD = "cis-PBd"
    PDMS = "PDMS*"
    A_PMMA = "a-PMMA"
    A_PS = "a-PS*"


class CoarseGrainModel(Enum):
    """Supported bead-spring coarse-graining parameterizations."""

    GAUSSIAN = ParameterType.GAUSSIAN
    KUHN = ParameterType.KUHN
    KG_LJ = ParameterType.KG_LJ


@dataclass(frozen=True)
class PolymerSpec:
    """Metadata and display properties for a selected polymer chemistry."""

    chemistry: PolymerChemistry
    display_name: str
    category: str
    color: str
    marker: str


POLYMER_SPECS: Tuple[PolymerSpec, ...] = (
    PolymerSpec(
        chemistry=PolymerChemistry.PE,
        display_name="PE (298 K)",
        category="Tightly entangled",
        color="#d73027",
        marker="o",
    ),
    PolymerSpec(
        chemistry=PolymerChemistry.CIS_PBD,
        display_name="cis-PBd (298 K)",
        category="Tightly entangled",
        color="#fc8d59",
        marker="s",
    ),
    PolymerSpec(
        chemistry=PolymerChemistry.PDMS,
        display_name="PDMS (298 K)",
        category="Loosely entangled",
        color="#2166ac",
        marker="^",
    ),
    PolymerSpec(
        chemistry=PolymerChemistry.A_PMMA,
        display_name="a-PMMA (413 K)",
        category="Stiff / glassy-melt",
        color="#1b7837",
        marker="D",
    ),
    PolymerSpec(
        chemistry=PolymerChemistry.A_PS,
        display_name="a-PS (413 K)",
        category="Stiff / glassy-melt",
        color="#762a83",
        marker="v",
    ),
)

COARSE_GRAIN_POLYMERS: Tuple[PolymerChemistry, ...] = (
    PolymerChemistry.CIS_PBD,
    PolymerChemistry.PDMS,
)

COARSE_GRAIN_MODELS: Tuple[CoarseGrainModel, ...] = (
    CoarseGrainModel.GAUSSIAN,
    CoarseGrainModel.KUHN,
    CoarseGrainModel.KG_LJ,
)


@dataclass(frozen=True)
class ModelConstants:
    """Dimensionally consistent simulation constants in a model's native distance units."""

    polymer_name: str
    parameter_type: ParameterType
    distance_unit_nm: float
    mw_g_per_mol: float
    mk_g_per_mol: float
    lk_nm: float
    nek: float
    rho_g_per_cm3: float
    temperature_k: float
    ge_mpa: float
    me_g_per_mol: float
    b_sim: float
    b2_sim: float
    r02_sim: float
    bead_density_sim: float
    ent_density_sim: float
    cutoff_sim: float
    gamma_conv_mpa: float


def load_model_constants(
    chemistry: PolymerChemistry,
    model: CoarseGrainModel,
) -> ModelConstants:
    """Load bead-spring parameters and extract native-unit simulation constants.

    Note that for ``ParameterType.KG_LJ``, ``distance_units`` is a Pint Quantity
    with magnitude ``sigma`` (in nm) and unit ``nanometer``. Converting lengths
    and densities by dividing/multiplying directly with ``distance_units`` guarantees
    exact dimensional consistency across ``GAUSSIAN``, ``KUHN``, and ``KG_LJ``.
    """
    polymer_name = chemistry.value
    parameter_type = model.value
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        params: Parameters = get_parameters_for_polymer(polymer_name, parameter_type)
    row = _get_relevant_everaers_row(polymer_name)

    du = params.get("distance_units")
    du_nm = float(du.to("nanometer").magnitude)
    mw = float(params.get("Mw").to("g/mol").magnitude)
    rho_g_cm3 = float(params.get("rho").to("g/cm^3").magnitude)
    temperature_k = float(params.get("T").to("kelvin").magnitude)
    ge_mpa = float(params.get("Ge").to("MPa").magnitude)
    # 1 MPa = 1 J/cm^3, so rho [g/cm^3] * R [J/(mol K)] * T [K] / Ge [MPa] -> g/mol
    me_g_per_mol = rho_g_cm3 * _GAS_CONSTANT_J_PER_MOL_K * temperature_k / ge_mpa

    b_sim = float((params.get("<b>") / du).to_reduced_units().magnitude)
    b2_sim = float((params.get("<b^2>") / (du**2)).to_reduced_units().magnitude)
    r02_sim = float((params.get("R02") / (du**2)).to_reduced_units().magnitude)
    bead_density_sim = float(params.get_bead_density())
    ent_density_sim = float(
        ((params.get("Ge") / (params.get("kb") * params.get("T"))) * (du**3))
        .to_reduced_units()
        .magnitude
    )
    cutoff_sim = float(params.get_sampling_cutoff())
    gamma_conv_mpa = float(params.get_gamma_conversion_factor().to("MPa").magnitude)

    return ModelConstants(
        polymer_name=polymer_name,
        parameter_type=parameter_type,
        distance_unit_nm=du_nm,
        mw_g_per_mol=mw,
        mk_g_per_mol=float(row["M_k"]),
        lk_nm=float(row["l_K"]) * 0.1,
        nek=float(row["N_ek"]),
        rho_g_per_cm3=rho_g_cm3,
        temperature_k=temperature_k,
        ge_mpa=ge_mpa,
        me_g_per_mol=me_g_per_mol,
        b_sim=b_sim,
        b2_sim=b2_sim,
        r02_sim=r02_sim,
        bead_density_sim=bead_density_sim,
        ent_density_sim=ent_density_sim,
        cutoff_sim=cutoff_sim,
        gamma_conv_mpa=gamma_conv_mpa,
    )


def _make_linear_solver_config() -> pt.MEHPLinearSolverConfig:
    """Construct the standard linear MEHP solver configuration."""
    cfg = pt.MEHPLinearSolverConfig()
    cfg.simplification_mode = pt.StructureSimplificationMode.INACTIVE_THEN_X2F
    return cfg


def _simulate_single_network(
    constants: ModelConstants,
    n_strands: int,
    n_beads: int,
    vol_sim: float,
    conversion: float,
    functionality: int,
    seed: int,
    ent_seed: str,
) -> Tuple[float, float, float, int, int]:
    """Build one MC network and evaluate phantom & entangled MEHP shear moduli."""
    n_xlinks = n_strands // (functionality // 2)
    box_l = vol_sim ** (1.0 / 3.0)

    gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(seed)
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(constants.b_sim)
    gen.set_mean_squared_bead_distance(constants.b2_sim)
    gen.add_crosslinkers(n_xlinks, functionality)
    gen.add_strands(n_strands, [n_beads] * n_strands, strand_atom_type=1)
    gen.link_strands_to_conversion(conversion)

    solver_cfg = _make_linear_solver_config()
    fb_ph = gen.get_mehp_simulator(
        entanglement_mode=pt.EntanglementMode.NONE,
        crosslinker_type=2,
    )
    fb_ph.run(solver_cfg)
    g_ph = float(
        constants.gamma_conv_mpa
        * np.sum(fb_ph.get_gamma_factors(constants.r02_sim))
        / vol_sim
    )
    w_act = float(fb_ph.get_active_weight_fraction())

    n_ent_target = int(round(constants.ent_density_sim * vol_sim))
    es = gen.sample_entanglements(
        n_ent_target,
        constants.cutoff_sim,
        0.0,
        0,
        3.0,
        ent_seed,
        pt.EntanglementFilteringMode.OMIT_INACTIVE,
    )
    n_ent_active = int(es.size())

    fb_en = gen.get_mehp_simulator(
        es,
        entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
        crosslinker_type=2,
    )
    fb_en.run(solver_cfg)
    g_en = float(
        constants.gamma_conv_mpa
        * np.sum(fb_en.get_gamma_factors(constants.r02_sim))
        / vol_sim
    )
    return g_ph, g_en, w_act, n_ent_active, n_ent_target


def _worker_part1_point(
    task: Tuple[int, int, int, str, float, int, float, int],
) -> Tuple[int, int, int, float, float, float, int, int]:
    """Worker for one (polymer, target Mc/Me, seed) point in Part 1."""
    (
        poly_idx,
        ratio_idx,
        seed_idx,
        poly_val,
        target_ratio,
        n_strands,
        conversion,
        seed,
    ) = task
    chemistry = PolymerChemistry(poly_val)
    constants = load_model_constants(chemistry, CoarseGrainModel.GAUSSIAN)
    target_mc = target_ratio * constants.me_g_per_mol
    n_beads = max(2, int(round(target_mc / constants.mw_g_per_mol)))
    n_xlinks = n_strands // 2
    vol_sim = (n_strands * n_beads + n_xlinks) / constants.bead_density_sim

    g_ph, g_en, w_act, n_ent_active, n_ent_target = _simulate_single_network(
        constants=constants,
        n_strands=n_strands,
        n_beads=n_beads,
        vol_sim=vol_sim,
        conversion=conversion,
        functionality=4,
        seed=seed,
        ent_seed=f"p1_{poly_idx}_{ratio_idx}_{seed_idx}",
    )
    return (
        poly_idx,
        ratio_idx,
        seed_idx,
        g_ph,
        g_en,
        w_act,
        n_ent_active,
        n_ent_target,
    )


def _worker_part2_point(
    task: Tuple[int, int, int, int, str, str, float, int, float, int],
) -> Tuple[int, int, int, int, float, float, float, int, int]:
    """Worker for one (polymer, model, target Mc/Mk, seed) point in Part 2."""
    (
        poly_idx,
        model_idx,
        ratio_idx,
        seed_idx,
        poly_val,
        model_name,
        mc_over_mk,
        n_strands,
        conversion,
        seed,
    ) = task
    chemistry = PolymerChemistry(poly_val)
    model = CoarseGrainModel[model_name]
    constants = load_model_constants(chemistry, model)

    target_mc = mc_over_mk * constants.mk_g_per_mol
    n_beads = max(2, int(round(target_mc / constants.mw_g_per_mol)))

    # Match the physical box volume in nm^3 identically across GAUSSIAN, KUHN, and KG_LJ
    rho_g_per_nm3 = constants.rho_g_per_cm3 / _CM3_TO_NM3
    vol_nm3 = (n_strands * target_mc) / (rho_g_per_nm3 * _AVOGADRO_NUMBER)
    vol_sim = vol_nm3 / (constants.distance_unit_nm**3)

    g_ph, g_en, w_act, n_ent_active, n_ent_target = _simulate_single_network(
        constants=constants,
        n_strands=n_strands,
        n_beads=n_beads,
        vol_sim=vol_sim,
        conversion=conversion,
        functionality=4,
        seed=seed,
        ent_seed=f"p2_{poly_idx}_{model_idx}_{ratio_idx}_{seed_idx}",
    )
    return (
        poly_idx,
        model_idx,
        ratio_idx,
        seed_idx,
        g_ph,
        g_en,
        w_act,
        n_ent_active,
        n_ent_target,
    )


def _compute_universal_mmt_reference(
    mc_over_me_grid: np.ndarray,
    conversion: float,
    functionality: int = 4,
    stoichiometry: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """Compute chemistry-independent normalized MMT crossover curves G/Ge vs. Mc/Me.

    Since ``nu * k_B * T = (rho / M_c) * R * T = G_e / (M_c / M_e)``, the
    dimensionless ratios ``G_ph / G_e`` and ``G_tot / G_e`` depend only on
    ``(r, p, f)`` and ``M_c / M_e``.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        params = get_parameters_for_polymer(
            PolymerChemistry.PDMS.value, ParameterType.GAUSSIAN
        )
    ureg = params.get_unit_registry()
    temperature = params.get("T")
    g_e = params.get("Ge")
    kb = params.get("kb")
    ge_over_kbt_nm3 = float((g_e / (kb * temperature)).to("nm^-3").magnitude)

    g_ph_norm = np.zeros_like(mc_over_me_grid, dtype=float)
    g_tot_norm = np.zeros_like(mc_over_me_grid, dtype=float)
    g_e_mpa = float(g_e.to("MPa").magnitude)

    for idx, ratio in enumerate(mc_over_me_grid):
        nu_nm3 = ge_over_kbt_nm3 / float(ratio)
        g_ph_q, g_en_q, _, _ = mmt.compute_modulus_decomposition(
            ureg=ureg,
            r=stoichiometry,
            p=conversion,
            f=functionality,
            nu=nu_nm3 * (ureg("nm") ** -3),
            temperature=temperature,
            g_e_1=g_e,
        )
        g_ph_val = float(g_ph_q.to("MPa").magnitude)
        g_en_val = float(g_en_q.to("MPa").magnitude)
        g_ph_norm[idx] = g_ph_val / g_e_mpa
        g_tot_norm[idx] = (g_ph_val + g_en_val) / g_e_mpa

    # Extract asymptotic phantom prefactor and trapping factor T_e
    phantom_prefactor = float(g_ph_norm[0] * mc_over_me_grid[0])
    trapping_factor = float(g_tot_norm[-1] - g_ph_norm[-1])
    return g_ph_norm, g_tot_norm, phantom_prefactor, trapping_factor


def run_multichemistry_study(
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Run Study 6: Multi-Chemistry Crossover & Coarse-Graining Invariance.

    Parameters
    ----------
    max_workers : int, optional
        Maximum worker processes for parallel execution. Defaults to
        ``resolve_max_workers(max_workers)``.

    Returns
    -------
    dict
        Dictionary of NumPy arrays suitable for ``np.savez_compressed``.
    """
    workers = resolve_max_workers(max_workers)
    configure_process_threading(workers=workers, verbose=True)

    n_repeats = 6
    n_strands = 1000
    conversion = 0.95
    stoichiometry = 1.0
    functionality = 4

    # ------------------------------------------------------------------
    # Part 1: Multi-Chemistry Crossover Across Mc / Me (GAUSSIAN)
    # ------------------------------------------------------------------
    p1_target_ratios = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0], dtype=float)
    n_polys = len(POLYMER_SPECS)
    n_p1_ratios = len(p1_target_ratios)

    p1_poly_names = np.array([spec.chemistry.value for spec in POLYMER_SPECS])
    p1_display_names = np.array([spec.display_name for spec in POLYMER_SPECS])
    p1_categories = np.array([spec.category for spec in POLYMER_SPECS])
    p1_colors = np.array([spec.color for spec in POLYMER_SPECS])
    p1_markers = np.array([spec.marker for spec in POLYMER_SPECS])

    p1_ge_mpa = np.zeros(n_polys, dtype=float)
    p1_me_g_mol = np.zeros(n_polys, dtype=float)
    p1_mk_g_mol = np.zeros(n_polys, dtype=float)
    p1_mw_g_mol = np.zeros(n_polys, dtype=float)
    p1_nek = np.zeros(n_polys, dtype=float)
    p1_temp_k = np.zeros(n_polys, dtype=float)
    p1_rho_g_cm3 = np.zeros(n_polys, dtype=float)
    p1_lk_nm = np.zeros(n_polys, dtype=float)

    p1_n_beads = np.zeros((n_polys, n_p1_ratios), dtype=int)
    p1_actual_mc = np.zeros((n_polys, n_p1_ratios), dtype=float)
    p1_actual_mc_over_me = np.zeros((n_polys, n_p1_ratios), dtype=float)

    p1_tasks: List[Tuple[int, int, int, str, float, int, float, int]] = []
    for p_idx, spec in enumerate(POLYMER_SPECS):
        c = load_model_constants(spec.chemistry, CoarseGrainModel.GAUSSIAN)
        p1_ge_mpa[p_idx] = c.ge_mpa
        p1_me_g_mol[p_idx] = c.me_g_per_mol
        p1_mk_g_mol[p_idx] = c.mk_g_per_mol
        p1_mw_g_mol[p_idx] = c.mw_g_per_mol
        p1_nek[p_idx] = c.nek
        p1_temp_k[p_idx] = c.temperature_k
        p1_rho_g_cm3[p_idx] = c.rho_g_per_cm3
        p1_lk_nm[p_idx] = c.lk_nm
        for r_idx, target_ratio in enumerate(p1_target_ratios):
            n_b = max(2, int(round(target_ratio * c.me_g_per_mol / c.mw_g_per_mol)))
            p1_n_beads[p_idx, r_idx] = n_b
            p1_actual_mc[p_idx, r_idx] = n_b * c.mw_g_per_mol
            p1_actual_mc_over_me[p_idx, r_idx] = (n_b * c.mw_g_per_mol) / c.me_g_per_mol
            for s_idx in range(n_repeats):
                seed = 61000 + p_idx * 1000 + r_idx * 50 + s_idx
                p1_tasks.append(
                    (
                        p_idx,
                        r_idx,
                        s_idx,
                        spec.chemistry.value,
                        float(target_ratio),
                        n_strands,
                        conversion,
                        seed,
                    )
                )

    p1_g_ph_raw = np.zeros((n_polys, n_p1_ratios, n_repeats), dtype=float)
    p1_g_en_raw = np.zeros((n_polys, n_p1_ratios, n_repeats), dtype=float)
    p1_w_act_raw = np.zeros((n_polys, n_p1_ratios, n_repeats), dtype=float)
    p1_n_ent_active_raw = np.zeros((n_polys, n_p1_ratios, n_repeats), dtype=float)
    p1_n_ent_target = np.zeros((n_polys, n_p1_ratios), dtype=int)

    # ------------------------------------------------------------------
    # Part 2: Coarse-Graining Discretization Invariance & Breakdown
    # ------------------------------------------------------------------
    p2_mc_over_mk = np.array([1.5, 2.0, 3.0, 5.0, 8.0, 15.0, 25.0, 40.0], dtype=float)
    n_p2_polys = len(COARSE_GRAIN_POLYMERS)
    n_p2_models = len(COARSE_GRAIN_MODELS)
    n_p2_ratios = len(p2_mc_over_mk)

    p2_poly_names = np.array([c.value for c in COARSE_GRAIN_POLYMERS])
    p2_model_names = np.array([m.name for m in COARSE_GRAIN_MODELS])
    p2_model_labels = np.array([m.value.value for m in COARSE_GRAIN_MODELS])

    p2_mk_g_mol = np.zeros(n_p2_polys, dtype=float)
    p2_me_g_mol = np.zeros(n_p2_polys, dtype=float)
    p2_ge_mpa = np.zeros(n_p2_polys, dtype=float)
    p2_mw_g_mol = np.zeros((n_p2_polys, n_p2_models), dtype=float)
    p2_du_nm = np.zeros((n_p2_polys, n_p2_models), dtype=float)
    p2_n_beads = np.zeros((n_p2_polys, n_p2_models, n_p2_ratios), dtype=int)
    p2_actual_mc = np.zeros((n_p2_polys, n_p2_models, n_p2_ratios), dtype=float)
    p2_target_mc = np.zeros((n_p2_polys, n_p2_ratios), dtype=float)
    p2_vol_nm3 = np.zeros((n_p2_polys, n_p2_ratios), dtype=float)

    p2_tasks: List[Tuple[int, int, int, int, str, str, float, int, float, int]] = []
    for p_idx, chem in enumerate(COARSE_GRAIN_POLYMERS):
        for m_idx, model in enumerate(COARSE_GRAIN_MODELS):
            c = load_model_constants(chem, model)
            p2_mk_g_mol[p_idx] = c.mk_g_per_mol
            p2_me_g_mol[p_idx] = c.me_g_per_mol
            p2_ge_mpa[p_idx] = c.ge_mpa
            p2_mw_g_mol[p_idx, m_idx] = c.mw_g_per_mol
            p2_du_nm[p_idx, m_idx] = c.distance_unit_nm
            rho_g_per_nm3 = c.rho_g_per_cm3 / _CM3_TO_NM3

            for r_idx, ratio in enumerate(p2_mc_over_mk):
                target_mc = float(ratio) * c.mk_g_per_mol
                p2_target_mc[p_idx, r_idx] = target_mc
                p2_vol_nm3[p_idx, r_idx] = (n_strands * target_mc) / (
                    rho_g_per_nm3 * _AVOGADRO_NUMBER
                )
                n_b = max(2, int(round(target_mc / c.mw_g_per_mol)))
                p2_n_beads[p_idx, m_idx, r_idx] = n_b
                p2_actual_mc[p_idx, m_idx, r_idx] = n_b * c.mw_g_per_mol
                for s_idx in range(n_repeats):
                    seed = 62000 + p_idx * 5000 + m_idx * 1000 + r_idx * 50 + s_idx
                    p2_tasks.append(
                        (
                            p_idx,
                            m_idx,
                            r_idx,
                            s_idx,
                            chem.value,
                            model.name,
                            float(ratio),
                            n_strands,
                            conversion,
                            seed,
                        )
                    )

    p2_g_ph_raw = np.zeros(
        (n_p2_polys, n_p2_models, n_p2_ratios, n_repeats), dtype=float
    )
    p2_g_en_raw = np.zeros(
        (n_p2_polys, n_p2_models, n_p2_ratios, n_repeats), dtype=float
    )
    p2_w_act_raw = np.zeros(
        (n_p2_polys, n_p2_models, n_p2_ratios, n_repeats), dtype=float
    )
    p2_n_ent_active_raw = np.zeros(
        (n_p2_polys, n_p2_models, n_p2_ratios, n_repeats), dtype=float
    )
    p2_n_ent_target = np.zeros((n_p2_polys, n_p2_models, n_p2_ratios), dtype=int)

    print(
        f"Running Study 6 ({len(p1_tasks)} Part-1 tasks + {len(p2_tasks)} Part-2 tasks) "
        f"with {workers} worker(s)..."
    )

    with ProcessPoolExecutor(
        max_workers=workers, initializer=worker_process_initializer
    ) as executor:
        f1 = [executor.submit(_worker_part1_point, t) for t in p1_tasks]
        f2 = [executor.submit(_worker_part2_point, t) for t in p2_tasks]

        for fut in as_completed(f1):
            p_idx, r_idx, s_idx, g_ph, g_en, w_act, n_act, n_tgt = fut.result()
            p1_g_ph_raw[p_idx, r_idx, s_idx] = g_ph
            p1_g_en_raw[p_idx, r_idx, s_idx] = g_en
            p1_w_act_raw[p_idx, r_idx, s_idx] = w_act
            p1_n_ent_active_raw[p_idx, r_idx, s_idx] = float(n_act)
            p1_n_ent_target[p_idx, r_idx] = n_tgt

        for fut in as_completed(f2):
            p_idx, m_idx, r_idx, s_idx, g_ph, g_en, w_act, n_act, n_tgt = fut.result()
            p2_g_ph_raw[p_idx, m_idx, r_idx, s_idx] = g_ph
            p2_g_en_raw[p_idx, m_idx, r_idx, s_idx] = g_en
            p2_w_act_raw[p_idx, m_idx, r_idx, s_idx] = w_act
            p2_n_ent_active_raw[p_idx, m_idx, r_idx, s_idx] = float(n_act)
            p2_n_ent_target[p_idx, m_idx, r_idx] = n_tgt

    # Ensemble statistics for Part 1
    p1_g_ph_stats = compute_ensemble_statistics(p1_g_ph_raw, axis=-1)
    p1_g_en_stats = compute_ensemble_statistics(p1_g_en_raw, axis=-1)
    p1_w_act_stats = compute_ensemble_statistics(p1_w_act_raw, axis=-1)
    p1_ent_act_stats = compute_ensemble_statistics(p1_n_ent_active_raw, axis=-1)

    ge_col = p1_ge_mpa[:, None, None]
    p1_g_ph_norm_stats = compute_ensemble_statistics(p1_g_ph_raw / ge_col, axis=-1)
    p1_g_en_norm_stats = compute_ensemble_statistics(p1_g_en_raw / ge_col, axis=-1)

    # Analytical universal MMT reference curve
    mmt_grid = np.geomspace(0.20, 10.0, 120)
    mmt_g_ph_norm, mmt_g_tot_norm, mmt_ph_pref, mmt_trapping = (
        _compute_universal_mmt_reference(
            mc_over_me_grid=mmt_grid,
            conversion=conversion,
            functionality=functionality,
            stoichiometry=stoichiometry,
        )
    )

    # Ensemble statistics for Part 2
    p2_g_ph_stats = compute_ensemble_statistics(p2_g_ph_raw, axis=-1)
    p2_g_en_stats = compute_ensemble_statistics(p2_g_en_raw, axis=-1)
    p2_w_act_stats = compute_ensemble_statistics(p2_w_act_raw, axis=-1)
    p2_ent_act_stats = compute_ensemble_statistics(p2_n_ent_active_raw, axis=-1)

    return {
        # Global metadata
        "n_repeats": np.array(n_repeats, dtype=int),
        "n_strands": np.array(n_strands, dtype=int),
        "conversion": np.array(conversion, dtype=float),
        "stoichiometry": np.array(stoichiometry, dtype=float),
        "functionality": np.array(functionality, dtype=int),
        # Part 1: Multi-Chemistry Crossover
        "p1_poly_names": p1_poly_names,
        "p1_display_names": p1_display_names,
        "p1_categories": p1_categories,
        "p1_colors": p1_colors,
        "p1_markers": p1_markers,
        "p1_target_ratios": p1_target_ratios,
        "p1_ge_mpa": p1_ge_mpa,
        "p1_me_g_mol": p1_me_g_mol,
        "p1_mk_g_mol": p1_mk_g_mol,
        "p1_mw_g_mol": p1_mw_g_mol,
        "p1_nek": p1_nek,
        "p1_temp_k": p1_temp_k,
        "p1_rho_g_cm3": p1_rho_g_cm3,
        "p1_lk_nm": p1_lk_nm,
        "p1_n_beads": p1_n_beads,
        "p1_actual_mc": p1_actual_mc,
        "p1_actual_mc_over_me": p1_actual_mc_over_me,
        "p1_g_ph_mean": p1_g_ph_stats["mean"],
        "p1_g_ph_ci": p1_g_ph_stats["ci"],
        "p1_g_en_mean": p1_g_en_stats["mean"],
        "p1_g_en_ci": p1_g_en_stats["ci"],
        "p1_g_ph_norm_mean": p1_g_ph_norm_stats["mean"],
        "p1_g_ph_norm_ci": p1_g_ph_norm_stats["ci"],
        "p1_g_en_norm_mean": p1_g_en_norm_stats["mean"],
        "p1_g_en_norm_ci": p1_g_en_norm_stats["ci"],
        "p1_w_act_mean": p1_w_act_stats["mean"],
        "p1_n_ent_active_mean": p1_ent_act_stats["mean"],
        "p1_n_ent_target": p1_n_ent_target,
        "p1_mmt_grid": mmt_grid,
        "p1_mmt_g_ph_norm": mmt_g_ph_norm,
        "p1_mmt_g_tot_norm": mmt_g_tot_norm,
        "p1_mmt_ph_prefactor": np.array(mmt_ph_pref, dtype=float),
        "p1_mmt_trapping": np.array(mmt_trapping, dtype=float),
        # Part 2: Coarse-Graining Discretization Invariance & Breakdown
        "p2_poly_names": p2_poly_names,
        "p2_model_names": p2_model_names,
        "p2_model_labels": p2_model_labels,
        "p2_mc_over_mk": p2_mc_over_mk,
        "p2_mk_g_mol": p2_mk_g_mol,
        "p2_me_g_mol": p2_me_g_mol,
        "p2_ge_mpa": p2_ge_mpa,
        "p2_mw_g_mol": p2_mw_g_mol,
        "p2_du_nm": p2_du_nm,
        "p2_n_beads": p2_n_beads,
        "p2_actual_mc": p2_actual_mc,
        "p2_target_mc": p2_target_mc,
        "p2_vol_nm3": p2_vol_nm3,
        "p2_g_ph_mean": p2_g_ph_stats["mean"],
        "p2_g_ph_ci": p2_g_ph_stats["ci"],
        "p2_g_en_mean": p2_g_en_stats["mean"],
        "p2_g_en_ci": p2_g_en_stats["ci"],
        "p2_w_act_mean": p2_w_act_stats["mean"],
        "p2_n_ent_active_mean": p2_ent_act_stats["mean"],
        "p2_n_ent_target": p2_n_ent_target,
    }


def main() -> None:
    """CLI entry point to compute or refresh the Study 6 cache."""
    cache_path = get_example_cache_path(
        "multichemistry_and_coarse_graining_cache.npz", "force_balance"
    )
    load_or_compute_cache(
        cache_path,
        run_multichemistry_study,
        recompute="--recompute" in sys.argv,
    )


if __name__ == "__main__":
    main()
