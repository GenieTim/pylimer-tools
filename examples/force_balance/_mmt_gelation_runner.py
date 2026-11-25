"""
Multi-seed simulation and analytical MMT runner for Study 1:
Mean-Field Theory (MMT) vs. Spatial MC + MEHP Across Gelation & Loop Regimes.

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
from typing import Any, Dict, Optional, Sequence, Tuple

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
    get_parameters_for_polymer,
)
import pylimer_tools_cpp as pt  # noqa: E402


class LoopRegime(Enum):
    """Loop formation regimes for Monte Carlo network synthesis."""

    IDEAL_NO_LOOPS = "ideal_no_loops"
    STANDARD_SPATIAL = "standard_spatial"
    DILUTED_5X = "diluted_5x"


@dataclass(frozen=True)
class _LoopRegimeConfig:
    """Configuration parameters for a loop formation regime."""

    regime: LoopRegime
    label: str
    primary_loop_prob: float
    secondary_loop_prob: float
    volume_multiplier: float


_LOOP_REGIME_CONFIGS: Tuple[_LoopRegimeConfig, ...] = (
    _LoopRegimeConfig(
        regime=LoopRegime.IDEAL_NO_LOOPS,
        label="Ideal (Loop-Suppressed MC)",
        primary_loop_prob=0.0,
        secondary_loop_prob=0.0,
        volume_multiplier=1.0,
    ),
    _LoopRegimeConfig(
        regime=LoopRegime.STANDARD_SPATIAL,
        label="Standard Spatial MC (Melt)",
        primary_loop_prob=1.0,
        secondary_loop_prob=1.0,
        volume_multiplier=1.0,
    ),
    _LoopRegimeConfig(
        regime=LoopRegime.DILUTED_5X,
        label="Diluted Precursor (5x Volume)",
        primary_loop_prob=1.0,
        secondary_loop_prob=1.0,
        volume_multiplier=5.0,
    ),
)


@dataclass(frozen=True)
class _PdmsConstants:
    """Numeric PDMS Gaussian bead-spring constants derived from ParameterProvider."""

    b: float
    b2_ms: float
    r02: float
    bead_density: float
    ent_density: float
    cutoff: float
    gamma_conv: float
    ge1_mpa: float
    temperature_k: float


def _load_pdms_constants() -> Tuple[Parameters, _PdmsConstants]:
    """Load PDMS Gaussian parameters and extract scalar simulation constants."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    dist_unit = params.get("distance_units")
    constants = _PdmsConstants(
        b=float(params.get("<b>").to(dist_unit).magnitude),
        b2_ms=float(params.get("<b^2>").to(dist_unit**2).magnitude),
        r02=float(params.get("R02").to(dist_unit**2).magnitude),
        bead_density=float(params.get_bead_density()),
        ent_density=float(params.get_entanglement_density()),
        cutoff=float(params.get_sampling_cutoff()),
        gamma_conv=float(params.get_gamma_conversion_factor().to("MPa").magnitude),
        ge1_mpa=float(params.get("Ge").to("MPa").magnitude),
        temperature_k=float(params.get("T").to("kelvin").magnitude),
    )
    return params, constants


def _make_linear_solver_config() -> pt.MEHPLinearSolverConfig:
    """Construct the standard linear MEHP solver configuration."""
    cfg = pt.MEHPLinearSolverConfig()
    cfg.simplification_mode = pt.StructureSimplificationMode.INACTIVE_THEN_X2F
    return cfg


def _build_mc_network(
    seed: int,
    box_l: float,
    constants: _PdmsConstants,
    n_xlinks: int,
    f: int,
    n_b2: int,
    strand_len: int,
    conversion: float,
    n_b1: int = 0,
    primary_loop_prob: float = 1.0,
    secondary_loop_prob: float = 1.0,
) -> pt.MCUniverseGenerator:
    """Build and link an end-linked MC polymer network."""
    gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(seed)
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(constants.b)
    gen.set_mean_squared_bead_distance(constants.b2_ms)
    gen.config_primary_loop_probability(primary_loop_prob)
    gen.config_secondary_loop_probability(secondary_loop_prob)
    gen.add_crosslinkers(n_xlinks, f)
    if n_b2 > 0:
        gen.add_strands(n_b2, [strand_len] * n_b2, strand_atom_type=1)
    if n_b1 > 0:
        gen.add_monofunctional_strands(
            n_b1, [strand_len] * n_b1, monofunctional_strand_atom_type=4
        )
    gen.link_strands_to_conversion(conversion)
    return gen


def _evaluate_mehp_phantom_and_entangled(
    gen: pt.MCUniverseGenerator,
    constants: _PdmsConstants,
    ent_seed: str,
    solver_cfg: pt.MEHPLinearSolverConfig,
    evaluate_entangled: bool = True,
) -> Tuple[float, float, float, float, float]:
    """Evaluate phantom and entangled MEHP moduli and topological weight fractions."""
    vol = float(gen.get_box().get_volume())

    fb_ph = gen.get_mehp_simulator(
        entanglement_mode=pt.EntanglementMode.NONE,
        crosslinker_type=2,
    )
    fb_ph.run(solver_cfg)
    g_ph = float(
        constants.gamma_conv * np.sum(fb_ph.get_gamma_factors(constants.r02)) / vol
    )
    w_sol = float(fb_ph.get_soluble_weight_fraction())
    w_dang = float(fb_ph.get_dangling_weight_fraction())
    w_act = float(fb_ph.get_active_weight_fraction())

    if not evaluate_entangled:
        return g_ph, math.nan, w_sol, w_dang, w_act

    n_ent_target = int(round(constants.ent_density * vol))
    es = gen.sample_entanglements(
        n_ent_target,
        constants.cutoff,
        0.0,
        0,
        3.0,
        ent_seed,
        pt.EntanglementFilteringMode.OMIT_INACTIVE,
    )
    fb_en = gen.get_mehp_simulator(
        es,
        entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
        crosslinker_type=2,
    )
    fb_en.run(solver_cfg)
    g_en = float(
        constants.gamma_conv * np.sum(fb_en.get_gamma_factors(constants.r02)) / vol
    )
    return g_ph, g_en, w_sol, w_dang, w_act


def _compute_mmt_state(
    params: Parameters,
    r: float,
    p: float,
    f: int,
    n_b2: int,
    n_b1: int,
    n_xlinks: int,
    strand_len: int,
    vol: float,
    b2: float = 1.0,
) -> Dict[str, float]:
    """Evaluate analytical MMT quantities for a monodisperse A_f + B_2/B_1 system."""
    ureg = params.get_unit_registry()
    temperature = params.get("T")
    ge1 = params.get("Ge")

    n_atoms = (n_b2 + n_b1) * strand_len + n_xlinks
    weight_fractions: Dict[int, float] = {
        1: (n_b2 * strand_len) / n_atoms,
        2: n_xlinks / n_atoms,
    }
    functionality_per_type: Dict[int, int] = {1: 2, 2: f}
    if n_b1 > 0:
        weight_fractions[4] = (n_b1 * strand_len) / n_atoms
        functionality_per_type[4] = 1

    p_gel = mmt.predict_gelation_point(r=r, f=f, b2=b2)
    xlink_conc = (n_xlinks / vol) * (ureg("nm") ** -3)
    nu_b2_q = (n_b2 / vol) * (ureg("nm") ** -3)

    if p <= p_gel:
        return {
            "g_ph": 0.0,
            "g_en": 0.0,
            "g_tot": 0.0,
            "g_junc_exact": 0.0,
            "g_ent_patel": 0.0,
            "g_tot_patel": 0.0,
            "g_ext_builtin": 0.0,
            "g_ext_patel": 0.0,
            "w_sol": 1.0,
            "w_dang": 0.0,
            "w_act": 0.0,
        }

    g_mmt_ph, g_mmt_en, _, _ = mmt.compute_modulus_decomposition(
        network=None,
        ureg=ureg,
        r=r,
        p=p,
        f=f,
        nu=nu_b2_q,
        temperature=temperature,
        g_e_1=ge1,
        b2=b2,
    )
    g_ph_val = float(g_mmt_ph.to("MPa").magnitude)
    g_en_val = float(g_mmt_en.to("MPa").magnitude)

    g_junc_exact = float(
        mmt.compute_junction_modulus(
            p=p,
            r=r,
            xlink_concentration_0=xlink_conc,
            ureg=ureg,
            f=f,
            temperature=temperature,
            b2=b2,
        )
        .to("MPa")
        .magnitude
    )

    w_sol = float(
        mmt.compute_weight_fraction_of_soluble_material(
            network=None,
            crosslinker_type=2,
            functionality_per_type=functionality_per_type,
            weight_fractions=weight_fractions,
            r=r,
            p=p,
            b2=b2,
        )
    )
    w_dang = float(
        mmt.compute_weight_fraction_of_dangling_chains(
            network=None,
            crosslinker_type=2,
            functionality_per_type=functionality_per_type,
            weight_fractions=weight_fractions,
            r=r,
            p=p,
            b2=b2,
        )
    )
    w_act = float(
        mmt.compute_weight_fraction_of_backbone(
            network=None,
            crosslinker_type=2,
            functionality_per_type=functionality_per_type,
            weight_fractions=weight_fractions,
            r=r,
            p=p,
            b2=b2,
        )
    )

    # Patel et al. (1992) trapping factor includes w_B2^2 weight dilution of B2 chains
    w_b2_polymer = n_b2 / (n_b2 + n_b1)
    g_ent_patel = (w_b2_polymer**2) * g_en_val

    g_ext_builtin = float(
        mmt.compute_extracted_modulus(
            p=p,
            r=r,
            f=f,
            g_e_1=ge1,
            w_sol=w_sol,
            xlink_concentration_0=xlink_conc,
            ureg=ureg,
            temperature=temperature,
            b2=b2,
        )
        .to("MPa")
        .magnitude
    )
    g_ext_patel = (1.0 - w_sol) ** (-1.0 / 3.0) * g_junc_exact + (1.0 - w_sol) ** (
        -2.0
    ) * g_ent_patel

    return {
        "g_ph": g_ph_val,
        "g_en": g_en_val,
        "g_tot": g_ph_val + g_en_val,
        "g_junc_exact": g_junc_exact,
        "g_ent_patel": g_ent_patel,
        "g_tot_patel": g_junc_exact + g_ent_patel,
        "g_ext_builtin": g_ext_builtin,
        "g_ext_patel": g_ext_patel,
        "w_sol": w_sol,
        "w_dang": w_dang,
        "w_act": w_act,
    }


def _simulate_single_seed(
    seed_idx: int,
    base_seed: int,
    constants: _PdmsConstants,
    strand_len: int,
    f: int,
    p_sweep: Sequence[float],
    n_strands_p: int,
    r_targets: Sequence[float],
    n_strands_r: int,
    loop_p_sweep: Sequence[float],
    n_strands_loop: int,
    b2_targets: Sequence[float],
    total_b_sites_b2: int,
) -> Dict[str, Any]:
    """Run all 4 sub-studies for a single Monte Carlo random seed."""
    solver_cfg = _make_linear_solver_config()
    seed_offset = base_seed + 1009 * seed_idx

    # --- Part 1a: Conversion p sweep across p_c at r = 1.0 ---
    n_p = len(p_sweep)
    p_g_ph = np.zeros(n_p)
    p_g_en = np.zeros(n_p)
    p_w_sol = np.zeros(n_p)
    p_w_dang = np.zeros(n_p)
    p_w_act = np.zeros(n_p)

    n_x_p = n_strands_p // 2
    vol_p = (n_strands_p * strand_len + n_x_p) / constants.bead_density
    box_l_p = vol_p ** (1.0 / 3.0)

    for i, p_val in enumerate(p_sweep):
        gen = _build_mc_network(
            seed=seed_offset + i,
            box_l=box_l_p,
            constants=constants,
            n_xlinks=n_x_p,
            f=f,
            n_b2=n_strands_p,
            strand_len=strand_len,
            conversion=float(p_val),
        )
        g_ph, g_en, w_sol, w_dang, w_act = _evaluate_mehp_phantom_and_entangled(
            gen,
            constants,
            ent_seed=f"psweep_{seed_idx}_{i}",
            solver_cfg=solver_cfg,
            evaluate_entangled=True,
        )
        p_g_ph[i] = g_ph
        p_g_en[i] = g_en
        p_w_sol[i] = w_sol
        p_w_dang[i] = w_dang
        p_w_act[i] = w_act

    # --- Part 1b: Stoichiometric imbalance r sweep at p = 0.95 * p_max(r, f) ---
    n_r = len(r_targets)
    r_g_ph = np.zeros(n_r)
    r_g_en = np.zeros(n_r)
    r_w_sol = np.zeros(n_r)
    r_w_dang = np.zeros(n_r)
    r_w_act = np.zeros(n_r)

    for i, r_val in enumerate(r_targets):
        n_x_r = int(round(r_val * n_strands_r / 2.0))
        r_actual = (f * n_x_r) / (2.0 * n_strands_r)
        p_max = mmt.predict_maximum_p(r_actual, f)
        assert p_max is not None
        p_target = 0.95 * p_max
        vol_r = (n_strands_r * strand_len + n_x_r) / constants.bead_density
        box_l_r = vol_r ** (1.0 / 3.0)

        gen = _build_mc_network(
            seed=seed_offset + 100 + i,
            box_l=box_l_r,
            constants=constants,
            n_xlinks=n_x_r,
            f=f,
            n_b2=n_strands_r,
            strand_len=strand_len,
            conversion=p_target,
        )
        g_ph, g_en, w_sol, w_dang, w_act = _evaluate_mehp_phantom_and_entangled(
            gen,
            constants,
            ent_seed=f"rsweep_{seed_idx}_{i}",
            solver_cfg=solver_cfg,
            evaluate_entangled=True,
        )
        r_g_ph[i] = g_ph
        r_g_en[i] = g_en
        r_w_sol[i] = w_sol
        r_w_dang[i] = w_dang
        r_w_act[i] = w_act

    # --- Part 2: Loop & dilution regime across p near p_c ---
    n_regimes = len(_LOOP_REGIME_CONFIGS)
    n_lp = len(loop_p_sweep)
    loop_g_ph = np.zeros((n_regimes, n_lp))
    loop_g_ph_norm = np.zeros((n_regimes, n_lp))
    loop_w_act = np.zeros((n_regimes, n_lp))
    loop_w_sol = np.zeros((n_regimes, n_lp))
    loop_primary = np.zeros((n_regimes, n_lp))
    loop_secondary = np.zeros((n_regimes, n_lp))

    n_x_l = n_strands_loop // 2
    vol_melt_l = (n_strands_loop * strand_len + n_x_l) / constants.bead_density
    prim_len = strand_len + 1
    sec_len = 2 * (strand_len + 1)

    for reg_idx, cfg in enumerate(_LOOP_REGIME_CONFIGS):
        vol_l = vol_melt_l * cfg.volume_multiplier
        box_l_l = vol_l ** (1.0 / 3.0)
        for j, p_val in enumerate(loop_p_sweep):
            gen = _build_mc_network(
                seed=seed_offset + 200 + reg_idx * 50 + j,
                box_l=box_l_l,
                constants=constants,
                n_xlinks=n_x_l,
                f=f,
                n_b2=n_strands_loop,
                strand_len=strand_len,
                conversion=float(p_val),
                primary_loop_prob=cfg.primary_loop_prob,
                secondary_loop_prob=cfg.secondary_loop_prob,
            )
            g_ph, _, w_sol, _, w_act = _evaluate_mehp_phantom_and_entangled(
                gen,
                constants,
                ent_seed="",
                solver_cfg=solver_cfg,
                evaluate_entangled=False,
            )
            u = gen.get_universe()
            loop_counts = u.count_loop_lengths(max_length=sec_len)
            loop_g_ph[reg_idx, j] = g_ph
            loop_g_ph_norm[reg_idx, j] = g_ph * cfg.volume_multiplier
            loop_w_act[reg_idx, j] = w_act
            loop_w_sol[reg_idx, j] = w_sol
            loop_primary[reg_idx, j] = float(loop_counts.get(prim_len, 0))
            loop_secondary[reg_idx, j] = float(loop_counts.get(sec_len, 0))

    # --- Part 3: Monofunctional chains (b2 < 1) & post-extraction deswelling ---
    n_b2_pts = len(b2_targets)
    b2_g_ph_pre = np.zeros(n_b2_pts)
    b2_g_en_pre = np.zeros(n_b2_pts)
    b2_g_ph_ext_rescaled = np.zeros(n_b2_pts)
    b2_g_en_ext_rescaled = np.zeros(n_b2_pts)
    b2_g_ph_ext_fixed = np.zeros(n_b2_pts)
    b2_w_sol = np.zeros(n_b2_pts)
    b2_w_dang = np.zeros(n_b2_pts)
    b2_w_act = np.zeros(n_b2_pts)

    n_x_b2 = total_b_sites_b2 // f  # keeps r = 1.0 exactly
    for i, b2_val in enumerate(b2_targets):
        n_b2 = int(round(b2_val * total_b_sites_b2 / 2.0))
        n_b1 = total_b_sites_b2 - 2 * n_b2
        vol_b2 = ((n_b2 + n_b1) * strand_len + n_x_b2) / constants.bead_density
        box_l_b2 = vol_b2 ** (1.0 / 3.0)
        net_seed = seed_offset + 500 + i

        gen_rescaled = _build_mc_network(
            seed=net_seed,
            box_l=box_l_b2,
            constants=constants,
            n_xlinks=n_x_b2,
            f=f,
            n_b2=n_b2,
            strand_len=strand_len,
            conversion=0.95,
            n_b1=n_b1,
        )
        g_ph_pre, g_en_pre, w_sol, w_dang, w_act = _evaluate_mehp_phantom_and_entangled(
            gen_rescaled,
            constants,
            ent_seed=f"b2pre_{seed_idx}_{i}",
            solver_cfg=solver_cfg,
            evaluate_entangled=True,
        )
        b2_g_ph_pre[i] = g_ph_pre
        b2_g_en_pre[i] = g_en_pre
        b2_w_sol[i] = w_sol
        b2_w_dang[i] = w_dang
        b2_w_act[i] = w_act

        # Post-extraction with deswelling (rescale_box=True)
        gen_rescaled.remove_soluble_fraction(rescale_box=True)
        g_ph_ext_r, g_en_ext_r, _, _, _ = _evaluate_mehp_phantom_and_entangled(
            gen_rescaled,
            constants,
            ent_seed=f"b2ext_{seed_idx}_{i}",
            solver_cfg=solver_cfg,
            evaluate_entangled=True,
        )
        b2_g_ph_ext_rescaled[i] = g_ph_ext_r
        b2_g_en_ext_rescaled[i] = g_en_ext_r

        # Post-extraction without box rescaling (rescale_box=False)
        gen_fixed = _build_mc_network(
            seed=net_seed,
            box_l=box_l_b2,
            constants=constants,
            n_xlinks=n_x_b2,
            f=f,
            n_b2=n_b2,
            strand_len=strand_len,
            conversion=0.95,
            n_b1=n_b1,
        )
        gen_fixed.remove_soluble_fraction(rescale_box=False)
        g_ph_ext_f, _, _, _, _ = _evaluate_mehp_phantom_and_entangled(
            gen_fixed,
            constants,
            ent_seed="",
            solver_cfg=solver_cfg,
            evaluate_entangled=False,
        )
        b2_g_ph_ext_fixed[i] = g_ph_ext_f

    return {
        "seed_idx": seed_idx,
        "p_g_ph": p_g_ph,
        "p_g_en": p_g_en,
        "p_w_sol": p_w_sol,
        "p_w_dang": p_w_dang,
        "p_w_act": p_w_act,
        "r_g_ph": r_g_ph,
        "r_g_en": r_g_en,
        "r_w_sol": r_w_sol,
        "r_w_dang": r_w_dang,
        "r_w_act": r_w_act,
        "loop_g_ph": loop_g_ph,
        "loop_g_ph_norm": loop_g_ph_norm,
        "loop_w_act": loop_w_act,
        "loop_w_sol": loop_w_sol,
        "loop_primary": loop_primary,
        "loop_secondary": loop_secondary,
        "b2_g_ph_pre": b2_g_ph_pre,
        "b2_g_en_pre": b2_g_en_pre,
        "b2_g_ph_ext_rescaled": b2_g_ph_ext_rescaled,
        "b2_g_en_ext_rescaled": b2_g_en_ext_rescaled,
        "b2_g_ph_ext_fixed": b2_g_ph_ext_fixed,
        "b2_w_sol": b2_w_sol,
        "b2_w_dang": b2_w_dang,
        "b2_w_act": b2_w_act,
    }


def _pack_stats(prefix: str, arr: np.ndarray, axis: int = -1) -> Dict[str, np.ndarray]:
    """Compute ensemble mean/std/ci and return flattened dict keys."""
    stats = compute_ensemble_statistics(arr, axis=axis)
    return {
        f"{prefix}_mean": stats["mean"],
        f"{prefix}_std": stats["std"],
        f"{prefix}_ci": stats["ci"],
    }


def run_mmt_vs_mehp_study(
    n_repeats: int = 8,
    strand_len: int = 50,
    n_strands: int = 1500,
    base_seed: int = 2026,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Run Study 1: MMT vs. Spatial MC + MEHP across gelation, stoichiometry, loops, and b2."""
    params, constants = _load_pdms_constants()
    f = 4

    p_sweep = np.array(
        [0.45, 0.50, 0.55, 0.577, 0.60, 0.63, 0.67, 0.72, 0.78, 0.85, 0.92, 0.98],
        dtype=float,
    )
    r_targets = np.array([0.45, 0.65, 0.80, 1.00, 1.25, 1.60, 2.00], dtype=float)
    loop_p_sweep = np.array(
        [0.50, 0.55, 0.58, 0.61, 0.65, 0.70, 0.78, 0.86, 0.95], dtype=float
    )
    b2_targets = np.array([0.55, 0.65, 0.75, 0.85, 0.95, 1.00], dtype=float)

    n_strands_p = n_strands
    n_strands_r = n_strands
    n_strands_loop = n_strands
    total_b_sites_b2 = (
        2 * n_strands
    )  # so b2=1.0 corresponds to n_strands bifunctional chains

    # --- Precompute fine-grid and point-matched analytical MMT ---
    p_fine = np.linspace(0.45, 0.98, 200)
    n_x_p = n_strands_p // 2
    vol_p = (n_strands_p * strand_len + n_x_p) / constants.bead_density
    p_mmt_fine_ph = np.zeros(len(p_fine))
    p_mmt_fine_tot = np.zeros(len(p_fine))
    p_mmt_fine_w_sol = np.zeros(len(p_fine))
    p_mmt_fine_w_dang = np.zeros(len(p_fine))
    p_mmt_fine_w_act = np.zeros(len(p_fine))
    for i, pv in enumerate(p_fine):
        res = _compute_mmt_state(
            params,
            r=1.0,
            p=float(pv),
            f=f,
            n_b2=n_strands_p,
            n_b1=0,
            n_xlinks=n_x_p,
            strand_len=strand_len,
            vol=vol_p,
        )
        p_mmt_fine_ph[i] = res["g_ph"]
        p_mmt_fine_tot[i] = res["g_tot"]
        p_mmt_fine_w_sol[i] = res["w_sol"]
        p_mmt_fine_w_dang[i] = res["w_dang"]
        p_mmt_fine_w_act[i] = res["w_act"]

    r_actual = np.array(
        [
            (f * int(round(rv * n_strands_r / 2.0))) / (2.0 * n_strands_r)
            for rv in r_targets
        ],
        dtype=float,
    )
    r_fine = np.linspace(0.40, 2.05, 150)
    r_mmt_fine_ph = np.zeros(len(r_fine))
    r_mmt_fine_tot = np.zeros(len(r_fine))
    r_mmt_fine_w_sol = np.zeros(len(r_fine))
    r_mmt_fine_w_dang = np.zeros(len(r_fine))
    r_mmt_fine_w_act = np.zeros(len(r_fine))
    for i, rv in enumerate(r_fine):
        p_max = mmt.predict_maximum_p(float(rv), f)
        assert p_max is not None
        p_val = 0.95 * p_max
        n_x_exact = rv * n_strands_r / 2.0
        vol_exact = (n_strands_r * strand_len + n_x_exact) / constants.bead_density
        res = _compute_mmt_state(
            params,
            r=float(rv),
            p=p_val,
            f=f,
            n_b2=n_strands_r,
            n_b1=0,
            n_xlinks=n_x_exact,
            strand_len=strand_len,
            vol=vol_exact,  # type: ignore[arg-type]
        )
        r_mmt_fine_ph[i] = res["g_ph"]
        r_mmt_fine_tot[i] = res["g_tot"]
        r_mmt_fine_w_sol[i] = res["w_sol"]
        r_mmt_fine_w_dang[i] = res["w_dang"]
        r_mmt_fine_w_act[i] = res["w_act"]

    loop_mmt_ph = np.zeros(len(loop_p_sweep))
    loop_mmt_w_act = np.zeros(len(loop_p_sweep))
    n_x_l = n_strands_loop // 2
    vol_l = (n_strands_loop * strand_len + n_x_l) / constants.bead_density
    for i, pv in enumerate(loop_p_sweep):
        res = _compute_mmt_state(
            params,
            r=1.0,
            p=float(pv),
            f=f,
            n_b2=n_strands_loop,
            n_b1=0,
            n_xlinks=n_x_l,
            strand_len=strand_len,
            vol=vol_l,
        )
        loop_mmt_ph[i] = res["g_ph"]
        loop_mmt_w_act[i] = res["w_act"]

    n_x_b2 = total_b_sites_b2 // f
    b2_actual = np.zeros(len(b2_targets))
    b2_mmt_ph_junc = np.zeros(len(b2_targets))
    b2_mmt_tot_builtin = np.zeros(len(b2_targets))
    b2_mmt_tot_patel = np.zeros(len(b2_targets))
    b2_mmt_ext_builtin = np.zeros(len(b2_targets))
    b2_mmt_ext_patel = np.zeros(len(b2_targets))
    b2_mmt_w_sol = np.zeros(len(b2_targets))
    b2_mmt_w_dang = np.zeros(len(b2_targets))
    b2_mmt_w_act = np.zeros(len(b2_targets))
    for i, b2v in enumerate(b2_targets):
        n_b2 = int(round(b2v * total_b_sites_b2 / 2.0))
        n_b1 = total_b_sites_b2 - 2 * n_b2
        b2_act = (2.0 * n_b2) / (2.0 * n_b2 + n_b1)
        b2_actual[i] = b2_act
        vol_b2 = ((n_b2 + n_b1) * strand_len + n_x_b2) / constants.bead_density
        res = _compute_mmt_state(
            params,
            r=1.0,
            p=0.95,
            f=f,
            n_b2=n_b2,
            n_b1=n_b1,
            n_xlinks=n_x_b2,
            strand_len=strand_len,
            vol=vol_b2,
            b2=b2_act,
        )
        b2_mmt_ph_junc[i] = res["g_junc_exact"]
        b2_mmt_tot_builtin[i] = res["g_junc_exact"] + res["g_en"]
        b2_mmt_tot_patel[i] = res["g_tot_patel"]
        b2_mmt_ext_builtin[i] = res["g_ext_builtin"]
        b2_mmt_ext_patel[i] = res["g_ext_patel"]
        b2_mmt_w_sol[i] = res["w_sol"]
        b2_mmt_w_dang[i] = res["w_dang"]
        b2_mmt_w_act[i] = res["w_act"]

    # --- Execute multi-seed Monte Carlo + MEHP ensemble ---
    n_p = len(p_sweep)
    n_r = len(r_targets)
    n_regimes = len(_LOOP_REGIME_CONFIGS)
    n_lp = len(loop_p_sweep)
    n_b2_pts = len(b2_targets)

    all_p_g_ph = np.zeros((n_p, n_repeats))
    all_p_g_en = np.zeros((n_p, n_repeats))
    all_p_w_sol = np.zeros((n_p, n_repeats))
    all_p_w_dang = np.zeros((n_p, n_repeats))
    all_p_w_act = np.zeros((n_p, n_repeats))

    all_r_g_ph = np.zeros((n_r, n_repeats))
    all_r_g_en = np.zeros((n_r, n_repeats))
    all_r_w_sol = np.zeros((n_r, n_repeats))
    all_r_w_dang = np.zeros((n_r, n_repeats))
    all_r_w_act = np.zeros((n_r, n_repeats))

    all_loop_g_ph = np.zeros((n_regimes, n_lp, n_repeats))
    all_loop_g_ph_norm = np.zeros((n_regimes, n_lp, n_repeats))
    all_loop_w_act = np.zeros((n_regimes, n_lp, n_repeats))
    all_loop_w_sol = np.zeros((n_regimes, n_lp, n_repeats))
    all_loop_prim = np.zeros((n_regimes, n_lp, n_repeats))
    all_loop_sec = np.zeros((n_regimes, n_lp, n_repeats))

    all_b2_g_ph_pre = np.zeros((n_b2_pts, n_repeats))
    all_b2_g_en_pre = np.zeros((n_b2_pts, n_repeats))
    all_b2_g_ph_ext_rescaled = np.zeros((n_b2_pts, n_repeats))
    all_b2_g_en_ext_rescaled = np.zeros((n_b2_pts, n_repeats))
    all_b2_g_ph_ext_fixed = np.zeros((n_b2_pts, n_repeats))
    all_b2_w_sol = np.zeros((n_b2_pts, n_repeats))
    all_b2_w_dang = np.zeros((n_b2_pts, n_repeats))
    all_b2_w_act = np.zeros((n_b2_pts, n_repeats))

    tasks = [
        (
            seed_idx,
            base_seed,
            constants,
            strand_len,
            f,
            p_sweep,
            n_strands_p,
            r_targets,
            n_strands_r,
            loop_p_sweep,
            n_strands_loop,
            b2_targets,
            total_b_sites_b2,
        )
        for seed_idx in range(n_repeats)
    ]

    workers = min(resolve_max_workers(max_workers), n_repeats)
    configure_process_threading(workers)

    def _store_seed_result(res: Dict[str, Any]) -> None:
        idx = int(res["seed_idx"])
        all_p_g_ph[:, idx] = res["p_g_ph"]
        all_p_g_en[:, idx] = res["p_g_en"]
        all_p_w_sol[:, idx] = res["p_w_sol"]
        all_p_w_dang[:, idx] = res["p_w_dang"]
        all_p_w_act[:, idx] = res["p_w_act"]

        all_r_g_ph[:, idx] = res["r_g_ph"]
        all_r_g_en[:, idx] = res["r_g_en"]
        all_r_w_sol[:, idx] = res["r_w_sol"]
        all_r_w_dang[:, idx] = res["r_w_dang"]
        all_r_w_act[:, idx] = res["r_w_act"]

        all_loop_g_ph[:, :, idx] = res["loop_g_ph"]
        all_loop_g_ph_norm[:, :, idx] = res["loop_g_ph_norm"]
        all_loop_w_act[:, :, idx] = res["loop_w_act"]
        all_loop_w_sol[:, :, idx] = res["loop_w_sol"]
        all_loop_prim[:, :, idx] = res["loop_primary"]
        all_loop_sec[:, :, idx] = res["loop_secondary"]

        all_b2_g_ph_pre[:, idx] = res["b2_g_ph_pre"]
        all_b2_g_en_pre[:, idx] = res["b2_g_en_pre"]
        all_b2_g_ph_ext_rescaled[:, idx] = res["b2_g_ph_ext_rescaled"]
        all_b2_g_en_ext_rescaled[:, idx] = res["b2_g_en_ext_rescaled"]
        all_b2_g_ph_ext_fixed[:, idx] = res["b2_g_ph_ext_fixed"]
        all_b2_w_sol[:, idx] = res["b2_w_sol"]
        all_b2_w_dang[:, idx] = res["b2_w_dang"]
        all_b2_w_act[:, idx] = res["b2_w_act"]

    if workers > 1 and n_repeats > 1:
        print(f"Executing {n_repeats} seeds in parallel across {workers} workers...")
        completed = 0
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = [executor.submit(_simulate_single_seed, *t) for t in tasks]
            for fut in as_completed(futures):
                _store_seed_result(fut.result())
                completed += 1
                print(f"  [MMT vs MEHP] Completed seed {completed}/{n_repeats}...")
    else:
        for t in tasks:
            _store_seed_result(_simulate_single_seed(*t))
            print(f"  [MMT vs MEHP] Completed seed {t[0] + 1}/{n_repeats}...")

    result: Dict[str, Any] = {
        "n_repeats": np.array(n_repeats),
        "strand_len": np.array(strand_len),
        "n_strands": np.array(n_strands),
        "f": np.array(f),
        "p_c_stoich": np.array(mmt.predict_gelation_point(r=1.0, f=f, b2=1.0)),
        "box_l_melt": np.array(vol_p ** (1.0 / 3.0)),
        # Part 1a: p sweep
        "p_sweep": p_sweep,
        "p_fine": p_fine,
        "p_mmt_fine_ph": p_mmt_fine_ph,
        "p_mmt_fine_tot": p_mmt_fine_tot,
        "p_mmt_fine_w_sol": p_mmt_fine_w_sol,
        "p_mmt_fine_w_dang": p_mmt_fine_w_dang,
        "p_mmt_fine_w_act": p_mmt_fine_w_act,
        # Part 1b: r sweep
        "r_targets": r_actual,
        "r_fine": r_fine,
        "r_mmt_fine_ph": r_mmt_fine_ph,
        "r_mmt_fine_tot": r_mmt_fine_tot,
        "r_mmt_fine_w_sol": r_mmt_fine_w_sol,
        "r_mmt_fine_w_dang": r_mmt_fine_w_dang,
        "r_mmt_fine_w_act": r_mmt_fine_w_act,
        # Part 2: loops & dilution
        "loop_p_sweep": loop_p_sweep,
        "loop_regime_keys": np.array([c.regime.value for c in _LOOP_REGIME_CONFIGS]),
        "loop_regime_labels": np.array([c.label for c in _LOOP_REGIME_CONFIGS]),
        "loop_mmt_ph": loop_mmt_ph,
        "loop_mmt_w_act": loop_mmt_w_act,
        # Part 3: b2 & extraction
        "b2_targets": b2_actual,
        "b2_mmt_ph_junc": b2_mmt_ph_junc,
        "b2_mmt_tot_builtin": b2_mmt_tot_builtin,
        "b2_mmt_tot_patel": b2_mmt_tot_patel,
        "b2_mmt_ext_builtin": b2_mmt_ext_builtin,
        "b2_mmt_ext_patel": b2_mmt_ext_patel,
        "b2_mmt_w_sol": b2_mmt_w_sol,
        "b2_mmt_w_dang": b2_mmt_w_dang,
        "b2_mmt_w_act": b2_mmt_w_act,
    }

    for prefix, arr in (
        ("p_g_ph", all_p_g_ph),
        ("p_g_en", all_p_g_en),
        ("p_w_sol", all_p_w_sol),
        ("p_w_dang", all_p_w_dang),
        ("p_w_act", all_p_w_act),
        ("r_g_ph", all_r_g_ph),
        ("r_g_en", all_r_g_en),
        ("r_w_sol", all_r_w_sol),
        ("r_w_dang", all_r_w_dang),
        ("r_w_act", all_r_w_act),
        ("loop_g_ph", all_loop_g_ph),
        ("loop_g_ph_norm", all_loop_g_ph_norm),
        ("loop_w_act", all_loop_w_act),
        ("loop_w_sol", all_loop_w_sol),
        ("loop_prim_frac", all_loop_prim / float(n_strands_loop)),
        ("loop_sec_frac", all_loop_sec / float(n_strands_loop)),
        ("b2_g_ph_pre", all_b2_g_ph_pre),
        ("b2_g_en_pre", all_b2_g_en_pre),
        ("b2_g_ph_ext_rescaled", all_b2_g_ph_ext_rescaled),
        ("b2_g_en_ext_rescaled", all_b2_g_en_ext_rescaled),
        ("b2_g_ph_ext_fixed", all_b2_g_ph_ext_fixed),
        ("b2_w_sol", all_b2_w_sol),
        ("b2_w_dang", all_b2_w_dang),
        ("b2_w_act", all_b2_w_act),
    ):
        result.update(_pack_stats(prefix, arr, axis=-1))

    return result


if __name__ == "__main__":
    cache_path = get_example_cache_path(
        "mmt_vs_mehp_gelation_and_loops_cache.npz", "force_balance"
    )
    load_or_compute_cache(
        cache_path,
        run_mmt_vs_mehp_study,
        recompute="--recompute" in sys.argv,
        max_workers=2 if "--workers" not in sys.argv else None,
    )
