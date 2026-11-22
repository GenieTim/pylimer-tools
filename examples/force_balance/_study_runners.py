"""
Multi-seed simulation and benchmark runners for the force_balance examples.

This module is prefixed with an underscore so Sphinx-Gallery ignores it when
building documentation, keeping the example scripts focused on idiomatic API
usage and resulting plots.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from enum import Enum
import json
import math
from pathlib import Path
import random
import subprocess
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

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
    get_parameters_for_polymer,
)
import pylimer_tools_cpp as pt  # noqa: E402
from pylimer_tools_cpp import (  # noqa: E402
    EntanglementFilteringMode,
    MCUniverseGenerator,
)


# --- Entanglement Link vs. Spring ---


def _simulate_link_vs_spring_single_repeat(
    r: int,
    seed: int,
    box_l: float,
    b: float,
    n_strands: int,
    strand_len: int,
    conversion: float,
    density_factors: Sequence[float],
    nominal_n_ent: int,
    cutoff: float,
    r02: float,
    vol: float,
    gamma_conv: float,
) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Execute a single Monte Carlo network realization across entanglement densities."""
    gen = MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(seed)
    gen.set_bead_distance(b)
    gen.set_mean_squared_bead_distance(b**2)
    gen.add_strands(n_strands, [strand_len] * n_strands)
    gen.add_crosslinkers(int(0.5 * n_strands), 4)
    gen.link_strands_to_conversion(conversion)

    n_factors = len(density_factors)
    g_links = np.zeros(n_factors)
    g_springs = np.zeros(n_factors)
    extra_links = np.zeros(n_factors)
    extra_springs = np.zeros(n_factors)

    for i, factor in enumerate(density_factors):
        n_target = int(factor * nominal_n_ent)

        if n_target == 0:
            # Bare unentangled crosslinker network
            fb_bare = gen.get_mehp_simulator()
            fb_bare.run()
            g_pa = (gamma_conv * fb_bare.get_gamma_factor(r02) / vol) * 1e6
            g_links[i] = g_pa
            g_springs[i] = g_pa
            extra_links[i] = 0
            extra_springs[i] = 0
        else:
            # Sample entanglements using seeded Brownian bridge
            es = gen.sample_entanglements(
                n_target,
                cutoff,
                0.0,
                0,
                3.0,
                f"ent_{seed}_{factor}",
                filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
            )

            # 1. Entanglement Links (merged junction nodes)
            fb_l = gen.get_mehp_simulator(
                es,
                entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
                crosslinker_type=2,
            )
            extra_links[i] = fb_l.get_nr_of_extra_atoms()
            fb_l.run()
            g_links[i] = (gamma_conv * fb_l.get_gamma_factor(r02) / vol) * 1e6

            # 2. Entanglement Springs (explicit Hookean coupling springs)
            fb_s = gen.get_mehp_simulator(
                es,
                entanglement_mode=pt.EntanglementMode.SPRING,
                crosslinker_type=2,
            )
            extra_springs[i] = fb_s.get_nr_of_extra_atoms()
            fb_s.run()
            g_springs[i] = (gamma_conv * fb_s.get_gamma_factor(r02) / vol) * 1e6

    return r, g_links, g_springs, extra_links, extra_springs


def run_entanglement_link_vs_spring_simulation(
    density_factors: Sequence[float] = (0.0, 0.5, 1.0, 1.5, 2.0),
    n_strands: int = 10000,
    strand_len: int = 50,
    conversion: float = 0.95,
    n_repeats: int = 10,
    base_seed: int = 123,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Evaluate shear modulus and extra degrees of freedom across entanglement densities.

    Averages over independent network realizations and stochastic entanglement samplings
    to provide statistically robust comparisons with Student's t 95% confidence intervals.
    """
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = params.get("<b>").to(params.get("distance_units")).magnitude
    bead_density = params.get_bead_density()
    cutoff = params.get_sampling_cutoff()
    nominal_ent_density = params.get_entanglement_density()
    r02 = params.get("R02").to(params.get("distance_units") ** 2).magnitude
    kbt = params.get("T") * params.get("kb")
    gamma_conv = ((kbt / ((params.get("distance_units")) ** 3)).to("MPa")).magnitude

    n_atoms = n_strands * strand_len + int(0.5 * n_strands)
    vol = n_atoms / bead_density
    box_l = vol ** (1.0 / 3.0)
    nominal_n_ent = int(nominal_ent_density * vol)

    all_g_links = np.zeros((len(density_factors), n_repeats))
    all_g_springs = np.zeros((len(density_factors), n_repeats))
    all_extra_links = np.zeros((len(density_factors), n_repeats))
    all_extra_springs = np.zeros((len(density_factors), n_repeats))

    tasks = [
        (
            r,
            base_seed + 101 * r,
            box_l,
            b,
            n_strands,
            strand_len,
            conversion,
            density_factors,
            nominal_n_ent,
            cutoff,
            r02,
            vol,
            gamma_conv,
        )
        for r in range(n_repeats)
    ]

    workers = min(resolve_max_workers(max_workers), n_repeats)
    configure_process_threading(workers)

    if workers > 1 and n_repeats > 1:
        print(
            f"Executing {n_repeats} repeats in parallel across {workers} worker processes..."
        )
        completed = 0
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = [
                executor.submit(_simulate_link_vs_spring_single_repeat, *task)
                for task in tasks
            ]
            for fut in as_completed(futures):
                r, g_l, g_s, ex_l, ex_s = fut.result()
                all_g_links[:, r] = g_l
                all_g_springs[:, r] = g_s
                all_extra_links[:, r] = ex_l
                all_extra_springs[:, r] = ex_s
                completed += 1
                print(
                    f"  [Link vs Spring] Completed repeat {completed}/{
                        n_repeats
                    } (repeat index {r + 1})..."
                )
    else:
        for task in tasks:
            r, g_l, g_s, ex_l, ex_s = _simulate_link_vs_spring_single_repeat(*task)
            all_g_links[:, r] = g_l
            all_g_springs[:, r] = g_s
            all_extra_links[:, r] = ex_l
            all_extra_springs[:, r] = ex_s
            print(f"  [Link vs Spring] Completed repeat {r + 1}/{n_repeats}...")

    # Compute statistical moments and 95% confidence intervals across seeds
    stat_links = compute_ensemble_statistics(all_g_links, axis=1)
    stat_springs = compute_ensemble_statistics(all_g_springs, axis=1)
    stat_extra_l = compute_ensemble_statistics(all_extra_links, axis=1)
    stat_extra_s = compute_ensemble_statistics(all_extra_springs, axis=1)

    return {
        "density_factors": np.array(density_factors),
        "n_repeats": np.array(n_repeats),
        "nominal_n_ent": np.array(nominal_n_ent),
        "vol": np.array(vol),
        "g_links_mean": stat_links["mean"],
        "g_links_std": stat_links["std"],
        "g_links_sem": stat_links["sem"],
        "g_links_ci": stat_links["ci"],
        "g_links_all": all_g_links,
        "g_springs_mean": stat_springs["mean"],
        "g_springs_std": stat_springs["std"],
        "g_springs_sem": stat_springs["sem"],
        "g_springs_ci": stat_springs["ci"],
        "g_springs_all": all_g_springs,
        "extra_atoms_links_mean": stat_extra_l["mean"],
        "extra_atoms_links_std": stat_extra_l["std"],
        "extra_atoms_springs_mean": stat_extra_s["mean"],
        "extra_atoms_springs_std": stat_extra_s["std"],
    }


# --- Slipping vs. Fixed Links ---


def _simulate_slipping_vs_fixed_single_repeat(
    r: int,
    seed: int,
    box_l: float,
    b: float,
    n_strands: int,
    strand_len: int,
    conversion: float,
    density_factors: Sequence[float],
    nominal_n_ent: int,
    cutoff: float,
    r02: float,
    vol: float,
    gamma_conv: float,
) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Execute a single Monte Carlo realization comparing fixed vs. slipping links."""
    gen = MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(seed)
    gen.set_bead_distance(b)
    gen.set_mean_squared_bead_distance(b**2)
    gen.add_strands(n_strands, [strand_len] * n_strands)
    gen.add_crosslinkers(int(0.5 * n_strands), 4)
    gen.link_strands_to_conversion(conversion)

    n_factors = len(density_factors)
    g_fixed = np.zeros(n_factors)
    g_slip = np.zeros(n_factors)
    iter_fixed = np.zeros(n_factors)
    iter_slip = np.zeros(n_factors)

    for i, factor in enumerate(density_factors):
        n_target = int(factor * nominal_n_ent)

        if n_target == 0:
            # Bare unentangled crosslinker network
            fb = gen.get_mehp_simulator()
            fb.run(
                pt.MEHPSlipLinkSolverConfig(
                    max_nr_of_steps=5000, slip_mode=pt.SlipMode.FIXED
                )
            )
            g_pa = (gamma_conv * np.sum(fb.get_gamma_factors(r02)) / vol) * 1e6
            it = fb.get_nr_of_iterations()
            g_fixed[i] = g_pa
            g_slip[i] = g_pa
            iter_fixed[i] = it
            iter_slip[i] = it
        else:
            # Sample entanglements on demand
            es = gen.sample_entanglements(
                n_target,
                cutoff,
                0.0,
                0,
                3.0,
                f"ent_{seed}_{factor}",
                filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
            )

            # 1. Fixed Entanglement Links (slip_mode=SlipMode.FIXED)
            fb_fixed = gen.get_mehp_simulator(es)
            fb_fixed.run(
                pt.MEHPSlipLinkSolverConfig(
                    max_nr_of_steps=5000, slip_mode=pt.SlipMode.FIXED
                )
            )
            g_fixed[i] = (
                gamma_conv * np.sum(fb_fixed.get_gamma_factors(r02)) / vol
            ) * 1e6
            iter_fixed[i] = fb_fixed.get_nr_of_iterations()

            # 2. Dynamic Slipping Links (slip_mode=SlipMode.SLIDING)
            fb_slip = gen.get_mehp_simulator(es)
            fb_slip.run(
                pt.MEHPSlipLinkSolverConfig(
                    max_nr_of_steps=5000, slip_mode=pt.SlipMode.SLIDING
                )
            )
            g_slip[i] = (
                gamma_conv * np.sum(fb_slip.get_gamma_factors(r02)) / vol
            ) * 1e6
            iter_slip[i] = fb_slip.get_nr_of_iterations()

    return r, g_fixed, g_slip, iter_fixed, iter_slip


def run_slipping_vs_fixed_simulation(
    density_factors: Sequence[float] = (0.0, 0.5, 1.0, 1.5, 2.0),
    n_strands: int = 10000,
    strand_len: int = 50,
    conversion: float = 0.95,
    n_repeats: int = 10,
    base_seed: int = 42,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Evaluate shear modulus and iteration count comparing slipping vs. fixed links."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = params.get("<b>").to(params.get("distance_units")).magnitude
    bead_density = params.get_bead_density()
    cutoff = params.get_sampling_cutoff()
    nominal_ent_density = params.get_entanglement_density()
    r02 = params.get("R02").to(params.get("distance_units") ** 2).magnitude
    kbt = params.get("T") * params.get("kb")
    gamma_conv = ((kbt / ((params.get("distance_units")) ** 3)).to("MPa")).magnitude

    n_atoms = n_strands * strand_len + int(0.5 * n_strands)
    vol = n_atoms / bead_density
    box_l = vol ** (1.0 / 3.0)
    nominal_n_ent = int(nominal_ent_density * vol)

    all_g_fixed = np.zeros((len(density_factors), n_repeats))
    all_g_slip = np.zeros((len(density_factors), n_repeats))
    all_iter_fixed = np.zeros((len(density_factors), n_repeats))
    all_iter_slip = np.zeros((len(density_factors), n_repeats))

    tasks = [
        (
            r,
            base_seed + 101 * r,
            box_l,
            b,
            n_strands,
            strand_len,
            conversion,
            density_factors,
            nominal_n_ent,
            cutoff,
            r02,
            vol,
            gamma_conv,
        )
        for r in range(n_repeats)
    ]

    workers = min(resolve_max_workers(max_workers), n_repeats)
    configure_process_threading(workers)

    if workers > 1 and n_repeats > 1:
        print(
            f"Executing {n_repeats} repeats in parallel across {workers} worker processes..."
        )
        completed = 0
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = [
                executor.submit(_simulate_slipping_vs_fixed_single_repeat, *task)
                for task in tasks
            ]
            for fut in as_completed(futures):
                r, g_fixed, g_slip, iter_fixed, iter_slip = fut.result()
                all_g_fixed[:, r] = g_fixed
                all_g_slip[:, r] = g_slip
                all_iter_fixed[:, r] = iter_fixed
                all_iter_slip[:, r] = iter_slip
                completed += 1
                print(
                    f"  [Slipping vs Fixed] Completed repeat {completed}/{
                        n_repeats
                    } (repeat index {r + 1})..."
                )
    else:
        for task in tasks:
            r, g_fixed, g_slip, iter_fixed, iter_slip = (
                _simulate_slipping_vs_fixed_single_repeat(*task)
            )
            all_g_fixed[:, r] = g_fixed
            all_g_slip[:, r] = g_slip
            all_iter_fixed[:, r] = iter_fixed
            all_iter_slip[:, r] = iter_slip
            print(f"  [Slipping vs Fixed] Completed repeat {r + 1}/{n_repeats}...")

    stat_g_fixed = compute_ensemble_statistics(all_g_fixed, axis=1)
    stat_g_slip = compute_ensemble_statistics(all_g_slip, axis=1)
    stat_it_fixed = compute_ensemble_statistics(all_iter_fixed, axis=1)
    stat_it_slip = compute_ensemble_statistics(all_iter_slip, axis=1)

    return {
        "density_factors": np.array(density_factors),
        "n_repeats": np.array(n_repeats),
        "nominal_n_ent": np.array(nominal_n_ent),
        "vol": np.array(vol),
        "g_fixed_mean": stat_g_fixed["mean"],
        "g_fixed_std": stat_g_fixed["std"],
        "g_fixed_sem": stat_g_fixed["sem"],
        "g_fixed_ci": stat_g_fixed["ci"],
        "g_fixed_all": all_g_fixed,
        "g_slip_mean": stat_g_slip["mean"],
        "g_slip_std": stat_g_slip["std"],
        "g_slip_sem": stat_g_slip["sem"],
        "g_slip_ci": stat_g_slip["ci"],
        "g_slip_all": all_g_slip,
        "iter_fixed_mean": stat_it_fixed["mean"],
        "iter_fixed_std": stat_it_fixed["std"],
        "iter_fixed_sem": stat_it_fixed["sem"],
        "iter_fixed_ci": stat_it_fixed["ci"],
        "iter_slip_mean": stat_it_slip["mean"],
        "iter_slip_std": stat_it_slip["std"],
        "iter_slip_sem": stat_it_slip["sem"],
        "iter_slip_ci": stat_it_slip["ci"],
    }


# --- Polydispersity Study ---


def _simulate_polydispersity_single_repeat(
    r: int,
    seed: int,
    pdi_list: Sequence[float],
    n_strands: int,
    avg_len: int,
    conversion: float,
    b: float,
    cutoff: float,
    ent_density: float,
    r02: float,
    bead_density: float,
    gamma_conv: float,
) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Execute a single Monte Carlo realization across the specified polydispersity range."""
    rng = random.Random(seed)
    n_pdi = len(pdi_list)
    g_total = np.zeros(n_pdi)
    g_phantom = np.zeros(n_pdi)
    w_active = np.zeros(n_pdi)
    w_dangling = np.zeros(n_pdi)

    for i, d in enumerate(pdi_list):
        s = math.sqrt(math.log(d)) if d > 1.0 else 0.0
        m = math.log(avg_len) - 0.5 * s**2 if d > 1.0 else math.log(avg_len)

        # Sample strand lengths from log-normal distribution (or monodisperse
        # for d == 1.0)
        strand_lengths = [
            max(1, int(round(rng.lognormvariate(m, s)))) if d > 1.0 else avg_len
            for _ in range(n_strands)
        ]

        n_atoms_total = sum(strand_lengths) + 0.5 * n_strands
        vol = n_atoms_total / bead_density

        gen = MCUniverseGenerator(
            vol ** (1.0 / 3.0),
            vol ** (1.0 / 3.0),
            vol ** (1.0 / 3.0),
        )
        gen.set_seed(seed + 31 * i)
        gen.config_nr_of_mc_steps(0)
        gen.set_bead_distance(b)
        gen.add_strands(nr_of_strands=n_strands, strand_lengths=strand_lengths)
        gen.add_crosslinkers(int(0.5 * n_strands), 4)
        gen.link_strands_to_conversion(crosslinker_conversion=conversion)

        # 1. Total Entangled Shear Modulus (with static entanglement links)
        fb_ent = gen.get_mehp_simulator(
            nr_of_entanglements_to_sample=int(ent_density * vol),
            upper_cutoff=cutoff,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            seed=str(seed + 31 * i),
        )
        fb_ent.run()
        g_total[i] = gamma_conv * np.sum(fb_ent.get_gamma_factors(r02)) / vol
        w_active[i] = fb_ent.get_active_weight_fraction()
        w_dangling[i] = fb_ent.get_dangling_weight_fraction()

        # 2. Phantom Shear Modulus (unentangled network baseline)
        fb_ph = gen.get_mehp_simulator()
        fb_ph.run()
        g_phantom[i] = gamma_conv * np.sum(fb_ph.get_gamma_factors(r02)) / vol

    return r, g_total, g_phantom, w_active, w_dangling


def run_polydispersity_simulation(
    pdi_values: Sequence[float] = (1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8),
    n_strands: int = 10000,
    avg_len: int = 50,
    conversion: float = 0.925,
    n_repeats: int = 8,
    base_seed: int = 42,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Execute multi-seed ensemble simulations of polymer network polydispersity."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = params.get("<b>").to(params.get("distance_units")).magnitude
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()
    r02 = params.get("R02").to(params.get("distance_units") ** 2).magnitude
    bead_density = params.get_bead_density()
    gamma_conv = params.get_gamma_conversion_factor().to("MPa").magnitude

    n_pdi = len(pdi_values)
    all_g_total = np.zeros((n_pdi, n_repeats))
    all_g_phantom = np.zeros((n_pdi, n_repeats))
    all_w_active = np.zeros((n_pdi, n_repeats))
    all_w_dangling = np.zeros((n_pdi, n_repeats))

    tasks = [
        (
            r,
            base_seed + 101 * r,
            pdi_values,
            n_strands,
            avg_len,
            conversion,
            b,
            cutoff,
            ent_density,
            r02,
            bead_density,
            gamma_conv,
        )
        for r in range(n_repeats)
    ]

    workers = min(resolve_max_workers(max_workers), n_repeats)
    configure_process_threading(workers)

    if workers > 1 and n_repeats > 1:
        print(
            f"Executing {n_repeats} polydispersity repeats in parallel across {workers} workers..."
        )
        completed = 0
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = [
                executor.submit(_simulate_polydispersity_single_repeat, *task)
                for task in tasks
            ]
            for fut in as_completed(futures):
                r, g_tot, g_ph, w_act, w_dang = fut.result()
                all_g_total[:, r] = g_tot
                all_g_phantom[:, r] = g_ph
                all_w_active[:, r] = w_act
                all_w_dangling[:, r] = w_dang
                completed += 1
                print(
                    f"  [Polydispersity Study] Completed repeat {completed}/{
                        n_repeats
                    } (repeat index {r + 1})..."
                )
    else:
        for task in tasks:
            r, g_tot, g_ph, w_act, w_dang = _simulate_polydispersity_single_repeat(
                *task
            )
            all_g_total[:, r] = g_tot
            all_g_phantom[:, r] = g_ph
            all_w_active[:, r] = w_act
            all_w_dangling[:, r] = w_dang
            print(f"  [Polydispersity Study] Completed repeat {r + 1}/{n_repeats}...")

    all_g_ent = all_g_total - all_g_phantom

    stat_g_total = compute_ensemble_statistics(all_g_total, axis=1)
    stat_g_phantom = compute_ensemble_statistics(all_g_phantom, axis=1)
    stat_g_ent = compute_ensemble_statistics(all_g_ent, axis=1)
    stat_w_active = compute_ensemble_statistics(all_w_active, axis=1)
    stat_w_dangling = compute_ensemble_statistics(all_w_dangling, axis=1)

    return {
        "pdi_values": np.array(pdi_values),
        "n_repeats": np.array(n_repeats),
        "n_strands": np.array(n_strands),
        "avg_len": np.array(avg_len),
        "conversion": np.array(conversion),
        "g_total_mean": stat_g_total["mean"],
        "g_total_std": stat_g_total["std"],
        "g_total_sem": stat_g_total["sem"],
        "g_total_ci": stat_g_total["ci"],
        "g_total_all": all_g_total,
        "g_phantom_mean": stat_g_phantom["mean"],
        "g_phantom_std": stat_g_phantom["std"],
        "g_phantom_sem": stat_g_phantom["sem"],
        "g_phantom_ci": stat_g_phantom["ci"],
        "g_phantom_all": all_g_phantom,
        "g_ent_mean": stat_g_ent["mean"],
        "g_ent_std": stat_g_ent["std"],
        "g_ent_sem": stat_g_ent["sem"],
        "g_ent_ci": stat_g_ent["ci"],
        "g_ent_all": all_g_ent,
        "w_active_mean": stat_w_active["mean"],
        "w_active_std": stat_w_active["std"],
        "w_active_sem": stat_w_active["sem"],
        "w_active_ci": stat_w_active["ci"],
        "w_active_all": all_w_active,
        "w_dangling_mean": stat_w_dangling["mean"],
        "w_dangling_std": stat_w_dangling["std"],
        "w_dangling_sem": stat_w_dangling["sem"],
        "w_dangling_ci": stat_w_dangling["ci"],
        "w_dangling_all": all_w_dangling,
    }


# --- Sampling Cutoff Effect ---


def _simulate_sampling_cutoff_single_repeat(
    repeat_idx: int,
    seed: int,
    cutoffs: Sequence[float],
    n_strands: int,
    b: float,
    rho: float,
    r02_mag: float,
    gamma_conv: float,
    ent_density: float,
) -> Dict[str, Any]:
    """Execute a single Monte Carlo realization across sampling cutoffs and network conditions."""
    n_cutoffs = len(cutoffs)

    # -------------------------------------------------------------
    # 1. Standard Network: N=50, r=1.0, p=0.95 (Links vs. Springs)
    # -------------------------------------------------------------
    N_std = 50
    vol_std = (n_strands * N_std + int(0.5 * n_strands)) / rho
    box_l_std = vol_std ** (1.0 / 3.0)
    n_ent_std = int(ent_density * vol_std)

    gen_std = pt.MCUniverseGenerator(box_l_std, box_l_std, box_l_std)
    gen_std.set_seed(seed)
    gen_std.config_nr_of_mc_steps(0)
    gen_std.set_bead_distance(b)
    gen_std.set_mean_squared_bead_distance(b**2)
    gen_std.add_strands(n_strands, [N_std] * n_strands)
    gen_std.add_crosslinkers(int(0.5 * n_strands), 4)
    gen_std.link_strands_to_conversion(0.95)

    g_links = np.zeros(n_cutoffs)
    g_springs = np.zeros(n_cutoffs)
    res_links = np.zeros(n_cutoffs)
    res_springs = np.zeros(n_cutoffs)
    iter_links = np.zeros(n_cutoffs)
    iter_springs = np.zeros(n_cutoffs)

    for i, c in enumerate(cutoffs):
        # Rigid Entanglement Links
        fb_l = gen_std.get_mehp_simulator(
            nr_of_entanglements_to_sample=n_ent_std,
            upper_cutoff=c,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            seed=f"link_{seed}_{c}",
        )
        fb_l.validate_network()
        res_links[i] = fb_l.get_residual_norm()
        fb_l.run()
        g_links[i] = gamma_conv * np.sum(fb_l.get_gamma_factors(r02_mag)) / vol_std
        iter_links[i] = fb_l.get_nr_of_iterations()

        # Compliant Entanglement Springs
        fb_s = gen_std.get_mehp_simulator(
            nr_of_entanglements_to_sample=n_ent_std,
            upper_cutoff=c,
            entanglement_mode=pt.EntanglementMode.SPRING,
            seed=f"spring_{seed}_{c}",
        )
        fb_s.validate_network()
        res_springs[i] = fb_s.get_residual_norm()
        fb_s.run()
        g_springs[i] = gamma_conv * np.sum(fb_s.get_gamma_factors(r02_mag)) / vol_std
        iter_springs[i] = fb_s.get_nr_of_iterations()

    # -------------------------------------------------------------
    # 2. Strand Length Variation: N=30 vs N=70 (Links, r=1.0, p=0.95)
    # -------------------------------------------------------------
    # N = 30 (Crosslink-dominated)
    vol_n30 = (n_strands * 30 + int(0.5 * n_strands)) / rho
    box_l_n30 = vol_n30 ** (1.0 / 3.0)
    n_ent_n30 = int(ent_density * vol_n30)
    gen_n30 = pt.MCUniverseGenerator(box_l_n30, box_l_n30, box_l_n30)
    gen_n30.set_seed(seed + 100)
    gen_n30.config_nr_of_mc_steps(0)
    gen_n30.set_bead_distance(b)
    gen_n30.set_mean_squared_bead_distance(b**2)
    gen_n30.add_strands(n_strands, [30] * n_strands)
    gen_n30.add_crosslinkers(int(0.5 * n_strands), 4)
    gen_n30.link_strands_to_conversion(0.95)

    g_n30 = np.zeros(n_cutoffs)
    for i, c in enumerate(cutoffs):
        fb = gen_n30.get_mehp_simulator(
            nr_of_entanglements_to_sample=n_ent_n30,
            upper_cutoff=c,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            seed=f"n30_{seed}_{c}",
        )
        fb.run()
        g_n30[i] = gamma_conv * np.sum(fb.get_gamma_factors(r02_mag)) / vol_n30

    # N = 70 (Entanglement-dominated)
    vol_n70 = (n_strands * 70 + int(0.5 * n_strands)) / rho
    box_l_n70 = vol_n70 ** (1.0 / 3.0)
    n_ent_n70 = int(ent_density * vol_n70)
    gen_n70 = pt.MCUniverseGenerator(box_l_n70, box_l_n70, box_l_n70)
    gen_n70.set_seed(seed + 200)
    gen_n70.config_nr_of_mc_steps(0)
    gen_n70.set_bead_distance(b)
    gen_n70.set_mean_squared_bead_distance(b**2)
    gen_n70.add_strands(n_strands, [70] * n_strands)
    gen_n70.add_crosslinkers(int(0.5 * n_strands), 4)
    gen_n70.link_strands_to_conversion(0.95)

    g_n70 = np.zeros(n_cutoffs)
    for i, c in enumerate(cutoffs):
        fb = gen_n70.get_mehp_simulator(
            nr_of_entanglements_to_sample=n_ent_n70,
            upper_cutoff=c,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            seed=f"n70_{seed}_{c}",
        )
        fb.run()
        g_n70[i] = gamma_conv * np.sum(fb.get_gamma_factors(r02_mag)) / vol_n70

    # -------------------------------------------------------------
    # 3. Conversion & Stoichiometry: p=0.82 vs 0.95, r=0.75 vs 1.00
    # -------------------------------------------------------------
    # Incomplete Conversion: p = 0.82 (r = 1.0, N = 50)
    gen_p82 = pt.MCUniverseGenerator(box_l_std, box_l_std, box_l_std)
    gen_p82.set_seed(seed + 300)
    gen_p82.config_nr_of_mc_steps(0)
    gen_p82.set_bead_distance(b)
    gen_p82.set_mean_squared_bead_distance(b**2)
    gen_p82.add_strands(n_strands, [N_std] * n_strands)
    gen_p82.add_crosslinkers(int(0.5 * n_strands), 4)
    gen_p82.link_strands_to_conversion(0.82)

    g_p82 = np.zeros(n_cutoffs)
    for i, c in enumerate(cutoffs):
        fb = gen_p82.get_mehp_simulator(
            nr_of_entanglements_to_sample=n_ent_std,
            upper_cutoff=c,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            seed=f"p82_{seed}_{c}",
        )
        fb.run()
        g_p82[i] = gamma_conv * np.sum(fb.get_gamma_factors(r02_mag)) / vol_std

    # Stoichiometric Imbalance: r = 0.75 (xlinks = round(0.75 * 0.5 *
    # n_strands) = 112, p = 0.90)
    n_xlinks_r75 = int(0.75 * 0.5 * n_strands)
    vol_r75 = (n_strands * N_std + n_xlinks_r75) / rho
    box_l_r75 = vol_r75 ** (1.0 / 3.0)
    n_ent_r75 = int(ent_density * vol_r75)
    gen_r75 = pt.MCUniverseGenerator(box_l_r75, box_l_r75, box_l_r75)
    gen_r75.set_seed(seed + 400)
    gen_r75.config_nr_of_mc_steps(0)
    gen_r75.set_bead_distance(b)
    gen_r75.set_mean_squared_bead_distance(b**2)
    gen_r75.add_strands(n_strands, [N_std] * n_strands)
    gen_r75.add_crosslinkers(n_xlinks_r75, 4)
    gen_r75.link_strands_to_conversion(0.90)

    g_r75 = np.zeros(n_cutoffs)
    for i, c in enumerate(cutoffs):
        fb = gen_r75.get_mehp_simulator(
            nr_of_entanglements_to_sample=n_ent_r75,
            upper_cutoff=c,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            seed=f"r75_{seed}_{c}",
        )
        fb.run()
        g_r75[i] = gamma_conv * np.sum(fb.get_gamma_factors(r02_mag)) / vol_r75

    # -------------------------------------------------------------
    # 4. Pair Separation Distances for c = 1.2, 2.0, 3.0
    # -------------------------------------------------------------
    sampled_dists = {}
    for c_dist in [1.2, 2.0, 3.0]:
        ent = gen_std.sample_entanglements(
            nr_of_entanglements_to_sample=n_ent_std,
            upper_cutoff=c_dist,
            seed=f"dist_{seed}_{c_dist}",
        )
        pos1 = np.array(ent.position_1).reshape(-1, 3)
        pos2 = np.array(ent.position_2).reshape(-1, 3)
        diff = pos2 - pos1
        diff = diff - np.round(diff / box_l_std) * box_l_std
        dists = np.linalg.norm(diff, axis=1)
        sampled_dists[f"dists_{c_dist:.1f}"] = np.array(dists, dtype=float)

    return {
        "repeat_idx": repeat_idx,
        "g_links": g_links,
        "g_springs": g_springs,
        "res_links": res_links,
        "res_springs": res_springs,
        "iter_links": iter_links,
        "iter_springs": iter_springs,
        "g_n30": g_n30,
        "g_n70": g_n70,
        "g_p82": g_p82,
        "g_r75": g_r75,
        "sampled_dists": sampled_dists,
    }


def run_sampling_cutoff_simulation(
    n_repeats: int = 8,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Execute multi-seed ensemble simulation across sampling cutoffs and system parameters."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = float(params.get("<b>").to(params.get("distance_units")).magnitude)
    rho = float(params.get_bead_density())
    r02_mag = float(params.get("R02").to(params.get("distance_units") ** 2).magnitude)
    kbt = params.get("T") * params.get("kb")
    gamma_conv = float((kbt / (params.get("distance_units") ** 3)).to("MPa").magnitude)
    ent_density = float(params.get_entanglement_density())
    default_cutoff = float(params.get_sampling_cutoff())

    cutoffs = np.array([0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5])
    n_cutoffs = len(cutoffs)
    n_strands = 300

    tasks = [
        (
            r,
            1001 + r * 37,
            cutoffs,
            n_strands,
            b,
            rho,
            r02_mag,
            gamma_conv,
            ent_density,
        )
        for r in range(n_repeats)
    ]

    all_g_links = np.zeros((n_cutoffs, n_repeats))
    all_g_springs = np.zeros((n_cutoffs, n_repeats))
    all_res_links = np.zeros((n_cutoffs, n_repeats))
    all_res_springs = np.zeros((n_cutoffs, n_repeats))
    all_iter_links = np.zeros((n_cutoffs, n_repeats))
    all_iter_springs = np.zeros((n_cutoffs, n_repeats))
    all_g_n30 = np.zeros((n_cutoffs, n_repeats))
    all_g_n70 = np.zeros((n_cutoffs, n_repeats))
    all_g_p82 = np.zeros((n_cutoffs, n_repeats))
    all_g_r75 = np.zeros((n_cutoffs, n_repeats))

    collected_dists = {"1.2": [], "2.0": [], "3.0": []}

    workers = min(resolve_max_workers(max_workers), n_repeats)
    configure_process_threading(workers)
    print(f"Executing {n_repeats} repeats in parallel across {workers} workers...")

    if workers > 1 and n_repeats > 1:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = [
                executor.submit(_simulate_sampling_cutoff_single_repeat, *task)
                for task in tasks
            ]
            completed = 0
            for fut in as_completed(futures):
                res = fut.result()
                r = res["repeat_idx"]
                all_g_links[:, r] = res["g_links"]
                all_g_springs[:, r] = res["g_springs"]
                all_res_links[:, r] = res["res_links"]
                all_res_springs[:, r] = res["res_springs"]
                all_iter_links[:, r] = res["iter_links"]
                all_iter_springs[:, r] = res["iter_springs"]
                all_g_n30[:, r] = res["g_n30"]
                all_g_n70[:, r] = res["g_n70"]
                all_g_p82[:, r] = res["g_p82"]
                all_g_r75[:, r] = res["g_r75"]
                for c_k in ["1.2", "2.0", "3.0"]:
                    collected_dists[c_k].append(res["sampled_dists"][f"dists_{c_k}"])
                completed += 1
                print(
                    f"  [Cutoff Study] Completed repeat {completed}/{
                        n_repeats
                    } (repeat {r + 1})..."
                )
    else:
        for task in tasks:
            res = _simulate_sampling_cutoff_single_repeat(*task)
            r = res["repeat_idx"]
            all_g_links[:, r] = res["g_links"]
            all_g_springs[:, r] = res["g_springs"]
            all_res_links[:, r] = res["res_links"]
            all_res_springs[:, r] = res["res_springs"]
            all_iter_links[:, r] = res["iter_links"]
            all_iter_springs[:, r] = res["iter_springs"]
            all_g_n30[:, r] = res["g_n30"]
            all_g_n70[:, r] = res["g_n70"]
            all_g_p82[:, r] = res["g_p82"]
            all_g_r75[:, r] = res["g_r75"]
            for c_k in ["1.2", "2.0", "3.0"]:
                collected_dists[c_k].append(res["sampled_dists"][f"dists_{c_k}"])
            print(f"  [Cutoff Study] Completed repeat {r + 1}/{n_repeats}...")

    # Compute statistical moments across repeats
    stat_g_links = compute_ensemble_statistics(all_g_links, axis=1)
    stat_g_springs = compute_ensemble_statistics(all_g_springs, axis=1)
    stat_res_links = compute_ensemble_statistics(all_res_links, axis=1)
    stat_res_springs = compute_ensemble_statistics(all_res_springs, axis=1)
    stat_it_links = compute_ensemble_statistics(all_iter_links, axis=1)
    stat_it_springs = compute_ensemble_statistics(all_iter_springs, axis=1)
    stat_g_n30 = compute_ensemble_statistics(all_g_n30, axis=1)
    stat_g_n70 = compute_ensemble_statistics(all_g_n70, axis=1)
    stat_g_p82 = compute_ensemble_statistics(all_g_p82, axis=1)
    stat_g_r75 = compute_ensemble_statistics(all_g_r75, axis=1)

    # Aggregate distance histograms
    bin_edges = np.linspace(0.0, 3.5, 36)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    hist_12, _ = np.histogram(
        np.concatenate(collected_dists["1.2"]), bins=bin_edges, density=True
    )
    hist_20, _ = np.histogram(
        np.concatenate(collected_dists["2.0"]), bins=bin_edges, density=True
    )
    hist_30, _ = np.histogram(
        np.concatenate(collected_dists["3.0"]), bins=bin_edges, density=True
    )

    return {
        "cutoffs": cutoffs,
        "n_repeats": np.array(n_repeats),
        "b": np.array(b),
        "default_cutoff": np.array(default_cutoff),
        # Study 1: Links vs Springs
        "g_links_mean": stat_g_links["mean"],
        "g_links_ci": stat_g_links["ci"],
        "g_links_all": all_g_links,
        "g_springs_mean": stat_g_springs["mean"],
        "g_springs_ci": stat_g_springs["ci"],
        "g_springs_all": all_g_springs,
        "res_links_mean": stat_res_links["mean"],
        "res_links_ci": stat_res_links["ci"],
        "res_springs_mean": stat_res_springs["mean"],
        "res_springs_ci": stat_res_springs["ci"],
        "iter_links_mean": stat_it_links["mean"],
        "iter_links_ci": stat_it_links["ci"],
        "iter_springs_mean": stat_it_springs["mean"],
        "iter_springs_ci": stat_it_springs["ci"],
        # Study 2: N=30 vs N=70
        "g_n30_mean": stat_g_n30["mean"],
        "g_n30_ci": stat_g_n30["ci"],
        "g_n70_mean": stat_g_n70["mean"],
        "g_n70_ci": stat_g_n70["ci"],
        # Study 3: p=0.82 vs p=0.95 & r=0.75 vs r=1.00
        "g_p82_mean": stat_g_p82["mean"],
        "g_p82_ci": stat_g_p82["ci"],
        "g_r75_mean": stat_g_r75["mean"],
        "g_r75_ci": stat_g_r75["ci"],
        # Study 4: Distance distributions
        "bin_centers": bin_centers,
        "hist_12": hist_12,
        "hist_20": hist_20,
        "hist_30": hist_30,
    }


# --- Inactive Removal & Sampling ---


def _simulate_inactive_removal_single_repeat(
    r: int,
    seed: int,
    box_l: float,
    b: float,
    n_strands: int,
    strand_len: int,
    conversions: Sequence[float],
    n_ent: int,
    cutoff: float,
    r02: float,
    vol: float,
    gamma_conv: float,
    inactive_cutoff: float,
) -> Tuple[
    int,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Execute a single Monte Carlo realization across conversions and simplification modes."""
    n_conv = len(conversions)

    # Modulus arrays [kPa]
    g_unf_nosimp = np.zeros(n_conv)
    g_unf_inact = np.zeros(n_conv)
    g_unf_x2f = np.zeros(n_conv)
    g_filt_nosimp = np.zeros(n_conv)
    g_filt_x2f = np.zeros(n_conv)

    # Spring count arrays
    sp_unf_nosimp = np.zeros(n_conv)
    sp_unf_inact = np.zeros(n_conv)
    sp_unf_x2f = np.zeros(n_conv)
    sp_filt_nosimp = np.zeros(n_conv)
    sp_filt_x2f = np.zeros(n_conv)

    # Link count arrays
    links_unf_nosimp = np.zeros(n_conv)
    links_unf_inact = np.zeros(n_conv)
    links_unf_x2f = np.zeros(n_conv)
    links_filt_nosimp = np.zeros(n_conv)
    links_filt_x2f = np.zeros(n_conv)

    # Iteration count arrays
    iter_unf_nosimp = np.zeros(n_conv)
    iter_unf_inact = np.zeros(n_conv)
    iter_unf_x2f = np.zeros(n_conv)
    iter_filt_nosimp = np.zeros(n_conv)
    iter_filt_x2f = np.zeros(n_conv)

    for i, conv in enumerate(conversions):
        gen = pt.MCUniverseGenerator(box_l, box_l, box_l)
        gen.set_seed(seed + i * 1000)
        gen.set_bead_distance(b)
        gen.set_mean_squared_bead_distance(b**2)
        gen.add_strands(n_strands, [strand_len] * n_strands)
        gen.add_crosslinkers(int(0.5 * n_strands), 4)
        gen.link_strands_to_conversion(conv)

        # Entanglement sampling: filtered vs. unfiltered
        es_filt = gen.sample_entanglements(
            n_ent,
            cutoff,
            0.0,
            0,
            3.0,
            f"s_filt_{seed}_{conv}",
            filtering_mode=pt.EntanglementFilteringMode.OMIT_INACTIVE,
        )
        es_unf = gen.sample_entanglements(
            n_ent,
            cutoff,
            0.0,
            0,
            3.0,
            f"s_unf_{seed}_{conv}",
            filtering_mode=pt.EntanglementFilteringMode.NONE,
        )

        # 1. Unfiltered + No Simplification
        fb1 = gen.get_mehp_simulator(es_unf)
        fb1.run(
            pt.MEHPLinearSolverConfig(
                simplification_mode=pt.StructureSimplificationMode.NO_SIMPLIFICATION
            )
        )
        g_unf_nosimp[i] = (gamma_conv * np.sum(fb1.get_gamma_factors(r02)) / vol) * 1e6
        sp_unf_nosimp[i] = fb1.network.nr_of_springs
        links_unf_nosimp[i] = fb1.network.nr_of_links
        iter_unf_nosimp[i] = fb1.get_nr_of_iterations()

        # 2. Unfiltered + INACTIVE_ONLY
        fb2 = gen.get_mehp_simulator(es_unf)
        fb2.run(
            pt.MEHPLinearSolverConfig(
                simplification_mode=pt.StructureSimplificationMode.INACTIVE_ONLY,
                inactive_removal_cutoff=inactive_cutoff,
            )
        )
        g_unf_inact[i] = (gamma_conv * np.sum(fb2.get_gamma_factors(r02)) / vol) * 1e6
        sp_unf_inact[i] = fb2.network.nr_of_springs
        links_unf_inact[i] = fb2.network.nr_of_links
        iter_unf_inact[i] = fb2.get_nr_of_iterations()

        # 3. Unfiltered + INACTIVE_THEN_X2F
        fb3 = gen.get_mehp_simulator(es_unf)
        fb3.run(
            pt.MEHPLinearSolverConfig(
                simplification_mode=pt.StructureSimplificationMode.INACTIVE_THEN_X2F,
                inactive_removal_cutoff=inactive_cutoff,
            )
        )
        g_unf_x2f[i] = (gamma_conv * np.sum(fb3.get_gamma_factors(r02)) / vol) * 1e6
        sp_unf_x2f[i] = fb3.network.nr_of_springs
        links_unf_x2f[i] = fb3.network.nr_of_links
        iter_unf_x2f[i] = fb3.get_nr_of_iterations()

        # 4. Filtered + No Simplification
        fb4 = gen.get_mehp_simulator(es_filt)
        fb4.run(
            pt.MEHPLinearSolverConfig(
                simplification_mode=pt.StructureSimplificationMode.NO_SIMPLIFICATION
            )
        )
        g_filt_nosimp[i] = (gamma_conv * np.sum(fb4.get_gamma_factors(r02)) / vol) * 1e6
        sp_filt_nosimp[i] = fb4.network.nr_of_springs
        links_filt_nosimp[i] = fb4.network.nr_of_links
        iter_filt_nosimp[i] = fb4.get_nr_of_iterations()

        # 5. Filtered + INACTIVE_THEN_X2F
        fb5 = gen.get_mehp_simulator(es_filt)
        fb5.run(
            pt.MEHPLinearSolverConfig(
                simplification_mode=pt.StructureSimplificationMode.INACTIVE_THEN_X2F,
                inactive_removal_cutoff=inactive_cutoff,
            )
        )
        g_filt_x2f[i] = (gamma_conv * np.sum(fb5.get_gamma_factors(r02)) / vol) * 1e6
        sp_filt_x2f[i] = fb5.network.nr_of_springs
        links_filt_x2f[i] = fb5.network.nr_of_links
        iter_filt_x2f[i] = fb5.get_nr_of_iterations()

    return (
        r,
        g_unf_nosimp,
        g_unf_inact,
        g_unf_x2f,
        g_filt_nosimp,
        g_filt_x2f,
        sp_unf_nosimp,
        sp_unf_inact,
        sp_unf_x2f,
        sp_filt_nosimp,
        sp_filt_x2f,
        links_unf_nosimp,
        links_unf_inact,
        links_unf_x2f,
        links_filt_nosimp,
        links_filt_x2f,
        iter_unf_nosimp,
        iter_unf_inact,
        iter_unf_x2f,
        iter_filt_nosimp,
        iter_filt_x2f,
    )


def run_inactive_removal_and_sampling_simulation(
    conversions: Sequence[float] = (0.70, 0.75, 0.80, 0.85, 0.90, 0.95),
    n_strands: int = 2000,
    strand_len: int = 50,
    box_l: float = 50.0,
    n_repeats: int = 8,
    base_seed: int = 42,
    inactive_cutoff: float = 1e-4,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Execute multi-seed ensemble simulation comparing inactive strand removal and sampling."""
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = params.get("<b>").to(params.get("distance_units")).magnitude
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()
    r02 = params.get("R02").to(params.get("distance_units") ** 2).magnitude
    gamma_conv = (
        (params.get("T") * params.get("kb") / (params.get("distance_units") ** 3))
        .to("MPa")
        .magnitude
    )

    vol = box_l**3
    n_ent = int(ent_density * vol)
    effective_workers = resolve_max_workers(max_workers)

    print(
        f"Simulating Inactive Removal & Sampling study: {n_repeats} repeats, "
        f"{len(conversions)} conversions, N_strands={n_strands}, N_ent={n_ent} across {effective_workers} workers..."
    )

    n_conv = len(conversions)
    all_g_unf_nosimp = np.zeros((n_conv, n_repeats))
    all_g_unf_inact = np.zeros((n_conv, n_repeats))
    all_g_unf_x2f = np.zeros((n_conv, n_repeats))
    all_g_filt_nosimp = np.zeros((n_conv, n_repeats))
    all_g_filt_x2f = np.zeros((n_conv, n_repeats))

    all_sp_unf_nosimp = np.zeros((n_conv, n_repeats))
    all_sp_unf_inact = np.zeros((n_conv, n_repeats))
    all_sp_unf_x2f = np.zeros((n_conv, n_repeats))
    all_sp_filt_nosimp = np.zeros((n_conv, n_repeats))
    all_sp_filt_x2f = np.zeros((n_conv, n_repeats))

    all_links_unf_nosimp = np.zeros((n_conv, n_repeats))
    all_links_unf_inact = np.zeros((n_conv, n_repeats))
    all_links_unf_x2f = np.zeros((n_conv, n_repeats))
    all_links_filt_nosimp = np.zeros((n_conv, n_repeats))
    all_links_filt_x2f = np.zeros((n_conv, n_repeats))

    all_iter_unf_nosimp = np.zeros((n_conv, n_repeats))
    all_iter_unf_inact = np.zeros((n_conv, n_repeats))
    all_iter_unf_x2f = np.zeros((n_conv, n_repeats))
    all_iter_filt_nosimp = np.zeros((n_conv, n_repeats))
    all_iter_filt_x2f = np.zeros((n_conv, n_repeats))

    tasks = [
        (
            r,
            base_seed + r * 100,
            box_l,
            b,
            n_strands,
            strand_len,
            conversions,
            n_ent,
            cutoff,
            r02,
            vol,
            gamma_conv,
            inactive_cutoff,
        )
        for r in range(n_repeats)
    ]

    configure_process_threading(effective_workers)

    if effective_workers > 1:
        with ProcessPoolExecutor(
            max_workers=effective_workers,
            initializer=worker_process_initializer,
        ) as executor:
            futures = [
                executor.submit(_simulate_inactive_removal_single_repeat, *t)
                for t in tasks
            ]
            completed = 0
            for fut in as_completed(futures):
                (
                    r,
                    g1,
                    g2,
                    g3,
                    g4,
                    g5,
                    sp1,
                    sp2,
                    sp3,
                    sp4,
                    sp5,
                    lk1,
                    lk2,
                    lk3,
                    lk4,
                    lk5,
                    it1,
                    it2,
                    it3,
                    it4,
                    it5,
                ) = fut.result()
                all_g_unf_nosimp[:, r] = g1
                all_g_unf_inact[:, r] = g2
                all_g_unf_x2f[:, r] = g3
                all_g_filt_nosimp[:, r] = g4
                all_g_filt_x2f[:, r] = g5

                all_sp_unf_nosimp[:, r] = sp1
                all_sp_unf_inact[:, r] = sp2
                all_sp_unf_x2f[:, r] = sp3
                all_sp_filt_nosimp[:, r] = sp4
                all_sp_filt_x2f[:, r] = sp5

                all_links_unf_nosimp[:, r] = lk1
                all_links_unf_inact[:, r] = lk2
                all_links_unf_x2f[:, r] = lk3
                all_links_filt_nosimp[:, r] = lk4
                all_links_filt_x2f[:, r] = lk5

                all_iter_unf_nosimp[:, r] = it1
                all_iter_unf_inact[:, r] = it2
                all_iter_unf_x2f[:, r] = it3
                all_iter_filt_nosimp[:, r] = it4
                all_iter_filt_x2f[:, r] = it5

                completed += 1
                print(
                    f"  [Inactive Removal] Completed repeat {completed}/{
                        n_repeats
                    } (repeat index {r + 1})..."
                )
    else:
        for t in tasks:
            (
                r,
                g1,
                g2,
                g3,
                g4,
                g5,
                sp1,
                sp2,
                sp3,
                sp4,
                sp5,
                lk1,
                lk2,
                lk3,
                lk4,
                lk5,
                it1,
                it2,
                it3,
                it4,
                it5,
            ) = _simulate_inactive_removal_single_repeat(*t)
            all_g_unf_nosimp[:, r] = g1
            all_g_unf_inact[:, r] = g2
            all_g_unf_x2f[:, r] = g3
            all_g_filt_nosimp[:, r] = g4
            all_g_filt_x2f[:, r] = g5

            all_sp_unf_nosimp[:, r] = sp1
            all_sp_unf_inact[:, r] = sp2
            all_sp_unf_x2f[:, r] = sp3
            all_sp_filt_nosimp[:, r] = sp4
            all_sp_filt_x2f[:, r] = sp5

            all_links_unf_nosimp[:, r] = lk1
            all_links_unf_inact[:, r] = lk2
            all_links_unf_x2f[:, r] = lk3
            all_links_filt_nosimp[:, r] = lk4
            all_links_filt_x2f[:, r] = lk5

            all_iter_unf_nosimp[:, r] = it1
            all_iter_unf_inact[:, r] = it2
            all_iter_unf_x2f[:, r] = it3
            all_iter_filt_nosimp[:, r] = it4
            all_iter_filt_x2f[:, r] = it5
            print(f"  [Inactive Removal] Completed repeat {r + 1}/{n_repeats}...")

    # Compute ensemble statistics
    stat_g_unf_nosimp = compute_ensemble_statistics(all_g_unf_nosimp, axis=1)
    stat_g_unf_inact = compute_ensemble_statistics(all_g_unf_inact, axis=1)
    stat_g_unf_x2f = compute_ensemble_statistics(all_g_unf_x2f, axis=1)
    stat_g_filt_nosimp = compute_ensemble_statistics(all_g_filt_nosimp, axis=1)
    stat_g_filt_x2f = compute_ensemble_statistics(all_g_filt_x2f, axis=1)

    stat_sp_unf_nosimp = compute_ensemble_statistics(all_sp_unf_nosimp, axis=1)
    stat_sp_unf_inact = compute_ensemble_statistics(all_sp_unf_inact, axis=1)
    stat_sp_unf_x2f = compute_ensemble_statistics(all_sp_unf_x2f, axis=1)
    stat_sp_filt_nosimp = compute_ensemble_statistics(all_sp_filt_nosimp, axis=1)
    stat_sp_filt_x2f = compute_ensemble_statistics(all_sp_filt_x2f, axis=1)

    stat_links_unf_nosimp = compute_ensemble_statistics(all_links_unf_nosimp, axis=1)
    stat_links_unf_inact = compute_ensemble_statistics(all_links_unf_inact, axis=1)
    stat_links_unf_x2f = compute_ensemble_statistics(all_links_unf_x2f, axis=1)
    stat_links_filt_nosimp = compute_ensemble_statistics(all_links_filt_nosimp, axis=1)
    stat_links_filt_x2f = compute_ensemble_statistics(all_links_filt_x2f, axis=1)

    stat_iter_unf_nosimp = compute_ensemble_statistics(all_iter_unf_nosimp, axis=1)
    stat_iter_unf_inact = compute_ensemble_statistics(all_iter_unf_inact, axis=1)
    stat_iter_unf_x2f = compute_ensemble_statistics(all_iter_unf_x2f, axis=1)
    stat_iter_filt_nosimp = compute_ensemble_statistics(all_iter_filt_nosimp, axis=1)
    stat_iter_filt_x2f = compute_ensemble_statistics(all_iter_filt_x2f, axis=1)

    return {
        "conversions": np.array(conversions),
        "n_repeats": np.array(n_repeats),
        "n_strands": np.array(n_strands),
        "strand_len": np.array(strand_len),
        "box_l": np.array(box_l),
        "vol": np.array(vol),
        "n_ent": np.array(n_ent),
        # Modulus
        "g_unf_nosimp_mean": stat_g_unf_nosimp["mean"],
        "g_unf_nosimp_ci": stat_g_unf_nosimp["ci"],
        "g_unf_nosimp_all": all_g_unf_nosimp,
        "g_unf_inact_mean": stat_g_unf_inact["mean"],
        "g_unf_inact_ci": stat_g_unf_inact["ci"],
        "g_unf_inact_all": all_g_unf_inact,
        "g_unf_x2f_mean": stat_g_unf_x2f["mean"],
        "g_unf_x2f_ci": stat_g_unf_x2f["ci"],
        "g_unf_x2f_all": all_g_unf_x2f,
        "g_filt_nosimp_mean": stat_g_filt_nosimp["mean"],
        "g_filt_nosimp_ci": stat_g_filt_nosimp["ci"],
        "g_filt_nosimp_all": all_g_filt_nosimp,
        "g_filt_x2f_mean": stat_g_filt_x2f["mean"],
        "g_filt_x2f_ci": stat_g_filt_x2f["ci"],
        "g_filt_x2f_all": all_g_filt_x2f,
        # Springs
        "sp_unf_nosimp_mean": stat_sp_unf_nosimp["mean"],
        "sp_unf_inact_mean": stat_sp_unf_inact["mean"],
        "sp_unf_x2f_mean": stat_sp_unf_x2f["mean"],
        "sp_filt_nosimp_mean": stat_sp_filt_nosimp["mean"],
        "sp_filt_x2f_mean": stat_sp_filt_x2f["mean"],
        # Links
        "links_unf_nosimp_mean": stat_links_unf_nosimp["mean"],
        "links_unf_inact_mean": stat_links_unf_inact["mean"],
        "links_unf_x2f_mean": stat_links_unf_x2f["mean"],
        "links_filt_nosimp_mean": stat_links_filt_nosimp["mean"],
        "links_filt_x2f_mean": stat_links_filt_x2f["mean"],
        # Iterations
        "iter_unf_nosimp_mean": stat_iter_unf_nosimp["mean"],
        "iter_unf_inact_mean": stat_iter_unf_inact["mean"],
        "iter_unf_x2f_mean": stat_iter_unf_x2f["mean"],
        "iter_filt_nosimp_mean": stat_iter_filt_nosimp["mean"],
        "iter_filt_x2f_mean": stat_iter_filt_x2f["mean"],
    }


# --- Entanglement Sampling Comparison (Discrete vs. On-Demand) ---


_SCALING_WORKER_CODE = """
import sys, time, resource, json
import pylimer_tools_cpp as pt

mode = sys.argv[1]
n_strands = int(sys.argv[2])
strand_len = int(sys.argv[3])
seed = int(sys.argv[4])
bead_density = float(sys.argv[5])
cutoff = float(sys.argv[6])
ent_density = float(sys.argv[7])

n_atoms = n_strands * strand_len + int(0.5 * n_strands)
vol = n_atoms / bead_density
L = vol ** (1.0 / 3.0)
n_ent = int(ent_density * vol)

filt = pt.EntanglementFilteringMode.OMIT_INACTIVE
samp_seed = f"seed_{seed}"

# Phase 0: network definition (identical for all methods)
t0 = time.perf_counter()
gen = pt.MCUniverseGenerator(L, L, L)
gen.config_nr_of_mc_steps(0)
gen.set_seed(seed)
gen.set_mean_squared_bead_distance(1.0)
gen.add_strands(n_strands, [strand_len] * n_strands)
gen.add_crosslinkers(int(0.5 * n_strands), 4)
gen.link_strands_to_conversion(0.92)
t_setup = time.perf_counter()

# Phase 1: bead realization (discrete only) and entanglement sampling,
# all with the same filtering (OMIT_INACTIVE) and target count.
if mode == "discrete":
    u = gen.get_universe()
    t_universe = time.perf_counter()
    ents = pt.randomly_sample_entanglements(
        u, n_ent, cutoff, 0.0, 0, 3.0, samp_seed, filtering_mode=filt
    )
    n_sampled = len(ents.pairs_of_atoms)
    t_sample = time.perf_counter()
    fb = pt.MEHPSimulator(
        u,
        ents,
        crosslinker_type=2,
        entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
        is_2d=False,
    )
else:
    t_universe = t_setup
    es = gen.sample_entanglements(
        n_ent, cutoff, 0.0, 0, 3.0, samp_seed, filt,
    )
    n_sampled = len(es.alpha_1)
    t_sample = time.perf_counter()
    # Phase 2: force balance construction
    fb = gen.get_mehp_simulator(
        es, entanglement_mode=pt.EntanglementMode.MERGED_LINKS, crosslinker_type=2
    )
fb.validate_network()
t1 = time.perf_counter()

maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
if sys.platform == "darwin":
    rss_mb = maxrss / (1024.0 * 1024.0)
else:
    rss_mb = maxrss / 1024.0

# "time" excludes the common network definition; it covers universe + sampling + FB build.
print(json.dumps({
    "mode": mode, "n_strands": n_strands, "strand_len": strand_len,
    "n_atoms": n_atoms,
    "time": t1 - t_setup, "t_setup": t_setup - t0,
    "t_universe": t_universe - t_setup,
    "t_sample": t_sample - t_universe, "t_build": t1 - t_sample,
    "n_sampled": n_sampled, "rss_mb": rss_mb,
}))
"""


def _run_scaling_grid(
    configs: Sequence[Tuple[int, int]],
    seed: int = 42,
    n_repeats: int = 5,
    prefix: str = "scaling",
) -> Dict[str, np.ndarray]:
    """Run subprocess scaling benchmarks for a sequence of ``(n_strands, strand_len)`` pairs.

    Each repeat uses a distinct random seed (``seed + 101 * r``) so that median and
    interquartile-range (25th–75th percentile) statistics capture both timing jitter and
    network realization variance.
    """
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    bead_density = params.get_bead_density()
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()

    modes = ["on_demand", "discrete"]
    records: Dict[str, List[Dict[str, float]]] = {m: [] for m in modes}
    atoms: List[int] = []
    strands_out: List[int] = []
    strand_lens_out: List[int] = []

    for n_strands, strand_len in configs:
        runs: Dict[str, List[Dict[str, Any]]] = {m: [] for m in modes}
        # Interleave modes across repeats so that transient machine load
        # affects all methods alike; vary seed per repeat to capture
        # realization variance.
        for r in range(n_repeats):
            run_seed = seed + 101 * r
            for mode in modes:
                cmd = [
                    sys.executable,
                    "-c",
                    _SCALING_WORKER_CODE,
                    mode,
                    str(n_strands),
                    str(strand_len),
                    str(run_seed),
                    str(bead_density),
                    str(cutoff),
                    str(ent_density),
                ]
                res = subprocess.run(cmd, capture_output=True, text=True, check=True)
                lines = [
                    line
                    for line in res.stdout.strip().split("\n")
                    if line.startswith("{")
                ]
                runs[mode].append(json.loads(lines[-1]))
        atoms.append(int(runs[modes[0]][0]["n_atoms"]))
        strands_out.append(int(n_strands))
        strand_lens_out.append(int(strand_len))
        for mode in modes:
            med: Dict[str, float] = {
                k: float(np.median([r[k] for r in runs[mode]]))
                for k in (
                    "time",
                    "t_universe",
                    "t_sample",
                    "t_build",
                    "rss_mb",
                    "n_sampled",
                )
            }
            for k in ("time", "rss_mb"):
                vals = [r[k] for r in runs[mode]]
                med[f"{k}_q25"] = float(np.percentile(vals, 25))
                med[f"{k}_q75"] = float(np.percentile(vals, 75))
            records[mode].append(med)
            print(
                f"{mode:10s} N_s={n_strands:5d} N_b={strand_len:3d} "
                f"N_atoms={atoms[-1]:7d} | Time={med['time']:6.2f}s "
                f"[{med['time_q25']:.2f}, {med['time_q75']:.2f}] "
                f"(univ={med['t_universe']:5.2f}s, "
                f"sample={med['t_sample']:5.2f}s, build={med['t_build']:5.2f}s) | "
                f"RSS={med['rss_mb']:6.1f} MB | n_sampled={med['n_sampled']:.0f}"
            )

    def col(mode: str, key: str) -> np.ndarray:
        return np.array([r[key] for r in records[mode]])

    out: Dict[str, np.ndarray] = {
        f"{prefix}_atoms": np.array(atoms),
        f"{prefix}_strands": np.array(strands_out),
        f"{prefix}_strand_lens": np.array(strand_lens_out),
        f"{prefix}_n_repeats": np.array(n_repeats),
    }
    for tag in ("discrete", "on_demand"):
        out[f"{prefix}_time_{tag}"] = col(tag, "time")
        out[f"{prefix}_time_q25_{tag}"] = col(tag, "time_q25")
        out[f"{prefix}_time_q75_{tag}"] = col(tag, "time_q75")
        out[f"{prefix}_rss_{tag}"] = col(tag, "rss_mb")
        out[f"{prefix}_rss_q25_{tag}"] = col(tag, "rss_mb_q25")
        out[f"{prefix}_rss_q75_{tag}"] = col(tag, "rss_mb_q75")
        out[f"{prefix}_time_universe_{tag}"] = col(tag, "t_universe")
        out[f"{prefix}_time_sample_{tag}"] = col(tag, "t_sample")
        out[f"{prefix}_time_build_{tag}"] = col(tag, "t_build")
        out[f"{prefix}_n_sampled_{tag}"] = col(tag, "n_sampled")
    return out


def run_scaling_benchmark(
    sizes: Sequence[int] = (200, 500, 1000, 2000, 4000, 7000, 10000, 15000),
    strand_len: int = 50,
    seed: int = 42,
    n_repeats: int = 5,
) -> Dict[str, np.ndarray]:
    """Measure peak memory and runtime (median & IQR of ``n_repeats`` seeds) across system sizes.

    All methods use the same filtering mode and sampling target; the network definition
    (identical for all) is excluded from the timings. Sampling and force balance
    construction are timed separately.
    """
    print("Running Scaling Benchmark across system sizes...")
    return _run_scaling_grid(
        [(n, strand_len) for n in sizes],
        seed=seed,
        n_repeats=n_repeats,
        prefix="scaling",
    )


def run_strand_length_benchmark(
    strand_lens: Sequence[int] = (15, 25, 50, 100, 200, 500),
    target_atoms: int = 250000,
    seed: int = 42,
    n_repeats: int = 5,
) -> Dict[str, np.ndarray]:
    """Measure peak memory and runtime across strand lengths at fixed total atom count.

    Varying ``strand_len`` at fixed ``target_atoms`` separates :math:`O(N_{\\text{strands}})`
    effects (on-demand LRU cache capacity and strand neighbor list) from
    :math:`O(N_{\\text{atoms}})` effects (discrete ``Universe`` construction and RAM).
    """
    print(f"Running Strand-Length Benchmark at ~{target_atoms} atoms...")
    configs = [(max(10, round(target_atoms / (sl + 0.5))), sl) for sl in strand_lens]
    return _run_scaling_grid(
        configs,
        seed=seed,
        n_repeats=n_repeats,
        prefix="strlen",
    )


# %%
# Benchmark 2: Conformational Statistics Comparison
# --------------------------------------------------
# We sample entanglements and compare the distribution of continuous contour positions
# :math:`\alpha` and spatial Euclidean distances :math:`r` against explicit discrete bead sampling.


def run_distribution_comparison(
    n_strands: int = 10000,
    strand_len: int = 50,
    conversion: float = 0.95,
    seed: int = 123,
) -> Dict[str, np.ndarray]:
    """Sample entanglements and compare contour position alpha and spatial distance distributions."""
    print("Running Conformational Distribution Comparison...")
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = params.get("<b>").to(params.get("distance_units")).magnitude
    bead_density = params.get_bead_density()
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()

    n_atoms = n_strands * strand_len + int(0.5 * n_strands)
    vol = n_atoms / bead_density
    L = vol ** (1.0 / 3.0)
    n_ent = int(ent_density * vol)

    gen = pt.MCUniverseGenerator(L, L, L)
    gen.config_nr_of_mc_steps(0)
    gen.set_seed(seed)
    gen.set_bead_distance(b)
    gen.set_mean_squared_bead_distance(b**2)
    gen.add_strands(n_strands, [strand_len] * n_strands)
    gen.add_crosslinkers(int(0.5 * n_strands), 4)
    gen.link_strands_to_conversion(conversion)

    # 1. On-Demand sampling (Bounding box & on-demand bead regeneration)
    es_od = gen.sample_entanglements(
        n_ent,
        cutoff,
        0.0,
        0,
        3.0,
        "od_seed",
        pt.EntanglementFilteringMode.NONE,
    )
    alpha_od = np.concatenate([np.array(es_od.alpha_1), np.array(es_od.alpha_2)])

    pos1_od = np.array(es_od.position_1).reshape(-1, 3)
    pos2_od = np.array(es_od.position_2).reshape(-1, 3)
    diff_od = pos2_od - pos1_od
    diff_od = diff_od - np.round(diff_od / L) * L
    dist_od = np.linalg.norm(diff_od, axis=1)

    # 3. Discrete sampling on explicit Universe beads
    u = gen.get_universe()
    cu = gen.get_crosslinker_universe()
    atom_pairs = pt.randomly_sample_entanglements(
        u, n_ent, cutoff, 0.0, 0, 3.0, "disc_seed"
    )

    n_xlinks = len(cu.xlink_x)
    offsets = [n_xlinks]
    for n_beads in cu.beads_in_strand:
        offsets.append(offsets[-1] + n_beads)
    offsets_arr = np.array(offsets)

    alpha_disc = []
    dist_disc = []
    box = u.get_box()

    for a1, a2 in atom_pairs.pairs_of_atoms:
        idx1 = a1 - 1
        idx2 = a2 - 1
        if idx1 >= n_xlinks:
            s1 = np.searchsorted(offsets_arr, idx1, side="right") - 1
            b1 = idx1 - offsets_arr[s1]
            alpha_disc.append((b1 + 1.0) / (strand_len + 1.0))
        if idx2 >= n_xlinks:
            s2 = np.searchsorted(offsets_arr, idx2, side="right") - 1
            b2 = idx2 - offsets_arr[s2]
            alpha_disc.append((b2 + 1.0) / (strand_len + 1.0))

        atom1 = u.get_atom(a1)
        atom2 = u.get_atom(a2)
        dist_disc.append(atom1.distance_to(atom2, box))

    print(f"Sampled {len(dist_od)} on-demand pairs, {len(dist_disc)} discrete pairs.")
    return {
        "alpha_on_demand": np.array(alpha_od),
        "dist_on_demand": np.array(dist_od),
        "alpha_discrete": np.array(alpha_disc),
        "dist_discrete": np.array(dist_disc),
    }


# %%
# Benchmark 3: Crosslinker Conversion Sweep & Multi-Seed Ensembles
# ----------------------------------------------------------------
# We evaluate how defective networks behave across crosslinker conversion :math:`p \in [0.75, 1.00]`.
# Multi-seed repetitions provide empirical variance and Student's t 95%
# confidence intervals.


def _count_active_entanglements(fb: pt.MEHPSimulator) -> int:
    """Count the number of active load-bearing entanglements in the network."""
    net = fb.network
    if hasattr(net, "spring_is_entanglement") and net.spring_is_entanglement.any():
        return int(net.spring_is_entanglement.sum())
    count = 0
    for i in range(net.nr_of_links):
        if net.link_is_entanglement[i] and len(net.strand_indices_of_link[i]) >= 2:
            count += 1
    return count


def _simulate_conversion_sweep_single_repeat(
    r: int,
    seed: int,
    L: float,
    b: float,
    n_strands: int,
    strand_len: int,
    conversions: Sequence[float],
    n_ent: int,
    cutoff: float,
    r02: float,
    vol: float,
    gamma_conv: float,
) -> Tuple[
    int,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Execute a single Monte Carlo repeat for the conversion sweep benchmark."""
    n_p = len(conversions)
    g_d = np.zeros(n_p)
    w_d = np.zeros(n_p)
    n_d = np.zeros(n_p)
    g_odf = np.zeros(n_p)
    w_odf = np.zeros(n_p)
    n_odf = np.zeros(n_p)
    g_odu = np.zeros(n_p)
    w_odu = np.zeros(n_p)
    n_odu = np.zeros(n_p)

    for i, p in enumerate(conversions):
        gen = pt.MCUniverseGenerator(L, L, L)
        gen.config_nr_of_mc_steps(0)
        gen.set_seed(seed)
        gen.set_bead_distance(b)
        gen.set_mean_squared_bead_distance(b**2)
        gen.add_strands(n_strands, [strand_len] * n_strands)
        gen.add_crosslinkers(int(0.5 * n_strands), 4)
        gen.link_strands_to_conversion(p)

        # 1. Discrete Full Universe with INACTIVE_THEN_X2F simplification
        # (this constructor filters entanglements by default, i.e. it is
        # equivalent to the "Filtered" variants below)
        u = gen.get_universe()
        fb_disc = pt.MEHPSimulator.construct_with_random_entanglements(
            u,
            nr_of_entanglements_to_sample=n_ent,
            upper_sampling_cutoff=cutoff,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
        )
        fb_disc.run(
            pt.MEHPLinearSolverConfig(
                simplification_mode=pt.StructureSimplificationMode.INACTIVE_THEN_X2F
            )
        )
        g_d[i] = gamma_conv * np.sum(fb_disc.get_gamma_factors(r02)) / vol
        w_d[i] = fb_disc.get_active_weight_fraction()
        n_d[i] = _count_active_entanglements(fb_disc)

        # 2. On-Demand Filtered with INACTIVE_THEN_X2F simplification
        es_od_filt = gen.sample_entanglements(
            n_ent,
            cutoff,
            0.0,
            0,
            3.0,
            f"seed_od_f_{seed}",
            pt.EntanglementFilteringMode.OMIT_INACTIVE,
        )
        fb_od_filt = gen.get_mehp_simulator(
            es_od_filt,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            crosslinker_type=2,
        )
        fb_od_filt.run(
            pt.MEHPLinearSolverConfig(
                simplification_mode=pt.StructureSimplificationMode.INACTIVE_THEN_X2F
            )
        )
        g_odf[i] = gamma_conv * np.sum(fb_od_filt.get_gamma_factors(r02)) / vol
        w_odf[i] = fb_od_filt.get_active_weight_fraction()
        n_odf[i] = _count_active_entanglements(fb_od_filt)

        # 3. On-Demand Unfiltered with INACTIVE_THEN_X2F simplification
        es_od_unfilt = gen.sample_entanglements(
            n_ent,
            cutoff,
            0.0,
            0,
            3.0,
            f"seed_od_u_{seed}",
            pt.EntanglementFilteringMode.NONE,
        )
        fb_od_unfilt = gen.get_mehp_simulator(
            es_od_unfilt,
            entanglement_mode=pt.EntanglementMode.MERGED_LINKS,
            crosslinker_type=2,
        )
        fb_od_unfilt.run(
            pt.MEHPLinearSolverConfig(
                simplification_mode=pt.StructureSimplificationMode.INACTIVE_THEN_X2F
            )
        )
        g_odu[i] = gamma_conv * np.sum(fb_od_unfilt.get_gamma_factors(r02)) / vol
        w_odu[i] = fb_od_unfilt.get_active_weight_fraction()
        n_odu[i] = _count_active_entanglements(fb_od_unfilt)

    return (
        r,
        g_d,
        w_d,
        n_d,
        g_odf,
        w_odf,
        n_odf,
        g_odu,
        w_odu,
        n_odu,
    )


def run_conversion_sweep(
    conversions: Sequence[float] = (0.75, 0.80, 0.85, 0.90, 0.95, 1.00),
    n_strands: int = 10000,
    strand_len: int = 50,
    n_repeats: int = 10,
    base_seed: int = 777,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Compare modulus, active fractions, and active entanglement count across conversions."""
    print("Running Conversion Sweep across multi-seed ensembles...")
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = params.get("<b>").to(params.get("distance_units")).magnitude
    bead_density = params.get_bead_density()
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()
    r02 = params.get("R02").to(params.get("distance_units") ** 2).magnitude
    kbt = params.get("T") * params.get("kb")
    gamma_conv = ((kbt / ((params.get("distance_units")) ** 3)).to("MPa")).magnitude

    n_atoms = n_strands * strand_len + int(0.5 * n_strands)
    vol = n_atoms / bead_density
    L = vol ** (1.0 / 3.0)
    n_ent = int(ent_density * vol)

    n_p = len(conversions)
    g_d_all = np.zeros((n_p, n_repeats))
    w_d_all = np.zeros((n_p, n_repeats))
    n_d_all = np.zeros((n_p, n_repeats))
    g_odf_all = np.zeros((n_p, n_repeats))
    w_odf_all = np.zeros((n_p, n_repeats))
    n_odf_all = np.zeros((n_p, n_repeats))
    g_odu_all = np.zeros((n_p, n_repeats))
    w_odu_all = np.zeros((n_p, n_repeats))
    n_odu_all = np.zeros((n_p, n_repeats))

    tasks = [
        (
            r,
            base_seed + 101 * r,
            L,
            b,
            n_strands,
            strand_len,
            conversions,
            n_ent,
            cutoff,
            r02,
            vol,
            gamma_conv,
        )
        for r in range(n_repeats)
    ]

    workers = min(resolve_max_workers(max_workers), n_repeats)
    configure_process_threading(workers)

    if workers > 1 and n_repeats > 1:
        print(
            f"Executing {n_repeats} conversion sweep repeats in parallel across {workers} worker processes..."
        )
        completed = 0
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = [
                executor.submit(_simulate_conversion_sweep_single_repeat, *task)
                for task in tasks
            ]
            for fut in as_completed(futures):
                (
                    r,
                    g_d,
                    w_d,
                    n_d,
                    g_odf,
                    w_odf,
                    n_odf,
                    g_odu,
                    w_odu,
                    n_odu,
                ) = fut.result()
                g_d_all[:, r] = g_d
                w_d_all[:, r] = w_d
                n_d_all[:, r] = n_d
                g_odf_all[:, r] = g_odf
                w_odf_all[:, r] = w_odf
                n_odf_all[:, r] = n_odf
                g_odu_all[:, r] = g_odu
                w_odu_all[:, r] = w_odu
                n_odu_all[:, r] = n_odu
                completed += 1
                print(
                    f"  [Conversion Sweep] Completed repeat {completed}/{
                        n_repeats
                    } (repeat index {r + 1})..."
                )
    else:
        for task in tasks:
            (
                r,
                g_d,
                w_d,
                n_d,
                g_odf,
                w_odf,
                n_odf,
                g_odu,
                w_odu,
                n_odu,
            ) = _simulate_conversion_sweep_single_repeat(*task)
            g_d_all[:, r] = g_d
            w_d_all[:, r] = w_d
            n_d_all[:, r] = n_d
            g_odf_all[:, r] = g_odf
            w_odf_all[:, r] = w_odf
            n_odf_all[:, r] = n_odf
            g_odu_all[:, r] = g_odu
            w_odu_all[:, r] = w_odu
            n_odu_all[:, r] = n_odu
            print(f"  [Conversion Sweep] Completed repeat {r + 1}/{n_repeats}...")

    stat_gd = compute_ensemble_statistics(g_d_all, axis=1)
    stat_wd = compute_ensemble_statistics(w_d_all, axis=1)
    stat_nd = compute_ensemble_statistics(n_d_all, axis=1)
    stat_godf = compute_ensemble_statistics(g_odf_all, axis=1)
    stat_wodf = compute_ensemble_statistics(w_odf_all, axis=1)
    stat_nodf = compute_ensemble_statistics(n_odf_all, axis=1)
    stat_godu = compute_ensemble_statistics(g_odu_all, axis=1)
    stat_wodu = compute_ensemble_statistics(w_odu_all, axis=1)
    stat_nodu = compute_ensemble_statistics(n_odu_all, axis=1)

    return {
        "conversions": np.array(conversions),
        "target_n_ent": np.array(n_ent),
        "n_repeats": np.array(n_repeats),
        "g_discrete": stat_gd["mean"],
        "g_discrete_ci": stat_gd["ci"],
        "w_active_discrete": stat_wd["mean"],
        "w_active_discrete_ci": stat_wd["ci"],
        "n_ent_discrete": stat_nd["mean"],
        "n_ent_discrete_ci": stat_nd["ci"],
        "g_od_filtered": stat_godf["mean"],
        "g_od_filtered_ci": stat_godf["ci"],
        "w_active_od_filtered": stat_wodf["mean"],
        "w_active_od_filtered_ci": stat_wodf["ci"],
        "n_ent_od_filtered": stat_nodf["mean"],
        "n_ent_od_filtered_ci": stat_nodf["ci"],
        "g_od_unfiltered": stat_godu["mean"],
        "g_od_unfiltered_ci": stat_godu["ci"],
        "w_active_od_unfiltered": stat_wodu["mean"],
        "w_active_od_unfiltered_ci": stat_wodu["ci"],
        "n_ent_od_unfiltered": stat_nodu["mean"],
        "n_ent_od_unfiltered_ci": stat_nodu["ci"],
    }


def run_all_benchmarks(max_workers: Optional[int] = None) -> Dict[str, Any]:
    """Execute all benchmark suites and compile into a single dataset dictionary."""
    combined: Dict[str, Any] = {}
    combined.update(run_scaling_benchmark())
    combined.update(run_strand_length_benchmark())
    combined.update(run_distribution_comparison())
    combined.update(run_conversion_sweep(max_workers=max_workers))
    return combined


# --- Representative Volume Element (RVE) Study for Realistic PDMS ---


class _RveStudyKind(Enum):
    """Identifier for sub-studies in the PDMS Representative Volume Element sweep."""

    CHAIN_LENGTH = "N"
    CONVERSION = "p"
    STOICHIOMETRY = "r"
    SAMPLES = "samples"


def _run_single_pdms_rve_point(
    box_l: float,
    strand_len: int,
    conversion: float,
    stoichiometry_r: float,
    seed: int,
    b: float,
    b2: float,
    r02: float,
    cutoff: float,
    ent_density: float,
    bead_density: float,
    gamma_conv: float,
) -> Tuple[float, float, float, int, int]:
    """Simulate a single PDMS network realization at a given box size."""
    vol = box_l**3
    n_beads_target = bead_density * vol
    n_strands = max(
        10, int(round(n_beads_target / (strand_len + 0.5 * stoichiometry_r)))
    )
    n_xlinks = max(5, int(math.floor(0.5 * stoichiometry_r * n_strands)))
    max_p = min(1.0, (2.0 * n_strands) / (4.0 * n_xlinks))
    effective_p = min(conversion, max_p)

    gen = MCUniverseGenerator(box_l, box_l, box_l)
    gen.set_seed(seed)
    gen.config_nr_of_mc_steps(0)
    gen.set_bead_distance(b)
    gen.set_mean_squared_bead_distance(b2)
    gen.add_strands(n_strands, [strand_len] * n_strands)
    gen.add_crosslinkers(n_xlinks, 4)
    gen.link_strands_to_conversion(effective_p)

    es = gen.sample_entanglements(
        nr_of_entanglements_to_sample=int(ent_density * vol),
        upper_cutoff=cutoff,
        seed=f"rve_{seed}",
        filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
    )
    fb_ent = gen.get_mehp_simulator(
        es, entanglement_mode=pt.EntanglementMode.MERGED_LINKS
    )
    fb_ent.run(
        pt.MEHPLinearSolverConfig(
            simplification_mode=pt.StructureSimplificationMode.INACTIVE_THEN_X2F
        )
    )
    g_total = float(gamma_conv * fb_ent.get_gamma_factor(r02, nr_of_chains=1) / vol)
    w_active = float(fb_ent.get_active_weight_fraction())

    fb_ph = gen.get_mehp_simulator()
    fb_ph.run(
        pt.MEHPLinearSolverConfig(
            simplification_mode=pt.StructureSimplificationMode.INACTIVE_THEN_X2F
        )
    )
    g_phantom = float(gamma_conv * fb_ph.get_gamma_factor(r02, nr_of_chains=1) / vol)

    n_atoms = n_strands * strand_len + n_xlinks
    return g_total, g_phantom, w_active, n_strands, n_atoms


def _simulate_rve_sweep_repeat(
    study_kind: _RveStudyKind,
    param_idx: int,
    repeat_idx: int,
    seed: int,
    box_lengths: Sequence[float],
    strand_len: int,
    conversion: float,
    stoichiometry_r: float,
    b: float,
    b2: float,
    r02: float,
    cutoff: float,
    ent_density: float,
    bead_density: float,
    gamma_conv: float,
) -> Tuple[
    _RveStudyKind,
    int,
    int,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Execute all box sizes for a single (study, parameter, repeat) combination."""
    n_boxes = len(box_lengths)
    g_total = np.zeros(n_boxes)
    g_phantom = np.zeros(n_boxes)
    w_active = np.zeros(n_boxes)
    n_strands_arr = np.zeros(n_boxes, dtype=int)
    n_atoms_arr = np.zeros(n_boxes, dtype=int)

    for j, box_l in enumerate(box_lengths):
        g_t, g_p, w_a, n_str, n_at = _run_single_pdms_rve_point(
            box_l=box_l,
            strand_len=strand_len,
            conversion=conversion,
            stoichiometry_r=stoichiometry_r,
            seed=seed + 37 * j,
            b=b,
            b2=b2,
            r02=r02,
            cutoff=cutoff,
            ent_density=ent_density,
            bead_density=bead_density,
            gamma_conv=gamma_conv,
        )
        g_total[j] = g_t
        g_phantom[j] = g_p
        w_active[j] = w_a
        n_strands_arr[j] = n_str
        n_atoms_arr[j] = n_at

    return (
        study_kind,
        param_idx,
        repeat_idx,
        g_total,
        g_phantom,
        w_active,
        n_strands_arr,
        n_atoms_arr,
    )


def run_representative_volume_pdms_study(
    box_lengths: Sequence[float] = (14.0, 20.0, 28.0, 40.0, 56.0, 80.0, 110.0),
    n_values: Sequence[int] = (50, 100, 200, 400),
    p_values: Sequence[float] = (0.75, 0.85, 0.925, 0.98),
    r_values: Sequence[float] = (0.65, 0.80, 1.00, 1.25),
    n_repeats: int = 12,
    n_samples_study: int = 24,
    base_seed: int = 1000,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """Execute the multi-seed Representative Volume Element (RVE) study for PDMS networks.

    Varies box length L across four sub-studies:
      1. Chain length N (at p=0.95, r=1.0)
      2. Cross-link conversion p (at N=100, r=1.0)
      3. Stoichiometric imbalance r (at N=100, p = 0.95 * min(1, 1/r))
      4. Sample count precision vs. box size (at N=100, p=0.925, r=1.0, 24 repeats)
    """
    params = get_parameters_for_polymer("PDMS", parameter_type=ParameterType.GAUSSIAN)
    dist_unit = params.get("distance_units")
    b = params.get("<b>").to(dist_unit).magnitude
    b2 = params.get("<b^2>").to(dist_unit**2).magnitude
    r02 = params.get("R02").to(dist_unit**2).magnitude
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()
    bead_density = params.get_bead_density()
    gamma_conv = params.get_gamma_conversion_factor().to("MPa").magnitude

    n_boxes = len(box_lengths)
    r_conversions = [0.95 * min(1.0, 1.0 / r_val) for r_val in r_values]

    # Preallocate result arrays
    n_g_tot = np.zeros((len(n_values), n_boxes, n_repeats))
    n_g_ph = np.zeros((len(n_values), n_boxes, n_repeats))
    n_w_act = np.zeros((len(n_values), n_boxes, n_repeats))
    n_strands_grid = np.zeros((len(n_values), n_boxes), dtype=int)
    n_atoms_grid = np.zeros((len(n_values), n_boxes), dtype=int)

    p_g_tot = np.zeros((len(p_values), n_boxes, n_repeats))
    p_g_ph = np.zeros((len(p_values), n_boxes, n_repeats))
    p_w_act = np.zeros((len(p_values), n_boxes, n_repeats))

    r_g_tot = np.zeros((len(r_values), n_boxes, n_repeats))
    r_g_ph = np.zeros((len(r_values), n_boxes, n_repeats))
    r_w_act = np.zeros((len(r_values), n_boxes, n_repeats))

    s_g_tot = np.zeros((n_boxes, n_samples_study))
    s_g_ph = np.zeros((n_boxes, n_samples_study))
    s_w_act = np.zeros((n_boxes, n_samples_study))
    s_atoms_grid = np.zeros(n_boxes, dtype=int)

    tasks = []
    # Study 1: Chain length N
    for idx, n_len in enumerate(n_values):
        for rep in range(n_repeats):
            tasks.append(
                (
                    _RveStudyKind.CHAIN_LENGTH,
                    idx,
                    rep,
                    base_seed + 10000 * idx + 503 * rep,
                    box_lengths,
                    int(n_len),
                    0.95,
                    1.0,
                    b,
                    b2,
                    r02,
                    cutoff,
                    ent_density,
                    bead_density,
                    gamma_conv,
                )
            )

    # Study 2: Conversion p
    for idx, p_val in enumerate(p_values):
        for rep in range(n_repeats):
            tasks.append(
                (
                    _RveStudyKind.CONVERSION,
                    idx,
                    rep,
                    base_seed + 100000 + 10000 * idx + 503 * rep,
                    box_lengths,
                    100,
                    float(p_val),
                    1.0,
                    b,
                    b2,
                    r02,
                    cutoff,
                    ent_density,
                    bead_density,
                    gamma_conv,
                )
            )

    # Study 3: Stoichiometric imbalance r
    for idx, (r_val, conv_val) in enumerate(zip(r_values, r_conversions)):
        for rep in range(n_repeats):
            tasks.append(
                (
                    _RveStudyKind.STOICHIOMETRY,
                    idx,
                    rep,
                    base_seed + 200000 + 10000 * idx + 503 * rep,
                    box_lengths,
                    100,
                    float(conv_val),
                    float(r_val),
                    b,
                    b2,
                    r02,
                    cutoff,
                    ent_density,
                    bead_density,
                    gamma_conv,
                )
            )

    # Study 4: Sample count study (N=100, p=0.925, r=1.0)
    for rep in range(n_samples_study):
        tasks.append(
            (
                _RveStudyKind.SAMPLES,
                0,
                rep,
                base_seed + 300000 + 503 * rep,
                box_lengths,
                100,
                0.925,
                1.0,
                b,
                b2,
                r02,
                cutoff,
                ent_density,
                bead_density,
                gamma_conv,
            )
        )

    workers = min(resolve_max_workers(max_workers), len(tasks))
    configure_process_threading(workers)
    total_tasks = len(tasks)

    def _store_result(
        res: Tuple[
            _RveStudyKind,
            int,
            int,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
            np.ndarray,
        ],
    ) -> None:
        kind, p_idx, rep_idx, gt, gp, wa, n_str, n_at = res
        if kind == _RveStudyKind.CHAIN_LENGTH:
            n_g_tot[p_idx, :, rep_idx] = gt
            n_g_ph[p_idx, :, rep_idx] = gp
            n_w_act[p_idx, :, rep_idx] = wa
            n_strands_grid[p_idx, :] = n_str
            n_atoms_grid[p_idx, :] = n_at
        elif kind == _RveStudyKind.CONVERSION:
            p_g_tot[p_idx, :, rep_idx] = gt
            p_g_ph[p_idx, :, rep_idx] = gp
            p_w_act[p_idx, :, rep_idx] = wa
        elif kind == _RveStudyKind.STOICHIOMETRY:
            r_g_tot[p_idx, :, rep_idx] = gt
            r_g_ph[p_idx, :, rep_idx] = gp
            r_w_act[p_idx, :, rep_idx] = wa
        elif kind == _RveStudyKind.SAMPLES:
            s_g_tot[:, rep_idx] = gt
            s_g_ph[:, rep_idx] = gp
            s_w_act[:, rep_idx] = wa
            s_atoms_grid[:] = n_at

    if workers > 1:
        print(
            f"Executing {total_tasks} PDMS RVE sweep tasks in parallel across {workers} workers..."
        )
        completed = 0
        with ProcessPoolExecutor(
            max_workers=workers, initializer=worker_process_initializer
        ) as executor:
            futures = [executor.submit(_simulate_rve_sweep_repeat, *t) for t in tasks]
            for fut in as_completed(futures):
                _store_result(fut.result())
                completed += 1
                if completed % 20 == 0 or completed == total_tasks:
                    print(
                        f"  [PDMS RVE Study] Completed {completed}/{total_tasks} tasks..."
                    )
    else:
        for idx, t in enumerate(tasks, start=1):
            _store_result(_simulate_rve_sweep_repeat(*t))
            if idx % 20 == 0 or idx == total_tasks:
                print(f"  [PDMS RVE Study] Completed {idx}/{total_tasks} tasks...")

    stat_n_gt = compute_ensemble_statistics(n_g_tot, axis=2)
    stat_n_gp = compute_ensemble_statistics(n_g_ph, axis=2)
    stat_n_wa = compute_ensemble_statistics(n_w_act, axis=2)

    stat_p_gt = compute_ensemble_statistics(p_g_tot, axis=2)
    stat_p_gp = compute_ensemble_statistics(p_g_ph, axis=2)
    stat_p_wa = compute_ensemble_statistics(p_w_act, axis=2)

    stat_r_gt = compute_ensemble_statistics(r_g_tot, axis=2)
    stat_r_gp = compute_ensemble_statistics(r_g_ph, axis=2)
    stat_r_wa = compute_ensemble_statistics(r_w_act, axis=2)

    stat_s_gt = compute_ensemble_statistics(s_g_tot, axis=1)
    stat_s_gp = compute_ensemble_statistics(s_g_ph, axis=1)
    stat_s_wa = compute_ensemble_statistics(s_w_act, axis=1)

    n_r0_nm = np.array([math.sqrt(n_len * b2) for n_len in n_values])

    return {
        "box_lengths": np.array(box_lengths, dtype=float),
        "n_repeats": np.array(n_repeats),
        "n_samples_study": np.array(n_samples_study),
        "bead_density": np.array(bead_density),
        "ent_density": np.array(ent_density),
        "r02_nm2": np.array(r02),
        # Study 1: Chain length N
        "n_values": np.array(n_values, dtype=int),
        "n_r0_nm": n_r0_nm,
        "n_strands_grid": n_strands_grid,
        "n_atoms_grid": n_atoms_grid,
        "n_g_total_raw": n_g_tot,
        "n_g_total_mean": stat_n_gt["mean"],
        "n_g_total_std": stat_n_gt["std"],
        "n_g_total_ci": stat_n_gt["ci"],
        "n_g_phantom_mean": stat_n_gp["mean"],
        "n_g_phantom_ci": stat_n_gp["ci"],
        "n_w_active_mean": stat_n_wa["mean"],
        # Study 2: Conversion p
        "p_values": np.array(p_values, dtype=float),
        "p_g_total_raw": p_g_tot,
        "p_g_total_mean": stat_p_gt["mean"],
        "p_g_total_std": stat_p_gt["std"],
        "p_g_total_ci": stat_p_gt["ci"],
        "p_g_phantom_mean": stat_p_gp["mean"],
        "p_g_phantom_ci": stat_p_gp["ci"],
        "p_w_active_mean": stat_p_wa["mean"],
        # Study 3: Stoichiometric imbalance r
        "r_values": np.array(r_values, dtype=float),
        "r_conversions": np.array(r_conversions, dtype=float),
        "r_g_total_raw": r_g_tot,
        "r_g_total_mean": stat_r_gt["mean"],
        "r_g_total_std": stat_r_gt["std"],
        "r_g_total_ci": stat_r_gt["ci"],
        "r_g_phantom_mean": stat_r_gp["mean"],
        "r_g_phantom_ci": stat_r_gp["ci"],
        "r_w_active_mean": stat_r_wa["mean"],
        # Study 4: Sample count study
        "s_atoms_grid": s_atoms_grid,
        "s_g_total_raw": s_g_tot,
        "s_g_total_mean": stat_s_gt["mean"],
        "s_g_total_std": stat_s_gt["std"],
        "s_g_total_ci": stat_s_gt["ci"],
        "s_g_phantom_mean": stat_s_gp["mean"],
        "s_w_active_mean": stat_s_wa["mean"],
    }
