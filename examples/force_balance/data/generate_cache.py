#!/usr/bin/env python
"""
Generate cached data for plot_entanglement_sampling_comparison.py
"""

import os
import sys
import time
import subprocess
import json
from pathlib import Path
import numpy as np

import pylimer_tools_cpp as pt
from pylimer_tools.io.bead_spring_parameter_provider import (
    ParameterType,
    get_parameters_for_polymer,
)


def run_scaling_benchmark(sizes, strand_len=50, seed=42):
    """Measure peak memory and runtime for discrete and analytical methods across system sizes."""
    print("Running Scaling Benchmark...")
    params = get_parameters_for_polymer(
        "PDMS", parameter_type=ParameterType.GAUSSIAN)
    bead_density = params.get_bead_density()
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()

    worker_code = """
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

t0 = time.perf_counter()
gen = pt.MCUniverseGenerator(L, L, L)
gen.set_seed(seed)
gen.set_mean_squared_bead_distance(1.0)
gen.add_strands(n_strands, [strand_len] * n_strands)
gen.add_crosslinkers(int(0.5 * n_strands), 4)
gen.link_strands_to_conversion(0.92)

if mode == "discrete":
    u = gen.get_universe()
    fb = pt.MEHPForceBalance2(
        u,
        nr_of_entanglements_to_sample=n_ent,
        upper_sampling_cutoff=cutoff,
        entanglements_as_springs=False,
    )
    fb.validate_network()
elif mode == "analytical":
    es = gen.sample_entanglements(n_ent, cutoff, 0.0, 0, 3.0, "seed", True)
    fb = gen.get_force_balance2(es, 2, False)
    fb.validate_network()

t1 = time.perf_counter()
maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
if sys.platform == "darwin":
    rss_mb = maxrss / (1024.0 * 1024.0)
else:
    rss_mb = maxrss / 1024.0

print(json.dumps({"mode": mode, "n_strands": n_strands, "n_atoms": n_atoms, "time": t1 - t0, "rss_mb": rss_mb}))
"""

    disc_atoms = []
    disc_time = []
    disc_rss = []

    ana_atoms = []
    ana_time = []
    ana_rss = []

    for n in sizes:
        for mode in ["analytical", "discrete"]:
            cmd = [
                sys.executable,
                "-c",
                worker_code,
                mode,
                str(n),
                str(strand_len),
                str(seed),
                str(bead_density),
                str(cutoff),
                str(ent_density),
            ]
            res = subprocess.run(
                cmd, capture_output=True, text=True, check=True)
            lines = [l for l in res.stdout.strip().split("\n")
                     if l.startswith("{")]
            data = json.loads(lines[-1])
            if mode == "discrete":
                disc_atoms.append(data["n_atoms"])
                disc_time.append(data["time"])
                disc_rss.append(data["rss_mb"])
                print(
                    f"Discrete:   N_atoms={data['n_atoms']:6d} | Time={
                        data['time']:6.2f}s | RSS={data['rss_mb']:6.1f} MB"
                )
            else:
                ana_atoms.append(data["n_atoms"])
                ana_time.append(data["time"])
                ana_rss.append(data["rss_mb"])
                print(
                    f"Analytical: N_atoms={data['n_atoms']:6d} | Time={
                        data['time']:6.2f}s | RSS={data['rss_mb']:6.1f} MB"
                )

    return {
        "scaling_atoms": np.array(ana_atoms),
        "scaling_time_discrete": np.array(disc_time),
        "scaling_rss_discrete": np.array(disc_rss),
        "scaling_time_analytical": np.array(ana_time),
        "scaling_rss_analytical": np.array(ana_rss),
    }


def run_distribution_comparison(
    n_strands=2000, strand_len=50, conversion=0.95, seed=123
):
    """Sample entanglements and compare contour position alpha and spatial distance distributions."""
    print("Running Distribution Comparison...")
    params = get_parameters_for_polymer(
        "PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = params.get("<b>").to(params.get("distance_units")).magnitude
    bead_density = params.get_bead_density()
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()

    n_atoms = n_strands * strand_len + int(0.5 * n_strands)
    vol = n_atoms / bead_density
    L = vol ** (1.0 / 3.0)
    n_ent = int(ent_density * vol)

    gen = pt.MCUniverseGenerator(L, L, L)
    gen.set_seed(seed)
    gen.set_bead_distance(b)
    gen.set_mean_squared_bead_distance(b**2)
    gen.add_strands(n_strands, [strand_len] * n_strands)
    gen.add_crosslinkers(int(0.5 * n_strands), 4)
    gen.link_strands_to_conversion(conversion)

    # 1. Analytical sampling
    es = gen.sample_entanglements(
        n_ent, cutoff, 0.0, 0, 3.0, "ana_seed", False)
    alpha_ana = np.concatenate([np.array(es.alpha_1), np.array(es.alpha_2)])

    pos1 = np.array(es.position_1).reshape(-1, 3)
    pos2 = np.array(es.position_2).reshape(-1, 3)
    diff = pos2 - pos1
    # minimum image convention
    diff = diff - np.round(diff / L) * L
    dist_ana = np.linalg.norm(diff, axis=1)

    # 2. Discrete sampling on Universe
    u = gen.get_universe()
    cu = gen.get_crosslinker_universe()
    atom_pairs = pt.randomly_sample_entanglements(
        u, n_ent, cutoff, 0.0, 0, 3.0, "disc_seed"
    )

    n_xlinks = len(cu.xlink_x)
    offsets = [n_xlinks]
    for n_beads in cu.beads_in_strand:
        offsets.append(offsets[-1] + n_beads)
    offsets = np.array(offsets)

    alpha_disc = []
    dist_disc = []
    box = u.get_box()

    for a1, a2 in atom_pairs.pairs_of_atoms:
        idx1 = a1 - 1
        idx2 = a2 - 1
        # only calculate contour alpha if both are strand beads (idx >=
        # n_xlinks)
        if idx1 >= n_xlinks:
            s1 = np.searchsorted(offsets, idx1, side="right") - 1
            b1 = idx1 - offsets[s1]
            alpha_disc.append((b1 + 1.0) / (strand_len + 1.0))
        if idx2 >= n_xlinks:
            s2 = np.searchsorted(offsets, idx2, side="right") - 1
            b2 = idx2 - offsets[s2]
            alpha_disc.append((b2 + 1.0) / (strand_len + 1.0))

        atom1 = u.get_atom(a1)
        atom2 = u.get_atom(a2)
        dist_disc.append(atom1.distance_to(atom2, box))

    print(
        f"Analytical samples: {
            len(dist_ana)}, Discrete samples: {
            len(dist_disc)}")
    return {
        "alpha_analytical": np.array(alpha_ana),
        "dist_analytical": np.array(dist_ana),
        "alpha_discrete": np.array(alpha_disc),
        "dist_discrete": np.array(dist_disc),
    }


def run_conversion_sweep(
    conversions=(0.75, 0.80, 0.85, 0.90, 0.95, 1.00),
    n_strands=1000,
    strand_len=50,
    seed=777,
):
    """Compare modulus, active fractions, and active entanglement count across conversions."""
    print("Running Crosslinker Conversion Sweep...")
    params = get_parameters_for_polymer(
        "PDMS", parameter_type=ParameterType.GAUSSIAN)
    b = params.get("<b>").to(params.get("distance_units")).magnitude
    bead_density = params.get_bead_density()
    cutoff = params.get_sampling_cutoff()
    ent_density = params.get_entanglement_density()
    r02 = params.get("R02").to(params.get("distance_units") ** 2).magnitude
    kbt = params.get("T") * params.get("kb")
    gamma_conv = (
        (kbt / ((params.get("distance_units")) ** 3)).to("MPa")).magnitude

    n_atoms = n_strands * strand_len + int(0.5 * n_strands)
    vol = n_atoms / bead_density
    L = vol ** (1.0 / 3.0)
    n_ent = int(ent_density * vol)

    res = {
        "conversions": np.array(conversions),
        "g_discrete": [],
        "w_active_discrete": [],
        "g_ana_filtered": [],
        "w_active_ana_filtered": [],
        "n_ent_ana_filtered": [],
        "g_ana_unfiltered": [],
        "w_active_ana_unfiltered": [],
        "n_ent_ana_unfiltered": [],
        "g_ana_springs": [],
    }

    for p in conversions:
        gen = pt.MCUniverseGenerator(L, L, L)
        gen.set_seed(seed)
        gen.set_bead_distance(b)
        gen.set_mean_squared_bead_distance(b**2)
        gen.add_strands(n_strands, [strand_len] * n_strands)
        gen.add_crosslinkers(int(0.5 * n_strands), 4)
        gen.link_strands_to_conversion(p)

        # 1. Discrete Full Universe
        u = gen.get_universe()
        fb_disc = pt.MEHPForceBalance2(
            u,
            nr_of_entanglements_to_sample=n_ent,
            upper_sampling_cutoff=cutoff,
            entanglements_as_springs=False,
        )
        fb_disc.run_force_relaxation()
        g_d = gamma_conv * fb_disc.get_gamma_factor(r02) / vol
        w_d = fb_disc.get_active_weight_fraction()
        res["g_discrete"].append(g_d)
        res["w_active_discrete"].append(w_d)

        # 2. Analytical Filtered (Default)
        es_filt = gen.sample_entanglements(
            n_ent, cutoff, 0.0, 0, 3.0, "seed_f", True)
        fb_filt = gen.get_force_balance2(es_filt, 2, False)
        fb_filt.run_force_relaxation()
        g_af = gamma_conv * fb_filt.get_gamma_factor(r02) / vol
        w_af = fb_filt.get_active_weight_fraction()
        res["g_ana_filtered"].append(g_af)
        res["w_active_ana_filtered"].append(w_af)
        res["n_ent_ana_filtered"].append(len(es_filt))

        # 3. Analytical Unfiltered
        es_unfilt = gen.sample_entanglements(
            n_ent, cutoff, 0.0, 0, 3.0, "seed_u", False
        )
        fb_unfilt = gen.get_force_balance2(es_unfilt, 2, False)
        fb_unfilt.run_force_relaxation()
        g_au = gamma_conv * fb_unfilt.get_gamma_factor(r02) / vol
        w_au = fb_unfilt.get_active_weight_fraction()
        res["g_ana_unfiltered"].append(g_au)
        res["w_active_ana_unfiltered"].append(w_au)
        res["n_ent_ana_unfiltered"].append(len(es_unfilt))

        # 4. Analytical with Hookean Springs
        fb_springs = gen.get_force_balance2(es_filt, 2, True)
        fb_springs.run_force_relaxation()
        g_as = gamma_conv * fb_springs.get_gamma_factor(r02) / vol
        res["g_ana_springs"].append(g_as)

        print(
            f"p={p:.2f} | Discrete G={g_d:.4f} MPa (w={w_d:.3f}) | "
            f"Ana-Filt G={g_af:.4f} MPa (w={w_af:.3f}, n={len(es_filt)}) | "
            f"Ana-Unfilt G={
                g_au:.4f} MPa (w={
                w_au:.3f}, n={
                len(es_unfilt)}) | "
            f"Ana-Springs G={g_as:.4f} MPa"
        )

    for k in [
        "g_discrete",
        "w_active_discrete",
        "g_ana_filtered",
        "w_active_ana_filtered",
        "n_ent_ana_filtered",
        "g_ana_unfiltered",
        "w_active_ana_unfiltered",
        "n_ent_ana_unfiltered",
        "g_ana_springs",
    ]:
        res[k] = np.array(res[k])
    res["target_n_ent"] = n_ent

    return res


def main():
    data_dir = Path(__file__).resolve().parent
    cache_file = data_dir / "entanglement_sampling_comparison_cache.npz"

    scaling_sizes = [200, 500, 1000, 2000, 4000]
    scaling_data = run_scaling_benchmark(scaling_sizes, strand_len=50)
    dist_data = run_distribution_comparison(
        n_strands=2000, strand_len=50, conversion=0.95
    )
    conv_data = run_conversion_sweep(
        conversions=[0.75, 0.80, 0.85, 0.90, 0.95, 1.00], n_strands=1000, strand_len=50
    )

    combined = {}
    combined.update(scaling_data)
    combined.update(dist_data)
    combined.update(conv_data)

    np.savez_compressed(cache_file, **combined)
    print(
        f"Successfully saved cache to {cache_file} ({
            os.path.getsize(cache_file) / 1024:.1f} KB)"
    )


if __name__ == "__main__":
    main()
