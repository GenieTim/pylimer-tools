#!/usr/bin/env python
"""
Utility script to recompute and update all precomputed cached datasets
for the force_balance examples.

Usage:
    python recompute_caches.py [--sliplink-only] [--sampling-only]
"""

import argparse
from pathlib import Path
import sys

# Ensure repo root, src, and examples directory are in sys.path
data_dir = Path(__file__).resolve().parent
force_balance_dir = data_dir.parent
repo_root = force_balance_dir.parents[1]
src_dir = repo_root / "src"

for p in [str(src_dir), str(force_balance_dir), str(repo_root)]:
    if p not in sys.path:
        sys.path.insert(0, p)

from pylimer_tools.utils.example_cache import get_example_cache_path, load_or_compute_cache


def recompute_sliplink_vs_slipspring():
    """Recompute the multi-seed ensemble dataset for plot_sliplink_vs_slipspring."""
    from plot_sliplink_vs_slipspring import run_sliplink_vs_slipspring_simulation

    cache_path = get_example_cache_path("sliplink_vs_slipspring_cache.npz", "force_balance")
    print(f"\n--- Recomputing Slip-Link vs. Slip-Spring cache ({cache_path.name}) ---")
    load_or_compute_cache(cache_path, run_sliplink_vs_slipspring_simulation, recompute=True)


def recompute_entanglement_sampling_comparison():
    """Recompute the multi-seed benchmark dataset for plot_entanglement_sampling_comparison."""
    from plot_entanglement_sampling_comparison import run_all_benchmarks

    cache_path = get_example_cache_path("entanglement_sampling_comparison_cache.npz", "force_balance")
    print(f"\n--- Recomputing Entanglement Sampling Comparison cache ({cache_path.name}) ---")
    load_or_compute_cache(cache_path, run_all_benchmarks, recompute=True)


def main():
    parser = argparse.ArgumentParser(description="Recompute force_balance example caches.")
    parser.add_argument("--sliplink-only", action="store_true", help="Recompute only slip-link vs slip-spring cache")
    parser.add_argument("--sampling-only", action="store_true", help="Recompute only entanglement sampling comparison cache")
    args = parser.parse_args()

    run_all = not (args.sliplink_only or args.sampling_only)
    if run_all or args.sliplink_only:
        recompute_sliplink_vs_slipspring()
    if run_all or args.sampling_only:
        recompute_entanglement_sampling_comparison()

    print("\nAll requested caches have been recomputed successfully!")


if __name__ == "__main__":
    main()
