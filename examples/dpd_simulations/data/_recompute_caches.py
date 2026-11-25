#!/usr/bin/env python
"""
Utility script to recompute and update precomputed cached datasets
for the dpd_simulations examples.

Usage:
    python _recompute_caches.py [-j WORKERS] [--rheology-only] [--curing-only]
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Optional

data_dir = Path(__file__).resolve().parent
dpd_dir = data_dir.parent
examples_dir = dpd_dir.parent
repo_root = examples_dir.parent
src_dir = repo_root / "src"

for p in (str(src_dir), str(examples_dir), str(dpd_dir), str(repo_root)):
    if p not in sys.path:
        sys.path.insert(0, p)

from _dpd_curing_runner import run_dpd_curing_study  # noqa: E402
from _dpd_rheology_runner import run_dpd_rheology_study  # noqa: E402
from _example_cache import get_example_cache_path, load_or_compute_cache  # noqa: E402


def recompute_dpd_rheology(max_workers: Optional[int] = None) -> None:
    """Recompute the multi-seed dataset for plot_dpd_slipspring_rheology."""
    cache_path = get_example_cache_path(
        "dpd_slipspring_rheology_cache.npz", "dpd_simulations"
    )
    print(f"\n--- Recomputing DPD Slip-Spring Rheology cache ({cache_path.name}) ---")
    load_or_compute_cache(
        cache_path,
        run_dpd_rheology_study,
        recompute=True,
        max_workers=max_workers,
    )


def recompute_dpd_curing(max_workers: Optional[int] = None) -> None:
    """Recompute the multi-seed dataset for plot_dpd_dynamic_curing_vs_mc."""
    cache_path = get_example_cache_path(
        "dpd_dynamic_curing_vs_mc_cache.npz", "dpd_simulations"
    )
    print(f"\n--- Recomputing DPD Dynamic Curing cache ({cache_path.name}) ---")
    load_or_compute_cache(
        cache_path,
        run_dpd_curing_study,
        recompute=True,
        max_workers=max_workers,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recompute cached datasets for dpd_simulations examples."
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=None,
        help="Number of parallel worker processes (default: all CPU cores).",
    )
    parser.add_argument(
        "--rheology-only",
        action="store_true",
        help="Recompute only dpd_slipspring_rheology_cache.npz.",
    )
    parser.add_argument(
        "--curing-only",
        action="store_true",
        help="Recompute only dpd_dynamic_curing_vs_mc_cache.npz.",
    )
    args = parser.parse_args()
    run_all = not (args.rheology_only or args.curing_only)
    if run_all or args.rheology_only:
        recompute_dpd_rheology(max_workers=args.workers)
    if run_all or args.curing_only:
        recompute_dpd_curing(max_workers=args.workers)


if __name__ == "__main__":
    main()
