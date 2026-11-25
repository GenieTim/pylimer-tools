#!/usr/bin/env python
"""
Utility script to recompute and update precomputed cached datasets
for the normal_mode_analysis examples.

Usage:
    python _recompute_caches.py [-j WORKERS]
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Optional

data_dir = Path(__file__).resolve().parent
nma_dir = data_dir.parent
examples_dir = nma_dir.parent
repo_root = examples_dir.parent
src_dir = repo_root / "src"

for p in (str(src_dir), str(examples_dir), str(nma_dir), str(repo_root)):
    if p not in sys.path:
        sys.path.insert(0, p)

from _example_cache import get_example_cache_path, load_or_compute_cache  # noqa: E402
from _kpm_limitations_runner import run_kpm_limitations_study  # noqa: E402
from _nma_study_runners import run_kpm_vs_exact_pdms_study  # noqa: E402


def recompute_kpm_vs_exact_pdms(max_workers: Optional[int] = None) -> None:
    """Recompute the multi-seed dataset for plot_kpm_vs_exact_pdms."""
    cache_path = get_example_cache_path(
        "kpm_vs_exact_pdms_cache.npz", "normal_mode_analysis"
    )
    print(f"\n--- Recomputing KPM vs. Exact PDMS cache ({cache_path.name}) ---")
    load_or_compute_cache(
        cache_path,
        run_kpm_vs_exact_pdms_study,
        recompute=True,
        max_workers=max_workers,
    )


def recompute_kpm_limitations_and_hybrid(max_workers: Optional[int] = None) -> None:
    """Recompute the multi-seed dataset for plot_kpm_limitations_and_hybrid."""
    cache_path = get_example_cache_path(
        "kpm_limitations_and_hybrid_cache.npz", "normal_mode_analysis"
    )
    print(f"\n--- Recomputing KPM Limitations & Hybrid cache ({cache_path.name}) ---")
    load_or_compute_cache(
        cache_path,
        run_kpm_limitations_study,
        recompute=True,
        max_workers=max_workers,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recompute cached datasets for normal_mode_analysis examples."
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=None,
        help="Number of parallel worker processes (default: all CPU cores).",
    )
    parser.add_argument(
        "--kpm-vs-exact-only",
        action="store_true",
        help="Recompute only kpm_vs_exact_pdms_cache.npz.",
    )
    parser.add_argument(
        "--kpm-limitations-only",
        action="store_true",
        help="Recompute only kpm_limitations_and_hybrid_cache.npz.",
    )
    args = parser.parse_args()
    run_all = not (args.kpm_vs_exact_only or args.kpm_limitations_only)
    if run_all or args.kpm_vs_exact_only:
        recompute_kpm_vs_exact_pdms(max_workers=args.workers)
    if run_all or args.kpm_limitations_only:
        recompute_kpm_limitations_and_hybrid(max_workers=args.workers)


if __name__ == "__main__":
    main()
