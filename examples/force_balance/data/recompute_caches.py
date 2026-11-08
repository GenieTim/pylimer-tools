#!/usr/bin/env python
"""
Utility script to recompute and update all precomputed cached datasets
for the force_balance examples.

Usage:
    python recompute_caches.py [-j WORKERS] [--parallel-caches]
                               [--link-vs-spring-only] [--slipping-vs-fixed-only] [--sampling-only]
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from pathlib import Path
import sys
from typing import Optional

# Ensure repo root, src, and examples directory are in sys.path
data_dir = Path(__file__).resolve().parent
force_balance_dir = data_dir.parent
examples_dir = force_balance_dir.parent
repo_root = examples_dir.parent
src_dir = repo_root / "src"

for p in [str(src_dir), str(examples_dir), str(
        force_balance_dir), str(repo_root)]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _example_cache import (  # noqa: E402
    get_example_cache_path,
    load_or_compute_cache,
    resolve_max_workers,
)


def recompute_entanglement_link_vs_spring(
        max_workers: Optional[int] = None) -> None:
    """Recompute the multi-seed ensemble dataset for plot_entanglement_link_vs_spring."""
    from plot_entanglement_link_vs_spring import (
        run_entanglement_link_vs_spring_simulation,
    )

    cache_path = get_example_cache_path(
        "entanglement_link_vs_spring_cache.npz", "force_balance"
    )
    print(
        f"\n--- Recomputing Entanglement Link vs. Spring cache ({
            cache_path.name}) ---"
    )
    load_or_compute_cache(
        cache_path,
        run_entanglement_link_vs_spring_simulation,
        recompute=True,
        max_workers=max_workers,
    )


def recompute_slipping_vs_fixed_links(
        max_workers: Optional[int] = None) -> None:
    """Recompute the multi-seed ensemble dataset for plot_slipping_vs_fixed_links."""
    from plot_slipping_vs_fixed_links import run_slipping_vs_fixed_simulation

    cache_path = get_example_cache_path(
        "slipping_vs_fixed_cache.npz", "force_balance")
    print(
        f"\n--- Recomputing Slipping vs. Fixed Links cache ({cache_path.name}) ---")
    load_or_compute_cache(
        cache_path,
        run_slipping_vs_fixed_simulation,
        recompute=True,
        max_workers=max_workers,
    )


def recompute_entanglement_sampling_comparison(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed benchmark dataset for plot_entanglement_sampling_comparison."""
    from plot_entanglement_sampling_comparison import run_all_benchmarks

    cache_path = get_example_cache_path(
        "entanglement_sampling_comparison_cache.npz", "force_balance"
    )
    print(
        f"\n--- Recomputing Entanglement Sampling Comparison cache ({
            cache_path.name
        }) ---"
    )
    load_or_compute_cache(
        cache_path,
        run_all_benchmarks,
        recompute=True,
        max_workers=max_workers,
    )


def recompute_inactive_removal_and_sampling(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed ensemble dataset for plot_inactive_removal_and_sampling."""
    from plot_inactive_removal_and_sampling import (
        run_inactive_removal_and_sampling_simulation,
    )

    cache_path = get_example_cache_path(
        "inactive_removal_and_sampling_cache.npz", "force_balance"
    )
    print(
        f"\n--- Recomputing Inactive Removal and Sampling cache ({
            cache_path.name}) ---"
    )
    load_or_compute_cache(
        cache_path,
        run_inactive_removal_and_sampling_simulation,
        recompute=True,
        max_workers=max_workers,
    )


def _run_cache_task(task_name: str, workers: Optional[int]) -> str:
    """Worker wrapper to execute an individual cache recomputation task."""
    if task_name == "link_vs_spring":
        recompute_entanglement_link_vs_spring(max_workers=workers)
    elif task_name == "slipping_vs_fixed":
        recompute_slipping_vs_fixed_links(max_workers=workers)
    elif task_name == "sampling":
        recompute_entanglement_sampling_comparison(max_workers=workers)
    elif task_name == "inactive_removal":
        recompute_inactive_removal_and_sampling(max_workers=workers)
    else:
        raise ValueError(f"Unknown task name: {task_name}")
    return task_name


def main():
    parser = argparse.ArgumentParser(
        description="Recompute force_balance example caches in parallel."
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=None,
        help="Number of worker processes for parallel simulation tasks (default: all available CPU cores).",
    )
    parser.add_argument(
        "--parallel-caches",
        action="store_true",
        help="Recompute independent cache files concurrently across separate processes.",
    )
    parser.add_argument(
        "--link-vs-spring-only",
        action="store_true",
        help="Recompute only entanglement link vs spring cache",
    )
    parser.add_argument(
        "--slipping-vs-fixed-only",
        action="store_true",
        help="Recompute only slipping vs fixed links cache",
    )
    parser.add_argument(
        "--sampling-only",
        action="store_true",
        help="Recompute only entanglement sampling comparison cache",
    )
    parser.add_argument(
        "--inactive-removal-only",
        action="store_true",
        help="Recompute only inactive removal and sampling cache",
    )
    args = parser.parse_args()

    total_workers = resolve_max_workers(args.workers)

    tasks_to_run = []
    run_all = not (
        args.link_vs_spring_only
        or args.slipping_vs_fixed_only
        or args.sampling_only
        or args.inactive_removal_only
    )
    if run_all or args.link_vs_spring_only:
        tasks_to_run.append("link_vs_spring")
    if run_all or args.slipping_vs_fixed_only:
        tasks_to_run.append("slipping_vs_fixed")
    if run_all or args.sampling_only:
        tasks_to_run.append("sampling")
    if run_all or args.inactive_removal_only:
        tasks_to_run.append("inactive_removal")

    print(
        f"Recomputing {len(tasks_to_run)} cache task(s) with total available workers = {
            total_workers
        }..."
    )

    if args.parallel_caches and len(tasks_to_run) > 1:
        workers_per_cache = max(1, total_workers // len(tasks_to_run))
        print(
            f"Executing {len(tasks_to_run)} caches in parallel ({
                workers_per_cache
            } worker(s) per cache)..."
        )
        with ProcessPoolExecutor(max_workers=len(tasks_to_run)) as pool:
            futures = [
                pool.submit(_run_cache_task, t, workers_per_cache) for t in tasks_to_run
            ]
            for fut in as_completed(futures):
                task_name = fut.result()
                print(f"Task '{task_name}' completed successfully.")
    else:
        for t in tasks_to_run:
            _run_cache_task(t, total_workers)

    print("\nAll requested caches have been recomputed successfully!")


if __name__ == "__main__":
    main()
