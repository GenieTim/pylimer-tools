#!/usr/bin/env python
"""
Utility script to recompute and update all precomputed cached datasets
for the force_balance examples.

Usage:
    python _recompute_caches.py [-j WORKERS] [--parallel-caches]
                               [--link-vs-spring-only] [--slipping-vs-fixed-only] [--sampling-only]
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from enum import Enum
from pathlib import Path
import sys
from typing import Callable, Dict, Optional, Tuple

# Ensure repo root, src, and examples directory are in sys.path
data_dir = Path(__file__).resolve().parent
force_balance_dir = data_dir.parent
examples_dir = force_balance_dir.parent
repo_root = examples_dir.parent
src_dir = repo_root / "src"

for p in [str(src_dir), str(examples_dir), str(force_balance_dir), str(repo_root)]:
    if p not in sys.path:
        sys.path.insert(0, p)

from _example_cache import (  # noqa: E402
    configure_process_threading,
    get_example_cache_path,
    load_or_compute_cache,
    resolve_max_workers,
    worker_process_initializer,
)


class CacheTask(Enum):
    """Enumeration of precomputed dataset caches for the force_balance examples."""

    LINK_VS_SPRING = "link_vs_spring"
    SLIPPING_VS_FIXED = "slipping_vs_fixed"
    SAMPLING = "sampling"
    INACTIVE_REMOVAL = "inactive_removal"
    SAMPLING_CUTOFF = "sampling_cutoff"
    POLYDISPERSITY = "polydispersity"
    RVE_PDMS = "rve_pdms"


def _recompute_cache(
    cache_filename: str,
    title: str,
    runner: Callable[..., Dict[str, object]],
    max_workers: Optional[int] = None,
) -> None:
    """Recompute and save a single example cache dataset."""
    cache_path = get_example_cache_path(cache_filename, "force_balance")
    print(f"\n--- Recomputing {title} cache ({cache_path.name}) ---")
    load_or_compute_cache(
        cache_path,
        runner,
        recompute=True,
        max_workers=max_workers,
    )


def recompute_entanglement_link_vs_spring(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed ensemble dataset for plot_entanglement_link_vs_spring."""
    from _study_runners import run_entanglement_link_vs_spring_simulation

    _recompute_cache(
        "entanglement_link_vs_spring_cache.npz",
        "Entanglement Link vs. Spring",
        run_entanglement_link_vs_spring_simulation,
        max_workers=max_workers,
    )


def recompute_slipping_vs_fixed_links(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed ensemble dataset for plot_slipping_vs_fixed_links."""
    from _study_runners import run_slipping_vs_fixed_simulation

    _recompute_cache(
        "slipping_vs_fixed_cache.npz",
        "Slipping vs. Fixed Links",
        run_slipping_vs_fixed_simulation,
        max_workers=max_workers,
    )


def recompute_entanglement_sampling_comparison(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed benchmark dataset for plot_entanglement_sampling_comparison."""
    from _study_runners import run_all_benchmarks

    _recompute_cache(
        "entanglement_sampling_comparison_cache.npz",
        "Entanglement Sampling Comparison",
        run_all_benchmarks,
        max_workers=max_workers,
    )


def recompute_inactive_removal_and_sampling(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed ensemble dataset for plot_inactive_removal_and_sampling."""
    from _study_runners import run_inactive_removal_and_sampling_simulation

    _recompute_cache(
        "inactive_removal_and_sampling_cache.npz",
        "Inactive Removal and Sampling",
        run_inactive_removal_and_sampling_simulation,
        max_workers=max_workers,
    )


def recompute_sampling_cutoff_effect(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed ensemble dataset for plot_sampling_cutoff_effect."""
    from _study_runners import run_sampling_cutoff_simulation

    _recompute_cache(
        "sampling_cutoff_effect_cache.npz",
        "Sampling Cutoff Effect",
        run_sampling_cutoff_simulation,
        max_workers=max_workers,
    )


def recompute_polydispersity_study(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed ensemble dataset for plot_polydispersity_study."""
    from _study_runners import run_polydispersity_simulation

    _recompute_cache(
        "polydispersity_study_cache.npz",
        "Polydispersity Study",
        run_polydispersity_simulation,
        max_workers=max_workers,
    )


def recompute_representative_volume_pdms(
    max_workers: Optional[int] = None,
) -> None:
    """Recompute the multi-seed ensemble dataset for plot_representative_volume_pdms."""
    from _study_runners import run_representative_volume_pdms_study

    _recompute_cache(
        "representative_volume_pdms_cache.npz",
        "Representative Volume PDMS",
        run_representative_volume_pdms_study,
        max_workers=max_workers,
    )


_TASK_HANDLERS: Dict[CacheTask, Callable[[Optional[int]], None]] = {
    CacheTask.LINK_VS_SPRING: recompute_entanglement_link_vs_spring,
    CacheTask.SLIPPING_VS_FIXED: recompute_slipping_vs_fixed_links,
    CacheTask.SAMPLING: recompute_entanglement_sampling_comparison,
    CacheTask.INACTIVE_REMOVAL: recompute_inactive_removal_and_sampling,
    CacheTask.SAMPLING_CUTOFF: recompute_sampling_cutoff_effect,
    CacheTask.POLYDISPERSITY: recompute_polydispersity_study,
    CacheTask.RVE_PDMS: recompute_representative_volume_pdms,
}


def _run_cache_task(task: CacheTask, workers: Optional[int]) -> CacheTask:
    """Worker wrapper to execute an individual cache recomputation task."""
    handler = _TASK_HANDLERS.get(task)
    if handler is None:
        raise ValueError(f"Unknown cache task: {task}")
    handler(workers)
    return task


def main() -> None:
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
    parser.add_argument(
        "--sampling-cutoff-only",
        action="store_true",
        help="Recompute only sampling cutoff effect cache",
    )
    parser.add_argument(
        "--polydispersity-only",
        action="store_true",
        help="Recompute only polydispersity study cache",
    )
    parser.add_argument(
        "--rve-pdms-only",
        action="store_true",
        help="Recompute only representative volume PDMS study cache",
    )
    args = parser.parse_args()

    total_workers = resolve_max_workers(args.workers)
    configure_process_threading(total_workers)

    flag_to_task: Tuple[Tuple[bool, CacheTask], ...] = (
        (args.link_vs_spring_only, CacheTask.LINK_VS_SPRING),
        (args.slipping_vs_fixed_only, CacheTask.SLIPPING_VS_FIXED),
        (args.sampling_only, CacheTask.SAMPLING),
        (args.inactive_removal_only, CacheTask.INACTIVE_REMOVAL),
        (args.sampling_cutoff_only, CacheTask.SAMPLING_CUTOFF),
        (args.polydispersity_only, CacheTask.POLYDISPERSITY),
        (args.rve_pdms_only, CacheTask.RVE_PDMS),
    )
    run_all = not any(flag for flag, _ in flag_to_task)
    tasks_to_run = [task for flag, task in flag_to_task if run_all or flag]

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
        with ProcessPoolExecutor(
            max_workers=len(tasks_to_run),
            initializer=worker_process_initializer,
        ) as pool:
            futures = [
                pool.submit(_run_cache_task, t, workers_per_cache) for t in tasks_to_run
            ]
            for fut in as_completed(futures):
                completed_task = fut.result()
                print(f"Task '{completed_task.value}' completed successfully.")
    else:
        for t in tasks_to_run:
            _run_cache_task(t, total_workers)

    print("\nAll requested caches have been recomputed successfully!")


if __name__ == "__main__":
    main()
