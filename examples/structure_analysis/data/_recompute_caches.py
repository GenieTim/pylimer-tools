#!/usr/bin/env python
"""
Utility script to recompute and update precomputed cached datasets
for the structure_analysis examples.

Usage:
    python _recompute_caches.py [-j WORKERS]
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Optional

data_dir = Path(__file__).resolve().parent
struct_dir = data_dir.parent
examples_dir = struct_dir.parent
repo_root = examples_dir.parent
src_dir = repo_root / "src"

for p in (str(src_dir), str(examples_dir), str(struct_dir), str(repo_root)):
    if p not in sys.path:
        sys.path.insert(0, p)

from _example_cache import get_example_cache_path, load_or_compute_cache  # noqa: E402
from _resistance_mehp_runner import run_resistance_mehp_study  # noqa: E402


def recompute_resistance_mehp(max_workers: Optional[int] = None) -> None:
    """Recompute the multi-seed dataset for plot_resistance_vs_mehp_stress."""
    cache_path = get_example_cache_path(
        "resistance_vs_mehp_stress_cache.npz", "structure_analysis"
    )
    print(f"\n--- Recomputing Resistance vs. MEHP Stress cache ({cache_path.name}) ---")
    load_or_compute_cache(
        cache_path,
        run_resistance_mehp_study,
        recompute=True,
        max_workers=max_workers,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recompute cached datasets for structure_analysis examples."
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=None,
        help="Number of parallel worker processes (default: all CPU cores).",
    )
    args = parser.parse_args()
    recompute_resistance_mehp(max_workers=args.workers)


if __name__ == "__main__":
    main()
