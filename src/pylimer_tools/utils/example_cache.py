"""
Caching and statistical utilities for Sphinx-Gallery examples.

Provides standardized functions to resolve cache file locations, compute ensemble
statistics (means, standard errors, and Student's t confidence intervals), and
transparently cache/load precomputed simulation data to accelerate documentation builds
and CI while enabling large-scale, statistically robust example studies.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
from typing import Any, Callable, Dict, Optional, Union

import numpy as np


def get_example_cache_path(filename: str, example_folder: str = "force_balance") -> Path:
    """Resolve the filesystem path to an example's cached data file (.npz).

    Searches intelligently from current working directory, caller file location,
    and repository root to find existing data files, or creates the target directory
    structure for new caches.

    Parameters
    ----------
    filename : str
        Name of the cache file (e.g. "sliplink_vs_slipspring_cache.npz").
    example_folder : str, default="force_balance"
        Subfolder within ``examples/`` where the cache resides.

    Returns
    -------
    Path
        Absolute path to the cache file.
    """
    target_path = Path(filename)
    if target_path.is_file():
        return target_path.resolve()

    # Search candidates from caller file, cwd, and parent directories
    candidates = []
    if "__file__" in globals() and globals()["__file__"]:
        candidates.append(Path(__file__).resolve().parent)
    candidates.append(Path.cwd().resolve())

    # Traverse upward from candidate directories
    for start_dir in list(candidates):
        cur = start_dir
        for p in [cur] + list(cur.parents):
            # Check p/data/filename
            if (p / "data" / filename).is_file():
                return (p / "data" / filename).resolve()
            # Check p/examples/<example_folder>/data/filename
            if (p / "examples" / example_folder / "data" / filename).is_file():
                return (p / "examples" / example_folder / "data" / filename).resolve()
            # Check p/<example_folder>/data/filename
            if (p / example_folder / "data" / filename).is_file():
                return (p / example_folder / "data" / filename).resolve()

    # If not found, resolve default location for creation
    for start_dir in list(candidates):
        cur = start_dir
        for p in [cur] + list(cur.parents):
            target_dir = p / "examples" / example_folder / "data"
            if target_dir.is_dir():
                return (target_dir / filename).resolve()
            if p.name == example_folder and (p / "data").is_dir():
                return (p / "data" / filename).resolve()

    # Fallback: create path relative to cwd or repo root
    fallback_dir = Path.cwd() / "examples" / example_folder / "data"
    fallback_dir.mkdir(parents=True, exist_ok=True)
    return (fallback_dir / filename).resolve()


def compute_ensemble_statistics(
    data: np.ndarray | list[Any],
    axis: int = -1,
    confidence: float = 0.95,
) -> Dict[str, np.ndarray]:
    """Compute ensemble mean, standard deviation, standard error, and Student's t confidence interval.

    Parameters
    ----------
    data : np.ndarray or array-like
        Array of observations where one axis corresponds to independent repetitions/seeds.
    axis : int, default=-1
        Axis representing independent repetitions.
    confidence : float, default=0.95
        Confidence level for the two-sided interval (e.g. 0.95 for 95% CI).

    Returns
    -------
    dict
        Dictionary containing:
        - ``"mean"``: Sample mean along ``axis``.
        - ``"std"``: Sample standard deviation (ddof=1) along ``axis``.
        - ``"sem"``: Standard error of the mean along ``axis``.
        - ``"ci"``: Half-width of the confidence interval (:math:`t_{\\text{crit}} \\times \\text{SEM}`).
        - ``"ci_lower"``: Lower bound (:math:`\\text{mean} - \\text{ci}`).
        - ``"ci_upper"``: Upper bound (:math:`\\text{mean} + \\text{ci}`).
        - ``"n"``: Number of independent repetitions.
    """
    arr = np.asarray(data, dtype=float)
    n = arr.shape[axis]

    mean = np.mean(arr, axis=axis)
    if n <= 1:
        zero_arr = np.zeros_like(mean)
        return {
            "mean": mean,
            "std": zero_arr,
            "sem": zero_arr,
            "ci": zero_arr,
            "ci_lower": mean,
            "ci_upper": mean,
            "n": np.array(n),
        }

    std = np.std(arr, axis=axis, ddof=1)
    sem = std / np.sqrt(n)

    # Calculate Student's t critical value
    try:
        from scipy import stats

        t_crit = float(stats.t.ppf((1.0 + confidence) / 2.0, df=n - 1))
    except Exception:
        # Fallback t-table approximations if scipy is unavailable
        t_table = {
            2: 12.706,
            3: 4.303,
            4: 3.182,
            5: 2.776,
            6: 2.571,
            7: 2.447,
            8: 2.365,
            9: 2.306,
            10: 2.262,
        }
        t_crit = t_table.get(n, 1.96 if confidence == 0.95 else 2.58)

    ci = t_crit * sem
    return {
        "mean": mean,
        "std": std,
        "sem": sem,
        "ci": ci,
        "ci_lower": mean - ci,
        "ci_upper": mean + ci,
        "n": np.array(n),
    }


def should_recompute(recompute_flag: bool = False) -> bool:
    """Determine whether cached data should be recomputed.

    Evaluates CLI arguments (``--recompute``), environment variables
    (``PYLIMER_RECOMPUTE_CACHE=1``), and explicit function arguments.
    """
    if recompute_flag:
        return True
    if "--recompute" in sys.argv:
        return True
    if os.environ.get("PYLIMER_RECOMPUTE_CACHE", "").lower() in ("1", "true", "yes"):
        return True
    return False


def load_or_compute_cache(
    cache_path: Union[str, Path],
    compute_fn: Callable[..., Dict[str, Any]],
    *args: Any,
    recompute: bool = False,
    verbose: bool = True,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Load precomputed data from cache or execute computation and save results.

    Parameters
    ----------
    cache_path : str or Path
        Target file path (.npz).
    compute_fn : callable
        Function returning a dictionary of numpy arrays or serializable values.
    *args : Any
        Positional arguments passed to ``compute_fn``.
    recompute : bool, default=False
        Force recomputation even if cache file exists.
    verbose : bool, default=True
        Print status messages about cache hits or computation.
    **kwargs : Any
        Keyword arguments passed to ``compute_fn``.

    Returns
    -------
    dict
        Dictionary containing cached or newly computed arrays.
    """
    path = Path(cache_path).resolve()
    recompute_needed = should_recompute(recompute) or not path.is_file()

    if not recompute_needed:
        if verbose:
            print(f"Loading cached example data from {path.name}...")
        raw_data = np.load(path)
        # Convert np.load archive into standard dictionary
        return {k: raw_data[k] for k in raw_data.files}

    if verbose:
        print(f"Computing example data (will save to {path})...")
    data = compute_fn(*args, **kwargs)

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **data)
    if verbose:
        size_kb = os.path.getsize(path) / 1024.0
        print(f"Saved example cache to {path} ({size_kb:.1f} KB)")

    return data
