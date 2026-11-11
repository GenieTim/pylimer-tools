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


def get_example_cache_path(
    filename: str, example_folder: str = "force_balance"
) -> Path:
    """Resolve the filesystem path to an example's cached data file (.npz).

    Searches intelligently from current working directory, caller file location,
    and repository root to find existing data files, or creates the target directory
    structure for new caches.

    Parameters
    ----------
    filename : str
        Name of the cache file (e.g. "entanglement_link_vs_spring_cache.npz").
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
                return (p / "examples" / example_folder /
                        "data" / filename).resolve()
            # Check p/<example_folder>/data/filename
            if (p / example_folder / "data" / filename).is_file():
                return (p / example_folder / "data" / filename).resolve()

    # If not found, resolve default location for creation
    for start_dir in list(candidates):
        cur = start_dir
        for p in [cur] + list(cur.parents):
            if (p / example_folder / "data").is_dir():
                return ((p / example_folder / "data") / filename).resolve()
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
    from scipy import stats

    t_crit = float(stats.t.ppf((1.0 + confidence) / 2.0, df=n - 1))
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


def resolve_max_workers(max_workers: Optional[int] = None) -> int:
    """Resolve the effective number of worker processes for parallel simulation tasks.

    Checks:
    1. Explicit ``max_workers`` argument (if provided and >= 1).
    2. CLI arguments: ``--workers N``, ``-j N``, or ``--workers=N`` in ``sys.argv``.
    3. Environment variables: ``PYLIMER_NUM_WORKERS`` or ``PYLIMER_MAX_WORKERS``.
    4. Default fallback: available CPU cores via ``os.cpu_count() or 1``.

    Parameters
    ----------
    max_workers : int, optional
        Explicit number of workers requested.

    Returns
    -------
    int
        Positive integer for the number of workers to use.
    """
    if max_workers is not None and max_workers >= 1:
        return max_workers

    # Check CLI arguments in sys.argv
    for i, arg in enumerate(sys.argv):
        if arg in ("--workers", "-j") and i + 1 < len(sys.argv):
            try:
                val = int(sys.argv[i + 1])
                if val >= 1:
                    return val
            except ValueError:
                pass
        elif arg.startswith("--workers="):
            try:
                val = int(arg.split("=", 1)[1])
                if val >= 1:
                    return val
            except ValueError:
                pass

    # Check environment variables
    for env_var in ("PYLIMER_NUM_WORKERS", "PYLIMER_MAX_WORKERS"):
        val_str = os.environ.get(env_var, "")
        if val_str:
            try:
                val = int(val_str)
                if val >= 1:
                    return val
            except ValueError:
                pass

    return max(1, os.cpu_count() or 1)


def configure_process_threading(
    workers: Optional[int] = None,
    force_single_thread: bool = False,
    verbose: bool = True,
) -> None:
    """Configure threading environment variables when running in parallel.

    When running simulation workloads across multiple worker processes (e.g. using
    ProcessPoolExecutor), each process running C++ extensions (Eigen, OpenMP,
    OpenBLAS) would otherwise default to using all available CPU cores. This causes
    massive CPU oversubscription (e.g. 12 processes x 12 OpenMP threads = 144 threads
    competing for 12 physical cores), severely degrading performance.

    If workers >= os.cpu_count() (running on all processors) or workers > 1,
    and OMP_NUM_THREADS is not already explicitly set by the user, this function sets
    OMP_NUM_THREADS=1 and related BLAS/LAPACK threading variables to 1.
    """
    if "OMP_NUM_THREADS" not in os.environ or force_single_thread:
        effective_workers = workers if workers is not None else resolve_max_workers()
        cpu_total = os.cpu_count() or 1
        if (
            effective_workers >= cpu_total
            or effective_workers > 1
            or force_single_thread
        ):
            for var in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            ):
                os.environ[var] = "1"
            try:
                import threadpoolctl

                threadpoolctl.threadpool_limits(1)
            except ImportError:
                pass
            if verbose:
                print(
                    f"Configured single-threaded OpenMP/BLAS for {effective_workers} parallel worker(s) (OMP_NUM_THREADS=1)."
                )


def worker_process_initializer() -> None:
    """Initializer for ProcessPoolExecutor worker processes to enforce single-threaded OpenMP/BLAS."""
    for var in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[var] = "1"
    try:
        import threadpoolctl

        threadpoolctl.threadpool_limits(1)
    except ImportError:
        pass


def should_recompute(recompute_flag: bool = False) -> bool:
    """Determine whether cached data should be recomputed.

    Evaluates CLI arguments (``--recompute``), environment variables
    (``PYLIMER_RECOMPUTE_CACHE=1``), and explicit function arguments.
    """
    if recompute_flag:
        return True
    if "--recompute" in sys.argv:
        return True
    if os.environ.get("PYLIMER_RECOMPUTE_CACHE",
                      "").lower() in ("1", "true", "yes"):
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
    configure_process_threading(verbose=verbose)
    data = compute_fn(*args, **kwargs)

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **data)
    if verbose:
        size_kb = os.path.getsize(path) / 1024.0
        print(f"Saved example cache to {path} ({size_kb:.1f} KB)")

    return data
