# Changelog

## main

### Apple Silicon Metal GPU Acceleration
- Add Apple Silicon GPU acceleration via Apple Metal framework (`USE_METAL` option enabled by default on macOS):
  - **Normal Mode Analysis**:
    - Implement `computeAllEigenvaluesMetal` using parallel Jacobi rotations on Metal compute shaders.
    - Implement `computeChebyshevMomentsKPMMetal` for Kernel Polynomial Method (KPM) with custom Metal shaders for CSR sparse matrix-vector multiplication (SpMV), Chebyshev polynomial recurrence, and parallel reduction dot products (delivering over 25x speedup for large networks).
    - Add `NormalModeSolver.METAL` enum value and string option (`"metal"`, `"gpu"`) to `NormalModeAnalyzer.find_all_eigenvalues`.
    - Expose `compute_chebyshev_moments_kpm_metal` in Python bindings and type stubs.
  - **MEHP Network Force Relaxation**:
    - Add GPU-accelerated sparse linear equation solvers in `MetalForceBalanceSolver`:
      - `METAL_CG`: Conjugate Gradient SLE solver using Metal compute shaders.
      - `METAL_BARZILAI_BORWEIN`: Barzilai-Borwein two-point step size gradient descent solver using Metal compute shaders.
    - Add `SLESolver.METAL_CG` and `SLESolver.METAL_BARZILAI_BORWEIN` to `MEHPForceBalance2` and Python bindings.
    - **Performance & Stability Enhancements**:
      - SIMDgroup-cooperative SpMV kernel (`csr_spmv_coop8`) with 8 threads per row, coalesced memory access, and `simd_shuffle_down` reductions.
      - Persistent `MTLBuffer` pooling in `MetalFBSBackend` to eliminate per-solve buffer allocation overhead.
      - GPU-resident batched iteration loops (20-25 steps per batch) evaluating step sizes and convergence on-device without CPU round-trips.
      - Translational zero-mode (center-of-mass) projection to eliminate rigid body drift in unanchored periodic networks.
      - Periodic gradient recomputation in Barzilai-Borwein batches to prevent single-precision drift.
      - Multi-RHS 3D Block-SpMV formulation: represent isotropic spring networks as an $N \times N$ scalar network Laplacian rather than a $3N \times 3N$ sparse matrix, achieving an exact $3\times$ matrix storage and memory bandwidth reduction with fused 3-RHS SpMM kernels (`csr_spmm_3rhs`, `csr_spmm_3rhs_coop8`).
      - GPU-resident Jacobi (diagonal scaling) preconditioning (`bufInvDiag`, `pcg_apply_jacobi_and_dot_3d`), accelerating Conjugate Gradient convergence by 20%–40% on heterogeneous polymer networks.
      - End-to-end GPU spring evaluation (`evaluateSpringStrainAndInactivityMetal`) to evaluate strand overextension and inactivity cutoffs directly on-device using GPU-resident coordinates, skipping outer relaxation loops when no bonds break.
      - Pure SpMV benchmarking utility (`benchmarkSpMVMetal`) and standalone 3D SpMV verification (`spmvMetal3D`).
- Fix `ADD_ITERATIONS` assertion in `MEHPForceBalance2` when residual vector is zero.
- Add comprehensive Catch2 accuracy benchmarks and tests for Metal KPM, Jacobi, Conjugate Gradient, Barzilai-Borwein solvers, buffer pooling, and zero-mode projection with cross-platform build guards.

## v0.4.1

### Build System & PyPI WebAssembly Toolchain
- Upgrade WebAssembly / Pyodide build toolchain to Pyodide 314 (`pyemscripten_2026_0_wasm32`), Python 3.14, and Emscripten 5.0.3 per PEP 783.
- Publish `pyemscripten` WebAssembly wheels directly to PyPI alongside native wheels and sdist, enabling direct in-browser installation (`micropip.install("pylimer-tools")`).
- Pass WebAssembly exception handling (`-fwasm-exceptions`) and forward toolchain configurations to vendor dependencies.
- Remove obsolete downstream `pyodide-recipes` sync automation now that PyPI directly hosts WebAssembly wheels.

### Normal Mode Analysis & Rheology
- Add `NormalModeSolver` enum (`AUTO`, `DENSE_LAPACK`, `DENSE_EIGEN`, `SPARSE_SLICING`) and solver string selection to `NormalModeAnalyzer.find_all_eigenvalues`.
- Add connected component decomposition (`decompose_components`) with analytical solutions for simple linear Rouse chains and ring polymers.
- Implement sparse slicing eigenvalue solver (`SPARSE_SLICING`) based on Sylvester's Law of Inertia with `Eigen::SimplicialLDLT` and `Spectra::SymEigsShiftSolver`.
- Add Kernel Polynomial Method (KPM) for $O(N)$ spectral calculations without full matrix diagonalization:
  - `compute_chebyshev_moments_kpm`: Chebyshev moment computation with Rademacher stochastic trace estimation and Jackson damping kernel.
  - `compute_density_of_states_kpm`: continuous density of states (DOS) reconstruction.
  - `evaluate_stress_autocorrelation_kpm`: direct evaluation of stress autocorrelation $G(t)$ via Gauss-Chebyshev quadrature.
  - `evaluate_storage_modulus_kpm`: direct evaluation of storage modulus $G'(\omega)$ via Gauss-Chebyshev quadrature.
  - `evaluate_loss_modulus_kpm`: direct evaluation of loss modulus $G''(\omega)$ via Gauss-Chebyshev quadrature.
- Add comprehensive unit tests in Python and C++ covering all normal mode solvers, multigraph parallel bonds, component decomposition, and KPM viscoelastic moduli.

## v0.4.0

### Simulation & Physics
- Fix inverse Langevin approximation (`langevin_inv`) symmetry for negative displacements, ensuring restorative forces have the correct sign.
- Fix stress evaluation in `NonGaussianSpringForceEvaluator` to properly account for non-linear Langevin chain tension.
- Prevent division by zero and NaN in `evaluateForceSetGradient` and `csch` when spring distance is zero.
- Ensure safe move semantics and delete copy operations for `MEHPForceRelaxation` to prevent solver handle double-free.
- Fix periodic boundary validation in `Box::validateOffset` for negative displacement offsets.
- Prevent mutating input vectors in `median(Eigen::DenseBase)` and improve its performance with `std::nth_element`.
- Support empty networks and direct crosslinker bonds in `compute_topological_factor`.
- Optimize DataFrame memory reduction utility to only downcast numeric column types.

### Build System & Performance
- Performance: switch to divide-and-conquer `LAPACKE_dsyevd` in `NormalModeAnalyzer` when computing both eigenvalues and eigenvectors, significantly speeding up dense decomposition.
- Feature: add GPU-accelerated dense eigenvalue calculation via cuSOLVER (`computeAllEigenvaluesGPU`) with fallback to CPU if compiled without CUDA.
- Improve vendor dependency builds (`nlopt`, `cereal`, `igraph`, `spectra`) with forwarded architecture, sysroot, and compiler flags (`vendorForwardArgs.cmake`).
- Remove `-ffast-math` from high-performance compiler flags to avoid NaN and precision inconsistencies.
- Detect Ninja build generator via `ninja.BIN_DIR` and fix parallel build core allocation in `setup.py`.
- Link standard filesystem library conditionally via `STD_FS_LIB` rather than hardcoding `stdc++fs`.
- Clean up build-time requirements in `pyproject.toml`.

### Documentation & Examples
- Document missing modules: `analyse_networks`, `basic_lammps_structure_stats`, `displace_randomly`, `generate_network`, `bead_spring_parameter_provider`, `polymer_data`, and `stoichiometric_relationships`.
- Fix Sphinx configuration and nitpicky reference warnings; update intersphinx to include `scipy`.
- Standardize fixture path discovery across all example scripts (`examples/**/plot_*.py`) for standalone runs and Sphinx Gallery generation.
- Add comprehensive unit tests and regression tests for MEHP simulation, box offsets, and Eigen utilities.

## v0.3.16-v0.3.19

- New attempts to fix webassembly deployment
- Deployment improvements, in particular for manylinux

## v0.3.15

- Fix an issue (#8) where the bond_types could be in an invalid state.
- Improve compatibility with newer versions of pandas

## v0.3.14

Thank you to @jasonmulderrig for the contribution, see PR #7:

- Add the characteristic ratio to the bead-spring parameter provider 
- Fix some issues in the documentation, volumes being computed incorrectly in some examples

## v0.3.13

- Implement new utility methods in `pylimer_tools.calc.stoichiometric_relationships` in particular to help with the parameter `b2`
- Restructure Python tests a bit to better reflect the structure of the Python code
- Finally re-implement possibility to use Python functions for determining the forces in `MEHPForceRelaxation`
- Migrate to pybind11 v3

## v0.3.12

- More build system improvements

## v0.3.11

- Improve performance of webassembly by compiling SIMD
- (Possibly) improve deployment

## v0.3.10

- Documentation improvements
- Upgrade igraph to 1.0
- Fix an issue in `predict_p_from_w_sol`, which returned incorrect results for `w_sol = 0.`

## v0.3.7 - v0.3.9

- Build system improvements (include missing files in sdist) with objective to get [pyodide recipe](https://github.com/pyodide/pyodide-recipes/pull/349) working

## v0.3.7

- Build system improvements (incl. missing files in tarball)
- Improve the documentation of the MCUniverseGenerator
- Fix an issue where serialization to text files might not work

## v0.3.6

- Next attempt to improve CI (Windows did not properly read the version from git)
- Fix an issue where certain unit styles could not be queried anymore

## v0.3.5

- More build & deploy system improvements (we got wasm now!)
- Less expensive imports, use lazy importing of scipy
- Add method to add regularly spaced cross-linker chains in `MCUniverseGenerator`

## v0.3.4

- Build & deploy system improvements

## v0.3.3

- Enable switching between the unit types in the `generate_network` command
- (Possibly) fix compilation if Eigen3 is not installed

## v0.3.2

- Add methods to the `MCUniverseGenerator` to spawn star-like crosslinkers (`add_star_crosslinkers`) as well as link those to each other (`link_strands_to_strands_to_conversion`)
- Implement Python binding for the `MCUniverseGenerator`'s `link_strands_callback`, i.e., allow the Python users to stear and stop the cross-linking procedure
- Switch from `MEHPForceRelaxation` to `MEHPForceBalance2` internally of the `MCUniverseGenerator` to improve performance
- Allow `MEHPForceRelaxation` and `MEHPForceBalance2` to compute soluble fractions even when the Universe is empty

## v0.3.1

This version mainly fixes some issues in the CI.

## v0.3.0

This has become the new baseline version, since it's the first public version since the very first try to publish many years ago.
There are only breaking changes, compared to previous versions, such that a CHANGELOG has become useless.
Sorry for that, please enjoy the many new features by taking a look at the revised documentation.
