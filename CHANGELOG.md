# Changelog

## main

### OpenMP Parallelization, Thread-Safety & Build Optimizations
- **OpenMP Configuration & macOS `libomp` Discovery (`vendor/openmp.CMakeLists.cmake`)**:
  - Automatically detect Homebrew `libomp` on macOS with Clang/AppleClang (`OpenMP_ROOT`, `brew --prefix libomp`, `/opt/homebrew/opt/libomp`, `/usr/local/opt/libomp`) and unify OpenMP configuration across `pylimer_tools` and `test_sources`.
  - Include `<omp.h>` in `ExtraEigenTypes.h` when `OPENMP_FOUND` is enabled and add Catch2 unit tests for OpenMP `Eigen::VectorXd` and `Eigen::Matrix3d` parallel reductions.
- **Parallelized Algorithms & Thread-Safety**:
  - Parallelize `NormalModeAnalyzer::computeChebyshevMomentsKPM` across Rademacher trace estimation vectors with deterministic pre-generated random vectors and OpenMP `Eigen::VectorXd` reductions.
  - Parallelize `ResistanceDistanceAnalysis::computeAllPairsEffectiveResistance` and `computeKirchhoffIndex`.
  - Parallelize `RouseBrownianSimulator::computeMSD`, `computeCenterOfMassMSD`, and `computeEndToEndAutocorrelation` over lag times $\tau$.
  - Parallelize `Universe::countAtomsInSkinDistance` and `UniverseSequence::computeDistanceFromToAtoms` / `computeVectorFromToAtoms`.
  - Add thread-private PRNGs (`threadRNGs`, `initThreadRNGs`) and persistent thread count configuration (`setNumOmpThreads`, `getNumOmpThreads`) in `DPDSimulator` to eliminate data races on `std::mt19937` during multi-threaded pair-force evaluation.
  - Fix missing `totalMass` reduction clause in `Molecule::computeRadiusOfGyration` and add loop size guards (`if (...)`) to avoid OpenMP overhead on small inputs.
- **Build System & CI Workflows**:
  - Propagate `EIGEN3_INCLUDE_DIRS` from `Eigen3::Eigen` target properties in `vendor/eigen.CMakeLists.cmake` and link `Eigen3::Eigen` directly in `tests/CMakeLists.txt`.
  - Optimize GitHub Actions workflows with path triggers, `ccache`, Ninja build generator, and conditional C++ test execution (`RUN_CPP_TESTS`) when only Python files are modified.

### Analytical Brownian-Bridge Entanglement Sampling & Universe-Free Force Balance
- **Analytical Entanglement Detection (`EntanglementDetector`)**:
  - Implement `randomlyFindEntanglementsAnalytical` (`randomly_sample_entanglements_analytical`) to probabilistically sample spatial slip-links directly on a `CrosslinkerUniverse` using conditional Brownian-bridge and random-walk Gaussian distributions between crosslinker anchors without instantiating explicit 3D bead coordinates or constructing an intermediate `Universe`.
  - Add topological leaf peeling to `randomlyFindEntanglementsAnalytical` (`filter_dangling_and_soluble=True`) to iteratively prune dangling chain ends, soluble components, and degree-$\le 1$ dead ends, ensuring slip-links are placed exclusively on load-bearing strands.
  - Introduce `StrandPairEntanglements` structure to store sampled strand indices, contour fractions ($\alpha_1, \alpha_2$), bead indices, and 3D coordinates, with `toAtomPairEntanglements` mapping directly to discrete atom pairs.
- **Universe-Free Network Mapping (`MCUniverseGenerator`)**:
  - Add `MCUniverseGenerator::sampleEntanglements` to sample entanglements directly from the generator's internal `simplifiedUniverse` (`CrosslinkerUniverse`).
  - Add `MCUniverseGenerator::getForceBalance` and `MCUniverseGenerator::getForceBalance2` overloads taking `StrandPairEntanglements` or sampling parameters directly, bypassing full bead generation and reducing peak resident memory (RSS) by 80% to over 90% on large networks.
  - Add `MCUniverseGenerator::convertToForceBalanceNetwork` and `MCUniverseGenerator::convertToForceBalance2Network` supporting both merged static entanglement links and split entanglement springs.
  - Expose `MCUniverseGenerator::getCrosslinkerUniverse` and `MCUniverseGenerator::getBox` in Python bindings.
- **Universe-Independent Force Balance Solvers (`MEHPForceBalance` & `MEHPForceBalance2`)**:
  - Extend `MEHPForceBalance::randomlyAddSliplinks` to support analytical entanglement sampling when initialized without a full `Universe` (`universe.getNrOfAtoms() == 0`).
  - Implement `MEHPForceBalance2::extractNetworkMetadata` to infer total atom count, cluster topology, and active/soluble weight fractions directly from network connectivity when initialized without a `Universe`.
  - Expose `get_active_weight_fraction` in Python bindings for both `MEHPForceBalance` and `MEHPForceBalance2`.
  - Expose `Universe::clearCaches` (`clear_caches`) in Python bindings to reclaim decomposed molecule and cluster cache memory.
- **Documentation & Benchmarks**:
  - Add comprehensive Sphinx-gallery example `examples/force_balance/plot_entanglement_sampling_comparison.py` comparing memory scaling, runtime, contour distributions, and network shear modulus across conversions between discrete and analytical methods, with accompanying precomputed benchmark cache and generator script.


### Automatic Solver Selection Optimizations
- **NormalModeAnalyzer**:
  - Fix `NormalModeSolver::AUTO` to only select `SPARSE_SLICING` when `!includeEigenvectors && n > 5000` (previously selected sparse slicing for $n > 2000$ even when eigenvectors were requested, triggering redundant eigenvalue bisection followed by dense fallback diagonalization).
  - Raise dense-to-sparse eigenvalue threshold to $n = 5000$ where dense LAPACK (`LAPACKE_dsyevr` MRRR) is up to ~29x faster than sparse slicing and fits in < 200 MB RAM, while preserving sparse slicing for large scales ($N > 55,000$, where dense matrix allocation requires $> 24$ GB RAM).
  - Automatically accelerate `computeEffectiveSpectrumKPM` and `computeDensityOfStatesKPM` via `computeChebyshevMomentsKPMMetal` when Metal GPU acceleration is available (> 25x speedup for large systems), falling back seamlessly to CPU when unavailable.
- **MEHPForceBalance2**:
  - Always default `SLESolver::DEFAULT` and `SLESolver::AUTO` to `SLESolver::CONJUGATE_GRADIENT` with automatic non-convergence and breakdown fallback to `BiCGSTAB`, ensuring robust and optimal iterative performance across all network topologies (including networks with more links than nodes).
  - Add `SLESolver::AUTO` enum (`"auto"`) and Python string aliases (`"auto"`, `"default"`, `"cg"`).
  - Add adaptive solver fallback from `ConjugateGradient` to `BiCGSTAB` when CG reaches a configurable iteration threshold without converging (`configCGFallbackMaxIterations` / `getCGFallbackMaxIterations`, `wasCGFallbackTriggered`, `getEffectiveSolverChoice`, default `20000` iterations), and increase the default `max_iterations` in Python `run_force_relaxation` to `50000`.
  - Replace full `Universe` member storage in `MEHPForceBalance2` with lightweight pre-extracted atom and cluster metadata (`extractUniverseMetadata`, `getNrOfAtomsInUniverse`) and `const Universe&` `pybind11` constructor bindings, eliminating deep `Universe` copies on instantiation and repeated `getClusters()` decompositions in `countActiveClusteredAtoms`.
  - Allow `sle_solver` in Python `run_force_relaxation` to accept both `SLESolver` enum instances and string descriptors (`"auto"`, `"metal_cg"`, `"cg"`, etc.).
  - Fix `exitReason` handling in `METAL_CG` and `METAL_BARZILAI_BORWEIN` when constant vector is zero to correctly record `ExitReason::X_TOLERANCE`.

### Apple Silicon Metal GPU Acceleration
- Add Apple Silicon GPU acceleration via Apple Metal framework (`USE_METAL` option enabled by default on macOS):
  - **Normal Mode Analysis**:
    - Implement `computeAllEigenvaluesMetal` using parallel Jacobi rotations on Metal compute shaders.
    - Implement `computeChebyshevMomentsKPMMetal` for Kernel Polynomial Method (KPM) with custom Metal shaders for CSR sparse matrix-vector multiplication (SpMV), Chebyshev polynomial recurrence, and parallel reduction dot products (delivering over 25x speedup for large networks).
    - Add multi-RHS block-batched KPM Metal kernels (`batch_size = 4` with `float4` and `batch_size = 8` with `Float8`) to `computeChebyshevMomentsKPMMetal` to amortize sparse matrix memory reads across multiple Rademacher trace estimation vectors.
    - Add `NormalModeSolver.METAL` enum value and string option (`"metal"`, `"gpu"`) to `NormalModeAnalyzer.find_all_eigenvalues`.
    - Expose `compute_chebyshev_moments_kpm_metal` in Python bindings and type stubs.
  - **MEHP Network Force Relaxation**:
    - Add GPU-accelerated sparse linear equation solvers in `MetalForceBalanceSolver`:
      - `METAL_CG`: Conjugate Gradient SLE solver using Metal compute shaders.
      - `METAL_BARZILAI_BORWEIN`: Barzilai-Borwein two-point step size gradient descent solver using Metal compute shaders.
    - Add `SLESolver.METAL_CG` and `SLESolver.METAL_BARZILAI_BORWEIN` to `MEHPForceBalance2` and Python bindings.
    - **Performance & Stability Enhancements**:
      - SIMDgroup-cooperative SpMV kernel (`csr_spmv_coop8`) with 8 threads per row, coalesced memory access, and `simd_shuffle_down` reductions.
      - Adaptive SIMDgroup-cooperative 3-RHS SpMM kernels with 2, 4, or 8 threads per row (`csr_spmm_3rhs_coop2`, `csr_spmm_3rhs_coop4`, `csr_spmm_3rhs_coop8`, `spmvMetal3DCoop`, and `benchmarkCoopWidthsMetal`), automatically selecting the optimal width based on average node degree.
      - Degree-aware hybrid SIMDgroup SpMM work-queue (`spmvMetal3DHybrid`, `benchmarkHybridSpMVMetal`, `csr_spmm_3rhs_hybrid_low_coop2`, `csr_spmm_3rhs_hybrid_high_coop8`, `csr_spmm_3rhs_hybrid_high_coop16_stateless`) dispatching low-degree chain beads (`degree <= 3`) to `coop2` and high-degree crosslinkers to `coop8`/`coop16`.
      - Persistent `MTLBuffer` pooling in `MetalFBSBackend` to eliminate per-solve buffer allocation overhead.
      - GPU-resident batched iteration loops (20-25 steps per batch) evaluating step sizes and convergence on-device without CPU round-trips.
      - Translational zero-mode (center-of-mass) projection and algebraically pinned reference node solver (`solveSLEMetalCG3DPinned`, `pinnedNode` in `solveSLEMetalCG3D`) to eliminate rigid body drift in unanchored periodic networks and improve conditioning without COM projection passes.
      - Periodic gradient recomputation in Barzilai-Borwein batches, mixed-precision iterative refinement in `solveSLEMetalCG3D` for tight tolerances (`< 1e-5`), and automatic CPU double-precision fallback (`enableCpuFallback` in `solveSLEMetalCG` / `solveSLEMetalCG3D`) for ill-conditioned networks (`tolerance < 1e-6`).
      - Multi-RHS 3D Block-SpMV formulation: represent isotropic spring networks as an $N \times N$ scalar network Laplacian rather than a $3N \times 3N$ sparse matrix, achieving an exact $3\times$ matrix storage and memory bandwidth reduction with fused 3-RHS SpMM kernels (`csr_spmm_3rhs`, `csr_spmm_3rhs_coop8`), and extend the scalar $N \times N$ + 3-RHS formulation to CPU direct and iterative solvers in `MEHPForceBalance2`.
      - In-place CSR system matrix and RHS updates when strands break (`breakTooLongStrands` with `BrokenSpringRecord`) under `NO_SIMPLIFICATION`, with a configurable cumulative broken-spring rebuild threshold (`configInPlaceRebuildThreshold` / `getInPlaceRebuildThreshold`, default 15%) to trigger a clean CSR rebuild and symbolic re-factorization once dead structural zeros accumulate.
      - GPU-resident Jacobi (diagonal scaling) preconditioning (`bufInvDiag`, `pcg_apply_jacobi_and_dot_3d`), accelerating Conjugate Gradient convergence by 20%–40% on heterogeneous polymer networks.
      - Zero-copy end-to-end GPU relaxation step (`runGpuRelaxationStep`, `evaluateSpringStrainAndInactivityMetal`, `getGpuResidentCoordinates`) combining the GPU SLE solve and spring strain/inactivity evaluation on GPU-resident coordinates, skipping outer relaxation loops when no bonds break.
      - Pure SpMV benchmarking utility (`benchmarkSpMVMetal`) and standalone 3D SpMV verification (`spmvMetal3D`).
- Fix `ADD_ITERATIONS` assertion in `MEHPForceBalance2` when residual vector is zero.
- Add comprehensive Catch2 accuracy benchmarks and tests for Metal KPM, Jacobi, Conjugate Gradient, Barzilai-Borwein solvers, buffer pooling, and zero-mode projection with cross-platform build guards.

### Topological Resistance Distance & Rouse Brownian Dynamics
- Add `ResistanceDistanceAnalysis` (`pylimer_tools_cpp.ResistanceDistanceAnalysis`) for graph-Laplacian topological analysis:
  - Single-pair, batched multi-RHS (`computeEffectiveResistanceBatched`), and all-pairs effective resistance distance calculations with disconnected component support.
  - Level-3 BLAS (`dtrsm`) multi-column chunked triangular solves (`solveIdentityChunked` / `solve_identity_chunked` and `chunk_size` parameter in `computeAllPairsEffectiveResistance` and `computeMoorePenrosePseudoInverse`, default `128`).
  - Kirchhoff index (direct and spectral), Moore-Penrose Laplacian pseudo-inverse, node resistance centrality, and random-walk commute times.
  - Canonical topology factories (`createLinearChain`, `createRing`, `createComplete`, `createStar`).
- Add `RouseBrownianSimulator` (`pylimer_tools_cpp.RouseBrownianSimulator`, `RouseSimulationResult`, `RouseSolverType`) for implicit backward-Euler Rouse and Brownian dynamics simulations of polymer chains and spring networks:
  - Single scalar shifted-Laplacian factorization ($A = I + \Delta t (k/\gamma) L$) reused across all 3 spatial RHS columns (`DIRECT_LDLT`, `DIRECT_LLT`, `CONJUGATE_GRADIENT`).
  - Support for thermal noise, deterministic relaxation, external forces, and trajectory observables (bead MSD, center-of-mass MSD, end-to-end vector autocorrelation, radius of gyration, and potential energy).

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
