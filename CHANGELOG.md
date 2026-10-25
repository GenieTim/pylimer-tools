# Changelog

## main

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
