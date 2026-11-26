# Agent Guidelines for `pylimer-tools`

## Python Environment

Whenever Python, `pip`, `unittest`, or Python-based linters/formatters are needed, activate the project virtual environment first:

```bash
source /Users/timbernhard/Privat/Programming/Studium/Doctorate/venv/bin/activate
```

## Building & Running Tests

- **Only run tests when making changes**: Do not run the test suite for review-only or git-only tasks unless changes were actually made.
- **Exclude `[long]` tests and benchmarks by default**: When running C++ Catch2 tests (`pylimer_tests`), always exclude `[long]` tagged tests and benchmarks unless your changes specifically affect them:
  ```bash
  ./build-debug/tests/pylimer_tests "~[long]" --skip-benchmarks
  ```
  - For faster iteration, filter by specific test name or tag when working on a particular module, then run the broader non-`[long]`, non-benchmark suite.
  - Build C++ debug/test targets via `./bin/build-fast.sh debug` (builds into `build-debug/`) or `cmake --build build-debug --parallel`.
- **Python tests**: After activating the virtualenv (and rebuilding the extension if C++ bindings changed), run:
  ```bash
  python -m unittest discover -v
  ```

## Code Formatting & Pre-Commit Hook

- **Automatic formatting on commit**: `.git/hooks/pre-commit` runs `./bin/format-code.sh` (or `pre-commit`), which formats C++ files (`.cpp`, `.h` in `src/pylimer_tools_cpp` and `tests/pylimer_tools_cpp` via `clang-format --style=file --fallback-style="Mozilla"`) and Python files (`ruff format` / `autopep8`) **in-place**.
- **Important**: Because `./bin/format-code.sh` modifies files in the working tree without re-staging them, running `git commit` on unformatted code will leave formatting diffs unstaged after the commit completes.
  - **Recommended workflow**: Run `./bin/format-code.sh` _before_ staging (`git add -u`), or check `git status -s -uno` immediately after `git commit` and fold any hook-formatted files into the commit with `git add -u && git commit --amend --no-edit`.

## Code Style & Implementation Guidelines

- **Use enums instead of raw strings**:
  - Prefer strongly-typed enums (`enum class` in C++, `enum.Enum` in Python) over raw strings, magic numbers, or boolean flags for modes, options, solver types, algorithms, and states.
  - When exposing C++ enums to Python, bind them via pybind11 (`py::enum_<...>`) to ensure consistent, type-safe usage across both languages.
  - If string parsing or serialization is needed (e.g., CLI arguments, file formats, or config inputs), isolate string-to-enum and enum-to-string mappings in dedicated helper functions/converters rather than branching on raw strings throughout core logic.
- **Prefer explicit types over `auto`**:
  - In C++, write explicit types (e.g., `size_t`, `double`, `std::vector<int>`, `Eigen::VectorXd`, `const Molecule&`) instead of `auto`.
  - Reserve `auto` strictly for cases where:
    - The type is explicitly repeated on the right-hand side (e.g., `auto ptr = std::make_shared<Universe>(...)`, `auto x = static_cast<size_t>(...)`).
    - The type is complex, verbose, or unnameable (e.g., complex iterator types, lambda signatures).
  - Explicit types improve readability, prevent unintended implicit conversions, and avoid accidentally dropping `const` or reference qualifiers.
  - In Python, provide comprehensive type annotations for all function signatures (arguments and return types) and class attributes.
- **Const-correctness**:
  - In C++, mark member functions, local variables, and parameters `const` by default whenever they are not modified.
  - Pass non-primitive parameters by `const &` (or `std::string_view` / `std::span` when appropriate) unless passing by value or transferring ownership.
- **Error handling & validation**:
  - Prefer throwing standard C++ exceptions (`std::invalid_argument`, `std::runtime_error`, `std::out_of_range`) with descriptive messages over returning error codes or silent failure; pybind11 automatically translates these into the corresponding Python exceptions.
  - Validate function arguments and invariants early.

## Utility Scripts (`bin/`)

The `bin/` directory contains standard helper scripts for builds, testing, formatting, and type stub generation:

- **Type stubs generation (`bin/make-stubs.sh`)**:
  - Automatically regenerates `src/pylimer_tools_cpp-stubs/__init__.pyi` using `pybind11-stubgen`.
  - **Never edit stubs by hand**: Run `./bin/make-stubs.sh` (or activate the virtualenv and run `pybind11-stubgen pylimer_tools_cpp -o src/pylimer_tools_cpp-stubs --numpy-array-remove-parameters && mv src/pylimer_tools_cpp-stubs/pylimer_tools_cpp.pyi src/pylimer_tools_cpp-stubs/__init__.pyi`) after modifying C++ bindings.
- **Code formatting (`bin/format-code.sh`)**:
  - In-place formatting of C++ (`clang-format`) and Python (`ruff format` / `autopep8`).
  - Always run before staging (`git add -u`).
- **Fast build (`bin/build-fast.sh [preset] [clean]`)**:
  - Builds CMake presets (`default`, `debug`, `high-performance`, `testing`) with optimal parallel jobs into `build-<preset>/`.
- **Test runners (`bin/run-tests-short.sh` / `bin/run-tests.sh`)**:
  - Convenient wrappers to run unit tests and check output.
- **Documentation (`bin/make-docs.sh`)**:
  - Builds Sphinx documentation.
- **Packaging (`bin/build-wheel.sh` / `bin/build-sdist.sh`)**:
  - Builds distribution wheels and source distributions.

## Git Workflow (Staging, Pulling & Pushing)

1. **Stage selectively**: The repository root often contains untracked local scripts, test outputs, or build directories (`build-debug/`, scratch `.py`/`.gml` files, etc.). Use `git add -u` or specify exact file paths—never blindly run `git add .` or `git add -A`.
2. **Ensure a clean tracked working tree before pulling**: Make sure no unstaged formatting changes remain from the pre-commit hook (`git status -s -uno` should be empty).
3. **Pull with rebase before pushing**: CI workflows automatically commit and push updated coverage badges (`Updated coverage badges [skip ci]`) to `main`. If `git push` is rejected as non-fast-forward, run:
   ```bash
   git pull --rebase
   git push
   ```

## System Sizes

This library is designed for realistic polymer simulations, which typically involve large systems with many chains and entanglements. When investigating performance, target system sizes should be on the order of approx 1e5 chains, each with 1e2 to 1e3 beads.
The stoichiometric imbalance $r$ and the degree of polymerization $p$ (see also `docs/nomenclature.rst`) can vary widely.

## Dangers of This Codebase

Polymer network topology and mechanics contain subtle physical and algorithmic edge cases. Read `docs/nomenclature.rst` and `docs/assumptions.rst` before modifying topology, network construction, or force-balance code, and always respect the following invariants:

### 1. Periodic Boundary Conditions (PBC) & `assume_box_large_enough`

PBC handling is one of the easiest places to introduce silent, catastrophic physics bugs (see `docs/assumptions.rst`):

- **Coarse-grained springs can exceed half the box length ($L/2$)**: In bead-spring MD, individual covalent bonds are much shorter than $L/2$. In `MEHPNetwork` / `MEHPSimulator` (and coarse-grained `DPDSimulator` setups), a single spring represents an entire polymer strand (or sub-strand) whose end-to-end vector can legitimately exceed $L/2$ or wrap multiple times across the periodic box.
- **Two PBC modes (`assume_box_large_enough`)**:
  - `assume_box_large_enough = False` (**default**, tracked periodic offsets): Walks the initial bead-by-bead contour to bake per-spring periodic box offsets (`springBoxOffset`, still assuming initial *atomistic* bonds are shorter than $L/2$). Strands longer than $L/2$ and PBC-wrapping loops stay stretched instead of collapsing. Any topology change (inserting, removing, or merging springs; sliding or swapping slip-links; flipping spring endpoint order `A <-> B`) or box deformation (`Box::adjustCoordinatesTo`) **must** preserve, split, negate, or scale `springBoxOffset` consistently.
  - `assume_box_large_enough = True` (minimum-image convention): Applies `Box::handlePBC` directly to endpoint coordinate differences on every evaluation, clamping vectors to $[-L/2, L/2]$ per axis. Any spring exceeding $L/2$ flips direction across the boundary, causing box-spanning strands and PBC-wrapping self-loops to collapse. Watch out for the API default asymmetry: constructors default to `assume_box_large_enough=False`, whereas calling `config_assume_box_large_enough()` with no arguments defaults its `assume` parameter to `True`.
- **Wrapped coordinates, image flags, and simple shear**: `Atom` stores wrapped coordinates alongside integer image flags (`nx`, `ny`, `nz`)—be careful to distinguish wrapped (`getCoordinates()`, `vectorTo()`) from unwrapped (`getUnwrappedCoordinates(box)`, `vectorToUnwrapped()`) queries. `Box` also supports simple shear (`isSheared()`, un-shearing coordinates before minimum-image wrapping in `Box::handlePBC`), though general triclinic tilt is not supported.

### 2. Primary Loops (Self-Loops) & Topological vs. Elastic Activity

- **Primary loops are NOT universally elastically inactive**: A primary loop (a strand starting and ending at the *same* crosslink) relaxes to zero tension only if it closes locally without entanglements. It becomes **elastically effective** in two critical situations:
  1. **Entangled primary loops**: When an entanglement link, entanglement spring, or slip-link is associated with the loop, coupling it to surrounding strands so it carries tension.
  2. **PBC-spanning self-loops**: Under `assume_box_large_enough = False`, even an entanglement-free self-loop (or crosslinker-less ring) that wraps through the periodic box carries a non-zero `springBoxOffset` and cannot collapse to zero length, contributing directly to network elasticity.
- **Never blindly strip self-loops or multi-edges**: Calling `Universe::simplify` (`igraph_simplify`) or passing `remove_primary_loops = True` deletes primary loops (and in `MEHPNetwork`, also removes any entanglement springs attached to them). Never assume self-loops or parallel strands (secondary loops) are disposable graph artifacts.
- **Purely topological pruning vs. PBC-aware relaxation**: `peel_topological_leaves` and `identify_active_atoms` perform 2-core peeling strictly on crosslinker connectivity (`strand_from` / `strand_to`), ignoring coordinates and PBC completely. Consequently, they intentionally **keep** all primary and multi-crosslinker loops (even local unentangled ones, because entanglements could activate them) and **prune** crosslinker-less solvent rings (even PBC-spanning ones). Distinguishing local zero-tension loops from PBC-spanning or entangled active loops requires PBC-aware force relaxation (`MEHPSimulator` with `assume_box_large_enough = False`).

### 3. Entanglements, Slip-Links, Strands & Springs

Several distinct physical models share similar names—follow `docs/nomenclature.rst` strictly:

- **Entanglement representation (`EntanglementMode`)**: `NONE` (phantom crosslink-only network), `MERGED_LINKS` (**entanglement links**: a single 4-functional merged node shared directly on two chains' contours), and `SPRING` (**entanglement springs**: chains remain separate and are coupled by a compliant Hookean spring flagged via `springIsEntanglement`).
- **Contour mobility (`SlipMode`)**: `FIXED` (contour fractions $\alpha$ are frozen—used by `LinearSystemMEHPSolver` and `NonlinearOptimizationSolver`) vs. `SLIDING` (**slip-links / slip-springs** in `SlipLinkRelaxationSolver`, where links slide along chain contours to equalize tension, redistributing `springPartitions` / `springContourLength` and updating `springBoxOffset` when links cross).
- **Strand vs. Spring & link indexing in `MEHPNetwork`**: A **strand** is a full polymer chain between chemical crosslinks (or an entanglement coupling strand in `SPRING` mode), whereas a **spring** is an individual sub-segment between two adjacent links (crosslinks or entanglement links) along that strand. All links (chemical crosslinks, chain ends, and entanglement links) live in `[0, nrOfLinks)`—always check `linkIsEntanglement` rather than assuming crosslinks and entanglements occupy separate index ranges.

### 4. Dual Indexing Schemes & Graph Assumptions

- **Atom IDs ("bonds") vs. Graph Vertex IDs ("edges")**: As detailed in `docs/nomenclature.rst`, methods referring to **bonds** or **atom IDs** use LAMMPS atom IDs (`Atom::getId()`, typically 1-based and potentially non-contiguous), whereas methods referring to **vertices** or **edges** use `igraph`'s contiguous `0`-based `[0, N)` indices (which shift whenever vertices or edges are deleted). Mixing the two silently corrupts or crashes lookups.
- **Bifunctional strands & cache invalidation**: Strand-decomposition routines assume interior chain beads have degree $\le 2$ and expect a single crosslinker atom type per call (`docs/assumptions.rst`). When mutating the underlying `igraph_t` of `Universe` / `Molecule` (`AtomGraphParent`), always use the `GraphMutator` RAII guard so cached properties are invalidated.
