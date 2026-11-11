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
The stoichiometric imbalance $r$ and the degree of polymerization $p$ (see also docs/nomenclature.rst) can vary widely.
