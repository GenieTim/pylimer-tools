# pylimer-tools

[![Run Tests](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml/badge.svg)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)
[![Test Coverage (Python)](https://github.com/GenieTim/pylimer-tools/blob/main/.github/coverage.svg?raw=true)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)
[![Test Coverage (C++)](https://github.com/GenieTim/pylimer-tools/blob/main/.github/cpp-coverage.svg?raw=true)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)
[![Total Coverage](https://codecov.io/gh/GenieTim/pylimer-tools/branch/main/graph/badge.svg?token=5ZE1VSDXJQ)](https://codecov.io/gh/GenieTim/pylimer-tools)
[![Docs](https://github.com/GenieTim/pylimer-tools/actions/workflows/publish-documentation-html.yml/badge.svg)](https://github.com/GenieTim/pylimer-tools/actions/workflows/publish-documentation-html.yml)
[![PyPI version](https://badge.fury.io/py/pylimer-tools.svg)](https://pypi.org/project/pylimer-tools/)
[![Downloads](https://img.shields.io/pypi/dm/pylimer-tools.svg)](https://pypi.org/project/pylimer-tools/)
[![License](https://img.shields.io/pypi/l/pylimer-tools.svg)](LICENSE)

Pronunciation: "pylimer” like "pü-limer”, with the "py” as in the German word "müde” (IPA: /ˈpyːlɪmɚ/).

`pylimer-tools` is a toolkit for simulation, analysis, and data handling of bead–spring polymer systems in Python and C++. It combines high-level Python utilities with performant C++ extensions for common tasks in computational polymer science.

## Table of Contents

- [pylimer-tools](#pylimer-tools)
  - [Table of Contents](#table-of-contents)
  - [1. Features](#1-features)
  - [2. Installation](#2-installation)
    - [Requirements](#requirements)
    - [Stable Release](#stable-release)
    - [Build from Source](#build-from-source)
      - [Requirements](#requirements-1)
      - [Build](#build)
  - [3. Test Installation](#3-test-installation)
  - [4. CLI Tools](#4-cli-tools)
  - [5. Documentation](#5-documentation)
  - [6. Development \& Testing](#6-development--testing)
    - [Adding Features](#adding-features)
  - [7. Citing](#7-citing)
  - [8. Contributing](#8-contributing)
  - [9. Code of Conduct](#9-code-of-conduct)
  - [10. Acknowledgements](#10-acknowledgements)
  - [11. License](#11-license)
  - [Pronunciation Note](#pronunciation-note)

## 1. Features

A selection of features includes:

- Monte Carlo structure (network "universe”) generation
- Dissipative Particle Dynamics (DPD) simulation with slip-springs
- Maximum Entropy Homogenization Procedure (MEHP) with and without slip-links
- LAMMPS output readers: data, dump, thermodynamic outputs
- Network analysis: loops, chain reconstruction, degree statistics
- Structural metrics: radius of gyration, end-to-end distance, distributions
- Normal mode analysis for stress autocorrelation, loss and storage modulus
- Command line interfaces for batch workflows

## 2. Installation

### Requirements

Python >= 3.9.

### Stable Release

Stable release from PyPI:

```
pip install pylimer-tools
```

### Build from Source

#### Requirements

Build requires CMake, a C++17 compiler, and (optionally) Ninja for faster builds.

Additionally, the system packages `flex` and `bison` are required (`winflexbison` on Windows using `choco`) to build the dependency `igraph`.

#### Build

To build from source (compiles C++ extension):

```
git clone https://github.com/GenieTim/pylimer-tools.git
cd pylimer-tools
pip install -e .
```

Optional: Use provided helper scripts in `./bin` (e.g. `./bin/build-wheel.sh`, `./bin/build-tests.sh`).

## 3. Test Installation

```python
import pylimer_tools
import pylimer_tools_cpp

print("Installed version: {} == {}".format(
  pylimer_tools.__version__,
  pylimer_tools_cpp.__version__
))

```

More examples: see the [examples in the documentation](https://genietim.github.io/pylimer-tools/auto_examples/index.html) and [their code in `examples/`](examples/).
Additionally, the [tests](tests/) and the [CLI Tools](src/pylimer_tools/) may serve as examples.

## 4. CLI Tools

Installed console scripts:

- `pylimer-generate-network` – generate random bead-spring networks using our MC generation procedure
- `pylimer-analyse-networks` – batch analysis / statistics given LAMMPS data (structure) files
- `pylimer-basic-lammps-stats` – quick structural stats from LAMMPS data (structure) file
- `pylimer-displace-randomly` – random displacement utility

Run any with `--help` for usage.

## 5. Documentation

Full documentation (API reference, tutorials, examples):
[https://genietim.github.io/pylimer-tools](https://genietim.github.io/pylimer-tools)

## 6. Development & Testing

Clone and install in editable mode (see Installation). Then:

```
./bin/run-tests.sh       # full test suite (Python + C++), includes benchmarks & tests that may fail
./bin/run-tests-short.sh # quicker subset, generates coverage, is what's run in the CI
./bin/format-code.sh     # apply formatting & style (run before PR)
```

Generate docs:

```
./bin/make-stubs.sh      # builds the stubs for the C++ module
./bin/make-docs.sh       # build Sphinx HTML docs
```

### Adding Features

- Add tests in `tests/` (unit or integration). New functionality without tests may be deferred.
- Keep public APIs documented in docstrings so they surface in Sphinx.
- If a change alters behavior, update existing tests rather than deleting them. Explain rationale in the PR.

## 7. Citing

If you use `pylimer-tools` in published work, please cite it. A minimal BibTeX example:

```bibtex
@software{pylimer-tools,
	title = {pylimer-tools},
	author = {Bernhard, Tim},
	url = {https://github.com/GenieTim/pylimer-tools},
	year = {2025},
	note = {See CITATION.cff for full metadata and related method references}
}
```

Also cite the specific theoretical / methodological papers corresponding to the components you use (listed in `CITATION.cff`).

## 8. Contributing

We welcome contributions on GitHub via Issues and Pull Requests.

1. Discuss larger ideas in an Issue first (helps align scope).
2. Fork, branch, implement, add tests & docs.
3. Run formatting and tests locally.
4. Submit PR referencing the Issue (if any).

See also: [Code of Conduct](CODE_OF_CONDUCT.md).

## 9. Code of Conduct

We strive for an inclusive, respectful environment. 
By participating you agree to uphold the [Code of Conduct](CODE_OF_CONDUCT.md). 
Report concerns to the maintainer email specified there.

## 10. Acknowledgements

The authors gratefully acknowledge financial support from the Swiss National Science Foundation (SNSF project 200021_204196).

## 11. License

GPL-3.0-or-later. See [LICENSE](LICENSE).

## Pronunciation Note

"pylimer” resembles "polymer”; the playful spelling emphasizes Python integration.
