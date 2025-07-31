# pylimer-tools

[![Test Coverage of Python Code](https://github.com/GenieTim/pylimer-tools/blob/main/.github/coverage.svg?raw=true)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)
[![Test Coverage of C++ Code](https://github.com/GenieTim/pylimer-tools/blob/main/.github/cpp-coverage.svg?raw=true)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)
[![Total Code Test Coverage](https://codecov.io/gh/GenieTim/pylimer-tools/branch/main/graph/badge.svg?token=5ZE1VSDXJQ)](https://codecov.io/gh/GenieTim/pylimer-tools)
[![Run Tests](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml/badge.svg)](https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml)
[![Publish Documentation](https://github.com/GenieTim/pylimer-tools/actions/workflows/publish-documentation-html.yml/badge.svg)](https://github.com/GenieTim/pylimer-tools/actions/workflows/publish-documentation-html.yml)[![PyPI version](https://badge.fury.io/py/pylimer-tools.svg)](https://badge.fury.io/py/pylimer-tools)
[![PyPI download month](https://img.shields.io/pypi/dm/pylimer-tools.svg)](https://pypi.python.org/pypi/pylimer-tools/)
[![PyPI license](https://img.shields.io/pypi/l/pylimer-tools.svg)](https://pypi.python.org/pypi/pylimer-tools/)

pylimer-tools (with "pylimer" pronounced like "pü-limer", with the "py" as in the German word "müde", IPA: /ˈpyːlɪmɚ/) is a collection of utility Python functions for handling handling bead-spring polymers in Python.

This toolbox provides an Monte Carlo (MC) structure generator, 
simulation methods for Dissipative Particle Dynamics (DPD) with slip-springs,
implementations of the Maximum Entropy Homogenization Procedure (MEHP), with and without slip-links,
and various means to read LAMMPS output: be it data, dump or thermodynamic data files.
Additionally, it provides various methods to calculate properties of the bead-spring networks, such as computing the radius of gyration, mean end to end distance, finding loops, or simply splitting a polymer network back up into its chains.

## Installation

Use pip:

`pip install pylimer-tools`

or refer to the [documentation](https://genietim.github.io/pylimer-tools) for instructions on how to compile the package yourself.

## Usage

See the [documentation](https://genietim.github.io/pylimer-tools) for a current list of all available functions.

Example usage can be found in the [documentation](https://genietim.github.io/pylimer-tools/auto_examples/index.html), the [tests](https://github.com/GenieTim/pylimer-tools/tree/main/tests), the [examples](https://github.com/GenieTim/pylimer-tools/tree/main/examples)
or the [CLI applications](https://github.com/GenieTim/pylimer-tools/tree/main//src/pylimer_tools/).

## Contributing

We welcome any contributions right here on [GitHub](https://github.com/GenieTim/pylimer-tools/).
If you have questions, open an [Issue](https://github.com/GenieTim/pylimer-tools/issues), and if you can contribute with code, open a [Pull Request](https://github.com/GenieTim/pylimer-tools/pulls).

If you add new methods or functionality, we kindly ask you to add appropriate [tests](./tests/).

If you change existing functionality, we require the existing tests to pass, or be modified accordingly if they turn out to have been wrong.

Be aware of the scripts provided in [`./bin`](./bin/) to improve your development experience; for example, run `./bin/format-code.sh` before submitting a PR.

## Acknowledgements

The authors gratefully acknowledge financial support from the
Swiss National Science Foundation (SNSF project 200021_204196).
