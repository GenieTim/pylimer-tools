#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

# Make sure you have the latest version of PyPA’s build installed:
# python3 -m pip install --upgrade build
# python3 -m pip install --upgrade twine

rm -rf dist/

pip install . --verbose
pybind11-stubgen pylimer_tools_cpp -o src
python -m build --sdist
