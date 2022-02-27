#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

# Make sure you have the latest version of PyPA’s build installed:
# python3 -m pip install --upgrade build
# python3 -m pip install --upgrade twine

rm -rf dist/

python -m build --wheel
# pybind11-stubgen pylimer_tools_cpp -o src
