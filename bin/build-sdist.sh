#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

# Make sure you have the latest version of PyPA’s build installed:
# python3 -m pip install --upgrade build
# python3 -m pip install --upgrade twine

rm -rf dist/

pip install . --verbose
pybind11-stubgen pylimer_tools_cpp -o src/pylimer_tools_cpp-stubs --numpy-array-remove-parameters
mv src/pylimer_tools_cpp-stubs/pylimer_tools_cpp.pyi src/pylimer_tools_cpp-stubs/__init__.pyi
rm -rf dist/
python -m build --sdist
