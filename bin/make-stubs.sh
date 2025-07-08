#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

python -m pip install . --verbose || {
    echo "Failed to install the package. Please check the output for errors."
    exit 1
}
pybind11-stubgen pylimer_tools_cpp -o src/pylimer_tools_cpp-stubs --numpy-array-remove-parameters
mv src/pylimer_tools_cpp-stubs/pylimer_tools_cpp.pyi src/pylimer_tools_cpp-stubs/__init__.pyi
python -m pip install .
