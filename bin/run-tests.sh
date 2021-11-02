#!/usr/bin/env bash

cd "$(dirname "$0")/../tests" || exit 2

# first, run cpp tests
mkdir -p build
cd build || exit 5
cmake ..
cmake --build .
./pylimer_tests || exit 6

cd "$(dirname "$0")/.." || exit 8

# then, build/install project for Python
python -m pip install --verbose . || exit 3

cd "$(dirname "$0")/../tests" || exit 4

# then, run Python tests
python -m coverage run -m unittest discover -v || exit 7

# generate coverage report
python -m coverage report --include="pylimer_tools/**/*.py"
# python -m coverage html --include="pylimer_tools/**/*.py" -d ../coverage.html
