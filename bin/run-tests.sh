#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 2

# build/install project
python -m pip install --verbose . || exit 3

cd "$(dirname "$0")/../tests" || exit 4

# run tests
python -m coverage run -m unittest discover -v || exit 5

# generate coverage report
python -m coverage report --include="pylimer_tools/**/*.py"
# python -m coverage html --include="pylimer_tools/**/*.py" -d ../coverage.html
