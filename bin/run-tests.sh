#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 2

python -m pip install --verbose . || exit 3

python -m coverage run -m unittest discover -v || exit 4

cd "$(dirname "$0")/../src" || exit 5

# python -m coverage report --include="pylimer_tools/**/*.py"
python -m coverage html --include="pylimer_tools/**/*.py" -d ../coverage.html
