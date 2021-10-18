#!/usr/bin/env bash

cd "$(dirname "$0")/../src" || exit 2

python -m coverage run -m unittest discover -v || exit 3

# python -m coverage report --include="pylimer_tools/**/*.py"
python -m coverage html --include="pylimer_tools/**/*.py" -d ../coverage.html
