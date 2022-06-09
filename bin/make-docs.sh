#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
# ROOT_DIR=$(pwd)

python -m pip install --verbose . || exit 7

# make sure you have sphinx installed:
# pip3 install sphinx
# and the template:
# pip install furo
sphinx-apidoc -o ./docs ./src || exit 2 # -f -P

sphinx-build -b html ./docs ./docs-html

touch ./docs-html/.nojekyll
