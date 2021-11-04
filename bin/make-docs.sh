#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
# ROOT_DIR=$(pwd)

pip install . || exit 7

# make sure you have sphinx installed:
# pip3 install sphinx
sphinx-apidoc -f -o ./docs ./src || exit 2

sphinx-build -b html ./docs ./docs-html
