#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

rm -rf ./docs
mkdir -p ./docs

# make sure you have pdoc3 installed: 
# pip3 install pdoc3
pdoc ./src/pylimer_tools/ --html --output-dir ./docs-html
