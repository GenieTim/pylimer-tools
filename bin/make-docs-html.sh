#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

rm -rf ./docs-html
mkdir -p ./docs-html

# make sure you have pdoc3 installed: 
# pip3 install pdoc3
pdoc ./src/pylimer_tools/ --template-dir ./docs-template --config latex_math=True --html --output-dir ./docs-html
