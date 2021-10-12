#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

rm -rf ./docs-html
mkdir -p ./docs-html

# make sure you have pdoc3 installed: 
# pip3 install pdoc3
pdoc ./src/pylimer_tools/ --math --template-directory ./docs-template --output-directory ./docs-html
