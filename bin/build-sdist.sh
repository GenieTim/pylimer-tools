#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

# Make sure you have the latest version of PyPA’s build installed:
# python3 -m pip install --upgrade build
# python3 -m pip install --upgrade twine

rm -rf dist/

pip install . --verbose || {
    echo "Failed to install the package. Please check the output for errors."
    exit 1
}

rm -rf dist/

python -m build --sdist 
