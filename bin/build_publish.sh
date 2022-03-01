#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 2

# Make sure you have the latest version of PyPA’s build installed:
# python3 -m pip install --upgrade build
# python3 -m pip install --upgrade twine

./bin/build-wheel.sh || exit 3

python3 -m twine upload dist/*
