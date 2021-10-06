#!/usr/bin/env bash

cd "$(dirname "$0")/../src" || exit

python -m unittest discover
