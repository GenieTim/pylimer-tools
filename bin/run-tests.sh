#!/usr/bin/env bash

cd "$(dirname "$0")/../src" || exit 2

coverage run -m unittest discover -v || exit 3

coverage report --include="pylimer_tools/**/*.py"
