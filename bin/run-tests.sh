#!/usr/bin/env bash

cd "$(dirname "$0")/../src" || exit

coverage run -m unittest discover -v

coverage report --include="pylimer_tools/**/*.py"
