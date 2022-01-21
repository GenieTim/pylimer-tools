#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
ROOT_DIR=$(pwd)

cd "$ROOT_DIR/tests" || exit 2

# first, run cpp tests
# rm -rf build
mkdir -p build
cd build || exit 5
cmake .. -DCMAKE_OSX_SYSROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.1.sdk || exit 1
cmake --build . || exit 9
ASAN_OPTIONS=detect_leaks=1 ./pylimer_tests || exit 6 # -s --durations yes 

cd "$ROOT_DIR" || exit 8

# then, build/install project for Python
python -m pip install --verbose --use-feature=in-tree-build . || exit 3

cd "$ROOT_DIR" || exit 4

# then, run Python tests
python -m coverage run -m unittest discover -v || exit 7

# generate coverage report
python -m coverage report --include="src/**/*.py"
# python -m coverage html --include="pylimer_tools/**/*.py" -d ../coverage.html
