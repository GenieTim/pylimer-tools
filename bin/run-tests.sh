#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
ROOT_DIR=$(pwd)


# pip/skbuild uses ninja as a generator, 
# however, it uses a bundled one in a virtual env
# therefore, we need to delete vendor caches
if [ -d "_skbuild" ]; then
  rm -rf ./_skbuild
  rm -rf ./vendor/igraph/src/igraphLib-build
fi

cd "$ROOT_DIR/tests" || exit 2

# first, run cpp tests
# rm -rf build
mkdir -p build
cd build || exit 5
# force use of g++ if available for coverage
CXXCOMPILER=$(which g++ || which clang)
CCOMPILER=$(which gcc || which clang)
ADDITIONALFLAGS=()
if command -v g++; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=ON -D LEAK_ANALYSIS=OFF -D CMAKE_C_COMPILER="$CCOMPILER" -D CMAKE_CXX_COMPILER="$CXXCOMPILER")
elif command -v clang; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=OFF -D LEAK_ANALYSIS=ON -D CMAKE_C_COMPILER="$CCOMPILER" -D CMAKE_CXX_COMPILER="$CXXCOMPILER")
else
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=OFF -D LEAK_ANALYSIS=OFF)
fi
if command -v ninja; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" "-GNinja")
fi
cmake .. -D CODE_COVERAGE=ON -D LEAK_ANALYSIS=OFF "${ADDITIONALFLAGS[@]}" || exit 1
cmake --build . || exit 9
echo "======== Starting tests ========"
ASAN_OPTIONS=detect_leaks=1 ./pylimer_tests || exit 6 # -s --durations yes
make pylimer_tests-gcov
make test_sources-gcov
make pylimer_tools-gcov
make header_tests-gcov

# TODO: the following is
# some sort of fix for header file coverage being assembled incorrectly
# but this fix does not really work
rm ./*"#src#pylimer_tools_cpp#calc#^#entities#Atom.h.gcov"
rm pylimer_tools.out/*"#entities#Atom.h.gcov"
rm ./*"#pylimer_tools_cpp#entities#Atom.h.gcov"
# ENDTODO

make pylimer_tests-geninfo
make pylimer_tools-geninfo
make header_tests-geninfo
make pylimer_tools-genhtml

cd "$ROOT_DIR" || exit 8

# copy outside such that pip installation does not remove it
# cp tests/build/lcov/data/capture/pylimer_tools.info pylimer_tools_lcoverage.info

if command -v npx; then
  npx -y lcov-badge2 -l "C++ Code Coverage" -o ".github/cpp-coverage.svg" pylimer_tools_lcoverage.info || echo "Failed to generate coverage badge"
fi

# then, build/install project for Python
python -m pip install --verbose --use-feature=in-tree-build . || exit 3

cd "$ROOT_DIR" || exit 4

# then, run Python tests
python -m coverage run -m unittest discover -v || exit 7

# generate coverage report
python -m coverage report --include="src/**/*.py"
python -m coverage xml --include="src/**/*.py"
# python -m coverage html --include="pylimer_tools/**/*.py" -d ../coverage.html

exit 0
