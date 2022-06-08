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
# rm -rf build; rm -rf vendor/igraph;
mkdir -p build
cd build || exit 5
rm ./**/*.gcda
GENERATOR_BIN="make"
# force use of g++ if available for coverage
# or clang on MacOS, as g++ leak analysis is not supported there
if [[ $OSTYPE == 'darwin'* ]]; then
  CXXCOMPILER=$(which clang++ || which g++)
  CCOMPILER=$(which clang || which gcc)
else
  CXXCOMPILER=$(which g++ || which clang++)
  CCOMPILER=$(which gcc || which clang)
fi
ADDITIONALFLAGS=()
if command -v clang++ || command -v g++; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=ON -D LEAK_ANALYSIS=ON -D CMAKE_C_COMPILER="$CCOMPILER" -D CMAKE_CXX_COMPILER="$CXXCOMPILER")
else
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=OFF -D LEAK_ANALYSIS=OFF)
fi
if command -v ninja; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" "-GNinja")
  GENERATOR_BIN="ninja"
fi
cmake .. "${ADDITIONALFLAGS[@]}" || exit 1
cmake --build . || exit 9
echo "======== Starting tests ========"
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1:detect_container_overflow=0:strict_string_checks=1:detect_stack_use_after_return=1:check_initialization_order=1:strict_init_order=1 LSAN_OPTIONS=suppressions=../lsan.supp ./pylimer_tests || exit 6 # -s --durations yes
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1:detect_container_overflow=0:strict_string_checks=1:detect_stack_use_after_return=1:check_initialization_order=1:strict_init_order=1 LSAN_OPTIONS=suppressions=../lsan.supp ./header_tests || exit 7
"$GENERATOR_BIN" pylimer_tests-gcov
"$GENERATOR_BIN" test_sources-gcov
"$GENERATOR_BIN" pylimer_tools-gcov
"$GENERATOR_BIN" header_tests-gcov

"$GENERATOR_BIN" pylimer_tests-geninfo
"$GENERATOR_BIN" pylimer_tools-geninfo
"$GENERATOR_BIN" header_tests-geninfo
"$GENERATOR_BIN" pylimer_tools-genhtml

cd "$ROOT_DIR" || exit 8

# copy outside such that pip installation does not remove it
# cp tests/build/lcov/data/capture/pylimer_tools.info pylimer_tools.info

if command -v npx; then
  npx -y lcov-badge2 -l "C++ Code Coverage" -o ".github/cpp-coverage.svg" tests/build/lcov/data/capture/pylimer_tools.info || echo "Failed to generate coverage badge"
fi

# then, build/install project for Python
python -m pip install --verbose . || exit 3

cd "$ROOT_DIR" || exit 4

# then, run Python tests
python -m coverage run -m unittest discover -v || exit 7

# generate coverage report
python -m coverage report --include="src/**/*.py"
python -m coverage xml --include="src/**/*.py"
# python -m coverage html --include="pylimer_tools/**/*.py" -d ../coverage.html

exit 0
