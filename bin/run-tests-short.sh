#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
ROOT_DIR=$(pwd)

# pip/skbuild uses ninja as a generator,
# however, it uses a bundled one in a virtual env
# therefore, we need to delete vendor caches
if [ -d "_skbuild" ]; then
  rm -rf ./_skbuild
fi

# first, run cpp tests
# (build them first)
"$ROOT_DIR/bin/build-tests.sh" || exit 3
cd "$ROOT_DIR/tests/build" || exit 2

GENERATOR_BIN="make"
if command -v ninja; then
  GENERATOR_BIN="ninja"
fi

ls ./*
echo "======== Starting gcov ========"
find . -name "*Universe.cpp.gcov" -exec cat {} \;
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1:detect_container_overflow=0:strict_string_checks=1:detect_stack_use_after_return=1:check_initialization_order=1:strict_init_order=1 LSAN_OPTIONS=suppressions=$ROOT_DIR/tests/lsan.supp time "$GENERATOR_BIN" pylimer_tests-gcov
find . -name "*Universe.cpp.gcov" -exec cat {} \;
echo "========== /ran gcov =========="
echo "======= Analyzing output ======"
for of in ./*_test_output.log.txt; do
  echo "==== $of ===="
  cat "$of"
  echo "==== /$of ===="
  if grep -q "FAILED:" "$of"; then
    echo "!!!!! TESTS FAILED !!!!!"
    exit 6
  fi
done
echo "====== /Analyzing output ======"
ls ./*

cd "$ROOT_DIR" || exit 8

# copy outside such that pip installation does not remove it
# cp tests/build/lcov/data/capture/pylimer_tools.info pylimer_tools.info

if command -v npx; then
  npx -y lcov-badge2 -l "C++ Code Coverage" -o ".github/cpp-coverage.svg" tests/build/lcov/data/capture/pylimer_tools.info || echo "Failed to generate coverage badge"
fi

# then, build/install project for Python
python -m pip install --verbose -e . || exit 3

cd "$ROOT_DIR" || exit 4

# then, run Python tests with XML output
python -m coverage run -m xmlrunner discover -v -o ./test-reports || exit 7

# generate coverage report
python -m coverage report
python -m coverage xml
# python -m coverage html --include="pylimer_tools/**/*.py" -d ../coverage.html

exit 0
