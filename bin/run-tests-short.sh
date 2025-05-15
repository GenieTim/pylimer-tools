#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
ROOT_DIR=$(pwd)

# Create arrays to store timing information
# (associative arrays are not supported in bash < 4, 
# and MacOS comes with 3, so let's go without)
timing_names=()
timing_values=()

# Function to run a command with timing and store the result
run_timed() {
  local name=$1
  shift
  local start_time
  start_time=$(date +%s.%N)
  time "$@" || return $?
  local end_time
  end_time=$(date +%s.%N)
  local elapsed
  elapsed=$(echo "$end_time - $start_time" | bc)
  
  # Store timing information in parallel arrays
  timing_names+=("$name")
  timing_values+=("$elapsed")

  return 0
}

# pip/skbuild uses ninja as a generator,
# however, it uses a bundled one in a virtual env
# therefore, we need to delete vendor caches
if [ -d "_skbuild" ]; then
  rm -rf ./_skbuild
fi

# first, run cpp tests
# (build them first)
run_timed "Build Tests" "$ROOT_DIR/bin/build-tests.sh" || exit 3
cd "$ROOT_DIR/tests/build" || exit 2

GENERATOR_BIN="make"
if command -v ninja; then
  GENERATOR_BIN="ninja"
fi

ls ./*
echo "======== Starting gcov ========"
find . -name "*Universe.cpp.gcov" -exec cat {} \;
MallocNanoZone=0 ASAN_OPTIONS=detect_leaks=1:detect_container_overflow=0:strict_string_checks=1:detect_stack_use_after_return=1:check_initialization_order=1:strict_init_order=1 LSAN_OPTIONS=suppressions=$ROOT_DIR/tests/lsan.supp run_timed "Run Tests (gcov)" "$GENERATOR_BIN" pylimer_tests-gcov
find . -name "*Universe.cpp.gcov" -exec cat {} \;
echo "========== /ran gcov =========="
echo "======= Analyzing output ======"
ls
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

if command -v npx; then
  run_timed "C++ Coverage Badge" npx -y lcov-badge2 -l "C++ Code Coverage" -o ".github/cpp-coverage.svg" tests/build/lcov/data/capture/pylimer_tools.info || echo "Failed to generate coverage badge"
fi

# then, build/install project for Python
echo "====== Installing with PIP ======"
run_timed "PIP Installation" python ./bin/dev_install.py || exit 3
echo "===== /Installing with PIP ======"

cd "$ROOT_DIR" || exit 4

# then, run Python tests with XML output
run_timed "Python Tests" python -m coverage run -m xmlrunner discover -v -o ./test-reports || exit 7

# generate coverage report
python -m coverage report
python -m coverage xml

# Print performance summary
echo ""
echo "====== PERFORMANCE SUMMARY ======"
for i in "${!timing_names[@]}"; do
  printf "%-22s: %7.2f s\n" "${timing_names[$i]}" "${timing_values[$i]}"
done
echo "================================="

# Calculate total runtime
total_time=0
for i in "${!timing_values[@]}"; do
  total_time=$(echo "$total_time + ${timing_values[$i]}" | bc)
done
printf "%-22s: %7.2f s\n" "Total Runtime" "$total_time"
echo "================================="


exit 0
