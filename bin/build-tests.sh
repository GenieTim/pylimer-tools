#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
ROOT_DIR=$(pwd)

cd "$ROOT_DIR/tests" || exit 2

# first, run cpp tests
# rm -rf build; rm -rf vendor/igraph;
mkdir -p build
cd build || exit 5
# rm ./**/*.gcda
GENERATOR_BIN="make"
# force use of g++ if available for coverage
# or clang on MacOS, as g++ leak analysis is not supported there
if [ -z "$CC" ] || [ -z "$CXX" ]; then
  if [[ $OSTYPE == 'darwin'* ]]; then
    CXXCOMPILER=$(which clang++ || which g++)
    CCOMPILER=$(which clang || which gcc)
  else
    CXXCOMPILER=$(which g++ || which clang++)
    CCOMPILER=$(which gcc || which clang)
  fi
else
  CXXCOMPILER=$CXX
  CCOMPILER=$CC
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
echo "======== Starting build ========"
cmake --build . || exit 9
