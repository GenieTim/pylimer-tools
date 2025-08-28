#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
ROOT_DIR=$(pwd)

mkdir -p "$ROOT_DIR/test-reports"

cd "$ROOT_DIR/tests" || exit 2

PATH="/usr/local/opt/llvm/bin:$PATH"

# first, run cpp tests
# rm -rf build; rm -rf vendor/igraph;
mkdir -p build
cd build || exit 5

echo "======== Starting configuration ========"

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
ADDITIONALFLAGS=(-D CMAKE_BUILD_TYPE=Debug)
if [ -n "$WITH_ERROR" ]; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D WITH_ERROR="$WITH_ERROR")
fi
if (command -v clang++ || command -v g++) && ! ([ -z "$CCOMPILER" ] || [ -z "$CXXCOMPILER" ]) && ([ -z "${DIABLE_COVERAGE}" ]); then
  echo "Using compiler: $CCOMPILER and $CXXCOMPILER"
  # Disable LEAK_ANALYSIS on macOS due to AddressSanitizer compatibility issues
  if [[ $OSTYPE == 'darwin'* ]]; then
    ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=ON -D LEAK_ANALYSIS=OFF -D CMAKE_C_COMPILER="$CCOMPILER" -D CMAKE_CXX_COMPILER="$CXXCOMPILER")
    ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CMAKE_CXX_FLAGS="-fno-stack-check" -D CMAKE_C_FLAGS="-fno-stack-check" -D CMAKE_EXE_LINKER_FLAGS="-fno-stack-check")
  else
    ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=ON -D LEAK_ANALYSIS=ON -D CMAKE_C_COMPILER="$CCOMPILER" -D CMAKE_CXX_COMPILER="$CXXCOMPILER")
  fi
else
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=OFF -D LEAK_ANALYSIS=OFF)
fi
if command -v ninja; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" "-GNinja")
  GENERATOR_BIN="ninja"
fi
CXX=$CXXCOMPILER CC=$CCOMPILER cmake .. "${ADDITIONALFLAGS[@]}" || exit 1
