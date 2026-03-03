#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
ROOT_DIR=$(pwd)

mkdir -p "$ROOT_DIR/test-reports"

cd "$ROOT_DIR/tests" || exit 2

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
    CXXCOMPILER=$(xcrun --find clang++)
    CCOMPILER=$(xcrun --find clang)
  else
    CXXCOMPILER=$(which g++ || which clang++)
    CCOMPILER=$(which gcc || which clang)
  fi
else
  CXXCOMPILER=$CXX
  CCOMPILER=$CC
fi
BUILD_TYPE=${CMAKE_BUILD_TYPE:-Debug}
ADDITIONALFLAGS=(-D CMAKE_BUILD_TYPE="$BUILD_TYPE")
if [ -n "$WITH_ERROR" ]; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D WITH_ERROR="$WITH_ERROR")
fi
LEAK_ANALYSIS_FLAG="OFF"
if [ "${ENABLE_LEAK_ANALYSIS:-1}" = "1" ]; then
  LEAK_ANALYSIS_FLAG="ON"
fi
if (command -v clang++ || command -v g++) && ! ([ -z "$CCOMPILER" ] || [ -z "$CXXCOMPILER" ]) && ([ -z "${DISABLE_COVERAGE}" ]); then
  echo "Using compiler: $CCOMPILER and $CXXCOMPILER"
  # Disable LEAK_ANALYSIS on macOS due to AddressSanitizer compatibility issues
  if [[ $OSTYPE == 'darwin'* ]]; then
    ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=ON -D LEAK_ANALYSIS=OFF -D CMAKE_C_COMPILER="$CCOMPILER" -D CMAKE_CXX_COMPILER="$CXXCOMPILER")
    ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CMAKE_CXX_FLAGS="-fno-stack-check" -D CMAKE_C_FLAGS="-fno-stack-check" -D CMAKE_EXE_LINKER_FLAGS="-fno-stack-check")
  else
    ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=ON -D LEAK_ANALYSIS="$LEAK_ANALYSIS_FLAG" -D CMAKE_C_COMPILER="$CCOMPILER" -D CMAKE_CXX_COMPILER="$CXXCOMPILER")
  fi
else
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CODE_COVERAGE=OFF -D LEAK_ANALYSIS=OFF)
fi
if command -v ccache >/dev/null 2>&1; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache)
fi
if command -v ninja; then
  ADDITIONALFLAGS=("${ADDITIONALFLAGS[@]}" "-GNinja")
  GENERATOR_BIN="ninja"
fi
CXX=$CXXCOMPILER CC=$CCOMPILER cmake .. "${ADDITIONALFLAGS[@]}" || exit 1
