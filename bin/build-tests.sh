#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
ROOT_DIR=$(pwd)

bash "$ROOT_DIR/bin/config-tests-build.sh"

echo "======== Configuration done ========"

mkdir -p "$ROOT_DIR/test-reports"

cd "$ROOT_DIR/tests/build" || exit 2

echo "======== Starting build ========"
cmake --build . || exit 9
echo "======== Tests built with $GENERATOR_BIN (" "${ADDITIONALFLAGS[@]}" ") ========"
