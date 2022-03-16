#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

find ./src \( -name "*.cpp" -o -name "*.h" \) -exec clang-format -i {} \;
