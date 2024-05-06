#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

find ./src \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} \;
find ./tests \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} \;
find ./src \(  -name "*.py" -o -name "*.pyi" \) -exec python -m autopep8 --in-place --aggressive {} \;
find ./tests \(  -name "*.py" -o -name "*.pyi" \) -exec python -m autopep8 --in-place --aggressive {} \;
