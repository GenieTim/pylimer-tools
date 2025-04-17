#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

find ./src/pylimer_tools_cpp \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} \;
find ./tests/pylimer_tools \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} \;
find ./src/pylimer_tools \(  -name "*.py" -o -name "*.pyi" \) -exec python -m autopep8 --in-place --aggressive {} \;
find ./tests \(  -name "*.py" -o -name "*.pyi" \) -exec python -m autopep8 --in-place --aggressive {} \;

./bin/build-tests.sh

find ./src/pylimer_tools_cpp \( -name "*.cpp" -o -name "*.h" \) -exec clang-tidy -p ./tests/build {} --fix --header-filter=".*" \;
find ./tests/pylimer_tools \( -name "*.cpp" -o -name "*.h" \) -exec clang-tidy -p ./tests/build {} --fix --header-filter=".*" \;
