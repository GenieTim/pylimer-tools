#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

find ./src \( -name "*.cpp" -o -name "*.h" \) -exec clang-format -i {} \;
find ./ \(  -name "*.py" -o -name "*.pyi" \)  -exec python -m autopep8 --in-place --aggressive {} \;

git submodule foreach git reset --hard
