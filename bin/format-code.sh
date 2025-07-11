#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

# Format C++ code
find ./src/pylimer_tools_cpp \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} \;
find ./tests/pylimer_tools \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} \;

# Format Python code
python -m ruff format ./src/pylimer_tools
python -m ruff format ./examples
python -m ruff format ./tests
find ./src/pylimer_tools \( -name "*.py" -o -name "*.pyi" \) -exec python -m autopep8 --in-place --aggressive {} \;
find ./tests \( -name "*.py" -o -name "*.pyi" \) -exec python -m autopep8 --in-place --aggressive {} \;
find ./examples \( -name "*.py" -o -name "*.pyi" \) -exec python -m autopep8 --in-place --aggressive {} \;
