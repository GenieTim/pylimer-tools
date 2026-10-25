#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

PYTHON="${PYTHON:-$(command -v python3 || command -v python)}"

# Format C++ code
if command -v clang-format &>/dev/null; then
  find ./src/pylimer_tools_cpp \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} +
  find ./tests/pylimer_tools_cpp \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} +
fi

# Format Python code
if [ -n "$PYTHON" ]; then
  $PYTHON -m ruff format ./src/pylimer_tools ./examples ./tests 2>/dev/null || true
  find ./src/pylimer_tools ./tests ./examples \( -name "*.py" -o -name "*.pyi" \) -exec $PYTHON -m autopep8 --in-place --aggressive {} + 2>/dev/null || true
fi
