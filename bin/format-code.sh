#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit

PYTHON="${PYTHON:-$(command -v python3 || command -v python)}"

# Format C++ code
if command -v clang-format &>/dev/null; then
  find ./src/pylimer_tools_cpp \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} +
  find ./tests/pylimer_tools_cpp \( -name "*.cpp" -o -name "*.h" \) -exec clang-format --style=file --fallback-style="Mozilla" -i {} +
fi

# Format Python code: ruff (matches .pre-commit-config.yaml) is authoritative;
# autopep8 is only a fallback. Running both makes them fight over line length
# and autopep8's E501 fixer can split f-string replacement fields across lines,
# which is a SyntaxError on Python < 3.12.
if [ -n "$PYTHON" ]; then
  if $PYTHON -m ruff --version &>/dev/null; then
    $PYTHON -m ruff format ./src/pylimer_tools ./examples ./tests
  else
    find ./src/pylimer_tools ./tests ./examples \( -name "*.py" -o -name "*.pyi" \) -exec $PYTHON -m autopep8 --in-place --aggressive --max-line-length=88 --ignore=E402,E501 {} + 2>/dev/null || true
  fi
fi
