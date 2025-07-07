#!/usr/bin/env bash

cd "$(dirname "$0")/.." || exit 10
# ROOT_DIR=$(pwd)

python -c "import pylimer_tools_cpp" 2>/dev/null ||
  python -m pip install --verbose . || {
  echo "Failed to install the package. Please check the output above."
  exit 7
}

# make sure you have sphinx installed:
# pip3 install sphinx
# and the template:
# pip install furo
sphinx-apidoc -o ./docs ./src || {
  echo "Failed to generate Sphinx API documentation. Please check the output above."
  exit 8
}

sphinx-build -b html ./docs ./docs-html

touch ./docs-html/.nojekyll
