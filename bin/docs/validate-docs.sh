#!/usr/bin/env bash

cd "$(dirname "$0")/../.." || exit 10
# ROOT_DIR=$(pwd)

# iterate over all *.rst files in the docs directory
for file in docs/**/*.rst; do
  # check if the file exists
  if [[ -f "$file" ]]; then
    # validate the file using rstcheck
    rstcheck "$file" --ignore-directives automodule,autosummary,autoclass,autofunction,autocode || {
      echo "Validation failed for $file. Please check the output above."
      echo "--------------------------"
      echo ""
    }
  else
    echo "File $file does not exist."
    exit 2
  fi
done
