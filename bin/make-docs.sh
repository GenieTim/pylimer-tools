#!/usr/bin/env bash

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Utility functions for colored output
print_error() {
    echo -e "${RED}❌ Error: $1${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

cd "$(dirname "$0")/.." || exit 10
# ROOT_DIR=$(pwd)

print_info "Checking and installing package dependencies..."
python -c "import pylimer_tools_cpp" 2>/dev/null ||
  python -m pip install --verbose . || {
  print_error "Failed to install the package. Please check the output above."
  exit 7
}

print_info "Generating thumbnails..."
python ./docs/generate-thumbnails.py || {
  print_error "Failed to generate thumbnails. Please check the output above."
  exit 6
}

# make sure you have sphinx installed:
# pip3 install sphinx
# and the template:
# pip install furo
print_info "Generating Sphinx API documentation..."
sphinx-apidoc -o ./docs ./src || {
  print_error "Failed to generate Sphinx API documentation. Please check the output above."
  exit 8
}

print_info "Building Sphinx documentation..."
sphinx-build -b html ./docs ./docs-html || {
  print_error "Failed to build Sphinx documentation. Please check the output above."
  exit 9
}

touch ./docs-html/.nojekyll

print_success "Documentation generated successfully in ./docs-html"
