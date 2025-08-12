#!/usr/bin/env bash

# Get version from latest git tag, with fallback strategies
# This script returns a clean version number suitable for Python/CMake

set -e

# Function to get version from git tag
get_git_version() {
    # Try to get the latest tag
    local latest_tag
    if latest_tag=$(git describe --tags --abbrev=0 2>/dev/null); then
        # Remove 'v' prefix if present
        echo "${latest_tag#v}"
        return 0
    fi
    
    # Fallback: if no tags exist, use 0.0.0-dev
    echo "0.0.0-dev"
    return 1
}

# Function to validate version format
validate_version() {
    local version="$1"
    if [[ ! "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+([.-].*)?$ ]]; then
        echo "Warning: Version '$version' does not match expected pattern" >&2
        return 1
    fi
    return 0
}

# Main execution
if [ $# -eq 0 ]; then
    # No arguments: return clean version
    version=$(get_git_version)
    validate_version "$version" || echo "0.0.0"  # Fallback to 0.0.0
    echo "$version"
elif [ "$1" = "--with-commit" ]; then
    # Include commit info for development builds
    version=$(get_git_version)
    if git_hash=$(git rev-parse --short HEAD 2>/dev/null); then
        echo "${version}+${git_hash}"
    else
        echo "$version"
    fi
elif [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
    echo "Usage: $0 [--with-commit] [--help]"
    echo "  (no args)      Return clean version from git tag"
    echo "  --with-commit  Include git commit hash"
    echo "  --help         Show this help"
else
    echo "Unknown option: $1" >&2
    exit 1
fi
