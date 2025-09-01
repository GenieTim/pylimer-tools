#!/usr/bin/env python3
"""
Version management script for pylimer-tools.

This script provides a single source of truth for version information by:
1. Attempting to extract version from git tags (when available)
2. Falling back to a version file (_version.txt) 
3. Writing/updating the version file when git version is available and differs

The version file is created during build time to ensure it's available in sdist/tarballs
where git information is not present.
"""

import argparse
import os
import re
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Optional


def get_project_root() -> Path:
    """Get the project root directory."""
    # Start from this script's directory and walk up to find project root
    current = Path(__file__).parent
    while current != current.parent:
        if (current / "setup.py").exists() or (current / "pyproject.toml").exists():
            return current
        current = current.parent
    # Fallback to parent of bin directory
    return Path(__file__).parent.parent


def normalize_version(version: str) -> str:
    """Normalize version string by removing 'v' prefix and validating format."""
    if not version:
        return "0.0.0"
    
    # Remove 'v' prefix if present
    if version.startswith('v'):
        version = version[1:]
    
    # Validate basic version format (X.Y.Z with optional suffix)
    if not re.match(r'^\d+\.\d+\.\d+', version):
        return "0.0.0"
    
    return version.strip()


def get_version_from_git() -> Optional[str]:
    """
    Attempt to get version from git tags.
    
    Returns:
        Version string if successful, None if git is unavailable or no tags found
    """
    try:
        # First, check if we're in a git repository
        result = subprocess.run(
            ["git", "rev-parse", "--git-dir"],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=get_project_root()
        )
        if result.returncode != 0:
            return None
        
        # Try to get the latest tag
        result = subprocess.run(
            ["git", "describe", "--tags", "--abbrev=0"],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=get_project_root()
        )
        
        if result.returncode == 0 and result.stdout.strip():
            version = normalize_version(result.stdout.strip())
            if version != "0.0.0":
                return version
        
        return None
        
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return None


def get_version_file_path() -> Path:
    """Get the path to the version file."""
    return get_project_root() / "_version.txt"


def read_version_from_file() -> Optional[str]:
    """
    Read version from the version file.
    
    Returns:
        Version string if file exists and is readable, None otherwise
    """
    version_file = get_version_file_path()
    try:
        if version_file.exists():
            version = version_file.read_text().strip()
            if version:
                return normalize_version(version)
    except (OSError, IOError):
        pass
    return None


def write_version_to_file(version: str) -> bool:
    """
    Write version to the version file.
    
    Args:
        version: Version string to write
        
    Returns:
        True if successful, False otherwise
    """
    try:
        version_file = get_version_file_path()
        version_file.write_text(version + "\n")
        return True
    except (OSError, IOError):
        return False


def get_version(update_file: bool = True) -> str:
    """
    Get the current version using the following priority:
    1. Git tags (if available)
    2. Version file
    3. Fallback to "0.0.0"
    
    Args:
        update_file: Whether to update the version file if git version is available
        
    Returns:
        Current version string
    """
    # Try to get version from git first
    git_version = get_version_from_git()
    
    if git_version:
        # Git version is available
        if update_file:
            # Check if we need to update the version file
            file_version = read_version_from_file()
            if file_version != git_version:
                if not write_version_to_file(git_version):
                    warnings.warn("Could not update version file", file=sys.stderr)
        return git_version
    
    # Fallback to version file
    file_version = read_version_from_file()
    if file_version:
        return file_version
    
    # Final fallback
    return "0.0.0"


def ensure_version_file() -> str:
    """
    Ensure the version file exists with the current version.
    This is useful for packaging where we want to guarantee the file exists.
    
    Returns:
        Current version string
    """
    version = get_version(update_file=True)
    
    # Ensure the file exists (even if it already has the right content)
    version_file = get_version_file_path()
    if not version_file.exists():
        write_version_to_file(version)
        # print(f"Created version file with version {version}", file=sys.stderr)
    
    return version


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Get or manage version information for pylimer-tools"
    )
    parser.add_argument(
        "--no-update", 
        action="store_true",
        help="Don't update the version file even if git version is available"
    )
    parser.add_argument(
        "--ensure-file",
        action="store_true", 
        help="Ensure version file exists (useful for packaging)"
    )
    parser.add_argument(
        "--file-only",
        action="store_true",
        help="Only read from version file, ignore git"
    )
    parser.add_argument(
        "--git-only", 
        action="store_true",
        help="Only read from git, ignore version file"
    )
    
    args = parser.parse_args()
    
    if args.file_only:
        version = read_version_from_file()
        if version is None:
            print("No version file found", file=sys.stderr)
            sys.exit(1)
        print(version)
    elif args.git_only:
        version = get_version_from_git()
        if version is None:
            print("No git version available", file=sys.stderr)
            sys.exit(1)
        print(version)
    elif args.ensure_file:
        version = ensure_version_file()
        print(version)
    else:
        version = get_version(update_file=not args.no_update)
        print(version)


if __name__ == "__main__":
    main()
