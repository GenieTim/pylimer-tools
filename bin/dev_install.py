#!/usr/bin/env python
"""
Development installation script that preserves the build directory
and enables faster incremental builds.
"""
import os
import subprocess
import sys

if __name__ == "__main__":
    # Run pip install with development flags
    cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "--editable",
        ".",
        "--no-deps",
        "--no-build-isolation",
    ]

    # Add any additional arguments
    if len(sys.argv) > 1:
        cmd.extend(sys.argv[1:])

    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=os.path.join(os.path.dirname(__file__), ".."))
