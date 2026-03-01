#!/usr/bin/env python

import os
import subprocess

"""
This is a custom wrapper around the test binary,
since that one does not necessarily exit
with the right error code when a test fails.
"""

test_bin = os.path.join(os.path.dirname(__file__), "build", "pylimer_tests")

os.chdir(os.path.dirname(test_bin))
all_tests = subprocess.run(
    [test_bin, "--list-tests"],
    capture_output=True,  # Python >= 3.7 only
    text=True,  # Python >= 3.7 only
)

tests_output = {}

n_tests = 0
n_tests_failed = 0

all_lines = all_tests.stdout.split("\n")
for i, line in enumerate(all_lines):
    if i == 0:
        continue
    line = line.strip()
    if line == "":
        continue

    if line.startswith("["):
        continue

    n_tests += 1
    test_out = subprocess.run(
        [test_bin, line],
        capture_output=True,  # Python >= 3.7 only
        text=True,  # Python >= 3.7 only
    )
    if test_out.returncode != 0:
        n_tests_failed += 1
    print("Test {}: {}".format(line, "PASS" if test_out.returncode == 0 else "FAIL"))
    tests_output[line] = test_out

print("Ran {} tests, {} failed.".format(n_tests, n_tests_failed))

for test, result in tests_output.items():
    if result.returncode != 0:
        print(
            "Test {} failed with output:\n{}\n{}\n\n".format(
                test, result.stderr, result.stdout
            )
        )

if n_tests_failed > 0:
    exit(1)
else:
    exit(0)
