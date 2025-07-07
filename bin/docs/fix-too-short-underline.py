#!/usr/bin/env python
"""
This script fixes too short underlines in reStructuredText (reST) files.
It ensures that underlines are as long as the text they underline.
"""
import os

docs_dir = os.path.join(os.path.dirname(__file__), "../..", "docs")

for root, _, files in os.walk(docs_dir):
    for file in files:
        if file.endswith(".rst"):
            file_path = os.path.join(root, file)
            with open(file_path, "r") as f:
                lines = f.readlines()

            prev_line_length = 0
            with open(file_path, "w") as f:
                for line in lines:
                    if (
                        line.startswith("===")
                        or line.startswith("---")
                        or line.startswith("~~~")
                        or line.startswith("^^^")
                    ):
                        # adjust the underline length if needed
                        if (
                            prev_line_length > 0
                            and len(line.rstrip()) < prev_line_length
                        ):
                            line = (
                                line.rstrip()
                                + line[0] * (prev_line_length - len(line.rstrip()))
                                + "\n"
                            )
                        f.write(line)
                    else:
                        # assume the line is a heading
                        prev_line_length = len(line.rstrip())
                        f.write(line)
            print(f"Processed file: {file_path}")
print("All reST files processed.")
# This script will ensure that all headings in reST files have underlines that match the length
