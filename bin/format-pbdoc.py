#!/usr/bin/env python

import os
import re
import sys
from pathlib import Path
import argparse
import textwrap


def find_pbdoc_strings(content):
    """Find all pbdoc raw strings in the content."""
    pattern = r'R"pbdoc\((.*?)\)pbdoc"'
    return re.finditer(pattern, content, re.DOTALL)


def format_pbdoc_string(pbdoc_string):
    """Format a pbdoc string with consistent indentation."""
    # Split the string into lines
    lines = pbdoc_string.split('\n')
    
    # Remove empty lines at the beginning and end
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    
    if not lines:
        return ""
    
    # Determine the minimum indentation (ignoring empty lines)
    non_empty_lines = [line for line in lines if line.strip()]
    if not non_empty_lines:
        return ""
    
    # Find leading whitespace in each non-empty line
    leading_spaces = [len(line) - len(line.lstrip()) for line in non_empty_lines]
    if leading_spaces:
        min_indent = min(leading_spaces)
    else:
        min_indent = 0
    
    # Special handling for notes, warnings, etc.
    special_blocks = []
    current_block = []
    in_special_block = False
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        # Check for special block markers
        if stripped.startswith(('.. note::', '.. warning::', '.. seealso::', '.. admonition::')):
            if current_block:
                special_blocks.append((current_block[0], len(current_block)))
            current_block = [i]
            in_special_block = True
        elif in_special_block:
            current_block.append(i)
            # Check if this line might end the special block
            if not stripped or i == len(lines) - 1:
                special_blocks.append((current_block[0], len(current_block)))
                current_block = []
                in_special_block = False
    
    # If we're still in a special block at the end
    if current_block:
        special_blocks.append((current_block[0], len(current_block)))
    
    # Create a set of line indices that are part of special blocks
    special_line_indices = set()
    for start, length in special_blocks:
        for i in range(start, start + length):
            special_line_indices.add(i)
    
    # Format each line
    formatted_lines = []
    for i, line in enumerate(lines):
        if i in special_line_indices:
            # For special blocks, preserve the original indentation
            formatted_lines.append(line)
        else:
            # For regular lines, normalize the indentation
            if line.strip():  # Non-empty line
                formatted_lines.append(' ' * 4 + line[min_indent:])
            else:  # Empty line
                formatted_lines.append('')
    
    return '\n'.join(formatted_lines)


def process_file(file_path, dry_run=False):
    """Process a single file to format pbdoc strings."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find all pbdoc strings
    matches = list(find_pbdoc_strings(content))
    if not matches:
        print(f"No pbdoc strings found in {file_path}")
        return False
    
    # Process the file content from end to beginning to avoid offset issues
    new_content = content
    for match in reversed(matches):
        original_string = match.group(1)
        formatted_string = format_pbdoc_string(original_string)
        
        # Replace the original string with the formatted one
        start, end = match.span(1)
        new_content = new_content[:start] + formatted_string + new_content[end:]
    
    # Check if content has changed
    if new_content == content:
        print(f"No changes needed for {file_path}")
        return False
    
    # Write the updated content back to the file
    if not dry_run:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(new_content)
        print(f"Updated {file_path}")
    else:
        print(f"Would update {file_path} (dry run)")
    
    return True


def main():
    parser = argparse.ArgumentParser(description='Format pbdoc strings in C++ files.')
    parser.add_argument('--path', default='src/pylimer_tools_cpp/pybind11', 
                        help='Path to the directory containing C++ files')
    parser.add_argument('--dry-run', action='store_true', 
                        help='Show what would be done without making changes')
    args = parser.parse_args()
    
    # Find the project root (assuming this script is in the bin directory)
    script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
    project_root = script_dir.parent
    
    # Get the full path to the target directory
    target_dir = project_root / args.path
    
    if not target_dir.exists() or not target_dir.is_dir():
        print(f"Error: Directory {target_dir} does not exist")
        return 1
    
    # Process all .cpp files in the directory
    files_processed = 0
    files_updated = 0
    
    for file_path in target_dir.glob('*.cpp'):
        files_processed += 1
        if process_file(file_path, args.dry_run):
            files_updated += 1
    
    print(f"\nProcessed {files_processed} files, updated {files_updated} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
