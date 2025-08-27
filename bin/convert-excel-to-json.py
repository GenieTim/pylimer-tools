#!/usr/bin/env python3
"""
Convert Everaers et al. Excel data to JSON format.

This script converts the polymer properties Excel file to JSON format,
eliminating the need for pandas and openpyxl as runtime dependencies.

Usage:
    python bin/convert-excel-to-json.py

This script should be run during development when the Excel file is updated.
"""

import json
import os
import sys
from pathlib import Path
import pandas as pd


def convert_excel_to_json():
    """Convert the Everaers et al. Excel file to JSON format."""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    excel_path = project_root / "src" / "pylimer_tools" / "data" / "everaers_et_al_unit_properties.xlsx"
    json_path = project_root / "src" / "pylimer_tools" / "data" / "everaers_et_al_unit_properties.json"
    
    if not excel_path.exists():
        print(f"Error: Excel file not found at {excel_path}")
        sys.exit(1)
    
    print(f"Converting {excel_path} to {json_path}")
    
    try:
        # Read the Excel file
        df = pd.read_excel(excel_path)
        
        # Convert to JSON format that preserves the structure
        data = {
            'polymers': df.to_dict('records'),
            'columns': list(df.columns),
            'metadata': {
                'source': 'Everaers et al. (2020)',
                'description': 'Polymer unit properties for LJ unit conversions',
                'generated_from': excel_path.name,
                'conversion_note': 'This file was generated from the Excel source. Do not edit manually.',
                'reference': 'https://doi.org/10.1021/acs.macromol.9b02428'
            }
        }
        
        # Save as JSON with nice formatting
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2, sort_keys=False)
        
        print(f"✓ Successfully converted Excel data to {json_path}")
        print(f"  Data contains {len(data['polymers'])} polymers with {len(data['columns'])} properties each")
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)


if __name__ == "__main__":
    convert_excel_to_json()
