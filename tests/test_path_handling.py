#!/usr/bin/env python3
"""
Test script to verify path handling fixes.
"""

import os
import sys
from pathlib import Path

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# Import the path conversion function
from analysis_pipeline.core.utils import convert_win_path_to_docker
from analysis_pipeline.docker.docker_utils import convert_win_path_to_docker as docker_convert_win_path_to_docker

def test_path_conversion():
    """Test the path conversion function."""
    test_paths = [
        r"C:\Users\spejo\Documents\Scripts\Bindra_Lab_CRISPR_analysis",
        r"C:\Program Files\Docker\Docker",
        r"D:\Data\CRISPR",
        r"\\server\share\data",
    ]
    
    print("Testing path conversion in utils.py:")
    for path in test_paths:
        try:
            converted = convert_win_path_to_docker(path)
            print(f"  Original: {path}")
            print(f"  Converted: {converted}")
        except Exception as e:
            print(f"  Error converting {path}: {e}")
    
    print("\nTesting path conversion in docker_utils.py:")
    for path in test_paths:
        try:
            converted = docker_convert_win_path_to_docker(path)
            print(f"  Original: {path}")
            print(f"  Converted: {converted}")
        except Exception as e:
            print(f"  Error converting {path}: {e}")

if __name__ == "__main__":
    test_path_conversion() 