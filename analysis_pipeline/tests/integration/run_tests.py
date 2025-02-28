#!/usr/bin/env python3
"""
Script to run the integration tests for the CRISPR analysis pipeline.
"""

import os
import sys
import subprocess
from pathlib import Path

def run_integration_tests():
    """Run the integration tests."""
    print("Running integration tests with fixed output directory paths...")
    
    # Get the path to the test runners
    script_dir = Path(__file__).parent
    test_runners_dir = script_dir / "test_runners"
    
    # List of test files to run
    test_files = [
        "test_full_workflow.py",
        "test_incremental_workflow.py",
        "test_error_handling.py"
    ]
    
    # Run each test file
    for test_file in test_files:
        test_path = test_runners_dir / test_file
        if not test_path.exists():
            print(f"Error: Test file {test_file} not found.")
            continue
            
        print(f"\nRunning test: {test_file}")
        cmd = ["python", "-m", "pytest", str(test_path), "-v"]
        
        # Run the test
        try:
            subprocess.run(cmd, check=True)
            print(f"Test {test_file} passed.")
        except subprocess.CalledProcessError:
            print(f"Test {test_file} failed.")
            
    print("\nIntegration tests completed.")

if __name__ == "__main__":
    run_integration_tests() 