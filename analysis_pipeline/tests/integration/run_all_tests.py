#!/usr/bin/env python3
"""
Run all integration tests for the CRISPR analysis pipeline.

This script runs all integration tests and reports a summary of results.
"""

import os
import pytest
import sys
import time
from pathlib import Path

def main():
    """
    Main function to run all integration tests and report results.
    """
    start_time = time.time()
    
    # Directory containing this script
    test_dir = Path(__file__).parent
    
    # Add integration test directory to Python path
    if str(test_dir) not in sys.path:
        sys.path.insert(0, str(test_dir))
    
    # Set up test arguments
    pytest_args = [
        "-xvs",              # Verbose, stop on first failure, show output
        "--no-header",       # Don't show pytest header
        "--tb=native",       # Use native traceback style
        str(test_dir / "test_runners")  # Test directory
    ]
    
    # Run the tests
    print("\n=== Running CRISPR Analysis Pipeline Integration Tests ===\n")
    result = pytest.main(pytest_args)
    
    # Calculate elapsed time
    elapsed_time = time.time() - start_time
    
    # Print summary
    print(f"\n=== Integration Test Summary ===")
    print(f"Time elapsed: {elapsed_time:.2f} seconds")
    print(f"Exit code: {result}")
    
    if result == 0:
        print("Status: All integration tests PASSED")
    else:
        print("Status: Some integration tests FAILED")
    
    return result

if __name__ == "__main__":
    sys.exit(main()) 