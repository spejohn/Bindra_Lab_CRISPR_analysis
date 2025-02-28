#!/usr/bin/env python3
"""
Utility script to run all tests for the analysis pipeline.

This script ensures the tests are run from the correct directory and 
with the proper Python path configuration.
"""

import os
import sys
import pytest
import argparse
from pathlib import Path


def run_tests(test_dir=None, verbose=False, coverage=False):
    """
    Run tests in the specified directory with pytest.
    
    Parameters
    ----------
    test_dir : str, optional
        Directory containing tests to run. If None, runs all tests.
    verbose : bool, optional
        Whether to run with verbose output.
    coverage : bool, optional
        Whether to run with coverage reporting.
    """
    # Get the repo root directory (parent of the current directory)
    repo_root = Path(__file__).parent.parent
    
    # Ensure we're in the repo root for proper imports
    os.chdir(repo_root)
    
    # Add the repo root to the Python path
    sys.path.insert(0, str(repo_root))
    
    # Build the pytest arguments
    pytest_args = []
    
    # Add the test directory
    if test_dir:
        test_path = repo_root / "tests" / test_dir
        if not test_path.exists():
            print(f"Error: Test directory '{test_dir}' not found")
            return 1
        pytest_args.append(str(test_path))
    else:
        pytest_args.append(str(repo_root / "tests"))
    
    # Add verbosity
    if verbose:
        pytest_args.append("-v")
    
    # Add coverage if requested
    if coverage:
        pytest_args.extend(["--cov=analysis_pipeline", "--cov-report=term"])
    
    # Run pytest
    return pytest.main(pytest_args)


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Run analysis pipeline tests")
    
    parser.add_argument(
        "--dir", "-d",
        help="Test directory to run (e.g., 'analysis' for tests/analysis)",
        default=None
    )
    
    parser.add_argument(
        "--verbose", "-v",
        help="Run with verbose output",
        action="store_true"
    )
    
    parser.add_argument(
        "--coverage", "-c",
        help="Run with coverage reporting",
        action="store_true"
    )
    
    args = parser.parse_args()
    
    sys.exit(run_tests(args.dir, args.verbose, args.coverage))


if __name__ == "__main__":
    main() 