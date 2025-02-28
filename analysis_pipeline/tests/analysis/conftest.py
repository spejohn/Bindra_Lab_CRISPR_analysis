"""
Conftest for analysis tests.
Contains shared fixtures and test setup/teardown.
"""

import os
import pytest
import shutil
from pathlib import Path


@pytest.fixture(scope="session")
def test_data_root():
    """Root fixture for all test data - this is shared across all tests"""
    # Create a test data directory at the module level
    test_dir = Path(__file__).parent / "test_data"
    os.makedirs(test_dir, exist_ok=True)
    
    # Create output directory
    output_dir = test_dir / "output"
    os.makedirs(output_dir, exist_ok=True)
    
    yield test_dir
    
    # Clean up after all tests are done
    # Comment out for debugging if needed
    # shutil.rmtree(test_dir)


@pytest.fixture(autouse=True)
def clean_output_dir(test_data_root):
    """Clean the output directory before each test"""
    output_dir = test_data_root / "output"
    
    # Ensure output dir exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Clean all files in the output directory but keep the directory
    for item in output_dir.iterdir():
        if item.is_file():
            item.unlink()
        elif item.is_dir():
            shutil.rmtree(item)
    
    yield
    
    # No teardown needed as the next test will clean it again 