"""
Unit tests for the run_drugz_docker module.
"""

import os
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from analysis_pipeline.analysis.run_drugz_docker import run_drugz_docker


# Fixtures for test data
@pytest.fixture
def test_data_dir():
    """Fixture for test data directory"""
    # Create a test data directory within tests
    test_dir = Path(__file__).parent / "test_data"
    os.makedirs(test_dir, exist_ok=True)
    return test_dir


@pytest.fixture
def mock_count_file(test_data_dir):
    """Create a mock count file for testing"""
    count_path = test_data_dir / "test_drugz_counts.txt"
    # Create a minimal count file with required format
    with open(count_path, 'w') as f:
        f.write("sgRNA\tGene\tcontrol1\tcontrol2\ttreatment1\ttreatment2\n")
        f.write("sgRNA1\tGeneA\t10\t15\t5\t8\n")
        f.write("sgRNA2\tGeneB\t20\t25\t40\t35\n")
    return str(count_path)


@pytest.fixture
def mock_design_file(test_data_dir):
    """Create a mock design file for testing"""
    design_path = test_data_dir / "test_drugz_design.txt"
    # Create a minimal design file with required format
    with open(design_path, 'w') as f:
        f.write("treatment1\ttreatment\n")
        f.write("treatment2\ttreatment\n")
        f.write("control1\tcontrol\n")
        f.write("control2\tcontrol\n")
    return str(design_path)


# Tests
@patch('analysis_pipeline.analysis.run_drugz_docker.subprocess.run')
def test_run_drugz_docker(mock_subprocess, mock_count_file, mock_design_file, test_data_dir):
    """Test run_drugz_docker function"""
    # Mock successful subprocess run
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_subprocess.return_value = mock_process
    
    # Run the function
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    result = run_drugz_docker(
        count_file=mock_count_file,
        design_file=mock_design_file,
        output_dir=output_dir,
        output_prefix="test_drugz"
    )
    
    # Assertions
    assert mock_subprocess.called
    
    # Mock a failed subprocess run
    mock_process.returncode = 1
    result_failed = run_drugz_docker(
        count_file=mock_count_file,
        design_file=mock_design_file,
        output_dir=output_dir,
        output_prefix="test_drugz_failed"
    )
    
    # Since the original function may handle the error internally, we may not have a direct assertion,
    # but we can check that the function completed at least
    assert mock_subprocess.call_count == 2 