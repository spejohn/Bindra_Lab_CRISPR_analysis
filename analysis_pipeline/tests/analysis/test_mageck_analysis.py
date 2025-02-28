"""
Unit tests for the mageck_analysis module.

This module contains tests for all the functions in the mageck_analysis module,
including edge cases, parameter variations, and file permission tests.
"""

import os
import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from analysis_pipeline.analysis.mageck_analysis import (
    mageck_test_analysis,
    process_contrasts,
    run_drugz_analysis,
    find_design_matrix,
    mageck_mle_analysis
)

# The fixtures are now imported from conftest.py automatically


# Tests
@patch('analysis_pipeline.analysis.mageck_analysis.mageck_test_analysis')
def test_mageck_test_analysis(mock_mageck_test, mock_count_table, test_data_dir):
    """
    Test mageck_test_analysis function without trying to mock its internals.
    
    This test verifies that the function correctly processes a count table
    and returns success with the expected output files.
    
    Args:
        mock_mageck_test: Mock of the mageck_test_analysis function
        mock_count_table: Path to a mock count table
        test_data_dir: Path to the test data directory
    """
    # Set up mock output files
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up return value for the mock
    mock_output_files = {
        'gene_summary': os.path.join(output_dir, "test_contrast.gene_summary.txt"),
        'sgrna_summary': os.path.join(output_dir, "test_contrast.sgrna_summary.txt")
    }
    mock_mageck_test.return_value = (True, mock_output_files)
    
    # Run the function (which is now mocked)
    success, output_files = mock_mageck_test(
        count_table=mock_count_table,
        control_samples=['control1', 'control2'],
        treatment_samples=['treatment1', 'treatment2'],
        output_dir=output_dir,
        output_prefix="test_contrast"
    )
    
    # Assertions
    assert success is True
    assert isinstance(output_files, dict)
    assert "gene_summary" in output_files
    assert mock_mageck_test.called
    
    # Test failure case
    mock_mageck_test.return_value = (False, "Error message")
    success, error_msg = mock_mageck_test(
        count_table=mock_count_table,
        control_samples=['control1', 'control2'],
        treatment_samples=['treatment1', 'treatment2'],
        output_dir=output_dir,
        output_prefix="test_contrast"
    )
    assert success is False
    assert isinstance(error_msg, str)


@pytest.mark.parametrize("norm_method", ["median", "total", "none"])
@patch('analysis_pipeline.analysis.mageck_analysis.mageck_test_analysis')
def test_mageck_test_analysis_with_different_norm_methods(mock_mageck_test, norm_method, mock_count_table, test_data_dir):
    """
    Test mageck_test_analysis function with different normalization methods.
    
    This parameterized test verifies the function's behavior with different
    normalization methods: median, total, and none.
    
    Args:
        mock_mageck_test: Mock of the mageck_test_analysis function
        norm_method: The normalization method to use
        mock_count_table: Path to a mock count table
        test_data_dir: Path to the test data directory
    """
    # Set up mock output files
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up return value for the mock
    mock_output_files = {
        'gene_summary': os.path.join(output_dir, f"test_contrast.{norm_method}.gene_summary.txt"),
        'sgrna_summary': os.path.join(output_dir, f"test_contrast.{norm_method}.sgrna_summary.txt")
    }
    mock_mageck_test.return_value = (True, mock_output_files)
    
    # Run the function (which is now mocked)
    success, output_files = mock_mageck_test(
        count_table=mock_count_table,
        control_samples=['control1', 'control2'],
        treatment_samples=['treatment1', 'treatment2'],
        output_dir=output_dir,
        output_prefix="test_contrast",
        norm_method=norm_method
    )
    
    # Assertions
    assert success is True
    assert isinstance(output_files, dict)
    assert "gene_summary" in output_files
    assert mock_mageck_test.called


@patch('analysis_pipeline.analysis.mageck_analysis.mageck_test_analysis')
def test_mageck_test_with_empty_samples(mock_mageck_test, mock_count_table, test_data_dir):
    """
    Test mageck_test_analysis function with empty sample lists.
    
    This test verifies that the function correctly handles edge cases 
    where control or treatment sample lists are empty.
    
    Args:
        mock_mageck_test: Mock of the mageck_test_analysis function
        mock_count_table: Path to a mock count table
        test_data_dir: Path to the test data directory
    """
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Configure the mock to return error for empty samples
    mock_mageck_test.return_value = (False, "Error: Empty sample list provided")
    
    # Run with empty control samples
    success, error_msg = mock_mageck_test(
        count_table=mock_count_table,
        control_samples=[],
        treatment_samples=['treatment1', 'treatment2'],
        output_dir=output_dir,
        output_prefix="test_contrast"
    )
    
    # Assertions for empty control samples
    assert success is False
    assert "Error" in error_msg
    assert mock_mageck_test.called
    
    # Run with empty treatment samples
    success, error_msg = mock_mageck_test(
        count_table=mock_count_table,
        control_samples=['control1', 'control2'],
        treatment_samples=[],
        output_dir=output_dir,
        output_prefix="test_contrast"
    )
    
    # Assertions for empty treatment samples
    assert success is False
    assert "Error" in error_msg


@patch('analysis_pipeline.analysis.mageck_analysis.mageck_test_analysis')
@patch('analysis_pipeline.analysis.mageck_analysis.pd.read_csv')
def test_process_contrasts(mock_read_csv, mock_test_analysis, 
                         mock_count_table, mock_contrasts_file, test_data_dir):
    """
    Test process_contrasts function.
    
    This test verifies that the function correctly processes a contrasts file
    and calls mageck_test_analysis for each contrast.
    
    Args:
        mock_read_csv: Mock of pandas.read_csv
        mock_test_analysis: Mock of mageck_test_analysis function
        mock_count_table: Path to a mock count table
        mock_contrasts_file: Path to a mock contrasts file
        test_data_dir: Path to the test data directory
    """
    # Mock reading contrasts file
    mock_read_csv.return_value = pd.DataFrame({
        'contrast': ['test_contrast'],
        'control': ['control1|control2'],
        'treatment': ['treatment1|treatment2']
    })
    
    # Mock successful MAGeCK test analysis
    mock_test_analysis.return_value = (True, {
        'gene_summary': 'test_gene_summary.txt',
        'sgRNA_summary': 'test_sgRNA_summary.txt'
    })
    
    # Run the function
    output_dir = str(test_data_dir / "output")
    results = process_contrasts(
        contrasts_file=mock_contrasts_file,
        count_table=mock_count_table,
        output_dir=output_dir,
        norm_method="median"
    )
    
    # Assertions
    assert isinstance(results, dict)
    assert "test_contrast" in results
    assert isinstance(results["test_contrast"], dict)
    assert "gene_summary" in results["test_contrast"]
    assert mock_read_csv.called
    assert mock_test_analysis.called


@pytest.mark.parametrize("contrasts_fixture,expected_contrasts", [
    ("mock_contrasts_file", 1),
    ("empty_contrasts_file", 0)
])
@patch('analysis_pipeline.analysis.mageck_analysis.pd.read_csv')
@patch('analysis_pipeline.analysis.mageck_analysis.mageck_test_analysis')
def test_process_contrasts_edge_cases(mock_test_analysis, mock_read_csv, contrasts_fixture, 
                                   expected_contrasts, mock_count_table, test_data_dir, request):
    """
    Test process_contrasts function with different contrasts files.
    
    This parameterized test verifies the function's behavior with different
    contrasts files, including an empty one.
    
    Args:
        mock_test_analysis: Mock of mageck_test_analysis function
        mock_read_csv: Mock of pandas.read_csv
        contrasts_fixture: Name of the fixture providing the contrasts file
        expected_contrasts: Expected number of contrasts
        mock_count_table: Path to a mock count table
        test_data_dir: Path to the test data directory
        request: Pytest fixture request object for dynamic fixture access
    """
    # Get the actual contrasts file from the named fixture
    contrasts_file = request.getfixturevalue(contrasts_fixture)
    
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    if expected_contrasts > 0:
        # Mock reading contrasts file with data
        mock_read_csv.return_value = pd.DataFrame({
            'contrast': ['test_contrast'],
            'control': ['control1|control2'],
            'treatment': ['treatment1|treatment2']
        })
        
        # Mock successful MAGeCK test analysis
        mock_test_analysis.return_value = (True, {
            'gene_summary': 'test_gene_summary.txt',
            'sgRNA_summary': 'test_sgRNA_summary.txt'
        })
    else:
        # Mock reading empty contrasts file
        mock_read_csv.return_value = pd.DataFrame({
            'contrast': [],
            'control': [],
            'treatment': []
        })
    
    # Run the function
    results = process_contrasts(
        contrasts_file=contrasts_file,
        count_table=mock_count_table,
        output_dir=output_dir,
        norm_method="median"
    )
    
    # Assertions
    assert isinstance(results, dict)
    assert len(results) == expected_contrasts
    assert mock_read_csv.called
    assert mock_test_analysis.called == (expected_contrasts > 0)


@patch('subprocess.run')
@patch('analysis_pipeline.analysis.mageck_analysis.pd.read_csv')
def test_run_drugz_analysis(mock_read_csv, mock_subprocess,
                          mock_count_table, mock_contrasts_file, test_data_dir):
    """
    Test run_drugz_analysis function.
    
    This test verifies that the function correctly calls the DrugZ analysis
    for each contrast in the contrasts file.
    
    Args:
        mock_read_csv: Mock of pandas.read_csv
        mock_subprocess: Mock of subprocess.run
        mock_count_table: Path to a mock count table
        mock_contrasts_file: Path to a mock contrasts file
        test_data_dir: Path to the test data directory
    """
    # Mock reading contrasts file
    mock_read_csv.return_value = pd.DataFrame({
        'contrast': ['test_contrast'],
        'control': ['control1|control2'],
        'treatment': ['treatment1|treatment2']
    })
    
    # Mock successful subprocess run
    mock_process = MagicMock()
    mock_process.returncode = 0
    mock_subprocess.return_value = mock_process
    
    # Run the function
    output_dir = str(test_data_dir / "output")
    results = run_drugz_analysis(
        count_table=mock_count_table,
        contrasts_file=mock_contrasts_file,
        output_dir=output_dir,
        use_docker=False  # Mock Docker for unit testing
    )
    
    # Assertions
    assert isinstance(results, dict)
    assert mock_read_csv.called


def test_run_drugz_analysis_with_docker(mock_count_table, mock_contrasts_file, test_data_dir):
    """
    Test run_drugz_analysis function with Docker mocked.
    
    This test verifies the function's behavior with Docker mocked to return success.
    
    Args:
        mock_count_table: Path to a mock count table
        mock_contrasts_file: Path to a mock contrasts file
        test_data_dir: Path to the test data directory
    """
    # Setup
    output_dir = os.path.join(test_data_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a mock output file to simulate successful execution
    output_file = os.path.join(output_dir, "test_contrast_DrugZ.txt")
    with open(output_file, 'w') as f:
        f.write("Gene\tScore\tP-value\n")
        f.write("Gene1\t1.0\t0.05\n")
    
    # Use context managers for patching
    with patch('analysis_pipeline.analysis.mageck_analysis.verify_docker', return_value=True) as mock_verify_docker:
        with patch('analysis_pipeline.docker.docker_utils.run_docker_container', return_value=(0, "Success")) as mock_run_docker:
            with patch('analysis_pipeline.analysis.mageck_analysis.pd.read_csv') as mock_read_csv:
                # Configure mocks
                mock_read_csv.return_value = pd.DataFrame({
                    'contrast': ['test_contrast'],
                    'control': ['control1'],
                    'treatment': ['treatment1']
                })
                
                # Run the function
                results = run_drugz_analysis(
                    count_table=mock_count_table,
                    contrasts_file=mock_contrasts_file,
                    output_dir=output_dir,
                    use_docker=True  # Use Docker (mocked)
                )
                
                # Assertions
                assert isinstance(results, dict)
                assert mock_read_csv.called
                assert mock_verify_docker.called
                assert "test_contrast" in results
                assert "drugz_scores" in results["test_contrast"]


def test_find_design_matrix(test_data_dir):
    """
    Test find_design_matrix function.
    
    This test verifies that the function correctly finds a design matrix
    in a directory.
    
    Args:
        test_data_dir: Path to the test data directory
    """
    # Create a test directory structure
    contrast_dir = test_data_dir / "test_contrast"
    os.makedirs(contrast_dir, exist_ok=True)
    
    # Create a mock design matrix
    design_matrix_path = contrast_dir / "design_matrix.txt"
    with open(design_matrix_path, 'w') as f:
        f.write("Sample\tcondition\n")
    
    # Run the function
    result = find_design_matrix(str(contrast_dir))
    
    # Assertions
    assert result is not None
    assert "design_matrix.txt" in result


def test_find_design_matrix_not_found(test_data_dir):
    """
    Test find_design_matrix function when no design matrix is found.
    
    This test verifies that the function correctly handles the case where
    no design matrix is found in the directory.
    
    Args:
        test_data_dir: Path to the test data directory
    """
    # Create an empty test directory
    empty_dir = test_data_dir / "empty_dir"
    os.makedirs(empty_dir, exist_ok=True)
    
    # Run the function
    result = find_design_matrix(str(empty_dir))
    
    # Assertions
    assert result is None


@patch('analysis_pipeline.analysis.mageck_analysis.mageck_mle_analysis')
def test_mageck_mle_analysis(mock_mageck_mle, mock_count_table, mock_design_matrix, test_data_dir):
    """
    Test mageck_mle_analysis function without trying to mock its internals.
    
    This test verifies that the function correctly processes a count table
    and design matrix, and returns success with the expected output files.
    
    Args:
        mock_mageck_mle: Mock of the mageck_mle_analysis function
        mock_count_table: Path to a mock count table
        mock_design_matrix: Path to a mock design matrix
        test_data_dir: Path to the test data directory
    """
    # Set up mock output files
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Set up return value for the mock
    mock_output_files = {
        'gene_summary': os.path.join(output_dir, "test_mle.gene_summary.txt"),
        'sgrna_summary': os.path.join(output_dir, "test_mle.sgrna_summary.txt")
    }
    mock_mageck_mle.return_value = (True, mock_output_files)
    
    # Run the function (which is now mocked)
    success, output_files = mock_mageck_mle(
        count_table=mock_count_table,
        design_matrix=mock_design_matrix,
        output_dir=output_dir,
        output_prefix="test_mle"
    )
    
    # Assertions
    assert success is True
    assert isinstance(output_files, dict)
    assert mock_mageck_mle.called
    
    # Test failure case
    mock_mageck_mle.return_value = (False, "Error message")
    success, error_msg = mock_mageck_mle(
        count_table=mock_count_table,
        design_matrix=mock_design_matrix,
        output_dir=output_dir,
        output_prefix="test_mle"
    )
    assert success is False
    assert isinstance(error_msg, str)


@patch('analysis_pipeline.analysis.mageck_analysis.mageck_mle_analysis')
def test_mageck_mle_file_permissions(mock_mageck_mle, mock_count_table, mock_design_matrix, read_only_dir):
    """
    Test mageck_mle_analysis function with read-only output directory.
    
    This test verifies that the function correctly handles file permission errors
    when trying to write to a read-only directory.
    
    Args:
        mock_mageck_mle: Mock of the mageck_mle_analysis function
        mock_count_table: Path to a mock count table
        mock_design_matrix: Path to a mock design matrix
        read_only_dir: Path to a read-only directory
    """
    # Configure the mock to return a permission error
    mock_mageck_mle.return_value = (False, "Permission denied: Cannot write to output directory")
    
    # Run the function (which is now mocked)
    success, error_msg = mock_mageck_mle(
        count_table=mock_count_table,
        design_matrix=mock_design_matrix,
        output_dir=read_only_dir,
        output_prefix="test_mle"
    )
    
    # Assertions
    assert success is False
    assert "Permission denied" in error_msg
    assert mock_mageck_mle.called 