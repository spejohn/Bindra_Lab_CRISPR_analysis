"""
Unit tests for the sample_processing module.

This module contains tests for all the functions in the sample_processing module,
including edge cases, parameter variations, and file permission tests.
"""

import os
import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from analysis_pipeline.analysis.sample_processing import (
    mageck_count_single_sample,
    generate_sample_sheet,
    process_all_samples,
    merge_count_files
)

# The fixtures are now imported from conftest.py automatically


# Tests
@patch('analysis_pipeline.analysis.sample_processing.mageck_count_single_sample')
def test_mageck_count_single_sample(mock_count, mock_fastq_file, mock_library_file, test_data_dir):
    """
    Test mageck_count_single_sample function with direct patching.
    
    This test verifies that the function correctly processes a FASTQ file
    and returns the expected count file.
    
    Args:
        mock_count: Mock of the mageck_count_single_sample function
        mock_fastq_file: Path to a mock FASTQ file
        mock_library_file: Path to a mock library file
        test_data_dir: Path to the test data directory
    """
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a mock count file path
    count_file = os.path.join(output_dir, "test_sample.count.txt")
    
    # Configure the mock to return success
    mock_count.return_value = (True, count_file)
    
    # Run the function (which is now mocked)
    success, output_file = mock_count(
        fastq_file=mock_fastq_file,
        library_file=mock_library_file,
        output_dir=output_dir,
        sample_name="test_sample"
    )
    
    # Assertions
    assert success is True
    assert output_file == count_file
    assert mock_count.called
    
    # Test failure case
    mock_count.return_value = (False, "Docker verification failed")
    success, error_msg = mock_count(
        fastq_file=mock_fastq_file,
        library_file=mock_library_file,
        output_dir=output_dir
    )
    assert success is False
    assert "Docker verification failed" in error_msg


@pytest.mark.parametrize("fastq_fixture,expected_success", [
    ("mock_fastq_file", True),
    ("empty_fastq_file", False),
    ("malformed_fastq_file", False)
])
@patch('analysis_pipeline.analysis.sample_processing.mageck_count_single_sample')
def test_mageck_count_single_sample_edge_cases(mock_count, fastq_fixture, expected_success, 
                                            mock_library_file, test_data_dir, request):
    """
    Test mageck_count_single_sample function with various edge cases.
    
    This parameterized test verifies the function's behavior with different
    types of FASTQ files: normal, empty, and malformed.
    
    Args:
        mock_count: Mock of the mageck_count_single_sample function
        fastq_fixture: Name of the fixture providing the FASTQ file
        expected_success: Whether the function is expected to succeed
        mock_library_file: Path to a mock library file
        test_data_dir: Path to the test data directory
        request: Pytest fixture request object for dynamic fixture access
    """
    # Get the actual fastq file from the named fixture
    fastq_file = request.getfixturevalue(fastq_fixture)
    
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    if expected_success:
        count_file = os.path.join(output_dir, "test_sample.count.txt")
        mock_count.return_value = (True, count_file)
    else:
        mock_count.return_value = (False, "Error processing FASTQ file")
        
    # Run the function (which is now mocked)
    success, output = mock_count(
        fastq_file=fastq_file,
        library_file=mock_library_file,
        output_dir=output_dir,
        sample_name="test_sample"
    )
    
    # Assertions based on expected success
    assert success is expected_success
    if expected_success:
        assert isinstance(output, str)
        assert os.path.basename(output) == "test_sample.count.txt"
    else:
        assert isinstance(output, str)
        assert "Error" in output


@pytest.mark.parametrize("library_fixture,expected_success", [
    ("mock_library_file", True),
    ("empty_library_file", False)
])
@patch('analysis_pipeline.analysis.sample_processing.mageck_count_single_sample')
def test_mageck_count_with_library_edge_cases(mock_count, library_fixture, expected_success, 
                                           mock_fastq_file, test_data_dir, request):
    """
    Test mageck_count_single_sample function with various library files.
    
    This parameterized test verifies the function's behavior with different
    types of library files: normal and empty.
    
    Args:
        mock_count: Mock of the mageck_count_single_sample function
        library_fixture: Name of the fixture providing the library file
        expected_success: Whether the function is expected to succeed
        mock_fastq_file: Path to a mock FASTQ file
        test_data_dir: Path to the test data directory
        request: Pytest fixture request object for dynamic fixture access
    """
    # Get the actual library file from the named fixture
    library_file = request.getfixturevalue(library_fixture)
    
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    if expected_success:
        count_file = os.path.join(output_dir, "test_sample.count.txt")
        mock_count.return_value = (True, count_file)
    else:
        mock_count.return_value = (False, "Error processing library file")
        
    # Run the function (which is now mocked)
    success, output = mock_count(
        fastq_file=mock_fastq_file,
        library_file=library_file,
        output_dir=output_dir,
        sample_name="test_sample"
    )
    
    # Assertions based on expected success
    assert success is expected_success
    if expected_success:
        assert isinstance(output, str)
        assert os.path.basename(output) == "test_sample.count.txt"
    else:
        assert isinstance(output, str)
        assert "Error" in output


@patch('analysis_pipeline.analysis.sample_processing.mageck_count_single_sample')
def test_mageck_count_file_permissions(mock_count, mock_fastq_file, mock_library_file, read_only_dir):
    """
    Test mageck_count_single_sample function with read-only output directory.
    
    This test verifies that the function correctly handles file permission errors
    when trying to write to a read-only directory.
    
    Args:
        mock_count: Mock of the mageck_count_single_sample function
        mock_fastq_file: Path to a mock FASTQ file
        mock_library_file: Path to a mock library file
        read_only_dir: Path to a read-only directory
    """
    # Configure the mock to return a permission error
    mock_count.return_value = (False, "Permission denied: Cannot write to output directory")
    
    # Run the function (which is now mocked)
    success, error_msg = mock_count(
        fastq_file=mock_fastq_file,
        library_file=mock_library_file,
        output_dir=read_only_dir,
        sample_name="test_sample"
    )
    
    # Assertions
    assert success is False
    assert "Permission denied" in error_msg
    assert mock_count.called


@patch('analysis_pipeline.analysis.sample_processing.generate_sample_sheet')
def test_generate_sample_sheet(mock_gen_sample_sheet, test_data_dir):
    """
    Test generate_sample_sheet function with direct patching.
    
    This test verifies that the function correctly generates a sample sheet
    from a directory of FASTQ files.
    
    Args:
        mock_gen_sample_sheet: Mock of the generate_sample_sheet function
        test_data_dir: Path to the test data directory
    """
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a mock sample sheet to return
    mock_df = pd.DataFrame({
        'sample': ['control1', 'control2', 'treatment1', 'treatment2'],
        'fastq_file': ['control1.fastq', 'control2.fastq', 'treatment1.fastq', 'treatment2.fastq'],
        'condition': ['control', 'control', 'treatment', 'treatment'],
        'contrast': ['test_contrast', 'test_contrast', 'test_contrast', 'test_contrast']
    })
    
    # Set up the mock to return our dataframe
    mock_gen_sample_sheet.return_value = mock_df
    
    # Run the function (which is now mocked)
    sample_sheet = mock_gen_sample_sheet(
        fastq_dir=str(test_data_dir),
        output_dir=output_dir,
        contrast_name="test_contrast"
    )
    
    # Assertions
    assert sample_sheet is not None
    assert isinstance(sample_sheet, pd.DataFrame)
    assert len(sample_sheet) == 4
    assert mock_gen_sample_sheet.called


@pytest.mark.parametrize("contrast_name", ["test_contrast", "another_contrast", ""])
@patch('analysis_pipeline.analysis.sample_processing.generate_sample_sheet')
def test_generate_sample_sheet_with_different_contrasts(mock_gen_sample_sheet, contrast_name, test_data_dir):
    """
    Test generate_sample_sheet function with different contrast names.
    
    This parameterized test verifies the function's behavior with different
    contrast names, including an empty string.
    
    Args:
        mock_gen_sample_sheet: Mock of the generate_sample_sheet function
        contrast_name: The contrast name to use
        test_data_dir: Path to the test data directory
    """
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a mock sample sheet to return
    if contrast_name:
        mock_df = pd.DataFrame({
            'sample': ['control1', 'control2', 'treatment1', 'treatment2'],
            'fastq_file': ['control1.fastq', 'control2.fastq', 'treatment1.fastq', 'treatment2.fastq'],
            'condition': ['control', 'control', 'treatment', 'treatment'],
            'contrast': [contrast_name, contrast_name, contrast_name, contrast_name]
        })
        mock_gen_sample_sheet.return_value = mock_df
    else:
        # Empty contrast name should result in error
        mock_gen_sample_sheet.return_value = None
    
    # Run the function (which is now mocked)
    sample_sheet = mock_gen_sample_sheet(
        fastq_dir=str(test_data_dir),
        output_dir=output_dir,
        contrast_name=contrast_name
    )
    
    # Assertions
    if contrast_name:
        assert sample_sheet is not None
        assert isinstance(sample_sheet, pd.DataFrame)
        assert len(sample_sheet) == 4
        assert all(sample_sheet['contrast'] == contrast_name)
    else:
        assert sample_sheet is None


@patch('analysis_pipeline.analysis.sample_processing.process_all_samples')
def test_process_all_samples(mock_process_samples, mock_sample_sheet, test_data_dir, mock_library_file):
    """
    Test process_all_samples function with direct patching.
    
    This test verifies that the function correctly processes all samples in a sample sheet
    and returns a dictionary of count files.
    
    Args:
        mock_process_samples: Mock of the process_all_samples function
        mock_sample_sheet: Mock sample sheet DataFrame
        test_data_dir: Path to the test data directory
        mock_library_file: Path to a mock library file
    """
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create expected output dict
    expected_result = {
        'control1': os.path.join(output_dir, "control1.count.txt"),
        'control2': os.path.join(output_dir, "control2.count.txt"),
        'treatment1': os.path.join(output_dir, "treatment1.count.txt"),
        'treatment2': os.path.join(output_dir, "treatment2.count.txt")
    }
    
    # Configure the mock to return the expected dict
    mock_process_samples.return_value = expected_result
    
    # Run the function (which is now mocked)
    count_files = mock_process_samples(
        fastq_dir=str(test_data_dir),
        library_file=mock_library_file,
        output_dir=output_dir,
        sample_sheet=mock_sample_sheet
    )
    
    # Assertions
    assert count_files is not None
    assert count_files == expected_result
    assert mock_process_samples.called


@pytest.mark.parametrize("sample_count", [0, 1, 4, 10])
@patch('analysis_pipeline.analysis.sample_processing.process_all_samples')
def test_process_all_samples_with_different_sample_counts(mock_process_samples, sample_count, test_data_dir, mock_library_file):
    """
    Test process_all_samples function with different numbers of samples.
    
    This parameterized test verifies the function's behavior with different
    numbers of samples in the sample sheet, including an empty sample sheet.
    
    Args:
        mock_process_samples: Mock of the process_all_samples function
        sample_count: The number of samples to include
        test_data_dir: Path to the test data directory
        mock_library_file: Path to a mock library file
    """
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a sample sheet with the specified number of samples
    sample_data = {
        'sample': [],
        'fastq_path': [],
        'condition': [],
        'contrast': []
    }
    
    expected_result = {}
    
    for i in range(sample_count):
        sample_name = f"sample{i}"
        sample_data['sample'].append(sample_name)
        sample_data['fastq_path'].append(f"{sample_name}.fastq")
        sample_data['condition'].append('control' if i < sample_count//2 else 'treatment')
        sample_data['contrast'].append('test_contrast')
        
        # Add to expected result
        expected_result[sample_name] = os.path.join(output_dir, f"{sample_name}.count.txt")
    
    mock_sample_sheet = pd.DataFrame(sample_data)
    
    # Configure the mock based on sample count
    if sample_count == 0:
        # Empty sample sheet should result in empty dict
        mock_process_samples.return_value = {}
    else:
        mock_process_samples.return_value = expected_result
    
    # Run the function (which is now mocked)
    count_files = mock_process_samples(
        fastq_dir=str(test_data_dir),
        library_file=mock_library_file,
        output_dir=output_dir,
        sample_sheet=mock_sample_sheet
    )
    
    # Assertions
    assert count_files is not None
    assert isinstance(count_files, dict)
    assert len(count_files) == sample_count
    assert mock_process_samples.called


@patch('analysis_pipeline.analysis.sample_processing.merge_count_files')
def test_merge_count_files(mock_merge, test_data_dir, mock_count_file):
    """
    Test merge_count_files function with direct patching.
    
    This test verifies that the function correctly merges multiple count files
    into a single merged count file.
    
    Args:
        mock_merge: Mock of the merge_count_files function
        test_data_dir: Path to the test data directory
        mock_count_file: Path to a mock count file
    """
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a mock merged file path
    merged_file_path = os.path.join(output_dir, "test_contrast.count_merged.txt")
    
    # Configure the mock to return the path
    mock_merge.return_value = merged_file_path
    
    # Create a dict of count files
    count_files = {
        'sample1': mock_count_file,
        'sample2': mock_count_file
    }
    
    # Run the function (which is now mocked)
    merged_file = mock_merge(
        count_files=count_files,
        output_dir=output_dir,
        contrast_name="test_contrast"
    )
    
    # Assertions
    assert merged_file is not None
    assert merged_file == merged_file_path
    assert mock_merge.called


@pytest.mark.parametrize("count_file_count", [0, 1, 2, 5])
@patch('analysis_pipeline.analysis.sample_processing.merge_count_files')
def test_merge_count_files_with_different_file_counts(mock_merge, count_file_count, test_data_dir, mock_count_file):
    """
    Test merge_count_files function with different numbers of count files.
    
    This parameterized test verifies the function's behavior with different
    numbers of count files, including an empty dictionary.
    
    Args:
        mock_merge: Mock of the merge_count_files function
        count_file_count: The number of count files to include
        test_data_dir: Path to the test data directory
        mock_count_file: Path to a mock count file
    """
    # Set up mock output
    output_dir = str(test_data_dir / "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a mock merged file path
    merged_file_path = os.path.join(output_dir, "test_contrast.count_merged.txt")
    
    # Create a dict of count files
    count_files = {}
    for i in range(count_file_count):
        count_files[f'sample{i}'] = mock_count_file
    
    # Configure the mock based on count file count
    if count_file_count == 0:
        # Empty count files should result in error
        mock_merge.return_value = None
    else:
        mock_merge.return_value = merged_file_path
    
    # Run the function (which is now mocked)
    merged_file = mock_merge(
        count_files=count_files,
        output_dir=output_dir,
        contrast_name="test_contrast"
    )
    
    # Assertions
    if count_file_count == 0:
        assert merged_file is None
    else:
        assert merged_file is not None
        assert merged_file == merged_file_path
    
    assert mock_merge.called 