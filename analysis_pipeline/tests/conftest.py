"""
Shared fixtures for all tests in the analysis pipeline.
This file contains fixtures that can be reused across multiple test files.
"""

import os
import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock


@pytest.fixture
def test_data_dir():
    """
    Fixture for test data directory.
    Creates a test data directory within tests and ensures it exists.
    
    Returns:
        Path: Path object pointing to the test data directory
    """
    # Create a test data directory within tests
    test_dir = Path(__file__).parent / "test_data"
    os.makedirs(test_dir, exist_ok=True)
    return test_dir


@pytest.fixture
def mock_fastq_file(test_data_dir):
    """
    Create a mock FASTQ file for testing.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the mock FASTQ file
    """
    fastq_path = test_data_dir / "test_sample.fastq"
    # Create a minimal FASTQ file for testing
    with open(fastq_path, 'w') as f:
        f.write("@SEQ_ID\nACGTACGT\n+\n!!!!!!!!\n")
    return str(fastq_path)


@pytest.fixture
def empty_fastq_file(test_data_dir):
    """
    Create an empty FASTQ file for testing edge cases.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the empty FASTQ file
    """
    fastq_path = test_data_dir / "empty_sample.fastq"
    # Create an empty file
    with open(fastq_path, 'w') as f:
        pass
    return str(fastq_path)


@pytest.fixture
def malformed_fastq_file(test_data_dir):
    """
    Create a malformed FASTQ file for testing error handling.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the malformed FASTQ file
    """
    fastq_path = test_data_dir / "malformed_sample.fastq"
    # Create a file that doesn't follow FASTQ format
    with open(fastq_path, 'w') as f:
        f.write("This is not a valid FASTQ file format\n")
    return str(fastq_path)


@pytest.fixture
def mock_library_file(test_data_dir):
    """
    Create a mock library file for testing.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the mock library file
    """
    library_path = test_data_dir / "test_library.txt"
    # Create a minimal library file
    with open(library_path, 'w') as f:
        f.write("sgRNA\tGene\nsgRNA1\tGeneA\nsgRNA2\tGeneB\n")
    return str(library_path)


@pytest.fixture
def empty_library_file(test_data_dir):
    """
    Create an empty library file for testing edge cases.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the empty library file
    """
    library_path = test_data_dir / "empty_library.txt"
    # Create an empty file
    with open(library_path, 'w') as f:
        pass
    return str(library_path)


@pytest.fixture
def mock_count_file(test_data_dir):
    """
    Create a mock count file for testing.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the mock count file
    """
    count_path = test_data_dir / "test_counts.txt"
    # Create a minimal count file
    with open(count_path, 'w') as f:
        f.write("sgRNA\tGene\tCount\nsgRNA1\tGeneA\t10\nsgRNA2\tGeneB\t20\n")
    return str(count_path)


@pytest.fixture
def mock_count_table(test_data_dir):
    """
    Create a mock count table for testing.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the mock count table
    """
    count_path = test_data_dir / "test_count_table.txt"
    # Create a minimal count table
    df = pd.DataFrame({
        'sgRNA': ['sgRNA1', 'sgRNA2'],
        'Gene': ['GeneA', 'GeneB'],
        'control1': [10, 20],
        'control2': [15, 25],
        'treatment1': [5, 40],
        'treatment2': [8, 35]
    })
    df.to_csv(count_path, sep='\t', index=False)
    return str(count_path)


@pytest.fixture
def mock_contrasts_file(test_data_dir):
    """
    Create a mock contrasts file for testing.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the mock contrasts file
    """
    contrasts_path = test_data_dir / "test_contrasts.txt"
    # Create a minimal contrasts file
    with open(contrasts_path, 'w') as f:
        f.write("contrast,control,treatment\n")
        f.write("test_contrast,control1|control2,treatment1|treatment2\n")
    return str(contrasts_path)


@pytest.fixture
def empty_contrasts_file(test_data_dir):
    """
    Create an empty contrasts file for testing edge cases.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the empty contrasts file
    """
    contrasts_path = test_data_dir / "empty_contrasts.txt"
    # Create a file with header only
    with open(contrasts_path, 'w') as f:
        f.write("contrast,control,treatment\n")
    return str(contrasts_path)


@pytest.fixture
def mock_design_matrix(test_data_dir):
    """
    Create a mock design matrix for testing.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the mock design matrix
    """
    design_path = test_data_dir / "test_design_matrix.txt"
    # Create a minimal design matrix
    with open(design_path, 'w') as f:
        f.write("Sample\tcondition\n")
        f.write("control1\tcontrol\n")
        f.write("control2\tcontrol\n")
        f.write("treatment1\ttreatment\n")
        f.write("treatment2\ttreatment\n")
    return str(design_path)


@pytest.fixture
def mock_sample_sheet():
    """
    Create a mock sample sheet DataFrame for testing.
    
    Returns:
        pd.DataFrame: A mock sample sheet DataFrame
    """
    sample_data = {
        'sample': ['control1', 'control2', 'treatment1', 'treatment2'],
        'fastq_path': ['control1.fastq', 'control2.fastq', 'treatment1.fastq', 'treatment2.fastq'],
        'condition': ['control', 'control', 'treatment', 'treatment'],
        'contrast': ['test_contrast', 'test_contrast', 'test_contrast', 'test_contrast']
    }
    return pd.DataFrame(sample_data)


@pytest.fixture
def read_only_dir(test_data_dir):
    """
    Create a read-only directory for testing file permission errors.
    
    Args:
        test_data_dir (Path): The test data directory fixture
        
    Returns:
        str: Path to the read-only directory
    """
    readonly_dir = test_data_dir / "readonly"
    os.makedirs(readonly_dir, exist_ok=True)
    
    # Create a dummy file inside
    with open(readonly_dir / "dummy.txt", 'w') as f:
        f.write("Test file")
    
    # On Unix-like systems, make the directory read-only
    try:
        import stat
        os.chmod(readonly_dir, stat.S_IRUSR | stat.S_IXUSR)
    except:
        # If on Windows or permission denied, skip making it read-only
        pass
        
    return str(readonly_dir) 