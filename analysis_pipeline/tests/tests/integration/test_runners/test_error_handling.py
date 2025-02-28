#!/usr/bin/env python3
"""
Integration test for error handling in the CRISPR analysis workflow.

This test verifies that the workflow correctly handles problematic input.
"""

import os
import pytest
import stat
from pathlib import Path


# Test error handling for missing library file
def test_missing_library_file(setup_test_workflow, run_snakemake):
    """
    Test that the workflow correctly handles a missing library file.
    
    Args:
        setup_test_workflow: Fixture that sets up the test workflow
        run_snakemake: Fixture that provides a function to run Snakemake
    """
    # Get the test workflow setup
    workflow_setup = setup_test_workflow
    
    # Create test samples but intentionally DON'T create a library file
    create_test_samples_without_library(workflow_setup["input_dir"])
    
    # Run the workflow
    result = run_snakemake(
        workflow_dir=workflow_setup["workflow_dir"],
        configfile=workflow_setup["config_file"],
        cores=1
    )
    
    # The workflow should fail or produce expected error output
    # We're looking for a clear error message, not just any failure
    assert not result["success"], "Workflow should fail when library file is missing"
    assert "library file" in result["error"].lower() or "library file" in result["output"].lower(), \
        "Error message should mention the missing library file"


# Test error handling for empty files
def test_empty_fastq_file(setup_test_workflow, run_snakemake):
    """
    Test that the workflow correctly handles empty FASTQ files.
    
    Args:
        setup_test_workflow: Fixture that sets up the test workflow
        run_snakemake: Fixture that provides a function to run Snakemake
    """
    # Get the test workflow setup
    workflow_setup = setup_test_workflow
    
    # Create a valid library file
    create_library_file(workflow_setup["input_dir"])
    
    # Create an empty FASTQ file
    fastq_path = Path(workflow_setup["input_dir"]) / "empty_sample.fastq"
    with open(fastq_path, 'w') as f:
        pass  # Create an empty file
    
    # Run the workflow
    result = run_snakemake(
        workflow_dir=workflow_setup["workflow_dir"],
        configfile=workflow_setup["config_file"],
        cores=1
    )
    
    # The workflow should either:
    # 1. Fail with a clear error message about the empty file
    # 2. Process the file but produce a log message about empty data
    # Either behavior is acceptable as long as it's handled gracefully
    error_text = result["error"].lower() + result["output"].lower()
    assert ("empty" in error_text or "no reads" in error_text or "zero reads" in error_text), \
        "Workflow should report issues with empty FASTQ file"


# Test error handling for malformed input
def test_malformed_fastq_file(setup_test_workflow, run_snakemake):
    """
    Test that the workflow correctly handles malformed FASTQ files.
    
    Args:
        setup_test_workflow: Fixture that sets up the test workflow
        run_snakemake: Fixture that provides a function to run Snakemake
    """
    # Get the test workflow setup
    workflow_setup = setup_test_workflow
    
    # Create a valid library file
    create_library_file(workflow_setup["input_dir"])
    
    # Create a malformed FASTQ file (not following FASTQ format)
    fastq_path = Path(workflow_setup["input_dir"]) / "malformed_sample.fastq"
    with open(fastq_path, 'w') as f:
        f.write("This is not a valid FASTQ format\n")
        f.write("It's missing the @ header, sequence, + line, and quality scores\n")
    
    # Run the workflow
    result = run_snakemake(
        workflow_dir=workflow_setup["workflow_dir"],
        configfile=workflow_setup["config_file"],
        cores=1
    )
    
    # The workflow should report issues with the malformed file
    error_text = result["error"].lower() + result["output"].lower()
    assert ("format" in error_text or "invalid" in error_text or "error" in error_text), \
        "Workflow should report issues with malformed FASTQ file"


# Test error handling for file permission issues
def test_read_only_output_dir(setup_test_workflow, run_snakemake):
    """
    Test that the workflow correctly handles read-only output directories.
    
    Args:
        setup_test_workflow: Fixture that sets up the test workflow
        run_snakemake: Fixture that provides a function to run Snakemake
    """
    # This test may not work on all platforms due to permission handling differences
    # Skip on Windows as chmod doesn't work the same way
    if os.name == 'nt':
        pytest.skip("Read-only permission test not supported on Windows")
    
    # Get the test workflow setup
    workflow_setup = setup_test_workflow
    
    # Create test samples and library
    create_test_samples(workflow_setup["input_dir"])
    
    # Make the output directory read-only
    try:
        os.chmod(workflow_setup["output_dir"], stat.S_IRUSR | stat.S_IXUSR)
    except:
        pytest.skip("Unable to set directory to read-only, skipping test")
    
    # Run the workflow
    result = run_snakemake(
        workflow_dir=workflow_setup["workflow_dir"],
        configfile=workflow_setup["config_file"],
        cores=1
    )
    
    # Reset permissions for cleanup
    try:
        os.chmod(workflow_setup["output_dir"], stat.S_IRWXU)
    except:
        pass
    
    # The workflow should fail with a permission error
    assert not result["success"], "Workflow should fail with read-only output directory"
    error_text = result["error"].lower() + result["output"].lower()
    assert ("permission" in error_text or "access" in error_text), \
        "Error message should mention permission issues"


def create_test_samples_without_library(input_dir):
    """
    Create test FASTQ files but no library file.
    
    Args:
        input_dir (Path): Path to the input directory
    """
    # Create test FASTQ files
    sample_names = ["sample1", "sample2"]
    for sample_name in sample_names:
        fastq_path = Path(input_dir) / f"{sample_name}.fastq"
        with open(fastq_path, 'w') as f:
            f.write("@SEQ_ID_1\n")
            f.write("ACGTACGTACGT\n")
            f.write("+\n")
            f.write("!!!!!!!!!!!!4\n")


def create_library_file(input_dir):
    """
    Create a valid library file for testing.
    
    Args:
        input_dir (Path): Path to the input directory
    """
    library_path = Path(input_dir) / "test_library.txt"
    with open(library_path, 'w') as f:
        f.write("sgRNA\tGene\tSequence\n")
        f.write("sgRNA1\tGeneA\tACGTACGTACGT\n")
        f.write("sgRNA2\tGeneB\tGCTAGCTAGCTA\n")


def create_test_samples(input_dir):
    """
    Create test FASTQ files and library file.
    
    Args:
        input_dir (Path): Path to the input directory
    """
    # Create test FASTQ files
    sample_names = ["sample1", "sample2"]
    for sample_name in sample_names:
        fastq_path = Path(input_dir) / f"{sample_name}.fastq"
        with open(fastq_path, 'w') as f:
            f.write("@SEQ_ID_1\n")
            f.write("ACGTACGTACGT\n")
            f.write("+\n")
            f.write("!!!!!!!!!!!!4\n")
    
    # Create a library file
    create_library_file(input_dir)


if __name__ == "__main__":
    # Allow running this test directly
    pytest.main(["-xvs", __file__]) 