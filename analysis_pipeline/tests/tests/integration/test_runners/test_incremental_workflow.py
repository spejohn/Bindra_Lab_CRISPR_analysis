#!/usr/bin/env python3
"""
Integration test for the incremental workflow functionality.

This test verifies that Snakemake correctly processes new samples added to an existing analysis directory.
"""

import os
import time
import pytest
import pandas as pd
from pathlib import Path

# Test incremental workflow (adding new samples)
def test_incremental_workflow(setup_test_workflow, run_snakemake, verify_output_files):
    """
    Test that Snakemake correctly processes only new samples.
    
    This test:
    1. Sets up a test workflow with initial samples
    2. Runs the workflow to process these samples
    3. Adds a new sample
    4. Runs the workflow again
    5. Verifies that only the new sample is processed
    
    Args:
        setup_test_workflow: Fixture that sets up the test workflow
        run_snakemake: Fixture that provides a function to run Snakemake
        verify_output_files: Fixture that provides a function to verify outputs
    """
    # Get the test workflow setup
    workflow_setup = setup_test_workflow
    
    # Create initial test data (2 samples)
    initial_samples = ["sample1", "sample2"]
    create_test_samples(workflow_setup["input_dir"], initial_samples)
    
    # Run the workflow for initial samples
    result1 = run_snakemake(
        workflow_dir=workflow_setup["workflow_dir"],
        configfile=workflow_setup["config_file"],
        cores=1
    )
    
    # Check that the initial workflow ran successfully
    assert result1["success"], f"Initial workflow failed with error: {result1['error']}"
    
    # Wait a moment to ensure file timestamps will be different
    time.sleep(1)
    
    # Capture file modification times after the first run
    initial_files = {}
    for sample in initial_samples:
        count_file = Path(workflow_setup["output_dir"]) / "counts" / f"{sample}.count.txt"
        if count_file.exists():
            initial_files[sample] = count_file.stat().st_mtime
    
    # Add a new sample
    new_sample = "sample3"
    create_test_samples(workflow_setup["input_dir"], [new_sample])
    
    # Update the contrasts file to include the new sample
    update_contrasts_file(workflow_setup["input_dir"], initial_samples, [new_sample])
    
    # Run the workflow again
    result2 = run_snakemake(
        workflow_dir=workflow_setup["workflow_dir"],
        configfile=workflow_setup["config_file"],
        cores=1
    )
    
    # Check that the second workflow ran successfully
    assert result2["success"], f"Incremental workflow failed with error: {result2['error']}"
    
    # Verify that all output files exist
    all_samples = initial_samples + [new_sample]
    expected_files = []
    for sample in all_samples:
        expected_files.append({
            "path": f"counts/{sample}.count.txt",
            "min_size": 10,
            "contains": "sgRNA\tGene\tCount"
        })
    
    expected_files.extend([
        {
            "path": "merged/count_table.txt",
            "min_size": 10,
            "contains": "sgRNA\tGene"
        },
        {
            "path": "analysis/mageck_test/test_contrast.gene_summary.txt",
            "min_size": 10,
            "contains": "Gene\tsgRNA\tp-value"
        },
        {
            "path": "qc/qc_report.html",
            "min_size": 100,
            "contains": "CRISPR Analysis QC Report"
        }
    ])
    
    success, message = verify_output_files(workflow_setup["output_dir"], expected_files)
    assert success, message
    
    # Check that initial files weren't reprocessed (unchanged timestamps)
    for sample in initial_samples:
        count_file = Path(workflow_setup["output_dir"]) / "counts" / f"{sample}.count.txt"
        assert count_file.exists(), f"Count file for {sample} is missing"
        
        # Check if timestamp matches the initial one
        # We only assert that the individual sample count files weren't regenerated
        # Merged files will naturally be regenerated
        current_mtime = count_file.stat().st_mtime
        assert current_mtime == initial_files[sample], \
            f"File for {sample} was regenerated (mtime changed from {initial_files[sample]} to {current_mtime})"
    
    # Check that the new sample was processed
    new_count_file = Path(workflow_setup["output_dir"]) / "counts" / f"{new_sample}.count.txt"
    assert new_count_file.exists(), f"Count file for new sample {new_sample} is missing"


def create_test_samples(input_dir, sample_names):
    """
    Create test FASTQ files for the specified samples.
    
    Args:
        input_dir (Path): Path to the input directory
        sample_names (list): List of sample names to create
    """
    # Create test FASTQ files (minimal content)
    for sample_name in sample_names:
        fastq_path = Path(input_dir) / f"{sample_name}.fastq"
        with open(fastq_path, 'w') as f:
            # Write minimal FASTQ content
            f.write("@SEQ_ID_1\n")
            f.write("ACGTACGTACGT\n")
            f.write("+\n")
            f.write("!!!!!!!!!!!!4\n")
            
            # Add a second read
            f.write("@SEQ_ID_2\n")
            f.write("GCTAGCTAGCTA\n")
            f.write("+\n")
            f.write("!!!!!!!!!!!!4\n")
    
    # Create library file if it doesn't exist
    library_path = Path(input_dir) / "test_library.txt"
    if not library_path.exists():
        with open(library_path, 'w') as f:
            f.write("sgRNA\tGene\tSequence\n")
            f.write("sgRNA1\tGeneA\tACGTACGTACGT\n")
            f.write("sgRNA2\tGeneB\tGCTAGCTAGCTA\n")


def update_contrasts_file(input_dir, control_samples, treatment_samples):
    """
    Update the contrasts file to include all samples.
    
    Args:
        input_dir (Path): Path to the input directory
        control_samples (list): List of control sample names
        treatment_samples (list): List of treatment sample names
    """
    contrasts_path = Path(input_dir) / "test_contrasts.txt"
    
    with open(contrasts_path, 'w') as f:
        f.write("contrast,control,treatment\n")
        control_str = "|".join(control_samples)
        treatment_str = "|".join(treatment_samples)
        f.write(f"test_contrast,{control_str},{treatment_str}\n")


if __name__ == "__main__":
    # Allow running this test directly
    pytest.main(["-xvs", __file__]) 