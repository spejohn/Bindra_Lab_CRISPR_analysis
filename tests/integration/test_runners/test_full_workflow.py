#!/usr/bin/env python3
"""
Integration test for the complete CRISPR analysis workflow.

This test verifies that the full workflow runs successfully with the test data
and produces the expected output files.
"""

import os
import pytest
import pandas as pd
from pathlib import Path

# Test the full workflow
def test_full_workflow_execution(setup_test_workflow, run_snakemake, verify_output_files):
    """
    Test that the full CRISPR analysis workflow runs successfully.
    
    This test:
    1. Sets up a test workflow directory with test data
    2. Runs the full Snakemake workflow
    3. Verifies that all expected output files are produced
    
    Args:
        setup_test_workflow: Fixture that sets up the test workflow
        run_snakemake: Fixture that provides a function to run Snakemake
        verify_output_files: Fixture that provides a function to verify outputs
    """
    # Get the test workflow setup
    workflow_setup = setup_test_workflow
    
    # Create minimal test FASTQ and library files
    create_test_data(workflow_setup["input_dir"])
    
    # Run the full workflow
    result = run_snakemake(
        workflow_dir=workflow_setup["workflow_dir"],
        configfile=workflow_setup["config_file"],
        cores=1
    )
    
    # Check that the workflow ran successfully
    assert result["success"], f"Workflow failed with error: {result['error']}"
    
    # Define expected output files
    expected_files = [
        {
            "path": "counts/sample1.count.txt",
            "min_size": 10,
            "contains": "sgRNA\tGene\tCount"
        },
        {
            "path": "counts/sample2.count.txt",
            "min_size": 10,
            "contains": "sgRNA\tGene\tCount"
        },
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
    ]
    
    # Verify the output files
    success, message = verify_output_files(workflow_setup["output_dir"], expected_files)
    assert success, message
    
    # Additional checks on specific output files
    check_output_file_content(workflow_setup["output_dir"])


def create_test_data(input_dir):
    """
    Create minimal test data for the workflow.
    
    Args:
        input_dir (Path): Path to the input directory
    """
    # Create test FASTQ files (minimal content)
    for sample_name in ["sample1", "sample2"]:
        fastq_path = Path(input_dir) / f"{sample_name}.fastq"
        with open(fastq_path, 'w') as f:
            # Write minimal FASTQ content
            f.write("@SEQ_ID_1\n")
            f.write("ACGTACGTACGT\n")
            f.write("+\n")
            f.write("!!!!!!!!!!!!4\n")
            
            f.write("@SEQ_ID_2\n")
            f.write("GCTAGCTAGCTA\n")
            f.write("+\n")
            f.write("!!!!!!!!!!!!4\n")
    
    # Create a test library file
    library_path = Path(input_dir) / "test_library.txt"
    with open(library_path, 'w') as f:
        f.write("sgRNA\tGene\tSequence\n")
        f.write("sgRNA1\tGeneA\tACGTACGTACGT\n")
        f.write("sgRNA2\tGeneB\tGCTAGCTAGCTA\n")
        
    # Create a test contrasts file
    contrasts_path = Path(input_dir) / "test_contrasts.txt"
    with open(contrasts_path, 'w') as f:
        f.write("contrast,control,treatment\n")
        f.write("test_contrast,sample1,sample2\n")


def check_output_file_content(output_dir):
    """
    Perform additional checks on specific output files.
    
    Args:
        output_dir (Path): Path to the output directory
    """
    # Example: Check the gene summary file has the expected structure
    gene_summary_file = Path(output_dir) / "analysis/mageck_test/test_contrast.gene_summary.txt"
    if gene_summary_file.exists():
        try:
            df = pd.read_csv(gene_summary_file, sep='\t')
            # Check that required columns are present
            required_columns = ["Gene", "p-value", "FDR", "LFC"]
            for col in required_columns:
                assert col in df.columns, f"Column {col} not found in gene summary file"
                
            # Check that there is at least one row of data
            assert len(df) > 0, "Gene summary file is empty"
        except Exception as e:
            pytest.fail(f"Error checking gene summary file: {str(e)}")
            
    # Could add more checks for other output files as needed


if __name__ == "__main__":
    # Allow running this test directly
    pytest.main(["-xvs", __file__]) 