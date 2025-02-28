#!/usr/bin/env python3
"""
Debugging script for integration test workflow issues.

This script sets up a test workflow similar to the tests and runs it
with detailed logging to diagnose issues.
"""

import os
import sys
import shutil
import tempfile
import subprocess
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def setup_test_workflow():
    """
    Set up a test workflow in a temporary directory.
    
    Returns:
        dict: Dictionary with workflow paths
    """
    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp(prefix="crispr_test_debug_"))
    logger.info(f"Created temporary directory: {temp_dir}")
    
    # Create input and output directories
    input_dir = temp_dir / "input"
    output_dir = temp_dir / "results"
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Get the source directories
    script_dir = Path(__file__).parent
    test_workflow_dir = script_dir / "test_workflow"
    
    # Copy the test Snakefile
    shutil.copy(
        test_workflow_dir / "Snakefile",
        temp_dir / "Snakefile"
    )
    logger.info(f"Copied Snakefile to {temp_dir}")
    
    # Copy the test config
    shutil.copy(
        test_workflow_dir / "config.yaml",
        temp_dir / "config.yaml"
    )
    logger.info(f"Copied config.yaml to {temp_dir}")
    
    # Create test data
    create_test_data(input_dir)
    logger.info(f"Created test data in {input_dir}")
    
    return {
        "workflow_dir": temp_dir,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "config_file": temp_dir / "config.yaml"
    }

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

def run_snakemake(workflow_dir, configfile=None, cores=1, dryrun=False):
    """
    Run Snakemake with detailed logging.
    
    Args:
        workflow_dir: Directory containing the Snakefile
        configfile: Path to config file
        cores: Number of cores to use
        dryrun: Whether to do a dry run
        
    Returns:
        dict: Results of the Snakemake run
    """
    workflow_dir = Path(workflow_dir)
    
    # Build the Snakemake command
    cmd = ["snakemake", "--cores", str(cores)]
    
    if configfile:
        cmd.extend(["--configfile", str(configfile)])
        
    if dryrun:
        cmd.append("--dryrun")
    
    # Add verbose flag for debugging
    cmd.append("--verbose")
    
    # Log the command
    logger.info(f"Running Snakemake with command: {' '.join(map(str, cmd))}")
    
    # Run Snakemake
    try:
        result = subprocess.run(
            cmd,
            cwd=str(workflow_dir),
            capture_output=True,
            text=True,
            check=False
        )
        
        success = result.returncode == 0
        
        # Log the results
        if success:
            logger.info("Snakemake run successful")
        else:
            logger.error(f"Snakemake run failed with return code {result.returncode}")
            
        if result.stdout:
            logger.debug(f"Snakemake stdout:\n{result.stdout}")
        if result.stderr:
            logger.error(f"Snakemake stderr:\n{result.stderr}")
        
        return {
            "success": success,
            "output": result.stdout,
            "error": result.stderr,
            "returncode": result.returncode
        }
    except Exception as e:
        logger.exception(f"Exception running Snakemake: {e}")
        return {
            "success": False,
            "output": "",
            "error": str(e),
            "returncode": -1
        }

def verify_outputs(workflow_dir):
    """
    Check for expected output files and print their status.
    
    Args:
        workflow_dir: Workflow directory containing the results
    """
    output_dir = Path(workflow_dir) / "results"
    
    expected_files = [
        "counts/sample1.count.txt",
        "counts/sample2.count.txt",
        "merged/count_table.txt",
        "analysis/mageck_test/test_contrast.gene_summary.txt",
        "qc/qc_report.html"
    ]
    
    logger.info("Verifying output files:")
    for file_path in expected_files:
        full_path = output_dir / file_path
        if full_path.exists():
            size = full_path.stat().st_size
            logger.info(f"✓ {file_path} exists (size: {size} bytes)")
        else:
            logger.error(f"✗ {file_path} does not exist")
            
            # Check if parent directory exists
            parent_dir = full_path.parent
            if parent_dir.exists():
                logger.info(f"  Parent directory {parent_dir} exists")
                # List files in parent directory
                files = list(parent_dir.iterdir())
                if files:
                    logger.info(f"  Files in {parent_dir}:")
                    for f in files:
                        logger.info(f"    - {f.name}")
                else:
                    logger.info(f"  No files in {parent_dir}")
            else:
                logger.error(f"  Parent directory {parent_dir} does not exist")

def main():
    """Run the debugging workflow"""
    logger.info("Starting debug workflow")
    
    # Setup the test workflow
    workflow = setup_test_workflow()
    
    # First try a dry run
    logger.info("Running Snakemake in dry-run mode first")
    dryrun_result = run_snakemake(
        workflow_dir=workflow["workflow_dir"], 
        configfile=workflow["config_file"], 
        cores=1, 
        dryrun=True
    )
    
    # Then run for real
    logger.info("Running Snakemake for real")
    result = run_snakemake(
        workflow_dir=workflow["workflow_dir"], 
        configfile=workflow["config_file"], 
        cores=1
    )
    
    # Verify outputs
    verify_outputs(workflow["workflow_dir"])
    
    # Report results
    if result["success"]:
        logger.info("Debug workflow completed successfully")
    else:
        logger.error("Debug workflow failed")
    
    logger.info(f"Temporary directory: {workflow['workflow_dir']}")
    logger.info("You can inspect this directory for debugging")

if __name__ == "__main__":
    main() 