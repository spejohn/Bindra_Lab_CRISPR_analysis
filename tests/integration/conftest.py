"""
Fixtures for integration testing with Snakemake.

This module provides pytest fixtures specifically for integration testing
the CRISPR analysis pipeline with Snakemake.
"""

import os
import shutil
import tempfile
import subprocess
import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def integration_test_data_dir():
    """
    Fixture providing the path to the integration test data directory.
    
    Returns:
        Path: Path to the test data directory for integration tests
    """
    test_dir = Path(__file__).parent / "test_data"
    os.makedirs(test_dir, exist_ok=True)
    return test_dir


@pytest.fixture(scope="function")
def temp_workflow_dir():
    """
    Creates a temporary directory for running a workflow test.
    
    This fixture creates a clean temporary directory for each test,
    which is automatically cleaned up after the test completes.
    
    Returns:
        Path: Path to a temporary directory for workflow testing
    """
    temp_dir = Path(tempfile.mkdtemp(prefix="crispr_test_"))
    yield temp_dir
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture(scope="function")
def setup_test_workflow(integration_test_data_dir, temp_workflow_dir):
    """
    Sets up a test workflow by copying test data and configuration.
    
    Args:
        integration_test_data_dir: Fixture providing the test data directory
        temp_workflow_dir: Fixture providing a temporary directory
        
    Returns:
        dict: A dictionary containing paths and configuration for the test
    """
    # Create input and output directories
    input_dir = temp_workflow_dir / "input"
    os.makedirs(input_dir, exist_ok=True)
    
    # The output directory should match what's specified in config.yaml
    output_dir = temp_workflow_dir / "crispr_analysis_pipeline_results"
    
    # Copy test data and workflow files
    test_workflow_dir = Path(__file__).parent / "test_workflow"
    
    # Copy the test Snakefile
    shutil.copy(
        test_workflow_dir / "Snakefile",
        temp_workflow_dir / "Snakefile"
    )
    
    # Copy the test config
    shutil.copy(
        test_workflow_dir / "config.yaml",
        temp_workflow_dir / "config.yaml"
    )
    
    # Copy test data
    test_data_src = integration_test_data_dir
    for item in test_data_src.glob("*"):
        if item.is_dir():
            shutil.copytree(item, input_dir / item.name)
        else:
            shutil.copy(item, input_dir)
    
    return {
        "workflow_dir": temp_workflow_dir,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "config_file": temp_workflow_dir / "config.yaml"
    }


@pytest.fixture(scope="function")
def run_snakemake():
    """
    Fixture that provides a function to run Snakemake with specified parameters.
    
    Returns:
        callable: A function that runs Snakemake with the provided arguments
    """
    def _run_snakemake(workflow_dir, configfile=None, cores=1, dryrun=False, targets=None):
        """
        Run Snakemake with the specified parameters.
        
        Args:
            workflow_dir (str or Path): Directory containing the Snakefile
            configfile (str or Path, optional): Path to config file
            cores (int, optional): Number of cores to use. Defaults to 1.
            dryrun (bool, optional): Whether to do a dry run. Defaults to False.
            targets (list, optional): List of target rules to run. Defaults to None.
            
        Returns:
            dict: Dictionary with 'success' (bool) and 'output' (str) keys
        """
        workflow_dir = Path(workflow_dir)
        
        # Build the Snakemake command
        cmd = ["snakemake", "--cores", str(cores)]
        
        if configfile:
            cmd.extend(["--configfile", str(configfile)])
            
        if dryrun:
            cmd.append("--dryrun")
            
        if targets:
            cmd.extend(targets)
        
        # Print the command for debugging
        print(f"Running command: {' '.join(cmd)} in directory {workflow_dir}")

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
            
            # Print error info if the command failed
            if not success:
                print(f"Command failed with return code: {result.returncode}")
                print(f"Error output: {result.stderr}")
                
            return {
                "success": success,
                "output": result.stdout,
                "error": result.stderr,
                "returncode": result.returncode
            }
        except Exception as e:
            return {
                "success": False,
                "output": "",
                "error": str(e),
                "returncode": -1
            }
    
    return _run_snakemake


@pytest.fixture(scope="function")
def verify_output_files():
    """
    Fixture that provides a function to verify the existence and properties of output files.
    
    Returns:
        callable: A function that verifies output files match expected properties
    """
    def _verify_output_files(output_dir, expected_files):
        """
        Verify that output files exist and match expected properties.
        
        Args:
            output_dir (str or Path): Directory containing output files
            expected_files (list): List of dictionaries with file specifications
                Each dict should have at least a 'path' key and optionally:
                'min_size': minimum file size in bytes
                'contains': text that should be in the file
        
        Returns:
            tuple: (bool success, str message)
        """
        output_dir = Path(output_dir)
        missing_files = []
        invalid_files = []
        
        for file_spec in expected_files:
            file_path = output_dir / file_spec['path']
            
            # Check if file exists
            if not file_path.exists():
                missing_files.append(file_spec['path'])
                continue
                
            # Check file size if specified
            if 'min_size' in file_spec and file_path.stat().st_size < file_spec['min_size']:
                invalid_files.append(f"{file_spec['path']} (too small)")
                continue
                
            # Check file contents if specified
            if 'contains' in file_spec:
                with open(file_path, 'r') as f:
                    content = f.read()
                    if file_spec['contains'] not in content:
                        invalid_files.append(f"{file_spec['path']} (missing expected content)")
                        continue
        
        success = not (missing_files or invalid_files)
        messages = []
        
        if missing_files:
            messages.append(f"Missing files: {', '.join(missing_files)}")
        
        if invalid_files:
            messages.append(f"Invalid files: {', '.join(invalid_files)}")
            
        return success, "\n".join(messages) if messages else "All output files verified successfully"
    
    return _verify_output_files 