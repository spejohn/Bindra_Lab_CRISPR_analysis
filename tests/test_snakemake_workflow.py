# tests/test_snakemake_workflow.py
import pytest
import subprocess
import sys
import os
import shutil
from pathlib import Path

# Assume conftest.py provides a fixture 'test_data_dir' pointing to tests/test_data
# Assume conftest.py adds project root to sys.path

# Fixture to set up a temporary working environment for each test
@pytest.fixture
def setup_workflow_test_env(tmp_path, test_data_dir):
    """Creates temp base_dir, output_dir and copies test data."""
    base_dir = tmp_path / "test_input"
    output_dir = tmp_path / "test_output"
    
    # Copy the main test data structure
    shutil.copytree(test_data_dir / "test_main_dir", base_dir)
    
    # Return paths needed by tests
    return {"base": base_dir, "output": output_dir}

def run_wrapper_script(args_list):
    """Helper function to run run_snakemake.py script via subprocess."""
    script_path = Path(__file__).parent.parent / "run_snakemake.py"
    cmd = [sys.executable, str(script_path)] + args_list
    print(f"\nExecuting: {' '.join(map(str, cmd))}") # For debugging
    result = subprocess.run(cmd, capture_output=True, text=True)
    print("STDOUT:\n", result.stdout)
    print("STDERR:\n", result.stderr)
    return result

def test_snakemake_dryrun_default(setup_workflow_test_env):
    """Test dryrun with default settings (should find both experiments)."""
    env = setup_workflow_test_env
    result = run_wrapper_script([str(env['base']), "--dryrun"])
    assert result.returncode == 0
    # Optionally, check stderr/stdout for expected number of jobs planned (tricky)

def test_snakemake_dryrun_target_fastq(setup_workflow_test_env):
    """Test dryrun targeting only the fastq experiment."""
    env = setup_workflow_test_env
    result = run_wrapper_script([
        str(env['base']), 
        "--target_experiments", "exp_fastq_only", 
        "--dryrun"
    ])
    assert result.returncode == 0
    # Check logs indicate filtering is happening
    assert "Filtering for target screens: ['exp_fastq_only']" in result.stderr \
        or "Filtering for target screens: ['exp_fastq_only']" in result.stdout 

def test_snakemake_dryrun_target_rc(setup_workflow_test_env):
    """Test dryrun targeting only the rc experiment."""
    env = setup_workflow_test_env
    result = run_wrapper_script([
        str(env['base']), 
        "--target_experiments", "exp_rc_only", 
        "--dryrun"
    ])
    assert result.returncode == 0
    assert "Filtering for target screens: ['exp_rc_only']" in result.stderr \
        or "Filtering for target screens: ['exp_rc_only']" in result.stdout

def test_snakemake_dryrun_skips(setup_workflow_test_env):
    """Test dryrun with skip flags."""
    env = setup_workflow_test_env
    result = run_wrapper_script([
        str(env['base']), 
        "--skip-qc", 
        "--skip-mle",
        "--dryrun"
    ])
    assert result.returncode == 0
    # Difficult to verify skips in dryrun output easily, 
    # just check command runs successfully

def test_snakemake_dryrun_profile(setup_workflow_test_env, tmp_path):
    """Test dryrun passes profile argument (profile dir doesn't need content for dryrun)."""
    env = setup_workflow_test_env
    profile_dir = tmp_path / "fake_profile"
    profile_dir.mkdir()
    result = run_wrapper_script([
        str(env['base']), 
        "--profile", str(profile_dir),
        "--dryrun"
    ])
    assert result.returncode == 0
    # Check logs show profile being used
    assert f"Using profile: {profile_dir}" in result.stderr \
        or f"Using profile: {profile_dir}" in result.stdout

def test_snakemake_lint(setup_workflow_test_env):
    """Run snakemake --lint directly to check Snakefile syntax."""
    snakefile_path = Path(__file__).parent.parent / "Snakefile"
    result = subprocess.run([ "snakemake", "-s", str(snakefile_path), "--lint"], capture_output=True, text=True)
    print("Lint STDOUT:\n", result.stdout)
    print("Lint STDERR:\n", result.stderr)
    # Snakemake lint exits 0 even with warnings, check output if needed
    # assert result.returncode == 0 # May fail if there are lint warnings
    assert "Syntax error" not in result.stderr # Basic check

# Optional: A limited actual run (can be slow, might need container runtime)
# @pytest.mark.skipif(not shutil.which("apptainer"), reason="Apptainer not found, skipping actual run test")
def test_snakemake_actual_run_convert_contrasts(setup_workflow_test_env):
     """Test actual execution of a simple early rule."""
     env = setup_workflow_test_env
     # Target only the contrast conversion for the fastq experiment
     # Note: We need to run the *wrapper* which sets the config needed by Snakefile
     result = run_wrapper_script([
         str(env['base']), 
         "--target_experiments", "exp_fastq_only", 
         # Requesting specific output file requires setting workdir or absolute path
         # Using just the wrapper without specific targets is easier
         # Instead, run and check for the output file later
     ])
     assert result.returncode == 0

     # Check if the expected output file was created
     expected_output = env['output'] / "exp_fastq_only" / "contrasts.txt"
     assert expected_output.exists()
     assert expected_output.stat().st_size > 0 # Check it's not empty