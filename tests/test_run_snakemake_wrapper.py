# tests/test_run_snakemake_wrapper.py
import pytest
import subprocess
from unittest.mock import patch, MagicMock
import sys
import os
from pathlib import Path

# Add project root to sys.path to import run_snakemake
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import the functions to test AFTER modifying sys.path
import run_snakemake

# Fixture to create dummy directories for testing
@pytest.fixture
def mock_dirs(tmp_path):
    base_dir = tmp_path / "mock_base"
    profile_dir = tmp_path / "mock_profile"
    base_dir.mkdir()
    profile_dir.mkdir()
    # Create dummy Snakefile in parent of tmp_path (simulating script location)
    (tmp_path.parent / "Snakefile").touch() 
    yield {"base": str(base_dir), "profile": str(profile_dir)}
    (tmp_path.parent / "Snakefile").unlink()


# --- Tests for argument parsing ---
def test_parse_args_defaults(mock_dirs):
    """Test default argument parsing."""
    with patch.object(sys, 'argv', ['run_snakemake.py', mock_dirs['base']]):
        args = run_snakemake.parse_args()
        assert args.base_dir == mock_dirs['base']
        assert args.output_dir is None
        assert args.target_experiments is None
        assert args.profile is None
        assert args.cores is None # Will be auto-detected later
        assert not args.dryrun
        assert not args.skip_qc
        assert not args.skip_mle
        assert not args.skip_drugz

def test_parse_args_custom(mock_dirs):
    """Test custom arguments."""
    test_args = [
        'run_snakemake.py', mock_dirs['base'],
        '-o', 'custom_out',
        '--target_experiments', 'exp1', 'exp2',
        '-j', '16',
        '--profile', mock_dirs['profile'],
        '--skip-qc',
        '--dryrun'
    ]
    with patch.object(sys, 'argv', test_args):
        args = run_snakemake.parse_args()
        assert args.base_dir == mock_dirs['base']
        assert args.output_dir == 'custom_out'
        assert args.target_experiments == ['exp1', 'exp2']
        assert args.cores == 16
        assert args.profile == mock_dirs['profile']
        assert args.skip_qc
        assert args.dryrun
        assert not args.skip_mle

# --- Tests for snakemake command construction ---
@patch('run_snakemake.subprocess.Popen')
def test_run_snakemake_command_basic(mock_popen, mock_dirs):
    """Test the basic snakemake command structure."""
    # Mock Popen and communicate
    mock_proc = MagicMock()
    mock_proc.communicate.return_value = ('stdout', 'stderr')
    mock_proc.returncode = 0
    mock_popen.return_value = mock_proc

    # Setup args
    args = MagicMock()
    args.base_dir = mock_dirs['base']
    args.output_dir = None # Use default logic
    args.target_experiments = None
    args.cores = 8
    args.profile = None
    args.dryrun = False
    args.skip_qc = False
    args.skip_mle = False
    args.skip_drugz = False
    
    # Expected default output dir calculation
    expected_output_dir = Path(mock_dirs['base']).resolve().parent / "crispr_analysis_pipeline_results"

    run_snakemake.run_snakemake(args)

    mock_popen.assert_called_once()
    call_args, call_kwargs = mock_popen.call_args
    cmd_str = call_args[0] # The command string

    # Check essential parts (use 'in' because exact quoting/spacing might vary slightly)
    assert 'snakemake' in cmd_str
    assert f'-s {Path(__file__).parent.parent / "Snakefile"}' in cmd_str # Check Snakefile path
    assert '-j 8' in cmd_str
    assert f'base_dir={mock_dirs["base"]}' in cmd_str
    assert f'output_dir={expected_output_dir}' in cmd_str
    assert 'skip_qc=false' in cmd_str
    assert 'use_apptainer=true' in cmd_str
    assert 'target_screens=None' in cmd_str # Check None is passed
    assert '--profile' not in cmd_str
    assert '--dryrun' not in cmd_str

@patch('run_snakemake.subprocess.Popen')
def test_run_snakemake_command_targets_profile_dryrun(mock_popen, mock_dirs):
    """Test command with targets, profile, and dryrun."""
    # Mock Popen
    mock_proc = MagicMock()
    mock_proc.communicate.return_value = ('', '')
    mock_proc.returncode = 0
    mock_popen.return_value = mock_proc
    
    # Setup args
    args = MagicMock()
    args.base_dir = mock_dirs['base']
    args.output_dir = 'test_out'
    args.target_experiments = ['exp_A', 'exp_B']
    args.cores = None # Test auto-detect path (mock get_available_cores if needed)
    args.profile = mock_dirs['profile']
    args.dryrun = True
    args.skip_qc = True
    args.skip_mle = False
    args.skip_drugz = False

    run_snakemake.run_snakemake(args)

    mock_popen.assert_called_once()
    call_args, call_kwargs = mock_popen.call_args
    cmd_str = call_args[0]

    assert 'snakemake' in cmd_str
    assert f'base_dir={mock_dirs["base"]}' in cmd_str
    assert 'output_dir=test_out' in cmd_str
    assert 'skip_qc=true' in cmd_str
    # Check list formatting for target_screens
    assert 'target_screens=[\'exp_A\',\'exp_B\']' in cmd_str.replace(" ", "") # Remove spaces for robust check
    assert f'--profile {mock_dirs["profile"]}' in cmd_str
    assert '--dryrun' in cmd_str

# --- Test error handling ---
@patch('run_snakemake.os.path.isdir')
def test_main_missing_basedir(mock_isdir):
    """Test main function exits if base_dir doesn't exist."""
    mock_isdir.return_value = False
    with patch.object(sys, 'argv', ['run_snakemake.py', 'nonexistent_dir']):
        with pytest.raises(SystemExit) as e:
            run_snakemake.main()
        assert e.value.code == 1 # Check exit code

@patch('run_snakemake.os.path.isdir')
def test_main_missing_profiledir(mock_isdir, mock_dirs):
    """Test main function exits if profile dir doesn't exist."""
    # Make base_dir exist, but profile dir not exist
    mock_isdir.side_effect = lambda x: x == mock_dirs['base'] 
    with patch.object(sys, 'argv', ['run_snakemake.py', mock_dirs['base'], '--profile', 'nonexistent_profile']):
        with pytest.raises(SystemExit) as e:
            run_snakemake.main()
        assert e.value.code == 1
