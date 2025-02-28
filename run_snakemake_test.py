#!/usr/bin/env python3
"""
Script to test Snakemake functionality with the CRISPR analysis pipeline.
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

# Ensure analysis_pipeline is importable 
project_root = Path(__file__).parent.absolute()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Try to import analysis_pipeline
try:
    import analysis_pipeline
    print(f"Using analysis_pipeline version {analysis_pipeline.__version__}")
except ImportError:
    print("Warning: Could not import analysis_pipeline module.")
    print("You may need to run `python install_dev.py` first.")
    print("Continuing anyway, but the test may fail.")

# Set up paths
config_file = project_root / "test_output/test_config.yaml"
workflow_dir = project_root

# Create test output directory
test_output_dir = project_root / "test_output"
test_output_dir.mkdir(exist_ok=True)

# Clean any previous test data
for item in test_output_dir.iterdir():
    if item.is_dir():
        shutil.rmtree(item)
    else:
        os.unlink(item)

# Create test directory structure
input_dir = test_output_dir / "input"
library_dir = input_dir / "library"
fastq_dir = input_dir / "fastq"
results_dir = test_output_dir / "results"

input_dir.mkdir(exist_ok=True)
library_dir.mkdir(exist_ok=True)
fastq_dir.mkdir(exist_ok=True)
results_dir.mkdir(exist_ok=True)

# Create a minimal test config
with open(config_file, 'w') as f:
    f.write("""# Test config for Snakemake workflow
input_dir: "test_output/input"
output_dir: "test_output/results"
experiment_name: "test_run"
library_file: "test_output/input/library/test_library.txt"
contrasts_file: "test_output/input/test_contrasts.txt"
norm_method: "median"
fdr_threshold: 0.05
threads: 1
memory_gb: 2
is_test: true
skip_long_steps: true
""")

# Create minimal test library and data files
library_file = library_dir / "test_library.txt"
with open(library_file, 'w') as f:
    f.write("sgRNA\tGene\tSequence\n")
    f.write("sgRNA1\tGeneA\tACGTACGTACGT\n")
    f.write("sgRNA2\tGeneB\tGCTAGCTAGCTA\n")

# Create test fastq files
for sample in ["sample1", "sample2"]:
    fastq_file = fastq_dir / f"{sample}.fastq"
    with open(fastq_file, 'w') as f:
        f.write("@SEQ_ID_1\n")
        f.write("ACGTACGTACGT\n")
        f.write("+\n")
        f.write("IIIIIIIIIIII\n")
        f.write("@SEQ_ID_2\n")
        f.write("GCTAGCTAGCTA\n")
        f.write("+\n")
        f.write("IIIIIIIIIIII\n")

# Create a test contrasts file
contrasts_file = input_dir / "test_contrasts.txt"
with open(contrasts_file, 'w') as f:
    f.write("contrast,control,treatment\n")
    f.write("test_contrast,sample1,sample2\n")

print(f"Testing Snakemake workflow in directory: {workflow_dir}")
print(f"Using config file: {config_file}")
print(f"Output will be in: {test_output_dir}")
print(f"Test input data created in: {input_dir}")

# Run Snakemake in dry-run mode first
try:
    print("\n=== Running Snakemake in dry-run mode ===\n")
    dryrun_cmd = [
        "snakemake",
        "--cores", "1",
        "--configfile", str(config_file),
        "--dryrun",
        "-p"  # print shell commands that will be executed
    ]
    
    dryrun_result = subprocess.run(
        dryrun_cmd,
        check=False,
        cwd=project_root,
        capture_output=True,
        text=True
    )
    
    print(f"Dry run exit code: {dryrun_result.returncode}")
    print("Output:")
    print(dryrun_result.stdout)
    
    if dryrun_result.stderr:
        print("Errors:")
        print(dryrun_result.stderr)
        
    # If dry run successful, offer to run for real
    if dryrun_result.returncode == 0:
        print("\nDry run successful! Ready to run the actual workflow.")
        
        run_actual = input("\nDo you want to run the actual workflow? (y/n): ").strip().lower()
        if run_actual == 'y':
            print("\n=== Running actual Snakemake workflow ===\n")
            run_cmd = [
                "snakemake",
                "--cores", "1", 
                "--configfile", str(config_file)
            ]
            
            run_result = subprocess.run(
                run_cmd, 
                check=False,
                cwd=project_root
            )
            
            print(f"\nWorkflow run completed with exit code: {run_result.returncode}")
        else:
            print("\nSkipping actual workflow run.")
            print(f"To run manually, execute:")
            print(f"  snakemake --cores 1 --configfile {config_file}")
    else:
        print("\nDry run failed. Fix the errors before running the actual workflow.")
        
except Exception as e:
    print(f"Error running Snakemake: {e}")
    sys.exit(1) 