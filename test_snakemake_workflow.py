#!/usr/bin/env python3
"""
Script to run a test of the Snakemake workflow with test data.
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path

# Set up paths
project_root = Path(__file__).parent.absolute()
test_data_dir = project_root / "tests/test_data"
snakefile = project_root / "Snakefile"
output_dir = project_root / "test_output"

# Create test directory structure
input_dir = output_dir / "input"
library_dir = input_dir / "library"
fastq_dir = input_dir / "fastq"

# Clean up previous test outputs
if output_dir.exists():
    print(f"Cleaning up previous test output: {output_dir}")
    # Just clean subdirectories, leaving the directory itself
    for item in output_dir.iterdir():
        if item.is_dir():
            shutil.rmtree(item)
        else:
            item.unlink()
            
# Create directory structure
output_dir.mkdir(exist_ok=True)
input_dir.mkdir(exist_ok=True)
library_dir.mkdir(exist_ok=True)
fastq_dir.mkdir(exist_ok=True)

# Create a minimal test config
config_path = output_dir / "test_config.yaml"
with open(config_path, 'w') as f:
    f.write("""# Test config for Snakemake workflow
input_dir: "input"
output_dir: "results"
experiment_name: "test_run"
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

print(f"Test data prepared in: {output_dir}")
print(f"- Library file: {library_file}")
print(f"- Contrasts file: {contrasts_file}")
print(f"- Config file: {config_path}")

# Run Snakemake in dry-run mode
print("\n=== Running Snakemake in dry-run mode ===\n")
try:
    snakemake_cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--directory", str(output_dir),
        "--configfile", str(config_path),
        "--cores", "1",
        "--dryrun",
        "-p"  # print shell commands that would be executed
    ]
    
    result = subprocess.run(
        snakemake_cmd,
        check=False,
        capture_output=True,
        text=True
    )
    
    print(f"\nDry run exit code: {result.returncode}")
    
    if result.stdout:
        print("\nOutput:")
        print(result.stdout)
    
    if result.stderr:
        print("\nErrors:")
        print(result.stderr)
        
    if result.returncode == 0:
        print("\nDry run successful!")
        print("\nTo run the actual workflow, execute:")
        print(f"cd {output_dir} && snakemake --snakefile {snakefile} --configfile {config_path} --cores 1")
    else:
        print("\nDry run failed. Check errors above.")
        sys.exit(1)
        
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1) 