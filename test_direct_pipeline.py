#!/usr/bin/env python3
"""
Script to test the analysis_pipeline directly, bypassing Snakemake.
"""

import os
import sys
import shutil
from pathlib import Path

# Ensure the analysis_pipeline is importable
project_root = Path(__file__).parent.absolute()
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Import from the installed package
try:
    import analysis_pipeline
    from analysis_pipeline.pipeline import run_pipeline
    print(f"Using analysis_pipeline version {analysis_pipeline.__version__}")
except ImportError as e:
    print(f"Error importing analysis_pipeline: {e}")
    print("Please run `python install_dev.py` first.")
    sys.exit(1)

# Set up test data paths
test_output_dir = project_root / "test_output"
input_dir = test_output_dir / "input"
library_dir = input_dir / "library"
fastq_dir = input_dir / "fastq"
results_dir = test_output_dir / "results"

# Clean up any previous test results (keep input data)
if results_dir.exists():
    shutil.rmtree(results_dir)
results_dir.mkdir(exist_ok=True)

# Create test directories if they don't exist
input_dir.mkdir(exist_ok=True)
library_dir.mkdir(exist_ok=True)
fastq_dir.mkdir(exist_ok=True)

# Create minimal test library file if it doesn't exist
library_file = library_dir / "test_library.txt"
if not library_file.exists():
    with open(library_file, 'w') as f:
        f.write("sgRNA\tGene\tSequence\n")
        f.write("sgRNA1\tGeneA\tACGTACGTACGT\n")
        f.write("sgRNA2\tGeneB\tGCTAGCTAGCTA\n")

# Create test fastq files if they don't exist
for sample in ["sample1", "sample2"]:
    fastq_file = fastq_dir / f"{sample}.fastq"
    if not fastq_file.exists():
        with open(fastq_file, 'w') as f:
            f.write(f"@{sample}_SEQ_1\n")
            f.write("ACGTACGTACGT\n")
            f.write("+\n")
            f.write("IIIIIIIIIIII\n")
            f.write(f"@{sample}_SEQ_2\n")
            f.write("GCTAGCTAGCTA\n")
            f.write("+\n")
            f.write("IIIIIIIIIIII\n")

# Create a test contrasts file if it doesn't exist
contrasts_file = input_dir / "test_contrasts.txt"
if not contrasts_file.exists():
    with open(contrasts_file, 'w') as f:
        f.write("contrast,control,treatment\n")
        f.write("test_contrast,sample1,sample2\n")

print(f"Test data prepared in: {test_output_dir}")
print(f"- Input directory: {input_dir}")
print(f"- Library file: {library_file}")
print(f"- Contrasts file: {contrasts_file}")
print(f"- FASTQ files: {list(fastq_dir.glob('*.fastq'))}")

# Run the pipeline
print("\n=== Running CRISPR Analysis Pipeline ===\n")
try:
    # Call the pipeline function directly with the correct parameters
    result = run_pipeline(
        input_dir=str(input_dir),
        output_dir=str(results_dir),
        experiment_name="test_run",
        library_file=str(library_file),
        contrasts_file=str(contrasts_file),
        norm_method="median",
        fdr_threshold=0.05,
        overwrite=True,
        skip_drugz=True,  # Skip DrugZ analysis for testing
        skip_mle=True,    # Skip MLE analysis for testing
        use_docker=False, # Skip Docker for testing
        parallel=False    # Disable parallelism for simpler debugging
    )
    
    print(f"\nPipeline run completed successfully!")
    print(f"Results are in: {results_dir}")
    
except Exception as e:
    print(f"Error running pipeline: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1) 