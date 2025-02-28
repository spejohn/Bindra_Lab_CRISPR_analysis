# CRISPR Analysis Pipeline Testing Guide

This document provides instructions for testing the CRISPR Analysis Pipeline with both test data and real datasets.

## Setup and Installation

Before running any tests, you need to install the package in development mode:

```bash
python install_dev.py
```

This will install the `analysis_pipeline` package in development mode, making it importable from anywhere in your Python environment.

## Testing Options

There are several ways to test the pipeline:

### 1. Direct Pipeline Testing

The most straightforward way to test the pipeline is to use the `test_direct_pipeline.py` script:

```bash
python test_direct_pipeline.py
```

This script:
- Creates minimal test data (FASTQ files, library file, contrasts file)
- Runs the pipeline directly using the `run_pipeline` function
- Bypasses Docker for testing purposes
- Disables parallelism for easier debugging

### 2. Integration Tests

To run the full suite of integration tests:

```bash
python run_integration_tests.py
```

These tests verify:
- Error handling for missing or malformed files
- Full workflow execution with test data
- Incremental workflow execution

### 3. Snakemake Workflow Testing

To test the Snakemake workflow:

```bash
python run_snakemake_test.py
```

This script:
- Creates test data and configuration
- Runs Snakemake in dry-run mode first
- Optionally runs the full workflow

## Testing with Real Data

When testing with real data:

1. **Start Small**: Begin with a small subset of your data (2-3 samples)
2. **Check Paths**: Ensure all file paths are correct and accessible
3. **Docker Configuration**: Make sure Docker is properly configured if using Docker-based tools
4. **Memory Requirements**: Be aware of memory requirements for larger datasets

Example command for real data:

```bash
python pipeline.py \
    --input-dir /path/to/data \
    --output-dir /path/to/results \
    --experiment-name my_experiment \
    --library-file /path/to/library.txt \
    --contrasts-file /path/to/contrasts.txt
```

## Troubleshooting

### Import Errors

If you encounter import errors:

1. Make sure you've run `python install_dev.py`
2. Check that the package is properly installed: `python -c "import analysis_pipeline; print(analysis_pipeline.__version__)"`
3. If still having issues, try adding the project root to your PYTHONPATH:
   ```bash
   export PYTHONPATH=/path/to/project:$PYTHONPATH  # Linux/Mac
   set PYTHONPATH=C:\path\to\project;%PYTHONPATH%  # Windows
   ```

### Docker Issues

If Docker-related errors occur:

1. Ensure Docker is installed and running
2. Try updating the Docker Python package: `pip install -U docker`
3. For testing without Docker, use the `--use-docker false` flag or set `use_docker=False` in code

### Snakemake Issues

For Snakemake-related problems:

1. Check the Snakefile for rule conflicts
2. Ensure your config file has all required parameters
3. Run with `--dryrun` first to validate the workflow
4. Check the `.snakemake/log` directory for detailed error logs

## Conclusion

By following these testing procedures, you can ensure the CRISPR Analysis Pipeline works correctly with your data before running full-scale analyses. 