# CRISPR Analysis Pipeline

A modular pipeline for CRISPR screening data analysis, supporting both MAGeCK and DrugZ methods.

## Features

- Individual sample processing with Docker for better parallelization
- Comprehensive quality control and visualization
- Robust error handling and diagnostics
- Support for both read count tables and raw FASTQ files
- Flexible input formats with automatic validation and fixing

## Directory Structure

The pipeline is organized into the following modules:

- `core/` - Core pipeline functionality
  - `config.py` - Configuration and constants
  - `logging_setup.py` - Logging configuration
  - `validation.py` - Input validation functions
  - `file_handling.py` - File and directory management utilities

- `docker/` - Docker-related functionality
  - `docker_utils.py` - Docker container management utilities
  - `mageck_docker.py` - MAGeCK Docker functions
  - `drugz_utils.py` - DrugZ execution utilities

- `analysis/` - Analysis functions
  - `mageck_count.py` - Functions for MAGeCK count operations
  - `mageck_test.py` - Functions for MAGeCK test operations
  - `drugz.py` - Functions for DrugZ analysis

- `qc/` - Quality control components
  - `plot_qa.py` - QA/QC plotting functions
  - `correlation_plots.py` - Guide- and gene-level correlation plots
  - `statistics.py` - Statistical analysis utilities

- `utils/` - General utility functions
  - `sample_utils.py` - Sample sheet generation and management
  - `diagnostic_utils.py` - Diagnostic and troubleshooting tools

## Main Entrypoints

- `run_pipeline.py` - Main pipeline script
- `run_individual_sample.py` - Process a single sample (for parallelization)
- `run_qc.py` - Standalone QC analysis

## Usage

```bash
# Run full pipeline
python run_pipeline.py --input /path/to/input --output /path/to/output --library /path/to/library.csv

# Process individual sample
python run_individual_sample.py --fastq /path/to/sample.fastq --library /path/to/library.csv --output /path/to/output

# Run QC on existing results
python run_qc.py --input /path/to/results --output /path/to/qc_output
```

## Requirements

- Python 3.8+
- Docker
- Required Python packages are listed in requirements.txt 