# CRISPR Analysis Pipeline

A modular pipeline for CRISPR screening data analysis, supporting MAGeCK (RRA and MLE methods) and DrugZ analysis.

## Features

- Individual sample processing with Docker for better parallelization
- Comprehensive quality control and visualization
- Robust error handling and diagnostics
- Support for both read count tables and raw FASTQ files
- Flexible input formats with automatic validation and fixing
- Conditional MAGeCK MLE analysis based on design matrix availability
- Standardized output directory structure
- **NEW: Snakemake workflow management for reproducible analyses**
- **NEW: Cross-platform compatibility (Windows, macOS, Linux)**
- **NEW: Enhanced Docker volume handling and path resolution**

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
  - `sample_processing.py` - Sample processing and count aggregation
  - `mageck_analysis.py` - MAGeCK RRA and MLE analysis

- `qc/` - Quality control components
  - `plot_quality.py` - QA/QC plotting functions
  - `correlation_plots.py` - Guide- and gene-level correlation plots
  - `statistics.py` - Statistical analysis utilities

- `utils/` - General utility functions
  - `sample_utils.py` - Sample sheet generation and management
  - `diagnostic_utils.py` - Diagnostic and troubleshooting tools

## Output Directory Structure

For each experiment, the pipeline creates a standardized directory structure:

```
<output_dir>/
├── <experiment_name>/
│   ├── logs/                  # Log files
│   ├── counts/                # Count files 
│   ├── samplesheets/          # Generated sample sheets
│   ├── library/               # Validated library files
│   ├── RRA/                   # MAGeCK RRA analysis results
│   ├── MLE/                   # MAGeCK MLE analysis results (when applicable)
│   ├── drugz/                 # DrugZ analysis results 
│   ├── qc/                    # Quality control plots and reports
│   └── figures/               # Publication-ready figures
```

## Main Entrypoints

- `pipeline.py` - Main pipeline script
- `run_individual_sample.py` - Process a single sample (for parallelization)
- `run_qc.py` - Standalone QC analysis

## Usage

```bash
# Run full pipeline
python pipeline.py --input /path/to/input --output /path/to/output --library /path/to/library.csv

# Run pipeline without MLE analysis
python pipeline.py --input /path/to/input --output /path/to/output --library /path/to/library.csv --skip-mle

# Process individual sample
python run_individual_sample.py --fastq /path/to/sample.fastq --library /path/to/library.csv --output /path/to/output

# Run QC on existing results
python run_qc.py --input /path/to/results --output /path/to/qc_output
```

## MAGeCK MLE Analysis

The pipeline automatically detects design matrix files in the same directory as contrast tables. If a design matrix is found, it will run MAGeCK MLE analysis in addition to the standard RRA analysis.

Requirements for design matrix files:
- Must be named following patterns like `design_matrix.txt`, `design.txt`, etc.
- Must be in the same directory as the contrast table
- Must follow the MAGeCK MLE design matrix format with columns for Sample and condition variables

## Multi-Screen Analysis

The pipeline now supports automatic detection and analysis of multiple screens:

```bash
# Run analysis on all screens in a directory
python analysis_pipeline/run_multiple_screens.py /path/to/screens_directory
```

This will:
1. Automatically detect all directories containing CRISPR screen data
2. Find library and contrast files for each screen
3. Determine which screens have already been analyzed
4. Only process screens that need analysis
5. Generate a summary report of all screens

## Smart Path Detection

The pipeline can now automatically find necessary files:

```bash
# Simplified usage with auto-detection of outputs and library
python analysis_pipeline/pipeline.py /path/to/screen_directory
```

The pipeline will:
- Use `<input_dir>/library.csv` as the default library file if not specified
- Use `<input_dir>/contrasts.csv` as the default contrasts file if not specified
- Create output in `<input_dir>/results` by default
- Detect if the screen has already been analyzed and skip if appropriate

## Cross-Platform Compatibility

This pipeline is designed to work seamlessly across Windows, macOS, and Linux operating systems. Special attention has been given to path handling to ensure compatibility with Docker volumes on all platforms.

### Windows Path Handling for Docker

When running on Windows, the pipeline automatically converts Windows-style paths to Docker-compatible paths:

- Windows path: `C:\Users\username\data`
- Docker path: `/c/Users/username/data`

This conversion happens transparently in the background, so you can use regular Windows paths throughout the pipeline.

### Path Conversion Function

The pipeline provides a utility function for path conversion in `core/utils.py`:

```python
from analysis_pipeline.core.utils import convert_win_path_to_docker

# Convert a Windows path to Docker-compatible format
windows_path = r"C:\Users\username\data"
docker_path = convert_win_path_to_docker(windows_path)
# Result: "/c/Users/username/data"
```

### Testing Path Conversion

You can verify path conversion is working correctly by running the included test script:

```bash
python analysis_pipeline/test_docker_paths.py
```

This will output the conversion results for various test paths and your current working directory.

### Docker Volume Mounting

The pipeline handles Docker volume mounting differently depending on the operating system:

- **Windows**: Automatically converts paths for Docker compatibility
- **macOS/Linux**: Uses paths directly with minimal modification

This ensures that file access within Docker containers works correctly regardless of the host operating system.

## Requirements

- Python 3.8+
- Docker
- Required Python packages are listed in requirements.txt

## Docker Containers

The pipeline now supports containerized execution of analysis tools, which provides several advantages:

- Consistent execution environment across different systems
- No need to install tools and dependencies locally
- Easier deployment and reproducibility

### Available Containers

- **MAGeCK Container**: For MAGeCK RRA and MLE analyses
- **DrugZ Container**: For DrugZ analysis

### Building Docker Images

To build the Docker images, run:

```bash
python analysis_pipeline/docker/build_images.py
```

Or manually:

```bash
cd analysis_pipeline/docker
docker-compose build
```

### Using Docker Containers

Docker containers are used automatically when running the pipeline if Docker is available. To enable Docker usage:

```bash
python analysis_pipeline/pipeline.py \
    --input-dir /path/to/fastq \
    --output-dir /path/to/output \
    --library-file /path/to/library.csv \
    --experiment-name my_experiment \
    --use-docker
```

You can also use the DrugZ Docker container directly:

```bash
python analysis_pipeline/analysis/run_drugz_docker.py \
    --input /path/to/count_file.csv \
    --design /path/to/design_file.txt \
    --output /path/to/output \
    --prefix my_experiment
```

### Docker Directory Structure

```
analysis_pipeline/docker/
├── docker-compose.yml           # Docker Compose configuration
├── build_images.py              # Script to build Docker images
└── dockerfiles/
    ├── Dockerfile.mageck        # Dockerfile for MAGeCK
    └── Dockerfile.drugz         # Dockerfile for DrugZ
```

## Snakemake Workflow

The pipeline now includes a Snakemake workflow for improved reproducibility and dependency management:

```bash
# Run the workflow with default settings
python analysis_pipeline/run_snakemake.py /path/to/screens_directory

# Run with multiple cores
python analysis_pipeline/run_snakemake.py /path/to/screens_directory -j 8

# Do a dry run to see what would be executed
python analysis_pipeline/run_snakemake.py /path/to/screens_directory --dryrun

# Specify custom output directory
python analysis_pipeline/run_snakemake.py /path/to/screens_directory -o /path/to/output
```

### Advantages of Snakemake

- **Automatic dependency tracking**: Only runs necessary steps based on file timestamps
- **Parallel execution**: Process multiple samples simultaneously
- **Resume functionality**: Pick up where you left off if a run is interrupted
- **Workflow visualization**: Generate DAG visualizations of your workflow
- **Reproducibility**: Ensures the same steps are run in the same order

### Configuration

The workflow can be configured using a `config.yaml` file or command-line parameters:

```yaml
# Example config.yaml
input_dir: "input"
output_dir: "results"
threads: 8
skip_drugz: false
use_docker: true
```

To use a custom config file:

```bash
python analysis_pipeline/run_snakemake.py /path/to/screens_directory --configfile /path/to/config.yaml
```

## Installation

### Dependencies

Install the required dependencies:

```bash
pip install -r requirements.txt
```

For Snakemake to work properly on Windows, you may need to install dependencies using Conda:

```bash
conda install -c bioconda -c conda-forge snakemake
```

### Windows Compatibility

The pipeline now includes enhanced Windows support:

- **Path Handling**: Automatic conversion of Windows paths to Docker-compatible format
- **Volume Mounting**: Improved Docker volume mounting for Windows host paths
- **Debugging**: Additional logging for Docker command execution and path conversion
- **Docker Desktop**: Requires Docker Desktop on Windows with WSL 2 backend for optimal performance

To verify Docker is properly configured on Windows:

```bash
python -c "from analysis_pipeline.docker.docker_utils import verify_docker; verify_docker()"
```

### Note on Linter Errors in Snakefile

The Snakefile may show linter errors in some IDEs due to Snakemake's unique syntax, which extends Python. These errors can be safely ignored as they don't affect the functionality of the workflow. Common false positive errors include:

- Triple-quoted shell commands ("""...""")
- The use of rule blocks and wildcards
- Snakemake-specific directives like configfile, rule, and input/output sections

If modifying the Snakefile, it's recommended to test changes with a dry run:

```bash
python analysis_pipeline/run_snakemake.py --dryrun /path/to/screen_directory
```